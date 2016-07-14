package fastx

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"regexp"
	"syscall"

	"github.com/brentp/xopen"
	"github.com/shenwei356/bio/seq"
)

// ErrClosed means that the reading process is canceled
var ErrClosed = errors.New("bio.seqio.fastx: reading canceled")

// ErrNotFASTXFormat means that the file is not FASTA/Q
var ErrNotFASTXFormat = errors.New("bio.seqio.fastx: not FASTA/Q format")

// ErrBadFASTQFormat means bad fastq format
var ErrBadFASTQFormat = errors.New("bio.seqio.fastx: bad fastq format")

// Reader asynchronously parse FASTX file with buffer,
// each buffer contain a chunk of multiple fastx records (RecordChunk).
// Reader also support safe cancellation.
type Reader struct {
	fh                 *xopen.Reader // file handle, xopen is such a wonderful package
	lastPart           bool
	buf                []byte
	buffer             *bytes.Buffer
	needMoreCheckOfBuf bool
	lastByte           byte

	t        *seq.Alphabet  // alphabet
	IDRegexp *regexp.Regexp // regexp for parsing record id

	checkSeqType    bool
	firstseq        bool // for guess alphabet by the first seq
	r               int
	head, seq, qual []byte
	seqBuffer       *bytes.Buffer
	record          *Record

	delim   byte
	IsFastq bool

	Err      error
	finished bool
}

// regexp for checking idRegexp string.
// The regular expression must contain "(" and ")" to capture matched ID
var reCheckIDregexpStr = regexp.MustCompile(`\(.+\)`)

// DefaultIDRegexp is the default ID parsing regular expression
var DefaultIDRegexp = `^([^\s]+)\s?`

// NewDefaultReader automaticly recognizes sequence type and parses id with default manner
func NewDefaultReader(file string) (*Reader, error) {
	return NewReader(nil, file, "")
}

// NewReader is constructor of FASTX Reader.
//
// Parameters:
//
//        t            sequence alphabet
//                     if nil is given, it will guess alphabet by the first record
//        file         file name, "-" for stdin
//        idRegexp     id parsing regular expression string, must contains "(" and ")" to capture matched ID
//                     "" for default value: `^([^\s]+)\s?`
//                     if record head does not match the idRegxp, whole name will be the id
//
func NewReader(t *seq.Alphabet, file string, idRegexp string) (*Reader, error) {
	var r *regexp.Regexp
	if idRegexp == "" {
		r = regexp.MustCompile(DefaultIDRegexp)
	} else {
		if !reCheckIDregexpStr.MatchString(idRegexp) {
			return nil, fmt.Errorf(`bio.seqio.fastx: regular expression must contain "(" and ")" to capture matched ID. default: %s`, DefaultIDRegexp)
		}
		var err error
		r, err = regexp.Compile(idRegexp)
		if err != nil {
			return nil, fmt.Errorf("bio.seqio.fastx: fail to compile regexp: %s", idRegexp)
		}
	}

	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, fmt.Errorf("bio.seqio.fastx: %s", err)
	}

	fastxReader := &Reader{
		fh:           fh,
		buf:          make([]byte, syscall.Getpagesize()),
		t:            t,
		IDRegexp:     r,
		firstseq:     true,
		checkSeqType: true,
	}
	fastxReader.buffer = bytes.NewBuffer(make([]byte, 0, defaultBytesBufferSize))
	fastxReader.seqBuffer = bytes.NewBuffer(make([]byte, 0, defaultBytesBufferSize))

	// fastxReader.read()

	return fastxReader, nil
}

func (fastxReader *Reader) close() {
	fastxReader.fh.Close()
}

// Next reads and return the record
func (fastxReader *Reader) Read() (*Record, error) {
	fastxReader.read()
	return fastxReader.record, fastxReader.Err
}

func (fastxReader *Reader) read() {
	if fastxReader.lastPart && fastxReader.finished {
		// fmt.Println("shit---------------------------------")
		fastxReader.Err = io.EOF
		return
	}
	if fastxReader.Err != nil {
		return
	}

	var n int
	var err error

	for {
		if !fastxReader.needMoreCheckOfBuf && !fastxReader.lastPart {
			// fmt.Println("----- read data ------")
			n, err = fastxReader.fh.Read(fastxReader.buf)
			if err != nil { // Error
				if err == io.EOF {
					// fmt.Println("----- eof of reading ------")
					fastxReader.lastPart = true
				} else {
					fastxReader.Err = err
					fastxReader.close()
					return
				}
			}

			if n < len(fastxReader.buf) { // last part of file
				// fmt.Println("----- eof ------", n, len(fastxReader.buf))
				fastxReader.lastPart = true
				fastxReader.buf = fastxReader.buf[0:n] // very important!
			}
			fastxReader.r = 0 /// TO CHECK
		}

		if fastxReader.checkSeqType {
			pn := 0
		FORCHECK:
			for i := range fastxReader.buf {
				switch fastxReader.buf[i] {
				case '>':
					fastxReader.checkSeqType = false
					fastxReader.IsFastq = false
					fastxReader.delim = '>'
					fastxReader.r = i + 1
					break FORCHECK
				case '@':
					fastxReader.IsFastq = true
					fastxReader.delim = '@'
					fastxReader.r = i + 1
					break FORCHECK
				case '\n': // allow some line
					pn++
					if pn > 100 {
						if i > 10240 { // ErrNotFASTXFormat
							fastxReader.Err = ErrNotFASTXFormat
							fastxReader.close()
							return
						}
					}
					// break FORCHECK
				default: // not typical FASTA/Q
					if i > 10240 { // ErrNotFASTXFormat
						fastxReader.Err = ErrNotFASTXFormat
						fastxReader.close()
						return
					}
				}
			}
			fastxReader.checkSeqType = false
		}

	FORFOUND:
		for {
			// fmt.Printf("======r:%d, buffer: =%s=, buf: =%s=\n", fastxReader.r,
			// 	fastxReader.buffer.Bytes(), fastxReader.buf[fastxReader.r:])
			if i := bytes.IndexByte(fastxReader.buf[fastxReader.r:], fastxReader.delim); i >= 0 {
				if i > 0 {
					fastxReader.lastByte = fastxReader.buf[fastxReader.r+i-1]
				}
				if fastxReader.lastByte == '\n' { // yes!
					// fmt.Println("----- found ------")
					if i > 0 {
						fastxReader.buffer.Write(dropCR(fastxReader.buf[fastxReader.r : fastxReader.r+i-1]))
					}

					fastxReader.parseRecord()
					fastxReader.buffer.Reset()
					fastxReader.needMoreCheckOfBuf = true
					fastxReader.r += i + 1
					return
				}
				// fmt.Println("----- inline ------")
				// inline ">"
				fastxReader.buffer.Write(fastxReader.buf[fastxReader.r : fastxReader.r+i+1])
				fastxReader.r += i + 1
				fastxReader.needMoreCheckOfBuf = true
				continue FORFOUND
			}

			fastxReader.buffer.Write(fastxReader.buf[fastxReader.r:])
			if fastxReader.lastPart {
				// fmt.Println("----- last part ------")
				fastxReader.parseRecord()
				fastxReader.buffer.Reset()
				fastxReader.close()
				fastxReader.finished = true
				return
			}
			// fmt.Println("----- need more data ------")
			fastxReader.needMoreCheckOfBuf = false
			break FORFOUND
		}
	}
}

func (fastxReader *Reader) parseRecord() {
	fastxReader.seqBuffer.Reset()

	var p = fastxReader.buffer.Bytes()
	if j := bytes.IndexByte(p, '\n'); j > 0 {
		fastxReader.head = dropCR(p[0:j])

		r := j + 1
		for {
			if k := bytes.IndexByte(p[r:], '\n'); k >= 0 {
				fastxReader.seqBuffer.Write(dropCR(p[r : r+k]))
				r += k + 1
				continue
			}
			fastxReader.seqBuffer.Write(dropCR(p[r:]))
			break
		}
		fastxReader.seq = fastxReader.seqBuffer.Bytes()
	} else {
		fastxReader.head = p
		fastxReader.seq = []byte{}
	}

	// fmt.Printf("head:%s\nseq:=%s=\nseq len:%d\n", fastxReader.head, fastxReader.seq, len(fastxReader.seq))

	// guess alphabet
	if fastxReader.firstseq {
		if fastxReader.t == nil {
			fastxReader.t = seq.GuessAlphabetLessConservatively(fastxReader.seq)
		}
		fastxReader.firstseq = false
	}

	// new record
	fastxReader.record, fastxReader.Err = NewRecord(fastxReader.t,
		fastxReader.parseHeadID(fastxReader.head), fastxReader.head, fastxReader.seq)
	if fastxReader.Err != nil {
		fastxReader.close()
	}
}

func (fastxReader *Reader) parseHeadID(head []byte) []byte {
	return ParseHeadID(fastxReader.IDRegexp, head)
}

// ParseHeadID parse ID from head by IDRegexp
func ParseHeadID(idRegexp *regexp.Regexp, head []byte) []byte {
	found := idRegexp.FindAllSubmatch(head, -1)
	if found == nil { // not match
		return head
	}
	return found[0][1]
}

// Cancel method cancel the reading process
func (fastxReader *Reader) Cancel() {
	fastxReader.close()
}

// Alphabet returns Alphabet of the file
func (fastxReader *Reader) Alphabet() *seq.Alphabet {
	return fastxReader.t
}

func dropCR(data []byte) []byte {
	if len(data) > 0 && data[len(data)-1] == '\r' {
		return data[0 : len(data)-1]
	}
	return data
}

func dropLF(data []byte) []byte {
	if len(data) > 0 && data[len(data)-1] == '\n' {
		return data[0 : len(data)-1]
	}
	return data
}

// -------------------------------------------------

// RecordChunk  is
type RecordChunk struct {
	ID   uint64
	Data []*Record
	Err  error
}

func (fastxReader *Reader) ChunkChan(bufferSize int, chunkSize int) chan RecordChunk {
	var ch chan RecordChunk
	if bufferSize <= 0 {
		ch = make(chan RecordChunk)
	} else {
		ch = make(chan RecordChunk, bufferSize)
	}
	if chunkSize < 1 {
		chunkSize = 1
	}

	go func() {
		var i int
		var id uint64
		chunkData := make([]*Record, chunkSize)

		for {
			record, err := fastxReader.Read()
			if err != nil {
				if err == io.EOF {
					break
				}
				ch <- RecordChunk{id, chunkData[0:i], err}
				close(ch)
				return
			}
			chunkData[i] = record.Clone()
			i++

			if i == chunkSize {
				ch <- RecordChunk{id, chunkData[0:i], nil}
				id++
				i = 0
				chunkData = make([]*Record, chunkSize)
			}
		}

		ch <- RecordChunk{id, chunkData[0:i], nil}
		close(ch)
		return
	}()

	return ch
}
