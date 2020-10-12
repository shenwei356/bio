package fastx

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"regexp"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/xopen"
)

// ErrNotFASTXFormat means that the file is not FASTA/Q
var ErrNotFASTXFormat = errors.New("fastx: invalid FASTA/Q format")

// ErrBadFASTQFormat means bad fastq format
var ErrBadFASTQFormat = errors.New("fastx: bad fastq format")

// ErrUnequalSeqAndQual means unequal sequence and quality
var ErrUnequalSeqAndQual = errors.New("fastx: unequal sequence and quality")

// ErrNoContent means nothing in the file or stream
var ErrNoContent = errors.New("fastx: no content found")

var pageSize = os.Getpagesize() * 2

// Reader seamlessly parse both FASTA and FASTQ formats
type Reader struct {
	fh *xopen.Reader // file handle, xopen is such a wonderful package

	buf                []byte // for store readed data from fh
	r                  int
	buffer             *bytes.Buffer // buffer of a record
	needMoreCheckOfBuf bool
	lastByte           byte
	checkSeqType       bool
	lastPart           bool
	finished           bool

	firstseq bool // for guess alphabet by the first seq
	delim    byte
	IsFastq  bool // if the file is fastq format

	t        *seq.Alphabet  // alphabet
	IDRegexp *regexp.Regexp // regexp for parsing record id

	head, seq, qual []byte
	seqBuffer       *bytes.Buffer
	qualBuffer      *bytes.Buffer
	record          *Record

	Err error // Current error
}

// regexp for checking idRegexp string.
// The regular expression must contain "(" and ")" to capture matched ID
var reCheckIDregexpStr = regexp.MustCompile(`\(.+\)`)

// DefaultIDRegexp is the default ID parsing regular expression
var DefaultIDRegexp = `^(\S+)\s?`
var isUsingDefaultIDRegexp bool

// NewDefaultReader automaticlly recognizes sequence type and parses id with default manner
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
		isUsingDefaultIDRegexp = true
	} else {
		if !reCheckIDregexpStr.MatchString(idRegexp) {
			return nil, fmt.Errorf(`fastx: regular expression must contain "(" and ")" to capture matched ID. default: %s`, DefaultIDRegexp)
		}
		var err error
		r, err = regexp.Compile(idRegexp)
		if err != nil {
			return nil, fmt.Errorf("fastx: fail to compile regexp: %s", idRegexp)
		}
		if idRegexp == DefaultIDRegexp {
			isUsingDefaultIDRegexp = true
		}
	}

	fh, err := xopen.Ropen(file)
	if err != nil {
		if err == xopen.ErrNoContent {
			return &Reader{
					fh:       nil,
					finished: true,
					lastPart: true,
					Err:      io.EOF,
				},
				nil
		}
		return nil, fmt.Errorf("fastx: %s", err)
	}

	fastxReader := &Reader{
		fh:           fh,
		buf:          make([]byte, pageSize),
		t:            t,
		IDRegexp:     r,
		firstseq:     true,
		checkSeqType: true,
	}
	fastxReader.buffer = bytes.NewBuffer(make([]byte, 0, defaultBytesBufferSize))
	fastxReader.seqBuffer = bytes.NewBuffer(make([]byte, 0, defaultBytesBufferSize))
	fastxReader.qualBuffer = bytes.NewBuffer(make([]byte, 0, defaultBytesBufferSize))

	return fastxReader, nil
}

// Close closes the reader
func (fastxReader *Reader) Close() {
	if fastxReader.fh != nil {
		fastxReader.fh.Close()
	}
}

// Read reads and return one FASTA/Q record.
// Note that, similar to bytes.Buffer.Bytes() method,
// the current record will change after your another call of this method.
// So, you could use record.Clone() to make a copy.
func (fastxReader *Reader) Read() (*Record, error) {
	fastxReader.read()
	return fastxReader.record, fastxReader.Err
}

func (fastxReader *Reader) read() {
	if fastxReader.lastPart && fastxReader.finished {
		fastxReader.Err = io.EOF
		return
	}
	if fastxReader.Err != nil {
		return
	}

	var n int
	var err error
	var p []byte

	for {
		if !fastxReader.needMoreCheckOfBuf && !fastxReader.lastPart {
			n, err = fastxReader.fh.Read(fastxReader.buf)
			if err != nil {
				if err == io.EOF {
					fastxReader.lastPart = true
				} else {
					fastxReader.Err = err
					fastxReader.Close()
					return
				}
			} else if n == 0 {
				fastxReader.Err = io.ErrUnexpectedEOF
				fastxReader.Close()
				return
			}

			// last part of file OR just because reader not fulfill the buf,
			// like reading from stdin
			if n < len(fastxReader.buf) {
				// fastxReader.lastPart = true
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
				case '\n': // allow some lines
					pn++
					if pn > 100 {
						if i > 10240 { // ErrNotFASTXFormat
							fastxReader.Err = ErrNotFASTXFormat
							fastxReader.Close()
							return
						}
					}
					// break FORCHECK
				default: // not typical FASTA/Q
					// if i > 10240 || fastxReader.lastPart { // ErrNotFASTXFormat
					fastxReader.Err = ErrNotFASTXFormat
					fastxReader.Close()
					return
					// }
				}
			}
			fastxReader.checkSeqType = false
		}

		var shorterQual bool
		var i int
	FORSEARCH:
		for {
			if i = bytes.IndexByte(fastxReader.buf[fastxReader.r:], fastxReader.delim); i >= 0 {
				if i > 0 {
					fastxReader.lastByte = fastxReader.buf[fastxReader.r+i-1]
				} else {
					p = fastxReader.buffer.Bytes()
					if len(p) == 0 {
						fastxReader.lastByte = '\x00'
					} else {
						fastxReader.lastByte = p[len(p)-1]
					}
				}
				if fastxReader.lastByte == '\n' { // yes!
					if i > 0 {
						fastxReader.buffer.Write(dropCR(fastxReader.buf[fastxReader.r : fastxReader.r+i-1]))
					} else {
						fastxReader.buffer.WriteByte('\n')
					}

					// we have to avoid the case of quality line starts with "@"
					shorterQual, err = fastxReader.parseRecord()
					if fastxReader.IsFastq && err != nil && err == ErrUnequalSeqAndQual {
						if shorterQual {
							fastxReader.buffer.WriteByte('\n')
							fastxReader.buffer.WriteByte(fastxReader.delim)
							fastxReader.needMoreCheckOfBuf = true
							fastxReader.r += i + 1
							continue FORSEARCH
						}
						fastxReader.Err = ErrBadFASTQFormat
						return
					}
					fastxReader.buffer.Reset()
					fastxReader.needMoreCheckOfBuf = true
					fastxReader.r += i + 1
					return
				}
				// inline >/@
				fastxReader.buffer.Write(fastxReader.buf[fastxReader.r : fastxReader.r+i+1])
				fastxReader.r += i + 1
				fastxReader.needMoreCheckOfBuf = true
				continue FORSEARCH
			}

			fastxReader.buffer.Write(fastxReader.buf[fastxReader.r:])
			if fastxReader.lastPart {
				_, err = fastxReader.parseRecord()
				if err != nil { // no any chance
					fastxReader.Err = err
					fastxReader.Close()
					return
				}
				fastxReader.buffer.Reset()
				fastxReader.Close()
				fastxReader.finished = true
				return
			}
			fastxReader.needMoreCheckOfBuf = false
			break FORSEARCH
		}
	}
}

// parseRecord parse a FASTA/Q record from the fastxReader.buffer
func (fastxReader *Reader) parseRecord() (bool, error) {
	fastxReader.seqBuffer.Reset()
	if fastxReader.IsFastq {
		fastxReader.qualBuffer.Reset()
	}

	var p = fastxReader.buffer.Bytes()
	if j := bytes.IndexByte(p, '\n'); j > 0 {
		fastxReader.head = dropCR(p[0:j])
		r := j + 1

		if !fastxReader.IsFastq { // FASTA
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
		} else { // FASTQ
			var isQual bool
			for {
				if k := bytes.IndexByte(p[r:], '\n'); k >= 0 {
					if k > 0 && p[r] == '+' && !isQual {
						isQual = true
					} else if isQual {
						fastxReader.qualBuffer.Write(dropCR(p[r : r+k]))
					} else {
						fastxReader.seqBuffer.Write(dropCR(p[r : r+k]))
					}
					r += k + 1
					continue
				}
				if isQual {
					fastxReader.qualBuffer.Write(dropCR(p[r:]))
				}
				break
			}

			// may be the case of quality line starts with "@"
			if fastxReader.seqBuffer.Len() != fastxReader.qualBuffer.Len() {
				return fastxReader.seqBuffer.Len() > fastxReader.qualBuffer.Len(), ErrUnequalSeqAndQual
			}

			fastxReader.seq = fastxReader.seqBuffer.Bytes()
			fastxReader.qual = fastxReader.qualBuffer.Bytes()
		}

	} else {
		fastxReader.head = dropCR(dropLF(p))
		fastxReader.seq = []byte{}
		fastxReader.qual = []byte{}
	}

	// guess alphabet
	if fastxReader.firstseq {
		if fastxReader.t == nil {
			fastxReader.t = seq.GuessAlphabetLessConservatively(fastxReader.seq)
		}
		fastxReader.firstseq = false
	}

	if len(fastxReader.head) == 0 && len(fastxReader.seq) == 0 {
		return false, io.EOF
	}
	// new record
	if fastxReader.IsFastq {
		fastxReader.record, fastxReader.Err = NewRecordWithQual(fastxReader.t,
			fastxReader.parseHeadID(fastxReader.head), fastxReader.head,
			fastxReader.parseHeadDesc(fastxReader.head),
			fastxReader.seq, fastxReader.qual)
	} else {
		fastxReader.record, fastxReader.Err = NewRecord(fastxReader.t,
			fastxReader.parseHeadID(fastxReader.head), fastxReader.head,
			fastxReader.parseHeadDesc(fastxReader.head),
			fastxReader.seq)
	}

	if fastxReader.Err != nil {
		fastxReader.Close()
	}

	return false, fastxReader.Err
}

func (fastxReader *Reader) parseHeadID(head []byte) []byte {
	return parseHeadID(fastxReader.IDRegexp, head)
}

func (fastxReader *Reader) parseHeadDesc(head []byte) []byte {
	return parseHeadDesc(fastxReader.IDRegexp, head)
}

// ParseHeadID parse ID from head by IDRegexp
func ParseHeadID(idRegexp *regexp.Regexp, head []byte) []byte {
	found := idRegexp.FindSubmatch(head)
	if found == nil { // not match
		return head
	}
	return found[1]
}

// parseHeadID parse ID from head by IDRegexp
func parseHeadID(idRegexp *regexp.Regexp, head []byte) []byte {
	if isUsingDefaultIDRegexp {
		if i := bytes.IndexByte(head, ' '); i > 0 {
			return head[0:i]
		}
		if i := bytes.IndexByte(head, '\t'); i > 0 {
			return head[0:i]
		}
		return head
	}

	found := idRegexp.FindSubmatch(head)
	if found == nil { // not match
		return head
	}
	return found[1]
}

// parseHeadDesc returns description of header
func parseHeadDesc(idRegexp *regexp.Regexp, head []byte) []byte {
	if isUsingDefaultIDRegexp {
		if i := bytes.IndexByte(head, ' '); i > 0 {
			return bytes.TrimLeft(head[i:], " \t")
		}
		if i := bytes.IndexByte(head, '\t'); i > 0 {
			return bytes.TrimLeft(head[i:], " \t")
		}
		return []byte{}
	}
	return []byte{}
}

// Alphabet returns Alphabet of the file
func (fastxReader *Reader) Alphabet() *seq.Alphabet {
	if fastxReader.t == nil {
		return seq.Unlimit
	}
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

// RecordChunk is chunk for records
type RecordChunk struct {
	ID   uint64
	Data []*Record
	Err  error
}

// ChunkChan asynchronously reads FASTA/Q records, and returns a channel of
// Record Chunk, from which you can easily access the records.
// bufferSize is the number of buffered chunks, and chunkSize is the size
// of records in a chunk.
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
					if i == 0 { // no any seqs
						close(ch)
						return
					}
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
