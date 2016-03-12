package fasta

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"regexp"

	"github.com/brentp/xopen"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/byteutil"
)

// FastaRecord struct
type FastaRecord struct {
	ID   []byte   // id
	Name []byte   // full name
	Seq  *seq.Seq // seq
}

func (fastRecord FastaRecord) String() string {
	return fmt.Sprintf(">%s\n%s", fastRecord.Name,
		byteutil.WrapByteSlice(fastRecord.Seq.Seq, 70))
}

// NewFastaRecord is constructor of type FastaRecord
func NewFastaRecord(t *seq.Alphabet, id, name, s []byte) (*FastaRecord, error) {
	seq, err := seq.NewSeq(t, s)
	if err != nil {
		return nil, fmt.Errorf("error when parsing seq: %s (%s)", id, err)
	}
	return &FastaRecord{id, name, seq}, nil
}

// FormatSeq formats output, i.e. wrap string into fixed width
func (fastRecord *FastaRecord) FormatSeq(width int) []byte {
	return byteutil.WrapByteSlice(fastRecord.Seq.Seq, width)
}

// FastaRecordChunk  is
type FastaRecordChunk struct {
	ID   uint64
	Data []*FastaRecord
	Err  error
}

// FastaReader asynchronously parse FASTA file with buffer,
// each buffer contain a chunk of multiple fasta records (FastaRecordChunk).
// FastaReader also support safe cancellation.
type FastaReader struct {
	t  *seq.Alphabet // alphabet
	fh *xopen.Reader // file handle, xopen is such a wonderful package

	BufferSize int                   // buffer size
	ChunkSize  int                   // chunk size
	Ch         chan FastaRecordChunk // chanel for output of records chunks
	IDRegexp   *regexp.Regexp        // regexp ofr parsing record id

	firstseq  bool          // for guess alphabet by the first seq
	done      chan struct{} // for cancellation
	finished  bool          // for cancellation
	cancelled bool          // for cancellation
}

// regexp for checking idRegexp string.
// The regular expression must contains "(" and ")" to capture matched ID
var reCheckIDregexpStr = regexp.MustCompile(`\(.+\)`)

// DefaultIDRegexp is the default ID parsing regular expression
var DefaultIDRegexp = `^([^\s]+)\s?`

// NewFastaReader is constructor of FastaReader.
//
// Parameters:
//
//        t            sequence alphabet
//                     if nil is given, it will guess alphabet by the first record
//        file         file name, "-" for stdin
//        bufferSize   buffer size
//        chunkSize    chunk size
//        idRegexp     id parsing regular expression string, must contains "(" and ")" to capture matched ID
//                     "" for default value: `^([^\s]+)\s?`
//                     if record head does not match the idRegxp, whole name will be the id
//
func NewFastaReader(t *seq.Alphabet, file string, bufferSize int, chunkSize int, idRegexp string) (*FastaReader, error) {
	if bufferSize < 0 {
		bufferSize = 0
	}
	if chunkSize < 1 {
		chunkSize = 1
	}

	var r *regexp.Regexp
	if idRegexp == "" {
		r = regexp.MustCompile(DefaultIDRegexp)
	} else {
		if !reCheckIDregexpStr.MatchString(idRegexp) {
			return nil, fmt.Errorf(`regular expression must contains "(" and ")" to capture matched ID. default: %s`, DefaultIDRegexp)
		}
		var err error
		r, err = regexp.Compile(idRegexp)
		if err != nil {
			return nil, fmt.Errorf("fail to compile regexp: %s", idRegexp)
		}
	}

	fh, err := xopen.Ropen(file)
	if err != nil {
		return nil, err
	}

	fastaReader := &FastaReader{
		t:          t,
		fh:         fh,
		BufferSize: bufferSize,
		ChunkSize:  chunkSize,
		Ch:         make(chan FastaRecordChunk, bufferSize),
		IDRegexp:   r,
		firstseq:   true,
		done:       make(chan struct{}),
		finished:   false,
		cancelled:  false,
	}

	fastaReader.read()

	return fastaReader, nil
}

var reTrimRightSpace = regexp.MustCompile(`[\r\n]+$`)
var reTrimSpace = regexp.MustCompile(`[\r\n\s]+`)

// ErrorCanceled means that the reading process is canceled
var ErrorCanceled = errors.New("reading canceled")

// SplitFASTA is split function for FASTA format
func SplitFASTA(data []byte, atEOF bool) (advance int, token []byte, err error) {
	if atEOF && len(data) == 0 {
		return 0, nil, nil
	}
	i := bytes.IndexByte(data, '>')
	if i == 0 { // first '>'
		return i + 1, nil, nil
	}
	if i > 0 {
		return i + 1, dropCR(data[0:i]), nil
	}
	if atEOF {
		return len(data), dropCR(data), nil
	}
	return 0, nil, nil
}

// dropCR drops a terminal \r from the data.
func dropCR(data []byte) []byte {
	if len(data) > 0 && data[len(data)-1] == '\r' {
		return data[0 : len(data)-1]
	}
	return data
}

func (fastaReader *FastaReader) read() {
	go func() {
		var text []byte
		var i int
		var id uint64
		chunkData := make([]*FastaRecord, fastaReader.ChunkSize)

		scanner := bufio.NewScanner(fastaReader.fh)
		scanner.Split(SplitFASTA)
		for scanner.Scan() {
			select {
			case <-fastaReader.done:
				if !fastaReader.finished {
					fastaReader.Ch <- FastaRecordChunk{id, chunkData[0:i], ErrorCanceled}

					fastaReader.finished = true
					fastaReader.fh.Close()
					close(fastaReader.Ch)
					return
				}
			default:
			}

			text = scanner.Bytes()
			j := bytes.IndexByte(text, '\n')

			name := []byte(string(text[0:j]))
			sequence := []byte(string(reTrimSpace.ReplaceAll(text[j+1:], []byte(""))))

			if fastaReader.firstseq {
				if fastaReader.t == nil {
					fastaReader.t = seq.GuessAlphabetLessConservatively(sequence)
				}
				fastaReader.firstseq = false
			}

			fastaRecord, err := NewFastaRecord(fastaReader.t, fastaReader.parseHeadID(name), name, sequence)
			if err != nil {
				fastaReader.Ch <- FastaRecordChunk{id, chunkData[0:i], err}
				fastaReader.fh.Close()
				close(fastaReader.Ch)
				return
			}

			chunkData[i] = fastaRecord
			i++
			if i == fastaReader.ChunkSize {
				fastaReader.Ch <- FastaRecordChunk{id, chunkData[0:i], nil}
				id++
				i = 0
				chunkData = make([]*FastaRecord, fastaReader.ChunkSize)
			}

		}

		fastaReader.Ch <- FastaRecordChunk{id, chunkData[0:i], scanner.Err()}
		fastaReader.fh.Close()
		close(fastaReader.Ch)
	}()
}

func (fastaReader *FastaReader) parseHeadID(head []byte) []byte {
	found := fastaReader.IDRegexp.FindAllSubmatch(head, -1)
	if found == nil { // not match
		return head
	}
	return found[0][1]
}

// Cancel method cancel the reading process
func (fastaReader *FastaReader) Cancel() {
	if !fastaReader.finished && !fastaReader.cancelled {
		close(fastaReader.done)
		fastaReader.cancelled = true
	}
}

// Alphabet returns Alphabet of the file
func (fastaReader *FastaReader) Alphabet() *seq.Alphabet {
	return fastaReader.t
}
