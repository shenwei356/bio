package fastx

import (
	"bytes"
	"fmt"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/xopen"
)

// Record is a struct for FASTA/Q
type Record struct {
	ID   []byte   // id
	Name []byte   // full name
	Desc []byte   // Description
	Seq  *seq.Seq // seq
}

// Clone of a Record
func (record *Record) Clone() *Record {
	return &Record{
		[]byte(string(record.ID)),
		[]byte(string(record.Name)),
		[]byte(string(record.Desc)),
		record.Seq.Clone(),
	}
}

func (record *Record) String() string {
	return string(record.Format(60))
}

// NewRecord is constructor of type Record for FASTA
func NewRecord(t *seq.Alphabet, id, name, desc, s []byte) (*Record, error) {
	seq, err := seq.NewSeq(t, s)
	if err != nil {
		return nil, fmt.Errorf("error when parsing seq: %s (%s)", id, err)
	}
	return &Record{id, name, desc, seq}, nil
}

// NewRecordWithoutValidation is constructor of type Record for FASTA
// without validation of the sequence
func NewRecordWithoutValidation(t *seq.Alphabet, id, name, desc, s []byte) (*Record, error) {
	seq, err := seq.NewSeqWithoutValidation(t, s)
	if err != nil {
		return nil, err
	}
	return &Record{id, name, desc, seq}, nil
}

// NewRecordWithSeq is constructor of type Record
// for FASTA with a existed seq.Seq object
func NewRecordWithSeq(id, name, desc []byte, s *seq.Seq) (*Record, error) {
	return &Record{id, name, desc, s}, nil
}

// NewRecordWithQual is constructor of type Record for FASTQ
func NewRecordWithQual(t *seq.Alphabet, id, name, desc, s, q []byte) (*Record, error) {
	seq, err := seq.NewSeqWithQual(t, s, q)
	if err != nil {
		return nil, fmt.Errorf("error when parsing seq: %s (%s)", id, err)
	}
	return &Record{id, name, desc, seq}, nil
}

// NewRecordWithQualWithoutValidation is constructor of type Record for FASTQ
func NewRecordWithQualWithoutValidation(t *seq.Alphabet, id, name, desc, s, q []byte) (*Record, error) {
	seq, err := seq.NewSeqWithQualWithoutValidation(t, s, q)
	if err != nil {
		return nil, err
	}
	return &Record{id, name, desc, seq}, nil
}

// ForcelyOutputFastq means outputing record as fastq even if it has no quality (zero-length fastq)
var ForcelyOutputFastq bool

// Format returns formated (wrapped with fixed length of) sequence record
func (record *Record) Format(width int) []byte {
	var buf bytes.Buffer

	if len(record.Seq.Qual) > 0 || ForcelyOutputFastq {
		buf.Write(_mark_fastq)
		buf.Write(record.Name)
		buf.Write(_mark_newline)

		buf.Write(record.Seq.Seq)
		buf.Write(_mark_newline_plus_newline)

		buf.Write(record.Seq.Qual)
		buf.Write(_mark_newline)

		return buf.Bytes()
	}

	buf.Write(_mark_fasta)
	buf.Write(record.Name)
	buf.Write(_mark_newline)

	if width < 1 {
		buf.Write(record.Seq.Seq)
	} else {
		var text []byte
		buffer := poolBuffer.Get().(*bytes.Buffer)
		text, buffer = wrapByteSlice(record.Seq.Seq, width, buffer)
		buf.Write(text)
		poolBuffer.Put(buffer)
	}

	buf.Write(_mark_newline)

	return buf.Bytes()
}

// FormatToWriter formats and directly writes to writer
func (record *Record) FormatToWriter(outfh *xopen.Writer, width int) {
	if len(record.Seq.Qual) > 0 || ForcelyOutputFastq {
		outfh.Write(_mark_fastq)
		outfh.Write(record.Name)
		outfh.Write(_mark_newline)

		outfh.Write(record.Seq.Seq)
		outfh.Write(_mark_newline_plus_newline)

		outfh.Write(record.Seq.Qual)
		outfh.Write(_mark_newline)

		return
	}

	outfh.Write(_mark_fasta)
	outfh.Write(record.Name)
	outfh.Write(_mark_newline)

	if width < 1 {
		outfh.Write(record.Seq.Seq)
	} else {
		var text []byte
		buffer := poolBuffer.Get().(*bytes.Buffer)
		text, buffer = wrapByteSlice(record.Seq.Seq, width, buffer)
		outfh.Write(text)
		poolBuffer.Put(buffer)
	}

	outfh.Write(_mark_newline)
}

// It's unsafe for concurrency
// var buffer *bytes.Buffer

var poolBuffer = &sync.Pool{New: func() interface{} {
	return bytes.NewBuffer(make([]byte, 0, 1024))
}}

var _mark_fasta = []byte{'>'}
var _mark_fastq = []byte{'@'}
var _mark_newline_plus_newline = []byte{'\n', '+', '\n'}
var _mark_newline = []byte{'\n'}

func wrapByteSlice(s []byte, width int, buffer *bytes.Buffer) ([]byte, *bytes.Buffer) {
	if width < 1 {
		return s, buffer
	}
	l := len(s)
	if l == 0 {
		return s, buffer
	}

	var lines int
	if l%width == 0 {
		lines = l/width - 1
	} else {
		lines = int(l / width)
	}

	if buffer == nil {
		buffer = bytes.NewBuffer(make([]byte, 0, l+lines))
	} else {
		buffer.Reset()
	}

	var start, end int
	for i := 0; i <= lines; i++ {
		start = i * width
		end = (i + 1) * width
		if end > l {
			end = l
		}

		buffer.Write(s[start:end])
		if i < lines {
			buffer.Write(_mark_newline)
		}
	}
	return buffer.Bytes(), buffer
}
