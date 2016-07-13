package fastx

import (
	"fmt"

	"github.com/brentp/xopen"
	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/byteutil"
)

// Record is a struct for FASTA/Q
type Record struct {
	ID   []byte   // id
	Name []byte   // full name
	Seq  *seq.Seq // seq
}

// Clone of a Record
func (record *Record) Clone() *Record {
	return &Record{
		[]byte(string(record.ID)),
		[]byte(string(record.Name)),
		record.Seq.Clone(),
	}
}

func (record *Record) String() string {
	return string(record.Format(70))
}

// NewRecord is constructor of type Record for FASTA
func NewRecord(t *seq.Alphabet, id, name, s []byte) (*Record, error) {
	seq, err := seq.NewSeq(t, s)
	if err != nil {
		return nil, fmt.Errorf("error when parsing seq: %s (%s)", id, err)
	}
	return &Record{id, name, seq}, nil
}

// NewRecordWithoutValidation is constructor of type Record for FASTA
// without validation of the sequence
func NewRecordWithoutValidation(t *seq.Alphabet, id, name, s []byte) (*Record, error) {
	seq, _ := seq.NewSeqWithoutValidation(t, s)
	return &Record{id, name, seq}, nil
}

// NewRecordWithSeq is constructor of type Record
// for FASTA with a existed seq.Seq object
func NewRecordWithSeq(id, name []byte, s *seq.Seq) (*Record, error) {
	return &Record{id, name, s}, nil
}

// NewRecordWithQual is constructor of type Record for FASTQ
func NewRecordWithQual(t *seq.Alphabet, id, name, s, q []byte) (*Record, error) {
	seq, err := seq.NewSeqWithQual(t, s, q)
	if err != nil {
		return nil, fmt.Errorf("error when parsing seq: %s (%s)", id, err)
	}
	return &Record{id, name, seq}, nil
}

// NewRecordWithQualWithoutValidation is constructor of type Record for FASTQ
func NewRecordWithQualWithoutValidation(t *seq.Alphabet, id, name, s, q []byte) (*Record, error) {
	seq, _ := seq.NewSeqWithQualWithoutValidation(t, s, q)
	return &Record{id, name, seq}, nil
}

// Format returns formated (wrapped with fixed length of) sequence record
func (record *Record) Format(width int) []byte {
	if len(record.Seq.Qual) > 0 {
		return append(append(append(append([]byte(fmt.Sprintf("@%s\n", record.Name)),
			byteutil.WrapByteSlice(record.Seq.Seq, width)...), []byte("\n+\n")...),
			byteutil.WrapByteSlice(record.Seq.Qual, width)...), []byte("\n")...)
	}
	return append(append([]byte(fmt.Sprintf(">%s\n", record.Name)),
		byteutil.WrapByteSlice(record.Seq.Seq, width)...), []byte("\n")...)
}

var defaultBytesBufferSize = 1 << 20

// bufferedByteSliceWrapper is for function FormatToWriter
var bufferedByteSliceWrapper *byteutil.BufferedByteSliceWrapper

func init() {
	bufferedByteSliceWrapper = byteutil.NewBufferedByteSliceWrapper(1, defaultBytesBufferSize)
}

// FormatToWriter formats and directly writes to writer
func (record *Record) FormatToWriter(outfh *xopen.Writer, width int) {
	if len(record.Seq.Qual) > 0 {
		outfh.Write([]byte(fmt.Sprintf("@%s\n", record.Name)))

		// outfh.Write(byteutil.WrapByteSlice(record.Seq.Seq, width))
		text, b := bufferedByteSliceWrapper.Wrap(record.Seq.Seq, width)
		outfh.Write(text)
		outfh.Flush()
		bufferedByteSliceWrapper.Recycle(b)

		outfh.Write([]byte("\n+\n"))

		// outfh.Write(byteutil.WrapByteSlice(record.Seq.Qual, width))
		text, b = bufferedByteSliceWrapper.Wrap(record.Seq.Qual, width)
		outfh.Write(text)
		outfh.Write([]byte("\n"))
		outfh.Flush()
		bufferedByteSliceWrapper.Recycle(b)

		return
	}
	outfh.Write([]byte(fmt.Sprintf(">%s\n", record.Name)))

	// outfh.Write(byteutil.WrapByteSlice(record.Seq.Seq, width))
	text, b := bufferedByteSliceWrapper.Wrap(record.Seq.Seq, width)
	outfh.Write(text)
	outfh.Write([]byte("\n"))
	outfh.Flush()
	bufferedByteSliceWrapper.Recycle(b)
}
