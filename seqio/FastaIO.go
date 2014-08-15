package seqio

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"regexp"

	"github.com/shenwei356/bio/seq"
)

// Type *FastaRecord*
type FastaRecord struct {
	Id  []byte
	Seq *seq.Seq
}

// Constructor of type *FastaRecord*
func NewFastaRecord(t *seq.Alphabet, id, str []byte) (*FastaRecord, error) {
	sequence, err := seq.NewSeq(t, str)
	if err != nil {
		return nil, errors.New(fmt.Sprintf("%s: %s", err, id))
	}
	return &FastaRecord{id, sequence}, nil
}

// Format output
func (record *FastaRecord) FormatSeq(width int) []byte {
	return seq.FormatSeq(record.Seq.Seq, width)
}

// Fasta Writter
type FastaWriter struct {
	filename  string
	lineWidth int
}

// Constructor of type *FastaWriter*
func NewFastaWriter(filename string, lineWidth int) *FastaWriter {
	return &FastaWriter{filename, lineWidth}
}

// Write list of FastaRecords to file
func (writer FastaWriter) Write(records []*FastaRecord) (int, error) {
	fout, err := os.Create(writer.filename)
	defer fout.Close()
	if err != nil {
		return 0, err
	}

	w := bufio.NewWriter(fout)
	n := 0
	for _, record := range records {
		fmt.Fprintf(w, ">%s\n%s", record.Id, record.FormatSeq(writer.lineWidth))
		n++
	}
	w.Flush()
	return n, nil
}

/* FastaReader. Usage:

   // Sequence type must be specified.
   fasta, err := seqio.NewFastaReader(seq.RNAredundant, "hairpin.fa")

   // check err, usually caused by a wrong file path.
   if err != nil {
       fmt.Println(err)
       return
   }

   // read and check if more record existed
   for fasta.HasNext() {
        record, err := fasta.NextSeq()
        // check err, the record may contain invalid sequence !!!
        if err != nil {
            fmt.Println(err)
            continue
        }

        // deal with the fasta record
        // format output
        fmt.Printf(">%s\n%s", record.Id, record.FormatSeq(70))
   }

*/
type FastaReader struct {
	t        *seq.Alphabet
	filename string
	nextseq  *FastaRecord

	fh                *os.File
	reader            *bufio.Reader
	buffer            bytes.Buffer
	secondLastHead    []byte // report when NextSeq() return a invalid record
	lastHead          []byte
	hasSeq            bool
	fileHandlerClosed bool
}

// Constructor of FastaReader. Sequence Type (Alphabet) must be given,
// to validate the sequence.
func NewFastaReader(t *seq.Alphabet, filename string) (*FastaReader, error) {
	fh, err := os.Open(filename)
	if err != nil {
		return nil, errors.New("failed to open file: " + filename)
	}

	fasta := new(FastaReader)
	fasta.t = t
	fasta.filename = filename
	fasta.fh = fh

	fasta.reader = bufio.NewReader(fh)
	fasta.hasSeq = false
	fasta.fileHandlerClosed = false

	return fasta, nil
}

// return the next record
func (fasta *FastaReader) NextSeq() (*FastaRecord, error) {
	if fasta.nextseq == nil {
		return nil, errors.New("invalid " + fasta.t.Type() +
			" sequence: " + string(fasta.secondLastHead))
	}
	return fasta.nextseq, nil
}

// parsing fasta file to check if more record existed,
// if existed, store it.
func (fasta *FastaReader) HasNext() bool {
	if fasta.fileHandlerClosed {
		return false
	}

	trimLeftSpace := regexp.MustCompile(`^\s+`)
	trimSpace := regexp.MustCompile(`[\r\n\s]+`)
	var seq, head []byte
	for {
		str, err := fasta.reader.ReadBytes('\n')
		str = trimLeftSpace.ReplaceAll(str, []byte("")) // spaces before > are allowed

		if err == io.EOF {
			// do not forget the last line,
			// even which does not ends with "\n"
			fasta.buffer.Write(trimSpace.ReplaceAll(str, []byte("")))

			seq = fasta.buffer.Bytes()
			fasta.buffer.Reset()
			fasta.secondLastHead = fasta.lastHead
			head = fasta.lastHead
			fasta.fh.Close()
			fasta.fileHandlerClosed = true
			err = nil

			fasta.nextseq, _ = NewFastaRecord(fasta.t, head, seq)
			return true
		}

		if bytes.HasPrefix(str, []byte(">")) {
			fasta.hasSeq = true

			thisHead := trimSpace.ReplaceAll(str[1:], []byte(""))
			if fasta.buffer.Len() > 0 { // no-first seq head
				seq = fasta.buffer.Bytes()
				fasta.buffer.Reset()
				fasta.secondLastHead = fasta.lastHead
				head, fasta.lastHead = fasta.lastHead, thisHead

				fasta.nextseq, _ = NewFastaRecord(fasta.t, head, seq)
				return true
			} else { // first sequence head
				fasta.secondLastHead = thisHead
				fasta.lastHead = thisHead
			}
		} else if fasta.hasSeq { // append sequence
			fasta.buffer.Write(trimSpace.ReplaceAll(str, []byte("")))
		} else {
			// some line before the first ">"
		}
	}
}

// Return an fasta record iterator.
// It's Go-ish by using channel!
//
// The arguments buffersize specifies the buffer size of channel.
// If buffersize == 1, it reads one fasta record and waits the
// record being deliveried.
// When buffersize > 1, records are well prepared in the backgound.
// It's useful when reading large sequences (human chromosome) and
// processing of sequence is also time-consuming.
//
// Usage:
//
//		fasta, err := NewFastaReader(seq.Unlimit, "test.fa")
//		if err != nil {
//			t.Error(err)
//			return
//		}
//
//		for record = range fasta.Iterator(10) {
//			fmt.Printf(">%s\n%s", record.Id, record.FormatSeq(70))
//		}
//
//
func (fasta *FastaReader) Iterator(buffersize int) chan *FastaRecord {
	if buffersize < 1 {
		buffersize = 1
	}

	ch := make(chan *FastaRecord, buffersize)
	go func() {
		for fasta.HasNext() {
			record, err := fasta.NextSeq()
			if err != nil {
				fmt.Fprint(os.Stderr, err)
				continue
			}
			ch <- record
		}
		close(ch)
	}()
	return ch
}
