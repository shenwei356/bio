package fai

import (
	"fmt"
	"io"
	"os"

	"github.com/edsrzf/mmap-go"
)

// Faidx is
type Faidx struct {
	reader io.ReadSeeker
	Index  Index
	mmap   mmap.MMap
}

// New try to get Faidx from fasta file
func New(fileSeq string) (*Faidx, error) {
	fileFai := fileSeq + ".fai"
	return NewWithCustomExt(fileSeq, fileFai)
}

// NewWithCustomExt try to get Faidx from fasta file, and .fai is specified
func NewWithCustomExt(fileSeq, fileFai string) (*Faidx, error) {
	var index Index
	if _, err := os.Stat(fileFai); os.IsNotExist(err) {
		index, err = Create(fileSeq, fileFai)
		if err != nil {
			return nil, err
		}
	} else {
		index, err = Read(fileFai)
		if err != nil {
			return nil, err
		}
	}

	return NewWithIndex(fileSeq, index)
}

// NewWithIndex return faidx from file and readed Index.
// Useful for using custom IDRegexp
func NewWithIndex(file string, index Index) (*Faidx, error) {
	reader, err := os.Open(file)
	if err != nil {
		return nil, err
	}

	m, err := mmap.Map(reader, mmap.RDONLY, 0)
	if err != nil {
		return nil, err
	}

	return &Faidx{reader, index, m}, nil
}

// p is 0-based
func position(r Record, p int) int64 {
	if p < 0 {
		p = 0
	}
	if p > r.Length {
		p = r.Length
	}
	return r.Start + int64(p/r.BasesPerLine*r.BytesPerLine+p%r.BasesPerLine)
}

// ErrSeqNotExists means that sequence not exists
var ErrSeqNotExists = fmt.Errorf("sequence not exists")

// Seq returns sequence of chr
func (f *Faidx) Seq(chr string) ([]byte, error) {
	index, ok := f.Index[chr]
	if !ok {
		return nil, ErrSeqNotExists
	}

	pstart := position(index, 0)
	pend := position(index, index.Length)
	return []byte(string(cleanSeq(f.mmap[pstart:pend]))), nil
}

// SubSeq returns subsequence of chr from start to end. start and end are 1-based.
func (f *Faidx) SubSeq(chr string, start int, end int) ([]byte, error) {
	index, ok := f.Index[chr]
	if !ok {
		return nil, ErrSeqNotExists
	}

	start, end, ok = SubLocation(index.Length, start, end)
	if !ok {
		return []byte{}, nil
	}

	pstart := position(index, start-1)
	pend := position(index, end)
	return []byte(string(cleanSeq(f.mmap[pstart:pend]))), nil
}

// Base returns base in postion pos. pos is 1 based
func (f *Faidx) Base(chr string, pos int) (byte, error) {
	index, ok := f.Index[chr]
	if !ok {
		return '*', ErrSeqNotExists
	}

	ppos := position(index, pos-1)
	return f.mmap[ppos], nil
}

// Close the readers
func (f *Faidx) Close() {
	f.reader.(io.Closer).Close()
	f.mmap.Unmap()
}

/*SubLocation is my sublocation strategy,
start, end and returned start and end are all 1-based

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

*/
func SubLocation(length, start, end int) (int, int, bool) {
	if start < 1 {
		if start == 0 {
			start = 1
		} else if start < 0 {
			if end < 0 && start > end {
				return start, end, false
			}
		}
		start = length + start + 1
	}
	if start > length {
		return start, end, false
	}

	if end > length {
		end = length
	}
	if end < 1 {
		if end == 0 {
			end = -1
		}
		end = length + end + 1
	}

	if start-1 > end {
		return start - 1, end, false
	}
	return start, end, true
}
