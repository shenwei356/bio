package seqio

import (
	"testing"

	"github.com/shenwei356/bio/seq"
)

func ReadSeqs(file string) (int, int, error) {
	fasta, err := NewFastaReader(seq.Unlimit, file)
	if err != nil {
		return -1, -1, err
	}

	n := 0
	baseNum := 0
	for fasta.HasNext() {
		record, err := fasta.NextSeq()
		if err != nil {
			return -1, -1, err
		}

		n++
		baseNum += len(record.Seq.Seq)
	}
	return n, baseNum, nil
}

func TestFastaReaderIterator(t *testing.T) {
	seqfile := "../example/hairpin.fa"
	n, baseNum, err := ReadSeqs(seqfile)
	if err != nil {
		t.Error(err)
		return
	}

	fasta, err := NewFastaReader(seq.Unlimit, seqfile)
	if err != nil {
		t.Error(err)
		return
	}

	n2 := 0
	baseNum2 := 0
	for record := range fasta.Iterator(100) {
		n2++
		baseNum2 += len(record.Seq.Seq)
	}

	if n != n2 {
		t.Error("amount of sequences mismatched!")
		return
	}
	if baseNum != baseNum2 {
		t.Error("amount of bases mismatched!")
		return
	}
}
