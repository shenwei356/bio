package seqio

import (
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestFastaReader(t *testing.T) {
	fasta, err := NewFastaReader(seq.DNA, "test.fa")
	if err != nil {
		t.Error(err)
		return
	}

	n := 0
	for fasta.HasNext() {
		fasta.NextSeq()
		n++
	}
	if n != 3 {
		t.Error("amount of sequences mismatched!")
		return
	}
}
