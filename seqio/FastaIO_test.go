package seqio

import (
	"fmt"
	"os"
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestFastaReader(t *testing.T) {
	fasta, err := NewFastaReader(seq.Unlimit, "test.fa")
	if err != nil {
		t.Error(err)
		return
	}

	n := 0
	for fasta.HasNext() {
		_, err := fasta.NextSeq()
		if err != nil {
			fmt.Fprintln(os.Stderr, err)
			continue
		}

		n++
	}
	if n != 3 {
		t.Error("amount of sequences mismatched!")
		return
	}
}

func TestFastaReaderIterator(t *testing.T) {
	fasta, err := NewFastaReader(seq.Unlimit, "test.fa")
	if err != nil {
		t.Error(err)
		return
	}

	n := 0
	for _ = range fasta.Iterator(2) {

		n++
	}

	if n != 3 {
		t.Error("amount of sequences mismatched!")
		return
	}
}
