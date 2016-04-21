package fai

import (
	"testing"
)

func TestFastaReader(t *testing.T) {
	file := "seq.fa"
	idx, err := New(file)
	if err != nil {
		t.Error(err)
	}

	s, err := idx.At("cel-let-7", 1)
	if err != nil {
		t.Error(err)
	}
	if s != 'U' {
		t.Errorf("unmatched sequences: cel-let-7")
	}

	seq, err := idx.Get("cel-mir-2", 15, 20)
	if err != nil {
		t.Error(err)
	}
	if seq != "AAAGC" {
		t.Errorf("unmatched sequences cel-mir-2")
	}
}
