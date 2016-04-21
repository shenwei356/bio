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
	defer idx.Close()

	s, err := idx.Base("cel-let-7", 1)
	if err != nil {
		t.Error(err)
	}
	if s != 'U' {
		t.Errorf("unmatched sequences: cel-let-7")
	}

	seq, err := idx.SubSeq("cel-mir-2", 15, 19)
	if err != nil {
		t.Error(err)
	}
	if seq != "AAAGC" {
		t.Errorf("unmatched sequences cel-mir-2")
	}

	seq, err = idx.Seq("cel-mir-2")
	if err != nil {
		t.Error(err)
	}
	if seq != "UAAACAGUAUACAGAAAGCCAUCAAAGC" {
		t.Errorf("unmatched sequences cel-mir-2: %s", seq)
	}
}
