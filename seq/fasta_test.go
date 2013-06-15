package seq

import (
	"testing"
)

func TestFastaReader(t *testing.T) {
	NextSeq, err := FastaReader("test.fa")
	if err != nil {
		t.Error(err)
		return
	}

	n := 0
	for {
		_, _, err := NextSeq()
		if err != nil {
			break
		}
		n++
	}
	if n != 3 {
		t.Error("amount of sequences mismatched!")
		return
	}
}
