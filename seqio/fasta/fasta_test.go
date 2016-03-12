package fasta

import (
	"fmt"
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestFastaReader(t *testing.T) {
	file := "test.fa"
	fastaReader, err := NewFastaReader(seq.Unlimit, file, 1, 0, "")
	if err != nil {
		t.Error(t)
	}
	n := 0
	for chunk := range fastaReader.Ch {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		for _, record := range chunk.Data {
			fmt.Println(record)
		}
	}
	if n != 6 {
		t.Errorf("seq number mismatch %d != %d", 6, n)
	}

}
