package fastx

import (
	"fmt"
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestFastaReader2(t *testing.T) {
	file := "test.fa"
	fastxReader, err := NewReader(seq.Unlimit, file, 0, 1, "")
	if err != nil {
		t.Error(t)
	}
	n := 0
	for chunk := range fastxReader.Ch {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		for _, record := range chunk.Data {
			fmt.Print(record)
		}
	}
	fmt.Println(n)
}
