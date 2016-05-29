package fastx

import (
	"fmt"
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestFastaReader(t *testing.T) {
	file := "test.fa"
	fastxReader, err := NewReader(seq.Unlimit, file, 1, 0, "")
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
			fmt.Println(record)
		}
	}
	if n != 6 {
		t.Errorf("seq number mismatch %d != %d", 6, n)
	}
}

func TestFastqReader(t *testing.T) {
	file := "test.fq"
	fastxReader, err := NewReader(seq.Unlimit, file, 1, 0, "")
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
			fmt.Println(record)
		}
	}
	if n != 5 {
		t.Errorf("seq number mismatch %d != %d", 5, n)
	}
}
func TestFastqReader2(t *testing.T) {
	file := "test2.fq"
	fastxReader, err := NewReader(seq.Unlimit, file, 1, 0, "")
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
			fmt.Println(record)
		}
	}
	if n != 5 {
		t.Errorf("seq number mismatch %d != %d", 5, n)
	}
}

func TestFastqReader3(t *testing.T) {
	file := "test3.fq"
	fastxReader, err := NewReader(seq.Unlimit, file, 1, 0, "")
	if err != nil {
		t.Error(t)
	}
	n := 0
	l := -1
	for chunk := range fastxReader.Ch {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		for _, record := range chunk.Data {
			fmt.Println(record)
			if l == -1 {
				l = len(record.Seq.Seq)
			} else {
				if l != len(record.Seq.Seq) {
					t.Errorf("parse error")
					return
				}
			}
		}
	}
	if n != 3 {
		t.Errorf("seq number mismatch %d != %d", 5, n)
	}
}

func TestBlankFile(t *testing.T) {
	file := "blank.fx"
	fastxReader, err := NewReader(seq.Unlimit, file, 1, 0, "")
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
			fmt.Println(record)
		}
	}
	if n != 0 {
		t.Errorf("seq number mismatch %d != %d", 0, n)
	}
}
