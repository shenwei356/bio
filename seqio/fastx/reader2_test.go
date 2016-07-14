package fastx

import (
	"fmt"
	"io"
	"testing"
)

func TestFastaReader2(t *testing.T) {
	file := "test.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(t)
	}
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Error(err)
			break
		}
		fmt.Print(record)
	}
}

func TestFastaReader3(t *testing.T) {
	file := "test.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(t)
	}

	ch := reader.ChunkChan(0, 1)
	for chunk := range ch {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		for _, record := range chunk.Data {
			fmt.Print(record)
		}
	}
}
