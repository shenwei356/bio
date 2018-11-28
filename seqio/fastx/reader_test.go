package fastx

import (
	// "fmt"

	"io"
	"testing"
)

func TestFastaReader2(t *testing.T) {
	file := "test.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	for {
		_, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Error(err)
			break
		}
		// fmt.Print(record)
	}
}

func TestFastaReader3(t *testing.T) {
	file := "test.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}

	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		// for _, record := range chunk.Data {
		// 	// fmt.Print(record)
		// }
	}
}

func TestFastqReadern(t *testing.T) {
	file := "test.fq"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}

	for {
		_, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			t.Error(err)
			break
		}
		// fmt.Print(record)
	}
}

// -----------------------------

func TestFastaReader(t *testing.T) {
	file := "test.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		// for _, record := range chunk.Data {
		// fmt.Println(record)
		// }
	}
	if n != 6 {
		t.Errorf("seq number mismatch %d != %d", 6, n)
	}
}

func TestFastqReader(t *testing.T) {
	file := "test.fq"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		// for _, record := range chunk.Data {
		// fmt.Println(record)
		// }
	}
	if n != 8 {
		t.Errorf("seq number mismatch %d != %d", 8, n)
	}
}
func TestFastqReader2(t *testing.T) {
	file := "test2.fq"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		// for _, record := range chunk.Data {
		// fmt.Println(record)
		// }
	}
	if n != 5 {
		t.Errorf("seq number mismatch %d != %d", 5, n)
	}
}

func TestFastqReader3(t *testing.T) {
	file := "test3.fq"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	l := -1
	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			t.Error(chunk.Err)
		}
		n += len(chunk.Data)
		for _, record := range chunk.Data {
			// fmt.Println(record)
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
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	for chunk := range reader.ChunkChan(0, 1) {
		if chunk.Err != nil {
			if chunk.Err != ErrNotFASTXFormat {
				t.Error(chunk.Err)
			}
		}
	}
}

func TestBlankFile2(t *testing.T) {
	file := "blank1.fx"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	for chunk := range reader.ChunkChan(0, 1) {
		// should not reach here
		t.Errorf("should not reach here. error: %s", chunk.Err)
		return
	}
}
func TestEmptyFile(t *testing.T) {
	file := "empty.fx"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	for chunk := range reader.ChunkChan(0, 1) {
		// should not reach here
		t.Errorf("should not reach here. error: %s", chunk.Err)
		return
	}
}
