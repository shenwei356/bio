package fastx

import (
	"bytes"
	"io"
	"testing"
)

func TestReaderIDRegexpStateIsPerReader(t *testing.T) {
	customReader, err := NewReaderFromIO(nil, bytes.NewBufferString(">id description\nACGT\n"), `^(.+)$`)
	if err != nil {
		t.Fatal(err)
	}
	defer customReader.Close()

	defaultReader, err := NewReaderFromIO(nil, bytes.NewBufferString(">other description\nACGT\n"), "")
	if err != nil {
		t.Fatal(err)
	}
	defer defaultReader.Close()

	defaultRecord, err := defaultReader.Read()
	if err != nil {
		t.Fatal(err)
	}
	if got, want := string(defaultRecord.ID), "other"; got != want {
		t.Fatalf("default ID parsing: got %q, want %q", got, want)
	}

	record, err := customReader.Read()
	if err != nil {
		t.Fatal(err)
	}
	if got, want := string(record.ID), "id description"; got != want {
		t.Fatalf("custom ID regexp affected by another reader: got %q, want %q", got, want)
	}
}

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
		// 	fmt.Println(record)
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
	file := "blank2.fx"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	for chunk := range reader.ChunkChan(0, 1) {
		// should not reach here
		// r := chunk.Data[0]
		// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)

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

func TestEmptyID(t *testing.T) {
	file := "empty_id.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		n++

		r := chunk.Data[0]
		// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)
		if !(len(r.ID) == 0 && len(r.Seq.Seq) > 0) {
			t.Errorf("sequence id should be empty and the seq should be non-empty")
		}
	}

	if n != 1 {
		t.Errorf("exact one record should be returned")
	}
}

func TestEmptySeq(t *testing.T) {
	file := "empty_seq.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		n++

		r := chunk.Data[0]
		// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)
		if !(len(r.ID) > 0 && len(r.Seq.Seq) == 0) {
			t.Errorf("sequence id should be non-empty and the seq should be empty")
		}
	}

	if n != 1 {
		t.Errorf("exact one record should be returned")
	}
}

func TestEmptySeq2(t *testing.T) {
	file := "empty_seq2.fa"
	reader, err := NewDefaultReader(file)
	if err != nil {
		t.Error(err)
	}
	n := 0
	for chunk := range reader.ChunkChan(0, 1) {
		n++

		r := chunk.Data[0]
		// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)
		if !(len(r.ID) > 0 && len(r.Seq.Seq) == 0) {
			t.Errorf("sequence id should be non-empty and the seq should be empty")
		}
	}

	if n != 2 {
		t.Errorf("exact one record should be returned")
	}
}

func TestEmptyIDAndSeq(t *testing.T) {
	files := []string{"empty_id_and_seq.fa", "empty_id_and_seq.fa"}
	for _, file := range files {
		reader, err := NewDefaultReader(file)
		if err != nil {
			t.Error(err)
		}
		n := 0
		for chunk := range reader.ChunkChan(0, 1) {
			n++

			r := chunk.Data[0]
			// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)
			if !(len(r.ID) == 0 && len(r.Seq.Seq) == 0) {
				t.Errorf("file: %s, both sequence id and seq should be empty", file)
			}
		}

		if n != 1 {
			t.Errorf("file: %s, exact one record should be returned", file)
		}
	}
}

func TestEmptyIDAndSeq2(t *testing.T) {
	files := []string{"empty_id_and_seq_one_record.fa", "empty_id_and_seq_one_record2.fa"}
	for _, file := range files {
		reader, err := NewDefaultReader(file)
		if err != nil {
			t.Error(err)
		}
		n := 0
		for chunk := range reader.ChunkChan(0, 1) {
			n++

			if n == 2 {
				r := chunk.Data[0]
				// fmt.Printf("id:(%s)\nseq:(%s)\n", r.ID, r.Seq.Seq)
				if !(len(r.ID) == 0 && len(r.Seq.Seq) == 0) {
					t.Errorf("file: %s, both sequence id and seq should be empty", file)
				}
			}
		}

		if n != 3 {
			t.Errorf("file: %s, exact one record should be returned", file)
		}
	}
}
