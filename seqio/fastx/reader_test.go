package fastx

import (
	"bytes"
	"io"
	"regexp"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seq"
)

type trackingReadCloser struct {
	*bytes.Reader
	closed bool
}

func (r *trackingReadCloser) Close() error {
	r.closed = true
	return nil
}

func TestReaderKeepsLargeInitialBuffers(t *testing.T) {
	reader, err := NewReaderFromIO(nil, strings.NewReader(">id\nACGT\n"), "")
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	for name, capacity := range map[string]int{
		"record": reader.buffer.Cap(),
		"seq":    reader.seqBuffer.Cap(),
		"qual":   reader.qualBuffer.Cap(),
	} {
		if capacity != defaultBytesBufferSize {
			t.Fatalf("%s buffer capacity = %d, want %d", name, capacity, defaultBytesBufferSize)
		}
	}
}

func TestReaderCloseClosesInput(t *testing.T) {
	source := &trackingReadCloser{Reader: bytes.NewReader([]byte(">id\nACGT\n"))}
	reader, err := NewReaderFromIO(nil, source, "")
	if err != nil {
		t.Fatal(err)
	}

	reader.Close()
	if !source.closed {
		t.Fatal("closing the FASTX reader did not close its input")
	}
}

func TestReaderFromEmptyIO(t *testing.T) {
	source := &trackingReadCloser{Reader: bytes.NewReader(nil)}
	reader, err := NewReaderFromIO(seq.DNA, source, "")
	if err != nil {
		t.Fatal(err)
	}

	if reader.Alphabet() != seq.DNA {
		t.Fatal("empty reader did not retain its configured alphabet")
	}
	if _, err := reader.Read(); err != io.EOF {
		t.Fatalf("empty reader returned %v, want io.EOF", err)
	}
	reader.Close()
	if !source.closed {
		t.Fatal("closing an empty reader did not close its input")
	}
}

func TestReaderRejectsIDRegexpWithoutCapture(t *testing.T) {
	for _, pattern := range []string{`^\S+`, `^\(literal-parentheses\)$`} {
		if _, err := NewReaderFromIO(nil, strings.NewReader(">id\nACGT\n"), pattern); err == nil {
			t.Fatalf("accepted regexp without a capturing group: %q", pattern)
		}
	}

	if got := ParseHeadID(regexp.MustCompile(`^\S+`), []byte("id description")); string(got) != "id description" {
		t.Fatalf("ParseHeadID without capture returned %q", got)
	}
}

func TestReaderHandlesLeadingNewlinesAcrossBuffers(t *testing.T) {
	input := strings.Repeat("\n", bufSize+1) + ">id\nACGT\n"
	reader, err := NewReaderFromIO(nil, strings.NewReader(input), "")
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	record, err := reader.Read()
	if err != nil {
		t.Fatal(err)
	}
	if string(record.ID) != "id" || string(record.Seq.Seq) != "ACGT" {
		t.Fatalf("unexpected record: id=%q seq=%q", record.ID, record.Seq.Seq)
	}
}

func TestReaderGrowsBuffersForLargeRecords(t *testing.T) {
	sequence := strings.Repeat("ACGT", bufSize/2)
	reader, err := NewReaderFromIO(nil, strings.NewReader(">large\n"+sequence+"\n"), "")
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	record, err := reader.Read()
	if err != nil {
		t.Fatal(err)
	}
	if string(record.Seq.Seq) != sequence {
		t.Fatalf("large record length = %d, want %d", len(record.Seq.Seq), len(sequence))
	}
}

func TestFastaReaderCompactsWrappedSequenceInplace(t *testing.T) {
	input := ">id description\r\nAC\r\nG\rT\nA\r"
	reader, err := NewReaderFromIO(nil, strings.NewReader(input), "")
	if err != nil {
		t.Fatal(err)
	}
	defer reader.Close()

	record, err := reader.Read()
	if err != nil {
		t.Fatal(err)
	}
	if got, want := string(record.ID), "id"; got != want {
		t.Fatalf("ID = %q, want %q", got, want)
	}
	if got, want := string(record.Seq.Seq), "ACG\rTA"; got != want {
		t.Fatalf("sequence = %q, want %q", got, want)
	}
	if reader.seqBuffer.Len() != 0 {
		t.Fatalf("FASTA parsing copied sequence into the secondary buffer: %d bytes", reader.seqBuffer.Len())
	}
}

func TestRemoveLineBreaksInplace(t *testing.T) {
	for _, test := range []struct {
		input, want string
	}{
		{"ACGT", "ACGT"},
		{"AC\nGT\n", "ACGT"},
		{"AC\r\nGT\r\n", "ACGT"},
		{"A\rC\nG\r", "A\rCG"},
	} {
		data := []byte(test.input)
		if got := string(removeLineBreaksInplace(data)); got != test.want {
			t.Errorf("removeLineBreaksInplace(%q) = %q, want %q", test.input, got, test.want)
		}
	}
}

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
