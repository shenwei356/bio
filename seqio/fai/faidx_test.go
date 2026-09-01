package fai

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/shenwei356/bio/seqio/fastx"
)

func TestFastaReader(t *testing.T) {
	file := "seq.fa"
	idx, err := New(file)
	if err != nil {
		t.Errorf("failed to create faidx for %s: %s", file, err)
		return
	}

	// all sequences
	seqs, err := fastx.GetSeqs(file, nil, 4, 10, fastx.DefaultIDRegexp)
	if err != nil {
		t.Errorf("failed to read seqs: %v", err)
	}
	for _, rec := range seqs {
		seq, err := idx.Seq(string(rec.ID))
		checkErr(t, err)

		if !bytes.Equal(seq, rec.Seq.Seq) {
			t.Errorf("unmatched sequences %s: %s", rec.ID, seq)
		}
	}

	//

	chr := "cel-let-7"
	s, err := idx.Base(chr, 1)
	checkErr(t, err)
	if s != 'U' {
		t.Errorf("unmatched sequences %s: %s", chr, []byte{s})
	}

	chr = "blank"
	seq, err := idx.Seq(chr)
	checkErr(t, err)
	if string(seq) != "" {
		t.Errorf("unmatched sequences %s: %s", chr, seq)
	}

	chr = "cel-mir-2"
	seq, err = idx.Seq(chr)
	checkErr(t, err)
	if string(seq) != "UAAACAGUAUACAGAAAGCCAUCAAAGC" {
		t.Errorf("unmatched sequences %s: %s", chr, seq)
	}

	start, end := 15, 19
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "AAAGC" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	start, end = -3, -1
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "AGC" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	start, end = 1, -1
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "UAAACAGUAUACAGAAAGCCAUCAAAGC" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	start, end = 2, -2
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "AAACAGUAUACAGAAAGCCAUCAAAG" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	start, end = 20, -2
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "CAUCAAAG" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	start, end = 50, -2
	seq, err = idx.SubSeq(chr, start, end)
	checkErr(t, err)
	if string(seq) != "" {
		t.Errorf("unmatched sequences %s from %d to %d: %s", chr, start, end, seq)
	}

	chr = "seq"
	seq, err = idx.Seq(chr)
	checkErr(t, err)
	if string(seq) != "ACTGACTG" {
		t.Errorf("unmatched sequences %s: %s", chr, seq)
	}

	err = idx.Close()
	if err != nil {
		t.Errorf("fail to close faidx: %v", err)
	}
}

func TestFastaReaderNotMapWholeFile(t *testing.T) {
	MapWholeFile = false
	TestFastaReader(t)
	MapWholeFile = true
}

func TestFastaReader2NotMapWholeFile(t *testing.T) {
	MapWholeFile = false
	TestFastaReader2(t)
	MapWholeFile = true
}

func checkErr(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}

func TestFastaReader2(t *testing.T) {
	file := "seq2.fa"
	idx, err := New(file)
	if err != nil {
		t.Errorf("failed to create faidx for %s: %s", file, err)
		return
	}

	seqs, err := fastx.GetSeqs(file, nil, 4, 10, fastx.DefaultIDRegexp)
	if err != nil {
		t.Errorf("failed to read seqs: %v", err)
	}
	for _, rec := range seqs {
		seq, err := idx.Seq(string(rec.ID))
		checkErr(t, err)

		if !bytes.Equal(seq, rec.Seq.Seq) {
			t.Errorf("unmatched sequences %s: %s", rec.ID, seq)
		}
	}

	err = idx.Close()
	if err != nil {
		t.Errorf("fail to close faidx: %v", err)
	}
}

func TestBaseOutOfRangeReturnsError(t *testing.T) {
	idx, err := New("seq.fa")
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()

	if _, err := idx.Base("cel-let-7", 1<<20); !errors.Is(err, ErrPositionOutOfRange) {
		t.Fatalf("Base returned %v, want ErrPositionOutOfRange", err)
	}
}

func TestCloseIsIdempotent(t *testing.T) {
	idx, err := New("seq.fa")
	if err != nil {
		t.Fatal(err)
	}
	if err := idx.Close(); err != nil {
		t.Fatal(err)
	}
	if err := idx.Close(); err != nil {
		t.Fatalf("second Close returned %v", err)
	}
}

func TestFaidxKeepsItsConfiguredAccessMode(t *testing.T) {
	previous := MapWholeFile
	MapWholeFile = false
	idx, err := New("seq.fa")
	MapWholeFile = true
	defer func() { MapWholeFile = previous }()
	if err != nil {
		t.Fatal(err)
	}
	defer idx.Close()

	base, err := idx.Base("cel-let-7", 1)
	if err != nil {
		t.Fatal(err)
	}
	if base != 'U' {
		t.Fatalf("Base returned %q, want U", base)
	}
}

func TestCreateWithIDRegexpDoesNotChangeDefault(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "input.fa")
	if err := os.WriteFile(fasta, []byte(">id description\nACGT\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	custom, err := CreateWithFullHead(fasta, filepath.Join(tmpDir, "custom.fai"))
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := custom["id description"]; !ok {
		t.Fatalf("custom index did not contain the full header: %v", custom)
	}

	standard, err := Create(fasta, filepath.Join(tmpDir, "standard.fai"))
	if err != nil {
		t.Fatal(err)
	}
	if _, ok := standard["id"]; !ok {
		t.Fatalf("custom regexp leaked into default Create: %v", standard)
	}
}

func TestCreateRejectsRegexpWithoutCapture(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "input.fa")
	if err := os.WriteFile(fasta, []byte(">id\nACGT\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	if _, err := CreateWithIDRegexp(fasta, filepath.Join(tmpDir, "input.fai"), `^\S+`); err == nil {
		t.Fatal("accepted ID regexp without a capturing group")
	}
}

func TestCreateRejectsEmptyFASTA(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "empty.fa")
	if err := os.WriteFile(fasta, nil, 0o600); err != nil {
		t.Fatal(err)
	}

	if _, err := Create(fasta, filepath.Join(tmpDir, "empty.fa.fai")); err == nil {
		t.Fatal("created an index for an empty FASTA file")
	}
}

func TestCreateHandlesCRLFAndTrailingEmptyLines(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "input.fa")
	content := []byte(">id description\r\nACGT\r\nAC\r\n\r\n>next\nTT\n")
	if err := os.WriteFile(fasta, content, 0o600); err != nil {
		t.Fatal(err)
	}

	index, err := Create(fasta, filepath.Join(tmpDir, "input.fai"))
	if err != nil {
		t.Fatal(err)
	}
	record := index["id"]
	if record.Length != 6 || record.Start != 17 || record.BasesPerLine != 4 || record.BytesPerLine != 6 {
		t.Fatalf("unexpected CRLF index record: %+v", record)
	}
	next := index["next"]
	if next.Length != 2 || next.Start != 35 || next.BasesPerLine != 2 || next.BytesPerLine != 3 {
		t.Fatalf("unexpected second index record: %+v", next)
	}
}

func TestCreateKeepsBlankRecordLineWidthsCompatible(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "input.fa")
	if err := os.WriteFile(fasta, []byte(">id\nACGT\n>blank\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	index, err := Create(fasta, filepath.Join(tmpDir, "input.fai"))
	if err != nil {
		t.Fatal(err)
	}
	blank := index["blank"]
	if blank.Length != 0 || blank.BasesPerLine != 4 || blank.BytesPerLine != 5 {
		t.Fatalf("unexpected blank index record: %+v", blank)
	}
}

func TestCreatePreservesLineWidthValidation(t *testing.T) {
	for _, test := range []struct {
		name    string
		content string
		valid   bool
	}{
		{"repeated short final width", ">id\nAAAA\nAA\nAA\n", true},
		{"increasing width", ">id\nAA\nAAAA\n", false},
		{"two width transitions", ">id\nAAAA\nAA\nA\n", false},
		{"trailing empty lines", ">id\nAAAA\nAA\n\n\n", true},
	} {
		t.Run(test.name, func(t *testing.T) {
			tmpDir := t.TempDir()
			fasta := filepath.Join(tmpDir, "input.fa")
			if err := os.WriteFile(fasta, []byte(test.content), 0o600); err != nil {
				t.Fatal(err)
			}
			_, err := Create(fasta, filepath.Join(tmpDir, "input.fai"))
			if test.valid && err != nil {
				t.Fatalf("valid line widths returned %v", err)
			}
			if !test.valid && (err == nil || !strings.Contains(err.Error(), "different line length")) {
				t.Fatalf("invalid line widths returned %v", err)
			}
		})
	}
}

func TestCreateHandlesLinesLargerThanReaderBuffer(t *testing.T) {
	tmpDir := t.TempDir()
	fasta := filepath.Join(tmpDir, "input.fa")
	description := strings.Repeat("x", 70<<10)
	sequence := strings.Repeat("ACGT", 20<<10)
	header := ">id " + description + "\n"
	if err := os.WriteFile(fasta, []byte(header+sequence+"\nACGT\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	index, err := Create(fasta, filepath.Join(tmpDir, "input.fai"))
	if err != nil {
		t.Fatal(err)
	}
	record := index["id"]
	if record.Length != len(sequence)+4 || record.Start != int64(len(header)) || record.BasesPerLine != len(sequence) || record.BytesPerLine != len(sequence)+1 {
		t.Fatalf("unexpected long-line index record: %+v", record)
	}
}
