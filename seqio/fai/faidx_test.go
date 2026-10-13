package fai

import (
	"bytes"
	"errors"
	"os"
	"path/filepath"
	"reflect"
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

// TestFaidxEmptyAndLeadingInput checks creation, index reloading, and sequence access.
func TestFaidxEmptyAndLeadingInput(t *testing.T) {
	previous := MapWholeFile
	defer func() { MapWholeFile = previous }()
	for _, test := range []struct {
		name, input, id, sequence string
		count                     int
		invalid                   bool
	}{
		{name: "zero bytes"},
		{name: "one newline", input: "\n"},
		{name: "multiple newlines", input: "\n\n\n"},
		{name: "newlines across buffers", input: strings.Repeat("\n", faiReaderSize+1)},
		{name: "space", input: " ", invalid: true},
		{name: "space and newline", input: " \n", invalid: true},
		{name: "tab", input: "\t", invalid: true},
		{name: "tab and newline", input: "\t\n", invalid: true},
		{name: "CRLF", input: "\r\n", invalid: true},
		{name: "carriage return", input: "\r", invalid: true},
		{name: "space after newlines", input: "\n\n \n", invalid: true},
		{name: "space across buffers", input: strings.Repeat("\n", faiReaderSize) + " ", invalid: true},
		{name: "sequence without header", input: "ACGT\n", invalid: true},
		{name: "BOM only", input: "\xef\xbb\xbf"},
		{name: "BOM and newlines", input: "\xef\xbb\xbf\n\n"},
		{name: "BOM and bare header", input: "\xef\xbb\xbf>", count: 1},
		{name: "BOM", input: "\xef\xbb\xbf>id\nACGT\n", id: "id", sequence: "ACGT", count: 1},
		{name: "BOM after newline", input: "\n\xef\xbb\xbf>id\nACGT\n", invalid: true},
		{name: "indented header", input: " >id\nACGT\n", invalid: true},
		{name: "CRLF before header", input: "\r\n>id\nACGT\n", invalid: true},
		{name: "FASTQ", input: "@id\nACGT\n+\n!!!!\n", invalid: true},
		{name: "FASTQ after newlines", input: "\n\n@id\nACGT\n+\n!!!!\n", invalid: true},
		{name: "bare header", input: ">", count: 1},
		{name: "empty record", input: ">\n", count: 1},
		{name: "duplicate empty IDs", input: ">\n>\n", count: 1},
		{name: "empty record after newlines", input: "\n\n>", count: 1},
		{name: "header without final newline", input: ">id", id: "id", count: 1},
		{name: "last header without final newline", input: ">first\nAA\n>id", id: "id", count: 2},
		{name: "empty sequence", input: ">id\n", id: "id", count: 1},
		{name: "empty ID", input: ">\nACGT\n", sequence: "ACGT", count: 1},
		{name: "leading newlines", input: "\n\n>id\nACGT\n", id: "id", sequence: "ACGT", count: 1},
		{name: "header across buffers", input: strings.Repeat("\n", faiReaderSize+1) + ">id\nACGT\n", id: "id", sequence: "ACGT", count: 1},
		{name: "sequence without final newline", input: ">id\nACGT", id: "id", sequence: "ACGT", count: 1},
		{name: "sequence ending in CR", input: ">id\nACGT\r", id: "id", sequence: "ACGT", count: 1},
		{name: "wrapped sequence without final newline", input: ">id\nAC\nGT", id: "id", sequence: "ACGT", count: 1},
		{name: "trailing empty lines", input: ">id\nACGT\n\n\n", id: "id", sequence: "ACGT", count: 1},
		{name: "CRLF record", input: ">id\r\nACGT\r\n", id: "id", sequence: "ACGT", count: 1},
	} {
		t.Run(test.name, func(t *testing.T) {
			for _, mapWholeFile := range []bool{true, false} {
				MapWholeFile = mapWholeFile
				fasta := filepath.Join(t.TempDir(), "input.fa")
				if err := os.WriteFile(fasta, []byte(test.input), 0o600); err != nil {
					t.Fatal(err)
				}
				// The second New must behave the same when an index already exists.
				for attempt := 0; attempt < 2; attempt++ {
					idx, err := New(fasta)
					if test.invalid {
						if err == nil {
							idx.Close()
							t.Fatal("accepted invalid FASTA input")
						}
						if _, err := os.Stat(fasta + ".fai"); !os.IsNotExist(err) {
							t.Fatalf("failed creation left an index: %v", err)
						}
						continue
					}
					if err != nil {
						t.Fatal(err)
					}
					func() {
						defer idx.Close()
						if len(idx.Index) != test.count {
							t.Fatalf("got %d records, want %d", len(idx.Index), test.count)
						}
						reloaded, err := Read(fasta + ".fai")
						if err != nil {
							t.Fatal(err)
						}
						if !reflect.DeepEqual(idx.Index, reloaded) {
							t.Fatalf("reloaded index differs: %v != %v", reloaded, idx.Index)
						}
						if test.count == 0 {
							if _, err := idx.Seq("missing"); !errors.Is(err, ErrSeqNotExists) {
								t.Fatalf("empty index query returned %v", err)
							}
							return
						}
						sequence, err := idx.Seq(test.id)
						if err != nil || string(sequence) != test.sequence {
							t.Fatalf("got sequence=%q error=%v, want %q", sequence, err, test.sequence)
						}
						if idx.Index[test.id].Length != len(test.sequence) {
							t.Fatalf("incorrect sequence length: %+v", idx.Index[test.id])
						}
						if len(test.sequence) > 0 {
							base, err := idx.Base(test.id, 1)
							if err != nil || base != test.sequence[0] {
								t.Fatalf("got first base=%q error=%v", base, err)
							}
						}
					}()
				}
			}
		})
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
		{"invalid width after a complete record", ">first\nAA\n>id\nAA\nAAAA\n", false},
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
			if !test.valid {
				if _, err := os.Stat(filepath.Join(tmpDir, "input.fai")); !os.IsNotExist(err) {
					t.Fatalf("invalid line widths left an incomplete index: %v", err)
				}
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
