package fai

import (
	"bytes"
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
