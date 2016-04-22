package fai

import (
	"testing"
)

func TestFastaReader(t *testing.T) {
	file := "seq.fa"
	idx, err := New(file)
	checkErr(t, err)
	defer idx.Close()

	chr := "cel-let-7"
	s, err := idx.Base(chr, 1)
	checkErr(t, err)
	if s != 'U' {
		t.Errorf("unmatched sequences %s: %s", chr, []byte{s})
	}

	chr = "cel-mir-2"
	seq, err := idx.Seq(chr)
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
}

func checkErr(t *testing.T, err error) {
	if err != nil {
		t.Error(err)
	}
}
