package seq

import (
	"fmt"
	"testing"
)

func TestValidateSequence(t *testing.T) {
	dna := []byte("acgt")
	dna2 := []byte("ACGTRYMKSWHBVDN")
	rna := []byte("acgu")
	fake := []byte("acgturymkswhbvdnz")

	ok := DNA.IsValid(dna) == nil &&
		DNAredundant.IsValid(dna2) == nil &&
		RNA.IsValid(rna) == nil &&
		RNA.IsValid(fake) != nil

	if !ok {
		t.Error("validate sequence failed.")
		return
	}
}

func TestRevCom(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtccn-"))
	if string(dna.RevCom().Seq) != "-nggacgt" {
		t.Error("revcom sequence failed.")
		return
	}

	rna, _ := NewSeq(RNA, []byte("auguccn-"))
	if string(rna.RevCom().Seq) != "-nggacau" {
		t.Error("revcom sequence failed.")
		return
	}
}

func TestBaseContent(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtACGT"))
	content := dna.BaseContent("gc")
	wanted := 0.5
	if content != wanted {
		t.Error(fmt.Printf("compution of base content failed: %f != %f", content, wanted))
		return
	}
}

func TestSubSeq(t *testing.T) {
	s, _ := NewSeqWithoutValidate(DNA, []byte("ACGTNacgtn"))
	ok := string(s.SubSeq(1, 1).Seq) == "A" &&
		string(s.SubSeq(2, 4).Seq) == "CGT" &&
		string(s.SubSeq(-4, -2).Seq) == "cgt" &&
		string(s.SubSeq(-1, -1).Seq) == "n" &&
		string(s.SubSeq(2, -2).Seq) == "CGTNacgt" &&
		string(s.SubSeq(1, -1).Seq) == "ACGTNacgtn"

	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	ok = string(s.SubSeq(-4, 2).Seq) == "" &&
		string(s.SubSeq(-3, -4).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}
}
