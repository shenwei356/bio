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
	s, _ := NewSeqWithoutValidation(DNA, []byte("ACGTNacgtn"))
	ok := string(s.SubSeq(1, 1).Seq) == "A" &&
		string(s.SubSeq(2, 4).Seq) == "CGT" &&
		string(s.SubSeq(-4, -2).Seq) == "cgt" &&
		string(s.SubSeq(-4, -1).Seq) == "cgtn" &&
		string(s.SubSeq(-1, -1).Seq) == "n" &&
		string(s.SubSeq(2, -2).Seq) == "CGTNacgt" &&
		string(s.SubSeq(1, -1).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(12, 14).Seq) == "" &&
		string(s.SubSeq(-10, -1).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(-10, -3).Seq) == "ACGTNacg" &&
		string(s.SubSeq(1, 10).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(3, 12).Seq) == "GTNacgtn" &&
		string(s.SubSeq(3, 100).Seq) == "GTNacgtn"

	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	ok = string(s.SubSeq(-4, 2).Seq) == "" &&
		string(s.SubSeq(-3, -4).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	s, _ = NewSeqWithoutValidation(DNA, []byte(""))
	ok = string(s.SubSeq(1, 4).Seq) == "" &&
		string(s.SubSeq(2, 4).Seq) == "" &&
		string(s.SubSeq(1, -1).Seq) == "" &&
		string(s.SubSeq(-4, -1).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}
}

func TestSubSeqInplace(t *testing.T) {
	s, _ := NewSeqWithoutValidation(DNA, []byte("ACGTNacgtn"))
	ok := string(s.Clone().SubSeqInplace(1, 1).Seq) == "A" &&
		string(s.Clone().SubSeqInplace(2, 4).Seq) == "CGT" &&
		string(s.Clone().SubSeqInplace(-4, -2).Seq) == "cgt" &&
		string(s.Clone().SubSeqInplace(-4, -1).Seq) == "cgtn" &&
		string(s.Clone().SubSeqInplace(-1, -1).Seq) == "n" &&
		string(s.Clone().SubSeqInplace(2, -2).Seq) == "CGTNacgt" &&
		string(s.Clone().SubSeqInplace(1, -1).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(-10, -1).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(-10, -3).Seq) == "ACGTNacg" &&
		string(s.Clone().SubSeqInplace(1, 10).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(3, 10).Seq) == "GTNacgtn" &&
		string(s.Clone().SubSeqInplace(3, 100).Seq) == "GTNacgtn"

	if !ok {
		t.Error(fmt.Printf("SubSeqInplace error"))
	}

	ok = string(s.Clone().SubSeqInplace(-4, 2).Seq) == "" &&
		string(s.Clone().SubSeqInplace(-3, -4).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("SubSeqInplace error"))
	}

	s, _ = NewSeqWithoutValidation(DNA, []byte(""))
	ok = string(s.Clone().SubSeqInplace(1, 4).Seq) == "" &&
		string(s.Clone().SubSeqInplace(2, 4).Seq) == "" &&
		string(s.Clone().SubSeqInplace(1, -1).Seq) == "" &&
		string(s.Clone().SubSeqInplace(-4, -1).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}
}
