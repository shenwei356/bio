package seq

import (
	"testing"
)

func TestValidateSequence(t *testing.T) {
	dna := []byte("acgt")
	dna2 := []byte("ACGTRYMKSWHBVDN")
	rna := []byte("acgu")
	fake := []byte("acgturymkswhbvdnz")

	ok := DNA.IsValid(dna) && DNAredundant.IsValid(dna2) &&
		RNA.IsValid(rna) && !RNA.IsValid(fake)

	if !ok {
		t.Error("validate sequence failed.")
		return
	}
}

func TestRevcom(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtccn-"))
	if string(dna.Revcom()) != "-nggacgt" {
		t.Error("Revcom sequence failed.")
		return
	}

	rna, _ := NewSeq(RNA, []byte("auguccn-"))
	if string(rna.Revcom()) != "-nggacau" {
		t.Error("Revcom sequence failed.")
		return
	}
}

func TestBaseContent(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtACGT"))

	if dna.BaseContent([]byte("gc")) == 0.50 {
		t.Error("Compution of base content failed.")
		return
	}
}
