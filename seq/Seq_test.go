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
		t.Error("revcom sequence failed.")
		return
	}

	rna, _ := NewSeq(RNA, []byte("auguccn-"))
	if string(rna.Revcom()) != "-nggacau" {
		t.Error("revcom sequence failed.")
		return
	}
}

func TestBaseContent(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtACGT"))
	content := dna.BaseContent([]byte("gc"))
	wanted := 0.5
	if content != wanted {
		t.Error(fmt.Printf("compution of base content failed: %f != %f", content, wanted))
		return
	}
}
