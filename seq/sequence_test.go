package seq

import (
	"testing"
)

func TestValidateSequence(t *testing.T) {
	dna := "acgt"
	dna2 := "ACGTURYMKSWHBVDN"
	rna := "acgu"
	fake := "acgturymkswhbvdnz"

	ok := ValidateSequence(dna) && ValidateSequence(dna2) &&
		ValidateSequence(rna) && !ValidateSequence(fake)
	if !ok {
		t.Error("Validate sequence failed.")
		return
	}
}

func TestRevcom(t *testing.T) {
	dna := "acgt"
	rna := "acgu"

	ok := Revcom(dna, DNA) == "ACGT" && Revcom(rna, RNA) == "ACGU"
	if !ok {
		t.Error("Revcom sequence failed.")
		return
	}
}

func TestBaseContent(t *testing.T) {
	dna := "acgtACGT"

	ok := BaseContent(dna, "a") == 0.25 && BaseContent(dna, "c") == 0.25 &&
		BaseContent(dna, "g") == 0.25 && BaseContent(dna, "t") == 0.25 &&
		BaseContent(dna, "u") == 0
	if !ok {
		t.Error("Compution of base content failed.")
		return
	}
}
