// Copyright © 2018
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
package seq

import (
	"bytes"
	"regexp"
	"testing"
)

func TestDNAToProteinInvalidInputs(t *testing.T) {
	dna := []byte("ATGGGCAAAGGCAAAGGC")
	seq, err := NewSeq(DNA, dna)
	if err != nil {
		t.Error("creating DNA sequence failed.")
		return
	}
	transl_table := 99
	_, err = DNAToProtein(seq, transl_table)
	if err == nil {
		t.Error("invalid transl table number should have been detected")
		return
	}
}

func testDNAToProteinHelper(t *testing.T, dna []byte, protein []byte, transl_table int) {
	dnaValid := DNA.IsValid(dna) == nil
	if !dnaValid {
		t.Error("validating DNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	seq, err := NewSeq(DNA, dna)
	if err != nil {
		t.Error("creating DNA sequence failed.")
		return
	}
	translatedProtein, err := DNAToProtein(seq, transl_table)
	if err != nil {
		t.Error("error converting DNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting DNA to protein doesnt match")
		return
	}
}
func TestDNAToProteinTransl1(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	testDNAToProteinHelper(t, dna, protein, 1)
}

func TestDNAToProteinTransl2(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`**MW`, ""))
	testDNAToProteinHelper(t, dna, protein, 2)
}

func TestDNAToProteinTransl3(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`ATACTTCTCCTACTGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`MTTTTW`, ""))
	testDNAToProteinHelper(t, dna, protein, 3)
}

func TestDNAToProteinTransl4(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACCTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FASTW`, ""))
	testDNAToProteinHelper(t, dna, protein, 4)
}

func TestDNAToProteinTransl5(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SSMW`, ""))
	testDNAToProteinHelper(t, dna, protein, 5)
}

func TestDNAToProteinTransl6(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`QQ`, ""))
	testDNAToProteinHelper(t, dna, protein, 6)
}

func TestDNAToProteinTransl9(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSW`, ""))
	testDNAToProteinHelper(t, dna, protein, 9)
}

func TestDNAToProteinTransl10(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`C`, ""))
	testDNAToProteinHelper(t, dna, protein, 10)
}

func TestDNAToProteinTransl11(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	testDNAToProteinHelper(t, dna, protein, 11)
}

func TestDNAToProteinTransl12(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CTG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`S`, ""))
	testDNAToProteinHelper(t, dna, protein, 12)
}

func TestDNAToProteinTransl13(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGMW`, ""))
	testDNAToProteinHelper(t, dna, protein, 13)
}

func TestDNAToProteinTransl14(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGTAATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSYW`, ""))
	testDNAToProteinHelper(t, dna, protein, 14)
}

func TestDNAToProteinTransl16(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`L`, ""))
	testDNAToProteinHelper(t, dna, protein, 16)
}

func TestDNAToProteinTransl21(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGAATAAGAAGGAAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`WMSSN`, ""))
	testDNAToProteinHelper(t, dna, protein, 21)
}

func TestDNAToProteinTransl22(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TCATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`*L`, ""))
	testDNAToProteinHelper(t, dna, protein, 22)
}

func TestDNAToProteinTransl23(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	testDNAToProteinHelper(t, dna, protein, 23)
}

func TestDNAToProteinTransl24(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SKW`, ""))
	testDNAToProteinHelper(t, dna, protein, 24)
}

func TestDNAToProteinTransl25(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`G`, ""))
	testDNAToProteinHelper(t, dna, protein, 25)
}

func TestDNAToProteinTransl26(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CTG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`A`, ""))
	testDNAToProteinHelper(t, dna, protein, 26)
}

func TestDNAToProteinTransl27(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAGTAATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG*`, ""))
	testDNAToProteinHelper(t, dna, protein, 27)
}

func TestDNAToProteinTransl28(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGW`, ""))
	testDNAToProteinHelper(t, dna, protein, 28)
}

func TestDNAToProteinTransl29(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YY`, ""))
	testDNAToProteinHelper(t, dna, protein, 29)
}

func TestDNAToProteinTransl30(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG`, ""))
	testDNAToProteinHelper(t, dna, protein, 30)
}

func TestDNAToProteinTransl31(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGATAGTAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YGG`, ""))
	testDNAToProteinHelper(t, dna, protein, 31)
}
