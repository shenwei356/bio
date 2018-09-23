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
	transl_table := 99
	_, err := DNAToProtein(dna, transl_table)
	if err == nil {
		t.Error("invalid transl table number should have been detected")
		return
	}
}

func TestDNAToProteinTransl1(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 1)
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

func TestDNAToProteinTransl2(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`**MW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 2)
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

func TestDNAToProteinTransl3(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`ATACTTCTCCTACTGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`MTTTTW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 3)
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

func TestDNAToProteinTransl4(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACCTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FASTW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 4)
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

func TestDNAToProteinTransl5(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SSMW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 5)
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

func TestDNAToProteinTransl6(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`QQ`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 6)
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

func TestDNAToProteinTransl9(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 9)
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

func TestDNAToProteinTransl10(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`C`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 10)
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

func TestDNAToProteinTransl11(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 11)
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

func TestDNAToProteinTransl12(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CTG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`S`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 12)
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

func TestDNAToProteinTransl13(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGATATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGMW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 13)
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

func TestDNAToProteinTransl14(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGTAATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSYW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 14)
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

func TestDNAToProteinTransl16(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`L`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 16)
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

func TestDNAToProteinTransl21(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGAATAAGAAGGAAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`WMSSN`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 21)
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

func TestDNAToProteinTransl22(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TCATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`*L`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 22)
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

func TestDNAToProteinTransl23(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TTTGCTTCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 23)
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

func TestDNAToProteinTransl24(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SKW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 24)
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

func TestDNAToProteinTransl25(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`G`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 25)
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

func TestDNAToProteinTransl26(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CTG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`A`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 26)
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

func TestDNAToProteinTransl27(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAGTAATGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG*`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 27)
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

func TestDNAToProteinTransl28(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAGTGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGW`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 28)
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

func TestDNAToProteinTransl29(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YY`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 29)
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

func TestDNAToProteinTransl30(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TAATAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 30)
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

func TestDNAToProteinTransl31(t *testing.T) {
	dna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`TGATAGTAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YGG`, ""))
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
	translatedProtein, err := DNAToProtein(dna, 31)
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
