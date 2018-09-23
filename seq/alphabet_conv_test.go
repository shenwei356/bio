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

func TestRNAToProteinInvalidInputs(t *testing.T) {
	rna := []byte("ATGGGCAAAGGCAAAGGC")
	transl_table := 99
	_, err := RNAToProtein(rna, transl_table)
	if err == nil {
		t.Error("invalid transl table number should have been detected")
		return
	}
}

func TestRNAToProteinTransl1(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UUUGCUUCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 1)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}

}

func TestRNAToProteinTransl2(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGAUAUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`**MW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 2)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl3(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AUACUUCUCCUACUGUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`MTTTTW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 3)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl4(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UUUGCUUCAACCUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FASTW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 4)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl5(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGAUAUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SSMW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 5)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl6(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAAUAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`QQ`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 6)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl9(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 9)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl10(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`C`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 10)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl11(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UUUGCUUCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 11)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl12(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CUG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`S`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 12)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl13(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGAUAUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGMW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 13)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl14(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AAAAGAAGGUAAUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`NSSYW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 14)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl16(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`L`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 16)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl21(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UGAAUAAGAAGGAAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`WMSSN`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 21)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl22(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UCAUAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`*L`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 22)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl23(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UUUGCUUCAACC`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`FAST`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 23)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl24(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`AGAAGGUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`SKW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 24)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl25(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`G`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 25)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl26(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`CUG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`A`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 26)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl27(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAGUAAUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG*`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 27)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl28(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAAUAGUGA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GGW`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 28)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl29(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAAUAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YY`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 29)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl30(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UAAUAG`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`GG`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 30)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}

func TestRNAToProteinTransl31(t *testing.T) {
	rna := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`UGAUAGUAA`, ""))
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`YGG`, ""))
	rnaValid := RNA.IsValid(rna) == nil
	if !rnaValid {
		t.Error("validating RNA sequence failed.")
		return
	}
	proteinValid := Protein.IsValid(protein) == nil
	if !proteinValid {
		t.Error("validating protein sequence failed.")
		return
	}
	translatedProtein, err := RNAToProtein(rna, 31)
	if err != nil {
		t.Error("error converting RNA to protein")
		return
	}
	tranlatedProteinValid := Protein.IsValid(translatedProtein) == nil
	if !tranlatedProteinValid {
		t.Error("validating translated protein sequence failed.")
		return
	}
	if !bytes.Equal(protein, translatedProtein) {
		t.Error("converting RNA to protein doesnt match")
		return
	}
}
