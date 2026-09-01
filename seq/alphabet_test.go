package seq

import (
	"bytes"
	"fmt"
	"regexp"
	"strings"
	"testing"
)

func TestAlphabet(t *testing.T) {
	dna := []byte("acgtACGT")
	dna2 := []byte("ACGTRYSWKMBDHV")
	rna := []byte("AUGCUUCCGGCCUGUUCCCUGAGACCUCAAGUGUGAGUGUACUAU" +
		"UGAUGCUUCACACCUGGGCUCUCCGGGUACCAGGACGGUUUGAGCAGAU")
	rna2 := []byte("ACGURYSWKMBDHV")
	protein := []byte(regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`MGLNRFMRAMMVVFITANCITINPDIIFAATDSEDSSLNTDEWEEEKTEEQPSEVNTGPR
YETAREVSSRDIKELEKSNKVRNTNKADLIAMLKEKAEKGPNINNNNSEQTENAAINEEA
SGADRPAIQVERRHPGLPSDSAAEIKKRRKAIASSDSELESLTYPDKPTKVNKKKVAKES
VADASESDLDSSMQSADESSPQPLKANQQPFFPKVFKKIKDAGKWVRDKIDENPEVKKAI
VDKSAGLIDQLLTKKKSEEVNASDFPPPPTDEELRLALPETPMLLGFNAPATSEPSSFEF
PPPPTDEELRLALPETPMLLGFNAPATSEPSSFEFPPPPTEDELEIIRETASSLDSSFTR
GDLASLRNAINRHSQNFSDFPPIPTEEELNGRGGRPTSEEFSSLNSGDFTDDENSETTEE
EIDRLADLRDRGTGKHSRNAGFLPLNPFASSPVPSLSPKVSKISAPALISDITKKTPFKN
PSQPLNVFNKKTTTKTVTKKPTPVKTAPKLAELPATKPQETVLRENKTPFIEKQAETNKQ
SINMPSLPVIQKEATESDKEEMKPQTEEKMVEESESANNANGKNRSAGIEEGKLIAKSAE
DEKAKEEPGNHTTLILAMLAIGVFSLGAFIKIIQLRKNN`, ""))

	ok := DNA.IsValid(dna) == nil &&
		DNAredundant.IsValid(dna2) == nil &&
		RNA.IsValid(rna) == nil &&
		RNAredundant.IsValid(rna2) == nil &&
		Protein.IsValid(protein) == nil
	if !ok {
		t.Error("validating sequence failed.")
		return
	}
	// fmt.Println("protein", GuessAlphabet(protein))
	// fmt.Println("dna2", GuessAlphabet(dna2))
	// fmt.Println("dna", GuessAlphabet(dna))
	// fmt.Println("rna2", GuessAlphabet(rna2))
	// fmt.Println("rna", GuessAlphabet(rna))
	ok = GuessAlphabet(dna) == DNA && GuessAlphabet(dna2) == DNAredundant && GuessAlphabet(rna) == RNA && GuessAlphabet(rna2) == RNAredundant && GuessAlphabet(protein) == Protein
	if !ok {
		t.Error("guessing alphabet error")
		return
	}
}

func TestAlphabet2(t *testing.T) {
	s := regexp.MustCompile(`\r?\n|\s`).ReplaceAllString(
		`MGLNRFMRAMMVVFITANCITINPDIIFAATDSEDSSLNTDEWEEEKTEEQPSEVNTGPR
YETAREVSSRDIKELEKSNKVRNTNKADLIAMLKEKAEKGPNINNNNSEQTENAAINEEA
SGADRPAIQVERRHPGLPSDSAAEIKKRRKAIASSDSELESLTYPDKPTKVNKKKVAKES
VADASESDLDSSMQSADESSPQPLKANQQPFFPKVFKKIKDAGKWVRDKIDENPEVKKAI
VDKSAGLIDQLLTKKKSEEVNASDFPPPPTDEELRLALPETPMLLGFNAPATSEPSSFEF
PPPPTDEELRLALPETPMLLGFNAPATSEPSSFEFPPPPTEDELEIIRETASSLDSSFTR
GDLASLRNAINRHSQNFSDFPPIPTEEELNGRGGRPTSEEFSSLNSGDFTDDENSETTEE
EIDRLADLRDRGTGKHSRNAGFLPLNPFASSPVPSLSPKVSKISAPALISDITKKTPFKN
PSQPLNVFNKKTTTKTVTKKPTPVKTAPKLAELPATKPQETVLRENKTPFIEKQAETNKQ
SINMPSLPVIQKEATESDKEEMKPQTEEKMVEESESANNANGKNRSAGIEEGKLIAKSAE
DEKAKEEPGNHTTLILAMLAIGVFSLGAFIKIIQLRKNN>`, "")

	protein := []byte(strings.Repeat(s, 1000))

	ok := Protein.IsValid(protein) == nil
	if ok {
		t.Error("validating sequence failed.")
		return
	}

}

func TestAlphabetParallelValidationReportsFirstInvalidLetter(t *testing.T) {
	oldThreshold, oldThreads := ValidSeqLengthThreshold, ValidSeqThreads
	defer func() {
		ValidSeqLengthThreshold = oldThreshold
		ValidSeqThreads = oldThreads
	}()
	ValidSeqLengthThreshold = 1
	ValidSeqThreads = 4

	sequence := bytes.Repeat([]byte{'A'}, (64<<10)*4)
	sequence[123] = '!'
	sequence[len(sequence)-10] = '?'

	err := DNA.IsValid(sequence)
	if err == nil {
		t.Fatal("invalid sequence passed validation")
	}
	if !strings.Contains(err.Error(), "at 123") {
		t.Fatalf("reported the wrong invalid position: %v", err)
	}
}

func TestAlphabetShortParallelValidationAvoidsWorkers(t *testing.T) {
	oldThreshold, oldThreads := ValidSeqLengthThreshold, ValidSeqThreads
	defer func() {
		ValidSeqLengthThreshold = oldThreshold
		ValidSeqThreads = oldThreads
	}()
	ValidSeqLengthThreshold = 64 << 10
	ValidSeqThreads = 16

	sequence := bytes.Repeat([]byte{'A'}, 10<<10)
	if allocs := testing.AllocsPerRun(1000, func() {
		if err := DNA.IsValid(sequence); err != nil {
			panic(err)
		}
	}); allocs != 0 {
		t.Fatalf("short validation allocated %.1f times per call", allocs)
	}
}

func BenchmarkAlphabetIsValid(b *testing.B) {
	oldThreshold, oldThreads := ValidSeqLengthThreshold, ValidSeqThreads
	defer func() {
		ValidSeqLengthThreshold = oldThreshold
		ValidSeqThreads = oldThreads
	}()
	ValidSeqLengthThreshold = 64 << 10
	ValidSeqThreads = 16

	for _, size := range []int{10 << 10, 32 << 10, 64 << 10, 100 << 10, 1 << 20} {
		sequence := bytes.Repeat([]byte{'A'}, size)
		b.Run(fmt.Sprintf("%d", size), func(b *testing.B) {
			b.SetBytes(int64(size))
			b.ReportAllocs()
			for i := 0; i < b.N; i++ {
				if err := DNA.IsValid(sequence); err != nil {
					b.Fatal(err)
				}
			}
		})
	}
}
