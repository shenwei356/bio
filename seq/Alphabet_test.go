package seq

import (
	"regexp"
	"testing"
)

func Test(t *testing.T) {
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
