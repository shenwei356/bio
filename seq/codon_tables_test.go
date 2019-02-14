package seq

import (
	"regexp"
	"testing"
)

type codonTableTest struct {
	table int
	nt    string
	aa    string
	frame int
	trim  bool
	clean bool

	allowUnknownCodon bool
	markInitCodonAsM  bool
}

var codonTableTests []codonTableTest

func init() {
	re := regexp.MustCompile(`\r?\n|\s`)
	codonTableTests = make([]codonTableTest, 0, 10)

	// https://www.ncbi.nlm.nih.gov/nuccore/AB021961.1
	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 1,
			nt: re.ReplaceAllString(`atgactgccatggaggagtcacagtcggatatcagcctcgagctccctctgagccaggag
acattttcaggcttatggaaactacttcctccagaagatatcctgccatcacctcactgc
atggacgatctgttgctgccccaggatgttgaggagttttttgaaggcccaagtgaagcc
ctccgagtgtcaggagctcctgcagcacaggaccctgtcaccgagacccctgggccagtg
gcccctgccccagccactccatggcccctgtcatcttttgtcccttctcaaaaaacttac
cagggcaactatggcttccacctgggcttcctgcagtctgggacagccaagtctgttatg
tgcacgtactctcctcccctcaataagctattctgccagctggcgaagacgtgccctgtg
cagttgtgggtcagcgccacacctccagctgggagccgtgtccgcgccatggccatctac
aagaagtcacagcacatgacggaggtcgtgagacgctgcccccaccatgagcgctgctcc
gatggtgatggcctggctcctccccagcatcgtatccgggtggaaggaaatttgtatccc
gagtatctggaagacaggcagacttttcgccacagcgtggtggtaccttatgagccaccc
gaggccggctctgagtataccaccatccactacaagtacatgtgtaatagctcctgcatg
gggggcatgaaccgccgacctatccttaccatcatcacactggaagactccagtgggaac
cttctgggacgggacagctttgaggttcgtgtttgtgcctgccctgggagagaccgccgt
acagaagaagaaaatttccgcaaaaaggaagtcctttgccctgaactgcccccagggagc
gcaaagagagcgctgcccacctgcacaagcgcctctcccccgcaaaagaaaaaaccactt
gatggagagtatttcaccctcaagatccgcgggcgtaaacgcttcgagatgttccgggag
ctgaatgaggccttagagttaaaggatgcccatgctacagaggagtctggagacagcagg
gctcactccagctacctgaagaccaagaagggccagtctacttcccgccataaaaaaaca
atggtcaagaaagtggggcctgactcagactga`, ""),
			aa: re.ReplaceAllString(`MTAMEESQSDISLELPLSQETFSGLWKLLPPEDILPSPHCMDDL
LLPQDVEEFFEGPSEALRVSGAPAAQDPVTETPGPVAPAPATPWPLSSFVPSQKTYQG
NYGFHLGFLQSGTAKSVMCTYSPPLNKLFCQLAKTCPVQLWVSATPPAGSRVRAMAIY
KKSQHMTEVVRRCPHHERCSDGDGLAPPQHRIRVEGNLYPEYLEDRQTFRHSVVVPYE
PPEAGSEYTTIHYKYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRDSFEVRVCACPG
RDRRTEEENFRKKEVLCPELPPGSAKRALPTCTSASPPQKKKPLDGEYFTLKIRGRKR
FEMFRELNEALELKDAHATEESGDSRAHSSYLKTKKGQSTSRHKKTMVKKVGPDSD`, ""),
			frame: 1,
			trim:  true,
			clean: false,
		})

	// ORF25 from https://www.ncbi.nlm.nih.gov/nuccore/AJ245616.1
	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 11,
			nt: re.ReplaceAllString(`atggaggaacaagcatggcgagaagtcctcgaacgtttagctcga
attgaaacaaagttagataactatgaaacagttcgagataaagcagaacgagcgctccta
atagctcaatcaaatgcgaaacttatagaaaaaatggaagctaataataagtgggcttgg
ggctttatgcttactcttgccgtaactgttattggttatttattcactaaaattagattc
tga`, ""),
			aa: re.ReplaceAllString(`MEEQAWREVLERLARIETKLDNYETVRDKAERALLIAQSNAKLI
EKMEANNKWAWGFMLTLAVTVIGYLFTKIRF`, ""),
			frame: 1,
			trim:  true,
			clean: false,
		})

	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 11,
			nt: re.ReplaceAllString(`tcagaatctaattttagtgaataaataaccaataacagttacggcaagagtaagcataaa
gccccaagcccacttattattagcttccattttttctataagtttcgcatttgattgagc
tattaggagcgctcgttctgctttatctcgaactgtttcatagttatctaactttgtttc
aattcgagctaaacgttcgaggacttctcgccatgcttgttcctccat`, ""),
			aa: re.ReplaceAllString(`MEEQAWREVLERLARIETKLDNYETVRDKAERALLIAQSNAKLI
EKMEANNKWAWGFMLTLAVTVIGYLFTKIRF`, ""),
			frame: -1,
			trim:  true,
			clean: false,
		})

	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 11,
			nt: re.ReplaceAllString(`tcagaatctaattttagtgaataaataaccaataacagttacggcaagagtaagcataaa
gccccaagcccacttattattagcttccattttttctataagtttcgcatttgattgagc
tattaggagcgctcgttctgctttatctcgaactgtttcatagttatctaactttgtttc
aattcgagctaaacgttcgaggacttctcgccatgcttgttcctccatc`, ""),
			aa: re.ReplaceAllString(`MEEQAWREVLERLARIETKLDNYETVRDKAERALLIAQSNAKLI
EKMEANNKWAWGFMLTLAVTVIGYLFTKIRF`, ""),
			frame: -2,
			trim:  true,
			clean: false,
		})

	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 11,
			nt: re.ReplaceAllString(`tcagaatctaattttagtgaataaataaccaataacagttacggcaagagtaagcataaa
gccccaagcccacttattattagcttccattttttctataagtttcgcatttgattgagc
tattaggagcgctcgttctgctttatctcgaactgtttcatagttatctaactttgtttc
aattcgagctaaacgttcgaggacttctcgccatgcttgttcctccatcc`, ""),
			aa: re.ReplaceAllString(`MEEQAWREVLERLARIETKLDNYETVRDKAERALLIAQSNAKLI
EKMEANNKWAWGFMLTLAVTVIGYLFTKIRF`, ""),
			frame: -3,
			trim:  true,
			clean: false,
		})

	codonTableTests = append(codonTableTests,
		codonTableTest{
			table: 11,
			nt: re.ReplaceAllString(`atNgaggaacaagcatggcgagaagtcctcgaacgtttagctcga
attgaaacaaagttagataactatgaaacagttcgagataaagcagaacgagcgctccta
atagctcaatcaaatgcgaaacttatagaaaaaatggaagctaataataagtgggcttgg
ggctttatgcttactcttgccgtaactgttattggttatttattcactaaaattagattc
tga`, ""),
			aa: re.ReplaceAllString(`XEEQAWREVLERLARIETKLDNYETVRDKAERALLIAQSNAKLI
EKMEANNKWAWGFMLTLAVTVIGYLFTKIRF*`, ""),
			frame:             1,
			trim:              false,
			clean:             false,
			allowUnknownCodon: true,
		})
}

func TestCodonTableStranslation(t *testing.T) {
	var aa []byte
	var err error
	for i, test := range codonTableTests {
		aa, err = CodonTables[test.table].Translate([]byte(test.nt), test.frame, test.trim, test.clean, test.allowUnknownCodon, test.markInitCodonAsM)
		if err != nil {
			t.Errorf("test %d err: %s", i, err)
		}
		if len(aa) != len(test.aa) {
			t.Errorf("test %d err: len unequal, answer: %d, result: %d", i, len(aa), len(test.aa))
		}
		if string(aa) != test.aa {
			t.Errorf("test %d err: result not right", i)
		}
	}
}
