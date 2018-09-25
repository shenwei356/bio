// https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
package seq

import (
	"bytes"
	"errors"
	"fmt"
	"sort"
	"strings"
)

// ErrInvalidCodon means the codon contains bases expcet for A C T G U.
var ErrInvalidCodon = errors.New("seq: invalid codon")

// CodonTable represents a codon table
type CodonTable struct {
	ID         int
	Name       string
	InitCodons map[string]struct{} // codon as string, map for fast quering
	StopCodons map[string]struct{}
	table      [4][4][4]byte // matrix is much faster than map for quering
}

// NewCodonTable contructs a CodonTable with ID and Name,
// you need to set the detailed codon table by calling Set or Set2.
func NewCodonTable(id int, name string) *CodonTable {
	t := &CodonTable{ID: id, Name: name}
	t.InitCodons = make(map[string]struct{}, 1)
	t.StopCodons = make(map[string]struct{}, 1)
	t.table = [4][4][4]byte{}
	return t
}

// String returns a short description of the CodonTable.
func (t CodonTable) String() string {
	var b bytes.Buffer
	b.WriteString(fmt.Sprintf("%s (transl_table=%d). Initiation Codon: ", t.Name, t.ID))

	codons := make([]string, len(t.InitCodons))
	i := 0
	for codon := range t.InitCodons {
		codons[i] = codon
		i++
	}
	sort.Strings(codons)
	b.WriteString(strings.Join(codons, ", "))

	b.WriteString(". Stop Codon: ")

	codons = make([]string, len(t.StopCodons))
	i = 0
	for codon := range t.StopCodons {
		codons[i] = codon
		i++
	}
	sort.Strings(codons)
	b.WriteString(strings.Join(codons, ", "))

	b.WriteString(".")

	return b.String()
}

// base2idx converts the base to integer using for indexing the condon table (matrix).
func base2idx(base byte) (int, error) {
	var v int
	switch base {
	case 'A', 'a':
		v = 0
	case 'C', 'c':
		v = 1
	case 'G', 'g':
		v = 2
	case 'T', 't', 'U', 'u':
		v = 3
	default:
		return 0, ErrInvalidCodon
	}
	return v, nil
}

// codon2idx returns the location of a codon in the matrix.
func codon2idx(codon []byte) (int, int, int, error) {
	if len(codon) != 3 {
		return 0, 0, 0, ErrInvalidCodon
	}
	var err error
	var i, j, k int
	i, err = base2idx(codon[0])
	if err != nil {
		return 0, 0, 0, err
	}
	j, err = base2idx(codon[1])
	if err != nil {
		return 0, 0, 0, err
	}
	k, err = base2idx(codon[2])
	if err != nil {
		return 0, 0, 0, err
	}
	return i, j, k, nil
}

// Set sets a codon of byte slice.
func (t *CodonTable) Set(codon []byte, aminoAcid byte) error {
	i, j, k, err := codon2idx(codon)
	if err != nil {
		return err
	}
	t.table[i][j][k] = aminoAcid
	return nil
}

// Set2 sets a codon of string.
func (t *CodonTable) Set2(codon string, aminoAcid byte) error {
	return t.Set([]byte(codon), aminoAcid)
}

// Get returns the amino acid of the codon ([]byte), codon can be DNA or RNA.
func (t *CodonTable) Get(codon []byte) (byte, error) {
	i, j, k, err := codon2idx(codon)
	if err != nil {
		return 0, err
	}
	aa := t.table[i][j][k]
	if aa == 0 {
		aa = 'X'
	}
	return aa, nil
}

// Get2 returns the amino acid of the codon (string), codon can be DNA or RNA.
func (t *CodonTable) Get2(codon string) (byte, error) {
	return t.Get([]byte(codon))
}

// Clone returns a deep copy of the CodonTable.
func (t *CodonTable) Clone() CodonTable {
	initCodons := make(map[string]struct{}, len(t.InitCodons))
	for k, v := range t.InitCodons {
		initCodons[k] = v
	}
	stopCodons := make(map[string]struct{}, len(t.StopCodons))
	for k, v := range t.StopCodons {
		stopCodons[k] = v
	}

	table := [4][4][4]byte{}
	for i := 0; i < 4; i++ {
		for j := 0; j < 4; j++ {
			for k := 0; k < 4; k++ {
				table[i][j][k] = t.table[i][j][k]
			}
		}
	}
	return CodonTable{ID: t.ID, Name: t.Name, InitCodons: initCodons, StopCodons: stopCodons, table: table}
}

// Translate translates a DNA/RNA sequence to amino acid sequences.
// Available frame: 1, 2, 3; for -1, -2 ,-3, please call the Seq.Translate().
// If option trim is true, it removes all 'X' and '*' characters from the right end of the translation.
// If option clean is true, it changes all STOP codon positions from the '*' character to 'X' (an unknown residue).
func (t *CodonTable) Translate(sequence []byte, frame int, trim bool, clean bool) ([]byte, error) {
	if len(sequence) < 3 {
		return nil, fmt.Errorf("seq: sequence too short to translate: %d", len(sequence))
	}
	if frame < 1 || frame > 3 {
		return nil, fmt.Errorf("seq: invalid frame: %d. available: 1, 2, 3", frame)
	}
	aas := make([]byte, 0, int((len(sequence)+2)/3))
	var aa byte
	var err error
	for i := frame - 1; i < len(sequence)-2; i += 3 {
		aa, err = t.Get(sequence[i : i+3])
		if err != nil {
			return nil, err
		}

		if trim && (aa == 'X' || aa == '*') {
			break
		}
		if clean && aa == '*' {
			aa = 'X'
		}

		aas = append(aas, aa)
	}
	return aas, nil
}

// CodonTables contains all the codon tables from NCBI:
//
//     1: The Standard Code
//     2: The Vertebrate Mitochondrial Code
//     3: The Yeast Mitochondrial Code
//     4: The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
//     5: The Invertebrate Mitochondrial Code
//     6: The Ciliate, Dasycladacean and Hexamita Nuclear Code
//     9: The Echinoderm and Flatworm Mitochondrial Code
//     10: The Euplotid Nuclear Code
//     11: The Bacterial, Archaeal and Plant Plastid Code
//     12: The Alternative Yeast Nuclear Code
//     13: The Ascidian Mitochondrial Code
//     14: The Alternative Flatworm Mitochondrial Code
//     16: Chlorophycean Mitochondrial Code
//     21: Trematode Mitochondrial Code
//     22: Scenedesmus obliquus Mitochondrial Code
//     23: Thraustochytrium Mitochondrial Code
//     24: Pterobranchia Mitochondrial Code
//     25: Candidate Division SR1 and Gracilibacteria Code
//     26: Pachysolen tannophilus Nuclear Code
//     27: Karyorelict Nuclear
//     28: Condylostoma Nuclear
//     29: Mesodinium Nuclear
//     30: Peritrich Nuclear
//     31: Blastocrithidia Nuclear
//
var CodonTables map[int]*CodonTable

// https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=cgencodes#SG1
func codonTableFromText(id int, name string, text string) *CodonTable {
	t := NewCodonTable(id, name)
	data := strings.Split(text, "\n")
	if !(len(data) == 5 && len(data[0]) == 64 && len(data[1]) == 64 &&
		len(data[2]) == 64 && len(data[3]) == 64 && len(data[4]) == 64) {
		panic("please paste the right text from NCBI")
	}
	var aa, start byte
	var codon []byte
	for i := 0; i < 64; i++ {
		aa = data[0][i]
		start = data[1][i]
		codon = []byte{data[2][i], data[3][i], data[4][i]}

		t.Set(codon, aa)
		if start == 'M' {
			t.InitCodons[string(codon)] = struct{}{}
		} else if start == '*' {
			t.StopCodons[string(codon)] = struct{}{}
		}
	}
	return t
}

func init() {
	CodonTables = make(map[int]*CodonTable, 31)

	CodonTables[1] = codonTableFromText(1,
		"The Standard Code",
		`FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
---M------**--*----M---------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[2] = codonTableFromText(2,
		"The Vertebrate Mitochondrial Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG
----------**--------------------MMMM----------**---M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[3] = codonTableFromText(3,
		"The Yeast Mitochondrial Code",
		`FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**----------------------MM----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[4] = codonTableFromText(4,
		"The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--MM------**-------M------------MMMM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[5] = codonTableFromText(5,
		"The Invertebrate Mitochondrial Code ",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG
---M------**--------------------MMMM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[6] = codonTableFromText(6,
		"The Ciliate, Dasycladacean and Hexamita Nuclear Code",
		`FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--------------*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[9] = codonTableFromText(9,
		"The Echinoderm and Flatworm Mitochondrial Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
----------**-----------------------M---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[10] = codonTableFromText(10,
		"The Euplotid Nuclear Code",
		`FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**-----------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[11] = codonTableFromText(11,
		"The Bacterial, Archaeal and Plant Plastid Code",
		`FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
---M------**--*----M------------MMMM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[12] = codonTableFromText(12,
		"The Alternative Yeast Nuclear Code",
		`FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**--*----M---------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[13] = codonTableFromText(13,
		"The Ascidian Mitochondrial Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG
---M------**----------------------MM---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[14] = codonTableFromText(14,
		"The Alternative Flatworm Mitochondrial Code",
		`FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG
-----------*-----------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[16] = codonTableFromText(16,
		"Chlorophycean Mitochondrial Code",
		`FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------*---*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[21] = codonTableFromText(21,
		"Trematode Mitochondrial Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG
----------**-----------------------M---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[22] = codonTableFromText(22,
		"Scenedesmus obliquus Mitochondrial Code",
		`FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
------*---*---*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[23] = codonTableFromText(23,
		"Thraustochytrium Mitochondrial Code",
		`FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--*-------**--*-----------------M--M---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[24] = codonTableFromText(24,
		"Pterobranchia Mitochondrial Code",
		`FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG
---M------**-------M---------------M---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[25] = codonTableFromText(25,
		"Candidate Division SR1 and Gracilibacteria Code",
		`FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
---M------**-----------------------M---------------M------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[26] = codonTableFromText(26,
		"Pachysolen tannophilus Nuclear Code",
		`FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**--*----M---------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[27] = codonTableFromText(27,
		"Karyorelict Nuclear",
		`FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--------------*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[28] = codonTableFromText(28,
		"Condylostoma Nuclear",
		`FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**--*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[29] = codonTableFromText(29,
		"Mesodinium Nuclear",
		`FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--------------*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[30] = codonTableFromText(30,
		"Peritrich Nuclear",
		`FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
--------------*--------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	CodonTables[31] = codonTableFromText(31,
		"Blastocrithidia Nuclear",
		`FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG
----------**-----------------------M----------------------------
TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG
TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG
TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG`)

	// ks := make([]int, len(CodonTables))
	// i := 0
	// for k := range CodonTables {
	// 	ks[i] = k
	// 	i++
	// }
	// sort.Ints(ks)
	// for _, i = range ks {
	// 	// fmt.Println(CodonTables[i])
	// 	fmt.Printf("%d: %s\n", CodonTables[i].ID, CodonTables[i].Name)
	// }

}
