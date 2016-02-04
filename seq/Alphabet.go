/*Package seq balabala

This package defines a *Seq* type, and provides some basic operations of sequence,
like validation of DNA/RNA/Protein sequence and getting reverse complement sequence.

This package was inspired by
[biogo](https://code.google.com/p/biogo/source/browse/#git%2Falphabet).

IUPAC nucleotide code: ACGTURYSWKMBDHVN

http://droog.gs.washington.edu/parc/images/iupac.html


	code	base	Complement
	A	A	T
	C	C	G
	G	G	C
	T/U	T	A

	R	A/G	Y
	Y	C/T	R
	S	C/G	S
	W	A/T	W
	K	G/T	M
	M	A/C	K

	B	C/G/T	V
	D	A/G/T	H
	H	A/C/T	D
	V	A/C/G	B

	X/N	A/C/G/T	X
	.	not A/C/G/T
	 or-	gap

IUPAC amino acid code: `ACGTRYSWKMBDHV`

	A	Ala	Alanine
	B	Asx	Aspartic acid or Asparagine [2]
	C	Cys	Cysteine
	D	Asp	Aspartic Acid
	E	Glu	Glutamic Acid
	F	Phe	Phenylalanine
	G	Gly	Glycine
	H	His	Histidine
	I	Ile	Isoleucine
	J		Isoleucine or Leucine [4]
	K	Lys	Lysine
	L	Leu	Leucine
	M	Met	Methionine
	N	Asn	Asparagine
	P	Pro	Proline
	Q	Gln	Glutamine
	R	Arg	Arginine
	S	Ser	Serine
	T	Thr	Threonine
	V	Val	Valine
	W	Trp	Tryptophan
	Y	Tyr	Tyrosine
	Z	Glx	Glutamine or Glutamic acid [2]

Other links:

	1. http://www.bioinformatics.org/sms/iupac.html
	2. http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	3. http://www.bioinformatics.org/sms2/iupac.html
	4. http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html

*/
package seq

import (
	"errors"
	"fmt"
)

/*Alphabet could be defined. Attention that,
**the letters are case sensitive**.

For exmaple, DNA:

	DNA, _ = NewAlphabet(
		"DNA",
		[]byte("acgtACGT"),
		[]byte("tgcaTGCA"),
		[]byte(" -"),
		[]byte("nN"))

*/
type Alphabet struct {
	t         string
	isUnlimit bool

	letters   []byte
	pairs     []byte
	gap       []byte
	ambiguous []byte

	pairLetters map[byte]byte
}

// NewAlphabet is Constructor for type *Alphabet*
func NewAlphabet(
	t string,
	isUnlimit bool,
	letters []byte,
	pairs []byte,
	gap []byte,
	ambiguous []byte,
) (*Alphabet, error) {

	a := &Alphabet{t, isUnlimit, letters, pairs, gap, ambiguous, nil}

	if isUnlimit {
		return a, nil
	}

	if len(letters) != len(pairs) {
		return a, errors.New("mismatch of length of letters and pairs")
	}

	a.pairLetters = make(map[byte]byte, len(letters))
	for i := 0; i < len(letters); i++ {
		a.pairLetters[letters[i]] = pairs[i]
	}

	// add gap and ambiguous code
	for _, v := range gap {
		a.pairLetters[v] = v
	}
	for _, v := range ambiguous {
		a.pairLetters[v] = v
	}

	return a, nil
}

// Type return type of the alphabet
func (a *Alphabet) Type() string {
	return a.t
}

// Return type of the alphabet
func (a *Alphabet) String() string {
	return a.t
}

// IsValidLetter is used to validate a letter
func (a *Alphabet) IsValidLetter(b byte) bool {
	if a.isUnlimit {
		return true
	}
	_, ok := a.pairLetters[b]
	return ok
}

// IsValid is used to validate a byte slice
func (a *Alphabet) IsValid(s []byte) bool {
	if a.isUnlimit {
		return true
	}

	for _, b := range s {
		if !a.IsValidLetter(b) {
			return false
		}
	}
	return true
}

// PairLetter returns the Pair Letter
func (a *Alphabet) PairLetter(b byte) (byte, error) {
	if a.isUnlimit {
		return b, nil
	}

	if !a.IsValidLetter(b) {
		return b, fmt.Errorf("invalid letter: %c", b)
	}
	v, _ := a.pairLetters[b]
	return v, nil
}

/*Four types of alphabets are pre-defined:

  DNA           Deoxyribonucleotide code
  DNAredundant  DNA + Ambiguity Codes
  RNA           Oxyribonucleotide code
  RNAredundant  RNA + Ambiguity Codes
  Protein       Amino Acide single-letter Code
  Unlimit       Self-defined, including all 26 English letters

*/
var (
	DNA          *Alphabet
	DNAredundant *Alphabet
	RNA          *Alphabet
	RNAredundant *Alphabet
	Protein      *Alphabet
	Unlimit      *Alphabet
)

func init() {
	DNA, _ = NewAlphabet(
		"DNA",
		false,
		[]byte("acgtACGT."),
		[]byte("tgcaTGCA."),
		[]byte(" -"),
		[]byte("nxNX"))

	DNAredundant, _ = NewAlphabet(
		"DNAredundant",
		false,
		[]byte("acgtryswkmbdhvACGTRYSWKMBDHV."),
		[]byte("tgcayrswmkvhdbTGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	RNA, _ = NewAlphabet(
		"RNA",
		false,
		[]byte("acguACGU"),
		[]byte("ugcaUGCA"),
		[]byte(" -"),
		[]byte("nxNX"))

	RNAredundant, _ = NewAlphabet(
		"RNAredundant",
		false,
		[]byte("acguryswkmbdhvACGURYSWKMBDHV."),
		[]byte("ugcayrswmkvhdbUGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	Protein, _ = NewAlphabet(
		"Protein",
		false,
		[]byte("abcdefghijklmnpqrstvwyzABCDEFGHIJKLMNPQRSTVWYZ*_."),
		[]byte("abcdefghijklmnpqrstvwyzABCDEFGHIJKLMNPQRSTVWYZ*_."),
		[]byte(" -"),
		[]byte("xX"))

	Unlimit, _ = NewAlphabet(
		"Unlimit",
		true,
		nil,
		nil,
		nil,
		nil)
}
