/*

This package defines a *Seq* type, and provides some basic operations of sequence,
like validation of DNA/RNA/Protein sequence and getting reverse complement sequence.

This package was inspired by
[biogo](https://code.google.com/p/biogo/source/browse/#git%2Falphabet).

[IUPAC nucleotide code](http://droog.gs.washington.edu/parc/images/iupac.html): ACGTURYSWKMBDHVN


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

[IUPAC amino acid code]j(http://www.bioinformatics.org/sms/iupac.html): ACGTRYSWKMBDHV


	A	Ala	Alanine
	C	Cys	Cysteine
	D	Asp	Aspartic Acid
	E	Glu	Glutamic Acid
	F	Phe	Phenylalanine
	G	Gly	Glycine
	H	His	Histidine
	I	Ile	Isoleucine
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

Other links:

- (http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html)

*/
package seq

import (
	"errors"
	"fmt"
)

/*
Four types of alphabets are pre-defined:

	DNA           Deoxyribonucleotide code
	DNAredundant  DNA + Ambiguity Codes
	RNA           Oxyribonucleotide code
	RNAredundant  RNA + Ambiguity Codes
	Protein       Amino Acide single-letter Code
	Unlimit       Self-defined, to including all 26 English letters

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
		[]byte("acgtACGT."),
		[]byte("tgcaTGCA."),
		[]byte(" -"),
		[]byte("nxNX"))

	DNAredundant, _ = NewAlphabet(
		"DNAredundant",
		[]byte("acgtryswkmbdhvACGTRYSWKMBDHV."),
		[]byte("tgcayrswmkvhdbTGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	RNA, _ = NewAlphabet(
		"RNA",
		[]byte("acguACGU"),
		[]byte("ugcaUGCA"),
		[]byte(" -"),
		[]byte("nxNX"))

	RNAredundant, _ = NewAlphabet(
		"RNAredundant",
		[]byte("acguryswkmbdhvACGuRYSWKMBDHV."),
		[]byte("ugcayrswmkvhdbUGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	Protein, _ = NewAlphabet(
		"Protein",
		[]byte("acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY*_"),
		[]byte("acdefghiklmnpqrstvwyACDEFGHIKLMNPQRSTVWY*_"),
		[]byte(" -"),
		[]byte("xX"))

	Unlimit, _ = NewAlphabet(
		"Unlimit",
		[]byte("abcdefghijklmnopqrstuvwxyz*_."),
		[]byte("abcdefghijklmnopqrstuvwxyz*_."),
		[]byte(" -"),
		[]byte(""))
}

/*
type *Alphabet*, you can defined yours. Attention that,
**the letter is case sensitive**.

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
	letters   []byte
	pairs     []byte
	gap       []byte
	ambiguous []byte

	pairLetters map[byte]byte
}

/*
Constructor for type *Alphabet*
*/
func NewAlphabet(
	t string,
	letters []byte,
	pairs []byte,
	gap []byte,
	ambiguous []byte,
) (*Alphabet, error) {

	a := &Alphabet{t, letters, pairs, gap, ambiguous, nil}

	if len(letters) != len(pairs) {
		return a, errors.New("mismarch of length of letters and pairs")
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

// Return type of the alphabet
func (a *Alphabet) Type() string {
	return a.t
}

// Return type of the alphabet
func (a *Alphabet) String() string {
	return a.t
}

// Validate a letter
func (a *Alphabet) IsvalidLetter(b byte) bool {
	_, ok := a.pairLetters[b]
	return ok
}

// Validate a byte slice
func (a *Alphabet) IsValid(s []byte) bool {
	for _, b := range s {
		if !a.IsvalidLetter(b) {
			return false
		}
	}
	return true
}

// Return the Pair Letter
func (a *Alphabet) PairLetter(b byte) (byte, error) {
	if !a.IsvalidLetter(b) {
		return b, errors.New(fmt.Sprintf("invalid letter: %c", b))
	}
	v, _ := a.pairLetters[b]
	return v, nil
}
