package seq

import (
	"errors"
	"fmt"
)

var (
	DNA          *Alphabet
	DNAredundant *Alphabet
	RNA          *Alphabet
	RNAredundant *Alphabet
	Protein      *Alphabet
)

func init() {
	DNA, _ = NewAlphabet(
		"DNA",
		[]byte("acgtACGT"),
		[]byte("tgcaTGCA"),
		'-',
		[]byte("nN"))

	/*	A	A
		T	T
		C	C
		G	G

		R	AG
		Y	CT
		M	AC
		K	GT
		S	CG
		W	AT

		H	ACT
		B	CGT
		V	ACG
		D	ATG

		N	ACGT

		ATCGRYMKSWHBVD
	*/
	DNAredundant, _ = NewAlphabet(
		"DNAredundant",
		[]byte("acgtrymkswhbvdACGTRYMKSWHBVD"),
		[]byte("tgcayrkmswdvbhTGCAYRKMSWDVBH"),

		// from https://code.google.com/p/biogo/source/browse/alphabet/alphabet.go
		// acmgrsvtwyhkdbnxACMGRSVTWYHKDB
		// tgkcysbawrdmhvnxTGKCYSBAWRDMHV
		'-',
		[]byte("nN"))

	RNA, _ = NewAlphabet(
		"RNA",
		[]byte("acguACGU"),
		[]byte("ugcaUGCA"),
		'-',
		[]byte("nN"))

	RNAredundant, _ = NewAlphabet(
		"RNAredundant",
		[]byte("acgurymkswhbvdACGURYMKSWHBVD"),
		[]byte("ugcayrkmswdvbhUGCAYRKMSWDVBH"),
		'-',
		[]byte("nN"))

	Protein, _ = NewAlphabet(
		"Protein",
		[]byte("abcdefghijklmnpqrstvwyz*ABCDEFGHIJKLMNPQRSTVWYZ"),
		[]byte("abcdefghijklmnpqrstvwyz*ABCDEFGHIJKLMNPQRSTVWYZ"),
		'-',
		[]byte("xX"))

}

type Alphabet struct {
	t         string
	letters   []byte
	pairs     []byte
	gap       byte
	ambiguous []byte

	letterPairs map[byte]byte
}

func NewAlphabet(
	t string,
	letters []byte,
	pairs []byte,
	gap byte,
	ambiguous []byte,
) (*Alphabet, error) {
	a := &Alphabet{t, letters, pairs, gap, ambiguous, nil}

	if len(letters) != len(pairs) {
		return a, errors.New("mismarch of length of letters and pairs")
	}

	a.letterPairs = make(map[byte]byte, len(letters))
	for i := 0; i < len(letters); i++ {
		a.letterPairs[letters[i]] = pairs[i]
	}

	a.letterPairs[gap] = gap
	for _, v := range ambiguous {
		a.letterPairs[v] = v
	}

	return a, nil
}

func (a *Alphabet) Type() string {
	return a.t
}

func (a *Alphabet) IsvalidLetter(b byte) bool {
	_, ok := a.letterPairs[b]
	return ok
}

func (a *Alphabet) IsValid(s []byte) bool {
	for _, b := range s {
		if !a.IsvalidLetter(b) {
			return false
		}
	}
	return true
}

func (a *Alphabet) LetterPair(b byte) (byte, error) {
	if !a.IsvalidLetter(b) {
		return b, errors.New(fmt.Sprintf("invalid letter: %c", b))
	}
	v, _ := a.letterPairs[b]
	return v, nil
}
