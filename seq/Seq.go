package seq

import (
	"bytes"
	"errors"
)

/*
*Seq* has three attributes, alphabet, seq,
and length (avoid compute repeatly)
*/
type Seq struct {
	Alphabet *Alphabet
	Seq      []byte

	Len int // avoid repeat computation
}

/*
Constructor for type *Seq*
*/
func NewSeq(t *Alphabet, s []byte) (*Seq, error) {
	// check sequene first
	if !t.IsValid(s) {
		return nil, errors.New("invalid " + t.Type() + " sequence")
	}

	seq := &Seq{t, s, 0}
	seq.Len = len(s)
	return seq, nil
}

// Return reverse complement sequence
func (seq *Seq) Revcom() []byte {
	return ReverseByteSlice(seq.Complement())
}

// Reverse sequence
func (seq *Seq) Reverse() []byte {
	return ReverseByteSlice(seq.Seq)
}

// Return complement sequence
func (seq *Seq) Complement() []byte {
	s := make([]byte, seq.Len)
	var p byte
	for i := 0; i < seq.Len; i++ {
		p, _ = seq.Alphabet.pairLetters[seq.Seq[i]]
		s[i] = p
	}
	return s
}

/* Compute base content. For example:

   seq.BaseContent("gc")

*/
func (seq *Seq) BaseContent(list string) float64 {
	sum := 0
	for _, b := range []byte(list) {
		up := bytes.ToUpper([]byte{b})
		lo := bytes.ToLower([]byte{b})
		sum += bytes.Count(seq.Seq, up) + bytes.Count(seq.Seq, lo)
	}

	return float64(sum) / float64(seq.Len)
}
