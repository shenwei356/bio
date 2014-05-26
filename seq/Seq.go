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

	length int // avoid repeat computation. Only set in constructor
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
	seq.length = len(s)
	return seq, nil
}

// Return the lenght of sequence
func (seq *Seq) Length() int {
	return seq.length
}

// Return reverse complement sequence
func (seq *Seq) Revcom() *Seq {
	return seq.Complement().Reverse()
}

// Reverse sequence
func (seq *Seq) Reverse() *Seq {
	s := ReverseByteSlice(seq.Seq)

	newseq, _ := NewSeq(seq.Alphabet, s)
	return newseq
}

// Return complement sequence
func (seq *Seq) Complement() *Seq {
	s := make([]byte, seq.length)
	var p byte
	for i := 0; i < seq.length; i++ {
		p, _ = seq.Alphabet.pairLetters[seq.Seq[i]]
		s[i] = p
	}

	newseq, _ := NewSeq(seq.Alphabet, s)
	return newseq
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

	return float64(sum) / float64(seq.length)
}
