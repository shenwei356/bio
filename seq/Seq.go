package seq

import (
	"bytes"
	"errors"
)

// Seq struct has two attributes, alphabet, seq,
type Seq struct {
	Alphabet *Alphabet
	Seq      []byte
}

// NewSeq is constructor for type *Seq*
func NewSeq(t *Alphabet, s []byte) (*Seq, error) {
	// check sequene first
	if !t.IsValid(s) {
		return nil, errors.New("invalid " + t.Type() + " sequence")
	}

	seq := &Seq{t, s}
	return seq, nil
}

// Length returns the lenght of sequence
func (seq *Seq) Length() int {
	return len(seq.Seq)
}

// Revcom returns reverse complement sequence
func (seq *Seq) Revcom() *Seq {
	return seq.Complement().Reverse()
}

// Reverse a sequence
func (seq *Seq) Reverse() *Seq {
	s := ReverseByteSlice(seq.Seq)

	newseq, _ := NewSeq(seq.Alphabet, s)
	return newseq
}

// Complement returns complement sequence
func (seq *Seq) Complement() *Seq {
	if seq.Alphabet == Unlimit {
		newseq, _ := NewSeq(seq.Alphabet, []byte(""))
		return newseq
	}

	s := make([]byte, len(seq.Seq))
	var p byte
	for i := 0; i < len(seq.Seq); i++ {
		p, _ = seq.Alphabet.pairLetters[seq.Seq[i]]
		s[i] = p
	}

	newseq, _ := NewSeq(seq.Alphabet, s)
	return newseq
}

/*BaseContent returns base content for given bases. For example:

  seq.BaseContent("gc")

*/
func (seq *Seq) BaseContent(list string) float64 {
	if len(seq.Seq) == 0 {
		return float64(0)
	}

	sum := 0
	for _, b := range []byte(list) {
		up := bytes.ToUpper([]byte{b})
		lo := bytes.ToLower([]byte{b})
		if string(up) == string(lo) {
			sum += bytes.Count(seq.Seq, up)
		} else {
			sum += bytes.Count(seq.Seq, up) + bytes.Count(seq.Seq, lo)
		}
	}

	return float64(sum) / float64(len(seq.Seq))
}
