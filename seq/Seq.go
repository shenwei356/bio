package seq

import (
	"bytes"
	"errors"
)

type Seq struct {
	Alphabet *Alphabet
	Seq      []byte

	Len int
}

func NewSeq(t *Alphabet, s []byte) (*Seq, error) {
	if !t.IsValid(s) {
		return nil, errors.New("invalid " + t.Type() + " sequence")
	}
	seq := &Seq{t, s, 0}
	seq.Len = len(s)
	return seq, nil
}

func (seq *Seq) Revcom() []byte {
	return ReverseByteSlice(seq.Complement())
}

func (seq *Seq) Reverse() []byte {
	return ReverseByteSlice(seq.Seq)
}

func (seq *Seq) Complement() []byte {
	s := make([]byte, seq.Len)
	var p byte
	for i := 0; i < seq.Len; i++ {
		p, _ = seq.Alphabet.letterPairs[seq.Seq[i]]
		s[i] = p
	}
	return s
}

func (seq *Seq) BaseContent(list []byte) float32 {
	sum := 0
	for _, b := range list {
		up := bytes.ToUpper([]byte{b})
		lo := bytes.ToLower([]byte{b})
		sum += bytes.Count(seq.Seq, up) + bytes.Count(seq.Seq, lo)
	}

	return float32(sum) / float32(len(seq.Seq))
}
