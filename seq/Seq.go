package seq

import (
	"bytes"
	"errors"
	"strings"

	"github.com/shenwei356/util/byteutil"
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

// SubSeq returns a sub seq. start and end is 1-based
func (seq *Seq) SubSeq(start int, end int) *Seq {
	if start < 1 {
		start = 1
	}
	if end > len(seq.Seq) {
		end = len(seq.Seq)
	}
	newseq, _ := NewSeq(seq.Alphabet, seq.Seq[start-1:end])
	return newseq
}

// RemoveGaps remove gaps
func (seq *Seq) RemoveGaps(letters string) *Seq {
	m := make(map[byte]bool)
	for i := 0; i < len(letters); i++ {
		m[letters[i]] = true
	}

	s := []byte{}
	for _, b := range seq.Seq {
		if _, ok := m[b]; !ok {
			s = append(s, b)
		}
	}
	newseq, _ := NewSeq(seq.Alphabet, s)
	return newseq
}

// RevCom returns reverse complement sequence
func (seq *Seq) RevCom() *Seq {
	return seq.Complement().Reverse()
}

// Reverse a sequence
func (seq *Seq) Reverse() *Seq {
	s := byteutil.ReverseByteSlice(seq.Seq)

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

// FormatSeq wrap seq
func (seq *Seq) FormatSeq(width int) []byte {
	return byteutil.WrapByteSlice(seq.Seq, width)
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

// DegenerateBaseMapNucl mappings nucleic acid degenerate base to
// regular expression
var DegenerateBaseMapNucl = map[byte]string{
	'A': "A",
	'T': "T",
	'U': "U",
	'C': "C",
	'G': "G",
	'R': "[AG]",
	'Y': "[CT]",
	'M': "[AC]",
	'K': "[GT]",
	'S': "[CG]",
	'W': "[AT]",
	'H': "[ACT]",
	'B': "[CGT]",
	'V': "[ACG]",
	'D': "[AGT]",
	'N': "[ACGT]",
	'a': "a",
	't': "t",
	'u': "u",
	'c': "c",
	'g': "g",
	'r': "[ag]",
	'y': "[ct]",
	'm': "[ac]",
	'k': "[gt]",
	's': "[cg]",
	'w': "[at]",
	'h': "[act]",
	'b': "[cgt]",
	'v': "[acg]",
	'd': "[agt]",
	'n': "[acgt]",
}

// DegenerateBaseMapProt mappings protein degenerate base to
// regular expression
var DegenerateBaseMapProt = map[byte]string{
	'A': "A",
	'B': "[DN]",
	'C': "C",
	'D': "D",
	'E': "E",
	'F': "F",
	'G': "G",
	'H': "H",
	'I': "I",
	'J': "[IL]",
	'K': "K",
	'L': "L",
	'M': "M",
	'N': "N",
	'P': "P",
	'Q': "Q",
	'R': "R",
	'S': "S",
	'T': "T",
	'V': "V",
	'W': "W",
	'Y': "Y",
	'Z': "[QE]",
	'a': "a",
	'b': "[dn]",
	'c': "c",
	'd': "d",
	'e': "e",
	'f': "f",
	'g': "g",
	'h': "h",
	'i': "i",
	'j': "[il]",
	'k': "k",
	'l': "l",
	'm': "m",
	'n': "n",
	'p': "p",
	'q': "q",
	'r': "r",
	's': "s",
	't': "t",
	'v': "v",
	'w': "w",
	'y': "y",
	'z': "[qe]",
}

// Degenerate2Regexp transform seqs containing degenrate base to regular expression
func (seq *Seq) Degenerate2Regexp() string {
	var m map[byte]string
	if seq.Alphabet == Protein {
		m = DegenerateBaseMapProt
	} else {
		m = DegenerateBaseMapNucl
	}

	s := make([]string, len(seq.Seq))
	for i, base := range seq.Seq {
		if _, ok := m[base]; ok {
			s[i] = m[base]
		} else {
			s[i] = string(base)
		}
	}
	return strings.Join(s, "")
}
