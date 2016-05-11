package seq

import (
	"bytes"
	"fmt"
	"runtime"
	"strings"
	"sync"

	"github.com/shenwei356/util/byteutil"
	"github.com/shenwei356/util/stringutil"
)

// Seq struct has two attributes, alphabet, seq,
type Seq struct {
	Alphabet  *Alphabet
	Seq       []byte
	Qual      []byte
	QualValue []int
}

// ValidateSeq decides whether check sequence or not
var ValidateSeq = true

// NewSeq is constructor for type *Seq*
func NewSeq(t *Alphabet, s []byte) (*Seq, error) {
	if ValidateSeq {
		// check sequene first
		if err := t.IsValid(s); err != nil {
			return nil, err
		}
	}

	seq := &Seq{Alphabet: t, Seq: s}
	return seq, nil
}

// NewSeqWithQual is used to store fastq sequence
func NewSeqWithQual(t *Alphabet, s []byte, q []byte) (*Seq, error) {
	if len(s) != len(q) {
		return nil, fmt.Errorf("unmatched length of sequence (%d) and quality (%d)", len(s), len(q))
	}
	seq, err := NewSeq(t, s)
	if err != nil {
		return nil, err
	}
	seq.Qual = q
	return seq, nil
}

// NewSeqWithoutValidate create Seq without check the sequences
func NewSeqWithoutValidate(t *Alphabet, s []byte) (*Seq, error) {
	seq := &Seq{Alphabet: t, Seq: s}
	return seq, nil
}

// NewSeqWithQualWithoutValidate create Seq with quality without check the sequences
func NewSeqWithQualWithoutValidate(t *Alphabet, s []byte, q []byte) (*Seq, error) {
	if len(s) != len(q) {
		return nil, fmt.Errorf("unmatched length of sequence (%d) and quality (%d)", len(s), len(q))
	}
	seq := &Seq{Alphabet: t, Seq: s, Qual: q}
	return seq, nil
}

// Length returns the lenght of sequence
func (seq *Seq) Length() int {
	return len(seq.Seq)
}

/*SubSeq returns a sub seq. start and end is 1-based.

Examples:

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n
*/
func (seq *Seq) SubSeq(start int, end int) *Seq {
	var newseq *Seq
	start, end, ok := SubLocation(len(seq.Seq), start, end)
	if ok {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, stringutil.Str2Bytes(string(seq.Seq[start-1:end])))
		if len(seq.Qual) > 0 {
			newseq.Qual = stringutil.Str2Bytes(string(seq.Qual[start-1 : end]))
		}
		if len(seq.QualValue) > 0 {
			qv := make([]int, end-start+1)
			for i, v := range seq.QualValue[start-1 : end] {
				qv[i] = v
			}
			newseq.QualValue = qv
		}
	} else {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, []byte(""))
	}

	return newseq
}

// SubSeqInplace return subseq inplace
func (seq *Seq) SubSeqInplace(start int, end int) *Seq {
	var newseq *Seq
	start, end, ok := SubLocation(len(seq.Seq), start, end)
	if ok {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, seq.Seq[start-1:end])
		if len(seq.Qual) > 0 {
			newseq.Qual = seq.Qual[start-1 : end]
		}
		if len(seq.QualValue) > 0 {
			newseq.QualValue = seq.QualValue[start-1 : end]
		}
	} else {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, []byte(""))
	}

	return newseq
}

/*SubLocation is my sublocation strategy,
start, end and returned start and end are all 1-based

 1-based index    1 2 3 4 5 6 7 8 9 10
negative index    0-9-8-7-6-5-4-3-2-1
           seq    A C G T N a c g t n
           1:1    A
           2:4      C G T
         -4:-2                c g t
         -4:-1                c g t n
         -1:-1                      n
          2:-2      C G T N a c g t
          1:-1    A C G T N a c g t n

*/
func SubLocation(length, start, end int) (int, int, bool) {
	if start < 1 {
		if start == 0 {
			start = 1
		} else if start < 0 {
			if end < 0 && start > end {
				return start, end, false
			}
		}
		start = length + start + 1
	}
	if start > length {
		return start, end, false
	}

	if end > length {
		end = length
	}
	if end < 1 {
		if end == 0 {
			end = -1
		}
		end = length + end + 1
	}

	if start-1 > end {
		return start - 1, end, false
	}
	return start, end, true
}

// RemoveGaps remove gaps in place
func (seq *Seq) RemoveGaps(letters string) *Seq {
	if len(letters) == 0 {
		newseq, _ := NewSeqWithQualWithoutValidate(seq.Alphabet, stringutil.Str2Bytes(string(seq.Seq)), stringutil.Str2Bytes(string(seq.Qual)))
		return newseq
	}

	// do not use map
	querySlice := make([]byte, 256)
	for i := 0; i < len(letters); i++ {
		querySlice[int(letters[i])] = letters[i]
	}

	s := make([]byte, len(seq.Seq))
	q := make([]byte, len(seq.Qual))
	var b, g byte
	var j int
	for i := 0; i < len(seq.Seq); i++ {
		b = seq.Seq[i]

		g = querySlice[int(b)]
		if g == 0 { // not gap
			s[j] = b
			if len(seq.Qual) > 0 {
				q[j] = seq.Qual[i]
			}
			j++
		}
	}
	var newSeq *Seq
	if len(seq.Qual) > 0 {
		newSeq, _ = NewSeqWithQualWithoutValidate(seq.Alphabet, s[0:j], q[0:j])
	} else {
		newSeq, _ = NewSeqWithoutValidate(seq.Alphabet, s[0:j])
	}
	return newSeq
}

// RevCom returns reverse complement sequence
func (seq *Seq) RevCom() *Seq {
	return seq.Reverse().Complement()
}

// RevComInplace reverses complement sequence in place
func (seq *Seq) RevComInplace() *Seq {
	return seq.ReverseInplace().ComplementInplace()
}

// Reverse a sequence
func (seq *Seq) Reverse() *Seq {
	if len(seq.Qual) > 0 {
		s := byteutil.ReverseByteSlice(seq.Seq)
		newseq, _ := NewSeqWithQualWithoutValidate(seq.Alphabet, s, byteutil.ReverseByteSlice(seq.Qual))
		return newseq
	}
	s := byteutil.ReverseByteSlice(seq.Seq)
	newseq, _ := NewSeqWithoutValidate(seq.Alphabet, s)
	return newseq
}

// ReverseInplace reverses the sequence content
func (seq *Seq) ReverseInplace() *Seq {
	if len(seq.Qual) > 0 {
		byteutil.ReverseByteSliceInplace(seq.Qual)
	}
	byteutil.ReverseByteSliceInplace(seq.Seq)
	return seq
}

// ComplementSeqLenThreshold is the threshold of sequence length that
// needed to  parallelly complement sequence
var ComplementSeqLenThreshold = 1000

// ComplementThreads is the threads number of parallelly complement sequence
var ComplementThreads = runtime.NumCPU()

// Complement returns complement sequence.
func (seq *Seq) Complement() *Seq {
	var newseq *Seq
	if seq.Alphabet == Unlimit {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, stringutil.Str2Bytes(string(seq.Seq)))
		return newseq
	}

	s := stringutil.Str2Bytes(string(seq.Seq))
	if len(seq.Qual) > 0 {
		newseq, _ = NewSeqWithQualWithoutValidate(seq.Alphabet, s, stringutil.Str2Bytes(string(seq.Qual)))
	} else {
		newseq, _ = NewSeqWithoutValidate(seq.Alphabet, s)
	}

	newseq = newseq.ComplementInplace()
	return newseq
}

// ComplementInplace returns complement sequence.
func (seq *Seq) ComplementInplace() *Seq {
	if seq.Alphabet == Unlimit {
		return seq
	}

	l := len(seq.Seq)
	if l < ComplementSeqLenThreshold {
		var p byte
		for i := 0; i < len(seq.Seq); i++ {
			p, _ = seq.Alphabet.PairLetter(seq.Seq[i])
			seq.Seq[i] = p
		}
		return seq
	}

	chunkSize, start, end := int(l/ComplementThreads), 0, 0
	var wg sync.WaitGroup
	tokens := make(chan int, ComplementThreads)
	for i := 0; i < ComplementThreads; i++ {
		start = i * chunkSize
		end = (i + 1) * chunkSize
		if end > l {
			end = l
		}
		tokens <- 1
		wg.Add(1)

		go func(alphabet *Alphabet, start, end int) {
			defer func() {
				<-tokens
				wg.Done()
			}()

			var p byte
			for i := start; i < end; i++ {
				p = alphabet.pairLetters[seq.Seq[i]-'\x00']
				if p != 0 {
					seq.Seq[i] = p
				}
			}
		}(seq.Alphabet, start, end)
	}
	wg.Wait()

	return seq
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

// GC returns the GC content
func (seq *Seq) GC() float64 {
	return seq.BaseContent("gc")
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
