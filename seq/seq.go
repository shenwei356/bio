package seq

import (
	"bytes"
	"fmt"
	"math"
	"runtime"
	"strings"
	"sync"

	"github.com/shenwei356/util/byteutil"
)

var defaultBytesBufferSize = 10 << 20

var bufferedByteSliceWrapper *byteutil.BufferedByteSliceWrapper

var QUAL_MAP [256]float64

func initQualMap() {
	for i, _ := range QUAL_MAP {
		QUAL_MAP[i] = math.Pow(10, float64(i)/-10)
	}
	QUAL_MAP[255] = 1.0

}

func init() {
	bufferedByteSliceWrapper = byteutil.NewBufferedByteSliceWrapper(1, defaultBytesBufferSize)
	initQualMap()
}

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
		return nil, fmt.Errorf("seq: unmatched length of sequence (%d) and quality (%d)", len(s), len(q))
	}
	seq, err := NewSeq(t, s)
	if err != nil {
		return nil, err
	}
	seq.Qual = q
	return seq, nil
}

// NewSeqWithoutValidation create Seq without check the sequences
func NewSeqWithoutValidation(t *Alphabet, s []byte) (*Seq, error) {
	seq := &Seq{Alphabet: t, Seq: s}
	return seq, nil
}

// NewSeqWithQualWithoutValidation create Seq with quality without check the sequences
func NewSeqWithQualWithoutValidation(t *Alphabet, s []byte, q []byte) (*Seq, error) {
	if len(s) != len(q) {
		return nil, fmt.Errorf("seq: unmatched length of sequence (%d) and quality (%d)", len(s), len(q))
	}
	seq := &Seq{Alphabet: t, Seq: s, Qual: q}
	return seq, nil
}

// Length returns the length of sequence
func (seq *Seq) Length() int {
	return len(seq.Seq)
}

func (seq *Seq) String() string {
	return fmt.Sprintf("%s, len:%d, seq:%s, qual:%s", seq.Alphabet.String(), len(seq.Seq), seq.Seq, seq.Qual)
}

// Clone of a Seq
func (seq *Seq) Clone() *Seq {
	s := make([]byte, len(seq.Seq))
	copy(s, seq.Seq)

	q := make([]byte, len(seq.Qual))
	copy(q, seq.Qual)

	qv := make([]int, len(seq.QualValue))
	copy(qv, seq.QualValue)

	return &Seq{
		Alphabet:  seq.Alphabet.Clone(),
		Seq:       s,
		Qual:      q,
		QualValue: qv,
	}
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
		  1:12    A C G T N a c g t n
		-12:-1    A C G T N a c g t n
*/
func (seq *Seq) SubSeq(start int, end int) *Seq {
	var newseq *Seq
	start, end, ok := SubLocation(len(seq.Seq), start, end)
	if ok {
		s := make([]byte, end-start+1)
		copy(s, seq.Seq[start-1:end])
		newseq, _ = NewSeqWithoutValidation(seq.Alphabet, s)
		if len(seq.Qual) > 0 {
			newseq.Qual = []byte(string(seq.Qual[start-1 : end]))
		}
		if len(seq.QualValue) > 0 {
			qv := make([]int, end-start+1)
			for i, v := range seq.QualValue[start-1 : end] {
				qv[i] = v
			}
			newseq.QualValue = qv
		}
	} else {
		newseq, _ = NewSeqWithoutValidation(seq.Alphabet, []byte(""))
	}

	return newseq
}

// SubSeqInplace return subseq inplace
func (seq *Seq) SubSeqInplace(start int, end int) *Seq {
	start, end, ok := SubLocation(len(seq.Seq), start, end)
	if ok {
		if len(seq.Seq) > 0 {
			seq.Seq = seq.Seq[start-1 : end]
		}
		if len(seq.Qual) > 0 {
			seq.Qual = seq.Qual[start-1 : end]
		}
		if len(seq.QualValue) > 0 {
			seq.QualValue = seq.QualValue[start-1 : end]
		}
	} else {
		seq.Seq = []byte{}
		seq.Qual = []byte{}
		seq.QualValue = []int{}
	}

	return seq
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
		  1:12    A C G T N a c g t n
		-12:-1    A C G T N a c g t n

*/
func SubLocation(length, start, end int) (int, int, bool) {
	if length == 0 {
		return 0, 0, false
	}
	if start < 1 {
		if start == 0 {
			start = 1
		} else if start < 0 {
			if end < 0 && start > end {
				return start, end, false
			}

			if -start > length {
				return start, end, false
			}
			start = length + start + 1
		}
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

// RemoveGaps return a new seq without gaps
func (seq *Seq) RemoveGaps(letters string) *Seq {
	if len(letters) == 0 {
		return seq.Clone()
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
		newSeq, _ = NewSeqWithQualWithoutValidation(seq.Alphabet, s[0:j], q[0:j])
	} else {
		newSeq, _ = NewSeqWithoutValidation(seq.Alphabet, s[0:j])
	}
	return newSeq
}

// RemoveGapsInplace removes gaps in place
func (seq *Seq) RemoveGapsInplace(letters string) *Seq {
	if len(letters) == 0 {
		return seq
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
	seq.Seq = s[0:j]
	if len(seq.Qual) > 0 {
		seq.Qual = q[0:j]
	}
	return seq
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
	return seq.Clone().ReverseInplace()
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
	return seq.Clone().ComplementInplace()
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
		if i == ComplementThreads-1 && end < l {
			end = l
		}

		if end > l {
			end = l
		}
		tokens <- 1
		wg.Add(1)

		go func(alphabet *Alphabet, start, end int) {
			var p byte
			for i := start; i < end; i++ {
				p, _ = seq.Alphabet.PairLetter(seq.Seq[i])
				seq.Seq[i] = p
			}
			<-tokens
			wg.Done()
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

// BaseContentCaseSensitive returns base content for given case sensitive bases.
func (seq *Seq) BaseContentCaseSensitive(list string) float64 {
	if len(seq.Seq) == 0 {
		return float64(0)
	}

	sum := 0
	for _, b := range []byte(list) {
		sum += bytes.Count(seq.Seq, []byte{b})
	}

	return float64(sum) / float64(len(seq.Seq))
}

// GC returns the GC content
func (seq *Seq) GC() float64 {
	return seq.BaseContent("gcs")
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

// DegenerateBaseMapNucl2 mappings nucleic acid degenerate base to all bases.
var DegenerateBaseMapNucl2 = map[byte]string{
	'A': "A",
	'T': "T",
	'U': "U",
	'C': "C",
	'G': "G",
	'R': "AG",
	'Y': "CT",
	'M': "AC",
	'K': "GT",
	'S': "CG",
	'W': "AT",
	'H': "ACT",
	'B': "CGT",
	'V': "ACG",
	'D': "AGT",
	'N': "ACGT",
	'a': "a",
	't': "t",
	'u': "u",
	'c': "c",
	'g': "g",
	'r': "ag",
	'y': "ct",
	'm': "ac",
	'k': "gt",
	's': "cg",
	'w': "at",
	'h': "act",
	'b': "cgt",
	'v': "acg",
	'd': "agt",
	'n': "acgt",
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
	'X': "[A-Z]",
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
	'x': "[a-z]",
	'y': "y",
	'z': "[qe]",
}

// Degenerate2Regexp transforms seqs containing degenrate base to regular expression
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

// Degenerate2Seqs transforms seqs containing degenrate bases to all possible sequences.
func Degenerate2Seqs(s []byte) (dseqs [][]byte, err error) {
	dseqs = [][]byte{[]byte{}}
	var i, j, k int
	var ok bool
	var dbases string
	var dbase byte
	for _, base := range s {
		if dbases, ok = DegenerateBaseMapNucl2[base]; ok {
			if len(dbases) == 1 {
				dbase = dbases[0]
				for i = 0; i < len(dseqs); i++ {
					dseqs[i] = append(dseqs[i], dbase)
				}
			} else {
				// 2nd
				more := make([][]byte, len(dseqs)*(len(dbases)-1))
				k = 0
				for i = 1; i < len(dbases); i++ {
					for j = 0; j < len(dseqs); j++ {
						more[k] = []byte(string(append(dseqs[j], dbases[i])))
						k++
					}
				}

				// 1th
				for i = 0; i < len(dseqs); i++ {
					dseqs[i] = append(dseqs[i], dbases[0])
				}

				dseqs = append(dseqs, more...)
			}

		} else {
			return dseqs, fmt.Errorf("seq: invalid letter: %c", base)
		}
	}
	return dseqs, nil
}

// Translate translates the RNA/DNA to amino acid sequence.
// Available frame: 1, 2, 3, -1, -2 ,-3.
// If option trim is true, it removes all 'X' and '*' characters from the right end of the translation.
// If option clean is true, it changes all STOP codon positions from the '*' character to 'X' (an unknown residue).
// If option allowUnknownCodon is true, codons not in the codon table will be translated to 'X'.
// If option markInitCodonAsM is true, initial codon at beginning will be represented as 'M'.
func (seq *Seq) Translate(transl_table int, frame int, trim bool, clean bool, allowUnknownCodon bool, markInitCodonAsM bool) (*Seq, error) {
	if !(seq.Alphabet == DNA || seq.Alphabet == DNAredundant || seq.Alphabet == RNA || seq.Alphabet == RNAredundant) {
		return nil, fmt.Errorf("seq: only DNA/RNA sequence can all method Translate, the alphabet is %s", seq.Alphabet.String())
	}
	var codonTable *CodonTable
	var ok bool
	if codonTable, ok = CodonTables[transl_table]; !ok {
		return nil, fmt.Errorf("seq: invalid codon table: %d", transl_table)
	}
	if !(frame == 1 || frame == 2 || frame == 3 || frame == -1 || frame == -2 || frame == -3) {
		return nil, fmt.Errorf("seq: invalid frame: %d. available: 1, 2, 3, -1, -2, -3", frame)
	}

	aa, err := codonTable.Translate(seq.Seq, frame, trim, clean, allowUnknownCodon, markInitCodonAsM)
	if err != nil {
		return nil, err
	}

	t, err := NewSeqWithoutValidation(Protein, aa)
	if err != nil {
		return nil, err
	}
	return t, nil
}

// ParseQual parses sequence quality, asciiBase = 33 for Phred+33.
func (seq *Seq) ParseQual(asciiBase int) {
	if len(seq.Qual) == 0 {
		return
	}
	qv := make([]int, len(seq.Qual))
	for i, q := range seq.Qual {
		qv[i] = int(q) - asciiBase
	}
	seq.QualValue = qv
}

// Calculate average quality value.
func (seq *Seq) AvgQual(asciiBase int) float64 {
	if len(seq.Qual) > 0 && len(seq.QualValue) == 0 {
		seq.ParseQual(asciiBase)
	}
	if len(seq.QualValue) == 0 {
		return 0.0
	}
	var sum float64
	for _, q := range seq.QualValue {
		sum += QUAL_MAP[q]
	}
	return -10 * math.Log10(sum/float64(len(seq.QualValue)))
}
