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

	M	A/C	K
	R	A/G	Y
	W	A/T	W
	S	C/G	S
	Y	C/T	R
	K	G/T	M

	V	A/C/G	B
	H	A/C/T	D
	D	A/G/T	H
	B	C/G/T	V

	X/N	A/C/G/T	X
	.	not A/C/G/T
	 or-	gap

IUPAC amino acid code

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
	O		pyrrolysine [6]
	P	Pro	Proline
	Q	Gln	Glutamine
	R	Arg	Arginine
	S	Ser	Serine
	T	Thr	Threonine
	U	Sec	selenocysteine [5,6]
	V	Val	Valine
	W	Trp	Tryptophan
	Y	Tyr	Tyrosine
	Z	Glx	Glutamine or Glutamic acid [2]

	X	unknown amino acid
	.	gaps
	*	End

Reference:

	1. http://www.bioinformatics.org/sms/iupac.html
	2. http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	3. http://www.bioinformatics.org/sms2/iupac.html
	4. http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html
	5. http://www.sbcs.qmul.ac.uk/iupac/AminoAcid/A2021.html#AA21
	6. https://en.wikipedia.org/wiki/Amino_acid

*/
package seq

import (
	"errors"
	"fmt"
	"runtime"
	"sync"

	"github.com/shenwei356/util/byteutil"
)

/*Alphabet could be defined. Attention that,
**the letters are case sensitive**.

For example, DNA:

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

	allLetters []byte

	pairLetters []byte
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

	a := &Alphabet{t, isUnlimit, letters, pairs, gap, ambiguous, []byte{}, []byte{}}

	if isUnlimit {
		return a, nil
	}

	if len(letters) != len(pairs) {
		return a, errors.New("mismatch of length of letters and pairs")
	}

	for i := 0; i < len(letters); i++ {
		a.allLetters = append(a.allLetters, letters[i])
	}
	// add gap and ambiguous code
	for _, v := range gap {
		a.allLetters = append(a.allLetters, v)
	}
	for _, v := range ambiguous {
		a.allLetters = append(a.allLetters, v)
	}

	// construct special slice.
	// index are the integer of a byte, and value is the original byte.
	// it's faster than map!!!!
	max := -1
	for i := 0; i < len(a.allLetters); i++ {
		b := int(a.allLetters[i])
		if max < b {
			max = b
		}
	}
	a.pairLetters = make([]byte, max+1)
	for i := 0; i < len(letters); i++ {
		a.pairLetters[letters[i]-'\x00'] = pairs[i]
	}
	for _, v := range gap {
		a.pairLetters[v-'\x00'] = v
	}
	for _, v := range ambiguous {
		a.pairLetters[v-'\x00'] = v
	}

	return a, nil
}

// Clone of a Alphabet
func (a *Alphabet) Clone() *Alphabet {
	return &Alphabet{
		t:           a.t,
		isUnlimit:   a.isUnlimit,
		letters:     []byte(string(a.letters)),
		pairs:       []byte(string(a.pairs)),
		gap:         []byte(string(a.gap)),
		ambiguous:   []byte(string(a.ambiguous)),
		allLetters:  []byte(string(a.allLetters)),
		pairLetters: []byte(string(a.pairLetters)),
	}

}

// Type returns type of the alphabet
func (a *Alphabet) Type() string {
	return a.t
}

// Letters returns letters
func (a *Alphabet) Letters() []byte {
	return a.letters
}

// Gaps returns gaps
func (a *Alphabet) Gaps() []byte {
	return a.gap
}

// AmbiguousLetters returns AmbiguousLetters
func (a *Alphabet) AmbiguousLetters() []byte {
	return a.ambiguous
}

// AllLetters return all letters
func (a *Alphabet) AllLetters() []byte {
	return a.allLetters
}

// String returns type of the alphabet
func (a *Alphabet) String() string {
	return a.t
}

// IsValidLetter is used to validate a letter
func (a *Alphabet) IsValidLetter(b byte) bool {
	if a.isUnlimit {
		return true
	}
	i := int(b)
	if i >= len(a.pairLetters) {
		return false
	}
	return a.pairLetters[i] != 0
}

// ValidSeqLengthThreshold is the threshold of sequence length that
// needed to  parallelly checking sequence
var ValidSeqLengthThreshold = 10000

// ValidateWholeSeq is used to determin whether validate all bases of a seq
var ValidateWholeSeq = true

// ValidSeqThreads is the threads number of parallelly checking sequence
var ValidSeqThreads = runtime.NumCPU()

type seqCheckStatus struct {
	err error
}

// IsValid is used to validate a byte slice
func (a *Alphabet) IsValid(s []byte) error {
	if len(s) == 0 {
		return nil
	}
	if a == nil || a.isUnlimit {
		return nil
	}

	l := len(s)
	var i int
	if l < ValidSeqLengthThreshold {
		for _, b := range s {
			i = int(b)
			if i >= len(a.pairLetters) || a.pairLetters[i] == 0 {
				return fmt.Errorf("seq: invalid %s letter: %s", a, []byte{b})
			}
		}
		return nil
	}

	if ValidateWholeSeq || ValidSeqThreads == 0 {
		ValidSeqThreads = len(s)
	}
	chunkSize, start, end := int(l/ValidSeqThreads), 0, 0

	var wg sync.WaitGroup
	tokens := make(chan int, ValidSeqThreads)
	ch := make(chan seqCheckStatus, ValidSeqThreads)
	done := make(chan struct{})
	finished := false
	for i := 0; i < ValidSeqThreads; i++ {
		start = i * chunkSize
		end = (i + 1) * chunkSize
		if end > l {
			end = l
		}
		tokens <- 1
		wg.Add(1)
		go func(start, end int) {
			defer func() {
				<-tokens
				wg.Done()
			}()

			select {
			case <-done:
				if !finished {
					finished = true
					close(ch)
					return
				}
			default:

			}

			var j int
			for i := start; i < end; i++ {
				j = int(s[i])
				if j >= len(a.pairLetters) || a.pairLetters[j] == 0 {
					ch <- seqCheckStatus{fmt.Errorf("seq: invalid %s lebtter: %s at %d", a, []byte{s[i]}, i)}
					close(done)
					return
				}
			}
			ch <- seqCheckStatus{nil}
		}(start, end)
	}
	wg.Wait()
	close(ch)
	for status := range ch {
		if status.err != nil {
			return status.err
		}
	}

	return nil
}

// PairLetter return the Pair Letter
func (a *Alphabet) PairLetter(b byte) (byte, error) {
	if a.isUnlimit {
		return b, nil
	}
	if int(b) >= len(a.pairLetters) {
		return b, fmt.Errorf("seq: invalid letter: %c", b)
	}
	p := a.pairLetters[b-'\x00']
	if p == 0 {
		return b, fmt.Errorf("seq: invalid letter: %c", b)
	}
	return p, nil
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

	abProtein      map[byte]bool
	abDNAredundant map[byte]bool
	abDNA          map[byte]bool
	abRNAredundant map[byte]bool
	abRNA          map[byte]bool
)

func init() {
	DNA, _ = NewAlphabet(
		"DNA",
		false,
		[]byte("acgtACGT"),
		[]byte("tgcaTGCA"),
		[]byte(" -."),
		[]byte("nN."))

	DNAredundant, _ = NewAlphabet(
		"DNAredundant",
		false,
		[]byte("acgtryswkmbdhvACGTRYSWKMBDHV"),
		[]byte("tgcayrswmkvhdbTGCAYRSWMKVHDB"),
		[]byte(" -."),
		[]byte("nN."))

	RNA, _ = NewAlphabet(
		"RNA",
		false,
		[]byte("acguACGU"),
		[]byte("ugcaUGCA"),
		[]byte(" -."),
		[]byte("nN"))

	RNAredundant, _ = NewAlphabet(
		"RNAredundant",
		false,
		[]byte("acguryswkmbdhvACGURYSWKMBDHV"),
		[]byte("ugcayrswmkvhdbUGCAYRSWMKVHDB"),
		[]byte(" -."),
		[]byte("nN"))

	Protein, _ = NewAlphabet(
		"Protein",
		false,
		[]byte("abcdefghijklmnopqrstuvwyzABCDEFGHIJKLMNOPQRSTUVWYZ"),
		[]byte("abcdefghijklmnopqrstuvwyzABCDEFGHIJKLMNOPQRSTUVWYZ"),
		[]byte(" -"),
		[]byte("xX*_."))

	Unlimit, _ = NewAlphabet(
		"Unlimit",
		true,
		nil,
		nil,
		nil,
		nil)

	abProtein = slice2map(byteutil.Alphabet(Protein.AllLetters()))
	abDNAredundant = slice2map(byteutil.Alphabet(DNAredundant.AllLetters()))
	abDNA = slice2map(byteutil.Alphabet(DNA.AllLetters()))
	abRNAredundant = slice2map(byteutil.Alphabet(RNAredundant.AllLetters()))
	abRNA = slice2map(byteutil.Alphabet(RNA.AllLetters()))
}

// AlphabetGuessSeqLengthThreshold is the length of sequence prefix of the first FASTA record
// based which FastaRecord guesses the sequence type. 0 for whole seq
var AlphabetGuessSeqLengthThreshold = 10000

// GuessAlphabet guesses alphabet by given
func GuessAlphabet(seqs []byte) *Alphabet {
	if len(seqs) == 0 {
		return Unlimit
	}
	var alphabetMap map[byte]bool
	if AlphabetGuessSeqLengthThreshold == 0 || len(seqs) <= AlphabetGuessSeqLengthThreshold {
		alphabetMap = slice2map(byteutil.Alphabet(seqs))
	} else { // reduce guessing time
		alphabetMap = slice2map(byteutil.Alphabet(seqs[0:AlphabetGuessSeqLengthThreshold]))
	}
	if isSubset(alphabetMap, abDNA) {
		return DNA
	}
	if isSubset(alphabetMap, abRNA) {
		return RNA
	}
	if isSubset(alphabetMap, abDNAredundant) {
		return DNAredundant
	}
	if isSubset(alphabetMap, abRNAredundant) {
		return RNAredundant
	}
	if isSubset(alphabetMap, abProtein) {
		return Protein
	}

	return Unlimit
}

// GuessAlphabetLessConservatively change DNA to DNAredundant and RNA to RNAredundant
func GuessAlphabetLessConservatively(seqs []byte) *Alphabet {
	ab := GuessAlphabet(seqs)
	if ab == DNA {
		return DNAredundant
	}
	if ab == RNA {
		return RNAredundant
	}
	return ab
}

func isSubset(query, subject map[byte]bool) bool {
	for b := range query {
		if _, ok := subject[b]; !ok {
			return false
		}
	}
	return true
}

func slice2map(s []byte) map[byte]bool {
	m := make(map[byte]bool)
	for _, b := range s {
		m[b] = true
	}
	return m
}
