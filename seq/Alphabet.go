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

IUPAC amino acid code: `ACGTRYSWKMBDHV`

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
	P	Pro	Proline
	Q	Gln	Glutamine
	R	Arg	Arginine
	S	Ser	Serine
	T	Thr	Threonine
	V	Val	Valine
	W	Trp	Tryptophan
	Y	Tyr	Tyrosine
	Z	Glx	Glutamine or Glutamic acid [2]

Other links:

	1. http://www.bioinformatics.org/sms/iupac.html
	2. http://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
	3. http://www.bioinformatics.org/sms2/iupac.html
	4. http://www.matrixscience.com/blog/non-standard-amino-acid-residues.html

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
		a.pairLetters[int(letters[i])] = pairs[i]
	}
	for _, v := range gap {
		a.pairLetters[int(v)] = v
	}
	for _, v := range ambiguous {
		a.pairLetters[int(v)] = v
	}

	return a, nil
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
var ValidSeqLengthThreshold = 100000

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
	if l < ValidSeqLengthThreshold {
		for _, b := range s {
			if !a.IsValidLetter(b) {
				return fmt.Errorf("invalid %s letter: %s", a, []byte{b})
			}
		}
	} else if ValidateWholeSeq {
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

				for i := start; i < end; i++ {
					if !a.IsValidLetter(s[i]) {
						ch <- seqCheckStatus{fmt.Errorf("invalid %s lebtter: %s at %d", a, []byte{s[i]}, i)}
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
	}
	return nil
}

// PairLetter return the Pair Letter
func (a *Alphabet) PairLetter(b byte) (byte, error) {
	if a.isUnlimit {
		return b, nil
	}

	p := a.pairLetters[int(b)]
	if p == 0 {
		return b, fmt.Errorf("invalid letter: %c", b)
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
		[]byte("acgtACGT."),
		[]byte("tgcaTGCA."),
		[]byte(" -"),
		[]byte("nxNX"))

	DNAredundant, _ = NewAlphabet(
		"DNAredundant",
		false,
		[]byte("acgtryswkmbdhvACGTRYSWKMBDHV."),
		[]byte("tgcayrswmkvhdbTGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	RNA, _ = NewAlphabet(
		"RNA",
		false,
		[]byte("acguACGU"),
		[]byte("ugcaUGCA"),
		[]byte(" -"),
		[]byte("nxNX"))

	RNAredundant, _ = NewAlphabet(
		"RNAredundant",
		false,
		[]byte("acguryswkmbdhvACGURYSWKMBDHV."),
		[]byte("ugcayrswmkvhdbUGCAYRSWMKVHDB."),
		[]byte(" -"),
		[]byte("nxNX"))

	Protein, _ = NewAlphabet(
		"Protein",
		false,
		[]byte("abcdefghijklmnpqrstvwyzABCDEFGHIJKLMNPQRSTVWYZ*_."),
		[]byte("abcdefghijklmnpqrstvwyzABCDEFGHIJKLMNPQRSTVWYZ*_."),
		[]byte(" -"),
		[]byte("xX"))

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

// AlphabetGuessSeqLenghtThreshold is the length of sequence prefix of the first FASTA record
// based which FastaRecord guesses the sequence type. 0 for whole seq
var AlphabetGuessSeqLenghtThreshold = 10000

// GuessAlphabet guesses alphabet by given
func GuessAlphabet(seqs []byte) *Alphabet {
	var alphabetMap map[byte]bool
	if AlphabetGuessSeqLenghtThreshold == 0 || len(seqs) <= AlphabetGuessSeqLenghtThreshold {
		alphabetMap = slice2map(byteutil.Alphabet(seqs))
	} else { // reduce guessing time
		alphabetMap = slice2map(byteutil.Alphabet(seqs[0:AlphabetGuessSeqLenghtThreshold]))
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
