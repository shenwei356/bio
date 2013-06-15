package seq

import (
	"bytes"
	"fmt"
	"regexp"
	"strings"
)

type SeqType int

const (
	DNA = iota
	RNA = iota
)

var (
	BasesDNA     []byte
	BasesRNA     []byte
	BasesFull    []byte
	BasesDNA2Com map[byte]byte
	BasesRNA2Com map[byte]byte
	BasesFullMap map[byte][]byte
)

func init() {
	BasesDNA = []byte{'A', 'C', 'G', 'T'}
	BasesRNA = []byte{'A', 'C', 'G', 'U'}
	BasesFull = []byte{'A', 'C', 'G', 'T', 'U', 'R', 'Y', 'M', 'K', 'S', 'W',
		'H', 'B', 'V', 'D', 'N'}
	BasesDNA2Com = map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	BasesRNA2Com = map[byte]byte{'A': 'U', 'U': 'A', 'C': 'G', 'G': 'C'}
	BasesFullMap = map[byte][]byte{'A': {'A'}, 'T': {'T'}, 'C': {'C'},
		'G': {'G'}, 'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'M': {'A', 'C'},
		'K': {'G', 'T'}, 'S': {'C', 'G'}, 'W': {'A', 'T'},
		'H': {'A', 'C', 'T'}, 'B': {'C', 'G', 'T'}, 'V': {'A', 'C', 'G'},
		'D': {'A', 'T', 'G'}, 'N': {'A', 'C', 'G', 'T'}, 'U': {'U'}}

}

/*
Validate nucleotide sequence.
*/
func ValidateSequence(seq string) bool {
	var buffer bytes.Buffer
	for _, v := range BasesFull {
		buffer.WriteByte(v)
	}
	validSeq, err := regexp.Compile(fmt.Sprintf("(?i)[^%s]", buffer.String()))
	if err != nil {
		recover()
	}
	buffer.Reset()
	return !validSeq.MatchString(seq)
}

/*
Return the base content of a nucleotide sequence.
*/
func BaseContent(seq string, base string) float32 {
	up := strings.ToUpper(base)
	lo := strings.ToLower(base)
	return float32((strings.Count(seq, up) +
		strings.Count(seq, lo))) / float32(len(seq))
}

/*
Return the reverse-complement counterpart of a nucleotide sequence.
*/
func Revcom(seq string, seqtype SeqType) string {
	var Bases2Com map[byte]byte
	switch seqtype {
	case DNA:
		Bases2Com = BasesDNA2Com
	case RNA:
		Bases2Com = BasesRNA2Com
	}
	seq = strings.ToUpper(seq)
	var buffer bytes.Buffer
	for _, v := range []byte(seq) {
		buffer.WriteByte(Bases2Com[v])
	}
	return Reverse(buffer.String())
}

/*
Reverse a string
*/
func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}
