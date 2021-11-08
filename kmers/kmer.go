// Copyright © 2018-2021 Wei Shen <shenwei356@gmail.com>
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//b
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package kmers

import (
	"bytes"
	"errors"
)

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("kmers: illegal base")

// ErrKOverflow means K > 32.
var ErrKOverflow = errors.New("kmers: k-mer size (1-32) overflow")

// ErrCodeOverflow means the encode interger is bigger than 4^k.
var ErrCodeOverflow = errors.New("kmers: code value overflow")

// ErrKMismatch means K size mismatch.
var ErrKMismatch = errors.New("kmers: K mismatch")

// slice is much faster than switch and map.
var base2bit = [256]uint64{
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 0, 1, 1, 0, 4, 4, 2, 0, 4, 4, 2, 4, 0, 0, 4,
	4, 4, 0, 1, 3, 3, 0, 0, 4, 1, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
	4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
}

// var base2bit []uint64

// MaxCode is the maxinum interger for all Ks.
var MaxCode []uint64

func init() {
	MaxCode = make([]uint64, 33)
	for i := 1; i <= 32; i++ {
		MaxCode[i] = 1<<uint(i*2) - 1
	}
}

// Encode converts byte slice to bits.
//
// Codes:
//
// 	  A    0b00
// 	  C    0b01
// 	  G    0b10
// 	  T    0b11
//
// For degenerate bases, only the first base is kept.
//
//     M       AC     A
//     V       ACG    A
//     H       ACT    A
//     R       AG     A
//     D       AGT    A
//     W       AT     A
//     S       CG     C
//     B       CGT    C
//     Y       CT     C
//     K       GT     G
//     N       ACGT   A
//
func Encode(kmer []byte) (code uint64, err error) {
	if len(kmer) == 0 || len(kmer) > 32 {
		return 0, ErrKOverflow
	}

	var v uint64
	for _, b := range kmer {
		code <<= 2
		v = base2bit[b]
		// if v > 3 {
		if v == 4 {
			return code, ErrIllegalBase
		}
		code |= v
	}
	return code, nil
}

// ErrNotConsecutiveKmers means the two k-mers are not consecutive.
var ErrNotConsecutiveKmers = errors.New("kmers: not consecutive k-mers")

// MustEncodeFromFormerKmer encodes from former the k-mer,
// assuming the k-mer and leftKmer are both OK.
func MustEncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	v := base2bit[kmer[len(kmer)-1]]
	// if v > 3 {
	if v == 4 {
		return leftCode, ErrIllegalBase
	}
	// retrieve (k-1)*2 bits and << 2, and then add v
	return leftCode&((1<<(uint(len(kmer)-1)<<1))-1)<<2 | v, nil
}

// EncodeFromFormerKmer encodes from the former k-mer, inspired by ntHash
func EncodeFromFormerKmer(kmer []byte, leftKmer []byte, leftCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(leftKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(kmer[0:len(kmer)-1], leftKmer[1:]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromFormerKmer(kmer, leftKmer, leftCode)
}

// MustEncodeFromLatterKmer encodes from the latter k-mer,
// assuming the k-mer and rightKmer are both OK.
func MustEncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	v := base2bit[kmer[0]]
	// if v > 3 {
	if v == 4 {
		return rightCode, ErrIllegalBase
	}

	return v<<(uint(len(kmer)-1)<<1) | rightCode>>2, nil
}

// EncodeFromLatterKmer encodes from the former k-mer.
func EncodeFromLatterKmer(kmer []byte, rightKmer []byte, rightCode uint64) (uint64, error) {
	if len(kmer) == 0 {
		return 0, ErrKOverflow
	}
	if len(kmer) != len(rightKmer) {
		return 0, ErrKMismatch
	}
	if !bytes.Equal(rightKmer[0:len(kmer)-1], kmer[1:len(rightKmer)]) {
		return 0, ErrNotConsecutiveKmers
	}
	return MustEncodeFromLatterKmer(kmer, rightKmer, rightCode)
}

// Reverse returns code of the reversed sequence.
func Reverse(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code & 3)
	// 	code >>= 2
	// }
	// return

	// https: //www.biostars.org/p/113640, with a little modification
	c = code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	return (c >> (2 * (32 - k)))
}

// MustReverse is similar to Reverse, but does not check k.
func MustReverse(code uint64, k int) (c uint64) {
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code & 3)
	// 	code >>= 2
	// }
	// return

	// https: //www.biostars.org/p/113640, with a little modification
	c = code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	return (c >> (2 * (32 - k)))
}

// Complement returns code of complement sequence.
func Complement(code uint64, k int) uint64 {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	return code ^ (1<<uint(k<<1) - 1)
}

// MustComplement is similar to Complement, but does not check k.
func MustComplement(code uint64, k int) uint64 {
	return code ^ (1<<uint(k<<1) - 1)
}

// RevComp returns code of reverse complement sequence.
func RevComp(code uint64, k int) (c uint64) {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code&3 ^ 3)
	// 	code >>= 2
	// }
	// return

	// https://www.biostars.org/p/113640/#9474334
	c = ^code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	return (c >> (2 * (32 - k)))
}

// MustRevComp is similar to RevComp, but does not check k.
func MustRevComp(code uint64, k int) (c uint64) {
	// for i := 0; i < k; i++ {
	// 	c = (c << 2) | (code&3 ^ 3)
	// 	code >>= 2
	// }
	// return

	// https://www.biostars.org/p/113640/#9474334
	c = ^code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	return (c >> (2 * (32 - k)))
}

// Canonical returns code of its canonical kmer.
func Canonical(code uint64, k int) uint64 {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}

	var rc uint64
	// c := code
	// for i := 0; i < k; i++ {
	// 	rc = (rc << 2) | (c&3 ^ 3)
	// 	c >>= 2
	// }

	// https://www.biostars.org/p/113640/#9474334
	c := ^code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	rc = (c >> (2 * (32 - k)))

	if rc < code {
		return rc
	}
	return code
}

// MustCanonical is similar to Canonical, but does not check k.
func MustCanonical(code uint64, k int) uint64 {
	var rc uint64
	// c := code
	// for i := 0; i < k; i++ {
	// 	rc = (rc << 2) | (c&3 ^ 3)
	// 	c >>= 2
	// }

	// https://www.biostars.org/p/113640/#9474334
	c := ^code
	c = ((c >> 2 & 0x3333333333333333) | (c&0x3333333333333333)<<2)
	c = ((c >> 4 & 0x0F0F0F0F0F0F0F0F) | (c&0x0F0F0F0F0F0F0F0F)<<4)
	c = ((c >> 8 & 0x00FF00FF00FF00FF) | (c&0x00FF00FF00FF00FF)<<8)
	c = ((c >> 16 & 0x0000FFFF0000FFFF) | (c&0x0000FFFF0000FFFF)<<16)
	c = ((c >> 32 & 0x00000000FFFFFFFF) | (c&0x00000000FFFFFFFF)<<32)
	rc = (c >> (2 * (32 - k)))

	if rc < code {
		return rc
	}
	return code
}

// bit2base is for mapping bit to base.
var bit2base = [4]byte{'A', 'C', 'G', 'T'}

// bit2str is for output bits string
var bit2str = [4]string{"00", "01", "10", "11"}

// Decode converts the code to original seq
func Decode(code uint64, k int) []byte {
	if k <= 0 || k > 32 {
		panic(ErrKOverflow)
	}
	if code > MaxCode[k] {
		panic(ErrCodeOverflow)
	}
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]
		code >>= 2
	}
	return kmer
}

// MustDecode is similar to Decode, but does not check k and code.
func MustDecode(code uint64, k int) []byte {
	kmer := make([]byte, k)
	for i := 0; i < k; i++ {
		kmer[k-1-i] = bit2base[code&3]
		code >>= 2
	}
	return kmer
}

// KmerCode is a struct representing a k-mer in 64-bits.
type KmerCode struct {
	Code uint64
	K    int
}

// NewKmerCode returns a new KmerCode struct from byte slice.
func NewKmerCode(kmer []byte) (KmerCode, error) {
	code, err := Encode(kmer)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeFromFormerOne computes KmerCode from the Former consecutive k-mer.
func NewKmerCodeFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := EncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// NewKmerCodeMustFromFormerOne computes KmerCode from the Former consecutive k-mer,
// assuming the k-mer and leftKmer are both OK.
func NewKmerCodeMustFromFormerOne(kmer []byte, leftKmer []byte, preKcode KmerCode) (KmerCode, error) {
	code, err := MustEncodeFromFormerKmer(kmer, leftKmer, preKcode.Code)
	if err != nil {
		return KmerCode{}, err
	}
	return KmerCode{code, len(kmer)}, err
}

// Equal checks wether two KmerCodes are the same.
func (kcode KmerCode) Equal(kcode2 KmerCode) bool {
	return kcode.K == kcode2.K && kcode.Code == kcode2.Code
}

// Rev returns KmerCode of the reverse sequence.
func (kcode KmerCode) Rev() KmerCode {
	return KmerCode{MustReverse(kcode.Code, kcode.K), kcode.K}
}

// Comp returns KmerCode of the complement sequence.
func (kcode KmerCode) Comp() KmerCode {
	return KmerCode{MustComplement(kcode.Code, kcode.K), kcode.K}
}

// RevComp returns KmerCode of the reverse complement sequence.
func (kcode KmerCode) RevComp() KmerCode {
	return KmerCode{MustRevComp(kcode.Code, kcode.K), kcode.K}
}

// Canonical returns its canonical kmer
func (kcode KmerCode) Canonical() KmerCode {
	rcKcode := kcode.RevComp()
	if rcKcode.Code < kcode.Code {
		return rcKcode
	}
	return kcode
}

// Bytes returns k-mer in []byte.
func (kcode KmerCode) Bytes() []byte {
	return Decode(kcode.Code, kcode.K)
}

// String returns k-mer in string
func (kcode KmerCode) String() string {
	return string(Decode(kcode.Code, kcode.K))
}

// BitsString returns code to string
func (kcode KmerCode) BitsString() string {
	var buf bytes.Buffer
	for _, b := range Decode(kcode.Code, kcode.K) {
		buf.WriteString(bit2str[base2bit[b]])
	}
	return buf.String()
}
