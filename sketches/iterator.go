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
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package sketches

import (
	"errors"
	"fmt"
	"math"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/kmers"
	"github.com/will-rowe/nthash"
)

// ErrInvalidK means k < 1.
var ErrInvalidK = fmt.Errorf("sketches: invalid k-mer size")

// ErrEmptySeq sequence is empty.
var ErrEmptySeq = fmt.Errorf("sketches: empty sequence")

// ErrShortSeq means the sequence is shorter than k
var ErrShortSeq = fmt.Errorf("sketches: sequence too short")

// ErrIllegalBase means that base beyond IUPAC symbols are  detected.
var ErrIllegalBase = errors.New("sketches: illegal base")

// ErrKTooLarge means that the k-mer size is too large.
var ErrKTooLarge = fmt.Errorf("sketches: k-mer size is too large")

// ErrInvalidM means that the m-mer size is too large or too small, should be in range of [4, k].
var ErrInvalidM = fmt.Errorf("sketches: invalid m-mer size, should be in range of [4, k]")

// ErrInvalidScale means
var ErrInvalidScale = fmt.Errorf("sketches: invalid scale, should be in range of [1, k-m+1]")

var poolIterator = &sync.Pool{New: func() interface{} {
	return &Iterator{}
}}

// Iterator is a kmer code (k<=32) or hash iterator.
type Iterator struct {
	s         *seq.Seq // only used for KmerIterator
	k         int
	kUint     uint // uint(k)
	kP1       int  // k -1
	kP1Uint   uint // uint(k-1)
	canonical bool
	circular  bool

	hash bool

	finished     bool
	revcomStrand bool
	idx          int

	// for KmerIterator
	length    int
	end, e    int
	first     bool
	kmer      []byte
	codeBase  uint64
	preCode   uint64
	preCodeRC uint64
	codeRC    uint64

	mask1 uint64 // (1<<(kP1Uint*2))-1
	mask2 uint   // iter.kP1Uint*2

	// for HashIterator
	hasher *nthash.NTHi

	// for SimHash
	simhash  bool
	m        int       // size of m-mer
	em       int       // end of idxJ in current k-mer
	hashes   []uint64  // cycle buffer of previous hashes of m-mer
	sum      [64]int16 // the vector
	nPosHash int16     // number of vector value > 0
	nThsld   int16     // the threshold to decide should a vector value be 1 (x >= nThsld) or not
	idxJ     int       // tmp index
	preHashI int       // index of the m-mer just getting out of the current k-mer
	preHash  uint64    //
	_hash    uint64    // ntHash

	fracMinHash bool   // scale > 1
	maxHash     uint64 // the max hash value for a scale, i.e. max uint64 / scale
}

// NewSimHashIterator returns a SimHash Iterator. Main parameters:
//
//	k:     k-mer size.
//	m:     size of m-mer in a k-mer. range: [4, k]
//	scale: scale of FracMinHash of m-mers. range: [1, k-m+1]
func NewSimHashIterator(s *seq.Seq, k int, m int, scale int, canonical bool, circular bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if k >= 65535 {
		return nil, ErrKTooLarge
	}

	if m < 4 || m > k {
		return nil, ErrInvalidM
	}
	if scale < 1 || scale > k-m+1 {
		return nil, ErrInvalidScale
	}

	if len(s.Seq) < k {
		return nil, ErrShortSeq
	}

	var s2 *seq.Seq
	if circular {
		s2 = s.Clone() // do not edit original sequence
		s2.Seq = append(s2.Seq, s.Seq[0:k-1]...)
	} else {
		s2 = s
	}

	// iter := &Iterator{s: s2, k: k, canonical: canonical, circular: circular}
	iter := poolIterator.Get().(*Iterator)
	iter.s = s2
	iter.k = k
	iter.canonical = canonical
	iter.circular = circular
	iter.finished = false
	iter.revcomStrand = false
	iter.idx = 0

	iter.length = len(s2.Seq)
	iter.end = iter.length - k + 1
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)

	iter.first = true

	var err error
	iter.hasher, err = nthash.NewHasher(&s2.Seq, uint(m))
	if err != nil {
		return nil, err
	}

	for i := 0; i < 64; i++ {
		iter.sum[i] = 0
	}
	iter.nPosHash = 0

	iter.simhash = true
	iter.m = m
	iter.em = iter.k - iter.m
	iter.nThsld = 0                          // int16((float64(iter.k-iter.m+1) + 1) / float64(2))
	if len(iter.hashes) >= iter.k-iter.m+1 { // resuse objects
		iter.hashes = iter.hashes[0 : iter.k-iter.m+1]
	} else {
		iter.hashes = make([]uint64, iter.k-iter.m+1)
	}
	iter.preHashI = 0
	iter.idxJ = 0
	iter.fracMinHash = scale > 1
	if iter.fracMinHash {
		iter.maxHash = math.MaxUint64 / uint64(scale)
	} else {
		iter.maxHash = math.MaxUint64
	}

	return iter, nil
}

// NextSimHash returns next SimHash.
func (iter *Iterator) NextSimHash() (code uint64, ok bool) {
	if iter.finished {
		return 0, false
	}

	if iter.idx == iter.end {
		iter.finished = true
		poolIterator.Put(iter)
		return 0, false
	}

	if !iter.first {
		// subtract the previous  hash, which is out of the window of current k-mer

		iter.preHash = iter.hashes[iter.preHashI]
		// fmt.Printf("[%d] prevHashI: %d\n", iter.idx, iter.preHashI)

		if iter.preHash > 0 {
			iter.nPosHash--

			iter.sum[0] -= int16(iter.preHash >> 63 & 1)
			iter.sum[1] -= int16(iter.preHash >> 62 & 1)
			iter.sum[2] -= int16(iter.preHash >> 61 & 1)
			iter.sum[3] -= int16(iter.preHash >> 60 & 1)
			iter.sum[4] -= int16(iter.preHash >> 59 & 1)
			iter.sum[5] -= int16(iter.preHash >> 58 & 1)
			iter.sum[6] -= int16(iter.preHash >> 57 & 1)
			iter.sum[7] -= int16(iter.preHash >> 56 & 1)
			iter.sum[8] -= int16(iter.preHash >> 55 & 1)
			iter.sum[9] -= int16(iter.preHash >> 54 & 1)
			iter.sum[10] -= int16(iter.preHash >> 53 & 1)
			iter.sum[11] -= int16(iter.preHash >> 52 & 1)
			iter.sum[12] -= int16(iter.preHash >> 51 & 1)
			iter.sum[13] -= int16(iter.preHash >> 50 & 1)
			iter.sum[14] -= int16(iter.preHash >> 49 & 1)
			iter.sum[15] -= int16(iter.preHash >> 48 & 1)
			iter.sum[16] -= int16(iter.preHash >> 47 & 1)
			iter.sum[17] -= int16(iter.preHash >> 46 & 1)
			iter.sum[18] -= int16(iter.preHash >> 45 & 1)
			iter.sum[19] -= int16(iter.preHash >> 44 & 1)
			iter.sum[20] -= int16(iter.preHash >> 43 & 1)
			iter.sum[21] -= int16(iter.preHash >> 42 & 1)
			iter.sum[22] -= int16(iter.preHash >> 41 & 1)
			iter.sum[23] -= int16(iter.preHash >> 40 & 1)
			iter.sum[24] -= int16(iter.preHash >> 39 & 1)
			iter.sum[25] -= int16(iter.preHash >> 38 & 1)
			iter.sum[26] -= int16(iter.preHash >> 37 & 1)
			iter.sum[27] -= int16(iter.preHash >> 36 & 1)
			iter.sum[28] -= int16(iter.preHash >> 35 & 1)
			iter.sum[29] -= int16(iter.preHash >> 34 & 1)
			iter.sum[30] -= int16(iter.preHash >> 33 & 1)
			iter.sum[31] -= int16(iter.preHash >> 32 & 1)
			iter.sum[32] -= int16(iter.preHash >> 31 & 1)
			iter.sum[33] -= int16(iter.preHash >> 30 & 1)
			iter.sum[34] -= int16(iter.preHash >> 29 & 1)
			iter.sum[35] -= int16(iter.preHash >> 28 & 1)
			iter.sum[36] -= int16(iter.preHash >> 27 & 1)
			iter.sum[37] -= int16(iter.preHash >> 26 & 1)
			iter.sum[38] -= int16(iter.preHash >> 25 & 1)
			iter.sum[39] -= int16(iter.preHash >> 24 & 1)
			iter.sum[40] -= int16(iter.preHash >> 23 & 1)
			iter.sum[41] -= int16(iter.preHash >> 22 & 1)
			iter.sum[42] -= int16(iter.preHash >> 21 & 1)
			iter.sum[43] -= int16(iter.preHash >> 20 & 1)
			iter.sum[44] -= int16(iter.preHash >> 19 & 1)
			iter.sum[45] -= int16(iter.preHash >> 18 & 1)
			iter.sum[46] -= int16(iter.preHash >> 17 & 1)
			iter.sum[47] -= int16(iter.preHash >> 16 & 1)
			iter.sum[48] -= int16(iter.preHash >> 15 & 1)
			iter.sum[49] -= int16(iter.preHash >> 14 & 1)
			iter.sum[50] -= int16(iter.preHash >> 13 & 1)
			iter.sum[51] -= int16(iter.preHash >> 12 & 1)
			iter.sum[52] -= int16(iter.preHash >> 11 & 1)
			iter.sum[53] -= int16(iter.preHash >> 10 & 1)
			iter.sum[54] -= int16(iter.preHash >> 9 & 1)
			iter.sum[55] -= int16(iter.preHash >> 8 & 1)
			iter.sum[56] -= int16(iter.preHash >> 7 & 1)
			iter.sum[57] -= int16(iter.preHash >> 6 & 1)
			iter.sum[58] -= int16(iter.preHash >> 5 & 1)
			iter.sum[59] -= int16(iter.preHash >> 4 & 1)
			iter.sum[60] -= int16(iter.preHash >> 3 & 1)
			iter.sum[61] -= int16(iter.preHash >> 2 & 1)
			iter.sum[62] -= int16(iter.preHash >> 1 & 1)
			iter.sum[63] -= int16(iter.preHash & 1)
		}

		// add newly added m-mer

		iter._hash, _ = iter.hasher.Next(iter.canonical)

		if iter.fracMinHash && iter._hash > iter.maxHash { // discard it, see sourmash paper for FacMinHash
			iter._hash = 0
		} else if iter._hash > 0 {
			iter.nPosHash++
		}
		// fmt.Printf("[%d] new: %064b\n", iter.idx, iter._hash)

		iter.hashes[iter.preHashI] = iter._hash // update the hash value

		if iter._hash > 0 {
			iter.sum[0] += int16(iter._hash >> 63 & 1)
			iter.sum[1] += int16(iter._hash >> 62 & 1)
			iter.sum[2] += int16(iter._hash >> 61 & 1)
			iter.sum[3] += int16(iter._hash >> 60 & 1)
			iter.sum[4] += int16(iter._hash >> 59 & 1)
			iter.sum[5] += int16(iter._hash >> 58 & 1)
			iter.sum[6] += int16(iter._hash >> 57 & 1)
			iter.sum[7] += int16(iter._hash >> 56 & 1)
			iter.sum[8] += int16(iter._hash >> 55 & 1)
			iter.sum[9] += int16(iter._hash >> 54 & 1)
			iter.sum[10] += int16(iter._hash >> 53 & 1)
			iter.sum[11] += int16(iter._hash >> 52 & 1)
			iter.sum[12] += int16(iter._hash >> 51 & 1)
			iter.sum[13] += int16(iter._hash >> 50 & 1)
			iter.sum[14] += int16(iter._hash >> 49 & 1)
			iter.sum[15] += int16(iter._hash >> 48 & 1)
			iter.sum[16] += int16(iter._hash >> 47 & 1)
			iter.sum[17] += int16(iter._hash >> 46 & 1)
			iter.sum[18] += int16(iter._hash >> 45 & 1)
			iter.sum[19] += int16(iter._hash >> 44 & 1)
			iter.sum[20] += int16(iter._hash >> 43 & 1)
			iter.sum[21] += int16(iter._hash >> 42 & 1)
			iter.sum[22] += int16(iter._hash >> 41 & 1)
			iter.sum[23] += int16(iter._hash >> 40 & 1)
			iter.sum[24] += int16(iter._hash >> 39 & 1)
			iter.sum[25] += int16(iter._hash >> 38 & 1)
			iter.sum[26] += int16(iter._hash >> 37 & 1)
			iter.sum[27] += int16(iter._hash >> 36 & 1)
			iter.sum[28] += int16(iter._hash >> 35 & 1)
			iter.sum[29] += int16(iter._hash >> 34 & 1)
			iter.sum[30] += int16(iter._hash >> 33 & 1)
			iter.sum[31] += int16(iter._hash >> 32 & 1)
			iter.sum[32] += int16(iter._hash >> 31 & 1)
			iter.sum[33] += int16(iter._hash >> 30 & 1)
			iter.sum[34] += int16(iter._hash >> 29 & 1)
			iter.sum[35] += int16(iter._hash >> 28 & 1)
			iter.sum[36] += int16(iter._hash >> 27 & 1)
			iter.sum[37] += int16(iter._hash >> 26 & 1)
			iter.sum[38] += int16(iter._hash >> 25 & 1)
			iter.sum[39] += int16(iter._hash >> 24 & 1)
			iter.sum[40] += int16(iter._hash >> 23 & 1)
			iter.sum[41] += int16(iter._hash >> 22 & 1)
			iter.sum[42] += int16(iter._hash >> 21 & 1)
			iter.sum[43] += int16(iter._hash >> 20 & 1)
			iter.sum[44] += int16(iter._hash >> 19 & 1)
			iter.sum[45] += int16(iter._hash >> 18 & 1)
			iter.sum[46] += int16(iter._hash >> 17 & 1)
			iter.sum[47] += int16(iter._hash >> 16 & 1)
			iter.sum[48] += int16(iter._hash >> 15 & 1)
			iter.sum[49] += int16(iter._hash >> 14 & 1)
			iter.sum[50] += int16(iter._hash >> 13 & 1)
			iter.sum[51] += int16(iter._hash >> 12 & 1)
			iter.sum[52] += int16(iter._hash >> 11 & 1)
			iter.sum[53] += int16(iter._hash >> 10 & 1)
			iter.sum[54] += int16(iter._hash >> 9 & 1)
			iter.sum[55] += int16(iter._hash >> 8 & 1)
			iter.sum[56] += int16(iter._hash >> 7 & 1)
			iter.sum[57] += int16(iter._hash >> 6 & 1)
			iter.sum[58] += int16(iter._hash >> 5 & 1)
			iter.sum[59] += int16(iter._hash >> 4 & 1)
			iter.sum[60] += int16(iter._hash >> 3 & 1)
			iter.sum[61] += int16(iter._hash >> 2 & 1)
			iter.sum[62] += int16(iter._hash >> 1 & 1)
			iter.sum[63] += int16(iter._hash & 1)
		}

		// decoding vector to hash

		code = 0
		iter.nThsld = (iter.nPosHash + 1) / 2 // the threshold
		// fmt.Printf("nThsld: %d\n", iter.nThsld)

		if iter.nPosHash > 0 {
			// if n >= N/2 ? 1 : 0
			code |= uint64((iter.sum[0]-iter.nThsld)>>15&1^1) << 63
			code |= uint64((iter.sum[1]-iter.nThsld)>>15&1^1) << 62
			code |= uint64((iter.sum[2]-iter.nThsld)>>15&1^1) << 61
			code |= uint64((iter.sum[3]-iter.nThsld)>>15&1^1) << 60
			code |= uint64((iter.sum[4]-iter.nThsld)>>15&1^1) << 59
			code |= uint64((iter.sum[5]-iter.nThsld)>>15&1^1) << 58
			code |= uint64((iter.sum[6]-iter.nThsld)>>15&1^1) << 57
			code |= uint64((iter.sum[7]-iter.nThsld)>>15&1^1) << 56
			code |= uint64((iter.sum[8]-iter.nThsld)>>15&1^1) << 55
			code |= uint64((iter.sum[9]-iter.nThsld)>>15&1^1) << 54
			code |= uint64((iter.sum[10]-iter.nThsld)>>15&1^1) << 53
			code |= uint64((iter.sum[11]-iter.nThsld)>>15&1^1) << 52
			code |= uint64((iter.sum[12]-iter.nThsld)>>15&1^1) << 51
			code |= uint64((iter.sum[13]-iter.nThsld)>>15&1^1) << 50
			code |= uint64((iter.sum[14]-iter.nThsld)>>15&1^1) << 49
			code |= uint64((iter.sum[15]-iter.nThsld)>>15&1^1) << 48
			code |= uint64((iter.sum[16]-iter.nThsld)>>15&1^1) << 47
			code |= uint64((iter.sum[17]-iter.nThsld)>>15&1^1) << 46
			code |= uint64((iter.sum[18]-iter.nThsld)>>15&1^1) << 45
			code |= uint64((iter.sum[19]-iter.nThsld)>>15&1^1) << 44
			code |= uint64((iter.sum[20]-iter.nThsld)>>15&1^1) << 43
			code |= uint64((iter.sum[21]-iter.nThsld)>>15&1^1) << 42
			code |= uint64((iter.sum[22]-iter.nThsld)>>15&1^1) << 41
			code |= uint64((iter.sum[23]-iter.nThsld)>>15&1^1) << 40
			code |= uint64((iter.sum[24]-iter.nThsld)>>15&1^1) << 39
			code |= uint64((iter.sum[25]-iter.nThsld)>>15&1^1) << 38
			code |= uint64((iter.sum[26]-iter.nThsld)>>15&1^1) << 37
			code |= uint64((iter.sum[27]-iter.nThsld)>>15&1^1) << 36
			code |= uint64((iter.sum[28]-iter.nThsld)>>15&1^1) << 35
			code |= uint64((iter.sum[29]-iter.nThsld)>>15&1^1) << 34
			code |= uint64((iter.sum[30]-iter.nThsld)>>15&1^1) << 33
			code |= uint64((iter.sum[31]-iter.nThsld)>>15&1^1) << 32
			code |= uint64((iter.sum[32]-iter.nThsld)>>15&1^1) << 31
			code |= uint64((iter.sum[33]-iter.nThsld)>>15&1^1) << 30
			code |= uint64((iter.sum[34]-iter.nThsld)>>15&1^1) << 29
			code |= uint64((iter.sum[35]-iter.nThsld)>>15&1^1) << 28
			code |= uint64((iter.sum[36]-iter.nThsld)>>15&1^1) << 27
			code |= uint64((iter.sum[37]-iter.nThsld)>>15&1^1) << 26
			code |= uint64((iter.sum[38]-iter.nThsld)>>15&1^1) << 25
			code |= uint64((iter.sum[39]-iter.nThsld)>>15&1^1) << 24
			code |= uint64((iter.sum[40]-iter.nThsld)>>15&1^1) << 23
			code |= uint64((iter.sum[41]-iter.nThsld)>>15&1^1) << 22
			code |= uint64((iter.sum[42]-iter.nThsld)>>15&1^1) << 21
			code |= uint64((iter.sum[43]-iter.nThsld)>>15&1^1) << 20
			code |= uint64((iter.sum[44]-iter.nThsld)>>15&1^1) << 19
			code |= uint64((iter.sum[45]-iter.nThsld)>>15&1^1) << 18
			code |= uint64((iter.sum[46]-iter.nThsld)>>15&1^1) << 17
			code |= uint64((iter.sum[47]-iter.nThsld)>>15&1^1) << 16
			code |= uint64((iter.sum[48]-iter.nThsld)>>15&1^1) << 15
			code |= uint64((iter.sum[49]-iter.nThsld)>>15&1^1) << 14
			code |= uint64((iter.sum[50]-iter.nThsld)>>15&1^1) << 13
			code |= uint64((iter.sum[51]-iter.nThsld)>>15&1^1) << 12
			code |= uint64((iter.sum[52]-iter.nThsld)>>15&1^1) << 11
			code |= uint64((iter.sum[53]-iter.nThsld)>>15&1^1) << 10
			code |= uint64((iter.sum[54]-iter.nThsld)>>15&1^1) << 9
			code |= uint64((iter.sum[55]-iter.nThsld)>>15&1^1) << 8
			code |= uint64((iter.sum[56]-iter.nThsld)>>15&1^1) << 7
			code |= uint64((iter.sum[57]-iter.nThsld)>>15&1^1) << 6
			code |= uint64((iter.sum[58]-iter.nThsld)>>15&1^1) << 5
			code |= uint64((iter.sum[59]-iter.nThsld)>>15&1^1) << 4
			code |= uint64((iter.sum[60]-iter.nThsld)>>15&1^1) << 3
			code |= uint64((iter.sum[61]-iter.nThsld)>>15&1^1) << 2
			code |= uint64((iter.sum[62]-iter.nThsld)>>15&1^1) << 1
			code |= uint64((iter.sum[63]-iter.nThsld)>>15&1 ^ 1)
		} else {
			code = 0
		}

		// update the preHashI
		if iter.preHashI == iter.k-iter.m {
			iter.preHashI = 0
		} else {
			iter.preHashI++
		}
	} else {
		iter.nPosHash = 0
		for iter.idxJ = 0; iter.idxJ <= iter.em; iter.idxJ++ {
			iter._hash, _ = iter.hasher.Next(iter.canonical)

			// fmt.Println(iter.fracMinHash, iter._hash, iter.maxHash, iter._hash > iter.maxHash)
			if iter.fracMinHash && iter._hash > iter.maxHash {
				iter.hashes[iter.idxJ] = 0
				continue
			}
			iter.hashes[iter.idxJ] = iter._hash

			if iter._hash == 0 {
				continue
			}

			iter.nPosHash++

			// fmt.Printf("[%02d] %02d: %064b %d\n", iter.idx, iter.idxJ, iter._hash, iter._hash)

			iter.sum[0] += int16(iter._hash >> 63 & 1)
			iter.sum[1] += int16(iter._hash >> 62 & 1)
			iter.sum[2] += int16(iter._hash >> 61 & 1)
			iter.sum[3] += int16(iter._hash >> 60 & 1)
			iter.sum[4] += int16(iter._hash >> 59 & 1)
			iter.sum[5] += int16(iter._hash >> 58 & 1)
			iter.sum[6] += int16(iter._hash >> 57 & 1)
			iter.sum[7] += int16(iter._hash >> 56 & 1)
			iter.sum[8] += int16(iter._hash >> 55 & 1)
			iter.sum[9] += int16(iter._hash >> 54 & 1)
			iter.sum[10] += int16(iter._hash >> 53 & 1)
			iter.sum[11] += int16(iter._hash >> 52 & 1)
			iter.sum[12] += int16(iter._hash >> 51 & 1)
			iter.sum[13] += int16(iter._hash >> 50 & 1)
			iter.sum[14] += int16(iter._hash >> 49 & 1)
			iter.sum[15] += int16(iter._hash >> 48 & 1)
			iter.sum[16] += int16(iter._hash >> 47 & 1)
			iter.sum[17] += int16(iter._hash >> 46 & 1)
			iter.sum[18] += int16(iter._hash >> 45 & 1)
			iter.sum[19] += int16(iter._hash >> 44 & 1)
			iter.sum[20] += int16(iter._hash >> 43 & 1)
			iter.sum[21] += int16(iter._hash >> 42 & 1)
			iter.sum[22] += int16(iter._hash >> 41 & 1)
			iter.sum[23] += int16(iter._hash >> 40 & 1)
			iter.sum[24] += int16(iter._hash >> 39 & 1)
			iter.sum[25] += int16(iter._hash >> 38 & 1)
			iter.sum[26] += int16(iter._hash >> 37 & 1)
			iter.sum[27] += int16(iter._hash >> 36 & 1)
			iter.sum[28] += int16(iter._hash >> 35 & 1)
			iter.sum[29] += int16(iter._hash >> 34 & 1)
			iter.sum[30] += int16(iter._hash >> 33 & 1)
			iter.sum[31] += int16(iter._hash >> 32 & 1)
			iter.sum[32] += int16(iter._hash >> 31 & 1)
			iter.sum[33] += int16(iter._hash >> 30 & 1)
			iter.sum[34] += int16(iter._hash >> 29 & 1)
			iter.sum[35] += int16(iter._hash >> 28 & 1)
			iter.sum[36] += int16(iter._hash >> 27 & 1)
			iter.sum[37] += int16(iter._hash >> 26 & 1)
			iter.sum[38] += int16(iter._hash >> 25 & 1)
			iter.sum[39] += int16(iter._hash >> 24 & 1)
			iter.sum[40] += int16(iter._hash >> 23 & 1)
			iter.sum[41] += int16(iter._hash >> 22 & 1)
			iter.sum[42] += int16(iter._hash >> 21 & 1)
			iter.sum[43] += int16(iter._hash >> 20 & 1)
			iter.sum[44] += int16(iter._hash >> 19 & 1)
			iter.sum[45] += int16(iter._hash >> 18 & 1)
			iter.sum[46] += int16(iter._hash >> 17 & 1)
			iter.sum[47] += int16(iter._hash >> 16 & 1)
			iter.sum[48] += int16(iter._hash >> 15 & 1)
			iter.sum[49] += int16(iter._hash >> 14 & 1)
			iter.sum[50] += int16(iter._hash >> 13 & 1)
			iter.sum[51] += int16(iter._hash >> 12 & 1)
			iter.sum[52] += int16(iter._hash >> 11 & 1)
			iter.sum[53] += int16(iter._hash >> 10 & 1)
			iter.sum[54] += int16(iter._hash >> 9 & 1)
			iter.sum[55] += int16(iter._hash >> 8 & 1)
			iter.sum[56] += int16(iter._hash >> 7 & 1)
			iter.sum[57] += int16(iter._hash >> 6 & 1)
			iter.sum[58] += int16(iter._hash >> 5 & 1)
			iter.sum[59] += int16(iter._hash >> 4 & 1)
			iter.sum[60] += int16(iter._hash >> 3 & 1)
			iter.sum[61] += int16(iter._hash >> 2 & 1)
			iter.sum[62] += int16(iter._hash >> 1 & 1)
			iter.sum[63] += int16(iter._hash & 1)
		}

		// fmt.Printf("   sum: %v\n", iter.sum)

		code = 0
		iter.nThsld = (iter.nPosHash + 1) / 2
		// fmt.Printf("nThsld: %d\n", iter.nThsld)
		if iter.nPosHash > 0 {
			// if n >= N/2 ? 1 : 0
			code |= uint64((iter.sum[0]-iter.nThsld)>>15&1^1) << 63
			code |= uint64((iter.sum[1]-iter.nThsld)>>15&1^1) << 62
			code |= uint64((iter.sum[2]-iter.nThsld)>>15&1^1) << 61
			code |= uint64((iter.sum[3]-iter.nThsld)>>15&1^1) << 60
			code |= uint64((iter.sum[4]-iter.nThsld)>>15&1^1) << 59
			code |= uint64((iter.sum[5]-iter.nThsld)>>15&1^1) << 58
			code |= uint64((iter.sum[6]-iter.nThsld)>>15&1^1) << 57
			code |= uint64((iter.sum[7]-iter.nThsld)>>15&1^1) << 56
			code |= uint64((iter.sum[8]-iter.nThsld)>>15&1^1) << 55
			code |= uint64((iter.sum[9]-iter.nThsld)>>15&1^1) << 54
			code |= uint64((iter.sum[10]-iter.nThsld)>>15&1^1) << 53
			code |= uint64((iter.sum[11]-iter.nThsld)>>15&1^1) << 52
			code |= uint64((iter.sum[12]-iter.nThsld)>>15&1^1) << 51
			code |= uint64((iter.sum[13]-iter.nThsld)>>15&1^1) << 50
			code |= uint64((iter.sum[14]-iter.nThsld)>>15&1^1) << 49
			code |= uint64((iter.sum[15]-iter.nThsld)>>15&1^1) << 48
			code |= uint64((iter.sum[16]-iter.nThsld)>>15&1^1) << 47
			code |= uint64((iter.sum[17]-iter.nThsld)>>15&1^1) << 46
			code |= uint64((iter.sum[18]-iter.nThsld)>>15&1^1) << 45
			code |= uint64((iter.sum[19]-iter.nThsld)>>15&1^1) << 44
			code |= uint64((iter.sum[20]-iter.nThsld)>>15&1^1) << 43
			code |= uint64((iter.sum[21]-iter.nThsld)>>15&1^1) << 42
			code |= uint64((iter.sum[22]-iter.nThsld)>>15&1^1) << 41
			code |= uint64((iter.sum[23]-iter.nThsld)>>15&1^1) << 40
			code |= uint64((iter.sum[24]-iter.nThsld)>>15&1^1) << 39
			code |= uint64((iter.sum[25]-iter.nThsld)>>15&1^1) << 38
			code |= uint64((iter.sum[26]-iter.nThsld)>>15&1^1) << 37
			code |= uint64((iter.sum[27]-iter.nThsld)>>15&1^1) << 36
			code |= uint64((iter.sum[28]-iter.nThsld)>>15&1^1) << 35
			code |= uint64((iter.sum[29]-iter.nThsld)>>15&1^1) << 34
			code |= uint64((iter.sum[30]-iter.nThsld)>>15&1^1) << 33
			code |= uint64((iter.sum[31]-iter.nThsld)>>15&1^1) << 32
			code |= uint64((iter.sum[32]-iter.nThsld)>>15&1^1) << 31
			code |= uint64((iter.sum[33]-iter.nThsld)>>15&1^1) << 30
			code |= uint64((iter.sum[34]-iter.nThsld)>>15&1^1) << 29
			code |= uint64((iter.sum[35]-iter.nThsld)>>15&1^1) << 28
			code |= uint64((iter.sum[36]-iter.nThsld)>>15&1^1) << 27
			code |= uint64((iter.sum[37]-iter.nThsld)>>15&1^1) << 26
			code |= uint64((iter.sum[38]-iter.nThsld)>>15&1^1) << 25
			code |= uint64((iter.sum[39]-iter.nThsld)>>15&1^1) << 24
			code |= uint64((iter.sum[40]-iter.nThsld)>>15&1^1) << 23
			code |= uint64((iter.sum[41]-iter.nThsld)>>15&1^1) << 22
			code |= uint64((iter.sum[42]-iter.nThsld)>>15&1^1) << 21
			code |= uint64((iter.sum[43]-iter.nThsld)>>15&1^1) << 20
			code |= uint64((iter.sum[44]-iter.nThsld)>>15&1^1) << 19
			code |= uint64((iter.sum[45]-iter.nThsld)>>15&1^1) << 18
			code |= uint64((iter.sum[46]-iter.nThsld)>>15&1^1) << 17
			code |= uint64((iter.sum[47]-iter.nThsld)>>15&1^1) << 16
			code |= uint64((iter.sum[48]-iter.nThsld)>>15&1^1) << 15
			code |= uint64((iter.sum[49]-iter.nThsld)>>15&1^1) << 14
			code |= uint64((iter.sum[50]-iter.nThsld)>>15&1^1) << 13
			code |= uint64((iter.sum[51]-iter.nThsld)>>15&1^1) << 12
			code |= uint64((iter.sum[52]-iter.nThsld)>>15&1^1) << 11
			code |= uint64((iter.sum[53]-iter.nThsld)>>15&1^1) << 10
			code |= uint64((iter.sum[54]-iter.nThsld)>>15&1^1) << 9
			code |= uint64((iter.sum[55]-iter.nThsld)>>15&1^1) << 8
			code |= uint64((iter.sum[56]-iter.nThsld)>>15&1^1) << 7
			code |= uint64((iter.sum[57]-iter.nThsld)>>15&1^1) << 6
			code |= uint64((iter.sum[58]-iter.nThsld)>>15&1^1) << 5
			code |= uint64((iter.sum[59]-iter.nThsld)>>15&1^1) << 4
			code |= uint64((iter.sum[60]-iter.nThsld)>>15&1^1) << 3
			code |= uint64((iter.sum[61]-iter.nThsld)>>15&1^1) << 2
			code |= uint64((iter.sum[62]-iter.nThsld)>>15&1^1) << 1
			code |= uint64((iter.sum[63]-iter.nThsld)>>15&1 ^ 1)
		} else {
			code = 0
		}

		iter.preHashI = 0
		iter.first = false
	}

	// code &= 0xffffffff

	// fmt.Printf("kmer: %03d-%s, %064b, %d\n\n", iter.idx, iter.s.Seq[iter.idx:iter.idx+iter.k], code, code)

	iter.preCode = code
	iter.idx++

	return code, true
}

// NewHashIterator returns ntHash Iterator.
func NewHashIterator(s *seq.Seq, k int, canonical bool, circular bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s.Seq) < k {
		return nil, ErrShortSeq
	}

	// iter := &Iterator{s: s, k: k, canonical: canonical, circular: circular}
	iter := poolIterator.Get().(*Iterator)
	iter.s = s
	iter.k = k
	iter.canonical = canonical
	iter.circular = circular
	iter.finished = false
	iter.revcomStrand = false
	iter.idx = 0

	iter.hash = true
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)
	// iter.mask1 = (1 << (iter.kP1Uint << 1)) - 1
	// iter.mask2 = iter.kP1Uint << 1

	var err error
	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(s.Seq), len(s.Seq)+k-1)
		copy(seq2, s.Seq) // do not edit original sequence
		seq2 = append(seq2, s.Seq[0:k-1]...)
	} else {
		seq2 = s.Seq
	}
	iter.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	return iter, nil
}

// NextHash returns next ntHash.
func (iter *Iterator) NextHash() (code uint64, ok bool) {
	code, ok = iter.hasher.Next(iter.canonical)
	if !ok {
		poolIterator.Put(iter)
	}
	iter.idx++
	return code, ok
}

// NewKmerIterator returns k-mer code iterator.
func NewKmerIterator(s *seq.Seq, k int, canonical bool, circular bool) (*Iterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s.Seq) < k {
		return nil, ErrShortSeq
	}

	var s2 *seq.Seq
	if circular {
		s2 = s.Clone() // do not edit original sequence
		s2.Seq = append(s2.Seq, s.Seq[0:k-1]...)
	} else {
		s2 = s
	}

	// iter := &Iterator{s: s2, k: k, canonical: canonical, circular: circular}
	iter := poolIterator.Get().(*Iterator)
	iter.s = s2
	iter.k = k
	iter.canonical = canonical
	iter.circular = circular
	iter.finished = false
	iter.revcomStrand = false
	iter.idx = 0

	iter.length = len(s2.Seq)
	iter.end = iter.length - k + 1
	iter.kUint = uint(k)
	iter.kP1 = k - 1
	iter.kP1Uint = uint(k - 1)
	iter.mask1 = (1 << (iter.kP1Uint << 1)) - 1
	iter.mask2 = iter.kP1Uint << 1

	iter.first = true

	return iter, nil
}

// NextKmer returns next k-mer code.
func (iter *Iterator) NextKmer() (code uint64, ok bool, err error) {
	if iter.finished {
		return 0, false, nil
	}

	if iter.idx == iter.end {
		if iter.canonical || iter.revcomStrand {
			iter.finished = true
			poolIterator.Put(iter)
			return 0, false, nil
		}
		iter.s.RevComInplace()
		iter.idx = 0
		iter.revcomStrand = true
		iter.first = true
	}

	iter.e = iter.idx + iter.k
	iter.kmer = iter.s.Seq[iter.idx:iter.e]

	if !iter.first {
		iter.codeBase = base2bit[iter.kmer[iter.kP1]]
		if iter.codeBase == 4 {
			err = ErrIllegalBase
		}

		// compute code from previous one
		//  code = iter.preCode&((1<<(iter.kP1Uint<<1))-1)<<2 | iter.codeBase
		code = (iter.preCode&iter.mask1)<<2 | iter.codeBase

		// compute code of revcomp kmer from previous one
		// iter.codeRC = (iter.codeBase^3)<<(iter.kP1Uint<<1) | (iter.preCodeRC >> 2)
		iter.codeRC = (iter.codeBase^3)<<(iter.mask2) | (iter.preCodeRC >> 2)
	} else {
		code, err = kmers.Encode(iter.kmer)
		iter.codeRC = kmers.MustRevComp(code, iter.k)
		iter.first = false
	}
	if err != nil {
		return 0, false, fmt.Errorf("encode %s: %s", iter.kmer, err)
	}

	iter.preCode = code
	iter.preCodeRC = iter.codeRC
	iter.idx++

	if iter.canonical && code > iter.codeRC {
		code = iter.codeRC
	}

	return code, true, nil
}

// Next is a wrapter for NextHash and NextKmer.
func (iter *Iterator) Next() (code uint64, ok bool, err error) {
	if iter.hash {
		if iter.simhash {
			code, ok = iter.NextSimHash()
			return
		}
		code, ok = iter.NextHash()
		return
	}
	code, ok, err = iter.NextKmer()
	return
}

// Index returns current 0-baesd index.
func (iter *Iterator) Index() int {
	return iter.idx - 1
}
