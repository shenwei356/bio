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
	"sort"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/zeebo/wyhash"
)

// ProteinMinimizerSketch is a protein k-mer minimizer iterator
type ProteinMinimizerSketch struct {
	s *seq.Seq // amino acid

	k    int
	end0 int
	idx  int

	// ----------------------

	skip bool
	w    int

	end int
	r   int // L-s

	i, mI     int
	mV        uint64
	preMinIdx int

	buf     []IdxValue
	i2v     IdxValue
	flag    bool
	t, b, e int
}

var poolProteinMinimizerSketch = &sync.Pool{New: func() interface{} {
	return &ProteinMinimizerSketch{}
}}

// NewProteinMinimizerSketch returns a ProteinMinimizerSketch
func NewProteinMinimizerSketch(S *seq.Seq, k int, codonTable int, frame int, w int) (*ProteinMinimizerSketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(S.Seq) < k*3 {
		return nil, ErrShortSeq
	}

	if w < 1 || w > (1<<31)-1 {
		return nil, ErrInvalidW
	}
	if len(S.Seq) < k*3+w-1 {
		return nil, ErrShortSeq
	}

	// s := &ProteinMinimizerSketch{s0: S, k: k, w: w}
	s := poolProteinMinimizerSketch.Get().(*ProteinMinimizerSketch)
	s.k = k
	s.w = w

	var err error
	if S.Alphabet != seq.Protein {
		s.s, err = S.Translate(codonTable, frame, false, false, true, false)
		if err != nil {
			return nil, err
		}
	} else {
		s.s = S
	}

	s.idx = 0
	s.end0 = len(s.s.Seq) - k

	s.skip = w == 1
	s.end = len(s.s.Seq) - 1
	s.r = w - 1 // L-k

	s.buf = make([]IdxValue, 0, w)
	s.preMinIdx = -1

	return s, nil
}

// Next returns next hash value
func (s *ProteinMinimizerSketch) Next() (code uint64, ok bool) {
	for {
		// if s.idx > s.end {
		// 	return 0, false
		// }

		if s.idx > s.end0 {
			poolProteinIterator.Put(s)
			return 0, false
		}

		code = wyhash.Hash(s.s.Seq[s.idx:s.idx+s.k], 1)

		if s.skip {
			s.mI = s.idx
			s.idx++
			return code, true
		}

		// in window
		if s.idx < s.r {
			s.buf = append(s.buf, IdxValue{Idx: s.idx, Val: code})

			s.idx++
			continue
		}

		// end of w
		if s.idx == s.r {
			s.buf = append(s.buf, IdxValue{Idx: s.idx, Val: code})
			sort.Sort(idxValues(s.buf)) // sort

			s.i2v = s.buf[0]

			s.mI, s.mV = s.i2v.Idx, s.i2v.Val
			s.preMinIdx = s.mI

			s.idx++
			return s.i2v.Val, true
		}

		// find min k-mer
		// remove k-mer not in this window.
		// have to check position/index one by one
		for s.i, s.i2v = range s.buf {
			if s.i2v.Idx == s.idx-s.w {
				if s.i < s.r {
					copy(s.buf[s.i:s.r], s.buf[s.i+1:])
				} // happen to be at the end
				s.buf = s.buf[:s.r]
				break
			}
		}

		// add new k-mer
		s.flag = false
		// using binary search, faster han linear search
		s.b, s.e = 0, s.r-1
		for {
			s.t = s.b + (s.e-s.b)/2
			if code < s.buf[s.t].Val {
				s.e = s.t - 1 // end search here
				if s.e <= s.b {
					s.flag = true
					s.i = s.b
					break
				}
			} else {
				s.b = s.t + 1 // start here
				if s.b >= s.r {
					s.flag = false
					break
				}
				if s.b >= s.e {
					s.flag = true
					s.i = s.e // right here
					break
				}
			}
		}
		if !s.flag { // it's the biggest one, append to the end
			s.buf = append(s.buf, IdxValue{s.idx, code})
		} else {
			if code >= s.buf[s.i].Val { // have to check again
				s.i++
			}
			s.buf = append(s.buf, blankI2V)     // append one element
			copy(s.buf[s.i+1:], s.buf[s.i:s.r]) // move right
			s.buf[s.i] = IdxValue{s.idx, code}
		}

		s.i2v = s.buf[0]
		if s.i2v.Idx == s.preMinIdx { // deduplicate
			s.idx++
			continue
		}

		s.mI, s.mV = s.i2v.Idx, s.i2v.Val
		s.preMinIdx = s.mI

		s.idx++
		return s.i2v.Val, true
	}
}

// Index returns current 0-baesd index.
func (s *ProteinMinimizerSketch) Index() int {
	return s.mI
}
