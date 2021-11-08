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

package kmers

import (
	"fmt"
	"sort"
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/will-rowe/nthash"
	// hasher "github.com/zeebo/wyhash"
)

// ErrInvalidS means s >= k.
var ErrInvalidS = fmt.Errorf("kmers: invalid s-mer size")

// ErrInvalidW means w < 2 or w > (1<<32)-1
var ErrInvalidW = fmt.Errorf("kmers: invalid minimimzer window")

// ErrBufNil means the buffer is nil
var ErrBufNil = fmt.Errorf("kmers: buffer slice is nil")

// ErrBufNotEmpty means the buffer has some elements
var ErrBufNotEmpty = fmt.Errorf("kmers: buffer has elements")

// Sketch is a k-mer sketch iterator
type Sketch struct {
	S        []byte
	k        int
	s        int
	circular bool
	hasher   *nthash.NTHi

	kMs int // k-s, just for syncmer
	r   int // L-s

	idx int // current location, 0-based
	end int

	i, mI     int
	v, mV     uint64
	preMinIdx int

	buf     []IdxValue
	i2v     IdxValue
	flag    bool
	t, b, e int

	// ------ just for syncmer -------
	hasherS           *nthash.NTHi
	bsyncmerIdx       int
	lateOutputThisOne bool
	preMinIdxs        []int

	// ------ just for minimizer -----
	skip      bool
	minimizer bool
	w         int
}

var poolSketch = &sync.Pool{New: func() interface{} {
	return &Sketch{}
}}

// NewMinimizerSketch returns a SyncmerSketch Iterator.
// It returns the minHashes in all windows of w (w>=1) bp.
func NewMinimizerSketch(S *seq.Seq, k int, w int, circular bool) (*Sketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if w < 1 || w > (1<<31)-1 {
		return nil, ErrInvalidW
	}
	if len(S.Seq) < k+w-1 {
		return nil, ErrShortSeq
	}

	// sketch := &Sketch{S: S.Seq, w: w, k: k, circular: circular}
	sketch := poolSketch.Get().(*Sketch)

	sketch.minimizer = true
	sketch.k = k
	sketch.w = w
	sketch.circular = circular
	sketch.skip = w == 1

	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(S.Seq), len(S.Seq)+k-1)
		copy(seq2, S.Seq) // do not edit original sequence
		seq2 = append(seq2, S.Seq[0:k-1]...)
		sketch.S = seq2
	} else {
		seq2 = S.Seq
	}

	sketch.idx = 0
	sketch.end = len(seq2) - 1
	sketch.r = w - 1 // L-k

	var err error
	sketch.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	if sketch.buf == nil {
		sketch.buf = make([]IdxValue, 0, w)
	} else {
		sketch.buf = sketch.buf[:0]
	}
	if sketch.preMinIdxs == nil {
		sketch.preMinIdxs = make([]int, 0, 8)
	} else {
		sketch.preMinIdxs = sketch.preMinIdxs[:0]
	}
	sketch.preMinIdx = -1

	return sketch, nil
}

// NewSyncmerSketch returns a SyncmerSketch Iterator.
// 1<=s<=k.
func NewSyncmerSketch(S *seq.Seq, k int, s int, circular bool) (*Sketch, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if s > k || s == 0 {
		return nil, ErrInvalidS
	}
	if len(S.Seq) < k*2-s-1 {
		return nil, ErrShortSeq
	}

	// sketch := &Sketch{S: S.Seq, s: s, k: k, circular: circular}
	sketch := poolSketch.Get().(*Sketch)

	sketch.minimizer = false
	sketch.k = k
	sketch.s = s
	sketch.circular = circular
	sketch.skip = s == k

	var seq2 []byte
	if circular {
		seq2 = make([]byte, len(S.Seq), len(S.Seq)+k-1)
		copy(seq2, S.Seq) // do not edit original sequence
		seq2 = append(seq2, S.Seq[0:k-1]...)
		sketch.S = seq2
	} else {
		seq2 = S.Seq
	}

	sketch.idx = 0
	sketch.end = len(seq2) - 2*k + s + 1 // len(sequence) - L (2*k - s - 1)
	sketch.r = 2*k - s - 1 - s           // L-s
	sketch.kMs = k - s                   // k-s
	sketch.w = k - s

	var err error
	sketch.hasher, err = nthash.NewHasher(&seq2, uint(k))
	if err != nil {
		return nil, err
	}

	sketch.hasherS, err = nthash.NewHasher(&seq2, uint(s))
	if err != nil {
		return nil, err
	}

	if sketch.buf == nil {
		sketch.buf = make([]IdxValue, 0, (k-s)<<1)
	} else {
		sketch.buf = sketch.buf[:0]
	}
	if sketch.preMinIdxs == nil {
		sketch.preMinIdxs = make([]int, 0, 8)
	} else {
		sketch.preMinIdxs = sketch.preMinIdxs[:0]
	}
	sketch.preMinIdx = -1

	return sketch, nil
}

// NextMinimizer returns next minimizer.
func (s *Sketch) NextMinimizer() (code uint64, ok bool) {
	for {
		if s.idx > s.end {
			return 0, false
		}

		// nthash of current k-mer
		code, ok = s.hasher.Next(true)
		if !ok {
			poolSketch.Put(s)
			return code, false
		}

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

// NextSyncmer returns next syncmer.
func (s *Sketch) NextSyncmer() (code uint64, ok bool) {
	for {
		if s.idx > s.end {
			return 0, false
		}

		// nthash of current k-mer
		code, ok = s.hasher.Next(true)
		if !ok {
			poolSketch.Put(s)
			return code, false
		}

		// fmt.Printf("\nidx: %d, %s, %d\n", s.idx, s.S[s.idx:s.idx+s.s], code)
		// fmt.Printf("idx: %d, pres: %v, pre: %d\n", s.idx, s.preMinIdxs, s.preMinIdx)

		if s.skip {
			s.idx++
			return code, true
		}

		if len(s.preMinIdxs) > 0 && s.idx == s.preMinIdxs[0] {
			// we will output this one in this round
			s.lateOutputThisOne = true
		} else {
			s.lateOutputThisOne = false
		}

		// find min s-mer
		if s.idx == 0 {
			for s.i = s.idx; s.i <= s.idx+s.r; s.i++ {
				// fmt.Printf("s: %d\n", s.i)
				s.v, ok = s.hasherS.Next(true)
				if !ok {
					return code, false
				}
				s.buf = append(s.buf, IdxValue{Idx: s.i, Val: s.v})
			}
			sort.Sort(idxValues(s.buf))
		} else {
			// remove s-mer not in this window.
			// have to check position/index one by one
			for s.i, s.i2v = range s.buf {
				if s.i2v.Idx == s.idx-1 {
					if s.i < s.r {
						copy(s.buf[s.i:s.r], s.buf[s.i+1:])
					} // happen to be at the end
					s.buf = s.buf[:s.r]
					break
				}
			}

			// add new s-mer
			// fmt.Printf("s: %d\n", s.idx+s.r)
			s.v, ok = s.hasherS.Next(true)
			if !ok {
				return code, false
			}
			s.flag = false
			// using binary search, faster han linear search
			s.b, s.e = 0, s.r-1
			for {
				s.t = s.b + (s.e-s.b)/2
				if s.v < s.buf[s.t].Val {
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
				s.buf = append(s.buf, IdxValue{s.idx + s.r, s.v})
			} else {
				if s.v >= s.buf[s.i].Val { // have to check again
					s.i++
				}
				s.buf = append(s.buf, blankI2V)     // append one element
				copy(s.buf[s.i+1:], s.buf[s.i:s.r]) // move right
				s.buf[s.i] = IdxValue{s.idx + s.r, s.v}
			}
		}

		s.i2v = s.buf[0]
		s.mI, s.mV = s.i2v.Idx, s.i2v.Val

		// fmt.Printf("  smer: %d: %d\n", s.mI, s.mV)

		// find the location of bounded syncmer
		if s.mI-s.idx < s.w { // syncmer at the beginning of kmer
			s.bsyncmerIdx = s.mI
			// fmt.Printf("  bIdx: start: %d\n", s.bsyncmerIdx)
		} else { // at the end
			s.bsyncmerIdx = s.mI - s.kMs
			// fmt.Printf("  bIdx:   end: %d\n", s.bsyncmerIdx)
		}

		// ----------------------------------

		// duplicated
		if len(s.preMinIdxs) > 0 && s.bsyncmerIdx == s.preMinIdxs[0] {
			// fmt.Printf("  duplicated:  %d\n", s.bsyncmerIdx)
			if s.lateOutputThisOne {
				// remove the first element
				copy(s.preMinIdxs[0:len(s.preMinIdxs)-1], s.preMinIdxs[1:])
				s.preMinIdxs = s.preMinIdxs[0 : len(s.preMinIdxs)-1]

				s.idx++
				s.preMinIdx = s.bsyncmerIdx
				return code, true
			}

			s.idx++
			// s.preMinIdx = s.bsyncmerIdx
			continue
		}

		if s.lateOutputThisOne {
			// remove the first element
			copy(s.preMinIdxs[0:len(s.preMinIdxs)-1], s.preMinIdxs[1:])
			s.preMinIdxs = s.preMinIdxs[0 : len(s.preMinIdxs)-1]

			if s.preMinIdx != s.bsyncmerIdx {
				s.preMinIdxs = append(s.preMinIdxs, s.bsyncmerIdx)
			}
			// fmt.Printf("    late2: %d\n", s.preMinIdxs[0])

			s.idx++
			s.preMinIdx = s.bsyncmerIdx
			return code, true
		}

		// is it current kmer?
		if s.bsyncmerIdx == s.idx {
			// fmt.Printf("  current: %d\n", s.bsyncmerIdx)
			if len(s.preMinIdxs) > 0 {
				// remove the first element
				copy(s.preMinIdxs[0:len(s.preMinIdxs)-1], s.preMinIdxs[1:])
				s.preMinIdxs = s.preMinIdxs[0 : len(s.preMinIdxs)-1]
			}
			s.idx++
			s.preMinIdx = s.bsyncmerIdx
			return code, true
		}

		if s.preMinIdx != s.bsyncmerIdx {
			s.preMinIdxs = append(s.preMinIdxs, s.bsyncmerIdx)
		}
		// fmt.Printf("  return it later: %d\n", s.bsyncmerIdx)
		s.idx++
		s.preMinIdx = s.bsyncmerIdx
	}
}

// Next returns next sketch
func (s *Sketch) Next() (uint64, bool) {
	if s.minimizer {
		return s.NextMinimizer()
	}
	return s.NextSyncmer()
}

// Index returns current  0-baesd index
func (s *Sketch) Index() int {
	if s.minimizer {
		return s.mI
	}
	return s.idx - 1
}

// IdxValue is for storing k-mer hash and it's location when computing k-mer sketches.
type IdxValue struct {
	Idx int    // index
	Val uint64 // hash
}

var blankI2V = IdxValue{0, 0}

type idxValues []IdxValue

func (l idxValues) Len() int               { return len(l) }
func (l idxValues) Less(i int, j int) bool { return l[i].Val < l[j].Val }
func (l idxValues) Swap(i int, j int)      { l[i], l[j] = l[j], l[i] }
