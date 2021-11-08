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
	"sync"

	"github.com/shenwei356/bio/seq"
	"github.com/zeebo/wyhash"
)

var poolProteinIterator = &sync.Pool{New: func() interface{} {
	return &ProteinIterator{}
}}

// ProteinIterator is a iterator for protein sequence.
type ProteinIterator struct {
	s0 *seq.Seq // only used for KmerProteinIterator
	s  *seq.Seq // amino acid

	k        int
	finished bool
	end      int
	idx      int
}

// NewProteinIterator returns an iterator for hash of amino acids
func NewProteinIterator(s *seq.Seq, k int, codonTable int, frame int) (*ProteinIterator, error) {
	if k < 1 {
		return nil, ErrInvalidK
	}
	if len(s.Seq) < k*3 {
		return nil, ErrShortSeq
	}

	// iter := &ProteinIterator{s0: s, k: k}
	iter := poolProteinIterator.Get().(*ProteinIterator)
	iter.s0 = s
	iter.k = k
	iter.finished = false
	iter.idx = 0

	var err error
	if s.Alphabet != seq.Protein {
		iter.s, err = s.Translate(codonTable, frame, false, false, true, false)
		if err != nil {
			return nil, err
		}
	} else {
		iter.s = s
	}
	iter.end = len(iter.s.Seq) - k

	return iter, nil
}

// Next return's a hash
func (iter *ProteinIterator) Next() (code uint64, ok bool) {
	if iter.finished {
		return 0, false
	}

	if iter.idx > iter.end {
		iter.finished = true
		poolProteinIterator.Put(iter)
		return 0, false
	}

	code = wyhash.Hash(iter.s.Seq[iter.idx:iter.idx+iter.k], 1)
	iter.idx++
	return code, true
}

// Index returns current 0-baesd index.
func (iter *ProteinIterator) Index() int {
	return iter.idx - 1
}
