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
	"slices"
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestProteinMinimizer(t *testing.T) {
	_s := "AAGTTTGAATCATTCAACTATCTAGTTTTCAGAGAACAATGTTCTCTAAAGAATAGAAAAGAGTCATTGTGCGGTGATGATGGCGGGAAGGATCCACCTG"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 10
	w := 3

	sketch, err := NewProteinMinimizerSketch(sequence, k, 1, 1, w)
	if err != nil {
		t.Errorf("fail to create minizimer sketch")
	}

	var code uint64
	var ok bool
	// var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok = sketch.Next()
		if !ok {
			break
		}

		// idx = sketch.Index()
		// fmt.Printf("aa: %d-%s, %d\n", idx, sketch.s.Seq[idx:idx+k], code)

		codes = append(codes, code)
	}
	if !sketch.finished {
		t.Fatal("exhausted protein minimizer sketch was not marked finished")
	}
	for i := 0; i < 3; i++ {
		if code, ok = sketch.Next(); ok || code != 0 {
			t.Fatalf("finished protein minimizer sketch returned (%d, %t)", code, ok)
		}
	}

	reused, err := NewProteinMinimizerSketch(sequence, k, 1, 1, w)
	if err != nil {
		t.Fatal(err)
	}
	var reusedCodes []uint64
	for {
		code, ok = reused.Next()
		if !ok {
			break
		}
		reusedCodes = append(reusedCodes, code)
	}
	if !slices.Equal(codes, reusedCodes) {
		t.Fatal("protein minimizer output changed after reusing a pooled sketch")
	}

	// A protein minimizer must never contaminate the ProteinIterator pool.
	iter, err := NewProteinIterator(sequence, k, 1, 1)
	if err != nil {
		t.Fatal(err)
	}
	for {
		_, ok = iter.Next()
		if !ok {
			break
		}
	}
}
