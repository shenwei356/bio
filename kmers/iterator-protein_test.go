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
	"testing"

	"github.com/shenwei356/bio/seq"
)

func TestProteinIterator(t *testing.T) {
	_s := "AAGTTTGAATCATTCAACTATCTAGTTTTCAGAGAACAATGTTCTCTAAAGAATAGAAAAGAGTCATTGTGCGGTGATGATGGCGGGAAGGATCCACCTG"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 10

	iter, err := NewProteinIterator(sequence, k, 1, 1)
	if err != nil {
		t.Errorf("fail to create aa iter rator")
	}

	var code uint64
	var ok bool
	// var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok = iter.Next()
		if !ok {
			break
		}

		// idx = iter.Index()
		// fmt.Printf("aa: %d-%s, %d\n", idx, iter.s.Seq[idx:idx+k], code)

		codes = append(codes, code)
	}

	if len(codes) != len(_s)/3-k+1 {
		t.Errorf("k-mer hashes number error")
	}

}
