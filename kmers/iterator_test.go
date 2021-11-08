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
	"math/rand"
	"testing"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/util/bytesize"
)

func TestKmerIterator(t *testing.T) {
	_s := "AAGTTTGAATCATTCAACTATCTAGTTTTCAGAGAACAATGTTCTCTAAAGAATAGAAAAGAGTCATTGTGCGGTGATGATGGCGGGAAGGATCCACCTG"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 10

	iter, err := NewKmerIterator(sequence, k, true, false)
	if err != nil {
		t.Errorf("fail to create aa iter rator")
	}

	var code uint64
	var ok bool
	// var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok, err = iter.Next()
		if err != nil {
			t.Error(err)
		}
		if !ok {
			break
		}

		// idx = iter.Index()
		// fmt.Printf("kmer: %d-%s, %d\n", idx, iter.s.Seq[idx:idx+k], code)

		codes = append(codes, code)
	}

	if len(codes) != len(_s)-k+1 {
		t.Errorf("k-mers number error")
	}
}

func TestHashIterator(t *testing.T) {
	_s := "AAGTTTGAATCATTCAACTATCTAGTTTTCAGAGAACAATGTTCTCTAAAGAATAGAAAAGAGTCATTGTGCGGTGATGATGGCGGGAAGGATCCACCTG"
	sequence, err := seq.NewSeq(seq.DNA, []byte(_s))
	if err != nil {
		t.Errorf("fail to create sequence: %s", _s)
	}
	k := 10

	iter, err := NewHashIterator(sequence, k, true, false)
	if err != nil {
		t.Errorf("fail to create aa iter rator")
	}

	var code uint64
	var ok bool
	// var idx int
	codes := make([]uint64, 0, 1024)
	for {
		code, ok, err = iter.Next()
		if err != nil {
			t.Error(err)
		}
		if !ok {
			break
		}

		// idx = iter.Index()
		// fmt.Printf("kmer: %d-%s, %d\n", idx, iter.s.Seq[idx:idx+k], code)

		codes = append(codes, code)
	}

	if len(codes) != len(_s)-k+1 {
		t.Errorf("k-mer hashes number error")
	}
}

var benchSeqs []*seq.Seq
var _code uint64

func init() {
	rand.Seed(11)

	sizes := []int{1 << 10} //, 1 << 20, 10 << 20}
	benchSeqs = make([]*seq.Seq, len(sizes))
	var err error
	for i, size := range sizes {
		sequence := make([]byte, size)

		// fmt.Printf("generating pseudo DNA with length of %s ...\n", bytesize.ByteSize(size))
		for j := 0; j < size; j++ {
			sequence[j] = bit2base[rand.Intn(4)]
		}
		benchSeqs[i], err = seq.NewSeq(seq.DNA, sequence)
		if err != nil {
			panic("should not happen")
		}
		// fmt.Println(benchSeqs[i])
	}
	// fmt.Printf("%d DNA sequences generated\n", len(sizes))
}

func BenchmarkKmerIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool

			for j := 0; j < b.N; j++ {
				iter, err := NewKmerIterator(benchSeqs[i], 31, true, false)
				if err != nil {
					b.Errorf("fail to create hash iterator. seq length: %d", size)
				}
				for {
					code, ok, err = iter.NextKmer()
					if err != nil {
						b.Errorf("fail to get kmer code: %d-%s", iter.Index(),
							benchSeqs[i].Seq[iter.Index():iter.Index()+31])
					}

					if !ok {
						break
					}

					_code = code
				}
			}
		})
	}
}

func BenchmarkHashIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool

			for j := 0; j < b.N; j++ {
				iter, err := NewHashIterator(benchSeqs[i], 31, true, false)
				if err != nil {
					b.Errorf("fail to create hash iterator. seq length: %d", size)
				}

				for {
					code, ok = iter.NextHash()
					if !ok {
						break
					}

					_code = code
				}
			}
		})
	}
}

func BenchmarkProteinIterator(b *testing.B) {
	for i := range benchSeqs {
		size := len(benchSeqs[i].Seq)
		b.Run(bytesize.ByteSize(size).String(), func(b *testing.B) {
			var code uint64
			var ok bool

			for j := 0; j < b.N; j++ {
				iter, err := NewProteinIterator(benchSeqs[i], 10, 1, 1)
				if err != nil {
					b.Errorf("fail to create hash iterator. seq length: %d", size)
				}

				for {
					code, ok = iter.Next()
					if !ok {
						break
					}

					_code = code
				}
			}
		})
	}
}
