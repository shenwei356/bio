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
	"bytes"
	"fmt"
	"math/rand"
	"testing"
)

var randomMers [][]byte
var randomMersN = 100000

var benchMer = []byte("ACTGactgGTCAgtcaactgGTCAACTGGTCA")
var codeBenchMer uint64 = 2170370756141391540
var benchMer2 = []byte("CTGactgGTCAgtcaactgGTCAACTGGTCAC")
var codeBenchMer2 uint64 = 8681483024565566161
var benchCode uint64
var benchKmerCode KmerCode

func init() {
	randomMers = make([][]byte, randomMersN)
	for i := 0; i < randomMersN; i++ {
		randomMers[i] = make([]byte, rand.Intn(32)+1)
		for j := range randomMers[i] {
			randomMers[i][j] = bit2base[rand.Intn(4)]
		}
	}

	// for benchmark
	var err error
	benchCode, err = Encode(benchMer)
	if err != nil {
		panic(fmt.Sprintf("init: fail to encode %s", benchMer))
	}

	benchKmerCode, err = NewKmerCode(benchMer)
	if err != nil {
		panic(fmt.Sprintf("init: fail to create KmerCode from %s", benchMer))
	}
}

// TestEncodeDecode tests encode and decode
func TestEncodeDecode(t *testing.T) {
	var kcode KmerCode
	var err error
	for _, mer := range randomMers {
		kcode, err = NewKmerCode(mer) // encode
		if err != nil {
			t.Errorf("Encode error: %s", mer)
		}

		if !bytes.Equal(mer, kcode.Bytes()) { // decode
			t.Errorf("Decode error: %s != %s ", mer, kcode.Bytes())
		}
	}
}

// TestEncodeFromFormerKmer tests TestEncodeFromFormerKmer
func TestEncodeFromFormerKmer(t *testing.T) {
	var err error
	k := 5
	first := true
	var code, code0, pCode uint64
	var kmer, pKmer []byte
	for i := 0; i < len(benchMer)-k; i++ {
		kmer = benchMer[i : i+k]
		if first {
			code, err = Encode(kmer)
			if err != nil {
				t.Errorf("Encode error: %s", kmer)
			}

			pCode = code
			first = false
			continue
		}
		pKmer = benchMer[i-1 : i+k-1]
		code, err = EncodeFromFormerKmer(kmer, pKmer, pCode)
		if err != nil {
			t.Errorf("Encode error: %s", kmer)
		}

		code0, err = Encode(kmer)
		if err != nil {
			t.Errorf("Encode error: %s", kmer)
		}
		if code0 != code {
			t.Errorf("EncodeFromFormerKmer error for %s: wrong %d != right %d", kmer, code, code0)
		}

		pCode = code
	}
}

func TestEncodeFromLatterKmer(t *testing.T) {
	var err error
	k := 5
	first := true
	var code, code0, pCode uint64
	var kmer, pKmer []byte
	for i := len(benchMer) - k - 1; i >= 0; i-- {
		kmer = benchMer[i : i+k]
		if first {
			code, err = Encode(kmer)
			if err != nil {
				t.Errorf("Encode error: %s", kmer)
			}

			pCode = code
			first = false
			continue
		}
		pKmer = benchMer[i+1 : i+k+1]
		code, err = EncodeFromLatterKmer(kmer, pKmer, pCode)
		if err != nil {
			t.Errorf("Encode error: %s", kmer)
		}

		code0, err = Encode(kmer)
		if err != nil {
			t.Errorf("Encode error: %s", kmer)
		}
		if code0 != code {
			t.Errorf("EncodeFromLatterKmer error for %s: wrong %d != right %d", kmer, code, code0)
		}

		pCode = code
	}
}

// TestRevComp tests revcomp
func TestRevComp(t *testing.T) {
	var kcode KmerCode
	for _, mer := range randomMers {
		kcode, _ = NewKmerCode(mer)

		// fmt.Printf("%s, rev:%s\n", kcode, kcode.Rev())

	}

	for _, mer := range randomMers {
		kcode, _ = NewKmerCode(mer)

		if !kcode.Rev().Rev().Equal(kcode) {
			t.Errorf("Rev() error: %s, Rev(): %s", kcode, kcode.Rev())
		}

		if !kcode.Comp().Comp().Equal(kcode) {
			t.Errorf("Comp() error: %s, Comp(): %s", kcode, kcode.Comp())
		}

		if !kcode.Comp().Rev().Equal(kcode.RevComp()) {
			t.Errorf("Rev().Comp() error: %s, Rev(): %s, Comp(): %s, RevComp: %s", kcode, kcode.Rev(), kcode.Comp(), kcode.RevComp())
		}
	}
}

var result uint64

// BenchmarkEncode tests speed of Encode()
func BenchmarkEncodeK32(b *testing.B) {
	var code uint64
	var err error
	for i := 0; i < b.N; i++ {
		code, err = Encode(benchMer)
		if err != nil {
			b.Errorf("Encode error: %s", benchMer)
		}
		if code != codeBenchMer {
			b.Errorf("wrong result: %s", benchMer)
		}
	}
	result = code
}

// BenchmarkEncode tests speed of EncodeFromFormerKmer
func BenchmarkEncodeFromFormerKmerK32(b *testing.B) {
	var code uint64
	var err error
	for i := 0; i < b.N; i++ {
		code, err = EncodeFromFormerKmer(benchMer2, benchMer, benchCode)
		if err != nil {
			b.Errorf("Encode error: %s", benchMer)
		}
		if code != codeBenchMer2 {
			b.Errorf("wrong result: %s", benchMer)
		}
	}
	result = code
}

// BenchmarkEncode tests speed of MustEncodeFromFormerKmer
func BenchmarkMustEncodeFromFormerKmerK32(b *testing.B) {
	var code uint64
	var err error
	for i := 0; i < b.N; i++ {
		code, err = MustEncodeFromFormerKmer(benchMer2, benchMer, benchCode)
		if err != nil {
			b.Errorf("Encode error: %s", benchMer)
		}
		if code != codeBenchMer2 {
			b.Errorf("wrong result: %s", benchMer)
		}
	}
	result = code
}

var result2 []byte

// BenchmarkDecode tests speed of decode
func BenchmarkDecodeK32(b *testing.B) {
	var r []byte
	for i := 0; i < b.N; i++ {
		r = Decode(benchCode, len(benchMer))
	}
	result2 = r
}

func BenchmarkMustDecodeK32(b *testing.B) {
	var r []byte
	for i := 0; i < b.N; i++ {
		r = MustDecode(benchCode, len(benchMer))
	}
	result2 = r
}

var result3 KmerCode

// BenchmarkRevK32 tests speed of rev
func BenchmarkRevK32(b *testing.B) {
	var r KmerCode
	for i := 0; i < b.N; i++ {
		r = benchKmerCode.Rev()
	}
	result3 = r
}

// BenchmarkRevK32 tests speed of comp
func BenchmarkCompK32(b *testing.B) {
	var r KmerCode
	for i := 0; i < b.N; i++ {
		r = benchKmerCode.Comp()
	}
	result3 = r
}

// BenchmarkRevCompK32 tests speed of revcomp
func BenchmarkRevCompK32(b *testing.B) {
	var r KmerCode
	for i := 0; i < b.N; i++ {
		r = benchKmerCode.RevComp()
	}
	result3 = r
}

func BenchmarkCannonalK32(b *testing.B) {
	var r KmerCode
	for i := 0; i < b.N; i++ {
		r = benchKmerCode.Canonical()
	}
	result3 = r
}
