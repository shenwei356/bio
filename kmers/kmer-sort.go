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

// KmerCodeSlice is a slice of KmerCode, for sorting
type KmerCodeSlice []KmerCode

// Len return length of the slice
func (codes KmerCodeSlice) Len() int {
	return len(codes)
}

// Swap swaps two elements
func (codes KmerCodeSlice) Swap(i, j int) {
	codes[i], codes[j] = codes[j], codes[i]
}

// Less simply compare two KmerCode
func (codes KmerCodeSlice) Less(i, j int) bool {
	return codes[i].Code < codes[j].Code
}

// func splitKmer(code uint64, k int) (uint64, uint64, uint64, uint64) {
// 	// -====, k = 4:  ---, -, =, ===
// 	return code >> 2, code & 3, code >> (uint(k-1) << 1) & 3, code & ((1 << (uint(k-1) << 1)) - 1)
// }

// CodeSlice is a slice of Kmer code (uint64), for sorting
type CodeSlice []uint64

// Len return length of the slice
func (codes CodeSlice) Len() int {
	return len(codes)
}

// Swap swaps two elements
func (codes CodeSlice) Swap(i, j int) {
	codes[i], codes[j] = codes[j], codes[i]
}

// Less simply compare two KmerCode
func (codes CodeSlice) Less(i, j int) bool {
	return codes[i] < codes[j]
}
