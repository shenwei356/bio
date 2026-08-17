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

package taxdump

import (
	"os"
	"path/filepath"
	"testing"
)

func TestPackTwoTaxids(t *testing.T) {
	type Test struct {
		a, b uint32
		c    uint64
	}
	tests := []Test{
		{0, 0, 0},
		{1, 1, 1<<32 + 1},
		{2, 1, 1<<32 + 2},
	}

	for _, test := range tests {
		c := pack2uint32(test.a, test.b)
		if c != test.c {
			t.Errorf("pack2uint32 error: %d != %d ", c, test.c)
		}
	}
}

func TestLoadNamesWithTaxIDInHighestColumn(t *testing.T) {
	file := filepath.Join(t.TempDir(), "names.tsv")
	if err := os.WriteFile(file, []byte("scientific name\twanted\tunused\t42\n"), 0o600); err != nil {
		t.Fatal(err)
	}

	taxonomy := &Taxonomy{}
	if err := taxonomy.LoadNames(file, 4, 1, 2, "wanted"); err != nil {
		t.Fatal(err)
	}
	if got := taxonomy.Name(42); got != "scientific name" {
		t.Fatalf("loaded name %q, want %q", got, "scientific name")
	}
}
