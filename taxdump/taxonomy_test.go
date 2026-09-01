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
	"math/rand"
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

func TestLCA(t *testing.T) {
	taxonomy := &Taxonomy{
		rootNode: 1,
		Nodes: map[uint32]uint32{
			1: 1,
			2: 1, 3: 1,
			4: 2, 5: 2, 6: 3, 7: 3,
			8: 4,
		},
		MergeNodes: map[uint32]uint32{80: 8, 90: 9},
	}
	tests := []struct {
		a, b uint32
		want uint32
	}{
		{0, 8, 0},
		{8, 8, 8},
		{99, 99, 99},
		{4, 5, 2},
		{8, 5, 2},
		{8, 2, 2},
		{2, 8, 2},
		{8, 6, 1},
		{8, 99, 0},
		{80, 8, 4},
		{8, 80, 4},
		{80, 4, 4},
		{80, 5, 2},
		{80, 90, 0},
	}
	for _, test := range tests {
		if got := taxonomy.LCA(test.a, test.b); got != test.want {
			t.Errorf("LCA(%d, %d) = %d, want %d", test.a, test.b, got, test.want)
		}
	}

	taxonomy.CacheLCA()
	for _, test := range tests {
		for i := 0; i < 2; i++ {
			if got := taxonomy.LCA(test.a, test.b); got != test.want {
				t.Errorf("cached LCA(%d, %d) = %d, want %d", test.a, test.b, got, test.want)
			}
		}
	}
}

func TestLCAAgainstReference(t *testing.T) {
	rng := rand.New(rand.NewSource(7))
	for tree := 0; tree < 20; tree++ {
		nodes := map[uint32]uint32{1: 1}
		for id := uint32(2); id <= 500; id++ {
			nodes[id] = uint32(rng.Intn(int(id-1))) + 1
		}
		merged := make(map[uint32]uint32, 50)
		for id := uint32(1001); id <= 1050; id++ {
			merged[id] = uint32(rng.Intn(500)) + 1
		}
		taxonomy := &Taxonomy{rootNode: 1, Nodes: nodes, MergeNodes: merged}

		for query := 0; query < 5000; query++ {
			a := randomTaxid(rng)
			b := randomTaxid(rng)
			got := taxonomy.LCA(a, b)
			want := referenceLCA(taxonomy, a, b)
			if got != want {
				t.Fatalf("tree %d: LCA(%d, %d) = %d, want %d", tree, a, b, got, want)
			}
		}
	}
}

func TestLCAFallsBackForDeepLineages(t *testing.T) {
	nodes := map[uint32]uint32{1: 1}
	for id := uint32(2); id <= lcaPathBufferSize*2; id++ {
		nodes[id] = id - 1
	}
	taxonomy := &Taxonomy{rootNode: 1, Nodes: nodes}
	if got, want := taxonomy.LCA(lcaPathBufferSize*2, lcaPathBufferSize), uint32(lcaPathBufferSize); got != want {
		t.Fatalf("deep LCA = %d, want %d", got, want)
	}
	if allocs := testing.AllocsPerRun(1000, func() {
		taxonomy.LCA(lcaPathBufferSize*2, lcaPathBufferSize)
	}); allocs != 0 {
		t.Fatalf("deep LCA allocated %.1f times per call", allocs)
	}
}

func TestLCAPreservesRootFallbackForForests(t *testing.T) {
	taxonomy := &Taxonomy{
		rootNode: 2,
		Nodes: map[uint32]uint32{
			1: 1, 2: 2,
			3: 1, 4: 2,
		},
	}
	if got, want := taxonomy.LCA(3, 4), uint32(2); got != want {
		t.Fatalf("forest LCA = %d, want %d", got, want)
	}

	deepNodes := map[uint32]uint32{1: 1, 2: 2}
	for id := uint32(3); id <= lcaPathBufferSize+2; id++ {
		deepNodes[id] = id - 1
	}
	deep := &Taxonomy{rootNode: 2, Nodes: deepNodes}
	if got, want := deep.LCA(lcaPathBufferSize+2, 1), uint32(2); got != want {
		t.Fatalf("deep forest LCA = %d, want %d", got, want)
	}
}

func BenchmarkLCA(b *testing.B) {
	nodes := make(map[uint32]uint32, 1<<20)
	nodes[1] = 1
	for id := uint32(2); id <= 1<<20; id++ {
		nodes[id] = id >> 1
	}
	pairs := make([][2]uint32, 1024)
	rng := rand.New(rand.NewSource(7))
	for i := range pairs {
		pairs[i] = [2]uint32{uint32(rng.Intn(1<<20)) + 1, uint32(rng.Intn(1<<20)) + 1}
	}

	b.Run("uncached", func(b *testing.B) {
		taxonomy := &Taxonomy{rootNode: 1, Nodes: nodes}
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			pair := pairs[i&(len(pairs)-1)]
			_ = taxonomy.LCA(pair[0], pair[1])
		}
	})
	b.Run("cached", func(b *testing.B) {
		taxonomy := &Taxonomy{rootNode: 1, Nodes: nodes}
		taxonomy.CacheLCA()
		for _, pair := range pairs {
			taxonomy.LCA(pair[0], pair[1])
		}
		b.ReportAllocs()
		b.ResetTimer()
		for i := 0; i < b.N; i++ {
			pair := pairs[i&(len(pairs)-1)]
			_ = taxonomy.LCA(pair[0], pair[1])
		}
	})
}

func randomTaxid(rng *rand.Rand) uint32 {
	switch rng.Intn(10) {
	case 0:
		return 0
	case 1:
		return uint32(rng.Intn(50)) + 1001
	case 2:
		return uint32(rng.Intn(50)) + 2001
	default:
		return uint32(rng.Intn(500)) + 1
	}
}

func referenceLCA(t *Taxonomy, a uint32, b uint32) uint32 {
	if a == 0 || b == 0 {
		return 0
	}
	if a == b {
		return a
	}

	ancestors := make(map[uint32]struct{}, 16)
	child := a
	for {
		parent, ok := t.Nodes[child]
		if !ok {
			newTaxid, merged := t.MergeNodes[child]
			if !merged {
				return 0
			}
			child = newTaxid
			parent, ok = t.Nodes[child]
			if !ok {
				return 0
			}
		}
		if parent == child {
			ancestors[parent] = struct{}{}
			break
		}
		if parent == b {
			return b
		}
		ancestors[parent] = struct{}{}
		child = parent
	}

	child = b
	for {
		parent, ok := t.Nodes[child]
		if !ok {
			newTaxid, merged := t.MergeNodes[child]
			if !merged {
				return 0
			}
			child = newTaxid
			parent, ok = t.Nodes[child]
			if !ok {
				return 0
			}
		}
		if parent == child {
			break
		}
		if parent == a {
			return a
		}
		if _, ok := ancestors[parent]; ok {
			return parent
		}
		child = parent
	}
	return t.rootNode
}
