package seq

import (
	"bytes"
	"sort"
	"testing"

	"github.com/cznic/sortutil"
)

func TestBases2AmbBase(t *testing.T) {
	type Test struct {
		bases  []byte
		ambase byte
		err    error
	}

	tests := []Test{
		{[]byte{'A'}, 'A', nil},
		{[]byte{'c'}, 'C', nil},
		{[]byte{'t'}, 'T', nil},
		{[]byte{'u'}, 'T', nil},
		{[]byte{'g'}, 'G', nil},

		{[]byte{'a', 'c'}, 'M', nil},
		{[]byte{'c', 'a'}, 'M', nil},
		{[]byte{'a', 'g'}, 'R', nil},
		{[]byte{'g', 'a'}, 'R', nil},
		{[]byte{'a', 't'}, 'W', nil},
		{[]byte{'t', 'a'}, 'W', nil},
		{[]byte{'c', 'g'}, 'S', nil},
		{[]byte{'g', 'c'}, 'S', nil},
		{[]byte{'c', 't'}, 'Y', nil},
		{[]byte{'t', 'c'}, 'Y', nil},
		{[]byte{'g', 't'}, 'K', nil},
		{[]byte{'t', 'g'}, 'K', nil},

		{[]byte{'a', 'c', 'g'}, 'V', nil},
		{[]byte{'a', 'c', 't'}, 'H', nil},
		{[]byte{'a', 'G', 't'}, 'D', nil},
		{[]byte{'C', 'g', 't'}, 'B', nil},
		{[]byte{'C', 'g', 't', 'a'}, 'N', nil},

		{[]byte{'j'}, 'N', ErrInvalidDNABase},
	}
	var amb byte
	var err error
	for _, test := range tests {
		amb, err = Bases2AmbBase(test.bases)
		if err != nil {
			if err != ErrInvalidDNABase {
				t.Errorf("error should be: %s, got: %s", ErrInvalidDNABase, err)
			} else {
				continue
			}
		}
		if amb != test.ambase {
			t.Errorf("Test Bases2AmbBase Err: %s, need: %c, got: %c\n", test.bases, test.ambase, amb)
		}
	}
}

type _AmbBase2BasesTest struct {
	ambase byte
	bases  []byte
	err    error
}

var tests4AmbBase2Bases = []_AmbBase2BasesTest{
	{'A', []byte{'A'}, nil},
	{'C', []byte{'C'}, nil},
	{'G', []byte{'G'}, nil},
	{'T', []byte{'T'}, nil},

	{'M', []byte{'A', 'C', 'M'}, nil},
	{'R', []byte{'A', 'G', 'R'}, nil},
	{'W', []byte{'A', 'T', 'W'}, nil},
	{'S', []byte{'C', 'G', 'S'}, nil},
	{'Y', []byte{'C', 'T', 'Y'}, nil},
	{'K', []byte{'G', 'T', 'K'}, nil},

	{'V', []byte{'A', 'C', 'G', 'M', 'R', 'S', 'V'}, nil},
	{'H', []byte{'A', 'C', 'T', 'M', 'W', 'Y', 'H'}, nil},
	{'D', []byte{'A', 'G', 'T', 'R', 'W', 'K', 'D'}, nil},
	{'B', []byte{'C', 'G', 'T', 'S', 'Y', 'K', 'B'}, nil},

	{'N', []byte{'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}, nil},
}

func TestAmbBase2Bases(t *testing.T) {
	var bases []byte
	var err error
	var b1, b2 sortutil.ByteSlice
	for _, test := range tests4AmbBase2Bases {
		bases, err = AmbBase2Bases0(test.ambase)
		if err != nil {
			if err != ErrInvalidDNABase {
				t.Errorf("error should be: %s, got: %s", ErrInvalidDNABase, err)
			} else {
				continue
			}
		}
		b1 = sortutil.ByteSlice(bases)
		b2 = sortutil.ByteSlice(test.bases)
		sort.Sort(b1)
		sort.Sort(b2)
		if bytes.Compare(b1, b2) != 0 {
			t.Errorf("Test AmbBase2Bases0 Err: %c, need: %s, got: %s\n", test.ambase, b2, b1)
		}
	}
}

func BenchmarkAmbBase2Bases0(b *testing.B) {
	for i := 0; i < b.N; i++ {
		var err error
		result := [][]byte{}

		for _, test := range tests4AmbBase2Bases {
			var bases []byte
			bases, err = AmbBase2Bases0(test.ambase)
			if err != nil {
				b.Errorf("Bench AmbBase2Bases0 Err")
			}
			result = append(result, bases)
		}
	}
}

func BenchmarkAmbBase2Bases(b *testing.B) {
	for i := 0; i < b.N; i++ {
		var ok bool
		result := [][]byte{}

		for _, test := range tests4AmbBase2Bases {
			var bases []byte
			if bases, ok = AmbBase2Bases[test.ambase]; ok {
				result = append(result, bases)
			}
		}
	}
}
