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
		Test{[]byte{'A'}, 'A', nil},
		Test{[]byte{'c'}, 'C', nil},
		Test{[]byte{'t'}, 'T', nil},
		Test{[]byte{'u'}, 'T', nil},
		Test{[]byte{'g'}, 'G', nil},

		Test{[]byte{'a', 'c'}, 'M', nil},
		Test{[]byte{'c', 'a'}, 'M', nil},
		Test{[]byte{'a', 'g'}, 'R', nil},
		Test{[]byte{'g', 'a'}, 'R', nil},
		Test{[]byte{'a', 't'}, 'W', nil},
		Test{[]byte{'t', 'a'}, 'W', nil},
		Test{[]byte{'c', 'g'}, 'S', nil},
		Test{[]byte{'g', 'c'}, 'S', nil},
		Test{[]byte{'c', 't'}, 'Y', nil},
		Test{[]byte{'t', 'c'}, 'Y', nil},
		Test{[]byte{'g', 't'}, 'K', nil},
		Test{[]byte{'t', 'g'}, 'K', nil},

		Test{[]byte{'a', 'c', 'g'}, 'V', nil},
		Test{[]byte{'a', 'c', 't'}, 'H', nil},
		Test{[]byte{'a', 'G', 't'}, 'D', nil},
		Test{[]byte{'C', 'g', 't'}, 'B', nil},
		Test{[]byte{'C', 'g', 't', 'a'}, 'N', nil},

		Test{[]byte{'j'}, 'N', ErrInvalidDNABase},
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
	_AmbBase2BasesTest{'A', []byte{'A'}, nil},
	_AmbBase2BasesTest{'C', []byte{'C'}, nil},
	_AmbBase2BasesTest{'G', []byte{'G'}, nil},
	_AmbBase2BasesTest{'T', []byte{'T'}, nil},

	_AmbBase2BasesTest{'M', []byte{'A', 'C', 'M'}, nil},
	_AmbBase2BasesTest{'R', []byte{'A', 'G', 'R'}, nil},
	_AmbBase2BasesTest{'W', []byte{'A', 'T', 'W'}, nil},
	_AmbBase2BasesTest{'S', []byte{'C', 'G', 'S'}, nil},
	_AmbBase2BasesTest{'Y', []byte{'C', 'T', 'Y'}, nil},
	_AmbBase2BasesTest{'K', []byte{'G', 'T', 'K'}, nil},

	_AmbBase2BasesTest{'V', []byte{'A', 'C', 'G', 'M', 'R', 'S', 'V'}, nil},
	_AmbBase2BasesTest{'H', []byte{'A', 'C', 'T', 'M', 'W', 'Y', 'H'}, nil},
	_AmbBase2BasesTest{'D', []byte{'A', 'G', 'T', 'R', 'W', 'K', 'D'}, nil},
	_AmbBase2BasesTest{'B', []byte{'C', 'G', 'T', 'S', 'Y', 'K', 'B'}, nil},

	_AmbBase2BasesTest{'N', []byte{'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}, nil},
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
