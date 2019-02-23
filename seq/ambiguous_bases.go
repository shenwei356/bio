package seq

import "errors"

/*
	base    bases   code
	A       A       1
	C       C       2
	G       G       4
	T/U     T       8
	M       A/C     3
	R       A/G     5
	W       A/T     9
	S       C/G     6
	Y       C/T     10
	K       G/T     12
	V       A/C/G   7
	H       A/C/T   11
	D       A/G/T   13
	B       C/G/T   14
	N       A/C/G/T 15
*/
var code2base = [16]byte{'-', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'}

// ErrInvalidDNABase means invalid DNA base
var ErrInvalidDNABase = errors.New("seq: invalid DNA base")

func base2code(b byte) (int, error) {
	var c int
	switch b {
	case 'A', 'a':
		c = 1
	case 'C', 'c':
		c = 2
	case 'G', 'g':
		c = 4
	case 'T', 't', 'U', 'u':
		c = 8
	case 'N', 'n':
		c = 15
	case 'M', 'm':
		c = 3
	case 'R', 'r':
		c = 5
	case 'W', 'w':
		c = 9
	case 'S', 's':
		c = 6
	case 'Y', 'y':
		c = 10
	case 'K', 'k':
		c = 12
	case 'V', 'v':
		c = 7
	case 'H', 'h':
		c = 11
	case 'D', 'd':
		c = 13
	case 'B', 'b':
		c = 14
	case ' ', '*', '-':
		c = 0
	default:
		return 0, ErrInvalidDNABase
	}
	return c, nil
}

// Bases2AmbBase converts list of bases to ambiguous base
func Bases2AmbBase(bs []byte) (byte, error) {
	var code, c int
	var err error
	for _, b := range bs {
		c, err = base2code(b)
		if err != nil {
			return 0, err
		}
		code |= c
	}
	return code2base[code], nil
}

// Codes2AmbCode converts list of codes of bases to code of ambiguous base
func Codes2AmbCode(codes []int) (int, error) {
	var code int
	for _, c := range codes {
		code |= c
	}
	return code, nil
}

// AmbBase2Bases0 converts ambiguous base to bases it represents, slower than AmbBase2Bases
func AmbBase2Bases0(b byte) ([]byte, error) {
	code, err := base2code(b)
	if err != nil {
		return nil, err
	}

	bases := []byte{}
	var c int
	for comb := 1; comb <= code; comb++ {
		c = 0
		for i := uint(0); i < 4; i++ {
			if (comb>>i)&1 == 1 && // this combination needs this bit being 1
				code&(1<<i) > 0 { // this bit is 1
				c += 1 << i
			}
		}
		if c == comb { // all bits are set
			bases = append(bases, code2base[c])
		}
	}
	return bases, nil
}

// AmbBase2Bases holds relationship of ambiguous base and bases it represents, faster than AmbBase2Bases0
var AmbBase2Bases = map[byte][]byte{
	'A': []byte{'A'},
	'a': []byte{'A'},
	'C': []byte{'C'},
	'c': []byte{'C'},
	'G': []byte{'G'},
	'g': []byte{'G'},
	'T': []byte{'T'},
	't': []byte{'T'},
	'U': []byte{'T'},
	'u': []byte{'T'},

	'M': []byte{'A', 'C', 'M'},
	'm': []byte{'A', 'C', 'M'},
	'R': []byte{'A', 'G', 'R'},
	'r': []byte{'A', 'G', 'R'},
	'W': []byte{'A', 'T', 'W'},
	'w': []byte{'A', 'T', 'W'},
	'S': []byte{'C', 'G', 'S'},
	's': []byte{'C', 'G', 'S'},
	'Y': []byte{'C', 'T', 'Y'},
	'y': []byte{'C', 'T', 'Y'},
	'K': []byte{'G', 'T', 'K'},
	'k': []byte{'G', 'T', 'K'},

	'V': []byte{'A', 'C', 'G', 'M', 'R', 'S', 'V'},
	'v': []byte{'A', 'C', 'G', 'M', 'R', 'S', 'V'},
	'H': []byte{'A', 'C', 'T', 'M', 'W', 'Y', 'H'},
	'h': []byte{'A', 'C', 'T', 'M', 'W', 'Y', 'H'},
	'D': []byte{'A', 'G', 'T', 'R', 'W', 'K', 'D'},
	'd': []byte{'A', 'G', 'T', 'R', 'W', 'K', 'D'},
	'B': []byte{'C', 'G', 'T', 'S', 'Y', 'K', 'B'},
	'b': []byte{'C', 'G', 'T', 'S', 'Y', 'K', 'B'},

	'N': []byte{'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'},
	'n': []byte{'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'},
}

// AmbCodes2Codes is code version of AmbBase2Bases
var AmbCodes2Codes = map[int][]int{
	1: []int{1},
	2: []int{2},
	4: []int{4},
	8: []int{8},

	3:  []int{1, 2, 3},
	5:  []int{1, 4, 5},
	9:  []int{1, 8, 9},
	6:  []int{2, 4, 6},
	10: []int{2, 8, 10},
	12: []int{4, 8, 12},

	7:  []int{1, 2, 4, 3, 5, 6, 7},
	11: []int{1, 2, 8, 3, 9, 10, 11},
	13: []int{1, 4, 8, 5, 9, 12, 13},
	14: []int{2, 4, 8, 6, 10, 12, 14},

	15: []int{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
}
