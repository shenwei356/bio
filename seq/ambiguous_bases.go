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
	'A': {'A'},
	'a': {'A'},
	'C': {'C'},
	'c': {'C'},
	'G': {'G'},
	'g': {'G'},
	'T': {'T'},
	't': {'T'},
	'U': {'T'},
	'u': {'T'},

	'M': {'A', 'C', 'M'},
	'm': {'A', 'C', 'M'},
	'R': {'A', 'G', 'R'},
	'r': {'A', 'G', 'R'},
	'W': {'A', 'T', 'W'},
	'w': {'A', 'T', 'W'},
	'S': {'C', 'G', 'S'},
	's': {'C', 'G', 'S'},
	'Y': {'C', 'T', 'Y'},
	'y': {'C', 'T', 'Y'},
	'K': {'G', 'T', 'K'},
	'k': {'G', 'T', 'K'},

	'V': {'A', 'C', 'G', 'M', 'R', 'S', 'V'},
	'v': {'A', 'C', 'G', 'M', 'R', 'S', 'V'},
	'H': {'A', 'C', 'T', 'M', 'W', 'Y', 'H'},
	'h': {'A', 'C', 'T', 'M', 'W', 'Y', 'H'},
	'D': {'A', 'G', 'T', 'R', 'W', 'K', 'D'},
	'd': {'A', 'G', 'T', 'R', 'W', 'K', 'D'},
	'B': {'C', 'G', 'T', 'S', 'Y', 'K', 'B'},
	'b': {'C', 'G', 'T', 'S', 'Y', 'K', 'B'},

	'N': {'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'},
	'n': {'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'},
}

// AmbCodes2Codes is code version of AmbBase2Bases
var AmbCodes2Codes = map[int][]int{
	1: {1},
	2: {2},
	4: {4},
	8: {8},

	3:  {1, 2, 3},
	5:  {1, 4, 5},
	9:  {1, 8, 9},
	6:  {2, 4, 6},
	10: {2, 8, 10},
	12: {4, 8, 12},

	7:  {1, 2, 4, 3, 5, 6, 7},
	11: {1, 2, 8, 3, 9, 10, 11},
	13: {1, 4, 8, 5, 9, 12, 13},
	14: {2, 4, 8, 6, 10, 12, 14},

	15: {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
}
