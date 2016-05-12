package seq

import (
	"fmt"
	"github.com/shenwei356/util/byteutil"
	"sort"
)

var (
	// ErrInvalideLetter means invalid letter
	ErrInvalideLetter = fmt.Errorf("DNACoder: invalid letter")
	// ErrInvalideCode means invalid code
	ErrInvalideCode = fmt.Errorf("DNACoder: invalid code")
)

// DNACoder is used to convert betweeen byte and int
type DNACoder struct {
	Alphabet           []byte
	alphabetQuerySlice []byte
	dna2int            []int
	int2dna            []byte
}

// NewDNACoder Create a DNACoder type
func NewDNACoder(alphabet []byte) (*DNACoder, error) {
	if len(alphabet) == 0 {
		return nil, fmt.Errorf("dnacoder: alphabet should not be empty")
	}

	m := make(map[byte]struct{}, len(alphabet))
	for _, a := range alphabet {
		m[a] = struct{}{}
	}

	max := -1
	var b int
	for a := range m {
		b = int(a)
		if max < b {
			max = b
		}
	}

	alphabet2 := make([]byte, len(m))
	slice := make([]byte, max+1)
	i := 0
	for a := range m {
		slice[a-'\x00'] = a

		alphabet2[i] = a
		i++
	}

	sort.Sort(byteutil.ByteSlice(alphabet2))

	dna2int := make([]int, max+1)
	int2dna := make([]byte, len(m))
	for i, a := range alphabet2 {
		dna2int[a-'\x00'] = i
		int2dna[i] = a
	}

	return &DNACoder{Alphabet: alphabet2, alphabetQuerySlice: slice,
		dna2int: dna2int, int2dna: int2dna}, nil
}

func (coder *DNACoder) String() string {
	return fmt.Sprintf(`DNACoder: alphabet:"%s" num:%d`, coder.Alphabet, len(coder.Alphabet))
}

// Encode converts []byte to []int
func (coder *DNACoder) Encode(s []byte) ([]int, error) {
	code := make([]int, len(s))
	for i, b := range s {
		if int(b) > len(coder.alphabetQuerySlice) {
			return nil, ErrInvalideLetter
		}
		v := coder.alphabetQuerySlice[b-'\x00']
		if v == 0 {
			return nil, ErrInvalideLetter
		}
		code[i] = coder.dna2int[v]
	}
	return code, nil
}

// Decode convert []int to []byte
func (coder *DNACoder) Decode(code []int) ([]byte, error) {
	dna := make([]byte, len(code))
	for i, b := range code {
		if b >= len(coder.int2dna) {
			return nil, ErrInvalideCode
		}
		v := coder.int2dna[b]
		if v == 0 {
			return nil, ErrInvalideCode
		}
		dna[i] = v
	}
	return dna, nil
}
