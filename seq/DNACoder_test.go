package seq

import (
	"testing"
)

func TestDNACoder(t *testing.T) {
	coder, err := NewDNACoder(DNAredundant.Letters())
	if err != nil {
		t.Error(err)
	}

	dna2int, err := coder.Encode([]byte("Jj"))
	if err != ErrInvalideLetter {
		t.Error(err)
	}

	dna2int, err = coder.Encode([]byte("acTg"))
	if err != nil {
		t.Error(err)
	}
	int2dna, err := coder.Decode(dna2int)
	if err != nil {
		t.Error(err)
	}

	if string(int2dna) != "acTg" {
		t.Errorf("DNACoder test error")
	}

}
