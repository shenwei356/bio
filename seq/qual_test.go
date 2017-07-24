package seq

import (
	"math"
	"testing"
)

func TestPhredAndSolexa(t *testing.T) {
	var q2, q3 float64
	var err error
	for q := 2.0; q <= 40; q++ {
		q2, err = Phred2Solexa(q)
		if err != nil {
			t.Error(err)
		}
		q3, err = Solexa2Phred(q2)
		if err != nil {
			t.Error(err)
		}

		if math.Abs(q3-q) > 0.01 {
			t.Errorf("%.2f, %.2f, %.2f", q, q2, q3)
		}
	}
}
