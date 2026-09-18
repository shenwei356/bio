package util

import (
	"math/rand"
	"testing"
)

type testCaseLengthStats struct {
	data []uint64

	q1, median, q3 float64
	n50            uint64
	l50            int
	min, max       uint64
}

var cases = []testCaseLengthStats{
	testCaseLengthStats{
		data:   []uint64{},
		median: 0,
		q1:     0,
		q3:     0,
		n50:    0,
		l50:    0,
		min:    0,
		max:    0,
	},
	testCaseLengthStats{
		data:   []uint64{2},
		median: 2,
		q1:     2,
		q3:     2,
		n50:    2,
		l50:    1,
		min:    2,
		max:    2,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2},
		median: 1.5,
		q1:     1,
		q3:     2,
		n50:    2,
		l50:    1,
		min:    1,
		max:    2,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2, 3},
		median: 2,
		q1:     1.5,
		q3:     2.5,
		n50:    3,
		l50:    1,
		min:    1,
		max:    3,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2, 3, 4},
		median: 2.5,
		q1:     1.5,
		q3:     3.5,
		n50:    3,
		l50:    2,
		min:    1,
		max:    4,
	},

	testCaseLengthStats{
		data:   []uint64{2, 3, 4, 5, 6, 7, 8, 9},
		median: 5.5,
		q1:     3.5,
		q3:     7.5,
		n50:    7,
		l50:    3,
		min:    2,
		max:    9,
	},

	// All values are the same - odd count
	testCaseLengthStats{
		data:   []uint64{3, 3, 3},
		median: 3,
		q1:     3,
		q3:     3,
		n50:    3,
		l50:    2,
		min:    3,
		max:    3,
	},

	// All values are the same - even count
	testCaseLengthStats{
		data:   []uint64{5, 5, 5, 5},
		median: 5,
		q1:     5,
		q3:     5,
		n50:    5,
		l50:    2,
		min:    5,
		max:    5,
	},

	// All values are the same - larger count
	testCaseLengthStats{
		data:   []uint64{7, 7, 7, 7, 7, 7, 7, 7},
		median: 7,
		q1:     7,
		q3:     7,
		n50:    7,
		l50:    4,
		min:    7,
		max:    7,
	},
}

func Test(t *testing.T) {
	for i, _case := range cases {
		rand.Shuffle(len(_case.data), func(i, j int) {
			_case.data[i], _case.data[j] = _case.data[j], _case.data[i]
		})

		stats := NewLengthStats()
		for _, l := range _case.data {
			stats.Add(l)
		}
		if stats.Count() != uint64(len(_case.data)) {
			t.Errorf("case %d: count mismatch", i)
		}

		min := stats.Min()
		if min != _case.min {
			t.Errorf("case %d: min mismatch: %d != %d", i, min, _case.min)
		}

		max := stats.Max()
		if max != _case.max {
			t.Errorf("case %d: max mismatch: %d != %d", i, max, _case.max)
		}

		median := stats.Median()
		if median != _case.median {
			t.Errorf("case %d: median mismatch: %f != %f", i, median, _case.median)
		}

		q1 := stats.Q1()
		if q1 != _case.q1 {
			t.Errorf("case %d: q1 mismatch: %f != %f", i, q1, _case.q1)
		}

		q3 := stats.Q3()
		if q3 != _case.q3 {
			t.Errorf("case %d: q3 mismatch: %f != %f", i, q3, _case.q3)
		}

		n50 := stats.N50()
		if n50 != _case.n50 {
			t.Errorf("case %d: n50 mismatch: %d != %d", i, n50, _case.n50)
		}

		l50 := stats.L50()
		if l50 != _case.l50 {
			t.Errorf("case %d: l50 mismatch: %d != %d", i, l50, _case.l50)
		}
	}
}

func TestNXLX(t *testing.T) {
	// Test NX and LX with all same values
	stats := NewLengthStats()
	for i := 0; i < 8; i++ {
		stats.Add(10)
	}

	// Sum = 80, N50 boundary = 40
	// Need 4 sequences of length 10 to reach 40
	n50 := stats.N50()
	if n50 != 10 {
		t.Errorf("N50 should be 10, got %d", n50)
	}
	l50 := stats.L50()
	if l50 != 4 {
		t.Errorf("L50 should be 4, got %d", l50)
	}

	// Test N75 (75% of 80 = 60)
	// Need 6 sequences of length 10 to reach 60
	n75 := stats.NX(75)
	if n75 != 10 {
		t.Errorf("N75 should be 10, got %d", n75)
	}
	l75 := stats.LX(75)
	if l75 != 6 {
		t.Errorf("L75 should be 6, got %d", l75)
	}

	// Test N90 (90% of 80 = 72)
	// Need 8 sequences (ceiling of 72/10 = 7.2)
	n90 := stats.NX(90)
	if n90 != 10 {
		t.Errorf("N90 should be 10, got %d", n90)
	}
	l90 := stats.LX(90)
	if l90 != 8 {
		t.Errorf("L90 should be 8, got %d", l90)
	}
}
