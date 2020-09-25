package util

import (
	"math/rand"
	"testing"
)

type testCaseLengthStats struct {
	data []uint64

	q1, median, q3 float64
	n50            uint64
	min, max       uint64
}

var cases = []testCaseLengthStats{
	testCaseLengthStats{
		data:   []uint64{},
		median: 0,
		q1:     0,
		q3:     0,
		n50:    0,
		min:    0,
		max:    0,
	},
	testCaseLengthStats{
		data:   []uint64{2},
		median: 2,
		q1:     1,
		q3:     1,
		n50:    2,
		min:    2,
		max:    2,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2},
		median: 1.5,
		q1:     1,
		q3:     2,
		n50:    2,
		min:    1,
		max:    2,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2, 3},
		median: 2,
		q1:     1.5,
		q3:     2.5,
		n50:    3,
		min:    1,
		max:    3,
	},
	testCaseLengthStats{
		data:   []uint64{1, 2, 3, 4},
		median: 2.5,
		q1:     1.5,
		q3:     3.5,
		n50:    3,
		min:    1,
		max:    4,
	},

	testCaseLengthStats{
		data:   []uint64{2, 3, 4, 5, 6, 7, 8, 9},
		median: 5.5,
		q1:     3.5,
		q3:     7.5,
		n50:    7,
		min:    2,
		max:    9,
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
		if q1 != _case.q1 {
			t.Errorf("case %d: q3 mismatch: %f != %f", i, q3, _case.q3)
		}

		n50 := stats.N50()
		if n50 != _case.n50 {
			t.Errorf("case %d: n50 mismatch: %d != %d", i, n50, _case.n50)
		}
	}
}
