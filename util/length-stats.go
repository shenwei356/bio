package util

import (
	"math"

	"github.com/twotwotwo/sorts"
)

// LengthStats accepts sequence lengths and calculate N50 et al..
type LengthStats struct {
	lens     map[uint64]uint64
	sum      uint64
	count    uint64
	min, max uint64

	counts    [][2]uint64 // length and it's number
	accCounts [][2]uint64 // accumulative count
	// accLens   [][2]uint64 // accumulative length
	sorted bool

	n50Calculated bool
	l50           int

	nXCalculated bool
	lX           int
}

// NewLengthStats initializes a LengthStats
func NewLengthStats() *LengthStats {
	return &LengthStats{lens: make(map[uint64]uint64, 256), min: 1<<64 - 1}
}

// Add adds a new length
func (stats *LengthStats) Add(length uint64) {
	stats.count++
	stats.sum += length

	stats.lens[length]++

	if length > stats.max {
		stats.max = length
	}
	if length < stats.min {
		stats.min = length
	}

	stats.sorted = false
}

type lengthCount [][2]uint64

func (lc lengthCount) Len() int           { return len(lc) }
func (lc lengthCount) Less(i, j int) bool { return lc[i][0] < lc[j][0] }
func (lc lengthCount) Swap(i, j int)      { lc[i], lc[j] = lc[j], lc[i] }

func (stats *LengthStats) sort() {
	if len(stats.lens) == 0 {
		stats.sorted = true
		return
	}

	stats.counts = make([][2]uint64, 0, len(stats.lens))
	for length, count := range stats.lens {
		stats.counts = append(stats.counts, [2]uint64{length, count})
	}
	// sort.Slice(stats.counts, func(i, j int) bool { return stats.counts[i][0] < stats.counts[j][0] })
	sorts.Quicksort(lengthCount(stats.counts))

	stats.accCounts = make([][2]uint64, len(stats.lens))

	stats.accCounts[0] = [2]uint64{stats.counts[0][0], stats.counts[0][1]}
	if len(stats.counts) > 1 {
		for i, data := range stats.counts[1:] {
			stats.accCounts[i+1] = [2]uint64{data[0], data[1] + stats.accCounts[i][1]}
		}
	}

	// stats.accLens = make([][2]uint64, len(stats.lens))
	// for i, data := range stats.counts {
	// 	if i == 0 {
	// 		stats.accLens[i] = [2]uint64{data[0], data[0] * data[1]}
	// 	} else {
	// 		stats.accLens[i] = [2]uint64{data[0], data[0]*data[1] + stats.accLens[i-1][1]}
	// 	}
	// }
	stats.sorted = true
}

// Count returns number of elements
func (stats *LengthStats) Count() uint64 {
	return stats.count
}

// Min returns the minimum length
func (stats *LengthStats) Min() uint64 {
	if stats.count == 0 {
		return 0
	}
	return stats.min
}

// Max returns the maxinimum length
func (stats *LengthStats) Max() uint64 {
	return stats.max
}

// Sum returns the length sum
func (stats *LengthStats) Sum() uint64 {
	return stats.sum
}

// Mean returns mean
func (stats *LengthStats) Mean() float64 {
	return float64(stats.sum) / float64(stats.count)
}

// Q2 returns Q2
func (stats *LengthStats) Q2() float64 {
	return stats.Median()
}

// Median returns median
func (stats *LengthStats) Median() float64 {
	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		return float64(stats.counts[0][0])
	}

	even := stats.count&1 == 0    // %2 == 0
	var iMedianL, iMedianR uint64 // 0-based
	if even {
		iMedianL = uint64(stats.count/2) - 1 // 3
		iMedianR = uint64(stats.count / 2)   // 4
	} else {
		iMedianL = uint64(stats.count / 2)
	}

	return stats.getValue(even, iMedianL, iMedianR)
}

// Q1 returns Q1
func (stats *LengthStats) Q1() float64 {
	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		return float64(stats.counts[0][0])
	}

	even := stats.count&1 == 0    // %2 == 0
	var iMedianL, iMedianR uint64 // 0-based
	var n uint64
	if even {
		n = stats.count / 2
	} else {
		n = (stats.count + 1) / 2
	}

	even = n%2 == 0
	if even {
		iMedianL = uint64(n/2) - 1
		iMedianR = uint64(n / 2)
	} else {
		iMedianL = uint64(n / 2)
	}

	return stats.getValue(even, iMedianL, iMedianR)
}

// Q3 returns Q3
func (stats *LengthStats) Q3() float64 {
	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		return float64(stats.counts[0][0])
	}

	even := stats.count&1 == 0    // %2 == 0
	var iMedianL, iMedianR uint64 // 0-based
	var mean, n uint64
	if even {
		n = stats.count / 2
		mean = n
	} else {
		n = (stats.count + 1) / 2
		mean = stats.count / 2
	}

	even = n%2 == 0
	if even {
		iMedianL = uint64(n/2) - 1 + mean
		iMedianR = uint64(n/2) + mean
	} else {
		iMedianL = uint64(n/2) + mean
	}

	return stats.getValue(even, iMedianL, iMedianR)
}

func (stats *LengthStats) getValue(even bool, iMedianL uint64, iMedianR uint64) float64 {

	var accCount uint64

	var flag bool
	var prev uint64

	for _, data := range stats.accCounts {
		accCount = data[1]

		if flag {
			// the middle two having different length.
			// example: 1, 2, 3, 4 or 1, 2
			return float64(data[0]+prev) / 2
		}

		if accCount >= iMedianL+1 {
			if even {
				if accCount >= iMedianR+1 {
					// having >=2 of same length in the middle.
					// example: 2, 2, 2, 3, 3, 4, 8, 8
					return float64(data[0])
				}
				flag = true
				prev = data[0]
			} else {
				// right here
				return float64(data[0])
			}
		}
	}

	// never happen
	// panic("bio/util: should never happen")
	return 0
}

func (stats *LengthStats) Percentile(percent float64) float64 {
	if percent <= 0 || percent > 100 {
		panic("invalid percentile")
	}
	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		return float64(stats.counts[0][0])
	}

	i0 := float64(stats.count) * percent / 100
	i := math.Floor(i0)

	even := math.Abs(i0-i) > 0.001
	var iMedianL, iMedianR uint64 // 0-based
	if even {
		iMedianL = uint64(i) - 1
		iMedianR = uint64(i)
	} else {
		iMedianL = uint64(i - 1)
	}

	return stats.getValue(even, iMedianL, iMedianR)
}

// N50 returns N50
func (stats *LengthStats) N50() uint64 {
	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		stats.l50 = 1
		return stats.counts[0][0]
	}

	var sumLen float64
	var data [2]uint64
	var halfSum = float64(stats.sum) / 2
	for i := len(stats.counts) - 1; i >= 0; i-- {
		data = stats.counts[i]

		sumLen += float64(data[0] * data[1])
		if sumLen >= halfSum {
			stats.l50 = len(stats.counts) - i
			stats.n50Calculated = true
			return data[0]
		}
	}

	// never happen
	// panic("bio/util: should never happen")
	return 0
}

// NX returns something like N50, where X could be a number in the range of [0, 100]
func (stats *LengthStats) NX(n float64) uint64 {
	if n < 0 || n > 100 {
		panic("NX: where X should be in range of [0, 100]")
	}

	if !stats.sorted {
		stats.sort()
	}
	if len(stats.counts) == 0 {
		return 0
	}

	if len(stats.counts) == 1 {
		return stats.counts[0][0]
	}

	var sumLen float64
	var data [2]uint64
	var boundary = float64(stats.sum) * n / 100
	for i := len(stats.counts) - 1; i >= 0; i-- {
		data = stats.counts[i]

		sumLen += float64(data[0] * data[1])
		if sumLen >= boundary {
			stats.lX = i + 1
			stats.nXCalculated = true
			return data[0]
		}
	}

	// never happen
	// panic("bio/util: should never happen")
	return 0
}

// L50 returns L50
func (stats *LengthStats) L50() int {
	if !stats.sorted {
		stats.N50()
	}
	if !stats.n50Calculated {
		stats.N50()
	}
	return stats.l50
}

// LX returns something like L50, where X could be a number in the range of [0, 100]
func (stats *LengthStats) LX(n float64) int {
	if !stats.sorted {
		stats.NX(n)
	}
	if !stats.nXCalculated {
		stats.NX(n)
	}
	return stats.lX
}
