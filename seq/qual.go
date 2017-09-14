package seq

import (
	"errors"
	"math"

	"github.com/shenwei356/util/byteutil"
)

// ErrInvalidPhredQuality occurs for phred quality less than 0.
var ErrInvalidPhredQuality = errors.New("seq: invalid Phred quality")

// ErrInvalidSolexaQuality occurs for solexa quality less than -5.
var ErrInvalidSolexaQuality = errors.New("seq: invalid Solexa quality")

// Phred2Solexa converts Phred quality to Solexa quality.
func Phred2Solexa(q float64) (float64, error) {
	if q == 0 {
		return -5, nil
	}
	if q < 0 {
		return -5, ErrInvalidPhredQuality
	}
	return max(-5, 10*math.Log10(math.Pow(10, q/10)-1)), nil
}

// Solexa2Phred converts Solexa quality to Phred quality.
func Solexa2Phred(q float64) (float64, error) {
	if q < -5 {
		return 0, ErrInvalidSolexaQuality
	}
	return 10 * math.Log10(math.Pow(10, q/10)+1), nil
}

func max(a float64, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

// QualityEncoding is the type of quality encoding
type QualityEncoding int

// NQualityEncoding is the number of QualityEncoding + 1: 5 + 1 = 6
const NQualityEncoding int = 6

const (
	// Unknown quality encoding
	Unknown QualityEncoding = iota
	// Sanger format can encode a Phred quality score from 0 to 93 using
	// ASCII 33 to 126 (although in raw read data the Phred quality score
	// rarely exceeds 60, higher scores are possible in assemblies or read maps).
	Sanger
	// Solexa /Illumina 1.0 format can encode a Solexa/Illumina quality score
	// from -5 to 62 using ASCII 59 to 126 (although in raw read data Solexa
	// scores from -5 to 40 only are expected).
	Solexa
	// Illumina1p3 means Illumina 1.3+.
	// Starting with Illumina 1.3 and before Illumina 1.8, the format
	// encoded a Phred quality score from 0 to 62 using ASCII 64 to 126
	// (although in raw read data Phred scores from 0 to 40 only are expected).
	Illumina1p3
	// Illumina1p5 means Illumina 1.5+.
	// Starting in Illumina 1.5 and before Illumina 1.8, the Phred scores
	//  0 to 2 have a slightly different meaning. The values 0 and 1 are
	// no longer used and the value 2, encoded by ASCII 66 "B", is used
	// also at the end of reads as a Read Segment Quality Control Indicator.
	Illumina1p5
	// Illumina1p8 means Illumina 1.8+.
	// Starting in Illumina 1.8, the quality scores have basically
	// returned to the use of the Sanger format (Phred+33)
	Illumina1p8
)

func (qe QualityEncoding) String() string {
	switch qe {
	case Sanger:
		return "Sanger"
	case Solexa:
		return "Solexa"
	case Illumina1p3:
		return "Illumina-1.3+"
	case Illumina1p5:
		return "Illumina-1.5+"
	case Illumina1p8:
		return "Illumina-1.8+"
	}
	return "Unknown"
}

// QualityRange is the typical quality range
func (qe QualityEncoding) QualityRange() []int {
	switch qe {
	case Sanger:
		return []int{33, 73}
	case Solexa:
		return []int{59, 104}
	case Illumina1p3:
		return []int{64, 104}
	case Illumina1p5:
		return []int{66, 105}
	case Illumina1p8:
		return []int{33, 74}
	}
	return []int{127, 256}
}

// Offset is the ASCII offset
func (qe QualityEncoding) Offset() int {
	switch qe {
	case Sanger:
		return 33
	case Solexa:
		return 64
	case Illumina1p3:
		return 64
	case Illumina1p5:
		return 64
	case Illumina1p8:
		return 33
	}
	return 0
}

// IsSolexa tells whether the encoding is Solexa
func (qe QualityEncoding) IsSolexa() bool {
	switch qe {
	case Solexa:
		return true
	}
	return false
}

// ErrUnknownQualityEncoding is error for Unknown quality encoding type
var ErrUnknownQualityEncoding = errors.New("unkown quality encoding")

// QualityValue returns quality value for given encoding and quality string
func QualityValue(encoding QualityEncoding, quality []byte) ([]int, error) {
	offset := encoding.Offset()
	if offset == 0 {
		return nil, ErrUnknownQualityEncoding
	}

	qv := make([]int, len(quality))
	for i, q := range quality {
		qv[i] = int(q) - offset
	}
	return qv, nil
}

// QualityConvert convert quality from one encoding to another encoding.
// Force means forcely truncate scores > 40 to 40 when converting Illumina-1.8+
// to Sanger.
func QualityConvert(from, to QualityEncoding, quality []byte, force bool) ([]byte, error) {
	if from == to || from == Unknown || to == Unknown {
		return quality, nil
	}
	qv, _ := QualityValue(from, quality)
	isSolexaFrom := from.IsSolexa()
	isSolexaTo := to.IsSolexa()
	offsetTo := to.Offset()
	var err error
	var q2 float64

	qualityNew := make([]byte, len(quality))
	for i, q := range qv {

		if force { // Illumina -> Sanger
			if from == Illumina1p8 && to == Sanger && q > 40 {
				q = 40
			}
		}

		if from == Illumina1p5 && q == 2 { // special case of Illumina 1.5
			q = 0
		}

		if isSolexaFrom == isSolexaTo {
			qualityNew[i] = byte(q + offsetTo)
		} else if isSolexaFrom {
			q2, err = Solexa2Phred(float64(q))
			if err != nil {
				return nil, err
			}
			qualityNew[i] = byte(int(q2) + offsetTo)
		} else {
			q2, err = Phred2Solexa(float64(q))
			if err != nil {
				return nil, err
			}
			qualityNew[i] = byte(int(q2) + offsetTo)
		}
	}
	return qualityNew, nil
}

// NMostCommonThreshold is the threshold of 'B' in
// top N most common quality for guessing Illumina 1.5.
var NMostCommonThreshold = 4

// GuessQualityEncoding returns potential quality encodings.
func GuessQualityEncoding(quality []byte) []QualityEncoding {
	var encodings []QualityEncoding
	min, max := qualRange(quality)
	var encoding QualityEncoding
	var r []int
	var count map[byte]int
	var countSorted byteutil.ByteCountList
	var BEnriched bool
	for i := 1; i < NQualityEncoding; i++ {
		encoding = QualityEncoding(i)
		r = encoding.QualityRange()
		if min >= r[0] && max <= r[1] {
			if encoding == Illumina1p5 {
				if count == nil {
					count = byteutil.CountOfByte(quality)
				}
				if count['@'] > 0 || count['A'] > 0 { // exclude Illumina 1.5
					continue
				} else { //
					countSorted = byteutil.SortCountOfByte(count, true)
					BEnriched = false
					for i := 0; i < len(countSorted) && i < NMostCommonThreshold; i++ {
						if countSorted[i].Key == 'B' {
							BEnriched = true
							break
						}
					}
					if BEnriched {
						return []QualityEncoding{Illumina1p5}
					}
				}
			}
			encodings = append(encodings, encoding)
		}
	}
	return encodings
}

func qualRange(quality []byte) (int, int) {
	var min, max int = 126, 0
	var v int
	for _, q := range quality {
		v = int(q)
		if v > max {
			max = v
		}
		if v < min {
			min = v
		}
	}
	return min, max
}
