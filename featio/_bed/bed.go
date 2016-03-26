//Package bed is used to read bed features.
// ref: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
// ref: https://github.com/biogo/biogo/blob/master/io/featio/bed/bed.go
package bed

import (
	"errors"
	"fmt"
	"os"
	"reflect"
	"regexp"
	"runtime"
	"strconv"
	"strings"

	"github.com/shenwei356/breader"
)

var (
	// ErrBadBEDType error
	ErrBadBEDType = errors.New("bad BED type")
	// ErrBadBrowserLine error
	ErrBadBrowserLine = errors.New("bad browser line")
)

// Feature is the interface of BED feature
type Feature interface {
	Chr() string
	Start() int
	End() int
	Len() int
	Strand() *string
	String() string
}

// BED3 struct
type BED3 struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
}

// Chr returns chromosome
func (b *BED3) Chr() string { return b.Chrom }

// Start returns start
func (b *BED3) Start() int { return b.ChromStart }

// End returns end. (the site is not included)
func (b *BED3) End() int { return b.ChromEnd }

// Length returns length
func (b *BED3) Length() int { return b.ChromEnd - b.ChromStart }

// Strand returns strand
func (b *BED3) Strand() *string { return nil }

// Name returns name
func (b *BED3) Name() string { return b.String() }
func (b *BED3) String() string {
	return fmt.Sprintf("BED3 %s:[%d,%d)", b.Chrom, b.ChromStart, b.ChromEnd)
}

func parseBED3(line string) (interface{}, bool, error) {
	if len(line) == 0 {
		return nil, false, nil
	}
	items := strings.Split(strings.TrimRight(line, "\n"), "\t")
	if len(items) < 3 {
		return nil, false, ErrBadBEDType
	}
	start, err := strconv.Atoi(items[1])
	if err != nil {
		return nil, false, fmt.Errorf("bad start: %s", items[1])
	}
	end, err := strconv.Atoi(items[2])
	if err != nil {
		return nil, false, fmt.Errorf("bad end: %s", items[2])
	}
	return BED3{items[0], start, end}, true, nil
}

// BED4 struct
type BED4 struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	FeatName   string
}

// Chr returns chromosome
func (b *BED4) Chr() string { return b.Chrom }

// Start returns start
func (b *BED4) Start() int { return b.ChromStart }

// End returns end. (the site is not included)
func (b *BED4) End() int { return b.ChromEnd }

// Length returns length
func (b *BED4) Length() int { return b.ChromEnd - b.ChromStart }

// Strand returns strand
func (b *BED4) Strand() *string { return nil }

// Name returns name
func (b *BED4) Name() string { return b.FeatName }
func (b *BED4) String() string {
	return fmt.Sprintf("BED4 %s:[%d,%d) %s", b.Chrom, b.ChromStart, b.ChromEnd, b.FeatName)
}

func parseBED4(line string) (interface{}, bool, error) {
	if len(line) == 0 {
		return nil, false, nil
	}
	items := strings.Split(strings.TrimRight(line, "\n"), "\t")
	if len(items) < 4 {
		return nil, false, ErrBadBEDType
	}
	start, err := strconv.Atoi(items[1])
	if err != nil {
		return nil, false, fmt.Errorf("bad start: %s", items[1])
	}
	end, err := strconv.Atoi(items[2])
	if err != nil {
		return nil, false, fmt.Errorf("bad end: %s", items[2])
	}
	return BED4{items[0], start, end, items[3]}, true, nil
}

// BED5 struct
type BED5 struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	FeatName   string
	FeatScore  int
}

// Chr returns chromosome
func (b *BED5) Chr() string { return b.Chrom }

// Start returns start
func (b *BED5) Start() int { return b.ChromStart }

// End returns end. (the site is not included)
func (b *BED5) End() int { return b.ChromEnd }

// Length returns length
func (b *BED5) Length() int { return b.ChromEnd - b.ChromStart }

// Strand returns strand
func (b *BED5) Strand() *string { return nil }

// Name returns name
func (b *BED5) Name() string { return b.FeatName }
func (b *BED5) String() string {
	return fmt.Sprintf("BED5 %s:[%d,%d) %s (score: %d)", b.Chrom, b.ChromStart, b.ChromEnd, b.FeatName, b.FeatScore)
}

func parseBED5(line string) (interface{}, bool, error) {
	if len(line) == 0 {
		return nil, false, nil
	}
	items := strings.Split(strings.TrimRight(line, "\n"), "\t")
	if len(items) < 5 {
		return nil, false, ErrBadBEDType
	}
	start, err := strconv.Atoi(items[1])
	if err != nil {
		return nil, false, fmt.Errorf("bad start: %s", items[1])
	}
	end, err := strconv.Atoi(items[2])
	if err != nil {
		return nil, false, fmt.Errorf("bad end: %s", items[2])
	}
	score, err := strconv.Atoi(items[4])
	if err != nil {
		return nil, false, fmt.Errorf("bad score: %s", items[4])
	}
	return BED5{items[0], start, end, items[3], score}, true, nil
}

// BED6 struct
type BED6 struct {
	Chrom      string
	ChromStart int
	ChromEnd   int
	FeatName   string
	FeatScore  int
	FeatStrand *string
}

// Chr returns chromosome
func (b *BED6) Chr() string { return b.Chrom }

// Start returns start
func (b *BED6) Start() int { return b.ChromStart }

// End returns end. (the site is not included)
func (b *BED6) End() int { return b.ChromEnd }

// Length returns length
func (b *BED6) Length() int { return b.ChromEnd - b.ChromStart }

// Strand returns strand
func (b *BED6) Strand() *string { return b.FeatStrand }

// Name returns name
func (b *BED6) Name() string { return b.FeatName }
func (b *BED6) String() string {
	return fmt.Sprintf("BED6 %s:[%d,%d)%s %s (score: %d)", b.Chrom, b.ChromStart, b.ChromEnd, *b.FeatStrand, b.FeatName, b.FeatScore)
}

func parseBED6(line string) (interface{}, bool, error) {
	if len(line) == 0 {
		return nil, false, nil
	}
	items := strings.Split(strings.TrimRight(line, "\n"), "\t")
	if len(items) < 6 {
		return nil, false, ErrBadBEDType
	}
	start, err := strconv.Atoi(items[1])
	if err != nil {
		return nil, false, fmt.Errorf("bad start: %s", items[1])
	}
	end, err := strconv.Atoi(items[2])
	if err != nil {
		return nil, false, fmt.Errorf("bad end: %s", items[2])
	}
	score, err := strconv.Atoi(items[4])
	if err != nil {
		return nil, false, fmt.Errorf("bad score: %s", items[4])
	}
	var strand *string
	if items[5] != "." {
		strand = &items[5]
	}
	return BED6{items[0], start, end, items[3], score, strand}, true, nil
}

// BED12 struct
type BED12 struct {
	Chrom       string
	ChromStart  int
	ChromEnd    int
	FeatName    string
	FeatScore   int
	FeatStrand  *string
	ThickStart  int
	ThickEnd    int
	RGB         string
	BlockCount  int
	BlockSizes  []int
	BlockStarts []int
}

// Chr returns chromosome
func (b *BED12) Chr() string { return b.Chrom }

// Start returns start
func (b *BED12) Start() int { return b.ChromStart }

// End returns end. (the site is not included)
func (b *BED12) End() int { return b.ChromEnd }

// Length returns length
func (b *BED12) Length() int { return b.ChromEnd - b.ChromStart }

// Strand returns strand
func (b *BED12) Strand() *string { return b.FeatStrand }

// Name returns name
func (b *BED12) Name() string { return b.FeatName }
func (b *BED12) String() string {
	return fmt.Sprintf("BED12 %s:[%d,%d)%s %s (score: %d)", b.Chrom, b.ChromStart, b.ChromEnd, *b.FeatStrand, b.FeatName, b.FeatScore)
}

func parseBED12(line string) (interface{}, bool, error) {
	if len(line) == 0 {
		return nil, false, nil
	}
	items := strings.Split(strings.TrimRight(line, "\n"), "\t")
	if len(items) < 12 {
		return nil, false, ErrBadBEDType
	}
	start, err := strconv.Atoi(items[1])
	if err != nil {
		return nil, false, fmt.Errorf("bad start: %s", items[1])
	}
	end, err := strconv.Atoi(items[2])
	if err != nil {
		return nil, false, fmt.Errorf("bad end: %s", items[2])
	}
	score, err := strconv.Atoi(items[4])
	if err != nil {
		return nil, false, fmt.Errorf("bad score: %s", items[4])
	}
	var strand *string
	if items[5] != "." {
		strand = &items[5]
	}
	thickStart, err := strconv.Atoi(items[6])
	if err != nil {
		return nil, false, fmt.Errorf("bad thick start: %s", items[6])
	}
	thickEnd, err := strconv.Atoi(items[7])
	if err != nil {
		return nil, false, fmt.Errorf("bad thick end: %s", items[7])
	}
	blockCount, err := strconv.Atoi(items[9])
	if err != nil {
		return nil, false, fmt.Errorf("bad block count: %s", items[9])
	}
	blockSizes := []int{}
	if blockCount > 0 {
		items2 := strings.Split(items[10], ",")
		for _, size := range items2[0:blockCount] {
			s, err := strconv.Atoi(size)
			if err != nil {
				return nil, false, fmt.Errorf("bad block size: %s", size)
			}
			blockSizes = append(blockSizes, s)
		}
	}
	blockStarts := []int{}
	if blockCount > 0 {
		items2 := strings.Split(items[11], ",")
		for _, start := range items2[0:blockCount] {
			s, err := strconv.Atoi(start)
			if err != nil {
				return nil, false, fmt.Errorf("bad block start: %s", start)
			}
			blockStarts = append(blockStarts, s)
		}
	}

	return BED12{items[0], start, end, items[3], score, strand,
		thickStart, thickEnd, items[8], blockCount, blockSizes, blockStarts}, true, nil
}

type meta struct {
	name    string
	details map[string]string
}

// TrackItemRegexp is regular expression for parsing track items
var TrackItemRegexp = regexp.MustCompile(`(\w+)="?[^#]+?"?`)

// ReadFeatures returns bed data of a file, availabe type values are 3,4,5,6
func ReadFeatures(file string, n int) ([]Feature, map[string]map[string]string, error) {
	if _, err := os.Stat(file); os.IsNotExist(err) {
		return nil, nil, err
	}
	fn := func(line string) (interface{}, bool, error) {
		if line[0] == '#' {
			return nil, false, nil
		}
		if string(line[0:7]) == "browser" {
			items := strings.Split(strings.TrimRight(line, "\n"), " ")
			if len(items) < 3 {
				return nil, false, ErrBadBrowserLine
			}
			details := make(map[string]string)
			details[items[1]] = items[2]
			return meta{"browser", details}, true, nil
		}
		if string(line[0:5]) == "track" {
			details := make(map[string]string)
			found := TrackItemRegexp.FindAllStringSubmatch(line, -1)
			for _, sub := range found {
				details[sub[0]] = sub[1]
			}
			return meta{"track", details}, true, nil
		}
		return nil, false, nil
	}

	reader, err := breader.NewBufferedReader(file, runtime.NumCPU(), 100, fn)
	if err != nil {
		return nil, nil, err
	}
	features := []Feature{}
	for chunk := range reader.Ch {
		if chunk.Err != nil {
			return nil, nil, chunk.Err
		}
		for _, data := range chunk.Data {
			fmt.Println(reflect.TypeOf(data).Kind())
		}
	}
	return features, nil, nil
}
