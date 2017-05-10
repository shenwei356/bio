//Package gtf is used to read gtf features.
// ref: http://mblab.wustl.edu/GTF22.html
package gtf

import (
	"fmt"
	"os"
	"runtime"
	"strconv"
	"strings"

	"github.com/shenwei356/breader"
)

// Version is the GTF version
const Version = 2.2

// Feature is the gff feature struct
type Feature struct {
	SeqName    string
	Source     string
	Feature    string
	Start      int
	End        int
	Score      *float64
	Strand     *string
	Frame      *int
	Attributes []Attribute
}

// Attribute is the attribute
type Attribute struct {
	Tag, Value string
}

// Threads for bread.NewBufferedReader()
var Threads = runtime.NumCPU()

// ReadFeatures returns gtf features of a file
func ReadFeatures(file string) ([]Feature, error) {
	return ReadFilteredFeatures(file, []string{}, []string{}, []string{})
}

// ReadFilteredFeatures returns gtf features of specific chrs in a file
func ReadFilteredFeatures(file string, chrs []string, feats []string, attrs []string) ([]Feature, error) {
	if _, err := os.Stat(file); os.IsNotExist(err) {
		return nil, err
	}
	chrsMap := make(map[string]struct{}, len(chrs))
	for _, chr := range chrs {
		chrsMap[strings.ToLower(chr)] = struct{}{}
	}

	featsMap := make(map[string]struct{}, len(feats))
	for _, f := range feats {
		featsMap[strings.ToLower(f)] = struct{}{}
	}

	attrsMap := make(map[string]struct{}, len(attrs))
	for _, f := range attrs {
		attrsMap[strings.ToLower(f)] = struct{}{}
	}

	fn := func(line string) (interface{}, bool, error) {
		if len(line) == 0 || line[0] == '#' {
			return nil, false, nil
		}
		line = strings.TrimRight(line, "\r\n")
		items := strings.Split(line, "\t")

		if len(items) != 9 {
			return nil, false, nil
		}

		if len(chrs) > 0 { // selected chrs
			if _, ok := chrsMap[strings.ToLower(items[0])]; !ok {
				return nil, false, nil
			}
		}

		if len(feats) > 0 { // selected features
			if _, ok := featsMap[strings.ToLower(items[2])]; !ok {
				return nil, false, nil
			}
		}
		var err error

		start, err := strconv.Atoi(items[3])
		if err != nil {
			return nil, false, fmt.Errorf("%s: bad start: %s", items[0], items[3])
		}

		end, err := strconv.Atoi(items[4])
		if err != nil {
			return nil, false, fmt.Errorf("%s: bad end: %s", items[0], items[4])
		}

		if start > end {
			return nil, false, fmt.Errorf("%s: start (%d) must be < end (%d)", items[0], start, end)
		}

		var score *float64
		if items[5] != "." {
			s, err := strconv.ParseFloat(items[5], 64)
			if err != nil {
				return nil, false, fmt.Errorf("%s: bad score: %s", items[0], items[5])
			}
			score = &s
		}

		var strand *string
		if items[6] != "." {
			s := items[6]
			if !(s == "+" || s == "-") {
				return nil, false, fmt.Errorf("%s: illigal strand: %s", items[0], s)
			}
			strand = &s
		}

		var frame *int
		if items[7] != "." {
			f, err := strconv.Atoi(items[7])
			if err != nil {
				return nil, false, fmt.Errorf("%s: bad frame: %s", items[0], items[7])
			}
			if !(f == 0 || f == 1 || f == 2) {
				return nil, false, fmt.Errorf("%s: illigal frame: %d", items[0], f)
			}
			frame = &f
		}

		feature := Feature{items[0], items[1], items[2], start, end, score, strand, frame, nil}

		tagValues := strings.Split(items[8], "; ")
		if len(tagValues) > 0 {
			var ok bool
			feature.Attributes = []Attribute{}
			for _, tagValue := range tagValues[0 : len(tagValues)-1] {
				items2 := strings.SplitN(tagValue, " ", 2)
				tag := items2[0]
				if _, ok = attrsMap[tag]; !ok {
					continue
				}
				value := items2[1]
				// if value[len(value)-1] == ';' {
				// 	value = value[0 : len(value)-1]
				// }
				if len(value) > 2 {
					value = value[1 : len(value)-1]
				} else {
					value = ""
				}
				feature.Attributes = append(feature.Attributes, Attribute{tag, value})
			}
		}
		return feature, true, nil
	}
	reader, err := breader.NewBufferedReader(file, Threads, 100, fn)
	if err != nil {
		return nil, err
	}
	features := []Feature{}
	for chunk := range reader.Ch {
		if chunk.Err != nil {
			return nil, chunk.Err
		}
		for _, data := range chunk.Data {
			features = append(features, data.(Feature))
		}
	}
	return features, nil
}
