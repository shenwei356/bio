// Package gtf is used to read gtf features.
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

var strandPositive = "+"
var strandNegative = "-"
var strandNotspecified = "."

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

		// if start > end {
		// 	return nil, false, fmt.Errorf("%s: start (%d) must be < end (%d)", items[0], start, end)
		// }

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
			if s == "+" {
				strand = &strandPositive
			} else if s == "-" {
				strand = &strandNegative
			} else {
				return nil, false, fmt.Errorf("%s: illegal strand: %s", items[0], s)
			}
			strand = &s

			if start > end {
				if s == "+" {
					return nil, false, fmt.Errorf(`%s: start (%d) should be < end (%d) when the strand is "+"`, items[0], start, end)
				} else {
					strand = &strandNegative
					tmp := start
					start = end
					end = tmp
				}
			}
		} else {
			strand = &strandNotspecified
		}

		var frame *int
		if items[7] != "." {
			f, err := strconv.Atoi(items[7])
			if err != nil {
				return nil, false, fmt.Errorf("%s: bad frame: %s", items[0], items[7])
			}
			if !(f == 0 || f == 1 || f == 2) {
				return nil, false, fmt.Errorf("%s: illegal frame: %d", items[0], f)
			}
			frame = &f
		}

		feature := Feature{items[0], items[1], items[2], start, end, score, strand, frame, nil}

		attributes := items[8]
		feature.Attributes = []Attribute{}
		var attributeStart int // Start of the current attribute in column 9.
		var quoted bool        // Semicolons inside quoted values are not separators.
		for i := 0; i <= len(attributes); i++ {
			if i < len(attributes) {
				if quoted && attributes[i] == '\\' && i+1 < len(attributes) {
					i++
					continue
				}
				if attributes[i] == '"' {
					quoted = !quoted
				}
				if attributes[i] != ';' || quoted {
					continue
				}
			}
			// Also parse the final attribute when there is no trailing semicolon.
			tagValue := strings.TrimSpace(attributes[attributeStart:i])
			attributeStart = i + 1
			separator := strings.IndexAny(tagValue, " \t")
			if separator < 0 {
				continue
			}
			tag := tagValue[:separator]
			if _, ok := attrsMap[tag]; !ok {
				continue
			}
			value := strings.TrimSpace(tagValue[separator:])
			if value == "" {
				continue
			}
			if value[0] == '"' {
				if len(value) < 2 || value[len(value)-1] != '"' {
					continue
				}
				value = value[1 : len(value)-1]
			}
			feature.Attributes = append(feature.Attributes, Attribute{tag, value})
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
