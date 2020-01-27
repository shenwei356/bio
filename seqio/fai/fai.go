package fai

import (
	"bufio"
	"bytes"
	"fmt"
	"os"
	"regexp"
	"strconv"
	"strings"
)

// Record is FASTA index record
type Record struct {
	Name         string
	Length       int
	Start        int64
	BasesPerLine int
	BytesPerLine int
}

// Index is FASTA index
type Index map[string]Record

// Read faidx from .fai file
func Read(fileFai string) (Index, error) {
	fh, err := os.Open(fileFai)
	if err != nil {
		return nil, fmt.Errorf("read faidx: %s", err)
	}
	defer fh.Close()

	index := make(map[string]Record)

	reader := bufio.NewReader(fh)
	var items []string
	var line, name string
	var length int
	var start int64
	var BasesPerLine, bytesPerLine int
	for {
		line, err = reader.ReadString('\n')
		if line != "" {
			line = string(dropCR([]byte(line[0 : len(line)-1])))
			items = strings.Split(line, "\t")
			if len(items) != 5 {
				return nil, fmt.Errorf("invalid fai records: %s", line)
			}
			name = items[0]

			length, err = strconv.Atoi(items[1])
			if err != nil {
				return nil, fmt.Errorf("invalid fai records: %s", line)
			}

			start, err = strconv.ParseInt(items[2], 10, 64)
			if err != nil {
				return nil, fmt.Errorf("invalid fai records: %s", line)
			}

			BasesPerLine, err = strconv.Atoi(items[3])
			if err != nil {
				return nil, fmt.Errorf("invalid fai records: %s", line)
			}

			bytesPerLine, err = strconv.Atoi(items[4])
			if err != nil {
				return nil, fmt.Errorf("invalid fai records: %s", line)
			}

			index[name] = Record{
				Name:         name,
				Length:       length,
				Start:        start,
				BasesPerLine: BasesPerLine,
				BytesPerLine: bytesPerLine,
			}
		}

		if err != nil {
			break
		}
	}

	return index, nil
}

// CreateWithFullHead uses full head instead of just sequence ID
func CreateWithFullHead(fileSeq, fileFai string) (Index, error) {
	return CreateWithIDRegexp(fileSeq, fileFai, `^(.+)$`)
}

// CreateWithIDRegexp uses custom regular expression to get sequence ID
func CreateWithIDRegexp(fileSeq, fileFai string, idRegexp string) (Index, error) {
	if idRegexp != defaultIDRegexp {
		if !reCheckIDregexpStr.MatchString(idRegexp) {
			return nil, fmt.Errorf(`regular expression must contain "(" and ")" to capture matched ID. default: %s`, `^([^\s]+)\s?`)
		}
		var err error
		IDRegexp, err = regexp.Compile(idRegexp)
		if err != nil {
			return nil, fmt.Errorf("fail to Compile idRegexp: %s", err)
		}
		isUsingDefaultIDRegexp = false
	}
	return Create(fileSeq, fileFai)
}

// Create .fai for file
func Create(fileSeq, fileFai string) (Index, error) {
	fh, err := os.Open(fileSeq)
	if err != nil {
		return nil, fmt.Errorf("fail to open seq file: %s", err)
	}
	defer fh.Close()

	outfh, err := os.Create(fileFai)
	if err != nil {
		return nil, fmt.Errorf("fail to write fai file: %s", err)
	}
	defer outfh.Close()

	index := make(map[string]Record)

	reader := bufio.NewReader(fh)
	checkSeqType := true
	seqLen := 0
	var hasSeq bool
	var lastName, thisName []byte
	var id string
	var lastStart, thisStart int64
	var lineWidths, seqWidths []int
	var lastLineWidth, lineWidth, seqWidth int
	var chances int
	for {
		line, err := reader.ReadBytes('\n')
		if err != nil {
			id = string(parseHeadID(lastName))

			// check lineWidths
			lastLineWidth, chances = -2, 2
			for i := len(lineWidths) - 1; i >= 0; i-- {
				if lastLineWidth == -2 {
					lastLineWidth = lineWidths[i]
					continue
				}
				if lineWidths[i] != lastLineWidth {
					chances--
					if chances == 0 {
						return nil, fmt.Errorf("different line length in sequence: %s", id)
					}
				}
				lastLineWidth = lineWidths[i]
			}
			// lineWidth = 0
			if len(lineWidths) > 0 {
				lineWidth = lineWidths[0]
			}
			// seqWidth = 0
			if len(seqWidths) > 0 {
				seqWidth = seqWidths[0]
			}

			if len(line) > 0 && line[len(line)-1] != '\n' {
				fmt.Fprintln(os.Stderr, `[WARNING]: newline character ('\n') not detected at end of file, truncated file?`)
			}

			seqLen += len(line)

			if _, ok := index[id]; ok {
				// return index, fmt.Errorf(`ignoring duplicate sequence "%s" at byte offset %d`, id, lastStart)
				os.Stderr.WriteString(fmt.Sprintf("[fai warning] ignoring duplicate sequence \"%s\" at byte offset %d\n", id, lastStart))
			} else {
				outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\n", id, seqLen, lastStart, seqWidth, lineWidth))
				index[id] = Record{
					Name:         id,
					Length:       seqLen,
					Start:        lastStart,
					BasesPerLine: seqWidth,
					BytesPerLine: lineWidth,
				}
			}

			seqLen = 0

			break
		}

		if checkSeqType {
			if line[0] == '@' {
				os.Remove(fileFai)
				return nil, fmt.Errorf("FASTQ format not supported")
			}
			checkSeqType = false
		}
		if line[0] == '>' {
			hasSeq = true
			thisName = dropCR(line[1 : len(line)-1])

			if lastName != nil {
				id = string(parseHeadID(lastName))

				// check lineWidths
				lastLineWidth, chances = -1, 2
				for i := len(lineWidths) - 1; i >= 0; i-- {
					if lastLineWidth == -1 {
						lastLineWidth = lineWidths[i]
						continue
					}
					if lineWidths[i] != lastLineWidth {
						chances--
						if chances == 0 {
							return nil, fmt.Errorf("different line length in sequence: %s", id)
						}
					}
					lastLineWidth = lineWidths[i]
				}
				// lineWidth = 0
				if len(lineWidths) > 0 {
					lineWidth = lineWidths[0]
				}
				// seqWidth = 0
				if len(seqWidths) > 0 {
					seqWidth = seqWidths[0]
				}

				if _, ok := index[id]; ok {
					// return index, fmt.Errorf(`ignoring duplicate sequence "%s" at byte offset %d`, id, lastStart)
					os.Stderr.WriteString(fmt.Sprintf("[fai warning] ignoring duplicate sequence \"%s\" at byte offset %d\n", id, lastStart))
				} else {
					outfh.WriteString(fmt.Sprintf("%s\t%d\t%d\t%d\t%d\n", id, seqLen, lastStart, seqWidth, lineWidth))
					index[id] = Record{
						Name:         id,
						Length:       seqLen,
						Start:        lastStart,
						BasesPerLine: seqWidth,
						BytesPerLine: lineWidth,
					}
				}

				seqLen = 0
			}
			lineWidths = []int{}
			seqWidths = []int{}
			thisStart += int64(len(line))
			lastStart = thisStart
			lastName = thisName
		} else if hasSeq {
			seqLen += len((dropCR(line[0 : len(line)-1])))
			thisStart += int64(len(line))

			lineWidths = append(lineWidths, len(line))
			seqWidths = append(seqWidths, len(dropCR(line[0:len(line)-1])))
		}
	}

	return index, nil
}

// ------------------------------------------------------------

var reCheckIDregexpStr = regexp.MustCompile(`\(.+\)`)

var defaultIDRegexp = `^(\S+)\s?`

// IDRegexp is regexp for parsing record id
var IDRegexp = regexp.MustCompile(defaultIDRegexp)
var isUsingDefaultIDRegexp = true

func parseHeadID(head []byte) []byte {
	if isUsingDefaultIDRegexp {
		if i := bytes.IndexByte(head, ' '); i > 0 {
			return head[0:i]
		}
		if i := bytes.IndexByte(head, '\t'); i > 0 {
			return head[0:i]
		}
		return head
	}

	found := IDRegexp.FindSubmatch(head)
	if found == nil { // not match
		return head
	}
	return found[1]
}

func dropCR(data []byte) []byte {
	if len(data) > 0 && data[len(data)-1] == '\r' {
		return data[0 : len(data)-1]
	}
	return data
}
