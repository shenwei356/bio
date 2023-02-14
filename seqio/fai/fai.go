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

	scanner := bufio.NewScanner(fh)
	items := make([]string, 5)
	var line, name string
	var length int
	var start int64
	var BasesPerLine, bytesPerLine int
	for scanner.Scan() {
		line = scanner.Text()
		if line != "" {
			line = dropCRStr(line)
			stringSplitNByByte(line, '\t', 5, &items)
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

		if err := scanner.Err(); err != nil {
			return nil, err
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
	var line, lineDropCR []byte
	var seenSeqs bool
	for {
		line, err = reader.ReadBytes('\n')
		if err != nil { // end of file
			id = string(parseHeadID(lastName))

			// check lineWidths
			lastLineWidth, chances = -2, 2
			seenSeqs = false
			for i := len(lineWidths) - 1; i >= 0; i-- {
				if !seenSeqs && seqWidths[i] == 0 { // skip empty lines in the end
					continue
				}
				seenSeqs = true

				if lastLineWidth == -2 {
					lastLineWidth = lineWidths[i]
					continue
				}
				if lineWidths[i] != lastLineWidth {
					chances--
					if chances == 0 || lineWidths[i] < lastLineWidth {
						return nil, fmt.Errorf("different line length in sequence: %s. Please format the file with 'seqkit seq'", id)
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
				fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\n", id, seqLen, lastStart, seqWidth, lineWidth)
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

			if lastName != nil { // not the first record
				id = string(parseHeadID(lastName))

				// check lineWidths
				lastLineWidth, chances = -1, 2
				seenSeqs = false
				for i := len(lineWidths) - 1; i >= 0; i-- {
					if !seenSeqs && seqWidths[i] == 0 { // skip empty lines in the end
						continue
					}
					seenSeqs = true

					if lastLineWidth == -1 {
						lastLineWidth = lineWidths[i]
						continue
					}
					if lineWidths[i] != lastLineWidth {
						chances--
						if chances == 0 || lineWidths[i] < lastLineWidth {
							return nil, fmt.Errorf("different line length in sequence: %s. Please format the file with 'seqkit seq'", id)
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
					fmt.Fprintf(outfh, "%s\t%d\t%d\t%d\t%d\n", id, seqLen, lastStart, seqWidth, lineWidth)
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
			lineDropCR = dropCR(line[0 : len(line)-1])
			seqLen += len(lineDropCR)
			thisStart += int64(len(line))

			lineWidths = append(lineWidths, len(line))
			seqWidths = append(seqWidths, len(lineDropCR))
		} else {
			return nil, fmt.Errorf("invalid fasta file: %s", fileSeq)
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

func dropCRStr(data string) string {
	if len(data) > 0 && data[len(data)-1] == '\r' {
		return data[0 : len(data)-1]
	}
	return data
}

func stringSplitNByByte(s string, sep byte, n int, a *[]string) {
	if a == nil {
		tmp := make([]string, n)
		a = &tmp
	}

	n--
	i := 0
	for i < n {
		m := strings.IndexByte(s, sep)
		if m < 0 {
			break
		}
		(*a)[i] = s[:m]
		s = s[m+1:]
		i++
	}
	(*a)[i] = s

	(*a) = (*a)[:i+1]
}
