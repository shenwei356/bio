package fai

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"regexp"
	"strconv"
	"strings"
)

const faiReaderSize = 64 << 10

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
	idRegexpCompiled, err := regexp.Compile(idRegexp)
	if err != nil {
		return nil, fmt.Errorf("fail to compile idRegexp %q: %w", idRegexp, err)
	}
	if idRegexpCompiled.NumSubexp() == 0 {
		return nil, fmt.Errorf("regular expression must contain a capturing group for the ID; default: %s", defaultIDRegexp)
	}
	return create(fileSeq, fileFai, idRegexpCompiled, idRegexp == defaultIDRegexp)
}

// Create .fai for file. Empty files and leading LF-only blank lines are allowed.
func Create(fileSeq, fileFai string) (Index, error) {
	idRegexp := IDRegexp
	if idRegexp == nil {
		idRegexp = regexp.MustCompile(defaultIDRegexp)
	}
	return create(fileSeq, fileFai, idRegexp, idRegexp.String() == defaultIDRegexp)
}

// create builds an index while preserving the original FASTA byte offsets.
func create(fileSeq, fileFai string, idRegexp *regexp.Regexp, isUsingDefaultIDRegexp bool) (Index, error) {
	fh, err := os.Open(fileSeq)
	if err != nil {
		return nil, fmt.Errorf("fail to open seq file: %s", err)
	}
	defer fh.Close()

	outfh, err := os.Create(fileFai)
	if err != nil {
		return nil, fmt.Errorf("fail to write fai file: %s", err)
	}
	// Remove incomplete indexes so a failed creation cannot succeed on retry.
	complete := false
	defer func() {
		outfh.Close()
		if !complete {
			os.Remove(fileFai)
		}
	}()

	index := make(map[string]Record)

	reader := bufio.NewReaderSize(fh, faiReaderSize)
	checkSeqType := true
	var hasSeq bool
	var state faiSequenceState
	var offset int64
	var scratch []byte
	var tail []byte // final sequence line without a newline
	for {
		var line []byte
		line, err = readFAILine(reader, &scratch)
		if err != nil && err != io.EOF {
			return nil, fmt.Errorf("read fasta file: %w", err)
		}
		// Match fastx/xopen's handling of a UTF-8 BOM, retaining byte offsets.
		if offset == 0 && bytes.HasPrefix(line, []byte("\xef\xbb\xbf")) {
			line = line[3:]
			offset = 3
		}
		if len(line) == 0 {
			break
		}

		if checkSeqType {
			// Only LF-only blank lines may precede the first FASTA header.
			if len(line) == 1 && line[0] == '\n' {
				offset++
				continue
			}
			if line[0] == '@' {
				return nil, fmt.Errorf("FASTQ format not supported")
			}
			checkSeqType = false
		}

		if line[0] == '>' {
			if hasSeq {
				if err := state.finish(index, outfh, idRegexp, isUsingDefaultIDRegexp, nil, false); err != nil {
					return nil, err
				}
			}
			hasSeq = true
			head := line[1:]
			if err != io.EOF {
				head = head[:len(head)-1] // strip the terminating LF
			}
			state.reset(dropCR(head), offset+int64(len(line)))
		} else if hasSeq {
			if err == io.EOF {
				tail = dropCR(line)
				break
			}
			state.addLine(line)
		} else {
			return nil, fmt.Errorf("invalid fasta file: %s", fileSeq)
		}
		offset += int64(len(line))
		if err == io.EOF {
			break
		}
	}

	if hasSeq {
		if err := state.finish(index, outfh, idRegexp, isUsingDefaultIDRegexp, tail, true); err != nil {
			return nil, err
		}
	}
	complete = true
	return index, nil
}

type faiLineWidths struct {
	last        int
	transitions int
	started     bool
	invalid     bool
}

func (w *faiLineWidths) add(width int) {
	if !w.started {
		w.last = width
		w.started = true
		return
	}
	if width == w.last {
		return
	}
	// Preserve the previous rule: line widths may change once, and only to a
	// shorter final width.
	w.transitions++
	if w.transitions > 1 || width > w.last {
		w.invalid = true
	}
	w.last = width
}

type faiSequenceState struct {
	name          []byte
	start         int64
	length        int
	firstSeqWidth int
	firstWidth    int
	hasLine       bool

	widths              faiLineWidths
	pendingEmptyWidths  faiLineWidths
	hasPendingEmptyLine bool
}

func (s *faiSequenceState) reset(name []byte, start int64) {
	s.name = append(s.name[:0], name...)
	s.start = start
	s.length = 0
	s.hasLine = false
	s.widths = faiLineWidths{}
	s.pendingEmptyWidths = faiLineWidths{}
	s.hasPendingEmptyLine = false
}

func (s *faiSequenceState) addLine(line []byte) {
	lineWidth := len(line)
	seqWidth := len(dropCR(line[:lineWidth-1]))
	s.length += seqWidth
	if !s.hasLine {
		s.firstSeqWidth = seqWidth
		s.firstWidth = lineWidth
		s.hasLine = true
	}

	if seqWidth == 0 {
		// Trailing empty lines are ignored by the existing FAI validation.
		// Keep their tentative state until another non-empty line appears.
		if !s.hasPendingEmptyLine {
			s.pendingEmptyWidths = s.widths
			s.hasPendingEmptyLine = true
		}
		s.pendingEmptyWidths.add(lineWidth)
		return
	}

	if s.hasPendingEmptyLine {
		s.widths = s.pendingEmptyWidths
		s.hasPendingEmptyLine = false
	}
	s.widths.add(lineWidth)
}

// finish validates and writes a sequence, including an optional unterminated tail.
func (s *faiSequenceState) finish(index Index, out io.Writer, idRegexp *regexp.Regexp, isUsingDefaultIDRegexp bool, tail []byte, warnNoFinalNewline bool) error {
	id := string(parseHeadID(idRegexp, isUsingDefaultIDRegexp, s.name))
	if strings.Contains(id, "\t") {
		id = reTabs.ReplaceAllString(id, " ")
	}
	if s.widths.invalid {
		return fmt.Errorf("different line length in sequence: %s. Please format the file with 'seqkit seq'", id)
	}
	if warnNoFinalNewline && len(tail) > 0 {
		fmt.Fprintln(os.Stderr, `[WARNING]: newline character ('\n') not detected at end of file, truncated file?`)
	}

	length := s.length + len(tail)
	if !s.hasLine && len(tail) > 0 {
		// A single unterminated sequence line still needs widths for random access.
		s.firstSeqWidth = len(tail)
		s.firstWidth = len(tail)
	}
	if _, ok := index[id]; ok {
		os.Stderr.WriteString(fmt.Sprintf("[fai warning] ignoring duplicate sequence \"%s\" at byte offset %d\n", id, s.start))
		return nil
	}
	if _, err := fmt.Fprintf(out, "%s\t%d\t%d\t%d\t%d\n", id, length, s.start, s.firstSeqWidth, s.firstWidth); err != nil {
		return err
	}
	index[id] = Record{
		Name:         id,
		Length:       length,
		Start:        s.start,
		BasesPerLine: s.firstSeqWidth,
		BytesPerLine: s.firstWidth,
	}
	return nil
}

func readFAILine(reader *bufio.Reader, scratch *[]byte) ([]byte, error) {
	line, err := reader.ReadSlice('\n')
	if err != bufio.ErrBufferFull {
		return line, err
	}

	*scratch = append((*scratch)[:0], line...)
	for err == bufio.ErrBufferFull {
		line, err = reader.ReadSlice('\n')
		*scratch = append(*scratch, line...)
	}
	return *scratch, err
}

// ------------------------------------------------------------

var defaultIDRegexp = `^(\S+)\s?`

var reTabs = regexp.MustCompile(`\t+`)

// IDRegexp is regexp for parsing record id
var IDRegexp = regexp.MustCompile(defaultIDRegexp)

func parseHeadID(idRegexp *regexp.Regexp, isUsingDefaultIDRegexp bool, head []byte) []byte {
	if isUsingDefaultIDRegexp {
		if i := bytes.IndexByte(head, ' '); i > 0 {
			return head[0:i]
		}
		if i := bytes.IndexByte(head, '\t'); i > 0 {
			return head[0:i]
		}
		return head
	}

	found := idRegexp.FindSubmatch(head)
	if len(found) < 2 { // not match or no capturing group
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
