package fastx

import (
	"runtime"
	"strings"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/breader"
)

// Threads number
var Threads = runtime.NumCPU()

// GetSeqNames returns the names of a fasta/q file
func GetSeqNames(file string) ([]string, error) {
	names := []string{}
	seq.ValidateSeq = false
	fastxReader, err := NewReader(nil, file, 1, 1, "")
	if err != nil {
		return nil, nil
	}
	for chunk := range fastxReader.Ch {
		if chunk.Err != nil {
			return nil, chunk.Err
		}

		for _, record := range chunk.Data {
			names = append(names, string(record.Name))
		}
	}
	return names, nil
}

// GetSeqNumber returns the sequences number of FASTA/Q files
func GetSeqNumber(file string) (int, error) {
	n := 0
	seq.ValidateSeq = false
	fastxReader, err := NewReader(nil, file, 1, 1, "")
	if err != nil {
		return 0, nil
	}
	for chunk := range fastxReader.Ch {
		if chunk.Err != nil {
			return n, chunk.Err
		}
		n += len(chunk.Data)
	}
	return n, nil
}

// EstimateSeqNumber estimates sequences number of FASTA/Q files.
// It may over count for FASTQ file.
func EstimateSeqNumber(file string) (int, error) {
	fn := func(line string) (interface{}, bool, error) {
		line = strings.TrimRight(line, "\r\n")
		if len(line) == 0 {
			return 0, false, nil
		}

		if line[0] == '>' || line[0] == '@' {
			return 1, true, nil
		}
		return 0, false, nil
	}
	reader, err := breader.NewBufferedReader(file, Threads, 100, fn)
	if err != nil {
		return 0, err
	}

	n := 0
	for chunk := range reader.Ch {
		if chunk.Err != nil {
			return n, err
		}
		n += len(chunk.Data)
	}
	return n, nil
}

// GetSeqs return fastx records of a file.
// when alphabet is nil or seq.Unlimit, it will automaticlly detect the alphabet.
// when idRegexp is "", default idRegexp ( ^([^\s]+)\s? ) will be used.
func GetSeqs(file string, alphabet *seq.Alphabet, bufferSize int, chunkSize int, idRegexp string) ([]*Record, error) {
	records := []*Record{}
	if alphabet == nil || alphabet == seq.Unlimit {
		alphabet = nil
	}
	fastxReader, err := NewReader(alphabet, file, bufferSize, chunkSize, idRegexp)
	if err != nil {
		return records, err
	}
	for chunk := range fastxReader.Ch {
		if err != nil {
			return records, err
		}

		for _, record := range chunk.Data {
			records = append(records, record.Clone())
		}
	}
	return records, nil
}

// GetSeqsMap returns all seqs as a map for fasta file
func GetSeqsMap(file string, alphabet *seq.Alphabet, bufferSize int, chunkSize int, idRegexp string) (map[string]*Record, error) {
	m := make(map[string]*Record)
	records, err := GetSeqs(file, alphabet, bufferSize, chunkSize, idRegexp)
	if err != nil {
		return m, err
	}
	for _, record := range records {
		m[string(record.Name)] = record
	}
	return m, nil
}

// GuessAlphabet guess the alphabet of the file by the first maxLen bases
func GuessAlphabet(file string) (*seq.Alphabet, bool, error) {
	var isFastq bool
	var alphabet *seq.Alphabet
	fastxReader, err := NewReader(nil, file, 1, 1, "")
	if err != nil {
		return nil, false, err
	}
LOOP:
	for {
		select {
		case chunk := <-fastxReader.Ch:
			if chunk.Err != nil {
				return nil, false, chunk.Err
			}

			isFastq = fastxReader.IsFastq
			alphabet = fastxReader.Alphabet()

			fastxReader.Cancel()
			break LOOP
		default:
		}
	}
	return alphabet, isFastq, nil
}
