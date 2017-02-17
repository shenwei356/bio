package fastx

import (
	"io"

	"github.com/shenwei356/bio/seq"
)

// GetSeqNames returns the names of a fasta/q file
func GetSeqNames(file string) ([]string, error) {
	names := []string{}
	seq.ValidateSeq = false
	reader, err := NewDefaultReader(file)
	if err != nil {
		return nil, nil
	}
	for {
		record, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return nil, err
		}
		names = append(names, string(record.Name))
	}
	return names, nil
}

// GetSeqNumber returns the sequences number of FASTA/Q files
func GetSeqNumber(file string) (int, error) {
	n := 0
	seq.ValidateSeq = false
	reader, err := NewDefaultReader(file)
	if err != nil {
		return 0, nil
	}
	for {
		_, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			return 0, err
		}
		n++
	}
	return n, nil
}

// GetSeqs return fastx records of a file.
// when alphabet is nil or seq.Unlimit, it will automaticlly detect the alphabet.
// when idRegexp is "", default idRegexp ( ^([^\s]+)\s? ) will be used.
func GetSeqs(file string, alphabet *seq.Alphabet, bufferSize int, chunkSize int, idRegexp string) ([]*Record, error) {
	records := []*Record{}

	reader, err := NewReader(alphabet, file, idRegexp)
	if err != nil {
		return records, err
	}
	for chunk := range reader.ChunkChan(bufferSize, chunkSize) {
		if err != nil {
			return records, err
		}

		for _, record := range chunk.Data {
			records = append(records, record)
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
	reader, err := NewDefaultReader(file)
	if err != nil {
		return nil, false, err
	}
	for {
		_, err := reader.Read()
		if err != nil {
			if err == io.EOF {
				return reader.Alphabet(), false, io.EOF
			}
			return nil, false, err
		}
		return reader.Alphabet(), reader.IsFastq, nil
	}
}
