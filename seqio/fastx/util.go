package fastx

import (
	"github.com/shenwei356/bio/seq"
)

// GetSeqNames returns the names of a fasta file
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