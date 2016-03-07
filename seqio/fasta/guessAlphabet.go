package fasta

import (
	"errors"

	"github.com/shenwei356/bio/seq"
)

// GuessAlphabet guess the alphabet of the file
func GuessAlphabet(file string) (*seq.Alphabet, error) {
	fastaReader, err := NewFastaReader(seq.Unlimit, file, 1, 1, "")
	if err != nil {
		return seq.Unlimit, err
	}

	for {
		select {
		case chunk := <-fastaReader.Ch:
			if chunk.Err != nil {
				return seq.Unlimit, chunk.Err
			}
			if len(chunk.Data) == 0 {
				return seq.Unlimit, errors.New("no fasta records found in file: " + file)
			}
			firstRecord := chunk.Data[0]
			fastaReader.Cancel()
			return seq.GuessAlphabet(firstRecord.Seq.Seq), nil
		default:
		}
	}
}
