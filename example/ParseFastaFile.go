package main

import (
	"fmt"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio"
)

func main() {
	fasta, err := seqio.NewFastaReader(seq.RNAredundant, "hairpin.fa")
	if err != nil {
		fmt.Println(err)
		return
	}

	records := make([]*seqio.FastaRecord, 0)

	for fasta.HasNext() {
		record, err := fasta.NextSeq()
		if err != nil { // invalid sequence
			fmt.Println(err)
			continue
		}

		s := record.Seq

		// format output
		fmt.Printf(">%s\n%s", record.Id, record.FormatSeq(70))

		// reverse complement
		fmt.Printf("\nrevcom\n%s", seq.FormatSeq(s.Revcom().Seq, 70))

		// length
		fmt.Printf("Seq length: %d\n", s.Length())

		// base content
		fmt.Printf("GC content: %.2f\n", s.BaseContent("gc"))
		fmt.Println()

		records = append(records, record)
	}

	// write to file
	writer := seqio.NewFastaWriter("tmp.fasta", 70)
	writer.Write(records)
}
