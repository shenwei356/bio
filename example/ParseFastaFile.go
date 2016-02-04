package main

import (
	"fmt"
	"os"

	"github.com/shenwei356/bio/seq"
	"github.com/shenwei356/bio/seqio"
)

func main() {
	if len(os.Args) == 1 {
		fmt.Fprintln(os.Stderr, "Please feed me one fasta file!")
		os.Exit(0)
	}
	fasta, err := seqio.NewFastaReader(seq.Unlimit, os.Args[1])
	//fasta, err := seqio.NewFastaReader(seq.Unlimit, "test.fa")
	if err != nil {
		fmt.Println(err)
		return
	}

	var records []*seqio.FastaRecord

	// You can use Iterator or HasNext() - NextSeq() pair.
	// for record := range fasta.Iterator(10) {
	//	fmt.Printf(">%s\n%s", record.Id, record.FormatSeq(70))
	// }

	// for fasta.HasNext() {
	// 	record, err := fasta.NextSeq()
	// 	if err != nil { // invalid sequence
	// 		fmt.Println(err)
	// 		continue
	// 	}

	for record := range fasta.Iterator(100) {
		s := record.Seq

		// fmt.Printf(">%s\n%s\n", record.Id, record.Seq.Seq)

		// format output
		fmt.Printf(">%s\n%s", record.ID, record.FormatSeq(70))

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
