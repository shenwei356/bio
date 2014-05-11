seqio
=====

This package provides operations of reading and writing sequence file, mainly FASTA file right now.

Documentation
-------------
[Documentation on gowalker](http://gowalker.org/github.com/shenwei356/bio/seqio).

Example
-------

    import (
        "fmt"

        "github.com/shenwei356/bio/seq"
        "github.com/shenwei356/bio/seqio"
    )

    func main() {
       // Sequence type must be specified.
       fasta, err := seqio.NewFastaReader(seq.RNAredundant, "hairpin.fa")
    
       // check err, usually caused by a wrong file path.
       if err != nil {
           fmt.Println(err)
           return
       }

       // read and check if more record existed
       for fasta.HasNext() {
            record, err := fasta.NextSeq()
            // check err, the record may contain invalid sequence !!!
            if err != nil {
                fmt.Println(err)
                continue
            }
            // format output
            fmt.Printf(">%s\n%s", record.Id, record.FormatSeq(70))

            // reverse complement
            fmt.Printf("\nrevcom\n%s", seq.FormatSeq(record.Seq.Revcom(), 70))

            // length
            fmt.Printf("Seq length: %d\n", record.Seq.Len)

            // base content
            fmt.Printf("GC content: %.2f\n", record.Seq.BaseContent([]byte("gc")))
            fmt.Println()

            records = append(records, record)
        }

        // write to file
        writer := seqio.NewFastaWriter("tmp.fasta", 70)
        writer.Write(records)
    }
