/*Package fai implements fasta sequence file index handling, including creating
, reading and random accessing.

Package fai implements fasta sequence file index handling, including creating
, reading and random accessing.

Code of fai data structure were copied and edited from [1].

But I wrote the code of creating and reading fai, and so did test code.

Code of random accessing subsequences were copied from [2], but I extended them.

Reference:

1. https://github.com/biogo/biogo/blob/master/io/seqio/fai/fai.go
2. https://github.com/brentp/faidx/blob/master/faidx.go

Examples:

    import "github.com/shenwei356/bio/seqio/fai"

    file := "seq.fa"
    idx, err := New(file)
    checkErr(err)
    defer idx.Close()

    // single base
    s, err := idx.Base("cel-let-7", 1)
    checkErr(err)

    // subsequence
    seq, err := idx.SubSeq("cel-mir-2", 15, 19)
    checkErr(err)

    // whole sequence
    seq, err := idx.Seq("cel-mir-2")
    checkErr(err)

*/
package fai
