# fai

Package fai implements fasta sequence file index handling, including creating
, reading and random accessing.

Code of fai data structure were copied and edited from [1].

But I wrote the code of creating and reading fai, and so did test code.

Code of random accessing subsequences were copied from [2], but I extended them.

Reference:

[1]. https://github.com/biogo/biogo/blob/master/io/seqio/fai/fai.go

[2]. https://github.com/brentp/faidx/blob/master/faidx.go

## Usage

    import "github.com/shenwei356/bio/seqio/fai"

    file := "seq.fa"
    idx, err := New(file)
    checkErr(err)
    defer idx.Close()

    // single base
    s, err := idx.Base("cel-let-7", 1)
    checkErr(err)

    // subsequence. start and end are all 1-based
    seq, err := idx.SubSeq("cel-mir-2", 15, 19)
    checkErr(err)

    // whole sequence
    seq, err := idx.Seq("cel-mir-2")
    checkErr(err)

Extended SubSeq

start and end are all 1-based.


     1-based index    1 2 3 4 5 6 7 8 9 10
    negative index    0-9-8-7-6-5-4-3-2-1
               seq    A C G T N a c g t n
               1:1    A
               2:4      C G T
             -4:-2                c g t
             -4:-1                c g t n
             -1:-1                      n
              2:-2      C G T N a c g t
              1:-1    A C G T N a c g t n



## Documentation

[Documentation on godoc](https://godoc.org/github.com/shenwei356/bio/seqio/fai).
