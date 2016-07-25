# fai

[![GoDoc](https://godoc.org/github.com/shenwei356/bio?status.svg)](https://godoc.org/github.com/shenwei356/bio/seqio/fai)


Package fai implements fasta sequence file index handling, including creating
, reading and random accessing.

Code of fai data structure were copied and edited from [1].

But I wrote the code of creating and reading fai, and so did test code.

Code of random accessing subsequences were copied from [2], but I extended them a lot.

Reference:

[1]. https://github.com/biogo/biogo/blob/master/io/seqio/fai/fai.go

[2]. https://github.com/brentp/faidx/blob/master/faidx.go

## General Usage

    import "github.com/shenwei356/bio/seqio/fai"

    file := "seq.fa"
    faidx, err := fai.New(file)
    checkErr(err)
    defer func() {
        checkErr(faidx.Close())
    }()

    // whole sequence
    seq, err := faidx.Seq("cel-mir-2")
    checkErr(err)

    // single base
    s, err := faidx.Base("cel-let-7", 1)
    checkErr(err)

    // subsequence. start and end are all 1-based
    seq, err := faidx.SubSeq("cel-mir-2", 15, 19)
    checkErr(err)


## Extended SubSeq


For extended SubSeq, negative position is allowed.


This is my custom locating strategy. Start and end are all 1-based.
To better understand the locating strategy, see examples below:


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

Examples:

    // last 12 bases
    seq, err := faidx.SubSeq("cel-mir-2", -12, -1)
    checkErr(err)

## Advanced Usage

Function `fai.New(file string)` is a wraper to simplify the process of
creating and reading FASTA index . Let's see what's happend inside:

    func New(file string) (*Faidx, error) {
            fileFai := file + ".fai"
            var index Index
            if _, err := os.Stat(fileFai); os.IsNotExist(err) {
                    index, err = Create(file)
                    if err != nil {
                            return nil, err
                    }
            } else {
                    index, err = Read(fileFai)
                    if err != nil {
                            return nil, err
                    }
            }

            return NewWithIndex(file, index)
    }

By default, sequence ID is used as key in FASTA index file.
Inside the package, a regular expression is used to get sequence ID from
full head. The default value is `^([^\s]+)\s?`, i.e. getting
first non-space characters of head.
So you can just use `fai.Create(file string)` to create .fai file.

If you want to use full head instead of sequence ID (first non-space characters of head),
you could use `fai.CreateWithIDRegexp(file string, idRegexp string)` to create faidx.
Here, the `idRegexp` should be `^(.+)$`. For convenience, you can use another function
`CreateWithFullHead`.


## More Advanced Usages

Note that, ***by default, whole file is mapped into shared memory***,
which is OK for small files (smaller than your RAM).
For very big files, you should disable that.
Instead, file seeking is used.

    // change the global variable
    fai.MapWholeFile = false

    // then do other things


## Documentation

[Documentation on godoc](https://godoc.org/github.com/shenwei356/bio/seqio/fai).
