# fai

Package fai implements fasta sequence file index handling, including creating
, reading and random accessing.

Code of fai data structure were copied and edited from [1].

But I wrote the code of creating and reading fai, and so did test code.

Code of random accessing subsequences are copied and editted from [2]

Reference:

[1]. https://github.com/biogo/biogo/blob/master/io/seqio/fai/fai.go

[2]. https://github.com/brentp/faidx/blob/master/faidx.go

## Usage

    import "github.com/shenwei356/bio/seqio/fai"

    file := "seq.fa"
	idx, err := New(file)
	if err != nil {
		t.Error(err)
	}

	s, err := idx.At("cel-let-7", 1)
	if err != nil {
		t.Error(err)
	}
	if s != 'U' {
		t.Errorf("unmatched sequences: cel-let-7")
	}

	seq, err := idx.Get("cel-mir-2", 15, 20)
	if err != nil {
		t.Error(err)
	}
	if seq != "AAAGC" {
		t.Errorf("unmatched sequences cel-mir-2")
	}


## Documentation

[Documentation on godoc](https://godoc.org/github.com/shenwei356/bio/seqio/fai).
