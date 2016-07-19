bio
===
[![GoDoc](https://godoc.org/github.com/shenwei356/bio?status.svg)](https://godoc.org/github.com/shenwei356/bio)

A lightweight and high-performance
(see [fakit](https://github.com/shenwei356/fakit#benchmark) benchmark of FASTA parsing)
bioinformatics package.

Compare to kseq.h (klib)
-----------------------

***This package has high performance close to the famous C lib
[`kseq.h`](https://github.com/attractivechaos/klib/blob/master/kseq.h).***

To test the performance, three datasets are used:

- dataset_A, bacteria genomes, 2.7G
- dataset_B, human genome,     2.9G
- dataset_C, Illumina reads,   2.2G

summary by [`fakit`](https://github.com/shenwei356/fakit):

    file           seq_format   seq_type   num_seqs   min_len        avg_len       max_len
    dataset_A.fa   FASTA        DNA          67,748        56       41,442.5     5,976,145
    dataset_B.fa   FASTA        DNA             194       970   15,978,096.5   248,956,422
    dataset_C.fq   FASTQ        DNA       9,186,045       100            100           100

[`seqtk`](https://github.com/lh3/seqtk/)
(Version [1.1-r92-dirty](https://github.com/lh3/seqtk/tree/fb85aad4ce1fc7b3d4543623418a1ae88fe1cea6),
using `kseq.h`)
and [`fakit`](https://github.com/shenwei356/fakit)
(Version [v0.2.8](https://github.com/shenwei356/fakit/releases/tag/v0.2.8),
using this package) were used to test.
**Note** that `seqtk` does not support wrapped (fixed line width) ouputing, so `fakit` uses
`-w 0` to disable outputing wrapping.
Script [`memusg`](https://github.com/shenwei356/memusg) is used to assess running time
and peak memory usage.

Commands:

    for f in dataset_*.f{a,q}; do
        echo $f;
        cat $f > t; /bin/rm t; # warm up

        echo seqtk
        memusg -t seqtk seq $f      > r.seqtk.fx
        echo fakit
        memusg -t fakit seq $f -w 0 > r.fakit.fx

        # md5sum r.seqtk.fx r.fakit.fx
        /bin/rm r.seqtk.fx r.fakit.fx
    done

Results:

dataset     |software|second|memory  |realative.speed
:-----------|:-------|-----:|-------:|--------------:
dataset_A.fa|seqtk   |9.377 |7.34MB  |*1.00*
dataset_A.fa|fakit   |8.568 |38.07MB |***1.09***
dataset_B.fa|seqtk   |10.983|239.03MB|***1.00***
dataset_B.fa|fakit   |11.536|768.63MB|*0.95*
dataset_C.fa|seqtk   |9.572 |956.0KB |***1.00***
dataset_C.fa|fakit   |15.580|13.11MB |*0.61*

For the memory usage, the real-time memory is actually about 300 MB.
However, it's higher due to the limitation of garbage collection mechanism in
 Go programming language, and it may be solved in the future.

Install
-------
This package is "go-gettable", just:

    go get -u github.com/shenwei356/bio


More
----
See the README of sub package.

Documentation
-------------
[See documentation on godoc for more detail](https://godoc.org/github.com/shenwei356/bio/).

Copyright (c) 2013-2016, Wei Shen (shenwei356@gmail.com)

[MIT License](https://github.com/shenwei356/bio/blob/master/LICENSE)
