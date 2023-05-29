

## FASTA/Q reading and writing

***bio/seqio/fastx has a high performance close to the famous C lib
[`kseq.h`](https://github.com/attractivechaos/klib/blob/master/kseq.h).***

To test the performance, three datasets and their gzip-compressed file are used:

- dataset_A, bacteria genomes, 2.7G
- dataset_B, human genome,     2.9G
- dataset_C, Illumina reads,   2.2G

Summary by [`seqkit`](https://github.com/shenwei356/seqkit):

    file           seq_format   seq_type   num_seqs   min_len        avg_len       max_len
    dataset_A.fa   FASTA        DNA          67,748        56       41,442.5     5,976,145
    dataset_B.fa   FASTA        DNA             194       970   15,978,096.5   248,956,422
    dataset_C.fq   FASTQ        DNA       9,186,045       100            100           100

[`seqtk`](https://github.com/lh3/seqtk/)
(Version [1.3-r119-dirty](https://github.com/lh3/seqtk/commit/f6ea81cc30b9232e244dffa94187114275389132),
using `kseq.h`)
and [`seqkit`](https://github.com/shenwei356/seqkit)
(Version [v2.4.0](https://github.com/shenwei356/seqkit/releases/tag/v2.4.0),
using this package) were used to test.
**Note** that `seqtk` does not support wrapped (fixed line width) ouputing, so `seqkit` uses
`-w 0` to disable outputing wrapping.
Script [`memusg`](https://github.com/shenwei356/memusg) is used to assess running time
and peak memory usage.

[Commands](https://github.com/shenwei356/bio/blob/master/benchmark/)

Tests were repeated 4 times and average time and memory usage were computed.

Results:

<img src="benchmark.tsv.png" alt="" width="700" align="center" />

Notes:

- `seqkit` uses 4 threads by default.
- `seqkit_t1` uses 1 thread.
- `seqtk` is single-threaded.
- `seqtk+gzip`: `seqtk` pipes data to the single-threaded `gzip`.
- `seqtk+pigz`: `seqtk` pipes data to the multithreaded `pigz` which uses 4 threads here.

 
## Run

    ./run.pl -n 4 run_benchmark_*.sh --outfile benchmark.tsv

    # PLOT
    ./plot.sh
