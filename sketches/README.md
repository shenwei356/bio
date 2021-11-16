# sketches

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/bio/sketches.svg)](https://pkg.go.dev/github.com/shenwei356/bio/sketches)


This package provides iterators for k-mer and k-mer sketches 
([Minimizer](https://academic.oup.com/bioinformatics/article/20/18/3363/202143),
 [Scaled MinHash](https://f1000research.com/articles/8-1006),
 [Closed Syncmers](https://peerj.com/articles/10805/)).
K-mers are either encoded (k<=32) or hashed (arbitrary k, using [ntHash](https://github.com/will-rowe/nthash)) into `uint64`.

Related projects:

- [kmers](https://github.com/shenwei356/kmers) provides manipulations for bit-packed k-mers (k<=32, encoded in `uint64`).
- [kmcp](https://github.com/shenwei356/kmcp) uses this package.

## Benchmark

CPU: AMD Ryzen 7 2700X Eight-Core Processor, 3.7 GHz

    $ go test . -bench=Bench* -benchmem \
        | grep Bench \
        | perl -pe 's/\s\s+/\t/g' \
        | csvtk cut -Ht -f 1,3-5 \
        | csvtk add-header -t -n test,time,memory,allocs \
        | csvtk pretty -t -r
 
                                          test           time     memory        allocs
    ------------------------------------------   ------------   --------   -----------
              BenchmarkKmerIterator/1.00_KB-16    11292 ns/op     0 B/op   0 allocs/op
              BenchmarkHashIterator/1.00_KB-16     7146 ns/op    24 B/op   1 allocs/op
           BenchmarkProteinIterator/1.00_KB-16    13985 ns/op   432 B/op   2 allocs/op

           BenchmarkMinimizerSketch/1.00_KB-16    58062 ns/op    48 B/op   2 allocs/op
             BenchmarkSyncmerSketch/1.00_KB-16   102475 ns/op   977 B/op   7 allocs/op
    BenchmarkProteinMinimizerSketch/1.00_KB-16    21617 ns/op   733 B/op   5 allocs/op


## History

This package was originally maintained in [unikmer](https://github.com/shenwei356/unikmer).
