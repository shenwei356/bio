# kmers

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/bio/kmers.svg)](https://pkg.go.dev/github.com/shenwei356/bio/kmers)


This package provides manipulations of bit-packed k-mers and iterators for k-mer sketches 
([Minimizer](https://academic.oup.com/bioinformatics/article/20/18/3363/202143),
 [Scaled MinHash](https://f1000research.com/articles/8-1006),
 [Closed Syncmers](https://peerj.com/articles/10805/)).
K-mers are either encoded (k<=32) or hashed (arbitrary k, using [ntHash](https://github.com/will-rowe/nthash)) into `uint64`.

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
                         BenchmarkEncodeK32-16    19.67 ns/op     0 B/op   0 allocs/op
           BenchmarkEncodeFromFormerKmerK32-16    7.692 ns/op     0 B/op   0 allocs/op
       BenchmarkMustEncodeFromFormerKmerK32-16    2.008 ns/op     0 B/op   0 allocs/op
                         BenchmarkDecodeK32-16    80.73 ns/op    32 B/op   1 allocs/op
                     BenchmarkMustDecodeK32-16    76.93 ns/op    32 B/op   1 allocs/op
                            BenchmarkRevK32-16    3.617 ns/op     0 B/op   0 allocs/op
                           BenchmarkCompK32-16   0.7999 ns/op     0 B/op   0 allocs/op
                        BenchmarkRevCompK32-16    3.814 ns/op     0 B/op   0 allocs/op
                       BenchmarkCannonalK32-16    4.147 ns/op     0 B/op   0 allocs/op

              BenchmarkKmerIterator/1.00_KB-16    11292 ns/op     0 B/op   0 allocs/op
              BenchmarkHashIterator/1.00_KB-16     7146 ns/op    24 B/op   1 allocs/op
           BenchmarkProteinIterator/1.00_KB-16    13985 ns/op   432 B/op   2 allocs/op

           BenchmarkMinimizerSketch/1.00_KB-16    58062 ns/op    48 B/op   2 allocs/op
             BenchmarkSyncmerSketch/1.00_KB-16   102475 ns/op   977 B/op   7 allocs/op
    BenchmarkProteinMinimizerSketch/1.00_KB-16    21617 ns/op   733 B/op   5 allocs/op


## History

This package was originally maintained in [unikmer](https://github.com/shenwei356/unikmer).
