# Changelog

### v0.9.1 - 2023-11-08

- seqio/fastx.Reader: recycle []byte buffer to save memory for reading a large number of sequences.

### v0.9.0 - 2023-06-25

- util/LengthStats: new method to compute N50 for any number
- seq/alphabet: faster with asciiset

### v0.8.4 - 2023-02-14

- seqio/fai: report error for non-fasta files

### v0.8.3 - 2022-12-02

- util.LengthStats: fix computing Q1 and Q3 for one element.

### v0.8.2 - 2022-11-16

- faidx: allow empty lines at the end of sequences

### v0.8.1 - 2022-09-06

- fastx: fix concurrency bug of `record.FormatToWriter()`.

### v0.8.0 - 2022-09-06

- sketches: added an iterator of SimHash.

### v0.7.1 - 2022-04-19

- taxdump: allow reading empty merged.dmp and delnodes.dmp

### v0.6.4 - 2022-03-13

- update xopen version

### v0.6.3 - 2022-03-12

- use new versin of xopen which support .xz and .zst

### v0.6.2 - 2021-12-01

- taxdump: more robust

### v0.6.1 - 2021-11-16

- taxdump: fix pkg name in errors

### v0.6.0 - 2021-11-08

- move sketches and taxdump packages from unikmer to here
