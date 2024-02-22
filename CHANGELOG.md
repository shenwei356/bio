# Changelog

### v0.13.1 - 2024-02-22

- util: fix computation of L50.

### v0.13.0 - 2024-02-19

- seq: remove the global variable: `seq.ValidateWholeSeq`.

### v0.12.1 - 2024-01-03

- seqio/fai: when using the whole FASTA header as the sequence ID, replace possible tabs in FASTA header with spaces.

### v0.12.0 - 2023-12-04

- seqio/fastx.Reader: reuse the reader with an object pool, this requires users to call reader.Close() after using.
  The benefit is that it reduces memory when handling of a lot of sequences.

### v0.11.0 - 2023-12-03

Do not use this version!

- seqio/fastx.Reader: delete reader.Recycle() to avoid API changes.

### v0.10.0 - 2023-12-03

- seqio/fastx.Reader: reuse reader with object pool, this requires users to call reader.Recycle() after using.

### v0.9.3 - 2023-11-11

- seqio/fastx.Reader: fix a panic of nil pointer when some files has file size of 0.

### v0.9.2 - 2023-11-10

- seq/alphabet.IsValid: fix panic: close of closed channel

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
