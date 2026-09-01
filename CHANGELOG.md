# Changelog

### v0.15.1 - 2026-09-01

- **fastx: reduce time and memory when reading large FASTA records**.
- fai: reduce time and memory when creating indexes for large FASTA files.
- seq: add `seq.AvgQualOfRegion` to compute the average quality of a region.
- seq: avoid allocating parsed quality values in `seq.AvgQual`.
- seq: reuse sequence and quality buffers in `seq.RemoveGapsInplace`.
- seq: increase `seq.ValidSeqLengthThreshold` to 64 KiB to reduce parallel validation overhead.
- seq: speed up copying parsed quality values in `seq.SubSeq`.
- sketches: safely return exhausted sketches to the correct object pools.
- sketches: reuse sketch buffers to reduce allocations.
- taxdump: reduce time and allocations in `Taxonomy.LCA`.

### v0.15.0 - 2026-08-07

- fastx: avoid data race when using different id regexps.
- fastx: close underlying readers before returning readers to the pool.
- fastx: handle empty IO readers without panicking.
- fastx: better handle leading blank lines across input buffers.
- fastx: restore the input buffer size after short reads.
- seq: avoid data race when validating long sequences in parallel.
- fai: avoid leaking custom ID regexp state between index creations.
- fai: return an error for out-of-range base queries.
- fai: make closing indexes idempotent.
- fai: close files when memory mapping fails.
- fai: preserve the access mode of each opened index.
- fai: reject empty FASTA files and invalid ID regexps.
- sketches: reset pooled iterator mode state.
- sketches: avoid returning hash iterators to the pool more than once.
- sketches: avoid mutating input sequences during non-canonical k-mer iteration.
- taxdump: fix loading names when the taxid column has the highest index.

### v0.14.0 - 2026-05-13

- fastx: allow wrapping FASTQ sequences.

### v0.13.9 - 2025-10-30

- fix reading records with both empty id and sequence.

### v0.13.8 - 2025-08-29

- update dependencies.

### v0.13.7 - 2025-08-12

- fastx: fix a bug that it did not report error when the invalid record is not the last one.

### v0.13.6 - 2024-09-12

- fix sequence ID parsing with the default regular expression (in this case, we actually use bytes.Index instead) for a rare case: "xxx\tyyy zzz" was wrongly parsed as "xxx\tyyy".

### v0.13.5 - 2024-09-05

- fix sequence ID parsing with the default regular expression (in this case, we actually use bytes.Index instead) for a rare case: "xxx\tyyy zzz" was wrongly parsed as "xxx\tyyy".

### v0.13.4 - 2024-07-15

- fix seq.SubSeq for (-start,-end), where start > len(seq).

### v0.13.3 - 2024-03-11

- fix dependency.

### v0.13.2 - 2024-03-11

- seq: increase default value of `seq.ComplementSeqLenThreshold` from 1000 to 1000000, which means only parallelize complement computation for sequences longer than 1Mb.

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
