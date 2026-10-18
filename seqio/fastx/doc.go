/*
Package fastx reads and writes FASTA and FASTQ records. NewDefaultReader
detects the format and guesses the alphabet from the first record. It uses the
first space- or tab-delimited field of the header as the record ID.

# Reading and writing records

This program reads standard input and writes records to standard output without
line wrapping. It reports read and output errors and closes both streams.

	package main

	import (
		"errors"
		"fmt"
		"io"
		"os"

		"github.com/shenwei356/bio/seqio/fastx"
		"github.com/shenwei356/xopen"
	)

	func main() {
		if err := run(); err != nil {
			fmt.Fprintln(os.Stderr, err)
			os.Exit(1)
		}
	}

	// run returns errors after deferred cleanup has finished.
	func run() (err error) {
		reader, err := fastx.NewDefaultReader("-")
		if err != nil {
			return err
		}
		defer reader.Close()

		outfh, err := xopen.Wopen("-")
		if err != nil {
			return err
		}
		defer func() {
			// Close flushes buffered output and reports write errors.
			err = errors.Join(err, outfh.Close())
		}()

		for {
			record, err := reader.Read()
			if err == io.EOF {
				return nil
			}
			if err != nil {
				return err
			}
			// Keep zero-length FASTQ records in FASTQ format.
			fastx.ForcelyOutputFastq = reader.IsFastq
			record.FormatToWriter(outfh, 0)
		}
	}

`Read` reuses the record and its backing buffers. Call `record.Clone()` before
the next `Read` if you need to retain a record. Always call `reader.Close()`
when finished, including after an error or an early stop; do not use the reader
after closing it.

`reader.IsFastq` is determined by the first successful `Read`.
`ForcelyOutputFastq` is a package-level setting used by the formatting methods;
this example changes it while processing a single input synchronously. Do not
change it while other goroutines are formatting records.

Sequence validation is enabled by default. For known DNA input, pass an explicit
alphabet to `NewReader` (see below); automatic alphabet detection can fall back
to `seq.Unlimit`, which accepts any sequence bytes. Setting `seq.ValidateSeq = false`
disables alphabet validation for all readers. Configure it before starting any
reader goroutines.

# Asynchronous parsing

`ChunkChan(bufferSize, chunkSize)` reads records in a background goroutine.
`bufferSize` is the number of chunks the channel can buffer; `chunkSize` is the
maximum number of records per chunk. A chunk can contain valid records followed
by a read error, so process `chunk.Data` before handling `chunk.Err`.

This function counts bases, returning the count read so far if parsing fails:

	package main

	import "github.com/shenwei356/bio/seqio/fastx"

	func countBases(file string) (int64, error) {
		reader, err := fastx.NewDefaultReader(file)
		if err != nil {
			return 0, err
		}
		defer reader.Close()

		var total int64
		var readErr error
		for chunk := range reader.ChunkChan(2, 100) {
			for _, record := range chunk.Data {
				total += int64(len(record.Seq.Seq))
			}
			if chunk.Err != nil {
				readErr = chunk.Err
			}
		}
		return total, readErr
	}

Records from `ChunkChan` are already cloned and can be retained without another
`Clone`. Consume the channel until it closes before closing the reader. Stopping
early without draining the channel can leave its producer blocked; use `Read`
for workflows that need to stop early. Do not call `Read` concurrently with
`ChunkChan` on the same reader.

# Custom alphabet and identifier regular expression

Pass `nil` as the alphabet to guess it from the first record, or use an explicit
alphabet such as `seq.DNAredundant` for DNA with IUPAC ambiguity codes.
The ID regular expression must contain a capturing group; the first group is
used as the ID. An empty pattern uses the default `^(\S+)\s?`. If a custom
pattern does not match, the entire header becomes the ID.

This function reads the first ID from a DNA file, taking the text before the
first `|`, space, or tab. For example, `>seq1|annotation description` yields `seq1`.
Use a raw string literal (backticks) for regular expressions containing backslashes.

	package main

	import (
		"github.com/shenwei356/bio/seq"
		"github.com/shenwei356/bio/seqio/fastx"
	)

	func firstID(file string) (string, error) {
		reader, err := fastx.NewReader(seq.DNAredundant, file, `^([^|\s]+)`)
		if err != nil {
			return "", err
		}
		defer reader.Close()

		record, err := reader.Read()
		if err != nil {
			return "", err
		}
		return string(record.ID), nil
	}

# Leading FASTA comments

Lines starting with ';' before the first FASTA header are skipped, including
FHR metadata lines starting with ';~'. These comments are discarded and are not
preserved when records are written. Comment-only input returns io.EOF without
a record. Leading comments before FASTQ are rejected. Lines starting with ';'
after the first FASTA header are not skipped.
*/
package fastx
