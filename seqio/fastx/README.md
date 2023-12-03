# fastx

[![Go Reference](https://pkg.go.dev/badge/github.com/shenwei356/bio/seqio/fastx.svg)](https://pkg.go.dev/github.com/shenwei356/bio/seqio/fastx)

This package seamlessly parses both FASTA and FASTQ formats.


## Examples


### Common operation

    package main

    import (
    	"fmt"
    	"io"
    	"os"

    	// "github.com/shenwei356/bio/seq"
    	"github.com/shenwei356/bio/seqio/fastx"
    	"github.com/shenwei356/xopen"
    )

    func main() {
    	// use buffered out stream for output
    	outfh, err := xopen.Wopen("-") // "-" for STDOUT
    	checkError(err)
    	defer outfh.Close()

    	// disable sequence validation could reduce time when reading large sequences
    	// seq.ValidateSeq = false

    	reader, err := fastx.NewDefaultReader("-")
    	checkError(err)

    	var record *fastx.Record
    	for {
    		record, err = reader.Read()
    		if err != nil {
    			if err == io.EOF {
    				break
    			}
    			checkError(err)
    			break
    		}

    		// fmt is slow for output, because it's not buffered
    		// fmt.Printf("%s", record.Format(0))

    		record.FormatToWriter(outfh, 0)
    	}

    	reader.Recycle() // remember to recycle the reader.
    }

    func checkError(err error) {
    	if err != nil {
    		fmt.Fprintln(os.Stderr, err)
    		os.Exit(1)
    	}
    }


***Note that***, similar with `bytes.Buffer.Bytes()` method,
the current record will change after your another call of this method.
You may use `record.Clone()` to make a copy.

### Asynchronously parsing

`ChunkChan` asynchronously reads FASTA/Q records, and returns a channel of
Record Chunk, from which you can easily access the records.
`bufferSize` is the number of buffered chunks, and `chunkSize` is the size
of records in a chunk.

    reader, err := fastx.NewDefaultReader(file)
    checkError(err)

    for chunk := range reader.ChunkChan(bufferSize, chunkSize) {
        checkError(chunk.Err)

        for _, record := range chunk.Data {
            fmt.Print(record)
        }
    }

***Note that***, these's no need to clone the record by `record.Clone()` here.

### Custom alphabet and identifier regular expression

    import (
        "github.com/shenwei356/bio/seq"
        "github.com/shenwei356/bio/seqio/fastx"
    )

    reader, err := fastx.NewReader(seq.DNA, file, "^([^\s]+)\s?")
