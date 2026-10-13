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
    	reader.Close() // Please remember to call this !!!
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

## Input handling

In the examples below, `\n` denotes LF and `\r\n` denotes CRLF. Empty inputs
return `io.EOF` on the first call to `Read()` without returning a record.

| Input | Handling |
| --- | --- |
| Zero-byte file or stream | Accepted; 0 records. |
| One or more LF line breaks, e.g. `\n` or `\n\n\n` | Accepted; 0 records, including inputs spanning multiple read buffers. |
| Leading LF-only blank lines before a FASTA or FASTQ header | Accepted; records are parsed normally. |
| A space or tab before any header, alone or mixed with LF line breaks, e.g. ` `, ` \n`, `\t`, or `\n\t\n` | Rejected with `ErrNotFASTXFormat`. |
| Leading CRLF blank lines, alone or before a header, e.g. `\r\n` or `\r\n>id\nACGT\n` | Rejected with `ErrNotFASTXFormat`. |
| A header with leading spaces or tabs, e.g. ` >id\nACGT\n` | Rejected with `ErrNotFASTXFormat`. |
| Sequence data without a header, e.g. `ACGT\n` | Rejected with `ErrNotFASTXFormat`. |
| An initial UTF-8 BOM | Ignored; a BOM alone or followed only by LF line breaks yields 0 records. |
| A UTF-8 BOM after leading line breaks | Rejected with `ErrNotFASTXFormat`; only a BOM at the start is ignored. |
| `>`, `>\n`, or `\n\n>` | Accepted; 1 record with an empty ID and sequence. |
| `>id\n` | Accepted; 1 record with ID `id` and an empty sequence. |
| `>\nACGT\n` | Accepted; 1 record with an empty ID and sequence `ACGT`. |
| FASTA records using CRLF line endings, e.g. `>id\r\nACGT\r\n` | Accepted; line endings are removed from the parsed header and sequence. |
| A FASTA header or sequence without a final LF, e.g. `>id` or `>id\nACGT` | Accepted; 1 record. |
| Repeated IDs, including empty IDs, e.g. `>\n>\n` | Every record is returned; IDs are not deduplicated. |
| A valid gzip file whose decompressed content is empty or contains only LF line breaks | Accepted; 0 records. Whitespace handling is the same as for uncompressed input. |

