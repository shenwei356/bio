/*Package fastx seamlessly parses FASTA and FASTQ format file
This package seamlessly parses both FASTA and FASTQ formats.

## Examples


### Common operation

    import (
        "github.com/shenwei356/bio/seqio/fastx"
        "io"
    )

    reader, err := NewDefaultReader(file)
    checkErr(err)
    for {
        record, err := reader.Read()
        if err != nil {
            if err == io.EOF {
                break
            }
            checkErr(err)
            break
        }

        fmt.Print(record)
    }

***Note that***, similar to `bytes.Buffer.Bytes()` method,
the current record will change after your another call of this method.
So, you could use `record.Clone()` to make a copy.

### Asynchronously parsing

ChunkChan asynchronously reads FASTA/Q records, and returns a channel of
Record Chunk, from which you can easily access the records.
bufferSize is the number of buffered chunks, and chunkSize is the size
of records in a chunk.

    reader, err := NewDefaultReader(file)
    if err != nil {
        t.Error(t)
    }

    for chunk := range reader.ChunkChan(bufferSize, chunkSize) {
        if chunk.Err != nil {
            t.Error(chunk.Err)
        }

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

    reader, err := NewReader(seq.DNA, file, "^([^\s]+)\s?")

*/
package fastx
