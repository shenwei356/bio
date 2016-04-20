## Example


### Common operation


    import (
        "github.com/shenwei356/bio/seq"
        "github.com/shenwei356/bio/seqio/fastx"
    )

    fastxReader, err := fasta.NewReader(seq.DNAredundant, file, threads, chunkSize, "")
    checkError(err)

    for chunk := range fastxReader.Ch {
        checkError(chunk.Err)

        for _, record := range chunk.Data {
            fmt.Printf(">%s\n%s\n", record.Name, record.FormatSeq(lineWidth))
        }
    }


### Cancellation

Reading the first record.

**Note that `range chanel` is buffered, therefore `for-select-case` is used.**


    fastxReader, err := fasta.NewReader(seq.Unlimit, file, 1, 1, "")
    checkError(err)

    // note that range is bufferd. using range will be failed
    // for chunk := range fastxReader.Ch {
    for {
        select {
        case chunk := <-fastxReader.Ch:
            checkError(chunk.Err)

            // do some thing

            reader.Cancel()
        default:
        }
    }

### Guessing Alphabet

If alphabet is nil, it will guess alphabet by the first 
`seq.AlphabetGuessSeqLenghtThreshold` (default 10000, 0 for whole seq)
letters of first record.
    
    
    fastxReader, err := fasta.NewReader(nil, file, threads, chunkSize, "")
    checkError(err)
