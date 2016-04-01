/*Package fasta parses fasta format file

Examples:

1. Common operation

    import (
        "github.com/shenwei356/bio/seq"
        "github.com/shenwei356/bio/seqio/fasta"
    )

    fastaReader, err := fasta.NewFastaReader(seq.DNAredundant, file, threads, chunkSize, "")
    checkError(err)

    for chunk := range fastaReader.Ch {
        checkError(chunk.Err)

        for _, record := range chunk.Data {
            fmt.Printf(">%s\n%s\n", record.Name, record.FormatSeq(lineWidth))
        }
    }

2. Cancellation

Reading the first record.

**Note that `range chanel` is buffered, therefore `for-select-case` is used.**


    fastaReader, err := fasta.NewFastaReader(seq.Unlimit, file, 1, 1, "")
    checkError(err)

    // note that range is bufferd. using range will be failed
    // for chunk := range fastaReader.Ch {
    for {
        select {
        case chunk := <-fastaReader.Ch:
            checkError(chunk.Err)

            // do some thing

            reader.Cancel()
        default:
        }
    }

3. Guessing Alphabet

If alphabet is nil, it will guess alphabet by the first record
    
    fastaReader, err := fasta.NewFastaReader(nil, file, threads, chunkSize, "")
    checkError(err)

*/
package fasta
