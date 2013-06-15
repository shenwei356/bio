/*
 * Package for reading sequence from fasta file.
 * by Wei Shen (shenwei356@gmail.com)
 */

package seq

import (
	"bufio"
	"bytes"
	"errors"
	"io"
	"os"
	"strings"
)

/*
FastaReader is a fasta file parser, which returns a function that
returns a pair of head and sequence when it was called.
Example:
	NextSeq, err := seq.FastaReader("test.fa")
	if err != nil {
		recover()
		fmt.Println(err)
		return
	}

	for {
		head, seq, err := NextSeq()
		if err != nil {
			// fmt.Println(err)
			break
		}

		fmt.Printf(">%s\n%s\n", head, seq)
	}
*/
func FastaReader(file string) (func() (string, string, error), error) {
	fh, err := os.Open(file)
	if err != nil {
		recover()
		return nil, errors.New("Failed to open file: " + file)
	}

	reader := bufio.NewReader(fh)
	var buffer bytes.Buffer
	var lastHead string
	var fileHandlerClosed bool = false

	return func() (head, seq string, err error) {
		if fileHandlerClosed {
			return "", "", io.EOF
		}

		var str string
		for {
			str, err = reader.ReadString('\n')
			if err == io.EOF {
				seq = buffer.String()
				buffer.Reset()
				head = lastHead
				fh.Close()
				fileHandlerClosed = true
				err = nil
				return
			}

			if strings.HasPrefix(str, ">") {
				thisHead := strings.TrimRight(str[1:], "\r?\n")
				if buffer.Len() > 0 { // no-first seq head
					seq = buffer.String()
					buffer.Reset()
					head, lastHead = lastHead, thisHead
					return
				} else { // first sequence head
					lastHead = thisHead
				}
			} else { // append sequence
				buffer.WriteString(strings.TrimRight(str, "\r?\n"))
			}
		}
	}, err
}
