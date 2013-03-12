/*
 * Package for reading sequence from fasta file.
 * by Wei Shen (shenwei356@gmail.com)
 */

package main

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
)

func main() {
	file := "test.fa"
	nextSeq, err := FastaReader(file)
	if err != nil {
		fmt.Println(err)
		return
	}

	var head, seq string
	for {
		head, seq = nextSeq()
		if head == "" {
			break
		}

		fmt.Printf(">%s\n%s\n", head, seq)
	}
}

func FastaReader(file string) (func() (string, string), error) {
	fh, err := os.Open(file)

	if err != nil {
		recover()
		return nil, errors.New("[Error] Failed to open file (" + file + ").")
	}

	reader := bufio.NewReader(fh)
	var buffer bytes.Buffer
	var lastHead string

	return func() (head, seq string) {
		var str string
		var err error
		for {
			str, err = reader.ReadString('\n')

			if err != nil {
				if err == io.EOF { // EOF
					if buffer.Len() > 0 {
						buffer.WriteString(strings.TrimRight(str, "\r?\n"))
						seq = buffer.String()
						buffer.Reset()
						return lastHead, seq
					} else {
						err = io.EOF
						break
					}
				}
			}

			if strings.HasPrefix(str, ">") {
				thisHead := strings.TrimRight(str[1:], "\r?\n")
				if buffer.Len() > 0 { // no-first seq head
					seq = buffer.String()
					buffer.Reset()
					head, lastHead = lastHead, thisHead
					return head, seq
				} else { // first sequence head
					lastHead = thisHead
				}
			} else { // append sequence
				buffer.WriteString(strings.TrimRight(str, "\r?\n"))
			}
		}
		return
	}, err
}
