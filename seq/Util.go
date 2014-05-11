package seq

import "bytes"

// reverse a byte slice
func ReverseByteSlice(s []byte) []byte {
	// make a copy of s
	l := len(s)
	t := make([]byte, l)
	for i := 0; i < l; i++ {
		t[i] = s[i]
	}

	// reverse
	for i, j := 0, len(t)-1; i < j; i, j = i+1, j-1 {
		t[i], t[j] = t[j], t[i]
	}
	return t
}

// format sequence for ouput
func FormatSeq(seq []byte, width int) []byte {
	var buffer bytes.Buffer
	l := len(seq)
	lines := int(l / width)
	var start, end int
	for i := 0; i <= lines; i++ {
		start = i * width
		end = (i + 1) * width
		if end > l {
			end = l
		}

		buffer.Write(seq[start:end])
		buffer.WriteString("\n")
	}
	return buffer.Bytes()
}
