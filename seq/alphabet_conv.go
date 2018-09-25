// Copyright © 2018
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

package seq

import (
	"bytes"
	"errors"
)

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes
var condonTables = map[int]map[string]byte{}
var frameShiftingOpts = map[string]bool{}

func init() {

	frameShiftingOpts = map[string]bool{
		"F1": true,
		"F2": true,
		"F3": true,
		"R1": true,
		"R2": true,
		"R3": true,
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1
	condonTables[1] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG2
	condonTables[2] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'M', "ACA": 'T', "AAA": 'K', "AGA": '*',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": '*',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG3
	condonTables[3] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'T', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'T', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'T', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'T', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'T', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'M', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG4
	condonTables[4] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'S',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG5
	condonTables[5] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'M', "ACA": 'T', "AAA": 'K', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'S',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG6
	condonTables[6] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'Q', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'Q', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG9
	condonTables[9] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'N', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'S',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG10
	condonTables[10] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'C',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG11
	condonTables[11] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG12
	condonTables[12] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'S', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG13
	condonTables[13] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'M', "ACA": 'T', "AAA": 'K', "AGA": 'G',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'G',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG14
	condonTables[14] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'Y', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'N', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'S',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG16
	condonTables[16] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'L', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG21
	condonTables[21] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'M', "ACA": 'T', "AAA": 'N', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'S',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG22
	condonTables[22] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": '*', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'L', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG23
	condonTables[23] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": '*', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG24
	condonTables[24] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'S',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'K',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG25
	condonTables[25] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": 'G',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG26
	condonTables[26] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": '*', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": '*', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'A', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG27
	condonTables[27] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'G', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'G', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG28
	condonTables[28] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'G', "TGA": 'W',
		"TTG": 'L', "TCG": 'S', "TAG": 'G', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG29
	condonTables[29] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'Y', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'Y', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG30
	condonTables[30] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'G', "TGA": '*',
		"TTG": 'L', "TCG": 'S', "TAG": 'G', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG31
	condonTables[31] = map[string]byte{
		"TTT": 'F', "TCT": 'S', "TAT": 'Y', "TGT": 'C',
		"TTC": 'F', "TCC": 'S', "TAC": 'Y', "TGC": 'C',
		"TTA": 'L', "TCA": 'S', "TAA": 'G', "TGA": 'Y',
		"TTG": 'L', "TCG": 'S', "TAG": 'G', "TGG": 'W',
		"CTT": 'L', "CCT": 'P', "CAT": 'H', "CGT": 'R',
		"CTC": 'L', "CCC": 'P', "CAC": 'H', "CGC": 'R',
		"CTA": 'L', "CCA": 'P', "CAA": 'Q', "CGA": 'R',
		"CTG": 'L', "CCG": 'P', "CAG": 'Q', "CGG": 'R',
		"ATT": 'I', "ACT": 'T', "AAT": 'N', "AGT": 'S',
		"ATC": 'I', "ACC": 'T', "AAC": 'N', "AGC": 'S',
		"ATA": 'I', "ACA": 'T', "AAA": 'K', "AGA": 'R',
		"ATG": 'M', "ACG": 'T', "AAG": 'K', "AGG": 'R',
		"GTT": 'V', "GCT": 'A', "GAT": 'D', "GGT": 'G',
		"GTC": 'V', "GCC": 'A', "GAC": 'D', "GGC": 'G',
		"GTA": 'V', "GCA": 'A', "GAA": 'E', "GGA": 'G',
		"GTG": 'V', "GCG": 'A', "GAG": 'E', "GGG": 'G',
	}

}

// DNAToProtein accepts a DNA Sequence, frameShiftingOptions, translation table from NIBH standard and returns a protein sequence
func DNAToProtein(dna *Seq, frameShiftingOptions string, transl_table int) ([]byte, error) {

	if dna.Seq == nil {
		return nil, errors.New("seq.DNAToProtein: input sequence is nil")
	}

	if frameShiftingOptions == "" || !frameShiftingOpts[frameShiftingOptions] {
		return nil, errors.New("seq.DNAToProtein: frameShiftingOptions is missing or invalid")
	}

	if _, found := condonTables[transl_table]; !found {
		return nil, errors.New("seq.DNAToProtein: unavailable condon table")
	}

	var i int = 0
	if frameShiftingOptions == "F2" || frameShiftingOptions == "R2" {
		i = 1
	}
	if frameShiftingOptions == "F3" || frameShiftingOptions == "R3" {
		i = 2
	}
	var ok bool
	var seq []byte = dna.Seq
	var protein bytes.Buffer
	if frameShiftingOptions == "R1" || frameShiftingOptions == "R2" || frameShiftingOptions == "R3" {
		seq = dna.Reverse().Seq
	}
	for i < len(seq)-2 {
		if _, ok = condonTables[transl_table][string(seq[i:i+3])]; ok {
			err := protein.WriteByte(condonTables[transl_table][string(seq[i:i+3])])
			if err != nil {
				return nil, err
			}
		}
		i += 3
	}
	return protein.Bytes(), nil
}
