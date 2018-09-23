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
var condonTables = map[int]map[string][]byte{}

func initializeCodonTables() {

	condonTables[1] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG2
	condonTables[2] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("M"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("*"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("*"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG3
	condonTables[3] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("T"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("T"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("T"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("T"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("T"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("M"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG4
	condonTables[4] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("S"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG5
	condonTables[5] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("M"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("S"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG6
	condonTables[6] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("Q"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("Q"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG9
	condonTables[9] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("N"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("S"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG10
	condonTables[10] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("C"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG11
	condonTables[11] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG12
	condonTables[12] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("S"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG13
	condonTables[13] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("M"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("G"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("G"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG14
	condonTables[14] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("Y"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("N"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("S"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG16
	condonTables[16] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("L"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG21
	condonTables[21] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("M"), "ACA": []byte("T"), "AAA": []byte("N"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("S"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG22
	condonTables[22] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("*"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("L"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG23
	condonTables[23] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("*"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG24
	condonTables[24] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("S"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("K"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG25
	condonTables[25] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("G"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG26
	condonTables[26] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("*"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("*"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("A"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG27
	condonTables[27] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("G"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("G"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG28
	condonTables[28] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("G"), "TGA": []byte("W"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("G"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG29
	condonTables[29] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("Y"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("Y"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG30
	condonTables[30] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("G"), "TGA": []byte("*"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("G"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

	//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG31
	condonTables[31] = map[string][]byte{
		"TTT": []byte("F"), "TCT": []byte("S"), "TAT": []byte("Y"), "TGT": []byte("C"),
		"TTC": []byte("F"), "TCC": []byte("S"), "TAC": []byte("Y"), "TGC": []byte("C"),
		"TTA": []byte("L"), "TCA": []byte("S"), "TAA": []byte("G"), "TGA": []byte("Y"),
		"TTG": []byte("L"), "TCG": []byte("S"), "TAG": []byte("G"), "TGG": []byte("W"),
		"CTT": []byte("L"), "CCT": []byte("P"), "CAT": []byte("H"), "CGT": []byte("R"),
		"CTC": []byte("L"), "CCC": []byte("P"), "CAC": []byte("H"), "CGC": []byte("R"),
		"CTA": []byte("L"), "CCA": []byte("P"), "CAA": []byte("Q"), "CGA": []byte("R"),
		"CTG": []byte("L"), "CCG": []byte("P"), "CAG": []byte("Q"), "CGG": []byte("R"),
		"ATT": []byte("I"), "ACT": []byte("T"), "AAT": []byte("N"), "AGT": []byte("S"),
		"ATC": []byte("I"), "ACC": []byte("T"), "AAC": []byte("N"), "AGC": []byte("S"),
		"ATA": []byte("I"), "ACA": []byte("T"), "AAA": []byte("K"), "AGA": []byte("R"),
		"ATG": []byte("M"), "ACG": []byte("T"), "AAG": []byte("K"), "AGG": []byte("R"),
		"GTT": []byte("V"), "GCT": []byte("A"), "GAT": []byte("D"), "GGT": []byte("G"),
		"GTC": []byte("V"), "GCC": []byte("A"), "GAC": []byte("D"), "GGC": []byte("G"),
		"GTA": []byte("V"), "GCA": []byte("A"), "GAA": []byte("E"), "GGA": []byte("G"),
		"GTG": []byte("V"), "GCG": []byte("A"), "GAG": []byte("E"), "GGG": []byte("G"),
	}

}

func validateTranslationTable(transl_table int) error {
	acceptable_transl_tables := []int{1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}
	var found bool = false
	for _, t := range acceptable_transl_tables {
		if transl_table == t {
			found = true
			break
		}
	}
	if !found {
		return errors.New("seq.DNAToProtein: translation table is not found")
	}
	return nil
}

// DNAToProtein accepts a DNA Sequence, translation table from NIBH standard and returns a protein sequence
func DNAToProtein(dna []byte, transl_table int) ([]byte, error) {

	if dna == nil {
		return nil, errors.New("seq.DNAToProtein: input sequence is nil")
	}

	err := validateTranslationTable(transl_table)
	if err != nil {
		return nil, err
	}

	initializeCodonTables()

	var protein bytes.Buffer
	for i := 0; i < len(dna)-2; i += 3 {
		if aminoacid, ok := condonTables[transl_table][string(dna[i:i+3])]; ok {
			_, err := protein.Write(aminoacid)
			if err != nil {
				return nil, err
			}
		}
	}
	return protein.Bytes(), nil
}
