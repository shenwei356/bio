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

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG1
var rnaCodonToAminoAcidTranslTable1 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG2
var rnaCodonToAminoAcidTranslTable2 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "M", "ACA": "T", "AAA": "K", "AGA": "*",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "*",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG3
var rnaCodonToAminoAcidTranslTable3 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "T", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "T", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "T", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "T", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "T", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "M", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG4
var rnaCodonToAminoAcidTranslTable4 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "S",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG5
var rnaCodonToAminoAcidTranslTable5 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "M", "ACA": "T", "AAA": "K", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "S",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG6
var rnaCodonToAminoAcidTranslTable6 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "Q", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "Q", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG9
var rnaCodonToAminoAcidTranslTable9 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "N", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "S",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG10
var rnaCodonToAminoAcidTranslTable10 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "C",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG11
var rnaCodonToAminoAcidTranslTable11 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG12
var rnaCodonToAminoAcidTranslTable12 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "S", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG13
var rnaCodonToAminoAcidTranslTable13 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "M", "ACA": "T", "AAA": "K", "AGA": "G",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "G",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG14
var rnaCodonToAminoAcidTranslTable14 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "Y", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "N", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "S",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG16
var rnaCodonToAminoAcidTranslTable16 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "L", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG21
var rnaCodonToAminoAcidTranslTable21 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "M", "ACA": "T", "AAA": "N", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "S",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG22
var rnaCodonToAminoAcidTranslTable22 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "*", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "L", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG23
var rnaCodonToAminoAcidTranslTable23 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "*", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG24
var rnaCodonToAminoAcidTranslTable24 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "S",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "K",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG25
var rnaCodonToAminoAcidTranslTable25 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "G",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG26
var rnaCodonToAminoAcidTranslTable26 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "*", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "*", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "A", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG27
var rnaCodonToAminoAcidTranslTable27 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "G", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "G", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG28
var rnaCodonToAminoAcidTranslTable28 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "G", "UGA": "W",
	"UUG": "L", "UCG": "S", "UAG": "G", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG29
var rnaCodonToAminoAcidTranslTable29 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "Y", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "Y", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG30
var rnaCodonToAminoAcidTranslTable30 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "G", "UGA": "*",
	"UUG": "L", "UCG": "S", "UAG": "G", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

//https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/index.cgi?chapter=tgencodes#SG31
var rnaCodonToAminoAcidTranslTable31 = map[string]string{
	"UUU": "F", "UCU": "S", "UAU": "Y", "UGU": "C",
	"UUC": "F", "UCC": "S", "UAC": "Y", "UGC": "C",
	"UUA": "L", "UCA": "S", "UAA": "G", "UGA": "Y",
	"UUG": "L", "UCG": "S", "UAG": "G", "UGG": "W",
	"CUU": "L", "CCU": "P", "CAU": "H", "CGU": "R",
	"CUC": "L", "CCC": "P", "CAC": "H", "CGC": "R",
	"CUA": "L", "CCA": "P", "CAA": "Q", "CGA": "R",
	"CUG": "L", "CCG": "P", "CAG": "Q", "CGG": "R",
	"AUU": "I", "ACU": "T", "AAU": "N", "AGU": "S",
	"AUC": "I", "ACC": "T", "AAC": "N", "AGC": "S",
	"AUA": "I", "ACA": "T", "AAA": "K", "AGA": "R",
	"AUG": "M", "ACG": "T", "AAG": "K", "AGG": "R",
	"GUU": "V", "GCU": "A", "GAU": "D", "GGU": "G",
	"GUC": "V", "GCC": "A", "GAC": "D", "GGC": "G",
	"GUA": "V", "GCA": "A", "GAA": "E", "GGA": "G",
	"GUG": "V", "GCG": "A", "GAG": "E", "GGG": "G",
}

// RNAToProtein accepts a RNA Sequence, translation table from NIBH standard and returns a protein sequence
func RNAToProtein(rna []byte, transl_table int) (protein []byte, err error) {

	if rna == nil {
		return nil, errors.New("seq.RNAToProtein: input sequence is nil")
	}
	var acceptable_transl_tables []int
	acceptable_transl_tables = append(acceptable_transl_tables, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31)
	var found bool = false
	for _, t := range acceptable_transl_tables {
		if transl_table == t {
			found = true
			break
		}
	}
	if !found {
		return nil, errors.New("seq.RNAToProtein: translation table is not found")
	}

	var last3NTs bytes.Buffer
	for _, nt := range rna {

		err := last3NTs.WriteByte(nt)
		if err != nil {
			return nil, err
		}

		if len(last3NTs.String()) == 3 {
			switch transl_table {
			case 1:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable1[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 2:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable2[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 3:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable3[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 4:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable4[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 5:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable5[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 6:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable6[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 9:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable9[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 10:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable10[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 11:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable11[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 12:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable12[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 13:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable13[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 14:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable14[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 16:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable16[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 21:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable21[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 22:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable22[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 23:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable23[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 24:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable24[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 25:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable25[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 26:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable26[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 27:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable27[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 28:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable28[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 29:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable29[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 30:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable30[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			case 31:
				if aminoacid, ok := rnaCodonToAminoAcidTranslTable31[last3NTs.String()]; ok {
					protein = append(protein, []byte(aminoacid)...)
				}
			}
			last3NTs.Reset()
		}
	}
	return protein, nil
}
