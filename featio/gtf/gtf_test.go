package gtf

import (
	"fmt"
	"os"
	"path/filepath"
	"slices"
	"testing"
)

// TestGTFAttributeParsing checks final attributes, quoted values, and missing values.
func TestGTFAttributeParsing(t *testing.T) {
	for _, test := range []struct {
		name, attributes string
		want             []Attribute
	}{
		{"final semicolon", `gene_id "g1"; gene_name "A";`, []Attribute{{"gene_id", "g1"}, {"gene_name", "A"}}},
		{"no final semicolon", `gene_name "A"; gene_id "g1"`, []Attribute{{"gene_name", "A"}, {"gene_id", "g1"}}},
		{"single attribute", `gene_id "g1"`, []Attribute{{"gene_id", "g1"}}},
		{"trailing whitespace", `gene_id "g1";  `, []Attribute{{"gene_id", "g1"}}},
		{"no separator space", `gene_id "g1";gene_name "A";`, []Attribute{{"gene_id", "g1"}, {"gene_name", "A"}}},
		{"multiple spaces", `  gene_id   "g1";   gene_name  "A";`, []Attribute{{"gene_id", "g1"}, {"gene_name", "A"}}},
		{"quoted spaces", `gene_name "alpha beta"; gene_id "g1"; `, []Attribute{{"gene_name", "alpha beta"}, {"gene_id", "g1"}}},
		{"quoted semicolon", `gene_name "A;B"; gene_id "g1"; `, []Attribute{{"gene_name", "A;B"}, {"gene_id", "g1"}}},
		{"quoted semicolon and space", `gene_name "A; B"; gene_id "g1";`, []Attribute{{"gene_name", "A; B"}, {"gene_id", "g1"}}},
		{"escaped quote", `gene_name "A\";B"; gene_id "g1";`, []Attribute{{"gene_name", `A\";B`}, {"gene_id", "g1"}}},
		{"empty quoted value", `gene_id ""; gene_name "A";`, []Attribute{{"gene_id", ""}, {"gene_name", "A"}}},
		{"missing value", `gene_id; gene_name "A";`, []Attribute{{"gene_name", "A"}}},
		{"blank value", `gene_id  ; gene_name "A";`, []Attribute{{"gene_name", "A"}}},
		{"only missing value", `gene_id`, nil},
		{"empty attributes", ``, nil},
		{"unspecified attributes", `.`, nil},
		{"unquoted number", `exon_number 1;`, []Attribute{{"exon_number", "1"}}},
		{"unselected attribute", `other "ignored"; gene_id "g1";`, []Attribute{{"gene_id", "g1"}}},
	} {
		t.Run(test.name, func(t *testing.T) {
			path := filepath.Join(t.TempDir(), "input.gtf")
			input := "chr1\tsrc\tgene\t1\t4\t.\t+\t.\t" + test.attributes + "\n"
			if err := os.WriteFile(path, []byte(input), 0600); err != nil {
				t.Fatal(err)
			}
			features, err := ReadFilteredFeatures(path, nil, nil, []string{"gene_id", "gene_name", "exon_number"})
			if err != nil || len(features) != 1 {
				t.Fatalf("got %d features, err %v; want one feature", len(features), err)
			}
			if got := features[0].Attributes; !slices.Equal(got, test.want) {
				t.Fatalf("attributes = %v, want %v", got, test.want)
			}
		})
	}
}

func TestGTF(t *testing.T) {
	file := "test1.gtf"
	features, err := ReadFeatures(file)
	if err != nil {
		t.Error(err)
		return
	}
	if len(features) != 14 {
		t.Error(err)
		return
	}
	for _, feature := range features {
		fmt.Println(feature)
	}
}

func TestGTF2(t *testing.T) {
	file := "test2.gtf"
	features, err := ReadFeatures(file)
	if err != nil {
		t.Error(err)
		return
	}
	if len(features) != 10 {
		t.Error(err)
		return
	}
	for _, feature := range features {
		fmt.Println(feature)
	}
}

func TestGTF3(t *testing.T) {
	file := "test3.gtf"
	features, err := ReadFeatures(file)
	if err != nil {
		t.Error(err)
		return
	}
	if len(features) != 6 {
		t.Error(err)
		return
	}
	for _, feature := range features {
		fmt.Println(feature)
	}
}
