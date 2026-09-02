package seq

import (
	"fmt"
	"math"
	"testing"
)

func TestValidateSequence(t *testing.T) {
	dna := []byte("acgt")
	dna2 := []byte("ACGTRYMKSWHBVDN")
	rna := []byte("acgu")
	fake := []byte("acgturymkswhbvdnz")

	ok := DNA.IsValid(dna) == nil &&
		DNAredundant.IsValid(dna2) == nil &&
		RNA.IsValid(rna) == nil &&
		RNA.IsValid(fake) != nil

	if !ok {
		t.Error("validate sequence failed.")
		return
	}
}

func TestRevCom(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtccn-"))
	if string(dna.RevCom().Seq) != "-nggacgt" {
		t.Error("revcom sequence failed.")
		return
	}

	rna, _ := NewSeq(RNA, []byte("auguccn-"))
	if string(rna.RevCom().Seq) != "-nggacau" {
		t.Error("revcom sequence failed.")
		return
	}
}

func TestBaseContent(t *testing.T) {
	dna, _ := NewSeq(DNA, []byte("acgtACGT"))
	content := dna.BaseContent("gc")
	wanted := 0.5
	if content != wanted {
		t.Error(fmt.Printf("compution of base content failed: %f != %f", content, wanted))
		return
	}
}

func TestRemoveGapsInplaceReusesBackingArrays(t *testing.T) {
	bases := []byte("A-C.G-T")
	quality := []byte("1234567")
	s, err := NewSeqWithQualWithoutValidation(DNA, bases, quality)
	if err != nil {
		t.Fatal(err)
	}

	basePtr := &bases[0]
	qualityPtr := &quality[0]
	s.RemoveGapsInplace("-.")

	if got, want := string(s.Seq), "ACGT"; got != want {
		t.Fatalf("unexpected sequence after removing gaps: %q != %q", got, want)
	}
	if got, want := string(s.Qual), "1357"; got != want {
		t.Fatalf("unexpected quality after removing gaps: %q != %q", got, want)
	}
	if &s.Seq[0] != basePtr || &s.Qual[0] != qualityPtr {
		t.Fatal("RemoveGapsInplace replaced a backing array")
	}

	if allocs := testing.AllocsPerRun(1000, func() {
		s.RemoveGapsInplace("-.")
	}); allocs != 0 {
		t.Fatalf("RemoveGapsInplace allocated %.1f times per call", allocs)
	}
}

func TestAvgQualDoesNotPopulateQualValue(t *testing.T) {
	s := &Seq{Qual: []byte("IIII")}
	if got, want := s.AvgQual(33), 40.0; math.Abs(got-want) > 1e-12 {
		t.Fatalf("unexpected average quality: %g != %g", got, want)
	}
	if s.QualValue != nil {
		t.Fatal("AvgQual unexpectedly populated QualValue")
	}

	s = &Seq{QualValue: []int{40, 40, 40, 40}}
	if got, want := s.AvgQual(33), 40.0; math.Abs(got-want) > 1e-12 {
		t.Fatalf("unexpected average parsed quality: %g != %g", got, want)
	}

	s = &Seq{Qual: []byte("!!!!"), QualValue: []int{40, 40, 40, 40}}
	if got, want := s.AvgQual(33), 40.0; math.Abs(got-want) > 1e-12 {
		t.Fatalf("AvgQual did not prefer QualValue: %g != %g", got, want)
	}
}

func TestAvgQualOfFASTA(t *testing.T) {
	s, err := NewSeq(DNA, []byte("ACGT"))
	if err != nil {
		t.Fatal(err)
	}

	if got := s.AvgQual(33); got != 0 {
		t.Errorf("AvgQual returned %g instead of 0", got)
	}
	if got := s.AvgQualOfRegion(33, 1, -1); got != 0 {
		t.Errorf("AvgQualOfRegion returned %g instead of 0", got)
	}
}

func TestAvgQualOfRegion(t *testing.T) {
	qualityValue := []int{10, 15, 20, 25, 30, 35, 40, 45, 50, 55}
	quality := make([]byte, len(qualityValue))
	for i, q := range qualityValue {
		quality[i] = byte(q + 33)
	}

	avg := func(values []int) float64 {
		var sum float64
		for _, q := range values {
			sum += math.Pow(10, float64(q)/-10)
		}
		return -10 * math.Log10(sum/float64(len(values)))
	}

	tests := []struct {
		name       string
		start, end int
		values     []int
	}{
		{"first base", 1, 1, qualityValue[0:1]},
		{"positive coordinates", 2, 4, qualityValue[1:4]},
		{"negative coordinates", -4, -2, qualityValue[6:9]},
		{"negative end", 2, -2, qualityValue[1:9]},
		{"clamped end", 3, 100, qualityValue[2:10]},
		{"clamped start", -12, -1, qualityValue[0:10]},
		{"zero coordinates", 0, 0, qualityValue[0:10]},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			want := avg(tt.values)
			for name, s := range map[string]*Seq{
				"QualValue": {QualValue: qualityValue},
				"Qual":      {Qual: quality},
			} {
				if got := s.AvgQualOfRegion(33, tt.start, tt.end); math.Abs(got-want) > 1e-12 {
					t.Errorf("%s: unexpected average quality: %g != %g", name, got, want)
				}
			}
		})
	}

	invalidRegions := [][2]int{{12, 14}, {-4, 2}, {-3, -4}, {2, 1}}
	for _, region := range invalidRegions {
		if got := (&Seq{QualValue: qualityValue}).AvgQualOfRegion(33, region[0], region[1]); got != 0 {
			t.Errorf("invalid region %v returned %g instead of 0", region, got)
		}
	}
	if got := (&Seq{}).AvgQualOfRegion(33, 1, -1); got != 0 {
		t.Errorf("empty quality returned %g instead of 0", got)
	}

	s := &Seq{Qual: []byte("!!!!"), QualValue: []int{40, 40, 40, 40}}
	if got, want := s.AvgQualOfRegion(33, 2, 3), 40.0; math.Abs(got-want) > 1e-12 {
		t.Fatalf("AvgQualOfRegion did not prefer QualValue: %g != %g", got, want)
	}
}

func BenchmarkAvgQual(b *testing.B) {
	quality := make([]byte, 1<<20)
	qualityValue := make([]int, len(quality))
	for i := range quality {
		quality[i] = byte(33 + i%42)
		qualityValue[i] = i % 42
	}
	b.Run("QualValue", func(b *testing.B) {
		s := &Seq{Qual: quality, QualValue: qualityValue}
		b.SetBytes(int64(len(quality)))
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			_ = s.AvgQual(33)
		}
	})
	b.Run("Qual", func(b *testing.B) {
		s := &Seq{Qual: quality}
		b.SetBytes(int64(len(quality)))
		b.ReportAllocs()
		for i := 0; i < b.N; i++ {
			_ = s.AvgQual(33)
		}
	})
}

func TestSubSeq(t *testing.T) {
	s, _ := NewSeqWithoutValidation(DNA, []byte("ACGTNacgtn"))
	ok := string(s.SubSeq(1, 1).Seq) == "A" &&
		string(s.SubSeq(2, 4).Seq) == "CGT" &&
		string(s.SubSeq(-4, -2).Seq) == "cgt" &&
		string(s.SubSeq(-4, -1).Seq) == "cgtn" &&
		string(s.SubSeq(-1, -1).Seq) == "n" &&
		string(s.SubSeq(2, -2).Seq) == "CGTNacgt" &&
		string(s.SubSeq(1, -1).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(12, 14).Seq) == "" &&
		string(s.SubSeq(-10, -1).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(-10, -3).Seq) == "ACGTNacg" &&
		string(s.SubSeq(1, 10).Seq) == "ACGTNacgtn" &&
		string(s.SubSeq(3, 12).Seq) == "GTNacgtn" &&
		string(s.SubSeq(3, 100).Seq) == "GTNacgtn"

	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	ok = string(s.SubSeq(-4, 2).Seq) == "" &&
		string(s.SubSeq(-3, -4).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	s, _ = NewSeqWithoutValidation(DNA, []byte(""))
	ok = string(s.SubSeq(1, 4).Seq) == "" &&
		string(s.SubSeq(2, 4).Seq) == "" &&
		string(s.SubSeq(1, -1).Seq) == "" &&
		string(s.SubSeq(-4, -1).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}

	s, _ = NewSeqWithoutValidation(DNA, []byte("ACGT"))
	s.QualValue = []int{10, 20, 30, 40}
	sub := s.SubSeq(2, 3)
	if got, want := fmt.Sprint(sub.QualValue), "[20 30]"; got != want {
		t.Fatalf("unexpected quality values: %s != %s", got, want)
	}
	sub.QualValue[0] = 0
	if s.QualValue[1] != 20 {
		t.Fatal("SubSeq shares its quality values with the original sequence")
	}
}

func TestSubSeqInplace(t *testing.T) {
	s, _ := NewSeqWithoutValidation(DNA, []byte("ACGTNacgtn"))
	ok := string(s.Clone().SubSeqInplace(1, 1).Seq) == "A" &&
		string(s.Clone().SubSeqInplace(2, 4).Seq) == "CGT" &&
		string(s.Clone().SubSeqInplace(-4, -2).Seq) == "cgt" &&
		string(s.Clone().SubSeqInplace(-4, -1).Seq) == "cgtn" &&
		string(s.Clone().SubSeqInplace(-1, -1).Seq) == "n" &&
		string(s.Clone().SubSeqInplace(2, -2).Seq) == "CGTNacgt" &&
		string(s.Clone().SubSeqInplace(1, -1).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(-10, -1).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(-10, -3).Seq) == "ACGTNacg" &&
		string(s.Clone().SubSeqInplace(1, 10).Seq) == "ACGTNacgtn" &&
		string(s.Clone().SubSeqInplace(3, 10).Seq) == "GTNacgtn" &&
		string(s.Clone().SubSeqInplace(3, 100).Seq) == "GTNacgtn"

	if !ok {
		t.Error(fmt.Printf("SubSeqInplace error"))
	}

	ok = string(s.Clone().SubSeqInplace(-4, 2).Seq) == "" &&
		string(s.Clone().SubSeqInplace(-3, -4).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("SubSeqInplace error"))
	}

	s, _ = NewSeqWithoutValidation(DNA, []byte(""))
	ok = string(s.Clone().SubSeqInplace(1, 4).Seq) == "" &&
		string(s.Clone().SubSeqInplace(2, 4).Seq) == "" &&
		string(s.Clone().SubSeqInplace(1, -1).Seq) == "" &&
		string(s.Clone().SubSeqInplace(-4, -1).Seq) == ""
	if !ok {
		t.Error(fmt.Printf("subseq error"))
	}
}
