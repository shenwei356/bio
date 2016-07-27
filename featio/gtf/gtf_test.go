package gtf

import (
	"fmt"
	"testing"
)

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
