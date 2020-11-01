package main

import (
	"flag"
	"fmt"
	"github.com/saskamegaprogrammist/bioinformatics/comparison"
)



func main() {
	var flags flag.FlagSet
	filenameInput, _, err := comparison.ArgsChecker(&flags)
	if err != nil {
		fmt.Println(err)
		return
	}
	proteinStrings, err := comparison.ReadFile(filenameInput)
	if err != nil {
		fmt.Println(err)
		return
	}
	err = comparison.Comparing(&flags, proteinStrings, false)
	if err != nil {
		fmt.Println(err)
	}
}
