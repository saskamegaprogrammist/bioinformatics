package main

import (
	"flag"
	"io/ioutil"
	"log"
	"os"
	"testing"
)

func Test_Compare_Simple_1(t *testing.T) {
	os.Args[1] = `-i=testFiles/prot_1.txt`
	os.Args[2] = `-o=testFiles/out_1.txt`
	var flags flag.FlagSet
	err := argsChecker(&flags)
	if err != nil {
		t.Errorf("Test_Compare_Simple_1 failed: %s", err)
	} else {
		proteinStrings, err := readFile(filenameInput)
		if err != nil {
			t.Errorf("Test_Compare_Simple_1 failed: %s", err)
		} else {
			err = comparing(&flags, proteinStrings)
			if err != nil {
				t.Errorf("Test_Compare_Simple_1 failed: %s", err)
			} else {
				outData, err := ioutil.ReadFile(filenameOutput)
				if err != nil {
					log.Fatal(err)
				}
				rightData, err := ioutil.ReadFile("testFiles/answer_1.txt")
				if err != nil {
					log.Fatal(err)
				}
				if string(outData) != string(rightData) {
					t.Errorf("Test_Compare_Simple_1 failed, result not match")
				}
			}
		}
	}
}

func Test_Compare_Simple_2(t *testing.T) {
	os.Args[1] = `-i=testFiles/prot_2.txt`
	os.Args[2] = `-o=testFiles/out_2.txt`
	var flags flag.FlagSet
	err := argsChecker(&flags)
	if err != nil {
		t.Errorf("Test_Compare_Simple_2 failed: %s", err)
	} else {
		proteinStrings, err := readFile(filenameInput)
		if err != nil {
			t.Errorf("Test_Compare_Simple_2 failed: %s", err)
		} else {
			err = comparing(&flags, proteinStrings)
			if err != nil {
				t.Errorf("Test_Compare_Simple_2 failed: %s", err)
			} else {
				outData, err := ioutil.ReadFile("testFiles/out_2.txt")
				if err != nil {
					log.Fatal(err)
				}
				rightData, err := ioutil.ReadFile("testFiles/answer_2.txt")
				if err != nil {
					log.Fatal(err)
				}
				if string(outData) != string(rightData) {
					t.Errorf("Test_Compare_Simple_2 failed, result not match")
				}
			}
		}
	}
}

func Test_Compare_Simple_3(t *testing.T) {
	os.Args[1] = `-i=testFiles/dna_1.txt`
	os.Args[2] = `-o=testFiles/out_3.txt`
	var flags flag.FlagSet
	err := argsChecker(&flags)
	if err != nil {
		t.Errorf("Test_Compare_Simple_3 failed: %s", err)
	} else {
		proteinStrings, err := readFile(filenameInput)
		if err != nil {
			t.Errorf("Test_Compare_Simple_3 failed: %s", err)
		} else {
			err = comparing(&flags, proteinStrings)
			if err != nil {
				t.Errorf("Test_Compare_Simple_3 failed: %s", err)
			} else {
				outData, err := ioutil.ReadFile(filenameOutput)
				if err != nil {
					log.Fatal(err)
				}
				rightData, err := ioutil.ReadFile("testFiles/answer_3.txt")
				if err != nil {
					log.Fatal(err)
				}
				if string(outData) != string(rightData) {
					t.Errorf("Test_Compare_Simple_3 failed, result not match")
				}
			}
		}
	}
}

func Test_Compare_Blosum_2(t *testing.T) {
	os.Args[1] = `-i=testFiles/prot_2.txt`
	os.Args[2] = `-o=testFiles/out_2_blosum.txt`
	os.Args[3] = `-type=AMINO`
	var flags flag.FlagSet
	err := argsChecker(&flags)
	if err != nil {
		t.Errorf("Test_Compare_Blosum_2 failed: %s", err)
	} else {
		proteinStrings, err := readFile(filenameInput)
		if err != nil {
			t.Errorf("Test_Compare_Blosum_2 failed: %s", err)
		} else {
			err = comparing(&flags, proteinStrings)
			if err != nil {
				t.Errorf("Test_Compare_Blosum_2 failed: %s", err)
			} else {
				outData, err := ioutil.ReadFile(filenameOutput)
				if err != nil {
					log.Fatal(err)
				}
				rightData, err := ioutil.ReadFile("testFiles/answer_2_blosum.txt")
				if err != nil {
					log.Fatal(err)
				}
				if string(outData) != string(rightData) {
					t.Errorf("Test_Compare_Blosum_2 failed, result not match")
				}
			}
		}
	}
}

func Test_Compare_DNA_3(t *testing.T) {
	os.Args[1] = `-i=testFiles/dna_1.txt`
	os.Args[2] = `-o=testFiles/out_3_dna.txt`
	os.Args[3] = `-type=NUC`
	var flags flag.FlagSet
	err := argsChecker(&flags)
	if err != nil {
		t.Errorf("Test_Compare_Simple_3 failed: %s", err)
	} else {
		proteinStrings, err := readFile(filenameInput)
		if err != nil {
			t.Errorf("Test_Compare_Simple_3 failed: %s", err)
		} else {
			err = comparing(&flags, proteinStrings)
			if err != nil {
				t.Errorf("Test_Compare_Simple_3 failed: %s", err)
			} else {
				outData, err := ioutil.ReadFile(filenameOutput)
				if err != nil {
					log.Fatal(err)
				}
				rightData, err := ioutil.ReadFile("testFiles/answer_3_dna.txt")
				if err != nil {
					log.Fatal(err)
				}
				if string(outData) != string(rightData) {
					t.Errorf("Test_Compare_Simple_3 failed, result not match")
				}
			}
		}
	}
}
