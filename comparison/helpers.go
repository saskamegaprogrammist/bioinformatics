package comparison

import (
	"bufio"
	"flag"
	"fmt"
	"os"
	"strings"
)

func ReadFile(filename string) ([]string, error) {
	var proteinStrings []string
	_, err := os.Stat(filename)
	if err != nil {
		if os.IsNotExist(err) {
			return proteinStrings, fmt.Errorf("file does not exist")
		}
	}
	file, err := os.Open(filename)
	if err != nil {
		return proteinStrings, fmt.Errorf("error opening file")
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	scanner.Split(bufio.ScanLines)

	var currentProtein strings.Builder

	for scanner.Scan() {
		if scanner.Err() != nil {
			return proteinStrings, fmt.Errorf("error reading file")
		} else {
			text := scanner.Text()
			if text[0] == '>' {
				if len(currentProtein.String()) != 0 {
					proteinStrings = append(proteinStrings, currentProtein.String())
				}
				currentProtein.Reset()
			} else {
				currentProtein.WriteString(text)
			}
		}
	}
	if len(currentProtein.String()) != 0 {
		proteinStrings = append(proteinStrings, currentProtein.String())
	}
	return proteinStrings, nil
}

func createWriter(flags *flag.FlagSet, writer *bufio.Writer, fileOut *os.File) (*bufio.Writer, error) {
	fileOutputName := flags.Lookup("o").Value.String()
	if fileOutputName != flags.Lookup("o").DefValue {
		_, err := os.Stat(flags.Lookup("o").Value.String())
		if err != nil {
			if os.IsNotExist(err) {
				fileOut, err = os.Create(fileOutputName)
				if err != nil {
					return writer, fmt.Errorf("could not create file")
				}
			}
		}

		fileOut, err = os.OpenFile(fileOutputName, os.O_CREATE|os.O_TRUNC|os.O_WRONLY, 0666)
		if err != nil {
			return writer, fmt.Errorf("error opening file")
		}
		writer = bufio.NewWriter(fileOut)
	} else {
		writer = bufio.NewWriter(os.Stdout)

	}
	return writer, nil
}

func argsAssign(flags *flag.FlagSet) {
	flags.Int64("gap", 0, "gap for distance")
	flags.Int64("gap_open", 0, "gap opening for distance")
	flags.Int64("gap_ext", 0, "gap extension for distance")
	flags.Bool("emboss", false, "use emboss algorithm")
	flags.String("o", "", "output file")
	flags.String("i", "", "input file")
	flags.String("type", "", "type of proteins")
}

func ArgsChecker(flags *flag.FlagSet) (string, string, bool, error) {
	argsAssign(flags)
	var filenameInput string
	var filenameOutput string
	var embossBool bool

	err := flags.Parse(os.Args[1:])
	if err != nil {
		return filenameInput, filenameOutput, embossBool, fmt.Errorf("flag input error")
	}
	filenameInput = flags.Lookup("i").Value.String()
	if filenameInput == flags.Lookup("i").DefValue {
		return filenameInput, filenameOutput, embossBool, fmt.Errorf("enter input file name")
	} else {
		if !strings.Contains(filenameInput, ".txt") &&
			!strings.Contains(filenameInput, ".fasta") {
			return filenameInput, filenameOutput, embossBool, fmt.Errorf("wrong input file format")
		}
	}
	filenameOutput = flags.Lookup("o").Value.String()
	if filenameOutput != flags.Lookup("o").DefValue &&
		!strings.Contains(filenameOutput, ".txt") {
		return filenameInput, filenameOutput, embossBool, fmt.Errorf("wrong output file format")
	}
	if flags.Lookup("emboss").Value.String() == "true" {
		embossBool = true
	}
	return filenameInput, filenameOutput, embossBool, nil
}

func Comparing(flags *flag.FlagSet, proteinStrings []string, score bool, emboss bool) error {
	var writer *bufio.Writer
	var fileOut *os.File
	defer fileOut.Close()
	writer, err := createWriter(flags, writer, fileOut)
	if err != nil {
		return err
	}
	//for _, str := range proteinStrings{
	//	fmt.Println(str)
	//}
	if score {
		if emboss {
			err = ComparisonEmboss(flags, proteinStrings, writer)
			if err != nil {
				return err
			}
		} else {
			err = Comparison(flags, proteinStrings, writer)
			if err != nil {
				return err
			}
		}

	} else {
		err = ComparisonNoScore(flags, proteinStrings, writer)
		if err != nil {
			return err
		}
	}

	err = writer.Flush()
	if err != nil {
		return err
	}
	return nil
}
