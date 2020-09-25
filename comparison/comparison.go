package comparison

import (
	"bufio"
	"flag"
	"fmt"
	"strconv"
	"strings"
)

func getComparisonArguments(flags *flag.FlagSet, gap *int, typeCompare *int) error {
	var err error
	gapString := flags.Lookup("gap").Value.String()
	if gapString != "0" {
		*gap, err = strconv.Atoi(gapString)
		if err != nil {
			return fmt.Errorf("wrong gap format")
		}
	} else {
		*gap = DEFAULT_GAP
	}
	typeString := flags.Lookup("type").Value.String()
	if typeString != "" {
		if typeString == "NUC" {
			*typeCompare = NUC
		} else if  typeString == "AMINO" {
			*typeCompare = AMINO
		} else {
			return fmt.Errorf("wrong type format")
		}
	} else {
		*typeCompare = DEFAULT
	}
	return nil
}


func max(valueUp int, valueDiagonal int, valueLeft int, ) (int, int) {
	max := valueDiagonal
	from := DIAG
	maxLocal := valueUp
	fromLocal := UP
	if valueLeft > maxLocal {
		maxLocal = valueLeft
		fromLocal = LEFT
	}
	if maxLocal > max {
		max = maxLocal
		from = fromLocal
	}
	return max, from
}


func makeAlignmentStrings(table [][]Cell, firstProtein string, secondProtein string, writer *bufio.Writer) error {
	var alignmentFirst strings.Builder
	var alignmentSecond strings.Builder
	for i, j := len(firstProtein), len(secondProtein); j!=0 && i!=0; {
		from := table[i][j].from
		if from == UP {
			alignmentFirst.WriteString(string(firstProtein[i-1]))
			alignmentSecond.WriteString("_")
			i--
		} else if from == LEFT {
			alignmentSecond.WriteString(string(secondProtein[j-1]))
			alignmentFirst.WriteString("_")
			j--
		} else if from == DIAG {
			alignmentSecond.WriteString(string(secondProtein[j-1]))
			alignmentFirst.WriteString(string(firstProtein[i-1]))
			i--
			j--
		} else {
			return fmt.Errorf("wrong from value")
		}
	}
	finalStringFirst := reverseString(alignmentFirst.String())
	finalStringSecond := reverseString(alignmentSecond.String())
	if finalStringFirst == finalStringSecond {
		return nil
	}
	_, err := writer.WriteString( finalStringFirst + "\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}
	_, err = writer.WriteString(finalStringSecond + "\n\n\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}
	return nil
}

func compareProteins(firstProtein string, secondProtein string, gap int, comparisonType int, writer *bufio.Writer) (int, error) {
	var table [][]Cell
	tableRows := len(firstProtein) + 1
	tableColumns := len(secondProtein) + 1
	err := initializeMatrixCell(&table, tableRows, tableColumns)
	if err != nil {
		return 0, fmt.Errorf("couldn't initialize matrix")
	}
	for i:=0; i<tableRows; i++ {
		table[i][0] = Cell{value:i*DEFAULT_GAP, from:UP}
	}
	for j:=0; j<tableColumns; j++ {
		table[0][j] = Cell{value:j*DEFAULT_GAP, from:LEFT}
	}
	for i:=1; i<tableRows; i++ {
		for j:=1; j<tableColumns; j++ {
			diagValue := table[i-1][j-1].value
			if firstProtein[i-1] == secondProtein[j-1] {
				if comparisonType == NUC {
					diagValue += DNAFULL_MATCH
				} else if comparisonType == AMINO {
					aminoValue, err := CompareAminoLetter(string(firstProtein[i-1]))
					if err != nil {
						return 0, err
					}
					diagValue += aminoValue
				} else {
					diagValue += DEFAULT_MATCH
				}
			} else {
				if comparisonType == NUC {
					diagValue += DNAFULL_MISMATCH
				} else if comparisonType == AMINO {
					aminoValue, err := CompareAminoLetters(string(firstProtein[i-1]), string(secondProtein[j-1]))
					if err != nil {
						return 0, err
					}
					diagValue += aminoValue
				} else {
					diagValue += DEFAULT_MISMATCH
				}
			}
			max, from := max(table[i-1][j].value + gap, diagValue, table[i][j-1].value + gap)
			table[i][j].from = from
			table[i][j].value = max
		}
	}
	err = makeAlignmentStrings(table, firstProtein, secondProtein, writer)
	if err != nil {
		return table[tableRows-1][tableColumns-1].value, err
	}
	return table[tableRows-1][tableColumns-1].value, nil
}

func Comparison(flags *flag.FlagSet, proteinStrings []string, writer *bufio.Writer) error {
	var gap int
	var typeCompare int
	err := getComparisonArguments(flags, &gap, &typeCompare)
	if err != nil {
		return err
	}

	var table [][]int
	tableSize := len(proteinStrings)
	err = initializeMatrixInt(&table, tableSize, tableSize)
	if err != nil {
		return fmt.Errorf("couldn't initialize matrix")
	}

	for i:=0; i < tableSize; i++ {
		for j:=i; j<tableSize; j++ {
			score, err := compareProteins(proteinStrings[i], proteinStrings[j], gap, typeCompare, writer)
			if err != nil {
				return fmt.Errorf("error comparing proteins")
			}
			table[i][j] = score
		}
	}

	for i:=0; i < tableSize; i++ {
		for j:=0; j<tableSize; j++ {
			fmt.Printf("%d ", table[i][j])
		}
		fmt.Println()
	}

	return nil
}
