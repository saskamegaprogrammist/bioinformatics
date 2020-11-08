package comparison

import (
	"bufio"
	"flag"
	"fmt"
	"strconv"
	"strings"
)

func getProteinsType(flags *flag.FlagSet) (int, error) {
	var typeCompare int
	typeString := flags.Lookup("type").Value.String()
	if typeString != "" {
		if typeString == "NUC" {
			typeCompare = NUC
		} else if typeString == "AMINO" {
			typeCompare = AMINO
		} else {
			return typeCompare, fmt.Errorf("wrong type format")
		}
	} else {
		typeCompare = DEFAULT
	}
	return typeCompare, nil
}

func getGapValueSimple(flags *flag.FlagSet) (int, error) {
	var err error
	var gap int
	gapString := flags.Lookup("gap").Value.String()
	if gapString != "0" {
		gap, err = strconv.Atoi(gapString)
		if err != nil {
			return gap, fmt.Errorf("wrong gap format")
		}
	} else {
		gap = DEFAULT_GAP
	}
	return gap, nil
}

func max(valueUp int, valueDiagonal int, valueLeft int) (int, int) {
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

func writeAlignmentStringsToWriter(firstProtein string, secondProtein string, writer *bufio.Writer) error {
	if firstProtein == secondProtein {
		return nil
	}
	_, err := writer.WriteString(firstProtein + "\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}
	_, err = writer.WriteString(secondProtein + "\n\n\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}
	return nil
}

func makeAlignmentStrings(table [][]Cell, firstProtein string, secondProtein string) (string, string, error) {
	var alignmentFirst strings.Builder
	var alignmentSecond strings.Builder
	var err error
	i, j := len(firstProtein), len(secondProtein)
	for j != 0 && i != 0 {
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
			err = fmt.Errorf("wrong from value")
			break
		}
	}
	for i != 0 {
		alignmentFirst.WriteString(string(firstProtein[i-1]))
		alignmentSecond.WriteString("_")
		i--
	}
	for j != 0 {
		alignmentSecond.WriteString(string(secondProtein[j-1]))
		alignmentFirst.WriteString("_")
		j--
	}
	finalStringFirst := reverseString(alignmentFirst.String())
	finalStringSecond := reverseString(alignmentSecond.String())
	return finalStringFirst, finalStringSecond, err
}

func getScore(firstLetter byte, secondLetter byte, comparisonType int) (int, error) {
	var score int
	if firstLetter == secondLetter {
		if comparisonType == NUC {
			score = DNAFULL_MATCH
		} else if comparisonType == AMINO {
			aminoValue, err := CompareAminoLetter(string(firstLetter))
			if err != nil {
				return score, err
			}
			score = aminoValue
		} else {
			score = DEFAULT_MATCH
		}
	} else {
		if comparisonType == NUC {
			score = DNAFULL_MISMATCH
		} else if comparisonType == AMINO {
			aminoValue, err := CompareAminoLetters(string(firstLetter), string(secondLetter))
			if err != nil {
				return score, err
			}
			score = aminoValue
		} else {
			score = DEFAULT_MISMATCH
		}
	}
	return score, nil
}

func NW(firstProtein string, secondProtein string, gap int, comparisonType int) (int, string, string, error) {
	var table [][]Cell
	var firstAligned string
	var secondAligned string
	tableRows := len(firstProtein) + 1
	tableColumns := len(secondProtein) + 1
	err := initializeMatrixCell(&table, tableRows, tableColumns)
	if err != nil {
		return 0, firstAligned, secondAligned, fmt.Errorf("couldn't initialize matrix")
	}
	for i := 0; i < tableRows; i++ {
		table[i][0] = Cell{value: i * gap, from: UP}
	}
	for j := 0; j < tableColumns; j++ {
		table[0][j] = Cell{value: j * gap, from: LEFT}
	}
	for i := 1; i < tableRows; i++ {
		for j := 1; j < tableColumns; j++ {
			score, err := getScore(firstProtein[i-1], secondProtein[j-1], comparisonType)
			if err != nil {
				return 0, firstAligned, secondAligned, err
			}
			max, from := max(table[i-1][j].value+gap, table[i-1][j-1].value+score, table[i][j-1].value+gap)
			table[i][j].from = from
			table[i][j].value = max
		}
	}
	firstAligned, secondAligned, err = makeAlignmentStrings(table, firstProtein, secondProtein)
	if err != nil {
		return table[tableRows-1][tableColumns-1].value, firstAligned, secondAligned, err
	}
	return table[tableRows-1][tableColumns-1].value, firstAligned, secondAligned, nil
}

func compareProteinsNWScore(firstProtein string, secondProtein string, gap int, comparisonType int, writer *bufio.Writer) (int, error) {
	score, firstAligned, secondAligned, err := NW(firstProtein, secondProtein, gap, comparisonType)
	if err != nil {
		return score, err
	}
	err = writeAlignmentStringsToWriter(firstAligned, secondAligned, writer)
	return score, err
}

func compareProteinsNWNoScore(firstProtein string, secondProtein string, gap int, comparisonType int, writer *bufio.Writer) error {
	_, firstAligned, secondAligned, err := NW(firstProtein, secondProtein, gap, comparisonType)
	if err != nil {
		return err
	}
	return writeAlignmentStringsToWriter(firstAligned, secondAligned, writer)
}

func Comparison(flags *flag.FlagSet, proteinStrings []string, writer *bufio.Writer) error {
	gap, err := getGapValueSimple(flags)
	if err != nil {
		return err
	}
	typeCompare, err := getProteinsType(flags)
	if err != nil {
		return err
	}

	var table [][]int
	tableSize := len(proteinStrings)
	err = initializeMatrixInt(&table, tableSize, tableSize)
	if err != nil {
		return fmt.Errorf("couldn't initialize matrix")
	}

	for i := 0; i < tableSize; i++ {
		for j := i; j < tableSize; j++ {
			score, err := compareProteinsNWScore(proteinStrings[i], proteinStrings[j], gap, typeCompare, writer)
			if err != nil {
				return fmt.Errorf("error comparing proteins")
			}
			table[i][j] = score
		}
	}

	for i := 0; i < tableSize; i++ {
		for j := 0; j < tableSize; j++ {
			fmt.Printf("%d ", table[i][j])
		}
		fmt.Println()
	}

	return nil
}

func ComparisonNoScore(flags *flag.FlagSet, proteinStrings []string, writer *bufio.Writer) error {
	gap, err := getGapValueSimple(flags)
	if err != nil {
		return err
	}
	typeCompare, err := getProteinsType(flags)
	if err != nil {
		return err
	}
	proteinsSize := len(proteinStrings)

	_, err = writer.WriteString("alignment with Needleman-Wunsch" + "\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}

	for i := 0; i < proteinsSize; i++ {
		for j := i; j < proteinsSize; j++ {
			err := compareProteinsNWNoScore(proteinStrings[i], proteinStrings[j], gap, typeCompare, writer)
			if err != nil {
				return fmt.Errorf("error comparing proteins")
			}
		}
	}

	_, err = writer.WriteString("alignment with Hirschberg" + "\n")
	if err != nil {
		return fmt.Errorf("error writing string")
	}

	for i := 0; i < proteinsSize; i++ {
		for j := i; j < proteinsSize; j++ {
			err = compareProteinsHirshberg(proteinStrings[i], proteinStrings[j], gap, typeCompare, writer)
			if err != nil {
				return fmt.Errorf("error comparing proteins")
			}
		}
	}

	return nil
}
