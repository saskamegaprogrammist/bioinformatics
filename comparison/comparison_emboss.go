package comparison

import (
	"bufio"
	"flag"
	"fmt"
	"strconv"
	"strings"
)

func getGapValueEmboss(flags *flag.FlagSet) (int, int, error) {
	var err error
	var gapOpen int
	var gapExtension int
	gapOpenString := flags.Lookup("gap_open").Value.String()
	if gapOpenString != "0" {
		gapOpen, err = strconv.Atoi(gapOpenString)
		if err != nil {
			return gapOpen, gapExtension, fmt.Errorf("wrong gap open format")
		}
	} else {
		gapOpen = DEFAULT_GAP_OPEN
	}
	gapExtString := flags.Lookup("gap_ext").Value.String()
	if gapExtString != "0" {
		gapExtension, err = strconv.Atoi(gapExtString)
		if err != nil {
			return gapOpen, gapExtension, fmt.Errorf("wrong gap extension format")
		}
	} else {
		gapExtension = DEFAULT_GAP_EXTEND
	}
	return gapOpen, gapExtension, nil
}

func maxEmboss(valueM int, valueI int, valueD int, last int) (int, int) {
	max := valueM
	from := M
	maxLocal := valueI
	fromLocal := I
	if valueD > maxLocal {
		maxLocal = valueD
		fromLocal = D
	}
	if last == D && valueD == valueI {
		maxLocal = valueD
		fromLocal = D
	}
	if maxLocal > max {
		max = maxLocal
		from = fromLocal
	}
	if last != M && maxLocal == max {
		max = maxLocal
		from = fromLocal
	}
	return max, from
}

func maxSimple(valueM int, valueI int, valueD int) int {
	max := valueM
	maxLocal := valueI
	if valueD >= maxLocal {
		maxLocal = valueD
	}
	if maxLocal >= max {
		max = maxLocal
	}
	return max
}

func makeAlignmentStringsEmboss(table [][]CellEmboss, firstProtein string, secondProtein string) (string, string, error) {
	var alignmentFirst strings.Builder
	var alignmentSecond strings.Builder
	var err error
	var last = M
	i, j := len(firstProtein), len(secondProtein)
	for j != 0 && i != 0 {
		_, from := maxEmboss(table[i][j].valueM, table[i][j].valueI, table[i][j].valueD, last)
		last = from
		if from == D {
			alignmentFirst.WriteString(string(firstProtein[i-1]))
			alignmentSecond.WriteString("_")
			i--
		} else if from == I {
			alignmentSecond.WriteString(string(secondProtein[j-1]))
			alignmentFirst.WriteString("_")
			j--
		} else if from == M {
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

func NWEmboss(firstProtein string, secondProtein string, gapOpen int, gapExtension int, comparisonType int) (int, string, string, error) {
	var table [][]CellEmboss
	var firstAligned string
	var secondAligned string
	tableRows := len(firstProtein) + 1
	tableColumns := len(secondProtein) + 1
	INFINITY := 2*gapOpen + (tableColumns+tableRows)*gapExtension - 1
	err := initializeMatrixCellEmboss(&table, tableRows, tableColumns)
	if err != nil {
		return 0, firstAligned, secondAligned, fmt.Errorf("couldn't initialize matrix")
	}
	table[0][0] = CellEmboss{
		valueM: 0,
		valueI: INFINITY,
		valueD: INFINITY,
	}
	for i := 1; i < tableRows; i++ {
		table[i][0] = CellEmboss{
			valueM: INFINITY,
			valueI: INFINITY,
			valueD: gapOpen + (i-1)*gapExtension,
		}
	}
	for j := 1; j < tableColumns; j++ {
		table[0][j] = CellEmboss{
			valueM: INFINITY,
			valueI: gapOpen + (j-1)*gapExtension,
			valueD: INFINITY,
		}
	}
	for i := 1; i < tableRows; i++ {
		for j := 1; j < tableColumns; j++ {
			score, err := getScore(firstProtein[i-1], secondProtein[j-1], comparisonType)
			if err != nil {
				return 0, firstAligned, secondAligned, err
			}
			maxM := maxSimple(table[i-1][j-1].valueM+score, table[i-1][j-1].valueI+score, table[i-1][j-1].valueD+score)
			table[i][j].valueM = maxM
			maxI := maxSimple(table[i][j-1].valueM+gapOpen, table[i][j-1].valueI+gapExtension, table[i][j-1].valueD+gapOpen)
			table[i][j].valueI = maxI
			maxD := maxSimple(table[i-1][j].valueM+gapOpen, table[i-1][j].valueI+gapOpen, table[i-1][j].valueD+gapExtension)
			table[i][j].valueD = maxD
		}
	}

	finalScore := maxSimple(table[tableRows-1][tableColumns-1].valueM, table[tableRows-1][tableColumns-1].valueI, table[tableRows-1][tableColumns-1].valueD)

	firstAligned, secondAligned, err = makeAlignmentStringsEmboss(table, firstProtein, secondProtein)
	if err != nil {
		return finalScore, firstAligned, secondAligned, err
	}
	return finalScore, firstAligned, secondAligned, nil
}

func compareProteinsNWEmbossScore(firstProtein string, secondProtein string, gapOpen int, gapExtension int, comparisonType int, writer *bufio.Writer) (int, error) {
	score, firstAligned, secondAligned, err := NWEmboss(firstProtein, secondProtein, gapOpen, gapExtension, comparisonType)
	if err != nil {
		return score, err
	}
	err = writeAlignmentStringsToWriter(firstAligned, secondAligned, writer)
	return score, err
}

func ComparisonEmboss(flags *flag.FlagSet, proteinStrings []string, writer *bufio.Writer) error {
	gapOpen, gapExtend, err := getGapValueEmboss(flags)
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
			score, err := compareProteinsNWEmbossScore(proteinStrings[i], proteinStrings[j], gapOpen, gapExtend, typeCompare, writer)
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
