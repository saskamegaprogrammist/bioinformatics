package comparison

import (
	"bufio"
	"fmt"
	"strings"
)

func optimisedScore(firstProtein string, secondProtein string, gap int, comparisonType int) ([]int, error) {
	var table [][]int
	tableRows := 2
	tableColumns := len(secondProtein) + 1
	err := initializeMatrixInt(&table, tableRows, tableColumns)
	if err != nil {
		return table[1], fmt.Errorf("couldn't initialize matrix")
	}
	for j := 0; j < tableColumns; j++ {
		table[0][j] = j * gap
	}
	for i := 1; i < len(firstProtein)+1; i++ {
		table[1][0] = table[0][0] + gap
		for j := 1; j < tableColumns; j++ {
			score, err := getScore(firstProtein[i-1], secondProtein[j-1], comparisonType)
			if err != nil {
				return table[1], err
			}
			max, _ := max(table[0][j]+gap, table[0][j-1]+score, table[1][j-1]+gap)
			table[1][j] = max
		}
		for j := 0; j < tableColumns; j++ {
			table[0][j] = table[1][j]
		}
	}
	return table[1], nil
}

func getSplitIndex(firstScore []int, secondScore []int) int {
	length := len(firstScore)
	max := firstScore[0] + secondScore[length-1]
	var maxIndex int
	for i := 1; i < length; i++ {
		scoreLocal := firstScore[i] + secondScore[length-i-1]
		if scoreLocal > max {
			max = scoreLocal
			maxIndex = i
		}
	}
	return maxIndex
}

func Hirshberg(firstProtein string, secondProtein string, gap int, comparisonType int) (string, string) {
	var alignmentFirst strings.Builder
	var alignmentSecond strings.Builder
	var firstAligned string
	var secondAligned string
	firstLength := len(firstProtein)
	secondLength := len(secondProtein)
	if firstLength == 0 {
		for i := 0; i < secondLength; i++ {
			alignmentFirst.WriteString("_")
			alignmentSecond.WriteString(string(secondProtein[i]))
		}
		return alignmentFirst.String(), alignmentSecond.String()
	} else if secondLength == 0 {
		for i := 0; i < firstLength; i++ {
			alignmentFirst.WriteString(string(firstProtein[i]))
			alignmentSecond.WriteString("_")
		}
		return alignmentFirst.String(), alignmentSecond.String()
	} else if len(firstProtein) == 1 || len(secondProtein) == 1 {
		_, firstAligned, secondAligned, err := NW(firstProtein, secondProtein, gap, comparisonType)
		if err != nil {
			fmt.Println(err)
		}
		return firstAligned, secondAligned
	} else {
		firstMiddle := firstLength / 2
		firstHalfFirst := firstProtein[0:firstMiddle]
		secondHalfFirst := firstProtein[firstMiddle:firstLength]
		scoreFirst, err := optimisedScore(firstHalfFirst, secondProtein, gap, comparisonType)
		if err != nil {
			fmt.Println(err)
			return firstAligned, secondAligned
		}
		scoreSecond, err := optimisedScore(reverseString(secondHalfFirst), reverseString(secondProtein), gap, comparisonType)
		if err != nil {
			fmt.Println(err)
			return firstAligned, secondAligned
		}
		index := getSplitIndex(scoreFirst, scoreSecond)
		firstHalfSecond := secondProtein[0:index]
		secondHalfSecond := secondProtein[index:secondLength]
		firstSplitFirst, firstSplitSecond := Hirshberg(firstHalfFirst, firstHalfSecond, gap, comparisonType)
		secondSplitFirst, secondSplitSecond := Hirshberg(secondHalfFirst, secondHalfSecond, gap, comparisonType)
		return firstSplitFirst + secondSplitFirst, firstSplitSecond + secondSplitSecond
	}
}

func compareProteinsHirshberg(firstProtein string, secondProtein string, gap int, comparisonType int, writer *bufio.Writer) error {
	firstAligned, secondAligned := Hirshberg(firstProtein, secondProtein, gap, comparisonType)
	err := writeAlignmentStringsToWriter(firstAligned, secondAligned, writer)
	return err
}
