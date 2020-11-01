package comparison

import "fmt"

func initializeMatrixInt(array *[][]int, rows int, columns int) error {
	if !(rows > 0 && columns > 0) {
		return fmt.Errorf("invalid matrix sizes")
	}
	*array = make([][]int, 0)

	for i := 0; i < rows; i++ {
		inner := make([]int, 0)
		for j := 0; j < columns; j++ {
			inner = append(inner, 0)
		}
		*array = append(*array, inner)
	}
	return nil
}

func initializeMatrixCell(array *[][]Cell, rows int, columns int) error {
	if !(rows > 0 && columns > 0) {
		return fmt.Errorf("invalid matrix sizes")
	}
	*array = make([][]Cell, 0)

	for i := 0; i < rows; i++ {
		inner := make([]Cell, 0)
		for j := 0; j < columns; j++ {
			inner = append(inner, Cell{value: 0, from: 0})
		}
		*array = append(*array, inner)
	}
	return nil
}

func reverseString(alignment string) string {
	n := 0
	stringToRune := make([]rune, len(alignment))
	for _, r := range alignment {
		stringToRune[n] = r
		n++
	}
	stringToRune = stringToRune[0:n]

	for i := 0; i < n/2; i++ {
		stringToRune[i], stringToRune[n-1-i] = stringToRune[n-1-i], stringToRune[i]
	}

	return string(stringToRune)
}
