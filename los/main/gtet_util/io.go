package gtet_util

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func stdinLines() ([]string, error) {
	bs, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		return nil, fmt.Errorf(
			"Error reading stdin: %s.", err.Error(),
		)
	}

	text := string(bs)
	return strings.Split(text, "\n"), nil
}

// ParseStdin parses a stdin which consists of a series of columns. The first
// two columns must be included in all non=empty rows and must be ints. The
// first will be parsed at halo IDs and the second as halo snapshots. Aside
// from these two columns, each line may have any number of additional columns,
// which will be parsed as floats and returned line-by-line as inVals.
//
// An error is returned either on I/O errors or on parsing errors.
func ParseStdin() (ids, snaps []int, inVals [][]float64, err error) {
	ids, snaps, inVals = []int{}, []int{}, [][]float64{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, nil, err }
	for i, line := range lines {
		rawTokens := strings.Split(line, " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}

		var (
			id, snap int
			vals []float64
			err error
		)
		switch {
		case len(tokens) == 0:
			continue
		case len(tokens) == 1:
			if tokens[0] == "" { continue }
			return nil, nil, nil, fmt.Errorf(
				"Line %d of stdin has 1 token, but " +
					"at least 2 are required.", i + 1,
			)
		case len(tokens) >= 2:
			id, err = strconv.Atoi(tokens[0])
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[0],
				)
			} 
			snap, err = strconv.Atoi(tokens[1]) 
			if err != nil {
				return nil, nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[1],
				)
			}
			
			vals = make([]float64, len(tokens) - 2) 
			for i := range vals {
				vals[i], err = strconv.ParseFloat(tokens[i + 2], 64)
			}
		}

		ids = append(ids, id)
		snaps = append(snaps, snap)
		inVals = append(inVals, vals)
	}

	return ids, snaps, inVals, nil
}

// PrintRows takes a seqeunce of variable-length rows, each associated with
// a halo ID and snapshot, and prints them out as uniformly-spaced columns
// to stdout.
//
// PrintRows panics if the length of ids, snaps, and rows are not all the same.
func PrintRows(ids, snaps []int, rows [][]float64) {
	height := len(ids)
	if height != len(snaps) {
		panic("Height of ID column does not equal height of snapshot column.")
	} else if height != len(rows) {
		panic("Height of rows does not equal height of ID column.")
	}

	// Maximum number of items per row.
	maxWidth := 0
	for i := range rows {
		if len(rows[i]) > maxWidth { maxWidth = len(rows[i]) }
	}

	idW, snapW := 0, 0
	ws := make([]int, maxWidth)

	// Maximum character width needed within each row.
	for i := 0; i < height; i++ {
		idN := len(fmt.Sprintf("%d", ids[i]))
		if idN > idW { idW = idN }
		snapN := len(fmt.Sprintf("%d", snaps[i]))
		if snapN > snapW { snapW = snapN }
		for j := range rows[i] {
			valN := len(fmt.Sprintf("%.10g", rows[i][j]))
			if valN > ws[j] { ws[j] = valN }
		}
	}

	intFmt := fmt.Sprintf("%%%dd %%%dd", idW, snapW)
	floatFmts := make([]string, maxWidth)
	for i := 0; i < len(ws); i++ {
		floatFmts[i] = fmt.Sprintf(" %%%d.10g", ws[i])
	}

	for i := 0; i < height; i++ {
		fmt.Printf(intFmt, ids[i], snaps[i])
		for j := range rows[i] { fmt.Printf(floatFmts[j], rows[i][j]) }
		fmt.Println()
	}
}

// PrintCols prints out a sequence of columns to stdin such that each column
// has the same character width.
//
// PrintCols panics if all the columns are not the same height.
func PrintCols(ids, snaps []int, cols ...[]float64) {
	height := len(ids)
	width := len(cols)

	// Start with some consistency checks.

	if len(ids) != len(snaps) {
		panic("Height of ID column does not equal height of snapshot column.")
	}

	for i := range cols {
		if len(cols[i]) != height {
			panic("All printed columns must be the same height.")
		}
	}

	// Transpose and print.

	rows := make([][]float64, height)
	for i := range rows { rows[i] = make([]float64, width) }

	for y := 0; y < height; y++ {
		for x := 0; x < width; x++ {
			rows[y][x] = cols[x][y]
		}
	}

	PrintRows(ids, snaps, rows)
}
