package geom

import (
	"fmt"
	"testing"
)

func TestDeriv(t *testing.T) {
	low, width := [3]int{0, 0, 0}, [3]int{1, 10, 1}
	vals, out := make([]float64, 10), make([]float64, 10)
	for i := range vals { vals[i] = float64(i - 5) * float64(i - 5) }

	g := NewGrid(low, width)

	opt := &DerivOptions{ true, None, 4 }
	g.Deriv(vals, out, 1, 10.0, 10, opt)
	printSlice(vals)
	printSlice(out)
}

func TestCurl(t *testing.T) {
	cells := 11

	width, vecs, out, g := initPaddleWheel(cells)

	op := &DerivOptions{ false, Add, 2 }

	g.Curl(vecs, out, width, cells, op)

	printGrid(out[2], cells)
}

func initPaddleWheel(cells int) (width float64, vecs, out [3][]float64, g *Grid) {
	width = float64(cells)
	for i := 0; i < 3; i++ {
		vecs[i] = make([]float64, cells * cells)
		out[i] = make([]float64, cells * cells)
	}

	g = NewGrid([3]int{0, 0, 0}, [3]int{cells, cells, 1})

	for y := 0; y < cells; y++ {
		for x := 0; x < cells; x++ {
			
			vecs[0][x + y*cells] = float64(+(y - cells/2))
			vecs[1][x + y*cells] = float64(-(x - cells/2))
		}
	}

	return width, vecs, out, g
}

func printGrid(grid []float64, cells int) {
	for y := 0; y < cells; y++ {
		for x := 0; x < cells; x++ {
			fmt.Printf("%g ", grid[x + y * cells])
		}
		fmt.Println()
	}
}

func printSlice(xs []float64) {
	fmt.Print("[ ")
	for _, x := range xs { fmt.Printf("%g ", x) }
	fmt.Println("]")
}
