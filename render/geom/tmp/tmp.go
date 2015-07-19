package main

import (
	"fmt"
	"github.com/phil-mansfield/gotetra/geom"
	"github.com/gonum/matrix/mat64"
)

func main() {
	xs := make([]float64, 16)
	
	c1 := &geom.Vec{0, 0, 0}
	c2 := &geom.Vec{1, 0, 0}
	c3 := &geom.Vec{0, 1, 0}
	c4 := &geom.Vec{0, 0, 1}

	tet, _ := geom.NewTetra(c1, c2, c3, c4, 1e6)
	J := FromTetra(tet, xs)
	v0, v1, v2, v3 := 1.0, 1.0, 1.0, 1.0
	fmt.Println(
		DerivX(J, v0, v1, v2, v3),
		DerivY(J, v0, v1, v2, v3),
		DerivZ(J, v0, v1, v2, v3),
	)
}

func FromTetra(t *geom.Tetra, xs []float64) *mat64.Dense {
	for i := 0; i < 4; i++ {
		xs[i] = 1.0
		xs[4 + i] =  float64(t.Corners[i][0])
		xs[8 + i] =  float64(t.Corners[i][1])
		xs[12 + i] = float64(t.Corners[i][2])
	}

	JInv, _ := mat64.Inverse(mat64.NewDense(4, 4, xs))
	return JInv
}

func Deriv(inv *mat64.Dense, dir int, val0, val1, val2, val3 float64) float64 {
	d := dir + 1
	return inv.At(0, d) * val0 + inv.At(1, d) * val1 +
		inv.At(2, d) * val2 + inv.At(3, d) * val3
}

func DerivX(J *mat64.Dense, val0, val1, val2, val3 float64) float64 {
	return Deriv(J, 0, val0, val1, val2, val3)
}

func DerivY(J *mat64.Dense, val0, val1, val2, val3 float64) float64 {
	return Deriv(J, 1, val0, val1, val2, val3)
}

func DerivZ(J *mat64.Dense, val0, val1, val2, val3 float64) float64 {
	return Deriv(J, 2, val0, val1, val2, val3)
}

func PrintMat(mat *mat64.Dense) {
	rows, cols := mat.Dims()
	for r := 0; r < rows; r++ {
		fmt.Print("[")
		for c := 0; c < cols; c++ {
			fmt.Print(mat.At(r, c))
			if c != cols - 1 { fmt.Print(" ") }
		}
		fmt.Println("]")
	}
}
