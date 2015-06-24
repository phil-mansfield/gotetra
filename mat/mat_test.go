package mat

import (
	"github.com/gonum/matrix/mat64"
	"fmt"
	"testing"
)

func TestLU(t *testing.T) {
	M := NewMatrix([]float64{
		1, 3, 5,
		2, 4, 7,
		1, 1, 0,
	}, 3, 3)
	
	//M = NewMatrix([]float64{
	//	2, 1,
	//	-1, 0,
	//}, 2, 2)

	fmt.Println("M", M.Vals)
	luf := M.LU()
	fmt.Println("LU", luf.lu.Vals)
	fmt.Println("Pivot", luf.pivot)
	luf.Invert(M)
	fmt.Println("M^-1", M.Vals)
	fmt.Println("Det", luf.Determinant())

	fmt.Println("------------------------------")
	D := mat64.NewDense(3, 3, []float64{
		1, 3, 5,
		2, 4, 7,
		1, 1, 0,
	})

	goLU := mat64.LU(D)
	//inv, _ := mat64.Inverse(D)
	inv := goLU.Solve(mat64.NewDense(3, 3, []float64{
		1, 0, 0,
		0, 1, 0, 
		0, 0, 1,
	}))
	for i := 0; i < 3; i++ {
        row := inv.Row(nil, i)
        fmt.Println(row)
    }
}
