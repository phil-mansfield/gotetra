package density

import (
	"fmt"
	"math/rand"
	"testing"

	"github.com/phil-mansfield/gotetra/geom"
)


func PrintGrid(grid []float32, cells int) {
	idx := 0
	for z := 0; z < cells; z++ {
		for y := 0; y < cells; y++ {
			for x := 0; x < cells; x++ {
				fmt.Printf("%d ", int(grid[idx]))
				idx++
			}
			fmt.Println()
		}
		fmt.Println()
	}
}

func TestAddBuffer(t *testing.T) {
	cells := 4

	grid := make([]float32, cells * cells * cells)
	cb := &geom.CellBounds{[3]int{0, 0, 0}, [3]int{3, 2, 2}}
	buf := make([]float32, cb.Width[0] * cb.Width[1] * cb.Width[2])	

	for i := range buf { buf[i] = 1.0 }
	for i := range grid { grid[i] = 1.0 }

	AddBuffer(grid, buf, cb, cells)
	//PrintGrid(grid, cells)
}

func BenchmarkAddBuffer(b *testing.B) {
	grid := make([]float32, 512 * 512 * 512)
	buf := make([]float32, 64 * 64 * 64)
	for i := range buf { buf[i] = 0 }

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		cb := &geom.CellBounds{
			[3]int{rand.Intn(512), rand.Intn(512), rand.Intn(512)},
			[3]int{64, 64, 64},
		}

		AddBuffer(grid, buf, cb, 512)
	}
}
