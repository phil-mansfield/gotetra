package interpolate

import (
	"encoding/binary"
	"testing"

	"github.com/stretchr/testify/assert"
)

func value(x, y, z float64) float64 {
	return 2*x + 3*y + 5*z
}

func TestUniformTriLinear(t *testing.T) {
	minVal := 0.0
	n := 11
	step := 0.1
	vals := make([]float64, n*n*n)
	idx := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				vals[idx] = value(minVal+float64(i)*step, minVal+float64(j)*step, minVal+float64(k)*step)
				idx++
			}
		}
	}
	interp := NewUniformTriLinear(
		minVal, step, n,
		minVal, step, n,
		minVal, step, n,
		vals, binary.BigEndian,
	)
	// points on the grid should work
	assert.Equal(t, value(0.5, 0.5, 0.5), interp.Eval(0.5, 0.5, 0.5), "on grid")
	// points just off the grid should also work
	assert.Equal(t, value(0.51, 0.50, 0.50), interp.Eval(0.51, 0.50, 0.50), "nearby x")
	assert.Equal(t, value(0.50, 0.51, 0.50), interp.Eval(0.50, 0.51, 0.50), "nearby y")
	assert.Equal(t, value(0.50, 0.50, 0.51), interp.Eval(0.50, 0.50, 0.51), "nearby z")
	// points on the edge of the grid should work
	assert.Equal(t, value(0, 0, 0), interp.Eval(0, 0, 0), "grid edge")
	assert.Equal(t, value(0.01, 0, 0), interp.Eval(0.01, 0, 0), "grid edge nearby x")
}

func TestUniformTriLinearLittleEndian(t *testing.T) {
	minVal := 0.0
	n := 11
	step := 0.1
	vals := make([]float64, n*n*n)
	idx := 0
	for k := 0; k < n; k++ {
		for j := 0; j < n; j++ {
			for i := 0; i < n; i++ {
				vals[idx] = value(minVal+float64(i)*step, minVal+float64(j)*step, minVal+float64(k)*step)
				idx++
			}
		}
	}
	interp := NewUniformTriLinear(
		minVal, step, n,
		minVal, step, n,
		minVal, step, n,
		vals, binary.LittleEndian,
	)
	// points just off the grid should also work
	assert.Equal(t, value(0.51, 0.50, 0.50), interp.Eval(0.51, 0.50, 0.50), "nearby x")
	assert.Equal(t, value(0.50, 0.51, 0.50), interp.Eval(0.50, 0.51, 0.50), "nearby y")
	assert.Equal(t, value(0.50, 0.50, 0.51), interp.Eval(0.50, 0.50, 0.51), "nearby z")
}

func TestUniformTriLinearMulti(t *testing.T) {
	minVal := 0.0
	n := 11
	step := 0.1
	vals := make([][]float64, n*n*n)
	idx := 0
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			for k := 0; k < n; k++ {
				vals[idx] = make([]float64, 1)
				vals[idx][0] = value(minVal+float64(i)*step, minVal+float64(j)*step, minVal+float64(k)*step)
				idx++
			}
		}
	}
	interp := NewUniformTriLinearMulti(
		minVal, step, n,
		minVal, step, n,
		minVal, step, n,
		vals, binary.BigEndian,
	)
	// assert.Equal(t, 1, interp.valsDim)
	// points on the grid should work
	assert.Equal(t, value(0.5, 0.5, 0.5), interp.Eval(0.5, 0.5, 0.5)[0], "on grid")
	// points just off the grid should also work
	assert.Equal(t, value(0.51, 0.50, 0.50), interp.Eval(0.51, 0.50, 0.50)[0], "nearby x")
	assert.Equal(t, value(0.50, 0.51, 0.50), interp.Eval(0.50, 0.51, 0.50)[0], "nearby y")
	assert.Equal(t, value(0.50, 0.50, 0.51), interp.Eval(0.50, 0.50, 0.51)[0], "nearby z")

	// points on the edge of the grid should work
	assert.Equal(t, value(0, 0, 0), interp.Eval(0, 0, 0)[0], "grid edge")
	assert.Equal(t, value(0.01, 0, 0), interp.Eval(0.01, 0, 0)[0], "grid edge nearby x")
}
