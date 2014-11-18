/*package geom contains routines for computing geometrical quantities in a
box with periodic boundary conditions. */
package geom

import (
	"math"
)

// Vec represents a 3D vector. All vector methods which return vectors will
// also contain *Self() and *At() variants which compute the operation in-place
// and at the specified location, respectively. All output vectors are valid
// unless they overlap with an input vector but are not equal to that vector.
type Vec [3]float32

// Scale multiplies all components of a vector by a constant.
func (v *Vec) Scale(k float64) *Vec {
	return v.ScaleAt(k, &Vec{})
}

func (v *Vec) ScaleSelf(k float64) *Vec {
	return v.ScaleAt(k, v)
}

func (v *Vec) ScaleAt(k float64, out *Vec) *Vec {
	for i := 0; i < 3; i++ {
		out[i] = float32(float64(v[i]) * k)
	}

	return out
}

// Mod calculates the value of a vector within the fundamental domain of a box
// with period boundary condiitions of the given width.
func (v *Vec) Mod(width float64) *Vec {
	return v.ModAt(width, &Vec{})
}

func (v *Vec) ModSelf(width float64) *Vec {
	return v.ModAt(width, v)
}

func (v *Vec) ModAt(width float64, out *Vec) *Vec {
	for i := 0; i < 3; i++ {
		out[i] = float32(math.Mod(float64(v[i]), width))
		if out[i] < 0 {
			out[i] += float32(width)
		}
	}

	return out
}

// Add adds two vectors together.
func (v1 *Vec) Add(v2 *Vec) *Vec {
	return v1.AddAt(v2, &Vec{})
}

func (v1 *Vec) AddSelf(v2 *Vec) *Vec {
	return v1.AddAt(v2, v1)
}

func (v1 *Vec) AddAt(v2, out *Vec) *Vec {
	for i := 0; i < 3; i++ {
		out[i] = v1[i] + v2[i]
	}

	return out
}

// Sub caluculates the dispacement vector between two vectors,
// assuming a box of the given width with periodic boundary conditions.
func (v1 *Vec) Sub(v2 *Vec, width float64) *Vec {
	return v1.SubAt(v2, width, &Vec{})
}

func (v1 *Vec) SubSelf(v2 *Vec, width float64) *Vec {
	return v1.SubAt(v2, width, v1)
}

func (v1 *Vec) SubAt(v2 *Vec, width float64, out *Vec) *Vec {
	w2 := width / 2.0
	for i := 0; i < 3; i++ {
		diff := math.Mod(float64(v1[i]-v2[i]), width)

		if diff > w2 {
			diff -= width
		} else if diff < -w2 {
			diff += width
		}

		out[i] = float32(diff)
	}

	return out
}

// Norm computes the norm of a vector.
func (v *Vec) Norm() float64 {
	return math.Sqrt(v.Dot(v))
}

// Dot computes the dot product of two vectors.
func (v1 *Vec) Dot(v2 *Vec) float64 {
	sum := 0.0
	for i := 0; i < 3; i++ {
		sum += float64(v1[i] * v2[i])
	}
	return sum
}

// Cross computes the cross product of two vectors.
func (v1 *Vec) Cross(v2 *Vec) *Vec {
	return v1.CrossAt(v2, &Vec{})
}

func (v1 *Vec) CrossSelf(v2 *Vec) *Vec {
	return v1.CrossAt(v2, v1)
}

func (v2 *Vec) CrossAt(v1, out *Vec) *Vec {
	out0 := v1[1]*v2[2] - v1[2]*v2[1]
	out1 := v1[2]*v2[0] - v1[0]*v2[2]
	out2 := v1[0]*v2[1] - v1[1]*v2[0]
	out[0], out[1], out[2] = out0, out1, out2
	return out
}
