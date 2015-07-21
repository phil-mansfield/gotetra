package geom

import (
	. "math"

	"github.com/phil-mansfield/gotetra/math/mat"
)

// EulerMatrix creates a 3D rotation matrix based off the Euler angles phi,
// theta, and psi. These represent three consecutive rotations arounfd the x,
// y, and z axes, respectively.
func EulerMatrix(phi, theta, psi float64) *mat.Matrix32 {
	A := []float32{
		float32(Cos(theta)*Cos(psi)),
		float32(Cos(phi)*Sin(psi) + Sin(phi)*Sin(theta)*Cos(psi)),
		float32(Sin(phi)*Sin(psi) - Cos(phi)*Sin(theta)*Cos(psi)),
		float32(-Cos(theta)*Sin(psi)),
		float32(Cos(phi)*Cos(psi) - Sin(phi)*Sin(theta)*Sin(psi)),
		float32(Sin(phi)*Cos(psi) + Cos(phi)*Sin(theta)*Sin(psi)),
		float32(Sin(theta)),
		float32(-Sin(phi)*Cos(theta)),
		float32(Cos(phi)*Cos(theta)),
	}

	return mat.NewMatrix32(A, 3, 3)
}

// Rotate rotates a vector by the given rotation matrix.
func (v *Vec) Rotate(m *mat.Matrix32) {
	v0 := m.Vals[0]*v[0] + m.Vals[1]*v[1] + m.Vals[2]*v[2]
	v1 := m.Vals[3]*v[0] + m.Vals[4]*v[1] + m.Vals[5]*v[2]
	v2 := m.Vals[6]*v[0] + m.Vals[7]*v[1] + m.Vals[8]*v[2]
	v[0], v[1], v[2] = v0, v1, v2
}

// Rotate rotates a tetrahedron by the given rotation matrix.
func (t *Tetra) Rotate(m *mat.Matrix32){
	for i := 0; i < 4; i++ { t[i].Rotate(m) }
}
