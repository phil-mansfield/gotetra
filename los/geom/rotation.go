package geom

import (
	"math"

	"github.com/phil-mansfield/gotetra/math/mat"
)

// EulerMatrix creates a 3D rotation matrix based off the Euler angles phi,
// theta, and psi. These represent three consecutive rotations arounfd the z,
// x, and z axes, respectively.
//
// This is really slow right now, but I'm seriously confused by Euler angles :(
func EulerMatrix(phi, theta, psi float32) *mat.Matrix32 {
	c1, s1 := float32(math.Cos(float64(phi))), float32(math.Sin(float64(phi)))
	c2, s2:=float32(math.Cos(float64(theta))),float32(math.Sin(float64(theta)))
	c3, s3 := float32(math.Cos(float64(psi))), float32(math.Sin(float64(psi)))

	R1 := []float32{c1, -s1, 0, s1, c1, 0, 0, 0, 1}
	R2 := []float32{1, 0, 0, 0, c2, -s2, 0, s2, c2}
	R3 := []float32{c3, -s3, 0, s3, c3, 0, 0, 0, 1}
	m1 := mat.NewMatrix32(R1, 3, 3)
	m2 := mat.NewMatrix32(R2, 3, 3)
	m3 := mat.NewMatrix32(R3, 3, 3)
	return m3.Mult(m2.Mult(m1))
}

// EulerMatrixBetween creates a 3D rotation matrix which such that M * v1 = v2.
func EulerMatrixBetween(v1, v2 *Vec) *mat.Matrix32 {
	x1, y1, z1 := v1[0], v1[1], v1[2]
	phi1, theta1 := SphericalAngles(x1, y1, z1)
	x2, y2, z2 := v2[0], v2[1], v2[2]
	phi2, theta2 := SphericalAngles(x2, y2, z2)
	pi2 := float32(math.Pi / 2)
	return EulerMatrix(pi2 - phi1, theta1 - theta2, phi2 - pi2)
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

// SphericalAngles returns the azimuthal and polar angles (respectively) of
// the point specified by x, y, z.
func SphericalAngles(x, y, z float32) (phi, theta float32) {
	phi = float32(math.Atan2(float64(y), float64(x)))
	r := math.Sqrt(float64(x*x + y*y +z*z))
	theta = float32(math.Acos(float64(z) / r))
	return phi, theta
}

// UnitSphericalAngles returns the azimuthal and polar angles (respectively)
// of the point specified by x, y, z. It is assumed that
// Sqrt(x*x + y*y + z*z) = 1
func UnitSphericalAngles(x, y, z float32) (phi, theta float32) {
	phi = float32(math.Atan2(float64(y), float64(x)))
	theta = float32(math.Acos(float64(z)))
	return phi, theta
}

// PolarAngle returns the polar angle of the point as x, y.
func PolarAngle(x, y float32) float32 {
	return float32(math.Atan2(float64(y), float64(x)))
}

// AngularDistance computes the signed distance from phi1 to phi2 in the
// range [0, 2 pi).
func AngularDistance(phi1, phi2 float32) float32 {
	dPhi := phi2 - phi1
	if dPhi > math.Pi {
		return dPhi - 2*math.Pi
	} else if dPhi < -math.Pi {
		return dPhi + 2*math.Pi
	} else {
		return dPhi
	}
}

// AngleInRange takes two angles in the range [0, 2 pi) and an angular width in
// the range [0, 2 pi]. True is returned if phi is within the range inclusive
// range specified by low and width and false is returned otherwise.
func AngleInRange(phi, low, width float32) bool {
	high := low + width
	if high > 2*math.Pi {
		high -= 2*math.Pi
		return phi >= low || phi <= high
	} else {
		return phi >= low && phi <= high
	}
}
