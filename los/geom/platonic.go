package geom

import (
	"fmt"
	"math"
)

// PlatonicSolid is a polyhedron in which every face is identical.
type PlatonicSolid int

const (
	PlatonicTetrahedron PlatonicSolid = iota
	PlatonicHexahedron
	PlatonicOctahedron
	PlatonicDodecahedron
	PlatonicIcosahedron
)

// Sides returns the number of sides contianed by a Platonic solide.
func (solid PlatonicSolid) Sides() int {
	switch solid {
	case PlatonicTetrahedron:  return 4
	case PlatonicHexahedron:   return 6
	case PlatonicOctahedron:   return 8
	case PlatonicDodecahedron: return 12
	case PlatonicIcosahedron:  return 20
	default: panic(":3")
	}
}

// NewPlatonicSolid returns the flag corresponding to the Platonic solid
// with the specified number of sides. ok is returned as true if there exists
// a Platonic solid with that many sides and false otherwise.
//
// In case you forgot, valid side numbers are: 4, 6, 8, 12, 20.
func NewPlatonicSolid(sides int) (solid PlatonicSolid, ok bool) {
	switch sides {
	case 4:  return PlatonicTetrahedron,  true
	case 6:  return PlatonicHexahedron,   true
	case 8:  return PlatonicOctahedron,   true
	case 12: return PlatonicIcosahedron,  true
	case 20: return PlatonicDodecahedron, true
	default: return 0, false
	}
}

func NewUniquePlatonicSolid(sides int) (solid PlatonicSolid, ok bool) {
	switch sides {
	case 3:  return PlatonicHexahedron,    true
	case 4:  return PlatonicTetrahedron,  true
	case 6:  return PlatonicOctahedron,   true
	case 10: return PlatonicDodecahedron, true
	default: return 0, false
	}
}

var (
	platonicTetrahedronVertices [][]Vec
	platonicHexahedronVertices [][]Vec
	platonicOctahedronVertices [][]Vec
	platonicDodecahedronVertices [][]Vec
	platonicIcosahedronVertices [][]Vec
)

func init() {
	phi := float32((1 + math.Sqrt(5)) / 2)

	platonicTetrahedronVertices = [][]Vec{
		{{ 1, 1,  1}, {-1,  1, -1}, { 1, -1, -1}},
		{{-1, 1, -1}, {-1, -1,  1}, { 1, -1, -1}},
		{{ 1, 1,  1}, { 1, -1, -1}, {-1, -1,  1}},
		{{ 1, 1,  1}, {-1, -1,  1}, {-1,  1, -1}},
	}

	platonicHexahedronVertices = [][]Vec{
		{{-1, -1, -1}, { 1, -1, -1}, { 1, -1,  1}, {-1, -1,  1}},
		{{-1, -1, -1}, {-1, -1,  1}, {-1,  1,  1}, {-1,  1, -1}},
		{{-1, -1, -1}, {-1,  1, -1}, { 1,  1, -1}, { 1, -1, -1}},

		{{-1,  1, -1}, {-1,  1,  1}, { 1,  1,  1}, { 1,  1, -1}},
		{{ 1, -1, -1}, { 1,  1, -1}, { 1,  1,  1}, { 1, -1,  1}},
		{{-1, -1,  1}, { 1, -1,  1}, { 1,  1,  1}, {-1,  1,  1}},
	}

	a, b := 1 / (2 * float32(math.Sqrt(2))), float32(0.5)
	platonicOctahedronVertices = [][]Vec{
		{{-a, 0,  a}, {-a, 0, -a}, {0,  b, 0}},
		{{-a, 0, -a}, { a, 0, -a}, {0,  b, 0}},
		{{ a, 0, -a}, { a, 0,  a}, {0,  b, 0}},
		{{ a, 0,  a}, {-a, 0,  a}, {0,  b, 0}},

		{{ a, 0,  a}, { a, 0, -a}, {0, -b, 0}},
		{{-a, 0,  a}, { a, 0,  a}, {0, -b, 0}},
		{{-a, 0, -a}, {-a, 0,  a}, {0, -b, 0}},
		{{ a, 0, -a}, {-a, 0, -a}, {0, -b, 0}},
	}

	b = 1 / phi
	c := 2 - phi
	platonicDodecahedronVertices = [][]Vec{
		{{ c,  0,  1}, {-c,  0,  1}, {-b,  b,  b}, { 0,  1,  c}, { b,  b,  b}},
		{{-c,  0,  1}, { c,  0,  1}, { b, -b,  b}, { 0, -1,  c}, {-b, -b,  b}},
		{{ 0, -1, -c}, { 0, -1,  c}, {-b, -b,  b}, {-1, -c,  0}, {-b, -b, -b}},
		{{ 0, -1,  c}, { 0, -1, -c}, { b, -b, -b}, { 1, -c,  0}, { b, -b,  b}},
		{{ 1,  c,  0}, { 1, -c,  0}, { b, -b,  b}, { c,  0,  1}, { b,  b,  b}},
		{{-1, -c,  0}, {-1,  c,  0}, {-b,  b,  b}, {-c,  0,  1}, {-b, -b,  b}},

		{{ c,  0, -1}, {-c,  0, -1}, {-b, -b, -b}, { 0, -1, -c}, { b, -b, -b}},
		{{-c,  0, -1}, { c,  0, -1}, { b,  b, -b}, { 0,  1, -c}, {-b,  b, -b}},
		{{ 0,  1, -c}, { 0,  1,  c}, { b,  b,  b}, { 1,  c,  0}, { b,  b, -b}},
		{{ 0,  1,  c}, { 0,  1, -c}, {-b,  b, -b}, {-1,  c,  0}, {-b,  b,  b}},
		{{-1,  c,  0}, {-1, -c,  0}, {-b, -b, -b}, {-c,  0, -1}, {-b,  b, -b}},
		{{ 1, -c,  0}, { 1,  c,  0}, { b,  b, -b}, { c,  0, -1}, { b, -b, -b}},
	}

	a, b = float32(0.5), 1 / (2 * phi)
	platonicIcosahedronVertices = [][]Vec{
		{{ 0,  b,  a}, {-b,  a,  0}, { b,  a,  0}},
		{{ 0,  b, -a}, { b,  a,  0}, {-b,  a,  0}},
		{{ 0,  b,  a}, { 0, -b,  a}, {-a,  0,  b}},
		{{ 0,  b,  a}, { a,  0,  b}, { 0, -b,  a}},
		{{ 0,  b,  a}, {-a,  0,  b}, {-b,  a,  0}},
		{{ 0, -b, -a}, {-a,  0, -b}, {-b, -a,  0}},
		{{ 0,  b,  a}, { b,  a,  0}, { a,  0,  b}},
		{{ b,  a,  0}, { a,  0, -b}, { a,  0,  b}},
		{{ 0, -b,  a}, { a,  0,  b}, { b, -a,  0}},
		{{-b,  a,  0}, {-a,  0,  b}, {-a,  0, -b}},


		{{ 0, -b, -a}, {-b, -a,  0}, { b, -a,  0}},
		{{ 0, -b,  a}, { b, -a,  0}, {-b, -a,  0}},
		{{ 0,  b, -a}, { 0, -b, -a}, { a,  0, -b}},
		{{ 0,  b, -a}, {-a,  0, -b}, { 0, -b, -a}},
		{{ 0, -b, -a}, { b, -a,  0}, { a,  0, -b}},
		{{ 0,  b, -a}, { a,  0, -b}, { b,  a,  0}},
		{{ 0,  b, -a}, {-b,  a,  0}, {-a,  0, -b}},
		{{-b, -a,  0}, {-a,  0, -b}, {-a,  0,  b}},
		{{ 0, -b,  a}, {-b, -a,  0}, {-a,  0,  b}},
		{{ b, -a,  0}, { a,  0,  b}, { a,  0, -b}},
	}
}

// FaceVertices returns the coordinates of vertices of the specifed face.
//
// These coordinates are not normalized in any particular way.
func (solid PlatonicSolid) FaceVertices(i int) []Vec {
	if i >= solid.Sides() {
		panic(fmt.Sprintf(
			"Invalid side number %d for PlatonidSolid %d", i, solid,
		))
	}
	switch solid {
	case PlatonicTetrahedron:
		return platonicTetrahedronVertices[i]
	case PlatonicHexahedron:
		return platonicHexahedronVertices[i]
	case PlatonicOctahedron:
		return platonicOctahedronVertices[i]
	case PlatonicDodecahedron:
		return platonicDodecahedronVertices[i]
	case PlatonicIcosahedron: 
		return platonicIcosahedronVertices[i]
	default:
		panic(":3")
	}
	panic("NYI")
}

// Normals returns the normal vectors of all the faces of a Platonic solid.
//
// Note: (since I know what you're going to use this for) if you want to
// generate origin-anchored planes from these vectors, make sure to remove
// vectors which would result in duplicate planes (i.e. filter one of each pair
// of vectors that point in opposite directions.)
func (solid PlatonicSolid) Normals() []Vec {
	vs := make([]Vec, solid.Sides())
	for i := range vs {
		v := &vs[i]
		verts := solid.FaceVertices(i)
		for j := range verts {
			for dim := 0; dim < 3; dim++ {
				v[dim] += verts[j][dim]
			}
		}

		sum := float32(0)
		for dim := 0; dim < 3; dim++ { sum += v[dim]*v[dim] }
		sum = float32(math.Sqrt(float64(sum)))
		for dim := 0; dim < 3; dim++ { v[dim] /= sum }
	}

	return vs
}

// UniqueNormals returns all the face-centered normal vectors which specify
// unique origin-centered planes.
func (solid PlatonicSolid) UniqueNormals() []Vec {
	vs := solid.Normals()
	switch solid {
	case PlatonicTetrahedron:
		return vs
	default:
		return vs[:len(vs)/2]
	}
}
