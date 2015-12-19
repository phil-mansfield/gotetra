package loop

import (
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/los/geom"
)

// InterpolatorType is the type of interpolation used during the iteration.
type InterpolatorType int

const (
	Cubic InterpolatorType = iota
	Linear
)

// An object which will be passed to Loop(). It must be less than half the
// size of the simulation box (because delaing with larger objects is a
// nightmare).
type Object interface {
	// Transform does a coordinate transformation on a slice of vectors such
	// that they are all .
	Transform(vecs []rgeom.Vec, tw float64)
	// Contains returns true if a point is contained inside the Object.
	Contains(x, y, z float64) bool
	// IntersectBox returns true if a box intersects the Object. This function
	// is allowed to be slow.
	IntersectBox(origin, span [3]float64, tw float64) bool

	// ThreadCopy copies enough of the object's internal state to a new instance
	// of that object so that both Objects can maintain calls to Use*foo*
	// simultaneously.
	ThreadCopy(id, threads int) Object
	// ThreadMerge combines the data from a group of objects into the reciever.
	// The reciever will _not_ be in the object list, but its initial data
	// should be combined likewise.
	ThreadMerge(objs []Object)

	// UsesTetra returns true if the Object requires a loop over tetrahedra
	// and false if it requires a loop over points.
	UsesTetra() bool
	// UsePoint supplies a point to the object.
	UsePoint(x, y, z float64)
	// UseTetra supplies a tetrahedron to the object.
	UseTetra(t *geom.Tetra)
}

// Loop iterates over all the tetrahedra or all the points in the simulation
// box. A lot of tricks are used to minimize the number of operations done.
// pts is the number of 
func Loop(objs []Object, intr InterpolatorType, pts int) {
	panic("NYI")
}
