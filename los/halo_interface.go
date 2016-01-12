package los

import (
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/geom"
)

// Halo is a _very leaky_ abstraction around the different types of halos.
// Mainly provided as a convenience for the already terrible gtet_shell.go
// file.
type Halo interface {
	GetRs(buf []float64)
	GetRhos(ring, losIdx int, buf []float64)
	MeanProfile() []float64
	MedianProfile() []float64
	Phi(i int) float64
	LineSegment(ring, losIdx int, out *geom.LineSegment)
	SheetIntersect(hd *io.SheetHeader) bool
	PlaneToVolume(ring int, px, py float64) (x, y, z float64)
	RMax() float64
}

// typechecking
var _ Halo = &HaloProfiles{}
