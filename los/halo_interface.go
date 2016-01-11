package los

// Halo is a _very leaky_ abstraction around the different types of halos.
// Mainly provided as a convenience for the already terrible gtet_shell.go
// file.
type Halo interface {
	GetRs(buf []float64)
	GetRhos(ring, losIdx int, buf []float64)
	MeanProfile() []float64
	MedianProfile() []float64
}

// typechecking
var _ Halo = &HaloProfiles{}
