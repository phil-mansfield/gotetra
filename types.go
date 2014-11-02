package gotetra

type Particle struct {
	Xs [3]float64
	Vs [3]float64
	Id int64
}

type Header struct {
	Cosmo CosmologyHeader
	// particle info
	Mass float64 
	Count, TotalCount int64
	// binning info
	Idx, GridWidth int64
	Width, TotalWidth float64
}

type CosmologyHeader struct {
	Z float64
	OmegaM float64
	OmegaL float64
	H100 float64
}
