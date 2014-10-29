package gotetra

type Particle struct {
	Xs [3]float64
	Vs [3]float64
	Id int64
}

type Header struct {
	Cosmo CosmologyHeader
	Count, TotalCount int64
	Mass, BoxSize float64
}

type CosmologyHeader struct {
	Z float64
	OmegaM float64
	OmegaL float64
	H100 float64
}
