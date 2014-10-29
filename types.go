package gotetra

type Particle struct {
	Xs [3]float64
	Vs [3]float64
	Id uint64
}

type Header struct {
	Cosmology CosmologyHeader
	Particles int64
	ParticleMass float64
	BoxSize float64
}

type CosmologyHeader struct {
	OmegaM float64
	OmegaL float64
	H100 float64
}
