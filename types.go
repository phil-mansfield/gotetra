package gotetra

type Particle struct {
	Xs [3]float64
	Vs [3]float64
	Id int64
}

// Header describes meta-information about the current catalog.
type Header struct {
	Cosmo CosmologyHeader

	Mass float64  // Mass of one particle
	Count int64 // Number of particles in catalog
	TotalCount int64 // Number of particles in all catalogs
	CountWidth int64 // Number of particles "on one side": TotalCount^(1/3)

	Idx int64 // Index of catalog: x-major ordering is used
	GridWidth int64 // Number of gird cells "on one side"
	Width int64 // Width of the catalog's bounding box
	TotalWidth float64 // Width of the sim's bounding box
}

type CosmologyHeader struct {
	Z float64
	OmegaM float64
	OmegaL float64
	H100 float64
}
