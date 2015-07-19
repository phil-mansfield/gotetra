package halo

import (
	"math"
	"strings"

	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/cosmo"
)

type Radius int
const (
	RVir Radius = iota
	R200c
	R200m
	R500c
)

func RadiusFromString(s string) (r Radius, ok bool) {
	s = strings.ToLower(s)
	switch s {
	case "200m", "r200m":
		return R200m, true
	case "vir", "rvir":
		return RVir, true
	case "200c", "r200c":
		return R200c, true
	case "500c", "r500c":
		return R500c, true
	}
	return RVir, false
}

func (r Radius) String() string {
	switch r {
	case R200m:
		return "R200m"
	case R200c:
		return "R200c"
	case R500c:
		return "R500c"
	case RVir:
		return "RVir"
	}
	panic(":3")
}

func (r Radius) Radius(c *io.CosmologyHeader, ms, out []float64) {
	var rho float64
	h0 := c.H100 * 100

	switch r {
	case RVir:
		rho = 177.653 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	case R200c:
		rho = 200 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	case R200m:
		rho = 200 * cosmo.RhoAverage(h0, c.OmegaM, c.OmegaL, c.Z)
	case R500c:
		rho = 500 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	}

	a := 1 / (1 + c.Z)
	factor := rho * 4 * math.Pi / 3

	for i, m := range ms {
		out[i] = math.Pow(m / factor, 1.0/3) / a
	}
}

func (r Radius) Mass(c *io.CosmologyHeader, rs, out []float64) {
	var rho float64
	h0 := c.H100 * 100

	switch r {
	case RVir:
		rho = 177.653 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	case R200c:
		rho = 200 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	case R200m:
		rho = 200 * cosmo.RhoAverage(h0, c.OmegaM, c.OmegaL, c.Z)
	case R500c:
		rho = 500 * cosmo.RhoCritical(h0, c.OmegaM, c.OmegaL, c.Z)
	}

	a := 1 / (1 + c.Z)
	factor := rho * 4 * math.Pi / 3
	for i, r := range rs {
		r = r * a
		out[i] = factor * (r * r * r)
	}
}


func (r Radius) RockstarColumn() int {
	switch r {
	case RVir:
		return 11
	case R200c:
		return 37
	case R200m:
		return 36
	case R500c:
		return 38
	}
	panic(":3")
}

func (r Radius) RockstarMass() bool {
	switch r {
	case RVir:
		return false
	case R200c, R200m, R500c:
		return true
	}
	panic(":3")
}
