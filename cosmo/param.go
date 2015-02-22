package cosmo

import (
	"math"
)

// HubbleFrac calculates h(z) = H(z)/H0. Here H(z) is from Hubble's Law,
// H(z)**2 + k (c/a)**2 = H0**2 h100**2 (OmegaR a**-4 + OmegaM a**-3 + OmegaL).
// The hubble's constant in const.go is H0 = H(z = 0). An alternate
// formulation is h(a) = da/dt / (a H0). Assumes k, r = 0.
func HubbleFrac(omegaM, omegaL, z float64) float64 {
	return math.Sqrt(omegaM*math.Pow(1.0+z, 3.0) + omegaL)
}

// (And by "Mks", I mean "Mks/h".)
func rhoCriticalMks(H0, omegaM, omegaL, z float64) float64 {
	H0Mks := (H0 * 1000) / MpcMks
	H100 := H0 / 100
	// m = m * H100
	H0MksH := H0Mks / H100

	H := HubbleFrac(omegaM, omegaL, z) * H0MksH
	return 3.0 * H * H / (8.0 * math.Pi * GMks)
}

// RhoCritical calculates the critical density of the universe. This shows
// up (among other places) in halo definitions and in the definitions of
// the omages (OmegaFoo = pFoo / pCritical).  The returned value is in
// comsological units / h.
func RhoCritical(H0, omegaM, omegaL, z float64) float64 {
	return rhoCriticalMks(H0, omegaM, omegaL, z) * math.Pow(MpcMks, 3) / MSunMks
}

// RhoAverage calculates the average density of matter in the universe. The
// returned value is in cosmological units.
func RhoAverage(H0, omegaM, omegaL, z float64) float64 {
	return RhoCritical(H0, omegaM, omegaL, 0) * omegaM * math.Pow(1+z, 3.0)
}
