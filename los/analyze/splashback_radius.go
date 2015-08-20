package analyze


type splashbackRadiusParams struct { dLim float64 }
type internalSplashbackRadiusOption func(*splashbackRadiusParams)
type SplashbackRadiusOption internalSplashbackRadiusOption

// DLim sets limit for d ln(rho) / d ln(r) above which point cannot be the
// splashback radius. The default value is -5.
func DLim(dLim float64) SplashbackRadiusOption {
	return func(p *splashbackRadiusParams) { p.dLim = dLim }
}

func (p *splashbackRadiusParams) loadOptions(opts []SplashbackRadiusOption) {
	p.dLim = -5
	for _, opt := range opts { opt(p) }
}

func SplashbackRadius(
	rs, rhos, derivs []float64, opts ...SplashbackRadiusOption,
) (r float64, ok bool) {
	p := new(splashbackRadiusParams)
	p.loadOptions(opts)

	if len(rhos) != len(derivs) { panic("len(rhos) != len(derivs)") }
	if len(rhos) == 0 { return 0, false }

	rhoMin := rhos[0]
	dMin, iMin := p.dLim, -1
	for i := 1; i < len(rs) - 1; i++ {
		if rhos[i] < rhoMin {
			rhoMin = rhos[i]
			if isMinimum(derivs, i) && (derivs[i] < dMin) {
				dMin, iMin = derivs[i], i
			}
		} 
	}

	if iMin == -1 { return 0, false }
	return rs[iMin], true
}

func isMinimum(xs []float64, i int) bool {
	return xs[i] > xs[i+1] && xs[i] < xs[i-1]
}
