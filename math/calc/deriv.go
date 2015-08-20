package calc

type derivParams struct { out []float64 }
type internalDerivOption func(*derivParams)
type DerivOption internalDerivOption

// Out supplies a call to Deriv with a slice to write derivatives to.
func Out(out []float64) DerivOption {
	return func(p *derivParams) { p.out = out }
}

func (p *derivParams) loadOptions(opts []DerivOption) {
	for _, opt := range opts { opt(p) }
}

// Deriv computes the numerical derivative of a a sequence of (x, y) points.
//
// The only supported orders are 2 and 4.
func Deriv(xs, ys []float64, order int, opts ...DerivOption) []float64 {
	n := len(xs)

	p := new(derivParams)
	p.loadOptions(opts)
	out := p.out
	if out == nil { out = make([]float64, n) }

	if len(ys) != n {
		panic("Length of ys and xs are not the same.")
	} else if len(out) != n {
		panic("Length of out and xs are not the same.")
	}

	if order == 0 {
		for i := range xs { out[i] = xs[i] }
	} else if order == 2 {
		for i := 1; i < n - 1; i++ {
			out[i] = (ys[i+1] - ys[i-1]) / (xs[i+1] - xs[i-1])
		}
		out[0] = (-3*ys[0] + 4*ys[1] - ys[2]) / (xs[2] - xs[0])
		out[n-1] = -(-3*ys[n-1] + 4*ys[n-2] - ys[n-3]) / (xs[n-1] - xs[n-3])
	} else if order == 4 {
		for i := 2; i < n - 2; i++ {
			out[i] = (-ys[i+2] + 8*ys[i+1] - 8*ys[i-1] + ys[i-2]) /
				(3*(xs[i+2] - xs[i-1]))
		}

		out[0] = ((-3*ys[4] + 16*ys[3] - 36*ys[2] + 48*ys[1] -
			25*ys[0]) / (3*(xs[4] - xs[0])))
		out[n-2]= ((-3*ys[n-1] - 10*ys[n-2] + 18*ys[n-3] - 6*ys[n-4] +
			ys[n-5]) / (3*(xs[n-5] - xs[n-1])))
		out[1]= ((-3*ys[0] - 10*ys[1] + 18*ys[2] - 6*ys[3] +
			ys[4]) / (3*(xs[4] - xs[0])))
		out[n-1] = ((-3*ys[n-5] + 16*ys[n-4] - 36*ys[n-3] + 48*ys[n-2] -
			25*ys[n-1]) / (3*(xs[n-5] - xs[n-1])))
	} else {
		panic("Invalid order.")
	}
	return out
}
