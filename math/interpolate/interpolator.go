package interpolate

type Interpolator interface {
	Eval(x float64) float64
	EvalAll(xs []float64, out ...[]float64) []float64
}

var (
	_ Interpolator = &Spline{}
	_ Interpolator = &Linear{}
)

type BiInterpolator interface {
	Eval(x, y float64) float64
	EvalAll(xs, ys []float64, out ...[]float64) []float64
}

var (
	_ BiInterpolator = &BiLinear{}
	_ BiInterpolator = &BiCubic{}
)

type TriInterpolator interface {
	Eval(x, y, z float64) float64
	EvalAll(xs, ys, zs []float64, out ...[]float64) []float64
}

var (
	_ TriInterpolator = &TriLinear{}
	_ TriInterpolator = &TriCubic{}
)
