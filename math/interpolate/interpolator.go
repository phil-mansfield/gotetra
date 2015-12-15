package interpolate

type Interpolator interface {
	Eval(x float64) float64
	EvalAll(xs []float64, out ...[]float64) []float64
	Ref() Interpolator
}

var (
	_ Interpolator = &Spline{}
	_ Interpolator = &Linear{}
	_ Interpolator = &splineRef{}
	_ Interpolator = &linearRef{}
)

type BiInterpolator interface {
	Eval(x, y float64) float64
	EvalAll(xs, ys []float64, out ...[]float64) []float64
	Ref() BiInterpolator
}

var (
	_ BiInterpolator = &BiLinear{}
	_ BiInterpolator = &BiCubic{}
	_ BiInterpolator = &biLinearRef{}
	_ BiInterpolator = &biCubicRef{}
)

type TriInterpolator interface {
	Eval(x, y, z float64) float64
	EvalAll(xs, ys, zs []float64, out ...[]float64) []float64
	Ref() TriInterpolator
}

var (
	_ TriInterpolator = &TriLinear{}
	_ TriInterpolator = &TriCubic{}
	_ TriInterpolator = &triLinearRef{}
	_ TriInterpolator = &triCubicRef{}
)
