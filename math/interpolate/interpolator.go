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

	EvalAllX(x float64, ys []float64, out ...[]float64) []float64
	EvalAllY(xs []float64, y float64, out ...[]float64) []float64
}

var (
	_ BiInterpolator = &BiLinear{}
)

type TriInterpolator interface {
	Eval(x, y, z float64) float64
	EvalAll(xs, ys, zs []float64, out ...[]float64) []float64

	EvalAllXY(x, y float64, zs []float64, out ...[]float64) []float64
	EvalAllYZ(x float64, ys []float64, z float64,  out ...[]float64) []float64
	EvalAllZZ(xs []float64, y, z float64, out ...[]float64) []float64
}

var (
	_ TriInterpolator = &TriLinear{}
)
