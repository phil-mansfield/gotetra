package geom

import (
	"log"
)

type DerivOptions struct {
	Periodic bool
	Op DerivOp
	Order int
}

type DerivOp int
const (
	Add DerivOp = iota
	Subtract
	None
)

var DerivOptionsDefault = &DerivOptions{ false, None, 4 }
func (g *GridLocation) Deriv(
	vals, out []float32, axis int,
	opt *DerivOptions,
) {
	if opt == nil { opt = DerivOptionsDefault }

	if axis > 2 || axis < 0 {
		log.Fatalf("Unrecognized axis %d", axis)
	} else if g.Width[axis] < 3 {
		log.Fatalf("Width of array must be at least 3 for a " +
			"second-order derivative.")
	} else if opt.Order != 2 && opt.Order != 4 {
		log.Fatalf("Grid.Deriv() can only compute 2nd and 4th order " +
			"derivatives, not %d order.", opt.Order)
	}

	lBounds := [3]int{0, 0, 0}
	uBounds := g.Width
	lBounds[axis] += opt.Order / 2
	uBounds[axis] -= opt.Order / 2
	dx := g.BoxWidth / float64(g.Cells)

	d := newDerivTaker(opt.Order, vals, out, dx, opt.Op)

	di := 1
	if axis == 1 {
		di = g.Length
	} else if axis == 2 {
		di = g.Area
	}

	for x := lBounds[0]; x < uBounds[0]; x++ {
		for y := lBounds[1]; y < uBounds[1]; y++ {
			for z := lBounds[2]; z < uBounds[2]; z++ {
				idx := x + y * g.Length + z * g.Area
				d.flatDeriv(idx, di)
			}
		}
	}
	
	if !opt.Periodic || g.Cells !=  g.Width[axis]{
		g.edgesDeriv(axis, opt, lBounds, uBounds, d)
	} else {
		g.wrapDeriv(axis, opt, lBounds, uBounds, d)
	}
}

func (g *Grid) wrapDeriv(
	axis int, opt *DerivOptions, lBounds, uBounds [3]int, d derivTaker,
) {
	var aIdx int
	for pos := -opt.Order/2; pos <= opt.Order/2; pos++ {
		if pos == 0 { continue }
		if pos > 0 { aIdx = pos - 1}
		if pos < 0 { aIdx = g.Width[axis] + pos }
		switch axis {
		case 0:
			x := aIdx
			for z := lBounds[2]; z < uBounds[2]; z++ {
				for  y := lBounds[1]; y < uBounds[1]; y++ {
					idx := x + y * g.Length + z * g.Area
					d.wrapDeriv(idx, 1, pos, g.Width[0])
				}
			}
		case 1:
			y := aIdx
			for z := lBounds[2]; z < uBounds[2]; z++ {
				for  x := lBounds[0]; x < uBounds[0]; x++ {
					idx := x + y * g.Length + z * g.Area
					d.wrapDeriv(idx, g.Length, pos, g.Width[1])
				}
			}
		case 2:
			z := aIdx
			for y := lBounds[1]; y < uBounds[1]; y++ {
				for  x := lBounds[0]; x < uBounds[0]; x++ {
					idx := x + y * g.Length + z * g.Area
					d.wrapDeriv(idx, 1, pos, g.Width[2])
				}
			}
		}
	}
}

func (g *Grid) edgesDeriv(
	axis int, opt *DerivOptions, lBounds, uBounds [3]int, d derivTaker,
) {
	var aIdx int
	for pos := -opt.Order/2; pos <= opt.Order/2; pos++ {
		if pos == 0 { continue }
		if pos > 0 { aIdx = pos - 1}
		if pos < 0 { aIdx = g.Width[axis] + pos }
		switch axis {
		case 0:
			x := aIdx
			for z := lBounds[2]; z < uBounds[2]; z++ {
				for  y := lBounds[1]; y < uBounds[1]; y++ {
					idx := x + y * g.Length + z * g.Area
					d.edgeDeriv(idx, 1, pos)
				}
			}
		case 1:
			y := aIdx
			for z := lBounds[2]; z < uBounds[2]; z++ {
				for  x := lBounds[0]; x < uBounds[0]; x++ {
					idx := x + y * g.Length + z * g.Area
					d.edgeDeriv(idx, g.Length, pos)
				}
			}
		case 2:
			z := aIdx
			for y := lBounds[1]; y < uBounds[1]; y++ {
				for  x := lBounds[0]; x < uBounds[0]; x++ {
					idx := x + y * g.Length + z * g.Area
					d.edgeDeriv(idx, 1, pos)
				}
			}
		}
	}
}

type derivTaker interface {
	init(vals, out []float32, denom float32)
	flatDeriv(idx, di int)
	edgeDeriv(idx, di, pos int)
	wrapDeriv(idx, di, pos, width int)

	flat(idx, di int) float32
	edge(idx, di, pos int) float32
	wrap(idx, di, pos, width int) float32
}

type derivBase struct {
	vals, out []float32
	denom float32
}

func (d *derivBase) init(vals, out []float32, denom float32) {
	d.vals, d.out, d.denom = vals, out, denom
}

type deriv2Base struct { derivBase }
type deriv4Base struct { derivBase }

func (d *deriv2Base) flat(idx, di int) float32 {
	return (d.vals[idx + di] - d.vals[idx - di]) / d.denom
}

func (d *deriv4Base) flat(idx, di int) float32 {
	return (-d.vals[idx + 2 *di] + 8*d.vals[idx + di] +
		-8*d.vals[idx - di] + d.vals[idx - 2*di]) / d.denom
}

func (d *deriv2Base) edge(idx, di, pos int) float32 {
	switch pos {
	case +1:
		return  (-3 * d.vals[idx] + 4 * d.vals[idx + di] +
			-d.vals[idx + 2*di]) / d.denom
	case -1:
		return (-3 * d.vals[idx] + 4 * d.vals[idx - di] +
			-d.vals[idx - 2*di]) / -d.denom
	}
	panic(":3")
}

func (d *deriv2Base) wrap(idx, di, pos, width int) float32 {
	var lo, hi int
	switch pos {
	case +1:
		lo, hi = idx + (width - 1) * di, idx + di
	case -1:
		hi, lo = idx - (width - 1) * di, idx - di
	default:
		panic(":3")
	}
	return (d.vals[hi] - d.vals[lo]) / d.denom
}

func (d *deriv4Base) wrap(idx, di, pos, width int) float32 {
	var im2, im1, ip1, ip2 int
	switch pos {
	case +1:
		im2, im1 = idx + (width - 2) * di, idx + (width - 1) * di
		ip1, ip2 = idx + di, idx + 2 * di
	case +2:
		im2, im1 = idx + (width - 2) * di, idx - di
		ip1, ip2 = idx + di, idx + 2 * di
	case -2:
		ip2, ip1 = idx - (width - 2) * di, idx + di
		im1, im2 = idx - di, idx - 2 * di
	case -1:
		ip2, ip1 = idx - (width - 2)*di, idx - (width - 1)*di
		im1, im2 = idx - di, idx - 2*di
	default:
		panic(":3")
	}
 	return (-d.vals[ip2] + 8*d.vals[ip1] +
		-8*d.vals[im1] + d.vals[im2]) / d.denom
}

func (d *deriv4Base) edge(idx, di, pos int) float32 {
	v := d.vals
	switch pos {
	case +1:
		return (-3 * v[idx + 4*di] + 16 * v[idx + 3*di] - 36 * v[idx + 2*di] +
			48 * v[idx + 1*di] - 25 * v[idx]) / d.denom
	case +2:
		return (-3 * v[idx - di] - 10 * v[idx] + 18 * v[idx + di] + 
			-6 * v[idx + 2*di] + v[idx + 3*di]) / d.denom
	case -2:
		return (-3 * v[idx + di] - 10 * v[idx] + 18 * v[idx - di] + 
			-6 * v[idx - 2*di] + v[idx - 3*di]) / -d.denom
	case -1:
		return (-3 * v[idx - 4*di] + 16 * v[idx - 3*di] - 36 * v[idx - 2*di] +
			48 * v[idx - di] - 25 * v[idx]) / -d.denom
	}
	panic(":3")
}

type deriv2 struct { deriv2Base }
type deriv2Add struct { deriv2Base }
type deriv2Sub struct { deriv2Base }
type deriv4 struct { deriv4Base }
type deriv4Add struct { deriv4Base }
type deriv4Sub struct { deriv4Base }

/* This whole program structure is an abomination just to save a
conditional or two. Why. */

func (d *deriv2) flatDeriv(idx, di int) { d.out[idx] = d.flat(idx, di) }
func (d *deriv2Add) flatDeriv(idx, di int) { d.out[idx] += d.flat(idx, di) }
func (d *deriv2Sub) flatDeriv(idx, di int) { d.out[idx] -= d.flat(idx, di) }

func (d *deriv4) flatDeriv(idx, di int) { d.out[idx] = d.flat(idx, di) }
func (d *deriv4Add) flatDeriv(idx, di int) { d.out[idx] += d.flat(idx, di) }
func (d *deriv4Sub) flatDeriv(idx, di int) { d.out[idx] -= d.flat(idx, di) }

func (d *deriv2) edgeDeriv(idx, di, pos int) { d.out[idx] = d.edge(idx, di, pos) }
func (d *deriv2Add) edgeDeriv(idx, di, pos int) { d.out[idx] += d.edge(idx, di, pos) }
func (d *deriv2Sub) edgeDeriv(idx, di, pos int) { d.out[idx] -= d.edge(idx, di, pos) }

func (d *deriv4) edgeDeriv(idx, di, pos int) { d.out[idx] = d.edge(idx, di, pos) }
func (d *deriv4Add) edgeDeriv(idx, di, pos int) { d.out[idx] += d.edge(idx, di, pos) }
func (d *deriv4Sub) edgeDeriv(idx, di, pos int) { d.out[idx] -= d.edge(idx, di, pos) }

func (d *deriv2) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] = d.wrap(idx, di, pos, width)
}
func (d *deriv2Add) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] += d.wrap(idx, di, pos, width)
}
func (d *deriv2Sub) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] -= d.wrap(idx, di, pos, width)
}

func (d *deriv4) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] = d.wrap(idx, di, pos, width)
}
func (d *deriv4Add) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] += d.wrap(idx, di, pos, width)
}
func (d *deriv4Sub) wrapDeriv(idx, di, pos, width int) {
	d.out[idx] -= d.wrap(idx, di, pos, width)
}

func newDerivTaker(
	order int, vals, out []float32, dx float64, op DerivOp,
) derivTaker {
	var (
		d derivTaker
		denom float32
	)

	if order == 2 {
		denom = float32(dx * 2)
		switch op {
		case None: 
			d = new(deriv2)
		case Add:
			d = new(deriv2Add)
		case Subtract:
			d = new(deriv2Sub)
		}
	} else if order == 4 {
		denom = float32(dx * 12)
		switch op {
		case None:
			d = new(deriv4)
		case Add:
			d = new(deriv4Add)
		case Subtract:
			d = new(deriv4Sub)
		}
	} else {
		log.Fatalf("Unrecognized order, %d", order)
	}

	d.init(vals, out, denom)
	return d
}

func (g *GridLocation) Curl(
	vecs, out [3][]float32,
	op *DerivOptions,
) {
	// We can't change and revert due to thread safety.
	tmp := op
	if tmp == nil {
		tmp, op = DerivOptionsDefault, DerivOptionsDefault
	}
	*op = *tmp

	g.Deriv(vecs[2], out[0], 1, op)
	g.Deriv(vecs[0], out[1], 2, op)
	g.Deriv(vecs[1], out[2], 0, op)

	op.Op = Subtract
	g.Deriv(vecs[1], out[0], 2, op)
	g.Deriv(vecs[2], out[1], 0, op)
	g.Deriv(vecs[0], out[2], 1, op)
}

func (g *GridLocation) Divergence(
	vecs [3][]float32, out []float32,
	op *DerivOptions,
) {
	tmp := op
	if tmp == nil {
		tmp, op = DerivOptionsDefault, DerivOptionsDefault
	}
	*op = *tmp
	g.Deriv(vecs[0], out, 0, op)
	op.Op = Add
	g.Deriv(vecs[1], out, 1, op)
	g.Deriv(vecs[2], out, 2, op)
}

func (g *GridLocation) Gradient(
	vals []float32, out [3][]float32,
	op *DerivOptions,
) {
	if op == nil { op = DerivOptionsDefault }
	g.Deriv(vals, out[0], 0, op)
	g.Deriv(vals, out[1], 1, op)
	g.Deriv(vals, out[2], 2, op)
}
