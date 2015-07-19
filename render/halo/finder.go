package halo

import (
	"sort"
)

type SubhaloFinder struct {
	g *Grid
	subhaloStarts, subhalos []int
}

type Bounds struct {
	Origin, Span [3]int
}

func (b *Bounds) SphereBounds(pos [3]float64, r, cw, width float64) {
	for i := 0; i < 3; i++ {
		min, max := pos[i] - r, pos[i] + r
		if min < 0 { min += width }
		if min > max { max +=  width }
		minCell, maxCell := int(min / cw), int(max / cw)
		b.Origin[i] = minCell
		b.Span[i] = maxCell - minCell + 1
	}
}

func (b *Bounds) ConvertIndices(x, y, z, width int) (bx, by, bz int) {
	bx = x - b.Origin[0]
	if bx < 0 { bx += width }
	by = y - b.Origin[1]
	if by < 0 { by += width }
	bz = z - b.Origin[2]
	if bz < 0 { bz += width }
	return bx, by, bz
}

func (b *Bounds) Inside(val int, width int, dim int) bool {
	lo, hi := b.Origin[dim], b.Origin[dim] + b.Span[dim]
	if val >= hi {
		val -= width
	} else if val < lo {
		val += width 
	}
	return val < hi && val >= lo
}

func NewSubhaloFinder(g *Grid) *SubhaloFinder {
	i := &SubhaloFinder{
		g: g,
		subhaloStarts: make([]int, len(g.Next)),
		subhalos: make([]int, 0, int(2.5 * float64(len(g.Next)))),
	}

	return i
}

func (sf *SubhaloFinder) FindSubhalos(
	xs, ys, zs, rs []float64, mult float64,
) {
	b := &Bounds{}
	pos := [3]float64{}
	c := sf.g.Cells
	
	buf := make([]int, 0, sf.g.MaxLength())
	for i := range rs { rs[i] *= mult }
	cMax := cumulativeMax(rs)

	vol := 0.0
	for ih, rh := range rs {
		sf.startSubhalos(ih)
		maxSR := cMax[ih]

		pos[0], pos[1], pos[2] = xs[ih], ys[ih], zs[ih]
		b.SphereBounds(pos, rh + maxSR, sf.g.cw, sf.g.Width)
		vol += float64(b.Span[0] * b.Span[1] * b.Span[2])

		for dz := 0; dz < b.Span[2]; dz++ {
			z := b.Origin[2] + dz
			if z >= c { z -= c}
			zOff := z * c * c
			for dy := 0; dy < b.Span[1]; dy++ {
				y := b.Origin[1] + dy
				if y >= c { y -= c }
				yOff := y * c
				for dx := 0; dx < b.Span[0]; dx++ {
					x := b.Origin[0] + dx
					if x >= c { x -= c }
					
					idx := zOff + yOff + x

					buf = sf.g.ReadIndexes(idx, buf)
					sf.markSubhalos(ih, buf, xs, ys, zs, rs)
				}
			}
		}
	}

	sf.crossMatch()
	for i := range rs { rs[i] /= mult }
}

func (sf *SubhaloFinder) StartEnd(ih int) (start, end int) {
	if ih == len(sf.subhaloStarts) - 1 {
		return sf.subhaloStarts[ih], len(sf.subhalos)
	} else {
		return sf.subhaloStarts[ih], sf.subhaloStarts[ih+1]
	}
}

func (sf *SubhaloFinder) IntersectCount(ih int) int {
	start, end := sf.StartEnd(ih)
	return end - start
}

func (sf *SubhaloFinder) HostCount(ih int) int {
	is := sf.Intersects(ih)
	for i, sh := range is {
		if sh > ih { return i }
	}
	return len(is)
}

func (sf *SubhaloFinder) SubhaloCount(ih int) int {
	return sf.IntersectCount(ih) - sf.HostCount(ih)
}

func (sf *SubhaloFinder) Intersects(ih int) []int {
	start, end := sf.StartEnd(ih)
	return sf.subhalos[start: end]
}

func (sf *SubhaloFinder) Hosts(ih int) []int {
	return sf.Intersects(ih)[:sf.HostCount(ih)]
}

func (sf *SubhaloFinder) Subhalos(ih int) []int {
	return sf.Intersects(ih)[sf.HostCount(ih):]
}

// This is dumb for performance reasons, but this would be fast anyway, so
// I don't know why I wrote it like this.
func (sf *SubhaloFinder) crossMatch() {
	counts := make([]int, len(sf.subhaloStarts))
	for ih := range sf.subhaloStarts {
		start, end := sf.StartEnd(ih)
		for ish := start; ish < end; ish++ {
			sh := sf.subhalos[ish]
			counts[sh]++
		}
	}

	hostStarts := make([]int, len(counts))
	for i := 1; i < len(counts); i++ {
		hostStarts[i] = hostStarts[i - 1] + counts[i - 1]
	}

	hosts := make([]int, len(sf.subhalos))
	for ih := range sf.subhaloStarts {
		start, end := sf.StartEnd(ih)
		for ish := start; ish < end; ish++ {
			sh := sf.subhalos[ish]
			hosts[hostStarts[sh]] = ih
			hostStarts[sh]++
		}
	}

	for i := range counts { hostStarts[i] -= counts[i] }

	subStarts := sf.subhaloStarts
	subs := sf.subhalos
	newSubhalos := make([]int, len(hosts) + len(sf.subhalos))
	newStarts := make([]int, len(subStarts))
	j := 0
	for ih := 0; ih < len(counts); ih++ {
		newStarts[ih] = j
		sf.subhaloStarts = hostStarts
		start, end := sf.StartEnd(ih)
		for ish := start; ish < end; ish++ {
			newSubhalos[j] = hosts[ish]
			j++
		}
		sf.subhaloStarts = subStarts
		start, end = sf.StartEnd(ih)
		for ish := start; ish < end; ish++ {
			newSubhalos[j] = subs[ish]
			j++
		}
	}

	sf.subhalos = newSubhalos
	sf.subhaloStarts = newStarts

	for i := 0; i < len(sf.subhaloStarts); i++ {
		is := sort.IntSlice(sf.Intersects(i))
		sort.Sort(is)
	}
}

func (sf *SubhaloFinder) startSubhalos(i int) {
	sf.subhaloStarts[i] = len(sf.subhalos)
}

func (sf *SubhaloFinder) markSubhalos(
	ih int, idxs []int, xs, ys, zs, rs []float64,
) {
	hx, hy, hz, hr := xs[ih], ys[ih], zs[ih], rs[ih]
	for _, j := range idxs {
		if j <= ih { continue }
		sx, sy, sz, sr := xs[j], ys[j], zs[j], rs[j]
		dx, dy, dz, dr := hx - sx, hy - sy, hz - sz, hr + sr
		if dr*dr >= dx*dx + dy*dy + dz*dz {
			sf.subhalos = append(sf.subhalos, j)
		}
	}
}

func cumulativeMax(xs []float64) []float64 {
	ms := make([]float64, len(xs))
	max := xs[len(xs) - 1]
	for i := len(xs) - 1; i >= 0; i-- {
		x := xs[i]
		if x > max { max = x }
		ms[i] = max
	}
	return ms
}
