package halo

const (
	tail = -1
)

type Grid struct {
	Cells int
	cw, Width float64

	// Grid-sized
	Heads []int
	// Data-sized
	Next []int
}

func NewGrid(cells int, width float64, dataLen int) *Grid {
	g := &Grid{
		Cells: cells,
		cw: width / float64(cells),
		Width: width,
		Heads: make([]int, cells * cells * cells),
		Next: make([]int, dataLen),
	}

	for i := range g.Heads { g.Heads[i] = tail }

	return g
}

func (g *Grid) Length(idx int) int {
	next := g.Heads[idx]
	n := 0
	for next != tail {
		n++
		next = g.Next[next]
	}
	return n
}

func (g *Grid) Insert(xs, ys, zs []float64) {
	for i := range xs {
		x, y, z := xs[i], ys[i], zs[i]
		if x > g.Width {
			x -= g.Width
		} else if x < 0 {
			x += g.Width
		}
		if y > g.Width {
			y -= g.Width
		} else if y < 0 {
			y += g.Width
		}
		if z > g.Width {
			z -= g.Width
		} else if z < 0 {
			z += g.Width
		}
				

		ix, iy, iz := int(x/g.cw), int(y/g.cw), int(z/g.cw)
		idx := ix + iy * g.Cells + iz * g.Cells * g.Cells

		g.Next[i] = g.Heads[idx]
		g.Heads[idx] = i
	}
}

func (g *Grid) TotalCells() int {
	return len(g.Heads)
}

func (g *Grid) MaxLength() int {
	max := 0
	for i := 0; i < g.TotalCells(); i++ {
		l := g.Length(i)
		if l > max { max = l }
	}
	return max
}

func (g *Grid) AverageLength() int {
	return len(g.Next) / g.TotalCells()
}

func (g *Grid) ReadIndexes(idx int, buf []int) []int {
	buf = buf[:cap(buf)]
	
	next := g.Heads[idx]
	n := 0
	for next != tail {
		buf[n] = next
		n++
		next = g.Next[next]
	}

	return buf[:n]
}
