package geom

// Grid provides an interface for reasoning over a 1D slice as if it were a
// 3D grid.
type Grid struct {
	Origin  [3]int
	Width   int
	Area    int
	Volume  int
	uBounds [3]int
}

// CellBounds represents a bounding box aligned to grid cells.
type CellBounds struct {
	Origin, Width [3]int
}

// NewGrid returns a new Grid instance.
func NewGrid(origin *[3]int, width int) *Grid {
	g := &Grid{}
	g.Init(origin, width)
	return g
}

// Init initializes a Grid instance.
func (g *Grid) Init(origin *[3]int, width int) {
	g.Origin = *origin
	g.Width = width
	g.Area = width * width
	g.Volume = width * width * width

	for i := 0; i < 3; i++ {
		g.uBounds[i] = g.Origin[i] + g.Width
	}
}

// Bounds returns a Grid's bounding box.
func (g *Grid) CellBounds() *CellBounds {
	b := &CellBounds{}
	g.CellBoundsAt(b)
	return b
}

// BoundsAt puts a Grid' bounding box at the specified location.
func (g *Grid) CellBoundsAt(out *CellBounds) {
	out.Origin = g.Origin
	for d := 0; d < 3; d++ {
		out.Width[d] = g.Width
	}
}

// Idx returns the grid index corresponding to a set of coordinates.
func (g *Grid) Idx(x, y, z int) int {
	return ((x - g.Origin[0]) + (y-g.Origin[1])*g.Width +
		(z-g.Origin[2])*g.Area)
}

// IdxCheck returns an index and true if the given coordinate are valid and
// false otherwise.
func (g *Grid) IdxCheck(x, y, z int) (idx int, ok bool) {
	if !g.BoundsCheck(x, y, z) {
		return -1, false
	}

	return g.Idx(x, y, z), true
}

// TODO: make this correct. Doesn't account for the origin.
func (g *Grid) Wrap(x, y, z int) (wx, wy, wz int) {
	wx = x % g.Width
	if wx < 0 {
		wx += g.Width
	}
	wy = y % g.Width
	if wy < 0 {
		wy += g.Width
	}
	wz = z % g.Width
	if wz < 0 {
		wz += g.Width
	}
	return wx, wy, wz
}

// BoundsCheck returns true if the given coordinates are within the Grid and
// false otherwise.
func (g *Grid) BoundsCheck(x, y, z int) bool {
	return (g.Origin[0] <= x && g.Origin[1] <= y && g.Origin[2] <= z) &&
		(x < g.uBounds[0] && y < g.uBounds[1] &&
			z < g.uBounds[2])
}

// Coords returns the x, y, z coordinates of a point from its grid index.
func (g *Grid) Coords(idx int) (x, y, z int) {
	x = idx % g.Width
	y = (idx % g.Area) / g.Width
	z = idx / g.Area
	return x, y, z
}

// pMod computes the positive modulo x % y.
func pMod(x, y int) int {
	m := x % y
	if m < 0 {
		m += y
	}
	return m
}

func (cb1 *CellBounds) Intersect(cb2 *CellBounds, width int) bool {
	intersect := true
	w2 := width / 2

	for i := 0; intersect && i < 3; i++ {
		diff := cb1.Origin[i] - cb2.Origin[i]
		if diff > w2 { 
			diff -= width 
		} else if diff < -w2 {
			diff += width
		}

		intersect = intersect &&
			((diff >= 0 && diff < cb2.Width[i]) || 
			(diff < 0 && -diff < cb1.Width[i]))
	}
	return intersect
}
