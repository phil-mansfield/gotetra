package geom

// Grid provides an interface for reasoning over a 1D slice as if it were a
// 3D grid with periodic boundary conditions.
type Grid struct {
	Origin              [3]int
	Width, Area, Volume int
}

// CellBounds represents a bounding box aligned to grid cells.
type CellBounds struct {
	Min, Max [3]int
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
}

// Bounds returns a Grid's bounding box.
func (g *Grid) Bounds() *CellBounds {
	b := &CellBounds{}
	g.BoundsAt(b)
	return b
}

// BoundsAt puts a Grid' bounding box at the specified location.
func (g *Grid) BoundsAt(out *CellBounds) {
	out.Min = g.Origin
	for d := 0; d < 3; d++ {
		out.Max[d] = out.Min[d]
	}
}

// Idx returns the grid index corresponding to a set of coordinates. Periodic
// boundary conditions are assumed.
func (g *Grid) Idx(x, y, z int) int {
	return pMod(x, g.Width) + pMod(y, g.Width)*g.Width + pMod(z, g.Width)*g.Area
}

// Coords returns the x, y, z coordinates of a point from its grid index.
func (g *Grid) Coords(idx int) (x, y, z int) {
	x = idx % g.Width
	y = (idx % g.Area) / g.Width
	z = idx / g.Area
	return x, y, z
}

// pMod computes the positive module x % y.
func pMod(x, y int) int {
	m := x % y
	if m < 0 {
		m += y
	}
	return m
}
