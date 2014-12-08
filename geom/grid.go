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
	out.Min = g.Origin
	for d := 0; d < 3; d++ {
		out.Max[d] = out.Min[d]
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

// TODO: make this correct.
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

// pMod computes the positive module x % y.
func pMod(x, y int) int {
	m := x % y
	if m < 0 {
		m += y
	}
	return m
}

func (g *Grid) Intersect(cb *CellBounds, bounds *Grid) bool {
	allContained := true

	for i := 0; i < 3; i++ {
		allContained = allContained && wrapOverlap(cb.Min[i], cb.Max[i],
			g.Origin[i], g.uBounds[i], bounds.Width)
	}

	return allContained
}

func wrapOverlap(wMin, wMax, min, max, width int) bool {
	if overlap(wMin, wMax, min, max) {
		return true
	}

	if wMax >= width {
		wMin -= width
		wMax -= width
		return overlap(wMin, wMax, min, max)
	} else if wMin < 0 {
		wMin += width
		wMax += width
		return overlap(wMin, wMax, min, max)
	}

	return false
}

func overlap(min1, max1, min2, max2 int) bool {
	return !(min1 > max2 || min2 > max1)
}
