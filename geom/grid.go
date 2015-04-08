package geom

// Grid provides an interface for reasoning over a 1D slice as if it were a
// 3D grid.
type Grid struct {
	CellBounds
	Length, Area, Volume int
	uBounds [3]int
}

// CellBounds represents a bounding box aligned to grid cells.
type CellBounds struct {
	Origin, Width [3]int
}

// NewGrid returns a new Grid instance.
func NewGrid(origin [3]int, width [3]int) *Grid {
	g := &Grid{}
	g.Init(origin, width)
	return g
}

// Init initializes a Grid instance.
func (g *Grid) Init(origin [3]int, width [3]int) {
	g.Origin = origin
	g.Width = width

	g.Length = width[0]
	g.Area = width[0] * width[1]
	g.Volume = width[0] * width[1] * width[2]

	for i := 0; i < 3; i++ {
		g.uBounds[i] = g.Origin[i] + g.Width[i]
	}
}

// Idx returns the grid index corresponding to a set of coordinates.
func (g *Grid) Idx(x, y, z int) int {
	// Those subtractions are actually unneccessary.
	return ((x - g.Origin[0]) + (y-g.Origin[1])*g.Length +
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

// BoundsCheck returns true if the given coordinates are within the Grid and
// false otherwise.
func (g *Grid) BoundsCheck(x, y, z int) bool {
	return (g.Origin[0] <= x && g.Origin[1] <= y && g.Origin[2] <= z) &&
		(x < g.uBounds[0] && y < g.uBounds[1] &&
			z < g.uBounds[2])
}

// Coords returns the x, y, z coordinates of a point from its grid index.
func (g *Grid) Coords(idx int) (x, y, z int) {
	x = idx % g.Length
	y = (idx % g.Area) / g.Length
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

// Intersect retursn true if the two bounding boxes overlap and false otherwise.
func (cb1 *CellBounds) Intersect(cb2 *CellBounds, width int) bool {
	intr := true
	var ( 
		oSmall, oBig, wSmall, wBig int
	)
	for i := 0; intr && i < 3; i++ {
		if cb1.Width[i] < cb2.Width[i] {
			oSmall, wSmall = cb1.Origin[i], cb1.Width[i]
			oBig, wBig = cb2.Origin[i], cb2.Width[i]
		} else {
			oSmall, wSmall = cb2.Origin[i], cb2.Width[i]
			oBig, wBig = cb1.Origin[i], cb1.Width[i]
		}

		eSmall := oSmall + wSmall
		beSmall := bound(eSmall, oBig, width)
		boSmall := bound(oSmall, oBig, width)

		intr = intr && (beSmall < wBig || boSmall < wBig)
	}
	return intr
}

func bound(x, origin, width int) int {
	diff := x - origin
	if diff < 0 { return diff + width }
	if diff > width { return diff - width }
	return diff
}

func (cb *CellBounds) ScaleVecs(vs []Vec, cells int, boxWidth float64) {
	fCells := float64(cells)
	scale := fCells / boxWidth

	origin := &Vec{
		float32(cb.Origin[0]), float32(cb.Origin[1]), float32(cb.Origin[2]),
	}

	for i := range vs {
		vs[i].ScaleSelf(scale)
		vs[i].SubSelf(origin, fCells)
	}
}

func fMinMax(min, max, x float32) (float32, float32) {
	if x < min {
		return x, max
	}
	if x > max {
		return min, x
	}
	return min, max
}
