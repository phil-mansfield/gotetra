package geom

import (
	"testing"
	"math"
	"math/rand"

	"github.com/phil-mansfield/gotetra/render/io"
	rGeom "github.com/phil-mansfield/gotetra/render/geom"
)

var (
	hx= float32(58.68211)
	hy = float32(59.48198)
	hz = float32(8.70855)
	rMax = float32(0.6318755421407911)
	rMin = float32(0)

	L = Vec{0, 0, 1}

	hd io.SheetHeader
	xs []rGeom.Vec
	ts []Tetra
	pts []PluckerTetra
	mainSuccess = myMain()
)

func randomizeTetra(t *Tetra, low, high float32) {
	for v := 0; v < 4; v++ {
		for i := 0; i < 3; i++ {
			t[v][i] = (high - low) * rand.Float32() + low
		}
	}
	t.Orient(+1)
}

func BenchmarkZPlaneSliceHit(b *testing.B) {
	n := 1000
	ts := make([]Tetra, n)
	pts := make([]PluckerTetra, n)
	poly := new(TetraSlice)
	for i := range ts {
		randomizeTetra(&ts[i], -1, +1)
		ts[i].Orient(+1)
		pts[i].Init(&ts[i])
	}

	for i := 0; i < b.N; i++ {
		idx := i % n
		ts[idx].ZPlaneSlice(&pts[idx], 0, poly)
	}
}

func BenchmarkZPlaneSliceMIss(b *testing.B) {
	n := 1000
	ts := make([]Tetra, n)
	pts := make([]PluckerTetra, n)
	poly := new(TetraSlice)
	for i := range ts {
		randomizeTetra(&ts[i], -1, +1)
		ts[i].Orient(+1)
		pts[i].Init(&ts[i])
	}

	for i := 0; i < b.N; i++ {
		idx := i % n
		ts[idx].ZPlaneSlice(&pts[idx], 1, poly)
	}
}

func BenchmarkAngleRange(b *testing.B) {
	n := 1000
	polys := make([]TetraSlice, 1)
	t, pt := new(Tetra), new(PluckerTetra)
	for len(polys) <= n {
		randomizeTetra(t, -1, +1)
		t.Orient(+1)
		pt.Init(t)

		if t.ZPlaneSlice(pt, 0, &polys[len(polys) - 1]) {
			polys = append(polys, TetraSlice{})
		}
	}

	for i := 0; i < b.N; i++ {
		idx := i % n
		polys[idx].AngleRange()
	}
}

func randomAnchoredPluckerVec() *AnchoredPluckerVec {
	x := rand.Float32()
	y := rand.Float32()
	z := rand.Float32()
	norm := float32(math.Sqrt(float64(x*x + y*y + z*z)))

	P := &Vec{0, 0, 0}
	L := &Vec{ x/norm, y/norm, z/norm }
	p := new(AnchoredPluckerVec)
	p.Init(P, L)
	return p
}


func BenchmarkPluckerTetraInit(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
	}

	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		pts[i].Init(&ts[i])
	}
}


func BenchmarkIntersectionSheet(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }

	for idx := range ts {
		pts[idx].Init(&ts[idx])
	}
	ap := new(AnchoredPluckerVec)
	P := Vec{float32(hx), float32(hy), float32(hz)}
	ap.Init(&P, &L)
	w := new(IntersectionWorkspace)

	b.ResetTimer()

	for i := 0; i < b.N; i++ {
		for idx := range ts {
			w.IntersectionDistance(&pts[idx], &ts[idx], ap)
		}
	}
}

func BenchmarkIntersectionIntersectOnly(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }

	for idx := range ts {
		pts[idx].Init(&ts[idx])
	}
	ap := new(AnchoredPluckerVec)
	P := Vec{float32(hx), float32(hy), float32(hz)}
	ap.Init(&P, &L)
	w := new(IntersectionWorkspace)

	valid := make([]bool, len(ts))
	for idx := range ts {
		re, rl, ok := w.IntersectionDistance(&pts[idx], &ts[idx], ap)
		if ok && ((re < rMax && re > rMin) || (rl < rMax && rl > rMin)) {
			valid[idx] = true
		}
	}

	b.ResetTimer()
	for i := 0; i < b.N; i++ {
		for idx := range ts {
			if valid[idx] {
				w.IntersectionDistance(&pts[idx], &ts[idx], ap)
			}
		}
	}
}

func BenchmarkSheetPreparation(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }
	vecs, phis := vectorRing(1000)
	t := new(Tetra)
	rot := EulerMatrix(0, 0, 0)

	for idx := range ts {
		*t = ts[idx]
		t.Rotate(rot)
		t.Orient(+1)
	}

	_, _ = vecs, phis

}

func BenchmarkSheetHaloIntersect(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }
	for idx := range ts {
		ts[idx].Orient(+1)
	}

	b.ResetTimer()

	var n int
	for m := 0; m < b.N; m++ {
		n = 0
		rSqr := rMax*rMax
		for i := range ts {
			for j := 0; j < 4; j++ {
				x, y, z := ts[i][j][0]-hx, ts[i][j][1]-hy, ts[i][j][2]-hz
				if rSqr >= x*x + y*y + z*z {
					n++
					break
				}
			}
		}
	}
}

func BenchmarkSheetRingIntersectionDistance(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }

	valid := make([]bool, len(ts))
	rMaxSqr, rMinSqr := rMax*rMax, rMin*rMin
	_ = rMinSqr
	for i := range ts {
		for j := 0; j < 4; j++ {
			x, y, z := ts[i][j][0]-hx, ts[i][j][1]-hy, ts[i][j][2]-hz
			rSqr := x*x + y*y + z*z
			if rSqr < rMaxSqr {
				valid[i] = true
				break
			}
		}
	}

	b.ResetTimer()

	t := new(Tetra)
	pt := new(PluckerTetra)
	poly := new(TetraSlice)
	w := new(IntersectionWorkspace)

	rot := EulerMatrix(0, 0, 0)
	dr := &Vec{-hx, -hy, -hz}
	vecs, _ := vectorRing(1000)

	var m int
	for n := 0; n < b.N; n++ {
		m = 0
		for i, ok := range valid {
			if ok {
				*t = ts[i]
				t.Translate(dr)
				t.Rotate(rot)
				pt.Init(t)
				ok := t.ZPlaneSlice(pt, 0, poly)
				if ok {
					lowPhi, phiWidth := poly.AngleRange()
					                    
					lowIdx, idxWidth := AngleBinRange(
						lowPhi, phiWidth, len(vecs),
					)

					for idx := lowIdx; idx < lowIdx + idxWidth; idx++ {
						m++
						j := idx
						if j >= len(vecs) { j -= len(vecs) }
						w.IntersectionDistance(pt, t, &vecs[j])
					}
				}
			}
		}
	}
	println(m)
}

func BenchmarkSheetRingLineSolve(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }

	y01, m1 := 2.0, -7.0
	y02, m2 := 5.0, 1.0

	valid := make([]bool, len(ts))
	rMaxSqr, rMinSqr := rMax*rMax, rMin*rMin
	_ = rMinSqr
	for i := range ts {
		for j := 0; j < 4; j++ {
			x, y, z := ts[i][j][0]-hx, ts[i][j][1]-hy, ts[i][j][2]-hz
			rSqr := x*x + y*y + z*z
			if rSqr < rMaxSqr {
				valid[i] = true
				break
			}
		}
	}

	b.ResetTimer()

	t := new(Tetra)
	pt := new(PluckerTetra)
	poly := new(TetraSlice)

	rot := EulerMatrix(0, 0, 0)
	dr := &Vec{-hx, -hy, -hz}
	vecs, _ := vectorRing(1000)

	for n := 0; n < b.N; n++ {
		for i, ok := range valid {
			if ok {
				*t = ts[i]
				t.Translate(dr)
				t.Rotate(rot)
				pt.Init(t)
				ok := t.ZPlaneSlice(pt, 0, poly)
				if ok {
					lowPhi, phiWidth := poly.AngleRange()
					                    
					lowIdx, idxWidth := AngleBinRange(
						lowPhi, phiWidth, len(vecs),
					)

					for idx := lowIdx; idx < lowIdx + idxWidth; idx++ {
						j := idx
						if j >= len(vecs) { j -= len(vecs) }
						x := (y02 - y01) / (m2 - m1)
						y := y01  + x * m1
						math.Sqrt(x*x + y*y)
					}
				}
			}
		}
	}
}

func BenchmarkIntersectionBary(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
		pts[i].Init(&ts[i])
	}

	w := new(IntersectionWorkspace)
	ap := randomAnchoredPluckerVec()

	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		w.IntersectionBary(&pts[i], &ap.PluckerVec)
	}
}

func BenchmarkIntersectionDistance(b *testing.B) {
	ts := make([]Tetra, 1<<10)
	pts := make([]PluckerTetra, len(ts))
	for i := range ts {
		randomizeTetra(&ts[i], 0, 1)
		pts[i].Init(&ts[i])
	}

	w := new(IntersectionWorkspace)
	ap := randomAnchoredPluckerVec()

	b.ResetTimer()

	for n := 0; n < b.N; n++ {
		i := n % len(ts)
		w.IntersectionDistance(&pts[i], &ts[i], ap)
	}
}

func almostEq(x1, x2, eps float32) bool {
	return x1 + eps > x2 && x1 - eps < x2
}

func TestIntersectionDistance(t *testing.T) {
	P, L := Vec{-0.5, 1, 0.5}, Vec{1, 0, 0}
	eps := float32(1e-4)

	ap := new(AnchoredPluckerVec)
	pt := new(PluckerTetra)
	w := new(IntersectionWorkspace)

	ap.Init(&P, &L)

	table := []struct{
		t Tetra
		enter, exit float32 
		ok bool
	}{	
		{Tetra{Vec{1,0,4},Vec{1,4,0},Vec{1,0,0},Vec{5,0,0}},1.5,4,true},
		{Tetra{Vec{1,4,0},Vec{1,0,4},Vec{3,0,0},Vec{5, 0, 0}},2.75,4,true},
		{Tetra{Vec{1,0,4},Vec{1,4,0},Vec{1,0,0},Vec{5,0,0}},1.5,4,true},

		{Tetra{Vec{-1,4,0},Vec{-1,0,4},Vec{-1,0,0},Vec{3,0,0}},-0.5,2,true},
		{Tetra{Vec{9,4,0},Vec{9,0,4},Vec{9,0,0},Vec{13,0,0}},9.5,12,true},


		{Tetra{Vec{1,0,4},Vec{1,0,0},Vec{5,0,0},Vec{1,4,0}},1.5,4,true},
		{Tetra{Vec{1,0,0},Vec{5,0,0},Vec{1,4,0},Vec{1,0,4}},1.5,4,true},
		{Tetra{Vec{5,0,0},Vec{1,4,0},Vec{1,0,4},Vec{1,0,0},},1.5,4,true},
		
		{Tetra{Vec{1,6,0},Vec{1,2,4},Vec{1,2,0},Vec{5,2,0}},0,0,false},

	}

	for i, test := range table {
		test.t.Orient(+1)
		pt.Init(&test.t)
		enter, exit, ok := w.IntersectionDistance(pt, &test.t, ap)
		if ok != test.ok {
			t.Errorf("%d) Expected ok = %v, but got %v.", i+1, test.ok, ok)
		} else if !almostEq(enter, test.enter, eps) {
			t.Errorf(
				"%d) Expected enter = %g, but got %g", i+1, test.enter, enter,
			)
		} else if !almostEq(exit, test.exit, eps) {
			t.Errorf(
				"%d) Expected leave = %g, but got %g", i+1, test.exit, exit,
			)
		}
	}
}

func almostContains(x, y float32, xs, ys []float32) bool {
	for i := range xs {
		if almostEq(xs[i], x, 1e-4) && almostEq(ys[i], y, 1e-4) { return true }
	}
	return false
}

func almostEqTetraSlice(poly *TetraSlice, xs, ys []float32) bool {
	if len(xs) != poly.Points { return false }
	for i := range xs {
		if !almostContains(poly.Xs[i], poly.Ys[i], xs, ys) {
			return false
		}
	}
	return true
}

func TestZPlaneSlice(t *testing.T) {
	tet := Tetra{Vec{0, 4, 3}, Vec{0, 4, -1}, Vec{4, 4, -1}, Vec{0, 8, -1}}
	xs, ys := []float32{0, 0, 3}, []float32{7, 4, 4}

	tet.Orient(+1)
	pt := new(PluckerTetra)
	pt.Init(&tet)
	poly := new(TetraSlice)

	if !tet.ZPlaneSlice(pt, 0, poly) {
		t.Errorf("z=0 did not slice a tetrahedron it intersected.")
	} else if !almostEqTetraSlice(poly, xs, ys) {
		t.Errorf(
			"Expected xs = %v, ys = %v, but got xs = %v, ys = %v.",
			xs, ys, poly.Xs, poly.Ys,
		)
	}


	tet.Translate(&Vec{0, 0, 10})
	pt.Init(&tet)
	if tet.ZPlaneSlice(pt, 0, poly) {
		t.Errorf("z=0 sliced a tetrahedron it did not intersect.")
	}
}

func vectorRing(n int) (vecs []AnchoredPluckerVec, phis []float32) {
	phis = make([]float32, n)
	vecs = make([]AnchoredPluckerVec, n)
	P := Vec{0, 0, 0}

	for i := 0; i < n; i++ {
		phi := float32(i) * math.Pi * 2 / float32(n)
		phis[i] = phi
		x := float32(math.Cos(float64(phi)))
		y := float32(math.Sin(float64(phi)))
		L := Vec{x, y, 0}
		vecs[i].Init(&P, &L)
	}

	return vecs, phis
}

func tetXs(t *Tetra) []float32 {
	xs := make([]float32, 4)
	for i := 0; i < 4; i++ { xs[i] = t[i][0] }
	return xs
}

func tetYs(t *Tetra) []float32 {
	ys := make([]float32, 4)
	for i := 0; i < 4; i++ { ys[i] = t[i][1] }
	return ys
}

func TestAngleRange(t *testing.T) {
	n, m := 100, 1000
	vecs, phis := vectorRing(n)
	tet := new(Tetra)
	poly := new(TetraSlice)
	pt := new(PluckerTetra)
	w := new(IntersectionWorkspace)

	for i := 0; i < m; i++ {
		randomizeTetra(tet, -1, +1)
		pt.Init(tet)

		if !tet.ZPlaneSlice(pt, 0, poly) { continue }
		low, width := poly.AngleRange()

		for j, phi := range phis {

			_, exit, ok := w.IntersectionDistance(pt, tet, &vecs[j])
			ok = ok && exit > 0
			inRange := AngleInRange(phi, low, width)
			if inRange && !ok {
				t.Errorf("Expected vec %d to be in range for tet %d.",
					j+1, i+1)
			} else if !inRange && ok {
				t.Errorf("Expected vec %d to not be in range for tet %d.",
					j+1, i+1)
			}
		}
	}
}

func BenchmarkSolveLine(b *testing.B) {
	y01, m1 := 2.0, -7.0
	y02, m2 := 5.0, 1.0
	for i := 0; i < b.N; i++ {
		x := (y02 - y01) / (m2 - m1)
		y := y01  + x * m1
		math.Sqrt(x*x + y*y)
	}
}

func BenchmarkSphereLineSegmentIntersection(b *testing.B) {
	n := 1000
	ls := make([]LineSegment, n)
	for i := range ls {
		ls[i] = LineSegment{ Vec{ float32(rand.Float64()),
			float32(rand.Float64()),
			float32(rand.Float64())}, Vec{0, 0, 1}, 0, 1 }
	}
	// Want most of the lines to intersect.
	sphere := Sphere{ Vec{0.5, 0.5, 0.5}, 0.5 }

	b.ResetTimer()
	idx := 0
	for i := 0; i < b.N; i++ {
		sphere.LineSegmentIntersect(&ls[idx])
		idx++
		if idx == n { idx = 0 }
	}
}

func coords(idx, cells int64) (x, y, z int64) {
	x = idx % cells
	y = (idx % (cells * cells)) / cells
	z = idx / (cells * cells)
	return x, y, z
}

func index(x, y, z, cells int64) int64 {
	return x + y * cells + z * cells * cells
}

func readTetra(idxs *rGeom.TetraIdxs, xs []rGeom.Vec, t *Tetra) {
	for i := 0; i < 4; i++ {
		t[i] = Vec(xs[idxs[i]])
	}
}

func myMain() int {
	file := "/project/surph/mansfield/data/sheet_segments/" + 
		"Box_L0063_N1024_G0008_CBol/snapdir_100/sheet167.dat"
	if err := io.ReadSheetHeaderAt(file, &hd); err != nil {
		return 1
	}
	xs = make([]rGeom.Vec, hd.GridCount)
	if err := io.ReadSheetPositionsAt(file, xs); err != nil {
		panic(err.Error())
	}

	n := hd.SegmentWidth * hd.SegmentWidth * hd.SegmentWidth
	ts = make([]Tetra, n * 6)
	pts = make([]PluckerTetra, n * 6)

	idxBuf := &rGeom.TetraIdxs{}
	for writeIdx := int64(0); writeIdx < n; writeIdx++ {
		x, y, z := coords(writeIdx, hd.SegmentWidth)
		readIdx := index(x, y, z, hd.SegmentWidth)

		for dir := int64(0); dir < 6; dir++ {
			tIdx := 6 * writeIdx + dir
			idxBuf.Init(readIdx, hd.GridWidth, 1, int(dir))
			readTetra(idxBuf, xs, &ts[tIdx])
			ts[tIdx].Orient(+1)
		}
	}

	return 0
}
