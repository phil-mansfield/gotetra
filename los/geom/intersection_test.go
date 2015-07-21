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

func BenchmarkPluckerTetraInitSheet(b *testing.B) {
	if mainSuccess == 1 { b.FailNow() }
	for i := 0; i < b.N; i++ {
		for idx := range ts {
			pts[idx].Init(&ts[idx])
		}
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
		{Tetra{Vec{1,4,0},Vec{1,0,4},Vec{1,0,0},Vec{5,0,0}},1.5,2.5,true},
		{Tetra{Vec{1,0,4},Vec{1,4,0},Vec{1,0,0},Vec{5,0,0}},1.5,2.5,true},

		{Tetra{Vec{-1,4,0},Vec{-1,0,4},Vec{-1,0,0},Vec{3,0,0}},-0.5,0.5,true},
		{Tetra{Vec{9,4,0},Vec{9,0,4},Vec{9,0,0},Vec{13,0,0}},9.5,10.5,true},


		{Tetra{Vec{1,0,4},Vec{1,0,0},Vec{5,0,0},Vec{1,4,0}},1.5,2.5,true},
		{Tetra{Vec{1,0,0},Vec{5,0,0},Vec{1,4,0},Vec{1,0,4}},1.5,2.5,true},
		{Tetra{Vec{5,0,0},Vec{1,4,0},Vec{1,0,4},Vec{1,0,0},},1.5,2.5,true},
		
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
