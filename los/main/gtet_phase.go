package main

import (
	"flag"
	"log"
	"math"
	"runtime"
	
	"github.com/phil-mansfield/gotetra/render/io"
	"github.com/phil-mansfield/gotetra/los/geom"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/halo"
)

func main() {
	log.Println("gtet_phase")
	
	p := parseCmd()
	ids, snaps, _, err := util.ParseStdin()
	if err != nil { log.Fatal(err.Error()) }

	bounds, grid, err := phaseBounds(ids, snaps, p)
	if err != nil { log.Fatal(err.Error()) }
	
	printPhases(ids, snaps, bounds, grid)
}

type Param struct {
	R200mMult, VMaxMult float64
	VWidth, RWidth int
	SubsampleLength int
}

func parseCmd() *Param {
	p := &Param{}
	flag.Float64Var(&p.R200mMult, "R200mMult", 3.0, "Multiplier to R200m")
	flag.Float64Var(&p.VMaxMult, "VMaxMult", 2.5, "Multiplier to VMax")
	flag.IntVar(&p.VWidth, "VWidth", 100, "Number of velocity bins.")
	flag.IntVar(&p.RWidth, "RWidth", 100, "Number of radial bins.")
	flag.IntVar(&p.SubsampleLength, "SubsampleLength", 1,
		"Subsample length. Must be 2^n.")
	flag.Parse()

	return p
}

type Bound struct {
	VLow, VHigh, RLow, RHigh float64
	VWidth, RWidth int
	SubsampleLength int
	x, y, z, vx, vy, vz float64
}

func phaseBounds(ids, snaps []int, p *Param) ([]Bound, [][]int, error) {
	bs := make([]Bound, len(ids))
	grids := make([][]int, len(ids))

	var xs, vs []rgeom.Vec

	snapBins, idxBins := binBySnap(snaps, ids)
	for _, snap := range snaps {
		runtime.GC()
		
		snapIDs := snapBins[snap]
		idxs := idxBins[snap]

		hds, files, err := util.ReadHeaders(snap)
		n := hds[0].GridWidth*hds[0].GridWidth*hds[0].GridWidth
		if len(xs) == 0 { xs = make([]rgeom.Vec, n) }
		if len(vs) == 0 { vs = make([]rgeom.Vec, n) }


		if err != nil { log.Fatal(err.Error()) }
		vals, err := util.ReadRockstar(
			snap, snapIDs, halo.X, halo.Y, halo.Z, halo.Rad200b,
			halo.Vx, halo.Vy, halo.Vz, halo.VMax,
		)
		if err != nil { log.Fatal(err.Error()) }

		ss := []geom.Sphere{}
		for i := range snapIDs {
			x, y, z, r := vals[0][i], vals[1][i], vals[2][i], vals[3][i]
			vx, vy, vz, vmax := vals[4][i], vals[5][i], vals[6][i], vals[7][i]
			
			b := &bs[idxs[i]]
			b.VHigh = vmax * p.VMaxMult
			b.VLow = -vmax * p.VMaxMult
			b.RHigh = r * p.R200mMult
			b.RLow = 0
			b.VWidth = p.VWidth
			b.RWidth = p.RWidth

			b.x, b.y, b.z = x, y, z
			b.vx, b.vy, b.vz = vx, vy, vz
			s := geom.Sphere{}
			s.R = float32(r) * float32(p.R200mMult)
			s.C = geom.Vec{float32(x), float32(y), float32(z)}
			
			grids[idxs[i]] = make([]int, p.VWidth * p.RWidth)
			ss = append(ss, s)
		}

		intrBins := binIntersections(hds, ss)
		for i := range hds {
			if len(intrBins[i]) == 0 { continue }
			err := io.ReadSheetPositionsAt(files[i], xs)
			if err != nil { return nil, nil, err }
			err = io.ReadSheetVelocitiesAt(files[i], vs)
			if err != nil { return nil, nil, err }
			
			for _, hi := range intrBins[i] {
				idx := idxs[hi]
				binParticles(
					&hds[i], xs, vs, &bs[idx],
					grids[idx], p.SubsampleLength,
				)
			}
		}
	}

	return bs, grids, nil
}

func binParticles(
	hd *io.SheetHeader, xs, vs []rgeom.Vec, b *Bound, grid []int, skip int,
) {
	dr := (b.RHigh - b.RLow) / float64(b.RWidth)
	dv := (b.VHigh - b.VLow) / float64(b.VWidth)
	incr := skip*skip*skip
	tw := hd.TotalWidth

	sw := int(hd.SegmentWidth)
	for iz := 0; iz < sw; iz += skip {
		for iy := 0; iy < sw; iy += skip {
			for ix := 0; ix < sw; ix += skip {
				i := ix + iy*sw + iz*sw*sw
				xVec, vVec := xs[i], vs[i]
				x, y, z := float64(xVec[0]), float64(xVec[1]), float64(xVec[2])
				x, y, z = x - b.x, y - b.y, z - b.z
				if x > tw / 2 { x -= tw }
				if y > tw / 2 { y -= tw }
				if z > tw / 2 { z -= tw }
				if x < -tw / 2 { x += tw }
				if y < -tw / 2 { y += tw }
				if z < -tw / 2 { z += tw }
				
				r := math.Sqrt(x*x + y*y + z*z)
				ri := int(math.Floor((r - b.RLow) / dr))
				if ri >= b.RWidth || ri < 0 { continue }
				
				vx, vy, vz := float64(vVec[0]),float64(vVec[1]),float64(vVec[2])
				vx, vy, vz = vx - b.vx, vy - b.vy, vz - b.vz
				vrx, vry, vrz := vx * x / r, vy * y / r, vz * z / r
				vr := vrx + vry + vrz
				vri := int(math.Floor((vr - b.VLow) / dv))
				if vri >= b.VWidth || vri < 0 { continue }
			
				grid[vri*b.RWidth + ri] += incr
			}
		}
	}
}

func printPhases(ids, snaps []int, bounds []Bound, grids [][]int) {
	rows := [][]float64{}
	for i := range ids {
		floatGrid := make([]float64, len(grids[i]))
		for j := range grids[i] { floatGrid[j] = float64(grids[i][j]) }
		
		b := bounds[i]
		boundsRow := []float64{
			b.VLow, b.VHigh, b.RLow, b.RHigh,
			float64(b.VWidth), float64(b.RWidth),
		}
		
		row := append(boundsRow, floatGrid...)
		rows = append(rows, row)
	}
	util.PrintRows(ids, snaps, rows)
}

func binBySnap(snaps, ids []int) (snapBins, idxBins map[int][]int) {
	snapBins = make(map[int][]int)
	idxBins = make(map[int][]int)
	for i, snap := range snaps {
		snapBins[snap] = append(snapBins[snap], ids[i])
		idxBins[snap] = append(idxBins[snap], i)
	}
	return snapBins, idxBins

}

func binIntersections(
	hds []io.SheetHeader, spheres []geom.Sphere,
) [][]int {
	bins := make([][]int, len(hds))
	for i := range hds {
		for si := range spheres {
			if sheetIntersect(spheres[si], &hds[i]) {
				bins[i] = append(bins[i], si)
			}
		}
	}
	return bins
}

func sheetIntersect(s geom.Sphere, hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(s.C[0], s.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(s.C[1], s.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(s.C[2], s.R, hd.Origin[2], hd.Width[2], tw)
}

func inRange(x, r, low, width, tw float32) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

func wrapDist(x1, x2, width float32) float32 {
	dist := x1 - x2
	if dist > width / 2 {
		return dist - width
	} else if dist < width / -2 {
		return dist + width
	} else {
		return dist
	}
}
