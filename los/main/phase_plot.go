package main

import (
	"log"
	"os"
	"strconv"
	"strings"
	"io/ioutil"
	"math"
	"fmt"

	"github.com/phil-mansfield/gotetra/los/geom"
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/halo"
	"github.com/phil-mansfield/gotetra/render/io"
	util "github.com/phil-mansfield/gotetra/los/main/gtet_util"
)

const G = 6.144e-9

func main() {
	ids, snaps, err := parseStdin()
	if err != nil { log.Fatal(err.Error()) }

	// Screw it, I don't care about speed.
	for i, snap := range snaps {
		hds, files, err := util.ReadHeaders(snap)
		if err != nil { log.Fatal(err.Error()) }
		vals, err := util.ReadRockstar(
			snap, ids[i:i+1], halo.X, halo.Y, halo.Z, halo.Rad200b,
			halo.Vx, halo.Vy, halo.Vz, halo.VMax,
		)
		if err != nil { log.Fatal(err.Error()) }
		x, y, z, r := vals[0][0], vals[1][0], vals[2][0], vals[3][0]
		vx, vy, vz, vmax := vals[4][0], vals[5][0], vals[6][0], vals[7][0]
		
		s := &geom.Sphere{}
		s.R = float32(r) * 3
		s.C = geom.Vec{float32(x), float32(y), float32(z)}
		v := geom.Vec{float32(vx), float32(vy), float32(vz)}

		vBins, rBins := 100, 100
		vLim := 2.5 * vmax
		ns := make([]int, vBins * rBins)

		xs := make([]rgeom.Vec, hds[0].GridCount)
		vs := make([]rgeom.Vec, hds[0].GridCount)
		for j := range hds {
			if !sheetIntersect(s, &hds[j]) { continue }
			err = io.ReadSheetPositionsAt(files[j], xs)
			if err != nil { log.Fatal(err.Error()) }
			err = io.ReadSheetVelocitiesAt(files[j], vs)
			if err != nil { log.Fatal(err.Error()) }

			bin(xs, vs, ns, rBins, vBins, s, v, vLim)
		}
		writeBins(ns, snap, ids[i], rBins, vBins)
	}
}

func bin(xs, vs []rgeom.Vec, ns []int, rBins, vBins int,
	s *geom.Sphere, v geom.Vec, vLim float64) {
	dr, dv := s.R / float32(rBins), float32(vLim / float64(vBins / 2))
	minv, maxv := float32(0.0), float32(0.0)
	minvx, minvy, minvz := vs[0][0], vs[0][1], vs[0][2]
	maxvx, maxvy, maxvz := vs[0][0], vs[0][1], vs[0][2]

	for i := range xs {
		dx, dy, dz := xs[i][0] - s.C[0], xs[i][1] - s.C[1], xs[i][2] - s.C[2]
		vx, vy, vz := vs[i][0] - v[0], vs[i][1] - v[1], vs[i][2] - v[2]
		r := float32(math.Sqrt(float64(dx*dx + dy*dy + dz*dz)))
		vr := (dx*vx + dy*vy + dz*vz) / r

		if vr > maxv {
			maxv = vr
		} else if vr < minv {
			minv = vr
		}

		if vs[i][0] > maxvx {
			maxvx = vs[i][0]
		} else if vs[i][0] < minvx {
			minvx = vs[i][0]
		}

		if vs[i][1] > maxvy {
			maxvy = vs[i][1]
		} else if vs[i][1] < minvy {
			minvy = vs[i][1]
		}
		
		if vs[i][2] > maxvz {
			maxvz = vs[i][2]
		} else if vs[i][2] < minvz {
			minvz = vs[i][2]
		}
		
		ri := int(r / dr)
		vi := int(vr / dv) + (vBins / 2)
		
		if ri >= rBins || vi >= vBins || vi < 0 { continue }
		ns[ri + vi * rBins]++
	}
}

func writeBins(ns []int, snap, id, rBins, vBins int) error {
	fname := fmt.Sprintf("s%d_id%d.dat", snap, id)
	f, err := os.Create(fname)
	defer f.Close()
	if err != nil { return err }

	bs := []byte{}
	for row := 0; row < vBins; row++ {
		for i := 0; i < rBins; i++ {
			bs = strconv.AppendInt(bs, int64(ns[i + row*rBins]), 10)
			bs = append(bs, ' ')
		}
		bs = append(bs, '\n')
	}
	_, err = f.Write(bs)
	return err
}

func parseStdin() (ids, snaps []int, err error) {
	ids, snaps = []int{}, []int{}
	lines, err := stdinLines()
	if err != nil { return nil, nil, err }
	for i, line := range lines {
		rawTokens := strings.Split(line, " ")
		tokens := make([]string, 0, len(rawTokens))
		for _, tok := range rawTokens {
			if len(tok) != 0 { tokens = append(tokens, tok) }
		}

		var (
			id, snap int
			err error
		)
		switch len(tokens) {
		case 0:
			continue
		case 2:
			id, err = strconv.Atoi(tokens[0])
			if err != nil {
				return nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[0],
				)
			} 
			snap, err = strconv.Atoi(tokens[1]) 
			if err != nil {
				return nil, nil, fmt.Errorf(
					"One line %d of stdin, %s does not parse as an int.",
					i + 1, tokens[1],
				)
			} 
		case 1:
			if tokens[0] == "" { continue }
			fallthrough
		default:
			return nil, nil, fmt.Errorf(
				"Line %d of stdin has %d tokens, but 2 are required.",
				i + 1, len(tokens),
			)
		}

		ids = append(ids, id)
		snaps = append(snaps, snap)
	}

	return ids, snaps, nil
}
	
func stdinLines() ([]string, error) {
		bs, err := ioutil.ReadAll(os.Stdin)
	if err != nil {
		return nil, fmt.Errorf(
			"Error reading stdin: %s.", err.Error(),
		)
	}

	text := string(bs)
	return strings.Split(text, "\n"), nil
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

func inRange(x, r, low, width, tw float32) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

// SheetIntersect returns true if the given halo and sheet intersect one another
// and false otherwise.
func sheetIntersect(s *geom.Sphere, hd *io.SheetHeader) bool {
	tw := float32(hd.TotalWidth)
	return inRange(s.C[0], s.R, hd.Origin[0], hd.Width[0], tw) &&
		inRange(s.C[1], s.R, hd.Origin[1], hd.Width[1], tw) &&
		inRange(s.C[2], s.R, hd.Origin[2], hd.Width[2], tw)
}
