package loop

import (
	rgeom "github.com/phil-mansfield/gotetra/render/geom"
)

type Sphere struct {
	origin [3]float64
	rMin, rMax, rMin2, rMax2, tw float64
}

func (s *Sphere) Init(origin [3]float64, rMin, rMax float64) {
	s.rMin = rMin
	s.rMax = rMax
	s.rMin2 = rMin * rMin
	s.rMax2 = rMax * rMax
	s.origin = origin
}

func (s *Sphere) Transform(vecs []rgeom.Vec, totalWidth float64) {
    x0 := float32(s.origin[0])
    y0 := float32(s.origin[1])
    z0 := float32(s.origin[2])
    tw := float32(totalWidth)
    tw2 := tw / 2

    for i, vec := range vecs {
		x, y, z := vec[0], vec[1], vec[2]
		dx, dy, dz := x - x0, y - y0, z - z0

        if dx > tw2 {
            vecs[i][0] -= tw
        } else if dx < -tw2 {
            vecs[i][0] += tw
        }

        if dy > tw2 {
            vecs[i][1] -= tw
        } else if dy < -tw2 {
            vecs[i][1] += tw
        }

        if dz > tw2 {
            vecs[i][2] -= tw
        } else if dz < -tw2 {
            vecs[i][2] += tw
        }
	}
}

func (s *Sphere) Contains(x, y, z float64) bool {
    x0, y0, z0 := s.origin[0], s.origin[1], s.origin[2]
    dx, dy, dz := x - x0, y - y0, z - z0
    r2 :=  dx*dx + dy*dy + dz*dz
	return s.rMin2 < r2 && r2 < s.rMax2
}

func (s *Sphere) IntersectBox(origin, span [3]float64, tw float64) bool {
	return inRange(s.origin[0], s.rMax, origin[0], span[0], tw) &&
		inRange(s.origin[1], s.rMax, origin[1], span[1], tw) &&
		inRange(s.origin[2], s.rMax, origin[2], span[2], tw)
}

func inRange(x, r, low, width, tw float64) bool {
	return wrapDist(x, low, tw) > -r && wrapDist(x, low + width, tw) < r
}

func wrapDist(x1, x2, width float64) float64 {
	dist := x1 - x2
	if dist > width / 2 {
		return dist - width
	} else if dist < width / -2 {
		return dist + width
	} else {
		return dist
	}
}

type Box struct {
	origin, span [3]float64
}

func (b *Box) Init(origin, span [3]float64) {
	b.origin = origin
	b.span = span
}

func (b *Box) Transform(vecs []rgeom.Vec, totalWidth float64) {
    x0 := float32(b.origin[0])
    y0 := float32(b.origin[1])
    z0 := float32(b.origin[2])
    tw := float32(totalWidth)
    tw2 := tw / 2

    for i, vec := range vecs {
		x, y, z := vec[0], vec[1], vec[2]
		dx, dy, dz := x - x0, y - y0, z - z0

        if dx > tw2 {
            vecs[i][0] -= tw
        } else if dx < -tw2 {
            vecs[i][0] += tw
        }

        if dy > tw2 {
            vecs[i][1] -= tw
        } else if dy < -tw2 {
            vecs[i][1] += tw
        }

        if dz > tw2 {
            vecs[i][2] -= tw
        } else if dz < -tw2 {
            vecs[i][2] += tw
        }
	}
}

func (b *Box) Contains(x, y, z float64) bool {
	lowX, highX := b.origin[0], b.origin[0] + b.span[0]
	lowY, highY := b.origin[1], b.origin[1] + b.span[1]
	lowZ, highZ := b.origin[2], b.origin[2] + b.span[2]
	return lowX < x && x < highX && 
		lowY < y && y < highY && 
		lowZ < z && z < highZ
}

func (b *Box) IntersectBox(origin, span [3]float64, tw float64) bool {
	s2x := b.span[0] / 2
	s2y := b.span[1] / 2
	s2z := b.span[2] / 2

	return inRange(b.origin[0] + s2x, s2x, origin[0], span[0], tw) &&
		inRange(b.origin[1] + s2y, s2y, origin[1], span[1], tw) &&
		inRange(b.origin[2] + s2z, s2z, origin[2], span[2], tw)
}
