package main

import (
	"encoding/binary"
	"log"
	"io"
	"math"
	"os"

	"github.com/phil-mansfield/gotetra/cosmo"

	"unsafe"
)

var end = binary.LittleEndian

type GridHeader struct {
	Type TypeInfo
    Cosmo CosmoInfo
    Render RenderInfo
	Loc LocationInfo
}

type TypeInfo struct {
    Endianness int64
	HeaderSize int64
    GridType int64
    IsVectorGrid int64
}

type CosmoInfo struct {
    Redshift, ScaleFactor float64
    OmegaM, OmegaL, Hubble float64
    RhoMean, RhoCritical float64
	BoxWidth float64
}

type RenderInfo struct {
    Particles int64
	TotalPixels int64
    SubsampleLength int64
    MinProjectionDepth int64
}

type LocationInfo struct {
    Origin, Span Vector
    PixelOrigin, PixelSpan IntVector
    PixelWidth float64
}

type Vector [3]float64
type IntVector [3]int64

type GridFlag int64
const (
	Density GridFlag = iota
	DensityGradient
	Velocity
	VelocityDivergence
	VelocityCurl
)
	
func NewCosmoInfo(H0, omegaM, omegaL, z, boxWidth float64) CosmoInfo {
	a := 1 / (1 + z)
	rhoC := cosmo.RhoCritical(H0, omegaM, omegaL, z)
	rhoM := cosmo.RhoAverage(H0, omegaM, omegaL, z)
	
	ci := CosmoInfo{z, a, omegaM, omegaL, H0, rhoM, rhoC, boxWidth}
	return ci
}
	
func NewRenderInfo(particles, totalCells, skip int) RenderInfo {
	proj := projectionDepth(particles, totalCells, skip)
	ri := RenderInfo{
		int64(particles), int64(totalCells), int64(skip), int64(proj),
	}
	return ri
}
	
func projectionDepth(particles, totalCells, skip int) int {
	proj := (20000.0 / float64(particles)) *
		math.Pow(float64(totalCells) / 5000, 3) /
		math.Pow(float64(skip), 3)
	
	return int(math.Ceil(proj))
}

func NewLocationInfo(origin, span [3]int, cellWidth float64) LocationInfo {
	loc := LocationInfo{ }

	for i := 0; i < 3; i++ {
		loc.Origin[i] = float64(origin[i]) * cellWidth
		loc.Span[i] = float64(span[i]) * cellWidth
		loc.PixelOrigin[i] = int64(origin[i])
		loc.PixelSpan[i] = int64(span[i])
	}

	loc.PixelWidth = cellWidth

	return loc
}

func WriteDensity(
	rhos []float32, cos CosmoInfo, render RenderInfo,
	loc LocationInfo, wr io.Writer,
) {
	WriteGrid(Density, rhos, cos, render, loc,  wr)
}

func WriteVelocityDivergence(
	divs []float32, cos CosmoInfo, render RenderInfo,
	loc LocationInfo, wr io.Writer,
) {
	WriteGrid(VelocityDivergence, divs, cos, render, loc, wr)
}

func WriteGrid(
	flag GridFlag, xs []float32, cosmo CosmoInfo, render RenderInfo,
	loc LocationInfo, wr io.Writer,
) {
	var endFlag int64
	if end == binary.LittleEndian {
		endFlag = -1
	} else {
		endFlag = 0
	}

	hd := GridHeader{ }
	hd.Type.Endianness = endFlag
	hd.Type.HeaderSize = int64(unsafe.Sizeof(hd))
	hd.Type.GridType = int64(flag)
	if flag != Density && flag != VelocityDivergence {
		hd.Type.IsVectorGrid = 1
	} else {
		hd.Type.IsVectorGrid = 0
	}		

	hd.Cosmo = cosmo
	hd.Render = render
	hd.Loc = loc

	if hd.Loc.PixelSpan[0]*hd.Loc.PixelSpan[1]*hd.Loc.PixelSpan[2] != 
		int64(len(xs)) {

		log.Fatalf(
			"PixelSpan %v is not the same saize as slice of length %d.",
			hd.Loc.PixelSpan, len(xs),
		)
	}

	binary.Write(wr, end, &hd)
	binary.Write(wr, end, xs)
}

func WriteVelocity() { }
func WriteDensityGradient() { }
func WriteVelocityCurl() { }
func WriteVectorGrid() { }

func main() {
	f, err := os.Create(os.Args[1])
	if err != nil { log.Fatal(err.Error()) }
	defer f.Close()
	
	span := [3]int{5, 3, 3}
	origin := [3]int{3, 4, 5}
	cellWidth := 2.0
	loc := NewLocationInfo(origin, span, cellWidth)

	H0, omegaM, omegaL := 0.3, 0.7, 70.0
	z := 0.0
	boxWidth := 100.0
	cos := NewCosmoInfo(H0, omegaM, omegaL, z, boxWidth)

	particles := 20000
	totalCells := 5000
	skip := 1
	render := NewRenderInfo(particles, totalCells, skip)

	rhos := make([]float32, span[0] * span[1] * span[2])
	for i := range rhos { rhos[i] = float32(i) }

	println(
		unsafe.Sizeof(GridHeader{}),
		unsafe.Sizeof(TypeInfo{}),
		unsafe.Sizeof(CosmoInfo{}),
		unsafe.Sizeof(RenderInfo{}),
		unsafe.Sizeof(LocationInfo{}),
	)
	
	WriteDensity(rhos, cos, render, loc, f)
}
