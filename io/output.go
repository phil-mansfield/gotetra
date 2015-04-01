package io

import (
	"encoding/binary"
	"io"
	"math"

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
	ProjectionAxis int64
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
	
func NewRenderInfo(particles, totalCells, skip int, axisStr string) RenderInfo {
	projDepth := projectionDepth(particles, totalCells, skip)

	axis := -1
	if axisStr == "X" { axis = 0 }
	if axisStr == "Y" { axis = 1 }
	if axisStr == "Z" { axis = 2 }

	ri := RenderInfo{
		int64(particles), int64(totalCells), int64(skip),
		int64(projDepth), int64(axis),
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

	binary.Write(wr, end, &hd)
	binary.Write(wr, end, xs)
}

func WriteVelocity() {
	panic("Not yet implemented.")
}

func WriteDensityGradient() {
	panic("Not yet implemented.")
}

func WriteVelocityCurl() {
	panic("Not yet implemented.")
}

func WriteVectorGrid() {
	panic("NOt yet implemented.")
}
