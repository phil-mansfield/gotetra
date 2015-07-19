package io

import (
	"encoding/binary"
	"fmt"
	"io"
	"math"
	"os"

	"github.com/phil-mansfield/gotetra/cosmo"
	"github.com/phil-mansfield/gotetra/density"

	"unsafe"
)

var (
	end = binary.LittleEndian
)
const (
	Version = uint64(2)
)

type GridHeader struct {
	EndiannessVersion uint64
	Type TypeInfo
    Cosmo CosmoInfo
    Render RenderInfo
	Loc LocationInfo
	Vel LocationInfo
}

type TypeInfo struct {
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

func ReadGridHeader(fname string) (*GridHeader, error) {
	f, err := os.Open(fname)
	defer f.Close()
	if err != nil { return nil, err }
	hd := &GridHeader{}
	err = binary.Read(f, end, hd)
	if err != nil { return nil, err }
	return hd, nil
}

func ReadGrid(fname string) ([]float64, error) {
	f, err := os.Open(fname)
	defer f.Close()
	if err != nil { return nil, err }
	hd := &GridHeader{}
	err = binary.Read(f, end, hd)
	if err != nil { return nil, err }

	if hd.Type.IsVectorGrid == 1 {
		return nil, fmt.Errorf("io.ReadGrid() can only read scalar grids.")
	}

	sp := hd.Loc.PixelSpan
	val32s := make([]float32, sp[0] * sp[1] * sp[2])
	err = binary.Read(f, end, val32s)
	if err != nil { return nil, err }

	vals := make([]float64, len(val32s))
	for i, x := range val32s { vals[i] = float64(x) }
	return vals, nil
}

type Vector [3]float64
type IntVector [3]int64

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

func WriteBuffer(
	buf density.Buffer,
	cosmo CosmoInfo, render RenderInfo, loc LocationInfo,
	wr io.Writer,
) {
	hd := GridHeader{}
	hd.EndiannessVersion = EndiannessVersionFlag(end)
	hd.Type.HeaderSize = int64(unsafe.Sizeof(hd))
	hd.Type.GridType = int64(buf.Quantity())

	hd.Cosmo = cosmo
	hd.Render = render
	hd.Loc = loc

	binary.Write(wr, end, &hd)
	if xs, ok := buf.FinalizedScalarBuffer(); ok {
		hd.Type.IsVectorGrid = 0
		binary.Write(wr, end, xs)
	} else if xs, ys, zs, ok := buf.FinalizedVectorBuffer(); ok {
		hd.Type.IsVectorGrid = 1
		binary.Write(wr, end, xs)
		binary.Write(wr, end, ys)
		binary.Write(wr, end, zs)
	} else {
		panic("Buffer is neither scalar nor vector.")
	}
}

func EndiannessVersionFlag(end binary.ByteOrder) uint64 {
	// If the version number ever gets larger than IntMax32, we have bigger
	// problems on our hands.

	flag := reverse(Version-1) | (Version-1)
	if end == binary.LittleEndian { return ^flag }
	return flag
}

func reverse(x uint64) uint64 {
	out := uint64(0)
	for i := uint(0); i < 64; i++ {
		if 1 & (x >> i) == 1 { out |= 1 << (63 - i) }
	}
	return out
}
