package main

import (
	"encoding/binary"
	"fmt"
	"math"
	"os"
	"strconv"

	"github.com/phil-mansfield/gotetra"
)

func fRange(min, max, x float64) (outMin, outMax float64) {
	if x > max {
		return min, x
	} else if x < min {
		return x, max
	}
	return min, max
}

func main() {
	ps, hd := gotetra.ReadGadget(os.Args[1], binary.LittleEndian)

	boxCount, err := strconv.Atoi(os.Args[2])
    if err != nil { panic(err) }

    boxVolume := math.Pow(hd.BoxSize , 3.0) / float64(boxCount)

	xMin, yMin, zMin := math.MaxFloat64, math.MaxFloat64, math.MaxFloat64
	xMax, yMax, zMax := 0.0, 0.0, 0.0

	for i := range ps {
		x, y, z := ps[i].Xs[0], ps[i].Xs[1], ps[i].Xs[2]
		xMin, xMax = fRange(xMin, xMax, x)
		yMin, yMax = fRange(yMin, yMax, y)
		zMin, zMax = fRange(zMin, zMax, z)
	}

	fmt.Printf("%11s: [%10g %10g]\n", "X Range", xMin, xMax)
	fmt.Printf("%11s: [%10g %10g]\n", "Y Range", yMin, yMax)
	fmt.Printf("%11s: [%10g %10g]\n", "Z Range", zMin, zMax)
	fmt.Printf("%11s:  %10g\n", "Volume", (xMax - xMin) * (yMax - yMin) * (zMax - zMin))
	fmt.Printf("%11s:  %10g\n", "Exp. Volume", boxVolume)
}
