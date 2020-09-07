package main
import (
	"flag"
	"fmt"
	"path"
	"encoding/binary"
	"strings"
	"math"
	"log"
	"runtime"
	"runtime/pprof"
	"os"
	"io/ioutil"

	"gopkg.in/gcfg.v1"

	"github.com/phil-mansfield/gotetra/render"
	"github.com/phil-mansfield/gotetra/render/density"
	"github.com/phil-mansfield/gotetra/render/geom"
	"github.com/phil-mansfield/gotetra/render/io"

	"unsafe"
)

const (
	catalogBufLen = 1<<12
)

// FileGroup contains utility files for logging and writing profiles to.
type FileGroup struct {
	log, prof *os.File
}

// Close closes the files inside FileGroup.
func (fg *FileGroup) Close() {
	if fg.log != nil {
		err := fg.log.Close()
		if err != nil { log.Fatal(err.Error()) }
	}

	if fg.prof != nil {
		pprof.StopCPUProfile()
		err := fg.prof.Close()
		if err != nil { log.Fatal(err.Error()) }
	}
}

var (
	// All Gadget files are little endian.
	gadgetEndianness = binary.LittleEndian
)

func main() {
	// The main function manages input sanitization and calls the secondary
	// main functions for each mode. The code tries to fail gracefully if the
	// user provides incorrect input. This takes a lot more code than you would
	// think.

	var (
		renderStr, convertSnapshot string
		exampleConfig string
	)
	vars := map[string]*string {
		"Render": &renderStr,
		"ConvertSnapshot": &convertSnapshot,
		"ExampleConfig": &exampleConfig,
	}

	flag.IntVar(
		&render.NumCores, "Threads", runtime.NumCPU(),
		"Number of threads used. Default is the number of logical cores.",
	)
	flag.StringVar(
		&renderStr, "Render", "",
		"Configuration file for [Render] mode, along with at least one " + 
			"Bounds file that specifies the render range.",
	)
	flag.StringVar(
		&convertSnapshot, "ConvertSnapshot", "",
		"Configuration file for [ConvertSnapshot] mode.",
	)
	flag.StringVar(
		&exampleConfig,
		"ExampleConfig", "", "Prints an example configuration file of the " + 
			"specified type to stdout. Accepted arguments are 'Render', " +
			"'ConvertSnapshot', and 'Bounds'.",
	)

	flag.Parse()

	// Figure out the mode and fail with a descriptive error is the user gave 
	// incorrect flags.
	modeName, err := getModeName(vars)
	if err != nil { log.Fatal(err.Error()) }
	
	// Based on the mode being run, do some basic input sanitization and 
	// run the econdary main function.
	switch modeName {
	case "Render":
		wrap := io.DefaultRenderWrapper()
		err := gcfg.ReadFileInto(wrap, renderStr)

		if err != nil { log.Fatal(err.Error()) }
		con := &wrap.Render

		if !con.ValidInput() {
			log.Fatal("Invalid/non-existent 'Input' value.")
		} else if !con.ValidOutput() {
			log.Fatal("Invalid/non-existent 'Output' value.")
		} else if !con.ValidSubsampleLength() {
			log.Fatal("Invalid 'SubsampleLength' value.")
		}

		if !con.ValidImagePixels() && !con.ValidTotalPixels() {
			log.Fatal(
				"You must set either a valid 'ImagePixels' " + 
					"or a valid 'TotalPixels'.",
			)
		} else if !con.ValidParticles() &&
			!con.ValidProjectionDepth() &&
			!con.AutoParticles {
			log.Fatal(
				"You must set either a valid 'Particles' or a valid " + 
					"'ProjectionDepth' or must set 'AutoParticles' to " +
					"true.",
			)
		}

		bounds := flag.Args()
		if len(bounds) < 1 {
			log.Fatal("Must supply at least one bounds file.")
		}
		renderMain(con, bounds)

	case "ConvertSnapshot":
		wrap := io.DefaultConvertSnapshotWrapper()
		err := gcfg.ReadFileInto(wrap, convertSnapshot)
		if err != nil { log.Fatal(err.Error()) }
		con := &wrap.ConvertSnapshot

		if !con.ValidInput() && !con.ValidIteratedInput() {
			log.Fatal("Invalid/non-existent 'Input' value.")
		} else if !con.ValidOutput() && !con.ValidIteratedOutput() {
			log.Fatal("Invalid/non-existent 'Output' value.")
		} else if !con.ValidCells() {
			log.Fatal("Invalid/non-existent 'Cells' value.")
		} else if con.ValidIteratedInput() != con.ValidIteratedOutput() {
			log.Fatal("Only one of IteratedInput and IteratedOutput is set.")
		}

		switch con.InputFormat {
		case "LGadget-2":
			lGadget2Main(con)
		default:
			log.Fatalf("Only LGadget-2 snapshots can be read at this time.")
		}

	case "ExampleConfig":
		switch exampleConfig {
		case "ConvertSnapshot":
			fmt.Println(io.ExampleConvertSnapshotFile)
		case "Render":
			fmt.Println(io.ExampleRenderFile)
		case "Bounds":
			fmt.Println(io.ExampleBoundsFile)
		default:
			log.Fatal(
				"Unrecognized 'ExampleConfig' argument. Only recognized " +
					"arguments are 'Render', 'Bounds', and 'ConvertSnapshot'.",
			)
		}
	default:
		panic("Impossible")
	}
}

// getModeName returns the name of the mode and fails with a descriptive error 
// if the user provided less or more than one mode flag.
func getModeName(vars map[string]*string) (string, error) {
	setNames := []string{}

	for name, varPtr := range vars {
		if *varPtr != "" { setNames = append(setNames, name) }
	}

	if len(setNames) == 0 {
		return "", fmt.Errorf("No flags have been set.")
	}
	
	if len(setNames) > 1 {
		return "", fmt.Errorf(
			"The following flags were set: %s, but gotetra_cmd " + 
				"only accepts one flag at a time.", 
			strings.Join(setNames, ", "),
		)
	}

	return setNames[0], nil
}

// lGadget2Main converts a set of LGadget-2 snapshots to gotetra files
// based on the input config file.
func lGadget2Main(con *io.ConvertSnapshotConfig) {
	if !con.ValidIteratedInput() {
		con.IterationStart = 0
		con.IterationEnd = 0
	}

	// If the user supplied format strings for iterated output, loop over the
	// iteration range. Otherwise, just use the supplied snapshot.
	for i := con.IterationStart; i <= con.IterationEnd; i++ {
		input, output := con.Input, con.Output
		if con.ValidIteratedInput() {
			input = fmt.Sprintf(con.IteratedInput, i)
			output = fmt.Sprintf(con.IteratedOutput, i)
		}

		infos, err := ioutil.ReadDir(input)
		if err != nil { log.Fatal(err.Error()) }
		files := make([]string, len(infos))
		
		for i, info := range infos {
			files[i] = path.Join(input, info.Name())
		}

		// Part 1: read data into memory and put it into a single in-memory
		// grid.
		hd, xs, vs := createGrids(files)

		if err = os.MkdirAll(output, 0777); err != nil {
			log.Fatalf(err.Error())
		}

		// Part 2: write that grid into gtet files.
		writeGrids(output, hd, con.Cells, xs, vs)
	}
}

// createGrids reads reads snapshot data into memory.
func createGrids(
	catalogs []string,
) (hd *io.CatalogHeader, xs, vs []geom.Vec) {
	hs := make([]io.CatalogHeader, len(catalogs))
	for i := range hs {
		hs[i] = *io.ReadGadgetHeader(catalogs[i], gadgetEndianness)
	}

	// Monitor memory usage.
	ms := &runtime.MemStats{}
	runtime.ReadMemStats(ms)
	log.Printf("Allocated %d MB %d B (sys: %d MB)\n",
		ms.Alloc / 1000000, ms.Alloc, ms.Sys / 1000000,
	)
	log.Printf("About to allocate another %d MB\n",
		(hs[0].TotalCount * int64(unsafe.Sizeof(geom.Vec{}))) / 1000000,
	)
	runtime.GC()

	// Allocate position/velocity grids and sinle-file buffers.
	xs = make([]geom.Vec, hs[0].TotalCount)
	vs = make([]geom.Vec, hs[0].TotalCount)
	
	maxLen := int64(0)
	for _, h := range hs {
		if h.Count > maxLen {
			maxLen = h.Count
		}
	}
	idBuf := make([]int64, maxLen)
	xBuf := make([]geom.Vec, maxLen)
	vBuf := make([]geom.Vec, maxLen)
	
	buf := io.NewParticleBuffer(xs, vs, catalogBufLen)

	for i, cat := range catalogs {
		// Read in a single snapsot file.
		if i % 25 == 0 {
			log.Printf("Read %d/%d catalogs", i, len(catalogs))
		}

		N := hs[i].Count

		idBuf = idBuf[0: N]
		xBuf = xBuf[0: N]
		vBuf = vBuf[0: N]

		io.ReadGadgetParticlesAt(
			cat, gadgetEndianness, xBuf, vBuf, idBuf,
		)

		runtime.GC()
		buf.Append(xBuf, vBuf, idBuf)
	}
	buf.Flush()

	if len(catalogs) % 25 != 0 {
		log.Printf("Read %d/%d catalogs", len(catalogs), len(catalogs))
	}

	hs[0].CountWidth = intCubeRoot(hs[0].TotalCount)

	return &hs[0], xs, vs
}

// intCubeRoot returns the cube root of a perfect cube.
func intCubeRoot(x int64) int64 {
	cr := int64(round(math.Pow(float64(x), 1.0 / 3.0)))
	if cr * cr * cr != x { panic("You gave a non-cube to intCubeRoot") }
	return cr
}

// writeGrids writes the in-memory grids to disk as gtet files.
func writeGrids(outDir string, hd *io.CatalogHeader,
	cells int, xs, vs []geom.Vec) {

	log.Println("Writing to directory", outDir)

	segmentWidth := int(hd.CountWidth) / cells
	gridWidth := segmentWidth + 1

	xsSeg := make([]geom.Vec, gridWidth * gridWidth * gridWidth)
	vsSeg := make([]geom.Vec, gridWidth * gridWidth * gridWidth)

	// Write the header
	shd := &io.SheetHeader{}
	shd.Cosmo = hd.Cosmo
	shd.CountWidth = hd.CountWidth
	shd.Count = hd.Count
	shd.Mass = hd.Mass
	shd.TotalWidth = hd.TotalWidth

	shd.SegmentWidth = int64(segmentWidth)
	shd.GridWidth = int64(gridWidth)
	shd.GridCount = int64(shd.GridWidth * shd.GridWidth * shd.GridWidth)
	shd.Cells = int64(cells)

	// Write the data. (x, y, z) indexes over the files.
	for z := int64(0); z < shd.Cells; z++ {
		for y := int64(0); y < shd.Cells; y++ {
			for x := int64(0); x < shd.Cells; x++ {
				copyToSegment(shd, xs, vs, xsSeg, vsSeg)
				file := path.Join(outDir, fmt.Sprintf(
					"sheet%d%d%d.dat", x, y, z,
				))
				io.WriteSheet(file, shd, xsSeg, vsSeg)
				runtime.GC()

				if shd.Idx % 25 == 0 {
					log.Printf("Wrote %d/%d sheet segments.",
						shd.Idx, shd.Cells * shd.Cells * shd.Cells,
					)
				}
				shd.Idx++
			}
		}
	}
	if shd.Idx % 25 != 0 {
		log.Printf("Wrote %d/%d sheet segments.",
			shd.Idx, shd.Cells * shd.Cells * shd.Cells,
		)
	}
}

// boundingBox is the smallest bounding box which contains a group of particles
// in a periodic simulation box. The algorithm I used here is not perfect, but
// avoids needing to sort the particles. It is only garuanteed to work if the
// particles span less than half the width of the box.
type boundingBox struct {
	Width float64
	Center geom.Vec
	ToMax, ToMin, ToPt geom.Vec
}

// Init greates a boundingBox with an initial center, pt, within a box of a
// given width.
func (box *boundingBox) Init(pt *geom.Vec, width float64) {
	box.Width = width
	box.Center = *pt
}

// Add inserts a point into the boundingBox.
func (box *boundingBox) Add(pt *geom.Vec) {
	pt.SubAt(&box.Center, box.Width, &box.ToPt)

	for i := 0; i < 3; i++ {
		box.ToMin[i], box.ToMax[i] = minMax(
			box.ToMin[i], box.ToMax[i], box.ToPt[i],
		)
	}

	//ToPt is now a buffer for the neccesary shift in the center.
	// I don't really understand how this works or why this is so complicated.
	box.ToMax.AddAt(&box.ToMin, &box.ToPt)
	box.ToPt.ScaleSelf(0.5)
	box.Center.AddSelf(&box.ToPt)
    box.ToMax.SubSelf(&box.ToPt, box.Width)
    box.ToMin.SubSelf(&box.ToPt, box.Width)
}

// minMax returns the minimum and maximum of (min, max, x), given the min < max.
func minMax(min, max, x float32) (outMin, outMax float32) {
	if x > max {
		return min, x
	} else if x < min {
		return x, max
	} else {
		return min, max
	}
}

// copyToSegment copies the x and v values in xs and vs into Lagrangian sheet
// segments, xsSeg and vsSeg. Information about the sheet is given by
// SheetHeader.
func copyToSegment(shd *io.SheetHeader, xs, vs, xsSeg, vsSeg []geom.Vec) {
	xStart := shd.SegmentWidth * (shd.Idx % shd.Cells)
	yStart := shd.SegmentWidth * ((shd.Idx / shd.Cells) % shd.Cells)
	zStart := shd.SegmentWidth * (shd.Idx / (shd.Cells * shd.Cells))

	N, N2 := shd.CountWidth, shd.CountWidth * shd.CountWidth

	box := &boundingBox{}
	box.Init(&xs[xStart + N * yStart + N2 * zStart], shd.TotalWidth)

	// smallidx is the index within the segment.
	smallIdx := 0
	vMin := vs[0]
	vMax := vs[0]

	for z := zStart; z < zStart + shd.GridWidth; z++ {
		zIdx := z
		if zIdx == shd.CountWidth { zIdx = 0 }
		for y := yStart; y < yStart + shd.GridWidth; y++ {
			yIdx := y
			if yIdx == shd.CountWidth { yIdx = 0 }
			for x := xStart; x < xStart + shd.GridWidth; x++ {
				xIdx := x
				if xIdx == shd.CountWidth { xIdx = 0 }

				// largeIdx is the index within the overall sheet.
				largeIdx := xIdx + yIdx * N + zIdx * N2
				
				// Copy.
				xsSeg[smallIdx] = xs[largeIdx]
				vsSeg[smallIdx] = vs[largeIdx]

				// Update x and v bounding boxes.
				box.Add(&xsSeg[smallIdx])
				for dim, v := range vs[largeIdx] {
					if v < vMin[dim] {
						vMin[dim] = v
					} else if v > vMax[dim] {
						vMax[dim] = v
					}
				}

				smallIdx++
			}
		}	
	}
	
	box.Center.AddAt(&box.ToMin, &shd.Origin)
	shd.Origin.ModSelf(shd.TotalWidth)
	box.ToMax.ScaleAt(2.0, &shd.Width)

	shd.VelocityOrigin = vMin
	for dim := range vMax { shd.VelocityWidth[dim] = vMax[dim] - vMin[dim] }
}

// validCellNum returns true if cells is a power of two.
func validCellNum(cells int) bool {
	for cells > 1 { cells /= 2 }
	return cells == 1
}

// renderMain is the main function for rendering density fields.
func renderMain(con *io.RenderConfig, bounds []string) {
	// Get the I/O essentials.
	fileNames, hd, fg := densitySetupIO(con)
	defer fg.Close()

	// Generate bounds files.
	configBoxes := make([]io.BoxConfig, 0)
	for _, boundsFile := range bounds {
		boxes, err := io.ReadBoundsConfig(boundsFile, hd.TotalWidth)
		if err != nil { log.Fatal(err.Error()) }
		configBoxes = append(configBoxes, boxes...)
	}

	// Figure out which quantity we're rendering.
	q, ok := density.QuantityFromString(con.Quantity)
	if !ok {
		log.Fatalf("Invalid quantity, '%s'", con.Quantity)
	}
	
	// Figure out cell sizes, pixel sizes, and particle counts.
	boxes := make([]render.Box, len(configBoxes))
	for i := range boxes {
		cells := totalPixels(con, &configBoxes[i], hd.TotalWidth)
		pts := particles(con, &configBoxes[i], hd.TotalWidth)
		boxes[i] = render.NewBox(
			hd.TotalWidth, pts, cells, q, &configBoxes[i],
		)
		log.Println(
			"Rendering to box:", boxes[i].CellSpan(),
			"pixels,", pts, "particles per tetrahedron",
		)
	}

	// Interpolate.
	man, err := render.NewManager(fileNames, boxes, true, q)
	if err != nil { log.Fatal(err.Error()) }
	man.Subsample(con.SubsampleLength)
	man.RenderDensity() // This actually means "render quantity"

	// Write output.
	for i, cBox := range configBoxes {
		box := boxes[i]

		out := path.Join(con.Output, fmt.Sprintf("%s%s%s.gtet",
			con.PrependName, cBox.Name, con.AppendName))

		log.Printf("Writing to %s", out)
		f, err := os.Create(out)
		defer f.Close()
		if err != nil { log.Fatalf("Could not create %s.", out) }
		
		loc := io.NewLocationInfo(
			box.CellOrigin(), box.CellSpan(), box.CellWidth(),
		)
		cos := io.NewCosmoInfo(
			hd.Cosmo.H100 * 100, hd.Cosmo.OmegaM,
			hd.Cosmo.OmegaL, hd.Cosmo.Z, hd.TotalWidth,
		)

		renderInfo := io.NewRenderInfo(
			con.Particles, con.TotalPixels, con.SubsampleLength,
			cBox.ProjectionAxis,
		)

		io.WriteBuffer(box.Vals(), cos, renderInfo, loc, f)
	}
}

// toFloat32 converts a float64 array to a float32 array.
func toFloat32(xs []float64) []float32 {
	ys := make([]float32, len(xs))
	for i := range xs { ys[i] = float32(xs[i]) }
	return ys
}

// densitySetupIO sets up the I/O for a routine which reads gotetra sheet files.
// It gets the relevant filenames, a header for the sheet (i.e. the header of
// an arbitrary file to be read), and a FileGroup struct which handles logging
// and profiling. The filenames have not been filtered by location yet.
func densitySetupIO(con *io.RenderConfig) (
	files []string,
	hd *io.SheetHeader,
	fg *FileGroup,
) {
	var err error
	fg = new(FileGroup)
	hd = new(io.SheetHeader)

	// Set up log file.
	if con.ValidLogFile() {
		fg.log, err = os.Create(con.LogFile)
		if err != nil { log.Fatal(err.Error()) }
		log.SetOutput(fg.log)
	}

	// Set up profile file.
	if con.ValidProfileFile() {
		fg.prof, err = os.Create(con.ProfileFile)
		if err != nil { log.Fatal(err.Error()) }
		err = pprof.StartCPUProfile(fg.prof)
		if err != nil { log.Fatal(err.Error()) }
	}

	// Get file names.
	infos, err := ioutil.ReadDir(con.Input)
	if err != nil { log.Fatal(err.Error()) }
	files = make([]string, len(infos))
	for i := range infos { files[i] = path.Join(con.Input, infos[i].Name()) }

	io.ReadSheetHeaderAt(files[0], hd)

	return files, hd, fg
}

// round rounds x to the nearest integer
func round(x float64) int {
	floor, ceil := math.Floor(x), math.Ceil(x)
	if ceil - x > x - floor {
		return int(ceil)
	} else {
		return int(floor)
	}
}

// totalPixels counts the number of pixels implied by the 
func totalPixels(
	con *io.RenderConfig, box *io.BoxConfig, boxWidth float64,
) int {
	if con.ValidImagePixels() {
		w := maxWidth(box)
		if w > boxWidth {
			log.Fatalf(
				"Requested dimensions of '%s' are larger than the " + 
					"simulation box.", box.Name,
			) 
		}

		return round(boxWidth / w * float64(con.ImagePixels))
	}
	return con.TotalPixels
}


func maxWidth(box *io.BoxConfig) float64 {
	var max float64
	if box.XWidth > box.YWidth {
		max = box.XWidth
	} else {
		max = box.YWidth
	}

	if max > box.ZWidth {
		return max
	} else {
		return box.ZWidth
	}
}

func particles(con *io.RenderConfig, box *io.BoxConfig, boxWidth float64) int {

	if con.AutoParticles {
		cells := totalPixels(con, box, boxWidth)
		cellWidth := boxWidth / float64(cells)

		if box.ProjectionAxis == "X" {
			con.ProjectionDepth = int(math.Ceil(box.XWidth / cellWidth))
		} else if box.ProjectionAxis == "Y" {
			con.ProjectionDepth = int(math.Ceil(box.YWidth / cellWidth))
		} else if box.ProjectionAxis == "Z" {
			con.ProjectionDepth = int(math.Ceil(box.ZWidth / cellWidth))
		} else {
			con.ProjectionDepth = 1
		}
	}

	if con.ValidProjectionDepth() {
		// We're garuanteed that con.ValidImagePixels() is also true.
		pixels := totalPixels(con, box, boxWidth)

		refPixels := 500.0
		refParticles := 5.0
		refDepth := 1.0
		refSubsample := 1.0
		return int(math.Ceil(
			refParticles * math.Pow(float64(pixels) / refPixels, 3) * 
			math.Pow(float64(con.SubsampleLength) / refSubsample, 3) *
			(refDepth / float64(con.ProjectionDepth)),
		))
	}

	return con.Particles
}
