package io

import (
	"fmt"
	"strings"

	"gopkg.in/gcfg.v1"
	"github.com/phil-mansfield/gotetra/render/density"
)


const (
	ExampleConvertSnapshotFile = `[ConvertSnapshot]

#######################
# Required Parameters #
#######################

# Directory containing the input files.
Input = path/to/input/dir
# Directory which output files will be written to.
Output = path/to/output/dir

# The format of the input files. LGadget-2 is currently the only format
# supported. I think vanilla Gadget-2 files should also work here, but keep
# an eye out for any bugs. Support for other file formats will come in the
# future.
InputFormat = LGadget-2

# Specifies the geometry of the output files. It's unlikely that you will want
# to change this.
Cells = 8 

#######################
# Optional Parameters #
#######################

# Gadget2IDSize specifies the size in bits of the IDs used in your Gadget-2
# files. Default is 64.
# Gadget2IDSize = 64

# Output files which are useful for profiling and debugging. Generally, there
# isn't a reason to use these unless something goes wrong.
# ProfileFile = prof.out
# LogFile = log.out

# In some cases, you might want to convert several snapshots at once (e.g.
# sims/output/snapdir001, sims/output/snapdir002, etc.) This can be done with
# IteratedInput and IteratedOutput, which are printf format strings that
# describe your output format (e.g. sims/output/snapdir%03d). You must create
# all the directories which IteratedOutput will use before running, which is
# easiest to do with a Python loop over os.system("mkdir ...") calls.
# IteratedInput = path/to/iterated/input_with_single_%d_format
# IteratedOutput = path/to/iterated/input_with_single_%d_format

# (Inclusive) range for this iteration. If IterationEnd isn't set, folders will
# be iterated through until an invalid one is found.
# IterationStart = 0 
# IterationEnd = 100`
	ExampleRenderFile = `[Render]

# Render specified 

######################
# RequiredParameters #
######################

# Quantity can be set to one of:
# [ Density ]
# In the future, 
# [ Velocity | DensityGradient | VelocityCurl | VelocityDivergence ]
# will be supported.
Quantity = Density

# Directory containing the input files.
Input  = path/to/input/dir
# Directory where the output .gtet files will be written to.
Output = path/to/output/dir

# Default way of specifying pixel size. TotalPixels gives the number of pixels
# across one side of the entire simulation box. This is useful if you are
# rendering many images with various sizes and want a fixed pixel size. You
# do not need to specify this if you use ImagePixels instead.
TotalPixels = 500

# Alternative way of specifying pixel size: the number of pixels required to
# render the longest axis of eaxh boxing box. This is a more natural way to
# specify pixl size if you are only rendering a single image at a time.
# ImagePixels = 100

# Default way of specifying redering resolution. Rendering is performed by
# populating tetrahedra with points, and this variable to controls the number
# of points used. You can think of this as "increasing" the number of particles
# in the simulation by a factor of 6*Particles.
# Expect to rerun the rendering a couple times to get this number right.
Particles = 25

#####################
# OptionalParamters #
#####################

# SubsampleLength allows you to render the image/volume using a subselection of
# the particles in the snapshot files. This is useful for resolution tests and
# some science applications. SubsampleLength must be a power of 2.
# SubsampleLength = 2

# Rendering output files are named after the bounding box. For example, a
# bounding box with the header [Box "halo_1"] will be written to halo_1.gtet.
# You can add leading and ending text to these files names using the following
# two variables (e.g. pre_halo_1_app.gtet). This is useful if you want to
# annotate file names with rendering information, like pixel size or particle
# count.
# PrependName = pre_
# AppendName  = _app

# Output files which are useful for profiling and debugging. Generally, there
# isn't a reason to use these unless something goes wrong.
# ProfileFile = prof.out
# LogFile = log.out

# Use these next two options with caution: I need to do much more testing.

# Alternative way of specifying the particle count. gotetra will (attempt to)
# automatically calculate how many particles are needed so that all projections
# rendered from the resulting grid will have enough particles-per-tetrahedron
# to avoid artifacts.
# AutoParticles = true

# Alternative way of specifying the particle count. Identical to specifying
# AutoParticles, except that projections with a depth below ProjectionDepth
# may contain artifacts. 
# ProjectionDepth = 3`
	ExampleBoxFile = `[Box "my_z_slice"]
# This file creates a bounding box which specified a volume or image that will
# be rendered. It is paired with a Render config file. If ProjectionAxis is set,
# an image will be rendered. If ProjectionAxis is not set, a 3D volume will
# be rendered.

# Location of lowermost corner:
X = 107.9
Y = 79
Z = 78.5

# Width of the box in each dimension:
XWidth = 42.14
YWidth = 42.14
ZWidth = 4.21

#######################
# Optional Parameters #
#######################

# Projection axis must be one of [ X | Y | Z ].
# ProjectionAxis = Z`
	ExampleBallFile = `[Ball "my_halo"]
# This file creates a bounding box defined by a sphere with a given radius.
# This is an alternative to using Box to specify rendering regions and is
# convenient when you're looking at individual haloes.

X = 4.602
Y = 100.7
Z = 80.7

Radius = 2.17

#######################
# Optional Parameters #
#######################

# Creates an image instead of a volume rendering when set. It must be one of
# [ X | Y | Z ].
# ProjectionAxis = Z

# Multiplies Radius by a constant.
# RadiusMultiplier = 3`
	ExampleTetraHistFile = `[TetraHist]

#######################
# Required Parameters #
#######################

# The property which a distribution is being measured for. Supported properties:
# Density: The mass-weighted density.
Quantity = Density

HistMin = 1e-2
HistMax = 1e4
HistBins = 200

# Must be "Log" or "Linear"
HistScale = Log

Particles = 50

Input = path/to/input/dir/
Output = path/to/outut/dir/
GridFile = path/to/gtet/grid.gtet

#######################
# Optional Parameters #
#######################

# SubsampleLength = 1

# Will result in files named pre_*foo*_app.txt:
# PrependName = pre_
# AppendName  = _app

# ProfileFile = prof.out
# LogFile = log.out`
)

type SharedConfig struct {
	// Required
	Input, Output string
	// Optional
	LogFile, ProfileFile string
}

func (con *SharedConfig) ValidInput() bool {
	return con.Input != ""
}
func (con *SharedConfig) ValidOutput() bool {
	return con.Output != ""
}
func (con *SharedConfig) ValidLogFile() bool {
	return con.LogFile != ""
}
func (con *SharedConfig) ValidProfileFile() bool {
	return con.ProfileFile != ""
}

type ConvertSnapshotConfig struct {
	SharedConfig
	// Required
	Cells int
	Gadget2IDSize int
	InputFormat string

	// Optional
	IteratedInput, IteratedOutput string
	IterationStart, IterationEnd int
}

func DefaultConvertSnapshotWrapper() *ConvertSnapshotWrapper {
	con := ConvertSnapshotConfig{}
	con.IterationStart = 0
	con.IterationEnd = -1
	con.Gadget2IDSize = 64
	return &ConvertSnapshotWrapper{con}
}


func (con *ConvertSnapshotConfig) ValidCells() bool {
	return con.Cells > 0
}
func (con *ConvertSnapshotConfig) ValidGadget2IDSize() bool {
	return con.Gadget2IDSize == 64 || con.Gadget2IDSize == 32
}
func (con *ConvertSnapshotConfig) ValidInputFormat() bool {
	return con.InputFormat != ""
}
func (con *ConvertSnapshotConfig) ValidIteratedInput() bool {
	return con.IteratedInput != ""
}
func (con *ConvertSnapshotConfig) ValidIteratedOutput() bool {
	return con.IteratedOutput != ""
}
func (con *ConvertSnapshotConfig) ValidIterationStart() bool {
	return con.IterationStart >= 0
}
func (con *ConvertSnapshotConfig) ValidIterationEnd() bool {
	return con.IterationEnd >= 0
}

type RenderConfig struct {
	SharedConfig
	
	// Required
	Quantity string
	TotalPixels, Particles int

	// Optional
	AutoParticles bool
	ImagePixels, ProjectionDepth int
	SubsampleLength int
	AppendName, PrependName string
}

func DefaultRenderWrapper() *RenderWrapper {
	rc := RenderConfig{ }
	rc.SubsampleLength = 1
	return &RenderWrapper{rc}
}

func (b *RenderConfig) ValidQuantity() bool {
	var q density.Quantity
	for q = 0; q < density.EndQuantity; q++ {
		if strings.ToLower(q.String()) == strings.ToLower(b.Quantity) {
			return true
		}
	}
	return false
}

func (con *RenderConfig) ValidTotalPixels() bool {
	return con.TotalPixels > 0
}

func (con *RenderConfig) ValidParticles() bool {
	return con.Particles > 0
}
func (con *RenderConfig) ValidSubsampleLength() bool {
	return con.SubsampleLength > 0
}
func (con *RenderConfig) ValidImagePixels() bool {
	return con.ImagePixels > 0
}
func (con *RenderConfig) ValidProjectionDepth() bool {
	return con.ProjectionDepth > 0
}

type ConvertSnapshotWrapper struct {
	ConvertSnapshot ConvertSnapshotConfig
}

type RenderWrapper struct {
	Render RenderConfig
}

type BallConfig struct {
	// Required
	X, Y, Z, Radius float64

	// Optional
	RadiusMultiplier float64
	Name string
}

func (ball *BallConfig) CheckInit(name string, totalWidth float64) error {
	if ball.Radius == 0 {
		return fmt.Errorf(
			"Need to specify a positive radius for Ball '%s'.", name,
		)
	}

	if ball.X >= totalWidth || ball.X < 0 {
		return fmt.Errorf(
			"X center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.X,
		)
	} else if ball.Y >= totalWidth || ball.Y < 0 {
		return fmt.Errorf(
			"Y center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.Y,
		)
	} else if ball.Z >= totalWidth || ball.Z < 0 {
		return fmt.Errorf(
			"Z center of Ball '%s' must be in range [0, %g), but is %g",
			name, totalWidth, ball.Z,
		)
	}

	ball.Name = name
	if ball.RadiusMultiplier == 0 {
		ball.RadiusMultiplier = 1
	} else if ball.RadiusMultiplier < 0 {
		return fmt.Errorf(
			"Ball '%s' given a negative radius multiplier, %g.",
			name, ball.RadiusMultiplier,
		)
	}

	return nil
}

func (ball *BallConfig) Box(totalWidth float64) *BoxConfig {
	box := &BoxConfig{}
	rad := ball.Radius * ball.RadiusMultiplier

	box.XWidth, box.YWidth, box.ZWidth = 2 * rad, 2 * rad, 2 * rad

	if ball.X > rad {
		box.X = ball.X - rad
	} else {
		box.X = ball.X - rad + totalWidth
	}

	if ball.Y > rad {
		box.Y = ball.Y - rad
	} else {
		box.Y = ball.Y - rad + totalWidth
	}

	if ball.Z > rad {
		box.Z = ball.Z - rad
	} else {
		box.Z = ball.Z - rad + totalWidth
	}

	box.Name = ball.Name

	return box
}

type BoxConfig struct {
	// Required
	X, Y, Z float64
	XWidth, YWidth, ZWidth float64

	// Optional
	ProjectionAxis string

	// Optional, "undocumented"
	Name string
}

func (box *BoxConfig) CheckInit(name string, totalWidth float64) error {
	if box.XWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive XWidth for Box '%s'", name,
		)
	} else if box.YWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive YWidth for Box '%s'", name,
		)
	} else if box.ZWidth <= 0 {
		return fmt.Errorf(
			"Need to specify a positive ZWidth for Box '%s'", name,
		)
	}

	if box.X >= totalWidth || box.X < 0 {
		return fmt.Errorf(
			"X origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.X,
		)
	} else if box.Y >= totalWidth || box.Y < 0 {
		return fmt.Errorf(
			"Y origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.Y,
		)
	} else if box.Z >= totalWidth || box.Z < 0 {
		return fmt.Errorf(
			"Z origin of Box '%s' must be in range [0, %g), but is %g",
			name, totalWidth, box.Z,
		)
	}

	tmp := box.ProjectionAxis
	box.ProjectionAxis = strings.Trim(strings.ToUpper(box.ProjectionAxis), " ")
	if box.ProjectionAxis != "" &&
		box.ProjectionAxis != "X" &&
		box.ProjectionAxis != "Y" &&
		box.ProjectionAxis != "Z" {

		return fmt.Errorf(
			"ProjectionAxis of Box '%s' must be one of [X | Y | Z]. '%s' is " + 
				"not recognized.", name, tmp,
		)
	}

	box.Name = name

	return nil
}

func (box *BoxConfig) IsProjection() bool { return box.ProjectionAxis != "" }

type BoundsConfig struct {
	Ball map[string]*BallConfig
	Box  map[string]*BoxConfig
}

func ReadBoundsConfig(fname string, totalWidth float64) ([]BoxConfig, error) {
	bc := BoundsConfig{}

	if err := gcfg.ReadFileInto(&bc, fname); err != nil {
		return nil, err
	}

	boxes := []BoxConfig{}
	for name, ball := range bc.Ball {
		if err := ball.CheckInit(name, totalWidth); err != nil {
			return nil, err
		}
		boxes = append(boxes, *ball.Box(totalWidth))
	}
	for name, box := range bc.Box {
		if err := box.CheckInit(name, totalWidth); err != nil {
			return nil, err
		}
		boxes = append(boxes, *box)
	}

	return boxes, nil
}

type TetraHistConfig struct {
	SharedConfig

	Quantity string

	HistMin, HistMax float64
	HistBins int
	HistScale string

	Particles int
	
	SubsampleLength int

	GridFile string
	AppendName, PrependName string
}

type TetraHistWrapper struct {
	TetraHist TetraHistConfig
}

func DefaultTetraHistWrapper() *TetraHistWrapper {
	cfg := TetraHistConfig{ SubsampleLength: 1 }
	return &TetraHistWrapper{ cfg }
}

func (con *TetraHistConfig) ValidInput() bool {
	return con.Input != ""
}

func (con *TetraHistConfig) ValidOutput() bool {
	return con.Output != ""
}

func (con *TetraHistConfig) ValidGridFile() bool {
	return con.GridFile != ""
}

func (con *TetraHistConfig) ValidQuantity() bool {
	switch con.Quantity {
	case "Density":
		return true
	}
	return false
}

func (con *TetraHistConfig) ValidHistMin() bool {
	return (con.HistScale == "Linear" ||
		(con.HistScale == "Log" && con.HistMin > 0)) && 
		con.HistMin < con.HistMax
}

func (con *TetraHistConfig) ValidHistMax() bool {
	return (con.HistScale == "Linear" ||
		(con.HistScale == "Log" && con.HistMax > 0)) &&
		con.HistMin < con.HistMax
}

func (con *TetraHistConfig) ValidHistBins() bool {
	return con.HistBins > 0
}

func (con *TetraHistConfig) ValidHistScale() bool {
	return con.HistScale == "Linear" || con.HistScale == "Log"
}

func (con *TetraHistConfig) ValidSubsampleLength() bool {
	s := con.SubsampleLength
	if s <= 0 { return false }

	for {
		if s == 1 { return true }
		if s % 2 == 1 { return false }
		s /= 2
	}
}
