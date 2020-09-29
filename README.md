# gotetra

`gotetra` is a Go-based package which uses phase-space tesselation techniques
to extract information about cosmological N-body simulations. The key
application this code is the generation of high resolution images of density fields.

This README describes how to install and use gotetra. Programmers interested in
modifying or using the source code can find documentation at https://godoc.org/github.com/phil-mansfield/gotetra.

![An example image](https://www.cfa.harvard.edu/~bdiemer/Visualizations/Images/Density_halo_L0063_s100_169074096_z_w15.0_t5.0_afmhot.png)

Further example images can be found at:
http://www.benediktdiemer.com/visualization/images/

## Installation

1. Install the Go compiler. The installation instructions can be found
here: https://golang.org/doc/install. If you are installing on a unix-based system
that you don't have root access to, follow the Linux instructions and replace
`/usr/local` with some directory that you do have access to (I use `~/local/go/`).
2. Test that the Go compiler works by running `$ go version`.
3. Create a directory that will host all your Go code. (I use `~/code/go/`).
4. Add a line to your `.profile`/`.bashrc`/`.bash_profile` which tells the Go compiler where
this code will be. If you your code directory is `~/code/go/`, write
`export GOPATH=$HOME/code/go`.
5. Download this package and dependencies into right directories by running
```
$ go get github.com/phil-mansfield/gotetra
$ go get gopkg.in/gcfg.v1
```
6. Go to `$GOPATH/src/github.com/phil-mansfield/gotetra/render/main` and run the command
`go build main.go`. This will create the main gotetra binary.
7. Test this binary was created correctly by running `./main -help`.

## Running Gotetra

Gotetra create images in two steps. First, it converts particle snapshots into a new format.
This is done using a config file which describes the location of your particle snapshots.
Second, it renders either a 2D image or a 3D density grid using a pair of config files which
describe rendering parameters and the portion of the simulation being rendered, respectively.
This will output a `.gtet` file, which can be read using either Gotetra or a Python pacakge.

### Step 1: Converting Particle Snapshots

To convert snapshots, you must run Gotetra on a machine with enough memory to store the entire
simulation in RAM and much have enough disk space to store a second copy of your snapshots.

Generate an example configuration file, `convert.cfg`, by running
`$ ./main -ExampleConfig ConvertSnapshot > convert.cfg`. Go through that example configuration
file and change the variables to match the simulation you are working with.

Start converting files by running `$ ./main -ConvertSnapshot convert.cfg`. This will use all the
threads on your machine unless you add an additional flag telling it how many to use
(e.g. `-Threads 8`).

### Step 2: Rendering Images

Rendering an image requires a rendering config file, `render.cfg` which specifies rendering
resolution and a bounds file, `box.cfg`, which specified the image/volume being rendered.

As in step 1, create example configuration files by running
```
$ ./main -ExampleConfig Render > convert.cfg
$ ./main -ExampleConfig Box > box.cfg
```
and change the variables to match your rendering targets. If you are rendering the area around a
halo, it may be more convenient to specify the rendering target using `./main -ExampleConfig Ball > ball.cfg`.

Render the image by running the command `$ ./main -Render render.cfg box.cfg`. This will create
a `.gtet` file at the directory given in `render.cfg`. By default, this command will use all the
threads on your machine, but you can change that with the `-Threads` flag (e.g. `-Threads 8`). You can
render multiple images at once by chaining them together at the end of the command
(e.g. `$ ./main -Render render.cfg box1.cfg ball1.cfg box2.gfc`).

## Python

Python code for interfacing with `gotetra` output is provided in the `python/`
directory. `gotetra.py` is both a python library and a command line utility.
Running it as `$ python gotetra.py my_gotetra_file.gtet` will print out
information about the file, and importing it will give you access to functions
which can read in `gotetra` headers and arrays. `gotetra.py` describes the functions
and data structures in more depth, but `read_header()` returns bit a bunch of 
information about the rendering and `read_grid()` returns the grid corresponding to the
image or volume being rendered.

`example.py` contains some example Python code that uses `gotetra.py`.

`render.py` is an incomprehensible blob of Python code that I use to generate images from
`.gtet` files. I don't plan to document or maintain this, but you are free to use it
if you'd like.

## Version

0.2

This project may experience breaking changes.

## License

MIT
