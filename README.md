gotetra
=======

`gotetra` is a Go-based package for calculating field densities in cosmological
N-body simulations. It is based on the method outlined in Abel, Hahn, and
Kaehler's numerous recent papers on the subject and is based on projecting the
simulation's phase sheet onto position space. This project is serving as a
prototype for an eventual C-based package, `tetra` which will be a more heavily
optimized and documented version of the same project.

Extremely unfinished function documentation can be found at 
http://godoc.org/github.com/phil-mansfield/gotetra

Project Structure
-----------------

 This project contains the following subpackages:

- Package `io` handles all catalog and phase sheet io.

- Package `geom` provides functions for computing properties of 3D objects
within domains with periodic boundary conditions.

- Package `rand` contains functions for computing random floating point
numbers.

- Package `density` interpolates particle sequences onto density grids.

- Package `main` compiles to a very hacky binary which can be run with a host of command line options.

Installing This Package on Midway
---------------------------------

First you must install Go. This can be done by putting the following in your `.bashrc`

    module load go/1.3
    export GOPATH=$HOME/code/go

(the exact location of `$GOPATH` can be whatever you want). All go packages will live
within a subdirectory of `$GOPATH/src/`. Since packages rely on knowing their exact location
with respect to `$GOPATH` you should either install this package by calling

    $ go get github.com/phil-mansfield/gotetra

or by manually installing via git in `$GOPATH/src/github.com/phil-mansfield/gotetra`.

To create a usable binary, go to `$GOPATH/src/github.com/phil-mansfield/gotetra/main` and
type

    $ go build

This will create a binary called `main`. You can see what command line argumetns `main` takes
by running

    $ `./main --help`

You will almost certainly want to run some permutation of the command line arguments found in
`main/gotetra.sh`.

Command Line Arguments to `main`
--------------------------------

calls to `main` will look like this:

    $ ./main *Mode Flag* *Other Flags* path/to/input/directory path/to/output/directory

Since I made several extremely poor design decisions back when I first made `main.go` a few
months ago, `main.go` can run in multiple modes. The most useful of which is `BoundedDensity`
which allows you to obtain density fields for small zoomed in boxes. This is what `gotetra.sh`
is set to do by default. You will probably want to change some of the other flags, though.

- `Cells` controls how many cells would be needed to span the entire width of the box. If you are
zooming in on a halo, this will generally mean that the full grid would be orders of magnitude
larger than what could be held in physical memory. (Increasing this number increases the memory footprint
of the program. Every running thread has its own buffers, so you can trade some performance for some
memory by reducing the number of threads.)

- `Points` the number of interpolation points per tetrahedron. (Increasing this number has the obvious
effects, but also increases the memory footprint of the program, since interpolation points are
heavily cached for performance reasons. If this becomes a problem you can reduce `UnitBufferLen` on line
62 of `gotetra.go`. If you reduce this number too far, your density fields will start to obtain stripes,
so exercise caution.)

- `Skip` the number of particles to skip when subsampling.

- `Bounds file` the file which contains the bounding boxes of regions that will be sampled. (This is a
text file which should consist of rows which each contain six space-separated numbers. The first three
will be the coordinates of the lower corner of a bounding box and the next three will correspond to the
width of the box in each direction.)

Version
-------

0.0.1

This project is not even close to stable yet.

License
----

MIT
