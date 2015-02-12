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

    $ go get go get github.com/phil-mansfield/gotetra

or by manually installing via git in `$GOPATH/src/github.com/phil-mansfield/gotetra`.

To create a usable binary, go to `$GOPATH/src/github.com/phil-mansfield/gotetra/main` and
type

    $ go build

This will create a binary called `main`. You can see what command line argumetns `main` takes
by running

    $ `./main --help`

Version
-------

0.0.1

This project is not even close to stable yet.

License
----

MIT