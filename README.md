gotetra
=======

`gotetra` is a Go-based package for calculating field densities in cosmological
N-body simulations. It is based on the method outlined in Abel, Hahn, and
Kaehler's numerous recent papers on the subject and is based on projecting the
simulation's phase sheet onto position space. This project is serving as a
prototype for an eventual C-based package, `tetra` which will be a more heavily
optimized and documented version of the same project.

Unfinished documentation can be found at 
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

Version
-------

0.0.1

This project is not even close to stable yet.

License
----

MIT