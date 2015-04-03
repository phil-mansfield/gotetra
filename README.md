gotetra
=======

`gotetra` is a Go-based package for calculating field densities in cosmological
N-body simulations. It exists as both a command-line utility and a library for
other Go code.

Extremely unfinished function documentation can be found at 
http://godoc.org/github.com/phil-mansfield/gotetra

Running Gotetra
---------------

`gotetra`'s command line utility can be buily by running the command

    $ go build

in the directory `main/`. Some Go compilers [read: the Go compiler on U Chicago's
Midway cluster] support the following build command:

    $ go build -gcflags=-B

which will speed up execution through the removal of bounds checks. Either of these
commands will create a `main` binary in that directory. This binary can be run in several
different modes using the following pattern:

    $ ./main -ModeName mode_config_file.txt [bounds_file_1.txt bounds_file_2.txt ...]

The supported modes are currently `-Density` and `-ConvertSnapshot`.

I plan to make detailed documentation for different types of config files eventually,
but for now, you can create example config files with minimal parameter descriptions
my running `$ ./main -Config Density`, `$ ./main -Config ConvertSnapshot`, or
`$ ./main -Config Bounds`. Hopefully, I've made the parameters simple enough
that the reader can figure it out from there.

Python
------

Python code for interfacing with `gotetra` output is provided in the `python/`
directory. `gotetra.py` is both a python library and a command line utility. Running
it as `$ python gotetra.py my_gotetra_file.gtet` will print out information about
the file, and importing it will give you access to functions which can read in
`gotetra` headers and arrays.

The documentation for this is not fantastic, and I plan to make it better later.

`render.py` is a terrible hacky script that I wrote which is completely undocumented
and renders `gotetra` files into pngs. I've included it in case the reader would find
it helpful.

Bugs
----

I garuantee there are still bugs lurking around. They should mostly be gone from the
actually rendering code, but I haven't tested the actual interface code nearly as
thoroughly. In particular, the free parameters that I chose for the auto-particle count
choosers are completely arbitrary.

Let me know if you find anything that doesn't work right.

Performance and Memory Footprint
--------------------------------

The exact runtime of a `gotetra` job will depend upon a number of factors (which I
will go into more detail on when I get more time). Although the program is relatively
well-behaved in its scaling behavior, it is always a good idea to gradually scale up
to the largest pixel and particle counts to ensure that memory and time constaints won't
be.

To give the reader an order of magnitude estimate, creating a 16 megapixel image of a
63 Mpc/h simulation box with 1500 particles pertetrahedron on a 16-core
Midway node consumes less than 600 MB, and finishes execution in 20 minutes.

Version
-------

0.0.1

This project is not even close to stable yet.

License
----

MIT
