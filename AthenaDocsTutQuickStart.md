---
title: Quick Start
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[Tutorial]({{site.baseurl}}/AthenaDocsTut)/QuickStart


To install, configure, compile, and run a test of the code do the following.

 1. Download, then uncompress and untar the source code distribution file.

	% gunzip athena4.0.tar.gz
	% tar xf athena4.0.tar


 2. Create the configure script by running `autoconf` in the `./athena` directory.

	% cd athena4.0
	% autoconf


 3. Test the install by running `configure`; `make all`; `make test`.

	% configure
	% make all
	% make test
	(cd tst/1D-mhd; ./run.test)
	zone-cycles/cpu-second = 4.826353e+05
	zone-cycles/wall-second = 2.407434e+05
	L1 norm for density: 6.333383e-11

 `make test` runs a linear wave convergence test on a grid of 512 zones, then computes and outputs the L1 error norm in the fast magnetosonic wave compared
 to the analytic solution. If the benchmark fails to run, or if the resulting error norm is large, then something has gone wrong in the installation
 or in the compilation of the code. 

 The zone-cycles/cpu-second is a useful measure of the code performance, it corresponds to the number of grid cells updated per cpu second. This of course depends on the physics included, the
 geometry of the problem, as well as the processor used for the calculation. The values reported above are for a 2.6Ghz Intel Nehalum processor. The error norm is absolute, 
 and should be very small.

 4. The steps above will have created an executable called `athena` in a new directory `./athena/bin`.  There are various command line arguments you can use
 with the executable, one of the most useful is

	% cd bin
	% athena -c

 This will output information about the configuration (physics and algorithm options) of this executable.  It is extremely useful if you running many problems at once
 and forget which executable does what.  A list of the valid options is given by the `-h` option.

If there are no errors in any of the above steps
the code can now be re-configured, re-compiled, and used to run any
of the test problems in `./athena/src/prob`.  For example, the next step for new users would be to run a 1D hydrodynamics problem and
plot the result (see the next section in the Tutorial).
