---
title: Static and Self-Gravity
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Gravity

Acceleration due to gravitational forces can be added in several ways.

Static gravitational potential
=============================

For gravitational potentials that are fixed in time, the appropriate source terms in the momentum and energy
equations are added using function pointers.  To set the gravitational potential function pointer (called `StaticGravPot`),
add the line

	StaticGravPot = myfunc;

anywhere in the problem generator.  *myfunc* must be of type `Real`, and take the `x1`, `x2`, and `x3` coordinate
positions as arguments, for example

	static Real myfunc(const Real x1, const Real x2, const Real x3)
	{
	  ...
	}

The function should return the gravitational potential.

Static potentials work with both the CTU and VL integrators.  In both cases, the momentum and energy source term updates are
*fully second-order*, and moreover, the total momentum and energy (including the gravitational potential energy) are 
*conserved exactly*.

For an example of how to include a static potential in a problem generator, see `/src/prob/rt.c`.

**With cylindrical coordinates, a function to compute the gravitational acceleration given the coordinates is required
in addition to the potential function above.**

Self-gravity
============

Adding self-gravity is far more complicated than adding a fixed potential; e.g. Poisson's equation must be solved
for the potential.  In Athena, there are two algorithms implemented for the solution of Poisson's equation, one based on
FFTs, the other using multigrid.

To use the FFT-based algorithm, configure the code with

	% configure --with-gravity=fft --enable-fft

The code uses the [FFTW](http://www.fftw.org/) library, which must be pre-installed in the system and linked during the compile step.
This requires setting default paths to the FFTW libraries, or editing the `Makeoptions.in` file to 
set paths using the `FFTWLIB` and `FFTWINC` macros.  (This latter step is similar to setting paths
for the MPI libraries for parallel calculations.)  **This method only works with periodic boundary conditions.**

To use the multigrid-based algorithm, configure the code with

	% configure --with-gravity=multigrid

No external libraries are needed, and any type of boundary condition is allowed.

In both cases, there are two parameter values that must be set in the problem generator; generally these are read from the `<problem>`
block in the input file using the lines

	four_pi_G     = par_getd("problem","four_pi_G");
	grav_mean_rho = par_getd("problem","grav_mean_rho");

added anywhere in the input file.  The first parameter effectively sets the gravitational constant G, and therefore sets the units used in the
simulation.  The second parameter is the mean density in the computational volume, which is required when FFTs are used
with periodic boundary conditions.  This parameter can be calculated from the initial conditions, rather than input, if desired.
If values for these parameters are not set, the code should exit with an error message.

In both cases, the gravitational accelerations are *second-order*, and the algorithm applies the accelerations as
the divergence of the gravitational stress tensor, so that the *total momentum is conserved exactly*.  Total energy
(including the gravitational potential energy) is not conserved, however.

**Self-gravity does not work with SMR.  Multigrid does not work with MPI.**
