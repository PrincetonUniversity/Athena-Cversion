---
title: Obtaining and Installing Athena
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Installing

The source code for Athena can be downloaded as a gzipped tar file from this site.
To install and run the code requires only a C compiler.
After downloading the tar file, uncompress and untar it:

	% gunzip athena4.0.tar.gz
	% tar xf athena4.0.tar

The following directory structure should be created:

	/athena
	       /apps           problem and input files used in previous applications
	       /doc            selected additional documentation
	       /src            source code, and include files
	           /fftsrc          FFT interface
	           /gravity         self-gravity 
	           /integrators     unsplit integrators
	           /microphysics    resistivity, thermal conduction, viscosity, etc.
	           /particles
	           /prob            problem files (/src/problem.c is a symbolic link to one file in this dir)
	           /reconstruction  spatial reconstruction algorithms
	           /rsolvers        Riemann solvers
	       /tst            various input files for tests
	           /1D-hydro
	           /1D-mhd
	           /1D-sr-hydro
	           /1D-sr-mhd
	           /2D-hydro
	           /2D-mhd
	           /2D-sr-hydro
	           /2D-sr-MHD
	           /3D-hydro
	           /3D-mhd
	           /3D-sr-hydro
	           /3D-sr-MHD
	           /cylindrical
	           /particle
	       /vis            visualization tools and scripts
	           /idl        IDL (rsinc.com) scripts
	           /matlab     MatLab .m files   
	           /particle   
	           /sm         Super Mongo scripts
	           /vtk        code to join VTK legacy files



In addition to those above, another directory will be created by the Makefile
when Athena is compiled for the first time, using `make all`.

	/athena
	       /bin           contains executable, created by Makefile

(Trying to compile the code with `make compile` before
`make all` will result in an error since the `bin` directory will
not yet exist.)
