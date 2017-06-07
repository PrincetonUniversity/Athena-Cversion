---
title: Visualizing Outputs from Athena
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Visualization

Athena does not come with a default graphics package.  Instead,
the user must decide which visualization package is best suited to their
needs, output the data in a format which can be read by this package,
and then proceed.  As a start, rudimentary scripts for several different
graphics packages are supplied with the source code; future versions
may incorporate more sophisticated visualization tools.  The following
subsections describe useful visualization packages for Athena data files
(the discussion assumes the code has already been run to produce output).

Supermongo
==========

A popular package for making publication-quality one-dimensional plots
is [SM](http://www.astro.princeton.edu/~rhl/sm). 
A simple SM macro that can read tabular output from Athena is provided in 
`./athena/vis/sm`.

IDL procedures
==============

[IDL (Interactive Data Language)](http://www.rsinc.com)
procedures that can read both the binary and VTK dump files and make plots
are included in `./athena/vis/idl`.  To run these procedures
IDL must be installed on the system.  From the `./athena/bin` directory, use
the following to read a binary file and make some one-dimensional plots:

        % idl
        IDL> .run ../vis/idl/pltath.pro 
        IDL> nine_plot,'Brio-Wu.0040.bin',1

A variety of potentially useful
procedures are included in the `pltath.pro` file.

MatLab .m files
===============

[Matlab](http://www.mathworks.com/) .m files can be used to read Athena .tab, .bin, and .vtk files.
A number of examples are provided in the `./athena/vis/matlab` directory.

Using VisIt to read VTK files
=============================

We have found the [VisIt](http://www.llnl.gov/visit) package, now supported by the DOE
[VACET](http://www.vacet.org) center,
useful for plotting 3D data sets.  The VTK legacy
format produced by Athena can be read by VisIt.  For data created with
executables parallelized with MPI, a C code that joins multiple files into
one is provided in `./athena/vis/vtk`.  This is useful for jobs run
on massively parallel clusters.

2D animations
=============

Athena can output 2D images that can be displayed directly, or easily
turned into animations.  For example, if a series of ppm images of a single
variable have been created, they can be displayed either using
[ImageMagick](http://www.imagemagick.org):

        % animate *.ppm

or, alternatively, `xanim` (which requires converting the ppm images to FLI
format):

        % ls Wind*ppm > list1
        % ppm2fli -g80x80 list1 Wind.fli
        % xanim Wind.fli
