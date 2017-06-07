---
title: Data Output with Cylindrical Coordinates
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Data Output with cylindrical

Most output formats work seamlessly with cylindrical coordinates.  Visualization
packages such as [VisIt](http://www.llnl.gov/visit) can read Athena `.vtk` files and plot them in cylindrical coordinates automatically.

We have found it useful to use certain data post-processing functions 
with binary dumps.  Thus we have added a coordinate system flag to the
`dump_binary.c` function.  Currently, this flag is either -1 for Cartesian
coordinates or -2 for cylindrical (we made these negative to avoid
confusion with the various Boolean flags).  As usual, the default output
formats are enrolled in the input files.  This allows packages such as 
[IDL](http://www.rsinc.com) to
recognize which coordinate system is being used.

A large number of
[Matlab](http://www.mathworks.com/) scripts for data
analysis are now included in the `/vis/matlab` directory, some
of which are useful in cylindrical coordinates.
