---
title: Specifying Output Slices
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Specifying Slices

Often it is useful to output a *slice* of data from the computational domain, rather than the entire array.  This
can be achieved by specifying values for the **x1, x2, x3** parameters in the appropriate input block (see
[Output Blocks]({{site.baseurl}}/AthenaDocsUGOutputBlck) in the User Guide).

Either a single value (which will output a slice),
or a range of values over
which data is to be averaged, can be specified.  In either case, the data is
reduced in dimension by one for each value specified.

In the following example, the density at x2=0.0 will be output in .tab format.  The output will be a 1D vector in a 2D calculation,
and a 2D slice in the x1-x3 plane in a 3D calculation.

        <output2>
        out_fmt = tab
        out     = d
        x2      = 0.0


In the following example, the density will be output in .vtk format averaged over the entire x1 axis.  This is particularly
useful, e.g., for making images of the column density.

        <output2>
        out_fmt = vtk
        out     = d
        x1      = :


In the following example, the density will be output in .vtk format averaged over the range 0.0 < x1 < 5.0.

        <output2>
        out_fmt = vtk
        out     = d
        x1      = 0.0:5.0


In the following example, the density will be output in .vtk format averaged over the range x1min < x1 < 5.0.

        <output2>
        out_fmt = vtk
        out     = d
        x1      = :5.0


In the following example, the density will be output from a 3D grid in .tab format as a 1D vector at x1=2.0, x2=4.0.

        <output2>
        out_fmt = tab
        out     = d
        x1      = 2.0
        x2      = 4.0


If the value specified lies outside
the computational domain, no output will occur, and no errors will be generated.  This is important to remember with
MPI and/or SMR.  In particular, if a Grid on a given processor (with MPI), or the Domain at some level (with SMR) does not contain data within
the specified range or slice, no output will be generated.
