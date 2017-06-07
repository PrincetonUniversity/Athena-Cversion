---
title: The <Output> Block
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Output Block(s)

Parameters in these blocks control the writing of *outputs*,
e.g. restart dumps, images, history files, etc.  There can be any number of output blocks in an input file, and they
can be listed in any order.  Example


        <output1>
        out_fmt = hst                # History data dump
        dt      = 62.831853          # time increment between outputs
        
        <output2>
        out_fmt = rst                # Restart dump
        dt      = 6.2831853e3
        
        <output3>
        out_fmt = bin                # Binary data dump
        dt      = 628.31853          # time increment between outputs
        
        <output4> 
        out_fmt = ppm                # ppm image dump
        out     = dVy                # output variable
        id      = dVy                # id added to file names
        usr_expr_flag = 1            # indicates dVy is user defined quantity
        palette = jh_colors          # color palette for images
        dt      = 62.831853          # time step between output of delta V3
        dmin    = -0.0006            # min value for imaging delta V3
        dmax    =  0.0006            # max value for imaging delta V3
        x2      = 0.0                # slice in X-Z plane at x2=0.0


**out:** Variable(s) to be output.  Currently accepted values are: `cons, prim, d, M1, M2, M3, E, B1c, B2c, B3c, ME, V1, V2, V3, P, S, cs2, G`.  If `out=cons` (or
`out=prim`), then `out_fmt` must be one of `bin, tab` or `vtk`, and output variables are the conserved (or primitive) variables.

**out_fmt:** Output format, e.g.
`bin, hst, tab, rst, vtk, pdf, pgm, ppm`.  See the 
[User Guide]({{site.baseurl}}/AthenaDocsUG) for more
about these formats.

**dat_fmt:** Optional field for controlling the format string used
to write tabular output files, e.g. `%12.5e`.  This value should
not appear in quotes and no white space should be present.

**dt:** Time increment between outputs (in problem time).

**time:** Time of next output (in problem time).  If not set,
the default will be the initial problem time (for new runs), or the
current problem time (for restarts).

**id:** Any string, added to label output filenames.

**dmin/dmax:** max/min applied to output (useful for images).

**palette:** Color palette for images.  Currently available palettes
are `rainbow, jh_colors, idl1, idl2, step8, step32, heat`.

**x1, x2, x3:** Value of x1, x2, or x3 at which data is to be output.  Either a single value (which
will produce a slice of the data at that value),
or a range of values over
which data is to be averaged, can be specified.  In either case, the data is
reduced in dimension by one for each value specified.  If the value specified lies outside
the computational domain, no output will occur, and no errors will be generated.  See [Specifying Slices]({{site.baseurl}}/AthenaDocsUGSlices)
for more details.

**usr_expr_flag:** Set to 1 if a user-defined expression is to
be used to compute output quantity, see [User Expressions]({{site.baseurl}}/AthenaDocsUGUserExpress) for more details.

**level:** Level of mesh to output; root level is specified by `level=0`.  Set to -1 to produce output
on *all* levels.  Default is -1.

**domain:** Integer index of Domain at each level of the Mesh to output.  Counting starts at zero, so the root Domain is `level=0, domain=0`.
Useful if output on only one Domain in a SMR
hierarchy is desired.  Set to -1 to produce output on *all* Domains at each level.  Default is -1.
