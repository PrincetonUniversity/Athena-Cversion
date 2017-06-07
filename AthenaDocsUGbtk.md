---
title: vtk (Legacy) Files
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/vtk Files

Much like binary dumps, vtk files contain an unformatted write of selected variables over all zones, with the data written
in [VTK](http://www.vtk.org) legacy format.  Each file contains data
at a single time.
A new file is created whenever the integration time exceeds an integer multiple of <output>/dt.  At the end of execution, the lesser of
`tlim/dt` or `<time>/nlim`
sequentially numbered files will be created.
Since they are unformatted, the outputs are compact, and are most useful for 3D data.  The
[VisIt](http://visitusers.org) visualization tool can read Athena vtk files directly.

vtk **dumps** contain all the primitive or conserved variables.  vtk **outputs** contain only a single
variable.  Each is described in more detail below.

Both vtk dumps and outputs are hardwired to be **single precision only**.  If data is needed in double precision, then the `dump_vtk.c` function would
need to be modified appropriately, and used as a new user-defined output function.

vtk Dumps
=========

The following example shows an <output> block in an input file that generates a vtk dump of the primitive variables:

	<output1>
	out_fmt = vtk               # vtk data dump
	out     = prim              # variables to be dumped
	dt      = 0.1               # time increment between outputs

vtk files consist of an ASCII header containing information about the dimensions and variables contained in the file,
followed by the unformatted data itself.  Since the header is ASCII format, vtk files can be edited to read the
header information.

vtk Outputs
===========

Single variables can also be output as vtk files.  In addition, new (user-defined) variables can be output as 
vtk files by following the steps described in
[User-defined Outputs]({{site.baseurl}}/AthenaDocsUGUserExpress).  For example,
the following shows an <output> block in an input file that generates a vtk dump of a user defined variable:

	<output1>
	out_fmt = vtk               # vtk data dump
	out     = dVy               # user defined variable: y-velocity fluctuations
	id      = dVy               # file name string identifier
	dt      = 0.1               # time increment between outputs
	usr_expr_flag = 1           # user defined variable

vtk output files are created by the function output_vtk.c.
