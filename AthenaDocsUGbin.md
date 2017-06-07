---
title: Binary Files
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Binary Files

Contain unformatted write of all (primitive or conserved) dependent variables over all active zones.  A new sequentially numbered
file is created whenever the integration time exceeds an integer multiple of `<output>/dt`. At the end of execution, the lesser of tlim/dt or `<time>/nlim`
sequentially numbered files will be created.

An example of an `<output>` block which creates binary dumps is given below.

	<output2>
	out_fmt = bin                # Binary data dump
	dt      = 628.31853          # time increment between outputs

Since the `out=prim` parameter is not specified in the above example, by default the conserved variables are output.

Binary files created by Athena contain several words of data in a header that describe the dimensions and physics of the problem, followed
by the data itself.  The data must be read by special scripts that read this header information, see 
[Visualizing Outputs]({{site.baseurl}}/AthenaDocsUGVis).
Created by the function dump_binary.c; consult this function for more information about the header data.

Binary dumps are compact, and can be used effectively for storing 3D data.  However, the data is not portable between machines with
different endianness.

Binary dumps are hardwired to be **single precision only**.  If data is needed in double precision, then the `dump_binary.c` function would
need to be modified appropriately, and used as a new user-defined output function.
