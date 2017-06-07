---
title: Overview of Data Output
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Output Overview

Data output in the Athena code is
controlled by the `<output>` blocks in the input file (see [Output Block(s)]({{site.baseurl}}/AthenaDocsUGOutputBlck) for more
details).  There should
be one block for each type of data output required.  There is no limit on
the total number of outputs.  

Output filenames follow the convention

	basename-id#.dumpid.outid.type

*basename* is inherited
from the `problem_id` parameter in the `<job>` block
(see [Job Block]({{site.baseurl}}/AthenaDocsUGJobBlck))

*-id#* labels the
processor id for jobs run with MPI with *#* an integer equal to the rank
of the MPI process (the root process does not contain an
*-id0*, nor is it present for serial jobs)

*dumpid* is a zero
filled unsigned integer of length `num_digit` (currently
fixed at four, this value can be changed in `./athena4.0/src/defs.h.in`)

*outid* is the string specified by the `id` parameter in the `<output>` block.
(if not specified, the string will be `out-#`, where `#`
denotes the block number in the input
file which generated the output)

*type* denotes the output format (`bin, tab, hst, vtk, rst, pdf, pgm, ppm`).

Note that history filenames do not include a *dumpid* or
*outid*.  Also note that if the output is a *dump* (a `.bin`, `.vtk`, or `.tab` output of all the
conserved or primitive variables), the filename does not
include the *outid* string. 

It is important to note that **output files in Athena will always be
silently overwritten**, and **history files will always be appended**. 

More information about each of the output file formats is provided by the subsections in the User Guide.
