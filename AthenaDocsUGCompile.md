---
title: The Compile Step
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Compiling

After running `configure`, the next step is to compile the code
using `make`.  In the same directory as the configure script, use:

	% make all

This automatically creates
the directory `./athena/bin`, which will contain the executable,
and runs make in `./athena/src` and its subdirectories to compile and link
the code.

The top-level Makefile contains the following targets

| TARGET | Comments |
|--------|----------|
| all | creates `./athena/bin` directory and compiles code |
| compile | compiles code |
| clean | removes .o files in `./athena/src` |
| help | prints help message |
| test | runs install test (see the [Tutorials/QuickStart]({{site.baseurl}}/AthenaDocsTutQuickStart) page) |

The top-level Makefile uses an additional macro `MACHINE` that can be used
to select compiler options and library paths specific to a certain computer.  For example,


	make all MACHINE=ophir


will compile with `icc -O3 -xW -ipo -i static` as well as use 
paths for the FFT and MPI libraries specific to this machine.  Values for this macro are
set in `Makeoptions.in`, and new values (machines, libraries, and paths) can be added
to this file.

Remember that the makefile is re-generated from the template `Makefile.in` every time configure is
run.  Thus, any changes made to the makefile will be lost if configure is run again.  **Any permanent modifications 
should be made to `Makefile.in`, and not the `Makefile`**.
