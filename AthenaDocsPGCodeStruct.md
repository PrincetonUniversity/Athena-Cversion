---
title: Code Structure
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Code Structure

The main integration loop is contained in `main.c`, which calls functions to
initialize the Mesh, Domains, and Grids, and orchestrates I/O.  Which
integrator is actually used is determined by a function pointer
set at run time by `integrate_init()` in the file `/athena/src/integrators/integrate.c`.

The image below shows a schematic flow chart for Athena in a multidimensional calculation
![alt Flow Chart]({{site.baseurl}}/images/flow_chart.png)

The 1D, 2D, and 3D integrators are very different in structure.
However, much effort has been put into keeping them as modular as
possible.  Thus, they all share the same reconstruction functions (e.g.,
`lr_states_prim2.c`, etc.) and flux functions (e.g., the Roe Riemann solver,
etc.).  Eigensystems in the conserved and primitive variables, needed by
some of the Riemann solvers and reconstruction algorithms respectively,
are contained in the `esystem_roe.c` and `esystem_prim.c` files.

The input file is parsed by the functions in `par.c`, and I/O
is coordinated by the functions in `output.c`.

Users who wish to make substantial modification or extension of
the algorithms should start by understanding the 1D integrator
(`integrate_1d.c`) and the functions it calls, before moving on to the 2D
and 3D integrators.
