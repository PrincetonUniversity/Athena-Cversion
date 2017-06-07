---
title: Grid Based Particle Output
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Grid Based Output

Grid-binned particle output is essentially the same as that for gas quantities. One example is given as below


	<output1>
	  out_fmt	= ppm		/* ppm image */
	  out		= M1par	 	/* particle x1-momentum */
	  id		= M1par
	  x2		= 1.0
	  dt		= 1.0		/* time increment between outputs */


Here the variable `M1par` denotes the particle x1-momentum. All built-in quantities for the particle output include **dpar, M1par, M2par, M3par, V1par, V2par, V3par**, which denote particle density, momentum and velocity respectively. The item `pargrid` is automatically set to 1 for these built-in output variables. By default, [particle property selection function]({{site.baseurl}}/AthenaDocsParOutOverview) is set to `all`, but the user is free to define their own property selection functions.

The built-in variables can be retrieved in the grid from `pG->Coup[k][j][i]`, where `pG` is the pointer to the grid, and `Coup` is of type `GPCouple` (see `/athena/src/athena.h`). It is mainly used for storing the intermediate step gas quantities for the particle integrators, but is also used in storing the grid-binned particle quantities for output purposes. Here its elements include `grid_d` (particle density), `grid_v1`, `grid_v2`, `grid_v3` (particle momentum). The user can define their own particle-related output variables based on these grid-binned quantities in the same way as for gas variables (see [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress)), but in such cases `pargrid` in the [output block]({{site.baseurl}}/AthenaDocsParOutOverview) must be set to 1 explicitly.

For [binary]({{site.baseurl}}/AthenaDocsUGbin)/[vtk]({{site.baseurl}}/AthenaDocsUGbtk)/[tab]({{site.baseurl}}/AthenaDocsUGtab) data dumps, they contain binned particle variables `dpar`, `M1par`, `M2par` and `M3par` with property selection function `all`, which can **NOT** be changed.
