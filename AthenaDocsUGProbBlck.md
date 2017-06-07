---
title: The <Problem> Block
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Problem Block


Parameters in this block are accessed by the problem-generator,
and therefore depend on the problem being run.  For example, the
problem-generator `./athena4.0/src/prob/shkset1d.c`
(which is used for the Brio & Wu test problem)
requires the following parameter values in the `<problem>` block:


	<problem>
	gamma           = 2.0       # gamma = C_p/C_v
	shk_dir         = 1         # Shock Direction -- (1,2,3) = (x1,x2,x3)
	
	dl              = 1.0       # density on left half of grid
	pl              = 1.0       # pressure
	v1l             = 0.0       # X-velocity
	v2l             = 0.0       # Y-velocity
	v3l             = 0.0       # Z-velocity
	b1l             = 0.75      # X-magnetic field
	b2l             = 1.0       # Y-magnetic field
	b3l             = 0.0       # Z-magnetic field
	s0l             = 1.0       # "color"
	
	dr              = 0.125     # density on right half of grid
	pr              = 0.1       # pressure
	v1r             = 0.0       # X-velocity
	v2r             = 0.0       # Y-velocity
	v3r             = 0.0       # Z-velocity
	b1r             = 0.75      # X-magnetic field
	b2r             = -1.0      # Y-magnetic field
	b3r             = 0.0       # Z-magnetic field
	s0r             = 0.0       # "color"

The only parameter values in the `<problem>` block that are common to all input files are:

**gamma:** ratio of specific heats, only used for adiabatic EOS.

**iso_csound:** isothermal sound speed, only used for isothermal EOS.

All of the remaining parameter names in the example above are specific to the Brio &
Wu shocktube problem.  In general, the input file for other problems
will have different variable names in the problem block.
