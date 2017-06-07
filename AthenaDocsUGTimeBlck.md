---
title: The <time> Block
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Time Block

Parameters in this block control times in a job (such as ending time).

        <time>
        cour_no         = 0.4       # The Courant, Friedrichs, & Lewy (CFL) Number
        nlim            = 100000    # cycle limit
        tlim            = 6.2832e4  # time limit (10 orbits)

**tlim:** Time to stop integration, in units defined by problem.

**nlim:** Maximum number of cycles of the main loop before stopping.
Set to -1 to stop only on time limit `tlim`.

**cour_no:**  CFL number, must be less than 1.0 for 1D and 2D CTU integrator,
and 0.5 for 2D VL integrator, or 3D.
