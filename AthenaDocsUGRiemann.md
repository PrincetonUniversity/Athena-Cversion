---
title: Riemann Solvers
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Riemann Solvers

The Riemann solver is the method by which time-averaged fluxes of all conserved quantities
are calculated at cell interfaces, see section 4.3 in the 
[ApJS Method Paper](http://adsabs.harvard.edu/abs/2008ApJS..178..137S).  There are entire
monographs written on exact and approximate Riemann solvers for hydrodynamics and MHD
(e.g. E.F. Toro, *Riemann Solvers and Numerical Methods for Fluid Dynamics*, 1999).
References to "Toro" below refer to this book.

Along with [Reconstruction]({{site.baseurl}}/AthenaDocsUGReconstruct) and the 
[Integrator]({{site.baseurl}}/AthenaDocsUGInt), the Riemann solver is one of the most important
algorithmic elements of a Godunov scheme.  For these reason, a variety of choices for the
solver are implemented in Athena.

The Riemann solvers are implemented in functions in the directory `/athena/src/rsolvers/`

To specify the Riemann solver in Athena, configure the code with

	% configure --with-flux=choice

where the table below gives the valid choices implemented in Athena.

| Choice | Comment | Physics | Reference |
|:------ | ------- | ------- | --------- |
| force  | Toro's FORCE flux | Hydro and MHD | Toro section 7.4.2 |
| two-shock | Two-shock approximation | Hydro | Toro section 9.4.2 |
| exact | eact solver | Hydro | Toro chapter 4 |
| hlle   | Harten-Lax-van Leer with Einfeldt fix | Hydro and MHD | Toro section 10.3 |
| hllc   | Harten-Lax-van Leer with contact | Hydro | Toro section 10.4 |
| hlld   | Harten-Lax-van Leer with contact and Alfven mode | MHD | Miyoshi & Kusano, JCP, 208, 305 |
| Roe    | Roe's linearized solver | Hydro and MHD | Toro chapter 11 |

For **hydrodynamics**, use of the **Roe or HLLC solver** is strongly recommended.

For **MHD**, use of the **Roe or HLLD solver** is strongly recommended.

The exact solvers are useful for testing, but are generally too slow for applications, and often do not increase
the accuracy of solutions in any case.
