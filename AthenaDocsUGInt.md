---
title: Integrators
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Integrators

The third important algorithmic element of a Godunov scheme is the Integrator: the way in which time-averaged fluxes
are calculated at each interface in multidimensions.  **Directional splitting cannot be used in MHD**, at least
if the equations are solved in the conservative form.  Thus, some kind of unsplit integrator is required.  In Athena,
two different unsplit integrators are implemented.

Each of the integrators are implemented in separate files in the directory `/athena/src/integrators/`

Corner Transport Upwind (CTU) Integrator (default)
==================================================

Specify during configure using

	% configure --with-integrator=ctu

The most complicated, and therefore the most accurate, algorithm.  Not all physics options
will work with CTU.  However, if possible, **use of CTU over VL is always recommended**.

See the 
[3D JCP Method Paper](http://adsabs.harvard.edu/abs/2008JCoPh.227.4123G) for complete details of our extension
of the CTU integrator to MHD using Constrained Transport.

MUSCL-Hancock (VL) Integrator
=============================

Specify during configure using

	% configure --with-integrator=vl

A much simpler, and more robust, algorithm.  Recommended for problems in extreme regimes, e.g. highly
supersonic an/or strongly magnetized turbulence,
where stability of the algorithm is a limiting factor.

Some physics (e.g. special relativity) only works with the VL integrator.  An error message should be generated
during the configure step if physics that does not work with the chosen integrator is selected.

See the [van Leer Integrator Method Paper](http://adsabs.harvard.edu/abs/2009NewA...14..139S) for complete details,
and comparisons of the accuracy of the VL and CTU integrators on various test problems.

### First-order flux correction

The stability of the VL integrator can be improved by enabling this option with configure, using

	% configure --with-integrator=vl --enable-fofc

Whenever the unsplit update produces a cell with negative density or pressure, the calculation *for this cell only*
is repeated using first-order (donor cell) fluxes.  The extra diffusion this introduces is often enough to 
prevent unphysical states, without compromising the accuracy of the solution anywhere else in the domain.
