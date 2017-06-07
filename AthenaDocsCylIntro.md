---
title: Introduction and Algorithms with Cylindrical Coordinates
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Introduction to cylindrical

Athena can compute solutions in cylindrical coordinates, using algorithms implemented by A. Skinner and E. Ostriker.  See the
[cylindrical coordinates method paper](http://adsabs.harvard.edu/abs/2010ApJS..188..290S)for details of the method and tests.

Configure
---------

To enable the cylindrical integrators, configure with

	% configure --with-coord=cylindrical

Riemann Solvers
---------------

The integration algorithm in cylindrical coordinates requires the time-centered pressure at cell
interfaces for second-order accuracy.  This is too expensive to compute using the conventional CTU
method, so instead we return it directly from the internal calculations
of the Riemann solver.  In the basic Riemann fan, this could refer to the
pressure of the left or right state in the case of supersonic flow, or
to the intermediate “star”-state.  The way this is computed varies
according to the particular Riemann solver.  We have appropriately
altered the Roe, HLLE, HLLC, and HLLD solvers to return this pressure;
to use any other Riemann solver with the cylindrical integrators, one
must alter them accordingly.  These are configured in the usual way:

	% configure --with-flux=roe
	% configure --with-flux=hlle
	% configure --with-flux=hllc
	% configure --with-flux=hlld

Spatial Reconstruction
----------------------

Reconstruction and characteristic evolution in the R-direction are
quite a bit different than in Cartesian coordinates.  Currently, we have
implemented modules for piecewise linear (2nd order) and quadratic (3rd
order) reconstruction.  To use other reconstructions with the cylindrical
integrators, one must rewrite them for the R-direction accordingly.
(Note that we make use of the Cartesian implementation for the $\phi-$
and z-directions.)  These are configured in the usual way:

	% configure --with-order=2
	% configure --with-order=3

With the CTU integrators, there is an additional half-timestep evolution
of the reconstructed states.  This is also significantly different in
the R-direction.  


Limitatations
-------------

 * Only a uniformally-spaced grid in all three coordinates $(r, \phi, z)$ is allowed.
 * SMR is not currently implemented with cylindrical coordinates
