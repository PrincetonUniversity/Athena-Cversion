---
title: Shearing Box and Orbital Advection
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Shearing Box

Shearing Box Approximation
==========================

To study the dynamics of accretion flows in the local shearing-box approximation, Coriolis and tidal
gravity source terms must be added to the equations of motion.  See the 
[Shearing Box Method Paper](http://adsabs.harvard.edu/abs/2010ApJS..189..142S) for details.

To use the shearing-box approximation, configure with

	% configure --enable-shearing-box

Currently only the CTU integrator can be used with the shearing-box.

Orbital Advection
=================

Orbital advection (sometimes called the FARGO algorithm) integrates the equations of motion with the
background orbital motion removed.  This can greatly increase the efficiency of the algorithm by increasing
the time step.  Orbital advection can be used both with the shearing-box approximation for local studies,
and with a cylindrical grid for global simulations.

A new algorithm for orbital advection with MHD, based on the Constrained Transport method, has been
implemented in Athena.  See the [Shearing Box Method Paper](http://adsabs.harvard.edu/abs/2010ApJS..189..142S) for details.

Configure using

	% configure --enable-fargo

The boundary conditions must be periodic in the direction of the shearing motion.
