---
title: Introduction to Athena
layout: page
---

Introduction
============

Athena is a grid-based code for astrophysical magnetohydrodynamics (MHD).
It was developed primarily for studies
of the interstellar medium, star formation, and accretion flows.
Athena has been made freely available to the community
in the hope that others may find it useful.

The current version (v4.2) implements algorithms for the following physics:
* compressible hydrodynamics and MHD in 1D, 2D, and 3D,
* special relativistic hydrodynamics and MHD,
* ideal gas equation of state with arbitrary \\(γ\\) (including \\(γ = 1\\), an isothermal EOS),
* an arbitrary number of passive scalars advected with the flow,
* self-gravity, and/or a static gravitational potential,
* Ohmic resistivity, ambipolar diffusion, and the Hall effect,
* both Navier-Stokes and anisotropic (Braginskii) viscosity,
* both isotropic and anisotropic thermal conduction,
* optically-thin radiative cooling.

In addition, Athena allows for the following grid and parallelization options:
* Cartesian or cylindrical coordinates,
* static (fixed) mesh refinement,
* shearing-box source terms, and an orbital advection algorithm for MHD,
* parallelization using domain decomposition and [MPI](http://www.mcs.anl.gov/mpi). 
A variety of choices are also available for the numerical algorithms, such as different Riemann solvers
and spatial reconstruction methods. 

The code has been developed using the [GNU](http://www.gnu.org) development tools,
maintaining strict adherence to ANSI standards, thus it should be possible
to configure, compile, and run the code on any platform that supports these
standards.  Athena has been run on everything from a Mac laptop to a 25,000 processor Cray XT4.

Learn More
==========
[Equations Solved:](/AthenaEqns) What system of equations does the code actually solve?

[Documentation:](/AthenaDocs) Tutorial, User Guide, Programmer Guide, and research
papers describing the numerical algorithms in Athena.

[Gallery:](http://www.astro.princeton.edu/~jstone/Athena/athena-apps/index.html) Images and links to research papers that include results from Athena.

[Tests:](http://www.astro.princeton.edu/~jstone/Athena/tests/index.html) Suite of 1D, 2D, and 3D hydrodynamic and MHD problems used to validate
Athena.  Potentially useful to anyone developing algorithms for MHD.

[Parallel Performance:](/AthenaDocsScaling) Results from strong and weak scaling tests of Athena, up to 10^5^ cores.

[People:](/AthenaDocsPeople) The list of collaborators developing Athena.

[Is Athena Galilean invariant?](/AthenaDocsGalileo)  Well, for resolved solutions, yes.

Download
========

[Downloads:](/AthenaDocsDownLd) Get the complete source code distribution for the latest public version of Athena.

What's New?
===========

* **Jan 2011:** v4.1 has been released, with algorithm extensions for special relativistic hydrodynamics and MHD, as well as various bug fixes.
* **Jul 2013:** v4.2 has been released, updates are mostly bug fixes from v4.1

Other Links
===========

[Athena3D in Fortran](http://www.astro.virginia.edu/VITA/athena.php) written by John Hawley and Jake Simon.

[Minerva](http://www.astro.umd.edu/~askinner/minerva), the original cylindrical coordinate version of Athena written by Aaron Skinner and
Eve Ostriker.
