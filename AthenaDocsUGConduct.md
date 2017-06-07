---
title: Thermal Conduction
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Thermal Conduction

Both isotropic and anisotropic thermal conduction can be added.  In either case,
configure the code with

        % configure --enable-conduction

Both are added at first-order via operator splitting.

The update is explicit in time, so that a very restrictive CFL constraint on the time step will be used.

Isotropic thermal conduction
============================

Enable by setting a value for the coefficient of thermal diffusion `kappa_iso` in the problem generator.  For example, to read a value
from the `<problem>` block in the input file, or to set a default value of zero if a value is not specified in the input file, add the line

        kappa_iso = par_getd_def("problem","kappa_iso",0.0);

anywhere in the problem generator.
If the code detects `kappa_iso > 0`, then thermal conduction will be applied.  Currently `kappa_iso` must be a constant.

If both `kappa_iso` and `kappa_aniso` (see below) are zero, the code assumes an error was made in initializing these constants, and aborts.

Anisotropic thermal conduction
==============================

In this case, the heat flux is parallel to the magnetic field lines.  **MHD must be enabled**
to use anisotropic thermal conduction.

Enable by setting a value for the coefficient of parallel thermal diffusion `kappa_aniso` in the problem generator.  For example, to read a value
from the `<problem>` block in the input file, or to set a default value of zero if a value is not specified in the input file, add the line

        kappa_aniso = par_getd_def("problem",kappa_aniso",0.0);

anywhere in the problem generator.
If the code detects `kappa_aniso > 0`, then thermal conduction will be added.  Currently `kappa_aniso` must be a constant.

To add *both* non-zero parallel and isotropic thermal conduction, simply set values for both `kappa_iso` and `kappa_aniso` in the
problem generator.

If both `kappa_iso` and `kappa_aniso` are zero, and conduction is enabled, the code will terminate with an error message.

Units
=====

Note the kappa's are **diffusivities**, not **conductivities**.  Also note the current implementation
uses "dimensionless units" in that the factor (mbar/k_B) is not
included in calculating the temperature (instead, T=P/d is used).  For cgs
units, kappa must be entered in units of [cm^2^/s], and the heat fluxes calculated in the counduction functions would
need to be multiplied by (k_B/mbar).
