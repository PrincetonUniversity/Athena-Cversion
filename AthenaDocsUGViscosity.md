---
title: Viscosity
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Viscosity

Both isotropic (Navier-Stokes) and anisotropic (Braginskii) viscosity can be added.  In either case,
configure the code with

        % configure --enable-viscosity

Both are added at first-order via operator splitting.

The update is explicit in time, so that a very restrictive CFL constraint on the time step will be used.

Navier-Stokes viscosity
=======================

Enable by setting a value for the coefficient of kinematic viscosity `nu_iso` in the problem generator.  For example, to read a value
from the `<problem>` block in the input file, or to set a default value of zero if a value is not specified in the input file, add the line

        nu_iso = par_getd_def("problem","nu_iso",0.0);

anywhere in the problem generator.
If the code detects `nu_iso > 0`, then the viscous stresses will be added.  Currently `nu_iso` must be a constant.

Braginskii viscosity
====================

In this case, viscous stresses are applied only in the direction parallel to the magnetic field lines.  **MHD must be enabled**
to use Braginskii viscosity.

Enable by setting a value for the coefficient of parallel kinematic viscosity `nu_aniso` in the problem generator.  For example, to read a value
from the `<problem>` block in the input file, or to set a default value of zero if a value is not specified in the input file, add the line

        nu_aniso = par_getd_def("problem","nu_aniso",0.0);

anywhere in the problem generator.
If the code detects `nu_aniso > 0`, then the viscous stresses will be added.  Currently `nu_aniso` must be a constant.

To add *both* non-zero parallel and isotropic viscosities, simply set values for both `nu_iso` and `nu_aniso` in the
problem generator.

If both `nu_iso` and `nu_aniso` are zero, and viscosity is enabled, the code will terminate with an error message.
