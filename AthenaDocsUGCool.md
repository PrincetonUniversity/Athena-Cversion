---
title: Optically-thin Cooling
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Cooling

Any arbitrary cooling and/or heating function can be added as a source term in the energy equation.  The implementation uses
function pointers which are set in the problem generator, providing the user flexibility in the form of the cooling function.
To set the cooling function pointer (called `CoolingFunc`), add the line

	CoolingFunc = myfunc;

anywhere in the problem generator.  *myfunc* must be a function of type `Real` with the density, pressure, and time step
as arguments, for example:

	Real myfunc(const Real dens, const Real Press, const Real dt)
	{
	  ...
	}

The function should return the *net cooling (heating) rate per unit volume*.  Any function of the density and pressure (temperature) is
allowed.  The time step is included as an argument so that the total change in energy can be limited to prevent negative pressures.

Some simple examples are given in the file `/athena/src/microphysics/cool.c`.

**Optically-thin cooling only works with the CTU integrator.**  The cooling source terms are added to the reconstruction and interface-state
correction steps in the integrator, so that the *cooling terms are fully second-order*.
