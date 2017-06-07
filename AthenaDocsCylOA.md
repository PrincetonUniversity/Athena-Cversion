---
title: Orbital Advection in Cylindrical Coordinates
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Orbital Advection in cylindrical


The implementation of cylindrical coordinates in Athena includes an orbital
advection, or FARGO, algorithm which greatly increases the time step in problems consisting
of a supersonic background flow, for example geometrically thin accretion disks.
Use of orbital advection can yield a
speed-up of up to the maximum Keplerian Mach number of the simulation;
in practice generally an order of magnitude for thin disks.  Orbital advection
is turned on by using the configure flag `--enable-fargo`.  Currently,
orbital advection is limited to the
CTU integrator.  The algorithm supports
both second- and third-order reconstruction.

Properly setting up orbital advection in a cylindrical simulation
is more difficult than the implementation for shearing
box (Cartesian) simulations.  First of all, two user-specified functions of the
cylindrical radius are required, one representing the Keplerian angular
velocity and the other representing the shearing parameter.  The shear
parameter, \\(q\\), is given as \\(q = -0.5 * dln(\Omega^2)/dlnR\\), where \\(\Omega =
\Omega(R)\\) is the Keplerian angular velocity profile.  For example, given
an unstratified Newtonian potential \\(Phi = -1/R\\), \\(\Omega(R) = R^{-1.5}\\)
and \\(q(R) = 1.5\\).  These functions would be defined similarly to a static
gravitational potential, for instance:

	static Real Omega(const Real R) {
	  Real Arg;
	  Arg = pow(R,1.5);
	  return (1.0/Arg);
	}
	static Real Shear(const Real R) {
	  return 1.5;
	}

These functions need to be enrolled in the main body of the problem file as:

	OrbitalProfile = Omega;
	ShearProfile = Shear;


Once the function pointers `OrbitalProfile` and `ShearProfile` have been
defined, the only remaining change from the standard cylindrical geometry
is in the definition of the x2-momentum.  When using orbital advection
the x2-momentum is not meant to represent the full azimuthal velocity,
but only the velocity in the locally-rotating frame.  For a precisely
Keplerian disk this means that the x2-momentum is identically zero.
An example of proper initialization is given below.

	#ifdef FARGO
	       pG->U[k][j][i].M2 = 0.0;
	#else  
	       pG->U[k][j][i].M2 = pG->U[k][j][i].d*avg1d(vphi,pG,i,j,k);
	#endif

A useful example of orbital advection is the problem file `cylnewtmri.c`
in the standard Athena distribution.  Worth noting, is that because the
x2-momentum is stored as the value in the rotating frame it will also be
output that way.  This can sometimes require that the Keplerian velocity
to be added back in using postprocessing scripts.
