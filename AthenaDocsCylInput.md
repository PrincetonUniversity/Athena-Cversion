---
title: Input Files and Problem Generators with Cylindrical Coordinates
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Input and Generators with cylindrical

To use
the cylindrical integrators, one must remember that the coordinate order
is always $(x1,x2,x3) = (R,\phi,z)$.  Thus, in the level 0 domain, x1min
and x1max specify the minimum and maximum values of R and so forth.  The
current implementation of cylindrical coordinates expects the coordinate
singularity at the origin to be excised from the computational domain
and unexpected behavior may result from a lack of proper regularization
if one uses a domain which includes it.  Also, note that in Athena,
the $\phi$ coordinate, hence also x2min and x2max, is measured in radians.
The standard boundary conditions behave as one would expect in the $\phi-$ and z-directions.  For example, a domain with

	x2min  = 0
	x2max  = 3.141592653589793
	bc_ix2 = 4
	bc_ox2 = 4

describes an annular region periodic in $\phi$.

The standard boundary conditions do not apply in the R-direction.
Typically one must write and enroll their own function in a problem
generator according to some particular application.  We have written a
simple Dirichlet boundary condition function called `do_nothing_bc`, which
leaves the initialized ghost zones intact throughout the integration.

Note that the nature of constrained transport in Athena causes the
divergence of the magnetic field to be preserved to within machine
precision of its initial value.  It is the user's responsibility to
ensure that this initial value is 0, recalling that the divergence
operator includes extra geometric scale factors in non-Cartesian
coordinate systems.  We have included a routine `compute_div_b()` located in
`utils.c` to help the user with debugging.

One typically refers to interface-centered and cell-centered quantities in
problem files.  However, it must be noted that the conserved quantities in
Athena are the area-averages of the magnetic fields over cell interfaces,
and the volume-averages of the density, momenta, and sometimes total
energy.  The easiest way to obtain the correct area-averages of the
magnetic fields, which also ensures that the divergence vanish, is to use
a vector potential.  Additionally, one must specify the volume-averages of
the magnetic field.  The easiest way to do this is to set the interface
fields first, then use the following second-order accurate averages to
set the cell-centered fields:

	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	Rc = x1;  // R at cell center
	Rl = x1 – 0.5*pG->dx1;  // R at left interface
	Rr = x1 + 0.5*pG->dx1;  // R at right interface
	
	pG->U[k][j][i].B1c = (Rl*pG->B1i[k][j][i] + Rr*pG->B1i[k][j][i+1])/Rc;
	pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
	pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);

It can be tricky to initialize the remaining variables with the correct
volume-averages, but in some applications this is especially important.
Consider, for example, a centrifugally balanced system with solid-body
rotation profile (i.e. no shearing).  Long term stability of this system
relies on careful balance between the inward gravitational and outward
centrifugal forces, both of which are derived from volume-averaged
source terms.  If these forces are unbalanced, the system will move
radially inward or outward and begin to shear.

To aid the user in producing highly accurate initial volume-averages,
we have included the numerical quadrature functions  `avg1d()`, `avg2d()`,
and `avg3d()` located in `utils.c`.  These functions take a function of three spatial
coordinates and produce the corresponding volume-averages in one, two,
and three spatial dimensions, respectively.  These functions do not
scale well with increasing dimension, but the user may specify an error
tolerance and maximum iteration number, and they need only be run to
initialize a problem.  Typically, one can get away with averaging in
the R-direction only.  For example, to study the Rayleigh stability
of a system with shear parameter q, we could compute $M_\phi$ by hand
as a function of q, or we could simply define a function `M2()` for the
azimuthal momentum,

	Real M2(const Real x1, const Real x2, const Real x3) {
	  return rho0*omega0*pow(x1,1.0-q);
	}

Then to initialize the momentum in the $\phi$-direction, we set

	pG->U[k][j][i].M2 = avg1d(M2,pG,i,j,k);

Note that this averaging is not required!  In most cases, it is sufficient
(and second-order accurate) to equate a volume-average with the pointwise
value of a function at the grid cell center.
