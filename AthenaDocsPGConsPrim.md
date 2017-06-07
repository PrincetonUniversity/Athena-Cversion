---
title: Structures for the Conserved and Primitive Variables
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Variable Structures

Conserved Variables (ConsS structure)
====================================

The cell-centered, volume-averaged values of each of the
conserved variables are stored in the *ConsS* structure, which is defined in `/athena/src/athena.h` and
reproduced below.

	typedef struct Cons_s{
	  Real d;                       /* density */
	  Real M1;                      /* momentum density in 1,2,3 directions */
	  Real M2;
	  Real M3;
	#ifndef BAROTROPIC
	  Real E;                       /* total energy density */
	#endif /* BAROTROPIC */
	#ifdef MHD
	  Real B1c;                     /* cell centered magnetic fields in 1,2,3 */
	  Real B2c;
	  Real B3c;
	#endif /* MHD */
	#if (NSCALARS > 0)
	  Real s[NSCALARS];             /* passively advected scalars */
	#endif
	}ConsS;


The vector **M** denotes the linear momentum.  The ConsS structure also stores the face-centered field, but this is
not the primary description of the field.  For adiabatic MHD, the structure contains 8 variables, plus passively advected scalars.

Primitive Variables (PrimS structure)
====================================

The corresponding structure containing the primitive variables is *PrimS*,
also defined in `/athena/src/athena.h`.

	typedef struct Prim_s{
	  Real d;                       /* density */
	  Real V1;                      /* velocity in 1,2,3 */
	  Real V2;
	  Real V3;
	#ifndef BAROTROPIC
	  Real P;                       /* pressure */
	#endif /* BAROTROPIC */
	#ifdef MHD
	  Real B1c;                     /* cell centered magnetic fields in 1,2,3 */
	  Real B2c;
	  Real B3c;
	#endif /* MHD */
	}PrimS;


Conserved and Primitive Variables in 1D: the Cons1D and Prim1D structures
=========================================================================

The longitudinal component of the magnetic field is dropped in the one-dimensional
system of equations.  Thus, we define new structures for the dependent variable
in one-dimension, each of which only contains 7 variables (for adiabatic MHD), plus passively advected scalars.
For the conserved variables, the structure is Cons1DS:

	typedef struct Cons1D_s{
	  Real d;                       /* density */
	  Real Mx;                      /* momentum density in X,Y,Z; where X is     */
	  Real My;                      /* direction longitudinal to 1D slice; which */
	  Real Mz;                      /* can be in any dimension: 1,2,or 3         */
	#ifndef BAROTROPIC
	  Real E;                       /* total energy density */
	#endif /* BAROTROPIC */
	#ifdef MHD
	  Real By;                      /* cell centered magnetic fields in Y,Z */
	  Real Bz;
	#endif /* MHD */
	#if (NSCALARS > 0)
	  Real s[NSCALARS];             /* passively advected scalars */
	#endif
	}Cons1DS;

Note that we have used the subscripts x,y, and z to denote the components
of vectors in the one-dimensional variables contained in Cons1DS and
Prim1DS, whereas we used the subscripts 1,2, and 3 to denote the components
of vectors in the ConsS and PrimS structures.  In the latter, the variables
are fixed with respect to the coordinates of the grid (e.g., `M1`
corresponds to the momentum in the 1-direction).  However, since the
one-dimensional variables may represent a slice in any direction, the
x,y, and z components are not fixed with respect to the grid (that is,
`Mx` corresponds to `M1` along a slice in the 1-direction,
but `Mx` corresponds to `M2` along a slice in the 2-direction).

The order of the variables in these structures is
extremely important and
cannot be changed for several reasons.  Firstly, this order determines the
order of the elements used in the eigensystem of the linearized equations
(computed in the functions contained in the files `esystem_prim.c`
and `esystem_roe.c`).  Secondly, in several functions we set a
pointer to the first element in the structure, and then address successive
elements (variables) by incrementing the pointer, as in the
following code fragment from the function `lr_states_prim2.c`

	    pWl = (Real *) &(Wl[i+1]);
	
	    qx = 0.5*MAX(ev[NWAVE-1],0.0)*dtodx;
	    for (n=0; n<NWAVE; n++) {
	      pWl[n] = Wrv[n] - qx*dW[n];
	    }

Here, `Wl`, `Wrv`, and `dW` are all structures of type Prim1DS,
and `NWAVE` is the number of components of these structures (which
depends on whether the problem is hydrodynamic or MHD, and adiabatic or
isothermal).  Using pointers in this way leads to more compact coding,
is more efficient (since it allows vectorization of the loop), and by
using structures to ensure the components of the vectors are stored
contiguously, it is more cache efficient.  However, *be warned*.
It also means the order of the variables is hardwired in the code,
and cannot be changed, which has the potential to be the source of some
nasty bugs.
