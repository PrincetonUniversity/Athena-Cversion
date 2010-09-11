#include "../copyright.h"
/*==============================================================================
 * FILE: radHD_hlle.c
 *
 * PURPOSE: Computes 1D fluxes using the Harten-Lax-van Leer (HLLE) Riemann
 *   solver.  This flux is very diffusive, especially for contacts, and so it
 *   is not recommended for use in applications.  However, as shown by Einfeldt
 *   et al.(1991), it is positively conservative, so it is a useful option when
 *   other approximate solvers fail and/or when extra dissipation is needed.
 *
 * REFERENCES:
 *   E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999) chpt. 10.
 *
 *   Einfeldt et al., "On Godunov-type methods near low densities",
 *   JCP, 92, 273 (1991)
 *
 *   A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and
 *   Godunov-type schemes for hyperbolic conservation laws",
 *   SIAM Review 25, 35-61 (1983).
 *
 *   B. Einfeldt, "On Godunov-type methods for gas dynamics",
 *   SIAM J. Numerical Analysis 25, 294-318 (1988).
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifndef SPECIAL_RELATIVITY


#ifdef RADHD_HLLE_FLUX
/*----------------------------------------------------------------------------*/
/* 
 *   Input Arguments:
 *     Bxi = B in direction of 1D slice at cell interface
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void rad_fluxes(const Cons1DS Ul, const Cons1DS Ur,
                   const Prim1DS Wl, const Prim1DS Wr,
                   const Real Bxi, Cons1DS *pFlux, const Real dt)
{
	Real aeffl, aeffr, bp, bm, al, ar, tmp;
	int n;
	/* Flux of the left and right state for the conserved quantity */
	Cons1DS Fl,Fr;
	Real *pFl, *pFr, *pF;

	aeffl = eff_sound(Ul,dt);
	aeffr = eff_sound(Ur,dt);


/*--- Step 1. ------------------------------------------------------------------
 * We may use  Roe-averaged data from left- and right-states. Not for now
 */

	al = Wl.Vx - aeffl;
	ar = Wr.Vx + aeffr;

	bp = MAX(ar, 0.0);
 	bm = MIN(al, 0.0);


/*--- Step 5. ------------------------------------------------------------------
 * Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
 */

 	Fl.d  = Ul.Mx - bm*Ul.d;
 	Fr.d  = Ur.Mx - bp*Ur.d;

 	Fl.Mx = Ul.Mx*(Wl.Vx - bm);
 	Fr.Mx = Ur.Mx*(Wr.Vx - bp);

 	Fl.My = Ul.My*(Wl.Vx - bm);
  	Fr.My = Ur.My*(Wr.Vx - bp);

  	Fl.Mz = Ul.Mz*(Wl.Vx - bm);
  	Fr.Mz = Ur.Mz*(Wr.Vx - bp);


	Fl.Mx += Wl.P;
	Fr.Mx += Wr.P;

  	Fl.E  = Ul.E*(Wl.Vx - bm) + Wl.P*Wl.Vx;
	Fr.E  = Ur.E*(Wr.Vx - bp) + Wr.P*Wr.Vx;


/*--- Step 6. ------------------------------------------------------------------
 * Compute the HLLE flux at interface.
 */

  	pFl = (Real *)&(Fl);
  	pFr = (Real *)&(Fr);
  	pF  = (Real *)pFlux;
  	tmp = 0.5*(bp + bm)/(bp - bm);
 	 for (n=0; n<(NWAVE+NSCALARS); n++){
	    pF[n] = 0.5*(pFl[n] + pFr[n]) + (pFl[n] - pFr[n])*tmp;
	 }
	
  

  return;
}
#endif /* RADIHD_HLLE_FLUX */
#endif
