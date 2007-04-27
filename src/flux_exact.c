#include "copyright.h"
/*==============================================================================
 * FILE: flux_exact.c
 *
 * PURPOSE: Computes 1D fluxes using exact nonlinear Riemann solver.
 *   Currently only isothermal hydrodynamics has been implemented.  
 *
 * REFERENCES:
 *   R.J. LeVeque, "Numerical Methods for Conservation Laws", 2nd ed.,
 *   Birkhauser Verlag, Basel, (1992).
 *
 *   E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics",
 *   2nd ed., Springer-Verlag, Berlin, (1999).
 *
 * HISTORY:
 *   dec-2006  Isothermal hydro version written by Nicole Lemaster
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   flux_exact()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MHD
#error : The exact flux only works for hydro.
#endif /* MHD */

#ifndef ISOTHERMAL
#error : The exact flux only works for isothermal equation of state.
#endif /* ISOTHERMAL */

static void srder(double dm, double vl, double vr, double dmin, double dmax, 
			double *y, double *dydx);
static double rtsafe(void (*funcd)(double, double, double, double, double,
			double *, double *), double x1, double x2, double xacc, 
			double vl, double vr, double dmin, double dmax);

/*----------------------------------------------------------------------------*/
/* flux_exact:
 *   Input Arguments:
 *     Bxi = B in direction of slice at cell interface
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void flux_exact(const Real Bxi, const Cons1D Ul, const Cons1D Ur, Cons1D *pF)
{
  Prim1D Wl, Wr;
  Real zl, zr, zm, dm, Vxm, Mxm, tmp, pbl, pbr, dmin, dmax;
  Real sl, sr;    /* Left and right going shock velocity */
  Real hdl, hdr;  /* Left and right going rarefaction head velocity */
  Real tll, tlr;  /* Left and right going rarefaction tail velocity */
  char soln;      /* two bits: 0=shock, 1=raref */

  if(!(Ul.d > 0.0)||!(Ur.d > 0.0))
    ath_error("[flux_exact]: Non-positive densities: dl = %e  dr = %e\n",
	      Ul.d, Ur.d);

/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);

/*--- Step 2. ------------------------------------------------------------------
 * Compute the density and momentum of the intermediate state
 */

  zl = sqrt((double)Wl.d);
  zr = sqrt((double)Wr.d);

  /* --- 1-shock and 2-shock --- */
  soln = 0;

  /* Start by finding density if shocks on both left and right.
   * This will only be the case if dm > Wl.d and dm > Wr.d */
  tmp = zl*zr*(Wl.Vx - Wr.Vx)/(2.0*Iso_csound*(zl + zr));
  zm = tmp + sqrt((double)(tmp*tmp + zl*zr));
  dm = zm*zm;

  /* Get velocity from 1-shock formula */
  Vxm = Wl.Vx - Iso_csound*(dm-Wl.d)/(zm*zl);

  /* If left or right density is greater than intermediate density,
   * then at least one side has rarefaction instead of shock */
  dmin = MIN(Wl.d, Wr.d);
  dmax = MAX(Wl.d, Wr.d);
  if (dm < dmax) {
    /* --- 1-rarefaction and 2-rarefaction --- */
    soln = 3;

    /* Try rarefactions on both left and right, since it's a quicker
     * calculation than 1-shock+2-raref or 1-raref+2-shock */
    dm = zl*zr*exp((Wl.Vx-Wr.Vx)/(2.0*Iso_csound));

    /* Get velocity from 1-rarefaction formula */
    Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);

    /* If left or right density is smaller than intermediate density,
     * then we must instead have a combination of shock and rarefaction */
    if (dm > dmin) {
      /* --- EITHER 1-rarefaction and 2-shock ---
       * --- OR     1-shock and 2-rarefaction --- */

      /* Solve iteratively equation for shock and rarefaction
       * If Wl.d > Wr.d ==> 1-rarefaction and 2-shock
       * If Wr.d > Wl.d ==> 1-shock and 2-rarefaction */
      if (Wl.d > Wr.d) soln = 2; else soln = 1;

      dm = rtsafe(&srder,dmin,dmax, 2.0*DBL_EPSILON, Wl.Vx, Wr.Vx, dmin, dmax);

      /* Don't be foolish enough to take ln of zero */
      if ((dm > dmin) && (dm <= dmax)) {
        if (Wl.d > Wr.d) {
          /* Get velocity from 1-rarefaction formula */
          Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);
        } else {
          /* Get velocity from 2-rarefaction formula */
          Vxm = Wr.Vx + Iso_csound*log(dm/Wr.d);
        }
      } else {
        /* --- DEFAULT 1-rarefaction and 2-rarefaction --- */
        soln = 3;

        /* In the event that the intermediate density fails to fall between
         * the left and right densities (should only happen when left and
         * right densities differ only slightly and intermediate density
         * calculated in any step has significant truncation and/or roundoff
         * errors), default to rarefactions on both left and right */
        dm = zl*zr*exp((Wl.Vx-Wr.Vx)/(2.0*Iso_csound));

        /* Get velocity from 1-rarefaction formula */
        Vxm = Wl.Vx - Iso_csound*log(dm/Wl.d);
      }
    }
  }

  if (dm < 0.0)
    ath_error("[flux_exact]: Solver finds negative density %5.4e\n", dm);

/*--- Step 3. ------------------------------------------------------------------
 * Calculate the Interface Flux if the wave speeds are such that we aren't
 * actually in the intermediate state
 */

  if (soln & 2) { /* left rarefaction */
    /* The L-going rarefaction head/tail velocity */
    hdl = Wl.Vx - Iso_csound;
    tll = Vxm - Iso_csound;

    if (hdl >= 0.0) {
      /* To left of rarefaction */
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.d*Iso_csound2;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);
      return;
    } else if (tll >= 0.0) {
      /* Inside rarefaction fan */
      dm = Ul.d*exp(hdl/Iso_csound);
      Mxm = Ul.d*Iso_csound*exp(hdl/Iso_csound);
      Vxm = (dm == 0.0 ? 0.0 : Mxm / dm);

      pF->d  = Mxm;
      pF->Mx = Mxm*Vxm + dm*Iso_csound2;
      pF->My = Mxm*Wl.Vy;
      pF->Mz = Mxm*Wl.Vz;
      return;
    }
  } else { /* left shock */
    /* The L-going shock velocity */
    sl = Wl.Vx - Iso_csound*sqrt(dm)/zl;

    if(sl >= 0.0) {
      /* To left of shock */
      pF->d  = Ul.Mx;
      pF->Mx = Ul.Mx*(Wl.Vx) + Wl.d*Iso_csound2;
      pF->My = Ul.My*(Wl.Vx);
      pF->Mz = Ul.Mz*(Wl.Vx);
      return;
    }
  }

  if (soln & 1) { /* right rarefaction */
    /* The R-going rarefaction head/tail velocity */
    hdr = Wr.Vx + Iso_csound;
    tlr = Vxm + Iso_csound;

    if (hdr <= 0.0) {
      /* To right of rarefaction */
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.d*Iso_csound2;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return;
    } else if (tlr <= 0.0) {
      /* Inside rarefaction fan */
      tmp = dm;
      dm = tmp*exp(-tlr/Iso_csound);
      Mxm = -tmp*Iso_csound*exp(-tlr/Iso_csound);
      Vxm = (dm == 0.0 ? 0.0 : Mxm / dm);

      pF->d  = Mxm;
      pF->Mx = Mxm*Vxm + dm*Iso_csound2;
      pF->My = Mxm*Wr.Vy;
      pF->Mz = Mxm*Wr.Vz;
      return;
    }
  } else { /* right shock */
    /* The R-going shock velocity */
    sr = Wr.Vx + Iso_csound*sqrt(dm)/zr;

    if(sr <= 0.0) {
      /* To right of shock */
      pF->d  = Ur.Mx;
      pF->Mx = Ur.Mx*(Wr.Vx) + Wr.d*Iso_csound2;
      pF->My = Ur.My*(Wr.Vx);
      pF->Mz = Ur.Mz*(Wr.Vx);
      return;
    }
  }

/* If we make it this far, then we're in the intermediate state */

/*--- Step 4. ------------------------------------------------------------------
 * Calculate the Interface Flux */

  if(Vxm >= 0.0){
    pF->d  = dm*Vxm;
    pF->Mx = dm*Vxm*Vxm + dm*Iso_csound2;
    pF->My = dm*Vxm*Wl.Vy;
    pF->Mz = dm*Vxm*Wl.Vz;
  }
  else{
    pF->d  = dm*Vxm;
    pF->Mx = dm*Vxm*Vxm + dm*Iso_csound2;
    pF->My = dm*Vxm*Wr.Vy;
    pF->Mz = dm*Vxm*Wr.Vz;
  }

  return;
}

/* Equation to solve iteratively for shock and rarefaction as well
 * as its derivative.  Used by rtsafe() */

static void srder(double dm, double vl, double vr, double dmin, double dmax, 
			double *y, double *dydx)
{
  *y = (vr - vl) + Iso_csound*(log(dm/dmax) + (dm-dmin)/sqrt(dm*dmin));
  *dydx = Iso_csound/dm*(1.0 + 0.5*(dm+dmin)/sqrt(dm*dmin));

  return;
}

/* Numerical Recipes function rtsafe modified to use doubles and take
 * extra parameters to pass to the derivative function above */

#define MAXIT 100

static double rtsafe(void (*funcd)(double, double, double, double, double,
			double *, double *), double x1, double x2, double xacc, 
			double vl, double vr, double dmin, double dmax)
{
	int j;
	double df,dx,dxold,f,fh,fl;
	double temp,xh,xl,rts;

	(*funcd)(x1,vl,vr,dmin,dmax,&fl,&df);
	(*funcd)(x2,vl,vr,dmin,dmax,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
		return 0.0;
	}
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,vl,vr,dmin,dmax,&f,&df);
	for (j=1;j<=MAXIT;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
				if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
		(*funcd)(rts,vl,vr,dmin,dmax,&f,&df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
        /* Cap the number of iterations but don't fail */
	return rts;
}

#undef MAXIT
