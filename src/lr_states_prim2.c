#include "copyright.h"
/*==============================================================================
 * FILE: lr_states_prim2.c
 *
 * PURPOSE: Second order (piecewise linear) spatial reconstruction using
 *   characteristic interpolation in the primitive variables.  A time-evolution
 *   (characteristic tracing) step is used to interpolate interface values to
 *   the half time level {n+1/2}, unless the unsplit integrator in 3D is VL. 
 *
 * NOTATION: 
 *   U_{L,i-1/2} is reconstructed value on the left-side of interface at i-1/2
 *   U_{R,i-1/2} is reconstructed value on the right-side of interface at i-1/2
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 *   U_{L,i-1/2} is denoted by Ul[i  ];   U_{R,i-1/2} is denoted by Ur[i  ]
 *   U_{L,i+1/2} is denoted by Ul[i+1];   U_{R,i+1/2} is denoted by Ur[i+1]
 *
 *   Internally, in this routine, Wlv and Wrv are the reconstructed values on
 *   the left-and right-side of cell center.  Thus (see Step 9),
 *     U_{L,i-1/2} = Wrv(i-1);  U_{R,i-1/2} = Wlv(i)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   lr_states()          - computes L/R states
 *   lr_states_init()     - initializes memory for static global arrays
 *   lr_states_destruct() - frees memory for static global arrays
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef SECOND_ORDER

static Real **pW=NULL;


/*----------------------------------------------------------------------------*/
/* lr_states:
 * Input Arguments:
 *   U1d = CONSERVED variables at cell centers along 1-D slice
 *   Bx{c/i} = B in direction of slice at cell {center/interface}
 *   dt = timestep;   dtodx = dt/dx
 *   il,iu = lower and upper indices of zone centers in slice
 * U1d and Bxc must be initialized over [il-2:iu+2]
 * Bxi must be initialized over [il:iu+1]
 *
 * Output Arguments:
 *   Ul,Ur = L/R-states of CONSERVED variables at interfaces over [il:iu+1]
 */

void lr_states(const Prim1D W[], const Real Bxc[],
               const Real dt, const Real dtodx, const int il, const int iu,
               Prim1D Wl[], Prim1D Wr[])
{
  int i,n,m;
  Real pb,lim_slope1,lim_slope2,qa,qb,qc,qx;
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real dWc[NWAVE+NSCALARS],dWl[NWAVE+NSCALARS];
  Real dWr[NWAVE+NSCALARS],dWg[NWAVE+NSCALARS];
  Real dac[NWAVE+NSCALARS],dal[NWAVE+NSCALARS];
  Real dar[NWAVE+NSCALARS],dag[NWAVE+NSCALARS],da[NWAVE+NSCALARS];
  Real Wlv[NWAVE+NSCALARS],Wrv[NWAVE+NSCALARS];
  Real dW[NWAVE+NSCALARS],dWm[NWAVE+NSCALARS];
  Real *pWl, *pWr;

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

/*--- Step 1. ------------------------------------------------------------------
 * Transform to primitive variables over 1D slice, W=(d,Vx,Vy,Vz,[P],[By,Bz])
 */

/*
  for (i=il-2; i<=iu+2; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }
*/

  for (i=il-2; i<=iu+2; i++) pW[i] = (Real*)&(W[i]);
/*=============== START BIG LOOP OVER i ===============*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 2. ------------------------------------------------------------------
 * Compute eigensystem in primitive variables.  */

#ifdef HYDRO
#ifdef ISOTHERMAL
    esys_prim_iso_hyd(W[i].d,W[i].Vx,       ev,rem,lem);
#else
    esys_prim_adb_hyd(W[i].d,W[i].Vx,W[i].P,ev,rem,lem);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#ifdef ISOTHERMAL
    esys_prim_iso_mhd(W[i].d,W[i].Vx,       Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#else
    esys_prim_adb_mhd(W[i].d,W[i].Vx,W[i].P,Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#endif /* ISOTHERMAL */
#endif /* MHD */

/*--- Step 3. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dWc[n] = pW[i+1][n] - pW[i-1][n];
      dWl[n] = pW[i][n]   - pW[i-1][n];
      dWr[n] = pW[i+1][n] - pW[i][n];
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 4. ------------------------------------------------------------------
 * Project differences in primitive variables along characteristics */

    for (n=0; n<NWAVE; n++) {
      dac[n] = lem[n][0]*dWc[0];
      dal[n] = lem[n][0]*dWl[0];
      dar[n] = lem[n][0]*dWr[0];
      dag[n] = lem[n][0]*dWg[0];
      for (m=1; m<NWAVE; m++) {
	dac[n] += lem[n][m]*dWc[m];
	dal[n] += lem[n][m]*dWl[m];
	dar[n] += lem[n][m]*dWr[m];
	dag[n] += lem[n][m]*dWg[m];
      }
    }

/* Advected variables are treated differently; for them the right and left
 * eigenmatrices are simply the identitiy matrix.
 */
#if (NSCALARS > 0)
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      dac[n] = dWc[n];
      dal[n] = dWl[n];
      dar[n] = dWr[n];
      dag[n] = dWg[n];
    }
#endif

/*--- Step 5. ------------------------------------------------------------------
 * Apply monotonicity constraints to characteristic projections */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      da[n] = 0.0;
      if (dal[n]*dar[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dal[n]),fabs(dar[n]));
        lim_slope2 = MIN(0.5*fabs(dac[n]),fabs(dag[n]));
        da[n] = SIGN(dac[n])*MIN(2.0*lim_slope1,lim_slope2);
      }
    }

/*--- Step 6. ------------------------------------------------------------------
 * Project monotonic slopes in characteristic back to primitive variables  */

    for (n=0; n<NWAVE; n++) {
      dWm[n] = da[0]*rem[n][0];
      for (m=1; m<NWAVE; m++) {
        dWm[n] += da[m]*rem[n][m];
      }
    }

#if (NSCALARS > 0)
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      dWm[n] = da[n];
    }
#endif

/*--- Step 7. ------------------------------------------------------------------
 * Limit velocity difference to sound speed
 * Limit velocity so momentum is always TVD (using only minmod limiter)
 * CURRENTLY NOT USED.  Was added to make code more robust for turbulence
 * simulations, but found it added noise to Noh shocktube.
 */

#ifdef H_CORRECTION
/*
#ifdef ISOTHERMAL
    qa = Iso_csound;
#else
    qa = sqrt(Gamma*W[i].P/W[i].d);
#endif
    dWm[1] = SIGN(dWm[1])*MIN(fabs(dWm[1]),qa);
*/
#endif /* H_CORRECTION */
/*
    qa = W[i  ].Vx*W[i  ].d - W[i-1].Vx*W[i-1].d;
    qb = W[i+1].Vx*W[i+1].d - W[i  ].Vx*W[i  ].d;
    qc = W[i+1].Vx*W[i+1].d - W[i-1].Vx*W[i-1].d;
    qx = SIGN(qc)*MIN(2.0*MIN(fabs(qa),fabs(qb)), 0.5*fabs(qc));

    if ((-W[i].Vx*dWm[0]) > 0.0) {
      qa = 0.0;
      qb = -W[i].Vx*dWm[0];
    } else {
      qa = -W[i].Vx*dWm[0];
      qb = 0.0;
    }
    if (qx > 0.0) {
      qb += qx;
    } else {
      qa += qx;
    }
    qa = qa/W[i].d;
    qb = qb/W[i].d;

    dWm[1] = MIN(dWm[1],qb);
    dWm[1] = MAX(dWm[1],qa);
*/

/*--- Step 8. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = pW[i][n] - 0.5*dWm[n];
      Wrv[n] = pW[i][n] + 0.5*dWm[n];
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
    }

    for (n=0; n<(NWAVE+NSCALARS); n++) {
      dW[n] = Wrv[n] - Wlv[n];
    }

/*--- Step 9. ------------------------------------------------------------------
 * Integrate linear interpolation function over domain of dependence defined by
 * max(min) eigenvalue
 */

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

    qx = 0.5*MAX(ev[NWAVE-1],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWl[n] = Wrv[n] - qx*dW[n];
    }

    qx = -0.5*MIN(ev[0],0.0)*dtodx;
    for (n=0; n<(NWAVE+NSCALARS); n++) {
      pWr[n] = Wlv[n] + qx*dW[n];
    }

#ifndef THREED_VL /* include step 10 only if not using VL 3D integrator */
/* NB: If THREED_VL is defined, cannot run 1D or 2D problems; see integrate.c */

/*--- Step 10. -----------------------------------------------------------------
 * Then subtract amount of each wave n that does not reach the interface
 * during timestep (CW eqn 3.5ff).  For HLL fluxes, must subtract waves that
 * move in both directions.
 */

    for (n=0; n<NWAVE; n++) {
      if (ev[n] > 0.) {
	qa  = 0.0;
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*0.5*dtodx*(ev[NWAVE-1]-ev[n])*dW[m];
	}
	for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
        qa = 0.0;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*0.5*dtodx*(ev[n]-ev[0])*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] -= qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

    for (n=0; n<NWAVE; n++) {
      if (ev[n] < 0.) {
        qa = 0.0;
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*0.5*dtodx*(ev[0]-ev[n])*dW[m];
        }
        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
/* For HLL fluxes, subtract wave moving away from interface as well. */
#if defined(HLLE_FLUX) || defined(HLLC_FLUX) || defined(HLLD_FLUX)
	qa  = 0.0;
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*0.5*dtodx*(ev[n]-ev[NWAVE-1])*dW[m];
	}
	for (m=0; m<NWAVE; m++) pWl[m] -= qa*rem[m][n];
#endif /* HLL_FLUX */
      }
    }

/* Wave subtraction for advected quantities */
    for (n=NWAVE; n<(NWAVE+NSCALARS); n++) {
      if (W[i].Vx > 0.) {
        pWl[n] += 0.5*dtodx*(ev[NWAVE-1]-W[i].Vx)*dW[n];
      } else if (W[i].Vx < 0.) {
        pWr[n] += 0.5*dtodx*(ev[0]-W[i].Vx)*dW[n];
      }
    }

#endif /* THREED_VL */

  } /*=============== END BIG LOOP OVER i ===============*/

/*--- Step 11. -----------------------------------------------------------------
 * Convert back to conserved variables, and done.  */

/*
  for (i=il; i<=iu+1; i++) {
    pb = Prim1D_to_Cons1D(&Ul[i],&Wl[i],&Bxi[i]);
    pb = Prim1D_to_Cons1D(&Ur[i],&Wr[i],&Bxi[i]);
  }
*/

  return;
}

/*----------------------------------------------------------------------------*/
/* lr_states_init:  Allocate enough memory for work arrays */

void lr_states_init(int nx1, int nx2, int nx3)
{
  int i, nmax;
  nmax =  nx1 > nx2  ? nx1 : nx2;
  nmax = (nx3 > nmax ? nx3 : nmax) + 2*nghost;

  if ((pW = (Real**)malloc(nmax*sizeof(Real*))) == NULL) goto on_error;

  return;
  on_error:
    lr_states_destruct();
    ath_error("[lr_states_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* lr_states_destruct:  Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (pW != NULL) free(pW);
  return;
}

#endif /* SECOND_ORDER */
