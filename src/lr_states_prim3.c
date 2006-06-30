#include "copyright.h"
/*==============================================================================
 * FILE: lr_states_prim3.c
 *
 * PURPOSE: Third order (piecewise parabolic) spatial reconstruction using
 *   characteristic interpolation in the primitive variables.  A time-evolution
 *   (characteristic tracing) step is used to interpolate interface values to
 *   the half time level {n+1/2}, unless the unsplit integrator in 3D is VL.
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 *   U_{L,i-1/2} is denoted by Ul[i  ];   U_{R,i-1/2} is denoted by Ur[i  ]
 *   U_{L,i+1/2} is denoted by Ul[i+1];   U_{R,i+1/2} is denoted by Ur[i+1]
 *
 * SOURCE TERMS: are added to the primitive variables
 *
 * REFERENCE:
 *   P. Colella & P. Woodward, "The piecewise parabolic method (PPM) for
 *   gas-dynamical simulations", JCP, 54, 174 (1984).
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
#include "prototypes.h"

#ifdef THIRD_ORDER

static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;
static Real **pW=NULL, **dWm=NULL, **Wim1h=NULL;

/*----------------------------------------------------------------------------*/
/* lr_states:
 * Input Arguments:
 *   U1d = CONSERVED variables at cell centers along 1-D slice
 *   Bxc = B in direction of slice at cell centers
 *   dt = timestep
 *   dtodx = dt/dx
 *   is,ie = starting and ending indices of zone centers in slice
 * U1d and Bxc must be initialized over [is-nghost:ie+nghost]
 *
 * Output Arguments:
 *   Ul,Ur = L/R-states of CONSERVED variables at interfaces over [is:ie+1]
 */

void lr_states(const Cons1D U1d[], const Real Bxc[], const Real Bxi[],
	       const Real dt, const Real dtodx, const int is, const int ie,
	       const Prim1D *Wsrc, Cons1D Ul[], Cons1D Ur[])
{
  int i,n,m;
  Real maxevlr=0.0, pb,lim_slope1,lim_slope2,qa,qb,qc,qx;
  Real ev    [NWAVE],rem    [NWAVE][NWAVE],lem    [NWAVE][NWAVE];
  Real ev_ip1[NWAVE],rem_ip1[NWAVE][NWAVE],lem_ip1[NWAVE][NWAVE];
  Real dWc[NWAVE],dWl[NWAVE],dWr[NWAVE],dWg[NWAVE];
  Real dac[NWAVE],dal[NWAVE],dar[NWAVE],dag[NWAVE],da[NWAVE];
  Real Wlv[NWAVE],Wrv[NWAVE],dW[NWAVE],W6[NWAVE];
  Real *pWl, *pWr, *pWsrc;

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
      rem_ip1[n][m] = 0.0;
      lem_ip1[n][m] = 0.0;
    }
  }

/*--- Step 1. ------------------------------------------------------------------
 * Transform to primitive variables over 1D slice, W=(d,Vx,Vy,Vz,[P],[By,Bz])
 */

  for (i=is-3; i<=ie+3; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }

/*=============== START LOOP OVER is-2:is-1 ===============*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute eigensystem in primitive variables.  */

  for (i=is-2; i<=is-1; i++) {

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

    for (n=0; n<NWAVE; n++) {
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

/*--- Step 5. ------------------------------------------------------------------
 * Apply monotonicity constraints to characteristic projections */

    for (n=0; n<NWAVE; n++) {
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
      dWm[i][n] = da[0]*rem[n][0];
      for (m=1; m<NWAVE; m++) {
        dWm[i][n] += da[m]*rem[n][m];
      }
    }

/*--- Step 7. ------------------------------------------------------------------
 * Limit velocity difference across cell to sound speed, limit velocity slope
 * so momentum is always TVD (using only minmod limiter) */

#ifdef ISOTHERMAL
    qa = Iso_csound;
#else
    qa = sqrt(Gamma*W[i].P/W[i].d);
#endif
    dWm[i][1] = SIGN(dWm[i][1])*MIN(fabs(dWm[i][1]),qa);

    qa = U1d[i  ].Mx - U1d[i-1].Mx;
    qb = U1d[i+1].Mx - U1d[i  ].Mx;
    qc = U1d[i+1].Mx - U1d[i-1].Mx;
    qx = SIGN(qc)*MIN(2.0*MIN(fabs(qa),fabs(qb)), 0.5*fabs(qc));

    if ((-W[i].Vx*dWm[i][0]) > 0.0) {
      qa = 0.0;
      qb = -W[i].Vx*dWm[i][0];
    } else {
      qa = -W[i].Vx*dWm[i][0];
      qb = 0.0;
    }
    if (qx > 0.0) {
      qb += qx;
    } else {
      qa += qx;
    }
    qa = qa/W[i].d;
    qb = qb/W[i].d;

    dWm[i][1] = MIN(dWm[i][1],qb);
    dWm[i][1] = MAX(dWm[i][1],qa);

  }
/*=============== END LOOP OVER is-2:is-1 ===============*/


/*--- Step 8. ------------------------------------------------------------------
 * Construct parabolic interpolant in primitive variables at left-interface
 * of cell is-1 ("W[is-3/2]", CW eqn 1.6) using linear TVD slopes at is-2 and
 * is-1 computed in Steps 2-7.
 */

  for (n=0; n<NWAVE; n++) {
    Wim1h[is-1][n] =.5*(pW[is-1][n]+pW[is-2][n])-(dWm[is-1][n]-dWm[is-2][n])/6.;
  }

/* Loop over is-2:is-1 in Steps 2-7 was necessary to bootstrap method by
 * computing Wim1h[is-1].  Now repeat these steps for rest of grid.  Splitting
 * into two loops avoids calculating eigensystems twice per cell.
 *
 * At the start of the loop, rem and lem still store values at i=is-1 computed
 * at the end of Step 2.  For each i, the eigensystem at i+1 is stored in
 * rem_ip1 and lem_ip1.  At the end of the loop rem[lem] is then set to
 * rem_ip1[lem_ip1] in preparation for the next iteration.
 */

/*=============== START BIG LOOP OVER i ===============*/
/* Steps 9-15 below are identical to steps 2-8 above */
  for (i=is-1; i<=ie+1; i++) {

/*--- Step 9. ------------------------------------------------------------------
 * Compute eigensystem in primitive variables.  */

#ifdef HYDRO
#ifdef ISOTHERMAL
    esys_prim_iso_hyd(W[i+1].d,W[i+1].Vx,         ev_ip1,rem_ip1,lem_ip1);
#else
    esys_prim_adb_hyd(W[i+1].d,W[i+1].Vx,W[i+1].P,ev_ip1,rem_ip1,lem_ip1);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#ifdef ISOTHERMAL
    esys_prim_iso_mhd(W[i+1].d,W[i+1].Vx,
		       Bxc[i+1],W[i+1].By,W[i+1].Bz,ev_ip1,rem_ip1,lem_ip1);
#else
    esys_prim_adb_mhd(W[i+1].d,W[i+1].Vx,W[i+1].P,
		       Bxc[i+1],W[i+1].By,W[i+1].Bz,ev_ip1,rem_ip1,lem_ip1);
#endif /* ISOTHERMAL */
#endif /* MHD */

/*--- Step 10. -----------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables */

    for (n=0; n<NWAVE; n++) {
      dWc[n] = pW[i+2][n] - pW[i][n];
      dWl[n] = pW[i+1][n] - pW[i][n];
      dWr[n] = pW[i+2][n] - pW[i+1][n];
      if (dWl[n]*dWr[n] > 0.0) {
        dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      } else {
        dWg[n] = 0.0;
      }
    }

/*--- Step 11. -----------------------------------------------------------------
 * Project differences in primitive variables along characteristics */

    for (n=0; n<NWAVE; n++) {
      dac[n] = lem_ip1[n][0]*dWc[0];
      dal[n] = lem_ip1[n][0]*dWl[0];
      dar[n] = lem_ip1[n][0]*dWr[0];
      dag[n] = lem_ip1[n][0]*dWg[0];
      for (m=1; m<NWAVE; m++) {
        dac[n] += lem_ip1[n][m]*dWc[m];
        dal[n] += lem_ip1[n][m]*dWl[m];
        dar[n] += lem_ip1[n][m]*dWr[m];
        dag[n] += lem_ip1[n][m]*dWg[m];
      }
    }

/*--- Step 12. -----------------------------------------------------------------
 * Apply monotonicity constraints to characteristic projections */

    for (n=0; n<NWAVE; n++) {
      da[n] = 0.0;
      if (dal[n]*dar[n] > 0.0) {
        lim_slope1 = MIN(    fabs(dal[n]),fabs(dar[n]));
        lim_slope2 = MIN(0.5*fabs(dac[n]),fabs(dag[n]));
        da[n] = SIGN(dac[n])*MIN(2.0*lim_slope1,lim_slope2);
      }
    }

/*--- Step 13. -----------------------------------------------------------------
 * Project monotonic slopes in characteristic back to primitive variables  */

    for (n=0; n<NWAVE; n++) {
      dWm[i+1][n] = da[0]*rem_ip1[n][0];
      for (m=1; m<NWAVE; m++) {
        dWm[i+1][n] += da[m]*rem_ip1[n][m];
      }
    }

/*--- Step 14. -----------------------------------------------------------------
 * Limit velocity difference across cell to sound speed, limit velocity slope
 * so momentum is always TVD (using only minmod limiter) */

#ifdef ISOTHERMAL
    qa = Iso_csound;
#else
    qa = sqrt(Gamma*W[i+1].P/W[i+1].d);
#endif
    dWm[i+1][1] = SIGN(dWm[i+1][1])*MIN(fabs(dWm[i+1][1]),qa);

    qa = U1d[i+1].Mx - U1d[i  ].Mx;
    qb = U1d[i+2].Mx - U1d[i+1].Mx;
    qc = U1d[i+2].Mx - U1d[i  ].Mx;
    qx = SIGN(qc)*MIN(2.0*MIN(fabs(qa),fabs(qb)), 0.5*fabs(qc));

    if ((-W[i+1].Vx*dWm[i+1][0]) > 0.0) {
      qa = 0.0;
      qb = -W[i+1].Vx*dWm[i+1][0];
    } else {
      qa = -W[i+1].Vx*dWm[i+1][0];
      qb = 0.0;
    }
    if (qx > 0.0) {
      qb += qx;
    } else {
      qa += qx;
    }
    qa = qa/W[i+1].d;
    qb = qb/W[i+1].d;

    dWm[i+1][1] = MIN(dWm[i+1][1],qb);
    dWm[i+1][1] = MAX(dWm[i+1][1],qa);

/*--- Step 15. -----------------------------------------------------------------
 * Construct parabolic interpolant in primitive variables at left-interface
 * of cell i+1 ("W[i+1/2]", CW eqn 1.6) using linear TVD slopes at i and
 * i+1 computed in Steps 2-7.
 */

    for (n=0; n<NWAVE; n++) {
      Wim1h[i+1][n] = 0.5*(pW[i+1][n]+pW[i][n]) - (dWm[i+1][n]-dWm[i][n])/6.0;
    }

/*--- Step 16. -----------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

    for (n=0; n<NWAVE; n++) {
      Wlv[n] = Wim1h[i  ][n];
      Wrv[n] = Wim1h[i+1][n];
    }

    for (n=0; n<NWAVE; n++) {
      Wlv[n] = MAX(MIN(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wlv[n] = MIN(MAX(pW[i][n],pW[i-1][n]),Wlv[n]);
      Wrv[n] = MAX(MIN(pW[i][n],pW[i+1][n]),Wrv[n]);
      Wrv[n] = MIN(MAX(pW[i][n],pW[i+1][n]),Wrv[n]);
    }

/*--- Step 17. -----------------------------------------------------------------
 * Monotonize again (CW eqn 1.10) */

    for (n=0; n<NWAVE; n++) {
      qa = (Wrv[n]-pW[i][n])*(pW[i][n]-Wlv[n]);
      qb = Wrv[n]-Wlv[n];
      qc = 6.0*(pW[i][n] - 0.5*(Wlv[n]+Wrv[n]));
      if (qa <= 0.0) {
        Wlv[n] = pW[i][n];
        Wrv[n] = pW[i][n];
      } else if (qb*(qb - qc) < 0.0) {
        Wlv[n] = 3.0*pW[i][n] - 2.0*Wrv[n];
      } else if (qb*(qb + qc) < 0.0) {
        Wrv[n] = 3.0*pW[i][n] - 2.0*Wlv[n];
      }

/*--- Step 18. -----------------------------------------------------------------
 * Compute coefficients of interpolation parabolae (CW eqn 1.5) */

      dW[n] = Wrv[n] - Wlv[n];
      W6[n] = 6.0*(pW[i][n] - 0.5*(Wlv[n] + Wrv[n]));
    }

#ifdef CTU_INTEGRATOR /* include steps 19-21 only if using CTU 3D integrator */
/*--- Step 19. -----------------------------------------------------------------
 * Integrate linear interpolation function over domain of dependence defined by
 * max(min) eigenvalue (CW eqn 1.12)
 */

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

    qx = TWO_3RDS*MAX(ev[NWAVE-1],0.0)*dtodx;
    for (n=0; n<NWAVE; n++) {
      pWl[n] = Wrv[n] - 0.75*qx*(dW[n] - (1.0 - qx)*W6[n]);
    }

    qx = -TWO_3RDS*MIN(ev[0],0.0)*dtodx;
    for (n=0; n<NWAVE; n++) {
      pWr[n] = Wlv[n] + 0.75*qx*(dW[n] + (1.0 - qx)*W6[n]);
    }

/*--- Step 20. -----------------------------------------------------------------
 * Then subtract amount of each wave m that does not reach the interface
 * during timestep (CW eqn 3.5ff)
 */

    for (n=0; n<NWAVE-1; n++) {
      if (ev[n] > 0.) {
	qa  = 0.0;
        qb = 0.5*dtodx*(ev[NWAVE-1]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[NWAVE-1]*ev[NWAVE-1] - ev[n]*ev[n]);
	for (m=0; m<NWAVE; m++) {
	  qa += lem[n][m]*(qb*(dW[m]-W6[m]) + qc*W6[m]);
	}
	for (m=0; m<NWAVE; m++) pWl[m] += qa*rem[m][n];
      }
    }

    for (n=1; n<NWAVE; n++) {
      if (ev[n] < 0.) {
        qa = 0.0;
        qb = 0.5*dtodx*(ev[0]-ev[n]);
        qc = 0.5*dtodx*dtodx*TWO_3RDS*(ev[0]*ev[0] - ev[n]*ev[n]);
        for (m=0; m<NWAVE; m++) {
          qa += lem[n][m]*(qb*(dW[m]+W6[m]) + qc*W6[m]);
        }
        for (m=0; m<NWAVE; m++) pWr[m] += qa*rem[m][n];
      }
    }

/*--- Step 21. -----------------------------------------------------------------
 * Add the source terms to the appropriate interface states */

    if(Wsrc != NULL){
      pWsrc = (Real *)&(Wsrc[i]);
      for (n=0; n<NWAVE; n++) {
	pWl[n] += 0.5*dt*pWsrc[n];
	pWr[n] += 0.5*dt*pWsrc[n];
      }
    }
#endif /* CTU_INTEGRATOR */

/*--- Step 22. -----------------------------------------------------------------
 * Save eigenvalues and eigenmatrices at i+1 for use in next iteration */

    for (m=0; m<NWAVE; m++) {
      ev[m] = ev_ip1[m];
      for (n=0; n<NWAVE; n++) {
        rem[m][n] = rem_ip1[m][n];
        lem[m][n] = lem_ip1[m][n];
      }
    }

  } /*=============== END BIG LOOP OVER i ===============*/

/*--- Step 23. -----------------------------------------------------------------
 * Convert back to conserved variables, and done
 */

  for (i=is; i<=ie+1; i++) {
    pb = Prim1D_to_Cons1D(&Ul[i],&Wl[i],&Bxi[i]);
    pb = Prim1D_to_Cons1D(&Ur[i],&Wr[i],&Bxi[i]);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* lr_states_init:  Allocate enough memory for work arrays */

void lr_states_init(int nx1, int nx2, int nx3)
{
  int i, nmax;
  nmax =  nx1 > nx2  ? nx1 : nx2;
  nmax = (nx3 > nmax ? nx3 : nmax) + 2*nghost;

  if ((W  = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;

  if ((pW = (Real**)malloc(nmax*sizeof(Real*))) == NULL) goto on_error;
  for (i=0; i<nmax; i++) pW[i] = (Real*)&(W[i]);

  if ((dWm = (Real**)calloc_2d_array(nmax, NWAVE, sizeof(Real))) == NULL)
    goto on_error;

  if ((Wim1h = (Real**)calloc_2d_array(nmax, NWAVE, sizeof(Real))) == NULL)
    goto on_error;

  return;
  on_error:
    lr_states_destruct();
    ath_error("[lr_states_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* lr_states_destruct:  Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (W != NULL) free(W);
  if (Wl != NULL) free(Wl);
  if (Wr != NULL) free(Wr);
  if (pW != NULL) free(pW);
  if (dWm != NULL) free_2d_array((void**)dWm);
  if (Wim1h != NULL) free_2d_array((void**)Wim1h);
  return;
}

#endif /* THIRD_ORDER */
