#include "copyright.h"
/*==============================================================================
 * FILE: lr_states_tag.c
 *
 * PURPOSE: Third order (piecewise parabolic) spatial reconstruction using
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
 *   the left-and right-side of cell center.  Thus (see Step 19),
 *     U_{L,i-1/2} = Wrv(i-1);  U_{R,i-1/2} = Wlv(i)
 *
 * REFERENCE:
 *   Ask T. A. Gardiner
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

#if defined(THIRD_ORDER) && defined(TAG_ORDER)

#define NSCALARS 0
#define num_total (NWAVE + NSCALARS)

static Real Wlv[num_total], Wrv[num_total];
static Real dWc[num_total], d2Wc[num_total], dWl[num_total], dWr[num_total];
static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;

#undef num_total

#ifndef THREED_VL

#define num_total (NWAVE + NSCALARS)

static Real ev[num_total];
static Real rem[NWAVE][NWAVE], lem[NWAVE][NWAVE];
static Real dac[num_total], da2c[num_total], dal[num_total], dar[num_total];

#undef num_total

void lr_states(const Cons1D U1d[], const Real Bxc[], const Real Bxi[],
               const Real dt, const Real dtodx, const int il, const int iu,
               Cons1D Ul[], Cons1D Ur[]){
/*
void lr_states(const Prim1D W[], const Real Bxc[],
               const Real dt, const Real dtodx, const int il, const int iu,
               Prim1D Wl[], Prim1D Wr[]){
*/
  const int num_total = NWAVE + NSCALARS;
  Real *pWl=NULL, *pWc=NULL, *pWr=NULL;
  Real a, s, t, tt, abs_b, eps, pb;
/*  Real a, s, t, tt, abs_b, eps; */
  int i, n, m;

  for (n=0; n<NWAVE; ++n) {
    for (m=0; m<NWAVE; ++m) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

/*--- Step 1. ------------------------------------------------------------------
 * Transform to primitive variables over 1D slice, W=(d,Vx,Vy,Vz,[P],[By,Bz])
 */

  for (i=il-2; i<=iu+2; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }

/*=============== START BIG LOOP OVER i ===============*/

  for(i=il-1; i<=iu+1; ++i){

/*--- Step 2. ------------------------------------------------------------------
 * Initialize the array of primative variables for the left, center, and right
 * states */

    pWl = (Real*)&(W[i-1]);
    pWc = (Real*)&(W[i]);
    pWr = (Real*)&(W[i+1]);

/*--- Step 3. ------------------------------------------------------------------
 * Compute eigensystem in primitive variables.  */

#ifdef HYDRO
#ifdef ISOTHERMAL
  esys_prim_iso_hyd(W[i].d,W[i].Vx,       ev,rem,lem);
#else /* ADIABATIC */
  esys_prim_adb_hyd(W[i].d,W[i].Vx,W[i].P,ev,rem,lem);
#endif /* ISOTHERMAL || ADIABATIC */
#endif /* HYDRO */

#ifdef MHD
#ifdef ISOTHERMAL
  esys_prim_iso_mhd(W[i].d,W[i].Vx,       Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#else /* ADIABATIC */
  esys_prim_adb_mhd(W[i].d,W[i].Vx,W[i].P,Bxc[i],W[i].By,W[i].Bz,ev,rem,lem);
#endif /* ISOTHERMAL || ADIABATIC */
#endif /* MHD */

#if (NSCALARS > 0)
  /* Species modes propagate at the contact speed */
  for(n=NWAVE; n<num_total; ++n) ev[n] = pWc[1];
#endif /* NSCALARS > 0 */

/*--- Step 4. ------------------------------------------------------------------
 * Compute quadratic reconstruction in primitive variables and apply a first-
 * order limiter */

  for(n=0; n<num_total; ++n){
    /* Calculate a quadratic reconstruction in primitive variables */
    dWl[n] = (pWc[n] - pWl[n]);
    dWr[n] = (pWr[n] - pWc[n]);
    dWc[n] = 0.5*(pWr[n] - pWl[n]);
    d2Wc[n] = 0.5*(pWr[n] - 2.0*pWc[n] + pWl[n]);

    /* Now apply a first order limiter */
    eps = 1.0;

    /* Right Interface */
    s = dWc[n] + d2Wc[n]/3.0;

    /* The average slope is bounded between 0 and twice the right slope */
    t = 2.0*(pWr[n] - pWc[n]);

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Left Interface */
    s = dWc[n] - d2Wc[n]/3.0;

    /* The average slope is bounded between 0 and twice the left slope */
    t = 2.0*(pWc[n] - pWl[n]);

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Now limit dWc[n] & d2Wc[n] */
    dWc[n] *= eps;
    d2Wc[n] *= eps;
  }

/*--- Step 5. ------------------------------------------------------------------
 * In this block of code we ensure that intrinsically non-negative quantities
 * (density, pressure, species fraction) have a non-negative reconstruction.
 *   -- T. A. Gardiner -- April 27, 2007 */

  /* Check the density for a positive reconstruction */
  abs_b = fabs(dWc[0]);
  t = 1.5*abs_b - d2Wc[0]*13.0/6.0; /* x = +/- 3dx/2 constraint */

  if(3.0*fabs(d2Wc[0]) > abs_b){
    /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
    tt = 0.25*dWc[0]*dWc[0]/d2Wc[0] + d2Wc[0]/12.0;
    t = tt > t ? tt : t;
  }

  if(pWc[0] > 0.0){
    if(t > pWc[0]){
      eps = pWc[0]/t;
      fprintf(stderr,"Density positivity constraint: eps = %e\n",eps);
      dWc[0]  *= eps;
      d2Wc[0] *= eps;
    }
  }
  else
    dWc[0] = d2Wc[0] = 0.0; /* eps = 0.0 */

#ifndef ISOTHERMAL
  /* Check the pressure for a positive reconstruction */
  abs_b = fabs(dWc[4]);
  t = 1.5*abs_b - d2Wc[4]*13.0/6.0; /* x = +/- 3dx/2 constraint
*/

  if(3.0*fabs(d2Wc[4]) > abs_b){
    /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
    tt = 0.25*dWc[4]*dWc[4]/d2Wc[4] + d2Wc[4]/12.0;
    t = tt > t ? tt : t;
  }

  if(pWc[4] > 0.0){
    if(t > pWc[4]){
      eps = pWc[4]/t;
      fprintf(stderr,"Pressure positivity constraint: eps = %e\n",eps);
      dWc[4]  *= eps;
      d2Wc[4] *= eps;
    }
  }
  else
    dWc[4] = d2Wc[4] = 0.0; /* eps = 0.0 */
#endif /* ISOTHERMAL */

#if (NSCALARS > 0)
  /* Species Mass Density / Total Density */
  for(n=NWAVE; n<num_total; ++n){
    /* Check the species fraction for a positive reconstruction */
    abs_b = fabs(dWc[n]);
    t = 1.5*abs_b - d2Wc[n]*13.0/6.0; /* x = +/- 3dx/2 constraint */

    if(3.0*fabs(d2Wc[n]) > abs_b){
      /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
      tt = 0.25*dWc[n]*dWc[n]/d2Wc[n] + d2Wc[n]/12.0;
      t = tt > t ? tt : t;
    }

    if(pWc[n] > 0.0){
      if(t > pWc[n]){
	eps = pWc[n]/t;
        fprintf(stderr,"Species Fraction[%d] positivity constraint: eps = %e\n",
		 n-NWAVE,eps);
	dWc[n]  *= eps;
	d2Wc[n] *= eps;
      }
    }
    else
      dWc[n] = d2Wc[n] = 0.0; /* eps = 0.0 */
  }
#endif /* NSCALARS > 0 */

/*--- Step 6. ------------------------------------------------------------------
 * In trying to simulate driven turbulence I learned that the velocity field
 * can become relatively smooth and supersonic!  Allowing this results in BIG
 * velocity errors which kills the calculation.  Limiting the x-gradient of Vx
 * to be sub-sonic (compressions and rarefactions) solves this problem.
     -- T. A. Gardiner -- March 7, 2006 */

  /* Compute the sound speed */
#ifdef ADIABATIC
  a = sqrt((double)(Gamma*W[i].P/W[i].d));
#elif defined ISOTHERMAL
  a = Iso_csound;
#else
#error : [lr_states] : Unknown equation of state.
#endif

  a += a; /* Double the sound speeed */

  if(dWc[1] > a){ /* No transonic rarefactions interior to a cell */
    eps = a/dWc[1];
    /* printf("[lr_states Rarefaction]: a = %e, dWc[1] = %e, eps = %e\n",
       a,dWc[1],eps); */
    dWc[1] = a; /* dWc[1] *= eps; */
    d2Wc[1] *= eps;
  }
  else if(dWc[1] < -a){ /* No transonic compressions interior to a cell */
    eps = a/-dWc[1];
    /* printf("[lr_states Compression]: a = %e, dWc[1] = %e, eps = %e\n",
       a,dWc[1],eps); */
    dWc[1] = -a; /* dWc[1] *= eps; */
    d2Wc[1] *= a/-dWc[1];
  }

/*--- Step 7. ------------------------------------------------------------------
 * Calculate the associated characteristic reconstruction */

  for(n=0; n<NWAVE; ++n){
    dac[n] = lem[n][0]*dWc[0];
    da2c[n] = lem[n][0]*d2Wc[0];
    dal[n] = lem[n][0]*dWl[0];
    dar[n] = lem[n][0]*dWr[0];
    for (m=1; m<NWAVE; m++) {
      dac[n] += lem[n][m]*dWc[m];
      da2c[n] += lem[n][m]*d2Wc[m];
      dal[n] += lem[n][m]*dWl[m];
      dar[n] += lem[n][m]*dWr[m];
    }

    /* Now apply a first order limiter */
    eps = 1.0;

    /* Right Interface */
    s = dac[n] + da2c[n]/3.0;

    /* The average slope is bounded between 0 and twice the right slope */
    t = 2.0*dar[n];

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
        t /= s;
        eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
        t /= s;
        eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Left Interface */
    s = dac[n] - da2c[n]/3.0;

    /* The average slope is bounded between 0 and twice the left slope */
    t = 2.0*dal[n];

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
        t /= s;
        eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
        t /= s;
        eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Now limit dWc[n] & d2Wc[n] */
    dac[n] *= eps;
    da2c[n] *= eps;
  }

#if (NSCALARS > 0)
  /* The primitive and characteristic differences are identical for the
     species density fraction, pWc[num_waves+n] = (q.s[n]/q.d). */
  for(n=NWAVE; n<num_total; ++n){
    dac[n] = dWc[n];
    da2c[n] = d2Wc[n];
  }
#endif /* NSCALARS > 0 */

/*--- Step 8. ------------------------------------------------------------------
 * Average over the domain of dependence of this wave. */

  for(n=0; n<num_total; ++n){
    /* right state */
    dar[n] = 0.5*( 1.0 - dtodx*ev[n])*dac[n] +
      (1.0/6.0 - 0.5*dtodx*ev[n] + ev[n]*ev[n]*dtodx*dtodx/3.0)*da2c[n];

    /* left state */
    dal[n] = -0.5*(1.0 + dtodx*ev[n])*dac[n] +
      (1.0/6.0 + 0.5*dtodx*ev[n] + ev[n]*ev[n]*dtodx*dtodx/3.0)*da2c[n];
  }

/*--- Step 9. ------------------------------------------------------------------
 * Return to primitive variables */

  for(n=0; n<NWAVE; ++n){
    /* Calculate the left and right primitive, interface states,
       composed of each wave averaged over its respective domain of
       dependence. */
    Wlv[n] = pWc[n];
    Wrv[n] = pWc[n];
    for (m=0; m<NWAVE; m++) {
      Wlv[n] += dal[m]*rem[n][m];
      Wrv[n] += dar[m]*rem[n][m];
    }
  }

#if (NSCALARS > 0)
  /* The primitive and characteristic differences are identical for the
     species density fraction, pWc[num_waves+n] = (q.s[n]/q.d). */
  for(n=NWAVE; n<num_total; ++n){
    Wlv[n] = pWc[n] + dal[n];
    Wrv[n] = pWc[n] + dar[n];
  }
#endif /* NSCALARS > 0 */

/*--- Step 10. -----------------------------------------------------------------
 * Ensure that the density and pressure of the interface states is positive */

  if(!(Wlv[0] > 0.0) || !(Wrv[0] > 0.0)
#ifndef ISOTHERMAL
     || !(Wlv[4] > 0.0) || !(Wrv[4] > 0.0)
#endif /* ISOTHERMAL */
     ){
    fprintf(stderr,"[lr_states]: Non-physical states predicted\n");

    pWr = (Real*)&(Wr[i  ]); /* Wlv */
    pWl = (Real*)&(Wl[i+1]); /* Wrv */
    for(n=0; n<num_total; ++n){
      /* left state */
      pWr[n] = Wlv[n];

      /* right state */
      pWl[n] = Wrv[n];
    }

    /* If you want to dump some debug info, do it here */


    ath_error("[lr_states]: Terminating\n");
  }

/*--- Step 11. -----------------------------------------------------------------
 * Convert LR states arrays back to primitive data type */

  pWr = (Real*)&(Wr[i  ]); /* Wlv */
  pWl = (Real*)&(Wl[i+1]); /* Wrv */
  for(n=0; n<num_total; ++n){
    /* left state */
    pWr[n] = Wlv[n];

    /* right state */
    pWl[n] = Wrv[n];
  }

  } /*=============== END BIG LOOP OVER i ===============*/

/*--- Step 12. -----------------------------------------------------------------
 * Convert back to conserved variables, and done */

  for (i=il; i<=iu+1; i++) {
    pb = Prim1D_to_Cons1D(&Ul[i],&Wl[i],&Bxi[i]);
    pb = Prim1D_to_Cons1D(&Ur[i],&Wr[i],&Bxi[i]);
  }

  return;
}

#else /* THREED_VL */

void lr_states(const Cons1D U1d[], const Real Bxc[], const Real Bxi[],
               const Real dt, const Real dtodx, const int il, const int iu,
               Cons1D Ul[], Cons1D Ur[]){
/*
void lr_states(const Prim1D W[], const Real Bxc[],
               const Real dt, const Real dtodx, const int il, const int iu,
               Prim1D Wl[], Prim1D Wr[]){
*/

  const int num_total = NWAVE + NSCALARS;
  Real *pWl=NULL, *pWc=NULL, *pWr=NULL;
  Real a, s, t, tt, abs_b, eps, pb;
/*  Real a, s, t, tt, abs_b, eps; */
  int i, n, m;

/*--- Step 1. ------------------------------------------------------------------
 * Transform to primitive variables over 1D slice, W=(d,Vx,Vy,Vz,[P],[By,Bz])
 */

  for (i=il-2; i<=iu+2; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }

/*=============== START BIG LOOP OVER i ===============*/

  for(i=il-1; i<=iu+1; ++i){

/*--- Step 2. ------------------------------------------------------------------
 * Initialize the array of primative variables for the left, center, and right
 * states */

    pWl = (Real*)&(W[i-1]);
    pWc = (Real*)&(W[i]);
    pWr = (Real*)&(W[i+1]);

/*--- Step 3. ------------------------------------------------------------------
 * NO NEED TO Compute eigensystem in primitive variables.  */

/*--- Step 4. ------------------------------------------------------------------
 * Compute quadratic reconstruction in primitive variables and apply a first-
 * order limiter */

  for(n=0; n<num_total; ++n){
    /* Calculate a quadratic reconstruction in primitive variables */
    dWl[n] = (pWc[n] - pWl[n]);
    dWr[n] = (pWr[n] - pWc[n]);
    dWc[n] = 0.5*(pWr[n] - pWl[n]);
    d2Wc[n] = 0.5*(pWr[n] - 2.0*pWc[n] + pWl[n]);

    /* Now apply a first order limiter */
    eps = 1.0;

    /* Right Interface */
    s = dWc[n] + d2Wc[n]/3.0;

    /* The average slope is bounded between 0 and twice the right slope */
    t = 2.0*(pWr[n] - pWc[n]);

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Left Interface */
    s = dWc[n] - d2Wc[n]/3.0;

    /* The average slope is bounded between 0 and twice the left slope */
    t = 2.0*(pWc[n] - pWl[n]);

    if(s > 0.0){
      if(t > 0.0){ /* bounds are 0 <= s <= t */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }
    else if(s < 0.0){
      if(t < 0.0){ /* bounds are t <= s <= 0 */
	t /= s;
	eps = eps < t ? eps : t;
      }
      else eps = 0.0;
    }

    /* Now limit dWc[n] & d2Wc[n] */
    dWc[n] *= eps;
    d2Wc[n] *= eps;
  }

/*--- Step 5. ------------------------------------------------------------------
 * In this block of code we ensure that intrinsically non-negative quantities
 * (density, pressure, species fraction) have a non-negative reconstruction.
 *   -- T. A. Gardiner -- April 27, 2007 */

  /* Check the density for a positive reconstruction */
  abs_b = fabs(dWc[0]);
  t = 1.5*abs_b - d2Wc[0]*13.0/6.0; /* x = +/- 3dx/2 constraint */

  if(3.0*fabs(d2Wc[0]) > abs_b){
    /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
    tt = 0.25*dWc[0]*dWc[0]/d2Wc[0] + d2Wc[0]/12.0;
    t = tt > t ? tt : t;
  }

  if(pWc[0] > 0.0){
    if(t > pWc[0]){
      eps = pWc[0]/t;
      fprintf(stderr,"Density positivity constraint: eps = %e\n",eps);
      dWc[0]  *= eps;
      d2Wc[0] *= eps;
    }
  }
  else
    dWc[0] = d2Wc[0] = 0.0; /* eps = 0.0 */

#ifndef ISOTHERMAL
  /* Check the pressure for a positive reconstruction */
  abs_b = fabs(dWc[4]);
  t = 1.5*abs_b - d2Wc[4]*13.0/6.0; /* x = +/- 3dx/2 constraint
*/

  if(3.0*fabs(d2Wc[4]) > abs_b){
    /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
    tt = 0.25*dWc[4]*dWc[4]/d2Wc[4] + d2Wc[4]/12.0;
    t = tt > t ? tt : t;
  }

  if(pWc[4] > 0.0){
    if(t > pWc[4]){
      eps = pWc[4]/t;
      fprintf(stderr,"Pressure positivity constraint: eps = %e\n",eps);
      dWc[4]  *= eps;
      d2Wc[4] *= eps;
    }
  }
  else
    dWc[4] = d2Wc[4] = 0.0; /* eps = 0.0 */
#endif /* ISOTHERMAL */

#if (NSCALARS > 0)
  /* Species Mass Density / Total Density */
  for(n=NWAVE; n<num_total; ++n){
    /* Check the species fraction for a positive reconstruction */
    abs_b = fabs(dWc[n]);
    t = 1.5*abs_b - d2Wc[n]*13.0/6.0; /* x = +/- 3dx/2 constraint */

    if(3.0*fabs(d2Wc[n]) > abs_b){
      /* there's an extremum in the domain -3dx/2 < x < 3dx/2 */
      tt = 0.25*dWc[n]*dWc[n]/d2Wc[n] + d2Wc[n]/12.0;
      t = tt > t ? tt : t;
    }

    if(pWc[n] > 0.0){
      if(t > pWc[n]){
	eps = pWc[n]/t;
        fprintf(stderr,"Species Fraction[%d] positivity constraint: eps = %e\n",
		 n-NWAVE,eps);
	dWc[n]  *= eps;
	d2Wc[n] *= eps;
      }
    }
    else
      dWc[n] = d2Wc[n] = 0.0; /* eps = 0.0 */
  }
#endif /* NSCALARS > 0 */

/*--- Step 6. ------------------------------------------------------------------
 * In trying to simulate driven turbulence I learned that the velocity field
 * can become relatively smooth and supersonic!  Allowing this results in BIG
 * velocity errors which kills the calculation.  Limiting the x-gradient of Vx
 * to be sub-sonic (compressions and rarefactions) solves this problem.
     -- T. A. Gardiner -- March 7, 2006 */

  /* Compute the sound speed */
#ifdef ADIABATIC
  a = sqrt((double)(Gamma*W[i].P/W[i].d));
#elif defined ISOTHERMAL
  a = Iso_csound;
#else
#error : [lr_states] : Unknown equation of state.
#endif

  a += a; /* Double the sound speeed */

  if(dWc[1] > a){ /* No transonic rarefactions interior to a cell */
    eps = a/dWc[1];
    /* printf("[lr_states Rarefaction]: a = %e, dWc[1] = %e, eps = %e\n",
       a,dWc[1],eps); */
    dWc[1] = a; /* dWc[1] *= eps; */
    d2Wc[1] *= eps;
  }
  else if(dWc[1] < -a){ /* No transonic compressions interior to a cell */
    eps = a/-dWc[1];
    /* printf("[lr_states Compression]: a = %e, dWc[1] = %e, eps = %e\n",
       a,dWc[1],eps); */
    dWc[1] = -a; /* dWc[1] *= eps; */
    d2Wc[1] *= a/-dWc[1];
  }

/*--- Step 7. ------------------------------------------------------------------
 * NO NEED TO Calculate the associated characteristic reconstruction */

/*--- Step 8. ------------------------------------------------------------------
 * NO NEED TO Average over the domain of dependence of this wave. */

/*--- Step 9. ------------------------------------------------------------------
 * Return to primitive variables */

  for(n=0; n<num_total; ++n){
    /* Calculate the left and right primitive, interface states,
       composed of each wave averaged over its respective domain of
       dependence. */
    /* right state */
    Wrv[n] = pWc[n] + 0.5*dWc[n] + (1.0/6.0)*d2Wc[n];

    /* left state */
    Wlv[n] = pWc[n] - 0.5*dWc[n] + (1.0/6.0)*d2Wc[n];
  }

/*--- Step 10. -----------------------------------------------------------------
 * Ensure that the density and pressure of the interface states is positive */

  if(!(Wlv[0] > 0.0) || !(Wrv[0] > 0.0)
#ifndef ISOTHERMAL
     || !(Wlv[4] > 0.0) || !(Wrv[4] > 0.0)
#endif /* ISOTHERMAL */
     ){
    fprintf(stderr,"[lr_states]: Non-physical states predicted\n");

    pWr = (Real*)&(Wr[i  ]); /* Wlv */
    pWl = (Real*)&(Wl[i+1]); /* Wrv */
    for(n=0; n<num_total; ++n){
      /* left state */
      pWr[n] = Wlv[n];

      /* right state */
      pWl[n] = Wrv[n];
    }

    /* If you want to dump some debug info, do it here */


    ath_error("[lr_states]: Terminating\n");
  }

/*--- Step 11. -----------------------------------------------------------------
 * Convert LR states arrays back to primative data type */

  pWr = (Real*)&(Wr[i  ]); /* Wlv */
  pWl = (Real*)&(Wl[i+1]); /* Wrv */
  for(n=0; n<num_total; ++n){
    /* left state */
    pWr[n] = Wlv[n];

    /* right state */
    pWl[n] = Wrv[n];
  }

  } /*=============== END BIG LOOP OVER i ===============*/

/*--- Step 12. -----------------------------------------------------------------
 * Convert back to conserved variables, and done */

  for (i=il; i<=iu+1; i++) {
    pb = Prim1D_to_Cons1D(&Ul[i],&Wl[i],&Bxi[i]);
    pb = Prim1D_to_Cons1D(&Ur[i],&Wr[i],&Bxi[i]);
  }

  return;
}

#endif /* THREED_VL */

/*----------------------------------------------------------------------------*/
/* lr_states_init:  Allocate enough memory for work arrays */

void lr_states_init(int nx1, int nx2, int nx3)
{
  int i, nmax;
  nmax =  nx1 > nx2  ? nx1 : nx2;
  nmax = (nx3 > nmax ? nx3 : nmax) + 2*nghost;

  if ((W = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr = (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;

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
  return;
}



#endif /* THIRD_ORDER */
