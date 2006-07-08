#include "copyright.h"
/*==============================================================================
 * FILE: flux_hlld.c
 *
 * PURPOSE: Computes 1D fluxes using the HLLD Riemann solver, and extension of
 *   the HLLE solver to MHD.  Only works for MHD problems.
 *
 * REFERENCES:
 *   T. Miyoshi & K. Kusano, "A multi-state HLL approximate Riemann solver
 *   for ideal MHD", JCP, 208, 315 (2005)
 *
 * HISTORY: written by Brian Biskeborn, May 8, 2006, COS sophmore project.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   flux_hlld()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error : The HLLD flux only works with adiabatic EOS.
#endif /* ISOTHERMAL */
#ifndef MHD
#error : The HLLD flux only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* flux_hlld:
 * Input Arguments:
 *   Bxi = B in direction of slice at cell interface
 *   Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *
 * Output Arguments:
 *   Flux = fluxes of CONSERVED variables at cell interface
 */

void flux_hlld(const Real Bxi, const Cons1D Ul, const Cons1D Ur, Cons1D *pFlux)
{
  Cons1D Ulst,Uldst,Urdst,Urst;       /* Conserved variable for all states */
  Prim1D Wl,Wlst,Wldst,Wrdst,Wrst,Wr; /* Primitive variables for all states */
  Cons1D Fl,Fr;                         /* Fluxes for left & right states */
  Real spd[5],maxspd;                 /* signal speeds, left to right */
  Real sdl,sdr,sdml,sdmr;             /* S_i-u_i, S_i-S_M (i=L or R) */
  Real pbl,pbr,pbstl,pbstr,pbdstl,pbdstr; /* Magnetic pressures */
  Real cfl,cfr,cfmax;                 /* Cf (left & right), max(cfl,cfr) */
  Real gpl,gpr,gpbl,gpbr;             /* gamma*P, gamma*P + B */
  Real sqrtdl,sqrtdr;                 /* sqrt of the L* & R* densities */
  Real invsumd;                       /* 1/(sqrtdl + sqrtdr) */
  Real ptl,ptr,ptst;                  /* total pressures */
  Real vbstl,vbstr;                   /* v_i* dot B_i* for i=L or R */
  Real Bxsig;                         /* sign(Bx) = 1 for Bx>0, -1 for Bx<0 */
  Real Bxsq = SQR(Bxi);               /* Bx^2 */
  Real tmp;                      /* Temp variable for repeated calculations */

/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);

/*--- Step 2. ------------------------------------------------------------------
 * Compute left & right wave speeds according to Miyoshi & Kusano, eqn. (67)
 */

  gpl  = Gamma * Wl.P;
  gpr  = Gamma * Wr.P;
  gpbl = gpl + 2.0*pbl;
  gpbr = gpr + 2.0*pbr;

  cfl = sqrt((gpbl + sqrt(SQR(gpbl)-4*gpl*Bxsq))/(2.0*Wl.d));
  cfr = sqrt((gpbr + sqrt(SQR(gpbr)-4*gpr*Bxsq))/(2.0*Wr.d));
  cfmax = MAX(cfl,cfr);

  if(Wl.Vx <= Wr.Vx) {
    spd[0] = Wl.Vx - cfmax;
    spd[4] = Wr.Vx + cfmax;
  }
  else {
    spd[0] = Wr.Vx - cfmax;
    spd[4] = Wl.Vx + cfmax;
  }

  maxspd = MAX(fabs(spd[0]),fabs(spd[4]));

  if(Ul.d  == Ur.d &&
     Ul.Mx == Ur.Mx &&
     Ul.My == Ur.My &&
     Ul.Mz == Ur.Mz &&
     Ul.E  == Ur.E &&
     Ul.By == Ur.By &&
     Ul.Bz == Ur.Bz) {
/* return Fl (= Fr) */
    pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
    ptl = Wl.P + pbl;
    pFlux->d  = Ul.Mx;
    pFlux->Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
    pFlux->My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
    pFlux->Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
    pFlux->E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
    pFlux->By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
    pFlux->Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;
    return;
  }

  if(spd[0] >= 0.0) {
/* return Fl */
    ptl = Wl.P + pbl;
    pFlux->d  = Ul.Mx;
    pFlux->Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
    pFlux->My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
    pFlux->Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
    pFlux->E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
    pFlux->By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
    pFlux->Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;
    return;
  }

  if(spd[4] <= 0.0) {
/* return Fr */
    ptr = Wr.P + pbr;
    pFlux->d  = Ur.Mx;
    pFlux->Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
    pFlux->My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
    pFlux->Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
    pFlux->E  = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
    pFlux->By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
    pFlux->Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;
    return;
  }

/*--- Step 3. ------------------------------------------------------------------
 * Compute middle and Alfven wave speeds
 */

  sdl = spd[0] - Wl.Vx;
  sdr = spd[4] - Wr.Vx;

  ptl = Wl.P + pbl;
  ptr = Wr.P + pbr;

  spd[2] = (sdr*Wr.d*Wr.Vx - sdl*Wl.d*Wl.Vx - ptr + ptl) /
           (sdr*Wr.d-sdl*Wl.d);

  sdml   = spd[0] - spd[2];
  sdmr   = spd[4] - spd[2];
  Ulst.d = Ul.d * sdl/sdml;
  Urst.d = Ur.d * sdr/sdmr;
  sqrtdl = sqrt(Ulst.d);
  sqrtdr = sqrt(Urst.d);

  spd[1] = spd[2] - fabs(Bxi)/sqrtdl;
  spd[3] = spd[2] + fabs(Bxi)/sqrtdr;

/*--- Step 4. ------------------------------------------------------------------
 * Compute intermediate states
 */

  ptst = ptl + Ul.d*sdl*(sdl-sdml);

/* Ul* */
  Ulst.Mx = Ulst.d * spd[2];
  if(spd[2] == Wl.Vx) {
    Ulst.My = Ulst.d * Wl.Vy;
    Ulst.Mz = Ulst.d * Wl.Vz;
  }
  else {
    tmp = Bxi*(sdl-sdml)/(Ul.d*sdl*sdml-Bxsq);
    Ulst.My = Ulst.d * (Wl.Vy - Ul.By*tmp);
    Ulst.Mz = Ulst.d * (Wl.Vz - Ul.Bz*tmp);
  }
  if(Ul.By == 0.0 && Ul.Bz == 0.0) {
    Ulst.By = 0.0;
    Ulst.Bz = 0.0;
  }
  else {
    tmp = (Ul.d*SQR(sdl)-Bxsq)/(Ul.d*sdl*sdml - Bxsq);
    Ulst.By = Ul.By * tmp;
    Ulst.Bz = Ul.Bz * tmp;
  }
  vbstl = (Ulst.Mx*Bxi+Ulst.My*Ulst.By+Ulst.Mz*Ulst.Bz)/Ulst.d;
  Ulst.E = (sdl*Ul.E - ptl*Wl.Vx + ptst*spd[2] +
            Bxi*(Wl.Vx*Bxi+Wl.Vy*Ul.By+Wl.Vz*Ul.Bz - vbstl))/sdml;
  pbstl = Cons1D_to_Prim1D(&Ulst,&Wlst,&Bxi);


/* Ur* */
  Urst.Mx = Urst.d * spd[2];
  if(spd[2] == Wr.Vx) {
    Urst.My = Urst.d * Wr.Vy;
    Urst.Mz = Urst.d * Wr.Vz;
  }
  else {
    tmp = Bxi*(sdr-sdmr)/(Ur.d*sdr*sdmr-Bxsq);
    Urst.My = Urst.d * (Wr.Vy - Ur.By*tmp);
    Urst.Mz = Urst.d * (Wr.Vz - Ur.Bz*tmp);
  }
  if(Ur.By == 0.0 && Ur.Bz == 0.0) {
    Urst.By = 0.0;
    Urst.Bz = 0.0;
  }
  else {
    tmp = (Ur.d*SQR(sdr)-Bxsq)/(Ur.d*sdr*sdmr - Bxsq);
    Urst.By = Ur.By * tmp;
    Urst.Bz = Ur.Bz * tmp;
  }
  vbstr = (Urst.Mx*Bxi+Urst.My*Urst.By+Urst.Mz*Urst.Bz)/Urst.d;
  Urst.E = (sdr*Ur.E - ptr*Wr.Vx + ptst*spd[2] +
            Bxi*(Wr.Vx*Bxi+Wr.Vy*Ur.By+Wr.Vz*Ur.Bz - vbstr))/sdmr;
  pbstr = Cons1D_to_Prim1D(&Urst,&Wrst,&Bxi);


/* Ul** and Ur** - if Bx is zero, same as *-states */
  if(Bxi == 0.0) {
    Uldst = Ulst;
    Urdst = Urst;
    Wldst = Wlst;
    Wrdst = Wrst;
  }
  else {
    invsumd = 1.0/(sqrtdl + sqrtdr);
    if(Bxi > 0) Bxsig =  1;
    else        Bxsig = -1;

    Uldst.d = Ulst.d;
    Urdst.d = Urst.d;

    Uldst.Mx = Ulst.Mx;
    Urdst.Mx = Urst.Mx;

    tmp = invsumd*(sqrtdl*Wlst.Vy + sqrtdr*Wrst.Vy + Bxsig*(Urst.By-Ulst.By));
    Uldst.My = Uldst.d * tmp;
    Urdst.My = Urdst.d * tmp;

    tmp = invsumd*(sqrtdl*Wlst.Vz + sqrtdr*Wrst.Vz + Bxsig*(Urst.Bz-Ulst.Bz));
    Uldst.Mz = Uldst.d * tmp;
    Urdst.Mz = Urdst.d * tmp;

    tmp = invsumd*(sqrtdl*Urst.By + sqrtdr*Ulst.By +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vy-Wlst.Vy));
    Uldst.By = Urdst.By = tmp;

    tmp = invsumd*(sqrtdl*Urst.Bz + sqrtdr*Ulst.Bz +
                   Bxsig*sqrtdl*sqrtdr*(Wrst.Vz-Wlst.Vz));
    Uldst.Bz = Urdst.Bz = tmp;

    tmp = spd[2]*Bxi + (Uldst.My*Uldst.By + Uldst.Mz*Uldst.Bz)/Uldst.d;
    Uldst.E = Ulst.E - sqrtdl*Bxsig*(vbstl - tmp);
    Urdst.E = Urst.E + sqrtdr*Bxsig*(vbstr - tmp);
  }
  pbdstl = Cons1D_to_Prim1D(&Uldst,&Wldst,&Bxi);
  pbdstr = Cons1D_to_Prim1D(&Urdst,&Wrdst,&Bxi);

/*--- Step 5. ------------------------------------------------------------------
 * Compute flux
 */

  Fl.d  = Ul.Mx;
  Fl.Mx = Ul.Mx*Wl.Vx + ptl - Bxsq;
  Fl.My = Ul.d*Wl.Vx*Wl.Vy - Bxi*Ul.By;
  Fl.Mz = Ul.d*Wl.Vx*Wl.Vz - Bxi*Ul.Bz;
  Fl.E  = Wl.Vx*(Ul.E + ptl - Bxsq) - Bxi*(Wl.Vy*Ul.By + Wl.Vz*Ul.Bz);
  Fl.By = Ul.By*Wl.Vx - Bxi*Wl.Vy;
  Fl.Bz = Ul.Bz*Wl.Vx - Bxi*Wl.Vz;

  Fr.d  = Ur.Mx;
  Fr.Mx = Ur.Mx*Wr.Vx + ptr - Bxsq;
  Fr.My = Ur.d*Wr.Vx*Wr.Vy - Bxi*Ur.By;
  Fr.Mz = Ur.d*Wr.Vx*Wr.Vz - Bxi*Ur.Bz;
  Fr.E  = Wr.Vx*(Ur.E + ptr - Bxsq) - Bxi*(Wr.Vy*Ur.By + Wr.Vz*Ur.Bz);
  Fr.By = Ur.By*Wr.Vx - Bxi*Wr.Vy;
  Fr.Bz = Ur.Bz*Wr.Vx - Bxi*Wr.Vz;

  if(spd[1] >= 0) {
/* return Fl* */
    pFlux->d  = Fl.d  + spd[0]*(Ulst.d  - Ul.d);
    pFlux->Mx = Fl.Mx + spd[0]*(Ulst.Mx - Ul.Mx);
    pFlux->My = Fl.My + spd[0]*(Ulst.My - Ul.My);
    pFlux->Mz = Fl.Mz + spd[0]*(Ulst.Mz - Ul.Mz);
    pFlux->E  = Fl.E  + spd[0]*(Ulst.E  - Ul.E);
    pFlux->By = Fl.By + spd[0]*(Ulst.By - Ul.By);
    pFlux->Bz = Fl.Bz + spd[0]*(Ulst.Bz - Ul.Bz);
  }
  else if(spd[2] >= 0) {
/* return Fl** */
    tmp = spd[1] - spd[0];
    pFlux->d  = Fl.d  - spd[0]*Ul.d  - tmp*Ulst.d  + spd[1]*Uldst.d;
    pFlux->Mx = Fl.Mx - spd[0]*Ul.Mx - tmp*Ulst.Mx + spd[1]*Uldst.Mx;
    pFlux->My = Fl.My - spd[0]*Ul.My - tmp*Ulst.My + spd[1]*Uldst.My;
    pFlux->Mz = Fl.Mz - spd[0]*Ul.Mz - tmp*Ulst.Mz + spd[1]*Uldst.Mz;
    pFlux->E  = Fl.E  - spd[0]*Ul.E  - tmp*Ulst.E  + spd[1]*Uldst.E;
    pFlux->By = Fl.By - spd[0]*Ul.By - tmp*Ulst.By + spd[1]*Uldst.By;
    pFlux->Bz = Fl.Bz - spd[0]*Ul.Bz - tmp*Ulst.Bz + spd[1]*Uldst.Bz;
  }
  else if(spd[3] > 0) {
/* return Fr** */
    tmp = spd[3] - spd[4];
    pFlux->d  = Fr.d  - spd[4]*Ur.d  - tmp*Urst.d  + spd[3]*Urdst.d;
    pFlux->Mx = Fr.Mx - spd[4]*Ur.Mx - tmp*Urst.Mx + spd[3]*Urdst.Mx;
    pFlux->My = Fr.My - spd[4]*Ur.My - tmp*Urst.My + spd[3]*Urdst.My;
    pFlux->Mz = Fr.Mz - spd[4]*Ur.Mz - tmp*Urst.Mz + spd[3]*Urdst.Mz;
    pFlux->E  = Fr.E  - spd[4]*Ur.E  - tmp*Urst.E  + spd[3]*Urdst.E;
    pFlux->By = Fr.By - spd[4]*Ur.By - tmp*Urst.By + spd[3]*Urdst.By;
    pFlux->Bz = Fr.Bz - spd[4]*Ur.Bz - tmp*Urst.Bz + spd[3]*Urdst.Bz;
  }
  else {
/* return Fr* */
    pFlux->d  = Fr.d  + spd[4]*(Urst.d  - Ur.d);
    pFlux->Mx = Fr.Mx + spd[4]*(Urst.Mx - Ur.Mx);
    pFlux->My = Fr.My + spd[4]*(Urst.My - Ur.My);
    pFlux->Mz = Fr.Mz + spd[4]*(Urst.Mz - Ur.Mz);
    pFlux->E  = Fr.E  + spd[4]*(Urst.E  - Ur.E);
    pFlux->By = Fr.By + spd[4]*(Urst.By - Ur.By);
    pFlux->Bz = Fr.Bz + spd[4]*(Urst.Bz - Ur.Bz);
  }

  return;
}
