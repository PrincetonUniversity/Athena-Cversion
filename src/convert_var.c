#include "copyright.h"
/*==============================================================================
 * FILE: convert_var.c
 *
 * PURPOSE: Functions to convert conservative->primitive variables, and
 *   vice-versa. Routines for both Newtonian and special relativistic physics
 *   are included. Also contains function to compute fast magnetosonic speed for
 *   Newtonian flows only.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   Cons_to_Prim()     - converts Cons type to Prim type
 *   Cons1D_to_Prim1D() - converts 1D vector (Bx passed through arguments)
 *   Prim1D_to_Cons1D() - converts 1D vector (Bx passed through arguments)
 *   cfast()            - computes fast magnetosonic speed
 * For special relativity, there are two versions of the Cons1D_to_Prim1D
 * functions, one for HYDRO and one for MHD.
 *
 * REFERENCE:
 *   S. Noble et al., "Primitive Variable Solvers for Conservative General
 *   Relativistic Magnetohydrodynamics", ApJ. 641, 626 (2006)
 *
 *   A. Mignone & J. McKinney, "Equation of state in relativistic MHD: variable
 *   versus constant adiabatic index", MNRAS, 378, 1118 (2007)
 *
 * HISTORY: (functions for special relativity)
 *   First version written by Kevin Tian, summer 2007.
 *   Rewritten and revised by Jonathan Fulton, Senior Thesis 2009.
 *   Rewritten to use M&McK (2007) by Kris Beckwith, Fall 2009.
 *============================================================================*/

#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#if defined(SPECIAL_RELATIVITY) && defined(MHD)
/* prototypes for private functions needed with SR MHD */
Real Gamma_1overGamma;

static Real calc_func (Real Q, Real E, Real Bsq, Real Ssq, Real Vsq, Real pgas);
static Real calc_dfunc (Real Q, Real Bsq, Real Msq, Real Ssq, Real d, Real Vsq,
  Real Gsq, Real Chi);
static Real calc_vsq (Real Bsq, Real Msq, Real Ssq, Real Q);
static Real calc_chi (Real d, Real Vsq, Real Gsq, Real Q);
#endif /* SPECIAL_RELATIVITY && MHD */


/*----------------------------------------------------------------------------*/
/* Cons_to_Prim: wrapper for the Cons1D_to_Prim1D function, works for both
 *   NEWTONIAN and SPECIAL_RELATIVITY
 * conserved variables = (d,M1,M2,M3,[E],[B1c,B2c,B3c],[s(n)])
 * primitive variables = (d,V1,V2,V3,[P],[B1c,B2c,B3c],[r(n)])
 */

PrimS Cons_to_Prim(const ConsS *pCons)
{
  Cons1DS U;
  Prim1DS W;
  PrimS Prim;
#if (NSCALARS > 0)
  int n;
#endif
  Real Bx=0.0;

  U.d  = pCons->d;
  U.Mx = pCons->M1;
  U.My = pCons->M2;
  U.Mz = pCons->M3;
#ifndef ISOTHERMAL
  U.E  = pCons->E;
#endif /* ISOTHERMAL */
#ifdef MHD
  Bx = pCons->B1c;
  U.By = pCons->B2c;
  U.Bz = pCons->B3c;
#endif /* MHD */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) U.s[n] = pCons->s[n];
#endif

  W = Cons1D_to_Prim1D(&U, &Bx);

  Prim.d  = W.d;
  Prim.V1 = W.Vx;
  Prim.V2 = W.Vy;
  Prim.V3 = W.Vz;
#ifndef ISOTHERMAL
  Prim.P = W.P;
#endif /* ISOTHERMAL */
#ifdef MHD
  Prim.B1c = Bx;
  Prim.B2c = W.By;
  Prim.B3c = W.Bz;
#endif /* MHD */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) Prim.r[n] = W.r[n];
#endif

  return Prim;
}

#ifndef SPECIAL_RELATIVITY /* Following versions for Newtonian dynamics */
/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: NEWTONIAN VERSION
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz],[s(n)])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz],[r(n)])
 * Bx is passed in through the argument list.
 */

Prim1DS Cons1D_to_Prim1D(const Cons1DS *pU, const Real *pBx)
{
  Prim1DS Prim1D;
#if (NSCALARS > 0)
  int n;
#endif
  Real di = 1.0/pU->d;

  Prim1D.d  = pU->d;
  Prim1D.Vx = pU->Mx*di;
  Prim1D.Vy = pU->My*di;
  Prim1D.Vz = pU->Mz*di;

#ifndef ISOTHERMAL
  Prim1D.P = pU->E - 0.5*(SQR(pU->Mx)+SQR(pU->My)+SQR(pU->Mz))*di;
#ifdef MHD
  Prim1D.P -= 0.5*(SQR(*pBx) + SQR(pU->By) + SQR(pU->Bz));
#endif /* MHD */
  Prim1D.P *= Gamma_1;
  Prim1D.P = MAX(Prim1D.P,TINY_NUMBER);
#endif /* ISOTHERMAL */

#ifdef MHD
  Prim1D.By = pU->By;
  Prim1D.Bz = pU->Bz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) Prim1D.r[n] = pU->s[n]*di;
#endif

  return Prim1D;
}

/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D: NEWTONIAN VERSION
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz],[r(n)])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz],[s(n)])
 * Bx is passed in through the argument list.
 */

Cons1DS Prim1D_to_Cons1D(const Prim1DS *pW, const Real *pBx)
{
  Cons1DS Cons1D;
#if (NSCALARS > 0)
  int n;
#endif

  Cons1D.d  = pW->d;
  Cons1D.Mx = pW->d*pW->Vx;
  Cons1D.My = pW->d*pW->Vy;
  Cons1D.Mz = pW->d*pW->Vz;

#ifndef ISOTHERMAL
  Cons1D.E = pW->P/Gamma_1 + 0.5*pW->d*(SQR(pW->Vx) +SQR(pW->Vy) +SQR(pW->Vz));
#ifdef MHD
  Cons1D.E += 0.5*(SQR(*pBx) + SQR(pW->By) + SQR(pW->Bz));
#endif /* MHD */
#endif /* ISOTHERMAL */

#ifdef MHD
  Cons1D.By = pW->By;
  Cons1D.Bz = pW->Bz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) Cons1D.s[n] = pW->r[n]*pW->d;
#endif

  return Cons1D;
}

/*----------------------------------------------------------------------------*/
/* cfast: returns fast magnetosonic speed given input 1D vector of conserved
 *   variables and Bx.   NEWTONIAN PHYSICS ONLY.
 */

Real cfast(const Cons1DS *U, const Real *Bx)
{
  Real asq;
#ifndef ISOTHERMAL
  Real p,pb=0.0;
#endif
#ifdef MHD
  Real ctsq,casq,tmp,cfsq;
#endif

#ifdef ISOTHERMAL
  asq = Iso_csound2;
#else
#ifdef MHD
  pb = 0.5*(SQR(*Bx) + SQR(U->By) + SQR(U->Bz));
#endif /* MHD */
  p = Gamma_1*(U->E - pb - 0.5*(SQR(U->Mx)+SQR(U->My)+SQR(U->Mz))/U->d);
  asq = Gamma*p/U->d;
#endif /* ISOTHERMAL */

#ifndef MHD
  return sqrt(asq);
#else
  ctsq = (SQR(U->By) + SQR(U->Bz))/U->d;
  casq = SQR(*Bx)/U->d;
  tmp = casq + ctsq - asq;
  cfsq = 0.5*((asq+ctsq+casq) + sqrt(tmp*tmp + 4.0*asq*ctsq));
  return sqrt(cfsq);
#endif
}
#endif /* not SPECIAL_RELATIVITY */

#if defined(SPECIAL_RELATIVITY) && defined(HYDRO) /* special relativity only */
/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: SPECIAL RELATIVISTIC HYDRODYNAMICS VERSION
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 * Bx is passed in through the argument list.
 */

Prim1DS Cons1D_to_Prim1D(const Cons1DS *U, const Real *Bx)
{
  Prim1DS Prim1D;
  Real Msq, M, ME, Dsq, Gamma_1sq, denom;
  Real a3, a2, a1, a0;
  Real i1, i2, i3;
  Real iR, iS, iT;
  Real ix1=0.0, iB, iC, v, vOverM;

/* Step One: Reduce equations to a quartic polynomial in v */

  Msq = SQR(U->Mx) + SQR(U->My) + SQR(U->Mz);
  M = sqrt(Msq);

  if (fabs(M) < TINY_NUMBER) {

    v = 0.0;   /* if M is zero, then v must be identically zero */

  } else {

    ME = M * U->E;
    Dsq = SQR(U->d);

    Gamma_1sq = SQR(Gamma_1);

    denom = 1.0 / (Gamma_1sq * (Msq + Dsq));

    a3 = (-2.0 * Gamma * Gamma_1 * ME) * denom;
    a2 = (SQR(Gamma) * SQR(U->E) + 2.0*Gamma_1*Msq - Gamma_1sq*Dsq)*denom;
    a1 = (-2.0 * Gamma * ME) * denom;
    a0 = Msq * denom;

/* Step Two: Solve the polynomial analytically */

    i1 = -a2;
    i2 = a3 * a1 - 4.0 * a0;
    i3 = 4.0 * a2 * a0 - SQR(a1) - SQR(a3) * a0;

    iR = (9.0 * i1 * i2 - 27.0 * i3 - 2.0 * SQR(i1) * i1) / 54.0;
    iS = (3.0 * i2 - SQR(a2)) / 9.0;
    iT = SQR(iR) + SQR(iS) * iS;

/* iT may be negative, but ix1 then adds a complex number and its conjugate,
 * giving us a real value, so ix1 is real no matter what */

    if (iT < 0) {
      ix1 = 2.0*pow(sqrt(iR*iR + iT),(ONE_3RD))*cos(atan2(sqrt(-iT),iR)/3.0)
             - i1/3.0;
    } else {
      ix1 = pow((iR + sqrt(iT)),(ONE_3RD)) + pow((iR - sqrt(iT)),(ONE_3RD)) 
             - i1/3.0;
    }

    iB = 0.5*(a3 + sqrt(SQR(a3) - 4.0*a2 + 4.0*ix1));
    iC = 0.5*(ix1 - sqrt(SQR(ix1) - 4.0*a0));
    v = (-iB + sqrt(SQR(iB) - 4.0*iC))/2.0;

  }

/* Step Three: Resolve primitives with the value of v */

  /* ensure that v is physical */
  v = MAX(v, 0.0);
  v = MIN(v, 1.0 - 1.0e-15);

  vOverM = v / M;

  if (fabs(M) < TINY_NUMBER) {
    vOverM = 0.0;
  }

  Prim1D.d = sqrt(1.0 - SQR(v)) * U->d;

  Prim1D.Vx = U->Mx * vOverM;
  Prim1D.Vy = U->My * vOverM;
  Prim1D.Vz = U->Mz * vOverM;

  Prim1D.P = Gamma_1*
    ((U->E - U->Mx*Prim1D.Vx - U->My*Prim1D.Vy - U->Mz*Prim1D.Vz) - Prim1D.d);

  return Prim1D;
}
#endif /* SPECIAL_RELATIVITY && HYDRO */

#if defined(SPECIAL_RELATIVITY) && defined(MHD) /* special relativity only */
/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: SPECIAL RELATIVISTIC MHD VERSION 
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 * Bx is passed in through the argument list.
 *
 * IMPORTANT: This algorithm uses an iterative (Newton-Raphson) root-finding
 * step, which requires an initial guess for W. This is provided by solving
 * a cubic equation for W based on passed values of U->E & U->d and assuming
 * that v^2 = 1, as in Appendix A3 of Mignone & McKinney. Note that the
 * conserved quantity is the total energy, 
 * E = D*h*\gamma - p + 0.5*B^2 + 0.5*(v^2*B^2 - v \dot B)
 */

Prim1DS Cons1D_to_Prim1D(const Cons1DS *U, const Real *Bx)
{
  Prim1DS Prim1D;
  Real Bsq = 0.0, Msq = 0.0, S = 0.0, Ssq = 0.0;
  Real Qp = 0.0, Q = 0.0, Ep = 0.0, E = 0.0, d = 0.0;
  Real scrh1, scrh2, tmp1, tmp2;
  Real Usq, Vsq, Gsq, rho, Chi, pgas, ent;
  Real fQ, dfQ, dQstep;
  Real tol = 1.0e-3;

  int nr_success;
  int q_incs;

  Gamma_1overGamma = Gamma_1/Gamma;

  /* Calculate Bsq = B^2, Msq = M^2, S = M \dot B, Ssq = S^2 */
  Bsq = SQR((*Bx)) + SQR(U->By) + SQR(U->Bz);
  Msq = SQR(U->Mx) + SQR(U->My) + SQR(U->Mz);
  S = U->Mx * (*Bx) + U->My * U->By + U->Mz * U->Bz;
  Ssq = SQR(S);

  /* Assign input energy & density to local variables */
  E = U->E;
  d = U->d;
/*printf ("Start pars %10.4e %10.4e %10.4e %10.4e %10.4e\n",Bsq,Msq,Ssq,E,d);*/

  /* Starting guess for W, based on taking the +ve */
  /* root of Eqn. A27, guarantees that p is +ve    */
  scrh1 = -4.0*(E - Bsq);
  scrh2 = Msq - 2.0*E*Bsq + Bsq*Bsq;
  Q = ( - scrh1 + sqrt(fabs(scrh1*scrh1 - 12.0*scrh2)))/6.0;
  /*printf ("Initial guess %10.4e\n",Q);*/

  if (Q < 0.0) {
    printf ("Startup parameters %10.4e %10.4e %10.4e %10.4e %10.4e\n",
      Bsq,Msq,Ssq,E,d);
    printf("[Cons1D_to_Prim1D]: Initial guess has Q < 0 %10.4e\n",Q);
    Q = d;
  } else if (Q != Q) {
    printf ("Startup parameters %10.4e %10.4e %10.4e %10.4e %10.4e\n",
      Bsq,Msq,Ssq,E,d);
    ath_error("[Cons1D_to_Prim1D]: Initial guess has Q = NaN %10.4e\n",Q);
  }

  /* 1d NR algorithm to find W from A1, converges */
  /* when W'(k+1) / W(k) < tol, or max iterations reached */
  dQstep = 1.0;
  nr_success = 0;
  q_incs = 0;
  while (nr_success == 0 && q_incs < 100) {

    if (fabs(dQstep) <= tol) nr_success = 1;

    /* Obtain scalar quantities based on current solution estimate */
    Vsq = calc_vsq (Bsq,Msq,Ssq,Q);
    if (Vsq != Vsq){
      printf("Vsq variables %10.4e %10.4e %10.4e %10.4e\n",Bsq,Msq,Ssq,Q);
    }
    Gsq = 1.0/(1.0-Vsq);
    Chi = calc_chi (d,Vsq,Gsq,Q);
    rho = d / sqrt(fabs(Gsq));
    pgas = Gamma_1*Chi/Gamma;

    /* Evaluate Eqn. A24, A25 for the total energy density */ 
    fQ = calc_func (Q, E, Bsq, Ssq, Vsq, pgas);
    dfQ = calc_dfunc (Q, Bsq, Msq, Ssq, d, Vsq, Gsq, Chi);

    /* Check that we didn't get a NaN in the process */
    /*printf ("Function evaluations %12.6e %12.6e\n",fQ,dfQ);*/
    if (fQ != fQ) {
      printf("Soln %10.4e %10.4e %10.4e %10.4e %10.4e\n",rho,pgas,Chi,Vsq,Gsq);
      ath_error("[Cons1D_to_Prim1D]: Got a NaN in fQ\n");
    }
    if (dfQ != dfQ) {
      printf("Soln %10.4e %10.4e %10.4e %10.4e %10.4e\n",rho,pgas,Chi,Vsq,Gsq);
      ath_error("[Cons1D_to_Prim1D]: Got a NaN in dfQ\n");
    }

/* If we start out very close to the solution try and make sure we don't
 * overshoot.--- Not actually clear that this does anything, so it can probably
 * be removed */
    if (fabs(fQ) < 0.1 && q_incs == 0) {
      Q *= 10;
      Vsq = calc_vsq (Bsq,Msq,Ssq,Q);
      Gsq = 1.0/(1.0-Vsq);
      Chi = calc_chi (d,Vsq,Gsq,Q);
      rho = d / sqrt(fabs(Gsq));
      pgas = Gamma_1*Chi/Gamma;

      fQ = calc_func (Q, E, Bsq, Ssq, Vsq, pgas);
      dfQ = calc_dfunc (Q, Bsq, Msq, Ssq, d, Vsq, Gsq, Chi);
    }
    
/* Calculate step for NR update, check that it's not a Nan, update/check Q */
    dQstep = fQ / dfQ;
    if (dQstep != dQstep){
      printf("NR variables %10.4e %10.4e %10.4e\n",fQ,dfQ,dQstep);
      ath_error("[Cons1D_to_Prim1D]: Got a NaN in dQstep\n");
    }
    Q -= dQstep;
    if (Q != Q){
      printf("NR variables %10.4e %10.4e %10.4e\n",fQ,dfQ,dQstep);
      ath_error("[Cons1D_to_Prim1D]: Got a NaN in Q\n");
    }
    q_incs ++;

/* Rinse and repeat */

/*printf ("Cycle parameters %10.4e %10.4e %10.4e %1.1i %1.1i \n",
 *Q,fQ,dQstep,nr_success,q_incs);*/
  }

/* If we convered (indicated by nr_success = 1) then check solution */
  if (nr_success == 1){
    /*Q = MAX(Q, d*(1.0e0+1.0e-6));*/

    Vsq = calc_vsq (Bsq,Msq,Ssq,Q);
    Gsq = 1.0/(1.0-Vsq);
    Chi = calc_chi (d,Vsq,Gsq,Q);
    rho = d / sqrt(fabs(Gsq));
    pgas = Gamma_1*Chi/Gamma;
/*printf("Soln %10.4e %10.4e %10.4e %10.4e %10.4e\n",rho,pgas,Chi,Vsq,Gsq);*/

    if (pgas < 0.0) nr_success = 2;
    if (Vsq > 1.0) nr_success =  3;
    if (Vsq < 0.0) nr_success =  4;
  }

  if (nr_success == 2) {
/* Case of p < 0, which triggers an exit currently */
/* Should be amended to either solve the cold MHD eqn's, fall back on the
 * entropy or tell the integrator to switch to a more diffusive update */
/*    printf("Soln %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",
      rho,pgas,Chi,Vsq,Gsq,Q);
    ath_error("[Cons1D_to_Prim1D]: Got a negative pressure\n");
*/
    tmp1 = 1.0 / Q;
    tmp2 = 1.0 / (Q + Bsq);
    Prim1D.d = MAX(rho,TINY_NUMBER);
    Prim1D.P = MAX(pgas,1.0e-5);
    Prim1D.Vx = (U->Mx + S*(*Bx)*tmp1)*tmp2;
    Prim1D.Vy = (U->My + S*U->By*tmp1)*tmp2;
    Prim1D.Vz = (U->Mz + S*U->Bz*tmp1)*tmp2;

    Prim1D.By = U->By;
    Prim1D.Bz = U->Bz;

  } else if (nr_success == 3) {
/* Case of V^2 > 0, in which case we try again with some rescalings */
    printf("[Cons1D_to_Prim1D]: Got a superluminal velocity, fixing\n");

    tmp1 = 1.0 / Q;
    tmp2 = 1.0 / (Q + Bsq);
    Prim1D.Vx = (U->Mx + S*(*Bx)*tmp1)*tmp2;
    Prim1D.Vy = (U->My + S*U->By*tmp1)*tmp2;
    Prim1D.Vz = (U->Mz + S*U->Bz*tmp1)*tmp2;

    scrh1 = SQR(Prim1D.Vx) + SQR(Prim1D.Vy) + SQR(Prim1D.Vz);

    Prim1D.Vx *= 0.9999/scrh1;
    Prim1D.Vy *= 0.9999/scrh1;
    Prim1D.Vz *= 0.9999/scrh1;
    Vsq = SQR(Prim1D.Vx) + SQR(Prim1D.Vy) + SQR(Prim1D.Vz);

    Gsq = 1.0/(1.0-Vsq);
    Chi = calc_chi (d,Vsq,Gsq,Q);
    rho = d / sqrt(fabs(Gsq));
    pgas = Gamma_1*Chi/Gamma;
    printf("Soln %10.4e %10.4e %10.4e %10.4e %10.4e\n",rho,pgas,Chi,Vsq,Gsq);

    Prim1D.d = MAX(rho,TINY_NUMBER);
    Prim1D.P = MAX(pgas,TINY_NUMBER);
    Prim1D.By = U->By;
    Prim1D.Bz = U->Bz;

  } else if (nr_success == 4) {
/* Case of v^2 < 0, causes an exit */
    ath_error("[Cons1D_to_Prim1D]: Got an imaginary velocity\n");
  } else if (nr_success == 1) {
/* It worked!!! Should have a valid solution, so now set up primitives */
    tmp1 = 1.0 / Q;
    tmp2 = 1.0 / (Q + Bsq);
    Prim1D.d = MAX(rho,TINY_NUMBER);
    Prim1D.P = MAX(pgas,TINY_NUMBER);
    Prim1D.Vx = (U->Mx + S*(*Bx)*tmp1)*tmp2;
    Prim1D.Vy = (U->My + S*U->By*tmp1)*tmp2;
    Prim1D.Vz = (U->Mz + S*U->Bz*tmp1)*tmp2;

    Prim1D.By = U->By;
    Prim1D.Bz = U->Bz;
  } else {
/* NR iteration failed to converge, causes an exit */
    printf ("Start pars %10.4e %10.4e %10.4e %10.4e %10.4e\n",Bsq,Msq,Ssq,E,d);
    scrh1 = -4.0*(E - Bsq);
    scrh2 = Msq - 2.0*E*Bsq + Bsq*Bsq;
    Qp = ( - scrh1 + sqrt(scrh1*scrh1 - 12.0*scrh2))/6.0;
    Qp = MAX(Q, d*(1.0e0+1.0e-6));
    printf ("Initial guess %10.4e\n",Qp);
    ath_error("[Cons1D_to_Prim1D]: NR iteration failed to converge %10.4e %10.4e %6.6i %6.6i \n",Q,dQstep,nr_success,q_incs); 
  }

  return Prim1D;
}
#endif /* SPECIAL_RELATIVITY && MHD */

#ifdef SPECIAL_RELATIVITY /* special relativity only */
/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D: SPECIAL RELATIVITY VERSION
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 * Bx is passed in through the argument list.
 */

Cons1DS Prim1D_to_Cons1D(const Prim1DS *W, const Real *Bx)
{
  Cons1DS Cons1D;
  Real vsq, U0, Bsq = 0.0, vDotB = 0.0, wU0sq;

  vsq = SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz);
  U0 = 1.0 / (1.0 - vsq);

#ifdef MHD
  Bsq = (*Bx)*(*Bx) + SQR(W->By) + SQR(W->Bz);
  vDotB = (*Bx)*W->Vx + W->By*W->Vy + W->Bz*W->Vz;
#endif

  wU0sq = (W->d + Gamma/Gamma_1 * W->P)*U0;

  Cons1D.d  = sqrt(U0) * W->d;
  Cons1D.Mx = wU0sq * W->Vx;
  Cons1D.My = wU0sq * W->Vy;
  Cons1D.Mz = wU0sq * W->Vz;

  Cons1D.E  = wU0sq - W->P;

#ifdef MHD
  Cons1D.Mx += Bsq*W->Vx - vDotB*(*Bx);
  Cons1D.My += Bsq*W->Vy - vDotB*W->By;
  Cons1D.Mz += Bsq*W->Vz - vDotB*W->Bz;

  Cons1D.E  += (1.0 + vsq) * Bsq/2.0 - SQR(vDotB) / 2.0;

  Cons1D.By = W->By;
  Cons1D.Bz = W->Bz;
#endif /* MHD */

  return Cons1D;
}
#endif /* SPECIAL_RELATIVITY */

#if defined(SPECIAL_RELATIVITY) && defined(MHD)
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
static Real calc_func (Real Q, Real E, Real Bsq, Real Ssq, Real Vsq, Real pgas)
{
  Real tmp1;
  /* Evaluate equation A25 for E, rather than E' */
  return Q - pgas + 0.5*(1.0+Vsq)*Bsq - (0.5*Ssq/Q/Q) - E;
}

static Real calc_dfunc(Real Q, Real Bsq, Real Msq, Real Ssq, Real d, Real Vsq,
       Real Gsq, Real Chi)
{
  /* Evaluate A8 for E & Q, rahter than E', W' */
  Real dp_dQ, dp_dchi, dchi_dQ, dp_drho, drho_dQ, dVsq_dQ;
  Real G, Qsq, Qth, scrh1;

  Qsq = Q*Q;
  Qth = Qsq*Q;
  G = sqrt(fabs(Gsq));

  scrh1 = Q + Bsq;
  dVsq_dQ  = Ssq*(3.0*Q*scrh1 + Bsq*Bsq) + Msq*Qth;
  dVsq_dQ *= -2.0/Qth/(scrh1*scrh1*scrh1);

  /* -- kinematical terms -- */
  /* see eqn. A14 - A16 */
  dchi_dQ =  1.0 - Vsq - 0.5*G*(d + 2.0*Chi*G)*dVsq_dQ;
  drho_dQ = -0.5*d*G*dVsq_dQ;

  /* -- thermo terms, change for different EOS --*/
  /* see Section A2 of M&M */
  dp_dchi = (Gamma - 1.0)/Gamma;
  dp_drho = 0.0;

  dp_dQ = dp_dchi*dchi_dQ + dp_drho*drho_dQ;

  return 1.0 - dp_dQ +0.5*Bsq*dVsq_dQ + Ssq / Qth;
}

static Real calc_vsq (Real Bsq, Real Msq, Real Ssq, Real Q)
{
  /* Obtain v^2 via eqn. A3 */
  Real Qsq, Ssq_Qsq, scrh1;
  Qsq = Q*Q;
  Ssq_Qsq = Ssq / Qsq;
  scrh1 = Q + Bsq;
  return (Msq + Ssq_Qsq*(scrh1 + Q))/(scrh1*scrh1);
}

static Real calc_chi (Real d, Real Vsq, Real Gsq, Real Q)
{
  /* Evaluate \Chi from eqn. A11 using Q, rathher than Q' */
  Real G, tmp1, tmp2;
  G = sqrt(fabs(Gsq));
  tmp1 = 1.0 - Vsq;
  tmp2 = Q - d*G;
  return tmp2*tmp1;
}
#endif /* SPECIAL_RELATIVITY && MHD */
