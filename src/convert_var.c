#include "copyright.h"
/*==============================================================================
 * FILE: convert_var.c
 *
 * PURPOSE: Functions to convert conservative to primitive vars, and vice versa.
 *   Also contains function to compute fast magnetosonic speed.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   Cons1D_to_Prim1D() - converts 1D vector (Bx passed through arguments)
 *   Prim1D_to_Cons1D() - converts 1D vector (Bx passed through arguments)
 *   cfast()            - compute fast speed given input Cons1D, Bx
 *============================================================================*/

#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: 
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz],[s(n)])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz],[r(n)])
 * Bx is passed in through the argument list.
 */

void Cons1D_to_Prim1D(const Cons1D *pU, Prim1D *pW, const Real *pBx)
{
#if (NSCALARS > 0)
  int n;
#endif
  Real di = 1.0/pU->d;

  pW->d  = pU->d;
  pW->Vx = pU->Mx*di;
  pW->Vy = pU->My*di;
  pW->Vz = pU->Mz*di;

#ifndef ISOTHERMAL
  pW->P = pU->E - 0.5*(SQR(pU->Mx)+SQR(pU->My)+SQR(pU->Mz))*di;
#ifdef MHD
  pW->P -= 0.5*(SQR(*pBx) + SQR(pU->By) + SQR(pU->Bz));
#endif /* MHD */
  pW->P *= Gamma_1;
  pW->P = MAX(pW->P,TINY_NUMBER);
#endif /* ISOTHERMAL */

#ifdef MHD
  pW->By = pU->By;
  pW->Bz = pU->Bz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) pW->r[n] = pU->s[n]*di;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D: 
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz],[r(n)])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz],[s(n)])
 * Bx is passed in through the argument list.
 */

void Prim1D_to_Cons1D(Cons1D *pU, const Prim1D *pW, const Real *pBx)
{
#if (NSCALARS > 0)
  int n;
#endif

  pU->d  = pW->d;
  pU->Mx = pW->d*pW->Vx;
  pU->My = pW->d*pW->Vy;
  pU->Mz = pW->d*pW->Vz;

#ifndef ISOTHERMAL
  pU->E = pW->P/Gamma_1 + 0.5*pW->d*(SQR(pW->Vx) + SQR(pW->Vy) + SQR(pW->Vz));
#ifdef MHD
  pU->E += 0.5*(SQR(*pBx) + SQR(pW->By) + SQR(pW->Bz));
#endif /* MHD */
#endif /* ISOTHERMAL */

#ifdef MHD
  pU->By = pW->By;
  pU->Bz = pW->Bz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) pU->s[n] = pW.r[n]*pW->d;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* cfast: returns fast magnetosonic speed given input 1D vector of conserved
 *   variables and Bx. 
 */

Real cfast(const Cons1D *U, const Real *Bx)
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
