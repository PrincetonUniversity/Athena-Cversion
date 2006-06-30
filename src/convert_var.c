#include "copyright.h"
/*==============================================================================
 * FILE: convert_var.c
 *
 * PURPOSE: Functions to convert conservative to primitive vars, and vice versa.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   Cons1D_to_Prim1D() - converts 1D vector (lacks Bx)
 *   Prim1D_to_Cons1D() - converts 1D vector (lacks Bx)
 *   Gas_to_Prim()      - converts Gas structure (includes B1c,B2c,B3c)
 *   Prim_to_Gas()      - converts Gas structure (includes B1c,B2c,B3c)
 *============================================================================*/

#include <math.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* Cons1D_to_Prim1D: 
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 * Bx is passed in through the argument list.  Returns the magnetic pressure.
 */

Real Cons1D_to_Prim1D(const Cons1D *U, Prim1D *W, const Real *Bx)
{
  Real pb=0.0,di;
  di = 1.0/U->d;

  W->d  = U->d;
  W->Vx = U->Mx*di;
  W->Vy = U->My*di;
  W->Vz = U->Mz*di;

#ifdef MHD
  W->By = U->By;
  W->Bz = U->Bz;
  pb = 0.5*(SQR(*Bx) + SQR(U->By) + SQR(U->Bz));
#endif /* MHD */

#ifndef ISOTHERMAL
  W->P = Gamma_1*(U->E - pb - 0.5*(SQR(U->Mx)+SQR(U->My)+SQR(U->Mz))*di);
  W->P = MAX(W->P,TINY_NUMBER);
#endif /* ISOTHERMAL */

  return(pb);
}

/*----------------------------------------------------------------------------*/
/* Prim1D_to_Cons1D: 
 *   primitive variables = (d,Vx,Vy,Vz,[P],[By,Bz])
 *   conserved variables = (d,Mx,My,Mz,[E],[By,Bz])
 * Bx is passed in through the argument list.  Returns the magnetic pressure.
 */

Real Prim1D_to_Cons1D(Cons1D *U, const Prim1D *W, const Real *Bx)
{
  Real pb=0.0;

  U->d = W->d;
  U->Mx = W->Vx*W->d;
  U->My = W->Vy*W->d;
  U->Mz = W->Vz*W->d;

#ifdef MHD
  pb = 0.5*(SQR(*Bx) + SQR(W->By) + SQR(W->Bz));
  U->By = W->By;
  U->Bz = W->Bz;
#endif /* MHD */

#ifndef ISOTHERMAL
  U->E = W->P/Gamma_1 + pb + 0.5*W->d*(SQR(W->Vx) + SQR(W->Vy) + SQR(W->Vz));
#endif /* ISOTHERMAL */

  return(pb);
}

/*----------------------------------------------------------------------------*/
/* Gas_to_Prim: 
 *   conserved (Gas) variables = (d,M1,M2,M3,[E],[B1c,B2c,B3c])
 *   primitive       variables = (d,V1,V2,V3,[P],[B1c,B2c,B3c])
 * Returns the magnetic pressure.
 */


Real Gas_to_Prim(const Gas *U, Prim *W)
{
  Real pb=0.0,di;
  di = 1.0/U->d;

  W->d  = U->d;
  W->V1 = U->M1*di;
  W->V2 = U->M2*di;
  W->V3 = U->M3*di;

#ifdef MHD
  W->B1c = U->B1c;
  W->B2c = U->B2c;
  W->B3c = U->B3c;
  pb = 0.5*(SQR(U->B1c) + SQR(U->B2c) + SQR(U->B3c));
#endif /* MHD */

#ifndef ISOTHERMAL
  W->P = Gamma_1*(U->E - pb - 0.5*(SQR(U->M1)+SQR(U->M2)+SQR(U->M3))*di);
  W->P = MAX(W->P,TINY_NUMBER);
#endif /* ISOTHERMAL */

  return(pb);
}

/*----------------------------------------------------------------------------*/
/* Prim_to_Gas: 
 *   primitive       variables = (d,V1,V2,V3,[P],[B1c,B2c,B3c])
 *   conserved (Gas) variables = (d,M1,M2,M3,[E],[B1c,B2c,B3c])
 * Returns the magnetic pressure.
 */

Real Prim_to_Gas(Gas *U, const Prim *W)
{
  Real pb=0.0;

  U->d = W->d;
  U->M1 = W->V1*W->d;
  U->M2 = W->V2*W->d;
  U->M3 = W->V3*W->d;

#ifdef MHD
  pb = 0.5*(SQR(W->B1c) + SQR(W->B2c) + SQR(W->B3c));
  U->B1c = W->B1c;
  U->B2c = W->B2c;
  U->B3c = W->B3c;
#endif /* MHD */

#ifndef ISOTHERMAL
  U->E = W->P/Gamma_1 + pb + 0.5*W->d*(SQR(W->V1) + SQR(W->V2) + SQR(W->V3));
#endif /* ISOTHERMAL */

  return(pb);
}
