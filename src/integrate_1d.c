#include "copyright.h"
/*==============================================================================
 * FILE: integrate_1d.c
 *
 * PURPOSE: Integrate MHD equations in 1D.  Updates U.[d,M1,M2,M3,E,B2c,B3c]
 *   in Grid structure, where U is of type Gas.  Adds gravitational source
 *   terms.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_1d()
 *   integrate_init_1d()
 *   integrate_destruct_1d()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real *Bxc=NULL, *Bxi=NULL;
static Cons1D *Ul_x1Face=NULL, *Ur_x1Face=NULL, *U1d=NULL, *x1Flux=NULL;
static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;

/*----------------------------------------------------------------------------*/
/* integrate_1d:   */

void integrate_1d(Grid *pGrid)
{
  Real dtodx1 = pGrid->dt/pGrid->dx1, dt = pGrid->dt;
  int i, is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks;
#if (NSCALARS > 0)
  int n;
#endif
  Real pb,x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;

/*--- Step 1 -------------------------------------------------------------------
 * Load 1D vector of conserved variables;  
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i].d  = pGrid->U[ks][js][i].d;
    U1d[i].Mx = pGrid->U[ks][js][i].M1;
    U1d[i].My = pGrid->U[ks][js][i].M2;
    U1d[i].Mz = pGrid->U[ks][js][i].M3;
#ifndef ISOTHERMAL
    U1d[i].E  = pGrid->U[ks][js][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
    U1d[i].By = pGrid->U[ks][js][i].B2c;
    U1d[i].Bz = pGrid->U[ks][js][i].B3c;
    Bxc[i] = pGrid->U[ks][js][i].B1c;
    Bxi[i] = pGrid->B1i[ks][js][i];
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pGrid->U[ks][js][i].s[n];
#endif
  }

/*--- Step 2 -------------------------------------------------------------------
 * Convert to primitive variables, compute L and R states at X1-interfaces.
 */

  for (i=is-2; i<=ie+2; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }
  lr_states(W,Bxc,dt,dtodx1,is,ie,Wl,Wr);

/*--- Step 3 -------------------------------------------------------------------
 * Add gravitational source terms for dt/2 from static potential to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie+1; i++) {

/* Calculate the potential at i [phi-center-right],i-1 [phi-center-left],
 * and i-1/2 [phi-face] */
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phicr = (*StaticGravPot)( x1                ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pGrid->dx1),x2,x3);
      phifc = (*StaticGravPot)((x1-0.5*pGrid->dx1),x2,x3);

/* Apply gravitational source terms to velocity using gradient of potential. */
      Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);
    }
  }

/*--- Step 4 -------------------------------------------------------------------
 * Convert back to conserved variables, and compute 1D fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    pb = Prim1D_to_Cons1D(&Ul_x1Face[i],&Wl[i],&Bxi[i]);
    pb = Prim1D_to_Cons1D(&Ur_x1Face[i],&Wr[i],&Bxi[i]);
    
  }

  for (i=is; i<=ie+1; i++) {
    GET_FLUXES(Bxi[i],Ul_x1Face[i],Ur_x1Face[i],&x1Flux[i]);
  }

/*--- Step 5 -------------------------------------------------------------------
 * To keep the gravitational source terms 2nd order, add 0.5 the gravitational
 * acceleration to the momentum equation now (using d^{n}), before the update
 * of the cell-centered variables due to flux gradient.
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phir = (*StaticGravPot)((x1+0.5*pGrid->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pGrid->dx1),x2,x3);

      pGrid->U[ks][js][i].M1 -= 0.5*dtodx1*(phir-phil)*pGrid->U[ks][js][i].d;
    }
  }

/*--- Step 6 -------------------------------------------------------------------
 * Update cell-centered variables in pGrid using 1D-fluxes
 */

  for (i=is; i<=ie; i++) {
    pGrid->U[ks][js][i].d  -= dtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    pGrid->U[ks][js][i].M1 -= dtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    pGrid->U[ks][js][i].M2 -= dtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    pGrid->U[ks][js][i].M3 -= dtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef ISOTHERMAL
    pGrid->U[ks][js][i].E  -= dtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
    pGrid->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    pGrid->U[ks][js][i].B3c -= dtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pGrid->U[ks][js][i].s[n] -= dtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
  }

/*--- Step 7 -------------------------------------------------------------------
 * Complete the gravitational source terms by adding 0.5 the gravitational
 * acceleration to the momentum equation (using d^{n+1}), and the energy source
 * terms constructed so total energy (E + \rho\Phi)  is strictly conserved
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1               ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pGrid->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pGrid->dx1),x2,x3);

      pGrid->U[ks][js][i].M1 -= 0.5*dtodx1*(phir-phil)*pGrid->U[ks][js][i].d;
#ifndef ISOTHERMAL
      pGrid->U[ks][js][i].E += dtodx1*(x1Flux[i  ].d*(phil - phic) +
                                       x1Flux[i+1].d*(phic - phir));
#endif
    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(int nx1)
{
  int Nx1;
  Nx1 = nx1 + 2*nghost;
  if ((Bxc = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  if ((U1d       = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul_x1Face = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur_x1Face = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((W  = (Prim1D*)malloc(Nx1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl = (Prim1D*)malloc(Nx1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr = (Prim1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((x1Flux    = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;

  return;

  on_error:
    ath_error("[integrate_init_1d]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_1d: Free temporary integration arrays  */

void integrate_destruct_1d(void)
{
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
  if (U1d != NULL) free(U1d);
  if (Ul_x1Face != NULL) free(Ul_x1Face);
  if (Ur_x1Face != NULL) free(Ur_x1Face);
  if (W  != NULL) free(W);
  if (Wl != NULL) free(Wl);
  if (Wr != NULL) free(Wr);
  if (x1Flux != NULL) free(x1Flux);

  return;
}
