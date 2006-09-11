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

/*----------------------------------------------------------------------------*/
/* integrate_1d:   */

void integrate_1d(Grid *pGrid)
{
  Real dtodx1 = pGrid->dt/pGrid->dx1, hdt = 0.5*pGrid->dt, dt = pGrid->dt;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks;
  int i,il = pGrid->is - 2,iu = pGrid->ie + 2;
  Real x1, x2, x3, g;

/*--- Step 1 -------------------------------------------------------------------
 * Load 1D vector of conserved variables;  U1d = (d, M1, M2, M3, E, B2c, B3c)
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
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute L and R states at X1-interfaces.
 */

  lr_states(U1d,Bxc,Bxi,dt,dtodx1,il,iu,Ul_x1Face,Ur_x1Face);

/*--- Step 3 -------------------------------------------------------------------
 * Add gravitational source terms for dt/2 from static potential to L/R states
 */

  if (x1GravAcc != NULL){
    for (i=il; i<=iu+1; i++) {

/* Calculate the face-centered acceleration */
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      g = (*x1GravAcc)((x1-0.5*pGrid->dx1),x2,x3);

/* Apply gravitational source terms to total energy and momentum.  Note E must
 * be updated first before Mx  */
#ifndef ISOTHERMAL
      Ul_x1Face[i].E += hdt*Ul_x1Face[i].Mx*g;
      Ur_x1Face[i].E += hdt*Ur_x1Face[i].Mx*g;
#endif
      Ul_x1Face[i].Mx += hdt*Ul_x1Face[i].d*g;
      Ur_x1Face[i].Mx += hdt*Ur_x1Face[i].d*g;
    }
  }

/*--- Step 4 -------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction
 */

  for (i=il; i<=iu+1; i++) {
    GET_FLUXES(Bxi[i],Ul_x1Face[i],Ur_x1Face[i],&x1Flux[i]);
  }

/*--- Step 5 -------------------------------------------------------------------
 * To keep the gravitational source terms 2nd order, add 0.5 the gravitational
 * acceleration to the momentum equation now (using d^{n}), before the update
 * of the cell-centered variables due to flux gradient.
 */

  if (x1GravAcc != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      g = (*x1GravAcc)(x1,x2,x3);

      pGrid->U[ks][js][i].M1 += hdt*pGrid->U[ks][js][i].d*g;
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
  }

/*--- Step 7 -------------------------------------------------------------------
 * Complete the gravitational source terms by adding 0.5 the acceleration at
 * time level n+1, and the energy source term at time level {n+1/2}.
 */

  if (x1GravAcc != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      g = (*x1GravAcc)(x1,x2,x3);

      pGrid->U[ks][js][i].M1 += hdt*pGrid->U[ks][js][i].d*g;
#ifndef ISOTHERMAL
      pGrid->U[ks][js][i].E += hdt*(x1Flux[i].d + x1Flux[i+1].d)*g;
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
  if (x1Flux != NULL) free(x1Flux);

  return;
}
