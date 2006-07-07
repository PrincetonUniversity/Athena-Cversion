#include "copyright.h"
/*==============================================================================
 * FILE: integrate_1d.c
 *
 * PURPOSE: Functions to integrate MHD equations in 1D.
 * The variables in Grid which are updated are:
 *    U.[d,M1,M2,M3,E,B2c,B3c] -- where U is of type Gas
 *    time,dt,nstep
 * source terms for momentum M1 can be added through a user-specified
 * time-independent conserved potential.  In this case, energy is conserved
 * exactly, and the algorithm exploits this fact.
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
#include "prototypes.h"

static Real *Bxc=NULL, *Bxi=NULL;
static Cons1D *Ul_x1Face=NULL, *Ur_x1Face=NULL, *U1d=NULL, *x1Flux=NULL;

/*----------------------------------------------------------------------------*/
/* integrate_1d:   */

void integrate_1d(Grid *pGrid)
{
  Real dtodx1 = pGrid->dt/pGrid->dx1;
  Real hdt = 0.5*pGrid->dt, dt = pGrid->dt;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks;
  int i,il = pGrid->is - 2,iu = pGrid->ie + 2;
  Real x1, x2, x3;
  Real phic_i, phic_im1, phirx1, philx1;

/*--- Step 1 ------------------------------------------------------------------
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

/*--- Step 2 ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces.
 */

  lr_states(U1d,Bxc,Bxi,dt,dtodx1,il,iu,Ul_x1Face,Ur_x1Face);

/*--- Step 3 ------------------------------------------------------------------
 * Add source terms from time-independent conservative potential to L/R states
 */

/*
  if (cons_pot_fun != NULL){
    for (i=il; i<=iu+1; i++) {
*/

/* Calculate the cell-centered potential at i, i-1 */
/*
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phic_i   = (*cons_pot_fun)( x1           ,x2,x3);
      phic_im1 = (*cons_pot_fun)((x1-pGrid->dx),x2,x3);

      Ul_x1Face[i].M1 += hdt*Ul_x1Face[i].d*(phic_i - phic_im1)/pGrid->dx1;
      Ur_x1Face[i].M1 += hdt*Ur_x1Face[i].d*(phic_i - phic_im1)/pGrid->dx1;
*/

/* THIS MUST BE WRONG */
#ifndef ISOTHERMAL
/*
      Ul_x1Face[i].E += hdt*Ul_x1Face[i].M1*(phic_i - phic_im1);
      Ur_x1Face[i].E += hdt*Ur_x1Face[i].M1*(phic_i - phic_im1);
*/
#endif
/*
    }
  }
*/

/*--- Step 4 ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction
 */

  for (i=il; i<=iu+1; i++) {
    GET_FLUXES(Bxi[i],Ul_x1Face[i],Ur_x1Face[i],&x1Flux[i]);
  }

/*--- Step 5 ----------------------------------------------------------------
 * Update cell-centered variables using source terms from time-independent
 * conservative potential.  This step must be before the update of
 * pGrid->U[ks][js][i].d
 */

/*
  if (cons_pot_fun != NULL){
    for (i=is; i<=ie; i++) {
*/

/* Calculate the potential at cell center, left- and right-interfaces */
/*
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phic_i = (*cons_pot_fun)(x1,x2,x3);
      phirx1 = (*cons_pot_fun)(x1 + 0.5*pGrid->dx1,x2,x3);
      philx1 = (*cons_pot_fun)(x1 - 0.5*pGrid->dx1,x2,x3);

      d_nph = pGrid->U[ks][js][i].d - 0.5*dtodx1*(x1Flux[i+1].d - x1Flux[i].d);
      pGrid->U[ks][js][i].M1 += d_nph*(philx1 - phirx1)/pGrid->dx1;

*/
/* THIS MUST BE WRONG */
#ifndef ISOTHERMAL
/*
      pGrid->U[ks][js][i].E += (x1Flux[i  ].d*(philx1 - phic_i) +
                                x1Flux[i+1].d*(phic_i - phirx1))/pGrid->dx1;
*/
#endif
/*
    }
  }
*/

/*--- Step 6 ----------------------------------------------------------------
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

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(int nx1)
{
  int Nx1 = nx1 + 2*nghost;
#ifdef MHD
  if ((Bxc = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
#endif
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
#ifdef MHD
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#endif
  if (U1d != NULL) free(Uld);
  if (Ul_x1Face != NULL) free(Ul_x1Face);
  if (Ur_x1Face != NULL) free(Ur_x1Face);
  if (x1Flux != NULL) free(x1Flux);

  return;
}
