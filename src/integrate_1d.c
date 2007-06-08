#include "copyright.h"
/*==============================================================================
 * FILE: integrate_1d.c
 *
 * PURPOSE: Integrate MHD equations in 1D.  Updates U.[d,M1,M2,M3,E,B2c,B3c,s]
 *   in Grid structure, where U is of type Gas.  Adds gravitational source
 *   terms, and self-gravity.
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
  Real dtodx1 = pGrid->dt/pGrid->dx1, hdtodx1 = 0.5*pGrid->dt/pGrid->dx1;
  int i, is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks;
#if (NSCALARS > 0)
  int n;
#endif
  Real pb,x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,dhalf;
#ifdef SELF_GRAVITY
  Real gxl,gxr,flm1,frm1;
#endif

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

  for (i=is-nghost; i<=ie+nghost; i++) {
    pb = Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }
  lr_states(W,Bxc,pGrid->dt,dtodx1,is,ie,Wl,Wr);

/*--- Step 3a ------------------------------------------------------------------
 * Add gravitational source terms for a static potential for dt/2 to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie+1; i++) {

/* Calculate the potential at i [phi-center-right],i-1 [phi-center-left],
 * and i-1/2 [phi-face] */
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phicr = (*StaticGravPot)( x1                ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pGrid->dx1),x2,x3);
      phifc = (*StaticGravPot)((x1-0.5*pGrid->dx1),x2,x3);

/* Apply gravitational source terms to velocity using gradient of potential
 * for (dt/2).   S_{V} = -Grad(Phi) */
      Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);
    }
  }

/*--- Step 3b ------------------------------------------------------------------
 * Add gravitational source terms for self-gravity for dt/2 to L/R states
 */

#ifdef SELF_GRAVITY
  for (i=is; i<=ie+1; i++) {
    Wl[i].Vx -= hdtodx1*(pGrid->Phi[ks][js][i] - pGrid->Phi[ks][js][i-1]);
    Wr[i].Vx -= hdtodx1*(pGrid->Phi[ks][js][i] - pGrid->Phi[ks][js][i-1]);
  }
#endif


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
 * Add the gravitational source terms at second order.  To improve conservation
 * of total energy, we average the energy source term computed at cell faces. 
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1               ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pGrid->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pGrid->dx1),x2,x3);

      dhalf = pGrid->U[ks][js][i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
      pGrid->U[ks][js][i].M1 -= dtodx1*dhalf*(phir-phil);
#ifndef ISOTHERMAL
      pGrid->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                       x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 6 -------------------------------------------------------------------
 * Add divergence of gravitational stress tensor and energy source terms for
 * self-gravity, using Phi^{n}, the current value of the potential at t^{n}.
 * This is going to require a flux correction later using Phi^{n+1} to make
 * update 2nd order.
 *    dM/dt = -Div(G);  S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */
#ifdef SELF_GRAVITY
  for (i=is; i<=ie; i++) {

/*  Add the source term using Phi^{n} for a full timestep */

      phic = pGrid->Phi[ks][js][i];
      phil = 0.5*(pGrid->Phi[ks][js][i-1] + pGrid->Phi[ks][js][i  ]);
      phir = 0.5*(pGrid->Phi[ks][js][i  ] + pGrid->Phi[ks][js][i+1]);

      gxl = (pGrid->Phi[ks][js][i-1] - pGrid->Phi[ks][js][i  ])/(pGrid->dx1);
      gxr = (pGrid->Phi[ks][js][i  ] - pGrid->Phi[ks][js][i+1])/(pGrid->dx1);

/* 1-momentum flux at cell faces.  grav_mean_rho is non-zero if Jeans swindle
   was used (potential computed using (\rho - grav_mean_rho))   */
      flm1 = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      frm1 = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

      pGrid->U[ks][js][i].M1 -= dtodx1*(frm1-flm1);
#ifndef ISOTHERMAL
      pGrid->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                       x1Flux[i+1].d*(phir - phic));
#endif
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */
  for (i=is; i<=ie+1; i++) {
    pGrid->x1MassFlux[ks][js][i] = x1Flux[i].d;
  }
#endif

/*--- Step 7 -------------------------------------------------------------------
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
