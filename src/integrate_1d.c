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

void integrate_1d(Grid *pG, Domain *pD)
{
  Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
#if (NSCALARS > 0)
  int n;
#endif
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,dhalf;
#ifdef SELF_GRAVITY
  Real gxl,gxr,flux_m1l,flux_m1r;
#endif

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;  
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i].d  = pG->U[ks][js][i].d;
    U1d[i].Mx = pG->U[ks][js][i].M1;
    U1d[i].My = pG->U[ks][js][i].M2;
    U1d[i].Mz = pG->U[ks][js][i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = pG->U[ks][js][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    U1d[i].By = pG->U[ks][js][i].B2c;
    U1d[i].Bz = pG->U[ks][js][i].B3c;
    Bxc[i] = pG->U[ks][js][i].B1c;
    Bxi[i] = pG->B1i[ks][js][i];
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif
  }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces.
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    Cons1D_to_Prim1D(&U1d[i],&W[i],&Bxc[i]);
  }
  lr_states(W,Bxc,pG->dt,dtodx1,is,ie,Wl,Wr);

/*--- Step 1c -----------------------------------------------------------------
 * Add gravitational source terms for a static potential for dt/2 to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie+1; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phicr = (*StaticGravPot)( x1                ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
      phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* Apply gravitational source terms to velocity using gradient of potential
 * for (dt/2).   S_{V} = -Grad(Phi) */
      Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);
    }
  }

/*--- Step 1d ------------------------------------------------------------------
 * Add gravitational source terms for self-gravity for dt/2 to L/R states
 */

#ifdef SELF_GRAVITY
  for (i=is; i<=ie+1; i++) {
    Wl[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
    Wr[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
  }
#endif


/*--- Step 1e ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    Prim1D_to_Cons1D(&Ul_x1Face[i],&Wl[i],&Bxi[i]);
    Prim1D_to_Cons1D(&Ur_x1Face[i],&Wr[i],&Bxi[i]);

    GET_FLUXES(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i]);
  }

/*--- Step 2a ------------------------------------------------------------------
 * Source term update (static gravitational potential).  To improve
 * conservation of E, average the energy source term computed at cell faces. 
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1               ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      dhalf = pG->U[ks][js][i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
      pG->U[ks][js][i].M1 -= dtodx1*dhalf*(phir-phil);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                       x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Source term update (self-gravity).
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */
#ifdef SELF_GRAVITY
  for (i=is; i<=ie; i++) {
      phic = pG->Phi[ks][js][i];
      phil = 0.5*(pG->Phi[ks][js][i-1] + pG->Phi[ks][js][i  ]);
      phir = 0.5*(pG->Phi[ks][js][i  ] + pG->Phi[ks][js][i+1]);

      gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* 1-momentum flux.  2nd term is needed only if Jean's swindle used */
      flux_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flux_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

      pG->U[ks][js][i].M1 -= dtodx1*(flux_m1r - flux_m1l);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                       x1Flux[i+1].d*(phir - phic));
#endif
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */
  for (i=is; i<=ie+1; i++) {
    pG->x1MassFlux[ks][js][i] = x1Flux[i].d;
  }
#endif

/*--- Step 3 -------------------------------------------------------------------
 * Update with flux differences.
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].d  -= dtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    pG->U[ks][js][i].M1 -= dtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    pG->U[ks][js][i].M2 -= dtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    pG->U[ks][js][i].M3 -= dtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E  -= dtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    pG->U[ks][js][i].B3c -= dtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pG->U[ks][js][i].s[n] -= dtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
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
