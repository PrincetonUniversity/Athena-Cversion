#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_1d_vl.c
 *
 * PURPOSE: Integrate MHD equations using 1D version of MUSCL-Hancock (VL)
 *   integrator.  Updates U.[d,M1,M2,M3,E,B2c,B3c,s] in Grid structure.
 *   Adds gravitational source terms, self-gravity.
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
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#if defined(VL_INTEGRATOR) && defined(CARTESIAN)

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1D *Wl_x1Face=NULL, *Wr_x1Face=NULL;
static Cons1D *x1Flux=NULL;

/* The interface magnetic field */
#ifdef MHD
static Real *B1_x1Face=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
#ifdef MHD
static Real *Bxc=NULL;
#endif /* MHD */
static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1D *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* conserved and primitive variables at t^{n+1/2} computed in predict step */
static Gas *Uhalf=NULL;
static Prim *Whalf=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_1d: 1D version of van Leer unsplit integrator for MHD. 
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 1D.
 */

void integrate_1d(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, hdtodx1=0.5*pG->dt/pG->dx1;
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,flx_m1l,flx_m1r;
#endif

  int il=is-(nghost-1), iu=ie+(nghost-1);

  for (i=is-nghost; i<=ie+nghost; i++) {
    Uhalf[i] = pG->U[ks][js][i];
#ifdef MHD
    B1_x1Face[i] = pG->B1i[ks][js][i];
#endif /* MHD */
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

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
    Bxc[i]    = pG->U[ks][js][i].B1c;
    U1d[i].By = pG->U[ks][js][i].B2c;
    U1d[i].Bz = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif

/* Convert to primitive variables -- for SR use primitive variables stored
 * in Grid structure */

#ifndef SPECIAL_RELATIVITY
    Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));
#else
    W[i].d  = pG->W[ks][js][i].d;
    W[i].Vx = pG->W[ks][js][i].V1;
    W[i].Vy = pG->W[ks][js][i].V2;
    W[i].Vz = pG->W[ks][js][i].V3;
#ifndef BAROTROPIC
    W[i].P  = pG->W[ks][js][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    W[i].By = pG->W[ks][js][i].B2c;
    W[i].Bz = pG->W[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W[i].r[n] = pG->W[ks][js][i].r[n];
#endif /* NSCALARS */
#endif /* SPECIAL_RELATIVITY */
  }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states in U and W */

  for (i=il; i<=ie+nghost; i++) {
    Wl[i] = W[i-1];
    Wr[i] = W[i  ];

    Ul[i] = U1d[i-1];
    Ur[i] = U1d[i  ];
  }

/*--- Step 1c ------------------------------------------------------------------
 * Compute flux in x1-direction */

  for (i=il; i<=ie+nghost; i++) {
    fluxes(Ul[i],Ur[i],Wl[i],Wr[i], MHDARG( B1_x1Face[i] , ) &x1Flux[i]);
  }

/*=== STEPS 2-4: Not needed in 1D ===*/

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables (including B2c and B3c) to half-timestep
 */

  for (i=il; i<=iu; i++) {
    Uhalf[i].d   -= hdtodx1*(x1Flux[i+1].d  - x1Flux[i].d );
    Uhalf[i].M1  -= hdtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
    Uhalf[i].M2  -= hdtodx1*(x1Flux[i+1].My - x1Flux[i].My);
    Uhalf[i].M3  -= hdtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
    Uhalf[i].E   -= hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif /* BAROTROPIC */
#ifdef MHD
    Uhalf[i].B2c -= hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    Uhalf[i].B3c -= hdtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      Uhalf[i].s[n] -= hdtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
  }

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=il; i<=iu; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      Uhalf[i].M1 -= hdtodx1*pG->U[ks][js][i].d*(phir-phil);
#ifndef BAROTROPIC
      Uhalf[i].E -= hdtodx1*(x1Flux[i  ].d*(phic - phil) +
                             x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 6b ------------------------------------------------------------------
 * Add source terms for self gravity for 0.5*dt to predict step.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (i=il; i<=iu; i++) {
    phic = pG->Phi[ks][js][i];
    phir = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i+1]);
    phil = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i-1]);

    Uhalf[i].M1 -= hdtodx1*pG->U[ks][js][i].d*(phir-phil);
#ifndef BAROTROPIC
    Uhalf[i].E -= hdtodx1*(x1Flux[i  ].d*(phic - phil) +
                           x1Flux[i+1].d*(phir - phic));
#endif
  }
#endif /* SELF_GRAVITY */

/*=== STEP 7: Compute primitive variables at half timestep  ===================*/

/*--- Step 7a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=il; i<=iu; i++) {
    U1d[i].d  = Uhalf[i].d;
    U1d[i].Mx = Uhalf[i].M1;
    U1d[i].My = Uhalf[i].M2;
    U1d[i].Mz = Uhalf[i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = Uhalf[i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    Bxc[i]    = Uhalf[i].B1c; /* although in 1D Bx does not evolve */
    U1d[i].By = Uhalf[i].B2c;
    U1d[i].Bz = Uhalf[i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = Uhalf[i].s[n];
#endif

/*--- Step 7b ------------------------------------------------------------------
 * With special relativity, load primitive variables at last timestep into W
 * as first guess in iterative conversion method */

#ifdef SPECIAL_RELATIVITY
    W[i].d  = pG->W[ks][js][i].d;
    W[i].Vx = pG->W[ks][js][i].V1;
    W[i].Vy = pG->W[ks][js][i].V2;
    W[i].Vz = pG->W[ks][js][i].V3;
#ifndef BAROTROPIC
    W[i].P  = pG->W[ks][js][i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    W[i].By = pG->W[ks][js][i].B2c;
    W[i].Bz = pG->W[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W[i].r[n] = pG->W[ks][js][i].r[n];
#endif /* NSCALARS */
#endif /* SPECIAL_RELATIVITY */

/*--- Step 7c ------------------------------------------------------------------
 * convert conserved to primitive, and store in array Whalf */

    Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));

    Whalf[i].d  = W[i].d;
    Whalf[i].V1 = W[i].Vx;
    Whalf[i].V2 = W[i].Vy;
    Whalf[i].V3 = W[i].Vz;
#ifndef BAROTROPIC
    Whalf[i].P  = W[i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    Whalf[i].B1c = Bxc[i];
    Whalf[i].B2c = W[i].By;
    Whalf[i].B3c = W[i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) Whalf[i].r[n] = W[i].r[n];
#endif /* NSCALARS */
  }

/*=== STEP 8: Compute second-order L/R x1-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of primitive variables;
 * W = (d, V1, V2, V3, E, B2c, B3c, s[n])
 */

  for (i=il; i<=iu; i++) {
    W[i].d  = Whalf[i].d;
    W[i].Vx = Whalf[i].V1;
    W[i].Vy = Whalf[i].V2;
    W[i].Vz = Whalf[i].V3;
#ifndef BAROTROPIC
    W[i].P  = Whalf[i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    Bxc[i]  = Whalf[i].B1c;
    W[i].By = Whalf[i].B2c;
    W[i].Bz = Whalf[i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W[i].r[n] = Whalf[i].r[n];
#endif /* NSCALARS */
  }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L/R states on x1-interfaces, store into arrays
 */

  lr_states(W, MHDARG( Bxc , ) 0.0,is,ie,Wl,Wr);

  for (i=is; i<=ie+1; i++) {
    Wl_x1Face[i] = Wl[i];
    Wr_x1Face[i] = Wr[i];
  }

/*=== STEPS 9-10: Not needed in 1D ===*/

/*=== STEP 11: Compute x1-Flux ===============================================*/

/*--- Step 11b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (i=is; i<=ie+1; i++) {
    Prim1D_to_Cons1D(&Ul[i],&Wl_x1Face[i] MHDARG( , &B1_x1Face[i]));
    Prim1D_to_Cons1D(&Ur[i],&Wr_x1Face[i] MHDARG( , &B1_x1Face[i]));
    fluxes(Ul[i],Ur[i],Wl_x1Face[i],Wr_x1Face[i],
      MHDARG(B1_x1Face[i], ) &x1Flux[i]);
  }

/*=== STEP 12: Not needed in 1D ===*/
        
/*=== STEP 13: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 13a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      pG->U[ks][js][i].M1 -= dtodx1*Uhalf[i].d*(phir-phil);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 13b -----------------------------------------------------------------
 * Add gravitational source terms for self-gravity.
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */

#ifdef SELF_GRAVITY
/* Add fluxes and source terms due to (d/dx1) terms  */

  for (i=is; i<=ie; i++){
    phic = pG->Phi[ks][js][i];
    phil = 0.5*(pG->Phi[ks][js][i-1] + pG->Phi[ks][js][i  ]);
    phir = 0.5*(pG->Phi[ks][js][i  ] + pG->Phi[ks][js][i+1]);

/* gx, gy and gz centered at L and R x1-faces */
    gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
    gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
    flx_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
    flx_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx1 terms  */
    pG->U[ks][js][i].M1 -= dtodx1*(flx_m1r - flx_m1l);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                  x1Flux[i+1].d*(phir - phic));
#endif /* BAROTROPIC */
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */

  for (i=is; i<=ie+1; i++) {
    pG->x1MassFlux[ks][js][i] = x1Flux[i].d;
  }
#endif /* SELF_GRAVITY */

/*=== STEP 14: Update cell-centered values for a full timestep ===============*/

/*--- Step 14a -----------------------------------------------------------------
 * Update cell-centered variables in pG (including B2c and B3c) using x1-Fluxes
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
/* For consistency, set B2i and B3i to cell-centered values.  */
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pG->U[ks][js][i].s[n] -= dtodx1*(x1Flux[i+1].s[n] - x1Flux[i].s[n]);
#endif
  }

/*--- Step 14b -----------------------------------------------------------------
 * With special relativity, compute primitive variables at new timestep and
 * store in Grid structure
 */

#ifdef SPECIAL_RELATIVITY

/* Load 1D vector of conserved variables */

  for (i=is; i<=ie; i++) {
    U1d[i].d  = pG->U[ks][js][i].d;
    U1d[i].Mx = pG->U[ks][js][i].M1;
    U1d[i].My = pG->U[ks][js][i].M2;
    U1d[i].Mz = pG->U[ks][js][i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = pG->U[ks][js][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    Bxc[i]    = pG->U[ks][js][i].B1c;
    U1d[i].By = pG->U[ks][js][i].B2c;
    U1d[i].Bz = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif

/* Load 1D vector of primitive variables from Whalf as initial guess */

    W[i].d  = Whalf[i].d;
    W[i].Vx = Whalf[i].V1;
    W[i].Vy = Whalf[i].V2;
    W[i].Vz = Whalf[i].V3;
#ifndef BAROTROPIC
    W[i].P  = Whalf[i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    W[i].By = Whalf[i].B2c;
    W[i].Bz = Whalf[i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) W[i].r[n] = Whalf[i].r[n];
#endif /* NSCALARS */

/* Convert variables and store result into Grid */

    Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));

    pG->W[ks][js][i].d  = W[i].d;
    pG->W[ks][js][i].V1 = W[i].Vx;
    pG->W[ks][js][i].V2 = W[i].Vy;
    pG->W[ks][js][i].V3 = W[i].Vz;
#ifndef BAROTROPIC
    pG->W[ks][js][i].P = W[i].P;
#endif /* BAROTROPIC */
#ifdef MHD
    pG->W[ks][js][i].B1c = Bxc[i];
    pG->W[ks][js][i].B2c = W[i].By;
    pG->W[ks][js][i].B3c = W[i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) pG->W[ks][js][i].r[n] = W[i].r[n];
#endif /* NSCALARS */
  }

#endif /* SPECIAL_RELATIVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(MeshS *pM)
{
  int size1=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->N[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;

  if ((Wl_x1Face = (Prim1D*)malloc(size1*sizeof(Prim1D))) ==NULL) goto on_error;
  if ((Wr_x1Face = (Prim1D*)malloc(size1*sizeof(Prim1D))) ==NULL) goto on_error;
  if ((x1Flux    = (Cons1D*)malloc(size1*sizeof(Cons1D))) ==NULL) goto on_error;

#ifdef MHD
  if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((B1_x1Face = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
#endif /* MHD */

  if ((U1d =      (Cons1D*)malloc(size1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul  =      (Cons1D*)malloc(size1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur  =      (Cons1D*)malloc(size1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((W   =      (Prim1D*)malloc(size1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl  =      (Prim1D*)malloc(size1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr  =      (Prim1D*)malloc(size1*sizeof(Prim1D))) == NULL) goto on_error;

  if ((Uhalf = (Gas*)malloc(size1*sizeof(Gas))) == NULL) goto on_error;
  if ((Whalf = (Prim*)malloc(size1*sizeof(Prim))) == NULL) goto on_error;

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_1d:  Free temporary integration arrays */

void integrate_destruct_1d(void)
{
  if (Wl_x1Face != NULL) free(Wl_x1Face);
  if (Wr_x1Face != NULL) free(Wr_x1Face);
  if (x1Flux    != NULL) free(x1Flux);

#ifdef MHD
  if (Bxc != NULL) free(Bxc);
  if (B1_x1Face != NULL) free(B1_x1Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Uhalf    != NULL) free(Uhalf);
  if (Whalf    != NULL) free(Whalf);

  return;
}
#endif /* VL_INTEGRATOR */
