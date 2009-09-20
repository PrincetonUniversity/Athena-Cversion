#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_1d_ctu.c
 *
 * PURPOSE: Integrate MHD equations using 1D version of the CTU integrator.
 *   Updates U.[d,M1,M2,M3,E,B2c,B3c,s] in Grid structure, where U is of type
 *   Gas.  Adds gravitational source terms, self-gravity, and optically-thin
 *   cooling.
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
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(CTU_INTEGRATOR) && defined(CARTESIAN)

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1D *Ul_x1Face=NULL, *Ur_x1Face=NULL, *x1Flux=NULL;

/* 1D scratch vectors used by lr_states and flux functions */
#ifdef MHD
static Real *Bxc=NULL, *Bxi=NULL;
#endif /* MHD */
static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1D *U1d=NULL;

/* Variables at t^{n+1/2} used for source terms */
static Real *dhalf = NULL, *phalf = NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_1d: 1D version of CTU unsplit integrator for MHD
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 1D.
 */

void integrate_1d(Grid *pG, Domain *pD)
{
  Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
  int i,il,iu, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h;
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolf,Eh=0.0;
#endif
#if defined(MHD) && !defined(BAROTROPIC)
  Real B1ch,B2ch,B3ch;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,flux_m1l,flux_m1r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif

/* With particles, one more ghost cell must be updated in predict step */
#ifdef PARTICLES
  Real d1;
  il = is - 3;
  iu = ie + 3;
#else
  il = is - 1;
  iu = ie + 1;
#endif

/* Compute predictor feedback from particle drag */
#ifdef FEEDBACK
  feedback_predictor(pG);
#endif

/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

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
    Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));
  }

  lr_states(W, MHDARG( Bxc , ) dtodx1,il+1,iu-1,Wl,Wr);

/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=il+1; i<=iu; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phicr = (*StaticGravPot)( x1                ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
      phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);
    }
  }

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
  for (i=il+1; i<=iu; i++) {
    Wl[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
    Wr[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
  }
#endif

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (i=il+1; i<=iu; i++) {
      coolfl = (*CoolingFunc)(Wl[i].d,Wl[i].P,(0.5*pG->dt));
      coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));

      Wl[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
      Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;
    }
  }
#endif /* BAROTROPIC */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {
      d1 = 1.0/W[i-1].d;
      Wl[i].Vx -= pG->feedback[ks][js][i-1].x1*d1;
      Wl[i].Vy -= pG->feedback[ks][js][i-1].x2*d1;
      Wl[i].Vz -= pG->feedback[ks][js][i-1].x3*d1;

      d1 = 1.0/W[i].d;
      Wr[i].Vx -= pG->feedback[ks][js][i].x1*d1;
      Wr[i].Vy -= pG->feedback[ks][js][i].x2*d1;
      Wr[i].Vz -= pG->feedback[ks][js][i].x3*d1;

#ifndef BAROTROPIC
      Wl[i].P += pG->Eloss[ks][js][i-1]*Gamma_1;
      Wr[i].P += pG->Eloss[ks][js][i]*Gamma_1;
#endif
    }
#endif /* FEEDBACK */


/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction
 */

  for (i=il+1; i<=iu; i++) {
    Prim1D_to_Cons1D(&Ul_x1Face[i],&Wl[i] MHDARG( , &Bxi[i]));
    Prim1D_to_Cons1D(&Ur_x1Face[i],&Wr[i] MHDARG( , &Bxi[i]));

    fluxes(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],
               MHDARG( Bxi[i] , ) &x1Flux[i]);
  }

/*=== STEPS 2-7: Not needed in 1D ===*/

/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 8a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, or cooling)
 */

#ifndef PARTICLES
  if ((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
  {
    for (i=il+1; i<=iu-1; i++) {
      dhalf[i] = pG->U[ks][js][i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
#ifdef PARTICLES
      grid_d[ks][js][i] = dhalf[i];
#endif
    }
  }

/*--- Step 8b ------------------------------------------------------------------
 * Calculate P^{n+1/2} (needed with cooling)
 */

#ifndef PARTICLES
  if (CoolingFunc != NULL)
#endif /* PARTICLES */
  {
    for (i=il+1; i<=iu-1; i++) {
      M1h = pG->U[ks][js][i].M1 - hdtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
      M2h = pG->U[ks][js][i].M2 - hdtodx1*(x1Flux[i+1].My - x1Flux[i].My);
      M3h = pG->U[ks][js][i].M3 - hdtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
      Eh  = pG->U[ks][js][i].E  - hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif

/* Add source terms for fixed gravitational potential */
      if (StaticGravPot != NULL){
        cc_pos(pG,i,js,ks,&x1,&x2,&x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
        M1h -= hdtodx1*(phir-phil)*pG->U[ks][js][i].d;
      }

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      phir = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i+1]);
      phil = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i-1]);
      M1h -= hdtodx1*(phir-phil)*pG->U[ks][js][i].d;
#endif /* SELF_GRAVITY */

/* Add the particle feedback terms */
#ifdef FEEDBACK
      M1h -= pG->feedback[ks][js][i].x1;
      M2h -= pG->feedback[ks][js][i].x2;
      M3h -= pG->feedback[ks][js][i].x3;
#endif /* FEEDBACK */

#ifndef BAROTROPIC
      phalf[i] = Eh - 0.5*(M1h*M1h + M2h*M2h + M3h*M3h)/dhalf[i];

#ifdef MHD
      B1ch = pG->U[ks][js][i].B1c;
      B2ch = pG->U[ks][js][i].B2c - hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
      B3ch = pG->U[ks][js][i].B3c - hdtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
      phalf[i] -= 0.5*(B1ch*B1ch + B2ch*B2ch + B3ch*B3ch);
#endif /* MHD */

      phalf[i] *= Gamma_1;
#endif /* BAROTROPIC */

#ifdef PARTICLES
      d1 = 1.0/dhalf[i];
      grid_v[ks][js][i].x1 = M1h*d1;
      grid_v[ks][js][i].x2 = M2h*d1;
      grid_v[ks][js][i].x3 = M3h*d1;
#ifndef ISOTHERMAL
      grid_cs[ks][js][i] = sqrt(Gamma*phalf[i]*d1);
#endif  /* ISOTHERMAL */
#endif /* PARTICLES */

    }
  }

/*=== STEP 8.5: Integrate the particles, compute the feedback ================*/

#ifdef PARTICLES
  Integrate_Particles(pG,pD);
#ifdef FEEDBACK
  exchange_feedback(pG, pD);
#endif
#endif

/*=== STEPS 9-10: Not needed in 1D ===*/

/*=== STEP 11: Add source terms for a full timestep using n+1/2 states =======*/

/*--- Step 11a -----------------------------------------------------------------
 * Add gravitational source terms as a Static Potential.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      pG->U[ks][js][i].M1 -= dtodx1*dhalf[i]*(phir-phil);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 11b -----------------------------------------------------------------
 * Add source terms for self-gravity.
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
#endif /* SELF_GRAVITY */

/*--- Step 11c -----------------------------------------------------------------
 * Add source terms for optically thin cooling
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (i=is; i<=ie; i++) {
      coolf = (*CoolingFunc)(dhalf[i],phalf[i],pG->dt);
      pG->U[ks][js][i].E -= pG->dt*coolf;
    }
  }
#endif /* BAROTROPIC */

/*--- Step 11d -----------------------------------------------------------------
 * Add source terms for particle feedback
 */

#ifdef FEEDBACK
  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].M1 -= pG->feedback[ks][js][i].x1;
    pG->U[ks][js][i].M2 -= pG->feedback[ks][js][i].x2;
    pG->U[ks][js][i].M3 -= pG->feedback[ks][js][i].x3;
#ifndef BAROTROPIC
    pG->U[ks][js][i].E += pG->Eloss[ks][js][i];
    pG->Eloss[ks][js][i] *= dt1;
#endif
  }
#endif

/*=== STEP 12: Update cell-centered values for a full timestep ===============*/

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 1D x1-fluxes
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

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(int nx1)
{
  int Nx1;
  Nx1 = nx1 + 2*nghost;

#ifdef MHD
  if ((Bxc = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
#endif /* MHD */

  if ((U1d       = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul_x1Face = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur_x1Face = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;
  if ((x1Flux    = (Cons1D*)malloc(Nx1*sizeof(Cons1D))) == NULL) goto on_error;

  if ((W  = (Prim1D*)malloc(Nx1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl = (Prim1D*)malloc(Nx1*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr = (Prim1D*)malloc(Nx1*sizeof(Prim1D))) == NULL) goto on_error;

  if((StaticGravPot != NULL) || (CoolingFunc != NULL)){
    if ((dhalf  = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  }
  if(CoolingFunc != NULL){
    if ((phalf  = (Real*)malloc(Nx1*sizeof(Real))) == NULL) goto on_error;
  }

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
#endif /* MHD */

  if (U1d != NULL) free(U1d);
  if (Ul_x1Face != NULL) free(Ul_x1Face);
  if (Ur_x1Face != NULL) free(Ur_x1Face);
  if (x1Flux != NULL) free(x1Flux);

  if (W  != NULL) free(W);
  if (Wl != NULL) free(Wl);
  if (Wr != NULL) free(Wr);

  if (dhalf != NULL) free(dhalf);
  if (phalf != NULL) free(phalf);

  return;
}
#endif /* CTU_INTEGRATOR */
