#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_2d_ctu.c
 *
 * PURPOSE: Integrate MHD equations in 2D using the directionally unsplit CTU
 *   method of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type Gas
 *      B1i, B2i -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically-thin cooling,
 *   and H-correction of Sanders et al.
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport", JCP, 205, 509 (2005)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_2d()
 *   integrate_init_2d()
 *   integrate_destruct_2d()
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

#ifdef CTU_INTEGRATOR

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1D **Ul_x1Face=NULL, **Ur_x1Face=NULL;
static Cons1D **Ul_x2Face=NULL, **Ur_x2Face=NULL;
static Cons1D **x1Flux=NULL, **x2Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
#ifdef MHD
static Real *Bxc=NULL, *Bxi=NULL;
#endif /* MHD */
static Prim1D *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1D *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real **dhalf = NULL,**phalf = NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method in Gardiner & Stone (2005) 
 *============================================================================*/

#ifdef MHD
static void integrate_emf3_corner(Grid *pG);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_2d:  CTU integrator in 2D.
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */

void integrate_2d(Grid *pG, Domain *pD)
{
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  Real hdt = 0.5*pG->dt;
  int i,il,iu,is=pG->is, ie=pG->ie;
  int j,jl,ju,js=pG->js, je=pG->je;
  int ks=pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,d1;
  Real coolfl,coolfr,coolf,M1h,M2h,M3h,Eh=0.0;
#ifdef MHD
  Real MHD_src,dbx,dby,B1,B2,B3,V3;
  Real B1ch, B2ch, B3ch;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flux_m1l,flux_m1r,flux_m2l,flux_m2r;
#endif
#ifdef SHEARING_BOX
  Real M1n, dM3n;   /* M1, dM3=(My+d*q*Omega_0*x) at time n */
  Real M1e, dM3e;   /* M1, dM3 evolved by dt/2  */
  Real flx1_dM3, frx1_dM3, flx2_dM3, frx2_dM3;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */

/* With particles, one more ghost cell must be updated in predict step */
#ifdef PARTICLES
  il = is - 3;
  iu = ie + 3;
  jl = js - 3;
  ju = je + 3;
#else
  il = is - 2;
  iu = ie + 2;
  jl = js - 2;
  ju = je + 2;
#endif

/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (j=jl; j<=ju; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      U1d[i].d  = pG->U[ks][j][i].d;
      U1d[i].Mx = pG->U[ks][j][i].M1;
      U1d[i].My = pG->U[ks][j][i].M2;
      U1d[i].Mz = pG->U[ks][j][i].M3;
#ifndef BAROTROPIC
      U1d[i].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[i].By = pG->U[ks][j][i].B2c;
      U1d[i].Bz = pG->U[ks][j][i].B3c;
      Bxc[i] = pG->U[ks][j][i].B1c;
      Bxi[i] = pG->B1i[ks][j][i];
      B1_x1Face[j][i] = pG->B1i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (i=is-nghost; i<=ie+nghost; i++) {
      Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));
    }

    lr_states(W, MHDARG( Bxc , ) pG->dt,dtodx1,il+1,iu-1,Wl,Wr);

#ifdef MHD
    for (i=il+1; i<=iu; i++) {
      MHD_src = (pG->U[ks][j][i-1].M2/pG->U[ks][j][i-1].d)*
               (pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1])/pG->dx1;
      Wl[i].By += hdt*MHD_src;

      MHD_src = (pG->U[ks][j][i].M2/pG->U[ks][j][i].d)*
               (pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i])/pG->dx1;
      Wr[i].By += hdt*MHD_src;
    }
#endif

/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

    if (StaticGravPot != NULL){
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*StaticGravPot)( x1             ,x2,x3);
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
      Wl[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      Wr[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
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
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
 * The tidal gravity terms are added through the StaticGravPot
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
 */

#ifdef SHEARING_BOX
    for (i=il+1; i<=iu; i++) {
      Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vz;
#ifdef FARGO
      Wl[i].Vz += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
#else
      Wl[i].Vz -= pG->dt*Omega_0*W[i-1].Vx;
#endif

      Wr[i].Vx += pG->dt*Omega_0*W[i].Vz; 
#ifdef FARGO
      Wr[i].Vz += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
      Wr[i].Vz -= pG->dt*Omega_0*W[i].Vx;
#endif
    }
#endif /* SHEARING_BOX */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {
      d1 = 1.0/W[i-1].d;
      Wl[i].Vx -= pG->feedback[ks][j][i-1].x1*d1;
      Wl[i].Vy -= pG->feedback[ks][j][i-1].x2*d1;
      Wl[i].Vz -= pG->feedback[ks][j][i-1].x3*d1;

      d1 = 1.0/W[i].d;
      Wr[i].Vx -= pG->feedback[ks][j][i].x1*d1;
      Wr[i].Vy -= pG->feedback[ks][j][i].x2*d1;
      Wr[i].Vz -= pG->feedback[ks][j][i].x3*d1;
    }
#endif /* FEEDBACK */

/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 2D array
 */

    for (i=il+1; i<=iu; i++) {
      Prim1D_to_Cons1D(&Ul_x1Face[j][i],&Wl[i] MHDARG( , &Bxi[i]));
      Prim1D_to_Cons1D(&Ur_x1Face[j][i],&Wr[i] MHDARG( , &Bxi[i]));

      fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],
                 MHDARG( B1_x1Face[j][i] , ) &x1Flux[j][i]);
    }

  }

/*=== STEP 2: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (i=il; i<=iu; i++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      U1d[j].d  = pG->U[ks][j][i].d;
      U1d[j].Mx = pG->U[ks][j][i].M2;
      U1d[j].My = pG->U[ks][j][i].M3;
      U1d[j].Mz = pG->U[ks][j][i].M1;
#ifndef BAROTROPIC
      U1d[j].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[j].By = pG->U[ks][j][i].B3c;
      U1d[j].Bz = pG->U[ks][j][i].B1c;
      Bxc[j] = pG->U[ks][j][i].B2c;
      Bxi[j] = pG->B2i[ks][j][i];
      B2_x2Face[j][i] = pG->B2i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (j=js-nghost; j<=je+nghost; j++) {
      Cons1D_to_Prim1D(&U1d[j],&W[j] MHDARG( , &Bxc[j]));
    }

    lr_states(W, MHDARG( Bxc , ) pG->dt,dtodx2,jl+1,ju-1,Wl,Wr);

#ifdef MHD
    for (j=jl+1; j<=ju; j++) {
      MHD_src = (pG->U[ks][j-1][i].M1/pG->U[ks][j-1][i].d)*
        (pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i])/pG->dx2;
      Wl[j].Bz += hdt*MHD_src;

      MHD_src = (pG->U[ks][j][i].M1/pG->U[ks][j][i].d)*
        (pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i])/pG->dx2;
      Wr[j].Bz += hdt*MHD_src;
    }
#endif

/*--- Step 2c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
  
    if (StaticGravPot != NULL){
      for (j=jl+1; j<=ju; j++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*StaticGravPot)(x1, x2             ,x3);
        phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
        phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Wl[j].Vx -= dtodx2*(phifc - phicl);
        Wr[j].Vx -= dtodx2*(phicr - phifc);
      }
    }

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    for (j=jl+1; j<=ju; j++) {
      Wl[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      Wr[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
    }
#endif

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
    if (CoolingFunc != NULL){
      for (j=jl+1; j<=ju; j++) {
        coolfl = (*CoolingFunc)(Wl[j].d,Wl[j].P,(0.5*pG->dt));
        coolfr = (*CoolingFunc)(Wr[j].d,Wr[j].P,(0.5*pG->dt));

        Wl[j].P -= 0.5*pG->dt*Gamma_1*coolfl;
        Wr[j].P -= 0.5*pG->dt*Gamma_1*coolfr;
      }
    }
#endif /* BAROTROPIC */

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
   for (j=jl+1; j<=ju; j++) {
      d1 = 1.0/W[j-1].d;
      Wl[j].Vx -= pG->feedback[ks][j-1][i].x2*d1;
      Wl[j].Vy -= pG->feedback[ks][j-1][i].x3*d1;
      Wl[j].Vz -= pG->feedback[ks][j-1][i].x1*d1;

      d1 = 1.0/W[j].d;
      Wr[j].Vx -= pG->feedback[ks][j][i].x2*d1;
      Wr[j].Vy -= pG->feedback[ks][j][i].x3*d1;
      Wr[j].Vz -= pG->feedback[ks][j][i].x1*d1;
    }
#endif /* FEEDBACK */

/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 2D array
 */

    for (j=jl+1; j<=ju; j++) {
      Prim1D_to_Cons1D(&Ul_x2Face[j][i],&Wl[j] MHDARG( , &Bxi[j]));
      Prim1D_to_Cons1D(&Ur_x2Face[j][i],&Wr[j] MHDARG( , &Bxi[j]));

      fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],
                 MHDARG( B2_x2Face[j][i] , ) &x2Flux[j][i]);
    }
  }

/*=== STEP 3: Not needed in 2D ===*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf_3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      emf3_cc[j][i] =
	(pG->U[ks][j][i].B1c*pG->U[ks][j][i].M2 -
	 pG->U[ks][j][i].B2c*pG->U[ks][j][i].M1 )/pG->U[ks][j][i].d;
    }
  }
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    B1_x1Face[j][iu] -= hdtodx2*(emf3[j+1][iu] - emf3[j][iu]);
  }
  for (i=il+1; i<=iu-1; i++) {
    B2_x2Face[ju][i] += hdtodx1*(emf3[ju][i+1] - emf3[ju][i]);
  }
#endif

/*=== STEP 5: Correct x1-interface states with transverse flux gradients =====*/

/*--- Step 5a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
      Ul_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i-1].d  - x2Flux[j][i-1].d );
      Ul_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i-1].Mz - x2Flux[j][i-1].Mz);
      Ul_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i-1].Mx - x2Flux[j][i-1].Mx);
      Ul_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i-1].My - x2Flux[j][i-1].My);
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i-1].E  - x2Flux[j][i-1].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ul_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i-1].By - x2Flux[j][i-1].By);
#endif

      Ur_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i  ].d  - x2Flux[j][i  ].d );
      Ur_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i  ].Mz - x2Flux[j][i  ].Mz);
      Ur_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i  ].Mx - x2Flux[j][i  ].Mx);
      Ur_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i  ].My - x2Flux[j][i  ].My);
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i  ].E  - x2Flux[j][i  ].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ur_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i  ].By - x2Flux[j][i  ].By);
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) {
      Ul_x1Face[j][i].s[n]-=hdtodx2*(x2Flux[j+1][i-1].s[n]-x2Flux[j][i-1].s[n]);
      Ur_x1Face[j][i].s[n]-=hdtodx2*(x2Flux[j+1][i  ].s[n]-x2Flux[j][i  ].s[n]);
      }
#endif
    }
  }

/*--- Step 5b: Not needed in 2D ---
/*--- Step 5c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2-flux gradient to the conservative
 * variables on the x1Face..
 */

#ifdef MHD
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
      dbx = pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1];
      B1 = pG->U[ks][j][i-1].B1c;
      B2 = pG->U[ks][j][i-1].B2c;
      B3 = pG->U[ks][j][i-1].B3c;
      V3 = pG->U[ks][j][i-1].M3/pG->U[ks][j][i-1].d;

      Ul_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ul_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ul_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ul_x1Face[j][i].Bz += hdtodx1*V3*dbx;
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;
#endif /* BAROTROPIC */

      dbx = pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i];
      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ur_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ur_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ur_x1Face[j][i].Bz += hdtodx1*V3*dbx;
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;
#endif /* BAROTROPIC */
    }
  }
#endif /* MHD */

/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x2-Flux
 * gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)(x1, x2             ,x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                      x2Flux[j+1][i  ].d*(phir - phic));
#endif

        phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);
        
        Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;
#ifndef BAROTROPIC
        Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                      x2Flux[j+1][i-1].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 5d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x2-Flux gradient
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);

      Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                    x2Flux[j+1][i  ].d*(phir - phic));
#endif

      phic = pG->Phi[ks][j][i-1];
      phir = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j+1][i-1]);
      phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j-1][i-1]);

      Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;
#ifndef BAROTROPIC
      Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                    x2Flux[j+1][i-1].d*(phir - phic));
#endif
    }
  }
#endif /* SELF_GRAVITY */

/*=== STEP 6: Correct x2-interface states with transverse flux gradients =====*/

/*--- Step 6a ------------------------------------------------------------------
 * Correct x2-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */

  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      Ul_x2Face[j][i].d  -= hdtodx1*(x1Flux[j-1][i+1].d  - x1Flux[j-1][i].d );
      Ul_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j-1][i+1].My - x1Flux[j-1][i].My);
      Ul_x2Face[j][i].My -= hdtodx1*(x1Flux[j-1][i+1].Mz - x1Flux[j-1][i].Mz);
      Ul_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j-1][i+1].Mx - x1Flux[j-1][i].Mx);
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E  -= hdtodx1*(x1Flux[j-1][i+1].E  - x1Flux[j-1][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ul_x2Face[j][i].By -= hdtodx1*(x1Flux[j-1][i+1].Bz - x1Flux[j-1][i].Bz);
#endif

      Ur_x2Face[j][i].d  -= hdtodx1*(x1Flux[j  ][i+1].d  - x1Flux[j  ][i].d );
      Ur_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j  ][i+1].My - x1Flux[j  ][i].My);
      Ur_x2Face[j][i].My -= hdtodx1*(x1Flux[j  ][i+1].Mz - x1Flux[j  ][i].Mz);
      Ur_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j  ][i+1].Mx - x1Flux[j  ][i].Mx);
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E  -= hdtodx1*(x1Flux[j  ][i+1].E  - x1Flux[j  ][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      Ur_x2Face[j][i].By -= hdtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz);
#endif
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) {
      Ul_x2Face[j][i].s[n]-=hdtodx1*(x1Flux[j-1][i+1].s[n]-x1Flux[j-1][i].s[n]);
      Ur_x2Face[j][i].s[n]-=hdtodx1*(x1Flux[j  ][i+1].s[n]-x1Flux[j  ][i].s[n]);
      }
#endif
    }
  }

/*--- Step 6b: Not needed in 2D ---
/*--- Step 6c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1-flux-gradients to the
 * conservative variables on the x2Face.
 */

#ifdef MHD
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      dby = pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i];
      B1 = pG->U[ks][j-1][i].B1c;
      B2 = pG->U[ks][j-1][i].B2c;
      B3 = pG->U[ks][j-1][i].B3c;
      V3 = pG->U[ks][j-1][i].M3/pG->U[ks][j-1][i].d;

      Ul_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ul_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ul_x2Face[j][i].My += hdtodx2*B3*dby;
      Ul_x2Face[j][i].By += hdtodx2*V3*dby;
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E  += hdtodx2*B3*V3*dby;
#endif /* BAROTROPIC */

      dby = pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i];
      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ur_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ur_x2Face[j][i].My += hdtodx2*B3*dby;
      Ur_x2Face[j][i].By += hdtodx2*V3*dby;
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E  += hdtodx2*B3*V3*dby;
#endif /* BAROTROPIC */
    }
  }
#endif /* MHD */

/*--- Step 6d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * gradients. To improve conservation of total energy, average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                      x1Flux[j  ][i+1].d*(phir - phic));
#endif

        phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

        Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#ifndef BAROTROPIC
        Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                      x1Flux[j-1][i+1].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x1-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);

      Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                    x1Flux[j  ][i+1].d*(phir - phic));
#endif

      phic = pG->Phi[ks][j-1][i];
      phir = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i+1]);
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i-1]);

      Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#ifndef BAROTROPIC
      Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                    x1Flux[j-1][i+1].d*(phir - phic));
#endif
    }
  }

#endif /* SELF_GRAVITY */

/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
 * x1-Flux gradient.  The tidal gravity terms are added via StaticGravPot
 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in shearing sheet
 */

#ifdef SHEARING_BOX
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      Ur_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j][i].M3;
#ifdef FARGO
      Ur_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
#else
      Ur_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
#endif

      Ul_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j-1][i].M3;
#ifdef FARGO
      Ul_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j-1][i].M1;
#else
      Ul_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j-1][i].M1;
#endif
    }
  }
#endif /* SHEARING_BOX */

/*=== STEP 7: Not needed in 2D ===*/

/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 8a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, cooling, or MHD)
 */

#ifndef MHD
#ifndef PARTICLES
  if ((StaticGravPot != NULL) || (CoolingFunc != NULL)) 
#endif
#endif
  {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {
        dhalf[j][i] = pG->U[ks][j][i].d
          - hdtodx1*(x1Flux[j  ][i+1].d - x1Flux[j][i].d)
          - hdtodx2*(x2Flux[j+1][i  ].d - x2Flux[j][i].d);
#ifdef PARTICLES
        grid_d[ks][j][i] = dhalf[j][i];
#endif
      }
    }
  }

/*--- Step 8b ------------------------------------------------------------------
 * Calculate P^{n+1/2} (needed with cooling), and cell centered emf3_cc^{n+1/2}
 */

#ifndef MHD
#ifndef PARTICLES
  if (CoolingFunc != NULL) 
#endif /* PARTICLES */
#endif /* MHD */
  {
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      M1h = pG->U[ks][j][i].M1
        - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
        - hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);

      M2h = pG->U[ks][j][i].M2
        - hdtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My)
        - hdtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);

      M3h = pG->U[ks][j][i].M3
        - hdtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz)
        - hdtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);

#ifndef BAROTROPIC
      Eh = pG->U[ks][j][i].E
        - hdtodx1*(x1Flux[j][i+1].E - x1Flux[j][i].E)
        - hdtodx2*(x2Flux[j+1][i].E - x2Flux[j][i].E);
#endif

/* Add source terms for fixed gravitational potential */
      if (StaticGravPot != NULL){
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
        M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
        M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
      }

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);
      M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);
      M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#endif /* SELF_GRAVITY */

/* Add the Coriolis terms for shearing box. */
#ifdef SHEARING_BOX
      M1h += pG->dt*Omega_0*pG->U[ks][j][i].M3;
#ifdef FARGO
      M3h += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
#else
      M3h -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
#endif

#endif /* SHEARING_BOX */

/* Add the particle feedback terms */
#ifdef FEEDBACK
      M1h -= pG->feedback[ks][j][i].x1;
      M2h -= pG->feedback[ks][j][i].x2;
      M3h -= pG->feedback[ks][j][i].x3;
#endif /* FEEDBACK */

#ifndef BAROTROPIC
      phalf[j][i] = Eh - 0.5*(M1h*M1h + M2h*M2h + M3h*M3h)/dhalf[j][i];
#endif

#ifdef MHD
      B1ch = 0.5*(B1_x1Face[j][i] + B1_x1Face[j][i+1]);
      B2ch = 0.5*(B2_x2Face[j][i] + B2_x2Face[j+1][i]);
      B3ch = pG->U[ks][j][i].B3c 
        - hdtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz)
        - hdtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
      emf3_cc[j][i] = (B1ch*M2h - B2ch*M1h)/dhalf[j][i];
#ifndef BAROTROPIC
      phalf[j][i] -= 0.5*(B1ch*B1ch + B2ch*B2ch + B3ch*B3ch);
#endif
#endif /* MHD */

#ifndef BAROTROPIC
      phalf[j][i] *= Gamma_1;
#endif

#ifdef PARTICLES
      d1 = 1.0/dhalf[j][i];
      grid_v[ks][j][i].x1 = M1h*d1;
      grid_v[ks][j][i].x2 = M2h*d1;
      grid_v[ks][j][i].x3 = M3h*d1;
#ifndef ISOTHERMAL
  #ifdef ADIABATIC
      grid_cs[ks][j][i] = sqrt(Gamma*phalf[j][i]*d1);
  #else
      ath_error("[get_gasinfo] can not calculate the sound speed!\n");
  #endif /* ADIABATIC */
#endif  /* ISOTHERMAL */
#endif /* PARTICLES */
    }
  }
  }

/*=== STEP 9: Compute 2D x1-Flux, x2-Flux ====================================*/

/*--- Step 9a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

#ifdef H_CORRECTION
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+2; i++) {
      cfr = cfast(&(Ur_x1Face[j][i]) MHDARG( , &(B1_x1Face[j][i])));
      cfl = cfast(&(Ul_x1Face[j][i]) MHDARG( , &(B1_x1Face[j][i])));
      lambdar = Ur_x1Face[j][i].Mx/Ur_x1Face[j][i].d + cfr;
      lambdal = Ul_x1Face[j][i].Mx/Ul_x1Face[j][i].d - cfl;
      eta1[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }

  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+1; i++) {
      cfr = cfast(&(Ur_x2Face[j][i]) MHDARG( , &(B2_x2Face[j][i])));
      cfl = cfast(&(Ul_x2Face[j][i]) MHDARG( , &(B2_x2Face[j][i])));
      lambdar = Ur_x2Face[j][i].Mx/Ur_x2Face[j][i].d + cfr;
      lambdal = Ul_x2Face[j][i].Mx/Ul_x2Face[j][i].d - cfl;
      eta2[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }
#endif /* H_CORRECTION */

/*--- Step 9b ------------------------------------------------------------------
 * Compute 2D x1-fluxes from corrected L/R states.
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta2[j][i-1],eta2[j][i]);
      etah = MAX(etah,eta2[j+1][i-1]);
      etah = MAX(etah,eta2[j+1][i  ]);
      etah = MAX(etah,eta1[j  ][i  ]);
#endif /* H_CORRECTION */
      Cons1D_to_Prim1D(&Ul_x1Face[j][i],&Wl[i] MHDARG( , &B1_x1Face[j][i]));
      Cons1D_to_Prim1D(&Ur_x1Face[j][i],&Wr[i] MHDARG( , &B1_x1Face[j][i]));

      fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],
                 MHDARG( B1_x1Face[j][i] , ) &x1Flux[j][i]);
    }
  }

/*--- Step 9c ------------------------------------------------------------------
 * Compute 2D x2-fluxes from corrected L/R states.
 */

  for (j=js; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta1[j-1][i],eta1[j][i]);
      etah = MAX(etah,eta1[j-1][i+1]);
      etah = MAX(etah,eta1[j  ][i+1]);
      etah = MAX(etah,eta2[j  ][i  ]);
#endif /* H_CORRECTION */
      Cons1D_to_Prim1D(&Ul_x2Face[j][i],&Wl[i] MHDARG( , &B2_x2Face[j][i]));
      Cons1D_to_Prim1D(&Ur_x2Face[j][i],&Wr[i] MHDARG( , &B2_x2Face[j][i]));

      fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[i],Wr[i],
                 MHDARG( B2_x2Face[j][i] , )&x2Flux[j][i]);
    }
  }

/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef MHD
  integrate_emf3_corner(pG);

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i] -= dtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      pG->B2i[ks][j][i] += dtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    pG->B1i[ks][j][ie+1] -= dtodx2*(emf3[j+1][ie+1] - emf3[j][ie+1]);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf3[je+1][i+1] - emf3[je+1][i]);
  }
#endif

/*=== STEP 11: Add source terms for a full timestep using n+1/2 states =======*/

/*--- Step 11a -----------------------------------------------------------------
 * Add gravitational (or shearing box) source terms as a Static Potential.
 * Remember with the 2D shearing box (1,2,3) = (x,z,y)
 *   A Crank-Nicholson update is used for shearing box terms.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;
  for(j=js; j<=je; ++j) {
    for(i=is; i<=ie; ++i) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);

/* Store the current state */
      M1n  = pG->U[ks][j][i].M1;
#ifdef FARGO
      dM3n = pG->U[ks][j][i].M3;
#else
      dM3n = pG->U[ks][j][i].M3 + qom*x1*pG->U[ks][j][i].d;
#endif

/* Calculate the flux for the y-momentum fluctuation (M3 in 2D) */
      frx1_dM3 = x1Flux[j][i+1].Mz;
      flx1_dM3 = x1Flux[j][i  ].Mz;
      frx2_dM3 = x2Flux[j+1][i].My;
      flx2_dM3 = x2Flux[j  ][i].My;
#ifndef FARGO
      frx1_dM3 += qom*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
      flx1_dM3 += qom*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
      frx2_dM3 += qom*(x1            )*x2Flux[j+1][i].d;
      flx2_dM3 += qom*(x1            )*x2Flux[j  ][i].d;
#endif

/* evolve M1n and dM3n by dt/2 using flux gradients */
      M1e = M1n - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
                - hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);

      dM3e = dM3n - hdtodx1*(frx1_dM3 - flx1_dM3)
                  - hdtodx2*(frx2_dM3 - flx2_dM3);

/* Update the 1- and 3-momentum (X and Y in 2D shearing box) for the Coriolis
 * and tidal potential source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

      pG->U[ks][j][i].M1 += (4.0*dM3e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
      pG->U[ks][j][i].M3 += 2.0*(qshear-2.)*(M1e + om_dt*dM3e)*fact;
#ifndef FARGO
      pG->U[ks][j][i].M3 -= 0.5*qshear*om_dt*(x1Flux[j][i].d+x1Flux[j][i+1].d);
#endif

/* Update the energy for a fixed potential, and add the Z-component (M2)
 * of the gravitational acceleration.
 * This update is identical to non-SHEARING_BOX below  */

      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic));
#endif

      phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
      phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

      pG->U[ks][j][i].M2 -= dtodx2*dhalf[j][i]*(phir-phil);

#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic));
#endif
    }
  }

#else /* ! SHEARING_BOX */

  if (StaticGravPot != NULL){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        pG->U[ks][j][i].M1 -= dtodx1*dhalf[j][i]*(phir-phil);

#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                     x1Flux[j][i+1].d*(phir - phic));
#endif
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[ks][j][i].M2 -= dtodx2*dhalf[j][i]*(phir-phil);

#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                     x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }

#endif /* SHEARING_BOX */

/*--- Step 11b -----------------------------------------------------------------
 * Add source terms for self-gravity.
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */

#ifdef SELF_GRAVITY
/* Add fluxes and source terms due to (d/dx1) terms  */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][j][i];
      phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j][i  ]);
      phir = 0.5*(pG->Phi[ks][j][i  ] + pG->Phi[ks][j][i+1]);

/* gx and gy centered at L and R x1-faces */
      gxl = (pG->Phi[ks][j][i-1] - pG->Phi[ks][j][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][j][i  ] - pG->Phi[ks][j][i+1])/(pG->dx1);

      gyl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j+1][i-1]) +
                  (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) )/(pG->dx2);
      gyr = 0.25*((pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) +
                  (pG->Phi[ks][j-1][i+1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx2);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
      flux_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
      flux_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

      flux_m2l = gxl*gyl/four_pi_G;
      flux_m2r = gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
      pG->U[ks][j][i].M1 -= dtodx1*(flux_m1r - flux_m1l);
      pG->U[ks][j][i].M2 -= dtodx1*(flux_m2r - flux_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic));
#endif
    }
  }

/* Add fluxes and source terms due to (d/dx2) terms  */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][j][i];
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j  ][i]);
      phir = 0.5*(pG->Phi[ks][j  ][i] + pG->Phi[ks][j+1][i]);

/* gx and gy centered at L and R x2-faces */
      gxl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j-1][i+1]) +
                  (pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) )/(pG->dx1);
      gxr = 0.25*((pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) +
                  (pG->Phi[ks][j+1][i-1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx1);

      gyl = (pG->Phi[ks][j-1][i] - pG->Phi[ks][j  ][i])/(pG->dx2);
      gyr = (pG->Phi[ks][j  ][i] - pG->Phi[ks][j+1][i])/(pG->dx2);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
      flux_m1l = gyl*gxl/four_pi_G;
      flux_m1r = gyr*gxr/four_pi_G;

      flux_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flux_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx2 terms  */
      pG->U[ks][j][i].M1 -= dtodx2*(flux_m1r - flux_m1l);
      pG->U[ks][j][i].M2 -= dtodx2*(flux_m2r - flux_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic));
#endif
    }
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pG->x1MassFlux[ks][j][i] = x1Flux[j][i].d;
      pG->x2MassFlux[ks][j][i] = x2Flux[j][i].d;
    }
  }
#endif /* SELF_GRAVITY */

/*--- Step 11c -----------------------------------------------------------------
 * Add source terms for optically thin cooling
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        coolf = (*CoolingFunc)(dhalf[j][i],phalf[j][i],pG->dt);
        pG->U[ks][j][i].E -= pG->dt*coolf;
      }
    }
  }
#endif /* BAROTROPIC */

/*=== STEP 12: Update cell-centered values for a full timestep ===============*/

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x1-fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d  -= dtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      pG->U[ks][j][i].M2 -= dtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      pG->U[ks][j][i].M3 -= dtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B2c -= dtodx1*(x1Flux[j][i+1].By - x1Flux[j][i].By);
      pG->U[ks][j][i].B3c -= dtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx1*(x1Flux[j][i+1].s[n] 
                                         - x1Flux[j][i  ].s[n]);
#endif

    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x2-fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d  -= dtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      pG->U[ks][j][i].M2 -= dtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      pG->U[ks][j][i].M3 -= dtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B3c -= dtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
      pG->U[ks][j][i].B1c -= dtodx2*(x2Flux[j+1][i].Bz - x2Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx2*(x2Flux[j+1][i].s[n] 
                                         - x2Flux[j  ][i].s[n]);
#endif
    }
  }

/*--- Step 12c: Not needed in 2D ---*/
/*--- Step 12d -----------------------------------------------------------------
 * LAST STEP!
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef MHD
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {

      pG->U[ks][j][i].B1c =0.5*(pG->B1i[ks][j][i]+pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c =0.5*(pG->B2i[ks][j][i]+pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
      pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
    }
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct_2d:  Free temporary integration arrays */

void integrate_destruct_2d(void)
{
#ifdef MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_2d_array(eta1);
  if (eta2 != NULL) free_2d_array(eta2);
#endif /* H_CORRECTION */
#ifdef MHD
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_2d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_2d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_2d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_2d_array(Ur_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (dhalf     != NULL) free_2d_array(dhalf);
  if (phalf     != NULL) free_2d_array(phalf);

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_init_2d:    Allocate temporary integration arrays */

void integrate_init_2d(int nx1, int nx2)
{
  int nmax;
  int Nx1 = nx1 + 2*nghost;
  int Nx2 = nx2 + 2*nghost;
  nmax = MAX(Nx1,Nx2);

#ifdef MHD
  if ((emf3 = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((eta2 = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* H_CORRECTION */

#ifdef MHD
  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

  if ((B1_x1Face = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((W  =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;

  if ((Ul_x1Face = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;
  if ((Ur_x1Face = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;
  if ((Ul_x2Face = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;
  if ((Ur_x2Face = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((x1Flux    = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;
  if ((x2Flux    = (Cons1D**)calloc_2d_array(Nx2, Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
  {
  if ((dhalf = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  }

  return;

  on_error:
  integrate_destruct();
  ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf3_corner:  */

#ifdef MHD
static void integrate_emf3_corner(Grid *pG)
{
  int i,is,ie,j,js,je;
  Real emf_l1, emf_r1, emf_l2, emf_r2;

  is = pG->is;   ie = pG->ie;
  js = pG->js;   je = pG->je;

/* NOTE: The x1-Flux of B2 is -E3.  The x2-Flux of B1 is +E3. */
  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
      if (x1Flux[j-1][i].d > 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1]);
      }
      else if (x1Flux[j-1][i].d < 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i].Bz - emf3_cc[j-1][i]);

      } else {
	emf_l2 = -x1Flux[j-1][i].By
	  + 0.5*(x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1] + 
		 x2Flux[j][i  ].Bz - emf3_cc[j-1][i  ] );
      }

      if (x1Flux[j][i].d > 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i-1].Bz - emf3_cc[j][i-1]);
      }
      else if (x1Flux[j][i].d < 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i].Bz - emf3_cc[j][i]);

      } else {
	emf_r2 = -x1Flux[j][i].By 
	  + 0.5*(x2Flux[j][i-1].Bz - emf3_cc[j][i-1] + 
		 x2Flux[j][i  ].Bz - emf3_cc[j][i  ] );
      }

      if (x2Flux[j][i-1].d > 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]);
      }
      else if (x2Flux[j][i-1].d < 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz 
          + (-x1Flux[j][i].By - emf3_cc[j][i-1]);
      } else {
	emf_l1 = x2Flux[j][i-1].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i-1] );
      }

      if (x2Flux[j][i].d > 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i]);
      }
      else if (x2Flux[j][i].d < 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j][i].By - emf3_cc[j][i]);
      } else {
	emf_r1 = x2Flux[j][i].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i] );
      }

      emf3[j][i] = 0.25*(emf_l1 + emf_r1 + emf_l2 + emf_r2);
    }
  }

  return;
}
#endif /* MHD */

#endif /* CTU_INTEGRATOR */
