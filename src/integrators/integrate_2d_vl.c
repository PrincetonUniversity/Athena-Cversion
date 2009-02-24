#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_2d_vl.c
 *
 * PURPOSE: Integrate MHD equations using 2D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type Gas
 *      B1i, B2i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *
 * REFERENCE: J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensional dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_2d
 *   integrate_destruct_2d()
 *   integrate_init_2d()
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

#ifdef VL_INTEGRATOR

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1D **Wl_x1Face=NULL, **Wr_x1Face=NULL;
static Prim1D **Wl_x2Face=NULL, **Wr_x2Face=NULL;
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

/* conserved variables at t^{n+1/2} computed in predict step */
static Gas **Uhalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method of Gardiner & Stone (2005) 
 *============================================================================*/

#ifdef MHD
static void integrate_emf3_corner(const Grid *pG);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_2d: van Leer unsplit integrator in 2D. 
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */

void integrate_2d(Grid *pG, Domain *pD)
{
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  Real dt = pG->dt, hdt = 0.5*pG->dt;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
#ifdef MHD
  Real d, M1, M2, M3, B1c, B2c, B3c;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flx_m1l,flx_m1r,flx_m2l,flx_m2r;
#endif

  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      Uhalf[j][i] = pG->U[ks][j][i];
    }
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-(nghost-1); i<=ie+nghost; i++) {
      Ul[i].d  = pG->U[ks][j][i-1].d;
      Ul[i].Mx = pG->U[ks][j][i-1].M1;
      Ul[i].My = pG->U[ks][j][i-1].M2;
      Ul[i].Mz = pG->U[ks][j][i-1].M3;
#ifndef BAROTROPIC
      Ul[i].E  = pG->U[ks][j][i-1].E;
#endif /* BAROTROPIC */
#ifdef MHD
      B1_x1Face[j][i] = pG->B1i[ks][j][i];
      Ul[i].By = pG->U[ks][j][i-1].B2c;
      Ul[i].Bz = pG->U[ks][j][i-1].B3c;
      Bxc[i] = pG->U[ks][j][i-1].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) Ul[i].s[n] = pG->U[ks][j][i-1].s[n];
#endif

      Cons1D_to_Prim1D(&Ul[i],&Wl[i] MHDARG( , &Bxc[i]));

      Ur[i].d  = pG->U[ks][j][i].d;
      Ur[i].Mx = pG->U[ks][j][i].M1;
      Ur[i].My = pG->U[ks][j][i].M2;
      Ur[i].Mz = pG->U[ks][j][i].M3;
#ifndef BAROTROPIC
      Ur[i].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      Ur[i].By = pG->U[ks][j][i].B2c;
      Ur[i].Bz = pG->U[ks][j][i].B3c;
      Bxc[i] = pG->U[ks][j][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) Ur[i].s[n] = pG->U[ks][j][i  ].s[n];
#endif

      Cons1D_to_Prim1D(&Ur[i],&Wr[i] MHDARG( , &Bxc[i]));

/* Compute flux in x1-direction  */
      fluxes(Ul[i],Ur[i],Wl[i],Wr[i],MHDARG(B1_x1Face[j][i] ,) &x1Flux[j][i]);
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

  for (i=is-nghost; i<=ie+nghost; i++) {
    for (j=js-(nghost-1); j<=je+nghost; j++) {
      Ul[j].d  = pG->U[ks][j-1][i].d;
      Ul[j].Mx = pG->U[ks][j-1][i].M2;
      Ul[j].My = pG->U[ks][j-1][i].M3;
      Ul[j].Mz = pG->U[ks][j-1][i].M1;
#ifndef BAROTROPIC
      Ul[j].E  = pG->U[ks][j-1][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      B2_x2Face[j][i] = pG->B2i[ks][j][i];
      Ul[j].By = pG->U[ks][j-1][i].B3c;
      Ul[j].Bz = pG->U[ks][j-1][i].B1c;
      Bxc[j] = pG->U[ks][j-1][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) Ul[j].s[n] = pG->U[ks][j-1][i].s[n];
#endif

      Cons1D_to_Prim1D(&Ul[j],&Wl[j] MHDARG( , &B2_x2Face[j][i]));

      Ur[j].d  = pG->U[ks][j][i].d;
      Ur[j].Mx = pG->U[ks][j][i].M2;
      Ur[j].My = pG->U[ks][j][i].M3;
      Ur[j].Mz = pG->U[ks][j][i].M1;
#ifndef BAROTROPIC
      Ur[j].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      Ur[j].By = pG->U[ks][j][i].B3c;
      Ur[j].Bz = pG->U[ks][j][i].B1c;
      Bxc[j] = pG->U[ks][j][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) Ur[j].s[n] = pG->U[ks][j][i].s[n];
#endif

      Cons1D_to_Prim1D(&Ur[j],&Wr[j] MHDARG( , &B2_x2Face[j][i]));

/* Compute flux in x2-direction */
      fluxes(Ul[j],Ur[j],Wl[j],Wr[j],MHDARG(B2_x2Face[j][i] ,) &x2Flux[j][i]);
    }
  }

/*=== STEP 3: Not needed in 2D ===*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (j=js-(nghost-1); j<=je+(nghost-1); j++) {
    for (i=is-(nghost-1); i<=ie+(nghost-1); i++) {
      emf3_cc[j][i] =
        (pG->U[ks][j][i].B1c*pG->U[ks][j][i].M2 -
         pG->U[ks][j][i].B2c*pG->U[ks][j][i].M1)/pG->U[ks][j][i].d;
    }
  }
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    B1_x1Face[j][iu] -= hdtodx2*(emf3[j+1][iu]-emf3[j][iu]);
  }
  for (i=il+1; i<=iu-1; i++) {
    B2_x2Face[ju][i] += hdtodx1*(emf3[ju][i+1]-emf3[ju][i]);
  }

/*--- Step 4c ------------------------------------------------------------------
 * Compute cell-centered magnetic fields at half-timestep from average of
 * face-centered fields.
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      Uhalf[j][i].B1c = 0.5*(B1_x1Face[j][i] + B1_x1Face[j][i+1]);
      Uhalf[j][i].B2c = 0.5*(B2_x2Face[j][i] + B2_x2Face[j+1][i]);
    }
  }
#endif /* MHD */

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x1-fluxes
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      Uhalf[j][i].d  -= hdtodx1*(x1Flux[j][i+1].d  - x1Flux[j][i].d );
      Uhalf[j][i].M1 -= hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx);
      Uhalf[j][i].M2 -= hdtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My);
      Uhalf[j][i].M3 -= hdtodx1*(x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      Uhalf[j][i].E  -= hdtodx1*(x1Flux[j][i+1].E  - x1Flux[j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Uhalf[j][i].s[n] -= hdtodx1*(x1Flux[j][i+1].s[n] - x1Flux[j][i].s[n]);
#endif
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x2-fluxes
 */

  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      Uhalf[j][i].d  -= hdtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      Uhalf[j][i].M1 -= hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      Uhalf[j][i].M2 -= hdtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      Uhalf[j][i].M3 -= hdtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      Uhalf[j][i].E  -= hdtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Uhalf[j][i].s[n] -= hdtodx2*(x2Flux[j+1][i].s[n] - x2Flux[j][i].s[n]);
#endif
    }
  }

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)( x1,             x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        Uhalf[j][i].M1 -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Uhalf[j][i].E -= hdtodx1*(x1Flux[j][i  ].d*(phic - phil)
                           + x1Flux[j][i+1].d*(phir - phic));
#endif

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Uhalf[j][i].M2 -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
        Uhalf[j][i].E -= hdtodx2*(x2Flux[j  ][i].d*(phic - phil)
                                + x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 6b ------------------------------------------------------------------
 * Add source terms for self gravity for 0.5*dt to predict step.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);

      Uhalf[j][i].M1 -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Uhalf[j][i].E -= hdtodx1*(x1Flux[j][i  ].d*(phic - phil)
                              + x1Flux[j][i+1].d*(phir - phic));
#endif

      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);

      Uhalf[j][i].M2 -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
      Uhalf[j][i].E -= hdtodx2*(x2Flux[j  ][i].d*(phic - phil)
                              + x2Flux[j+1][i].d*(phir - phic));
#endif
    }
  }
#endif /* SELF_GRAVITY */

/*=== STEP 7: Compute second-order L/R x1-interface states ===================*/

/*--- Step 7a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      U1d[i].d  = Uhalf[j][i].d;
      U1d[i].Mx = Uhalf[j][i].M1;
      U1d[i].My = Uhalf[j][i].M2;
      U1d[i].Mz = Uhalf[j][i].M3;
#ifndef BAROTROPIC
      U1d[i].E  = Uhalf[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[i].By = Uhalf[j][i].B2c;
      U1d[i].Bz = Uhalf[j][i].B3c;
      Bxc[i] = Uhalf[j][i].B1c;
      Bxi[i] = B1_x1Face[j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[i].s[n] = Uhalf[j][i].s[n];
#endif

      Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));
    }

/*--- Step 7b ------------------------------------------------------------------
 * Compute L and R states at x1-interfaces, store in 2D array
 */

    lr_states(W, MHDARG( Bxc , ) 0.0,0.0,il+1,iu-1,Wl,Wr);

    for (i=il+1; i<=iu; i++) {
      Wl_x1Face[j][i] = Wl[i];
      Wr_x1Face[j][i] = Wr[i];
    }
  }

/*=== STEP 8: Compute second-order L/R x2-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (i=il; i<=iu; i++) {
    for (j=jl; j<=ju; j++) {
      U1d[j].d  = Uhalf[j][i].d;
      U1d[j].Mx = Uhalf[j][i].M2;
      U1d[j].My = Uhalf[j][i].M3;
      U1d[j].Mz = Uhalf[j][i].M1;
#ifndef BAROTROPIC
      U1d[j].E  = Uhalf[j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[j].By = Uhalf[j][i].B3c;
      U1d[j].Bz = Uhalf[j][i].B1c;
      Bxc[j] = Uhalf[j][i].B2c;
      Bxi[j] = B2_x2Face[j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = Uhalf[j][i].s[n];
#endif
      Cons1D_to_Prim1D(&U1d[j],&W[j] MHDARG( , &Bxc[j]));
    }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L and R states at x2-interfaces, store in 2D array
 */

    lr_states(W, MHDARG( Bxc , ) 0.0,0.0,jl+1,ju-1,Wl,Wr);

    for (j=jl+1; j<=ju; j++) {
      Wl_x2Face[j][i] = Wl[j];
      Wr_x2Face[j][i] = Wr[j];
    }
  }

/*=== STEP 9: Not needed in 2D ===*/

/*=== STEP 10: Compute 2D x1-Flux, x2-Flux, x3-Flux ==========================*/

/*--- Step 10a -----------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction, if needed.
 */

#ifdef H_CORRECTION
  for (j=js-1; j<=je+1; j++) {
    for (i=is-11; i<=iu; i++) {
      cfr = cfast(&(Ur_x1Face[j][i]) MHDARG( , &(B1_x1Face[j][i])));
      cfl = cfast(&(Ul_x1Face[j][i]) MHDARG( , &(B1_x1Face[j][i])));
      lambdar = Ur_x1Face[j][i].Mx/Ur_x1Face[j][i].d + cfr;
      lambdal = Ul_x1Face[j][i].Mx/Ul_x1Face[j][i].d - cfl;
      eta1[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }

  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
      cfr = cfast(&(Ur_x2Face[j][i]) MHDARG( , &(B2_x2Face[j][i])));
      cfl = cfast(&(Ul_x2Face[j][i]) MHDARG( , &(B2_x2Face[j][i])));
      lambdar = Ur_x2Face[j][i].Mx/Ur_x2Face[j][i].d + cfr;
      lambdal = Ul_x2Face[j][i].Mx/Ul_x2Face[j][i].d - cfl;
      eta2[j][i] = 0.5*fabs(lambdar - lambdal);
    }
  }
#endif /* H_CORRECTION */

/*--- Step 10b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta2[j][i-1],eta2[j][i]);
      etah = MAX(etah,eta2[j+1][i-1]);
      etah = MAX(etah,eta2[j+1][i  ]);
      etah = MAX(etah,eta1[j  ][i  ]);
#endif /* H_CORRECTION */

      Prim1D_to_Cons1D(&Ul[i],&Wl_x1Face[j][i] MHDARG( , &B1_x1Face[j][i]));
      Prim1D_to_Cons1D(&Ur[i],&Wr_x1Face[j][i] MHDARG( , &B1_x1Face[j][i]));
      fluxes(Ul[i],Ur[i],Wl_x1Face[j][i],Wr_x1Face[j][i],
                 MHDARG( B1_x1Face[j][i] , ) &x1Flux[j][i]);
    }
  }

/*--- Step 10c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

  for (j=js; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
      etah = MAX(eta1[j-1][i],eta1[j][i]);
      etah = MAX(etah,eta1[j-1][i+1]);
      etah = MAX(etah,eta1[j  ][i+1]);
      etah = MAX(etah,eta2[j  ][i]);
#endif /* H_CORRECTION */
      Prim1D_to_Cons1D(&Ul[i],&Wl_x2Face[j][i] MHDARG( , &B2_x2Face[j][i]));
      Prim1D_to_Cons1D(&Ur[i],&Wr_x2Face[j][i] MHDARG( , &B2_x2Face[j][i]));
      fluxes(Ul[i],Ur[i],Wl_x2Face[j][i],Wr_x2Face[j][i],
                 MHDARG( B2_x2Face[j][i] , ) &x2Flux[j][i]);
    }
  }

/*=== STEP 11: Update face-centered B for a full timestep ====================*/
        
/*--- Step 11a -----------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at the half-time-step.
 */

#ifdef MHD
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      d  = Uhalf[j][i].d ;
      M1 = Uhalf[j][i].M1;
      M2 = Uhalf[j][i].M2;
      B1c = Uhalf[j][i].B1c;
      B2c = Uhalf[j][i].B2c;

      emf3_cc[j][i] = (B1c*M2 - B2c*M1)/d;
    }
  }
#endif

/*--- Step 11b -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef MHD
  integrate_emf3_corner(pG);

/*--- Step 11c -----------------------------------------------------------------
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

/*=== STEP 12: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 12a -----------------------------------------------------------------
 * Add gravitational source terms due to a Static Potential
 * To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)( x1,             x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        pG->U[ks][j][i].M1 -= dtodx1*(phir-phil)*Uhalf[j][i].d;
#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil)
                                   + x1Flux[j][i+1].d*(phir - phic));
#endif

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[ks][j][i].M2 -= dtodx2*(phir-phil)*Uhalf[j][i].d;
#ifndef BAROTROPIC
        pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil)
                                   + x2Flux[j+1][i].d*(phir - phic));
#endif
      }
    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Add gravitational source terms for self-gravity.
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

/* gx, gy and gz centered at L and R x1-faces */
      gxl = (pG->Phi[ks][j][i-1] - pG->Phi[ks][j][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][j][i  ] - pG->Phi[ks][j][i+1])/(pG->dx1);

      gyl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j+1][i-1]) +
                  (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]))/(pG->dx2);
      gyr = 0.25*((pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) +
                  (pG->Phi[ks][j-1][i+1] - pG->Phi[ks][j+1][i+1]))/(pG->dx2);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
      flx_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
      flx_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

      flx_m2l = gxl*gyl/four_pi_G;
      flx_m2r = gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
      pG->U[ks][j][i].M1 -= dtodx1*(flx_m1r - flx_m1l);
      pG->U[ks][j][i].M2 -= dtodx1*(flx_m2r - flx_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic));
#endif /* BAROTROPIC */
    }
  }

/* Add fluxes and source terms due to (d/dx2) terms  */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      phic = pG->Phi[ks][j][i];
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j  ][i]);
      phir = 0.5*(pG->Phi[ks][j  ][i] + pG->Phi[ks][j+1][i]);

/* gx, gy and gz centered at L and R x2-faces */
      gxl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j-1][i+1]) +
                  (pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]))/(pG->dx1);
      gxr = 0.25*((pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) +
                  (pG->Phi[ks][j+1][i-1] - pG->Phi[ks][j+1][i+1]))/(pG->dx1);

      gyl = (pG->Phi[ks][j-1][i] - pG->Phi[ks][j  ][i])/(pG->dx2);
      gyr = (pG->Phi[ks][j  ][i] - pG->Phi[ks][j+1][i])/(pG->dx2);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
      flx_m1l = gyl*gxl/four_pi_G;
      flx_m1r = gyr*gxr/four_pi_G;

      flx_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flx_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx2 terms  */
      pG->U[ks][j][i].M1 -= dtodx2*(flx_m1r - flx_m1l);
      pG->U[ks][j][i].M2 -= dtodx2*(flx_m2r - flx_m2l);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic));
#endif /* BAROTROPIC */
    }
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pG->x1MassFlux[j][i] = x1Flux[j][i].d;
      pG->x2MassFlux[j][i] = x2Flux[j][i].d;
      pG->x3MassFlux[j][i] = x3Flux[j][i].d;
    }
  }
#endif /* SELF_GRAVITY */

/*=== STEP 13: Update cell-centered values for a full timestep ===============*/

/*--- Step 13a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x1-Fluxes
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
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx1*(x1Flux[j][i+1].s[n]
                                      - x1Flux[j][i  ].s[n]);
#endif
    }
  }

/*--- Step 13b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 2D x2-Fluxes
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
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx2*(x2Flux[j+1][i].s[n]
                                      - x2Flux[j  ][i].s[n]);
#endif
    }
  }

/*--- Step 13e -----------------------------------------------------------------
 * LAST STEP!
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef MHD
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i]+pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i]+pG->B2i[ks][j+1][i]);
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

  if (Wl_x1Face != NULL) free_2d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_2d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_2d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_2d_array(Wr_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);

  if (Uhalf    != NULL) free_2d_array(Uhalf);

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_init_2d: Allocate temporary integration arrays */

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
  if ((B1_x1Face = (Real**)calloc_2d_array(Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((W   =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wl  =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;
  if ((Wr  =      (Prim1D*)malloc(nmax*sizeof(Prim1D))) == NULL) goto on_error;

  if ((Wl_x1Face = (Prim1D**)calloc_2d_array(Nx2,Nx1, sizeof(Prim1D)))
    == NULL) goto on_error;
  if ((Wr_x1Face = (Prim1D**)calloc_2d_array(Nx2,Nx1, sizeof(Prim1D)))
    == NULL) goto on_error;
  if ((Wl_x2Face = (Prim1D**)calloc_2d_array(Nx2,Nx1, sizeof(Prim1D)))
    == NULL) goto on_error;
  if ((Wr_x2Face = (Prim1D**)calloc_2d_array(Nx2,Nx1, sizeof(Prim1D)))
    == NULL) goto on_error;
  if ((x1Flux    = (Cons1D**)calloc_2d_array(Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;
  if ((x2Flux    = (Cons1D**)calloc_2d_array(Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;

  if ((Uhalf = (Gas**)calloc_2d_array(Nx2,Nx1, sizeof(Gas))) == NULL)
    goto on_error;

  return;

  on_error:
  integrate_destruct();
  ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf3_corner()
 *   Integrates face centered B-fluxes to compute corner EMFs.  Note:
 *   x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *   x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *   x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *   x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 */ 

#ifdef MHD
static void integrate_emf3_corner(const Grid *pG)
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

#endif /* VL_INTEGRATOR */
