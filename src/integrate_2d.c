#include "copyright.h"
/*==============================================================================
 * FILE: integrate_2d.c
 *
 * PURPOSE: Updates the input Grid structure pointed to by *pG by one 
 *   timestep using directionally unsplit CTU method of Colella (1990).  The
 *   variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type Gas
 *      B1i, B2i -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and H-correction
 *   of Sanders et al.
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
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

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

/* density at t^{n+1/2} needed by both MHD and to make gravity 2nd order */
static Real **dhalf = NULL;

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
/* integrate_2d:  CTU integrator in 2D  */

void integrate_2d(Grid *pG, Domain *pD)
{
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,il,iu;
  int j,jl,ju;
#ifdef MHD
  Real MHD_src,dbx,dby,B1,B2,B3,V3;
  Real d, M1, M2, B1c, B2c;
  Real hdt = 0.5*pG->dt;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flux_m1l,flux_m1r,flux_m2l,flux_m2r;
#endif
#ifdef SHEARING_BOX
  Real M1n, dM3n;   /* M1, dM3=(My+d*1.5*Omega*x) at time n */
  Real M1e, dM3e;   /* M1, dM3 evolved by dt/2  */
  Real flx1_dM3, frx1_dM3, flx2_dM3, frx2_dM3;
  Real fact, TH_om, om_dt = Omega*pG->dt;
#endif /* SHEARING_BOX */

  il = is - 2;
  iu = ie + 2;

  jl = js - 2;
  ju = je + 2;

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
 * Compute L and R states at X1-interfaces.
 */

    for (i=is-nghost; i<=ie+nghost; i++) {
      Cons1D_to_Prim1D(&U1d[i],&W[i] MHDARG( , &Bxc[i]));
    }
    lr_states(W, MHDARG( Bxc , ) pG->dt,dtodx1,is-1,ie+1,Wl,Wr);

/* Add "MHD source terms" for 0.5*dt */

#ifdef MHD
    for (i=is-1; i<=iu; i++) {
      MHD_src = (pG->U[ks][j][i-1].M2/pG->U[ks][j][i-1].d)*
               (pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1])/pG->dx1;
      Wl[i].By += hdt*MHD_src;

      MHD_src = (pG->U[ks][j][i].M2/pG->U[ks][j][i].d)*
               (pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i])/pG->dx1;
      Wr[i].By += hdt*MHD_src;
    }
#endif

/*--- Step 1c ------------------------------------------------------------------
 * Add gravitational source terms from static potential for 0.5*dt to L/R states
 */

    if (StaticGravPot != NULL){
      for (i=is-1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*StaticGravPot)( x1             ,x2,x3);
        phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
        phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* Apply gravitational source terms to velocity using gradient of potential.
 * for (dt/2).   S_{V} = -Grad(Phi) */
        Wl[i].Vx -= dtodx1*(phifc - phicl);
        Wr[i].Vx -= dtodx1*(phicr - phifc);
      }
    }

/*--- Step 1d ------------------------------------------------------------------
 * Add gravitational source terms for self-gravity for dt/2 to L/R states
 */

#ifdef SELF_GRAVITY
    for (i=is-1; i<=iu; i++) {
      Wl[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      Wr[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
    }
#endif

/*--- Step 1d (cont) -----------------------------------------------------------
 * Shearing box source terms (Coriolis forces) in 2D X-Z plane.
 *  (x1,x2,x3) in code = (X,Z,Y) in shearing sheet
 */

#ifdef SHEARING_BOX
      for (i=is-1; i<=iu; i++) {
        Wl[i].Vx += pG->dt*Omega*W[i-1].Vz; /* (dt/2)*( 2 Omega Vy) */
        Wl[i].Vz -= pG->dt*Omega*W[i-1].Vx; /* (dt/2)*(-2 Omega Vx) */

        Wr[i].Vx += pG->dt*Omega*W[i].Vz; /* (dt/2)*( 2 Omega Vy) */
        Wr[i].Vz -= pG->dt*Omega*W[i].Vx; /* (dt/2)*(-2 Omega Vx) */
    }
#endif /* SHEARING_BOX */

/*--- Step 1e ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 2D array
 */

    for (i=is-1; i<=iu; i++) {
      Prim1D_to_Cons1D(&Ul_x1Face[j][i],&Wl[i] MHDARG( , &Bxi[i]));
      Prim1D_to_Cons1D(&Ur_x1Face[j][i],&Wr[i] MHDARG( , &Bxi[i]));

      GET_FLUXES(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],
                 MHDARG( B1_x1Face[j][i] , ) &x1Flux[j][i]);
    }
  }

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
 * Compute L and R states at X2-interfaces.
 */

    for (j=js-nghost; j<=je+nghost; j++) {
      Cons1D_to_Prim1D(&U1d[j],&W[j] MHDARG( , &Bxc[j]));
    }
    lr_states(W, MHDARG( Bxc , ) pG->dt,dtodx2,js-1,je+1,Wl,Wr);

/* Add "MHD source terms" for 0.5*dt */

#ifdef MHD
    for (j=js-1; j<=ju; j++) {
      MHD_src = (pG->U[ks][j-1][i].M1/pG->U[ks][j-1][i].d)*
        (pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i])/pG->dx2;
      Wl[j].Bz += hdt*MHD_src;

      MHD_src = (pG->U[ks][j][i].M1/pG->U[ks][j][i].d)*
        (pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i])/pG->dx2;
      Wr[j].Bz += hdt*MHD_src;
    }
#endif

/*--- Step 2c ------------------------------------------------------------------
 * Add gravitational source terms from static potential for 0.5*dt to L/R states
 */
  
    if (StaticGravPot != NULL){
      for (j=js-1; j<=ju; j++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*StaticGravPot)(x1, x2             ,x3);
        phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
        phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

/* Apply gravitational source terms to velocity using gradient of potential.
 * for (dt/2).   S_{V} = -Grad(Phi) */
        Wl[j].Vx -= dtodx2*(phifc - phicl);
        Wr[j].Vx -= dtodx2*(phicr - phifc);
      }
    }

/*--- Step 2d ------------------------------------------------------------------
 * Add gravitational source terms for self-gravity for dt/2 to L/R states
 */

#ifdef SELF_GRAVITY
      for (j=js-1; j<=ju; j++) {
        Wl[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
        Wr[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      }
#endif

/*--- Step 2e ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 2D array
 */

    for (j=js-1; j<=ju; j++) {
      Prim1D_to_Cons1D(&Ul_x2Face[j][i],&Wl[j] MHDARG( , &Bxi[j]));
      Prim1D_to_Cons1D(&Ur_x2Face[j][i],&Wr[j] MHDARG( , &Bxi[j]));

      GET_FLUXES(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],
                 MHDARG( B2_x2Face[j][i] , ) &x2Flux[j][i]);
    }
  }

/*--- Step 3 ------------------------------------------------------------------
 * Calculate the cell centered value of emf_3 at t^{n}
 */

#ifdef MHD
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      emf3_cc[j][i] =
	(pG->U[ks][j][i].B1c*pG->U[ks][j][i].M2 -
	 pG->U[ks][j][i].B2c*pG->U[ks][j][i].M1 )/pG->U[ks][j][i].d;
    }
  }

/*--- Step 4 ------------------------------------------------------------------
 * Integrate emf3 to the grid cell corners and then update the 
 * interface magnetic fields using CT for a half time step.
 */

  integrate_emf3_corner(pG);

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    B1_x1Face[j][iu] -= hdtodx2*(emf3[j+1][iu] - emf3[j][iu]);
  }
  for (i=is-1; i<=ie+1; i++) {
    B2_x2Face[ju][i] += hdtodx1*(emf3[ju][i+1] - emf3[ju][i]);
  }
#endif

/*--- Step 5a ------------------------------------------------------------------
 * Correct the L/R states at x1-interfaces using transverse flux-gradients in
 * the x2-direction for 0.5*dt using x2-fluxes computed in Step 2e.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=iu; i++) {
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

/*--- Step 5b ------------------------------------------------------------------
 * Add the "MHD source terms" to the x2 (conservative) flux gradient.
 */

#ifdef MHD
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=iu; i++) {
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

/*--- Step 5c ------------------------------------------------------------------
 * Add source terms for a Static Gravitational Potential to L/R states.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=iu; i++) {
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

/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for self gravity to L/R states.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=iu; i++) {
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

/*--- Step 6a ------------------------------------------------------------------
 * Correct the L/R states at x2-interfaces using transverse flux-gradients in
 * the x1-direction for 0.5*dt using x1-fluxes computed in Step 1e.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (y,z,x) on LHS */

  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
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

/*--- Step 6b ------------------------------------------------------------------
 * Add the "MHD source terms" to the x1 (conservative) flux gradient.
 */

#ifdef MHD
  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
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

/*--- Step 6c ------------------------------------------------------------------
 * Add source terms for a Static Gravitational Potential to L/R states.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=js-1; j<=ju; j++) {
      for (i=is-1; i<=ie+1; i++) {
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

/*--- Step 6d ------------------------------------------------------------------
 * Add source terms for self gravity to L/R states.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
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

/*--- Step 6e ------------------------------------------------------------------
 * Add the tidal potential and Coriolis terms in 2D X-Z plane.
 *  (x1,x2,x3) in code = (X,Z,Y) in shearing sheet
 *    Vx source term is (dt/2)( 2 Omega V y); Mx on x2Face is Mz 
 *    Vy source term is (dt/2)(-2 Omega V x); My on x2Face is Mx
 */
#ifdef SHEARING_BOX
  for (j=js-1; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
      Ur_x2Face[j][i].Mz += pG->dt*Omega*pG->U[ks][j][i].M3;
      Ur_x2Face[j][i].Mx -= pG->dt*Omega*pG->U[ks][j][i].M1;

      Ul_x2Face[j][i].Mz += pG->dt*Omega*pG->U[ks][j-1][i].M3;
      Ul_x2Face[j][i].Mx -= pG->dt*Omega*pG->U[ks][j-1][i].M1;
    }
  }
#endif /* SHEARING_BOX */

/*--- Step 7 ------------------------------------------------------------------
 * Calculate the cell centered value of emf_3 at t^{n+1/2}, needed by CT 
 * algorithm to integrate emf to corner in step 10
 */

  if (dhalf != NULL){
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        dhalf[j][i] = pG->U[ks][j][i].d
          - hdtodx1*(x1Flux[j  ][i+1].d - x1Flux[j][i].d)
          - hdtodx2*(x2Flux[j+1][i  ].d - x2Flux[j][i].d);
      }
    }
  }

#ifdef MHD
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);

      d  = dhalf[j][i];

      M1 = pG->U[ks][j][i].M1
        - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
        - hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);

      M2 = pG->U[ks][j][i].M2
        - hdtodx1*(x1Flux[j][i+1].My - x1Flux[j][i].My)
        - hdtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);

/* Add source terms for fixed gravitational potential */
      if (StaticGravPot != NULL){
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
        M1 -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
        M2 -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
      }

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);
      M1 -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);
      M2 -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
#endif /* SELF_GRAVITY */

/* Add the tidal potential and Coriolis terms for shearing box. */
#ifdef SHEARING_BOX
      M1 += pG->dt*Omega*pG->U[ks][j][i].M3;
#endif /* SHEARING_BOX */

      B1c = 0.5*(B1_x1Face[j][i] + B1_x1Face[j][i+1]);
      B2c = 0.5*(B2_x2Face[j][i] + B2_x2Face[j+1][i]);

      emf3_cc[j][i] = (B1c*M2 - B2c*M1)/d;
    }
  }
#endif /* MHD */

/*--- Step 8a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

#ifdef H_CORRECTION
  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=iu; i++) {
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

/*--- Step 8b ------------------------------------------------------------------
 * Compute x1-fluxes from corrected L/R states.
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

      GET_FLUXES(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],
                 MHDARG( B1_x1Face[j][i] , ) &x1Flux[j][i]);
    }
  }

/*--- Step 8c ------------------------------------------------------------------
 * Compute x2-fluxes from corrected L/R states.
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

      GET_FLUXES(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[i],Wr[i],
                 MHDARG( B2_x2Face[j][i] , )&x2Flux[j][i]);
    }
  }

/*--- Step 9 -------------------------------------------------------------------
 * Integrate emf3^{n+1/2} to the grid cell corners and then update the 
 * interface magnetic fields using CT for a full time step.
 */

#ifdef MHD
  integrate_emf3_corner(pG);

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


/*--- Step 10a -----------------------------------------------------------------
 * Add the gravitational (or shearing box) source terms expressed as a Static
 * Potential.  Remember with the 2D shearing box (1,2,3) = (x,z,y)
 *   A Crank-Nicholson update is used for shearing box terms.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

#ifdef SHEARING_BOX
  fact = om_dt/(1.0 + 0.25*om_dt*om_dt);
  TH_om = 1.5*Omega; /* Three-Halves Omega */
  for(j=js; j<=je; ++j) {
    for(i=is; i<=ie; ++i) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);

/* Store the current state */
      M1n  = pG->U[ks][j][i].M1;
      dM3n = pG->U[ks][j][i].M3 + pG->U[ks][j][i].d*TH_om*x1;

/* Calculate the flux for the y-momentum fluctuation (M3 in 2D) */
      frx1_dM3 = x1Flux[j][i+1].Mz + TH_om*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
      flx1_dM3 = x1Flux[j][i  ].Mz + TH_om*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
      frx2_dM3 = x2Flux[j+1][i].My + TH_om*(x1            )*x2Flux[j+1][i].d;
      flx2_dM3 = x2Flux[j  ][i].My + TH_om*(x1            )*x2Flux[j  ][i].d;

/* evolve M1n and dM3n by dt/2 using Forward Euler */
      M1e = M1n - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
                - hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);

      dM3e = dM3n - hdtodx1*(frx1_dM3 - flx1_dM3)
                  - hdtodx2*(frx2_dM3 - flx2_dM3);

/* Update the 1- and 3-momentum (X and Y in 2D shearing box) for the Coriolis
 * and tidal potential momentum source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

      pG->U[ks][j][i].M1 += (2.0*dM3e - 0.5*om_dt*M1e)*fact;
      pG->U[ks][j][i].M3 -= (0.5*(M1e + om_dt*dM3e)*fact + 
           0.75*om_dt*(x1Flux[j][i].d + x1Flux[j][i+1].d));

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

/*--- Step 10b -----------------------------------------------------------------
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

/*--- Step 11a -----------------------------------------------------------------
 * Update cell-centered variables in pG using x1-fluxes
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

/*--- Step 11b -----------------------------------------------------------------
 * Update cell-centered variables in pG using x2-fluxes
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

/*--- Step 13 ------------------------------------------------------------------
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

#if defined MHD || defined SHEARING_BOX
  if ((dhalf = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#else
  if(StaticGravPot != NULL){
    if ((dhalf = (Real**)calloc_2d_array(Nx2, Nx1, sizeof(Real))) == NULL)
      goto on_error;
  }
#endif

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
