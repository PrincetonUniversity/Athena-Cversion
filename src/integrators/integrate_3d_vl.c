#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_3d_vl.c
 *  \brief Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit MUSCL-Hancock (VL) integrator.  The variables updated are:
 *    - U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *    - B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, and the H-correction
 *   of Sanders et al.
 *   - For adb hydro, requires (9*Cons1DS + 3*Real + 1*ConsS) = 53 3D arrays
 *   - For adb mhd, requires   (9*Cons1DS + 9*Real + 1*ConsS) = 80 3D arrays
 *
 * REFERENCE: 
 * - J.M Stone & T.A. Gardiner, "A simple, unsplit Godunov method
 *   for multidimensional MHD", NewA 14, 139 (2009)
 *
 * - R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_3d_vl()
 * - integrate_destruct_3d()
 * - integrate_init_3d() */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifndef SPECIAL_RELATIVITY

#if defined(VL_INTEGRATOR)

/* The L/R states of primitive variables and fluxes at each cell face */
static Prim1DS ***Wl_x1Face=NULL, ***Wr_x1Face=NULL;
static Prim1DS ***Wl_x2Face=NULL, ***Wr_x2Face=NULL;
static Prim1DS ***Wl_x3Face=NULL, ***Wr_x3Face=NULL;
static Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
#ifdef FIRST_ORDER_FLUX_CORRECTION
static Cons1DS ***x1FluxP=NULL, ***x2FluxP=NULL, ***x3FluxP=NULL;
static Real ***emf1P=NULL, ***emf2P=NULL, ***emf3P=NULL;
static Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1, x3FD_k, x3FD_kp1;
#ifdef MHD
static Real emf1D_kj, emf1D_kp1j, emf1D_kjp1, emf1D_kp1jp1;
static Real emf2D_ki, emf2D_kip1, emf2D_kp1i, emf2D_kp1ip1;
static Real emf3D_ji, emf3D_jip1, emf3D_jp1i, emf3D_jp1ip1;
#endif
#endif

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
static Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W1d=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* conserved variables at t^{n+1/2} computed in predict step */
static ConsS ***Uhalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real ***eta1=NULL, ***eta2=NULL, ***eta3=NULL;
#endif

/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static ConsS **Flxiib=NULL, **Flxoib=NULL;
static ConsS **rFlxiib=NULL, **rFlxoib=NULL;
#endif
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - upwind CT method of GS (2005) for emf1
 *   integrate_emf2_corner() - upwind CT method of GS (2005) for emf2 
 *   integrate_emf3_corner() - upwind CT method of GS (2005) for emf3
 *   FixCell() - apply first-order correction to one cell
 *============================================================================*/
#ifdef MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
static void FixCell(GridS *pG, Int3Vect);
static void ApplyCorr(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_3d_vl(DomainS *pD)
 *  \brief 3D van Leer unsplit integrator for MHD. 
 */
void integrate_3d_vl(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  PrimS W,Whalf;
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  Real dt = pG->dt, hdt = 0.5*pG->dt, dx2=pG->dx2;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,Bx;
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
  Real M1n, dM2n; /* M1, dM2=(My+d*q*Omega_0*x) at time n */
  Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif
#ifdef FIRST_ORDER_FLUX_CORRECTION
  int flag_cell=0,negd=0,negP=0,superl=0,NaNFlux=0;
  Real Vsq;
  Int3Vect BadCell;
#endif
  int il=is-(nghost-1), iu=ie+(nghost-1);
  int jl=js-(nghost-1), ju=je+(nghost-1);
  int kl=ks-(nghost-1), ku=ke+(nghost-1);

#ifdef CYLINDRICAL
  Real Ekin,Emag,Ptot,B2sq;
  const Real *r=pG->r,*ri=pG->ri;
#ifdef FARGO
  Real Om,qshear;
#endif
#endif
  Real lsf=1.0,rsf=1.0,g;
  Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;

#if defined(CYLINDRICAL) && defined(FARGO)
  if (OrbitalProfile==NULL || ShearProfile==NULL)
    ath_error("[integrate_3d_vl]:  OrbitalProfile() and ShearProfile() *must* be defined.\n");
#endif

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        Uhalf[k][j][i] = pG->U[k][j][i];
#ifdef MHD
        B1_x1Face[k][j][i] = pG->B1i[k][j][i];
        B2_x2Face[k][j][i] = pG->B2i[k][j][i];
        B3_x3Face[k][j][i] = pG->B3i[k][j][i];
#endif /* MHD */
      }
    }
  }

/*=== STEP 1: Compute first-order fluxes at t^{n} in x1-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	U1d[i].d  = pG->U[k][j][i].d;
	U1d[i].Mx = pG->U[k][j][i].M1;
	U1d[i].My = pG->U[k][j][i].M2;
	U1d[i].Mz = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
	U1d[i].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
	U1d[i].By = pG->U[k][j][i].B2c;
	U1d[i].Bz = pG->U[k][j][i].B3c;
        Bxc[i] = pG->U[k][j][i].B1c;
        Bxi[i] = pG->B1i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[k][j][i].s[n];
#endif
      }

/*--- Step 1b ------------------------------------------------------------------
 * Compute first-order L/R states */

    for (i=is-nghost; i<=ie+nghost; i++) {
      W1d[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
    }

    for (i=il; i<=ie+nghost; i++) {
      Wl[i] = W1d[i-1];
      Wr[i] = W1d[i  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
      Ul[i] = Prim1D_to_Cons1D(&Wl[i], &Bxc[i-1]);
      Ur[i] = Prim1D_to_Cons1D(&Wr[i], &Bxc[i  ]);
    }

/*--- Step 1c ------------------------------------------------------------------
 * No source terms needed.
 */

/*--- Step 1d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (i=il; i<=ie+nghost; i++) {
        fluxes(Ul[i],Ur[i],Wl[i],Wr[i],Bxi[i],&x1Flux[k][j][i]);
      }
    }
  }

/*=== STEP 2: Compute first-order fluxes at t^{n} in x2-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
	U1d[j].d  = pG->U[k][j][i].d;
	U1d[j].Mx = pG->U[k][j][i].M2;
	U1d[j].My = pG->U[k][j][i].M3;
	U1d[j].Mz = pG->U[k][j][i].M1;
#ifndef BAROTROPIC
	U1d[j].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
	U1d[j].By = pG->U[k][j][i].B3c;
	U1d[j].Bz = pG->U[k][j][i].B1c;
        Bxc[j] = pG->U[k][j][i].B2c;
        Bxi[j] = pG->B2i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[k][j][i].s[n];
#endif
      }

/*--- Step 2b ------------------------------------------------------------------
 * Compute first-order L/R states */

      for (j=js-nghost; j<=je+nghost; j++) {
        W1d[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
      }

      for (j=jl; j<=je+nghost; j++) {
        Wl[j] = W1d[j-1];
        Wr[j] = W1d[j  ];

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
        Ul[j] = Prim1D_to_Cons1D(&Wl[j], &Bxc[j-1]);
        Ur[j] = Prim1D_to_Cons1D(&Wr[j], &Bxc[j  ]);
      }

/*--- Step 2c ------------------------------------------------------------------
 * No source terms needed
 */

/*--- Step 2d ------------------------------------------------------------------
 * Compute flux in x2-direction */

      for (j=jl; j<=je+nghost; j++) {
        fluxes(Ul[j],Ur[j],Wl[j],Wr[j],Bxi[j],&x2Flux[k][j][i]);
      }
    }
  }

/*=== STEP 3: Compute first-order fluxes at t^{n} in x3-direction ============*/
/* No source terms are needed since there is no temporal evolution */

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M3, M1, M2, E, B1c, B2c, s[n])
 */

  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
	U1d[k].d  = pG->U[k][j][i].d;
	U1d[k].Mx = pG->U[k][j][i].M3;
	U1d[k].My = pG->U[k][j][i].M1;
	U1d[k].Mz = pG->U[k][j][i].M2;
#ifndef BAROTROPIC
	U1d[k].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
	U1d[k].By = pG->U[k][j][i].B1c;
	U1d[k].Bz = pG->U[k][j][i].B2c;
        Bxc[k] = pG->U[k][j][i].B3c;
        Bxi[k] = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[k].s[n] = pG->U[k][j][i].s[n];
#endif
      }

/*--- Step 3b ------------------------------------------------------------------
 * Compute first-order L/R states */      
        
      for (k=ks-nghost; k<=ke+nghost; k++) {
        W1d[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
      }

      for (k=kl; k<=ke+nghost; k++) { 
        Wl[k] = W1d[k-1];
        Wr[k] = W1d[k  ]; 

/* Compute U from W in case Pfloor used in Cons1D_to_Prim1D */
        Ul[k] = Prim1D_to_Cons1D(&Wl[k], &Bxc[k-1]);
        Ur[k] = Prim1D_to_Cons1D(&Wr[k], &Bxc[k  ]);
      }

/*--- Step 3c ------------------------------------------------------------------
 * No source terms needed. 
 */

/*--- Step 3d ------------------------------------------------------------------
 * Compute flux in x1-direction */

      for (k=kl; k<=ke+nghost; k++) {
        fluxes(Ul[k],Ur[k],Wl[k],Wr[k],Bxi[k],&x3Flux[k][j][i]);
      }
    }
  }

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef MHD
  for (k=ks-nghost; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        Whalf = Cons_to_Prim(&pG->U[k][j][i]);
        emf1_cc[k][j][i] = (Whalf.B2c*Whalf.V3 - Whalf.B3c*Whalf.V2);
        emf2_cc[k][j][i] = (Whalf.B3c*Whalf.V1 - Whalf.B1c*Whalf.V3);
        emf3_cc[k][j][i] = (Whalf.B1c*Whalf.V2 - Whalf.B2c*Whalf.V1);
      }
    }
  }
  integrate_emf1_corner(pG);
  integrate_emf2_corner(pG);
  integrate_emf3_corner(pG);

/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
        q2 = hdt/(ri[i]*pG->dx2);
#endif
        B1_x1Face[k][j][i] += q3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                              q2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        B2_x2Face[k][j][i] += q1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                              q3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        q2 = hdt/(r[i]*pG->dx2);
#endif
        B3_x3Face[k][j][i] += q2*(    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                              q1*(rsf*emf2[k  ][j  ][i+1] - lsf*emf2[k][j][i]);
      }
#ifdef CYLINDRICAL
      q2 = hdt/(ri[iu+1]*pG->dx2);
#endif
      B1_x1Face[k][j][iu+1] += q3*(emf2[k+1][j  ][iu+1]-emf2[k][j][iu+1]) -
                               q2*(emf3[k  ][j+1][iu+1]-emf3[k][j][iu+1]);
    }
    for (i=il; i<=iu; i++) {
      B2_x2Face[k][ju+1][i] += q1*(emf3[k  ][ju+1][i+1]-emf3[k][ju+1][i]) -
                               q3*(emf1[k+1][ju+1][i  ]-emf1[k][ju+1][i]);
    }
  }
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      q2 = hdt/(r[i]*pG->dx2);
#endif
      B3_x3Face[ku+1][j][i] += q2*(    emf1[ku+1][j+1][i  ]-    emf1[ku+1][j][i]) -
                               q1*(rsf*emf2[ku+1][j  ][i+1]-lsf*emf2[ku+1][j][i]);
    }
  }

/*--- Step 4c ------------------------------------------------------------------
 * Compute cell-centered magnetic fields at half-timestep from average of
 * face-centered fields.
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        Uhalf[k][j][i].B1c = 0.5*(lsf*B1_x1Face[k][j][i] + rsf*B1_x1Face[k][j][i+1]);
        Uhalf[k][j][i].B2c = 0.5*(    B2_x2Face[k][j][i] +     B2_x2Face[k][j+1][i]);
        Uhalf[k][j][i].B3c = 0.5*(    B3_x3Face[k][j][i] +     B3_x3Face[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/*=== STEP 5: Update cell-centered variables to half-timestep ================*/

/*--- Step 5a ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x1-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        Uhalf[k][j][i].d   -= q1*(rsf*x1Flux[k][j][i+1].d  - lsf*x1Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q1*(rsf*x1Flux[k][j][i+1].Mx - lsf*x1Flux[k][j][i].Mx);
        Uhalf[k][j][i].M2  -= q1*(SQR(rsf)*x1Flux[k][j][i+1].My - SQR(lsf)*x1Flux[k][j][i].My);
        Uhalf[k][j][i].M3  -= q1*(rsf*x1Flux[k][j][i+1].Mz - lsf*x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q1*(rsf*x1Flux[k][j][i+1].E  - lsf*x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -= 
            q1*(rsf*x1Flux[k][j][i+1].s[n] - lsf*x1Flux[k][j][i  ].s[n]);
#endif
      }
    }
  }

/*--- Step 5b ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x2-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
        q2 = hdt/(r[i]*pG->dx2);
#endif
        Uhalf[k][j][i].d   -= q2*(x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz);
        Uhalf[k][j][i].M2  -= q2*(x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx);
        Uhalf[k][j][i].M3  -= q2*(x2Flux[k][j+1][i].My - x2Flux[k][j][i].My);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q2*(x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -= 
            q2*(x2Flux[k][j+1][i].s[n] - x2Flux[k][j  ][i].s[n]);
#endif
      }
    }
  }

/*--- Step 5c ------------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x3-fluxes
 */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].d   -= q3*(x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d );
        Uhalf[k][j][i].M1  -= q3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);
        Uhalf[k][j][i].M2  -= q3*(x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz);
        Uhalf[k][j][i].M3  -= q3*(x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx);
#ifndef BAROTROPIC
        Uhalf[k][j][i].E   -= q3*(x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          Uhalf[k][j][i].s[n] -=
            q3*(x3Flux[k+1][j][i].s[n] - x3Flux[k  ][j][i].s[n]);
#endif
      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*--- Step 5d ------------------------------------------------------------------
 * With first-order flux correction, save predict fluxes and emf3
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        x1FluxP[k][j][i] = x1Flux[k][j][i];
        x2FluxP[k][j][i] = x2Flux[k][j][i];
        x3FluxP[k][j][i] = x3Flux[k][j][i];
#ifdef MHD
        emf1P[k][j][i] = emf1[k][j][i];
        emf2P[k][j][i] = emf2[k][j][i];
        emf3P[k][j][i] = emf3[k][j][i];
#endif
      }
    }
  }

#if defined (MHD) && defined (SHEARING_BOX)
  get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* initialize remapped Fluxes */

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je+1; j++){
        Flxiib[k][j].d  = x1FluxP[k][j][is].d;
        Flxiib[k][j].M1 = x1FluxP[k][j][is].Mx;
        Flxiib[k][j].M2 = x1FluxP[k][j][is].My;
        Flxiib[k][j].M3 = x1FluxP[k][j][is].Mz;

        Flxoib[k][j].d  = x1FluxP[k][j][ie+1].d;
        Flxoib[k][j].M1 = x1FluxP[k][j][ie+1].Mx;
        Flxoib[k][j].M2 = x1FluxP[k][j][ie+1].My;
        Flxoib[k][j].M3 = x1FluxP[k][j][ie+1].Mz;
#ifndef BAROTROPIC
        Flxiib[k][j].E = x1FluxP[k][j][is].E;
        Flxoib[k][j].E = x1FluxP[k][j][ie+1].E;
#endif
#ifdef MHD
        Flxiib[k][j].B1c = emf1P[k][j][is];
        Flxiib[k][j].B2c = emf2P[k][j][is];
        Flxiib[k][j].B3c = emf3P[k][j][is];

        Flxoib[k][j].B1c = emf1P[k][j][ie+1];
        Flxoib[k][j].B2c = emf2P[k][j][ie+1];
        Flxoib[k][j].B3c = emf3P[k][j][ie+1];
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          Flxiib[k][j].s[n] = x1FluxP[k][j][is].s[n];
          Flxoib[k][j].s[n] = x1FluxP[k][j][ie+1].s[n];
        }
#endif
      }
    }

/* compute remapped Fluxes from opposite side of grid */

    if (my_iproc == 0) {
      RemapFlx_ix1(pD, Flxiib, Flxoib, rFlxiib);
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapFlx_ox1(pD, Flxiib, Flxoib, rFlxoib);
    }

/* Now average fluxes and remapped fluxes */

    if (my_iproc == 0) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
          emf2P[k][j][is] = 0.5*(emf2P[k][j][is] + rFlxiib[k][j].B2c);
          x1FluxP[k][j][is].d  = 0.5*(x1FluxP[k][j][is].d  + rFlxiib[k][j].d);
          x1FluxP[k][j][is].Mx = 0.5*(x1FluxP[k][j][is].Mx + rFlxiib[k][j].M1);
          x1FluxP[k][j][is].My = 0.5*(x1FluxP[k][j][is].My + rFlxiib[k][j].M2);
          x1FluxP[k][j][is].Mz = 0.5*(x1FluxP[k][j][is].Mz + rFlxiib[k][j].M3);
        }
      }
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
          emf2P[k][j][ie+1] = 0.5*(emf2P[k][j][ie+1] + rFlxoib[k][j].B2c);
          x1FluxP[k][j][ie+1].d =0.5*(x1FluxP[k][j][ie+1].d  + rFlxoib[k][j].d);
          x1FluxP[k][j][ie+1].Mx=0.5*(x1FluxP[k][j][ie+1].Mx + rFlxoib[k][j].M1);
          x1FluxP[k][j][ie+1].My=0.5*(x1FluxP[k][j][ie+1].My + rFlxoib[k][j].M2);
          x1FluxP[k][j][ie+1].Mz=0.5*(x1FluxP[k][j][ie+1].Mz + rFlxoib[k][j].M3);
        }
      }
    }
#endif /* MHD & SHEARING_BOX */


#endif

/*=== STEP 6: Add source terms to predict values at half-timestep ============*/

/*--- Step 6a ------------------------------------------------------------------
 * Add source terms from a static gravitational potential for 0.5*dt to predict
 * step.  To improve conservation of total energy, we average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phic = (*StaticGravPot)(x1,x2,x3);
          phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          g = (phir-phil)*dx1i;
#ifdef CYLINDRICAL
          rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
          q2 = hdt/(r[i]*pG->dx2);
#ifdef FARGO
          g -= r[i]*SQR((*OrbitalProfile)(r[i]));
#endif
#endif /* CYLINDRICAL */
          Uhalf[k][j][i].M1 -= hdt*pG->U[k][j][i].d*g;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q1*(lsf*x1Flux[k][j][i  ].d*(phic - phil)
                                + rsf*x1Flux[k][j][i+1].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          Uhalf[k][j][i].M2 -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q2*(x2Flux[k][j  ][i].d*(phic - phil)
                                + x2Flux[k][j+1][i].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          Uhalf[k][j][i].M3 -= q3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
          Uhalf[k][j][i].E -= q3*(x3Flux[k  ][j][i].d*(phic - phil)
                                + x3Flux[k+1][j][i].d*(phir - phic));
#endif
        }
      }
    }
  }

/*--- Step 6b ------------------------------------------------------------------
 * Add source terms for self gravity for 0.5*dt to predict step.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        phic = pG->Phi[k][j][i];
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i+1]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i-1]);

        Uhalf[k][j][i].M1 -= q1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q1*(x1Flux[k][j][i  ].d*(phic - phil)
                              + x1Flux[k][j][i+1].d*(phir - phic));
#endif
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j+1][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j-1][i]);

        Uhalf[k][j][i].M2 -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q2*(x2Flux[k][j  ][i].d*(phic - phil)
                              + x2Flux[k][j+1][i].d*(phir - phic));
#endif
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k+1][j][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k-1][j][i]);

        Uhalf[k][j][i].M3 -= q3*(phir-phil)*pG->U[k][j][i].d; 
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q3*(x3Flux[k  ][j][i].d*(phic - phil)
                              + x3Flux[k+1][j][i].d*(phir - phic));
#endif
      }
    }
  }
#endif /* SELF_GRAVITY */

/*--- Step 6c -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
 * x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
 */

#ifdef SHEARING_BOX
  if (ShearingBoxPot != NULL){
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        phic = (*ShearingBoxPot)((x1            ),x2,x3);
        phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

        Uhalf[k][j][i].M1 -= q1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q1*(x1Flux[k][j][i  ].d*(phic - phil)
                              + x1Flux[k][j][i+1].d*(phir - phic));
#endif
        phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);

        Uhalf[k][j][i].M2 -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q2*(x2Flux[k][j  ][i].d*(phic - phil)
                              + x2Flux[k][j+1][i].d*(phir - phic));
#endif
        phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
        phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));

        Uhalf[k][j][i].M3 -= q3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Uhalf[k][j][i].E -= q3*(x3Flux[k  ][j][i].d*(phic - phil)
                              + x3Flux[k+1][j][i].d*(phir - phic));
#endif
      }
    }
  }}

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Uhalf[k][j][i].M1 += pG->dt*Omega_0*pG->U[k][j][i].M2;
#ifdef FARGO
        Uhalf[k][j][i].M2 += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
#else
        Uhalf[k][j][i].M2 -= pG->dt*Omega_0*pG->U[k][j][i].M1;
#endif
      }
    }
  }
#endif /* SHEARING_BOX */

#if defined(CYLINDRICAL) && defined(FARGO)
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Om = (*OrbitalProfile)(r[i]);
        qshear = (*ShearProfile)(r[i]);
        /* This *is* a half-timestep update below (see note above) */
        Uhalf[k][j][i].M1 += pG->dt*Om*pG->U[k][j][i].M2;
        Uhalf[k][j][i].M2 += hdt*(qshear - 2.0)*Om*pG->U[k][j][i].M1;
      }
    }
  }
#endif

/*--- Step 6d ------------------------------------------------------------------
 * Add the geometric source-term now using cell-centered conserved variables
 *   at time t^n 
 */

#ifdef CYLINDRICAL
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        Ekin = 0.5*(SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#ifdef MHD
        B2sq = SQR(pG->U[k][j][i].B2c);
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + B2sq + SQR(pG->U[k][j][i].B3c));
#else
        B2sq = 0.0;
        Emag = 0.0;
#endif

#ifdef ISOTHERMAL
        Ptot = Iso_csound2*pG->U[k][j][i].d;
#else
        Ptot = Gamma_1*(pG->U[k][j][i].E - Ekin - Emag);
#endif
        Ptot = MAX(Ptot,TINY_NUMBER);
        Ptot += Emag;

        Uhalf[k][j][i].M1 += hdt*(SQR(pG->U[k][j][i].M2)/pG->U[k][j][i].d - B2sq + Ptot)/r[i];
      }
    }
  }
#endif /* CYLINDRICAL */

/*=== STEP 7: Compute second-order L/R x1-interface states ===================*/

/*--- Step 7a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=il; i<=iu; i++) {
        U1d[i].d  = Uhalf[k][j][i].d;
        U1d[i].Mx = Uhalf[k][j][i].M1;
        U1d[i].My = Uhalf[k][j][i].M2;
        U1d[i].Mz = Uhalf[k][j][i].M3;
#ifndef BAROTROPIC
        U1d[i].E  = Uhalf[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[i].By = Uhalf[k][j][i].B2c;
        U1d[i].Bz = Uhalf[k][j][i].B3c;
        Bxc[i] = Uhalf[k][j][i].B1c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[i].s[n] = Uhalf[k][j][i].s[n];
#endif
      }

/*--- Step 7b ------------------------------------------------------------------
 * Compute L and R states at x1-interfaces, store in 3D array
 */

      for (i=il; i<=iu; i++) {
        W1d[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
      }

      lr_states(pG,W1d,Bxc,pG->dt,pG->dx1,is,ie,Wl,Wr,1);

      for (i=is; i<=ie+1; i++) {
        Wl_x1Face[k][j][i] = Wl[i];
        Wr_x1Face[k][j][i] = Wr[i];
      }
    }
  }

/*=== STEP 8: Compute second-order L/R x2-interface states ===================*/

/*--- Step 8a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (i=is-1; i<=ie+1; i++) {
      for (j=jl; j<=ju; j++) {
        U1d[j].d  = Uhalf[k][j][i].d;
        U1d[j].Mx = Uhalf[k][j][i].M2;
        U1d[j].My = Uhalf[k][j][i].M3;
        U1d[j].Mz = Uhalf[k][j][i].M1;
#ifndef BAROTROPIC
        U1d[j].E  = Uhalf[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[j].By = Uhalf[k][j][i].B3c;
        U1d[j].Bz = Uhalf[k][j][i].B1c;
        Bxc[j] = Uhalf[k][j][i].B2c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[j].s[n] = Uhalf[k][j][i].s[n];
#endif
      }

/*--- Step 8b ------------------------------------------------------------------
 * Compute L and R states at x2-interfaces, store in 3D array
 */

      for (j=jl; j<=ju; j++) {
        W1d[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
      }

#ifdef CYLINDRICAL
      dx2 = r[i]*pG->dx2;
#endif
      lr_states(pG,W1d,Bxc,pG->dt,dx2,js,je,Wl,Wr,2);

      for (j=js; j<=je+1; j++) {
        Wl_x2Face[k][j][i] = Wl[j];
        Wr_x2Face[k][j][i] = Wr[j];
      }
    }
  }

/*=== STEP 9: Compute second-order L/R x3-interface states ===================*/

/*--- Step 9a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M3, M1, M2, E, B1c, B2c, s[n])
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      for (k=kl; k<=ku; k++) {
        U1d[k].d  = Uhalf[k][j][i].d;
        U1d[k].Mx = Uhalf[k][j][i].M3;
        U1d[k].My = Uhalf[k][j][i].M1;
        U1d[k].Mz = Uhalf[k][j][i].M2;
#ifndef BAROTROPIC
        U1d[k].E  = Uhalf[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[k].By = Uhalf[k][j][i].B1c;
        U1d[k].Bz = Uhalf[k][j][i].B2c;
        Bxc[k] = Uhalf[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[k].s[n] = Uhalf[k][j][i].s[n];
#endif
      }

/*--- Step 9b ------------------------------------------------------------------
 * Compute L and R states at x3-interfaces, store in 3D array
 */

      for (k=kl; k<=ku; k++) {
        W1d[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
      }

      lr_states(pG,W1d,Bxc,pG->dt,pG->dx3,ks,ke,Wl,Wr,3);

      for (k=ks; k<=ke+1; k++) {
        Wl_x3Face[k][j][i] = Wl[k];
        Wr_x3Face[k][j][i] = Wr[k];
      }
    }
  }

/*=== STEP 10: Compute 3D x1-Flux, x2-Flux, x3-Flux ==========================*/

/*--- Step 10a -----------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction, if needed.
 */

#ifdef H_CORRECTION
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=iu; i++) {
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x1Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x1Face[k][j][i]),&Bx);
        lambdar = Ur_x1Face[k][j][i].Mx/Ur_x1Face[k][j][i].d + cfr;
        lambdal = Ul_x1Face[k][j][i].Mx/Ul_x1Face[k][j][i].d - cfl;
        eta1[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=ju; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x2Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x2Face[k][j][i]),&Bx);
        lambdar = Ur_x2Face[k][j][i].Mx/Ur_x2Face[k][j][i].d + cfr;
        lambdal = Ul_x2Face[k][j][i].Mx/Ul_x2Face[k][j][i].d - cfl;
        eta2[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ku; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x3Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x3Face[k][j][i]),&Bx);
        lambdar = Ur_x3Face[k][j][i].Mx/Ur_x3Face[k][j][i].d + cfr;
        lambdal = Ul_x3Face[k][j][i].Mx/Ul_x3Face[k][j][i].d - cfl;
        eta3[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }
#endif /* H_CORRECTION */

/*--- Step 10b -----------------------------------------------------------------
 * Compute second-order fluxes in x1-direction
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta2[k][j][i-1],eta2[k][j][i]);
        etah = MAX(etah,eta2[k][j+1][i-1]);
        etah = MAX(etah,eta2[k][j+1][i  ]);

        etah = MAX(etah,eta3[k  ][j][i-1]);
        etah = MAX(etah,eta3[k  ][j][i  ]);
        etah = MAX(etah,eta3[k+1][j][i-1]);
        etah = MAX(etah,eta3[k+1][j][i  ]);

        etah = MAX(etah,eta1[k  ][j][i  ]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x1Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x1Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x1Face[k][j][i],Wr_x1Face[k][j][i],Bx,
               &x1Flux[k][j][i]);

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux Nan'ed */
        if ((x1Flux[k][j][i].d  != x1Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x1Flux[k][j][i].E  != x1Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x1Flux[k][j][i].By != x1Flux[k][j][i].By) ||
            (x1Flux[k][j][i].Bz != x1Flux[k][j][i].Bz) ||
#endif
            (x1Flux[k][j][i].Mx != x1Flux[k][j][i].Mx) ||
            (x1Flux[k][j][i].My != x1Flux[k][j][i].My) ||
            (x1Flux[k][j][i].Mz != x1Flux[k][j][i].Mz)) {
          x1Flux[k][j][i] = x1FluxP[k][j][i];
          NaNFlux++;
        }
#endif

      }
    }
  }

/*--- Step 10c -----------------------------------------------------------------
 * Compute second-order fluxes in x2-direction
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k][j-1][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k][j-1][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta3[k  ][j-1][i]);
        etah = MAX(etah,eta3[k  ][j  ][i]);
        etah = MAX(etah,eta3[k+1][j-1][i]);
        etah = MAX(etah,eta3[k+1][j  ][i]);

        etah = MAX(etah,eta2[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x2Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x2Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x2Face[k][j][i],Wr_x2Face[k][j][i],Bx,
               &x2Flux[k][j][i]);

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x2Flux[k][j][i].d  != x2Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x2Flux[k][j][i].E  != x2Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x2Flux[k][j][i].By != x2Flux[k][j][i].By) ||
            (x2Flux[k][j][i].Bz != x2Flux[k][j][i].Bz) ||
#endif
            (x2Flux[k][j][i].Mx != x2Flux[k][j][i].Mx) ||
            (x2Flux[k][j][i].My != x2Flux[k][j][i].My) ||
            (x2Flux[k][j][i].Mz != x2Flux[k][j][i].Mz)) {
          x2Flux[k][j][i] = x2FluxP[k][j][i];
          NaNFlux++;
        }
#endif

      }
    }
  }

/*--- Step 10d -----------------------------------------------------------------
 * Compute second-order fluxes in x3-direction
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k-1][j][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k-1][j][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta2[k-1][j  ][i]);
        etah = MAX(etah,eta2[k  ][j  ][i]);
        etah = MAX(etah,eta2[k-1][j+1][i]);
        etah = MAX(etah,eta2[k  ][j+1][i]);

        etah = MAX(etah,eta3[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        Ul[i] = Prim1D_to_Cons1D(&Wl_x3Face[k][j][i],&Bx);
        Ur[i] = Prim1D_to_Cons1D(&Wr_x3Face[k][j][i],&Bx);

        fluxes(Ul[i],Ur[i],Wl_x3Face[k][j][i],Wr_x3Face[k][j][i],Bx,
               &x3Flux[k][j][i]);

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x3Flux[k][j][i].d  != x3Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x3Flux[k][j][i].E  != x3Flux[k][j][i].E)  ||
#endif
#ifdef MHD
            (x3Flux[k][j][i].By != x3Flux[k][j][i].By) ||
            (x3Flux[k][j][i].Bz != x3Flux[k][j][i].Bz) ||
#endif
            (x3Flux[k][j][i].Mx != x3Flux[k][j][i].Mx) ||
            (x3Flux[k][j][i].My != x3Flux[k][j][i].My) ||
            (x3Flux[k][j][i].Mz != x3Flux[k][j][i].Mz)) {
          x3Flux[k][j][i] = x3FluxP[k][j][i];
          NaNFlux++;
        }
#endif

      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step10] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif

/*=== STEP 11: Update face-centered B for a full timestep ====================*/
        
/*--- Step 11a -----------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at the half-time-step.
 */

#ifdef MHD
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        Whalf = Cons_to_Prim(&Uhalf[k][j][i]);
        emf1_cc[k][j][i] = (Whalf.B2c*Whalf.V3 - Whalf.B3c*Whalf.V2);
        emf2_cc[k][j][i] = (Whalf.B3c*Whalf.V1 - Whalf.B1c*Whalf.V3);
        emf3_cc[k][j][i] = (Whalf.B1c*Whalf.V2 - Whalf.B2c*Whalf.V1);
      }
    }
  }
#endif

/*--- Step 11b -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef MHD
  integrate_emf1_corner(pG);
  integrate_emf2_corner(pG);
  integrate_emf3_corner(pG);
#endif

/* Remap Fluxes at is and ie+1 to conserve quantities in shearing box */

#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* initialize remapped Fluxes */

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je+1; j++){
        Flxiib[k][j].d = x1Flux[k][j][is].d;
        Flxiib[k][j].M1 = x1Flux[k][j][is].Mx;
        Flxiib[k][j].M2 = x1Flux[k][j][is].My;
        Flxiib[k][j].M3 = x1Flux[k][j][is].Mz;

        Flxoib[k][j].d = x1Flux[k][j][ie+1].d;
        Flxoib[k][j].M1 = x1Flux[k][j][ie+1].Mx;
        Flxoib[k][j].M2 = x1Flux[k][j][ie+1].My;
        Flxoib[k][j].M3 = x1Flux[k][j][ie+1].Mz;
#ifndef BAROTROPIC
        Flxiib[k][j].E = x1Flux[k][j][is].E;
        Flxoib[k][j].E = x1Flux[k][j][ie+1].E;
#endif
#ifdef MHD
        Flxiib[k][j].B1c = emf1[k][j][is];
        Flxiib[k][j].B2c = emf2[k][j][is];
        Flxiib[k][j].B3c = emf3[k][j][is];

        Flxoib[k][j].B1c = emf1[k][j][ie+1];
        Flxoib[k][j].B2c = emf2[k][j][ie+1];
        Flxoib[k][j].B3c = emf3[k][j][ie+1];
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          Flxiib[k][j].s[n] = x1Flux[k][j][is].s[n];
          Flxoib[k][j].s[n] = x1Flux[k][j][ie+1].s[n];
        }
#endif
      }
    }

/* compute remapped Fluxes from opposite side of grid */

    if (my_iproc == 0) {
      RemapFlx_ix1(pD, Flxiib, Flxoib, rFlxiib);
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapFlx_ox1(pD, Flxiib, Flxoib, rFlxoib);
    }
  
/* Now average fluxes and remapped fluxes */

    if (my_iproc == 0) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
#ifdef MHD
          emf2[k][j][is] = 0.5*(emf2[k][j][is] + rFlxiib[k][j].B2c);
#endif
          x1Flux[k][j][is].d  = 0.5*(x1Flux[k][j][is].d  + rFlxiib[k][j].d);
          x1Flux[k][j][is].Mx = 0.5*(x1Flux[k][j][is].Mx + rFlxiib[k][j].M1);
          x1Flux[k][j][is].My = 0.5*(x1Flux[k][j][is].My + rFlxiib[k][j].M2);
          x1Flux[k][j][is].Mz = 0.5*(x1Flux[k][j][is].Mz + rFlxiib[k][j].M3);
        }
      }
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
#ifdef MHD
          emf2[k][j][ie+1] = 0.5*(emf2[k][j][ie+1] + rFlxoib[k][j].B2c);
#endif
          x1Flux[k][j][ie+1].d =0.5*(x1Flux[k][j][ie+1].d  + rFlxoib[k][j].d);
          x1Flux[k][j][ie+1].Mx=0.5*(x1Flux[k][j][ie+1].Mx + rFlxoib[k][j].M1);
          x1Flux[k][j][ie+1].My=0.5*(x1Flux[k][j][ie+1].My + rFlxoib[k][j].M2);
          x1Flux[k][j][ie+1].Mz=0.5*(x1Flux[k][j][ie+1].Mz + rFlxoib[k][j].M3);
        }
      }
    }
#endif /* SHEARING_BOX */


/*--- Step 11c -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        dtodx2 = pG->dt/(ri[i]*pG->dx2);
#endif
        pG->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                            dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        pG->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                            dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
        pG->B3i[k][j][i] += dtodx2*(    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                            dtodx1*(rsf*emf2[k  ][j  ][i+1] - lsf*emf2[k][j][i]);
      }
#ifdef CYLINDRICAL
      dtodx2 = pG->dt/(ri[ie+1]*pG->dx2);
#endif
      pG->B1i[k][j][ie+1] +=
        dtodx3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
        dtodx2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[k][je+1][i] +=
        dtodx1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
        dtodx3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
      pG->B3i[ke+1][j][i] += 
        dtodx2*(    emf1[ke+1][j+1][i  ] -     emf1[ke+1][j][i]) -
        dtodx1*(rsf*emf2[ke+1][j  ][i+1] - lsf*emf2[ke+1][j][i]);
    }
  }

/*--- Step 11d -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
      }
    }
  }
#endif

/*=== STEP 12: Add source terms for a full timestep using n+1/2 states =======*/
       
/*--- Step 12a -----------------------------------------------------------------
 * Add gravitational (or shearing box) source terms as a Static Potential.
 *   A Crank-Nicholson update is used for shearing box terms.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;
  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++) {
      for(i=is; i<=ie; i++) {
	cc_pos(pG,i,j,k,&x1,&x2,&x3);

/* Store the current state */
	M1n  = pG->U[k][j][i].M1;
#ifdef FARGO
	dM2n = pG->U[k][j][i].M2;
#else
	dM2n = pG->U[k][j][i].M2 + qom*x1*pG->U[k][j][i].d;
#endif

/* Calculate the flux for the y-momentum fluctuation */
	frx1_dM2 = x1Flux[k][j][i+1].My;
	flx1_dM2 = x1Flux[k][j][i  ].My;
	frx2_dM2 = x2Flux[k][j+1][i].Mx;
	flx2_dM2 = x2Flux[k][j  ][i].Mx;
	frx3_dM2 = x3Flux[k+1][j][i].Mz;
	flx3_dM2 = x3Flux[k  ][j][i].Mz;
#ifndef FARGO
	frx1_dM2 += qom*(x1+0.5*pG->dx1)*x1Flux[k][j][i+1].d;
	flx1_dM2 += qom*(x1-0.5*pG->dx1)*x1Flux[k][j][i  ].d;
	frx2_dM2 += qom*(x1            )*x2Flux[k][j+1][i].d;
	flx2_dM2 += qom*(x1            )*x2Flux[k][j  ][i].d;
	frx3_dM2 += qom*(x1            )*x3Flux[k+1][j][i].d;
	flx3_dM2 += qom*(x1            )*x3Flux[k  ][j][i].d;
#endif

/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
	M1e = M1n - q1*(x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx)
	          - q2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz)
	          - q3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);

	dM2e = dM2n - q1*(frx1_dM2 - flx1_dM2)
	            - q2*(frx2_dM2 - flx2_dM2) 
                    - q3*(frx3_dM2 - flx3_dM2);

/* Update the 1- and 2-momentum for the Coriolis and tidal
 * potential momentum source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

	pG->U[k][j][i].M1 += (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
	pG->U[k][j][i].M2 += 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;

#ifndef FARGO
	pG->U[k][j][i].M2 -= 0.5*qshear*om_dt*
           (x1Flux[k][j][i].d + x1Flux[k][j][i+1].d);
#endif

/* Update the energy for a fixed potential.
 * This update is identical to non-SHEARING_BOX below  */

	phic = (*ShearingBoxPot)((x1            ),x2,x3);
	phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
	phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
	pG->U[k][j][i].E -= dtodx1*(x1Flux[k][j][i  ].d*(phic - phil) +
                                    x1Flux[k][j][i+1].d*(phir - phic));
#endif

	phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
	phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
#ifndef BAROTROPIC
	pG->U[k][j][i].E -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil) +
                                    x2Flux[k][j+1][i].d*(phir - phic));
#endif

	phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
	phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
#ifndef BAROTROPIC
	pG->U[k][j][i].E -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil) +
                                    x3Flux[k+1][j][i].d*(phir - phic));
#endif
      }
    }
  }
#endif /* SHEARING_BOX */


  if (StaticGravPot != NULL){
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phic = (*StaticGravPot)(x1,x2,x3);
          phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          g = (phir-phil)*dx1i;
#ifdef CYLINDRICAL
          rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
          dtodx2 = pG->dt/(r[i]*pG->dx2);
#ifdef FARGO
          g -= r[i]*SQR((*OrbitalProfile)(r[i]));
#endif
#endif
          pG->U[k][j][i].M1 -= pG->dt*Uhalf[k][j][i].d*g;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx1*(lsf*x1Flux[k][j][i  ].d*(phic - phil)
                                    + rsf*x1Flux[k][j][i+1].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          pG->U[k][j][i].M2 -= dtodx2*(phir-phil)*Uhalf[k][j][i].d;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil)
                                    + x2Flux[k][j+1][i].d*(phir - phic));
#endif
          phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          pG->U[k][j][i].M3 -= dtodx3*(phir-phil)*Uhalf[k][j][i].d;
#ifndef BAROTROPIC
          pG->U[k][j][i].E -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil)
                                    + x3Flux[k+1][j][i].d*(phir - phic));
#endif
        }
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

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j][i  ]);
        phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+1]);

/* gx, gy and gz centered at L and R x1-faces */
        gxl = (pG->Phi[k][j][i-1] - pG->Phi[k][j][i  ])/(pG->dx1);
        gxr = (pG->Phi[k][j][i  ] - pG->Phi[k][j][i+1])/(pG->dx1);

        gyl = 0.25*((pG->Phi[k][j-1][i-1] - pG->Phi[k][j+1][i-1]) +
                    (pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]))/(pG->dx2);
        gyr = 0.25*((pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]) +
                    (pG->Phi[k][j-1][i+1] - pG->Phi[k][j+1][i+1]))/(pG->dx2);

        gzl = 0.25*((pG->Phi[k-1][j][i-1] - pG->Phi[k+1][j][i-1]) +
                    (pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]))/(pG->dx3);
        gzr = 0.25*((pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) +
                    (pG->Phi[k-1][j][i+1] - pG->Phi[k+1][j][i+1]))/(pG->dx3);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
        flx_m1l = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        flx_m1r = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        flx_m2l = gxl*gyl/four_pi_G;
        flx_m2r = gxr*gyr/four_pi_G;

        flx_m3l = gxl*gzl/four_pi_G;
        flx_m3r = gxr*gzr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
        pG->U[k][j][i].M1 -= dtodx1*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= dtodx1*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= dtodx1*(flx_m3r - flx_m3l);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -= dtodx1*(x1Flux[k][j][i  ].d*(phic - phil) +
                                    x1Flux[k][j][i+1].d*(phir - phic));
#endif /* BAROTROPIC */
      }
    }
  }

/* Add fluxes and source terms due to (d/dx2) terms  */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k][j  ][i]);
        phir = 0.5*(pG->Phi[k][j  ][i] + pG->Phi[k][j+1][i]);

/* gx, gy and gz centered at L and R x2-faces */
        gxl = 0.25*((pG->Phi[k][j-1][i-1] - pG->Phi[k][j-1][i+1]) +
                    (pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]))/(pG->dx1);
        gxr = 0.25*((pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]) +
                    (pG->Phi[k][j+1][i-1] - pG->Phi[k][j+1][i+1]))/(pG->dx1);

        gyl = (pG->Phi[k][j-1][i] - pG->Phi[k][j  ][i])/(pG->dx2);
        gyr = (pG->Phi[k][j  ][i] - pG->Phi[k][j+1][i])/(pG->dx2);

        gzl = 0.25*((pG->Phi[k-1][j-1][i] - pG->Phi[k+1][j-1][i]) +
                    (pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]) )/(pG->dx3);
        gzr = 0.25*((pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]) +
                    (pG->Phi[k-1][j+1][i] - pG->Phi[k+1][j+1][i]) )/(pG->dx3);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
        flx_m1l = gyl*gxl/four_pi_G;
        flx_m1r = gyr*gxr/four_pi_G;

        flx_m2l = 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        flx_m2r = 0.5*(gyr*gyr-gxr*gxr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        flx_m3l = gyl*gzl/four_pi_G;
        flx_m3r = gyr*gzr/four_pi_G;

/* Update momenta and energy with d/dx2 terms  */
        pG->U[k][j][i].M1 -= dtodx2*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= dtodx2*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= dtodx2*(flx_m3r - flx_m3l);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil) +
                                    x2Flux[k][j+1][i].d*(phir - phic));
#endif /* BAROTROPIC */
      }
    }
  }

/* Add fluxes and source terms due to (d/dx3) terms  */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        phic = pG->Phi[k][j][i];
        phil = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k  ][j][i]);
        phir = 0.5*(pG->Phi[k  ][j][i] + pG->Phi[k+1][j][i]);

/* gx, gy and gz centered at L and R x3-faces */
        gxl = 0.25*((pG->Phi[k-1][j][i-1] - pG->Phi[k-1][j][i+1]) +
                    (pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]))/(pG->dx1);
        gxr = 0.25*((pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]) +
                    (pG->Phi[k+1][j][i-1] - pG->Phi[k+1][j][i+1]))/(pG->dx1);

        gyl = 0.25*((pG->Phi[k-1][j-1][i] - pG->Phi[k-1][j+1][i]) +
                    (pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]))/(pG->dx2);
        gyr = 0.25*((pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]) +
                    (pG->Phi[k+1][j-1][i] - pG->Phi[k+1][j+1][i]))/(pG->dx2);

        gzl = (pG->Phi[k-1][j][i] - pG->Phi[k  ][j][i])/(pG->dx3);
        gzr = (pG->Phi[k  ][j][i] - pG->Phi[k+1][j][i])/(pG->dx3);

/* momentum fluxes in x3.  2nd term is needed only if Jean's swindle used */
        flx_m1l = gzl*gxl/four_pi_G;
        flx_m1r = gzr*gxr/four_pi_G;

        flx_m2l = gzl*gyl/four_pi_G;
        flx_m2r = gzr*gyr/four_pi_G;

        flx_m3l = 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
        flx_m3r = 0.5*(gzr*gzr-gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx3 terms  */
        pG->U[k][j][i].M1 -= dtodx3*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= dtodx3*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= dtodx3*(flx_m3r - flx_m3l);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil) +
                                    x3Flux[k+1][j][i].d*(phir - phic));
#endif /* BAROTROPIC */
      }
    }
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        pG->x1MassFlux[k][j][i] = x1Flux[k][j][i].d;
        pG->x2MassFlux[k][j][i] = x2Flux[k][j][i].d;
        pG->x3MassFlux[k][j][i] = x3Flux[k][j][i].d;
      }
    }
  }
#endif /* SELF_GRAVITY */

/*--- Step 12c -----------------------------------------------------------------
 * Add the geometric source-term now using cell-centered conserved variables
 *   at time t^{n+1/2}
 */

#ifdef CYLINDRICAL
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {

        Ekin = 0.5*(SQR(Uhalf[k][j][i].M1)+SQR(Uhalf[k][j][i].M2)+SQR(Uhalf[k][j][i].M3))/Uhalf[k][j][i].d;
#ifdef MHD
        B2sq = SQR(Uhalf[k][j][i].B2c);
        Emag = 0.5*(SQR(Uhalf[k][j][i].B1c) + B2sq + SQR(Uhalf[k][j][i].B3c));
#else
        B2sq = 0.0;
        Emag = 0.0;
#endif

#ifdef ISOTHERMAL
        Ptot = Iso_csound2*Uhalf[k][j][i].d;
#else
        Ptot = Gamma_1*(Uhalf[k][j][i].E - Ekin - Emag);
#endif
        Ptot = MAX(Ptot,TINY_NUMBER);
        Ptot += Emag;

        pG->U[k][j][i].M1 += pG->dt*(SQR(Uhalf[k][j][i].M2)/Uhalf[k][j][i].d - B2sq + Ptot)/r[i];

#ifdef FARGO
        /* Use average values to apply source terms for full time-step */
        Om = (*OrbitalProfile)(r[i]);
        qshear = (*ShearProfile)(r[i]);
        pG->U[k][j][i].M1 += pG->dt*(2.0*Om*Uhalf[k][j][i].M2);
        pG->U[k][j][i].M2 += pG->dt*(Om*(qshear-2.0)*Uhalf[k][j][i].M1);
#endif /* FARGO */
      }
    }
  }
#endif /* CYLINDRICAL */

/*=== STEP 13: Update cell-centered values for a full timestep ===============*/

/*--- Step 13a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x1-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        pG->U[k][j][i].d  -= dtodx1*(rsf*x1Flux[k][j][i+1].d  - lsf*x1Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx1*(rsf*x1Flux[k][j][i+1].Mx - lsf*x1Flux[k][j][i].Mx);
        pG->U[k][j][i].M2 -= dtodx1*(SQR(rsf)*x1Flux[k][j][i+1].My - SQR(lsf)*x1Flux[k][j][i].My);
        pG->U[k][j][i].M3 -= dtodx1*(rsf*x1Flux[k][j][i+1].Mz - lsf*x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  -= dtodx1*(rsf*x1Flux[k][j][i+1].E  - lsf*x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx1*(rsf*x1Flux[k][j][i+1].s[n]
                                       - lsf*x1Flux[k][j][i  ].s[n]);
#endif
      }
    }
  }

/*--- Step 13b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x2-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
        pG->U[k][j][i].d  -= dtodx2*(x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz);
        pG->U[k][j][i].M2 -= dtodx2*(x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx);
        pG->U[k][j][i].M3 -= dtodx2*(x2Flux[k][j+1][i].My - x2Flux[k][j][i].My);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  -= dtodx2*(x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx2*(x2Flux[k][j+1][i].s[n]
                                       - x2Flux[k][j  ][i].s[n]);
#endif
      }
    }
  }

/*--- Step 13c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d  -= dtodx3*(x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);
        pG->U[k][j][i].M2 -= dtodx3*(x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz);
        pG->U[k][j][i].M3 -= dtodx3*(x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  -= dtodx3*(x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx3*(x3Flux[k+1][j][i].s[n]
                                       - x3Flux[k  ][j][i].s[n]);
#endif
      }
    }
  }

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*=== STEP 14: First-order flux correction ===================================*/

/* If cell-centered d or P have gone negative, or if v^2 > 1 in SR, correct
 * by using 1st order predictor fluxes */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        W = Cons_to_Prim(&(pG->U[k][j][i]));
        if (W.d < 0.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negd++;
        }
#ifndef ISOTHERMAL
        if (W.P < 0.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negP++;
        }
#endif
        if (flag_cell != 0) {
          FixCell(pG, BadCell);
          flag_cell=0;
        }
      }
    }
  }

  if (negd > 0 || negP > 0)
    printf("[Step14]: at t=%g; %i cells had d<0; %i cells had P<0\n",pG->time,negd,negP);
#endif /* FIRST_ORDER_FLUX_CORRECTION */

#ifdef STATIC_MESH_REFINEMENT
/*=== STEP 15: With SMR, store fluxes at fine/coarse boundaries ==============*/

/*--- Step 15a -----------------------------------------------------------------
 * store fluxes at boundaries of child grids.
 */

  for (ncg=0; ncg<pG->NCGrid; ncg++) {

/* x1-boundaries of child Grids (interior to THIS Grid) */

    for (dim=0; dim<2; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->CGrid[ncg].ijks[0];
        if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;
        jcs = pG->CGrid[ncg].ijks[1];
        jce = pG->CGrid[ncg].ijke[1];
        kcs = pG->CGrid[ncg].ijks[2];
        kce = pG->CGrid[ncg].ijke[2];

        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (j=jcs, jj=0; j<=jce; j++, jj++){
            pG->CGrid[ncg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            pG->CGrid[ncg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][kk][jj].B1c = 0.0;
            pG->CGrid[ncg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By;
            pG->CGrid[ncg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          for (j=jcs, jj=0; j<=jce; j++, jj++){
            pG->CGrid[ncg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          }
        }
        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (j=jcs, jj=0; j<=jce+1; j++, jj++){
            pG->CGrid[ncg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of child Grids (interior to THIS Grid) */

    for (dim=2; dim<4; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        ics = pG->CGrid[ncg].ijks[0];
        ice = pG->CGrid[ncg].ijke[0];
        if (dim==2) j = pG->CGrid[ncg].ijks[1];
        if (dim==3) j = pG->CGrid[ncg].ijke[1] + 1;
        kcs = pG->CGrid[ncg].ijks[2];
        kce = pG->CGrid[ncg].ijke[2];

        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            pG->CGrid[ncg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz;
            pG->CGrid[ncg].myFlx[dim][kk][ii].B2c = 0.0;
            pG->CGrid[ncg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          }
        }
        for (k=kcs, kk=0; k<=kce; k++, kk++){
          for (i=ics, ii=0; i<=ice+1; i++, ii++){
            pG->CGrid[ncg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x3-boundaries of child Grids (interior to THIS Grid) */

    for (dim=4; dim<6; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        ics = pG->CGrid[ncg].ijks[0];
        ice = pG->CGrid[ncg].ijke[0];
        jcs = pG->CGrid[ncg].ijks[1];
        jce = pG->CGrid[ncg].ijke[1];
        if (dim==4) k = pG->CGrid[ncg].ijks[2];
        if (dim==5) k = pG->CGrid[ncg].ijke[2] + 1;

        for (j=jcs, jj=0; j<=jce; j++, jj++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            pG->CGrid[ncg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx;
#ifndef BAROTROPIC
            pG->CGrid[ncg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->CGrid[ncg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By;
            pG->CGrid[ncg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz;
            pG->CGrid[ncg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->CGrid[ncg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (j=jcs, jj=0; j<=jce+1; j++, jj++){
          for (i=ics, ii=0; i<=ice; i++, ii++){
            pG->CGrid[ncg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          }
        }
        for (j=jcs, jj=0; j<=jce; j++, jj++){
          for (i=ics, ii=0; i<=ice+1; i++, ii++){
            pG->CGrid[ncg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          }
        }
#endif /* MHD */
      }
    }
  } /* end loop over child Grids */

/*--- Step 15b -----------------------------------------------------------------
 * store fluxes at boundaries with parent grids.
 */

  for (npg=0; npg<pG->NPGrid; npg++) {

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=0; dim<2; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->PGrid[npg].ijks[0];
        if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;
        jps = pG->PGrid[npg].ijks[1];
        jpe = pG->PGrid[npg].ijke[1];
        kps = pG->PGrid[npg].ijks[2];
        kpe = pG->PGrid[npg].ijke[2];

        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (j=jps, jj=0; j<=jpe; j++, jj++){
            pG->PGrid[npg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx;
            pG->PGrid[npg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            pG->PGrid[npg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][kk][jj].B1c = 0.0;
            pG->PGrid[npg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By;
            pG->PGrid[npg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          for (j=jps, jj=0; j<=jpe; j++, jj++){
            pG->PGrid[npg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          }
        }
        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (j=jps, jj=0; j<=jpe+1; j++, jj++){
            pG->PGrid[npg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=2; dim<4; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        ips = pG->PGrid[npg].ijks[0];
        ipe = pG->PGrid[npg].ijke[0];
        if (dim==2) j = pG->PGrid[npg].ijks[1];
        if (dim==3) j = pG->PGrid[npg].ijke[1] + 1;
        kps = pG->PGrid[npg].ijks[2];
        kpe = pG->PGrid[npg].ijke[2];

        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz;
            pG->PGrid[npg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            pG->PGrid[npg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz;
            pG->PGrid[npg].myFlx[dim][kk][ii].B2c = 0.0;
            pG->PGrid[npg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          }
        }
        for (k=kps, kk=0; k<=kpe; k++, kk++){
          for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            pG->PGrid[npg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          }
        }
#endif /* MHD */
      }
    }

/* x3-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=4; dim<6; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        ips = pG->PGrid[npg].ijks[0];
        ipe = pG->PGrid[npg].ijke[0];
        jps = pG->PGrid[npg].ijks[1];
        jpe = pG->PGrid[npg].ijke[1];
        if (dim==4) k = pG->PGrid[npg].ijks[2];
        if (dim==5) k = pG->PGrid[npg].ijke[2] + 1;

        for (j=jps, jj=0; j<=jpe; j++, jj++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d;
            pG->PGrid[npg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My;
            pG->PGrid[npg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            pG->PGrid[npg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx;
#ifndef BAROTROPIC
            pG->PGrid[npg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
            pG->PGrid[npg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By;
            pG->PGrid[npg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz;
            pG->PGrid[npg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++)
              pG->PGrid[npg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n];
#endif
          }
        }
#ifdef MHD
        for (j=jps, jj=0; j<=jpe+1; j++, jj++){
          for (i=ips, ii=0; i<=ipe; i++, ii++){
            pG->PGrid[npg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          }
        }
        for (j=jps, jj=0; j<=jpe; j++, jj++){
          for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            pG->PGrid[npg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          }
        }
#endif /* MHD */
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays */
void integrate_init_3d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

#ifdef MHD
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;

  if ((emf1_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* H_CORRECTION */
  if ((Wl_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x1Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x2Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wl_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;
  if ((Wr_x3Face=(Prim1DS***)calloc_3d_array(size3,size2,size1,sizeof(Prim1DS)))
    == NULL) goto on_error;

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef MHD
  if ((B1_x1Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
#endif /* MHD */

  if ((U1d = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur  = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W1d = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux = (Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
#ifdef MHD
  if ((emf1P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf2P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf3P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */

  if ((Uhalf = (ConsS***)calloc_3d_array(size3,size2,size1,sizeof(ConsS)))
    == NULL) goto on_error;

#ifdef SHEARING_BOX
  if ((Flxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((Flxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
#endif


  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_3d(void)
 *  \brief Free temporary integration arrays */
void integrate_destruct_3d(void)
{
#ifdef MHD
  if (emf1 != NULL) free_3d_array(emf1);
  if (emf2 != NULL) free_3d_array(emf2);
  if (emf3 != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_3d_array(eta1);
  if (eta2 != NULL) free_3d_array(eta2);
  if (eta3 != NULL) free_3d_array(eta3);
#endif /* H_CORRECTION */
  if (Wl_x1Face != NULL) free_3d_array(Wl_x1Face);
  if (Wr_x1Face != NULL) free_3d_array(Wr_x1Face);
  if (Wl_x2Face != NULL) free_3d_array(Wl_x2Face);
  if (Wr_x2Face != NULL) free_3d_array(Wr_x2Face);
  if (Wl_x3Face != NULL) free_3d_array(Wl_x3Face);
  if (Wr_x3Face != NULL) free_3d_array(Wr_x3Face);

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_3d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array(B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array(B3_x3Face);
#endif /* MHD */

  if (U1d != NULL) free(U1d);
  if (Ul  != NULL) free(Ul);
  if (Ur  != NULL) free(Ur);
  if (W1d != NULL) free(W1d);
  if (Wl  != NULL) free(Wl);
  if (Wr  != NULL) free(Wr);

  if (x1Flux  != NULL) free_3d_array(x1Flux);
  if (x2Flux  != NULL) free_3d_array(x2Flux);
  if (x3Flux  != NULL) free_3d_array(x3Flux);
#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxP != NULL) free_3d_array(x1FluxP);
  if (x2FluxP != NULL) free_3d_array(x2FluxP);
  if (x3FluxP != NULL) free_3d_array(x3FluxP);
#ifdef MHD
  if (emf1P   != NULL) free_3d_array(emf1P);
  if (emf2P   != NULL) free_3d_array(emf2P);
  if (emf3P   != NULL) free_3d_array(emf3P);
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */

  if (Uhalf  != NULL) free_3d_array(Uhalf);

#ifdef SHEARING_BOX
  if (Flxiib != NULL) free_2d_array(Flxiib);
  if (Flxoib != NULL) free_2d_array(Flxoib);
  if (rFlxiib != NULL) free_2d_array(rFlxiib);
  if (rFlxoib != NULL) free_2d_array(rFlxoib);
#endif


  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/

#ifdef MHD
/*! \fn static void integrate_emf1_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *   Note:
 *-  x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *-  x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *-  x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *-  x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *-  x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *-  x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */
static void integrate_emf1_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf2_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *   Note:
 *-  x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *-  x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *-  x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *-  x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *-  x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *-  x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */
static void integrate_emf2_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku+1; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu+1; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf3_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *   Note:
 *-  x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *-  x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *-  x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *-  x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *-  x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *-  x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX
 *- first_order_correction() - Added by N. Lemaster to run supersonic turbulence
 */
static void integrate_emf3_corner(const GridS *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku;
  Real de3_l1, de3_r1, de3_l2, de3_r2;
  Real rsf=1.0,lsf=1.0;

  il = pG->is-(nghost-1);   iu = pG->ie+(nghost-1);
  jl = pG->js-(nghost-1);   ju = pG->je+(nghost-1);
  kl = pG->ks-(nghost-1);   ku = pG->ke+(nghost-1);

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju+1; j++) {
      for (i=il; i<=iu+1; i++) {
#ifdef CYLINDRICAL
        rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf +
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf +
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}
#endif /* MHD */

#ifdef FIRST_ORDER_FLUX_CORRECTION
/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect ix)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect ix)
{
  int n;
#ifdef MHD
  int i,j,k;
#endif
  Real rsf=1.0,lsf=1.0;

/* Compute difference of predictor and corrector fluxes at cell faces */

  x1FD_i.d = x1Flux[ix.k][ix.j][ix.i].d - x1FluxP[ix.k][ix.j][ix.i].d;
  x2FD_j.d = x2Flux[ix.k][ix.j][ix.i].d - x2FluxP[ix.k][ix.j][ix.i].d;
  x3FD_k.d = x3Flux[ix.k][ix.j][ix.i].d - x3FluxP[ix.k][ix.j][ix.i].d;

  x1FD_ip1.d = x1Flux[ix.k][ix.j][ix.i+1].d - x1FluxP[ix.k][ix.j][ix.i+1].d;
  x2FD_jp1.d = x2Flux[ix.k][ix.j+1][ix.i].d - x2FluxP[ix.k][ix.j+1][ix.i].d;
  x3FD_kp1.d = x3Flux[ix.k+1][ix.j][ix.i].d - x3FluxP[ix.k+1][ix.j][ix.i].d;

  x1FD_i.Mx = x1Flux[ix.k][ix.j][ix.i].Mx - x1FluxP[ix.k][ix.j][ix.i].Mx;
  x2FD_j.Mx = x2Flux[ix.k][ix.j][ix.i].Mx - x2FluxP[ix.k][ix.j][ix.i].Mx;
  x3FD_k.Mx = x3Flux[ix.k][ix.j][ix.i].Mx - x3FluxP[ix.k][ix.j][ix.i].Mx;

  x1FD_ip1.Mx = x1Flux[ix.k][ix.j][ix.i+1].Mx - x1FluxP[ix.k][ix.j][ix.i+1].Mx;
  x2FD_jp1.Mx = x2Flux[ix.k][ix.j+1][ix.i].Mx - x2FluxP[ix.k][ix.j+1][ix.i].Mx;
  x3FD_kp1.Mx = x3Flux[ix.k+1][ix.j][ix.i].Mx - x3FluxP[ix.k+1][ix.j][ix.i].Mx;

  x1FD_i.My = x1Flux[ix.k][ix.j][ix.i].My - x1FluxP[ix.k][ix.j][ix.i].My;
  x2FD_j.My = x2Flux[ix.k][ix.j][ix.i].My - x2FluxP[ix.k][ix.j][ix.i].My;
  x3FD_k.My = x3Flux[ix.k][ix.j][ix.i].My - x3FluxP[ix.k][ix.j][ix.i].My;

  x1FD_ip1.My = x1Flux[ix.k][ix.j][ix.i+1].My - x1FluxP[ix.k][ix.j][ix.i+1].My;
  x2FD_jp1.My = x2Flux[ix.k][ix.j+1][ix.i].My - x2FluxP[ix.k][ix.j+1][ix.i].My;
  x3FD_kp1.My = x3Flux[ix.k+1][ix.j][ix.i].My - x3FluxP[ix.k+1][ix.j][ix.i].My;

  x1FD_i.Mz = x1Flux[ix.k][ix.j][ix.i].Mz - x1FluxP[ix.k][ix.j][ix.i].Mz;
  x2FD_j.Mz = x2Flux[ix.k][ix.j][ix.i].Mz - x2FluxP[ix.k][ix.j][ix.i].Mz;
  x3FD_k.Mz = x3Flux[ix.k][ix.j][ix.i].Mz - x3FluxP[ix.k][ix.j][ix.i].Mz;

  x1FD_ip1.Mz = x1Flux[ix.k][ix.j][ix.i+1].Mz - x1FluxP[ix.k][ix.j][ix.i+1].Mz;
  x2FD_jp1.Mz = x2Flux[ix.k][ix.j+1][ix.i].Mz - x2FluxP[ix.k][ix.j+1][ix.i].Mz;
  x3FD_kp1.Mz = x3Flux[ix.k+1][ix.j][ix.i].Mz - x3FluxP[ix.k+1][ix.j][ix.i].Mz;

#ifndef BAROTROPIC
  x1FD_i.E = x1Flux[ix.k][ix.j][ix.i].E - x1FluxP[ix.k][ix.j][ix.i].E;
  x2FD_j.E = x2Flux[ix.k][ix.j][ix.i].E - x2FluxP[ix.k][ix.j][ix.i].E;
  x3FD_k.E = x3Flux[ix.k][ix.j][ix.i].E - x3FluxP[ix.k][ix.j][ix.i].E;

  x1FD_ip1.E = x1Flux[ix.k][ix.j][ix.i+1].E - x1FluxP[ix.k][ix.j][ix.i+1].E;
  x2FD_jp1.E = x2Flux[ix.k][ix.j+1][ix.i].E - x2FluxP[ix.k][ix.j+1][ix.i].E;
  x3FD_kp1.E = x3Flux[ix.k+1][ix.j][ix.i].E - x3FluxP[ix.k+1][ix.j][ix.i].E;
#endif /* BAROTROPIC */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++){
    x1FD_i.s[n] = x1Flux[ix.k][ix.j][ix.i].s[n] - x1FluxP[ix.k][ix.j][ix.i].s[n];
    x2FD_j.s[n] = x2Flux[ix.k][ix.j][ix.i].s[n] - x2FluxP[ix.k][ix.j][ix.i].s[n];
    x3FD_k.s[n] = x3Flux[ix.k][ix.j][ix.i].s[n] - x3FluxP[ix.k][ix.j][ix.i].s[n];
    x1FD_ip1.s[n] = x1Flux[ix.k][ix.j][ix.i+1].s[n] - x1FluxP[ix.k][ix.j][ix.i+1].s[n];
    x2FD_jp1.s[n] = x2Flux[ix.k][ix.j+1][ix.i].s[n] - x2FluxP[ix.k][ix.j+1][ix.i].s[n];
    x3FD_kp1.s[n] = x3Flux[ix.k+1][ix.j][ix.i].s[n] - x3FluxP[ix.k+1][ix.j][ix.i].s[n];
  }
#endif

#ifdef MHD
  emf1D_kj     = emf1[ix.k  ][ix.j  ][ix.i] - emf1P[ix.k  ][ix.j  ][ix.i];
  emf1D_kjp1   = emf1[ix.k  ][ix.j+1][ix.i] - emf1P[ix.k  ][ix.j+1][ix.i];
  emf1D_kp1j   = emf1[ix.k+1][ix.j  ][ix.i] - emf1P[ix.k+1][ix.j  ][ix.i];
  emf1D_kp1jp1 = emf1[ix.k+1][ix.j+1][ix.i] - emf1P[ix.k+1][ix.j+1][ix.i];

  emf2D_ki     = emf2[ix.k  ][ix.j][ix.i  ] - emf2P[ix.k  ][ix.j][ix.i  ];
  emf2D_kip1   = emf2[ix.k  ][ix.j][ix.i+1] - emf2P[ix.k  ][ix.j][ix.i+1];
  emf2D_kp1i   = emf2[ix.k+1][ix.j][ix.i  ] - emf2P[ix.k+1][ix.j][ix.i  ];
  emf2D_kp1ip1 = emf2[ix.k+1][ix.j][ix.i+1] - emf2P[ix.k+1][ix.j][ix.i+1];

  emf3D_ji     = emf3[ix.k][ix.j  ][ix.i  ] - emf3P[ix.k][ix.j  ][ix.i  ];
  emf3D_jip1   = emf3[ix.k][ix.j  ][ix.i+1] - emf3P[ix.k][ix.j  ][ix.i+1];
  emf3D_jp1i   = emf3[ix.k][ix.j+1][ix.i  ] - emf3P[ix.k][ix.j+1][ix.i  ];
  emf3D_jp1ip1 = emf3[ix.k][ix.j+1][ix.i+1] - emf3P[ix.k][ix.j+1][ix.i+1];
#endif /* MHD */

/* Use flux differences to correct bad cell-centered quantities */
  ApplyCorr(pG,ix.i,ix.j,ix.k,1,1,1,1,1,1);

#ifdef SELF_GRAVITY
/* Save mass fluxes in Grid structure for source term correction in main loop */
  pG->x1MassFlux[ix.k  ][ix.j  ][ix.i  ] = x1FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x2MassFlux[ix.k  ][ix.j  ][ix.i  ] = x2FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x3MassFlux[ix.k  ][ix.j  ][ix.i  ] = x3FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x1MassFlux[ix.k  ][ix.j  ][ix.i+1] = x1FluxP[ix.k  ][ix.j  ][ix.i+1].d;
  pG->x2MassFlux[ix.k  ][ix.j+1][ix.i  ] = x2FluxP[ix.k  ][ix.j+1][ix.i  ].d;
  pG->x3MassFlux[ix.k+1][ix.j  ][ix.i  ] = x3FluxP[ix.k+1][ix.j  ][ix.i  ].d;
#endif /* SELF_GRAVITY */


/* Use flux differences to correct cell-centered values at i-1 and i+1 */
  
  if (ix.i > pG->is) ApplyCorr(pG,ix.i-1,ix.j,ix.k,-1,0,0,0,0,0);
  if (ix.i < pG->ie) ApplyCorr(pG,ix.i+1,ix.j,ix.k,0,-1,0,0,0,0);

/* Use flux differences to correct cell-centered values at j-1 and j+1 */
  
  if (ix.j > pG->js) ApplyCorr(pG,ix.i,ix.j-1,ix.k,0,0,-1,0,0,0);
  if (ix.j < pG->je) ApplyCorr(pG,ix.i,ix.j+1,ix.k,0,0,0,-1,0,0);

/* Use flux differences to correct cell-centered values at k-1 and k+1 */
 
  if (ix.k > pG->ks) ApplyCorr(pG,ix.i,ix.j,ix.k-1,0,0,0,0,-1,0);
  if (ix.k < pG->ke) ApplyCorr(pG,ix.i,ix.j,ix.k+1,0,0,0,0,0,-1);

#ifdef MHD
/* Compute new cell-centered fields */
  for (k=(ix.k-1); k<=(ix.k+1); k++) {
  for (j=(ix.j-1); j<=(ix.j+1); j++) {
  for (i=(ix.i-1); i<=(ix.i+1); i++) {
#ifdef CYLINDRICAL
    rsf = pG->ri[i+1]/pG->r[i];  lsf = pG->ri[i]/pG->r[i];
#endif
    pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
    pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
    pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
  }}}
#endif /* MHD */

// /* With SMR, replace higher-order fluxes with predict fluxes in case they are
//  * used at fine/coarse grid boundaries */
// #ifdef STATIC_MESH_REFINEMENT
//   x1Flux[ix.k][ix.j][ix.i] = x1FluxP[ix.k][ix.j][ix.i];
//   x2Flux[ix.k][ix.j][ix.i] = x2FluxP[ix.k][ix.j][ix.i];
//   x3Flux[ix.k][ix.j][ix.i] = x3FluxP[ix.k][ix.j][ix.i];
// #endif /* STATIC_MESH_REFINEMENT */

/* Must replace higher-order fluxes/emfs with predict fluxes/emfs in case two
 * adjacent cells are corrected.  Otherwise, flux/emf differences may get
 * corrected more than once.  By overwriting the higher-order fluxes/emfs,
 * the second time through, the differences will be zero. */
  x1Flux[ix.k  ][ix.j  ][ix.i  ] = x1FluxP[ix.k  ][ix.j  ][ix.i  ];
  x2Flux[ix.k  ][ix.j  ][ix.i  ] = x2FluxP[ix.k  ][ix.j  ][ix.i  ];
  x3Flux[ix.k  ][ix.j  ][ix.i  ] = x3FluxP[ix.k  ][ix.j  ][ix.i  ];
  x1Flux[ix.k  ][ix.j  ][ix.i+1] = x1FluxP[ix.k  ][ix.j  ][ix.i+1];
  x2Flux[ix.k  ][ix.j+1][ix.i  ] = x2FluxP[ix.k  ][ix.j+1][ix.i  ];
  x3Flux[ix.k+1][ix.j  ][ix.i  ] = x3FluxP[ix.k+1][ix.j  ][ix.i  ];
#ifdef MHD
  emf1[ix.k  ][ix.j  ][ix.i  ] = emf1P[ix.k  ][ix.j  ][ix.i  ];
  emf1[ix.k  ][ix.j+1][ix.i  ] = emf1P[ix.k  ][ix.j+1][ix.i  ];
  emf1[ix.k+1][ix.j  ][ix.i  ] = emf1P[ix.k+1][ix.j  ][ix.i  ];
  emf1[ix.k+1][ix.j+1][ix.i  ] = emf1P[ix.k+1][ix.j+1][ix.i  ];
  emf2[ix.k  ][ix.j  ][ix.i  ] = emf2P[ix.k  ][ix.j  ][ix.i  ];
  emf2[ix.k  ][ix.j  ][ix.i+1] = emf2P[ix.k  ][ix.j  ][ix.i+1];
  emf2[ix.k+1][ix.j  ][ix.i  ] = emf2P[ix.k+1][ix.j  ][ix.i  ];
  emf2[ix.k+1][ix.j  ][ix.i+1] = emf2P[ix.k+1][ix.j  ][ix.i+1];
  emf3[ix.k  ][ix.j  ][ix.i  ] = emf3P[ix.k  ][ix.j  ][ix.i  ];
  emf3[ix.k  ][ix.j  ][ix.i+1] = emf3P[ix.k  ][ix.j  ][ix.i+1];
  emf3[ix.k  ][ix.j+1][ix.i  ] = emf3P[ix.k  ][ix.j+1][ix.i  ];
  emf3[ix.k  ][ix.j+1][ix.i+1] = emf3P[ix.k  ][ix.j+1][ix.i+1];
#endif
}

static void ApplyCorr(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3)
{
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real x1,x2,x3,phil,phir,phic;
  int n;
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */
  Real lsf=1.0,rsf=1.0;

#ifdef CYLINDRICAL
  rsf = (rx1>0) ? pG->ri[i+1]/pG->r[i] : pG->ri[i  ]/pG->r[i];
  lsf = (lx1>0) ? pG->ri[i  ]/pG->r[i] : pG->ri[i+1]/pG->r[i];
  dtodx2 = pG->dt/(pG->r[i]*pG->dx2);
#endif
  pG->U[k][j][i].d  += dtodx1*(rsf*rx1*x1FD_ip1.d  - lsf*lx1*x1FD_i.d );
  pG->U[k][j][i].M1 += dtodx1*(rsf*rx1*x1FD_ip1.Mx - lsf*lx1*x1FD_i.Mx);
  pG->U[k][j][i].M2 += dtodx1*(SQR(rsf)*rx1*x1FD_ip1.My - SQR(lsf)*lx1*x1FD_i.My);
  pG->U[k][j][i].M3 += dtodx1*(rsf*rx1*x1FD_ip1.Mz - lsf*lx1*x1FD_i.Mz);
#ifndef BAROTROPIC
  pG->U[k][j][i].E  += dtodx1*(rsf*rx1*x1FD_ip1.E  - lsf*lx1*x1FD_i.E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++)
    pG->U[k][j][i].s[n] += dtodx1*(rsf*rx1*x1FD_ip1.s[n] - lsf*lx1*x1FD_i.s[n]);
#endif

  pG->U[k][j][i].d  += dtodx2*(rx2*x2FD_jp1.d  - lx2*x2FD_j.d );
  pG->U[k][j][i].M1 += dtodx2*(rx2*x2FD_jp1.Mz - lx2*x2FD_j.Mz);
  pG->U[k][j][i].M2 += dtodx2*(rx2*x2FD_jp1.Mx - lx2*x2FD_j.Mx);
  pG->U[k][j][i].M3 += dtodx2*(rx2*x2FD_jp1.My - lx2*x2FD_j.My);
#ifndef BAROTROPIC                                   
  pG->U[k][j][i].E  += dtodx2*(rx2*x2FD_jp1.E  - lx2*x2FD_j.E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++)
    pG->U[k][j][i].s[n] += dtodx2*(rx2*x2FD_jp1.s[n] - lx2*x2FD_j.s[n]);
#endif

  pG->U[k][j][i].d  += dtodx3*(rx3*x3FD_kp1.d  - lx3*x3FD_k.d );
  pG->U[k][j][i].M1 += dtodx3*(rx3*x3FD_kp1.My - lx3*x3FD_k.My);
  pG->U[k][j][i].M2 += dtodx3*(rx3*x3FD_kp1.Mz - lx3*x3FD_k.Mz);
  pG->U[k][j][i].M3 += dtodx3*(rx3*x3FD_kp1.Mx - lx3*x3FD_k.Mx);
#ifndef BAROTROPIC                                   
  pG->U[k][j][i].E  += dtodx3*(rx3*x3FD_kp1.E  - lx3*x3FD_k.E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++)
    pG->U[k][j][i].s[n] += dtodx3*(rx3*x3FD_kp1.s[n] - lx3*x3FD_k.s[n]);
#endif

#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

/* Calculate the flux for the y-momentum fluctuation */
  frx1_dM2 = x1FD_ip1.My;
  flx1_dM2 = x1FD_i.My;
  frx2_dM2 = x2FD_jp1.Mx;
  flx2_dM2 = x2FD_j.Mx;
  frx3_dM2 = x3FD_kp1.Mz;
  flx3_dM2 = x3FD_k.Mz;
#ifndef FARGO
  frx1_dM2 += qom*(x1+rx1*0.5*pG->dx1)*x1FD_ip1.d;
  flx1_dM2 += qom*(x1-lx1*0.5*pG->dx1)*x1FD_i.d;
  frx2_dM2 += qom*(x1            )*x2FD_jp1.d;
  flx2_dM2 += qom*(x1            )*x2FD_j.d;
  frx3_dM2 += qom*(x1            )*x3FD_kp1.d;
  flx3_dM2 += qom*(x1            )*x3FD_k.d;
#endif

/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
  M1e = q1*(rx1*x1FD_ip1.Mx - lx1*x1FD_i.Mx)
      + q2*(rx2*x2FD_jp1.Mz - lx2*x2FD_j.Mz)
      + q3*(rx3*x3FD_kp1.My - lx3*x3FD_k.My);

  dM2e = q1*(rx1*frx1_dM2 - lx1*flx1_dM2)
       + q2*(rx2*frx2_dM2 - lx2*flx2_dM2) 
       + q3*(rx3*frx3_dM2 - lx3*flx3_dM2);

/* Update the 1- and 2-momentum for the Coriolis and tidal
 * potential momentum source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

  pG->U[k][j][i].M1 += (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
  pG->U[k][j][i].M2 += 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;

#ifndef FARGO
  pG->U[k][j][i].M2 += 0.5*qshear*om_dt*(ABS(lx1)*x1FD_i.d + ABS(rx1)*x1FD_ip1.d);
#endif

/* Update the energy for a fixed potential.
 * This update is identical to non-SHEARING_BOX below  */

  phic = (*ShearingBoxPot)((x1            ),x2,x3);
  phir = (*ShearingBoxPot)((x1+rx1*0.5*pG->dx1),x2,x3);
  phil = (*ShearingBoxPot)((x1-lx1*0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx1*(lx1*x1FD_i.d*(phic - phil) +
                              rx1*x1FD_ip1.d*(phir - phic));
#endif

  phir = (*ShearingBoxPot)(x1,(x2+rx2*0.5*pG->dx2),x3);
  phil = (*ShearingBoxPot)(x1,(x2-lx2*0.5*pG->dx2),x3);
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                              rx2*x2FD_jp1.d*(phir - phic));
#endif

  phir = (*ShearingBoxPot)(x1,x2,(x3+rx3*0.5*pG->dx3));
  phil = (*ShearingBoxPot)(x1,x2,(x3-lx3*0.5*pG->dx3));
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                              rx3*x3FD_kp1.d*(phir - phic));
#endif
#endif /* SHEARING_BOX */

  if (StaticGravPot != NULL){
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    phic = (*StaticGravPot)(x1,x2,x3);
    phir = (*StaticGravPot)((x1+rx1*0.5*pG->dx1),x2,x3);
    phil = (*StaticGravPot)((x1-lx1*0.5*pG->dx1),x2,x3);

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx1*(lsf*lx1*x1FD_i.d*(phic - phil) +
                                rsf*rx1*x1FD_ip1.d*(phir - phic));
#endif

    phir = (*StaticGravPot)(x1,(x2+rx2*0.5*pG->dx2),x3);
    phil = (*StaticGravPot)(x1,(x2-lx2*0.5*pG->dx2),x3);

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                                rx2*x2FD_jp1.d*(phir - phic));
#endif

    phir = (*StaticGravPot)(x1,x2,(x3+rx3*0.5*pG->dx3));
    phil = (*StaticGravPot)(x1,x2,(x3-lx3*0.5*pG->dx3));

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                                rx3*x3FD_kp1.d*(phir - phic));
#endif
  }

#ifdef SELF_GRAVITY
  phic = pG->Phi[k][j][i];
  phil = 0.5*(pG->Phi[k][j][i-lx1] + pG->Phi[k][j][i  ]);
  phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+rx1]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx1*(lx1*x1FD_i.d*(phic - phil) +
                              rx1*x1FD_ip1.d*(phir - phic));
#endif

  phil = 0.5*(pG->Phi[k][j-lx2][i] + pG->Phi[k][j  ][i]);
  phir = 0.5*(pG->Phi[k][j  ][i] + pG->Phi[k][j+rx2][i]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                              rx2*x2FD_jp1.d*(phir - phic));
#endif

  phil = 0.5*(pG->Phi[k-lx3][j][i] + pG->Phi[k  ][j][i]);
  phir = 0.5*(pG->Phi[k  ][j][i] + pG->Phi[k+rx3][j][i]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                              rx3*x3FD_kp1.d*(phir - phic));
#endif

#endif /* SELF_GRAVITY */


/* Use emf differences to correct face-centered fields on edges of bad cell */
#ifdef MHD
#ifdef CYLINDRICAL
  dtodx2 = pG->dt/(pG->ri[i]*pG->dx2);
#endif
  pG->B1i[k][j][i  ] -= dtodx3*(rx3*emf2D_kp1i   - lx3*emf2D_ki) -
                        dtodx2*(rx2*emf3D_jp1i   - lx2*emf3D_ji);
#ifdef CYLINDRICAL
  dtodx2 = pG->dt/(pG->ri[i+1]*pG->dx2);
#endif
  pG->B1i[k][j][i+1] -= dtodx3*(rx3*emf2D_kp1ip1 - lx3*emf2D_kip1) -
                        dtodx2*(rx2*emf3D_jp1ip1 - lx2*emf3D_jip1);
  pG->B2i[k][j  ][i] -= dtodx1*(rx1*emf3D_jip1   - lx1*emf3D_ji) -
                        dtodx3*(rx3*emf1D_kp1j   - lx3*emf1D_kj);
  pG->B2i[k][j+1][i] -= dtodx1*(rx1*emf3D_jp1ip1 - lx1*emf3D_jp1i) -
                        dtodx3*(rx3*emf1D_kp1jp1 - lx3*emf1D_kjp1);
#ifdef CYLINDRICAL
  rsf = (rx1>0) ? pG->ri[i+1]/pG->r[i] : pG->ri[i  ]/pG->r[i];
  lsf = (lx1>0) ? pG->ri[i  ]/pG->r[i] : pG->ri[i+1]/pG->r[i];
  dtodx2 = pG->dt/(pG->r[i]*pG->dx2);
#endif
  pG->B3i[k  ][j][i] -= dtodx2*(rx2*emf1D_kjp1   - lx2*emf1D_kj) -
                        dtodx1*(rsf*rx1*emf2D_kip1   - lsf*lx1*emf2D_ki);
  pG->B3i[k+1][j][i] -= dtodx2*(rx2*emf1D_kp1jp1 - lx2*emf1D_kp1j) -
                        dtodx1*(rsf*rx1*emf2D_kp1ip1 - lsf*lx1*emf2D_kp1i);
#endif 

}

#endif /* FIRST_ORDER_FLUX_CORRECTION */

#endif /* VL_INTEGRATOR */

#endif /* SPECIAL_RELATIVITY */
