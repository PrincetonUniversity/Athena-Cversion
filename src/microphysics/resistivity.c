#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Adds explicit resistivity terms to the induction and energy eqns,
 *      dB/dt = -Curl(E)
 *      dE/dt = Div(B X E)
 *   where E = eta_Ohm J + eta_Hall(J X B)/B + eta_AD J_perp = (emf) 
 *         J = Curl(B) = current density
 *         eta_Ohm = Ohmic resistivity
 *         eta_Hall = Hall diffusion coefficient
 *         eta_AD = ambipolar diffusion coefficient
 *   The induction equation is updated using CT to keep div(B)=0.  The total
 *   electric field (resistive EMF) is computed from calls to the EField_*
 *   functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  resistivity() - updates induction and energy eqns with resistive term.
 *  resistivity_init() - allocates memory needed
 *  resistivity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef RESISTIVITY

#ifdef HYDRO
#error : resistivity only works for MHD.
#endif /* HYDRO */

/* current, emf, and energy flux, contained in 3D vector structure */
Real3Vect ***J=NULL, ***emf=NULL, ***EnerFlux=NULL;

/* emf and intermediate B and J for Hall MHD */
static Real3Vect ***emfh=NULL, ***Bcor=NULL, ***Jcor=NULL;

/* for 3D shearing box, variables needed to conserve net Bz */
#ifdef SHEARING_BOX
static Real ***emf2=NULL;
static Real **remapEyiib=NULL, **remapEyoib=NULL;
static Real ***J2=NULL;
static Real ***remapJyiib=NULL,***remapJyoib=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   EField_Ohm  - computes electric field due to Ohmic dissipation
 *   EField_Hall - computes electric field due to Hall effect
 *   EField_AD   - computes electric field due to ambipolar diffusion
 *   hyper_diffusion? - add hyper-resistivity to help stabilize the Hall term
 *============================================================================*/

void EField_Ohm(DomainS *pD);
void EField_Hall(DomainS *pD);
void EField_AD(DomainS *pD);

void hyper_diffusion4(DomainS *pD, Real prefac);
void hyper_diffusion6(DomainS *pD, Real prefac);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* resistivity:
 */

void resistivity(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;
#ifdef SHEARING_BOX
  int my_iproc, my_jproc, my_kproc;
  int nlayer=1;
#endif
  int ndim=1;
#ifdef STS
  Real my_dt = STS_dt;
#else
  Real my_dt = pG->dt;
#endif
  Real dtodx1 = my_dt/pG->dx1, dtodx2 = 0.0, dtodx3 = 0.0;

#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif
  Real dx1i=1.0/pG->dx1, dx2i=0.0, dx3i=0.0;
  Real lsf=1.0,rsf=1.0;

  if (pG->Nx[1] > 1){
    jl = js - 4;
    ju = je + 4;
    dtodx2 = my_dt/pG->dx2;
    dx2i=1.0/pG->dx2;
    ndim++;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 4;
    ku = ke + 4;
    dtodx3 = my_dt/pG->dx3;
    dx3i=1.0/pG->dx3;
    ndim++;
  } else {
    kl = ks;
    ku = ke;
  }

/* zero fluxes (electric fields) */

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
    for (i=is-4; i<=ie+4; i++) {
      emf[k][j][i].x1 = 0.0;
      emf[k][j][i].x2 = 0.0;
      emf[k][j][i].x3 = 0.0;
    }
  }}

  if (Q_Hall > 0.0) {
    for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-4; i<=ie+4; i++) {
        emfh[k][j][i].x1 = 0.0;
        emfh[k][j][i].x2 = 0.0;
        emfh[k][j][i].x3 = 0.0;
        Bcor[k][j][i].x1 = pG->B1i[k][j][i];
        Bcor[k][j][i].x2 = pG->B2i[k][j][i];
        Bcor[k][j][i].x3 = pG->B3i[k][j][i];
      }
    }}
  }

#ifdef SHEARING_BOX
  if (Q_AD > 0.0) {
    nlayer = 2;
  }
  if (Q_Hall > 0.0) {
    nlayer = 4;
  }
#endif

/*--- Step 1. Compute currents.------------------------------------------------
 * Note:  J1 = (dB3/dx2 - dB2/dx3)
 *        J2 = (dB1/dx3 - dB3/dx1)
 *        J3 = (dB2/dx1 - dB1/dx2) */

/* 1D PROBLEM */
  if (ndim == 1){
    for (i=is-3; i<=ie+4; i++) {
      J[ks][js][i].x1 = 0.0;
      J[ks][js][i].x2 = -(pG->U[ks][js][i].B3c - pG->U[ks][js][i-1].B3c)/pG->dx1;
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      J[ks][js][i].x3 =  (rsf*pG->U[ks][js][i].B2c - lsf*pG->U[ks][js][i-1].B2c)/pG->dx1;
    }
    J[ks][js][0].x1 = 0.0;
  }

/* 2D PROBLEM */
  if (ndim == 2){
    for (j=js-3; j<=je+4; j++) {
      for (i=is-3; i<=ie+4; i++) {
#ifdef CYLINDRICAL
        dx2i = 1.0/(r[i]*pG->dx2);
#endif
        J[ks][j][i].x1=  dx2i*(pG->U[ks][j][i].B3c - pG->U[ks][j-1][i  ].B3c);
        J[ks][j][i].x2= -dx1i*(pG->U[ks][j][i].B3c - pG->U[ks][j  ][i-1].B3c);
#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
        dx2i = 1.0/(ri[i]*pG->dx2);
#endif
        J[ks][j][i].x3=  dx1i*(rsf*pG->B2i[ks][j][i] - lsf*pG->B2i[ks][j  ][i-1]) -
                         dx2i*(    pG->B1i[ks][j][i] -     pG->B1i[ks][j-1][i  ]);
      }
      i = is-4;
#ifdef CYLINDRICAL
      dx2i = 1.0/(r[i]*pG->dx2);
#endif
      J[ks][j][i].x1=  dx2i*(pG->U[ks][j][i].B3c - pG->U[ks][j-1][i  ].B3c);
    }
    j = js-4;
    for (i=is-3; i<=ie+4; i++) {
      J[ks][j][i].x2= -dx1i*(pG->U[ks][j][i].B3c - pG->U[ks][j  ][i-1].B3c);
    }
  }

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=ks-3; k<=ke+4; k++) {
     for (j=js-3; j<=je+4; j++) {
      for (i=is-3; i<=ie+4; i++) {
#ifdef CYLINDRICAL
        dx2i=1.0/(r[i]*pG->dx2);
#endif
        J[k][j][i].x1 = dx2i*(pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ]) -
                        dx3i*(pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ]);
        J[k][j][i].x2 = dx3i*(pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ]) -
                        dx1i*(pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1]);
#ifdef CYLINDRICAL
        dx2i=1.0/(ri[i]*pG->dx2);
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        J[k][j][i].x3 = dx1i*(rsf*pG->B2i[k][j][i] - lsf*pG->B2i[k  ][j  ][i-1]) -
                        dx2i*(    pG->B1i[k][j][i] -     pG->B1i[k  ][j-1][i  ]);
      }
      i = is-4;
#ifdef CYLINDRICAL
      dx2i=1.0/(r[i]*pG->dx2);
#endif
      J[k][j][i].x1 = dx2i*(pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ]) -
                      dx3i*(pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ]);
     }
     j = js-4;
     for (i=is-3; i<=ie+4; i++) {
        J[k][j][i].x2 = dx3i*(pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ]) -
                        dx1i*(pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1]);
     }
    }
    k = ks-4;
    for (j=js-3; j<=je+4; j++) {
    for (i=is-3; i<=ie+4; i++) {
#ifdef CYLINDRICAL
      dx2i=1.0/(ri[i]*pG->dx2);
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      J[k][j][i].x3 = dx1i*(rsf*pG->B2i[k][j][i] - lsf*pG->B2i[k  ][j  ][i-1]) -
                      dx2i*(    pG->B1i[k][j][i] -     pG->B1i[k  ][j-1][i  ]);
   }}
  }

/* Remap Jy at the inner and outer shearing box boundaries */
#ifdef SHEARING_BOX
  if (pG->Nx[2] > 1){

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* compute remapped Jy from opposite side of grid */
    for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
    for(j=js; j<=je; j++) {
    for(i=is; i<=ie+1; i++) {
      J2[k][j][i]   = J[k][j][i].x2;
    }}}

    if (my_iproc == 0) {
      RemapJy_ix1(pD, J2, remapJyiib, nlayer);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapJy_ox1(pD, J2, remapJyoib, nlayer);
    }

/* Now average Jy and remapped Jy */
    if (my_iproc == 0) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=is-nlayer+1; i<=is; i++) {
        J[k][j][i].x2  = 0.5*(J[k][j][i].x2 + remapJyiib[i-is+nlayer-1][k][j]);
      }}}
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=ie+1; i<=ie+nlayer; i++) {
        J[k][j][i].x2 = 0.5*(J[k][j][i].x2 + remapJyoib[i-ie-1][k][j]);
      }}}
    }
  }
#endif /* SHEARING_BOX */


/*--- Step 2.  Call functions to compute resistive EMFs ------------------------
 * including Ohmic dissipation, the Hall effect, and ambipolar diffusion.
 * Current density (J) and emfs are global variables in this file. */

  if (eta_Ohm > 0.0) EField_Ohm(pD);
  if (Q_Hall > 0.0)  EField_Hall(pD);
  if (Q_AD > 0.0)    EField_AD(pD);

/* Remap Ey at is and ie+1 to conserve Bz in shearing box */
#ifdef SHEARING_BOX
  if (pG->Nx[2] > 1){

/* compute remapped Ey from opposite side of grid */
    for(k=ks; k<=ke+1; k++) {
    for(j=js; j<=je; j++)   {
      emf2[k][j][is]   = emf[k][j][is].x2;
      emf2[k][j][ie+1] = emf[k][j][ie+1].x2;
    }}

    if (my_iproc == 0) {
      RemapEy_ix1(pD, emf2, remapEyiib);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapEy_ox1(pD, emf2, remapEyoib);
    }

/* Now average Ey and remapped Ey */
    if (my_iproc == 0) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++) {
          emf[k][j][is].x2 = 0.5*(emf2[k][j][is] + remapEyiib[k][j]);
    }}}

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++) {
          emf[k][j][ie+1].x2 = 0.5*(emf2[k][j][ie+1] + remapEyoib[k][j]);
    }}}

  }
#endif /* SHEARING_BOX */

#ifndef BAROTROPIC
/*--- Step 3.  Compute energy fluxes -------------------------------------------
 * flux of total energy due to resistive diffusion = B X emf
 *  EnerFlux.x =  By*emf.z - Bz*emf.y
 *  EnerFlux.y =  Bz*emf.x - Bx*emf.z
 *  EnerFlux.z =  Bx*emf.y - By*emf.x
 */

/* 1D PROBLEM */
  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      EnerFlux[ks][js][i].x1 =
         0.5*(rsf*pG->U[ks][js][i].B2c + lsf*pG->U[ks][js][i-1].B2c)*emf[ks][js][i].x3
       - 0.5*(rsf*pG->U[ks][js][i].B3c + lsf*pG->U[ks][js][i-1].B3c)*emf[ks][js][i].x2;
    }
  } 
      
/* 2D PROBLEM */
  if (ndim == 2){
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      EnerFlux[ks][j][i].x1 = 0.25*(rsf*pG->U[ks][j][i].B2c + lsf*pG->U[ks][j][i-1].B2c)*
                            (emf[ks][j][i].x3 + emf[ks][j+1][i].x3)
         - 0.5*(rsf*pG->U[ks][j][i].B3c + lsf*pG->U[ks][j][i-1].B3c)*emf[ks][j][i].x2;
    }}
    
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      EnerFlux[ks][j][i].x2 =
         0.5*(rsf*pG->U[ks][j][i].B3c + lsf*pG->U[ks][j-1][i].B3c)*emf[ks][j][i].x1;
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      EnerFlux[ks][j][i].x2 -=
         0.25*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c)*
                (lsf*emf[ks][j][i].x3 + rsf*emf[ks][j][i+1].x3);
    }}
  }   

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        EnerFlux[k][j][i].x1 = 0.25*(rsf*pG->U[k][j][i].B2c + lsf*pG->U[k][j][i-1].B2c)*
                             (emf[k][j][i].x3 + emf[k][j+1][i].x3)
                            - 0.25*(rsf*pG->U[k][j][i].B3c + lsf*pG->U[k][j][i-1].B3c)*
                             (emf[k][j][i].x2 + emf[k+1][j][i].x2);
      }
    }}

    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        EnerFlux[k][j][i].x2 = 0.25*(rsf*pG->U[k][j][i].B3c + lsf*pG->U[k][j-1][i].B3c)*
                             (emf[k][j][i].x1 + emf[k+1][j][i].x1);
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        EnerFlux[k][j][i].x2-= 0.25*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c)*
                             (lsf*emf[k][j][i].x3 + rsf*emf[k][j][i+1].x3);
      }
    }}

    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        EnerFlux[k][j][i].x3 = 0.25*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c)*
                             (lsf*emf[k][j][i].x2 + emf[k][j][i+1].x2)
                            - 0.25*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c)*
                             (lsf*emf[k][j][i].x1 + rsf*emf[k][j+1][i].x1);
      }
    }}
  }

/*--- Step 4.  Update total energy ---------------------------------------------
 * Update energy using x1-fluxes */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[k][j][i].E += dtodx1*(rsf*EnerFlux[k][j][i+1].x1 - lsf*EnerFlux[k][j][i].x1);
    }
  }}

/* Update energy using x2-fluxes */

  if (pG->Nx[1] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        dtodx2 = my_dt/(r[i]*pG->dx2);
#endif
        pG->U[k][j][i].E += dtodx2*(EnerFlux[k][j+1][i].x2 -EnerFlux[k][j][i].x2);
      }
    }}
  }

/* Update energy using x3-fluxes */

  if (pG->Nx[2] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(EnerFlux[k+1][j][i].x3 -EnerFlux[k][j][i].x3);
      }
    }}
  }
#endif /* BAROTROPIC */

/*--- Step 5. CT update of magnetic field -------------------------------------
 * using total resistive EMFs.  This is identical
 * to the CT update in the integrators: dB/dt = -Curl(E) */

/* 1D PROBLEM: centered differences for B2c and B3c */
  if (ndim == 1){
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[ks][js][i].B2c += dtodx1*(    emf[ks][js][i+1].x3 -     emf[ks][js][i].x3);
      pG->U[ks][js][i].B3c -= dtodx1*(rsf*emf[ks][js][i+1].x2 - lsf*emf[ks][js][i].x2);
/* For consistency, set B2i and B3i to cell-centered values. */
      pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
      pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
    }
  }

/* 2D PROBLEM: CT +  centered differences for B3c */
  if (ndim == 2){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        dtodx2 = my_dt/(ri[i]*pG->dx2);
#endif
        pG->B1i[ks][j][i] -= dtodx2*(emf[ks][j+1][i  ].x3 - emf[ks][j][i].x3);
        pG->B2i[ks][j][i] += dtodx1*(emf[ks][j  ][i+1].x3 - emf[ks][j][i].x3);

#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
        dtodx2 = my_dt/(r[i]*pG->dx2);
#endif
        pG->U[ks][j][i].B3c += dtodx2*(    emf[ks][j+1][i  ].x1 -     emf[ks][j][i].x1) -
                               dtodx1*(rsf*emf[ks][j  ][i+1].x2 - lsf*emf[ks][j][i].x2);
      }
#ifdef CYLINDRICAL
      dtodx2 = my_dt/(ri[ie+1]*pG->dx2);
#endif
      pG->B1i[ks][j][ie+1] -= dtodx2*(emf[ks][j+1][ie+1].x3 -emf[ks][j][ie+1].x3);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[ks][je+1][i] += dtodx1*(emf[ks][je+1][i+1].x3 -emf[ks][je+1][i].x3);
    } 
/* Set cell centered magnetic fields to average of face centered */
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        pG->U[ks][j][i].B1c = 0.5*(lsf*pG->B1i[ks][j][i] + rsf*pG->B1i[ks][j][i+1]);
        pG->U[ks][j][i].B2c = 0.5*(    pG->B2i[ks][j][i] +     pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
        pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
      }
    }
  }

/* 3D PROBLEM: CT */
  if (ndim == 3){
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
          dtodx2 = my_dt/(ri[i]*pG->dx2);
#endif
          pG->B1i[k][j][i] += dtodx3*(emf[k+1][j  ][i  ].x2 - emf[k][j][i].x2) -
                              dtodx2*(emf[k  ][j+1][i  ].x3 - emf[k][j][i].x3);
          pG->B2i[k][j][i] += dtodx1*(emf[k  ][j  ][i+1].x3 - emf[k][j][i].x3) -
                              dtodx3*(emf[k+1][j  ][i  ].x1 - emf[k][j][i].x1);
#ifdef CYLINDRICAL
          rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
          dtodx2 = my_dt/(r[i]*pG->dx2);
#endif
          pG->B3i[k][j][i] += dtodx2*(    emf[k  ][j+1][i  ].x1 -     emf[k][j][i].x1) -
                              dtodx1*(rsf*emf[k  ][j  ][i+1].x2 - lsf*emf[k][j][i].x2);
        }
#ifdef CYLINDRICAL
        dtodx2 = my_dt/(ri[ie+1]*pG->dx2);
#endif  
        pG->B1i[k][j][ie+1] +=
          dtodx3*(emf[k+1][j  ][ie+1].x2 - emf[k][j][ie+1].x2) -
          dtodx2*(emf[k  ][j+1][ie+1].x3 - emf[k][j][ie+1].x3);
      }
      for (i=is; i<=ie; i++) {
        pG->B2i[k][je+1][i] +=
          dtodx1*(emf[k  ][je+1][i+1].x3 - emf[k][je+1][i].x3) -
          dtodx3*(emf[k+1][je+1][i  ].x1 - emf[k][je+1][i].x1);
      }
    }
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      dtodx2 = my_dt/(r[i]*pG->dx2);
#endif
      pG->B3i[ke+1][j][i] +=
        dtodx2*(    emf[ke+1][j+1][i  ].x1 -     emf[ke+1][j][i].x1) -
        dtodx1*(rsf*emf[ke+1][j  ][i+1].x2 - lsf*emf[ke+1][j][i].x2);
    }}
/* Set cell centered magnetic fields to average of face centered */
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
      pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
      pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
    }}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* EField_Ohm:  Resistive EMF from Ohmic dissipation.   E = \eta_Ohm J
 */

void EField_Ohm(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
  Real eta_O;

#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif
  Real lsf=1.0,rsf=1.0;

  if (pG->Nx[1] > 1)    ndim++;
  if (pG->Nx[2] > 1)    ndim++;

/* For Ohmic resistivity, E = \eta_Ohm J  */

/* 1D PROBLEM: */
  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {

#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      eta_O = 0.5*(rsf*pG->eta_Ohm[ks][js][i] + lsf*pG->eta_Ohm[ks][js][i-1]);

      emf[ks][js][i].x2 += eta_O * J[ks][js][i].x2;
      emf[ks][js][i].x3 += eta_O * J[ks][js][i].x3;
    }
  }

/* 2D PROBLEM: */
  if (ndim == 2){
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      eta_O = 0.5*(pG->eta_Ohm[ks][j][i] + pG->eta_Ohm[ks][j-1][i]);

      emf[ks][j][i].x1 += eta_O * J[ks][j][i].x1;

#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      eta_O = 0.5*(rsf*pG->eta_Ohm[ks][j][i] + lsf*pG->eta_Ohm[ks][j][i-1]);

      emf[ks][j][i].x2 += eta_O * J[ks][j][i].x2; 

      eta_O = 0.25*(rsf*pG->eta_Ohm[ks][j][i  ] + rsf*pG->eta_Ohm[ks][j-1][i  ] +
                    lsf*pG->eta_Ohm[ks][j][i-1] + lsf*pG->eta_Ohm[ks][j-1][i-1]);

      emf[ks][j][i].x3 += eta_O * J[ks][j][i].x3;
    }}
  }


/* 3D PROBLEM: */
  
  if (ndim == 3){

    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {

        eta_O = 0.25*(pG->eta_Ohm[k][j  ][i] + pG->eta_Ohm[k-1][j  ][i] +
                      pG->eta_Ohm[k][j-1][i] + pG->eta_Ohm[k-1][j-1][i]);

        emf[k][j][i].x1 += eta_O * J[k][j][i].x1;

#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        eta_O = 0.25*(rsf*pG->eta_Ohm[k][j][i  ] + rsf*pG->eta_Ohm[k-1][j][i  ] +
                      lsf*pG->eta_Ohm[k][j][i-1] + lsf*pG->eta_Ohm[k-1][j][i-1]);

        emf[k][j][i].x2 += eta_O * J[k][j][i].x2;

        eta_O = 0.25*(rsf*pG->eta_Ohm[k][j][i  ] + rsf*pG->eta_Ohm[k][j-1][i  ] +
                      lsf*pG->eta_Ohm[k][j][i-1] + lsf*pG->eta_Ohm[k][j-1][i-1]);

        emf[k][j][i].x3 += eta_O * J[k][j][i].x3;
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* EField_Hall:  Resistive EMF from Hall effect.  E = Q_H (J X B)
 */

void EField_Hall(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, il,iu, is = pG->is, ie = pG->ie;
  int j, jl,ju, js = pG->js, je = pG->je;
  int k, kl,ku, ks = pG->ks, ke = pG->ke;
#ifdef SHEARING_BOX
  int my_iproc, my_jproc, my_kproc;
  int nlayer;
#endif
  int ndim=1;
  Real eta_H, Bmag;
#ifdef STS
  Real my_dt = STS_dt;
#else
  Real my_dt = pG->dt;
#endif
  Real dtodx1 = my_dt/pG->dx1, dtodx2 = 0.0, dtodx3 = 0.0;

  il = is - 4;  iu = ie + 4;

  if (pG->Nx[1] > 1){
    jl = js - 4;    ju = je + 4;
    dtodx2 = my_dt/pG->dx2;
    ndim++;
  } else {
    jl = js;        ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 4;    ku = ke + 4;
    dtodx3 = my_dt/pG->dx3;
    ndim++;
  } else {
    kl = ks;        ku = ke;
  }

/* Preliminary: hyper-diffusion */

//  hyper_diffusion4(pD, 0.1);

/* Preliminary: divide eta_Hall by B for convenience */
  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=il; i<=iu; i++) {

      Bmag = sqrt(SQR(pG->U[k][j][i].B1c)
                + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));

      pG->eta_Hall[k][j][i] /= Bmag+TINY_NUMBER;
  }}}

#ifdef SHEARING_BOX
  get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
#endif

/* 1D PROBLEM:
 *   emf.x =  0.0
 *   emf.y =  Jz*Bx
 *   emf.z = -Jy*Bx   */

  if (ndim == 1){
    /* x2-sweep */
    for (i=is-2; i<=ie+3; i++) {
      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);

      emfh[ks][js][i].x2 = eta_H*J[ks][js][i].x3 * pG->B1i[ks][js][i];
    }

    /* update magnetic field */
    for (i=is-2; i<=ie+2; i++) {
      Bcor[ks][js][i].x3 -= dtodx1*(emfh[ks][js][i+1].x2 - emfh[ks][js][i].x2);
    }

    /* update current */
    for (i=is-1; i<=ie+2; i++) {
      Jcor[ks][js][i].x2 = -(Bcor[ks][js][i].x3 - Bcor[ks][js][i-1].x3)/pG->dx1;
    }

    /* x3-sweep */
    for (i=is-1; i<=ie+2; i++) {
      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);

      emfh[ks][js][i].x3 = -eta_H*Jcor[ks][js][i].x2 * pG->B1i[ks][js][i];
    }

    /* update global emf */
    for (i=is-1; i<=ie+1; i++) {
      emf[ks][js][i].x2 += emfh[ks][js][i].x2;
      emf[ks][js][i].x3 += emfh[ks][js][i].x3;
    }

  }

/* 2D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 2){
    for (j=js-3; j<=je+4; j++) {
    for (i=is-3; i<=ie+3; i++) {

      /* x1-sweep */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j-1][i]);

      emfh[ks][j][i].x1  = eta_H*(
        0.125*(   J[ks][j  ][i].x2 +    J[ks][j  ][i+1].x2
             +    J[ks][j-1][i].x2 +    J[ks][j-1][i+1].x2)
             *(Bcor[ks][j  ][i].x3 + Bcor[ks][j-1][i  ].x3) -
         0.5*((   J[ks][j  ][i].x3 +    J[ks][j  ][i+1].x3)*Bcor[ks][j][i].x2));
      }}

    /* update the magnetic field */
    for (j=js-3; j<=je+3; j++) {
    for (i=is-3; i<=ie+3; i++) {
        Bcor[ks][j][i].x3 += dtodx2*(emfh[ks][j+1][i].x1 - emfh[ks][j][i].x1);
    }}

    /* update the current density */
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+3; i++) {
        Jcor[ks][j][i].x1=  (Bcor[ks][j][i].x3 - Bcor[ks][j-1][i  ].x3)/pG->dx2;
        Jcor[ks][j][i].x2= -(Bcor[ks][j][i].x3 - Bcor[ks][j  ][i-1].x3)/pG->dx1;
        Jcor[ks][j][i].x3=  J[ks][j][i].x3;
      }
      i = is-3;
      Jcor[ks][j][i].x1=  (Bcor[ks][j][i].x3 - Bcor[ks][j-1][i  ].x3)/pG->dx2;
    }
    j = js-3;
    for (i=is-2; i<=ie+3; i++) {
      Jcor[ks][j][i].x2= -(Bcor[ks][j][i].x3 - Bcor[ks][j  ][i-1].x3)/pG->dx1;
    }

    /* x2-sweep */
    for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+3; i++) {
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j][i-1]);

      emfh[ks][j][i].x2 = eta_H*(
         0.5*((Jcor[ks][j][i  ].x3 + Jcor[ks][j+1][i  ].x3)*Bcor[ks][j][i].x1) -
        0.125*(Jcor[ks][j][i  ].x1 + Jcor[ks][j+1][i  ].x1
             + Jcor[ks][j][i-1].x1 + Jcor[ks][j+1][i-1].x1)
             *(Bcor[ks][j][i  ].x3 + Bcor[ks][j  ][i-1].x3));
    }}

    /* update the magnetic field */
    for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      Bcor[ks][j][i].x3 -= dtodx1*(emfh[ks][j][i+1].x2 - emfh[ks][j][i].x2);
    }}

    /* update the current density */
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
        Jcor[ks][j][i].x1=  (Bcor[ks][j][i].x3 - Bcor[ks][j-1][i  ].x3)/pG->dx2;
        Jcor[ks][j][i].x2= -(Bcor[ks][j][i].x3 - Bcor[ks][j  ][i-1].x3)/pG->dx1;
      }
      i = is-2;
      Jcor[ks][j][i].x1=  (Bcor[ks][j][i].x3 - Bcor[ks][j-1][i  ].x3)/pG->dx2;
    }
    j = js-2;
    for (i=is-1; i<=ie+2; i++) {
      Jcor[ks][j][i].x2= -(Bcor[ks][j][i].x3 - Bcor[ks][j  ][i-1].x3)/pG->dx1;
    }

    /* x3-sweep */
    for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
      eta_H = 0.25*(pG->eta_Hall[ks][j][i  ] + pG->eta_Hall[ks][j-1][i  ] +
                    pG->eta_Hall[ks][j][i-1] + pG->eta_Hall[ks][j-1][i-1]);

      emfh[ks][j][i].x3 = eta_H*(
        0.25*(Jcor[ks][j][i].x1 + Jcor[ks][j][i-1].x1)
            *(Bcor[ks][j][i].x2 + Bcor[ks][j][i-1].x2) -
        0.25*(Jcor[ks][j][i].x2 + Jcor[ks][j-1][i].x2)
            *(Bcor[ks][j][i].x1 + Bcor[ks][j-1][i].x1));
    }}

    /* update the global emf */
    for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) { 
      emf[ks][j][i].x1 += emfh[ks][j][i].x1;
      emf[ks][j][i].x2 += emfh[ks][j][i].x2;
      emf[ks][j][i].x3 += emfh[ks][j][i].x3;
    }}

  }

/* 3D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 3){

    /* x1-sweep */
    for (k=ks-3; k<=ke+4; k++) {
    for (j=js-3; j<=je+4; j++) {
      for (i=is-3; i<=ie+3; i++) {

        eta_H = 0.25*(pG->eta_Hall[k][j  ][i] + pG->eta_Hall[k-1][j  ][i] +
                      pG->eta_Hall[k][j-1][i] + pG->eta_Hall[k-1][j-1][i]);

        emfh[k][j][i].x1 = 0.125*eta_H*(
                (J[k  ][j  ][i].x2    + J[k  ][j  ][i+1].x2
               + J[k  ][j-1][i].x2    + J[k  ][j-1][i+1].x2)
               *(Bcor[k  ][j  ][i].x3 + Bcor[k  ][j-1][i  ].x3)-
                (J[k  ][j  ][i].x3    + J[k  ][j  ][i+1].x3
               + J[k-1][j  ][i].x3    + J[k-1][j  ][i+1].x3)
               *(Bcor[k  ][j  ][i].x2 + Bcor[k-1][j  ][i  ].x2));
      }
    }}

    /* update the magnetic field */
    for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        Bcor[k][j][i].x2 -= dtodx3*(emfh[k+1][j  ][i  ].x1 - emfh[k][j][i].x1);
        Bcor[k][j][i].x3 += dtodx2*(emfh[k  ][j+1][i  ].x1 - emfh[k][j][i].x1);
      }
    }}

    /* update the current density */
    for (k=ks-2; k<=ke+3; k++) {
      for (j=js-2; j<=je+3; j++) {
        for (i=is-2; i<=ie+3; i++) {
          Jcor[k][j][i].x1 = (Bcor[k][j][i].x3 - Bcor[k  ][j-1][i  ].x3)/pG->dx2 -
                             (Bcor[k][j][i].x2 - Bcor[k-1][j  ][i  ].x2)/pG->dx3;
          Jcor[k][j][i].x2 = (Bcor[k][j][i].x1 - Bcor[k-1][j  ][i  ].x1)/pG->dx3 -
                             (Bcor[k][j][i].x3 - Bcor[k  ][j  ][i-1].x3)/pG->dx1;
          Jcor[k][j][i].x3 = (Bcor[k][j][i].x2 - Bcor[k  ][j  ][i-1].x2)/pG->dx1 -
                             (Bcor[k][j][i].x1 - Bcor[k  ][j-1][i  ].x1)/pG->dx2;
        }
        i = is-3;
        Jcor[k][j][i].x1 = (Bcor[k][j][i].x3 - Bcor[k  ][j-1][i  ].x3)/pG->dx2 -
                           (Bcor[k][j][i].x2 - Bcor[k-1][j  ][i  ].x2)/pG->dx3;
      }
      j = js-3;
      for (i=is-2; i<=ie+3; i++) {
        Jcor[k][j][i].x2 = (Bcor[k][j][i].x1 - Bcor[k-1][j  ][i  ].x1)/pG->dx3 -
                           (Bcor[k][j][i].x3 - Bcor[k  ][j  ][i-1].x3)/pG->dx1;
      }
    }
    k = ks-3;
    for (j=js-2; j<=je+3; j++) {
    for (i=is-2; i<=ie+3; i++) {
      Jcor[k][j][i].x3 = (Bcor[k][j][i].x2 - Bcor[k  ][j  ][i-1].x2)/pG->dx1 -
                         (Bcor[k][j][i].x1 - Bcor[k  ][j-1][i  ].x1)/pG->dx2;
    }}

#ifdef SHEARING_BOX
    nlayer = 4;

/* remap Jy */
    for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
    for(j=js; j<=je; j++) {
    for(i=is; i<=ie+1; i++) {
      J2[k][j][i]   = Jcor[k][j][i].x2;
    }}}

    if (my_iproc == 0) {
      RemapJy_ix1(pD, J2, remapJyiib, nlayer);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapJy_ox1(pD, J2, remapJyoib, nlayer);
    }

/* Now average Ey and remapped Ey */
    if (my_iproc == 0) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=is-nlayer+1; i<=is; i++) {
        Jcor[k][j][i].x2  = 0.5*(Jcor[k][j][i].x2 + remapJyiib[i-is+nlayer-1][k][j]);
      }}}
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=ie+1; i<=ie+nlayer; i++) {
        Jcor[k][j][i].x2 = 0.5*(Jcor[k][j][i].x2 + remapJyoib[i-ie-1][k][j]);
      }}}
    }
#endif /* SHEARING_BOX */

    /* x2-sweep */
    for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k-1][j][i  ] +
                      pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k-1][j][i-1]);

        emfh[k][j][i].x2 += 0.125*eta_H*(
                (Jcor[k  ][j][i  ].x3 + Jcor[k  ][j+1][i  ].x3
               + Jcor[k-1][j][i  ].x3 + Jcor[k-1][j+1][i  ].x3)
               *(Bcor[k  ][j][i  ].x1 + Bcor[k-1][j  ][i  ].x1)-
                (Jcor[k  ][j][i  ].x1 + Jcor[k  ][j+1][i  ].x1
               + Jcor[k  ][j][i-1].x1 + Jcor[k  ][j+1][i-1].x1)
               *(Bcor[k  ][j][i  ].x3 + Bcor[k  ][j  ][i-1].x3));
      }
    }}

    /* update the magnetic field */
    for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Bcor[k][j][i].x1 += dtodx3*(emfh[k+1][j  ][i  ].x2 - emfh[k][j][i].x2);
        Bcor[k][j][i].x3 -= dtodx1*(emfh[k  ][j  ][i+1].x2 - emfh[k][j][i].x2);
      }
    }}

    /* update the current density */
    for (k=ks-1; k<=ke+2; k++) {
      for (j=js-1; j<=je+2; j++) {
        for (i=is-1; i<=ie+2; i++) {
          Jcor[k][j][i].x1 = (Bcor[k][j][i].x3 - Bcor[k  ][j-1][i  ].x3)/pG->dx2 -
                             (Bcor[k][j][i].x2 - Bcor[k-1][j  ][i  ].x2)/pG->dx3;
          Jcor[k][j][i].x2 = (Bcor[k][j][i].x1 - Bcor[k-1][j  ][i  ].x1)/pG->dx3 -
                             (Bcor[k][j][i].x3 - Bcor[k  ][j  ][i-1].x3)/pG->dx1;
          Jcor[k][j][i].x3 = (Bcor[k][j][i].x2 - Bcor[k  ][j  ][i-1].x2)/pG->dx1 -
                             (Bcor[k][j][i].x1 - Bcor[k  ][j-1][i  ].x1)/pG->dx2;
        }
        i = is-2;
        Jcor[k][j][i].x1 = (Bcor[k][j][i].x3 - Bcor[k  ][j-1][i  ].x3)/pG->dx2 -
                           (Bcor[k][j][i].x2 - Bcor[k-1][j  ][i  ].x2)/pG->dx3;
      }
      j = js-2;
      for (i=is-1; i<=ie+2; i++) {
        Jcor[k][j][i].x2 = (Bcor[k][j][i].x1 - Bcor[k-1][j  ][i  ].x1)/pG->dx3 -
                           (Bcor[k][j][i].x3 - Bcor[k  ][j  ][i-1].x3)/pG->dx1;
      }
    }
    k = ks-2;
    for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
      Jcor[k][j][i].x3 = (Bcor[k][j][i].x2 - Bcor[k  ][j  ][i-1].x2)/pG->dx1 -
                         (Bcor[k][j][i].x1 - Bcor[k  ][j-1][i  ].x1)/pG->dx2;
    }}

#ifdef SHEARING_BOX
    nlayer = 3;

/* remap Jy */
    for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
    for(j=js; j<=je; j++) {
    for(i=is; i<=ie+1; i++) {
      J2[k][j][i]   = Jcor[k][j][i].x2;
    }}}

    if (my_iproc == 0) {
      RemapJy_ix1(pD, J2, remapJyiib, nlayer);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapJy_ox1(pD, J2, remapJyoib, nlayer);
    }

/* Now average Ey and remapped Ey */
    if (my_iproc == 0) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=is-nlayer+1; i<=is; i++) {
        Jcor[k][j][i].x2  = 0.5*(Jcor[k][j][i].x2 + remapJyiib[i-is+nlayer-1][k][j]);
      }}}
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks-nlayer+1; k<=ke+nlayer; k++) {
      for(j=js-nlayer; j<=je+nlayer; j++) {
      for(i=ie+1; i<=ie+nlayer; i++) {
        Jcor[k][j][i].x2 = 0.5*(Jcor[k][j][i].x2 + remapJyoib[i-ie-1][k][j]);
      }}}
    }
#endif /* SHEARING_BOX */

    /* x3-sweep */
    for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
        eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k][j-1][i  ] +
                      pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k][j-1][i-1]);

        emfh[k][j][i].x3 = 0.125*eta_H*(
                (Jcor[k][j  ][i  ].x1 + Jcor[k+1][j  ][i  ].x1
               + Jcor[k][j  ][i-1].x1 + Jcor[k+1][j  ][i-1].x1)
               *(Bcor[k][j  ][i  ].x2 + Bcor[k  ][j  ][i-1].x2)-
                (Jcor[k][j  ][i  ].x2 + Jcor[k+1][j  ][i  ].x2
               + Jcor[k][j-1][i  ].x2 + Jcor[k+1][j-1][i  ].x2)
               *(Bcor[k][j  ][i  ].x1 + Bcor[k  ][j-1][i  ].x1));
      }
    }}

    /* update the global emf */
    for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        emf[k][j][i].x1 += emfh[k][j][i].x1;
        emf[k][j][i].x2 += emfh[k][j][i].x2;
        emf[k][j][i].x3 += emfh[k][j][i].x3;
      }
    }}

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* EField_AD:  Resistive EMF from ambipolar diffusion.  E = Q_AD (J X B) X B
 */

void EField_AD(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
  Real eta_A;
  Real intBx,intBy,intBz,intJx,intJy,intJz,Bsq,JdotB;

#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif
  Real lsf=1.0,rsf=1.0;

  if (pG->Nx[1] > 1) ndim++;
  if (pG->Nx[2] > 1) ndim++;

/* 1D PROBLEM:
 *   emf.x = 0
 *   emf.y = (J_perp)_y
 *   emf.z = (J_perp)_z  */

  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      eta_A = 0.5*(rsf*pG->eta_AD[ks][js][i] + lsf*pG->eta_AD[ks][js][i-1]);

      intBx = pG->B1i[ks][js][i];
      intBy = 0.5*(rsf*pG->U[ks][js][i].B2c + lsf*pG->U[ks][js][i-1].B2c);
      intBz = 0.5*(rsf*pG->U[ks][js][i].B3c + lsf*pG->U[ks][js][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
      JdotB = J[ks][js][i].x2*intBy + J[ks][js][i].x3*intBz;

      emf[ks][js][i].x2 += eta_A*(J[ks][js][i].x2 - JdotB*intBy/Bsq);
      emf[ks][js][i].x3 += eta_A*(J[ks][js][i].x3 - JdotB*intBz/Bsq);
    }
  }

/* 2D PROBLEM:
 *   emf.x = (J_perp)_x
 *   emf.y = (J_perp)_y
 *   emf.z = (J_perp)_z  */

  if (ndim == 2){
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* emf.x */
#ifdef CYLINDRICAL
      lsf = ri[i]/r[i];  rsf = ri[i+1]/r[i];
#endif
      eta_A = 0.5*(pG->eta_AD[ks][j][i] + pG->eta_AD[ks][j-1][i]);

      intJx = J[ks][j][i].x1;
      intJy = 0.25*(lsf*J[ks][j  ][i].x2 + rsf*J[ks][j  ][i+1].x2
                  + lsf*J[ks][j-1][i].x2 + rsf*J[ks][j-1][i+1].x2);
      intJz = 0.5 *(lsf*J[ks][j][i].x3   + rsf*J[ks][j][i+1].x3);

      intBx = 0.5*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c);
      intBy = pG->B2i[ks][j][i];
      intBz = 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].x1 += eta_A*(J[ks][j][i].x1 - JdotB*intBx/Bsq);

      /* emf.y */
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      eta_A = 0.5*(rsf*pG->eta_AD[ks][j][i] + lsf*pG->eta_AD[ks][j][i-1]);

      intJx = 0.25*(rsf*J[ks][j][i  ].x1 + rsf*J[ks][j+1][i  ].x1
                  + lsf*J[ks][j][i-1].x1 + lsf*J[ks][j+1][i-1].x1);
      intJy = J[ks][j][i].x2;
      intJz = 0.5 *(    J[ks][j][i].x3   +     J[ks][j+1][i].x3);

      intBx = pG->B1i[ks][j][i];
      intBy = 0.5*(rsf*pG->U[ks][j][i].B2c + lsf*pG->U[ks][j][i-1].B2c);
      intBz = 0.5*(lsf*pG->U[ks][j][i].B3c + lsf*pG->U[ks][j][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].x2 += eta_A*(J[ks][j][i].x2 - JdotB*intBy/Bsq);

     /* emf.z */
#ifdef CYLINDRICAL
      rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
      eta_A = 0.25*(rsf*pG->eta_AD[ks][j  ][i] + lsf*pG->eta_AD[ks][j  ][i-1]
                  + rsf*pG->eta_AD[ks][j-1][i] + lsf*pG->eta_AD[ks][j-1][i-1]);

      intJx = 0.5*(rsf*J[ks][j][i].x1 + lsf*J[ks][j][i-1].x1);
      intJy = 0.5*(    J[ks][j][i].x2 +     J[ks][j-1][i].x2);
      intJz = J[ks][j][i].x3;

      intBx = 0.5*(    pG->B1i[ks][j][i] +     pG->B1i[ks][j-1][i]);
      intBy = 0.5*(rsf*pG->B2i[ks][j][i] + lsf*pG->B2i[ks][j][i-1]);
      intBz = 0.25*(rsf*pG->U[ks][j  ][i].B3c + rsf*pG->U[ks][j  ][i-1].B3c
                  + lsf*pG->U[ks][j-1][i].B3c + lsf*pG->U[ks][j-1][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].x3 += eta_A*(J[ks][j][i].x3 - JdotB*intBz/Bsq);
    }}
  }

/* 3D PROBLEM:
 *   emf.x = (J_perp)_x
 *   emf.y = (J_perp)_y
 *   emf.z = (J_perp)_z  */

  if (ndim == 3){
    for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {

        /* emf.x */
#ifdef CYLINDRICAL
        lsf = ri[i]/r[i];  rsf = ri[i+1]/r[i];
#endif
        eta_A = 0.25*(pG->eta_AD[k][j  ][i] + pG->eta_AD[k-1][j  ][i] +
                      pG->eta_AD[k][j-1][i] + pG->eta_AD[k-1][j-1][i]);

        intJx = J[k][j][i].x1;
        intJy = 0.25*(lsf*J[k][j  ][i].x2 + rsf*J[k][j  ][i+1].x2
                    + lsf*J[k][j-1][i].x2 + rsf*J[k][j-1][i+1].x2);
        intJz = 0.25*(lsf*J[k  ][j][i].x3 + rsf*J[k  ][j][i+1].x3
                    + lsf*J[k-1][j][i].x3 + rsf*J[k-1][j][i+1].x3);

        intBx = 0.25*(pG->U[k][j  ][i].B1c + pG->U[k-1][j  ][i].B1c +
                      pG->U[k][j-1][i].B1c + pG->U[k-1][j-1][i].B1c);
        intBy = 0.5*(pG->B2i[k][j][i] + pG->B2i[k-1][j][i]);
        intBz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k][j-1][i]);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].x1 += eta_A*(J[k][j][i].x1 - JdotB*intBx/Bsq);

        /* emf.y */
#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        eta_A = 0.25*(rsf*pG->eta_AD[k][j][i  ] + rsf*pG->eta_AD[k-1][j][i  ] +
                      lsf*pG->eta_AD[k][j][i-1] + lsf*pG->eta_AD[k-1][j][i-1]);

        intJx = 0.25*(rsf*J[k][j][i  ].x1 + rsf*J[k][j+1][i  ].x1
                    + lsf*J[k][j][i-1].x1 + lsf*J[k][j+1][i-1].x1);;
        intJy = J[k][j][i].x2;
        intJz = 0.25*(    J[k  ][j][i].x3 +     J[k  ][j+1][i].x3
                    +     J[k-1][j][i].x3 +     J[k-1][j+1][i].x3);

        intBx = 0.5*(    pG->B1i[k][j][i] +     pG->B1i[k-1][j][i]);
        intBy = 0.25*(rsf*pG->U[k][j][i  ].B2c + rsf*pG->U[k-1][j][i  ].B2c +
                      lsf*pG->U[k][j][i-1].B2c + lsf*pG->U[k-1][j][i-1].B2c);
        intBz = 0.5*(rsf*pG->B3i[k][j][i] + lsf*pG->B3i[k][j][i-1]);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].x2 += eta_A*(J[k][j][i].x2 - JdotB*intBy/Bsq);

        /* emf.z */
#ifdef CYLINDRICAL
        rsf = r[i]/ri[i];  lsf = r[i-1]/ri[i];
#endif
        eta_A = 0.25*(rsf*pG->eta_AD[k][j][i  ] + rsf*pG->eta_AD[k][j-1][i  ] +
                      lsf*pG->eta_AD[k][j][i-1] + lsf*pG->eta_AD[k][j-1][i-1]);

        intJx = 0.25*(rsf*J[k][j][i  ].x1 + rsf*J[k+1][j][i  ].x1
                    + lsf*J[k][j][i-1].x1 + lsf*J[k+1][j][i-1].x1);;
        intJy = 0.25*(    J[k][j  ][i].x2 +     J[k+1][j  ][i].x2
                    +     J[k][j-1][i].x2 +     J[k+1][j-1][i].x2);
        intJz = J[k][j][i].x3;

        intBx = 0.5*(    pG->B1i[k][j][i] +     pG->B1i[k][j-1][i]);
        intBy = 0.5*(rsf*pG->B2i[k][j][i] + lsf*pG->B2i[k][j][i-1]);
        intBz = 0.25*(rsf*pG->U[k][j  ][i].B3c + rsf*pG->U[k][j  ][i-1].B3c +
                      lsf*pG->U[k][j-1][i].B3c + lsf*pG->U[k][j-1][i-1].B3c);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz) + TINY_NUMBER;
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].x3 += eta_A*(J[k][j][i].x3 - JdotB*intBz/Bsq);
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* hyper_diffusion: calculate the higher-order derivatives of J
 */  

/* 4th order diffusion */
void hyper_diffusion4(DomainS *pD, Real prefac)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
#ifdef STS
  Real my_dt = STS_dt;
#else
  Real my_dt = pG->dt;
#endif
  Real eta_H,eta_4,dx21,dy21=0.0,dz21=0.0;

  dx21 = 1.0/SQR(pG->dx1);
  if (pG->Nx[1]>1) {
    ndim++;
    dy21 = 1.0/SQR(pG->dx2);
  }
  if (pG->Nx[2]>1) {
    ndim++;
    dz21 = 1.0/SQR(pG->dx3);
  }

  /* 1D */
  if (ndim == 1) {
    for (i=is; i<=ie+1; i++) {
      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[ks][js][i].x2 -= eta_4 * (J[ks][js][i-1].x2 - 2.0*J[ks][js][i].x2
                                 + J[ks][js][i+1].x2) * dx21;
      emf[ks][js][i].x3 -= eta_4 * (J[ks][js][i-1].x3 - 2.0*J[ks][js][i].x3
                                 + J[ks][js][i+1].x3) * dx21;
    }
  }

  /* 2D */
  if (ndim == 2) {
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* x1 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j-1][i]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[ks][j][i].x1 -= eta_4 * ((J[ks][j][i-1].x1 - 2.0*J[ks][j][i].x1
                                 + J[ks][j][i+1].x1) * dx21
                                 +(J[ks][j-1][i].x1 - 2.0*J[ks][j][i].x1
                                 + J[ks][j+1][i].x1) * dy21);

      /* x2 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j][i-1]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[ks][j][i].x2 -= eta_4 * ((J[ks][j][i-1].x2 - 2.0*J[ks][j][i].x2
                                 + J[ks][j][i+1].x2) * dx21
                                 +(J[ks][j-1][i].x2 - 2.0*J[ks][j][i].x2
                                 + J[ks][j+1][i].x2) * dy21);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[ks][j][i  ] + pG->eta_Hall[ks][j-1][i  ] +
                    pG->eta_Hall[ks][j][i-1] + pG->eta_Hall[ks][j-1][i-1]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[ks][j][i].x3 -= eta_4 * ((J[ks][j][i-1].x3 - 2.0*J[ks][j][i].x3
                                 + J[ks][j][i+1].x3) * dx21
                                 +(J[ks][j-1][i].x3 - 2.0*J[ks][j][i].x3
                                 + J[ks][j+1][i].x3) * dy21);
    }}
  }

  /* 3D */
  if (ndim == 3) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* x1 */
      eta_H = 0.25*(pG->eta_Hall[k][j  ][i] + pG->eta_Hall[k-1][j  ][i] +
                    pG->eta_Hall[k][j-1][i] + pG->eta_Hall[k-1][j-1][i]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[k][j][i].x1 -= eta_4 * ((J[k][j][i-1].x1 - 2.0*J[k][j][i].x1
                                + J[k][j][i+1].x1) * dx21
                                +(J[k][j-1][i].x1 - 2.0*J[k][j][i].x1
                                + J[k][j+1][i].x1) * dy21
                                +(J[k-1][j][i].x1 - 2.0*J[k][j][i].x1
                                + J[k+1][j][i].x1) * dz21);

      /* x2 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k-1][j][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k-1][j][i-1]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[k][j][i].x2 -= eta_4 * ((J[k][j][i-1].x2 - 2.0*J[k][j][i].x2
                                + J[k][j][i+1].x2) * dx21
                                +(J[k][j-1][i].x2 - 2.0*J[k][j][i].x2
                                + J[k][j+1][i].x2) * dy21
                                +(J[k-1][j][i].x2 - 2.0*J[k][j][i].x2
                                + J[k+1][j][i].x2) * dz21);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k][j-1][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k][j-1][i-1]);
      eta_4 = prefac * SQR(eta_H) * my_dt;

      emf[k][j][i].x3 -= eta_4 * ((J[k][j][i-1].x3 - 2.0*J[k][j][i].x3
                                + J[k][j][i+1].x3) * dx21
                                +(J[k][j-1][i].x3 - 2.0*J[k][j][i].x3
                                + J[k][j+1][i].x3) * dy21
                                +(J[k-1][j][i].x3 - 2.0*J[k][j][i].x3
                                + J[k+1][j][i].x3) * dz21);
    }}}
  }

  return;
}

/* 6th order diffusion */
void hyper_diffusion6(DomainS *pD, Real prefac)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
#ifdef STS
  Real my_dt = STS_dt;
#else
  Real my_dt = pG->dt;
#endif
  Real eta_H,eta_6,dx41,dy41=0.0,dz41=0.0;
  Real fac,fac2,fac3;

  dx41 = 1.0/SQR(SQR(pG->dx1)); 
  if (pG->Nx[1]>1) {
    ndim++;
    dy41 = 1.0/SQR(SQR(pG->dx2));
    fac2 = SQR(pG->dx1/pG->dx2);
  }   
  if (pG->Nx[2]>1) {
    ndim++;
    dz41 = 1.0/SQR(SQR(pG->dx3)); 
    fac3 = SQR(pG->dx1/pG->dx3);
  }
  fac = prefac*SQR(my_dt/pG->dx1)*my_dt;

  /* 1D */
  if (ndim == 1) {
    for (i=is; i<=ie+1; i++) {

      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][js][i].x2 += eta_6 *  (J[ks][js][i-2].x2 - 4.0*J[ks][js][i-1].x2
         + 6.0*J[ks][js][i].x2 - 4.0*J[ks][js][i+1].x2 + J[ks][js][i+2].x2) * dx41;
      emf[ks][js][i].x3 += eta_6 *  (J[ks][js][i-2].x3 - 4.0*J[ks][js][i-1].x3
         + 6.0*J[ks][js][i].x3 - 4.0*J[ks][js][i+1].x3 + J[ks][js][i+2].x3) * dx41;
    }
  }

  /* 2D */
  if (ndim == 2) {
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* x1 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j-1][i]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].x1 += eta_6 * ((J[ks][j][i-2].x1 - 4.0*J[ks][j][i-1].x1
         + 6.0*J[ks][j][i].x1 - 4.0*J[ks][j][i+1].x1 + J[ks][j][i+2].x1)*dx41
                         + fac2 * (J[ks][j-2][i].x1 - 4.0*J[ks][j-1][i].x1
         + 6.0*J[ks][j][i].x1 - 4.0*J[ks][j+1][i].x1 + J[ks][j+2][i].x1) * dy41);

      /* x2 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].x2 += eta_6 * ((J[ks][j][i-2].x2 - 4.0*J[ks][j][i-1].x2
         + 6.0*J[ks][j][i].x2 - 4.0*J[ks][j][i+1].x2 + J[ks][j][i+2].x2) * dx41
                         + fac2 * (J[ks][j-2][i].x2 - 4.0*J[ks][j-1][i].x2
         + 6.0*J[ks][j][i].x2 - 4.0*J[ks][j+1][i].x2 + J[ks][j+2][i].x2) * dy41);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[ks][j][i  ] + pG->eta_Hall[ks][j-1][i  ] +
                    pG->eta_Hall[ks][j][i-1] + pG->eta_Hall[ks][j-1][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].x3 += eta_6 * ((J[ks][j][i-2].x3 - 4.0*J[ks][j][i-1].x3
         + 6.0*J[ks][j][i].x3 - 4.0*J[ks][j][i+1].x3 + J[ks][j][i+2].x3) * dx41
                         + fac2 * (J[ks][j-2][i].x3 - 4.0*J[ks][j-1][i].x3
         + 6.0*J[ks][j][i].x3 - 4.0*J[ks][j+1][i].x3 + J[ks][j+2][i].x3) * dy41);
    }}
  }

  /* 3D */
  if (ndim == 3) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* x1 */
      eta_H = 0.25*(pG->eta_Hall[k][j  ][i] + pG->eta_Hall[k-1][j  ][i] +
                    pG->eta_Hall[k][j-1][i] + pG->eta_Hall[k-1][j-1][i]);
      eta_6 = SQR(SQR(eta_H)) * fac;
      emf[k][j][i].x1 += eta_6 * ((J[k][j][i-2].x1 - 4.0*J[k][j][i-1].x1
         + 6.0*J[k][j][i].x1 - 4.0*J[k][j][i+1].x1 + J[k][j][i+2].x1) * dx41
                        + fac2 * (J[k][j-2][i].x1 - 4.0*J[k][j-1][i].x1
         + 6.0*J[k][j][i].x1 - 4.0*J[k][j+1][i].x1 + J[k][j+2][i].x1) * dy41
                        + fac3 * (J[k-2][j][i].x1 - 4.0*J[k-1][j][i].x1
         + 6.0*J[k][j][i].x1 - 4.0*J[k+1][j][i].x1 + J[k+2][j][i].x1) * dz41);

      /* x2 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k-1][j][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k-1][j][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[k][j][i].x2 += eta_6 * ((J[k][j][i-2].x2 - 4.0*J[k][j][i-1].x2
         + 6.0*J[k][j][i].x2 - 4.0*J[k][j][i+1].x2 + J[k][j][i+2].x2) * dx41
                        + fac2 * (J[k][j-2][i].x2 - 4.0*J[k][j-1][i].x2
         + 6.0*J[k][j][i].x2 - 4.0*J[k][j+1][i].x2 + J[k][j+2][i].x2) * dy41
                        + fac3 * (J[k-2][j][i].x2 - 4.0*J[k-1][j][i].x2
         + 6.0*J[k][j][i].x2 - 4.0*J[k+1][j][i].x2 + J[k+2][j][i].x2) * dz41);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k][j-1][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k][j-1][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[k][j][i].x3 += eta_6 * ((J[k][j][i-2].x3 - 4.0*J[k][j][i-1].x3
         + 6.0*J[k][j][i].x3 - 4.0*J[k][j][i+1].x3 + J[k][j][i+2].x3) * dx41
                        + fac2 * (J[k][j-2][i].x3 - 4.0*J[k][j-1][i].x3
         + 6.0*J[k][j][i].x3 - 4.0*J[k][j+1][i].x3 + J[k][j+2][i].x3) * dy41
                        + fac3 * (J[k-2][j][i].x3 - 4.0*J[k-1][j][i].x3
         + 6.0*J[k][j][i].x3 - 4.0*J[k+1][j][i].x3 + J[k+2][j][i].x3) * dz41);
    }}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* resistivity_init: Allocate temporary arrays
 */

void resistivity_init(MeshS *pM)
{
  int nl,nd,size1=0,size2=0,size3=0,Nx1,Nx2,Nx3;
  int mycase;

/* Assign the function pointer for diffusivity calculation */
  mycase = par_geti_def("problem","CASE",1);

  if (mycase == 1)
    /* standard (no small grain) prescription with constant coefficients */
    get_myeta = eta_standard; 
  else
    /* general prescription with user defined diffusivities */
    get_myeta = get_eta_user;

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

  Nx1 = size1 + 2*nghost;

  if (pM->Nx[1] > 1){
    Nx2 = size2 + 2*nghost;
  } else {
    Nx2 = size2;
  }

  if (pM->Nx[2] > 1){
    Nx3 = size3 + 2*nghost;
  } else {
    Nx3 = size3;
  }

  if ((J = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
    goto on_error;
  if ((emf=(Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
    goto on_error;
#ifndef BAROTROPIC
  if ((EnerFlux = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real3Vect)))
    == NULL) goto on_error;
#endif
  if (Q_Hall > 0.0) {
    if ((Bcor = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
    if ((Jcor = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
    if ((emfh = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
  }
#ifdef SHEARING_BOX
  if (pM->Nx[2] > 1){
    if ((emf2 = (Real***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real)))==NULL)
      goto on_error;
    if ((remapEyiib = (Real**)calloc_2d_array(Nx3,Nx2,sizeof(Real)))==NULL)
      goto on_error;
    if ((remapEyoib = (Real**)calloc_2d_array(Nx3,Nx2,sizeof(Real)))==NULL)
      goto on_error;
    if ((J2   = (Real***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real)))==NULL)
      goto on_error;
    if ((remapJyiib = (Real***)calloc_3d_array(nghost,Nx3,Nx2,sizeof(Real)))==NULL)
      goto on_error;
    if ((remapJyoib = (Real***)calloc_3d_array(nghost,Nx3,Nx2,sizeof(Real)))==NULL)
      goto on_error;
  }
#endif

  return;

  on_error:
  resistivity_destruct();
  ath_error("[resistivity_init]: malloc returned a NULL pointer\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* resistivity_destruct: Free temporary arrays
 */

void resistivity_destruct()
{
  get_myeta = NULL;

  if (J != NULL) free_3d_array(J);
  if (emf != NULL) free_3d_array(emf);
#ifndef BAROTROPIC
  if (EnerFlux != NULL) free_3d_array(EnerFlux);
#endif

  if (Bcor != NULL) free_3d_array(Bcor);
  if (Jcor != NULL) free_3d_array(Jcor);
  if (emfh != NULL) free_3d_array(emfh);

#ifdef SHEARING_BOX
  if (emf2 != NULL) free_3d_array(emf2);
  if (remapEyiib != NULL) free_2d_array(remapEyiib);
  if (remapEyoib != NULL) free_2d_array(remapEyoib);
  if (J2 != NULL)   free_3d_array(J2);
  if (remapJyiib != NULL) free_3d_array(remapJyiib);
  if (remapJyoib != NULL) free_3d_array(remapJyoib);
#endif

  return;
}

#endif /* RESISTIVITY */
