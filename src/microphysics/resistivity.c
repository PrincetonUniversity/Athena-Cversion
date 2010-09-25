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
static Real3Vect ***J=NULL, ***emf=NULL, ***EnerFlux=NULL;

/* current, Bi in the corrector step & emf in the predict step for Hall MHD */
static Real3Vect ***Jcor=NULL, ***Bcor=NULL, ***emfp=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   EField_Ohm  - computes electric field due to Ohmic dissipation
 *   EField_Hall - computes electric field due to Hall effect
 *   EField_AD   - computes electric field due to ambipolar diffusion
 *   hyper_diffusion? -- compute higher-order derivatives of current
 *============================================================================*/

void EField_Ohm(DomainS *pD);
void EField_Hall(DomainS *pD);
void EField_AD(DomainS *pD);

void EField_Hall_sub(DomainS *pD, Real3Vect ***Bs, Real3Vect ***Js,
                                  Real3Vect ***emfs, int noff);
void hyper_diffusion6(DomainS *pD);

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
  int ndim=1;
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = 0.0, dtodx3 = 0.0;

  if (pG->Nx[1] > 1){
    jl = js - 3;
    ju = je + 4;
    dtodx2 = pG->dt/pG->dx2;
    ndim++;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 3;
    ku = ke + 4;
    dtodx3 = pG->dt/pG->dx3;
    ndim++;
  } else {
    kl = ks;
    ku = ke;
  }

/* zero fluxes (electric fields) */

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
    for (i=is-1; i<=ie+1; i++) {
      emf[k][j][i].x = 0.0;
      emf[k][j][i].y = 0.0;
      emf[k][j][i].z = 0.0;
    }
  }}

  if (Q_Hall > 0.0) {
    for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-3; i<=ie+4; i++) {
        emfp[k][j][i].x = 0.0;
        emfp[k][j][i].y = 0.0;
        emfp[k][j][i].z = 0.0;
      }
    }}
  }

/*--- Step 1. Compute currents.------------------------------------------------
 * Note:  J1 = (dB3/dx2 - dB2/dx3)
 *        J2 = (dB1/dx3 - dB3/dx1)
 *        J3 = (dB2/dx1 - dB1/dx2) */

/* 1D PROBLEM */
  if (ndim == 1){
    for (i=is-3; i<=ie+4; i++) {
      J[ks][js][i].x = 0.0;
      J[ks][js][i].y = -(pG->U[ks][js][i].B3c - pG->U[ks][js][i-1].B3c)/pG->dx1;
      J[ks][js][i].z =  (pG->U[ks][js][i].B2c - pG->U[ks][js][i-1].B2c)/pG->dx1;
    }
  }

/* 2D PROBLEM */
  if (ndim == 2){
    for (j=js-3; j<=je+4; j++) {
      for (i=is-3; i<=ie+4; i++) {
        J[ks][j][i].x=  (pG->U[ks][j][i].B3c - pG->U[ks][j-1][i  ].B3c)/pG->dx2;
        J[ks][j][i].y= -(pG->U[ks][j][i].B3c - pG->U[ks][j  ][i-1].B3c)/pG->dx1;
        J[ks][j][i].z=  (pG->B2i[ks][j][i] - pG->B2i[ks][j  ][i-1])/pG->dx1 -
                        (pG->B1i[ks][j][i] - pG->B1i[ks][j-1][i  ])/pG->dx2;
      }
    }
  }

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=ks-3; k<=ke+4; k++) {
     for (j=js-3; j<=je+4; j++) {
      for (i=is-3; i<=ie+4; i++) {
        J[k][j][i].x = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                       (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        J[k][j][i].y = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                       (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        J[k][j][i].z = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                       (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
      }
     }
    }
  }

/*--- Step 2.  Call functions to compute resistive EMFs ------------------------
 * including Ohmic dissipation, the Hall effect, and ambipolar diffusion.
 * Current density (J) and emfs are global variables in this file. */

  if (eta_Ohm > 0.0)  EField_Ohm(pD);
  if (Q_Hall > 0.0) EField_Hall(pD);
  if (Q_AD > 0.0)   EField_AD(pD);

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
      EnerFlux[ks][js][i].x =
         0.5*(pG->U[ks][js][i].B2c + pG->U[ks][js][i-1].B2c)*emf[ks][js][i].z
       - 0.5*(pG->U[ks][js][i].B3c + pG->U[ks][js][i-1].B3c)*emf[ks][js][i].y;
    }
  } 
      
/* 2D PROBLEM */
  if (ndim == 2){
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      EnerFlux[ks][j][i].x = 0.25*(pG->U[ks][j][i].B2c + pG->U[ks][j][i-1].B2c)*
                            (emf[ks][j][i].z + emf[ks][j+1][i].z)
         - 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j][i-1].B3c)*emf[ks][j][i].y;
    }}
    
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      EnerFlux[ks][j][i].y =
         0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c)*emf[ks][j][i].x
       - 0.25*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c)*
                (emf[ks][j][i].z + emf[ks][j][i+1].z);
    }}
  }   

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        EnerFlux[k][j][i].x = 0.25*(pG->U[k][j][i].B2c + pG->U[k][j][i-1].B2c)*
                             (emf[k][j][i].z + emf[k][j+1][i].z)
                            - 0.25*(pG->U[k][j][i].B3c + pG->U[k][j][i-1].B3c)*
                             (emf[k][j][i].y + emf[k+1][j][i].y);
      }
    }}

    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        EnerFlux[k][j][i].y = 0.25*(pG->U[k][j][i].B3c + pG->U[k][j-1][i].B3c)*
                             (emf[k][j][i].x + emf[k+1][j][i].x)
                            - 0.25*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c)*
                             (emf[k][j][i].z + emf[k][j][i+1].z);
      }
    }}

    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        EnerFlux[k][j][i].z = 0.25*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c)*
                             (emf[k][j][i].y + emf[k][j][i+1].y)
                            - 0.25*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c)*
                             (emf[k][j][i].x + emf[k][j+1][i].x);
      }
    }}
  }

/*--- Step 4.  Update total energy ---------------------------------------------
 * Update energy using x1-fluxes */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[k][j][i].E += dtodx1*(EnerFlux[k][j][i+1].x - EnerFlux[k][j][i].x);
    }
  }}

/* Update energy using x2-fluxes */

  if (pG->Nx[1] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx2*(EnerFlux[k][j+1][i].y -EnerFlux[k][j][i].y);
      }
    }}
  }

/* Update energy using x3-fluxes */

  if (pG->Nx[2] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(EnerFlux[k+1][j][i].z -EnerFlux[k][j][i].z);
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
      pG->U[ks][js][i].B2c += dtodx1*(emf[ks][js][i+1].z - emf[ks][js][i].z);
      pG->U[ks][js][i].B3c -= dtodx1*(emf[ks][js][i+1].y - emf[ks][js][i].y);
/* For consistency, set B2i and B3i to cell-centered values. */
      pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
      pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
    }
  }

/* 2D PROBLEM: CT +  centered differences for B3c */
  if (ndim == 2){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->B1i[ks][j][i] -= dtodx2*(emf[ks][j+1][i  ].z - emf[ks][j][i].z);
        pG->B2i[ks][j][i] += dtodx1*(emf[ks][j  ][i+1].z - emf[ks][j][i].z);

        pG->U[ks][j][i].B3c += dtodx2*(emf[ks][j+1][i  ].x - emf[ks][j][i].x) -
                               dtodx1*(emf[ks][j  ][i+1].y - emf[ks][j][i].y);
      }
      pG->B1i[ks][j][ie+1] -= dtodx2*(emf[ks][j+1][ie+1].z -emf[ks][j][ie+1].z);
    }
    for (i=is; i<=ie; i++) {
      pG->B2i[ks][je+1][i] += dtodx1*(emf[ks][je+1][i+1].z -emf[ks][je+1][i].z);
    } 
/* Set cell centered magnetic fields to average of face centered */
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j][i+1]);
        pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i] + pG->B2i[ks][j+1][i]);
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
          pG->B1i[k][j][i] += dtodx3*(emf[k+1][j  ][i  ].y - emf[k][j][i].y) -
                              dtodx2*(emf[k  ][j+1][i  ].z - emf[k][j][i].z);
          pG->B2i[k][j][i] += dtodx1*(emf[k  ][j  ][i+1].z - emf[k][j][i].z) -
                              dtodx3*(emf[k+1][j  ][i  ].x - emf[k][j][i].x);
          pG->B3i[k][j][i] += dtodx2*(emf[k  ][j+1][i  ].x - emf[k][j][i].x) -
                              dtodx1*(emf[k  ][j  ][i+1].y - emf[k][j][i].y);
        }
        pG->B1i[k][j][ie+1] +=
          dtodx3*(emf[k+1][j  ][ie+1].y - emf[k][j][ie+1].y) -
          dtodx2*(emf[k  ][j+1][ie+1].z - emf[k][j][ie+1].z);
      }
      for (i=is; i<=ie; i++) {
        pG->B2i[k][je+1][i] +=
          dtodx1*(emf[k  ][je+1][i+1].z - emf[k][je+1][i].z) -
          dtodx3*(emf[k+1][je+1][i  ].x - emf[k][je+1][i].x);
      }
    }
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B3i[ke+1][j][i] +=
        dtodx2*(emf[ke+1][j+1][i  ].x - emf[ke+1][j][i].x) -
        dtodx1*(emf[ke+1][j  ][i+1].y - emf[ke+1][j][i].y);
    }}
/* Set cell centered magnetic fields to average of face centered */
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
      pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
      pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
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

  if (pG->Nx[1] > 1)    ndim++;
  if (pG->Nx[2] > 1)    ndim++;

/* For Ohmic resistivity, E = \eta_Ohm J  */

/* 1D PROBLEM: */
  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {

      eta_O = 0.5*(pG->eta_Ohm[ks][js][i] + pG->eta_Ohm[ks][js][i-1]);

      emf[ks][js][i].y += eta_O * J[ks][js][i].y;
      emf[ks][js][i].z += eta_O * J[ks][js][i].z;
    }
  }

/* 2D PROBLEM: */
  if (ndim == 2){
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      eta_O = 0.5*(pG->eta_Ohm[ks][j][i] + pG->eta_Ohm[ks][j-1][i]);

      emf[ks][j][i].x += eta_O * J[ks][j][i].x;

      eta_O = 0.5*(pG->eta_Ohm[ks][j][i] + pG->eta_Ohm[ks][j][i-1]);

      emf[ks][j][i].y += eta_O * J[ks][j][i].y; 

      eta_O = 0.25*(pG->eta_Ohm[ks][j][i  ] + pG->eta_Ohm[ks][j-1][i  ] +
                    pG->eta_Ohm[ks][j][i-1] + pG->eta_Ohm[ks][j-1][i-1]);

      emf[ks][j][i].z += eta_O * J[ks][j][i].z;
    }}
  }


/* 3D PROBLEM: */
  
  if (ndim == 3){
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {

        eta_O = 0.25*(pG->eta_Ohm[k][j  ][i] + pG->eta_Ohm[k-1][j  ][i] +
                      pG->eta_Ohm[k][j-1][i] + pG->eta_Ohm[k-1][j-1][i]);

        emf[k][j][i].x += eta_O * J[k][j][i].x;

        eta_O = 0.25*(pG->eta_Ohm[k][j][i  ] + pG->eta_Ohm[k-1][j][i  ] +
                      pG->eta_Ohm[k][j][i-1] + pG->eta_Ohm[k-1][j][i-1]);

        emf[k][j][i].y += eta_O * J[k][j][i].y;

        eta_O = 0.25*(pG->eta_Ohm[k][j][i  ] + pG->eta_Ohm[k][j-1][i  ] +
                      pG->eta_Ohm[k][j][i-1] + pG->eta_Ohm[k][j-1][i-1]);

        emf[k][j][i].z += eta_O * J[k][j][i].z;
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
  int ndim=1;
  Real eta_H, Bmag;
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = 0.0, dtodx3 = 0.0;

  il = is - 3;  iu = ie + 4;

  if (pG->Nx[1] > 1){
    jl = js - 3;    ju = je + 4;
    dtodx2 = pG->dt/pG->dx2;
    ndim++;
  } else {
    jl = js;        ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 3;    ku = ke + 4;
    dtodx3 = pG->dt/pG->dx3;
    ndim++;
  } else {
    kl = ks;        ku = ke;
  }

/* Preliminary: add hyper-diffusion */

  hyper_diffusion6(pD);

/* Preliminary: divide eta_Hall by B for convenience */
  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=il; i<=iu; i++) {

      Bmag = sqrt(SQR(pG->U[k][j][i].B1c)
                + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));

      pG->eta_Hall[k][j][i] /= Bmag+TINY_NUMBER;

      Bcor[k][j][i].x = pG->B1i[k][j][i];
      Bcor[k][j][i].y = pG->B2i[k][j][i];
      Bcor[k][j][i].z = pG->B3i[k][j][i];
  }}}

/* Step 1: Predictor step for emf calculation at t=t(n)  */
  EField_Hall_sub(pD, Bcor, J, emfp, 3);

/* Step 2: Calculate Bcor and Jcor at t=t(n+1/2) */

  /* 1D PROBLEM: centered differences for B2c and B3c */
  if (ndim == 1){
    for (i=is-2; i<=ie+2; i++) {
      Bcor[ks][js][i].y += 0.5*dtodx1*(emfp[ks][js][i+1].z - emfp[ks][js][i].z);
      Bcor[ks][js][i].z -= 0.5*dtodx1*(emfp[ks][js][i+1].y - emfp[ks][js][i].y);
    }
  }

  /* 2D PROBLEM: CT +  centered differences for B3c */
  if (ndim == 2){
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Bcor[ks][j][i].x -= 0.5*dtodx2*(emfp[ks][j+1][i  ].z - emfp[ks][j][i].z);
        Bcor[ks][j][i].y += 0.5*dtodx1*(emfp[ks][j  ][i+1].z - emfp[ks][j][i].z);

        Bcor[ks][j][i].z += pG->eta_Hall[ks][j][i]*(
                            0.5*dtodx2*(emfp[ks][j+1][i  ].x - emfp[ks][j][i].x) -
                            0.5*dtodx1*(emfp[ks][j  ][i+1].y - emfp[ks][j][i].y));
      }
    }
  }
  
  if (ndim == 3){
    for (k=ks-2; k<=ke+2; k++) {
      for (j=js-2; j<=je+2; j++) {
        for (i=is-2; i<=ie+2; i++) {
          Bcor[k][j][i].x += 0.5*dtodx3*(emfp[k+1][j  ][i  ].y - emfp[k][j][i].y) -
                             0.5*dtodx2*(emfp[k  ][j+1][i  ].z - emfp[k][j][i].z);
          Bcor[k][j][i].y += 0.5*dtodx1*(emfp[k  ][j  ][i+1].z - emfp[k][j][i].z) -
                             0.5*dtodx3*(emfp[k+1][j  ][i  ].x - emfp[k][j][i].x);
          Bcor[k][j][i].z += 0.5*dtodx2*(emfp[k  ][j+1][i  ].x - emfp[k][j][i].x) -
                             0.5*dtodx1*(emfp[k  ][j  ][i+1].y - emfp[k][j][i].y);
        }
      }
    }
  }

/* Step 3: compute current for the corrector step */

/* 1D PROBLEM */
  if (ndim == 1){
    for (i=is-1; i<=ie+2; i++) {
      Jcor[ks][js][i].x = 0.0;
      Jcor[ks][js][i].y = -(Bcor[ks][js][i].z - Bcor[ks][js][i-1].z)/pG->dx1;
      Jcor[ks][js][i].z =  (Bcor[ks][js][i].y - Bcor[ks][js][i-1].y)/pG->dx1;
    }
  }

/* 2D PROBLEM */
  if (ndim == 2){
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
        Jcor[ks][j][i].x=  (Bcor[ks][j][i].z - Bcor[ks][j-1][i  ].z)/pG->dx2;
        Jcor[ks][j][i].y= -(Bcor[ks][j][i].z - Bcor[ks][j  ][i-1].z)/pG->dx1;
        Jcor[ks][j][i].z=  (Bcor[ks][j][i].y - Bcor[ks][j  ][i-1].y)/pG->dx1 -
                           (Bcor[ks][j][i].x - Bcor[ks][j-1][i  ].x)/pG->dx2;
      }
    }
  }

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
        Jcor[k][j][i].x = (Bcor[k][j][i].z - Bcor[k  ][j-1][i  ].z)/pG->dx2 -
                          (Bcor[k][j][i].y - Bcor[k-1][j  ][i  ].y)/pG->dx3;
        Jcor[k][j][i].y = (Bcor[k][j][i].x - Bcor[k-1][j  ][i  ].x)/pG->dx3 -
                          (Bcor[k][j][i].z - Bcor[k  ][j  ][i-1].z)/pG->dx1;
        Jcor[k][j][i].z = (Bcor[k][j][i].y - Bcor[k  ][j  ][i-1].y)/pG->dx1 -
                          (Bcor[k][j][i].x - Bcor[k  ][j-1][i  ].x)/pG->dx2;
      }
    }}
  }

/* Step 4: Corrector step fir emf  calculation at t=t(n+1/2) */
  EField_Hall_sub(pD, Bcor, Jcor, emf, 1);

  return;

}

/* Update the Hall emfs (subroutine for the Hall calculation) */
void EField_Hall_sub(DomainS *pD, Real3Vect ***Bs, Real3Vect ***Js,
                                  Real3Vect ***emfs, int noff)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
  Real eta_H;

  if (pG->Nx[1] > 1)  ndim++;
  if (pG->Nx[2] > 1)  ndim++;

/* 1D PROBLEM:
 *   emf.x =  0.0
 *   emf.y =  Jz*Bx
 *   emf.z = -Jy*Bx   */

  if (ndim == 1){
    for (i=is-noff+1; i<=ie+noff; i++) {
      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);

      emfs[ks][js][i].y += eta_H*Js[ks][js][i].z * Bs[ks][js][i].x;
      emfs[ks][js][i].z -= eta_H*Js[ks][js][i].y * Bs[ks][js][i].x;
    }
  }

/* 2D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 2){
    for (j=js-noff+1; j<=je+noff; j++) {
    for (i=is-noff+1; i<=ie+noff; i++) {

      /* x1 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j-1][i]);

      emfs[ks][j][i].x += eta_H*(
        0.125*(Js[ks][j  ][i].y + Js[ks][j  ][i+1].y
             + Js[ks][j-1][i].y + Js[ks][j-1][i+1].y)
             *(Bs[ks][j  ][i].z + Bs[ks][j-1][i  ].z) -
         0.5*((Js[ks][j  ][i].z + Js[ks][j  ][i+1].z)*Bs[ks][j  ][i].y));

      /* x2 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j][i-1]);

      emfs[ks][j][i].y += eta_H*(
         0.5*((Js[ks][j][i  ].z + Js[ks][j+1][i  ].z)*Bs[ks][j][i  ].x) -
        0.125*(Js[ks][j][i  ].x + Js[ks][j+1][i  ].x
             + Js[ks][j][i-1].x + Js[ks][j+1][i-1].x)
             *(Bs[ks][j][i  ].z + Bs[ks][j  ][i-1].z));
    
      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[ks][j][i  ] + pG->eta_Hall[ks][j-1][i  ] +
                    pG->eta_Hall[ks][j][i-1] + pG->eta_Hall[ks][j-1][i-1]);

      emfs[ks][j][i].z += eta_H*(
        0.25*(Js[ks][j][i].x + Js[ks][j][i-1].x)
            *(Bs[ks][j][i].y + Bs[ks][j][i-1].y) -
        0.25*(Js[ks][j][i].y + Js[ks][i-1][i].y)
            *(Bs[ks][j][i].x + Bs[ks][j-1][i].x));
    }}
  }

/* 3D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 3){
    for (k=ks-noff+1; k<=ke+noff; k++) {
    for (j=js-noff+1; j<=je+noff; j++) {
      for (i=is-noff+1; i<=ie+noff; i++) {

        /* x1 */
        eta_H = 0.25*(pG->eta_Hall[k][j  ][i] + pG->eta_Hall[k-1][j  ][i] +
                      pG->eta_Hall[k][j-1][i] + pG->eta_Hall[k-1][j-1][i]);

        emfs[k][j][i].x += 0.125*eta_H*(
                (Js[k  ][j  ][i].y + Js[k  ][j  ][i+1].y
               + Js[k  ][j-1][i].y + Js[k  ][j-1][i+1].y)
               *(Bs[k  ][j  ][i].z + Bs[k  ][j-1][i  ].z)-
                (Js[k  ][j  ][i].z + Js[k  ][j  ][i+1].z
               + Js[k-1][j  ][i].z + Js[k-1][j  ][i+1].z)
               *(Bs[k  ][j  ][i].y + Bs[k-1][j  ][i  ].y));

        /* x2 */
        eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k-1][j][i  ] +
                      pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k-1][j][i-1]);

        emfs[k][j][i].y += 0.125*eta_H*(
                (Js[k  ][j][i  ].z + Js[k  ][j+1][i  ].z
               + Js[k-1][j][i  ].z + Js[k-1][j+1][i  ].z)
               *(Bs[k  ][j][i  ].x + Bs[k-1][j  ][i  ].x)-
                (Js[k  ][j][i  ].x + Js[k  ][j+1][i  ].x
               + Js[k  ][j][i-1].x + Js[k  ][j+1][i-1].x)
               *(Bs[k  ][j][i  ].z + Bs[k  ][j  ][i-1].z));

        /* x3 */
        eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k][j-1][i  ] +
                      pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k][j-1][i-1]);

        emfs[k][j][i].z += 0.125*eta_H*(
                (Js[k][j  ][i  ].x + Js[k+1][j  ][i  ].x
               + Js[k][j  ][i-1].x + Js[k+1][j  ][i-1].x)
               *(Bs[k][j  ][i  ].y + Bs[k  ][j  ][i-1].y)-
                (Js[k][j  ][i  ].y + Js[k+1][j  ][i  ].y
               + Js[k][j-1][i  ].y + Js[k+1][j-1][i  ].y)
               *(Bs[k][j  ][i  ].x + Bs[  k][j-1][i  ].x));
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

  if (pG->Nx[1] > 1) ndim++;
  if (pG->Nx[2] > 1) ndim++;

/* 1D PROBLEM:
 *   emf.x = 0
 *   emf.y = (J_perp)_y
 *   emf.z = (J_perp)_z  */

  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {
      eta_A = 0.5*(pG->eta_AD[ks][js][i] + pG->eta_AD[ks][js][i-1]);

      intBx = pG->B1i[ks][js][i];
      intBy = 0.5*(pG->U[ks][js][i].B2c + pG->U[ks][js][i-1].B2c);
      intBz = 0.5*(pG->U[ks][js][i].B3c + pG->U[ks][js][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
      JdotB = J[ks][js][i].y*intBy + J[ks][js][i].z*intBz;

      emf[ks][js][i].y += eta_A*(J[ks][js][i].y - JdotB*intBy/Bsq);
      emf[ks][js][i].z += eta_A*(J[ks][js][i].z - JdotB*intBz/Bsq);
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
      eta_A = 0.5*(pG->eta_AD[ks][j][i] + pG->eta_AD[ks][j-1][i]);

      intJx = J[ks][j][i].x;
      intJy = 0.25*(J[ks][j  ][i].y + J[ks][j  ][i+1].y
                  + J[ks][j-1][i].y + J[ks][j-1][i+1].y);
      intJz = 0.5 *(J[ks][j][i].z   + J[ks][j][i+1].z);

      intBx = 0.5*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c);
      intBy = pG->B2i[ks][j][i];
      intBz = 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].x += eta_A*(J[ks][j][i].x - JdotB*intBx/Bsq);

      /* emf.y */
      eta_A = 0.5*(pG->eta_AD[ks][j][i] + pG->eta_AD[ks][j][i-1]);

      intJx = 0.25*(J[ks][j][i  ].x + J[ks][j+1][i  ].x
                  + J[ks][j][i-1].x + J[ks][j+1][i-1].x);
      intJy = J[ks][j][i].y;
      intJz = 0.5 *(J[ks][j][i].z   + J[ks][j+1][i].z);

      intBx = pG->B1i[ks][j][i];
      intBy = 0.5*(pG->U[ks][j][i].B2c + pG->U[ks][j][i-1].B2c);
      intBz = 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].y += eta_A*(J[ks][j][i].y - JdotB*intBy/Bsq);

     /* emf.z */
      eta_A = 0.25*(pG->eta_AD[ks][j  ][i] + pG->eta_AD[ks][j  ][i-1]
                  + pG->eta_AD[ks][j-1][i] + pG->eta_AD[ks][j-1][i-1]);

      intJx = 0.5*(J[ks][j][i].x + J[ks][j][i-1].x);
      intJy = 0.5*(J[ks][j][i].y + J[ks][j-1][i].y);
      intJz = J[ks][j][i].z;

      intBx = 0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j-1][i]);
      intBy = 0.5*(pG->B2i[ks][j][i] + pG->B2i[ks][j][i-1]);
      intBz = 0.25*(pG->U[ks][j  ][i].B3c + pG->U[ks][j  ][i-1].B3c
                  + pG->U[ks][j-1][i].B3c + pG->U[ks][j-1][i-1].B3c);

      Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
      JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

      emf[ks][j][i].z += eta_A*(J[ks][j][i].z - JdotB*intBz/Bsq);

    }}
  }

/* 3D PROBLEM:
 *   emf.x = (J_perp)_x
 *   emf.y = (J_perp)_y
 *   emf.z = (J_perp)_z  */

  if (ndim == 3){
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {

        /* emf.x */
        eta_A = 0.25*(pG->eta_AD[k][j  ][i] + pG->eta_AD[k-1][j  ][i] +
                      pG->eta_AD[k][j-1][i] + pG->eta_AD[k-1][j-1][i]);

        intJx = J[k][j][i].x;
        intJy = 0.25*(J[k][j  ][i].y + J[k][j  ][i+1].y
                    + J[k][j-1][i].y + J[k][j-1][i+1].y);
        intJz = 0.25*(J[k  ][j][i].z + J[k  ][j][i+1].z
                    + J[k-1][j][i].z + J[k-1][j][i+1].z);

        intBx = 0.25*(pG->U[k][j  ][i].B1c + pG->U[k-1][j  ][i].B1c +
                      pG->U[k][j-1][i].B1c + pG->U[k-1][j-1][i].B1c);
        intBy = 0.5*(pG->B2i[k][j][i] + pG->B2i[k-1][j][i]);
        intBz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k][j-1][i]);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].x += eta_A*(J[k][j][i].x - JdotB*intBx/Bsq);

        /* emf.y */
        eta_A = 0.25*(pG->eta_AD[k][j][i  ] + pG->eta_AD[k-1][j][i  ] +
                      pG->eta_AD[k][j][i-1] + pG->eta_AD[k-1][j][i-1]);

        intJx = 0.25*(J[k][j][i  ].x + J[k][j+1][i  ].x
                    + J[k][j][i-1].x + J[k][j+1][i-1].x);;
        intJy = J[k][j][i].y;
        intJz = 0.25*(J[k  ][j][i].z + J[k  ][j+1][i].z
                    + J[k-1][j][i].z + J[k-1][j+1][i].z);

        intBx = 0.5*(pG->B1i[k][j][i] + pG->B1i[k-1][j][i]);
        intBy = 0.25*(pG->U[k][j][i  ].B2c + pG->U[k-1][j][i  ].B2c +
                      pG->U[k][j][i-1].B2c + pG->U[k-1][j][i-1].B2c);
        intBz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k][j][i-1]);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].y += eta_A*(J[k][j][i].y - JdotB*intBy/Bsq);

        /* emf.z */
        eta_A = 0.25*(pG->eta_AD[k][j][i  ] + pG->eta_AD[k][j-1][i  ] +
                      pG->eta_AD[k][j][i-1] + pG->eta_AD[k][j-1][i-1]);

        intJx = 0.25*(J[k][j][i  ].x + J[k+1][j][i  ].x
                    + J[k][j][i-1].x + J[k+1][j][i-1].x);;
        intJy = 0.25*(J[k][j  ][i].y + J[k+1][j  ][i].y
                    + J[k][j-1][i].y + J[k+1][j-1][i].y);
        intJz = J[k][j][i].z;

        intBx = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j-1][i]);
        intBy = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j][i-1]);
        intBz = 0.25*(pG->U[k][j  ][i].B3c + pG->U[k][j  ][i-1].B3c +
                      pG->U[k][j-1][i].B3c + pG->U[k][j-1][i-1].B3c);

        Bsq = SQR(intBx) + SQR(intBy) + SQR(intBz);
        JdotB = intJx*intBx + intJy*intBy + intJz*intBz;

        emf[k][j][i].z += eta_A*(J[k][j][i].z - JdotB*intBz/Bsq);
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* hyper_diffusion: calculate the higher-order derivatives of J
 */

/* 6th order diffusion */
void hyper_diffusion6(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
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
  fac = 2.0*SQR(pG->dt/pG->dx1)*pG->dt;
  
  /* 1D */ 
  if (ndim == 1) {
    for (i=is; i<=ie+1; i++) {

      eta_H = 0.5*(pG->eta_Hall[ks][js][i] + pG->eta_Hall[ks][js][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][js][i].y += eta_6 *  (J[ks][js][i-2].y - 4.0*J[ks][js][i-1].y
         + 6.0*J[ks][js][i].y - 4.0*J[ks][js][i+1].y + J[ks][js][i+2].y) * dx41;
      emf[ks][js][i].z += eta_6 *  (J[ks][js][i-2].z - 4.0*J[ks][js][i-1].z
         + 6.0*J[ks][js][i].z - 4.0*J[ks][js][i+1].z + J[ks][js][i+2].z) * dx41;
    }
  }
  
  /* 2D */ 
  if (ndim == 2) {
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      /* x1 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j-1][i]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].x += eta_6 * ((J[ks][j][i-2].x - 4.0*J[ks][j][i-1].x
         + 6.0*J[ks][j][i].x - 4.0*J[ks][j][i+1].x + J[ks][j][i+2].x)*dx41
                         + fac2 * (J[ks][j-2][i].x - 4.0*J[ks][j-1][i].x
         + 6.0*J[ks][j][i].x - 4.0*J[ks][j+1][i].x + J[ks][j+2][i].x) * dy41);

      /* x2 */
      eta_H = 0.5*(pG->eta_Hall[ks][j][i] + pG->eta_Hall[ks][j][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].y += eta_6 * ((J[ks][j][i-2].y - 4.0*J[ks][j][i-1].y
         + 6.0*J[ks][j][i].y - 4.0*J[ks][j][i+1].y + J[ks][j][i+2].y) * dx41
                         + fac2 * (J[ks][j-2][i].y - 4.0*J[ks][j-1][i].y
         + 6.0*J[ks][j][i].y - 4.0*J[ks][j+1][i].y + J[ks][j+2][i].y) * dy41);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[ks][j][i  ] + pG->eta_Hall[ks][j-1][i  ] +
                    pG->eta_Hall[ks][j][i-1] + pG->eta_Hall[ks][j-1][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[ks][j][i].z += eta_6 * ((J[ks][j][i-2].z - 4.0*J[ks][j][i-1].z
         + 6.0*J[ks][j][i].z - 4.0*J[ks][j][i+1].z + J[ks][j][i+2].z) * dx41
                         + fac2 * (J[ks][j-2][i].z - 4.0*J[ks][j-1][i].z
         + 6.0*J[ks][j][i].z - 4.0*J[ks][j+1][i].z + J[ks][j+2][i].z) * dy41);
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

      emf[k][j][i].x += eta_6 * ((J[k][j][i-2].x - 4.0*J[k][j][i-1].x
         + 6.0*J[k][j][i].x - 4.0*J[k][j][i+1].x + J[k][j][i+2].x) * dx41
                        + fac2 * (J[k][j-2][i].x - 4.0*J[k][j-1][i].x
         + 6.0*J[k][j][i].x - 4.0*J[k][j+1][i].x + J[k][j+2][i].x) * dy41
                        + fac3 * (J[k-2][j][i].x - 4.0*J[k-1][j][i].x
         + 6.0*J[k][j][i].x - 4.0*J[k+1][j][i].x + J[k+2][j][i].x) * dz41);

      /* x2 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k-1][j][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k-1][j][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[k][j][i].y += eta_6 * ((J[k][j][i-2].y - 4.0*J[k][j][i-1].y
         + 6.0*J[k][j][i].y - 4.0*J[k][j][i+1].y + J[k][j][i+2].y) * dx41
                        + fac2 * (J[k][j-2][i].y - 4.0*J[k][j-1][i].y
         + 6.0*J[k][j][i].y - 4.0*J[k][j+1][i].y + J[k][j+2][i].y) * dy41
                        + fac3 * (J[k-2][j][i].y - 4.0*J[k-1][j][i].y
         + 6.0*J[k][j][i].y - 4.0*J[k+1][j][i].y + J[k+2][j][i].y) * dz41);

      /* x3 */
      eta_H = 0.25*(pG->eta_Hall[k][j][i  ] + pG->eta_Hall[k][j-1][i  ] +
                    pG->eta_Hall[k][j][i-1] + pG->eta_Hall[k][j-1][i-1]);
      eta_6 = SQR(SQR(eta_H)) * fac;

      emf[k][j][i].z += eta_6 * ((J[k][j][i-2].z - 4.0*J[k][j][i-1].z
         + 6.0*J[k][j][i].z - 4.0*J[k][j][i+1].z + J[k][j][i+2].z) * dx41
                        + fac2 * (J[k][j-2][i].z - 4.0*J[k][j-1][i].z
         + 6.0*J[k][j][i].z - 4.0*J[k][j+1][i].z + J[k][j+2][i].z) * dy41
                        + fac3 * (J[k-2][j][i].z - 4.0*J[k-1][j][i].z
         + 6.0*J[k][j][i].z - 4.0*J[k+1][j][i].z + J[k+2][j][i].z) * dz41);
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

  switch (mycase)
  {
    /* single-ion prescription with constant coefficients */
    case 1:  get_myeta = eta_single_const; break;

    /* single-ion prescription with user defined diffusiviteis */
    case 2:  get_myeta = eta_single_user;  break;

    /* general prescription with user defined diffusivities */
    case 3:  get_myeta = eta_general_user; break;

    default: ath_error("[resistivity_init]: CASE must equal to 1, 2 or 3!\n");
             return;
  }

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
    if ((Jcor = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
    if ((Bcor = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
    if ((emfp = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1,sizeof(Real3Vect)))==NULL)
      goto on_error;
  }

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

  if (Jcor != NULL) free_3d_array(Jcor);
  if (Bcor != NULL) free_3d_array(Bcor);
  if (emfp != NULL) free_3d_array(emfp);

  return;
}

#endif /* RESISTIVITY */
