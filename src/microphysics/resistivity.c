#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Adds explicit resistivity terms to the induction and energy eqns,
 *      dB/dt = -Curl(E)
 *      dE/dt = Div(B X E)
 *   where E = eta_Ohm J +  = electric field (emf) 
 *         J = Curl(B) = current density
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

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   EField_Ohm  - computes electric field due to Ohmic dissipation
 *   EField_Hall - computes electric field due to Hall effect
 *   EField_AD   - computes electric field due to ambipolar diffusion
 *============================================================================*/

void EField_Ohm(DomainS *pD);
void EField_Hall(DomainS *pD);
void EField_AD(DomainS *pD);

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
    jl = js - 1;
    ju = je + 2;
    dtodx2 = pG->dt/pG->dx2;
    ndim++;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;
    ku = ke + 2;
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

/*--- Step 1. Compute currents.------------------------------------------------
 * Note:  J1 = (dB3/dx2 - dB2/dx3)
 *        J2 = (dB1/dx3 - dB3/dx1)
 *        J3 = (dB2/dx1 - dB1/dx2) */

/* 1D PROBLEM */
  if (ndim == 1){
    for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-1; i<=ie+2; i++) {
        J[k][j][i].x = 0.0;
        J[k][j][i].y = -(pG->U[k][j][i].B3c - pG->U[k][j][i-1].B3c)/pG->dx1;
        J[k][j][i].z =  (pG->U[k][j][i].B2c - pG->U[k][j][i-1].B2c)/pG->dx1;
      }
    }}
  }

/* 2D PROBLEM */
  if (ndim == 2){
    for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-1; i<=ie+2; i++) {
        J[k][j][i].x =  (pG->U[k][j][i].B3c - pG->U[k][j-1][i  ].B3c)/pG->dx2;
        J[k][j][i].y = -(pG->U[k][j][i].B3c - pG->U[k][j  ][i-1].B3c)/pG->dx1;
        J[k][j][i].z =  (pG->B2i[k][j][i] - pG->B2i[k][j  ][i-1])/pG->dx1 -
                        (pG->B1i[k][j][i] - pG->B1i[k][j-1][i  ])/pG->dx2;
      }
    }}
  }

/* 3D PROBLEM */
  if (ndim == 3){
    for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-1; i<=ie+2; i++) {
        J[k][j][i].x = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                       (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        J[k][j][i].y = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                       (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        J[k][j][i].z = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                       (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
      }
    }}
  }

/*--- Step 2.  Call functions to compute resistive EMFs ------------------------
 * including Ohmic dissipation, the Hall effect, and ambipolar diffusion.
 * Current density (J) and emfs are global variables in this file. */

  if (eta_Ohm > 0.0)  EField_Ohm(pD);
  if (eta_Hall > 0.0) EField_Hall(pD);
  if (eta_AD > 0.0)   EField_AD(pD);

/*--- Step 3. CT update of magnetic field -------------------------------------
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

#ifndef BAROTROPIC
/*--- Step 4.  Compute energy fluxes -------------------------------------------
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

/*--- Step 5.  Update total energy ---------------------------------------------
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

  return;
}

/*----------------------------------------------------------------------------*/
/* EField_Ohm:  Resistive EMF from Ohmic dissipation.   E = \eta_Ohm J
 */

void EField_Ohm(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;

  if (pG->Nx[1] > 1){
    jl = js - 1;
    ju = je + 2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
    kl = ks - 1;
    ku = ke + 2;
  } else {
    kl = ks;
    ku = ke;
  }

/* For Ohmic resistivity, E = \eta_Ohm J  */
  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
    for (i=is; i<=ie+1; i++) {
      emf[k][j][i].x += eta_Ohm*J[k][j][i].x;
      emf[k][j][i].y += eta_Ohm*J[k][j][i].y;
      emf[k][j][i].z += eta_Ohm*J[k][j][i].z;
    }
  }}

  return;
}

/*----------------------------------------------------------------------------*/
/* EField_Hall:  Resistive EMF from Hall effect.  E = Q_H (J X B)
 */

void EField_Hall(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ndim=1;
  Real Q_H;

  if (pG->Nx[1] > 1) ndim++;
  if (pG->Nx[2] > 1) ndim++;

/* 1D PROBLEM:
 *   emf.x =  0.0
 *   emf.y =  Jz*Bx 
 *   emf.z = -Jy*Bx   */

  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {
      Q_H = 2.0*eta_Hall/(pG->U[ks][js][i].d + pG->U[ks][js][i-1].d);
      emf[ks][js][i].x += 0.0;
      emf[ks][js][i].y += Q_H*J[ks][js][i].z*pG->B1i[ks][js][i];
      emf[ks][js][i].z -= Q_H*J[ks][js][i].y*pG->B1i[ks][js][i];
    }
  }

/* 2D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 2){
    for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {

      Q_H = 2.0*eta_Hall/(pG->U[ks][j][i].d + pG->U[ks][j-1][i].d);

      emf[ks][j][i].x += Q_H*(  
        0.25*((J[ks][j  ][i].y + J[ks][j  ][i+1].y)*pG->U[ks][j  ][i].B3c +
              (J[ks][j-1][i].y + J[ks][j-1][i+1].y)*pG->U[ks][j-1][i].B3c) -
         0.5*((J[ks][j  ][i].z + J[ks][j  ][i+1].z)*pG->B2i[ks][j  ][i]) );

      Q_H = 2.0*eta_Hall/(pG->U[ks][j][i].d + pG->U[ks][j][i-1].d);

      emf[ks][j][i].y +=  Q_H*(
         0.5*((J[ks][j][i  ].z + J[ks][j+1][i  ].z)*pG->B1i[ks][j][i  ]) -
        0.25*((J[ks][j][i  ].x + J[ks][j+1][i  ].x)*pG->U[ks][j][i  ].B3c +
              (J[ks][j][i-1].x + J[ks][j+1][i-1].x)*pG->U[ks][j][i-1].B3c) );
    
      Q_H = 4.0*eta_Hall/(pG->U[ks][j  ][i].d + pG->U[ks][j  ][i-1].d
                        + pG->U[ks][j-1][i].d + pG->U[ks][j-1][i-1].d);

      emf[ks][j][i].z += Q_H*(
        0.5*(J[ks][j  ][i  ].x*pG->B2i[ks][j  ][i  ] +
             J[ks][j  ][i-1].x*pG->B2i[ks][j  ][i-1]) -
        0.5*(J[ks][j  ][i  ].y*pG->B1i[ks][j  ][i  ] +
             J[ks][j-1][i  ].y*pG->B1i[ks][j-1][i  ]) );
    }}
  }

/* 3D PROBLEM:
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx  */

  if (ndim == 3){
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {

        Q_H = 4.0*eta_Hall/(pG->U[k  ][j][i].d + pG->U[k  ][j-1][i].d
                          + pG->U[k-1][j][i].d + pG->U[k-1][j-1][i].d);

        emf[k][j][i].x += Q_H*(
          0.25*((J[k  ][j  ][i].y + J[k  ][j  ][i+1].y)*pG->B3i[k  ][j  ][i] +
                (J[k  ][j-1][i].y + J[k  ][j-1][i+1].y)*pG->B3i[k  ][j-1][i])-
          0.25*((J[k  ][j  ][i].z + J[k  ][j  ][i+1].z)*pG->B2i[k  ][j  ][i] +
                (J[k-1][j  ][i].z + J[k-1][j  ][i+1].z)*pG->B2i[k-1][j  ][i]) );

        Q_H = 4.0*eta_Hall/(pG->U[k  ][j][i].d + pG->U[k  ][j][i-1].d
                          + pG->U[k-1][j][i].d + pG->U[k-1][j][i-1].d);

        emf[k][j][i].y += Q_H*(
          0.25*((J[k  ][j][i  ].z + J[k  ][j+1][i  ].z)*pG->B1i[k  ][j][i  ] +
                (J[k-1][j][i  ].z + J[k-1][j+1][i  ].z)*pG->B1i[k-1][j][i  ])-
          0.25*((J[k  ][j][i  ].x + J[k  ][j+1][i  ].x)*pG->B3i[k  ][j][i  ] +
                (J[k  ][j][i-1].x + J[k  ][j+1][i-1].x)*pG->B3i[k  ][j][i-1]) );

        Q_H = 4.0*eta_Hall/(pG->U[k][j  ][i].d + pG->U[k][j  ][i-1].d
                          + pG->U[k][j-1][i].d + pG->U[k][j-1][i-1].d);

        emf[k][j][i].z += Q_H*(
          0.25*((J[k][j  ][i  ].x + J[k+1][j  ][i  ].x)*pG->B2i[k][j  ][i  ] +
                (J[k][j  ][i-1].x + J[k+1][j  ][i-1].x)*pG->B2i[k][j  ][i-1])-
          0.25*((J[k][j  ][i  ].y + J[k+1][j  ][i  ].y)*pG->B1i[k][j  ][i  ] +
                (J[k][j-1][i  ].y + J[k+1][j-1][i  ].y)*pG->B1i[k][j-1][i  ]) );
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
  Real Q_AD;

  if (pG->Nx[1] > 1) ndim++;
  if (pG->Nx[2] > 1) ndim++;

/* 1D PROBLEM: */

  if (ndim == 1){
    for (i=is; i<=ie+1; i++) {
      Q_AD = eta_AD;
      emf[ks][js][i].x += 0.0;
      emf[ks][js][i].y += Q_AD*J[ks][js][i].z*pG->B1i[ks][js][i];
      emf[ks][js][i].z -= Q_AD*J[ks][js][i].y*pG->B1i[ks][js][i];
    }
  }

/* 2D PROBLEM: */

/* 3D PROBLEM: */

  return;
}

/*----------------------------------------------------------------------------*/
/* resistivity_init: Allocate temporary arrays
 */

void resistivity_init(MeshS *pM)
{
  int nl,nd,size1=0,size2=0,size3=0,Nx1,Nx2,Nx3;

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
  return;

  on_error:
  resistivity_destruct();
  ath_error("[resisticvity_init]: malloc returned a NULL pointer\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* resistivity_destruct: Free temporary arrays
 */

void resistivity_destruct(void)
{
  if (J != NULL) free_3d_array(J);
  if (emf != NULL) free_3d_array(emf);
#ifndef BAROTROPIC
  if (EnerFlux != NULL) free_3d_array(EnerFlux);
#endif
  return;
}
#endif /* RESISTIVITY */
