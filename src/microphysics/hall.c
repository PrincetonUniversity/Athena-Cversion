#include "../copyright.h"
/*==============================================================================
 * FILE: hall.c
 *
 * PURPOSE: Implements explicit Hall and Ohmic resistivity, that is
 *      dB/dt = -Curl(Q_H[J X B] + [\eta J])    where J=Curl(B).
 *      dE/dt = Div(B X [Q_H(JXB)+\eta J])
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  hall_resistivity_1d()
 *  hall_resistivity_2d()
 *  hall_resistivity_3d()
 *  hall_resistivity_init() - allocates memory needed
 *  hall_resistivity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef HALL_MHD
#ifdef HYDRO
#error : Hall resistivity only works for MHD.
#endif /* HYDRO */
#endif

#ifdef HALL_MHD
/* The Hall emfs, etc., contained in special structure */
static Real3Vect ***J=NULL, ***emf=NULL, ***EFlux=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* hall_resistivity_1d:
 */

void hall_resistivity_1d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1, Q_H;

/*--- Step 1 -------------------------------------------------------------------
 * Compute currents. Note:
 *   Jx = 0
 *   Jy = (-dB3/dx1)
 *   Jz = (dB2/dx1)
 * Jy and Jz use B3c and B2c respectively, and are centered at x1-interfaces
 */

  for (i=is-1; i<=ie+2; i++) {
    J[ks][js][i].x = 0.0;
    J[ks][js][i].y = -(pG->U[ks][js][i].B3c - pG->U[ks][js][i-1].B3c)/pG->dx1;
    J[ks][js][i].z =  (pG->U[ks][js][i].B2c - pG->U[ks][js][i-1].B2c)/pG->dx1;
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute Hall EMFs.  emf.x is never used, and so is not computed.
 * emf.y and emf.z are face-ctrd.
 *   emf.x =  0.0
 *   emf.y =  Jz*Bx 
 *   emf.z = -Jy*Bx
 */

  for (i=is; i<=ie+1; i++) {
    emf[ks][js][i].x = 0.0;
    emf[ks][js][i].y =  J[ks][js][i].z*pG->B1i[ks][js][i];
    emf[ks][js][i].z = -J[ks][js][i].y*pG->B1i[ks][js][i];

    Q_H = 2.0*eta_Hall/(pG->U[ks][js][i].d + pG->U[ks][js][i-1].d);
    emf[ks][js][i].y *= Q_H;
    emf[ks][js][i].z *= Q_H;
  }

/*--- Step 3 -------------------------------------------------------------------
 * Add Ohmic resistivity (\eta J) to Hall EMFs, if required.
 */

  if (eta_Ohm != 0.0) {
    for (i=is; i<=ie+1; i++) {
      emf[ks][js][i].y += eta_Ohm*J[ks][js][i].y;
      emf[ks][js][i].z += eta_Ohm*J[ks][js][i].z;
    }
  }

#ifndef BAROTROPIC
/*--- Step 4 -------------------------------------------------------------------
 * Compute fluxes of energy
 */

  for (i=is; i<=ie+1; i++) {
    EFlux[ks][js][i].x = 
       0.5*(pG->U[ks][js][i].B2c + pG->U[ks][js][i-1].B2c)*emf[ks][js][i].z
     - 0.5*(pG->U[ks][js][i].B3c + pG->U[ks][js][i-1].B3c)*emf[ks][js][i].y;
  }
#endif

/*--- Step 5 -------------------------------------------------------------------
 * CT update of magnetic field using Hall [+Ohmic] EMFs. In 1D, this reduces to
 * centered differences for the resistive fluxes of B2c and B3c
 */
  
  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].B2c += dtodx1*(emf[ks][js][i+1].z - emf[ks][js][i].z);
    pG->U[ks][js][i].B3c -= dtodx1*(emf[ks][js][i+1].y - emf[ks][js][i].y);
  }

/* For consistency, set B2i and B3i to cell-centered values. */
  
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
  }

#ifndef BAROTROPIC
/*--- Step 6 -------------------------------------------------------------------
 * Update energy using resistive fluxes in each dimension (dE/dt = Div(F))
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].E  += dtodx1*(EFlux[ks][js][i+1].x - EFlux[ks][js][i].x);
  }
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* hall_resistivity_2d:
 */

void hall_resistivity_2d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real Q_H;

/*--- Step 1 -------------------------------------------------------------------
 * Compute currents. Note:
 *   J1 = (dB3/dx2)
 *   J2 = (-dB3/dx1)
 *   J3 = (dB2/dx1 - dB1/dx2)
 * J1 and J2 use B3c, and in 2D are centered at x2- and x1-interfaces
 */

  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
      J[ks][j][i].x =  (pG->U[ks][j][i].B3c - pG->U[ks][j-1][i  ].B3c)/pG->dx2;
      J[ks][j][i].y = -(pG->U[ks][j][i].B3c - pG->U[ks][j  ][i-1].B3c)/pG->dx1;

      J[ks][j][i].z = (pG->B2i[ks][j][i] - pG->B2i[ks][j  ][i-1])/pG->dx1 -
                       (pG->B1i[ks][j][i] - pG->B1i[ks][j-1][i  ])/pG->dx2;
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute Hall EMFs
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx
 */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      emf[ks][j][i].x = 
        0.25*((J[ks][j  ][i].y + J[ks][j  ][i+1].y)*pG->U[ks][j  ][i].B3c +
              (J[ks][j-1][i].y + J[ks][j-1][i+1].y)*pG->U[ks][j-1][i].B3c) -
         0.5*((J[ks][j  ][i].z + J[ks][j  ][i+1].z)*pG->B2i[ks][j  ][i]);

      emf[ks][j][i].y = 
         0.5*((J[ks][j][i  ].z + J[ks][j+1][i  ].z)*pG->B1i[ks][j][i  ]) -
        0.25*((J[ks][j][i  ].x + J[ks][j+1][i  ].x)*pG->U[ks][j][i  ].B3c +
              (J[ks][j][i-1].x + J[ks][j+1][i-1].x)*pG->U[ks][j][i-1].B3c);

      emf[ks][j][i].z = 
        0.5*(J[ks][j  ][i  ].x*pG->B2i[ks][j  ][i  ] +
             J[ks][j  ][i-1].x*pG->B2i[ks][j  ][i-1]) -
        0.5*(J[ks][j  ][i  ].y*pG->B1i[ks][j  ][i  ] +
             J[ks][j-1][i  ].y*pG->B1i[ks][j-1][i  ]);

      Q_H = 2.0*eta_Hall/(pG->U[ks][j][i].d + pG->U[ks][j-1][i].d);
      emf[ks][j][i].x *= Q_H;

      Q_H = 2.0*eta_Hall/(pG->U[ks][j][i].d + pG->U[ks][j][i-1].d);
      emf[ks][j][i].y *= Q_H;

      Q_H = 4.0*eta_Hall/(pG->U[ks][j  ][i].d + pG->U[ks][j  ][i-1].d 
                        + pG->U[ks][j-1][i].d + pG->U[ks][j-1][i-1].d);
      emf[ks][j][i].z *= Q_H;
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Add Ohmic resistivity (\eta J) to Hall EMFs, if required.
 */

  if (eta_Ohm != 0.0) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf[ks][j][i].x += eta_Ohm*J[ks][j][i].x;
        emf[ks][j][i].y += eta_Ohm*J[ks][j][i].y;
        emf[ks][j][i].z += eta_Ohm*J[ks][j][i].z;
      }
    }
  }

#ifndef BAROTROPIC
/*--- Step 4 -------------------------------------------------------------------
 * Compute fluxes of energy
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      EFlux[ks][j][i].x = 
         0.25*(pG->U[ks][j][i].B2c + pG->U[ks][j][i-1].B2c)*
                (emf[ks][j][i].z + emf[ks][j+1][i].z)
       - 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j][i-1].B3c)*emf[ks][j][i].y;
    }
  }

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[ks][j][i].y =   
         0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c)*emf[ks][j][i].x
       - 0.25*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c)*
                (emf[ks][j][i].z + emf[ks][j][i+1].z);
    }
  }
#endif

/*--- Step 5 -------------------------------------------------------------------
 * CT update of magnetic field using Hall [+Ohmic] EMFs.  This is identical to
 * the CT update in the integrators: dB/dt = -Curl(E)
 */
  
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

/* Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].B1c = 0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c = 0.5*(pG->B2i[ks][j][i] + pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
      pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
    }
  }

#ifndef BAROTROPIC
/*--- Step 6 -------------------------------------------------------------------
 * Update energy using resistive fluxes in each dimension (dE/dt = Div(F))
 */
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].E  += dtodx1*(EFlux[ks][j][i+1].x - EFlux[ks][j][i].x);
      pG->U[ks][j][i].E  += dtodx2*(EFlux[ks][j+1][i].y - EFlux[ks][j][i].y);
    }
  }
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* hall_resistivity_3d:
 */

void hall_resistivity_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real Q_H;

/*--- Step 1 -------------------------------------------------------------------
 * Compute currents. Note:
 *   J1 = (dB3/dx2 - dB2/dx3)
 *   J2 = (dB1/dx3 - dB3/dx1)
 *   J3 = (dB2/dx1 - dB1/dx2)
 */

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
        J[k][j][i].x = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                        (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        J[k][j][i].y = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                        (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        J[k][j][i].z = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                        (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
      }
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute Hall EMFs
 *  emf.x =  Jy*Bz - Jz*By
 *  emf.y =  Jz*Bx - Jx*Bz
 *  emf.z =  Jx*By - Jy*Bx
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf[k][j][i].x = 
          0.25*((J[k  ][j  ][i].y + J[k  ][j  ][i+1].y)*pG->B3i[k  ][j  ][i] +
                (J[k  ][j-1][i].y + J[k  ][j-1][i+1].y)*pG->B3i[k  ][j-1][i])-
          0.25*((J[k  ][j  ][i].z + J[k  ][j  ][i+1].z)*pG->B2i[k  ][j  ][i] +
                (J[k-1][j  ][i].z + J[k-1][j  ][i+1].z)*pG->B2i[k-1][j  ][i]);

        emf[k][j][i].y = 
          0.25*((J[k  ][j][i  ].z + J[k  ][j+1][i  ].z)*pG->B1i[k  ][j][i  ] +
                (J[k-1][j][i  ].z + J[k-1][j+1][i  ].z)*pG->B1i[k-1][j][i  ])-
          0.25*((J[k  ][j][i  ].x + J[k  ][j+1][i  ].x)*pG->B3i[k  ][j][i  ] +
                (J[k  ][j][i-1].x + J[k  ][j+1][i-1].x)*pG->B3i[k  ][j][i-1]);

        emf[k][j][i].z = 
          0.25*((J[k][j  ][i  ].x + J[k+1][j  ][i  ].x)*pG->B2i[k][j  ][i  ] +
                (J[k][j  ][i-1].x + J[k+1][j  ][i-1].x)*pG->B2i[k][j  ][i-1])-
          0.25*((J[k][j  ][i  ].y + J[k+1][j  ][i  ].y)*pG->B1i[k][j  ][i  ] +
                (J[k][j-1][i  ].y + J[k+1][j-1][i  ].y)*pG->B1i[k][j-1][i  ]);

        Q_H = 4.0*eta_Hall/(pG->U[k  ][j][i].d + pG->U[k  ][j-1][i].d 
                          + pG->U[k-1][j][i].d + pG->U[k-1][j-1][i].d);
        emf[k][j][i].x *= Q_H;

        Q_H = 4.0*eta_Hall/(pG->U[k  ][j][i].d + pG->U[k  ][j][i-1].d 
                          + pG->U[k-1][j][i].d + pG->U[k-1][j][i-1].d);
        emf[k][j][i].y *= Q_H;

        Q_H = 4.0*eta_Hall/(pG->U[k][j  ][i].d + pG->U[k][j  ][i-1].d 
                          + pG->U[k][j-1][i].d + pG->U[k][j-1][i-1].d);
        emf[k][j][i].z *= Q_H;
      }
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Add Ohmic resistivity (\eta J) to Hall EMFs, if required.
 */

  if (eta_Ohm != 0.0) {
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie+1; i++) {
          emf[k][j][i].x += eta_Ohm*J[k][j][i].x;
          emf[k][j][i].y += eta_Ohm*J[k][j][i].y;
          emf[k][j][i].z += eta_Ohm*J[k][j][i].z;
        }
      }
    }
  }

#ifndef BAROTROPIC
/*--- Step 4 -------------------------------------------------------------------
 * Compute fluxes of energy
 */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      EFlux[k][j][i].x = 
         0.25*(pG->U[k][j][i].B2c + pG->U[k][j][i-1].B2c)*
                (emf[k][j][i].z + emf[k][j+1][i].z)
       - 0.25*(pG->U[k][j][i].B3c + pG->U[k][j][i-1].B3c)*
                (emf[k][j][i].y + emf[k+1][j][i].y);
    }
  }}

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[k][j][i].y =   
         0.25*(pG->U[k][j][i].B3c + pG->U[k][j-1][i].B3c)*
                (emf[k][j][i].x + emf[k+1][j][i].x)
       - 0.25*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c)*
                (emf[k][j][i].z + emf[k][j][i+1].z);
    }
  }}

  for (k=ks; k<=ke+1; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[k][j][i].z =   
         0.25*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c)*
                (emf[k][j][i].y + emf[k][j][i+1].y)
       - 0.25*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c)*
                (emf[k][j][i].x + emf[k][j+1][i].x);
    }
  }}
#endif

/*--- Step 5 -------------------------------------------------------------------
 * CT update of magnetic field using Hall [+Ohmic] EMFs.  This is identical to
 * the CT update in the integrators: dB/dt = -Curl(E)
 */

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
    }
  }

/* Set cell centered magnetic fields to average of updated face centered fields.
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
      }
    }
  }

#ifndef BAROTROPIC
/*--- Step 6 -------------------------------------------------------------------
 * Update energy using resistive fluxes in each dimension (dE/dt = Div(F))
 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E  += dtodx1*(EFlux[k][j][i+1].x - EFlux[k][j][i].x);
        pG->U[k][j][i].E  += dtodx2*(EFlux[k][j+1][i].y - EFlux[k][j][i].y);
        pG->U[k][j][i].E  += dtodx3*(EFlux[k+1][j][i].z - EFlux[k][j][i].z);
      }
    }
  }
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* hall_resistivity_init: Allocate temporary arrays
 */

void hall_resistivity_init(MeshS *pM)
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
  
  if ((J = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real3Vect)))
    == NULL) goto on_error;
  if ((emf = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real3Vect)))
    == NULL) goto on_error;
#ifndef BAROTROPIC
  if ((EFlux = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real3Vect)))
    == NULL) goto on_error;
#endif 
  return;

  on_error:
  hall_resistivity_destruct();
  ath_error("[hall_resisticvity_init]: malloc returned a NULL pointer\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* hall_resistivity_destruct: Free temporary arrays
 */

void hall_resistivity_destruct(void)
{
  if (J != NULL) free_3d_array(J);
  if (emf != NULL) free_3d_array(emf);
#ifndef BAROTROPIC
  if (EFlux != NULL) free_3d_array(EFlux);
#endif
  return;
}
#endif /* HALL_MHD */
