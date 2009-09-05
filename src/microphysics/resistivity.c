#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Implements explicit Ohmic resistivity, that is
 *      dB/dt = -Curl(\eta J)    where J=Curl(B).
 *      dE/dt = Div(B X \eta J)
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  ohmic_resistivity_1d()
 *  ohmic_resistivity_2d()
 *  ohmic_resistivity_3d()
 *  ohmic_resistivity_init() - allocates memory needed
 *  ohmic_resistivity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef OHMIC
#ifdef HYDRO
#error : Ohmic resistivity only works for MHD.
#endif /* HYDRO */
#endif

#ifdef OHMIC
/* The resistive emfs and energy flux, contained in 3D vector structure */
static Vector ***emf=NULL, ***EFlux=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_1d:
 */

void ohmic_resistivity_1d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMF.  Note:
 *   emf.x1 = eta*J1 = 0
 *   emf.x2 = eta*J2 = eta_Ohm*(-dB3/dx1)
 *   emf.x3 = eta*J3 = eta_Ohm*(dB2/dx1)
 * emf.x2 and emf.x3 use B3c and B2c respectively, and are centered at x1-faces
 */

  for (i=is; i<=ie+1; i++) {
    emf[ks][js][i].x1 = 0.0;
    emf[ks][js][i].x2 = -(pG->U[ks][js][i].B3c -pG->U[ks][js][i-1].B3c)/pG->dx1;
    emf[ks][js][i].x3 =  (pG->U[ks][js][i].B2c -pG->U[ks][js][i-1].B2c)/pG->dx1;
/* Multiple components by constant \eta_Ohm */
    emf[ks][js][i].x2 *= eta_Ohm;
    emf[ks][js][i].x3 *= eta_Ohm;
  }

#ifndef BAROTROPIC
/*--- Step 2 -------------------------------------------------------------------
 * Compute flux of total energy due to resistive diffusion = B X emf
 *  EFlux.x1 =  By*emf.x3 - Bz*emf.y
 */   
      
  for (i=is; i<=ie+1; i++) {
    EFlux[ks][js][i].x1 =
       0.5*(pG->U[ks][js][i].B2c + pG->U[ks][js][i-1].B2c)*emf[ks][js][i].x3
     - 0.5*(pG->U[ks][js][i].B3c + pG->U[ks][js][i-1].B3c)*emf[ks][js][i].x2;
  }
#endif

/*--- Step 3 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMF.  In 1D, this reduces to
 * centered differences for the resistive fluxes of B2c and B3c
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].B2c += dtodx1*(emf[ks][js][i+1].x3 - emf[ks][js][i].x3);
    pG->U[ks][js][i].B3c -= dtodx1*(emf[ks][js][i+1].x2 - emf[ks][js][i].x2);
  }

/* For consistency, set B2i and B3i to cell-centered values. */

  for (i=is; i<=ie; i++) {
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
  }

#ifndef BAROTROPIC
/*--- Step 4 -------------------------------------------------------------------
 * Update energy using resistive fluxes (dE/dt = Div(F))
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].E  += dtodx1*(EFlux[ks][js][i+1].x1 - EFlux[ks][js][i].x1);
  }                       
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_2d:
 */

void ohmic_resistivity_2d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMF.  Note:
 *   emf.x1 = eta*J1 = eta_Ohm*(dB3/dx2)
 *   emf.x2 = eta*J2 = eta_Ohm*(-dB3/dx1)
 *   emf.x3 = eta*J3 = eta_Ohm*(dB2/dx1 - dB1/dx2)
 * emf.x1 and emf.x2 use B3c, and in 2D are centered at x2- and x1-interfaces
 */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      emf[ks][j][i].x1 =  (pG->U[ks][j][i].B3c - pG->U[ks][j-1][i].B3c)/pG->dx2;
      emf[ks][j][i].x2 = -(pG->U[ks][j][i].B3c - pG->U[ks][j][i-1].B3c)/pG->dx1;

      emf[ks][j][i].x3 = (pG->B2i[ks][j][i] - pG->B2i[ks][j  ][i-1])/pG->dx1 -
                         (pG->B1i[ks][j][i] - pG->B1i[ks][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta_Ohm */
      emf[ks][j][i].x1 *= eta_Ohm;
      emf[ks][j][i].x2 *= eta_Ohm;
      emf[ks][j][i].x3 *= eta_Ohm;
    }
  }

#ifndef BAROTROPIC
/*--- Step 2 -------------------------------------------------------------------
 * Compute flux of total energy due to resistive diffusion = B X emf
 *  EFlux.x1 =  By*emf.x3 - Bz*emf.x2
 *  EFlux.x2 =  Bz*emf.x1 - Bx*emf.x3
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      EFlux[ks][j][i].x1 =
         0.25*(pG->U[ks][j][i].B2c + pG->U[ks][j][i-1].B2c)*
                (emf[ks][j][i].x3 + emf[ks][j+1][i].x3)
       - 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j][i-1].B3c)*emf[ks][j][i].x2;
    }
  }
  
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[ks][j][i].x2 =   
         0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c)*emf[ks][j][i].x1
       - 0.25*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c)*
                (emf[ks][j][i].x3 + emf[ks][j][i+1].x3);
    }
  }
#endif

/*--- Step 3 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMF.  This is identical to the
 * CT update in the 2D integrators: dB/dt = -Curl(E).  For B3, the CT formula
 * reduces to centered differences for the diffusive (resistive) fluxes of B3c.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i] -= dtodx2*(emf[ks][j+1][i  ].x3 - emf[ks][j][i].x3);
      pG->B2i[ks][j][i] += dtodx1*(emf[ks][j  ][i+1].x3 - emf[ks][j][i].x3);

      pG->U[ks][j][i].B3c += dtodx2*(emf[ks][j+1][i  ].x1 - emf[ks][j][i].x1) -
                             dtodx1*(emf[ks][j  ][i+1].x2 - emf[ks][j][i].x2);
    }
    pG->B1i[ks][j][ie+1] -= dtodx2*(emf[ks][j+1][ie+1].x3 -emf[ks][j][ie+1].x3);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf[ks][je+1][i+1].x3 -emf[ks][je+1][i].x3);
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
/*--- Step 4 -------------------------------------------------------------------
 * Update energy using resistive fluxes in each dimension (dE/dt = Div(F))
 */
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].E  += dtodx1*(EFlux[ks][j][i+1].x1 - EFlux[ks][j][i].x1);
      pG->U[ks][j][i].E  += dtodx2*(EFlux[ks][j+1][i].x2 - EFlux[ks][j][i].x2);
    } 
  }                       
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_3d:
 */

void ohmic_resistivity_3d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMFs.  Note:
 *   emf.x1 = eta*J1 = eta_Ohm*(dB3/dx2 - dB2/dx3)
 *   emf.x2 = eta*J2 = eta_Ohm*(dB1/dx3 - dB3/dx1)
 *   emf.x3 = eta*J3 = eta_Ohm*(dB2/dx1 - dB1/dx2)
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf[k][j][i].x1 = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                          (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        emf[k][j][i].x2 = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                          (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        emf[k][j][i].x3 = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                          (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta_Ohm */
        emf[k][j][i].x1 *= eta_Ohm;
        emf[k][j][i].x2 *= eta_Ohm;
        emf[k][j][i].x3 *= eta_Ohm;
      }
    }
  }

#ifndef BAROTROPIC
/*--- Step 2 -------------------------------------------------------------------
 * THIS WILL USE OLD B
 * Compute flux of total energy due to resistive diffusion = B X emf
 *  EFlux.x1 =  By*emf.x3 - Bz*emf.x2
 *  EFlux.x2 =  Bz*emf.x1 - Bx*emf.x3
 *  EFlux.x3 =  Bx*emf.x2 - By*emf.x1
 */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      EFlux[k][j][i].x1 = 
         0.25*(pG->U[k][j][i].B2c + pG->U[k][j][i-1].B2c)*
                (emf[k][j][i].x3 + emf[k][j+1][i].x3)
       - 0.25*(pG->U[k][j][i].B3c + pG->U[k][j][i-1].B3c)*
                (emf[k][j][i].x2 + emf[k+1][j][i].x2);
    }
  }}

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[k][j][i].x2 =   
         0.25*(pG->U[k][j][i].B3c + pG->U[k][j-1][i].B3c)*
                (emf[k][j][i].x1 + emf[k+1][j][i].x1)
       - 0.25*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c)*
                (emf[k][j][i].x3 + emf[k][j][i+1].x3);
    }
  }}

  for (k=ks; k<=ke+1; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      EFlux[k][j][i].x3 =   
         0.25*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c)*
                (emf[k][j][i].x2 + emf[k][j][i+1].x2)
       - 0.25*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c)*
                (emf[k][j][i].x1 + emf[k][j+1][i].x1);
    }
  }}
#endif

/*--- Step 3 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMFs.  This is identical to the
 * CT update in the integrators: dB/dt = -Curl(E)
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->B1i[k][j][i] += dtodx3*(emf[k+1][j  ][i  ].x2 - emf[k][j][i].x2) -
                            dtodx2*(emf[k  ][j+1][i  ].x3 - emf[k][j][i].x3);
        pG->B2i[k][j][i] += dtodx1*(emf[k  ][j  ][i+1].x3 - emf[k][j][i].x3) -
                            dtodx3*(emf[k+1][j  ][i  ].x1 - emf[k][j][i].x1);
        pG->B3i[k][j][i] += dtodx2*(emf[k  ][j+1][i  ].x1 - emf[k][j][i].x1) -
                            dtodx1*(emf[k  ][j  ][i+1].x2 - emf[k][j][i].x2);
      }
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
      pG->B3i[ke+1][j][i] +=
        dtodx2*(emf[ke+1][j+1][i  ].x1 - emf[ke+1][j][i].x1) -
        dtodx1*(emf[ke+1][j  ][i+1].x2 - emf[ke+1][j][i].x2);
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
/*--- Step 4 -------------------------------------------------------------------
 * Update energy using resistive fluxes in each dimension (dE/dt = Div(F))
 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E  += dtodx1*(EFlux[k][j][i+1].x1 - EFlux[k][j][i].x1);
        pG->U[k][j][i].E  += dtodx2*(EFlux[k][j+1][i].x2 - EFlux[k][j][i].x2);
        pG->U[k][j][i].E  += dtodx3*(EFlux[k+1][j][i].x3 - EFlux[k][j][i].x3);
      }
    }
  }
#endif /* BAROTROPIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_init: Allocate temporary arrays
 */

void ohmic_resistivity_init(int nx1, int nx2, int nx3)
{
  int Nx1 = nx1 + 2*nghost, Nx2, Nx3;
  if (nx2 > 1){
    Nx2 = nx2 + 2*nghost;
  } else {
    Nx2 = nx2;
  }
  if (nx3 > 1){
    Nx3 = nx3 + 2*nghost;
  } else {
    Nx3 = nx3;
  }
  
  if ((emf = (Vector***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Vector)))
    == NULL) goto on_error;
#ifndef BAROTROPIC
  if ((EFlux = (Vector***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Vector)))
    == NULL) goto on_error;
#endif 
  return;

  on_error:
  ohmic_resistivity_destruct();
  ath_error("[ohmic_resisticvity_init]: malloc returned a NULL pointer\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_destruct: Free temporary arrays
 */

void ohmic_resistivity_destruct(void)
{
  if (emf != NULL) free_3d_array(emf);
#ifndef BAROTROPIC
  if (EFlux != NULL) free_3d_array(EFlux);
#endif
  return;
}
#endif /* OHMIC */
