#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Implements explicit Ohmic resistivity, that is
 *      dB/dt = -Curl(\eta J)    where J=Curl(B).
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
#ifdef ADIABATIC
#error : resistivity only works for isothermal EOS.
#endif /* ADIABATIC */
#ifdef HYDRO
#error : Ohmic resistivity only works for MHD.
#endif /* HYDRO */
#endif

/* The resistive emfs, contained in special structure */
typedef struct ThreeDVect_t{
  Real x;
  Real y;
  Real z;
}ThreeDVect;
static ThreeDVect ***emf=NULL;

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_1d:
 */

void ohmic_resistivity_1d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_2d:
 */

void ohmic_resistivity_2d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_3d:
 */

void ohmic_resistivity_3d(Grid *pG, Domain *pD)
{
#ifdef OHMIC
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMFs.  Note:
 *   emf.x = eta*J1 = eta_R*(dB3/dx2 - dB2/dx3)
 *   emf.y = eta*J2 = eta_R*(dB1/dx3 - dB3/dx1)
 *   emf.z = eta*J3 = eta_R*(dB2/dx1 - dB1/dx2)
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf[k][j][i].x = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                         (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        emf[k][j][i].y = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                         (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        emf[k][j][i].z = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                         (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta_R */
        emf[k][j][i].x *= eta_R;
        emf[k][j][i].y *= eta_R;
        emf[k][j][i].z *= eta_R;
      }
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * CT update of magnetic field using resistive EMFs.  This is identical to the
 * CT update in the integrators: dB/dt = -Curl(E)
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

/*--- Step 3 -------------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields. */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* OHMIC */

  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_init: Allocate temporary arrays
 */

void ohmic_resistivity_init(int nx1, int nx2, int nx3)
{
#ifdef OHMIC
  int Nx1 = nx1 + 2;
  int Nx2 = nx2 + 2;
  int Nx3 = nx3 + 2;
/* Calculate the dimensions  */
  dim=0;
  if(Nx1 > 1) dim++;
  if(Nx2 > 1) dim++;
  if(Nx3 > 1) dim++;
  
  switch(dim){
  case 1:
    break;
  case 2:
    break;
  case 3:
    if ((emf = (ThreeDVect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ThreeDVect)))
      == NULL) goto on_error;
    break;
  }
  return;

  on_error:
  ohmic_resistivity_destruct();
  ath_error("[ohmic_resisticvity_init]: malloc returned a NULL pointer\n");
#endif /* OHMIC */
  return;
}

/*----------------------------------------------------------------------------*/
/* ohmic_resistivity_destruct: Free temporary arrays
 */

void ohmic_resistivity_destruct(void)
{
#ifdef OHMIC
/* dim set in ohmic_resistivity_init() at begnning of run */
  switch(dim){
  case 1:
    break;
  case 2:
    break;
  case 3:
    if (emf != NULL) free_3d_array(emf);
    break;
  }
#endif /* OHMIC */

  return;
}
