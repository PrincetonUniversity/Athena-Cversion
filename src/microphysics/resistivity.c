#include "../copyright.h"
/*==============================================================================
 * FILE: resistivity.c
 *
 * PURPOSE: Implements explicit resistivity using operator splitting, that is
 *      dB/dt = -Curl(\eta J)    where J=Curl(B).
 *   This function should be called after the integrate step in the main loop.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef RESISTIVITY
#ifdef ADIABATIC
#error : resistivity only works for isothermal EOS.
#endif /* ADIABATIC */
#endif


/* The resistive emfs, uses arrays allocated in integrator to save memory */
#ifdef MHD
extern Real ***emf1, ***emf2, ***emf3;
#endif /* MHD */


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* resistivity_3d:
 */

void resistivity_3d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real eta;

  eta = par_getd("problem","eta");

#ifdef RESISTIVITY

/*--- Step 1 -------------------------------------------------------------------
 * Compute resistive EMFs.  Note:
 *   emf1 = eta*J1 = eta*(dB3/dx2 - dB2/dx3)
 *   emf2 = eta*J2 = eta*(dB1/dx3 - dB3/dx1)
 *   emf3 = eta*J3 = eta*(dB2/dx1 - dB1/dx2)
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        emf1[k][j][i] = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                        (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;
        emf2[k][j][i] = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                        (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;
        emf3[k][j][i] = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                        (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;
/* Multiple components by constant \eta */
        emf1[k][j][i] *= eta;
        emf2[k][j][i] *= eta;
        emf3[k][j][i] *= eta;
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
        pG->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                            dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        pG->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                            dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        pG->B3i[k][j][i] += dtodx2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                            dtodx1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
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
      pG->B3i[ke+1][j][i] +=
        dtodx2*(emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
        dtodx1*(emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields. */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i]+pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* RESISTIVITY */

  return;
}
