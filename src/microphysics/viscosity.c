#include "../copyright.h"
/*==============================================================================
 * FILE: viscosity.c
 *
 * PURPOSE: Implements explicit viscosity using operator splitting, that is
 *      dM/dt = Div(T)    where T= viscous stress tensor.
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

#ifdef VISCOSITY

#ifdef ADIABATIC
#error : Navier-Stokes viscosity only works for isothermal EOS.
#endif /* ADIABATIC */

#ifdef BRAGINSKII
#error : Navier-Stokes viscosity cannot be used with Braginskii viscosity.
#endif /* BRAGINSKII */

#endif /* VISCOSITY */

/* The viscous fluxes, uses arrays allocated in integrator to save memory */
#ifdef VISCOSITY
extern Cons1D ***x1Flux, ***x2Flux, ***x3Flux;
Real ***Vx=NULL, ***Vy=NULL, ***Vz=NULL, ***divv=NULL;
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* viscous_3d:
 */

void viscosity_3d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real nu, nud;

  nu = par_getd("problem","nu");

#ifdef VISCOSITY
/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity and div(V) at cell centers (needed everywhere below)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Vx[k][j][i] = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        Vy[k][j][i] = pG->U[k][j][i].M2/pG->U[k][j][i].d;
        Vz[k][j][i] = pG->U[k][j][i].M3/pG->U[k][j][i].d;
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        divv[k][j][i] = ((Vx[k][j][i+1] - Vx[k][j][i-1])/(2.0*pG->dx1) +
                         (Vy[k][j+1][i] - Vy[k][j-1][i])/(2.0*pG->dx2) +
                         (Vz[k+1][j][i] - Vz[k-1][j][i])/(2.0*pG->dx3));
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        x1Flux[k][j][i].Mx = 2.0*(Vx[k][j][i] - Vx[k][j][i-1])/pG->dx1
           - ONE_3RD*(divv[k][j][i] + divv[k][j][i-1]);
        x1Flux[k][j][i].My = (Vy[k][j][i] - Vy[k][j][i-1])/pG->dx1
          + ((Vx[k][j+1][i]+Vx[k][j+1][i-1]) - (Vx[k][j-1][i]+Vx[k][j-1][i-1]))
            /(4.0*pG->dx2); 
        x1Flux[k][j][i].Mz = (Vz[k][j][i] - Vz[k][j][i-1])/pG->dx1
          + ((Vx[k+1][j][i]+Vx[k+1][j][i-1]) - (Vx[k-1][j][i]+Vx[k-1][j][i-1]))
            /(4.0*pG->dx3);
#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu*0.5*(pG->U[k][j][i].d + pG->U[k][j][i-1].d);
        x1Flux[k][j][i].Mx *= nud;
        x1Flux[k][j][i].My *= nud;
        x1Flux[k][j][i].Mz *= nud;
      }
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute viscous fluxes in 2-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        x2Flux[k][j][i].Mx = (Vx[k][j][i] - Vx[k][j-1][i])/pG->dx2
          + ((Vy[k][j][i+1]+Vy[k][j-1][i+1]) - (Vy[k][j][i-1]+Vy[k][j-1][i-1]))
            /(4.0*pG->dx1);
        x2Flux[k][j][i].My = 2.0*(Vy[k][j][i] - Vy[k][j-1][i])/pG->dx2
           - ONE_3RD*(divv[k][j][i] + divv[k][j-1][i]);
        x2Flux[k][j][i].Mz = (Vz[k][j][i] - Vz[k][j-1][i])/pG->dx2
          + ((Vy[k+1][j][i]+Vy[k+1][j-1][i]) - (Vy[k-1][j][i]+Vy[k-1][j-1][i]))
             /(4.0*pG->dx3);
#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu*0.5*(pG->U[k][j][i].d + pG->U[k][j-1][i].d);
        x2Flux[k][j][i].Mx *= nud;
        x2Flux[k][j][i].My *= nud;
        x2Flux[k][j][i].Mz *= nud;
      }
    }
  }

/*--- Step 2c ------------------------------------------------------------------
 * Compute viscous fluxes in 3-direction
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        x3Flux[k][j][i].Mx = (Vx[k][j][i] - Vx[k-1][j][i])/pG->dx3
          + ((Vz[k][j][i+1]+Vz[k-1][j][i+1]) - (Vz[k][j][i-1]+Vz[k-1][j][i-1]))
            /(4.0*pG->dx1);
        x3Flux[k][j][i].My = (Vy[k][j][i] - Vy[k-1][j][i])/pG->dx3
          + ((Vz[k][j+1][i]+Vz[k-1][j+1][i]) - (Vz[k][j-1][i]+Vz[k-1][j-1][i]))
            /(4.0*pG->dx2);
        x3Flux[k][j][i].Mz = 2.0*(Vz[k][j][i] - Vz[k-1][j][i])/pG->dx3
           - ONE_3RD*(divv[k][j][i] + divv[k-1][j][i]);
#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu*0.5*(pG->U[k][j][i].d + pG->U[k-1][j][i].d);
        x3Flux[k][j][i].Mx *= nud;
        x3Flux[k][j][i].My *= nud;
        x3Flux[k][j][i].Mz *= nud;
      }
    }
  }

/*--- Step 3a ------------------------------------------------------------------
 * Update momentum and energy using x1-fluxes (dM/dt = Div(T))
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].M1 += dtodx1*(x1Flux[k][j][i+1].Mx-x1Flux[k][j][i].Mx);
        pG->U[k][j][i].M2 += dtodx1*(x1Flux[k][j][i+1].My-x1Flux[k][j][i].My);
        pG->U[k][j][i].M3 += dtodx1*(x1Flux[k][j][i+1].Mz-x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  += dtodx1*(x1Flux[k][j][i+1].E -x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
      }
    }
  }

/*--- Step 3b ------------------------------------------------------------------
 * Update momentum and energy using x2-fluxes (dM/dt = Div(T))
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].M1 += dtodx2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        pG->U[k][j][i].M2 += dtodx2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
        pG->U[k][j][i].M3 += dtodx2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  +=dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
      }
    }
  }

/*--- Step 3c ------------------------------------------------------------------
 * Update momentum and energy using x3-fluxes (dM/dt = Div(T))
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].M1 += dtodx3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
        pG->U[k][j][i].M2 += dtodx3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        pG->U[k][j][i].M3 += dtodx3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  += dtodx3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
      }
    }
  }

#endif /* VISCOSITY */
  return;
}

/*----------------------------------------------------------------------------*/
/* viscous_init_3d: Allocate temporary integration arrays
*/

void viscosity_init_3d(int nx1, int nx2, int nx3)
{
  int Nx1 = nx1 + 2*nghost;
  int Nx2 = nx2 + 2*nghost;
  int Nx3 = nx3 + 2*nghost;

#ifdef VISCOSITY
  if ((Vx = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Vy = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Vz = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((divv = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* VISCOSITY */

  return;

#ifdef VISCOSITY
  on_error:
  integrate_destruct();
  ath_error("[viscosity_init]: malloc returned a NULL pointer\n");
#endif /* VISCOSITY */
}
