#include "../copyright.h"
/*==============================================================================
 * FILE: braginskii.c
 *
 * PURPOSE: Implements anisotropic (Braginskii) viscosity using operator
 *   splitting, that is
 *      dM/dt = Div(T)    where T=Braginski stress tensor.
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

#ifdef BRAGINSKII

#ifdef ADIABATIC
#error : Braginskii viscosity only works for isothermal EOS.
#endif /* ADIABATIC */

#ifdef HYDRO
#error : Braginskii viscosity only works for MHD.
#endif /* HYDRO */

#ifdef VISCOSITY
#error : Braginskii viscosity cannot be used with Navier-Stokes viscosity
#endif /* VISCOSITY */

#endif /* BRAGINSKII */

/* The viscous fluxes, uses arrays allocated in integrator to save memory */
#ifdef BRAGINSKII
extern Cons1D ***x1Flux, ***x2Flux, ***x3Flux;
Real ***Vx=NULL, ***Vy=NULL, ***Vz=NULL, ***divv=NULL, ***BBdV=NULL;
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* viscous_3d:
 */

void braginskii_3d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real B02,Bx,By,Bz,qa,nu,nud;
#ifdef FARGO
  Real x1,x2,x3;
#endif

  nu = par_getd("problem","nu");


#ifdef BRAGINSKII
#ifdef MHD
/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity, div(V), and BBdV=B_{m}B_{k}d_{k}B_{m} at cell centers
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Vx[k][j][i] = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        Vy[k][j][i] = pG->U[k][j][i].M2/pG->U[k][j][i].d;
#ifdef FARGO
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        Vy[k][j][i] -= 1.5*Omega*x1;
#endif
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
        B02 = pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
              pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
              pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;
        B02 = MAX(B02,TINY_NUMBER);
        BBdV[k][j][i] = (pG->U[k][j][i].B1c/B02)*
           (pG->U[k][j][i].B1c*(Vx[k][j][i+1]-Vx[k][j][i-1])/(2.0*pG->dx1)
          + pG->U[k][j][i].B2c*(Vx[k][j+1][i]-Vx[k][j-1][i])/(2.0*pG->dx2)
          + pG->U[k][j][i].B3c*(Vx[k+1][j][i]-Vx[k-1][j][i])/(2.0*pG->dx3))
                      + (pG->U[k][j][i].B2c/B02)*
           (pG->U[k][j][i].B1c*(Vy[k][j][i+1]-Vy[k][j][i-1])/(2.0*pG->dx1)
          + pG->U[k][j][i].B2c*(Vy[k][j+1][i]-Vy[k][j-1][i])/(2.0*pG->dx2)
          + pG->U[k][j][i].B3c*(Vy[k+1][j][i]-Vy[k-1][j][i])/(2.0*pG->dx3))
                      + (pG->U[k][j][i].B3c/B02)*
           (pG->U[k][j][i].B1c*(Vz[k][j][i+1]-Vz[k][j][i-1])/(2.0*pG->dx1)
          + pG->U[k][j][i].B2c*(Vz[k][j+1][i]-Vz[k][j-1][i])/(2.0*pG->dx2)
          + pG->U[k][j][i].B3c*(Vz[k+1][j][i]-Vz[k-1][j][i])/(2.0*pG->dx3));
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        qa = 0.5*((BBdV[k][j][i  ] - ONE_3RD*divv[k][j][i  ])
                + (BBdV[k][j][i-1] - ONE_3RD*divv[k][j][i-1]));
        Bx = pG->B1i[k][j][i];
        By = 0.5*(pG->U[k][j][i].B2c + pG->U[k][j][i-1].B2c);
        Bz = 0.5*(pG->U[k][j][i].B3c + pG->U[k][j][i-1].B3c);
        B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);
        x1Flux[k][j][i].Mx = qa*(3.0*Bx*Bx/B02 - 1.0);
        x1Flux[k][j][i].My = qa*(3.0*By*Bx/B02);
        x1Flux[k][j][i].Mz = qa*(3.0*Bz*Bx/B02);
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
        qa = 0.5*((BBdV[k][j  ][i] - ONE_3RD*divv[k][j  ][i])
                + (BBdV[k][j-1][i] - ONE_3RD*divv[k][j-1][i]));
        Bx = 0.5*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c);
        By = pG->B2i[k][j][i];
        Bz = 0.5*(pG->U[k][j][i].B3c + pG->U[k][j-1][i].B3c);
        B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);
        x2Flux[k][j][i].Mx = qa*(3.0*Bx*By/B02);
        x2Flux[k][j][i].My = qa*(3.0*By*By/B02 - 1.0);
        x2Flux[k][j][i].Mz = qa*(3.0*Bz*By/B02);
#ifndef BAROTROPIC
        x2Flux[k][j][i].E  =
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
        qa = 0.5*((BBdV[k  ][j][i] - ONE_3RD*divv[k  ][j][i])
                + (BBdV[k-1][j][i] - ONE_3RD*divv[k-1][j][i]));
        Bx = 0.5*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c);
        By = 0.5*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c);
        Bz = pG->B3i[k][j][i];
        B02 = Bx*Bx + By*By + Bz*Bz;
        B02 = MAX(B02,TINY_NUMBER);
        x3Flux[k][j][i].Mx = qa*(3.0*Bx*Bz/B02);
        x3Flux[k][j][i].My = qa*(3.0*By*Bz/B02);
        x3Flux[k][j][i].Mz = qa*(3.0*Bz*Bz/B02 - 1.0);
#ifndef BAROTROPIC
        x3Flux[k][j][i].E  =
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

#endif /* MHD */
#endif /* BRAGINSKII */
  return;
}

/*----------------------------------------------------------------------------*/
/* viscous_init_3d: Allocate temporary integration arrays
*/

void braginskii_init_3d(int nx1, int nx2, int nx3)
{
  int Nx1 = nx1 + 2*nghost;
  int Nx2 = nx2 + 2*nghost;
  int Nx3 = nx3 + 2*nghost;

#ifdef BRAGINSKII
  if ((Vx = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Vy = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((Vz = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((divv = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((BBdV = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* BRAGINSKII */

  return;

#ifdef BRAGINSKII
  on_error:
  integrate_destruct();
  ath_error("[braginskii_init]: malloc returned a NULL pointer\n");
#endif /* BRAGINSKII */
}
