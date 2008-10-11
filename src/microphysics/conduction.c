#include "../copyright.h"
/*==============================================================================
 * FILE: conduction.c
 *
 * PURPOSE: Implements explicit thermal conduction using operator splitting,
 *      dE/dt = Div(kappa_T Grad(T))    where T=(P/d)*(mbar/k)=temperature
 *
 *   Functions for both isotropic and anisotropic conduction are included.
 *   These are selected in configure using 
 *      --with-conduction=[isotropic,anisotropic]
 *
 *   Isotropic conduction works in 1D, 2D and 3D.
 *   Anisotropic conduction works in 2D and 3D.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   isoconduct() - isotropic conduction in 1D/2D/3D
 *   anisoconduct_2d() - anisotropic conduction in 2D
 *   anisoconduct_3d() - anisotropic conduction in 3D
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef CONDUCTION
#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif
#endif

/* Arrays for the temperature and heat fluxes */
static Real ***Temp=NULL, ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* viscous_3d:
 */

void isoconduct(Grid *pG, Domain *pD)
{
#ifndef BAROTROPIC
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1, dtodx2=0.0, dtodx3=0.0;

  if (pG->Nx2 > 1){
    jl = js - 1;
    ju = je + 1;
    dtodx2 = pG->dt/pG->dx2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx3 > 1){
    kl = ks - 1;
    ku = ke + 1;
    dtodx3 = pG->dt/pG->dx3;
  } else {
    kl = ks;
    ku = ke;
  }

/*--- Step 1 -------------------------------------------------------------------
 * Compute temperature at cell centers. Note temperature is dimensionless
 * (in units of [mbar/k]).  For cgs units, the value of kappa_T passed to this
 * routine must be multiplied by (mbar/k).
 */

  for (k=kl; k<=ku; k++) {
  for (j=jl; j<=ju; j++) {
  for (i=is-1; i<=ie+1; i++) {
    Temp[k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
      (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));
#ifdef MHD
    Temp[k][j][i] -= (0.5)*
      (SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
#endif
    Temp[k][j][i] *= (Gamma_1/pG->U[k][j][i].d);
  }}}

/*--- Step 2a ------------------------------------------------------------------
 * Compute heat fluxes in 1-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        x1Flux[k][j][i] = (Temp[k][j][i] - Temp[k][j][i-1])/pG->dx1;
        x1Flux[k][j][i] *= kappa_T;
      }
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute heat fluxes in 2-direction
 */

  if (pG->Nx2 > 1) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        x2Flux[k][j][i] = (Temp[k][j][i] - Temp[k][j-1][i])/pG->dx2;
        x2Flux[k][j][i] *= kappa_T;
      }
    }}
  }

/*--- Step 2c ------------------------------------------------------------------
 * Compute heat fluxes in 3-direction
 */

  if (pG->Nx3 > 1) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        x3Flux[k][j][i] = (Temp[k][j][i] - Temp[k-1][j][i])/pG->dx3;
        x3Flux[k][j][i] *= kappa_T;
      }
    }}
  }

/*--- Step 3a ------------------------------------------------------------------
 * Update energy using x1-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx1*(x1Flux[k][j][i+1] - x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 3b ------------------------------------------------------------------
 * Update energy using x2-fluxes
 */

  if (pG->Nx2 > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx2*(x2Flux[k][j+1][i] - x2Flux[k][j][i]);
      }
    }}
  }

/*--- Step 3c ------------------------------------------------------------------
 * Update energy using x3-fluxes
 */

  if (pG->Nx3 > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(x3Flux[k+1][j][i] - x3Flux[k][j][i]);
      }
    }}
  }
#endif /* BAROTROPIC */

  return;
}


/*----------------------------------------------------------------------------*/
/* isoconduct_init_3d: Allocate temporary integration arrays
*/

void isoconduct_init(int nx1, int nx2, int nx3)
{
  int Nx1 = nx1 + 2;
  int Nx2 = nx2 + 2;
  int Nx3 = nx3 + 2;

  if ((Temp = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x1Flux = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x2Flux = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((x3Flux = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;

  return;

  on_error:
  integrate_destruct();
  ath_error("[isoconduct_init]: malloc returned a NULL pointer\n");
}
