#include "../copyright.h"
/*==============================================================================
 * FILE: isotropic_conduction.c
 *
 * PURPOSE: Implements explicit isotropic thermal conduction, that is
 *      dE/dt = Div(kappa_T Grad(T))    where T=(P/d)*(mbar/k)=temperature
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   isoconduct() - isotropic conduction in 1D/2D/3D
 *   isoconduct_init() - allocates memory needed
 *   isoconduct_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef ISOTROPIC_CONDUCTION
#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif
#endif

#ifdef ISOTROPIC_CONDUCTION
/* Arrays for the temperature and heat fluxes */
static Real ***Temp=NULL;
static Real3Vect ***EFlux=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* isoconduct: isotropic conduction in 1d/2d/3d.
 */

void isoconduct(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, jl, ju, js = pG->js, je = pG->je;
  int k, kl, ku, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1, dtodx2=0.0, dtodx3=0.0;

  if (pG->Nx[1] > 1){
    jl = js - 1;
    ju = je + 1;
    dtodx2 = pG->dt/pG->dx2;
  } else {
    jl = js;
    ju = je;
  }
  if (pG->Nx[2] > 1){
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
        EFlux[k][j][i].x = (Temp[k][j][i] - Temp[k][j][i-1])/pG->dx1;
        EFlux[k][j][i].x *= kappa_T;
      }
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute heat fluxes in 2-direction
 */

  if (pG->Nx[1] > 1) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        EFlux[k][j][i].y = (Temp[k][j][i] - Temp[k][j-1][i])/pG->dx2;
        EFlux[k][j][i].y *= kappa_T;
      }
    }}
  }

/*--- Step 2c ------------------------------------------------------------------
 * Compute heat fluxes in 3-direction
 */

  if (pG->Nx[2] > 1) {
    for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        EFlux[k][j][i].z = (Temp[k][j][i] - Temp[k-1][j][i])/pG->dx3;
        EFlux[k][j][i].z *= kappa_T;
      }
    }}
  }

/*--- Step 3a ------------------------------------------------------------------
 * Update energy using x1-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx1*(EFlux[k][j][i+1].x - EFlux[k][j][i].x);
      }
    }
  }

/*--- Step 3b ------------------------------------------------------------------
 * Update energy using x2-fluxes
 */

  if (pG->Nx[1] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx2*(EFlux[k][j+1][i].y - EFlux[k][j][i].y);
      }
    }}
  }

/*--- Step 3c ------------------------------------------------------------------
 * Update energy using x3-fluxes
 */

  if (pG->Nx[2] > 1){
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx3*(EFlux[k+1][j][i].z - EFlux[k][j][i].z);
      }
    }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* isoconduct_init: Allocate temporary arrays
*/

void isoconduct_init(MeshS *pM)
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
  if ((Temp = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((EFlux = (Real3Vect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real3Vect)))
    == NULL) goto on_error;
  return;

  on_error:
  isoconduct_destruct();
  ath_error("[isoconduct_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* isoconduct_destruct: Free temporary arrays
 */

void isoconduct_destruct(void)
{
  if (Temp != NULL) free_3d_array(Temp);
  if (EFlux != NULL) free_3d_array(EFlux);
  return;
}
#endif /* ISOTROPIC_CONDUCTION */
