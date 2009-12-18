#include "../copyright.h"
/*==============================================================================
 * FILE: anisotropic_conduction.c
 *
 * PURPOSE: Implements explicit anisotropic thermal conduction, that is
 *      dE/dt = Div(chi_C[b Dot Grad(T)]b)  where T=(P/d)*(mbar/k)=temperature
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 *   Based on the original version written by Ian Parrish 10/11/2008
 *   Reference: Parrish & Stone (2005), ApJ, 633, 334
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   anisoconduct_2d() - anisotropic conduction in 2D
 *   anisoconduct_3d() - anisotropic conduction in 3D
 *   anisoconduct_init() - allocates memory needed
 *   anisoconduct_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef ANISOTROPIC_CONDUCTION
#ifdef BAROTROPIC
#error : Thermal conduction requires an adiabatic EOS
#endif /* BAROTROPIC */
#ifdef HYDRO
#error : Anisotropic conduction requires MHD
#endif /* HYDRO */
#endif /* ANISOTROPIC_CONDUCTION */

/* Arrays for the temperature and heat fluxes */
#ifdef ANISOTROPIC_CONDUCTION
static Real ***Temp=NULL;
static Real3Vect ***EFlux=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* anisoconduct_2d: anisotropic conduction in 2d.
 */

void anisoconduct_2d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real Bx,By,B02,dTc,dTl,dTr,lim_slope,dTdx,dTdy,bDotGradT;

/*--- Step 1 -------------------------------------------------------------------
 * Compute temperature at cell centers. Note temperature is dimensionless
 * (in units of [mbar/k]).  For cgs units, the value of chi_C passed to this
 * routine must be multiplied by (mbar/k).
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      Temp[ks][j][i] = pG->U[ks][j][i].E - (0.5/pG->U[ks][j][i].d)*
      (SQR(pG->U[ks][j][i].M1)+SQR(pG->U[ks][j][i].M2)+SQR(pG->U[ks][j][i].M3));
#ifdef MHD
      Temp[ks][j][i] -= (0.5)*(SQR(pG->U[ks][j][i].B1c) +
         SQR(pG->U[ks][j][i].B2c) + SQR(pG->U[ks][j][i].B3c));
#endif
      Temp[ks][j][i] *= (Gamma_1/pG->U[ks][j][i].d);
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute heat fluxes in 1-direction
 */
  
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
/* Compute normalized field at x1-interface */
      By = 0.5*(pG->U[ks][j][i-1].B2c + pG->U[ks][j][i].B2c);
      B02 = SQR(pG->B1i[ks][j][i]) + SQR(By);
      B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized temperature difference dT/dy */
      dTr = 0.5*((Temp[ks][j+1][i-1] + Temp[ks][j+1][i]) -
                 (Temp[ks][j  ][i-1] + Temp[ks][j  ][i]));
      dTl = 0.5*((Temp[ks][j  ][i-1] + Temp[ks][j  ][i]) -
                 (Temp[ks][j-1][i-1] + Temp[ks][j-1][i]));
      dTc = dTr + dTl;

      dTdy = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
      }

      bDotGradT = pG->B1i[ks][j][i]*(Temp[ks][j][i]-Temp[ks][j][i-1])/pG->dx1
         + By*dTdy;
      EFlux[ks][j][i].x = pG->B1i[ks][j][i]*bDotGradT/B02;
      EFlux[ks][j][i].x *= chi_C;
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute heat fluxes in 2-direction
 */
  
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
/* Compute normalized field at x2-interface */
      Bx = 0.5*(pG->U[ks][j-1][i].B1c + pG->U[ks][j][i].B1c);
      B02 = SQR(Bx) + SQR(pG->B2i[ks][j][i]);
      B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized temperature difference dT/dx */
      dTr = 0.5*((Temp[ks][j-1][i+1] + Temp[ks][j][i+1]) -
                 (Temp[ks][j-1][i  ] + Temp[ks][j][i  ]));
      dTl = 0.5*((Temp[ks][j-1][i  ] + Temp[ks][j][i  ]) -
                 (Temp[ks][j-1][i-1] + Temp[ks][j][i-1]));
      dTc = dTr + dTl;

      dTdx = 0.0;
      if (dTl*dTr > 0.0) {
        lim_slope = MIN(fabs(dTl),fabs(dTr));
        dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
      }

      bDotGradT = pG->B2i[ks][j][i]*(Temp[ks][j][i]-Temp[ks][j-1][i])/pG->dx2
         + Bx*dTdx;
      EFlux[ks][j][i].y = pG->B2i[ks][j][i]*bDotGradT/B02;
      EFlux[ks][j][i].y *= chi_C;
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Update energy using heat-fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].E += dtodx1*(EFlux[ks][j][i+1].x - EFlux[ks][j][i].x);
      pG->U[ks][j][i].E += dtodx2*(EFlux[ks][j+1][i].y - EFlux[ks][j][i].y);
    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/* anisoconduct_3d: anisotropic conduction in 3d.
 */

void anisoconduct_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real Bx,By,Bz,B02,dTc,dTl,dTr,lim_slope,dTdx,dTdy,dTdz,bDotGradT;

/*--- Step 1 -------------------------------------------------------------------
 * Compute temperature at cell centers. Note temperature is dimensionless
 * (in units of [mbar/k]).  For cgs units, the value of chi_C passed to this
 * routine must be multiplied by (mbar/k).
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        Temp[k][j][i] = pG->U[k][j][i].E - (0.5/pG->U[k][j][i].d)*
         (SQR(pG->U[k][j][i].M1)+SQR(pG->U[k][j][i].M2)+SQR(pG->U[k][j][i].M3));
#ifdef MHD
        Temp[k][j][i] -= (0.5)*(SQR(pG->U[k][j][i].B1c) +
           SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
        Temp[k][j][i] *= (Gamma_1/pG->U[k][j][i].d);
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute heat fluxes in 1-direction
 */
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
/* Compute normalized field at x1-interface */
        By = 0.5*(pG->U[k][j][i-1].B2c + pG->U[k][j][i].B2c);
        Bz = 0.5*(pG->U[k][j][i-1].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(pG->B1i[k][j][i]) + SQR(By) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized temperature difference dT/dy */
        dTr = 0.5*((Temp[k][j+1][i-1] + Temp[k][j+1][i]) -
                   (Temp[k][j  ][i-1] + Temp[k][j  ][i]));
        dTl = 0.5*((Temp[k][j  ][i-1] + Temp[k][j  ][i]) -
                   (Temp[k][j-1][i-1] + Temp[k][j-1][i]));
        dTc = dTr + dTl;

        dTdy = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
        }

/* Monotonized temperature difference dT/dz */
        dTr = 0.5*((Temp[k+1][j][i-1] + Temp[k+1][j][i]) -
                   (Temp[k  ][j][i-1] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j][i-1] + Temp[k  ][j][i]) -
                   (Temp[k-1][j][i-1] + Temp[k-1][j][i]));
        dTc = dTr + dTl;

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }

        bDotGradT = pG->B1i[k][j][i]*(Temp[k][j][i]-Temp[k][j][i-1])/pG->dx1
           + By*dTdy + Bz*dTdz;
        EFlux[k][j][i].x = pG->B1i[k][j][i]*bDotGradT/B02;
        EFlux[k][j][i].x *= chi_C;
      }
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute heat fluxes in 2-direction
 */
  
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
/* Compute normalized field at x2-interface */
        Bx = 0.5*(pG->U[k][j-1][i].B1c + pG->U[k][j][i].B1c);
        Bz = 0.5*(pG->U[k][j-1][i].B3c + pG->U[k][j][i].B3c);
        B02 = SQR(Bx) + SQR(pG->B2i[k][j][i]) + SQR(Bz);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized temperature difference dT/dx */
        dTr = 0.5*((Temp[k][j-1][i+1] + Temp[k][j][i+1]) -
                   (Temp[k][j-1][i  ] + Temp[k][j][i  ]));
        dTl = 0.5*((Temp[k][j-1][i  ] + Temp[k][j][i  ]) -
                   (Temp[k][j-1][i-1] + Temp[k][j][i-1]));
        dTc = dTr + dTl;

        dTdx = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
        }

/* Monotonized temperature difference dT/dz */
        dTr = 0.5*((Temp[k+1][j-1][i] + Temp[k+1][j][i]) -
                   (Temp[k  ][j-1][i] + Temp[k  ][j][i]));
        dTl = 0.5*((Temp[k  ][j-1][i] + Temp[k  ][j][i]) -
                   (Temp[k-1][j-1][i] + Temp[k-1][j][i]));
        dTc = dTr + dTl;

        dTdz = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdz = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx3;
        }

        bDotGradT = pG->B2i[k][j][i]*(Temp[k][j][i]-Temp[k][j-1][i])/pG->dx2
           + Bx*dTdx + Bz*dTdz;
        EFlux[k][j][i].y = pG->B2i[k][j][i]*bDotGradT/B02;
        EFlux[k][j][i].y *= chi_C;
      }
    }
  }

/*--- Step 2c ------------------------------------------------------------------
 * Compute heat fluxes in 3-direction
 */
  
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
/* Compute normalized field at x3-interface */
        Bx = 0.5*(pG->U[k-1][j][i].B1c + pG->U[k][j][i].B1c);
        By = 0.5*(pG->U[k-1][j][i].B2c + pG->U[k][j][i].B2c);
        B02 = SQR(Bx) + SQR(By) + SQR(pG->B3i[k][j][i]);
        B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized temperature difference dT/dx */
        dTr = 0.5*((Temp[k-1][j][i+1] + Temp[k][j][i+1]) -
                   (Temp[k-1][j][i  ] + Temp[k][j][i  ]));
        dTl = 0.5*((Temp[k-1][j][i  ] + Temp[k][j][i  ]) -
                   (Temp[k-1][j][i-1] + Temp[k][j][i-1]));
        dTc = dTr + dTl;

        dTdx = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdx = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx1;
        }

/* Monotonized temperature difference dT/dy */
        dTr = 0.5*((Temp[k-1][j+1][i] + Temp[k][j+1][i]) -
                   (Temp[k-1][j  ][i] + Temp[k][j  ][i]));
        dTl = 0.5*((Temp[k-1][j  ][i] + Temp[k][j  ][i]) -
                   (Temp[k-1][j-1][i] + Temp[k][j-1][i]));
        dTc = dTr + dTl;

        dTdy = 0.0;
        if (dTl*dTr > 0.0) {
          lim_slope = MIN(fabs(dTl),fabs(dTr));
          dTdy = SIGN(dTc)*MIN(0.5*fabs(dTc),2.0*lim_slope)/pG->dx2;
        }

        bDotGradT = pG->B3i[k][j][i]*(Temp[k][j][i]-Temp[k-1][j][i])/pG->dx3
           + Bx*dTdx + By*dTdy;
        EFlux[k][j][i].z = pG->B3i[k][j][i]*bDotGradT/B02;
        EFlux[k][j][i].z *= chi_C;
      }
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Update energy using heat-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].E += dtodx1*(EFlux[k][j][i+1].x - EFlux[k][j][i].x);
        pG->U[k][j][i].E += dtodx2*(EFlux[k][j+1][i].y - EFlux[k][j][i].y);
        pG->U[k][j][i].E += dtodx3*(EFlux[k+1][j][i].z - EFlux[k][j][i].z);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* anisoconduct_init: Allocate temporary arrays
*/

void anisoconduct_init(MeshS *pM)
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
  anisoconduct_destruct();
  ath_error("[anisoconduct_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* anisoconduct_destruct: Free temporary arrays
 */

void anisoconduct_destruct(void)
{
  if (Temp != NULL) free_3d_array(Temp);
  if (EFlux != NULL) free_3d_array(EFlux);

  return;
}
#endif /* ANISOTROPIC_CONDUCTION */
