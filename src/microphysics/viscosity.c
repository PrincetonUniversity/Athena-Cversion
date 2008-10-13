#include "../copyright.h"
/*==============================================================================
 * FILE: viscosity.c
 *
 * PURPOSE: Implements explicit Navier-Stokes viscosity, that is
 *      dM/dt = Div(T)    where T = \nu Grad(V) = viscous stress tensor.
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  ns_viscosity_1d()
 *  ns_viscosity_2d()
 *  ns_viscosity_3d()
 *  ns_viscosity_init() - allocates memory needed
 *  ns_viscosity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef NAVIER_STOKES
#ifdef ADIABATIC
#error : Navier-Stokes viscosity only works for isothermal EOS.
#endif /* ADIABATIC */
#endif /* NAVIER_STOKES */

/* The viscous fluxes and velocities, contained in special structures */
typedef struct ViscFlux_t{
  Real Mx;
  Real My;
  Real Mz;
}ViscFlux;
typedef struct ThreeDVect_t{
  Real x;
  Real y;
  Real z;
}ThreeDVect;

static ViscFlux ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
static ThreeDVect ***Vel=NULL;
static Real ***divv=NULL;

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ns_viscosity_1d: Navier-Stokes viscosity in 1d
 */

void ns_viscosity_1d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_2d: Navier-Stokes viscosity in 2d
 */

void ns_viscosity_2d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_3d: Navier-Stokes viscosity in 3d
 */

void ns_viscosity_3d(Grid *pG, Domain *pD)
{
#ifdef NAVIER_STOKES
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real nud;

/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity and div(V) at cell centers (needed everywhere below)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Vel[k][j][i].x = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        Vel[k][j][i].y = pG->U[k][j][i].M2/pG->U[k][j][i].d;
        Vel[k][j][i].z = pG->U[k][j][i].M3/pG->U[k][j][i].d;
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        divv[k][j][i] = ((Vel[k][j][i+1].x - Vel[k][j][i-1].x)/(2.0*pG->dx1) +
                         (Vel[k][j+1][i].y - Vel[k][j-1][i].y)/(2.0*pG->dx2) +
                         (Vel[k+1][j][i].z - Vel[k-1][j][i].z)/(2.0*pG->dx3));
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        x1Flux[k][j][i].Mx = 2.0*(Vel[k][j][i].x - Vel[k][j][i-1].x)/pG->dx1
           - ONE_3RD*(divv[k][j][i] + divv[k][j][i-1]);

        x1Flux[k][j][i].My = (Vel[k][j][i].y - Vel[k][j][i-1].y)/pG->dx1
          + ((Vel[k][j+1][i].x + Vel[k][j+1][i-1].x) - 
             (Vel[k][j-1][i].x + Vel[k][j-1][i-1].x))/(4.0*pG->dx2); 

        x1Flux[k][j][i].Mz = (Vel[k][j][i].z - Vel[k][j][i-1].z)/pG->dx1
          + ((Vel[k+1][j][i].x + Vel[k+1][j][i-1].x) - 
             (Vel[k-1][j][i].x + Vel[k-1][j][i-1].x))/(4.0*pG->dx3);

#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j][i-1].d);
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
        x2Flux[k][j][i].Mx = (Vel[k][j][i].x - Vel[k][j-1][i].x)/pG->dx2
          + ((Vel[k][j][i+1].y + Vel[k][j-1][i+1].y) - 
             (Vel[k][j][i-1].y + Vel[k][j-1][i-1].y))/(4.0*pG->dx1);

        x2Flux[k][j][i].My = 2.0*(Vel[k][j][i].y - Vel[k][j-1][i].y)/pG->dx2
           - ONE_3RD*(divv[k][j][i] + divv[k][j-1][i]);

        x2Flux[k][j][i].Mz = (Vel[k][j][i].z - Vel[k][j-1][i].z)/pG->dx2
          + ((Vel[k+1][j][i].y + Vel[k+1][j-1][i].y) - 
             (Vel[k-1][j][i].y + Vel[k-1][j-1][i].y))/(4.0*pG->dx3);
#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j-1][i].d);
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
        x3Flux[k][j][i].Mx = (Vel[k][j][i].x - Vel[k-1][j][i].x)/pG->dx3
          + ((Vel[k][j][i+1].z + Vel[k-1][j][i+1].z) -
             (Vel[k][j][i-1].z + Vel[k-1][j][i-1].z))/(4.0*pG->dx1);

        x3Flux[k][j][i].My = (Vel[k][j][i].y - Vel[k-1][j][i].y)/pG->dx3
          + ((Vel[k][j+1][i].z + Vel[k-1][j+1][i].z) -
             (Vel[k][j-1][i].z + Vel[k-1][j-1][i].z))/(4.0*pG->dx2);

        x3Flux[k][j][i].Mz = 2.0*(Vel[k][j][i].z - Vel[k-1][j][i].z)/pG->dx3
           - ONE_3RD*(divv[k][j][i] + divv[k-1][j][i]);
#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
#endif /* BAROTROPIC */

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k-1][j][i].d);
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

#endif /* NAVIER_STOKES */
  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_init: Allocate temporary arrays
 */

void ns_viscosity_init(int nx1, int nx2, int nx3)
{
#ifdef NAVIER_STOKES
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
    if ((x1Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
      == NULL) goto on_error;
    if ((x2Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
      == NULL) goto on_error;
    if ((x3Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
      == NULL) goto on_error;
    if ((Vel = (ThreeDVect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ThreeDVect)))
      == NULL) goto on_error;
    if ((divv = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
      goto on_error;
    break;
  }
  return;

  on_error:
  ns_viscosity_destruct();
  ath_error("[ns_viscosity_init]: malloc returned a NULL pointer\n");
#endif /* NAVIER_STOKES */
  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_destruct: Free temporary arrays
 */      

void ns_viscosity_destruct(void)
{   
#ifdef NAVIER_STOKES
/* dim set in ns_viscosity_init() at begnning of run */
  switch(dim){
  case 1:
    break;
  case 2:
    break;
  case 3:
    if (x1Flux != NULL) free_3d_array(x1Flux);
    if (x2Flux != NULL) free_3d_array(x2Flux);
    if (x3Flux != NULL) free_3d_array(x3Flux);
    if (Vel != NULL) free_3d_array(Vel);
    if (divv != NULL) free_3d_array(divv);
    break;
  }
#endif /* NAVIER_STOKES */

  return;
}
