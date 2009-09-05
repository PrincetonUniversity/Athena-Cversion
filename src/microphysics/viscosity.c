#include "../copyright.h"
/*==============================================================================
 * FILE: viscosity.c
 *
 * PURPOSE: Implements explicit Navier-Stokes viscosity, that is
 *      dM/dt = Div(T)    where T = \nu Grad(V) = viscous stress tensor.
 *      dE/dt = Div(v.T)
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
/* The viscous fluxes and velocities, contained in special structures */
typedef struct ViscFlux_t{
  Real Mx;
  Real My;
  Real Mz;
#ifndef BAROTROPIC
  Real E;
#endif
}ViscFlux;

static ViscFlux ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
static Vector ***Vel=NULL;
static Real ***divv=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ns_viscosity_1d: Navier-Stokes viscosity in 1d
 */

void ns_viscosity_1d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real nud;

/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity and div(V) at cell centers (needed everywhere below)
 */

  for (i=is-2; i<=ie+2; i++) {
    Vel[ks][js][i].x1 = pG->U[ks][js][i].M1/pG->U[ks][js][i].d;
    Vel[ks][js][i].x2 = pG->U[ks][js][i].M2/pG->U[ks][js][i].d;
    Vel[ks][js][i].x3 = pG->U[ks][js][i].M3/pG->U[ks][js][i].d;
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (i=is; i<=ie+1; i++) {
    x1Flux[ks][js][i].Mx = (Vel[ks][js][i].x1 - Vel[ks][js][i-1].x1)/pG->dx1;
    x1Flux[ks][js][i].My = (Vel[ks][js][i].x2 - Vel[ks][js][i-1].x2)/pG->dx1;
    x1Flux[ks][js][i].Mz = (Vel[ks][js][i].x3 - Vel[ks][js][i-1].x3)/pG->dx1;

    nud = nu_V*0.5*(pG->U[ks][js][i].d + pG->U[ks][js][i-1].d);
    x1Flux[ks][js][i].Mx *= (4./3.)*nud;
    x1Flux[ks][js][i].My *= nud;
    x1Flux[ks][js][i].Mz *= nud;

#ifndef BAROTROPIC
    x1Flux[ks][js][i].E  =
       0.5*(Vel[ks][js][i-1].x1 + Vel[ks][js][i].x1)*x1Flux[ks][js][i].Mx +
       0.5*(Vel[ks][js][i-1].x2 + Vel[ks][js][i].x2)*x1Flux[ks][js][i].My +
       0.5*(Vel[ks][js][i-1].x3 + Vel[ks][js][i].x3)*x1Flux[ks][js][i].Mz;
#endif /* BAROTROPIC */
  }

/*--- Step 3a ------------------------------------------------------------------
 * Update momentum and energy using x1-fluxes (dM/dt = Div(T))
 */

  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].M1 += dtodx1*(x1Flux[ks][js][i+1].Mx-x1Flux[ks][js][i].Mx);
    pG->U[ks][js][i].M2 += dtodx1*(x1Flux[ks][js][i+1].My-x1Flux[ks][js][i].My);
    pG->U[ks][js][i].M3 += dtodx1*(x1Flux[ks][js][i+1].Mz-x1Flux[ks][js][i].Mz);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E  += dtodx1*(x1Flux[ks][js][i+1].E -x1Flux[ks][js][i].E );
#endif /* BAROTROPIC */
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_2d: Navier-Stokes viscosity in 2d
 */

void ns_viscosity_2d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real nud;

/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity and div(V) at cell centers (needed everywhere below)
 */

  for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      Vel[ks][j][i].x1 = pG->U[ks][j][i].M1/pG->U[ks][j][i].d;
      Vel[ks][j][i].x2 = pG->U[ks][j][i].M2/pG->U[ks][j][i].d;
      Vel[ks][j][i].x3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;
    }
  }

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      divv[ks][j][i] = ((Vel[ks][j][i+1].x1 -Vel[ks][j][i-1].x1)/(2.0*pG->dx1) +
                        (Vel[ks][j+1][i].x2 -Vel[ks][j-1][i].x2)/(2.0*pG->dx2));
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {
      x1Flux[ks][j][i].Mx = 2.0*(Vel[ks][j][i].x1 - Vel[ks][j][i-1].x1)/pG->dx1
         - ONE_3RD*(divv[ks][j][i] + divv[ks][j][i-1]);

      x1Flux[ks][j][i].My = (Vel[ks][j][i].x2 - Vel[ks][j][i-1].x2)/pG->dx1
        + ((Vel[ks][j+1][i].x1 + Vel[ks][j+1][i-1].x1) - 
           (Vel[ks][j-1][i].x1 + Vel[ks][j-1][i-1].x1))/(4.0*pG->dx2); 

      x1Flux[ks][j][i].Mz = (Vel[ks][j][i].x3 - Vel[ks][j][i-1].x3)/pG->dx1;

      nud = nu_V*0.5*(pG->U[ks][j][i].d + pG->U[ks][j][i-1].d);
      x1Flux[ks][j][i].Mx *= nud;
      x1Flux[ks][j][i].My *= nud;
      x1Flux[ks][j][i].Mz *= nud;

#ifndef BAROTROPIC
      x1Flux[ks][j][i].E  =
         0.5*(Vel[ks][j][i-1].x1 + Vel[ks][j][i].x1)*x1Flux[ks][j][i].Mx +
         0.5*(Vel[ks][j][i-1].x2 + Vel[ks][j][i].x2)*x1Flux[ks][j][i].My +
         0.5*(Vel[ks][j][i-1].x3 + Vel[ks][j][i].x3)*x1Flux[ks][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute viscous fluxes in 2-direction
 */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      x2Flux[ks][j][i].Mx = (Vel[ks][j][i].x1 - Vel[ks][j-1][i].x1)/pG->dx2
        + ((Vel[ks][j][i+1].x2 + Vel[ks][j-1][i+1].x2) - 
           (Vel[ks][j][i-1].x2 + Vel[ks][j-1][i-1].x2))/(4.0*pG->dx1);

      x2Flux[ks][j][i].My = 2.0*(Vel[ks][j][i].x2 - Vel[ks][j-1][i].x2)/pG->dx2
         - ONE_3RD*(divv[ks][j][i] + divv[ks][j-1][i]);

      x2Flux[ks][j][i].Mz = (Vel[ks][j][i].x3 - Vel[ks][j-1][i].x3)/pG->dx2;

      nud = nu_V*0.5*(pG->U[ks][j][i].d + pG->U[ks][j-1][i].d);
      x2Flux[ks][j][i].Mx *= nud;
      x2Flux[ks][j][i].My *= nud;
      x2Flux[ks][j][i].Mz *= nud;

#ifndef BAROTROPIC
      x2Flux[ks][j][i].E  =
         0.5*(Vel[ks][j-1][i].x1 + Vel[ks][j][i].x1)*x2Flux[ks][j][i].Mx +
         0.5*(Vel[ks][j-1][i].x2 + Vel[ks][j][i].x2)*x2Flux[ks][j][i].My +
         0.5*(Vel[ks][j-1][i].x3 + Vel[ks][j][i].x3)*x2Flux[ks][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }

/*--- Step 3a ------------------------------------------------------------------
 * Update momentum and energy using x1-fluxes (dM/dt = Div(T))
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].M1 += dtodx1*(x1Flux[ks][j][i+1].Mx-x1Flux[ks][j][i].Mx);
      pG->U[ks][j][i].M2 += dtodx1*(x1Flux[ks][j][i+1].My-x1Flux[ks][j][i].My);
      pG->U[ks][j][i].M3 += dtodx1*(x1Flux[ks][j][i+1].Mz-x1Flux[ks][j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  += dtodx1*(x1Flux[ks][j][i+1].E -x1Flux[ks][j][i].E );
#endif /* BAROTROPIC */
    }
  }

/*--- Step 3b ------------------------------------------------------------------
 * Update momentum and energy using x2-fluxes (dM/dt = Div(T))
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].M1 += dtodx2*(x2Flux[ks][j+1][i].Mx-x2Flux[ks][j][i].Mx);
      pG->U[ks][j][i].M2 += dtodx2*(x2Flux[ks][j+1][i].My-x2Flux[ks][j][i].My);
      pG->U[ks][j][i].M3 += dtodx2*(x2Flux[ks][j+1][i].Mz-x2Flux[ks][j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  += dtodx2*(x2Flux[ks][j+1][i].E -x2Flux[ks][j][i].E );
#endif /* BAROTROPIC */
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_3d: Navier-Stokes viscosity in 3d
 */

void ns_viscosity_3d(Grid *pG, Domain *pD)
{
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
        Vel[k][j][i].x1 = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        Vel[k][j][i].x2 = pG->U[k][j][i].M2/pG->U[k][j][i].d;
        Vel[k][j][i].x3 = pG->U[k][j][i].M3/pG->U[k][j][i].d;
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        divv[k][j][i] = ((Vel[k][j][i+1].x1 - Vel[k][j][i-1].x1)/(2.0*pG->dx1) +
                         (Vel[k][j+1][i].x2 - Vel[k][j-1][i].x2)/(2.0*pG->dx2) +
                         (Vel[k+1][j][i].x3 - Vel[k-1][j][i].x3)/(2.0*pG->dx3));
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+1; i++) {
        x1Flux[k][j][i].Mx = 2.0*(Vel[k][j][i].x1 - Vel[k][j][i-1].x1)/pG->dx1
           - ONE_3RD*(divv[k][j][i] + divv[k][j][i-1]);

        x1Flux[k][j][i].My = (Vel[k][j][i].x2 - Vel[k][j][i-1].x2)/pG->dx1
          + ((Vel[k][j+1][i].x1 + Vel[k][j+1][i-1].x1) - 
             (Vel[k][j-1][i].x1 + Vel[k][j-1][i-1].x1))/(4.0*pG->dx2); 

        x1Flux[k][j][i].Mz = (Vel[k][j][i].x3 - Vel[k][j][i-1].x3)/pG->dx1
          + ((Vel[k+1][j][i].x1 + Vel[k+1][j][i-1].x1) - 
             (Vel[k-1][j][i].x1 + Vel[k-1][j][i-1].x1))/(4.0*pG->dx3);

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j][i-1].d);
        x1Flux[k][j][i].Mx *= nud;
        x1Flux[k][j][i].My *= nud;
        x1Flux[k][j][i].Mz *= nud;

#ifndef BAROTROPIC
        x1Flux[k][j][i].E  =
           0.5*(Vel[k][j][i-1].x1 + Vel[k][j][i].x1)*x1Flux[k][j][i].Mx +
           0.5*(Vel[k][j][i-1].x2 + Vel[k][j][i].x2)*x1Flux[k][j][i].My +
           0.5*(Vel[k][j][i-1].x3 + Vel[k][j][i].x3)*x1Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
      }
    }
  }

/*--- Step 2b ------------------------------------------------------------------
 * Compute viscous fluxes in 2-direction
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie; i++) {
        x2Flux[k][j][i].Mx = (Vel[k][j][i].x1 - Vel[k][j-1][i].x1)/pG->dx2
          + ((Vel[k][j][i+1].x2 + Vel[k][j-1][i+1].x2) - 
             (Vel[k][j][i-1].x2 + Vel[k][j-1][i-1].x2))/(4.0*pG->dx1);

        x2Flux[k][j][i].My = 2.0*(Vel[k][j][i].x2 - Vel[k][j-1][i].x2)/pG->dx2
           - ONE_3RD*(divv[k][j][i] + divv[k][j-1][i]);

        x2Flux[k][j][i].Mz = (Vel[k][j][i].x3 - Vel[k][j-1][i].x3)/pG->dx2
          + ((Vel[k+1][j][i].x2 + Vel[k+1][j-1][i].x2) - 
             (Vel[k-1][j][i].x2 + Vel[k-1][j-1][i].x2))/(4.0*pG->dx3);

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j-1][i].d);
        x2Flux[k][j][i].Mx *= nud;
        x2Flux[k][j][i].My *= nud;
        x2Flux[k][j][i].Mz *= nud;

#ifndef BAROTROPIC
        x2Flux[k][j][i].E  =
           0.5*(Vel[k][j-1][i].x1 + Vel[k][j][i].x1)*x2Flux[k][j][i].Mx +
           0.5*(Vel[k][j-1][i].x2 + Vel[k][j][i].x2)*x2Flux[k][j][i].My +
           0.5*(Vel[k][j-1][i].x3 + Vel[k][j][i].x3)*x2Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
      }
    }
  }

/*--- Step 2c ------------------------------------------------------------------
 * Compute viscous fluxes in 3-direction
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        x3Flux[k][j][i].Mx = (Vel[k][j][i].x1 - Vel[k-1][j][i].x1)/pG->dx3
          + ((Vel[k][j][i+1].x3 + Vel[k-1][j][i+1].x3) -
             (Vel[k][j][i-1].x3 + Vel[k-1][j][i-1].x3))/(4.0*pG->dx1);

        x3Flux[k][j][i].My = (Vel[k][j][i].x2 - Vel[k-1][j][i].x2)/pG->dx3
          + ((Vel[k][j+1][i].x3 + Vel[k-1][j+1][i].x3) -
             (Vel[k][j-1][i].x3 + Vel[k-1][j-1][i].x3))/(4.0*pG->dx2);

        x3Flux[k][j][i].Mz = 2.0*(Vel[k][j][i].x3 - Vel[k-1][j][i].x3)/pG->dx3
           - ONE_3RD*(divv[k][j][i] + divv[k-1][j][i]);

        nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k-1][j][i].d);
        x3Flux[k][j][i].Mx *= nud;
        x3Flux[k][j][i].My *= nud;
        x3Flux[k][j][i].Mz *= nud;

#ifndef BAROTROPIC
        x3Flux[k][j][i].E  =
           0.5*(Vel[k-1][j][i].x1 + Vel[k][j][i].x1)*x3Flux[k][j][i].Mx +
           0.5*(Vel[k-1][j][i].x2 + Vel[k][j][i].x2)*x3Flux[k][j][i].My +
           0.5*(Vel[k-1][j][i].x3 + Vel[k][j][i].x3)*x3Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
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
        pG->U[k][j][i].E  += dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
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

  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_init: Allocate temporary arrays
 */

void ns_viscosity_init(int nx1, int nx2, int nx3)
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

  if ((x1Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
    == NULL) goto on_error;
  if ((x2Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
    == NULL) goto on_error;
  if ((x3Flux = (ViscFlux***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ViscFlux)))
    == NULL) goto on_error;
  if ((Vel = (Vector***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Vector)))
    == NULL) goto on_error;
  if ((divv = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  return;

  on_error:
  ns_viscosity_destruct();
  ath_error("[ns_viscosity_init]: malloc returned a NULL pointer\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* ns_viscosity_destruct: Free temporary arrays
 */      

void ns_viscosity_destruct(void)
{   
  if (x1Flux != NULL) free_3d_array(x1Flux);
  if (x2Flux != NULL) free_3d_array(x2Flux);
  if (x3Flux != NULL) free_3d_array(x3Flux);
  if (Vel != NULL) free_3d_array(Vel);
  if (divv != NULL) free_3d_array(divv);
  return;
}
#endif /* NAVIER_STOKES */
