#include "../copyright.h"
/*==============================================================================
 * FILE: braginskii.c
 *
 * PURPOSE: Implements anisotropic (Braginskii) viscosity, that is
 *      dM/dt = Div(T)    where T=Braginski stress tensor.
 *      dE/dt = Div(v.T)
 *   Functions are called by integrate_diffusion() in the main loop, which
 *   coordinates adding all diffusion operators (viscosity, resistivity, thermal
 *   conduction) using operator splitting.
 *
 *   An explicit timestep limit must be applied if these routines are used.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *  brag_viscosity_2d()
 *  brag_viscosity_3d()
 *  brag_viscosity_init() - allocates memory needed
 *  brag_viscosity_destruct() - frees memory used
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef BRAGINSKII
#ifdef HYDRO
#error : Braginskii viscosity only works for MHD.
#endif /* HYDRO */
#endif /* BRAGINSKII */

/* The viscous fluxes, contained in special structures */
typedef struct ViscFlux_t{
  Real Mx;
  Real My;
  Real Mz;
#ifndef BAROTROPIC
  Real E;
#endif
}ViscFlux;
typedef struct ThreeDVect_t{
  Real x;
  Real y;
  Real z;
}ThreeDVect;

static ViscFlux ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
static ThreeDVect ***Vel=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* brag_viscosity_2d: Braginskii viscosity in 2d
 */

void brag_viscosity_2d(Grid *pG, Domain *pD)
{
#ifdef BRAGINSKII
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real Bx,By,Bz,B02,dVc,dVl,dVr,lim_slope,BBdV,divV,qa,nud;
  Real dVxdx,dVydx,dVzdx,dVxdy,dVydy,dVzdy;
#ifdef FARGO
  Real x1,x2,x3;
#endif

/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity. In shearing box Vy must include background shear, because
 * coefficients of dVy/dx terms are not constant.
 */

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      Vel[ks][j][i].x = pG->U[ks][j][i].M1/pG->U[ks][j][i].d;
      Vel[ks][j][i].y = pG->U[ks][j][i].M2/pG->U[ks][j][i].d;
#ifdef FARGO
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);
      Vel[ks][j][i].y -= qshear*Omega_0*x1;
#endif
      Vel[ks][j][i].z = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction, centered at X1-Faces
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {

/* Compute field components at x1-interface */

      Bx = pG->B1i[ks][j][i];
      By = 0.5*(pG->U[ks][j][i].B2c + pG->U[ks][j][i-1].B2c);
      Bz = 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j][i-1].B3c);
      B02 = Bx*Bx + By*By + Bz*Bz;
      B02 = MAX(B02,TINY_NUMBER);  /* limit in case B=0 */

/* Monotonized Velocity gradient dVx/dy */
      dVr = 0.5*((Vel[ks][j+1][i-1].x + Vel[ks][j+1][i].x) -
                 (Vel[ks][j  ][i-1].x + Vel[ks][j  ][i].x));
      dVl = 0.5*((Vel[ks][j  ][i-1].x + Vel[ks][j  ][i].x) -
                 (Vel[ks][j-1][i-1].x + Vel[ks][j-1][i].x));
      dVc = dVr + dVl;
  
      dVxdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVy/dy */
      dVr = 0.5*((Vel[ks][j+1][i-1].y + Vel[ks][j+1][i].y) -
                 (Vel[ks][j  ][i-1].y + Vel[ks][j  ][i].y));
      dVl = 0.5*((Vel[ks][j  ][i-1].y + Vel[ks][j  ][i].y) -
                 (Vel[ks][j-1][i-1].y + Vel[ks][j-1][i].y));
      dVc = dVr + dVl;
  
      dVydy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVz/dy */
      dVr = 0.5*((Vel[ks][j+1][i-1].z + Vel[ks][j+1][i].z) -
                 (Vel[ks][j  ][i-1].z + Vel[ks][j  ][i].z));
      dVl = 0.5*((Vel[ks][j  ][i-1].z + Vel[ks][j  ][i].z) -
                 (Vel[ks][j-1][i-1].z + Vel[ks][j-1][i].z));
      dVc = dVr + dVl;
  
      dVzdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* compute BBdV and div(V) */

      BBdV = 
        Bx/B02*(Bx*(Vel[ks][j][i].x-Vel[ks][j][i-1].x)/pG->dx1 + By*dVxdy) +
        By/B02*(Bx*(Vel[ks][j][i].y-Vel[ks][j][i-1].y)/pG->dx1 + By*dVydy) +
        Bz/B02*(Bx*(Vel[ks][j][i].z-Vel[ks][j][i-1].z)/pG->dx1 + By*dVzdy);

      divV = (Vel[ks][j][i].x-Vel[ks][j][i-1].x)/pG->dx1 + dVydy;
      nud = nu_V*0.5*(pG->U[ks][j][i].d + pG->U[ks][j][i-1].d);

      qa = nud*(BBdV - ONE_3RD*divV);

/* Compute fluxes */

      x1Flux[ks][j][i].Mx = qa*(3.0*Bx*Bx/B02 - 1.0);
      x1Flux[ks][j][i].My = qa*(3.0*By*Bx/B02);
      x1Flux[ks][j][i].Mz = qa*(3.0*Bz*Bx/B02);

#ifndef BAROTROPIC
      x1Flux[ks][j][i].E =
         0.5*(Vel[ks][j][i-1].x + Vel[ks][j][i].x)*x1Flux[ks][j][i].Mx +
         0.5*(Vel[ks][j][i-1].y + Vel[ks][j][i].y)*x1Flux[ks][j][i].My +
         0.5*(Vel[ks][j][i-1].z + Vel[ks][j][i].z)*x1Flux[ks][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }

/*--- Step 3 -------------------------------------------------------------------
 * Compute viscous fluxes in 2-direction, centered at X2-Faces
 */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {

/* Compute field components at x2-interface */

      Bx = 0.5*(pG->U[ks][j][i].B1c + pG->U[ks][j-1][i].B1c);
      By = pG->B2i[ks][j][i];
      Bz = 0.5*(pG->U[ks][j][i].B3c + pG->U[ks][j-1][i].B3c);
      B02 = Bx*Bx + By*By + Bz*Bz;
      B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized Velocity gradient dVx/dx */
      dVr = 0.5*((Vel[ks][j-1][i+1].x + Vel[ks][j][i+1].x) -
                 (Vel[ks][j-1][i  ].x + Vel[ks][j][i  ].x));
      dVl = 0.5*((Vel[ks][j-1][i  ].x + Vel[ks][j][i  ].x) -
                 (Vel[ks][j-1][i-1].x + Vel[ks][j][i-1].x));
      dVc = dVr + dVl;

      dVxdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVy/dx */
      dVr = 0.5*((Vel[ks][j-1][i+1].y + Vel[ks][j][i+1].y) -
                 (Vel[ks][j-1][i  ].y + Vel[ks][j][i  ].y));
      dVl = 0.5*((Vel[ks][j-1][i  ].y + Vel[ks][j][i  ].y) -
                 (Vel[ks][j-1][i-1].y + Vel[ks][j][i-1].y));
      dVc = dVr + dVl;

      dVydx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVz/dx */
      dVr = 0.5*((Vel[ks][j-1][i+1].z + Vel[ks][j][i+1].z) -
                 (Vel[ks][j-1][i  ].z + Vel[ks][j][i  ].z));
      dVl = 0.5*((Vel[ks][j-1][i  ].z + Vel[ks][j][i  ].z) -
                 (Vel[ks][j-1][i-1].z + Vel[ks][j][i-1].z));
      dVc = dVr + dVl;

      dVzdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* compute BBdV and div(V) */

      BBdV =
        Bx/B02*(Bx*dVxdx + By*(Vel[ks][j][i].x-Vel[ks][j-1][i].x)/pG->dx2) +
        By/B02*(Bx*dVydx + By*(Vel[ks][j][i].y-Vel[ks][j-1][i].y)/pG->dx2) +
        Bz/B02*(Bx*dVzdx + By*(Vel[ks][j][i].z-Vel[ks][j-1][i].z)/pG->dx2);

      divV = dVxdx + (Vel[ks][j][i].y-Vel[ks][j-1][i].y)/pG->dx2;
      nud = nu_V*0.5*(pG->U[ks][j][i].d + pG->U[ks][j-1][i].d);

      qa = nud*(BBdV - ONE_3RD*divV);

/* Compute fluxes */

      x2Flux[ks][j][i].Mx = qa*(3.0*Bx*By/B02);
      x2Flux[ks][j][i].My = qa*(3.0*By*By/B02 - 1.0);
      x2Flux[ks][j][i].Mz = qa*(3.0*Bz*By/B02);

#ifndef BAROTROPIC
        x2Flux[ks][j][i].E =
           0.5*(Vel[ks][j-1][i].x + Vel[ks][j][i].x)*x2Flux[ks][j][i].Mx +
           0.5*(Vel[ks][j-1][i].y + Vel[ks][j][i].y)*x2Flux[ks][j][i].My +
           0.5*(Vel[ks][j-1][i].z + Vel[ks][j][i].z)*x2Flux[ks][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }

/*--- Step 4a ------------------------------------------------------------------
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

/*--- Step 4b ------------------------------------------------------------------
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

#endif /* BRAGINSKII */
  return;
}


/*----------------------------------------------------------------------------*/
/* brag_viscosity_3d: Braginskii viscosity in 3d
 */

void brag_viscosity_3d(Grid *pG, Domain *pD)
{
#ifdef BRAGINSKII
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real Bx,By,Bz,B02,dVc,dVl,dVr,lim_slope,BBdV,divV,qa,nud;
  Real dVxdx,dVydx,dVzdx,dVxdy,dVydy,dVzdy,dVxdz,dVydz,dVzdz;
#ifdef FARGO
  Real x1,x2,x3;
#endif

/*--- Step 1 -------------------------------------------------------------------
 * Compute velocity. In shearing box Vy must include background shear, because
 * coefficients of dVy/dx terms are not constant.
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
        Vel[k][j][i].x = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        Vel[k][j][i].y = pG->U[k][j][i].M2/pG->U[k][j][i].d;
#ifdef FARGO
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        Vel[k][j][i].y -= qshear*Omega_0*x1;
#endif
        Vel[k][j][i].z = pG->U[k][j][i].M3/pG->U[k][j][i].d;
      }
    }
  }

/*--- Step 2 -------------------------------------------------------------------
 * Compute viscous fluxes in 1-direction, centered at x1-Faces
 */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie+1; i++) {

/* Compute field components at x1-interface */

      Bx = pG->B1i[k][j][i];
      By = 0.5*(pG->U[k][j][i].B2c + pG->U[k][j][i-1].B2c);
      Bz = 0.5*(pG->U[k][j][i].B3c + pG->U[k][j][i-1].B3c);
      B02 = Bx*Bx + By*By + Bz*Bz;
      B02 = MAX(B02,TINY_NUMBER);  /* limit in case B=0 */

/* Monotonized Velocity gradient dVx/dy */
      dVr = 0.5*((Vel[k][j+1][i-1].x + Vel[k][j+1][i].x) -
                 (Vel[k][j  ][i-1].x + Vel[k][j  ][i].x));
      dVl = 0.5*((Vel[k][j  ][i-1].x + Vel[k][j  ][i].x) -
                 (Vel[k][j-1][i-1].x + Vel[k][j-1][i].x));
      dVc = dVr + dVl;

      dVxdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVy/dy */
      dVr = 0.5*((Vel[k][j+1][i-1].y + Vel[k][j+1][i].y) -
                 (Vel[k][j  ][i-1].y + Vel[k][j  ][i].y));
      dVl = 0.5*((Vel[k][j  ][i-1].y + Vel[k][j  ][i].y) -
                 (Vel[k][j-1][i-1].y + Vel[k][j-1][i].y));
      dVc = dVr + dVl;

      dVydy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVz/dy */
      dVr = 0.5*((Vel[k][j+1][i-1].z + Vel[k][j+1][i].z) -
                 (Vel[k][j  ][i-1].z + Vel[k][j  ][i].z));
      dVl = 0.5*((Vel[k][j  ][i-1].z + Vel[k][j  ][i].z) -
                 (Vel[k][j-1][i-1].z + Vel[k][j-1][i].z));
      dVc = dVr + dVl;

      dVzdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVx/dz */
      dVr = 0.5*((Vel[k+1][j][i-1].x + Vel[k+1][j][i].x) -
                 (Vel[k  ][j][i-1].x + Vel[k  ][j][i].x));
      dVl = 0.5*((Vel[k  ][j][i-1].x + Vel[k  ][j][i].x) -
                 (Vel[k-1][j][i-1].x + Vel[k-1][j][i].x));
      dVc = dVr + dVl;

      dVxdz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* Monotonized Velocity gradient dVy/dz */
      dVr = 0.5*((Vel[k+1][j][i-1].y + Vel[k+1][j][i].y) -
                 (Vel[k  ][j][i-1].y + Vel[k  ][j][i].y));
      dVl = 0.5*((Vel[k  ][j][i-1].y + Vel[k  ][j][i].y) -
                 (Vel[k-1][j][i-1].y + Vel[k-1][j][i].y));
      dVc = dVr + dVl;

      dVydz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* Monotonized Velocity gradient dVz/dz */
      dVr = 0.5*((Vel[k+1][j][i-1].z + Vel[k+1][j][i].z) -
                 (Vel[k  ][j][i-1].z + Vel[k  ][j][i].z));
      dVl = 0.5*((Vel[k  ][j][i-1].z + Vel[k  ][j][i].z) -
                 (Vel[k-1][j][i-1].z + Vel[k-1][j][i].z));
      dVc = dVr + dVl;

      dVzdz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* compute BBdV and div(V) */

      BBdV =
        Bx*(Bx*(Vel[k][j][i].x-Vel[k][j][i-1].x)/pG->dx1 + By*dVxdy + Bz*dVxdz)+
        By*(Bx*(Vel[k][j][i].y-Vel[k][j][i-1].y)/pG->dx1 + By*dVydy + Bz*dVydz)+
        Bz*(Bx*(Vel[k][j][i].z-Vel[k][j][i-1].z)/pG->dx1 + By*dVzdy + Bz*dVzdz);
      BBdV /= B02;

      divV = (Vel[k][j][i].x-Vel[k][j][i-1].x)/pG->dx1 + dVydy + dVzdz;
      nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j][i-1].d);

      qa = nud*(BBdV - ONE_3RD*divV);

/* Compute fluxes */

      x1Flux[k][j][i].Mx = qa*(3.0*Bx*Bx/B02 - 1.0);
      x1Flux[k][j][i].My = qa*(3.0*By*Bx/B02);
      x1Flux[k][j][i].Mz = qa*(3.0*Bz*Bx/B02);

#ifndef BAROTROPIC
      x1Flux[k][j][i].E =
         0.5*(Vel[k][j][i-1].x + Vel[k][j][i].x)*x1Flux[k][j][i].Mx +
         0.5*(Vel[k][j][i-1].y + Vel[k][j][i].y)*x1Flux[k][j][i].My +
         0.5*(Vel[k][j][i-1].z + Vel[k][j][i].z)*x1Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }}

/*--- Step 3 -------------------------------------------------------------------
 * Compute viscous fluxes in 2-direction, centered at X2-Faces
 */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {

/* Compute field components at x2-interface */

      Bx = 0.5*(pG->U[k][j][i].B1c + pG->U[k][j-1][i].B1c);
      By = pG->B2i[k][j][i];
      Bz = 0.5*(pG->U[k][j][i].B3c + pG->U[k][j-1][i].B3c);
      B02 = Bx*Bx + By*By + Bz*Bz;
      B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized Velocity gradient dVx/dx */
      dVr = 0.5*((Vel[k][j-1][i+1].x + Vel[k][j][i+1].x) -
                 (Vel[k][j-1][i  ].x + Vel[k][j][i  ].x));
      dVl = 0.5*((Vel[k][j-1][i  ].x + Vel[k][j][i  ].x) -
                 (Vel[k][j-1][i-1].x + Vel[k][j][i-1].x));
      dVc = dVr + dVl;

      dVxdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVy/dx */
      dVr = 0.5*((Vel[k][j-1][i+1].y + Vel[k][j][i+1].y) -
                 (Vel[k][j-1][i  ].y + Vel[k][j][i  ].y));
      dVl = 0.5*((Vel[k][j-1][i  ].y + Vel[k][j][i  ].y) -
                 (Vel[k][j-1][i-1].y + Vel[k][j][i-1].y));
      dVc = dVr + dVl;

      dVydx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVz/dx */
      dVr = 0.5*((Vel[k][j-1][i+1].z + Vel[k][j][i+1].z) -
                 (Vel[k][j-1][i  ].z + Vel[k][j][i  ].z));
      dVl = 0.5*((Vel[k][j-1][i  ].z + Vel[k][j][i  ].z) -
                 (Vel[k][j-1][i-1].z + Vel[k][j][i-1].z));
      dVc = dVr + dVl;

      dVzdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVx/dz */
      dVr = 0.5*((Vel[k+1][j-1][i].x + Vel[k+1][j][i].x) -
                 (Vel[k  ][j-1][i].x + Vel[k  ][j][i].x));
      dVl = 0.5*((Vel[k  ][j-1][i].x + Vel[k  ][j][i].x) -
                 (Vel[k-1][j-1][i].x + Vel[k-1][j][i].x));
      dVc = dVr + dVl;

      dVxdz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* Monotonized Velocity gradient dVy/dz */
      dVr = 0.5*((Vel[k+1][j-1][i].y + Vel[k+1][j][i].y) -
                 (Vel[k  ][j-1][i].y + Vel[k  ][j][i].y));
      dVl = 0.5*((Vel[k  ][j-1][i].y + Vel[k  ][j][i].y) -
                 (Vel[k-1][j-1][i].y + Vel[k-1][j][i].y));
      dVc = dVr + dVl;

      dVydz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* Monotonized Velocity gradient dVz/dz */
      dVr = 0.5*((Vel[k+1][j-1][i].z + Vel[k+1][j][i].z) -
                 (Vel[k  ][j-1][i].z + Vel[k  ][j][i].z));
      dVl = 0.5*((Vel[k  ][j-1][i].z + Vel[k  ][j][i].z) -
                 (Vel[k-1][j-1][i].z + Vel[k-1][j][i].z));
      dVc = dVr + dVl;

      dVzdz = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdz = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx3;
      }

/* compute BBdV and div(V) */

      BBdV =
        Bx*(Bx*dVxdx + By*(Vel[k][j][i].x-Vel[k][j-1][i].x)/pG->dx2 + Bz*dVxdz)+
        By*(Bx*dVydx + By*(Vel[k][j][i].y-Vel[k][j-1][i].y)/pG->dx2 + Bz*dVydz)+
        Bz*(Bx*dVzdx + By*(Vel[k][j][i].z-Vel[k][j-1][i].z)/pG->dx2 + Bz*dVzdz);
      BBdV /= B02;

      divV = dVxdx + (Vel[k][j][i].y-Vel[k][j-1][i].y)/pG->dx2 + dVzdz;
      nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k][j-1][i].d);

      qa = nud*(BBdV - ONE_3RD*divV);

/* Compute fluxes */

      x2Flux[k][j][i].Mx = qa*(3.0*Bx*By/B02);
      x2Flux[k][j][i].My = qa*(3.0*By*By/B02 - 1.0);
      x2Flux[k][j][i].Mz = qa*(3.0*Bz*By/B02);

#ifndef BAROTROPIC
        x2Flux[k][j][i].E =
           0.5*(Vel[k][j-1][i].x + Vel[k][j][i].x)*x2Flux[k][j][i].Mx +
           0.5*(Vel[k][j-1][i].y + Vel[k][j][i].y)*x2Flux[k][j][i].My +
           0.5*(Vel[k][j-1][i].z + Vel[k][j][i].z)*x2Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
    }
  }}

/*--- Step 4 -------------------------------------------------------------------
 * Compute viscous fluxes in 3-direction, centered at x3-Faces
 */

  for (k=ks; k<=ke+1; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {

/* Compute field components at x3-interface */

      Bx = 0.5*(pG->U[k][j][i].B1c + pG->U[k-1][j][i].B1c);
      By = 0.5*(pG->U[k][j][i].B2c + pG->U[k-1][j][i].B2c);
      Bz = pG->B3i[k][j][i];
      B02 = Bx*Bx + By*By + Bz*Bz;
      B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

/* Monotonized Velocity gradient dVx/dx */
      dVr = 0.5*((Vel[k-1][j][i+1].x + Vel[k][j][i+1].x) -
                 (Vel[k-1][j][i  ].x + Vel[k][j][i  ].x));
      dVl = 0.5*((Vel[k-1][j][i  ].x + Vel[k][j][i  ].x) -
                 (Vel[k-1][j][i-1].x + Vel[k][j][i-1].x));
      dVc = dVr + dVl;

      dVxdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVy/dx */
      dVr = 0.5*((Vel[k-1][j][i+1].y + Vel[k][j][i+1].y) -
                 (Vel[k-1][j][i  ].y + Vel[k][j][i  ].y));
      dVl = 0.5*((Vel[k-1][j][i  ].y + Vel[k][j][i  ].y) -
                 (Vel[k-1][j][i-1].y + Vel[k][j][i-1].y));
      dVc = dVr + dVl;

      dVydx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVz/dx */
      dVr = 0.5*((Vel[k-1][j][i+1].z + Vel[k][j][i+1].z) -
                 (Vel[k-1][j][i  ].z + Vel[k][j][i  ].z));
      dVl = 0.5*((Vel[k-1][j][i  ].z + Vel[k][j][i  ].z) -
                 (Vel[k-1][j][i-1].z + Vel[k][j][i-1].z));
      dVc = dVr + dVl;

      dVzdx = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdx = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx1;
      }

/* Monotonized Velocity gradient dVx/dy */
      dVr = 0.5*((Vel[k-1][j+1][i].x + Vel[k][j+1][i].x) -
                 (Vel[k-1][j  ][i].x + Vel[k][j  ][i].x));
      dVl = 0.5*((Vel[k-1][j  ][i].x + Vel[k][j  ][i].x) -
                 (Vel[k-1][j-1][i].x + Vel[k][j-1][i].x));
      dVc = dVr + dVl;

      dVxdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVxdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVy/dy */
      dVr = 0.5*((Vel[k-1][j+1][i].y + Vel[k][j+1][i].y) -
                 (Vel[k-1][j  ][i].y + Vel[k][j  ][i].y));
      dVl = 0.5*((Vel[k-1][j  ][i].y + Vel[k][j  ][i].y) -
                 (Vel[k-1][j-1][i].y + Vel[k][j-1][i].y));
      dVc = dVr + dVl;

      dVydy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVydy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* Monotonized Velocity gradient dVz/dy */
      dVr = 0.5*((Vel[k-1][j+1][i].z + Vel[k][j+1][i].z) -
                 (Vel[k-1][j  ][i].z + Vel[k][j  ][i].z));
      dVl = 0.5*((Vel[k-1][j  ][i].z + Vel[k][j  ][i].z) -
                 (Vel[k-1][j-1][i].z + Vel[k][j-1][i].z));
      dVc = dVr + dVl;

      dVzdy = 0.0;
      if (dVl*dVr > 0.0) {
        lim_slope = MIN(fabs(dVl),fabs(dVr));
        dVzdy = SIGN(dVc)*MIN(0.5*fabs(dVc),2.0*lim_slope)/pG->dx2;
      }

/* compute BBdV and div(V) */

      BBdV =
        Bx*(Bx*dVxdx + By*dVxdy + Bz*(Vel[k][j][i].x-Vel[k-1][j][i].x)/pG->dx3)+
        By*(Bx*dVydx + By*dVydy + Bz*(Vel[k][j][i].y-Vel[k-1][j][i].y)/pG->dx3)+
        Bz*(Bx*dVzdx + By*dVzdy + Bz*(Vel[k][j][i].z-Vel[k-1][j][i].z)/pG->dx3);
      BBdV /= B02;

      divV = dVxdx + dVydy + (Vel[k][j][i].z-Vel[k-1][j][i].z)/pG->dx3;
      nud = nu_V*0.5*(pG->U[k][j][i].d + pG->U[k-1][j][i].d);

      qa = nud*(BBdV - ONE_3RD*divV);

/* Compute fluxes */

        x3Flux[k][j][i].Mx = qa*(3.0*Bx*Bz/B02);
        x3Flux[k][j][i].My = qa*(3.0*By*Bz/B02);
        x3Flux[k][j][i].Mz = qa*(3.0*Bz*Bz/B02 - 1.0);

#ifndef BAROTROPIC
        x3Flux[k][j][i].E  =
           0.5*(Vel[k-1][j][i].x + Vel[k][j][i].x)*x3Flux[k][j][i].Mx +
           0.5*(Vel[k-1][j][i].y + Vel[k][j][i].y)*x3Flux[k][j][i].My +
           0.5*(Vel[k-1][j][i].z + Vel[k][j][i].z)*x3Flux[k][j][i].Mz;
#endif /* BAROTROPIC */
      }
    }
  }

/*--- Step 5a ------------------------------------------------------------------
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

/*--- Step 5b ------------------------------------------------------------------
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

/*--- Step 5c ------------------------------------------------------------------
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
#endif /* BRAGINSKII */

  return;
}

/*----------------------------------------------------------------------------*/
/* brag_viscosity_init: Allocate temporary arrays
 */

void brag_viscosity_init(int nx1, int nx2, int nx3)
{
#ifdef BRAGINSKII
  int Nx1 = nx1 + 2*nghost, Nx2, Nx3;
  if (nx2 > 1){
    Nx2 = nx2 + 2*nghost;
  } else {
    ath_error("[brag_viscosity_init]:Braginskii viscosity requires 2D or 3D\n");
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
  if ((Vel = (ThreeDVect***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(ThreeDVect)))
    == NULL) goto on_error;
  return;

  on_error:
  brag_viscosity_destruct();
  ath_error("[brag_viscosity_init]: malloc returned a NULL pointer\n");
#endif /* BRAGINSKII */
  return;
}

/*----------------------------------------------------------------------------*/
/* brag_viscosity_destruct: Free temporary arrays
 */

void brag_viscosity_destruct(void)
{
#ifdef BRAGINSKII
  if (x1Flux != NULL) free_3d_array(x1Flux);
  if (x2Flux != NULL) free_3d_array(x2Flux);
  if (x3Flux != NULL) free_3d_array(x3Flux);
  if (Vel != NULL) free_3d_array(Vel);
#endif /* BRAGINSKII */

  return;
}
