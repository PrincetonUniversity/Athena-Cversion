#include "copyright.h"
/*==============================================================================
 * FILE: self_gravity.c
 *
 * PURPOSE: Contains functions to solve Poisson's equation for self-gravity in
 *   1D, 2D and 3D.  Contains algorithms based on both multigrid and FFTs
 *
 * HISTORY:
 *   may-2007 - Initial framework and 1D solver written by J. Stone
 *   may-2007 - Multigrid solvers written by Irene Balmes
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   self_gravity_init()  - sets pointer to appropriate self-gravity function
 *   selfg_by_multig_1d() - multigrid algorithm in 1D
 *   
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* selfg_flux_correction: Corrects source terms added to momentum and energy
 *   equations in integrator (which used old potential at t^{n}) with terms
 *   computed using the new potential (at t^{n+1}) to keep terms 2nd order.
 *   Requires knowing the mass fluxes at cell faces, stored in Grid.
 */

void selfg_flux_correction(Grid *pG)
{
#ifdef SELF_GRAVITY
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int dim=0;
  Real dtodx1 = pG->dt/pG->dx1;
  Real dtodx2 = pG->dt/pG->dx2;
  Real dtodx3 = pG->dt/pG->dx3;
  Real phic,phil,phir;
  Real gxl,gxr,gyl,gyr,gzl,gzr,frm1,flm1,frm2,flm2,frm3,flm3;
  
/* Calculate the dimensions  */
  if(pG->Nx1 > 1) dim++;
  if(pG->Nx2 > 1) dim++;
  if(pG->Nx3 > 1) dim++;


/* Since the form of the source terms depend on dimensions of problem, there are
 * different loops for 1D, 2D, and 3D 
 */

  switch(dim){
/*------------------------- 1D problem ---------------------------------------*/
  case 1:
    for (i=is; i<=ie; i++){

/*  Subtract 1/2 the source term at the old time level (Phi^{n}) */

      phic = pG->Phi_old[ks][js][i];
      phil = 0.5*(pG->Phi_old[ks][js][i-1] + pG->Phi_old[ks][js][i  ]);
      phir = 0.5*(pG->Phi_old[ks][js][i  ] + pG->Phi_old[ks][js][i+1]);

      gxl = (pG->Phi_old[ks][js][i-1] - pG->Phi_old[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi_old[ks][js][i  ] - pG->Phi_old[ks][js][i+1])/(pG->dx1);

/* 1-momentum flux */
/*
      flm1 = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      frm1 = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;
*/
      flm1 = 0.5*(gxl*gxl)/four_pi_G;
      frm1 = 0.5*(gxr*gxr)/four_pi_G;

      pG->U[ks][js][i].M1 -= 0.5*dtodx1*(frm1-flm1);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E -=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(phil - phic) +
                     pG->x1MassFlux[ks][js][i+1]*(phic - phir));
#endif

/*  Add 1/2 the source term at the new time level (Phi^{n+1}) */

      phic = pG->Phi[ks][js][i];
      phil = 0.5*(pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ]);
      phir = 0.5*(pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1]);

      gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* 1-momentum flux */
/*
      flm1 = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      frm1 = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;
*/
      flm1 = 0.5*(gxl*gxl)/four_pi_G;
      frm1 = 0.5*(gxr*gxr)/four_pi_G;

      pG->U[ks][js][i].M1 += 0.5*dtodx1*(frm1-flm1);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E +=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(phil - phic) +
                     pG->x1MassFlux[ks][js][i+1]*(phic - phir));
#endif
    }
    break;

/*------------------------- 2D problem ---------------------------------------*/
  case 2:
    for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){

/*  Subtract 1/2 the source term at the old time level (Phi^{n}) */

      phic = pG->Phi_old[ks][js][i];
      phil = 0.5*(pG->Phi_old[ks][js][i-1] - pG->Phi_old[ks][js][i  ]);
      phir = 0.5*(pG->Phi_old[ks][js][i  ] - pG->Phi_old[ks][js][i+1]);

      gxl = (pG->Phi_old[k][j][i-1] - pG->Phi_old[k][j][i  ])/(pG->dx1);
      gxr = (pG->Phi_old[k][j][i  ] - pG->Phi_old[k][j][i+1])/(pG->dx1);

/* 1-momentum flux */
      flm1 = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      frm1 = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

      pG->U[ks][js][i].M1 += 0.5*dtodx1*(frm1-flm1);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E -=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(phil - phic) +
                     pG->x1MassFlux[ks][js][i+1]*(phic - phir));
#endif

/*  Add 1/2 the source term at the new time level (Phi^{n+1}) */

      phic = pG->Phi[ks][js][i];
      phil = 0.5*(pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ]);
      phir = 0.5*(pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1]);

      gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* 1-momentum flux */
      flm1 = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      frm1 = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

      pG->U[ks][js][i].M1 -= 0.5*dtodx1*(frm1-flm1);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E +=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(phil - phic) +
                     pG->x1MassFlux[ks][js][i+1]*(phic - phir));
#endif
    }}
    break;

/*------------------------- 3D problem ---------------------------------------*/
  case 3:
    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
/* Calculate an effective phi at the left and right x-interface */
        phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j][i  ]);
        phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+1]);

        gxl = (pG->Phi[k][j][i-1] - pG->Phi[k][j][i  ])/(pG->dx1);
        gxr = (pG->Phi[k][j][i  ] - pG->Phi[k][j][i+1])/(pG->dx1);

        gyl = 0.25*((pG->Phi[k][j-1][i-1] - pG->Phi[k][j+1][i-1]) +
                    (pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]) )/(pG->dx2);
        gyr = 0.25*((pG->Phi[k][j-1][i ] - pG->Phi[k][j+1][i  ]) +
                    (pG->Phi[k][j-1][i+1] - pG->Phi[k][j+1][i+1]) )/(pG->dx2);

        gzl = 0.25*((pG->Phi[k-1][j][i-1] - pG->Phi[k+1][j][i-1]) +
                    (pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) )/(pG->dx3);
        gzr = 0.25*((pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) +
                    (pG->Phi[k-1][j][i+1] - pG->Phi[k+1][j][i+1]) )/(pG->dx3);

/* 1-momentum flux */
/*
        flm1 = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        frm1 = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;
*/
        flm1 = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G;
        frm1 = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G;

/* 2-momentum flux */
        flm2 = gxl*gyl/four_pi_G;
        frm2 = gxr*gyr/four_pi_G;

/* 3-momentum flux */
        flm3 = gxl*gzl/four_pi_G;
        frm3 = gxr*gzr/four_pi_G;

/* Update the momenta */
        pG->U[k][j][i].M1 += dtodx1*(flm1 - frm1);
        pG->U[k][j][i].M2 += dtodx1*(flm2 - frm2);
        pG->U[k][j][i].M3 += dtodx3*(flm3 - frm3);

/* Update the total energy, if needed */
#ifdef ADIABATIC
        pG->U[k][j][i].E += 0.5*dtodx1*
          (pG->x1MassFlux[k][j][i-1]*(pG->Phi[k][j][i-1] - pG->Phi[k][j][i  ]) +
           pG->x1MassFlux[k][j][i  ]*(pG->Phi[k][j][i  ] - pG->Phi[k][j][i+1]));
#endif /* ADIABATIC */
      }
    }}
    break;

  } /* end of switch statement */
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_FEBS_1d:  Uses forward elimination - back substituion.
 *   Only works for uniform grid, periodic boundary conditions 
 *   This algorithm taken from pp.35-38 of Hockney & Eastwood
*/

void selfg_by_FEBS_1d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real drho;

/* Copy current potential into old */

  for (i=is-nghost; i<=ie+nghost; i++){
    pG->Phi_old[ks][js][i] = pG->Phi[ks][js][i];
  }

/* Compute new potential */

  pG->Phi[ks][js][is] = 0.0;
  for (i=is; i<=ie; i++) {
    drho = pG->U[ks][js][i].d - grav_mean_rho;
    pG->Phi[ks][js][is] += (float)(i-is+1)*four_pi_G*drho;
  }
  pG->Phi[ks][js][is] /= (float)(pG->Nx1);

  drho = pG->U[ks][js][is].d - grav_mean_rho;
  pG->Phi[ks][js][is+1] = 2.0*pG->Phi[ks][js][is] + four_pi_G*drho;
  for (i=is+2; i<=ie; i++) {
    drho = pG->U[ks][js][i-1].d - grav_mean_rho;
    pG->Phi[ks][js][i] = four_pi_G*drho 
      + 2.0*pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i-2];
  }

/* Apply periodic boundary conditions */

  for (i=1; i<=nghost; i++) {
    pG->Phi[ks][js][is-i] =  pG->Phi[ks][js][ie-(i-1)];
    pG->Phi[ks][js][ie+i] =  pG->Phi[ks][js][is+(i-1)];
  }

#endif
}

void selfg_by_multig_2d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;

/* Copy current potential into old */

  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[ks][j][i] = pG->Phi[ks][j][i];
    }
  }

#endif
  return;
}

void selfg_by_multig_3d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* self_gravity_init: initialize pointer to appropriate self-gravity f'n
 *   VGDFun_t is a function of type void which takes a Grid and a Domain as 
 *   arguments.
 */

VGDFun_t self_gravity_init(int Nx1, int Nx2, int Nx3)
{
  int dim = 0;

/* Calculate the dimensions  */
  if(Nx1 > 1) dim++;
  if(Nx2 > 1) dim++;
  if(Nx3 > 1) dim++;

/* set function pointer based on dimensions and algorithm */
  switch(dim){
#ifdef SELF_GRAVITY_USING_MULTIGRID
  case 1:
    if(Nx1 <= 1) break;
    return selfg_by_FEBS_1d;
  case 2:
    if(Nx3 > 1) break;
    return selfg_by_multig_2d;
  case 3:
    return selfg_by_multig_3d;
#endif
#ifdef SELF_GRAVITY_USING_FFT
  case 1:
    if(Nx1 <= 1) break;
    return selfg_by_fft_1d;
  case 2:
    if(Nx3 > 1) break;
    return selfg_by_fft_2d;
  case 3:
    return selfg_by_fft_3d;
#endif
  }

/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}
