#include "copyright.h"
/*==============================================================================
 * FILE: selfg.c
 *
 * PURPOSE: Contains functions to control solution of Poisson's equation for
 *   self-gravity.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_flux_correction() - 2nd order corrections for self-gravity terms
 *   selfg_init()            - sets pointer to appropriate self-gravity function
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real ***dPhi=NULL;

/*----------------------------------------------------------------------------*/
/* selfg_flux_correction: Adds second-order correction to source terms for the
 *   momentum and energy equations.  This requires subtracting 1/2 the source
 *   terms computed with the old potential, and adding 1/2 the source terms
 *   computed with the new potential.
 *
 *   The source terms for momentum are computed using the divergence of the
 *   gravitational stress tensor to conserve momentum exactly.
 *     dM/dt = -Div(G);   G = (gg - 0.5g^2)/4\piG;   g=-Grad(Phi);
 *
 *   The source terms for the energy are added using the mass fluxes at cell
 *   faces, to improve conservation.
 *     S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
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
  Real dphic,dphil,dphir,gxl,gxr,gyl,gyr,gzl,gzr;
  Real flx_m1r,flx_m1l,flx_m2r,flx_m2l,flx_m3r,flx_m3l;
  
/* Calculate the dimensions  */
  if(pG->Nx1 > 1) dim++;
  if(pG->Nx2 > 1) dim++;
  if(pG->Nx3 > 1) dim++;


/* The divergence of the gravitational stress tensor depends on the dimensions
 * of the problem.
 */

  switch(dim){
/*------------------------- 1D problem ---------------------------------------*/
  case 1:
    for (i=is-1; i<=ie+1; i++){
      dPhi[ks][js][i] = pG->Phi[ks][js][i] - pG->Phi_old[ks][js][i];
    }

/* Step 1 for 1D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (i=is; i<=ie; i++){
      dphic = dPhi[ks][js][i];
      dphil = 0.5*(dPhi[ks][js][i-1] + dPhi[ks][js][i  ]);
      dphir = 0.5*(dPhi[ks][js][i  ] + dPhi[ks][js][i+1]);

/* gx centered at L and R x1-faces */
      gxl = (dPhi[ks][js][i-1] - dPhi[ks][js][i  ])/(pG->dx1);
      gxr = (dPhi[ks][js][i  ] - dPhi[ks][js][i+1])/(pG->dx1);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
      flx_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*dphil;
      flx_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*dphir;

/* Update momenta and energy with d/dx1 terms  */
      pG->U[ks][js][i].M1 -= 0.5*dtodx1*(flx_m1r-flx_m1l);
#ifndef ISOTHERMAL
      pG->U[ks][js][i].E -=
         0.5*dtodx1*(pG->x1MassFlux[ks][js][i  ]*(dphic - dphil) +
                     pG->x1MassFlux[ks][js][i+1]*(dphir - dphic));
#endif
    }
    break;

/*------------------------- 2D problem ---------------------------------------*/
  case 2:
    for (j=js-1; j<=je+1; j++){
      for (i=is-1; i<=ie+1; i++){
        dPhi[ks][j][i] = pG->Phi[ks][j][i] - pG->Phi_old[ks][j][i];
      }
    }

/* Step 1 for 2D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        dphic = dPhi[ks][j][i];
        dphil = 0.5*(dPhi[ks][j][i-1] + dPhi[ks][j][i  ]);
        dphir = 0.5*(dPhi[ks][j][i  ] + dPhi[ks][j][i+1]);

/* gx and gy centered at L and R x1-faces */
        gxl = (dPhi[ks][j][i-1] - dPhi[ks][j][i  ])/(pG->dx1);
        gxr = (dPhi[ks][j][i  ] - dPhi[ks][j][i+1])/(pG->dx1);

        gyl = 0.25*((dPhi[ks][j-1][i-1] - dPhi[ks][j+1][i-1]) +
                    (dPhi[ks][j-1][i  ] - dPhi[ks][j+1][i  ]) )/(pG->dx2);
        gyr = 0.25*((dPhi[ks][j-1][i  ] - dPhi[ks][j+1][i  ]) +
                    (dPhi[ks][j-1][i+1] - dPhi[ks][j+1][i+1]) )/(pG->dx2);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
        flx_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*dphil;
        flx_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*dphir;

        flx_m2l = gxl*gyl/four_pi_G;
        flx_m2r = gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
        pG->U[ks][j][i].M1 -= 0.5*dtodx1*(flx_m1r - flx_m1l);
        pG->U[ks][j][i].M2 -= 0.5*dtodx1*(flx_m2r - flx_m2l);
#ifndef ISOTHERMAL
        pG->U[ks][j][i].E -=
           0.5*dtodx1*(pG->x1MassFlux[ks][j][i  ]*(dphic - dphil) +
                       pG->x1MassFlux[ks][j][i+1]*(dphir - dphic));
#endif
      }
    }

/* Step 2 for 2D.  Add fluxes and source terms due to (d/dx2) terms  */

    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        dphic = dPhi[ks][j][i];
        dphil = 0.5*(dPhi[ks][j-1][i] + dPhi[ks][j  ][i]);
        dphir = 0.5*(dPhi[ks][j  ][i] + dPhi[ks][j+1][i]);

/* gx and gy centered at L and R x2-faces */
        gxl = 0.25*((dPhi[ks][j-1][i-1] - dPhi[ks][j-1][i+1]) +
                    (dPhi[ks][j  ][i-1] - dPhi[ks][j  ][i+1]) )/(pG->dx1);
        gxr = 0.25*((dPhi[ks][j  ][i-1] - dPhi[ks][j  ][i+1]) +
                    (dPhi[ks][j+1][i-1] - dPhi[ks][j+1][i+1]) )/(pG->dx1);

        gyl = (dPhi[ks][j-1][i] - dPhi[ks][j  ][i])/(pG->dx2);
        gyr = (dPhi[ks][j  ][i] - dPhi[ks][j+1][i])/(pG->dx2);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
        flx_m1l = gyl*gxl/four_pi_G;
        flx_m1r = gyr*gxr/four_pi_G;

        flx_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*dphil;
        flx_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*dphir;

/* Update momenta and energy with d/dx2 terms  */
        pG->U[ks][j][i].M1 -= 0.5*dtodx2*(flx_m1r - flx_m1l);
        pG->U[ks][j][i].M2 -= 0.5*dtodx2*(flx_m2r - flx_m2l);
#ifndef ISOTHERMAL
        pG->U[ks][j][i].E -=
           0.5*dtodx2*(pG->x2MassFlux[ks][j  ][i]*(dphic - dphil) +
                       pG->x2MassFlux[ks][j+1][i]*(dphir - dphic));
#endif
      }
    }

    break;

/*------------------------- 3D problem ---------------------------------------*/
  case 3:
    for (k=ks-1; k<=ke+1; k++){
    for (j=js-1; j<=je+1; j++){
      for (i=is-1; i<=ie+1; i++){
        dPhi[k][j][i] = pG->Phi[k][j][i] - pG->Phi_old[k][j][i];
      }
    }}

/* Step 1 for 3D.  Add fluxes and source terms due to (d/dx1) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        dphic = dPhi[k][j][i];
        dphil = 0.5*(dPhi[k][j][i-1] + dPhi[k][j][i  ]);
        dphir = 0.5*(dPhi[k][j][i  ] + dPhi[k][j][i+1]);

/* gx, gy and gz centered at L and R x1-faces */
        gxl = (dPhi[k][j][i-1] - dPhi[k][j][i  ])/(pG->dx1);
        gxr = (dPhi[k][j][i  ] - dPhi[k][j][i+1])/(pG->dx1);

        gyl = 0.25*((dPhi[k][j-1][i-1] - dPhi[k][j+1][i-1]) +
                    (dPhi[k][j-1][i  ] - dPhi[k][j+1][i  ]))/(pG->dx2);
        gyr = 0.25*((dPhi[k][j-1][i  ] - dPhi[k][j+1][i  ]) +
                    (dPhi[k][j-1][i+1] - dPhi[k][j+1][i+1]))/(pG->dx2);

        gzl = 0.25*((dPhi[k-1][j][i-1] - dPhi[k+1][j][i-1]) +
                    (dPhi[k-1][j][i  ] - dPhi[k+1][j][i  ]))/(pG->dx3);
        gzr = 0.25*((dPhi[k-1][j][i  ] - dPhi[k+1][j][i  ]) +
                    (dPhi[k-1][j][i+1] - dPhi[k+1][j][i+1]))/(pG->dx3);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
        flx_m1l = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G + grav_mean_rho*dphil;
        flx_m1r = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G + grav_mean_rho*dphir;

        flx_m2l = gxl*gyl/four_pi_G;
        flx_m2r = gxr*gyr/four_pi_G;

        flx_m3l = gxl*gzl/four_pi_G;
        flx_m3r = gxr*gzr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx1*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx1*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx1*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx1*
          (pG->x1MassFlux[k][j][i  ]*(dphic - dphil) +
           pG->x1MassFlux[k][j][i+1]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

/* Step 2 for 3D.  Add fluxes and source terms due to (d/dx2) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        dphic = dPhi[k][j][i];
        dphil = 0.5*(dPhi[k][j-1][i] + dPhi[k][j  ][i]);
        dphir = 0.5*(dPhi[k][j  ][i] + dPhi[k][j+1][i]);

/* gx, gy and gz centered at L and R x2-faces */
        gxl = 0.25*((dPhi[k][j-1][i-1] - dPhi[k][j-1][i+1]) +
                    (dPhi[k][j  ][i-1] - dPhi[k][j  ][i+1]))/(pG->dx1);
        gxr = 0.25*((dPhi[k][j  ][i-1] - dPhi[k][j  ][i+1]) +
                    (dPhi[k][j+1][i-1] - dPhi[k][j+1][i+1]))/(pG->dx1);

        gyl = (dPhi[k][j-1][i] - dPhi[k][j  ][i])/(pG->dx2);
        gyr = (dPhi[k][j  ][i] - dPhi[k][j+1][i])/(pG->dx2);

        gzl = 0.25*((dPhi[k-1][j-1][i] - dPhi[k+1][j-1][i]) +
                    (dPhi[k-1][j  ][i] - dPhi[k+1][j  ][i]))/(pG->dx3);
        gzr = 0.25*((dPhi[k-1][j  ][i] - dPhi[k+1][j  ][i]) +
                    (dPhi[k-1][j+1][i] - dPhi[k+1][j+1][i]))/(pG->dx3);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
        flx_m1l = gyl*gxl/four_pi_G;
        flx_m1r = gyr*gxr/four_pi_G;

        flx_m2l = 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G + grav_mean_rho*dphil;
        flx_m2r = 0.5*(gyr*gyr-gxr*gxr-gzr*gzr)/four_pi_G + grav_mean_rho*dphir;

        flx_m3l = gyl*gzl/four_pi_G;
        flx_m3r = gyr*gzr/four_pi_G;

/* Update momenta and energy with d/dx2 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx2*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx2*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx2*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx2*
          (pG->x2MassFlux[k][j  ][i]*(dphic - dphil) +
           pG->x2MassFlux[k][j+1][i]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

/* Step 3 for 3D.  Add fluxes and source terms due to (d/dx3) terms  */

    for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        dphic = dPhi[k][j][i];
        dphil = 0.5*(dPhi[k-1][j][i] + dPhi[k  ][j][i]);
        dphir = 0.5*(dPhi[k  ][j][i] + dPhi[k+1][j][i]);

/* gx, gy and gz centered at L and R x3-faces */
        gxl = 0.25*((dPhi[k-1][j][i-1] - dPhi[k-1][j][i+1]) +
                    (dPhi[k  ][j][i-1] - dPhi[k  ][j][i+1]) )/(pG->dx1);
        gxr = 0.25*((dPhi[k  ][j][i-1] - dPhi[k  ][j][i+1]) +
                    (dPhi[k+1][j][i-1] - dPhi[k+1][j][i+1]) )/(pG->dx1);

        gyl = 0.25*((dPhi[k-1][j-1][i] - dPhi[k-1][j+1][i]) +
                    (dPhi[k  ][j-1][i] - dPhi[k  ][j+1][i]) )/(pG->dx2);
        gyr = 0.25*((dPhi[k  ][j-1][i] - dPhi[k  ][j+1][i]) +
                    (dPhi[k+1][j-1][i] - dPhi[k+1][j+1][i]) )/(pG->dx2);

        gzl = (dPhi[k-1][j][i] - dPhi[k  ][j][i])/(pG->dx3);
        gzr = (dPhi[k  ][j][i] - dPhi[k+1][j][i])/(pG->dx3);

/* momentum fluxes in x3.  2nd term is needed only if Jean's swindle used */
        flx_m1l = gzl*gxl/four_pi_G;
        flx_m1r = gzr*gxr/four_pi_G;

        flx_m2l = gzl*gyl/four_pi_G;
        flx_m2r = gzr*gyr/four_pi_G;

        flx_m3l = 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*dphil;
        flx_m3r = 0.5*(gzr*gzr-gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*dphir;

/* Update momenta and energy with d/dx3 terms  */
        pG->U[k][j][i].M1 -= 0.5*dtodx3*(flx_m1r - flx_m1l);
        pG->U[k][j][i].M2 -= 0.5*dtodx3*(flx_m2r - flx_m2l);
        pG->U[k][j][i].M3 -= 0.5*dtodx3*(flx_m3r - flx_m3l);
#ifdef ADIABATIC
        pG->U[k][j][i].E -= 0.5*dtodx3*
          (pG->x3MassFlux[k  ][j][i]*(dphic - dphil) +
           pG->x3MassFlux[k+1][j][i]*(dphir - dphic));
#endif /* ADIABATIC */
      }
    }}

    break;

  } /* end of switch statement */
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_init: initialize pointer to appropriate self-gravity f'n, 
 *   allocates memory for dPhi array used for flux correction.
 *
 *   VGDFun_t is a function of type void which takes a Grid and a Domain as 
 *   arguments.
 */

VGDFun_t selfg_init(Grid *pG, Domain *pD)
{
  int nx1,nx2,nx3,dim = 0;

  nx1 = pG->Nx1 + 2*nghost;
  nx2 = pG->Nx2 + 2*nghost;
  nx3 = pG->Nx3 + 2*nghost;
  if ((dPhi = (Real***)calloc_3d_array(nx3,nx2,nx1,sizeof(Real))) == NULL)
    ath_error("[self_gravity_init]: malloc returned a NULL pointer\n");

/* Calculate the dimensions  */
  if(pG->Nx1 > 1) dim++;
  if(pG->Nx2 > 1) dim++;
  if(pG->Nx3 > 1) dim++;

/* Return function pointer based on dimensions and algorithm */

  switch(dim){
#ifdef SELF_GRAVITY_USING_MULTIGRID
  case 1:
    return selfg_by_multig_1d;
  case 2:
    return selfg_by_multig_2d;
  case 3:
    return selfg_by_multig_3d;
#endif

/* for gravity using FFTs, also initialize plans and data for FFTW */
#ifdef SELF_GRAVITY_USING_FFT
  case 1:
    return selfg_by_fft_1d;
  case 2:
    selfg_by_fft_2d_init(pG, pD);
    return selfg_by_fft_2d;
  case 3:
    selfg_by_fft_3d_init(pG, pD);
    return selfg_by_fft_3d;
#endif
  }

  return NULL;
}
