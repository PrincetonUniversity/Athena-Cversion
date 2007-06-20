#include "copyright.h"
/*==============================================================================
 * FILE: selfg_fft.c
 *
 * PURPOSE: Contains functions to solve Poisson's equation for self-gravity in
 *   1D, 2D and 3D using FFTs (actually, the 1D algorithm uses Forward 
 *   Elimination followed by Back Substitution: FEBS).
 *
 *   These functions require PERIODIC BCs and use the Jeans swindle.
 *
 *   The 2D and 3D f'ns use FFTW3.x, and for MPI parallel use Steve Plimpton's
 *   block decomposition routines added by N. Lemaster to /athena/fftsrc.
 *   This means to use these fns the code must be
 *      (1) configured with --enable-fft
 *      (2) compiled with links to FFTW libraries
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_by_fft_1d() - actually uses FEBS
 *   selfg_by_fft_2d() - 2D Poisson solver using FFTs
 *   selfg_by_fft_3d() - 3D Poisson solver using FFTs
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef FFT_ENABLED

#include "ath_fft.h"

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_2d_fft_plan *fplan2d, *bplan2d;
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL;

#endif /* FFT_ENABLED */

/*----------------------------------------------------------------------------*/
/* selfg_by_fft_1d:  Actually uses forward elimination - back substituion!!
 *   Only works for uniform grid, periodic boundary conditions 
 *   This algorithm taken from pp.35-38 of Hockney & Eastwood
 */

void selfg_by_fft_1d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY_USING_FFT
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real total_Phi=0.0,drho,dx_sq = (pG->dx1*pG->dx1);

/* Copy current potential into old */

  for (i=is-nghost; i<=ie+nghost; i++){
    pG->Phi_old[ks][js][i] = pG->Phi[ks][js][i];
  }

/* Compute new potential */

  pG->Phi[ks][js][is] = 0.0;
  for (i=is; i<=ie; i++) {
    drho = (pG->U[ks][js][i].d - grav_mean_rho);
    pG->Phi[ks][js][is] += ((float)(i-is+1))*four_pi_G*dx_sq*drho;
  }
  pG->Phi[ks][js][is] /= (float)(pG->Nx1);

  drho = (pG->U[ks][js][is].d - grav_mean_rho);
  pG->Phi[ks][js][is+1] = 2.0*pG->Phi[ks][js][is] + four_pi_G*dx_sq*drho;
  for (i=is+2; i<=ie; i++) {
    drho = (pG->U[ks][js][i-1].d - grav_mean_rho);
    pG->Phi[ks][js][i] = four_pi_G*dx_sq*drho 
      + 2.0*pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i-2];
  }

/* Normalize so mean Phi is zero */

  for (i=is; i<=ie; i++) {
    total_Phi += pG->Phi[ks][js][i];
  }
  total_Phi /= (float)(pG->Nx1);

  for (i=is; i<=ie; i++) {
    pG->Phi[ks][js][i] -= total_Phi;
  }


/* Apply periodic boundary conditions */

  for (i=1; i<=nghost; i++) {
    pG->Phi[ks][js][is-i] =  pG->Phi[ks][js][ie-(i-1)];
    pG->Phi[ks][js][ie+i] =  pG->Phi[ks][js][is+(i-1)];
  }

#endif /* SELF_GRAVITY_USING_FFT */
}

#ifdef FFT_ENABLED
/*----------------------------------------------------------------------------*/
/* selfg_by_fft_2d:
 *   Only works for uniform grid, periodic boundary conditions, and dx1=dx2
 */

void selfg_by_fft_2d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY_USING_FFT
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dkx,dky;

/* Copy current potential into old */

  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[ks][j][i] = pG->Phi[ks][j][i];
    }
  }

/* Forward FFT of 4\piG*(d-d0) */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][0] =
        four_pi_G*(pG->U[ks][j][i].d - grav_mean_rho);
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][1] = 0.0;
    }
  }

  ath_2d_fft(fplan2d, work);

/* Compute potential in Fourier space.  Multiple loops are used to avoid divide
 * by zero at i=is,j=js   */

  dkx = 2.0*PI/(double)(pD->ixe - pD->ixs + 1);
  dky = 2.0*PI/(double)(pD->jxe - pD->jxs + 1);

  work[F2DI(0,0,pG->Nx1,pG->Nx2)][0] = 0.0;
  work[F2DI(0,0,pG->Nx1,pG->Nx2)][1] = 0.0;

  for (j=js+1; j<=je; j++){
    work[F2DI(0,j-js,pG->Nx1,pG->Nx2)][0] *= SQR(pG->dx1)/
     (2.0*cos((j-js)*dky) - 2.0);
    work[F2DI(0,j-js,pG->Nx1,pG->Nx2)][1] *= SQR(pG->dx1)/
     (2.0*cos((j-js)*dky) - 2.0);
  }

  for (i=is+1; i<=ie; i++){
    for (j=js; j<=je; j++){
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][0] *= SQR(pG->dx1)/
       (2.0*cos((i-is)*dkx) + 2.0*cos((j-js)*dky) - 4.0);
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][1] *= SQR(pG->dx1)/
       (2.0*cos((i-is)*dkx) + 2.0*cos((j-js)*dky) - 4.0);
    }
  }

/* Backward FFT and set potential in real space */

  ath_2d_fft(bplan2d, work);

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[ks][j][i] = work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][0]/
        bplan2d->cnt;
    }
  }

/* Set periodic boundary conditions */

  for (j=js; j<=je; j++){
    for (i=1; i<=nghost; i++) {
      pG->Phi[ks][j][is-i] = pG->Phi[ks][j][ie-(i-1)];
      pG->Phi[ks][j][ie+i] = pG->Phi[ks][j][is+(i-1)];
    }
  }

  for (j=1; j<=nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pG->Phi[ks][js-j][i] = pG->Phi[ks][je-(j-1)][i];
      pG->Phi[ks][je+j][i] = pG->Phi[ks][js+(j-1)][i];
    }
  }

#endif /* SELF_GRAVITY_USING_FFT */
  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_fft_3d:
 *   Only works for uniform grid, periodic boundary conditions, and dx1=dx2=dx3
 */

void selfg_by_fft_3d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY_USING_FFT
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dkx,dky,dkz;

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
    }
  }}

/* Forward FFT of 4\piG*(d-d0) */

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] = 
        four_pi_G*(pG->U[k][j][i].d - grav_mean_rho);
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] = 0.0;
    }
  }}

  ath_3d_fft(fplan3d, work);

/* Compute potential in Fourier space.  Multiple loops are used to avoid divide
 * by zero at i=is,j=js   */

  dkx = 2.0*PI/(double)(pD->ixe - pD->ixs + 1);
  dky = 2.0*PI/(double)(pD->jxe - pD->jxs + 1);
  dkz = 2.0*PI/(double)(pD->kxe - pD->kxs + 1);

  work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][0] = 0.0;
  work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][1] = 0.0;

  for (k=ks+1; k<=ke; k++){
    work[F3DI(0,0,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= SQR(pG->dx1)/
     (2.0*cos((k-ks)*dkz) - 2.0);
    work[F3DI(0,0,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= SQR(pG->dx1)/
     (2.0*cos((k-ks)*dkz) - 2.0);
  }

  for (j=js+1; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(0,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= SQR(pG->dx1)/
       (2.0*cos((j-js)*dky) + 2.0*cos((k-ks)*dkz) - 4.0);
      work[F3DI(0,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= SQR(pG->dx1)/
       (2.0*cos((j-js)*dky) + 2.0*cos((k-ks)*dkz) - 4.0);
    }
  }

  for (i=is+1; i<=ie; i++){
  for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= SQR(pG->dx1)/
       (2.0*cos((i-is)*dkx) + 2.0*cos((j-js)*dky) + 2.0*cos((k-ks)*dkz) - 6.0);
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= SQR(pG->dx1)/
       (2.0*cos((i-is)*dkx) + 2.0*cos((j-js)*dky) + 2.0*cos((k-ks)*dkz) - 6.0);
    }
  }}

/* Backward FFT and set potential in real space */

  ath_3d_fft(bplan3d, work);

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] = work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0]/
        bplan3d->cnt;
    }
  }}

/* Set periodic boundary conditions */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=1; i<=nghost; i++) {
        pG->Phi[k][j][is-i] = pG->Phi[k][j][ie-(i-1)];
        pG->Phi[k][j][ie+i] = pG->Phi[k][j][is+(i-1)];
      }
    }
  }

  for (k=ks; k<=ke; k++){
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->Phi[k][js-j][i] = pG->Phi[k][je-(j-1)][i];
        pG->Phi[k][je+j][i] = pG->Phi[k][js+(j-1)][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pG->Phi[ks-k][j][i] = pG->Phi[ke-(k-1)][j][i];
        pG->Phi[ke+k][j][i] = pG->Phi[ks+(k-1)][j][i];
      }
    }
  }

#endif /* SELF_GRAVITY_USING_FFT */
  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_fft_2d_init:
 *   Initializes plans for forward/backward FFTs, and allocates memory needed
 *   by FFTW.
 */

void selfg_by_fft_2d_init(Grid *pG, Domain *pD)
{
  fplan2d = ath_2d_fft_quick_plan(pG, pD, NULL, ATH_FFT_FORWARD);
  bplan2d = ath_2d_fft_quick_plan(pG, pD, NULL, ATH_FFT_BACKWARD);
  work = ath_2d_fft_malloc(fplan2d);
}

/*----------------------------------------------------------------------------*/
/* selfg_by_fft_3d_init:
 *   Initializes plans for forward/backward FFTs, and allocates memory needed
 *   by FFTW.
 */

void selfg_by_fft_3d_init(Grid *pG, Domain *pD)
{
  fplan3d = ath_3d_fft_quick_plan(pG, pD, NULL, ATH_FFT_FORWARD);
  bplan3d = ath_3d_fft_quick_plan(pG, pD, NULL, ATH_FFT_BACKWARD);
  work = ath_3d_fft_malloc(fplan3d);
}

#endif /* FFT_ENABLED */
