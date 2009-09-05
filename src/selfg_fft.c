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
 *     (1) configured with --with-gravity=fft --enable-fft
 *     (2) compiled with links to FFTW libraries (may need to edit Makeoptions)
 *
 *   For NON-PERIODIC BCs, use selfg_multig() functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_by_fft_1d() - actually uses FEBS
 *   selfg_by_fft_2d() - 2D Poisson solver using FFTs
 *   selfg_by_fft_3d() - 3D Poisson solver using FFTs
 *   selfg_by_fft_2d_init() - initializes FFT plans for 2D
 *   selfg_by_fft_3d_init() - initializes FFT plans for 3D
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef FFT_ENABLED

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

#ifdef SELF_GRAVITY_USING_FFT
void selfg_by_fft_1d(Grid *pG, Domain *pD)
{
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
}
#endif /* SELF_GRAVITY_USING_FFT */

#ifdef FFT_ENABLED
#ifdef SELF_GRAVITY_USING_FFT
/*----------------------------------------------------------------------------*/
/* selfg_by_fft_2d:
 *   Only works for uniform grid, periodic boundary conditions
 */

void selfg_by_fft_2d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2);
  Real dkx,dky,pcoeff;

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
 * by zero at i=is,j=js, and to avoid if statement in loop   */
/* To compute kx,ky note that indices relative to whole Domain are needed */

  dkx = 2.0*PI/(double)(pD->ide - pD->ids + 1);
  dky = 2.0*PI/(double)(pD->jde - pD->jds + 1);

  if ((js+pG->jdisp)==0 && (is+pG->idisp)==0) {
    work[F2DI(0,0,pG->Nx1,pG->Nx2)][0] = 0.0;
    work[F2DI(0,0,pG->Nx1,pG->Nx2)][1] = 0.0;
  } else {
    pcoeff = 1.0/(((2.0*cos((is+pG->idisp)*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((js+pG->jdisp)*dky)-2.0)/dx2sq));
    work[F2DI(0,0,pG->Nx1,pG->Nx2)][0] *= pcoeff;
    work[F2DI(0,0,pG->Nx1,pG->Nx2)][1] *= pcoeff;
  }

  for (j=js+1; j<=je; j++){
    pcoeff = 1.0/(((2.0*cos((is+pG->idisp)*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((j +pG->jdisp)*dky)-2.0)/dx2sq));
    work[F2DI(0,j-js,pG->Nx1,pG->Nx2)][0] *= pcoeff;
    work[F2DI(0,j-js,pG->Nx1,pG->Nx2)][1] *= pcoeff;
  }

  for (i=is+1; i<=ie; i++){
    for (j=js; j<=je; j++){
      pcoeff = 1.0/(((2.0*cos((i+pG->idisp)*dkx)-2.0)/dx1sq) +
                    ((2.0*cos((j+pG->jdisp)*dky)-2.0)/dx2sq));
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][0] *= pcoeff;
      work[F2DI(i-is,j-js,pG->Nx1,pG->Nx2)][1] *= pcoeff;
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

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_fft_3d:
 *   Only works for uniform grid, periodic boundary conditions
 */

void selfg_by_fft_3d(Grid *pG, Domain *pD)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3);
  Real dkx,dky,dkz,pcoeff;

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
 * by zero at i=is,j=js,k=ks, and to avoid if statement in loop   */
/* To compute kx,ky,kz, note that indices relative to whole Domain are needed */

  dkx = 2.0*PI/(double)(pD->ide - pD->ids + 1);
  dky = 2.0*PI/(double)(pD->jde - pD->jds + 1);
  dkz = 2.0*PI/(double)(pD->kde - pD->kds + 1);

  if ((ks+pG->kdisp)==0 && (js+pG->jdisp)==0 && (is+pG->idisp)==0) {
    work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][0] = 0.0;
    work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][1] = 0.0;
  } else {
    pcoeff = 1.0/(((2.0*cos((is+pG->idisp)*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((js+pG->jdisp)*dky)-2.0)/dx2sq) +
                  ((2.0*cos((ks+pG->kdisp)*dkz)-2.0)/dx3sq));
    work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= pcoeff;
    work[F3DI(0,0,0,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= pcoeff;
  }


  for (k=ks+1; k<=ke; k++){
    pcoeff = 1.0/(((2.0*cos((is+pG->idisp)*dkx)-2.0)/dx1sq) +
                  ((2.0*cos((js+pG->jdisp)*dky)-2.0)/dx2sq) +
                  ((2.0*cos((k +pG->kdisp)*dkz)-2.0)/dx3sq));
    work[F3DI(0,0,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= pcoeff;
    work[F3DI(0,0,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= pcoeff;
  }

  for (j=js+1; j<=je; j++){
    for (k=ks; k<=ke; k++){
      pcoeff = 1.0/(((2.0*cos((is+pG->idisp)*dkx)-2.0)/dx1sq) +
                    ((2.0*cos((j +pG->jdisp)*dky)-2.0)/dx2sq) +
                    ((2.0*cos((k +pG->kdisp)*dkz)-2.0)/dx3sq));
      work[F3DI(0,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= pcoeff;
      work[F3DI(0,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= pcoeff;
    }
  }

  for (i=is+1; i<=ie; i++){
  for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      pcoeff = 1.0/(((2.0*cos((i+pG->idisp)*dkx)-2.0)/dx1sq) +
                    ((2.0*cos((j+pG->jdisp)*dky)-2.0)/dx2sq) +
                    ((2.0*cos((k+pG->kdisp)*dkz)-2.0)/dx3sq));
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0] *= pcoeff;
      work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][1] *= pcoeff;
    }
  }}

/* Backward FFT and set potential in real space.  Normalization of Phi is over
 * total number of cells in Domain */

  ath_3d_fft(bplan3d, work);

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] = work[F3DI(i-is,j-js,k-ks,pG->Nx1,pG->Nx2,pG->Nx3)][0]/
        bplan3d->gcnt;
    }
  }}

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

#endif /* SELF_GRAVITY_USING_FFT */
#endif /* FFT_ENABLED */
