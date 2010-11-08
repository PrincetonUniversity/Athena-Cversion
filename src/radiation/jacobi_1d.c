#include "../copyright.h"
/*==============================================================================
 * FILE: jacobi_1d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 1D grid using jacobi's method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution_1d.c()
 *   formal_solution_1d_destruct()
 *   formal_solution_1d_init()
 *   jacobi_pass_pointers_to_mg_1d()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION
#ifdef JACOBI

#ifdef RAD_MULTIG
static Real ******psi = NULL;
#else
static Real *****psi = NULL;
#endif
static Real *mu1 = NULL;
static Real ***imuo = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_1d()     - computes a single sweep in one direction (up or down)
 *   update_sfunc() - updates source function after compute of mean intensity
 *============================================================================*/
static void sweep_1d(RadGridS *pRG, int sx);
static void update_sfunc(RadS *R);

void formal_solution_1d(RadGridS *pRG)
{
  int i, l;
  int nmu = pRG->nmu, nmu2 = pRG->nmu / 2;
  int ng = pRG->ng, ng2 = pRG->ng / 2;
  int ks = pRG->ks, js = pRG->js;
  int is = pRG->is, ie = pRG->ie; 
  int ifr, nf = pRG->nf;

/* initialize mean intensities at all depths to zero */
  for(i=is-1; i<=ie+1; i++)
    for(ifr=0; ifr<nf; ifr++)
      pRG->R[ks][js][i][ifr].J = 0.0;

/* Compute formal solution for all downward going rays in 
 * each vertical gridzone */

/* Account for ix1 boundary condition */
  for(ifr=0; ifr<nf; ifr++)
    for(l=nmu2; l<nmu; l++) 
      imuo[ifr][l][0] = pRG->l1imu[ks][js][ifr][l][0];

/* sweep forward in x1 */
  sweep_1d(pRG, 1);

/* update outward emission at ox1 boundary */
  for(ifr=0; ifr<nf; ifr++)
    for(l=nmu2; l<nmu; l++) 
      pRG->r1imu[ks][js][ifr][l][0] = imuo[ifr][l][0];

/* Compute formal solution for all upward going rays in 
 * each vertical gridzone */

/* Account for ox1 boundary condition */
  for(ifr=0; ifr<nf; ifr++)
    for(l=0; l<nmu2; l++) 
      imuo[ifr][l][0] = pRG->r1imu[ks][js][ifr][l][0];

/* sweep upward */
  sweep_1d(pRG, -1);
  
/* update outward emission at ix1 boundary */
  for(ifr=0; ifr<nf; ifr++)
    for(l=0; l<nmu2; l++) 
      pRG->l1imu[ks][js][ifr][l][0] = imuo[ifr][l][0];

/* update source function */
  for(i=is; i<=ie; i++) 
    for(ifr=0; ifr<nf; ifr++)
      update_sfunc(&(pRG->R[ks][js][i][ifr]));

  return;
}

static void sweep_1d(RadGridS *pRG, int sx)
{
  int it0, i, l;
  int nmu = pRG->nmu, nmu2 = pRG->nmu/2;
  int ng = pRG->ng, ng2 = pRG->ng/2;
  int js = pRG->js, ks = pRG->ks;
  int is = pRG->is, ie = pRG->ie;
  int ifr, nf = pRG->nf;
  Real imu, S0, S2;
  Real am, am1, bm, bm1;
  Real mus, mue, jacobi, ge;
  Real deltas;

  if (sx == 1) {
    mus = nmu2; mue = nmu-1;
  } else {
    mus = 0; mue = nmu2-1;
  }

  for(it0=is; it0<=ie; it0++) {
     if (sx == 1) 
       i = it0;
     else
       i = ie + is - it0;

     for(ifr=0; ifr<nf; ifr++) {
       for(l=mus; l<=mue; l++) {

	 S0 = pRG->R[ks][js][i-sx][ifr].S;
	 S2 = pRG->R[ks][js][i+sx][ifr].S;

#ifdef RAD_MULTIG
	 imu = psi[img][i][ifr][l][0][0] * imuo[ifr][l][0] + 
	       psi[img][i][ifr][l][0][1] * S0 +
	       psi[img][i][ifr][l][0][2] * pRG->R[ks][js][i][ifr].S +
	       psi[img][i][ifr][l][0][3] * S2;
#else
	 imu = psi[i][ifr][l][0][0] * imuo[ifr][l][0] + 
	       psi[i][ifr][l][0][1] * S0 +
	       psi[i][ifr][l][0][2] * pRG->R[ks][js][i][ifr].S +
	       psi[i][ifr][l][0][3] * S2;
#endif 
/* Add to mean intensity and save for next iteration */
	 pRG->R[ks][js][i][ifr].J += pRG->w[l][0] * imu;
	 imuo[ifr][l][0] = imu;
       }
     }
  }
  return;
}

static void update_sfunc(RadS *R)
{
  Real snew, deltas;
  
  snew = (1.0 - R->eps) * R->J + R->eps * R->B;
  deltas = (snew - R->S) / (1.0 - (1.0 - R->eps) * R->lamstr);
  R->S += deltas;
  return;
}

void formal_solution_1d_destruct(void)
{

  int i;

#ifdef RAD_MULTIG
  if (psi != NULL) {
    for(i=0; i<nmgrid; i++)
      if (psi[i] != NULL) free_5d_array(psi[i]);
    free(psi);
  }
#else
  if (psi != NULL) free_5d_array(psi);
#endif
  if (imuo != NULL) free_3d_array(imuo);
  if (mu1 != NULL) free(mu1);

  return;
}

void formal_solution_1d_init(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0];
  int nmu = pRG->nmu, nmu2 = pRG->nmu/2, ng = pRG->ng;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, ks = pRG->ks;
  int ifr, nf = pRG->nf;
  int i, l;
  int sx;
  Real dx = pRG->dx1;
  Real edtau, a0, a1, a2;
  Real chim, chio, chip, dtaum, dtaup;

  if ((imuo = (Real ***)calloc_3d_array(nf,nmu,ng,sizeof(Real))) == NULL) 
    goto on_error;

  if ((mu1 = (Real *)calloc(nmu,sizeof(Real))) == NULL)
    goto on_error;

  for(l=0; l<nmu; l++) 
    mu1[l] = fabs(1.0 / pRG->mu[l]);

#ifdef RAD_MULTIG

  if ((psi = (Real ******)calloc(nmgrid,sizeof(Real*****))) == NULL) 
    goto on_error;
  for(i=0; i<nmgrid; i++)
    psi[i] = NULL;

  if ((psi[0] = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  for(i=is; i<=ie; i++) 
    for(ifr=0; ifr<nf; ifr++) {
      chio = pRG->R[ks][js][i][ifr].chi;
      for(l=0; l<nmu; l++) {
	if(l < nmu2) sx = -1; else sx = 1;

	chim = pRG->R[ks][js][i-sx][ifr].chi;
	chip = pRG->R[ks][js][i+sx][ifr].chi;
      
	dtaum = 0.5 * (chim + chio) * dx * mu1[l]; 
	dtaup = 0.5 * (chip + chio) * dx * mu1[l]; 
	get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);	

	psi[0][i][ifr][l][0][0] = edtau;
	psi[0][i][ifr][l][0][1] = a0;
	psi[0][i][ifr][l][0][2] = a1;
	psi[0][i][ifr][l][0][3] = a2;
	pRG->R[ks][js][i][ifr].lamstr += pRG->w[l][0] * a1;

      }
    }
#else
  if ((psi = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  for(i=is; i<=ie; i++) 
    for(ifr=0; ifr<nf; ifr++) {
      chio = pRG->R[ks][js][i][ifr].chi;
      for(l=0; l<nmu; l++) {
	if(l < nmu2) sx = -1; else sx = 1;

	chim = pRG->R[ks][js][i-sx][ifr].chi;
	chip = pRG->R[ks][js][i+sx][ifr].chi;
      
	dtaum = 0.5 * (chim + chio) * dx * mu1[l]; 
	dtaup = 0.5 * (chip + chio) * dx * mu1[l]; 
	get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);	

	psi[i][ifr][l][0][0] = edtau;
	psi[i][ifr][l][0][1] = a0;
	psi[i][ifr][l][0][2] = a1;
	psi[i][ifr][l][0][3] = a2;
	pRG->R[ks][js][i][ifr].lamstr += pRG->w[l][0] * a1;

      }
    }
#endif /* RAD_MULTIG */
  
  return;

  on_error:
  formal_solution_1d_destruct();
  ath_error("[formal_solution__1d_init]: Error allocating memory\n");
  return;

}

#ifdef RAD_MULTIG
void jacobi_pass_pointers_to_mg_1d(Real *******psi0, Real **mu10)
{
  *psi0 = psi;
  *mu10 = mu1;

  return;
}
#endif /* RAD_MULTIG */

#endif /* JACOBI */
#endif /* RADIATION */
