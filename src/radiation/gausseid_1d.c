#include "../copyright.h"
/*==============================================================================
 * FILE: gausseid_1d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a grid using Gauss-Seidel method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *          SOR method is also implemented via this function
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution_1d.c()
 *   formal_solution_1d_destruct()
 *   formal_solution_1d_init()
 *   gausseid_pass_pointers_to_mg_1d()
 *============================================================================*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER
#ifdef GAUSSEID

#ifdef RAD_MULTIG
static Real ******psi = NULL;
static Real ***psiint = NULL;
#else
static Real *****psi = NULL;
static Real **psiint = NULL;
#endif
static Real *mu1 = NULL, *mu2 = NULL;
static Real ***imuo = NULL;
static Real sorw, dsrold, dsrmx;
static int isor;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_1d()     - computes a single sweep in one direction (up or down)
 *   update_sfunc() - updates source function after compute of mean intensity
 *============================================================================*/
static void sweep_1d(RadGridS *pRG, int sz);
static void update_sfunc(RadS *R, Real *deltas);

void formal_solution_1d(RadGridS *pRG, Real *dSrmax)
{
  int i, l;
  int nmu = pRG->nmu, nmu2 = pRG->nmu / 2;
  int ng = pRG->ng, ng2 = pRG->ng / 2;
  int ks = pRG->ks, js = pRG->js;
  int is = pRG->is, ie = pRG->ie; 
  int ifr, nf = pRG->nf;

/* Initialize dsrmx */
  dsrmx = 0.0;

/* initialize mean intensities at all depths to zero */
  for(i=is-1; i<ie+1; i++)
    for(ifr=0; ifr<nf; ifr++) {
      pRG->R[ks][js][i][ifr].J = 0.0;
      pRG->R[ks][js][i][ifr].K[0] = 0.0;
    }
/* Compute formal solution for all downward going rays in 
 *  each vertical gridzone */

/* Account for ix1 boundary condition */
  for(ifr=0; ifr<nf; ifr++)
    for(l=nmu2; l<nmu; l++) 
      imuo[ifr][l][0] = pRG->l1imu[ks][js][ifr][l][0];

/* sweep downward z */
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

/* Evaluate relative change if SOR */
  if ((dsrmx < 0.03) && (isor == 1) ) {
    sorw = 2.0 / (1.0 + sqrt(1.0 - dsrmx/dsrold));
    isor = 0;
  }
  dsrold = dsrmx;
  (*dSrmax) = dsrmx;

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
  Real mus, mue, gs, ge;
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
	 imu = psi[img][i][ifr][l][0][1] * S0 +
	       psi[img][i][ifr][l][0][2] * pRG->R[ks][js][i][ifr].S +
	       psi[img][i][ifr][l][0][3] * S2;
	 if (imu < 0.0) imu=0.0;
	 imu += psi[img][i][ifr][l][0][0] * imuo[ifr][l][0]; 
#else
	 imu = psi[i][ifr][l][0][1] * S0 +
	       psi[i][ifr][l][0][2] * pRG->R[ks][js][i][ifr].S +
	       psi[i][ifr][l][0][3] * S2;
	 if (imu < 0.0) imu=0.0;
	 imu += psi[i][ifr][l][0][0] * imuo[ifr][l][0];	
#endif      
/* Add imu to mean intensity and save for next iteration */
	 pRG->R[ks][js][i][ifr].J += pRG->w[l][0] * imu;
	 pRG->R[ks][js][i][ifr].K[0] += mu2[l] * pRG->w[l][0] * imu;

	 imuo[ifr][l][0] = imu;
       }

       if (sx == -1) {
	 update_sfunc(&(pRG->R[ks][js][i][ifr]), &deltas);
#ifdef RAD_MULTIG
	 for(l=mus; l<=mue; l++)
	   imuo[ifr][l][0] += deltas * pRG->w[l][0] * 
	     psi[img][i][ifr][l][0][2];       
	 pRG->R[ks][js][i-1][ifr].J += deltas * psiint[img][i-1][ifr];
#else
	 for(l=mus; l<=mue; l++)
	   imuo[ifr][l][0] += deltas * pRG->w[l][0] * 
	     psi[i][ifr][l][0][2];       
	 pRG->R[ks][js][i-1][ifr].J += deltas * psiint[i-1][ifr];
#endif
       }
     }
  }
  return;
}

static void update_sfunc(RadS *R, Real *deltas)
{
  Real snew, deltasr;
  
  snew = (1.0 - R->eps) * R->J + R->eps * R->B;

  (*deltas) = (snew - R->S) / (1.0 - (1.0 - R->eps) * R->lamstr);
  deltasr = fabs(*deltas) / R->S;
  (*deltas) *= sorw;
  R->S += (*deltas);

  if (deltasr > dsrmx) dsrmx = deltasr; 
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
  if (psiint != NULL) {
    for(i=0; i<nmgrid; i++)
      if (psiint[i] != NULL) free_2d_array(psiint[i]);
    free(psiint);
  }
#else
  if (psi != NULL) free_5d_array(psi);
  if (psiint != NULL) free_2d_array(psiint);
#endif
  if (imuo != NULL) free_3d_array(imuo);
  if (mu1 != NULL) free(mu1);
  if (mu2 != NULL) free(mu2);

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

/* Initialize variables for sor*/
  isor = par_geti_def("problem","isor",0);
  sorw = 1.0;
  dsrold = 0.0;

  if ((mu1 = (Real *)calloc(nmu,sizeof(Real))) == NULL)
    goto on_error;

  if ((mu2 = (Real *)calloc(nmu,sizeof(Real))) == NULL)
    goto on_error;

  for(l=0; l<nmu; l++) {
    mu1[l] = fabs(1.0 / pRG->mu[l]);
    mu2[l] = pRG->mu[l] * pRG->mu[l];
  }

#ifdef RAD_MULTIG
  if ((psi = (Real ******)calloc(nmgrid,sizeof(Real*****))) == NULL) 
    goto on_error;
  for(i=0; i<nmgrid; i++)
    psi[i] = NULL;

  if ((psiint = (Real ***)calloc(nmgrid,sizeof(Real**))) == NULL) 
    goto on_error;
  for(i=0; i<nmgrid; i++)
    psiint[i] = NULL;

  if ((psi[0] = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psiint[0] = (Real **)calloc_2d_array(nx1+2,nf,sizeof(Real))) == NULL) 
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
	if (sx == 1) 
	  psiint[0][i][ifr] += pRG->w[l][0] * a2;
      }
    }
#else
  if ((psi = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psiint = (Real **)calloc_2d_array(nx1+2,nf,sizeof(Real))) == NULL) 
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
	if (sx == 1) 
	  psiint[i][ifr] += pRG->w[l][0] * a2; 
      }
    }
#endif

  if ((imuo = (Real ***)calloc_3d_array(nf,nmu,ng,sizeof(Real))) == NULL) 
    goto on_error;

  return;

  on_error:
  formal_solution_1d_destruct();
  ath_error("[formal_solution_1d_init]: Error allocating memory\n");
  return;

}

#ifdef RAD_MULTIG
void gs_pass_pointers_to_mg_1d(Real *******psi0, Real **mu10, Real ****psiint0)
{
  *psi0 = psi;
  *mu10 = mu1;
  *psiint0 = psiint;

  return;
}
#endif /* RAD_MULTIG */

#endif /* GAUSSEID */
#endif /* RADIATION_TRANSFER */
