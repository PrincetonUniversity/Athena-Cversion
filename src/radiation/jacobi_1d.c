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
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER
#if defined(JACOBI)

static Real *****psi = NULL;
static Real **lamstr = NULL;
static Real *muinv = NULL, *mu2 = NULL;
static Real ***imuo = NULL;
static Real **Jold = NULL;
static int svwght;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_1d()     - computes a single sweep in one direction (up or down)
 *   update_sfunc() - updates source function after compute of mean intensity
 *============================================================================*/
static void sweep_1d(RadGridS *pRG, int sx, int ifr);
static void update_sfunc(RadS *R, Real *dSr, Real lamstr);

void formal_solution_1d(RadGridS *pRG, Real *dSrmax, int ifr)
{
  int i, m;
  int ks = pRG->ks, js = pRG->js;
  int is = pRG->is, ie = pRG->ie; 
  int nang = pRG->nang;
  Real dSr, dJ, dJmax;

/* if LTE then store J values from previous iteration */
  if(lte != 0)
    for(i=is; i<=ie; i++) 
      Jold[ifr][i] = pRG->R[ifr][ks][js][i].J;

/* initialize mean intensities at all depths to zero */
  for(i=is; i<=ie; i++) {
    pRG->R[ifr][ks][js][i].J = 0.0;
    pRG->R[ifr][ks][js][i].H[0] = 0.0;
    pRG->R[ifr][ks][js][i].K[0] = 0.0;
    if (svwght == 0) lamstr[ifr][i] = 0.0;
  }

/* Account for ix1 boundary condition */
  for(m=0; m<nang; m++) 
    imuo[ifr][0][m] = pRG->Ghstl1i[ifr][ks][js][0][m];

/* sweep forward in x1 */
  sweep_1d(pRG, 1, ifr);

/* update outward emission at ox1 boundary */
  for(m=0; m<nang; m++) 
    pRG->r1imu[ifr][ks][js][0][m] = imuo[ifr][0][m];

/* Account for ox1 boundary condition */
  for(m=0; m<nang; m++) 
    imuo[ifr][1][m] = pRG->Ghstr1i[ifr][ks][js][1][m];

/* sweep upward */
  sweep_1d(pRG, -1, ifr);
  
/* update outward emission at ix1 boundary */
  for(m=0; m<nang; m++) 
    pRG->l1imu[ifr][ks][js][1][m] = imuo[ifr][1][m];

  if(lte == 0) {
/* update source function */
    (*dSrmax) = 0.0;
    for(i=is; i<=ie; i++) {
      update_sfunc(&(pRG->R[ifr][ks][js][i]),&dSr,lamstr[ifr][i]);
      if( dSr > (*dSrmax)) (*dSrmax) = dSr;
    }
  } else {
    (*dSrmax) = 0.0;
    dJmax = 0.0;
    for(i=is; i<=ie; i++) {
      dJ = fabs(pRG->R[ifr][ks][js][i].J - Jold[ifr][i]);
      if(dJ > dJmax) dJmax = dJ;
      if (Jold[ifr][i] > 0.0)
	dSr = dJ / Jold[ifr][i];
      else
	dSr = 0;
      if( dSr > (*dSrmax)) (*dSrmax) = dSr;	 
    }
    if(((*dSrmax) == 0.0) && (dJmax > 0.0)) (*dSrmax) = 1.0;
  }

  return;
}

static void sweep_1d(RadGridS *pRG, int sx, int ifr)
{
  int it0, i, l, m;
  int js = pRG->js, ks = pRG->ks;
  int is = pRG->is, ie = pRG->ie;
  int nang = pRG->nang;
  Real imu, wimu, S0, S2;
  Real chio, chim, chip;
  Real dtaum, dtaup, dx = pRG->dx1;
  Real edtau, a0, a1, a2, a3;

  if (sx == 1) l = 0; else l = 1;

  for(it0=is; it0<=ie; it0++) {
     if (sx == 1) 
       i = it0;
     else
       i = ie + is - it0;

     S0 = pRG->R[ifr][ks][js][i-sx].S;
     S2 = pRG->R[ifr][ks][js][i+sx].S;
     if(svwght == 0) {
       chio = pRG->R[ifr][ks][js][i].chi;
       chim = pRG->R[ifr][ks][js][i-sx].chi;
       chip = pRG->R[ifr][ks][js][i+sx].chi;
	 /*dtaum = 0.5 * (chim + chio);
	   dtaup = 0.5 * (chip + chio);*/
       interp_quad_chi(chim,chio,chip,&dtaum);
       interp_quad_chi(chip,chio,chim,&dtaup);
       dtaum *= dx; 
       dtaup *= dx;
     }
     for(m=0; m<nang; m++) {
       if(svwght == 1) {
	 imu = psi[ifr][i][l][m][1] * S0 +
	   psi[ifr][i][l][m][2] * pRG->R[ifr][ks][js][i].S +
	   psi[ifr][i][l][m][3] * S2;
	 if (imu < 0.0) imu=0.0;
	 imu += psi[ifr][i][l][m][0] * imuo[ifr][l][m];	
       } else {
	 /*interp_quad_source(dtaum*muinv[m],dtaup*muinv[m], &edtau, &a0, &a1, &a2,
	   S0, pRG->R[ifr][ks][js][i].S, S2);*/
	 interp_quad_source_slope_lim(dtaum*muinv[m],dtaup*muinv[m], &edtau, &a0, &a1, &a2,
				      S0, pRG->R[ifr][ks][js][i].S, S2);
	 imu = a0 * S0 + a1 * pRG->R[ifr][ks][js][i].S + a2 * S2;
	 imu += edtau * imuo[ifr][l][m];
	 lamstr[ifr][i] += pRG->wmu[m] * a1;
       } 
/* Add to mean intensity and save for next iteration */
       wimu = pRG->wmu[m] * imu;
       pRG->R[ifr][ks][js][i].J += wimu;
       pRG->R[ifr][ks][js][i].H[0] += pRG->mu[l][m][0] * wimu;
       pRG->R[ifr][ks][js][i].K[0] += mu2[m] * wimu;
       imuo[ifr][l][m] = imu;
     }
  }
  return;
}

static void update_sfunc(RadS *R, Real *dSr, Real lamstr)
{
  Real Snew, deltaS;
  
  Snew = (1.0 - R->eps) * R->J + R->eps * R->B;
  deltaS = (Snew - R->S) / (1.0 - (1.0 - R->eps) * lamstr);
  if (R->S > 0.0) (*dSr) = fabs(deltaS/R->S);
  R->S += deltaS;
  return;
}

void formal_solution_1d_destruct(void)
{

  int i;

  if (psi != NULL) free_5d_array(psi);
  if (lamstr != NULL) free_2d_array(lamstr);
  if (imuo != NULL) free_3d_array(imuo);
  if (muinv != NULL) free(muinv);
  if (mu2 != NULL) free(mu2);
  if (Jold != NULL) free_2d_array(Jold);

  return;
}

void formal_solution_1d_init(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0];
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;
  int ifr, i, l, m;
  int sx;
  Real dx = pRG->dx1;
  Real edtau, a0, a1, a2;
  Real chim, chio, chip, dtaum, dtaup;
  Real S0, S2;

  svwght = par_geti("radiation","svwght");

  if ((imuo = (Real ***)calloc_3d_array(nf,2,nang,sizeof(Real))) == NULL) 
    goto on_error;

  if ((muinv = (Real *)calloc(nang,sizeof(Real))) == NULL)
    goto on_error;

  if ((mu2 = (Real *)calloc(nang,sizeof(Real))) == NULL)
    goto on_error;

  for(m=0; m<nang; m++) {
    muinv[m] = fabs(1.0 / pRG->mu[0][m][0]);
    mu2[m] = pRG->mu[0][m][0] * pRG->mu[0][m][0];
  }

  if ((lamstr = (Real **)calloc_2d_array(nf,nx1+2,sizeof(Real))) == NULL) 
    goto on_error;

  if (lte != 0) 
    if ((Jold = (Real **)calloc_2d_array(nf,nx1+2,sizeof(Real))) == NULL)
      goto on_error;

  if(svwght == 1) {
/* compute weights once and save for next iteration */
    if ((psi = (Real *****)calloc_5d_array(nf,nx1+2,2,nang,4,sizeof(Real))) == NULL) 
      goto on_error;

    for(i=is; i<=ie; i++) 
      for(ifr=0; ifr<nf; ifr++) {
	lamstr[ifr][i] = 0.0;
	chio = pRG->R[ifr][ks][js][i].chi;
	for(l=0; l<2; l++) {
	  if(l == 0) sx = 1; else sx = -1;
	  S0 = pRG->R[ifr][ks][js][i-sx].S;
	  S2 = pRG->R[ifr][ks][js][i+sx].S;
	  for(m=0; m<nang; m++) { 
	    chim = pRG->R[ifr][ks][js][i-sx].chi;
	    chip = pRG->R[ifr][ks][js][i+sx].chi;
	    
	    /*dtaum = 0.5 * (chim + chio) * dx * muinv[m]; 
	      dtaup = 0.5 * (chip + chio) * dx * muinv[m];*/ 
	    interp_quad_chi(chim,chio,chip,&dtaum);
	    interp_quad_chi(chip,chio,chim,&dtaup);
	    dtaum *= dx * muinv[m]; 
	    dtaup *= dx * muinv[m];
	    interp_quad_source(dtaum, dtaup, &edtau, &a0, &a1, &a2,
	    		       S0, pRG->R[ifr][ks][js][i].S, S2);

	    psi[ifr][i][l][m][0] = edtau;
	    psi[ifr][i][l][m][1] = a0;
	    psi[ifr][i][l][m][2] = a1;
	    psi[ifr][i][l][m][3] = a2;
	    lamstr[ifr][i] += pRG->wmu[m] * a1;
	  }
	}
      }
  }
  
  return;

  on_error:
  formal_solution_1d_destruct();
  ath_error("[formal_solution__1d_init]: Error allocating memory\n");
  return;

}

#endif /* JACOBI */
#endif /* RADIATION_TRANSFER */
