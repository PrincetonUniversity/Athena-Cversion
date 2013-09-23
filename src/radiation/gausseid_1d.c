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

/* Working arrays used in formal solution */
static Real **psiint = NULL;
static Real **lamstr = NULL;
static Real ***imuo = NULL;

/* Variables used for succesive over-relaxation */
static Real sorw, dsrold, dsrmx;
static int isor;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_1d()     - computes a single sweep in one direction (up or down)
 *   update_sfunc() - updates source function after compute of mean intensity
 *============================================================================*/
static void sweep_1d(RadGridS *pRG, int sx, int ifr);
static void update_sfunc(RadS *R, Real *deltas, Real lamstr);

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_1d(RadGridS *pRG, Real *dSrmax, int ifr)
 *  \brief formal solution for single freq. in 1D using Gauss-Seidel method */
void formal_solution_1d(RadGridS *pRG, Real *dSrmax, int ifr)
{
  int i, m;
  int ks = pRG->ks, js = pRG->js;
  int is = pRG->is, ie = pRG->ie; 
  int nf = pRG->nf, nang = pRG->nang;

/* Initialize dsrmx */
  dsrmx = 0.0;

/* initialize mean intensities at all depths to zero */
  for(i=is-1; i<=ie+1; i++) {
    pRG->R[ifr][ks][js][i].J = 0.0;
    pRG->R[ifr][ks][js][i].H[0] = 0.0;
    pRG->R[ifr][ks][js][i].K[0] = 0.0;
    lamstr[ifr][i] = 0.0;
    psiint[ifr][i] = 0.0;
  }

/* Account for ix1 boundary condition */
  for(m=0; m<nang; m++) 
    imuo[ifr][0][m] = pRG->Ghstl1i[ifr][ks][js][0][m];

/* sweep forward in x1 */
  sweep_1d(pRG, 1, ifr);

/* update outward emission at ox1 boundary */
  for(m=0; m<nang; m++) 
    pRG->r1imu[ifr][ks][js][0][m] = imuo[ifr][0][m];

/* Compute formal solution for all upward going rays in 
 * each vertical gridzone */

/* Account for ox1 boundary condition */
  for(m=0; m<nang; m++) 
    imuo[ifr][1][m] = pRG->Ghstr1i[ifr][ks][js][1][m];

/* sweep upward */
  sweep_1d(pRG, -1, ifr);
  
/* update outward emission at ix1 boundary */
  for(m=0; m<nang; m++) 
    pRG->l1imu[ifr][ks][js][1][m] = imuo[ifr][1][m];

/* Evaluate relative change if SOR */
  if ((dsrmx < 0.03) && (isor == 1) ) {
    sorw = 2.0 / (1.0 + sqrt(1.0 - dsrmx/dsrold));
    isor = 0;
  }
  dsrold = dsrmx;
  (*dSrmax) = dsrmx;

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_1d_destruct(void)
 *  \brief free temporary working arrays */
void formal_solution_1d_destruct(void)
{

  if (psiint != NULL) free_2d_array(psiint);
  if (lamstr != NULL) free_2d_array(lamstr);
  if (imuo != NULL)   free_3d_array(imuo);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_1d_init(DomainS *pD)
 *  \brief Allocate memory for working arrays */
void formal_solution_1d_init(DomainS *pD)
{
  RadGridS *pRG, *pROG;
  int nx1;
  int nf, nang;

  if (radt_mode == 0) { /* integration only */
    pRG = pD->RadGrid;
    nx1 = pRG->Nx[0];
    nf = pRG->nf; 
    nang = pRG->nang;
  } else if (radt_mode == 1) { /* output only */
    pRG = pD->RadOutGrid;
    nx1 = pRG->Nx[0];
    nf = pRG->nf;
    nang = pRG->nang;
  } else if (radt_mode == 2) { /* integration and output */
    pRG = pD->RadGrid;
    pROG = pD->RadOutGrid;
    nx1 = pRG->Nx[0];
    nf = MAX(pRG->nf,pROG->nf);
    nang = MAX(pRG->nang,pROG->nang);
  }

/* Initialize variables for sor*/
  isor = par_geti_def("radiation","isor",0);
  sorw = 1.0;
  dsrold = 0.0;

  if ((imuo = (Real ***)calloc_3d_array(nf,2,nang,sizeof(Real))) == NULL) 
    goto on_error;

  if ((lamstr = (Real **)calloc_2d_array(nf,nx1+2,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psiint = (Real **)calloc_2d_array(nf,nx1+2,sizeof(Real))) == NULL) 
    goto on_error;
  
  return;

  on_error:
  formal_solution_1d_destruct();
  ath_error("[formal_solution__1d_init]: Error allocating memory\n");
  return;

}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_1d(RadGridS *pRG, int sx, int ifr)
 *  \brief Perform Gauss-Seidel sweep in one direction */
static void sweep_1d(RadGridS *pRG, int sx, int ifr)
{
  int it0, i, l, m;
  int js = pRG->js, ks = pRG->ks;
  int is = pRG->is, ie = pRG->ie;
  int nf = pRG->nf, nang = pRG->nang;
  Real imu, wimu, S0, S2;
  Real deltas;
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
     chio = pRG->R[ifr][ks][js][i].chi;
     for(m=0; m<nang; m++) {
       chim = pRG->R[ifr][ks][js][i-sx].chi;
       chip = pRG->R[ifr][ks][js][i+sx].chi;
       interp_quad_chi(chim,chio,chip,&dtaum);
       interp_quad_chi(chip,chio,chim,&dtaup);
       dtaum *= dx / pRG->mu[0][m][0]; 
       dtaup *= dx / pRG->mu[0][m][0];
/* ---------  compute intensity at grid center and add to mean intensity ------- */
       interp_quad_source_slope_lim(dtaum, dtaup, &edtau, &a0, &a1, &a2,
				    S0, pRG->R[ifr][ks][js][i].S, S2);
       imu = a0 * S0 + a1 * pRG->R[ifr][ks][js][i].S + a2 * S2;
       imu += edtau * imuo[ifr][l][m];
       lamstr[ifr][i] += pRG->wmu[m] * a1;
       if (sx == 1) psiint[ifr][i] += pRG->wmu[m] * a2; 
     
/* Add to mean intensity and save for next iteration */
       wimu = pRG->wmu[m] * imu;
       pRG->R[ifr][ks][js][i].J += wimu;
       pRG->R[ifr][ks][js][i].H[0] += pRG->mu[l][m][0] * wimu;
       pRG->R[ifr][ks][js][i].K[0] += pRG->mu[l][m][0] * pRG->mu[l][m][0] * wimu;
       imuo[ifr][l][m] = imu;
     }
     if (sx == -1) {
       update_sfunc(&(pRG->R[ifr][ks][js][i]), &deltas, lamstr[ifr][i]);
       for(m=0; m<nang; m++) 
	 imuo[ifr][l][m] += deltas * pRG->wmu[m] * a1;	          
       pRG->R[ifr][ks][js][i-1].J += deltas * psiint[ifr][i-1];
     }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_sfunc(RadS *R, Real *dSr, Real lamstr)
 *  \brief Gauss-Seidel update of source function with new mean intensity */
static void update_sfunc(RadS *R, Real *deltas, Real lamstr)
{
  Real snew, deltasr;
  
  snew = (1.0 - R->eps) * R->J + R->eps * R->B + R->Snt;

  (*deltas) = (snew - R->S) / (1.0 - (1.0 - R->eps) * lamstr);
  deltasr = fabs(*deltas) / R->S;
  (*deltas) *= sorw;
  R->S += (*deltas);

  if (deltasr > dsrmx) dsrmx = deltasr; 
  return;
}


#endif /* GAUSSEID */
#endif /* RADIATION_TRANSFER */
