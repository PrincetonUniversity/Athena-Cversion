#include "../copyright.h"
/*==============================================================================
 * FILE: multig_1d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 1D grid using multigrid method.  Works with Jacobi
 *          or Gauss-Seidel formal solution algorithms.  Calls private
 *          multigrid_1d() function which handles bulk of the computation.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   fomal_solution_mg_1d.c()
 *============================================================================*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER
#ifdef RAD_MULTIG

static int it0 = 1;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   multigrid_1d()   - recursive function which contains multigrid algorithm 
 *   restrict_1d()    - 1d implementation of restriction (fine-to-coarse)
 *   prolongate_1d()  - 1d implementation of prolongation (coarse-to-fine)
 *   formal_solution_mg_1d_psi() - computes psi matrix for each grid in multigrid
 *============================================================================*/
void multigrid_1d(RadGridS *mG, Real *dSrmax);
static void restrict_1d(RadGridS *fmG, RadGridS *cmG);
static void prolongate_1d(RadGridS *cmG, RadGridS *fmG);
static void formal_solution_mg_1d_psi(RadGridS *pRG);

void formal_solution_mg_1d(RadGridS *pRG, Real *dSrmax)
{
  //RadGridS *pRG=(pD->RadGrid);
  int i;
  
  //for(i=0; i < 0; i++) 
  //  formal_solution_1d(pRG);

  multigrid_1d(pRG,dSrmax);

  for(i=0; i < 40; i++) {
    formal_solution_1d(pRG,dSrmax);
    bvals_rad(pRG);
    //output_mean_intensity_1d(pRG,it0);
    it0++;
  }

  return;
}

void multigrid_1d(RadGridS *mG, Real *dSrmax)
{
  RadGridS cmG;
  int i, j, k;
  int is, ie, js, ks;
  int nfsf = 8, nfs = 1;

  if (mG->Nx[0] <= 16) {
    for(i=0; i < nfsf; i++) {
      formal_solution_1d(mG,dSrmax);
      bvals_rad(mG);
      //output_mean_intensity_1d(mG,it0);
      it0++;
    }
  } else {
    
    for(i=0; i < nfs; i++) {
      formal_solution_1d(mG,dSrmax);
      bvals_rad(mG);
      //output_mean_intensity_1d(mG,it0);
      it0++;
    }

    cmG.Nx[2] = 1;
    cmG.Nx[1] = 1;
    cmG.Nx[0] = mG->Nx[0]/2;
    cmG.dx1 = 2.0 * mG->dx1; 
    cmG.nf = mG->nf;
    cmG.nmu = mG->nmu;
    cmG.ng = 1;
    cmG.ks = 0;
    cmG.js = 0;
    cmG.is = 1;
    cmG.ie = cmG.Nx[0];
    cmG.mu = NULL;
    cmG.gamma = NULL;  
    cmG.r3imu = NULL;
    cmG.l3imu = NULL;
    cmG.r2imu = NULL;
    cmG.l2imu = NULL;
    cmG.ix1_RBCFun = mG->ix1_RBCFun;
    cmG.ox1_RBCFun = mG->ox1_RBCFun;
    cmG.rx1_id = mG->rx1_id;
    cmG.lx1_id = mG->lx1_id;

    js = cmG.js; ks = cmG.ks;

    if ((cmG.R = (RadS ****)calloc_4d_array(1,1,cmG.Nx[0]+2,cmG.nf,sizeof(RadS))) == NULL) 
      ath_error("[multigrid_1d]: Error allocating memory for multigrid radiation.\n");
    
    if ((cmG.w = (Real **)calloc_2d_array(cmG.nmu,1,sizeof(Real))) == NULL) 
      ath_error("Error allocating memory for w.\n");
    
    if ((cmG.r1imu = (Real *****)calloc_5d_array(1,1,cmG.nf,cmG.nmu,1,sizeof(Real))) == NULL) 
      ath_error("[multigrid_1d]: Error allocating memory for r3imu.\n");
    
    if ((cmG.l1imu = (Real *****)calloc_5d_array(1,1,cmG.nf,cmG.nmu,1,sizeof(Real))) == NULL)
      ath_error("[multigrid_1d]: Error allocating memory for l3imu.\n");

    for(i=0; i<cmG.nmu; i++) {
      cmG.w[i][0] = mG->w[i][0];
      cmG.r1imu[ks][js][0][i][0] = mG->r1imu[ks][js][0][i][0];
      cmG.l1imu[ks][js][0][i][0] = mG->l1imu[ks][js][0][i][0];
    }

    restrict_1d(mG, &cmG);    
    formal_solution_mg_1d_psi(&cmG);
    multigrid_1d(&cmG);
    prolongate_1d(&cmG, mG);
  
    for(i=0; i < nfs; i++) {
      formal_solution_1d(mG,dSrmax);
      bvals_rad(mG);
      //      output_mean_intensity_1d(mG,it0);
      it0++;
    }
    radgrid_destruct(&cmG);
  }
  
  return;
}

static void restrict_1d(RadGridS *fmG, RadGridS *cmG) 
{
  int i;
  int ifr = 0;
  int isc = cmG->is, iec = cmG->ie;
  int isf = fmG->is, ief = fmG->ie;
  int js = cmG->js, ks = cmG->ks;

  // uses boundary values in restriction
  for (i=cmG->is; i<=cmG->ie; i++) {
    cmG->R[ks][js][ifr][i].J = 0.0;
    cmG->R[ks][js][i][ifr].B = 
      0.25 * (fmG->R[ks][js][2*i  ][ifr].B   + fmG->R[ks][js][2*i-2][ifr].B) +
      0.5  *  fmG->R[ks][js][2*i-1][ifr].B;
    cmG->R[ks][js][i][ifr].S = 
      0.25 * (fmG->R[ks][js][2*i  ][ifr].S   + fmG->R[ks][js][2*i-2][ifr].S) +
      0.5  *  fmG->R[ks][js][2*i-1][ifr].S;    
    cmG->R[ks][js][i][ifr].eps = 
      0.25 * (fmG->R[ks][js][2*i  ][ifr].eps + fmG->R[ks][js][2*i-2][ifr].eps) +
      0.5  *  fmG->R[ks][js][2*i-1][ifr].eps;
    cmG->R[ks][js][ifr][i].chi = 
      0.25 * (fmG->R[ks][js][2*i  ][ifr].chi + fmG->R[ks][js][2*i-2][ifr].chi) +
      0.5  *  fmG->R[ks][js][2*i-1][ifr].chi;
  }
  cmG->R[ks][js][isc-1][ifr].B   = fmG->R[ks][js][isf-1][ifr].B;
  cmG->R[ks][js][isc-1][ifr].S   = fmG->R[ks][js][isf-1][ifr].S;
  cmG->R[ks][js][isc-1][ifr].eps = fmG->R[ks][js][isf-1][ifr].eps;
  cmG->R[ks][js][isc-1][ifr].chi = fmG->R[ks][js][isf-1][ifr].chi;
  cmG->R[ks][js][iec+1][ifr].B   = fmG->R[ks][js][ief+1][ifr].B;
  cmG->R[ks][js][iec+1][ifr].S   = fmG->R[ks][js][ief+1][ifr].S;
  cmG->R[ks][js][iec+1][ifr].eps = fmG->R[ks][js][ief+1][ifr].eps;
  cmG->R[ks][js][iec+1][ifr].chi = fmG->R[ks][js][ief+1][ifr].chi;

  img++;
  return;
}

static void prolongate_1d(RadGridS *cmG, RadGridS *fmG) 
{
  int i;
  int ifr = 0;
  int ks = fmG->ks, js = fmG->js;

  for (i=fmG->is; i<=fmG->ie; i+=2) 
    fmG->R[ks][js][i][ifr].S = cmG->R[ks][js][i/2+1][ifr].S;

  for (i=fmG->is+1; i<=fmG->ie-1; i+=2) 
    fmG->R[ks][js][i    ][ifr].S = 0.5 * (cmG->R[ks][js][i/2][ifr].S + 
    cmG->R[ks][js][i/2+1][ifr].S);
  
  img--;
  return;
}

void formal_solution_mg_1d_psi(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0];
  int nmu = pRG->nmu, nmu2 = pRG->nmu/2, ng = pRG->ng;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, ks = pRG->ks;
  int ifr, nf = pRG->nf;
  int i, it4;
  int sx;
  Real dx = pRG->dx1;
  Real edtau, a0, a1, a2;
  Real chim, chio, chip, dtaum, dtaup;
  Real ******psi = NULL, *mu1 = NULL;
  Real ***psiint = NULL;

#if defined(JACOBI)
  jacobi_pass_pointers_to_mg_1d(&psi, &mu1);
#elif defined(GAUSSEID)
  gs_pass_pointers_to_mg_1d(&psi, &mu1, &psiint);
  if ((psiint[img] = (Real **)calloc_2d_array(nx1+2,nf,sizeof(Real))) == NULL) 
    goto on_error;
#endif

  if ((psi[img] = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  for(i=is; i<=ie; i++) 
    for(ifr=0; ifr<nf; ifr++) {
      chio = pRG->R[ks][js][i][ifr].chi;
      for(it4=0; it4<nmu; it4++) {
	if(it4 < nmu2) sx = -1; else sx = 1;

	chim = pRG->R[ks][js][i-sx][ifr].chi;
	chip = pRG->R[ks][js][i+sx][ifr].chi;
      
	dtaum = 0.5 * (chim + chio) * dx * mu1[it4]; 
	dtaup = 0.5 * (chip + chio) * dx * mu1[it4]; 
	get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);	
    
	psi[img][i][ifr][it4][0][0] = edtau;
	psi[img][i][ifr][it4][0][1] = a0;
	psi[img][i][ifr][it4][0][2] = a1;
	psi[img][i][ifr][it4][0][3] = a2;
	pRG->R[ks][js][i][ifr].lamstr += pRG->w[it4][0] * a1;
#ifdef GAUSSEID
	if (sx == 1) 
	  psiint[img][i][ifr] += pRG->w[it4][0] * a2;
#endif
      }
    }

  return;

  on_error:
  formal_solution_1d_destruct();
  ath_error("[formal_solution_mg_1d_psi]: Error allocating memory\n");
  return;

}

#endif /* RAD_MULTIG */
#endif /* RADIATION_TRANSFER */
