#include "../copyright.h"
/*==============================================================================
 * FILE: multig_2d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 2D grid using multigrid method.  Works with Jacobi
 *          or Gauss-Seidel formal solution algorithms.  Calls private
 *          multigrid_2d() function which handles bulk of the computation.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   fomal_solution_mg_2d.c()
 *============================================================================*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define OUT_DIAG

#ifdef RADIATION
#ifdef RAD_MULTIG

static int it0 = 1;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   multigrid_2d()   - recursive function which contains multigrid algorithm 
 *   restrict_2d()    - 2d implementation of restriction (fine-to-coarse)
 *   prolongate_2d()  - 2d implementation of prolongation (coarse-to-fine)
 *   formal_solution_mg_2d_psi() - computes psi matrix for each grid in multigrid
 *============================================================================*/
static void multigrid_2d(RadGridS *mG);
static void restrict_2d(RadGridS *fmG, RadGridS *cmG);
static void prolongate_2d(RadGridS *cmG, RadGridS *fmG);
static void formal_solution_mg_2d_psi(RadGridS *pRG);

void formal_solution_mg_2d(RadGridS *pRG)
{
  //RadGridS *pRG=(pD->RadGrid);
  int i;

  //for(i=0; i < 0; i++) {
  //  formal_solution_2d(pRG);
  //    bvals_rad(pRG);
  //}
  multigrid_2d(pRG);

  for(i=0; i < 30; i++) {
    formal_solution_2d(pRG);
    bvals_rad(pRG);
#ifdef OUT_DIAG
    output_mean_intensity_2d(pRG,it0);
    //dumpvtk2d(pRG,it0);
    it0++;
#endif
  }
  return;
}

static void multigrid_2d(RadGridS *mG)
{
  RadGridS cmG;
  int i, j;
  int nmu = mG->nmu, ng = mG->ng, nf = mG->nf;
  int nfsf = 10, nfs = 1;

  if (mG->Nx[1] <= 16) {
    for(i=0; i < nfsf; i++) {
      formal_solution_2d(mG);
      bvals_rad(mG);
#ifdef OUT_DIAG
      output_mean_intensity_2d(mG,it0);
      it0++;
#endif
    }
  } else {
    
    for(i=0; i < nfs; i++) {
      formal_solution_2d(mG);
      bvals_rad(mG);
#ifdef OUT_DIAG
      output_mean_intensity_2d(mG,it0);
      it0++;
#endif
    }

    cmG.Nx[2] = 1;
    cmG.Nx[1] = mG->Nx[1]/2;
    cmG.Nx[0] = mG->Nx[0]/2;
    cmG.dx2 = 2.0 * mG->dx2; 
    cmG.dx1 = 2.0 * mG->dx1;
    cmG.nf = nf;
    cmG.nmu = nmu;
    cmG.ng = ng;
    cmG.ks = 0;
    cmG.ke = 0;
    cmG.js = 1;
    cmG.je = cmG.Nx[1];
    cmG.is = 1;
    cmG.ie = cmG.Nx[0];
    cmG.mu = NULL;
    cmG.gamma = NULL;  
    cmG.r3imu = NULL;
    cmG.l3imu = NULL;
    cmG.ix1_RBCFun = mG->ix1_RBCFun;
    cmG.ox1_RBCFun = mG->ox1_RBCFun;
    cmG.ix2_RBCFun = mG->ix2_RBCFun;
    cmG.ox2_RBCFun = mG->ox2_RBCFun;
    cmG.rx1_id = mG->rx1_id;
    cmG.lx1_id = mG->lx1_id;
    cmG.rx2_id = mG->rx2_id;
    cmG.lx2_id = mG->lx2_id;

    if ((cmG.R = (RadS ****)calloc_4d_array(1,cmG.Nx[1]+2,cmG.Nx[0]+2,nf,sizeof(RadS))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for multigrid radiation.\n");
    
    if ((cmG.w = (Real **)calloc_2d_array(nmu,ng,sizeof(Real))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for w.\n");

    if ((cmG.r2imu = (Real ******)calloc_6d_array(1,cmG.Nx[0]+2,nf,nmu,ng,2,sizeof(Real))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for r2imu.\n");

    if ((cmG.l2imu = (Real ******)calloc_6d_array(1,cmG.Nx[0]+2,nf,nmu,ng,2,sizeof(Real))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for l2imu.\n");

    if ((cmG.l1imu = (Real *****)calloc_5d_array(1,cmG.Nx[1]+2,nf,nmu,ng,sizeof(Real))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for l1imu.\n");

    if ((cmG.r1imu = (Real *****)calloc_5d_array(1,cmG.Nx[1]+2,nf,nmu,ng,sizeof(Real))) == NULL)
      ath_error("[multig_2d]: Error allocating memory for r1imu.\n");
  
    for(i=0; i<nmu; i++) 
       for(j=0; j<ng; j++) 
	cmG.w[i][j] = mG->w[i][j];

    restrict_2d(mG, &cmG);
    formal_solution_mg_2d_psi(&cmG);
    multigrid_2d(&cmG);
    prolongate_2d(&cmG, mG);

    for(i=0; i < 1; i++) {
      formal_solution_2d(mG);
      bvals_rad(mG);
#ifdef OUT_DIAG
      output_mean_intensity_2d(mG,it0);
      //dumpvtk2d(mG,it0);
      it0++;
#endif
    }
    radgrid_destruct(&cmG);
  }
  return;
}

static void restrict_2d(RadGridS *fmG, RadGridS *cmG) 
{
  int i, j, l, m, n;
  int ifr = 0;
  int isc = cmG->is, iec = cmG->ie;
  int isf = fmG->is, ief = fmG->ie;
  int jsc = cmG->js, jec = cmG->je;
  int jsf = fmG->js, jef = fmG->je;
  int ks = cmG->ks;

  // uses boundary values in restriction
  for (j=jsc; j<=jec; j++) 
    for (i=isc; i<=iec; i++) {  
      cmG->R[ks][j][i][ifr].J = 0.0;
      cmG->R[ks][j][i][ifr].B = 
        0.0625 * (fmG->R[ks][2*j-2][2*i-2][ifr].B   + fmG->R[ks][2*j-2][2*i  ][ifr].B +
                  fmG->R[ks][2*j  ][2*i-2][ifr].B   + fmG->R[ks][2*j  ][2*i  ][ifr].B) +
	0.125  * (fmG->R[ks][2*j-1][2*i  ][ifr].B   + fmG->R[ks][2*j-1][2*i-2][ifr].B +
                  fmG->R[ks][2*j  ][2*i-1][ifr].B   + fmG->R[ks][2*j-2][2*i-1][ifr].B) +
	0.25   *  fmG->R[ks][2*j-1][2*i-1][ifr].B;
      cmG->R[ks][j][i][ifr].S = 
        0.0625 * (fmG->R[ks][2*j-2][2*i-2][ifr].S   + fmG->R[ks][2*j-2][2*i  ][ifr].S +
                  fmG->R[ks][2*j  ][2*i-2][ifr].S   + fmG->R[ks][2*j  ][2*i  ][ifr].S) +
	0.125  * (fmG->R[ks][2*j-1][2*i  ][ifr].S   + fmG->R[ks][2*j-1][2*i-2][ifr].S +
                  fmG->R[ks][2*j  ][2*i-1][ifr].S   + fmG->R[ks][2*j-2][2*i-1][ifr].S) +
	0.25   *  fmG->R[ks][2*j-1][2*i-1][ifr].S;
      cmG->R[ks][j][i][ifr].eps = 
        0.0625 * (fmG->R[ks][2*j-2][2*i-2][ifr].eps + fmG->R[ks][2*j-2][2*i  ][ifr].eps +
                  fmG->R[ks][2*j  ][2*i-2][ifr].eps + fmG->R[ks][2*j  ][2*i  ][ifr].eps) +
	0.125  * (fmG->R[ks][2*j-1][2*i  ][ifr].eps + fmG->R[ks][2*j-1][2*i-2][ifr].eps +
                  fmG->R[ks][2*j  ][2*i-1][ifr].eps + fmG->R[ks][2*j-2][2*i-1][ifr].eps) +
	0.25   *  fmG->R[ks][2*j-1][2*i-1][ifr].eps;      
      cmG->R[ks][j][i][ifr].chi = 
        0.0625 * (fmG->R[ks][2*j-2][2*i-2][ifr].chi + fmG->R[ks][2*j-2][2*i  ][ifr].chi +
                  fmG->R[ks][2*j  ][2*i-2][ifr].chi + fmG->R[ks][2*j  ][2*i  ][ifr].chi) +
	0.125  * (fmG->R[ks][2*j-1][2*i  ][ifr].chi + fmG->R[ks][2*j-1][2*i-2][ifr].chi +
                  fmG->R[ks][2*j  ][2*i-1][ifr].chi + fmG->R[ks][2*j-2][2*i-1][ifr].chi) +
	0.25   *  fmG->R[ks][2*j-1][2*i-1][ifr].chi;
    }
  for (i=isc; i<=iec; i++) {
    cmG->R[ks][jsc-1][i][ifr].B   =       0.5 * fmG->R[ks][jsf-1][2*i-1][ifr].B +
      0.25 * (fmG->R[ks][jsf-1][2*i][ifr].B   + fmG->R[ks][jsf-1][2*i-2][ifr].B);
    cmG->R[ks][jsc-1][i][ifr].S   =       0.5 * fmG->R[ks][jsf-1][2*i-1][ifr].S +
      0.25 * (fmG->R[ks][jsf-1][2*i][ifr].S   + fmG->R[ks][jsf-1][2*i-2][ifr].S);
    cmG->R[ks][jsc-1][i][ifr].eps =       0.5 * fmG->R[ks][jsf-1][2*i-1][ifr].eps +
      0.25 * (fmG->R[ks][jsf-1][2*i][ifr].eps + fmG->R[ks][jsf-1][2*i-2][ifr].eps);
    cmG->R[ks][jsc-1][i][ifr].chi =       0.5 * fmG->R[ks][jsf-1][2*i-1][ifr].chi +
      0.25 * (fmG->R[ks][jsf-1][2*i][ifr].chi + fmG->R[ks][jsf-1][2*i-2][ifr].chi);
    cmG->R[ks][jec+1][i][ifr].B   =       0.5 * fmG->R[ks][jef+1][2*i-1][ifr].B +
      0.25 * (fmG->R[ks][jef+1][2*i][ifr].B   + fmG->R[ks][jef+1][2*i-2][ifr].B);
    cmG->R[ks][jec+1][i][ifr].S   =       0.5 * fmG->R[ks][jef+1][2*i-1][ifr].S +
      0.25 * (fmG->R[ks][jef+1][2*i][ifr].S   + fmG->R[ks][jef+1][2*i-2][ifr].S);
    cmG->R[ks][jec+1][i][ifr].eps =       0.5 * fmG->R[ks][jef+1][2*i-1][ifr].eps +
      0.25 * (fmG->R[ks][jef+1][2*i][ifr].eps + fmG->R[ks][jef+1][2*i-2][ifr].eps);
    cmG->R[ks][jec+1][i][ifr].chi =       0.5 * fmG->R[ks][jef+1][2*i-1][ifr].chi +
      0.25 * (fmG->R[ks][jef+1][2*i][ifr].chi + fmG->R[ks][jef+1][2*i-2][ifr].chi);
    for(l=0; l<fmG->nmu; l++) 
      for(m=0; m<fmG->ng; m++) {
	cmG->r2imu[ks][i][ifr][l][m][n] =       0.5 * fmG->r2imu[ks][2*i-1][ifr][l][m][n] +
	  0.25 * (fmG->r2imu[ks][2*i][ifr][l][m][n] + fmG->r2imu[ks][2*i-2][ifr][l][m][n]);
	cmG->l2imu[ks][i][ifr][l][m][n] =       0.5 * fmG->l2imu[ks][2*i-1][ifr][l][m][n] +
	  0.25 * (fmG->l2imu[ks][2*i][ifr][l][m][n] + fmG->l2imu[ks][2*i-2][ifr][l][m][n]);
     }
  }
  for (j=cmG->is; j<=cmG->ie; j++) {
   cmG->R[ks][j][isc-1][ifr].B   =       0.5 * fmG->R[ks][2*j-1][isf-1][ifr].B +
     0.25 * (fmG->R[ks][2*j][isf-1][ifr].B   + fmG->R[ks][2*j-2][isf-1][ifr].B);
   cmG->R[ks][j][isc-1][ifr].S   =       0.5 * fmG->R[ks][2*j-1][isf-1][ifr].S +
     0.25 * (fmG->R[ks][2*j][isf-1][ifr].S   + fmG->R[ks][2*j-2][isf-1][ifr].S);
   cmG->R[ks][j][isc-1][ifr].eps =       0.5 * fmG->R[ks][2*j-1][isf-1][ifr].eps +
     0.25 * (fmG->R[ks][2*j][isf-1][ifr].eps + fmG->R[ks][2*j-2][isf-1][ifr].eps);
   cmG->R[ks][j][isc-1][ifr].chi =       0.5 * fmG->R[ks][2*j-1][isf-1][ifr].chi +
     0.25 * (fmG->R[ks][2*j][isf-1][ifr].chi + fmG->R[ks][2*j-2][isf-1][ifr].chi);

   cmG->R[ks][j][iec+1][ifr].B   =       0.5 * fmG->R[ks][2*j-1][ief+1][ifr].B +
     0.25 * (fmG->R[ks][2*j][ief+1][ifr].B   + fmG->R[ks][2*j-2][ief+1][ifr].B);
   cmG->R[ks][j][iec+1][ifr].S   =       0.5 * fmG->R[ks][2*j-1][ief+1][ifr].S +
     0.25 * (fmG->R[ks][2*j][ief+1][ifr].S   + fmG->R[ks][2*j-2][ief+1][ifr].S);
   cmG->R[ks][j][iec+1][ifr].eps =       0.5 * fmG->R[ks][2*j-1][ief+1][ifr].eps +
     0.25 * (fmG->R[ks][2*j][ief+1][ifr].eps + fmG->R[ks][2*j-2][ief+1][ifr].eps);
   cmG->R[ks][j][iec+1][ifr].chi =       0.5 * fmG->R[ks][2*j-1][ief+1][ifr].chi +
     0.25 * (fmG->R[ks][2*j][ief+1][ifr].chi + fmG->R[ks][2*j-2][ief+1][ifr].chi);

   for(l=0; l<fmG->nmu; l++) 
     for(m=0; m<fmG->ng; m++) {
       cmG->r1imu[ks][j][ifr][l][m] =       0.5 * fmG->r1imu[ks][2*j-1][ifr][l][m] +
	 0.25 * (fmG->r1imu[ks][2*j][ifr][l][m] + fmG->r1imu[ks][2*j-2][ifr][l][m]);
       cmG->l1imu[ks][j][ifr][l][m] =       0.5 * fmG->l1imu[ks][2*j-1][ifr][l][m] +
	 0.25 * (fmG->l1imu[ks][2*j][ifr][l][m] + fmG->l1imu[ks][2*j-2][ifr][l][m]);
     }
  }

  img++;
  return;
}

static void prolongate_2d(RadGridS *cmG, RadGridS *fmG) 
{
  int i, j, l, m;
  int ifr = 0;
  int isc = cmG->is, iec = cmG->ie;
  int isf = fmG->is, ief = fmG->ie;
  int jsc = cmG->js, jec = cmG->je;
  int jsf = fmG->js, jef = fmG->je;
  int ks = cmG->ks;

  for (j=jsf; j<=jef; j+=2) {
    for (i=isf; i<=ief; i+=2) 
      fmG->R[ks][j][i][ifr].S = cmG->R[ks][j/2+1][i/2+1][ifr].S;
    for (i=isf+1; i<=ief; i+=2) 
      fmG->R[ks][j][i][ifr].S = 0.5 * 
	(cmG->R[ks][j/2+1][i/2][ifr].S + cmG->R[ks][j/2+1][i/2+1][ifr].S);
  }
  for (j=jsf+1; j<=jef; j+=2) {
    for (i=fmG->is; i<=ief; i+=2) 
      fmG->R[ks][j][i][ifr].S = 0.5 * 
	(cmG->R[ks][j/2][i/2+1][ifr].S + cmG->R[ks][j/2+1][i/2+1][ifr].S);
    for (i=isf+1; i<=ief; i+=2) 
      fmG->R[ks][j][i][ifr].S = 0.25 * 
	(cmG->R[ks][j/2  ][i/2  ][ifr].S + cmG->R[ks][j/2+1][i/2  ][ifr].S +
	 cmG->R[ks][j/2  ][i/2+1][ifr].S + cmG->R[ks][j/2+1][i/2+1][ifr].S);
  }
  for (j=jsf; j<=jef; j+=2) {
    fmG->R[ks][j][isf-1][ifr].S = cmG->R[ks][j/2+1][isc-1][ifr].S;
    fmG->R[ks][j][ief+1][ifr].S = cmG->R[ks][j/2+1][iec+1][ifr].S;
  }
  for (j=jsf+1; j<=jef; j+=2) {
    fmG->R[ks][j][isf-1][ifr].S = 
      0.5 * (cmG->R[ks][j/2][isc-1][ifr].S + cmG->R[ks][j/2+1][isc-1][ifr].S);
    fmG->R[ks][j][ief+1][ifr].S = 
      0.5 * (cmG->R[ks][j/2][iec+1][ifr].S + cmG->R[ks][j/2+1][iec+1][ifr].S);  
  }
  for (j=jsf; j<=jef; j+=2) 
    for(l=0; l<fmG->nmu; l++) 
      for(m=0; m<fmG->ng; m++) {
	fmG->r1imu[ks][j][ifr][l][m] = cmG->r1imu[ks][j/2+1][ifr][l][m];
	fmG->l1imu[ks][j][ifr][l][m] = cmG->l1imu[ks][j/2+1][ifr][l][m];
      }   
  for (j=jsf+1; j<=jef; j+=2) 
    for(l=0; l<fmG->nmu; l++) 
      for(m=0; m<fmG->ng; m++) {
	fmG->r1imu[ks][j][ifr][l][m] = 
	  0.5 * (cmG->r1imu[ks][j/2][ifr][l][m] + cmG->r1imu[ks][j/2+1][ifr][l][m]);
	fmG->l1imu[ks][j][ifr][l][m] = 
	  0.5 * (cmG->l1imu[ks][j/2][ifr][l][m] + cmG->l1imu[ks][j/2+1][ifr][l][m]);
      }
  img--;
  return;
}

static void formal_solution_mg_2d_psi(RadGridS *pRG)
{

  int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1];
  int nmu = pRG->nmu, ng = pRG->ng;
  int nmu2 = pRG->nmu/2, ng2 = pRG->ng/2;
  int ifr, nf = pRG->nf;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks; 
  Real dx = pRG->dx1, dy = pRG->dx2;
  int i, j, l, m;
  int sy, sx;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;
  Real am, am1, bm, bm1;
  Real *******psi = NULL;
  Real *mu1 = NULL, *gamma1 =NULL, **am0 = NULL;
  Real *****psiint = NULL;

#if defined(JACOBI)
  jacobi_pass_pointers_to_mg_2d(&psi, &mu1, &gamma1, &am0);
#elif defined(GAUSSEID)
  gs_pass_pointers_to_mg_2d(&psi, &mu1, &gamma1, &am0, &psiint);
  if ((psiint[img] = (Real ****)calloc_4d_array(nx2+2,nx1+2,nf,4,sizeof(Real))) == NULL) 
    goto on_error;
#endif

  if ((psi[img] = (Real ******)calloc_6d_array(nx2+2,nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

/* ---------  Compute weights and exponentials ------------------------- */
  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++) 
      for(ifr=0; ifr<nf; ifr++) {
	chi1 = pRG->R[ks][j][i][ifr].chi;
	for(l=0; l<nmu; l++) {
	  if(l < nmu2) sy = -1; else sy = 1;
	  for(m=0; m<ng; m++) {
	    if(m < ng2) sx = -1; else sx = 1;
	    am = am0[l][m];
	    if (am <= 1.0) {
	      am1 = 1.0 - am;
	      chi0 = am  * pRG->R[ks][j-sy][i-sx][ifr].chi + 
                     am1 * pRG->R[ks][j-sy][i   ][ifr].chi;
	      chi2 = am  * pRG->R[ks][j+sy][i+sx][ifr].chi + 
                     am1 * pRG->R[ks][j+sy][i   ][ifr].chi;
	      dtaum = 0.5 * (chi0 + chi1) * dy * mu1[l]; 
	      dtaup = 0.5 * (chi2 + chi1) * dy * mu1[l]; 
	    } else {
	      bm = 1.0 / am;
	      bm1 = 1.0 - bm;
	      chi0 = bm  * pRG->R[ks][j-sy][i-sx][ifr].chi + 
	             bm1 * pRG->R[ks][j   ][i-sx][ifr].chi;
	      chi2 = bm  * pRG->R[ks][j+sy][i+sx][ifr].chi +
                     bm1 * pRG->R[ks][j   ][i+sx][ifr].chi;
	      dtaum = 0.5 * (chi0 + chi1) * dx * gamma1[m]; 
	      dtaup = 0.5 * (chi2 + chi1) * dx * gamma1[m]; 
	    }

	    get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);
	    psi[img][j][i][ifr][l][m][0] = edtau;
	    psi[img][j][i][ifr][l][m][1] = a0;
	    psi[img][j][i][ifr][l][m][2] = a1;
	    psi[img][j][i][ifr][l][m][3] = a2;
	    pRG->R[ks][j][i][ifr].lamstr += pRG->w[l][m] * a1;

#ifdef GAUSSEID
	    if (sy == 1) 
	      if (sx == 1) 
		if (am <= 1.0) {
		  psiint[img][j][i][ifr][1] += am  * pRG->w[l][m] * a2; 
		  psiint[img][j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		} else {
		  psiint[img][j][i][ifr][0] += bm1 * pRG->w[l][m] * a2; 
		  psiint[img][j][i][ifr][1] += bm  * pRG->w[l][m] * a2; 
		}
	      else
		if (am <= 1.0) {
		  psiint[img][j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		  psiint[img][j][i][ifr][3] += am  * pRG->w[l][m] * a2; 
		} else 
		  psiint[img][j][i][ifr][3] += bm  * pRG->w[l][m] * a2;
#endif
	  }
	}
      }

  return;

  on_error:
  formal_solution_2d_destruct();
  ath_error("[formal_solution_mg_2d_psi]: Error allocating memory\n");
  return;

}

#endif /* RAD_MULTIG */
#endif /* RADIATION */
