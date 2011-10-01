#include "../copyright.h"
/*==============================================================================
 * FILE: gausseid_2d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 2D grid using Gauss-Seidel method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *          SOR is incompletely implemented for 2D and 3D grids with periodic
 *          boundary conditions, and generally not suitable.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution_2d.c()
 *   formal_solution_2d_destruct()
 *   formal_solution_2d_init()
 *   gs_pass_pointers_to_mg_2d()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define INTERP_2D

#ifdef RADIATION_TRANSFER
#ifdef GAUSSEID

static Real ***lamstr = NULL, **psi0= NULL;
static Real ******psi = NULL, ****psiint = NULL;
static Real *****imu1 = NULL, *****imu2 = NULL;
static Real **muinv = NULL, *am0 = NULL, ***mu2 = NULL;
static Real dSrmx;
static int ntot;
static int iter;
static int svwght;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_2d()     - computes a single sweep in one direction (right or left)
 *   update_sfunc() - updates source function after compute of mean intensity
 *   set_bvals_imu_y()      - set imu array at vertical boundary
 *   set_bvals_imu_y_j()    - set imu array at horizontal boundary
 *   update_bvals_imu_y()   - update outgoing radiation at vertical boundary
 *   update_bvals_imu_y_j() - update outgoing radiation at horizontal boundary
 *============================================================================*/

static void update_sfunc(RadS *R, Real *dS, Real lamstr);
static void update_cell(RadGridS *pRG, Real *****imuo, int ifr, int k, int j, int i, 
			int l, int m);
static void sweep_2d_forward(RadGridS *pRG);
static void sweep_2d_backward(RadGridS *pRG);

void formal_solution_2d(RadGridS *pRG, Real *dSrmax, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie; 
  int js = pRG->js, je = pRG->je; 
  int ks = pRG->ks; 
  int nf = pRG->nf;  

  for(ifr=0; ifr<nf; ifr++) 
    for(m=0; m<pRG->nang; m++) {
      pRG->l2imu[ifr][ks][is-1][0][m] = pRG->l1imu[ifr][ks][js-1][0][m];
      pRG->l2imu[ifr][ks][is-1][0][m] = pRG->l1imu[ifr][ks][js-1][0][m];
      pRG->l2imu[ifr][ks][ie+1][1][m] = pRG->r1imu[ifr][ks][js-1][1][m];
      pRG->l2imu[ifr][ks][ie+1][1][m] = pRG->r1imu[ifr][ks][js-1][1][m];

      pRG->r2imu[ifr][ks][is-1][2][m] = pRG->l1imu[ifr][ks][je+1][2][m];
      pRG->r2imu[ifr][ks][is-1][2][m] = pRG->l1imu[ifr][ks][je+1][2][m];
      pRG->r2imu[ifr][ks][ie+1][3][m] = pRG->r1imu[ifr][ks][je+1][3][m];
      pRG->r2imu[ifr][ks][ie+1][3][m] = pRG->r1imu[ifr][ks][je+1][3][m];
    }

/* Initialize dSrmx */
  dSrmx = 0.0;

/* initialize mean intensities at all depths to zero */
  for(j=js-1; j<+je+1; j++)
    for(i=is-1; i<=ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++) {
	pRG->R[ks][j][i][ifr].J = 0.0;
	pRG->R[ks][j][i][ifr].K[0] = 0.0;
	pRG->R[ks][j][i][ifr].K[1] = 0.0;
	pRG->R[ks][j][i][ifr].K[2] = 0.0;
	if (svwght == 0) {
	  lamstr[j][i][ifr] = 0.0;
	  for(l=0; l<4; l++)
	    psiint[j][i][ifr][l] = 0.0;
	}
      }

/* Compute formal solution for all upward and rightward going rays in 
 * each gridzone */
  sweep_2d_forward(pRG);

/* Compute formal solution for all downward and leftward going rays in 
 * each vertical gridzone */
  sweep_2d_backward(pRG);

/* Return maximum relative change to test convergence*/
  (*dSrmax) = dSrmx;

  return;
}

static void update_cell(RadGridS *pRG, Real *****imuo, int ifr, int k, int j, int i,
			int l, int m)
{

  int im, ip, jm, jp, imm, imp;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real am, am1, bm, bm1;
  Real w0, w1, w2;
  Real maxint, minint;
  Real dx = pRG->dx1, dy = pRG->dx2;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2, wa2;

/* initialize stencil base on quadrant*/  
  if(l == 0) {
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
    imm = 1; imp = 0;
  } else if (l == 1) {
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
    imm = 0; imp = 1;
  } else if (l == 2) {
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
    imm = 0; imp = 1;
  } else {
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
    imm = 1; imp = 0;
  }

 if(svwght == 0) chi1 = pRG->R[k][j][i][ifr].chi;
/* --------- Interpolate intensity and source functions at endpoints --------- 
 * --------- of characteristics                                      --------- */
  am = am0[m];
  if (am <= 1.0) {
    am1 = 1.0 - am;
    /* Use linear interpolation for source functions */
    S0 = am  * pRG->R[k][jm][im][ifr].S +
         am1 * pRG->R[k][jm][i ][ifr].S;
    S2 = am  * pRG->R[k][jp][ip][ifr].S +
         am1 * pRG->R[k][jp][i ][ifr].S;

#ifdef INTERP_2D /* Use parabolic interpolation for intensity */
    w0 = 0.5 * am * (1.0 + am);
    w1 = am1 * (1.0 + am);
    w2 = -0.5 * am * am1;
    imu0 = w0 * imuo[ifr][im][l][m][imm] + w1 * imuo[ifr][i][l][m][0] +
           w2 * imuo[ifr][ip][l][m][imp];
    maxint = MAX(imuo[ifr][im][l][m][imm],imuo[ifr][i][l][m][0]);
    minint = MIN(imuo[ifr][im][l][m][imm],imuo[ifr][i][l][m][0]);
    if(imu0 > maxint) imu0 = maxint;
    if(imu0 < minint) imu0 = minint;
#else     /* Use linear interpolation for intensity */
    imu0 = am  * imuo[ifr][im][l][m][imm] + am1 * imuo[ifr][i][l][m][0];
#endif
  } else {
    bm = 1.0 / am;
    bm1 = 1.0 - bm;

    /* Use linear interpolation for source functions */
    S0 = bm  * pRG->R[k][jm][im][ifr].S +
         bm1 * pRG->R[k][j ][im][ifr].S;
    S2 = bm  * pRG->R[k][jp][ip][ifr].S +
         bm1 * pRG->R[k][j ][ip][ifr].S;

#ifdef INTERP_2D /* Use parabolic interpolation for intensity */
    w0 = 0.5 * bm1 * (1.0 + bm1);  /* Modify to compute only once? */
    w1 = bm * (1.0 + bm1);
    w2 = -0.5 * bm * bm1;

    imu0 = w0 * imuo[ifr][jm][l][m][imm] + w1 * imuo[ifr][j][l][m][0] +
           w2 * imuo[ifr][jp][l][m][imp];
    maxint = MAX(imuo[ifr][jm][l][m][imm],imuo[ifr][j][l][m][0]);
    minint = MIN(imuo[ifr][jm][l][m][imm],imuo[ifr][j][l][m][0]);
    if(imu0 > maxint) imu0 = maxint;
    if(imu0 < minint) imu0 = minint;
#else     /* Use linear interpolation for intensity */
    imu0 = bm  * imuo[ifr][jm][l][m][imm] + bm1 * imuo[ifr][j][l][m][0];
#endif
  }
/* ---------  compute intensity at grid center and add to mean intensity ------- */
  if(svwght == 1) {
    imu = psi[j][i][ifr][l][m][1] * S0 +
          psi[j][i][ifr][l][m][2] * pRG->R[k][j][i][ifr].S +
          psi[j][i][ifr][l][m][3] * S2;	
    if (imu < 0.0) imu=0.0;
    imu += psi[j][i][ifr][l][m][0] * imu0;
  } else {
    if (am <= 1.0) {
      chi0 = am  * pRG->R[k][jm][im][ifr].chi + 
	     am1 * pRG->R[k][jm][i ][ifr].chi;
      chi2 = am  * pRG->R[k][jp][ip][ifr].chi + 
	     am1 * pRG->R[k][jp][i ][ifr].chi;
      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dy * muinv[m][1]; 
      dtaup *= dy * muinv[m][1]; 
    } else {
      chi0 = bm  * pRG->R[k][jm][im][ifr].chi + 
	     bm1 * pRG->R[k][j ][im][ifr].chi;
      chi2 = bm  * pRG->R[k][jp][ip][ifr].chi +
	     bm1 * pRG->R[k][j ][ip][ifr].chi;
      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dx * muinv[m][0]; 
      dtaup *= dx * muinv[m][0]; 
    }
    interp_quad_source_slope_lim(dtaum, dtaup, &edtau, &a0, &a1, &a2,
    S0, pRG->R[k][j][i][ifr].S, S2);

    imu = a0 * S0 + a1 * pRG->R[k][j][i][ifr].S + a2 * S2;
    imu += edtau * imu0;
    lamstr[j][i][ifr] += pRG->wmu[m] * a1;
    psi0[l][m] = a1;

    if (l == 0) {
      if (am <= 1.0) {
	psiint[j][i][ifr][1] += am  * pRG->wmu[m] * a2; 
	psiint[j][i][ifr][2] += am1 * pRG->wmu[m] * a2; 
      } else {
	psiint[j][i][ifr][0] += bm1 * pRG->wmu[m] * a2; 
	psiint[j][i][ifr][1] += bm  * pRG->wmu[m] * a2; 
      }
    } else if (l == 1) {
      if (am <= 1.0) {
	psiint[j][i][ifr][2] += am1 * pRG->wmu[m] * a2; 
	psiint[j][i][ifr][3] += am  * pRG->wmu[m] * a0; 
      }
    } else if (l == 2) {
      if (am <= 1.0) {
	psiint[j][i][ifr][3] += am * pRG->wmu[m] * a2; 
      } else {
	psiint[j][i][ifr][0] += bm1 * pRG->wmu[m] * a2; 
	psiint[j][i][ifr][3] += bm  * pRG->wmu[m] * a2;
      }
    }
  }

  /* Add to radiation moments and save for next iteration */
  wimu = pRG->wmu[m] * imu;
  pRG->R[k][j][i][ifr].J += wimu;
  pRG->R[k][j][i][ifr].K[0] += mu2[l][m][0] * wimu;
  pRG->R[k][j][i][ifr].K[1] += mu2[l][m][1] * wimu;
  pRG->R[k][j][i][ifr].K[2] += mu2[l][m][2] * wimu;
  /* Update intensity workspace */
  if (am <= 1.0) {
    imuo[ifr][i][l][m][1] = imuo[ifr][i][l][m][0];
    imuo[ifr][i][l][m][0] = imu;
  } else {
    imuo[ifr][j][l][m][1] = imuo[ifr][j][l][m][0];
    imuo[ifr][j][l][m][0] = imu;
  }

  return;
}

static void sweep_2d_forward(RadGridS *pRG)
{
  int ifr, l, m, n;
  int itot;
  int i, il, iu;
  int j, jl, ju;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;   
  int nf = pRG->nf, nang = pRG->nang;


  for(ifr=0; ifr<nf; ifr++)  {
/* Account for ix2 boundary condition */
    for(i=is-1; i<=ie+1; i++)
      for(m=0; m<nang; m++) 
	for(n=0; n<2; n++) {
	  imu2[ifr][i][0][m][n] = pRG->l2imu[ifr][ks][i][0][m];
	  imu2[ifr][i][1][m][n] = pRG->l2imu[ifr][ks][i][1][m];
	}
/* Account for ix1 boundary condition */
    for(j=js-1; j<=je+1; j++)
      for(m=0; m<nang; m++)
	for(n=0; n<2; n++) {
	  imu1[ifr][j][0][m][n] = pRG->l1imu[ifr][ks][j][0][m];
	  imu1[ifr][j][2][m][n] = pRG->l1imu[ifr][ks][j][2][m];
	}
/* Peform diagonal sweep through 2D grid */
    il = is; iu = is;  jl = js; ju = js;
    for(itot=0; itot < ntot; itot++) {
      j = jl;  i = iu;
      while ((j <= ju) && (i >= il)) {
	for(m=0; m<nang; m++)
	  if(am0[m] <= 1.0) {
/* modify imu2 to include left and right boundary intensity */
	    if(i == is) {
	      imu2[ifr][i-1][0][m][1] = pRG->l1imu[ifr][ks][j-1][0][m];
	      imu2[ifr][i-1][1][m][1] = pRG->l1imu[ifr][ks][j-1][1][m];
	    }
	    if(i == ie) {
	      imu2[ifr][i+1][0][m][0] = pRG->r1imu[ifr][ks][j-1][0][m];
	      imu2[ifr][i+1][1][m][0] = pRG->r1imu[ifr][ks][j-1][1][m];
	    }
/* compute imu in current cell */
	    update_cell(pRG,imu2,ifr,ks,j,i,0,m);
	    update_cell(pRG,imu2,ifr,ks,j,i,1,m);
/* update left and right boundary intensity */
	    if(i == is)
	      pRG->l1imu[ifr][ks][j][1][m] = imu2[ifr][i][1][m][0];	    
	    if(i == ie)
	      pRG->r1imu[ifr][ks][j][0][m] = imu2[ifr][i][0][m][0];	    
	  }
	j++; i--;
      }
      j = ju;  i = il;
      while ((j >= jl) && (i <= iu)) {
	for(m=0; m<nang; m++)
	  if(am0[m] > 1.0) {
/* modify imu1 to include top and boundar boundary intensity */
	    /* if(j == js) { */
	    if((j == js) && (i != is)) {
	      imu1[ifr][j-1][0][m][1] = pRG->l2imu[ifr][ks][i-1][0][m];
	      imu1[ifr][j-1][2][m][1] = pRG->l2imu[ifr][ks][i-1][2][m];
	    }
	    /* if(j == je) { */
	    if((j == je) && (i != is)) {
	      imu1[ifr][j+1][0][m][0] = pRG->r2imu[ifr][ks][i-1][0][m];
	      imu1[ifr][j+1][2][m][0] = pRG->r2imu[ifr][ks][i-1][2][m];
	    }
/* compute imu in current cell */
	    update_cell(pRG,imu1,ifr,ks,j,i,0,m);
	    update_cell(pRG,imu1,ifr,ks,j,i,2,m);
/* update top and bottom boundary intensity */
	    if(j == js)
	      pRG->l2imu[ifr][ks][i][2][m] = imu1[ifr][j][2][m][0];
	    if(j == je)
	      pRG->r2imu[ifr][ks][i][0][m] = imu1[ifr][j][0][m][0];
	  }
	j--; i++;
      }
      if(iu == ie) jl++;
      if(ju == je) il++;      
      if(iu < ie) iu++;
      if(ju < je) ju++;
    }

/* Update ox2 boundary condition */
    for(i=is; i<=ie; i++)
      for(m=0; m<nang; m++)
	if(am0[m] <= 1.0) {
	  pRG->r2imu[ifr][ks][i][0][m] = imu2[ifr][i][0][m][0];
	  pRG->r2imu[ifr][ks][i][1][m] = imu2[ifr][i][1][m][0];
	}

/* Update ox1 boundary condition */
    for(j=js; j<=je; j++)
      for(m=0; m<nang; m++) 
	if(am0[m] > 1.0) {
	  pRG->r1imu[ifr][ks][j][0][m] = imu1[ifr][j][0][m][0];
	  pRG->r1imu[ifr][ks][j][2][m] = imu1[ifr][j][2][m][0];
	}    
  }

  return;
}

static void sweep_2d_backward(RadGridS *pRG)
{
  int ifr, l, m, n;
  int itot;
  int i, il, iu;
  int j, jl, ju;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;   
  int nf = pRG->nf, nang = pRG->nang;
  Real dS;

  iter=0;

  for(ifr=0; ifr<nf; ifr++)  {
/* Account for ox2 boundary condition */
    for(i=is-1; i<=ie+1; i++)
      for(m=0; m<nang; m++) {
	imu2[ifr][i][2][m][0] = pRG->r2imu[ifr][ks][i][2][m];
	imu2[ifr][i][3][m][0] = pRG->r2imu[ifr][ks][i][3][m];
      }
/* Account for ox1 boundary condition */
    for(j=js-1; j<=je+1; j++)
      for(m=0; m<nang; m++) {
	imu1[ifr][j][1][m][0] = pRG->r1imu[ifr][ks][j][1][m];
	imu1[ifr][j][3][m][0] = pRG->r1imu[ifr][ks][j][3][m];
      }
/* Peform diagonal sweep through 2D grid */
    il = ie; iu = ie;  jl = je; ju = je;
    for(itot=0; itot < ntot; itot++) {
      j = ju;  i = il;
      while ((j >= jl) && (i <= iu)) {
	for(m=0; m<nang; m++)
	  if(am0[m] <= 1.0) {	    
/* modify imu2 to include left and right boundary intensity */
	    if(i == is) {
	      imu2[ifr][i-1][2][m][0] = pRG->l1imu[ifr][ks][j+1][2][m];
	      imu2[ifr][i-1][3][m][0] = pRG->l1imu[ifr][ks][j+1][3][m];
	    }
	    if(i == ie) {
	      imu2[ifr][i+1][2][m][1] = pRG->r1imu[ifr][ks][j+1][2][m];
	      imu2[ifr][i+1][3][m][1] = pRG->r1imu[ifr][ks][j+1][3][m];
	    }
/* compute imu in current cell */
	    update_cell(pRG,imu2,ifr,ks,j,i,2,m);
	    update_cell(pRG,imu2,ifr,ks,j,i,3,m);
/* update left and right boundary intensity */
	    if(i == is)
	      pRG->l1imu[ifr][ks][j][3][m] = imu2[ifr][i][3][m][0];	    
	    if(i == ie)
	      pRG->r1imu[ifr][ks][j][2][m] = imu2[ifr][i][2][m][0];	    
	  }
	j--; i++;
      }
      j = jl;  i = iu;
      while ((j <= ju) && (i >= il)) {
	for(m=0; m<nang; m++) {
	  if(am0[m] > 1.0) {
/* modify imu1 to include top and boundar boundary intensity */
	    /* if(j == js) { */
	    if((j == js) && (i != ie)) {
	      imu1[ifr][j-1][1][m][0] = pRG->l2imu[ifr][ks][i+1][1][m];
	      imu1[ifr][j-1][3][m][0] = pRG->l2imu[ifr][ks][i+1][3][m];
	    }
	    /* if(j == je) { */
	    if((j == je) && (i != ie)) {
	      imu1[ifr][j+1][1][m][1] = pRG->r2imu[ifr][ks][i+1][1][m];
	      imu1[ifr][j+1][3][m][1] = pRG->r2imu[ifr][ks][i+1][3][m];
	    }
/* compute imu in current cell */
	    update_cell(pRG,imu1,ifr,ks,j,i,1,m);
	    update_cell(pRG,imu1,ifr,ks,j,i,3,m);
/* update top and bottom boundary intensity */
	    if(j == js)
	      pRG->l2imu[ifr][ks][i][3][m] = imu1[ifr][j][3][m][0];
	    if(j == je)
	      pRG->r2imu[ifr][ks][i][1][m] = imu1[ifr][j][1][m][0];
	  }
	}
/* update source function and relevant quantities for GS integration */
	update_sfunc(&(pRG->R[ks][j][i]), &dS, lamstr[j][i][ifr]);
/* correct intensties */
	if (svwght == 1) {
	  for(m=0; m<nang; m++) {
	    if(am0[m] > 1.0) {
	      imu1[ifr][j][1][m][0] += dS * pRG->wmu[m] *
	                               psi[j][i][ifr][1][m][2];
	      imu1[ifr][j][3][m][0] += dS * pRG->wmu[m] *
	                               psi[j][i][ifr][3][m][2];
	    } else {
	      imu2[ifr][i][2][m][0] += dS * pRG->wmu[m] *
	                               psi[j][i][ifr][2][m][2];
	      imu2[ifr][i][3][m][0] += dS * pRG->wmu[m] *
	                               psi[j][i][ifr][3][m][2];
	    }
	  }
	} else {
	  for(m=0; m<nang; m++) {
	    if(am0[m] > 1.0) {
	      imu1[ifr][j][1][m][0] += dS * pRG->wmu[m] *
	                               psi0[1][m];
	      imu1[ifr][j][3][m][0] += dS * pRG->wmu[m] *
	                               psi0[3][m];
	    } else {
	      imu2[ifr][i][2][m][0] += dS * pRG->wmu[m] *
	                               psi0[2][m];
	      imu2[ifr][i][3][m][0] += dS * pRG->wmu[m] *
	                               psi0[3][m];
	    }
	  }
	} 
/* Correct J w/ updated S from "new" neighbors, but not in ghostzones */	
	if(i != is) {
	  pRG->R[ks][j][i-1][ifr].J += dS * psiint[j][i-1][ifr][0];
	  if (j != js)
	    pRG->R[ks][j-1][i-1][ifr].J += dS * psiint[j-1][i-1][ifr][1];
	  if (j != je)
	    pRG->R[ks][j+1][i-1][ifr].J += dS * psiint[j+1][i-1][ifr][3];
	} 
	if(j != js)
	  pRG->R[ks][j-1][i  ][ifr].J += dS * psiint[j-1][i  ][ifr][2];

	j++; i--;
      }
      if(il == is) ju--;
      if(jl == js) iu--;  
      if(il > is) il--;
      if(jl > js) jl--;
    }    
/* Update ix2 boundary condition */
    for(i=is; i<=ie; i++)
      for(m=0; m<nang; m++)
	if(am0[m] <= 1.0) {
	  pRG->l2imu[ifr][ks][i][2][m] = imu2[ifr][i][2][m][0];
	  pRG->l2imu[ifr][ks][i][3][m] = imu2[ifr][i][3][m][0];
	}
/* Update ix1 boundary condition */
    for(j=js; j<=je; j++)
      for(m=0; m<nang; m++) 
	if(am0[m] > 1.0) {
	  pRG->l1imu[ifr][ks][j][1][m] = imu1[ifr][j][1][m][0];
	  pRG->l1imu[ifr][ks][j][3][m] = imu1[ifr][j][3][m][0];
	}    
  }

  return;
}

static void update_sfunc(RadS *R, Real *dS, Real lamstr)
{
  Real snew, dSr;
  
  snew = (1.0 - R->eps) * R->J + R->eps * R->B;

  (*dS) = (snew - R->S) / (1.0 - (1.0 - R->eps) * lamstr);
  dSr = fabs((*dS) / R->S);
  R->S += (*dS);
  if (dSr > dSrmx) dSrmx = dSr; 
  return;
}

void formal_solution_2d_destruct(void)
{
  if (psi    != NULL) free_6d_array(psi);
  if (psiint != NULL) free_4d_array(psiint);
  if (psi0   != NULL) free_2d_array(psi0);
  if (lamstr != NULL) free_3d_array(lamstr);
  if (imu1   != NULL) free_5d_array(imu1);
  if (imu2   != NULL) free_5d_array(imu2);
  if (muinv  != NULL) free_2d_array(muinv);
  if (am0    != NULL) free_1d_array(am0);
  if (mu2    != NULL) free_3d_array(mu2);

  return;
}

void formal_solution_2d_init(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
  int nf = pRG->nf, nang = pRG->nang;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks; 
  Real dx = pRG->dx1, dy = pRG->dx2;
  int ifr, i, j, l, m;
  int sy, sx;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;
  Real am, am1, bm, bm1;

 svwght = par_geti("radiation","svwght");

  ntot = je + ie - (js + is) + 1;

/* Initialize variables for sor */

  if ((imu1 = (Real *****)calloc_5d_array(nf,nx2+2,4,nang,2,sizeof(Real))) == NULL)
    goto on_error;

  if ((imu2 = (Real *****)calloc_5d_array(nf,nx1+2,4,nang,2,sizeof(Real))) == NULL)
    goto on_error;

  if ((muinv = (Real **)calloc_2d_array(nang,2,sizeof(Real))) == NULL)
    goto on_error;

  if ((am0 = (Real *)calloc_1d_array(nang,sizeof(Real))) == NULL)
    goto on_error;

  if ((mu2 = (Real ***)calloc_3d_array(4,nang,3,sizeof(Real))) == NULL)
    goto on_error;

  for(i=0; i<nang; i++)  
    for(j=0; j<2; j++) 
      muinv[i][j] = fabs(1.0 / pRG->mu[0][i][j]);

  for(i=0; i<nang; i++)
    am0[i]   = fabs(dy * muinv[i][1] / (dx * muinv[i][0]));

  for(i=0; i<4; i++) 
    for(j=0; j<nang; j++)  {
      mu2[i][j][0] = pRG->mu[i][j][0] * pRG->mu[i][j][0];
      mu2[i][j][1] = pRG->mu[i][j][0] * pRG->mu[i][j][1];
      mu2[i][j][2] = pRG->mu[i][j][1] * pRG->mu[i][j][1];
    }

  if ((lamstr = (Real ***)calloc_3d_array(nx2+1,nx1+2,nf,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psi0 = (Real **)calloc_2d_array(4,nang,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psiint = (Real ****)calloc_4d_array(nx2+2,nx1+2,nf,5,sizeof(Real))) == NULL) 
    goto on_error;

  if(svwght == 1) {

    if ((psi = (Real ******)calloc_6d_array(nx2+2,nx1+2,nf,4,nang,4,sizeof(Real))) == NULL) 
      goto on_error;
/* ---------  Compute weights and exponentials ------------------------- */
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) 
	for(ifr=0; ifr<nf; ifr++) {
	  chi1 = pRG->R[ks][j][i][ifr].chi;
	  for(l=0; l<4; l++) {
	    if(l < 2) sy = 1; else sy = -1;
	    if((l == 0) || (l == 2)) sx = 1; else sx = -1;
	    for(m=0; m<nang; m++) {
	      am = am0[m];
	      if (am <= 1.0) {
		am1 = 1.0 - am;
		chi0 = am  * pRG->R[ks][j-sy][i-sx][ifr].chi + 
                       am1 * pRG->R[ks][j-sy][i   ][ifr].chi;
		chi2 = am  * pRG->R[ks][j+sy][i+sx][ifr].chi + 
                       am1 * pRG->R[ks][j+sy][i   ][ifr].chi;
		dtaum = 0.5 * (chi0 + chi1) * dy * muinv[m][1]; 
		dtaup = 0.5 * (chi2 + chi1) * dy * muinv[m][1]; 
	      } else {
		bm = 1.0 / am;
		bm1 = 1.0 - bm;
		chi0 = bm  * pRG->R[ks][j-sy][i-sx][ifr].chi + 
	               bm1 * pRG->R[ks][j   ][i-sx][ifr].chi;
		chi2 = bm  * pRG->R[ks][j+sy][i+sx][ifr].chi +
                       bm1 * pRG->R[ks][j   ][i+sx][ifr].chi;
		dtaum = 0.5 * (chi0 + chi1) * dx * muinv[m][0]; 
		dtaup = 0.5 * (chi2 + chi1) * dx * muinv[m][0];
	      }
	      get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);
	      psi[j][i][ifr][l][m][0] = edtau;
	      psi[j][i][ifr][l][m][1] = a0;
	      psi[j][i][ifr][l][m][2] = a1;
	      psi[j][i][ifr][l][m][3] = a2;
	      lamstr[j][i][ifr] += pRG->wmu[m] * a1;

	      if (l == 0) {
		if (am <= 1.0) {
		  psiint[j][i][ifr][1] += am  * pRG->wmu[m] * a2; 
		  psiint[j][i][ifr][2] += am1 * pRG->wmu[m] * a2; 
		} else {
		  psiint[j][i][ifr][0] += bm1 * pRG->wmu[m] * a2; 
		  psiint[j][i][ifr][1] += bm  * pRG->wmu[m] * a2; 
		}
	      } else if (l == 1) {
		if (am <= 1.0) {
		  psiint[j][i][ifr][2] += am1 * pRG->wmu[m] * a2; 
		  psiint[j][i][ifr][3] += am  * pRG->wmu[m] * a0; 
		}
	      } else if (l == 2) {
		if (am <= 1.0) {
		  psiint[j][i][ifr][3] += am * pRG->wmu[m] * a2; 
		} else {
		  psiint[j][i][ifr][0] += bm1 * pRG->wmu[m] * a2; 
		  psiint[j][i][ifr][3] += bm  * pRG->wmu[m] * a2;
		}
	      }
	    }
	  }
	}
  }

  return;

  on_error:
  formal_solution_2d_destruct();
  ath_error("[formal_solution__2d_init]: Error allocating memory\n");
  return;

}

#endif /* GAUSSEID */
#endif /* RADIATION_TRANSFER */
