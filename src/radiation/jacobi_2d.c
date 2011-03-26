 #include "../copyright.h"
/*==============================================================================
 * FILE: jacobi_2d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 2D grid using jacobi's method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution_2d.c()
 *   formal_solution_2d_destruct()
 *   formal_solution_2d_init()
 *   jacobi_pass_pointers_to_mg_2d()
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
#ifdef JACOBI


static Real ******psi = NULL;
static Real *****imu1 = NULL, *****imu2 = NULL;
static Real **muinv = NULL, *am0 = NULL, ***mu2 = NULL;
static Real ***Jold = NULL;
static int ntot, lte;
static int svwght, frstit;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_2d()     - computes a single sweep in one direction (right or left)
 *   update_sfunc() - updates source function after compute of mean intensity
 *   set_bvals_imu_y()      - set imu array at vertical boundary
 *   set_bvals_imu_y_j()    - set imu array at horizontal boundary
 *   update_bvals_imu_y()   - update outgoing radiation at vertical boundary
 *   update_bvals_imu_y_j() - update outgoing radiation at horizontal boundary
 *============================================================================*/

static void update_sfunc(RadS *R, Real *dSr);
static void update_cell(RadGridS *pRG, Real *****imuo, int ifr, int k, int j, int i, int l, int m);
static void sweep_2d_forward(RadGridS *pRG);
static void sweep_2d_backward(RadGridS *pRG);

void formal_solution_2d(RadGridS *pRG, Real *dSrmax)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie; 
  int js = pRG->js, je = pRG->je; 
  int ks = pRG->ks; 
  int ifr, nf = pRG->nf;  
  Real dSr, dJ, dJmax;

/* if LTE then store J values from previous iteration */
  if(lte != 0)
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) 
	for(ifr=0; ifr<nf; ifr++) 
	  Jold[j][i][ifr] = pRG->R[ks][j][i][ifr].J;

/* initialize mean intensities at all depths to zero */
  for(j=js-1; j<+je+1; j++)
    for(i=is-1; i<=ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++) {
	pRG->R[ks][j][i][ifr].J = 0.0;
	pRG->R[ks][j][i][ifr].K[0] = 0.0;
	pRG->R[ks][j][i][ifr].K[1] = 0.0;
	pRG->R[ks][j][i][ifr].K[2] = 0.0;
      }

/* Compute formal solution for all upward and rightward going rays in 
 * each gridzone */
  sweep_2d_forward(pRG);

/* Compute formal solution for all downward and leftward going rays in 
 * each vertical gridzone */
  sweep_2d_backward(pRG);

  if(lte == 0) {
/* Update source function */
    (*dSrmax) = 0.0;
    for(j=js; j<=je; j++) 
      for(i=is; i<=ie; i++) 
	for(ifr=0; ifr<nf; ifr++) {
	  update_sfunc(&(pRG->R[ks][j][i][ifr]),&dSr);
	  if( dSr > (*dSrmax)) (*dSrmax) = dSr;
	}
  } else {
/* Use delta J / J as convergence criterion */
    (*dSrmax) = 0.0;
    dJmax = 0.0;
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) 
	for(ifr=0; ifr<nf; ifr++) {
	  dJ = fabs(pRG->R[ks][j][i][ifr].J - Jold[j][i][ifr]);
	  if(dJ > dJmax) dJmax = dJ;
	  if (Jold[j][i][ifr] > 0.0)
	    dSr = dJ / Jold[j][i][ifr];
	  else
	    dSr = 0;
	  if( dSr > (*dSrmax)) (*dSrmax) = dSr;	 
	}
    if(((*dSrmax) == 0.0) && (dJmax > 0.0)) (*dSrmax) = 1.0;
  }

/* flag for storing lamstr on first iteration */
  if(svwght == 0) frstit = 0;

  return;
}

static void update_cell(RadGridS *pRG, Real *****imuo, int ifr, int k, int j, int i, int l, int m)
{

  int im, ip, jm, jp, imm, imp;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real am, am1, bm, bm1;
  Real w0, w1, w2;
  Real maxint, minint;
  Real dx = pRG->dx1, dy = pRG->dx2;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;

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
      dtaum = 0.5 * (chi0 + chi1) * dy * muinv[m][1]; 
      dtaup = 0.5 * (chi2 + chi1) * dy * muinv[m][1]; 
    } else {
      chi0 = bm  * pRG->R[k][jm][im][ifr].chi + 
	     bm1 * pRG->R[k][j ][im][ifr].chi;
      chi2 = bm  * pRG->R[k][jp][ip][ifr].chi +
	     bm1 * pRG->R[k][j ][ip][ifr].chi;
      dtaum = 0.5 * (chi0 + chi1) * dx * muinv[m][0]; 
      dtaup = 0.5 * (chi2 + chi1) * dx * muinv[m][0];
    }
    get_weights_parabolic(dtaum, dtaup, &edtau, &a0, &a1, &a2);
    imu = a0 * S0 + a1 * pRG->R[k][j][i][ifr].S + a2 * S2;
    if (imu < 0.0) imu=0.0;
    imu += edtau * imu0;
/* only compute lamstr on first iteration */
    if(frstit == 1) pRG->R[k][j][i][ifr].lamstr += pRG->wmu[m] * a1;
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
	    if(j == js) {
	      imu1[ifr][j-1][0][m][1] = pRG->l2imu[ifr][ks][i-1][0][m];
	      imu1[ifr][j-1][2][m][1] = pRG->l2imu[ifr][ks][i-1][2][m];
	    }
	    if(j == je) {
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
	for(m=0; m<nang; m++)
	  if(am0[m] > 1.0) {
/* modify imu1 to include top and boundar boundary intensity */
	    if(j == js) {
	      imu1[ifr][j-1][1][m][0] = pRG->l2imu[ifr][ks][i+1][1][m];
	      imu1[ifr][j-1][3][m][0] = pRG->l2imu[ifr][ks][i+1][3][m];
	    }
	    if(j == je) {
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


static void update_sfunc(RadS *R, Real *dSr)
{
  Real Snew, dS;
  
  Snew = (1.0 - R->eps) * R->J + R->eps * R->B;
  dS = (Snew - R->S) / (1.0 - (1.0 - R->eps) * R->lamstr);
  if (R->S > 0.0) (*dSr) = fabs(dS / R->S);
  R->S += dS;

  return;
}

void formal_solution_2d_destruct(void)
{
  int i;

  if (psi != NULL) free_6d_array(psi);
  if (imu1  != NULL) free_5d_array(imu1);
  if (imu2  != NULL) free_5d_array(imu2);
  if (muinv != NULL) free_2d_array(muinv);
  if (am0   != NULL) free_1d_array(am0);
  if (mu2   != NULL) free_3d_array(mu2);
  if (Jold  != NULL) free_3d_array(Jold);

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
  if (svwght == 0) frstit = 1;
  lte = par_geti("radiation","lte");

  ntot = je + ie - (js + is) + 1;

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

  if (lte != 0) 
    if ((Jold = (Real ***)calloc_3d_array(nx2+2,nx1+2,nf,sizeof(Real))) == NULL)
      goto on_error;
  
  if ((psi = (Real ******)calloc_6d_array(nx2+2,nx1+2,nf,4,nang,4,sizeof(Real))) == NULL) 
    goto on_error;

  if(svwght == 1) {
/* compute weights once and save for next iteration */
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
	      pRG->R[ks][j][i][ifr].lamstr += pRG->wmu[m] * a1;
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

#ifdef RAD_MULTIG
void jacobi_pass_pointers_to_mg_2d(Real ********psi0, Real **muinv0, Real **am00)
{
  *psi0 = psi;
  *muinv0 = muinv;
  *am00 = am0;
  return;
}
#endif /* RAD_MULTIG */

#endif /* JACOBI */
#endif /* RADIATION_TRANSFER */
