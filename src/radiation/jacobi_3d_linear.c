
#include "../copyright.h"
/*==============================================================================
 * FILE: jacobi_3d_linear.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 3D grid using jacobi's method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution_3d.c()
 *   formal_solution_3d_destruct()
 *   formal_solution_3d_init()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"


#ifdef RADIATION_TRANSFER
#ifdef JACOBI_LINEAR

static Real ****lamstr = NULL;
static Real ******imuo = NULL;
static int *face = NULL;
static Real **muinv = NULL, **coeff = NULL, ***mu2 = NULL;
static Real ****Jold = NULL;


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_3d()     - computes a single sweep in one direction (right or left)
 *   update_sfunc() - updates source function after compute of mean intensity
 *   set_bvals_imu_y()      - set imu array at vertical boundary
 *   set_bvals_imu_y_j()    - set imu array at horizontal boundary
 *   update_bvals_imu_y()   - update outgoing radiation at vertical boundary
 *   update_bvals_imu_y_j() - update outgoing radiation at horizontal boundary
 *============================================================================*/

static void update_sfunc(RadS *R, Real *dSr, Real lam);
static void sweep_3d_forward(RadGridS *pRG);
static void sweep_3d_backward(RadGridS *pRG);
static void update_cell(RadGridS *pRG, Real ******imuo, int k, int j, int i, int l);

void formal_solution_3d(RadGridS *pRG, Real *dSrmax)
{
  int i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie; 
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int ifr, nf = pRG->nf;
  int ismx, jsmx, ksmx;
  Real dSr, dJ, dJmax;

  /*for (m=0;m<pRG->nang;m++) {
    printf("0: %d %g %g %g \n",m,pRG->l1imu[0][ks-1][js-1][0][m],
	   pRG->l2imu[0][ks-1][is-1][0][m], pRG->l3imu[0][js-1][is-1][0][m]);
    printf("1: %d %g %g %g \n",m,pRG->r1imu[0][ks-1][js-1][1][m],
	   pRG->l2imu[0][ks-1][ie+1][1][m], pRG->l3imu[0][js-1][ie+1][1][m]);
    printf("2: %d %g %g %g \n",m,pRG->l1imu[0][ks-1][je+1][2][m],
	   pRG->r2imu[0][ks-1][is-1][2][m], pRG->l3imu[0][je+1][is-1][2][m]);
    printf("3: %d %g %g %g \n",m,pRG->r1imu[0][ks-1][je+1][3][m],
	   pRG->r2imu[0][ks-1][ie+1][3][m], pRG->l3imu[0][je+1][ie+1][3][m]);
    printf("4: %d %g %g %g \n",m,pRG->l1imu[0][ke+1][js-1][4][m],
	   pRG->l2imu[0][ke+1][is-1][4][m], pRG->r3imu[0][js-1][is-1][4][m]);
    printf("5: %d %g %g %g \n",m,pRG->r1imu[0][ke+1][js-1][5][m],
	   pRG->l2imu[0][ke+1][ie+1][5][m], pRG->r3imu[0][js-1][ie+1][5][m]);
    printf("6: %d %g %g %g \n",m,pRG->l1imu[0][ke+1][je+1][6][m],
	   pRG->r2imu[0][ke+1][is-1][6][m], pRG->r3imu[0][je+1][is-1][6][m]);
    printf("7: %d %g %g %g \n",m,pRG->r1imu[0][ke+1][je+1][7][m],
	   pRG->r2imu[0][ke+1][ie+1][7][m], pRG->r3imu[0][je+1][ie+1][7][m]);
    printf("---\n");
  }
  printf("\n");*/
/* if LTE then store J values from previous iteration */
  if(lte != 0) {
    for(k=ks; k<=ke; k++) 
      for(j=js; j<=je; j++)
	for(i=is; i<=ie; i++) 
	  for(ifr=0; ifr<nf; ifr++) 
	    Jold[k][j][i][ifr] = pRG->R[k][j][i][ifr].J;
  }

/* initialize mean intensities at all depths to zero */
  for(k=ks-1; k<=ke+1; k++) 
    for(j=js-1; j<=je+1; j++)
      for(i=is-1; i<=ie+1; i++) 
	for(ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][i][ifr].J = 0.0;
	  pRG->R[k][j][i][ifr].K[0] = 0.0;
	  pRG->R[k][j][i][ifr].K[1] = 0.0;
	  pRG->R[k][j][i][ifr].K[2] = 0.0;
	  pRG->R[k][j][i][ifr].K[3] = 0.0;
	  pRG->R[k][j][i][ifr].K[4] = 0.0;
	  pRG->R[k][j][i][ifr].K[5] = 0.0;
	  lamstr[k][j][i][ifr] = 0.0;
	}

/* Compute formal solution and for all rays in each gridzone and 
 * update boundary emission*/
  sweep_3d_forward(pRG);

  sweep_3d_backward(pRG);

  if(lte == 0) {
/* Update source function */
    (*dSrmax) = 0.0;
    for(k=ks; k<=ke; k++) 
      for(j=js; j<=je; j++) 
	for(i=is; i<=ie; i++) 
	  for(ifr=0; ifr<nf; ifr++) {
	    update_sfunc(&(pRG->R[k][j][i][ifr]),&dSr,lamstr[k][j][i][ifr]);
	    if( dSr > (*dSrmax)) {
	      (*dSrmax) = dSr; ismx=i; jsmx=j; ksmx=k;
	    }
	  }
  } else {
/* Use delta J / J as convergence criterion */
    (*dSrmax) = 0.0;
    dJmax = 0.0;
    for(k=ks; k<=ke; k++) 
      for(j=js; j<=je; j++)
	for(i=is; i<=ie; i++) 
	  for(ifr=0; ifr<nf; ifr++) {
	    dJ = fabs(pRG->R[k][j][i][ifr].J - Jold[k][j][i][ifr]);
	    if(dJ > dJmax) dJmax = dJ;
	    if (Jold[k][j][i][ifr] > 0.0)
	      dSr = dJ / Jold[k][j][i][ifr];
	    else
	      dSr = 0;
	    if( dSr > (*dSrmax)) {
	      (*dSrmax) = dSr; ismx=i; jsmx=j; ksmx=k;
	    }	 
	  }
    if(((*dSrmax) == 0.0) && (dJmax > 0.0)) (*dSrmax) = 1.0;
  }

  return;
}

static void sweep_3d_forward(RadGridS *pRG)
{
  int ifr, i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ix3 boundary intensities */
  for(ifr=0; ifr<nf; ifr++) {
    for(j=js-1; j<=je+1; j++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=0; l<=3; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][j][i][l][m][0] = pRG->l3imu[ifr][j][i][l][m];
	  }}}}}

  /* sweep forward in x3 */
  for(k=ks; k<=ke; k++) {

    /* Account for ix2 boundary intensities.  Note that this uses
     * l2imu to initialize imuo on edge */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=0; l<=1; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][js-1][i][l][m][1] = imuo[ifr][js-1][i][l][m][0];
	    imuo[ifr][js-1][i][l][m][0] = pRG->l2imu[ifr][k][i][l][m];
	  }}}}

    /* Sweep forward in x2 */
    for(j=js; j<=je; j++) {

      /* Account for ix1 boundary intensities */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++) {
	  /* ix1/ox1 boundary conditions*/
	  imuo[ifr][j][is-1][0][m][1] = imuo[ifr][j][is-1][0][m][0];
	  imuo[ifr][j][ie+1][1][m][1] = imuo[ifr][j][ie+1][1][m][0];
	  imuo[ifr][j][is-1][0][m][0] = pRG->l1imu[ifr][k][j][0][m];
	  imuo[ifr][j][ie+1][1][m][0] = pRG->r1imu[ifr][k][j][1][m];
	}}

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++) 
	update_cell(pRG,imuo,k,j,i,0);

      /* Update intensity at the ox1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->r1imu[ifr][k][j][0][m] = imuo[ifr][j][ie][0][m][0];
	}}

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--) 
	update_cell(pRG,imuo,k,j,i,1);

      /* Update intensity at the ix1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->l1imu[ifr][k][j][1][m] = imuo[ifr][j][is][1][m][0];
	}}

    }

    /* Update intensity at the ox2 boundary */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is; i<=ie; i++) { 
	for(l=0; l<=1; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->r2imu[ifr][k][i][l][m] = imuo[ifr][je][i][l][m][0];
	  }}}}

/* ----------------  Start of reverse sweep ---------------------- */

    /* Account for ox2 boundary intensities */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=2; l<=3; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][je+1][i][l][m][1] = imuo[ifr][je+1][i][l][m][0];
	    imuo[ifr][je+1][i][l][m][0] = pRG->r2imu[ifr][k][i][l][m];
	  }}}}

    /* sweep backward in x2 */
    for(j=je; j>=js; j--) {

      /* Account for ix1 boundary intensities */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++) {
	  /* ix1/ox1 boundary conditions*/
	  imuo[ifr][j][is-1][2][m][1] = imuo[ifr][j][is-1][2][m][0];
	  imuo[ifr][j][ie+1][3][m][1] = imuo[ifr][j][ie+1][3][m][0];
	  imuo[ifr][j][is-1][2][m][0] = pRG->l1imu[ifr][k][j][2][m];
	  imuo[ifr][j][ie+1][3][m][0] = pRG->r1imu[ifr][k][j][3][m];
	}}

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++) 
	update_cell(pRG,imuo,k,j,i,2);

      /* Update intensity at the ox1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->r1imu[ifr][k][j][2][m] = imuo[ifr][j][ie][2][m][0];
	}}

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--) 
	update_cell(pRG,imuo,k,j,i,3);

      /* Update intensity at the ix1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->l1imu[ifr][k][j][3][m] = imuo[ifr][j][is][3][m][0];
	}}
    }

    /* Update intensity at the ix2 boundary */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is; i<=ie; i++) { 
	for(l=2; l<=3; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->l2imu[ifr][k][i][l][m] = imuo[ifr][js][i][l][m][0];
	  }}}}
  }

   /* Update intensity at the ox3 boundary */
  for(ifr=0; ifr<nf; ifr++) {
    for(j=js; j<=je; j++) {
      for(i=is; i<=ie; i++) { 
	for(l=0; l<=3; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->r3imu[ifr][j][i][l][m] = imuo[ifr][j][i][l][m][0];
	  }}}}}

  return;
}

static void sweep_3d_backward(RadGridS *pRG)
{
  int ifr, i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ix3 boundary intensities */
  for(ifr=0; ifr<nf; ifr++) {
    for(j=js-1; j<=je+1; j++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=4; l<=7; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][j][i][l][m][0] = pRG->r3imu[ifr][j][i][l][m];
	  }}}}}

  /* sweep forward in x3 */
  for(k=ke; k>=ks; k--) {

    /* Account for ix2 boundary intensities.  Note that this uses
     * l2imu to initialize imuo on edge */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=4; l<=5; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][js-1][i][l][m][1] = imuo[ifr][js-1][i][l][m][0];
	    imuo[ifr][js-1][i][l][m][0] = pRG->l2imu[ifr][k][i][l][m];
	  }}}}

    /* Sweep forward in x2 */
    for(j=js; j<=je; j++) {

      /* Account for ix1 boundary intensities */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++) {
	  /* ix1/ox1 boundary conditions*/
	  imuo[ifr][j][is-1][4][m][1] = imuo[ifr][j][is-1][4][m][0];
	  imuo[ifr][j][ie+1][5][m][1] = imuo[ifr][j][ie+1][5][m][0];
	  imuo[ifr][j][is-1][4][m][0] = pRG->l1imu[ifr][k][j][4][m];
	  imuo[ifr][j][ie+1][5][m][0] = pRG->r1imu[ifr][k][j][5][m];
	}}

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++) 
	update_cell(pRG,imuo,k,j,i,4);

      /* Update intensity at the ox1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->r1imu[ifr][k][j][4][m] = imuo[ifr][j][ie][4][m][0];
	}}

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--) 
	update_cell(pRG,imuo,k,j,i,5);

      /* Update intensity at the ix1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->l1imu[ifr][k][j][5][m] = imuo[ifr][j][is][5][m][0];
	}}
    }

    /* Update intensity at the ox2 boundary */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is; i<=ie; i++) { 
	for(l=4; l<=5; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->r2imu[ifr][k][i][l][m] = imuo[ifr][je][i][l][m][0];
	  }}}}

/* ----------------  Start of reverse sweep ---------------------- */

    /* Account for ox2 boundary intensities */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is-1; i<=ie+1; i++) {
	for(l=6; l<=7; l++)  {
	  for(m=0; m<nang; m++) {
	    imuo[ifr][je+1][i][l][m][1] = imuo[ifr][je+1][i][l][m][0];
	    imuo[ifr][je+1][i][l][m][0] = pRG->r2imu[ifr][k][i][l][m];
	  }}}}

    /* sweep backward in x2 */
    for(j=je; j>=js; j--) {

      /* Account for ix1 boundary intensities */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++) {
	  /* ix1/ox1 boundary conditions*/
	  imuo[ifr][j][is-1][6][m][1] = imuo[ifr][j][is-1][6][m][0];
	  imuo[ifr][j][ie+1][7][m][1] = imuo[ifr][j][ie+1][7][m][0];
	  imuo[ifr][j][is-1][6][m][0] = pRG->l1imu[ifr][k][j][6][m];
	  imuo[ifr][j][ie+1][7][m][0] = pRG->r1imu[ifr][k][j][7][m];
	}}

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++) 
	update_cell(pRG,imuo,k,j,i,6);

      /* Update intensity at the ox1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->r1imu[ifr][k][j][6][m] = imuo[ifr][j][ie][6][m][0];
	}}

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--) 
	update_cell(pRG,imuo,k,j,i,7);

      /* Update intensity at the ix1 boundary */
      for(ifr=0; ifr<nf; ifr++) {
	for(m=0; m<nang; m++)  {
	  pRG->l1imu[ifr][k][j][7][m] = imuo[ifr][j][is][7][m][0];
	}}
    }

    /* Update intensity at the ix2 boundary */
    for(ifr=0; ifr<nf; ifr++) {
      for(i=is; i<=ie; i++) { 
	for(l=6; l<=7; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->l2imu[ifr][k][i][l][m] = imuo[ifr][js][i][l][m][0];
	  }}}}
  }

   /* Update intensity at the ox3 boundary */
  for(ifr=0; ifr<nf; ifr++) {
    for(j=js; j<=je; j++) {
      for(i=is; i<=ie; i++) { 
	for(l=4; l<=7; l++) { 
	  for(m=0; m<nang; m++) { 
	    pRG->l3imu[ifr][j][i][l][m] = imuo[ifr][j][i][l][m][0];
	  }}}}}

  return;
}


static void update_cell(RadGridS *pRG, Real ******imuo, int k, int j, int i, int l)

{

  int im, ip, jm, jp, km, kp;
  int ifr, m, nf = pRG->nf, nang = pRG->nang;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real w0, w1, w2;
  Real dx = pRG->dx1, dy = pRG->dx2, dz = pRG->dx3;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;

/* initialize stencil base on quadrant*/  
  if(l == 0) {
    kp = k + 1;  km = k - 1;
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 1) {
    kp = k + 1;  km = k - 1;
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 2) {
    kp = k + 1;  km = k - 1;
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 3) {
    kp = k + 1;  km = k - 1;
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 4) {
    kp = k - 1;  km = k + 1;
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 5) {
    kp = k - 1;  km = k + 1;
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 6) {
    kp = k - 1;  km = k + 1;
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else {
    kp = k - 1;  km = k + 1;
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  }  


  for(ifr=0; ifr<nf; ifr++) {
    for(m=0; m<nang; m++) {
      chi1 = pRG->R[k][j][i][ifr].chi;
/* --------- Interpolate intensity and source functions at endpoints --------- 
 * --------- of characteristics                                      --------- */
      if (face[m] == 0) {
	/* interpolation in x2-x3 plane */
	S0 = coeff[m][0] * pRG->R[k ][j ][im][ifr].S +
	     coeff[m][1] * pRG->R[km][j ][im][ifr].S +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].S +
	     coeff[m][3] * pRG->R[k ][jm][im][ifr].S;
	S2 = coeff[m][0] * pRG->R[k ][j ][ip][ifr].S +
	     coeff[m][1] * pRG->R[kp][j ][ip][ifr].S +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].S +
	     coeff[m][3] * pRG->R[k ][jp][ip][ifr].S;
      imu0 = coeff[m][0] * imuo[ifr][j ][im][l][m][0] +
	     coeff[m][1] * imuo[ifr][j ][im][l][m][1] +
	     coeff[m][2] * imuo[ifr][jm][im][l][m][1] +
	     coeff[m][3] * imuo[ifr][jm][im][l][m][0];
      chi0 = coeff[m][0] * pRG->R[k ][j ][im][ifr].chi +
	     coeff[m][1] * pRG->R[km][j ][im][ifr].chi +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].chi +
	     coeff[m][3] * pRG->R[k ][jm][im][ifr].chi;
      chi2 = coeff[m][0] * pRG->R[k ][j ][ip][ifr].chi +
	     coeff[m][1] * pRG->R[kp][j ][ip][ifr].chi +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].chi +
	     coeff[m][3] * pRG->R[k ][jp][ip][ifr].chi;
        interp_quad_chi(chi0,chi1,chi2,&dtaum);
	interp_quad_chi(chi2,chi1,chi0,&dtaup);
	dtaum *= dx * muinv[m][0]; 
	dtaup *= dx * muinv[m][0]; 
      } else if(face[m] == 1) {
	/* interpolation in x1-x3 plane */
	S0 = coeff[m][0] * pRG->R[k ][jm][i ][ifr].S +
	     coeff[m][1] * pRG->R[km][jm][i ][ifr].S +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].S +
	     coeff[m][3] * pRG->R[k ][jm][im][ifr].S;
	S2 = coeff[m][0] * pRG->R[k ][jp][i ][ifr].S +
	     coeff[m][1] * pRG->R[kp][jp][i ][ifr].S +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].S +
	     coeff[m][3] * pRG->R[k ][jp][ip][ifr].S;
      imu0 = coeff[m][0] * imuo[ifr][jm][i ][l][m][0] +
	     coeff[m][1] * imuo[ifr][jm][i ][l][m][1] +
	     coeff[m][2] * imuo[ifr][jm][im][l][m][1] +
	     coeff[m][3] * imuo[ifr][jm][im][l][m][0];
      chi0 = coeff[m][0] * pRG->R[k ][jm][i ][ifr].chi +
	     coeff[m][1] * pRG->R[km][jm][i ][ifr].chi +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].chi +
	     coeff[m][3] * pRG->R[k ][jm][im][ifr].chi;
      chi2 = coeff[m][0] * pRG->R[k ][jp][i ][ifr].chi +
	     coeff[m][1] * pRG->R[kp][jp][i ][ifr].chi +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].chi +
	     coeff[m][3] * pRG->R[k ][jp][ip][ifr].chi;
        interp_quad_chi(chi0,chi1,chi2,&dtaum);
	interp_quad_chi(chi2,chi1,chi0,&dtaup);
	dtaum *= dy * muinv[m][1]; 
	dtaup *= dy * muinv[m][1]; 
      } else  {
	/* interpolation in x1-x2 plane */
	S0 = coeff[m][0] * pRG->R[km][j ][i ][ifr].S +
	     coeff[m][1] * pRG->R[km][jm][i ][ifr].S +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].S +
	     coeff[m][3] * pRG->R[km][j ][im][ifr].S;
	S2 = coeff[m][0] * pRG->R[kp][j ][i ][ifr].S +
	     coeff[m][1] * pRG->R[kp][jp][i ][ifr].S +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].S +
	     coeff[m][3] * pRG->R[kp][j ][ip][ifr].S;	
      imu0 = coeff[m][0] * imuo[ifr][j ][i ][l][m][0] +
	     coeff[m][1] * imuo[ifr][jm][i ][l][m][0] +
	     coeff[m][2] * imuo[ifr][jm][im][l][m][0] +
	     coeff[m][3] * imuo[ifr][j ][im][l][m][0];
      /*if((l ==6) || (l ==7)) {
	if((m == 0) && (j == pRG->je) && (i == pRG->ie) && (k == pRG->ks)) {
	  printf("%d %d %d %d %d %d\n",km,k,jm,j,im,i);
	  printf("%g %g %g %g\n",imuo[ifr][j ][i ][l][m][0],imuo[ifr][jm][i ][l][m][0],
		 imuo[ifr][jm][im][l][m][0],imuo[ifr][j ][im][l][m][0]);
	  //	  printf("%g %g\n",pRG->r1imu[ifr][k][jm][l][m],pRG->r1imu[ifr][k][j][l][m]);
	  }}*/
      chi0 = coeff[m][0] * pRG->R[km][j ][i ][ifr].chi +
	     coeff[m][1] * pRG->R[km][jm][i ][ifr].chi +
	     coeff[m][2] * pRG->R[km][jm][im][ifr].chi +
	     coeff[m][3] * pRG->R[km][j ][im][ifr].chi;
      chi2 = coeff[m][0] * pRG->R[kp][j ][i ][ifr].chi +
	     coeff[m][1] * pRG->R[kp][jp][i ][ifr].chi +
	     coeff[m][2] * pRG->R[kp][jp][ip][ifr].chi +
	     coeff[m][3] * pRG->R[kp][j ][ip][ifr].chi;
        interp_quad_chi(chi0,chi1,chi2,&dtaum);
	interp_quad_chi(chi2,chi1,chi0,&dtaup);
	dtaum *= dz * muinv[m][2]; 
	dtaup *= dz * muinv[m][2]; 
      }
/* ---------  compute intensity at grid center and add to mean intensity ------- */
      interp_quad_source(dtaum, dtaup, &edtau, &a0, &a1, &a2,
			 S0, pRG->R[k][j][i][ifr].S, S2);
      imu = a0 * S0 + a1 * pRG->R[k][j][i][ifr].S + a2 * S2 + edtau * imu0;
      lamstr[k][j][i][ifr] += pRG->wmu[m] * a1;
    
/* Add to radiation moments and save for next iteration */
      wimu = pRG->wmu[m] * imu;
      pRG->R[k][j][i][ifr].J += wimu;
      //printf("%d %d %d %d %d %g %g %g %g\n",k,j,i,l,m,imu0,imu,wimu,pRG->R[k][j][i][ifr].J);
      /*if((l ==6) || (l ==7)) {
	if((j == pRG->je) && (i == pRG->ie) && (k == pRG->ks)) {
	  printf("%d %d %g %g %g %g\n",l,m,imu0,imu,wimu,pRG->R[k][j][i][ifr].J);
	  printf("%d %d %g %g\n",l,m,pRG->r1imu[ifr][k][j][6][m],pRG->r1imu[ifr][k][j][7][m]);
	
	  printf("%d %d %g %g\n",l,m,pRG->r1imu[ifr][k][j][6][m],pRG->l1imu[ifr][k][j][6][m]);
	  printf("%d %d %g %g\n",l,m,pRG->r1imu[ifr][k][j][7][m],pRG->l1imu[ifr][k][j][7][m]);
	  }}*/
      pRG->R[k][j][i][ifr].K[0] += mu2[l][m][0] * wimu;
      pRG->R[k][j][i][ifr].K[1] += mu2[l][m][1] * wimu;
      pRG->R[k][j][i][ifr].K[2] += mu2[l][m][2] * wimu;
      pRG->R[k][j][i][ifr].K[3] += mu2[l][m][3] * wimu;
      pRG->R[k][j][i][ifr].K[4] += mu2[l][m][4] * wimu;
      pRG->R[k][j][i][ifr].K[5] += mu2[l][m][5] * wimu;
/* Update intensity workspace */
      imuo[ifr][j][i][l][m][1] = imuo[ifr][j][i][l][m][0];
      imuo[ifr][j][i][l][m][0] = imu;
    }}
  
  return;
}


static void update_sfunc(RadS *R, Real *dSr, Real lam)
{
  Real Snew, dS;
  
  Snew = (1.0 - R->eps) * R->J + R->eps * R->B;
  dS = (Snew - R->S) / (1.0 - (1.0 - R->eps) * lam);
  if (R->S > 0.0) (*dSr) = fabs(dS / R->S);
  R->S += dS;

  return;
}

void formal_solution_3d_destruct(void)
{
  int i;

  if (lamstr != NULL) free_4d_array(lamstr);
  if (imuo   != NULL) free_6d_array(imuo);
  if (muinv  != NULL) free_2d_array(muinv);
  if (face   != NULL) free_1d_array(face);
  if (coeff  != NULL) free_2d_array(coeff);
  if (mu2    != NULL) free_3d_array(mu2);
  if (Jold   != NULL) free_4d_array(Jold);

  return;
}

void formal_solution_3d_init(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
  int nf = pRG->nf, nang = pRG->nang;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks; 
  Real dx = pRG->dx1, dy = pRG->dx2, dz = pRG->dx3;;
  int ifr, i, j, l, m;
  int sy, sx;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;
  Real am, bm;
  Real lx, ly, lz, lmin;

  if ((lamstr = (Real ****)calloc_4d_array(nx3+2,nx2+2,nx1+2,nf,sizeof(Real))) == NULL) 
   goto on_error;

  if ((imuo = (Real ******)calloc_6d_array(nf,nx2+2,nx1+2,8,nang,2,sizeof(Real))) == NULL)
    goto on_error;

  if ((muinv = (Real **)calloc_2d_array(nang,3,sizeof(Real))) == NULL)
    goto on_error;

  if ((face = (int *)calloc_1d_array(nang,sizeof(int))) == NULL)
    goto on_error;

  if ((coeff = (Real **)calloc_2d_array(nang,4,sizeof(Real))) == NULL)
    goto on_error;

  if ((mu2 = (Real ***)calloc_3d_array(8,nang,6,sizeof(Real))) == NULL)
    goto on_error;

  for(i=0; i<nang; i++)  
    for(j=0; j<3; j++) 
      muinv[i][j] = fabs(1.0 / pRG->mu[0][i][j]);

  for(i=0; i<nang; i++) {
    lx = muinv[i][0] * dx;
    ly = muinv[i][1] * dy;
    lz = muinv[i][2] * dz;
    lmin = MIN(MIN(lx,ly),lz);
    if (lz == lmin) {
      face[i]=2;
      am = lmin/lx; bm = lmin/ly;
    } else if (ly == lmin) {
      face[i]=1;
      am = lmin/lx; bm = lmin/lz;
    } else {
      face[i]=0;
      am = lmin/ly; bm = lmin/lz;
    }      
    coeff[i][0] = (1.0 - am)*(1.0 - bm);
    coeff[i][1] = (1.0 - am)*       bm;
    coeff[i][2] =        am *       bm;
    coeff[i][3] =        am *(1.0 - bm);
  }

  for(i=0; i<8; i++) 
    for(j=0; j<nang; j++)  {
      mu2[i][j][0] = pRG->mu[i][j][0] * pRG->mu[i][j][0];
      mu2[i][j][1] = pRG->mu[i][j][0] * pRG->mu[i][j][1];
      mu2[i][j][2] = pRG->mu[i][j][1] * pRG->mu[i][j][1];
      mu2[i][j][3] = pRG->mu[i][j][0] * pRG->mu[i][j][2];
      mu2[i][j][4] = pRG->mu[i][j][1] * pRG->mu[i][j][2];
      mu2[i][j][5] = pRG->mu[i][j][2] * pRG->mu[i][j][2];
    }

 
  if (lte != 0) 
    if ((Jold = (Real ****)calloc_4d_array(nx3+2,nx2+2,nx1+2,nf,sizeof(Real))) == NULL)
      goto on_error;

  return;

  on_error:
  formal_solution_3d_destruct();
  ath_error("[formal_solution__3d_init]: Error allocating memory\n");
  return;

}

#endif /* JACOBI_LINEAR */
#endif /* RADIATION_TRANSFER */
