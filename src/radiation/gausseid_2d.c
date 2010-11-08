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
 *   gausseid_pass_pointers_to_mg_2d()
 *============================================================================*/


#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define INTERP_2D

#ifdef RADIATION
#ifdef GAUSSEID

#ifdef RAD_MULTIG
static Real *******psi = NULL;
static Real *****psiint = NULL;
#else
static Real ******psi = NULL;
static Real ****psiint = NULL;
#endif
static Real *****imuo = NULL;
static Real *mu1 = NULL, *gamma1 = NULL, **am0 = NULL;
static Real sorw, dsrold, dSrmx;
static int isor;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   sweep_2d()     - computes a single sweep in one direction (right or left)
 *   update_sfunc() - updates source function after compute of mean intensity
 *   set_bvals_imu_y()      - set imu array at vertical boundary
 *   set_bvals_imu_y_j()    - set imu array at horizontal boundary
 *   update_bvals_imu_y()   - update outgoing radiation at vertical boundary
 *   update_bvals_imu_y_j() - update outgoing radiation at horizontal boundary
 *============================================================================*/
static void sweep_2d(RadGridS *pRG, int j, int sx, int sy);
static void update_sfunc(RadS *R, Real *deltas);
static void set_bvals_imu_y(RadGridS *pRG, int sy);
static void set_bvals_imu_y_j(RadGridS *pRG, int j, int sy);
static void update_bvals_imu_y(RadGridS *pRG, int sy);
static void update_bvals_imu_y_j(RadGridS *pRG, int j, int sx, int sy);

void formal_solution_2d(RadGridS *pRG, Real *dSrmax)
{
  int i, j, l, m;
  int nmu = pRG->nmu, nmu2 = pRG->nmu / 2;
  int ng = pRG->ng, ng2 = pRG->ng / 2;
  int is = pRG->is, ie = pRG->ie; 
  int js = pRG->js, je = pRG->je; 
  int ks = pRG->ks; 
  int ifr, nf = pRG->nf;

/* Initialize dSrmx */
  dSrmx = 0.0;

/* initialize mean intensities at all depths to zero */
  for(j=js-1; j<je+1; j++)
    for(i=is-1; i<ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++)
	pRG->R[ks][j][i][ifr].J = 0.0;

/* Compute formal solution for all upward going rays in 
 * each vertical gridzone */

/* set ix2 boundary condition */
  set_bvals_imu_y(pRG, 1);

  for(j=js; j<=je; j++) {
/* set ix1/ox1 boundary conditions*/
    set_bvals_imu_y_j(pRG, j, 1);

/* Sweep forward and upward */
    sweep_2d(pRG, j, 1, 1);

/* Update ox1 boundary intensities */
    update_bvals_imu_y_j(pRG, j, 1, 1);

/* Sweep backward and upward */
    sweep_2d(pRG, j, -1, 1);

/* Update ix1 boundary intensities */
    update_bvals_imu_y_j(pRG, j, -1, 1);
  }

/* Update ox2 boundary intensities */
  update_bvals_imu_y(pRG, 1);
  
/* Compute formal solution for all downward going rays in 
 * each vertical gridzone */

/* set ox2 boundary condition */
  set_bvals_imu_y(pRG, -1);

  for(j=je; j>=js; j--) {
/* set ix1/ox1 boundary conditions*/
    set_bvals_imu_y_j(pRG, j, -1);

/* Sweep forward and downward */
    sweep_2d(pRG, j, 1, -1);

/* Update ox1 boundary intensities */   
    update_bvals_imu_y_j(pRG, j, 1, -1);

/* Sweep backward and downward */
    sweep_2d(pRG, j, -1, -1);

/* Update ix1 boundary intensities */
    update_bvals_imu_y_j(pRG, j, -1, -1);
  }
/* Update ix2 boundary intensities */
  update_bvals_imu_y(pRG, -1);

/* Evaluate relative change if SOR */
  if ((dSrmx < 0.03) && (isor == 1) ) {
    sorw = 2.0 / (1.0 + sqrt(1.0 - dSrmx/dsrold));
    isor = 0;
  }
  dsrold = dSrmx;
/* Return maximum relative change to test convergence*/
  (*dSrmax) = dSrmx;

  return;
}


static void sweep_2d(RadGridS *pRG, int j, int sx, int sy)
{
  int it0, i, l, m;
  int nmu = pRG->nmu, nmu2 = pRG->nmu/2;
  int ng = pRG->ng, ng2 = pRG->ng/2;
  int is = pRG->is, ie = pRG->ie;
  int ks = pRG->ks;
  int ifr, nf = pRG->nf;
  Real imu, imu0;
  Real S0, S2;
  Real am, am1, bm, bm1;
  Real w0, w1, w2;
  Real mus, mue, gs, ge;
  Real maxint, minint;
  Real deltas;

  if (sy == 1) {
    mus = nmu2;
    mue = nmu-1;
  } else {
    sy = -1;
    mus = 0;
    mue = nmu2-1;
  }

  if (sx == 1) {
    gs = ng2;
    ge = ng-1;
  } else {
    sx = -1;
    gs = 0;
    ge = ng2-1;
  }

  for(it0=is; it0<=ie; it0++) {
     if (sx == 1) 
       i = it0;
     else
       i = ie + is - it0;

     for(ifr=0; ifr<nf; ifr++) {
       for(l=mus; l<=mue; l++)
	 for(m=gs; m<=ge; m++) {

/* --------- Interpolate intensity and source functions at endpoints --------- 
 * --------- of characteristics                                      --------- */

	   am = am0[l][m];
	   if (am <= 1.0) {
	     am1 = 1.0 - am;
/* Use linear interpolation for source functions */
	     S0 = am  * pRG->R[ks][j-sy][i-sx][ifr].S +
                  am1 * pRG->R[ks][j-sy][i   ][ifr].S;
	     S2 = am  * pRG->R[ks][j+sy][i+sx][ifr].S +
                  am1 * pRG->R[ks][j+sy][i   ][ifr].S;
#ifdef INTERP_2D
/* Use parabolic interpolation for intensity */
	     w0 = 0.5 * am * (1.0 + am);
	     w1 = am1 * (1.0 + am);
	     w2 = -0.5 * am * am1;
	     imu0 = w0 * imuo[i-sx][ifr][l][m][1] + w1 * imuo[i][ifr][l][m][0] +
	            w2 * imuo[i+sx][ifr][l][m][0];
	     maxint = MAX(imuo[i-sx][ifr][l][m][1],imuo[i][ifr][l][m][0]);
	     minint = MIN(imuo[i-sx][ifr][l][m][1],imuo[i][ifr][l][m][0]);
	     if(imu0 > maxint) imu0 = maxint;
	     if(imu0 < minint) imu0 = minint;
#else
/* Use linear interpolation for intensity */
	     imu0 = am  * imuo[i-sx][ifr][l][m][1] +
                    am1 * imuo[i   ][ifr][l][m][0];
#endif
	   } else {
	     bm = 1.0 / am;
	     bm1 = 1.0 - bm;	     
/* Use linear interpolation for source functions */
	     S0 = bm  * pRG->R[ks][j-sy][i-sx][ifr].S +
                  bm1 * pRG->R[ks][j   ][i-sx][ifr].S;
	     S2 = bm  * pRG->R[ks][j+sy][i+sx][ifr].S +
                  bm1 * pRG->R[ks][j   ][i+sx][ifr].S;
#ifdef INTERP_2D 
/* Use parabolic interpolation for intensity */
	     w0 = 0.5 * bm1 * (1.0 + bm1);  // Modify to compute only once
	     w1 = bm * (1.0 + bm1);
	     w2 = -0.5 * bm * bm1;

	     imu0 = w0 * imuo[i-sx][ifr][l][m][0] + w1 * imuo[i-sx][ifr][l][m][1] +
	            w2 * imuo[i-sx][ifr][l][m][2];
	     maxint = MAX(imuo[i-sx][ifr][l][m][0],imuo[i-sx][ifr][l][m][1]);
	     minint = MIN(imuo[i-sx][ifr][l][m][0],imuo[i-sx][ifr][l][m][1]);
	     if(imu0 > maxint) imu0 = maxint;
	     if(imu0 < minint) imu0 =minint;
#else
/* Use linear interpolation for intensity */
	     imu0 = bm  * imuo[i-sx][ifr][l][m][1] +
                    bm1 * imuo[i-sx][ifr][l][m][0];
#endif
	}
/* ---------  compute intensity at grid center and add to mean intensity ------- */
#ifdef RAD_MULTIG
	   imu = psi[img][j][i][ifr][l][m][0] * imu0 + 
	         psi[img][j][i][ifr][l][m][1] * S0 +
	         psi[img][j][i][ifr][l][m][2] * pRG->R[ks][j][i][ifr].S +
	         psi[img][j][i][ifr][l][m][3] * S2;	
#else
	   imu = psi[j][i][ifr][l][m][0] * imu0 + 
	         psi[j][i][ifr][l][m][1] * S0 +
	         psi[j][i][ifr][l][m][2] * pRG->R[ks][j][i][ifr].S +
	         psi[j][i][ifr][l][m][3] * S2;	
#endif	   
/* Add to mean intensity and save for next iteration */
	   pRG->R[ks][j][i][ifr].J += pRG->w[l][m] * imu;
/* Update intensity workspace */
#ifdef INTERP_2D
	   imuo[i][ifr][l][m][2] = imuo[i][ifr][l][m][1];
#endif
	   imuo[i][ifr][l][m][1] = imuo[i][ifr][l][m][0];
	   imuo[i][ifr][l][m][0] = imu;
	 }
       if ((sy == -1) && (sx == -1)) {
	 update_sfunc(&(pRG->R[ks][j][i][ifr]), &deltas);
#ifdef RAD_MULTIG
	 for(l=mus; l<=mue; l++)
	   for(m=0; m<=ng-1; m++) 
	     imuo[i][ifr][l][m][0] += pRG->w[l][m] * 
	       psi[img][j][i][ifr][l][m][2] * deltas;
	 
	 pRG->R[ks][j  ][i-1][ifr].J += deltas * psiint[img][j  ][i-1][ifr][0];
	 pRG->R[ks][j-1][i-1][ifr].J += deltas * psiint[img][j-1][i-1][ifr][1];
	 pRG->R[ks][j-1][i  ][ifr].J += deltas * psiint[img][j-1][i  ][ifr][2];
	 pRG->R[ks][j-1][i+1][ifr].J += deltas * psiint[img][j-1][i+1][ifr][3];
#else
	 for(l=mus; l<=mue; l++)
	   for(m=0; m<=ng-1; m++) 
	     imuo[i][ifr][l][m][0] += pRG->w[l][m] * 
	       psi[j][i][ifr][l][m][2] * deltas;
	 
	 pRG->R[ks][j  ][i-1][ifr].J += deltas * psiint[j  ][i-1][ifr][0];
	 pRG->R[ks][j-1][i-1][ifr].J += deltas * psiint[j-1][i-1][ifr][1];
	 pRG->R[ks][j-1][i  ][ifr].J += deltas * psiint[j-1][i  ][ifr][2];
	 pRG->R[ks][j-1][i+1][ifr].J += deltas * psiint[j-1][i+1][ifr][3];
#endif
       }
     }
  }
  return;
}

static void update_sfunc(RadS *R, Real *deltas)
{
  Real snew, deltasr;
  
  snew = (1.0 - R->eps) * R->J + R->eps;

  (*deltas) = (snew - R->S) / (1.0 - (1.0 - R->eps) * R->lamstr);
  deltasr = fabs((*deltas) / R->S);
  (*deltas) *= sorw;
  R->S += (*deltas);
  if (deltasr > dSrmx) dSrmx = deltasr; 
  return;
}


void formal_solution_2d_destruct(void)
{
  int i;

#ifdef RAD_MULTIG
  if (psi != NULL) {
    for(i=0; i<nmgrid; i++)
      if (psi[i] != NULL) free_6d_array(psi[i]);
    free(psi);
  }
  if (psiint != NULL) {
    for(i=0; i<nmgrid; i++)
      if (psiint[i] != NULL) free_4d_array(psiint[i]);
    free(psiint);
  }
#else
  if (psi != NULL) free_6d_array(psi);
  if (psiint != NULL) free_4d_array(psiint);
#endif
  if (imuo != NULL) free_5d_array(imuo);
  if (mu1 != NULL) free(mu1);
  if (gamma1 != NULL) free(gamma1);
  if (am0 != NULL) free_2d_array(am0);

  return;
}


void formal_solution_2d_init(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
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

/* Initialize variables for sor */
  isor = par_geti_def("problem","isor",0);
  sorw = 1.0;
  dsrold = 0.0;

  if ((imuo = (Real *****)calloc_5d_array(nx1+2,nf,nmu,ng,3,sizeof(Real))) == NULL)
    goto on_error;

  if ((mu1 = (Real *)calloc(nmu,sizeof(Real))) == NULL)
    goto on_error;

  if ((gamma1 = (Real *)calloc(ng,sizeof(Real))) == NULL)
    goto on_error;

  if ((am0 = (Real **)calloc_2d_array(nmu,ng,sizeof(Real))) == NULL)
    goto on_error;

  for(i=0; i<nmu; i++) 
    mu1[i] = fabs(1.0 / pRG->mu[i]);

  for(i=0; i<ng; i++) 
    gamma1[i] = fabs(1.0 / pRG->gamma[i]);

  for(i=0; i<nmu; i++) 
    for(j=0; j<ng; j++) 
      am0[i][j] = fabs(dy * mu1[i] / (dx * gamma1[j]));
    

#ifdef RAD_MULTIG
  if ((psi = (Real *******)calloc(nmgrid,sizeof(Real******))) == NULL) 
    goto on_error;
  if ((psiint = (Real *****)calloc(nmgrid,sizeof(Real****))) == NULL) 
    goto on_error;

  for(i=0; i<nmgrid; i++) {
    psi[i] = NULL;
    psiint[i] = NULL;
  }

  if ((psi[0] = (Real ******)calloc_6d_array(nx2+2,nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL)
    goto on_error;

  if ((psiint[0] = (Real ****)calloc_4d_array(nx2+2,nx1+2,nf,4,sizeof(Real))) == NULL) 
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
	    psi[0][j][i][ifr][l][m][0] = edtau;
	    psi[0][j][i][ifr][l][m][1] = a0;
	    psi[0][j][i][ifr][l][m][2] = a1;
	    psi[0][j][i][ifr][l][m][3] = a2;
	    pRG->R[ks][j][i][ifr].lamstr += pRG->w[l][m] * a1;

	    if (sy == 1) 
	      if (sx == 1) 
		if (am <= 1.0) {
		  psiint[0][j][i][ifr][1] += am  * pRG->w[l][m] * a2; 
		  psiint[0][j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		} else {
		  psiint[0][j][i][ifr][0] += bm1 * pRG->w[l][m] * a2; 
		  psiint[0][j][i][ifr][1] += bm  * pRG->w[l][m] * a2; 
		}
	      else
		if (am <= 1.0) {
		  psiint[0][j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		  psiint[0][j][i][ifr][3] += am  * pRG->w[l][m] * a2; 
		} else 
		  psiint[0][j][i][ifr][3] += bm  * pRG->w[l][m] * a2; 
	  }
	}
      }
#else
  if ((psi = (Real ******)calloc_6d_array(nx2+2,nx1+2,nf,nmu,ng,4,sizeof(Real))) == NULL) 
    goto on_error;

  if ((psiint = (Real ****)calloc_4d_array(nx2+2,nx1+2,nf,4,sizeof(Real))) == NULL) 
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
	    psi[j][i][ifr][l][m][0] = edtau;
	    psi[j][i][ifr][l][m][1] = a0;
	    psi[j][i][ifr][l][m][2] = a1;
	    psi[j][i][ifr][l][m][3] = a2;
	    pRG->R[ks][j][i][ifr].lamstr += pRG->w[l][m] * a1;
	    
	    if (sy == 1) 
	      if (sx == 1) 
		if (am <= 1.0) {
		  psiint[j][i][ifr][1] += am  * pRG->w[l][m] * a2; 
		  psiint[j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		} else {
		  psiint[j][i][ifr][0] += bm1 * pRG->w[l][m] * a2; 
		  psiint[j][i][ifr][1] += bm  * pRG->w[l][m] * a2; 
		}
	      else
		if (am <= 1.0) {
		  psiint[j][i][ifr][2] += am1 * pRG->w[l][m] * a2; 
		  psiint[j][i][ifr][3] += am  * pRG->w[l][m] * a2; 
		} else 
		  psiint[j][i][ifr][3] += bm  * pRG->w[l][m] * a2; 	 
	  }
	}
      }
#endif /* RAD_MULTIG */

  return;

  on_error:
  formal_solution_2d_destruct();
  ath_error("[formal_solution_2d_init]: Error allocating memory\n");
  return;

}

#ifdef RAD_MULTIG
void gs_pass_pointers_to_mg_2d(Real ********psi0, Real **mu10, Real **gamma10, 
			       Real ***am00, Real ******psiint0)
{
  *psi0 = psi;
  *mu10 = mu1;
  *gamma10 = gamma1;
  *am00 = am0;
  *psiint0 = psiint;

  return;
}
#endif /* RAD_MULTIG */

static void set_bvals_imu_y(RadGridS *pRG, int sy)
{
  int is = pRG->is, ie = pRG->ie;
  int ks = pRG->ks;
  int nmu = pRG->nmu, nmu2 = pRG->nmu / 2;
  int ng = pRG->ng;
  int ifr, nf = pRG->nf;
  int i, l, m;

  if(sy == 1) {
/* Account for ix2 boundary condition */
    for(i=is-1; i<=ie+1; i++)
      for(ifr=0; ifr<nf; ifr++)
	for(l=nmu2; l<nmu; l++) 
	  for(m=0; m<ng; m++) {
	    imuo[i][ifr][l][m][0] = pRG->l2imu[ks][i][ifr][l][m][0];
	    imuo[i][ifr][l][m][1] = pRG->l2imu[ks][i][ifr][l][m][1];
#ifdef INTERP_2D
	    imuo[i][ifr][l][m][2] = pRG->l2imu[ks][i][ifr][l][m][1];
#endif
	  }

  } else {

/* Account for ox2 boundary condition */
    for(i=is-1; i<=ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++)
	for(l=0; l<nmu2; l++) 
	  for(m=0; m<ng; m++) {
	    imuo[i][ifr][l][m][0] = pRG->r2imu[ks][i][ifr][l][m][0];
	    imuo[i][ifr][l][m][1] = pRG->r2imu[ks][i][ifr][l][m][1];
#ifdef INTERP_2D
	    imuo[i][ifr][l][m][2] = pRG->r2imu[ks][i][ifr][l][m][1];
#endif
	  }
  }


  return;
}

static void set_bvals_imu_y_j(RadGridS *pRG, int j, int sy)
{
  int is = pRG->is, ie = pRG->ie;
  int ks = pRG->ks;
  int ng = pRG->ng;
  int ifr, nf = pRG->nf;
  int mus, mue;
  int l, m;

  if (sy == 1) {
    mus = pRG->nmu / 2;
    mue = pRG->nmu - 1;
  } else {
    mus = 0;
    mue = pRG->nmu / 2 - 1;
  }

  for(ifr=0; ifr<nf; ifr++)
    for(l=mus; l<=mue; l++) 
      for(m=0; m<ng; m++) {
/* set ix1/ox1 boundary conditions*/
#ifdef INTERP_2D
	if ((j >= pRG->js+1) && (j <= pRG->je-1)) {
	  imuo[is-1][ifr][l][m][2] = imuo[is-1][ifr][l][m][1];
	  imuo[ie+1][ifr][l][m][2] = imuo[ie+1][ifr][l][m][1];
	}
#endif
	imuo[is-1][ifr][l][m][1] = imuo[is-1][ifr][l][m][0];
	imuo[ie+1][ifr][l][m][1] = imuo[ie+1][ifr][l][m][0];
	imuo[is-1][ifr][l][m][0] = pRG->l1imu[ks][j][ifr][l][m];
	imuo[ie+1][ifr][l][m][0] = pRG->r1imu[ks][j][ifr][l][m];
      }

  return;
}

static void update_bvals_imu_y_j(RadGridS *pRG, int j, int sx, int sy)
{
  int is = pRG->is, ie = pRG->ie;
  int ks = pRG->ks;
  int ng = pRG->ng, nmu = pRG->nmu;
  int ifr, nf = pRG->nf;
  int mus, mue;
  int gs, ge;
  int l, m;
  
  if (sy == 1) {
    mus = nmu / 2;
    mue = nmu - 1;
  } else {
    mus = 0;
    mue = nmu / 2 - 1;
  }

  if(sx == 1) {
/* Update intensity at the ox1 boundary */
    gs = ng / 2;
    ge = ng - 1;
    for(ifr=0; ifr<nf; ifr++)
      for(l=mus; l<=mue; l++) 
	for(m=gs; m<=ge; m++) {
	  pRG->r1imu[ks][j][ifr][l][m] = imuo[pRG->ie][ifr][l][m][0];
	}
  } else {
/* Update intensity at the ix1 boundary */
    gs = 0;
    ge = ng / 2 - 1;
    for(ifr=0; ifr<nf; ifr++)
      for(l=mus; l<=mue; l++) 
	for(m=gs; m<=ge; m++) {
	  pRG->l1imu[ks][j][ifr][l][m] = imuo[pRG->is][ifr][l][m][0];
	}
  }

  return;
}

static void update_bvals_imu_y(RadGridS *pRG, int sy)
{
  int is = pRG->is, ie = pRG->ie;
  int ks = pRG->ks;
  int nmu2 = pRG->nmu / 2, nmu = pRG->nmu;
  int ng = pRG->ng;
  int ifr, nf = pRG->nf;
  int i, l, m, n;


  if(sy == 1)
/* Update intensity at the ox2 boundary */
    for(i=is-1; i<=ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++)
	for(l=nmu2; l<nmu; l++) 
	  for(m=0; m<ng; m++) 
	    for(n=0; n<2; n++)
	      pRG->r2imu[ks][i][ifr][l][m][n] = imuo[i][ifr][l][m][n];
  else
/* Update intensity at the ix2 boundary */
    for(i=is-1; i<=ie+1; i++) 
      for(ifr=0; ifr<nf; ifr++)
	for(l=0; l<nmu2; l++) 
	  for(m=0; m<ng; m++)
 	    for(n=0; n<2; n++)
	      pRG->l2imu[ks][i][ifr][l][m][n] = imuo[i][ifr][l][m][n];

  return;
}

#endif /* GAUSSEID */
#endif /* RADIATION */
