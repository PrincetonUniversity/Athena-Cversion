#include "copyright.h"

/*==============================================================================
 * FILE: rad2d.c
 *
 * PURPOSE:  Problem generator for a 2D test of radiative transfer routine
 *
 * Initial conditions available:
 *  iang = 1 - use Gauss-Legendre quadrature
 *  iang = 2 - use Eddington approximation
 *
 *============================================================================*/

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

static Real eps0;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * gauleg()           - gauss-legendre quadrature from NR
 *============================================================================*/
static void gauleg(Real x1, Real x2,  Real *x, Real *w, int n);
static Real const_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
static Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);

void problem(DomainS *pDomain)
{
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int nmu=pRG->nmu, ng=pRG->ng, nf=pRG->nf;
  int i, j, k, ifr, l;
  int iang;
  Real x2, ytop, ybtm;  
  Real den = 1.0;
  Real wedd[2] = {0.5, 0.5};
  Real muedd[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  Real stheta;
  Real *mul = NULL, *wl = NULL, *tau = NULL, *gl = NULL, *wgl = NULL;

/* Read problem parameters. */

  iang = par_geti_def("problem","iang", 1);
  eps0 = par_getd("problem","eps");

/* Setup density structure */ 
/* tau is used to initialize density grid */
  if ((tau = (Real *)calloc_1d_array(pG->Nx[1]+nghost+2,sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory");
  }

  ytop = pDomain->RootMaxX[1];
  ybtm = pDomain->RootMinX[1];

  for(j=pG->js; j<=pG->je+1; j++) {
    x2 = pG->MinX[1] + (Real)(j-js)*pG->dx2;
    tau[j] = pow(10.0,-3.0 + 10.0 * ((x2-ybtm)/(ytop-ybtm)));
  }

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      pG->U[ks][j][i].d  = (tau[j+1] - tau[j]) / pG->dx2;
    }
  }

/* Initialize arrays of angles and weights for angular quadrature */

  if (iang == 1) { 
/* angles and weight determined by gauss-legendre */
    pRG->r1ls=0; pRG->r1le=nmu-1; pRG->l1ls=0; pRG->l1le=nmu-1;
    pRG->r1ms=ng/2; pRG->r1me=ng-1; pRG->l1ms=0; pRG->l1me=ng/2-1;
    pRG->r2ls=nmu/2; pRG->r2le=nmu-1; pRG->l2ls=0; pRG->l2le=nmu/2-1;
    pRG->r2ms=0; pRG->r2me=ng-1; pRG->l2ms=0; pRG->l2me=ng-1;
    if ((mul = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((wl = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      free_1d_array(mul);
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((gl = (Real *)calloc_1d_array(ng,sizeof(Real))) == NULL) {
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((wgl = (Real *)calloc_1d_array(ng,sizeof(Real))) == NULL) {
      free_1d_array(mul);
      ath_error("[rad1d]: Error allocating memory");
    }
    gauleg(0.0, 1.0, mul, wl, nmu/2);
    gauleg(0.0, 2.0 * PI, gl, wgl, ng);

    for(i=0; i<nmu/2; i++) {
      pRG->mu[i] = -mul[nmu/2-i-1];
      stheta = sqrt(1.0 - pRG->mu[i] * pRG->mu[i]);
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = 0.25 * wl[nmu/2-i-1] * wgl[j] / PI; 
	pRG->gamma[i][j] = cos(gl[j]) * stheta;  
      } 
    }
    for(i=nmu/2; i<nmu; i++) {
      pRG->mu[i] = mul[i-nmu/2];   
      stheta = sqrt(1.0 - pRG->mu[i] * pRG->mu[i]);
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = 0.25 * wl[i-nmu/2] * wgl[j] / PI; 
	pRG->gamma[i][j] = cos(gl[j]) * stheta; 
      } 
    }

/* free up memory */
    free_1d_array(mul);
    free_1d_array(wl);
    free_1d_array(gl);
    free_1d_array(wgl);     
  } else if (iang == 2) {
/* Angles reproduce Eddington approximation */
    nmu = pRG->nmu = 2; ng = pRG->ng = 2;
    pRG->r1ls=0; pRG->r1le=1; pRG->l1ls=0; pRG->l1le=1;
    pRG->r1ms=1; pRG->r1me=1; pRG->l1ms=0; pRG->l1me=0;
    pRG->r2ls=1; pRG->r2le=1; pRG->l2ls=0; pRG->l2le=0;
    pRG->r2ms=0; pRG->r2me=1; pRG->l2ms=0; pRG->l2me=1;
    for(i=0; i<nmu; i++) {
      pRG->mu[i] = muedd[i];
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = wedd[i] * wedd[j];
	pRG->gamma[i][j] = muedd[j];
      }
    }
  }
/* ------- Initialize boundary emission ---------------------------------- */

  for(j=0; j<nmu; j++) 
    for(k=0; k<ng; k++) 
      for(ifr=0; ifr<nf; ifr++) {
	pRG->r1imu[pRG->ks][0][ifr][j][k] = 0.0;
	pRG->l1imu[pRG->ks][0][ifr][j][k] = 0.0;
      }

  for(j=pRG->js; j<=pRG->je+1; j++) {
    for(i=pRG->r1ls; i<=pRG->r1le; i++) 
      for(k=pRG->r1ms; k<=pRG->r1me; k++)
	for(ifr=0; ifr<nf; ifr++)
	  pRG->l1imu[pRG->ks][j][ifr][i][k] = 1.0;

    for(i=pRG->l1ls; i<=pRG->l1le; i++) 
      for(k=pRG->l1ms; k<=pRG->l1me; k++) 
	for(ifr=0; ifr<nf; ifr++)
	  pRG->r1imu[pRG->ks][j][ifr][i][k] = 1.0;
  }

  for(i=pRG->is-1; i<=pRG->ie+1; i++) {
/* incident radiation at lower boundary */
    for(j=pRG->r2ls; j<=pRG->r2le; j++) 
      for(k=pRG->r2ms; k<=pRG->r2me; k++)
	for(ifr=0; ifr<nf; ifr++)
	  for(l=0; l<2; l++)
/* lower boundary is tau=0, no irradiation */
	    pRG->l2imu[pRG->ks][i][ifr][j][k][l] = 0.0;

/* incident radiation at upper boundary */
    for(j=pRG->l2ls; j<=pRG->l2le; j++) 
      for(k=pRG->l2ms; k<=pRG->l2me; k++) 
	for(ifr=0; ifr<nf; ifr++)
	  for(l=0; l<2; l++)
/* upper boundary is large tau, eps=1 */
	    pRG->r2imu[pRG->ks][i][ifr][j][k][l] = 1.0;
  }

/* enroll radiation specification functions */
get_thermal_source = const_B;
get_thermal_fraction = const_eps;
get_total_opacity = const_opacity;

/* Free up memory */
  free_1d_array(tau);

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}


void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

static Real const_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
  return 1.0;
}

static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k)
{

  return eps0;
  
}

static Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k)
{

  return pG->U[k][j][i].d;
  
}

static void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
{

  Real eps = 3.0e-14;
  Real xm, xl, z, z1;
  Real p1, p2, p3, pp;
  int m, i, j;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=1; i<=m; i++) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2.0 * (Real)j - 1.0) * z * p2 - ((Real)j - 1.0) * p3) / (Real)j;
      }
      pp = (Real)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}

