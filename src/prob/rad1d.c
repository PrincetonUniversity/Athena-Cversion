#include "copyright.h"

/*==============================================================================
 * FILE: rad1d.c
 *
 * PURPOSE:  Problem generator for a 1D test of radiative transfer routine
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
  int nmu = pRG->nmu, nf=pRG->nf;
  int i, j, k, ifr;
  int iang;
  Real x1, xtop, xbtm;  
  Real den = 1.0;
  Real wedd[2] = {0.5, 0.5};
  Real muedd[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  Real *mul = NULL, *wl = NULL, *tau = NULL;

/* Read problem parameters. */

  iang = par_geti_def("problem","iang", 1);
  eps0 = par_getd("problem","eps");
/* Setup density structure */ 
/* tau is used to initialize density grid */
  if ((tau = (Real *)calloc_1d_array(pG->Nx[0]+nghost+2,sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory");
  }

  xtop = pDomain->RootMaxX[0];
  xbtm = pDomain->RootMinX[0];

  for(i=pG->is; i<=pG->ie+2; i++) {
    x1 = pG->MinX[0] + (Real)(i-is)*pG->dx1;
    tau[i] = pow(10.0,-3.0 + 10.0 * ((x1-xbtm)/(xtop-xbtm)));
  }

  //for(i=pG->is; i<=pG->ie+1; i++) 
  //  tau[i] = pow(10.0,-3.0 + pG->dx1 * (Real)(i-pG->is-0.5) * 10.0);

  for (i=is-1; i<=ie+1; i++) {
    pG->U[ks][js][i].d  = (tau[i+1] - tau[i]) / pG->dx1;
  }

/* Initialize arrays of angles and weights for angular quadrature */
  pRG->ng=1;
  pRG->r1ms=0; pRG->r1me=0; pRG->l1ms=0; pRG->l1me=0;
  if (iang == 1) { 
/* angles and weight determined by gauss-legendre */
    pRG->r1ls=nmu/2; pRG->r1le=nmu-1; pRG->l1ls=0; pRG->l1le=nmu/2-1;

    if ((mul = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      ath_error("[rad1d]: Error allocating memory");
    }

    if ((wl = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      free_1d_array(mul);
      ath_error("[rad1d]: Error allocating memory");
    }
    gauleg(0.0, 1.0, mul, wl, nmu/2);
    for(i=0; i<nmu/2; i++) {
      pRG->mu[i] = -mul[nmu/2-i-1];   
      pRG->w[i][0] = 0.5 * wl[nmu/2-i-1];  
    }
    for(i=nmu/2; i<nmu; i++) {
      pRG->mu[i] = mul[i-nmu/2];   
      pRG->w[i][0] = 0.5 * wl[i-nmu/2];  
    } 
/* free up memory */
    free_1d_array(mul);
    free_1d_array(wl);
     
  } else if (iang == 2) {
/* Angles reproduce Eddington approximation */
    nmu = pRG->nmu = 2;
    pRG->r1ls=1; pRG->r1le=1; pRG->l1ls=0; pRG->l1le=0;
    for(i=0; i<nmu; i++) {
      pRG->mu[i] = muedd[i];   
      pRG->w[i][0] = wedd[i];  
    }
  }

/* ------- Initialize boundary emission ---------------------------------- */

  /* incident radiation at lower (iz=0)  boundary */
  /* lower boundary is tau=0, no irradiation */
  for(j=nmu/2; j<nmu; j++) 
    for(ifr=0; ifr<nf; ifr++) 
    pRG->l1imu[pRG->ks][pRG->js][ifr][j][0] = 0.0; 

  /* incident radiation at upper (iz=nx) boundary */
  /* upper boundary is large tau, eps=1 */
  for(j=0; j<nmu/2; j++) 
    for(ifr=0; ifr<nf; ifr++) 
      pRG->r1imu[pRG->ks][pRG->js][ifr][j][0] = 1.0;

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

  z1 = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5)) + 2.0 * eps;
  for (i=1; i<=m; i++) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    while(fabs(z - z1) > eps) {
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
    }
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}

