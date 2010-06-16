#include "../copyright.h"
/*==============================================================================
 * FILE: rad_trans.c
 *
 * PURPOSE: Called from main loop before integrator.  Identifies and calls
 *          the appropriate formal solution algorithm.  If non-LTE will
 *          iteratively call formal_solution until convergence criterion
 *          is met.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   rad_trans()  - interate formal solution until convergence
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION

Real ***sol;

void get_solution_1d(RadGridS *pRG);
void get_solution_2d(RadGridS *pRG);
void max_dev(RadGridS *pRG, Real *dsm, int *ism);
void output_diag(Real *dsarr, int *isarr, int ns);
char *construct_filename(char *basename,char *key,int dump,char *ext);
void output_mean_intensity_2d(RadGridS *pRG, int itr);

void rad_trans(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  int i, niter, ndim;
  Real *dsmax, dsm;
  int  *isarr, ism;

  if ((dsmax = (Real *)calloc(10000,sizeof(Real))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }
  if ((isarr = (int *)calloc(10000,sizeof(int))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }

/* number of dimensions in Grid. */
  ndim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) ndim++;

/* maximum number of iterations */
  niter = par_geti("problem","niter");

/* compute formal solution */
  if (ndim == 1) {
    get_solution_1d(pRG);
    bvals_rad_1d_init(pRG);
    formal_solution_1d_init(pRG);

    for(i=0; i<niter; i++) {
#ifdef RAD_MULTIG
      formal_solution_mg_1d(pRG);
#else
      formal_solution_1d(pRG);
#endif
      max_dev(pRG,&dsm,&ism);
      dsmax[i] = dsm;
      isarr[i] = ism;
    }
    formal_solution_1d_destruct();
  } else if (ndim == 2) {
    get_solution_2d(pRG);
    bvals_rad_2d_init(pRG);
    formal_solution_2d_init(pRG);

    for(i=0; i<niter; i++) {
#ifdef RAD_MULTIG
      formal_solution_mg_2d(pRG);
#else
      formal_solution_2d(pRG);
#endif
      bvals_rad_2d(pRG);
      max_dev(pRG,&dsm,&ism);
      dsmax[i] = dsm;
      isarr[i] = ism;
    }
    output_mean_intensity_2d(pRG,0);
    formal_solution_2d_destruct();
  }

  output_diag(dsmax,isarr,niter);
  if (sol != NULL) free_3d_array(sol);

  return;
}


/* for testing purposes */
void get_solution_1d(RadGridS *pRG)
{
  int i, j, k;
  Real dtau, tau, jmean, eps;
  
  eps = par_getd("problem","eps");

  if ((sol = (Real ***)calloc_3d_array(pRG->Nx[2],pRG->Nx[1],pRG->Nx[0],
				       sizeof(Real))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }

  tau=0.0;
  for(i=pRG->is; i<=pRG->ie; i++) {
    dtau= 0.5 * pRG->dx1 * (pRG->R[0][0][i-1][0].chi + 
			   pRG->R[0][0][i  ][0].chi);
    tau += dtau;
    //    tau = pow(10.0,-3.0 + pRG->dx3 * (Real)(i-pRG->ks-0.5) * 10.0);
    jmean = 1.0 - sqrt(3.0) * exp(-sqrt(3.0 * eps) * tau) /
      (sqrt(3.0) + sqrt(3.0 * eps));
    for(k=pRG->ks; k<=pRG->ke; k++)
      for(j=pRG->js; j<=pRG->je; j++) 
	sol[k-pRG->ks][j-pRG->js][i-pRG->is] = eps + (1.0-eps) * jmean;
    //sol[pRG->nx3-2-i][j-1] = eps + (1.0-eps) * jmean;
  }

  return;

}

void get_solution_2d(RadGridS *pRG)
{
  int i, j, k;
  Real dtau, tau, jmean, eps;
  
  eps = par_getd("problem","eps");

  if ((sol = (Real ***)calloc_3d_array(pRG->Nx[2],pRG->Nx[1],pRG->Nx[0],
				       sizeof(Real))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }

  tau=0.0;
  for(j=pRG->js; j<=pRG->je; j++) {
     dtau= 0.5 * pRG->dx2 * (pRG->R[0][j-1][1][0].chi + 
			    pRG->R[0][j  ][1][0].chi);
    tau += dtau;
    //    tau = pow(10.0,-3.0 + pRG->dx3 * (Real)(i-pRG->ks-0.5) * 10.0);
    jmean = 1.0 - sqrt(3.0) * exp(-sqrt(3.0 * eps) * tau) /
      (sqrt(3.0) + sqrt(3.0 * eps));
     for(k=pRG->ks; k<=pRG->ke; k++)
      for(i=pRG->is; i<=pRG->ie; i++) 
	sol[k-pRG->ks][j-pRG->js][i-pRG->is] = eps + (1.0-eps) * jmean;
    //sol[pRG->nx3-2-i][j-1] = eps + (1.0-eps) * jmean;
  }

  return;

}

void max_dev(RadGridS *pRG, Real *dsm, int *ism)
{
  int i,j,k;
  int ixmax, iymax, izmax;
  Real ds, dst;

  ds = 0;
  for(k=pRG->ks; k<=pRG->ke; k++) 
    for(j=pRG->js; j<=pRG->je; j++)
      for(i=pRG->is; i<=pRG->ie; i++) {
	dst = fabs(pRG->R[k][j][i][0].S - sol[k-pRG->ks][j-pRG->js][i-pRG->is]) / 
	  sol[k-pRG->ks][j-pRG->js][i-pRG->is];
	//printf("%d %d %g %g\n",i,j,pRG->R[k][j][i][0].S,sol[k-pRG->ks][j-pRG->js][i-pRG->is]);
	if (dst > ds) {ds = dst; ixmax = i; iymax = j; izmax = k;}
      }

  (*dsm) = ds;
  (*ism) = iymax;
}

void output_diag(Real *dsarr, int *isarr, int ns)
{
  FILE *fp;
  int it1;

  if ((fp = fopen("diag.out", "w")) == NULL) {
    ath_error("## Error opening diag.out\n");
  }
  
  for(it1=0; it1<ns; it1++) 
    fprintf(fp,"%d %g %d\n",it1+1,dsarr[it1],isarr[it1]);

  fclose(fp);
  return;
}

void output_mean_intensity_2d(RadGridS *pRG, int itr)
{
  FILE *fp;
  char *fname;
  int i, j;
  Real tau, dtau;

  fname = construct_filename("imean", NULL, itr, "out");

  if ((fp = fopen(fname, "w")) == NULL) {
    ath_error("## Error opening imean.out\n");
  }

  tau = 0.0;
  for(i=pRG->js; i<=pRG->je; i++) {
    //tau = pow(10.0,-3.0 + pRG->dx3 * (Real)(i-1) * 10.0);
    dtau= 0.5 * pRG->dx1 * (pRG->R[0][i-1][0][0].chi + 
			    pRG->R[0][i  ][0][0].chi);
    tau += dtau;
    fprintf(fp,"%g %g %g\n",tau,pRG->R[pRG->ks][i][pRG->is][0].J,
	                        pRG->R[pRG->ks][i][pRG->is][0].S);
  }

  fclose(fp);
  return;
}

char *construct_filename(char *basename,char *key,int dump,char *ext)
{
  char *fname = NULL;

  int namelen = strlen(basename)+1+4+1+strlen(ext)+1;
  if (key != NULL) namelen += 1+strlen(key);

  if ((fname = (char *)calloc(namelen,sizeof(char))) == NULL)
    ath_error("[problem]: Error allocating memory for filename\n");

  if (key != NULL) {
    sprintf(fname, "%s-%s.%04d.%s", basename, key, dump, ext);
  } else {
    sprintf(fname, "%s.%04d.%s", basename, dump, ext);
  }

  return fname;
}


#endif /* RADIATION */
