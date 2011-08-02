#include "../copyright.h"
/*==============================================================================
 * FILE: formal_solution.c
 *
 * PURPOSE: Called from main loop before integrator.  Identifies and calls
 *          the appropriate formal solution algorithm.  If non-LTE will
 *          iteratively call formal_solution until convergence criterion
 *          is met.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution()  - interate formal solution until convergence
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER
Real ***sol;

void output_diag(Real *dsarr, int *isarr, int ns);
void get_solution_2d(RadGridS *pRG);
void max_dev(RadGridS *pRG, Real *dsm, int *ism);

void formal_solution(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  int i, ndim;
  Real *dSmax, dsm, dSmin, dSrmax;
  int  *isarr, ism, sflag;
#ifdef MPI_PARALLEL
  Real gdSrmax;
#endif
  if ((dSmax = (Real *)calloc(10000,sizeof(Real))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }
  if ((isarr = (int *)calloc(10000,sizeof(int))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }
  //if(lte == 0) sflag = 1; else sflag = 0;
  sflag = 1;
/* number of dimensions in Grid. */
  ndim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) ndim++;

  if (ndim == 1) {
/* compute formal solution with 1D method*/
    formal_solution_1d_init(pRG);
    for(i=0; i<niter; i++) {
      if(i > 0) bvals_rad(pD,sflag);
#ifdef RAD_MULTIG
      formal_solution_mg_1d(pRG,&dSrmax);
#else
      formal_solution_1d(pRG,&dSrmax);
#endif

/* User work (defined in problem()) */
    Userwork_in_formal_solution(pD);
   
/* Check whether convergence criterion is met. */
#ifdef MPI_PARALLEL
      MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      dSrmax=gdSrmax;
#endif
      if(dSrmax <= dScnv) {
	i++;
	break;
      }
    }
    formal_solution_1d_destruct();
  } else if (ndim == 2) {
    /*get_solution_2d(pRG); used for testing */
/* compute formal solution with 2D method*/
    formal_solution_2d_init(pRG);
    for(i=0; i<niter; i++) {
      if(i > 0) bvals_rad(pD,sflag);
#ifdef RAD_MULTIG
      formal_solution_mg_2d(pRG,&dSrmax);
#else
      formal_solution_2d(pRG,&dSrmax);
#endif

/* User work (defined in problem()) */
    Userwork_in_formal_solution(pD);

/* Check whether convergence criterion is met. */
#ifdef MPI_PARALLEL
      MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      dSrmax=gdSrmax;
#endif
      if(dSrmax <= dScnv) {
	i++;
	break;
      }
/* temporary output for debugging purpose */
/*      max_dev(pRG,&dsm,&ism);
      dSmax[i] = dsm;
      isarr[i] = ism;*/
    }
    formal_solution_2d_destruct();
  } else if (ndim == 3) {
/* compute formal solution with 3D method*/
    formal_solution_3d_init(pRG);
    for(i=0; i<niter; i++) {
      if(i > 0) bvals_rad(pD,sflag);
      formal_solution_3d(pRG,&dSrmax);

/* User work (defined in problem()) */
    Userwork_in_formal_solution(pD);

/* Check whether convergence criterion is met. */
#ifdef MPI_PARALLEL
      MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      dSrmax=gdSrmax;
#endif
      if(dSrmax <= dScnv) {
	i++;
	break;
      }
    }
    formal_solution_3d_destruct();
  }

  if (myID_Comm_world == 0) printf("iterations: %d, dSrmax: %g\n",i,dSrmax);
  if((i == niter) && (myID_Comm_world == 0)) 
    printf("Maximum number of iterations: niter=%d exceeded.\n",niter);

  /*  Used for testing purposes in old version of code */
  /*output_diag(dSmax,isarr,niter);
    if (sol != NULL) free_3d_array(sol); */

  /* free up memory */
  free(dSmax);
  free(isarr);

  return;
}


/* The following functions are for code testing purposes only */

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
	/*printf("%d %d %g %g\n",i,j,pRG->R[k][j][i][0].S,sol[k-pRG->ks][j-pRG->js][i-pRG->is]);
	*/
	if (dst > ds) {ds = dst; ixmax = i; iymax = j; izmax = k;}
      }

  (*dsm) = ds;
  (*ism) = iymax;
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

  tau=0.5 * pRG->dx2 * pRG->R[0][pRG->js-1][1][0].chi;
  for(j=pRG->js; j<=pRG->je; j++) {
    dtau= 0.5 * pRG->dx2 * (pRG->R[0][j-1][1][0].chi + 
			    pRG->R[0][j  ][1][0].chi);
    tau += dtau;
    /*    tau = pow(10.0,-3.0 + pRG->dx3 * (Real)(i-pRG->ks-0.5) * 10.0);
    */
    jmean = 1.0 - sqrt(3.0) * exp(-sqrt(3.0 * eps) * tau) /
      (sqrt(3.0) + sqrt(3.0 * eps));
     for(k=pRG->ks; k<=pRG->ke; k++)
      for(i=pRG->is; i<=pRG->ie; i++) 
	sol[k-pRG->ks][j-pRG->js][i-pRG->is] = eps + (1.0-eps) * jmean;
    /*sol[pRG->nx3-2-i][j-1] = eps + (1.0-eps) * jmean;
    */
  }

  return;

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

#endif /* RADIATION_TRANSFER */
