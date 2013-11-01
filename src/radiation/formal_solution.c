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

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn formal_solution(DomainS *pD)
 * Controls computation of formal solution 
 */
void formal_solution(DomainS *pD, const int outflag, const int fstflag)
{

  RadGridS *pRG;
  int niter, ndim, nf;
  int i, ifr, nfc = 0;
  Real *dSmaxa, dSmax = 0.0, dsm, dSmin, dSrmax;
  Real dScnv;
  int  ism, *iconv;
  
#ifdef MPI_PARALLEL
  Real gdSrmax;
#endif


  if (outflag == 0)  {
    pRG = pD->RadGrid;    /* set ptr to RadGrid */
    niter = par_geti("radiation","niter");
    dScnv = par_getd("radiation","dScnv");
    if (fstflag == 1)  {
      niter = par_geti_def("radiation","niter0",niter);
      dScnv = par_getd_def("radiation","dScnv0",dScnv);
    }
  } else { 
    pRG = pD->RadOutGrid; /* set ptr to RadOutGrid */
    niter = par_geti("radiation_output","niter");
    dScnv = par_getd("radiation_output","dScnv");
    printf("fs: %d %g\n",niter,dScnv);
    if (fstflag == 1)  {
      niter = par_geti_def("radiation_output","niter0",niter);
      dScnv = par_getd_def("radiation_output","dScnv0",dScnv);
    }
  }

/* allocate memory for convergence arrays */
  nf = pRG->nf;
  if ((dSmaxa = (Real *)calloc(nf,sizeof(Real))) == NULL) {
     ath_error("[formal_solution]: Error allocating memory\n");
   }

  if ((iconv = (int *)calloc(nf,sizeof(int))) == NULL) {
    ath_error("[formal_solution]: Error allocating memory\n");
  }

/* set iconv = 0 indicating non-convergence */
  for(ifr=0; ifr<nf; ifr++) 
    iconv[ifr] = 0;

/* number of dimensions in Grid. */
  ndim=1;
  for (i=1; i<3; i++) if (pD->Nx[i]>1) ndim++;

/* ----------------------- Ray Tracing Block --------------------- */
#ifdef RAY_TRACING
/* if ray tracing opacity functions not defined in problem generator
 * assume they are equivalent to the radiation transfer routine
 * functions */
  if (get_raytrace_thermal_fraction == NULL)
    get_raytrace_thermal_fraction = get_thermal_fraction;
  if (get_raytrace_opacity == NULL)
    get_raytrace_opacity = get_total_opacity;
/* frequency conversion not set, use default function defined in raytrace.c */
  if (raytrace_to_radtrans == NULL)
    raytrace_to_radtrans = raytrace_to_radtrans_default;
/* If enabled, compute contribution of external radiation. Note that we pass
 * the Domain pD as pGrid data is needed by the ray trace function to compute
 * opacities (which are stored for standard radiative transfer). */
  ray_trace(pD,outflag);
#endif

  if (ndim == 1) {
/* compute formal solution with 1D method*/
    for(i=0; i<niter; i++) {
/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;
      
      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,outflag);
	  formal_solution_1d(pRG,&dSrmax,ifr);
/* find global maximum relative change */
#ifdef MPI_PARALLEL
	  MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  dSrmax=gdSrmax;
#endif
	  dSmaxa[ifr] = dSrmax;
/* Check whether convergence criterion is met. */   
	  if(dSrmax <= dScnv) {
	    i++;
	    nfc++;
	    iconv[ifr] = 1;
	  }
	}
      }
/* User work (defined in problem()) */
      Userwork_in_formal_solution(pD);
    }
  } else if (ndim == 2) {
/* compute formal solution with 2D method*/
    for(i=0; i<niter; i++) { 
/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;

      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,outflag);
	  formal_solution_2d(pRG,&dSrmax,ifr);
/* find global maximum relative change */
#ifdef MPI_PARALLEL
	  MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  dSrmax=gdSrmax;
#endif
	  dSmaxa[ifr] = dSrmax;
/* Check whether convergence criterion is met. */
	  if(dSrmax <= dScnv) {
	    i++;
	    nfc++;
	    iconv[ifr] = 1;
	  }
	}
      }
  
/* User work (defined in problem()) */
      Userwork_in_formal_solution(pD);
    }
  } else if (ndim == 3) {
/* compute formal solution with 3D method*/
    for(i=0; i<niter; i++) {

/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;

      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,outflag);
	  formal_solution_3d(pRG,&dSrmax,ifr);
/* find global maximum relative change */
#ifdef MPI_PARALLEL
	  MPI_Allreduce(&dSrmax, &gdSrmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	  dSrmax=gdSrmax;
#endif	  
	  if (myID_Comm_world == 0) printf("iter: %d %g\n",i,dSrmax);
	  dSmaxa[ifr] = dSrmax;
/* Check whether convergence criterion is met. */
	  if(dSrmax <= dScnv) {
	    i++;
	    nfc++;
	    iconv[ifr] = 1;
	  }
	}
      }
/* User work (defined in problem()) */
      Userwork_in_formal_solution(pD);
    }
  }

  if (myID_Comm_world == 0) {
    for(ifr=0; ifr<nf; ifr++)
      if(dSmax < dSmaxa[ifr]) dSmax=dSmaxa[ifr];
    if (outflag == 0)
      printf("integration grid: ");
    else 
      printf("output grid: ");
    printf("iterations=%d, dSmax=%g\n",i,dSmax);
  }
  /* free up memory */
  free(iconv);
  free(dSmaxa);

/* User work (defined in problem()) */
  Userwork_after_formal_solution(pD);

  return;
}



#endif /* RADIATION_TRANSFER */
