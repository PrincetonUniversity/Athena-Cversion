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

void formal_solution(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  int i, ndim;
  int ifr, nf = pRG->nf, nfc = 0;
  Real *dSmaxa, dSmax = 0.0, dsm, dSmin, dSrmax;
  int  ism, *iconv;
#ifdef MPI_PARALLEL
  Real gdSrmax;
#endif
  if ((dSmaxa = (Real *)calloc(nf,sizeof(Real))) == NULL) {
     ath_error("[get_solution]: Error allocating memory\n");
   }

  if ((iconv = (int *)calloc(nf,sizeof(int))) == NULL) {
    ath_error("[get_solution]: Error allocating memory\n");
  }

/* set iconv = 0 indicating non-convergence */
  for(ifr=0; ifr<nf; ifr++) 
    iconv[ifr] = 0;

/* number of dimensions in Grid. */
  ndim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) ndim++;

  if (ndim == 1) {
/* compute formal solution with 1D method*/
    formal_solution_1d_init(pRG);
    for(i=0; i<niter; i++) {
/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;
      
      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,ifr);
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
    formal_solution_1d_destruct();
  } else if (ndim == 2) {
/* compute formal solution with 2D method*/
    formal_solution_2d_init(pRG);
    for(i=0; i<niter; i++) {
/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;

      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,ifr);
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
    formal_solution_2d_destruct();
  } else if (ndim == 3) {
/* compute formal solution with 3D method*/
    formal_solution_3d_init(pRG);
    for(i=0; i<niter; i++) {
/* break out of loop if all frequencies are converged */
      if (nfc == nf) break;

      for(ifr=0; ifr<nf; ifr++) {
/* perform iteration only if frequency has not converged */
	if(iconv[ifr] != 1) {
	  bvals_rad(pD,ifr,ifr);
	  formal_solution_3d(pRG,&dSrmax,ifr);
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
    formal_solution_3d_destruct();
  }

  if (myID_Comm_world == 0) {
    for(ifr=0; ifr<nf; ifr++)
      if(dSmax < dSmaxa[ifr]) dSmax=dSmaxa[ifr];
    printf("iterations=%d, dSmax=%g\n",i,dSmax);
  }
  /* free up memory */
  free(iconv);
  free(dSmaxa);

  return;
}



#endif /* RADIATION_TRANSFER */
