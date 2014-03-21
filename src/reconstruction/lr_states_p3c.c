#include "../copyright.h"
/*============================================================================*/
/*! \file lr_states_p3c.c
 *  \brief Third-order compact reconstruction in the primitive variables 
 *         using the Koren slope limiter or Cada & Torrilhon LimO3 limiter.
 *         Currently this works with the VL integrator in Cartesian coordinates. 
 *         For LimO3, a parameter r (RLIM) must be set and currently it is 0.1
 *         (more tests are required). Note that even with these limiters, 
 *         overall accuracy is 2nd-order as the time-integrator is 2nd-order.
 *         It is notable that LimO3 can preserve some local extrema.
 *         See Cada & Torrilhon, 2009, JCoPh, 228, 4118
 *         Configuration options: --with-order=3pck: Koren
 *                                --with-order=3pcl: LimO3
 *
 * PURPOSE: Third order (piecewise parabolic) compact spatial reconstruction in
 *   the primitive variables. Limiting is performed in the primitive variables.
 *
 * NOTATION: 
 * - W_{L,i-1/2} is reconstructed value on the left-side of interface at i-1/2
 * - W_{R,i-1/2} is reconstructed value on the right-side of interface at i-1/2
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 * - W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 * - W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - lr_states()          - computes L/R states
 * - lr_states_init()     - initializes memory for static global arrays
 * - lr_states_destruct() - frees memory for static global arrays	      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#if defined(THIRD_ORDER_COMPACT_PRIM_KOREN) || defined(THIRD_ORDER_COMPACT_PRIM_LIMO3)
#ifdef SPECIAL_RELATIVITY
#error: Third-prder compact reconstruction is incompatible with SR.
#endif
#ifdef CTU_INTEGRATOR
#error: Third-prder compact reconstruction is incompatible with CTU.
#endif
#if defined(SPHERICAL) || defined(CYLINDRICAL)
#error: Third-prder compact reconstruction can be used only in Cartesian coordinates.
#endif

/* Parameter for the LIMO3 parameter */
#define RLIM (0.1)

static Real **pW=NULL;

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states(const GridS *pG, const Prim1DS W[], const Real Bxc[],
 *               const Real dt, const Real dx, const int il, const int iu,
 *               Prim1DS Wl[], Prim1DS Wr[], const int dir)
 *  \brief Computes L/R states
 * Input Arguments:
 * - W = PRIMITIVE variables at cell centers along 1-D slice
 * - Bxc = B in direction of slice at cell center
 * - il,iu = lower and upper indices of zone centers in slice
 * W and Bxc must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 * - Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states(const GridS *pG __attribute((unused)),
               const Prim1DS W[], const Real Bxc[],
               const Real dt, const Real dx, const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[], 
               const int dir __attribute__((unused)))
{
  int i,n;
  Real dWl,dWr,dWm;
  Real pl,pr,pc;
  Real *pWl, *pWr;
  Real lmr,lml,theta,thetai,phir,phil,eta;

/* Set pointer to primitive variables */
  for (i=il-2; i<=iu+2; i++) pW[i] = (Real*)&(W[i]);


/*========================== START BIG LOOP OVER i =======================*/
  for (i=il-1; i<=iu+1; i++) {

    pWl = (Real *) &(Wl[i+1]);
    pWr = (Real *) &(Wr[i]);

    for (n=0; n<(NWAVE+NSCALARS); n++) {
/*--- Step 1. ------------------------------------------------------------------
 * Compute L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */
      pl = pW[i-1][n];
      pc = pW[i][n];
      pr = pW[i+1][n];
      dWl = pc - pl;
      dWr = pr - pc;
      theta=dWl/(fabs(dWr)+1e-40)*SIGN(dWr);
      thetai=dWr/(fabs(dWl)+1e-40)*SIGN(dWl);
#ifdef THIRD_ORDER_COMPACT_PRIM_KOREN
      phir=MIN(theta,(2.0+theta)/6.0);
      phir=MIN(phir,1.0);
      lmr=MAX(phir,0.0);
      phil=MIN(thetai,(2.0+thetai)/6.0);
      phil=MIN(phil,1.0);
      lml=MAX(phil,0.0);
#endif
#ifdef THIRD_ORDER_COMPACT_PRIM_LIMO3
      eta=(SQR(dWl)+SQR(dWr))/SQR(RLIM*dx);
      lmr=(2.0+theta)/6.0;
      lml=(2.0+thetai)/6.0;
      if(eta >= 1.0)
      {
        phir=MIN(theta,0.8);
        phir=MIN(lmr,phir);
        phir=MAX(-0.25*theta,phir);
        phir=MIN(lmr,phir);
        lmr=MAX(phir,0.0);
        phil=MIN(thetai,0.8);
        phil=MIN(lml,phil);
        phil=MAX(-0.25*thetai,phil);
        phil=MIN(lml,phil);
        lml=MAX(phil,0.0);
      }
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Set L/R values */
      pWl[n] = pc + lmr*dWr;
      pWr[n] = pc - lml*dWl;
    }

  } /*===================== END BIG LOOP OVER i ===========================*/

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_init(MeshS *pM)
 *  \brief Allocate enough memory for work arrays */

void lr_states_init(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd,n4v=4;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

  if ((pW = (Real**)malloc(nmax*sizeof(Real*))) == NULL) goto on_error;

  return;
  on_error:
    lr_states_destruct();
    ath_error("[lr_states_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void lr_states_destruct(void)
 *  \brief Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (pW != NULL) free(pW);
  return;
}

#endif /* THIRD_ORDER_COMPACT_PRIM */
