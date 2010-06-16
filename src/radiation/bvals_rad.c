#include "../copyright.h"
/*=============================================================================
 * FILE: bvals_rad.c
 *
 * PURPOSE: Sets boundary condistions for radiative transfer on each edge of
 *          the grid.  It closely follows the methods and conventions and
 *          methods used for the hydro integration (see e.g. bvals_mhd.c).
 *
 *          Current version only rudimentary functions which handle specific
 *          problems and is a place-holder for more comprehensive methods
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   bvals_rad_1d_init()
 *   bvals_rad_2d_init() 
 *   bvals_rad_2d()
 */

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../prototypes.h"

#ifdef RADIATION

void bvals_rad_1d_init(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int ifr;

  for(ifr=0; ifr<pRG->nf; ifr++) {
    pRG->R[ks][js][is-1][ifr].S = 0.0;
    pRG->R[ks][js][is-1][ifr].eps = 
    pRG->R[ks][js][is  ][ifr].eps;
    pRG->R[ks][js][is-1][ifr].B = 1.0;
    pRG->R[ks][js][is-1][ifr].chi = 0.0;
    pRG->R[ks][js][ie+1][ifr].S = 1.0;
    pRG->R[ks][js][ie+1][ifr].eps = pRG->R[ks][js][ie][ifr].eps;
    pRG->R[ks][js][ie+1][ifr].B = 1.0;
    pRG->R[ks][js][ie+1][ifr].chi = pRG->R[ks][js][ie][ifr].chi;
  }

  return;
}

// FIX FOR MPI!
void bvals_rad_2d_init(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int i,j,ifr;

  for(ifr=0; ifr<pRG->nf; ifr++) {
/* set x2 boundaries */
    for(i=is-1; i<=ie+1; i++) {
      pRG->R[ks][js-1][i][ifr].S   = 0.0;
      pRG->R[ks][js-1][i][ifr].B   = 1.0;
      pRG->R[ks][js-1][i][ifr].eps = pRG->R[ks][js][i][ifr].eps;
      pRG->R[ks][js-1][i][ifr].chi = 0.0;
      pRG->R[ks][je+1][i][ifr].S   = 1.0;
      pRG->R[ks][je+1][i][ifr].B   = 1.0;
      pRG->R[ks][je+1][i][ifr].eps = pRG->R[ks][je][i][ifr].eps;
      pRG->R[ks][je+1][i][ifr].chi = pRG->R[ks][je][i][ifr].chi;
    }
/* set x1 boundaries */
    for(j=js; j<=je; j++) {
      pRG->R[ks][j][is-1][ifr].S   = pRG->R[ks][j][ie][ifr].S;
      pRG->R[ks][j][is-1][ifr].B   = pRG->R[ks][j][ie][ifr].B;
      pRG->R[ks][j][is-1][ifr].eps = pRG->R[ks][j][ie][ifr].eps;
      pRG->R[ks][j][is-1][ifr].chi = pRG->R[ks][j][ie][ifr].chi;
      pRG->R[ks][j][ie+1][ifr].S   = pRG->R[ks][j][is][ifr].S;
      pRG->R[ks][j][ie+1][ifr].B   = pRG->R[ks][j][is][ifr].B;
      pRG->R[ks][j][ie+1][ifr].eps = pRG->R[ks][j][is][ifr].eps;
      pRG->R[ks][j][ie+1][ifr].chi = pRG->R[ks][j][is][ifr].chi;
    }

  }
  return;
}

void bvals_rad_2d(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int i,j,ifr;

  for(ifr=0; ifr<pRG->nf; ifr++) 
/* update S at x1 boundaries */
    for(j=js; j<=je; j++) {
      pRG->R[ks][j][is-1][ifr].S = pRG->R[ks][j][ie][ifr].S;
      pRG->R[ks][j][ie+1][ifr].S = pRG->R[ks][j][is][ifr].S;
    }

  return;
}

#endif /* RADIATION */
