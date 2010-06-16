#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_rad.c
 *
 * PURPOSE: Computes radiation variables (thermal source function,
 *          thermalization parameter, and total opacity) based on
 *          state of conserved variables after last integration
 *          step.  Currently uses simple constructions for testing
 *          purposes.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   hydro_to_rad()
 *============================================================================*/

#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"


#ifdef RADIATION

/*----------------------------------------------------------------------------*/
/* hydro_to_rad:  */

void hydro_to_rad(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k, ifr=0;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int ig,jg,kg,ioff,joff,koff;
  Real eps;

  if (pG->Nx[0] > 1) ioff = nghost - 1; else ioff = 0;
  if (pG->Nx[1] > 1) joff = nghost - 1; else joff = 0;
  if (pG->Nx[2] > 1) koff = nghost - 1; else koff = 0;

  eps = par_getd("problem","eps");
/* Compute radiation variables from conserved variables */
  for (k=ks; k<=ke; k++) {
    kg = k + koff;
  for (j=js; j<=je; j++) {
    jg = j + joff;
    for (i=is; i<=ie; i++) {
      ig = i + ioff;
      pRG->R[k][j][i][ifr].S = 1.0;
      pRG->R[k][j][i][ifr].B = 1.0;
      pRG->R[k][j][i][ifr].eps = eps;
      pRG->R[k][j][i][ifr].chi = pG->U[kg][jg][ig].d;
    }
  }
  }

  return;
}

#endif /* RADIATION */
