#include "copyright.h"
/*==============================================================================
 * FILE: integrate.c
 *
 * PURPOSE: Contains public functions to set integrator and source terms.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_init()        - set pointer to integrate function based on dim
 *   integrate_destruct()    - call destruct integrate function based on dim
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*----------------------------------------------------------------------------*/
/* integrate_init: initialize integrator; VGDFun_t defined in athena.h   */

VGDFun_t integrate_init(int Nx1, int Nx2, int Nx3)
{

/* Calculate the dimensions  */
  dim = 0;
  if(Nx1 > 1) dim++;
  if(Nx2 > 1) dim++;
  if(Nx3 > 1) dim++;

/* Since algorithm in lr_states is modified when VL integrator is defined, then
 * problem must be 3d */
#ifdef THREED_VL
  if(dim < 3) ath_error(
    "[integrate_init]: Problem must be 3D when VL integrator is defined");
#endif

/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){
  case 1:
    if(Nx1 <= 1) break;
    integrate_init_1d(Nx1);
    return integrate_1d;
  case 2:
    if(Nx3 > 1) break;
    integrate_init_2d(Nx1,Nx2);
    return integrate_2d;
  case 3:
    integrate_init_3d(Nx1,Nx2,Nx3);
    return THREE_D_INTEGRATOR;
  }

  if (dim == 1)
    ath_error("[integrate_init]: 1D problem must have Nx1 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",Nx1,Nx2,Nx3);
  if (dim == 2)
     ath_error("[integrate_init]: 2D problem must have Nx1 and Nx2 > 1: Nx1=%d, Nx2=%d, Nx3=%d\n",Nx1,Nx2,Nx3);

/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct:  */

void integrate_destruct()
{
  switch(dim){
  case 1:
    integrate_destruct_1d();
    return;
  case 2:
    integrate_destruct_2d();
    return;
  case 3:
    integrate_destruct_3d();
    return;
  }

  ath_error("[integrate_destruct]: Grid dimension = %d\n",dim);
}
