#include "../copyright.h"
/*==============================================================================
 * FILE: integrate.c
 *
 * PURPOSE: Contains public functions to set integrator.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_init()        - set pointer to integrate function based on dim
 *   integrate_destruct()    - call destruct integrate function based on dim
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"

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

/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){

  case 1:
    if(Nx1 <= 1) break;
    integrate_init_1d(Nx1);
#if defined(CTU_INTEGRATOR) || defined(VL_INTEGRATOR)
    return integrate_1d;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 1D problem");
#endif

  case 2:
    if(Nx3 > 1) break;
    integrate_init_2d(Nx1,Nx2);
#if defined(CTU_INTEGRATOR) || defined(VL_INTEGRATOR)
    return integrate_2d;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 2D problem");
#endif

  case 3:
    integrate_init_3d(Nx1,Nx2,Nx3);
#if defined(CTU_INTEGRATOR) || defined(VL_INTEGRATOR)
    return integrate_3d;
#else
    ath_err("[integrate_init]: Invalid integrator defined for 3D problem");
#endif

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
