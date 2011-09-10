#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Contains public functions to set BackEuler.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   BackEuler_init()        - set pointer to integrate function based on dim
 *   BackEuler_destruct()    - call destruct integrate function based on dim
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)

/* dimension of calculation (determined at runtime) */
static int dim=0;


/*----------------------------------------------------------------------------*/
/* BackEuler_init: initialize BackEuler; VMFun_t defined in athena.h   */

VMFun_t BackEuler_init(MeshS *pM)
{
  int i;
 
/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pM->Nx[i] > 1) dim++;



/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){

  case 1:
    if(pM->Nx[0] <= 1) break;


    BackEuler_init_1d(pM);
    return BackEuler_1d;

  case 2:
    if(pM->Nx[2] > 1) break;

    BackEuler_init_2d(pM);

    return BackEuler_2d;	

  case 3:
    
    BackEuler_init_3d(pM);

    return BackEuler_3d;

  }

 
/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/* BackEuler_destruct:  */

void BackEuler_destruct()
{
  switch(dim){
  case 1:
    BackEuler_destruct_1d();
    return;
  case 2:
    BackEuler_destruct_2d();
    return;
  case 3:
    BackEuler_destruct_3d();
    return;
  }

  ath_error("[integrate_destruct]: Grid dimension = %d\n",dim);
}


#endif

