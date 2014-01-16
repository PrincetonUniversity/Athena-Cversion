#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_cooling.c
 *
 * PURPOSE: Contains public functions to integrate implicit cooling term
 *   using operator splitting.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_cooling() - calls functions for cooling
 *   integrate_cooling_init() - allocates memory for diff functions
 *   integrate_cooling_destruct() - frees memory for diff functions
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "prototypes.h"


/*----------------------------------------------------------------------------*/
/* integrate_cooling:  called in main loop, sets timestep and calls cooling functions
 */

void integrate_cooling(MeshS *pM)
{
  GridS *pG;
  int nl,nd;

/* Limit timestep by minimum for explicit update of diffusion operators */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
#ifdef OPERATOR_SPLIT_COOLING
        pG=pM->Domain[nl][nd].Grid;
        pG->dt = pM->dt;

/* Call cooling solver Mesh hierarchy. */
        cooling_solver(pG);
        pM->dt = MIN(pM->dt, pG->dt);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_cooling_init: Set function pointers for cooling operators */

void integrate_cooling_init(MeshS *pM)
{   
#ifdef OPERATOR_SPLIT_COOLING
#ifdef VL_INTEGRATOR
  if(pM->Nx[2]>1 && CourNo > 0.33333) ath_error("[integrate_cooling] CourNo should be smaller than 1/3 for 3D VL integrator with operator split cooling.\n");
#endif
  cooling_solver_init(pM);
#endif
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct:  Frees memory associated with cooling funcs  */

void integrate_cooling_destruct()
{
#ifdef OPERATOR_SPLIT_COOLING
  cooling_solver_destruct();
#endif
}
