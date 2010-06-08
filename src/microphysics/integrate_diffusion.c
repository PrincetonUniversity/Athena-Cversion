#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_diffusion.c
 *
 * PURPOSE: Contains public functions to integrate explicit diffusion terms
 *   using operator splitting.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_diff() - calls functions for each diffusion operator
 *   integrate_diff_init() - allocates memory for diff functions
 *   integrate_diff_destruct() - frees memory for diff functions
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* integrate_diff:  called in main loop, sets timestep and/or orchestrates
 * subcycling, calls appropriate functions for each diffusion operator
 */

void integrate_diff(MeshS *pM)
{
  GridS *pG;
  int nl,nd;
  Real dtmin_expl;

  dtmin_expl = diff_dt(pM);

/* Limit timestep by minimum for explicit update of diffusion operators.
 * Currently, subcycling is not implemented!! */

  pM->dt = MIN(pM->dt, dtmin_expl);;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;
        pG->dt = pM->dt;

/* Call diffusion operators across Mesh hierarchy.
 * Conduction must be called first to avoid an extra call to bval_mhd().  */

#ifdef THERMAL_CONDUCTION
        conduction(&(pM->Domain[nl][nd]));
#endif

#ifdef RESISTIVITY
        resistivity(&(pM->Domain[nl][nd]));
#endif

#ifdef VISCOSITY
        viscosity(&(pM->Domain[nl][nd]));
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_diff_init: call functions to allocate memory
 */

void integrate_diff_init(MeshS *pM)
{   
/* Check that diffusion coefficients were set in problem generator, call memory
 * allocation routines.  */

#ifdef THERMAL_CONDUCTION
  if ((kappa_iso + kappa_aniso) <= 0.0) 
    ath_error("[diff_init] coefficents of thermal conduction not set\n");
  conduction_init(pM);
#endif

#ifdef VISCOSITY
  if ((nu_iso + nu_aniso) <= 0.0) 
    ath_error("[diff_init] coefficents of viscosity not set\n");
  viscosity_init(pM);
#endif

#ifdef RESISTIVITY
  if ((eta_Ohm + eta_Hall + eta_AD) <= 0.0) 
    ath_error("[diff_init] coefficents of resistivity not set\n");
  resistivity_init(pM);
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct:  Frees memory associated with diffusion funcs  */

void integrate_diff_destruct()
{
#ifdef THERMAL_CONDUCTION
  conduction_destruct();
#endif
#ifdef RESISTIVTY
  resistivity_destruct();
#endif
#ifdef VISCOSITY
  viscosity_destruct();
#endif
}
