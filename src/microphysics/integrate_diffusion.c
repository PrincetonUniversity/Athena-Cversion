#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_diffusion.c
 *  \brief Contains public functions to integrate explicit diffusion terms
 *   using operator splitting.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_diff() - calls functions for each diffusion operator
 * - integrate_diff_init() - allocates memory for diff functions
 * - integrate_diff_destruct() - frees memory for diff functions */
/*============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_diff(MeshS *pM)
 *  \brief Called in main loop, sets timestep and/or orchestrates
 * subcycling, calls appropriate functions for each diffusion operator
 */

void integrate_diff(MeshS *pM)
{
  GridS *pG;
  int nl,nd;

/* Call diffusion operators across Mesh hierarchy.
 * Conduction must be called first to avoid an extra call to bval_mhd().  */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;

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
/*! \fn void integrate_diff_init(MeshS *pM)
 *  \brief Call functions to allocate memory
 */

void integrate_diff_init(MeshS *pM __attribute__ ((unused)))
{   
/* Check that diffusion coefficients were set in problem generator, call memory
 * allocation routines.  */

#ifdef THERMAL_CONDUCTION
  if ((kappa_iso+kappa_aniso)==0.0 || kappa_iso<0.0 || kappa_aniso<0.0) 
    ath_error("[diff_init] problem with coefficents of thermal conduction\n");
  conduction_init(pM);
#endif

#ifdef VISCOSITY
  if ((nu_iso+nu_aniso)==0.0 || nu_iso<0.0 || nu_aniso<0.0) 
    ath_error("[diff_init] problem with coefficents of viscosity\n");
  viscosity_init(pM);
#endif

#ifdef RESISTIVITY
  if ((eta_Ohm+Q_Hall+Q_AD)==0.0 || eta_Ohm<0.0 || Q_Hall<0.0 || Q_AD<0.0) 
    ath_error("[diff_init] problem with coefficents of resistivity\n");
  resistivity_init(pM);
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_diff_destruct()
 *  \brief Frees memory associated with diffusion funcs  */
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
