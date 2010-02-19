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

#ifdef EXPLICIT_DIFFUSION
/* function pointers for diffusion operators (determined at runtime) */
static VDFun_t TConductF=NULL;
static VDFun_t ViscosityF=NULL;
static VDFun_t ResistivityF=NULL;

/* minimum timestep for explicit integration of diffusion operators */
static Real dtmin_diffusion;

/*----------------------------------------------------------------------------*/
/* integrate_diff:  called in main loop, sets timestep and calls
 *   appropriate functions for each diffusion operator
 */

void integrate_diff(MeshS *pM)
{
  GridS *pG;
  int nl,nd;

/* Limit timestep by minimum for explicit update of diffusion operators */

  pM->dt = MIN(pM->dt, dtmin_diffusion);;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;
        pG->dt = pM->dt;

/* Call diffusion operators across Mesh hierarchy.  Function pointers set in
 * integrate_diff_init().  TConductF must be called first to avoid an extra call
 * to set_bval_mhd().  */

        if (TConductF != NULL) (*TConductF)(&(pM->Domain[nl][nd]));

        if (ViscosityF != NULL) (*ViscosityF)(&(pM->Domain[nl][nd]));

        if (ResistivityF != NULL) (*ResistivityF)(&(pM->Domain[nl][nd]));
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_diff_init: Set function pointers for diffusion operators, set time
 * step minimum, call initialization routines to allocate memory
 */

void integrate_diff_init(MeshS *pM)
{   
  GridS *pG;
  int dim=0,nl,nd;
  Real dxmin, qa;

/* Calculate the dimensions  */
  if(pM->Nx[0] > 1) dim++;
  if(pM->Nx[1] > 1) dim++;
  if(pM->Nx[2] > 1) dim++;

  dtmin_diffusion = 0.0;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG=pM->Domain[nl][nd].Grid;
        dxmin = pG->dx1;
        if (pG->Nx[1] > 1) dxmin = MIN(dxmin,pG->dx2);
        if (pG->Nx[2] > 1) dxmin = MIN(dxmin,pG->dx3);
        qa = CourNo*(dxmin*dxmin)/4.0;

#if defined(ISOTROPIC_CONDUCTION)
        dtmin_diffusion = MIN(dtmin_diffusion,(qa/kappa_T));
#endif
#if defined(ANISOTROPIC_CONDUCTION)
        dtmin_diffusion = MIN(dtmin_diffusion,(qa/chi_C));
#endif

#if defined(NAVIER_STOKES) || defined(BRAGINSKII)
        dtmin_diffusion = MIN(dtmin_diffusion,(qa/nu_V));
#endif

#ifdef OHMIC
        dtmin_diffusion = MIN(dtmin_diffusion,(qa/eta_Ohm));
#endif        
/* I think the Hall timestep limit needs density */
#ifdef HALL_MHD
        dtmin_diffusion = MIN(dtmin_diffusion,(qa/(eta_Ohm + eta_Hall)));
#endif        
      }
    }
  }

/* Set function pointers for thermal conduciton, viscosity, and resistivity
 * based on dimension of problem, and macros set by configure.  Also check that
 * diffusion coefficients were set in problem generator.
 */

/* For isotropic thermal conduction the same function handles 1d/2d/3d */

#ifdef ISOTROPIC_CONDUCTION
  TConductF = isoconduct;
  isoconduct_init(pM);
  if (kappa_T <= 0.0) 
    ath_error("[diff_init] coefficent of conduction kappa_T was not set\n");
#endif
#ifdef ANISOTROPIC_CONDUCTION
  switch(dim){
    ath_error("[diff_init] anisotropic conduction requires 2D or 3D\n");
    break;
  case 2:
    TConductF = anisoconduct_2d;
    anisoconduct_init(pM);
    break;
  case 3:
    TConductF = anisoconduct_3d;
    anisoconduct_init(pM);
    break;
  }
  if (chi_C <= 0.0) 
    ath_error("[diff_init] coefficent of conduction chi_C was not set\n");
#endif

/* viscosity */

#ifdef NAVIER_STOKES
  switch(dim){
  case 1:
    ViscosityF = ns_viscosity_1d;
    ns_viscosity_init(pM);
    break;
  case 2:
    ViscosityF = ns_viscosity_2d;
    ns_viscosity_init(pM);
    break;
  case 3:
    ViscosityF = ns_viscosity_3d;
    ns_viscosity_init(pM);
    break;
  }
  if (nu_V <= 0.0) 
    ath_error("[diff_init] coefficent of viscosity nu_V was not set\n");
#endif
#ifdef BRAGINSKII
  switch(dim){
  case 1:
    ath_error("[diff_init] Braginskii viscosity requires 2D or 3D\n");
    break;
  case 2:
    ViscosityF = brag_viscosity_2d;
    brag_viscosity_init(pM);
    break;
  case 3:
    ViscosityF = brag_viscosity_3d;
    brag_viscosity_init(pM);
    break;
  }
  if (nu_V <= 0.0) 
    ath_error("[diff_init] coefficent of viscosity nu_V was not set\n");
#endif

/* resistivity */

#ifdef OHMIC
  switch(dim){
  case 1:
    ResistivityF = ohmic_resistivity_1d;
    ohmic_resistivity_init(pM);
    break;
  case 2:
    ResistivityF = ohmic_resistivity_2d;
    ohmic_resistivity_init(pM);
    break;
  case 3:
    ResistivityF = ohmic_resistivity_3d;
    ohmic_resistivity_init(pM);
    break;
  }
  if (eta_Ohm <= 0.0) 
    ath_error("[diff_init] coefficent of resistivity eta_Ohm was not set\n");
#endif
#ifdef HALL_MHD
  switch(dim){
  case 1:
    ResistivityF = hall_resistivity_1d;
    hall_resistivity_init(pM);
    break;
  case 2:
    ResistivityF = hall_resistivity_2d;
    hall_resistivity_init(pM);
    break;
  case 3:
    ResistivityF = hall_resistivity_3d;
    hall_resistivity_init(pM);
    break;
  }
  if (eta_Hall <= 0.0) 
    ath_error("[diff_init] coefficent of resistivity eta_Hall was not set\n");
#endif

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct:  Frees memory associated with diffusion funcs  */

void integrate_diff_destruct()
{

#ifdef NAVIER_STOKES
  ns_viscosity_destruct();
#endif
#ifdef BRAGINSKII
  brag_viscosity_destruct();
#endif
#ifdef OHMIC
  ohmic_resistivity_destruct();
#endif
#ifdef HALL_MHD
  hall_resistivity_destruct();
#endif
#ifdef ISOTROPIC_CONDUCTION
  isoconduct_destruct();
#endif
#ifdef ANISOTROPIC_CONDUCTION
  anisoconduct_destruct();
#endif
}
#endif /* EXPLICIT_DIFFUSION */
