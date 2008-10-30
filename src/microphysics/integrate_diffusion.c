#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_diffusion.c
 *
 * PURPOSE: Contains public functions to integrate explicit diffusion terms
 *   using operator splitting.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_explicit_diff() - calls functions for each diffusion operator
 *   integrate_explicit_diff_init() - allocates memory for diff functions
 *   integrate_explicit_diff_destruct() - frees memory for diff functions
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#include "prototypes.h"

/* function pointers for diffusion operators (determined at runtime) */
static VGDFun_t ApplyViscosity=NULL;
static VGDFun_t ApplyResistivity=NULL;
static VGDFun_t ApplyThermalConduct=NULL;

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*----------------------------------------------------------------------------*/
/* integrate_explicit_diff:  called in main loop, sets timestep and calls
 *   appropriate functions for each diffusion operator
 */

void integrate_explicit_diff(Grid *pGrid, Domain *pDomain)
{
  Real min_dx,min_dt=0.0;

/* Calculate explicit diffusion timestep */

  min_dx = pGrid->dx1;
  if (pGrid->Nx2 > 1) min_dx = MIN(min_dx,pGrid->dx2);
  if (pGrid->Nx3 > 1) min_dx = MIN(min_dx,pGrid->dx3);
#ifdef OHMIC
  pGrid->dt = MIN(pGrid->dt,CourNo*((min_dx*min_dx)/(4.0*eta_Ohm)));
#endif        
/* I think the Hall timestep limit needs density */
#ifdef HALL_MHD
  pGrid->dt = MIN(pGrid->dt,CourNo*((min_dx*min_dx)/(4.0*(eta_Ohm + eta_Hall))));
#endif        
#if defined(NAVIER_STOKES) || defined(BRAGINSKII)
  pGrid->dt = MIN(pGrid->dt,CourNo*((min_dx*min_dx)/(4.0*nu_V)));
#endif
#if defined(ISOTROPIC_CONDUCTION)
  pGrid->dt = MIN(pGrid->dt,CourNo*((min_dx*min_dx)/(4.0*kappa_T)));
#endif
#if defined(ANISOTROPIC_CONDUCTION)
  pGrid->dt = MIN(pGrid->dt,CourNo*((min_dx*min_dx)/(4.0*chi_C)));
#endif

/* Call diffusion operators.  Function pointers set in 
 * integrate_explicit_diff_init() */

  if (ApplyViscosity != NULL) (*ApplyViscosity)(pGrid, pDomain);

  if (ApplyResistivity != NULL) (*ApplyResistivity)(pGrid, pDomain);

  if (ApplyThermalConduct != NULL) (*ApplyThermalConduct)(pGrid, pDomain);

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_explicit_diff_init: Set function pointers for diffusion 
 *   operators, call initialization routines to allocate memory
 */

void integrate_explicit_diff_init(Grid *pGrid, Domain *pDomain)
{   
/* Calculate the dimensions  */
  dim=0;
  if(pGrid->Nx1 > 1) dim++;
  if(pGrid->Nx2 > 1) dim++;
  if(pGrid->Nx3 > 1) dim++;

/* Set function pointers for viscosity, resistivity, and anisotropic conduction
 * based on dimension of problem, and macros set by configure 
 */

  switch(dim){

/* 1D (braginskii viscosity and anisotropic conduction not allowed) */
  case 1:
#ifdef NAVIER_STOKES
    ApplyViscosity = ns_viscosity_1d;
    ns_viscosity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef OHMIC
    ApplyResistivity = ohmic_resistivity_1d;
    ohmic_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef HALL_MHD
    ApplyResistivity = hall_resistivity_1d;
    hall_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
    break;

/* 2D */
  case 2:
#ifdef NAVIER_STOKES
    ApplyViscosity = ns_viscosity_2d;
    ns_viscosity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef BRAGINSKII
    ApplyViscosity = brag_viscosity_2d;
    brag_viscosity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef OHMIC
    ApplyResistivity = ohmic_resistivity_2d;
    ohmic_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef HALL_MHD
    ApplyResistivity = hall_resistivity_2d;
    hall_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef ANISOTROPIC_CONDUCTION
    ApplyThermalConduct = anisoconduct_2d;
    anisoconduct_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
    break;

/* 3D */
  case 3:
#ifdef NAVIER_STOKES
    ApplyViscosity = ns_viscosity_3d;
    ns_viscosity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef BRAGINSKII
    ApplyViscosity = brag_viscosity_3d;
    brag_viscosity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef OHMIC
    ApplyResistivity = ohmic_resistivity_3d;
    ohmic_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef HALL_MHD
    ApplyResistivity = hall_resistivity_3d;
    hall_resistivity_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
#ifdef ANISOTROPIC_CONDUCTION
    ApplyThermalConduct = anisoconduct_3d;
    anisoconduct_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif
    break;
  }

/* Set function pointer for isotropic thermal diffusion operator (the same
 * function handles 1d/2d/3d)
 */

#ifdef ISOTROPIC_CONDUCTION
  ApplyThermalConduct = isoconduct;
  isoconduct_init(pGrid->Nx1, pGrid->Nx2, pGrid->Nx3);
#endif

  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct:  Frees memory associated with diffusion funcs  */

void integrate_explicit_diff_destruct()
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
