#include "../copyright.h"
/*==============================================================================
 * FILE: diff_dt.c
 *
 * PURPOSE: Computes diffusion timestep using CFL condition, for all diffusive
 *   processes currently implemented in code.  These include:
 *     * Ohmic dissipation, Hall effect, ambipolar diffusion
 *     * Navier-Stokes and Braginskii viscosity
 *     * isotropic and anisotropic thermal conduction
 *   With MPI parallel jobs, finds minimum dt across all processors.
 *   Function returns minimum diffusion dt.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   diff_dt()  - computes dt
 *============================================================================*/

#include <stdio.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/* diff_dt:  */

Real diff_dt(MeshS *pM)
{
  int irefine, ir;
  Real dtmin_diffusion=(HUGE_NUMBER);
  Real dxmin,qa;
#ifdef MPI_PARALLEL
  double my_dt, dt;
  integer ierr;
#endif

/* Calculate minimum dx.  Always given by Grid on highest level of refinement */

  irefine = 1;
  for (ir=1; ir<(pM->NLevels); ir++) irefine *= 2;

  dxmin = pM->dx[0]/(Real)(irefine);
  if (pM->Nx[1] > 1) dxmin = MIN( dxmin, (pM->dx[1]/(Real)(irefine)) );
  if (pM->Nx[2] > 1) dxmin = MIN( dxmin, (pM->dx[2]/(Real)(irefine)) );

  qa = CourNo*(dxmin*dxmin)/2.0;
  if (pM->Nx[1] > 1) qa /= 2.0;
  if (pM->Nx[2] > 1) qa /= 2.0;

#ifdef THERMAL_CONDUCTION
  dtmin_diffusion = MIN(dtmin_diffusion,(qa/(kappa_iso + kappa_aniso)));
#endif
#ifdef VISCOSITY
  dtmin_diffusion = MIN(dtmin_diffusion,(qa/(nu_iso + nu_aniso)));
#endif

#ifdef RESISTIVITY
  if (eta_Ohm > 0.0) dtmin_diffusion = MIN(dtmin_diffusion,(qa/eta_Ohm));

/* FIX NEEDED: Hall timestep limit needs density and B */
  if (eta_Hall > 0.0) dtmin_diffusion = MIN(dtmin_diffusion,(qa/eta_Hall));

/* FIX NEEDED: AD timestep limit needs B */
  if (eta_AD > 0.0) dtmin_diffusion = MIN(dtmin_diffusion,(qa/eta_AD));
#endif

/* Find minimum timestep over all processors */
#ifdef MPI_PARALLEL
  my_dt = dtmin_diffusion;
  ierr = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  dtmin_diffusion = dt;
#endif /* MPI_PARALLEL */

  return dtmin_diffusion;
}
