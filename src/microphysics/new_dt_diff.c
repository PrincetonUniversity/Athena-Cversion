#include "../copyright.h"
/*============================================================================*/
/*! \file new_dt_diff.c
 *  \brief Computes stability constraint on timestep for all diffusive
 *   processes currently implemented in code.
 *  
 *  These include:
 *     - Ohmic dissipation, Hall effect, ambipolar diffusion
 *     - Navier-Stokes and Braginskii viscosity
 *     - isotropic and anisotropic thermal conduction
 *  The function returns maximum inverse of dt for all Domains at all Levels.
 *  With MPI, this value is calculated only for the Grids being updated on
 *  this processor.  The calling function new_dt() is responsible for finding
 *  the maximum over all processors.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - new_dt_diff()  - computes maximum inverse of dt */
/*============================================================================*/

#include <stdio.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/*----------------------------------------------------------------------------*/
/*! \fn Real diff_dt(MeshS *pM)
 *  \brief Computes diffusion timestep */
Real new_dt_diff(MeshS *pM)
{
  Real max_dti_diff=(TINY_NUMBER);
  Real dxmin,qa;
#ifdef RESISTIVITY
  int i,j,k,nl,nd;
  GridS *pG;

/* Calculate the magnetic diffusivity array */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {

        pG=pM->Domain[nl][nd].Grid;

        get_eta(pG);
      }
    }
  }
#endif

/* Calculate minimum dx.  Always given by Grid on highest level of refinement */

  dxmin = pM->dx[0]/pow(2,((pM->NLevels)-1));
  if (pM->Nx[1] > 1) dxmin = MIN( dxmin, (pM->dx[1]/pow(2,((pM->NLevels)-1))) );
  if (pM->Nx[2] > 1) dxmin = MIN( dxmin, (pM->dx[2]/pow(2,((pM->NLevels)-1))) );

  qa = (dxmin*dxmin)/4.0;
  if (pM->Nx[1] > 1) qa = (dxmin*dxmin)/8.0;
  if (pM->Nx[2] > 1) qa = (dxmin*dxmin)/6.0;

#ifdef THERMAL_CONDUCTION
  max_dti_diff = MAX( max_dti_diff, ((kappa_iso + kappa_aniso)/qa) );
#endif
#ifdef VISCOSITY
  max_dti_diff = MAX( max_dti_diff, ((nu_iso + nu_aniso)/qa) );
#endif

#ifdef RESISTIVITY
/* Since resistivities can vary from cell to cell, must loop over all cells */
  for (nl=pM->NLevels-1; nl>=0; nl--){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pG = pM->Domain[nl][nd].Grid;

        dxmin = pG->dx1;
        if (pG->Nx[1] > 1) dxmin = MIN( dxmin, (pG->dx2) );
        if (pG->Nx[2] > 1) dxmin = MIN( dxmin, (pG->dx3) );

        qa = (dxmin*dxmin)/4.0;
        if (pG->Nx[1] > 1) qa = (dxmin*dxmin)/8.0;
        if (pG->Nx[2] > 1) qa = (dxmin*dxmin)/6.0;

        for (k=pG->ks; k<=pG->ke; k++) {
        for (j=pG->js; j<=pG->je; j++) {
        for (i=pG->is; i<=pG->ie; i++) {

          max_dti_diff = MAX( max_dti_diff, ((pG->eta_Ohm[k][j][i] +
                              pG->eta_AD[k][j][i])/qa) );
  
        }}}
        if (Q_Hall > 0.0) {
          for (k=pG->ks; k<=pG->ke; k++) {
          for (j=pG->js; j<=pG->je; j++) { 
          for (i=pG->is; i<=pG->ie; i++) {

            max_dti_diff = MAX( max_dti_diff, fabs(pG->eta_Hall[k][j][i])/qa);

          }}}
        }
      }
    }
  }
#endif

  return max_dti_diff;
}
