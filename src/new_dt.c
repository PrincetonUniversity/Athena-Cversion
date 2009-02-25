#include "copyright.h"
/*==============================================================================
 * FILE: new_dt.c
 *
 * PURPOSE: Computes timestep using CFL condition on cell-centered velocities
 *   and sound speed, and Alfven speed from face-centered B.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   new_dt()  - computes dt
 *   sync_dt() - synchronizes dt across all MPI patches
 *============================================================================*/

#include <stdio.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* new_dt:  */

void new_dt(Grid *pGrid)
{
  int i,j,k;
  Real di,v1,v2,v3,qsq,p,asq,cf1sq,cf2sq,cf3sq;
  Real max_v1, max_v2, max_v3, max_dti;
#ifdef MHD
  Real b1,b2,b3,bsq,tsum,tdif;
#endif /* MHD */
#ifdef PARTICLES
  long q;
#endif /* PARTICLES */
  Real nu, eta, min_dx;

  max_v1=0.0;	max_v2=0.0;	max_v3=0.0;
  max_dti = 0.0;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
  for (j=pGrid->js; j<=pGrid->je; j++) {
    for (i=pGrid->is; i<=pGrid->ie; i++) {
      di = 1.0/(pGrid->U[k][j][i].d);
      v1 = pGrid->U[k][j][i].M1*di;
      v2 = pGrid->U[k][j][i].M2*di;
      v3 = pGrid->U[k][j][i].M3*di;
      qsq = v1*v1 + v2*v2 + v3*v3;

#ifdef MHD

/* Use maximum of face-centered fields (always larger than cell-centered B) */
      b1 = pGrid->U[k][j][i].B1c 
        + fabs((double)(pGrid->B1i[k][j][i] - pGrid->U[k][j][i].B1c));
      b2 = pGrid->U[k][j][i].B2c 
        + fabs((double)(pGrid->B2i[k][j][i] - pGrid->U[k][j][i].B2c));
      b3 = pGrid->U[k][j][i].B3c 
        + fabs((double)(pGrid->B3i[k][j][i] - pGrid->U[k][j][i].B3c));
      bsq = b1*b1 + b2*b2 + b3*b3;
/* compute sound speed squared */
#ifdef ADIABATIC
      p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq
              - 0.5*bsq), TINY_NUMBER);
      asq = Gamma*p*di;
#elif defined ISOTHERMAL
      asq = Iso_csound2;
#endif /* EOS */
/* compute fast magnetosonic speed squared in each direction */
      tsum = bsq*di + asq;
      tdif = bsq*di - asq;
      cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3)*di));
      cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3)*di));
      cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2)*di));

#else /* MHD */

/* compute sound speed squared */
#ifdef ADIABATIC
      p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq),
              TINY_NUMBER);
      asq = Gamma*p*di;
#elif defined ISOTHERMAL
      asq = Iso_csound2;
#endif /* EOS */
/* compute fast magnetosonic speed squared in each direction */
      cf1sq = asq;
      cf2sq = asq;
      cf3sq = asq;

#endif /* MHD */

/* compute maximum cfl velocity (corresponding to minimum dt) */
      if (pGrid->Nx1 > 1)
        max_v1 = MAX(max_v1,fabs(v1)+sqrt((double)cf1sq));
      if (pGrid->Nx2 > 1)
        max_v2 = MAX(max_v2,fabs(v2)+sqrt((double)cf2sq));
      if (pGrid->Nx3 > 1)
        max_v3 = MAX(max_v3,fabs(v3)+sqrt((double)cf3sq));

    }
  }}

/* compute maximum cfl velocity with particles */
#ifdef PARTICLES
  for (q=0; q<pGrid->nparticle; q++) {
    if (pGrid->Nx1 > 1)
      max_v1 = MAX(max_v1, pGrid->particle[q].v1);
    if (pGrid->Nx2 > 1)
      max_v2 = MAX(max_v2, pGrid->particle[q].v2);
    if (pGrid->Nx3 > 1)
      max_v3 = MAX(max_v3, pGrid->particle[q].v3);
  }
#endif /* PARTICLES */

/* compute maximum inverse of dt (corresponding to minimum dt) */
  if (pGrid->Nx1 > 1)
    max_dti = MAX(max_dti, max_v1/pGrid->dx1);
  if (pGrid->Nx2 > 1)
    max_dti = MAX(max_dti, max_v2/pGrid->dx2);
  if (pGrid->Nx3 > 1)
    max_dti = MAX(max_dti, max_v3/pGrid->dx3);

/* new timestep.  Limit increase to 2x old value */
  if (pGrid->nstep == 0) {
    pGrid->dt = CourNo/max_dti;
  } else {
    pGrid->dt = MIN(2.0*pGrid->dt, CourNo/max_dti);
  }

#ifdef MPI_PARALLEL
  sync_dt(pGrid);
#endif /* MPI_PARALLEL */
  return;
}

/*----------------------------------------------------------------------------*/
/* sync_dt: uses MPI_Allreduce to compute minumum timestep over all MPI patches
 */

#ifdef MPI_PARALLEL
void sync_dt(Grid *pG)
{
  double dt, my_dt;
  int err;

  my_dt = pG->dt;

  err = MPI_Allreduce(&my_dt, &dt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[sync_dt]: MPI_Allreduce returned error code %d\n",err);

  pG->dt = dt;

  return;
}
#endif /* MPI_PARALLEL */
