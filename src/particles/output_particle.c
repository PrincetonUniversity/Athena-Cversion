#include "../copyright.h"
/*=============================================================================
FILE: output_particle.c
PURPOSE: contains all the routines necessary for outputting particles. There are
  two basic formats:
  1. Bin particles to the grid, then output particles as a grid.
  2. Dump the particle list directly.

  For particle binning, there can be many choices since particles may have different
  properties. We provide a default (and trivial) particle selection function, which
  select all the particles with different properties. The user can define their own
  particle selection functions in the problem generator and pass them to the main code.

  The output quantities include, density, momentum density and velocity of the selected
  particles averaged in one grid cell. The binned data are saved in arrays dpar, M1par,
  M2par and M3par, the expression function expr_??? are used to pick the relevant
  quantities, which is part of the output data structure. The way to output these binned
  particle quantities are then exactly the same as other gas quantities.

  Dumping particle list has not been developed yet since we need to figure out how to
  do visualization.

CONTAINS PUBLIC FUNCTIONS:
  void init_output_particle(Grid *pG);
  void particle_to_grid(Grid *pG, Domain *pD, Output *pout);
  void destruct_particle_grid();
  void dump_particle_binary(Grid *pG, Domain *pD);

History:
  Created:	Xuening Bai		Mar. 2009

==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "../globals.h"

#ifdef PARTICLES         /* endif at the end of the file */

float ***dpar;
float ***M1par;
float ***M2par;
float ***M3par;
int il,iu, jl,ju, kl,ku;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   expr_*
 *   property_all
 *============================================================================*/

Real expr_dpar (const Grid *pG, const int i, const int j, const int k);
Real expr_M1par(const Grid *pG, const int i, const int j, const int k);
Real expr_M2par(const Grid *pG, const int i, const int j, const int k);
Real expr_M3par(const Grid *pG, const int i, const int j, const int k);
Real expr_V1par(const Grid *pG, const int i, const int j, const int k);
Real expr_V2par(const Grid *pG, const int i, const int j, const int k);
Real expr_V3par(const Grid *pG, const int i, const int j, const int k);
int property_all(const int property);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/* Initialize particle output (for unbinned output format) */
void init_output_particle(Grid *pG)
{
  return;
}

/* Bin the particles to grid cells */
void particle_to_grid(Grid *pG, Domain *pD, Output *pout)
{
  int i, j, k, m1, m2, m3, Nx1T, Nx2T, Nx3T;
  long p;
  Real dx11, dx21, dx31, cellvol1, drho, a;
  Grain *gr;

  /* Get grid limit related quantities */
  cellvol1 = 1.0;

  if (pG->Nx1 > 1)  { dx11 = 1.0/pG->dx1; cellvol1 *= dx11; m1 = 1; }
  else { dx11 = 0.0; m1 = 0; }

  if (pG->Nx2 > 1)  { dx21 = 1.0/pG->dx2; cellvol1 *= dx21; m2 = 1; }
  else { dx21 = 0.0; m2 = 0; }

  if (pG->Nx3 > 1)  { dx31 = 1.0/pG->dx3; cellvol1 *= dx31; m3 = 1; }
  else { dx31 = 0.0; m3 = 0; }

  il = pG->is - m1*nghost;
  iu = pG->ie + m1*nghost;

  jl = pG->js - m2*nghost;
  ju = pG->je + m2*nghost;

  kl = pG->ks - m3*nghost;
  ku = pG->ke + m3*nghost;

  Nx1T = iu - il + 1;
  Nx2T = ju - jl + 1;
  Nx3T = ku - kl + 1;

  /* allocate memory for particle interpolation */
  dpar = (float***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(float));
  if (dpar == NULL) goto on_error;

  M1par = (float***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(float));
  if (M1par == NULL) goto on_error;

  M2par = (float***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(float));
  if (M1par == NULL) goto on_error;

  M3par = (float***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(float));
  if (M1par == NULL) goto on_error;

  /* initialization */
  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (i=il; i<=iu; i++) {
        dpar[k][j][i] = 0.0;
        M1par[k][j][i] = 0.0;
        M2par[k][j][i] = 0.0;
        M3par[k][j][i] = 0.0;
      }

  /* bin the particles */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);
    /* judge if the particle should be selected */
    if ((*(pout->par_prop))(gr->property)) {

      /* get grid index */
      if (dx11 > 0.0)
        celli(pG, gr->x1, dx11, &i, &a);
      else /* x1 dimension collapses */
        i = pG->is;
      if (dx21 > 0.0)
        cellj(pG, gr->x2, dx21, &j, &a);
      else /* x2 dimension collapses */
        j = pG->js;
      if (dx31 > 0.0)
        k = cellk(pG, gr->x3, dx31, &k, &a);
      else /* x3 dimension collapses */
        k = pG->ks;

      /* bin the particles to the grid */
#ifdef FEEDBACK
      drho = pG->grproperty[gr->property].m * cellvol1;
#else
      drho = cellvol1;
#endif
      dpar[k][j][i] += drho;
      M1par[k][j][i] += drho*gr->v1;
      M2par[k][j][i] += drho*gr->v2;
      M3par[k][j][i] += drho*gr->v3;
    }
  }

  return;

  on_error:
    ath_error("[init_output_particle]: Error allocating memory\n");
}

/* release particle grid memory */
void destruct_particle_grid()
{
  free_3d_array(dpar);
  free_3d_array(M1par);
  free_3d_array(M2par);
  free_3d_array(M3par);

  return;
}

/* dump unbinned particles in binary format */
void dump_particle_binary(Grid *pG, Domain *pD)
{
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* expr_*: where * are variables d,M1,M2,M3,V1,V2,V3 for particles */

Real expr_dpar(const Grid *pG, const int i, const int j, const int k) {
  if (dpar == NULL) ath_error("[expr_dpar]: Particles have not been binned for output, please set pargrid to 1.\n");
  return dpar[k][j][i];
}

Real expr_M1par(const Grid *pG, const int i, const int j, const int k) {
  if (M1par == NULL) ath_error("[expr_M1par]: Particles have not been binned for output, please set pargrid to 1.\n");
  return M1par[k][j][i];
}

Real expr_M2par(const Grid *pG, const int i, const int j, const int k) {
  if (M2par == NULL) ath_error("[expr_M2par]: Particles have not been binned for output, please set pargrid to 1.\n");
  return M2par[k][j][i];
}

Real expr_M3par(const Grid *pG, const int i, const int j, const int k) {
  if (M3par == NULL) ath_error("[expr_M3par]: Particles have not been binned for output, please set pargrid to 1.\n");
  return M3par[k][j][i];
}

Real expr_V1par(const Grid *pG, const int i, const int j, const int k) {
  if (dpar == NULL) ath_error("[expr_V1par]: Particles have not been binned for output, please set pargrid to 1.\n");
  if (dpar[k][j][i]>0.0)
    return M1par[k][j][i]/dpar[k][j][i];
  else return 0.0;
}

Real expr_V2par(const Grid *pG, const int i, const int j, const int k) {
  if (dpar == NULL) ath_error("[expr_V2par]: Particles have not been binned for output, please set pargrid to 1.\n");
  if (dpar[k][j][i]>0.0)
    return M2par[k][j][i]/dpar[k][j][i];
  else return 0.0;
}

Real expr_V3par(const Grid *pG, const int i, const int j, const int k) {
  if (dpar == NULL) ath_error("[expr_V3par]: Particles have not been binned for output, please set pargrid to 1.\n");
  if (dpar[k][j][i]>0.0)
    return M3par[k][j][i]/dpar[k][j][i];
  else return 0.0;
}

/* default choice for binning particles to the grid: 
   All the particles are binned, return true for any value.
*/
int property_all(const int property)
{
  return 1;  /* always true */
}

#endif /*PARTICLES*/
