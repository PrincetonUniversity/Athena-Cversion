#include "copyright.h"
/*==============================================================================
 * FILE: init_grid_block.c
 *
 * PURPOSE: Initializes most variables in the Grid structure:
 *      time,nstep,[ijk]s,[ijk]e,dx[123],[ijk]disp,x[123]_0
 *   For MPI jobs, calls domain partitioning function.  Allocates memory for
 *   gas arrays and interface B
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_grid_block()
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS: none
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

void init_grid_block(Grid *pGrid)
{
  int Nx1T,Nx2T,Nx3T;    /* Total Number of grid cells in x1,x2,x3 direction */
  Real x1min,x1max,x2min,x2max,x3min,x3max;   /* read from input file */

  pGrid->time = 0.0;
  pGrid->nstep = 0;

/* ---------------------  Intialize grid in 1-direction --------------------- */
/* Initialize is,ie */

  pGrid->Nx1 = par_geti("grid","Nx1");
  if(pGrid->Nx1 > 1) {
    pGrid->is = nghost;
    pGrid->ie = pGrid->Nx1 + nghost - 1;
  }
  else if(pGrid->Nx1 == 1) {
    pGrid->is = pGrid->ie = 0;
  }
  else {
    ath_error("[init_grid_block]: Nx1 = %d must be >= 1\n",pGrid->Nx1);
  }

/* Compute dx1 */

  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  if(x1max < x1min) {
    ath_error("[init_grid_block]: x1max = %g < x1min = %g\n",x1max,x1min);
  }

  pGrid->dx1 = (x1max - x1min)/(Real)pGrid->Nx1;

/* Initialize i-displacement, and the x1-position of coordinate ix = 0. */

  pGrid->idisp = -pGrid->is;   /* Index i = is has Coordinate ix = 0 */
  pGrid->x1_0 = x1min;         /* Coordinate ix = 0 is at position x1min */

/* ---------------------  Intialize grid in 2-direction --------------------- */
/* Initialize js,je */

  pGrid->Nx2 = par_geti("grid","Nx2");
  if(pGrid->Nx2 > 1) {
    pGrid->js = nghost;
    pGrid->je = pGrid->Nx2 + nghost - 1;
  }
  else if(pGrid->Nx2 == 1) {
    pGrid->js = pGrid->je = 0;
  }
  else {
    ath_error("[init_grid_block]: Nx2 = %d must be >= 1\n",pGrid->Nx2);
  }

/* Compute dx2 */

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  if(x2max < x2min) {
    ath_error("[init_grid_block]: x2max = %g < x2min = %g\n",x2max,x2min);
  }
  pGrid->dx2 = (x2max - x2min)/(Real)pGrid->Nx2;

/* Initialize j-displacement, and the x2-position of coordinate jx = 0. */

  pGrid->jdisp = -pGrid->js;   /* Index j = js has Coordinate jx = 0 */
  pGrid->x2_0 = x2min;         /* Coordinate jx = 0 is at position x2min */

/* ---------------------  Intialize grid in 3-direction --------------------- */
/* Initialize ks,ke */

  pGrid->Nx3 = par_geti("grid","Nx3");
  if(pGrid->Nx3 > 1) {
    pGrid->ks = nghost;
    pGrid->ke = pGrid->Nx3 + nghost - 1;
  }
  else if(pGrid->Nx3 == 1) {
    pGrid->ks = pGrid->ke = 0;
  }
  else {
    ath_error("[init_grid_block]: Nx3 = %d must be >= 1\n",pGrid->Nx3);
  }

/* Compute dx3 */

  x3min = par_getd("grid","x3min");
  x3max = par_getd("grid","x3max");
  if(x3max < x3min) {
    ath_error("[init_grid_block]: x3max = %g < x3min = %g\n",x3max,x3min);
  }
  pGrid->dx3 = (x3max - x3min)/(Real)pGrid->Nx3;

/* Initialize k-displacement, and the x3-position of coordinate kx = 0. */

  pGrid->kdisp = -pGrid->ks;   /* Index k = ks has Coordinate kx = 0 */
  pGrid->x3_0 = x3min;         /* Coordinate kx = 0 is at position x3min */

/* -------- For MPI Parallel Calculations Partition the Grid Domain -------- */

#ifdef MPI_PARALLEL
  domain_partition(pGrid);
#endif /* MPI_PARALLEL */

/* ---------  Allocate 3D arrays to hold Gas based on size of grid --------- */

  if (pGrid->Nx1 > 1)
    Nx1T = pGrid->Nx1 + 2*nghost;
  else
    Nx1T = 1;

  if (pGrid->Nx2 > 1)
    Nx2T = pGrid->Nx2 + 2*nghost;
  else
    Nx2T = 1;

  if (pGrid->Nx3 > 1)
    Nx3T = pGrid->Nx3 + 2*nghost;
  else
    Nx3T = 1;

/* Build a 3D array of type Gas */

  pGrid->U = (Gas***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(Gas));
  if (pGrid->U == NULL) goto on_error;

/* Build 3D arrays to hold interface field */

#ifdef MHD
  pGrid->B1i = (Real***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(Real));
  if (pGrid->B1i == NULL) {
    free_3d_array((void***)pGrid->U);
    goto on_error;
  }

  pGrid->B2i = (Real***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(Real));
  if (pGrid->B2i == NULL) {
    free_3d_array((void***)pGrid->U);
    free_3d_array((void***)pGrid->B1i);
    goto on_error;
  }

  pGrid->B3i = (Real***)calloc_3d_array(Nx3T, Nx2T, Nx1T, sizeof(Real));
  if (pGrid->B3i == NULL) {
    free_3d_array((void***)pGrid->U);
    free_3d_array((void***)pGrid->B1i);
    free_3d_array((void***)pGrid->B2i);
    goto on_error;
  }
#endif /* MHD */

  return;

  on_error:
    ath_error("[init_grid_block]: Error allocating memory\n");
}
