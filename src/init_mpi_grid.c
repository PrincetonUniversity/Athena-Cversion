#include "copyright.h"
/*==============================================================================
 * FILE: init_mpi_grid.c
 *
 * PURPOSE: Functions to decompose a grid level into an arbitrary nmuber of
 *   MPI blocks.  The grid can be decomposed in any direction.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   domain_partition()
 *   domain_destruct()
 *   domain_ijk()
 *============================================================================*/

#include <stdlib.h>
#include "athena.h"
#include "copyright.h"
#include "defs.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL

/* Min and Max bounding coordinates of the Complete Computational Grid. */
int cg_ixs, cg_jxs, cg_kxs;
int cg_ixe, cg_jxe, cg_kxe;

/* Number of grids in each direction, read from inputfile in domain_partition()*/
int NGrid_x1, NGrid_x2, NGrid_x3; 

/* 3D array of type Domain with dimensions = number of grids in 1,2,3 */
Domain ***grid_domain=NULL;


/*----------------------------------------------------------------------------*/
/* domain_partition:  */

void domain_partition(Grid *pG)
{
  int i, j, k;
  div_t x1div, x2div, x3div;        /* A divisor with quot and rem members */
  int id;

/* Initialize the min / max coordinate of the Complete comp. Grid */
  cg_ixs = pG->is + pG->idisp;
  cg_jxs = pG->js + pG->jdisp;
  cg_kxs = pG->ks + pG->kdisp;

  cg_ixe = pG->ie + pG->idisp;
  cg_jxe = pG->je + pG->jdisp;
  cg_kxe = pG->ke + pG->kdisp;

/* Get the number of grids in each direction */
  NGrid_x1 = par_geti("parallel","NGrid_x1");
  NGrid_x2 = par_geti("parallel","NGrid_x2");
  NGrid_x3 = par_geti("parallel","NGrid_x3");

  if(NGrid_x1*NGrid_x2*NGrid_x3 != pG->nproc)
    ath_error("[domain_partition]: There are %d Grids and %d processes\n",
	      NGrid_x1*NGrid_x2*NGrid_x3, pG->nproc);

/* test for dimensionality conflicts */
  if(NGrid_x1 > 1 && pG->Nx1 <= 1)
    ath_error("[domain_partition]: NGrid_x = %d and Nx = %d\n",
	      NGrid_x1,pG->Nx1);

  if(NGrid_x2 > 1 && pG->Nx2 <= 1)
    ath_error("[domain_partition]: NGrid_y = %d and Ny = %d\n",
	      NGrid_x2,pG->Nx2);

  if(NGrid_x3 > 1 && pG->Nx3 <= 1)
    ath_error("[domain_partition]: NGrid_z = %d and Nz = %d\n",
	      NGrid_x3,pG->Nx3);

/* Allocate the grid_domain[][][] array */
  grid_domain =
    (Domain***)calloc_3d_array(NGrid_x3, NGrid_x2, NGrid_x1, sizeof(Domain));
  if(grid_domain == NULL)
    ath_error("[domain_partition]: malloc returned a NULL pointer\n");

/* initialize the my_id field of the domain array */
  id = 0;
  for(k=0; k<NGrid_x3; k++){
    for(j=0; j<NGrid_x2; j++){
      for(i=0; i<NGrid_x1; i++){
	grid_domain[k][j][i].my_id = id++;
	/* printf("grid_domain[%d][%d][%d].my_id = %d\n",
	   k,j,i,grid_domain[k][j][i].my_id); */
      }
    }
  }

  x1div = div(pG->Nx1, NGrid_x1);
  x2div = div(pG->Nx2, NGrid_x2);
  x3div = div(pG->Nx3, NGrid_x3);

/* Initialize the domain of the first (0,0,0) Grid Domain */
  grid_domain[0][0][0].ixs = cg_ixs;
  grid_domain[0][0][0].jxs = cg_jxs;
  grid_domain[0][0][0].kxs = cg_kxs;

  grid_domain[0][0][0].ixe = cg_ixs + x1div.quot - 1;
  grid_domain[0][0][0].jxe = cg_jxs + x2div.quot - 1;
  grid_domain[0][0][0].kxe = cg_kxs + x3div.quot - 1;

/* If the domain is not evenly divisible put the extra cells on the
 * inner Grid domains */
  if(x1div.rem > 0){
    grid_domain[0][0][0].ixe++;
    x1div.rem--;
  }

  if(x2div.rem > 0){
    grid_domain[0][0][0].jxe++;
    x2div.rem--;
  }

  if(x3div.rem > 0){
    grid_domain[0][0][0].kxe++;
    x3div.rem--;
  }

/* Initialize the 1D domains along each direction independently */
/* Along the x1-direction */
  for(i=1; i<NGrid_x1; i++){
    grid_domain[0][0][i].ixs = grid_domain[0][0][i-1].ixe + 1;
    if(x1div.rem > 0){
      grid_domain[0][0][i].ixe =
	grid_domain[0][0][i].ixs + x1div.quot;
      x1div.rem--;
    }
    else{
      grid_domain[0][0][i].ixe =
	grid_domain[0][0][i].ixs + x1div.quot - 1;
    }
  }

/* Along the x2-direction */
  for(j=1; j<NGrid_x2; j++){
    grid_domain[0][j][0].jxs = grid_domain[0][j-1][0].jxe + 1;
    if(x2div.rem > 0){
      grid_domain[0][j][0].jxe =
	grid_domain[0][j][0].jxs + x2div.quot;
      x2div.rem--;
    }
    else{
      grid_domain[0][j][0].jxe =
	grid_domain[0][j][0].jxs + x2div.quot - 1;
    }
  }

/* Along the x3-direction */
  for(k=1; k<NGrid_x3; k++){
    grid_domain[k][0][0].kxs = grid_domain[k-1][0][0].kxe + 1;
    if(x3div.rem > 0){
      grid_domain[k][0][0].kxe =
	grid_domain[k][0][0].kxs + x3div.quot;
      x3div.rem--;
    }
    else{
      grid_domain[k][0][0].kxe =
	grid_domain[k][0][0].kxs + x3div.quot - 1;
    }
  }

/* Finally fill in the rest */
  for(k=0; k<NGrid_x3; k++){
    for(j=0; j<NGrid_x2; j++){
      for(i=0; i<NGrid_x1; i++){
/* The x-coordinates */
	if(j != 0 || k != 0){
	  grid_domain[k][j][i].ixs = grid_domain[0][0][i].ixs;
	  grid_domain[k][j][i].ixe = grid_domain[0][0][i].ixe;
	}
/* The y-coordinates */
	if(i != 0 || k != 0){
	  grid_domain[k][j][i].jxs = grid_domain[0][j][0].jxs;
	  grid_domain[k][j][i].jxe = grid_domain[0][j][0].jxe;
	}
/* The z-coordinates */
	if(i != 0 || j != 0){
	  grid_domain[k][j][i].kxs = grid_domain[k][0][0].kxs;
	  grid_domain[k][j][i].kxe = grid_domain[k][0][0].kxe;
	}
      }
    }
  }

/* Reset the Grid domain to match the domain of the sub-Grid just calculated. */
  domain_ijk(pG->my_id, &i, &j, &k);

/* Initialize Nx1, is, ie, idisp */
  pG->Nx1 = grid_domain[k][j][i].ixe - grid_domain[k][j][i].ixs + 1;

  if(pG->Nx1 > 1) {
    pG->is = nghost;
    pG->ie = pG->Nx1 + nghost - 1;
  }
  else
    pG->is = pG->ie = 0;

  pG->idisp = grid_domain[k][j][i].ixs - pG->is;

/* Initialize Nx2, js, je, jdisp */
  pG->Nx2 = grid_domain[k][j][i].jxe - grid_domain[k][j][i].jxs + 1;

  if(pG->Nx2 > 1) {
    pG->js = nghost;
    pG->je = pG->Nx2 + nghost - 1;
  }
  else
    pG->js = pG->je = 0;

  pG->jdisp = grid_domain[k][j][i].jxs - pG->js;

/* Initialize Nx3, ks, ke, kdisp */
  pG->Nx3 = grid_domain[k][j][i].kxe - grid_domain[k][j][i].kxs + 1;

  if(pG->Nx3 > 1) {
    pG->ks = nghost;
    pG->ke = pG->Nx3 + nghost - 1;
  }
  else
    pG->ks = pG->ke = 0;

  pG->kdisp = grid_domain[k][j][i].kxs - pG->ks;

/* Left-x1 */
  if(i > 0) pG->lx1_id = grid_domain[k][j][i-1].my_id;
  else pG->lx1_id = -1;

/* Right-x1 */
  if(i < (NGrid_x1 - 1)) pG->rx1_id = grid_domain[k][j][i+1].my_id;
  else pG->rx1_id = -1;

/* Left-x2 */
  if(j > 0) pG->lx2_id = grid_domain[k][j-1][i].my_id;
  else pG->lx2_id = -1;

/* Right-x2 */
  if(j < (NGrid_x2 - 1)) pG->rx2_id = grid_domain[k][j+1][i].my_id;
  else pG->rx2_id = -1;

/* Left-x3 */
  if(k > 0) pG->lx3_id = grid_domain[k-1][j][i].my_id;
  else pG->lx3_id = -1;

/* Right-x3 */
  if(k < (NGrid_x3 - 1)) pG->rx3_id = grid_domain[k+1][j][i].my_id;
  else pG->rx3_id = -1;

  return;
}

/*----------------------------------------------------------------------------*/
/* domain_destruct:  Free the grid_domain[][][] array */

void domain_destruct(void)
{
  if(grid_domain != NULL){
    free_3d_array((void***)grid_domain);
    grid_domain = NULL;
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* domain_ijk:
 * This is not optimal by any stretch, but I also don't expect to use
 * it more than a few times at the start of a simulation.   */

void domain_ijk(const int my_id, int *pi, int *pj, int *pk)
{
  int i, j, k;

  for(k=0; k<NGrid_x3; k++){
    for(j=0; j<NGrid_x2; j++){
      for(i=0; i<NGrid_x1; i++){
	if(grid_domain[k][j][i].my_id == my_id){
	  *pi = i;  *pj = j;  *pk = k;
	  return;
	}
      }
    }
  }

  ath_error("[domain_ijk]: Can't find (i,j,k) in the grid_domain array");
}

#endif
