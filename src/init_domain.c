#include "copyright.h"
/*==============================================================================
 * FILE: init_domain.c
 *
 * PURPOSE: Initialize a single Domain.  A Domain is the entire region over
 *   which we integrate the solution.  For a single processor job, the Domain
 *   and the Grid are identical.  For MPI parallel jobs, each Grid is one
 *   tile in the entire Domain; and init_domain() handles the distribution of
 *   the Domain over multiple processors.  For a nested grid, each level is a
 *   Domain, each of which may contain multiple Grids (MPI blocks).
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_domain()
 *   domain_destruct()
 *   domain_ijk()
 *============================================================================*/

#include <stdlib.h>
#include "athena.h"
#include "copyright.h"
#include "defs.h"
#include "prototypes.h"


/* Number of grids in each direction, read from inputfile in domain_partition()*/
int NGrid_x1, NGrid_x2, NGrid_x3; 


/*----------------------------------------------------------------------------*/
/* init_domain:  */

void init_domain(Grid *pG, Domain *pD)
{
  int i,j,k,ib,jb,kb;
  div_t x1div, x2div, x3div;        /* A divisor with quot and rem members */
  int id;

/* Initialize the min/max coordinates of this Domain, using parameters read
 * from <grid> block in input file.
 * Currently ixs,jxs,kxs hardwired to be zero; needs change for nested grids */

  nx1 = par_geti("grid","Nx1");
  if (nx1 > 1) {
    pD->ixs = 0;
    pD->ixe = nx1 - 1;
  }
  else if (nx1 == 1) {
    pD->ixs = pD->ixe = 0;
  }
  else {
    ath_error("[init_domain]: Nx1 = %d must be >= 1\n",nx1);
  }

  nx2 = par_geti("grid","Nx2");
  if (nx2 > 1) {
    pD->jxs = 0;
    pD->jxe = nx2 - 1;
  }
  else if (nx2 == 1) {
    pD->jxs = pD->jxe = 0;
  }
  else {
    ath_error("[init_domain]: Nx2 = %d must be >= 1\n",nx2);
  }

  nx3 = par_geti("grid","Nx3");
  if (nx3 > 1) {
    pD->kxs = 0;
    pD->kxe = nx3 - 1;
  }
  else if (nx3 == 1) {
    pD->kxs = pD->kxe = 0;
  }
  else {
    ath_error("[init_domain]: Nx3 = %d must be >= 1\n",nx3);
  }

/* Get the number of grids in each direction */

#ifdef MPI_PARALLEL
  NGrid_x1 = par_geti("parallel","NGrid_x1");
  NGrid_x2 = par_geti("parallel","NGrid_x2");
  NGrid_x3 = par_geti("parallel","NGrid_x3");

  if(NGrid_x1*NGrid_x2*NGrid_x3 != pG->nproc)
    ath_error("[domain_partition]: There are %d Grids and %d processes\n",
	      NGrid_x1*NGrid_x2*NGrid_x3, pG->nproc);
#endif

/* test for dimensionality conflicts */

  if(NGrid_x1 > 1 && nx1 <= 1)
    ath_error("[domain_partition]: NGrid_x = %d and Nx = %d\n",
	      NGrid_x1,nx1);

  if(NGrid_x2 > 1 && nx2 <= 1)
    ath_error("[domain_partition]: NGrid_y = %d and Ny = %d\n",
	      NGrid_x2,nx2);

  if(NGrid_x3 > 1 && nx3 <= 1)
    ath_error("[domain_partition]: NGrid_z = %d and Nz = %d\n",
	      NGrid_x3,nx3);

/* Build a 3D array of type Grid_Block */

  pD->grid_block = (Grid_Block***)calloc_3d_array(NGrid_x3, NGrid_x2, NGrid_x1,
     sizeof(Grid_Block));
  if(pD->grid_block == NULL)
    ath_error("[domain_partition]: malloc returned a NULL pointer\n");

/* Divide the domain into blocks */

  x1div = div(nx1, NGrid_x1);
  x2div = div(nx2, NGrid_x2);
  x3div = div(nx3, NGrid_x3);

/* Assign processor IDs to each grid block in domain.  We use a regular grid,
 * with neighboring IDs in the x1-direction.  This could be changed if other
 * layouts are more efficient on certain machines
 */

  id = 0;
  for(k=0; k<NGrid_x3; k++){
    for(j=0; j<NGrid_x2; j++){
      for(i=0; i<NGrid_x1; i++){
	pD->grid_block[k][j][i].my_id = id++;
      }
    }
  }

/* Initialize the min/max of the coordinate in the first (0,0,0) Grid Block */

  pD->grid_block[0][0][0].ixs = pD->ixs;
  pD->grid_block[0][0][0].jxs = pD->jxs;
  pD->grid_block[0][0][0].kxs = pD->kxs;

  pD->grid_block[0][0][0].ixe = pD->ixs + x1div.quot - 1;
  pD->grid_block[0][0][0].jxe = pD->jxs + x2div.quot - 1;
  pD->grid_block[0][0][0].kxe = pD->kxs + x3div.quot - 1;

/* If the domain is not evenly divisible put the extra cells on the inner 
 * Grid blocks */

  if(x1div.rem > 0){
    pD->grid_block[0][0][0].ixe++;
    x1div.rem--;
  }

  if(x2div.rem > 0){
    pD->grid_block[0][0][0].jxe++;
    x2div.rem--;
  }

  if(x3div.rem > 0){
    pD->grid_block[0][0][0].kxe++;
    x3div.rem--;
  }

/* Initialize the 1D domains along each direction independently */
/* Along the x1-direction */

  for (i=1; i<NGrid_x1; i++){
    pD->grid_block[0][0][i].ixs = pD->grid_block[0][0][i-1].ixe + 1;
    if(x1div.rem > 0){
      pD->grid_block[0][0][i].ixe =
	pD->grid_block[0][0][i].ixs + x1div.quot;
      x1div.rem--;
    } else{
      pD->grid_block[0][0][i].ixe =
	pD->grid_block[0][0][i].ixs + x1div.quot - 1;
    }
  }

/* Along the x2-direction */

  for (j=1; j<NGrid_x2; j++){
    pD->grid_block[0][j][0].jxs = pD->grid_block[0][j-1][0].jxe + 1;
    if(x2div.rem > 0){
      pD->grid_block[0][j][0].jxe =
	pD->grid_block[0][j][0].jxs + x2div.quot;
      x2div.rem--;
    } else{
      pD->grid_block[0][j][0].jxe =
	pD->grid_block[0][j][0].jxs + x2div.quot - 1;
    }
  }

/* Along the x3-direction */

  for (k=1; k<NGrid_x3; k++){
    pD->grid_block[k][0][0].kxs = pD->grid_block[k-1][0][0].kxe + 1;
    if(x3div.rem > 0){
      pD->grid_block[k][0][0].kxe =
	pD->grid_block[k][0][0].kxs + x3div.quot;
      x3div.rem--;
    } else{
      pD->grid_block[k][0][0].kxe =
	pD->grid_block[k][0][0].kxs + x3div.quot - 1;
    }
  }

/* Finally fill in the rest */

  for (k=0; k<NGrid_x3; k++){
    for (j=0; j<NGrid_x2; j++){
      for (i=0; i<NGrid_x1; i++){
/* The x-coordinates */
	if(j != 0 || k != 0){
	  pD->grid_block[k][j][i].ixs = pD->grid_block[0][0][i].ixs;
	  pD->grid_block[k][j][i].ixe = pD->grid_block[0][0][i].ixe;
	}
/* The y-coordinates */
	if(i != 0 || k != 0){
	  pD->grid_block[k][j][i].jxs = pD->grid_block[0][j][0].jxs;
	  pD->grid_block[k][j][i].jxe = pD->grid_block[0][j][0].jxe;
	}
/* The z-coordinates */
	if(i != 0 || j != 0){
	  pD->grid_block[k][j][i].kxs = pD->grid_block[k][0][0].kxs;
	  pD->grid_block[k][j][i].kxe = pD->grid_block[k][0][0].kxe;
	}
      }
    }
  }

/* Now get the i,j,k coordinates (in the 3D array of grid blocks that tile the
 * compuational domain) of the grid block being updated on this processor
 * (ib,jb,kb)  */

  get_block_ijk(pG->my_id, &ib, &jb, &kb);

/* Get IDs of neighboring grid blocks.  If edge of a grid block is at a
 * physical boundary, ID is set to -1
 */

/* Left-x1 */
  if(ib > 0) pG->lx1_id = pD->grid_block[kb][jb][ib-1].my_id;
  else pG->lx1_id = -1;

/* Right-x1 */
  if(ib < (NGrid_x1 - 1)) pG->rx1_id = pD->grid_block[kb][jb][ib+1].my_id;
  else pG->rx1_id = -1;

/* Left-x2 */
  if(jb > 0) pG->lx2_id = pD->grid_block[kb][jb-1][ib].my_id;
  else pG->lx2_id = -1;

/* Right-x2 */
  if(jb < (NGrid_x2 - 1)) pG->rx2_id = pD->grid_block[kb][jb+1][ib].my_id;
  else pG->rx2_id = -1;

/* Left-x3 */
  if(kb > 0) pG->lx3_id = pD->grid_block[kb-1][jb][ib].my_id;
  else pG->lx3_id = -1;

/* Right-x3 */
  if(kb < (NGrid_x3 - 1)) pG->rx3_id = pD->grid_block[kb+1][jb][ib].my_id;
  else pG->rx3_id = -1;

  return;
}

/*----------------------------------------------------------------------------*/
/* domain_destruct:  Free the grid_domain[][][] array */

void domain_destruct(void)
{
  if(pD->grid_block != NULL){
    free_3d_array((void***)pD->grid_block);
    pD->grid_block = NULL;
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
	if(pD->grid_block[k][j][i].my_id == my_id){
	  *pi = i;  *pj = j;  *pk = k;
	  return;
	}
      }
    }
  }

  ath_error("[domain_ijk]: Can't find (i,j,k) in the grid_domain array");
}

#endif
