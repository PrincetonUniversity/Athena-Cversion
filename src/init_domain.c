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
 * Should be called directly before init_grid(), which initializes data on an
 * individual Grid block.  Needed even for single processor jobs to set
 * indices in Domain structure, and IDs of neighboring grids in Grid structure
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_domain()
 *   get_myblock_ijk()
 *============================================================================*/

#include <math.h>
#include <stdlib.h>
#include "athena.h"
#include "globals.h"
#include "copyright.h"
#include "defs.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
static int dom_decomp(const int Nx, const int Ny, const int Nz,
		      const int Np, int *pNGx, int *pNGy, int *pNGz);
#endif

/* Number of grids in each direction, read from inputfile in init_domain()*/
int NGrid_x1, NGrid_x2, NGrid_x3; 

/*----------------------------------------------------------------------------*/
/* init_domain:  */

void init_domain(Grid *pG, Domain *pD)
{
  int i,j,k,ib,jb,kb,nx1,nx2,nx3,id;
  div_t x1div, x2div, x3div;        /* A divisor with quot and rem members */

/* Initialize the min/max coordinates of this Domain, using parameters read
 * from <grid> block in input file.
 * Currently ixs,jxs,kxs hardwired to be zero */

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
  if(par_geti_def("parallel","auto",0)){
    if(dom_decomp(nx1, nx2, nx3, pG->nproc,
		  &NGrid_x1, &NGrid_x2, &NGrid_x3))
      ath_error("[init_domain]: Auto. Domain Decomposition Error\n");

    /* Store the domain decomposition in the par database */
    par_seti("parallel","NGrid_x1","%d",NGrid_x1,"x1-dir. domain decomposition");
    par_seti("parallel","NGrid_x2","%d",NGrid_x2,"x2-dir. domain decomposition");
    par_seti("parallel","NGrid_x3","%d",NGrid_x3,"x3-dir. domain decomposition");
  }
  else{
    NGrid_x1 = par_geti("parallel","NGrid_x1");
    NGrid_x2 = par_geti("parallel","NGrid_x2");
    NGrid_x3 = par_geti("parallel","NGrid_x3");
  }
#else
  NGrid_x1 = 1;
  NGrid_x2 = 1;
  NGrid_x3 = 1;
#endif

  if(NGrid_x1*NGrid_x2*NGrid_x3 != pG->nproc)
    ath_error("[init_domain]: There are %d Grids and %d processes\n",
	      NGrid_x1*NGrid_x2*NGrid_x3, pG->nproc);

/* test for dimensionality conflicts */

  if(NGrid_x1 > 1 && nx1 <= 1)
    ath_error("[init_domain]: NGrid_x = %d and Nx = %d\n",
	      NGrid_x1,nx1);

  if(NGrid_x2 > 1 && nx2 <= 1)
    ath_error("[init_domain]: NGrid_y = %d and Ny = %d\n",
	      NGrid_x2,nx2);

  if(NGrid_x3 > 1 && nx3 <= 1)
    ath_error("[init_domain]: NGrid_z = %d and Nz = %d\n",
	      NGrid_x3,nx3);

/* Build a 3D array of type Grid_Block */

  pD->grid_block = (Grid_Block***)calloc_3d_array(NGrid_x3, NGrid_x2, NGrid_x1,
     sizeof(Grid_Block));
  if(pD->grid_block == NULL)
    ath_error("[init_domain]: malloc returned a NULL pointer\n");

/* Divide the domain into blocks */

  x1div = div(nx1, NGrid_x1);
  x2div = div(nx2, NGrid_x2);
  x3div = div(nx3, NGrid_x3);

/* Assign each grid block in Domain to one of the processor IDs created during
 * MPI_Init() in main.c.  We rely on the IDs being numbered sequentially from 0
 * to nproc-1.  We use a regular grid, with neighboring IDs in the x1-direction.
 * For single-processor jobs, the grid_block array has only one member.  */

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
 * computational domain) of the grid block being updated on this processor
 * (ib,jb,kb)  */

  get_myblock_ijk(pD, pG->my_id, &ib, &jb, &kb);

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
/* get_myblock_ijk: searches grid_block[][][] array to find i,j,k components
 *   of block being updated on this processor.  */

void get_myblock_ijk(Domain *pD, const int my_id, int *pi, int *pj, int *pk)
{
  int i, j, k;

  for (k=0; k<NGrid_x3; k++){
    for (j=0; j<NGrid_x2; j++){
      for (i=0; i<NGrid_x1; i++){
	if (pD->grid_block[k][j][i].my_id == my_id) {
	  *pi = i;  *pj = j;  *pk = k;
	  return;
	}
      }
    }
  }

  ath_error("[get_myblock_ijk]: Can't find my_id=%i in the grid_block array\n",
    my_id);
}


/* ========================================================================== */


/* Note, the TOTAL amount of data communicated (summed over all
   processes and all INTERNAL boundaries) divided by 2*nghost where
   the 2 is for two messages per internal interface is computed and
   stored in the variable I, the minimum of which is I0. */
/* The assumption here is that in the x-direction we communicate only
   computational cells, while in the y-direction we comunicate the
   computational cells and x-direction ghost cells. */
/* Total x-communication is (rx - 1)*Ny*(2*nghost) */
/* Total y-communication is (ry - 1)*(Nx + rx*(2*nghost))*(2*nghost) */

static int dom_decomp_2d(const int Nx, const int Ny,
			 const int Np, int *pNGx, int *pNGy){

  int rx, ry, I, rxs, rx_min, rx_max;
  int rx0=1, ry0=1, I0=0, init=1;
  div_t dv;

  /* Compute the ideal decomposition, truncated to an integer, which
     minimizes the amount of communication. */
  rxs = (int)sqrt((double)(Nx*Np)/(double)(Ny > 2*nghost ? Ny - 2*nghost : 1));

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1 and rx <= Nx */
  if(Ny < Np){ /* Require rx >= 1 and ry <= Ny */
    dv = div(Np, Ny);
    /* rx_min = the smallest integer >= Np/Ny */
    rx_min = dv.quot + (dv.rem > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  /* Constrain rxs to fall in this domain */
  rxs = rxs > rx_min ? rxs : rx_min;
  rxs = rxs < rx_max ? rxs : rx_max;

  /* Search down for a factor of Np */
  for(rx=rxs; rx>=rx_min; rx--){
    dv = div(Np,rx);
    if(dv.rem == 0){
      rx0 = rx;
      ry0 = dv.quot;
      I0 = (rx0 - 1)*Ny + (ry0 - 1)*(Nx + 2*nghost*rx0);
      init = 0;
      break;
    }
  }

  /* Search up for a factor of Np */
  for(rx=rxs+1; rx<=rx_max; rx++){
    dv = div(Np,rx);
    if(dv.rem == 0){
      ry = dv.quot;
      I = (rx - 1)*Ny + (ry - 1)*(Nx + 2*nghost*rx);

      if(init || I < I0){
	rx0 = rx;
	ry0 = ry;
	I0  = I;
	init = 0;
      }
      break;
    }
  }

  if(init) return 1; /* Error locating a solution */

  /* printf("Minimum messaging decomposition has: rx = %d, ry = %d, I = %d\n",
     rx0, ry0, I0); */

  *pNGx = rx0;
  *pNGy = ry0;

  return 0;
}


/* ========================================================================== */


static int dom_decomp_3d(const int Nx, const int Ny, const int Nz,
			 const int Np, int *pNGx, int *pNGy, int *pNGz){

  div_t dv;
  int rx_min, rx_max, rx, ry, rz, I;
  int rx0=1, ry0=1, rz0=1, I0=0, init=1;
  int err, t, Npt;

  /* Constrain the decomposition */
  rx_max = Nx < Np ? Nx : Np; /* Require ry >= 1, rz >= 1 and rx <= Nx */
  /* Compute a global minimum constraint on rx. */
  t = (Ny < Np ? Ny : Np)*(Nz < Np ? Nz : Np); /* t = Max(ry)*Max(rz) */
  if(t < Np){ /* Require rx >= 1, ry <= Ny and rz <= Nz */
    dv = div(Np, t);
    /* rx_min = the smallest integer >= Np/t */
    rx_min = dv.quot + (dv.rem > 0 ? 1 : 0);
  }
  else rx_min = 1;

  /* printf("rx_min = %d, rx_max = %d\n",rx_min, rx_max); */

  for(rx = rx_min; rx <= rx_max; rx++){
    dv = div(Np, rx);
    if(dv.rem == 0){
      Npt = dv.quot; /* Np for transverse (y,z) decomposition */

      err = dom_decomp_2d(Ny, Nz, Npt, &ry, &rz);
      if(err == 0){
	/* Now compute the amount of messaging */
	I = (rx - 1)*Ny*Nz + (ry - 1)*(Nx + 2*nghost*rx)*Nz
	  + (rz - 1)*(Nx + 2*nghost*rx)*(Ny + 2*nghost*ry);

	if(I < 0){ /* Integer Overflow */
	  /* printf("[3d new] I = %d\n",I); */
	  continue;
	}

	if(init || I < I0){
	  rx0 = rx;
	  ry0 = ry;
	  rz0 = rz;
	  I0  = I;
	  init = 0;
	  /* printf("I(rx = %d, ry = %d, rz = %d) = %d\n",rx,ry,rz,I); */
	}
      }
    }
  }

  if(init) return 1; /* Error locating a solution */

  *pNGx = rx0;
  *pNGy = ry0;
  *pNGz = rz0;

  return 0;
}


/* ========================================================================== */


static int dom_decomp(const int Nx, const int Ny, const int Nz,
		      const int Np, int *pNGx, int *pNGy, int *pNGz){

  if(Nx > 1 && Ny == 1 && Nz == 1){ /* 1-D */
    if(Np > Nx) return 1; /* Too many procs. */
    *pNGx = Np;
    *pNGy = 1;
    *pNGz = 1;
    return 0;
  }
  else if(Nx > 1 && Ny > 1 && Nz == 1){ /* 2-D */
    *pNGy = 1;
    return dom_decomp_2d(Nx, Nz, Np, pNGx, pNGz);
  }
  else if(Nx > 1 && Ny > 1 && Nz > 1){ /* 3-D */
    return dom_decomp_3d(Nx, Ny, Nz, Np, pNGx, pNGy, pNGz);
  }

  return 1; /* Error - particular case not expected */
}


/* ========================================================================== */
