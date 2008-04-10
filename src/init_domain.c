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
 *   get_myGridIndex()
 *============================================================================*/

#include <math.h>
#include <stdlib.h>
#include "athena.h"
#include "globals.h"
#include "copyright.h"
#include "defs.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   dom_decomp()    - calls auto domain decomposition functions 
 *   dom_decomp_2d() - finds optimum domain decomposition in 2D 
 *   dom_decomp_3d() - finds optimum domain decomposition in 3D 
 *============================================================================*/
#ifdef MPI_PARALLEL
static int dom_decomp(const int Nx, const int Ny, const int Nz,
		      const int Np, int *pNGx, int *pNGy, int *pNGz);
static int dom_decomp_2d(const int Nx, const int Ny,
                         const int Np, int *pNGx, int *pNGy);
static int dom_decomp_3d(const int Nx, const int Ny, const int Nz,
                         const int Np, int *pNGx, int *pNGy, int *pNGz);
#endif

/*----------------------------------------------------------------------------*/
/* init_domain:  */

void init_domain(Grid *pG, Domain *pD)
{
  int i,j,k,ib,jb,kb,id,nproc=0;
  div_t x1div, x2div, x3div;        /* A divisor with quot and rem members */

/* Initialize the min/max coordinates of this Domain, using parameters read
 * from <grid> block in input file. */

  pD->Nx1 = par_geti("grid","Nx1");
  if (pD->Nx1 > 1) {
    pD->ids = 0;
    pD->ide = pD->Nx1 - 1;
  }
  else if (pD->Nx1 == 1) {
    pD->ids = pD->ide = 0;
  }
  else {
    ath_error("[init_domain]: Nx1 = %d must be >= 1\n",pD->Nx1);
  }

  pD->Nx2 = par_geti("grid","Nx2");
  if (pD->Nx2 > 1) {
    pD->jds = 0;
    pD->jde = pD->Nx2 - 1;
  }
  else if (pD->Nx2 == 1) {
    pD->jds = pD->jde = 0;
  }
  else {
    ath_error("[init_domain]: Nx2 = %d must be >= 1\n",pD->Nx2);
  }

  pD->Nx3 = par_geti("grid","Nx3");
  if (pD->Nx3 > 1) {
    pD->kds = 0;
    pD->kde = pD->Nx3 - 1;
  }
  else if (pD->Nx3 == 1) {
    pD->kds = pD->kde = 0;
  }
  else {
    ath_error("[init_domain]: Nx3 = %d must be >= 1\n",pD->Nx3);
  }

/* Get the number of grids in each direction */

#ifdef MPI_PARALLEL
/* Get the number of processes */
  if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&nproc))
    ath_error("[init_domain]: Error on calling MPI_Comm_size\n");

/* Auto decompose the Domain into Grids */
  if(par_geti_def("parallel","auto",0)){
    if(dom_decomp(pD->Nx1,pD->Nx2,pD->Nx3,nproc,
                 &(pD->NGrid_x1),&(pD->NGrid_x2),&(pD->NGrid_x3)))
        ath_error("[init_domain]: Auto. Domain Decomposition Error\n");

    /* Store the domain decomposition in the par database */
    par_seti("parallel","NGrid_x1","%d",pD->NGrid_x1,"x1 domain decomposition");
    par_seti("parallel","NGrid_x2","%d",pD->NGrid_x2,"x2 domain decomposition");
    par_seti("parallel","NGrid_x3","%d",pD->NGrid_x3,"x3 domain decomposition");
  }

/* Else read number of Grids from input file */
  else{
    pD->NGrid_x1 = par_geti("parallel","NGrid_x1");
    pD->NGrid_x2 = par_geti("parallel","NGrid_x2");
    pD->NGrid_x3 = par_geti("parallel","NGrid_x3");
  }

  if((pD->NGrid_x1)*(pD->NGrid_x2)*(pD->NGrid_x3) != nproc)
    ath_error("[init_domain]: There are %d Grids and %d processes\n",
    (pD->NGrid_x1)*(pD->NGrid_x2)*(pD->NGrid_x3), nproc);

#else
  pD->NGrid_x1 = 1;
  pD->NGrid_x2 = 1;
  pD->NGrid_x3 = 1;
#endif /* MPI_PARALLEL */

/* test for dimensionality conflicts */

  if(pD->NGrid_x1 > 1 && pD->Nx1 <= 1)
    ath_error("[init_domain]: NGrid_x = %d and Nx = %d\n",pD->NGrid_x1,pD->Nx1);

  if(pD->NGrid_x2 > 1 && pD->Nx2 <= 1)
    ath_error("[init_domain]: NGrid_y = %d and Ny = %d\n",pD->NGrid_x2,pD->Nx2);

  if(pD->NGrid_x3 > 1 && pD->Nx3 <= 1)
    ath_error("[init_domain]: NGrid_z = %d and Nz = %d\n",pD->NGrid_x3,pD->Nx3);

/* Build a 3D array of type Grid_Indices */

  pD->GridArray = (Grid_Indices***)calloc_3d_array(pD->NGrid_x3, pD->NGrid_x2, 
     pD->NGrid_x1, sizeof(Grid_Indices));
  if(pD->GridArray == NULL)
    ath_error("[init_domain]: malloc returned a NULL pointer\n");

/* Divide the domain into blocks */

  x1div = div(pD->Nx1, pD->NGrid_x1);
  x2div = div(pD->Nx2, pD->NGrid_x2);
  x3div = div(pD->Nx3, pD->NGrid_x3);

/* Assign each Grid in Domain to one of the processor IDs created during
 * MPI_Init() in main.c.  We rely on the IDs being numbered sequentially from 0
 * to nproc-1.  We use a regular grid, with neighboring IDs in the x1-direction.
 * For single-processor jobs, the GridArray has only one member.  */

  id = 0;
  for(k=0; k<(pD->NGrid_x3); k++){
    for(j=0; j<(pD->NGrid_x2); j++){
      for(i=0; i<(pD->NGrid_x1); i++){
	pD->GridArray[k][j][i].id = id++;
      }
    }
  }

/* Initialize the min/max of the coordinate in the first (0,0,0) Grid */

  pD->GridArray[0][0][0].igs = pD->ids;
  pD->GridArray[0][0][0].jgs = pD->jds;
  pD->GridArray[0][0][0].kgs = pD->kds;

  pD->GridArray[0][0][0].ige = pD->ids + x1div.quot - 1;
  pD->GridArray[0][0][0].jge = pD->jds + x2div.quot - 1;
  pD->GridArray[0][0][0].kge = pD->kds + x3div.quot - 1;

/* If the Domain is not evenly divisible put the extra cells on the inner Grid*/

  if(x1div.rem > 0){
    pD->GridArray[0][0][0].ige++;
    x1div.rem--;
  }

  if(x2div.rem > 0){
    pD->GridArray[0][0][0].jge++;
    x2div.rem--;
  }

  if(x3div.rem > 0){
    pD->GridArray[0][0][0].kge++;
    x3div.rem--;
  }

/* Initialize the GridArray indices along each direction independently */
/* Along the x1-direction */

  for (i=1; i<(pD->NGrid_x1); i++){
    pD->GridArray[0][0][i].igs = pD->GridArray[0][0][i-1].ige + 1;
    if(x1div.rem > 0){
      pD->GridArray[0][0][i].ige = pD->GridArray[0][0][i].igs + x1div.quot;
      x1div.rem--;
    } else{
      pD->GridArray[0][0][i].ige = pD->GridArray[0][0][i].igs + x1div.quot - 1;
    }
  }

/* Along the x2-direction */

  for (j=1; j<(pD->NGrid_x2); j++){
    pD->GridArray[0][j][0].jgs = pD->GridArray[0][j-1][0].jge + 1;
    if(x2div.rem > 0){
      pD->GridArray[0][j][0].jge = pD->GridArray[0][j][0].jgs + x2div.quot;
      x2div.rem--;
    } else{
      pD->GridArray[0][j][0].jge = pD->GridArray[0][j][0].jgs + x2div.quot - 1;
    }
  }

/* Along the x3-direction */

  for (k=1; k<(pD->NGrid_x3); k++){
    pD->GridArray[k][0][0].kgs = pD->GridArray[k-1][0][0].kge + 1;
    if(x3div.rem > 0){
      pD->GridArray[k][0][0].kge = pD->GridArray[k][0][0].kgs + x3div.quot;
      x3div.rem--;
    } else{
      pD->GridArray[k][0][0].kge = pD->GridArray[k][0][0].kgs + x3div.quot - 1;
    }
  }

/* Finally fill in the rest */

  for (k=0; k<(pD->NGrid_x3); k++){
    for (j=0; j<(pD->NGrid_x2); j++){
      for (i=0; i<(pD->NGrid_x1); i++){
/* The x-coordinates */
	if(j != 0 || k != 0){
	  pD->GridArray[k][j][i].igs = pD->GridArray[0][0][i].igs;
	  pD->GridArray[k][j][i].ige = pD->GridArray[0][0][i].ige;
	}
/* The y-coordinates */
	if(i != 0 || k != 0){
	  pD->GridArray[k][j][i].jgs = pD->GridArray[0][j][0].jgs;
	  pD->GridArray[k][j][i].jge = pD->GridArray[0][j][0].jge;
	}
/* The z-coordinates */
	if(i != 0 || j != 0){
	  pD->GridArray[k][j][i].kgs = pD->GridArray[k][0][0].kgs;
	  pD->GridArray[k][j][i].kge = pD->GridArray[k][0][0].kge;
	}
      }
    }
  }

/* Now get the i,j,k coordinates (in the 3D array of Grid_Indices that tile the
 * computational domain) of the Grid being updated on this processor (ib,jb,kb)
 */

  get_myGridIndex(pD, pG->my_id, &ib, &jb, &kb);

/* Get IDs of neighboring Grids.  If edge of a Grid is at a physical boundary,
 * ID is set to -1
 */

/* Left-x1 */
  if(ib > 0) pG->lx1_id = pD->GridArray[kb][jb][ib-1].id;
  else pG->lx1_id = -1;

/* Right-x1 */
  if(ib < (pD->NGrid_x1 - 1)) pG->rx1_id = pD->GridArray[kb][jb][ib+1].id;
  else pG->rx1_id = -1;

/* Left-x2 */
  if(jb > 0) pG->lx2_id = pD->GridArray[kb][jb-1][ib].id;
  else pG->lx2_id = -1;

/* Right-x2 */
  if(jb < (pD->NGrid_x2 - 1)) pG->rx2_id = pD->GridArray[kb][jb+1][ib].id;
  else pG->rx2_id = -1;

/* Left-x3 */
  if(kb > 0) pG->lx3_id = pD->GridArray[kb-1][jb][ib].id;
  else pG->lx3_id = -1;

/* Right-x3 */
  if(kb < (pD->NGrid_x3 - 1)) pG->rx3_id = pD->GridArray[kb+1][jb][ib].id;
  else pG->rx3_id = -1;

  return;
}

/*----------------------------------------------------------------------------*/
/* get_myGridIndex: searches GridArray[][][] array to find i,j,k components
 *   of block being updated on this processor.  */

void get_myGridIndex(Domain *pD, const int my_id, int *pi, int *pj, int *pk)
{
  int i, j, k;

  for (k=0; k<(pD->NGrid_x3); k++){
    for (j=0; j<(pD->NGrid_x2); j++){
      for (i=0; i<(pD->NGrid_x1); i++){
	if (pD->GridArray[k][j][i].id == my_id) {
	  *pi = i;  *pj = j;  *pk = k;
	  return;
	}
      }
    }
  }

  ath_error("[get_myGridIndex]: Can't find id=%i in the GridArray\n", my_id);
}

#ifdef MPI_PARALLEL
/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* dom_decomp: calls apropriate 2D or 3D auto decomposition routines
 *   Functions written by T.A.G., added May 2007
 */

static int dom_decomp(const int Nx, const int Ny, const int Nz,
                      const int Np, int *pNGx, int *pNGy, int *pNGz)
{
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

/*----------------------------------------------------------------------------*/
/* dom_decomp_2d: optimizes domain decomposition in 2D.  The TOTAL amount of
 *   data communicated (summed over all processes and all INTERNAL boundaries)
 *   divided by 2*nghost (where the 2 is for two messages per internal
 *   interface) is computed and stored in the variable I, the minimum of which
 *   is I0.  This assumes that in the x-direction we communicate only
 *   computational cells, while in the y-direction we comunicate the
 *   computational cells and x-direction ghost cells.  Then
 *     Total x-communication is (rx - 1)*Ny*(2*nghost) 
 *     Total y-communication is (ry - 1)*(Nx + rx*(2*nghost))*(2*nghost)
 */

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

/*----------------------------------------------------------------------------*/
/* dom_decomp_3d: optimizes domain decomposition in 3D.  See the comments for
 *   dom_decomp_2d() for more about the algorithm
 */

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

#endif /* MPI_PARALLEL */
