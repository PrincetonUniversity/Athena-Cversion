#include "copyright.h"
/*==============================================================================
 * FILE: set_bvals.c
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) on each edge
 *   of a Grid.  Each edge of a Grid represents either:
 *    (1) the physical boundary of computational Domain; in which case BCs are
 *        specified by an integer flag input by user (or by user-defined BC
 *        function in the problem file)
 *    (2) the boundary between Grids resulting from decomposition of a larger
 *        computational domain using MPI; in which case BCs are handled
 *        by MPI calls
 *    (3) an internal boundary between fine/coarse grid levels in a nested grid;
 *        in which case BCs require prolongation and restriction operators (and
 *        possibly some combination of (1) and (2) as well!)
 *   This file contains functions called in the main loop that can handle each
 *   of these cases.  The naming convention for BCs is:
 *       ibc_x1 = Inner Boundary Condition for x1
 *       obc_x1 = Outer Boundary Condition for x1
 *   similarly for ibc_x2; obc_x2; ibc_x3; obc_x3
 *
 * For case (1) -- PHYSICAL BOUNDARIES
 *   The values of the integer flags are:
 *       1 = reflecting; 2 = outflow; 4 = periodic
 *   Following ZEUS conventions, 3 would be flow-in (ghost zones held at
 *   pre-determined fixed values), however in Athena instead we use pointers to
 *   user-defined BC functions for flow-in.
 *
 * For case (2) -- MPI BOUNDARIES
 *   We do the parallel synchronization by having every grid:
 *     1) Pack and send data to the grid on right  [send_ox1()]
 *     2) Listen to the left, unpack and set data  [receive_ix1()]
 *     3) Pack and send data to the grid on left   [send_ix1()]
 *     4) Listen to the right, unpack and set data [receive_ox1()]
 *   If the grid is at the edge of the Domain, we set BCs as in case (1) or (3).
 *   Currently the code uses NON-BLOCKING sends (MPI_Isend) and BLOCKING
 *   receives (MPI_Recv).  Some optimization could be achieved by interleaving
 *   non-blocking sends (MPI_Isend) and computations.
 *
 * For case (3) -- INTERNAL GRID LEVEL BOUNDARIES
 *
 * The type of BC is unchanged throughout a calculation.  Thus, during setup
 * we determine the BC type, and set a pointer to the appropriate BC function
 * using set_bvals_init().  MPI calls are used if the grid ID number to the
 * left or right is >= 0.
 * 
 * With SELF-GRAVITY: BCs for Phi have to be set independently of the Gas.  The
 *   argument list of the BC functions contain an integer flag which determines
 *   whether BCs for Phi (flag=1) or the Gas structure (flag=0) are set.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   set_bvals()      - calls appropriate functions to set ghost cells
 *   set_bvals_init() - sets function pointers used by set_bvals() based on flag
 *   set_bvals_fun()  - enrolls a pointer to a user-defined BC function
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL

/* MPI send and receive buffer, size dynamically determined near end of
 * set_bvals_init() based on number of zones in each grid */
static double *send_buf = NULL, *recv_buf = NULL;

/* Maximim number of variables passed in any one MPI message = variables in the
 * Gas structure, plus 3 extra for the interface magnetic fields.
 *
 * With self-gravity, BCs for Phi have to be set before the source term
 * correction in the main loop.  So Phi is always passed on its own, not
 * with the Gas variables.   */
#ifdef MHD
#define NVAR_SHARE (NVAR + 3)
#else
#define NVAR_SHARE NVAR
#endif

#endif /* MPI_PARALLEL */

/* boundary condition function pointers */
static VBCFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VBCFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VBCFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - apply reflecting BCs at boundary ???
 *   outflow_???()  - apply outflow BCs at boundary ???
 *   periodic_???() - apply periodic BCs at boundary ???
 *   send_???()     - MPI send of data at ??? boundary
 *   receive_???()  - MPI receive of data at ??? boundary
 *============================================================================*/

static void reflect_ix1(Grid *pG, int var_flag);
static void reflect_ox1(Grid *pG, int var_flag);
static void reflect_ix2(Grid *pG, int var_flag);
static void reflect_ox2(Grid *pG, int var_flag);
static void reflect_ix3(Grid *pG, int var_flag);
static void reflect_ox3(Grid *pG, int var_flag);

static void outflow_ix1(Grid *pG, int var_flag);
static void outflow_ox1(Grid *pG, int var_flag);
static void outflow_ix2(Grid *pG, int var_flag);
static void outflow_ox2(Grid *pG, int var_flag);
static void outflow_ix3(Grid *pG, int var_flag);
static void outflow_ox3(Grid *pG, int var_flag);

static void periodic_ix1(Grid *pG, int var_flag);
static void periodic_ox1(Grid *pG, int var_flag);
static void periodic_ix2(Grid *pG, int var_flag);
static void periodic_ox2(Grid *pG, int var_flag);
static void periodic_ix3(Grid *pG, int var_flag);
static void periodic_ox3(Grid *pG, int var_flag);

#ifdef MPI_PARALLEL
static void send_ix1(Grid *pG, int var_flag);
static void send_ox1(Grid *pG, int var_flag);
static void send_ix2(Grid *pG, int var_flag);
static void send_ox2(Grid *pG, int var_flag);
static void send_ix3(Grid *pG, int var_flag);
static void send_ox3(Grid *pG, int var_flag);

static void receive_ix1(Grid *pG, int var_flag, MPI_Request *prq);
static void receive_ox1(Grid *pG, int var_flag, MPI_Request *prq);
static void receive_ix2(Grid *pG, int var_flag, MPI_Request *prq);
static void receive_ox2(Grid *pG, int var_flag, MPI_Request *prq);
static void receive_ix3(Grid *pG, int var_flag, MPI_Request *prq);
static void receive_ox3(Grid *pG, int var_flag, MPI_Request *prq);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* set_bvals:  calls appropriate functions to set ghost zones.  The function
 *   pointers (*apply_???) are set during initialization by set_bvals_init()
 *   to be either a user-defined function, or one of the functions corresponding
 *   to reflecting, periodic, or outflow.  If the left- or right-Grid ID numbers
 *   are >= 1 (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void set_bvals(Grid *pGrid, Domain *pDomain, int var_flag)
{
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, err;
  MPI_Request rq;
#endif /* MPI_PARALLEL */
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx1 > 1){

#ifdef MPI_PARALLEL
    cnt2 = pGrid->Nx2 > 1 ? pGrid->Nx2 + 1 : 1;
    cnt3 = pGrid->Nx3 > 1 ? pGrid->Nx3 + 1 : 1;
    if (var_flag == 0) {
      cnt = nghost*cnt2*cnt3*NVAR_SHARE;
    } else {
      cnt = nghost*cnt2*cnt3;
    }

/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox1   (pGrid, var_flag);       /* send R */
      receive_ix1(pGrid, var_flag, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix1   (pGrid, var_flag);       /* send L */
      receive_ox1(pGrid, var_flag, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox1    (pGrid, var_flag);       /* send R */
      (*apply_ix1)(pGrid, var_flag);
      receive_ox1 (pGrid, var_flag, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix1    (pGrid, var_flag);       /* send L */
      (*apply_ox1)(pGrid, var_flag);
      receive_ix1 (pGrid, var_flag, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*apply_ix1)(pGrid, var_flag);
      (*apply_ox1)(pGrid, var_flag);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx2 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx1 > 1 ? pGrid->Nx1 + 2*nghost : 1;
    cnt3 = pGrid->Nx3 > 1 ? pGrid->Nx3 + 1 : 1;
    if (var_flag == 0) {
      cnt = nghost*cnt1*cnt3*NVAR_SHARE;
    } else {
      cnt = nghost*cnt1*cnt3;
    }

/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox2   (pGrid, var_flag);       /* send R */
      receive_ix2(pGrid, var_flag, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix2   (pGrid, var_flag);       /* send L */
      receive_ox2(pGrid, var_flag, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox2    (pGrid, var_flag);       /* send R */
      (*apply_ix2)(pGrid, var_flag);
      receive_ox2 (pGrid, var_flag, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix2    (pGrid, var_flag);       /* send L */
      (*apply_ox2)(pGrid, var_flag);
      receive_ix2 (pGrid, var_flag, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*apply_ix2)(pGrid, var_flag);
      (*apply_ox2)(pGrid, var_flag);
    }

/* shearing sheet BCs; function defined in problem generator */
#ifdef SHEARING_BOX
    get_myGridIndex(pDomain, pGrid->my_id, &my_iproc, &my_jproc, &my_kproc);
    if (my_iproc == 0) {
      ShearingSheet_ix1(pGrid, pDomain);
    }
    if (my_iproc == (pDomain->NGrid_x1-1)) {
      ShearingSheet_ox1(pGrid, pDomain);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx3 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx1 > 1 ? pGrid->Nx1 + 2*nghost : 1;
    cnt2 = pGrid->Nx2 > 1 ? pGrid->Nx2 + 2*nghost : 1;
    if (var_flag == 0) {
      cnt = nghost*cnt1*cnt2*NVAR_SHARE;
    } else {
      cnt = nghost*cnt1*cnt2;
    }

/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox3   (pGrid, var_flag);       /* send R */
      receive_ix3(pGrid, var_flag, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix3   (pGrid, var_flag);       /* send L */
      receive_ox3(pGrid, var_flag, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ox3    (pGrid, var_flag);       /* send R */
      (*apply_ix3)(pGrid, var_flag);
      receive_ox3 (pGrid, var_flag, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      send_ix3    (pGrid, var_flag);       /* send L */
      (*apply_ox3)(pGrid, var_flag);
      receive_ix3 (pGrid, var_flag, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*apply_ix3)(pGrid, var_flag);
      (*apply_ox3)(pGrid, var_flag);
    }

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void set_bvals_init(Grid *pG, Domain *pD)
{
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
#ifdef MPI_PARALLEL
  int i,j,k,ib,jb,kb;
  int my_id = pG->my_id;
  int x1cnt, x2cnt, x3cnt; /* Number of Gas passed in x1-, x2-, x3-dir. */
  int nx1t, nx2t, nx3t, size;
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction */

  if(pG->Nx1 > 1) {
    if(apply_ix1 == NULL){

      ibc_x1 = par_geti("grid","ibc_x1");
      switch(ibc_x1){

      case 1: /* Reflecting */
	apply_ix1 = reflect_ix1;
	break;

      case 2: /* Outflow */
	apply_ix1 = outflow_ix1;
	break;

      case 4: /* Periodic */
	apply_ix1 = periodic_ix1;
#ifdef MPI_PARALLEL
	if(pG->lx1_id < 0 && pD->NGrid_x1 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx1_id = pD->GridArray[kb][jb][pD->NGrid_x1-1].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: ibc_x1 = %d unknown\n",ibc_x1);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox1 == NULL){

      obc_x1 = par_geti("grid","obc_x1");
      switch(obc_x1){

      case 1: /* Reflecting */
	apply_ox1 = reflect_ox1;
	break;

      case 2: /* Outflow */
	apply_ox1 = outflow_ox1;
	break;

      case 4: /* Periodic */
	apply_ox1 = periodic_ox1;
#ifdef MPI_PARALLEL
	if(pG->rx1_id < 0 && pD->NGrid_x1 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx1_id = pD->GridArray[kb][jb][0].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: obc_x1 = %d unknown\n",obc_x1);
	exit(EXIT_FAILURE);
      }

    }
  }

/* Set function pointers for physical boundaries in x2-direction */

  if(pG->Nx2 > 1) {
    if(apply_ix2 == NULL){

      ibc_x2 = par_geti("grid","ibc_x2");
      switch(ibc_x2){

      case 1: /* Reflecting */
	apply_ix2 = reflect_ix2;
	break;

      case 2: /* Outflow */
	apply_ix2 = outflow_ix2;
	break;

      case 4: /* Periodic */
	apply_ix2 = periodic_ix2;
#ifdef MPI_PARALLEL
	if(pG->lx2_id < 0 && pD->NGrid_x2 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx2_id = pD->GridArray[kb][pD->NGrid_x2-1][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: ibc_x2 = %d unknown\n",ibc_x2);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox2 == NULL){

      obc_x2 = par_geti("grid","obc_x2");
      switch(obc_x2){

      case 1: /* Reflecting */
	apply_ox2 = reflect_ox2;
	break;

      case 2: /* Outflow */
	apply_ox2 = outflow_ox2;
	break;

      case 4: /* Periodic */
	apply_ox2 = periodic_ox2;
#ifdef MPI_PARALLEL
	if(pG->rx2_id < 0 && pD->NGrid_x2 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx2_id = pD->GridArray[kb][0][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: obc_x2 = %d unknown\n",obc_x2);
	exit(EXIT_FAILURE);
      }

    }
  }

/* Set function pointers for physical boundaries in x3-direction */

  if(pG->Nx3 > 1) {
    if(apply_ix3 == NULL){

      ibc_x3 = par_geti("grid","ibc_x3");
      switch(ibc_x3){

      case 1: /* Reflecting */
	apply_ix3 = reflect_ix3;
	break;

      case 2: /* Outflow */
	apply_ix3 = outflow_ix3;
	break;

      case 4: /* Periodic */
	apply_ix3 = periodic_ix3;
#ifdef MPI_PARALLEL
	if(pG->lx3_id < 0 && pD->NGrid_x3 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx3_id = pD->GridArray[pD->NGrid_x3-1][jb][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: ibc_x3 = %d unknown\n",ibc_x3);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox3 == NULL){

      obc_x3 = par_geti("grid","obc_x3");
      switch(obc_x3){

      case 1: /* Reflecting */
	apply_ox3 = reflect_ox3;
	break;

      case 2: /* Outflow */
	apply_ox3 = outflow_ox3;
	break;

      case 4: /* Periodic */
	apply_ox3 = periodic_ox3;
#ifdef MPI_PARALLEL
	if(pG->rx3_id < 0 && pD->NGrid_x3 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx3_id = pD->GridArray[0][jb][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_init]: obc_x3 = %d unknown\n",obc_x3);
	exit(EXIT_FAILURE);
      }

    }
  }

/* allcoate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  x1cnt = x2cnt = x3cnt = 0;

  for (k=0; k<(pD->NGrid_x3); k++){
    for (j=0; j<(pD->NGrid_x2); j++){
      for (i=0; i<(pD->NGrid_x1); i++){
	if(pD->NGrid_x1 > 1){
	  nx2t = pD->GridArray[k][j][i].jge - pD->GridArray[k][j][i].jgs + 1;
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GridArray[k][j][i].kge - pD->GridArray[k][j][i].kgs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
	}

	if(pD->NGrid_x2 > 1){
	  nx1t = pD->GridArray[k][j][i].ige - pD->GridArray[k][j][i].igs + 1;
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx3t = pD->GridArray[k][j][i].kge - pD->GridArray[k][j][i].kgs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
	}


	if(pD->NGrid_x3 > 1){
	  nx1t = pD->GridArray[k][j][i].ige - pD->GridArray[k][j][i].igs + 1;
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx2t = pD->GridArray[k][j][i].jge - pD->GridArray[k][j][i].jgs + 1;
	  if(nx2t > 1) nx2t += 2*nghost;

	  x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
	}
      }
    }
  }

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= nghost; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate receive buffer\n");
  }
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_fun:  sets function pointers for user-defined BCs in problem file
 */

void set_bvals_fun(enum Direction dir, VBCFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    apply_ix1 = prob_bc;
    break;
  case right_x1:
    apply_ox1 = prob_bc;
    break;
  case left_x2:
    apply_ix2 = prob_bc;
    break;
  case right_x2:
    apply_ox2 = prob_bc;
    break;
  case left_x3:
    apply_ix3 = prob_bc;
    break;
  case right_x3:
    apply_ox3 = prob_bc;
    break;
  default:
    ath_perr(-1,"[set_bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???
 *   outflow_???
 *   periodic_???
 *   send_???
 *   receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ix1(Grid *pGrid, int var_flag)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) return;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ox1(Grid *pGrid, int var_flag)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) return;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie-(i-2)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ix2(Grid *pGrid, int var_flag)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i]    =  pGrid->U[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].M2 = -pGrid->U[k][js-j][i].M2; /* reflect 2-mom. */
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js+(j-1)][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ox2(Grid *pGrid, int var_flag)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = -pGrid->U[k][je+j][i].M2; /* reflect 2-mom. */
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ix3(Grid *pGrid, int var_flag)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1)
 * Phi set by multipole expansion external to this function
 */

static void reflect_ox3(Grid *pGrid, int var_flag)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke-(k-2)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditionss, Inner x1 boundary (ibc_x1=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ix1(Grid *pGrid, int var_flag)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) return;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x1 boundary (obc_x1=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ox1(Grid *pGrid, int var_flag)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) return;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x2 boundary (ibc_x2=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ix2(Grid *pGrid, int var_flag)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][js][i];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js][i];
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x2 boundary (obc_x2=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ox2(Grid *pGrid, int var_flag)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][je][i];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je][i];
      }
    }
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x3 boundary (ibc_x3=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ix3(Grid *pGrid, int var_flag)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U  [ks-k][j][i] = pGrid->U  [ks][j][i];
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x3 boundary (obc_x3=2)
 * Phi set by multipole expansion external to this function
 */

static void outflow_ox3(Grid *pGrid, int var_flag)
{
  int ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (var_flag == 1) return;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_ix1(Grid *pGrid, int var_flag)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=1; i<=nghost; i++) {
            pGrid->Phi[k][j][is-i] = pGrid->Phi[k][j][ie-(i-1)];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][ie-(i-1)];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (obc_x1=4)
 */

static void periodic_ox1(Grid *pGrid, int var_flag)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=1; i<=nghost; i++) {
          pGrid->Phi[k][j][ie+i] = pGrid->Phi[k][j][is+(i-1)];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][is+(i-1)];
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2(Grid *pGrid, int var_flag)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=ks; k<=ke; k++) {
      for (j=1; j<=nghost; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->Phi[k][js-j][i] = pGrid->Phi[k][je-(j-1)][i];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][js-j][i] = pGrid->U[k][je-(j-1)][i];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][js-j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js-j][i] = pGrid->B2i[k][je-(j-1)][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (obc_x2=4)
 */

static void periodic_ox2(Grid *pGrid, int var_flag)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=ks; k<=ke; k++) {
      for (j=1; j<=nghost; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->Phi[k][je+j][i] = pGrid->Phi[k][js+(j-1)][i];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i] = pGrid->U[k][js+(j-1)][i];
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][js+(j-1)][i];
      }
    }
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][js+(j-1)][i];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (ibc_x3=4)
 */

static void periodic_ix3(Grid *pGrid, int var_flag)
{
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->Phi[ks-k][j][i] = pGrid->Phi[ke-(k-1)][j][i];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ke-(k-1)][j][i];
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (obc_x3=4)
 */

static void periodic_ox3(Grid *pGrid, int var_flag)
{
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (var_flag == 1) {
#ifdef SELF_GRAVITY
    for (k=1; k<=nghost; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->Phi[ke+k][j][i] = pGrid->Phi[ks+(k-1)][j][i];
        }
      }
    }
#endif /* SELF_GRAVITY */
    return;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ks+(k-1)][j][i];
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~760 lines */

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x1 boundary -- send left
 */

static void send_ix1(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  il = pG->is;
  iu = pG->is + nghost - 1;

  if(pG->Nx2 > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 1 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */
#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif

/* send contents of buffer to the neighboring grid on L-x1 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id,
		 boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x1 boundary -- send right
 */

static void send_ox1(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  il = pG->ie - nghost + 1;
  iu = pG->ie;

  if(pG->Nx2 > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 1 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */
#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif

/* send contents of buffer to the neighboring grid on R-x1 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id,
		 boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x2 boundary -- send left
 */

static void send_ix2(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->js;
  ju = pG->js + nghost - 1;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 2 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

/* send contents of buffer to the neighboring grid on L-x2 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
		 boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix2]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x2 boundary -- send right
 */

static void send_ox2(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->je - nghost + 1;
  ju = pG->je;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 2 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

/* send contents of buffer to the neighboring grid on R-x2 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
		 boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox2]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x3 boundary -- send left
 */

static void send_ix3(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx2 > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ks;
  ku = pG->ks + nghost - 1;

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 3 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

/* send contents of buffer to the neighboring grid on L-x3 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx3_id,
		  boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix3]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x3 boundary -- send right
 */

static void send_ox3(Grid *pG, int var_flag)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  double *pd = send_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx2 > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ke - nghost + 1;
  ku = pG->ke;

/* Pack data in Gas structure into send buffer */

  if (var_flag == 0) {
    /* Following expression gives same cnt as in Step 3 in set_bvals()  */
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          *(pd++) = pq->d;
          *(pd++) = pq->M1;
          *(pd++) = pq->M2;
          *(pd++) = pq->M3;
#ifdef MHD
          *(pd++) = pq->B1c;
          *(pd++) = pq->B2c;
          *(pd++) = pq->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pq->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pq->s[n];
#endif
        }
      }
    }
  }

/* Pack only Phi into send buffer */

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pG->Phi[k][j][i];
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

/* send contents of buffer to the neighboring grid on R-x3 */

  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx3_id,
		  boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox3]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x1 boundary -- listen left
 */

static void receive_ix1(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  il = pG->is - nghost;
  iu = pG->is - 1;

  if(pG->Nx2 > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the left grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix1]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x1 boundary -- listen right
 */

static void receive_ox1(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  il = pG->ie + 1;
  iu = pG->ie + nghost;

  if(pG->Nx2 > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the right grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox1]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GTRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x2 boundary -- listen left
 */

static void receive_ix2(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->js - nghost;
  ju = pG->js - 1;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the left grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix2]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x2 boundary -- listen right
 */

static void receive_ox2(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->je + 1;
  ju = pG->je + nghost;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the right grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox2]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x3 boundary -- listen left
 */

static void receive_ix3(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx2 > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ks - nghost;
  ku = pG->ks - 1;

/* Wait to receive the input data from the left grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix3]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x3 boundary -- listen right
 */

static void receive_ox3(Grid *pG, int var_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
#if (NSCALARS > 0)
  int n;
#endif
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx2 > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ke + 1;
  ku = pG->ke + nghost;

/* Wait to receive the input data from the right grid */

  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox3]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */

  if (var_flag == 0) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          /* Get a pointer to the Gas cell */
          pq = &(pG->U[k][j][i]);

          pq->d = *(pd++);
          pq->M1 = *(pd++);
          pq->M2 = *(pd++);
          pq->M3 = *(pd++);
#ifdef MHD
          pq->B1c = *(pd++);
          pq->B2c = *(pd++);
          pq->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pq->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pq->s[n] = *(pd++);
#endif
        }
      }
    }
  }

#ifdef SELF_GRAVITY
  if (var_flag == 1) {
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }
#endif /* SELF_GRAVITY */

  return;
}

#endif /* MPI_PARALLEL */
