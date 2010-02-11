#include "copyright.h"
/*==============================================================================
 * FILE: set_bvals_grav.c
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) for the
 *   gravitational potential on each edge of a Grid.  See comments at
 *   start of set_bvals_mhd.c for more details.
 * The only BC functions implemented here are for:
 *  1 = reflecting, 4 = periodic, and MPI boundaries
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   set_bvals_grav()      - calls appropriate functions to set ghost cells
 *   set_bvals_grav_init() - sets function pointers used by set_bvals_grav()
 *   set_bvals_grav_fun()  - enrolls a pointer to a user-defined BC function
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/* The functions in this file will only work with SELF_GRAVITY */
#ifdef SELF_GRAVITY

#ifdef MPI_PARALLEL
/* MPI send and receive buffer, size dynamically determined near end of
 * set_bvals_grav_init() based on number of zones in each grid */
static double *send_buf = NULL, *recv_buf = NULL;
#endif /* MPI_PARALLEL */

/* boundary condition function pointers, local to this file */
static VBCFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VBCFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VBCFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - apply reflecting BCs at boundary ???
 *   periodic_???() - apply periodic BCs at boundary ???
 *   send_???()     - MPI send of data at ??? boundary
 *   receive_???()  - MPI receive of data at ??? boundary
 *============================================================================*/

static void reflect_ix1(Grid *pG);
static void reflect_ox1(Grid *pG);
static void reflect_ix2(Grid *pG);
static void reflect_ox2(Grid *pG);
static void reflect_ix3(Grid *pG);
static void reflect_ox3(Grid *pG);

static void periodic_ix1(Grid *pG);
static void periodic_ox1(Grid *pG);
static void periodic_ix2(Grid *pG);
static void periodic_ox2(Grid *pG);
static void periodic_ix3(Grid *pG);
static void periodic_ox3(Grid *pG);

#ifdef MPI_PARALLEL
static void send_ix1(Grid *pG);
static void send_ox1(Grid *pG);
static void send_ix2(Grid *pG);
static void send_ox2(Grid *pG);
static void send_ix3(Grid *pG);
static void send_ox3(Grid *pG);

static void receive_ix1(Grid *pG, MPI_Request *prq);
static void receive_ox1(Grid *pG, MPI_Request *prq);
static void receive_ix2(Grid *pG, MPI_Request *prq);
static void receive_ox2(Grid *pG, MPI_Request *prq);
static void receive_ix3(Grid *pG, MPI_Request *prq);
static void receive_ox3(Grid *pG, MPI_Request *prq);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* set_bvals_grav: calls appropriate functions to set ghost zones.  The
 *   function pointers (*apply_???) are set during initialization by
 *   set_bvals_grav_init() to be either a user-defined function, or one of the
 *   functions corresponding to reflecting or periodic.  If the left-
 *   or right-Grid ID numbers are >= 1 (neighboring grids exist), then MPI calls
 *   are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void set_bvals_grav(Grid *pGrid, Domain *pDomain)
{
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr;
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
    cnt = nghost*cnt2*cnt3;

/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox1   (pGrid);       /* send R */
      receive_ix1(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix1   (pGrid);       /* send L */
      receive_ox1(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox1    (pGrid);       /* send R */
      (*apply_ix1)(pGrid);
      receive_ox1 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix1    (pGrid);       /* send L */
      (*apply_ox1)(pGrid);
      receive_ix1 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*apply_ix1)(pGrid);
      (*apply_ox1)(pGrid);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx2 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx1 > 1 ? pGrid->Nx1 + 2*nghost : 1;
    cnt3 = pGrid->Nx3 > 1 ? pGrid->Nx3 + 1 : 1;
    cnt = nghost*cnt1*cnt3;

/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox2   (pGrid);       /* send R */
      receive_ix2(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix2   (pGrid);       /* send L */
      receive_ox2(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox2    (pGrid);       /* send R */
      (*apply_ix2)(pGrid);
      receive_ox2 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix2    (pGrid);       /* send L */
      (*apply_ox2)(pGrid);
      receive_ix2 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*apply_ix2)(pGrid);
      (*apply_ox2)(pGrid);
    }

/* shearing sheet BCs; function defined in problem generator */
#ifdef SHEARING_BOX
    get_myGridIndex(pDomain, pGrid->my_id, &my_iproc, &my_jproc, &my_kproc);
    if (my_iproc == 0) {
      ShearingSheet_grav_ix1(pGrid, pDomain);
    }
    if (my_iproc == (pDomain->NGrid_x1-1)) {
      ShearingSheet_grav_ox1(pGrid, pDomain);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx3 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx1 > 1 ? pGrid->Nx1 + 2*nghost : 1;
    cnt2 = pGrid->Nx2 > 1 ? pGrid->Nx2 + 2*nghost : 1;
    cnt = nghost*cnt1*cnt2;

/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox3   (pGrid);       /* send R */
      receive_ix3(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix3   (pGrid);       /* send L */
      receive_ox3(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ox3    (pGrid);       /* send R */
      (*apply_ix3)(pGrid);
      receive_ox3 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);

      send_ix3    (pGrid);       /* send L */
      (*apply_ox3)(pGrid);
      receive_ix3 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*apply_ix3)(pGrid);
      (*apply_ox3)(pGrid);
    }

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_grav_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void set_bvals_grav_init(Grid *pG, Domain *pD)
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x1 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x1 = %d unknown\n",ibc_x1);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x1 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x1 = %d unknown\n",obc_x1);
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x2 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x2 = %d unknown\n",ibc_x2);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x2 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x2 = %d unknown\n",obc_x2);
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x3 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: ibc_x3 = %d unknown\n",ibc_x3);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x3 = 2 not defined\n");
	exit(EXIT_FAILURE);
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
	ath_perr(-1,"[set_bvals_grav_init]: obc_x3 = %d unknown\n",obc_x3);
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
    if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
      ath_error("[set_bvals_grav_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
      ath_error("[set_bvals_grav_init]: Failed to allocate receive buffer\n");
  }
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_grav_fun: sets function pointers for user-defined BCs in problem 
 */

void set_bvals_grav_fun(enum Direction dir, VBCFun_t prob_bc)
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
    ath_perr(-1,"[set_bvals_grav_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???
 *   periodic_???
 *   send_???
 *   receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1)
 */

static void reflect_ix1(Grid *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][is-i] = pGrid->Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1)
 */

static void reflect_ox1(Grid *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][ie+i] = pGrid->Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 */

static void reflect_ix2(Grid *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

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
        pGrid->Phi[k][js-j][i]    =  pGrid->Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1)
 */

static void reflect_ox2(Grid *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

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
        pGrid->Phi[k][je+j][i] = pGrid->Phi[k][je-(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 */

static void reflect_ix3(Grid *pGrid)
{
  int ks = pGrid->ks;
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

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->Phi[ks-k][j][i] = pGrid->Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1)
 */

static void reflect_ox3(Grid *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

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
        pGrid->Phi[ke+k][j][i] = pGrid->Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_ix1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][is-i] = pGrid->Phi[k][j][ie-(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (obc_x1=4)
 */

static void periodic_ox1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->Phi[k][j][ie+i] = pGrid->Phi[k][j][is+(i-1)];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2(Grid *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

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
        pGrid->Phi[k][js-j][i] = pGrid->Phi[k][je-(j-1)][i];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (obc_x2=4)
 */

static void periodic_ox2(Grid *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */

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
        pGrid->Phi[k][je+j][i] = pGrid->Phi[k][js+(j-1)][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (ibc_x3=4)
 */

static void periodic_ix3(Grid *pGrid)
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

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->Phi[ks-k][j][i] = pGrid->Phi[ke-(k-1)][j][i];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (obc_x3=4)
 */

static void periodic_ox3(Grid *pGrid)
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

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->Phi[ke+k][j][i] = pGrid->Phi[ks+(k-1)][j][i];
      }
    }
  }

  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~550 lines */

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x1 boundary -- send left
 */

static void send_ix1(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x1 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id,
		 boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x1 boundary -- send right
 */

static void send_ox1(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x1 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id,
		 boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x2 boundary -- send left
 */

static void send_ix2(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */

  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x2 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
		 boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x2 boundary -- send right
 */

static void send_ox2(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */

  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x2 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
		 boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x3 boundary -- send left
 */

static void send_ix3(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */

  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x3 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx3_id,
		  boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x3 boundary -- send right
 */

static void send_ox3(Grid *pG)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
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

/* Pack only Phi into send buffer */

  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1);
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        *(pd++) = pG->Phi[k][j][i];
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x3 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx3_id,
		  boundary_cells_tag, MPI_COMM_WORLD);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x1 boundary -- listen left
 */

static void receive_ix1(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x1 boundary -- listen right
 */

static void receive_ox1(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x2 boundary -- listen left
 */

static void receive_ix2(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x2 boundary -- listen right
 */

static void receive_ox2(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x3 boundary -- listen left
 */

static void receive_ix3(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x3 boundary -- listen right
 */

static void receive_ox3(Grid *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
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

  ierr = MPI_Wait(prq, &stat);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        pG->Phi[k][j][i] = *(pd++);
      }
    }
  }

  return;
}

#endif /* MPI_PARALLEL */

#endif /* SELF_GRAVITY */
