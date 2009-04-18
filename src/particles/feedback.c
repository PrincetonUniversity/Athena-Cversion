#include "../copyright.h"
/*=============================================================================
FILE: feedback.c
PURPOSE: Exchange particle feedback between boundary cells. The procedure is
  opposite to setting boundary conditions. Extra feedback forces are exterted
  to cells outside the grid boundary, by feedback exchange, these forces are
  properly added to cells inside the grid boundary, according to various
  boundary conditions.
    Currently, for shearing box simulation in 3D, FARGO must be enabled.

CONTAINS PUBLIC FUNCTIONS:
  void exchange_feedback(Grid *pG, Domain *pD);
  void exchange_feedback_init(Grid *pG, Domain *pD);
  void exchange_feedback_fun(enum Direction dir, VBCFun_t prob_bc);
  void exchange_feedback_destruct(Grid *pG, Domain *pD);

History:
  Written by Xuening Bai, Apr. 2009
==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"


#ifdef FEEDBACK         /* endif at the end of the file */

/* number of quantities to be exchanged for each cell */
#define NVAR_F 3
/* number of boundary cells to be exchanged */
#define NGF    2

#ifdef MPI_PARALLEL

/* MPI send and receive buffer, size dynamically determined near end of
 * set_bvals_init() based on number of zones in each grid */
static double *send_buf = NULL, *recv_buf = NULL;

#endif /* MPI_PARALLEL */

static int my_iproc, my_jproc, my_kproc;	/* processor indices in the computational domain */
static int il,iu, jl,ju, kl,ku;		/* grid index limit for feedback exchange */
static Real x1min,x1max,x2min,x2max,x3min,x3max;/* min and max coordinate limits of the computational domain */
static Real Lx1, Lx2, Lx3;			/* domain size in x1, x2, x3 direction */

/* boundary condition function pointers. local to this function  */
static VBCFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VBCFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VBCFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*====================== PROTOTYPE OF PRIVATE FUNCTIONS ======================*/
/*----------------------------------------------------------------------------*/

static void reflect_ix1_feedback(Grid *pG);
static void reflect_ox1_feedback(Grid *pG);
static void reflect_ix2_feedback(Grid *pG);
static void reflect_ox2_feedback(Grid *pG);
static void reflect_ix3_feedback(Grid *pG);
static void reflect_ox3_feedback(Grid *pG);

static void outflow_feedback(Grid *pG);

static void periodic_ix1_feedback(Grid *pG);
static void periodic_ox1_feedback(Grid *pG);
static void periodic_ix2_feedback(Grid *pG);
static void periodic_ox2_feedback(Grid *pG);
static void periodic_ix3_feedback(Grid *pG);
static void periodic_ox3_feedback(Grid *pG);

#ifdef MPI_PARALLEL
static void send_ix1_feedback(Grid *pG);
static void send_ox1_feedback(Grid *pG);
static void send_ix2_feedback(Grid *pG);
static void send_ox2_feedback(Grid *pG);
static void send_ix3_feedback(Grid *pG);
static void send_ox3_feedback(Grid *pG);

static void recv_ix1_feedback(Grid *pG, MPI_Request *prq);
static void recv_ox1_feedback(Grid *pG, MPI_Request *prq);
static void recv_ix2_feedback(Grid *pG, MPI_Request *prq);
static void recv_ox2_feedback(Grid *pG, MPI_Request *prq);
static void recv_ix3_feedback(Grid *pG, MPI_Request *prq);
static void recv_ox3_feedback(Grid *pG, MPI_Request *prq);
#endif /* MPI_PARALLEL */

#ifdef SHEARING_BOX
static void shearingbox_ix1_feedback(Grid *pG, Domain *pD);
static void shearingbox_ox1_feedback(Grid *pG, Domain *pD);
#endif /* SHEARING_BOX */


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* exchange_feedback: calls appropriate functions to set ghost zones.  The function
 *   pointers (*apply_???) are set during initialization by set_bvals_particle_init()
 *   to be either a user-defined function, or one of the functions corresponding
 *   to reflecting, periodic, or outflow.  If the left- or right-Grid ID numbers
 *   are >= 1 (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void exchange_feedback(Grid *pG, Domain *pD)
{
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, err;
  MPI_Request rq;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Feedback exchange in x3-direction */

  if (pG->Nx3 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pG->Nx1 > 1 ? pG->Nx1 + 1 : 1;
    cnt2 = pG->Nx2 > 1 ? pG->Nx2 + 1 : 1;
    cnt = NGF*cnt1*cnt2*NVAR_F;

/* MPI blocks to both left and right */
    if (pG->rx3_id >= 0 && pG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx3_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox3_feedback(pG);       /* send R */
      recv_ix3_feedback(pG, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx3_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix3_feedback(pG);       /* send L */
      recv_ox3_feedback(pG, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx3_id >= 0 && pG->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx3_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox3_feedback(pG);       /* send R */
      (*apply_ix3)     (pG);
      recv_ox3_feedback(pG, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx3_id < 0 && pG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx3_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix3_feedback(pG);       /* send L */
      (*apply_ox3)     (pG);
      recv_ix3_feedback(pG, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx3_id < 0 && pG->lx3_id < 0) {
      (*apply_ix3)(pG);
      (*apply_ox3)(pG);
    }

  }

/*--- Step 2. ------------------------------------------------------------------
 * Feedback exchange in x2-direction */

  if (pG->Nx2 > 1){

#ifdef MPI_PARALLEL
    cnt1 = pG->Nx1 > 1 ? pG->Nx1 + 1 : 1;
    cnt3 = pG->Nx3 > 1 ? pG->Nx3 + 2*NGF : 1;
    cnt = NGF*cnt1*cnt3*NVAR_F;

/* MPI blocks to both left and right */
    if (pG->rx2_id >= 0 && pG->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox2_feedback(pG);       /* send R */
      recv_ix2_feedback(pG, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix2_feedback(pG);       /* send L */
      recv_ox2_feedback(pG, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx2_id >= 0 && pG->lx2_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox2_feedback(pG);       /* send R */
      (*apply_ix2)     (pG);
      recv_ox2_feedback(pG, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx2_id < 0 && pG->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix2_feedback(pG);       /* send L */
      (*apply_ox2)     (pG);
      recv_ix2_feedback(pG, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx2_id < 0 && pG->lx2_id < 0) {
      (*apply_ix2)(pG);
      (*apply_ox2)(pG);
    }

  }

/*--- Step 2.cont.  -------------------------------------------------------------
 * shift the feedback in x1 direction */
#ifdef SHEARING_BOX
#ifndef FARGO
  if (pG->Nx3 > 1) /* 3D shearing box */
  {
    if (my_iproc == 0)
      shearingbox_ix1_feedback(pG, pD);

    if (my_iproc == pD->NGrid_x1-1)
      shearingbox_ox1_feedback(pG, pD);
  }
#endif
#endif

/*--- Step 3. ------------------------------------------------------------------
 * Feedback exchange in x1-direction */

  if (pG->Nx1 > 1){

#ifdef MPI_PARALLEL
    cnt2 = pG->Nx2 > 1 ? pG->Nx2 + 2*NGF : 1;
    cnt3 = pG->Nx3 > 1 ? pG->Nx3 + 2*NGF : 1;
    cnt = NGF*cnt2*cnt3*NVAR_F;

/* MPI blocks to both left and right */
    if (pG->rx1_id >= 0 && pG->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx1_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox1_feedback(pG);       /* send R */
      recv_ix1_feedback(pG, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx1_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix1_feedback(pG);       /* send L */
      recv_ox1_feedback(pG, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx1_id >= 0 && pG->lx1_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx1_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ox1_feedback(pG);       /* send R */
      (*apply_ix1)     (pG);
      recv_ox1_feedback(pG, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx1_id < 0 && pG->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx1_id, boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[exchange_feedback]: MPI_Irecv error = %d\n",err);

      send_ix1_feedback(pG);       /* send L */
      (*apply_ox1)     (pG);
      recv_ix1_feedback(pG, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx1_id < 0 && pG->lx1_id < 0) {
      (*apply_ix1)(pG);
      (*apply_ox1)(pG);
    } 

  }

  return;

}

/*----------------------------------------------------------------------------*/
/* exchange_feedback_init:  sets function pointers for feedback exchange during
 *   initialization, allocates memory for send/receive buffers with MPI
 */
void exchange_feedback_init(Grid *pG, Domain *pD)
{
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
#ifdef MPI_PARALLEL
  int i,j,k;
  int x1cnt, x2cnt, x3cnt; /* Number of Gas passed in x1-, x2-, x3-dir. */
  int nx1t, nx2t, nx3t, size;
#endif /* MPI_PARALLEL */

/* set left and right grid indices */
  if (pG->Nx1 > 1) {
    il = pG->is - NGF;		iu = pG->ie + NGF;
  } else {
    il = pG->is;		iu = pG->ie;
  }

  if (pG->Nx2 > 1) {
    jl = pG->js - NGF;		ju = pG->je + NGF;
  } else {
    jl = pG->js;		ju = pG->je;
  }

  if (pG->Nx3 > 1) {
    kl = pG->ks - NGF;		ku = pG->ke + NGF;
  } else {
    kl = pG->ks;		ku = pG->ke;
  }

/* calculate distances of the computational domain and shear velocity */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx1 = x1max - x1min;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Lx2 = x2max - x2min;

  x3min = par_getd("grid","x3min");
  x3max = par_getd("grid","x3max");
  Lx3 = x3max - x3min;

  get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

/* Set function pointers for physical boundaries in x1-direction */

  if(pG->Nx1 > 1) {
    if(apply_ix1 == NULL){

      ibc_x1 = par_geti("grid","ibc_x1");
      switch(ibc_x1){

      case 1: /* Reflecting */
	apply_ix1 = reflect_ix1_feedback;
	break;
      case 5: /* Reflecting */
	apply_ix1 = reflect_ix1_feedback;
	break;

      case 2: /* Outflow */
	apply_ix1 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ix1 = periodic_ix1_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: ibc_x1 = %d unknown\n",ibc_x1);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox1 == NULL){

      obc_x1 = par_geti("grid","obc_x1");
      switch(obc_x1){

      case 1: /* Reflecting */
	apply_ox1 = reflect_ox1_feedback;
	break;
      case 5: /* Reflecting */
	apply_ox1 = reflect_ox1_feedback;
	break;

      case 2: /* Outflow */
	apply_ox1 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ox1 = periodic_ox1_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: obc_x1 = %d unknown\n",obc_x1);
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
	apply_ix2 = reflect_ix2_feedback;
	break;
      case 5: /* Reflecting */
	apply_ix2 = reflect_ix2_feedback;
	break;

      case 2: /* Outflow */
	apply_ix2 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ix2 = periodic_ix2_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: ibc_x2 = %d unknown\n",ibc_x2);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox2 == NULL){

      obc_x2 = par_geti("grid","obc_x2");
      switch(obc_x2){

      case 1: /* Reflecting */
	apply_ox2 = reflect_ox2_feedback;
	break;
      case 5: /* Reflecting */
	apply_ox2 = reflect_ox2_feedback;
	break;

      case 2: /* Outflow */
	apply_ox2 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ox2 = periodic_ox2_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: obc_x2 = %d unknown\n",obc_x2);
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
	apply_ix3 = reflect_ix3_feedback;
	break;
      case 5: /* Reflecting */
	apply_ix3 = reflect_ix3_feedback;
	break;

      case 2: /* Outflow */
	apply_ix3 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ix3 = periodic_ix3_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: ibc_x3 = %d unknown\n",ibc_x3);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox3 == NULL){

      obc_x3 = par_geti("grid","obc_x3");
      switch(obc_x3){

      case 1: /* Reflecting */
	apply_ox3 = reflect_ox3_feedback;
	break;
      case 5: /* Reflecting */
	apply_ox3 = reflect_ox3_feedback;
	break;

      case 2: /* Outflow */
	apply_ox3 = outflow_feedback;
	break;

      case 4: /* Periodic */
	apply_ox3 = periodic_ox3_feedback;
	break;

      default:
	ath_perr(-1,"[exchange_feedback_init]: obc_x3 = %d unknown\n",obc_x3);
	exit(EXIT_FAILURE);
      }
    }
  }

#ifdef MPI_PARALLEL
  x1cnt = x2cnt = x3cnt = 0;

  for (k=0; k<(pD->NGrid_x3); k++){
    for (j=0; j<(pD->NGrid_x2); j++){
      for (i=0; i<(pD->NGrid_x1); i++){
	if(pD->NGrid_x3 > 1){
	  nx1t = pD->GridArray[k][j][i].ige - pD->GridArray[k][j][i].igs + 1;
	  if(nx1t > 1) nx1t += 2*NGF;

	  nx2t = pD->GridArray[k][j][i].jge - pD->GridArray[k][j][i].jgs + 1;
	  if(nx2t > 1) nx2t += 2*NGF;

	  x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
	}

	if(pD->NGrid_x2 > 1){
	  nx1t = pD->GridArray[k][j][i].ige - pD->GridArray[k][j][i].igs + 1;
	  if(nx1t > 1) nx1t += 2*NGF;

	  nx3t = pD->GridArray[k][j][i].kge - pD->GridArray[k][j][i].kgs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
	}

	if(pD->NGrid_x1 > 1){
	  nx2t = pD->GridArray[k][j][i].jge - pD->GridArray[k][j][i].jgs + 1;
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GridArray[k][j][i].kge - pD->GridArray[k][j][i].kgs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
	}

      }
    }
  }

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= NGF; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double*)malloc(size*NVAR_F*sizeof(double))) == NULL)
      ath_error("[exchange_feedback_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*NVAR_F*sizeof(double))) == NULL)
      ath_error("[exchange_feedback_init]: Failed to allocate receive buffer\n");
  }
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* exchange_feedback_fun:  sets function pointers for user-defined feedback exchange 
   in problem file
 */

void exchange_feedback_fun(enum Direction dir, VBCFun_t prob_bc)
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
    ath_perr(-1,"[set_bvals_particle_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/* finalize feedback exchange */
void exchange_feedback_destruct(Grid *pG, Domain *pD)
{
  apply_ix1 = NULL;
  apply_ox1 = NULL;
  apply_ix2 = NULL;
  apply_ox2 = NULL;
  apply_ix3 = NULL;
  apply_ox3 = NULL;
#ifdef MPI_PARALLEL
  free(send_buf);
  free(recv_buf);
#endif
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
/* REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1,5)
 */

static void reflect_ix3_feedback(Grid *pG)
{
  int kr;
  int i,j,k;

  for (k=nghost-NGF; k<nghost; k++) {
    kr = 2*nghost-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[kr][j][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[kr][j][i].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[kr][j][i].x3 -= pG->feedback[k][j][i].x3; /* reflect 3-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1,5)
 */

static void reflect_ox3_feedback(Grid *pG)
{
  int kr;
  int i,j,k;

  for (k=pG->ke+1; k<=pG->ke+NGF; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[kr][j][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[kr][j][i].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[kr][j][i].x3 -= pG->feedback[k][j][i].x3; /* reflect 3-mom. */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1,5)
 */

static void reflect_ix2_feedback(Grid *pG)
{
  int jr;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=nghost-NGF; j<nghost; j++) {
      jr = 2*nghost-j-1;
      for (i=il; i<=iu; i++) {
        pG->feedback[k][jr][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[k][jr][i].x2 -= pG->feedback[k][j][i].x2; /* reflect 2-mom. */
        pG->feedback[k][jr][i].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1,5)
 */

static void reflect_ox2_feedback(Grid *pG)
{
  int jr;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->je+1; j<=pG->je+NGF; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        pG->feedback[k][jr][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[k][jr][i].x2 -= pG->feedback[k][j][i].x2; /* reflect 2-mom. */
        pG->feedback[k][jr][i].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1,5)
 */

static void reflect_ix1_feedback(Grid *pG)
{
  int ir;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=nghost-NGF; i<nghost; i++) {
        ir = 2*nghost-i-1;
        pG->feedback[k][j][ir].x1 -= pG->feedback[k][j][i].x1; /* reflect 1-mom. */
        pG->feedback[k][j][ir].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[k][j][ir].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1,5)
 */

static void reflect_ox1_feedback(Grid *pG)
{
  int ir;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=pG->ie+1; i<=pG->ie+NGF; i++) {
        ir = 2*pG->ie-i+1;
        pG->feedback[k][j][ir].x1 -= pG->feedback[k][j][i].x1; /* reflect 1-mom. */
        pG->feedback[k][j][ir].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[k][j][ir].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions (ibc=1,5), essentially do nothing
 */

static void outflow_feedback(Grid *pG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x3 boundary (ibc_x3=4)
 */

static void periodic_ix3_feedback(Grid *pG)
{
  int dk = pG->ke - pG->ks + 1;
  int i,j,k;

  for (k=nghost; k<nghost+NGF; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[k][j][i].x1 += pG->feedback[k+dk][j][i].x1;
        pG->feedback[k][j][i].x2 += pG->feedback[k+dk][j][i].x2;
        pG->feedback[k][j][i].x3 += pG->feedback[k+dk][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Outer x3 boundary (obc_x3=4)
 */

static void periodic_ox3_feedback(Grid *pG)
{
  int dk = pG->ke - pG->ks + 1;
  int i,j,k;

  for (k=nghost-NGF; k<nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[k+dk][j][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[k+dk][j][i].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[k+dk][j][i].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2_feedback(Grid *pG)
{
  int dj = pG->je - pG->js + 1;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=nghost; j<nghost+NGF; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[k][j][i].x1 += pG->feedback[k][j+dj][i].x1;
        pG->feedback[k][j][i].x2 += pG->feedback[k][j+dj][i].x2;
        pG->feedback[k][j][i].x3 += pG->feedback[k][j+dj][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions,Outer x2 boundary (obc_x2=4)
 */

static void periodic_ox2_feedback(Grid *pG)
{
  int dj = pG->je - pG->js + 1;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=nghost-NGF; j<nghost; j++) {
      for (i=il; i<=iu; i++) {
        pG->feedback[k][j+dj][i].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[k][j+dj][i].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[k][j+dj][i].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_ix1_feedback(Grid *pG)
{
  int di = pG->ie - pG->is + 1;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=nghost; i<nghost+NGF; i++) {
        pG->feedback[k][j][i].x1 += pG->feedback[k][j][i+di].x1;
        pG->feedback[k][j][i].x2 += pG->feedback[k][j][i+di].x2;
        pG->feedback[k][j][i].x3 += pG->feedback[k][j][i+di].x3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Outer x1 boundary (obc_x1=4)
 */

static void periodic_ox1_feedback(Grid *pG)
{
  int di = pG->ie - pG->is + 1;
  int i,j,k;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=nghost-NGF; i<nghost; i++) {
        pG->feedback[k][j][i+di].x1 += pG->feedback[k][j][i].x1;
        pG->feedback[k][j][i+di].x2 += pG->feedback[k][j][i].x2;
        pG->feedback[k][j][i+di].x3 += pG->feedback[k][j][i].x3;
      }
    }
  }

  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~400 lines */

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x3 boundary -- send left
 */

static void send_ix3_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (iu-il+1)*(ju-jl+1)*NGF*NVAR_F;

  for (k=nghost-NGF; k<nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x3 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx3_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix3_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x3 boundary -- send right
 */

static void send_ox3_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (iu-il+1)*(ju-jl+1)*NGF*NVAR_F;

  for (k=pG->ke+1; k<=pG->ke+NGF; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x3 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx3_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox3_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x2 boundary -- send left
 */

static void send_ix2_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (iu-il+1)*(pG->ke-pG->ks+1)*NGF*NVAR_F;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=nghost-NGF; j<nghost; j++) {
      for (i=il; i<=iu; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x2 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix2_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x2 boundary -- send right
 */

static void send_ox2_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (iu-il+1)*(pG->ke-pG->ks+1)*NGF*NVAR_F;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->je+1; j<=pG->je+NGF; j++) {
      for (i=il; i<=iu; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x2 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox2_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x1 boundary -- send left
 */

static void send_ix1_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (pG->je-pG->js+1)*(pG->ke-pG->ks+1)*NGF*NVAR_F;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=nghost-NGF; i<nghost; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x1 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x1 boundary -- send right
 */

static void send_ox1_feedback(Grid *pG)
{
  int i,j,k,cnt,err;
  Vector *pq;
  double *pd = send_buf;

/* Pack feedback data into send buffer */
  cnt = (pG->je-pG->js+1)*(pG->ke-pG->ks+1)*NGF*NVAR_F;

  for (k=pG->ks; k<=pG->ke; k++) {
    for (j=pG->js; j<=pG->je; j++) {
      for (i=pG->ie+1; i<=pG->ie+NGF; i++) {

        pq = &(pG->feedback[k][j][i]);

        *(pd++) = pq->x1;
        *(pd++) = pq->x2;
        *(pd++) = pq->x3;
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x1 */
  err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id, boundary_cells_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1_feedback]: MPI_Send error = %d\n",err);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x3 boundary -- listen left
 */

static void recv_ix3_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the left grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix3_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ks; k<pG->ks+NGF; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outerer x3 boundary -- listen right
 */

static void recv_ox3_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the right grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox3_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ke-NGF+1; k<=pG->ke; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x2 boundary -- listen left
 */

static void recv_ix2_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the left grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix2_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ks; k<=pG->ke; k++){
    for (j=pG->js; k<pG->js+NGF; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x2 boundary -- listen right
 */

static void recv_ox2_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the right grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox2_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ks; k<=pG->ke; k++){
    for (j=pG->je-NGF+1; k<=pG->je; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x1 boundary -- listen left
 */

static void recv_ix1_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the left grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ix1_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ks; k<=pG->ke; k++){
    for (j=pG->js; k<=pG->je; j++){
      for (i=pG->is; i<pG->is+NGF; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x1 boundary -- listen right
 */

static void recv_ox1_feedback(Grid *pG, MPI_Request *prq)
{
  int i,j,k,err;
  Vector *pq;
  MPI_Status stat;
  double *pd = recv_buf;

/* Wait to receive the input data from the left grid */
  err = MPI_Wait(prq, &stat);
  if(err) ath_error("[receive_ox1_feedback]: MPI_Wait error = %d\n",err);

/* Manually unpack the data from the receive buffer */
  for (k=pG->ks; k<=pG->ke; k++){
    for (j=pG->js; k<=pG->je; j++){
      for (i=pG->ie-NGF+1; i<=pG->ie; i++){
        /* Get a pointer to the Gas cell */
        pq = &(pG->feedback[k][j][i]);

        pq->x1 += *(pd++);
        pq->x2 += *(pd++);
        pq->x3 += *(pd++);

      }
    }
  }

  return;
}

#endif /* MPI_PARALLEL */

#ifdef SHEARING_BOX

/*----------------------------------------------------------------------------*/
/* exchange feedback for 3D shearing box. Inner x1.
   Note that here we only need to shear the feedback array by an integer number
   of cells.
 */
static void shearingbox_ix1_feedback(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* exchange feedback for 3D shearing box. Outer x1.
   Note that here we only need to shear the feedback array by an integer number
   of cells.
 */
static void shearingbox_ox1_feedback(Grid *pG, Domain *pD)
{
  return;
}

#endif /* SHEARING_BOX */

#endif /* FEEDBACK */
