#include "copyright.h"
/*==============================================================================
 * FILE: set_bvals.c
 *
 * PURPOSE: Sets boundary conditions on each edge of the computational volume
 *   using integer flags.  The naming convention for the flags is:
 *     ibc_x1 = Inner Boundary Condition for x1
 *     obc_x1 = Outer Boundary Condition for x1
 *     similarly for ibc_x2; obc_x2; ibc_x3; obc_x3
 *   The values of the integer flags are:
 *     1 = reflecting; 2 = outflow; 4 = periodic
 *   Following ZEUS conventions, 3 would be flow-in (ghost zones held at
 *   pre-determined fixed values), however this is currently not implemented.
 *   Instead, use the boundary value function pointers for flow-in.
 *
 *   Routines to handle swapping of ghost zones between MPI patches are also
 *   included here.  
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   set_bvals()
 *   set_bvals_init()
 *   set_bvals_fun()
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#define MPI_USE_SEND_RECV

#ifdef MPI_PARALLEL

/* MPI send and receive buffer, size dynamically determined based on number of
 * zones in each grid */
static double *send_buf = NULL, *recv_buf = NULL;
static int size;

#ifdef MHD
/* Number of variables passed with MPI, including 3 extra for the interface
 * magnetic fields */
#define NVAR_SHARE (NVAR + 3)
#else
#define NVAR_SHARE NVAR
#endif
#endif /* MPI_PARALLEL */

/* boundary condition function pointers for edges of grid */
static VGFun_t apply_ix1_bc = NULL, apply_ox1_bc = NULL;
static VGFun_t apply_ix2_bc = NULL, apply_ox2_bc = NULL;
static VGFun_t apply_ix3_bc = NULL, apply_ox3_bc = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   pack_send()
 *   unpack_set()
 *   bc_send_recv()
 *   reflect_*()    - apply reflecting BCs at boundary *
 *   in_out_*
 *   periodic_*()
 *============================================================================*/

#ifdef MPI_PARALLEL
static void pack_send(Grid *pG, int to_id, 
            int is, int ie, int js, int je, int ks, int ke, MPI_Request *prq);
static void unpack_set(Grid *pG, int fr_id,
		       int is, int ie, int js, int je, int ks, int ke);
static void bc_send_recv(Grid *pG, int to_id, int s_is, int s_ie,
			 int s_js, int s_je, int s_ks, int s_ke,
			 int fr_id, int r_is, int r_ie,
			 int r_js, int r_je, int r_ks, int r_ke);
#endif


static void reflect_ix1(Grid *pG);
static void reflect_ox1(Grid *pG);
static void reflect_ix2(Grid *pG);
static void reflect_ox2(Grid *pG);
static void reflect_ix3(Grid *pG);
static void reflect_ox3(Grid *pG);

static void in_out_ix1(Grid *pG);
static void in_out_ox1(Grid *pG);
static void in_out_ix2(Grid *pG);
static void in_out_ox2(Grid *pG);
static void in_out_ix3(Grid *pG);
static void in_out_ox3(Grid *pG);

static void periodic_ix1(Grid *pG);
static void periodic_ox1(Grid *pG);
static void periodic_ix2(Grid *pG);
static void periodic_ox2(Grid *pG);
static void periodic_ix3(Grid *pG);
static void periodic_ox3(Grid *pG);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* set_bvals:  */

#ifndef MPI_PARALLEL

/* Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly */

void set_bvals(Grid *pGrid)
{

  if (pGrid->Nx1 > 1){
    (*apply_ix1_bc)(pGrid);
    (*apply_ox1_bc)(pGrid);
  }

  if (pGrid->Nx2 > 1){
    (*apply_ix2_bc)(pGrid);
    (*apply_ox2_bc)(pGrid);
  }

  if (pGrid->Nx3 > 1){
    (*apply_ix3_bc)(pGrid);
    (*apply_ox3_bc)(pGrid);
  }

  return;
}

#else

/* ========================================================================== */



static void pack_send(Grid *pG, int to_id, 
		      int is, int ie, int js, int je, int ks, int ke,
		      MPI_Request *prq){
  int i, j, k;
  Gas *pq;
  double *pd = send_buf;
  int cnt = (ie-is+1)*(je-js+1)*(ke-ks+1)*NVAR_SHARE;
  int err;

  for(k = ks; k <= ke; k++){
    for(j = js; j <= je; j++){
      for(i = is; i <= ie; i++){
	/* Get a pointer to the Gas cell */
	pq = &(pG->U[k][j][i]);

	/* Manually pack the data into the buffer send_buf */
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
#ifndef ISOTHERMAL
	*(pd++) = pq->E;
#endif /* ISOTHERMAL */
      }
    }
  }

  /* send this info to the neighboring grid */
  err = MPI_Isend(send_buf, cnt, MPI_DOUBLE, to_id,
		  boundary_cells_tag, MPI_COMM_WORLD, prq);
  if(err) ath_error("[set_bvals/pack_send]: MPI_Isend error = %d\n",err);

  return;
}


/* ========================================================================== */


static void unpack_set(Grid *pG, int fr_id,
		       int is, int ie, int js, int je, int ks, int ke){
  int i, j, k;
  Gas *pq;
  MPI_Status stat;
  double *pd = recv_buf;
  int cnt = (ie-is+1)*(je-js+1)*(ke-ks+1)*NVAR_SHARE;
  int err;

  /* Wait to receive the input data from the left grid */
  err = MPI_Recv(recv_buf, cnt, MPI_DOUBLE, fr_id,
		 boundary_cells_tag, MPI_COMM_WORLD, &stat);
  if(err) ath_error("[set_bvals/unpack_set]: MPI_Recv error = %d\n",err);

  for(k = ks; k <= ke; k++){
    for(j = js; j <= je; j++){
      for(i = is; i <= ie; i++){
	/* Get a pointer to the Gas cell */
	pq = &(pG->U[k][j][i]);

	/* Manually unpack the data from the buffer recv_buf */
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
#ifndef ISOTHERMAL
	pq->E = *(pd++);
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}


/* ========================================================================== */

#ifdef MPI_USE_SEND_RECV

static void bc_send_recv(Grid *pG, int to_id, int s_is, int s_ie,
			 int s_js, int s_je, int s_ks, int s_ke,
			 int fr_id, int r_is, int r_ie,
			 int r_js, int r_je, int r_ks, int r_ke){
  int i, j, k;
  Gas *pq;
  MPI_Status stat;
  int s_cnt = (s_ie-s_is+1)*(s_je-s_js+1)*(s_ke-s_ks+1)*NVAR_SHARE;
  int r_cnt = (r_ie-r_is+1)*(r_je-r_js+1)*(r_ke-r_ks+1)*NVAR_SHARE;
  double *pd;
  int err;

  pd = send_buf;
  for(k = s_ks; k <= s_ke; k++){
    for(j = s_js; j <= s_je; j++){
      for(i = s_is; i <= s_ie; i++){
	/* Get a pointer to the Gas cell */
	pq = &(pG->U[k][j][i]);

	/* Manually pack the data into the buffer send_buf */
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
#ifndef ISOTHERMAL
	*(pd++) = pq->E;
#endif /* ISOTHERMAL */
      }
    }
  }

  /* Both send the data to one neighbor and receive data from the other */
  err = MPI_Sendrecv(send_buf, s_cnt, MPI_DOUBLE, to_id, boundary_cells_tag,
		     recv_buf, r_cnt, MPI_DOUBLE, fr_id, boundary_cells_tag,
		     MPI_COMM_WORLD, &stat);
  if(err)
    ath_error("[set_bvals/bc_send_recv]: MPI_Sendrecv error = %d\n",err);

  pd = recv_buf;
  for(k = r_ks; k <= r_ke; k++){
    for(j = r_js; j <= r_je; j++){
      for(i = r_is; i <= r_ie; i++){
	/* Get a pointer to the Gas cell */
	pq = &(pG->U[k][j][i]);

	/* Manually unpack the data from the buffer recv_buf */
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
#ifndef ISOTHERMAL
	pq->E = *(pd++);
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

#endif /* MPI_USE_SEND_RECV */

/* ========================================================================== */


/* Apply boundary conditions by setting the ghost cells */
/* We do the parallel synchronization by having every grid:
   1) Pack up and send the appropriate data to the grid to its right
   2) Listen to the left, unpack the data, and set the values
   3) Pack up and send the appropriate data to the grid to its left
   4) Listen to the right, unpack the data, and set the values

   If there is no grid to either side we set our boundary conditions
   in the standard way on the listen step.
*/
void set_bvals(Grid *pG){

  int il,iu,jl,ju,kl,ku;
#ifdef MPI_USE_SEND_RECV
  int sl, su, rl, ru; /* The prefix s and r -> send and receive */
#endif
  MPI_Status stat;
  MPI_Request rq;
  int err;

  /* x_dir_bc: */
  if(pG->Nx1 == 1) goto y_dir_bc;

  if(pG->Nx2 > 1){
    jl = pG->js;
    ju = pG->je + 1;
  }
  else jl = ju = pG->js;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  }
  else kl = ku = pG->ks;

  /* ======================================================================== */

#ifdef MPI_USE_SEND_RECV
  if(pG->rx1_id >= 0 && pG->lx1_id >= 0 ){ /* We use send_recv */
    /* Pack up and send right & Listen left, unpack and set values */
    sl = pG->ie - nghost + 1;
    su = pG->ie;

    rl = pG->is - nghost;
    ru = pG->is - 1;

    /* printf("my_id = %d, calling bc_send_recv (right-x)\n",pG->my_id); */
    bc_send_recv(pG, pG->rx1_id, sl, su, jl, ju, kl, ku,
		     pG->lx1_id, rl, ru, jl, ju, kl, ku);

    /* Pack up and send left & Listen right, unpack and set values */
    sl = pG->is;
    su = pG->is + nghost - 1;

    rl = pG->ie + 1;
    ru = pG->ie + nghost;

    /* printf("my_id = %d, calling bc_send_recv (left-x)\n",pG->my_id); */

    bc_send_recv(pG, pG->lx1_id, sl, su, jl, ju, kl, ku,
		     pG->rx1_id, rl, ru, jl, ju, kl, ku);
  }
  else{
#endif
    /* Pack up and send right */
    if(pG->rx1_id >= 0){
      il = pG->ie - nghost + 1;
      iu = pG->ie;
      pack_send(pG,pG->rx1_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen left, unpack and set values */
    if(pG->lx1_id < 0){ /* There is no grid to the left */
      (*apply_ix1_bc)(pG);
    }
    else{ /* There is a grid to the left */
      il = pG->is - nghost;
      iu = pG->is - 1;
      unpack_set(pG,pG->lx1_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->rx1_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }

    /* ==================================================================== */

    /* Pack up and send left */
    if(pG->lx1_id >= 0){
      il = pG->is;
      iu = pG->is + nghost - 1;
      pack_send(pG,pG->lx1_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen right, unpack and set values */
    if(pG->rx1_id < 0){ /* There is no grid to the right */
      (*apply_ox1_bc)(pG);
    }
    else{ /* There is a grid to the right */
      il = pG->ie + 1;
      iu = pG->ie + nghost;
      unpack_set(pG,pG->rx1_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->lx1_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }
#ifdef MPI_USE_SEND_RECV
  }
#endif

  /* ======================================================================== */

 y_dir_bc:
  if(pG->Nx2 == 1) goto z_dir_bc;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  }
  else il = iu = pG->is;

  if(pG->Nx3 > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  }
  else kl = ku = pG->ks;

  /* ======================================================================== */

#ifdef MPI_USE_SEND_RECV
  if(pG->rx2_id >= 0 && pG->lx2_id >= 0 ){ /* We use send_recv */
    /* Pack up and send right & Listen left, unpack and set values */
    sl = pG->je - nghost + 1;
    su = pG->je;

    rl = pG->js - nghost;
    ru = pG->js - 1;

    /* printf("my_id = %d, calling bc_send_recv (right-y)\n",pG->my_id); */
    bc_send_recv(pG, pG->rx2_id, il, iu, sl, su, kl, ku,
		     pG->lx2_id, il, iu, rl, ru, kl, ku);

    /* Pack up and send left & Listen right, unpack and set values */
    sl = pG->js;
    su = pG->js + nghost - 1;

    rl = pG->je + 1;
    ru = pG->je + nghost;

    /* printf("my_id = %d, calling bc_send_recv (left-y)\n",pG->my_id); */
    bc_send_recv(pG, pG->lx2_id, il, iu, sl, su, kl, ku,
		     pG->rx2_id, il, iu, rl, ru, kl, ku);
  }
  else{
#endif
    /* Pack up and send right */
    if(pG->rx2_id >= 0){
      jl = pG->je - nghost + 1;
      ju = pG->je;
      pack_send(pG,pG->rx2_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen left, unpack and set values */
    if(pG->lx2_id < 0){ /* There is no grid to the left */
      (*apply_ix2_bc)(pG);
    }
    else{ /* There is a grid to the left */
      jl = pG->js - nghost;
      ju = pG->js - 1;
      unpack_set(pG,pG->lx2_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->rx2_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }

    /* ==================================================================== */

    /* Pack up and send left */
    if(pG->lx2_id >= 0){
      jl = pG->js;
      ju = pG->js + nghost - 1;
      pack_send(pG,pG->lx2_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen right, unpack and set values */
    if(pG->rx2_id < 0){ /* There is no grid to the right */
      (*apply_ox2_bc)(pG);
    }
    else{ /* There is a grid to the right */
      jl = pG->je + 1;
      ju = pG->je + nghost;
      unpack_set(pG,pG->rx2_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->lx2_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }
#ifdef MPI_USE_SEND_RECV
  }
#endif

  /* ======================================================================== */

 z_dir_bc:
  if(pG->Nx3 == 1) return;

  if(pG->Nx1 > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  }
  else il = iu = pG->is;

  if(pG->Nx2 > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  }
  else jl = ju = pG->js;

  /* ======================================================================== */

#ifdef MPI_USE_SEND_RECV
  if(pG->rx3_id >= 0 && pG->lx3_id >= 0 ){ /* We use send_recv */
    /* Pack up and send right & Listen left, unpack and set values */
    sl = pG->ke - nghost + 1;
    su = pG->ke;

    rl = pG->ks - nghost;
    ru = pG->ks - 1;

    bc_send_recv(pG, pG->rx3_id, il, iu, jl, ju, sl, su,
		     pG->lx3_id, il, iu, jl, ju, rl, ru);

    /* Pack up and send left & Listen right, unpack and set values */
    sl = pG->ks;
    su = pG->ks + nghost - 1;

    rl = pG->ke + 1;
    ru = pG->ke + nghost;

    bc_send_recv(pG, pG->lx3_id, il, iu, jl, ju, sl, su,
		     pG->rx3_id, il, iu, jl, ju, rl, ru);
  }
  else{
#endif
    /* Pack up and send right */
    if(pG->rx3_id >= 0){
      kl = pG->ke - nghost + 1;
      ku = pG->ke;
      pack_send(pG,pG->rx3_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen left, unpack and set values */
    if(pG->lx3_id < 0){ /* There is no grid to the left */
      (*apply_ix3_bc)(pG);
    }
    else{ /* There is a grid to the left */
      kl = pG->ks - nghost;
      ku = pG->ks - 1;
      unpack_set(pG,pG->lx3_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->rx3_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }

    /* ==================================================================== */

    /* Pack up and send left */
    if(pG->lx3_id >= 0){
      kl = pG->ks;
      ku = pG->ks + nghost - 1;
      pack_send(pG,pG->lx3_id,il,iu,jl,ju,kl,ku,&rq);
    }

    /* Listen right, unpack and set values */
    if(pG->rx3_id < 0){ /* There is no grid to the right */
      (*apply_ox3_bc)(pG);
    }
    else{ /* There is a grid to the right */
      kl = pG->ke + 1;
      ku = pG->ke + nghost;
      unpack_set(pG,pG->rx3_id,il,iu,jl,ju,kl,ku);
    }

    if(pG->lx3_id >= 0){ /* Ensure send is complete */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[set_bvals]: MPI_Wait error = %d\n",err);
    }
#ifdef MPI_USE_SEND_RECV
  }
#endif

  return;
}


/* ========================================================================== */


#endif /* MPI_PARALLEL */

/* ========================================================================== */



void set_bvals_fun(enum Direction dir, VGFun_t prob_bc){

  switch(dir){
  case left_x1:
    apply_ix1_bc = prob_bc;
    break;
  case right_x1:
    apply_ox1_bc = prob_bc;
    break;
  case left_x2:
    apply_ix2_bc = prob_bc;
    break;
  case right_x2:
    apply_ox2_bc = prob_bc;
    break;
  case left_x3:
    apply_ix3_bc = prob_bc;
    break;
  case right_x3:
    apply_ox3_bc = prob_bc;
    break;
  default:
    fprintf(stderr,"[set_bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}


/* ========================================================================== */


void set_bvals_init(Grid *pG)
{
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
#ifdef MPI_PARALLEL
  int i, j, k;
  int my_id = pG->my_id;
  int x1cnt, x2cnt, x3cnt; /* Number of Gas passed in x1-, x2-, x3-dir. */
  int nx1t, nx2t, nx3t;
#endif /* MPI_PARALLEL */

  if(pG->Nx1 > 1) {
    if(apply_ix1_bc == NULL){
      ibc_x1 = par_geti("grid","ibc_x1");
      switch(ibc_x1){
      case 1: /* Reflecting */
	apply_ix1_bc = reflect_ix1;
	break;
      case 2: /* Inflow / Outflow */
	apply_ix1_bc = in_out_ix1;
	break;
      case 4: /* Periodic */
	apply_ix1_bc = periodic_ix1;
#ifdef MPI_PARALLEL
	if(pG->lx1_id < 0 && NGrid_x1 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->lx1_id = grid_domain[k][j][NGrid_x1-1].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: ibc_x1 = %d unknown\n",ibc_x1);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox1_bc == NULL){
      obc_x1 = par_geti("grid","obc_x1");
      switch(obc_x1){
      case 1: /* Reflecting */
	apply_ox1_bc = reflect_ox1;
	break;
      case 2: /* Inflow / Outflow */
	apply_ox1_bc = in_out_ox1;
	break;
      case 4: /* Periodic */
	apply_ox1_bc = periodic_ox1;
#ifdef MPI_PARALLEL
	if(pG->rx1_id < 0 && NGrid_x1 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->rx1_id = grid_domain[k][j][0].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: obc_x1 = %d unknown\n",obc_x1);
	exit(EXIT_FAILURE);
      }
    }
  }

  if(pG->Nx2 > 1) {
    if(apply_ix2_bc == NULL){
      ibc_x2 = par_geti("grid","ibc_x2");
      switch(ibc_x2){
      case 1: /* Reflecting */
	apply_ix2_bc = reflect_ix2;
	break;
      case 2: /* Inflow / Outflow */
	apply_ix2_bc = in_out_ix2;
	break;
      case 4: /* Periodic */
	apply_ix2_bc = periodic_ix2;
#ifdef MPI_PARALLEL
	if(pG->lx2_id < 0 && NGrid_x2 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->lx2_id = grid_domain[k][NGrid_x2-1][i].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: ibc_x2 = %d unknown\n",ibc_x2);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox2_bc == NULL){
      obc_x2 = par_geti("grid","obc_x2");
      switch(obc_x2){
      case 1: /* Reflecting */
	apply_ox2_bc = reflect_ox2;
	break;
      case 2: /* Inflow / Outflow */
	apply_ox2_bc = in_out_ox2;
	break;
      case 4: /* Periodic */
	apply_ox2_bc = periodic_ox2;
#ifdef MPI_PARALLEL
	if(pG->rx2_id < 0 && NGrid_x2 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->rx2_id = grid_domain[k][0][i].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: obc_x2 = %d unknown\n",obc_x2);
	exit(EXIT_FAILURE);
      }
    }
  }

  if(pG->Nx3 > 1) {
    if(apply_ix3_bc == NULL){
      ibc_x3 = par_geti("grid","ibc_x3");
      switch(ibc_x3){
      case 1: /* Reflecting */
	apply_ix3_bc = reflect_ix3;
	break;
      case 2: /* Inflow / Outflow */
	apply_ix3_bc = in_out_ix3;
	break;
      case 4: /* Periodic */
	apply_ix3_bc = periodic_ix3;
#ifdef MPI_PARALLEL
	if(pG->lx3_id < 0 && NGrid_x3 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->lx3_id = grid_domain[NGrid_x3-1][j][i].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: ibc_x3 = %d unknown\n",ibc_x3);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox3_bc == NULL){
      obc_x3 = par_geti("grid","obc_x3");
      switch(obc_x3){
      case 1: /* Reflecting */
	apply_ox3_bc = reflect_ox3;
	break;
      case 2: /* Inflow / Outflow */
	apply_ox3_bc = in_out_ox3;
	break;
      case 4: /* Periodic */
	apply_ox3_bc = periodic_ox3;
#ifdef MPI_PARALLEL
	if(pG->rx3_id < 0 && NGrid_x3 > 1){
	  domain_ijk(my_id, &i, &j, &k);
	  pG->rx3_id = grid_domain[0][j][i].my_id;
	}
#endif /* MPI_PARALLEL */
	break;
      default:
	fprintf(stderr,"[set_bvals_init]: obc_x3 = %d unknown\n",obc_x3);
	exit(EXIT_FAILURE);
      }
    }
  }

#ifdef MPI_PARALLEL
  x1cnt = x2cnt = x3cnt = 0;
  for(k=0; k<NGrid_x3; k++){
    for(j=0; j<NGrid_x2; j++){
      for(i=0; i<NGrid_x1; i++){
	/* x1-dir. data pass */
	if(NGrid_x1 > 1){
	  nx2t = grid_domain[k][j][i].jxe - grid_domain[k][j][i].jxs + 1;
	  if(nx2t > 1) nx2t += 1;

	  nx3t = grid_domain[k][j][i].kxe - grid_domain[k][j][i].kxs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
	  /* printf("nx2t = %d, nx3t = %d, x1cnt = %d\n",nx2t, nx3t, x1cnt); */
	}

	/* x2-dir. data pass */
	if(NGrid_x2 > 1){
	  nx1t = grid_domain[k][j][i].ixe - grid_domain[k][j][i].ixs + 1;
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx3t = grid_domain[k][j][i].kxe - grid_domain[k][j][i].kxs + 1;
	  if(nx3t > 1) nx3t += 1;

	  x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
	  /* printf("nx1t = %d, nx3t = %d, x2cnt = %d\n",nx1t, nx3t, x2cnt); */
	}


	/* x3-dir. data pass */
	if(NGrid_x3 > 1){
	  nx1t = grid_domain[k][j][i].ixe - grid_domain[k][j][i].ixs + 1;
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx2t = grid_domain[k][j][i].jxe - grid_domain[k][j][i].jxs + 1;
	  if(nx2t > 1) nx2t += 2*nghost;

	  x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
	  /* printf("nx1t = %d, nx2t = %d, x3cnt = %d\n",nx1t, nx2t, x3cnt); */
	}
      }
    }
  }

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= nghost; /* Multiply by the third dimension */

  if(size > 0){
    /* printf("Allocating two vectors of size %d doubles\n",size*NVAR_SHARE); */
    if((send_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate receive buffer\n");
  }
#endif /* MPI_PARALLEL */

  return;
}


/* ====================== Start with Inner x1 boundary ====================== */


static void periodic_ix1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;

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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ix1(Grid *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;

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
#endif

  return;
}


/* ========================================================================== */


static void reflect_ix1(Grid *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;

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
#endif

  return;
}


/* ===================== Continue with Outer x1 boundary ==================== */


static void periodic_ox1(Grid *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;

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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ox1(Grid *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;

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
#endif

  return;
}


/* ========================================================================== */


static void reflect_ox1(Grid *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ju,ku;


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
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie-(i-1)];
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
#endif

  return;
}


/* ==================== Next comes the Inner x2 boundary ==================== */


static void periodic_ix2(Grid *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ix2(Grid *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
#endif

  return;
}


/* ========================================================================== */


static void reflect_ix2(Grid *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
#endif

  return;
}


/* ==================== Next comes the Outer x2 boundary ==================== */


static void periodic_ox2(Grid *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ox2(Grid *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
#endif

  return;
}


/* ========================================================================== */


static void reflect_ox2(Grid *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,ku; /* i-lower/upper;  k-upper */

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
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je-(j-1)][i];
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
#endif

  return;
}

/* ==================== Next comes the Inner x3 boundary ==================== */


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
        pGrid->U  [ks-k][j][i] = pGrid->U  [ke-(k-1)][j][i];
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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ix3(Grid *pGrid)
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
#endif

  return;
}


/* ========================================================================== */


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
#endif

  return;
}


/* ==================== Next comes the Outer x3 boundary ==================== */


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
#endif

  return;
}


/* ========================================================================== */


static void in_out_ox3(Grid *pGrid)
{
  int ke = pGrid->ke;
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
#endif

  return;
}


/* ========================================================================== */


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
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

