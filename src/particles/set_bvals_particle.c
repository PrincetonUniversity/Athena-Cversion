/*=============================================================================
FILE: set_bvals_particles.c
PURPOSE: Sets boundary conditions for particles. The basic algorithm is similar to
  the MHD boundary conditions. In setting B.C., particles are always packed to the
  send buffer (for both MPI and non-MPI cases), then unpacked with certain shifts
  in position or velocity. In the case of MPI, two communications are needed for
  sending and receiving particle, first on the number of particles to be sent/received,
  then to send/receive real data. The shearing box B.C. is first treated as periodic
  B.C., then move the particles in the ghost cells accordingly. Advection of particles
  when FARGO is turned on is also included.

CONTAINS PUBLIC FUNCTIONS:
  void set_bvals_particle(Grid *pG, Domain *pD);
  void advect_particles(Grid *pG, Domain *pD);
  void set_bvals_particle_init(Grid *pG, Domain *pD);
  void set_bvals_particle_fun(enum Direction dir, VBCFun_t prob_bc);
  void set_bvals_particle_destruct(Grid *pG);

History:
 Created:	Emmanuel Jacquet	May 2008
 Rewritten:	Xuening Bai		Feb. 2009
==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"


#ifdef PARTICLES         /* endif at the end of the file */

/* particle structure size */
#define NVAR_P 7

/* send and receive buffer, size dynamically determined
 * They are mainly used for MPI, and shearing box.
 */
static double *send_buf = NULL, *recv_buf = NULL;
static long NBUF;		/* buffer size unit (in number of particle) */
static long send_bufsize;	/* size of the send buffer (in unit of particles) */
static long recv_bufsize;	/* size of the recv buffer (in unit of particles) */

int my_iproc, my_jproc, my_kproc;	/* processor indices in the computational domain */
Real x1min,x1max,x2min,x2max,x3min,x3max;/* min and max coordinate limits of the computational domain */
Real Lx1, Lx2, Lx3;			/* domain size in x1, x2, x3 direction */

#ifdef SHEARING_BOX
Real vshear;				/* shear velocity */
#endif /* SHEARING_BOX */

/* boundary condition function pointers. local to this function  */
static VBCFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VBCFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VBCFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*====================== PROTOTYPE OF PRIVATE FUNCTIONS ======================*/
/*----------------------------------------------------------------------------*/
static void realloc_sendbuf();
static void realloc_recvbuf();

static void reflect_ix1_particle(Grid *pG);
static void reflect_ox1_particle(Grid *pG);
static void reflect_ix2_particle(Grid *pG);
static void reflect_ox2_particle(Grid *pG);
static void reflect_ix3_particle(Grid *pG);
static void reflect_ox3_particle(Grid *pG);

static void outflow_particle(Grid *pG);

static void periodic_ix1_particle(Grid *pG);
static void periodic_ox1_particle(Grid *pG);
static void periodic_ix2_particle(Grid *pG);
static void periodic_ox2_particle(Grid *pG);
static void periodic_ix3_particle(Grid *pG);
static void periodic_ox3_particle(Grid *pG);

static long packing_ix1_particle(Grid *pG, int nlayer);
static long packing_ox1_particle(Grid *pG, int nlayer);
static long packing_ix2_particle(Grid *pG, int nlayer);
static long packing_ox2_particle(Grid *pG, int nlayer);
static long packing_ix3_particle(Grid *pG, int nlayer);
static long packing_ox3_particle(Grid *pG, int nlayer);

static void shift_packed_particle(double *buf, long n, int index, double shift);
static void unpack_particle(Grid *pG, double *buf, long n);

#ifdef SHEARING_BOX
static void shearingbox_ix1_particle(Grid *pG, Domain *pD, long numpar);
static void shearingbox_ox1_particle(Grid *pG, Domain *pD, long numpar);

static long packing_ix1_particle_shear(Grid *pG, int reg, long numpar);
static long packing_ox1_particle_shear(Grid *pG, int reg, long numpar);

#ifdef FARGO
static long packing_particle_fargo(Grid *pG, Real yl, Real yu);
static int gridshift(Real shift);
#endif /* FARGO */

#endif /* SHEARING_BOX */


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* set_bvals_particle: calls appropriate functions to set ghost zones.  The function
 *   pointers (*apply_???) are set during initialization by set_bvals_particle_init()
 *   to be either a user-defined function, or one of the functions corresponding
 *   to reflecting, periodic, or outflow.  If the left- or right-Grid ID numbers
 *   are >= 1 (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void set_bvals_particle(Grid *pG, Domain *pD)
{
#ifdef MPI_PARALLEL
  int err;
  long cnt_send, cnt_recv;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */
#ifdef SHEARING_BOX
  long numpar;  /* number of particles before applying B.C. */
#endif

/*--- Step 1. ------------------------------------------------------------------
 * shuffle if necessary */

  /* shuffle every NSHUFFLE steps */
  if ((NShuffle>0) && (fmod(pG->nstep, NShuffle)<0.1))
    shuffle(pG);

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pG->Nx1 > 1){

#ifdef SHEARING_BOX
  numpar = pG->nparticle;
#endif

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx1_id >= 0 && pG->lx1_id >= 0) {
      /*-------- sub-step 1: send to right and receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox1_particle(pG, nghost);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (my_iproc == pD->NGrid_x1-1) {
        /* physical boundary on the rignt in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 1, -Lx1);
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
#ifndef FARGO
        /* velocity shift for shearing box */
        shift_packed_particle(send_buf, cnt_send, 6, vshear);
#endif
#endif /* SHEARING_BOX */
      }

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }

      /*-------- sub-step 2: send to left and receive from right --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix1_particle(pG, nghost);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (my_iproc == 0) {
        /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 1, Lx1);
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
#ifndef FARGO
        /* velocity shift for shearing box */
        shift_packed_particle(send_buf, cnt_send, 6, -vshear);
#endif
#endif /* SHEARING_BOX */
      }

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox1_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox1_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* Physical boundary on left, MPI block on right */
    if ((pG->rx1_id >= 0) && (pG->lx1_id < 0)) {
      /*-------- sub-step 1: send to right  --------*/
      /* packing particle on the right to send buffer */
      cnt_send = packing_ox1_particle(pG, nghost);
      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      }

     /*-------- sub-step 2: apply left boundary condition --------*/
      (*apply_ix1)(pG);

      /*-------- sub-step 3: receive from right  --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
     /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);
      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox1_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* MPI block on left, Physical boundary on right */
    if ((pG->lx1_id >= 0) && (pG->rx1_id < 0)) {
      /*-------- sub-step 1: receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);
      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      /*-------- sub-step 2: apply right boundary condition --------*/
      (*apply_ox1)(pG);

      /*-------- sub-step 3: send to left --------*/
      /* packing particle on the left to send buffer */
      cnt_send = packing_ix1_particle(pG, nghost);
      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx1_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox1_particle]: MPI_Send error = %d\n",err);
      }

      /*-------- sub-step 1: receive from left (cont) --------*/
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx1_id < 0 && pG->lx1_id < 0) {
      (*apply_ix1)(pG);
      (*apply_ox1)(pG);
    }

#ifdef SHEARING_BOX
    if (pG->Nx3>0) {
    /* For 3D shearing box boundary conditions */
      if (my_iproc == 0) /* inner boundary */
        shearingbox_ix1_particle(pG, pD, numpar);
      if (my_iproc == (pD->NGrid_x1-1)) /* outer boundary */
        shearingbox_ox1_particle(pG, pD, numpar);
    }
#endif /* SHEARING_BOX */

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pG->Nx2 > 1) {

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx2_id >= 0 && pG->lx2_id >= 0) {
      /*-------- sub-step 1: send to right and receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox2_particle(pG, nghost);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      if (my_jproc == pD->NGrid_x2-1) /* physical boundary on the rignt in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 2, -Lx2);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }

      /*-------- sub-step 2: send to left and receive from right --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix2_particle(pG, nghost);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      if (my_jproc == 0) /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 2, Lx2);

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox2_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox2_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx2_id >= 0 && pG->lx2_id < 0) {
      /*-------- sub-step 1: send to right --------*/
      /* packing particle on the right to send buffer */
      cnt_send = packing_ox2_particle(pG, nghost);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      }

      /*-------- sub-step 2: apply left boundary condition --------*/
      (*apply_ix2)(pG);

      /*-------- sub-step 3: receive from right --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
     /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox2_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx2_id < 0 && pG->lx2_id >= 0) {
      /*-------- sub-step 1: receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);
      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      /*-------- sub-step 2: apply right boundary condition --------*/
      (*apply_ox2)(pG);

      /*-------- sub-step 3: send to left --------*/
      /* packing particle on the left to send buffer */
      cnt_send = packing_ix2_particle(pG, nghost);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx2_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox2_particle]: MPI_Send error = %d\n",err);
      }

      /*-------- sub-step 1: receive from left (cont) --------*/
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx2_id < 0 && pG->lx2_id < 0) {
      (*apply_ix2)(pG);
      (*apply_ox2)(pG);
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pG->Nx3 > 1){

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx3_id >= 0 && pG->lx3_id >= 0) {
      /*-------- sub-step 1: send to right and receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox3_particle(pG, nghost);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      if (my_kproc == pD->NGrid_x3-1) /* physical boundary on the rignt in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 3, -Lx3);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }

      /*-------- sub-step 2: send to left and receive from right --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix3_particle(pG, nghost);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      if (my_kproc == 0) /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 3, Lx3);

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox3_particle]: MPI_Send error = %d\n",err);
      }
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox3_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx3_id >= 0 && pG->lx3_id < 0) {
      /*-------- sub-step 1: send to right --------*/
      /* packing particle on the right to send buffer */
      cnt_send = packing_ox3_particle(pG, nghost);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);
      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);
      }

      /*-------- sub-step 2: apply left boundary condition --------*/
      (*apply_ix3)(pG);

      /*-------- sub-step 3: receive from right --------*/
      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ox3_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx3_id < 0 && pG->lx3_id >= 0) {
      /*-------- sub-step 1: receive from left --------*/
      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);
      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf();
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      /*-------- sub-step 2: apply right boundary condition --------*/
      (*apply_ox3)(pG);

      /*-------- sub-step 3: send to left --------*/
      /* packing particle on the left to send buffer */
      cnt_send = packing_ix3_particle(pG, nghost);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);
      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx3_id, boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ox3_particle]: MPI_Send error = %d\n",err);
      }

      /*-------- sub-step 1: receive from left (cont) --------*/
      if (cnt_recv > 0) {
        /* receive buffer from the left grid */
        err = MPI_Wait(&rq, &stat);
        if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);
        /* unpack the received particle */
        unpack_particle(pG, recv_buf, cnt_recv);
      }
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx3_id < 0 && pG->lx3_id < 0) {
      (*apply_ix3)(pG);
      (*apply_ox3)(pG);
    } 
  }

  return;
}

#ifdef FARGO
/* Advect particles by 1.5*Omega*x1*dt for the FARGO algorithm.
*/
void advect_particles(Grid *pG, Domain *pD)
{
  Grain *cur;
  long p;
  Real x1l, x1u;
#ifdef MPI_PARALLEL
  long cnt_recv, n;
  int ishl, ishu, i;
  int inds, indr, ids, idr;
  Real x3len, yl, yu;
  int err;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */

  /* lower and upper bound of the grid in x1 direction */
  x1l = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;
  x1u = pG->x1_0 + (pG->ie+1 + pG->idisp)*pG->dx1;
  /* shift the particles */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    cur->x3 = x3min + fmod(cur->x3 + cur->shift - x3min + Lx3, Lx3);
  }

#ifdef MPI_PARALLEL
  /* calculate the farthest grid that advection can reach */
  x3len = pG->dx3*pG->Nx3;
  ishl = gridshift(1.5*Omega*(x1l-pG->dx1)*pG->dt/x3len);
  ishu = gridshift(1.5*Omega*(x1u+pG->dx1)*pG->dt/x3len);

  /* loop over all the possible destination grids */
  for (i=ishl; i<=ishu; i++)
  if (i != 0) { /* avoid moving particles to itself */
    /* find the processor id to send/receive data */
    inds = my_kproc + i;
    if (inds < 0) inds += pD->NGrid_x3;
    if (inds > (pD->NGrid_x3-1)) inds -= pD->NGrid_x3;
    ids = pD->GridArray[inds][my_jproc][my_iproc].id;	/* send to */

    indr = my_kproc - i;
    if (indr < 0) indr += pD->NGrid_x3;
    if (indr > (pD->NGrid_x3-1)) indr -= pD->NGrid_x3;
    idr = pD->GridArray[indr][my_jproc][my_iproc].id;	/* receive from */

    /* Post a non-blocking receive for the data size */
    err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, idr, boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

    /* packing particles */
    yl = x3min + inds*x3len;	yu = x3min + (inds+1)*x3len;
    n = packing_particle_fargo(pG, yl, yu);

    /* send buffer size */
    err = MPI_Send(&n, 1, MPI_LONG, ids, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
    /* receive buffer size */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);

    /* check space for the receive buffer */
    if ((cnt_recv+1) >= recv_bufsize)
      realloc_recvbuf();

    /* Post a non-blocking receive for data */
    if (cnt_recv > 0) {
      err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, idr, boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err); }
    /* send buffer */
    if (n > 0) {
      err = MPI_Send(send_buf, n*NVAR_P, MPI_DOUBLE, ids, boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
    }
    /* receive buffer */
    if (cnt_recv > 0) {
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);
      /* unpack the received particle */
      unpack_particle(pG, recv_buf, cnt_recv);
    }
  }
#endif /* MPI_PARALLEL */

  return;
}
#endif /* FARGO */

/*----------------------------------------------------------------------------*/
/* set_bvals_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */
void set_bvals_particle_init(Grid *pG, Domain *pD)
{
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
#ifdef MPI_PARALLEL
  int ib,jb,kb;
  int my_id = pG->my_id;
#endif /* MPI_PARALLEL */

/* initialize buffers */
  NBUF = (long)(0.15*pG->arrsize);

  send_bufsize = NBUF;
  recv_bufsize = NBUF;
  send_buf = (double*)calloc_1d_array(NVAR_P*send_bufsize, sizeof(double));
  recv_buf = (double*)calloc_1d_array(NVAR_P*recv_bufsize, sizeof(double));

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

  /* get the number of time steps for shuffle */
  NShuffle = par_geti("particle","nshuf");

#ifdef SHEARING_BOX
  vshear = 1.5 * Omega * Lx1;	/* shear velocity between inner and outer x1 boundaries */
#endif /* SHEARING_BOX */

/* Set function pointers for physical boundaries in x1-direction */

  if(pG->Nx1 > 1) {
    if(apply_ix1 == NULL){

      ibc_x1 = par_geti("grid","ibc_x1");
      switch(ibc_x1){

      case 1: /* Reflecting */
	apply_ix1 = reflect_ix1_particle;
	break;
      case 5: /* Reflecting */
	apply_ix1 = reflect_ix1_particle;
	break;

      case 2: /* Outflow */
	apply_ix1 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ix1 = periodic_ix1_particle;
#ifdef MPI_PARALLEL
	if(pG->lx1_id < 0 && pD->NGrid_x1 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx1_id = pD->GridArray[kb][jb][pD->NGrid_x1-1].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: ibc_x1 = %d unknown\n",ibc_x1);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox1 == NULL){

      obc_x1 = par_geti("grid","obc_x1");
      switch(obc_x1){

      case 1: /* Reflecting */
	apply_ox1 = reflect_ox1_particle;
	break;
      case 5: /* Reflecting */
	apply_ox1 = reflect_ox1_particle;
	break;

      case 2: /* Outflow */
	apply_ox1 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ox1 = periodic_ox1_particle;
#ifdef MPI_PARALLEL
	if(pG->rx1_id < 0 && pD->NGrid_x1 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx1_id = pD->GridArray[kb][jb][0].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: obc_x1 = %d unknown\n",obc_x1);
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
	apply_ix2 = reflect_ix2_particle;
	break;
      case 5: /* Reflecting */
	apply_ix2 = reflect_ix2_particle;
	break;

      case 2: /* Outflow */
	apply_ix2 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ix2 = periodic_ix2_particle;
#ifdef MPI_PARALLEL
	if(pG->lx2_id < 0 && pD->NGrid_x2 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx2_id = pD->GridArray[kb][pD->NGrid_x2-1][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: ibc_x2 = %d unknown\n",ibc_x2);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox2 == NULL){

      obc_x2 = par_geti("grid","obc_x2");
      switch(obc_x2){

      case 1: /* Reflecting */
	apply_ox2 = reflect_ox2_particle;
	break;
      case 5: /* Reflecting */
	apply_ox2 = reflect_ox2_particle;
	break;

      case 2: /* Outflow */
	apply_ox2 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ox2 = periodic_ox2_particle;
#ifdef MPI_PARALLEL
	if(pG->rx2_id < 0 && pD->NGrid_x2 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx2_id = pD->GridArray[kb][0][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: obc_x2 = %d unknown\n",obc_x2);
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
	apply_ix3 = reflect_ix3_particle;
	break;
      case 5: /* Reflecting */
	apply_ix3 = reflect_ix3_particle;
	break;

      case 2: /* Outflow */
	apply_ix3 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ix3 = periodic_ix3_particle;
#ifdef MPI_PARALLEL
	if(pG->lx3_id < 0 && pD->NGrid_x3 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->lx3_id = pD->GridArray[pD->NGrid_x3-1][jb][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: ibc_x3 = %d unknown\n",ibc_x3);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox3 == NULL){

      obc_x3 = par_geti("grid","obc_x3");
      switch(obc_x3){

      case 1: /* Reflecting */
	apply_ox3 = reflect_ox3_particle;
	break;
      case 5: /* Reflecting */
	apply_ox3 = reflect_ox3_particle;
	break;

      case 2: /* Outflow */
	apply_ox3 = outflow_particle;
	break;

      case 4: /* Periodic */
	apply_ox3 = periodic_ox3_particle;
#ifdef MPI_PARALLEL
	if(pG->rx3_id < 0 && pD->NGrid_x3 > 1){
	  get_myGridIndex(pD, my_id, &ib, &jb, &kb);
	  pG->rx3_id = pD->GridArray[0][jb][ib].id;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: obc_x3 = %d unknown\n",obc_x3);
	exit(EXIT_FAILURE);
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_fun:  sets function pointers for user-defined BCs in problem file
 */

void set_bvals_particle_fun(enum Direction dir, VBCFun_t prob_bc)
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

/* finalize boundary condition */
void set_bvals_particle_destruct(Grid *pG, Domain *pD)
{
  apply_ix1 = NULL;
  apply_ox1 = NULL;
  apply_ix2 = NULL;
  apply_ox2 = NULL;
  apply_ix3 = NULL;
  apply_ox3 = NULL;
  free(send_buf);
  free(recv_buf);
  send_bufsize = 0;
  recv_bufsize = 0;
  return;
}

/*----------------------------------------------------------------------------------------------*/
/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------------------------*/
/* Following are the functions:
 *   realloc_sendbuf & realloc_recvbuf
 *   reflecting_???_particle
 *   outflow_???_particle
 *   periodic_???_particle
 *   packing_???_particle
 *   unpack__particle
 *   shearing box related functions
 *   fargo related functions
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/* reallocate memory to send buffer */
static void realloc_sendbuf()
{
  send_bufsize += NBUF;
  if ((send_buf = (double*)realloc(send_buf, NVAR_P*(send_bufsize)*sizeof(double))) == NULL)
    ath_error("[set_bvals_prticles]: failed to allocate memory for buffer.\n");
  return;
}

/* reallocate memory to receive buffer */
static void realloc_recvbuf()
{
  recv_bufsize += NBUF;
  if ((recv_buf = (double*)realloc(recv_buf, NVAR_P*(recv_bufsize)*sizeof(double))) == NULL)
    ath_error("[set_bvals_prticles]: failed to allocate memory for buffer.\n");
  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1)
 */

static void reflect_ix1_particle(Grid *pG)
{
  Real x1b2;	/* 2*(x1 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x1b2 = 2.0*(pG->x1_0 + (pG->is + pG->idisp)*pG->dx1);
  /* pack boundary particles */
  n = packing_ix1_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x1 = x1b2 - pG->particle[p].x1;
    pG->particle[p].v1 = -pG->particle[p].v1;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (ibc_x1=1)
 */

static void reflect_ox1_particle(Grid *pG)
{
  Real x1b2;	/* 2*(x1 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x1b2 = 2.0*(pG->x1_0 + (pG->ie+1 + pG->idisp)*pG->dx1);
  /* pack boundary particles */
  n = packing_ox1_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x1 = x1b2 - pG->particle[p].x1;
    pG->particle[p].v1 = -pG->particle[p].v1;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 */

static void reflect_ix2_particle(Grid *pG)
{
  Real x2b2;	/* 2*(x2 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x2b2 = 2.0*(pG->x2_0 + (pG->js + pG->jdisp)*pG->dx2);
  /* pack boundary particles */
  n = packing_ix2_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x2 = x2b2 - pG->particle[p].x2;
    pG->particle[p].v2 = -pG->particle[p].v2;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (ibc_x2=1)
 */

static void reflect_ox2_particle(Grid *pG)
{
  Real x2b2;	/* 2*(x2 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x2b2 = 2.0*(pG->x2_0 + (pG->je+1 + pG->jdisp)*pG->dx2);
  /* pack boundary particles */
  n = packing_ox2_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x2 = x2b2 - pG->particle[p].x2;
    pG->particle[p].v2 = -pG->particle[p].v2;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 */

static void reflect_ix3_particle(Grid *pG)
{
  Real x3b2;	/* 2*(x3 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x3b2 = 2.0*(pG->x3_0 + (pG->ks + pG->kdisp)*pG->dx3);
  /* pack boundary particles */
  n = packing_ix3_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x3 = x3b2 - pG->particle[p].x3;
    pG->particle[p].v3 = -pG->particle[p].v3;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (ibc_x3=1)
 */

static void reflect_ox3_particle(Grid *pG)
{
  Real x3b2;	/* 2*(x3 border coordinate of grid cell pG->ie) */
  long n, n0, p;

  /* get lower and upper coordinate limit in x1 outer boundary */
  x3b2 = 2.0*(pG->x3_0 + (pG->ke+1 + pG->kdisp)*pG->dx3);
  /* pack boundary particles */
  n = packing_ox3_particle(pG, nghost);
  /* get the rear of the particle list */
  n0 = pG->nparticle;
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x3 = x3b2 - pG->particle[p].x3;
    pG->particle[p].v3 = -pG->particle[p].v3;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions
 * For particles, outflow B.C. = No B.C.  We only remove particles in the outermost
 * layer of the ghost cells, which is done in remove_ghost_particle(), see particle.c.
 */

static void outflow_particle(Grid *pG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 * Note: 2D shearing box B.C. is considered here!
 */

static void periodic_ix1_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ox1_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 1, -Lx1);
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
#ifndef FARGO
  /* velocity shift for shearing box */
  shift_packed_particle(send_buf, n, 6, vshear);
#endif
#endif /* SHEARING_BOX */
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Outer x1 boundary (ibc_x1=4)
 * Note: 2D shearing box B.C. is considered here!
 */

static void periodic_ox1_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ix1_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 1, Lx1);
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
#ifndef FARGO
  /* velocity shift for shearing box */
  shift_packed_particle(send_buf, n, 6, -vshear);
#endif
#endif /* SHEARING_BOX */
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ox2_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 2, -Lx2);
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Outer x2 boundary (ibc_x2=4)
 */

static void periodic_ox2_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ix2_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 2, Lx2);
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x3 boundary (ibc_x3=4)
 */

static void periodic_ix3_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ox3_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 3, -Lx3);
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Outer x3 boundary (ibc_x3=4)
 */

static void periodic_ox3_particle(Grid *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ix3_particle(pG, nghost);
  /* shift the particles */
  shift_packed_particle(send_buf, n, 3, Lx3);
  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/* packing the particle inside the inner x1 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ix1_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x1l,x1u;	/* lower and upper coordinate limit in x1 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x1l = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;
  x1u = pG->x1_0 + (pG->is+nlayer + pG->idisp)*pG->dx1;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x1 < x1u) && (cur->x1 >= x1l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }

  return n;
}

/* packing the particle inside the outer x1 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ox1_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x1l,x1u;	/* lower and upper coordinate limit in x1 outer boundary */
  long p, n = 0;
  double *pd = send_buf;
long i;
  /* get lower and upper coordinate limit in x1 inner boundary */
  x1l = pG->x1_0 + (pG->ie-nlayer+1 + pG->idisp)*pG->dx1;
  x1u = pG->x1_0 + (pG->ie+1 + pG->idisp)*pG->dx1;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x1 < x1u) && (cur->x1 >= x1l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }

  return n;
}

/* packing the particle inside the inner x2 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ix2_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x2l,x2u;	/* lower and upper coordinate limit in x2 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x2l = pG->x2_0 + (pG->js + pG->jdisp)*pG->dx2;
  x2u = pG->x2_0 + (pG->js+nlayer + pG->jdisp)*pG->dx2;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x2 < x2u) && (cur->x2 >= x2l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }

  return n;
}

/* packing the particle inside the outer x2 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ox2_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x2l,x2u;	/* lower and upper coordinate limit in x2 outer boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x2l = pG->x2_0 + (pG->je-nlayer+1 + pG->jdisp)*pG->dx2;
  x2u = pG->x2_0 + (pG->je+1 + pG->jdisp)*pG->dx2;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x2 < x2u) && (cur->x2 >= x2l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }

  return n;
}

/* packing the particle inside the inner x3 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ix3_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x3l,x3u;	/* lower and upper coordinate limit in x3 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x3l = pG->x3_0 + (pG->ks + pG->kdisp)*pG->dx3;
  x3u = pG->x3_0 + (pG->ks+nlayer + pG->kdisp)*pG->dx3;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x3 < x3u) && (cur->x3 >= x3l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }
  return n;
}

/* packing the particle inside the outer x3 boundary
   Input: pG: grid;
     nlayer: number of layers to be packed in the boundary
   Output:
     send_buf: buffer to save packed particle
     return: number of packed particle
*/
static long packing_ox3_particle(Grid *pG, int nlayer)
{
  Grain *cur;	/* current pointer */
  Real x3l,x3u;	/* lower and upper coordinate limit in x3 outer boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x3l = pG->x3_0 + (pG->ke-nlayer+1 + pG->kdisp)*pG->dx3;
  x3u = pG->x3_0 + (pG->ke+1 + pG->kdisp)*pG->dx3;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    cur = &(pG->particle[p]);
    if ((cur->x3 < x3u) && (cur->x3 >= x3l)) /* if the particle is the regime */
    { /* add one particle to buffer */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
    }
  }
  return n;
}

/* shift the coordinate/velocity of the packed particles by a constant amount
   Input:
     buf: buffer;	n: number of particles in the buffer
     index: 1: x1; 2: x2; 3: x3; 4: v1; 5: v2; 6: v3;
     shift: amount of change.
   Output:
     buf: buffer with shifted particles
*/
static void shift_packed_particle(double *buf, long n, int index, double shift)
{
  double *pd = buf;
  long i;

  if (n>0) pd += index-1;
  for (i=0; i<n-1; i++) {
    *pd += shift;
    pd += NVAR_P;
  }
  *pd += shift;

  return;
}

/* unpack received particle
   Input:
     pG: grid;
     buf: received buffer
     n: number of particle in the buffer
   Output:
     pG: grid with new particle added.
*/
static void unpack_particle(Grid *pG, double *buf, long n)
{
  Grain *cur;		/* current pointer */
  Grain *newgr;		/* space for new particle */
  double *pd = buf;
  long p, i;

  /* initialization */
  p = pG->nparticle;
  pG->nparticle += n;
  if (pG->nparticle >= pG->arrsize)
    particle_realloc(pG, pG->nparticle+1);
  /* unpacking */
  for (i=p; i<pG->nparticle; i++) {
    cur = &(pG->particle[i]);
    cur->x1 = *(pd++);
    cur->x2 = *(pd++);
    cur->x3 = *(pd++);
    cur->v1 = *(pd++);
    cur->v2 = *(pd++);
    cur->v3 = *(pd++);
    cur->property = (int)(*(pd++));
    pG->grproperty[cur->property].num += 1;
  }

  return;
}

#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
/* Shearing box boundary condition, Inner x1 boundary (ibc_x1=4)
 * This routine works for only 3D.
 * Input: numpar: for the packing routine, the array index to start with.
 */
static void shearingbox_ix1_particle(Grid *pG, Domain *pD, long numpar)
{
#ifdef MPI_PARALLEL
  Real yshear, yshift;		/* amount of shear, whole (yshear) and fractional (yshift) */
  Real x3len;			/* length in x3 (y) direction of a grid */
  int id1s, id1r, id2s, id2r;	/* processor ids to send and receive data */
  int ind1, ind2, ind3;		/* x1,x2,x3 indices of grid in the domain */
  long n1, n2;			/* number of packed particles */
  long cnt_recv;		/* received number of particles */
  int err;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */

#ifndef MPI_PARALLEL
  packing_ix1_particle_shear(pG, 0, numpar);
#else
/*-------------------------- MPI case: Step 1. Find locations ---------------------------------*/
  /* compute the distance the computational domain has sheared in y */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx3);
  x3len = pG->dx3*pG->Nx3;
  /* obtain processor ids to send and receive */
  ind1 = 0;		ind2 = my_jproc;
  ind3 = my_kproc + (int)(yshift/x3len) + 1;
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  id1s = pD->GridArray[ind3][ind2][ind1].id;	/* for region I (send to) */

  ind3 = my_kproc + (int)(yshift/x3len);
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  id2s = pD->GridArray[ind3][ind2][ind1].id;	/* for region II (send to) */

  ind3 = my_kproc - (int)(yshift/x3len) - 1;
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  id1r = pD->GridArray[ind3][ind2][ind1].id;	/* for region I (receive from) */

  ind3 = my_kproc - (int)(yshift/x3len);
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  id2r = pD->GridArray[ind3][ind2][ind1].id;	/* for region II (receive from) */

/*----------------------MPI case: Step 2. Exchange particles -------------------------------*/

  /* packing particles at region I */
  n1 = packing_ix1_particle_shear(pG, 1, numpar);

  /*-------- sub-step 1: exchange rigion I (id1) --------*/
  /* send and receive buffer size to/from region I (id1) */
  /* Post a non-blocking receive for the data size */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id1r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  /* send buffer size */
  err = MPI_Send(&n1, 1, MPI_LONG, id1s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  /* receive buffer size from */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf();

  /* send and receive buffer to/from region I (id1) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from outer region I */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id1r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }
  if (n1 > 0) {
    /* send buffer */
    err = MPI_Send(send_buf, n1*NVAR_P, MPI_DOUBLE, id1s, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  }
  /* packing particles at region II (N.B.!) */
  n2 = packing_ix1_particle_shear(pG, 2, numpar);
  if (cnt_recv > 0) {
    /* receive buffer */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);
    /* unpack the received particles */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

  /*-------- sub-step 2: exchange region II (id2)  --------*/
  /* send and receive buffer size to/from region II (id2) */
  /* Post a non-blocking receive for the data size */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id2r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  /* send buffer size */
  err = MPI_Send(&n2, 1, MPI_LONG, id2s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  /* receive buffer size */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf();

  /* send and receive buffer to/from region II (id2) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id2r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }
  if (n2 > 0) {
    /* send buffer */
    err = MPI_Send(send_buf, n2*NVAR_P, MPI_DOUBLE, id2s, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  }
  if (cnt_recv > 0) {
    /* receive buffer */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);
    /* unpack the received particles */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

#endif /* MPI_PARALLEL */

  return;
}

/* Shearing box boundary condition, Outer x1 boundary (obc_x1=4)
 * This routine works for only 3D.
 * Input: numpar: for the packing routine, the array index to start with.
 */
static void shearingbox_ox1_particle(Grid *pG, Domain *pD, long numpar)
{
#ifdef MPI_PARALLEL
  Real yshear, yshift;		/* amount of shear, whole (yshear) and fractional (yshift) */
  Real x3len;			/* length in x3 (y) direction of a grid */
  int id1s, id1r, id2s, id2r;	/* processor ids to send and receive data */
  int ind1, ind2, ind3;		/* x1,x2,x3 indices of grid in the domain */
  long n1, n2;			/* number of packed particles */
  long cnt_recv;		/* received number of particles */
  int err;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */

#ifndef MPI_PARALLEL
  packing_ox1_particle_shear(pG, 0, numpar);
#else
/*-------------------------- MPI case: Step 1. Find locations ---------------------------------*/
  /* compute the distance the computational domain has sheared in y */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx3);
  x3len = pG->dx3*pG->Nx3;
  /* obtain processor ids to send and receive */
  ind1 = pD->NGrid_x1-1;	ind2 = my_jproc;
  ind3 = my_kproc - (int)(yshift/x3len) - 1;
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  id1s = pD->GridArray[ind3][ind2][ind1].id;	/* for region I (send to) */

  ind3 = my_kproc - (int)(yshift/x3len);
  if (ind3 < 0) ind3 += pD->NGrid_x3;
  id2s = pD->GridArray[ind3][ind2][ind1].id;	/* for region II (send to) */

  ind3 = my_kproc + (int)(yshift/x3len) + 1;
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  id1r = pD->GridArray[ind3][ind2][ind1].id;	/* for region I (receive from) */

  ind3 = my_kproc + (int)(yshift/x3len);
  if (ind3 > (pD->NGrid_x3-1)) ind3 -= pD->NGrid_x3;
  id2r = pD->GridArray[ind3][ind2][ind1].id;	/* for region II (receive from) */

/*----------------------MPI case: Step 2. Exchange particles -------------------------------*/

  /* packing particles at region I */
  n1 = packing_ox1_particle_shear(pG, 1, numpar);

  /*-------- sub-step 1: outer region I with inner rigion I (id1) --------*/
  /* send and receive buffer size to/from inner region I (id1) */
  /* Post a non-blocking receive for the data size from inner regiion I */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id1r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  /* send buffer size to inner region I */
  err = MPI_Send(&n1, 1, MPI_LONG, id1s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  /* receive buffer size from the inner region I */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf();

  /* send and receive buffer to/from inner region I (id1) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from inner region I */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id1r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }
  if (n1 > 0) {
    /* send buffer to inner region I */
    err = MPI_Send(send_buf, n1*NVAR_P, MPI_DOUBLE, id1s, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  }
  /* packing particles at region II (N.B.!) */
  n2 = packing_ox1_particle_shear(pG, 2, numpar);
  if (cnt_recv > 0) {
    /* receive buffer inner region I */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);
    /* unpack the received particle */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

  /*-------- sub-step 2: outer region II with inner region II (id2)  --------*/
  /* send and receive buffer size to/from inner region II (id2) */
  /* Post a non-blocking receive for the data size from inner region II */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id2r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  /* send buffer size to inner region II */
  err = MPI_Send(&n2, 1, MPI_LONG, id2s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  /* receive buffer size from inner region II */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf();

  /* send and receive buffer to/from inner region II (id2) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from inner region II */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id2r, boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }
  if (n2 > 0) {
    /* send buffer to inner region II */
    err = MPI_Send(send_buf, n2*NVAR_P, MPI_DOUBLE, id2s, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  }
  if (cnt_recv > 0) {
    /* receive buffer from inner region II */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[receive_ox1_particle_shear]: MPI_Wait error = %d\n",err);
    /* unpack the received particle */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

#endif /* MPI_PARALLEL */

  return;
}

/* packing the particle inside the inner x1 boundary for shearing box
   Input: pG: grid; reg: region, 1 or 2 for mpi case, 0 for non-mpi case
          numpar: array index to start with.
   Return: number of packed particles in the specified region.
   Note: this routine serves for only 3D
*/
static long packing_ix1_particle_shear(Grid *pG, int reg, long numpar)
{
  Grain *cur;	/* current pointer */
  long n, p;
  Real ix1b;		/* coordinate limit in x1 inner boundary */
  Real yshear, yshift;	/* amount of shear, whole (yshear) and fractional (yshift) */
  Real x30, x3c;	/* x3c: y-coordinate marking the demarcation of the two regions */
  double *pd;

/*---------------- Step.1 -----------------------*/
  /* get the distance of shear */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx3);
  x30 = pG->x3_0+(pG->ke+pG->kdisp+1)*pG->dx3;
  x3c = x30 - fmod(yshear, pG->dx3*pG->Nx3);
  /* get coordinate limits for particles to be packed*/
  ix1b = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;
  /* buffer assignment */
  pd = send_buf;
  n = 0;

/*---------------- Step.2 -----------------------*/
  /* loop over all particle to pack particles in the boundary */
  p = numpar;
  while (p<pG->nparticle) {
    cur = &(pG->particle[p]);
    p += 1;
    if (cur->x1 < ix1b) { /* particles in inner boundary */
      if (((reg == 1) && (cur->x3 >= x3c)) || ((reg == 2) && (cur->x3 < x3c)))
      {         /* region I */                      /* region II */
        n += 1;
        if ((n+1) > send_bufsize) {
          realloc_sendbuf();
          pd = &send_buf[NVAR_P*(n-1)];
        }
        /* pack the particle and apply shift */
        *(pd++) = cur->x1;
        *(pd++) = cur->x2;
        *(pd++) = x3min + fmod(cur->x3 - x3min + yshift, Lx3);
        *(pd++) = cur->v1;
        *(pd++) = cur->v2;
        *(pd++) = cur->v3;
        *(pd++) = (double)(cur->property)+0.01;
        /* delete the particle */
        pG->nparticle -= 1;
        p -= 1;
        pG->particle[p] = pG->particle[pG->nparticle];
      }
      if (reg == 0)	/* non-mpi case, directly shift the particle positions */
        cur->x3 = x3min + fmod(cur->x3 - x3min + yshift, Lx3);
    }
  }

  return n;
}

/* packing the particle outside the outer x1 boundary for shearing box
   Input: pG: grid; reg: region, 1 or 2 for mpi case, 0 for non-mpi case
          numpar: array index to start with.
   Return: number of packed particles in the specified region.
   Note: this routine serves for only 3D
*/
static long packing_ox1_particle_shear(Grid *pG, int reg, long numpar)
{
  Grain *cur;	/* current pointer */
  long n, p;
  Real ox1b;		/* coordinate limit of x1 outer boundary */
  Real yshear, yshift;	/* amount of shear, whole (yshear) and fractional (yshift) */
  Real x30, x3c;	/* x3c: y-coordinate marking the demarcation of the two regions */
  double *pd;

/*---------------- Step.1 -----------------------*/
  /* get the distance of shear */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx3);
  x30 = pG->x3_0+(pG->ks+pG->kdisp)*pG->dx3;
  x3c = x30 + fmod(yshear, pG->dx3*pG->Nx3);
  /* get coordinate limits for particles to be packed*/
  ox1b = pG->x1_0 + (pG->ie + 1 + pG->idisp)*pG->dx1;
  /* buffer assignment */
  pd = send_buf;
  n = 0;

/*---------------- Step.2 -----------------------*/
  /* loop over all particle to pack particles in the boundary */
  p = numpar;
  while (p<pG->nparticle) {
    cur = &(pG->particle[p]);
    p += 1;
    if (cur->x1 >= ox1b) {/* particles in outer boundary */
      if (((reg == 1) && (cur->x3 < x3c)) || ((reg == 2) && (cur->x3 >= x3c)))
      {         /* region I */                      /* region II */
        n += 1;
        if ((n+1) > send_bufsize) {
          realloc_sendbuf();
          pd = &(send_buf[NVAR_P*(n-1)]);
	}
        /* pack the particle and apply shift */
        *(pd++) = cur->x1;
        *(pd++) = cur->x2;
        *(pd++) = x3min + fmod(cur->x3 - x3min + Lx3 - yshift, Lx3);
        *(pd++) = cur->v1;
        *(pd++) = cur->v2;
        *(pd++) = cur->v3;
        *(pd++) = (double)(cur->property)+0.01;
        /* delete the particle */
        pG->nparticle -= 1;
        p -= 1;
        pG->particle[p] = pG->particle[pG->nparticle];
      }
      if (reg == 0)	/* non-mpi case, directly shift the particle positions */
        cur->x3 = x3min + fmod(cur->x3 - x3min + Lx3 - yshift, Lx3);
    }
  }

  return n;
}

#ifdef FARGO
/* packing the particle for FARGO
   Input: pG: grid; 
      yl, yu: lower and upper limit of x3 (y) coordinate for particles to be packed
   Return: number of packed particles in the specified region.
*/
static long packing_particle_fargo(Grid *pG, Real yl, Real yu)
{
  Grain *cur;
  long p, n;
  double *pd;

  p = 0;
  n = 0;
  pd = send_buf;
  while (p<pG->nparticle) {
    cur = &(pG->particle[p]);
    p += 1;
    if ((cur->x3 >= yl) && (cur->x3 < yu))
    {/* particles in the range */
      n += 1;
      if ((n+1) > send_bufsize) {
        realloc_sendbuf();
        pd = &(send_buf[NVAR_P*(n-1)]);
      }
      /* pack the particle and apply shift */
      *(pd++) = cur->x1;
      *(pd++) = cur->x2;
      *(pd++) = cur->x3;
      *(pd++) = cur->v1;
      *(pd++) = cur->v2;
      *(pd++) = cur->v3;
      *(pd++) = (double)(cur->property)+0.01;
      /* delete the particle */
      pG->nparticle -= 1;
      p -= 1;
      pG->particle[p] = pG->particle[pG->nparticle];
    }
  }
  return n;
}

/* Get the number of shift of grid in the y direction */
static int gridshift(Real shift)
{
  if (shift>0)
    return (int)(shift)+1;
  else if (shift<0)
    return (int)(shift)-1;
  else return 0;
}
#endif /* FARGO */

#endif /* SHEARING_BOX */

#undef NBUF
#undef NVAR_P

#endif /*PARTICLES*/
