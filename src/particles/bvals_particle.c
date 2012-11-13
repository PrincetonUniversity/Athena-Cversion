#include "../copyright.h"
/*============================================================================*/
/*! \file bvals_particle.c
 *  \brief Sets boundary conditions for particles. 
 *
 * PURPOSE: Sets boundary conditions for particles. The basic algorithm is
 *   similar to the MHD boundary conditions. In setting B.C., particles are
 *   always packed to the send buffer (for both MPI and non-MPI cases), then
 *   unpacked with certain shifts in position or velocity. In the case of MPI,
 *   two communications are needed for sending and receiving particle, first on
 *   the number of particles to be sent/received, then to send/receive real
 *   data. The shearing box B.C. is first treated as periodic B.C., then move
 *   the particles in the ghost cells accordingly. Advection of particles when
 *   FARGO is turned on is also included.
 * 
 * CONTAINS PUBLIC FUNCTIONS:
 * - set_bvals_particle()
 * - advect_particles()
 * - set_bvals_particle_init()
 * - set_bvals_particle_fun()
 * - set_bvals_particle_destruct()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - realloc_???()            - reallocate send/recv buffer
 * - update_particle_status() - reset particle status (either ghost or grid)
 * - reflect_???_particle()   - apply reflecting BCs at boundary ???
 * - outflow_particle()       - apply outflow BCs at boundary ???
 * - periodic_???_particle()  - apply periodic BCs at boundary ???
 * - packing_???_particle()   - pack particles at boundary ???
 * - shift_packed_particle()  - shift packed particle (v or x) by certain amount
 * - unpack_particle()        - upack received particles
 * - shearingbox_???_particle()-shearing box BC at boundary ???
 * - packing_???_particle_shear()- pack particles for shearing box BC
 * - packing_particle_fargo() - pack particles for FARGO
 * - gridshift()              - calculate shift in grid in y for FARGO
 *
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"


#ifdef PARTICLES         /* endif at the end of the file */

/* particle structure size */
#ifdef MPI_PARALLEL
#define NVAR_P 10
#else
#define NVAR_P 9
#endif

/* send and receive buffer, size dynamically determined
 * They are mainly used for MPI, and shearing box.
 */
static double *send_buf = NULL, *recv_buf = NULL;
static long NBUF;	 /* buffer size unit (in number of particle) */
static long send_bufsize;/* size of the send buffer (in unit of particles) */
static long recv_bufsize;/* size of the recv buffer (in unit of particles) */
static int  nbc;	 /* number of boundary layers for particle BC */

/* processor indices in the computational domain */
static int my_iproc, my_jproc, my_kproc;
/* min and max coordinate limits of the computational domain */
static Real x1min,x1max,x2min,x2max,x3min,x3max;
static Real Lx1, Lx2, Lx3;/* domain size in x1, x2, x3 direction */
static Real TShuffle;	  /* number of time steps for resorting particles */

/* boundary condition function pointers. local to this function  */
static VGFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VGFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VGFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   realloc_???()            - reallocate send/recv buffer
 *   update_particle_status() - reset particle status (either ghost or grid)
 *   reflect_???_particle()   - apply reflecting BCs at boundary ???
 *   outflow_particle()       - apply outflow BCs at boundary ???
 *   periodic_???_particle()  - apply periodic BCs at boundary ???
 *   packing_???_particle()   - pack particles at boundary ???
 *   shift_packed_particle()  - shift packed particle (v or x) by certain amount
 *   unpack_particle()        - upack received particles
 *   shearingbox_???_particle()-shearing box BC at boundary ???
 *   packing_???_particle_shear()- pack particles for shearing box BC
 *   packing_particle_fargo() - pack particles for FARGO
 *   gridshift()              - calculate shift in grid in y for FARGO
 *============================================================================*/

static void realloc_sendbuf();
static void realloc_recvbuf(long newsize);

static void update_particle_status(GridS *pG);

static void reflect_ix1_particle(GridS *pG);
static void reflect_ox1_particle(GridS *pG);
static void reflect_ix2_particle(GridS *pG);
static void reflect_ox2_particle(GridS *pG);
static void reflect_ix3_particle(GridS *pG);
static void reflect_ox3_particle(GridS *pG);

static void outflow_particle(GridS *pG);

static void periodic_ix1_particle(GridS *pG);
static void periodic_ox1_particle(GridS *pG);
static void periodic_ix2_particle(GridS *pG);
static void periodic_ox2_particle(GridS *pG);
static void periodic_ix3_particle(GridS *pG);
static void periodic_ox3_particle(GridS *pG);

static long packing_ix1_particle(GridS *pG, int nlayer);
static long packing_ox1_particle(GridS *pG, int nlayer);
static long packing_ix2_particle(GridS *pG, int nlayer);
static long packing_ox2_particle(GridS *pG, int nlayer);
static long packing_ix3_particle(GridS *pG, int nlayer);
static long packing_ox3_particle(GridS *pG, int nlayer);
static void packing_one_particle(GrainS *gr, long n, short pos);
static void shift_packed_particle(double *buf, long n, int index, double shift);
static void unpack_particle(GridS *pG, double *buf, long n);

#ifdef SHEARING_BOX
static void shearingbox_ix1_particle(GridS *pG, DomainS *pD, long numpar);
static void shearingbox_ox1_particle(GridS *pG, DomainS *pD, long numpar);

static long packing_ix1_particle_shear(GridS *pG, int reg, long numpar);
static long packing_ox1_particle_shear(GridS *pG, int reg, long numpar);

#ifdef FARGO
static long packing_particle_fargo(GridS *pG, Real yl, Real yu);
static int gridshift(Real shift);
#endif /* FARGO */

#endif /* SHEARING_BOX */

extern void Delete_Ghost(GridS *pG);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void bvals_particle(GridS *pG, Domain *pD)
 *  \brief Calls appropriate functions to set particle BCs. 
 *
 *   The
 *   function pointers (*apply_???) are set during initialization by
 *   set_bvals_particle_init() to be either a user-defined function, or one of
 *   the functions corresponding to reflecting, periodic, or outflow.  If the
 *   left- or right-Grid ID numbers are >= 1 (neighboring grids exist), then
 *   MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void bvals_particle(DomainS *pD)
{
  GridS *pG = pD->Grid;
#ifdef MPI_PARALLEL
  int err;
  long cnt_send, cnt_recv;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */
#ifdef SHEARING_BOX
  long numpar;  /* number of particles before applying B.C. */
  int vyind;	/* index for vy, 5 for 3D (x1,x2,x3)=(X,Y,Z),
                                 6 for 2D (x1,x2,x3)=(X,Z,Y) */
#endif

/*--- Step 1. ------------------------------------------------------------------
 * shuffle if necessary */

  /* shuffle every time interval TShuffle */
  /* if TShuffle is not positive, don't shuffle */
  if ((TShuffle>0) && (fmod(pG->time, TShuffle)<pG->dt))
    shuffle(pG);

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pG->Nx[0] > 1){

#ifdef SHEARING_BOX
  numpar = pG->nparticle;
  if (ShBoxCoord == xy) /* (x1,x2,x3)=(X,Y,Z) */
    vyind = 5;
  else                  /* (x1,x2,x3)=(X,Z,Y) */
    vyind = 6;
#endif

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx1_id >= 0 && pG->lx1_id >= 0) {

    /*-------- sub-step 1: send to right and receive from left ------------*/

      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox1_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx1_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (my_iproc == pD->NGrid[0]-1) {
        /* physical boundary on the rignt in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 1, -Lx1);
#ifdef SHEARING_BOX
#ifndef FARGO
        /* velocity shift for shearing box */
        shift_packed_particle(send_buf, cnt_send, vyind, vshear);
#endif
#endif /* SHEARING_BOX */
      }

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx1_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix1_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx1_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (my_iproc == 0) {
        /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 1, Lx1);
#ifdef SHEARING_BOX
#ifndef FARGO
        /* velocity shift for shearing box */
        shift_packed_particle(send_buf, cnt_send, vyind, -vshear);
#endif
#endif /* SHEARING_BOX */
      }

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx1_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      cnt_send = packing_ox1_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx1_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      }

    /*-------- sub-step 2: apply left boundary condition --------*/
      (*apply_ix1)(pG);

    /*-------- sub-step 3: receive from right  --------*/

      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

     /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx1_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix1_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx1_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

    /*-------- sub-step 2: apply right boundary condition --------*/
      (*apply_ox1)(pG);

    /*-------- sub-step 3: send to left --------*/

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix1_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx1_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx1_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
   if (ShBoxCoord == xy) {
      if (my_iproc == 0) /* inner boundary */
        shearingbox_ix1_particle(pG, pD, numpar);

      if (my_iproc == (pD->NGrid[0]-1)) /* outer boundary */
        shearingbox_ox1_particle(pG, pD, numpar);
    }
#endif /* SHEARING_BOX */

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pG->Nx[1] > 1) {

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx2_id >= 0 && pG->lx2_id >= 0) {

    /*-------- sub-step 1: send to right and receive from left --------*/

      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox2_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx2_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      /* physical boundary on the rignt in periodic B.C. */
      if (my_jproc == pD->NGrid[1]-1)
        shift_packed_particle(send_buf, cnt_send, 2, -Lx2);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx2_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix2_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx2_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      if (my_jproc == 0) /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 2, Lx2);

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx2_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      cnt_send = packing_ox2_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx2_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);

      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix1_particle]: MPI_Send error = %d\n",err);
      }

    /*-------- sub-step 2: apply left boundary condition --------*/

      (*apply_ix2)(pG);

    /*-------- sub-step 3: receive from right --------*/

      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx2_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix2_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx2_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

    /*-------- sub-step 2: apply right boundary condition --------*/

      (*apply_ox2)(pG);

    /*-------- sub-step 3: send to left --------*/

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix2_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx2_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix2_particle]: MPI_Send error = %d\n",err);

      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx2_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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

  if (pG->Nx[2] > 1){

#ifdef MPI_PARALLEL

/* MPI blocks to both left and right */
    if (pG->rx3_id >= 0 && pG->lx3_id >= 0) {

    /*-------- sub-step 1: send to right and receive from left --------*/

      /* Post a non-blocking receive for the data size from the left grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the right to send buffer */
      cnt_send = packing_ox3_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx3_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      /* physical boundary on the rignt in periodic B.C. */
      if (my_kproc == pD->NGrid[2]-1) 
        shift_packed_particle(send_buf, cnt_send, 3, -Lx3);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the right and obtain recv_buf from the left */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx3_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix3_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx3_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      if (my_kproc == 0) /* physical boundary on the left in periodic B.C. */
        shift_packed_particle(send_buf, cnt_send, 3, Lx3);

      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      /* send send_buf to the left and obtain recv_buf from the right */
      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx3_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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
      cnt_send = packing_ox3_particle(pG, nbc);

      /* send buffer size to the right grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->rx3_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      /* send buffer to the right grid */
      if (cnt_send > 0) {
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->rx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
        if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);
      }

    /*-------- sub-step 2: apply left boundary condition --------*/

      (*apply_ix3)(pG);

    /*-------- sub-step 3: receive from right --------*/

      /* Post a non-blocking receive for the data size from the right grid */
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->rx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      /* receive buffer size from the right grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the right grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->rx3_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
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
      err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, pG->lx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

      /* receive buffer size from the left grid */
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[receive_ix3_particle]: MPI_Wait error = %d\n",err);

      /* check space for the receive buffer */
      if ((cnt_recv+1) >= recv_bufsize)
        realloc_recvbuf(cnt_recv+1);

      if (cnt_recv > 0) {
        /* Post a non-blocking receive for the input data from the left grid */
        err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, pG->lx3_id,
                                  boundary_particle_tag, MPI_COMM_WORLD, &rq);
        if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
      }

    /*-------- sub-step 2: apply right boundary condition --------*/

      (*apply_ox3)(pG);

    /*-------- sub-step 3: send to left --------*/

      /* packing particle on the left to send buffer */
      cnt_send = packing_ix3_particle(pG, nbc);

      /* send buffer size to the left grid */
      err = MPI_Send(&cnt_send, 1, MPI_LONG, pG->lx3_id,
                                boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ix3_particle]: MPI_Send error = %d\n",err);

      if (cnt_send > 0) {
        /* send buffer to the right grid */
        err = MPI_Send(send_buf, cnt_send*NVAR_P, MPI_DOUBLE, pG->lx3_id,
                                 boundary_particle_tag, MPI_COMM_WORLD);
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

/*--- Step 5. ------------------------------------------------------------------
 * Update the status of the crossing particles */
  update_particle_status(pG);

  Delete_Ghost(pG);

  return;
}


#ifdef FARGO
/*----------------------------------------------------------------------------*/
/*! \fn void advect_particles(DomainS *pD)
 *  \brief Advect particles by qshear*Omega_0*x1*dt for the FARGO algorithm. */
void advect_particles(DomainS *pD)
{
  GridS *pG = pD->Grid;
  GrainS *gr;
  long p;
  Real x1l, x1u;
#ifdef MPI_PARALLEL
  long cnt_recv, n;
  int ishl, ishu, i;
  int inds, indr, ids, idr;
  Real x2len, yl, yu;
  int err;
  MPI_Request rq;
  MPI_Status stat;
#endif /* MPI_PARALLEL */

  /* Do nothing if the azimuthal dimension is not present */
  if (ShBoxCoord != xy)
    return;

  /* lower and upper bound of the grid in x1 direction */
  x1l = pG->MinX[0];
  x1u = pG->MaxX[0];

  /* shift the particles */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);
    gr->x2 = x2min + fmod(gr->x2 + pG->parsub[p].shift - x2min + Lx2, Lx2);
  }

#ifdef MPI_PARALLEL
  /* calculate the farthest grid that advection can reach */
  x2len = pG->dx2*pG->Nx[1];
  ishl = gridshift((qshear*Omega_0*(x1l-pG->dx1)*pG->dt - pG->dx2)/x2len);
  ishu = gridshift((qshear*Omega_0*(x1u+pG->dx1)*pG->dt + pG->dx2)/x2len);

  /* loop over all the possible destination grids */
  for (i=ishl; i<=ishu; i++)
  if (i != 0) { /* avoid moving particles to the same grid */

    /* find the processor id to send/receive data */
    inds = my_jproc + i;
    if (inds < 0) inds += pD->NGrid[1];
    if (inds > (pD->NGrid[1]-1)) inds -= pD->NGrid[1];
    ids = pD->GData[my_kproc][inds][my_iproc].ID_Comm_Domain; /* send to */

    indr = my_jproc - i;
    if (indr < 0) indr += pD->NGrid[1];
    if (indr > (pD->NGrid[1]-1)) indr -= pD->NGrid[1];
    idr = pD->GData[my_kproc][indr][my_iproc].ID_Comm_Domain;/* receive from */

    /* Post a non-blocking receive for the data size */
    err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, idr, boundary_particle_tag,
                                                 MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

    /* packing particles */
    yl = x2min + inds*x2len;	yu = x2min + (inds+1)*x2len;
    n = packing_particle_fargo(pG, yl, yu);

    /* send buffer size */
    err = MPI_Send(&n, 1, MPI_LONG, ids, boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);

    /* receive buffer size */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[recv_ox1_particle_shear]: MPI_Wait error = %d\n",err);

    /* check space for the receive buffer */
    if ((cnt_recv+1) >= recv_bufsize)
      realloc_recvbuf(cnt_recv+1);

    /* Post a non-blocking receive for data */
    if (cnt_recv > 0) {
      err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, idr,
                                boundary_particle_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
    }

    /* send buffer */
    if (n > 0) {
      err = MPI_Send(send_buf, n*NVAR_P, MPI_DOUBLE, ids,
                               boundary_particle_tag, MPI_COMM_WORLD);
      if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
    }

    /* receive buffer */
    if (cnt_recv > 0) {
      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[recv_ox1_particle_shear]: MPI_Wait error = %d\n",err);

      /* unpack the received particle */
      unpack_particle(pG, recv_buf, cnt_recv);
    }
  }
#endif /* MPI_PARALLEL */

  return;
}
#endif /* FARGO */

/*----------------------------------------------------------------------------*/
/*! \fn void bvals_particle_init(MeshS *pM)
 *  \brief Sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */
void bvals_particle_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;

  if (pM->NLevels > 1)
    ath_error("[bval_particle_init]: particle module does not support SMR\n");

  pD = &(pM->Domain[0][0]);
  pG = pD->Grid;

#ifdef MPI_PARALLEL
  int ib,jb,kb;
#endif /* MPI_PARALLEL */

/* initialize buffers */
  NBUF = (long)(0.15*pG->arrsize);

  send_bufsize = NBUF;
  recv_bufsize = NBUF;
  send_buf = (double*)calloc_1d_array(NVAR_P*send_bufsize, sizeof(double));
  recv_buf = (double*)calloc_1d_array(NVAR_P*recv_bufsize, sizeof(double));

/* number of boundary layers to pack the particles */
#ifdef FEEDBACK
  nbc = 0;  /* need 4 layers of ghost particles for feedback predictor */
#else
  nbc = 0;  /* leave one layer for output purposes */
#endif

/* calculate distances of the computational domain and shear velocity */
  x1min = pD->RootMinX[0];
  x1max = pD->RootMaxX[0];
  Lx1 = x1max - x1min;

  x2min = pD->RootMinX[1];
  x2max = pD->RootMaxX[1];
  Lx2 = x2max - x2min;

  x3min = pD->RootMinX[2];
  x3max = pD->RootMaxX[2];
  Lx3 = x3max - x3min;

  get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

  /* get the number of time steps for shuffle */
  TShuffle = par_getd_def("particle","tshuf",0.0);/* by default, not shuffle */

#ifdef SHEARING_BOX
  /* shear velocity between inner and outer x1 boundaries */
  vshear = qshear * Omega_0 * Lx1;
#endif /* SHEARING_BOX */

/* Set function pointers for physical boundaries in x1-direction */

  if(pG->Nx[0] > 1) {
    if(apply_ix1 == NULL){

      switch(pM->BCFlag_ix1){

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
	if(pG->lx1_id < 0 && pD->NGrid[0] > 1){
	  pG->lx1_id =
          pD->GData[my_kproc][my_jproc][pD->NGrid[0]-1].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ix1 = %d unknown\n",
                    pM->BCFlag_ix1);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox1 == NULL){

      switch(pM->BCFlag_ox1){

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
	if(pG->rx1_id < 0 && pD->NGrid[0] > 1){
	  pG->rx1_id = pD->GData[my_kproc][my_jproc][0].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ox1 = %d unknown\n",
                    pM->BCFlag_ox1);
	exit(EXIT_FAILURE);
      }

    }
  }

/* Set function pointers for physical boundaries in x2-direction */

  if(pG->Nx[1] > 1) {
    if(apply_ix2 == NULL){

      switch(pM->BCFlag_ix2){

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
	if(pG->lx2_id < 0 && pD->NGrid[1] > 1){
	  pG->lx2_id =
          pD->GData[my_kproc][pD->NGrid[1]-1][my_iproc].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ix2 = %d unknown\n",
                    pM->BCFlag_ix2);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox2 == NULL){

      switch(pM->BCFlag_ox2){

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
	if(pG->rx2_id < 0 && pD->NGrid[1] > 1){
	  pG->rx2_id = pD->GData[my_kproc][0][my_iproc].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ox2 = %d unknown\n",
                                                         pM->BCFlag_ox2);
	exit(EXIT_FAILURE);
      }

    }
  }

/* Set function pointers for physical boundaries in x3-direction */

  if(pG->Nx[2] > 1) {
    if(apply_ix3 == NULL){

      switch(pM->BCFlag_ix3){

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
	if(pG->lx3_id < 0 && pD->NGrid[2] > 1){
	  pG->lx3_id =
          pD->GData[pD->NGrid[2]-1][my_jproc][my_iproc].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ix3 = %d unknown\n",
                    pM->BCFlag_ix3);
	exit(EXIT_FAILURE);
      }

    }

    if(apply_ox3 == NULL){

      switch(pM->BCFlag_ox3){

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
	if(pG->rx3_id < 0 && pD->NGrid[2] > 1){
	  pG->rx3_id = pD->GData[0][my_jproc][my_iproc].ID_Comm_Domain;
	}
#endif /* MPI_PARALLEL */
	break;

      default:
	ath_perr(-1,"[set_bvals_particle_init]: bc_ox3 = %d unknown\n",
                    pM->BCFlag_ox3);
	exit(EXIT_FAILURE);
      }

    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void set_bvals_particle_fun(enum BCDirection dir, VBCFun_t prob_bc)
 *  \brief Sets function pointers for user-defined BCs in problem file
 */

void set_bvals_particle_fun(enum BCDirection dir, VGFun_t prob_bc)
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

/*! \fn void bvals_particle_destruct(MeshS *pM)
 *  \brief Finalize boundary condition */
void bvals_particle_destruct(MeshS *pM)
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

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/* Following are the functions:
 *   realloc_sendbuf & realloc_recvbuf
 *   update_particle_status
 *   reflecting_???_particle
 *   outflow_???_particle
 *   periodic_???_particle
 *   packing_???_particle
 *   unpack__particle
 *   shearing box related functions
 *   fargo related functions
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 */

/*----------------------------------------------------------------------------*/
/*! \fn static void realloc_sendbuf()
 *  \brief Reallocate memory to send buffer */
static void realloc_sendbuf()
{
  send_bufsize += NBUF;
  ath_pout(1,"[set_bvals_prticles]: reallocating send buffer...");
  if ((send_buf = (double*)realloc(send_buf,
                           NVAR_P*(send_bufsize)*sizeof(double))) == NULL)
    ath_error("[set_bvals_prticles]: failed to allocate memory for buffer.\n");

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void realloc_recvbuf(long newsize)
 *  \brief Reallocate memory to receive buffer */
static void realloc_recvbuf(long newsize)
{
  recv_bufsize += NBUF;
  recv_bufsize = MAX(recv_bufsize, newsize);
  ath_pout(1,"[set_bvals_prticles]: reallocating receive buffer...");
  if ((recv_buf = (double*)realloc(recv_buf,
                           NVAR_P*(recv_bufsize)*sizeof(double))) == NULL)
    ath_error("[set_bvals_prticles]: failed to allocate memory for buffer.\n");

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_particle_status(GridS *pG)
 *  \brief Update the status of the particles after applying boundary conditions
 */
static void update_particle_status(GridS *pG)
{
  long p;
  GrainS *gr;

  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);
    if (gr->pos >= 10) /* crossing out/in particle from the previous step */
    {
      if ((gr->x1>=x1upar) || (gr->x1< x1lpar) || (gr->x2>=x2upar) ||
          (gr->x2< x2lpar) || (gr->x3>=x3upar) || (gr->x3< x3lpar))
        gr->pos = 0; /* ghost particle */

      else
        gr->pos = 1; /* grid particle */
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1) */
static void reflect_ix1_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ix1_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x1 = 2.0*pG->MinX[0] - pG->particle[p].x1;
    pG->particle[p].v1 = -pG->particle[p].v1;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (ibc_x1=1)
 */
static void reflect_ox1_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ox1_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x1 = 2.0*pG->MaxX[0] - pG->particle[p].x1;
    pG->particle[p].v1 = -pG->particle[p].v1;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1)
 */
static void reflect_ix2_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ix2_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x2 = 2.0*pG->MinX[1] - pG->particle[p].x2;
    pG->particle[p].v2 = -pG->particle[p].v2;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (ibc_x2=1)
 */
static void reflect_ox2_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ox2_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x2 = 2.0*pG->MaxX[1] - pG->particle[p].x2;
    pG->particle[p].v2 = -pG->particle[p].v2;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix3_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1)
 */
static void reflect_ix3_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ix3_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x3 = pG->MinX[2] - pG->particle[p].x3;
    pG->particle[p].v3 = -pG->particle[p].v3;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3_particle(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (ibc_x3=1)
 */
static void reflect_ox3_particle(GridS *pG)
{
  long n, n0, p;

  /* pack boundary particles */
  n = packing_ox3_particle(pG, nbc);

  /* get the rear of the particle list */
  n0 = pG->nparticle;

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);


  /* apply reflection boundary condition */
  for (p=n0; p<pG->nparticle; p++)
  {
    pG->particle[p].x3 = 2.0*pG->MaxX[2] - pG->particle[p].x3;
    pG->particle[p].v3 = -pG->particle[p].v3;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_particle(GridS *pG)
 *  \brief OUTFLOW boundary conditions
 *
 * For particles, outflow B.C. = No B.C.  We only remove particles in the
 * outermost layer of the ghost cells, which is done in remove_ghost_particle(),
 * see particle.c.
 */
static void outflow_particle(GridS *pG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 *
 * Note: 2D shearing box B.C. is considered here!
 */
static void periodic_ix1_particle(GridS *pG)
{
  long n = 0;
#ifdef SHEARING_BOX
  /* index for vy, 5 for 3D (x1,x2,x3)=(X,Y,Z), 6 for 2D (x1,x2,x3)=(X,Z,Y) */
  int vyind;

  if (ShBoxCoord == xy) /* (x1,x2,x3)=(X,Y,Z) */
    vyind = 5;
  else                  /* (x1,x2,x3)=(X,Z,Y) */
    vyind = 6;
#endif

  /* pack boundary particles */
  n = packing_ox1_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 1, -Lx1);

#ifdef SHEARING_BOX
#ifndef FARGO
  /* velocity shift for shearing box */
  shift_packed_particle(send_buf, n, vyind, vshear);
#endif
#endif /* SHEARING_BOX */

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Outer x1 boundary (ibc_x1=4)
 *
 * Note: 2D shearing box B.C. is considered here!
 */
static void periodic_ox1_particle(GridS *pG)
{
  long n = 0;
#ifdef SHEARING_BOX
  /* index for vy, 5 for 3D (x1,x2,x3)=(X,Y,Z), 6 for 2D (x1,x2,x3)=(X,Z,Y) */
  int vyind;

  if (ShBoxCoord == xy) /* (x1,x2,x3)=(X,Y,Z) */
    vyind = 5;
  else             /* (x1,x2,x3)=(X,Z,Y) */
    vyind = 6;
#endif

  /* pack boundary particles */
  n = packing_ix1_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 1, Lx1);

#ifdef SHEARING_BOX
#ifndef FARGO
  /* velocity shift for shearing box */
  shift_packed_particle(send_buf, n, vyind, -vshear);
#endif
#endif /* SHEARING_BOX */

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x2 boundary (ibc_x2=4)
 */
static void periodic_ix2_particle(GridS *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ox2_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 2, -Lx2);

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Outer x2 boundary (ibc_x2=4)
 */
static void periodic_ox2_particle(GridS *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ix2_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 2, Lx2);

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x3 boundary (ibc_x3=4)
 */
static void periodic_ix3_particle(GridS *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ox3_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 3, -Lx3);

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3_particle(GridS *pG)
 *  \brief PERIODIC boundary conditions, Outer x3 boundary (ibc_x3=4)
 */
static void periodic_ox3_particle(GridS *pG)
{
  long n = 0;

  /* pack boundary particles */
  n = packing_ix3_particle(pG, nbc);

  /* shift the particles */
  shift_packed_particle(send_buf, n, 3, Lx3);

  /* copy boudary particles to the rear */
  unpack_particle(pG, send_buf, n);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ix1_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the inner x1 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ix1_particle(GridS *pG, int nlayer)
{
  GrainS *gr;
  Real x1l,x1u;	/* lower and upper coordinate limit in x1 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x1l = pG->MinX[0];
  x1u = pG->MinX[0] + nlayer*pG->dx1;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x1 < x1u) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x1 >= x1l) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }

      else if (gr->pos != 21) /* it's not from ox1 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ix1 */
        packing_one_particle(gr, n, 11);

        n += 1;
      }
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ox1_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the outer x1 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ox1_particle(GridS *pG, int nlayer)
{
  GrainS *gr;
  Real x1l,x1u;	/* lower and upper coordinate limit in x1 outer boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x1l = pG->MaxX[0] - nlayer*pG->dx1;
  x1u = pG->MaxX[0];

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x1 >= x1l) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x1 < x1u) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }
      else if (gr->pos != 11) /* it's not from ix1 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ox1 */
        packing_one_particle(gr, n, 21);
        n += 1;
      }
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ix2_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the inner x2 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ix2_particle(GridS *pG, int nlayer)
{
  GrainS *gr;
  Real x2l,x2u;	/* lower and upper coordinate limit in x2 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x2l = pG->MinX[1];
  x2u = pG->MinX[1] + nlayer*pG->dx2;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x2 < x2u) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x2 >= x2l) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }
      else if (gr->pos != 22) /* it's not from ox2 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ox1 */
        packing_one_particle(gr, n, 12); 
        n += 1;
      }
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ox2_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the outer x2 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ox2_particle(GridS *pG, int nlayer)
{
  GrainS *gr;
  Real x2l,x2u;	/* lower and upper coordinate limit in x2 outer boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x2l = pG->MaxX[1] - nlayer*pG->dx2;
  x2u = pG->MaxX[1];

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x2 >= x2l) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x2 < x2u) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }
      else if (gr->pos != 12) /* it's not from ix2 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ox2 */
        packing_one_particle(gr, n, 22);
        n += 1;
      }
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ix3_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the inner x3 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ix3_particle(GridS *pG, int nlayer)
{
  GrainS *gr;	/* current pointer */
  Real x3l,x3u;	/* lower and upper coordinate limit in x3 inner boundary */
  long p, n = 0;
  double *pd = send_buf;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x3l = pG->MinX[2];
  x3u = pG->MinX[2] + nlayer*pG->dx3;

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x3 < x3u) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x3 >= x3l) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }
      else if (gr->pos != 23) /* it's not from ox3 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ix3 */
        packing_one_particle(gr, n, 13); 
        n += 1;
      }
    }
  }
  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ox3_particle(GridS *pG, int nlayer)
 *  \brief Packing the particle inside the outer x3 boundary
 *
 * Input: pG: grid;
 *   nlayer: number of layers to be packed in the boundary
 * Output:
 *   send_buf: buffer to save packed particle
 *   return: number of packed particle
 */
static long packing_ox3_particle(GridS *pG, int nlayer)
{
  GrainS *gr;	/* current pointer */
  Real x3l,x3u;	/* lower and upper coordinate limit in x3 outer boundary */
  long p, n = 0;

  /* get lower and upper coordinate limit in x1 inner boundary */
  x3l = pG->MaxX[2] - nlayer*pG->dx3;
  x3u = pG->MaxX[2];

  /* loop over all particle to pack the ones in the boundary */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    if (gr->x3 >= x3l) {
      if ((gr->pos == 0) || (gr->pos == 1))
      { /* ghost particle or grid particle */
        if (gr->x3 < x3u) {/* in the boundary */
          packing_one_particle(gr, n, 0); /* pack as ghost particle */
          n += 1;
        }
      }
      else if (gr->pos != 13) /* it's not from ix3 */
      {/* crossing particle in the boundary */
        /* pack as crossing particle from ox3 */
        packing_one_particle(gr, n, 23);
        n += 1;
      }
    }
  }
  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void packing_one_particle(GrainS *cur, long n, short pos)
 *  \brief Subroutine for packing one particle to send buffer
 *
 * Input:
 *   gr: particle pointer;
 *   p: starting index in the buffer
 *   pos: particle position (0: ghost; 1: grid; 2: cross in/out;
 * Output:
 *   one particle is added to the send buffer
 */
static void packing_one_particle(GrainS *gr, long n, short pos)
{
  double *pd;
  if ((n+2) > send_bufsize) {
    realloc_sendbuf();
  }
  pd = &(send_buf[NVAR_P*n]);

  /* pack the particle */
  *(pd++) = gr->x1;
  *(pd++) = gr->x2;
  *(pd++) = gr->x3;
  *(pd++) = gr->v1;
  *(pd++) = gr->v2;
  *(pd++) = gr->v3;
  *(pd++) = (double)(gr->property)+0.01;
  *(pd++) = (double)(pos)+0.01;
  *(pd++) = (double)(gr->my_id)+0.01;
#ifdef MPI_PARALLEL
  *(pd++) = (double)(gr->init_id)+0.01;
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void shift_packed_particle(double *buf, long n, int index, 
 *					  double shift)
 *  \brief shift the coordinate/velocity of the packed particles by a constant 
 *  amount
 *
 * Input:
 *   buf: buffer;	n: number of particles in the buffer
 *   index: 1: x1; 2: x2; 3: x3; 4: v1; 5: v2; 6: v3;
 *   shift: amount of change.
 * Output:
 *   buf: buffer with shifted particles
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

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_particle(GridS *pG, double *buf, long n)
 *  \brief Unpack received particle
 * Input:
 *   pG: grid;
 *   buf: received buffer
 *   n: number of particle in the buffer
 * Output:
 *   pG: grid with new particle added.
 */
static void unpack_particle(GridS *pG, double *buf, long n)
{
  GrainS *gr;		/* current pointer */
  double *pd = buf;
  long p, i;

  /* initialization */
  p = pG->nparticle;
  pG->nparticle += n;
  if (pG->nparticle >= pG->arrsize)
    particle_realloc(pG, pG->nparticle+1);

  /* unpacking */
  for (i=p; i<pG->nparticle; i++) {
    gr = &(pG->particle[i]);
    gr->x1 = *(pd++);
    gr->x2 = *(pd++);
    gr->x3 = *(pd++);
    gr->v1 = *(pd++);
    gr->v2 = *(pd++);
    gr->v3 = *(pd++);
    gr->property = (int)(*(pd++));
    grproperty[gr->property].num += 1;
    gr->pos = (short)(*(pd++));
    gr->my_id = (long)(*(pd++));
#ifdef MPI_PARALLEL
    gr->init_id = (int)(*(pd++));
#endif
  }

  return;
}

#ifdef SHEARING_BOX
/*----------------------------------------------------------------------------*/
/*! \fn static void shearingbox_ix1_particle(GridS *pG, DomainS *pD,long numpar)
 *  \brief Shearing box boundary condition, Inner x1 boundary (ibc_x1=4)
 *
 * This routine works only for ShBoxCoord = xy.
 * Input: numpar: for the packing routine, the array index to start with.
 */
static void shearingbox_ix1_particle(GridS *pG, DomainS *pD, long numpar)
{
#ifdef MPI_PARALLEL
  /* amount of shear, whole (yshear) and fractional (yshift) */
  Real yshear, yshift;
  Real x2len;			/* length in x2 (y) direction of a grid */
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
/*----------------- MPI case: Step 1. Find locations -------------------------*/

  /* compute the distance the computational domain has sheared in y */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx2);
  x2len = pG->dx2*pG->Nx[1];

  /* obtain processor ids to send and receive */
  ind1 = 0;		ind3 = my_kproc;
  ind2 = my_jproc + (int)(yshift/x2len) + 1;
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  id1s = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region I (send to) */

  ind2 = my_jproc + (int)(yshift/x2len);
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  id2s = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region II(send to) */

  ind2 = my_jproc - (int)(yshift/x2len) - 1;
  if (ind2 < 0) ind2 += pD->NGrid[1];
  if (ind2 < 0) ind2 += pD->NGrid[1];
  id1r = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region I (rcv from)*/

  ind2 = my_jproc - (int)(yshift/x2len);
  if (ind2 < 0) ind2 += pD->NGrid[1];
  id2r = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region II(rcv from)*/

/*------------------MPI case: Step 2. Exchange particles ---------------------*/

  /* packing particles at region I */
  n1 = packing_ix1_particle_shear(pG, 1, numpar);

/*-------- sub-step 1: exchange rigion I (id1) --------*/

  /* send and receive buffer size to/from region I (id1) */
  /* Post a non-blocking receive for the data size */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id1r,
                             boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

  /* send buffer size */
  err = MPI_Send(&n1, 1, MPI_LONG, id1s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);

  /* receive buffer size from */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf(cnt_recv+1);

  /* send and receive buffer to/from region I (id1) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from outer region I */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id1r,
                              boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }

  if (n1 > 0) {
    /* send buffer */
    err = MPI_Send(send_buf, n1*NVAR_P, MPI_DOUBLE, id1s,
                             boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  }

  /* packing particles at region II (N.B.!) */
  n2 = packing_ix1_particle_shear(pG, 2, numpar);
  if (cnt_recv > 0) {
    /* receive buffer */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[recv_ix1_particle_shear]: MPI_Wait error = %d\n",err);

    /* unpack the received particles */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

/*-------- sub-step 2: exchange region II (id2)  --------*/

  /* send and receive buffer size to/from region II (id2) */
  /* Post a non-blocking receive for the data size */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id2r,
                             boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

  /* send buffer size */
  err = MPI_Send(&n2, 1, MPI_LONG, id2s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);

  /* receive buffer size */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[receive_ix1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf(cnt_recv+1);

  /* send and receive buffer to/from region II (id2) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id2r,
                              boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }

  if (n2 > 0) {
    /* send buffer */
    err = MPI_Send(send_buf, n2*NVAR_P, MPI_DOUBLE, id2s,
                             boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ix1_particle_shear]: MPI_Send error = %d\n",err);
  }

  if (cnt_recv > 0) {
    /* receive buffer */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[recv_ix1_particle_shear]: MPI_Wait error = %d\n",err);

    /* unpack the received particles */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void shearingbox_ox1_particle(GridS *pG, DomainS *pD,long numpar)
 *  \brief Shearing box boundary condition, Outer x1 boundary (obc_x1=4)
 *
 * This routine works only for ShBoxCoord = xy.
 * Input: numpar: for the packing routine, the array index to start with.
 */
static void shearingbox_ox1_particle(GridS *pG, DomainS *pD, long numpar)
{
#ifdef MPI_PARALLEL
  /* amount of shear, whole (yshear) and fractional (yshift) */
  Real yshear, yshift;
  Real x2len;			/* length in x2 (y) direction of a grid */
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
/*--------------------- MPI case: Step 1. Find locations ---------------------*/

  /* compute the distance the computational domain has sheared in y */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx2);
  x2len = pG->dx2*pG->Nx[1];

  /* obtain processor ids to send and receive */
  ind1 = pD->NGrid[0]-1;	ind3 = my_kproc;
  ind2 = my_jproc - (int)(yshift/x2len) - 1;
  if (ind2 < 0) ind2 += pD->NGrid[1];
  if (ind2 < 0) ind2 += pD->NGrid[1];
  id1s = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region I (send to) */

  ind2 = my_jproc - (int)(yshift/x2len);
  if (ind2 < 0) ind2 += pD->NGrid[1];
  id2s = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region II(send to) */

  ind2 = my_jproc + (int)(yshift/x2len) + 1;
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  id1r = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region I (rcv from)*/

  ind2 = my_jproc + (int)(yshift/x2len);
  if (ind2 > (pD->NGrid[1]-1)) ind2 -= pD->NGrid[1];
  id2r = pD->GData[ind3][ind2][ind1].ID_Comm_Domain;/* for region II(rcv from)*/

/*--------------MPI case: Step 2. Exchange particles -------------------------*/

  /* packing particles at region I */
  n1 = packing_ox1_particle_shear(pG, 1, numpar);

/*-------- sub-step 1: outer region I with inner rigion I (id1) --------*/
  /* send and receive buffer size to/from inner region I (id1) */
  /* Post a non-blocking receive for the data size from inner regiion I */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id1r,
                             boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

  /* send buffer size to inner region I */
  err = MPI_Send(&n1, 1, MPI_LONG, id1s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);

  /* receive buffer size from the inner region I */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[recc_ox1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf(cnt_recv+1);

  /* send and receive buffer to/from inner region I (id1) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from inner region I */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id1r,
                              boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }

  if (n1 > 0) {
    /* send buffer to inner region I */
    err = MPI_Send(send_buf, n1*NVAR_P, MPI_DOUBLE, id1s,
                             boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  }

  /* packing particles at region II (N.B.!) */
  n2 = packing_ox1_particle_shear(pG, 2, numpar);
  if (cnt_recv > 0) {
    /* receive buffer inner region I */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[recv_ox1_particle_shear]: MPI_Wait error = %d\n",err);

    /* unpack the received particle */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

/*-------- sub-step 2: outer region II with inner region II (id2)  --------*/

  /* send and receive buffer size to/from inner region II (id2) */
  /* Post a non-blocking receive for the data size from inner region II */
  err = MPI_Irecv(&cnt_recv, 1, MPI_LONG, id2r,
                             boundary_particle_tag, MPI_COMM_WORLD, &rq);
  if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);

  /* send buffer size to inner region II */
  err = MPI_Send(&n2, 1, MPI_LONG, id2s, boundary_particle_tag, MPI_COMM_WORLD);
  if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);

  /* receive buffer size from inner region II */
  err = MPI_Wait(&rq, &stat);
  if(err) ath_error("[recv_ox1_particle_shear]: MPI_Wait error = %d\n",err);

  /* check space for the receive buffer */
  if ((cnt_recv+1) >= recv_bufsize)
    realloc_recvbuf(cnt_recv+1);

  /* send and receive buffer to/from inner region II (id2) */
  if (cnt_recv > 0) {
    /* Post a non-blocking receive for the data from inner region II */
    err = MPI_Irecv(recv_buf, cnt_recv*NVAR_P, MPI_DOUBLE, id2r,
                              boundary_particle_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals_particle]: MPI_Irecv error = %d\n",err);
  }

  if (n2 > 0) {
    /* send buffer to inner region II */
    err = MPI_Send(send_buf, n2*NVAR_P, MPI_DOUBLE, id2s,
                             boundary_particle_tag, MPI_COMM_WORLD);
    if(err) ath_error("[send_ox1_particle_shear]: MPI_Send error = %d\n",err);
  }

  if (cnt_recv > 0) {
    /* receive buffer from inner region II */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[recv_ox1_particle_shear]: MPI_Wait error = %d\n",err);

    /* unpack the received particle */
    unpack_particle(pG, recv_buf, cnt_recv);
  }

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ix1_particle_shear(GridS *pG, int reg, long numpar)
 *  \brief Packing the particle inside the inner x1 boundary for shearing box
 *
 * Input: pG: grid; reg: region, 1 or 2 for mpi case, 0 for non-mpi case
 *        numpar: array index to start with.
 * Return: number of packed particles in the specified region.
 * Note: this routine serves for only 3D
 */
static long packing_ix1_particle_shear(GridS *pG, int reg, long numpar)
{
  GrainS *gr;		/* current pointer */
  long n, p;
  Real ix1b;		/* coordinate limit in x1 inner boundary */
  /* amount of shear, whole (yshear) and fractional (yshift) */
  Real yshear, yshift;
  /* x2c: y-coordinate marking the demarcation of the two regions */
  Real x20, x2c;
  double *pd;

/*---------------- Step.1 -----------------------*/
  /* get the distance of shear */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx2);
  x20 = pG->MaxX[1];
  x2c = x20 - fmod(yshear, pG->dx2*pG->Nx[1]);

  /* get coordinate limits for particles to be packed*/
  ix1b = pG->MinX[0];

  /* buffer assignment */
  pd = send_buf;
  n = 0;

/*---------------- Step.2 -----------------------*/
  /* loop over all particle to pack particles in the boundary */
  p = numpar;
  while (p<pG->nparticle) {
    gr = &(pG->particle[p]);
    p += 1;
    if ((gr->pos == 21) || ( (gr->pos == 0) && (gr->x1 < ix1b)))
    { /* crossing particle or ghost particle from ox1 */

      if (((reg == 1) && (gr->x2 >= x2c)) || ((reg == 2) && (gr->x2 < x2c)))
      {         /* region I */                      /* region II */
        /* apply the shift */
        gr->x2 = x2min + fmod(gr->x2 - x2min + yshift, Lx2);

        /* pack the particle */
        packing_one_particle(gr, n, gr->pos);
        n += 1;

        /* delete the particle */
        pG->nparticle -= 1;
        grproperty[gr->property].num -= 1;
        p -= 1;
        pG->particle[p] = pG->particle[pG->nparticle];
      }

      if (reg == 0) /* non-mpi case, directly shift the particle positions */
        gr->x2 = x2min + fmod(gr->x2 - x2min + yshift, Lx2);
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static long packing_ox1_particle_shear(GridS *pG, int reg, long numpar)
 *  \brief Packing the particle outside the outer x1 boundary for shearing box
 *
 * Input: pG: grid; reg: region, 1 or 2 for mpi case, 0 for non-mpi case
 *        numpar: array index to start with.
 * Return: number of packed particles in the specified region.
 * Note: this routine serves for only 3D
 */
static long packing_ox1_particle_shear(GridS *pG, int reg, long numpar)
{
  GrainS *gr;		/* current pointer */
  long n, p;
  Real ox1b;		/* coordinate limit of x1 outer boundary */
  /* amount of shear, whole (yshear) and fractional (yshift) */
  Real yshear, yshift;
  /* x2c: y-coordinate marking the demarcation of the two regions */
  Real x20, x2c;
  double *pd;

/*---------------- Step.1 -----------------------*/
  /* get the distance of shear */
  yshear = vshear*pG->time;
  yshift = fmod(yshear, Lx2);
  x20 = pG->MinX[1];
  x2c = x20 + fmod(yshear, pG->dx2*pG->Nx[1]);

  /* get coordinate limits for particles to be packed*/
  ox1b = pG->MaxX[0];

  /* buffer assignment */
  pd = send_buf;
  n = 0;

/*---------------- Step.2 -----------------------*/
  /* loop over all particle to pack particles in the boundary */
  p = numpar;
  while (p<pG->nparticle) {
    gr = &(pG->particle[p]);
    p += 1;
    if ((gr->pos == 11) || ( (gr->pos == 0) && (gr->x1 >= ox1b)))
    { /* crossing particle or ghost particle from ix1 */

      if (((reg == 1) && (gr->x2 < x2c)) || ((reg == 2) && (gr->x2 >= x2c)))
      {         /* region I */                      /* region II */

        /* apply the shift */
        gr->x2 = x2min + fmod(gr->x2 - x2min + Lx2 - yshift, Lx2);

        /* pack the particle */
        packing_one_particle(gr, n, gr->pos);
        n += 1;

        /* delete the particle */
        pG->nparticle -= 1;
        grproperty[gr->property].num -= 1;
        p -= 1;
        pG->particle[p] = pG->particle[pG->nparticle];
      }
      if (reg == 0) /* non-mpi case, directly shift the particle positions */
        gr->x2 = x2min + fmod(gr->x2 - x2min + Lx2 - yshift, Lx2);
    }
  }

  return n;
}

#ifdef FARGO
/*----------------------------------------------------------------------------*/
/*! \fn static long packing_particle_fargo(GridS *pG, Real yl, Real yu)
 *  \brief Packing the particle for FARGO
 *
 * Input: pG: grid; 
 *        yl, yu: lower and upper limit of x2 (y) coordinate for particles to
 *                be packed
 * Return: number of packed particles in the specified region.
 */
static long packing_particle_fargo(GridS *pG, Real yl, Real yu)
{
  GrainS *gr;
  long p, n;
  double *pd;

  p = 0;
  n = 0;
  pd = send_buf;
  while (p<pG->nparticle) {
    gr = &(pG->particle[p]);
    p += 1;
    if ((gr->x2 >= yl) && (gr->x2 < yu))
    { /* pack the particle as it is */
      packing_one_particle(gr, n, gr->pos);
      n += 1;

      /* delete the particle */
      pG->nparticle -= 1;
      grproperty[gr->property].num -= 1;
      p -= 1;
      pG->particle[p] = pG->particle[pG->nparticle];
    }
  }

  return n;
}

/*----------------------------------------------------------------------------*/
/*! \fn static int gridshift(Real shift)
 *  \brief Get the number of grids to be shifted in the y direction */
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

#undef NVAR_P

#endif /*PARTICLES*/
