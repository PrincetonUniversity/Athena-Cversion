#include "../copyright.h"
/*============================================================================*/
/*! \file exchange.c
 *  \brief Exchange the gas-particle coupling array between ghost and
 *         boundary zones.
 *
 * PURPOSE: Particles near grid boundaries deposit their physical properties
 *   partially to the ghost zones. This part of the deposit is to be mapped
 *   to the grid zone. The procedure is opposite to setting boundary conditions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - exchange_gpcouple()
 * - exchange_gpcouple_init()
 * - exchange_gpcouple_fun()
 * - exchange_gpcouple_destruct()
 *
 * PRIVATE FUNCTIONS:
 * - reflecting_???
 * - outflow_???
 * - periodic_???
 * - send_???
 * - receive_???
 * where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
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

#ifdef PARTICLES /* endif at the end of the file */

/* Define the maximum number of variables in the gas-particle coupling array
 * to be exchanged at any step (0-2) */
#define NVar_Max 5
/* Define the maximum number of ghost zones to be exchanged */
#define NLayer_Max 6

/*! \struct GPExc
 *  \brief Define structure which holds variables for gas-particle exchange */
typedef struct GPExc_s{
  Real U[NVar_Max];
}GPExc;

/* number of variables to be exchanged for each cell in a specific step */
static int NVar;
/* number of ghost layers to be exchanged */
static int NExc;
/* nulber of grid layers (offset from boundary) to be copied
 * For pure exchange, set to 0;
 * To fill (copy) particle deposits to ghost zones, set to nghost
 * Total number of layers is thus NLayer = NExc + NOfst.  */
static int NOfst;

/* grid index limit for the exchange */
static int il,iu, jl,ju, kl,ku;
static int ib,it, jb,jt, kb,kt;

/* Temporary array where the exchange operation is executed */
static GPExc ***myCoup=NULL;

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

#ifdef SHEARING_BOX
static Real Delta; /* 0 (beginning) or 0.5 (middle) of a time step */
static Real *Flx=NULL;
static Real *UBuf=NULL;
static GPExc ***GhstZns_ix1=NULL;
static GPExc ***GhstZns_ox1=NULL;
static GPExc ***TempZns=NULL;
extern void RemapFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);
#endif

/* boundary condition function pointers. local to this function  */
static VGFun_t apply_ix1 = NULL, apply_ox1 = NULL;
static VGFun_t apply_ix2 = NULL, apply_ox2 = NULL;
static VGFun_t apply_ix3 = NULL, apply_ox3 = NULL;

/*====================== PROTOTYPE OF PRIVATE FUNCTIONS ======================*/
/*----------------------------------------------------------------------------*/

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - apply reflecting BCs at boundary ???
 *   outflow_???()  - apply outflow BCs at boundary ???
 *   periodic_???() - apply periodic BCs at boundary ???
 *   pack_???()     - pack data at ??? boundary
 *   unpack_???()   - unpack data at ??? boundary
 *============================================================================*/

static void reflect_ix1_exchange(GridS *pG);
static void reflect_ox1_exchange(GridS *pG);
static void reflect_ix2_exchange(GridS *pG);
static void reflect_ox2_exchange(GridS *pG);
static void reflect_ix3_exchange(GridS *pG);
static void reflect_ox3_exchange(GridS *pG);

static void outflow_exchange(GridS *pG);

static void periodic_ix1_exchange(GridS *pG);
static void periodic_ox1_exchange(GridS *pG);
static void periodic_ix2_exchange(GridS *pG);
static void periodic_ox2_exchange(GridS *pG);
static void periodic_ix3_exchange(GridS *pG);
static void periodic_ox3_exchange(GridS *pG);

#ifdef SHEARING_BOX
static void Remap_exchange_ix1(DomainS *pD);
static void Remap_exchange_ox1(DomainS *pD);
#endif

#ifdef MPI_PARALLEL
static void pack_ix1_exchange(GridS *pG);
static void pack_ox1_exchange(GridS *pG);
static void pack_ix2_exchange(GridS *pG);
static void pack_ox2_exchange(GridS *pG);
static void pack_ix3_exchange(GridS *pG);
static void pack_ox3_exchange(GridS *pG);

static void unpack_ix1_exchange(GridS *pG);
static void unpack_ox1_exchange(GridS *pG);
static void unpack_ix2_exchange(GridS *pG);
static void unpack_ox2_exchange(GridS *pG);
static void unpack_ix3_exchange(GridS *pG);
static void unpack_ox3_exchange(GridS *pG);

#ifdef SHEARING_BOX
static void pack_ix2_remap(GridS *pG, GPExc ***myZns);
static void pack_ox2_remap(GridS *pG, GPExc ***myZns);
static void unpack_ix2_remap(GridS *pG, GPExc ***myZns);
static void unpack_ox2_remap(GridS *pG, GPExc ***myZns);
#endif
#endif /* MPI_PARALLEL */


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void exchange_feedback(DomainS *pD, short step)
 *  \brief Calls appropriate functions to copy feedback in the ghost
 *    zones back to the grid.
 *
 *    The function pointers (*apply_???) are set during
 *    initialization by exchange_feedback_init() to be either a user-defined
 *    function, or one of the functions corresponding to reflecting, periodic,
 *    or outflow.  If the left- or right-Grid ID numbers are >= 1 (neighboring
 *    grids exist), then MPI calls are used.
 *
 *  Order for updating boundary conditions must always be x3-x2-x1 in order to
 *  fill the corner cells properly (opposite to setting MHD B.C.!)
 */

void exchange_gpcouple(DomainS *pD, short lab)
{
  GridS *pG = pD->Grid;
  int i,j,k;
#ifdef SHEARING_BOX
  int n,myL,myM,myN,BCFlag;
#ifndef FARGO
  Real Lx = pD->RootMaxX[0]-pD->RootMinX[0];
#endif
#endif
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr, mIndex;
#endif /* MPI_PARALLEL */
	
/*--- Step 1. ------------------------------------------------------------------	
 * Copy the information in the Gas-Particle coupling array into temporary array
 * This step depends on the parameter "lab", where for
 * lab = 0: particle binning for output purpose
 * lab = 1: predictor step of feedback exchange
 * lab = 2: corrector step of feedback exchange
 * All the operations in this routine are performed on the temporary array,
 * which will be copied back to the main array GPCoup at the end.
 *----------------------------------------------------------------------------*/
 
#ifdef SHEARING_BOX
  Delta = 0.0; 
#endif

  switch (lab) {
    case 0: /* particle binning for output purpose */
		  
      NVar = 4; NExc = 1; NOfst = 0;

      for (k=klp; k<=kup; k++) {
       for (j=jlp; j<=jup; j++) {
        for (i=ilp; i<=iup; i++) {
	  myCoup[k][j][i].U[0]=pG->Coup[k][j][i].grid_v1;
          myCoup[k][j][i].U[1]=pG->Coup[k][j][i].grid_v2;
          myCoup[k][j][i].U[2]=pG->Coup[k][j][i].grid_v3;
          myCoup[k][j][i].U[3]=pG->Coup[k][j][i].grid_d;
      }}}
      break;
		  
#ifdef FEEDBACK		  
    case 1: /* predictor step of feedback exchange */
		  
      NVar = 5; NExc = 1; NOfst = nghost;

      for (k=klp; k<=kup; k++) {
       for (j=jlp; j<=jup; j++) {
        for (i=ilp; i<=iup; i++) {
          myCoup[k][j][i].U[0]=pG->Coup[k][j][i].fb1;
          myCoup[k][j][i].U[1]=pG->Coup[k][j][i].fb2;
          myCoup[k][j][i].U[2]=pG->Coup[k][j][i].fb3;
          myCoup[k][j][i].U[3]=pG->Coup[k][j][i].FBstiff;
          myCoup[k][j][i].U[4]=pG->Coup[k][j][i].Eloss;
      }}}
      break;
		  
    case 2: /* corrector step of feedback exchange */
		  
      NVar = 4; NExc = 2; NOfst = 0;
#ifdef SHEARING_BOX
#ifndef FARGO
      Delta = 0.5; /* at the middle of a time step */
#endif
#endif
      for (k=klp; k<=kup; k++) {
       for (j=jlp; j<=jup; j++) {
        for (i=ilp; i<=iup; i++) {
          myCoup[k][j][i].U[0]=pG->Coup[k][j][i].fb1;
          myCoup[k][j][i].U[1]=pG->Coup[k][j][i].fb2;
          myCoup[k][j][i].U[2]=pG->Coup[k][j][i].fb3;
          myCoup[k][j][i].U[3]=pG->Coup[k][j][i].Eloss;		  
      }}}
      break;

    default:
      ath_perr(-1,"[exchange_GPCouple]: lab must be equal to 0, 1, or 2!\n");
#else
    default:
      ath_perr(-1,"[exchange_GPCouple]: lab must be equal to 0!\n");
#endif /* FEEDBACK */
  }

/* set left and right grid indices */
  if (pG->Nx[0] > 1) {
    il = pG->is - NExc;         iu = pG->ie + NExc;
    ib = pG->is - NOfst;        it = pG->ie + NOfst;
  } else {
    il = ib = pG->is;           iu = it = pG->ie;
  }   
      
  if (pG->Nx[1] > 1) {
    jl = pG->js - NExc;         ju = pG->je + NExc;
    jb = pG->js - NOfst;        jt = pG->je + NOfst;
  } else {
    jl = jb = pG->js;           ju = jt = pG->je;
  }   
  
  if (pG->Nx[2] > 1) {
    kl = pG->ks - NExc;         ku = pG->ke + NExc;
    kb = pG->ks - NOfst;        kt = pG->ke + NOfst;
  } else {
    kl = kb = pG->ks;           ku = kt = pG->ke;
  }

/*--- Step 2. ------------------------------------------------------------------
 * Feedback exchange in x3-direction */

  if (pG->Nx[2] > 1){

#ifdef MPI_PARALLEL
    cnt1 = pG->Nx[0] > 1 ? pG->Nx[0] + 2*NExc : 1;
    cnt2 = pG->Nx[1] > 1 ? pG->Nx[1] + 2*NExc : 1;
    cnt = (NExc+NOfst)*cnt1*cnt2*NVar;

/* MPI blocks to both left and right */
    if (pG->rx3_id >= 0 && pG->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_exchange(pG);
      if (mIndex == 1) unpack_ox3_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_exchange(pG);
      if (mIndex == 1) unpack_ox3_exchange(pG);

    }

    /* Physical boundary on left, MPI block on right */
    if (pG->rx3_id >= 0 && pG->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix3)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);
      
      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_exchange(pG);
      
    }

    /* MPI block on left, Physical boundary on right */
    if (pG->rx3_id < 0 && pG->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox3)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_exchange(pG);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx3_id < 0 && pG->lx3_id < 0) {
      (*apply_ix3)(pG);
      (*apply_ox3)(pG);
    }

  }

/*--- Step 3. ------------------------------------------------------------------
 * Feedback exchange in x2-direction */

  if (pG->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt1 = pG->Nx[0] > 1 ? pG->Nx[0] + 2*NExc : 1;
    cnt3 = pG->Nx[2] > 1 ? pG->Nx[2] + 2*NOfst: 1;
    cnt = (NExc+NOfst)*cnt1*cnt3*NVar;

/* MPI blocks to both left and right */
    if (pG->rx2_id >= 0 && pG->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_exchange(pG);
      if (mIndex == 1) unpack_ox2_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_exchange(pG);
      if (mIndex == 1) unpack_ox2_exchange(pG);

    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx2_id >= 0 && pG->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix2)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_exchange(pG);

    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx2_id < 0 && pG->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox2)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_exchange(pG);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx2_id < 0 && pG->lx2_id < 0) {
      (*apply_ix2)(pG);
      (*apply_ox2)(pG);
    }

/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);

    /* shearing-box in xy (2D and 3D) */
    if (pG->Nx[2] > 1 || ShBoxCoord==xy)
    {
      /* inner x1 boundary */
      BCFlag = par_geti_def("domain1","bc_ix1",0);
      if (myL == 0 && BCFlag == 4) {
        Remap_exchange_ix1(pD);
#ifndef FARGO
        if (lab == 0) { /* for output */
          for (k=kb; k<=kt; k++) {
           for (j=jb; j<=jt; j++) {
            for (i=0; i<NExc+NOfst; i++) {
              GhstZns_ix1[k][i][j].U[1] -= qshear*Omega_0*Lx;
          }}}
        }
#endif
      } else {
        for (k=kb; k<=kt; k++) {
         for (j=jb; j<=jt; j++) {
          for (i=0; i<NExc+NOfst; i++) {
           for (n=0;n<NVar; n++) {
            GhstZns_ix1[k][i][j].U[n] = myCoup[k][j][pG->is-NExc+i].U[n];
        }}}}
      }

      /* outer x1 boundary */
      BCFlag = par_geti_def("domain1","bc_ox1",0);
      if (myL == ((pD->NGrid[0])-1) && BCFlag == 4) {
        Remap_exchange_ox1(pD);
#ifndef FARGO
        if (lab == 0) { /* no fargo, for output */
          for (k=kb; k<=kt; k++) {
           for (j=jb; j<=jt; j++) {
            for (i=0; i<NExc+NOfst; i++) {
              GhstZns_ox1[k][i][j].U[1] += qshear*Omega_0*Lx;
          }}}
        }
#endif
      } else {
        for (k=kb; k<=kt; k++) {
         for (j=jb; j<=jt; j++) {
          for (i=0; i<NExc+NOfst; i++) {
           for (n=0;n<NVar; n++) {
            GhstZns_ox1[k][i][j].U[n] = myCoup[k][j][pG->ie-NOfst+i+1].U[n];
        }}}}
      }
    }
#ifndef FARGO
    /* 2D shearing-box in x-z, no fargo, for output */
    else if (lab == 0) 
    {
      if (myL == 0) {
        for (j=jb; j<=jt; j++) {
         for (i=pG->is-NExc; i<pG->is; i++) {
           myCoup[pG->ks][j][i].U[2] -= qshear*Omega_0*Lx;
        }}
      }
      if (myL == ((pD->NGrid[0])-1)) {
        for (j=jb; j<=jt; j++) {
         for (i=pG->ie+1; i<=pG->ie+NExc; i++) {
           myCoup[pG->ks][j][i].U[2] += qshear*Omega_0*Lx;
        }}
      }
    }
#endif
#endif /* SHEARING_BOX */
  }

/*--- Step 4. ------------------------------------------------------------------
 * Feedback exchange in x1-direction */

  if (pG->Nx[0] > 1){

#ifdef MPI_PARALLEL
    cnt2 = pG->Nx[1] > 1 ? pG->Nx[1] + 2*NOfst : 1;
    cnt3 = pG->Nx[2] > 1 ? pG->Nx[2] + 2*NOfst : 1;
    cnt = (NExc+NOfst)*cnt2*cnt3*NVar;

/* MPI blocks to both left and right */
    if (pG->rx1_id >= 0 && pG->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_exchange(pG);
      if (mIndex == 1) unpack_ox1_exchange(pG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_exchange(pG);
      if (mIndex == 1) unpack_ox1_exchange(pG);
    }

/* Physical boundary on left, MPI block on right */
    if (pG->rx1_id >= 0 && pG->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*apply_ix1)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_exchange(pG);

    }

/* MPI block on left, Physical boundary on right */
    if (pG->rx1_id < 0 && pG->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1_exchange(pG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*apply_ox1)(pG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_exchange(pG);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pG->rx1_id < 0 && pG->lx1_id < 0) {
      (*apply_ix1)(pG);
      (*apply_ox1)(pG);
    } 
  }

/*--- Step 5. ------------------------------------------------------------------	
 * Copy the variables from the temporary array where exchange has finished back
 * to the Gas-Particle coupling array. Again, for
 * lab = 0: particle binning for output purpose
 * lab = 1: predictor step of feedback exchange
 * lab = 2: corrector step of feedback exchange
 *----------------------------------------------------------------------------*/
	
  switch (lab) {
    case 0:	/* particle binning for output purpose */
      for (k=kb; k<=kt; k++) {
       for (j=jb; j<=jt; j++) {
        for (i=ib; i<=it; i++) {
          pG->Coup[k][j][i].grid_v1= myCoup[k][j][i].U[0];
          pG->Coup[k][j][i].grid_v2= myCoup[k][j][i].U[1];
          pG->Coup[k][j][i].grid_v3= myCoup[k][j][i].U[2];
          pG->Coup[k][j][i].grid_d = myCoup[k][j][i].U[3];
      }}}
      break;
			
#ifdef FEEDBACK		  
    case 1: /* predictor step of feedback exchange */
      for (k=kb; k<=kt; k++) {
       for (j=jb; j<=jt; j++) {
        for (i=ib; i<=it; i++) {
          pG->Coup[k][j][i].fb1    = myCoup[k][j][i].U[0];
          pG->Coup[k][j][i].fb2    = myCoup[k][j][i].U[1];
          pG->Coup[k][j][i].fb3    = myCoup[k][j][i].U[2];
          pG->Coup[k][j][i].FBstiff= myCoup[k][j][i].U[3];
          pG->Coup[k][j][i].Eloss  = myCoup[k][j][i].U[4];
      }}}
      break;
			
    case 2: /* corrector step of feedback exchange */
      for (k=kb; k<=kt; k++) {
       for (j=jb; j<=jt; j++) {
        for (i=ib; i<=it; i++) {
          pG->Coup[k][j][i].fb1  = myCoup[k][j][i].U[0];
          pG->Coup[k][j][i].fb2  = myCoup[k][j][i].U[1];
          pG->Coup[k][j][i].fb3  = myCoup[k][j][i].U[2];
          pG->Coup[k][j][i].Eloss= myCoup[k][j][i].U[3];		  
      }}}
      break;
			
    default:
      ath_perr(-1,"[exchange_GPCouple]: lab must be equal to 0, 1, or 2!\n");
#else
    default:
      ath_perr(-1,"[exchange_GPCouple]: lab must be equal to 0!\n");
#endif /* FEEDBACK */
  }
	
  return;

}

#ifdef SHEARING_BOX
/*----------------------------------------------------------------------------*/
/*! \fn void Remap_exchange_ix1(DomainS *pD);
 *  \brief Remap of the gas-particle coupling array for 3D shearing-box in ix1.  
 *
 * It applies a remap in Y for the radial ghost cells before they are sent to
 * the opposite side exchange_gpcouple.c
 *
 * This is a public function which is called by exchange_gpcouple() inside a
 * SHEARING_BOX macro.                                                        */
/*----------------------------------------------------------------------------*/
void Remap_exchange_ix1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int n,i,j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids,mIndex,Nx3;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ix1()  */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;
  
  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;
  
  qomL = qshear*Omega_0*Lx;
  yshear = qomL*(pG->time+Delta*pG->dt);
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;
  
#ifdef MPI_PARALLEL
  if (pG->Nx[2] > 1) {
    Nx3 = pG->Nx[2]+2*NOfst;
  } else {
    Nx3 = 1;
  }
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Copy myCoup to the ghost zone array where j is the fastest varying index.
 * Apply periodic BC in x2 to the ghost array.  Requires MPI calls if 
 * NGrid_x2 > 1 */
      
  for(k=kb; k<=kt; k++){
   for(j=jb; j<=jt; j++){
    for(i=0; i<NExc+NOfst; i++){
      for (n=0; n<NVar; n++){
	TempZns[k][i][j].U[n] = myCoup[k][j][is-NExc+i].U[n];
      }     
  }}}     

  if (NOfst < nghost) {
    if (pD->NGrid[1] == 1) {
  
      for(k=kb; k<=kt; k++){
       for(i=0; i<NExc+NOfst; i++){
        for(j=1; j<=nghost; j++){
	 for (n=0; n<NVar; n++){
	   TempZns[k][i][js-j].U[n] = TempZns[k][i][je-(j-1)].U[n];
	   TempZns[k][i][je+j].U[n] = TempZns[k][i][js+(j-1)].U[n];
	 }
      }}}     
        
#ifdef MPI_PARALLEL
    } else {

/* MPI calls to swap data */

      cnt = nghost*(NExc+NOfst)*Nx3*NVar;
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_remap(pG, TempZns);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_remap(pG, TempZns);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, TempZns);
      if (mIndex == 1) unpack_ox2_remap(pG, TempZns);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, TempZns);
      if (mIndex == 1) unpack_ox2_remap(pG, TempZns);

#endif /* MPI_PARALLEL */
    }
  } /* If NOfst < nghost */

/*--- Step 3. ------------------------------------------------------------------
 * Apply a conservative remap of solution over the fractional part of a cell */

  for(k=kb; k<=kt; k++){
   for(i=0; i<NExc+NOfst; i++){
    for (n=0; n<NVar; n++) {
      for(j=js-nghost; j<=je+nghost; j++){
	UBuf[j] = TempZns[k][i][j].U[n];
      }
      RemapFlux(UBuf,epso,js,je+1,Flx);
      for(j=js; j<=je; j++){
        TempZns[k][i][j].U[n] -= (Flx[j+1] - Flx[j]);
      }
  }}}

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from the ghost array back into myCoup.  */

  if (pD->NGrid[1] == 1) {

    for(k=kb; k<=kt; k++){
     for(i=0; i<NExc+NOfst; i++){
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx[1];
        for(n=0; n<NVar; n++)
         GhstZns_ix1[k][i][j].U[n] = TempZns[k][i][jremap].U[n];
       }
    }}

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from the ghost array back into myCoup.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(joverlap-1):je] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZns */

      cnt = (NExc+NOfst)*joverlap*Nx3*NVar;
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(&(recv_buf[0][0]), cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &(recv_rq[0]));

      pSnd = &(send_buf[0][0]);
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for (j=js; j<=js+(joverlap-1); j++) {
	  for (n=0; n<NVar; n++) {
            *(pSnd++) = TempZns[k][i][j].U[n];
	  }
      }}}
      ierr = MPI_Send(&(send_buf[0][0]),cnt,MPI_DOUBLE,sendto_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in myCoup
 */

      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);

      pRcv = &(recv_buf[0][0]);

      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
	for (j=je-(joverlap-1); j<=je; j++) {
	  for (n=0; n<NVar; n++) {
            GhstZns_ix1[k][i][j].U[n] = *(pRcv++);
	  }
      }}}
    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tJyBuf into tJy.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

    if (Ngrids == 0) {

      for(k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for(j=js; j<=je-joverlap; j++){
	  jremap = j+joverlap;
	  for (n=0; n<NVar; n++) {
            GhstZns_ix1[k][i][j].U[n]  = TempZns[k][i][jremap].U[n];
          }
      }}}

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZns (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = (NExc+NOfst)*(pG->Nx[1]-joverlap)*Nx3*NVar;
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]), cnt, MPI_DOUBLE, getfrom_id,
                       shearing_sheet_ix1_tag, pD->Comm_Domain, &(recv_rq[1]));

      pSnd = &(send_buf[1][0]);

      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
	for (j=js; j<=je-joverlap; j++) {
	  for (n=0; n<NVar; n++) {
            *(pSnd++) = TempZns[k][i][j].U[n];
	  }
      }}}
      ierr = MPI_Send(&(send_buf[1][0]),cnt,MPI_DOUBLE,sendto_id,
		      shearing_sheet_ix1_tag, pD->Comm_Domain);

/* unpack data sent from [js+overlap:je], and remap into cells in
 * [js:je-joverlap] in myCoup */

      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);

      pRcv = &(recv_buf[1][0]);
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
	for (j=js; j<=je-joverlap; j++) {
	  for (n=0; n<NVar; n++) {
	    GhstZns_ix1[k][i][j].U[n] = *(pRcv++);
	  }
      }}}
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*-----Step 6. Apply periodic BC in Y if necessary-----------------
 */
  if (NOfst > 0) {
    if (pD->NGrid[1] == 1) {
      
      for(k=kb; k<=kt; k++){
       for(i=0; i<NExc+NOfst; i++){
        for(j=1; j<=nghost; j++){
         for (n=0; n<NVar; n++){
           GhstZns_ix1[k][i][js-j].U[n] = GhstZns_ix1[k][i][je-(j-1)].U[n];
           GhstZns_ix1[k][i][je+j].U[n] = GhstZns_ix1[k][i][js+(j-1)].U[n];
         }
      }}}     
     
#ifdef MPI_PARALLEL
    } else {

/* MPI calls to swap data */

      cnt = nghost*(NExc+NOfst)*Nx3*NVar;
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_remap(pG, GhstZns_ix1);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_remap(pG, GhstZns_ix1);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, GhstZns_ix1);
      if (mIndex == 1) unpack_ox2_remap(pG, GhstZns_ix1);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, GhstZns_ix1);
      if (mIndex == 1) unpack_ox2_remap(pG, GhstZns_ix1);
    
#endif /* MPI_PARALLEL */
    }
  } /* If NOfst > 0 */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void Remap_exchange_ox1(DomainS *pD);
 *  \brief Remap of the gas-particle coupling array for 3D shearing-box in ox1.  
 *
 * It applies a remap in Y for the radial ghost cells before they are sent to
 * the opposite side exchange_gpcouple.c
 *
 * This is a public function which is called by exchange_gpcouple() inside a
 * SHEARING_BOX macro.                                                        */
/*----------------------------------------------------------------------------*/
void Remap_exchange_ox1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int n,i,j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids,mIndex,Nx3;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ix1()  */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*(pG->time+Delta*pG->dt);
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

#ifdef MPI_PARALLEL
  if (pG->Nx[2] > 1) {
    Nx3 = pG->Nx[2]+2*NOfst;
  } else {
    Nx3 = 1;
  }
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Copy myCoup to the ghost zone array where j is the fastest varying index.
 * Apply periodic BC in x2 to the ghost array.  Requires MPI calls if 
 * NGrid_x2 > 1 */
	
  for(k=kb; k<=kt; k++){
   for(j=jb; j<=jt; j++){
    for(i=0; i<NExc+NOfst; i++){
      for (n=0; n<NVar; n++){
	TempZns[k][i][j].U[n] = myCoup[k][j][ie-NOfst+i+1].U[n];
      }
  }}}

  if (NOfst < nghost) {
    if (pD->NGrid[1] == 1) {

      for(k=kb; k<=kt; k++){
       for(i=0; i<NExc+NOfst; i++){
        for(j=1; j<=nghost; j++){
	  for (n=0; n<NVar; n++){
	    TempZns[k][i][js-j].U[n] = TempZns[k][i][je-(j-1)].U[n];
            TempZns[k][i][je+j].U[n] = TempZns[k][i][js+(j-1)].U[n];
	  }
      }}}

#ifdef MPI_PARALLEL
    } else {

/* MPI calls to swap data */

      cnt = nghost*(NExc+NOfst)*Nx3*NVar;
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_remap(pG, TempZns);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_remap(pG, TempZns);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, TempZns);
      if (mIndex == 1) unpack_ox2_remap(pG, TempZns);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, TempZns);
      if (mIndex == 1) unpack_ox2_remap(pG, TempZns);

#endif /* MPI_PARALLEL */
    }
  } /* If NOfst < nghost */
  
/*--- Step 3. ------------------------------------------------------------------
 * Apply a conservative remap of solution over the fractional part of a cell */

  for(k=kb; k<=kt; k++){
   for(i=0; i<NExc+NOfst; i++){
    for (n=0; n<NVar; n++) {
      for(j=js-nghost; j<=je+nghost; j++){
	UBuf[j] = TempZns[k][i][j].U[n];
      }
      RemapFlux(UBuf,epsi,js,je+1,Flx);
      for(j=js; j<=je; j++){
        TempZns[k][i][j].U[n] -= (Flx[j+1] - Flx[j]);
      }
  }}}

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from the ghost array back into myCoup.  */

  if (pD->NGrid[1] == 1) {

    for(k=kb; k<=kt; k++){
     for(i=0; i<NExc+NOfst; i++){
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx[1];
	for (n=0; n<NVar; n++)
          GhstZns_ox1[k][i][j].U[n]  = TempZns[k][i][jremap].U[n];
    }}}

#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from the ghost array back into myCoup.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZns */

      cnt = (NExc+NOfst)*joverlap*Nx3*NVar;
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(&(recv_buf[0][0]), cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &(recv_rq[0]));

      pSnd = &(send_buf[0][0]);
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for (j=je-(joverlap-1); j<=je; j++) {
	 for (n=0; n<NVar; n++) {
          *(pSnd++) = TempZns[k][i][j].U[n];
	 }
     }}}
     ierr = MPI_Send(&(send_buf[0][0]),cnt,MPI_DOUBLE,sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in myCoup
 */

      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);

      pRcv = &(recv_buf[0][0]);
	
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for (j=js; j<=js+(joverlap-1); j++) {
          for (n=0; n<NVar; n++) {
	    GhstZns_ox1[k][i][j].U[n] = *(pRcv++);
	  }
      }}}
    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tJyBuf into tJy.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for(j=js+joverlap; j<=je; j++){
	  jremap = j-joverlap;
	  for (n=0; n<NVar; n++) {
            GhstZns_ox1[k][i][j].U[n]  = TempZns[k][i][jremap].U[n];
	  }
      }}}

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZns (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = (NExc+NOfst)*(pG->Nx[1]-joverlap)*Nx3*NVar;
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]), cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &(recv_rq[1]));

      pSnd = &(send_buf[1][0]);
      
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for (j=js; j<=je-joverlap; j++) {
	  for (n=0; n<NVar; n++) {
	    *(pSnd++) = TempZns[k][i][j].U[n];
          }
      }}}
      ierr = MPI_Send(&(send_buf[1][0]),cnt,MPI_DOUBLE,sendto_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in myCoup */

      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);

      pRcv = &(recv_buf[1][0]);
      for (k=kb; k<=kt; k++) {
       for(i=0; i<NExc+NOfst; i++){
        for (j=js+joverlap; j<=je; j++) {
	  for (n=0; n<NVar; n++) {
	    GhstZns_ox1[k][i][j].U[n] = *(pRcv++);
	  }
      }}}
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*-----Step 6. Apply periodic BC in Y if necessary-----------------
 */   
  if (NOfst > 0) {
    if (pD->NGrid[1] == 1) {
    
      for(k=kb; k<=kt; k++){
       for(i=0; i<NExc+NOfst; i++){
        for(j=1; j<=nghost; j++){
          for (n=0; n<NVar; n++){
            GhstZns_ox1[k][i][js-j].U[n] = GhstZns_ox1[k][i][je-(j-1)].U[n];
            GhstZns_ox1[k][i][je+j].U[n] = GhstZns_ox1[k][i][js+(j-1)].U[n];
          }
      }}}
  
#ifdef MPI_PARALLEL
    } else {

/* MPI calls to swap data */

      cnt = nghost*(NExc+NOfst)*Nx3*NVar;
      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_remap(pG, GhstZns_ox1);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_remap(pG, GhstZns_ox1);
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, GhstZns_ox1);
      if (mIndex == 1) unpack_ox2_remap(pG, GhstZns_ox1);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_remap(pG, GhstZns_ox1);
      if (mIndex == 1) unpack_ox2_remap(pG, GhstZns_ox1);

#endif /* MPI_PARALLEL */
    }
  } /* if NOfst > 0 */ 

  return;
}
#endif /* SHEARING_BOX */

/*----------------------------------------------------------------------------*/
/*! \fn void exchange_gpcouple_init(MeshS *pM)
 *  \brief Sets function pointers for the exchange of gas-particle coupling
 *   array, allocates memory, etc. for the initialization
 */
void exchange_gpcouple_init(MeshS *pM)
{
  int ibc_x1, obc_x1; /* x1 inner and outer boundary condition flag */
  int ibc_x2, obc_x2; /* x2 inner and outer boundary condition flag */
  int ibc_x3, obc_x3; /* x3 inner and outer boundary condition flag */
  int N1T,N2T,N3T;
#ifdef MPI_PARALLEL
  int i,j,k;
  int x1cnt, x2cnt, x3cnt; /* Number of Gas passed in x1-, x2-, x3-dir. */
  int nx1t, nx2t, nx3t, size;
#endif /* MPI_PARALLEL */

  GridS *pG;
  DomainS *pD;

  if (pM->NLevels > 1)
    ath_error("[exchange_init]: particle module does not suport SMR\n");

  pD = &(pM->Domain[0][0]);
  pG = pD->Grid;

/* Set function pointers for physical boundaries in x1-direction */

  if(pG->Nx[0] > 1) {
    if(apply_ix1 == NULL){

      switch(pM->BCFlag_ix1){

      case 1: /* Reflecting */
	apply_ix1 = reflect_ix1_exchange;
	break;
      case 5: /* Reflecting */
	apply_ix1 = reflect_ix1_exchange;
	break;

      case 2: /* Outflow */
	apply_ix1 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ix1 = periodic_ix1_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ix1 = %d unknown\n",
                    pM->BCFlag_ix1);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox1 == NULL){

      switch(pM->BCFlag_ox1){

      case 1: /* Reflecting */
	apply_ox1 = reflect_ox1_exchange;
	break;
      case 5: /* Reflecting */
	apply_ox1 = reflect_ox1_exchange;
	break;

      case 2: /* Outflow */
	apply_ox1 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ox1 = periodic_ox1_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ox1 = %d unknown\n",
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
	apply_ix2 = reflect_ix2_exchange;
	break;
      case 5: /* Reflecting */
	apply_ix2 = reflect_ix2_exchange;
	break;

      case 2: /* Outflow */
	apply_ix2 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ix2 = periodic_ix2_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ix2 = %d unknown\n",
                    pM->BCFlag_ix2);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox2 == NULL){

      switch(pM->BCFlag_ox2){

      case 1: /* Reflecting */
	apply_ox2 = reflect_ox2_exchange;
	break;
      case 5: /* Reflecting */
	apply_ox2 = reflect_ox2_exchange;
	break;

      case 2: /* Outflow */
	apply_ox2 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ox2 = periodic_ox2_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ox2 = %d unknown\n",
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
	apply_ix3 = reflect_ix3_exchange;
	break;
      case 5: /* Reflecting */
	apply_ix3 = reflect_ix3_exchange;
	break;

      case 2: /* Outflow */
	apply_ix3 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ix3 = periodic_ix3_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ix3 = %d unknown\n",
                    pM->BCFlag_ix3);
	exit(EXIT_FAILURE);
      }
    }

    if(apply_ox3 == NULL){

      switch(pM->BCFlag_ox3){

      case 1: /* Reflecting */
	apply_ox3 = reflect_ox3_exchange;
	break;
      case 5: /* Reflecting */
	apply_ox3 = reflect_ox3_exchange;
	break;

      case 2: /* Outflow */
	apply_ox3 = outflow_exchange;
	break;

      case 4: /* Periodic */
	apply_ox3 = periodic_ox3_exchange;
	break;

      default:
	ath_perr(-1,"[exchange_init]: bc_ox3 = %d unknown\n",
                    pM->BCFlag_ox3);
	exit(EXIT_FAILURE);
      }
    }
  }
	
  if (pG->Nx[0] > 1) {
      N1T = pG->Nx[0]+2*nghost;
  } else {
      N1T=1;
  }
  if (pG->Nx[1] > 1) {
      N2T = pG->Nx[1]+2*nghost;
  } else {
      N2T=1;
  }
  if (pG->Nx[2] > 1) {
      N3T = pG->Nx[2]+2*nghost;
  } else {
      N3T=1;
  }

  if ((myCoup = (GPExc***)calloc_3d_array(N3T,N2T,N1T, sizeof(GPExc))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the myCoup array.\n");
	
#ifdef SHEARING_BOX
  if ((Flx = (Real*)calloc_1d_array(N2T,sizeof(Real))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the Flux array.\n");
  if ((UBuf = (Real*)calloc_1d_array(N2T,sizeof(Real))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the UBuf array.\n");
  if ((GhstZns_ix1 = (GPExc***)calloc_3d_array(N3T,NLayer_Max,N2T,sizeof(GPExc))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the Ghost array.\n");
  if ((GhstZns_ox1 = (GPExc***)calloc_3d_array(N3T,NLayer_Max,N2T,sizeof(GPExc))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the Ghost array.\n");
  if ((TempZns = (GPExc***)calloc_3d_array(N3T,NLayer_Max,N2T,sizeof(GPExc))) == NULL)
    ath_error("[exchange_init]: Failed to allocate the Ghost array.\n");
#endif
	
#ifdef MPI_PARALLEL
  x1cnt = x2cnt = x3cnt = 0;

  for (k=0; k<(pD->NGrid[2]); k++){
    for (j=0; j<(pD->NGrid[1]); j++){
      for (i=0; i<(pD->NGrid[0]); i++){
        if(pD->NGrid[2] > 1){
          nx1t = pD->GData[k][j][i].Nx[0];
          if(nx1t > 1) nx1t += 2*NLayer_Max;

          nx2t = pD->GData[k][j][i].Nx[1];
          if(nx2t > 1) nx2t += 2*NLayer_Max;

          x3cnt = nx1t*nx2t > x3cnt ? nx1t*nx2t : x3cnt;
        }

        if(pD->NGrid[1] > 1){
          nx1t = pD->GData[k][j][i].Nx[0];
          if(nx1t > 1) nx1t += 2*NLayer_Max;

          nx3t = pD->GData[k][j][i].Nx[2];
          if(nx3t > 1) nx3t += 2*NLayer_Max;

          x2cnt = nx1t*nx3t > x2cnt ? nx1t*nx3t : x2cnt;
    	}

        if(pD->NGrid[0] > 1){
    	  nx2t = pD->GData[k][j][i].Nx[1];
          if(nx2t > 1) nx2t += 2*NLayer_Max;

          nx3t = pD->GData[k][j][i].Nx[2];
          if(nx3t > 1) nx3t += 2*NLayer_Max;

          x1cnt = nx2t*nx3t > x1cnt ? nx2t*nx3t : x1cnt;
    	}

      }
    }
  }

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= NLayer_Max*NVar_Max; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[exchange_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[exchange_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[exchange_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[exchange_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void exchange_gpcouple_fun(enum Direction dir, VGFun_t prob_bc)
 *  \brief Sets function pointers for user-defined exchange gas-particle
 *   coupling in the problem generator
 */

void exchange_gpcouple_fun(enum BCDirection dir, VGFun_t prob_bc)
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

/*! \fn void exchange_gpcouple_destruct(GridS *pG, Domain *pD)
 *  \brief Finalize the exchange of gas-particle coupling */
void exchange_gpcouple_destruct(MeshS *pM)
{
  apply_ix1 = NULL;
  apply_ox1 = NULL;
  apply_ix2 = NULL;
  apply_ox2 = NULL;
  apply_ix3 = NULL;
  apply_ox3 = NULL;
  free_3d_array(myCoup);
#ifdef SHEARING_BOX
  free_1d_array(Flx);
  free_1d_array(UBuf);
  free_3d_array(GhstZns_ix1);
  free_3d_array(GhstZns_ox1);
  free_3d_array(TempZns);
#endif
#ifdef MPI_PARALLEL
  free_2d_array(send_buf);
  free_2d_array(recv_buf);
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
/*! \fn static void reflect_ix3_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x3 boundary (ibc_x3=1,5)
 */

static void reflect_ix3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;

  for (k=kl; k<pG->ks; k++) {
    kr = 2*pG->ks-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[kr][j][i].U[2] -= myCoup[k][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }

  for (k=kb; k<pG->ks; k++) {
    kr = 2*pG->ks-k-1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[2] = -myCoup[kr][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
      } 
    }   
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox3_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x3 boundary (obc_x3=1,5)
 */

static void reflect_ox3_exchange(GridS *pG)
{
  int kr;
  int i,j,k,n;

  for (k=pG->ke+1; k<=ku; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[kr][j][i].U[2] -= myCoup[k][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[kr][j][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }

  for (k=pG->ke+1; k<=kt; k++) {
    kr = 2*pG->ke-k+1;
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[2] = -myCoup[kr][j][i].U[2];
        for (n=0;n<2; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
        for (n=3;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[kr][j][i].U[n];
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix2_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x2 boundary (ibc_x2=1,5)
 */
static void reflect_ix2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;

  for (k=kb; k<=kt; k++) {
    for (j=jl; j<pG->js; j++) {
      jr = 2*pG->js-j-1;
      for (i=il; i<=iu; i++) {
        myCoup[k][jr][i].U[0] += myCoup[k][j][i].U[0];
        myCoup[k][jr][i].U[1] -= myCoup[k][j][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    } 
  } 
    
  for (k=kb; k<=kt; k++) {
    for (j=jb; j<pG->js; j++) {
      jr = 2*pG->js-j-1;
      for (i=il; i<=iu; i++) { 
        myCoup[k][j][i].U[0] = myCoup[k][jr][i].U[0];
        myCoup[k][j][i].U[1] = -myCoup[k][jr][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][jr][i].U[n];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox2_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x2 boundary (obc_x2=1,5)
 */

static void reflect_ox2_exchange(GridS *pG)
{
  int jr;
  int i,j,k,n;

  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=ju; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        myCoup[k][jr][i].U[0] += myCoup[k][j][i].U[0];
        myCoup[k][jr][i].U[1] -= myCoup[k][j][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][jr][i].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }

  for (k=kb; k<=kt; k++) {
    for (j=pG->je+1; j<=jt; j++) {
      jr = 2*pG->je-j+1;
      for (i=il; i<=iu; i++) {
        myCoup[k][j][i].U[0] = myCoup[k][jr][i].U[0];
        myCoup[k][j][i].U[1] = -myCoup[k][jr][i].U[1];
        for (n=2;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][jr][i].U[n];
      }
    } 
  } 

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ix1_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Inner x1 boundary (ibc_x1=1,5)
 */

static void reflect_ix1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;

  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=jl; i<pG->js; i++) {
        ir = 2*pG->js-i-1;
        myCoup[k][j][ir].U[0] -= myCoup[k][j][i].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }

  for (k=kb; k<=kt; k++) {
    for (j=jb; j<=jt; j++) {
      for (i=jb; i<pG->js; i++) {
        ir = 2*pG->js-i-1;
        myCoup[k][j][i].U[0] = -myCoup[k][j][ir].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][ir].U[n];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void reflect_ox1_exchange(GridS *pG)
 *  \brief REFLECTING boundary conditions, Outer x1 boundary (obc_x1=1,5)
 */

static void reflect_ox1_exchange(GridS *pG)
{
  int ir;
  int i,j,k,n;

  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=iu; i++) {
        ir = 2*pG->ie-i+1;
        myCoup[k][j][ir].U[0] -= myCoup[k][j][i].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][ir].U[n] += myCoup[k][j][i].U[n];
      }
    }
  }
  
  for (k=kb; k<=kt; k++) { 
    for (j=jb; j<=jt; j++) {
      for (i=pG->ie+1; i<=it; i++) {
        ir = 2*pG->ie-i+1;
        myCoup[k][j][i].U[0] = -myCoup[k][j][ir].U[0];
        for (n=1;n<NVar; n++)
          myCoup[k][j][i].U[n] = myCoup[k][j][ir].U[n];
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void outflow_exchange(GridS *pG)
 *  \brief OUTFLOW boundary conditions (ibc=1,5), essentially do nothing
 */

static void outflow_exchange(GridS *pG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix3_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x3 boundary (ibc_x3=4)
 */
static void periodic_ix3_exchange(GridS *pG)
{
  int dk = pG->Nx[2];
  int i,j,k,k1,n;

  for (k=pG->ks; k<pG->ks+NExc; k++) {
   k1 = k+dk-NExc;
   for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
     for (n=0;n<NVar; n++) {
       myCoup[k ][j][i].U[n] += myCoup[k +dk][j][i].U[n];
       myCoup[k1][j][i].U[n] += myCoup[k1-dk][j][i].U[n];
  }}}}

  for (k=pG->ks-NOfst; k<pG->ks; k++) {
   k1 = k+dk+NOfst;
   for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
     for (n=0;n<NVar; n++) {
       myCoup[k ][j][i].U[n] = myCoup[k +dk][j][i].U[n];
       myCoup[k1][j][i].U[n] = myCoup[k1-dk][j][i].U[n];
  }}}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox3_exchange(GridS *pG) 
 *  \brief PERIODIC boundary conditions, Outer x3 boundary (obc_x3=4)
 */

static void periodic_ox3_exchange(GridS *pG)
{
/* Do nothing. All are handled by periodic_ix3_exchange */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix2_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x2 boundary (ibc_x2=4)
 */

static void periodic_ix2_exchange(GridS *pG)
{
  int dj = pG->Nx[1];
  int i,j,j1,k,n;

  for (k=kb; k<=kt; k++) {
   for (j=pG->js; j<pG->js+NExc; j++) {
    j1 = j+dj-NExc;
    for (i=il; i<=iu; i++) {
     for (n=0;n<NVar; n++) {
       myCoup[k][j ][i].U[n] += myCoup[k][j +dj][i].U[n];
       myCoup[k][j1][i].U[n] += myCoup[k][j1-dj][i].U[n];
  }}}}

  for (k=kb; k<=kt; k++) {
   for (j=pG->js-NOfst; j<pG->js; j++) {
    j1 = j+dj+NOfst;
    for (i=il; i<=iu; i++) {
     for (n=0;n<NVar; n++) {
       myCoup[k][j ][i].U[n] = myCoup[k][j +dj][i].U[n];
       myCoup[k][j1][i].U[n] = myCoup[k][j1-dj][i].U[n];
  }}}   }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox2_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions,Outer x2 boundary (obc_x2=4)
 */

static void periodic_ox2_exchange(GridS *pG)
{
/* Do nothing. All are handled by periodic_ix2_exchange */
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ix1_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Inner x1 boundary (ibc_x1=4)
 */

static void periodic_ix1_exchange(GridS *pG)
{
  int di = pG->Nx[0];
  int i,i1,j,k,n;

#ifdef SHEARING_BOX
  if (ShBoxCoord == xy) {
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=0; i<NOfst+NExc; i++) {
       i1 = i+di-NExc;
       for (n=0;n<NVar; n++) {
         myCoup[k][j][pG->is-NOfst +i].U[n] += GhstZns_ox1[k][i][j].U[n];
         myCoup[k][j][pG->ie-NExc+1+i].U[n] += GhstZns_ix1[k][i][j].U[n];
    }}}}

  } else {
#endif
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=pG->is; i<pG->is+NExc; i++) {
       i1 = i+di-NExc;
       for (n=0;n<NVar; n++) {
         myCoup[k][j][i ].U[n] += myCoup[k][j][i +di].U[n];
         myCoup[k][j][i1].U[n] += myCoup[k][j][i1-di].U[n];
    }}}}

    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=pG->is-NOfst; i<pG->is; i++) {
       i1 = i+di+NOfst;
       for (n=0;n<NVar; n++) {
         myCoup[k][j][i ].U[n] = myCoup[k][j][i +di].U[n];
         myCoup[k][j][i1].U[n] = myCoup[k][j][i1-di].U[n];
    }}}}
#ifdef SHEARING_BOX
  }
#endif
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void periodic_ox1_exchange(GridS *pG)
 *  \brief PERIODIC boundary conditions, Outer x1 boundary (obc_x1=4)
 */

static void periodic_ox1_exchange(GridS *pG)
{
/* Do nothing. All are handled by periodic_ix1_exchange */
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 16 funs; ~400 lines */

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix3_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x3 boundary
 */

static void pack_ix3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[0][0]);

  for (k=kl; k<=pG->ks+NOfst-1; k++) {
   for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      for (n=0;n<NVar; n++) {
        *(pd++) = myCoup[k][j][i].U[n];
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox3_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x3 boundary
 */

static void pack_ox3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[1][0]);

  for (k=pG->ke-NOfst+1; k<=ku; k++) {
   for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      for (n=0;n<NVar; n++) {
        *(pd++) = myCoup[k][j][i].U[n];
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix2_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x2 boundary
 */

static void pack_ix2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[0][0]);

  for (k=kb; k<=kt; k++) {
   for (j=jl; j<=pG->js+NOfst-1; j++) {
    for (i=il; i<=iu; i++) {
      for (n=0;n<NVar; n++) {
        *(pd++) = myCoup[k][j][i].U[n];
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox2_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x2 boundary
 */

static void pack_ox2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(send_buf[1][0]);

  for (k=kb; k<=kt; k++) {
   for (j=pG->je-NOfst+1; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      for (n=0;n<NVar; n++) {
        *(pd++) = myCoup[k][j][i].U[n];
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix1_exchange(GridS *pG)
 *  \brief pack the coupling array to bufer at inner x1 boundary
 */

static void pack_ix1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[0][0]);

#ifdef SHEARING_BOX
  if (ShBoxCoord == xy) {
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=0; i<NExc+NOfst; i++) {
       for (n=0;n<NVar; n++) {
         *(pd++) = GhstZns_ix1[k][i][j].U[n];
    }}}}
  } else {
#endif
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=il; i<=pG->is+NOfst-1; i++) {
       for (n=0;n<NVar; n++) {
         *(pd++) = myCoup[k][j][i].U[n];
    }}}}
#ifdef SHEARING_BOX
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox1_exchange(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x1 boundary
 */

static void pack_ox1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[1][0]);

#ifdef SHEARING_BOX
  if (ShBoxCoord == xy) {
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=0; i<NExc+NOfst; i++) {
       for (n=0;n<NVar; n++) {
         *(pd++) = GhstZns_ox1[k][i][j].U[n];
    }}}}
  } else {
#endif
    for (k=kb; k<=kt; k++) {
     for (j=jb; j<=jt; j++) {
      for (i=pG->ie-NOfst+1; i<=iu; i++) {
       for (n=0;n<NVar; n++) {
         *(pd++) = myCoup[k][j][i].U[n];
    }}}}
#ifdef SHEARING_BOX
  }
#endif

  return;
}

#ifdef SHEARING_BOX
/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ix2_remap(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x2 boundary
 */

static void pack_ix2_remap(GridS *pG, GPExc ***myZns)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[0][0]);
	
  for (k=kb; k<=kt; k++) {
   for (i=0; i<NExc+NOfst; i++) {
    for (j=pG->js; j<pG->js+nghost; j++) {
      for (n=0;n<NVar; n++) {
	*(pd++) = myZns[k][i][j].U[n];
      }
  }}}
	
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void pack_ox2_remap(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x2 boundary
 */

static void pack_ox2_remap(GridS *pG, GPExc ***myZns)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(send_buf[1][0]);
	
  for (k=kb; k<=kt; k++) {
   for (i=0; i<NExc+NOfst; i++) {
    for (j=pG->je-nghost+1; j<=pG->je; j++) {
      for (n=0;n<NVar; n++) {
        *(pd++) = myZns[k][i][j].U[n];
      }
  }}}

  return;
}

#endif /* SHEARING_BOX */

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix3_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x3 boundary
 */

static void unpack_ix3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);

  for (k=pG->ks-NOfst; k<pG->ks+NExc; k++){
   for (j=jl; j<=ju; j++){
    for (i=il; i<=iu; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox3_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x3 boundary
 */

static void unpack_ox3_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);

  for (k=pG->ke-NExc+1; k<=pG->ke+NOfst; k++){
   for (j=jl; j<=ju; j++){
    for (i=il; i<=iu; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix2_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x2 boundary
 */

static void unpack_ix2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);

  for (k=kb; k<=kt; k++){
   for (j=pG->js-NOfst; j<pG->js+NExc; j++){
    for (i=il; i<=iu; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox2_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x2 boundary
 */

static void unpack_ox2_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);

  for (k=kb; k<=kt; k++){
   for (j=pG->je-NExc+1; j<=pG->je+NOfst; j++){
    for (i=il; i<=iu; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix1_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at inner x1 boundary
 */

static void unpack_ix1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[0][0]);

  for (k=kb; k<=kt; k++){
   for (j=jb; j<=jt; j++){
    for (i=pG->is-NOfst; i<pG->is+NExc; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox1_exchange(GridS *pG)
 *  \brief unpack the coupling array from buffer at outer x1 boundary
 */

static void unpack_ox1_exchange(GridS *pG)
{
  int i,j,k,n;
  double *pd;
  pd = (double*)&(recv_buf[1][0]);

  for (k=kb; k<=kt; k++){
   for (j=jb; j<=jt; j++){
    for (i=pG->ie-NExc+1; i<=pG->ie+NOfst; i++){
      for (n=0;n<NVar; n++) {
        myCoup[k][j][i].U[n] += *(pd++);
      }
  }}}

  return;
}

#ifdef SHEARING_BOX
/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ix2_remap(GridS *pG)
 *  \brief pack the coupling array to buffer at inner x2 boundary
 */

static void unpack_ix2_remap(GridS *pG, GPExc ***myZns)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(recv_buf[0][0]);
	
  for (k=kb; k<=kt; k++) {
   for (i=0; i<NExc+NOfst; i++) {
    for (j=pG->js-nghost; j<pG->js; j++) {
      for (n=0;n<NVar; n++) {
        myZns[k][i][j].U[n] = *(pd++);
      }
   }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void unpack_ox2_remap(GridS *pG)
 *  \brief pack the coupling array to buffer at outer x2 boundary
 */

static void unpack_ox2_remap(GridS *pG, GPExc ***myZns)
{
  int i,j,k,n;
  double *pd; 
  pd = (double*)&(recv_buf[1][0]);

  for (k=kb; k<=kt; k++) {
   for (i=0; i<NExc+NOfst; i++) {
    for (j=pG->je+1; j<=pG->je+nghost; j++) {
      for (n=0;n<NVar; n++) {
        myZns[k][i][j].U[n] = *(pd++);
      }
  }}}
	
  return;
}
#endif /* SHEARING_BOX */

#endif /* MPI_PARALLEL */

#endif /* PARTICLES */
