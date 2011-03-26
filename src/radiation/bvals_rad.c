#include "../copyright.h"
/*=============================================================================
 * FILE: bvals_rad.c
 *
 * PURPOSE: Sets boundary condistions for radiative transfer on each edge of
 *          the grid.  It closely follows the methods and conventions and
 *          methods used for the hydro integration (see e.g. bvals_mhd.c).
 *
 *          Current version only rudimentary functions which handle specific
 *          problems and is a place-holder for more comprehensive methods
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   bvals_rad_init()
 *   bvals_rad()
 *   bvals_rad_1d_init()
 *   bvals_rad_2d_init() 
 */

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   periodic_??_rad?() - periodic BCs at boundary ???
 *   pack_???_rad()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???_rad()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/
static void periodic_ix1_rad(RadGridS *pRG);
static void periodic_ox1_rad(RadGridS *pRG);
static void periodic_ix2_rad(RadGridS *pRG);
static void periodic_ox2_rad(RadGridS *pRG);
static void periodic_ix3_rad(RadGridS *pRG);
static void periodic_ox3_rad(RadGridS *pRG);

static void ProlongateLater(RadGridS *pRG);
static void const_incident_rad(RadGridS *pRG);

#ifdef MPI_PARALLEL
static void pack_ix1_rad(RadGridS *pRG);
static void pack_ox1_rad(RadGridS *pRG);
static void pack_ix2_rad(RadGridS *pRG);
static void pack_ox2_rad(RadGridS *pRG);
static void pack_ix3_rad(RadGridS *pRG);
static void pack_ox3_rad(RadGridS *pRG);

static void unpack_ix1_rad(RadGridS *pRG);
static void unpack_ox1_rad(RadGridS *pRG);
static void unpack_ix2_rad(RadGridS *pRG);
static void unpack_ox2_rad(RadGridS *pRG);
static void unpack_ix3_rad(RadGridS *pRG);
static void unpack_ox3_rad(RadGridS *pRG);
#endif /* MPI_PARALLEL */

void bvals_rad(DomainS *pD)
{
  RadGridS *pRG=(pD->RadGrid);
#ifdef MPI_PARALLEL
  int cnt, ierr, mIndex;
  int cnt0 = pRG->nf * (1 + pRG->noct * pRG->nang);
  int cnt02 = pRG->nf * (1 + 2 * pRG->noct * pRG->nang);
#endif /* MPI_PARALLEL */
  int l;

  
/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pRG->Nx[0] > 1){
#ifdef MPI_PARALLEL
    cnt = cnt0*(pRG->Nx[1])*(pRG->Nx[2]);
/* MPI blocks to both left and right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1_rad(pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pRG);
      if (mIndex == 1) unpack_ox1_rad(pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pRG);
      if (mIndex == 1) unpack_ox1_rad(pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));
      /* set physical boundary */
      (*(pD->ix1_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_rad(pRG);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx1_id < 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_rad(pRG);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx1_id < 0 && pRG->lx1_id < 0) {
      (*(pD->ix1_RBCFun))(pRG);
      (*(pD->ox1_RBCFun))(pRG);
    }
  }
/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pRG->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt = cnt02*(pRG->Nx[0]+2)*(pRG->Nx[2]);
/* MPI blocks to both left and right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_rad(pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pRG);
      if (mIndex == 1) unpack_ox2_rad(pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pRG);
      if (mIndex == 1) unpack_ox2_rad(pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id < 0) {
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_rad(pRG);


    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx2_id < 0 && pRG->lx2_id >= 0) {
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_rad(pRG);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx2_id < 0 && pRG->lx2_id < 0) {
      (*(pD->ix2_RBCFun))(pRG);
      (*(pD->ox2_RBCFun))(pRG);
    } 
  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pRG->Nx[2] > 1){
#ifdef MPI_PARALLEL
    cnt = cnt0*(pRG->Nx[0]+2)*(pRG->Nx[1]+2);
/* MPI blocks to both left and right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3_rad(pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pRG);
      if (mIndex == 1) unpack_ox3_rad(pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pRG);
      if (mIndex == 1) unpack_ox3_rad(pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_rad(pRG);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx3_id < 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3_rad(pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_RBCFun))(pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_rad(pRG);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pRG->rx3_id < 0 && pRG->lx3_id < 0) {
      (*(pD->ix3_RBCFun))(pRG);
      (*(pD->ox3_RBCFun))(pRG);
    } 

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* bvals_rad_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI.
 *   Patterned closely after bvals_mhd_init().
 */

void bvals_rad_init(MeshS *pM)
{

  RadGridS *pRG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active RadGrids on this proc */

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pRG = pM->Domain[nl][nd].RadGrid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */
#ifdef MPI_PARALLEL
/* get (l,m,n) coordinates of Grid being updated on this processor */
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction -------------*/

    if(pRG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[0] != 0) {      
	pD->ix1_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix1){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ix1_RBCFun = periodic_ix1_rad;
#ifdef MPI_PARALLEL
	  if(pRG->lx1_id < 0 && pD->NGrid[0] > 1){
	    pRG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ix1_RBCFun = const_incident_rad;
	  break;

	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ix1=%d unknown\n",pM->RBCFlag_ix1);
	  exit(EXIT_FAILURE);
	}
      }


/*---- ox1 boundary ----------------------------------------------------------*/


/* Domain boundary is in interior of root */
      if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
	pD->ox1_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
      } else {
	switch(pM->RBCFlag_ox1){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ox1_RBCFun = periodic_ox1_rad;
#ifdef MPI_PARALLEL
	  if(pRG->rx1_id < 0 && pD->NGrid[0] > 1){
	    pRG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ox1_RBCFun = const_incident_rad;
	  break;
	    
	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ox1=%d unknown\n",pM->RBCFlag_ox1);
	  exit(EXIT_FAILURE);
	}
      }
    }

    if(pRG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[1] != 0) {      
	pD->ix2_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix2){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ix2_RBCFun = periodic_ix2_rad;
#ifdef MPI_PARALLEL
	  if(pRG->lx2_id < 0 && pD->NGrid[1] > 1){
	    pRG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ix2_RBCFun = const_incident_rad;
	  break;


	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ix2=%d unknown\n",pM->RBCFlag_ix2);
	  exit(EXIT_FAILURE);
	}
      }


/*---- ox2 boundary ----------------------------------------------------------*/


/* Domain boundary is in interior of root */
      if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
	pD->ox2_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
      } else {
	switch(pM->RBCFlag_ox2){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ox2_RBCFun = periodic_ox2_rad;
#ifdef MPI_PARALLEL
	  if(pRG->rx2_id < 0 && pD->NGrid[1] > 1){
	    pRG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ox2_RBCFun = const_incident_rad;
	  break;

	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ox2=%d unknown\n",pM->RBCFlag_ox2);
	  exit(EXIT_FAILURE);
	}
      }
    }

    if(pRG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[2] != 0) {      
	pD->ix3_RBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix3){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ix3_RBCFun = periodic_ix3_rad;
#ifdef MPI_PARALLEL
	  if(pRG->lx3_id < 0 && pD->NGrid[2] > 1){
	    pRG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ix3_RBCFun = const_incident_rad;
	  break;

	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ix3=%d unknown\n",pM->RBCFlag_ix3);
	  exit(EXIT_FAILURE);
	}
      }


/*---- ox3 boundary ----------------------------------------------------------*/


/* Domain boundary is in interior of root */
      if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
	pD->ox3_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
      } else {
	switch(pM->RBCFlag_ox3){

	case 1: /* Periodic. Handle with MPI calls for parallel jobs. */
	  pD->ox3_RBCFun = periodic_ox3_rad;
#ifdef MPI_PARALLEL
	  if(pRG->rx3_id < 0 && pD->NGrid[2] > 1){
	    pRG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
	  }
#endif /* MPI_PARALLEL */
	  break;

	case 2: /* Open boundary with fixed incoming radiation */
	  pD->ox3_RBCFun = const_incident_rad;
	  break;

	default:
	  ath_perr(-1,"[bvals_rad_init]:rbc_ox3=%d unknown\n",pM->RBCFlag_ox3);
	  exit(EXIT_FAILURE);
	}
      }
    }

/* Figure out largest size needed for send/receive buffers with MPI ----------*/
#ifdef MPI_PARALLEL

    for (n=0; n<(pD->NGrid[2]); n++){
    for (m=0; m<(pD->NGrid[1]); m++){
      for (l=0; l<(pD->NGrid[0]); l++){

/* x1cnt is surface area of x1 faces */
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
	}

/* x2cnt is surface area of x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
	}

/* x3cnt is surface area of x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 2*nghost;

          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
	}
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;


  size *= pRG->nf * (2.0 + 3.0 * pRG->noct * pRG->nang);

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (rbc_ix1=1) */

static void periodic_ix1_rad(RadGridS *pRG)
{

  int il = pRG->is-1, ie = pRG->ie;
  int jl = pRG->js  , ju = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++) {
	pRG->R[k][j][il][ifr].S    = pRG->R[k][j][ie][ifr].S;
	pRG->R[k][j][il][ifr].H[0] = pRG->R[k][j][ie][ifr].H[0];
      }
  /* if (pRG->Nx[1] > 1) { jl--; ju++; } */
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)
	for (m=0; m<nang; m++) {
	  pRG->l1imu[ifr][k][j][0][m] = pRG->r1imu[ifr][k][j][0][m];
	  pRG->l1imu[ifr][k][j][2][m] = pRG->r1imu[ifr][k][j][2][m];
	}

  if (pRG->Nx[1] > 1) {
    for (ifr=0; ifr<nf; ifr++)
      for (k=kl; k<=ku; k++)
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][0][m] = pRG->l2imu[ifr][k][ie][0][m];
	  pRG->l2imu[ifr][k][il][1][m] = pRG->l2imu[ifr][k][ie][1][m];
	  pRG->r2imu[ifr][k][il][2][m] = pRG->r2imu[ifr][k][ie][2][m];
	  pRG->r2imu[ifr][k][il][3][m] = pRG->r2imu[ifr][k][ie][3][m];
	}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (rbc_ox1=1) */

static void periodic_ox1_rad(RadGridS *pRG)
{
  int is = pRG->is, iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++) {
	pRG->R[k][j][iu][ifr].S = pRG->R[k][j][is][ifr].S;
	pRG->R[k][j][iu][ifr].H[0] = pRG->R[k][j][is][ifr].H[0];
      }

  /*if (pRG->Nx[1] > 1) { jl--; ju++; } */
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)
	for(m=0; m<nang; m++) {
	  pRG->r1imu[ifr][k][j][1][m] = pRG->l1imu[ifr][k][j][1][m];
	  pRG->r1imu[ifr][k][j][3][m] = pRG->l1imu[ifr][k][j][3][m];
	}

  if (pRG->Nx[1] > 1)
    for (ifr=0; ifr<nf; ifr++)
      for (k=kl; k<=ku; k++)
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][0][m] = pRG->l2imu[ifr][k][is][0][m];
	  pRG->l2imu[ifr][k][iu][1][m] = pRG->l2imu[ifr][k][is][1][m];
	  pRG->r2imu[ifr][k][iu][2][m] = pRG->r2imu[ifr][k][is][2][m];
	  pRG->r2imu[ifr][k][iu][3][m] = pRG->r2imu[ifr][k][is][3][m];
	}

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=1) */

static void periodic_ix2_rad(RadGridS *pRG)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, je = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][jl][i][ifr].S = pRG->R[k][je][i][ifr].S;

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][i][0][m] = pRG->r2imu[ifr][k][i][0][m];
	  pRG->l2imu[ifr][k][i][1][m] = pRG->r2imu[ifr][k][i][1][m];
	}

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][jl][0][m] = pRG->l1imu[ifr][k][je][0][m];
	pRG->l1imu[ifr][k][jl][2][m] = pRG->l1imu[ifr][k][je][2][m];
	pRG->r1imu[ifr][k][jl][1][m] = pRG->r1imu[ifr][k][je][1][m];
	pRG->r1imu[ifr][k][jl][3][m] = pRG->r1imu[ifr][k][je][3][m];
      }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=1) */

static void periodic_ox2_rad(RadGridS *pRG)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js  , ju = pRG->je+1;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][ju][i][ifr].S = pRG->R[k][js][i][ifr].S;

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  pRG->r2imu[ifr][k][i][2][m] = pRG->l2imu[ifr][k][i][2][m];
	  pRG->r2imu[ifr][k][i][3][m] = pRG->l2imu[ifr][k][i][3][m];
	}

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][ju][0][m] = pRG->l1imu[ifr][k][js][0][m];
	pRG->l1imu[ifr][k][ju][2][m] = pRG->l1imu[ifr][k][js][2][m];
	pRG->r1imu[ifr][k][ju][1][m] = pRG->r1imu[ifr][k][js][1][m];
	pRG->r1imu[ifr][k][ju][3][m] = pRG->r1imu[ifr][k][js][3][m];
      }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (rbc_ix3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ix3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf;
  int i, j, l, m, ifr;

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++) 
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[ks-1][j][i][ifr].S = pRG->R[ke][j][i][ifr].S;

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (rbc_ox3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ox3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf;
  int i, j, l, m, ifr;

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++) 
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[ke+1][j][i][ifr].S = pRG->R[ks][j][i][ifr].S;

  return;
}


/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(RadGridS *pRG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Time independent incident radiation boundary condition.  Nothing is done 
 * here, as the incident boundary radiaion is specified at initialization and
 * unchanged during computation. */

static void const_incident_rad(RadGridS *pRG)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs */
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1_rad(RadGridS *pRG)
{
  int is = pRG->is;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[k][j][is][ifr].S;      

  /*if (pRG->Nx[1] > 1) { jl--; ju++; }*/
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)      
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][1][m];
	    *(pSnd++) = pRG->l1imu[ifr][k][j][3][m];	
	  }

  if (pRG->Nx[1] > 1)
    for (ifr=0; ifr<nf; ifr++)
      for (k=kl; k<=ku; k++)
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][is][0][m];
	  *(pSnd++) = pRG->l2imu[ifr][k][is][1][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][is][2][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][is][3][m];
	}
  

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1_rad(RadGridS *pRG)
{
  int ie = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[k][j][ie][ifr].S;

  /*if (pRG->Nx[1] > 1) { jl--; ju++; } */
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)      
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->r1imu[ifr][k][j][0][m];
	  *(pSnd++) = pRG->r1imu[ifr][k][j][2][m];
	}

  if (pRG->Nx[1] > 1)
   for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {
	*(pSnd++) = pRG->l2imu[ifr][k][ie][0][m];
	*(pSnd++) = pRG->l2imu[ifr][k][ie][1][m];
	*(pSnd++) = pRG->r2imu[ifr][k][ie][2][m];
	*(pSnd++) = pRG->r2imu[ifr][k][ie][3][m];
      }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2_rad(RadGridS *pRG)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[k][js][i][ifr].S;

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][i][2][m];
	  *(pSnd++) = pRG->l2imu[ifr][k][i][3][m];
	}

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][k][js][0][m];
	*(pSnd++) = pRG->l1imu[ifr][k][js][2][m];
	*(pSnd++) = pRG->r1imu[ifr][k][js][1][m];
	*(pSnd++) = pRG->r1imu[ifr][k][js][3][m];
      }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2_rad(RadGridS *pRG)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int je = pRG->je;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[k][je][i][ifr].S;

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->r2imu[ifr][k][i][0][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][i][1][m];
	}

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][k][je][0][m];
	*(pSnd++) = pRG->l1imu[ifr][k][je][2][m];
	*(pSnd++) = pRG->r1imu[ifr][k][je][1][m];
	*(pSnd++) = pRG->r1imu[ifr][k][je][3][m];
      }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf;
  int i, j, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[ks][j][i][ifr].S;

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ke = pRG->ke;
  int nf = pRG->nf;
  int i, j, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++)
      for (ifr=0; ifr<nf; ifr++)
	*(pSnd++) = pRG->R[ke][j][i][ifr].S;

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1_rad(RadGridS *pRG)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][j][il][ifr].S = *(pRcv++);

  /*if (pRG->Nx[1] > 1) { jl--; ju++; }*/
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)
	for(m=0; m<nang; m++) {
	  pRG->l1imu[ifr][k][j][0][m] = *(pRcv++);
	  pRG->l1imu[ifr][k][j][2][m] = *(pRcv++);
	}


  if (pRG->Nx[1] > 1)
    for (ifr=0; ifr<nf; ifr++)
      for (k=kl; k<=ku; k++)
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][0][m] = *(pRcv++);
	  pRG->l2imu[ifr][k][il][1][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][il][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][il][3][m] = *(pRcv++);
	}

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1_rad(RadGridS *pRG)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int j, k, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][j][iu][ifr].S = *(pRcv++);


  /*if (pRG->Nx[1] > 1) { jl--; ju++; }*/
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)
	for(m=0; m<nang; m++) {
	  pRG->r1imu[ifr][k][j][1][m] = *(pRcv++);
	  pRG->r1imu[ifr][k][j][3][m] = *(pRcv++);
	}

  if (pRG->Nx[1] > 1)
    for (ifr=0; ifr<nf; ifr++)
      for (k=kl; k<=ku; k++)
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][0][m] = *(pRcv++);
	  pRG->l2imu[ifr][k][iu][1][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][iu][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][iu][3][m] = *(pRcv++);
	}

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2_rad(RadGridS *pRG)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][jl][i][ifr].S = *(pRcv++);

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][i][0][m] = *(pRcv++);
	  pRG->l2imu[ifr][k][i][1][m] = *(pRcv++);
	}
  
  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][jl][0][m] = *(pRcv++);
	pRG->l1imu[ifr][k][jl][2][m] = *(pRcv++);
	pRG->r1imu[ifr][k][jl][1][m] = *(pRcv++);
	pRG->r1imu[ifr][k][jl][3][m] = *(pRcv++);
      }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2_rad(RadGridS *pRG)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, k, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  for (k=kl; k<=ku; k++)
    for (i=il; i<=iu; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[k][ju][i][ifr].S = *(pRcv++);

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (i=il+1; i<=iu-1; i++)
	for(m=0; m<nang; m++) {
	  pRG->r2imu[ifr][k][i][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][i][3][m] = *(pRcv++);
	}

  for (ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][ju][0][m] = *(pRcv++);
	pRG->l1imu[ifr][k][ju][2][m] = *(pRcv++);
	pRG->r1imu[ifr][k][ju][1][m] = *(pRcv++);
	pRG->r1imu[ifr][k][ju][3][m] = *(pRcv++);
      }

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf;
  int i, j, ifr;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[ks][j][i][ifr].S = *(pRcv++);

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3_rad(RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ke = pRG->ke;
  int nf = pRG->nf;
  int i, j, ifr;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (j=js-1; j<=je+1; j++)
    for (i=is-1; i<=ie+1; i++)
      for (ifr=0; ifr<nf; ifr++)
	pRG->R[ke][j][i][ifr].S = *(pRcv++);

  return;
}
#endif /* MPI_PARALLEL */

#endif /* RADIATION_TRANSFER */
