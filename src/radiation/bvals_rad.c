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
 */

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

static int nDim;
#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
static  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */
#ifdef SHEARING_BOX
extern Real *****GhstZnsIntl, *****GhstZnsIntr;
#endif
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   periodic_??_rad?() - periodic BCs at boundary ???
 *   pack_???_rad()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???_rad()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/
static void periodic_ix1_rad(RadGridS *pRG, int sflag);
static void periodic_ox1_rad(RadGridS *pRG, int sflag);
static void periodic_ix2_rad(RadGridS *pRG, int sflag);
static void periodic_ox2_rad(RadGridS *pRG, int sflag);
static void periodic_ix3_rad(RadGridS *pRG, int sflag);
static void periodic_ox3_rad(RadGridS *pRG, int sflag);

static void ProlongateLater(RadGridS *pRG, int sflag);
static void const_incident_rad(RadGridS *pRG, int sflag);

#ifdef MPI_PARALLEL
static void pack_ix1_rad(RadGridS *pRG, int sflag);
static void pack_ox1_rad(RadGridS *pRG, int sflag);
static void pack_ix2_rad(RadGridS *pRG, int sflag);
static void pack_ox2_rad(RadGridS *pRG, int sflag);
static void pack_ix3_rad(RadGridS *pRG, int sflag);
static void pack_ox3_rad(RadGridS *pRG, int sflag);

static void unpack_ix1_rad(RadGridS *pRG, int sflag);
static void unpack_ox1_rad(RadGridS *pRG, int sflag);
static void unpack_ix2_rad(RadGridS *pRG, int sflag);
static void unpack_ox2_rad(RadGridS *pRG, int sflag);
static void unpack_ix3_rad(RadGridS *pRG, int sflag);
static void unpack_ox3_rad(RadGridS *pRG, int sflag);
#endif /* MPI_PARALLEL */

void bvals_rad(DomainS *pD, int sflag)
{
  RadGridS *pRG=(pD->RadGrid);
#ifdef SHEARING_BOX
  int myL,myM,myN,BCFlag;
#endif
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
    cnt = x1cnt;
/* MPI blocks to both left and right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1_rad(pRG,sflag);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox1_rad(pRG,sflag);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox1_rad(pRG,sflag);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));
      /* set physical boundary */
      (*(pD->ix1_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_rad(pRG,sflag);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx1_id < 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_rad(pRG,sflag);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx1_id < 0 && pRG->lx1_id < 0) {
      (*(pD->ix1_RBCFun))(pRG,sflag);
      (*(pD->ox1_RBCFun))(pRG,sflag);
    }
  }
/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pRG->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt = x2cnt;
/* MPI blocks to both left and right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2_rad(pRG,sflag);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox2_rad(pRG,sflag);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox2_rad(pRG,sflag);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id < 0) {
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_rad(pRG,sflag);


    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx2_id < 0 && pRG->lx2_id >= 0) {
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_rad(pRG,sflag);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx2_id < 0 && pRG->lx2_id < 0) {
      (*(pD->ix2_RBCFun))(pRG,sflag);
      (*(pD->ox2_RBCFun))(pRG,sflag);
    } 
 
/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */

#ifdef SHEARING_BOX 
    BCFlag = par_geti_def("domain1","rbc_ix1",0);
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0 && BCFlag == 4) {
      ShearingSheet_Rad_ix1(pD);
      } 
    BCFlag = par_geti_def("domain1","rbc_ox1",0);
    if (myL == ((pD->NGrid[0])-1) && BCFlag == 4) {
      ShearingSheet_Rad_ox1(pD);
      }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pRG->Nx[2] > 1){
#ifdef MPI_PARALLEL
    cnt = x3cnt;
/* MPI blocks to both left and right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3_rad(pRG,sflag);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox3_rad(pRG,sflag);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pRG,sflag);
      if (mIndex == 1) unpack_ox3_rad(pRG,sflag);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_rad(pRG,sflag);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx3_id < 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3_rad(pRG,sflag); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_RBCFun))(pRG,sflag);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_rad(pRG,sflag);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pRG->rx3_id < 0 && pRG->lx3_id < 0) {
      (*(pD->ix3_RBCFun))(pRG,sflag);
      (*(pD->ox3_RBCFun))(pRG,sflag);
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
  /*int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
  int nang, nf, noct, xcnt;
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
      nDim = 1;
/*---- ix1 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[0] != 0) {      
	pD->ix1_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix1){

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

      nDim = 2;
/*---- ix2 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[1] != 0) {      
	pD->ix2_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix2){

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

      nDim = 3;
/*---- ix3 boundary ----------------------------------------------------------*/

/* Domain boundary is in interior of root */
      if(pD->Disp[2] != 0) {      
	pD->ix3_RBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
      } else {                    
	switch(pM->RBCFlag_ix3){

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

	case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
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

    nang = pRG->nang;
    noct = pRG->noct;
    nf = pRG->nf;
    for (n=0; n<(pD->NGrid[2]); n++){
    for (m=0; m<(pD->NGrid[1]); m++){
      for (l=0; l<(pD->NGrid[0]); l++){

/* x1cnt is surface area of x1 faces */
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  nx3t = pD->GData[n][m][l].Nx[2];
	  xcnt =  nx2t * nx3t;
#ifdef SHEARING_BOX
	  xcnt += nx2t * nx3t * 2;
#endif
#ifdef JACOBI
	  if (noct == 4) xcnt += 2 * nx3t;
	  if (noct == 8) xcnt += 2 * nx2t;
	  xcnt *= noct * nang;
#else
	  if (noct == 4) xcnt += nx3t;
	  if (noct == 8) xcnt += nx2t;
	  xcnt *= noct * nang / 2;
#endif
	  xcnt += (nDim + 1) * nx2t * nx3t;	  
	  xcnt *= nf;
          if(xcnt > x1cnt) x1cnt = xcnt;
	}

/* x2cnt is surface area of x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2;
	  nx3t = pD->GData[n][m][l].Nx[2];
	  xcnt = nx3t * (nx1t + 2);
#ifdef SHEARING_BOX
	  xcnt +=  nx3t * 4; 
#endif
#ifdef JACOBI
	  if (noct == 8) xcnt += 2 * nx1t;
	  xcnt *= noct * nang;
#else
	  if (noct == 8) xcnt += nx1t;
	  xcnt *= noct * nang / 2;
#endif
	  xcnt += (nDim + 1) * nx3t * nx1t;
	  xcnt *= nf;
          if(xcnt > x2cnt) x2cnt = xcnt;
	}

/* x3cnt is surface area of x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2;
	  nx2t = pD->GData[n][m][l].Nx[1] + 2;
	  xcnt = nx1t * nx2t + 2 * (nx1t + nx2t);
#ifdef JACOBI
	  xcnt *= noct * nang;
#else
	  xcnt *= noct * nang / 2;
#endif
	  xcnt += (nDim + 1) * nx2t * nx1t;
	  xcnt *= nf;
          if(xcnt > x3cnt) x3cnt = xcnt;
	}
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  printf("size of buffers: %d %d %d %d\n",size,x1cnt,x2cnt,x3cnt);

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

#ifdef JACOBI
/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (rbc_ix1=1) */

static void periodic_ix1_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, ie = pRG->ie;
  int jl = pRG->js  , ju = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][il][ifr].S = pRG->R[k][j][ie][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][il][ifr].H[l] = pRG->R[k][j][ie][ifr].H[l];
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][j][l][m] = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++){
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    GhstZnsIntl[ifr][k][j][l][m] = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][il][l][m] = pRG->Ghstl1i[ifr][k][jl][l][m];
	    pRG->r2imu[ifr][k][il][l][m] = pRG->Ghstl1i[ifr][k][ju][l][m];
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][il][l][m] = pRG->Ghstl1i[ifr][kl][j][l][m];
	    pRG->r3imu[ifr][j][il][l][m] = pRG->Ghstl1i[ifr][ku][j][l][m];
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (rbc_ox1=1) */

static void periodic_ox1_rad(RadGridS *pRG, int sflag)
{
  int is = pRG->is, iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][iu][ifr].S = pRG->R[k][j][is][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][iu][ifr].H[l] = pRG->R[k][j][is][ifr].H[l];
	  }
	}}}
  }
  /* update r1imu using l1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][l][m] = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++){
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    GhstZnsIntr[ifr][k][j][l][m] = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#endif
 /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][iu][l][m] = pRG->Ghstr1i[ifr][k][jl][l][m];
	    pRG->r2imu[ifr][k][iu][l][m] = pRG->Ghstr1i[ifr][k][ju][l][m];
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][iu][l][m] = pRG->Ghstr1i[ifr][kl][j][l][m];
	    pRG->r3imu[ifr][j][iu][l][m] = pRG->Ghstr1i[ifr][ku][j][l][m];
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=1) */

static void periodic_ix2_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, je = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][jl][i][ifr].S = pRG->R[k][je][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][jl][i][ifr].H[l] = pRG->R[k][je][i][ifr].H[l];
	  }
	}}}
  }
  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones.  r2imu on the
   * corners/edges has already been updated (if necessary) by
   * x1 boundary routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
       for (i=il; i<=iu; i++) {
	 for (l=0; l<noct; l++) {
	   for(m=0; m<nang; m++) {
	     pRG->Ghstl2i[ifr][k][i][l][m] = pRG->r2imu[ifr][k][i][l][m];
	   }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][jl][l][m] = GhstZnsIntl[ifr][k][je][l][m];
	  GhstZnsIntr[ifr][k][jl][l][m] = GhstZnsIntr[ifr][k][je][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][jl][l][m] = pRG->Ghstl1i[ifr][k][je][l][m];
	  pRG->Ghstr1i[ifr][k][jl][l][m] = pRG->Ghstr1i[ifr][k][je][l][m];
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][jl][i][l][m] = pRG->Ghstl2i[ifr][kl][i][l][m];
	    pRG->r3imu[ifr][jl][i][l][m] = pRG->Ghstl2i[ifr][ku][i][l][m];
	  }}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=1) */

static void periodic_ox2_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js  , ju = pRG->je+1;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][ju][i][ifr].S = pRG->R[k][js][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][ju][i][ifr].H[l] = pRG->R[k][js][i][ifr].H[l];
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones.  l2imu on the
   * corners/edges has already been updated (if necessary) by
   * x1 boundary routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr2i[ifr][k][i][l][m] = pRG->l2imu[ifr][k][i][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][ju][l][m] = GhstZnsIntl[ifr][k][js][l][m];
	  GhstZnsIntr[ifr][k][ju][l][m] = GhstZnsIntr[ifr][k][js][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][ju][l][m] = pRG->Ghstl1i[ifr][k][js][l][m];
	  pRG->Ghstr1i[ifr][k][ju][l][m] = pRG->Ghstr1i[ifr][k][js][l][m];
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][ju][i][l][m] = pRG->Ghstr2i[ifr][kl][i][l][m];
	    pRG->r3imu[ifr][ju][i][l][m] = pRG->Ghstr2i[ifr][ku][i][l][m];
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (rbc_ix3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ix3_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[kl][j][i][ifr].S = pRG->R[ke][j][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[kl][j][i][ifr].H[l] = pRG->R[ke][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  r3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = pRG->r3imu[ifr][j][i][l][m];
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][kl][j][l][m] = pRG->Ghstl1i[ifr][ke][j][l][m];
	  pRG->Ghstr1i[ifr][kl][j][l][m] = pRG->Ghstr1i[ifr][ke][j][l][m];
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][kl][i][l][m] = pRG->Ghstl2i[ifr][ke][i][l][m];
	  pRG->Ghstr2i[ifr][kl][i][l][m] = pRG->Ghstr2i[ifr][ke][i][l][m];
	}}}}

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (rbc_ox3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks, ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[ku][j][i][ifr].S = pRG->R[ks][j][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[ku][j][i][ifr].H[l] = pRG->R[ks][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  l3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = pRG->l3imu[ifr][j][i][l][m];
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][ku][j][l][m] = pRG->Ghstl1i[ifr][ks][j][l][m];
	  pRG->Ghstr1i[ifr][ku][j][l][m] = pRG->Ghstr1i[ifr][ks][j][l][m];
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][ku][i][l][m] = pRG->Ghstl2i[ifr][ks][i][l][m];
	  pRG->Ghstr2i[ifr][ku][i][l][m] = pRG->Ghstr2i[ifr][ks][i][l][m];
	}}}}

  return;
}
#else
/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (rbc_ix1=1) */

static void periodic_ix1_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, ie = pRG->ie;
  int jl = pRG->js  , ju = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][il][ifr].S = pRG->R[k][j][ie][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][il][ifr].H[l] = pRG->R[k][j][ie][ifr].H[l];
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l1imu[ifr][k][j][0][m] = pRG->r1imu[ifr][k][j][0][m];
	  if(noct > 2) {
	    pRG->l1imu[ifr][k][j][2][m] = pRG->r1imu[ifr][k][j][2][m];
	    if(noct == 8) {
	      pRG->l1imu[ifr][k][j][4][m] = pRG->r1imu[ifr][k][j][4][m];
	      pRG->l1imu[ifr][k][j][6][m] = pRG->r1imu[ifr][k][j][6][m];
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++){
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    GhstZnsIntl[ifr][k][j][l][m] = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=0,2,4,6 are passed using r1imu */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][2][m] = pRG->r1imu[ifr][k][jl][2][m];
	  pRG->l2imu[ifr][k][il][3][m] = pRG->l2imu[ifr][k][ie][3][m];
	  pRG->r2imu[ifr][k][il][0][m] = pRG->r1imu[ifr][k][ju][0][m];
	  pRG->r2imu[ifr][k][il][1][m] = pRG->r2imu[ifr][k][ie][1][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][il][6][m] = pRG->r1imu[ifr][k][jl][6][m];
	      pRG->l2imu[ifr][k][il][7][m] = pRG->l2imu[ifr][k][ie][7][m];
	      pRG->r2imu[ifr][k][il][4][m] = pRG->r1imu[ifr][k][ju][4][m];
	      pRG->r2imu[ifr][k][il][5][m] = pRG->r2imu[ifr][k][ie][5][m];
	    }
	}}}
  }
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,2,4,6 are passed using r1imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][il][4][m] = pRG->r1imu[ifr][kl][j ][4][m];
	  pRG->l3imu[ifr][j][il][5][m] = pRG->l3imu[ifr][j ][ie][5][m];
	  pRG->l3imu[ifr][j][il][6][m] = pRG->r1imu[ifr][kl][j ][6][m];
	  pRG->l3imu[ifr][j][il][7][m] = pRG->l3imu[ifr][j ][ie][7][m];
	  pRG->r3imu[ifr][j][il][0][m] = pRG->r1imu[ifr][ku][j ][0][m];
	  pRG->r3imu[ifr][j][il][1][m] = pRG->r3imu[ifr][j ][ie][1][m];
	  pRG->r3imu[ifr][j][il][2][m] = pRG->r1imu[ifr][ku][j ][2][m];
	  pRG->r3imu[ifr][j][il][3][m] = pRG->r3imu[ifr][j ][ie][3][m];	   
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (rbc_ox1=1) */

static void periodic_ox1_rad(RadGridS *pRG, int sflag)
{
  int is = pRG->is, iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][iu][ifr].S = pRG->R[k][j][is][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][iu][ifr].H[l] = pRG->R[k][j][is][ifr].H[l];
	  }
	}}}
  }
  /* update r1imu using l1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(m=0; m<nang; m++) {
	  pRG->r1imu[ifr][k][j][1][m] = pRG->l1imu[ifr][k][j][1][m];
	  if(noct > 2) {
	    pRG->r1imu[ifr][k][j][3][m] = pRG->l1imu[ifr][k][j][3][m];
	    if(noct == 8) {
	      pRG->r1imu[ifr][k][j][5][m] = pRG->l1imu[ifr][k][j][5][m];
	      pRG->r1imu[ifr][k][j][7][m] = pRG->l1imu[ifr][k][j][7][m];
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++){
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    GhstZnsIntr[ifr][k][j][l][m] = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=1,3,5,7 are passed using l1imu */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][2][m] = pRG->l2imu[ifr][k][is][2][m];
	  pRG->l2imu[ifr][k][iu][3][m] = pRG->l1imu[ifr][k][jl][3][m];
	  pRG->r2imu[ifr][k][iu][0][m] = pRG->r2imu[ifr][k][is][0][m];
	  pRG->r2imu[ifr][k][iu][1][m] = pRG->l1imu[ifr][k][ju][1][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][iu][6][m] = pRG->l2imu[ifr][k][is][6][m];
	      pRG->l2imu[ifr][k][iu][7][m] = pRG->l1imu[ifr][k][jl][7][m];
	      pRG->r2imu[ifr][k][iu][4][m] = pRG->r2imu[ifr][k][is][4][m];
	      pRG->r2imu[ifr][k][iu][5][m] = pRG->l1imu[ifr][k][ju][5][m];
	    }
	}}}
  }
 /* pass l3imu/r3imu on the edges.  Note that values
   * l=1,3,5,7 are passed using r1imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][iu][4][m] = pRG->l3imu[ifr][j ][is][4][m];
	  pRG->l3imu[ifr][j][iu][5][m] = pRG->l1imu[ifr][kl][j ][5][m];
	  pRG->l3imu[ifr][j][iu][6][m] = pRG->l3imu[ifr][j ][is][6][m];
	  pRG->l3imu[ifr][j][iu][7][m] = pRG->l1imu[ifr][kl][j ][7][m];
	  pRG->r3imu[ifr][j][iu][0][m] = pRG->r3imu[ifr][j ][is][0][m];
	  pRG->r3imu[ifr][j][iu][1][m] = pRG->l1imu[ifr][ku][j ][1][m];
	  pRG->r3imu[ifr][j][iu][2][m] = pRG->r3imu[ifr][j ][is][2][m];	  
	  pRG->r3imu[ifr][j][iu][3][m] = pRG->l1imu[ifr][ku][j ][3][m];
	}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=1) */

static void periodic_ix2_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, je = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][jl][i][ifr].S = pRG->R[k][je][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][jl][i][ifr].H[l] = pRG->R[k][je][i][ifr].H[l];
	  }
	}}}
  }
  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones.  r2imu on the
   * corners/edges has already been updated (if necessary) by
   * x1 boundary routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][i][0][m] = pRG->r2imu[ifr][k][i][0][m];
	  pRG->l2imu[ifr][k][i][1][m] = pRG->r2imu[ifr][k][i][1][m];
	  if (noct == 8) {
	    pRG->l2imu[ifr][k][i][4][m] = pRG->r2imu[ifr][k][i][4][m];
	    pRG->l2imu[ifr][k][i][5][m] = pRG->r2imu[ifr][k][i][5][m];	    
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][jl][l][m] = GhstZnsIntl[ifr][k][je][l][m];
	  GhstZnsIntr[ifr][k][jl][l][m] = GhstZnsIntr[ifr][k][je][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][jl][0][m] = pRG->l1imu[ifr][k][je][0][m];
	pRG->l1imu[ifr][k][jl][2][m] = pRG->l1imu[ifr][k][je][2][m];
	pRG->r1imu[ifr][k][jl][1][m] = pRG->r1imu[ifr][k][je][1][m];
	pRG->r1imu[ifr][k][jl][3][m] = pRG->r1imu[ifr][k][je][3][m];
	if (noct == 8) {
	  pRG->l1imu[ifr][k][jl][4][m] = pRG->l1imu[ifr][k][je][4][m];
	  pRG->l1imu[ifr][k][jl][6][m] = pRG->l1imu[ifr][k][je][6][m];
	  pRG->r1imu[ifr][k][jl][5][m] = pRG->r1imu[ifr][k][je][5][m];
	  pRG->r1imu[ifr][k][jl][7][m] = pRG->r1imu[ifr][k][je][7][m];
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,1,4,5 are passed using r2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][jl][i][4][m] = pRG->r2imu[ifr][kl][i][4][m];
	  pRG->l3imu[ifr][jl][i][5][m] = pRG->r2imu[ifr][kl][i][5][m];
	  pRG->l3imu[ifr][jl][i][6][m] = pRG->l3imu[ifr][je][i][6][m];
	  pRG->l3imu[ifr][jl][i][7][m] = pRG->l3imu[ifr][je][i][7][m];
	  pRG->r3imu[ifr][jl][i][0][m] = pRG->r2imu[ifr][ku][i][0][m];
	  pRG->r3imu[ifr][jl][i][1][m] = pRG->r2imu[ifr][ku][i][1][m];
	  pRG->r3imu[ifr][jl][i][2][m] = pRG->r3imu[ifr][je][i][2][m];
	  pRG->r3imu[ifr][jl][i][3][m] = pRG->r3imu[ifr][je][i][3][m];
	}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=1) */

static void periodic_ox2_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js  , ju = pRG->je+1;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][ju][i][ifr].S = pRG->R[k][js][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[k][ju][i][ifr].H[l] = pRG->R[k][js][i][ifr].H[l];
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones.  l2imu on the
   * corners/edges has already been updated (if necessary) by
   * x1 boundary routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->r2imu[ifr][k][i][2][m] = pRG->l2imu[ifr][k][i][2][m];
	  pRG->r2imu[ifr][k][i][3][m] = pRG->l2imu[ifr][k][i][3][m];
	  if (noct == 8) {
	    pRG->r2imu[ifr][k][i][6][m] = pRG->l2imu[ifr][k][i][6][m];
	    pRG->r2imu[ifr][k][i][7][m] = pRG->l2imu[ifr][k][i][7][m];	    
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][ju][l][m] = GhstZnsIntl[ifr][k][js][l][m];
	  GhstZnsIntr[ifr][k][ju][l][m] = GhstZnsIntr[ifr][k][js][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][ju][0][m] = pRG->l1imu[ifr][k][js][0][m];
	pRG->l1imu[ifr][k][ju][2][m] = pRG->l1imu[ifr][k][js][2][m];
	pRG->r1imu[ifr][k][ju][1][m] = pRG->r1imu[ifr][k][js][1][m];
	pRG->r1imu[ifr][k][ju][3][m] = pRG->r1imu[ifr][k][js][3][m];
	if (noct == 8) {
	  pRG->l1imu[ifr][k][ju][4][m] = pRG->l1imu[ifr][k][js][4][m];
	  pRG->l1imu[ifr][k][ju][6][m] = pRG->l1imu[ifr][k][js][6][m];
	  pRG->r1imu[ifr][k][ju][5][m] = pRG->r1imu[ifr][k][js][5][m];
	  pRG->r1imu[ifr][k][ju][7][m] = pRG->r1imu[ifr][k][js][7][m];
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=2,3,6,7 are passed using l2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][ju][i][4][m] = pRG->l3imu[ifr][js][i][4][m];
	  pRG->l3imu[ifr][ju][i][5][m] = pRG->l3imu[ifr][js][i][5][m];
	  pRG->l3imu[ifr][ju][i][6][m] = pRG->l2imu[ifr][kl][i][6][m];
	  pRG->l3imu[ifr][ju][i][7][m] = pRG->l2imu[ifr][kl][i][7][m];
	  pRG->r3imu[ifr][ju][i][0][m] = pRG->r3imu[ifr][js][i][0][m];
	  pRG->r3imu[ifr][ju][i][1][m] = pRG->r3imu[ifr][js][i][1][m];
	  pRG->r3imu[ifr][ju][i][2][m] = pRG->l2imu[ifr][ku][i][2][m];
	  pRG->r3imu[ifr][ju][i][3][m] = pRG->l2imu[ifr][ku][i][3][m];
	}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (rbc_ix3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ix3_rad(RadGridS *pRG, int sflag)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[kl][j][i][ifr].S = pRG->R[ke][j][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[kl][j][i][ifr].H[l] = pRG->R[ke][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  r3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][i][l][m] = pRG->r3imu[ifr][j][i][l][m];
	  }}}}}
 /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][kl][j][0][m] = pRG->l1imu[ifr][ke][j][0][m];
	pRG->l1imu[ifr][kl][j][2][m] = pRG->l1imu[ifr][ke][j][2][m];
	pRG->l1imu[ifr][kl][j][4][m] = pRG->l1imu[ifr][ke][j][4][m];
	pRG->l1imu[ifr][kl][j][6][m] = pRG->l1imu[ifr][ke][j][6][m];
	pRG->r1imu[ifr][kl][j][1][m] = pRG->r1imu[ifr][ke][j][1][m];
	pRG->r1imu[ifr][kl][j][3][m] = pRG->r1imu[ifr][ke][j][3][m];
	pRG->r1imu[ifr][kl][j][5][m] = pRG->r1imu[ifr][ke][j][5][m];
	pRG->r1imu[ifr][kl][j][7][m] = pRG->r1imu[ifr][ke][j][7][m];
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->l2imu[ifr][kl][i][0][m] = pRG->l2imu[ifr][ke][i][0][m];
	pRG->l2imu[ifr][kl][i][1][m] = pRG->l2imu[ifr][ke][i][1][m];
	pRG->l2imu[ifr][kl][i][4][m] = pRG->l2imu[ifr][ke][i][4][m];
	pRG->l2imu[ifr][kl][i][5][m] = pRG->l2imu[ifr][ke][i][5][m];
	pRG->r2imu[ifr][kl][i][2][m] = pRG->r2imu[ifr][ke][i][2][m];
	pRG->r2imu[ifr][kl][i][3][m] = pRG->r2imu[ifr][ke][i][3][m];
	pRG->r2imu[ifr][kl][i][6][m] = pRG->r2imu[ifr][ke][i][6][m];
	pRG->r2imu[ifr][kl][i][7][m] = pRG->r2imu[ifr][ke][i][7][m];
      }}}

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (rbc_ox3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks, ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[ku][j][i][ifr].S = pRG->R[ks][j][i][ifr].S;
	  for(l=0; l < nDim; l++) {
	    pRG->R[ku][j][i][ifr].H[l] = pRG->R[ks][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  l3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->r3imu[ifr][j][i][l][m] = pRG->l3imu[ifr][j][i][l][m];
	  }}}}}
/* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][ku][j][0][m] = pRG->l1imu[ifr][ks][j][0][m];
	pRG->l1imu[ifr][ku][j][2][m] = pRG->l1imu[ifr][ks][j][2][m];
	pRG->l1imu[ifr][ku][j][4][m] = pRG->l1imu[ifr][ks][j][4][m];
	pRG->l1imu[ifr][ku][j][6][m] = pRG->l1imu[ifr][ks][j][6][m];
	pRG->r1imu[ifr][ku][j][1][m] = pRG->r1imu[ifr][ks][j][1][m];
	pRG->r1imu[ifr][ku][j][3][m] = pRG->r1imu[ifr][ks][j][3][m];
	pRG->r1imu[ifr][ku][j][5][m] = pRG->r1imu[ifr][ks][j][5][m];
	pRG->r1imu[ifr][ku][j][7][m] = pRG->r1imu[ifr][ks][j][7][m];
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->l2imu[ifr][ku][i][0][m] = pRG->l2imu[ifr][ks][i][0][m];
	pRG->l2imu[ifr][ku][i][1][m] = pRG->l2imu[ifr][ks][i][1][m];
	pRG->l2imu[ifr][ku][i][4][m] = pRG->l2imu[ifr][ks][i][4][m];
	pRG->l2imu[ifr][ku][i][5][m] = pRG->l2imu[ifr][ks][i][5][m];
	pRG->r2imu[ifr][ku][i][2][m] = pRG->r2imu[ifr][ks][i][2][m];
	pRG->r2imu[ifr][ku][i][3][m] = pRG->r2imu[ifr][ks][i][3][m];
	pRG->r2imu[ifr][ku][i][6][m] = pRG->r2imu[ifr][ks][i][6][m];
	pRG->r2imu[ifr][ku][i][7][m] = pRG->r2imu[ifr][ks][i][7][m];
      }}}


  return;
}
#endif /* JACOBI */

/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(RadGridS *pRG, int sflag)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Time independent incident radiation boundary condition.  Nothing is done 
 * here, as the incident boundary radiaion is specified at initialization and
 * unchanged during computation. */

static void const_incident_rad(RadGridS *pRG, int sflag)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs */
#ifdef JACOBI
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1_rad(RadGridS *pRG, int sflag)
{
  int is = pRG->is;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;
  
  pSnd = (double*)&(send_buf[0][0]);
  /* if LTE send the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][j][is][ifr].S;  
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][j][is][ifr].H[l];
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#endif
 /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][jl][l][m];
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][ju][l][m];
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][kl][j][l][m];
	    *(pSnd++) = pRG->Ghstr1i[ifr][ku][j][l][m];
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1_rad(RadGridS *pRG, int sflag)
{
  int ie = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);
  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][j][ie][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][j][ie][ifr].H[l];
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][jl][l][m];
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][ju][l][m];
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][kl][j][l][m];
	    *(pSnd++) = pRG->Ghstl1i[ifr][ku][j][l][m];
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][js][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][js][i][ifr].H[l];
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l2imu[ifr][k][i][l][m];
	  }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = GhstZnsIntl[ifr][k][js][l][m];
	  *(pSnd++) = GhstZnsIntr[ifr][k][js][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][js][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][js][l][m];
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr2i[ifr][kl][i][l][m];
	    *(pSnd++) = pRG->Ghstr2i[ifr][ku][i][l][m];
	  }}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int je = pRG->je;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][je][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][je][i][ifr].H[l];
	  }
	}}}
  }
  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
       for (i=il; i<=iu; i++) {
	 for (l=0; l<noct; l++) {
	   for(m=0; m<nang; m++) {
	     *(pSnd++) = pRG->r2imu[ifr][k][i][l][m];
	   }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = GhstZnsIntl[ifr][k][je][l][m];
	  *(pSnd++) = GhstZnsIntr[ifr][k][je][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][je][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][je][l][m];
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl2i[ifr][kl][i][l][m];
	    *(pSnd++) = pRG->Ghstl2i[ifr][ku][i][l][m];
	  }}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[ks][j][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[ks][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l3imu[ifr][j][i][l][m];
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][ks][j][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][ks][j][l][m];
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl2i[ifr][ks][i][l][m];
	  *(pSnd++) = pRG->Ghstr2i[ifr][ks][i][l][m];
	}}}}

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[ke][j][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[ke][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r3imu[ifr][j][i][l][m];
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][ke][j][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][ke][j][l][m];
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl2i[ifr][ke][i][l][m];
	  *(pSnd++) = pRG->Ghstr2i[ifr][ke][i][l][m];
	}}}}

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);
  /* if LTE receive the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][il][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][il][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][j][l][m] = *(pRcv++);
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    GhstZnsIntl[ifr][k][j][l][m]= *(pRcv++);
	  }}}}}
#endif
 /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][il][l][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][il][l][m] = *(pRcv++);
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][il][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][j][il][l][m] = *(pRcv++);
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1_rad(RadGridS *pRG, int sflag)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);
  /* if LTE receive the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][iu][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][iu][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update r1imu using l1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][l][m] = *(pRcv++);
	  }}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    GhstZnsIntr[ifr][k][j][l][m]= *(pRcv++);
	  }}}}}
#endif
 /* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][iu][l][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][iu][l][m] = *(pRcv++);
	  }}}}
  }
  /* pass l3imu/r3imu on the edges.*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][iu][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][j][iu][l][m] = *(pRcv++);
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][jl][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][jl][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
       for (i=il; i<=iu; i++) {
	 for (l=0; l<noct; l++) {
	   for(m=0; m<nang; m++) {
	     pRG->Ghstl2i[ifr][k][i][l][m] = *(pRcv++);
	   }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][jl][l][m] = *(pRcv++);
	  GhstZnsIntr[ifr][k][jl][l][m] = *(pRcv++);
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][jl][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][jl][l][m] = *(pRcv++);
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][jl][i][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][jl][i][l][m] = *(pRcv++);
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][ju][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][ju][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr2i[ifr][k][i][l][m] = *(pRcv++);
	  }}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][ju][l][m] = *(pRcv++);
	  GhstZnsIntr[ifr][k][ju][l][m] = *(pRcv++);
	}}}} 
#endif
   /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][ju][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][ju][l][m] = *(pRcv++);
	}}}}
  /* pass l3imu/r3imu on the edges. */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][ju][i][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][ju][i][l][m] = *(pRcv++);
	  }}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[kl][j][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[kl][j][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  r3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = *(pRcv++);
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][kl][j][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][kl][j][l][m] = *(pRcv++);
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][kl][i][l][m] = *(pRcv++);
	  pRG->Ghstr2i[ifr][kl][i][l][m] = *(pRcv++);
	}}}}
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[ku][j][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[ku][j][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.  l3imu on the
   * edges has already been updated (if necessary) by x1 and x2 boundary
   * routines */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = *(pRcv++);
	  }}}}}
  /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][ku][j][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][ku][j][l][m] = *(pRcv++);
	}}}}
  /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][ku][i][l][m] = *(pRcv++);
	  pRG->Ghstr2i[ifr][ku][i][l][m] = *(pRcv++);
	}}}}
  return;
}
#else
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1_rad(RadGridS *pRG, int sflag)
{
  int is = pRG->is;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  /* if LTE send the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][j][is][ifr].S;  
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][j][is][ifr].H[l];
	  }
	}}}
  }
 /* update r1imu using l1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l1imu[ifr][k][j][1][m];
	  if(noct > 2) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][3][m];
	    if(noct == 8) {
	      *(pSnd++) = pRG->l1imu[ifr][k][j][5][m];
	      *(pSnd++) = pRG->l1imu[ifr][k][j][7][m];
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=1,3,5,7 are passed using l1imu and set in unpack routine */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][is][2][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][is][0][m];
	  if(noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][is][6][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][is][4][m];
	  }
	}}}
  }
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=1,3,5,7 are passed using r1imu and set in unpack routine */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][j][is][4][m];
	  *(pSnd++) = pRG->l3imu[ifr][j][is][6][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][is][0][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][is][2][m];	  
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1_rad(RadGridS *pRG, int sflag)
{
  int ie = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);
  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][j][ie][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][j][ie][ifr].H[l];
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++)  = pRG->r1imu[ifr][k][j][0][m];
	  if(noct > 2) {
	    *(pSnd++)  = pRG->r1imu[ifr][k][j][2][m];
	    if(noct == 8) {
	      *(pSnd++)  = pRG->r1imu[ifr][k][j][4][m];
	      *(pSnd++)  = pRG->r1imu[ifr][k][j][6][m];
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r1imu[ifr][k][j][l][m];
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=0,2,4,6 are passed using r1imu and set in unpack routine */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][ie][3][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][ie][1][m];
	  if(noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][ie][7][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][ie][5][m];
	  }
	}}}
  }
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,2,4,6 are passed using r1imu and set in unpack routine */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][j][ie][5][m];
	  *(pSnd++) = pRG->l3imu[ifr][j][ie][7][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][ie][1][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][ie][3][m];	   
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][js][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][js][i][ifr].H[l];
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][i][2][m];
	  *(pSnd++) = pRG->l2imu[ifr][k][i][3][m];
	  if (noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][i][6][m];
	    *(pSnd++) = pRG->l2imu[ifr][k][i][7][m];
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = GhstZnsIntl[ifr][k][js][l][m];
	  *(pSnd++) = GhstZnsIntr[ifr][k][js][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][k][js][0][m];
	*(pSnd++) = pRG->l1imu[ifr][k][js][2][m];
	*(pSnd++) = pRG->r1imu[ifr][k][js][1][m];
	*(pSnd++) = pRG->r1imu[ifr][k][js][3][m];
	if (noct == 8) {
	  *(pSnd++) = pRG->l1imu[ifr][k][js][4][m];
	  *(pSnd++) = pRG->l1imu[ifr][k][js][6][m];
	  *(pSnd++) = pRG->r1imu[ifr][k][js][5][m];
	  *(pSnd++) = pRG->r1imu[ifr][k][js][7][m];
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=2,3,6,7 are passed using r2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][js][i][4][m];
	  *(pSnd++) = pRG->l3imu[ifr][js][i][5][m];
	  *(pSnd++) = pRG->r3imu[ifr][js][i][0][m];
	  *(pSnd++) = pRG->r3imu[ifr][js][i][1][m];
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int je = pRG->je;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[k][je][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[k][je][i][ifr].H[l];
	  }
	}}}
  }
  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->r2imu[ifr][k][i][0][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][i][1][m];
	  if (noct == 8) {
	    *(pSnd++) = pRG->r2imu[ifr][k][i][4][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][i][5][m];	    
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  *(pSnd++) = GhstZnsIntl[ifr][k][je][l][m];
	  *(pSnd++) = GhstZnsIntr[ifr][k][je][l][m];
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][k][je][0][m];
	*(pSnd++) = pRG->l1imu[ifr][k][je][2][m];
	*(pSnd++) = pRG->r1imu[ifr][k][je][1][m];
	*(pSnd++) = pRG->r1imu[ifr][k][je][3][m];
	if (noct == 8) {
	  *(pSnd++) = pRG->l1imu[ifr][k][je][4][m];
	  *(pSnd++) = pRG->l1imu[ifr][k][je][6][m];
	  *(pSnd++) = pRG->r1imu[ifr][k][je][5][m];
	  *(pSnd++) = pRG->r1imu[ifr][k][je][7][m];
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,1,4,5 are updated using l2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][je][i][6][m];
	  *(pSnd++) = pRG->l3imu[ifr][je][i][7][m];
	  *(pSnd++) = pRG->r3imu[ifr][je][i][2][m];
	  *(pSnd++) = pRG->r3imu[ifr][je][i][3][m];
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[ks][j][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[ks][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l3imu[ifr][j][i][l][m];
	  }}}}}
/* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][ks][j][0][m];
	*(pSnd++) = pRG->l1imu[ifr][ks][j][2][m];
	*(pSnd++) = pRG->l1imu[ifr][ks][j][4][m];
	*(pSnd++) = pRG->l1imu[ifr][ks][j][6][m];
	*(pSnd++) = pRG->r1imu[ifr][ks][j][1][m];
	*(pSnd++) = pRG->r1imu[ifr][ks][j][3][m];
	*(pSnd++) = pRG->r1imu[ifr][ks][j][5][m];
	*(pSnd++) = pRG->r1imu[ifr][ks][j][7][m];
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l2imu[ifr][ks][i][0][m];
	*(pSnd++) = pRG->l2imu[ifr][ks][i][1][m];
	*(pSnd++) = pRG->l2imu[ifr][ks][i][4][m];
	*(pSnd++) = pRG->l2imu[ifr][ks][i][5][m];
	*(pSnd++) = pRG->r2imu[ifr][ks][i][2][m];
	*(pSnd++) = pRG->r2imu[ifr][ks][i][3][m];
	*(pSnd++) = pRG->r2imu[ifr][ks][i][6][m];
	*(pSnd++) = pRG->r2imu[ifr][ks][i][7][m];
      }}}

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  *(pSnd++) = pRG->R[ke][j][i][ifr].S;
    	  for(l=0; l < nDim; l++) {
	    *(pSnd++) = pRG->R[ke][j][i][ifr].H[l];
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r3imu[ifr][j][i][l][m];
	  }}}}}
 /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l1imu[ifr][ke][j][0][m];
	*(pSnd++) = pRG->l1imu[ifr][ke][j][2][m];
	*(pSnd++) = pRG->l1imu[ifr][ke][j][4][m];
	*(pSnd++) = pRG->l1imu[ifr][ke][j][6][m];
	*(pSnd++) = pRG->r1imu[ifr][ke][j][1][m];
	*(pSnd++) = pRG->r1imu[ifr][ke][j][3][m];
	*(pSnd++) = pRG->r1imu[ifr][ke][j][5][m];
	*(pSnd++) = pRG->r1imu[ifr][ke][j][7][m];
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->l2imu[ifr][ke][i][0][m];
	*(pSnd++) = pRG->l2imu[ifr][ke][i][1][m];
	*(pSnd++) = pRG->l2imu[ifr][ke][i][4][m];
	*(pSnd++) = pRG->l2imu[ifr][ke][i][5][m];
	*(pSnd++) = pRG->r2imu[ifr][ke][i][2][m];
	*(pSnd++) = pRG->r2imu[ifr][ke][i][3][m];
	*(pSnd++) = pRG->r2imu[ifr][ke][i][6][m];
	*(pSnd++) = pRG->r2imu[ifr][ke][i][7][m];
      }}}

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);
  /* if LTE receive the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][il][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][il][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update l1imu using r1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l1imu[ifr][k][j][0][m] = *(pRcv++);
	  if(noct > 2) {
	    pRG->l1imu[ifr][k][j][2][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->l1imu[ifr][k][j][4][m] = *(pRcv++);
	      pRG->l1imu[ifr][k][j][6][m] = *(pRcv++);
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    GhstZnsIntl[ifr][k][j][l][m]= *(pRcv++);
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=0,2,4,6 are passed using r1imu and handled below */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][3][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][il][1][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][il][7][m] = *(pRcv++);
	      pRG->r2imu[ifr][k][il][5][m] = *(pRcv++);
	    }
	    }}}
  }
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,2,4,6 are passed using r1imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][il][5][m] = *(pRcv++);
	  pRG->l3imu[ifr][j][il][7][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][il][1][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][il][3][m] = *(pRcv++);
	}}}
  }
  /* copy duplicate edges/coners from corresponding l1imu */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][2][m] = pRG->l1imu[ifr][k][jl][2][m];
	  pRG->r2imu[ifr][k][il][0][m] = pRG->l1imu[ifr][k][ju][0][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][il][6][m] = pRG->l1imu[ifr][k][jl][6][m];
	      pRG->r2imu[ifr][k][il][4][m] = pRG->l1imu[ifr][k][ju][4][m];
	    }
	}}}
  }
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][il][4][m] = pRG->l1imu[ifr][kl][j][4][m];
	  pRG->l3imu[ifr][j][il][6][m] = pRG->l1imu[ifr][kl][j][6][m];
	  pRG->r3imu[ifr][j][il][0][m] = pRG->l1imu[ifr][ku][j][0][m];
	  pRG->r3imu[ifr][j][il][2][m] = pRG->l1imu[ifr][ku][j][2][m];
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1_rad(RadGridS *pRG, int sflag)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);
  /* if LTE receive the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][j][iu][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][j][iu][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update r1imu using l1imu */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(m=0; m<nang; m++) {
	  pRG->r1imu[ifr][k][j][1][m] = *(pRcv++);
	  if(noct > 2) {
	    pRG->r1imu[ifr][k][j][3][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->r1imu[ifr][k][j][5][m] = *(pRcv++);
	      pRG->r1imu[ifr][k][j][7][m] = *(pRcv++);
	    }
	  }
	}}}}
#ifdef SHEARING_BOX
  /* If shearing periodic pass all intensities, which are needed to set l2imu/r2imu
     in corners/edges after remapping */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    GhstZnsIntr[ifr][k][j][l][m]= *(pRcv++);
	  }}}}}
#endif
  /* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
   * l=1,3,5,7 are updated using r1imu below */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][iu][0][m] = *(pRcv++);
	  if(noct == 8) {
	    pRG->l2imu[ifr][k][iu][6][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][iu][4][m] = *(pRcv++);
	  }
	}}}
  }
 /* pass l3imu/r3imu on the edges.  Note that values
   * l=1,3,5,7 are update using r1imu below*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][iu][4][m] = *(pRcv++);
	  pRG->l3imu[ifr][j][iu][6][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][iu][0][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][iu][2][m] = *(pRcv++);
	}}}
  }
  /* copy duplicate edges/coners from corresponding r1imu */
  if (noct > 2) {
    for (ifr=0; ifr<nf; ifr++) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][3][m] = pRG->r1imu[ifr][k][jl][3][m];
	  pRG->r2imu[ifr][k][iu][1][m] = pRG->r1imu[ifr][k][ju][1][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][iu][7][m] = pRG->r1imu[ifr][k][jl][7][m];
	      pRG->r2imu[ifr][k][iu][5][m] = pRG->r1imu[ifr][k][ju][5][m];
	    }
	}}}
  }
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][iu][5][m] = pRG->r1imu[ifr][kl][j][5][m];
	  pRG->l3imu[ifr][j][iu][7][m] = pRG->r1imu[ifr][kl][j][7][m];
	  pRG->r3imu[ifr][j][iu][1][m] = pRG->r1imu[ifr][ku][j][1][m];
	  pRG->r3imu[ifr][j][iu][3][m] = pRG->r1imu[ifr][ku][j][3][m];
	}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][jl][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][jl][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }

  /* update l2imu using r2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][i][0][m] = *(pRcv++);
	  pRG->l2imu[ifr][k][i][1][m] = *(pRcv++);
	  if (noct == 8) {
	    pRG->l2imu[ifr][k][i][4][m] = *(pRcv++);
	    pRG->l2imu[ifr][k][i][5][m] = *(pRcv++);
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][jl][l][m] = *(pRcv++);
	  GhstZnsIntr[ifr][k][jl][l][m] = *(pRcv++);
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][jl][0][m] = *(pRcv++);
	pRG->l1imu[ifr][k][jl][2][m] = *(pRcv++);
	pRG->r1imu[ifr][k][jl][1][m] = *(pRcv++);
	pRG->r1imu[ifr][k][jl][3][m] = *(pRcv++);
	if (noct == 8) {
	  pRG->l1imu[ifr][k][jl][4][m] = *(pRcv++);
	  pRG->l1imu[ifr][k][jl][6][m] = *(pRcv++);
	  pRG->r1imu[ifr][k][jl][5][m] = *(pRcv++);
	  pRG->r1imu[ifr][k][jl][7][m] = *(pRcv++);
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=0,1,4,5 are updated using l2imu below*/
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][jl][i][6][m] = *(pRcv++);
	  pRG->l3imu[ifr][jl][i][7][m] = *(pRcv++);
	  pRG->r3imu[ifr][jl][i][2][m] = *(pRcv++);
	  pRG->r3imu[ifr][jl][i][3][m] = *(pRcv++);
	}}}
  }
  /* copy duplicate edges from corresponding l2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][jl][i][4][m] = pRG->l2imu[ifr][kl][i][4][m];
	  pRG->l3imu[ifr][jl][i][5][m] = pRG->l2imu[ifr][kl][i][5][m];
	  pRG->r3imu[ifr][jl][i][0][m] = pRG->l2imu[ifr][ku][i][0][m];
	  pRG->r3imu[ifr][jl][i][1][m] = pRG->l2imu[ifr][ku][i][1][m];
	}}}
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[k][ju][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[k][ju][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update r2imu using l2imu. Note that i runs from
   * is-1 to ie+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->r2imu[ifr][k][i][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][i][3][m] = *(pRcv++);
	  if (noct == 8) {
	    pRG->r2imu[ifr][k][i][6][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][i][7][m] = *(pRcv++);
	  }
	}}}}
#ifdef SHEARING_BOX
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  GhstZnsIntl[ifr][k][ju][l][m] = *(pRcv++);
	  GhstZnsIntr[ifr][k][ju][l][m] = *(pRcv++);
	}}}} 
#endif
  /* update r1imu/l1imu on corners[2d]/edges[3d] */
  for (ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][k][ju][0][m] = *(pRcv++);
	pRG->l1imu[ifr][k][ju][2][m] = *(pRcv++);
	pRG->r1imu[ifr][k][ju][1][m] = *(pRcv++);
	pRG->r1imu[ifr][k][ju][3][m] = *(pRcv++);
	if (noct == 8) {
	  pRG->l1imu[ifr][k][ju][4][m] = *(pRcv++);
	  pRG->l1imu[ifr][k][ju][6][m] = *(pRcv++);
	  pRG->r1imu[ifr][k][ju][5][m] = *(pRcv++);
	  pRG->r1imu[ifr][k][ju][7][m] = *(pRcv++);
	}
      }}}
  /* pass l3imu/r3imu on the edges.  Note that values
   * l=2,3,6,7 are updated using r2imu below */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][ju][i][4][m] = *(pRcv++);
	  pRG->l3imu[ifr][ju][i][5][m] = *(pRcv++);
	  pRG->r3imu[ifr][ju][i][0][m] = *(pRcv++);
	  pRG->r3imu[ifr][ju][i][1][m] = *(pRcv++);
	}}}
  }
  /* copy duplicate edges from corresponding l2imu */
  if (noct == 8) {
    for (ifr=0; ifr<nf; ifr++) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][ju][i][6][m] = pRG->r2imu[ifr][kl][i][6][m];
	  pRG->l3imu[ifr][ju][i][7][m] = pRG->r2imu[ifr][kl][i][7][m];
	  pRG->r3imu[ifr][ju][i][2][m] = pRG->r2imu[ifr][ku][i][2][m];
	  pRG->r3imu[ifr][ju][i][3][m] = pRG->r2imu[ifr][ku][i][3][m];
	}}}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[kl][j][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[kl][j][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
  /* update l3imu using r3imu. Note that i runs from  is-1 to ie+1 and 
   * j runs from js-1 to je+1 so loop includes ghost zones.*/
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][i][l][m] = *(pRcv++);
	  }}}}}
 /* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][kl][j][0][m] = *(pRcv++);
	pRG->l1imu[ifr][kl][j][2][m] = *(pRcv++);
	pRG->l1imu[ifr][kl][j][4][m] = *(pRcv++);
	pRG->l1imu[ifr][kl][j][6][m] = *(pRcv++);
	pRG->r1imu[ifr][kl][j][1][m] = *(pRcv++);
	pRG->r1imu[ifr][kl][j][3][m] = *(pRcv++);
	pRG->r1imu[ifr][kl][j][5][m] = *(pRcv++);
	pRG->r1imu[ifr][kl][j][7][m] = *(pRcv++);
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->l2imu[ifr][kl][i][0][m] = *(pRcv++);
	pRG->l2imu[ifr][kl][i][1][m] = *(pRcv++);
	pRG->l2imu[ifr][kl][i][4][m] = *(pRcv++);
	pRG->l2imu[ifr][kl][i][5][m] = *(pRcv++);
	pRG->r2imu[ifr][kl][i][2][m] = *(pRcv++);
	pRG->r2imu[ifr][kl][i][3][m] = *(pRcv++);
	pRG->r2imu[ifr][kl][i][6][m] = *(pRcv++);
	pRG->r2imu[ifr][kl][i][7][m] = *(pRcv++);
      }}}
 

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3_rad(RadGridS *pRG, int sflag)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  /* if LTE pass the source functions */
  if (sflag == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	for (ifr=0; ifr<nf; ifr++) {
	  pRG->R[ku][j][i][ifr].S = *(pRcv++);
    	  for(l=0; l < nDim; l++) {
	    pRG->R[ku][j][i][ifr].H[l] = *(pRcv++);
	  }
	}}}
  }
 /* update r3imu using l3imu. Note that i runs from  is-1 to ie+1 and 
  * j runs from js-1 to je+1 so loop includes ghost zones. */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->r3imu[ifr][j][i][l][m] = *(pRcv++);
	  }}}}}
/* update r1imu/l1imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->l1imu[ifr][ku][j][0][m] = *(pRcv++);
	pRG->l1imu[ifr][ku][j][2][m] = *(pRcv++);
	pRG->l1imu[ifr][ku][j][4][m] = *(pRcv++);
	pRG->l1imu[ifr][ku][j][6][m] = *(pRcv++);
	pRG->r1imu[ifr][ku][j][1][m] = *(pRcv++);
	pRG->r1imu[ifr][ku][j][3][m] = *(pRcv++);
	pRG->r1imu[ifr][ku][j][5][m] = *(pRcv++);
	pRG->r1imu[ifr][ku][j][7][m] = *(pRcv++);
      }}}
 /* update r2imu/l2imu on edges */
  for (ifr=0; ifr<nf; ifr++) {
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->l2imu[ifr][ku][i][0][m] = *(pRcv++);
	pRG->l2imu[ifr][ku][i][1][m] = *(pRcv++);
	pRG->l2imu[ifr][ku][i][4][m] = *(pRcv++);
	pRG->l2imu[ifr][ku][i][5][m] = *(pRcv++);
	pRG->r2imu[ifr][ku][i][2][m] = *(pRcv++);
	pRG->r2imu[ifr][ku][i][3][m] = *(pRcv++);
	pRG->r2imu[ifr][ku][i][6][m] = *(pRcv++);
	pRG->r2imu[ifr][ku][i][7][m] = *(pRcv++);
      }}}

  return;
}
#endif /* JACOBI */
#endif /* MPI_PARALLEL */

#endif /* RADIATION_TRANSFER */
