#include "../copyright.h"
/*=============================================================================
 * FILE: bvals_rad.c
 *
 * PURPOSE: Sets boundary condistions for radiative transfer on each edge of
 *          the grid.  It closely follows the methods and conventions and
 *          methods used for the hydro integration (see e.g. bvals_mhd.c).
 *          The radiation source function, radiative flux, and intensities
 *          on the boundaries of the RadGrid are copied to ghost zones.
 *          Hence, for the purposes of the formal_solution functions, all
 *          zones (including boundary zones) are treated identically.
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
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   periodic_??_rad?() - periodic BCs at boundary ???
 *   pack_???_rad()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???_rad()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/
static void periodic_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void periodic_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void periodic_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void periodic_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void periodic_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void periodic_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);

static void ProlongateLater(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_incident_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);

static void const_flux_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_flux_ox1(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_flux_ix2(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_flux_ox2(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_flux_ix3(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void const_flux_ox3(GridS *pG, RadGridS *pRG, int ifs, int ife);

static void constH_varJ_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife);

#ifdef MPI_PARALLEL
static void pack_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void pack_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void pack_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void pack_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void pack_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void pack_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);

static void unpack_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void unpack_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void unpack_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void unpack_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void unpack_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
static void unpack_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife);
#endif /* MPI_PARALLEL */

void bvals_rad(DomainS *pD, int ifs, int ife)
{
  RadGridS *pRG=(pD->RadGrid);
  GridS *pG=(pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN,BCFlag;
#endif
#ifdef MPI_PARALLEL
  int cnt, ierr, mIndex;
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
      pack_ix1_rad(pG,pRG,ifs,ife);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox1_rad(pG,pRG,ifs,ife);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox1_rad(pG,pRG,ifs,ife);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));
      /* set physical boundary */
      (*(pD->ix1_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_rad(pG,pRG,ifs,ife);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx1_id < 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_rad(pG,pRG,ifs,ife);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx1_id < 0 && pRG->lx1_id < 0) {
      (*(pD->ix1_RBCFun))(pG,pRG,ifs,ife);
      (*(pD->ox1_RBCFun))(pG,pRG,ifs,ife);
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
      pack_ix2_rad(pG,pRG,ifs,ife);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox2_rad(pG,pRG,ifs,ife);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox2_rad(pG,pRG,ifs,ife);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id < 0) {
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_rad(pG,pRG,ifs,ife);


    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx2_id < 0 && pRG->lx2_id >= 0) {
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_rad(pG,pRG,ifs,ife);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx2_id < 0 && pRG->lx2_id < 0) {
      (*(pD->ix2_RBCFun))(pG,pRG,ifs,ife);
      (*(pD->ox2_RBCFun))(pG,pRG,ifs,ife);
    } 
 
/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */


#ifdef SHEARING_BOX 
    BCFlag = par_geti_def("domain1","rbc_ix1",0);
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0 && BCFlag == 4) {
      ShearingSheet_Rad_ix1(pD,ifs,ife);
      } 
    BCFlag = par_geti_def("domain1","rbc_ox1",0);
    if (myL == ((pD->NGrid[0])-1) && BCFlag == 4) {
      ShearingSheet_Rad_ox1(pD,ifs,ife);
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
      pack_ix3_rad(pG,pRG,ifs,ife);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox3_rad(pG,pRG,ifs,ife);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_rad(pG,pRG,ifs,ife);
      if (mIndex == 1) unpack_ox3_rad(pG,pRG,ifs,ife);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_rad(pG,pRG,ifs,ife);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx3_id < 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3_rad(pG,pRG,ifs,ife); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_RBCFun))(pG,pRG,ifs,ife);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_rad(pG,pRG,ifs,ife);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pRG->rx3_id < 0 && pRG->lx3_id < 0) {
      (*(pD->ix3_RBCFun))(pG,pRG,ifs,ife);
      (*(pD->ox3_RBCFun))(pG,pRG,ifs,ife);
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
      
      if(pD->ix1_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if(pD->Disp[0] != 0) {      
	  pD->ix1_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
	} else {                    
	  switch(pM->RBCFlag_ix1){

	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ix1_RBCFun = const_incident_rad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix1_RBCFun = periodic_ix1_rad;
#ifdef MPI_PARALLEL
	    if(pRG->lx1_id < 0 && pD->NGrid[0] > 1){
	      pRG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix1_RBCFun = const_flux_ix1;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ix1=%d unknown\n",pM->RBCFlag_ix1);
	    exit(EXIT_FAILURE);
	  }
	}
      }


/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->ox1_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
	  pD->ox1_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
	} else {
	  switch(pM->RBCFlag_ox1){

	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ox1_RBCFun = const_incident_rad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox1_RBCFun = periodic_ox1_rad;
#ifdef MPI_PARALLEL
	    if(pRG->rx1_id < 0 && pD->NGrid[0] > 1){
	      pRG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;
	    
	  case 5: /* constant fixed flux at boundary */
	    pD->ox1_RBCFun = const_flux_ox1;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ox1=%d unknown\n",pM->RBCFlag_ox1);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }

    if(pRG->Nx[1] > 1) {

      nDim = 2;
/*---- ix2 boundary ----------------------------------------------------------*/

      if(pD->ix2_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if(pD->Disp[1] != 0) {      
	  pD->ix2_RBCFun = ProlongateLater;
/* Domain is at L-edge of root Domain */
	} else {                    
	  switch(pM->RBCFlag_ix2){

	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ix2_RBCFun = const_incident_rad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix2_RBCFun = periodic_ix2_rad;
#ifdef MPI_PARALLEL
	    if(pRG->lx2_id < 0 && pD->NGrid[1] > 1){
	      pRG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix2_RBCFun = const_flux_ix2;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ix2=%d unknown\n",pM->RBCFlag_ix2);
	    exit(EXIT_FAILURE);
	  }
	}
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(pD->ox2_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
	  pD->ox2_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
	} else {
	  switch(pM->RBCFlag_ox2){

	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ox2_RBCFun = const_incident_rad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox2_RBCFun = periodic_ox2_rad;
#ifdef MPI_PARALLEL
	    if(pRG->rx2_id < 0 && pD->NGrid[1] > 1){
	      pRG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ox2_RBCFun = const_flux_ox2;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ox2=%d unknown\n",pM->RBCFlag_ox2);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }

    if(pRG->Nx[2] > 1) {

      nDim = 3;
/*---- ix3 boundary ----------------------------------------------------------*/

      if(pD->ix3_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if(pD->Disp[2] != 0) {      
	  pD->ix3_RBCFun = ProlongateLater;

/* Domain is at L-edge of root Domain */
	} else {                    
	  switch(pM->RBCFlag_ix3){
	    
	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ix3_RBCFun = const_incident_rad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix3_RBCFun = periodic_ix3_rad;
#ifdef MPI_PARALLEL
	    if(pRG->lx3_id < 0 && pD->NGrid[2] > 1){
	      pRG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix3_RBCFun = const_flux_ix3;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ix3=%d unknown\n",pM->RBCFlag_ix3);
	    exit(EXIT_FAILURE);
	  }
	}
      }

/*---- ox3 boundary ----------------------------------------------------------*/

      if(pD->ox3_RBCFun == NULL) {    /* RBCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
	if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
	  pD->ox3_RBCFun = ProlongateLater;
/* Domain is at R-edge of root Domain */
	} else {
	  switch(pM->RBCFlag_ox3){

	  case 2: /* Open boundary with fixed incoming radiation */
	    pD->ox3_RBCFun = const_incident_rad;
	    break;
	    
	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox3_RBCFun = periodic_ox3_rad;
#ifdef MPI_PARALLEL
	    if(pRG->rx3_id < 0 && pD->NGrid[2] > 1){
	      pRG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ox3_RBCFun = const_flux_ox3;
	    break;

	  default:
	    ath_perr(-1,"[bvals_rad_init]:rbc_ox3=%d unknown\n",pM->RBCFlag_ox3);
	    exit(EXIT_FAILURE);
	  }
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

/* x1cnt is the number of Reals passed for x1 faces */
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  nx3t = pD->GData[n][m][l].Nx[2];
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	  xcnt = nx2t * nx3t * noct * nang;
#else
	  xcnt = nx2t * nx3t * noct * nang / 2;
#endif
#ifdef QUADRATIC_INTENSITY
	  if (noct > 2)  xcnt += nx3t * noct * nang * 2;
	  if (noct == 8) xcnt += nx2t * noct * nang * 2;
#else
	  if (noct > 2)  xcnt += nx3t * noct * nang / 2;
	  if (noct == 8) xcnt += nx2t * noct * nang / 2;
#endif
	  xcnt += nx2t * nx3t * (nDim + 2);	  
	  xcnt *= nf;
          if(xcnt > x1cnt) x1cnt = xcnt;
	}

/* x2cnt is the number of Reals passed for x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2;
	  nx3t = pD->GData[n][m][l].Nx[2];
#ifdef QUADRATIC_INTENSITY	  
	  xcnt = nx1t * nx3t * noct * nang;
#else
	  xcnt = nx1t * nx3t * noct * nang / 2;
#endif
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	  xcnt += nx3t * noct * nang * 2;
#else
	  xcnt += nx3t * noct * nang;
#endif
#ifdef QUADRATIC_INTENSITY
	  if (noct == 8) xcnt += nx1t * noct * nang * 2; 
#else
	  if (noct == 8) xcnt += nx1t * noct * nang / 2;   
#endif
	  xcnt += nx3t * nx1t * (nDim + 2);
	  xcnt *= nf;
          if(xcnt > x2cnt) x2cnt = xcnt;
	}

/* x3cnt is the number of Reals passed for x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2;
	  nx2t = pD->GData[n][m][l].Nx[1] + 2;
#ifdef QUADRATIC_INTENSITY
	  xcnt = nx1t * nx2t * noct * nang;
#else
	  xcnt = nx1t * nx2t * noct * nang / 2;
#endif
#ifdef QUADRATIC_INTENSITY
	  xcnt += (nx1t + nx2t) * noct * nang * 2; 
#else
	  xcnt += (nx1t + nx2t) * noct * nang; 
#endif
	  xcnt += nx2t * nx1t * (nDim + 2);
	  xcnt *= nf;
          if(xcnt > x3cnt) x3cnt = xcnt;
	}
	/*	printf("counts: %d %d %d %d\n",x1cnt,x2cnt,x3cnt,xcnt);*/
      }
    }}
#endif /* MPI_PARALLEL */

  }}}  /* End loop over all Domains with active Grids -----------------------*/

#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  if (size > 0) {
    if((send_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_rad_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_rad_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_rad_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_rad_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */


  return;
}

void bvals_rad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc)
{

  switch(dir){
  case left_x1:
    pD->ix1_RBCFun = prob_bc;
    break;
  case right_x1:
    pD->ox1_RBCFun = prob_bc;
    break;
  case left_x2:
    pD->ix2_RBCFun = prob_bc;
    break;
  case right_x2:
    pD->ox2_RBCFun = prob_bc;
    break;
  case left_x3:
    pD->ix3_RBCFun = prob_bc;
    break;
  case right_x3:
    pD->ox3_RBCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_rad_trans_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   periodic_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   const_flux_???
 *   pack_???
 *   unpack_???
 *  const_incident_rad
 */


/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (rbc_ix1=1) */

static void periodic_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, ie = pRG->ie;
  int jl = pRG->js  , ju = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {

/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	pRG->R[ifr][k][j][il].S = pRG->R[ifr][k][j][ie].S;
	pRG->R[ifr][k][j][il].J = pRG->R[ifr][k][j][ie].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][j][il].H[l] = pRG->R[ifr][k][j][ie].H[l];
	}
      }}
/* update Ghstl1i using r1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][j][l][m] = pRG->r1imu[ifr][k][j][l][m];
	  }}
#else
	for (m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = pRG->r1imu[ifr][k][j][0][m];
	  if(noct > 2) {
	    pRG->Ghstl1i[ifr][k][j][2][m] = pRG->r1imu[ifr][k][j][2][m];
	    if(noct == 8) {
	      pRG->Ghstl1i[ifr][k][j][4][m] = pRG->r1imu[ifr][k][j][4][m];
	      pRG->Ghstl1i[ifr][k][j][6][m] = pRG->r1imu[ifr][k][j][6][m];
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][il][l][m] = pRG->Ghstl1i[ifr][k][jl][l][m];
	    pRG->r2imu[ifr][k][il][l][m] = pRG->Ghstl1i[ifr][k][ju][l][m];
	  }}}
    }
/* pass l3imu/r3imu on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][il][l][m] = pRG->Ghstl1i[ifr][kl][j][l][m];
	    pRG->r3imu[ifr][j][il][l][m] = pRG->Ghstl1i[ifr][ku][j][l][m];
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=0,2,4,6 are passed using r1imu */
    if (noct > 2) {
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
	}}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=0,2,4,6 are passed using r1imu */
    if (noct == 8) {
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
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (rbc_ox1=1) */

static void periodic_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int is = pRG->is, iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	pRG->R[ifr][k][j][iu].S = pRG->R[ifr][k][j][is].S;
	pRG->R[ifr][k][j][iu].J = pRG->R[ifr][k][j][is].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][j][iu].H[l] = pRG->R[ifr][k][j][is].H[l];
	}
      }}
/* update Ghstr1i using l1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][l][m] = pRG->l1imu[ifr][k][j][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstr1i[ifr][k][j][1][m] = pRG->l1imu[ifr][k][j][1][m];
	  if(noct > 2) {
	    pRG->Ghstr1i[ifr][k][j][3][m] = pRG->l1imu[ifr][k][j][3][m];
	    if(noct == 8) {
	      pRG->Ghstr1i[ifr][k][j][5][m] = pRG->l1imu[ifr][k][j][5][m];
	      pRG->Ghstr1i[ifr][k][j][7][m] = pRG->l1imu[ifr][k][j][7][m];
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][iu][l][m] = pRG->Ghstr1i[ifr][k][jl][l][m];
	    pRG->r2imu[ifr][k][iu][l][m] = pRG->Ghstr1i[ifr][k][ju][l][m];
	  }}}
    }
/* pass l3imu/r3imu on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][iu][l][m] = pRG->Ghstr1i[ifr][kl][j][l][m];
	    pRG->r3imu[ifr][j][iu][l][m] = pRG->Ghstr1i[ifr][ku][j][l][m];
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=1,3,5,7 are passed using l1imu */
    if (noct > 2) {
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
	}}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=1,3,5,7 are passed using r1imu */
    if (noct == 8) {
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
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=1) */

static void periodic_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, je = pRG->je;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][k][jl][i].S = pRG->R[ifr][k][je][i].S;
	pRG->R[ifr][k][jl][i].J = pRG->R[ifr][k][je][i].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][jl][i].H[l] = pRG->R[ifr][k][je][i].H[l];
	}
      }}
/* update Ghstl2i using r2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones.  r2imu on the
 * corners/edges has already been updated (if necessary) by
 * x1 boundary routines */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl2i[ifr][k][i][l][m] = pRG->r2imu[ifr][k][i][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = pRG->r2imu[ifr][k][i][0][m];
	  pRG->Ghstl2i[ifr][k][i][1][m] = pRG->r2imu[ifr][k][i][1][m];
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = pRG->r2imu[ifr][k][i][4][m];
	    pRG->Ghstl2i[ifr][k][i][5][m] = pRG->r2imu[ifr][k][i][5][m];	    
	  }
	}
#endif
      }}
/* update Ghstl1i/Ghstr1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][jl][l][m] = pRG->Ghstl1i[ifr][k][je][l][m];
	  pRG->Ghstr1i[ifr][k][jl][l][m] = pRG->Ghstr1i[ifr][k][je][l][m];
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][k][jl][0][m] = pRG->Ghstl1i[ifr][k][je][0][m];
	pRG->Ghstl1i[ifr][k][jl][2][m] = pRG->Ghstl1i[ifr][k][je][2][m];
	pRG->Ghstr1i[ifr][k][jl][1][m] = pRG->Ghstr1i[ifr][k][je][1][m];
	pRG->Ghstr1i[ifr][k][jl][3][m] = pRG->Ghstr1i[ifr][k][je][3][m];
	if (noct == 8) {
	  pRG->Ghstl1i[ifr][k][jl][4][m] = pRG->Ghstl1i[ifr][k][je][4][m];
	  pRG->Ghstl1i[ifr][k][jl][6][m] = pRG->Ghstl1i[ifr][k][je][6][m];
	  pRG->Ghstr1i[ifr][k][jl][5][m] = pRG->Ghstr1i[ifr][k][je][5][m];
	  pRG->Ghstr1i[ifr][k][jl][7][m] = pRG->Ghstr1i[ifr][k][je][7][m];
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][jl][i][l][m] = pRG->Ghstl2i[ifr][kl][i][l][m];
	    pRG->r3imu[ifr][jl][i][l][m] = pRG->Ghstl2i[ifr][ku][i][l][m];
	  }}}
    }
#else
/* Note that values l=0,1,4,5 are passed using r2imu */
    if (noct == 8) {
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
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=1) */

static void periodic_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js  , ju = pRG->je+1;
  int kl = pRG->ks  , ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][k][ju][i].S = pRG->R[ifr][k][js][i].S;
	pRG->R[ifr][k][ju][i].J = pRG->R[ifr][k][js][i].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][ju][i].H[l] = pRG->R[ifr][k][js][i].H[l];
	}
      }}
/* update Ghstr2i using l2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones.  l2imu on the
 * corners/edges has already been updated (if necessary) by
 * x1 boundary routines */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr2i[ifr][k][i][l][m] = pRG->l2imu[ifr][k][i][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][k][i][2][m] = pRG->l2imu[ifr][k][i][2][m];
	  pRG->Ghstr2i[ifr][k][i][3][m] = pRG->l2imu[ifr][k][i][3][m];
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = pRG->l2imu[ifr][k][i][6][m];
	    pRG->Ghstr2i[ifr][k][i][7][m] = pRG->l2imu[ifr][k][i][7][m];	    
	  }
	}
#endif
      }}
/* update Ghstr1i/Ghstl1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][ju][l][m] = pRG->Ghstl1i[ifr][k][js][l][m];
	  pRG->Ghstr1i[ifr][k][ju][l][m] = pRG->Ghstr1i[ifr][k][js][l][m];
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][k][ju][0][m] = pRG->Ghstl1i[ifr][k][js][0][m];
	pRG->Ghstl1i[ifr][k][ju][2][m] = pRG->Ghstl1i[ifr][k][js][2][m];
	pRG->Ghstr1i[ifr][k][ju][1][m] = pRG->Ghstr1i[ifr][k][js][1][m];
	pRG->Ghstr1i[ifr][k][ju][3][m] = pRG->Ghstr1i[ifr][k][js][3][m];
	if (noct == 8) {
	  pRG->Ghstl1i[ifr][k][ju][4][m] = pRG->Ghstl1i[ifr][k][js][4][m];
	  pRG->Ghstl1i[ifr][k][ju][6][m] = pRG->Ghstl1i[ifr][k][js][6][m];
	  pRG->Ghstr1i[ifr][k][ju][5][m] = pRG->Ghstr1i[ifr][k][js][5][m];
	  pRG->Ghstr1i[ifr][k][ju][7][m] = pRG->Ghstr1i[ifr][k][js][7][m];
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][ju][i][l][m] = pRG->Ghstr2i[ifr][kl][i][l][m];
	    pRG->r3imu[ifr][ju][i][l][m] = pRG->Ghstr2i[ifr][ku][i][l][m];
	  }}}
    }
#else
/* Note that values l=2,3,6,7 are passed using l2imu */
    if (noct == 8) {
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
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (rbc_ix3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][kl][j][i].S = pRG->R[ifr][ke][j][i].S;
	pRG->R[ifr][kl][j][i].J = pRG->R[ifr][ke][j][i].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][kl][j][i].H[l] = pRG->R[ifr][ke][j][i].H[l];
	}
      }}
/* update Ghstl3i using r3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones.  r3imu on the
 * edges has already been updated (if necessary) by x1 and x2 boundary
 * routines */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = pRG->r3imu[ifr][j][i][l][m];
	  }}
#else
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = pRG->r3imu[ifr][j][i][l][m];
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][kl][j][l][m] = pRG->Ghstl1i[ifr][ke][j][l][m];
	  pRG->Ghstr1i[ifr][kl][j][l][m] = pRG->Ghstr1i[ifr][ke][j][l][m];
	}}}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][kl][i][l][m] = pRG->Ghstl2i[ifr][ke][i][l][m];
	  pRG->Ghstr2i[ifr][kl][i][l][m] = pRG->Ghstr2i[ifr][ke][i][l][m];
	}}}
#else
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][kl][j][0][m] = pRG->Ghstl1i[ifr][ke][j][0][m];
	pRG->Ghstl1i[ifr][kl][j][2][m] = pRG->Ghstl1i[ifr][ke][j][2][m];
	pRG->Ghstl1i[ifr][kl][j][4][m] = pRG->Ghstl1i[ifr][ke][j][4][m];
	pRG->Ghstl1i[ifr][kl][j][6][m] = pRG->Ghstl1i[ifr][ke][j][6][m];
	pRG->Ghstr1i[ifr][kl][j][1][m] = pRG->Ghstr1i[ifr][ke][j][1][m];
	pRG->Ghstr1i[ifr][kl][j][3][m] = pRG->Ghstr1i[ifr][ke][j][3][m];
	pRG->Ghstr1i[ifr][kl][j][5][m] = pRG->Ghstr1i[ifr][ke][j][5][m];
	pRG->Ghstr1i[ifr][kl][j][7][m] = pRG->Ghstr1i[ifr][ke][j][7][m];
      }}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl2i[ifr][kl][i][0][m] = pRG->Ghstl2i[ifr][ke][i][0][m];
	pRG->Ghstl2i[ifr][kl][i][1][m] = pRG->Ghstl2i[ifr][ke][i][1][m];
	pRG->Ghstl2i[ifr][kl][i][4][m] = pRG->Ghstl2i[ifr][ke][i][4][m];
	pRG->Ghstl2i[ifr][kl][i][5][m] = pRG->Ghstl2i[ifr][ke][i][5][m];
	pRG->Ghstr2i[ifr][kl][i][2][m] = pRG->Ghstr2i[ifr][ke][i][2][m];
	pRG->Ghstr2i[ifr][kl][i][3][m] = pRG->Ghstr2i[ifr][ke][i][3][m];
	pRG->Ghstr2i[ifr][kl][i][6][m] = pRG->Ghstr2i[ifr][ke][i][6][m];
	pRG->Ghstr2i[ifr][kl][i][7][m] = pRG->Ghstr2i[ifr][ke][i][7][m];
      }}
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (rbc_ox3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks, ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, l, m, ifr;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	pRG->R[ifr][ku][j][i].S = pRG->R[ifr][ks][j][i].S;
	pRG->R[ifr][ku][j][i].J = pRG->R[ifr][ks][j][i].J;
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][ku][j][i].H[l] = pRG->R[ifr][ks][j][i].H[l];
	}
      }}
/* update Ghstl3i using l3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones.  l3imu on the
 * edges has already been updated (if necessary) by x1 and x2 boundary
 * routines */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr3i[ifr][j][i][l][m] = pRG->l3imu[ifr][j][i][l][m];
	  }}
#else
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr3i[ifr][j][i][l][m] = pRG->l3imu[ifr][j][i][l][m];
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][ku][j][l][m] = pRG->Ghstl1i[ifr][ks][j][l][m];
	  pRG->Ghstr1i[ifr][ku][j][l][m] = pRG->Ghstr1i[ifr][ks][j][l][m];
	}}}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][ku][i][l][m] = pRG->Ghstl2i[ifr][ks][i][l][m];
	  pRG->Ghstr2i[ifr][ku][i][l][m] = pRG->Ghstr2i[ifr][ks][i][l][m];
	}}}
#else
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][ku][j][0][m] = pRG->Ghstl1i[ifr][ks][j][0][m];
	pRG->Ghstl1i[ifr][ku][j][2][m] = pRG->Ghstl1i[ifr][ks][j][2][m];
	pRG->Ghstl1i[ifr][ku][j][4][m] = pRG->Ghstl1i[ifr][ks][j][4][m];
	pRG->Ghstl1i[ifr][ku][j][6][m] = pRG->Ghstl1i[ifr][ks][j][6][m];
	pRG->Ghstr1i[ifr][ku][j][1][m] = pRG->Ghstr1i[ifr][ks][j][1][m];
	pRG->Ghstr1i[ifr][ku][j][3][m] = pRG->Ghstr1i[ifr][ks][j][3][m];
	pRG->Ghstr1i[ifr][ku][j][5][m] = pRG->Ghstr1i[ifr][ks][j][5][m];
	pRG->Ghstr1i[ifr][ku][j][7][m] = pRG->Ghstr1i[ifr][ks][j][7][m];
      }}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl2i[ifr][ku][i][0][m] = pRG->Ghstl2i[ifr][ks][i][0][m];
	pRG->Ghstl2i[ifr][ku][i][1][m] = pRG->Ghstl2i[ifr][ks][i][1][m];
	pRG->Ghstl2i[ifr][ku][i][4][m] = pRG->Ghstl2i[ifr][ks][i][4][m];
	pRG->Ghstl2i[ifr][ku][i][5][m] = pRG->Ghstl2i[ifr][ks][i][5][m];
	pRG->Ghstr2i[ifr][ku][i][2][m] = pRG->Ghstr2i[ifr][ks][i][2][m];
	pRG->Ghstr2i[ifr][ku][i][3][m] = pRG->Ghstr2i[ifr][ks][i][3][m];
	pRG->Ghstr2i[ifr][ku][i][6][m] = pRG->Ghstr2i[ifr][ks][i][6][m];
	pRG->Ghstr2i[ifr][ku][i][7][m] = pRG->Ghstr2i[ifr][ks][i][7][m];
      }}
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Requires a time independent constant flux on the boundary but allows 
 * mean intensity to vary. Inner x1 boundary  WORK IN PROGRESS!!! */

static void constH_varJ_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, is = pRG->is;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  RadS Rs;
  Real prod, dx, eps, J, H;
  
  dx = pRG->dx1;

  for(ifr=ifs; ifr<=ife; ifr++) {
/*  send the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
/* Compute ghost zone J using 2nd moment of RT eq. */
	Rs =pRG->R[ifr][k][j][is];
	H = pRG->R[ifr][k][j][il].H[0];
	J = Rs.J * (1.0 + H * Rs.chi * dx / Rs.K[0]);
	eps = pRG->R[ifr][k][j][il].eps;
	pRG->R[ifr][k][j][il].S = eps * pRG->R[ifr][k][j][il].B +
	  (1.0 - eps) * J;
	pRG->R[ifr][k][j][il].J = J;
/* update Ghstl1i using l1imu */
	prod = J * (J - H) / (Rs.J * (J + H) ); 
	for (m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = prod * pRG->l1imu[ifr][k][j][1][m];
	  if (noct > 2) {
	    pRG->Ghstl1i[ifr][k][j][2][m] = prod * pRG->l1imu[ifr][k][j][3][m];
	    if (noct == 8) {
	      pRG->Ghstl1i[ifr][k][j][4][m] = prod * pRG->l1imu[ifr][k][j][5][m];
	      pRG->Ghstl1i[ifr][k][j][6][m] = prod * pRG->l1imu[ifr][k][j][7][m];
	    }
	  }
	}
      }}
  }

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic.
 * Inner x1 boundary  */

static void const_flux_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  Real I0, Jm, H, Hm, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct > 2)
    if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstl1i using l1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	Hm = 0.0;
	Jm = 0.0;
	for (m=0; m<nang; m++) {
	  Hm += pRG->l1imu[ifr][k][j][0][m] * pRG->mu[0][m][0] * pRG->wmu[m];
	  Jm += pRG->l1imu[ifr][k][j][0][m] * pRG->wmu[m];
	  if (noct > 2) {
	    Hm += pRG->l1imu[ifr][k][j][2][m] * pRG->mu[2][m][0] * pRG->wmu[m];
	    Jm += pRG->l1imu[ifr][k][j][2][m] * pRG->wmu[m];	  
	    if (noct == 8) {
	      Hm += pRG->l1imu[ifr][k][j][4][m] * pRG->mu[4][m][0] * pRG->wmu[m];
	      Hm += pRG->l1imu[ifr][k][j][6][m] * pRG->mu[6][m][0] * pRG->wmu[m];
	      Jm += pRG->l1imu[ifr][k][j][4][m] * pRG->wmu[m];
	      Jm += pRG->l1imu[ifr][k][j][6][m] * pRG->wmu[m];
	    }
	  }
	}
	H = pRG->R[ifr][k][j][il].H[0];
	I0 = (H - Hm) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][k][j][il].J = 0.5 * I0 + Jm;
	for (m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = I0;
	  if (noct > 2) {
	    pRG->Ghstl1i[ifr][k][j][2][m] = I0;
	    if (noct == 8) {
	      pRG->Ghstl1i[ifr][k][j][4][m] = I0;
	      pRG->Ghstl1i[ifr][k][j][6][m] = I0;
	    }
	  }
	}
      }}
  }

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic.
 * Outer x1 boundary  */

static void const_flux_ox1(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  Real I0, Jp, H, Hp, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct > 2)
    if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstr1i using r1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	Hp = 0.0;
	Jp = 0.0;
	for (m=0; m<nang; m++) {
	  Hp += pRG->r1imu[ifr][k][j][1][m] * pRG->mu[1][m][0] * pRG->wmu[m];
	  Jp += pRG->r1imu[ifr][k][j][1][m] * pRG->wmu[m];
	  if (noct > 2) {
	    Hp += pRG->r1imu[ifr][k][j][3][m] * pRG->mu[3][m][0] * pRG->wmu[m];
	    Jp += pRG->r1imu[ifr][k][j][3][m] * pRG->wmu[m];
	    if (noct == 8) {
	      Hp += pRG->r1imu[ifr][k][j][5][m] * pRG->mu[5][m][0] * pRG->wmu[m];
	      Hp += pRG->r1imu[ifr][k][j][7][m] * pRG->mu[7][m][0] * pRG->wmu[m];
	      Jp += pRG->r1imu[ifr][k][j][5][m] * pRG->wmu[m];
	      Jp += pRG->r1imu[ifr][k][j][7][m] * pRG->wmu[m];
	    }
	  }
	}
	H = pRG->R[ifr][k][j][iu].H[0];
	I0 = (Hp - H) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][k][j][iu].J = 0.5 * I0 + Jp;
	for (m=0; m<nang; m++) {
	  pRG->Ghstr1i[ifr][k][j][1][m] = I0;
	  if (noct > 2) {
	    pRG->Ghstr1i[ifr][k][j][3][m] = I0;
	    if (noct == 8) {
	      pRG->Ghstr1i[ifr][k][j][5][m] = I0;
	      pRG->Ghstr1i[ifr][k][j][7][m] = I0;
	    }
	  }
	}
      }}}

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Inner x2 boundary  */

static void const_flux_ix2(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  Real I0, Jm, H, Hm, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstl2i using l2imu */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	Hm = 0.0;
	Jm = 0.0;
	for (m=0; m<nang; m++) {
	  Hm += pRG->l2imu[ifr][k][i][2][m] * pRG->mu[2][m][1] * pRG->wmu[m];
	  Hm += pRG->l2imu[ifr][k][i][3][m] * pRG->mu[3][m][1] * pRG->wmu[m];
	  Jm += pRG->l2imu[ifr][k][i][2][m] * pRG->wmu[m];
	  Jm += pRG->l2imu[ifr][k][i][3][m] * pRG->wmu[m];	  
	  if (noct == 8) {
	    Hm += pRG->l2imu[ifr][k][i][6][m] * pRG->mu[6][m][1] * pRG->wmu[m];
	    Hm += pRG->l2imu[ifr][k][i][7][m] * pRG->mu[7][m][1] * pRG->wmu[m];
	    Jm += pRG->l2imu[ifr][k][i][6][m] * pRG->wmu[m];
	    Jm += pRG->l2imu[ifr][k][i][7][m] * pRG->wmu[m];
	  }
	}	
	H = pRG->R[ifr][k][jl][i].H[1];
	I0 = (H - Hm) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][k][jl][i].J = 0.5 * I0 + Jm;
	for (m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = I0;
	  pRG->Ghstl2i[ifr][k][i][1][m] = I0;
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = I0;
	    pRG->Ghstl2i[ifr][k][i][5][m] = I0;
	  }
	}
      }}
  }

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Outer x2 boundary  */

static void const_flux_ox2(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  Real I0, Jp, H, Hp, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstr2i using r2imu */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	Hp = 0.0;
	Jp = 0.0;
	for (m=0; m<nang; m++) {
	  Hp += pRG->r2imu[ifr][k][i][0][m] * pRG->mu[0][m][1] * pRG->wmu[m];
	  Hp += pRG->r2imu[ifr][k][i][1][m] * pRG->mu[1][m][1] * pRG->wmu[m];
	  Jp += pRG->r2imu[ifr][k][i][0][m] * pRG->wmu[m];
	  Jp += pRG->r2imu[ifr][k][i][1][m] * pRG->wmu[m];
	  if (noct == 8) {
	    Hp += pRG->r2imu[ifr][k][i][4][m] * pRG->mu[4][m][1] * pRG->wmu[m];
	    Hp += pRG->r2imu[ifr][k][i][5][m] * pRG->mu[5][m][1] * pRG->wmu[m];
	    Jp += pRG->r2imu[ifr][k][i][4][m] * pRG->wmu[m];
	    Jp += pRG->r2imu[ifr][k][i][5][m] * pRG->wmu[m];
	  }
	}
	H = pRG->R[ifr][k][ju][i].H[1];
	I0 = (Hp - H) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][k][ju][i].J = 0.5 * I0 + Jp;
	for (m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][k][i][2][m] = I0;
	  pRG->Ghstr2i[ifr][k][i][3][m] = I0;
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = I0;
	    pRG->Ghstr2i[ifr][k][i][7][m] = I0;
	  }
	}
      }}
  }

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Inner x3 boundary  */

static void const_flux_ix3(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1;
  int nang = pRG->nang;
  int i, j, l, m, n, ifr;
  Real I0, Jm, H, Hm, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  gamma *= 4.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstl3i using l3imu */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	Hm = 0.0;
	Jm = 0.0;
	for (m=0; m<nang; m++) {
	  Hm += pRG->l3imu[ifr][j][i][0][m] * pRG->mu[0][m][2] * pRG->wmu[m];
	  Hm += pRG->l3imu[ifr][j][i][1][m] * pRG->mu[1][m][2] * pRG->wmu[m];
	  Hm += pRG->l3imu[ifr][j][i][2][m] * pRG->mu[2][m][2] * pRG->wmu[m];
	  Hm += pRG->l3imu[ifr][j][i][3][m] * pRG->mu[3][m][2] * pRG->wmu[m];
	  Jm += pRG->l3imu[ifr][j][i][0][m] * pRG->wmu[m];
	  Jm += pRG->l3imu[ifr][j][i][1][m] * pRG->wmu[m];	  
	  Jm += pRG->l3imu[ifr][j][i][2][m] * pRG->wmu[m];
	  Jm += pRG->l3imu[ifr][j][i][3][m] * pRG->wmu[m];
	}	
	H = pRG->R[ifr][kl][j][i].H[2];
	I0 = (H - Hm) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][kl][j][i].J = 0.5 * I0 + Jm;
	for (m=0; m<nang; m++) {
	  pRG->Ghstl3i[ifr][j][i][0][m] = I0;
	  pRG->Ghstl3i[ifr][j][i][1][m] = I0;
	  pRG->Ghstl3i[ifr][j][i][2][m] = I0;
	  pRG->Ghstl3i[ifr][j][i][3][m] = I0;
	}
      }}
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic.
 * Outer x3 boundary  */

static void const_flux_ox3(GridS *pG, RadGridS *pRG, int ifs, int ife)
{

  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ku = pRG->ke+1;
  int nang = pRG->nang;
  int i, j, l, m, n, ifr;
  Real I0, Jp, H, Hp, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  gamma *= 4.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstr3i using r3imu */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	Hp = 0.0;
	Jp = 0.0;
	for (m=0; m<nang; m++) {
	  Hp += pRG->r3imu[ifr][j][i][4][m] * pRG->mu[4][m][2] * pRG->wmu[m];
	  Hp += pRG->r3imu[ifr][j][i][5][m] * pRG->mu[5][m][2] * pRG->wmu[m];
	  Hp += pRG->r3imu[ifr][j][i][6][m] * pRG->mu[6][m][2] * pRG->wmu[m];
	  Hp += pRG->r3imu[ifr][j][i][7][m] * pRG->mu[7][m][2] * pRG->wmu[m];
	  Jp += pRG->r3imu[ifr][j][i][4][m] * pRG->wmu[m];
	  Jp += pRG->r3imu[ifr][j][i][5][m] * pRG->wmu[m];
	  Jp += pRG->r3imu[ifr][j][i][6][m] * pRG->wmu[m];
	  Jp += pRG->r3imu[ifr][j][i][7][m] * pRG->wmu[m];	  
	}

	H = pRG->R[ifr][ku][j][i].H[2];
	I0 = (Hp - H) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][ku][j][i].J = 0.5 * I0 + Jp;
	for (m=0; m<nang; m++) {
	  pRG->Ghstr3i[ifr][j][i][4][m] = I0;
	  pRG->Ghstr3i[ifr][j][i][5][m] = I0;
	  pRG->Ghstr3i[ifr][j][i][6][m] = I0;
	  pRG->Ghstr3i[ifr][j][i][7][m] = I0;
	}
      }}
  }

  return;

}

/*----------------------------------------------------------------------------*/
/* Time independent incident radiation boundary condition.  Nothing is done 
 * here, as the incident boundary radiaion is specified at initialization and
 * unchanged during computation. This means that any contribution from back
 * scattering of outgoing radiation is ignored.  */

static void const_incident_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  return;
}


#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs */

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int is = pRG->is;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);
  
  for(ifr=ifs; ifr<=ife; ifr++) {
/*  send the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	*(pSnd++) = pRG->R[ifr][k][j][is].S;  
	*(pSnd++) = pRG->R[ifr][k][j][is].J;  
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][j][is].H[l];
	}
      }}
/* update Ghstr1i using l1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l1imu[ifr][k][j][1][m];
	  if(noct > 2) {
	    *(pSnd++) = pRG->l1imu[ifr][k][j][3][m];
	    if(noct == 8) {
	      *(pSnd++) = pRG->l1imu[ifr][k][j][5][m];
	      *(pSnd++) = pRG->l1imu[ifr][k][j][7][m];
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass Ghstl2i/Ghstr2i on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][jl][l][m];
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][ju][l][m];
	  }}}
    }
/* pass Ghstl3i/Ghstr3i on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][kl][j][l][m];
	    *(pSnd++) = pRG->Ghstr1i[ifr][ku][j][l][m];
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=1,3,5,7 are passed using l1imu and set in unpack routine */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][is][2][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][is][0][m];
	  if(noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][is][6][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][is][4][m];
	  }
	}}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=1,3,5,7 are passed using r1imu and set in unpack routine */
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][j][is][4][m];
	  *(pSnd++) = pRG->l3imu[ifr][j][is][6][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][is][0][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][is][2][m];	  
	}}
    }
#endif
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int ie = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	*(pSnd++) = pRG->R[ifr][k][j][ie].S;
	*(pSnd++) = pRG->R[ifr][k][j][ie].J;
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][j][ie].H[l];
	}
      }}
/* update Ghstl1i using r1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r1imu[ifr][k][j][l][m];
	  }}
#else
	for (m=0; m<nang; m++) {
	  *(pSnd++)  = pRG->r1imu[ifr][k][j][0][m];
	  if(noct > 2) {
	    *(pSnd++)  = pRG->r1imu[ifr][k][j][2][m];
	    if(noct == 8) {
	      *(pSnd++)  = pRG->r1imu[ifr][k][j][4][m];
	      *(pSnd++)  = pRG->r1imu[ifr][k][j][6][m];
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass Ghstl2i/Ghstr2i on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][jl][l][m];
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][ju][l][m];
	  }}}
    }
/* pass Ghstl3i/Ghstr3i on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][kl][j][l][m];
	    *(pSnd++) = pRG->Ghstl1i[ifr][ku][j][l][m];
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=0,2,4,6 are passed using r1imu and set in unpack routine */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][ie][3][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][ie][1][m];
	  if(noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][ie][7][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][ie][5][m];
	  }
	}}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=0,2,4,6 are passed using r1imu and set in unpack routine */
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][j][ie][5][m];
	  *(pSnd++) = pRG->l3imu[ifr][j][ie][7][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][ie][1][m];
	  *(pSnd++) = pRG->r3imu[ifr][j][ie][3][m];	   
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int js = pRG->js;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	*(pSnd++) = pRG->R[ifr][k][js][i].S;
	*(pSnd++) = pRG->R[ifr][k][js][i].J;
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][js][i].H[l];
	}
      }}
/* update Ghstr2i using l2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones. */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l2imu[ifr][k][i][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l2imu[ifr][k][i][2][m];
	  *(pSnd++) = pRG->l2imu[ifr][k][i][3][m];
	  if (noct == 8) {
	    *(pSnd++) = pRG->l2imu[ifr][k][i][6][m];
	    *(pSnd++) = pRG->l2imu[ifr][k][i][7][m];
	  }
	}
#endif
      }}
/* update Ghstr1i/Ghstl1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][js][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][js][l][m];
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl1i[ifr][k][js][0][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][k][js][2][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][k][js][1][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][k][js][3][m];
	if (noct == 8) {
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][js][4][m];
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][js][6][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][js][5][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][js][7][m];
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr2i[ifr][kl][i][l][m];
	    *(pSnd++) = pRG->Ghstr2i[ifr][ku][i][l][m];
	  }}}
    }
#else
/* Note that values l=2,3,6,7 are passed using r2imu */
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][js][i][4][m];
	  *(pSnd++) = pRG->l3imu[ifr][js][i][5][m];
	  *(pSnd++) = pRG->r3imu[ifr][js][i][0][m];
	  *(pSnd++) = pRG->r3imu[ifr][js][i][1][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int je = pRG->je;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	*(pSnd++) = pRG->R[ifr][k][je][i].S;
	*(pSnd++) = pRG->R[ifr][k][je][i].J;
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][je][i].H[l];
	}
      }}
/* update l2imu using r2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones. */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r2imu[ifr][k][i][l][m];
	  }}
#else
	for(m=0; m<nang; m++) {
	  *(pSnd++) = pRG->r2imu[ifr][k][i][0][m];
	  *(pSnd++) = pRG->r2imu[ifr][k][i][1][m];
	  if (noct == 8) {
	    *(pSnd++) = pRG->r2imu[ifr][k][i][4][m];
	    *(pSnd++) = pRG->r2imu[ifr][k][i][5][m];	    
	  }
	}
#endif
      }}
/* update Ghstr1i/Ghstl1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][je][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][je][l][m];
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl1i[ifr][k][je][0][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][k][je][2][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][k][je][1][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][k][je][3][m];
	if (noct == 8) {
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][je][4][m];
	  *(pSnd++) = pRG->Ghstl1i[ifr][k][je][6][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][je][5][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][k][je][7][m];
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl2i[ifr][kl][i][l][m];
	    *(pSnd++) = pRG->Ghstl2i[ifr][ku][i][l][m];
	  }}}
    }
#else
/* Note that values l=0,1,4,5 are updated using l2imu */
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  *(pSnd++) = pRG->l3imu[ifr][je][i][6][m];
	  *(pSnd++) = pRG->l3imu[ifr][je][i][7][m];
	  *(pSnd++) = pRG->r3imu[ifr][je][i][2][m];
	  *(pSnd++) = pRG->r3imu[ifr][je][i][3][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[0][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	*(pSnd++) = pRG->R[ifr][ks][j][i].S;
	*(pSnd++) = pRG->R[ifr][ks][j][i].J;
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][ks][j][i].H[l];
	}
      }}
/* update Ghstr3i using l3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones. */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<8; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l3imu[ifr][j][i][l][m];
	  }}
#else
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->l3imu[ifr][j][i][l][m];
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][ks][j][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][ks][j][l][m];
	}}}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl2i[ifr][ks][i][l][m];
	  *(pSnd++) = pRG->Ghstr2i[ifr][ks][i][l][m];
	}}}
#else
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl1i[ifr][ks][j][0][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ks][j][2][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ks][j][4][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ks][j][6][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ks][j][1][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ks][j][3][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ks][j][5][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ks][j][7][m];
      }}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl2i[ifr][ks][i][0][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ks][i][1][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ks][i][4][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ks][i][5][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ks][i][2][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ks][i][3][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ks][i][6][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ks][i][7][m];
      }}
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pSnd;

  pSnd = (double*)&(send_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* pass the source function and flux */
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	*(pSnd++) = pRG->R[ifr][ke][j][i].S;
	*(pSnd++) = pRG->R[ifr][ke][j][i].J;
	for(l=0; l < nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][ke][j][i].H[l];
	}
      }}
/* update Ghstl3i using r3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones. */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<8; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r3imu[ifr][j][i][l][m];
	  }}
#else
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->r3imu[ifr][j][i][l][m];
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl1i[ifr][ke][j][l][m];
	  *(pSnd++) = pRG->Ghstr1i[ifr][ke][j][l][m];
	}}}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  *(pSnd++) = pRG->Ghstl2i[ifr][ke][i][l][m];
	  *(pSnd++) = pRG->Ghstr2i[ifr][ke][i][l][m];
	}}}
#else
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl1i[ifr][ke][j][0][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ke][j][2][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ke][j][4][m];
	*(pSnd++) = pRG->Ghstl1i[ifr][ke][j][6][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ke][j][1][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ke][j][3][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ke][j][5][m];
	*(pSnd++) = pRG->Ghstr1i[ifr][ke][j][7][m];
      }}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	*(pSnd++) = pRG->Ghstl2i[ifr][ke][i][0][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ke][i][1][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ke][i][4][m];
	*(pSnd++) = pRG->Ghstl2i[ifr][ke][i][5][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ke][i][2][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ke][i][3][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ke][i][6][m];
	*(pSnd++) = pRG->Ghstr2i[ifr][ke][i][7][m];
      }}
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	pRG->R[ifr][k][j][il].S = *(pRcv++);
	pRG->R[ifr][k][j][il].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][j][il].H[l] = *(pRcv++);
	}
      }}
/* update Ghstl1i using r1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][j][l][m] = *(pRcv++);
	  }}
#else
	for (m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = *(pRcv++);
	  if(noct > 2) {
	    pRG->Ghstl1i[ifr][k][j][2][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->Ghstl1i[ifr][k][j][4][m] = *(pRcv++);
	      pRG->Ghstl1i[ifr][k][j][6][m] = *(pRcv++);
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][il][l][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][il][l][m] = *(pRcv++);
	  }}}
    }
/* pass l3imu/r3imu on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][il][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][j][il][l][m] = *(pRcv++);
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=0,2,4,6 are passed using r1imu and handled below */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][3][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][il][1][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][il][7][m] = *(pRcv++);
	      pRG->r2imu[ifr][k][il][5][m] = *(pRcv++);
	    }
	    }}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=0,2,4,6 are passed using r1imu */
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][il][5][m] = *(pRcv++);
	  pRG->l3imu[ifr][j][il][7][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][il][1][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][il][3][m] = *(pRcv++);
	}}
    }
/* copy duplicate edges/coners from corresponding Ghstl1i */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][il][2][m] = pRG->Ghstl1i[ifr][k][jl][2][m];
	  pRG->r2imu[ifr][k][il][0][m] = pRG->Ghstl1i[ifr][k][ju][0][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][il][6][m] = pRG->Ghstl1i[ifr][k][jl][6][m];
	      pRG->r2imu[ifr][k][il][4][m] = pRG->Ghstl1i[ifr][k][ju][4][m];
	    }
	}}
    }
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][il][4][m] = pRG->Ghstl1i[ifr][kl][j][4][m];
	  pRG->l3imu[ifr][j][il][6][m] = pRG->Ghstl1i[ifr][kl][j][6][m];
	  pRG->r3imu[ifr][j][il][0][m] = pRG->Ghstl1i[ifr][ku][j][0][m];
	  pRG->r3imu[ifr][j][il][2][m] = pRG->Ghstl1i[ifr][ku][j][2][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int j, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	pRG->R[ifr][k][j][iu].S = *(pRcv++);
	pRG->R[ifr][k][j][iu].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][j][iu].H[l] = *(pRcv++);
	  }
      }}
/* update r1imu using l1imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][l][m] = *(pRcv++);
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstr1i[ifr][k][j][1][m] = *(pRcv++);
	  if(noct > 2) {
	    pRG->Ghstr1i[ifr][k][j][3][m] = *(pRcv++);
	    if(noct == 8) {
	      pRG->Ghstr1i[ifr][k][j][5][m] = *(pRcv++);
	      pRG->Ghstr1i[ifr][k][j][7][m] = *(pRcv++);
	    }
	  }
	}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* pass l2imu/r2imu on the corngers[2D]/edge[3D] */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l2imu[ifr][k][iu][l][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][iu][l][m] = *(pRcv++);
	  }}}
    }
/* pass l3imu/r3imu on the edges.*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][iu][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][j][iu][l][m] = *(pRcv++);
	  }}}
    }
#else
/* pass l2imu/r2imu on the corngers[2D]/edge[3D].  Note that values
 * l=1,3,5,7 are updated using r1imu below */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][2][m] = *(pRcv++);
	  pRG->r2imu[ifr][k][iu][0][m] = *(pRcv++);
	  if(noct == 8) {
	    pRG->l2imu[ifr][k][iu][6][m] = *(pRcv++);
	    pRG->r2imu[ifr][k][iu][4][m] = *(pRcv++);
	  }
	}}
    }
/* pass l3imu/r3imu on the edges.  Note that values
 * l=1,3,5,7 are update using r1imu below*/
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][iu][4][m] = *(pRcv++);
	  pRG->l3imu[ifr][j][iu][6][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][iu][0][m] = *(pRcv++);
	  pRG->r3imu[ifr][j][iu][2][m] = *(pRcv++);
	}}
    }
/* copy duplicate edges/coners from corresponding Ghstr1i */
    if (noct > 2) {
      for (k=kl; k<=ku; k++) {
	for (m=0; m<nang; m++) {
	  pRG->l2imu[ifr][k][iu][3][m] = pRG->Ghstr1i[ifr][k][jl][3][m];
	  pRG->r2imu[ifr][k][iu][1][m] = pRG->Ghstr1i[ifr][k][ju][1][m];
	    if(noct == 8) {
	      pRG->l2imu[ifr][k][iu][7][m] = pRG->Ghstr1i[ifr][k][jl][7][m];
	      pRG->r2imu[ifr][k][iu][5][m] = pRG->Ghstr1i[ifr][k][ju][5][m];
	    }
	}}
    }
    if (noct == 8) {
      for (j=jl; j<=ju; j++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][j][iu][5][m] = pRG->Ghstr1i[ifr][kl][j][5][m];
	  pRG->l3imu[ifr][j][iu][7][m] = pRG->Ghstr1i[ifr][kl][j][7][m];
	  pRG->r3imu[ifr][j][iu][1][m] = pRG->Ghstr1i[ifr][ku][j][1][m];
	  pRG->r3imu[ifr][j][iu][3][m] = pRG->Ghstr1i[ifr][ku][j][3][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][k][jl][i].S = *(pRcv++);
	pRG->R[ifr][k][jl][i].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][jl][i].H[l] = *(pRcv++);
	}
      }}
/* update Ghstl2i using r2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones. */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl2i[ifr][k][i][l][m] = *(pRcv++);
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = *(pRcv++);
	  pRG->Ghstl2i[ifr][k][i][1][m] = *(pRcv++);
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = *(pRcv++);
	    pRG->Ghstl2i[ifr][k][i][5][m] = *(pRcv++);
	  }
	}
#endif
      }}
/* update Ghstr1i/Ghstl1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][jl][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][jl][l][m] = *(pRcv++);
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][k][jl][0][m] = *(pRcv++);
	pRG->Ghstl1i[ifr][k][jl][2][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][k][jl][1][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][k][jl][3][m] = *(pRcv++);
	if (noct == 8) {
	  pRG->Ghstl1i[ifr][k][jl][4][m] = *(pRcv++);
	  pRG->Ghstl1i[ifr][k][jl][6][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][jl][5][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][jl][7][m] = *(pRcv++);
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][jl][i][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][jl][i][l][m] = *(pRcv++);
	  }}}
    }
#else
/* Note that values l=0,1,4,5 are updated using l2imu below*/
    if (noct == 8) {    
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][jl][i][6][m] = *(pRcv++);
	  pRG->l3imu[ifr][jl][i][7][m] = *(pRcv++);
	  pRG->r3imu[ifr][jl][i][2][m] = *(pRcv++);
	  pRG->r3imu[ifr][jl][i][3][m] = *(pRcv++);
	}}
  }
/* copy duplicate edges from corresponding l2imu */
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][jl][i][4][m] = pRG->Ghstl2i[ifr][kl][i][4][m];
	  pRG->l3imu[ifr][jl][i][5][m] = pRG->Ghstl2i[ifr][kl][i][5][m];
	  pRG->r3imu[ifr][jl][i][0][m] = pRG->Ghstl2i[ifr][ku][i][0][m];
	  pRG->r3imu[ifr][jl][i][1][m] = pRG->Ghstl2i[ifr][ku][i][1][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks,   ku = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][k][ju][i].S = *(pRcv++);
	pRG->R[ifr][k][ju][i].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][k][ju][i].H[l] = *(pRcv++);
	}
      }}
/* update Ghstr2i using l2imu. Note that i runs from
 * is-1 to ie+1 so loop includes ghost zones. */
    for (k=kl; k<=ku; k++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for (l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr2i[ifr][k][i][l][m] = *(pRcv++);
	  }}
#else
	for(m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][k][i][2][m] = *(pRcv++);
	  pRG->Ghstr2i[ifr][k][i][3][m] = *(pRcv++);
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = *(pRcv++);
	    pRG->Ghstr2i[ifr][k][i][7][m] = *(pRcv++);
	  }
	}
#endif
      }}
/* update Ghstr1i/Ghstl1i on corners[2d]/edges[3d] */
#if defined(QUADRATIC_INTENSITY) || defined(SHEARING_BOX)
    for (k=kl; k<=ku; k++) {
      for (l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][k][ju][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][ju][l][m] = *(pRcv++);
	}}}
#else
    for (k=kl; k<=ku; k++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][k][ju][0][m] = *(pRcv++);
	pRG->Ghstl1i[ifr][k][ju][2][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][k][ju][1][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][k][ju][3][m] = *(pRcv++);
	if (noct == 8) {
	  pRG->Ghstl1i[ifr][k][ju][4][m] = *(pRcv++);
	  pRG->Ghstl1i[ifr][k][ju][6][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][ju][5][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][k][ju][7][m] = *(pRcv++);
	}
      }}
#endif
/* pass l3imu/r3imu on the edges. */
#ifdef QUADRATIC_INTENSITY
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for(l=0; l<noct; l++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][ju][i][l][m] = *(pRcv++);
	    pRG->r3imu[ifr][ju][i][l][m] = *(pRcv++);
	  }}}
    }
#else
/* Note that values l=2,3,6,7 are updated using r2imu below */
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][ju][i][4][m] = *(pRcv++);
	  pRG->l3imu[ifr][ju][i][5][m] = *(pRcv++);
	  pRG->r3imu[ifr][ju][i][0][m] = *(pRcv++);
	  pRG->r3imu[ifr][ju][i][1][m] = *(pRcv++);
	}}
    }
/* copy duplicate edges from corresponding Ghstl2i */
    if (noct == 8) {
      for (i=il; i<=iu; i++) {
	for (m=0; m<nang; m++) {
	  pRG->l3imu[ifr][ju][i][6][m] = pRG->Ghstr2i[ifr][kl][i][6][m];
	  pRG->l3imu[ifr][ju][i][7][m] = pRG->Ghstr2i[ifr][kl][i][7][m];
	  pRG->r3imu[ifr][ju][i][2][m] = pRG->Ghstr2i[ifr][ku][i][2][m];
	  pRG->r3imu[ifr][ju][i][3][m] = pRG->Ghstr2i[ifr][ku][i][3][m];
	}}
    }
#endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int kl = pRG->ks-1;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[0][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (j=jl; j<=ju; j++) { 
      for (i=il; i<=iu; i++) {
	pRG->R[ifr][kl][j][i].S = *(pRcv++);
	pRG->R[ifr][kl][j][i].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][kl][j][i].H[l] = *(pRcv++);
	}
      }}
/* update Ghstl3i using r3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones.*/
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<8; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = *(pRcv++);
	  }}
#else
	for(l=0; l<=3; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][l][m] = *(pRcv++);
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][kl][j][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][kl][j][l][m] = *(pRcv++);
	}}}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][kl][i][l][m] = *(pRcv++);
	  pRG->Ghstr2i[ifr][kl][i][l][m] = *(pRcv++);
	}}}
#else
/* update Ghstr1i/Ghstl1i on edges */
    for (j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl1i[ifr][kl][j][0][m] = *(pRcv++);
	pRG->Ghstl1i[ifr][kl][j][2][m] = *(pRcv++);
	pRG->Ghstl1i[ifr][kl][j][4][m] = *(pRcv++);
	pRG->Ghstl1i[ifr][kl][j][6][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][kl][j][1][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][kl][j][3][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][kl][j][5][m] = *(pRcv++);
	pRG->Ghstr1i[ifr][kl][j][7][m] = *(pRcv++);
      }}
/* update Ghstr2i/Ghstl2i on edges */
    for (i=il; i<=iu; i++) {
      for(m=0; m<nang; m++) {	
	pRG->Ghstl2i[ifr][kl][i][0][m] = *(pRcv++);
	pRG->Ghstl2i[ifr][kl][i][1][m] = *(pRcv++);
	pRG->Ghstl2i[ifr][kl][i][4][m] = *(pRcv++);
	pRG->Ghstl2i[ifr][kl][i][5][m] = *(pRcv++);
	pRG->Ghstr2i[ifr][kl][i][2][m] = *(pRcv++);
	pRG->Ghstr2i[ifr][kl][i][3][m] = *(pRcv++);
	pRG->Ghstr2i[ifr][kl][i][6][m] = *(pRcv++);
	pRG->Ghstr2i[ifr][kl][i][7][m] = *(pRcv++);
      }}
 #endif
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3_rad(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int jl = pRG->js-1, ju = pRG->je+1;
  int ku = pRG->ke+1;
  int nf = pRG->nf, nang = pRG->nang;
  int i, j, l, m, ifr;
  double *pRcv;

  pRcv = (double*)&(recv_buf[1][0]);

  for(ifr=ifs; ifr<=ife; ifr++) {
/* receive the source function and flux */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) { 
	pRG->R[ifr][ku][j][i].S = *(pRcv++);
	pRG->R[ifr][ku][j][i].J = *(pRcv++);
	for(l=0; l < nDim; l++) {
	  pRG->R[ifr][ku][j][i].H[l] = *(pRcv++);
	}
      }}
/* update Ghstr3i using l3imu. Note that i runs from  is-1 to ie+1 and 
 * j runs from js-1 to je+1 so loop includes ghost zones. */
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
#ifdef QUADRATIC_INTENSITY
	for(l=0; l<8; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr3i[ifr][j][i][l][m] = *(pRcv++);
	  }}
#else
	for(l=4; l<=7; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr3i[ifr][j][i][l][m] = *(pRcv++);
	  }}
#endif
      }}
#ifdef QUADRATIC_INTENSITY
/* update r1imu/l1imu on edges */
    for (j=jl; j<=ju; j++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl1i[ifr][ku][j][l][m] = *(pRcv++);
	  pRG->Ghstr1i[ifr][ku][j][l][m] = *(pRcv++);
	}}}
/* update r2imu/l2imu on edges */
    for (i=il; i<=iu; i++) {
      for (l=0; l<8; l++) {
	for(m=0; m<nang; m++) {	
	  pRG->Ghstl2i[ifr][ku][i][l][m] = *(pRcv++);
	  pRG->Ghstr2i[ifr][ku][i][l][m] = *(pRcv++);
	}}}
#else
/* update r1imu/l1imu on edges */
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
      }}
/* update r2imu/l2imu on edges */
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
      }}
#endif
  }
  return;
}

#endif /* MPI_PARALLEL */

#endif /* RADIATION_TRANSFER */
