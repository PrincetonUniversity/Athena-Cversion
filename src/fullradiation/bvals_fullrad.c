#include "../copyright.h"
/*=============================================================================
 * FILE: bvals_fullrad.c
 * copy from file bvals_fullrad.c
 * PURPOSE: Sets boundary condistions for radiative transfer on each edge of
 *          the grid.  It closely follows the methods and conventions and
 *          methods used for the hydro integration (see e.g. bvals_mhd.c).
 *          The radiation source function, radiative flux, and intensities
 *          on the boundaries of the RadGrid are copied to ghost zones.
 *
 *  Need to update boundary conditions for each specific intensity
 *  the source terms heatcool and Scat 
 * CONTAINS PUBLIC FUNCTIONS: 
 *   bvals_fullrad_init()
 *   bvals_fullrad()
 */

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER

static int nDim;
static int Rotate_flag;
#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
static  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */

#endif /* MPI_PARALLEL */
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   periodic_??_fullrad?() - periodic BCs at boundary ???
 *   pack_???_fullrad()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???_fullrad()   - unpack data for MPI non-blocking receive at ??? boundary
 *============================================================================*/
static void periodic_ix1_fullrad(GridS *pG, RadGridS *pRG);
static void periodic_ox1_fullrad(GridS *pG, RadGridS *pRG);
static void periodic_ix2_fullrad(GridS *pG, RadGridS *pRG);
static void periodic_ox2_fullrad(GridS *pG, RadGridS *pRG);
static void periodic_ix3_fullrad(GridS *pG, RadGridS *pRG);
static void periodic_ox3_fullrad(GridS *pG, RadGridS *pRG);
static void Rotate180_ix2_fullrad(RadGridS *pRG);
static void Rotate180_ox2_fullrad(RadGridS *pRG);
static void Rotate90_ix2_fullrad(RadGridS *pRG);
static void Rotate90_ox2_fullrad(RadGridS *pRG);

static void outflow_ix1_fullrad(GridS *pG, RadGridS *pRG);
static void outflow_ox1_fullrad(GridS *pG, RadGridS *pRG);
static void outflow_ix2_fullrad(GridS *pG, RadGridS *pRG);
static void outflow_ox2_fullrad(GridS *pG, RadGridS *pRG);
static void outflow_ix3_fullrad(GridS *pG, RadGridS *pRG);
static void outflow_ox3_fullrad(GridS *pG, RadGridS *pRG);

static void vacuum_ix1_fullrad(GridS *pG, RadGridS *pRG);
static void vacuum_ox1_fullrad(GridS *pG, RadGridS *pRG);
static void vacuum_ix2_fullrad(GridS *pG, RadGridS *pRG);
static void vacuum_ox2_fullrad(GridS *pG, RadGridS *pRG);
static void vacuum_ix3_fullrad(GridS *pG, RadGridS *pRG);
static void vacuum_ox3_fullrad(GridS *pG, RadGridS *pRG);



static void ProlongateLater(GridS *pG, RadGridS *pRG);


static void const_flux_ix1(GridS *pG, RadGridS *pRG);
static void const_flux_ox1(GridS *pG, RadGridS *pRG);
static void const_flux_ix2(GridS *pG, RadGridS *pRG);
static void const_flux_ox2(GridS *pG, RadGridS *pRG);
static void const_flux_ix3(GridS *pG, RadGridS *pRG);
static void const_flux_ox3(GridS *pG, RadGridS *pRG);

#ifdef MPI_PARALLEL
static void pack_ix1_fullrad(GridS *pG, RadGridS *pRG);
static void pack_ox1_fullrad(GridS *pG, RadGridS *pRG);
static void pack_ix2_fullrad(GridS *pG, RadGridS *pRG);
static void pack_ox2_fullrad(GridS *pG, RadGridS *pRG);
static void pack_ix3_fullrad(GridS *pG, RadGridS *pRG);
static void pack_ox3_fullrad(GridS *pG, RadGridS *pRG);

static void unpack_ix1_fullrad(GridS *pG, RadGridS *pRG);
static void unpack_ox1_fullrad(GridS *pG, RadGridS *pRG);
static void unpack_ix2_fullrad(GridS *pG, RadGridS *pRG);
static void unpack_ox2_fullrad(GridS *pG, RadGridS *pRG);
static void unpack_ix3_fullrad(GridS *pG, RadGridS *pRG);
static void unpack_ox3_fullrad(GridS *pG, RadGridS *pRG);
#endif /* MPI_PARALLEL */

void bvals_fullrad(DomainS *pD)
{
  RadGridS *pRG=(pD->RadGrid);
  GridS *pG=(pD->Grid);
  int myL,myM,myN;
#ifdef SHEARING_BOX
  int BCFlag;
#endif
#ifdef MPI_PARALLEL
  int cnt, ierr, mIndex;
#endif /* MPI_PARALLEL */


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
      pack_ix1_fullrad(pG,pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox1_fullrad(pG,pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox1_fullrad(pG,pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx1_id >= 0 && pRG->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));
      /* set physical boundary */
      (*(pD->ix1_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1_fullrad(pG,pRG);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx1_id < 0 && pRG->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox1_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1_fullrad(pG,pRG);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx1_id < 0 && pRG->lx1_id < 0) {
      (*(pD->ix1_RBCFun))(pG,pRG);
      (*(pD->ox1_RBCFun))(pG,pRG);
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
      pack_ix2_fullrad(pG,pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox2_fullrad(pG,pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox2_fullrad(pG,pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx2_id >= 0 && pRG->lx2_id < 0) {
      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix2_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2_fullrad(pG,pRG);


    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx2_id < 0 && pRG->lx2_id >= 0) {
      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox2_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2_fullrad(pG,pRG);

    }
#endif /* MPI_PARALLEL */


/* Physical boundaries on both left and right */
    if (pRG->rx2_id < 0 && pRG->lx2_id < 0) {
      (*(pD->ix2_RBCFun))(pG,pRG);
      (*(pD->ox2_RBCFun))(pG,pRG);
    } 
	  

/*-------------------------------------------*/
/* Check whether we need to rotate the angles or not */
/* Only do this if Rotate_flag > 0 && the domain is near the J physical boundary */	  
/* Only rotate the angles for the ghost zones */

	/* Rotate_flag==1, rotate 90 degree */
	/* Rotate_flag==2, rotate 180 degree */

	get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
	if(myM == 0){
		if(Rotate_flag == 1)
			Rotate90_ix2_fullrad(pRG);
		else if(Rotate_flag == 2)
			Rotate180_ix2_fullrad(pRG);
	}
	  
	if (myM == ((pD->NGrid[1])-1)) {
		if(Rotate_flag == 1)
			Rotate90_ox2_fullrad(pRG);
		else if(Rotate_flag == 2)
			Rotate180_ox2_fullrad(pRG);
	}


	  
/*-------------------------------------------*/
	  
 
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
      pack_ix3_fullrad(pG,pRG);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox3_fullrad(pG,pRG);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3_fullrad(pG,pRG);
      if (mIndex == 1) unpack_ox3_fullrad(pG,pRG);

    }

/* Physical boundary on left, MPI block on right */
    if (pRG->rx3_id >= 0 && pRG->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pRG->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->ix3_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3_fullrad(pG,pRG);

    }

/* MPI block on left, Physical boundary on right */
    if (pRG->rx3_id < 0 && pRG->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3_fullrad(pG,pRG); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pRG->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->ox3_RBCFun))(pG,pRG);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3_fullrad(pG,pRG);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pRG->rx3_id < 0 && pRG->lx3_id < 0) {
      (*(pD->ix3_RBCFun))(pG,pRG);
      (*(pD->ox3_RBCFun))(pG,pRG);
    } 

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* bvals_fullrad_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI.
 *   Patterned closely after bvals_mhd_init().
 */

void bvals_fullrad_init(MeshS *pM)
{

  RadGridS *pRG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  /*int x1cnt=0, x2cnt=0, x3cnt=0; */ 
/* Number of words passed in x1/x2/x3-dir. */
  int nang, nf, noct, xcnt;
#endif /* MPI_PARALLEL */
	
	/* Default, no rotation for the angles */
	/* this flag only works for j boundary flag > 6 */
	Rotate_flag = 0;

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
	    pD->ix1_RBCFun = outflow_ix1_fullrad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix1_RBCFun = periodic_ix1_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->lx1_id < 0 && pD->NGrid[0] > 1){
	      pRG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix1_RBCFun = const_flux_ix1;
	    break;

	  case 6:
	    pD->ix1_RBCFun = vacuum_ix1_fullrad;
	    break;

	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ix1=%d unknown\n",pM->RBCFlag_ix1);
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
	    pD->ox1_RBCFun = outflow_ox1_fullrad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox1_RBCFun = periodic_ox1_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->rx1_id < 0 && pD->NGrid[0] > 1){
	      pRG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;
	    
	  case 5: /* constant fixed flux at boundary */
	    pD->ox1_RBCFun = const_flux_ox1;
	    break;
	
	  case 6:
	    pD->ox1_RBCFun = vacuum_ox1_fullrad;
	    break;

	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ox1=%d unknown\n",pM->RBCFlag_ox1);
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
	    pD->ix2_RBCFun = outflow_ix2_fullrad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix2_RBCFun = periodic_ix2_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->lx2_id < 0 && pD->NGrid[1] > 1){
	      pRG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix2_RBCFun = const_flux_ix2;
	    break;

	  case 6:
	    pD->ix2_RBCFun = vacuum_ix2_fullrad;
	    break;
			  
	  case 7: /* Rotate 90 degree */
		pD->ix2_RBCFun = periodic_ix2_fullrad;
	    Rotate_flag = 1;
#ifdef MPI_PARALLEL
	    if(pRG->lx2_id < 0 && pD->NGrid[1] > 1){
		  pRG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
		}
#endif /* MPI_PARALLEL */
	  case 8: /* Rotate 180 degree */
                pD->ix2_RBCFun = periodic_ix2_fullrad;
            Rotate_flag = 2;
#ifdef MPI_PARALLEL
            if(pRG->lx2_id < 0 && pD->NGrid[1] > 1){
                  pRG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
                }
#endif /* MPI_PARALLEL */
		  
	    break;
	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ix2=%d unknown\n",pM->RBCFlag_ix2);
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
	    pD->ox2_RBCFun = outflow_ox2_fullrad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox2_RBCFun = periodic_ox2_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->rx2_id < 0 && pD->NGrid[1] > 1){
	      pRG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ox2_RBCFun = const_flux_ox2;
	    break;

	  case 6:
	    pD->ox2_RBCFun = vacuum_ox2_fullrad;
	    break;
			  
	  case 7:/* Rotate 90 degree */
	    pD->ox2_RBCFun = periodic_ox2_fullrad;
	    Rotate_flag = 1;
			  
#ifdef MPI_PARALLEL
		if(pRG->rx2_id < 0 && pD->NGrid[1] > 1){
		  pRG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
		}
#endif /* MPI_PARALLEL */
	  case 8: /* Rotate 180 degree */
            pD->ox2_RBCFun = periodic_ox2_fullrad;
            Rotate_flag = 2;

#ifdef MPI_PARALLEL
                if(pRG->rx2_id < 0 && pD->NGrid[1] > 1){
                  pRG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
                }
#endif /* MPI_PARALLEL */

	  break;

	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ox2=%d unknown\n",pM->RBCFlag_ox2);
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
	    pD->ix3_RBCFun = outflow_ix3_fullrad;
	    break;

	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ix3_RBCFun = periodic_ix3_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->lx3_id < 0 && pD->NGrid[2] > 1){
	      pRG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ix3_RBCFun = const_flux_ix3;
	    break;

	  case 6:
	    pD->ix3_RBCFun = vacuum_ix3_fullrad;
	    break;

	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ix3=%d unknown\n",pM->RBCFlag_ix3);
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
	    pD->ox3_RBCFun = outflow_ox3_fullrad;
	    break;
	    
	  case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
	    pD->ox3_RBCFun = periodic_ox3_fullrad;
#ifdef MPI_PARALLEL
	    if(pRG->rx3_id < 0 && pD->NGrid[2] > 1){
	      pRG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
	    }
#endif /* MPI_PARALLEL */
	    break;

	  case 5: /* constant fixed flux at boundary */
	    pD->ox3_RBCFun = const_flux_ox3;
	    break;

	  case 6:
	    pD->ox3_RBCFun = vacuum_ox3_fullrad;
	    break;

	  default:
	    ath_perr(-1,"[bvals_fullrad_init]:rbc_ox3=%d unknown\n",pM->RBCFlag_ox3);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }

#ifdef MPI_PARALLEL 
    if(myID_Comm_world == 0){
#endif
        /* Print out the parameters for consistency check */
       printf("Azimuthal boundary flag: Rotateflag:  %d\n",Rotate_flag);
#ifdef MPI_PARALLEL 
    }
#endif


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

	  if(nx2t > 1){
		nx2t += 2 * Radghost;
	  }
	  if(nx3t > 1){
		nx3t += 2 * Radghost;
	  }

	  /* need to set each array along each direction */
	  /* space for J H, K and Sigma */
	  xcnt = Radghost * nx2t * nx3t * nf * (1 + 3 + 6 + 4);

	  xcnt += Radghost * nx2t * nx3t * noct * nang * nf;

	  if(xcnt > x1cnt) x1cnt = xcnt;
	}

/* x2cnt is the number of Reals passed for x2 faces */
	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2 * Radghost;
	  nx3t = pD->GData[n][m][l].Nx[2];

	  if(nx3t > 1){
		nx3t += 2 * Radghost;
	  }

	   /* space for J H, K and Sigma */
	  xcnt = Radghost * nx1t * nx3t * nf * (1 + 3 + 6 + 4);

	  xcnt += Radghost * nx1t * nx3t * noct * nang * nf;

          if(xcnt > x2cnt) x2cnt = xcnt;
	}

/* x3cnt is the number of Reals passed for x3 faces */
	if(pD->NGrid[2] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0] + 2 * Radghost;
	  nx2t = pD->GData[n][m][l].Nx[1] + 2 * Radghost;


	   /* space for J H, K and Sigma */
	  xcnt = Radghost * nx1t * nx2t * nf * (1 + 3 + 6 + 4);

	  xcnt += Radghost * nx1t * nx2t * noct * nang * nf;

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
      ath_error("[bvals_fullrad_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double**)calloc_2d_array(2,size,sizeof(double))) == NULL)
      ath_error("[bvals_fullrad_init]: Failed to allocate recv buffer\n");
  }

  if((recv_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_fullrad_init]: Failed to allocate recv MPI_Request array\n");
  if((send_rq = (MPI_Request*) calloc_1d_array(2,sizeof(MPI_Request))) == NULL)
    ath_error("[bvals_fullrad_init]: Failed to allocate send MPI_Request array\n");

#endif /* MPI_PARALLEL */


  return;
}


void bvals_fullrad_destruct()
{


#ifdef MPI_PARALLEL
	if(send_buf != NULL)
		free_2d_array(send_buf);

	if(recv_buf != NULL)
		free_2d_array(recv_buf);

	if(recv_rq != NULL)
		free_1d_array(recv_rq);

	if(send_rq != NULL)
		free_1d_array(send_rq);

#endif

}


void bvals_fullrad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc)
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
    ath_perr(-1,"[bvals_fullrad_trans_fun]: Unknown direction = %d\n",dir);
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
 *  
 */


/*--------------------------------------------------------------------------*/
/* OUTFLOW boundary condition */


static void outflow_ix1_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;

  /* set the angle independent values */

for(ifr=0; ifr<nf; ifr++) {

  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][is-i] = pRG->R[ifr][k][j][is];
	}
     }
  }
  
    for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
	
			pRG->imu[ifr][l][n][k][j][is-i] = pRG->imu[ifr][l][n][k][j][is];			
			pRG->heatcool[ifr][l][n][k][j][is-i] = pRG->heatcool[ifr][l][n][k][j][is];

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */
   }/* end ifr */
	  


  return;
}




static void outflow_ox1_fullrad(GridS *pG, RadGridS *pRG)
{

  int ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;

  /* set the angle independent values */

for(ifr=0; ifr<nf; ifr++) {

  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][ie+i] = pRG->R[ifr][k][j][ie];
	}
     }
  }
  
    for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
	
			pRG->imu[ifr][l][n][k][j][ie+i] = pRG->imu[ifr][l][n][k][j][ie];			
			pRG->heatcool[ifr][l][n][k][j][ie+i] = pRG->heatcool[ifr][l][n][k][j][ie];

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */
   }/* end ifr */
	  


  return;
}



static void outflow_ix2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			/* set the angle independent values */
			pRG->R[ifr][k][js-j][i] = pRG->R[ifr][k][js][i];
     	   }/* end i */
       }/* end J */
    } /* End k */

   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			pRG->imu[ifr][l][n][k][js-j][i] = pRG->imu[ifr][l][n][k][js][i];			
			pRG->heatcool[ifr][l][n][k][js-j][i] = pRG->heatcool[ifr][l][n][k][js][i]; 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */

  return;
}



static void outflow_ox2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			/* set the angle independent values */
			pRG->R[ifr][k][je+j][i] = pRG->R[ifr][k][je][i];
     	   }/* end i */
       }/* end J */
    } /* End k */

   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			pRG->imu[ifr][l][n][k][je+j][i] = pRG->imu[ifr][l][n][k][je][i];			
			pRG->heatcool[ifr][l][n][k][je+j][i] = pRG->heatcool[ifr][l][n][k][je][i]; 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */

  return;
}


static void outflow_ix3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=1; k<=Radghost; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][ks-k][j][i] = pRG->R[ifr][ks][j][i];
	   }/* end i */
       }/* end J */
    } /* End k */
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			pRG->imu[ifr][l][n][ks-k][j][i] = pRG->imu[ifr][l][n][ks][j][i];			
			pRG->heatcool[ifr][l][n][ks-k][j][i] = pRG->heatcool[ifr][l][n][ks][j][i];
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 

  return;
}


static void outflow_ox3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=1; k<=Radghost; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][ke+k][j][i] = pRG->R[ifr][ke][j][i];
	   }/* end i */
       }/* end J */
    } /* End k */
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			pRG->imu[ifr][l][n][ke+k][j][i] = pRG->imu[ifr][l][n][ke][j][i];			
			pRG->heatcool[ifr][l][n][ke+k][j][i] = pRG->heatcool[ifr][l][n][ke][j][i];
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 

  return;
}

/*--------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (rbc_ix1=1) */
/* We need to set the values in the corners of the boundary , which is needed */

static void periodic_ix1_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;

  /* set the angle independent values */

for(ifr=0; ifr<nf; ifr++) {

  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][is-i] = pRG->R[ifr][k][j][ie-(i-1)];
	}
     }
  }
  
    for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
	
			pRG->imu[ifr][l][n][k][j][is-i] = pRG->imu[ifr][l][n][k][j][ie-(i-1)];			
			pRG->heatcool[ifr][l][n][k][j][is-i] = pRG->heatcool[ifr][l][n][k][j][ie-(i-1)];

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */
   }/* end ifr */
	  


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (rbc_ox1=1) */

static void periodic_ox1_fullrad(GridS *pG, RadGridS *pRG)
{


  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


  

for(ifr=0; ifr<nf; ifr++) {

  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][ie+i] = pRG->R[ifr][k][j][is+(i-1)];
	}
     }
  }
  
    for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
	
			pRG->imu[ifr][l][n][k][j][ie+i] = pRG->imu[ifr][l][n][k][j][is+(i-1)];			
			pRG->heatcool[ifr][l][n][k][j][ie+i] = pRG->heatcool[ifr][l][n][k][j][is+(i-1)];

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */
   }/* end ifr */
	  
  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=1) */
/* We need to set the corner values */

static void periodic_ix2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			/* set the angle independent values */
			pRG->R[ifr][k][js-j][i] = pRG->R[ifr][k][je-(j-1)][i];
     	   }/* end i */
       }/* end J */
    } /* End k */

   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			pRG->imu[ifr][l][n][k][js-j][i] = pRG->imu[ifr][l][n][k][je-(j-1)][i];			
			pRG->heatcool[ifr][l][n][k][js-j][i] = pRG->heatcool[ifr][l][n][k][je-(j-1)][i]; 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=1) */

static void periodic_ox2_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;

for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			/* set the angle independent values */
			pRG->R[ifr][k][je+j][i] = pRG->R[ifr][k][js+(j-1)][i];
     	   }/* end i */
       }/* end J */
    } /* End k */

   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			pRG->imu[ifr][l][n][k][je+j][i] = pRG->imu[ifr][l][n][k][js+(j-1)][i];			
			pRG->heatcool[ifr][l][n][k][je+j][i] = pRG->heatcool[ifr][l][n][k][js+(j-1)][i]; 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (rbc_ix3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ix3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=1; k<=Radghost; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][ks-k][j][i] = pRG->R[ifr][ke-(k-1)][j][i];
	   }/* end i */
       }/* end J */
    } /* End k */
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			pRG->imu[ifr][l][n][ks-k][j][i] = pRG->imu[ifr][l][n][ke-(k-1)][j][i];			
			pRG->heatcool[ifr][l][n][ks-k][j][i] = pRG->heatcool[ifr][l][n][ke-(k-1)][j][i];
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (rbc_ox3=1) */
/* MUST BE MODIFIED TO INCLUDE RADIATION ONCE 3D RAD IS IMPLEMENTED */

static void periodic_ox3_fullrad(GridS *pG, RadGridS *pRG)
{
  
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    for(k=1; k<=Radghost; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][ke+k][j][i] = pRG->R[ifr][ks+(k-1)][j][i];
	   }/* end i */
       }/* end J */
    } /* End k */
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			pRG->imu[ifr][l][n][ke+k][j][i] = pRG->imu[ifr][l][n][ks+(k-1)][j][i];			
			pRG->heatcool[ifr][l][n][ke+k][j][i] = pRG->heatcool[ifr][l][n][ks+(k-1)][j][i];
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 



  return;
}


/*----------------------------------------------------------------------------*/
/* Rotate PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=7) */
/* Angles for specific intensities rotates for 90 degrees  */

static void Rotate90_ix2_fullrad(RadGridS *pRG)
{
	int is = pRG->is, ie = pRG->ie;
	int js = pRG->js;
	int ks = pRG->ks, ke = pRG->ke;
	int nang = pRG->nang;
	int nf = pRG->nf;
	int i, j, k, n, ifr;
	Real swap;
	
	
	for(ifr=0; ifr<nf; ifr++) {		
		for(n=0; n<nang; n++){
			for(k=ks; k<=ke; k++){
				for(j=1; j<=Radghost; j++){
					for(i=is-Radghost; i<=ie+Radghost; i++){
						/* Swap l = 0 && l =1, l=2 and l=3 */
						swap = pRG->imu[ifr][0][n][k][js-j][i];
						pRG->imu[ifr][0][n][k][js-j][i] = pRG->imu[ifr][1][n][k][js-j][i];
						pRG->imu[ifr][1][n][k][js-j][i] = swap;
						
						swap = pRG->imu[ifr][2][n][k][js-j][i];
						pRG->imu[ifr][2][n][k][js-j][i] = pRG->imu[ifr][3][n][k][js-j][i];
						pRG->imu[ifr][3][n][k][js-j][i] = swap;
						
						/* If for 3D */ 
						/* swap l ==4 && l==5, l==6 and l==7 */
						if(pRG->noct > 4){
							swap = pRG->imu[ifr][4][n][k][js-j][i];
							pRG->imu[ifr][4][n][k][js-j][i] = pRG->imu[ifr][5][n][k][js-j][i];
							pRG->imu[ifr][5][n][k][js-j][i] = swap;
							
							swap = pRG->imu[ifr][6][n][k][js-j][i];
							pRG->imu[ifr][6][n][k][js-j][i] = pRG->imu[ifr][7][n][k][js-j][i];
							pRG->imu[ifr][7][n][k][js-j][i] = swap;
						}
						
					}/* end i */
				}/* end J */
			} /* End k */
		}/* end nang */ 
	}/* end ifr */
	
	return;
}

/*----------------------------------------------------------------------------*/
/* Rotate PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=7) */

static void Rotate90_ox2_fullrad(RadGridS *pRG)
{
	
	int is = pRG->is, ie = pRG->ie;
	int je = pRG->je;
	int ks = pRG->ks, ke = pRG->ke;
	int nang = pRG->nang;
	int nf = pRG->nf;
	int i, j, k, n, ifr;
	Real swap;
	
	
	for(ifr=0; ifr<nf; ifr++) {		
		for(n=0; n<nang; n++){
			for(k=ks; k<=ke; k++){
				for(j=1; j<=Radghost; j++){
					for(i=is-Radghost; i<=ie+Radghost; i++){
						/* Swap l = 0 && l =1 , l==2 and l ==3*/
						swap = pRG->imu[ifr][0][n][k][je+j][i];
						pRG->imu[ifr][0][n][k][je+j][i] = pRG->imu[ifr][1][n][k][je+j][i];
						pRG->imu[ifr][1][n][k][je+j][i] = swap;
						
						swap = pRG->imu[ifr][2][n][k][je+j][i];
						pRG->imu[ifr][2][n][k][je+j][i] = pRG->imu[ifr][3][n][k][je+j][i];
						pRG->imu[ifr][3][n][k][je+j][i] = swap;
						
						/* swap l ==4 && l==5, l==6 and l==7 */
						if(pRG->noct > 4){
							swap = pRG->imu[ifr][4][n][k][je+j][i];
							pRG->imu[ifr][4][n][k][je+j][i] = pRG->imu[ifr][5][n][k][je+j][i];
							pRG->imu[ifr][5][n][k][je+j][i] = swap;
							
							swap = pRG->imu[ifr][6][n][k][je+j][i];
							pRG->imu[ifr][6][n][k][je+j][i] = pRG->imu[ifr][7][n][k][je+j][i];
							pRG->imu[ifr][7][n][k][je+j][i] = swap;
							
						}
						
					}/* end i */
				}/* end J */
			} /* End k */
		}/* end nang */ 
	}/* end ifr */
	
	
	
	return;
}





/*----------------------------------------------------------------------------*/
/* Rotate PERIODIC boundary conditions (cont), Inner x2 boundary (rbc_ix2=8) */
/* Angles for specific intensities rotates for 180 degrees  */

static void Rotate180_ix2_fullrad(RadGridS *pRG)
{
	int is = pRG->is, ie = pRG->ie;
	int js = pRG->js;
	int ks = pRG->ks, ke = pRG->ke;
	int nang = pRG->nang;
	int nf = pRG->nf;
	int i, j, k, n, ifr;
	Real swap;
	
	
	for(ifr=0; ifr<nf; ifr++) {		
			for(n=0; n<nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=1; j<=Radghost; j++){
						for(i=is-Radghost; i<=ie+Radghost; i++){
							/* Swap l = 0 && l =2, l=1 and l=3 */
							swap = pRG->imu[ifr][0][n][k][js-j][i];
							pRG->imu[ifr][0][n][k][js-j][i] = pRG->imu[ifr][2][n][k][js-j][i];
							pRG->imu[ifr][2][n][k][js-j][i] = swap;

							swap = pRG->imu[ifr][1][n][k][js-j][i];
                                                        pRG->imu[ifr][1][n][k][js-j][i] = pRG->imu[ifr][3][n][k][js-j][i];
                                                        pRG->imu[ifr][3][n][k][js-j][i] = swap;
							
							/* If for 3D */ 
							/* swap l ==4 && l==6, l==5 and l==7 */
							if(pRG->noct > 4){
								swap = pRG->imu[ifr][4][n][k][js-j][i];
								pRG->imu[ifr][4][n][k][js-j][i] = pRG->imu[ifr][6][n][k][js-j][i];
								pRG->imu[ifr][6][n][k][js-j][i] = swap;

								swap = pRG->imu[ifr][5][n][k][js-j][i];
                                                                pRG->imu[ifr][5][n][k][js-j][i] = pRG->imu[ifr][7][n][k][js-j][i];
                                                                pRG->imu[ifr][7][n][k][js-j][i] = swap;
							}
							
						}/* end i */
					}/* end J */
				} /* End k */
			}/* end nang */ 
	}/* end ifr */
	
	return;
}

/*----------------------------------------------------------------------------*/
/* Rotate PERIODIC boundary conditions (cont), Outer x2 boundary (rbc_ox2=8) */

static void Rotate180_ox2_fullrad(RadGridS *pRG)
{
	
	int is = pRG->is, ie = pRG->ie;
	int je = pRG->je;
	int ks = pRG->ks, ke = pRG->ke;
	int nang = pRG->nang;
	int nf = pRG->nf;
	int i, j, k, n, ifr;
	Real swap;
	
	
	for(ifr=0; ifr<nf; ifr++) {		
		for(n=0; n<nang; n++){
			for(k=ks; k<=ke; k++){
				for(j=1; j<=Radghost; j++){
					for(i=is-Radghost; i<=ie+Radghost; i++){
						/* Swap l = 0 && l =2 , l==1 and l ==3*/
						swap = pRG->imu[ifr][0][n][k][je+j][i];
						pRG->imu[ifr][0][n][k][je+j][i] = pRG->imu[ifr][2][n][k][je+j][i];
						pRG->imu[ifr][2][n][k][je+j][i] = swap;

						swap = pRG->imu[ifr][1][n][k][je+j][i];
                                                pRG->imu[ifr][1][n][k][je+j][i] = pRG->imu[ifr][3][n][k][je+j][i];
                                                pRG->imu[ifr][3][n][k][je+j][i] = swap;

						/* swap l ==4 && l==6, l==5 and l==7 */
						if(pRG->noct > 4){
							swap = pRG->imu[ifr][4][n][k][je+j][i];
							pRG->imu[ifr][4][n][k][je+j][i] = pRG->imu[ifr][6][n][k][je+j][i];
							pRG->imu[ifr][6][n][k][je+j][i] = swap;

							swap = pRG->imu[ifr][5][n][k][je+j][i];
                                                        pRG->imu[ifr][5][n][k][je+j][i] = pRG->imu[ifr][7][n][k][je+j][i];
                                                        pRG->imu[ifr][7][n][k][je+j][i] = swap;

						}
						
					}/* end i */
				}/* end J */
			} /* End k */
		}/* end nang */ 
	}/* end ifr */

	
	
	return;
}



/*---------------------------------------------------------------------------*/
/* The vacuum boundary conditions means outgoing specifiy intensity is continuous,*/
/* But incoming specific intensity is zero */

static void vacuum_ix1_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;
#ifdef CYLINDRICAL
  Real  miux, miuy;
#endif

  /* set the angle independent values */

for(ifr=0; ifr<nf; ifr++) {


 for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
#ifdef CYLINDRICAL
			miux = pRG->Rphimu[l][n][k][j][is-i][0]; 
			miuy = pRG->Rphimu[l][n][k][j][is-i][1];

			if(miux > 0.0){
				pRG->imu[ifr][l][n][k][j][is-i] = 0.0;  
			}else {
				pRG->imu[ifr][l][n][k][j][is-i] = pRG->imu[ifr][l][n][k][j][is];
			}
				
#else
			  
			if((l == 1) || (l == 3) || (l == 5) || (l == 7)){
				pRG->imu[ifr][l][n][k][j][is-i] = pRG->imu[ifr][l][n][k][j][is];
			}	
			else{
				pRG->imu[ifr][l][n][k][j][is-i] = 0.0;
			}	
#endif			  
		/*	This needs to be re-calculated based on local specific intensity */
		/*	pRG->heatcool[ifr][l][n][k][j][is-i] = pRG->heatcool[ifr][l][n][k][j][is];
		*/

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */

/*
  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][is-i] = pRG->R[ifr][k][j][is];
	}
     }
  }
* This needs to be re-calculated 
*/
  
   
   }/* end ifr */
	  


  return;
}


static void vacuum_ox1_fullrad(GridS *pG, RadGridS *pRG)
{

  int ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;
#ifdef CYLINDRICAL
	Real miux, miuy;
#endif	

  /* set the angle independent values */

for(ifr=0; ifr<nf; ifr++) {


 for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
        for(k=ks; k<=ke; k++){
     	    for(j=js; j<=je; j++){
	      for(i=1; i<=Radghost; i++){
#ifdef CYLINDRICAL
			  miux = pRG->Rphimu[l][n][k][j][ie+i][0]; 
			  miuy = pRG->Rphimu[l][n][k][j][ie+i][1];
			  if(miux < 0.0){
				  pRG->imu[ifr][l][n][k][j][ie+i] = 0.0;  
			  }else {
				  pRG->imu[ifr][l][n][k][j][ie+i] = pRG->imu[ifr][l][n][k][j][ie];
			  }
			  
#else			  
			if((l == 0) || (l == 2) || (l == 4) || (l == 6)){
				pRG->imu[ifr][l][n][k][j][ie+i] = pRG->imu[ifr][l][n][k][j][ie];
			}	
			else{
				pRG->imu[ifr][l][n][k][j][ie+i] = 0.0;
			}
			  
#endif
		/*	This needs to be re-calculated based on local specific intensity */
		/*	pRG->heatcool[ifr][l][n][k][j][is-i] = pRG->heatcool[ifr][l][n][k][j][is];
		*/

	      }/* end i */
       	    }/* end J */
         } /* End k */
       }  /* end nang */
     }/* end noctant */

/*
  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){
	for(i=1; i<=Radghost; i++){
	    pRG->R[ifr][k][j][is-i] = pRG->R[ifr][k][j][is];
	}
     }
  }
* This needs to be re-calculated 
*/
  
   
   }/* end ifr */
	  


  return;
}


static void vacuum_ix2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
   /* for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			
			pRG->R[ifr][k][js-j][i] = pRG->R[ifr][k][js][i];
     	   }
       }
    } 
*/
   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			if((l == 2) || (l == 3) || (l == 6) || (l == 7)){
				pRG->imu[ifr][l][n][k][js-j][i] = pRG->imu[ifr][l][n][k][js][i];
			}
			else{
				pRG->imu[ifr][l][n][k][js-j][i] = 0.0;
			}			
		/*	pRG->heatcool[ifr][l][n][k][js-j][i] = pRG->heatcool[ifr][l][n][k][js][i];
		*/ 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */

  return;
}




static void vacuum_ox2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;
  int je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
   /* for(k=ks; k<=ke; k++){
       for(j=1; j<=Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){	      
			
			pRG->R[ifr][k][js-j][i] = pRG->R[ifr][k][js][i];
     	   }
       }
    } 
*/
   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
       	    for(j=1; j<=Radghost; j++){
	       for(i=is-Radghost; i<=ie+Radghost; i++){	
			if((l == 0) || (l == 1) || (l == 4) || (l == 5)){
				pRG->imu[ifr][l][n][k][je+j][i] = pRG->imu[ifr][l][n][k][je][i];
			}
			else{
				pRG->imu[ifr][l][n][k][je+j][i] = 0.0;
			}			
		/*	pRG->heatcool[ifr][l][n][k][js-j][i] = pRG->heatcool[ifr][l][n][k][js][i];
		*/ 
	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */

  return;
}



static void vacuum_ix3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			if((l == 4) || (l == 5) || (l == 6) || (l == 7)){
				pRG->imu[ifr][l][n][ks-k][j][i] = pRG->imu[ifr][l][n][ks][j][i];
			}
			else{			
				pRG->imu[ifr][l][n][ks-k][j][i] = 0.0;
			}
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 

  return;
}




static void vacuum_ox3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;


for(ifr=0; ifr<nf; ifr++) {
    
 
    for(l=0; l<noct; l++){
	for(n=0; n<nang; n++){
	   for(k=1; k<=Radghost; k++){
       	      for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
			if((l == 0) || (l == 1) || (l == 2) || (l == 3)){
				pRG->imu[ifr][l][n][ke+k][j][i] = pRG->imu[ifr][l][n][ke][j][i];
			}
			else{			
				pRG->imu[ifr][l][n][ke+k][j][i] = 0.0;
			}
		 }/* end i */
       	      }/* end J */
    	   } /* End k */
    }/* end nang */
  }/* end noctant */
}/* end ifr */
	 

  return;
}





/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pG, RadGridS *pRG)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic.
 * Inner x1 boundary  */

static void const_flux_ix1(GridS *pG, RadGridS *pRG)
{

  return;
}

static void const_flux_ox1(GridS *pG, RadGridS *pRG)
{

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Inner x2 boundary  */

static void const_flux_ix2(GridS *pG, RadGridS *pRG)
{
  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Outer x2 boundary  */

static void const_flux_ox2(GridS *pG, RadGridS *pRG)
{

  return;
}
/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic. 
 * Inner x3 boundary  */

static void const_flux_ix3(GridS *pG, RadGridS *pRG)
{


  return;
}

/*----------------------------------------------------------------------------*/
/* Enforces a time independent constant flux on the boundary assuming the
 * the incoming radiation field is isotropic.
 * Outer x3 boundary  */

static void const_flux_ox3(GridS *pG, RadGridS *pRG)
{
  return;

}

/*----------------------------------------------------------------------------*/
/* Time independent incident radiation boundary condition.  Nothing is done 
 * here, as the incident boundary radiaion is specified at initialization and
 * unchanged during computation. This means that any contribution from back
 * scattering of outgoing radiation is ignored.  */



#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs */

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=js; j<=je; j++){
	   for(i=is; i<=is+(Radghost-1); i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
       	      for(j=js; j<=je; j++){
	   	 for(i=is; i<=is+(Radghost-1); i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1_fullrad(GridS *pG, RadGridS *pRG)
{


  int ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);


  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=js; j<=je; j++){
	   for(i=ie-(Radghost-1); i<=ie; i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
       	      for(j=js; j<=je; j++){
	   	 for(i=ie-(Radghost-1); i<=ie; i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);



  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=js; j<=js+(Radghost-1); j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
       	     for(j=js; j<=js+(Radghost-1); j++){
	        for(i=is-Radghost; i<=ie+Radghost; i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   
  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);



  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=je-(Radghost-1); j<=je; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
       	      for(j=je-(Radghost-1); j<=je; j++){
	   	 for(i=is-Radghost; i<=ie+Radghost; i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3_fullrad(GridS *pG, RadGridS *pRG)
{


  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;;
  int ks = pRG->ks;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);



  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ks+(Radghost-1); k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ks+(Radghost-1); k++){
              for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3_fullrad(GridS *pG, RadGridS *pRG)
{

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;;
  int ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);




  for(ifr=0; ifr<nf; ifr++) {
    for(k=ke-(Radghost-1); k<=ke; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
			/* set the angle independent values */
			*(pSnd++) = pRG->R[ifr][k][j][i].J;

			/* flux */
			for(m=0; m<3; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].H[m];

			/* rad pressure */
			for(m=0; m<6; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].K[m];

			/* Opacity */
			for(m=0; m<4; m++)
				*(pSnd++) = pRG->R[ifr][k][j][i].Sigma[m];

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ke-(Radghost-1); k<=ke; k++){
              for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){
		    *(pSnd++) = pRG->imu[ifr][l][n][k][j][i];
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1_fullrad(GridS *pG, RadGridS *pRG)
{


  int is = pRG->is;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);



  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=js; j<=je; j++){
	   for(i=is-Radghost; i<=is-1; i++){
	     
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
       	      for(j=js; j<=je; j++){
	         for(i=is-Radghost; i<=is-1; i++){
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1_fullrad(GridS *pG, RadGridS *pRG)
{

  int ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);



  for(ifr=0; ifr<nf; ifr++) {
    for(k=ks; k<=ke; k++){
       for(j=js; j<=je; j++){
	   for(i=ie+1; i<=ie+Radghost; i++){
	     
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
              for(j=js; j<=je; j++){
	         for(i=ie+1; i<=ie+Radghost; i++){
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   
  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;;
  int js = pRG->js;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);



  for(ifr=0; ifr<nf; ifr++) {
     for(k=ks; k<=ke; k++){
        for(j=js-Radghost; j<=js-1; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
	     
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ks; k<=ke; k++){
	      for(j=js-Radghost; j<=js-1; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){	       
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;;
  int je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);






  for(ifr=0; ifr<nf; ifr++) {
     for(k=ks; k<=ke; k++){
       for(j=je+1; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
	     
	     
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	  for(k=ks; k<=ke; k++){
             for(j=je+1; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){       
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	   


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);


  for(ifr=0; ifr<nf; ifr++) {
     for(k=ks-Radghost; k<=ks-1; k++){
       for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	  for(k=ks-Radghost; k<=ks-1; k++){
             for(j=js-Radghost; j<=je+Radghost; j++){
	        for(i=is-Radghost; i<=ie+Radghost; i++){      
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	

  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3_fullrad(GridS *pG, RadGridS *pRG)
{
  int is = pRG->is, ie = pRG->ie;;
  int js = pRG->js, je = pRG->je;
  int ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, m, n, ifr;
  
  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);




  for(ifr=0; ifr<nf; ifr++) {
     for(k=ke+1; k<=ke+Radghost; k++){
         for(j=js-Radghost; j<=je+Radghost; j++){
	   for(i=is-Radghost; i<=ie+Radghost; i++){
			/* set the angle independent values */
			pRG->R[ifr][k][j][i].J = *(pRcv++);

			/* flux */
			for(m=0; m<3; m++)
				pRG->R[ifr][k][j][i].H[m] = *(pRcv++);

			/* rad pressure */
			for(m=0; m<6; m++)
				pRG->R[ifr][k][j][i].K[m] = *(pRcv++);

			/* Opacity */
			for(m=0; m<4; m++)
				pRG->R[ifr][k][j][i].Sigma[m] = *(pRcv++);

	   }/* end i */
        }/* end J */
     } /* End k */	

			
    for(l=0; l<noct; l++){
       for(n=0; n<nang; n++){
	   for(k=ke+1; k<=ke+Radghost; k++){
              for(j=js-Radghost; j<=je+Radghost; j++){
	         for(i=is-Radghost; i<=ie+Radghost; i++){      
		    pRG->imu[ifr][l][n][k][j][i] = *(pRcv++);
		/*	*(pSnd++) = pRG->Scat[k][j][i][ifr][l][n];
			*(pSnd++) = pRG->heatcool[k][j][i][ifr][l][n];
		*/	

	      }/* end i */
       	   }/* end J */
    	} /* End k */		

       }/* end nang */
    }/* end noctant */
}/* end ifr */
	




  return;
}

#endif /* MPI_PARALLEL */

#endif /* FULL_RADIATION_TRANSFER */
