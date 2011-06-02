#include "copyright.h"
/*==============================================================================
 * FILE: bvals_radMHD.c
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) on each edge
 *   of a Grid for the MHD variables.  Each edge of a Grid represents either:
 *    (1) a physical boundary at the edge of the Mesh; in which case BCs are
 *        specified by an integer flag input by user (or by user-defined BC
 *        function in the problem file)
 *    (2) the boundary between Grids resulting from decomposition of a larger
 *        Domain using MPI; in which case BCs are handled by MPI calls
 *    (3) an internal boundary between fine/coarse grid levels in a nested Mesh;
 *        in which case the BCs require a prolongation operator from the parent
 *   This file contains functions that can handle the first two cases.  Case (3)
 *   is handled in the Prolongate function called from main loop.
 *   The naming convention of the integer flags for BCs is:
 *       bc_ix1 = Boundary Condition for Inner x1 (at i=is)
 *       bc_ox1 = Boundary Condition for Outer x1 (at i=ie)
 *   similarly for bc_ix2; bc_ox2; bc_ix3; bc_ox3

 *
 *
 * For case (1) -- PHYSICAL BOUNDARIES
 *   The values of the integer flags (bc_ix1, etc.) are:
 *     1 = reflecting; 2 = outflow; 4 = periodic; 5 = conductor
 *   For flow-in bondaries (ghost zones set to pre-determined values), pointers
 *   to user-defined functions in the problem file are used. 
 *
 * For case (2) -- MPI BOUNDARIES
 *   We do the parallel synchronization by having every grid:
 *     1) Post non-blocking receives for data from both L and R Grids
 *     2) Pack and send data to the Grids on both L and R
 *     3) Check for receives and unpack data in order of first to finish
 *   If the Grid is at the edge of the Domain, we set BCs as in case (1) or (3).
 *
 * For case (3) -- INTERNAL GRID LEVEL BOUNDARIES
 *   This step is complicated and must be handled separately, in the function
 *   Prolongate() called from the main loop.  In the algorithm below, nothing is
 *   done for this case; the BCs are left to be set later.
 *
 * With SELF-GRAVITY: BCs for Phi are set independently of the MHD variables
 *   in a separate function bvals_grav(). 
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   bvals_radMHD()      - calls appropriate functions to set ghost cells
 *   bvals_radMHD_init() - sets function pointers used by bvals_mhd()
 *   bvals_rad_fun()  - enrolls a pointer to a user-defined BC function for radiation quantities
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
static int Nrad = 12;
#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - reflecting BCs at boundary ???
 *   outflow_???()  - outflow BCs at boundary ???
 *   periodic_???() - periodic BCs at boundary ???
 *   conduct_???()  - conducting BCs at boundary ???
 *   pack_???()     - pack data for MPI non-blocking send at ??? boundary
 *   unpack_???()   - unpack data for MPI non-blocking receive at ??? boundary
 * NOTICE private functions are the same as functions in bvals_mhd 
 *============================================================================*/

static void reflect_ix1(GridS *pG);
static void reflect_ox1(GridS *pG);
static void reflect_ix2(GridS *pG);
static void reflect_ox2(GridS *pG);
static void reflect_ix3(GridS *pG);
static void reflect_ox3(GridS *pG);

static void outflow_ix1(GridS *pG);
static void outflow_ox1(GridS *pG);
static void outflow_ix2(GridS *pG);
static void outflow_ox2(GridS *pG);
static void outflow_ix3(GridS *pG);
static void outflow_ox3(GridS *pG);

static void periodic_ix1(GridS *pG);
static void periodic_ox1(GridS *pG);
static void periodic_ix2(GridS *pG);
static void periodic_ox2(GridS *pG);
static void periodic_ix3(GridS *pG);
static void periodic_ox3(GridS *pG);

static void conduct_ix1(GridS *pG);
static void conduct_ox1(GridS *pG);
static void conduct_ix2(GridS *pG);
static void conduct_ox2(GridS *pG);
static void conduct_ix3(GridS *pG);
static void conduct_ox3(GridS *pG);

static void ProlongateLater(GridS *pG);

#ifdef MPI_PARALLEL
static void pack_ix1(GridS *pG);
static void pack_ox1(GridS *pG);
static void pack_ix2(GridS *pG);
static void pack_ox2(GridS *pG);
static void pack_ix3(GridS *pG);
static void pack_ox3(GridS *pG);

static void unpack_ix1(GridS *pG);
static void unpack_ox1(GridS *pG);
static void unpack_ix2(GridS *pG);
static void unpack_ox2(GridS *pG);
static void unpack_ix3(GridS *pG);
static void unpack_ox3(GridS *pG);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* bvals_radMHD: calls appropriate functions to set ghost zones.  The function
 *   pointers (*(pD->???_BCFun)) are set by bvals_init() to be either a
 *   user-defined function, or one of the functions corresponding to reflecting,
 *   periodic, or outflow.  If the left- or right-Grid ID numbers are >= 1
 *   (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 * This function is form domain now. Should call for every domain in mesh 
 * Only 12 radiation related variables are set here. Other variables in ghost
 * zones are set in function bvals_MHD 
 */

void bvals_radMHD(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN,BCFlag;
#endif
#ifdef MPI_PARALLEL
  int cnt, cnt2, cnt3, ierr, mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL
    cnt = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(Nrad);


/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox1(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pGrid);
      if (mIndex == 1) unpack_ox1(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pGrid);
      if (mIndex == 1) unpack_ox1(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx1_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->rad_ix1_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx1_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->rad_ox1_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*(pD->rad_ix1_BCFun))(pGrid);
      (*(pD->rad_ox1_BCFun))(pGrid);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt = (pGrid->Nx[0] + 2*nghost)*nghost*(pGrid->Nx[2])*(Nrad);


/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox2(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pGrid);
      if (mIndex == 1) unpack_ox2(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pGrid);
      if (mIndex == 1) unpack_ox2(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx2_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->rad_ix2_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx2_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->rad_ox2_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*(pD->rad_ix2_BCFun))(pGrid);
      (*(pD->rad_ox2_BCFun))(pGrid);
    } 

/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */
/* This is done after periodic boundary condition has been applied */
#ifdef SHEARING_BOX
    BCFlag = par_geti_def("domain1","bc_ix1",0);
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0 && BCFlag == 4) {
      ShearingSheet_radMHD_ix1(pD);
    }
    BCFlag = par_geti_def("domain1","bc_ox1",0);
    if (myL == ((pD->NGrid[0])-1) && BCFlag == 4) {
      ShearingSheet_radMHD_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL
    cnt = (pGrid->Nx[0] + 2*nghost)*(pGrid->Nx[1] + 2*nghost)*nghost*(Nrad);


/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3(pGrid);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      pack_ox3(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pGrid);
      if (mIndex == 1) unpack_ox3(pGrid);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pGrid);
      if (mIndex == 1) unpack_ox3(pGrid);

    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,RtoL_tag,
        pD->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3(pGrid); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pGrid->rx3_id,LtoR_tag,
        pD->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(pD->rad_ix3_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3(pGrid);

    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,LtoR_tag,
        pD->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3(pGrid); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pGrid->lx3_id,RtoL_tag,
        pD->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(pD->rad_ox3_BCFun))(pGrid);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3(pGrid);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*(pD->rad_ix3_BCFun))(pGrid);
      (*(pD->rad_ox3_BCFun))(pGrid);
    } 

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* bvals_radMHD_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void bvals_radMHD_init(MeshS *pM)
{
  GridS *pG;
  DomainS *pD;
  int i,nl,nd,irefine;
#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active Grids on this proc */

  for (nl=0; nl<(pM->NLevels); nl++){
  for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
  if (pM->Domain[nl][nd].Grid != NULL) {
    pD = (DomainS*)&(pM->Domain[nl][nd]);  /* ptr to Domain */
    pG = pM->Domain[nl][nd].Grid;          /* ptr to Grid */
    irefine = 1;
    for (i=1;i<=nl;i++) irefine *= 2;   /* C pow fn only takes doubles !! */
#ifdef MPI_PARALLEL
/* get (l,m,n) coordinates of Grid being updated on this processor */
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
#endif /* MPI_PARALLEL */

/* Set function pointers for physical boundaries in x1-direction -------------*/

    if(pG->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

      if(pD->rad_ix1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {      
          pD->rad_ix1_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
/* BCFlag_??? is the same for both radiation and hydro variables      */
        } else {
          if(((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) && 
               pM->BCFlag_ix1 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix1b but not ox1b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {                    
            switch(pM->BCFlag_ix1){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ix1_BCFun = reflect_ix1;
            break;

            case 2: /* Outflow */
              pD->rad_ix1_BCFun = outflow_ix1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ix1_BCFun = periodic_ix1;
#ifdef MPI_PARALLEL
              if(pG->lx1_id < 0 && pD->NGrid[0] > 1){
                pG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ix1_BCFun = conduct_ix1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix1=%d unknown\n",pM->BCFlag_ix1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->rad_ox1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          pD->rad_ox1_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[0] != 0) && (pM->BCFlag_ox1 == 4)) {      
            ath_error("[bvals_init]:level=%d Domain touching ox1b but not ix1b and periodic BC not allowed\n",nl); 


/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox1){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ox1_BCFun = reflect_ox1;
            break;

            case 2: /* Outflow */
              pD->rad_ox1_BCFun = outflow_ox1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ox1_BCFun = periodic_ox1;
#ifdef MPI_PARALLEL
              if(pG->rx1_id < 0 && pD->NGrid[0] > 1){
                pG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ox1_BCFun = conduct_ox1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox1=%d unknown\n",pM->BCFlag_ox1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x2-direction -------------*/

    if(pG->Nx[1] > 1) {

/*---- ix2 boundary ----------------------------------------------------------*/

      if(pD->rad_ix2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          pD->rad_ix2_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) &&
               pM->BCFlag_ix2 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix2b but not ox2b and periodic BC not allowed\n",nl); 


/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix2){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ix2_BCFun = reflect_ix2;
            break;

            case 2: /* Outflow */
              pD->rad_ix2_BCFun = outflow_ix2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ix2_BCFun = periodic_ix2;
#ifdef MPI_PARALLEL
              if(pG->lx2_id < 0 && pD->NGrid[1] > 1){
                pG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;
  
            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ix2_BCFun = conduct_ix2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix2=%d unknown\n",pM->BCFlag_ix2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(pD->rad_ox2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          pD->rad_ox2_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[1] != 0) && (pM->BCFlag_ox2 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox2b but not ix2b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox2){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ox2_BCFun = reflect_ox2;
            break;

            case 2: /* Outflow */
              pD->rad_ox2_BCFun = outflow_ox2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ox2_BCFun = periodic_ox2;
#ifdef MPI_PARALLEL
              if(pG->rx2_id < 0 && pD->NGrid[1] > 1){
                pG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ox2_BCFun = conduct_ox2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox2=%d unknown\n",pM->BCFlag_ox2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }
    }

/* Set function pointers for physical boundaries in x3-direction -------------*/

    if(pG->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

      if(pD->rad_ix3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          pD->rad_ix3_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) &&
               pM->BCFlag_ix3 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix3b but not ox3b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix3){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ix3_BCFun = reflect_ix3;
            break;

            case 2: /* Outflow */
              pD->rad_ix3_BCFun = outflow_ix3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ix3_BCFun = periodic_ix3;
#ifdef MPI_PARALLEL
              if(pG->lx3_id < 0 && pD->NGrid[2] > 1){
                pG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ix3_BCFun = conduct_ix3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix3=%d unknown\n",pM->BCFlag_ix3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

      if(pD->rad_ox3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          pD->rad_ox3_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[2] != 0) && (pM->BCFlag_ox3 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox3b but not ix3b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox3){

            case 1: /* Reflecting, B_normal=0 */
              pD->rad_ox3_BCFun = reflect_ox3;
            break;

            case 2: /* Outflow */
              pD->rad_ox3_BCFun = outflow_ox3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->rad_ox3_BCFun = periodic_ox3;
#ifdef MPI_PARALLEL
              if(pG->rx3_id < 0 && pD->NGrid[2] > 1){
                pG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->rad_ox3_BCFun = conduct_ox3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox3=%d unknown\n",pM->BCFlag_ox3);
              exit(EXIT_FAILURE);
            }
          }
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
 /* Here we only need to send radiation variables *
  * NVAR is the number of hydro variables *
  * total is Er, Fr???, Sigma_??, Edd_?? 12 variables 
  */

  size *= nghost*Nrad;


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
/* bvals_rad_fun:  sets function ptrs for user-defined BCs for radiation quantities
 */

void bvals_rad_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    pD->rad_ix1_BCFun = prob_bc;
    break;
  case right_x1:
    pD->rad_ox1_BCFun = prob_bc;
    break;
  case left_x2:
    pD->rad_ix2_BCFun = prob_bc;
    break;
  case right_x2:
    pD->rad_ox2_BCFun = prob_bc;
    break;
  case left_x3:
    pD->rad_ix3_BCFun = prob_bc;
    break;
  case right_x3:
    pD->rad_ox3_BCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   outflow_???
 *   periodic_???
 *   conduct_???
 *   pack_???
 *   unpack_???
 */

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (bc_ix1=1) */

static void reflect_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
		pGrid->U[k][j][is-i].Er  	=  pGrid->U[k][j][is+(i-1)].Er;
		pGrid->U[k][j][is-i].Fr1 	= -pGrid->U[k][j][is+(i-1)].Fr1; /* reflect 1-flux. */
		pGrid->U[k][j][is-i].Fr2  	=  pGrid->U[k][j][is+(i-1)].Fr2;
		pGrid->U[k][j][is-i].Fr3  	=  pGrid->U[k][j][is+(i-1)].Fr3;
		pGrid->U[k][j][is-i].Edd_11  	=  pGrid->U[k][j][is+(i-1)].Edd_11;
		pGrid->U[k][j][is-i].Edd_21  	=  pGrid->U[k][j][is+(i-1)].Edd_21;
		pGrid->U[k][j][is-i].Edd_22  	=  pGrid->U[k][j][is+(i-1)].Edd_22;
		pGrid->U[k][j][is-i].Edd_31  	=  pGrid->U[k][j][is+(i-1)].Edd_31;
		pGrid->U[k][j][is-i].Edd_32  	=  pGrid->U[k][j][is+(i-1)].Edd_32;
		pGrid->U[k][j][is-i].Edd_33  	=  pGrid->U[k][j][is+(i-1)].Edd_33;
		pGrid->U[k][j][is-i].Sigma_t  	=  pGrid->U[k][j][is+(i-1)].Sigma_t;
		pGrid->U[k][j][is-i].Sigma_a  	=  pGrid->U[k][j][is+(i-1)].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1) */

static void reflect_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      		pGrid->U[k][j][ie+i].Er    	=  pGrid->U[k][j][ie-(i-1)].Er;
		pGrid->U[k][j][ie+i].Fr1 	= -pGrid->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */
		pGrid->U[k][j][ie+i].Fr2    	=  pGrid->U[k][j][ie-(i-1)].Fr2;
		pGrid->U[k][j][ie+i].Fr3    	=  pGrid->U[k][j][ie-(i-1)].Fr3;
		pGrid->U[k][j][ie+i].Edd_11  	=  pGrid->U[k][j][ie-(i-1)].Edd_11;
		pGrid->U[k][j][ie+i].Edd_21  	=  pGrid->U[k][j][ie-(i-1)].Edd_21;
		pGrid->U[k][j][ie+i].Edd_22  	=  pGrid->U[k][j][ie-(i-1)].Edd_22;
		pGrid->U[k][j][ie+i].Edd_31  	=  pGrid->U[k][j][ie-(i-1)].Edd_31;
		pGrid->U[k][j][ie+i].Edd_32  	=  pGrid->U[k][j][ie-(i-1)].Edd_32;
		pGrid->U[k][j][ie+i].Edd_33  	=  pGrid->U[k][j][ie-(i-1)].Edd_33;				
		pGrid->U[k][j][ie+i].Sigma_t    =  pGrid->U[k][j][ie-(i-1)].Sigma_t;
		pGrid->U[k][j][ie+i].Sigma_a    =  pGrid->U[k][j][ie-(i-1)].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1) */

static void reflect_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[k][js-j][i].Er   	=  pGrid->U[k][js+(j-1)][i].Er;
		pGrid->U[k][js-j][i].Fr1  	=  pGrid->U[k][js+(j-1)][i].Fr1; 
		pGrid->U[k][js-j][i].Fr2  	= -pGrid->U[k][js+(j-1)][i].Fr2; /* reflect 2-flux. */
		pGrid->U[k][js-j][i].Fr3  	=  pGrid->U[k][js+(j-1)][i].Fr3; 
		pGrid->U[k][js-j][i].Edd_11  	=  pGrid->U[k][js+(j-1)][i].Edd_11;
		pGrid->U[k][js-j][i].Edd_21  	=  pGrid->U[k][js+(j-1)][i].Edd_21;
		pGrid->U[k][js-j][i].Edd_22  	=  pGrid->U[k][js+(j-1)][i].Edd_22;
		pGrid->U[k][js-j][i].Edd_31  	=  pGrid->U[k][js+(j-1)][i].Edd_31;
		pGrid->U[k][js-j][i].Edd_32  	=  pGrid->U[k][js+(j-1)][i].Edd_32;
		pGrid->U[k][js-j][i].Edd_33  	=  pGrid->U[k][js+(j-1)][i].Edd_33;
		pGrid->U[k][js-j][i].Sigma_t  	=  pGrid->U[k][js+(j-1)][i].Sigma_t;
		pGrid->U[k][js-j][i].Sigma_a  	=  pGrid->U[k][js+(j-1)][i].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1) */

static void reflect_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[k][je+j][i].Er   	=  pGrid->U[k][je-(j-1)][i].Er;
		pGrid->U[k][je+j][i].Fr1  	=  pGrid->U[k][je-(j-1)][i].Fr1; 
		pGrid->U[k][je+j][i].Fr2  	= -pGrid->U[k][je-(j-1)][i].Fr2; /* reflect 2-flux. */
		pGrid->U[k][je+j][i].Fr3  	=  pGrid->U[k][je-(j-1)][i].Fr3; 
		pGrid->U[k][je+j][i].Edd_11  	=  pGrid->U[k][je-(j-1)][i].Edd_11;
		pGrid->U[k][je+j][i].Edd_21  	=  pGrid->U[k][je-(j-1)][i].Edd_21;
		pGrid->U[k][je+j][i].Edd_22  	=  pGrid->U[k][je-(j-1)][i].Edd_22;
		pGrid->U[k][je+j][i].Edd_31  	=  pGrid->U[k][je-(j-1)][i].Edd_31;
		pGrid->U[k][je+j][i].Edd_32  	=  pGrid->U[k][je-(j-1)][i].Edd_32;
		pGrid->U[k][je+j][i].Edd_33  	=  pGrid->U[k][je-(j-1)][i].Edd_33;
		pGrid->U[k][je+j][i].Sigma_t  	=  pGrid->U[k][je-(j-1)][i].Sigma_t;
		pGrid->U[k][je+j][i].Sigma_a  	=  pGrid->U[k][je-(j-1)][i].Sigma_a;

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1) */

static void reflect_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[ks-k][j][i].Er   	=  pGrid->U[ks+(k-1)][j][i].Er;
		pGrid->U[ks-k][j][i].Fr1  	=  pGrid->U[ks+(k-1)][j][i].Fr1; 
		pGrid->U[ks-k][j][i].Fr2  	=  pGrid->U[ks+(k-1)][j][i].Fr2; 
		pGrid->U[ks-k][j][i].Fr3  	= -pGrid->U[ks+(k-1)][j][i].Fr3; /* reflect 3-flux. */
		pGrid->U[ks-k][j][i].Edd_11  	=  pGrid->U[ks+(k-1)][j][i].Edd_11;
		pGrid->U[ks-k][j][i].Edd_21  	=  pGrid->U[ks+(k-1)][j][i].Edd_21;
		pGrid->U[ks-k][j][i].Edd_22  	=  pGrid->U[ks+(k-1)][j][i].Edd_22;
		pGrid->U[ks-k][j][i].Edd_31  	=  pGrid->U[ks+(k-1)][j][i].Edd_31;
		pGrid->U[ks-k][j][i].Edd_32  	=  pGrid->U[ks+(k-1)][j][i].Edd_32;
		pGrid->U[ks-k][j][i].Edd_33  	=  pGrid->U[ks+(k-1)][j][i].Edd_33;
		pGrid->U[ks-k][j][i].Sigma_t  	=  pGrid->U[ks+(k-1)][j][i].Sigma_t;
		pGrid->U[ks-k][j][i].Sigma_a  	=  pGrid->U[ks+(k-1)][j][i].Sigma_a;

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1) */

static void reflect_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
    		pGrid->U[ke+k][j][i].Er   	=  pGrid->U[ke-(k-1)][j][i].Er;	
		pGrid->U[ke+k][j][i].Fr1  	=  pGrid->U[ke-(k-1)][j][i].Fr1; 
		pGrid->U[ke+k][j][i].Fr2  	=  pGrid->U[ke-(k-1)][j][i].Fr2; 
		pGrid->U[ke+k][j][i].Fr3  	= -pGrid->U[ke-(k-1)][j][i].Fr3; /* reflect 3-flux. */
		pGrid->U[ke+k][j][i].Edd_11  	=  pGrid->U[ke-(k-1)][j][i].Edd_11;
		pGrid->U[ke+k][j][i].Edd_21  	=  pGrid->U[ke-(k-1)][j][i].Edd_21;
		pGrid->U[ke+k][j][i].Edd_22  	=  pGrid->U[ke-(k-1)][j][i].Edd_22;
		pGrid->U[ke+k][j][i].Edd_31  	=  pGrid->U[ke-(k-1)][j][i].Edd_31;
		pGrid->U[ke+k][j][i].Edd_32  	=  pGrid->U[ke-(k-1)][j][i].Edd_32;
		pGrid->U[ke+k][j][i].Edd_33  	=  pGrid->U[ke-(k-1)][j][i].Edd_33;
		pGrid->U[ke+k][j][i].Sigma_t  	=  pGrid->U[ke-(k-1)][j][i].Sigma_t;
		pGrid->U[ke+k][j][i].Sigma_a  	=  pGrid->U[ke-(k-1)][j][i].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        	pGrid->U[k][j][is-i].Er  	= pGrid->U[k][j][is].Er;
		pGrid->U[k][j][is-i].Fr1 	= pGrid->U[k][j][is].Fr1;
		pGrid->U[k][j][is-i].Fr2  	= pGrid->U[k][j][is].Fr2;
		pGrid->U[k][j][is-i].Fr3  	= pGrid->U[k][j][is].Fr3;
		pGrid->U[k][j][is-i].Edd_11  	= pGrid->U[k][j][is].Edd_11;
		pGrid->U[k][j][is-i].Edd_21  	= pGrid->U[k][j][is].Edd_21;
		pGrid->U[k][j][is-i].Edd_22  	= pGrid->U[k][j][is].Edd_22;
		pGrid->U[k][j][is-i].Edd_31  	= pGrid->U[k][j][is].Edd_31;
		pGrid->U[k][j][is-i].Edd_32  	= pGrid->U[k][j][is].Edd_32;
		pGrid->U[k][j][is-i].Edd_33  	= pGrid->U[k][j][is].Edd_33;
		pGrid->U[k][j][is-i].Sigma_t  	= pGrid->U[k][j][is].Sigma_t;
		pGrid->U[k][j][is-i].Sigma_a  	= pGrid->U[k][j][is].Sigma_a;
      }
    }
  }



  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        	pGrid->U[k][j][ie+i].Er  	= pGrid->U[k][j][ie].Er;
		pGrid->U[k][j][ie+i].Fr1 	= pGrid->U[k][j][ie].Fr1;
		pGrid->U[k][j][ie+i].Fr2  	= pGrid->U[k][j][ie].Fr2;
		pGrid->U[k][j][ie+i].Fr3  	= pGrid->U[k][j][ie].Fr3;
		pGrid->U[k][j][ie+i].Edd_11  	= pGrid->U[k][j][ie].Edd_11;
		pGrid->U[k][j][ie+i].Edd_21  	= pGrid->U[k][j][ie].Edd_21;
		pGrid->U[k][j][ie+i].Edd_22  	= pGrid->U[k][j][ie].Edd_22;
		pGrid->U[k][j][ie+i].Edd_31  	= pGrid->U[k][j][ie].Edd_31;
		pGrid->U[k][j][ie+i].Edd_32  	= pGrid->U[k][j][ie].Edd_32;
		pGrid->U[k][j][ie+i].Edd_33  	= pGrid->U[k][j][ie].Edd_33;
		pGrid->U[k][j][ie+i].Sigma_t  	= pGrid->U[k][j][ie].Sigma_t;
		pGrid->U[k][j][ie+i].Sigma_a  	= pGrid->U[k][j][ie].Sigma_a;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[k][js-j][i].Er  	=  pGrid->U[k][js][i].Er;
		pGrid->U[k][js-j][i].Fr1 	=  pGrid->U[k][js][i].Fr1;
		pGrid->U[k][js-j][i].Fr2  	=  pGrid->U[k][js][i].Fr2;
		pGrid->U[k][js-j][i].Fr3  	=  pGrid->U[k][js][i].Fr3;
		pGrid->U[k][js-j][i].Edd_11  	=  pGrid->U[k][js][i].Edd_11;
		pGrid->U[k][js-j][i].Edd_21  	=  pGrid->U[k][js][i].Edd_21;
		pGrid->U[k][js-j][i].Edd_22  	=  pGrid->U[k][js][i].Edd_22;
		pGrid->U[k][js-j][i].Edd_31  	=  pGrid->U[k][js][i].Edd_31;
		pGrid->U[k][js-j][i].Edd_32  	=  pGrid->U[k][js][i].Edd_32;
		pGrid->U[k][js-j][i].Edd_33 	=  pGrid->U[k][js][i].Edd_33;
		pGrid->U[k][js-j][i].Sigma_t 	=  pGrid->U[k][js][i].Sigma_t;
		pGrid->U[k][js-j][i].Sigma_a 	=  pGrid->U[k][js][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[k][je+j][i].Er  	=  pGrid->U[k][je][i].Er;
		pGrid->U[k][je+j][i].Fr1 	=  pGrid->U[k][je][i].Fr1;
		pGrid->U[k][je+j][i].Fr2  	=  pGrid->U[k][je][i].Fr2;
		pGrid->U[k][je+j][i].Fr3  	=  pGrid->U[k][je][i].Fr3;
		pGrid->U[k][je+j][i].Edd_11  	=  pGrid->U[k][je][i].Edd_11;
		pGrid->U[k][je+j][i].Edd_21  	=  pGrid->U[k][je][i].Edd_21;
		pGrid->U[k][je+j][i].Edd_22  	=  pGrid->U[k][je][i].Edd_22;
		pGrid->U[k][je+j][i].Edd_31  	=  pGrid->U[k][je][i].Edd_31;
		pGrid->U[k][je+j][i].Edd_32  	=  pGrid->U[k][je][i].Edd_32;
		pGrid->U[k][je+j][i].Edd_33 	=  pGrid->U[k][je][i].Edd_33;
		pGrid->U[k][je+j][i].Sigma_t 	=  pGrid->U[k][je][i].Sigma_t;
		pGrid->U[k][je+j][i].Sigma_a 	=  pGrid->U[k][je][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[ks-k][j][i].Er  		=  pGrid->U[ks][j][i].Er;
		pGrid->U[ks-k][j][i].Fr1 		=  pGrid->U[ks][j][i].Fr1;
		pGrid->U[ks-k][j][i].Fr2  		=  pGrid->U[ks][j][i].Fr2;
		pGrid->U[ks-k][j][i].Fr3  		=  pGrid->U[ks][j][i].Fr3;
		pGrid->U[ks-k][j][i].Edd_11  		=  pGrid->U[ks][j][i].Edd_11;
		pGrid->U[ks-k][j][i].Edd_21  		=  pGrid->U[ks][j][i].Edd_21;
		pGrid->U[ks-k][j][i].Edd_22  		=  pGrid->U[ks][j][i].Edd_22;
		pGrid->U[ks-k][j][i].Edd_31  		=  pGrid->U[ks][j][i].Edd_31;
		pGrid->U[ks-k][j][i].Edd_32  		=  pGrid->U[ks][j][i].Edd_32;
		pGrid->U[ks-k][j][i].Edd_33 		=  pGrid->U[ks][j][i].Edd_33;
		pGrid->U[ks-k][j][i].Sigma_t 		=  pGrid->U[ks][j][i].Sigma_t;
		pGrid->U[ks-k][j][i].Sigma_a 		=  pGrid->U[ks][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        	pGrid->U[ke+k][j][i].Er  		=  pGrid->U[ke][j][i].Er;
		pGrid->U[ke+k][j][i].Fr1 		=  pGrid->U[ke][j][i].Fr1;
		pGrid->U[ke+k][j][i].Fr2  		=  pGrid->U[ke][j][i].Fr2;
		pGrid->U[ke+k][j][i].Fr3  		=  pGrid->U[ke][j][i].Fr3;	
		pGrid->U[ke+k][j][i].Edd_11  		=  pGrid->U[ke][j][i].Edd_11;
		pGrid->U[ke+k][j][i].Edd_21  		=  pGrid->U[ke][j][i].Edd_21;
		pGrid->U[ke+k][j][i].Edd_22  		=  pGrid->U[ke][j][i].Edd_22;
		pGrid->U[ke+k][j][i].Edd_31  		=  pGrid->U[ke][j][i].Edd_31;
		pGrid->U[ke+k][j][i].Edd_32  		=  pGrid->U[ke][j][i].Edd_32;
		pGrid->U[ke+k][j][i].Edd_33 		=  pGrid->U[ke][j][i].Edd_33;
		pGrid->U[ke+k][j][i].Sigma_t 		=  pGrid->U[ke][j][i].Sigma_t;
		pGrid->U[ke+k][j][i].Sigma_a 		=  pGrid->U[ke][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
       		pGrid->U[k][j][is-i].Er 	=  pGrid->U[k][j][ie-(i-1)].Er;
		pGrid->U[k][j][is-i].Fr1 	=  pGrid->U[k][j][ie-(i-1)].Fr1;
		pGrid->U[k][j][is-i].Fr2 	=  pGrid->U[k][j][ie-(i-1)].Fr2;
		pGrid->U[k][j][is-i].Fr3 	=  pGrid->U[k][j][ie-(i-1)].Fr3;	
		pGrid->U[k][j][is-i].Edd_11  	=  pGrid->U[k][j][ie-(i-1)].Edd_11;
		pGrid->U[k][j][is-i].Edd_21  	=  pGrid->U[k][j][ie-(i-1)].Edd_21;
		pGrid->U[k][j][is-i].Edd_22  	=  pGrid->U[k][j][ie-(i-1)].Edd_22;
		pGrid->U[k][j][is-i].Edd_31  	=  pGrid->U[k][j][ie-(i-1)].Edd_31;
		pGrid->U[k][j][is-i].Edd_32  	=  pGrid->U[k][j][ie-(i-1)].Edd_32;
		pGrid->U[k][j][is-i].Edd_33  	=  pGrid->U[k][j][ie-(i-1)].Edd_33;
		pGrid->U[k][j][is-i].Sigma_t 	=  pGrid->U[k][j][ie-(i-1)].Sigma_t;
		pGrid->U[k][j][is-i].Sigma_a 	=  pGrid->U[k][j][ie-(i-1)].Sigma_a;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        	pGrid->U[k][j][ie+i].Er 	=  pGrid->U[k][j][is+(i-1)].Er;
		pGrid->U[k][j][ie+i].Fr1 	=  pGrid->U[k][j][is+(i-1)].Fr1;
		pGrid->U[k][j][ie+i].Fr2 	=  pGrid->U[k][j][is+(i-1)].Fr2;
		pGrid->U[k][j][ie+i].Fr3 	=  pGrid->U[k][j][is+(i-1)].Fr3;
		pGrid->U[k][j][ie+i].Edd_11  	=  pGrid->U[k][j][is+(i-1)].Edd_11;
		pGrid->U[k][j][ie+i].Edd_21  	=  pGrid->U[k][j][is+(i-1)].Edd_21;
		pGrid->U[k][j][ie+i].Edd_22  	=  pGrid->U[k][j][is+(i-1)].Edd_22;
		pGrid->U[k][j][ie+i].Edd_31  	=  pGrid->U[k][j][is+(i-1)].Edd_31;
		pGrid->U[k][j][ie+i].Edd_32  	=  pGrid->U[k][j][is+(i-1)].Edd_32;
		pGrid->U[k][j][ie+i].Edd_33  	=  pGrid->U[k][j][is+(i-1)].Edd_33;
		pGrid->U[k][j][ie+i].Sigma_t 	=  pGrid->U[k][j][is+(i-1)].Sigma_t;
		pGrid->U[k][j][ie+i].Sigma_a 	=  pGrid->U[k][j][is+(i-1)].Sigma_a;
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[k][js-j][i].Er 	=  pGrid->U[k][je-(j-1)][i].Er;
		pGrid->U[k][js-j][i].Fr1 	=  pGrid->U[k][je-(j-1)][i].Fr1;
		pGrid->U[k][js-j][i].Fr2 	=  pGrid->U[k][je-(j-1)][i].Fr2;
		pGrid->U[k][js-j][i].Fr3 	=  pGrid->U[k][je-(j-1)][i].Fr3;
		pGrid->U[k][js-j][i].Edd_11  	=  pGrid->U[k][je-(j-1)][i].Edd_11;
		pGrid->U[k][js-j][i].Edd_21  	=  pGrid->U[k][je-(j-1)][i].Edd_21;
		pGrid->U[k][js-j][i].Edd_22  	=  pGrid->U[k][je-(j-1)][i].Edd_22;
		pGrid->U[k][js-j][i].Edd_31  	=  pGrid->U[k][je-(j-1)][i].Edd_31;
		pGrid->U[k][js-j][i].Edd_32  	=  pGrid->U[k][je-(j-1)][i].Edd_32;
		pGrid->U[k][js-j][i].Edd_33  	=  pGrid->U[k][je-(j-1)][i].Edd_33;
		pGrid->U[k][js-j][i].Sigma_t 	=  pGrid->U[k][je-(j-1)][i].Sigma_t;
		pGrid->U[k][js-j][i].Sigma_a 	=  pGrid->U[k][je-(j-1)][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
      		pGrid->U[k][je+j][i].Er 	=  pGrid->U[k][js+(j-1)][i].Er;
		pGrid->U[k][je+j][i].Fr1 	=  pGrid->U[k][js+(j-1)][i].Fr1;
		pGrid->U[k][je+j][i].Fr2 	=  pGrid->U[k][js+(j-1)][i].Fr2;
		pGrid->U[k][je+j][i].Fr3 	=  pGrid->U[k][js+(j-1)][i].Fr3;
		pGrid->U[k][je+j][i].Edd_11  	=  pGrid->U[k][js+(j-1)][i].Edd_11;
		pGrid->U[k][je+j][i].Edd_21  	=  pGrid->U[k][js+(j-1)][i].Edd_21;
		pGrid->U[k][je+j][i].Edd_22  	=  pGrid->U[k][js+(j-1)][i].Edd_22;
		pGrid->U[k][je+j][i].Edd_31  	=  pGrid->U[k][js+(j-1)][i].Edd_31;
		pGrid->U[k][je+j][i].Edd_32  	=  pGrid->U[k][js+(j-1)][i].Edd_32;
		pGrid->U[k][je+j][i].Edd_33  	=  pGrid->U[k][js+(j-1)][i].Edd_33;
		pGrid->U[k][je+j][i].Sigma_t 	=  pGrid->U[k][js+(j-1)][i].Sigma_t;
		pGrid->U[k][je+j][i].Sigma_a 	=  pGrid->U[k][js+(j-1)][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[ks-k][j][i].Er 	=  pGrid->U[ke-(k-1)][j][i].Er;
		pGrid->U[ks-k][j][i].Fr1 	=  pGrid->U[ke-(k-1)][j][i].Fr1;
		pGrid->U[ks-k][j][i].Fr2 	=  pGrid->U[ke-(k-1)][j][i].Fr2;
		pGrid->U[ks-k][j][i].Fr3 	=  pGrid->U[ke-(k-1)][j][i].Fr3;		
		pGrid->U[ks-k][j][i].Edd_11  	=  pGrid->U[ke-(k-1)][j][i].Edd_11;
		pGrid->U[ks-k][j][i].Edd_21  	=  pGrid->U[ke-(k-1)][j][i].Edd_21;
		pGrid->U[ks-k][j][i].Edd_22  	=  pGrid->U[ke-(k-1)][j][i].Edd_22;
		pGrid->U[ks-k][j][i].Edd_31  	=  pGrid->U[ke-(k-1)][j][i].Edd_31;
		pGrid->U[ks-k][j][i].Edd_32  	=  pGrid->U[ke-(k-1)][j][i].Edd_32;
		pGrid->U[ks-k][j][i].Edd_33  	=  pGrid->U[ke-(k-1)][j][i].Edd_33;
		pGrid->U[ks-k][j][i].Sigma_t 	=  pGrid->U[ke-(k-1)][j][i].Sigma_t;
		pGrid->U[ks-k][j][i].Sigma_a 	=  pGrid->U[ke-(k-1)][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[ke+k][j][i].Er 	=  pGrid->U[ks+(k-1)][j][i].Er;
		pGrid->U[ke+k][j][i].Fr1 	=  pGrid->U[ks+(k-1)][j][i].Fr1;
		pGrid->U[ke+k][j][i].Fr2 	=  pGrid->U[ks+(k-1)][j][i].Fr2;
		pGrid->U[ke+k][j][i].Fr3 	=  pGrid->U[ks+(k-1)][j][i].Fr3;
		pGrid->U[ke+k][j][i].Edd_11  	=  pGrid->U[ks+(k-1)][j][i].Edd_11;
		pGrid->U[ke+k][j][i].Edd_21  	=  pGrid->U[ks+(k-1)][j][i].Edd_21;
		pGrid->U[ke+k][j][i].Edd_22  	=  pGrid->U[ks+(k-1)][j][i].Edd_22;
		pGrid->U[ke+k][j][i].Edd_31  	=  pGrid->U[ks+(k-1)][j][i].Edd_31;
		pGrid->U[ke+k][j][i].Edd_32  	=  pGrid->U[ks+(k-1)][j][i].Edd_32;
		pGrid->U[ke+k][j][i].Edd_33  	=  pGrid->U[ks+(k-1)][j][i].Edd_33;
		pGrid->U[ke+k][j][i].Sigma_t 	=  pGrid->U[ks+(k-1)][j][i].Sigma_t;		
		pGrid->U[ke+k][j][i].Sigma_a 	=  pGrid->U[ks+(k-1)][j][i].Sigma_a;
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Inner x1 boundary (bc_ix1=5) */
/* conductor boundary condition is the same as reflect boundary condition for 
 * radiation variables *
 */

static void conduct_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
		pGrid->U[k][j][is-i].Er  	=  pGrid->U[k][j][is+(i-1)].Er;
		pGrid->U[k][j][is-i].Fr1 	= -pGrid->U[k][j][is+(i-1)].Fr1; /* reflect 1-flux. */
		pGrid->U[k][j][is-i].Fr2  	=  pGrid->U[k][j][is+(i-1)].Fr2;
		pGrid->U[k][j][is-i].Fr3  	=  pGrid->U[k][j][is+(i-1)].Fr3;
		pGrid->U[k][j][is-i].Edd_11  	=  pGrid->U[k][j][is+(i-1)].Edd_11;
		pGrid->U[k][j][is-i].Edd_21  	=  pGrid->U[k][j][is+(i-1)].Edd_21;
		pGrid->U[k][j][is-i].Edd_22  	=  pGrid->U[k][j][is+(i-1)].Edd_22;
		pGrid->U[k][j][is-i].Edd_31  	=  pGrid->U[k][j][is+(i-1)].Edd_31;
		pGrid->U[k][j][is-i].Edd_32  	=  pGrid->U[k][j][is+(i-1)].Edd_32;
		pGrid->U[k][j][is-i].Edd_33  	=  pGrid->U[k][j][is+(i-1)].Edd_33;
		pGrid->U[k][j][is-i].Sigma_t  	=  pGrid->U[k][j][is+(i-1)].Sigma_t;
		pGrid->U[k][j][is-i].Sigma_a  	=  pGrid->U[k][j][is+(i-1)].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x1 boundary (bc_ox1=5) */

static void conduct_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
		pGrid->U[k][j][ie+i].Er    	=  pGrid->U[k][j][ie-(i-1)].Er;
		pGrid->U[k][j][ie+i].Fr1 	= -pGrid->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */
		pGrid->U[k][j][ie+i].Fr2    	=  pGrid->U[k][j][ie-(i-1)].Fr2;
		pGrid->U[k][j][ie+i].Fr3    	=  pGrid->U[k][j][ie-(i-1)].Fr3;
		pGrid->U[k][j][ie+i].Edd_11  	=  pGrid->U[k][j][ie-(i-1)].Edd_11;
		pGrid->U[k][j][ie+i].Edd_21  	=  pGrid->U[k][j][ie-(i-1)].Edd_21;
		pGrid->U[k][j][ie+i].Edd_22  	=  pGrid->U[k][j][ie-(i-1)].Edd_22;
		pGrid->U[k][j][ie+i].Edd_31  	=  pGrid->U[k][j][ie-(i-1)].Edd_31;
		pGrid->U[k][j][ie+i].Edd_32  	=  pGrid->U[k][j][ie-(i-1)].Edd_32;
		pGrid->U[k][j][ie+i].Edd_33  	=  pGrid->U[k][j][ie-(i-1)].Edd_33;
		pGrid->U[k][j][ie+i].Sigma_t    =  pGrid->U[k][j][ie-(i-1)].Sigma_t;
		pGrid->U[k][j][ie+i].Sigma_a    =  pGrid->U[k][j][ie-(i-1)].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Inner x2 boundary (bc_ix2=5) */

static void conduct_ix2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
      		pGrid->U[k][js-j][i].Er   	=  pGrid->U[k][js+(j-1)][i].Er;
		pGrid->U[k][js-j][i].Fr1  	=  pGrid->U[k][js+(j-1)][i].Fr1; 
		pGrid->U[k][js-j][i].Fr2  	= -pGrid->U[k][js+(j-1)][i].Fr2; /* reflect 2-flux. */
		pGrid->U[k][js-j][i].Fr3  	=  pGrid->U[k][js+(j-1)][i].Fr3; 
		pGrid->U[k][js-j][i].Edd_11  	=  pGrid->U[k][js+(j-1)][i].Edd_11;
		pGrid->U[k][js-j][i].Edd_21  	=  pGrid->U[k][js+(j-1)][i].Edd_21;
		pGrid->U[k][js-j][i].Edd_22  	=  pGrid->U[k][js+(j-1)][i].Edd_22;
		pGrid->U[k][js-j][i].Edd_31  	=  pGrid->U[k][js+(j-1)][i].Edd_31;
		pGrid->U[k][js-j][i].Edd_32  	=  pGrid->U[k][js+(j-1)][i].Edd_32;
		pGrid->U[k][js-j][i].Edd_33  	=  pGrid->U[k][js+(j-1)][i].Edd_33;
		pGrid->U[k][js-j][i].Sigma_t  	=  pGrid->U[k][js+(j-1)][i].Sigma_t;
		pGrid->U[k][js-j][i].Sigma_a  	=  pGrid->U[k][js+(j-1)][i].Sigma_a;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x2 boundary (bc_ox2=5) */

static void conduct_ox2(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[k][je+j][i].Er   	=  pGrid->U[k][je-(j-1)][i].Er;
		pGrid->U[k][je+j][i].Fr1  	=  pGrid->U[k][je-(j-1)][i].Fr1; 
		pGrid->U[k][je+j][i].Fr2  	= -pGrid->U[k][je-(j-1)][i].Fr2; /* reflect 2-flux. */
		pGrid->U[k][je+j][i].Fr3  	=  pGrid->U[k][je-(j-1)][i].Fr3; 
		pGrid->U[k][je+j][i].Edd_11  	=  pGrid->U[k][je-(j-1)][i].Edd_11;
		pGrid->U[k][je+j][i].Edd_21  	=  pGrid->U[k][je-(j-1)][i].Edd_21;
		pGrid->U[k][je+j][i].Edd_22  	=  pGrid->U[k][je-(j-1)][i].Edd_22;
		pGrid->U[k][je+j][i].Edd_31  	=  pGrid->U[k][je-(j-1)][i].Edd_31;
		pGrid->U[k][je+j][i].Edd_32  	=  pGrid->U[k][je-(j-1)][i].Edd_32;
		pGrid->U[k][je+j][i].Edd_33  	=  pGrid->U[k][je-(j-1)][i].Edd_33;
		pGrid->U[k][je+j][i].Sigma_t  	=  pGrid->U[k][je-(j-1)][i].Sigma_t;	
		pGrid->U[k][je+j][i].Sigma_a  	=  pGrid->U[k][je-(j-1)][i].Sigma_a;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Inner x3 boundary (bc_ix3=5) */

static void conduct_ix3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[ks-k][j][i].Er   	=  pGrid->U[ks+(k-1)][j][i].Er;
		pGrid->U[ks-k][j][i].Fr1  	=  pGrid->U[ks+(k-1)][j][i].Fr1; 
		pGrid->U[ks-k][j][i].Fr2  	=  pGrid->U[ks+(k-1)][j][i].Fr2; 
		pGrid->U[ks-k][j][i].Fr3  	= -pGrid->U[ks+(k-1)][j][i].Fr3; /* reflect 3-flux. */
		pGrid->U[ks-k][j][i].Edd_11  	=  pGrid->U[ks+(k-1)][j][i].Edd_11;
		pGrid->U[ks-k][j][i].Edd_21  	=  pGrid->U[ks+(k-1)][j][i].Edd_21;
		pGrid->U[ks-k][j][i].Edd_22  	=  pGrid->U[ks+(k-1)][j][i].Edd_22;
		pGrid->U[ks-k][j][i].Edd_31  	=  pGrid->U[ks+(k-1)][j][i].Edd_31;
		pGrid->U[ks-k][j][i].Edd_32  	=  pGrid->U[ks+(k-1)][j][i].Edd_32;
		pGrid->U[ks-k][j][i].Edd_33  	=  pGrid->U[ks+(k-1)][j][i].Edd_33;
		pGrid->U[ks-k][j][i].Sigma_t  	=  pGrid->U[ks+(k-1)][j][i].Sigma_t;
		pGrid->U[ks-k][j][i].Sigma_a  	=  pGrid->U[ks+(k-1)][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x3 boundary (bc_ox3=5) */

static void conduct_ox3(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ke = pGrid->ke;
  int i,j,k;

  for (k=1; k<=nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		pGrid->U[ke+k][j][i].Er   	=  pGrid->U[ke-(k-1)][j][i].Er;
		pGrid->U[ke+k][j][i].Fr1  	=  pGrid->U[ke-(k-1)][j][i].Fr1; 
		pGrid->U[ke+k][j][i].Fr2  	=  pGrid->U[ke-(k-1)][j][i].Fr2; 
		pGrid->U[ke+k][j][i].Fr3  	= -pGrid->U[ke-(k-1)][j][i].Fr3; /* reflect 3-flux. */
		pGrid->U[ke+k][j][i].Edd_11  	=  pGrid->U[ke-(k-1)][j][i].Edd_11;
		pGrid->U[ke+k][j][i].Edd_21  	=  pGrid->U[ke-(k-1)][j][i].Edd_21;
		pGrid->U[ke+k][j][i].Edd_22  	=  pGrid->U[ke-(k-1)][j][i].Edd_22;
		pGrid->U[ke+k][j][i].Edd_31  	=  pGrid->U[ke-(k-1)][j][i].Edd_31;
		pGrid->U[ke+k][j][i].Edd_32  	=  pGrid->U[ke-(k-1)][j][i].Edd_32;
		pGrid->U[ke+k][j][i].Edd_33  	=  pGrid->U[ke-(k-1)][j][i].Edd_33;
		pGrid->U[ke+k][j][i].Sigma_t  	=  pGrid->U[ke-(k-1)][j][i].Sigma_t;
		pGrid->U[ke+k][j][i].Sigma_a  	=  pGrid->U[ke-(k-1)][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(GridS *pGrid)
{
  return;
}

#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~800 lines */
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(nghost-1); i++){
        *(pSnd++) = pG->U[k][j][i].Er;
        *(pSnd++) = pG->U[k][j][i].Fr1;
        *(pSnd++) = pG->U[k][j][i].Fr2;
        *(pSnd++) = pG->U[k][j][i].Fr3;
        *(pSnd++) = pG->U[k][j][i].Edd_11;
	*(pSnd++) = pG->U[k][j][i].Edd_21;
	*(pSnd++) = pG->U[k][j][i].Edd_22;
        *(pSnd++) = pG->U[k][j][i].Edd_31;
        *(pSnd++) = pG->U[k][j][i].Edd_32;
        *(pSnd++) = pG->U[k][j][i].Edd_33;
        *(pSnd++) = pG->U[k][j][i].Sigma_t;
        *(pSnd++) = pG->U[k][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(nghost-1); i<=ie; i++){
		*(pSnd++) = pG->U[k][j][i].Er;
        	*(pSnd++) = pG->U[k][j][i].Fr1;
        	*(pSnd++) = pG->U[k][j][i].Fr2;
        	*(pSnd++) = pG->U[k][j][i].Fr3;
        	*(pSnd++) = pG->U[k][j][i].Edd_11;
		*(pSnd++) = pG->U[k][j][i].Edd_21;
		*(pSnd++) = pG->U[k][j][i].Edd_22;
        	*(pSnd++) = pG->U[k][j][i].Edd_31;
        	*(pSnd++) = pG->U[k][j][i].Edd_32;
        	*(pSnd++) = pG->U[k][j][i].Edd_33;
        	*(pSnd++) = pG->U[k][j][i].Sigma_t;
        	*(pSnd++) = pG->U[k][j][i].Sigma_a;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=js+(nghost-1); j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
		*(pSnd++) = pG->U[k][j][i].Er;
        	*(pSnd++) = pG->U[k][j][i].Fr1;
        	*(pSnd++) = pG->U[k][j][i].Fr2;
        	*(pSnd++) = pG->U[k][j][i].Fr3;
        	*(pSnd++) = pG->U[k][j][i].Edd_11;
		*(pSnd++) = pG->U[k][j][i].Edd_21;
		*(pSnd++) = pG->U[k][j][i].Edd_22;
        	*(pSnd++) = pG->U[k][j][i].Edd_31;
        	*(pSnd++) = pG->U[k][j][i].Edd_32;
        	*(pSnd++) = pG->U[k][j][i].Edd_33;
        	*(pSnd++) = pG->U[k][j][i].Sigma_t;
        	*(pSnd++) = pG->U[k][j][i].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=je-(nghost-1); j<=je; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
		*(pSnd++) = pG->U[k][j][i].Er;
        	*(pSnd++) = pG->U[k][j][i].Fr1;
        	*(pSnd++) = pG->U[k][j][i].Fr2;
        	*(pSnd++) = pG->U[k][j][i].Fr3;
        	*(pSnd++) = pG->U[k][j][i].Edd_11;
		*(pSnd++) = pG->U[k][j][i].Edd_21;
		*(pSnd++) = pG->U[k][j][i].Edd_22;
        	*(pSnd++) = pG->U[k][j][i].Edd_31;
        	*(pSnd++) = pG->U[k][j][i].Edd_32;
        	*(pSnd++) = pG->U[k][j][i].Edd_33;
        	*(pSnd++) = pG->U[k][j][i].Sigma_t;
        	*(pSnd++) = pG->U[k][j][i].Sigma_a;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ks+(nghost-1); k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
     		*(pSnd++) = pG->U[k][j][i].Er;
        	*(pSnd++) = pG->U[k][j][i].Fr1;
        	*(pSnd++) = pG->U[k][j][i].Fr2;
        	*(pSnd++) = pG->U[k][j][i].Fr3;
        	*(pSnd++) = pG->U[k][j][i].Edd_11;
		*(pSnd++) = pG->U[k][j][i].Edd_21;
		*(pSnd++) = pG->U[k][j][i].Edd_22;
        	*(pSnd++) = pG->U[k][j][i].Edd_31;
        	*(pSnd++) = pG->U[k][j][i].Edd_32;
        	*(pSnd++) = pG->U[k][j][i].Edd_33;
        	*(pSnd++) = pG->U[k][j][i].Sigma_t;
        	*(pSnd++) = pG->U[k][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ke-(nghost-1); k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
     		*(pSnd++) = pG->U[k][j][i].Er;
        	*(pSnd++) = pG->U[k][j][i].Fr1;
        	*(pSnd++) = pG->U[k][j][i].Fr2;
        	*(pSnd++) = pG->U[k][j][i].Fr3;
        	*(pSnd++) = pG->U[k][j][i].Edd_11;
		*(pSnd++) = pG->U[k][j][i].Edd_21;
		*(pSnd++) = pG->U[k][j][i].Edd_22;
        	*(pSnd++) = pG->U[k][j][i].Edd_31;
        	*(pSnd++) = pG->U[k][j][i].Edd_32;
        	*(pSnd++) = pG->U[k][j][i].Edd_33;
        	*(pSnd++) = pG->U[k][j][i].Sigma_t;
        	*(pSnd++) = pG->U[k][j][i].Sigma_a;
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is-nghost; i<=is-1; i++){
        pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
 	pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=js-1; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=je+1; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks-nghost; k<=ks-1; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
	pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ke = pG->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ke+1; k<=ke+nghost; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
  	pG->U[k][j][i].Er  	= *(pRcv++);
        pG->U[k][j][i].Fr1 	= *(pRcv++);
        pG->U[k][j][i].Fr2 	= *(pRcv++);
        pG->U[k][j][i].Fr3 	= *(pRcv++);
        pG->U[k][j][i].Edd_11  	= *(pRcv++);
        pG->U[k][j][i].Edd_21 	= *(pRcv++);
        pG->U[k][j][i].Edd_22 	= *(pRcv++);
        pG->U[k][j][i].Edd_31 	= *(pRcv++);
	pG->U[k][j][i].Edd_32 	= *(pRcv++);
	pG->U[k][j][i].Edd_33 	= *(pRcv++);
	pG->U[k][j][i].Sigma_t 	= *(pRcv++);
	pG->U[k][j][i].Sigma_a 	= *(pRcv++);
      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */


#endif /* End radiation hydro and radiation MHD */



