#include "copyright.h"
/*==============================================================================
 * FILE: set_bvals_mhd.c
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
 *       bc_ix1 = Inner Boundary Condition for x1
 *       bc_ox1 = Outer Boundary Condition for x1
 *   similarly for bc_ix2; bc_ox2; bc_ix3; bc_ox3
 *
 * For case (1) -- PHYSICAL BOUNDARIES
 *   The values of the integer flags (bc_ix1, etc.) are:
 *       1 = reflecting, B_normal = 0; 2 = outflow; 4 = periodic
 *       5 = reflecting, B_normal != 0
 *   Following ZEUS conventions, 3 would be flow-in (ghost zones held at
 *   pre-determined fixed values), however in Athena instead we use pointers to
 *   user-defined BC functions for flow-in.
 *
 * For case (2) -- MPI BOUNDARIES
 *   We do the parallel synchronization by having every grid:
 *     1) Pack and send data to the grid on right  [send_ox1()]
 *     2) Listen to the left, unpack and set data  [receive_ix1()]
 *     3) Pack and send data to the grid on left   [send_ix1()]
 *     4) Listen to the right, unpack and set data [receive_ox1()]
 *   If the grid is at the edge of the Domain, we set BCs as in case (1) or (3).
 *   Currently the code uses NON-BLOCKING sends (MPI_Isend) and BLOCKING
 *   receives (MPI_Recv).  Some optimization could be achieved by interleaving
 *   non-blocking sends (MPI_Isend) and computations.
 *
 * For case (3) -- INTERNAL GRID LEVEL BOUNDARIES
 *   This step is complicated and must be handled separately, in the function
 *   Prolongate() called from the main loop.  In the algorithm below, nothing is
 *   done for this case; the BCs are left to be set later.
 *
 * Which case of BC is needed is unchanged throughout a calculation.  Thus,
 * during setup we determine which case we need, and set a pointer to the
 * appropriate BC function using set_bvals_init().  These pointers are stored in
 * the Domain structure containing the Grid, since the BCs can be different for
 * each different domain.  MPI calls are used when the proc ID number to the
 * left or right is >= 0.
 * 
 * With SELF-GRAVITY: BCs for Phi are set independently of the MHD variables
 *   in a separate function set_bvals_grav()
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   set_bvals_mhd()      - calls appropriate functions to set ghost cells
 *   set_bvals_mhd_init() - sets function pointers used by set_bvals_mhd()
 *   set_bvals_mhd_fun()  - enrolls a pointer to a user-defined BC function
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
/* MPI send and receive buffer, size dynamically determined near end of
 * set_bvals_init() based on number of zones in each grid */
static double *send_buf = NULL, *recv_buf = NULL;

/* NVAR_SHARE = maximim number of variables passed in any one MPI message =
 * variables in ConsVarS, plus 3 extra for interface magnetic fields. */
#ifdef MHD
#define NVAR_SHARE (NVAR + 3)
#else
#define NVAR_SHARE NVAR
#endif /* MHD */

#endif /* MPI_PARALLEL */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   reflect_???()  - reflecting BCs at boundary ???
 *   outflow_???()  - outflow BCs at boundary ???
 *   periodic_???() - periodic BCs at boundary ???
 *   send_???()     - MPI send of data at ??? boundary
 *   receive_???()  - MPI receive of data at ??? boundary
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

static void ProlongateLater(GridS *pG);

#ifdef MPI_PARALLEL
static void send_ix1(DomainS *pD);
static void send_ox1(DomainS *pD);
static void send_ix2(DomainS *pD);
static void send_ox2(DomainS *pD);
static void send_ix3(DomainS *pD);
static void send_ox3(DomainS *pD);

static void receive_ix1(GridS *pG, MPI_Request *prq);
static void receive_ox1(GridS *pG, MPI_Request *prq);
static void receive_ix2(GridS *pG, MPI_Request *prq);
static void receive_ox2(GridS *pG, MPI_Request *prq);
static void receive_ix3(GridS *pG, MPI_Request *prq);
static void receive_ox3(GridS *pG, MPI_Request *prq);
#endif /* MPI_PARALLEL */

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* set_bvals_mhd: calls appropriate functions to set ghost zones.  The function
 *   pointers (*(pD->???_BCFun)) are set by set_bvals_init() to be either a
 *   user-defined function, or one of the functions corresponding to reflecting,
 *   periodic, or outflow.  If the left- or right-Grid ID numbers are >= 1
 *   (neighboring grids exist), then MPI calls are used.
 *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

void set_bvals_mhd(DomainS *pD)
{
  GridS *pGrid = (pD->Grid);
#ifdef SHEARING_BOX
  int myL,myM,myN;
#endif
#ifdef MPI_PARALLEL
  int cnt1, cnt2, cnt3, cnt, ierr;
  MPI_Request rq;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pGrid->Nx[0] > 1){

#ifdef MPI_PARALLEL
    cnt2 = pGrid->Nx[1] > 1 ? pGrid->Nx[1] + 1 : 1;
    cnt3 = pGrid->Nx[2] > 1 ? pGrid->Nx[2] + 1 : 1;
    cnt = nghost*cnt2*cnt3*NVAR_SHARE;

/* MPI blocks to both left and right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox1(pD);             /* send R */
      receive_ix1(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix1(pD);             /* send L */
      receive_ox1(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx1_id >= 0 && pGrid->lx1_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx1_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox1(pD);              /* send R */
      (*(pD->ix1_BCFun))(pGrid);
      receive_ox1 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx1_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix1(pD);              /* send L */
      (*(pD->ox1_BCFun))(pGrid);
      receive_ix1 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx1_id < 0 && pGrid->lx1_id < 0) {
      (*(pD->ix1_BCFun))(pGrid);
      (*(pD->ox1_BCFun))(pGrid);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pGrid->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx[0] > 1 ? pGrid->Nx[0] + 2*nghost : 1;
    cnt3 = pGrid->Nx[2] > 1 ? pGrid->Nx[2] + 1 : 1;
    cnt = nghost*cnt1*cnt3*NVAR_SHARE;

/* MPI blocks to both left and right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox2(pD);             /* send R */
      receive_ix2(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix2(pD);             /* send L */
      receive_ox2(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx2_id >= 0 && pGrid->lx2_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx2_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox2(pD);              /* send R */
      (*(pD->ix2_BCFun))(pGrid);
      receive_ox2 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx2_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix2(pD);              /* send L */
      (*(pD->ox2_BCFun))(pGrid);
      receive_ix2 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx2_id < 0 && pGrid->lx2_id < 0) {
      (*(pD->ix2_BCFun))(pGrid);
      (*(pD->ox2_BCFun))(pGrid);
    }

/* shearing sheet BCs; function defined in problem generator */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
    if (myL == 0) {
      ShearingSheet_ix1(pD);
    }
    if (myL == ((pD->NGrid[0])-1)) {
      ShearingSheet_ox1(pD);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pGrid->Nx[2] > 1){

#ifdef MPI_PARALLEL
    cnt1 = pGrid->Nx[0] > 1 ? pGrid->Nx[0] + 2*nghost : 1;
    cnt2 = pGrid->Nx[1] > 1 ? pGrid->Nx[1] + 2*nghost : 1;
    cnt = nghost*cnt1*cnt2*NVAR_SHARE;

/* MPI blocks to both left and right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox3(pD);             /* send R */
      receive_ix3(pGrid, &rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix3(pD);             /* send L */
      receive_ox3(pGrid, &rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pGrid->rx3_id >= 0 && pGrid->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->rx3_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ox3(pD);              /* send R */
      (*(pD->ix3_BCFun))(pGrid);
      receive_ox3 (pGrid, &rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pGrid->lx3_id,
		      boundary_cells_tag, pD->Comm_Domain, &rq);

      send_ix3(pD);              /* send L */
      (*(pD->ox3_BCFun))(pGrid);
      receive_ix3 (pGrid, &rq);  /* listen L */
    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pGrid->rx3_id < 0 && pGrid->lx3_id < 0) {
      (*(pD->ix3_BCFun))(pGrid);
      (*(pD->ox3_BCFun))(pGrid);
    }

  }

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_init:  sets function pointers for physical boundaries during
 *   initialization, allocates memory for send/receive buffers with MPI
 */

void set_bvals_mhd_init(MeshS *pM)
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

      if(pD->ix1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[0] != 0) {      
          pD->ix1_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) && 
               pM->BCFlag_ix1 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix1b but not ox1b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {                    
            switch(pM->BCFlag_ix1){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix1_BCFun = reflect_ix1;
            break;

            case 2: /* Outflow */
              pD->ix1_BCFun = outflow_ix1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix1_BCFun = periodic_ix1;
#ifdef MPI_PARALLEL
              if(pG->lx1_id < 0 && pD->NGrid[0] > 1){
                pG->lx1_id = pD->GData[myN][myM][pD->NGrid[0]-1].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ix1_BCFun = reflect_ix1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix1=%d unknown\n",pM->BCFlag_ix1);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox1 boundary ----------------------------------------------------------*/

      if(pD->ox1_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[0] + pD->Nx[0])/irefine != pM->Nx[0]) {
          pD->ox1_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[0] != 0) && (pM->BCFlag_ox1 == 4)) {      
            ath_error("[bvals_init]:level=%d Domain touching ox1b but not ix1b and periodic BC not allowed\n",nl); 


/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox1){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox1_BCFun = reflect_ox1;
            break;

            case 2: /* Outflow */
              pD->ox1_BCFun = outflow_ox1;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox1_BCFun = periodic_ox1;
#ifdef MPI_PARALLEL
              if(pG->rx1_id < 0 && pD->NGrid[0] > 1){
                pG->rx1_id = pD->GData[myN][myM][0].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox1_BCFun = reflect_ox1;
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

      if(pD->ix2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[1] != 0) {
          pD->ix2_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) &&
               pM->BCFlag_ix2 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix2b but not ox2b and periodic BC not allowed\n",nl); 


/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix2){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix2_BCFun = reflect_ix2;
            break;

            case 2: /* Outflow */
              pD->ix2_BCFun = outflow_ix2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix2_BCFun = periodic_ix2;
#ifdef MPI_PARALLEL
              if(pG->lx2_id < 0 && pD->NGrid[1] > 1){
                pG->lx2_id = pD->GData[myN][pD->NGrid[1]-1][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;
  
            case 5: /* Reflecting, B_normal!=0 */
              pD->ix2_BCFun = reflect_ix2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix2=%d unknown\n",pM->BCFlag_ix2);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox2 boundary ----------------------------------------------------------*/

      if(pD->ox2_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[1] + pD->Nx[1])/irefine != pM->Nx[1]) {
          pD->ox2_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[1] != 0) && (pM->BCFlag_ox2 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox2b but not ix2b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox2){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox2_BCFun = reflect_ox2;
            break;

            case 2: /* Outflow */
              pD->ox2_BCFun = outflow_ox2;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox2_BCFun = periodic_ox2;
#ifdef MPI_PARALLEL
              if(pG->rx2_id < 0 && pD->NGrid[1] > 1){
                pG->rx2_id = pD->GData[myN][0][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox2_BCFun = reflect_ox2;
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

      if(pD->ix3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if(pD->Disp[2] != 0) {
          pD->ix3_BCFun = ProlongateLater;

/* Domain is at L-edge of root Domain, but not R-edge and periodic BC  */
        } else {
          if(((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) &&
               pM->BCFlag_ix3 == 4) {
            ath_error("[bvals_init]:level=%d Domain touching ix3b but not ox3b and periodic BC not allowed\n",nl); 

/* Domain is at L-edge of root Domain */
          } else {
            switch(pM->BCFlag_ix3){

            case 1: /* Reflecting, B_normal=0 */
              pD->ix3_BCFun = reflect_ix3;
            break;

            case 2: /* Outflow */
              pD->ix3_BCFun = outflow_ix3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ix3_BCFun = periodic_ix3;
#ifdef MPI_PARALLEL
              if(pG->lx3_id < 0 && pD->NGrid[2] > 1){
                pG->lx3_id = pD->GData[pD->NGrid[2]-1][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ix3_BCFun = reflect_ix3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix3=%d unknown\n",pM->BCFlag_ix3);
              exit(EXIT_FAILURE);
            }
          }
        }
      }

/*---- ox3 boundary ----------------------------------------------------------*/

      if(pD->ox3_BCFun == NULL){    /* BCFun ptr was not set in prob gen */

/* Domain boundary is in interior of root */
        if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2]) {
          pD->ox3_BCFun = ProlongateLater;

/* Domain is at R-edge of root Domain, but not L-edge and periodic BC */
        } else {
          if((pD->Disp[2] != 0) && (pM->BCFlag_ox3 == 4)) {
            ath_error("[bvals_init]:level=%d Domain touching ox3b but not ix3b and periodic BC not allowed\n",nl); 

/* Domain is at R-edge of root Domain */
          } else {
            switch(pM->BCFlag_ox3){

            case 1: /* Reflecting, B_normal=0 */
              pD->ox3_BCFun = reflect_ox3;
            break;

            case 2: /* Outflow */
              pD->ox3_BCFun = outflow_ox3;
            break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              pD->ox3_BCFun = periodic_ox3;
#ifdef MPI_PARALLEL
              if(pG->rx3_id < 0 && pD->NGrid[2] > 1){
                pG->rx3_id = pD->GData[0][myM][myL].ID_Comm_Domain;
              }
#endif /* MPI_PARALLEL */
            break;

            case 5: /* Reflecting, B_normal!=0 */
              pD->ox3_BCFun = reflect_ox3;
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
	if(pD->NGrid[0] > 1){
	  nx2t = pD->GData[n][m][l].Nx[1];
	  if(nx2t > 1) nx2t += 1;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;
	}

	if(pD->NGrid[1] > 1){
	  nx1t = pD->GData[n][m][l].Nx[0];
	  if(nx1t > 1) nx1t += 2*nghost;

	  nx3t = pD->GData[n][m][l].Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
	}


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
/* Allocate memory for send/receive buffers in MPI parallel calculations */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;

  size *= nghost; /* Multiply by the third dimension */

  if (size > 0) {
    if((send_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate send buffer\n");

    if((recv_buf = (double*)malloc(size*NVAR_SHARE*sizeof(double))) == NULL)
      ath_error("[set_bvals_init]: Failed to allocate receive buffer\n");
  }
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_fun:  sets function pointers for user-defined BCs in problem file
 */

void set_bvals_mhd_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    pD->ix1_BCFun = prob_bc;
    break;
  case right_x1:
    pD->ox1_BCFun = prob_bc;
    break;
  case left_x2:
    pD->ix2_BCFun = prob_bc;
    break;
  case right_x2:
    pD->ox2_BCFun = prob_bc;
    break;
  case left_x3:
    pD->ix3_BCFun = prob_bc;
    break;
  case right_x3:
    pD->ox3_BCFun = prob_bc;
    break;
  default:
    ath_perr(-1,"[set_bvals_fun]: Unknown direction = %d\n",dir);
    exit(EXIT_FAILURE);
  }
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/* Following are the functions:
 *   reflecting_???:   where ???=[ix1,ox1,ix2,ox2,ix3,ox3]
 *   outflow_???
 *   periodic_???
 *   send_???
 *   receive_???
 */

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x1 boundary (bc_ix1=1,5)
 */

static void reflect_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int bc_ix1,ju,ku; /* j-upper, k-upper */
  Real qa;
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i]    =  pGrid->U[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].M1 = -pGrid->U[k][j][is-i].M1; /* reflect 1-mom. */
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][is-i]    =  pGrid->W[k][j][is+(i-1)];
        pGrid->W[k][j][is-i].V1 = -pGrid->W[k][j][is-i].V1; /* reflect 1-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ix1=1) */
  bc_ix1 = par_geti("grid","bc_ix1");
  if (bc_ix1 == 1) qa = -1.0;
  if (bc_ix1 == 5) qa =  1.0;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      if (bc_ix1 == 1) pGrid->B1i[k][j][is] = 0.0;
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i]   = qa*pGrid->B1i[k][j][is+i];
        pGrid->U[k][j][is-i].B1c = qa*pGrid->U[k][j][is+(i-1)].B1c;
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i]   = -qa*pGrid->B2i[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].B2c = -qa*pGrid->U[k][j][is+(i-1)].B2c;
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i]   = -qa*pGrid->B3i[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].B3c = -qa*pGrid->U[k][j][is+(i-1)].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1,5)
 */

static void reflect_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int bc_ox1,ju,ku; /* j-upper, k-upper */
  Real qa;
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i]    =  pGrid->U[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].M1 = -pGrid->U[k][j][ie+i].M1; /* reflect 1-mom. */
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][ie+i]    =  pGrid->W[k][j][ie-(i-1)];
        pGrid->W[k][j][ie+i].V1 = -pGrid->W[k][j][ie+i].V1; /* reflect 1-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ox1=1) */
  bc_ox1 = par_geti("grid","bc_ox1");
  if (bc_ox1 == 1) qa = -1.0;
  if (bc_ox1 == 5) qa =  1.0;

/* i=ie+1 is not set for the interface field B1i, except bc_ox1=1 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      if (bc_ox1 == 1 ) pGrid->B1i[k][j][ie+1] = 0.0;
      pGrid->U[k][j][ie+1].B1c = qa*pGrid->U[k][j][ie].B1c;
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i]   = qa*pGrid->B1i[k][j][ie-(i-2)];
        pGrid->U[k][j][ie+i].B1c = qa*pGrid->U[k][j][ie-(i-1)].B1c;
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i]   = -qa*pGrid->B2i[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].B2c = -qa*pGrid->U[k][j][ie-(i-1)].B2c;
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i]   = -qa*pGrid->B3i[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].B3c = -qa*pGrid->U[k][j][ie-(i-1)].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1,5)
 */

static void reflect_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int bc_ix2,ku; /* k-upper */
  Real qa;
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][js-j][i]    =  pGrid->W[k][js+(j-1)][i];
        pGrid->W[k][js-j][i].V2 = -pGrid->W[k][js-j][i].V2; /* reflect 2-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ix2=1) */
  bc_ix2 = par_geti("grid","bc_ix2");
  if (bc_ix2 == 1) qa = -1.0;
  if (bc_ix2 == 5) qa =  1.0;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][js-j][i]   = -qa*pGrid->B1i[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].B1c = -qa*pGrid->U[k][js+(j-1)][i].B1c;
      }
    }
  }

  for (k=ks; k<=ke; k++) {
    if (bc_ix2 == 1) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js][i] = 0.0;
      }
    }
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][js-j][i]   = qa*pGrid->B2i[k][js+j][i];
        pGrid->U[k][js-j][i].B2c = qa*pGrid->U[k][js+(j-1)][i].B2c;
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i]   = -qa*pGrid->B3i[k][js+(j-1)][i];
        pGrid->U[k][js-j][i].B3c = -qa*pGrid->U[k][js+(j-1)][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1,5)
 */

static void reflect_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int bc_ox2,ku; /* k-upper */
  Real qa;
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][je+j][i]    =  pGrid->W[k][je-(j-1)][i];
        pGrid->W[k][je+j][i].V2 = -pGrid->W[k][je+j][i].V2; /* reflect 2-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ox2=1) */
  bc_ox2 = par_geti("grid","bc_ox2");
  if (bc_ox2 == 1) qa = -1.0;
  if (bc_ox2 == 5) qa =  1.0;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i]   = -qa*pGrid->B1i[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].B1c = -qa*pGrid->U[k][je-(j-1)][i].B1c;
      }
    }
  }

/* j=je+1 is not set for the interface field B2i, except bc_ox2=1 */
  for (k=ks; k<=ke; k++) {
    for (i=il; i<=iu; i++) {
      if (bc_ox2 == 1) pGrid->B2i[k][je+1][i] = 0.0;
      pGrid->U[k][je+1][i].B2c = qa*pGrid->U[k][je][i].B2c;
    }
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i]   = qa*pGrid->B2i[k][je-(j-2)][i];
        pGrid->U[k][je+j][i].B2c = qa*pGrid->U[k][je-(j-1)][i].B2c;
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i]   = -qa*pGrid->B3i[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].B3c = -qa*pGrid->U[k][je-(j-1)][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1,5)
 */

static void reflect_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
#ifdef MHD
  int bc_ix3;
  Real qa;
#endif

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ks-k][j][i]    =  pGrid->W[ks+(k-1)][j][i];
        pGrid->W[ks-k][j][i].V3 = -pGrid->W[ks-k][j][i].V3; /* reflect 3-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ix3=1) */
  bc_ix3 = par_geti("grid","bci_x3");
  if (bc_ix3 == 1) qa = -1.0;
  if (bc_ix3 == 5) qa =  1.0;

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i]   = -qa*pGrid->B1i[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].B1c = -qa*pGrid->U[ks+(k-1)][j][i].B1c;
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i]   = -qa*pGrid->B2i[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].B2c = -qa*pGrid->U[ks+(k-1)][j][i].B2c;
      }
    }
  }

  if (bc_ix3 == 1) {
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      pGrid->B3i[ks][j][i] = 0.0;
    }
  }}
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i]   = qa*pGrid->B3i[ks+k][j][i];
        pGrid->U[ks-k][j][i].B3c = qa*pGrid->U[ks+(k-1)][j][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1,5)
 */

static void reflect_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
#ifdef MHD
  int bc_ox3;
  Real qa;
#endif

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ke+k][j][i]    =  pGrid->W[ke-(k-1)][j][i];
        pGrid->W[ke+k][j][i].V3 = -pGrid->W[ke+k][j][i].V3; /* reflect 3-vel. */
#endif
      }
    }
  }

#ifdef MHD
/* The multiplier qa=-1 if B_normal=0 (bc_ox3=1) */
  bc_ox3 = par_geti("grid","bc_ox3");
  if (bc_ox3 == 1) qa = -1.0;
  if (bc_ox3 == 5) qa =  1.0;

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i]   = -qa*pGrid->B1i[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].B1c = -qa*pGrid->U[ke-(k-1)][j][i].B1c;
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i]   = -qa*pGrid->B2i[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].B2c = -qa*pGrid->U[ke-(k-1)][j][i].B2c;
      }
    }
  }

/* k=ke+1 is not set for the interface field B3i, except bc_ox3=1 */
  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      if (bc_ox3 == 1) pGrid->B3i[ke+1][j][i] = 0.0;
      pGrid->U[ke+1][j][i].B3c = qa*pGrid->U[ke][j][i].B3c;
    }
  }
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i]   = qa*pGrid->B3i[ke-(k-2)][j][i];
        pGrid->U[ke+k][j][i].B3c = qa*pGrid->U[ke-(k-1)][j][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditionss, Inner x1 boundary (bc_ix1=2)
 */

static void outflow_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][is-i] = pGrid->W[k][j][is];
#endif
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

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2)
 */

static void outflow_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][ie+i] = pGrid->W[k][j][ie];
#endif
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

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2)
 */

static void outflow_ix2(GridS *pGrid)
{
  int js = pGrid->js;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][js-j][i] = pGrid->W[k][js][i];
#endif
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

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2)
 */

static void outflow_ox2(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][je+j][i] = pGrid->W[k][je][i];
#endif
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

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2)
 */

static void outflow_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ks][j][i];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ks-k][j][i] = pGrid->W[ks][j][i];
#endif
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
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2)
 */

static void outflow_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ke+k][j][i] = pGrid->W[ke][j][i];
#endif
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
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4)
 */

static void periodic_ix1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][ie-(i-1)];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][is-i] = pGrid->W[k][j][ie-(i-1)];
#endif
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

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][ie-(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][ie-(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4)
 */

static void periodic_ox1(GridS *pGrid)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][is+(i-1)];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][j][ie+i] = pGrid->W[k][j][is+(i-1)];
#endif
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

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][is+(i-1)];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][is+(i-1)];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4)
 */

static void periodic_ix2(GridS *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][js-j][i] = pGrid->W[k][je-(j-1)][i];
#endif
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

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][js-j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4)
 */

static void periodic_ox2(GridS *pGrid)
{
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
#ifdef MHD
  int ku; /* k-upper */
#endif

  if (pGrid->Nx[0] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[k][je+j][i] = pGrid->W[k][js+(j-1)][i];
#endif
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

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][js+(j-1)][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4)
 */

static void periodic_ix3(GridS *pGrid)
{
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i] = pGrid->U[ke-(k-1)][j][i];
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ks-k][j][i] = pGrid->W[ke-(k-1)][j][i];
#endif
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
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4)
 */

static void periodic_ox3(GridS *pGrid)
{
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx[0] > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx[1] > 1){
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
#ifdef SPECIAL_RELATIVITY
        pGrid->W[ke+k][j][i] = pGrid->W[ks+(k-1)][j][i];
#endif
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
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp function.
 */

static void ProlongateLater(GridS *pGrid)
{
  return;
}


#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~760 lines */

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x1 boundary -- send left
 */

static void send_ix1(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  il = pG->is;
  iu = pG->is + nghost - 1;

  if(pG->Nx[1] > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in ConsVarS structure into send buffer */

  /* Following expression gives same cnt as in Step 1 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x1 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id,
		 boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x1 boundary -- send right
 */

static void send_ox1(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  il = pG->ie - nghost + 1;
  iu = pG->ie;

  if(pG->Nx[1] > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in ConsVarS structure into send buffer */

  /* Following expression gives same cnt as in Step 1 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x1 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id,
		 boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x2 boundary -- send left
 */

static void send_ix2(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->js;
  ju = pG->js + nghost - 1;

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in ConsVarS structure into send buffer */

  /* Following expression gives same cnt as in Step 2 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x2 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
		 boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x2 boundary -- send right
 */

static void send_ox2(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->je - nghost + 1;
  ju = pG->je;

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Pack data in ConsVarS structure into send buffer */

  /* Following expression gives same cnt as in Step 2 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x2 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
		 boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Inner x3 boundary -- send left
 */

static void send_ix3(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx[1] > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ks;
  ku = pG->ks + nghost - 1;

/* Pack data in ConsVarS structure into send buffer */

  /* Following expression gives same cnt as in Step 3 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on L-x3 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx3_id,
		  boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SEND of boundary conditions, Outer x3 boundary -- send right
 */

static void send_ox3(DomainS *pD)
{
  int i,il,iu,j,jl,ju,k,kl,ku,cnt,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pSnd = send_buf;
  GridS *pG=pD->Grid;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx[1] > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ke - nghost + 1;
  ku = pG->ke;

/* Pack data in ConsVarS structure into send buffer */

    /* Following expression gives same cnt as in Step 3 in set_bvals()  */
  cnt = (iu-il+1)*(ju-jl+1)*(ku-kl+1)*NVAR_SHARE;
  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        *(pSnd++) = pCons->d;
        *(pSnd++) = pCons->M1;
        *(pSnd++) = pCons->M2;
        *(pSnd++) = pCons->M3;
#ifndef BAROTROPIC
        *(pSnd++) = pCons->E;
#endif /* BAROTROPIC */
#ifdef MHD
        *(pSnd++) = pCons->B1c;
        *(pSnd++) = pCons->B2c;
        *(pSnd++) = pCons->B3c;
        *(pSnd++) = pG->B1i[k][j][i];
        *(pSnd++) = pG->B2i[k][j][i];
        *(pSnd++) = pG->B3i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) *(pSnd++) = pCons->s[n];
#endif
      }
    }
  }

/* send contents of buffer to the neighboring grid on R-x3 */

  ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx3_id,
		  boundary_cells_tag, pD->Comm_Domain);

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x1 boundary -- listen left
 */

static void receive_ix1(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  il = pG->is - nghost;
  iu = pG->is - 1;

  if(pG->Nx[1] > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the left grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
        pG->B1i[k][j][i] = *(pRcv++);
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x1 boundary -- listen right
 */

static void receive_ox1(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  il = pG->ie + 1;
  iu = pG->ie + nghost;

  if(pG->Nx[1] > 1){
    jl = pG->js;
    ju = pG->je + 1;
  } else {
    jl = ju = pG->js;
  }

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the right grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
/* Do not set B1i[ie+1] for shearing sheet boundary conditions */
#ifdef SHEARING_BOX
        if (i>il) {pG->B1i[k][j][i] = *(pRcv++);}
        else {pRcv++;}
#else
        pG->B1i[k][j][i] = *(pRcv++);
#endif /* SHEARING_BOX */
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x2 boundary -- listen left
 */

static void receive_ix2(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->js - nghost;
  ju = pG->js - 1;

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the left grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
        pG->B1i[k][j][i] = *(pRcv++);
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x2 boundary -- listen right
 */

static void receive_ox2(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  jl = pG->je + 1;
  ju = pG->je + nghost;

  if(pG->Nx[2] > 1){
    kl = pG->ks;
    ku = pG->ke + 1;
  } else {
    kl = ku = pG->ks;
  }

/* Wait to receive the input data from the right grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
        pG->B1i[k][j][i] = *(pRcv++);
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Inner x3 boundary -- listen left
 */

static void receive_ix3(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx[1] > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ks - nghost;
  ku = pG->ks - 1;

/* Wait to receive the input data from the left grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
        pG->B1i[k][j][i] = *(pRcv++);
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_RECEIVE of boundary conditions, Outer x3 boundary -- listen right
 */

static void receive_ox3(GridS *pG, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,ierr;
#if (NSCALARS > 0)
  int n;
#endif
  ConsVarS *pCons;
  double *pRcv = recv_buf;

  if(pG->Nx[0] > 1){
    il = pG->is - nghost;
    iu = pG->ie + nghost;
  } else {
    il = iu = pG->is;
  }

  if(pG->Nx[1] > 1){
    jl = pG->js - nghost;
    ju = pG->je + nghost;
  } else {
    jl = ju = pG->js;
  }

  kl = pG->ke + 1;
  ku = pG->ke + nghost;

/* Wait to receive the input data from the right grid */

  ierr = MPI_Wait(prq, MPI_STATUS_IGNORE);

/* Manually unpack the data from the receive buffer */

  for (k=kl; k<=ku; k++){
    for (j=jl; j<=ju; j++){
      for (i=il; i<=iu; i++){
        /* Get a pointer to the ConsVarS cell */
        pCons = &(pG->U[k][j][i]);

        pCons->d  = *(pRcv++);
        pCons->M1 = *(pRcv++);
        pCons->M2 = *(pRcv++);
        pCons->M3 = *(pRcv++);
#ifndef BAROTROPIC
        pCons->E  = *(pRcv++);
#endif /* BAROTROPIC */
#ifdef MHD
        pCons->B1c = *(pRcv++);
        pCons->B2c = *(pRcv++);
        pCons->B3c = *(pRcv++);
        pG->B1i[k][j][i] = *(pRcv++);
        pG->B2i[k][j][i] = *(pRcv++);
        pG->B3i[k][j][i] = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) pCons->s[n] = *(pRcv++);
#endif
      }
    }
  }

  return;
}

#endif /* MPI_PARALLEL */
