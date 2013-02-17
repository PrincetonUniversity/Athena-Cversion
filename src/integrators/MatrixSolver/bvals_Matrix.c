#include "../../copyright.h"
/*==============================================================================
 * FILE: bvals_Matrix.c
 * Set boundary condition for matrix 
 * 
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../../prototypes.h"

#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)

#if defined(MATRIX_MULTIGRID) || defined(MATRIX_HYPRE)

#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double **send_buf = NULL, **recv_buf = NULL;
static MPI_Request *recv_rq, *send_rq;
static int Nrad = 4 * Matghost; /* Only need to update boundary condition for Er, Fr? for matrix solver, ghost zone numbers can be larger than 1 */
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

static void reflect_ix1(MatrixS *pMat);
static void reflect_ox1(MatrixS *pMat);
static void reflect_ix2(MatrixS *pMat);
static void reflect_ox2(MatrixS *pMat);
static void reflect_ix3(MatrixS *pMat);
static void reflect_ox3(MatrixS *pMat);

static void outflow_ix1(MatrixS *pMat);
static void outflow_ox1(MatrixS *pMat);
static void outflow_ix2(MatrixS *pMat);
static void outflow_ox2(MatrixS *pMat);
static void outflow_ix3(MatrixS *pMat);
static void outflow_ox3(MatrixS *pMat);

static void periodic_ix1(MatrixS *pMat);
static void periodic_ox1(MatrixS *pMat);
static void periodic_ix2(MatrixS *pMat);
static void periodic_ox2(MatrixS *pMat);
static void periodic_ix3(MatrixS *pMat);
static void periodic_ox3(MatrixS *pMat);

static void conduct_ix1(MatrixS *pMat);
static void conduct_ox1(MatrixS *pMat);
static void conduct_ix2(MatrixS *pMat);
static void conduct_ox2(MatrixS *pMat);
static void conduct_ix3(MatrixS *pMat);
static void conduct_ox3(MatrixS *pMat);

static void ProlongateLater(MatrixS *pMat);

#ifdef MPI_PARALLEL
static void pack_ix1(MatrixS *pMat);
static void pack_ox1(MatrixS *pMat);
static void pack_ix2(MatrixS *pMat);
static void pack_ox2(MatrixS *pMat);
static void pack_ix3(MatrixS *pMat);
static void pack_ox3(MatrixS *pMat);

static void unpack_ix1(MatrixS *pMat);
static void unpack_ox1(MatrixS *pMat);
static void unpack_ix2(MatrixS *pMat);
static void unpack_ox2(MatrixS *pMat);
static void unpack_ix3(MatrixS *pMat);
static void unpack_ox3(MatrixS *pMat);
#endif /* MPI_PARALLEL */


/*===========================*/
/* Function pointer for this Matrix object */
/* Initialize them to be NULL */

static VMatFun_t Mat_ix1_BCFun = NULL;
static VMatFun_t Mat_ox1_BCFun = NULL;  /* ix1/ox1 BC function pointers for this Dom for radiation quantities*/
static VMatFun_t Mat_ix2_BCFun = NULL;
static VMatFun_t Mat_ox2_BCFun = NULL;  /* ix1/ox1 BC function pointers for this Dom for radiation quantities*/
static VMatFun_t Mat_ix3_BCFun = NULL;
static VMatFun_t Mat_ox3_BCFun = NULL;  /* ix1/ox1 BC function pointers for this Dom for radiation quantities*/

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 * This function is form domain now. Should call for every domain in mesh 
 * Only 16 radiation related variables are set here. Other variables in Matghost
 * zones are set in function bvals_MHD 
 */
/* Initialize the boundary condition functions according to the functions in Mesh */

void bvals_Matrix(MatrixS *pMat)
{

#ifdef MPI_PARALLEL
  int cnt, cnt2, cnt3, ierr, mIndex;
#endif /* MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  if (pMat->Nx[0] > 1){

#ifdef MPI_PARALLEL
    cnt = Matghost*(pMat->Nx[1])*(pMat->Nx[2])*(Nrad);


/* MPI blocks to both left and right */
    if (pMat->rx1_id >= 0 && pMat->lx1_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix1(pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox1(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pMat);
      if (mIndex == 1) unpack_ox1(pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(pMat);
      if (mIndex == 1) unpack_ox1(pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx1_id >= 0 && pMat->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatRtoL_tag,
       pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(Mat_ix1_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1(pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx1_id < 0 && pMat->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1(pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(Mat_ox1_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1(pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx1_id < 0 && pMat->lx1_id < 0) {
      (*(Mat_ix1_BCFun))(pMat);
      (*(Mat_ox1_BCFun))(pMat);
    } 

  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  if (pMat->Nx[1] > 1){

#ifdef MPI_PARALLEL
    cnt = (pMat->Nx[0] + 2*Matghost)*Matghost*(pMat->Nx[2])*(Nrad);


/* MPI blocks to both left and right */
    if (pMat->rx2_id >= 0 && pMat->lx2_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix2(pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox2(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pMat);
      if (mIndex == 1) unpack_ox2(pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(pMat);
      if (mIndex == 1) unpack_ox2(pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx2_id >= 0 && pMat->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(Mat_ix2_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2(pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx2_id < 0 && pMat->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2(pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(Mat_ox2_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2(pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx2_id < 0 && pMat->lx2_id < 0) {
      (*(Mat_ix2_BCFun))(pMat);
      (*(Mat_ox2_BCFun))(pMat);
    } 

/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */
/* This is done after periodic boundary condition has been applied */
#ifdef SHEARING_BOX
    if (pMat->my_iproc == 0 && pMat->BCFlag_ix1 == 4) {

      ShearingSheet_Matrix_ix1(pMat);
    }
   
    if (pMat->my_iproc == ((pMat->NGrid[0])-1) && pMat->BCFlag_ox1 == 4) {
      ShearingSheet_Matrix_ox1(pMat);
    }
#endif

  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pMat->Nx[2] > 1){

#ifdef MPI_PARALLEL
    cnt = (pMat->Nx[0] + 2*Matghost)*(pMat->Nx[1] + 2*Matghost)*Matghost*(Nrad);


/* MPI blocks to both left and right */
    if (pMat->rx3_id >= 0 && pMat->lx3_id >= 0) {

      /* Post non-blocking receives for data from L and R Grids */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data L and R */
      pack_ix3(pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox3(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pMat);
      if (mIndex == 1) unpack_ox3(pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(pMat);
      if (mIndex == 1) unpack_ox3(pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx3_id >= 0 && pMat->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3(pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      (*(Mat_ix3_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3(pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx3_id < 0 && pMat->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3(pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      (*(Mat_ox3_BCFun))(pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3(pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx3_id < 0 && pMat->lx3_id < 0) {
      (*(Mat_ix3_BCFun))(pMat);
      (*(Mat_ox3_BCFun))(pMat);
    } 

  }

  return;
}



void bvals_Matrix_init(MatrixS *pMat)
{

/* MPI flag, lx1, rx1, rx2, lx2, rx3, lx3 are all set when */
/* Matrix object is created */ 

#ifdef MPI_PARALLEL
  int nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */

/* Cycle through all the Domains that have active Grids on this proc */


/* Set function pointers for physical boundaries in x1-direction -------------*/

    if(pMat->Nx[0] > 1) {

/*---- ix1 boundary ----------------------------------------------------------*/

            switch(pMat->BCFlag_ix1){
	
	    case 0: /* In case the grid does not touch the boundary of the root domain, we do not need to */ 
		    /* change the boundary at each relaxation as the ghost zones are fixed to be zero	 */
		    /* It will get updated every cycle */
	      Mat_ix1_BCFun = ProlongateLater;
	    break;

            case 1: /* Reflecting, B_normal=0 */
              Mat_ix1_BCFun = reflect_ix1;
            break;

            case 2: /* Outflow */
              Mat_ix1_BCFun = outflow_ix1;
            break;
		
	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ix1(&(Mat_ix1_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ix1_BCFun = periodic_ix1;
            break;

            case 5: /* Reflecting, B_normal!=0 */
              Mat_ix1_BCFun = conduct_ix1;
            break;

            default:
              ath_perr(-1,"[bvals_Matrix_init]:bc_ix1=%d unknown\n",pMat->BCFlag_ix1);
              exit(EXIT_FAILURE);
            }
	
         

/*---- ox1 boundary ----------------------------------------------------------*/

            switch(pMat->BCFlag_ox1){

	    case 0:
	      Mat_ox1_BCFun = ProlongateLater;
	    break;

            case 1: /* Reflecting, B_normal=0 */
              Mat_ox1_BCFun = reflect_ox1;
            break;

            case 2: /* Outflow */
              Mat_ox1_BCFun = outflow_ox1;
            break;

	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ox1(&(Mat_ox1_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ox1_BCFun = periodic_ox1;
            break;

            case 5: /* Reflecting, B_normal!=0 */
              Mat_ox1_BCFun = conduct_ox1;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox1=%d unknown\n",pMat->BCFlag_ox1);
              exit(EXIT_FAILURE);
 
          }
        
 	}

/* Set function pointers for physical boundaries in x2-direction -------------*/

    if(pMat->Nx[1] > 1) {
/*---- ix2 boundary ----------------------------------------------------------*/

      

            switch(pMat->BCFlag_ix2){

	    case 0:
	      Mat_ix2_BCFun = ProlongateLater;
	    break;

            case 1: /* Reflecting, B_normal=0 */
              Mat_ix2_BCFun = reflect_ix2;
            break;

            case 2: /* Outflow */
              Mat_ix2_BCFun = outflow_ix2;
            break;
		
	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ix2(&(Mat_ix2_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ix2_BCFun = periodic_ix2;
            break;
  
            case 5: /* Reflecting, B_normal!=0 */
              Mat_ix2_BCFun = conduct_ix2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix2=%d unknown\n",pMat->BCFlag_ix2);
              exit(EXIT_FAILURE);
            }
          
        
      

/*---- ox2 boundary ----------------------------------------------------------*/


            switch(pMat->BCFlag_ox2){

	    case 0: 
	      Mat_ox2_BCFun = ProlongateLater;
	    break;

            case 1: /* Reflecting, B_normal=0 */
              Mat_ox2_BCFun = reflect_ox2;
            break;

            case 2: /* Outflow */
              Mat_ox2_BCFun = outflow_ox2;
            break;

	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ox2(&(Mat_ox2_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ox2_BCFun = periodic_ox2;
            break;

            case 5: /* Reflecting, B_normal!=0 */
              Mat_ox2_BCFun = conduct_ox2;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox2=%d unknown\n",pMat->BCFlag_ox2);
              exit(EXIT_FAILURE);
            }
          
    }

/* Set function pointers for physical boundaries in x3-direction -------------*/

    if(pMat->Nx[2] > 1) {

/*---- ix3 boundary ----------------------------------------------------------*/

            switch(pMat->BCFlag_ix3){

 	    case 0:
	      Mat_ix3_BCFun = ProlongateLater;
	    break;        

            case 1: /* Reflecting, B_normal=0 */
              Mat_ix3_BCFun = reflect_ix3;
            break;

            case 2: /* Outflow */
              Mat_ix3_BCFun = outflow_ix3;
            break;

	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ix3(&(Mat_ix3_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ix3_BCFun = periodic_ix3;
            break;

            case 5: /* Reflecting, B_normal!=0 */
              Mat_ix3_BCFun = conduct_ix3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ix3=%d unknown\n",pMat->BCFlag_ix3);
              exit(EXIT_FAILURE);
            }
          
	

/*---- ox3 boundary ----------------------------------------------------------*/

            switch(pMat->BCFlag_ox3){

	    case 0:
	      Mat_ox3_BCFun = ProlongateLater;
	    break;

            case 1: /* Reflecting, B_normal=0 */
              Mat_ox3_BCFun = reflect_ox3;
            break;

            case 2: /* Outflow */
              Mat_ox3_BCFun = outflow_ox3;
            break;

	    case 3: /* inflow, function is set in problem generator */
	      bvals_mat_fun_ox3(&(Mat_ox3_BCFun));	
	    break;

            case 4: /* Periodic. Handle with MPI calls for parallel jobs. */
              Mat_ox3_BCFun = periodic_ox3;
            break;

            case 5: /* Reflecting, B_normal!=0 */
              Mat_ox3_BCFun = conduct_ox3;
            break;

            default:
              ath_perr(-1,"[bvals_init]:bc_ox3=%d unknown\n",pMat->BCFlag_ox3);
              exit(EXIT_FAILURE);
            }
          
	}
/* Figure out largest size needed for send/receive buffers with MPI ----------*/

#ifdef MPI_PARALLEL

/* Assuming every grid has the same size */

/* x1cnt is surface area of x1 faces */	
	if(pMat->NGrid[0] > 1){

		nx2t = pMat->Nx[1];
		nx3t = pMat->Nx[2];

		if(nx2t > 1) nx2t += 1;
		if(nx3t > 1) nx3t += 1;

		if(nx2t*nx3t > x1cnt) x1cnt = nx2t*nx3t;

	}  

/* x2cnt is surface area of x2 faces */

	if(pMat->NGrid[1] > 1){
		  nx1t = pMat->Nx[0];
	  if(nx1t > 1) nx1t += 2*Matghost;

	  nx3t = pMat->Nx[2];
	  if(nx3t > 1) nx3t += 1;

          if(nx1t*nx3t > x2cnt) x2cnt = nx1t*nx3t;
	}

/* x3cnt is surface area of x3 faces */
	if(pMat->NGrid[2] > 1){
	  nx1t = pMat->Nx[0];
	  if(nx1t > 1) nx1t += 2*Matghost;

	  nx2t = pMat->Nx[1];
	  if(nx2t > 1) nx2t += 2*Matghost;

          if(nx1t*nx2t > x3cnt) x3cnt = nx1t*nx2t;
	}
     
#endif /* MPI_PARALLEL */

  


#ifdef MPI_PARALLEL
/* Allocate memory for send/receive buffers and MPI_Requests */

  size = x1cnt > x2cnt ? x1cnt : x2cnt;
  size = x3cnt >  size ? x3cnt : size;
 /* Here we only need to send Matrix related variables *
  * total is Er, Fr???, Sigma_??, Edd_?? V?, T4 14+NOPACITY variables 
  */

  size *= Matghost*Nrad;


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

void bvals_Matrix_destruct(MatrixS *pMat)
{

/* Only need to do this for MPI case to delete the buffar */
/* Because the matrix will be destroyed at each level. We need to clean this */

#ifdef MPI_PARALLEL
	if(send_buf != NULL){
		free_2d_array(send_buf);
		send_buf = NULL;
	}

	if(recv_buf != NULL){
		free_2d_array(recv_buf);
		recv_buf = NULL;
	}

	if(recv_rq != NULL){
		free_1d_array(recv_rq);
		recv_rq = NULL;
	}

	if(send_rq != NULL){
		free_1d_array(send_rq);
		send_rq = NULL;
	}


#endif

}


void bvals_Matrix_fun(MatrixS *pMat, enum BCDirection dir, VMatFun_t prob_bc)
{
  switch(dir){
  case left_x1:
    Mat_ix1_BCFun = prob_bc;
    break;
  case right_x1:
    Mat_ox1_BCFun = prob_bc;
    break;
  case left_x2:
    Mat_ix2_BCFun = prob_bc;
    break;
  case right_x2:
    Mat_ox2_BCFun = prob_bc;
    break;
  case left_x3:
    Mat_ix3_BCFun = prob_bc;
    break;
  case right_x3:
    Mat_ox3_BCFun = prob_bc;
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

static void reflect_ix1(MatrixS *pMat)
{
  int is = pMat->is;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		pMat->U[k][j][is-i]	=  pMat->U[k][j][is+(i-1)];
		/* reflect velocity and flux */
		pMat->U[k][j][is-i].Fr1 = -pMat->U[k][j][is+(i-1)].Fr1;		

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

				pMat->U[k][j][is-i].Er = 0.0;
				pMat->U[k][j][is-i].Fr1 = 0.0;
				pMat->U[k][j][is-i].Fr2 = 0.0;
				pMat->U[k][j][is-i].Fr3 = 0.0;
		}
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x1 boundary (bc_ox1=1) */

static void reflect_ox1(MatrixS *pMat)
{
  int ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
      		pMat->U[k][j][ie+i]	=  pMat->U[k][j][ie-(i-1)];
		/* reflect the velocity and flux */
		pMat->U[k][j][ie+i].Fr1 = -pMat->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */		

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

				pMat->U[k][j][ie+i].Er = 0.0;
				pMat->U[k][j][ie+i].Fr1 = 0.0;
				pMat->U[k][j][ie+i].Fr2 = 0.0;
				pMat->U[k][j][ie+i].Fr3 = 0.0;
		}

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x2 boundary (bc_ix2=1) */

static void reflect_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[k][js-j][i]	=  pMat->U[k][js+(j-1)][i];
		/* reflect velocity and flux */
		pMat->U[k][js-j][i].Fr2 = -pMat->U[k][js+(j-1)][i].Fr2;		

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

				pMat->U[k][js-j][i].Er = 0.0;
				pMat->U[k][js-j][i].Fr1 = 0.0;
				pMat->U[k][js-j][i].Fr2 = 0.0;
				pMat->U[k][js-j][i].Fr3 = 0.0;
		}

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x2 boundary (bc_ox2=1) */

static void reflect_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[k][je+j][i]	=  pMat->U[k][je-(j-1)][i];
		/* reflect the velocity and flux */
		pMat->U[k][je+j][i].Fr2 = -pMat->U[k][je-(j-1)][i].Fr2; /* reflect 1-flux. */
		

		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

				pMat->U[k][je+j][i].Er = 0.0;
				pMat->U[k][je+j][i].Fr1 = 0.0;
				pMat->U[k][je+j][i].Fr2 = 0.0;
				pMat->U[k][je+j][i].Fr3 = 0.0;
		}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Inner x3 boundary (bc_ix3=1) */

static void reflect_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ks-k][j][i]	=  pMat->U[ks+(k-1)][j][i];
		/* reflect velocity and flux */
		pMat->U[ks-k][j][i].Fr3 = -pMat->U[ks+(k-1)][j][i].Fr3;
		
		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

				pMat->U[ks-k][j][i].Er = 0.0;
				pMat->U[ks-k][j][i].Fr1 = 0.0;
				pMat->U[ks-k][j][i].Fr2 = 0.0;
				pMat->U[ks-k][j][i].Fr3 = 0.0;
		}

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* REFLECTING boundary conditions, Outer x3 boundary (bc_ox3=1) */

static void reflect_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ke = pMat->ke;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ke+k][j][i]	=  pMat->U[ke-(k-1)][j][i];
		/* reflect the velocity and flux */
		pMat->U[ke+k][j][i].Fr3 = -pMat->U[ke-(k-1)][j][i].Fr3; /* reflect 1-flux. */
		

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

				pMat->U[ke+k][j][i].Er = 0.0;
				pMat->U[ke+k][j][i].Fr1 = 0.0;
				pMat->U[ke+k][j][i].Fr2 = 0.0;
				pMat->U[ke+k][j][i].Fr3 = 0.0;
		}

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary condition, Inner x1 boundary (bc_ix1=2) */

static void outflow_ix1(MatrixS *pMat)
{
  int is = pMat->is;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
        	pMat->U[k][j][is-i]	=  pMat->U[k][j][is];

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

				pMat->U[k][j][is-i].Er = 0.0;
				pMat->U[k][j][is-i].Fr1 = 0.0;
				pMat->U[k][j][is-i].Fr2 = 0.0;
				pMat->U[k][j][is-i].Fr3 = 0.0;
		}

      }
    }
  }



  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x1 boundary (bc_ox1=2) */

static void outflow_ox1(MatrixS *pMat)
{
  int ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
        	pMat->U[k][j][ie+i]	=  pMat->U[k][j][ie];

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

				pMat->U[k][j][ie+i].Er = 0.0;
				pMat->U[k][j][ie+i].Fr1 = 0.0;
				pMat->U[k][j][ie+i].Fr2 = 0.0;
				pMat->U[k][j][ie+i].Fr3 = 0.0;
		}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x2 boundary (bc_ix2=2) */

static void outflow_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[k][js-j][i]	=  pMat->U[k][js][i];

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

				pMat->U[k][js-j][i].Er = 0.0;
				pMat->U[k][js-j][i].Fr1 = 0.0;
				pMat->U[k][js-j][i].Fr2 = 0.0;
				pMat->U[k][js-j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x2 boundary (bc_ox2=2) */

static void outflow_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[k][je+j][i]	=  pMat->U[k][je][i];

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

				pMat->U[k][je+j][i].Er = 0.0;
				pMat->U[k][je+j][i].Fr1 = 0.0;
				pMat->U[k][je+j][i].Fr2 = 0.0;
				pMat->U[k][je+j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Inner x3 boundary (bc_ix3=2) */

static void outflow_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ks-k][j][i]	=  pMat->U[ks][j][i];

		/* in case background state is subtracted or in coarse grid */
		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

				pMat->U[ks-k][j][i].Er = 0.0;
				pMat->U[ks-k][j][i].Fr1 = 0.0;
				pMat->U[ks-k][j][i].Fr2 = 0.0;
				pMat->U[ks-k][j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions, Outer x3 boundary (bc_ox3=2) */

static void outflow_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ke = pMat->ke;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ke+k][j][i]	=  pMat->U[ke][j][i];

		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

			pMat->U[ke+k][j][i].Er = 0.0;
			pMat->U[ke+k][j][i].Fr1 = 0.0;
			pMat->U[ke+k][j][i].Fr2 = 0.0;
			pMat->U[ke+k][j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
       		pMat->U[k][j][is-i] 	=  pMat->U[k][j][ie-(i-1)];
		
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
        	pMat->U[k][j][ie+i] 	=  pMat->U[k][j][is+(i-1)];
		
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		pMat->U[k][js-j][i] 	=  pMat->U[k][je-(j-1)][i];
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
      		pMat->U[k][je+j][i] 	=  pMat->U[k][js+(j-1)][i];
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		pMat->U[ks-k][j][i] 	=  pMat->U[ke-(k-1)][j][i];
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		pMat->U[ke+k][j][i] 	=  pMat->U[ks+(k-1)][j][i];	
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

static void conduct_ix1(MatrixS *pMat)
{
  int is = pMat->is;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		pMat->U[k][j][is-i]	=  pMat->U[k][j][is+(i-1)];
		/* reflect velocity and flux */
		pMat->U[k][j][is-i].Fr1 = -pMat->U[k][j][is+(i-1)].Fr1;
		

		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

			pMat->U[k][j][is-i].Er = 0.0;
			pMat->U[k][j][is-i].Fr1 = 0.0;
			pMat->U[k][j][is-i].Fr2 = 0.0;
			pMat->U[k][j][is-i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x1 boundary (bc_ox1=5) */

static void conduct_ox1(MatrixS *pMat)
{
  int ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		pMat->U[k][j][ie+i]	=  pMat->U[k][j][ie-(i-1)];
		/* reflect the velocity and flux */
		pMat->U[k][j][ie+i].Fr1 = -pMat->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */
		

		if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){

			pMat->U[k][j][ie+i].Er = 0.0;
			pMat->U[k][j][ie+i].Fr1 = 0.0;
			pMat->U[k][j][ie+i].Fr2 = 0.0;
			pMat->U[k][j][ie+i].Fr3 = 0.0;
		}

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Inner x2 boundary (bc_ix2=5) */

static void conduct_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
      		pMat->U[k][js-j][i]	=  pMat->U[k][js+(j-1)][i];
		/* reflect velocity and flux */
		pMat->U[k][js-j][i].Fr2 = -pMat->U[k][js+(j-1)][i].Fr2;
		
		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

			pMat->U[k][js-j][i].Er = 0.0;
			pMat->U[k][js-j][i].Fr1 = 0.0;
			pMat->U[k][js-j][i].Fr2 = 0.0;
			pMat->U[k][js-j][i].Fr3 = 0.0;
		}

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x2 boundary (bc_ox2=5) */

static void conduct_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;


  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[k][je+j][i]	=  pMat->U[k][je-(j-1)][i];
		/* reflect the velocity and flux */
		pMat->U[k][je+j][i].Fr2 = -pMat->U[k][je-(j-1)][i].Fr2; /* reflect 1-flux. */
		

		if((pMat->bgflag) || (pMat->Nx[1] < pMat->RootNx[1])){

			pMat->U[k][je+j][i].Er = 0.0;
			pMat->U[k][je+j][i].Fr1 = 0.0;
			pMat->U[k][je+j][i].Fr2 = 0.0;
			pMat->U[k][je+j][i].Fr3 = 0.0;
		}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Inner x3 boundary (bc_ix3=5) */

static void conduct_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ks-k][j][i]	=  pMat->U[ks+(k-1)][j][i];
		/* reflect velocity and flux */
		pMat->U[ks-k][j][i].Fr3 = -pMat->U[ks+(k-1)][j][i].Fr3;
		

		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

			pMat->U[ks-k][j][i].Er = 0.0;
			pMat->U[ks-k][j][i].Fr1 = 0.0;
			pMat->U[ks-k][j][i].Fr2 = 0.0;
			pMat->U[ks-k][j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* CONDUCTOR boundary conditions, Outer x3 boundary (bc_ox3=5) */

static void conduct_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ke = pMat->ke;
  int i,j,k;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
        	pMat->U[ke+k][j][i]	=  pMat->U[ke-(k-1)][j][i];
		/* reflect the velocity and flux */
		pMat->U[ke+k][j][i].Fr3 = -pMat->U[ke-(k-1)][j][i].Fr3; /* reflect 1-flux. */
		

		if((pMat->bgflag) || (pMat->Nx[2] < pMat->RootNx[2])){

			pMat->U[ke+k][j][i].Er = 0.0;
			pMat->U[ke+k][j][i].Fr1 = 0.0;
			pMat->U[ke+k][j][i].Fr2 = 0.0;
			pMat->U[ke+k][j][i].Fr3 = 0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PROLONGATION boundary conditions.  Nothing is actually done here, the
 * prolongation is actually handled in ProlongateGhostZones in main loop, so
 * this is just a NoOp Grid function.  */

static void ProlongateLater(MatrixS *pMat)
{
  return;
}



#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~800 lines */
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=is+(Matghost-1); i++){
        *(pSnd++) = pMat->U[k][j][i].Er;
        *(pSnd++) = pMat->U[k][j][i].Fr1;
        *(pSnd++) = pMat->U[k][j][i].Fr2;
        *(pSnd++) = pMat->U[k][j][i].Fr3;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;
  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=ie-(Matghost-1); i<=ie; i++){
		*(pSnd++) = pMat->U[k][j][i].Er;
        	*(pSnd++) = pMat->U[k][j][i].Fr1;
        	*(pSnd++) = pMat->U[k][j][i].Fr2;
        	*(pSnd++) = pMat->U[k][j][i].Fr3;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=js+(Matghost-1); j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		*(pSnd++) = pMat->U[k][j][i].Er;
        	*(pSnd++) = pMat->U[k][j][i].Fr1;
        	*(pSnd++) = pMat->U[k][j][i].Fr2;
        	*(pSnd++) = pMat->U[k][j][i].Fr3;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ks; k<=ke; k++){
    for (j=je-(Matghost-1); j<=je; j++){
      for (i=is-Matghost; i<=ie+Matghost; i++){
		*(pSnd++) = pMat->U[k][j][i].Er;
        	*(pSnd++) = pMat->U[k][j][i].Fr1;
        	*(pSnd++) = pMat->U[k][j][i].Fr2;
        	*(pSnd++) = pMat->U[k][j][i].Fr3;		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[0][0]);

  for (k=ks; k<=ks+(Matghost-1); k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
     		*(pSnd++) = pMat->U[k][j][i].Er;
        	*(pSnd++) = pMat->U[k][j][i].Fr1;
        	*(pSnd++) = pMat->U[k][j][i].Fr2;
        	*(pSnd++) = pMat->U[k][j][i].Fr3;

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pSnd;
  pSnd = (double*)&(send_buf[1][0]);

  for (k=ke-(Matghost-1); k<=ke; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
     		*(pSnd++) = pMat->U[k][j][i].Er;
        	*(pSnd++) = pMat->U[k][j][i].Fr1;
        	*(pSnd++) = pMat->U[k][j][i].Fr2;
        	*(pSnd++) = pMat->U[k][j][i].Fr3;
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1(MatrixS *pMat)
{
  int is = pMat->is;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is-Matghost; i<=is-1; i++){
        pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1(MatrixS *pMat)
{
  int ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=ie+1; i<=ie+Matghost; i++) {
 	pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);


      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks; k<=ke; k++) {
    for (j=js-Matghost; j<=js-1; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
	pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ks; k<=ke; k++) {
    for (j=je+1; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
	pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[0][0]);

  for (k=ks-Matghost; k<=ks-1; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
	pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);	

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3(MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ke = pMat->ke;
  int i,j,k;

  double *pRcv;
  pRcv = (double*)&(recv_buf[1][0]);

  for (k=ke+1; k<=ke+Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
  	pMat->U[k][j][i].Er  	= *(pRcv++);
        pMat->U[k][j][i].Fr1 	= *(pRcv++);
        pMat->U[k][j][i].Fr2 	= *(pRcv++);
        pMat->U[k][j][i].Fr3 	= *(pRcv++);	

      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */

#endif /* MATRIX_MULTIGRID */

#endif /* End radiation hydro and radiation MHD */



