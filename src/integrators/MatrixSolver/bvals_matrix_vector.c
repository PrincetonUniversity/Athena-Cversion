#include "../../copyright.h"
/*==============================================================================
 * FILE: bvals_matrix_vector.c
 * Set boundary condition for matrix vector product A.v
 * Only need to exchange values for MPI and periodic boundary conditions.
 * For other boundary conditions, ghost zones are set to be zero  
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
#endif /* MPI_PARALLEL */
static int Nrad = 4;

static void periodic_ix1(Real ****vector, MatrixS *pMat);
static void periodic_ox1(Real ****vector, MatrixS *pMat);
static void periodic_ix2(Real ****vector, MatrixS *pMat);
static void periodic_ox2(Real ****vector, MatrixS *pMat);
static void periodic_ix3(Real ****vector, MatrixS *pMat);
static void periodic_ox3(Real ****vector, MatrixS *pMat);

static void non_periodic_ix1(Real ****vector, MatrixS *pMat);
static void non_periodic_ox1(Real ****vector, MatrixS *pMat);
static void non_periodic_ix2(Real ****vector, MatrixS *pMat);
static void non_periodic_ox2(Real ****vector, MatrixS *pMat);
static void non_periodic_ix3(Real ****vector, MatrixS *pMat);
static void non_periodic_ox3(Real ****vector, MatrixS *pMat);



#ifdef MPI_PARALLEL
static void pack_ix1(Real ****vector, MatrixS *pMat);
static void pack_ox1(Real ****vector, MatrixS *pMat);
static void pack_ix2(Real ****vector, MatrixS *pMat);
static void pack_ox2(Real ****vector, MatrixS *pMat);
static void pack_ix3(Real ****vector, MatrixS *pMat);
static void pack_ox3(Real ****vector, MatrixS *pMat);

static void unpack_ix1(Real ****vector, MatrixS *pMat);
static void unpack_ox1(Real ****vector, MatrixS *pMat);
static void unpack_ix2(Real ****vector, MatrixS *pMat);
static void unpack_ox2(Real ****vector, MatrixS *pMat);
static void unpack_ix3(Real ****vector, MatrixS *pMat);
static void unpack_ox3(Real ****vector, MatrixS *pMat);
#endif /* MPI_PARALLEL */




/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* *
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 * This function is form domain now. Should call for every domain in mesh 
 */


void bvals_matrix_vector(Real ****vector, MatrixS *pMat)
{
	void bvals_matrix_vector_init(MatrixS *pMat);
	void bvals_matrix_vector_destruct(MatrixS *pMat);

#ifdef SHEARING_BOX
	extern void ShearingSheet_matrix_vector_ix1(Real ****vector, MatrixS *pMat);
	extern void ShearingSheet_matrix_vector_ox1(Real ****vector, MatrixS *pMat);
	extern void bvals_matrix_vector_shear_init(MatrixS *pMat);
	extern void bvals_matrix_vector_shear_destruct(void);
#endif


	/* First, initialize the memory */
	bvals_matrix_vector_init(pMat);
#ifdef SHEARING_BOX
  int myL,myM,myN,BCFlag;
	/* also need to initialize for shearing boundary condition */
	bvals_matrix_vector_shear_init(pMat);
#endif


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
      pack_ix1(vector, pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox1(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(vector, pMat);
      if (mIndex == 1) unpack_ox1(vector, pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix1(vector, pMat);
      if (mIndex == 1) unpack_ox1(vector, pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx1_id >= 0 && pMat->lx1_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox1(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      /* Only need to consider periodic or non-periodic bounadry condition */
	if(pMat->BCFlag_ix1 == 4)
		periodic_ix1(vector, pMat);
	else
		non_periodic_ix1(vector, pMat);


      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox1(vector, pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx1_id < 0 && pMat->lx1_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix1(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx1_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

     /* set physical boundary */
      /* Only need to consider periodic or non-periodic bounadry condition */
	if(pMat->BCFlag_ox1 == 4)
		periodic_ox1(vector, pMat);
	else
		non_periodic_ox1(vector, pMat);


      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix1(vector, pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx1_id < 0 && pMat->lx1_id < 0) {
      if(pMat->BCFlag_ix1 == 4)
		periodic_ix1(vector, pMat);
	else
		non_periodic_ix1(vector, pMat);

      if(pMat->BCFlag_ox1 == 4)
		periodic_ox1(vector, pMat);
	else
		non_periodic_ox1(vector, pMat);

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
      pack_ix2(vector, pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox2(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(vector, pMat);
      if (mIndex == 1) unpack_ox2(vector, pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix2(vector, pMat);
      if (mIndex == 1) unpack_ox2(vector, pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx2_id >= 0 && pMat->lx2_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox2(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
        if(pMat->BCFlag_ix2 == 4)
		periodic_ix2(vector, pMat);
	else
		non_periodic_ix2(vector, pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox2(vector, pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx2_id < 0 && pMat->lx2_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix2(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx2_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
        if(pMat->BCFlag_ox2 == 4)
		periodic_ox2(vector, pMat);
	else
		non_periodic_ox2(vector, pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix2(vector, pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx2_id < 0 && pMat->lx2_id < 0) {
       if(pMat->BCFlag_ix2 == 4)
		periodic_ix2(vector, pMat);
	else
		non_periodic_ix2(vector, pMat);

      if(pMat->BCFlag_ox2 == 4)
		periodic_ox2(vector, pMat);
	else
		non_periodic_ox2(vector, pMat);
    } 

/* shearing sheet BCs; function defined in problem generator.
 * Enroll outflow BCs if perdiodic BCs NOT selected.  This assumes the root
 * level grid is specified by the <domain1> block in the input file */
/* This is done after periodic boundary condition has been applied */
#ifdef SHEARING_BOX
    if (pMat->my_iproc == 0 && pMat->BCFlag_ix1 == 4) {

      ShearingSheet_matrix_vector_ix1(vector,pMat);
    }
   
    if (pMat->my_iproc == ((pMat->NGrid[0])-1) && pMat->BCFlag_ox1 == 4) {
      ShearingSheet_matrix_vector_ox1(vector,pMat);
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
      pack_ix3(vector, pMat);
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      pack_ox3(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* check non-blocking sends have completed. */
      ierr = MPI_Waitall(2, send_rq, MPI_STATUS_IGNORE);

      /* check non-blocking receives and unpack data in any order. */
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(vector, pMat);
      if (mIndex == 1) unpack_ox3(vector, pMat);
      ierr = MPI_Waitany(2,recv_rq,&mIndex,MPI_STATUS_IGNORE);
      if (mIndex == 0) unpack_ix3(vector, pMat);
      if (mIndex == 1) unpack_ox3(vector, pMat);

    }

/* Physical boundary on left, MPI block on right */
    if (pMat->rx3_id >= 0 && pMat->lx3_id < 0) {

      /* Post non-blocking receive for data from R Grid */
      ierr = MPI_Irecv(&(recv_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(recv_rq[1]));

      /* pack and send data R */
      pack_ox3(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[1][0]),cnt,MPI_DOUBLE,pMat->rx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(send_rq[1]));

      /* set physical boundary */
      if(pMat->BCFlag_ix3 == 4)
		periodic_ix3(vector, pMat);
	else
		non_periodic_ix3(vector, pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[1]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from R and unpack data */
      ierr = MPI_Wait(&(recv_rq[1]), MPI_STATUS_IGNORE);
      unpack_ox3(vector, pMat);

    }

/* MPI block on left, Physical boundary on right */
    if (pMat->rx3_id < 0 && pMat->lx3_id >= 0) {

      /* Post non-blocking receive for data from L grid */
      ierr = MPI_Irecv(&(recv_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatLtoR_tag,
        pMat->Comm_Domain, &(recv_rq[0]));

      /* pack and send data L */
      pack_ix3(vector, pMat); 
      ierr = MPI_Isend(&(send_buf[0][0]),cnt,MPI_DOUBLE,pMat->lx3_id,MatRtoL_tag,
        pMat->Comm_Domain, &(send_rq[0]));

      /* set physical boundary */
      if(pMat->BCFlag_ox3 == 4)
		periodic_ox3(vector, pMat);
	else
		non_periodic_ox3(vector, pMat);

      /* check non-blocking send has completed. */
      ierr = MPI_Wait(&(send_rq[0]), MPI_STATUS_IGNORE);

      /* wait on non-blocking receive from L and unpack data */
      ierr = MPI_Wait(&(recv_rq[0]), MPI_STATUS_IGNORE);
      unpack_ix3(vector, pMat);

    }
#endif /* MPI_PARALLEL */

/* Physical boundaries on both left and right */
    if (pMat->rx3_id < 0 && pMat->lx3_id < 0) {
      if(pMat->BCFlag_ix3 == 4)
		periodic_ix3(vector, pMat);
	else
		non_periodic_ix3(vector, pMat);

       if(pMat->BCFlag_ox3 == 4)
		periodic_ox3(vector, pMat);
	else
		non_periodic_ox3(vector, pMat);
    } 

  }



	/* Destroy the allocated memory */
	 bvals_matrix_vector_destruct(pMat);
#ifdef SHEARING_BOX
	bvals_matrix_vector_shear_destruct();
#endif

  return;
}



void bvals_matrix_vector_init(MatrixS *pMat)
{
/* MPI flag, lx1, rx1, rx2, lx2, rx3, lx3 are all set when */
/* Matrix object is created */ 

#ifdef MPI_PARALLEL
  int myL,myM,myN,l,m,n,nx1t,nx2t,nx3t,size;
  int x1cnt=0, x2cnt=0, x3cnt=0; /* Number of words passed in x1/x2/x3-dir. */
#endif /* MPI_PARALLEL */


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

void bvals_matrix_vector_destruct(MatrixS *pMat)
{

/* Only need to do this for MPI case to delete the buffar */
/* Because the matrix will be destroyed at each level. We need to clean this */

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


/*=========================== PRIVATE FUNCTIONS ==============================*/


/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions, Inner x1 boundary (bc_ix1=4) */

static void periodic_ix1(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		for(m=0; m<Nrad; m++){
       			vector[k][j][is-i][m] =  vector[k][j][ie-(i-1)][m];
		}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void periodic_ox1(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		for(m=0; m<Nrad; m++){
        		vector[k][j][ie+i][m] =  vector[k][j][is+(i-1)][m];
		}
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void periodic_ix2(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j, k, m;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[k][js-j][i][m] =  vector[k][je-(j-1)][i][m];
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void periodic_ox2(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
      		for(m=0; m<Nrad; m++){
			vector[k][je+j][i][m] =  vector[k][js+(j-1)][i][m];
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void periodic_ix3(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[ks-k][j][i][m] =  vector[ke-(k-1)][j][i][m];
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void periodic_ox3(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[ke+k][j][i][m] =  vector[ks+(k-1)][j][i][m];	
      		}
	}
    }
  }


  return;
}

/*---------------------------------------------------------------------*/
/*=======================================================================*/
/* non-periodic boundary condition. Just set it to be zero */

static void non_periodic_ix1(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		for(m=0; m<Nrad; m++){
       			vector[k][j][is-i][m] =  0.0;
		}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x1 boundary (bc_ox1=4) */

static void non_periodic_ox1(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;


  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=Matghost; i++) {
		for(m=0; m<Nrad; m++){
        		vector[k][j][ie+i][m] =  0.0;
		}
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x2 boundary (bc_ix2=4) */

static void non_periodic_ix2(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j, k, m;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[k][js-j][i][m] =  0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x2 boundary (bc_ox2=4) */

static void non_periodic_ox2(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
      		for(m=0; m<Nrad; m++){
			vector[k][je+j][i][m] =  0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Inner x3 boundary (bc_ix3=4) */

static void non_periodic_ix3(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[ks-k][j][i][m] =  0.0;
		}
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PERIODIC boundary conditions (cont), Outer x3 boundary (bc_ox3=4) */

static void non_periodic_ox3(Real ****vector, MatrixS *pMat)
{
  int is = pMat->is, ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,j,k, m;

  for (k=1; k<=Matghost; k++) {
    for (j=js-Matghost; j<=je+Matghost; j++) {
      for (i=is-Matghost; i<=ie+Matghost; i++) {
		for(m=0; m<Nrad; m++){
			vector[ke+k][j][i][m] =  0.0;	
      		}
	}
    }
  }


  return;
}




#ifdef MPI_PARALLEL  /* This ifdef wraps the next 12 funs; ~800 lines */
/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x1 boundary */

static void pack_ix1(Real ****vector, MatrixS *pMat)
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
        *(pSnd++) = vector[k][j][i][0];
        *(pSnd++) = vector[k][j][i][1];
        *(pSnd++) = vector[k][j][i][2];
        *(pSnd++) = vector[k][j][i][3];
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x1 boundary */

static void pack_ox1(Real ****vector, MatrixS *pMat)
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
		*(pSnd++) = vector[k][j][i][0];
        	*(pSnd++) = vector[k][j][i][1];
        	*(pSnd++) = vector[k][j][i][2];
        	*(pSnd++) = vector[k][j][i][3];
	}
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x2 boundary */

static void pack_ix2(Real ****vector, MatrixS *pMat)
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
		*(pSnd++) = vector[k][j][i][0];
        	*(pSnd++) = vector[k][j][i][1];
        	*(pSnd++) = vector[k][j][i][2];
        	*(pSnd++) = vector[k][j][i][3];		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x2 boundary */

static void pack_ox2(Real ****vector, MatrixS *pMat)
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
		*(pSnd++) = vector[k][j][i][0];
        	*(pSnd++) = vector[k][j][i][1];
        	*(pSnd++) = vector[k][j][i][2];
        	*(pSnd++) = vector[k][j][i][3];
	}
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Inner x3 boundary */

static void pack_ix3(Real ****vector, MatrixS *pMat)
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
     		*(pSnd++) = vector[k][j][i][0];
        	*(pSnd++) = vector[k][j][i][1];
        	*(pSnd++) = vector[k][j][i][2];
        	*(pSnd++) = vector[k][j][i][3];		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* PACK boundary conditions for MPI_Isend, Outer x3 boundary */

static void pack_ox3(Real ****vector, MatrixS *pMat)
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
     		*(pSnd++) = vector[k][j][i][0];
        	*(pSnd++) = vector[k][j][i][1];
        	*(pSnd++) = vector[k][j][i][2];
        	*(pSnd++) = vector[k][j][i][3];
		
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x1 boundary */

static void unpack_ix1(Real ****vector, MatrixS *pMat)
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
        vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x1 boundary */

static void unpack_ox1(Real ****vector, MatrixS *pMat)
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
 	vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);
	

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x2 boundary */

static void unpack_ix2(Real ****vector, MatrixS *pMat)
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
	vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x2 boundary */

static void unpack_ox2(Real ****vector, MatrixS *pMat)
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
	vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);
	
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Inner x3 boundary */

static void unpack_ix3(Real ****vector, MatrixS *pMat)
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
	vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);

      }
    }
  }


  return;
}

/*----------------------------------------------------------------------------*/
/* UNPACK boundary conditions after MPI_Irecv, Outer x3 boundary */

static void unpack_ox3(Real ****vector, MatrixS *pMat)
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
  	vector[k][j][i][0]  	= *(pRcv++);
        vector[k][j][i][1] 	= *(pRcv++);
        vector[k][j][i][2] 	= *(pRcv++);
        vector[k][j][i][3] 	= *(pRcv++);
	
      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */

#endif /* MATRIX_MULTIGRID */

#endif /* End radiation hydro and radiation MHD */



