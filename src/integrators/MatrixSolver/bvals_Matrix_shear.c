#include "../../copyright.h"
/*==============================================================================
 * FILE: bvals_Matrix_shear.c
 *
 * PURPOSE: Shearing sheet boundary conditions at ix1 and ox1
 *
 * Configure code with
 *   --enable-shearing-box to use.
 *

 *
 * CONTAINS PUBLIC FUNCTIONS:

 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../../prototypes.h"

/* The functions in this file will only work with the shearing box */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)

#if defined(MATRIX_MULTIGRID) || defined(MATRIX_HYPRE)

#ifdef SHEARING_BOX

/* Define number of variables to be remapped */
/* We need to remap Er, Fr1, Fr2, Fr3, V1, V2, V3, T4, Edd11, Edd21, Edd22, Edd31, Edd32 */
/* Edd33, sigma_t, sigma_a */
#define NREMAP (4*Matghost)



/* Define structure which holds variables remapped by shearing sheet BCs */
typedef struct Remap_s{
  Real U[NREMAP];
}Remap;

/*
  Remap.U[0] = Er;
  Remap.U[1] = Fr1;
  Remap.U[2] = Fr2;
  Remap.U[3] = Fr3;
*/


/* The memory for all the arrays below is allocated in bvals_shear_init */
/* Arrays of ghost zones containing remapped conserved quantities */
static Remap ***GhstZns=NULL, ***GhstZnsBuf=NULL;
/* 1D vectors for reconstruction in conservative remap step */


/* MPI send and receive buffers */
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif

/*==============================================================================
 * RADIATION quantieis are only 1st order accuracy
 *============================================================================*/


#endif

#ifdef SHEARING_BOX
/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ShearingSheet_radMHD_ix1: 3D shearing-sheet BCs in x1.  It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in bvals_mhd.c
 * NOTICE that the periodic copy has been applied before this function is called.
 *
 * This is a public function which is called by bvals_radMHD() inside a
 * SHEARING_BOX macro.
 *----------------------------------------------------------------------------*/
void ShearingSheet_Matrix_ix1(MatrixS *pMat)
{
  int is = pMat->is;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real Lx,Ly,qomL,yshear,deltay, dFry;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  RadMHDS *pCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  Lx = pMat->Lx;

  Ly = pMat->Ly;

  qomL = qshear*Omega_0*Lx;
/* When solve the implicit matrix, need to add dt in */
/* Add a factor of 4.0/3.0 */
/*
  yshear = qomL*(pMat->time+pMat->dt) * (1.0 + 1.0/3.0);
*/
  yshear = qomL*pMat->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pMat->Nx2 > joffset > pMat->Nx2   */
/* For radiation quantities, we only need the integer part */

  joffset = (int)(deltay/pMat->dx2);

  if((deltay - joffset * pMat->dx2) > 0.5 * pMat->dx2)
    joffset += 1;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * switch i and j just to speed up the code in the cyclic step
 * Steps 2-10 are for 3D or 2d xy runs.  Step 11 handles 2D xz separately */

  if (pMat->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

    if (pMat->Nx[2] > 1) ku=ke+1; else ku=ke;
    for(k=ks; k<=ku; k++) {
      for(j=js-Matghost; j<=je+Matghost; j++){
        for(i=0; i<Matghost; i++){
          ii = is-Matghost+i;
          GhstZns[k][i][j].U[0] = pMat->U[k][j][ii].Er;
          GhstZns[k][i][j].U[1] = pMat->U[k][j][ii].Fr1;

          GhstZns[k][i][j].U[2] = pMat->U[k][j][ii].Fr2;

/*      Multigrid solves for the residule, background shearing is subtracted */

/*      if(!(pMat->bgflag) && (pMat->Nx[2] == pMat->RootNx[2])){
 */
/* We still need this even if the background state is subtracted */
          dFry = qomL * pMat->U[k][j][ii].Er * (1.0 + pMat->Ugas[k][j][ii].Edd_22)/Crat;
          GhstZns[k][i][j].U[2] += dFry;


/*      }
 */

          GhstZns[k][i][j].U[3] = pMat->U[k][j][ii].Fr3;

        }
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

    for(k=ks; k<=ku; k++) {
      for(i=0; i<Matghost; i++){

        for (n=0; n<(NREMAP); n++) {

          for(j=js-Matghost; j<=je+Matghost; j++){
            GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n];
          }
        }

      }
    }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

    if (pMat->NGrid[1] == 1) {

      for(k=ks; k<=ku; k++) {
        for(j=js; j<=je; j++){
          jremap = j - joffset;
          if (jremap < (int)js) jremap += pMat->Nx[1];

          for(i=0; i<Matghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }

          }

        }
      }

#ifdef MPI_PARALLEL
    } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

      my_iproc = pMat->my_iproc;
      my_jproc = pMat->my_jproc;
      my_kproc = pMat->my_kproc;

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
      Ngrids = (int)(joffset/pMat->Nx[1]);
      joverlap = joffset - Ngrids*pMat->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0 */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

      if (joverlap != 0) {

        jproc = my_jproc + (Ngrids + 1);
        if (jproc > (pMat->NGrid[1]-1)) jproc -= pMat->NGrid[1];

/*pMat->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain; */
/* ID of my_k, j, my_i is */
        sendto_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;


        jproc = my_jproc - (Ngrids + 1);
        if (jproc < 0) jproc += pMat->NGrid[1];
/*       getfrom_id = pMat->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain; */
        getfrom_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZnsBuf */

        cnt = Matghost*joverlap*(ku-ks+1)*(NREMAP);
/* Post a non-blocking receive for the input data */
        ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                         shearing_sheet_ix1_tag, pMat->Comm_Domain, &rq);

        pSnd = send_buf;
        for (k=ks; k<=ku; k++) {
          for (j=je-(joverlap-1); j<=je; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZnsBuf[k][i][j]);

              for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
            }
          }
        }
        ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                        shearing_sheet_ix1_tag, pMat->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

        ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

        pRcv = recv_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js; j<=js+(joverlap-1); j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZns[k][i][j]);

              for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);

            }
          }
        }

      }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

      if (Ngrids == 0) {

        for(k=ks; k<=ku; k++) {
          for(j=js+joverlap; j<=je; j++){
            jremap = j-joverlap;
            for(i=0; i<Matghost; i++){
              for (n=0; n<(NREMAP); n++) {
                GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
              }

            }
          }
        }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZnsBuf (this step replaces 5d) */

      } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

        jproc = my_jproc + Ngrids;
        if (jproc > (pMat->NGrid[1]-1)) jproc -= pMat->NGrid[1];

        sendto_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;

        jproc = my_jproc - Ngrids;
        if (jproc < 0) jproc += pMat->NGrid[1];
        getfrom_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;

        cnt = Matghost*(pMat->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP);
/* Post a non-blocking receive for the input data from the left grid */
        ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                         shearing_sheet_ix1_tag, pMat->Comm_Domain, &rq);

        pSnd = send_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js; j<=je-joverlap; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZnsBuf[k][i][j]);
              for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

            }
          }
        }
        ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                        shearing_sheet_ix1_tag, pMat->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */

        ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

        pRcv = recv_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js+joverlap; j<=je; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZns[k][i][j]);
              for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);

            }
          }
        }
      } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
    } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

    for(k=ks; k<=ke; k++) {
      for(j=js; j<=je; j++){
        for(i=0; i<Matghost; i++){

          pMat->U[k][j][is-Matghost+i].Er  = GhstZns[k][i][j].U[0];
          pMat->U[k][j][is-Matghost+i].Fr1 = GhstZns[k][i][j].U[1];
          pMat->U[k][j][is-Matghost+i].Fr2 = GhstZns[k][i][j].U[2];
          pMat->U[k][j][is-Matghost+i].Fr3 = GhstZns[k][i][j].U[3];


        }
      }
    }



/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

/* Here we handle the corners */

    if (pMat->NGrid[1] == 1) {

      for(k=ks; k<=ke; k++) {
        for(j=1; j<=Matghost; j++){
          for(i=is-Matghost; i<is; i++){
            pMat->U[k][js-j][i] = pMat->U[k][je-(j-1)][i];
            pMat->U[k][je+j][i] = pMat->U[k][js+(j-1)][i];

          }
        }
      }


#ifdef MPI_PARALLEL
    } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
      cnt = Matghost*Matghost*(ku-ks+1)*(NREMAP);
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMat->lx2_id,
                       shearing_sheet_ix1_tag, pMat->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++){
        for (j=je-Matghost+1; j<=je; j++){
          for (i=is-Matghost; i<is; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            *(pSnd++) = pCons->Er;
            *(pSnd++) = pCons->Fr1;
            *(pSnd++) = pCons->Fr2;
            *(pSnd++) = pCons->Fr3;

          }
        }
      }

/* send contents of buffer to the neighboring grid on R-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMat->rx2_id,
                      shearing_sheet_ix1_tag, pMat->Comm_Domain);

/* Wait to receive the input data from the left grid */
      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++){
        for (j=js-Matghost; j<=js-1; j++){
          for (i=is-Matghost; i<is; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            pCons->Er  = *(pRcv++);
            pCons->Fr1 = *(pRcv++);
            pCons->Fr2 = *(pRcv++);
            pCons->Fr3 = *(pRcv++);

          }
        }
      }

/* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMat->rx2_id,
                       shearing_sheet_ix1_tag, pMat->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++){
        for (j=js; j<=js+Matghost-1; j++){
          for (i=is-Matghost; i<is; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            *(pSnd++) = pCons->Er;
            *(pSnd++) = pCons->Fr1;
            *(pSnd++) = pCons->Fr2;
            *(pSnd++) = pCons->Fr3;

          }
        }
      }

/* send contents of buffer to the neighboring grid on L-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMat->lx2_id,
                      shearing_sheet_ix1_tag, pMat->Comm_Domain);

/* Wait to receive the input data from the left grid */
      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++){
        for (j=je+1; j<=je+Matghost; j++){
          for (i=is-Matghost; i<is; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            pCons->Er  = *(pRcv++);
            pCons->Fr1 = *(pRcv++);
            pCons->Fr2 = *(pRcv++);
            pCons->Fr3 = *(pRcv++);

          }
        }
      }
#endif /* MPI_PARALLEL */

    } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */

  } /* end of if */

/*--- Step 11 ------------------------------------------------------------------
 * Shearing sheet BC in 2D xz.  Periodic BC already applied in x1 and x2 in
 * bvals_mhd (including for MPI parallel jobs).  Now just have to add offset
 * to azimuthal velocity when FARGO not defined */

/* Assume shearing sheet is always in xy plane for radiaiton MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/* ShearingSheet_ox1: 3D shearing-sheet BCs in x1.  It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in bvals_mhd.c
 *
 * This is a public function which is called by bvals_mhd() inside a
 * SHEARING_BOX macro.
 *----------------------------------------------------------------------------*/

void ShearingSheet_Matrix_ox1(MatrixS *pMat)
{
  int ie = pMat->ie;
  int js = pMat->js, je = pMat->je;
  int ks = pMat->ks, ke = pMat->ke;
  int i,ii,j,k,ku,n,joffset,jremap;
  Real Lx,Ly,qomL,yshear,deltay, dFry;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  RadMHDS *pCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  Lx = pMat->Lx;

  Ly = pMat->Ly;

  qomL = qshear*Omega_0*Lx;
/* add a factor of dt and 4/3 */
/*
  yshear = qomL*(pMat->time+pMat->dt) * (1.0 + 1.0/3.0);
*/
  yshear = qomL*pMat->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pMat->Nx2 > joffset > pMat->Nx2   */

  joffset = (int)(deltay/pMat->dx2);

  if((deltay - joffset * pMat->dx2) > 0.5 * pMat->dx2)
    joffset += 1;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D or 2D xy runs.  Step 11 handles 2D xz separately */

  if (pMat->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

    if (pMat->Nx[2] > 1) ku=ke+1; else ku=ke;
    for(k=ks; k<=ku; k++) {
      for(j=js-Matghost; j<=je+Matghost; j++){
        for(i=0; i<Matghost; i++){
          ii = ie+1+i;
          GhstZns[k][i][j].U[0] = pMat->U[k][j][ii].Er;
          GhstZns[k][i][j].U[1] = pMat->U[k][j][ii].Fr1;

          GhstZns[k][i][j].U[2] = pMat->U[k][j][ii].Fr2;


/*      if(!(pMat->bgflag) && (pMat->Nx[2] == pMat->RootNx[2])){
 */
/* Still need this even background state is subtracted */
          dFry = qomL * pMat->U[k][j][ii].Er * (1.0 + pMat->Ugas[k][j][ii].Edd_22)/Crat;
          GhstZns[k][i][j].U[2] -= dFry;
/*      }
 */

          GhstZns[k][i][j].U[3] = pMat->U[k][j][ii].Fr3;

        }
      }
    }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

    for(k=ks; k<=ku; k++) {
      for(i=0; i<Matghost; i++){

        for (n=0; n<(NREMAP); n++) {

          for(j=js; j<=je; j++){
            GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n];
          }
        }

      }
    }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

    if (pMat->NGrid[1] == 1) {

      for(k=ks; k<=ku; k++) {
        for(j=js; j<=je; j++){
          jremap = j + joffset;
          if (jremap > (int)je) jremap -= pMat->Nx[1];

          for(i=0; i<Matghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }

          }

        }
      }

#ifdef MPI_PARALLEL
    } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

      my_iproc = pMat->my_iproc;
      my_jproc = pMat->my_jproc;
      my_kproc = pMat->my_kproc;

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
      Ngrids = (int)(joffset/pMat->Nx[1]);
      joverlap = joffset - Ngrids*pMat->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(overlap-1):je] is received from.  Only execute if joverlap>0  */
/* This can result in send/receive to self -- we rely on MPI to handle this
 * properly */

      if (joverlap != 0) {

        jproc = my_jproc - (Ngrids + 1);
        if (jproc < 0) jproc += pMat->NGrid[1];
        sendto_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;;

        jproc = my_jproc + (Ngrids + 1);
        if (jproc > (pMat->NGrid[1]-1)) jproc -= pMat->NGrid[1];
        getfrom_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZnsBuf */

        cnt = Matghost*joverlap*(ku-ks+1)*(NREMAP);
/* Post a non-blocking receive for the input data */
        ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                         shearing_sheet_ox1_tag, pMat->Comm_Domain, &rq);

        pSnd = send_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js; j<=js+(joverlap-1); j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZnsBuf[k][i][j]);

              for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

            }
          }
        }
        ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                        shearing_sheet_ox1_tag, pMat->Comm_Domain);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

        ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

        pRcv = recv_buf;
        for (k=ks; k<=ku; k++) {
          for (j=je-(joverlap-1); j<=je; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZns[k][i][j]);

              for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);

            }
          }
        }

      }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

      if (Ngrids == 0) {

        for(k=ks; k<=ku; k++) {
          for(j=js; j<=je-joverlap; j++){
            jremap = j+joverlap;
            for(i=0; i<Matghost; i++){
              for (n=0; n<(NREMAP); n++) {
                GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
              }

            }
          }
        }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+joverlap:je]
 * from GhstZnsBuf (this step replaces 5d) */

      } else {

/* index of sendto and getfrom processors in GData are -/+1 from Step 5a */

        jproc = my_jproc - Ngrids;
        if (jproc < 0) jproc += pMat->NGrid[1];
        sendto_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;;

        jproc = my_jproc + Ngrids;
        if (jproc > (pMat->NGrid[1]-1)) jproc -= pMat->NGrid[1];
        getfrom_id = my_kproc * pMat->NGrid[0] * pMat->NGrid[1] + jproc * pMat->NGrid[0] + my_iproc;

        cnt = Matghost*(pMat->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP);
/* Post a non-blocking receive for the input data from the left grid */
        ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                         shearing_sheet_ox1_tag, pMat->Comm_Domain, &rq);

        pSnd = send_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js+joverlap; j<=je; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZnsBuf[k][i][j]);
              for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

            }
          }
        }
        ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                        shearing_sheet_ox1_tag, pMat->Comm_Domain);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

        ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

        pRcv = recv_buf;
        for (k=ks; k<=ku; k++) {
          for (j=js; j<=je-joverlap; j++) {
            for(i=0; i<Matghost; i++){
/* Get a pointer to the Remap structure */
              pRemap = &(GhstZns[k][i][j]);
              for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pRcv++);

            }
          }
        }
      } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
    } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells, except B1i[ie+1] */

    for(k=ks; k<=ke; k++) {
      for(j=js; j<=je; j++){
        for(i=0; i<Matghost; i++){
          pMat->U[k][j][ie+1+i].Er  = GhstZns[k][i][j].U[0];
          pMat->U[k][j][ie+1+i].Fr1 = GhstZns[k][i][j].U[1];
          pMat->U[k][j][ie+1+i].Fr2 = GhstZns[k][i][j].U[2];
          pMat->U[k][j][ie+1+i].Fr3 = GhstZns[k][i][j].U[3];

        }
      }
    }



/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B as average of remapped face centered B, except B1.
 * The value of B2c at j=je is incorrect since B2i[je+1] not yet set -- fix in
 * step 10 below */


/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

    if (pMat->NGrid[1] == 1) {

      for(k=ks; k<=ke; k++) {
        for(j=1; j<=Matghost; j++){
          for(i=ie+1; i<=ie+Matghost; i++){
            pMat->U[k][js-j][i] = pMat->U[k][je-(j-1)][i];
            pMat->U[k][je+j][i] = pMat->U[k][js+(j-1)][i];

          }
        }
      }


#ifdef MPI_PARALLEL
    } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
      cnt = Matghost*Matghost*(ku-ks+1)*(NREMAP);
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMat->lx2_id,
                       shearing_sheet_ox1_tag, pMat->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++){
        for (j=je-Matghost+1; j<=je; j++){
          for (i=ie+1; i<=ie+Matghost; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            *(pSnd++) = pCons->Er;
            *(pSnd++) = pCons->Fr1;
            *(pSnd++) = pCons->Fr2;
            *(pSnd++) = pCons->Fr3;

          }
        }
      }

/* send contents of buffer to the neighboring grid on R-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMat->rx2_id,
                      shearing_sheet_ox1_tag, pMat->Comm_Domain);

/* Wait to receive the input data from the left grid */
      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++){
        for (j=js-Matghost; j<=js-1; j++){
          for (i=ie+1; i<=ie+Matghost; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            pCons->Er  = *(pRcv++);
            pCons->Fr1 = *(pRcv++);
            pCons->Fr2 = *(pRcv++);
            pCons->Fr3 = *(pRcv++);
          }
        }
      }

/* Post a non-blocking receive for the input data from the right grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMat->rx2_id,
                       shearing_sheet_ox1_tag, pMat->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++){
        for (j=js; j<=js+Matghost-1; j++){
          for (i=ie+1; i<=ie+Matghost; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            *(pSnd++) = pCons->Er;
            *(pSnd++) = pCons->Fr1;
            *(pSnd++) = pCons->Fr2;
            *(pSnd++) = pCons->Fr3;

          }
        }
      }

/* send contents of buffer to the neighboring grid on L-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMat->lx2_id,
                      shearing_sheet_ox1_tag, pMat->Comm_Domain);

/* Wait to receive the input data from the left grid */
      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++){
        for (j=je+1; j<=je+Matghost; j++){
          for (i=ie+1; i<=ie+Matghost; i++){
/* Get a pointer to the ConsS cell */
            pCons = &(pMat->U[k][j][i]);

            pCons->Er  = *(pRcv++);
            pCons->Fr1 = *(pRcv++);
            pCons->Fr2 = *(pRcv++);
            pCons->Fr3 = *(pRcv++);
          }
        }
      }
#endif /* MPI_PARALLEL */

    } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */


  } /* end of if */

/*--- Step 11 ------------------------------------------------------------------
 * Shearing sheet BC in 2D.  Periodic BC already applied in x1 and x2 in
 * bvals_mhd (including for MPI parallel jobs).  Now just have to add offset
 * to azimuthal velocity when FARGO not defined */

  return;
}




/*----------------------------------------------------------------------------*/
/* bvals_shear_init: allocates memory for temporary arrays/buffers
 */

void bvals_Matrix_shear_init(MatrixS *pMat)
{

  int nx1,nx2,nx3,max1=0,max2=0,max3=0;
#ifdef MPI_PARALLEL
  int size1=0,size2=0,size;
#endif

/* Loop over all Grids on this processor to find maximum size of arrays */

/* set pointer to Grid */

  nx1 = pMat->Nx[0] + 2*Matghost;
  nx2 = pMat->Nx[1] + 2*Matghost;
  nx3 = pMat->Nx[2] + 2*Matghost;
  max1 = nx1;
  max2 = nx2;
  max3 = nx3;


/* Allocate memory for temporary arrays and vectors */

  if((GhstZns=(Remap***)calloc_3d_array(max3,Matghost,max2,sizeof(Remap)))==NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsBuf=(Remap***)calloc_3d_array(max3,Matghost,max2,sizeof(Remap))) ==
     NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");



/* estimate extra ghost cells needed by FARGO, allocate arrays accordingly */



/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size1 = Matghost*pMat->Nx[1]*(pMat->Nx[2]+1)*(NREMAP);

  size = size1;

  if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate send buffer\n");

  if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate receive buffer\n");
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* bvals_shear_destruct:  Free temporary arrays
 */

void bvals_Matrix_shear_destruct(void)
{
  if (GhstZns    != NULL) free_3d_array(GhstZns);
  if (GhstZnsBuf != NULL) free_3d_array(GhstZnsBuf);


#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

  return;
}


#endif /* SHEARING_BOX */

#endif /* MATRIX_MULTIGRID */
#endif /* END RADIATION HYDRO OR RADIATION MDH */
