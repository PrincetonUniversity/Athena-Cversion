#include "copyright.h"
/*==============================================================================
 * FILE: bvals_radMHD_shear.c
 *
 * PURPOSE: Shearing sheet boundary conditions at ix1 and ox1 for both 2D and 3D
 *   Called by bvals_radMHD.  Set shearing sheet boundary condition for 
 *
 * Configure code with
 *   --enable-shearing-box to use.
 *

 *
 * CONTAINS PUBLIC FUNCTIONS:
 * ShearingSheet_radMHD_ix1() - shearing sheet BCs on ix1
 * ShearingSheet_radMHD_ox1() - shearing sheet BCs on ox1
 * bvals_shear_radMHD_init() - allocates memory for arrays used here
 * bvals_shear_radMHD_destruct() - frees memory for arrays used here
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* The functions in this file will only work with the shearing box */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef SHEARING_BOX

/* Define number of variables to be remapped */
/* We only need to remap radiation quantities Er, Fr1, Fr2, Fr3 here */
#define NREMAP (10+NOPACITY)
#define NVAR_SHARE NVAR



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
void ShearingSheet_radMHD_ix1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap,m;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay, dFry;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  ConsS *pCons;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */
/* For radiation quantities, we only need the integer part */
	
  joffset = (int)(deltay/pG->dx2);
	
	if((deltay - joffset * pG->dx2) > 0.5 * pG->dx2)
		joffset += 1;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * switch i and j just to speed up the code in the cyclic step 
 * Steps 2-10 are for 3D or 2d xy runs.  Step 11 handles 2D xz separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = is-nghost+i;
        GhstZns[k][i][j].U[0] = pG->U[k][j][ii].Er;
        GhstZns[k][i][j].U[1] = pG->U[k][j][ii].Fr1;
        
		GhstZns[k][i][j].U[2] = pG->U[k][j][ii].Fr2;
		dFry = qomL * pG->U[k][j][ii].Er * (1.0 + pG->U[k][j][ii].Edd_22)/Crat;
		GhstZns[k][i][j].U[2] += dFry;
 
		  
		GhstZns[k][i][j].U[3] = pG->U[k][j][ii].Fr3;
		GhstZns[k][i][j].U[4] = pG->U[k][j][ii].Edd_11;
		GhstZns[k][i][j].U[5] = pG->U[k][j][ii].Edd_21;
		GhstZns[k][i][j].U[6] = pG->U[k][j][ii].Edd_22;
		GhstZns[k][i][j].U[7] = pG->U[k][j][ii].Edd_31;
		GhstZns[k][i][j].U[8] = pG->U[k][j][ii].Edd_32;
		GhstZns[k][i][j].U[9] = pG->U[k][j][ii].Edd_33;
		  for(m=0;m<NOPACITY;m++){
			  GhstZns[k][i][j].U[10+m] = pG->U[k][j][ii].Sigma[m];
			  
		  }
		  
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){

      for (n=0; n<(NREMAP); n++) {

        for(j=js-nghost; j<=je+nghost; j++){
          GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n];
        }
      }

    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx[1];

        for(i=0; i<nghost; i++){
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
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];
          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
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
          for(i=0; i<nghost; i++){
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
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
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
      for(i=0; i<nghost; i++){
        pG->U[k][j][is-nghost+i].Er  = GhstZns[k][i][j].U[0];
        pG->U[k][j][is-nghost+i].Fr1 = GhstZns[k][i][j].U[1];
        pG->U[k][j][is-nghost+i].Fr2 = GhstZns[k][i][j].U[2];
        pG->U[k][j][is-nghost+i].Fr3 = GhstZns[k][i][j].U[3];
		pG->U[k][j][is-nghost+i].Edd_11 = GhstZns[k][i][j].U[4];
		pG->U[k][j][is-nghost+i].Edd_21 = GhstZns[k][i][j].U[5];
		pG->U[k][j][is-nghost+i].Edd_22 = GhstZns[k][i][j].U[6];
		pG->U[k][j][is-nghost+i].Edd_31 = GhstZns[k][i][j].U[7];
		pG->U[k][j][is-nghost+i].Edd_32 = GhstZns[k][i][j].U[8];
		pG->U[k][j][is-nghost+i].Edd_33 = GhstZns[k][i][j].U[9];
		  for(m=0;m<NOPACITY;m++){
			 pG->U[k][j][is-nghost+i].Sigma[m] = GhstZns[k][i][j].U[10+m];
			  
		  }  
		  
	  }
    }
  }



/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */
	  
	  /* Here we handle the corners */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=is-nghost; i<is; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];

        }
      }
    }


#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks+1)*(NREMAP);
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
			pCons = &(pG->U[k][j][i]);

			*(pSnd++) = pCons->Er;
			*(pSnd++) = pCons->Fr1;
			*(pSnd++) = pCons->Fr2;
			*(pSnd++) = pCons->Fr3;
			*(pSnd++) = pCons->Edd_11;
			*(pSnd++) = pCons->Edd_21;
			*(pSnd++) = pCons->Edd_22;
			*(pSnd++) = pCons->Edd_31;
			*(pSnd++) = pCons->Edd_32;
			*(pSnd++) = pCons->Edd_33;
			for(m=0;m<NOPACITY;m++){
				*(pSnd++) = pCons->Sigma[m];
				
			}
			

        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
			pCons = &(pG->U[k][j][i]);

			pCons->Er  = *(pRcv++);
			pCons->Fr1 = *(pRcv++);
			pCons->Fr2 = *(pRcv++);
			pCons->Fr3 = *(pRcv++);
			pCons->Edd_11  = *(pRcv++);
			pCons->Edd_21 = *(pRcv++);
			pCons->Edd_22 = *(pRcv++);
			pCons->Edd_31 = *(pRcv++);
			pCons->Edd_32  = *(pRcv++);
			pCons->Edd_33 = *(pRcv++);
			for(m=0;m<NOPACITY;m++){
				pCons->Sigma[m] = *(pRcv++);
				
			}
			
			
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

			*(pSnd++) = pCons->Er;
			*(pSnd++) = pCons->Fr1;
			*(pSnd++) = pCons->Fr2;
			*(pSnd++) = pCons->Fr3;
			*(pSnd++) = pCons->Edd_11;
			*(pSnd++) = pCons->Edd_21;
			*(pSnd++) = pCons->Edd_22;
			*(pSnd++) = pCons->Edd_31;
			*(pSnd++) = pCons->Edd_32;
			*(pSnd++) = pCons->Edd_33;
			for(m=0;m<NOPACITY;m++){
				*(pSnd++) = pCons->Sigma[m];
				
			}


        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_ix1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the ConsS cell */
			pCons = &(pG->U[k][j][i]);

			pCons->Er  = *(pRcv++);
			pCons->Fr1 = *(pRcv++);
			pCons->Fr2 = *(pRcv++);
			pCons->Fr3 = *(pRcv++);
			pCons->Edd_11  = *(pRcv++);
			pCons->Edd_21 = *(pRcv++);
			pCons->Edd_22 = *(pRcv++);
			pCons->Edd_31 = *(pRcv++);
			pCons->Edd_32  = *(pRcv++);
			pCons->Edd_33 = *(pRcv++);
			
			for(m=0;m<NOPACITY;m++){
				pCons->Sigma[m] = *(pRcv++);
				
			}
			
			
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

void ShearingSheet_radMHD_ox1(DomainS *pD)
{
  GridS *pG = pD->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,ku,n,joffset,jremap,m;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay, dFry;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
  Remap *pRemap;
  ConsS *pCons;
  MPI_Request rq;
#endif /* end MPI_PARALLEL */

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  qomL = qshear*Omega_0*Lx;
  yshear = qomL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
	
	if((deltay - joffset * pG->dx2) > 0.5 * pG->dx2)
		joffset += 1;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched.
 * Steps 2-10 are for 3D or 2D xy runs.  Step 11 handles 2D xz separately */

  if (pG->Nx[2] > 1 || ShBoxCoord==xy) {  /* this if ends at end of step 10 */

  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
		  ii = ie+1+i;
		  GhstZns[k][i][j].U[0] = pG->U[k][j][ii].Er;
		  GhstZns[k][i][j].U[1] = pG->U[k][j][ii].Fr1;
		  
		  GhstZns[k][i][j].U[2] = pG->U[k][j][ii].Fr2;
		  dFry = qomL * pG->U[k][j][ii].Er * (1.0 + pG->U[k][j][ii].Edd_22)/Crat;
		  GhstZns[k][i][j].U[2] -= dFry;
		  
		  GhstZns[k][i][j].U[3] = pG->U[k][j][ii].Fr3;
		  GhstZns[k][i][j].U[4] = pG->U[k][j][ii].Edd_11;
		  GhstZns[k][i][j].U[5] = pG->U[k][j][ii].Edd_21;
		  GhstZns[k][i][j].U[6] = pG->U[k][j][ii].Edd_22;
		  GhstZns[k][i][j].U[7] = pG->U[k][j][ii].Edd_31;
		  GhstZns[k][i][j].U[8] = pG->U[k][j][ii].Edd_32;
		  GhstZns[k][i][j].U[9] = pG->U[k][j][ii].Edd_33;
		  for(m=0;m<NOPACITY;m++){
			  GhstZns[k][i][j].U[10+m] = pG->U[k][j][ii].Sigma[m]; 
		  }
  
		  
		  
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ku; k++) {
    for(i=0; i<nghost; i++){

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

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx[1];

        for(i=0; i<nghost; i++){
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

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx[1]);
    joverlap = joffset - Ngrids*pG->Nx[1];

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(overlap-1):je] is received from.  Only execute if joverlap>0  */
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
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZnsBuf */

      cnt = nghost*joverlap*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
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
          for(i=0; i<nghost; i++){
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
      if (jproc < 0) jproc += pD->NGrid[1];
      sendto_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid[1]-1)) jproc -= pD->NGrid[1];
      getfrom_id = pD->GData[my_kproc][jproc][my_iproc].ID_Comm_Domain;

      cnt = nghost*(pG->Nx[1]-joverlap)*(ku-ks+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pSnd++) = pRemap->U[n];

          }
        }
      }
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for (k=ks; k<=ku; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
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
      for(i=0; i<nghost; i++){
		  pG->U[k][j][ie+1+i].Er  = GhstZns[k][i][j].U[0];
		  pG->U[k][j][ie+1+i].Fr1 = GhstZns[k][i][j].U[1];
		  pG->U[k][j][ie+1+i].Fr2 = GhstZns[k][i][j].U[2];
		  pG->U[k][j][ie+1+i].Fr3 = GhstZns[k][i][j].U[3];
		  pG->U[k][j][ie+1+i].Edd_11 = GhstZns[k][i][j].U[4];
		  pG->U[k][j][ie+1+i].Edd_21 = GhstZns[k][i][j].U[5];
		  pG->U[k][j][ie+1+i].Edd_22 = GhstZns[k][i][j].U[6];
		  pG->U[k][j][ie+1+i].Edd_31 = GhstZns[k][i][j].U[7];
		  pG->U[k][j][ie+1+i].Edd_32 = GhstZns[k][i][j].U[8];
		  pG->U[k][j][ie+1+i].Edd_33 = GhstZns[k][i][j].U[9];
		  for(m=0;m<NOPACITY;m++){
			  pG->U[k][j][ie+1+i].Sigma[m] =  GhstZns[k][i][j].U[10+m];			  
		  }

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

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=ie+1; i<=ie+nghost; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];

        }
      }
    }


#ifdef MPI_PARALLEL
  } else {

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */


/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(ku-ks+1)*(NREMAP);
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
			pCons = &(pG->U[k][j][i]);

			*(pSnd++) = pCons->Er;
			*(pSnd++) = pCons->Fr1;
			*(pSnd++) = pCons->Fr2;
			*(pSnd++) = pCons->Fr3;
			*(pSnd++) = pCons->Edd_11;
			*(pSnd++) = pCons->Edd_21;
			*(pSnd++) = pCons->Edd_22;
			*(pSnd++) = pCons->Edd_31;
			*(pSnd++) = pCons->Edd_32;
			*(pSnd++) = pCons->Edd_33;
			for(m=0;m<NOPACITY;m++){
				*(pSnd++) = pCons->Sigma[m];
			}

			
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   shearing_sheet_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

			pCons->Er  = *(pRcv++);
			pCons->Fr1 = *(pRcv++);
			pCons->Fr2 = *(pRcv++);
			pCons->Fr3 = *(pRcv++);
			pCons->Edd_11  = *(pRcv++);
			pCons->Edd_21 = *(pRcv++);
			pCons->Edd_22 = *(pRcv++);
			pCons->Edd_31 = *(pRcv++);
			pCons->Edd_32  = *(pRcv++);
			pCons->Edd_33 = *(pRcv++);
			for(m=0;m<NOPACITY;m++){
				pCons->Sigma[m] = *(pRcv++);
				
			}
			
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);
			
			*(pSnd++) = pCons->Er;
			*(pSnd++) = pCons->Fr1;
			*(pSnd++) = pCons->Fr2;
			*(pSnd++) = pCons->Fr3;
			*(pSnd++) = pCons->Edd_11;
			*(pSnd++) = pCons->Edd_21;
			*(pSnd++) = pCons->Edd_22;
			*(pSnd++) = pCons->Edd_31;
			*(pSnd++) = pCons->Edd_32;
			*(pSnd++) = pCons->Edd_33;
			for(m=0;m<NOPACITY;m++){
				*(pSnd++) = pCons->Sigma[m];
				
			}

        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   shearing_sheet_ox1_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the ConsS cell */
          pCons = &(pG->U[k][j][i]);

			pCons->Er  = *(pRcv++);
			pCons->Fr1 = *(pRcv++);
			pCons->Fr2 = *(pRcv++);
			pCons->Fr3 = *(pRcv++);
			pCons->Edd_11  = *(pRcv++);
			pCons->Edd_21 = *(pRcv++);
			pCons->Edd_22 = *(pRcv++);
			pCons->Edd_31 = *(pRcv++);
			pCons->Edd_32  = *(pRcv++);
			pCons->Edd_33 = *(pRcv++);
			for(m=0;m<NOPACITY;m++){
				pCons->Sigma[m] = *(pRcv++);
				
			}

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

void bvals_radMHD_shear_init(MeshS *pM)
{
  GridS *pG;
  int nl,nd,nx1,nx2,nx3,max1=0,max2=0,max3=0;
#ifdef MPI_PARALLEL
  int size1=0,size2=0,size;
#endif

/* Loop over all Grids on this processor to find maximum size of arrays */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this proc */
        pG=pM->Domain[nl][nd].Grid;          /* set pointer to Grid */

        nx1 = pG->Nx[0] + 2*nghost;
        nx2 = pG->Nx[1] + 2*nghost;
        nx3 = pG->Nx[2] + 2*nghost;
        max1 = MAX(max1,nx1);
        max2 = MAX(max2,nx2);
        max3 = MAX(max3,nx3);
      }
    }
  }

/* Allocate memory for temporary arrays and vectors */

  if((GhstZns=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap)))==NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsBuf=(Remap***)calloc_3d_array(max3,nghost,max2,sizeof(Remap))) ==
    NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");



/* estimate extra ghost cells needed by FARGO, allocate arrays accordingly */



/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size1 = nghost*pG->Nx[1]*(pG->Nx[2]+1)*(NREMAP+NSCALARS);

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

void bvals_radMHD_shear_destruct(void)
{
  if (GhstZns    != NULL) free_3d_array(GhstZns);
  if (GhstZnsBuf != NULL) free_3d_array(GhstZnsBuf);


#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

  return;
}


#endif /* END RADIATION HYDRO OR RADIATION MDH */
#endif /* SHEARING_BOX */
