#include "../copyright.h"
/*==============================================================================
 * FILE: bvals_shear.c
 *
 * PURPOSE: Shearing sheet boundary conditions for radiative transfer at ix1 
 *          and ox1 for both 2D and 3D.  Patterned after bvals_shear.c.  
 *          Called by bvals_rad. 
 *
 *          Configure code with  --enable-shearing-box to use.
 *
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * ShearingSheet_Rad_ix1() - shearing sheet BCs on ix1
 * ShearingSheet_Rad_ox1() - shearing sheet BCs on ox1
 *============================================================================*/


#include <stdlib.h>
#include <math.h>
#include <float.h>

#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER
#ifdef SHEARING_BOX

/* The memory for all the arrays below is allocated in bvals_shear_init */
/* Arrays of ghost zones containing remapped conserved quantities */
static Real ****GhstZnsMom=NULL, ****GhstZnsMomBuf=NULL;
static Real *****GhstZnsIntBuf=NULL;

/* 1D vectors for reconstruction in conservative remap step */
static Real *U=NULL, *Flx=NULL;
/* MPI send and receive buffers */
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif

static void RemapFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);

void ShearingSheet_Rad_ix1(DomainS *pD, int ifs, int ife)
{
  RadGridS *pRG = pD->RadGrid;
  Real time = pD->Grid->time;
  int il = pRG->is-1;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang, noct = pRG->noct;
  int i,j,k,l,m,ifr,joffset,jremap,nDim;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
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
  yshear = qomL*time;
  
  if (pRG->Nx[2] > 1) nDim = 3; else nDim=2;
/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

  //printf("shear_ix1: %d %g %g %g\n",myID_Comm_world,yshear,qomL,deltay);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pRG->dx2);
  epsi = (fmod(deltay,pRG->dx2))/pRG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array. */

/* Copy moments to temporary arrays */
  for(k=ks; k<=ke; k++) {
    for(j=js-1; j<=je+1; j++){
      for (ifr=ifs; ifr<=ife; ifr++) {
	GhstZnsMom[ifr][k][j][0] = pRG->R[ifr][k][j][il].S;
	GhstZnsMom[ifr][k][j][1] = pRG->R[ifr][k][j][il].J;
	for(l=0; l<nDim; l++) {
	  GhstZnsMom[ifr][k][j][l+2] = pRG->R[ifr][k][j][il].H[l];
	}
      }}}

/* intensities at all angles are copied into pRG->Ghstl1i array in unpack_ix1 fucntion */

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke; k++) {
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(l=0; l<nDim+2; l++) {
	for (j=js-1; j<=je+1; j++) U[j] =  GhstZnsMom[ifr][k][j][l];
	RemapFlux(U,epsi,js,je+1,Flx);
	for(j=js; j<=je; j++){
	  GhstZnsMomBuf[ifr][k][j][l] = GhstZnsMom[ifr][k][j][l] -
	    (Flx[j+1]-Flx[j]);
	}
      }}}
  for (ifr=ifs; ifr<=ife; ifr++) {
    for(k=ks; k<=ke; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  for (j=js-1; j<=je+1; j++) U[j] =  pRG->Ghstl1i[ifr][k][j][l][m];
	  RemapFlux(U,epsi,js,je+1,Flx);
	  for(j=js; j<=je; j++){
	    GhstZnsIntBuf[ifr][k][j][l][m] = pRG->Ghstl1i[ifr][k][j][l][m] -
	      (Flx[j+1]-Flx[j]);
	  }
	}}}}

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pRG->Nx[1];
	for (ifr=ifs; ifr<=ife; ifr++) {
	  for(l=0; l<nDim+2; l++) {
	    GhstZnsMom[ifr][k][j][l] = GhstZnsMomBuf[ifr][k][jremap][l];
	  }}}}

    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(j=js; j<=je; j++){
	  jremap = j - joffset;
	  if (jremap < (int)js) jremap += pRG->Nx[1];
	  for(l=0; l<noct; l++) {
	    for(m=0; m<nang; m++) {
	      pRG->Ghstl1i[ifr][k][j][l][m] = GhstZnsIntBuf[ifr][k][jremap][l][m];
	    }}}}}
#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
  
/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pRG->Nx[1]);
    joverlap = joffset - Ngrids*pRG->Nx[1];
 
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

      cnt = joverlap*(ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for(k=ks; k<=ke; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      (*pSnd++) = GhstZnsMomBuf[ifr][k][j][l];
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=je-(joverlap-1); j<=je; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		(*pSnd++) = GhstZnsIntBuf[ifr][k][j][l][m];
	      }}}}}

      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for(k=ks; k<=ke; k++) {
	for (j=js; j<=js+(joverlap-1); j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = *(pRcv++);
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js; j<=js+(joverlap-1); j++) {	  
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstl1i[ifr][k][j][l][m] = *(pRcv++);
	      }}}}}

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */
    if (Ngrids == 0) {

      for(k=ks; k<=ke; k++) {
	for(j=js+joverlap; j<=je; j++){
	  jremap = j-joverlap;
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = GhstZnsMomBuf[ifr][k][jremap][l];
	    }}}}
      
      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for(j=js+joverlap; j<=je; j++){
	    jremap = j-joverlap;
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstl1i[ifr][k][j][l][m] = GhstZnsIntBuf[ifr][k][jremap][l][m];
	      }}}}}

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

      cnt = (pRG->Nx[1]-joverlap)*(ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for(k=ks; k<=ke; k++) {
	for (j=js; j<=je-joverlap; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      (*pSnd++) = GhstZnsMomBuf[ifr][k][j][l];
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js; j<=je-joverlap; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		(*pSnd++) = GhstZnsIntBuf[ifr][k][j][l][m];
	      }}}}}

      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ix1_tag, pD->Comm_Domain);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for(k=ks; k<=ke; k++) {
	for (j=js+joverlap; j<=je; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = *(pRcv++);
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js+joverlap; j<=je; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstl1i[ifr][k][j][l][m] = *(pRcv++);
	      }}}}}
  
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

  for(k=ks; k<=ke; k++) {
    for(j=js-1; j<=je+1; j++){
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][j][il].S = GhstZnsMom[ifr][k][j][0];
	pRG->R[ifr][k][j][il].J = GhstZnsMom[ifr][k][j][1];
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][j][il].H[l] = GhstZnsMom[ifr][k][j][l+2];
	}
      }}}
/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][js-1][il].S = pRG->R[ifr][k][je][il].S;
	pRG->R[ifr][k][je+1][il].S = pRG->R[ifr][k][js][il].S;
	pRG->R[ifr][k][js-1][il].J = pRG->R[ifr][k][je][il].J;
	pRG->R[ifr][k][je+1][il].J = pRG->R[ifr][k][js][il].J;
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][js-1][il].H[l] = pRG->R[ifr][k][je][il].H[l];
	  pRG->R[ifr][k][je+1][il].H[l] = pRG->R[ifr][k][js][il].H[l];
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][js-1][l][m] = pRG->Ghstl1i[ifr][k][je][l][m];
	    pRG->Ghstl1i[ifr][k][je+1][l][m] = pRG->Ghstl1i[ifr][k][js][l][m];
	  }}}}
  
#ifdef MPI_PARALLEL
  } else {
/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_rad.c */

/* Post a non-blocking receive for the input data from the left grid */
    cnt = (ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2); 

    if (pRG->lx2_id != -1) {
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pRG->lx2_id,
		       shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);
    }      
    pSnd = send_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	*(pSnd++) = pRG->R[ifr][k][je][il].S;
	*(pSnd++) = pRG->R[ifr][k][je][il].J;
	for(l=0; l<nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][je][il].H[l];
	}
      }}      
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][je][l][m];
	  }}}}

    if (pRG->rx2_id != -1) {
/* send contents of buffer to the neighboring grid on R-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pRG->rx2_id,
		      shearing_sheet_ix1_tag, pD->Comm_Domain);
    }

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][js-1][il].S = *(pRcv++);
	pRG->R[ifr][k][js-1][il].J = *(pRcv++);
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][js-1][il].H[l] = *(pRcv++);
	}
      }}      
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][js-1][l][m] = *(pRcv++);
	  }}}}
    

/* Post a non-blocking receive for the input data from the right grid */
    if (pRG->rx2_id != -1) {
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pRG->rx2_id,
		       shearing_sheet_ix1_tag, pD->Comm_Domain, &rq);
    }

    pSnd = send_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	*(pSnd++) = pRG->R[ifr][k][js][il].S;
	*(pSnd++) = pRG->R[ifr][k][js][il].J;
	for(l=0; l<nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][js][il].H[l];
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstl1i[ifr][k][js][l][m];
	  }}}}

/* send contents of buffer to the neighboring grid on L-x2 */
    if (pRG->lx2_id != -1) {
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pRG->lx2_id,
		      shearing_sheet_ix1_tag, pD->Comm_Domain);
    }

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][je+1][il].S = *(pRcv++); 
	pRG->R[ifr][k][je+1][il].J = *(pRcv++); 
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][je+1][il].H[l] = *(pRcv++);
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][je+1][l][m] = *(pRcv++);
	  }}}}
  
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10. ------------------------------------------------------------------
/* Update l/r2imu, and l/r3imu on corners using the remapped Ghstl1i 
 * values. */
  for (ifr=ifs; ifr<=ife; ifr++) {
    for (k=ks; k<=ke; k++) {
      for (m=0; m<nang; m++) {
	pRG->l2imu[ifr][k][il][0][m] = pRG->Ghstl1i[ifr][k][js-1][0][m];
	pRG->l2imu[ifr][k][il][1][m] = pRG->Ghstl1i[ifr][k][js-1][1][m];
	pRG->r2imu[ifr][k][il][2][m] = pRG->Ghstl1i[ifr][k][je+1][2][m];
	pRG->r2imu[ifr][k][il][3][m] = pRG->Ghstl1i[ifr][k][je+1][3][m];
	if(noct == 8) {
	  pRG->l2imu[ifr][k][il][4][m] = pRG->Ghstl1i[ifr][k][js-1][4][m];
	  pRG->l2imu[ifr][k][il][5][m] = pRG->Ghstl1i[ifr][k][js-1][5][m];
	  pRG->r2imu[ifr][k][il][6][m] = pRG->Ghstl1i[ifr][k][je+1][6][m];
	  pRG->r2imu[ifr][k][il][7][m] = pRG->Ghstl1i[ifr][k][je+1][7][m];
	}
      }}}
/* For l/r3imu we reset the top/bottom zones using pRG->Ghstl1i.  The iteration
 * over j runs from js-1 to je+1 (differs from "normal" periodic boundary
 * routine) because this function is called after the "normal" radiation
 * boundary condition has already been applied to the x2 face. */
    if (noct == 8) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	for (j=js-1; j<=je+1; j++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][il][4][m] = pRG->Ghstl1i[ifr][ks][j][4][m];
	    pRG->l3imu[ifr][j][il][5][m] = pRG->Ghstl1i[ifr][ks][j][5][m];
	    pRG->l3imu[ifr][j][il][6][m] = pRG->Ghstl1i[ifr][ks][j][6][m];
	    pRG->l3imu[ifr][j][il][7][m] = pRG->Ghstl1i[ifr][ks][j][7][m];
	    pRG->r3imu[ifr][j][il][0][m] = pRG->Ghstl1i[ifr][ke][j][0][m];
	    pRG->r3imu[ifr][j][il][1][m] = pRG->Ghstl1i[ifr][ke][j][1][m];
	    pRG->r3imu[ifr][j][il][2][m] = pRG->Ghstl1i[ifr][ke][j][2][m];
	    pRG->r3imu[ifr][j][il][3][m] = pRG->Ghstl1i[ifr][ke][j][3][m];	   
	  }}}
    }

}

void ShearingSheet_Rad_ox1(DomainS *pD, int ifs, int ife)
{
  RadGridS *pRG = pD->RadGrid;
  Real time = pD->Grid->time;
  int iu = pRG->ie+1;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang, noct = pRG->noct;
  int i,j,k,l,m,ifr,joffset,jremap,nDim;
  Real xmin,xmax,Lx,Ly,qomL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int ierr,sendto_id,getfrom_id;
  double *pSnd,*pRcv;
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
  yshear = qomL*time;

  if (pRG->Nx[2] > 1) nDim = 3; else nDim=2;
/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);
  //printf("shear_ox1: %d %g %g %g\n",myID_Comm_world,yshear,qomL,deltay);
/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epso < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pRG->Nx2   */

  joffset = (int)(deltay/pRG->dx2);
  epso = -(fmod(deltay,pRG->dx2))/pRG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array. */

/* Copy moments to temporary arrays */
  for(k=ks; k<=ke; k++) {
    for(j=js-1; j<=je+1; j++){
      for (ifr=ifs; ifr<=ife; ifr++) {
	GhstZnsMom[ifr][k][j][0] = pRG->R[ifr][k][j][iu].S;
	GhstZnsMom[ifr][k][j][1] = pRG->R[ifr][k][j][iu].J;
	for(l=0; l<nDim; l++) {
	  GhstZnsMom[ifr][k][j][l+2] = pRG->R[ifr][k][j][iu].H[l];
	}
      }}}

/* Intensities are copied into pRG->Ghstr1i array in unpack_ox1 fucntion */


/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke; k++) {
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(l=0; l<nDim+2; l++) {
	for (j=js-1; j<=je+1; j++) U[j] =  GhstZnsMom[ifr][k][j][l];
	RemapFlux(U,epso,js,je+1,Flx);
	for(j=js; j<=je; j++){
	  GhstZnsMomBuf[ifr][k][j][l] = GhstZnsMom[ifr][k][j][l] -
	    (Flx[j+1]-Flx[j]);
	}
      }}}
  for (ifr=ifs; ifr<=ife; ifr++) {
    for(k=ks; k<=ke; k++) {
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  for (j=js-1; j<=je+1; j++) U[j] =  pRG->Ghstr1i[ifr][k][j][l][m];
	  RemapFlux(U,epso,js,je+1,Flx);
	  for(j=js; j<=je; j++){
	    GhstZnsIntBuf[ifr][k][j][l][m] = pRG->Ghstr1i[ifr][k][j][l][m] -
	      (Flx[j+1]-Flx[j]);
	  }
	}}}}

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=js; j<=je; j++){
	jremap = j + joffset;
        if (jremap > (int)je) jremap -= pRG->Nx[1];
	for (ifr=ifs; ifr<=ife; ifr++) {
	  for(l=0; l<nDim+2; l++) {
	    GhstZnsMom[ifr][k][j][l] = GhstZnsMomBuf[ifr][k][jremap][l];
	  }}}}

    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(j=js; j<=je; j++){
	  jremap = j + joffset;
	  if (jremap > (int)je) jremap -= pRG->Nx[1];
	  for(l=0; l<noct; l++) {
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr1i[ifr][k][j][l][m] = GhstZnsIntBuf[ifr][k][jremap][l][m];
	    }}}}}
#ifdef MPI_PARALLEL
  } else {

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
  
/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pRG->Nx[1]);
    joverlap = joffset - Ngrids*pRG->Nx[1];

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

      cnt = joverlap*(ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2);
/* Post a non-blocking receive for the input data */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for(k=ks; k<=ke; k++) {
	for (j=js; j<=js+(joverlap-1); j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      (*pSnd++) = GhstZnsMomBuf[ifr][k][j][l];
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js; j<=js+(joverlap-1); j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		(*pSnd++) = GhstZnsIntBuf[ifr][k][j][l][m];
	      }}}}}

      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for(k=ks; k<=ke; k++) {
	for (j=je-(joverlap-1); j<=je; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = *(pRcv++);
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=je-(joverlap-1); j<=je; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstr1i[ifr][k][j][l][m] = *(pRcv++);
	      }}}}}

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = GhstZnsMomBuf[ifr][k][jremap][l];
	    }}}}
      
      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for(j=js; j<=je-joverlap; j++){
	    jremap = j+joverlap;
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstr1i[ifr][k][j][l][m] = GhstZnsIntBuf[ifr][k][jremap][l][m];
	      }}}}}

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

      cnt = (pRG->Nx[1]-joverlap)*(ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2);
/* Post a non-blocking receive for the input data from the left grid */
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

      pSnd = send_buf;
      for(k=ks; k<=ke; k++) {
        for (j=js+joverlap; j<=je; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      (*pSnd++) = GhstZnsMomBuf[ifr][k][j][l];
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js; j<=je-joverlap; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		(*pSnd++) = GhstZnsIntBuf[ifr][k][j][l][m];
	      }}}}}

      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     shearing_sheet_ox1_tag, pD->Comm_Domain);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

      ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

      pRcv = recv_buf;
      for(k=ks; k<=ke; k++) {
	for (j=js; j<=je-joverlap; j++) {
	  for (ifr=ifs; ifr<=ife; ifr++) {
	    for(l=0; l<nDim+2; l++) {
	      GhstZnsMom[ifr][k][j][l] = *(pRcv++);
	    }}}}

      for (ifr=ifs; ifr<=ife; ifr++) {
	for(k=ks; k<=ke; k++) {
	  for (j=js; j<=je-joverlap; j++) {
	    for(l=0; l<noct; l++) {
	      for(m=0; m<nang; m++) {
		pRG->Ghstr1i[ifr][k][j][l][m] = *(pRcv++);
	      }}}}}
  
   } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now copy remapped variables back into ghost cells */

  for(k=ks; k<=ke; k++) {
    for(j=js-1; j<=je+1; j++){
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][j][iu].S = GhstZnsMom[ifr][k][j][0];
	pRG->R[ifr][k][j][iu].J = GhstZnsMom[ifr][k][j][1];
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][j][iu].H[l] = GhstZnsMom[ifr][k][j][l+2];
	}
      }}}

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in bvals_mhd.c) */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][js-1][iu].S = pRG->R[ifr][k][je][iu].S;
	pRG->R[ifr][k][je+1][iu].S = pRG->R[ifr][k][js][iu].S;
	pRG->R[ifr][k][js-1][iu].J = pRG->R[ifr][k][je][iu].J;
	pRG->R[ifr][k][je+1][iu].J = pRG->R[ifr][k][js][iu].J;
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][js-1][iu].H[l] = pRG->R[ifr][k][je][iu].H[l];
	  pRG->R[ifr][k][je+1][iu].H[l] = pRG->R[ifr][k][js][iu].H[l];
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][js-1][l][m] = pRG->Ghstr1i[ifr][k][je][l][m];
	    pRG->Ghstr1i[ifr][k][je+1][l][m] = pRG->Ghstr1i[ifr][k][js][l][m];
	  }}}}
 
#ifdef MPI_PARALLEL
  } else {
/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */

/* Post a non-blocking receive for the input data from the left grid */
    cnt = (ke-ks+1)*(ife-ifs+1)*(noct*nang+nDim+2); 

    if (pRG->lx2_id != -1) {
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pRG->lx2_id,
		       shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);

    }

    pSnd = send_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	*(pSnd++) = pRG->R[ifr][k][je][iu].S;
	*(pSnd++) = pRG->R[ifr][k][je][iu].J;
	for(l=0; l<nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][je][iu].H[l];
	}
      }}      
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][je][l][m];
	  }}}}
    

    if (pRG->rx2_id != -1) {
/* send contents of buffer to the neighboring grid on R-x2 */
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pRG->rx2_id,
		      shearing_sheet_ox1_tag, pD->Comm_Domain);
    }

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][js-1][iu].S = *(pRcv++);
	pRG->R[ifr][k][js-1][iu].J = *(pRcv++);
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][js-1][iu].H[l] = *(pRcv++);
	}
      }}      
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][js-1][l][m] = *(pRcv++);
	  }}}}

/* Post a non-blocking receive for the input data from the right grid */
    if (pRG->rx2_id != -1) {
      ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pRG->rx2_id,
		       shearing_sheet_ox1_tag, pD->Comm_Domain, &rq);
    }

    pSnd = send_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	*(pSnd++) = pRG->R[ifr][k][js][iu].S;
	*(pSnd++) = pRG->R[ifr][k][js][iu].J;
	for(l=0; l<nDim; l++) {
	  *(pSnd++) = pRG->R[ifr][k][js][iu].H[l];
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    *(pSnd++) = pRG->Ghstr1i[ifr][k][js][l][m];
	  }}}}

/* send contents of buffer to the neighboring grid on L-x2 */
    if (pRG->lx2_id != -1) {
      ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pRG->lx2_id,
		      shearing_sheet_ox1_tag, pD->Comm_Domain);
    }

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for(k=ks; k<=ke; k++) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	pRG->R[ifr][k][je+1][iu].S = *(pRcv++); 
	pRG->R[ifr][k][je+1][iu].J = *(pRcv++); 
	for(l=0; l<nDim; l++) {
	  pRG->R[ifr][k][je+1][iu].H[l] = *(pRcv++);
	}
      }}
    for (ifr=ifs; ifr<=ife; ifr++) {
      for(k=ks; k<=ke; k++) {
	for(l=0; l<noct; l++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][je+1][l][m] = *(pRcv++);
	  }}}}
  
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10. ------------------------------------------------------------------
 * Update l/r2imu, and l/r3imu on corners using the remapped Ghstr1i 
 * values. */
   for (ifr=ifs; ifr<=ife; ifr++) {
     for (k=ks; k<=ke; k++) {
       for (m=0; m<nang; m++) {
	 pRG->l2imu[ifr][k][iu][0][m] = pRG->Ghstr1i[ifr][k][js-1][0][m];
	 pRG->l2imu[ifr][k][iu][1][m] = pRG->Ghstr1i[ifr][k][js-1][1][m];
	 pRG->r2imu[ifr][k][iu][2][m] = pRG->Ghstr1i[ifr][k][je+1][2][m];
	 pRG->r2imu[ifr][k][iu][3][m] = pRG->Ghstr1i[ifr][k][je+1][3][m];
	 if(noct == 8) {
	   pRG->l2imu[ifr][k][iu][4][m] = pRG->Ghstr1i[ifr][k][js-1][4][m];
	   pRG->l2imu[ifr][k][iu][5][m] = pRG->Ghstr1i[ifr][k][js-1][5][m];
	   pRG->r2imu[ifr][k][iu][6][m] = pRG->Ghstr1i[ifr][k][je+1][6][m];
	   pRG->r2imu[ifr][k][iu][7][m] = pRG->Ghstr1i[ifr][k][je+1][7][m];
	 }
       }}}
/* For l/r3imu we reset the top/bottom zones using pRG->Ghstr1i.  The iteration
 * over j runs from js-1 to je+1 (differs from "normal" periodic boundary
 * routine) because this function is called after the "normal" radiation
 * boundary condition has already been applied to the x2 face. */
    if (noct == 8) {
      for (ifr=ifs; ifr<=ife; ifr++) {
	for (j=js-1; j<=je+1; j++) {
	  for (m=0; m<nang; m++) {
	    pRG->l3imu[ifr][j][iu][4][m] = pRG->Ghstr1i[ifr][ks][j][4][m];
	    pRG->l3imu[ifr][j][iu][5][m] = pRG->Ghstr1i[ifr][ks][j][5][m];
	    pRG->l3imu[ifr][j][iu][6][m] = pRG->Ghstr1i[ifr][ks][j][6][m];
	    pRG->l3imu[ifr][j][iu][7][m] = pRG->Ghstr1i[ifr][ks][j][7][m];
	    pRG->r3imu[ifr][j][iu][0][m] = pRG->Ghstr1i[ifr][ke][j][0][m];
	    pRG->r3imu[ifr][j][iu][1][m] = pRG->Ghstr1i[ifr][ke][j][1][m];
	    pRG->r3imu[ifr][j][iu][2][m] = pRG->Ghstr1i[ifr][ke][j][2][m];
	    pRG->r3imu[ifr][j][iu][3][m] = pRG->Ghstr1i[ifr][ke][j][3][m];	   
	  }}}
    }
}

/*----------------------------------------------------------------------------*/
/* bvals_shear_init: allocates memory for temporary arrays/buffers
 */

void bvals_rad_shear_init(MeshS *pM)
{
  RadGridS *pRG;
  int nl,nd,nx1,nx2,nx3,max1=0,max2=0,max3=0;
  int nDim,nf,nang,noct,maxf=0,maxa=0;
#ifdef MPI_PARALLEL
  int size;
#endif

  if (pM->Nx[2] > 1) {
    nDim = 3;
    noct = 8;
  } else {
    nDim = 2;
    noct = 4;
  }
/* Loop over all Grids on this processor to find maximum size of arrays */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL) { /* there is a Grid on this proc */
        pRG=pM->Domain[nl][nd].RadGrid;          /* set pointer to Grid */
	nang = pRG->nang;  /* allows for variation of nang and nf across grids */
	nf = pRG->nf;
        nx1 = pRG->Nx[0] + 2;
        nx2 = pRG->Nx[1] + 2;
        nx3 = pRG->Nx[2] + 2;
	maxf = MAX(maxf,nf);
	maxa = MAX(maxa,nang);
        max1 = MAX(max1,nx1);
        max2 = MAX(max2,nx2);
        max3 = MAX(max3,nx3);
      }
    }
  }

/* Allocate memory for temporary arrays and vectors */

  if((GhstZnsMom=(Real****)calloc_4d_array(maxf,max3,max2,nDim+2,sizeof(Real)))==NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsMomBuf=(Real****)calloc_4d_array(maxf,max3,max2,nDim+2,sizeof(Real))) ==
    NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsIntBuf=(Real*****)calloc_5d_array(maxf,max3,max2,noct,maxa,sizeof(Real))) ==
    NULL) ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((U = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

  if((Flx = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[bvals_shear_init]: malloc returned a NULL pointer\n");

/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size = max3*max2*nf*(noct*nang+nDim+2);

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

void bvals_rad_shear_destruct(void)
{

  if (GhstZnsMom    != NULL) free_4d_array(GhstZnsMom);
  if (GhstZnsMomBuf != NULL) free_4d_array(GhstZnsMomBuf);
  if (GhstZnsIntBuf != NULL) free_5d_array(GhstZnsIntBuf);

  if (U   != NULL) free(U);
  if (Flx != NULL) free(Flx);

#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

  return;
}

/*------------------------------------------------------------------------------
 * RemapFlux(): second order reconstruction for conservative remap.
 *   using piecewise linear reconstruction and min/mod limiters
 */

static void RemapFlux(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real dUc,dUl,dUr,dUm,lim_slope;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dUc = U[j+1] - U[j-1];
      dUl = U[j  ] - U[j-1];
      dUr = U[j+1] - U[j  ];

      dUm = 0.0;
      if (dUl*dUr > 0.0) {
        lim_slope = MIN(fabs(dUl),fabs(dUr));
        dUm = SIGN(dUc)*MIN(0.5*fabs(dUc),2.0*lim_slope);
      }
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      Flux[j+1] = eps*(U[j] + 0.5*(1.0 - eps)*dUm);
    } else {         /* eps always < 0 for outer i boundary */
      Flux[j  ] = eps*(U[j] - 0.5*(1.0 + eps)*dUm);
    }
  }

  return;
}


#endif /* SHEARING_BOX */
#endif /* RADIATION_TRANSFER */
