#include "copyright.h"
/*==============================================================================
 * FILE: set_bvals_shear.c
 *
 * PURPOSE: 3D shearing sheet boundary conditions at ix1 and ox1.  Called by
 *   set_bvals_mhd.  Decomposition of the Domain into MPI grids in X,Y and/or Z
 *   is allowed.  The RemapEy() function (which applies the shearing sheet
 *   boundary conditions to the y-component of the EMF to keep <Bz>=const) is
 *   called directly by the 3D integrator.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * ShearingSheet_ix1() - shearing sheet BCs on ix1, called by set_bval().
 * ShearingSheet_ox1() - shearing sheet BCs on ox1, called by set_bval().
 * RemapEy_ix1()      - sets Ey at ix1 in integrator to keep <Bz>=const. 
 * RemapEy_ox1()      - sets Ey at ox1 in integrator to keep <Bz>=const. 
 * set_bvals_shear_init() - allocates memory for arrays used here
 * set_bvals_shear_destruct() - frees memory for arrays used here
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
#ifdef SHEARING_BOX

/* Define number of variables to be remapped */
#ifdef BAROTROPIC /* BAROTROPIC EOS */
#ifdef HYDRO
 enum {NREMAP = 4};
#endif
#ifdef MHD
 enum {NREMAP = 8};
#endif
#else /* ADIABATIC or other EOS */
#ifdef HYDRO
 enum {NREMAP = 5};
#endif
#ifdef MHD
 enum {NREMAP = 9};
#endif
#endif /* EOS */

#ifdef MHD
#define NVAR_SHARE (NVAR + 3)
#else
#define NVAR_SHARE NVAR
#endif

/* Define structure which holds variables remapped by shearing sheet BCs */
typedef struct Remap_s{
  Real U[NREMAP];
#if (NSCALARS > 0)
  Real s[NSCALARS];
#endif
}Remap;

/* The memory for all the arrays below is allocated in set_bvals_shear_init */
/* Arrays of ghost zones containing remapped conserved quantities */
static Remap ***GhstZns=NULL, ***GhstZnsBuf=NULL;
/* 1D vectors for reconstruction in conservative remap step */
static Real *U=NULL, *Flx=NULL;
/* Arrays of Ey remapped at ix1/ox1 edges of Domain */
#ifdef MHD
static Real **tEyBuf=NULL;
#endif
/* temporary vector needed for 3rd order reconstruction in ghost zones */
#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
static Real *Uhalf=NULL;
#endif
/* MPI send and receive buffers */
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * RemapFlux() - 2nd or 3rd order reconstruction for remap in ghost zones
 *============================================================================*/

void RemapFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* ShearingSheet_ix1: 3D shearing-sheet BCs in x1.  It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in set_bvals.c
 *
 * This is a public function which is called by set_bvals_mhd() inside a
 * SHEARING_BOX macro.
 *----------------------------------------------------------------------------*/

void ShearingSheet_ix1(Grid *pG, Domain *pD)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,TH_omL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int err,sendto_id,getfrom_id;
  double *pd;
  Remap *pRemap;
  Gas *pGas;
  MPI_Request rq;
  MPI_Status stat;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  xmin = par_getd("grid","x2min");
  xmax = par_getd("grid","x2max");
  Ly = xmax - xmin;

  TH_omL = 1.5*Omega*Lx;
  yshear = TH_omL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched. */

  for(k=ks; k<=ke+1; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = is-nghost+i; n=3;
        GhstZns[k][i][j].U[0] = pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[1] = pG->U[k][j][ii].M1;
        GhstZns[k][i][j].U[2] = pG->U[k][j][ii].M2 + TH_omL*pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[3] = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        n++;
        GhstZns[k][i][j].U[n] = pG->U[k][j][ii].E + (0.5/GhstZns[k][i][j].U[0])*
          (SQR(GhstZns[k][i][j].U[2]) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        n++;
        GhstZns[k][i][j].U[n] = pG->U[k][j][ii].B1c;
        n++;
        GhstZns[k][i][j].U[n] = pG->B1i[k][j][ii];
        n++;
        GhstZns[k][i][j].U[n] = pG->B2i[k][j][ii];
        n++;
        GhstZns[k][i][j].U[n] = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) GhstZns[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    for(i=0; i<nghost; i++){

      for (n=0; n<(NREMAP); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].U[n];
        RemapFlux(U,epsi,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n] - (Flx[j+1]-Flx[j]);
        }
      }

#if (NSCALARS > 0)
      for (n=0; n<(NSCALARS); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].s[n];
        RemapFlux(U,epsi,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].s[n] = GhstZns[k][i][j].s[n] - (Flx[j+1]-Flx[j]);
        }
      }
#endif

    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx2;

        for(i=0; i<nghost; i++){
          for (n=0; n<(NREMAP); n++) {
            GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
          }
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) { 
            GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
          }
#endif
        }

      }
    }

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

  } else {
#ifdef MPI_PARALLEL
    get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx2);
    joverlap = joffset - Ngrids*pG->Nx2;

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2; 
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid_x2; 
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from GhstZnsBuf */

      cnt = nghost*joverlap*(pG->Nx3+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Irecv error at 1 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pd++) = pRemap->s[n];
#endif
          }
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Send error at 1 = %d\n",err);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in GhstZns */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Wait error at 1 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
  
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pd++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pd++);
#endif
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js+joverlap; j<=je; j++){
          jremap = j-joverlap;
          for(i=0; i<nghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) { 
              GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
            }
#endif
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      cnt = nghost*(pG->Nx2-joverlap)*(pG->Nx3+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Irecv error at 2 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pd++) = pRemap->s[n];
#endif
          }
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Send error at 2 = %d\n",err);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in GhstZns */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[ShearingSheet_ix1]: MPI_Wait error at 2 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pd++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pd++);
#endif
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
        n=3;
        pG->U[k][j][is-nghost+i].d  = GhstZns[k][i][j].U[0];
        pG->U[k][j][is-nghost+i].M1 = GhstZns[k][i][j].U[1];
        pG->U[k][j][is-nghost+i].M2 = GhstZns[k][i][j].U[2];
        pG->U[k][j][is-nghost+i].M3 = GhstZns[k][i][j].U[3];
#ifdef ADIABATIC
        n++;
        pG->U[k][j][is-nghost+i].E  = GhstZns[k][i][j].U[n];
#endif /* ADIABATIC */
#ifdef MHD
        n++;
        pG->U[k][j][is-nghost+i].B1c = GhstZns[k][i][j].U[n];
        n++;
        pG->B1i[k][j][is-nghost+i] = GhstZns[k][i][j].U[n];
        n++;
        pG->B2i[k][j][is-nghost+i] = GhstZns[k][i][j].U[n];
        n++;
        pG->B3i[k][j][is-nghost+i] = GhstZns[k][i][j].U[n];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][is-nghost+i].s[n] = GhstZns[k][i][j].s[n];
        }
#endif
      }
    }
  }

/* Copy the face-centered B3 component of the field at k=ke+1 */
#ifdef MHD
  for(j=js; j<=je; j++){
    for(i=0; i<nghost; i++){
      pG->B3i[ke+1][j][is-nghost+i] = GhstZns[ke+1][i][j].U[NREMAP-1];
    }
  }
#endif /* MHD */

/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B as average of remapped face centered B, except B1.
 * The value of B2c at j=je is incorrect since B2i[je+1] not yet set -- fix in
 * step 10 below */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=is-nghost; i<is; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in set_bavls_mhd.c) */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=is-nghost; i<is; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
          pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
          pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
          pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

          pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
          pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
          pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
        }
      }
    }
#ifdef MHD
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<is; i++) {
        pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
        pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
      }
    }
#endif /* MHD */

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in set_bvals_mhd.c */

  } else {
#ifdef MPI_PARALLEL

/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(pG->Nx3 + 1)*NVAR_SHARE;
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[ShearingSheet_ix1]: MPI_Irecv error at 3 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          *(pd++) = pGas->d;
          *(pd++) = pGas->M1;
          *(pd++) = pGas->M2;
          *(pd++) = pGas->M3;
#ifdef MHD
          *(pd++) = pGas->B1c;
          *(pd++) = pGas->B2c;
          *(pd++) = pGas->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pGas->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pGas->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[ShearingSheet_ix1]: MPI_Send error at 3 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[Shearing_Sheet_ix1]: MPI_Wait error at 3 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          pGas->d = *(pd++);
          pGas->M1 = *(pd++);
          pGas->M2 = *(pd++);
          pGas->M3 = *(pd++);
#ifdef MHD
          pGas->B1c = *(pd++);
          pGas->B2c = *(pd++);
          pGas->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pGas->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pGas->s[n] = *(pd++);
#endif
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[ShearingSheet_ix1]: MPI_Irecv error at 4 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          *(pd++) = pGas->d;
          *(pd++) = pGas->M1;
          *(pd++) = pGas->M2;
          *(pd++) = pGas->M3;
#ifdef MHD
          *(pd++) = pGas->B1c;
          *(pd++) = pGas->B2c;
          *(pd++) = pGas->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pGas->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pGas->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[ShearingSheet_ix1]: MPI_Send error at 4 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[ShearingSheet_ix1]: MPI_Wait error at 4 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=is-nghost; i<is; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          pGas->d = *(pd++);
          pGas->M1 = *(pd++);
          pGas->M2 = *(pd++);
          pGas->M3 = *(pd++);
#ifdef MHD
          pGas->B1c = *(pd++);
          pGas->B2c = *(pd++);
          pGas->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pGas->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pGas->s[n] = *(pd++);
#endif
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for(i=is-nghost; i<is; i++){
      pG->U[k][je  ][i].B2c = 0.5*(pG->B2i[k][je+1][i]+pG->B2i[k][je][i]);
      pG->U[k][js-1][i].B2c = 0.5*(pG->B2i[k][js-1][i]+pG->B2i[k][js][i]);
    }
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/* ShearingSheet_ox1: 3D shearing-sheet BCs in x1.  It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in set_bvals.c
 *
 * This is a public function which is called by set_bvals_mhd() inside a
 * SHEARING_BOX macro.
 *----------------------------------------------------------------------------*/

void ShearingSheet_ox1(Grid *pG, Domain *pD)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,n,joffset,jremap;
  Real xmin,xmax,Lx,Ly,TH_omL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int err,sendto_id,getfrom_id;
  double *pd;
  Remap *pRemap;
  Gas *pGas;
  MPI_Request rq;
  MPI_Status stat;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Compute the distance the computational domain has sheared in y */

  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  xmin = par_getd("grid","x2min");
  xmax = par_getd("grid","x2max");
  Ly = xmax - xmin;

  TH_omL = 1.5*Omega*Lx;
  yshear = TH_omL*pG->time;

/* Split this into integer and fractional peices of the Domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */

  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0.  If Domain has MPI decomposition in Y,
 * then it is possible that:  pD->Nx2 > joffset > pG->Nx2   */

  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 2. ------------------------------------------------------------------
 * Copy data into GhstZns array.  Note i and j indices are switched. */

  for(k=ks; k<=ke+1; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = ie+1+i; n=3;
        GhstZns[k][i][j].U[0] = pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[1] = pG->U[k][j][ii].M1;
        GhstZns[k][i][j].U[2] = pG->U[k][j][ii].M2 - TH_omL*pG->U[k][j][ii].d;
        GhstZns[k][i][j].U[3] = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        n++;
        GhstZns[k][i][j].U[n] = pG->U[k][j][ii].E + (0.5/GhstZns[k][i][j].U[0])*
          (SQR(GhstZns[k][i][j].U[2]) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        n++;
        GhstZns[k][i][j].U[n] = pG->U[k][j][ii].B1c;
        n++;
        GhstZns[k][i][j].U[n] = pG->B1i[k][j][ii];
        n++;
        GhstZns[k][i][j].U[n] = pG->B2i[k][j][ii];
        n++;
        GhstZns[k][i][j].U[n] = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) GhstZns[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy GhstZns into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    for(i=0; i<nghost; i++){

      for (n=0; n<(NREMAP); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].U[n];
        RemapFlux(U,epso,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].U[n] = GhstZns[k][i][j].U[n] - (Flx[j+1]-Flx[j]);
        }
      }

#if (NSCALARS > 0)
      for (n=0; n<(NSCALARS); n++) {
        for (j=js-nghost; j<=je+nghost; j++) U[j] = GhstZns[k][i][j].s[n];
        RemapFlux(U,epsi,js,je+1,Flx);
        for(j=js; j<=je; j++){
          GhstZnsBuf[k][i][j].s[n] = GhstZns[k][i][j].s[n] - (Flx[j+1]-Flx[j]);
        }
      }
#endif

    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into GhstZns.  */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx2;

        for(i=0; i<nghost; i++){
          for (n=0; n<(NREMAP); n++) {
            GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
          }
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) { 
            GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
          }
#endif
        }

      }
    }

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into GhstZns.  */

  } else {
#ifdef MPI_PARALLEL
    get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx2);
    joverlap = joffset - Ngrids*pG->Nx2;

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(overlap-1):je] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid_x2; 
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2; 
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from GhstZnsBuf */

      cnt = nghost*joverlap*(pG->Nx3+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Irecv error at 1 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);

            for (n=0; n<NREMAP; n++) *(pd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pd++) = pRemap->s[n];
#endif
          }
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Send error at 1 = %d\n",err);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in GhstZns
 */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Wait error at 1 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
  
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pd++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pd++);
#endif
          }
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from GhstZnsBuf into GhstZns.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-joverlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
          for(i=0; i<nghost; i++){
            for (n=0; n<(NREMAP); n++) {
              GhstZns[k][i][j].U[n]  = GhstZnsBuf[k][i][jremap].U[n];
            }
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) { 
              GhstZns[k][i][j].s[n] = GhstZnsBuf[k][i][jremap].s[n];
            }
#endif
          }
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+joverlap:je]
 * from GhstZnsBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      cnt = nghost*(pG->Nx2-joverlap)*(pG->Nx3+1)*(NREMAP+NSCALARS);
/* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Irecv error at 2 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZnsBuf[k][i][j]);
            for (n=0; n<NREMAP; n++) *(pd++) = pRemap->U[n];
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) *(pd++) = pRemap->s[n];
#endif
          }
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Send error at 2 = %d\n",err);

/* unpack data sent from [js+joverlap:je], and remap into cells in
 * [js:je-joverlap] in GhstZns */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[ShearingSheet_ox1]: MPI_Wait error at 2 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          for(i=0; i<nghost; i++){
            /* Get a pointer to the Remap structure */
            pRemap = &(GhstZns[k][i][j]);
            for (n=0; n<NREMAP; n++) pRemap->U[n] = *(pd++);
#if (NSCALARS > 0)
            for (n=0; n<NSCALARS; n++) pRemap->s[n] = *(pd++);
#endif
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
        n=3;
        pG->U[k][j][ie+1+i].d  = GhstZns[k][i][j].U[0];
        pG->U[k][j][ie+1+i].M1 = GhstZns[k][i][j].U[1];
        pG->U[k][j][ie+1+i].M2 = GhstZns[k][i][j].U[2];
        pG->U[k][j][ie+1+i].M3 = GhstZns[k][i][j].U[3];
#ifdef ADIABATIC
        n++;
        pG->U[k][j][ie+1+i].E  = GhstZns[k][i][j].U[n];
#endif /* ADIABATIC */
#ifdef MHD
        n++;
        pG->U[k][j][ie+1+i].B1c = GhstZns[k][i][j].U[n];
        n++;
        if(i>0) pG->B1i[k][j][ie+1+i] = GhstZns[k][i][j].U[n];
        n++;
        pG->B2i[k][j][ie+1+i] = GhstZns[k][i][j].U[n];
        n++;
        pG->B3i[k][j][ie+1+i] = GhstZns[k][i][j].U[n];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][ie+1+i].s[n] = GhstZns[k][i][j].s[n];
        }
#endif
      }
    }
  }

/* Copy the face-centered B3 component of the field at k=ke+1 */
#ifdef MHD
  for(j=js; j<=je; j++){
    for(i=0; i<nghost; i++){
      pG->B3i[ke+1][j][ie+1+i] = GhstZns[ke+1][i][j].U[NREMAP-1];
    }
  }
#endif /* MHD */

/*--- Step 7. ------------------------------------------------------------------
 * compute cell-centered B as average of remapped face centered B, except B1.
 * The value of B2c at j=je is incorrect since B2i[je+1] not yet set -- fix in
 * step 10 below */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=ie+1; i<=ie+nghost; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/*--- Step 8. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y (similar to
 * periodic_ix2() and periodic_ox2() in set_bavls_mhd.c) */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke; k++) {
      for(j=1; j<=nghost; j++){
        for(i=ie+1; i<=ie+nghost; i++){
          pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
          pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
          pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
          pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
          pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

          pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
          pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
          pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
        }
      }
    }
#ifdef MHD
    for (j=1; j<=nghost; j++) {
      for (i=ie+1; i<=ie+nghost; i++) {
        pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
        pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
      }
    }
#endif /* MHD */

/*--- Step 9. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in set_bvals_mhd.c */

  } else {
#ifdef MPI_PARALLEL

/* Post a non-blocking receive for the input data from the left grid */
    cnt = nghost*nghost*(pG->Nx3 + 1)*NVAR_SHARE;
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Irecv error at 3 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          *(pd++) = pGas->d;
          *(pd++) = pGas->M1;
          *(pd++) = pGas->M2;
          *(pd++) = pGas->M3;
#ifdef MHD
          *(pd++) = pGas->B1c;
          *(pd++) = pGas->B2c;
          *(pd++) = pGas->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pGas->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pGas->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Send error at 3 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Wait error at 3 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          pGas->d = *(pd++);
          pGas->M1 = *(pd++);
          pGas->M2 = *(pd++);
          pGas->M3 = *(pd++);
#ifdef MHD
          pGas->B1c = *(pd++);
          pGas->B2c = *(pd++);
          pGas->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pGas->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pGas->s[n] = *(pd++);
#endif
        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Irecv error at 4 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          *(pd++) = pGas->d;
          *(pd++) = pGas->M1;
          *(pd++) = pGas->M2;
          *(pd++) = pGas->M3;
#ifdef MHD
          *(pd++) = pGas->B1c;
          *(pd++) = pGas->B2c;
          *(pd++) = pGas->B3c;
          *(pd++) = pG->B1i[k][j][i];
          *(pd++) = pG->B2i[k][j][i];
          *(pd++) = pG->B3i[k][j][i];
#endif /* MHD */
#ifndef BAROTROPIC
          *(pd++) = pGas->E;
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) *(pd++) = pGas->s[n];
#endif
        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Send error at 4 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[ShearingSheet_ox1]: MPI_Wait error at 4 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        for (i=ie+1; i<=ie+nghost; i++){
          /* Get a pointer to the Gas cell */
          pGas = &(pG->U[k][j][i]);

          pGas->d = *(pd++);
          pGas->M1 = *(pd++);
          pGas->M2 = *(pd++);
          pGas->M3 = *(pd++);
#ifdef MHD
          pGas->B1c = *(pd++);
          pGas->B2c = *(pd++);
          pGas->B3c = *(pd++);
          pG->B1i[k][j][i] = *(pd++);
          pG->B2i[k][j][i] = *(pd++);
          pG->B3i[k][j][i] = *(pd++);
#endif /* MHD */
#ifndef BAROTROPIC
          pGas->E = *(pd++);
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++) pGas->s[n] = *(pd++);
#endif
        }
      }
    }
#endif /* MPI_PARALLEL */

  } /* end of step 9 - periodic BC in Y with MPI */

/*--- Step 10 ------------------------------------------------------------------
 * Fix B2c at j=je,js-1, now that B2i[je+1] has been set properly  */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (i=ie+1; i<=ie+nghost; i++){
      pG->U[k][je  ][i].B2c = 0.5*(pG->B2i[k][je+1][i]+pG->B2i[k][je][i]);
      pG->U[k][js-1][i].B2c = 0.5*(pG->B2i[k][js-1][i]+pG->B2i[k][js][i]);
    }
  }
#endif /* MHD */

  return;
}


/*------------------------------------------------------------------------------
 * RemapEy_ix1() - Remaps Ey at [is] due to background shear, and then
 * averages remapped and original field.  This guarantees the sums of Ey
 * along the x1 boundaries at [is] and [ie+1] are identical -- thus net Bz is
 * conserved
 *
 * This is a public function which is called by integrator (inside a
 * SHEARING_BOX macro).
 *----------------------------------------------------------------------------*/

#ifdef MHD
void RemapEy_ix1(Grid *pG, Domain *pD, Real ***emfy, Real **tEy)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,TH_omL,yshear,deltay,epsi;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int err,sendto_id,getfrom_id;
  double *pd;
  Real *pEy;
  MPI_Request rq;
  MPI_Status stat;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ix1()  */

  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  xmin = par_getd("grid","x2min");
  xmax = par_getd("grid","x2max");
  Ly = xmax - xmin;

  TH_omL = 1.5*Omega*Lx;
  yshear = TH_omL*pG->time;
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 1. ------------------------------------------------------------------
 * Copy Ey from [ie+1] into temporary array, using periodic BC in x1.
 * Requires MPI calls if NGrid_x1 > 1   */

  if (pD->NGrid_x1 == 1) {
    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        tEy[k][j] = emfy[k][j][ie+1];
      }
    }
  } else {

/* MPI calls to swap data */

#ifdef MPI_PARALLEL
    cnt = pG->Nx2*(pG->Nx3+1);
/* Post a non-blocking receive for the input data from remapEy_ox1 (listen L) */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx1_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ix1]: MPI_Irecv error at 1 = %d\n",err);

/* send Ey at [is] to ox1 (send L) -- this data is needed by remapEy_ox1 */

    pd = send_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
        pEy = &(emfy[k][j][is]);
        *(pd++) = *pEy;
      }
    }
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx1_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ix1]: MPI_Send error at 1 = %d\n",err);

/* Listen for data from ox1 (listen L), unpack and set temporary array */

    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ix1]: MPI_Wait error at 1 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pd++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Apply periodic BC in x2 to temporary array.  Requires MPI calls if 
 * NGrid_x2 > 1 */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=1; j<=nghost; j++){
        tEy[k][js-j] = tEy[k][je-(j-1)];
        tEy[k][je+j] = tEy[k][js+(j-1)];
      }
    }

  } else {

/* MPI calls to swap data */

#ifdef MPI_PARALLEL
    cnt = nghost*(pG->Nx3 + 1);
/* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ix1]: MPI_Irecv error at 2 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        pEy = &(tEy[k][j]);
        *(pd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ix1]: MPI_Send error at 2 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ix1]: MPI_Wait error at 2 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pd++);
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ix1]: MPI_Irecv error at 3 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        pEy = &(tEy[k][j]);
        *(pd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                     boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ix1]: MPI_Send error at 3 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ix1]: MPI_Wait error at 3 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pd++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy tEy into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    RemapFlux(tEy[k],epsi,js,je+1,Flx);
    for(j=js; j<=je; j++){
      tEyBuf[k][j] = tEy[k][j] - (Flx[j+1] - Flx[j]);
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into tEy.  */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j - joffset;
        if (jremap < (int)js) jremap += pG->Nx2;
        tEy[k][j]  = tEyBuf[k][jremap];
      }
    }

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into tEy.  */

  } else {
#ifdef MPI_PARALLEL
    get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx2);
    joverlap = joffset - Ngrids*pG->Nx2;

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [je-(joverlap-1):je] is sent to, and
 * data in [js:js+(joverlap-1)] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [je-(joverlap-1):je] from tEyBuf */

      cnt = joverlap*(pG->Nx3+1);
/* Post a non-blocking receive for the input data */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[RemapEy_ix1]: MPI_Irecv error at 4 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pd++) = *pEy;
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[RemapEy_ix1]: MPI_Send error at 4 = %d\n",err);


/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [je-(joverlap-1):je], and remap into cells in
 * [js:js+(joverlap-1)] in tEy
 */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[RemapEy_ix1]: MPI_Wait error at 4 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
            pEy = &(tEy[k][j]);
            *pEy = *(pd++);
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tEyBuf into tEy.  Cells in [js:je-joverlap] are shifted by
 * joverlap into [js+joverlap:je] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js+joverlap; j<=je; j++){
          jremap = j-joverlap;
          tEy[k][j]  = tEyBuf[k][jremap];
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js:je-joverlap]
 * from tEyBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      cnt = nghost*(pG->Nx2-joverlap)*(pG->Nx3+1);
/* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[RemapEy_ix1]: MPI_Irecv error at 5 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pd++) = *pEy;
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[RemapEy_ix1]: MPI_Send error at 5 = %d\n",err);

/* unpack data sent from [js:je-overlap], and remap into cells in
 * [js+joverlap:je] in tEy */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[RemapEy_ix1]: MPI_Wait error at 5 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pd++);
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now return remapped Ey */

  return;
}
#endif /* MHD */


/*------------------------------------------------------------------------------
 * RemapEy_ox1() - Remaps Ey at [ie+1] due to background shear, and then
 * averages remapped and original field.  This guarantees the sums of Ey
 * along the x1 boundaries at [is] and [ie+1] are identical -- thus net Bz is
 * conserved
 *
 * This is a public function which is called by integrator (inside a
 * SHEARING_BOX macro).
 *----------------------------------------------------------------------------*/

#ifdef MHD
void RemapEy_ox1(Grid *pG, Domain *pD, Real ***emfy, Real **tEy)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int j,k,joffset,jremap;
  Real xmin,xmax,Lx,Ly,TH_omL,yshear,deltay,epso;
#ifdef MPI_PARALLEL
  int my_iproc,my_jproc,my_kproc,cnt,jproc,joverlap,Ngrids;
  int err,sendto_id,getfrom_id;
  double *pd;
  Real *pEy;
  MPI_Request rq;
  MPI_Status stat;
#endif

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheet_ox1()  */

  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  xmin = par_getd("grid","x2min");
  xmax = par_getd("grid","x2max");
  Ly = xmax - xmin;

  TH_omL = 1.5*Omega*Lx;
  yshear = TH_omL*pG->time;
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epso = -(fmod(deltay,pG->dx2))/pG->dx2;

/*--- Step 1. ------------------------------------------------------------------
 * Copy Ey from [is] into temporary array, using periodic BC in x1.
 * Requires MPI calls if NGrid_x1 > 1   */

  if (pD->NGrid_x1 == 1) {
    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        tEy[k][j] = emfy[k][j][is];
      }
    }
  } else {

/* MPI calls to swap data */

#ifdef MPI_PARALLEL
    cnt = pG->Nx2*(pG->Nx3+1);
/* Post a non-blocking receive for the input data from remapEy_ix1 (listen R) */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx1_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ox1]: MPI_Irecv error at 1 = %d\n",err);

/* send Ey at [ie+1] to ix1 (send R) -- this data is needed by remapEy_ix1 */

    pd = send_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
        pEy = &(emfy[k][j][ie+1]);
        *(pd++) = *pEy;
      }
    }
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx1_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ox1]: MPI_Send error at 1 = %d\n",err);

/* Listen for data from ix1 (listen R), unpack and set temporary array */

    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ox1]: MPI_Wait error at 1 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pd++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Apply periodic BC in x2 to temporary array.  Requires MPI calls if 
 * NGrid_x2 > 1 */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=1; j<=nghost; j++){
        tEy[k][js-j] = tEy[k][je-(j-1)];
        tEy[k][je+j] = tEy[k][js+(j-1)];
      }
    }

  } else {

/* MPI calls to swap data */

#ifdef MPI_PARALLEL
    cnt = nghost*(pG->Nx3 + 1);
/* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ox1]: MPI_Irecv error at 2 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je-nghost+1; j<=je; j++){
        pEy = &(tEy[k][j]);
        *(pd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ox1]: MPI_Send error at 2 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ox1]: MPI_Wait error at 2 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js-nghost; j<=js-1; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pd++);
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[RemapEy_ox1]: MPI_Irecv error at 3 = %d\n",err);

    pd = send_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=js; j<=js+nghost-1; j++){
        pEy = &(tEy[k][j]);
        *(pd++) = *pEy;
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                     boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[RemapEy_ox1]: MPI_Send error at 3 = %d\n",err);

/* Wait to receive the input data from the left grid */
    err = MPI_Wait(&rq, &stat);
    if(err) ath_error("[RemapEy_ox1]: MPI_Wait error at 3 = %d\n",err);

    pd = recv_buf;
    for (k=ks; k<=ke+1; k++){
      for (j=je+1; j<=je+nghost; j++){
        pEy = &(tEy[k][j]);
        *pEy = *(pd++);
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Copy tEy into buffer, at the same time apply a conservative remap of
 * solution over the fractional part of grid cell */

  for(k=ks; k<=ke+1; k++) {
    RemapFlux(tEy[k],epso,js,je+1,Flx);
    for(j=js; j<=je; j++){
      tEyBuf[k][j] = tEy[k][j] - (Flx[j+1] - Flx[j]);
    }
  }

/*--- Step 4. ------------------------------------------------------------------
 * If no MPI decomposition in Y, apply shift over integer number of
 * grid cells during copy from buffer back into tEy.  */

  if (pD->NGrid_x2 == 1) {

    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je; j++){
        jremap = j + joffset;
        if (jremap > (int)je) jremap -= pG->Nx2;
        tEy[k][j]  = tEyBuf[k][jremap];
      }
    }

/*--- Step 5. ------------------------------------------------------------------
 * If Domain contains MPI decomposition in Y, then MPI calls are required for
 * the cyclic shift needed to apply shift over integer number of grid cells
 * during copy from buffer back into tEy.  */

  } else {
#ifdef MPI_PARALLEL
    get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

/* Find integer and fractional number of grids over which offset extends.
 * This assumes every grid has same number of cells in x2-direction! */
    Ngrids = (int)(joffset/pG->Nx2);
    joverlap = joffset - Ngrids*pG->Nx2;

/*--- Step 5a. -----------------------------------------------------------------
 * Find ids of processors that data in [js:js+(joverlap-1)] is sent to, and
 * data in [je-(joverlap-1):je] is received from.  Only execute if joverlap>0  */

    if (joverlap != 0) {

      jproc = my_jproc - (Ngrids + 1);
      if (jproc < 0) jproc += pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc + (Ngrids + 1);
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

/*--- Step 5b. -----------------------------------------------------------------
 * Pack send buffer and send data in [js:js+(joverlap-1)] from tEyBuf */

      cnt = joverlap*(pG->Nx3+1);
/* Post a non-blocking receive for the input data */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[RemapEy_ox1]: MPI_Irecv error at 4 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=js+(joverlap-1); j++) {
          pEy = &(tEyBuf[k][j]);
          *(pd++) = *pEy;
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                     remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[RemapEy_ox1]: MPI_Send error at 4 = %d\n",err);

/*--- Step 5c. -----------------------------------------------------------------
 * unpack data sent from [js:js+(joverlap-1)], and remap into cells in
 * [je-(joverlap-1):je] in tEy
 */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[RemapEy_ox1]: MPI_Wait error at 4 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=je-(joverlap-1); j<=je; j++) {
            pEy = &(tEy[k][j]);
            *pEy = *(pd++);
        }
      }

    }

/*--- Step 5d. -----------------------------------------------------------------
 * If shear is less one full Grid, remap cells which remain on same processor
 * from tEyBuf into tEy.  Cells in [js+joverlap:je] are shifted by
 * joverlap into [js:je-overlap] */

    if (Ngrids == 0) {

      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je-joverlap; j++){
          jremap = j+joverlap;
          tEy[k][j]  = tEyBuf[k][jremap];
        }
      }

/*--- Step 5e. -----------------------------------------------------------------
 * If shear is more than one Grid, pack and send data from [js+overlap:je]
 * from tEyBuf (this step replaces 5d) */

    } else {

/* index of sendto and getfrom processors in GridArray are -/+1 from Step 5a */

      jproc = my_jproc - Ngrids;
      if (jproc < 0) jproc += pD->NGrid_x2;
      sendto_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      jproc = my_jproc + Ngrids;
      if (jproc > (pD->NGrid_x2-1)) jproc -= pD->NGrid_x2;
      getfrom_id = pD->GridArray[my_kproc][jproc][my_iproc].id;

      cnt = (pG->Nx2-joverlap)*(pG->Nx3+1);
/* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, getfrom_id,
                      remap0_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[RemapEy_ox1]: MPI_Irecv error at 5 = %d\n",err);

      pd = send_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js+joverlap; j<=je; j++) {
          pEy = &(tEyBuf[k][j]);
          *(pd++) = *pEy;
        }
      }
      err = MPI_Send(send_buf, cnt, MPI_DOUBLE, sendto_id,
                   remap0_tag, MPI_COMM_WORLD);
      if(err) ath_error("[RemapEy_ox1]: MPI_Send error at 5 = %d\n",err);

/* unpack data sent from [js+overlap:je], and remap into cells in
 * [js:je-joverlap] in tEy */

      err = MPI_Wait(&rq, &stat);
      if(err) ath_error("[RemapEy_ox1]: MPI_Wait error at 5 = %d\n",err);

      pd = recv_buf;
      for (k=ks; k<=ke+1; k++) {
        for (j=js; j<=je-joverlap; j++) {
          pEy = &(tEy[k][j]);
          *pEy = *(pd++);
        }
      }
    } /* end of step 5e - shear is more than one Grid */

#endif /* MPI_PARALLEL */
  } /* end of step 5 - MPI decomposition in Y */

/*--- Step 6. ------------------------------------------------------------------
 * Now return remapped Ey */

  return;
}
#endif /* MHD */


/*----------------------------------------------------------------------------*/
/* set_bvals_shear_init: allocates memory for temporary arrays/buffers
 */

void set_bvals_shear_init(Grid *pG, Domain *pD)
{
  int nx2,nx3,size;
  nx2 = pG->Nx2 + 2*nghost;
  nx3 = pG->Nx3 + 2*nghost;

/* Allocate memory for temporary arrays and vectors */

  if((GhstZns=(Remap***)calloc_3d_array(nx3,nghost,nx2,sizeof(Remap)))==NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");

  if((GhstZnsBuf=(Remap***)calloc_3d_array(nx3,nghost,nx2,sizeof(Remap)))==NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");

  if((U = (Real*)malloc(nx2*sizeof(Real))) == NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");

  if((Flx = (Real*)malloc(nx2*sizeof(Real))) == NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");

#ifdef MHD
  if ((tEyBuf=(Real**)calloc_2d_array(nx3,nx2,sizeof(Real))) == NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");
#endif /* MHD */

#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
  if ((Uhalf = (Real*)malloc(nx2*sizeof(Real))) == NULL)
    ath_error("[set_bvals_shear_init]: malloc returned a NULL pointer\n");
#endif

/* allocate memory for send/receive buffers in MPI parallel calculations */

#ifdef MPI_PARALLEL
  size = nghost*pG->Nx2*(pG->Nx3+1)*(NREMAP+NSCALARS);

  if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[set_bvals_shear_init]: Failed to allocate send buffer\n");

  if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[set_bvals_shear_init]: Failed to allocate receive buffer\n");
#endif /* MPI_PARALLEL */

  return;
}

/*----------------------------------------------------------------------------*/
/* set_bvals_shear_destruct:  Free temporary arrays
 */

void set_bvals_shear_destruct(void)
{
  if (GhstZns    != NULL) free_3d_array(GhstZns);
  if (GhstZnsBuf != NULL) free_3d_array(GhstZnsBuf);
  if (U   != NULL) free(U);
  if (Flx != NULL) free(Flx);

#ifdef MHD
  if (tEyBuf != NULL) free_2d_array(tEyBuf);
#endif

#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

  return;
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * RemapFlux: computes "fluxes" of conserved variables for conservative remap
 * Input Arguments:
 *   U = 1D vector of conserved variable at cell centers along 1-D slice
 *   eps = fraction of a cell to be remapped
 * Output Arguments:
 *   Flux = fluxes of conserved variable at interfaces over [jinner:jouter]
 */

#ifdef SECOND_ORDER
/*------------------------------------------------------------------------------
 * RemapFlux(): second order reconstruction for conservative remap.
 * SECOND ORDER REMAP: piecewise linear reconstruction and min/mod limiters
 */

void RemapFlux(const Real *U, const Real eps,
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

#endif /* SECOND_ORDER */

#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
/*------------------------------------------------------------------------------
 * RemapFlux(): third order reconstruction for conservative remap. 
 * THIRD ORDER REMAP: Colella & Sekora extremum preserving algorithm (PPME)
 */

void RemapFlux(const Real *U, const Real eps,
               const int jinner, const int jouter, Real *Flux)
{
  int j,jl,ju;
  Real d2Uc,d2Ul,d2Ur,d2U,d2Ulim,lim_slope,Ulv,Urv,dU,U6,qa,qb,qc,qx;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju+1; j++) {
    Uhalf[j]=(7.0*(U[j-1]+U[j]) - (U[j-2]+U[j+1]))/12.0;
    d2Uc = 3.0*(U[j-1] - 2.0*Uhalf[j] + U[j]);
    d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
    d2Ur = (U[j-1] - 2.0*U[j  ] + U[j+1]);
    d2Ulim = 0.0;
    lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
    if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0) {
      d2Ulim = SIGN(d2Uc)*MIN(1.25*lim_slope,fabs(d2Uc));
    }
    Uhalf[j] = 0.5*((U[j-1]+U[j]) - d2Ulim/3.0);
  }

  for (j=jl; j<=ju; j++) {
    Ulv = Uhalf[j  ];
    Urv = Uhalf[j+1];

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = (U[j-1]-U[j])*(U[j]-U[j+1]);
    if (qa <= 0.0 && qb <= 0.0) {
      qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
      d2U  = -2.0*qc;
      d2Uc = (U[j-1] - 2.0*U[j  ] + U[j+1]);
      d2Ul = (U[j-2] - 2.0*U[j-1] + U[j  ]);
      d2Ur = (U[j  ] - 2.0*U[j+1] + U[j+2]);
      d2Ulim = 0.0;
      lim_slope = MIN(fabs(d2Ul),fabs(d2Ur));
      lim_slope = MIN(fabs(d2Uc),lim_slope);
      if (d2Uc > 0.0 && d2Ul > 0.0 && d2Ur > 0.0 && d2U > 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2Uc < 0.0 && d2Ul < 0.0 && d2Ur < 0.0 && d2U < 0.0) {
        d2Ulim = SIGN(d2U)*MIN(1.25*lim_slope,fabs(d2U));
      }
      if (d2U == 0.0) {
        Ulv = U[j];
        Urv = U[j];
      } else {
        Ulv = U[j] + (Ulv - U[j])*d2Ulim/d2U;
        Urv = U[j] + (Urv - U[j])*d2Ulim/d2U;
      }
    }

    qa = (Urv-U[j])*(U[j]-Ulv);
    qb = Urv-Ulv;
    qc = 6.0*(U[j] - 0.5*(Ulv+Urv));
    if (qa <= 0.0) {
      Ulv = U[j];
      Urv = U[j];
    } else if ((qb*qc) > (qb*qb)) {
      Ulv = 3.0*U[j] - 2.0*Urv;
    } else if ((qb*qc) < -(qb*qb)) {
      Urv = 3.0*U[j] - 2.0*Ulv;
    }

    dU = Urv - Ulv;
    U6 = 6.0*(U[j] - 0.5*(Ulv + Urv));

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      qx = TWO_3RDS*eps;
      Flux[j+1] = eps*(Urv - 0.75*qx*(dU - (1.0 - qx)*U6));

    } else {         /* eps always < 0 for outer i boundary */
      qx = -TWO_3RDS*eps;
      Flux[j  ] = eps*(Ulv + 0.75*qx*(dU + (1.0 - qx)*U6));
    }
  }

  return;
}

#endif /* THIRD_ORDER_EXTREMA_PRESERVING */

#endif /* SHEARING_BOX */
