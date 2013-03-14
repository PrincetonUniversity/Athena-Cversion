#include "../copyright.h"
/*==============================================================================
 * FILE: Rad_Advection_Flux.c
 *
 * PURPOSE: Calculate the advective raiation flux
 *
 *
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"



#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
/********Public function****************/
/*-------Rad_Advection_flux(): 
 * Calculate the advection flux, which is at the right hand side of the equation *
 * To reduce the numerical diffusion of advection flux, we use piecewise linear 
 * interpolation to reconstruct the E_r and velocity 
 * This function only calculate the flux for one cell 
 * Fargo scheme is used, background shearing is subtracted 
 * With STATIC MESH REFINEMENT, we need to restrict and correct flux at the boundary *
 */
#ifdef SHEARING_BOX
#ifdef FARGO
#ifndef RADFARGO
#error : rad_fargo must be used when fargo is used as we update advection term explicitly!.
#endif /* fargo */
#endif /* rad_fargo */
#endif /* shearing box */




#ifdef STATIC_MESH_REFINEMENT

static double **send_bufRC=NULL;
 
#ifdef MPI_PARALLEL
static double ***recv_bufRC=NULL;
static MPI_Request ***recv_rq=NULL;
static MPI_Request  **send_rq=NULL;
#endif



#endif /* Static mesh refinement */







static void vanLeer_slope(const Real Er1, const Real Er2, const Real Er3, Real *slope);
#ifdef SHEARING_BOX
#ifdef RADFARGO
static Real ***FargoVars = NULL;
static Real *U=NULL, *Flx=NULL;
static int nfghost;
static void RemapRadFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif

#endif /* Shearing box */
#endif /* RadFARGO */

/* Advection flux only include the vE term, vP_r term is not included here */
void Rad_Advection_Flux1D(const DomainS *pD, const int i, const int j, const int k, const Real AdvFlag, Real *x1Flux)
{
/* The returned flux is Er^{n+1} = Er^{n} - x1Flux - x2Flux*/
	
	GridS *pG = pD->Grid;
	Real dt = pG->dt;
	Real dtodx1 = dt/pG->dx1;

  
	int m; 

	/* temporary variables for the velocity terms */
	Real vx, vxi0;
	Real vFx, vFxi0, meanvFx, meanEr, slope;
	Real AdvFx[2]; 
	/* Actual advection flux used at the cell interface, some of them will be zero */

	/* First, calculate the flux along x direction */
	for(m=i; m<=i+1; m++){
		vxi0 = pG->U[k][j][m-1].M1 / pG->U[k][j][m-1].d;				

		vx = pG->U[k][j][m].M1 / pG->U[k][j][m].d;		

		vFx   = vx ;
		vFxi0 = vxi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[k][j][m-2].Er, pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, &slope);
			meanEr = pG->U[k][j][m-1].Er + (1.0 - meanvFx * dtodx1) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, pG->U[k][j][m+1].Er, &slope);
			meanEr = pG->U[k][j][m].Er - (1.0 + meanvFx * dtodx1) * 0.5 * slope;
		}
		
		AdvFx[m-i] = meanvFx * meanEr;
	}

	/*	*x1Flux = -AdvFlag * dtodx1 * (AdvFx[1] - AdvFx[0]);
	*/
		x1Flux[0] = -AdvFlag * dtodx1 * AdvFx[0];
		x1Flux[1] = -AdvFlag * dtodx1 * AdvFx[1];
	


  return;	
	

}





void Rad_Advection_Flux2D(const DomainS *pD, const int i, const int j, const int k, const Real AdvFlag, Real *x1Flux, Real *x2Flux)
{
/* The returned flux is Er^{n+1} = Er^{n} - x1Flux - x2Flux*/
	
	
	GridS *pG = pD->Grid;
	Real dt = pG->dt;
	Real dtodx1 = dt/pG->dx1;
	Real dtodx2 = dt/pG->dx2;

  
	int m; 
#ifdef SHEARING_BOX
	int jj;
#ifndef RADFARGO	
	Real vshear, x1, x2, x3;
#endif
#endif
	
	/* temporary variables for the velocity terms */
	Real vx, vy, vxi0, vyi0;
	Real vFx, vFxi0, meanvFx, meanEr, slope;
	Real AdvFx[2]; 
	/* Actual advection flux used at the cell interface, some of them will be zero */	


	/* First, calculate the flux along x direction */
	for(m=i; m<=i+1; m++){
		vxi0 = pG->U[k][j][m-1].M1 / pG->U[k][j][m-1].d;

		vx = pG->U[k][j][m].M1 / pG->U[k][j][m].d;

		vFx   = vx;
		vFxi0 = vxi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[k][j][m-2].Er, pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, &slope);
			meanEr = pG->U[k][j][m-1].Er + (1.0 - meanvFx * dtodx1) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, pG->U[k][j][m+1].Er, &slope);
			meanEr = pG->U[k][j][m].Er - (1.0 + meanvFx * dtodx1) * 0.5 * slope;
		}
		
		AdvFx[m-i] = meanvFx * meanEr;
	}

/*		*x1Flux = -AdvFlag * dtodx1 * (AdvFx[1] - AdvFx[0]);
*/
		x1Flux[0] = -AdvFlag * dtodx1 * AdvFx[0];
		x1Flux[1] = -AdvFlag * dtodx1 * AdvFx[1];


	/*========================================================*/
	/* Second, calculate the flux along y direction */
	for(m=j; m<=j+1; m++){
#ifdef SHEARING_BOX
#ifndef RADFARGO
		cc_pos(pG,i,m,k,&x1,&x2,&x3);
		vshear = qshear * Omega_0 * x1;
#endif
#endif	
		
		vyi0 = pG->U[k][m-1][i].M2 / pG->U[k][m-1][i].d;		
		vy   = pG->U[k][m][i].M2 / pG->U[k][m][i].d;

#ifdef SHEARING_BOX
#ifndef RADFARGO
		/* include background shearing */
		vyi0 -= vshear;
		vy   -= vshear; 
#endif
#endif		
		

		vFx   = vy;
		vFxi0 = vyi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[k][m-2][i].Er, pG->U[k][m-1][i].Er, pG->U[k][m][i].Er, &slope);
			meanEr = pG->U[k][m-1][i].Er + (1.0 - meanvFx * dtodx2) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[k][m-1][i].Er, pG->U[k][m][i].Er, pG->U[k][m+1][i].Er, &slope);
			meanEr = pG->U[k][m][i].Er   - (1.0 + meanvFx * dtodx2) * 0.5 * slope;
		}
		
		AdvFx[m-j] = meanvFx * meanEr;
	}

/*		*x2Flux = -AdvFlag * dtodx2 * (AdvFx[1] - AdvFx[0]);
*/
		x2Flux[0] = -AdvFlag * dtodx2 * AdvFx[0];
		x2Flux[1] = -AdvFlag * dtodx2 * AdvFx[1];


#ifdef SHEARING_BOX
#ifdef RADFARGO
		if(Erflag){
			jj = j - pG->js + nfghost;
		/*	*x2Flux -= (FargoFlx[k][i][jj+1]-FargoFlx[k][i][jj]);
		*/
			x2Flux[0] += (-pG->RadFargoFlx[k][i][jj]);
			x2Flux[1] += (-pG->RadFargoFlx[k][i][jj+1]);

		}
#endif
#endif


  return;	
	

}




void Rad_Advection_Flux3D(const DomainS *pD, const int i, const int j, const int k, const Real AdvFlag, Real *x1Flux, Real *x2Flux, Real *x3Flux)
{
/* The returned flux is Er^{n+1} = Er^{n} + (x1Flux[1] - x1Flux[0]) + (x2Flux[1] - x2Flux[0] + (x3Flux[1] - x3Flux[0]) */
	
	
	GridS *pG = pD->Grid;
	Real dt = pG->dt;
	Real dtodx1 = dt/pG->dx1;
	Real dtodx2 = dt/pG->dx2;
	Real dtodx3 = dt/pG->dx3;

  
	int m;
#ifdef SHEARING_BOX
	int jj;
#ifndef RADFARGO	
	Real vshear, x1, x2, x3;
#endif
#endif

	/* temporary variables for the velocity terms */
	Real vx, vy, vz, vxi0, vyi0, vzi0;
	Real vFx, vFxi0, meanvFx, meanEr, slope;
	Real AdvFx[2]; 
	/* Actual advection flux used at the cell interface, some of them will be zero */
	

	/* First, calculate the flux along x direction */
	for(m=i; m<=i+1; m++){

		vxi0 = pG->U[k][j][m-1].M1 / pG->U[k][j][m-1].d;				

		vx = pG->U[k][j][m].M1 / pG->U[k][j][m].d;
		
		

		vFx   = vx;
		vFxi0 = vxi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[k][j][m-2].Er, pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, &slope);
			meanEr = pG->U[k][j][m-1].Er + (1.0 - meanvFx * dtodx1) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[k][j][m-1].Er, pG->U[k][j][m].Er, pG->U[k][j][m+1].Er, &slope);
			meanEr = pG->U[k][j][m].Er - (1.0 + meanvFx * dtodx1) * 0.5 * slope;
		}
		
		AdvFx[m-i] = meanvFx * meanEr;
	}

		x1Flux[0] = -AdvFlag * dtodx1 * AdvFx[0];
		x1Flux[1] = -AdvFlag * dtodx1 * AdvFx[1];

	/*========================================================*/
	/* Second, calculate the flux along y direction */
	for(m=j; m<=j+1; m++){

#ifdef SHEARING_BOX
#ifndef RADFARGO
		cc_pos(pG,i,m,k,&x1,&x2,&x3);
		vshear = qshear * Omega_0 * x1;
#endif
#endif
			
		vyi0 = pG->U[k][m-1][i].M2 / pG->U[k][m-1][i].d;	
		
		vy = pG->U[k][m][i].M2 / pG->U[k][m][i].d;

#ifdef SHEARING_BOX
#ifndef RADFARGO
		/* include background shearing */
		vyi0 -= vshear;
		vy   -= vshear; 
#endif	
#endif

		vFx   = vy;
		vFxi0 = vyi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[k][m-2][i].Er, pG->U[k][m-1][i].Er, pG->U[k][m][i].Er, &slope);
			meanEr = pG->U[k][m-1][i].Er + (1.0 - meanvFx * dtodx2) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[k][m-1][i].Er, pG->U[k][m][i].Er, pG->U[k][m+1][i].Er, &slope);
			meanEr = pG->U[k][m][i].Er   - (1.0 + meanvFx * dtodx2) * 0.5 * slope;
		}
		
		AdvFx[m-j] = meanvFx * meanEr;
	}
		
		x2Flux[0] = -AdvFlag * dtodx2 * AdvFx[0];
		x2Flux[1] = -AdvFlag * dtodx2 * AdvFx[1];


	/* Add flux due to background shearing */
#ifdef SHEARING_BOX
#ifdef RADFARGO
		if(Erflag){
			jj = j - pG->js + nfghost;
		/*	*x2Flux -= (FargoFlx[k][i][jj+1]-FargoFlx[k][i][jj]);
		*/
			x2Flux[0] += (-pG->RadFargoFlx[k][i][jj]);
			x2Flux[1] += (-pG->RadFargoFlx[k][i][jj+1]);

		}
#endif
#endif

	/*========================================================*/
	/* Third, calculate the flux along z direction */
	for(m=k; m<=k+1; m++){
		
		vzi0 = pG->U[m-1][j][i].M3 / pG->U[m-1][j][i].d;
		
		vz = pG->U[m][j][i].M3 / pG->U[m][j][i].d;	
		

		vFx   = vz;
		vFxi0 = vzi0;
		
		meanvFx = 0.5 * (vFx + vFxi0);
		if(meanvFx > 0.0){
			vanLeer_slope(pG->U[m-2][j][i].Er, pG->U[m-1][j][i].Er, pG->U[m][j][i].Er, &slope);
			meanEr = pG->U[m-1][j][i].Er + (1.0 - meanvFx * dtodx3) * 0.5 * slope;
		}
		else{
			vanLeer_slope(pG->U[m-1][j][i].Er, pG->U[m][j][i].Er, pG->U[m+1][j][i].Er, &slope);
			meanEr = pG->U[m][j][i].Er   - (1.0 + meanvFx * dtodx3) * 0.5 * slope;
		}
		
		AdvFx[m-k] = meanvFx * meanEr;
	}

/*		*x3Flux = -AdvFlag * dtodx3 * (AdvFx[1] - AdvFx[0]);
*/
		x3Flux[0] = -AdvFlag * dtodx3 * AdvFx[0];
		x3Flux[1] = -AdvFlag * dtodx3 * AdvFx[1];


  return;	
	

}



#ifdef SHEARING_BOX
#ifdef RADFARGO
void Rad_Fargo_Pre(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int jfs = nfghost, jfe = pG->Nx[1] + nfghost - 1;
  int i,j,k,ku,jj,joffset;
  Real x1,x2,x3,yshear,eps;

#ifdef MPI_PARALLEL
  int ierr,cnt;
  double *pSnd,*pRcv;
  MPI_Request rq;
#endif

/*--- Step 1. ------------------------------------------------------------------
 * Copy data into FargoVars array.  Note i and j indices are switched.
 * Since the FargoVars array has extra ghost zones in y, it must be indexed
 * over [jfs:jfe] instead of [js:je]  */
  if (pG->Nx[2] > 1) ku=ke+1; else ku=ke;
  for(k=ks; k<=ku; k++) {
    for(j=jfs; j<=jfe+1; j++){
      for(i=is; i<=ie+1; i++){
        jj = j-(jfs-js);
        FargoVars[k][i][j] = pG->U[k][jj][i].Er; 
      }
    }
  }
/*--- Step 2. ------------------------------------------------------------------
 * With no MPI decomposition in Y, apply periodic BCs in Y to FargoVars array
 * (method is similar to periodic_ix2() and periodic_ox2() in bvals_mhd).
 * Note there are extra ghost zones in Y in the FargoVars array  */

  if (pD->NGrid[1] == 1) {

    for(k=ks; k<=ku; k++) {
      for(i=is; i<=ie+1; i++){
        for(j=1; j<=nfghost; j++){
          FargoVars[k][i][jfs-j] = FargoVars[k][i][jfe-(j-1)];
          FargoVars[k][i][jfe+j] = FargoVars[k][i][jfs+(j-1)];
        }
      }
    }

#ifdef MPI_PARALLEL
  } else {

/*--- Step 3. ------------------------------------------------------------------
 * With MPI decomposition in Y, use MPI calls to handle periodic BCs in Y (like
 * send_ox2/receive_ix1 and send_ix1/receive_ox2 pairs in bvals_mhd.c */

/* Post a non-blocking receive for the input data from the left grid */
/* There is only one variable */
    cnt = (pG->Nx[0]+1)*nfghost*(ku-ks+1);
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                    fargo_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfe-nfghost+1; j<=jfe; j++){        
          *(pSnd++) = FargoVars[k][i][j];
        }
      }
    }

/* send contents of buffer to the neighboring grid on R-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                   fargo_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfs-nfghost; j<=jfs-1; j++){
               FargoVars[k][i][j] = *(pRcv++);

        }
      }
    }

/* Post a non-blocking receive for the input data from the right grid */
    ierr = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pG->rx2_id,
                    fargo_tag, pD->Comm_Domain, &rq);

    pSnd = send_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfs; j<=jfs+nfghost-1; j++){
          *(pSnd++) = FargoVars[k][i][j];

        }
      }
    }

/* send contents of buffer to the neighboring grid on L-x2 */
    ierr = MPI_Send(send_buf, cnt, MPI_DOUBLE, pG->lx2_id,
                   fargo_tag, pD->Comm_Domain);

/* Wait to receive the input data from the left grid */
    ierr = MPI_Wait(&rq, MPI_STATUS_IGNORE);

    pRcv = recv_buf;
    for (k=ks; k<=ku; k++){
      for (i=is; i<=ie+1; i++){
        for (j=jfe+1; j<=jfe+nfghost; j++){
           FargoVars[k][i][j] = *(pRcv++);
        }
      }
    }
#endif /* MPI_PARALLEL */
  }

/*--- Step 4. ------------------------------------------------------------------
 * Compute fluxes, including both (1) the fractional part of grid cell, and
 * (2) the sum over integer number of cells  */
  for(k=ks; k<=ku; k++) {
    for(i=is; i<=ie+1; i++){

/* Compute integer and fractional peices of a cell covered by shear */
      cc_pos(pG, i, js, ks, &x1,&x2,&x3);

	/* rad fargo only solve the vEr part, vPr part is not included */
      yshear = -qshear*Omega_0*x1*pG->dt;


      joffset = (int)(yshear/pG->dx2);
      if (abs(joffset) > (jfs-js))
        ath_error("[bvals_shear]: FARGO offset exceeded # of gh zns\n");
      eps = (fmod(yshear,pG->dx2))/pG->dx2;



        for (j=jfs-nfghost; j<=jfe+nfghost; j++) U[j] = FargoVars[k][i][j];
        RemapRadFlux(U,eps,(jfs-joffset),(jfe+1-joffset),Flx);

        for(j=jfs; j<=jfe+1; j++){
          pG->RadFargoFlx[k][i][j] = Flx[j-joffset];

/* Sum the flux from integer offset: +/- for +/- joffset */
          for (jj=1; jj<=joffset; jj++) {
            pG->RadFargoFlx[k][i][j] += FargoVars[k][i][j-jj];
          }
          for (jj=(joffset+1); jj<=0; jj++) {
            pG->RadFargoFlx[k][i][j] -= FargoVars[k][i][j-jj];
          }
        }
      }

    }

	/* update Er with advection flux */

if(!Erflag){
	
   	for(k=ks; k<=ke; k++) {
    		for(j=js; j<=je; j++){
      			jj = j-js+nfghost;
      			for(i=is; i<=ie; i++){

				pG->U[k][j][i].Er -= (pG->RadFargoFlx[k][i][jj+1]-pG->RadFargoFlx[k][i][jj]);

      			}
    		}
  	}

	bvals_radMHD(pD);

}

  return;
}



/* Copy all the data along y direction into this array */
void Rad_Fargo_init(MeshS *pM)
{
	GridS *pG;
  int nl,nd,nx1,nx2,nx3,max1=0,max2=0,max3=0;
  Real xmin, xmax;
#ifdef MPI_PARALLEL
  int size1=0,size2=0,size;
#endif

  xmin = pM->RootMinX[0];
  xmax = pM->RootMaxX[0];

/* estimate extra ghost cells needed by FARGO, allocate arrays accordingly */
  nfghost = nghost + ((int)(1.5*CourNo*MAX(fabs(xmin),fabs(xmax))) + 2);
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

   max2 = max2 + 2*nfghost;

  if((U = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[Rad_Fargo_init]: malloc returned a NULL pointer\n");

  if((Flx = (Real*)malloc(max2*sizeof(Real))) == NULL)
    ath_error("[Rad_Fargo_init]: malloc returned a NULL pointer\n");

  if((FargoVars=(Real***)calloc_3d_array(max3,max1,max2,sizeof(Real)))
    ==NULL) ath_error("[Rad_Fargo_init]: malloc returned a NULL pointer\n");



#ifdef MPI_PARALLEL
  size1 = nghost*(max2-2*nghost-2*nfghost)*(max3-2*nghost+1);
  size2 = nfghost*(max1-2*nghost+1)*(max3-2*nghost+1);
  size = MAX(size1,size2);

  if((send_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate send buffer\n");

  if((recv_buf = (double*)malloc(size*sizeof(double))) == NULL)
    ath_error("[bvals_shear_init]: Failed to allocate receive buffer\n");
#endif /* MPI_PARALLEL */


}


void Rad_Fargo_destruct(void )
{
	if(U != NULL)
		free(U);
	if(Flx != NULL)
		free(Flx);

	if (FargoVars != NULL) free_3d_array(FargoVars);
	
#ifdef MPI_PARALLEL
  	if (send_buf != NULL) free(send_buf);
  	if (recv_buf != NULL) free(recv_buf);
#endif /* MPI_PARALLEL */

}



void RemapRadFlux(const Real *U, const Real eps,
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

#endif /* end rad fargo */
#endif


void vanLeer_slope(const Real Er1, const Real Er2, const Real Er3, Real *slope)
{

  Real lim_slope1,lim_slope2;
  Real dWc, dWl, dWr, dWg, dWm;	




/*--- Step 1. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables */
	  dWc = Er3 - Er1;
	  dWl = Er2 - Er1;
	  dWr = Er3 - Er2;
	  if(dWl * dWr > 0.0)
		  dWg = 2.0*dWl*dWr/(dWl+dWr);
	  else dWg = 0.0;
	  

/*--- Step 2. ------------------------------------------------------------------
 * Apply monotonicity constraints to differences in primitive vars. */
  
      dWm = 0.0;
      if (dWl*dWr > 0.0) {
        lim_slope1 = MIN(    fabs(dWl),fabs(dWr));
        lim_slope2 = MIN(0.5*fabs(dWc),fabs(dWg));
        dWm = SIGN(dWc)*MIN(2.0*lim_slope1,lim_slope2);
      }
    
	*slope = dWm;
 
  return;
}


#ifdef STATIC_MESH_REFINEMENT
/*==========================================================*/
/* Prepare the advective radiation energy flux, including ghost zones */
/* Calculate the advection flux with Er at the beginning of the step */
/*=============================================================*/



void GetAdvErFlx(MeshS *pM)
{

	int nl, nd, nDim;
	int i, j, k, il, iu, jl, ju, kl, ku;
	DomainS *pD;
	GridS *pG;
	Real Advxtemp[2], Advytemp[2], Advztemp[2];


	/* Determine the dimentionality of the problem */
	nDim = 1;
	for(i=1; i<3; i++)
		if(pM->Nx[i]>1) 
			nDim++;


	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				pD = &(pM->Domain[nl][nd]);
				pG = pD->Grid;

				if(Matghost > 3)
					ath_error("[AdvErFlx_pre]: Matghost: %d is larger than 3!\n",Matghost);

				il = pG->is - (Matghost);
				iu = pG->ie + (Matghost);
				if(nDim > 1){
					jl = pG->js - (Matghost);
					ju = pG->je + (Matghost);
				} 
				else{
					jl = pG->js;
					ju = pG->je;
				}

				if(nDim > 2){
					kl = pG->ks - (Matghost);
					ku = pG->ke + (Matghost);
				} 
				else{
					kl = pG->ks;
					ku = pG->ke;
				}
				
				/* Now calculate the advection flux for each grid for three directions */
				/* AdvErFlx i, j, k stores the flux from the left side of cell i, j, k */
				for(k=kl; k<=ku; k++){ 
				for(j=jl; j<=ju; j++){
				for(i=il; i<=iu; i++){
					
					if(nDim == 1){
						Rad_Advection_Flux1D(pD, i, j, k, 1.0, Advxtemp);
						/* i,j, k store the flux from the left boundary */
						/* right flux of i is the same as left flux of i+1 */
						/* In this way, right flux of the last cell is kept */
						pG->AdvErFlx[0][k][j][i]   = Advxtemp[0];
						pG->AdvErFlx[0][k][j][i+1] = Advxtemp[1];						

					}
					else if(nDim == 2){
						Rad_Advection_Flux2D(pD, i, j, k, 1.0, Advxtemp,Advytemp);
						pG->AdvErFlx[0][k][j][i]   = Advxtemp[0];	
						pG->AdvErFlx[0][k][j][i+1] = Advxtemp[1]; 
						pG->AdvErFlx[1][k][j][i]   = Advytemp[0];	
						pG->AdvErFlx[1][k][j+1][i] = Advytemp[1];						

					}
					else{
						Rad_Advection_Flux3D(pD, i, j, k, 1.0, Advxtemp,Advytemp,Advztemp);
						pG->AdvErFlx[0][k][j][i]   = Advxtemp[0];	
						pG->AdvErFlx[0][k][j][i+1] = Advxtemp[1];

						pG->AdvErFlx[1][k][j][i]   = Advytemp[0];
						pG->AdvErFlx[1][k][j+1][i] = Advytemp[1];

						pG->AdvErFlx[2][k][j][i]   = Advztemp[0];		
						pG->AdvErFlx[2][k+1][j][i] = Advztemp[1];		
					}

					
				}/* end i */
				}/* end j */
				}/* end k */

			}/* Finish if this core doesn't work in this grid */
		}/* Finish the domain at level nl */
	}/* Finish level nl */



}




/*==========================================================*/
/* Update Er with advection flux, including restriction, correction part */
/*=============================================================*/


void AdvErFlx_pre(MeshS *pM)
{

	int nl, nd, nDim, dim;
	int il, iu, jl, ju, kl, ku;
	int i, ics, ice, ips, ipe;
	int j, jcs, jce, jps, jpe;
	int k, kcs, kce, kps, kpe;
	DomainS *pD;
	GridS *pG;
	int nZeroRC;
	

	int ncg, npg, rbufN, start_addr, cnt, nFlx, count;
	double *pRcv, *pSnd;
	GridOvrlpS *pCO, *pPO;
	Real factor, flxdir;
	int shift; /* shift the cell position */
	/* flxdir to judge left or right hand side */

#ifdef MPI_PARALLEL
  	int ierr,mAddress,mIndex,mCount;
#endif

	/* Determine the dimentionality of the problem */
	nDim = 1;
	for(i=1; i<3; i++)
		if(pM->Nx[i]>1) 
			nDim++;
	
	/*------------------ step 1, prepare the advecion flux ------------------------------------*/
	/* First, calculate the flux for each grid at each domain and level */
	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				pD = &(pM->Domain[nl][nd]);
				pG = pD->Grid;

				if(Matghost > 3)
					ath_error("[AdvErFlx_pre]: Matghost: %d is larger than 3!\n",Matghost);

				il = pG->is - (Matghost);
				iu = pG->ie + (Matghost);
				if(nDim > 1){
					jl = pG->js - (Matghost);
					ju = pG->je + (Matghost);
				} 
				else{
					jl = pG->js;
					ju = pG->je;
				}

				if(nDim > 2){
					kl = pG->ks - (Matghost);
					ku = pG->ke + (Matghost);
				} 
				else{
					kl = pG->ks;
					ku = pG->ke;
				}
				
				

				/* Now update the solution, including one ghost zones */
				for(k=kl; k<=ku; k++){ 
				for(j=jl; j<=ju; j++){
				for(i=il; i<=iu; i++){
					if(nDim == 1){
						pG->U[k][j][i].Er += (pG->AdvErFlx[0][k][j][i+1] - pG->AdvErFlx[0][k][j][i]);
					}
					else if(nDim == 2){
						pG->U[k][j][i].Er += (pG->AdvErFlx[0][k][j][i+1] - pG->AdvErFlx[0][k][j][i]) 
									+ (pG->AdvErFlx[1][k][j+1][i] - pG->AdvErFlx[1][k][j][i]);
					}
					else{
						pG->U[k][j][i].Er += (pG->AdvErFlx[0][k][j][i+1] - pG->AdvErFlx[0][k][j][i]) 
									+ (pG->AdvErFlx[1][k][j+1][i] - pG->AdvErFlx[1][k][j][i])
									+ (pG->AdvErFlx[2][k+1][j][i] - pG->AdvErFlx[2][k][j][i]);
					}
				}
				}
				}


			}/* End if Grid is not NULL */
		}/* Finish looping all domains at this level */
	}/* Finish looping all levels */



	/*----------------------Restrict the flux---------------------------- */
	/* We also restrict Er in the active regions */
	/* Fine level solution will help reduce numerical diffusion due to advection */
	/* Restrict the flux at boundaries of refinement levels */

	/* This is adopted from RestrictCorrect function in smr , but we only need to do this for advection flux */
	for(nl=(pM->NLevels)-1; nl>=0; nl--){


#ifdef MPI_PARALLEL
/* Post non-blocking receives at level nl-1 for data from child Grids at this
 * level (nl).  This data is sent in Step 3 below, and will be read in Step 1
 * at the next iteration of the loop. */ 

  		if (nl>0) {
    			for (nd=0; nd<(pM->DomainsPerLevel[nl-1]); nd++){
      				if (pM->Domain[nl-1][nd].Grid != NULL) {
        				pG=pM->Domain[nl-1][nd].Grid;

/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0.
 * First index alternates between 0 and 1 for even/odd values of nl, since if
 * there are Grids on multiple levels there may be 2 receives posted at once */
        				mAddress = 0;
					nZeroRC = 0;
        				rbufN = ((nl-1) % 2);
        				for (ncg=(pG->NmyCGrid); ncg<(pG->NCGrid); ncg++){
						if(pG->CGrid[ncg].nWordsAdvEr == 0){
							nZeroRC += 1;
						}
						else{
          						mIndex = ncg - pG->NmyCGrid - nZeroRC;
          						ierr = MPI_Irecv(&(recv_bufRC[rbufN][nd][mAddress]),
            						pG->CGrid[ncg].nWordsAdvEr, MPI_DOUBLE, pG->CGrid[ncg].ID,
            						pG->CGrid[ncg].DomN, pM->Domain[nl-1][nd].Comm_Children,
            							&(recv_rq[nl-1][nd][mIndex]));
          						mAddress += pG->CGrid[ncg].nWordsAdvEr;
						}
        				}/* Finish looping all child grid */

      				}/* if Grid != NULL */
    			}/* Finish loop all domain at level nl */
  		}/* if nl > 0 */
#endif /* MPI_PARALLEL */


	/* Get child solution */
		for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){

  			if (pM->Domain[nl][nd].Grid != NULL) { /* there is a Grid on this processor */
    				pG=pM->Domain[nl][nd].Grid;
    				rbufN = (nl % 2);
				nZeroRC = 0;

				for(i=pG->NmyCGrid; i< pG->NCGrid; i++)
					if(pG->CGrid[i].nWordsAdvEr == 0)
						nZeroRC++;

    				for (ncg=0; ncg<(pG->NCGrid-nZeroRC); ncg++){

/*--- Step 1a. Get restricted solution and fluxes. ---------------------------*/

/* If child Grid is on this processor, set pointer to start of send buffer
 * loaded by this child Grid in Step 3 below during last iteration of loop. */

      					if (ncg < pG->NmyCGrid) {
        					pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);
        					pRcv = (double*)&(send_bufRC[pCO->DomN][0]);
      					} else {

#ifdef MPI_PARALLEL
/* Check non-blocking receives posted above for restricted solution from child
 * Grids, sent in Step 3 during last iteration of loop over nl.  Accept messages
 * in any order. */

        				mCount = pG->NCGrid - pG->NmyCGrid - nZeroRC;
        				ierr = MPI_Waitany(mCount,recv_rq[nl][nd],&mIndex,MPI_STATUS_IGNORE);
        				if(mIndex == MPI_UNDEFINED){
          					ath_error("[RestCorr]: Invalid request index nl=%i nd=%i\n",nl,nd);
        				}
      
/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0 */
        				
        				mIndex += pG->NmyCGrid;

					for(i=pG->NmyCGrid; i<=mIndex; i++)
						if(pG->CGrid[i].nWordsAdvEr == 0)	mIndex++;

					mAddress = 0;
        				for (i=pG->NmyCGrid; i<mIndex; i++) mAddress += pG->CGrid[i].nWordsAdvEr;
        					pCO=(GridOvrlpS*)&(pG->CGrid[mIndex]);
        					pRcv = (double*)&(recv_bufRC[rbufN][nd][mAddress]);
#else
/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        					ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",nl,nd);
#endif /* MPI_PARALLEL */
      					}/* if ncg not on the same CPU */


					if(pCO->nWordsAdvEr > 0){


					ics = pCO->ijks[0];
					ice = pCO->ijke[0];
					jcs = pCO->ijks[1];
					jce = pCO->ijke[1];
					kcs = pCO->ijks[2];
					kce = pCO->ijke[2];


					/* First, get the accurate solution from fine level */
					for(k=kcs; k<=kce; k++){
					for(j=jcs; j<=jce; j++){
					for(i=ics; i<=ice; i++){
						pG->U[k][j][i].Er = *(pRcv++);

					}/* End i */
					}/* End j */
					}/* End k */


					/* Flux at x1-face, for both the left and right hand side boundries */
					/* AdvErFlx[0] includes all flux along x direction */
					for(dim=0; dim<2; dim++){
						if(pCO->AdvEr[dim]){
							if(dim == 0){
								/* The left side boundary */
								i = ics;
								flxdir = 1.0;
								shift = -1;
							}
							if(dim == 1){
								i = ice + 1;
								flxdir = -1.0;	
								shift = 0;
							}
							/* [i][j][k] stores the flux at left boundary */
							for(k=kcs; k<=kce; k++){
							for(j=jcs; j<=jce; j++){
								/* As dt/dx is alread included in the flux, * 
								 * We need to account for the change of cell size */
								/* 0.5 is the ratio of cell size between parent and child grids */
								pG->U[k][j][i+shift].Er -= flxdir * (pG->AdvErFlx[0][k][j][i] - 0.5 * (*(pRcv++)));
	
							}/* Finish j */
							}/* Finish k */

						}
					}/* finish two interfaces at x boundary */

					/* Flux at x2 - face, for both the left and right hand size */
					/* AdvErFlx[1] includes all flux along y direction */
					for(dim=2; dim<4; dim++){
						if(pCO->AdvEr[dim]){
							if(dim == 2){
								/* The left side boundary */
								j = jcs;
								flxdir = 1.0;
								shift = -1;
							}
							if(dim == 3){
								j = jce + 1;
								flxdir = -1.0;
								shift  = 0;
							}
							/* [i][j][k] stores the flux at left boundary */
							for(k=kcs; k<=kce; k++){
							for(i=ics; i<=ice; i++){
								/* As dt/dx is alread included in the flux, * 
								 * We need to account for the change of cell size */
								/* 0.5 is the ratio of cell size between parent and child grids */
								pG->U[k][j+shift][i].Er -= flxdir * (pG->AdvErFlx[1][k][j][i] - 0.5 * (*(pRcv++)));	
							}/* Finish j */
							}/* Finish k */

						}
					}/* finish two interfaces at y boundary */

					/* Flux at x3 - face, for both the left and right hand size */
					/* AdvErFlx[2] includes all flux along z direction */
					for(dim=4; dim<6; dim++){
						if(pCO->AdvEr[dim]){
							if(dim == 4){
								/* The left side boundary */
								k = kcs;
								flxdir = 1.0;
								shift = -1;
							}
							if(dim == 5){
								k = kce + 1;
								flxdir = -1.0;
								shift = 0;
							}
							/* [i][j][k] stores the flux at left boundary */
							for(j=jcs; j<=jce; j++){
							for(i=ics; i<=ice; i++){
								/* As dt/dx is alread included in the flux, * 
								 * We need to account for the change of cell size */
								/* 0.5 is the ratio of cell size between parent and child grids */
								pG->U[k+shift][j][i].Er -= flxdir * (pG->AdvErFlx[2][k][j][i] - 0.5 * (*(pRcv++)));	
							}/* Finish j */
							}/* Finish k */

						}
					}/* finish two interfaces at z boundary */

					}/* if this child has restriction */
				}/* Finish all child grid */
			}/* If this grid is not null */
		}/* Finish all domain at this level */


		/*------Now restrict the flux --------------------------------*/
		for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				pG = pM->Domain[nl][nd].Grid;
				start_addr = 0;

				nZeroRC = 0;
				

				for(npg=0; npg<(pG->NPGrid); npg++){
					if(pG->PGrid[npg].nWordsAdvEr == 0){
						if(npg >= pG->NmyPGrid)
							nZeroRC += 1;
					}
					else{

					pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
					cnt = 0;
					
					ips = pPO->ijks[0];
					ipe = pPO->ijke[0];
      					jps = pPO->ijks[1];
      					jpe = pPO->ijke[1];
      					kps = pPO->ijks[2];
      					kpe = pPO->ijke[2];

					/* First the fine level solution */
					pSnd = (double*)&(send_bufRC[nd][start_addr]);
					for(k=kps; k<=kpe; k+=2){
					for(j=jps; j<=jpe; j+=2){
					for(i=ips; i<=ipe; i+=2){
						*(pSnd++) = pG->U[k][j][i].Er + pG->U[k][j][i+1].Er;
					}/* end i*/
					}/* end j */
					}/* end k */
					
					factor = 0.5;
					count = (ipe-ips+1)/2;	

					/* 2D and 3D problem */
					if(nDim > 1){
						pSnd = (double*)&(send_bufRC[nd][start_addr]);
						for(k=kps; k<=kpe; k+=2){
						for(j=jps; j<=jpe; j+=2){
						for(i=ips; i<=ipe; i+=2){
							*(pSnd++) += (pG->U[k][j+1][i].Er + pG->U[k][j+1][i+1].Er);
						}/* end i*/
						}/* end j */
						}/* end k */

						factor = 0.25;
						count = (jpe-jps+1)*(ipe-ips+1)/4;
					}/* End 2 or 3D problem */
	
					/* 3D problem */
					if(nDim > 2){
						pSnd = (double*)&(send_bufRC[nd][start_addr]);
						for(k=kps; k<=kpe; k+=2){
						for(j=jps; j<=jpe; j+=2){
						for(i=ips; i<=ipe; i+=2){
							*(pSnd++) += ((pG->U[k+1][j][i].Er + pG->U[k+1][j][i+1].Er)
									+ (pG->U[k+1][j+1][i].Er + pG->U[k+1][j+1][i+1].Er));
						}/* end i*/
						}/* end j */
						}/* end k */

						factor = 0.125;
						count = (jpe-jps+1)*(ipe-ips+1)*(kpe-kps+1)/8;
					}
					pSnd = (double*)&(send_bufRC[nd][start_addr]);

					for(i=0; i<count; i++) (*(pSnd++)) *= factor;
					cnt = count;

					/*-----------Restrict fluxes at x1 faces --------------*/
				

					for(dim=0; dim<2; dim++){
						if(pPO->AdvEr[dim]){	
							/* restrcit the flux */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							if(dim == 0)	i = ips;
							if(dim == 1)	i = ipe+1;

							if(nDim == 1){
								j = jps;
								k = kps;								
								*(pSnd++) = pG->AdvErFlx[0][k][j][i];
							}else{/* For 2D or 3D problem, average x2 direction for x1 flux */
								for(k=kps; k<=kpe; k+=2){
								for(j=jps; j<=jpe; j+=2){
									*(pSnd++) = pG->AdvErFlx[0][k][j][i] + pG->AdvErFlx[0][k][j+1][i];
								}/* End j */
								}/* End k */


								factor = 0.5;
								nFlx = (jpe-jps+1)/2;

								/* For 3D case, we need to rewrite the buffer */
								if(nDim == 3){	
									/* reset the pointer */
									pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
									for(k=kps; k<=kpe; k+=2){
									for(j=jps; j<=jpe; j+=2){
										*(pSnd++) += (pG->AdvErFlx[0][k+1][j][i] + pG->AdvErFlx[0][k+1][j+1][i]);
									}
									}/* end k */

									factor = 0.25;
									nFlx = ((kpe-kps+1)*(jpe-jps+1)/4);

								}/* End 3D case */
								
								/* for 2D and 3D case, multiple the normalization */
								pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
								for(i=0; i<nFlx; i++)	*(pSnd++) *= factor;	

							}/* End 2D or 3D */
							cnt += nFlx;
						}/* End if we need restrict this flux */
					}/* ENd for dim 0 and 1 */

					/*----------------------------------------------*/
					/*------Average x2 Faces --------------------------*/
					/*-----------Restrict fluxes at x2 faces --------------*/
					for(dim=2; dim<4; dim++){
						if(pPO->AdvEr[dim]){	
							/* restrcit the flux */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							if(dim == 2)	j = jps;
							if(dim == 3)	j = jpe+1;

							
							for(k=kps; k<=kpe; k+=2){
							for(i=ips; i<=ipe; i+=2){
								*(pSnd++) = pG->AdvErFlx[1][k][j][i] + pG->AdvErFlx[1][k][j][i+1];
							}/* End j */
							}/* End k */


							factor = 0.5;
							nFlx = (ipe-ips+1)/2;

							/* For 3D case, we need to rewrite the buffer */
							if(nDim == 3){	
								/* reset the pointer */
								pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
								for(k=kps; k<=kpe; k+=2){
								for(i=ips; i<=ipe; i+=2){
									*(pSnd++) += (pG->AdvErFlx[1][k+1][j][i] + pG->AdvErFlx[1][k+1][j][i+1]);
								}
								}/* end k */

								factor = 0.25;
								nFlx = ((kpe-kps+1)*(ipe-ips+1)/4);

							}/* End 3D case */
								
							/* for 2D and 3D case, multiple the normalization */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							for(i=0; i<nFlx; i++)	*(pSnd++) *= factor;	

							
							cnt += nFlx;
						}/* End if we need restrict this flux */
					}/* ENd for dim 2 and 3 */

					/*------------------------------------------------*/
					/*------Average x3 Faces --------------------------*/

					for(dim=4; dim<6; dim++){
						if(pPO->AdvEr[dim]){	
							/* restrcit the flux */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							if(dim == 4)	k = kps;
							if(dim == 5)	k = kpe+1;
							
							for(j=jps; j<=jpe; j+=2){
							for(i=ips; i<=ipe; i+=2){
								*(pSnd++) = pG->AdvErFlx[2][k][j][i] + pG->AdvErFlx[2][k][j][i+1];
							}/* End j */
							}/* End k */


							/* For 3D case, we need to rewrite the buffer */

							/* reset the pointer */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							for(j=jps; j<=jpe; j+=2){
							for(i=ips; i<=ipe; i+=2){
								*(pSnd++) += (pG->AdvErFlx[2][k][j+1][i] + pG->AdvErFlx[2][k][j+1][i+1]);
							}
							}/* end k */

							factor = 0.25;
							nFlx = ((jpe-jps+1)*(ipe-ips+1)/4);
							
							/* for 2D and 3D case, multiple the normalization */
							pSnd = (double*)&(send_bufRC[nd][start_addr+cnt]);
							for(i=0; i<nFlx; i++)	*(pSnd++) *= factor;	

							
							cnt += nFlx;
						}/* End if we need restrict this flux */
					}/* ENd for dim 2 and 3 */
					
					/*-------------------------------------------------*/
					/* now send the data */
#ifdef MPI_PARALLEL

				/* non-blocking send with MPI, using Domain number as tag.  */

      					if (npg >= pG->NmyPGrid){
        					mIndex = npg - pG->NmyPGrid - nZeroRC;
        					ierr = MPI_Isend(&(send_bufRC[nd][start_addr]), pG->PGrid[npg].nWordsAdvEr,
          					MPI_DOUBLE, pG->PGrid[npg].ID, nd, pM->Domain[nl][nd].Comm_Parent,
          						&(send_rq[nd][mIndex]));
      					}
#endif /* MPI_PARALLEL */

      						start_addr += pG->PGrid[npg].nWordsAdvEr;


					/*--------------------------------------------------*/
					}/* End if this parent grid has data to send */
				}/* Loop over all parent grid  */
			}/* If the CPU works in this grid */
		}/* Loop all domains at this level */


#ifdef MPI_PARALLEL
/*--- Step 4. Check non-blocking sends completed. ----------------------------*/
/* For MPI jobs, wait for all non-blocking sends in Step 3e to complete.  This
 * is more efficient if there are multiple messages per Grid. */

  		for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    			if (pM->Domain[nl][nd].Grid != NULL) {
      				pG=pM->Domain[nl][nd].Grid;

				nZeroRC = 0;
				/* Do not double count when NmyPGrid is Rad_nWordsRC =0 */
				for(i=pG->NmyPGrid; i<pG->NPGrid; i++)
					if(pG->PGrid[i].nWordsAdvEr == 0)	nZeroRC++;

      				if (pG->NPGrid > pG->NmyPGrid) {
        				mCount = pG->NPGrid - pG->NmyPGrid - nZeroRC;
        				ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      				}
    			}
  		}
#endif /* MPI_PARALLEL */


	}/* Finish looping all levels */
}


void AdvErFlx_init(MeshS *pM)
{

	int nl, nd, sendRC, recvRC, npg, ncg, maxND;
	int max_sendRC = 1, max_recvRC=1;
	int maxCG=1;

	GridS *pG;

	maxND = 1;

	for (nl=0; nl<(pM->NLevels); nl++) maxND=MAX(maxND,pM->DomainsPerLevel[nl]);
  	
	
	/* Loop over all grids to find maximum number of words for communication */
	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			sendRC = 0;
			recvRC = 0;
			if((pM->Domain[nl][nd].Grid) != NULL){
				pG = pM->Domain[nl][nd].Grid;
				for(npg=0; npg<pG->NPGrid; npg++)	sendRC += pG->PGrid[npg].nWordsAdvEr;
				for(ncg=0; ncg<pG->NCGrid; ncg++)	recvRC += pG->CGrid[ncg].nWordsAdvEr;
				
				max_sendRC = MAX(max_sendRC, sendRC);
				max_recvRC = MAX(max_recvRC, recvRC);

				
				maxCG = MAX(maxCG, pG->NCGrid);

			}/* End Grid is NULL */
		}/* ENd nd */
	}/* End nl */



	/*===== Allocate the memory ==========*/
		if((send_bufRC =
    			(double**)calloc_2d_array(maxND,max_sendRC,sizeof(double))) == NULL)
    				ath_error("[AdvEr_init]:Failed to allocate send_bufRC\n");

#ifdef MPI_PARALLEL
  		if((recv_bufRC =
    			(double***)calloc_3d_array(2,maxND,max_recvRC,sizeof(double))) == NULL)
    				ath_error("[AdvEr_init]: Failed to allocate recv_bufRC\n");
  		if((recv_rq = (MPI_Request***)
    			calloc_3d_array(pM->NLevels,maxND,maxCG,sizeof(MPI_Request))) == NULL)
    			ath_error("[AdvEr_init]: Failed to allocate recv MPI_Request array\n");
  		if((send_rq = (MPI_Request**)
    			calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    			ath_error("[AdvEr_init]: Failed to allocate send MPI_Request array\n");
#endif /* MPI_PARALLEL */




}


void AdvErFlx_destruct()
{


	
	if(send_bufRC != NULL){
		free_2d_array(send_bufRC);
		send_bufRC = NULL;
	}

#ifdef MPI_PARALLEL	
	if(recv_bufRC != NULL){
		free_3d_array(recv_bufRC);
		recv_bufRC = NULL;
	}	

	if(recv_rq != NULL){
		free_3d_array(recv_rq);
		recv_rq = NULL;
	}

	if(send_rq != NULL){
		free_2d_array(send_rq);
		send_rq = NULL;
	}

#endif /* MPI_PARALLEL */

}


#endif /* End STATIC MESH REFINEMENT */







#endif /* radMHD_INTEGRATOR */


