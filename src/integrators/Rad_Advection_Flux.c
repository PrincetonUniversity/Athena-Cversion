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
 */
#ifdef SHEARING_BOX
#ifdef FARGO
#ifndef RADFARGO
#error : rad_fargo must be used when fargo is used as we update advection term explicitly!.
#endif /* fargo */
#endif /* rad_fargo */
#endif /* shearing box */



static void vanLeer_slope(const Real Er1, const Real Er2, const Real Er3, Real *slope);
#ifdef SHEARING_BOX
#ifdef RADFARGO
static Real ***FargoVars = NULL;
static Real ***FargoFlx = NULL;
static Real *U=NULL, *Flx=NULL;
static int nfghost;
static void RemapRadFlux(const Real *U,const Real eps,const int ji,const int jo, Real *F);
#ifdef MPI_PARALLEL
static double *send_buf = NULL, *recv_buf = NULL;
#endif

#endif
#endif

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

		*x1Flux = -AdvFlag * dtodx1 * (AdvFx[1] - AdvFx[0]);
	


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

		*x1Flux = -AdvFlag * dtodx1 * (AdvFx[1] - AdvFx[0]);

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

		*x2Flux = -AdvFlag * dtodx2 * (AdvFx[1] - AdvFx[0]);
#ifdef SHEARING_BOX
#ifdef RADFARGO
		if(Erflag){
			jj = j - pG->js + nfghost;
			*x2Flux -= (FargoFlx[k][i][jj+1]-FargoFlx[k][i][jj]);
		}
#endif
#endif


  return;	
	

}




void Rad_Advection_Flux3D(const DomainS *pD, const int i, const int j, const int k, const Real AdvFlag, Real *x1Flux, Real *x2Flux, Real *x3Flux)
{
/* The returned flux is Er^{n+1} = Er^{n} - x1Flux - x2Flux - x3Flux */
	
	
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

		*x1Flux = -AdvFlag * dtodx1 * (AdvFx[1] - AdvFx[0]);

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

		*x2Flux = -AdvFlag * dtodx2 * (AdvFx[1] - AdvFx[0]);

	/* Add flux due to background shearing */
#ifdef SHEARING_BOX
#ifdef RADFARGO
		if(Erflag){
			jj = j - pG->js + nfghost;
			*x2Flux -= (FargoFlx[k][i][jj+1]-FargoFlx[k][i][jj]);
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

		*x3Flux = -AdvFlag * dtodx3 * (AdvFx[1] - AdvFx[0]);

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
          FargoFlx[k][i][j] = Flx[j-joffset];

/* Sum the flux from integer offset: +/- for +/- joffset */
          for (jj=1; jj<=joffset; jj++) {
            FargoFlx[k][i][j] += FargoVars[k][i][j-jj];
          }
          for (jj=(joffset+1); jj<=0; jj++) {
            FargoFlx[k][i][j] -= FargoVars[k][i][j-jj];
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

				pG->U[k][j][i].Er -= (FargoFlx[k][i][jj+1]-FargoFlx[k][i][jj]);

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

  if((FargoFlx=(Real***)calloc_3d_array(max3,max1,max2,sizeof(Real)))==NULL)
    ath_error("[Rad_Fargo_init]: malloc returned a NULL pointer\n");


#ifdef MPI_PARALLEL
  size1 = nghost*pG->Nx[1]*(pG->Nx[2]+1);
  size2 = nfghost*(pG->Nx[0]+1)*(pG->Nx[2]+1);
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
        if (FargoFlx  != NULL) free_3d_array(FargoFlx);
	
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



#endif /* radMHD_INTEGRATOR */


