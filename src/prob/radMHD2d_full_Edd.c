#include "copyright.h"
/*==============================================================================
 * FILE: radMHD1d.c
 *
 * PURPOSE: Problem generator to test the radiation MHD code. 
 *  Development is underway. To be modified later 
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */
void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
void radMHD_inflow2(GridS *pGrid);
void radMHD_rad_inflow2(GridS *pGrid);


#ifdef RADIATION_TRANSFER

static Real eps0;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * gauleg()           - gauss-legendre quadrature from NR
 *============================================================================*/
void gauleg(Real x1, Real x2,  Real *x, Real *w, int n);
Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);
void Eddington_FUN(GridS *pG, RadGridS *pRG);
#endif


void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
#ifdef RADIATION_TRANSFER
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int nmu=pRG->nmu, ng=pRG->ng, nf=pRG->nf;
  int i, j, k, ifr, l;
  int iang;
  Real ytop, ybtm;  
  Real den = 1.0;
  Real wedd[2] = {0.5, 0.5};
  Real muedd[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  Real stheta;
  Real *mul = NULL, *wl = NULL, *tau = NULL, *gl = NULL, *wgl = NULL;

#endif

  int iu, il, ju, jl, ku, kl;
  Real d0, u0, T0, x1, x2, x3, temperature;
  int shifti, shiftj;

#ifndef RADIATION_TRANSFER
  
  int i, j, k;

#endif

  
/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }


/* Parse global variables of unit ratio */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");	
#endif

#ifdef RADIATION_TRANSFER
  
  iang = par_geti_def("problem","iang", 1);
  eps0 = par_getd("problem","eps");
/* Setup density structure */ 
/* tau is used to initialize density grid */
 /* if ((tau = (Real *)calloc_1d_array(pGrid->Nx[1]+nghost+2,sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory");
  }

  ytop = pDomain->RootMaxX[1];
  ybtm = pDomain->RootMinX[1];

  for(j=pGrid->js; j<=pGrid->je+1; j++) {
    x2 = pGrid->MinX[1] + (Real)(j-js)*pGrid->dx2;
    tau[j] = pow(10.0,-3.0 + 10.0 * ((x2-ybtm)/(ytop-ybtm)));
  }
*/
#endif


/* Initialize the grid including the ghost cells.  */
 /* Notice that ghost zones must be updated here. 
  * radiation transfer part need to use one zone 
  */

	d0 = 1.0;
	u0 = -20.0;
	T0 = 1.0;

	shiftj = (int)((ju-jl)*2/5);


/* First initialize the gas part */
	for (k=kl; k<=ku; k++) {
      	for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
	if((j>=jl+shiftj) && (j<=ju-shiftj)){

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

	/* xmax = 2.0; */
/*	  temperature = 1.0 + 7.5 * x1 / 2.0;
*/

#ifdef RADIATION_MHD
	  pGrid->B1i[k][j][i] = sqrt(fabs(u0));
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;

#endif


	  temperature = 1.0;	  
          pGrid->U[k][j][i].d  = d0;
	  if(fabs(x1)>0.0)
          pGrid->U[k][j][i].M1 = d0 * u0 * x1/fabs(x1);
	  else
	  pGrid->U[k][j][i].M1 = 0.0;
          
	  pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 0.5 * d0 * u0 * u0 + R_ideal * d0 * temperature /(Gamma - 1.0) + 0.5 * fabs(u0);

	  pGrid->U[k][j][i].Edd_11 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	  pGrid->U[k][j][i].Edd_22 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	  pGrid->U[k][j][i].Sigma_t = 10.0;
	  pGrid->U[k][j][i].Sigma_a = 10.0;



	}
	else{

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

	/* xmax = 2.0; */
/*	  temperature = 1.0 + 7.5 * x1 / 2.0;
*/
	   temperature = 1.0;		  
          pGrid->U[k][j][i].d  = d0;
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = R_ideal * d0 * temperature /(Gamma - 1.0) + 0.5* fabs(u0);

	  pGrid->U[k][j][i].Edd_11 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	  pGrid->U[k][j][i].Edd_22 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	  pGrid->U[k][j][i].Sigma_t = 10.0;
	  pGrid->U[k][j][i].Sigma_a = 10.0;
#ifdef RADIATION_MHD
	  pGrid->B1i[k][j][i] = sqrt(fabs(u0));
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;

#endif


	}


	}
	}
}

#ifdef RADIATION_MHD

 for (k=pGrid->ks; k<=pGrid->ke; k++) {
      for (j=pGrid->js; j<=pGrid->je; j++) {
        for (i=pGrid->is; i<=pGrid->ie; i++) {
	  pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
          pGrid->U[k][j][i].B3c = 0.0;


	}
	}
      }
#endif


#ifdef RADIATION_TRANSFER

	/* Second, initialize the specific intensity along the boundary */
	
/* Initialize arrays of angles and weights for angular quadrature */

  	
  if (iang == 1) { 
/* angles and weight determined by gauss-legendre */
     pRG->r1ls=0; pRG->r1le=nmu-1; pRG->l1ls=0; pRG->l1le=nmu-1;
    pRG->r1ms=ng/2; pRG->r1me=ng-1; pRG->l1ms=0; pRG->l1me=ng/2-1;
    pRG->r2ls=nmu/2; pRG->r2le=nmu-1; pRG->l2ls=0; pRG->l2le=nmu/2-1;
    pRG->r2ms=0; pRG->r2me=ng-1; pRG->l2ms=0; pRG->l2me=ng-1;
    if ((mul = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((wl = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      free_1d_array(mul);
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((gl = (Real *)calloc_1d_array(ng,sizeof(Real))) == NULL) {
      ath_error("[rad1d]: Error allocating memory");
    }
    if ((wgl = (Real *)calloc_1d_array(ng,sizeof(Real))) == NULL) {
      free_1d_array(mul);
      ath_error("[rad1d]: Error allocating memory");
    }
    gauleg(0.0, 1.0, mul, wl, nmu/2);
    gauleg(0.0, 2.0 * PI, gl, wgl, ng);

    for(i=0; i<nmu/2; i++) {
      pRG->mu[i] = -mul[nmu/2-i-1];
      stheta = sqrt(1.0 - pRG->mu[i] * pRG->mu[i]);
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = 0.25 * wl[nmu/2-i-1] * wgl[j] / PI; 
	pRG->gamma[i][j] = cos(gl[j]) * stheta;  
      } 
    }
    for(i=nmu/2; i<nmu; i++) {
      pRG->mu[i] = mul[i-nmu/2];   
      stheta = sqrt(1.0 - pRG->mu[i] * pRG->mu[i]);
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = 0.25 * wl[i-nmu/2] * wgl[j] / PI; 
	pRG->gamma[i][j] = cos(gl[j]) * stheta; 
      } 
    }
  

/* free up memory */
    free_1d_array(mul);
    free_1d_array(wl);
    free_1d_array(gl);
    free_1d_array(wgl);     
  } else if (iang == 2) {
/* Angles reproduce Eddington approximation */
     nmu = pRG->nmu = 2; ng = pRG->ng = 2;
    pRG->r1ls=0; pRG->r1le=1; pRG->l1ls=0; pRG->l1le=1;
    pRG->r1ms=1; pRG->r1me=1; pRG->l1ms=0; pRG->l1me=0;
    pRG->r2ls=1; pRG->r2le=1; pRG->l2ls=0; pRG->l2le=0;
    pRG->r2ms=0; pRG->r2me=1; pRG->l2ms=0; pRG->l2me=1;
    for(i=0; i<nmu; i++) {
      pRG->mu[i] = muedd[i];
      for(j=0; j<ng; j++) {
	pRG->w[i][j] = wedd[i] * wedd[j];
	pRG->gamma[i][j] = muedd[j];
      }
    }
  }

/* ------- Initialize boundary emission ---------------------------------- */

  /* incident radiation at lower (iz=0)  boundary */
  /* lower boundary is tau=0, no irradiation */
 /* set Ghost zones */
  for(j=0; j<nmu; j++) 
  	for(k=0; k<ng; k++) 
  		for(ifr=0; ifr<nf; ifr++) {
  		pRG->l1imu[pRG->ks][0][ifr][j][k] = 0.0; 
		pRG->r1imu[pRG->ks][0][ifr][j][k] = 0.0;
	}

 /* Left and right ghost zone at x direction */
 for(j=pRG->js; j<=pRG->je+1; j++) {
    for(i=pRG->r1ls; i<=pRG->r1le; i++) 
      for(k=pRG->r1ms; k<=pRG->r1me; k++)
	for(ifr=0; ifr<nf; ifr++)
	  pRG->r1imu[pRG->ks][j][ifr][i][k] = Thermal_B(pGrid, ifr, pGrid->ie+1, j+nghost-1, pGrid->ks);;

    for(i=pRG->l1ls; i<=pRG->l1le; i++) 
      for(k=pRG->l1ms; k<=pRG->l1me; k++) 
	for(ifr=0; ifr<nf; ifr++)
	   pRG->l1imu[pRG->ks][j][ifr][i][k] = Thermal_B(pGrid, ifr, pGrid->is-1, j+nghost-1, pGrid->ks);
  }


  for(i=pRG->is-1; i<=pRG->ie+1; i++) {
/* incident radiation at lower boundary */
    for(j=pRG->r2ls; j<=pRG->r2le; j++) 
      for(k=pRG->r2ms; k<=pRG->r2me; k++)
	for(ifr=0; ifr<nf; ifr++)
	  for(l=0; l<2; l++)
/* lower boundary is tau=0, no irradiation */
	    pRG->l2imu[pRG->ks][i][ifr][j][k][l] = 0.0;

/* incident radiation at upper boundary */
    for(j=pRG->l2ls; j<=pRG->l2le; j++) 
      for(k=pRG->l2ms; k<=pRG->l2me; k++) 
	for(ifr=0; ifr<nf; ifr++)
	  for(l=0; l<2; l++)
/* upper boundary is large tau, eps=1 */
	    pRG->r2imu[pRG->ks][i][ifr][j][k][l] = 0.0;
  }

/* enroll radiation specification functions */
	get_thermal_source = Thermal_B;
	get_thermal_fraction = const_eps;
	get_total_opacity = const_opacity;

/* Free up memory */
  /*	free_1d_array(tau);
*/

	/* Third, initialize Eddington tensor by calling Shane's function */

	hydro_to_rad(pDomain);  
/* solve radiative transfer */
	formal_solution(pDomain);
	/* Get the Eddington tensor */
	Eddington_FUN(pGrid, pRG);
	

#endif


     /* Now initialize the radiation energy density and flux */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

	  temperature = 1.0;
	
#if defined (RADIATION_HYDRO) || defined(RADIATION_MHD)
	  pGrid->U[k][j][i].Er = temperature * temperature * temperature * temperature ;
	  pGrid->U[k][j][i].Fr1 = 0.0;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;	 
	  
		
#endif
        }
      }
    }



	bvals_mhd_fun(pDomain, right_x1, radMHD_inflow);
	bvals_rad_fun(pDomain, right_x1, radMHD_rad_inflow);
	bvals_mhd_fun(pDomain, left_x1, radMHD_inflow2);
	bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow2);

  return;
}

void radMHD_inflow(GridS *pGrid)
{
  	int i, je,j, ju, jl, shiftj, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	Real d0, u0, T0;
	d0 = 1.0;
	u0 = -20.0;
	T0 = 1.0;

	
	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);
 
    for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
		if(j>=jl+shiftj && j<=ju-shiftj){
#ifdef RADIATION_MHD
		pGrid->B1i[ks][j][ie+i] = sqrt(fabs(u0));
		pGrid->U[ks][j][ie+i].B1c =  sqrt(fabs(u0));
#endif
      		pGrid->U[ks][j][ie+i].d  = d0;
      		pGrid->U[ks][j][ie+i].M1 = d0 * u0;
		pGrid->U[ks][j][ie+i].M2 = 0.0;
		pGrid->U[ks][j][ie+i].M3 = 0.0;
      		pGrid->U[ks][j][ie+i].E  = 0.5 * d0 * u0 * u0 + R_ideal * d0 * T0/(Gamma - 1.0) + 0.5 * fabs(u0);
      		}
	    else{
#ifdef RADIATION_MHD
		pGrid->B1i[ks][j][ie+i] = sqrt(fabs(u0));
		pGrid->U[ks][j][ie+i].B1c =  sqrt(fabs(u0));
#endif
		pGrid->U[ks][js][ie+i].d  = d0;
      		pGrid->U[ks][js][ie+i].M1 = 0.0;
		pGrid->U[ks][js][ie+i].M2 = 0.0;
		pGrid->U[ks][js][ie+i].M3 = 0.0;
      		pGrid->U[ks][js][ie+i].E  = R_ideal * d0 * T0/(Gamma - 1.0) + 0.5 * fabs(u0);
	  }
      }
    }
  
}


void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, je,j, ju, jl, shiftj, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);

	for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
/*      		if(j>=jl+shiftj && j<=ju-shiftj){
*/		pGrid->U[ks][j][ie+i].Er  = 1.0;
		pGrid->U[ks][j][ie+i].Sigma_t = 10.0;
		pGrid->U[ks][j][ie+i].Sigma_a = 10.0;
		pGrid->U[ks][j][ie+i].Edd_11 = 1.0/3.0;
		pGrid->U[ks][j][ie+i].Edd_22 = 1.0/3.0;
      		pGrid->U[ks][j][ie+i].Fr1 = 0.0;
		pGrid->U[ks][j][ie+i].Fr2 = 0.0;
		pGrid->U[ks][j][ie+i].Fr3 = 0.0;
/*		}
	   else{
		pGrid->U[ks][j][ie+i].Er  = 1.0;
		pGrid->U[ks][j][ie+i].Sigma_t = 10.0;
		pGrid->U[ks][j][ie+i].Sigma_a = 10.0;
		pGrid->U[ks][j][ie+i].Edd_11 = 1.0/3.0;
		pGrid->U[ks][j][ie+i].Edd_22 = 1.0/3.0;
      		pGrid->U[ks][j][ie+i].Fr1 = 0.0;
		pGrid->U[ks][j][ie+i].Fr2 = 0.0;
	
	 }
*/
    }
	}

  
}


void radMHD_inflow2(GridS *pGrid)
{
  	int i, je,j, ju, jl, shiftj, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	Real d0, u0, T0;
	d0 = 1.0;
	u0 = 20.0;
	T0 = 1.0;

	
	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);

 
    for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
#ifdef RADIATION_MHD
		pGrid->B1i[ks][j][is-i] = sqrt(fabs(u0));
		pGrid->U[ks][j][is-i].B1c =  sqrt(fabs(u0));
#endif
		if(j>=jl+shiftj && j<=ju-shiftj){
      		pGrid->U[ks][js][is-i].d  = d0;
      		pGrid->U[ks][js][is-i].M1 = d0 * u0;
      		pGrid->U[ks][js][is-i].E  = 0.5 * d0 * u0 * u0 + R_ideal * d0 * T0/(Gamma - 1.0) + 0.5 * fabs(u0);
      		}
	    else{
		
		pGrid->U[ks][js][is-i].d  = d0;
      		pGrid->U[ks][js][is-i].M1 = 0.0;
      		pGrid->U[ks][js][is-i].E  = R_ideal * d0 * T0/(Gamma - 1.0) +  0.5 * fabs(u0);
	  }
      }
    }
  
}


void radMHD_rad_inflow2(GridS *pGrid)
{
  	int i, je,j, ju, jl, shiftj, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);


	for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
      		if(j>jl+shiftj && j<ju-shiftj){
		pGrid->U[ks][j][is-i].Er  = 1.0;
		pGrid->U[ks][j][is-i].Sigma_t = 10.0;
		pGrid->U[ks][j][is-i].Sigma_a = 10.0;
		pGrid->U[ks][j][is-i].Edd_11 = 1.0/3.0;
		pGrid->U[ks][j][is-i].Edd_22 = 1.0/3.0;
      		pGrid->U[ks][j][is-i].Fr1 = 0.0;
		pGrid->U[ks][j][is-i].Fr2 = 0.0;
		}
	   else{
		pGrid->U[ks][j][is-i].Er  = 1.0;
		pGrid->U[ks][j][is-i].Sigma_t = 10.0;
		pGrid->U[ks][j][is-i].Sigma_a = 10.0;
		pGrid->U[ks][j][is-i].Edd_11 = 1.0/3.0;
		pGrid->U[ks][j][is-i].Edd_22 = 1.0/3.0;
      		pGrid->U[ks][j][is-i].Fr1 = 0.0;
		pGrid->U[ks][j][is-i].Fr2 = 0.0;
	
	 }

    }
	}

  
}
/*
void Scattering(const Real rho, const Real T, Real *Sigma_t, Real *Sigma_a){
	*Sigma_t = 10.0 * rho;
	*Sigma_a = *Sigma_t;

 return; 

}
*/
/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
	fprintf(fp,"%5.3e\n",Gamma);
	fprintf(fp,"%5.3e\n",Prat);
	fprintf(fp,"%5.3e\n",Crat);
	fprintf(fp,"%5.3e\n",R_ideal);
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{

	bvals_mhd_fun(&(pM->Domain[0][0]), right_x1, radMHD_inflow);
	bvals_rad_fun(&(pM->Domain[0][0]), right_x1, radMHD_rad_inflow);
	fscanf(fp,"%lf",&Gamma);
	fscanf(fp,"%lf\n",&Prat);
	fscanf(fp,"%lf\n",&Crat);
	fscanf(fp,"%lf\n",&R_ideal);


  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
 return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

#ifdef RADIATION_TRANSFER

Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
	Real density, momentum, energy, pressure, T, B;
	density = pG->U[k][j][i].d;
	momentum = pG->U[k][j][i].M1;
	energy = pG->U[k][j][i].E;
	pressure = (energy - 0.5 * momentum * momentum / density) * (Gamma - 1.0);
	T = pressure / (density * R_ideal);
	B = T * T * T * T / (4.0 * 3.141592653);

  return B;
}

Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k)
{

  return eps0;
  
}

Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k)
{

  return 10.0;
  
}

void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
{

  Real eps = 3.0e-14;
  Real xm, xl, z, z1;
  Real p1, p2, p3, pp;
  int m, i, j;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=1; i<=m; i++) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2.0 * (Real)j - 1.0) * z * p2 - ((Real)j - 1.0) * p3) / (Real)j;
      }
      pp = (Real)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}


#endif





