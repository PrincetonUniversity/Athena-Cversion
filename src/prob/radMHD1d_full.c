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
  RadGridS *pRG = pDomain->RadGrid;
  int nmu = pRG->nmu, nf=pRG->nf;
  int ifr;
  int iang;
  Real xtop, xbtm;  
  Real den = 1.0;
  Real wedd[2] = {0.5, 0.5};
  Real muedd[2] = {-1.0 / sqrt(3.0), 1.0 / sqrt(3.0)};
  Real *mul = NULL, *wl = NULL, *tau = NULL;
#endif

  int i, j, k, iu, il, ju, jl, ku, kl;
  Real d0, u0, T0, x1, x2, x3, temperature;

  
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
#ifdef RADIATION_HYDRO
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");	
#endif

#ifdef RADIATION_TRANSFER
  
  iang = par_geti_def("problem","iang", 1);
  eps0 = par_getd("problem","eps");
/* Setup density structure */ 
/* tau is used to initialize density grid */
  if ((tau = (Real *)calloc_1d_array(pGrid->Nx[0]+nghost+2,sizeof(Real))) == NULL) {
    ath_error("[problem]: Error allocating memory");
  }

  xtop = pDomain->RootMaxX[0];
  xbtm = pDomain->RootMinX[0];

  for(i=pGrid->is; i<=pGrid->ie+2; i++) {
    cc_pos(pGrid, i, j,k, &x1, &x2, &x3);
    tau[i] = const_opacity(pGrid, 0,i,pGrid->js,pGrid->ks) * (xtop - x1);
  }


#endif


/* Initialize the grid including the ghost cells.  */

	d0 = 1.0;
	u0 = -20.0;
	T0 = 8.5;

/* First initialize the gas part */
	for (k=kl; k<=ku; k++) {
      	for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

	/* xmax = 2.0; */
	  temperature = 1.0 + 7.5 * x1 / 2.0;
          pGrid->U[k][j][i].d  = d0;
          pGrid->U[k][j][i].M1 = d0 * u0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 0.5 * d0 * u0 * u0 + R_ideal * d0 * temperature /(Gamma - 1.0);

	  pGrid->U[k][j][i].Edd_11 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	  pGrid->U[k][j][i].Sigma_t = 10.0;
	  pGrid->U[k][j][i].Sigma_a = 10.0;

#ifdef MHD
          pGrid->B1i[k][j][i] = 0.0;
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;
          pGrid->U[k][j][i].B1c = 0.0;
          pGrid->U[k][j][i].B2c = 0.0;
          pGrid->U[k][j][i].B3c = 0.0;
#endif
	}
	}
}

#ifdef RADIATION_TRANSFER

	/* Second, initialize the specific intensity along the boundary */
	
/* Initialize arrays of angles and weights for angular quadrature */
  	pRG->ng=1;
  	pRG->r1ms=0; pRG->r1me=0; pRG->l1ms=0; pRG->l1me=0;
  	if (iang == 1) { 
/* angles and weight determined by gauss-legendre */
    	pRG->r1ls=nmu/2; pRG->r1le=nmu-1; pRG->l1ls=0; pRG->l1le=nmu/2-1;

    	if ((mul = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      		ath_error("[rad1d]: Error allocating memory");
    	}

    	if ((wl = (Real *)calloc_1d_array(nmu/2,sizeof(Real))) == NULL) {
      		free_1d_array(mul);
      		ath_error("[rad1d]: Error allocating memory");
    	}
    	gauleg(0.0, 1.0, mul, wl, nmu/2);
    	for(i=0; i<nmu/2; i++) {
      		pRG->mu[i] = -mul[nmu/2-i-1];   
      		pRG->w[i][0] = 0.5 * wl[nmu/2-i-1];  
    	}
    	for(i=nmu/2; i<nmu; i++) {
      		pRG->mu[i] = mul[i-nmu/2];   
      		pRG->w[i][0] = 0.5 * wl[i-nmu/2];  
    	} 
/* free up memory */
    	free_1d_array(mul);
    	free_1d_array(wl);
     
  	} else if (iang == 2) {
/* Angles reproduce Eddington approximation */
    	nmu = 2;
	pRG->nmu = 2;
    	pRG->r1ls=1; pRG->r1le=1; pRG->l1ls=0; pRG->l1le=0;
    	for(i=0; i<nmu; i++) {
      		pRG->mu[i] = muedd[i];   
      		pRG->w[i][0] = wedd[i];  
    	}
  	}

/* ------- Initialize boundary emission ---------------------------------- */

  /* incident radiation at lower (iz=0)  boundary */
  /* lower boundary is tau=0, no irradiation */
  	for(j=nmu/2; j<nmu; j++) 
    	for(ifr=0; ifr<nf; ifr++) 
    		pRG->l1imu[pRG->ks][pRG->js][ifr][j][0] = Thermal_B(pGrid, ifr, pGrid->is, pGrid->js, pGrid->ks); 

  /* incident radiation at upper (iz=nx) boundary */
  /* upper boundary is large tau, eps=1 */
/* Should use reflected boundary condition */
  	for(j=0; j<nmu/2; j++) 
    	for(ifr=0; ifr<nf; ifr++) 
      		pRG->r1imu[pRG->ks][pRG->js][ifr][j][0] = Thermal_B(pGrid, ifr, pGrid->ie, pGrid->js, pGrid->ks);

/* enroll radiation specification functions */
	get_thermal_source = Thermal_B;
	get_thermal_fraction = const_eps;
	get_total_opacity = const_opacity;

/* Free up memory */
  	free_1d_array(tau);


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

	 temperature = 1.0 + 7.5 * x1 / 2.0;
	
#ifdef RADIATION_HYDRO
	  pGrid->U[k][j][i].Er = temperature * temperature * temperature * temperature ;
	  pGrid->U[k][j][i].Fr1 = -30.0 * pGrid->U[k][j][i].Edd_11 * temperature * temperature * temperature /                                           pGrid->U[k][j][i].Sigma_t;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;	 		
#endif
        }
      }
    }

	bvals_mhd_fun(pDomain, right_x1, radMHD_inflow);
	bvals_rad_fun(pDomain, right_x1, radMHD_rad_inflow);

  return;
}

void radMHD_inflow(GridS *pGrid)
{
  	int i, ie;
	int ks, js;
	ie = pGrid->ie;
  	ks = pGrid->ks;
	js = pGrid->js;

	Real d0, u0, T0;
	d0 = 1.0;
	u0 = -20.0;
	T0 = 8.5;
 
    for (i=1;  i<=nghost;  i++) {
      pGrid->U[ks][js][ie+i].d  = d0;
      pGrid->U[ks][js][ie+i].M1 = d0 * u0;
      pGrid->U[ks][js][ie+i].E  = 0.5 * d0 * u0 * u0 + R_ideal * d0 * T0/(Gamma - 1.0);
    }
  
}


void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, ie;
	int ks, js;
	ie = pGrid->ie;
  	ks = pGrid->ks;
	js = pGrid->js;

	for (i=1;  i<=nghost;  i++) {
      	pGrid->U[ks][js][ie+i].Er  = 1.0;
	pGrid->U[ks][js][ie+i].Sigma_t = 10.0;
	pGrid->U[ks][js][ie+i].Sigma_a = 10.0;
	pGrid->U[ks][js][ie+i].Edd_11 = 1.0/3.0;
      	pGrid->U[ks][js][ie+i].Fr1 = -30.0 * pGrid->U[ks][js][ie+i].Edd_11 * 8.5 * 8.5 * 8.5/pGrid->U[ks][js][ie+i].Sigma_t;
    }
	

  
}

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
/* During the first step, to update the raidative flux 
 * after the Eddington tensor is updated 
 *
 */
/*
#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)
#ifdef RADIATION_TRANSFER

	int nl, nd;
	GridS *pGrid;

	int i, j, k, iu, il, ju, jl, ku, kl;
	double temperature, energy, momentum, density, pressure;
	

	for (nl=0; nl<(pM->NLevels); nl++){ 
      		for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			pGrid = pM->Domain[nl][nd].Grid;
	
	
  			iu = pGrid->ie;
  			il = pGrid->is;

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


			if(pM->nstep == 0){

			 for (k=kl; k<=ku; k++) {
      				for (j=jl; j<=ju; j++) {
        				for (i=il; i<=iu; i++) {
						energy = pGrid->U[k][j][i].E;
						momentum = pGrid->U[k][j][i].M1;
						density = pGrid->U[k][j][i].d;	
						pressure = (energy - 0.5 * momentum * momentum / density) * (Gamma - 1.0);
						temperature = pressure /(density * R_ideal);
						pGrid->U[k][j][i].Fr1 = -30.0 * pGrid->U[k][j][i].Edd_11 * temperature * temperature * temperature / pGrid->U[k][j][i].Sigma_t;
					}
				}
			}



			}
		}
	}
#endif
#endif

  return;
*/
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





