#include "copyright.h"
/*==============================================================================
 * FILE: radMHD2d.c
 *
 * PURPOSE: Problem generator to test the radiation MHD code. 
 *  Development is underway. To be modified later 
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */
/*
void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
*/

static double kappa=1.0;
void constopa(const Real rho, const Real T, Real *Sigma_t, Real *Sigma_a, Real dSigma[4]);
static double ran2(long int *idum);

static Real vorticity_p(const GridS *pG, const int i, const int j, const int k);
static Real vorticity_n(const GridS *pG, const int i, const int j, const int k);


static Real baronn(const GridS *pG, const int i, const int j, const int k);
static Real baronp(const GridS *pG, const int i, const int j, const int k);

static Real baron(const GridS *pG, const int i, const int j, const int k);

static Real dTdotvp(const GridS *pG, const int i, const int j, const int k);
static Real dTdotvn(const GridS *pG, const int i, const int j, const int k);


static Real Ntemperature(const GridS *pG, const int i, const int j, const int k);

static Real Ptemperature(const GridS *pG, const int i, const int j, const int k);

static Real baronnEr(const GridS *pG, const int i, const int j, const int k);
static Real baronpEr(const GridS *pG, const int i, const int j, const int k);


static Real vortSFr_n(const GridS *pG, const int i, const int j, const int k);
static Real vortSFr_p(const GridS *pG, const int i, const int j, const int k);

static Real vorticityFr_n(const GridS *pG, const int i, const int j, const int k);
static Real vorticityFr_p(const GridS *pG, const int i, const int j, const int k);


static Real dErTv_p(const GridS *pG, const int i, const int j, const int k);
static Real dErTv_n(const GridS *pG, const int i, const int j, const int k);

static Real Radwork_p(const GridS *pG, const int i, const int j, const int k);
static Real Radwork_n(const GridS *pG, const int i, const int j, const int k);



static Real vortS_n(const GridS *pG, const int i, const int j, const int k);

static Real vortS_p(const GridS *pG, const int i, const int j, const int k);

static Real divV_p(const GridS *pG, const int i, const int j, const int k);
static Real divV_n(const GridS *pG, const int i, const int j, const int k);

static Real dTdotv(const GridS *pG, const int i, const int j, const int k);

/* For transfer module */
#ifdef RADIATION_TRANSFER

static Real eps0;

static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);




#endif


void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  int shifti, shiftj;

/* Parse global variables of unit ratio */

  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
	
  long int iseed = -1;

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





/* Initialize the grid including the ghost cells.  */

	Real d0, u0, T0, x1, x2, x3, temperature, t, theta, Omegaimg, Omegareal;
	d0 = 1.0;
	u0 = 0.0;
	T0 = 1.0;
	t = pGrid->time;
	Real factor = 1.0;

	Omegareal = 0.0;
	Omegaimg = -9.50765;

	Real v0, flux, flag, radius;
	v0 = 0.0;
	flux = (4.0/3.0)*pow(T0,4.0)*v0/Crat;
	flag = 0.0;
	Gamma = 5.0/3.0;

	/* Incline the wave with 30 degree */

  for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	theta = - 2.0 * 3.1415926 * x1 + Omegareal * t;	
	radius = sqrt(x1*x1+x2*x2);

	temperature = T0;
        pGrid->U[k][j][i].d  = 1.0 ;
/*          pGrid->U[k][j][i].M1 = v0 + flag * 1.0 * (ran2(&iseed) - 0.5) * 1.e-2;
          pGrid->U[k][j][i].M2 = flag * 1.0 * (ran2(&iseed) - 0.5) * 1.e-2;
  
*/
	if(radius < 0.3 && radius > 0.28){
        	pGrid->U[k][j][i].M1 = v0 - x2;
          	pGrid->U[k][j][i].M2 =  x1;

	}
	else{
	 pGrid->U[k][j][i].M1 = v0;
	 pGrid->U[k][j][i].M2 = 0.0; 
	}


        pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 1.0/(Gamma - 1.0)+0.5*((pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d + (pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d);

	 pGrid->U[k][j][i].Edd_11 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	 pGrid->U[k][j][i].Edd_22 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	 pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
	 pGrid->U[k][j][i].Sigma_t = kappa *  pGrid->U[k][j][i].d;
	 pGrid->U[k][j][i].Sigma_a = kappa *  pGrid->U[k][j][i].d;



	  pGrid->U[k][j][i].Er = 1.0;
	  pGrid->U[k][j][i].Fr1 = flux;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;

	
	 	
        }
      }
    }
	dump_history_enroll(vorticity_p,"vor_p");
	dump_history_enroll(vorticity_n,"vor_n");

	dump_history_enroll(vorticityFr_p,"vorFr_p");
	dump_history_enroll(vorticityFr_n,"vorFr_n");

	
	dump_history_enroll(vortSFr_p,"vorFrS_p");
	dump_history_enroll(vortSFr_n,"vorFrS_n");
	

	dump_history_enroll(dErTv_p,"dErTv_p");
	dump_history_enroll(dErTv_n,"dErTv_n");

	dump_history_enroll(Radwork_p,"Radwork_p");
	dump_history_enroll(Radwork_n,"Radwork_n");

	dump_history_enroll(vortS_p,"vortS_p");
	dump_history_enroll(vortS_n,"vortS_n");

	dump_history_enroll(divV_p,"divV_p");
	dump_history_enroll(divV_n,"divV_n");

	dump_history_enroll(baronp,"baron_p");
	dump_history_enroll(baronn,"baron_n");

	dump_history_enroll(dTdotv,"dTdotv");



	Opacity = constopa;

/* data for radiation transfer method */
#ifdef RADIATION_TRANSFER
  RadGridS *pRG = (pDomain->RadGrid);

  int nf=pRG->nf, nang=pRG->nang;
  int ifr, l, m;
  
/* We do not need to iniatialize tau */

/* Read problem parameters. */

  eps0 = par_getd("problem","eps");

/* ------- Initialize boundary emission ---------------------------------- */

  for(ifr=0; ifr<nf; ifr++)
    for(l=0; l<4; l++) 
      for(m=0; m<nang; m++) {
	pRG->r1imu[ifr][pRG->ks][0][l][m] = Thermal_B(pGrid, ifr, pGrid->ie+1, nghost-1, pGrid->ks);
	pRG->l1imu[ifr][pRG->ks][0][l][m] = Thermal_B(pGrid, ifr, pGrid->is-1, nghost-1, pGrid->ks);
      }
  for(j=pRG->js; j<=pRG->je+1; j++) {
/* incident radiation at left boundary */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<nang; m++) {
	  pRG->l1imu[ifr][pRG->ks][j][0][m] = Thermal_B(pGrid, ifr, pGrid->is-1, j+nghost-1, pGrid->ks);
	  pRG->l1imu[ifr][pRG->ks][j][2][m] = Thermal_B(pGrid, ifr, pGrid->is-1, j+nghost-1, pGrid->ks);
      }
/* incident radiation at right boundary */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<=nang; m++) {
	  pRG->r1imu[ifr][pRG->ks][j][1][m] = Thermal_B(pGrid, ifr, pGrid->ie+1, j+nghost-1, pGrid->ks);
	  pRG->r1imu[ifr][pRG->ks][j][3][m] = Thermal_B(pGrid, ifr, pGrid->ie+1, j+nghost-1, pGrid->ks);
      }
  }

  for(i=pRG->is-1; i<=pRG->ie+1; i++) {
/* incident radiation at upper and lower boundaries */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<nang; m++) {
/* lower boundary is tau=0, no irradiation */
	pRG->l2imu[ifr][pRG->ks][i][0][m] = Thermal_B(pGrid, ifr, i+nghost-1, pGrid->js-1, pGrid->ks);
	pRG->l2imu[ifr][pRG->ks][i][1][m] = Thermal_B(pGrid, ifr, i+nghost-1, pGrid->js-1, pGrid->ks);
/* upper boundary is large tau, eps=1 */
	pRG->r2imu[ifr][pRG->ks][i][2][m] = Thermal_B(pGrid, ifr, i+nghost-1, pGrid->je+1, pGrid->ks);
	pRG->r2imu[ifr][pRG->ks][i][3][m] = Thermal_B(pGrid, ifr, i+nghost-1, pGrid->je+1, pGrid->ks);
      }
  }

/* enroll radiation specification functions */
get_thermal_source = Thermal_B;
get_thermal_fraction = const_eps;
get_total_opacity = transfer_opacity;

#endif




/* INItialize the Eddington tensor */
/* If initial data depends on Eddington tensor, we need to use the Eddington tensor calculated here */
#ifdef RADIATION_TRANSFER

	hydro_to_rad(pDomain);  
/* solve radiative transfer */
	formal_solution(pDomain);
	/* Get the Eddington tensor */
	Eddington_FUN(pGrid, pRG);


#endif


  return;
}



void constopa(const Real rho, const Real T, Real *Sigma_t, Real *Sigma_a, Real dSigma[4]){
	if(Sigma_t != NULL)
		*Sigma_t = kappa * rho;
	
	if(Sigma_a != NULL)
		*Sigma_a = *Sigma_t;

	if(dSigma != NULL){
		dSigma[0] = kappa;
		dSigma[1] = kappa;
		dSigma[2] = 0.0;
		dSigma[3] = 0.0;
	}
	

 return; 

}
/*
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
	T0 = 1.0 + 7.5;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);

 	for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
      		pGrid->U[ks][j][ie+i].d  = d0;
		if((j>=jl+shiftj)&&(j<=ju-shiftj))
		pGrid->U[ks][j][ie+i].M1 = d0 * u0;
		else 
		pGrid->U[ks][j][ie+i].M1 = 0.0;

	      	pGrid->U[ks][j][ie+i].E  = 0.5 * d0 * (pGrid->U[ks][j][ie+i].M1) * (pGrid->U[ks][j][ie+i].M1)  + d0 * T0/(Gamma - 1.0);
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

	Real d0, u0, T0;
	d0 = 1.0;
	u0 = -20.0;
	T0 = 1.0 + 7.5;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	shiftj = (int)((ju-jl)*2/5);

 	for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
      		pGrid->U[ks][j][ie+i].Er  = 8.8 * 8.5 * 8.5 * 8.5;
	      	pGrid->U[ks][j][ie+i].Fr1 = - 30.0 * 0.3333333 * 8.5 * 8.5 * 8.5 / Sigma_t;
    		}
	}



  
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
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
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



Real vorticity_p(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);
	
	if(vorticity > 0.0)
		return vorticity;
	else
		return 0.0; 
}

Real vorticity_n(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);
	
	if(vorticity < 0.0)
		return vorticity;
	else
		return 0.0; 
}


Real Ptemperature(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real pressure, temperature;

	pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
			+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2)/pG->U[k][j][i].d) * (Gamma - 1.0);

/* if MHD - 0.5 * Bx * Bx   */

#ifdef RADIATION_MHD

		pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

    		temperature = pressure / (pG->U[k][j][i].d * R_ideal);

	
	if(vorticity > 0.0)
		return temperature;
	else
		return 0.0; 
}


Real Ntemperature(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real pressure, temperature;

	pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
			+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2)/pG->U[k][j][i].d) * (Gamma - 1.0);

/* if MHD - 0.5 * Bx * Bx   */

#ifdef RADIATION_MHD

		pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

    		temperature = pressure / (pG->U[k][j][i].d * R_ideal);

	
	if(vorticity < 0.0)
		return temperature;
	else
		return 0.0; 
}

Real dTdotvp(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real pressure, temperaturex1, temperaturex2, temperaturey1, temperaturey2;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	temperaturex1 = pressure / (pG->U[k][j][i-1].d * R_ideal);

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);

    	temperaturex2 = pressure / (pG->U[k][j][i+1].d * R_ideal);


	
	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

    	temperaturey1 = pressure / (pG->U[k][j-1][i].d * R_ideal);

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);

    	temperaturey2 = pressure / (pG->U[k][j+1][i].d * R_ideal);

	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;

	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	result = vx * (temperaturex2 - temperaturex1)/(2.0*dx) + vy * (temperaturey2 - temperaturey1)/(2.0*dy); 
	

	
	if(vorticity > 0.0)
		return result;
	else
		return 0.0; 
}

Real dTdotvn(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real pressure, temperaturex1, temperaturex2, temperaturey1, temperaturey2;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	temperaturex1 = pressure / (pG->U[k][j][i-1].d * R_ideal);

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);

    	temperaturex2 = pressure / (pG->U[k][j][i+1].d * R_ideal);


	
	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

    	temperaturey1 = pressure / (pG->U[k][j-1][i].d * R_ideal);

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);

    	temperaturey2 = pressure / (pG->U[k][j+1][i].d * R_ideal);

	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;

	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	result = vx * (temperaturex2 - temperaturex1)/(2.0*dx) + vy * (temperaturey2 - temperaturey1)/(2.0*dy); 
	

	
	if(vorticity < 0.0)
		return result;
	else
		return 0.0; 
}



Real dTdotv(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real pressure, temperaturex1, temperaturex2, temperaturey1, temperaturey2;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	temperaturex1 = pressure / (pG->U[k][j][i-1].d * R_ideal);

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);

    	temperaturex2 = pressure / (pG->U[k][j][i+1].d * R_ideal);


	
	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

    	temperaturey1 = pressure / (pG->U[k][j-1][i].d * R_ideal);

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);

    	temperaturey2 = pressure / (pG->U[k][j+1][i].d * R_ideal);

	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;

	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	result = vx * (temperaturex2 - temperaturex1)/(2.0*dx) + vy * (temperaturey2 - temperaturey1)/(2.0*dy); 
	

	return result;
	
}




Real baronp(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real drhodx, drhody, dpdx, dpdy, pressurex1, pressurex2, pressurey1, pressurey2;

	pressurex1 = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	

	pressurex2 = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);



	
	pressurey1 = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

	pressurey2 = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);


	drhodx = (pG->U[k][j][i+1].d - pG->U[k][j][i-1].d)/(2.0*dx);
	drhody = (pG->U[k][j+1][i].d - pG->U[k][j-1][i].d)/(2.0*dy);

	dpdx  = (pressurex2 - pressurex1)/(2.0*dx);
	dpdy  = (pressurey2 - pressurey1)/(2.0*dy);

	result = (drhodx*dpdy - drhody*dpdx)/(pG->U[k][j][i].d * pG->U[k][j][i].d);


	
	if(vorticity > 0.0)
		return result;
	else
		return 0.0; 
}



Real baronn(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real drhodx, drhody, dpdx, dpdy, pressurex1, pressurex2, pressurey1, pressurey2;

	pressurex1 = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	

	pressurex2 = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);



	
	pressurey1 = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

	pressurey2 = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);


	drhodx = (pG->U[k][j][i+1].d - pG->U[k][j][i-1].d)/(2.0*dx);
	drhody = (pG->U[k][j+1][i].d - pG->U[k][j-1][i].d)/(2.0*dy);

	dpdx  = (pressurex2 - pressurex1)/(2.0*dx);
	dpdy  = (pressurey2 - pressurey1)/(2.0*dy);

	result = (drhodx*dpdy - drhody*dpdx)/(pG->U[k][j][i].d * pG->U[k][j][i].d);

	
	if(vorticity < 0.0)
		return result;
	else
		return 0.0; 
}





Real baron(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real drhodx, drhody, dpdx, dpdy, pressurex1, pressurex2, pressurey1, pressurey2;

	pressurex1 = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);

    	

	pressurex2 = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);



	
	pressurey1 = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);

	pressurey2 = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);


	drhodx = (pG->U[k][j][i+1].d - pG->U[k][j][i-1].d)/(2.0*dx);
	drhody = (pG->U[k][j+1][i].d - pG->U[k][j-1][i].d)/(2.0*dy);

	dpdx  = (pressurex2 - pressurex1)/(2.0*dx);
	dpdy  = (pressurey2 - pressurey1)/(2.0*dy);

	result = (drhodx*dpdy - drhody*dpdx)/(pG->U[k][j][i].d * pG->U[k][j][i].d);

	result *= vorticity;

	if(fabs(result) > 0.0) result = result /fabs(result);


	return result;
	
}




Real baronpEr(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real drhodx, drhody, dpdx, dpdy, pressurex1, pressurex2, pressurey1, pressurey2;

	pressurex1 =Prat* pG->U[k][j][i-1].Er/3.0;

    	pressurex2 =Prat* pG->U[k][j][i+1].Er/3.0;

	pressurey1 =Prat* pG->U[k][j-1][i].Er/3.0;

    	pressurey2 =Prat* pG->U[k][j+1][i].Er/3.0;


	drhodx = (pG->U[k][j][i+1].d - pG->U[k][j][i-1].d)/(2.0*dx);
	drhody = (pG->U[k][j+1][i].d - pG->U[k][j-1][i].d)/(2.0*dy);

	dpdx  = (pressurex2 - pressurex1)/(2.0*dx);
	dpdy  = (pressurey2 - pressurey1)/(2.0*dy);

	result = drhodx*dpdy - drhody*dpdx;


	
	if(vorticity > 0.0)
		return result;
	else
		return 0.0; 
}



Real baronnEr(const GridS *pG, const int i, const int j, const int k)
{  
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity, vx, vy, result;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0*dx) - (vx2 - vx1)/(2.0*dy);

	Real drhodx, drhody, dpdx, dpdy, pressurex1, pressurex2, pressurey1, pressurey2;

	pressurex1 =Prat* pG->U[k][j][i-1].Er/3.0;

    	pressurex2 =Prat* pG->U[k][j][i+1].Er/3.0;

	pressurey1 =Prat* pG->U[k][j-1][i].Er/3.0;

    	pressurey2 =Prat* pG->U[k][j+1][i].Er/3.0;


	drhodx = (pG->U[k][j][i+1].d - pG->U[k][j][i-1].d)/(2.0*dx);
	drhody = (pG->U[k][j+1][i].d - pG->U[k][j-1][i].d)/(2.0*dy);

	dpdx  = (pressurex2 - pressurex1)/(2.0*dx);
	dpdy  = (pressurey2 - pressurey1)/(2.0*dy);

	result = drhodx*dpdy - drhody*dpdx;


	
	if(vorticity < 0.0)
		return result;
	else
		return 0.0; 
}



Real vorticityFr_p(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);
	
	if(vorticity > 0.0)
		return vorticity;
	else
		return 0.0; 
}


Real vorticityFr_n(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);
	
	if(vorticity < 0.0)
		return vorticity;
	else
		return 0.0; 
}


Real vortSFr_p(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);




	Real pressurex1, temperaturex1, Erx1, vx1, pressurex2, temperaturex2, Erx2, vx2;
	Real pressurey1, temperaturey1, Ery1, vy1, pressurey2, temperaturey2, Ery2, vy2;

	pressurex1 = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);
	temperaturex1 = pressurex1 / (pG->U[k][j-1][i].d * R_ideal);

	vx1 = pG->U[k][j-1][i].M1 /pG->U[k][j-1][i].d;
	Erx1 = pG->U[k][j-1][i].Er;

	pressurex2 = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);
	temperaturex2 = pressurex2 / (pG->U[k][j+1][i].d * R_ideal);

	vx2 = pG->U[k][j+1][i].M1 /pG->U[k][j+1][i].d;
	Erx2 = pG->U[k][j+1][i].Er;


	pressurey1 = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);
	temperaturey1 = pressurey1 / (pG->U[k][j][i-1].d * R_ideal);

	vy1 = pG->U[k][j][i-1].M2 /pG->U[k][j][i-1].d;
	Ery1 = pG->U[k][j][i-1].Er;

	pressurey2 = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);
	temperaturey2 = pressurey2 / (pG->U[k][j][i+1].d * R_ideal);

	vy2 = pG->U[k][j][i+1].M2 /pG->U[k][j][i+1].d;
	Ery2 = pG->U[k][j][i+1].Er;
	

	Real Sourcevort;

	Sourcevort = (vy2 * (pow(temperaturey2, 4.0) + Ery2/3.0) - vy1 * (pow(temperaturey1, 4.0) + Ery1/3.0))/(2.0 * dx)
			- (vx2 * (pow(temperaturex2, 4.0) + Erx2/3.0) - vx1 * (pow(temperaturex1, 4.0) + Erx1/3.0))/(2.0 * dy);

	
	if(vorticity > 0.0)
		return (-Crat* pG->U[k][j][i].Sigma_a * vorticity + pG->U[k][j][i].Sigma_a * Sourcevort);
	else
		return 0.0; 
}



Real vortSFr_n(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);




	Real pressurex1, temperaturex1, Erx1, vx1, pressurex2, temperaturex2, Erx2, vx2;
	Real pressurey1, temperaturey1, Ery1, vy1, pressurey2, temperaturey2, Ery2, vy2;

	pressurex1 = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);
	temperaturex1 = pressurex1 / (pG->U[k][j-1][i].d * R_ideal);

	vx1 = pG->U[k][j-1][i].M1 /pG->U[k][j-1][i].d;
	Erx1 = pG->U[k][j-1][i].Er;

	pressurex2 = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);
	temperaturex2 = pressurex2 / (pG->U[k][j+1][i].d * R_ideal);

	vx2 = pG->U[k][j+1][i].M1 /pG->U[k][j+1][i].d;
	Erx2 = pG->U[k][j+1][i].Er;


	pressurey1 = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);
	temperaturey1 = pressurey1 / (pG->U[k][j][i-1].d * R_ideal);

	vy1 = pG->U[k][j][i-1].M2 /pG->U[k][j][i-1].d;
	Ery1 = pG->U[k][j][i-1].Er;

	pressurey2 = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);
	temperaturey2 = pressurey2 / (pG->U[k][j][i+1].d * R_ideal);

	vy2 = pG->U[k][j][i+1].M2 /pG->U[k][j][i+1].d;
	Ery2 = pG->U[k][j][i+1].Er;
	

	Real Sourcevort;

	Sourcevort = (vy2 * (pow(temperaturey2, 4.0) + Ery2/3.0) - vy1 * (pow(temperaturey1, 4.0) + Ery1/3.0))/(2.0 * dx)
			- (vx2 * (pow(temperaturex2, 4.0) + Erx2/3.0) - vx1 * (pow(temperaturex1, 4.0) + Erx1/3.0))/(2.0 * dy);

	
	if(vorticity < 0.0)
		return (-Crat* pG->U[k][j][i].Sigma_a * vorticity + pG->U[k][j][i].Sigma_a * Sourcevort);
	else
		return 0.0; 
}





Real dErTv_p(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);

	Real pressure, temperature1, Er1, temperature2, Er2, dTdx, dTdy, vx, vy, result;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);	
	temperature1 = pressure / pG->U[k][j][i-1].d;

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);	
	temperature2 = pressure / pG->U[k][j][i+1].d;

	Er1 = pG->U[k][j][i-1].Er;
	Er2 = pG->U[k][j][i+1].Er;

	dTdx = ((pow(temperature2, 4.0) + Er2/3.0) - (pow(temperature1, 4.0) + Er1/3.0))/(2.0 * dx);


	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);	
	temperature1 = pressure / pG->U[k][j-1][i].d;

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);	
	temperature2 = pressure / pG->U[k][j+1][i].d;

	Er1 = pG->U[k][j-1][i].Er;
	Er2 = pG->U[k][j+1][i].Er;

	dTdy = ((pow(temperature2, 4.0) + Er2/3.0) - (pow(temperature1, 4.0) + Er1/3.0))/(2.0 * dy);


	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	result = dTdx * vy - dTdy * vx;

	
	if(vorticity > 0.0)
		return result;
	else
		return 0.0; 
}


Real dErTv_n(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);

	Real pressure, temperature1, Er1, temperature2, Er2, dTdx, dTdy, vx, vy, result;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d) * (Gamma - 1.0);	
	temperature1 = pressure / pG->U[k][j][i-1].d;

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/pG->U[k][j][i+1].d) * (Gamma - 1.0);	
	temperature2 = pressure / pG->U[k][j][i+1].d;

	Er1 = pG->U[k][j][i-1].Er;
	Er2 = pG->U[k][j][i+1].Er;

	dTdx = ((pow(temperature2, 4.0) + Er2/3.0) - (pow(temperature1, 4.0) + Er1/3.0))/(2.0 * dx);


	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d) * (Gamma - 1.0);	
	temperature1 = pressure / pG->U[k][j-1][i].d;

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/pG->U[k][j+1][i].d) * (Gamma - 1.0);	
	temperature2 = pressure / pG->U[k][j+1][i].d;

	Er1 = pG->U[k][j-1][i].Er;
	Er2 = pG->U[k][j+1][i].Er;

	dTdy = ((pow(temperature2, 4.0) + Er2/3.0) - (pow(temperature1, 4.0) + Er1/3.0))/(2.0 * dy);


	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	result = dTdx * vy - dTdy * vx;

	
	if(vorticity < 0.0)
		return result;
	else
		return 0.0; 
}



Real Radwork_p(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);

	Real vx, vy, Er, Fx, Fy, result;

	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	Fx = pG->U[k][j][i].Fr1;
	Fy = pG->U[k][j][i].Fr2;

	Er = pG->U[k][j][i].Er;
	
	result = vx * (Fx - 4.0 * vx * Er / (3.0 * Crat)) + vy * (Fy - 4.0 * vy * Er/(3.0 * Crat));
	
	if(vorticity > 0.0)
		return result;
	else
		return 0.0; 
}

Real Radwork_n(const GridS *pG, const int i, const int j, const int k)
{  
	Real Frx1, Frx2, Fry1, Fry2, dx, dy, vorticity;
	Frx1 = pG->U[k][j-1][i].Fr1;
	Frx2 = pG->U[k][j+1][i].Fr1;
	Fry1 = pG->U[k][j][i-1].Fr2;
	Fry2 = pG->U[k][j][i+1].Fr2;
	dx = pG->dx1;
	dy = pG->dx2;
	
	
	vorticity = (Fry2 - Fry1)/(2.0 * dx) - (Frx2 - Frx1)/(2.0 * dy);

	Real vx, vy, Er, Fx, Fy, result;

	vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
	vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;

	Fx = pG->U[k][j][i].Fr1;
	Fy = pG->U[k][j][i].Fr2;

	Er = pG->U[k][j][i].Er;
	
	result = vx * (Fx - 4.0 * vx * Er / (3.0 * Crat)) + vy * (Fy - 4.0 * vy * Er/(3.0 * Crat));
	
	if(vorticity < 0.0)
		return result;
	else
		return 0.0; 
}



Real vortS_p(const GridS *pG, const int i, const int j, const int k)
{  
	
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0 * dx) - (vx2 - vx1)/(2.0 * dy);


	Real pressure, temperature, Fr, v, d, Er, result, Sigma;
	Real Sx1, Sx2, Sy1, Sy2;

	d = pG->U[k][j-1][i].d;
	v = pG->U[k][j-1][i].M1/d;
	Fr = pG->U[k][j-1][i].Fr1;
	Er = pG->U[k][j-1][i].Er;
	Sigma = pG->U[k][j-1][i].Sigma_a;

	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sx1 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j+1][i].d;
	v = pG->U[k][j+1][i].M1/d;
	Fr = pG->U[k][j+1][i].Fr1;
	Er = pG->U[k][j+1][i].Er;
	Sigma = pG->U[k][j+1][i].Sigma_a;

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sx2 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j][i-1].d;
	v = pG->U[k][j][i-1].M2/d;
	Fr = pG->U[k][j][i-1].Fr2;
	Er = pG->U[k][j][i-1].Er;
	Sigma = pG->U[k][j][i-1].Sigma_a;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sy1 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j][i+1].d;
	v = pG->U[k][j][i+1].M2/d;
	Fr = pG->U[k][j][i+1].Fr2;
	Er = pG->U[k][j][i+1].Er;
	Sigma = pG->U[k][j][i+1].Sigma_a;

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sy2 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	result = (Sy2 - Sy1)/(2.0 * dx) - (Sx2 - Sx1)/(2.0 * dy);


	
	if(vorticity > 0.0)
		return -Prat * result;
	else
		return 0.0; 
}



Real vortS_n(const GridS *pG, const int i, const int j, const int k)
{  
	
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0 * dx) - (vx2 - vx1)/(2.0 * dy);


	Real pressure, temperature, Fr, v, d, Er, result, Sigma;
	Real Sx1, Sx2, Sy1, Sy2;

	d = pG->U[k][j-1][i].d;
	v = pG->U[k][j-1][i].M1/d;
	Fr = pG->U[k][j-1][i].Fr1;
	Er = pG->U[k][j-1][i].Er;
	Sigma = pG->U[k][j-1][i].Sigma_a;

	pressure = (pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 
			+ pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sx1 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j+1][i].d;
	v = pG->U[k][j+1][i].M1/d;
	Fr = pG->U[k][j+1][i].Fr1;
	Er = pG->U[k][j+1][i].Er;
	Sigma = pG->U[k][j+1][i].Sigma_a;

	pressure = (pG->U[k][j+1][i].E - 0.5 * (pG->U[k][j+1][i].M1 * pG->U[k][j+1][i].M1 
			+ pG->U[k][j+1][i].M2 * pG->U[k][j+1][i].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sx2 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j][i-1].d;
	v = pG->U[k][j][i-1].M2/d;
	Fr = pG->U[k][j][i-1].Fr2;
	Er = pG->U[k][j][i-1].Er;
	Sigma = pG->U[k][j][i-1].Sigma_a;

	pressure = (pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 
			+ pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sy1 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	d = pG->U[k][j][i+1].d;
	v = pG->U[k][j][i+1].M2/d;
	Fr = pG->U[k][j][i+1].Fr2;
	Er = pG->U[k][j][i+1].Er;
	Sigma = pG->U[k][j][i+1].Sigma_a;

	pressure = (pG->U[k][j][i+1].E - 0.5 * (pG->U[k][j][i+1].M1 * pG->U[k][j][i+1].M1 
			+ pG->U[k][j][i+1].M2 * pG->U[k][j][i+1].M2)/d) * (Gamma - 1.0);
	temperature = pressure / (d * R_ideal);

	Sy2 = (-Sigma * (Fr - 4.0 * v * Er/(3.0 * Crat)) + Sigma * v * (pow(temperature, 4.0) - Er)/Crat) / d;

	result = (Sy2 - Sy1)/(2.0 * dx) - (Sx2 - Sx1)/(2.0 * dy);


	
	if(vorticity < 0.0)
		return -Prat * result;
	else
		return 0.0; 
}



Real divV_p(const GridS *pG, const int i, const int j, const int k)
{  
	
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0 * dx) - (vx2 - vx1)/(2.0 * dy);

	vx1 = pG->U[k][j][i-1].M1 / pG->U[k][j][i-1].d;
	vx2 = pG->U[k][j][i+1].M1 / pG->U[k][j][i+1].d;
	vy1 = pG->U[k][j-1][i].M2 / pG->U[k][j-1][i].d;
	vy2 = pG->U[k][j+1][i].M2 / pG->U[k][j+1][i].d;

	Real result;

	result = -((vx2 - vx1)/(2.0 * dx) + (vy2 - vy1)/(2.0 * dy));

	
	if(vorticity > 0.0)
		return result * vorticity;
	else
		return 0.0; 
	
}

Real divV_n(const GridS *pG, const int i, const int j, const int k)
{  
	
	Real vx1, vx2, vy1, vy2, dx, dy, vorticity;
	vx1 = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;
	vx2 = pG->U[k][j+1][i].M1 / pG->U[k][j+1][i].d;
	vy1 = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
	vy2 = pG->U[k][j][i+1].M2 / pG->U[k][j][i+1].d;
	dx = pG->dx1;
	dy = pG->dx2;
	
	vorticity = (vy2 - vy1)/(2.0 * dx) - (vx2 - vx1)/(2.0 * dy);

	vx1 = pG->U[k][j][i-1].M1 / pG->U[k][j][i-1].d;
	vx2 = pG->U[k][j][i+1].M1 / pG->U[k][j][i+1].d;
	vy1 = pG->U[k][j-1][i].M2 / pG->U[k][j-1][i].d;
	vy2 = pG->U[k][j+1][i].M2 / pG->U[k][j+1][i].d;

	Real result;

	result = -((vx2 - vx1)/(2.0 * dx) + (vy2 - vy1)/(2.0 * dy));

	
	if(vorticity < 0.0)
		return result * vorticity;
	else
		return 0.0; 
	
}




/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


/* Function for transfer module */
#ifdef RADIATION_TRANSFER

static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
	Real density, vx, vy, vz, energy, pressure, T, B;
	density = pG->U[k][j][i].d;
	vx = pG->U[k][j][i].M1/density;
	vy = pG->U[k][j][i].M2/density;
	vz = pG->U[k][j][i].M3/density;
	energy = pG->U[k][j][i].E;
	pressure = (energy - 0.5 * density * (vx* vx + vy * vy + vz * vz)) * (Gamma - 1.0);
#ifdef RADIATION_MHD

	pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c) * (Gamma - 1.0);
#endif
	T = pressure / (density * R_ideal);
	B = T * T * T * T / (4.0 * PI);

  return B;
}

static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k)
{

  return eps0;
  
}

static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k)
{

  return (pG->U[k][j][i].d * kappa);
  
}

#endif
