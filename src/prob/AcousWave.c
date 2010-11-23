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


void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  int shift;

/* Parse global variables of unit ratio */
#ifdef RADIATION_HYDRO
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
	
#endif

/* Set up the index bounds for initializing the grid */
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

/* Initialize the grid including the ghost cells.  */
	Real d0, u0, T0, x1, x2, x3, temperature, t, theta, Omegaimg, Omegareal;
	d0 = 1.0;
	u0 = 0.0;
	T0 = 1.0;
	t = pGrid->time;
	Real flag = 1.0;
	Real factor = 1.e-3;

	Omegareal = 62.4774; 
	Omegaimg = 6.66773;


    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	theta = -10.0 * 2 * 3.1415926  * x1 + Omegareal * t;	

	  temperature = T0;
          pGrid->U[k][j][i].d  = 1.0 + flag * factor * (1.e-3 * cos(theta) + 0.0000 * sin(theta));
          pGrid->U[k][j][i].M1 = flag * factor * (9.94358e-4 * cos(theta) - 1.0612e-4 * sin(theta));
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 1.0/(Gamma - 1.0) + flag * factor * (1.5e-3 * cos(theta) - 2.3764e-8 * sin(theta));

	 pGrid->U[k][j][i].Edd_11 = 0.33333333; /* Set to be a constant in 1D. To be modified later */
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
#ifdef RADIATION_HYDRO
	  pGrid->U[k][j][i].Er = 1.0 + flag * factor * (-6.75337e-9 * cos(theta) - 6.33082e-8 * sin(theta));
	  pGrid->U[k][j][i].Fr1 = flag * factor *  (-1.1287e-11 * cos(theta) - 5.16247e-12 * sin(theta));
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;
	
#endif
        }
      }
    }
/*
	bvals_mhd_fun(pDomain, left_x1, radMHD_inflow);
	bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow);
	bvals_mhd_fun(pDomain, right_x1, radMHD_inflow2);
	bvals_rad_fun(pDomain, right_x1, radMHD_rad_inflow2);
*/
  return;
}

void radMHD_inflow(GridS *pGrid)
{
  	int i, is;
	int ks, js;
	is = pGrid->is;
  	ks = pGrid->ks;
	js = pGrid->js;

	Real t, x1, x2, x3,theta, Kimg;
	t = pGrid->time;
	int zones=100;
	Kimg = -0.00920423;
	Real factor = 0.00001;
 
    for (i=1;  i<=nghost+zones;  i++) {
	cc_pos(pGrid, is-i+zones, js,ks, &x1, &x2, &x3);
	theta = -0.409717 * x1 + 1.0 * t;

      pGrid->U[ks][js][is-i+zones].d  = 1.0 + exp(Kimg * x1) * factor * (0.001 * cos(theta) + 0.0000 * sin(theta));
      pGrid->U[ks][js][is-i+zones].M1 = exp(Kimg * x1) * factor * (0.00243948 * cos(theta) - 0.0000548025 * sin(theta));
      pGrid->U[ks][js][is-i+zones].E  = 1.0/(Gamma - 1.0) + factor * exp(Kimg * x1) * (0.00201782 * cos(theta) - 0.0000279815 * sin(theta));
    }
  
}

void radMHD_inflow2(GridS *pGrid)
{
  	int i, ie;
	int ks, js;
	ie = pGrid->ie;
  	ks = pGrid->ks;
	js = pGrid->js;

	Real t, x1, x2, x3,theta, Kimg;
	t = pGrid->time;
	int zones=0;
	Kimg = -0.01252;
 
    for (i=1;  i<=nghost+zones;  i++) {
	cc_pos(pGrid, ie+i-zones, js,ks, &x1, &x2, &x3);
	theta = -0.999993 * x1 + 1.0 * t;

      pGrid->U[ks][js][ie+i-zones].d  = 1.0 + exp(Kimg * x1) * (0.001 * cos(theta) + 0.0000 * sin(theta));
      pGrid->U[ks][js][ie+i-zones].M1 = exp(Kimg * x1) * (0.00100001 * cos(theta) - 1.25202e-6 * sin(theta));
      pGrid->U[ks][js][ie+i-zones].E  = 1.0/(Gamma - 1.0) + exp(Kimg * x1) * (0.001500001 * cos(theta) - 3.75257e-6 * sin(theta));
    }
  
}



void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, is;
	int ks, js;
	is = pGrid->is;
  	ks = pGrid->ks;
	js = pGrid->js;

	Real t, x1, x2, x3,theta, Kimg,dt;
	t = pGrid->time;
	dt = pGrid->time;
	int zones=100;
	Kimg = -0.00920423;
	Real factor = 0.001;

	for (i=1;  i<=nghost+zones;  i++) {
		cc_pos(pGrid, is-i+zones, js,ks, &x1, &x2, &x3);
		theta = -0.409717 * x1 + 1.0 * (t + dt);
	      	pGrid->U[ks][js][is-i+zones].Er  = 1.0 + factor * exp(Kimg * x1) * (0.00138085 * cos(theta) - 0.0000746174 * sin(theta));
		pGrid->U[ks][js][is-i+zones].Sigma_t = 10000;
		pGrid->U[ks][js][is-i+zones].Sigma_a = 10000;
		pGrid->U[ks][js][is-i+zones].Edd_11 = 0.33333;
      		pGrid->U[ks][js][is-i+zones].Fr1 =factor *  exp(Kimg * x1) * (3.24668e-7 * cos(theta) - 2.61885e-8 * sin(theta));
    }
}




void radMHD_rad_inflow2(GridS *pGrid)
{
  	int i, ie;
	int ks, js;
	ie = pGrid->ie;
  	ks = pGrid->ks;
	js = pGrid->js;

	Real t, x1, x2, x3,theta, Kimg,dt;
	t = pGrid->time;
	dt = pGrid->time;
	int zones=0;
	Kimg = -0.001252;

	for (i=1;  i<=nghost+zones;  i++) {
		cc_pos(pGrid, ie+i-zones, js,ks, &x1, &x2, &x3);
		theta = -0.999993 * x1 + 1.0 * (t + dt);
	      	pGrid->U[ks][js][ie+i-zones].Er  = 1.0 + exp(Kimg * x1) * (-3.62496e-11 * cos(theta) - 7.00009e-9 * sin(theta));
		pGrid->U[ks][js][ie+i-zones].Sigma_t = 0.01;
		pGrid->U[ks][js][ie+i-zones].Sigma_a = 0.01;
		pGrid->U[ks][js][ie+i-zones].Edd_11 = 0.33333;
      		pGrid->U[ks][js][ie+i-zones].Fr1 = exp(Kimg * x1) * (-9.99996e-8 * cos(theta) - 2.50759e-10 * sin(theta));
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
