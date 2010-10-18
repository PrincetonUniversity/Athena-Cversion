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

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  int shift;

/* Parse global variables of unit ratio */
#ifdef rad_hydro
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  Sigma_t = par_getd("problem","Sigma_t");
  Sigma_a = par_getd("problem","Sigma_a");
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
	Real d0, u0, T0, x1, x2, x3, temperature;
	d0 = 1.0;
	u0 = 0.0;
	T0 = 1.0;


    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	
	  temperature = T0;
          pGrid->U[k][j][i].d  = d0;
          pGrid->U[k][j][i].M1 = d0 * u0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 0.5 * d0 * u0 * u0 + d0 * temperature /(Gamma - 1.0);

	 pGrid->U[k][j][i].Edd_11 = 0.33333333; /* Set to be a constant in 1D. To be modified later */

#ifdef MHD
          pGrid->B1i[k][j][i] = 0.0;
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;
          pGrid->U[k][j][i].B1c = 0.0;
          pGrid->U[k][j][i].B2c = 0.0;
          pGrid->U[k][j][i].B3c = 0.0;
#endif
#ifdef rad_hydro
	  pGrid->U[k][j][i].Er = temperature * temperature * temperature * temperature ;
	  pGrid->U[k][j][i].Fr1 = 0.0;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;

	 		
#endif
        }
      }
    }

	bvals_mhd_fun(pDomain, left_x1, radMHD_inflow);
	bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow);

  return;
}

void radMHD_inflow(GridS *pGrid)
{
  	int i, is;
	int ks, js;
	is = pGrid->is;
  	ks = pGrid->ks;
	js = pGrid->js;

	double t=pGrid->time;
	

	Real vamp, omega;
	Real temp;
	
	vamp = 1.e-3;
	omega = 2.0 * 3.1415926;
 
    for (i=1;  i<=nghost;  i++) {
      pGrid->U[ks][js][is-i].d  = 1.0;
      pGrid->U[ks][js][is-i].M1 =vamp * cos(omega * t);
      temp = cos(omega*t);
      pGrid->U[ks][js][is-i].E  = 0.5 * vamp * vamp * cos(omega * t) * cos(omega * t) + 1.0/(Gamma - 1.0);
    }
  
}

void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, is;
	int ks, js;
	is = pGrid->is;
  	ks = pGrid->ks;
	js = pGrid->js;

	double t=pGrid->time;
	

	Real Eramp, omega;
	
	Eramp = 1.e-3;
	omega = 2.0 * 3.1415926;
 
    for (i=1;  i<=nghost;  i++) {
      pGrid->U[ks][js][is-i].Er  = 1.0 + Eramp * cos(omega * t);
      pGrid->U[ks][js][is-i].Fr1 = 0.0;
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
	fprintf(fp,"%5.3e\n",Sigma_t);
	fprintf(fp,"%5.3e\n",Sigma_a);
	fprintf(fp,"%5.3e\n",R_ideal);
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{

	bvals_mhd_fun(&(pM->Domain[0][0]), right_x1, radMHD_inflow);
	fscanf(fp,"%lf",&Gamma);
	fscanf(fp,"%lf\n",&Prat);
	fscanf(fp,"%lf\n",&Crat);
	fscanf(fp,"%lf\n",&Sigma_t);
	fscanf(fp,"%lf\n",&Sigma_a);
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
