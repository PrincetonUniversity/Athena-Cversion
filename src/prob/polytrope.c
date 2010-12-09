#include "copyright.h"
/*============================================================================*/
/*! \file polytrope.c
 *  \brief Problem generator for a self-gravitating polytrope.
 *
 * PURPOSE: Problem generator for a self-gravitating polytrope.
 *
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

double theta_interpol(double rad, double dr, double* theta, int nelem)
{
  int i = rad/dr;
  if(i+1 >= nelem) 
    ath_error("out of bounds array access in theta_interpol \n");
  Real delta = rad/dr - i;
  return (1-delta)*theta[i] + delta*theta[i+1];
}

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real pressure,rad,pa,da,x1,x2,x3;
  double npoly;
  Real b0=0.0,Bx=0.0,d_c,p_c;
  //step in xi, radius of star in xi units, conversion xi to r
  Real d_xi, dr, xi_star, xi_last; 
  const int buffersize = 100000;
  Real theta[buffersize];

  //read in variables from input file
  Real rstar = par_getd("problem","rstar");
  pa  = par_getd("problem","pamb"); //pa is a fraction of p_c
  da  = par_getd("problem","damb"); //da is a fraction of d_c
  d_c  = par_getd("problem","d_c");  
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = par_getd("problem","grav_mean_rho");
  //char* input_poly = "polytrope.txt";
  char* input_poly = par_gets("problem","input_poly");
  Real p_frac = par_getd("problem","p_frac"); //fraction of hydrostatic pres.
  Real u0 = par_getd("problem","u0"); //velocity of sphere across grid
  
  printf("rstar: %lf \n pa: %lf \n da: %lf \n npoly: %lf \n d_c: %lf \n 4piG: %lf \n d0: %lf \n input_poly: %s, p_frac: %lf, u0: %lf \n", rstar, pa, da, npoly, d_c, four_pi_G, grav_mean_rho, input_poly, p_frac, u0);

  //read in polytrope data
  FILE *infile;
  int nelem = 0;
  char dummy_s[80]; //dummy variable which is a string
  double dummy_d; //dummy variable which is a double
  infile = fopen(input_poly,"r");
  if (infile == NULL) 
    ath_error("problem reading input file for polytrope \n");
  fscanf(infile, "%s", dummy_s);
  fscanf(infile, "%s", dummy_s);
  fscanf(infile, "%s", dummy_s);
  fscanf(infile, "%lf", &npoly);
  fscanf(infile, "%s", dummy_s);
  fscanf(infile, "%lf", &d_xi);
  fscanf(infile, "%s", dummy_s);
  fscanf(infile, "%s", dummy_s);
  printf("last: %s \n", dummy_s);
  while(1) { //while true
    fscanf(infile, "%lf", &dummy_d);
    if(feof(infile)) break;
    if(nelem >= buffersize)
      ath_error("maximum buffer size for polytrope array exceeded \n");
    fscanf(infile, "%lf", &(theta[nelem]));
    xi_last = dummy_d;
    nelem ++;
  }
  
  //set d_xi, xi_star
  if(nelem < 2)
    ath_error("error polytrope file should contain more than one entry \n");
  if(!(theta[nelem-1] <= 0 && theta[nelem-2] > 0))
    ath_error("error format for theta data is incorrect. Should be theta[nelem-1] <= 0 && theta[nelem-2] > 0 \n");
  xi_star = xi_last - d_xi*theta[nelem-1]/(theta[nelem-1] - theta[nelem-2]);
  dr = d_xi*rstar/xi_star;

  //set the value of the central pressure
  p_c = p_frac*four_pi_G*d_c*d_c*dr*dr/(d_xi*d_xi*(npoly+1));

  //set the ambient pressure and density to be fractions of the central pressure and density
  pa*=p_c;
  da*=d_c;

  printf("d_xi: %lf, xi_star: %lf, p_c: %lf, pa: %lf, da: %lf \n", d_xi, xi_star, p_c, pa, da);

/* setup uniform ambient medium with central polytropic sphere. 
   Velocities are everywhere initialized to zero to simplify energy
   and momentum calculations.*/

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2);
	pGrid->U[k][j][i].d = da;
	if (rad < rstar) 
	  pGrid->U[k][j][i].d += d_c*pow(theta_interpol(rad,dr,theta,nelem),npoly);
	pGrid->U[k][j][i].M1 = u0*(pGrid->U[k][j][i].d);
	pGrid->U[k][j][i].M2 = 0;
	pGrid->U[k][j][i].M3 = 0;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = pa/(Gamma-1) + .5*u0*u0*(pGrid->U[k][j][i].d);
	if (rad < rstar) 
	  pGrid->U[k][j][i].E += p_c/(Gamma-1)*pow(theta_interpol(rad,dr,theta,nelem),npoly+1.);
#endif
 
#ifdef MHD
	pGrid->B1i[k][j][i] = 0.;
	pGrid->B2i[k][j][i] = 0.;
	pGrid->B3i[k][j][i] = 0.;
	pGrid->U[k][j][i].B1c = 0.;
	pGrid->U[k][j][i].B2c = 0.;
	pGrid->U[k][j][i].B3c = 0.;
#endif /* MHD */
      }
    }
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
}

void Userwork_after_loop(MeshS *pM)
{
}
