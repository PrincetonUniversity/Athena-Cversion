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

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static ConsS ***RootSoln=NULL;



void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
    int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3,Nx1,Nx2;
  int shift;

/* Parse global variables of unit ratio */
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
#endif
 is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  Nx1 = pDomain->Nx[0];
  Nx2 = pDomain->Nx[1];

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie+nghost;
  il = pGrid->is-nghost;

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

  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }


/* Initialize the grid including the ghost cells.  */
	Real d0, u0, T0, x1, x2, x3, temperature, t, theta, Omegareal,Omegaimg, rho0, P0, E0, B0, angle, thetaB, tlim;
	d0 = 1.0;
	u0 = 0.0;
	T0 = 1.0;
	rho0 = 1.0;
	P0 = 1.0;
	E0 = P0/(Gamma - 1.0)+d0*u0*u0/2.0;
#ifdef RADIATION_MHD
	B0 = sqrt(10.0/3.0);
	E0 += B0*B0/2.0;
	angle=45.0*PI/180.0;
#endif

	t = pGrid->time;
	Real flag = 1.0;
	Real factor = 1.e-3;

	Omegareal = 7.94937;
	Omegaimg = 0.587018;

	tlim = 2.0 * PI / Omegareal;


    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	theta = - 2.0 * PI * x1  + Omegareal * t;		

	  temperature = T0;
          pGrid->U[k][j][i].d  = rho0 + flag * factor * (1.0e-3 * cos(theta) + 0.0000 * sin(theta));
          pGrid->U[k][j][i].M1 = flag * factor * 1.0 * (1.2651809667779195e-3 * cos(theta) - 9.342678338278762e-5 * sin(theta));
          pGrid->U[k][j][i].M2 = flag * factor * 0.0 * (1.2651809667779195e-3 * cos(theta) - 9.342678338278762e-5 * sin(theta));
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = E0 + flag * factor * (2.387826854962222e-3 * cos(theta) - 3.545635964979445e-4 * sin(theta));

	 pGrid->U[k][j][i].Edd_11 = 1.0/3.0; /* Set to be a constant in 1D. To be modified later */
	 pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
	 pGrid->U[k][j][i].Sigma_t = 10.0;
	 pGrid->U[k][j][i].Sigma_a = 10.0;

#ifdef RADIATION_MHD
	 /* interface magnetic field */
/*	  thetaB = -2.0 * 3.1415926  * (x1 - pGrid->dx1/2.0) + Omegareal * t;	
*/
	  thetaB = theta;
          pGrid->B1i[k][j][i] = B0*cos(angle);
          pGrid->B2i[k][j][i] = B0 * sin(angle) + flag*factor*(-7.97878e-4*cos(thetaB)+1.98765e-8*sin(thetaB));
          pGrid->B3i[k][j][i] = 0.0;
        
#endif
	  pGrid->U[k][j][i].Er = 1.0 + flag * factor * (2.0922674468559182e-3 * cos(theta) - 8.35518866314462e-4 * sin(theta));
	  pGrid->U[k][j][i].Fr1 = flag * factor * 1.0 * (-1.7478816069974772e-4 * cos(theta) - 4.3823239331717244e-4 * sin(theta));
	  pGrid->U[k][j][i].Fr2 = flag * factor * 0.0 * (-1.7478816069974772e-4 * cos(theta) - 4.3823239331717244e-4 * sin(theta));
	  pGrid->U[k][j][i].Fr3 = 0.0;	
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

/* save solution on root grid */

  if (pDomain->Level == 0) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      RootSoln[k][j][i].d  = exp(-Omegaimg * tlim) * (pGrid->U[k][j][i].d - d0) + d0;
      RootSoln[k][j][i].M1 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].M1;
      RootSoln[k][j][i].M2 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].M2;
      RootSoln[k][j][i].M3 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].M3;
#ifndef ISOTHERMAL
      RootSoln[k][j][i].E  = exp(-Omegaimg * tlim) * (pGrid->U[k][j][i].E - E0) + E0;
#endif /* ISOTHERMAL */
#ifdef RADIATION_MHD
      RootSoln[k][j][i].B1c = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].B1c;
      RootSoln[k][j][i].B2c = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].B2c;
      RootSoln[k][j][i].B3c = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].B3c;
#endif
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
      RootSoln[k][j][i].Er =  exp(-Omegaimg * tlim) * (pGrid->U[k][j][i].Er - 1.0) + 1.0;
      RootSoln[k][j][i].Fr1 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].Fr1;
      RootSoln[k][j][i].Fr2 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].Fr2;
      RootSoln[k][j][i].Fr3 = exp(-Omegaimg * tlim) * pGrid->U[k][j][i].Fr3;
#endif


    }}}
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
  GridS *pGrid;
  int i=0,j=0,k=0;
 

  int is,ie,js,je,ks,ke;
  Real rms_error=0.0, t;
  ConsS error,total_error;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count;
#if defined MPI_PARALLEL
  double err[12+NSCALARS], tot_err[12+NSCALARS];
  int ierr,myID;
#endif

  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
#ifdef RADIATION_MHD
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;
#endif /* MHD */
  total_error.Er = 0.0;
  total_error.Fr1 = 0.0;
  total_error.Fr2 = 0.0;
  total_error.Fr3 = 0.0;
#ifndef ISOTHERMAL
  total_error.E = 0.0;
#endif /* ISOTHERMAL */


/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

  t = pGrid->time;

/* compute L1 error in each variable, and rms total error */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
#ifdef RADIATION_MHD
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
    error.E = 0.0;
#endif /* ISOTHERMAL */
    error.Er = 0.0;
    error.Fr1 = 0.0;
    error.Fr2 = 0.0;
    error.Fr3 = 0.0;


    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[k][j][i].d   - RootSoln[k][j][i].d );
      error.M1  += fabs(pGrid->U[k][j][i].M1  - RootSoln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - RootSoln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
#ifdef RADIATION_MHD
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[k][j][i].E   - RootSoln[k][j][i].E );
#endif /* ISOTHERMAL */
      error.Er += fabs(pGrid->U[k][j][i].Er -   RootSoln[k][j][i].Er);
      error.Fr1 += fabs(pGrid->U[k][j][i].Fr1 - RootSoln[k][j][i].Fr1);
      error.Fr2 += fabs(pGrid->U[k][j][i].Fr2 - RootSoln[k][j][i].Fr2);
      error.Fr3 += fabs(pGrid->U[k][j][i].Fr3 - RootSoln[k][j][i].Fr3);
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
#ifdef RADIATION_MHD
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E += error.E;
#endif /* ISOTHERMAL */
    total_error.Er += error.Er;
    total_error.Fr1 += error.Fr1;
    total_error.Fr2 += error.Fr2;
    total_error.Fr3 += error.Fr3;

  }}

#ifdef MPI_PARALLEL
  Nx1 = pM->Domain[0][0].Nx[0];
  Nx2 = pM->Domain[0][0].Nx[1];
  Nx3 = pM->Domain[0][0].Nx[2];
#else
  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;
#endif
  count = Nx1*Nx2*Nx3;

#ifdef MPI_PARALLEL 
/* Now we have to use an All_Reduce to get the total error over all the MPI
 * grids.  Begin by copying the error into the err[] array */

  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
#ifdef RADIATION_MHD
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */
  err[8] = total_error.Er;
  err[9] = total_error.Fr1;
  err[10] = total_error.Fr2;
  err[11] = total_error.Fr3;

/* Sum up the Computed Error */
  ierr = MPI_Reduce(err,tot_err,(8+NSCALARS),MPI_DOUBLE,MPI_SUM,0,
    pM->Domain[0][0].Comm_Domain);

/* If I'm the parent, copy the sum back to the total_error variable */

  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
#ifdef RADIATION_MHD
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = err[8+n];
#endif

  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
#ifdef RADIATION_MHD
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */
   rms_error += SQR(total_error.Er) + SQR(total_error.Fr1) + SQR(total_error.Fr2)
                + SQR(total_error.Fr3);

  rms_error = sqrt(rms_error)/(double)count;


/* Print error to file "LinWave-errors.#.dat", where #=wave_flag  */

#ifdef MPI_PARALLEL
  fname = ath_fname("../","LinWave-errors",NULL,NULL,1,NULL,NULL,"dat");
#else
  fname = ath_fname(NULL,"LinWave-errors",NULL,NULL,1,NULL,NULL,"dat");
#endif

/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      free(fname);
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
#ifdef RADIATION_MHD
    fprintf(fp,"  B1c  B2c  B3c");
#endif /* MHD */
    fprintf(fp,"  Er  Fr1  Fr2  Fr3");

    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double) count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count));

#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(total_error.E/(double)count));
#endif /* ISOTHERMAL */

#ifdef RADIATION_MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
#endif /* MHD */
   fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.Er/(double)count),
	  (total_error.Fr1/(double)count),
	  (total_error.Fr2/(double)count),
	(total_error.Fr3/(double)count));	



  fprintf(fp,"\n");

  fclose(fp);
  free(fname);

  return;
}
