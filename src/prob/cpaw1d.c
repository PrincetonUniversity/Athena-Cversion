#include "copyright.h"
/*==============================================================================
 * FILE: cpaw1d.c
 *
 * PURPOSE: Problem generator for 1-D circularly polarized Alfven wave (CPAW)
 *   test.  Only works in 1D (wavevector in x).  Tests in 2D and 3D are 
 *   initialized with different functions.
 *
 *   Can be used for either standing (problem/v_par=1.0) or travelling
 *   (problem/v_par=0.0) waves.
 *
 * USERWORK_AFTER_LOOP function computes L1 error norm in solution by comparing
 *   to initial conditions.  Problem must be evolved for an integer number of
 *   wave periods for this to work.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   problem - problem generator
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifndef MHD
#error : The cpaw1d test only works for mhd.
#endif

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
static Gas *Soln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pgrid)
{
  int i, is = pgrid->is, ie = pgrid->ie;
  int j, js = pgrid->js;
  int k, ks = pgrid->ks;
  Real x1,x2,x3,cs,sn,b_par,b_perp,lambda,k_par,v_par,v_perp,den,pres;
  Soln = (Gas*)malloc(((ie-is+1)+2*nghost)*sizeof(Gas));
  if (Soln == NULL) ath_error("[cpaw1d] Error initializing solution array");

  if (pGrid->Nx2 > 1 || pGrid->Nx3 > 1) {
    ath_error("[cpaw1d] grid must be 1D");
  }

/* Put one wavelength on the grid, and initialize k_parallel */

  lambda = pgrid->Nx1*pgrid->dx1; 
  k_par = 2.0*PI/lambda;

  b_par = par_getd("problem","b_par");
  den = 1.0;
  b_perp = par_getd("problem","b_perp");
  v_perp = b_perp/sqrt((double)den);

/* The gas pressure and parallel velocity are free parameters. */
  pres = par_getd("problem","pres");
  v_par = par_getd("problem","v_par");

/* Setup circularily polarized AW solution  */

  for (i=is; i<=ie+1; i++) {
    cc_pos(pgrid,i,js,ks,&x1,&x2,&x3);

    sn = sin(k_par*x1);
    cs = cos(k_par*x1);

    Soln[i].d = den;

    Soln[i].M1 = den*v_par;
    Soln[i].M2 = den*v_perp*sn;
    Soln[i].M3 = den*v_perp*cs;
 
    Soln[i].B1c = b_par;
    Soln[i].B2c = b_perp*sn;
    Soln[i].B3c = b_perp*cs;
#ifndef ISOTHERMAL
    Soln[i].E = pres/Gamma_1 
      + 0.5*den*(v_par*v_par + v_perp*v_perp)
      + 0.5*(b_par*b_par + b_perp*b_perp);
#endif
  }

/* set code variables to solution */

  for (i=is; i<=ie+1; i++) {
    pgrid->U[ks][js][i].d  = Soln[i].d;
    pgrid->U[ks][js][i].M1 = Soln[i].M1;
    pgrid->U[ks][js][i].M2 = Soln[i].M2;
    pgrid->U[ks][js][i].M3 = Soln[i].M3;

    pgrid->U[ks][js][i].B1c = pgrid->B1i[ks][js][i] = Soln[i].B1c;
    pgrid->U[ks][js][i].B2c = pgrid->B2i[ks][js][i] = Soln[i].B2c;
    pgrid->U[ks][js][i].B3c = pgrid->B3i[ks][js][i] = Soln[i].B3c;
#ifndef ISOTHERMAL
    pgrid->U[ks][js][i].E = Soln[i].E;
#endif
  }

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(Grid *pG, FILE *fp){
  return;
}

void problem_read_restart(Grid *pG, FILE *fp){
  return;
}

Gasfun_t get_usr_expr(const char *expr){
  return NULL;
}

void Userwork_in_loop(Grid *pGrid)
{
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in CPAW,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(Grid *pGrid)
{
  int i=0,is,ie,js,ks,Nx1;
  Real rms_error=0.0;
  Gas error;
  FILE *fp;
  char *fname;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js;
  ks = pGrid->ks;
  Nx1 = (ie-is+1);
  error.d = 0.0;
  error.M1 = 0.0;
  error.M2 = 0.0;
  error.M3 = 0.0;
  error.B1c = 0.0;
  error.B2c = 0.0;
  error.B3c = 0.0;
#ifndef ISOTHERMAL
  error.E = 0.0;
#endif /* ISOTHERMAL */

/* compute L1 error in each variable, and rms total error */

  for (i=is; i<=ie; i++) {
    error.d   += fabs(pGrid->U[ks][js][i].d   - Soln[i].d  );
    error.M1  += fabs(pGrid->U[ks][js][i].M1  - Soln[i].M1 );
    error.M2  += fabs(pGrid->U[ks][js][i].M2  - Soln[i].M2 );
    error.M3  += fabs(pGrid->U[ks][js][i].M3  - Soln[i].M3 );
    error.B1c += fabs(pGrid->U[ks][js][i].B1c - Soln[i].B1c);
    error.B2c += fabs(pGrid->U[ks][js][i].B2c - Soln[i].B2c);
    error.B3c += fabs(pGrid->U[ks][js][i].B3c - Soln[i].B3c);
#ifndef ISOTHERMAL
    error.E   += fabs(pGrid->U[ks][js][i].E   - Soln[i].E  );
#endif /* ISOTHERMAL */
  }

/* Compute RMS error over all variables */

  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
#ifndef ISOTHERMAL
  rms_error += SQR(error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)Nx1;

/* Print error to file "cpaw1d-errors.dat" */

  fname = "cpaw1d-errors.dat";
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      return;
    }
/* write out some header information */
    fprintf(fp,"# Nx1 RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %e  %e  %e  %e  %e",Nx1,rms_error,
          (error.d/(double)Nx1),
          (error.M1/(double)Nx1),
          (error.M2/(double)Nx1),
          (error.M3/(double)Nx1));
#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(error.E/(double)Nx1));
#endif /* ISOTHERMAL */
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)Nx1),
          (error.B2c/(double)Nx1),
          (error.B3c/(double)Nx1));
  fprintf(fp,"\n");

  fclose(fp);

  return;
}
