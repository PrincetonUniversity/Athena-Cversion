#include "copyright.h"
/*==============================================================================
 * FILE: cpaw2d.c
 *
 * PURPOSE: Problem generator for 2-D circularly polarized Alfven wave (CPAW)
 *   test.  Works for any arbitrary wavevector in the x1-x2 plane.  The wave is
 *   defined with reference to a coordinate system (x,y,z) with transformation
 *   rules to the code coordinate system (x1,x2,x3)
 *      x =  x1*cos(alpha) + x2*sin(alpha)
 *      y = -x1*sin(alpha) + x2*cos(alpha)
 *      z = x3
 *   The magnetic field is given by:
 *     B_x = b_par
 *     B_y = b_perp*sin(k*x)
 *     B_z = b_perp*cos(k*x)   where k = 2.0*PI/lambda
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
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The cpaw2d test only works for MHD.
#endif


/* Initial solution, shared with Userwork_after_loop to compute L1 error.
 * Vector potential A3 defined in prvate function below.  B_par, etc. must
 * be defined as globals to be used by A3() */
 
static Gas **Soln=NULL;
static Real A3(const Real x1, const Real x2);
Real sin_a, cos_a, b_par, b_perp;
Real k_par;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int nx1, nx2;
  Real angle;    /* Angle the wave direction makes with the x1-direction */
  Real dx1 = pGrid->dx1;
  Real dx2 = pGrid->dx2;
  Real x1,x2,x3,cs,sn;
  Real v_par, v_perp, den, pres;
  Real lambda; /* Wavelength */
 
  nx1 = (ie - is + 1) + 2*nghost;
  nx2 = (je - js + 1) + 2*nghost;

  if (pGrid->Nx1 == 1) {
    ath_error("[cpaw2d] Grid must be 2D with Nx1>1");
  }

  if ((Soln = (Gas**)calloc_2d_array(nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[pgflow]: Error allocating memory for Soln\n");

/* An angle =  0.0 is a wave aligned with the x1-direction. */
/* An angle = 90.0 is a wave aligned with the x2-direction. */

  angle = par_getd("problem","angle");

/* Compute the sin and cos of the angle and the wavelength. */

  if (angle == 0.0) {
    sin_a = 0.0;
    cos_a = 1.0;
    lambda = pGrid->Nx1*pGrid->dx1; /* Put one wavelength in the grid */
  }
  else if (angle == 90.0) {
    sin_a = 1.0;
    cos_a = 0.0;
    lambda = pGrid->Nx2*pGrid->dx2; /* Put one wavelength in the grid */
  }
  else {

/* We put 1 wavelength in each direction.  Hence the wavelength
 *     lambda = pGrid->Nx1*pGrid->dx1*cos_a;
 *     AND  lambda = pGrid->Nx2*pGrid->dx2*sin_a;
 *     are both satisfied. */

    if((pGrid->Nx1*pGrid->dx1) == (pGrid->Nx2*pGrid->dx2)){
      cos_a = sin_a = sqrt(0.5);
    }
    else{
      angle = atan((double)(pGrid->Nx1*pGrid->dx1)/(pGrid->Nx2*pGrid->dx2));
      sin_a = sin(angle);
      cos_a = cos(angle);
    }
/* Use the larger angle to determine the wavelength */
    if (cos_a >= sin_a) {
      lambda = pGrid->Nx1*pGrid->dx1*cos_a;
    } else {
      lambda = pGrid->Nx2*pGrid->dx2*sin_a;
    }
  }

/* Initialize k_parallel */

  k_par = 2.0*PI/lambda;
  b_par = par_getd("problem","b_par");
  den = 1.0;

  printf("va_parallel = %g\nlambda = %g\n",b_par/sqrt(den),lambda);

  b_perp = par_getd("problem","b_perp");
  v_perp = b_perp/sqrt((double)den);

/* The gas pressure and parallel velocity are free parameters. */

  pres = par_getd("problem","pres");
  v_par = par_getd("problem","v_par");

/* Use the vector potential to initialize the interface magnetic fields
 * The iterface fields are located at the left grid cell face normal */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        cs = cos(k_par*(x1*cos_a + x2*sin_a));

        x1 -= 0.5*pGrid->dx1;
        x2 -= 0.5*pGrid->dx2;

        pGrid->B1i[k][j][i] = -(A3(x1,x2+dx2) - A3(x1,x2))/dx2;
        pGrid->B2i[k][j][i] =  (A3(x1+dx1,x2) - A3(x1,x2))/dx1;
        pGrid->B3i[k][j][i] = b_perp*cs;
      }
    }
  }
  if (pGrid->Nx3 > 1) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        cs = cos(k_par*(x1*cos_a + x2*sin_a));
        pGrid->B3i[ke+1][j][i] = b_perp*cs;
      }
    }
  }

/* Now initialize the cell centered quantities */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

        sn = sin(k_par*(x1*cos_a + x2*sin_a));
        cs = cos(k_par*(x1*cos_a + x2*sin_a));

        Soln[j][i].d  = den;
        Soln[j][i].M1 = den*(v_par*cos_a - v_perp*sn*sin_a);
        Soln[j][i].M2 = den*(v_par*sin_a + v_perp*sn*cos_a);
        Soln[j][i].M3 = den*v_perp*cs;
        pGrid->U[k][j][i].d  = Soln[j][i].d;
        pGrid->U[k][j][i].M1 = Soln[j][i].M1;
        pGrid->U[k][j][i].M2 = Soln[j][i].M2;
        pGrid->U[k][j][i].M3 = Soln[j][i].M3;

        Soln[j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + pGrid->B1i[k][j][i+1]);
        Soln[j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + pGrid->B2i[k][j+1][i]);
        Soln[j][i].B3c = b_perp*cs;
        pGrid->U[k][j][i].B1c = Soln[j][i].B1c;
        pGrid->U[k][j][i].B2c = Soln[j][i].B2c;
        pGrid->U[k][j][i].B3c = Soln[j][i].B3c;

#ifndef ISOTHERMAL
        Soln[j][i].E = pres/Gamma_1 + 0.5*(SQR(pGrid->U[k][j][i].B1c) +
                 SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c) )
          + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) +
                 SQR(pGrid->U[k][j][i].M3) )/den;
        pGrid->U[k][j][i].E = Soln[j][i].E;
#endif
      }
    }
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
 * A3() - computes vector potential to initialize fields
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
  int i,j,is,ie,js,je,ks,Nx1,Nx2;
  Real rms_error=0.0;
  Gas error;
  FILE *fp;
  char *fname;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  Nx1 = (ie-is+1);
  Nx2 = (je-js+1);
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

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[ks][j][i].d   - Soln[j][i].d  );
      error.M1  += fabs(pGrid->U[ks][j][i].M1  - Soln[j][i].M1 );
      error.M2  += fabs(pGrid->U[ks][j][i].M2  - Soln[j][i].M2 );
      error.M3  += fabs(pGrid->U[ks][j][i].M3  - Soln[j][i].M3 );
      error.B1c += fabs(pGrid->U[ks][j][i].B1c - Soln[j][i].B1c);
      error.B2c += fabs(pGrid->U[ks][j][i].B2c - Soln[j][i].B2c);
      error.B3c += fabs(pGrid->U[ks][j][i].B3c - Soln[j][i].B3c);
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[ks][j][i].E   - Soln[j][i].E  );
#endif /* ISOTHERMAL */
    }
  }

/* Compute RMS error over all variables */

  rms_error = SQR(error.d) + SQR(error.M1) + SQR(error.M2) + SQR(error.M3);
  rms_error += SQR(error.B1c) + SQR(error.B2c) + SQR(error.B3c);
#ifndef ISOTHERMAL
  rms_error += SQR(error.E);
#endif /* ISOTHERMAL */
  rms_error = sqrt(rms_error)/(double)(Nx1*Nx2);

/* Print error to file "cpaw2d-errors.dat" */

  fname = "cpaw2d-errors.dat";
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
    fprintf(fp,"# Nx1   Nx2  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
    fprintf(fp,"  B1c  B2c  B3c");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d %e  %e  %e  %e  %e",Nx1,Nx2,rms_error,
          (error.d/(double)(Nx1*Nx2)),
          (error.M1/(double)(Nx1*Nx2)),
          (error.M2/(double)(Nx1*Nx2)),
          (error.M3/(double)(Nx1*Nx2)));
#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(error.E/(double)(Nx1*Nx2)));
#endif /* ISOTHERMAL */
  fprintf(fp,"  %e  %e  %e",
          (error.B1c/(double)(Nx1*Nx2)),
          (error.B2c/(double)(Nx1*Nx2)),
          (error.B3c/(double)(Nx1*Nx2)));
  fprintf(fp,"\n");

  fclose(fp);

  return;
}

/*---------------------------------------------------------------------------
 * A3: Define a scalar potential A3 such that:
 *     B_x = - $\partial A3 / \partial y$
 *     B_y =   $\partial A3 / \partial x$
 *   Then A3 is given in the function below.  */

static Real A3(const Real x1, const Real x2)
{
  return b_par*(x1*sin_a - x2*cos_a) 
     - (b_perp/k_par)*cos(k_par*(x1*cos_a + x2*sin_a));
}
