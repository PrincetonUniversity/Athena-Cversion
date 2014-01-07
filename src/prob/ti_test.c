#include "copyright.h"
/*==============================================================================
 * FILE: ti_test.c
 *
 * PURPOSE: Problem generator for thermal instability with operator split cooling.
 *          iprob=1 for eigenmode test
 *          iprob=2 for random perturbation
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef OPERATOR_SPLIT_COOLING
#error This problem only works with --enable-cooling option
#endif

#ifdef MHD
#error This problem doesn't work with MHD
#endif


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()
 *============================================================================*/

static double ran2(long int *idum);
static Real logd(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  int iprob;
  Real n0,T0,kappa,P0;
  long int iseed = -1-myID_Comm_world;

  ConstS consts;

  Real x1,x2,x3;
  Real x1min,x1max,Lx,kwx,amp;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Get box size and set wavenumber */
  x1min = pDomain->MinX[0];
  x1max = pDomain->MaxX[0];
  Lx = x1max - x1min;
  kwx = 2.0*PI/Lx;
  amp   = par_getd("problem","amp");

/* Set Units */
  init_consts(&consts);
  pGrid->units.Dcode = 1.4*consts.mH;
  pGrid->units.Vcode = consts.kms;
  pGrid->units.Lcode = consts.pc;
  pGrid->units.Tcode = pGrid->units.Lcode/pGrid->units.Vcode;
//  pGrid->units.Pcode = pGrid->units.Dcode*SQR(pGrid->units.Vcode);
  pGrid->heat0 = 2.e-26;
  pGrid->heat_ratio = 1.0;

/* Read problem parameters */

  T0    = par_getd("problem","T0");
  n0    = par_getd("problem","n0");
  P0    = 1.1*n0*T0*consts.kB/pGrid->units.Dcode/SQR(pGrid->units.Vcode);
  iprob = par_geti_def("problem","iprob",2);
  if(myID_Comm_world == 0) printf("P0/k0 = %g, P0 in code unit = %g\n",1.1*n0*T0,P0);

#ifdef THERMAL_CONDUCTION
  kappa = par_getd("problem","kappa");
  kappa_iso = kappa/(consts.kB*pGrid->units.Lcode*pGrid->units.Vcode);
  if(myID_Comm_world == 0) printf("kappa in c.g.s. = %g, in code unit = %g\n",kappa,kappa_iso);
#endif

/* Constant density and temperature initially */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

      if(iprob == 1) pGrid->U[k][j][i].d = n0*(1+amp*sin(kwx*x1));
      if(iprob == 2) pGrid->U[k][j][i].d = n0*(1+amp*(ran2(&iseed)-0.5));
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
      pGrid->U[k][j][i].M3 = 0.0;
#ifndef BAROTROPIC
      pGrid->U[k][j][i].E = P0/Gamma_1
             + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
             + SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif
    }
  }}

}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
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

static Real logd(const GridS *pG, const int i, const int j, const int k)
{
  return log10(pG->U[k][j][i].d);
}


ConsFun_t get_usr_expr(const char *expr)
{ 
  if(strcmp(expr,"logd")==0) return logd;
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
double ran2(long int *idum){
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
