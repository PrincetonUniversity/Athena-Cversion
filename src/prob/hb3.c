#include "copyright.h"
/*==============================================================================
 * FILE: hb3.c
 *
 * PURPOSE: Problem generator for 2D MRI simulations using the shearing sheet
 *   based on "A powerful local shear instability in weakly magnetized disks.
 *   III - Long-term evolution in a shearing sheet" by Hawley & Balbus.  This
 *   is the third of the HB papers on the MRI, thus hb3.
 * Several different field configurations and perturbations are possible:
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ipert = 1 - random perturbations to P [default, used by HB]
 *  ipert = 2 - uniform Vx=amp
 *
 * REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef SHEARING_BOX
#error : The HB3 problem requires shearing-box to be enabled.
#endif /* HYDRO */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2() - random number generator from NR
 * ShearingBoxPot() - tidal potential in 2D shearing box
 * expr_dV3() - computes delta(Vy)
 * hst_rho_Vx_dVy () - new history variable
 * hst_E_total() - new history variable
 * hst_dEk() - new history variable
 * hst_Bx()  - new history variable
 * hst_By()  - new history variable
 * hst_Bz()  - new history variable
 * hst_BxBy() - new history variable
 *============================================================================*/

static double ran2(long int *idum);
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Real expr_dV3(const Grid *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j,const int k);
static Real hst_dEk(const Grid *pG, const int i, const int j, const int k);
static Real hst_E_total(const Grid *pG, const int i, const int j, const int k);
#ifdef MHD
static Real hst_Bx(const Grid *pG, const int i, const int j, const int k);
static Real hst_By(const Grid *pG, const int i, const int j, const int k);
static Real hst_Bz(const Grid *pG, const int i, const int j, const int k);
static Real hst_BxBy(const Grid *pG, const int i, const int j, const int k);
#endif

/* boxsize, made a global variable so can be accessed by bval, etc. routines */
static Real Lx;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks;
  int i,j,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1, x2, x3, x1min, x1max;
  Real den = 1.0, pres = 1.0e-5, rd, rp, rvx;
  Real beta,B0,kx,amp;
  double rval;

  if (pGrid->Nx2 == 1){
    ath_error("[problem]: HB3 only works on a 2D grid\n");
  }

  if (pGrid->Nx3 > 1){
    ath_error("[problem]: HB3 does not work on 3D grid\n");
  }

/* Initialize boxsize */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;
  kx = 2.0*PI/Lx;

/* Read problem parameters */
  Omega = par_getd_def("problem","omega",1.0e-3);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  B0 = sqrt((double)(2.0*pres/beta));
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - isentropic perturbations to P & d [default]
 *  ipert = 2 - uniform Vx=amp, sinusoidal density
 *  ipert = 3 - random perturbations to P [used by HB]
 */
      if (ipert == 1) {
        rval = 1.0 + amp*(ran2(&iseed) - 0.5);
#ifdef ADIABATIC
        rp = rval*pres;
        rd = pow(rval,1.0/Gamma)*den;
#else
        rd = rval*den;
#endif
        rvx = 0.0;
      }
      if (ipert == 2) {
        rp = pres;
        rd = den*(1.0 + 0.1*sin((double)kx*x1));
        rvx = amp*sqrt(Gamma*pres/den);
      }
      if (ipert == 3) {
        rval = 1.0 + amp*(ran2(&iseed) - 0.5);
#ifdef ADIABATIC
        rp = rval*pres;
        rd = den;
#else
        rd = rval*den;
#endif
        rvx = 0.0;
      }

/* Initialize d, M, and P.  For 2D shearing box M1=Vx, M2=Vz, M3=Vy */ 

      pGrid->U[ks][j][i].d  = rd;
      pGrid->U[ks][j][i].M1 = rd*rvx;
      pGrid->U[ks][j][i].M2 = 0.0;
      pGrid->U[ks][j][i].M3 = -rd*1.5*Omega*x1;
#ifdef ADIABATIC
      pGrid->U[ks][j][i].E = rp/Gamma_1
        + 0.5*(SQR(pGrid->U[ks][j][i].M1) + SQR(pGrid->U[ks][j][i].M3))/rd;
#endif

/* Initialize magnetic field.  For 2D shearing box B1=Bx, B2=Bz, B3=By
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 */
#ifdef MHD
      if (ifield == 1) {
        pGrid->U[ks][j][i].B1c = 0.0;
        pGrid->U[ks][j][i].B2c = B0*(sin((double)kx*x1));
        pGrid->U[ks][j][i].B3c = 0.0;
        pGrid->B1i[ks][j][i] = 0.0;
        pGrid->B2i[ks][j][i] = B0*(sin((double)kx*x1));
        pGrid->B3i[ks][j][i] = 0.0;
        if (i==ie) pGrid->B1i[ks][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[ks][je+1][i] = B0*(sin((double)kx*x1));
      }
      if (ifield == 2) {
        pGrid->U[ks][j][i].B1c = 0.0;
        pGrid->U[ks][j][i].B2c = B0;
        pGrid->U[ks][j][i].B3c = 0.0;
        pGrid->B1i[ks][j][i] = 0.0;
        pGrid->B2i[ks][j][i] = B0;
        pGrid->B3i[ks][j][i] = 0.0;
        if (i==ie) pGrid->B1i[ks][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[ks][je+1][i] = B0;
      }
#ifdef ADIABATIC
      pGrid->U[ks][j][i].E += 0.5*(SQR(pGrid->U[ks][j][i].B1c)
         + SQR(pGrid->U[ks][j][i].B2c) + SQR(pGrid->U[ks][j][i].B3c));
#endif
#endif /* MHD */
    }
  }

/* enroll gravitational potential function, shearing sheet BC functions */

  StaticGravPot = ShearingBoxPot;

/* enroll new history variables */

  dump_history_enroll(hst_dEk, "<0.5rho(Vx^2+4dVy^2)>");
  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */


  return;
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

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll special boundary value functions,
 *    and initialize gravity on restarts
 */

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  Real x1min, x1max;

  Omega = par_getd_def("problem","omega",1.0e-3);

/* Must recompute global variable Lx needed by BC routines */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;

  StaticGravPot = ShearingBoxPot;

  return;
}

/* Get_user_expression computes dVy */
Gasfun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV3;
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

GVDFun_t get_usr_gasvshift(const char *name)
{
  return NULL;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
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

/*------------------------------------------------------------------------------
 * ShearingBoxPot: 
 */

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3){
  return -1.5*Omega*Omega*x1*x1;  
}

/*------------------------------------------------------------------------------
 * expr_dV3: computes delta(Vy) 
 */

static Real expr_dV3(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return (pG->U[k][j][i].M3/pG->U[k][j][i].d + 1.5*Omega*x1);
}

/*------------------------------------------------------------------------------
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 */

static Real hst_rho_Vx_dVy(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d + 1.5*Omega*x1);
}

/*------------------------------------------------------------------------------
 * hst_dEk: computes 0.5*(Vx^2 + 4(\delta Vy)^2), which for epicyclic motion
 *   is a constant, added as history variable 
 */

static Real hst_dEk(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  Real dMy, dE;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  dMy = (pG->U[k][j][i].M2 + 1.5*Omega*x1*pG->U[k][j][i].d);
  dE = 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 + 4.0*dMy*dMy)/pG->U[k][j][i].d;

  return dE;
}

/*------------------------------------------------------------------------------
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_E_total(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = ShearingBoxPot(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}

/*------------------------------------------------------------------------------
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

#ifdef MHD
static Real hst_Bx(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][pG->is].B1c;
}

static Real hst_By(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][pG->js][i].B2c;
}

static Real hst_Bz(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[pG->ks][j][i].B3c;
}

static Real hst_BxBy(const Grid *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c*pG->U[k][j][i].B2c;
}

#endif /* MHD */
