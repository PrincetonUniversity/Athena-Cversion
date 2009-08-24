#include "copyright.h"
/*==============================================================================
 * FILE: strat.c
 *
 * PURPOSE:  Problem generator for stratified 3D shearing sheet.  Based on the 
 *   initial conditions described in "Three-dimensional Magnetohydrodynamic
 *   Simulations of Vertically Stratified Accretion Disks" by Stone, Hawley,
 *   Gammie & Balbus.
 *
 * Several different field configurations and perturbations are possible:
 *
 *  ifield = 0 - uses field set by choice of ipert flag
 *  ifield = 1 - Bz=B0sin(kx*x1) field with zero-net-flux [default] (kx input)
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Stone, J., Hawley, J. & Balbus, S. A., ApJ 463, 656-673 (1996)
 *            Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992)
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()           - random number generator from NR
 * ShearingBoxPot() - tidal potential in 3D shearing box
 * expr_*()         - computes new output variables
 * hst_*            - adds new history variables
 *============================================================================*/

static double ran2(long int *idum);
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const Grid *pG, const int i, const int j, const int k);
static Real expr_beta(const Grid *pG, const int i, const int j, const int k);
static Real expr_ME(const Grid *pG, const int i, const int j, const int k);
static Real expr_KE(const Grid *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const Grid *pG, const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const Grid *pG, const int i, const int j, const int k);
#endif
#ifdef MHD
static Real hst_Bx(const Grid *pG, const int i, const int j, const int k);
static Real hst_By(const Grid *pG, const int i, const int j, const int k);
static Real hst_Bz(const Grid *pG, const int i, const int j, const int k);
static Real hst_BxBy(const Grid *pG, const int i, const int j, const int k);
#endif /* MHD */
void output_1d(Grid *pGrid, Domain *pD, Output *pOut);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  FILE *fp;
  Real xFP[160],dFP[160],vxFP[160],vyFP[160];
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,i,j,k,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax,Lx,Ly,Lz;
  Real den = 1.0, pres, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real beta,B0,kx,ky,kz,amp;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  double rval;

  if (pGrid->Nx2 == 1){
    ath_error("[problem]: HGB only works on a 2D or 3D grid\n");
  }

/* Read problem parameters.  Note Omega set to 10^{-3} by default */
  pres=den*Iso_csound2;
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  B0 = sqrt((double)(2.0*pres/beta));
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->is + pGrid->idisp;
  jxs = pGrid->js + pGrid->jdisp;
  kxs = pGrid->ks + pGrid->kdisp;
  iseed = -1 - (ixs + pDomain->Nx1*(jxs + pDomain->Nx2*kxs));

/* Initialize boxsize */
  xmin = par_getd("grid","x1min");
  xmax = par_getd("grid","x1max");
  Lx = xmax - xmin;

  xmin = par_getd("grid","x2min");
  xmax = par_getd("grid","x2max");
  Ly = xmax - xmin;

  xmin = par_getd("grid","x3min");
  xmax = par_getd("grid","x3max");
  Lz = xmax - xmin;

/* initialize wavenumbers, given input number of waves per L */
  nwx = par_geti_def("problem","nwx",1);
  nwy = par_geti_def("problem","nwy",1);
  nwz = par_geti_def("problem","nwz",1);
  kx = (2.0*PI/Lx)*((double)nwx);  /* nxw should be -ve for leading wave */
  ky = (2.0*PI/Ly)*((double)nwy);
  kz = (2.0*PI/Lz)*((double)nwz);



  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 */
      if (ipert == 1) {
        rval = amp*(ran2(&iseed) - 0.5);
	rd = den*exp(-x3*x3);
#ifdef ADIABATIC
        rp = pres*rd*(1.0 + 2.0*rval);
#else
        rd *= (1.0 + 2.0*rval);
#endif
/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = 0.4*rval*sqrt(pres/den);
      }

/* Initialize d, M, and P.  For 3D shearing box M1=Vx, M2=Vy, M3=Vz
 * With FARGO do not initialize the background shear */ 

      pGrid->U[k][j][i].d  = rd;
      pGrid->U[k][j][i].M1 = rd*rvx;
      pGrid->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pGrid->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = rp/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/rd;
#endif

/* Initialize magnetic field.  For 3D shearing box B1=Bx, B2=By, B3=Bz
 *  ifield = 0 - 
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 */
#ifdef MHD
      if (ifield == 0) {
        pGrid->B1i[k][j][i] = rbx;
        pGrid->B2i[k][j][i] = rby;
        pGrid->B3i[k][j][i] = rbz;
        if (i==ie) pGrid->B1i[k][j][ie+1] =  pGrid->B1i[k][j][is];
        if (j==je) pGrid->B2i[k][je+1][i] =  pGrid->B2i[k][js][i];
        if (k==ke) pGrid->B3i[ke+1][j][i] =  pGrid->B3i[ks][j][i];
      }
      if (ifield == 1) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }
      if (ifield == 2) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0;
      }
      if (ifield == 3) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0*(cos((double)kx*x1));
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0*(cos((double)kx*x1));
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }

#ifdef ADIABATIC
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
         + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
#endif /* MHD */
    }
  }}
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

/* enroll gravitational potential function */

  StaticGravPot = ShearingBoxPot;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */

  data_output_enroll(pGrid->time,6.2831853e2,0,output_1d,NULL,NULL,0,0.0,0.0,0,0);

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
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
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

/* enroll gravitational potential function */

  StaticGravPot = ShearingBoxPot;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */

  data_output_enroll(pGrid->time,6.2831853e2,0,output_1d,NULL,NULL,0,0.0,0.0,0,0);

  return;
}

/* Get_user_expression computes dVy */
Gasfun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  else if(strcmp(expr,"beta")==0) return expr_beta;
  else if(strcmp(expr,"ME")==0) return expr_ME;
  else if(strcmp(expr,"KE")==0) return expr_KE;
  else if(strcmp(expr,"BxBy")==0) return hst_BxBy;
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

void gasvshift(const Real x1, const Real x2, const Real x3, Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3, Real *w1, Real *w2, Real *w3)
{
  return;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

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
 * ShearingBoxPot */

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0,z;
  static Real zt,zb;
  static int frst=1;

#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  /* Ensure vertical periodicity in ghost zones */
  if (frst == 1) {
    zb = par_getd("grid","x3min");
    zt = par_getd("grid","x3max");
    frst=0;
  }
  if(x3 > zt) 
    z=x3-zt+zb;
  else if (x3 < zb)
    z=x3-zb+zt;
  else
    z=x3;
/*
  phi += 0.5*Omega*Omega*z*z;
*/

  phi += 0.5*Omega_0*Omega_0*(SQR(2.0-sqrt(SQR(2.0-fabs(z)) + 0.01)));

  return phi;
}

/*------------------------------------------------------------------------------
 * expr_dV2: computes delta(Vy) 
 */

static Real expr_dV2(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

/*------------------------------------------------------------------------------
 * expr_beta: computes beta=P/(B^2/8pi)  
 */

static Real expr_beta(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;
  return (2.0*Iso_csound2*pG->U[k][j][i].d/B2);
#else
  return NULL;
#endif
}

/*------------------------------------------------------------------------------
 * expr_ME: computes B^2/8pi
 */

static Real expr_ME(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;
  return (B2/2.0);
#else
  return NULL;
#endif
}
/*------------------------------------------------------------------------------
 * expr_KE: computes dens*(Vx^2+Vy^2+Vz^2)/2
 */

static Real expr_KE(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,Vy,Vx,Vz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  Vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  Vz = pG->U[k][j][i].M3/pG->U[k][j][i].d;

  return pG->U[k][j][i].d*(Vx*Vx + Vy*Vy + Vz*Vz)/2.0;

}

/*------------------------------------------------------------------------------
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

static Real hst_rho_dVy2(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  return pG->U[k][j][i].d*dVy*dVy;
}

#ifdef ADIABATIC
static Real hst_E_total(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = ShearingBoxPot(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*------------------------------------------------------------------------------
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

#ifdef MHD
static Real hst_Bx(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

static Real hst_By(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

static Real hst_Bz(const Grid *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}

static Real hst_BxBy(const Grid *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c*pG->U[k][j][i].B2c;
}
#endif /* MHD */

void output_1d(Grid *pGrid, Domain *pD, Output *pOut)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int tot1d,i1d,nzmx,kg,kdisp;
  int dnum = pOut->num;
  double darea,**out1d;

  FILE *p_1dfile;
  char *fname;
  int area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int my_area_rat; /* My (Grid Volume)/(dx1*dx2*dx3) */
  int err;
#endif

#ifdef MHD
  tot1d=5;
#else
  tot1d=3;
#endif /* MHD */
 
  nzmx=pD->kde+1;
  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->kdisp;
  /* Compute 1d averaged variables */
  for (k=ks; k<=ke; k++) {
    kg=k+kdisp;
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	i1d=0;
	out1d[kg][i1d] += pGrid->U[k][j][i].d;
	i1d++;
	out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
	i1d++;
	out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#ifdef MHD
	i1d++;
	out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
	i1d++;
	out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif
      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = (pD->ide - pD->ids + 1)*(pD->jde - pD->jds + 1); 

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    err = MPI_Reduce(my_out1d, g_out1d, nzmx,
		     MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if(err)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",err);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }

  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(pGrid->my_id == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pGrid->outfilename,num_digit,dnum,NULL,"1d");
#else
  fname = ath_fname(NULL,pGrid->outfilename,num_digit,dnum,NULL,"1d");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef MHD
    fprintf(p_1dfile,"%G %G %G %G %G\n",out1d[k][0],out1d[k][1],out1d[k][2],
	    out1d[k][3],out1d[k][4]);
#else
    fprintf(p_1dfile,"%G %G %G\n",out1d[k][0],out1d[k][1],out1d[k][2]);
#endif /* MHD */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
}


