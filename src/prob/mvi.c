#include "copyright.h"
/*==============================================================================
 * FILE: mvi.c
 *
 * PURPOSE:  Problem generator for 3D shearing sheet studies of MVI.  Based on
 *   hgb.c.  New diagnostics added which are relevant only with anisotropic
 *   viscosity, thus the motivation for a new problem generator.
 *
 * Several different field configurations and perturbations are possible:
 *
 *  ifield = 0 - uses field set by choice of ipert flag
 *  ifield = 1 - Bz=B0sin(kx*x1) field with zero-net-flux [default] (kx input)
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
 *  ifield = 5 - uniform By
 *
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *
 * Code must be configured using --enable-shearing-box --enable-viscosity
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : MVI only works for MHD.
#endif /* MHD */
#ifndef VISCOSITY
#error : MVI only works with VISCOSITY defined.
#endif /* VISCOSITY */

Real Lx,Ly,Lz; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);

static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
static Real hst_Visc_flx(const GridS *pG,const int i, const int j, const int k);
static Real hst_BBdV(const GridS *pG,const int i, const int j, const int k);
static Real hst_BrBpdOmega(const GridS *pG,const int i,const int j,const int k);
static Real hst_ThirddivV(const GridS *pG, const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  FILE *fp;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,i,j,k,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real beta=1.0,B0,kx,amp;
  int nwx;  /* input number of waves per Lx [default=1] */
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */

  if (pGrid->Nx[1] == 1){
    ath_error("[problem]: MVI only works on a 2D or 3D grid\n");
  }

/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
  Omega_0 = par_getd_def("problem","Omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);

/* Compute field strength based on beta.  */
#ifdef ISOTHERMAL
  pres = Iso_csound2;
#else
  pres = par_getd("problem","pres");
#endif
  B0 = sqrt((double)(2.0*pres/beta));

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* initialize wavenumbers, given input number of waves per L */
  nwx = par_geti_def("problem","nwx",1);
  kx = (2.0*PI/Lx)*((double)nwx);  /* nxw should be -ve for leading wave */

/* Rescale amp to sound speed for ipert 2,3 */
#ifdef ADIABATIC
  if (ipert == 2) amp *= sqrt(Gamma*pres/den);
#else
  if (ipert == 2) amp *= Iso_csound;
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 */
      if (ipert == 1) {
        rval = amp*(ran2(&iseed) - 0.5);
#ifdef ADIABATIC
        rp = pres*(1.0 + 2.0*rval);
        rd = den;
#else
        rd = den*(1.0 + 2.0*rval);
#endif
/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = 0.4*rval*sqrt(pres/den);
      }
      if (ipert == 2) {
        rp = pres;
        rd = den;
        rvx = amp;
        rvy = 0.0;
        rvz = 0.0;
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
 *  ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
 *  ifield = 5 - uniform By
 */
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
      if (ifield == 4) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0/sqrt(2);
        pGrid->B3i[k][j][i] = B0/sqrt(2);
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0/sqrt(2);
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0/sqrt(2);
      }
      if (ifield == 5) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0;
        pGrid->B3i[k][j][i] = 0.0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;
      }
    }
  }}
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
#ifdef ADIABATIC
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
         + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
      }
    }
  }

/* enroll gravitational potential function */

  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables, only once with SMR  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
    dump_history_enroll(hst_Visc_flx, "<ViscFlux>");
    dump_history_enroll(hst_BBdV, "<BBdV>");
    dump_history_enroll(hst_BrBpdOmega, "<BrBpdO>");
    dump_history_enroll(hst_ThirddivV, "<divV/3>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BxBy, "<-Bx By>");
    frst = 0;
  }

/* With viscosity and/or resistivity, read eta_Ohm and nu */
#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
  Q_AD    = par_getd_def("problem","Q_A",0.0);
#endif
  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);

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

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
/* Read Omega, and with viscosity and/or resistivity, read eta_Ohm and nu */

  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);

#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",0.0);
  Q_Hall  = par_getd_def("problem","Q_H",0.0);
  Q_AD    = par_getd_def("problem","Q_A",0.0);
#endif

  nu_iso = par_getd_def("problem","nu_iso",0.0);
  nu_aniso = par_getd_def("problem","nu_aniso",0.0);

/* enroll gravitational potential function */

  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
  dump_history_enroll(hst_Visc_flx, "<ViscFlux>");
  dump_history_enroll(hst_BBdV, "<BBdV>");
  dump_history_enroll(hst_BrBpdOmega, "<BrBpdO>");
  dump_history_enroll(hst_ThirddivV, "<divV/3>");
#ifdef ADIABATIC
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");

  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{
  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
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
 * UnstratifiedDisk:
 */

static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

/*------------------------------------------------------------------------------
 * expr_dV2: computes delta(Vy) 
 */

static Real expr_dV2(const GridS *pG, const int i, const int j, const int k)
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
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_Visc_flx: r-Phi viscous flux
 * hst_BBdV: BB \cdot dV (anisotropy parameter)
 * hst_BrBpdOmega: estimate of anistropy parameter
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
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

static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k)
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

static Real hst_BrBpdOmega(const GridS *pG,const int i,const int j,const int k)
{
  Real Bx,By,Bz,B02,BBdV;

  Bx = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
  By = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
  Bz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
  B02 = Bx*Bx + By*By + Bz*Bz;
  B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

  BBdV = Bx*By*qshear*Omega_0/B02;
  return BBdV;
}

static Real hst_BBdV(const GridS *pG, const int i, const int j, const int k)
{
  Real Bx,By,Bz,B02,BBdV,Vyp,Vym,x1,x2,x3;

  Bx = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
  By = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
  Bz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
  B02 = Bx*Bx + By*By + Bz*Bz;
  B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

  BBdV  = Bx*Bx*((pG->U[k][j][i+1].M1/pG->U[k][j][i+1].d) -
                 (pG->U[k][j][i-1].M1/pG->U[k][j][i-1].d))/(2.0*pG->dx1);
  BBdV += Bx*By*((pG->U[k][j+1][i].M1/pG->U[k][j+1][i].d) -
                 (pG->U[k][j-1][i].M1/pG->U[k][j-1][i].d))/(2.0*pG->dx2);
  BBdV += Bx*Bz*((pG->U[k+1][j][i].M1/pG->U[k+1][j][i].d) -
                 (pG->U[k-1][j][i].M1/pG->U[k-1][j][i].d))/(2.0*pG->dx3);

  Vyp = pG->U[k][j][i+1].M2/pG->U[k][j][i+1].d;
  Vym = pG->U[k][j][i-1].M2/pG->U[k][j][i-1].d;
#ifdef FARGO
  cc_pos(pG,(i+1),j,k,&x1,&x2,&x3);
  Vyp -= qshear*Omega_0*x1;
  cc_pos(pG,(i-1),j,k,&x1,&x2,&x3);
  Vym -= qshear*Omega_0*x1;
#endif
  BBdV += By*Bx*(Vyp - Vym)/(2.0*pG->dx1);

  Vyp = pG->U[k][j+1][i].M2/pG->U[k][j+1][i].d;
  Vym = pG->U[k][j-1][i].M2/pG->U[k][j-1][i].d;
#ifdef FARGO
  cc_pos(pG,i,(j+1),k,&x1,&x2,&x3);
  Vyp -= qshear*Omega_0*x1;
  cc_pos(pG,i,(j-1),k,&x1,&x2,&x3);
  Vym -= qshear*Omega_0*x1;
#endif
  BBdV += By*By*(Vyp - Vym)/(2.0*pG->dx2);

  Vyp = pG->U[k+1][j][i].M2/pG->U[k+1][j][i].d;
  Vym = pG->U[k-1][j][i].M2/pG->U[k-1][j][i].d;
#ifdef FARGO
  cc_pos(pG,i,j,(k+1),&x1,&x2,&x3);
  Vyp -= qshear*Omega_0*x1;
  cc_pos(pG,i,j,(k-1),&x1,&x2,&x3);
  Vym -= qshear*Omega_0*x1;
#endif
  BBdV += By*Bz*(Vyp - Vym)/(2.0*pG->dx3);

  BBdV += Bz*Bx*((pG->U[k][j][i+1].M3/pG->U[k][j][i+1].d) -
                 (pG->U[k][j][i-1].M3/pG->U[k][j][i-1].d))/(2.0*pG->dx1);
  BBdV += Bz*By*((pG->U[k][j+1][i].M3/pG->U[k][j+1][i].d) -
                 (pG->U[k][j-1][i].M3/pG->U[k][j-1][i].d))/(2.0*pG->dx2);
  BBdV += Bz*Bz*((pG->U[k+1][j][i].M3/pG->U[k+1][j][i].d) -
                 (pG->U[k-1][j][i].M3/pG->U[k-1][j][i].d))/(2.0*pG->dx3);

  BBdV /= B02;
  return BBdV;
}

static Real hst_Visc_flx(const GridS *pG, const int i, const int j, const int k)
{
  Real Bx,By,Bz,B02,BBdV,nud,qa;
  Real divV = 0.0;

  Bx = 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]);
  By = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
  Bz = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
  B02 = Bx*Bx + By*By + Bz*Bz;
  B02 = MAX(B02,TINY_NUMBER); /* limit in case B=0 */

  BBdV = hst_BBdV(pG, i, j, k);

  nud = nu_aniso*pG->U[k][j][i].d;
  qa = nud*(BBdV - ONE_3RD*divV);
  return qa*(3.0*By*Bx/B02);
}

static Real hst_ThirddivV(const GridS *pG, const int i, const int j, const int k)
{
  Real divV;

  divV  = ((pG->U[k][j][i+1].M1/pG->U[k][j][i+1].d)-
           (pG->U[k][j][i-1].M1/pG->U[k][j][i-1].d))/(2.0*pG->dx1);
  divV += ((pG->U[k][j+1][i].M2/pG->U[k][j+1][i].d)-
           (pG->U[k][j-1][i].M2/pG->U[k][j-1][i].d))/(2.0*pG->dx2);
  divV += ((pG->U[k+1][j][i].M3/pG->U[k+1][j][i].d)-
           (pG->U[k-1][j][i].M3/pG->U[k-1][j][i].d))/(2.0*pG->dx3);

  return ONE_3RD*divV;
}

#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*------------------------------------------------------------------------------
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

static Real hst_By(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}

static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c*pG->U[k][j][i].B2c;
}
