#include "copyright.h"
/*==============================================================================
 * FILE: hgb.c
 *
 * PURPOSE: Problem generator for 3D MRI simulations using the shearing sheet
 *   based on "Local Three-dimensional Magnetohydrodynamic Simulations of
 *   Accretion Disks" by Hawley, Gammie & Balbus, or HGB.
 * Several different field configurations and perturbations are possible:
 *  ifield = 1 - Bz=B0sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp
 *
 * To run simulations of stratified disks (include vertical gravity),
 * un-comment the macro STRATIFIED_DISK below.  NOT YET IMPLEMENTED.
 *
 * This file also contains shear_ix1_ox1(), a public function called by
 * set_bvals() which implements the 3D shearing sheet boundary conditions.
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

/* #define STRATIFIED_DISK */

/* prototype for shearing sheet BC function (called by set_bvals) */
void shear_ix1_ox1(Grid *pG, int var_flag);

/* Remapped conserved quantities in ghost zones, and their fluxes */
static Gas **RemapGas=NULL, *Flx=NULL;
Real **pU=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * CompRemapFlux() - 2nd order reconstruction for remap in ghost zones
 * ran2()          - random number generator from NR
 * ShearingBoxPot() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * no_op_VGfun() - no operation void grid function, replaces ix1/ox1 bval fns
 *============================================================================*/

void CompRemapFlux(const Gas U[], const Real eps,
                   const int il, const int iu, Gas Flux[]);
static double ran2(long int *idum);
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const Grid *pG, const int i, const int j, const int k);
static void no_op_VGfun(Grid *pGrid, int swap_phi);

/* boxsize, made a global variable so can be accessed by bval, etc. routines */
static Real Lx,Ly;

/*=========================== PUBLIC FUNCTIONS =================================
 * Contains the usual, plus:
 * shear_ix1_ox1() - shearing sheet BCs in 3D, called by set_bval().  Ensures
 *   ix1/ox1 bval routines are executaed simultanesouly in MPI parallel jobs
 * EvolveEy() - sets Ey at [is] and [ie+1] to keep <Bz>=const.  Called by
 *   integrator.
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ipert,ifield,Nx2;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1, x2, x3, x1min, x1max, x2min, x2max;
  Real den = 1.0, pres = 1.0e-6, rd, rp, rvx, rvy, rvz;
  Real beta,B0,kx,amp;
  double rval;

  if (pGrid->Nx2 == 1){
    ath_error("[problem]: HGB only works on a 2D or 3D grid\n");
  }

/* Initialize boxsize */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;
  kx = 2.0*PI/Lx;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Ly = x2max - x2min;

/* Read problem parameters.  Note Omega set to 10^{-3} by default */
  Omega = par_getd_def("problem","omega",1.0e-3);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  B0 = sqrt((double)(2.0*pres/beta));
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp
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
        rd = den*(1.0 + 0.1*sin((double)kx*x1));
        rvx = amp*sqrt(Gamma*pres/den);
        rvy = 0.0;
        rvz = 0.0;
      }

/* Initialize d, M, and P.  For 3D shearing box M1=Vx, M2=Vy, M3=Vz */ 

      pGrid->U[k][j][i].d  = rd;
      pGrid->U[k][j][i].M1 = rd*rvx;
      pGrid->U[k][j][i].M2 = rd*(rvy - 1.5*Omega*x1);
      pGrid->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = rp/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/rd;
#endif

/* Initialize magnetic field.  For 3D shearing box B1=Bx, B2=By, B3=Bz
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 */
#ifdef MHD
      if (ifield == 1) {
        pGrid->U[k][j][i].B1c = 0.0;
        pGrid->U[k][j][i].B2c = 0.0;
        pGrid->U[k][j][i].B3c = B0*(sin((double)kx*x1));
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }
      if (ifield == 2) {
        pGrid->U[k][j][i].B1c = 0.0;
        pGrid->U[k][j][i].B2c = 0.0;
        pGrid->U[k][j][i].B3c = B0;
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0;
      }
#ifdef ADIABATIC
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
         + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
#endif /* MHD */
    }
  }}

/* enroll gravitational potential function, set x1-bval functions to NoOp
 * function, since public function shear_ix1_ox1() below is called by
 * set_bvals() instead */

  StaticGravPot = ShearingBoxPot;
  set_bvals_fun(left_x1,  no_op_VGfun);
  set_bvals_fun(right_x1, no_op_VGfun);

/* Allocate memory for remapped variables in ghost zones */

  Nx2 = pGrid->Nx2 + 2*nghost;
  if ((Flx = (Gas*)malloc(Nx2*sizeof(Gas))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");
  if ((RemapGas=(Gas**)calloc_2d_array(nghost, Nx2, sizeof(Gas))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");
  if ((pU=(Real**)malloc(Nx2*sizeof(Real*))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
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
  int Nx2;
  Real x1min, x1max, x2min, x2max;

  Omega = par_getd_def("problem","omega",1.0e-3);

/* Must recompute global variable Lx needed by BC routines */
  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  Lx = x1max - x1min;

  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  Ly = x2max - x2min;

  StaticGravPot = ShearingBoxPot;
  set_bvals_fun(left_x1,  no_op_VGfun);
  set_bvals_fun(right_x1, no_op_VGfun);

/* Allocate memory for remapped variables in ghost zones */

  Nx2 = pG->Nx2 + 2*nghost;
  if ((Flx = (Gas*)malloc(Nx2*sizeof(Gas))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");
  if ((RemapGas=(Gas**)calloc_2d_array(nghost, Nx2, sizeof(Gas))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");
  if ((pU=(Real**)malloc(Nx2*sizeof(Real*))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");

  return;
}

/* Get_user_expression computes dVy */
Gasfun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  return NULL;
}

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*------------------------------------------------------------------------------
 * shear_ix1_ox1() - shearing-sheet BCs in x1 for 3D sims, does the ix1/ox1
 *   boundaries simultaneously which is necessary in MPI parallel jobs since
 *   ix1 must receive data sent by ox1, and vice versa
 * This is a public function which is called by set_bvals() (inside a
 * SHEARING_BOX macro).  The hgb problem generator enrolls NoOp functions for
 * the x1 bval routines, so that set_bvals() uses MPI to handle the internal
 * boundaries between grids properly, and this routine to handle the shearing
 * sheet BCs.
 *
 * RemapGas and Flx are defined as global arrays, memory allocated in hgb
 *----------------------------------------------------------------------------*/

void shear_ix1_ox1(Grid *pG, int var_flag)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k,j_offset,j_remap;
  Real yshear, deltay, epsi, epso;

  if (var_flag == 1) return;  /* BC for Phi with self-gravity not set here */

/* Compute the distance the computational domain has sheared in y */
  yshear = 1.5*Omega*Lx*pG->time;

/* Split this into integer and fractional peices of the domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */
  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0   */
  j_offset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;
  epso = -epsi;

/*================ START REMAP FOR IX1 BOUNDARY =================*/
/*------- START BIG LOOP OVER k in ix1 remap (remap on ij-slices) ---------*/
  for (k=ks; k<=ke; k++) {

/* Copy data from ox1 side of grid into RemapGas array, using integer offset
 * to address appropriate elements at ox1 */

    for(j=js; j<=je; j++){
      j_remap = j - j_offset;
      if(j_remap < (int)js) j_remap += pG->Nx2;

/* Note RemapGas array has j as fastest index, for use as 1D pencils below */
      for(i=0; i<nghost; i++){
        RemapGas[i][j].d  = pG->U[k][j_remap][ie-nghost+1+i].d;
        RemapGas[i][j].M1 = pG->U[k][j_remap][ie-nghost+1+i].M1;
        RemapGas[i][j].M2 = pG->U[k][j_remap][ie-nghost+1+i].M2;
        RemapGas[i][j].M2 += 1.5*Omega*Lx*RemapGas[i][j].d;
        RemapGas[i][j].M3 = pG->U[k][j_remap][ie-nghost+1+i].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        RemapGas[i][j].E  = pG->U[k][j_remap][ie-nghost+1+i].E;
        RemapGas[i][j].E += (0.5/RemapGas[i][j].d)*
          (SQR(RemapGas[i][j].M2) - SQR(pG->U[k][j_remap][ie-nghost+1+i].M2));
#endif /* ADIABATIC */
      }
    }

/* Apply y-periodicity to RemapGas array */
    for(i=0; i<nghost; i++){
      for(j=1; j<=nghost; j++){
        RemapGas[i][js-j] = RemapGas[i][je+1-j];
        RemapGas[i][je+j] = RemapGas[i][js+j-1];
      }
    }

/* Compute "fluxes" of conserved quantities associated with fractional offset
 * of a cell (epsi), and perform conservative remap for ix1 ghost zones */

/* RemapGas is passed as 1D pencil in y-direction; Flx is returned as 1D arr */
    for(i=0; i<nghost; i++){
      CompRemapFlux(RemapGas[i],epsi,js,je+1,Flx);

      for(j=js; j<=je; j++){
        pG->U[k][j][is-nghost+i].d  = RemapGas[i][j].d -(Flx[j+1].d -Flx[j].d );
        pG->U[k][j][is-nghost+i].M1 = RemapGas[i][j].M1-(Flx[j+1].M1-Flx[j].M1);
        pG->U[k][j][is-nghost+i].M2 = RemapGas[i][j].M2-(Flx[j+1].M2-Flx[j].M2);
        pG->U[k][j][is-nghost+i].M3 = RemapGas[i][j].M3-(Flx[j+1].M3-Flx[j].M3);
#ifdef ADIABATIC
        pG->U[k][j][is-nghost+i].E  = RemapGas[i][j].E -(Flx[j+1].E -Flx[j].E );
#endif /* ADIABATIC */
      }
    }

  } /*----------- END BIG LOOP OVER k in ix1 remap ---------------*/

/*================ START REMAP FOR OX1 BOUNDARY =================*/
/*------- START BIG LOOP OVER k in ox1 remap (remap on ij-slices) ---------*/
  for (k=ks; k<=ke; k++) {

/* Copy data from ix1 side of grid into RemapGas array, using integer offset
 * to address appropriate elements at ix1 */

    for(j=js; j<=je; j++){
      j_remap = j + j_offset;
      if(j_remap > (int)je) j_remap -= pG->Nx2;

/* Note RemapGas array has j as fastest index, for use as 1D pencils below */
      for(i=0; i<nghost; i++){
        RemapGas[i][j].d  = pG->U[k][j_remap][is+i].d;
        RemapGas[i][j].M1 = pG->U[k][j_remap][is+i].M1;
        RemapGas[i][j].M2 = pG->U[k][j_remap][is+i].M2;
        RemapGas[i][j].M2 -= 1.5*Omega*Lx*RemapGas[i][j].d;
        RemapGas[i][j].M3 = pG->U[k][j_remap][is+i].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        RemapGas[i][j].E  = pG->U[k][j_remap][is+i].E;
        RemapGas[i][j].E += (0.5/RemapGas[i][j].d)*
          (SQR(RemapGas[i][j].M2) - SQR(pG->U[k][j_remap][is+i].M2));
#endif /* ADIABATIC */
      }
    }

/* Apply y-periodicity to RemapGas array */
    for(i=0; i<nghost; i++){
      for(j=1; j<=nghost; j++){
        RemapGas[i][js-j] = RemapGas[i][je+1-j];
        RemapGas[i][je+j] = RemapGas[i][js+j-1];
      }
    }

/* Compute "fluxes" of conserved quantities associated with fractional offset
 * of a cell (epso), and perform conservative remap */

/* RemapGas is passed as 1D pencil in y-direction; Flx is returned as 1D arr */
    for(i=0; i<nghost; i++){
      CompRemapFlux(RemapGas[i],epso,js,je+1,Flx);

      for(j=js; j<=je; j++){
        pG->U[k][j][ie+1+i].d  = RemapGas[i][j].d  - (Flx[j+1].d  - Flx[j].d );
        pG->U[k][j][ie+1+i].M1 = RemapGas[i][j].M1 - (Flx[j+1].M1 - Flx[j].M1);
        pG->U[k][j][ie+1+i].M2 = RemapGas[i][j].M2 - (Flx[j+1].M2 - Flx[j].M2);
        pG->U[k][j][ie+1+i].M3 = RemapGas[i][j].M3 - (Flx[j+1].M3 - Flx[j].M3);
#ifdef ADIABATIC
        pG->U[k][j][ie+1+i].E  = RemapGas[i][j].E  - (Flx[j+1].E  - Flx[j].E );
#endif /* ADIABATIC */
      }
    }

  } /*----------- END BIG LOOP OVER k in ox1 remap ---------------*/

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * CompRemapFlux: computes "fluxes" of conserved variables through y-interfaces
 * in remapped ghost zones for 3D shearing sheet BCs needed in shear_ix1_ox1().
 * Uses second-order reconstruction in conserved variables.
 * Input Arguments:
 *   U = Conserved variables at cell centers along 1-D slice
 *   eps = fraction of a cell to be remapped
 *   il,iu = lower and upper indices of zone centers in slice
 * U must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 *   Flux = fluxes of conserved variables at interfaces over [il:iu+1]
 */

void CompRemapFlux(const Gas U[], const Real eps,
                   const int il, const int iu, Gas Flux[])
{
  int i,n;
  Real dUc[NVAR],dUl[NVAR],dUr[NVAR],dUm[NVAR];
  Real lim_slope;
  Real *pFlux;

/*--- Step 1.
 * Set pointer to array elements of input conserved variables */

  for (i=il-2; i<=iu+2; i++) pU[i] = (Real*)&(U[i]);

/*--- Step 2.
 * Compute centered, L/R, and van Leer differences of conserved variables
 * Note we access contiguous array elements by indexing pointers for speed */

  for (i=il-1; i<=iu+1; i++) {
    for (n=0; n<(NVAR); n++) {
      dUc[n] = pU[i+1][n] - pU[i-1][n];
      dUl[n] = pU[i  ][n] - pU[i-1][n];
      dUr[n] = pU[i+1][n] - pU[i  ][n];
    }

/*--- Step 3.
 * Apply monotonicity constraints to characteristic projections */

    for (n=0; n<(NVAR); n++) {
      dUm[n] = 0.0;
      if (dUl[n]*dUr[n] > 0.0) {
        lim_slope = MIN(fabs(dUl[n]),fabs(dUr[n]));
        dUm[n] = SIGN(dUc[n])*MIN(0.5*fabs(dUc[n]),2.0*lim_slope);
      }
    }

/*--- Step 4.
 * Integrate linear interpolation function over eps */
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      pFlux = (Real *) &(Flux[i+1]);
      for (n=0; n<(NVAR); n++) {
        pFlux[n] = eps*(pU[i][n] + 0.5*(1.0 - eps)*dUm[n]);
      }

    } else {         /* eps always < 0 for outer i boundary */
      pFlux = (Real *) &(Flux[i]);
      for (n=0; n<(NVAR); n++) {
        pFlux[n] = eps*(pU[i][n] - 0.5*(1.0 + eps)*dUm[n]);
      }
    }

  }  /* end loop over [il-1,iu+1] */

  return;
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
 * ShearingBoxPot: includes vertical gravity if macro STRATIFIED_DISK is
 *   defined above.
 */

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3){
#ifdef STRATIFIED_DISK
  return 0.5*Omega*Omega*(x3*x3 - 3.0*x1*x1);
#else
  return -1.5*Omega*Omega*x1*x1;  
#endif
}

/*------------------------------------------------------------------------------
 * expr_dV2: computes delta(Vy) 
 */

static Real expr_dV2(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + 1.5*Omega*x1);
}

/*------------------------------------------------------------------------------
 * no_op_VGfun: replaces x1-bval routines
 */

static void no_op_VGfun(Grid *pGrid, int phi_flag)
{
  return;
}
