#include "copyright.h"
/*==============================================================================
 * FILE: hgb.c
 *
 * PURPOSE:  Problem generator for 3D shearing sheet.  Based on the initial
 *   conditions described in "Local Three-dimensional Magnetohydrodynamic
 *   Simulations of Accretion Disks" by Hawley, Gammie & Balbus, or HGB.
 *
 * Several different field configurations and perturbations are possible:
 *
 *  ifield = 1 - Bz=B0sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *  ipert = 3 - vortical shwave (hydro test)
 *
 * To run simulations of stratified disks (including vertical gravity),
 * un-comment the macro VERTICAL_GRAVITY below.
 *
 * This file also contains ShearingSheetBC(), a public function called by
 * set_bvals() which implements the 3D shearing sheet boundary conditions.
 *
 * REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* #define VERTICAL_GRAVITY */

/* Define structure which holds variables remapped by shearing sheet BCs */
typedef struct RVars_s{
  Real d,M1,M2,M3;
#ifndef ISOTHERMAL
  Real E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  Real s[NSCALARS];
#endif
#ifdef MHD
  Real B1c,B1i,B2i,B3i;
#endif /* MHD */
}RVars;

/* Define number of variables to be remapped */
#ifdef MHD
 enum {NREMAP = NVAR + 1};
#else 
 enum {NREMAP = NVAR};
#endif /* MHD */

/* prototypes for shearing sheet BC function (called by set_bvals) and 
 * remap function for Ey (called by integrator) */
void ShearingSheetBC(Grid *pG);
#ifdef MHD
void RemapEy(Grid *pG, Real ***emfy);
#endif /* MHD */

/* Remapped conserved quantities in ghost zones, and their fluxes */
static RVars ***RemapVar=NULL, *Flx=NULL;
static Real **pU=NULL;
#ifdef MHD
static Real **tEyiib=NULL, **tEyoib=NULL, *FlxEy=NULL;
#endif
#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
static Real **Uhalf=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * CompRemapFlux() - 2nd or 3rd order reconstruction for remap in ghost zones
 * ran2()          - random number generator from NR
 * ShearingBoxPot() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

void CompRemapFlux(const RVars U[], const Real eps,
     const int jl, const int ju, RVars Flux[]);
static double ran2(long int *idum);
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const Grid *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const Grid *pG, const int i, const int j, const int k);
static Real hst_E_total(const Grid *pG, const int i, const int j, const int k);
#ifdef MHD
void CompEyFlux(const Real *E, const Real eps,
                const int jinner, const int jouter, Real *FluxE);
static Real hst_Bx(const Grid *pG, const int i, const int j, const int k);
static Real hst_By(const Grid *pG, const int i, const int j, const int k);
static Real hst_Bz(const Grid *pG, const int i, const int j, const int k);
static Real hst_BxBy(const Grid *pG, const int i, const int j, const int k);
#endif /* MHD */


/* boxsize, made a global variable so can be accessed by bval, etc. routines */
static Real Lx,Ly;

/*=========================== PUBLIC FUNCTIONS =================================
 * Contains the usual, plus:
 * ShearingSheetBC() - shearing sheet BCs in 3D, called by set_bval().
 * EvolveEy()      - sets Ey in integrator to keep <Bz>=const. 
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,ipert,ifield,Nx2m,Nx3m;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1, x2, x3, x1min, x1max, x2min, x2max;
  Real den = 1.0, pres = 1.0e-6, rd, rp, rvx, rvy, rvz;
  Real beta,B0,kx,ky,amp;
  Real fkx,fky; /* wavenumber; only used for shwave tests */
  int nwx,nwy;  /* input number of waves per Lx, Ly -- only used for shwave */
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
  ky = 2.0*PI/Ly;

/* For shwave test, initialize wavenumber */
  nwx = par_geti_def("problem","nwx",1);
  nwy = par_geti_def("problem","nwy",1);
  fkx = kx*((double)nwx);  /* nxw should be input as -ve for leading wave */
  fky = ky*((double)nwy);

/* Read problem parameters.  Note Omega set to 10^{-3} by default */
  Omega = par_getd_def("problem","omega",1.0e-3);
  amp = par_getd("problem","amp");
  beta = par_getd("problem","beta");
  B0 = sqrt((double)(2.0*pres/beta));
  ifield = par_geti_def("problem","ifield", 1);
  ipert = par_geti_def("problem","ipert", 1);

/* Rescale amp to sound speed for ipert 2,3 */
#ifdef ADIABATIC
  if (ipert == 2 || ipert == 3) amp *= sqrt(Gamma*pres/den);
#else
  if (ipert == 2 || ipert == 3) amp *= Iso_csound;
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *  ipert = 3 - vortical shwave (hydro test)
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
        rvx = amp;
        rvy = 0.0;
        rvz = 0.0;
      }
      if (ipert == 3) {
        rp = pres;
        rd = den;
        rvx = amp*sin((double)(fkx*x1 + fky*x2));
        rvy = -amp*(fkx/fky)*sin((double)(fkx*x1 + fky*x2));
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

/* enroll gravitational potential function */

  StaticGravPot = ShearingBoxPot;

/* enroll new history variables */

  dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
  dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
  dump_history_enroll(hst_E_total, "<E + rho Phi>");
#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_By, "<By>");
  dump_history_enroll(hst_Bz, "<Bz>");
  dump_history_enroll(hst_BxBy, "<-Bx By>");
#endif /* MHD */

/* Allocate memory for remapped variables in ghost zones */

  Nx2m = pGrid->Nx2 + 2*nghost;
  Nx3m = pGrid->Nx2 + 2*nghost;

  if ((Flx = (RVars*)malloc(Nx2m*sizeof(RVars))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if((RemapVar=(RVars***)calloc_3d_array(Nx3m,nghost,Nx2m,sizeof(RVars)))==NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if ((pU=(Real**)malloc(Nx2m*sizeof(Real*))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

#ifdef MHD
  if ((tEyiib=(Real**)calloc_2d_array(Nx3m,Nx2m,sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if ((tEyoib=(Real**)calloc_2d_array(Nx3m,Nx2m,sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if ((FlxEy=(Real*)malloc(Nx2m*sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");
#endif /* MHD */

#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
  if ((Uhalf = (Real**)calloc_2d_array(nmax, NREMAP, sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");
#endif

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
  int Nx2m,Nx3m;
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

/* Allocate memory for remapped variables in ghost zones */

  Nx2m = pG->Nx2 + 2*nghost;
  Nx3m = pG->Nx2 + 2*nghost;

  if((Flx = (RVars*)malloc(Nx2m*sizeof(RVars))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");

  if((RemapVar=(RVars***)calloc_3d_array(Nx3m,nghost,Nx2m,sizeof(RVars)))==NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");

  if((pU=(Real**)malloc(Nx2m*sizeof(Real*))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");

#ifdef MHD
  if ((tEyiib=(Real**)calloc_2d_array(Nx3m,Nx2m,sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if ((tEyoib=(Real**)calloc_2d_array(Nx3m,Nx2m,sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");

  if ((FlxEy=(Real*)malloc(Nx2m*sizeof(Real))) == NULL)
    ath_error("[hgb]: malloc returned a NULL pointer\n");
#endif /* MHD */

#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
  if((Uhalf = (Real**)calloc_2d_array(Nx2m, NREMAP, sizeof(Real))) == NULL)
    ath_error("[read_restart]: malloc returned a NULL pointer\n");
#endif

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
 * ShearingSheetBC() - 3D shearing-sheet BCs in x1.  It applies a remap
 * in Y after the ghost cells have been set by the usual periodic BCs in X and
 * Y implemented in set_bvals.c
 *
 * This is a public function which is called by set_bvals() (inside a
 * SHEARING_BOX macro).
 *
 * RemapVar and Flx are defined as global arrays, memory allocated in problem
 *----------------------------------------------------------------------------*/

void ShearingSheetBC(Grid *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,ii,j,k,joffset,jremap;
#if (NSCALARS > 0)
  int n;
#endif
  Real yshear, deltay, epsi, epso, qa;

/* Compute the distance the computational domain has sheared in y */
  yshear = 1.5*Omega*Lx*pG->time;

/* Split this into integer and fractional peices of the domain in y.  Ignore
 * the integer piece because the Grid is periodic in y */
  deltay = fmod(yshear, Ly);

/* further decompose the fractional peice into integer and fractional pieces of
 * a grid cell.  Note 0.0 <= epsi < 1.0   */
  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;
  epso = -epsi;

/*=============== START REMAP ON IX1 BOUNDARY  ================*/
/* When this fun is called, the ix1 ghost zones have been set by periodic BCs
 * in X and Y.  This routine now applies the remap. */

/* Copy data into RemapVar array.  Note i and j indices are switched. */

  qa = 1.5*Omega*Lx;
  for(k=ks; k<=ke+1; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = is-nghost+i;
        RemapVar[k][i][j].d  = pG->U[k][j][ii].d;
        RemapVar[k][i][j].M1 = pG->U[k][j][ii].M1;
        RemapVar[k][i][j].M2 = pG->U[k][j][ii].M2 + qa*pG->U[k][j][ii].d;
        RemapVar[k][i][j].M3 = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        RemapVar[k][i][j].E  = pG->U[k][j][ii].E + (0.5/RemapVar[k][i][j].d)*
          (SQR(RemapVar[k][i][j].M2) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        RemapVar[k][i][j].B1c = pG->U[k][j][ii].B1c;
        RemapVar[k][i][j].B1i = pG->B1i[k][j][ii];
        RemapVar[k][i][j].B2i = pG->B2i[k][j][ii];
        RemapVar[k][i][j].B3i = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) RemapVar[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/* Remap over (ks:ke+1)(js:je+1)(is-nghost:is-1) */

  for(k=ks; k<=ke+1; k++) {
    for(i=0; i<nghost; i++){
      CompRemapFlux(RemapVar[k][i],epsi,js,je+2,Flx);

      for(j=js; j<=je+1; j++){
        RemapVar[k][i][j].d  -= (Flx[j+1].d  - Flx[j].d );
        RemapVar[k][i][j].M1 -= (Flx[j+1].M1 - Flx[j].M1);
        RemapVar[k][i][j].M2 -= (Flx[j+1].M2 - Flx[j].M2);
        RemapVar[k][i][j].M3 -= (Flx[j+1].M3 - Flx[j].M3);
#ifdef ADIABATIC
        RemapVar[k][i][j].E  -= (Flx[j+1].E  - Flx[j].E );
#endif /* ADIABATIC */
#ifdef MHD
        RemapVar[k][i][j].B1c -= (Flx[j+1].B1c - Flx[j].B1c);
        RemapVar[k][i][j].B1i -= (Flx[j+1].B1i - Flx[j].B1i);
        RemapVar[k][i][j].B2i -= (Flx[j+1].B2i - Flx[j].B2i);
        RemapVar[k][i][j].B3i -= (Flx[j+1].B3i - Flx[j].B3i);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          RemapVar[k][i][j].s[n] -= (Flx[j+1].s[n] - Flx[j].s[n]);
        }
#endif
      }
    }
  }

/* Copy remapped variables back into appropriate ghost cells with integer
 * shift.  There are two different cases:
 *   (1) no MPI decomposition in Y, just do copy with shift
 *   (2) MPI decomposition in Y, do circular shift in Y using MPI calls
 */

/* CASE 1: no MPI decomposition in Y ---------------------------------------- */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      jremap = j - joffset;
      if (jremap < (int)js) jremap += pG->Nx2;
      for(i=0; i<nghost; i++){
        pG->U[k][j][is-nghost+i].d  = RemapVar[k][i][jremap].d ;
        pG->U[k][j][is-nghost+i].M1 = RemapVar[k][i][jremap].M1;
        pG->U[k][j][is-nghost+i].M2 = RemapVar[k][i][jremap].M2;
        pG->U[k][j][is-nghost+i].M3 = RemapVar[k][i][jremap].M3;
#ifdef ADIABATIC
        pG->U[k][j][is-nghost+i].E  = RemapVar[k][i][jremap].E ;
#endif /* ADIABATIC */
#ifdef MHD
        pG->U[k][j][is-nghost+i].B1c = RemapVar[k][i][jremap].B1c;
        pG->B1i[k][j][is-nghost+i] = RemapVar[k][i][jremap].B1i;
        pG->B2i[k][j][is-nghost+i] = RemapVar[k][i][jremap].B2i;
        pG->B3i[k][j][is-nghost+i] = RemapVar[k][i][jremap].B3i;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][is-nghost+i].s[n] = RemapVar[k][i][jremap].s[n];
        }
#endif
      }
    }

/* Copy the face-centered B2 component of the field at j=je+1 */
#ifdef MHD
    jremap = je+1 - joffset;
    if (jremap < (int)js) jremap += pG->Nx2;
    for(i=0; i<nghost; i++){
      pG->B2i[k][je+1][is-nghost+i] = RemapVar[k][i][jremap].B2i;
    }
#endif /* MHD */
  }

/* Copy the face-centered B3 component of the field at k=ke+1 */
#ifdef MHD
  for(j=js; j<=je; j++){
    jremap = j - joffset;
    if (jremap < (int)js) jremap += pG->Nx2;
    for(i=0; i<nghost; i++){
      pG->B3i[ke+1][j][is-nghost+i] = RemapVar[ke+1][i][jremap].B3i;
    }
  }
#endif /* MHD */

/* compute cell-centered B as average of remapped face centered B, except B1 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=is-nghost; i<is; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/* periodic BCs in Y (compare to periodic_ix2() and perdiodic_ox2() in
 * set_bavls.c) */

  for(k=ks; k<=ke; k++) {
    for(j=1; j<=nghost; j++){
      for(i=is-nghost; i<is; i++){
        pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
        pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
        pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
        pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
        pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

        pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
        pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
        pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
      }
    }
  }
#ifdef MHD
  for (j=1; j<=nghost; j++) {
    for (i=is-nghost; i<is; i++) {
      pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
      pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
    }
  }
#endif /* MHD */

/* CASE 2: MPI decomposition in Y -- requires cshift in Y --------------------*/


/*=============== START REMAP FOR OX1 BOUNDARY  ================*/
/* When this fun is called, the ox1 ghost zones have been set by periodic BCs
 * in X and Y.  This routine now applies the remap. */

/* Copy data into RemapVar array.  Note i and j indices are switched. */

  qa = 1.5*Omega*Lx;
  for (k=ks; k<=ke+1; k++) {
    for(j=js-nghost; j<=je+nghost; j++){
      for(i=0; i<nghost; i++){
        ii = ie+1+i;
        RemapVar[k][i][j].d  = pG->U[k][j][ii].d;
        RemapVar[k][i][j].M1 = pG->U[k][j][ii].M1;
        RemapVar[k][i][j].M2 = pG->U[k][j][ii].M2 - qa*pG->U[k][j][ii].d;
        RemapVar[k][i][j].M3 = pG->U[k][j][ii].M3;
#ifdef ADIABATIC
/* No change in the internal energy */
        RemapVar[k][i][j].E  = pG->U[k][j][ii].E + (0.5/RemapVar[k][i][j].d)*
          (SQR(RemapVar[k][i][j].M2) - SQR(pG->U[k][j][ii].M2));
#endif /* ADIABATIC */
#ifdef MHD
        RemapVar[k][i][j].B1c = pG->U[k][j][ii].B1c;
        RemapVar[k][i][j].B1i = pG->B1i[k][j][ii];
        RemapVar[k][i][j].B2i = pG->B2i[k][j][ii];
        RemapVar[k][i][j].B3i = pG->B3i[k][j][ii];
#endif /* MHD */
#if (NSCALARS > 0)
        for(n=0; n<NSCALARS; n++) RemapVar[k][i][j].s[n] = pG->U[k][j][ii].s[n];
#endif
      }
    }
  }

/* Remap over (ks:ke+1)(js:je+1)(ie+1:ie+nghost).
 * Note that remap of B1c is needed at i=ie+nghost */

  for(k=ks; k<=ke+1; k++) {
    for(i=0; i<nghost; i++){
      CompRemapFlux(RemapVar[k][i],epso,js,je+2,Flx);

      for(j=js; j<=je+1; j++){
        RemapVar[k][i][j].d  -= (Flx[j+1].d  - Flx[j].d );
        RemapVar[k][i][j].M1 -= (Flx[j+1].M1 - Flx[j].M1);
        RemapVar[k][i][j].M2 -= (Flx[j+1].M2 - Flx[j].M2);
        RemapVar[k][i][j].M3 -= (Flx[j+1].M3 - Flx[j].M3);
#ifdef ADIABATIC
        RemapVar[k][i][j].E  -= (Flx[j+1].E  - Flx[j].E );
#endif /* ADIABATIC */
#ifdef MHD
        RemapVar[k][i][j].B1c -= (Flx[j+1].B1c - Flx[j].B1c);
        RemapVar[k][i][j].B1i -= (Flx[j+1].B1i - Flx[j].B1i);
        RemapVar[k][i][j].B2i -= (Flx[j+1].B2i - Flx[j].B2i);
        RemapVar[k][i][j].B3i -= (Flx[j+1].B3i - Flx[j].B3i);
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          RemapVar[k][i][j].s[n] -= (Flx[j+1].s[n] - Flx[j].s[n]);
        }
#endif
      }
    }
  }

/* Copy remapped variables back into appropriate ghost cells with integer
 * shift.  There are two different cases:
 *   (1) no MPI decomposition in Y, just do copy with shift
 *   (2) MPI decomposition in Y, do circular shift in Y using MPI calls
 */

/* CASE 1: no MPI decomposition in Y ---------------------------------------- */

  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++){
      jremap = j + joffset;
      if(jremap > (int)je) jremap -= pG->Nx2;
      for(i=0; i<nghost; i++){
        pG->U[k][j][ie+1+i].d  = RemapVar[k][i][jremap].d ;
        pG->U[k][j][ie+1+i].M1 = RemapVar[k][i][jremap].M1;
        pG->U[k][j][ie+1+i].M2 = RemapVar[k][i][jremap].M2;
        pG->U[k][j][ie+1+i].M3 = RemapVar[k][i][jremap].M3;
#ifdef ADIABATIC
        pG->U[k][j][ie+1+i].E  = RemapVar[k][i][jremap].E ;
#endif /* ADIABATIC */
#ifdef MHD
        pG->U[k][j][ie+1+i].B1c = RemapVar[k][i][jremap].B1c;
        pG->B2i[k][j][ie+1+i] = RemapVar[k][i][jremap].B2i;
        pG->B3i[k][j][ie+1+i] = RemapVar[k][i][jremap].B3i;
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          pG->U[k][j][ie+1+i].s[n] = RemapVar[k][i][jremap].s[n];
        }
#endif
      }

/* B1i is not remapped at i=ie+1, this requires a special loop that skips i=0 */
#ifdef MHD
      for(i=1; i<nghost; i++){
        pG->B1i[k][j][ie+1+i] = RemapVar[k][i][jremap].B1i;
      }
#endif /* MHD */
    }

/* Copy the face-centered B2 component of the field at j=je+1 */
#ifdef MHD
    jremap = je+1 + joffset;
    if(jremap > (int)je) jremap -= pG->Nx2;
    for(i=0; i<nghost; i++){
      pG->B2i[k][je+1][ie+1+i] = RemapVar[k][i][jremap].B2i;
    }
#endif /* MHD */
  }

/* Copy the face-centered B3 component of the field at k=ke+1 */
#ifdef MHD
  for(j=js; j<=je; j++){
    jremap = j + joffset;
    if(jremap > (int)je) jremap -= pG->Nx2;
    for(i=0; i<nghost; i++){
      pG->B3i[ke+1][j][ie+1+i] = RemapVar[ke+1][i][jremap].B3i;
    }
  }
#endif /* MHD */

/* Compute cell-centered B as average of remapped face centered B, except B1 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for(i=ie+1; i<=ie+nghost; i++){
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/* periodic BCs in Y (compare to periodic_ix2() and perdiodic_ox2() in
 * set_bavls.c) */

  for (k=ks; k<=ke; k++) {
    for(j=1; j<=nghost; j++){
      for(i=ie+1; i<=ie+nghost; i++){
        pG->U[k][js-j][i] = pG->U[k][je-(j-1)][i];
        pG->U[k][je+j][i] = pG->U[k][js+(j-1)][i];
#ifdef MHD
        pG->B1i[k][js-j][i] = pG->B1i[k][je-(j-1)][i];
        pG->B2i[k][js-j][i] = pG->B2i[k][je-(j-1)][i];
        pG->B3i[k][js-j][i] = pG->B3i[k][je-(j-1)][i];

        pG->B1i[k][je+j][i] = pG->B1i[k][js+(j-1)][i];
        pG->B2i[k][je+j][i] = pG->B2i[k][js+(j-1)][i];
        pG->B3i[k][je+j][i] = pG->B3i[k][js+(j-1)][i];
#endif /* MHD */
      }
    }
  }
#ifdef MHD
  for (j=1; j<=nghost; j++) {
    for (i=ie+1; i<=ie+nghost; i++) {
      pG->B3i[ke+1][js-j][i] = pG->B3i[ke+1][je-(j-1)][i];
      pG->B3i[ke+1][je+j][i] = pG->B3i[ke+1][js+(j-1)][i];
    }
  }
#endif /* MHD */

/* CASE 2: MPI decomposition in Y -- requires cshift in Y --------------------*/

  return;
}

/*------------------------------------------------------------------------------
 * RemapEy() - Remaps Ey at is and ie+1 due to background shear, and then
 * averages remapped and original field.  This guarantees the sums of Ey
 * along the x1 boundaries at is and ie+1 are identical -- thus net Bz is
 * conserved
 *
 * This is a public function which is called by integrator (inside a
 * SHEARING_BOX macro).
 *----------------------------------------------------------------------------*/

#ifdef MHD
void RemapEy(Grid *pG, Real ***emfy)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int j,k,joffset,jremap;
  Real yshear, deltay, epsi, epso;
  Real sum;

/* Compute the distance the computational domain has sheared in y in integer
 * and fractional pieces of a cell.  Same code as in ShearingSheetBC()  */

  yshear = 1.5*Omega*Lx*pG->time;  /* SHOULD THIS INCLUDE 1/2 dt **/
  deltay = fmod(yshear, Ly);
  joffset = (int)(deltay/pG->dx2);
  epsi = (fmod(deltay,pG->dx2))/pG->dx2;
  epso = -epsi;

/* Copy Ey at inner-i and outer-i boundary into temporary arrays, using
 * periodic BC in x1 */

  for(k=ks; k<=ke+1; k++) {
    for(j=js; j<=je; j++){
      tEyiib[k][j] = emfy[k][j][ie+1];
      tEyoib[k][j] = emfy[k][j][is  ];
    }
  }

/* Apply periodic BC in x2 to remporary arrays */

  for(k=ks; k<=ke+1; k++) {
    for(j=1; j<=nghost; j++){
      tEyiib[k][js-j] = tEyiib[k][je-(j-1)];
      tEyoib[k][js-j] = tEyoib[k][je-(j-1)];
      tEyiib[k][je+j] = tEyiib[k][js+(j-1)];
      tEyoib[k][je+j] = tEyoib[k][js+(j-1)];
    }
  }

/*============== Remap of Ey ix1 ===================*/
/* Remap over (ks:ke+1)(js:je) at is */

  for(k=ks; k<=ke+1; k++) {
    CompEyFlux(tEyiib[k],epsi,js,je+1,FlxEy);
    for(j=js; j<=je; j++){
      tEyiib[k][j] -= (FlxEy[j+1] - FlxEy[j]);
    }
  }

/* Average remapped Ey back into appropriate ghost cells with integer
 * shift.  There are two different cases:
 *   (1) no MPI decomposition in Y, just do copy with shift
 *   (2) MPI decomposition in Y, do circular shift in Y using MPI calls
 */

/* CASE 1: no MPI decomposition in Y ---------------------------------------- */

  for(k=ks; k<=ke+1; k++) {
    for(j=js; j<=je; j++){
      jremap = j - joffset;
      if (jremap < (int)js) jremap += pG->Nx2;
      emfy[k][j][is]  = 0.5*(emfy[k][j][is] + tEyiib[k][jremap]);
    }
  }

/* CASE 2: MPI decomposition in Y -- requires cshift in Y --------------------*/


/*============== Remap of Ey ox1 ===================*/
/* Remap over (ks:ke+1)(js:je) at ie+1 */

  for(k=ks; k<=ke+1; k++) {
    CompEyFlux(tEyoib[k],epso,js,je+1,FlxEy);
    for(j=js; j<=je; j++){
      tEyoib[k][j] -= (FlxEy[j+1] - FlxEy[j]);
    }
  }

/* Average remapped Ey back into appropriate ghost cells with integer
 * shift.  There are two different cases:
 *   (1) no MPI decomposition in Y, just do copy with shift
 *   (2) MPI decomposition in Y, do circular shift in Y using MPI calls
 */

/* CASE 1: no MPI decomposition in Y ---------------------------------------- */

  for(k=ks; k<=ke+1; k++) {
    for(j=js; j<=je; j++){
      jremap = j + joffset;
      if (jremap > (int)je) jremap -= pG->Nx2;
      emfy[k][j][ie+1]  = 0.5*(emfy[k][j][ie+1] + tEyoib[k][jremap]);
    }
  }

/* CASE 2: MPI decomposition in Y -- requires cshift in Y --------------------*/


  return;
}
#endif /* MHD */

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * CompRemapFlux: computes "fluxes" of conserved variables through y-interfaces
 * in ShearingSheetBC().  
 * Input Arguments:
 *   U = Conserved variables at cell centers along 1-D slice
 *   eps = fraction of a cell to be remapped
 *   il,iu = lower and upper indices of zone centers in slice
 * Output Arguments:
 *   Flux = fluxes of conserved variables at interfaces over [il:iu+1]
 */


/* SECOND ORDER REMAP: piecewise linear reconstruction and min/mod limiters
 * U must be initialized over [il-2:iu+2] */

#ifdef SECOND_ORDER
void CompRemapFlux(const RVars U[], const Real eps,
                   const int jl, const int ju, RVars Flux[])
{
  int j,n;
  Real dUc[NREMAP],dUl[NREMAP],dUr[NREMAP],dUm[NREMAP];
  Real lim_slope;
  Real *pFlux;

/*--- Step 1.
 * Set pointer to array elements of input conserved variables */

  for (j=jl-2; j<=ju+2; j++) pU[j] = (Real*)&(U[j]);

/*--- Step 2.
 * Compute centered, L/R, and van Leer differences of conserved variables
 * Note we access contiguous array elements by indexing pointers for speed */

  for (j=jl-1; j<=ju+1; j++) {
    for (n=0; n<(NREMAP); n++) {
      dUc[n] = pU[j+1][n] - pU[j-1][n];
      dUl[n] = pU[j  ][n] - pU[j-1][n];
      dUr[n] = pU[j+1][n] - pU[j  ][n];
    }

/*--- Step 3.
 * Apply monotonicity constraints */

    for (n=0; n<(NREMAP); n++) {
      dUm[n] = 0.0;
      if (dUl[n]*dUr[n] > 0.0) {
        lim_slope = MIN(fabs(dUl[n]),fabs(dUr[n]));
        dUm[n] = SIGN(dUc[n])*MIN(0.5*fabs(dUc[n]),2.0*lim_slope);
      }
    }

/*--- Step 4.
 * Integrate linear interpolation function over eps */
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      pFlux = (Real *) &(Flux[j+1]);
      for (n=0; n<(NREMAP); n++) {
        pFlux[n] = eps*(pU[j][n] + 0.5*(1.0 - eps)*dUm[n]);
      }

    } else {         /* eps always < 0 for outer i boundary */
      pFlux = (Real *) &(Flux[j]);
      for (n=0; n<(NREMAP); n++) {
        pFlux[n] = eps*(pU[j][n] - 0.5*(1.0 + eps)*dUm[n]);
      }
    }

  }  /* end loop over [jl-1,ju+1] */

  return;
}

/*------------------------------------------------------------------------------
 * CompEyFlux(): second order reconstruction for Ey, nearly identical to
 *  CompeRemapFlux() above.  jinner/jouter are range of indices over which flux
 *  is needed in calling function */

#ifdef MHD
void CompEyFlux(const Real *E, const Real eps,
                const int jinner, const int jouter, Real *FluxE)
{
  int j,jl,ju;
  Real dEc,dEl,dEr,dEm,lim_slope;

/* jinner,jouter are index range over which flux must be returned.  Set loop
 * limits depending on direction of upwind differences  */

  if (eps > 0.0) { /* eps always > 0 for inner i boundary */
    jl = jinner-1;
    ju = jouter-1;
  } else {         /* eps always < 0 for outer i boundary */
    jl = jinner;
    ju = jouter;
  }

  for (j=jl; j<=ju; j++) {
      dEc = E[j+1] - E[j-1];
      dEl = E[j  ] - E[j-1];
      dEr = E[j+1] - E[j  ];

      dEm = 0.0;
      if (dEl*dEr > 0.0) {
        lim_slope = MIN(fabs(dEl),fabs(dEr));
        dEm = SIGN(dEc)*MIN(0.5*fabs(dEc),2.0*lim_slope);
      }
 
    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      FluxE[j+1] = eps*(E[j] + 0.5*(1.0 - eps)*dEm);
    } else {         /* eps always < 0 for outer i boundary */
      FluxE[j  ] = eps*(E[j] - 0.5*(1.0 + eps)*dEm);
    }
  }

  return;
}
#endif /* MHD */

#endif /* SECOND_ORDER */

/* THIRD ORDER REMAP: Colella & Sekora extremum preserving algorithm (PPME)
 * U must be initialized over [il-3:iu+3] */

#if defined(THIRD_ORDER) || defined(THIRD_ORDER_EXTREMA_PRESERVING)
void CompRemapFlux(const Gas U[], const Real eps,
                   const int jl, const int ju, Gas Flux[])
{
  int j,n;
  Real lim_slope,qa,qb,qc,qx;
  Real d2Uc[NREMAP],d2Ul[NREMAP],d2Ur[NREMAP],d2U [NREMAP],d2Ulim[NREMAP];
  Real Ulv[NREMAP],Urv[NREMAP],dU[NREMAP],U6[NREMAP];
  Real *pFlux;

/*--- Step 1.
 * Set pointer to array elements of input conserved variables */

  for (i=il-3; i<=iu+3; i++) pU[i] = (Real*)&(U[i]);

/*--- Step 2. 
 * Compute interface states (CS eqns 12-15) over entire 1D pencil.  Using usual
 * Athena notation that index i for face-centered quantities denotes L-edge
 * (interface i-1/2), then Uhalf[i] = U[i-1/2]. */

  for (i=il-1; i<=iu+2; i++) {
    for (n=0; n<(NREMAP); n++) {
      Uhalf[i][n]=(7.0*(pU[i-1][n]+pU[i][n]) - (pU[i-2][n]+pU[i+1][n]))/12.0;
    }
    for (n=0; n<(NREMAP); n++) {
      d2Uc[n] = 3.0*(pU[i-1][n] - 2.0*Uhalf[i][n] + pU[i][n]);
      d2Ul[n] = (pU[i-2][n] - 2.0*pU[i-1][n] + pU[i  ][n]);
      d2Ur[n] = (pU[i-1][n] - 2.0*pU[i  ][n] + pU[i+1][n]);
      d2Ulim[n] = 0.0;
      lim_slope = MIN(fabs(d2Ul[n]),fabs(d2Ur[n]));
      if (d2Uc[n] > 0.0 && d2Ul[n] > 0.0 && d2Ur[n] > 0.0) {
        d2Ulim[n] = SIGN(d2Uc[n])*MIN(1.25*lim_slope,fabs(d2Uc[n]));
      }
      if (d2Uc[n] < 0.0 && d2Ul[n] < 0.0 && d2Ur[n] < 0.0) {
        d2Ulim[n] = SIGN(d2Uc[n])*MIN(1.25*lim_slope,fabs(d2Uc[n]));
      }
    }
    for (n=0; n<(NREMAP); n++) {
      Uhalf[i][n] = 0.5*((pU[i-1][n]+pU[i][n]) - d2Ulim[n]/3.0);
    }
  }

/*--- Step 3.
 * Compute L/R values
 * Ulv = U at left  side of cell-center = U[i-1/2] = a_{j,-} in CS
 * Urv = U at right side of cell-center = U[i+1/2] = a_{j,+} in CS
 */

  for (i=il-1; i<=iu+1; i++) {
    for (n=0; n<(NREMAP); n++) {
      Ulv[n] = Uhalf[i  ][n];
      Urv[n] = Uhalf[i+1][n];
    }

/*--- Step 4.
 * Construct parabolic interpolant (CS eqn 16-19) */

    for (n=0; n<(NREMAP); n++) {
      qa = (Urv[n]-pU[i][n])*(pU[i][n]-Ulv[n]);
      qb = (pU[i-1][n]-pU[i][n])*(pU[i][n]-pU[i+1][n]);
      if (qa <= 0.0 && qb <= 0.0) {
        qc = 6.0*(pU[i][n] - 0.5*(Ulv[n]+Urv[n]));
        d2U [n] = -2.0*qc;
        d2Uc[n] = (pU[i-1][n] - 2.0*pU[i  ][n] + pU[i+1][n]);
        d2Ul[n] = (pU[i-2][n] - 2.0*pU[i-1][n] + pU[i  ][n]);
        d2Ur[n] = (pU[i  ][n] - 2.0*pU[i+1][n] + pU[i+2][n]);
        d2Ulim[n] = 0.0;
        lim_slope = MIN(fabs(d2Ul[n]),fabs(d2Ur[n]));
        lim_slope = MIN(fabs(d2Uc[n]),lim_slope);
        if (d2Uc[n] > 0.0 && d2Ul[n] > 0.0 && d2Ur[n] > 0.0 && d2U[n] > 0.0) {
          d2Ulim[n] = SIGN(d2U[n])*MIN(1.25*lim_slope,fabs(d2U[n]));
        }
        if (d2Uc[n] < 0.0 && d2Ul[n] < 0.0 && d2Ur[n] < 0.0 && d2U[n] < 0.0) {
          d2Ulim[n] = SIGN(d2U[n])*MIN(1.25*lim_slope,fabs(d2U[n]));
        }
        if (d2U[n] == 0.0) {
          Ulv[n] = pU[i][n];
          Urv[n] = pU[i][n];
        } else {
          Ulv[n] = pU[i][n] + (Ulv[n] - pU[i][n])*d2Ulim[n]/d2U[n];
          Urv[n] = pU[i][n] + (Urv[n] - pU[i][n])*d2Ulim[n]/d2U[n];
        }
      }
    }

/*--- Step 5.
 * Monotonize again (CW eqn 1.10), ensure they lie between neighboring
 * cell-centered vals */

    for (n=0; n<(NREMAP); n++) {
      qa = (Urv[n]-pU[i][n])*(pU[i][n]-Ulv[n]);
      qb = Urv[n]-Ulv[n];
      qc = 6.0*(pU[i][n] - 0.5*(Ulv[n]+Urv[n]));
      if (qa <= 0.0) {
        Ulv[n] = pU[i][n];
        Urv[n] = pU[i][n];
      } else if ((qb*qc) > (qb*qb)) {
        Ulv[n] = 3.0*pU[i][n] - 2.0*Urv[n];
      } else if ((qb*qc) < -(qb*qb)) {
        Urv[n] = 3.0*pU[i][n] - 2.0*Ulv[n];
      }
    }

/*
    for (n=0; n<(NREMAP); n++) {
      Ulv[n] = MAX(MIN(pU[i][n],pU[i-1][n]),Ulv[n]);
      Ulv[n] = MIN(MAX(pU[i][n],pU[i-1][n]),Ulv[n]);
      Urv[n] = MAX(MIN(pU[i][n],pU[i+1][n]),Urv[n]);
      Urv[n] = MIN(MAX(pU[i][n],pU[i+1][n]),Urv[n]);
    }
*/

/*--- Step 6.
 * Compute coefficients of interpolation parabolae (CW eqn 1.5) */

    for (n=0; n<(NREMAP); n++) {
      dU[n] = Urv[n] - Ulv[n];
      U6[n] = 6.0*(pU[i][n] - 0.5*(Ulv[n] + Urv[n]));
    }

/*--- Step 7.
 * Integrate parabolic interpolation function over eps */

    if (eps > 0.0) { /* eps always > 0 for inner i boundary */
      pFlux = (Real *) &(Flux[i+1]);
      qx = TWO_3RDS*eps;
      for (n=0; n<(NREMAP); n++) {
        pFlux[n] = eps*(Urv[n] - 0.75*qx*(dU[n] - (1.0 - qx)*U6[n]));
      }

    } else {         /* eps always < 0 for outer i boundary */
      pFlux = (Real *) &(Flux[i]);
      qx = -TWO_3RDS*eps;
      for (n=0; n<(NREMAP); n++) {
        pFlux[n] = eps*(Ulv[n] + 0.75*qx*(dU[n] + (1.0 - qx)*U6[n]));
      }
    }
  }

  return;
}
#endif /* THIRD_ORDER_EXTREMA_PRESERVING */

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
 * ShearingBoxPot: includes vertical gravity if macro VERTICAL_GRAVITY is
 *   defined above.
 */

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3){
#ifdef VERTICAL_GRAVITY
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
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_rho_Vx_dVy(const Grid *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d + 1.5*Omega*x1);
}

static Real hst_rho_dVy2(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + 1.5*Omega*x1);
  return pG->U[k][j][i].d*dVy*dVy;
}

static Real hst_E_total(const Grid *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = ShearingBoxPot(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}

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

