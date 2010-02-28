#include "copyright.h"
/*==============================================================================
 * FILE: planet-disk.c
 *
 * PURPOSE:  Problem generator for planet embedded in a disk, using the
 *   shearing sheet approximation.
 *
 * Code must be configured using --enable-shearing-box
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
 * ShearingBoxPot() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static Real Mp,Rsoft,Xplanet,Yplanet,Zplanet;
static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real den = 1.0, pres = 1.0e-6;
  static int frst=1;  /* flag so new history variables enrolled only once */

/* specify xy (r-phi) plane */
  ShBoxCoord = xy;

/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);

/* Compute field strength based on beta.  */
#ifdef ISOTHERMAL
  pres = Iso_csound2;
#endif

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize d, M, and P.  With FARGO do not initialize the background shear */

      pGrid->U[k][j][i].d  = den;
      pGrid->U[k][j][i].M1 = 0.0;
      pGrid->U[k][j][i].M2 = 0.0;
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= den*(qshear*Omega_0*x1);
#endif
      pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = pres/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/den;
#endif

    }
  }}

/* enroll gravitational potential function */

  StaticGravPot = ShearingBoxPot;

/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
    frst = 0;
  }

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif

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

/* 'problem_read_restart' must enroll gravity on restarts */

void problem_read_restart(MeshS *pM, FILE *fp)
{
/* Read Omega, and with viscosity and/or resistivity, read eta_Ohm and nu_V */

  Omega_0 = par_getd_def("problem","omega",1.0e-3);
  qshear  = par_getd_def("problem","qshear",1.5);
  Mp      = par_getd_def("problem","Mplanet",0.0);
  Xplanet = par_getd_def("problem","Xplanet",0.0);
  Yplanet = par_getd_def("problem","Yplanet",0.0);
  Zplanet = par_getd_def("problem","Zplanet",0.0);
  Rsoft   = par_getd_def("problem","Rsoft",0.1);
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

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * ShearingBoxPot:
 */

static Real ShearingBoxPot(const Real x1, const Real x2, const Real x3)
{
  Real rad,phi=0.0;
  rad = sqrt(SQR(x1-Xplanet) + SQR(x2-Yplanet) + SQR(x3-Zplanet));
  phi = -Mp/(rad+Rsoft);
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

#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = ShearingBoxPot(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */
