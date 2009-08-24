#include "copyright.h"
/*==============================================================================
 * FILE: cyladvect.c
 *
 * A simple square-pulse advection test in cylindrical coordinates
 *
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "cyl.h"

static Real omega,rho0,bz0,vz0,Ptotal,amp,r0,phi0,rad;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega);
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega);
}

static Gas ***Soln=NULL;


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3; 
  Real x1,x2,x3,y1,y2,y3;
  Real x2min,x2max,R1,R2,PHI1,PHI2,alpha,X0,Y0;
  Real Eint,Emag,Ekin;
  Real **tmpl;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cyladvect]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cyladvect]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cyladvect]: Only (R,phi) can be used in 2D!\n");
  }

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cyladvect]: Error allocating memory for solution\n");

  if ((tmpl = (Real**)calloc_2d_array(nx2,nx1,sizeof(Real))) == NULL) {
    ath_error("[cyladvect]: Error allocating memory for advection template\n");
  }

  /* PARSE INPUT FILE */
  x2min = par_getd("grid", "x2min");
  x2max = par_getd("grid", "x2max");
  omega  = par_getd("problem", "omega");
  rho0   = par_getd("problem", "rho0");
  bz0    = par_getd("problem", "bz0");
  vz0    = par_getd("problem", "vz0");
  Ptotal = par_getd("problem", "Ptotal");
  amp    = par_getd("problem", "amp");
  r0     = par_getd("problem", "r0");
  phi0   = par_getd("problem", "phi0");
  rad    = par_getd("problem", "rad");
  iprob  = par_geti("problem", "iprob");

#ifndef MHD
  if (iprob==2) ath_error("[cyladvect]:  This problem can only be run in MHD!\n");
#endif

  /* INITIALIZE TEMPLATE FOR ADVECTION PROBLEM */
  R1   = r0   - 0.5*rad;
  R2   = r0   + 0.5*rad;
  PHI1 = phi0 - 0.5*rad/r0;
  PHI2 = phi0 + 0.5*rad/r0;

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      tmpl[j][i] = (Real)((x1>=R1) && (x1<=R2) && (x2>=PHI1) && (x2<=PHI2));
    }
  }

  /* GIVES ONE SINUSOIDAL PERIOD IN PHI */
  alpha = 2.0*PI/(x2max-x2min);
  X0 = r0*cos(phi0);
  Y0 = r0*sin(phi0);

  /* SET DENSITY AND PHI-MOMENTUM (MUST USE VOLUME CENTER) */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));

        switch(iprob) {
          case 1:  // SQUARE-PULSE, HYDRO ONLY
                   pG->U[k][j][i].d   = rho0*(1.0 + amp*tmpl[j][i]);
                   break;
          case 2:  // SQUARE-PULSE, UNIF. B_z IN PULSE
                   pG->U[k][j][i].d   = rho0*(1.0 + amp*tmpl[j][i]);
#ifdef MHD
                   pG->U[k][j][i].B3c = bz0*tmpl[j][i];
                   pG->B3i[k][j][i]   = bz0*tmpl[j][i];
#endif
                   break;
          case 3:  // SINE-WAVE, HYDRO ONLY
                   pG->U[k][j][i].d   = rho0*(1.0 + amp*sin(alpha*x2));
                   break;
          case 4:  // GAUSSIAN PULSE, HYDRO ONLY
                   R1 = sqrt(SQR(x1*cos(x2)-X0)+SQR(x1*sin(x2)-Y0));
                   pG->U[k][j][i].d   = rho0*(1.0+exp(-0.5*SQR(3.0*R1/rad)));
                   break;
          default:
                   ath_error("[cyladvect]:  Not an accepted problem number\n");
        }

        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz0;

#ifndef ISOTHERMAL
        // INITIALIZE TOTAL ENERGY
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
        Eint = (Ptotal-Emag)/Gamma_1;
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;

        pG->U[k][j][i].E = Eint + Emag + Ekin;
#endif /* ISOTHERMAL */

        /* SAVE SOLUTION */
        Soln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  free_2d_array((void**)tmpl);

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_mhd_fun(left_x1,do_nothing_bc);
  set_bvals_mhd_fun(right_x1,do_nothing_bc);

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

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
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
  compute_l1_error("CylAdvect", pGrid, pDomain, Soln, 0);
}
