#include "copyright.h"
/*==============================================================================
 * FILE: cyl_blast.c
 *
 * PURPOSE: Problem generator for blast wave in cylindrical coords.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *      problem/radius = radius of field initial overpressured region
 *      problem/pamb   = ambient pressure
 *      problem/prat   = ratio of interior to ambient pressure
 *      problem/b0     = initial azimuthal magnetic field (units sqrt(Pamb))
 *      problem/omega  = initial azimuthal flow angular velocity
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

static Real radius,pamb,prat,b0,omega;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega);
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega);
}


/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real r0,phi0,x0,y0,z0,angle;
  Real x1,x2,x3,y1,y2,y3,x1i,x2i,x3i;
  Real x,y,z,pressure;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylblast]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylblast]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylblast]: Only (R,phi) can be used in 2D!\n");
  }

/* Read initial conditions */

  radius = par_getd("problem","radius");
  pamb = par_getd("problem","pamb");
  prat = par_getd("problem","prat");
  omega = par_getd("problem","omega");
#ifdef MHD
  b0 = par_getd("problem","b0");
#endif
  /* placement of center of blast */
  r0 = par_getd("problem","r0");
  phi0 = par_getd("problem","phi0");
  z0 = par_getd("problem","z0");

  angle = (PI/180.0)*par_getd("problem","angle");

  x0 = r0*cos(phi0);
  y0 = r0*sin(phi0);

/* set up uniform ambient medium with circular over-pressured region */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&x1,&x2,&x3);
	x1i = x1 - 0.5*pG->dx1;
	x2i = x2 - 0.5*pG->dx2;
	x3i = x3 - 0.5*pG->dx3;

        pG->U[k][j][i].d = 1.0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = 0.0;
#ifdef MHD
        // THIS IS A PLANAR MAGNETIC FIELD IN THE R-PHI PLANE
        pG->B1i[k][j][i] = b0*(cos(angle)*cos(x2)+sin(angle)*sin(x2)); 
        pG->B2i[k][j][i] = b0*(-cos(angle)*sin(x2i)+sin(angle)*cos(x2i));
        pG->B3i[k][j][i] = 0.0;
        pG->U[k][j][i].B1c = b0*(cos(angle)*cos(x2)+sin(angle)*sin(x2));
        pG->U[k][j][i].B2c = b0*(-cos(angle)*sin(x2)+sin(angle)*cos(x2));
        pG->U[k][j][i].B3c = 0.0;
#endif
        /* Cartesian position */
        x = x1*cos(x2);
        y = x1*sin(x2);
        z = x3;
        /* initialize total energy */
        pressure=pamb;
        if (SQR(x-x0) + SQR(y-y0) + SQR(z-z0) < SQR(radius)) {
          pressure=pamb*prat;
        }
#ifndef ISOTHERMAL
        pG->U[k][j][i].E = pressure/Gamma_1 
#ifdef MHD
          + 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
          + SQR(pG->U[k][j][i].B3c))
#endif
          + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
          + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#else
        if (SQR(x-x0) + SQR(y-y0) + SQR(z-z0) < SQR(radius)) {
          pG->U[k][j][i].d = prat;
        }
#endif /* ISOTHERMAL */

      }
    }
  }

  /* Enroll the gravitational function and radial BC */
  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;  
  set_bvals_mhd_fun(left_x1, do_nothing_bc);
  set_bvals_mhd_fun(right_x1,do_nothing_bc);
  set_bvals_mhd_fun(left_x2, do_nothing_bc);
  set_bvals_mhd_fun(right_x2,do_nothing_bc);

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
}
