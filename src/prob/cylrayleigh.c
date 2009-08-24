#define SEED 661979

#include "copyright.h"
/*==============================================================================
 * FILE: cylrayleigh.c
 *
 *  A simple magnetostatic test of the Rayleigh instability using 
 *  omega = 1/r^K.
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

static Real bphi0,omega0,rho0,pgas0,kappa,noise_level;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  Real omega = omega0/pow(x1,kappa);
  return 0.5*SQR(x1*omega)/(1.0-kappa);
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  Real omega = omega0/pow(x1,kappa);
  return x1*SQR(omega);
}


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3,myid;
  Real x1,x2,x3,R1,R2;
  Real r,noise,omega;
  Real Eint,Emag,Ekin;

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  il = is-nghost;  iu = ie+nghost;
  jl = js-nghost;  ju = je+nghost;
  if (ke-ks == 0) {
    kl = ks;  ku = ke; 
  } else {
    kl = ks-nghost;  ku = ke+nghost;
  }

  nx1 = iu-il+1;
  nx2 = ju-jl+1;
  nx3 = ku-kl+1;

  if ((nx2 == 1) && (nx3 == 1)) {
    ath_error("[cyl_rayleigh]: This problem can only be run in 2D or 3D!\n");
  }


  /* SEED THE RANDOM NUMBER GENERATOR */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  srand(SEED+myid);

  omega0      = par_getd("problem", "omega0");
  bphi0       = par_getd("problem", "bphi0");
  rho0        = par_getd("problem", "rho0");
  pgas0       = par_getd("problem", "pgas0");
  kappa       = par_getd("problem", "kappa");
  noise_level = par_getd("problem", "noise_level");


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));
	R1 = x1 - 0.5*pG->dx1;
	R2 = x1 + 0.5*pG->dx1;

        // RANDOM NUMBER BETWEEN 0 AND 1
        r = ((double) rand()/((double)RAND_MAX + 1.0)); 
	// RANDOM NUMBER BETWEEN +/- noise_level
        noise = noise_level*(2.0*r-1.0);

        pG->U[k][j][i].d   = rho0;
	omega = omega0/pow(x1,kappa);
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega;
//        pG->U[k][j][i].M2 = pG->U[k][j][i].d*omega0*(pow(R2,3-kappa)-pow(R1,3-kappa))/(x1*pG->dx1*(3-kappa));
        // NOW PERTURB M2
	if ((i>=is) && (i<=ie)) {
          pG->U[k][j][i].M2 *= (1.0 + noise);
        }

#ifdef MHD
        pG->U[k][j][i].B2c = bphi0/x1;
        pG->B2i[k][j][i]   = bphi0/x1;
#endif

#ifndef ISOTHERMAL
        Eint = pgas0/Gamma_1;
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1 ) + SQR(pG->U[k][j][i].M2 ) + SQR(pG->U[k][j][i].M3 ))/pG->U[k][j][i].d;
        pG->U[k][j][i].E = Eint + Emag + Ekin;
#endif
      }
    }
  }

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_fun(left_x1,do_nothing_bc);
  set_bvals_fun(right_x1,do_nothing_bc);
//  set_bvals_fun(left_x1,strict_outflow_ix1);
//  set_bvals_fun(right_x1,strict_outflow_ox1);

  return;
}



/*==============================================================================
 * PROBLEM USER FUNCTIONS:
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

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  omega0      = par_getd("problem", "omega0");
  bphi0       = par_getd("problem", "bphi0");
  rho0        = par_getd("problem", "rho0");
  pgas0       = par_getd("problem", "pgas0");
  kappa       = par_getd("problem", "kappa");

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_fun(left_x1,do_nothing_bc);
  set_bvals_fun(right_x1,do_nothing_bc);
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
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

void Userwork_in_loop(Grid *pG, Domain *pDomain)
{
//   printf("Max divB = %1.10e\n", compute_div_b(pG));
}

void Userwork_after_loop(Grid *pG, Domain *pDomain)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
