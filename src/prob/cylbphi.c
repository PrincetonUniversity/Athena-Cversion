#include "copyright.h"
/*==============================================================================
 * FILE: cylbphi.c
 *
 * A simple magnetostatic test of pressure balance using a B-field with uniform
 * phi-component.  
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

#if !defined MHD || !defined ADIABATIC
#error This problem only works for adiabatic MHD...
#endif 

static Real bphi0, omega, vz0, rho0, pgas0;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return 0.5*SQR(x1*omega) - (SQR(bphi0)/rho0)*log(x1);
              break;
    case 2:   return 0.5*SQR(x1*omega);
              break;
    default:  return 0.0;
  }
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return x1*SQR(omega) - SQR(bphi0)/(rho0*x1);
              break;
    case 2:   return x1*SQR(omega);
              break;
    default:  return 0.0;
  }
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
  Real x1,x2,x3;
  Real y1,y2,y3;
//   Real x1min, x1max, x2min, x2max, x3min, x3max; 
  Real Eint, Emag, Ekin;

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
    ath_error("[cylbphi]: This problem can only be run in 2D or 3D!\n");
    exit(EXIT_FAILURE);
  }

//   x1min = par_getd("grid","x1min");
//   x1max = par_getd("grid","x1max");
//   x2min = par_getd("grid","x2min");
//   x2max = par_getd("grid","x2max");
//   x3min = par_getd("grid","x3min");
//   x3max = par_getd("grid","x3max");

  omega  = par_getd("problem", "omega");
  vz0    = par_getd("problem", "vz0");
  bphi0  = par_getd("problem", "bphi0");
  rho0   = par_getd("problem", "rho0");
  pgas0  = par_getd("problem", "pgas0");
  iprob  = par_geti("problem", "iprob");


  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylbphi]: Error allocating memory for solution!\n");


  /* SET DENSITY, MOMENTUM, AND MAGNETIC FIELDS
     iprob = 1, CONSTANT B-PHI, CONSTANT PRESSURE
     iprob = 2, B-PHI GOES AS 1/R, CONSTANT PRESSURE
  */

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M1 = 0.0;
//         pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = vz0;
        switch (iprob) {
          case 1:   // CONSTANT B_phi
                    pG->B2i[k][j][i]   = bphi0;
                    pG->U[k][j][i].B2c = bphi0;
                    break;
          case 2:   // B_phi GOES AS 1/R
                    pG->B2i[k][j][i]   = bphi0/x1;
                    pG->U[k][j][i].B2c = bphi0/x1;
                    break;
          default:  printf("[cylbphi]:  Not an accepted problem number\n");
        }

        /* INITIALIZE TOTAL ENERGY */
        Eint = pgas0/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
                  + SQR(pG->U[k][j][i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2)
                  + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
        pG->U[k][j][i].E = Eint + Emag + Ekin;
      }
    }
  }


/* SAVE SOLUTION */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Soln[k][j][i] = pG->U[ks][j][i];
      }
    }
  }


  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_fun(left_x1,do_nothing_bc);
  set_bvals_fun(right_x1,do_nothing_bc);

  return;
}




/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets poEintr to expression for special output data
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

void Userwork_in_loop(Grid *pG, Domain *pDomain)
{
//   printf("Max divB = %1.10e\n", compute_div_b(pG));
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN EintGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(Grid *pG, Domain *pDomain)
{
  compute_l1_error("CylBPhi", pG, pDomain, Soln, 0);
}


