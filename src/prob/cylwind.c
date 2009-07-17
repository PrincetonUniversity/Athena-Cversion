#include "copyright.h"
/*==============================================================================
 * FILE: cylwind.c
 *
 *  A test of rotationless Bondi accretion (Parker wind) problem.  
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

static Real b0,lambda_s;
static int iprob;

Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -1.0/x1;
}

Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return 1.0/SQR(x1);
}

Real myfunc(const Real x, const Real v);

static Gas ***Soln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k,converged;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3; 
  Real x1,x2,x3,a,b;
  Real xs,vs,v,p0,pgas,beta;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylwind]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylwind]: Only R can be used in 1D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylwind]: Only (R,phi) can be used in 2D!\n");
  }

  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylwind]: Error allocating memory for solution\n");

  b0     = par_getd("problem", "b0");
  iprob  = par_geti("problem", "iprob");

  beta = 2*Gamma_1/(Gamma+1);
  xs = (3.0-Gamma)/2.0;
  lambda_s = pow(xs,(beta-1)/beta);
//   vs = pow(lambda_s/xs,beta/2);
  vs = sqrt(2.0/(3.0-Gamma));
  printf("xs = %f, \tlambda_s = %f, \tvs = %f\n", xs, lambda_s, vs);

  p0 = 1.0/Gamma;

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(Gas));

        switch(iprob) {
          case 1: /* PARKER WIND */
                  if (x1 < xs) { 
                    a = TINY_NUMBER;
                    b = vs;
                  } else { 
                    a = vs;
                    b = HUGE_NUMBER;
                  }
                  break;
          case 2: /* BONDI ACCRETION */
                  if (x1 < xs) { 
                    a = vs;
                    b = HUGE_NUMBER;
                  } else { 
                    a = TINY_NUMBER;
                    b = vs;
                  }
                  break;
          default:  ath_error("[cylwind]:  Not an accepted problem number!\n");
        }

        converged = bisection(myfunc,a,b,x1,&v);
        if (!converged) {
          ath_error("[cylwind]:  Bisection did not converge!\n");
        }

        pG->U[k][j][i].d   = lambda_s/(x1*v);

        switch(iprob) {
          case 1: /* PARKER WIND */
                  pG->U[k][j][i].M1  = lambda_s/x1;
                  break; 
          case 2: /* BONDI ACCRETION */
                  pG->U[k][j][i].M1  = -lambda_s/x1;
                  break;
          default:  ath_error("[cylwind]:  Not an accepted problem number!\n");
        }

#ifdef MHD
        pG->U[k][j][i].B1c = b0/x1;
        pG->B1i[k][j][i]   = b0/(x1-0.5*pG->dx1);
#endif /* MHD */

        /* INITIALIZE TOTAL ENERGY */
#ifndef ISOTHERMAL
        pgas = p0*pow(pG->U[k][j][i].d,Gamma);
        pG->U[k][j][i].E = pgas/Gamma_1 
#ifdef MHD
          + 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c))
#endif /* MHD */
          + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#endif /* ISOTHERMAL */

        /* SAVE SOLUTION */
        Soln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  StaticGravPot = grav_pot;  
  x1GravAcc = grav_acc;
  set_bvals_mhd_fun(left_x1,do_nothing_bc);
  set_bvals_mhd_fun(right_x1,do_nothing_bc);
//   set_bvals_fun(right_x1,grad_or_bc);

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
  compute_l1_error("CylWind", pGrid, pDomain, Soln, 0);
}

/*============================ PRIVATE FUNCTIONS ===============================
/*------------------------------------------------------------------------------
 *  FUNCTION myfunc
 *
 *  THIS FUNCITON IS USED TO CALCULATE VELOCITY v AS A FUNCTION OF POSITION x
 *  USING lambda_c, THE CRITICAL VALUE OF THE DIMENSIONLESS MASS WIND/ACCRETION
 *  RATE.  STANDARD BISECTION IS USED TO FIND THE ROOT(S).
 *----------------------------------------------------------------------------*/
Real myfunc(const Real x, const Real v)
{
  return Gamma_1*(1/x + 1/Gamma_1 - 0.5*SQR(v))*pow(v*x,Gamma_1)
    - pow(lambda_s,Gamma_1);
}
