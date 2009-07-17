#include "copyright.h"
/*==============================================================================
 * FILE: cylwindrot.c
 *
 *  A test of rotational Bondi accretion (Parker wind) problem.  
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

static Real br0, ang_mom, c_infty, lambda_s, vz0;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -SQR(c_infty)/x1;
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return SQR(c_infty/x1);
}

Real myfunc(Real x, Real v);

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
  Real xs,vs,v,pgas,alpha,beta,a,b,converged;

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


  if (((ie-is) == 0) || ((je-js) == 0)) {
    fprintf(stderr,"[cylwindrot]: This problem can only be run in 2D or 3D\n");
    exit(EXIT_FAILURE);
  }

  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylwindrot]: Error allocating memory for solution\n");

  ang_mom = par_getd("problem","ang_mom");
  c_infty = par_getd("problem","c_infty");
  br0      = par_getd("problem","br0");
  vz0      = par_getd("problem","vz0");
  printf("gamma = %f,\t ang_mom = %f,\t c_infty = %f\n", Gamma, ang_mom, c_infty);

  beta = 2.0*Gamma_1/(Gamma+1.0);
  xs = (3.0-Gamma+sqrt(SQR(Gamma-3.0)-SQR(4.0*ang_mom)))/4.0;
  lambda_s = 1.0/Gamma_1*pow(xs,beta)+pow(xs,beta-1.0)-0.5*SQR(ang_mom)*pow(xs,beta-2.0);
  lambda_s = pow(lambda_s/(0.5+1.0/Gamma_1),1.0/beta);
  vs = pow(lambda_s/xs,0.5*beta);
  printf("xs = %13.10f,\t lambda_s = %13.10f,\t vs = %13.10f\n", xs, lambda_s, vs);




  for (i=il; i<=iu; i++) {
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    memset(&(pG->U[ks][js][i]),0.0,sizeof(Gas));
    vs = pow(lambda_s/x1,0.5*beta);

    switch(iprob) {
      case 0: // PARKER WIND
              if (x1 < xs) { 
                a = TINY_NUMBER;
                b = vs;
              } else { 
                a = vs;
                b = HUGE_NUMBER;
              }
              break;
      case 1: // BONDI ACCRETION
              if (x1 < xs) { 
                a = vs;
                b = HUGE_NUMBER;
              } else { 
                a = TINY_NUMBER;
                b = vs;
              }
              break;
      default:  ath_error("[cylwindrot]:  Not an accepted problem number\n");
    }

    converged = bisection(myfunc,a,b,x1,&v);
    if (!converged) {
      ath_error("[cylwindrot]:  Bisection did not converge!\n");
    }

    pG->U[ks][js][i].d = lambda_s/(x1*v);

    switch(iprob) {
      case 0: // PARKER WIND
              pG->U[ks][js][i].M1  = lambda_s/x1;
              break; 
      case 1: // BONDI ACCRETION
              pG->U[ks][js][i].M1  = -lambda_s/x1;
              break;
      default:  ath_error("[cylwindrot]:  Not an accepted problem number\n");
    }

    pG->U[ks][js][i].M2  = pG->U[ks][js][i].d*ang_mom/x1; 
    pG->U[ks][js][i].M3  = pG->U[ks][js][i].d*vz0; 
#ifdef MHD
    pG->U[ks][js][i].B1c = br0/x1; 
    pG->B1i[ks][js][i]   = br0/(x1-0.5*pG->dx1); 
#endif /* MHD */

    /* INITIALIZE TOTAL ENERGY */
#ifndef ISOTHERMAL
    pgas = (1.0/Gamma)*pow(pG->U[ks][js][i].d,Gamma);
    pG->U[ks][js][i].E = pgas/Gamma_1
#ifdef MHD
      + 0.5*SQR(pG->U[ks][js][i].B1c)
#endif /* MHD */ 
      + 0.5*(SQR(pG->U[ks][js][i].M1) + SQR(pG->U[ks][js][i].M2) + SQR(pG->U[ks][js][i].M3))/pG->U[ks][js][i].d;
#endif /* ISOTHERMAL */
  }


  // COPY WIND SOLUTION AND SAVE
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
	  for (i=il; i<=iu; i++) {
	    pG->U[k][j][i] = pG->U[ks][js][i];
#ifdef MHD
        pG->B1i[k][j][i] = pG->B1i[ks][js][i];
#endif /* MHD */


        // SAVE SOLUTION
        Soln[k][j][i] = pG->U[ks][js][i];
      }
    }
  }


  StaticGravPot = grav_pot;  
  x1GravAcc = grav_acc;
  set_bvals_fun(left_x1,do_nothing_bc);
  set_bvals_fun(right_x1,do_nothing_bc);
//  set_bvals_fun(right_x1,grad_or_bc);

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
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
  return NULL;
}

void Userwork_in_loop(Grid *pG, Domain *pDomain)
{
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */
void Userwork_after_loop(Grid *pG, Domain *pDomain)
{
  compute_l1_error("CylWindRot", pG, pDomain, Soln, 0);
}



/*-----------------------------------------------------------------------------
 * Function func
 *
 * This function is used to calculate v as a function of x, gamma, and lambda 
 * using the bisection method.  
 *
 */
Real myfunc(Real x, Real v) 
{
  return Gamma_1*(1/x + 1/Gamma_1 - 0.5*(SQR(v)+SQR(ang_mom/x)))*pow(v*x,Gamma_1) - pow(lambda_s,Gamma_1);
}

