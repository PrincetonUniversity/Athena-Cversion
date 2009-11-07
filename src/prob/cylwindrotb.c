#include "copyright.h"
/*==============================================================================
 * FILE: cylwindrotb.c
 *
 *  A test of the magnetized rotational wind problem.  
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

static Real theta, omega, eta, E, v_z;

const Real rho_A = 1.0;
const Real R_A   = 1.0;
const Real GM    = 1.0;


static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return -GM/x1;
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return GM/SQR(x1);
}

Real myfunc(const Real x, const Real y);
static Gas ***Soln=NULL;


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k,n,converged;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1, nx2, nx3;
  Real x1, x2, x3;
  Real a,b,c,d,xmin,xmax,ymin,ymax;
  Real x,y,xslow,yslow,xfast,yfast;
  Real R0,R1,R2,rho,Mdot,K,Omega,P_gas,beta,v_R,B_R,v_phi,B_phi;
  Gas *Wind=NULL;
  Real *pU,*pWl,*pWr;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylwindrotb]: This problem only works in cylindrical!\n");
#endif
#ifndef MHD
  ath_error("[cylwindrotb]: This problem only works in MHD!\n");
#endif

  if (nx1==1) {
    ath_error("[cylwind]: Only R can be used in 1D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylwind]: Only (R,phi) can be used in 2D!\n");
  }

  // ALLOCATE MEMORY FOR WIND SOLUTION 
  if ((Wind = (Gas*)calloc_1d_array(nx1,sizeof(Gas))) == NULL)
    ath_error("[cylwindrotb]: Error allocating memory\n");

  // ALLOCATE MEMORY FOR GRID SOLUTION 
  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylwindrotb]: Error allocating memory\n");

  theta = par_getd("problem","theta");
  omega = par_getd("problem","omega");
  v_z   = par_getd("problem","v_z");

  // NUMERICAL SOLUTION OBTAINED FROM MATLAB
  // HOPEFULLY WE CAN REPLACE THIS WITH A FUNCTION SOLVER LATER...
  xslow = 0.5243264128;
  yslow = 2.4985859152;
  xfast = 1.6383327831;
  yfast = 0.5373957134;
  E     = 7.8744739104;
  eta   = 2.3608500383;

  xmin = par_getd("grid","x1min")/R_A;
  xmax = par_getd("grid","x1max")/R_A;
  ymin = 0.45/rho_A;
  ymax = 2.6/rho_A;

  printf("theta = %f,\t omega = %f,\t eta = %f,\t E = %f\n", theta,omega,eta,E);
  printf("xslow = %f,\t yslow = %f,\t xfast = %f,\t yfast = %f\n", xslow,yslow,xfast,yfast);
  printf("xmin = %f,\t ymin = %f,\t xmax = %f,\t ymax = %f\n", xmin,ymin,xmax,ymax);


  // CALCULATE THE WIND SOLUTION AT CELL-CENTERS
  for (i=il; i<=iu; i++) {
    cc_pos(pG,i,jl,kl,&x1,&x2,&x3);
    R0 = x1;
    x = R0/R_A;

    if (x < xslow) {
      sign_change(myfunc,yslow,10*ymax,x,&a,&b);
      sign_change(myfunc,b,10*ymax,x,&a,&b);
    } else if (x < 1.0) {
      sign_change(myfunc,1.0+TINY_NUMBER,yslow,x,&a,&b);
    } else if (x < xfast) {
      sign_change(myfunc,yfast,1.0-TINY_NUMBER,x,&a,&b);
      if (!sign_change(myfunc,b,1.0-TINY_NUMBER,x,&a,&b)) {
        a = yfast;
        b = 1.0-TINY_NUMBER;
      }
    } else {
      sign_change(myfunc,0.5*ymin,yfast,x,&a,&b);
    }
    converged = bisection(myfunc,a,b,x,&y);
    if(!converged) {
      ath_error("[cylwindrotb]:  Bisection did not converge!\n");
    }

    rho = rho_A*y;
    Mdot = sqrt(R_A*SQR(rho_A)*GM*eta);
    Omega = sqrt((GM*omega)/pow(R_A,3));
    K = (GM*theta)/(Gamma*pow(rho_A,Gamma_1)*R_A);
    P_gas = K*pow(rho,Gamma);
    v_R = Mdot/(R0*rho);
    beta = sqrt(1.0/rho_A);
    B_R = beta*rho*v_R;
    v_phi = R0*Omega*(1.0/SQR(x)-y)/(1.0-y);
    B_phi = beta*rho*(v_phi-R0*Omega);

// printf("\n");
// printf("r = %f,\t phi = %f\n", y1,x2);
// printf("rho = %f\n", rho);
// printf("Phi = %f\n", Phi);
// printf("K = %f\n", K);
// printf("Omega = %f\n", Omega);
// printf("f = %f\n", f);
// printf("P = %f\n", P);
// printf("vr = %f\n", vr);
// printf("vphi = %f\n", vphi);
// printf("Br = %f\n", Br);
// printf("Bphi = %f\n", Bphi);
// printf("\n");

    Wind[i].d   = rho;
    Wind[i].M1  = rho*v_R;
    Wind[i].M2  = rho*v_phi;
    Wind[i].M3  = rho*v_z;

    Wind[i].B1c = B_R;
    Wind[i].B2c = B_phi;
    Wind[i].B3c = 0.0;

    Wind[i].E   = P_gas/Gamma_1 
      + 0.5*(SQR(Wind[i].B1c) + SQR(Wind[i].B2c) + SQR(Wind[i].B3c))
      + 0.5*(SQR(Wind[i].M1) + SQR(Wind[i].M2) + SQR(Wind[i].M3))/Wind[i].d;
  }

  // COPY THE WIND SOLUTION ACROSS THE GRID, SAVE FOR ERROR ANALYSIS
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i] = Wind[i];

        // SAVE SOLUTION
        Soln[k][j][i]  = Wind[i];
      }
    }
  }


// /* INITIALIZE INTERFACE B-FIELDS */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      pG->B1i[k][j][il] = Wind[il].B1c;
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,jl,kl,&x1,&x2,&x3);
        R1 = x1 - pG->dx1;
        R0 = x1 - 0.5*pG->dx1;
        R2 = x1;
        pG->B1i[k][j][i] = (R1*Wind[i-1].B1c + R2*Wind[i].B1c)/(2.0*R0);
        pG->B2i[k][j][i] = Wind[i].B2c;
        pG->B3i[k][j][i] = 0.0;
      }
    }
  }


  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_mhd_fun(left_x1,do_nothing_bc);
  set_bvals_mhd_fun(right_x1,do_nothing_bc);

  free_1d_array((void *)Wind);

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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
//   printf("Max divB = %1.10e\n", compute_div_b(pGrid));
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  compute_l1_error("CylWindRotB", pGrid, pDomain, Soln, 0);
}


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*-----------------------------------------------------------------------------
 * Function func
 *
 * This function is used to calculate y (ie. rho) as a function of x (ie. R),
 * gamma, eta, theta, omega, and E using the bisection method.
 *
 */

Real myfunc(const Real x, const Real y) 
{
  return eta/(2.0*pow(x,2)*pow(y,2)) + (theta/Gamma_1)*pow(y,Gamma_1) 
    - 1.0/x + 0.5*omega*(pow(x-1.0/x,2)/pow(y-1,2) - pow(x,2)) - E;
}
