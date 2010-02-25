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

static Real bphi0, omega0, vz0, rho0, pgas0;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return 0.5*SQR(x1*omega0) - (SQR(bphi0)/rho0)*log(x1);
              break;
    case 2:   return 0.5*SQR(x1*omega0);
              break;
    default:  return 0.0;
  }
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return x1*SQR(omega0) - SQR(bphi0)/(rho0*x1);
              break;
    case 2:   return x1*SQR(omega0);
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
  Real Eint, Emag, Ekin;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef MHD
  ath_error("[cylbphi]: This problem only works in MHD!\n");
#endif
#ifndef CYLINDRICAL
  ath_error("[cylbphi]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylbphi]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylbphi]: Only (R,phi) can be used in 2D!\n");
  }

  omega0 = par_getd("problem", "omega0");
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
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega0;
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

        /* SAVE SOLUTION */
        Soln[k][j][i] = pG->U[ks][j][i];
      }
    }
  }

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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  compute_l1_error("CylBPhi", pGrid, pDomain, Soln, 0);
}

