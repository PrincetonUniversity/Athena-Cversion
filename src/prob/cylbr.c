#include "copyright.h"
/*==============================================================================
 * FILE: cylbr.c
 *
 * A simple magnetostatic test of force balance using a B-field with uniform
 * R-component.  
 *
 *
 *
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

static Real br, omega, vz, rho, a;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return 0.5*SQR(x1*omega);
              break;
    case 2:   return 0.5*SQR(x1*omega) - 0.5*SQR(br/x1);
              break;
    case 3:   return 0.5*SQR(x1*omega) + 2.0*log(x1);
              break;
    default:  return 0.0;
  }
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return x1*SQR(omega);
              break;
    case 2:   return x1*SQR(omega) + SQR(br)/pow(x1,3);
              break;
    case 3:   return x1*SQR(omega) + 2.0/x1;
              break;
    default:  return 0.0;
  }
}

void cylbr_ir_bc(Grid *pG, int var_flag);
void cylbr_or_bc(Grid *pG, int var_flag);

static Gas ***Soln=NULL;



/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3;
  Real x1,x2,x3,x1i,x2i,x3i,r1,r2,phi1,phi2,pgas,magE,kinE,tlim;
  Real x1min, x1max, x2min, x2max, x3min, x3max;
  Real divB=0.0, maxdivB=0.0;
  Real y1,y2,y3;

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
    ath_error("[field_loop]: This problem can only be run in 2D or 3D\n");
  }

  x1min = par_getd("grid","x1min");
  x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
  x3min = par_getd("grid","x3min");
  x3max = par_getd("grid","x3max");
  a = 2.0*PI/(x2max-x2min);

  omega  = par_getd("problem", "omega");
  vz     = par_getd("problem", "vz");
  br     = par_getd("problem", "br");
  rho    = par_getd("problem", "rho");
  iprob  = par_geti("problem", "iprob");


  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[cylbr]: Error allocating memory for solution\n");


  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);

        r1 = x1 - 0.5*pG->dx1;
        r2 = x1 + 0.5*pG->dx1;
        phi1 = x2 - 0.5*pG->dx2;
        phi2 = x2 + 0.5*pG->dx2;

        switch (iprob) {
          case 1:   pG->U[k][j][i].d   = rho;
                    pG->B1i[k][j][i]   = br/r1;
                    pG->U[k][j][i].B1c = br/x1;
                    pgas = 1.0;
                    break;
          case 2:   pG->U[k][j][i].d   = SQR(sin(a*x2));
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = 1.0 + 0.5*SQR(br*sin(a*x2)/x1);
                    break;
          case 3:   pG->U[k][j][i].d   = (0.5*SQR(br)/SQR(x1))*(SQR(sin(a*x2))+1.0);
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = pG->U[k][j][i].d;
                    break;
          default:  ath_error("[cylbr]:  Not an accepted problem number\n");
        }

        pG->U[k][j][i].M1 = 0.0;
//         pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz;

        magE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        kinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
        pG->U[k][j][i].E = pgas/Gamma_1 + magE + kinE;

      }
    }
  }



/* SAVE SOLUTION */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        Soln[k][j][i] = pG->U[k][j][i];
      }
    }
  }


  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
//   set_bvals_fun(left_x1,do_nothing_bc);
//   set_bvals_fun(right_x1,do_nothing_bc);

  set_bvals_fun(left_x1,cylbr_ir_bc);
  set_bvals_fun(right_x1,cylbr_or_bc);

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
//   printf("Max divB = %1.10e\n", compute_div_b(pG));
}

/*---------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * ASSUMING WAVE HAS PROPAGATED AN INTEGER NUMBER OF PERIODS
 * Must set parameters in input file appropriately so that this is true
 */

void Userwork_after_loop(Grid *pG, Domain *pDomain)
{
 compute_l1_error("CylBR", pG, pDomain, Soln, 0);
}



/*=========================== PRIVATE FUNCTIONS ==============================*/

/*-----------------------------------------------------------------------------
 * Function cylbr_or_bc
 *
 * Time-dependent boundary condition for outer-R boundary.  The initial 
 * boundary values are simply advected in the phi-direction by solid-body 
 * rotation.
 */
void cylbr_or_bc(Grid *pG, int var_flag)
{
  int i,j,k;
  int is,ie,js,je,ks,ke,il,iu,jl,ju,kl,ku;
  Real x1,x2,x3,r1,r2,phi1,phi2,pgas,magE,kinE;
  Real y1,y2,y3;

  if (var_flag == 1) return;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  il = is-nghost;  iu = ie+nghost;
  jl = js-nghost;  ju = je+nghost;
  if (ke-ks == 0) {
    kl = ks;  ku = ke; 
  } else {
    kl = ks-nghost;  ku = ke+nghost;
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=ie+1; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);

        r1 = x1 - 0.5*pG->dx1;
        r2 = x1 + 0.5*pG->dx1;
        phi1 = x2 - 0.5*pG->dx2;
        phi2 = x2 + 0.5*pG->dx2;
        x2 = x2 - omega*pG->time;

        switch (iprob) {
          case 1:   pG->U[k][j][i].d   = rho;
                    pG->B1i[k][j][i]   = br/r1;
                    pG->U[k][j][i].B1c = br/x1;
                    pgas = 1.0;
                    break;
          case 2:   pG->U[k][j][i].d   = SQR(sin(a*x2));
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = 1.0 + 0.5*SQR(br*sin(a*x2)/x1);
                    break;
          case 3:   pG->U[k][j][i].d   = (0.5*SQR(br)/SQR(x1))*(SQR(sin(a*x2))+1.0);
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = pG->U[k][j][i].d;
                    break;
          default:  printf("[cylbr_or_bc]:  Not an accepted problem number\n");
        }
        pG->U[k][j][i].M1 = 0.0;
//         pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz;

        magE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        kinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
        pG->U[k][j][i].E = pgas/Gamma_1 + magE + kinE;
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * Function cylbr_ir_bc
 *
 * Time-dependent boundary condition for inner-R boundary.  The initial 
 * boundary values are simply advected in the phi-direction by solid-body 
 * rotation.
 */
void cylbr_ir_bc(Grid *pG, int var_flag)
{
  int i,j,k;
  int is,ie,js,je,ks,ke,il,iu,jl,ju,kl,ku;
  Real x1,x2,x3,r1,r2,phi1,phi2,pgas,magE,kinE;
  Real y1,y2,y3;

  if (var_flag == 1) return;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  il = is-nghost;  iu = ie+nghost;
  jl = js-nghost;  ju = je+nghost;
  if (ke-ks == 0) {
    kl = ks;  ku = ke; 
  } else {
    kl = ks-nghost;  ku = ke+nghost;
  }

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=0; i<=is-1; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);

        r1 = x1 - 0.5*pG->dx1;
        r2 = x1 + 0.5*pG->dx1;
        phi1 = x2 - 0.5*pG->dx2;
        phi2 = x2 + 0.5*pG->dx2;
        x2 = x2 - omega*pG->time;

        switch (iprob) {
          case 1:   pG->U[k][j][i].d   = rho;
                    pG->B1i[k][j][i]   = br/r1;
                    pG->U[k][j][i].B1c = br/x1;
                    pgas = 1.0;
                    break;
          case 2:   pG->U[k][j][i].d   = SQR(sin(a*x2));
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = 1.0 + 0.5*SQR(br*sin(a*x2)/x1);
                    break;
          case 3:   pG->U[k][j][i].d   = (0.5*SQR(br)/SQR(x1))*(SQR(sin(a*x2))+1.0);
                    pG->B1i[k][j][i]   = br*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br*cos(a*x2)/x1;
                    pgas = pG->U[k][j][i].d;
                    break;
          default:  printf("[cylbr_ir_bc]:  Not an accepted problem number\n");
        }
        pG->U[k][j][i].M1 = 0.0;
//         pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz;

        magE = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        kinE = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
        pG->U[k][j][i].E = pgas/Gamma_1 + magE + kinE;
      }
    }
  }
}

