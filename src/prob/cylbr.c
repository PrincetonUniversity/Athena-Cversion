#include "copyright.h"
/*==============================================================================
 * FILE: cylbr.c
 *
 * A simple magnetostatic test of force balance using a B-field with uniform
 * R-component.  
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

static Real br0, omega0, vz0, rho0, pgas0, a, x2min, x2max;
static int iprob;

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return 0.5*SQR(x1*omega0);
              break;
    case 2:   return 0.5*SQR(x1*omega0) - 0.5*SQR(br0/x1);
              break;
    case 3:   return 0.5*SQR(x1*omega0) - 0.5*SQR(br0)/SQR(x1);
              break;
    default:  return 0.0;
  }
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  switch (iprob) {
    case 1:   return x1*SQR(omega0);
              break;
    case 2:   return x1*SQR(omega0) + SQR(br0)/pow(x1,3);
              break;
    case 3:   return x1*SQR(omega0) + SQR(br0)/pow(x1,3);
              break;
    default:  return 0.0;
  }
}

static Gas ***Soln=NULL;
void cylbr_ix1(Grid *pG);
void cylbr_ox1(Grid *pG);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pG, Domain *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3;
  Real x1,x2,x3,x1i,x2i,x3i,r1,r2,phi1,phi2,tlim;
//   Real x1min,x1max,x2min,x2max,x3min,x3max;
  Real Eint,Ekin,Emag,pgas,vphi;
  Real y1,y2,y3;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylbr]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylbr]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylbr]: Only (R,phi) can be used in 2D!\n");
  }

//   x1min = par_getd("grid","x1min");
//   x1max = par_getd("grid","x1max");
  x2min = par_getd("grid","x2min");
  x2max = par_getd("grid","x2max");
//   x3min = par_getd("grid","x3min");
//   x3max = par_getd("grid","x3max");
  a = 2.0*PI/(x2max-x2min);

  omega0 = par_getd("problem", "omega0");
  vz0    = par_getd("problem", "vz0");
  br0    = par_getd("problem", "br0");
  rho0   = par_getd("problem", "rho0");
  pgas0  = par_getd("problem", "pgas0");
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
          case 1:   pG->U[k][j][i].d   = rho0;
                    pG->B1i[k][j][i]   = br0/r1;
                    pG->U[k][j][i].B1c = br0/x1;
                    pgas = 1.0;
                    break;
          case 2:   pG->U[k][j][i].d   = SQR(sin(a*x2));
                    pG->B1i[k][j][i]   = br0*cos(a*x2)/r1;
                    pG->U[k][j][i].B1c = br0*cos(a*x2)/x1;
                    pgas = 1.0 + 0.5*SQR(br0*sin(a*x2)/x1);
                    break;
         case 3:    vphi = omega0*y1;
                    pG->U[k][j][i].d   = rho0 + SQR(sin(a*(x2-x2min)));
                    pG->B1i[k][j][i]   = br0*(cos(a*(x2-x2min)))/r1;
                    pG->U[k][j][i].B1c = br0*(cos(a*(x2-x2min)))/x1;
                    pgas = pgas0 + 0.5*SQR(br0)*pG->U[k][j][i].d/SQR(x1);
                    break;
          default:  ath_error("[cylbr]:  Not an accepted problem number\n");
        }

        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*vphi;
//         pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega0;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz0;

        Eint = pgas/Gamma_1;
// printf("(%d,%d,%d) Pg = %f\n", i,j,k,pgas);
        Emag = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
        pG->U[k][j][i].E = Eint + Emag + Ekin;

        /* SAVE SOLUTION */
        Soln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  set_bvals_mhd_fun(left_x1,cylbr_ix1);
  set_bvals_mhd_fun(right_x1,cylbr_ox1);

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

void Userwork_in_loop(Grid *pG, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pG, Domain *pDomain)
{
  compute_l1_error("CylBR", pG, pDomain, Soln, 0);
}



/*============================================================================
 * BOUNDARY CONDITION FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* B_R = B_0/R boundary conditions, Inner x1 boundary
 */

void cylbr_ix1(Grid *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3,y1,y2,y3,r1,r2;
  Real phi0,Eint,Emag,Ekin,vphi,pgas,C;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      cc_pos(pG,is,j,k,&x1,&x2,&x3);
      phi0 = a*(x2 - omega0*pG->time - x2min);
      C = (x1-0.5*pG->dx1)*pG->B1i[k][j][is];
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,is-i,j,k,&x1,&x2,&x3);
        vc_pos(pG,is-i,j,k,&y1,&y2,&y3);
        r1 = x1 - 0.5*pG->dx1;
        r2 = x1 + 0.5*pG->dx1;

        vphi = omega0*y1;
//         pG->U[k][j][is-i].d   = rho0 + SQR(sin(phi0));
        pG->U[k][j][is-i].d   = pG->U[k][j][is].d;
//         pG->B1i[k][j][is-i]   = br0*(cos(phi0))/r1;
//         pG->U[k][j][is-i].B1c = br0*(cos(phi0))/x1;
        pG->B1i[k][j][is-i]   = C/r1;
        pG->U[k][j][is-i].B1c = C/x1;
        pgas = pgas0 + 0.5*SQR(br0)*pG->U[k][j][is-i].d/SQR(x1);

        pG->U[k][j][is-i].M1 = 0.0;
        pG->U[k][j][is-i].M2 = pG->U[k][j][is-i].d*vphi;
        pG->U[k][j][is-i].M3 = pG->U[k][j][is-i].d*vz0;

        Eint = pgas/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][is-i].B1c) + SQR(pG->U[k][j][is-i].B2c) + SQR(pG->U[k][j][is-i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][is-i].M1) + SQR(pG->U[k][j][is-i].M2) + SQR(pG->U[k][j][is-i].M3))/pG->U[k][j][is-i].d; 
        pG->U[k][j][is-i].E = Eint + Emag + Ekin;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* B_R = B_0/R boundary conditions, Outer x1 boundary
 */

void cylbr_ox1(Grid *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3,y1,y2,y3,r1,r2;
  Real phi0,Eint,Emag,Ekin,vphi,pgas,C;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      cc_pos(pG,ie,j,k,&x1,&x2,&x3);
      phi0 = a*(x2 - omega0*pG->time - x2min);
      C = x1*pG->U[k][j][ie].B1c;
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
        vc_pos(pG,ie+i,j,k,&y1,&y2,&y3);
        r1 = x1 - 0.5*pG->dx1;
        r2 = x1 + 0.5*pG->dx1;

        vphi = omega0*y1;
//         pG->U[k][j][ie+i].d   = rho0 + SQR(sin(phi0));
        pG->U[k][j][ie+i].d   = pG->U[k][j][ie].d;
//         pG->B1i[k][j][ie+i]   = br0*(cos(phi0))/r1;
//         pG->U[k][j][ie+i].B1c = br0*(cos(phi0))/x1;
        pG->B1i[k][j][ie+i]   = C/r1;
        pG->U[k][j][ie+i].B1c = C/x1;
        pgas = pgas0 + 0.5*SQR(br0)*pG->U[k][j][ie+i].d/SQR(x1);

        pG->U[k][j][ie+i].M1 = 0.0;
        pG->U[k][j][ie+i].M2 = pG->U[k][j][ie+i].d*vphi;
        pG->U[k][j][ie+i].M3 = pG->U[k][j][ie+i].d*vz0;

        Eint = pgas/Gamma_1;
        Emag = 0.5*(SQR(pG->U[k][j][ie+i].B1c) + SQR(pG->U[k][j][ie+i].B2c) + SQR(pG->U[k][j][ie+i].B3c));
        Ekin = 0.5*(SQR(pG->U[k][j][ie+i].M1) + SQR(pG->U[k][j][ie+i].M2) + SQR(pG->U[k][j][ie+i].M3))/pG->U[k][j][ie+i].d; 
        pG->U[k][j][ie+i].E = Eint + Emag + Ekin;
      }
    }
  }

  return;
}
