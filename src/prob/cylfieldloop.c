#include "copyright.h"
/*==============================================================================
 * FILE: cylfieldloop.c
 *
 * PURPOSE: Problem generator for advection of a field loop test in cylindrical
 *   coordinates.  Can only be run in 2D or 3D.  Input parameters are:
 *      problem/r0    = radial coordinate of loop center
 *      problem/phi0  = angular coordinate of loop center
 *      problem/rad   = radius of field loop
 *      problem/amp   = amplitude of vector potential (and therefore B)
 *      problem/omega = flow angular velocity
 *      problem/vz0   = flow vertical velocity
 *
 * REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD
 *   via constrined transport", JCP, 205, 509 (2005)
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


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * grav_pot() - gravitational potential
 * grav_acc() - gravitational acceleration
 *============================================================================*/

static Real grav_pot(const Real x1, const Real x2, const Real x3);
static Real grav_acc(const Real x1, const Real x2, const Real x3);

static Real omega;
static ConsS ***RootSoln=NULL;

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real X0,Y0,X,Y,R1,R2,dist;
  Real x1,x2,x3,x1i,x2i,x3i,y1,y2,y3;
  Real **a3;
  Real r0,phi0,amp,rad,vz0;

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  il = is-nghost;  iu = ie+nghost;
  jl = js-nghost;  ju = je+nghost;
  kl = ks-nghost*(ke>ks);  ku = ke+nghost*(ke>ks);

  nx1 = iu-il+1;
  nx2 = ju-jl+1;
  nx3 = ku-kl+1;

#ifndef MHD
  ath_error("[cylfieldloop]: This problem can only be run in MHD!\n");
#endif

  if ((ie-is)==0 || (je-js)==0) {
    ath_error("[cylfieldloop]: This problem can only be run in 2D or 3D\n");
  }

  if ((a3 = (Real**)calloc_2d_array(nx2+1, nx1+1, sizeof(Real))) == NULL) {
    ath_error("[cylfieldloop]: Error allocating memory for vector pot\n");
  }

  /* ALLOCATE MEMORY FOR SOLUTION */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cylfieldloop]: Error allocating memory for solution\n");

  /* READ INITIAL CONDITIONS */
  r0    = par_getd("problem","r0");
  phi0  = par_getd("problem","phi0");
  amp   = par_getd("problem","amp");
  rad   = par_getd("problem","rad");
  omega = par_getd("problem","omega");
  vz0   = par_getd("problem","vz0");

  /* (X0,Y0) IS THE CENTER OF THE LOOP */
  X0 = r0*cos(phi0);
  Y0 = r0*sin(phi0);

  /* INITIALIZE VECTOR POTENTIAL FOR FIELD LOOP */
  for (j=jl; j<=ju+1; j++) {
    for (i=il; i<=iu+1; i++) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);
      x1 -= 0.5*pG->dx1;
      x2 -= 0.5*pG->dx2;

      /* (X,Y) IS GRID CELL CORNER */
      X = x1*cos(x2);
      Y = x1*sin(x2);

      /* dist IS DISTANCE FROM CENTER OF LOOP TO GRID CELL CORNER */
      dist = sqrt(SQR(X-X0) + SQR(Y-Y0));
      if (dist < rad) {
        a3[j][i] = amp*(rad - dist);
      }
    }
  }

  /* INITIALIZE DENSITY, MOMENTA, INTERFACE FIELDS */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        vc_pos(pG,i,j,k,&y1,&y2,&y3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));
        R1 = x1 - 0.5*pG->dx1;

        pG->U[k][j][i].d = 1.0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*y1*omega;
        pG->U[k][j][i].M3 = vz0*pG->U[k][j][i].d;
        pG->B1i[k][j][i] = (a3[j+1][i] - a3[j][i])/(R1*pG->dx2);
        pG->B2i[k][j][i] = -(a3[j][i+1] - a3[j][i])/(pG->dx1);
        pG->B3i[k][j][i] = 0.0;
      }
    }
  }

  /* INITIALIZE TOTAL ENERGY AND CELL-CENTERED FIELDS */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        R1 = x1 - 0.5*pG->dx1;
        R2 = x1 + 0.5*pG->dx1;

        if (i < iu) 
          pG->U[k][j][i].B1c = 0.5*(R1*pG->B1i[k][j][i] + R2*pG->B1i[k][j][i+1])/x1;
        else
          pG->U[k][j][i].B1c = 0.5*(R1*pG->B1i[k][j][i] + (a3[j+1][i+1] - a3[j][i+1])/(x1*pG->dx2));

        if (j < ju) 
          pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]);
        else
          pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i] - (a3[j+1][i+1] - a3[j+1][i])/pG->dx1);

        if (ke > ks)
          if (k < ku)
            pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]);
          else
            pG->U[k][j][i].B3c = 0.0;
        else
          pG->U[k][j][i].B3c = pG->B3i[k][j][i];

#ifndef ISOTHERMAL
        pG->U[k][j][i].E = 1.0/Gamma_1
          + 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c))
          + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#endif /* ISOTHERMAL */

        // SAVE SOLUTION
        RootSoln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  free_2d_array((void**)a3);

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain, left_x1, do_nothing_bc);
  bvals_mhd_fun(pDomain, right_x1, do_nothing_bc);

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * Emag()                  - computes magnetic pressure (Bx^2 + By^2)
 * divB()                  - computes divergence of magnetic field
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

#ifdef MHD
static Real Emag(const GridS *pG, const int i, const int j, const int k)
{
  return (SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
}

static Real divB(const GridS *pG, const int i, const int j, const int k)
{
  Real qa,x1,x2,x3;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  qa = ((x1+0.5*pG->dx1)*pG->B1i[k][j][i+1] - (x1-0.5*pG->dx1)*pG->B1i[k][j][i])/(x1*pG->dx1)
        + (pG->B2i[k][j+1][i] - pG->B2i[k][j][i])/(x1*pG->dx2);
  if (pG->Nx[2] > 1)
    qa += (pG->B3i[k+1][j][i] - pG->B3i[k][j][i])/(pG->dx3);

  return qa;
}
#endif

ConsFun_t get_usr_expr(const char *expr)
{
#ifdef MHD
  if(strcmp(expr,"Emag")==0) return Emag;
  else if(strcmp(expr,"divB")==0) return divB;
#endif
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylFieldLoop", pM, RootSoln);
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega);
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  return x1*SQR(omega);
}
