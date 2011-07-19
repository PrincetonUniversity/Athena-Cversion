#include "copyright.h"
/*============================================================================*/
/*! \file cyladvect.c
 *  \brief A simple density-pulse advection test in cylindrical coordinates
 *
 * PURPOSE:  A simple density-pulse advection test in cylindrical coordinates
 *   with no pressure or tension forces.
 *
 * AUTHOR:  Aaron Skinner
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * grav_pot() - gravitational potential
 * cyladvect_ix1() - Inner-R boundary conditions for advection problem
 * cyladvect_ox1() - Outer-R boundary conditions for advection problem
 *============================================================================*/

static Real omega0,rho0,bz0,vz0,Pgas0,amp,R0,phi0,z0,rad,alpha,x2min,x2max;
static int iprob;
static Real grav_pot(const Real x1, const Real x2, const Real x3);
void cyladvect_ix1(GridS *pG);
void cyladvect_ox1(GridS *pG);
static ConsS ***RootSoln=NULL;
#ifdef FARGO
static Real Omega(const Real R) {return omega0;}
static Real Shear(const Real R) {return 0.0;}
#endif /* Fargo */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3; 
  Real x1,x2,x3;
  Real Eint,Emag,Ekin;
  int mask=0;
  Real x,y,z,x0,y0,r;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cyladvect]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cyladvect]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cyladvect]: Only (R,phi) can be used in 2D!\n");
  }

  /* allocate memory for solution */
  if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS))) == NULL)
    ath_error("[cyladvect]: Error allocating memory for solution\n");

  /* parse input file */
  x2min  = par_getd("domain1", "x2min");
  x2max  = par_getd("domain1", "x2max");
  omega0 = par_getd("problem", "omega0");
  rho0   = par_getd("problem", "rho0");
  bz0    = par_getd("problem", "bz0");
  vz0    = par_getd("problem", "vz0");
  Pgas0  = par_getd("problem", "Pgas0");
  amp    = par_getd("problem", "amp");
  R0     = par_getd("problem", "R0");
  phi0   = par_getd("problem", "phi0");
  z0     = par_getd("problem", "z0");
  rad    = par_getd("problem", "rad");
  iprob  = par_geti("problem", "iprob");

  /* Period is one grid-crossing */
  alpha = 2.0*PI/(x2max-x2min);

  /* Set density and phi-momentum */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));

        switch (iprob) {
          case 1:
                  mask = (int)((x1>=R0-rad) && (x1<=R0+rad) && (x2>=phi0-rad/R0) && (x2<=phi0+rad/R0));
                  pG->U[k][j][i].d  = rho0*(1.0 + amp*mask);
#ifndef FARGO
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega0;
#else
                  pG->U[k][j][i].M2 = 0;
#endif
#ifdef MHD
                  pG->U[k][j][i].B3c = bz0*mask;
                  pG->B3i[k][j][i]   = bz0*mask;
#endif
                  break;
          case 2:
                  pG->U[k][j][i].d  = rho0*(1.0 + amp*sin(alpha*x2-phi0));
#ifndef FARGO
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega0;
#else
                  pG->U[k][j][i].M2 = 0;
#endif
                  break;
          case 3:
                  x  = x1*cos(x2);    y  = x1*sin(x2);    z  = x3;
                  x0 = R0*cos(phi0);  y0 = R0*sin(phi0);
                  r = sqrt(SQR(x-x0) + SQR(y-y0) + SQR(z-z0));
                  pG->U[k][j][i].d  = rho0*(1.0 + exp(-0.5*SQR(3.0*r/rad)));
#ifndef FARGO
                  pG->U[k][j][i].M2 = pG->U[k][j][i].d*x1*omega0;
#else
                  pG->U[k][j][i].M2 = 0;
#endif
#ifdef MHD
                  pG->U[k][j][i].B3c = bz0;
                  pG->B3i[k][j][i]   = bz0;
#endif
                  break;
          default:
                  ath_error("[cyladvect]:  Not an accepted problem number\n");
        }

        pG->U[k][j][i].M3 = pG->U[k][j][i].d*vz0;

#ifndef ISOTHERMAL
        /* Initialize total energy */
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*SQR(pG->U[k][j][i].B3c);
#endif
        Eint = (Pgas0-Emag)/Gamma_1;
        Ekin = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;

        pG->U[k][j][i].E = Eint + Emag + Ekin;
#endif /* ISOTHERMAL */

        /* Save solution */
        RootSoln[k][j][i] = pG->U[k][j][i];
      }
    }
  }

  StaticGravPot = grav_pot;
  if (iprob==2) {
    bvals_mhd_fun(pDomain, left_x1,  cyladvect_ix1);
    bvals_mhd_fun(pDomain, right_x1, cyladvect_ox1);
  } else {
    bvals_mhd_fun(pDomain, left_x1,  do_nothing_bc);
    bvals_mhd_fun(pDomain, right_x1, do_nothing_bc);
  }
#ifdef FARGO
  OrbitalProfile = Omega;
  ShearProfile = Shear;
#endif

  return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
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

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
#ifdef MHD
  printf("Max divB = %1.10e\n", compute_div_b(pM->Domain[0][0].Grid));
#endif
}

void Userwork_after_loop(MeshS *pM)
{
  compute_l1_error("CylAdvect", pM, RootSoln, 0);
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  return 0.5*SQR(x1*omega0);
}

/*----------------------------------------------------------------------------*/
/*! \fn void cyladvect_ix1(GridS *pG)
 *  \brief Inner-R boundary conditions for advection problem
 */

void cyladvect_ix1(GridS *pG)
{
  int is = pG->is;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real Eint,Emag,Ekin;

  phi0 = alpha*(omega0*pG->time + x2min);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,is-i,j,k,&x1,&x2,&x3);
        pG->U[k][j][is-i].d   = rho0*(1.0 + amp*sin(alpha*x2-phi0));
#ifndef FARGO
        pG->U[k][j][is-i].M2  = pG->U[k][j][is-i].d*x1*omega0;
#else
        pG->U[k][j][is-i].M2  = 0.0;
#endif
        pG->U[k][j][is-i].M3  = pG->U[k][j][is-i].d*vz0;
#ifdef MHD
        pG->U[k][j][is-i].B3c = bz0;
        pG->B3i[k][j][is-i]   = bz0;
#endif

#ifndef ISOTHERMAL
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*SQR(pG->U[k][j][is-i].B3c);
#endif
        Eint = (Pgas0-Emag)/Gamma_1;
        Ekin = 0.5*(SQR(pG->U[k][j][is-i].M1) + SQR(pG->U[k][j][is-i].M2) + SQR(pG->U[k][j][is-i].M3))/pG->U[k][j][is-i].d;

        pG->U[k][j][is-i].E = Eint + Emag + Ekin;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void cyladvect_ox1(GridS *pG)
 *  \brief Outer-R boundary conditions for advection problem
 */

void cyladvect_ox1(GridS *pG)
{
  int ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3;
  Real Eint,Emag,Ekin;

  phi0 = alpha*(omega0*pG->time + x2min);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
        pG->U[k][j][ie+i].d   = rho0*(1.0 + amp*sin(alpha*x2-phi0));
#ifndef FARGO
        pG->U[k][j][ie+i].M2  = pG->U[k][j][ie+i].d*x1*omega0;
#else
        pG->U[k][j][ie+i].M2  = 0.0;
#endif
        pG->U[k][j][ie+i].M3  = pG->U[k][j][ie+i].d*vz0;
#ifdef MHD
        pG->U[k][j][ie+i].B3c = bz0;
        pG->B3i[k][j][ie+i]   = bz0;
#endif

#ifndef ISOTHERMAL
        Emag = 0.0;
#ifdef MHD
        Emag = 0.5*SQR(pG->U[k][j][ie+i].B3c);
#endif
        Eint = (Pgas0-Emag)/Gamma_1;
        Ekin = 0.5*(SQR(pG->U[k][j][ie+i].M1) + SQR(pG->U[k][j][ie+i].M2) + SQR(pG->U[k][j][ie+i].M3))/pG->U[k][j][ie+i].d;

        pG->U[k][j][ie+i].E = Eint + Emag + Ekin;
#endif /* ISOTHERMAL */
      }
    }
  }

  return;
}
