#include "copyright.h"
/*============================================================================*/
/*! \file cylspiral.c
 *  \brief A test of Global Spiral Shocks in a Galactic Disk.              
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

static Real rho0,vcir,vt;
static Real atime;
static Real n_arm,i_pitch,omega_p,F_arm,ephi,mtani;
static Real taper_i,taper_o,t_arm;
static Real grav_pot(const Real x1, const Real x2, const Real x3);

/* take a simple analytic form for the rotation curve */
static Real Omega(const Real R) {
  vt = vcir*(1.-1./(1+R*R));
  return vt/R;
}
static Real Shear(const Real R) {
 return 1.-2./(1.+R*R);
}

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int myID_Comm_world = 0;
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3;
  Real x1,x2,x3;
  Real Eint,Emag,Ekin;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
  ath_error("[cylrayleigh]: This problem only works in cylindrical!\n");
#endif

  if (nx1==1) {
    ath_error("[cylrayleigh]: This problem can only be run in 2D or 3D!\n");
  }
  else if (nx2==1 && nx3>1) {
    ath_error("[cylrayleigh]: Only (R,phi) can be used in 2D!\n");
  }

#ifdef MPI_PARALLEL
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD, &myID_Comm_world))
    ath_error("[cylrayleigh]: Error on calling MPI_Comm_rank\n");
#endif

#ifdef MHD
  bphi0       = par_getd("problem", "bphi0");
#endif
  rho0        = par_getd("problem", "rho0");
  vcir        = par_getd("problem", "vcir");

  /* read the spiral arm parmeters*/
  
  n_arm       = par_getd("problem","n_arm");
  i_pitch     = par_getd("problem","i_pitch");
  F_arm       = par_getd("problem", "F_arm");
  taper_i     = par_getd("problem", "taper_i");
  taper_o     = par_getd("problem", "taper_o");
  t_arm       = par_getd("problem", "t_arm");
  omega_p     = par_getd("problem", "omega_p");

  atime = 0.0;
  mtani = n_arm/tan(i_pitch*PI/180.);
  ephi  = F_arm*SQR(vcir)/mtani; 

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);

//        vt = vcir*(1.-1./(1+x1*x1));
        vt = x1*Omega(x1);

        pG->U[k][j][i].d  = rho0;
        pG->U[k][j][i].M1 = 0.0;

#ifndef FARGO
        pG->U[k][j][i].M2 = rho0*vt;
#else
        pG->U[k][j][i].M2 = 0.0;
#endif

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
  bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
  bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
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
#ifdef MHD
  bphi0       = par_getd("problem", "bphi0");
#endif
  rho0        = par_getd("problem", "rho0");

  StaticGravPot = grav_pot;
//   bvals_mhd_fun(pDomain,left_x1,do_nothing_bc);
//   bvals_mhd_fun(pDomain,right_x1,do_nothing_bc);
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
 atime = pM->time;
}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  Real arg,arm_amp,Phi_bg,Phi_sp;
  arg = x1*x1 + 1.;
  Phi_bg = SQR(vcir)*0.5*(1./arg + log(arg));

  arm_amp = ephi;
/* tapering near the inner and outer boundaries*/
  if(x1 < taper_i) arm_amp =ephi*exp(-5.*SQR(x1-taper_i));
  if(x1 > taper_o) arm_amp =ephi*exp(-5.*SQR(x1-taper_o));

/* slow increase of the arm strength */
  if(atime < t_arm) arm_amp=arm_amp*atime/t_arm;
  Phi_sp  = arm_amp*cos(n_arm*(x2-omega_p*atime) + mtani*log(x1));

  return  Phi_bg+Phi_sp;
}

