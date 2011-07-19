#define SEED 661979

#include "copyright.h"
/*============================================================================*/
/*! \file cylrayleigh.c
 *  \brief A test of the Rayleigh instability using omega(R) = omega_0/R^q.
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

static Real rho0,omega0,q,bphi0,pgas0;
static Real grav_pot(const Real x1, const Real x2, const Real x3);

static Real Omega(const Real R) {
  return omega0/pow(R,q);
}
static Real Shear(const Real R) {
  return q;
}

static Real vphi(const Real x1, const Real x2, const Real x3) {
  return x1*Omega(x1);
}

static void diode_outflow_ix1(GridS *pGrid);
static void diode_outflow_ox1(GridS *pGrid);

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
  Real randnum,noise,noise_level;
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

  /* Seed the random number generator */
  srand(SEED + myID_Comm_world);

  omega0      = par_getd("problem", "omega0");
#ifdef MHD
  bphi0       = par_getd("problem", "bphi0");
#endif
  rho0        = par_getd("problem", "rho0");
  pgas0       = par_getd("problem", "pgas0");
  q           = par_getd("problem", "q");
  noise_level = par_getd("problem", "noise_level");

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));

        // Random number between 0 and 1
        randnum = ((double) rand()/((double)RAND_MAX + 1.0));
        // Random number between +/- noise_level
        noise = noise_level*(2.0*randnum-1.0);

        pG->U[k][j][i].d  = rho0;
#ifndef FARGO
        pG->U[k][j][i].M2 = rho0*avg1d(vphi,pG,i,j,k);
        // Now perturb v_phi
         if ((i>=is) && (i<=ie)) {
           pG->U[k][j][i].M2 *= (1.0 + noise);
         }
#else
         pG->U[k][j][i].M2 = 0.0;
        // Now perturb v_phi
         if ((i>=is) && (i<=ie)) {
           pG->U[k][j][i].M2 = noise*rho0*avg1d(vphi,pG,i,j,k);
         }
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
//   bvals_mhd_fun(pDomain,left_x1,diode_outflow_ix1);
//   bvals_mhd_fun(pDomain,right_x1,diode_outflow_ox1);
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
  omega0      = par_getd("problem", "omega0");
#ifdef MHD
  bphi0       = par_getd("problem", "bphi0");
#endif
  rho0        = par_getd("problem", "rho0");
  pgas0       = par_getd("problem", "pgas0");
  q           = par_getd("problem", "q");

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
}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*! \fn static Real grav_pot(const Real x1, const Real x2, const Real x3) 
 *  \brief Gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  if (q == 1.0) {
    return SQR(omega0)*log(x1);
  }
  else {
    Real omega = omega0/pow(x1,q);
    return 0.5*SQR(x1*omega)/(1.0-q);
  }
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions (w/diode condition), Inner x1 boundary */
void diode_outflow_ix1(GridS *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  const Real *r=pGrid->r,*ri=pGrid->ri;
  Real x1,x2,x3;
  Real L,B1,B3;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      L = r[is]*pGrid->U[k][j][is].M2;
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
        // HOLD ANGULAR MOMENTUM CONSTANT
        pGrid->U[k][j][is-i].M2 = L/r[is-i];
        pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      B1 = ri[is]*pGrid->B1i[k][j][is];
      for (i=1; i<=nghost; i++) {
        // HOLD FLUX OF B_R CONSTANT
        pGrid->B1i[k][j][is-i] = B1/ri[is-i];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }


  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      B3 = r[is]*pGrid->B3i[k][j][is];
      for (i=1; i<=nghost; i++) {
        // HOLD FLUX OF B_z CONSTANT
        pGrid->B3i[k][j][is-i] = B3/r[is-i];
      }
    }
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions (w/diode condition), Outer x1 boundary */
void diode_outflow_ox1(GridS *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  const Real *r=pGrid->r,*ri=pGrid->ri;
  Real x1,x2,x3;
  Real L,B1,B3;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      L = r[ie]*pGrid->U[k][j][ie].M2;
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
        // HOLD ANGULAR MOMENTUM CONSTANT
        pGrid->U[k][j][ie+i].M2 = L/r[ie+i];
        pGrid->U[k][j][ie+i].M1 = MAX(pGrid->U[k][j][ie+i].M1,0.0);
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      B1 = ri[ie+1]*pGrid->B1i[k][j][ie+1];
      for (i=2; i<=nghost; i++) {
        // HOLD FLUX OF B_R CONSTANT
        pGrid->B1i[k][j][ie+i] = B1/ri[ie+i];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }
  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      B3 = r[ie]*pGrid->B3i[k][j][ie];
      for (i=1; i<=nghost; i++) {
        // HOLD FLUX OF B_z CONSTANT
        pGrid->B3i[k][j][ie+i] = B3/r[ie+i];
      }
    }
  } 
#endif /* MHD */

  return;
}
