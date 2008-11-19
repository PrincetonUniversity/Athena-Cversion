#include "copyright.h"
/*==============================================================================
 * FILE: bubble.c
 *
 * PURPOSE: Problem generator for bubble in spherical isothermal atmosphere.
 *
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2() - random number generator from NR
 * grav_pot() - gravitational potential for 2D problem (accn in Y)
 * outflow_ix1() - sets BCs on L-x1 (left edge) of grid used in 2D
 * outflow_ox1() - sets BCs on R-x1 (right edge) of grid used in 2D
 * outflow_ox2() - sets BCs on R-x2 (right edge) of grid used in 2D
 *============================================================================*/

static double ran2(long int *idum);
static Real grav_pot(const Real x1, const Real x2, const Real x3);
static void outflow_ix1(Grid *pGrid);
static void outflow_ox1(Grid *pGrid);
static void outflow_ox2(Grid *pGrid);
static void outflow_ix3(Grid *pGrid);
static void outflow_ox3(Grid *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  long int iseed = -1;
  Real x1,x2,x3,r2,rad_bubble=0.25;;
#ifdef MHD
  Real b0;
#endif
  int Nx1, Nx2, Nx3;
  int ixs, jxs, kxs;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  Nx1 = par_geti("grid","Nx1");
  Nx2 = par_geti("grid","Nx2");
  Nx3 = par_geti("grid","Nx3");

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->is + pGrid->idisp;
  jxs = pGrid->js + pGrid->jdisp;
  kxs = pGrid->ks + pGrid->kdisp;
  iseed = -1 - (ixs + Nx1*(jxs + Nx2*kxs));

/* Read magnetic field strength */
#ifdef MHD
  b0 = par_getd("problem","b0");
#endif

/* Initialize atmosphere, and bubble centered at y=0.25 and radius=rad_bubble
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
	pGrid->U[k][j][i].d = pow((1.0 + r2),-0.75);
	pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        pGrid->U[k][j][i].E = pow((1.0 + r2),-0.75)/(Gamma*Gamma_1);

        r2 = x1*x1 + (x2-0.25)*(x2-0.25) + x3*x3;
        if (r2 < (rad_bubble*rad_bubble)) {
	  pGrid->U[k][j][i].d = 0.01;
	}
#ifdef MHD
	pGrid->B1i[k][j][i] = b0;
	pGrid->U[k][j][i].B1c = b0;
        pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif
      }
#ifdef MHD
    pGrid->B1i[k][j][ie+1] = b0;
#endif
    }
  }

/* Enroll gravitational potential, and special BC fns */

  StaticGravPot = grav_pot;
  set_bvals_mhd_fun(left_x1,  outflow_ix1);
  set_bvals_mhd_fun(right_x1, outflow_ox1);
  set_bvals_mhd_fun(right_x2, outflow_ox2);
  set_bvals_mhd_fun(left_x3,  outflow_ix3);
  set_bvals_mhd_fun(right_x3, outflow_ox3);

/* With viscosity and/or resistivity, read eta_Ohm and nu_V */

#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif
#ifdef BRAGINSKII
  nu_V = par_getd("problem","nu");
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

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll special boundary value functions,
 *    and initialize gravity on restarts
 */

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  StaticGravPot = grav_pot;
  set_bvals_mhd_fun(left_x1,  outflow_ix1);
  set_bvals_mhd_fun(right_x1, outflow_ox1);
  set_bvals_mhd_fun(right_x2, outflow_ox2);
  set_bvals_mhd_fun(left_x3,  outflow_ix3);
  set_bvals_mhd_fun(right_x3, outflow_ox3);

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
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX

/* Static gravitational potential */
static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
  Real r2;
  r2 = x1*x1 + x2*x2 + x3*x3;
  return (0.75/Gamma)*log(1.0+r2);
}

/*------------------------------------------------------------------------------
 * outflow_ix1: special outflow boundary function
 */

static void outflow_ix1(Grid *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real r2,x1,x2,x3,d0,p0;
#ifdef MHD
  int ibc_x1,ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,is,j,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        d0 = pGrid->U[k][j][is].d/pow((1.0 + r2),-0.75);

        p0 = pGrid->U[k][j][is].E - 0.5*(SQR(pGrid->U[k][j][is].M1)
          + SQR(pGrid->U[k][j][is].M2) + SQR(pGrid->U[k][j][is].M3))
             /pGrid->U[k][j][is].d;
        p0 /= pow((1.0 + r2),-0.75);

        cc_pos(pGrid,(is-i),j,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        pGrid->U[k][j][is-i].d = d0*pow((1.0 + r2),-0.75);
        pGrid->U[k][j][is-i].M1 = 0.0;
        pGrid->U[k][j][is-i].M2 = 0.0;
        pGrid->U[k][j][is-i].M3 = 0.0;
        pGrid->U[k][j][is-i].E = p0*pow((1.0 + r2),-0.75);
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i]   = pGrid->B1i[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].B1c = pGrid->U[k][j][is+(i-1)].B1c;
      }
      if (ibc_x1 == 1) pGrid->B1i[k][j][is] = 0.0;
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i]   = pGrid->B2i[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].B2c = pGrid->U[k][j][is+(i-1)].B2c;
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i]   = pGrid->B3i[k][j][is+(i-1)];
        pGrid->U[k][j][is-i].B3c = pGrid->U[k][j][is+(i-1)].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*------------------------------------------------------------------------------
 * outflow_ox1: special outflow boundary function
 */

static void outflow_ox1(Grid *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real r2,x1,x2,x3,d0,p0;
#ifdef MHD
  int obc_x1,ju,ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,ie,j,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        d0 = pGrid->U[k][j][ie].d/pow((1.0 + r2),-0.75);

        p0 = pGrid->U[k][j][ie].E - 0.5*(SQR(pGrid->U[k][j][ie].M1)
          + SQR(pGrid->U[k][j][ie].M2) + SQR(pGrid->U[k][j][ie].M3))
             /pGrid->U[k][j][ie].d;
        p0 /= pow((1.0 + r2),-0.75);

        cc_pos(pGrid,(ie+i),j,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        pGrid->U[k][j][ie+i].d = d0*pow((1.0 + r2),-0.75);
        pGrid->U[k][j][ie+i].M1 = 0.0;
        pGrid->U[k][j][ie+i].M2 = 0.0;
        pGrid->U[k][j][ie+i].M3 = 0.0;
        pGrid->U[k][j][ie+i].E = p0*pow((1.0 + r2),-0.75);
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not set for the interface field B1i, except obc_x1=1 */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      if (obc_x1 == 1 ) pGrid->B1i[k][j][ie+1] = 0.0;
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i]   = pGrid->B1i[k][j][ie-(i-2)];
        pGrid->U[k][j][ie+i].B1c = pGrid->U[k][j][ie-(i-2)].B1c;
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i]   = pGrid->B2i[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].B2c = pGrid->U[k][j][ie-(i-1)].B2c;
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i]   = pGrid->B3i[k][j][ie-(i-1)];
        pGrid->U[k][j][ie+i].B3c = pGrid->U[k][j][ie-(i-1)].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*------------------------------------------------------------------------------
 * outflow_ox2: special outflow boundary function
 */

static void outflow_ox2(Grid *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k,il,iu; /* i-lower/upper */
  Real r2,x1,x2,x3,d0,p0;
#ifdef MHD
  int obc_x2,ku; /* k-upper */
#endif

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,je,k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        d0 = pGrid->U[k][je][i].d/pow((1.0 + r2),-0.75);

        p0 = pGrid->U[k][je][i].E - 0.5*(SQR(pGrid->U[k][je][i].M1)
          + SQR(pGrid->U[k][je][i].M2) + SQR(pGrid->U[k][je][i].M3))
             /pGrid->U[k][je][i].d;
        p0 /= pow((1.0 + r2),-0.75);

        cc_pos(pGrid,i,(je+j),k,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        pGrid->U[k][je+j][i].d = d0*pow((1.0 + r2),-0.75);
        pGrid->U[k][je+j][i].M1 = 0.0;
        pGrid->U[k][je+j][i].M2 = 0.0;
        pGrid->U[k][je+j][i].M3 = 0.0;
        pGrid->U[k][je+j][i].E = p0*pow((1.0 + r2),-0.75);
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i]   = pGrid->B1i[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].B1c = pGrid->U[k][je-(j-1)][i].B1c;
      }
    }
  }

/* j=je+1 is not set for the interface field B2i, except obc_x2=1 */
  for (k=ks; k<=ke; k++) {
    if (obc_x2 == 1) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+1][i] = 0.0;
      }
    }
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i]   = pGrid->B2i[k][je-(j-2)][i];
        pGrid->U[k][je+j][i].B2c = pGrid->U[k][je-(j-2)][i].B2c;
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i]   = pGrid->B3i[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].B3c = pGrid->U[k][je-(j-1)][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*------------------------------------------------------------------------------
 * outflow_ix3: special outflow boundary function
 */

static void outflow_ix3(Grid *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real r2,x1,x2,x3,d0,p0;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        d0 = pGrid->U[ks][j][i].d/pow((1.0 + r2),-0.75);

        p0 = pGrid->U[ks][j][i].E - 0.5*(SQR(pGrid->U[ks][j][i].M1)
          + SQR(pGrid->U[ks][j][i].M2) + SQR(pGrid->U[ks][j][i].M3))
             /pGrid->U[ks][j][i].d;
        p0 /= pow((1.0 + r2),-0.75);

        cc_pos(pGrid,i,j,(ks-1),&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        pGrid->U[ks-k][j][i].d = d0*pow((1.0 + r2),-0.75);
        pGrid->U[ks-k][j][i].M1 = 0.0;
        pGrid->U[ks-k][j][i].M2 = 0.0;
        pGrid->U[ks-k][j][i].M3 = 0.0;
        pGrid->U[ks-k][j][i].E = p0*pow((1.0 + r2),-0.75);
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i]   = pGrid->B1i[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].B1c = pGrid->U[ks+(k-1)][j][i].B1c;
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i]   = pGrid->B2i[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].B2c = pGrid->U[ks+(k-1)][j][i].B2c;
      }
    }
  }

  if (ibc_x3 == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks][j][i] = 0.0;
      }
    }
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i]   = pGrid->B3i[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].B3c = pGrid->U[ks+(k-1)][j][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}

/*------------------------------------------------------------------------------
 * outflow_ox3: special outflow boundary function
 */

static void outflow_ox3(Grid *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */
  Real r2,x1,x2,x3,d0,p0;

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pGrid,i,j,ke,&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        d0 = pGrid->U[ke][j][i].d/pow((1.0 + r2),-0.75);

        p0 = pGrid->U[ke][j][i].E - 0.5*(SQR(pGrid->U[ke][j][i].M1)
          + SQR(pGrid->U[ke][j][i].M2) + SQR(pGrid->U[ke][j][i].M3))
             /pGrid->U[ke][j][i].d;
        p0 /= pow((1.0 + r2),-0.75);

        cc_pos(pGrid,i,j,(ke+k),&x1,&x2,&x3);
        r2 = x1*x1 + x2*x2 + x3*x3;
        pGrid->U[ke+k][j][i].d = d0*pow((1.0 + r2),-0.75);
        pGrid->U[ke+k][j][i].M1 = 0.0;
        pGrid->U[ke+k][j][i].M2 = 0.0;
        pGrid->U[ke+k][j][i].M3 = 0.0;
        pGrid->U[ke+k][j][i].E = p0*pow((1.0 + r2),-0.75);
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i]   = pGrid->B1i[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].B1c = pGrid->U[ke-(k-1)][j][i].B1c;
      }
    }
  }


  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i]   = pGrid->B2i[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].B2c = pGrid->U[ke-(k-1)][j][i].B2c;
      }
    }
  }

/* k=ke+1 is not set for the interface field B3i, except obc_x3=1 */
  if (obc_x3 == 1) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+1][j][i] = 0.0;
      }
    }
  }
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i]   = pGrid->B3i[ke-(k-2)][j][i];
        pGrid->U[ke+k][j][i].B3c = pGrid->U[ke-(k-2)][j][i].B3c;
      }
    }
  }
#endif /* MHD */

  return;
}
