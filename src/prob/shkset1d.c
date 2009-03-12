#include "copyright.h"
/*==============================================================================
 * FILE: shkset1d.c
 *
 * PURPOSE: Problem generator for 1-D Riemann problems.  Initial discontinuity
 *   is located so there are equal numbers of cells to the left and right (at
 *   center of grid based on integer index).  Initializes plane-parallel shock
 *   along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i,il,iu,j,jl,ju,k,kl,ku,middle;
  int shk_dir; /* Shock direction: {1,2,3} -> {x1,x2,x3} */
  Prim1D W1d, Wl, Wr;
  Cons1D U1d, Ul, Ur;
#ifdef MHD
  Real Bxl, Bxr;
#endif

/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  Wl.d = par_getd("problem","dl");
#ifdef ADIABATIC
  Wl.P = par_getd("problem","pl");
#endif
  Wl.Vx = par_getd("problem","v1l");
  Wl.Vy = par_getd("problem","v2l");
  Wl.Vz = par_getd("problem","v3l");
#ifdef MHD
  Bxl = par_getd("problem","b1l");
  Wl.By = par_getd("problem","b2l");
  Wl.Bz = par_getd("problem","b3l");
#endif
#if (NSCALARS > 0)
  Wl.r[0] = par_getd("problem","r[0]l");
#endif

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  Wr.d = par_getd("problem","dr");
#ifdef ADIABATIC
  Wr.P = par_getd("problem","pr");
#endif
  Wr.Vx = par_getd("problem","v1r");
  Wr.Vy = par_getd("problem","v2r");
  Wr.Vz = par_getd("problem","v3r");
#ifdef MHD
  Bxr = par_getd("problem","b1r");
  Wr.By = par_getd("problem","b2r");
  Wr.Bz = par_getd("problem","b3r");
  if (Bxr != Bxl) ath_error(0,"[shkset1d] L/R values of Bx not the same\n");
#endif
#if (NSCALARS > 0)
  Wr.r[0] = par_getd("problem","r[0]r");
#endif

  Prim1D_to_Cons1D(&Ul,&Wl MHDARG( , &Bxl));
  Prim1D_to_Cons1D(&Ur,&Wr MHDARG( , &Bxr));

/* Parse shock direction */
  shk_dir = par_geti("problem","shk_dir");
  if (shk_dir != 1 && shk_dir != 2 && shk_dir != 3) {
    ath_error("[problem]: shk_dir = %d must be either 1,2 or 3\n",shk_dir);
  }

/* Set up the index bounds for initializing the grid */
  if (pGrid->Nx1 > 1) {
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  }
  else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  if (pGrid->Nx2 > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx3 > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Initialize the grid including the ghost cells.  Discontinuity is always
 * located in middle of grid (at zone index half way between is,ie),
 * regardless of value of x-coordinate */

  switch(shk_dir) {
/*--- shock in 1-direction ---------------------------------------------------*/
  case 1:  /* shock in 1-direction  */
    middle = il - 1 + (int)(0.5*(iu-il+1));
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

/* set primitive and conserved variables to be L or R state */
          if (i <= middle) {
            W1d = Wl;
            U1d = Ul;
          } else {
            W1d = Wr;
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.Mx;
          pGrid->U[k][j][i].M2 = U1d.My;
          pGrid->U[k][j][i].M3 = U1d.Mz;
#ifdef MHD
          pGrid->B1i[k][j][i] = Bxl;
          pGrid->B2i[k][j][i] = U1d.By;
          pGrid->B3i[k][j][i] = U1d.Bz;
          pGrid->U[k][j][i].B1c = Bxl;
          pGrid->U[k][j][i].B2c = U1d.By;
          pGrid->U[k][j][i].B3c = U1d.Bz;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif

#ifdef SPECIAL_RELATIVITY
          pGrid->W[k][j][i].d  = W1d.d;
          pGrid->W[k][j][i].V1 = W1d.Vx;
          pGrid->W[k][j][i].V2 = W1d.Vy;
          pGrid->W[k][j][i].V3 = W1d.Vz;
#ifdef MHD
          pGrid->W[k][j][i].B1c = Bxl;
          pGrid->W[k][j][i].B2c = W1d.By;
          pGrid->W[k][j][i].B3c = W1d.Bz;
#endif
#ifdef ADIABATIC
          pGrid->W[k][j][i].P = W1d.P;
#endif
#if (NSCALARS > 0)
          pGrid->W[k][j][i].r[0] = W1d.r[0];
#endif
#endif /* SPECIAL_RELATIVITY */
        }
      }
    }
    break;

/*--- shock in 2-direction ---------------------------------------------------*/
  case 2:  /* shock in 2-direction  */
    middle = jl - 1 + (int)(0.5*(ju-jl+1));
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

/* set primitive variables to be L or R state */
          if (j <= middle) {
            W1d = Wl;
            U1d = Ul;
          } else {
            W1d = Wr;
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.Mz;
          pGrid->U[k][j][i].M2 = U1d.Mx;
          pGrid->U[k][j][i].M3 = U1d.My;
#ifdef MHD
          pGrid->B1i[k][j][i] = U1d.Bz;
          pGrid->B2i[k][j][i] = Bxl;
          pGrid->B3i[k][j][i] = U1d.By;
          pGrid->U[k][j][i].B1c = U1d.Bz;
          pGrid->U[k][j][i].B2c = Bxl;
          pGrid->U[k][j][i].B3c = U1d.By;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif

#ifdef SPECIAL_RELATIVITY
          pGrid->W[k][j][i].d  = W1d.d;
          pGrid->W[k][j][i].V1 = W1d.Vz;
          pGrid->W[k][j][i].V2 = W1d.Vx;
          pGrid->W[k][j][i].V3 = W1d.Vy;
#ifdef MHD
          pGrid->W[k][j][i].B1c = W1d.Bz;
          pGrid->W[k][j][i].B2c = Bxl;
          pGrid->W[k][j][i].B3c = W1d.By;
#endif
#ifdef ADIABATIC
          pGrid->W[k][j][i].P = W1d.P;
#endif
#if (NSCALARS > 0)
          pGrid->W[k][j][i].r[0] = W1d.r[0];
#endif
#endif /* SPECIAL_RELATIVITY */
        }
      }
    }
    break;

/*--- shock in 3-direction ---------------------------------------------------*/
  case 3:  /* shock in 3-direction  */
    middle = kl - 1 + (int)(0.5*(ku-kl+1));
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

/* set primitive variables to be L or R state */
          if (k <= middle) {
            W1d = Wl;
            U1d = Ul;
          } else {
            W1d = Wr;
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.My;
          pGrid->U[k][j][i].M2 = U1d.Mz;
          pGrid->U[k][j][i].M3 = U1d.Mx;
#ifdef MHD
          pGrid->B1i[k][j][i] = U1d.By;
          pGrid->B2i[k][j][i] = U1d.Bz;
          pGrid->B3i[k][j][i] = Bxl;
          pGrid->U[k][j][i].B1c = U1d.By;
          pGrid->U[k][j][i].B2c = U1d.Bz;
          pGrid->U[k][j][i].B3c = Bxl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif

#ifdef SPECIAL_RELATIVITY
          pGrid->W[k][j][i].d  = W1d.d;
          pGrid->W[k][j][i].V1 = W1d.Vy;
          pGrid->W[k][j][i].V2 = W1d.Vz;
          pGrid->W[k][j][i].V3 = W1d.Vx;
#ifdef MHD
          pGrid->W[k][j][i].B1c = W1d.By;
          pGrid->W[k][j][i].B2c = W1d.Bz;
          pGrid->W[k][j][i].B3c = Bxl;
#endif
#ifdef ADIABATIC
          pGrid->W[k][j][i].P = W1d.P;
#endif
#if (NSCALARS > 0)
          pGrid->W[k][j][i].r[0] = W1d.r[0];
#endif
#endif /* SPECIAL_RELATIVITY */
        }
      }
    }
  break;
  default:
    ath_error("[shkset1d]: invalid shk_dir = %i\n",shk_dir);
  }

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

GVDFun_t get_usr_gasvshift(const char *name)
{
  return NULL;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}
