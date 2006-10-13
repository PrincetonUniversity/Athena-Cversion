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

void problem(Grid *pGrid)
{
  int i,il,iu,j,jl,ju,k,kl,ku,middle;
  int shk_dir; /* Shock direction: {1,2,3} -> {x1,x2,x3} */
  Real dl,pl,ul,vl,wl,dr,pr,ur,vr,wr;
#ifdef MHD
  Real bxl,byl,bzl,bxr,byr,bzr;
#endif /* MHD */

/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  dl = par_getd("problem","dl");
#ifdef ADIABATIC
  pl = par_getd("problem","pl");
#endif
  ul = par_getd("problem","v1l");
  vl = par_getd("problem","v2l");
  wl = par_getd("problem","v3l");
#ifdef MHD
  bxl = par_getd("problem","b1l");
  byl = par_getd("problem","b2l");
  bzl = par_getd("problem","b3l");
#endif

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  dr = par_getd("problem","dr");
#ifdef ADIABATIC
  pr = par_getd("problem","pr");
#endif
  ur = par_getd("problem","v1r");
  vr = par_getd("problem","v2r");
  wr = par_getd("problem","v3r");
#ifdef MHD
  bxr = par_getd("problem","b1r");
  byr = par_getd("problem","b2r");
  bzr = par_getd("problem","b3r");
#endif

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
  case 1:  /* shock in 1-direction  */
    middle = il - 1 + (int)(0.5*(iu-il+1));
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=middle; i++) {
          pGrid->U[k][j][i].d  = dl;
          pGrid->U[k][j][i].M1 = ul*dl;
          pGrid->U[k][j][i].M2 = vl*dl;
          pGrid->U[k][j][i].M3 = wl*dl;
#ifdef MHD
          pGrid->B1i[k][j][i] = bxl;
          pGrid->B2i[k][j][i] = byl;
          pGrid->B3i[k][j][i] = bzl;
          pGrid->U[k][j][i].B1c = bxl;
          pGrid->U[k][j][i].B2c = byl;
          pGrid->U[k][j][i].B3c = bzl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pl/Gamma_1 
#ifdef MHD
	  + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
	  + 0.5*dl*(ul*ul + vl*vl + wl*wl);
#endif
        }
        for (i=middle+1; i<=iu; i++) {
          pGrid->U[k][j][i].d  = dr;
          pGrid->U[k][j][i].M1 = ur*dr;
          pGrid->U[k][j][i].M2 = vr*dr;
          pGrid->U[k][j][i].M3 = wr*dr;
#ifdef MHD
          pGrid->B1i[k][j][i] = bxr;
          pGrid->B2i[k][j][i] = byr;
          pGrid->B3i[k][j][i] = bzr;
          pGrid->U[k][j][i].B1c = bxr;
          pGrid->U[k][j][i].B2c = byr;
          pGrid->U[k][j][i].B3c = bzr;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pr/Gamma_1
#ifdef MHD
	  + 0.5*(bxr*bxr + byr*byr + bzr*bzr) 
#endif
	  + 0.5*dr*(ur*ur + vr*vr + wr*wr);
#endif
        }
      }
    }
    break;
  case 2:  /* shock in 2-direction  */
    middle = jl - 1 + (int)(0.5*(ju-jl+1));
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=middle; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->U[k][j][i].d  = dl;
          pGrid->U[k][j][i].M1 = wl*dl;
          pGrid->U[k][j][i].M2 = ul*dl;
          pGrid->U[k][j][i].M3 = vl*dl;
#ifdef MHD
          pGrid->B1i[k][j][i] = bzl;
          pGrid->B2i[k][j][i] = bxl;
          pGrid->B3i[k][j][i] = byl;
          pGrid->U[k][j][i].B1c = bzl;
          pGrid->U[k][j][i].B2c = bxl;
          pGrid->U[k][j][i].B3c = byl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pl/Gamma_1 
#ifdef MHD
            + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
	  + 0.5*dl*(ul*ul + vl*vl + wl*wl);
#endif
        }
      }
      for (j=middle+1; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->U[k][j][i].d  = dr;
          pGrid->U[k][j][i].M1 = wr*dr;
          pGrid->U[k][j][i].M2 = ur*dr;
          pGrid->U[k][j][i].M3 = vr*dr;
#ifdef MHD
          pGrid->B1i[k][j][i] = bzr;
          pGrid->B2i[k][j][i] = bxr;
          pGrid->B3i[k][j][i] = byr;
          pGrid->U[k][j][i].B1c = bzr;
          pGrid->U[k][j][i].B2c = bxr;
          pGrid->U[k][j][i].B3c = byr;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pr/Gamma_1
#ifdef MHD
	  + 0.5*(bxr*bxr + byr*byr + bzr*bzr) 
#endif
	  + 0.5*dr*(ur*ur + vr*vr + wr*wr);
#endif
        }
      }
    }
    break;
  case 3:  /* shock in 3-direction  */
    middle = kl - 1 + (int)(0.5*(ku-kl+1));
    for (k=kl; k<=middle; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->U[k][j][i].d  = dl;
          pGrid->U[k][j][i].M1 = vl*dl;
          pGrid->U[k][j][i].M2 = wl*dl;
          pGrid->U[k][j][i].M3 = ul*dl;
#ifdef MHD
          pGrid->B1i[k][j][i] = byl;
          pGrid->B2i[k][j][i] = bzl;
          pGrid->B3i[k][j][i] = bxl;
          pGrid->U[k][j][i].B1c = byl;
          pGrid->U[k][j][i].B2c = bzl;
          pGrid->U[k][j][i].B3c = bxl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pl/Gamma_1 
#ifdef MHD
            + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
	  + 0.5*dl*(ul*ul + vl*vl + wl*wl);
#endif
        }
      }
    }
    for (k=middle+1; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          pGrid->U[k][j][i].d  = dr;
          pGrid->U[k][j][i].M1 = vr*dr;
          pGrid->U[k][j][i].M2 = wr*dr;
          pGrid->U[k][j][i].M3 = ur*dr;
#ifdef MHD
          pGrid->B1i[k][j][i] = byr;
          pGrid->B2i[k][j][i] = bzr;
          pGrid->B3i[k][j][i] = bxr;
          pGrid->U[k][j][i].B1c = byr;
          pGrid->U[k][j][i].B2c = bzr;
          pGrid->U[k][j][i].B3c = bxr;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = pr/Gamma_1
#ifdef MHD
	  + 0.5*(bxr*bxr + byr*byr + bzr*bzr) 
#endif
	  + 0.5*dr*(ur*ur + vr*vr + wr*wr);
#endif
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
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(Grid *pG, FILE *fp){
  return;
}

void problem_read_restart(Grid *pG, FILE *fp){
  return;
}

Gasfun_t get_usr_expr(const char *expr){
  return NULL;
}

void Userwork_in_loop(Grid *pGrid){
  return;
}

void Userwork_after_loop(Grid *pGrid){
  return;
}
