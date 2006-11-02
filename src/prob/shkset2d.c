#include "copyright.h"
/*==============================================================================
 * FILE: shkset2d.c
 *
 * PURPOSE: Sets up shock at angle to grid to test multidimensional algorithm.
 *   The grid angle atan(Ly/Lx) is fixed to be atan(0.5), or atan(1), and 
 *   Nx1/Nx2 must be the same ratio as Lx/Ly.  Uses the angle of the shock to
 *   remap ghost cells to the equivalent active grid cells, which requires
 *   that Nx1>32, using special function pointers.  The shock is initialized
 *   with reference to a coordinate system (x,y,z) with transformation rules to
 *   the code coordinate system (x1,x2,x3)
 *      x =  x1*cos(alpha) + x2*sin(alpha)
 *      y = -x1*sin(alpha) + x2*cos(alpha)
 *      z = x3
 *   This inverts to:
 *      x1 = x*cos(alpha) - y*sin(alpha)
 *      x2 = x*sin(alpha) + y*cos(alpha)
 *      x3 = z
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.
 *============================================================================*/

void shkset2d_iib(Grid *pGrid);
void shkset2d_oib(Grid *pGrid);
void shkset2d_ijb(Grid *pGrid);
void shkset2d_ojb(Grid *pGrid);

/* Make size of box and dimension of unit cell (r1 x r2) static globals so they
 * can be accessed by boundary value functions */
static Real Lx,Ly;
static int r1,r2;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,ix1,ix2,nx1,nx2,gcd;
  Real angle, sin_a, cos_a; /* Angle the shock makes with the x1-direction */
  Gas ql, qr;
  Real dl,pl,ul,vl,wl,dr,pr,ur,vr,wr;
#ifdef MHD
  Real bxl,byl,bzl,bxr,byr,bzr;
#endif /* MHD */
  div_t id;   /* structure containing remainder and quotient */

/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
  int dll, dlr, drr, drl;
  Real afl_lx, afr_lx, afl_rx, afr_rx;
  Real afl_ly, afr_ly, afl_ry, afr_ry;
  Real vfl, vfr, B1r, B2r;

  if (pGrid->Nx3 > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  nx1 = (ie-is)+1;
  nx2 = (je-js)+1;
  if ((nx1 == 1) || (nx2 == 1)) {
    fprintf(stderr,"[shkset2d]: This problem can only be run in 2D\n");
    exit(EXIT_FAILURE);
  }

/* Compute greatest common divisor of nx1,nx2.  The size of the "unit cell"
 * is nx1/gcd by nx2/gcd */

  if((gcd = ath_gcd(nx1,nx2)) < 10){
    fprintf(stderr,"[shkset2d]: Greatest Common Divisor (nx1,nx2) = %d\n",gcd);
    exit(EXIT_FAILURE);
  }

  id = div(nx1,gcd);
  r1 = id.quot;
  if(id.rem != 0){
    fprintf(stderr,"[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx1,gcd,id.rem);
    exit(EXIT_FAILURE);
  }

  id = div(nx2,gcd);
  r2 = id.quot;
  if(id.rem != 0){
    fprintf(stderr,"[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx2,gcd,id.rem);
    exit(EXIT_FAILURE);
  }

  printf("The unit cell is (%d,1,%d) grid cells in size\n",r1,r2);

/* Compute angle initial interface makes to the grid */

  Lx = pGrid->Nx1*pGrid->dx1;
  Ly = pGrid->Nx2*pGrid->dx2;
  if(Lx == Ly){
    cos_a = sin_a = sqrt(0.5);
  }
  else{
    angle = atan((double)(Lx/Ly));
    sin_a = sin(angle);
    cos_a = cos(angle);
  }

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

/* Initialize ql rotated to the (x1,x2,x3) coordinate system */
  ql.d   = dl;
  ql.M1  = dl*(ul*cos_a - vl*sin_a);
  ql.M2  = dl*(ul*sin_a + vl*cos_a);
  ql.M3  = dl*wl;
#ifdef MHD
  ql.B1c = bxl*cos_a - byl*sin_a;
  ql.B2c = bxl*sin_a + byl*cos_a;
  ql.B3c = bzl;
#endif
#ifdef ADIABATIC
  ql.E   = pl/Gamma_1
#ifdef MHD
    + 0.5*(bxl*bxl + byl*byl + bzl*bzl)
#endif
    + 0.5*(ul*ul + vl*vl + wl*wl)*dl;
#endif

/* Initialize qr rotated to the (x1,x2,x3) coordinate system */
  qr.d   = dr;
  qr.M1  = dr*(ur*cos_a - vr*sin_a);
  qr.M2  = dr*(ur*sin_a + vr*cos_a);
  qr.M3  = dr*wr;
#ifdef MHD
  qr.B1c = bxr*cos_a - byr*sin_a;
  qr.B2c = bxr*sin_a + byr*cos_a;
  qr.B3c = bzr;
#endif
#ifdef ADIABATIC
  qr.E   = pr/Gamma_1
#ifdef MHD
    + 0.5*(bxr*bxr + byr*byr + bzr*bzr)
#endif
    + 0.5*(ur*ur + vr*vr + wr*wr)*dr;
#endif

/* Initialize the grid */

  for (k=kl; k<=ku; k++) {
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->jdisp;
      for (i=0; i<=ie+nghost; i++) {
	ix1 = i + pGrid->idisp;

/* cell is completely in the left state */
	if((drr = r2*(ix1) + r1*(ix2) - gcd*r1*r2) <= 0){
	  pGrid->U[k][j][i] = ql;
#ifdef MHD
	  pGrid->B1i[k][j][i] = ql.B1c;
	  pGrid->B2i[k][j][i] = ql.B2c;
	  pGrid->B3i[k][j][i] = ql.B3c;
#endif /* MHD */
	}
/* cell is completely in the right state */
	else if((dll = r2*(ix1-1) + r1*(ix2-1) - gcd*r1*r2) >= 0){
	  pGrid->U[k][j][i] = qr;
#ifdef MHD
	  pGrid->B1i[k][j][i] = qr.B1c;
	  pGrid->B2i[k][j][i] = qr.B2c;
	  pGrid->B3i[k][j][i] = qr.B3c;
#endif /* MHD */
	}
/* The more complicated case of a cell  split by the interface boundary */
	else{
	  dlr = r2*(ix1-1) + r1*(ix2) - gcd*r1*r2;

	  if(dlr < 0){ /* The boundary hits the right y-face */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = (Real)(-dlr)/(Real)(r2);
	    afr_ry = 1.0 - afl_ry;
	  }
	  else if(dlr > 0){ /* The boundary hits the left x-face */
	    afl_lx = (Real)(-dll)/(Real)(r1);
	    afr_lx = 1.0 - afl_lx;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }
	  else{ /* dlr == 0.0, The boundary hits the grid cell corner */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }

	  drl = r2*(ix1) + r1*(ix2-1) - gcd*r1*r2;

	  if(drl < 0){ /* The boundary hits the right x-face */
	    afl_rx = (Real)(-drl)/(Real)(r1);
	    afr_rx = 1.0 - afl_rx;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }
	  else if(drl > 0){ /* The boundary hits the left y-face */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = (Real)(-dll)/(Real)(r2);
	    afr_ly = 1.0 - afl_ly;
	  }
	  else{ /* drl == 0.0, The boundary hits the grid cell corner */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }

/* The boundary hits both x-interfaces */
	  if(dlr > 0 && drl < 0){ 
	    vfl = 0.5*(afl_lx + afl_rx);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both y-interfaces */
	  else if(dlr < 0 && drl > 0){ 
	    vfl = 0.5*(afl_ly + afl_ry);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both grid cell corners */
	  else if(dlr == 0 && drl == 0){ 
	    vfl = vfr = 0.5;
	  }
/* The boundary hits the left x- and left y-interface */
	  else if(dlr > 0 && drl > 0){
	    vfl = 0.5*afl_lx*afl_ly;
	    vfr = 1.0 - vfl;
	  }
/* dlr<0 && drl<0:  The boundary hits the right x- and right y-interface */
	  else{ 
	    vfr = 0.5*afr_rx*afr_ry;
	    vfl = 1.0 - vfr;
	  }

/* Initialize the x- and y-interface magnetic fields */
#ifdef MHD
	  pGrid->B1i[k][j][i] = afl_lx*ql.B1c + afr_lx*qr.B1c;
	  B1r              = afl_rx*ql.B1c + afr_rx*qr.B1c;

	  pGrid->B2i[k][j][i] = afl_ly*ql.B2c + afr_ly*qr.B2c;
	  B2r              = afl_ry*ql.B2c + afr_ry*qr.B2c;

	  pGrid->B3i[k][j][i] = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */

/* Initialize the volume averaged quantities */
	  pGrid->U[k][j][i].d  = vfl*ql.d + vfr*qr.d;
	  pGrid->U[k][j][i].M1 = vfl*ql.M1 + vfr*qr.M1;
	  pGrid->U[k][j][i].M2 = vfl*ql.M2 + vfr*qr.M2;
	  pGrid->U[k][j][i].M3 = vfl*ql.M3 + vfr*qr.M3;
#ifdef MHD
	  pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + B1r);
	  pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + B2r);
	  pGrid->U[k][j][i].B3c = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E  = vfl*ql.E + vfr*qr.E;
#endif
	}
      }
    }
  }

/* Set boundary value function pointers */

  set_bvals_fun(left_x1,shkset2d_iib);
  set_bvals_fun(left_x2,shkset2d_ijb);
  set_bvals_fun(right_x1,shkset2d_oib);
  set_bvals_fun(right_x2,shkset2d_ojb);

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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*-----------------------------------------------------------------------------
 * Function shkset2d_iib: sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_iib(Grid *pGrid)
{
  const int is = pGrid->is;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost + r2;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx3 > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U  [k][j][is-i] = pGrid->U  [k][j-r2][is-i+r1];
#ifdef MHD
	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j-r2][is-i+r1];
	pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j-r2][is-i+r1];
	pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j-r2][is-i+r1];
#endif
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------------
 * Function shkset2d_oib: same for oib
 */

void shkset2d_oib(Grid *pGrid)
{
  const int ie = pGrid->ie;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost - r2;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx3 > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that i=ie+1 is not a boundary condition for the interface field B1i */

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U[k][j][ie+i] = pGrid->U[k][j+r2][ie+i-r1];
#ifdef MHD
	if(i>1) pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j+r2][ie+i-r1];
	pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j+r2][ie+i-r1];
	pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j+r2][ie+i-r1];
#endif
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------------
 * Function shkset2d_ijb: same for ijb
 */

void shkset2d_ijb(Grid *pGrid)
{
  const int js = pGrid->js;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost + r1;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  if (pGrid->Nx3 > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U  [k][js-j][i] = pGrid->U  [k][js-j+r2][i-r1];
#ifdef MHD
	pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js-j+r2][i-r1];
	pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js-j+r2][i-r1];
	pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js-j+r2][i-r1];
#endif
      }
    }
  }
  return;
}

/*-----------------------------------------------------------------------------
 * Function shkset2d_ojb: same for ojb
 */

void shkset2d_ojb(Grid *pGrid)
{
  const int je = pGrid->je;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost - r1;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }

  if (pGrid->Nx3 > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U[k][je+j][i] = pGrid->U[k][je+j-r2][i+r1];
#ifdef MHD
	pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-r2][i+r1];
	if(j>1) pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je+j-r2][i+r1];
	pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je+j-r2][i+r1];
#endif
      }
    }
  }
  return;
}
