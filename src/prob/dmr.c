#include "copyright.h"
/*==============================================================================
 * FILE: dmr.c
 *
 * PURPOSE: Problem generator for double Mach reflection test.  Only works for
 *   genuinely 2D problems in X1-X2 plane.
 *
 * REFERENCE: P. Woodward & P. Colella, "The numerical simulation of 
 *   two-dimensional fluid flow with strong shocks", JCP, 54, 115, sect. IVc.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error : The dmr test only works for adiabatic EOS.
#endif /* ISOTHERMAL */
#ifndef HYDRO
#error : The dmr test only works for hydro.
#endif /* HYDRO */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * dmrbv_iib() - sets BCs on L-x1 (left edge) of grid.  
 * dmrbv_ijb() - sets BCs on L-x2 (bottom edge) of grid.  
 * dmrbv_ojb() - sets BCs on R-x2 (top edge) of grid.  
 *============================================================================*/

void dmrbv_iib(Grid *pGrid);
void dmrbv_ijb(Grid *pGrid);
void dmrbv_ojb(Grid *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
int i=0,j=0;
int is,ie,js,je,ks;
Real d0,e0,u0,v0,x1_shock,x1,x2,x3;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  if (pGrid->Nx1 == 1 || pGrid->Nx2 == 1) {
    ath_error("[dmr]: this test only works with Nx1 & Nx2 > 1\n");
  }
  if (pGrid->Nx3 > 1) {
    ath_error("[dmr]: this test only works for 2D problems, with Nx3=1\n");
  }

/* Initialize shock using parameters defined in Woodward & Colella */

  d0 = 8.0;
  e0 = 291.25;
  u0 =  8.25*sqrt(3.0)/2.0;
  v0 = -8.25*0.5;
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
      x1_shock = 0.1666666666 + x2/sqrt((double)3.0);
/* upstream conditions */
      pGrid->U[ks][j][i].d = 1.4;
      pGrid->U[ks][j][i].E = 2.5;
      pGrid->U[ks][j][i].M1 = 0.0;
      pGrid->U[ks][j][i].M2 = 0.0;
/* downstream conditions */
      if (x1 < x1_shock) {
        pGrid->U[ks][j][i].d = d0;
        pGrid->U[ks][j][i].E = e0 + 0.5*d0*(u0*u0+v0*v0);
        pGrid->U[ks][j][i].M1 = d0*u0;
        pGrid->U[ks][j][i].M2 = d0*v0;
      }
    }
  }

/* Set boundary value function pointers */

  set_bvals_mhd_fun(left_x1,dmrbv_iib);
  set_bvals_mhd_fun(left_x2,dmrbv_ijb);
  set_bvals_mhd_fun(right_x2,dmrbv_ojb);
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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*-----------------------------------------------------------------------------
 * dmrbv_iib: sets boundary condition on left X boundary (iib) for dmr test
 * Note quantities at this boundary are held fixed at the downstream state
 */

void dmrbv_iib(Grid *pGrid)
{
int i=0,j=0;
int is,ie,js,je,ks,jl,ju;
Real d0,e0,u0,v0;

  d0 = 8.0;
  e0 = 291.25;
  u0 =  8.25*sqrt(3.0)/2.0;
  v0 = -8.25*0.5;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  for (j=jl; j<=ju;  j++) {
    for (i=1;  i<=nghost;  i++) {
      pGrid->U[ks][j][is-i].d  = d0;
      pGrid->U[ks][j][is-i].M1 = d0*u0;
      pGrid->U[ks][j][is-i].M2 = d0*v0;
      pGrid->U[ks][j][is-i].E  = e0 + 0.5*d0*(u0*u0+v0*v0);
    }
  }
}

/*-----------------------------------------------------------------------------
 * dmrbv_ijb: sets boundary condition on lower Y boundary (ijb) for dmr test
 * Note quantaties at this boundary are held fixed at the downstream state for
 * x1 < 0.16666666, and are reflected for x1 > 0.16666666
 */

void dmrbv_ijb(Grid *pGrid)
{
int i=0,j=0;
int is,ie,js,je,ks,il,iu;
Real d0,e0,u0,v0,x1,x2,x3;

  d0 = 8.0;
  e0 = 291.25;
  u0 =  8.25*sqrt(3.0)/2.0;
  v0 = -8.25*0.5;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (j=1;  j<=nghost;  j++) {
    for (i=il; i<=iu; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
      if (x1 < 0.1666666666) {
/* fixed at downstream state */
        pGrid->U[ks][js-j][i].d  = d0;
        pGrid->U[ks][js-j][i].M1 = d0*u0;
        pGrid->U[ks][js-j][i].M2 = d0*v0;
        pGrid->U[ks][js-j][i].E  = e0 + 0.5*d0*(u0*u0+v0*v0);
      } else {
/* reflected */
        pGrid->U[ks][js-j][i].d  = pGrid->U[ks][js+(j-1)][i].d;
        pGrid->U[ks][js-j][i].M1 = pGrid->U[ks][js+(j-1)][i].M1;
        pGrid->U[ks][js-j][i].M2 = -pGrid->U[ks][js+(j-1)][i].M2;
        pGrid->U[ks][js-j][i].E  = pGrid->U[ks][js+(j-1)][i].E;
      }
    }
  }
}

/*-----------------------------------------------------------------------------
 * dmrbv_ojb: sets TIME-DEPENDENT boundary condition on upper Y boundary (ojb)
 * for dmr test.  Quantaties at this boundary are held fixed at the downstream
 * state for x1 < 0.16666666+v1_shock*time, and at the upstream state for
 * x1 > 0.16666666+v1_shock*time
 */

void dmrbv_ojb(Grid *pGrid)
{
int i=0,j=0;
int is,ie,js,je,ks,il,iu;
Real d0,e0,u0,v0,x1_shock,x1,x2,x3;

  d0 = 8.0;
  e0 = 291.25;
  u0 =  8.25*sqrt(3.0)/2.0;
  v0 = -8.25*0.5;
  x1_shock = 0.1666666666 + (1.0 + 20.0*pGrid->time)/sqrt(3.0);

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks;
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  for (j=1;  j<=nghost;  j++) {
    for (i=il; i<=iu; i++) {
      cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);
      if (x1 < x1_shock) {
/* fixed at downstream state */
        pGrid->U[ks][je+j][i].d  = d0;
        pGrid->U[ks][je+j][i].M1 = d0*u0;
        pGrid->U[ks][je+j][i].M2 = d0*v0;
        pGrid->U[ks][je+j][i].E  = e0 + 0.5*d0*(u0*u0+v0*v0);
      } else {
/* fixed at upstream state */
        pGrid->U[ks][je+j][i].d  = 1.4;
        pGrid->U[ks][je+j][i].M1 = 0.0;
        pGrid->U[ks][je+j][i].M2 = 0.0;
        pGrid->U[ks][je+j][i].E  = 2.5;
      }
    }
  }
}
