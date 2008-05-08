#include "copyright.h"
/*==============================================================================
 * FILE: blast.c
 *
 * PURPOSE: Problem generator for spherical blast wave problem.
 *
 * REFERENCE: P. Londrillo & L. Del Zanna, "High-order upwind schemes for
 *   multidimensional MHD", ApJ, 530, 508 (2000), and references therein.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid, Domain *pDomain)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real pressure,prat,rad,da,pa,ua,va,wa,x1,x2,x3;
  Real bxa,bya,bza,b0=0.0;
  Real rin;
  double theta;

  rin = par_getd("problem","radius");
  pa  = par_getd("problem","pamb");
  prat = par_getd("problem","prat");
#ifdef MHD
  b0 = par_getd("problem","b0");
#endif
  theta = (PI/180.0)*par_getd("problem","angle");

/* setup uniform ambient medium with spherical over-pressured region */

  da = 1.0;
  ua = 0.0;
  va = 0.0;
  wa = 0.0;
  bxa = b0*cos(theta);
  bya = b0*sin(theta);
  bza = 0.0;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->U[k][j][i].d  = da;
	pGrid->U[k][j][i].M1 = da*ua;
	pGrid->U[k][j][i].M2 = da*va;
	pGrid->U[k][j][i].M3 = da*wa;
#ifdef MHD
	pGrid->B1i[k][j][i] = bxa;
	pGrid->B2i[k][j][i] = bya;
	pGrid->B3i[k][j][i] = bza;
	pGrid->U[k][j][i].B1c = bxa;
	pGrid->U[k][j][i].B2c = bya;
	pGrid->U[k][j][i].B3c = bza;
	if (i == ie && ie > is) pGrid->B1i[k][j][i+1] = bxa;
	if (j == je && je > js) pGrid->B2i[k][j+1][i] = bya;
	if (k == ke && ke > ks) pGrid->B3i[k+1][j][i] = bza;
#endif
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	rad = sqrt(x1*x1 + x2*x2 + x3*x3);
	pressure = pa;
	if (rad < rin) pressure = prat*pa;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E = pressure/Gamma_1 
#ifdef MHD
          + 0.5*(bxa*bxa + bya*bya + bza*bza)
#endif
          + 0.5*da*(ua*ua + va*va + wa*wa);
#else
	if (rad < rin) pGrid->U[k][j][i].d = prat*da;
#endif /* ISOTHERMAL */
      }
    }
  }
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
