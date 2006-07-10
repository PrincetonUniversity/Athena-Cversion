#include "copyright.h"
/*==============================================================================
 * FILE: shu-osher.c
 *
 * PURPOSE: Problem generator for Shu-Osher shocktube test, involving
 *   interaction of a strong shock with a sine wave density distribution.  
 *
 * REFERENCE: C.W. Shu & S. Osher, "Efficient implementation of essentially
 *   non-oscillatory shock-capturing schemes, II", JCP, 83, 32 (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   problem - 
 *
 * PROBLEM USER FUNCTIONS: Must be included in every problem file, even if they
 *   are NoOPs and never used.  They provide user-defined functionality.
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifdef ISOTHERMAL
#error : The shu-osher test only works for adiabatic EOS.
#endif /* ISOTHERMAL */
#ifndef HYDRO
#error : The shu-osher test only works for hydro.
#endif /* HYDRO */

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid){
  int i,il,iu,js,ks;
  Real dl,pl,ul,vl,wl,x1,x2,x3;

/* Set up the grid bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  js = pGrid->js;
  ks = pGrid->ks;

  if (pGrid->Nx2 > 1 || pGrid->Nx3 > 1) {
    printf("Nx1=%i Nx2=%i Nx3=%i \n",pGrid->Nx1,pGrid->Nx2,pGrid->Nx3);
    ath_error("Shu Osher test only works for 1D problem in x1-direction\n");
  }

/* setup dependent variables */
  dl = 3.857143;
  ul = 2.629369;
  pl = 10.33333;
  wl = 0.0;
  vl = 0.0;

  for (i=il; i<=iu; i++) {
/* Calculate the cell center position of the cell i,j */
    cc_pos(pGrid,i,js,ks,&x1,&x2,&x3);

    if (x1 < -0.8) {
      pGrid->U[ks][js][i].d = dl;
      pGrid->U[ks][js][i].M1 = ul*dl;
      pGrid->U[ks][js][i].M2 = vl*dl;
      pGrid->U[ks][js][i].M3 = wl*dl;
      pGrid->U[ks][js][i].E = pl/Gamma_1 + 0.5*dl*(ul*ul + vl*vl + wl*wl);
    }
    else {
      pGrid->U[ks][js][i].d = 1.0 + 0.2*sin(5.0*PI*x1);
      pGrid->U[ks][js][i].M1 = 0.0;
      pGrid->U[ks][js][i].M2 = 0.0;
      pGrid->U[ks][js][i].M3 = 0.0;
      pGrid->U[ks][js][i].E = 1.0/Gamma_1;
    }
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

void Userwork_in_loop(Grid *pGrid)
{
}

void Userwork_after_loop(Grid *pGrid)
{
}
