#include "copyright.h"
/*==============================================================================
 * FILE: lw_implode.c
 *
 * PURPOSE: Problem generator for square implosion problem.
 *
 * REFERENCE: R. Liska & B. Wendroff, SIAM J. Sci. Comput., 25, 995 (2003)
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MHD
#error : The lw_implode problem only works for hydro.
#endif

/* computes difference d{i,j}-d{j,i} to test if solution is symmetric */
static Real expr_diff_d(const Grid *pG, const int i, const int j, const int k);

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(Grid *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int a;
  Real d_in, p_in, d_out, p_out;

/* Set up the grid bounds for initializing the grid */
  if (pGrid->Nx1 <= 1 || pGrid->Nx2 <= 1) {
    ath_error("[problem]: This problem requires Nx1 > 1, Nx2 > 1\n");
  }

  d_in = par_getd("problem","d_in");
  p_in = par_getd("problem","p_in");

  d_out = par_getd("problem","d_out");
  p_out = par_getd("problem","p_out");

/* Initialize the grid */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;

	a = (2*(i + pGrid->idisp) - 1)*pGrid->Nx2
	  + (2*(j + pGrid->jdisp) - 1)*pGrid->Nx1;

	if(a > pGrid->Nx1*pGrid->Nx2) {
	  pGrid->U[k][j][i].d  = d_out;
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E = p_out/Gamma_1;
#endif
	} else {
	  pGrid->U[k][j][i].d  = d_in;
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E = p_in/Gamma_1;
#endif
	}
      }
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

void Userwork_in_loop(Grid *pGrid)
{
}

void Userwork_after_loop(Grid *pGrid)
{
}

void problem_write_restart(Grid *pG, FILE *fp){
  return;
}

void problem_read_restart(Grid *pG, FILE *fp){
  return;
}

Gasfun_t get_usr_expr(const char *expr){
  if(strcmp(expr,"diff_d")==0) return expr_diff_d;
  return NULL;
}

static Real expr_diff_d(const Grid *pG, const int i, const int j, const int k)
{
 return (pG->U[k][j][i].d - pG->U[k][i][j].d);
}
