#include "copyright.h"
/*==============================================================================
 * FILE: radMHD1d.c
 *
 * PURPOSE: Problem generator to test the radiation MHD code. 
 *  Development is underway. To be modified later 
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

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  int shift;

/* Parse global variables of unit ratio */
#ifdef RADIATION_HYDRO
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
	
#endif

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Initialize the grid including the ghost cells.  */


    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

/* Initialize conserved (and  the primitive) variables in Grid */
	if(i < il + (iu - il)/2){

          pGrid->U[k][j][i].d  = 1.0;
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = 1/(Gamma - 1.0);
	  pGrid->U[k][j][i].Er = 1.0;
	  pGrid->U[k][j][i].Fr1 = 1.0;
#endif
	}
	else {
          pGrid->U[k][j][i].d  = 0.125;
          pGrid->U[k][j][i].M1 = 0.0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = 0.1/(Gamma - 1.0);
#endif
	  pGrid->U[k][j][i].Er = 0.1 * 0.1 * 0.1 * 0.1;
	  pGrid->U[k][j][i].Fr1 = 0.1;

	}




#ifdef MHD
          pGrid->B1i[k][j][i] = 0.0;
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;
          pGrid->U[k][j][i].B1c = 0.0;
          pGrid->U[k][j][i].B2c = 0.0;
          pGrid->U[k][j][i].B3c = 0.0;
#endif
#ifdef RADIATION_HYDRO
/*	  pGrid->U[k][j][i].Er = 1.0;
	  pGrid->U[k][j][i].Fr1 = 0.0;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;
*/
	  pGrid->U[k][j][i].Edd_11 = 1.0; /* Set to be a constant in 1D. To be modified later */		
#endif
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
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
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
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}
