#include "copyright.h"
/*==============================================================================
 * FILE: integrate.c
 *
 * PURPOSE: Contains public functions to set integrator and source terms.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_init()        - set pointer to integrate function based on dim
 *   integrate_destruct()    - call destruct integrate function based on dim
 *   cons_pot_fun_enroll()   - enroll function that computes source term that
 *                             can be described by a conservative potential
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/* dimension of calculation (determined at runtime) */
static int dim=0;

/*----------------------------------------------------------------------------*/
/* integrate_init: initialize integrator; VGFun_t defined in athena.h   */

VGFun_t integrate_init(int Nx1, int Nx2, int Nx3)
{

/* Calculate the dimensions  */
  dim = 0;
  if(Nx1 > 1) dim++;
  if(Nx2 > 1) dim++;
  if(Nx3 > 1) dim++;

/* set function pointer to appropriate integrator based on dimensions */
  switch(dim){
  case 1:
    if(Nx1 <= 1) break;
    integrate_init_1d(Nx1);
    return integrate_1d;
  case 2:
    if(Nx3 > 1) break;
    integrate_init_2d(Nx1,Nx2);
    return integrate_2d;
  case 3:
    integrate_init_3d(Nx1,Nx2,Nx3);
    return integrate_3d;
  }

  ath_error("[integrate_init]: Unsupported dim. Nx1=%d, Nx2=%d, Nx3=%d\n",
	    Nx1,Nx2,Nx3);

/* This is never executed, but generates a warning on some compilers. */
  return NULL;
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct:  */

void integrate_destruct()
{
  switch(dim){
  case 1:
    integrate_destruct_1d();
    return;
  case 2:
    integrate_destruct_2d();
    return;
  case 3:
    integrate_destruct_3d();
    return;
  }

  ath_error("[integrate_destruct]: Grid dimension = %d\n",dim);
}

/*----------------------------------------------------------------------------*/
/* cons_pot_fun_enroll:  */

void cons_pot_fun_enroll(ConsPotFun_t pfun)
{
  switch(dim){
  case 1:
    cons_pot_fun_enroll_1d(pfun);
    return;
  case 2:
    cons_pot_fun_enroll_2d(pfun);
    return;
  case 3:
    cons_pot_fun_enroll_3d(pfun);
    return;
  }

  ath_error("[cons_pot_fun_enroll]: Grid dimension = %d\n",dim);
}
