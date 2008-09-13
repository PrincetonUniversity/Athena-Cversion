#include "../copyright.h"
#define IONRAD_C
/*==============================================================================
 * FILE: ionrad.c
 *
 * PURPOSE: Contains function to control ionizing radiative transfer routines.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ion_radtransfer_init() - sets pointer to appropriate ionizing
 *                            radiative transfer function
 *   ion_radtransfer_init_domain() - sets domain information for ionizing
 *                                   radiative transfer module
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "ionrad.h"

#ifdef ION_RADIATION
static int dim=0;

void ion_radtransfer_init_domain(Grid *pGrid, Domain *pDomain) {

  /* Calcualte the dimensionality and error check */
  dim = 0;
  if(pGrid->Nx1 > 1) dim++;
  if(pGrid->Nx2 > 1) dim++;
  if(pGrid->Nx3 > 1) dim++;

  switch(dim) {
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_domain_3d(pGrid, pDomain);
    return;
  }

  ath_error("[ion_radtransfer_init_domain]: Unsupported dim. Nx1=%d, Nx2=%d, Nx3=%d\n",
	    pGrid->Nx1,pGrid->Nx2,pGrid->Nx3);
}

VGFun_t ion_radtransfer_init(Grid *pGrid, Domain *pDomain, int ires){

  /* Calcualte the dimensionality and error check */
  dim = 0;
  if(pGrid->Nx1 > 1) dim++;
  if(pGrid->Nx2 > 1) dim++;
  if(pGrid->Nx3 > 1) dim++;

  switch(dim){
  case 1: break;
  case 2: break;
  case 3:
    ion_radtransfer_init_3d(pGrid, pDomain, ires);
    return ion_radtransfer_3d;
  }

  ath_error("[ion_radtransfer_init]: Unsupported dim. Nx1=%d, Nx2=%d, Nx3=%d\n",
	    pGrid->Nx1,pGrid->Nx2,pGrid->Nx3);

  /* This is never executed, but lack of a return statement generates
     a warning on some compilers. */
  return NULL;
}

#endif /* ION_RADIATION */
