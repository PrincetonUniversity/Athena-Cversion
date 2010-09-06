#include "../copyright.h"
/*==============================================================================
 * FILE: lr_states_dc_radiMHD.c
 *
 * PURPOSE: First order (donor cell, piecewise constant) spatial reconstruction.
 * This is noly used to resconstruct radiation energy density and flux
 *   The L/R-states at the left-interface in each cell are indexed i.
 *   W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 *   W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   lr_states()          - computes L/R states
 *   lr_states_init()     - NoOp function in this case
 *   lr_states_destruct() - NoOp function in this case
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef radiation_HD

/*----------------------------------------------------------------------------*/
/* This function must be and is only used for the radiation MHD code.

/* lr_states:
 * Input Arguments:
 *   W = PRIMITIVE variables at cell centers along 1-D slice
 *   Bxc = B in direction of slice at cell centers
 *   dtodx = dt/dx
 *   il,iu = lower and upper indices of zone centers in slice
 * W must be initialized over [il-1:iu+1]
 *
 * Output Arguments:
 *   Wl,Wr = L/R-states of PRIMITIVE variables at interfaces over [il:iu+1]
 */

void lr_states_dc_radiMHD(const GridS *pG, const Prim1DS W[], 
		const int il, const int iu,
               Prim1DS Wl[], Prim1DS Wr[])
{
  int i;
/* Only update the radiation part, NOT used for others */
  for (i=il; i<=iu+1; i++) {
    Wl[i].Er 	 = W[i-1].Er;
    Wl[i].Fluxr1 = W[i-1].Fluxr1;
    Wl[i].Fluxr2 = W[i-1].Fluxr2;
    Wl[i].Fluxr3 = W[i-1].Fluxr3;
    Wr[i].Er 	 = W[i  ].Er;
    Wr[i].Fluxr1 = W[i  ].Fluxr1;
    Wr[i].Fluxr2 = W[i  ].Fluxr2;
    Wr[i].Fluxr3 = W[i  ].Fluxr3;
  }

  return;
}


#endif /* radiation_MHD */
