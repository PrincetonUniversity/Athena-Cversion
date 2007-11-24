#include "copyright.h"
/*==============================================================================
 * FILE: lr_states_dc.c
 *
 * PURPOSE: First order (donor cell, piecewise constant) spatial reconstruction.
 *   The L/R-states at the left-interface in each cell are indexed i.
 *   U_{L,i-1/2} is denoted by Ul[i  ];   U_{R,i-1/2} is denoted by Ur[i  ]
 *   U_{L,i+1/2} is denoted by Ul[i+1];   U_{R,i+1/2} is denoted by Ur[i+1]
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   lr_states()          - computes L/R states
 *   lr_states_init()     - NoOp function in this case
 *   lr_states_destruct() - NoOp function in this case
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifdef FIRST_ORDER

/*----------------------------------------------------------------------------*/
/* lr_states:
 * Input Arguments:
 *   U1d = CONSERVED variables at cell centers along 1-D slice
 *   Bxc = B in direction of slice at cell centers
 *   dt = timestep
 *   dtodx = dt/dx
 *   il,iu = lower and upper indices of zone centers in slice
 * U1d must be initialized over [il-1:iu+1]
 *
 * Output Arguments:
 *   Ul,Ur = L/R-states of CONSERVED variables at interfaces over [il:iu+1]
 */

void lr_states(const Prim1D W[], const Real Bxc[],
	       const Real dt, const Real dtodx, const int il, const int iu,
	       Prim1D Wl[], Prim1D Wr[])
{
  int i;

  for (i=il; i<=iu+1; i++) {
    Wl[i] = W1d[i-1];
    Wr[i] = W1d[i  ];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* lr_states_init:  NoOp for first order, but included for compatibility
 *   with integrator (needed for 2nd and 3rd order). */

void lr_states_init(int nx1, int nx2, int nx3)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* lr_states_destruct:  NoOp for first order, but included for compatibility
 *   with integrator (needed for 2nd and 3rd order). */

void lr_states_destruct(void)
{
  return;
}

#endif /* FIRST_ORDER */
