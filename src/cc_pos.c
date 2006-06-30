#include "copyright.h"
/*==============================================================================
 * FILE: cc_pos.c
 *
 * PURPOSE: Function to compute (x1,x2,x3) positions of center of cell i,j,k.  
 *   In a nested grid, each Grid structure is a patch in a larger computational
 *   domain (with the exception of the level0 grid).  The displacement of the
 *   origin of the Grid from the origin of the computational domain (level0 
 *   grid) is x1_{disp} = idisp*dx1.  Furthermore, the origin of the level0
 *   grid can be displaced by a distance x1_{0} from the origin of x1.  Thus,
 *   the x1 position of the center of cell i (x1_{cc,i}) in any level Grid is
 *            x1_{cc,i} = x1_{0} + ((i + idisp) + 0.5)*dx1
 *   Similarly for x2 and x3.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   cc_pos() -
 *============================================================================*/

#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* cc_pos:  */

void cc_pos(const Grid *pG, const int i, const int j,const int k,
	    Real *px1, Real *px2, Real *px3)
{
  *px1 = pG->x1_0 + ((i + pG->idisp) + 0.5)*pG->dx1;
  *px2 = pG->x2_0 + ((j + pG->jdisp) + 0.5)*pG->dx2;
  *px3 = pG->x3_0 + ((k + pG->kdisp) + 0.5)*pG->dx3;
  return;
}
