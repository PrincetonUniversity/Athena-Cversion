#include "../copyright.h"
/*==============================================================================
 * FILE: utils_rad.c
 *
 * PURPOSE: contains misc. functions require for computation of rad. transfer
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   get_weights_linear()     - 
 *   get_weights_parabolic()  -
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../prototypes.h"


#ifdef RADIATION

/*----------------------------------------------------------------------------*/
/* get weights using parabolic interpolation of source function
 */
void get_weights_parabolic(Real dtaum, Real dtaup, Real *edtau, Real *a0,
		 Real *a1, Real *a2)
{
  Real c0, c1, c2;
  Real dtausum;

  dtausum = dtaum + dtaup;
  (*edtau) = exp(-dtaum);

  //printf("%g %g\n",dtaum,dtaup);
  c0 = 1.0 - (*edtau);
  c1 = dtaum - c0;
  c2 = dtaum * dtaum - 2.0 * c1;
  
  (*a0) = c0 + (c2 - (dtaup + 2.0 * dtaum) * c1) / (dtaum * dtausum);
  (*a1) = (dtausum * c1 - c2) / (dtaum * dtaup);
  (*a2) = (c2 - dtaum * c1) / (dtaup * dtausum);


}

/*----------------------------------------------------------------------------*/
/* get weights using linear interpolation of source function
 */
void get_weights_linear(Real dtaum, Real *edtau, Real *a0, Real *a1)
{
  Real c0, c1;

  (*edtau) = exp(-dtaum);

  c0 = 1.0 - (*edtau);
  c1 = dtaum - c0;
  
  (*a0) = c0 - c1 / dtaum;
  (*a1) = c1 / dtaum;

}

#endif /* RADIATION */
