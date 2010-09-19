#include "./copyright.h"
/*==============================================================================
 * FILE: rad_utils.c
 *
 * PURPOSE: A set of functions used to calculate rad_hydro quantities.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   eff_sound() - To calculate the effective sound speed
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "./defs.h"
#include "./athena.h"
#include "./globals.h"
#include "prototypes.h"

#ifdef rad_hydro
/*----------------------------------------------------------------------------*/
/* 
 *   Input Arguments:
 *     U = Conserved variable
 *    The effective sound speed is calculated as conserved variable formula
 */

Real eff_sound(const Prim1DS W, Real dt)
{
	Real aeff, temperature, SPP, Alpha;

	temperature = W.P / (W.d * R_ideal);
	

	SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a 
		* temperature * temperature * temperature / (W.d * R_ideal);

	if(fabs(SPP * dt * 0.5) > 0.001)
	Alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
	else 
	Alpha = 1.0 + 0.25 * SPP * dt;
	/* In case SPP * dt  is small, use expansion expression */	


	aeff = ((Gamma - 1.0) * Alpha + 1.0) * W.P / W.d;

	aeff = sqrt(aeff); 

	return aeff;
}
#endif 
