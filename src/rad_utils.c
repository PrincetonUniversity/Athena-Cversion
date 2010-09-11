#include "./copyright.h"
/*==============================================================================
 * FILE: rad_utils.c
 *
 * PURPOSE: A set of functions used to calculate radiation quantities.
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

#ifdef RADIATION
/*----------------------------------------------------------------------------*/
/* 
 *   Input Arguments:
 *     U = Conserved variable
 */

Real eff_sound(const Cons1DS U, Real dt)
{
	Real pressure, temperature, velocity, enthalpy, SEE, Alpha, TEnergy;
	Real aeff;


	TEnergy = U.E;
	pressure = (TEnergy - 0.5 * U.Mx * U.Mx / U.d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (U.d * Ridealgas);
	velocity = U.Mx / U.d;


	enthalpy = Gamma * TEnergy / U.d - (Gamma - 1.0) * velocity * velocity / 2.0;
/* The Source term */
	SEE = 4.0 * Sigmaa * temperature * temperature * temperature * (Gamma - 1.0)/ U.d;
	Alpha = (1.0 - exp(-Pratio * Cratio * SEE * dt/2.0))/(Pratio * Cratio* SEE * dt/2.0);
	aeff = -(Gamma - 1.0) * velocity * velocity/2.0 + Alpha * (Gamma - 1.0) * enthalpy 
			+ (1.0 - Alpha) * (temperature + (Gamma - 1.0) * velocity * velocity/2.0);
	aeff = sqrt(aeff); 

	return aeff;
}
#endif 
