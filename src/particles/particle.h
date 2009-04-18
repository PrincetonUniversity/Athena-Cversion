#ifndef PARTICLE_H
#define PARTICLE_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: particle.h
 *
 * PURPOSE: Global variables for all functionsin in the src/particles
 *   directory.
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/
#ifdef PARTICLES

/* 3D Arrays to store gas quantities at (n+1/2) step */
Real   ***grid_d;		/* gas density */
Vector ***grid_v;		/* gas velocities */
#ifndef ISOTHERMAL
Real   ***grid_cs;		/* gas sound speed */
#endif

/* coefficient for the calculation of a/lambda */
Real alamcoeff;

/* an array of particle solid density times particle size in normalized unit */
Real *grrhoa;

/* left and right limit of grid indices */
int ilp,iup, jlp,jup, klp,kup;

/* left and right limit of grid boundary */
Real x1lpar, x1upar, x2lpar, x2upar, x3lpar, x3upar;

/* number of neighbouring cells involved in 1D interpolation */
int ncell;

/* value of the fixed stopping time */
Real tstop0;

#ifdef SHEARING_BOX
/* shear velocity */
Real vshear;

/* left and right limit in x1 boundary, for feedback */
int is0, ie0;
#endif

#endif /* PARTICLES */

#endif /* PARTICLE_H */
