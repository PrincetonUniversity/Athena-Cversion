#include "../copyright.h"
/*==============================================================================
 * FILE: ray_tracing.c
 *
 * PURPOSE: Contatins functions for handling ray tracing on Cartesian,
 * grid-alinged rays.  Models attenuation of "external" radiation that can
 * be approximated by parallel rays (e.g. a distant point source).  When
 * scattering is present, computes source term for the "diffuse" emission.
 * Current implementation assumes rays aligned with x1 direction with source
 * of radiation at ix1 boundary.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ray_trace()  - compute external radiation
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RAY_TRACING


#endif /* RAY_TRACING */
