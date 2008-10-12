#ifndef MICROPHYS_PROTOTYPES_H
#define MICROPHYS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the following files:
 *   braginskii.c
 *   cool.c
 *   resistivity.c
 *   viscosity.c
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* braginskii.c */
void braginskii_3d(Grid *pG, Domain *pD);
void braginskii_init_3d(int nx1, int nx2, int nx3);

/* cool.c */
Real KoyInut(const Real dens, const Real Press, const Real dt);

/* resistivity.c */
void resistivity_3d(Grid *pG, Domain *pD);

/* viscosity.c */
void viscosity_3d(Grid *pG, Domain *pD);
void viscosity_init_3d(int nx1, int nx2, int nx3);

#endif /* MICROPHYS_PROTOTYPES_H */
