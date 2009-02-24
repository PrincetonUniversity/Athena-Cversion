#ifndef INTEGRATORS_PROTOTYPES_H
#define INTEGRATORS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the src/integrators directory
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/
/* integrate.c */

VGDFun_t integrate_init(int Nx1, int Nx2, int Nx3);
void integrate_destruct(void);

/*----------------------------------------------------------------------------*/
/* integrate_1d_ctu.c and integrate_1d_vl.c contain the same functions */

void integrate_destruct_1d(void);
void integrate_init_1d(int Nx1);
void integrate_1d(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* integrate_2d_ctu.c and integrate_2d_vl.c contain the same functions */

void integrate_destruct_2d(void);
void integrate_init_2d(int Nx1, int Nx2);
void integrate_2d(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* integrate_3d_ctu.c and integrate_3d_vl.c contain the same functions */

void integrate_destruct_3d(void);
void integrate_init_3d(int Nx1, int Nx2, int Nx3);
void integrate_3d(Grid *pG, Domain *pD);

#endif /* INTEGRATORS_PROTOTYPES_H */
