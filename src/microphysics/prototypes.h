#ifndef MICROPHYS_PROTOTYPES_H
#define MICROPHYS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/microphysics dir
 *============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* conduction.c */
#ifdef THERMAL_CONDUCTION
void conduction(DomainS *pD);
void conduction_init(MeshS *pM);
void conduction_destruct(void);
#endif

/* cool.c */
Real KoyInut(const Real dens, const Real Press, const Real dt);

/* diff_dt.c */
Real diff_dt(MeshS *pM);

/* integrate_diffusion.c */
void integrate_diff(MeshS *pM);
void integrate_diff_init(MeshS *pM);
void integrate_diff_destruct(void);

/* resistivity.c */
#ifdef RESISTIVTY
void resistivity(DomainS *pD);
void resistivity_init(MeshS *pM);
void resistivity_destruct(void);
#endif

/* viscosity.c */
#ifdef VISCOSITY
void viscosity(DomainS *pD);
void viscosity_init(MeshS *pM);
void viscosity_destruct(void);
#endif

#endif /* MICROPHYS_PROTOTYPES_H */
