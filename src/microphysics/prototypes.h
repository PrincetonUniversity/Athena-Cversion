#ifndef MICROPHYS_PROTOTYPES_H
#define MICROPHYS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the following files:
 *   braginskii.c
 *   cool.c
 *   integrate_diffusion.c
 *   isotropic_conduction.c
 *   resistivity.c
 *   viscosity.c
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* braginskii.c */
void brag_viscosity_2d(Grid *pG, Domain *pD);
void brag_viscosity_3d(Grid *pG, Domain *pD);
void brag_viscosity_init(int nx1, int nx2, int nx3);
void brag_viscosity_destruct(void);

/* cool.c */
Real KoyInut(const Real dens, const Real Press, const Real dt);

/* integrate_diffusion.c */
void integrate_explicit_diff(Grid *pGrid, Domain *pDomain);
void integrate_explicit_diff_init(Grid *pGrid, Domain *pDomain);
void integrate_explicit_diff_destruct(void);

/* isotropic_conduction.c */
void isoconduct(Grid *pG, Domain *pD);
void isoconduct_init(int nx1, int nx2, int nx3);
void isoconduct_destruct(void);

/* resistivity.c */
void ohmic_resistivity_1d(Grid *pG, Domain *pD);
void ohmic_resistivity_2d(Grid *pG, Domain *pD);
void ohmic_resistivity_3d(Grid *pG, Domain *pD);
void ohmic_resistivity_init(int nx1, int nx2, int nx3);
void ohmic_resistivity_destruct(void);

/* viscosity.c */
void ns_viscosity_1d(Grid *pG, Domain *pD);
void ns_viscosity_2d(Grid *pG, Domain *pD);
void ns_viscosity_3d(Grid *pG, Domain *pD);
void ns_viscosity_init(int nx1, int nx2, int nx3);
void ns_viscosity_destruct(void);

#endif /* MICROPHYS_PROTOTYPES_H */
