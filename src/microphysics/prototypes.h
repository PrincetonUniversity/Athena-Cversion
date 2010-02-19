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

/* anisotropic_conduction.c */
#ifdef ANISOTROPIC_CONDUCTION
void anisoconduct_2d(DomainS *pD);
void anisoconduct_3d(DomainS *pD);
void anisoconduct_init(MeshS *pM);
void anisoconduct_destruct(void);
#endif

/* braginskii.c */
#ifdef BRAGINSKII
void brag_viscosity_2d(DomainS *pD);
void brag_viscosity_3d(DomainS *pD);
void brag_viscosity_init(MeshS *pM);
void brag_viscosity_destruct(void);
#endif

/* cool.c */
Real KoyInut(const Real dens, const Real Press, const Real dt);

/* hall.c */
#ifdef HALL_MHD
void hall_resistivity_1d(DomainS *pD);
void hall_resistivity_2d(DomainS *pD);
void hall_resistivity_3d(DomainS *pD);
void hall_resistivity_init(MeshS *pM);
void hall_resistivity_destruct(void);
#endif

/* integrate_diffusion.c */
#ifdef EXPLICIT_DIFFUSION
void integrate_diff(MeshS *pM);
void integrate_diff_init(MeshS *pM);
void integrate_diff_destruct(void);
#endif

/* isotropic_conduction.c */
#ifdef ISOTROPIC_CONDUCTION
void isoconduct(DomainS *pD);
void isoconduct_init(MeshS *pM);
void isoconduct_destruct(void);
#endif

/* resistivity.c */
#ifdef OHMIC
void ohmic_resistivity_1d(DomainS *pD);
void ohmic_resistivity_2d(DomainS *pD);
void ohmic_resistivity_3d(DomainS *pD);
void ohmic_resistivity_init(MeshS *pM);
void ohmic_resistivity_destruct(void);
#endif

/* viscosity.c */
#ifdef NAVIER_STOKES
void ns_viscosity_1d(DomainS *pD);
void ns_viscosity_2d(DomainS *pD);
void ns_viscosity_3d(DomainS *pD);
void ns_viscosity_init(MeshS *pM);
void ns_viscosity_destruct(void);
#endif

#endif /* MICROPHYS_PROTOTYPES_H */
