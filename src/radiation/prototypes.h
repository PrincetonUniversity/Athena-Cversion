#ifndef RADIATION_PROTOTYPES_H
#define RADIATION_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/radiation dir
 *============================================================================*/

#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"
#include "../config.h"

#ifdef RADIATION

/* init_radiation.c */

void init_radiation(MeshS *pM);
void radiation_temp_array_init(DomainS *pD);
void radiation_destruct(MeshS *pM);
void radgrid_destruct(RadGridS *pRG);

/* hydro_to_rad.c */

void hydro_to_rad(DomainS *pD);
void rad_to_hydro(DomainS *pD);

/* bvals_rad.c */
void bvals_rad(DomainS *pD);
void bvals_rad_init(MeshS *pM);

/* formal_solution.c */
void formal_solution(DomainS *pD);
#ifdef RAD_MULTIG
void output_mean_intensity_2d(RadGridS *pRG, int itr);
#endif

/* utils_rad.c */
void get_weights_parabolic(Real dtaum, Real dtaup, Real *edtau,
                           Real *a0, Real *a1, Real *a2);

void get_weights_linear(Real dtaum, Real *edtau, Real *a0, Real *a1);

/* all formal solutions algorithms contain the following functions */
void formal_solution_1d(RadGridS *pRG);
void formal_solution_1d_init(RadGridS *pRG);
void formal_solution_1d_destruct(void);

#ifdef RAD_MULTIG
/* gausseid_1d.c */
void gs_pass_pointers_to_mg_1d(Real *******psi0, Real **mu10, 
			       Real ****psiint0);
/* jacobi_1d.c */
void jacobi_pass_pointers_to_mg_1d(Real *******psi0, Real **mu10);

/* multigrid_1d */
void formal_solution_mg_1d(RadGridS *pRG);

/* multigrid_2d */
void formal_solution_mg_2d(RadGridS *pRG);

#endif /* RAD_MULTIG */

#endif /* RADIATION */
#endif /* RADIATION_PROTOTYPES_H */
