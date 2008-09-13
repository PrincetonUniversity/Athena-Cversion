#ifndef IONRAD_PROTOTYPES_H
#define IONRAD_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions from the following files:
 *   ionrad.c
 *   ionrad_3d.c
 *   ionrad_chemistry.c
 *   ionradplane_3d.c
 *   ionradpoint_3d.c
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

#ifdef ION_RADIATION
/*----------------------------------------------------------------------------*/
/* ionrad.c */
void ion_radtransfer_init_domain(Grid *pG, Domain *pD);
VGFun_t ion_radtransfer_init(Grid *pG, Domain *pD, int ires);

/*----------------------------------------------------------------------------*/
/* ionrad_3d.c */
void ion_radtransfer_3d(Grid *pG);
void ion_radtransfer_init_3d(Grid *pG, Domain *pD, int ires);
void ion_radtransfer_init_domain_3d(Grid *pG, Domain *pD);

/*----------------------------------------------------------------------------*/
/* ionrad_chemistry.c */
Real recomb_rate_coef(Real T);
Real coll_ion_rate_coef(Real T);
Real recomb_cool_rate_coef(Real T);
Real dmc_cool_rate(Real x, Real T);
Real osterbrock_cool_rate(Real T);
Real ki_cool_rate(Real T);
Real ki_heat_rate(void);
#endif /* ION_RADIATION */

#ifdef ION_RADPLANE
/*----------------------------------------------------------------------------*/
/* ionradplane_3d.c */
void add_radplane_3d(Grid *pGrid, int dir, Real flux);
void ion_radplane_init_domain_3d(Grid *pGrid, Domain *pDomain);
void get_ph_rate_plane(Real initflux, int dir, Real ***ph_rate, 
		       Grid *pGrid);
#endif /* ION_RADPLANE */

/*----------------------------------------------------------------------------*/
/* ionradpoint_3d.c */
#ifdef ION_RADPOINT
void get_ph_rate_point(Real s, Ray_Tree *tree, Real ***ph_rate,
		       Grid *pGrid);
void add_radpoint_3d(Grid *pG, Real x1, Real x2, Real x3, Real s);
void restore_radpoint_3d(Grid *pG, Real x1, Real x2, Real x3, Real s, 
			 int rebuild_ctr, float rotation[3][3]);
void refresh_trees(Grid *pG);
void ion_radpoint_init_domain_3d(Grid *pGrid, Domain *pDomain);
Rad_Ran2_State ion_radpoint_get_ranstate(void);
void ion_radpoint_set_ranstate(Rad_Ran2_State newstate);
void ion_radpoint_init_ranstate(void);
#endif /* ION_RADPOINT */

#endif /* IONRAD_PROTOTYPES_H */
