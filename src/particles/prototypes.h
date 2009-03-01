#ifndef PARTICLES_PROTOTYPES_H
#define PARTICLES_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the src/particles
 *   directory.
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/
/* particle.c */

#ifdef PARTICLES
void integrate_particle(Grid *pG);
void init_particle(Grid *pG, Domain *pD);
void particle_destruct(Grid *pG);
void particle_realloc(Grid *pG, long n);
void remove_ghost_particle(Grid *pG);
#ifdef FEEDBACK
void feedback(Grid *pG);
#endif
void shuffle(Grid *pG);
#endif

/* set_bvals_particle.c */
#ifdef PARTICLES
void set_bvals_particle(Grid *pG, Domain *pD);
#ifdef FARGO
void advect_particles(Grid *pG, Domain *pD);
#endif
void set_bvals_particle_init(Grid *pG, Domain *pD);
void set_bvals_particle_fun(enum Direction dir, VBCFun_t prob_bc);
void set_bvals_final_particle(Grid *pG, Domain *pD);
#endif

/* output_particle.c */
void init_output_particle(Grid *pG);
void particle_to_grid(Grid *pG, Domain *pD, Output *pout);
void destruct_particle_grid();
void dump_particle_binary(Grid *pG, Domain *pD);

#endif /* PARTICLES_PROTOTYPES_H */
