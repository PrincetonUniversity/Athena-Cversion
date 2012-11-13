#ifndef PARTICLES_PROTOTYPES_H
#define PARTICLES_PROTOTYPES_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions in the /src/particles dir      */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

#ifdef PARTICLES

/* bvals_particle.c */
void bvals_particle(DomainS *pD);
#ifdef FARGO
void advect_particles(DomainS *pD);
#endif
void bvals_particle_init(MeshS *pM);
void bvals_particle_fun(enum BCDirection dir, VGFun_t prob_bc);
void bvals_final_particle(MeshS *pM);

/* dump_particle_history.c */
void dump_particle_history(MeshS *pM, OutputS *pOut);
void dump_parhistory_enroll();

/* exchange.c */
void exchange_gpcouple(DomainS *pD, short lab);
void exchange_gpcouple_init(MeshS *pM);
void exchange_gpcouple_fun(enum BCDirection dir, VGFun_t prob_bc);
void exchange_gpcouple_destruct(MeshS *pM);

/* init_particle.c */
void init_particle(MeshS *pM);
void particle_destruct(MeshS *pM);
void particle_realloc(GridS *pG, long n);

/* integrators_particle.c */
void Integrate_Particles(DomainS *pD);
void int_par_exp   (GridS *pG, GrainS *curG, Real3Vect cell1,
                              Real *dv1, Real *dv2, Real *dv3, Real *ts);
void int_par_semimp(GridS *pG, GrainS *curG, Real3Vect cell1,
                              Real *dv1, Real *dv2, Real *dv3, Real *ts);
void int_par_fulimp(GridS *pG, GrainS *curG, Real3Vect cell1,
                              Real *dv1, Real *dv2, Real *dv3, Real *ts);
#ifdef FEEDBACK
void feedback_predictor(DomainS *pD);
void feedback_corrector(GridS *pG, GrainS *gri, GrainS *grf, Real3Vect cell1,
                              Real dv1, Real dv2, Real dv3, Real ts);
#endif

/* output_particle.c */
void particle_to_grid(DomainS *pD, PropFun_t par_prop);
void dump_particle_binary(MeshS *pM, OutputS *pOut);
int  property_all(const GrainS *gr, const GrainAux *grsub);

/* utils_particle.c */
void get_gasinfo(GridS *pG);

void getwei_linear(GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);
void getwei_TSC   (GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);
void getwei_QP    (GridS *pG, Real x1, Real x2, Real x3, Real3Vect cell1,
                              Real weight[3][3][3], int *is, int *js, int *ks);

int getvalues(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
#ifndef FEEDBACK
                        Real *rho, Real *u1, Real *u2, Real *u3, Real *cs
#else
             Real *rho, Real *u1,  Real *u2, Real *u3, Real *cs, Real *stiff
#endif
);

Real get_ts_epstein(GridS *pG, int type, Real rho, Real cs, Real vd);
Real get_ts_general(GridS *pG, int type, Real rho, Real cs, Real vd);
Real get_ts_fixed  (GridS *pG, int type, Real rho, Real cs, Real vd);

#ifdef FEEDBACK
void feedback_clear(GridS *pG);
void distrFB_pred(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
#ifndef BAROTROPIC
                            Real3Vect fb, Real stiffness, Real Elosspar
#else
                            Real3Vect fb, Real stiffness
#endif
);
void distrFB_corr(GridS *pG, Real weight[3][3][3], int is, int js, int ks,
                                             Real3Vect fb, Real Elosspar);
#endif

void shuffle(GridS *pG);

#endif /* PARTICLES */
#endif /* PARTICLES_PROTOTYPES_H */
