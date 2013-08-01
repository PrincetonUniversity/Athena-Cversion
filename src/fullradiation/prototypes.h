#ifndef FULL_RADIATION_TRANSFER_PROTOTYPES_H
#define FULL_RADIATION_TRANSFER_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/fullradiation dir
 *============================================================================*/

#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"
#include "../config.h"

#ifdef FULL_RADIATION_TRANSFER

/* init_fullradiation.c */

VDFun_t init_fullradiation(MeshS *pM);
void fullradiation_destruct(MeshS *pM);


/* hydro_to_rad.c */

void hydro_to_fullrad(DomainS *pD);


/* bvals_fullrad.c */
void bvals_fullrad(DomainS *pD);
void bvals_fullrad_init(MeshS *pM);
void bvals_fullrad_destruct();
void bvals_fullrad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc);

/* utils_fullrad.c */
/* Update the momentums of specific intensity for each cell */

void UpdateRT(DomainS *pD);
void ReduceVelocity(const Real sigma, const Real ds, Real *alpha);
void SpecialMatrix(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N);
void SpecialMatrix2(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N);
void SpecialMatrix3(Real *Ma, Real *Mb, Real *Mc, Real *Md, Real *RHS, Real **lN,  Real **UN, const int N);
void convert_angle(const Real x2, const Real miux0, const Real miuy0, Real *miux, Real *miuy);
/* FullRT_flux.c */

/* piece linear flux */
void flux_PLM(Real r[3]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[3], Real imhalf[1]);
void flux_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[5], Real imhalf[1]);
void lrstate_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real ds __attribute__((unused)), Real imu[5], Real iLeft[1], Real iRight[1]);
void flux_AdvJ(Real r[3]  __attribute__((unused)), const int dir __attribute__((unused)), Real *tempJ, Real *tempV, int nstart, int nend, Real ds, Real dt, Real *tempAdv);

int permutation(int i, int j, int k, int **pl, int np);

/* fullRT_2d.c */
void fullRT_2d_init(RadGridS *pRG);
void fullRT_2d_destruct(void);

void fullRT_2d(DomainS *pD);

/* fullRT_3d.c */
void fullRT_3d(DomainS *pD);
void fullRT_3d_init(RadGridS *pRG);
void fullRT_3d_destruct(void);

/* output_spec.c */
void output_spec(MeshS *pM);

#endif /* RADIATION_TRANSFER */
#endif /* RADIATION_TRANSFER_PROTOTYPES_H */
