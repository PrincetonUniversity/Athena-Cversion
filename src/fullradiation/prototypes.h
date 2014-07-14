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


void Absorption(const int nelements, const Real rho, const Real dtsigma, const Real Told, const Real Jold, Real Vel[3], Real vdotn[], Real vvdotnn[], Real wmu[], Real imu[], Real **ABCm, Real *Tnew);
void RadAsource(const int i, const int j, const int k, const int N, const Real Tcoef, const Real Tnew, RadGridS *pRG, GridS *pG,  Real **Coefn, Real **sol);
void RadSsource(const int i, const int j, const int k, const int N, RadGridS *pRG, GridS *pG, Real **Coefn, Real **sol);
void UpdateOpacity(DomainS *pD);
void GetVelguess(DomainS *pD);
void EstimateVel(const int il, const int iu, const int jl, const int ju, const int kl, const int ku, const int flag, RadGridS *pRG, GridS *pG);
void GetSpeedfactor(DomainS *pD);
void ComptIntensity(DomainS *pD);
void UpdateHcomp(Real *Hnew, Real Hold[3], const Real dtsigma, Real Vel[3]);
void UpdateKcomp(Real *Knew, Real Kold[6], Real Hnew[3], Real Vel[3], const Real Jnew, const Real dtsigma, const Real Source);


/* bvals_fullrad.c */
void bvals_fullrad(DomainS *pD);
void bvals_fullrad_init(MeshS *pM);
void bvals_fullrad_destruct();
void bvals_fullrad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc);

/* utils_fullrad.c */
/* Update the momentums of specific intensity for each cell */

void UpdateRT(DomainS *pD);
void CalMoment(const int il, const int iu, const int jl, const int ju, const int kl, const int ku, RadGridS *pRG);
void ReduceVelocity(const Real sigma, const Real ds, Real *alpha);
void SpecialMatrix(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N);
void SpecialMatrix2(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N);
void SpecialMatrix3(const int N, Real **Ma, Real *Md, Real *RHS,  Real *lN1, Real *lN2, Real *lN3,  Real *UN1, Real *UN2, Real *UN3);
void AbsorptionMatrix(const int N, Real **Ma, Real *Md, Real *RHS);
void LinearAbsorptionMatrix(const int N, Real **Ma, Real *RHS);
void convert_angle(const Real x2, const Real miux0, const Real miuy0, Real *miux, Real *miuy);
/* FullRT_flux.c */

/* piece linear flux */
void flux_PLM(Real r[3]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[3], Real imhalf[1]);
void flux_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[5], Real imhalf[1]);
void lrstate_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real ds __attribute__((unused)), Real imu[5], Real iLeft[1], Real iRight[1]);
void flux_AdvJ(const int nf, const int N, Real *r  __attribute__((unused)), const int dir __attribute__((unused)), Real ***tempJ, Real ***tempV, int nstart,int nend,Real ds, Real dt, Real ***tempAdv);

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

/* output_intensity_vtk.c */
/* only output intensity at the faces */
void output_ix1_vtk(MeshS *pM, OutputS *pOut);
void output_ox1_vtk(MeshS *pM, OutputS *pOut);
void output_ix2_vtk(MeshS *pM, OutputS *pOut);
void output_ox2_vtk(MeshS *pM, OutputS *pOut);
void output_ix3_vtk(MeshS *pM, OutputS *pOut);
void output_ox3_vtk(MeshS *pM, OutputS *pOut);

/* dump_intensity_vtk.c */
/* Dump intensities for all angles at all grids */
void dump_intensity_vtk(MeshS *pM, OutputS *pOut);



#endif /* RADIATION_TRANSFER */
#endif /* RADIATION_TRANSFER_PROTOTYPES_H */
