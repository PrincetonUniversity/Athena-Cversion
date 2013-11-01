#ifndef RADIATION_TRANSFER_PROTOTYPES_H
#define RADIATION_TRANSFER_PROTOTYPES_H 
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

#ifdef RADIATION_TRANSFER

/* init_radiation.c */

void init_radiation(MeshS *pM);
void radiation_destruct(MeshS *pM);
void formal_solution_init(DomainS *pD);

/* ray_tracing.c */
#ifdef RAY_TRACING
void ray_trace(DomainS *pD, const int outflag);
void init_ray_tracing(RadGridS *pRG);
void destruct_ray_tracing(RadGridS *pRG);
Real raytrace_to_radtrans_default(RadGridS *pRG, int ifray, int ifrad);
#endif /* RAY_TRACING */

/* radtrans_dt */

Real radtrans_dt(DomainS *pD);

/* hydro_to_rad.c */

void hydro_to_rad(DomainS *pD, const int outflag);
void rad_to_hydro(DomainS *pD);

/* bvals_rad.c */
void bvals_rad(DomainS *pD, int ifr, const int outflag);
void bvals_rad_init(MeshS *pM);
void bvals_rad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc);

/* bvals_rad_shear.c */
void bvals_rad_shear_init(MeshS *pM);
void bvals_rad_shear_destruct(void);
void ShearingSheet_Rad_ix1(DomainS *pD, int ifr);
void ShearingSheet_Rad_ox1(DomainS *pD, int ifr);

/* output_intensity_vtk.c */
void output_ix1_vtk(MeshS *pM, OutputS *pOut);
void output_ox1_vtk(MeshS *pM, OutputS *pOut);
void output_ix2_vtk(MeshS *pM, OutputS *pOut);
void output_ox2_vtk(MeshS *pM, OutputS *pOut);
void output_ix3_vtk(MeshS *pM, OutputS *pOut);
void output_ox3_vtk(MeshS *pM, OutputS *pOut);

/* output_spec.c */
void output_rad_mesh(MeshS *pM);

/* formal_solution.c */
void formal_solution(DomainS *pD, const int outflag, const int fstflag);

/* utils_rad.c */
void interp_quad_chi(Real chi0, Real chi1, Real chi2, Real *chi);
void interp_quad_source_slope_lim(Real dtaum, Real dtaup, Real *edtau, Real *a0,
				  Real *a1, Real *a2, Real S0, Real S1, Real S2);
void interp_quad_source(Real dtaum, Real dtaup, Real *edtau, Real *a0,
			Real *a1, Real *a2);
void interp_linear_soource(Real dtaum, Real *edtau, Real *a0, Real *a1);

/* angles.c */
void init_angles(RadGridS *pRG, const int qmeth, const int outflag);

/* all 1D formal solutions algorithms contain the following functions */
void formal_solution_1d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_1d_init(DomainS *pD);
void formal_solution_1d_destruct(void);

/* all 2D formal solutions algorithms contain the following functions */
void formal_solution_2d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_2d_init(DomainS *pD);
void formal_solution_2d_destruct(void);

/* all 3D formal solutions algorithms contain the following functions */
void formal_solution_3d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_3d_init(DomainS *pD);
void formal_solution_3d_destruct(void);

#endif /* RADIATION_TRANSFER */
#endif /* RADIATION_TRANSFER_PROTOTYPES_H */
