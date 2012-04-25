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
void radgrid_destruct(RadGridS *pRG);

/* ray_tracing.c */
#ifdef RAY_TRACING
void ray_trace(DomainS *pD);
#ifdef MPI_PARALLEL
void allocate_working_array_ray_tracing(RayGridS *pRayG);
void destruct_working_array_ray_tracing(void);
#endif
#endif /* RAY_TRACING */

/* radtrans_dt */

Real radtrans_dt(DomainS *pD);

/* hydro_to_rad.c */

void hydro_to_rad(DomainS *pD, int iflag);
void rad_to_hydro(DomainS *pD);

/* bvals_rad.c */
void bvals_rad(DomainS *pD, int ifs, int ife);
void bvals_rad_init(MeshS *pM);
void bvals_rad_trans_fun(DomainS *pD, enum BCDirection dir, VRGIFun_t prob_bc);

/* bvals_rad_shear.c */
void ShearingSheet_Rad_ix1(DomainS *pD, int ifs, int ife);
void ShearingSheet_Rad_ox1(DomainS *pD, int ifs, int ife);

/* dump_intensity_vtk.c */
void dump_ix1_vtk(MeshS *pM, OutputS *pOut);
void dump_ox1_vtk(MeshS *pM, OutputS *pOut);
void dump_ix2_vtk(MeshS *pM, OutputS *pOut);
void dump_ox2_vtk(MeshS *pM, OutputS *pOut);
void dump_ix3_vtk(MeshS *pM, OutputS *pOut);
void dump_ox3_vtk(MeshS *pM, OutputS *pOut);

/* output_spec.c */
void output_spec(MeshS *pM);

/* formal_solution.c */
void formal_solution(DomainS *pD);

/* utils_rad.c */
void interp_quad_chi(Real chi0, Real chi1, Real chi2, Real *chi);
void interp_quad_source(Real dtaum, Real dtaup, Real *edtau, Real *a0,
			Real *a1, Real *a2, Real S0, Real S1, Real S2);
void interp_quad_source_slope_lim(Real dtaum, Real dtaup, Real *edtau, Real *a0,
				  Real *a1, Real *a2, Real S0, Real S1, Real S2);
void get_weights_parabolic(Real dtaum, Real dtaup, Real *edtau,
                           Real *a0, Real *a1, Real *a2);

void get_weights_linear(Real dtaum, Real *edtau, Real *a0, Real *a1);

int permutation(int i, int j, int k, int **pl, int np);

/* all 1D formal solutions algorithms contain the following functions */
void formal_solution_1d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_1d_init(RadGridS *pRG);
void formal_solution_1d_destruct(void);

/* all 2D formal solutions algorithms contain the following functions */
void formal_solution_2d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_2d_init(RadGridS *pRG);
void formal_solution_2d_destruct(void);

/* all 3D formal solutions algorithms contain the following functions */
void formal_solution_3d(RadGridS *pRG, Real *dSrmax, int ifr);
void formal_solution_3d_init(RadGridS *pRG);
void formal_solution_3d_destruct(void);

#endif /* RADIATION_TRANSFER */
#endif /* RADIATION_TRANSFER_PROTOTYPES_H */
