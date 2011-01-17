#ifndef RSOLVERS_PROTOTYPES_H
#define RSOLVERS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/rsolvers directory.
 *============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* esystem_*.c */
void esys_roe_iso_hyd(const Real v1, const Real v2, const Real v3,
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4]);

void esys_roe_adb_hyd(const Real v1, const Real v2, const Real v3,
  const Real h, Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5]);

void esys_roe_iso_mhd(const Real d, const Real v1, const Real v2,
  const Real v3, const Real b1, const Real b2, const Real b3,
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6]);

void esys_roe_adb_mhd(const Real d, const Real v1, const Real v2,
  const Real v3, const Real h, const Real b1, const Real b2, const Real b3,
  const Real x, const Real y, Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7]);

#ifdef RADIATION_HYDRO

void esys_roe_rad_hyd(const Real v1, const Real v2, const Real v3, const Real h, const Real dt,
  const Real pressure, const Real Er, const Real Fr[3], const Real Edd[6], const Real Sigma_a, const Real Sigma_t, const int flag,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5]);

void hlle_thick(const Cons1DS Ul, const Cons1DS Ur,
                   const Prim1DS Wl, const Prim1DS Wr,
                   const Real Bxi, Cons1DS *pFlux);
#endif

#ifdef RADIATION_MHD
void esys_roe_rad_mhd(const Real d, const Real v1, const Real v2, const Real v3, const Real dt, 
  const Real pressure, const Real Er, const Real Fr[3], const Real Edd[6], const Real Sigma_a, const Real Sigma_t, const int flag,
  const Real h, const Real b1, const Real b2, const Real b3, 
  Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7]);
#endif

/* All of the Riemann solvers in this directory contain the same function name
 */
void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pF);



#ifdef SPECIAL_RELATIVITY
void entropy_flux (const Cons1DS Ul, const Cons1DS Ur,
		   const Prim1DS Wl, const Prim1DS Wr,
		   const Real Bx, Real *pFlux);
#endif

#endif /* RSOLVERS_PROTOTYPES_H */
