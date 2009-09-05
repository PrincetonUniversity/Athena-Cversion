#ifndef RECONSTRUCTION_PROTOTYPES_H
#define RECONSTRUCTION_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the src/reconstruction
 *   directory.
 *============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/
/* esystem_prim.c */

#if defined(BAROTROPIC) && defined(HYDRO)
void esys_prim_iso_hyd(const Real d, const Real v1,
  Real eigenvalues[],
  Real right_eigenmatrix[][4], Real left_eigenmatrix[][4]);
#endif

#if !defined(BAROTROPIC) && defined(HYDRO)
void esys_prim_adb_hyd(const Real d, const Real v1, const Real p,
  Real eigenvalues[],
  Real right_eigenmatrix[][5], Real left_eigenmatrix[][5]);
#endif

#if defined(BAROTROPIC) && defined(MHD)
void esys_prim_iso_mhd(const Real d, const Real v1, const Real b1,
  const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][6], Real left_eigenmatrix[][6]);
#endif

#if !defined(BAROTROPIC) && defined(MHD)
void esys_prim_adb_mhd(const Real d, const Real v1, const Real rho_a2,
  const Real b1, const Real b2, const Real b3, Real eigenvalues[],
  Real right_eigenmatrix[][7], Real left_eigenmatrix[][7]);
#endif

/*----------------------------------------------------------------------------*/
/*  All of the lr_states_*.c files in this directory contain the same function
 *  names below */

void lr_states_destruct(void);
void lr_states_init(int nx1, int nx2, int nx3);
void lr_states(const Prim1D W[], MHDARG( const Real Bxc[] , )
               const Real dtodx, const int is, const int ie,
               Prim1D Wl[], Prim1D Wr[]);
#ifdef CYLINDRICAL
void lr_states_cyl(const Prim1D W[], MHDARG( const Real Bxc[] , )
                   const Real dtodx, const int is, const int ie,
                   Prim1D Wl[], Prim1D Wr[]);
#endif /* CYLINDRICAL */

#endif /* RECONSTRUCTION_PROTOTYPES_H */
