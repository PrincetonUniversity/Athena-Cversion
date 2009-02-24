#include "copyright.h"
/*==============================================================================
 * FILE: show_config.c
 *
 * PURPOSE: Outputs information on configuration of Athena.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   show_config     - prints diagnostic message showinf code configuration 
 *   show_config_par - adds configuration information to database used by par
 *============================================================================*/

#include <stdio.h>
#include "defs.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* show_config:  the packages and features reported on by this functin should
 *   be kept consistent with the optional packages and features added by the
 *   file configure.ac in the top-level directory   */

void show_config(void)
{
  int nscal;
  ath_pout(0,"\nConfiguration details:\n\n");
  ath_pout(0," Problem:                 %s\n",A_PROBLEM);

#if defined(HYDRO)
  ath_pout(0," Gas properties:          HYDRO\n");
#elif defined(MHD)
  ath_pout(0," Gas properties:          MHD\n");
#endif

#if defined(ADIABATIC)
  ath_pout(0," Equation of State:       ADIABATIC\n");
#elif defined(ISOTHERMAL)
  ath_pout(0," Equation of State:       ISOTHERMAL\n");
#else
  ath_pout(0," Equation of State:       " EOS_STR "\n");
#endif

  nscal = NSCALARS;
  ath_pout(0," Passive scalars:         %d\n",nscal);

#if defined(SELF_GRAVITY_USING_MULTIGRID)
  ath_pout(0," Self-gravity:            using multigrid\n");
#elif defined(SELF_GRAVITY_USING_FFT)
  ath_pout(0," Self-gravity:            using FFTs\n");
#else
  ath_pout(0," Self-gravity:            OFF\n");
#endif

#if defined(OHMIC)
  ath_pout(0," Resistivity:             Ohmic\n");
#elif defined(HALL_MHD)
  ath_pout(0," Resistivity:             Hall\n");
#else
  ath_pout(0," Ohmic resistivity:       OFF\n");
#endif

#if defined(NAVIER_STOKES)
  ath_pout(0," Viscosity:               Navier-Stokes\n");
#elif defined(BRAGINSKII)
  ath_pout(0," Viscosity:               Braginskii\n");
#else
  ath_pout(0," Viscosity:               OFF\n");
#endif

#if defined(ISOTROPIC_CONDUCTION)
  ath_pout(0," Thermal conduction:      isotropic\n");
#elif defined(ANISOTROPIC_CONDUCTION)
  ath_pout(0," Thermal conduction:      anisotropic\n");
#else
  ath_pout(0," Thermal conduction:      OFF\n");
#endif

#ifdef ION_RADIATION
  ath_pout(0," Ionizing radiation:      ON\n");
#else
  ath_pout(0," Ionizing radiation:      OFF\n");
#endif

#ifdef ION_RADPOINT
  ath_pout(0,"   Ionizing point sources:  ON\n");
#else
  ath_pout(0,"   Ionizing point sources:  OFF\n");
#endif

#ifdef ION_RADPLANE
  ath_pout(0,"   Ionizing plane sources:  ON\n");
#else
  ath_pout(0,"   Ionizing plane sources:  OFF\n");
#endif

#if defined(FIRST_ORDER)
  ath_pout(0," Order of Accuracy:       1 (FIRST_ORDER)\n");
#elif defined(SECOND_ORDER)
  ath_pout(0," Order of Accuracy:       2 (SECOND_ORDER)\n");
#elif defined(THIRD_ORDER)
  ath_pout(0," Order of Accuracy:       3 (THIRD_ORDER)\n");
#elif defined(THIRD_ORDER_EXTREMA_PRESERVING)
  ath_pout(0," Order of Accuracy:       3e (THIRD_ORDER_EXTREMA_PRESERVING)\n");
#endif

#if defined(ROE_FLUX)
  ath_pout(0," Flux:                    roe\n");
#elif defined(HLLE_FLUX)
  ath_pout(0," Flux:                    hlle\n");
#elif defined(HLLD_FLUX)
  ath_pout(0," Flux:                    hlld\n");
#elif defined(HLLC_FLUX)
  ath_pout(0," Flux:                    hllc\n");
#elif defined(FORCE_FLUX)
  ath_pout(0," Flux:                    force\n");
#elif defined(EXACT_FLUX)
  ath_pout(0," Flux:                    exact\n");
#elif defined(TWO_SHOCK_FLUX)
  ath_pout(0," Flux:                    two-shock\n");
#endif

#if defined(CTU_INTEGRATOR)
  ath_pout(0," Unsplit integrator:      ctu\n");
#elif defined(VL_INTEGRATOR)
  ath_pout(0," Unsplit integrator:      vl\n");
#endif

#if defined(SINGLE_PREC)
  ath_pout(0," Precision:               SINGLE_PREC\n");
#elif defined(DOUBLE_PREC)
  ath_pout(0," Precision:               DOUBLE_PREC\n");
#endif

#ifdef WRITE_GHOST_CELLS
  ath_pout(0," Ghost cell Output:       ON\n");
#else
  ath_pout(0," Ghost cell Output:       OFF\n");
#endif

#if defined(MPI_PARALLEL)
  ath_pout(0," Parallel Modes: MPI:     ON\n");
#else
  ath_pout(0," Parallel Modes: MPI:     OFF\n");
#endif

#ifdef H_CORRECTION
  ath_pout(0," H-correction:            ON\n");
#else
  ath_pout(0," H-correction:            OFF\n");
#endif

#ifdef FFT_ENABLED
  ath_pout(0," FFT:                     ON\n");
#else
  ath_pout(0," FFT:                     OFF\n");
#endif

#ifdef SHEARING_BOX
  ath_pout(0," Shearing Box:            ON\n");
#else
  ath_pout(0," Shearing Box:            OFF\n");
#endif

#ifdef FARGO
  ath_pout(0," FARGO:                   ON\n");
#else
  ath_pout(0," FARGO:                   OFF\n");
#endif
}

/*----------------------------------------------------------------------------*/
/*  show_config_par:  Add the configure block to the parameter database used
 *    by the functions in par.c.  */

void show_config_par(void)
{
  par_sets("configure","problem",A_PROBLEM,"Name of the problem file");

#if defined(HYDRO)
  par_sets("configure","gas","hydro","Hydrodynamic gas");
#elif defined(MHD)
  par_sets("configure","gas","mhd","Magnetohydrodynamic gas");
#endif

#if defined(ADIABATIC)
  par_sets("configure","eq_state","adiabatic","Equation of state");
#elif defined(ISOTHERMAL)
  par_sets("configure","eq_state","isothermal","Equation of state");
#else
  par_sets("configure","eq_state",EOS_STR,"Equation of state");
#endif

  par_seti("configure","nscalars","%d",NSCALARS,"Number of passive scalars");

#if defined(SELF_GRAVITY_USING_MULTIGRID)
  par_sets("configure","self-gravity","multigrid","Self-gravity algorithm");
#elif defined(SELF_GRAVITY_USING_FFT)
  par_sets("configure","self-gravity","FFT","Self-gravity algorithm");
#else
  par_sets("configure","self-gravity","OFF","Self-gravity algorithm");
#endif

#if defined(OHMIC)
  par_sets("configure","resistivity","Ohmic","resistivity algorithm");
#elif defined(HALL_MHD)
  par_sets("configure","resistivity","Hall","resistivity algorithm");
#else
  par_sets("configure","resistivity","OFF","resistivity algorithm");
#endif

#if defined(NAVIER_STOKES)
  par_sets("configure","viscosity","Navier-Stokes","viscosity algorithm");
#elif defined(BRAGINSKII)
  par_sets("configure","viscosity","Braginskii","viscosity algorithm");
#else
  par_sets("configure","viscosity","none","viscosity algorithm");
#endif

#if defined(ISOTROPIC_CONDUCTION)
  par_sets("configure","thermal conduction","isotropic","conduction algorithm");
#elif defined(ANISOTROPIC_CONDUCTION)
  par_sets("configure","thermal conduction","anisotropic","conduction algorithm");
#else
  par_sets("configure","thermal conduction","none","conduction algorithm");
#endif

#if defined(ION_RADIATION)
  par_sets("configure","ionizing radiation","yes","Ionizing rad transfer?");
#else
  par_sets("configure","ionizing radiation","no","Ionizing rad transfer?");
#endif

#if defined(ION_RADPOINT)
  par_sets("configure","point source","yes","Point source of radiation?");
#else
  par_sets("configure","point source","no","Point source of radiation?");
#endif

#if defined(ION_RADPLANE)
  par_sets("configure","plane source","yes","Plane source of radiation?");
#else
  par_sets("configure","plane source","no","Plane source of radiation?");
#endif

#if defined(FIRST_ORDER)
  par_seti("configure","order","%d",1,"Order of accuracy");
#elif defined(SECOND_ORDER)
  par_seti("configure","order","%d",2,"Order of accuracy");
#elif defined(THIRD_ORDER)
  par_seti("configure","order","%d",3,"Order of accuracy");
#endif

#if defined(ROE_FLUX)
  par_sets("configure","flux","roe","Flux function");
#elif defined(HLLE_FLUX)
  par_sets("configure","flux","hlle","Flux function");
#elif defined(HLLD_FLUX)
  par_sets("configure","flux","hlld","Flux function");
#elif defined(HLLC_FLUX)
  par_sets("configure","flux","hllc","Flux function");
#elif defined(FORCE_FLUX)
  par_sets("configure","flux","force","Flux function");
#elif defined(EXACT_FLUX)
  par_sets("configure","flux","exact","Flux function");
#elif defined(TWO_SHOCK_FLUX)
  par_sets("configure","flux","two-shock","Flux function");
#endif

#if defined(CTU_INTEGRATOR)
  par_sets("configure","integrator","ctu","Unsplit integrator");
#elif defined(VL_INTEGRATOR)
  par_sets("configure","integrator","vl","Unsplit integrator");
#endif

#if defined(SINGLE_PREC)
  par_sets("configure","precision","single","Type of Real variables");
#elif defined(DOUBLE_PREC)
  par_sets("configure","precision","double","Type of Real variables");
#endif

#ifdef WRITE_GHOST_CELLS
  par_sets("configure","write_ghost","yes","Ghost cells included in output?");
#else
  par_sets("configure","write_ghost","no","Ghost cells included in output?");
#endif

#if defined(MPI_PARALLEL)
  par_sets("configure","mpi","yes","Is code MPI parallel enabled?");
#else
  par_sets("configure","mpi","no","Is code MPI parallel enabled?");
#endif

#ifdef H_CORRECTION
  par_sets("configure","H-correction","yes","H-correction enabled?");
#else
  par_sets("configure","H-correction","no","H-correction enabled?");
#endif

#ifdef FFT_ENABLED
  par_sets("configure","FFT","yes","FFT enabled?");
#else
  par_sets("configure","FFT","no","FFT enabled?");
#endif

#ifdef SHEARING_BOX
  par_sets("configure","ShearingBox","yes","Shearing box enabled?");
#else
  par_sets("configure","ShearingBox","no","Shearing box enabled?");
#endif

#ifdef FARGO
  par_sets("configure","FARGO","yes","FARGO enabled?");
#else
  par_sets("configure","FARGO","no","FARGO enabled?");
#endif

  return;
}
