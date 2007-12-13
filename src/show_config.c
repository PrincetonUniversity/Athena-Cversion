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
#endif

  nscal = NSCALARS;
  ath_pout(0," Passive scalars:         %d\n",nscal);

#if defined(SELF_GRAVITY_USING_MULTIGRID)
  ath_pout(0," Self-gravity:            using multigrid\n");
#elif defined(SELF_GRAVITY_USING_FFT)
  ath_pout(0," Self-gravity:            using FFTs\n");
#else
  ath_pout(0," Self-gravity:            none\n");
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

  ath_pout(0," Flux:                    %s\n",FLUX_TYPE);
  ath_pout(0," Unsplit 3D integrator:   %s\n",UNSPLIT_INTEGRATOR);

#if defined(SINGLE_PREC)
  ath_pout(0," Precision:               SINGLE_PREC\n");
#elif defined(DOUBLE_PREC)
  ath_pout(0," Precision:               DOUBLE_PREC\n");
#endif

  ath_pout(0," Output Modes:\n");
#ifdef WRITE_GHOST_CELLS
  ath_pout(0,"   Ghost Cells:           enabled\n");
#else
  ath_pout(0,"   Ghost Cells:           disabled\n");
#endif

  ath_pout(0," Parallel Modes:\n");
#if defined(MPI_SERIAL)
  ath_pout(0,"   MPI:                   MPI_SERIAL\n");
#elif defined(MPI_PARALLEL)
  ath_pout(0,"   MPI:                   MPI_PARALLEL\n");
#else
  ath_pout(0,"   MPI:                   undefined\n");
#endif

#ifdef H_CORRECTION
  ath_pout(0," H-correction:            enabled\n");
#else
  ath_pout(0," H-correction:            disabled\n");
#endif

#ifdef FFT_ENABLED
  ath_pout(0," FFT:                     enabled\n");
#else
  ath_pout(0," FFT:                     disabled\n");
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
#endif

  par_seti("configure","nscalars","%d",NSCALARS,"Number of passive scalars");

#if defined(SELF_GRAVITY_USING_MULTIGRID)
  par_sets("configure","self-gravity","multigrid","Self-gravity algorithm");
#elif defined(SELF_GRAVITY_USING_FFT)
  par_sets("configure","self-gravity","FFT","Self-gravity algorithm");
#else
  par_sets("configure","self-gravity","none","Self-gravity algorithm");
#endif

#if defined(FIRST_ORDER)
  par_seti("configure","order","%d",1,"Order of accuracy");
#elif defined(SECOND_ORDER)
  par_seti("configure","order","%d",2,"Order of accuracy");
#elif defined(THIRD_ORDER)
  par_seti("configure","order","%d",3,"Order of accuracy");
#endif

  par_sets("configure","flux",FLUX_TYPE,"Flux function");
  par_sets("configure","integrator",UNSPLIT_INTEGRATOR,"Unsplit 3D integrator");

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

#if defined(MPI_SERIAL)
  par_sets("configure","mpi","no","Is code serial or MPI parallel enabled?");
#elif defined(MPI_PARALLEL)
  par_sets("configure","mpi","yes","Is code serial or MPI parallel enabled?");
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

  return;
}
