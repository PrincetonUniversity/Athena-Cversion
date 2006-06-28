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
  fprintf(stderr,"\nConfiguration details:\n\n");
  fprintf(stderr," Problem:                 %s\n",A_PROBLEM);

#if defined(HYDRO)
  fprintf(stderr," Gas properties:          HYDRO\n");
#elif defined(MHD)
  fprintf(stderr," Gas properties:          MHD\n");
#endif

#if defined(ADIABATIC)
  fprintf(stderr," Equation of State:       ADIABATIC\n");
#elif defined(ISOTHERMAL)
  fprintf(stderr," Equation of State:       ISOTHERMAL\n");
#endif

#if defined(FIRST_ORDER)
  fprintf(stderr," Order of Accuracy:       1 (FIRST_ORDER)\n");
#elif defined(SECOND_ORDER)
  fprintf(stderr," Order of Accuracy:       2 (SECOND_ORDER)\n");
#elif defined(THIRD_ORDER)
  fprintf(stderr," Order of Accuracy:       3 (THIRD_ORDER)\n");
#endif

  fprintf(stderr," Flux:                    %s\n",FLUX_NAME);
  fprintf(stderr," Unsplit 3D integrator:   %s\n",THREE_DIM_INT);

#if defined(SINGLE_PREC)
  fprintf(stderr," Precision:               SINGLE_PREC\n");
#elif defined(DOUBLE_PREC)
  fprintf(stderr," Precision:               DOUBLE_PREC\n");
#endif

  fprintf(stderr," Output Modes:\n");
#ifdef WRITE_GHOST_CELLS
  fprintf(stderr,"   Ghost Cells:           enabled\n");
#else
  fprintf(stderr,"   Ghost Cells:           disabled\n");
#endif

  fprintf(stderr," Parallel Modes:\n");
#if defined(MPI_SERIAL)
  fprintf(stderr,"   MPI:                   MPI_SERIAL\n");
#elif defined(MPI_PARALLEL)
  fprintf(stderr,"   MPI:                   MPI_PARALLEL\n");
#else
  fprintf(stderr,"   MPI:                   undefined\n");
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

#if defined(FIRST_ORDER)
  par_seti("configure","order","%d",1,"Order of accuracy");
#elif defined(SECOND_ORDER)
  par_seti("configure","order","%d",2,"Order of accuracy");
#elif defined(THIRD_ORDER)
  par_seti("configure","order","%d",3,"Order of accuracy");
#endif

  par_sets("configure","flux",FLUX_NAME,"Flux function");
  par_sets("configure","integrator",THREE_DIM_INT,"Unsplit 3D integrator");

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

  return;
}
