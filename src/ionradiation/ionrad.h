#ifndef IONRAD_H
#define IONRAD_H  

#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

#ifdef ION_RADIATION
/* This hearder includes only things used for both point or planar
   radiation */

#ifdef DOUBLE_PREC
#  define EPS   DBL_EPSILON          /* Arithmetic precision */
#  define LARGE DBL_MAX              /* A big number */
#  define SMALL DBL_MIN              /* A small number */
#  define MP_RL MPI_DOUBLE
#else
#  define EPS   FLT_EPSILON
#  define LARGE FLT_MAX
#  define SMALL FLT_MIN
#  define MP_RL MPI_FLOAT
#endif

#define MINFLUXFRAC 1.0e-3         /* Smallest flux fraction to keep. */
#define COOLFRAC 1.0e-2            /* Level of neutrality or ionization
				      required to allow molecular cooling */
#define MINOPTDEPTH 1.0e-6         /* Minimum allowed optical depth per cell
				      to ionizing photons. */
#define IONFRACFLOOR 1.0e-6        /* Minimum ionization fraction. Ion
				      fractions are floored only if
				      the fraction is < IONFRACFLOOR
				      and the optical depth per cell
				      is below MINOPTDEPTH */
#define CION 8.0e5                 /* Approximate sound speed in fully
				      ionized gas */
#define MAXCELLCOUNT 20            /* Number of cells required to exceed
				      threshold for maximum change to
				      trigger new iteration. */


/* ------------------------------------------------------------
 * Global variable definition block
 * ------------------------------------------------------------
 *
 * Here we allocate storage for the global variables that will
 * be shared by all radiation routines. The first part of
 * this section is included only in ionrad.c, and it defines the
 * variables. The second part is included in all the ionizing
 * radiation modules, and provides data sharing between them.
 */

#ifdef IONRAD_C
/* Physical constants and parameters read from inputs file */
Real sigma_ph;                     /* Photoionization cross section */
Real m_H;                          /* Mean mass per H nucleus */
Real mu;                           /* Mean particle mass in neutral
				      gas */
Real alpha_C;                      /* Carbon mass fraction */
Real e_gamma;                      /* Energy added to gas by a single
				      photoionization */
Real k_B;                          /* Boltzmann's constant */
Real time_unit;                    /* Number of seconds in one unit of
				      code time */
Real max_de_iter;                  /* Maximum change in total gas
				      energy per iteration */
Real max_de_therm_iter;            /* Maximum change in gas thermal
				      energy per iteration */
Real max_dx_iter;                  /* Maximum change in ion
				      fraction per iteration */
Real max_de_step;                  /* Maximum change in total gas
				      energy per hydro step */
Real max_de_therm_step;            /* Maximum change in gas thermal
				      energy per hydro step */
Real max_dx_step;                  /* Maximum change in ion
				      fraction per hydro step */
Real tfloor;                       /* Temperature floor -- needed because
				      hydro / mhd doesn't guarantee that
				      the temperature is positive */
Real tceil;       		   /* Temperature max -- needed because hydro
				      / mhd can produce unphysically large
				      temperatures */
int maxiter;                       /* Maximum number of sub-cycle
				      iterations allowed */

/* Global grid information */
Real min_area;                     /* Smallest cell face area */
Real d_nlo;                        /* "Low" neutral density, defined
				      as the value that gives an
				      optical depth of MINOPTDEPTH */
#ifdef MPI_PARALLEL
Domain *pD;
int NGrid_x1, NGrid_x2, NGrid_x3;
#endif

#if defined(MPI_PARALLEL) && defined(ION_RADPOINT)
/* MPI buffer variables, used for point source radiation only.
   A few notes on the buffer: 
   1. The buffer is attached and detached at the start of every call
   to ion_radtransfer_3d. This costs some extra overhead, but ensures
   that the buffer we use here stays compartmentalized and doesn't
   affect communication in any other intergrator.
   2. To minimize overhead, we include the definition here and do the
   attachment and detachment inside ionrad_3d.c rather than inside the
   radpoint routines. That way we can attach and detach once per call
   to radiation, rather than once per subcycle.
   3. If we run out of buffer, then we have to wait for all pending
   communication to complete, detach the buffer, and reattach a larger
   one. Since this is costly, we should make the buffer as large as
   possible without hogging too many resources.
*/
#define INIT_MPI_BUFSIZE 32768
int mpi_bufsize = 0;
void *mpi_buffer;
#endif /* MPI_PARALLEL and ION_RADPOINT */

#else /* IONRAD_C */

extern Real sigma_ph;              /* Photoionization cross section */
extern Real m_H;                   /* Mean mass per H nucleus */
extern Real mu;                    /* Mean particle mass in neutral
				      gas */
extern Real alpha_C;               /* Carbon mass fraction */
extern Real e_gamma;               /* Energy added to gas by a single
				      photoionization */
extern Real k_B;                   /* Boltzmann's constant */
extern Real time_unit;             /* Number of seconds in one unit of
				      code time */
extern Real max_de_iter;           /* Maximum change in total gas
				      energy per iteration */
extern Real max_de_therm_iter;     /* Maximum change in gas thermal
				      energy per iteration */
extern Real max_dx_iter;           /* Maximum change in ion
				      fraction per iteration */
extern Real max_de_step;           /* Maximum change in total gas
				      energy per hydro step */
extern Real max_de_therm_step;     /* Maximum change in gas thermal
				      energy per hydro step */
extern Real max_dx_step;           /* Maximum change in ion
				      fraction per hydro step */
extern Real tfloor;                /* Temperature floor -- needed because
				      hydro / mhd doesn't guarantee that
				      the temperature is positive */
extern Real tceil;		   /* Temperature max -- needed because hydro
				      / mhd can produce unphysically large
				      temperatures */
extern int maxiter;                /* Maximum number of sub-cycle
				      iterations allowed */
extern Real min_area;              /* Smallest cell face area */
extern Real d_nlo;                 /* "Low" neutral density, defined
				      as the value that gives an
				      optical depth of MINOPTDEPTH */
#ifdef MPI_PARALLEL
extern Domain *pD;
extern int NGrid_x1, NGrid_x2, NGrid_x3;
#endif
#if defined(MPI_PARALLEL) && defined(ION_RADPOINT)
#define INIT_MPI_BUFSIZE 32768
extern int mpi_bufsize ;
extern void *mpi_buffer;
#endif /* MPI_PARALLEL and ION_RADPOINT */

#endif /* IONRAD_C */

#endif /* ION_RADIATION */

#endif /* IONRAD_H */
