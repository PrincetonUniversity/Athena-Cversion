#include "../copyright.h"
/*=============================================================================
FILE: init_particle.c
PURPOSE: Initialize particle related structures and functions. Particle
  integrator is enrolled by calling integrate_particle_init(int type).
  In init_particle(Grid *pG, Domain *pD), all the particle related arrays
  are allocated, interpolation and stopping time calculation functions are
  enrolled, most global variables for particles are evaluated. Also contains
  functions for particle array reallocation and destruction.

CONTAINS PUBLIC FUNCTIONS:
  VGFun_t integrate_particle_init(int type);
  void init_particle(Grid *pG, Domain *pD);
  void particle_destruct(Grid *pG);
  void particle_realloc(Grid *pG, long n);

History:
  Written by  Xuening Bai      Apr. 2009

==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

#ifdef PARTICLES         /* endif at the end of the file */

/*=========================== PROTOTYPES OF PRIVATE FUNCTIONS ===============================*/

void grid_limit(Grid *pG, Domain *pD);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/* Initialization for particles.
 * We assume to have "partypes" types of particles, each type has "parnum" particles.
 * We enforce that each type has equal number of particles to ensure equal resolution.
 * Allocate memory for the gas velocity/sound speed array, feedback array.
 */
void init_particle(Grid *pG, Domain *pD)
{
  int i, N1T, N2T, N3T, integratortype, interp, tsmode;
  Grain *GrArray;
  long size = 1000, size1 = 1, size2 = 1;

  /* get coordinate limit */
  grid_limit(pG, pD);
  N1T = iup-ilp+1;
  N2T = jup-jlp+1;
  N3T = kup-klp+1;

  /* check particle types */
  pG->partypes = par_geti_def("particle","partypes",1);

  if (pG->partypes < 0)
    ath_error("[init_particle]: Particle types must not be negative!\n");

  /* initialize the particle array */
  if(par_exist("particle","parnumcell"))
  {
    /* if we consider number of particles per cell */
    size1 = N1T*N2T*N3T*(long)(pG->partypes*par_geti("particle","parnumcell"));
    if (size1 < 0)
      ath_error("[init_particle]: Particle number must not be negative!\n");
  }

  if(par_exist("particle","parnumgrid"))
  {
    /* if we consider number of particles per grid */
    size2 = (long)(pG->partypes*par_geti("particle","parnumgrid"));
    if (size2 < 0)
      ath_error("[init_particle]: Particle number must not be negative!\n");
    /* account for the ghost cells */
    size2 = (long)(size2/((double)(pG->Nx1*pG->Nx2*pG->Nx3))*N1T*N2T*N3T);
  }

  size = MAX(size, MAX(size1, size2));
  pG->arrsize = (long)(1.2*size);	/* account for number fluctuations */

  pG->particle = (Grain*)calloc_1d_array(pG->arrsize, sizeof(Grain));
  if (pG->particle == NULL) goto on_error;

  /* allocate memory for particle properties */
  pG->grproperty = (Grain_Property*)calloc_1d_array(pG->partypes, sizeof(Grain_Property));
  if (pG->grproperty == NULL) goto on_error;

  grrhoa = (Real*)calloc_1d_array(pG->partypes, sizeof(Real));
  if (grrhoa == NULL) goto on_error;

  tstop0 = (Real*)calloc_1d_array(pG->partypes, sizeof(Real));
  if (tstop0 == NULL) goto on_error;

  /* by default these global values are zero */
  for (i=0; i<pG->partypes; i++) {
    tstop0[i] = 0.0;
    grrhoa[i] = 0.0;
  }
  alamcoeff = 0.0;

  /* Set particle integrator to according to particle types */
  for (i=0; i<pG->partypes; i++)
    pG->grproperty[i].integrator = par_geti_def("particle","integrator",2);

  /* set the interpolation function pointer */
  interp = par_geti_def("particle","interp",2);
  if (interp == 1)
  { /* linear interpolation */
    getweight = getwei_linear;
    ncell = 2;
  }
  else if (interp == 2)
  { /* TSC interpolation */
    getweight = getwei_TSC;
    ncell = 3;
  }
  else if (interp == 3)
  { /* Quadratic polynomial interpolation */
    getweight = getwei_QP;
    ncell = 3;
  }
  else
    ath_error("[init_particle]: Invalid interp value (should equals 1 or 2)!\n");

  /* set the stopping time function pointer */
  tsmode = par_geti("particle","tsmode");
  if (tsmode == 1)
    get_ts = get_ts_general;
  else if (tsmode == 2)
    get_ts = get_ts_epstein;
  else if (tsmode == 3)
    get_ts = get_ts_fixed;
  else
    ath_error("[init_particle]: tsmode must be 1, 2 or 3!\n");

  /* allocate the memory for gas and feedback arrays */
  grid_d = (Real***)calloc_3d_array(N3T, N2T, N1T, sizeof(Real));
  if (grid_d == NULL) goto on_error;

  grid_v = (Vector***)calloc_3d_array(N3T, N2T, N1T, sizeof(Vector));
  if (grid_v == NULL) goto on_error;

#ifndef ISOTHERMAL
  grid_cs = (Real***)calloc_3d_array(N3T, N2T, N1T, sizeof(Real));
  if (grid_cs == NULL) goto on_error;
#endif

#ifdef FEEDBACK
  pG->feedback = (Vector***)calloc_3d_array(N3T, N2T, N1T, sizeof(Vector));
  if (pG->feedback == NULL) goto on_error;

  pG->Eloss = (Real***)calloc_3d_array(N3T, N2T, N1T, sizeof(Real));
  if (pG->Eloss == NULL) goto on_error;
#endif

  return;

  on_error:
    ath_error("[init_particle]: Error allocating memory.\n");
}


/* Finalization for particles */
void particle_destruct(Grid *pG)
{
  free_1d_array(pG->particle);

  free_1d_array(pG->grproperty);
  free_1d_array(grrhoa);

  /* free memory for gas and feedback arrays */
  if (grid_d != NULL) free_3d_array(grid_d);
  if (grid_v != NULL) free_3d_array(grid_v);

#ifndef ISOTHERMAL
  if (grid_cs != NULL) free_3d_array(grid_cs);
#endif

#ifdef FEEDBACK
  if (pG->feedback != NULL) free_3d_array(pG->feedback);
  if (pG->Eloss != NULL) free_3d_array(pG->Eloss);
#endif

  return;
}


/* Enlarge the particle array */
void particle_realloc(Grid *pG, long n)
{
  pG->arrsize = MAX((long)(1.2*pG->arrsize), n);

  if ((pG->particle = (Grain*)realloc(pG->particle, pG->arrsize*sizeof(Grain))) == NULL)
    ath_error("[init_particle]: Error re-allocating memory with array size %ld.\n", n);

  return;
}

/*=========================== PRIVATE FUNCTIONS ===============================*/

/* Calculate the left and right grid limit
   Input: pG: grid;
   Output: ilp,iup,jlp,jup,klp,kup: grid limit indices;
           x1lpar,x1upar,x2lpar,x2upar,x3lpar,x3upar: grid boundary coordinates
*/
void grid_limit(Grid *pG, Domain *pD)
{
  int m1, m2, m3;	/* dimension flags */
  int my_iproc, my_jproc, my_kproc;

  if (pG->Nx1 > 1) m1 = 1;
  else m1 = 0;

  if (pG->Nx2 > 1) m2 = 1;
  else m2 = 0;

  if (pG->Nx3 > 1) m3 = 1;
  else m3 = 0;

  /* set left and right grid indices */
  ilp = pG->is - m1*nghost;
  iup = pG->ie + m1*nghost;

  jlp = pG->js - m2*nghost;
  jup = pG->je + m2*nghost;

  klp = pG->ks - m3*nghost;
  kup = pG->ke + m3*nghost;

  /* set left and right boundary for removing particles */
  /* Note: for outflow B.C. (ibc=2), we only remove the particles in
   * the outermost layer of the ghost cells */
  get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

  if ((par_geti_def("grid","ibc_x1",4) == 2) && (my_iproc == 0))
    x1lpar = pG->x1_0 + (ilp+m1 + pG->idisp)*pG->dx1;
  else
    x1lpar = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;

  if ((par_geti_def("grid","obc_x1",4) == 2) && (my_iproc == pD->NGrid_x1-1))
    x1upar = pG->x1_0 + (iup + pG->idisp)*pG->dx1;
  else
    x1upar = pG->x1_0 + (pG->ie + 1 + pG->idisp)*pG->dx1;

  if ((par_geti_def("grid","ibc_x2",4) == 2) && (my_jproc == 0))
    x2lpar = pG->x2_0 + (jlp+m2 + pG->jdisp)*pG->dx2;
  else
    x2lpar = pG->x2_0 + (pG->js + pG->jdisp)*pG->dx2;

  if ((par_geti_def("grid","obc_x2",4) == 2) && (my_jproc == pD->NGrid_x2-1))
    x2upar = pG->x2_0 + (jup + pG->jdisp)*pG->dx2;
  else
    x2upar = pG->x2_0 + (pG->je + 1 + pG->jdisp)*pG->dx2;

  if ((par_geti_def("grid","ibc_x3",4) == 2) && (my_kproc == 0))
    x3lpar = pG->x3_0 + (klp+m3 + pG->kdisp)*pG->dx3;
  else
    x3lpar = pG->x3_0 + (pG->ks + pG->kdisp)*pG->dx3;

  if ((par_geti_def("grid","obc_x3",4) == 2) && (my_kproc == pD->NGrid_x3-1))
    x3upar = pG->x3_0 + (kup + pG->kdisp)*pG->dx3;
  else
    x3upar = pG->x3_0 + (pG->ke + 1 + pG->kdisp)*pG->dx3;

  if (pG->Nx1 == 1) x1upar += MAX(0.1*fabs(pG->x1_0), 1.);
  if (pG->Nx2 == 1) x2upar += MAX(0.1*fabs(pG->x2_0), 1.);
  if (pG->Nx3 == 1) x3upar += MAX(0.1*fabs(pG->x3_0), 1.);

  return;
}

#endif /*PARTICLES*/
