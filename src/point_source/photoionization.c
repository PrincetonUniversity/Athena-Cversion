#include "../copyright.h"
/*==============================================================================
 * FILE: photoionization.c
 *
 * PURPOSE: Contains functions for computing photoionization equilibrium.
 *          This version works with radiation fields from point sources, 
 *           which are computed using the adaptive ray tracing scheme
 *          implmented in point_source.c.  Both implementations follow
 *          Krumholz et al. 2007 and are reimplemented in modern Athena
 *          by S. Davis.          
 *          
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   photoionization()        --  compute ionization rates
 *   init_photoionization()
 *   new_dt_phot()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

//#define IONIZATION_ONLY

#ifdef PHOTOIONIZATION

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

#define MAXSIGNCOUNT 4
#define DAMPFACTOR 0.5

/* Units and constants */
#define EV  1.602e-12
#define KB  1.38e-16
#define IHI (13.6*EV)    /* H ionization potential */
#define GAMMAKI 2.0e-26  /* Koyama & Inutsuka (2002) heating rate, in erg/s */

/* MacDonald & Bailey cooling function table */
static Real xmat[] = {-0.133, 0.105, 0.452, 0.715, 0.901
		      , 1.030, 1.082, 1.174, 1.257, 1.362
		      , 1.448, 1.523, 1.569, 1.582, 1.539
		      , 1.430, 1.275, 1.168, 1.092, 1.019
		      , 1.000, 1.004, 1.008, 0.987, 0.905
		      , 0.738, 0.603, 0.555, 0.552, 0.554
		      , 0.552, 0.535, 0.425, 0.275, 0.251
		      , 0.232, 0.247, 0.283, 0.322, 0.363
		      , 0.397};

static Real ***edot;               /* Rate of change of energy */
static Real ***nHdot;              /* Rate of change of neutral
				      density */
static int  ***last_sign;          /* Last sign of nHdot -- keep track
				      of this to avoid oscillatory
				      overstability */
static int  ***sign_count;         /* Number of successive times the
				      sign of nHdot has flipped -- use
				      this to avoid oscillatory
				      overstability */
static Real ***e_init;             /* Total energies on entry to
				      routine */
static Real ***e_th_init;          /* Thermal energies on entry to
				      routine */

/* Physical constants and parameters read from inputs file */
static Real mH;                   /* Mean mass per H nucleus */
static Real mu;                    /* Mean particle mass in neutral
				      gas */
static Real alphaC;               /* Carbon mass fraction */
static Real egamma;               /* Energy added to gas by a single
				      photoionization */
static Real kB;                   /* Boltzmann's constant */
static Real time_unit;             /* Number of seconds in one unit of
				      code time */
static Real enorm;
static Real max_de;                /* Maximum change in total gas
				      energy per time step */
static Real max_de_therm;          /* Maximum change in gas thermal
				      energy per time step */
static Real tfloor;                /* Temperature floor -- needed because
				      hydro / mhd doesn't guarantee that
				      the temperature is positive */
static int maxiter;                /* Maximum number of sub-cycle
				      iterations allowed */


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * compute_rates() - compute thermal rates
 * apply_temp_floor()
 * save_energy()
 *============================================================================*/
Real compute_therm_rates(GridS *pGrid);
Real compute_chem_rates(GridS *pGrid);
void apply_temp_floor(GridS *pGrid);
void save_energy (GridS *pGrid);
void ionization_update(GridS *pGrid, Real dt);
Real recomb_rate_coef(Real T);
Real coll_ion_rate_coef(Real T);
Real recomb_cool_rate_coef(Real T);
Real dmc_cool_rate(Real x, Real T);
Real ki_cool_rate(Real T);
Real ki_heat_rate();

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/

void init_photoionization(MeshS *pM)
{
  GridS *pGrid;
  int i,j,k;

/* Read in parameters */
  max_de = par_getd_def("problem", "max_de",3.);
  max_de_therm = par_getd_def("problem", "max_de_therm",3.);
  tfloor = par_getd_def("problem", "tfloor",0.0);
  maxiter = par_geti_def("problem", "maxiter",10);
  mH = par_getd("problem","mH");
  mu = par_getd("problem", "mu");
  egamma = par_getd("problem", "egamma");
  alphaC = par_getd_def("problem", "alphaC",3e-3);
  kB = par_getd("problem", "kB");
  time_unit = par_getd("problem", "time_unit");
  enorm = par_getd("problem","enorm");
  
/* Get pointer to pGrid, assuming no mesh refinement */
  if (pM->Domain[0][0].Grid != NULL)
    pGrid = pM->Domain[0][0].Grid;
  else
    ath_error("[init_photoionization]: pM->Domain[0][0].Grid = NULL\n");

/* Allocate memory for rate arrays */

  edot = (Real***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(Real));
  nHdot = (Real***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(Real));
  last_sign = (int***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(int));
  sign_count = (int***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(int));
  e_init = (Real***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(Real));
  e_th_init = (Real***) 
    calloc_3d_array(pGrid->Nx[2], pGrid->Nx[1], pGrid->Nx[0], 
		    sizeof(Real));

  /* Offset pointers to account for ghost cells */
  
  edot -= pGrid->ks;
  nHdot -= pGrid->ks;
  last_sign -= pGrid->ks;
  sign_count -= pGrid->ks;
  e_init -= pGrid->ks;
  e_th_init -= pGrid->ks;
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    edot[k] -= pGrid->js;
    nHdot[k] -= pGrid->js;
    last_sign[k] -= pGrid->js;
    sign_count[k] -= pGrid->js;
    e_init[k] -= pGrid->js;
    e_th_init[k] -= pGrid->js;
    for (j=pGrid->js; j<=pGrid->je; j++) {
      edot[k][j] -= pGrid->is;
      nHdot[k][j] -= pGrid->is;
      last_sign[k][j] -= pGrid->is;
      sign_count[k][j] -= pGrid->is;
      e_init[k][j] -= pGrid->is;
      e_th_init[k][j] -= pGrid->is;
    }
  }

  /* Initialize counter arrays */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	last_sign[k][j][i] = 0;
	sign_count[k][j][i] = 0;
      }
    }
  }
}

void photoionization(MeshS *pM)
{

  GridS *pGrid;
  Real dt_chem, dt_therm, dt, dt_done;
  int n, niter, hydro_done, thermal_advance;
  int i,j,k;
  int nl,nd;

/* Get pointer to pGrid, assuming no mesh refinement */
  if (pM->Domain[0][0].Grid != NULL)
    pGrid = pM->Domain[0][0].Grid;
  else
    ath_error("[photoionization]: pM->Domain[0][0].Grid = NULL\n");

/* Trees are initialized in call of init_pointsource() in main.c 
 * this will need to be modified if/when point source origins move
 * between timesteps */

/* Set all temperatures below the floor to the floor */
  apply_temp_floor(pGrid);

/* Save total and thermal energies passed in -- we use these to
   compute time steps */
  save_energy(pGrid);

  /* Begin the radiation sub-cycle */
  dt_done = 0.0;
  hydro_done = 0;
  for (niter = 0; niter<maxiter; niter++) {
    //printf("niter: %d %g\n",niter,dt_done);
/* Initialize photoionization rate array */
    for (k=pGrid->ks; k<=pGrid->ke; k++) {
      for (j=pGrid->js; j<=pGrid->je; j++) {
	for (i=pGrid->is; i<=pGrid->ie; i++) {
	  pGrid->phrate[k][j][i] = 0.0;
	}}}
    
/* Compute photoionization rate with adaptive ray tracing */
    for (n=0; n<pGrid->nsource; n++)
      update_tree_radiation(n, pGrid);
    
/* Compute rates and time step for thermal energy update */
    dt_therm = compute_therm_rates(pGrid);

/* Compute rates and time step for chemistry update */
    dt_chem = compute_chem_rates(pGrid);

/* Set time step to smaller of thermal and chemical time steps. Also record 
 * if we are advancing on the chemical or thermal time step. */
    if (dt_therm < dt_chem) {
      dt = dt_therm;
      thermal_advance = 1;
    } else {
      dt = dt_chem;
      thermal_advance = 0;
    }

/* If necessary, scale back time step to avoid exceeding hydro time step. */
    if (dt_done + dt > pGrid->dt) {
      dt = pGrid->dt - dt_done;
      hydro_done = 1;
    }

/* update ionization fraction*/
    ionization_update(pGrid, dt);

    dt_done += dt;

/* Set all temperatures below the floor to the floor */
    apply_temp_floor(pGrid);

/* If we advanced the full hydro time step, exit loop. */
    if (hydro_done) 
      break;
    
/* Was the last update taken at the chemistry or the thermal time step? If it 
 * was taken at the thermal time step, then the energy must have changed as much
 * as was allowed. Scale back the hydro time step, and exit the loop. */
    if (thermal_advance) {
      pM->dt = dt_done;
      break;
    }
  }

/* If we exceeded the maximum number of iterations return to hydro with the time 
 * step re-set to what we managed to do. */
  if (niter==maxiter)
    pM->dt = dt_done;

  if (myID_Comm_world == 0) {
    if (thermal_advance && !hydro_done)
      printf("Radiation done in %d iterations, terminating on thermal timestep\n",
	     niter+1);
    else if (hydro_done)
      printf("Radiation done in %d iterations, terminating on hydro timestep\n",
	     niter+1);
    else
      printf("Radiation done in %d iterations, terminating on iteration number\n",
	     niter+1);
  }
/* Spread timestep across all Grid structures in all Domains */

  if(!hydro_done) {
    for (nl=0; nl<=(pM->NLevels)-1; nl++){
      for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
	if (pM->Domain[nl][nd].Grid != NULL) {
	  pM->Domain[nl][nd].Grid->dt = pM->dt;
	}
      }
    }
  }

}

/* Private functions */

Real compute_chem_rates(GridS *pGrid)
{
  int i, j, k, n;
  Real n_H, n_Hplus, n_e;
  Real e_sp, T, x;
  Real dt_chem, dt_chem_min;
#ifdef MPI_PARALLEL
  int err;
  Real dt_chem_min_glob;
#endif
  Real sw = 0.0;

  /* Initialize chemistry time step to large time step */
  dt_chem_min = HUGE_NUMBER;

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Get species abundances */
	n_H = pGrid->U[k][j][i].dn / mH;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].dn) / mH;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alphaC / (14.0 * mH);
	x = n_Hplus / (n_H + n_Hplus);

	/* Get gas temperature in K */
	e_sp = (pGrid->U[k][j][i].E -
		0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		       pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		       pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
		/ pGrid->U[k][j][i].d
#ifdef MHD
		- 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
			 pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
			 pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
		) / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*mH+(1.0-x)*mu)/ kB;
	if (T < tfloor) T = tfloor;

	/* Get rate of change of neutral density */
	nHdot[k][j][i] = 
	  recomb_rate_coef(T) * time_unit * n_e * n_Hplus
	  - pGrid->phrate[k][j][i] * n_H
	  - sw *coll_ion_rate_coef(T) * time_unit * n_e * n_H;

	/* Check if the sign has flipped -- oscillatory overstability
	   check */
	if (nHdot[k][j][i] < 0.0) {
	  if (last_sign[k][j][i] == 1) sign_count[k][j][i]++;
	  else if (sign_count[k][j][i] > 0) sign_count[k][j][i]--;
	  last_sign[k][j][i] = -1;
	} else if (nHdot[k][j][i] > 0.0) {
	  if (last_sign[k][j][i] == -1) sign_count[k][j][i]++;
	  else if (sign_count[k][j][i] > 0) sign_count[k][j][i]--;
	  last_sign[k][j][i] = 1;
	} else {
	  sign_count[k][j][i] = last_sign[k][j][i] = 0;
	}

	/* If sign has flipped too many times successively, this cell
	 * is probably experiencing oscillatory overstability. To
	 * combat this, decrease the update amount by a chosen factor
	 * for every repeat of the same sign past the trigger number
	 */
	for (n=MAXSIGNCOUNT; n<sign_count[k][j][i]; n++) {
	  edot[k][j][i] *= DAMPFACTOR;
	  nHdot[k][j][i] *= DAMPFACTOR;
	}

	/* Compute chemistry time step for this cell and find the min */
	dt_chem = 0.1 * n_e / fabs(nHdot[k][j][i]);
	if (dt_chem_min > dt_chem) dt_chem_min = dt_chem;
	dt_chem = 0.1 * n_H / fabs(nHdot[k][j][i]);
	if (dt_chem_min > dt_chem) dt_chem_min = dt_chem;

	//if (dt_chem < 1.0e3) {
	//  printf("Warning: dt_chem = %e, T = %e, x = %e, ijk = %d %d %d\n",
	//		 dt_chem, T, n_Hplus/(n_H+n_Hplus), i, j, k);
	  if (dt_chem <= 0.0) ath_error("[compute_chem_rates]: bad dt_chem\n");
	  //}
      }
    }
  }

#ifdef MPI_PARALLEL
  /* Sync chemistry timestep across processors */
  err = MPI_Allreduce(&dt_chem_min, &dt_chem_min_glob, 1, MP_RL, 
		      MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[compute_chem_rates]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt_chem_min = dt_chem_min_glob;
#endif

  return(dt_chem_min);
}
#undef MAXSIGNCOUNT
#undef DAMPFACTOR



Real compute_therm_rates(GridS *pGrid)
{
  int i, j, k;
  Real n_H, n_Hplus, n_e, e_thermal;
  Real e_sp, T, x, e_sp_min, e_th_min, e_min;
  Real dt_therm, dt_therm1, dt_therm2, dt_therm_min;
#ifdef MPI_PARALLEL
  int err;
  Real dt_therm_min_glob;
#endif
#if 0
  int isave, jsave, ksave;
#endif
  Real sw = 0.0;

  /* Initialize thermal time step to large value */
  dt_therm_min = HUGE_NUMBER;

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Get species abundances */
	n_H = pGrid->U[k][j][i].dn / mH;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].dn) / mH;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alphaC / (14.0 * mH);
	x = n_Hplus / (n_H + n_Hplus);

	/* Get gas temperature in K */
	e_thermal = pGrid->U[k][j][i].E - 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;
	e_sp = e_thermal / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*mH+(1.0-x)*mu)/ kB;
	//printf("%d %d %d %g\n",i,j,k,T);
	/* Check temperature floor. If cell is below temperature
	   floor, skip it. We'll fix it later */
	if (T < tfloor) {
	  edot[k][j][i] = 0.0;
	  continue;
	}

	/* Get rate of change of gas energy */
	edot[k][j][i] = pGrid->phrate[k][j][i] * egamma * n_H
	  - sw * dmc_cool_rate(n_Hplus/(n_H+n_Hplus), T)
	  - recomb_cool_rate_coef(T) * time_unit * n_Hplus * n_e
	  - sw * ki_cool_rate(T) * time_unit * n_H * n_H
	  + sw *ki_heat_rate() * time_unit * n_H;
	edot[k][j][i] /= enorm;
#if 0
	if (pGrid->phrate[k][j][i] > 0.0) {
	  printf("%d %d %d %g %g %g %g %g %g %g %g\n",i,j,k,T,e_thermal*enorm,edot[k][j][i],pGrid->phrate[k][j][i] * egamma * n_H,
		 dmc_cool_rate(n_Hplus/(n_H+n_Hplus), T),recomb_cool_rate_coef(T) * time_unit * n_Hplus * n_e,
		 ki_cool_rate(T) * time_unit * n_H * n_H,ki_heat_rate() * time_unit * n_H);
	}
#endif
	/* Compute thermal time step for this cell and find the
	   min. Note that, if we're cooling, we need to take into
	   account the effect of the floor. */
	if (edot[k][j][i] > 0.0) {

	  /* We're heating, no need to consider floor */
	  dt_therm1 = ((1.0+max_de)*e_init[k][j][i] - pGrid->U[k][j][i].E) 
	    / edot[k][j][i];
	  dt_therm2 = ((1.0+max_de_therm)*e_th_init[k][j][i] - e_thermal) 
	    / edot[k][j][i];

	} else {

#if 0
	  /* We're cooling. Start by computing the total and thermal
	     energy the gas would have if it were at the temperature
	     floor. */
	  e_sp_min = tfloor * kB / ((x*0.5*mH+(1.0-x)*mu) * Gamma_1);
	  e_th_min = e_sp_min * pGrid->U[k][j][i].d;
	  e_min = 
	    0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		   pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		   pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    + 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
	    + e_th_min;

	  /* If cooling to the temperature floor would not violate the
	     constraint on the maximum allowable change in either
	     total or thermal energy, there is no constraint, and we
	     can continue to the next cell */
	  if ((e_th_init[k][j][i]-e_th_min < max_de_therm*e_th_init[k][j][i]) &&
	      (e_init[k][j][i]-e_min < max_de*e_init[k][j][i])) continue;

	  /* If we're here, cooling to the temperature floor would
	     violate our time step constraint. Therefore compute the
	     time step normally. */
	  dt_therm1 = (e_init[k][j][i]/(1.0+max_de) - pGrid->U[k][j][i].E)
	    / edot[k][j][i];
	  dt_therm2 = (e_th_init[k][j][i]/(1.0+max_de_therm) - e_thermal) 
	    / edot[k][j][i];
#else
	  dt_therm1 = dt_therm2 = HUGE_NUMBER;
#endif

	}

	dt_therm = (dt_therm1 < dt_therm2) ? dt_therm1 : dt_therm2;

	/* Compare time step for this cell to minimum found so far */
	if (dt_therm_min > dt_therm) {
	  dt_therm_min = dt_therm;
#if 0
	  isave = i;
	  jsave = j;
	  ksave = k;
#endif
	}
      }
    }
  }

#ifdef MPI_PARALLEL
  /* Sync thermal timestep across processors */
  err = MPI_Allreduce(&dt_therm_min, &dt_therm_min_glob, 1, MP_RL, 
		      MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[compute_therm_rates]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt_therm_min = dt_therm_min_glob;
#endif

#if 0
  n_H = pGrid->U[ksave][jsave][isave].dn / mH;
  n_Hplus = (pGrid->U[ksave][jsave][isave].d - pGrid->U[ksave][jsave][isave].dn) / mH;
  n_e = n_Hplus + pGrid->U[ksave][jsave][isave].d * alphaC / (14.0 * mH);
  x = n_Hplus / (n_H + n_Hplus);
  e_thermal = pGrid->U[ksave][jsave][isave].E - 0.5 *
    (pGrid->U[ksave][jsave][isave].M1 * pGrid->U[ksave][jsave][isave].M1 +
     pGrid->U[ksave][jsave][isave].M2 * pGrid->U[ksave][jsave][isave].M2 +
     pGrid->U[ksave][jsave][isave].M3 * pGrid->U[ksave][jsave][isave].M3) 
    / pGrid->U[ksave][jsave][isave].d
#ifdef MHD
    - 0.5 * (pGrid->U[ksave][jsave][isave].B1c * pGrid->U[ksave][jsave][isave].B1c +
	     pGrid->U[ksave][jsave][isave].B2c * pGrid->U[ksave][jsave][isave].B2c +
	     pGrid->U[ksave][jsave][isave].B3c * pGrid->U[ksave][jsave][isave].B3c) 
#endif
    ;
  e_sp = e_thermal / pGrid->U[ksave][jsave][isave].d;
  T = Gamma_1 * e_sp * (x*0.5*mH+(1.0-x)*mu)/ kB;
  printf("dt_therm = %e, T = %f, T_init = %f, x = %f, ijk = %d %d %d\n",
	 dt_therm_min, T, T*e_th_init[ksave][jsave][isave]/e_thermal, 
	 n_Hplus/(n_H+n_Hplus), isave, jsave, ksave);
#endif  

  return(dt_therm_min);
}


/* Routine to floor temperatures */
void apply_temp_floor(GridS *pGrid) {
  int i,j,k;
  Real e_sp, e_thermal, T, x, n_H, n_Hplus;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute temperature */
	n_H = pGrid->U[k][j][i].dn / mH;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].dn) / mH;
	x = n_Hplus / (n_H + n_Hplus);
	e_thermal = pGrid->U[k][j][i].E - 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;
	e_sp = e_thermal / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*mH+(1.0-x)*mu)/ kB;

	if (T < tfloor) {
	  e_sp = tfloor * kB / ((x*0.5*mH+(1.0-x)*mu) * Gamma_1);
	  pGrid->U[k][j][i].E = 
	    0.5 * (pGrid->U[k][j][i].M1*pGrid->U[k][j][i].M1 +
		   pGrid->U[k][j][i].M2*pGrid->U[k][j][i].M2 +
		   pGrid->U[k][j][i].M3*pGrid->U[k][j][i].M3) 
	    / pGrid->U[k][j][i].d
#ifdef MHD
	    + 0.5 * (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].B1c +
		     pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].B2c +
		     pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].B3c)
#endif
	    + e_sp * pGrid->U[k][j][i].d;
	}
      }
    }
  }
}


/* Routine to save energy and thermal energy passed to routine -- used
   for time step constraint. */
void save_energy (GridS *pGrid)
{
  int i,j,k;
  Real e_thermal, x, n_H, n_Hplus;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	n_H = pGrid->U[k][j][i].dn / mH;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].dn) / mH;
	x = n_Hplus / (n_H + n_Hplus);
	e_thermal = pGrid->U[k][j][i].E - 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d
#ifdef MHD
	  - 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		   pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		   pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c) 
#endif
	  ;
	e_init[k][j][i] = pGrid->U[k][j][i].E;
	e_th_init[k][j][i] = e_thermal;	
      }
    }
  }
}

void ionization_update(GridS *pGrid, Real dt)
{
  int i, j, k;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Update gas energy */
	pGrid->U[k][j][i].E += edot[k][j][i] * dt;

	/* Update neutral density */
	pGrid->U[k][j][i].dn += nHdot[k][j][i] * dt * mH;

      }
    }
  }
}

Real recomb_rate_coef(Real T) {
  /* Osterbrock (1989), pg. 19, rough fit to table values. Note that
   * the table contains an error in the 20,000 K column for alpha_A or
   * alpha_B. Fit is taken from Rijkhorst, Plewa, Dewey, & Mellema
   * (2005), eqn. 18.
   */
#ifdef IONIZATION_ONLY
  return(0.0);
#else
  return(2.59e-13*pow(T/1.0e4, -0.7));
#endif
}

Real coll_ion_rate_coef(Real T) {
  /* Tenorio-Tagle et al. 1986, eqn 8 */
#ifdef IONIZATION_ONLY
  return(0.0);
#else
  return(5.84e-11*sqrt(T)*exp(-IHI/(KB*T)));
#endif
}

Real recomb_cool_rate_coef(Real T) {
  /* Recombination cooling rate from Osterbrock (1989),
     pg. 50-51. Numerical values are a linear fit (in log space) to
     Table 3.2, column 3. */
#ifdef IONIZATION_ONLY
  return(0.0);
#else
  if (T < 100.0) {
    return(0.0);
  } else {
    return(6.11e-10*pow(T,-0.89)*KB*T);
  }
#endif
}

#define SCALEFACTOR 1.0e-23 /* Scaling to physical units */
Real dmc_cool_rate(Real x, Real T) {
#ifndef IONIZATION_ONLY
  /* Equilibrium cooling rate from Dalgarno & McCray (1972) */
  /* Implementation follows Jim Stone's dmc routine */
  Real le, lh, u, u2, om, dom, p1;
  Real qq2, qt1, qt2, qt3, qt4, xu1, xu2, xu3, xu4, tlost, tcool;
  int ipps, jaug;

  if (x < 1.0e-3) x = 1.0e-3;

  /* Electron impact excitation luminosity (eqn 3-10) */
  le = 0.0;
  if (T > 10) {
    le += 2.96e-23/sqrt(T) * exp(-92.0/T);
    if (T > 50) {
      le += 6.08e-23/sqrt(T) * exp(-413.0/T)
	+ 3.52e-23/sqrt(T) *
	(exp(-554.0/T) + 1.3* exp(-961.0/T));
      if (T > 2.0e4) {
	le = le + 4.14e-26*sqrt(T) * exp(-22700.0/T)
	  +  7.13e-26*sqrt(T)*(1.0-2.7e-9*T*T)
	  *exp(-27700.0/T);
      }
    }
  }

  /* Hydrogen cooling */
  if (T > 50.0) {
    lh = 2.37e-27*exp(-413/T)
      + 3.52e-27*(exp(-554.0/T)+1.4*exp(-961.0/T));
  } else {
    lh = 0.0;
  }

  /* Neutral cooling */
  u = T/157890. < 3.16 ? T/157890. : 3.16;
  u2 = u*u;
  om=.6098+1.489*u+.50755*u2-.38145*u*u2+.10196*u2*u2
    -.01007*u*u2*u2;
  dom=(1.489+2.*.50755*u-3.*.38145*u2+4.*.10196*u2*u
       -5.*.01007*u2*u2)/157890.;
  p1 = 0.0;
  if (T > 1.0e4)
    p1 = 0.5*1.41e-16*om*exp(-118000./T) / sqrt(T);

  /* Cooling in various regimes: < 100 K */
  if (T < 100.0)
    /* return(0.0); */
    return(x*le + lh + 
	   (1.0-x)*p1);
  if (T < 1.0e4)
    return(SCALEFACTOR*x*2.8347e-10*pow(T - 1.0e+02, 2.3562)
	   + x*le + lh
	   +(1.0-x)*p1);
  if (T > 1.27717e8)
    return(x*2.3988e-04*sqrt(T));

  tlost = log(T)/log(10.0);
  ipps = (int) floor(10.0*tlost) - 38;
  if (ipps > 41) ipps = 41;
  jaug = 2 > ipps ? 2 : ipps;
  qq2 = 3.8 + 0.1*jaug;
  qt2 = tlost - qq2;
  qt3 = qt2 - 0.1;

  if ((jaug == 2) || (jaug == 41)) {
    tcool = (xmat[jaug-1]*qt2 - xmat[jaug-2]*qt3)*10.0;
    return(SCALEFACTOR*pow(10.0, tcool)*x + (1.0-x)*p1);
  }

  qt1 = qt2 + 0.1;
  qt4 = qt3 - 0.1;

  xu1 = qt2*qt3*qt4/6.0e-03;
  xu2 = qt1*qt3*qt4/2.0e-03;
  xu3 = qt1*qt2*qt4/2.0e-03;
  xu4 = qt1*qt2*qt3/6.0e-03;

  tcool = -xmat[jaug-3]*xu1 + xmat[jaug-2]*xu2 -
    xmat[jaug-1]*xu3 + xmat[jaug]*xu4;
  return(SCALEFACTOR*pow(10.0,tcool)*x + (1.0-x)*p1);
#else
  return(0.0);
#endif
}
#undef SCALEFACTOR

Real ki_cool_rate(Real T) {
  return(GAMMAKI * (1.0e7 * exp(-114800.0/(T + 1000.0)) +
	  14.0*sqrt(T)*exp(-92.0/T)));
}

Real ki_heat_rate() {
  return(GAMMAKI);
}

#endif /* PHOTOIONIZATION */
