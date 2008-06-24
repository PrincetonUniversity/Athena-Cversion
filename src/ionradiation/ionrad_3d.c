#include "copyright.h"
/*==============================================================================
 * FILE: ionrad_3d.c
 *
 * PURPOSE: Contains functions to compute an ionization radiative transfer
 *   from update, using the algorithm described in Krumholz, Stone, 
 *   & Gardiner (2007).
 *
 *   Use of these routines requires that --enable-ion-radiation be set
 *   at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   ion_radtransfer_3d             - does an ionizing radiative transfer
 *                                      update
 *   ion_radtransfer_init_3d        - handles internal initialization
 *   ion_radtransfer_init_domain_3d - handles internal initialization
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "ionrad.h"
#include "prototypes.h"

#ifdef ION_RADIATION

/* Global storage arrays */
static Real ***ph_rate;            /* Photoionization rate */
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
static Real ***x_init;             /* Ionization fraction on entry to
				      routine */

/* ------------------------------------------------------------
 * Photoionization routines
 * ------------------------------------------------------------
 *
 * These routines do the work of computing the photoionization and
 * chemisty update.
 *
 */

/* Routine to zero out initial photoionization rates */
void ph_rate_init(Grid *pGrid)
{
  int i,j,k;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) ph_rate[k][j][i] = 0.0;
    }
  }
}

/* Routine to floor temperatures */
void apply_temp_floor(Grid *pGrid) {
  int i,j,k;
  Real e_sp, e_thermal, ke, T, x, n_H, n_Hplus;
#ifdef MHD
  Real be;
#endif

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute temperature */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	x = n_Hplus / (n_H + n_Hplus);
	ke = 0.5 *
	  (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 +
	   pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 +
	   pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3) 
	  / pGrid->U[k][j][i].d;
#ifdef MHD
	be = 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c +
		    pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c +
		    pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c);
#endif
	e_thermal = pGrid->U[k][j][i].E - ke;
#ifdef MHD
	e_thermal -= be; 
#endif
	e_sp = e_thermal / pGrid->U[k][j][i].d;
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;

	if (T < tfloor) {
	  e_sp = tfloor * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
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

	if ((T > tceil) && (tceil > 0)) {
	  e_sp = tceil * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
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


/* Routine to keep d_n > floor and d_n < d. */
void apply_neutral_floor(Grid *pGrid) {
  int i,j,k;
  Real d_nlim;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;
	if (pGrid->U[k][j][i].s[0] < d_nlim) {
	  pGrid->U[k][j][i].s[0] = d_nlim;
	} else if (pGrid->U[k][j][i].s[0] > pGrid->U[k][j][i].d) {
	  pGrid->U[k][j][i].s[0] = pGrid->U[k][j][i].d;
	}
      }
    }
  }
}


/* Routine to save energy, thermal energy, and ionization fraction
   passed to routine -- used for time step constraint. */
void save_energy_and_x(Grid *pGrid)
{
  int i, j, k;
  Real e_thermal, n_H, n_Hplus, n_e, x;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Compute thermal energy */
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

	/* Compute ion fraction */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	x = n_Hplus / (n_H + n_Hplus);

	/* Save thermal and total energy, and neutral fraction */
	e_init[k][j][i] = pGrid->U[k][j][i].E;
	e_th_init[k][j][i] = e_thermal;
	x_init[k][j][i] = pGrid->U[k][j][i].s[0] / pGrid->U[k][j][i].d;

	/* Initialize last_sign and sign_count arrays */
	last_sign[k][j][i] = 0;
	sign_count[k][j][i] = 0;
      }
    }
  }
}


/* Routine to check if we have changed the total energy, thermal
   energy, or x_n as much as we are allowed. */
int check_range(Grid *pGrid) {
  int i, j, k;
  Real e_thermal, n_H, n_Hplus, n_e, x;
  long cellcount = 0;
#ifdef MPI_PARALLEL
  int err;
  long cellcount_glob = 0;
#endif

  /* Check thermal energy */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Check D type condition */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	if (ph_rate[k][j][i] / (min_area * n_H) > 2.0*CION) continue;

	/* Check thermal energy condition */
	if (max_de_therm_step > 0) {
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
	}
	if ((e_thermal / e_th_init[k][j][i] >= 1 + max_de_therm_step) ||
	    (e_th_init[k][j][i] / e_thermal >= 1 + max_de_therm_step)) {
	  cellcount++;
	  continue;
	}

	/* Check total energy condition */
	if (max_de_step > 0) {
	  if ((pGrid->U[k][j][i].E / e_init[k][j][i] >= 1 + max_de_step) ||
	      (e_init[k][j][i] / pGrid->U[k][j][i].E >= 1 + max_de_step)) {
	    cellcount++;
	    continue;
	  }
	}

	/* Check neutral fraction condition */
	if (max_dx_step > 0) {
	  n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	  n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
	  x = n_Hplus / (n_H + n_Hplus);
	  if ((x / x_init[k][j][i] >= 1 + max_dx_step) ||
	      (x_init[k][j][i] / x >= 1 + max_dx_step)) {
	    cellcount++;
	    continue;
	  }
	}
      }
    }
  }

#ifdef MPI_PARALLEL
  err = MPI_Allreduce(&cellcount, &cellcount_glob, 1, MPI_LONG, 
		      MPI_SUM, MPI_COMM_WORLD);
  if (err) ath_error("[check_range]: MPI_Allreduce returned error code %d\n"
		    ,err);
  cellcount = cellcount_glob;
#endif
  if (cellcount > MAXCELLCOUNT) return(1);
  else return(0);
}


#define MAXSIGNCOUNT 4
#define DAMPFACTOR 0.5
Real compute_chem_rates(Grid *pGrid)
{
  int i, j, k, n;
  Real n_H, n_Hplus, n_e, d_nlim;
  Real e_sp, T, x;
  Real dt_chem, dt_chem1, dt_chem2, dt_chem_min;
#ifdef MPI_PARALLEL
  int err;
  Real dt_chem_min_glob;
#endif

  /* Initialize chemistry time step to large time step */
  dt_chem_min = LARGE;

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Get species abundances */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
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
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;
	if (T < tfloor) T = tfloor;

	/* Get rate of change of neutral density */
	nHdot[k][j][i] = 
	  recomb_rate_coef(T) * time_unit * n_e * n_Hplus
	  - ph_rate[k][j][i] * n_H
	  - coll_ion_rate_coef(T) * time_unit * n_e * n_H;

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

	/* Get ionization fraction floor */
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;

	/* Compute chemistry time step for this cell and find the min */
	if (nHdot[k][j][i] == 0.0) {
	  dt_chem1 = dt_chem2 = LARGE;
	} else if (nHdot[k][j][i] > 0.0) {
	  dt_chem1 = max_dx_iter / (1+max_dx_iter) * n_e / nHdot[k][j][i];
	  dt_chem2 = max_dx_iter * n_H / nHdot[k][j][i];
	} else if (pGrid->U[k][j][i].s[0] > 1.0001*d_nlim) {
	  dt_chem1 = -max_dx_iter * n_e / nHdot[k][j][i];
	  dt_chem2 = -max_dx_iter / (1+max_dx_iter) * n_H / nHdot[k][j][i];
	} else {
	  dt_chem1 = dt_chem2 = LARGE;
	}
	dt_chem = (dt_chem1 < dt_chem2) ? dt_chem1 : dt_chem2;
	dt_chem_min = (dt_chem < dt_chem_min) ? dt_chem : dt_chem_min;

	if (dt_chem < 0) {
	  ath_error("[compute_chem_rates]: cell %d %d %d: dt_chem = %e, d = %e, d_n = %e, T = %e, nH = %e, nH+ = %e, ne = %e, nHdot = %e\n", i,j,k, dt_chem, pGrid->U[k][j][i].d, pGrid->U[k][j][i].s[0], T, n_H, n_Hplus, n_e, nHdot[k][j][i]);
	}
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

  /*fclose(fp);*/

  return(dt_chem_min);
}
#undef MAXSIGNCOUNT
#undef DAMPFACTOR


Real compute_therm_rates(Grid *pGrid)
{
  int i, j, k;
  Real n_H, n_Hplus, n_e, e_thermal;
  Real e_sp, T, x, e_sp_min, e_th_min, e_min, d_nlim;
  Real dt_therm, dt_therm1, dt_therm2, dt_therm_min;
#ifdef MPI_PARALLEL
  int err;
  Real dt_therm_min_glob;
#endif

  /* Initialize thermal time step to large value */
  dt_therm_min = LARGE;

  /* Loop over cells to get timestep */
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	/* Get species abundances */
	n_H = pGrid->U[k][j][i].s[0] / m_H;
	n_Hplus = (pGrid->U[k][j][i].d - pGrid->U[k][j][i].s[0]) / m_H;
	n_e = n_Hplus + pGrid->U[k][j][i].d * alpha_C / (14.0 * m_H);
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
	T = Gamma_1 * e_sp * (x*0.5*m_H+(1.0-x)*mu)/ k_B;

	/* Check temperature floor. If cell is below temperature
	   floor, skip it. We'll fix it later */
	if (T < tfloor) {
	  edot[k][j][i] = 0.0;
	  continue;
	}

	/* If we're at the ionization floor and trying to ionize
	   further, we won't update the temperature, so skip this
	   cell. */
	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;
	if ((nHdot[k][j][i] < 0) &&
	    (pGrid->U[k][j][i].s[0] < 1.0001*d_nlim)) {
	  edot[k][j][i] = 0.0;
	  continue;
	}

	/* Get rate of change of gas energy. Only use molecular cooling
	   in cells with < COOLFRAC or > 1 - COOLFRAC ionization fraction, to
	   avoid artificial over-cooling in mixed or transition cells. */
	edot[k][j][i] = ph_rate[k][j][i] * e_gamma * n_H
	  - osterbrock_cool_rate(T) * n_e*n_Hplus
	  + recomb_cool_rate_coef(T) * time_unit * n_Hplus * n_e;
	if ((n_Hplus / (n_H+n_Hplus) < COOLFRAC) || 
	    (n_Hplus / (n_H+n_Hplus) > 1.0-COOLFRAC)) {
	  edot[k][j][i] += ki_heat_rate() * time_unit * n_H
	    - ki_cool_rate(T) * time_unit * n_H * n_H;
	}

	/* Compute thermal time step for this cell and find the
	   min. Note that, if we're cooling, we need to take into
	   account the effect of the floor. */
	if (edot[k][j][i] == 0.0) {
	  dt_therm1 = dt_therm2 = LARGE;
	} else if (edot[k][j][i] > 0.0) {

	  /* We're heating, no need to consider floor */
	  dt_therm1 = max_de_iter * pGrid->U[k][j][i].E / edot[k][j][i];
	  dt_therm2 = max_de_therm_iter * e_thermal / edot[k][j][i];

	} else {

	  /* We're cooling. Start by computing the total and thermal
	     energy the gas would have if it were at the temperature
	     floor. */
	  e_sp_min = tfloor * k_B / ((x*0.5*m_H+(1.0-x)*mu) * Gamma_1);
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
	  if ((e_thermal/(1.0+max_de_therm_iter) < e_th_min) &&
	      (pGrid->U[k][j][i].E/(1.0+max_de_iter) < e_min))
	    continue;

	  /* If we're here, cooling to the temperature floor would
	     violate our time step constraint. Therefore compute the
	     time step normally. */
	  dt_therm1 = -max_de_iter / (1+max_de_iter) * pGrid->U[k][j][i].E 
	    / edot[k][j][i];
	  dt_therm2 = -max_de_therm_iter / (1+max_de_therm_iter) * e_thermal /
	    edot[k][j][i];

	}

	/* Set time step to minimum */
	dt_therm = (dt_therm1 < dt_therm2) ? dt_therm1 : dt_therm2;
	dt_therm_min = (dt_therm < dt_therm_min) ? dt_therm : dt_therm_min;
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

  return(dt_therm_min);
}


void ionization_update(Grid *pGrid, Real dt)
{
  int i, j, k;
  Real d_nlim;

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {

	d_nlim = pGrid->U[k][j][i].d*IONFRACFLOOR;
	d_nlim = d_nlim < d_nlo ? d_nlim : d_nlo;

	if ((nHdot[k][j][i] > 0) ||
	    (pGrid->U[k][j][i].s[0] > 1.0001*d_nlim)) {

	  /* Update gas energy */
	  pGrid->U[k][j][i].E += edot[k][j][i] * dt;

	  /* Update neutral density */
	  pGrid->U[k][j][i].s[0] += nHdot[k][j][i] * dt * m_H;

	}
      }
    }
  }
}


Real compute_dt_hydro(Grid *pGrid) {
  int i,j,k;
  Real di,v1,v2,v3,qsq,p,asq,cf1sq,cf2sq,cf3sq,max_dti=0.0,dt;
#ifdef MHD
  Real b1,b2,b3,bsq,tsum,tdif;
#endif /* MHD */
#ifdef MPI_PARALLEL
  Real dt_glob;
  int err;
#endif /* MPI_PARALLEL */

  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    for (j=pGrid->js; j<=pGrid->je; j++) {
      for (i=pGrid->is; i<=pGrid->ie; i++) {
	di = 1.0/(pGrid->U[k][j][i].d);
	v1 = pGrid->U[k][j][i].M1*di;
	v2 = pGrid->U[k][j][i].M2*di;
	v3 = pGrid->U[k][j][i].M3*di;
	qsq = v1*v1 + v2*v2 + v3*v3;

#ifdef MHD
	/* Use maximum of face-centered fields (always larger than
	   cell-centered B) */
	b1 = pGrid->U[k][j][i].B1c 
	  + fabs((double)(pGrid->B1i[k][j][i] - pGrid->U[k][j][i].B1c));
	b2 = pGrid->U[k][j][i].B2c 
	  + fabs((double)(pGrid->B2i[k][j][i] - pGrid->U[k][j][i].B2c));
	b3 = pGrid->U[k][j][i].B3c 
	  + fabs((double)(pGrid->B3i[k][j][i] - pGrid->U[k][j][i].B3c));
	bsq = b1*b1 + b2*b2 + b3*b3;
	/* compute sound speed squared */
#ifdef ADIABATIC
	p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq
			 - 0.5*bsq), TINY_NUMBER);
	asq = Gamma*p*di;
#else
	asq = Iso_csound2;
#endif /* ADIABATIC */
	/* compute fast magnetosonic speed squared in each direction */
	tsum = bsq*di + asq;
	tdif = bsq*di - asq;
	cf1sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b2*b2+b3*b3)*di));
	cf2sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b3*b3)*di));
	cf3sq = 0.5*(tsum + sqrt(tdif*tdif + 4.0*asq*(b1*b1+b2*b2)*di));
	
#else /* MHD */
	
	/* compute sound speed squared */
#ifdef ADIABATIC
	p = MAX(Gamma_1*(pGrid->U[k][j][i].E - 0.5*pGrid->U[k][j][i].d*qsq),
		TINY_NUMBER);
	asq = Gamma*p*di;
#else
	asq = Iso_csound2;
#endif /* ADIABATIC */
	/* compute fast magnetosonic speed squared in each direction */
	cf1sq = asq;
	cf2sq = asq;
	cf3sq = asq;
	
#endif /* MHD */

	/* compute maximum inverse of dt (corresponding to minimum dt) */
	if (pGrid->Nx1 > 1)
	  max_dti = MAX(max_dti,(fabs(v1)+sqrt((double)cf1sq))/pGrid->dx1);
	if (pGrid->Nx2 > 1)
	  max_dti = MAX(max_dti,(fabs(v2)+sqrt((double)cf2sq))/pGrid->dx2);
	if (pGrid->Nx3 > 1)
	  max_dti = MAX(max_dti,(fabs(v3)+sqrt((double)cf3sq))/pGrid->dx3);
      }
    }
  }

  /* get timestep. */
  dt = CourNo/max_dti;

#ifdef MPI_PARALLEL
  err = MPI_Allreduce(&dt, &dt_glob, 1, MP_RL, MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[compute_dt_hydro]: MPI_Allreduce returned error code %d\n"
		    ,err);
  dt = dt_glob;
#endif /* MPI_PARALLEL */
  return(dt);
}

/* ------------------------------------------------------------
 * Initialize and store routines
 * ------------------------------------------------------------
 *
 * These are called at problem setup to allocate memory and read
 * setup parameters, or when checkpointing to dump internal data
 * needed on restart.
 *
 */

void ion_radtransfer_init_3d(Grid *pGrid, Domain *pDomain, int ires) {
  Real area1, area2, area3, maxdx;
  int j, k;

  /* Read input values  */
  sigma_ph = par_getd("ionradiation", "sigma_ph");
  m_H = par_getd("ionradiation", "m_H");
  mu = par_getd("ionradiation", "mu");
  e_gamma = par_getd("ionradiation", "e_gamma");
  alpha_C = par_getd("ionradiation", "alpha_C");
  k_B = par_getd("ionradiation", "k_B");
  time_unit = par_getd("ionradiation", "time_unit");
  max_de_iter = par_getd("ionradiation", "max_de_iter");
  max_de_therm_iter = par_getd("ionradiation", "max_de_therm_iter");
  max_dx_iter = par_getd("ionradiation", "max_dx_iter");
  max_de_step = par_getd("ionradiation", "max_de_step");
  max_de_therm_step = par_getd("ionradiation", "max_de_therm_step");
  max_dx_step = par_getd("ionradiation", "max_dx_step");
  tfloor = par_getd("ionradiation", "tfloor");
  tceil = par_getd("ionradiation", "tceil");
  maxiter = par_getd("ionradiation", "maxiter");

  /* Allocate memory for rate arrays */
  ph_rate = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(Real));
  edot = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(Real));
  nHdot = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(Real));
  last_sign = (int***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(int));
  sign_count = (int***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(int));
  e_init = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1,  
		    sizeof(Real));
  e_th_init = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(Real));
  x_init = (Real***) 
    calloc_3d_array(pGrid->Nx3, pGrid->Nx2, pGrid->Nx1, 
		    sizeof(Real));

  /* Offset pointers to account for ghost cells */
  ph_rate -= pGrid->ks;
  edot -= pGrid->ks;
  nHdot -= pGrid->ks;
  last_sign -= pGrid->ks;
  sign_count -= pGrid->ks;
  e_init -= pGrid->ks;
  e_th_init -= pGrid->ks;
  x_init -= pGrid->ks;
  for (k=pGrid->ks; k<=pGrid->ke; k++) {
    ph_rate[k] -= pGrid->js;
    edot[k] -= pGrid->js;
    nHdot[k] -= pGrid->js;
    last_sign[k] -= pGrid->js;
    sign_count[k] -= pGrid->js;
    e_init[k] -= pGrid->js;
    e_th_init[k] -= pGrid->js;
    x_init[k] -= pGrid->js;
    for (j=pGrid->js; j<=pGrid->je; j++) {
      ph_rate[k][j] -= pGrid->is;
      edot[k][j] -= pGrid->is;
      nHdot[k][j] -= pGrid->is;
      last_sign[k][j] -= pGrid->is;
      sign_count[k][j] -= pGrid->is;
      e_init[k][j] -= pGrid->is;
      e_th_init[k][j] -= pGrid->is;
      x_init[k][j] -= pGrid->is;
    }
  }

#ifdef ION_RADPOINT
  /* Initialize random number generator if this is not a restart. If
     this is a restart, we have already set up the random number
     generator from the restart file. */
  if (ires==0) ion_radpoint_init_ranstate();
#endif /* ION_RADPOINT */

  /* What's the smallest area a cell face can have? */
  area1 = pGrid->dx1 * pGrid->dx2;
  area2 = pGrid->dx1 * pGrid->dx3;
  area3 = pGrid->dx2 * pGrid->dx3;
  if (area1 < area2) {
    if (area1 < area3) min_area = area1;
    else min_area = area3;
  } else {
    if (area2 < area3) min_area = area2;
    else min_area = area3;
  }

  /* What's the "low" neutral density, corresponding to the minimum
     optical depth we care about? */
  maxdx = pGrid->dx1 > pGrid->dx2 ? pGrid->dx1 : pGrid->dx2;
  maxdx = maxdx > pGrid->dx3 ? maxdx : pGrid->dx2;
  d_nlo = MINOPTDEPTH * m_H / (sigma_ph * maxdx);

  return;
}


void ion_radtransfer_init_domain_3d(Grid *pGrid, Domain *pDomain) {

  /* Store parallel grid information for internal use */
#ifdef MPI_PARALLEL
  pD = pDomain;
  NGrid_x1 = par_geti("parallel","NGrid_x1");
  NGrid_x2 = par_geti("parallel","NGrid_x2");
  NGrid_x3 = par_geti("parallel","NGrid_x3");
#endif

  /* Store information specific to point sources and planes */
#ifdef ION_RADPOINT
  ion_radpoint_init_domain_3d(pGrid, pDomain);
#endif
#ifdef ION_RADPLANE
  ion_radplane_init_domain_3d(pGrid, pDomain);
#endif
}


/* ------------------------------------------------------------
 * Main integration routine
 * ------------------------------------------------------------
 *
 * This is the driver routine for the radiation integration step.
 *
 */

void ion_radtransfer_3d(Grid *pGrid) 
{
  Real dt_chem, dt_therm, dt_hydro, dt, dt_done;
  int n, niter, hydro_done;
  int nchem, ntherm;

#ifdef ION_RADPOINT
  /* Build or rebuild trees */
  refresh_trees(pGrid);
#endif /* ION_RADPOINT */

  /* Set all temperatures below the floor to the floor */
  apply_temp_floor(pGrid);

  /* Set neutral densities below floor to sensible values. This is necessary
     because the hydro update can make the neutral density negative. */
  apply_neutral_floor(pGrid);

  /* Save total and thermal energies and neutral fractions passed in
     -- we use these to compute time steps. Also initialize the last_sign
     and sign_count arrays */
  save_energy_and_x(pGrid);

#if defined(MPI_PARALLEL) && defined(ION_RADPOINT)
  /* Allocate buffer for MPI communication and attach here. Note that
     we detatch below, to ensure that the buffer we use here doesn't
     get tangled up with buffers that may be used in other
     modules. This is necessary because only one MPI buffer is allowed
     per process. We only use the buffer for point sources, since it     
     is easier to handle plane fronts using non-buffered communication. */
  if (mpi_bufsize==0) {
    mpi_bufsize = INIT_MPI_BUFSIZE;
    mpi_buffer = malloc(mpi_bufsize);
  }
  MPI_Buffer_attach(mpi_buffer, mpi_bufsize);
#endif

  /* Begin the radiation sub-cycle */
  dt_done = 0.0;
  hydro_done = 0;
  nchem = ntherm = 0;
  for (niter = 0; niter<maxiter;) {

    /* Initialize photoionization rate array */
    ph_rate_init(pGrid);

    /* Compute photoionization rate from all sources */
#ifdef ION_RADPOINT
    for (n=0; n<pGrid->nradpoint; n++) 
      get_ph_rate_point(pGrid->radpointlist[n].s, 
			&(pGrid->radpointlist[n].tree),
			ph_rate, pGrid);
#endif
#ifdef ION_RADPLANE
    for (n=0; n<pGrid->nradplane; n++) 
      get_ph_rate_plane(pGrid->radplanelist[n].flux,
			pGrid->radplanelist[n].dir,
			ph_rate, pGrid);
#endif

    /* Compute rates and time step for chemistry update */
    dt_chem = compute_chem_rates(pGrid);

    /* Compute rates and time step for thermal energy update */
    dt_therm = compute_therm_rates(pGrid);

    /* Set time step to smaller of thermal and chemical time
       steps, and record whether this is a thermal or chemical step */
    if (dt_chem < dt_therm) nchem++;
    else ntherm++;
    dt = MIN(dt_therm, dt_chem);

    /* If necessary, scale back time step to avoid exceeding hydro
       time step. */
    if (dt_done + dt > pGrid->dt) {
      dt = pGrid->dt - dt_done;
      hydro_done = 1;
    }

    /* Do an update */
    ionization_update(pGrid, dt);
    dt_done += dt;
    niter++;

    /* Set all temperatures below the floor to the floor */
    apply_temp_floor(pGrid);
    apply_neutral_floor(pGrid);

    /* Check new energies and ionization fractions against initial
       values to see if we've changed them as much as possible. If so,
       exit loop. */
    if (check_range(pGrid)) {
      pGrid->dt = dt_done;
      break;
    }

    /* Have we advanced the full hydro time step? If so, exit loop. */
    if (hydro_done) {
      break;
    }

    /* Compute a new hydro time step based on the new temperature
       distribution. If it's smaller than the time step we've already 
       advanced, then exit. */
    dt_hydro = compute_dt_hydro(pGrid);
    if (dt_hydro < dt_done) {
      pGrid->dt = dt_done;
      break;
    }
  }

  /* Have we exceeded the maximum number of iterations? If so, return
     to hydro with the time step re-set to what we managed to do. */
  if (niter==maxiter)
    pGrid->dt = dt_done;

#if defined(MPI_PARALLEL) && defined(ION_RADPOINT)
  /* Deallocate MPI buffer */
  MPI_Buffer_detach(&mpi_buffer, &mpi_bufsize);
#endif

  /* Write status */
  ath_pout(0, "Radiation done in %d iterations: %d thermal, %d chemical; new dt = %e\n", 
	   niter, ntherm, nchem, pGrid->dt);

  /* Sanity check */
  if (pGrid->dt < 0) {
    ath_error("[ion_radtransfer_3d]: dt = %e, dt_chem = %e, dt_therm = %e, dt_hydro = %e, dt_done = %e\n", pGrid->dt, dt_chem, dt_therm, dt_hydro, dt_done);
  }
}

#endif /* ION_RADIATION */
