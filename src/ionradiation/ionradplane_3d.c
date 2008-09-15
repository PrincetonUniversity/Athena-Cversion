#include "../copyright.h"
/*==============================================================================
 * FILE: ionradplane_3d.c
 *
 * PURPOSE: Contains functions to compute an ionization radiative transfer
 *   from point sources update, using the algorithm described in
 *   Krumholz, Stone, & Gardiner (2007), or transfer from a
 *   plane-parallel radiation front.
 *
 *   Use of these routines requires that --enable-ion-radiation be set
 *   at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   add_radplane_3d             - adds a new radiation source
 *   ion_radplane_init_domain_3d - handles internal initialization
 *   get_ph_rate_plane           - computed photoionzation rate from
 *                                    a planar source
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "prototypes.h"
#include "../prototypes.h"
#include "../globals.h"
#include "ionrad.h"

#ifdef ION_RADPLANE

/* Initialized number of radiation planes to zero */
void ion_radplane_init_domain_3d(Grid *pGrid, Domain *pDomain) {
  pGrid->nradplane = 0;
  return;
}

/* --------------------------------------------------------------
 * Routine to add new radiator plane
 * --------------------------------------------------------------
 */

void add_radplane_3d(Grid *pGrid, int dir, Real flux) {
  int n;

  /* Add radiator to pgrid structure */
  pGrid->nradplane++;
  if (pGrid->nradplane==1) {
    if (!(pGrid->radplanelist = malloc(sizeof(Radplane))))
      goto on_error;
  } else {
    if (!(pGrid->radplanelist = realloc(pGrid->radplanelist, 
					pGrid->nradplane*sizeof(Radplane))))
      goto on_error;
  }
  n = pGrid->nradplane-1;
  pGrid->radplanelist[n].dir = dir;
  pGrid->radplanelist[n].flux = flux;

  return;

 on_error:
  ath_error("[add_radplane_3d]: malloc returned a NULL pointer\n");
}

/* --------------------------------------------------------------
 * Routine to compute photoionization rate from a plane radiation
 * source.
 * --------------------------------------------------------------
 */

void get_ph_rate_plane(Real initflux, int dir, Real ***ph_rate, 
		       Grid *pGrid) {
  int n;
  int lr;
  Real tau, n_H, kph, etau, cell_len;
  Real flux, flux_frac;
  int i, j, k;
  int s, e;
#ifdef MPI_PARALLEL
  int nGrid;
  int myrank, nextproc, prevproc, err;
  int planesize;
  Real *planeflux = NULL;
  Real max_flux_frac, max_flux_frac_glob;
  MPI_Status stat;
#endif

  /* Set lr based on whether radiation is left or right
     propagating. lr = 1 is for radiation going left to right. Also
     set up the start and end indices based on the direction, and
     store the cell length. */
  lr = (dir < 0) ? 1 : -1;
  switch(dir) {
  case -1: case 1: {
    if (lr > 0) {
      s=pGrid->is; e=pGrid->ie;
    } else {
      s=pGrid->ie; e=pGrid->is;
    }
    cell_len = pGrid->dx1;
    break;
  }
  case -2: case 2: {
    if (lr > 0) {
      s=pGrid->js; e=pGrid->je;
    } else {
      s=pGrid->je; e=pGrid->js;
    }
    cell_len = pGrid->dx2;
    break;
  }
  case -3: case 3: {
    if (lr > 0) {
      s=pGrid->ks; e=pGrid->ke;
    } else {
      s=pGrid->ke; e=pGrid->ks;
    cell_len = pGrid->dx3;
    }
    break;
  }
  }

#ifdef MPI_PARALLEL
  /* Figure out processor geometry: where am I, where are my neighbors
     upstream and downstream, how many processors are there in the
     direction of radiation propagation, how big is the interface with
     my neighbor? */
  for (k=0; k<NGrid_x3; k++) {
    for (j=0; j<NGrid_x2; j++) {
      for (i=0; i<NGrid_x1; i++) {
	if (pD->GridArray[k][j][i].id == pGrid->my_id) {
	  switch(dir) {
	  case -1: case 1: {
	    nGrid=NGrid_x1;
	    myrank = lr > 0 ? i : nGrid - i - 1;
	    planesize=pGrid->Nx2*pGrid->Nx3;
	    if ((i-lr >= 0) && (i-lr <= NGrid_x1-1))
	      prevproc = pD->GridArray[k][j][i-lr].id;
	    else
	      prevproc = -1;
	    if ((i+lr >= 0) && (i+lr <= NGrid_x1-1))
	      nextproc = pD->GridArray[k][j][i+lr].id;
	    else
	      nextproc = -1;
	    break;
	  }
	  case -2: case 2: {
	    nGrid=NGrid_x2;
	    myrank = lr > 0 ? j : nGrid - j - 1;
	    planesize=pGrid->Nx1*pGrid->Nx2;
	    if ((j-lr >= 0) && (j-lr <= NGrid_x2-1))
	      prevproc = pD->GridArray[k][j-lr][i].id;
	    else
	      prevproc = -1;
	    if ((j+lr >= 0) && (j+lr <= NGrid_x2-1))
	      nextproc = pD->GridArray[k][j+lr][i].id;
	    else
	      nextproc = -1;
	    break;
	  }
	  case -3: case 3: {
	    nGrid=NGrid_x3;
	    myrank = lr > 0 ? k : nGrid - k - 1;
	    planesize=pGrid->Nx1*pGrid->Nx2;
	    if ((k-lr >= 0) && (k-lr <= NGrid_x3-1))
	      prevproc = pD->GridArray[k-lr][j][i].id;
	    else
	      prevproc = -1;
	    if ((k+lr >= 0) && (k+lr <= NGrid_x2-1))
	      nextproc = pD->GridArray[k+lr][j][i].id;
	    else
	      nextproc = -1;
	    break;
	  }
	  default:
	    ath_error("[get_ph_rate_plane]: dir must be +-1, 2, or 3\n");
	  }
	}
      }
      if (nGrid != 0) break;
    }
      if (nGrid != 0) break;
  }

  /* Allocate memory for flux at interface on first pass */
  if (!(planeflux=calloc(planesize, sizeof(Real))))
    ath_error("[get_ph_rate_plane]: calloc returned a null pointer!\n");

  /* Loop over processors in the direction of propagation */
  for (n=0; n<nGrid; n++) {

    /* Set to 0 flux fraction remaining to start */
    max_flux_frac = 0.0;

    /* Am I the rank before the current one? If so, pass the flux on
       to the next processor. */
    if (myrank == n-1) {
      err = MPI_Send(planeflux, planesize, MP_RL, nextproc, n,
		     MPI_COMM_WORLD);
      if (err) ath_error("[get_ph_rate_plane]: MPI_Send error = %d\n", err);
    } 

    /* Is it my turn to compute the transfer now? */
    if (myrank == n) {

      /* If I am not the first processor, get the flux from the
	 previous one */
      if (prevproc != -1) {
	err = MPI_Recv(planeflux, planesize, MP_RL, prevproc, n,
		       MPI_COMM_WORLD, &stat);
	if (err) ath_error("[get_ph_rate_plane]: MPI_Send error = %d\n", err);
      }
#endif /* MPI_PARALLEL */

      /* Propograte the radiation */
      switch(dir) {
      case -1: case 1: {
	for (k=pGrid->ks; k<=pGrid->ke; k++) {
	  for (j=pGrid->js; j<=pGrid->je; j++) {
#ifdef MPI_PARALLEL
	    /* Get initial flux from passed information or boundary
	       conditions */
	    if (prevproc != -1) 
	      flux = planeflux[(k-pGrid->ks)*pGrid->Nx2+j-pGrid->js];
	    else
#endif /* MPI_PARALLEL */
	      flux = initflux;
	    for (i=s; i<=e; i+=lr) {
	      n_H = pGrid->U[k][j][i].s[0] / m_H;
	      tau = sigma_ph * n_H * pGrid->dx1;
	      etau = exp(-tau);
	      kph = flux * (1.0-etau) / (n_H*cell_len);
	      ph_rate[k][j][i] += kph;
	      flux *= etau;
	      flux_frac = flux / initflux;
	      if (flux_frac < MINFLUXFRAC) break;
	    }
#ifdef MPI_PARALLEL
	    /* Store final flux to pass to next processor, or 0 if we
	       ended the loop early because we were below the minimum
	       fraction. */
	    planeflux[(k-pGrid->ks)*pGrid->Nx2+j-pGrid->js] = 
	      flux_frac < MINFLUXFRAC ? 0.0 : flux;
	    max_flux_frac = (flux_frac > max_flux_frac) ?
	      flux_frac : max_flux_frac;
#endif /* MPI_PARALLEL */
	  }
	}
	break;
      }
      case -2: case 2: {
	for (k=pGrid->ks; k<=pGrid->ke; k++) {
	  for (i=pGrid->is; i<=pGrid->ie; i++) {
#ifdef MPI_PARALLEL
	    /* Get initial flux from passed information or boundary
	       conditions */
	    if (prevproc != -1) 
	      flux = planeflux[(k-pGrid->ks)*pGrid->Nx1+i-pGrid->is];
	    else
#endif /* MPI_PARALLEL */
	      flux = initflux;
	    for (j=s; j<=e; j+=lr) {
	      n_H = pGrid->U[k][j][i].s[0] / m_H;
	      tau = sigma_ph * n_H * pGrid->dx1;
	      etau = exp(-tau);
	      kph = flux * (1.0-etau) / (n_H*cell_len);
	      ph_rate[k][j][i] += kph;
	      flux *= etau;
	      flux_frac = flux / initflux;
	      if (flux_frac < MINFLUXFRAC) break;
	    }
#ifdef MPI_PARALLEL
	    /* Store final flux to pass to next processor */
	    planeflux[(k-pGrid->ks)*pGrid->Nx1+i-pGrid->is] = 
	      flux_frac < MINFLUXFRAC ? 0.0 : flux;
	    max_flux_frac = (flux_frac > max_flux_frac) ?
	      flux_frac : max_flux_frac;
#endif /* MPI_PARALLEL */
	  }
	}
	break;
      }
      case -3: case 3: {
	for (j=pGrid->ks; j<=pGrid->ke; j++) {
	  for (i=pGrid->is; i<=pGrid->ie; i++) {
#ifdef MPI_PARALLEL
	    /* Get initial flux from passed information or boundary
	       conditions */
	    if (prevproc != -1) 
	      flux = planeflux[(j-pGrid->js)*pGrid->Nx1+i-pGrid->is];
	    else
#endif /* MPI_PARALLEL */
	      flux = initflux;
	    for (k=s; k<=e; k+=lr) {
	      n_H = pGrid->U[k][j][i].s[0] / m_H;
	      tau = sigma_ph * n_H * pGrid->dx1;
	      etau = exp(-tau);
	      kph = flux * (1.0-etau) / (n_H*cell_len);
	      ph_rate[k][j][i] += kph;
	      flux *= etau;
	      flux_frac = flux / initflux;
	      if (flux_frac < MINFLUXFRAC) break;
	    }
#ifdef MPI_PARALLEL
	    /* Store final flux to pass to next processor */
	    planeflux[(j-pGrid->js)*pGrid->Nx1+i-pGrid->is] = 
	      flux_frac < MINFLUXFRAC ? 0.0 : flux;
	    max_flux_frac = (flux_frac > max_flux_frac) ?
	      flux_frac : max_flux_frac;
#endif /* MPI_PARALLEL */
	  }
	}
	break;
      }
      }

#ifdef MPI_PARALLEL
    }

    /* If we're parallel, get the maximum flux fraction left and see
       if we should continue to the next set of processors. */
    err = MPI_Allreduce(&max_flux_frac, &max_flux_frac_glob, 1, MP_RL,
			MPI_MAX, MPI_COMM_WORLD);
    if (err) ath_error("[get_ph_rate_plane]: MPI_Allreduce error = %d\n", err);
    if (max_flux_frac_glob < MINFLUXFRAC) break;
  }

  free(planeflux);
#endif /* MPI_PARALLEL */

  return;
}

#endif /* ION_RADPLANE */
