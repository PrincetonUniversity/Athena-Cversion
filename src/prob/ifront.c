/*
 * Function ifront.c
 *
 * Problem generator for ionization front propagating from domain edge
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#define TII 1.0e4
#define ALPHA4 2.59e-13

void problem(Grid *pGrid, Domain *pDomain)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real cs, n_H, m_H, flux;
  Real rho, pressure;
#ifdef MHD
  Real bx, by, bz;
#endif

/* Set up uniform ambinet medium with ionizing source at its
 * edge. The user-input parameters are n_H (initial density),
 * cs (initial sound speed in neutral gas), flux (ionizing flux,
 * in units of photons per unit time per unit area), and bx, by, and
 * bz (initial vector magnetic field).
 */
  m_H = par_getd("ionradiation", "m_H");
  n_H = par_getd("problem","n_H");
  cs = par_getd("problem","cs");
  flux = par_getd("problem","flux");
#ifdef MHD
  bx = par_getd("problem","Bx");
  by = par_getd("problem","By");
  bz = par_getd("problem","Bz");
#endif

  /* Power-law pressure and density */
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
	rho = n_H*m_H;
	pressure = rho*cs*cs;
	pGrid->U[k][j][i].d  = rho;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;
	pGrid->U[k][j][i].E  = pressure/Gamma_1;
	pGrid->U[k][j][i].s[0] = rho;
#ifdef MHD
	pGrid->U[k][j][i].B1c = bx;
	pGrid->U[k][j][i].B2c = by;
	pGrid->U[k][j][i].B3c = bz;
	pGrid->B1i[k][j][i] = bx;
	pGrid->B2i[k][j][i] = by;
	pGrid->B3i[k][j][i] = bz;
	pGrid->U[k][j][i].E += (bx*bx+by*by+bz*bz)/2.0;
#endif
      }
    }
  }

  /* Radiation originating from domain edge */
  add_radplane_3d(pGrid, -1, flux);

  return;
}

/*
 * 'Userwork_before_loop' and 'Userwork_after_loop' are executed
 * before and after the main integration loop.
 */

void Userwork_before_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}


void problem_write_restart(Grid *pG, Domain *pDomain, FILE *fp){
  return;
}


void problem_read_restart(Grid *pG, Domain *pDomain, FILE *fp){
  return;
}


Gasfun_t get_usr_expr(const char *expr){
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}
