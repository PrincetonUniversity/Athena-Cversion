/*
 * Function hiireg.c
 *
 * Problem generator for 3-D spherical HII region
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
  Real r0, krho, cs, rs, s, n_H0, m_H, alpha_B;
  Real x1, x2, x3, r, rho, pressure;
#ifdef MHD
  Real b0;
#endif

/* Set up uniform ambinet medium with ionizing source at its
 * center. The user-input parameters are n_H0 (central density), r0
 * (scale radius), krho (density power law index), cs (sound speed in
 * unionized gas), and s (ionizing luminosity, in units of photons per
 * unit time.
 */
  m_H = par_getd("rayradiation", "m_H");
  n_H0 = par_getd("problem","n_H0");
  r0 = par_getd("problem","r0");
  krho = par_getd("problem","krho");
  cs = par_getd("problem","cs");
  rs = par_getd("problem","rs");
#ifndef IONIZATION_ONLY
  alpha_B = recomb_rate_coef(TII);
  s = 4./3. * PI * rs*rs*rs * n_H0*n_H0 * alpha_B;
#else
  s = 4./3. * PI * rs*rs*rs * n_H0*n_H0 * ALPHA4;
#endif
#ifdef MHD
  b0 = par_getd("problem","B0");
#endif

  /* Power-law pressure and density */
  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	r = sqrt(x1*x1 + x2*x2 + x3*x3);
	rho = n_H0*m_H*pow(r/r0,-krho);
	pressure = rho*cs*cs;
	pGrid->U[k][j][i].d  = rho;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;
	pGrid->U[k][j][i].E  = pressure/Gamma_1;
	pGrid->U[k][j][i].s[0] = rho;
#ifdef MHD
	pGrid->U[k][j][i].B1c = b0;
	pGrid->U[k][j][i].B2c = 0.0;
	pGrid->U[k][j][i].B3c = 0.0;
	pGrid->B1i[k][j][i] = b0;
	pGrid->B2i[k][j][i] = 0.0;
	pGrid->B3i[k][j][i] = 0.0;
	pGrid->U[k][j][i].E += b0*b0/2.0;
#endif
      }
    }
  }

  /* Single radiator at origin */
  add_radpoint_3d(pGrid, 0, 0, 0, s);

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
