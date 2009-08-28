#include "../copyright.h"
/*=============================================================================
FILE: utils.c
PURPOSE: Contains most of the utilities for the particle code: all the
  interpolation functions, stopping time calculation functions, shuffle
  algorithms. Also contained are the default (and trivial) gas velocity
  shift function. The get_gasinfo(Grid *pG) routine is used for test
  purposes only.

CONTAINS PUBLIC FUNCTIONS:

  void getwei_linear(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
  void getwei_TSC   (Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
  void getwei_QP    (Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
  int  getvalues(Grid *pG, Real weight[3][3][3], int is, int js, int ks, Real *rho, Real *u1, Real *u2, Real *u3, Real *cs);

  Real get_ts_epstein(Grid *pG, int type, Real rho, Real cs, Real vd);
  Real get_ts_general(Grid *pG, int type, Real rho, Real cs, Real vd);
  Real get_ts_fixed  (Grid *pG, int type, Real rho, Real cs, Real vd);

  void get_gasinfo(Grid *pG);
  void feedback_clear(Grid *pG);
  void distrFB      (Grid *pG, Real weight[3][3][3], int is, int js, int ks, Vector fb);

  void shuffle(Grid* pG);

  void gasvshift_zero(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3);
  void User_ParticleForce_Zero(Vector *ft, Real x1, Real x2, Real x3, Real *w1, Real *w2, Real *w3);

History:
  Written by Xuening Bai, Apr. 2009

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

int compare_gr(Grid *pG, Vector cell1, Grain gr1, Grain gr2);
void quicksort_particle(Grid *pG, Vector cell1, long start, long end);


/*=================================== ALL FUNCTIONS ======================================*/
/*----------------------------------------------------------------------------------------
 * interpolation functions;
 * stopping time functions;
 * feedback related functions;
 * shuffle related functions;
 *----------------------------------------------------------------------------------------*/


/*============================================================================*/
/*-----------------------------INTERPOLATION----------------------------------
 *
 * void getwei_linear(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
 * void getwei_TSC(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
 * int getvalues(Grid *pG, Real weight[3][3][3], int is, int js, int ks, Real *rho, Real *u1, Real *u2, Real *u3, Real *cs);
 */
/*============================================================================*/

/* get weight using linear interpolation
   Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
   Output: weight: weight function; is,js,ks: starting cell indices in the grid.
   Note: this interpolation works in any 1-3 dimensions.
*/

void getwei_linear(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c;				/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[2], wei2[2], wei3[2];	/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    i = celli(pG, x1, cell1.x1, &i1, &a);	/* x1 index */
    i1 = i+i1-1;	*is = i1;		/* starting x1 index */
    wei1[1] = a - i1 - 0.5;			/* one direction weight */
    wei1[0] = 1.0 - wei1[1];			/* 0: left; 1: right */
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    j = cellj(pG, x2, cell1.x2, &j1, &b);	/* x2 index */
    j1 = j+j1-1;	*js = j1;		/* starting x2 index */
    wei2[1] = b - j1 - 0.5;			/* one direction weight */
    wei2[0] = 1.0 - wei2[1];			/* 0: left; 1: right */
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    k = cellk(pG, x3, cell1.x3, &k1, &c);	/* x3 index */
    k1 = k+k1-1;	*ks = k1;		/* starting x3 index */
    wei3[1] = c - k1 - 0.5;			/* one direction weight */
    wei3[0] = 1.0 - wei3[1];			/* 0: left; 1: right */
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<2; k++)
    for (j=0; j<2; j++)
      for (i=0; i<2; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


/* get weight using Triangular Shaped Cloud (TSC) interpolation 
   Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
   Output: weight: weight function; is,js,ks: starting cell indices in the grid.
   Note: this interpolation works in any 1-3 dimensions.
*/
void getwei_TSC(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c, d;			/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[3], wei2[3], wei3[3];	/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
    *is = i - 1;				/* starting x1 index, wei[0] */
    d = a - i;
    wei1[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei1[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei1[2] = 0.5*SQR(d);
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;	wei1[2] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
    *js = j - 1;				/* starting x2 index */
    d = b - j;
    wei2[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei2[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei2[2] = 0.5*SQR(d);
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;	wei2[2] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
    *ks = k - 1;				/* starting x3 index */
    d = c - k;
    wei3[0] = 0.5*SQR(1.0-d);			/* 0: left; 2: right */
    wei3[1] = 0.75-SQR(d-0.5);			/* one direction weight */
    wei3[2] = 0.5*SQR(d);
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;	wei3[2] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


/* get weight using quadratic polynomial interpolation 
   Input: pG: grid; x1,x2,x3: global coordinate; cell1: 1 over dx1,dx2,dx3
   Output: weight: weight function; is,js,ks: starting cell indices in the grid.
   Note: this interpolation works in any 1-3 dimensions.
*/
void getwei_QP(Grid *pG, Real x1, Real x2, Real x3, Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c, d;			/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[3], wei2[3], wei3[3];	/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (cell1.x1 > 0.0) {
    celli(pG, x1, cell1.x1, &i, &a);		/* x1 index */
    *is = i - 1;				/* starting x1 index, wei[0] */
    d = a - i;
    wei1[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei1[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei1[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x1 dimension collapses */
    *is = pG->is;
    wei1[1] = 0.0;	wei1[2] = 0.0;
    wei1[0] = 1.0;
  }

  /* x2 direction */
  if (cell1.x2 > 0.0) {
    cellj(pG, x2, cell1.x2, &j, &b);		/* x2 index */
    *js = j - 1;				/* starting x2 index */
    d = b - j;
    wei2[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei2[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei2[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x2 dimension collapses */
    *js = pG->js;
    wei2[1] = 0.0;	wei2[2] = 0.0;
    wei2[0] = 1.0;
  }

  /* x3 direction */
  if (cell1.x3 > 0.0) {
    cellk(pG, x3, cell1.x3, &k, &c);		/* x3 index */
    *ks = k - 1;				/* starting x3 index */
    d = c - k;
    wei3[0] = 0.5*(0.5-d)*(1.5-d);		/* 0: left; 2: right */
    wei3[1] = 1.0-SQR(d-0.5);			/* one direction weight */
    wei3[2] = 0.5*(d-0.5)*(d+0.5);
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;
    wei3[1] = 0.0;	wei3[2] = 0.0;
    wei3[0] = 1.0;
  }

  /* calculate 3D weight */
  for (k=0; k<3; k++)
    for (j=0; j<3; j++)
      for (i=0; i<3; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}


/* get interpolated value using the weight
   Input:
     pG: grid; weight: weight function;
     is,js,ks: starting cell indices in the grid.
   Output:
     interpolated values of density, velocity and sound speed of the fluid
   Return: 0: normal exit;  -1: particle lie out of the grid, cannot interpolate!
   Note: this interpolation works in any 1-3 dimensions.
*/
int getvalues(Grid *pG, Real weight[3][3][3], int is, int js, int ks, Real *rho, Real *u1, Real *u2, Real *u3, Real *cs)
{
  int i, j, k, i1, j1, k1;
  Real D, v1, v2, v3;		/* density and velocity of the fluid */
#ifndef ISOTHERMAL
  Real C = 0.0;			/* fluid sound speed */
#endif
  Real totwei, totwei1;		/* total weight (in case of edge cells) */

  /* linear interpolation */
  D = 0.0; v1 = 0.0; v2 = 0.0; v3 = 0.0;
  totwei = 0.0;		totwei1 = 1.0;
  /* Interpolate density, velocity and sound speed */
  /* Note: in lower dimensions only wei[0] is non-zero, which ensures the validity */
  for (k=0; k<ncell; k++) {
    k1 = k+ks;
    if ((k1 <= kup) && (k1 >= klp)) {
      for (j=0; j<ncell; j++) {
        j1 = j+js;
        if ((j1 <= jup) && (j1 >= jlp)) {
          for (i=0; i<ncell; i++) {
            i1 = i+is;
            if ((i1 <= iup) && (i1 >= ilp)) {
              D += weight[k][j][i] * grid_d[k1][j1][i1];
              v1 += weight[k][j][i] * grid_v[k1][j1][i1].x1;
              v2 += weight[k][j][i] * grid_v[k1][j1][i1].x2;
              v3 += weight[k][j][i] * grid_v[k1][j1][i1].x3;
#ifndef ISOTHERMAL
              C += weight[k][j][i] * grid_cs[k1][j1][i1];
#endif
              totwei += weight[k][j][i];
            }
          }
        }
      }
    }
  }
  if (totwei < TINY_NUMBER) /* particle lies out of the grid, warning! */
    return -1;

  totwei1 = 1.0/totwei;
  *rho = D*totwei1;
  *u1 = v1*totwei1;	*u2 = v2*totwei1;	*u3 = v3*totwei1;
#ifdef ISOTHERMAL
  *cs = Iso_csound;
#else
  *cs = C*totwei1;
#endif /* ISOTHERMAL */

  return 0;
}


/*============================================================================*/
/*-----------------------------STOPPING TIME----------------------------------
 *
 * Real get_ts_general(Grid *pG, int type, Real rho, Real cs, Real vd);
 * Real get_ts_epstein(Grid *pG, int type, Real rho, Real cs, Real vd);
 * Real get_ts_fixed(Grid *pG, int type, Real rho, Real cs, Real vd);
 */
/*============================================================================*/

/* Calculate the stopping time for the most general case
   The relavent scale to calculate is: 
   1. a/lambda_m == alam
   2. rho_s*a in normalized unit == rhoa
*/
Real get_ts_general(Grid *pG, int type, Real rho, Real cs, Real vd)
{
  Real tstop;		/* stopping time */
  Real a, rhos;		/* primitive properties: size, solid density in cgs unit */
  Real alam;		/* a/lambda: particle size/gas mean free path */
  Real rhoa;		/* rhoa: rho_s*a in normalized unit (density, velocity, time) */
  Real Re, CD;		/* Reynolds number and drag coefficient */

  /* particle properties */
  a = pG->grproperty[type].rad;
  rhos = pG->grproperty[type].rho;
  rhoa = grrhoa[type];

  /* calculate particle size/gas mean free path */
  alam = alamcoeff * a * rho;  /* alamcoeff is defined in global.h */

  /* calculate the stopping time */
  if (alam < 2.25) {		/* Epstein regime */
    tstop = rhoa/(rho*cs);
  }

  else {
    Re = 4.0*alam*vd/cs;	/* the Reynolds number */

    if (Re < 1.) CD = 24.0/Re;	/* Stokes regime */
    else if (Re < 800.0) CD = 24.0*exp(-0.6*log(Re));
    else CD = 0.44;

    tstop = rhoa/(rho*vd*CD);
  } /* endif */

  /* avoid zero stopping time */
  if (tstop < 1.0e-8*pG->dt)
    tstop = 1.0e-8*pG->dt;

  return tstop;
}

/* Calculate the stopping time in the Epstein regime */
/* Note grrhoa == rho_s*a in normalized unit */
Real get_ts_epstein(Grid *pG, int type, Real rho, Real cs, Real vd)
{
  Real tstop = grrhoa[type]/(rho*cs);

  /* avoid zero stopping time */
  if (tstop < 1.0e-8*pG->dt)
    tstop = 1.0e-8*pG->dt;

  return tstop;
}

/* Return the fixed stopping time */
Real get_ts_fixed(Grid *pG, int type, Real rho, Real cs, Real vd)
{
  Real tstop = tstop0;

  /* avoid zero stopping time */
  if (tstop < 1.0e-8*pG->dt)
    tstop = 1.0e-8*pG->dt;

  return tstop;
}



/*============================================================================*/
/*--------------------------------FEEDBACK------------------------------------
 *
 * void get_gasinfo(Grid *pG);
 * void feedback_clear(Grid *pG);
 * void distrFB      (Grid *pG, Real weight[3][3][3], int is, int js, int ks, Vector fb);
 */
/*============================================================================*/

/* Calculate the gas information from conserved variables for feedback_predictor
   Input: pG: grid (not evolved yet).
   Output: calculate 3D array grid_v/grid_cs in the grid structure.
           Calculated are gas velocity and sound speed.
*/
void get_gasinfo(Grid *pG)
{
  int i,j,k;
  Real rho1;
#ifdef ADIABATIC
  Real P;
#endif

  /* get gas information */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++)
      {
        rho1 = 1.0/(pG->U[k][j][i].d);
        grid_d[k][j][i]    = pG->U[k][j][i].d;
        grid_v[k][j][i].x1 = pG->U[k][j][i].M1 * rho1;
        grid_v[k][j][i].x2 = pG->U[k][j][i].M2 * rho1;
        grid_v[k][j][i].x3 = pG->U[k][j][i].M3 * rho1;

#ifndef ISOTHERMAL
  #ifdef ADIABATIC
        /* E = P/(gamma-1) + rho*v^2/2 + B^2/2 */
        P = pG->U[k][j][i].E - 0.5*pG->U[k][j][i].d*(SQR(grid_v[k][j][i].x1) \
             + SQR(grid_v[k][j][i].x2) + SQR(grid_v[k][j][i].x3));
    #ifdef MHD
        P = P - 0.5*(SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
    #endif /* MHD */
        P = MAX(Gamma_1*P, TINY_NUMBER);
        grid_cs[k][j][i] = sqrt(Gamma*P*rho1);
  #else
        ath_error("[get_gasinfo] can not calculate the sound speed!\n");
  #endif /* ADIABATIC */
#endif /* ISOTHERMAL */
      }

  return;
}

#ifdef FEEDBACK

/* clean the feedback array */
void feedback_clear(Grid *pG)
{
  int i,j,k;

  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        pG->feedback[k][j][i].x1 = 0.0;
        pG->feedback[k][j][i].x2 = 0.0;
        pG->feedback[k][j][i].x3 = 0.0;
      }

  return;
}

/* Distribute the feedback force to grid cells
   Input: 
     pG: grid;   weight: weight function; 
     is,js,ks: starting cell indices in the grid.
     f1, f2, f3: feedback force from one particle.
   Output:
     pG: feedback array is updated.
*/
void distrFB(Grid *pG, Real weight[3][3][3], int is, int js, int ks, Vector fb)
{
  int i,j,k,i1,j1,k1;

  /* distribute feedback force */
  for (k=0; k<ncell; k++) {
    k1 = k+ks;
    if ((k1 <= kup) && (k1 >= klp)) {
      for (j=0; j<ncell; j++) {
        j1 = j+js;
        if ((j1 <= jup) && (j1 >= jlp)) {
          for (i=0; i<ncell; i++) {
            i1 = i+is;
            if ((i1 <= iup) && (i1 >= ilp)) {
              pG->feedback[k1][j1][i1].x1 += weight[k][j][i] * fb.x1;
              pG->feedback[k1][j1][i1].x2 += weight[k][j][i] * fb.x2;
              pG->feedback[k1][j1][i1].x3 += weight[k][j][i] * fb.x3;
            }
          }
        }
      }
    }
  }

  return;
}

#endif /* FEEDBACK */


/*============================================================================*/
/*---------------------------------SHUFFLE------------------------------------
 *
 * void shuffle(Grid *pG);
 * int compare_gr(Grid *pG, Vector cell1, Grain gr1, Grain gr2);
 * void quicksort_particle(Grid *pG, Vector cell1, long start, long end);
 */
/*============================================================================*/

/* Shuffle the particles
   Input: pG: grid with particles;
   Output: pG: particles in the linked list are rearranged by the order of their
           locations that are consistent with grid cell storage.
*/
void shuffle(Grid *pG)
{
  Grain *cur, *allgr;
  Vector cell1;

  if (pG->Nx1 > 1) cell1.x1 = 1.0/pG->dx1;  else  cell1.x1 = 0.0;
  if (pG->Nx2 > 1) cell1.x2 = 1.0/pG->dx2;  else  cell1.x2 = 0.0;
  if (pG->Nx3 > 1) cell1.x3 = 1.0/pG->dx3;  else  cell1.x3 = 0.0;

  /* output status */
  ath_pout(0, "Resorting particles...\n");

  /* sort the particles according to their positions */
  quicksort_particle(pG, cell1, 0, pG->nparticle-1);

  return;
}

/* Compare the order of two particles according to their positions in the grid
   Input: pG: grid; 
          cell1: 1/dx1,1/dx2,1/dx3, or 0 if that dimension collapses.
          gr1,gr2: pointers of the two particles to be compared.
   Output: pointer of the particle that should be put in front of the other.
*/
int compare_gr(Grid *pG, Vector cell1, Grain gr1, Grain gr2)
{
  int i1,j1,k1, i2,j2,k2;

  k1 = (int)((gr1.x3 - pG->x3_0) * cell1.x3);	/* x3 index of gr1 */
  k2 = (int)((gr2.x3 - pG->x3_0) * cell1.x3);	/* x3 index of gr2 */
  if (k1 < k2) return 1;
  if (k1 > k2) return 2;

  j1 = (int)((gr1.x2 - pG->x2_0) * cell1.x2);	/* x2 index of gr1 */
  j2 = (int)((gr2.x2 - pG->x2_0) * cell1.x2);	/* x2 index of gr2 */
  if (j1 < j2) return 1;
  if (j1 > j2) return 2;

  i1 = (int)((gr1.x1 - pG->x1_0) * cell1.x1);	/* x1 index of gr1 */
  i2 = (int)((gr2.x1 - pG->x1_0) * cell1.x1);	/* x1 index of gr2 */
  if (i1 < i2) return 1;
  if (i1 > i2) return 2;

  /* if they have equal indices, arbitrarily choose gr1 */
  return 1;
}

/* Quick sort algorithm to shuffle the particles
   Input: pG, cell1: for compare_gr subroutine only. See above.
          head, rear: head and rear of the linked list.
                      They do not contain data, or equal the pivot in the recursion.
          length: length of the linked list (does not contain head or rear).
   Output: *head: linked list with shuffling finished.
*/
void quicksort_particle(Grid *pG, Vector cell1, long start, long end)
{
  long i, pivot;
  Grain gr;
  if (end <= start) return;	/* automatically sorted already */

  /* location of the pivot at half chain length */
  pivot = (long)((start+end+1)/2);

  /* move the pivot to the start */
  gr = pG->particle[pivot];
  pG->particle[pivot] = pG->particle[start];
  pG->particle[start] = gr;

  /* initial configuration */
  pivot = start;
  i = start + 1;

  /* move the particles that are "smaller" than the pivot before it */
  while (i <= end) {
    if (compare_gr(pG, cell1, pG->particle[pivot], pG->particle[i]) == 2)
    {/* the ith particle is smaller, move it before the pivot */
      gr = pG->particle[pivot];
      pG->particle[pivot] = pG->particle[i];
      pG->particle[i] = pG->particle[pivot+1];
      pG->particle[pivot+1] = gr;
      pivot += 1;
    }
    i += 1;
  }

  /* recursively call this routine to complete sorting */
  quicksort_particle(pG, cell1, start, pivot-1);
  quicksort_particle(pG, cell1, pivot+1, end);

  return;
}


/*============================================================================*/
/*-----------------------gas velocity shift function--------------------------
 *
 * gasvshift_zero(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3);
 */
/*============================================================================*/

/* Infer new gas velocity (u1,u2,u3) based on current position and velocity.
   This is the default routine, which applies no velocity shift.
   This routine is particularly useful for code test where user can generate
   arbitrary gas velocity field, which can be done in the problem generator
   using get_usr_par_prop().
*/
void gasvshift_zero(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3)
{
  return;
}

/* User defined particle force generator.
 * This is the default routine, which has no additional forces.
 */
void User_ParticleForce_Zero(Vector *ft, Real x1, Real x2, Real x3, Real *w1, Real *w2, Real *w3)
{
  return;
}

#endif /*PARTICLES*/
