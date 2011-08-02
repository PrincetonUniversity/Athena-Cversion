#include "copyright.h"

/*==============================================================================
 * FILE: rad2d.c
 *
 * PURPOSE:  Problem generator for a 2D test of radiative transfer routine
 *
 * Initial conditions available:
 *  iang = 1 - use Gauss-Legendre quadrature
 *  iang = 2 - use Eddington approximation
 *
 *============================================================================*/

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * gauleg()           - gauss-legendre quadrature from NR
 *============================================================================*/

static Real const_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
static Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);
Real chi0;

void problem(DomainS *pDomain)
{
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int nf=pRG->nf, nang=pRG->nang;
  int i, j, k, ifr, l, m;
  Real tau0;
  int iang, ihor1, ihor2;

/* Read problem parameters. */

  tau0 = par_getd("problem","tau0");
  iang = par_geti("problem","iang");
  ihor1 = par_geti("problem","ihor1");
  ihor2 = par_geti("problem","ihor2");
  printf("beam1: %g %g\n",pRG->mu[0][iang][0],pRG->mu[0][iang][1]);
  printf("beam2: %g %g\n",pRG->mu[1][iang][0],pRG->mu[1][iang][1]);
  chi0 = tau0 / pG->dx2;

/* Setup density structure */ 

  for (j=js-1; j<=je+1; j++) {
    for (i=is-1; i<=ie+1; i++) {
      pG->U[ks][j][i].d  = 1.0;
    }
  }

/* Initialize source function */
  for(ifr=0; ifr<nf; ifr++)
    for (j=js; j<=je+1; j++)
      for(i=is-1; i<=ie+1; i++)
	pRG->R[0][j][i][ifr].S = 0.;

/* ------- Initialize boundary emission ---------------------------------- */

  for(ifr=0; ifr<nf; ifr++) {
    for(j=pRG->js-1; j<=pRG->je+1; j++)
      for(l=0; l<4; l++) 
	for(m=0; m<nang; m++) {
	  pRG->r1imu[ifr][pRG->ks][j][l][m] = 0.0;
	  pRG->l1imu[ifr][pRG->ks][j][l][m] = 0.0;
	}
    for(i=pRG->is-1; i<=pRG->ie+1; i++)
      for(l=0; l<4; l++) 
	for(m=0; m<nang; m++) {
	  pRG->r2imu[ifr][pRG->ks][i][l][m] = 0.0;
	  pRG->l2imu[ifr][pRG->ks][i][l][m] = 0.0;
	}
    pRG->l2imu[ifr][pRG->ks][ihor1][0][iang] = 1.0;
    pRG->l2imu[ifr][pRG->ks][ihor2][1][iang] = 1.0;
  }
/* enroll radiation specification functions */
get_thermal_source = const_B;
get_thermal_fraction = const_eps;
get_total_opacity = const_opacity;

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * Userwork_in_formal_solution  - problem specific work in formal solution loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}


void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_formal_solution(DomainS *pD)
{
  return;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

static Real const_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
  return 0.0;
}

static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k)
{

  return 0.0;
  
}

static Real const_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k)
{

  return chi0;
  
}
