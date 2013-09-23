#include "copyright.h"

/*==============================================================================
 * FILE: beam2d.c
 *
 * PURPOSE:  Problem generator for a simple 2D beam test of radiative transfer routine
 *
 *
 *============================================================================*/

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

Real chi0;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
#ifdef SHEARING_BOX
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
#endif
static Real const_B(const GridS *pG, const int ifr, const int i,
		    const int j, const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);
static Real const_opacity(const GridS *pG, const int ifr, const int i,
			  const int j, const int k);

void problem(DomainS *pDomain)
{
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int nf=pRG->nf, nang=pRG->nang;
  int i, j, k, ifr, l, m;
  Real tau0, t0;
  int iang, ihor1, ihor2, ivert1, ivert2;

/* Read problem parameters. */

  tau0 = par_getd("problem","tau0");
  iang = par_geti("problem","iang");
  ihor1 = par_geti("problem","ihor1");
  ihor2 = par_geti("problem","ihor2");
  ivert1 = par_geti("problem","ivert1");
  ivert2 = par_geti("problem","ivert2");
  printf("bbeam1: %g %g\n",pRG->mu[0][iang][0],pRG->mu[0][iang][1]);
  printf("bbeam2: %g %g\n",pRG->mu[1][iang][0],pRG->mu[1][iang][1]);
  chi0 = tau0 / pG->dx2;
#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","omega",1.0);
  qshear  = par_getd_def("problem","qshear",1.5);
  t0 = par_getd_def("problem","t0",0.0);
  pG->time = t0;
#endif
/* Setup density structure */ 

  R_ideal = 1.0;
  for (j=js-1; j<=je+1; j++) {
    for (i=is; i<=ie; i++) {
      pG->U[ks][j][i].d  = 1.0;
      pG->U[ks][j][i].M1 = 0.0;
      pG->U[ks][j][i].M2 = 0.0;
      pG->U[ks][j][i].M3 = 0.0;
      pG->U[ks][j][i].E = 1.0/Gamma/Gamma_1;
    }
  }
	
/* ------- Initialize boundary emission ---------------------------------- */

  for(ifr=0; ifr<nf; ifr++) {
    for(j=pRG->js-1; j<=pRG->je+1; j++)
      for(l=0; l<4; l++) 
	for(m=0; m<nang; m++) {
	  pRG->Ghstr1i[ifr][pRG->ks][j][l][m] = 0.0;
	  pRG->Ghstl1i[ifr][pRG->ks][j][l][m] = 0.0;
	}
    for(i=pRG->is-1; i<=pRG->ie+1; i++)
      for(l=0; l<4; l++) 
	for(m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][pRG->ks][i][l][m] = 0.0;
	  pRG->Ghstl2i[ifr][pRG->ks][i][l][m] = 0.0;
	}

    if(pRG->Disp[1] == 0) {
      if((ihor1 > 0) && (ihor1 >= pRG->Disp[0]) && (ihor1 < pRG->Disp[0]+pRG->Nx[0]))
	pRG->Ghstl2i[ifr][pRG->ks][ihor1-pRG->Disp[0]][0][iang] = 1.0;
      if((ihor2 > 0) && (ihor2 >= pRG->Disp[0]) && (ihor2 < pRG->Disp[0]+pRG->Nx[0]))
	pRG->Ghstl2i[ifr][pRG->ks][ihor2-pRG->Disp[0]][1][iang] = 1.0;
    }
    if(pRG->Disp[0] == 0) {
      if((ivert1 > 0) && (ivert1 >= pRG->Disp[1]) && (ivert1 < pRG->Disp[1]+pRG->Nx[1]))
	pRG->Ghstl1i[ifr][pRG->ks][ivert1-pRG->Disp[1]][0][iang] = 1.0;
      if((ivert2 > 0) && (ivert2 >= pRG->Disp[1]) && (ivert2 < pRG->Disp[1]+pRG->Nx[1]))
	pRG->Ghstl1i[ifr][pRG->ks][ivert2-pRG->Disp[1]][2][iang] = 1.0;
    }
  }

/* enroll radiation specification functions */
get_thermal_source = const_B;
get_thermal_fraction = const_eps;
get_total_opacity = const_opacity;

#ifdef SHEARING_BOX
 /* enroll gravitational potential function */
  ShearingBoxPot = UnstratifiedDisk;
#endif

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
 * Userwork_after_formal_solution  - problem specific work after formal solution
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

void Userwork_after_formal_solution(DomainS *pD)
{
  return;
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

/*=========================== PRIVATE FUNCTIONS ==============================*/


/*------------------------------------------------------------------------------
 * UnstratifiedDisk:
 */
#ifdef SHEARING_BOX
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}
#endif

static Real const_B(const GridS *pG, const RadGridS *pRG, const int ifr, const int i,
		    const int j, const int k)
{
  return 0.0;
}

static Real const_eps(const GridS *pG, const int ifr, const int i, 
		      const int j, const int k)
{

  return 0.0;
  
}

static Real const_opacity(const GridS *pG, const int ifr, const int i,
			  const int j, const int k)
{

  return chi0;
  
}
