#include "copyright.h"
/*==============================================================================
 * FILE: field_loop.c
 *
 * PURPOSE: Problem generator for 2d advection of a field loop test.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *      problem/rad = radius of field loop
 *      problem/amp = amplitude of vector potential (and therefore B)
 *      problem/vflow = flow velocity
 *   The flow is automatically set to run along the diagonal. 
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   problem -  problem generator
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid)
{
   int i=0,j=0,k=0;
   int is,ie,js,je,ks,ke,nx1,nx2,nx3,iprob;
   Real x1,x2,x3;
   Real rad,amp,vflow;
   Real ***az,***ay,***ax;

   is = pGrid->is; ie = pGrid->ie;
   js = pGrid->js; je = pGrid->je;
   ks = pGrid->ks; ke = pGrid->ke;
   nx1 = (ie-is)+1 + 2*nghost;
   nx2 = (je-js)+1 + 2*nghost;
   nx3 = (ke-ks)+1 + 2*nghost;

   if (((ie-is) == 0) || ((je-js) == 0)) {
     ath_error("[field_loop]: This problem can only be run in 2D or 3D\n");
   }

   if ((ay = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
     ath_error("[field_loop]: Error allocating memory for vector pot\n");
   }
   if ((az = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
     ath_error("[field_loop]: Error allocating memory for vector pot\n");
   }
   if ((ax = (Real***)calloc_3d_array(nx3, nx2, nx1, sizeof(Real))) == NULL) {
     ath_error("[field_loop]: Error allocating memory for vector pot\n");
   }

/* Read initial conditions */

   rad = par_getd("problem","rad");
   amp = par_getd("problem","amp");
   vflow = par_getd("problem","vflow");
   iprob = par_getd("problem","iprob");

/* Use vector potential to initialize field loop in poloidal plan  */

   for (k=ks; k<=ke+1; k++) {
   for (j=js; j<=je+1; j++) {
   for (i=is; i<=ie+1; i++) {
     ax[k][j][i] = 0.0;
     ay[k][j][i] = 0.0;
     az[k][j][i] = 0.0;
     cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
     x1 -= 0.5*pGrid->dx1;
     x2 -= 0.5*pGrid->dx2;
     x3 -= 0.5*pGrid->dx3;
     
/* field loop in x1-x2 plane (cylinder in 3D) */

     if(iprob==1) {  
       if ((x1*x1 + x2*x2) < rad*rad) {
         az[k][j][i] = amp*(rad - sqrt(x1*x1 + x2*x2));
       }
     }

/* field loop in x2-x3 plane (cylinder in 3D) */

     if(iprob==2) {  
       if ((x2*x2 + x3*x3) < rad*rad) {
         ax[k][j][i] = amp*(rad - sqrt(x2*x2 + x3*x3));
       }
     }

/* field loop in x3-x1 plane (cylinder in 3D) */

     if(iprob==3) {  
       if ((x1*x1 + x3*x3) < rad*rad) {
         ay[k][j][i] = amp*(rad - sqrt(x1*x1 + x3*x3));
       }
     }

/* spherical field loop in rotated plane */

     if(iprob==4) { 
       if ((x1*x1 + x2*x2 + x3*x3) < rad*rad) {
         az[k][j][i] = amp*(rad - sqrt(x1*x1 + x2*x2 + 
                            (x3+0.5*pGrid->dx3)*(x3+0.5*pGrid->dx3)));
         ay[k][j][i] = amp*(rad - sqrt(x1*x1 + 
                       (x2+0.5*pGrid->dx2)*(x2+0.5*pGrid->dx2) + x3*x3));
       }
     }

   }}}

   x1 = pGrid->dx1*(Real)par_geti("grid","Nx1");
   x2 = pGrid->dx2*(Real)par_geti("grid","Nx2");
   x3 = pGrid->dx3*(Real)par_geti("grid","Nx3");
   rad = sqrt(x1*x1 + x2*x2 + x3*x3);
   for (k=ks; k<=ke; k++) {
   for (j=js; j<=je; j++) {
   for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = 1.0;
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vflow*x1/rad;
      pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*vflow*x2/rad;
      pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*vflow*x3/rad;
#ifdef MHD
      pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 -
                            (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
      pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 -
                            (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
      pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 -
                            (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
#endif
   }}}

/* boundary conditions on interface B */

#ifdef MHD
   i = ie+1;
   for (k=ks; k<=ke; k++) {
     for (j=js; j<=je; j++) {
       pGrid->B1i[k][j][i] = (az[k][j+1][i] - az[k][j][i])/pGrid->dx2 -
	                     (ay[k+1][j][i] - ay[k][j][i])/pGrid->dx3;
     }
   }
   j = je+1;
   for (k=ks; k<=ke; k++) {
     for (i=is; i<=ie; i++) {
       pGrid->B2i[k][j][i] = (ax[k+1][j][i] - ax[k][j][i])/pGrid->dx3 -
                             (az[k][j][i+1] - az[k][j][i])/pGrid->dx1;
     }
   }
   if (ke > ks) {
     k = ke+1;
     for (j=js; j<=je; j++) {
       for (i=is; i<=ie; i++) {
	 pGrid->B3i[k][j][i] = (ay[k][j][i+1] - ay[k][j][i])/pGrid->dx1 -
                               (ax[k][j+1][i] - ax[k][j][i])/pGrid->dx2;
       }
     }
   }
#endif

/* initialize total energy and cell-centered B */

#if defined MHD || !defined ISOTHERMAL
   for (k=ks; k<=ke; k++) {
     for (j=js; j<=je; j++) {
       for (i=is; i<=ie; i++) {
#ifdef MHD
	 pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i  ] + 
				      pGrid->B1i[k][j][i+1]);
	 pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j  ][i] +
				      pGrid->B2i[k][j+1][i]);
	 if (ke > ks)
	   pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k  ][j][i] + 
					pGrid->B3i[k+1][j][i]);
	 else
	   pGrid->U[k][j][i].B3c = pGrid->B3i[k][j][i];
#endif

#ifndef ISOTHERMAL
	 pGrid->U[k][j][i].E = 1.0/Gamma_1
#ifdef MHD
	   + 0.5*(SQR(pGrid->U[k][j][i].B1c) + SQR(pGrid->U[k][j][i].B2c)
		+ SQR(pGrid->U[k][j][i].B3c))
#endif
	   + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2)
		+ SQR(pGrid->U[k][j][i].M3))/pGrid->U[k][j][i].d;
#endif /* ISOTHERMAL */
       }
     }
   }
#endif

   free_3d_array((void***)az);
   free_3d_array((void***)ay);
   free_3d_array((void***)ax);
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * current() - computes x3-component of current
 * Bp2()     - computes magnetic pressure (Bx2 + By2)
 *----------------------------------------------------------------------------*/

void problem_write_restart(Grid *pG, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, FILE *fp)
{
  return;
}

#ifdef MHD
static Real current(const Grid *pG, const int i, const int j, const int k)
{
  return ((pG->B2i[k][j][i]-pG->B2i[k][j][i-1])/pG->dx1 - 
	  (pG->B1i[k][j][i]-pG->B1i[k][j-1][i])/pG->dx2);
}

static Real Bp2(const Grid *pG, const int i, const int j, const int k)
{
  return (pG->U[k][j][i].B1c*pG->U[k][j][i].B1c + 
	  pG->U[k][j][i].B2c*pG->U[k][j][i].B2c);
}
#endif

Gasfun_t get_usr_expr(const char *expr)
{
#ifdef MHD
  if(strcmp(expr,"J3")==0) return current;
  else if(strcmp(expr,"Bp2")==0) return Bp2;
#endif
  return NULL;
}

void Userwork_in_loop(Grid *pGrid)
{
}

void Userwork_after_loop(Grid *pGrid)
{
}
