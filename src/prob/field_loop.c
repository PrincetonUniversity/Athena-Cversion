#include "copyright.h"
/*==============================================================================
 * FILE: field_loop.c
 *
 * PURPOSE: Problem generator for advection of a field loop test.  Can only
 *   be run in 2D or 3D.  Input parameters are:
 *      problem/rad   = radius of field loop
 *      problem/amp   = amplitude of vector potential (and therefore B)
 *      problem/vflow = flow velocity
 *   The flow is automatically set to run along the diagonal. 
 *   Various test cases are possible:
 *     (iprob=1): field loop in x1-x2 plane (cylinder in 3D)
 *     (iprob=2): field loop in x2-x3 plane (cylinder in 3D)
 *     (iprob=3): field loop in x3-x1 plane (cylinder in 3D) 
 *     (iprob=4): rotated cylindrical field loop in 3D.
 *     (iprob=5): spherical field loop in rotated plane
 *
 * REFERENCE: T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD
 *   via constrined transport", JCP, 205, 509 (2005)
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(Grid *pGrid, Domain *pDomain)
{
   int i=0,j=0,k=0;
   int is,ie,js,je,ks,ke,nx1,nx2,nx3,iprob;
   Real x1c,x2c,x3c,x1f,x2f,x3f;       /* cell- and face-centered coordinates */
   Real x1size,x3size,lambda=0.0,ang_2=0.0,sin_a2=0.0,cos_a2=1.0,x,y;
   Real rad,amp,vflow;
   Real ***az,***ay,***ax;

   is = pGrid->is; ie = pGrid->ie;
   js = pGrid->js; je = pGrid->je;
   ks = pGrid->ks; ke = pGrid->ke;
   nx1 = (ie-is)+1 + 2*nghost;
   nx2 = (je-js)+1 + 2*nghost;
   nx3 = (ke-ks)+1 + 2*nghost;

   if (((je-js) == 0)) {
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

/* For (iprob=4) -- rotated cylinder in 3D -- set up rotation angle and
 * wavelength of cylinder */

  if(iprob == 4){
     x1size = par_getd("grid","x1max") - par_getd("grid","x1min");
     x3size = par_getd("grid","x3max") - par_getd("grid","x3min");

/* We put 1 wavelength in each direction.  Hence the wavelength
 *     lambda = x1size*cos_a;
 *     AND   lambda = x3size*sin_a;
 *     are both satisfied. */

     if(x1size == x3size){
       ang_2 = PI/4.0;
       cos_a2 = sin_a2 = sqrt(0.5);
     }
     else{
       ang_2 = atan(x1size/x3size);
       sin_a2 = sin(ang_2);
       cos_a2 = cos(ang_2);
     }
/* Use the larger angle to determine the wavelength */
     if (cos_a2 >= sin_a2) {
       lambda = x1size*cos_a2;
     } else {
       lambda = x3size*sin_a2;
     }
   }

/* Use vector potential to initialize field loop */

   for (k=ks; k<=ke+1; k++) {
   for (j=js; j<=je+1; j++) {
   for (i=is; i<=ie+1; i++) {
     cc_pos(pGrid,i,j,k,&x1c,&x2c,&x3c);
     x1f = x1c - 0.5*pGrid->dx1;
     x2f = x2c - 0.5*pGrid->dx2;
     x3f = x3c - 0.5*pGrid->dx3;
     
/* (iprob=1): field loop in x1-x2 plane (cylinder in 3D) */

     if(iprob==1) {  
       ax[k][j][i] = 0.0;
       ay[k][j][i] = 0.0;
       if ((x1f*x1f + x2f*x2f) < rad*rad) {
         az[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2f*x2f));
       }
     }

/* (iprob=2): field loop in x2-x3 plane (cylinder in 3D) */

     if(iprob==2) {  
       if ((x2f*x2f + x3f*x3f) < rad*rad) {
         ax[k][j][i] = amp*(rad - sqrt(x2f*x2f + x3f*x3f));
       }
       ay[k][j][i] = 0.0;
       az[k][j][i] = 0.0;
     }

/* (iprob=3): field loop in x3-x1 plane (cylinder in 3D) */

     if(iprob==3) {  
       ax[k][j][i] = 0.0;
       if ((x1f*x1f + x3f*x3f) < rad*rad) {
         ay[k][j][i] = amp*(rad - sqrt(x1f*x1f + x3f*x3f));
       }
       az[k][j][i] = 0.0;
     }

/* (iprob=4): rotated cylindrical field loop in 3D.  Similar to iprob=1
 * with a rotation about the x2-axis.  Define coordinate systems (x1,x2,x3)
 * and (x,y,z) with the following transformation rules:
 *    x =  x1*cos(ang_2) + x3*sin(ang_2)
 *    y =  x2
 *    z = -x1*sin(ang_2) + x3*cos(ang_2)
 * This inverts to:
 *    x1  = x*cos(ang_2) - z*sin(ang_2)
 *    x2  = y
 *    x3  = x*sin(ang_2) + z*cos(ang_2)
 */

     if(iprob==4) {
       x = x1c*cos_a2 + x3f*sin_a2;
       y = x2f;
/* shift x back to the domain -0.5*lambda <= x <= 0.5*lambda */
       while(x >  0.5*lambda) x -= lambda;
       while(x < -0.5*lambda) x += lambda;
       if ((x*x + y*y) < rad*rad) {
         ax[k][j][i] = amp*(rad - sqrt(x*x + y*y))*(-sin_a2);
       }

       ay[k][j][i] = 0.0;

       x = x1f*cos_a2 + x3c*sin_a2;
       y = x2f;
/* shift x back to the domain -0.5*lambda <= x <= 0.5*lambda */
       while(x >  0.5*lambda) x -= lambda;
       while(x < -0.5*lambda) x += lambda;
       if ((x*x + y*y) < rad*rad) {
         az[k][j][i] = amp*(rad - sqrt(x*x + y*y))*(cos_a2);
       }
     }

/* (iprob=5): spherical field loop in rotated plane */

     if(iprob==5) { 
       ax[k][j][i] = 0.0;
       if ((x1f*x1f + x2c*x2c + x3f*x3f) < rad*rad) {
         ay[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2c*x2c + x3f*x3f));
       }
       if ((x1f*x1f + x2f*x2f + x3c*x3c) < rad*rad) {
         az[k][j][i] = amp*(rad - sqrt(x1f*x1f + x2f*x2f + x3c*x3c));
       }
     }

   }}}

   x1c = pGrid->dx1*(Real)par_geti("grid","Nx1");
   x2c = pGrid->dx2*(Real)par_geti("grid","Nx2");
   x3c = pGrid->dx3*(Real)par_geti("grid","Nx3");
   rad = sqrt(x1c*x1c + x2c*x2c + x3c*x3c);
   for (k=ks; k<=ke; k++) {
   for (j=js; j<=je; j++) {
   for (i=is; i<=ie; i++) {
      pGrid->U[k][j][i].d = 1.0;
      pGrid->U[k][j][i].M1 = pGrid->U[k][j][i].d*vflow*x1c/rad;
      pGrid->U[k][j][i].M2 = pGrid->U[k][j][i].d*vflow*x2c/rad;
      pGrid->U[k][j][i].M3 = pGrid->U[k][j][i].d*vflow*x3c/rad;
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

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
}
