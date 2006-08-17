#include "copyright.h"
/*==============================================================================
 * FILE: integrate_3d-ctu.c
 *
 * PURPOSE: Updates the input Grid structure pointed to by *pGrid by one
 *   timestep using directionally unsplit CTU method of Colella (1990).  The
 *   variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c] -- where U is of type Gas
 *      B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, and H-correction of Sanders et al.
 *   For adb hydro, requires (9*Cons1D +  6*Real) = 51 3D arrays
 *   For adb mhd, requires   (9*Cons1D + 12*Real) = 75 3D arrays
 *   The H-correction of Sanders et al. adds another 3 arrays.
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrined transport", JCP, ***, *** (2007)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_3d()
 *   integrate_init_3d()
 *   integrate_destruct_3d()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

#ifdef THREED_CTU  /**** This macro defined in defs.h by configure ****/

static Real *Bxc=NULL, *Bxi=NULL;
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
static Cons1D ***Ul_x1Face=NULL, ***Ur_x1Face=NULL;
static Cons1D ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
static Cons1D ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
static Cons1D *U1d=NULL, *Ul=NULL, *Ur=NULL;
static Cons1D ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
#ifdef MHD
static Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */
/* variables needed for H-correction of Sanders et al (1998) */
Real etah=0.0;
#ifdef H_CORRECTION
static Real ***eta1=NULL, ***eta2=NULL, ***eta3=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - 
 *   integrate_emf2_corner() - 
 *   integrate_emf3_corner() - 
 *============================================================================*/

static void integrate_emf1_corner(const Grid *pGrid);
static void integrate_emf2_corner(const Grid *pGrid);
static void integrate_emf3_corner(const Grid *pGrid);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_3d: 3D CTU integrator for MHD using 6-solve method */

void integrate_3d(Grid *pGrid)
{

  Real dtodx1,dtodx2,dtodx3,q1,q2,q3;
  Real dt = pGrid->dt, hdt = 0.5*pGrid->dt;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
#ifdef MHD
  Real MHD_src_By,MHD_src_Bz,mdb1,mdb2,mdb3;
  Real db1,db2,db3,l1,l2,l3,B1,B2,B3,V1,V2,V3;
  Real d, M1, M2, M3, B1c, B2c, B3c;
#endif
  Real g, x1, x2, x3;
  dtodx1 = pGrid->dt/pGrid->dx1;
  dtodx2 = pGrid->dt/pGrid->dx2;
  dtodx3 = pGrid->dt/pGrid->dx3;

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;  U1d = (d, M1, M2, M3, E, B2c, B3c)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        U1d[i].d  = pGrid->U[k][j][i].d;
        U1d[i].Mx = pGrid->U[k][j][i].M1;
        U1d[i].My = pGrid->U[k][j][i].M2;
        U1d[i].Mz = pGrid->U[k][j][i].M3;
#ifndef ISOTHERMAL
        U1d[i].E  = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[i].By = pGrid->U[k][j][i].B2c;
        U1d[i].Bz = pGrid->U[k][j][i].B3c;
        Bxc[i] = pGrid->U[k][j][i].B1c;
        Bxi[i] = pGrid->B1i[k][j][i];
        B1_x1Face[k][j][i] = pGrid->B1i[k][j][i];
#endif /* MHD */
      }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces.
 */

     lr_states(U1d,Bxc,Bxi,dt,dtodx1,is-2,ie+2,Ul_x1Face[k][j],Ur_x1Face[k][j]);

/*--- Step 1c ------------------------------------------------------------------
 * Add "MHD source terms" for 0.5*dt
 */

#ifdef MHD
      for (i=is-nghost; i<ie+nghost; i++) {
        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;

	if(db1 >= 0.0){
	  l3 = db1 < -db3 ? db1 : -db3;
	  l3 = l3 > 0.0 ? l3 : 0.0;

	  l2 = db1 < -db2 ? db1 : -db2;
	  l2 = l2 > 0.0 ? l2 : 0.0;
	}
	else{
	  l3 = db1 > -db3 ? db1 : -db3;
	  l3 = l3 < 0.0 ? l3 : 0.0;

	  l2 = db1 > -db2 ? db1 : -db2;
	  l2 = l2 < 0.0 ? l2 : 0.0;
	}

        MHD_src_By = (pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d)*l2;
        MHD_src_Bz = (pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d)*l3;

        Ul_x1Face[k][j][i].By += hdt*MHD_src_By;
        Ur_x1Face[k][j][i].By += hdt*MHD_src_By;
        Ul_x1Face[k][j][i].Bz += hdt*MHD_src_Bz;
        Ur_x1Face[k][j][i].Bz += hdt*MHD_src_Bz;
      }
#endif

/*--- Step 1d ------------------------------------------------------------------
 * Add gravitational source terms from static potential for 0.5*dt to L/R states
 */

      if (x1GravAcc != NULL){
        for (i=is-nghost; i<ie+nghost; i++) {

/* Calculate the face-centered acceleration */
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          g = (*x1GravAcc)((x1-0.5*pGrid->dx1),x2,x3);

/* Apply gravitational source terms to momentum and total energy */
          Ul_x1Face[k][j][i].Mx += hdt*Ul_x1Face[k][j][i].d*g;
          Ur_x1Face[k][j][i].Mx += hdt*Ur_x1Face[k][j][i].d*g;
#ifndef ISOTHERMAL
          Ul_x1Face[k][j][i].E += hdt*Ul_x1Face[k][j][i].Mx*g;
          Ur_x1Face[k][j][i].E += hdt*Ur_x1Face[k][j][i].Mx*g;
#endif
        }
      }
    }
  }

/*--- Step 1e ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        GET_FLUXES(B1_x1Face[k][j][i],
          Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],&x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;  U1d = (d, M2, M3, M1, E, B3c, B1c)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (i=is-2; i<=ie+2; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
        U1d[j].d  = pGrid->U[k][j][i].d;
        U1d[j].Mx = pGrid->U[k][j][i].M2;
        U1d[j].My = pGrid->U[k][j][i].M3;
        U1d[j].Mz = pGrid->U[k][j][i].M1;
#ifndef ISOTHERMAL
        U1d[j].E  = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[j].By = pGrid->U[k][j][i].B3c;
        U1d[j].Bz = pGrid->U[k][j][i].B1c;
        Bxc[j] = pGrid->U[k][j][i].B2c;
        Bxi[j] = pGrid->B2i[k][j][i];
        B2_x2Face[k][j][i] = pGrid->B2i[k][j][i];
#endif /* MHD */
      }

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces.
 */

      lr_states(U1d,Bxc,Bxi,dt,dtodx2,js-2,je+2,Ul,Ur);

/*--- Step 2c ------------------------------------------------------------------
 * Add "MHD source terms"
 */

#ifdef MHD
      for (j=js-nghost; j<je+nghost; j++) {
        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;

	if(db2 >= 0.0){
	  l1 = db2 < -db1 ? db2 : -db1;
	  l1 = l1 > 0.0 ? l1 : 0.0;

	  l3 = db2 < -db3 ? db2 : -db3;
	  l3 = l3 > 0.0 ? l3 : 0.0;
	}
	else{
	  l1 = db2 > -db1 ? db2 : -db1;
	  l1 = l1 < 0.0 ? l1 : 0.0;

	  l3 = db2 > -db3 ? db2 : -db3;
	  l3 = l3 < 0.0 ? l3 : 0.0;
	}

	MHD_src_By = (pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d)*l3;
	MHD_src_Bz = (pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d)*l1;

        Ul[j].By += hdt*MHD_src_By;
        Ur[j].By += hdt*MHD_src_By;
        Ul[j].Bz += hdt*MHD_src_Bz;
        Ur[j].Bz += hdt*MHD_src_Bz;
      }
#endif

/*--- Step 2d ------------------------------------------------------------------
 * Add gravitational source terms from static potential for 0.5*dt to L/R states
 */

      if (x2GravAcc != NULL){
        for (j=js-nghost; j<je+nghost; j++) {

/* Calculate the face-centered acceleration */
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          g = (*x2GravAcc)(x1,(x2-0.5*pGrid->dx2),x3);

/* Apply gravitational source terms to momentum and total energy */
          Ul[j].Mx += hdt*Ul[j].d*g;
          Ur[j].Mx += hdt*Ur[j].d*g;
#ifndef ISOTHERMAL
          Ul[j].E += hdt*Ul[j].Mx*g;
          Ur[j].E += hdt*Ur[j].Mx*g;
#endif
        }
      }

      for (j=js-2; j<=je+3; j++) {
        Ul_x2Face[k][j][i] = Ul[j];
        Ur_x2Face[k][j][i] = Ur[j];
      }
    }
  }

/*--- Step 2e ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 3D array
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
        GET_FLUXES(B2_x2Face[k][j][i],
           Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],&x2Flux[k][j][i]);
      }
    }
  }


/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;  U1d = (d, M3, M1, M2, E, B1c, B2c)
 */

  for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
        U1d[k].d  = pGrid->U[k][j][i].d;
        U1d[k].Mx = pGrid->U[k][j][i].M3;
        U1d[k].My = pGrid->U[k][j][i].M1;
        U1d[k].Mz = pGrid->U[k][j][i].M2;
#ifndef ISOTHERMAL
        U1d[k].E  = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[k].By = pGrid->U[k][j][i].B1c;
        U1d[k].Bz = pGrid->U[k][j][i].B2c;
        Bxc[k] = pGrid->U[k][j][i].B3c;
        Bxi[k] = pGrid->B3i[k][j][i];
        B3_x3Face[k][j][i] = pGrid->B3i[k][j][i];
#endif /* MHD */
      }

/*--- Step 3b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces.
 */

      lr_states(U1d,Bxc,Bxi,dt,dtodx3,ks-2,ke+2,Ul,Ur);

/*--- Step 3c ------------------------------------------------------------------
 * Add "MHD source terms"
 */

#ifdef MHD
      for (k=ks-nghost; k<ke+nghost; k++) {
        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;

	if(db3 >= 0.0){
	  l1 = db3 < -db1 ? db3 : -db1;
	  l1 = l1 > 0.0 ? l1 : 0.0;

	  l2 = db3 < -db2 ? db3 : -db2;
	  l2 = l2 > 0.0 ? l2 : 0.0;
	}
	else{
	  l1 = db3 > -db1 ? db3 : -db1;
	  l1 = l1 < 0.0 ? l1 : 0.0;

	  l2 = db3 > -db2 ? db3 : -db2;
	  l2 = l2 < 0.0 ? l2 : 0.0;
	}

	MHD_src_By = (pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d)*l1;
	MHD_src_Bz = (pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d)*l2;

        Ul[k].By += hdt*MHD_src_By;
        Ur[k].By += hdt*MHD_src_By;
        Ul[k].Bz += hdt*MHD_src_Bz;
        Ur[k].Bz += hdt*MHD_src_Bz;
      }
#endif

/*--- Step 3d ------------------------------------------------------------------
 * Add gravitational source terms from static potential for 0.5*dt to L/R states
 */

      if (x3GravAcc != NULL){
        for (k=ks-nghost; k<ke+nghost; k++) {

/* Calculate the face-centered acceleration */
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          g = (*x3GravAcc)(x1,x2,(x3-0.5*pGrid->dx3));

/* Apply gravitational source terms to momentum and total energy */
          Ul[k].Mx += hdt*Ul[k].d*g;
          Ur[k].Mx += hdt*Ur[k].d*g;
#ifndef ISOTHERMAL
          Ul[k].E += hdt*Ul[k].Mx*g;
          Ur[k].E += hdt*Ur[k].Mx*g;
#endif
        }
      }

      for (k=ks-2; k<=ke+3; k++) {
        Ul_x3Face[k][j][i] = Ul[k];
        Ur_x3Face[k][j][i] = Ur[k];
      }
    }
  }

/*--- Step 3e ------------------------------------------------------------------
 * Compute 1D fluxes in x3-direction, storing into 3D array
 */

  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        GET_FLUXES(B3_x3Face[k][j][i],
           Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],&x3Flux[k][j][i]);
      }
    }
  }

/*--- Step 4 -------------------------------------------------------------------
 * Calculate the cell centered value of emf1 at t^{n} and integrate to corner.
 * Repeat for emf2 and emf3
 */

#ifdef MHD
/* emf1 */
  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        emf1_cc[k][j][i] = (pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].M3 -
			    pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].M2)
                              /pGrid->U[k][j][i].d;
      }
    }
  }
  integrate_emf1_corner(pGrid);

/* emf2 */
  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        emf2_cc[k][j][i] = (pGrid->U[k][j][i].B3c*pGrid->U[k][j][i].M1 -
			    pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].M3)
                              /pGrid->U[k][j][i].d;
      }
    }
  }
  integrate_emf2_corner(pGrid);

/* emf3 */
  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
        emf3_cc[k][j][i] = (pGrid->U[k][j][i].B1c*pGrid->U[k][j][i].M2 -
			    pGrid->U[k][j][i].B2c*pGrid->U[k][j][i].M1)
                              /pGrid->U[k][j][i].d;
      }
    }
  }
  integrate_emf3_corner(pGrid);
#endif /* MHD */


/*--- Step 5 ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

#ifdef MHD
  q1 = 0.5*dtodx1;
  q2 = 0.5*dtodx2;
  q3 = 0.5*dtodx3;
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        B1_x1Face[k][j][i] += q3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                              q2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        B2_x2Face[k][j][i] += q1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                              q3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        B3_x3Face[k][j][i] += q2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                              q1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
    }
  }
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      B1_x1Face[k][j][ie+3] += q3*(emf2[k+1][j  ][ie+3]-emf2[k][j][ie+3]) -
                               q2*(emf3[k  ][j+1][ie+3]-emf3[k][j][ie+3]);
    }
  }
  for (k=ks-2; k<=ke+2; k++) {
    for (i=is-2; i<=ie+2; i++) {
      B2_x2Face[k][je+3][i] += q1*(emf3[k  ][je+3][i+1]-emf3[k][je+3][i]) -
                               q3*(emf1[k+1][je+3][i  ]-emf1[k][je+3][i]);
    }
  }
  for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      B3_x3Face[ke+3][j][i] += q2*(emf1[ke+3][j+1][i  ]-emf1[ke+3][j][i]) -
                               q1*(emf2[ke+3][j  ][i+1]-emf2[ke+3][j][i]);
    }
  }
#endif

/*--- Step 6a ------------------------------------------------------------------
 * Correct the L/R states at x1-interfaces using x2-fluxes computed in Step 2e.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

  q2 = 0.5*dtodx2;
  q3 = 0.5*dtodx3;
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        Ul_x1Face[k][j][i].d -=q2*(x2Flux[k][j+1][i-1].d -x2Flux[k][j][i-1].d );
        Ul_x1Face[k][j][i].Mx-=q2*(x2Flux[k][j+1][i-1].Mz-x2Flux[k][j][i-1].Mz);
        Ul_x1Face[k][j][i].My-=q2*(x2Flux[k][j+1][i-1].Mx-x2Flux[k][j][i-1].Mx);
        Ul_x1Face[k][j][i].Mz-=q2*(x2Flux[k][j+1][i-1].My-x2Flux[k][j][i-1].My);
#ifndef ISOTHERMAL
        Ul_x1Face[k][j][i].E -=q2*(x2Flux[k][j+1][i-1].E -x2Flux[k][j][i-1].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B3 */
	Ul_x1Face[k][j][i].Bz+=q2*0.5*
	  ((emf1[k  ][j+1][i-1] - emf1[k  ][j][i-1]) +
	   (emf1[k+1][j+1][i-1] - emf1[k+1][j][i-1]));
#endif

        Ur_x1Face[k][j][i].d -=q2*(x2Flux[k][j+1][i  ].d -x2Flux[k][j][i  ].d );
        Ur_x1Face[k][j][i].Mx-=q2*(x2Flux[k][j+1][i  ].Mz-x2Flux[k][j][i  ].Mz);
        Ur_x1Face[k][j][i].My-=q2*(x2Flux[k][j+1][i  ].Mx-x2Flux[k][j][i  ].Mx);
        Ur_x1Face[k][j][i].Mz-=q2*(x2Flux[k][j+1][i  ].My-x2Flux[k][j][i  ].My);
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E -=q2*(x2Flux[k][j+1][i  ].E -x2Flux[k][j][i  ].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B3 */
	Ur_x1Face[k][j][i].Bz+=q2*0.5*
	  ((emf1[k  ][j+1][i] - emf1[k  ][j][i]) +
	   (emf1[k+1][j+1][i] - emf1[k+1][j][i]));
#endif


/*--- Step 6b ------------------------------------------------------------------
 * Correct the L/R states at x1-interfaces using x3-fluxes computed in Step 3e.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (y,z,x) on LHS */

        Ul_x1Face[k][j][i].d -=q3*(x3Flux[k+1][j][i-1].d -x3Flux[k][j][i-1].d );
        Ul_x1Face[k][j][i].Mx-=q3*(x3Flux[k+1][j][i-1].My-x3Flux[k][j][i-1].My);
        Ul_x1Face[k][j][i].My-=q3*(x3Flux[k+1][j][i-1].Mz-x3Flux[k][j][i-1].Mz);
        Ul_x1Face[k][j][i].Mz-=q3*(x3Flux[k+1][j][i-1].Mx-x3Flux[k][j][i-1].Mx);
#ifndef ISOTHERMAL
        Ul_x1Face[k][j][i].E -=q3*(x3Flux[k+1][j][i-1].E -x3Flux[k][j][i-1].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B2 */
	Ul_x1Face[k][j][i].By-=q3*0.5*
	  ((emf1[k+1][j  ][i-1] - emf1[k][j  ][i-1]) +
	   (emf1[k+1][j+1][i-1] - emf1[k][j+1][i-1]));
#endif

        Ur_x1Face[k][j][i].d -=q3*(x3Flux[k+1][j][i  ].d -x3Flux[k][j][i  ].d );
        Ur_x1Face[k][j][i].Mx-=q3*(x3Flux[k+1][j][i  ].My-x3Flux[k][j][i  ].My);
        Ur_x1Face[k][j][i].My-=q3*(x3Flux[k+1][j][i  ].Mz-x3Flux[k][j][i  ].Mz);
        Ur_x1Face[k][j][i].Mz-=q3*(x3Flux[k+1][j][i  ].Mx-x3Flux[k][j][i  ].Mx);
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E -=q3*(x3Flux[k+1][j][i  ].E -x3Flux[k][j][i  ].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B2 */
	Ur_x1Face[k][j][i].By-=q3*0.5*
	  ((emf1[k+1][j  ][i] - emf1[k][j  ][i]) +
	   (emf1[k+1][j+1][i] - emf1[k][j+1][i]));
#endif
      }
    }
  }

/*--- Step 6c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2- and x3-flux-gradients to the
 * conservative variables on the x1Face.  Limiting is used as in GS (2007)
 */
#ifdef MHD
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        db1 = (pGrid->B1i[k ][j ][i  ] - pGrid->B1i[k][j][i-1])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i-1] - pGrid->B2i[k][j][i-1])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i-1] - pGrid->B3i[k][j][i-1])/pGrid->dx3;
        B1 = pGrid->U[k][j][i-1].B1c;
        B2 = pGrid->U[k][j][i-1].B2c;
        B3 = pGrid->U[k][j][i-1].B3c;
        V2 = pGrid->U[k][j][i-1].M2/pGrid->U[k][j][i-1].d;
        V3 = pGrid->U[k][j][i-1].M3/pGrid->U[k][j][i-1].d;

/* Calculate mdb2 = min_mod(-db1,db2) */
        if(db1 > 0.0 && db2 < 0.0){
          mdb2 = db2 > -db1 ? db2 : -db1;
        }
        else if(db1 < 0.0 && db2 > 0.0){
          mdb2 = db2 < -db1 ? db2 : -db1;
        }
        else mdb2 = 0.0;

/* Calculate mdb3 = min_mod(-db1,db3) */
        if(db1 > 0.0 && db3 < 0.0){
          mdb3 = db3 > -db1 ? db3 : -db1;
        }
        else if(db1 < 0.0 && db3 > 0.0){
          mdb3 = db3 < -db1 ? db3 : -db1;
        }
        else mdb3 = 0.0;

        Ul_x1Face[k][j][i].Mx += hdt*B1*db1;
        Ul_x1Face[k][j][i].My += hdt*B2*db1;
        Ul_x1Face[k][j][i].Mz += hdt*B3*db1;
        Ul_x1Face[k][j][i].By += hdt*V2*(-mdb3);
        Ul_x1Face[k][j][i].Bz += hdt*V3*(-mdb2);
#ifndef ISOTHERMAL
        Ul_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );
#endif /* ISOTHERMAL */

        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;
        B1 = pGrid->U[k][j][i].B1c;
        B2 = pGrid->U[k][j][i].B2c;
        B3 = pGrid->U[k][j][i].B3c;
        V2 = pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d;
        V3 = pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d;

/* Calculate mdb2 = min_mod(-db1,db2) */
        if(db1 > 0.0 && db2 < 0.0){
          mdb2 = db2 > -db1 ? db2 : -db1;
        }
        else if(db1 < 0.0 && db2 > 0.0){
          mdb2 = db2 < -db1 ? db2 : -db1;
        }
        else mdb2 = 0.0;

/* Calculate mdb3 = min_mod(-db1,db3) */
        if(db1 > 0.0 && db3 < 0.0){
          mdb3 = db3 > -db1 ? db3 : -db1;
        }
        else if(db1 < 0.0 && db3 > 0.0){
          mdb3 = db3 < -db1 ? db3 : -db1;
        }
        else mdb3 = 0.0;

        Ur_x1Face[k][j][i].Mx += hdt*B1*db1;
        Ur_x1Face[k][j][i].My += hdt*B2*db1;
        Ur_x1Face[k][j][i].Mz += hdt*B3*db1;
        Ur_x1Face[k][j][i].By += hdt*V2*(-mdb3);
        Ur_x1Face[k][j][i].Bz += hdt*V3*(-mdb2);
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );
#endif /* ISOTHERMAL */
      }
    }
  }
#endif /* MHD */

/*--- Step 6d ------------------------------------------------------------------
 * Add gravitational source terms in x2-direction and x3-direction to
 * corrected L/R states on x1-faces.
 */

  if (x2GravAcc != NULL){
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x2GravAcc)(x1,x2,x3);
        Ur_x1Face[k][j][i].My += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E += 
          0.5*hdt*(x2Flux[k][j][i  ].d + x2Flux[k][j+1][i  ].d)*g;
#endif

        g = (*x2GravAcc)((x1-pGrid->dx1),x2,x3);
        Ul_x1Face[k][j][i].My += hdt*pGrid->U[k][j][i-1].d*g;
#ifndef ISOTHERMAL
        Ul_x1Face[k][j][i].E += 
          0.5*hdt*(x2Flux[k][j][i-1].d + x2Flux[k][j+1][i-1].d)*g;
#endif
      }
    }
  }}

  if (x3GravAcc != NULL){
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x3GravAcc)(x1,x2,x3);
        Ur_x1Face[k][j][i].Mz += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E +=
          0.5*hdt*(x3Flux[k][j][i  ].d + x3Flux[k+1][j][i  ].d)*g;
#endif

        g = (*x3GravAcc)((x1-pGrid->dx1),x2,x3);
        Ul_x1Face[k][j][i].Mz += hdt*pGrid->U[k][j][i-1].d*g;
#ifndef ISOTHERMAL
        Ul_x1Face[k][j][i].E +=
          0.5*hdt*(x3Flux[k][j][i-1].d + x3Flux[k+1][j][i-1].d)*g;
#endif
      }
    }
  }}


/*--- Step 7a ------------------------------------------------------------------
 * Correct the L/R states at x2-interfaces using x1-fluxes computed in Step 1e.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */

  q1 = 0.5*dtodx1;
  q3 = 0.5*dtodx3;
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Ul_x2Face[k][j][i].d -=q1*(x1Flux[k][j-1][i+1].d -x1Flux[k][j-1][i].d );
        Ul_x2Face[k][j][i].Mx-=q1*(x1Flux[k][j-1][i+1].My-x1Flux[k][j-1][i].My);
        Ul_x2Face[k][j][i].My-=q1*(x1Flux[k][j-1][i+1].Mz-x1Flux[k][j-1][i].Mz);
        Ul_x2Face[k][j][i].Mz-=q1*(x1Flux[k][j-1][i+1].Mx-x1Flux[k][j-1][i].Mx);
#ifndef ISOTHERMAL
        Ul_x2Face[k][j][i].E -=q1*(x1Flux[k][j-1][i+1].E -x1Flux[k][j-1][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B3 */
	Ul_x2Face[k][j][i].By-=q1*0.5*
	  ((emf2[k  ][j-1][i+1] - emf2[k  ][j-1][i]) + 
	   (emf2[k+1][j-1][i+1] - emf2[k+1][j-1][i]));
#endif

        Ur_x2Face[k][j][i].d -=q1*(x1Flux[k][j  ][i+1].d -x1Flux[k][j  ][i].d );
        Ur_x2Face[k][j][i].Mx-=q1*(x1Flux[k][j  ][i+1].My-x1Flux[k][j  ][i].My);
        Ur_x2Face[k][j][i].My-=q1*(x1Flux[k][j  ][i+1].Mz-x1Flux[k][j  ][i].Mz);
        Ur_x2Face[k][j][i].Mz-=q1*(x1Flux[k][j  ][i+1].Mx-x1Flux[k][j  ][i].Mx);
#ifndef ISOTHERMAL
        Ur_x2Face[k][j][i].E -=q1*(x1Flux[k][j  ][i+1].E -x1Flux[k][j  ][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B3 */
	Ur_x2Face[k][j][i].By-=q1*0.5*
	  ((emf2[k  ][j][i+1] - emf2[k  ][j][i]) + 
	   (emf2[k+1][j][i+1] - emf2[k+1][j][i]));
#endif

/*--- Step 7b ------------------------------------------------------------------
 * Correct the L/R states at x2-interfaces using x3-fluxes computed in Step 3e.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

        Ul_x2Face[k][j][i].d -=q3*(x3Flux[k+1][j-1][i].d -x3Flux[k][j-1][i].d );
        Ul_x2Face[k][j][i].Mx-=q3*(x3Flux[k+1][j-1][i].Mz-x3Flux[k][j-1][i].Mz);
        Ul_x2Face[k][j][i].My-=q3*(x3Flux[k+1][j-1][i].Mx-x3Flux[k][j-1][i].Mx);
        Ul_x2Face[k][j][i].Mz-=q3*(x3Flux[k+1][j-1][i].My-x3Flux[k][j-1][i].My);
#ifndef ISOTHERMAL
        Ul_x2Face[k][j][i].E -=q3*(x3Flux[k+1][j-1][i].E -x3Flux[k][j-1][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B1 */
	Ul_x2Face[k][j][i].Bz+=q3*0.5*
	  ((emf2[k+1][j-1][i  ] - emf2[k][j-1][i  ]) +
	   (emf2[k+1][j-1][i+1] - emf2[k][j-1][i+1]));
#endif

        Ur_x2Face[k][j][i].d -=q3*(x3Flux[k+1][j  ][i].d -x3Flux[k][j  ][i].d );
        Ur_x2Face[k][j][i].Mx-=q3*(x3Flux[k+1][j  ][i].Mz-x3Flux[k][j  ][i].Mz);
        Ur_x2Face[k][j][i].My-=q3*(x3Flux[k+1][j  ][i].Mx-x3Flux[k][j  ][i].Mx);
        Ur_x2Face[k][j][i].Mz-=q3*(x3Flux[k+1][j  ][i].My-x3Flux[k][j  ][i].My);
#ifndef ISOTHERMAL
        Ur_x2Face[k][j][i].E -=q3*(x3Flux[k+1][j  ][i].E -x3Flux[k][j  ][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B1 */
	Ur_x2Face[k][j][i].Bz+=q3*0.5*
	  ((emf2[k+1][j][i  ] - emf2[k][j][i  ]) +
	   (emf2[k+1][j][i+1] - emf2[k][j][i+1]));
#endif
      }
    }
  }

/*--- Step 7c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x3-flux-gradients to the
 * conservative variables on the x2Face.  Limiting is used as in GS (2007)
 */
#ifdef MHD
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
        db1 = (pGrid->B1i[k][j-1][i+1] - pGrid->B1i[k][j-1][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k ][j  ][i ] - pGrid->B2i[k][j-1][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j-1][i] - pGrid->B3i[k][j-1][i])/pGrid->dx3;
        B1 = pGrid->U[k][j-1][i].B1c;
        B2 = pGrid->U[k][j-1][i].B2c;
        B3 = pGrid->U[k][j-1][i].B3c;
        V1 = pGrid->U[k][j-1][i].M1/pGrid->U[k][j-1][i].d;
        V3 = pGrid->U[k][j-1][i].M3/pGrid->U[k][j-1][i].d;

/* Calculate mdb1 = min_mod(-db2,db1) */
        if(db2 > 0.0 && db1 < 0.0){
          mdb1 = db1 > -db2 ? db1 : -db2;
        }
        else if(db2 < 0.0 && db1 > 0.0){
          mdb1 = db1 < -db2 ? db1 : -db2;
        }
        else mdb1 = 0.0;

/* Calculate mdb3 = min_mod(-db2,db3) */
        if(db2 > 0.0 && db3 < 0.0){
          mdb3 = db3 > -db2 ? db3 : -db2;
        }
        else if(db2 < 0.0 && db3 > 0.0){
          mdb3 = db3 < -db2 ? db3 : -db2;
        }
        else mdb3 = 0.0;

        Ul_x2Face[k][j][i].Mz += hdt*B1*db2;
        Ul_x2Face[k][j][i].Mx += hdt*B2*db2;
        Ul_x2Face[k][j][i].My += hdt*B3*db2;
        Ul_x2Face[k][j][i].By += hdt*V3*(-mdb1);
        Ul_x2Face[k][j][i].Bz += hdt*V1*(-mdb3);
#ifndef ISOTHERMAL
        Ul_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );
#endif /* ISOTHERMAL */

        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;
        B1 = pGrid->U[k][j][i].B1c;
        B2 = pGrid->U[k][j][i].B2c;
        B3 = pGrid->U[k][j][i].B3c;
        V1 = pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d;
        V3 = pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d;

/* Calculate mdb1 = min_mod(-db2,db1) */
        if(db2 > 0.0 && db1 < 0.0){
          mdb1 = db1 > -db2 ? db1 : -db2;
        }
        else if(db2 < 0.0 && db1 > 0.0){
          mdb1 = db1 < -db2 ? db1 : -db2;
        }
        else mdb1 = 0.0;

/* Calculate mdb3 = min_mod(-db2,db3) */
        if(db2 > 0.0 && db3 < 0.0){
          mdb3 = db3 > -db2 ? db3 : -db2;
        }
        else if(db2 < 0.0 && db3 > 0.0){
          mdb3 = db3 < -db2 ? db3 : -db2;
        }
        else mdb3 = 0.0;

        Ur_x2Face[k][j][i].Mz += hdt*B1*db2;
        Ur_x2Face[k][j][i].Mx += hdt*B2*db2;
        Ur_x2Face[k][j][i].My += hdt*B3*db2;
        Ur_x2Face[k][j][i].By += hdt*V3*(-mdb1);
        Ur_x2Face[k][j][i].Bz += hdt*V1*(-mdb3);
#ifndef ISOTHERMAL
        Ur_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );
#endif /* ISOTHERMAL */
      }
    }
  }
#endif /* MHD */

/*--- Step 7d ------------------------------------------------------------------
 * Add gravitational source terms in x1-direction and x3-direction to
 * corrected L/R states on x2-faces.
 */

  if (x1GravAcc != NULL){
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x1GravAcc)(x1,x2,x3);
        Ur_x2Face[k][j][i].Mz += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x2Face[k][j][i].E += 
          0.5*hdt*(x1Flux[k][j  ][i].d + x1Flux[k][j  ][i+1].d)*g;
#endif

        g = (*x1GravAcc)(x1,(x2-pGrid->dx2),x3);
        Ul_x2Face[k][j][i].Mz += hdt*pGrid->U[k][j-1][i].d*g;
#ifndef ISOTHERMAL
        Ul_x2Face[k][j][i].E += 
          0.5*hdt*(x1Flux[k][j-1][i].d + x1Flux[k][j-1][i+1].d)*g;
#endif
      }
    }
  }}

  if (x3GravAcc != NULL){
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x3GravAcc)(x1,x2,x3);
        Ur_x2Face[k][j][i].My += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x2Face[k][j][i].E +=
          0.5*hdt*(x3Flux[k][j  ][i].d + x3Flux[k+1][j  ][i].d)*g;
#endif

        g = (*x3GravAcc)(x1,(x2-pGrid->dx2),x3);
        Ul_x2Face[k][j][i].My += hdt*pGrid->U[k][j-1][i].d*g;
#ifndef ISOTHERMAL
        Ul_x2Face[k][j][i].E +=
          0.5*hdt*(x3Flux[k][j-1][i].d + x3Flux[k+1][j-1][i].d)*g;
#endif
      }
    }
  }}

/*--- Step 8a ------------------------------------------------------------------
 * Correct the L/R states at x3-interfaces using x1-fluxes computed in Step 1e.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

  q1 = 0.5*dtodx1;
  q2 = 0.5*dtodx2;
  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        Ul_x3Face[k][j][i].d -=q1*(x1Flux[k-1][j][i+1].d -x1Flux[k-1][j][i].d );
        Ul_x3Face[k][j][i].Mx-=q1*(x1Flux[k-1][j][i+1].Mz-x1Flux[k-1][j][i].Mz);
        Ul_x3Face[k][j][i].My-=q1*(x1Flux[k-1][j][i+1].Mx-x1Flux[k-1][j][i].Mx);
        Ul_x3Face[k][j][i].Mz-=q1*(x1Flux[k-1][j][i+1].My-x1Flux[k-1][j][i].My);
#ifndef ISOTHERMAL
        Ul_x3Face[k][j][i].E -=q1*(x1Flux[k-1][j][i+1].E -x1Flux[k-1][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B2 */
	Ul_x3Face[k][j][i].Bz+=q1*0.5*
	  ((emf3[k-1][j  ][i+1] - emf3[k-1][j  ][i]) +
	   (emf3[k-1][j+1][i+1] - emf3[k-1][j+1][i]));
#endif

        Ur_x3Face[k][j][i].d -=q1*(x1Flux[k  ][j][i+1].d -x1Flux[k  ][j][i].d );
        Ur_x3Face[k][j][i].Mx-=q1*(x1Flux[k  ][j][i+1].Mz-x1Flux[k  ][j][i].Mz);
        Ur_x3Face[k][j][i].My-=q1*(x1Flux[k  ][j][i+1].Mx-x1Flux[k  ][j][i].Mx);
        Ur_x3Face[k][j][i].Mz-=q1*(x1Flux[k  ][j][i+1].My-x1Flux[k  ][j][i].My);
#ifndef ISOTHERMAL
        Ur_x3Face[k][j][i].E -=q1*(x1Flux[k  ][j][i+1].E -x1Flux[k  ][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B2 */
	Ur_x3Face[k][j][i].Bz+=q1*0.5*
	  ((emf3[k][j  ][i+1] - emf3[k][j  ][i]) +
	   (emf3[k][j+1][i+1] - emf3[k][j+1][i]));
#endif

/*--- Step 8b ------------------------------------------------------------------
 * Correct the L/R states at x3-interfaces using x2-fluxes computed in Step 2e.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (y,z,x) on LHS 
 */

        Ul_x3Face[k][j][i].d -=q2*(x2Flux[k-1][j+1][i].d -x2Flux[k-1][j][i].d );
        Ul_x3Face[k][j][i].Mx-=q2*(x2Flux[k-1][j+1][i].My-x2Flux[k-1][j][i].My);
        Ul_x3Face[k][j][i].My-=q2*(x2Flux[k-1][j+1][i].Mz-x2Flux[k-1][j][i].Mz);
        Ul_x3Face[k][j][i].Mz-=q2*(x2Flux[k-1][j+1][i].Mx-x2Flux[k-1][j][i].Mx);
#ifndef ISOTHERMAL
        Ul_x3Face[k][j][i].E -=q2*(x2Flux[k-1][j+1][i].E -x2Flux[k-1][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B1 */
	Ul_x3Face[k][j][i].By-=q2*0.5*
	  ((emf3[k-1][j+1][i  ] - emf3[k-1][j][i  ]) +
	   (emf3[k-1][j+1][i+1] - emf3[k-1][j][i+1]));
#endif

        Ur_x3Face[k][j][i].d -=q2*(x2Flux[k  ][j+1][i].d -x2Flux[k  ][j][i].d );
        Ur_x3Face[k][j][i].Mx-=q2*(x2Flux[k  ][j+1][i].My-x2Flux[k  ][j][i].My);
        Ur_x3Face[k][j][i].My-=q2*(x2Flux[k  ][j+1][i].Mz-x2Flux[k  ][j][i].Mz);
        Ur_x3Face[k][j][i].Mz-=q2*(x2Flux[k  ][j+1][i].Mx-x2Flux[k  ][j][i].Mx);
#ifndef ISOTHERMAL
        Ur_x3Face[k][j][i].E -=q2*(x2Flux[k  ][j+1][i].E -x2Flux[k  ][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
/* Update B1 */
	Ur_x3Face[k][j][i].By-=q2*0.5*
	  ((emf3[k][j+1][i  ] - emf3[k][j][i  ]) +
	   (emf3[k][j+1][i+1] - emf3[k][j][i+1]));
#endif
      }
    }
  }

/*--- Step 8c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x2-flux-gradients to the
 * conservative variables on the x3Face.  Limiting is used as in GS (2007)
 */
#ifdef MHD
  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        db1 = (pGrid->B1i[k-1][j][i+1] - pGrid->B1i[k-1][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k-1][j+1][i] - pGrid->B2i[k-1][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k  ][j ][i ] - pGrid->B3i[k-1][j][i])/pGrid->dx3;
        B1 = pGrid->U[k-1][j][i].B1c;
        B2 = pGrid->U[k-1][j][i].B2c;
        B3 = pGrid->U[k-1][j][i].B3c;
	V1 = pGrid->U[k-1][j][i].M1/pGrid->U[k-1][j][i].d;
	V2 = pGrid->U[k-1][j][i].M2/pGrid->U[k-1][j][i].d;

/* Calculate mdb1 = min_mod(-db3,db1) */
	if(db3 > 0.0 && db1 < 0.0){
	  mdb1 = db1 > -db3 ? db1 : -db3;
	}
	else if(db3 < 0.0 && db1 > 0.0){
	  mdb1 = db1 < -db3 ? db1 : -db3;
	}
	else mdb1 = 0.0;

/* Calculate mdb2 = min_mod(-db3,db2) */
	if(db3 > 0.0 && db2 < 0.0){
	  mdb2 = db2 > -db3 ? db2 : -db3;
	}
	else if(db3 < 0.0 && db2 > 0.0){
	  mdb2 = db2 < -db3 ? db2 : -db3;
	}
	else mdb2 = 0.0;

        Ul_x3Face[k][j][i].My += hdt*B1*db3;
        Ul_x3Face[k][j][i].Mz += hdt*B2*db3;
        Ul_x3Face[k][j][i].Mx += hdt*B3*db3;
	Ul_x3Face[k][j][i].By += hdt*V1*(-mdb2);
	Ul_x3Face[k][j][i].Bz += hdt*V2*(-mdb1);
#ifndef ISOTHERMAL
	Ul_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );
#endif /* ISOTHERMAL */

        db1 = (pGrid->B1i[k][j][i+1] - pGrid->B1i[k][j][i])/pGrid->dx1;
        db2 = (pGrid->B2i[k][j+1][i] - pGrid->B2i[k][j][i])/pGrid->dx2;
        db3 = (pGrid->B3i[k+1][j][i] - pGrid->B3i[k][j][i])/pGrid->dx3;
        B1 = pGrid->U[k][j][i].B1c;
        B2 = pGrid->U[k][j][i].B2c;
        B3 = pGrid->U[k][j][i].B3c;
	V1 = pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d;
	V2 = pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d;

/* Calculate mdb1 = min_mod(-db3,db1) */
	if(db3 > 0.0 && db1 < 0.0){
	  mdb1 = db1 > -db3 ? db1 : -db3;
	}
	else if(db3 < 0.0 && db1 > 0.0){
	  mdb1 = db1 < -db3 ? db1 : -db3;
	}
	else mdb1 = 0.0;

/* Calculate mdb2 = min_mod(-db3,db2) */
	if(db3 > 0.0 && db2 < 0.0){
	  mdb2 = db2 > -db3 ? db2 : -db3;
	}
	else if(db3 < 0.0 && db2 > 0.0){
	  mdb2 = db2 < -db3 ? db2 : -db3;
	}
	else mdb2 = 0.0;

        Ur_x3Face[k][j][i].My += hdt*B1*db3;
        Ur_x3Face[k][j][i].Mz += hdt*B2*db3;
        Ur_x3Face[k][j][i].Mx += hdt*B3*db3;
	Ur_x3Face[k][j][i].By += hdt*V1*(-mdb2);
	Ur_x3Face[k][j][i].Bz += hdt*V2*(-mdb1);
#ifndef ISOTHERMAL
	Ur_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );
#endif /* ISOTHERMAL */
      }
    }
  }
#endif /* MHD */

/*--- Step 8d ------------------------------------------------------------------
 * Add gravitational source terms in x1-direction and x2-direction to
 * corrected L/R states on x3-faces.
 */

  if (x1GravAcc != NULL){
  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x1GravAcc)(x1,x2,x3);
        Ur_x3Face[k][j][i].My += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x3Face[k][j][i].E +=
          0.5*hdt*(x1Flux[k  ][j][i].d + x1Flux[k  ][j][i+1].d)*g;
#endif

        g = (*x1GravAcc)(x1,x2,(x3-pGrid->dx3));
        Ul_x3Face[k][j][i].My += hdt*pGrid->U[k-1][j][i].d*g;
#ifndef ISOTHERMAL
        Ul_x3Face[k][j][i].E +=
          0.5*hdt*(x1Flux[k-1][j][i].d + x1Flux[k-1][j][i+1].d)*g;
#endif
      }
    }
  }}

  if (x2GravAcc != NULL){
  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x2GravAcc)(x1,x2,x3);
        Ur_x3Face[k][j][i].Mz += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        Ur_x1Face[k][j][i].E +=
          0.5*hdt*(x2Flux[k  ][j][i].d + x2Flux[k  ][j+1][i].d)*g;
#endif

        g = (*x2GravAcc)(x1,x2,(x3-pGrid->dx3));
        Ul_x3Face[k][j][i].Mz += hdt*pGrid->U[k-1][j][i].d*g;
#ifndef ISOTHERMAL
        Ul_x3Face[k][j][i].E +=
          0.5*hdt*(x2Flux[k-1][j][i].d + x2Flux[k-1][j+1][i].d)*g;
#endif
      }
    }
  }}

/*--- Step 9 -------------------------------------------------------------------
 * Calculate the cell centered value of emf1 at the half-time-step and
 * integrate to corner.  Repeat for emf2 and emf3
 */

#ifdef MHD
  q1 = 0.5*dtodx1;
  q2 = 0.5*dtodx2;
  q3 = 0.5*dtodx3;

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

        d  = pGrid->U[k][j][i].d 
           - q1*(x1Flux[k  ][j  ][i+1].d - x1Flux[k][j][i].d)
           - q2*(x2Flux[k  ][j+1][i  ].d - x2Flux[k][j][i].d)
           - q3*(x3Flux[k+1][j  ][i  ].d - x3Flux[k][j][i].d);

        M1 = pGrid->U[k][j][i].M1
           - q1*(x1Flux[k  ][j  ][i+1].Mx - x1Flux[k][j][i].Mx)
           - q2*(x2Flux[k  ][j+1][i  ].Mz - x2Flux[k][j][i].Mz)
           - q3*(x3Flux[k+1][j  ][i  ].My - x3Flux[k][j][i].My);
        if (x1GravAcc != NULL) {
          g = (*x1GravAcc)(x1,x2,x3);
          M1 += hdt*pGrid->U[k][j][i].d*g;
        }

        M2 = pGrid->U[k][j][i].M2
           - q1*(x1Flux[k  ][j  ][i+1].My - x1Flux[k][j][i].My)
           - q2*(x2Flux[k  ][j+1][i  ].Mx - x2Flux[k][j][i].Mx)
           - q3*(x3Flux[k+1][j  ][i  ].Mz - x3Flux[k][j][i].Mz);
        if (x2GravAcc != NULL) {
          g = (*x2GravAcc)(x1,x2,x3);
          M2 += hdt*pGrid->U[k][j][i].d*g;
        }

        M3 = pGrid->U[k][j][i].M3
           - q1*(x1Flux[k  ][j  ][i+1].Mz - x1Flux[k][j][i].Mz)
           - q2*(x2Flux[k  ][j+1][i  ].My - x2Flux[k][j][i].My)
           - q3*(x3Flux[k+1][j  ][i  ].Mx - x3Flux[k][j][i].Mx);
        if (x3GravAcc != NULL) {
          g = (*x3GravAcc)(x1,x2,x3);
          M3 += hdt*pGrid->U[k][j][i].d*g;
        }

        B1c = 0.5*(B1_x1Face[k][j][i] + B1_x1Face[k  ][j  ][i+1]);
        B2c = 0.5*(B2_x2Face[k][j][i] + B2_x2Face[k  ][j+1][i  ]);
        B3c = 0.5*(B3_x3Face[k][j][i] + B3_x3Face[k+1][j  ][i  ]);

        emf1_cc[k][j][i] = (B2c*M3 - B3c*M2)/d;

        emf2_cc[k][j][i] = (B3c*M1 - B1c*M3)/d;

        emf3_cc[k][j][i] = (B1c*M2 - B2c*M1)/d;
      }
    }
  }
#endif

/*--- Step 10a -----------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

#ifdef H_CORRECTION
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        cfr = cfast(&(Ur_x1Face[k][j][i]), &(B1_x1Face[k][j][i]));
        cfl = cfast(&(Ul_x1Face[k][j][i]), &(B1_x1Face[k][j][i]));
        ur = Ur_x1Face[k][j][i].Mx/Ur_x1Face[k][j][i].d;
        ul = Ul_x1Face[k][j][i].Mx/Ul_x1Face[k][j][i].d;
        eta1[k][j][i] = 0.5*(fabs(ur - ul) + fabs(cfr - cfl));

        cfr = cfast(&(Ur_x2Face[k][j][i]), &(B2_x2Face[k][j][i]));
        cfl = cfast(&(Ul_x2Face[k][j][i]), &(B2_x2Face[k][j][i]));
        ur = Ur_x2Face[k][j][i].Mx/Ur_x2Face[k][j][i].d;
        ul = Ul_x2Face[k][j][i].Mx/Ul_x2Face[k][j][i].d;
        eta2[k][j][i] = 0.5*(fabs(ur - ul) + fabs(cfr - cfl));

        cfr = cfast(&(Ur_x3Face[k][j][i]), &(B3_x3Face[k][j][i]));
        cfl = cfast(&(Ul_x3Face[k][j][i]), &(B3_x3Face[k][j][i]));
        ur = Ur_x3Face[k][j][i].Mx/Ur_x3Face[k][j][i].d;
        ul = Ul_x3Face[k][j][i].Mx/Ul_x3Face[k][j][i].d;
        eta3[k][j][i] = 0.5*(fabs(ur - ul) + fabs(cfr - cfl));
      }
    }
  }
#endif /* H_CORRECTION */

/*--- Step 10b -----------------------------------------------------------------
 * Compute x1-fluxes from corrected L/R states.
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+3; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta2[k][j][i-1],eta2[k][j][i]);
        etah = MAX(etah,eta2[k][j+1][i-1]);
        etah = MAX(etah,eta2[k][j+1][i  ]);

        etah = MAX(etah,eta3[k  ][j][i-1]);
        etah = MAX(etah,eta3[k  ][j][i  ]);
        etah = MAX(etah,eta3[k+1][j][i-1]);
        etah = MAX(etah,eta3[k+1][j][i  ]);

        etah = MAX(etah,eta1[k  ][j][i  ]);
#endif /* H_CORRECTION */
        GET_FLUXES(B1_x1Face[k][j][i],
          Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],&x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 10c -----------------------------------------------------------------
 * Compute x2-fluxes from corrected L/R states.
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+3; j++) {
      for (i=is-2; i<=ie+2; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k][j-1][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k][j-1][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta3[k  ][j-1][i]);
        etah = MAX(etah,eta3[k  ][j  ][i]);
        etah = MAX(etah,eta3[k+1][j-1][i]);
        etah = MAX(etah,eta3[k+1][j  ][i]);

        etah = MAX(etah,eta2[k  ][j  ][i]);
#endif /* H_CORRECTION */
        GET_FLUXES(B2_x2Face[k][j][i],
          Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],&x2Flux[k][j][i]);
      }
    }
  }

/*--- Step 10d -----------------------------------------------------------------
 * Compute x3-fluxes from corrected L/R states.
 */

  for (k=ks-2; k<=ke+3; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k-1][j][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k-1][j][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta2[k-1][j  ][i]);
        etah = MAX(etah,eta2[k  ][j  ][i]);
        etah = MAX(etah,eta2[k-1][j+1][i]);
        etah = MAX(etah,eta2[k  ][j+1][i]);

        etah = MAX(etah,eta3[k  ][j  ][i]);
#endif /* H_CORRECTION */
        GET_FLUXES(B3_x3Face[k][j][i],
          Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],&x3Flux[k][j][i]);
      }
    }
  }

/*--- Step 11 ------------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners and then update the
 * interface magnetic fields using CT for a full time step.
 */

#ifdef MHD
  integrate_emf1_corner(pGrid);
  integrate_emf2_corner(pGrid);
  integrate_emf3_corner(pGrid);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                               dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        pGrid->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                               dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);
        pGrid->B3i[k][j][i] += dtodx2*(emf1[k  ][j+1][i  ] - emf1[k][j][i]) -
                               dtodx1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      }
    }
  }
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] +=
        dtodx3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
        dtodx2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    }
  }
  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      pGrid->B2i[k][je+1][i] +=
        dtodx1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
        dtodx3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    }
  }
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pGrid->B3i[ke+1][j][i] += 
        dtodx2*(emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
        dtodx1*(emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    }
  }
#endif

/*--- Step 12 ------------------------------------------------------------------
 * To keep the gravitational source terms 2nd order, add 0.5 the gravitational
 * acceleration to the momentum equation now (using d^{n}), before the update
 * of the cell-centered variables due to flux gradients.
 */

  if (x1GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x1GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M1 += hdt*pGrid->U[k][j][i].d*g;
      }
    }
  }}

  if (x2GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x2GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M2 += hdt*pGrid->U[k][j][i].d*g;
      }
    }
  }}

  if (x3GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x3GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M3 += hdt*pGrid->U[k][j][i].d*g;
      }
    }
  }}

/*--- Step 13a -----------------------------------------------------------------
 * Update cell-centered variables in pGrid using x1-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d -=dtodx1*(x1Flux[k][j][i+1].d -x1Flux[k][j][i].d );
        pGrid->U[k][j][i].M1-=dtodx1*(x1Flux[k][j][i+1].Mx-x1Flux[k][j][i].Mx);
        pGrid->U[k][j][i].M2-=dtodx1*(x1Flux[k][j][i+1].My-x1Flux[k][j][i].My);
        pGrid->U[k][j][i].M3-=dtodx1*(x1Flux[k][j][i+1].Mz-x1Flux[k][j][i].Mz);
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E -=dtodx1*(x1Flux[k][j][i+1].E -x1Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        pGrid->U[k][j][i].B2c-=dtodx1*(x1Flux[k][j][i+1].By-x1Flux[k][j][i].By);
        pGrid->U[k][j][i].B3c-=dtodx1*(x1Flux[k][j][i+1].Bz-x1Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }

/*--- Step 13b -----------------------------------------------------------------
 * Update cell-centered variables in pGrid using x2-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d -=dtodx2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
        pGrid->U[k][j][i].M1-=dtodx2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
        pGrid->U[k][j][i].M2-=dtodx2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        pGrid->U[k][j][i].M3-=dtodx2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E -=dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        pGrid->U[k][j][i].B3c-=dtodx2*(x2Flux[k][j+1][i].By-x2Flux[k][j][i].By);
        pGrid->U[k][j][i].B1c-=dtodx2*(x2Flux[k][j+1][i].Bz-x2Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }

/*--- Step 13c -----------------------------------------------------------------
 * Update cell-centered variables in pGrid using x3-fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].d -=dtodx3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );
        pGrid->U[k][j][i].M1-=dtodx3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        pGrid->U[k][j][i].M2-=dtodx3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
        pGrid->U[k][j][i].M3-=dtodx3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E -=dtodx3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        pGrid->U[k][j][i].B1c-=dtodx3*(x3Flux[k+1][j][i].By-x3Flux[k][j][i].By);
        pGrid->U[k][j][i].B2c-=dtodx3*(x3Flux[k+1][j][i].Bz-x3Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }

/*--- Step 14 ------------------------------------------------------------------
 * Complete the gravitational source terms by adding 0.5 the acceleration at
 * time level n+1, and the energy source term at time level {n+1/2}.
 */

  if (x1GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x1GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M1 += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += hdt*(x1Flux[k][j][i].d + x1Flux[k][j][i+1].d)*g;
#endif
      }
    }
  }}

  if (x2GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x2GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M2 += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += hdt*(x2Flux[k][j][i].d + x2Flux[k][j+1][i].d)*g;
#endif
      }
    }
  }}

  if (x3GravAcc != NULL){
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        g = (*x3GravAcc)(x1,x2,x3);

        pGrid->U[k][j][i].M3 += hdt*pGrid->U[k][j][i].d*g;
#ifndef ISOTHERMAL
        pGrid->U[k][j][i].E += hdt*(x3Flux[k][j][i].d + x3Flux[k+1][j][i].d)*g;
#endif
      }
    }
  }}

/*--- Step 15 ------------------------------------------------------------------
 * LAST STEP!
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

  return;
}


/*----------------------------------------------------------------------------*/
/*  Free temporary integration arrays 
 */

void integrate_destruct_3d(void)
{

#ifdef MHD
  if (emf1    != NULL) free_3d_array((void***)emf1);
  if (emf2    != NULL) free_3d_array((void***)emf2);
  if (emf3    != NULL) free_3d_array((void***)emf3);
  if (emf1_cc != NULL) free_3d_array((void***)emf1_cc);
  if (emf2_cc != NULL) free_3d_array((void***)emf2_cc);
  if (emf3_cc != NULL) free_3d_array((void***)emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_3d_array((void***)eta1);
  if (eta2 != NULL) free_3d_array((void***)eta2);
  if (eta3 != NULL) free_3d_array((void***)eta3);
#endif /* H_CORRECTION */

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
  if (B1_x1Face != NULL) free_3d_array((void***)B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array((void***)B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array((void***)B3_x3Face);

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);

  if (Ul_x1Face != NULL) free_3d_array((void***)Ul_x1Face);
  if (Ur_x1Face != NULL) free_3d_array((void***)Ur_x1Face);
  if (Ul_x2Face != NULL) free_3d_array((void***)Ul_x2Face);
  if (Ur_x2Face != NULL) free_3d_array((void***)Ur_x2Face);
  if (Ul_x3Face != NULL) free_3d_array((void***)Ul_x3Face);
  if (Ur_x3Face != NULL) free_3d_array((void***)Ur_x3Face);
  if (x1Flux    != NULL) free_3d_array((void***)x1Flux);
  if (x2Flux    != NULL) free_3d_array((void***)x2Flux);
  if (x3Flux    != NULL) free_3d_array((void***)x3Flux);

  return;
}

/*----------------------------------------------------------------------------*/
/* Allocate temporary integration arrays 
*/

void integrate_init_3d(int nx1, int nx2, int nx3)
{
  int nmax;
  int Nx1 = nx1 + 2*nghost;
  int Nx2 = nx2 + 2*nghost;
  int Nx3 = nx3 + 2*nghost;
  nmax = MAX(MAX(Nx1,Nx2),Nx3);

#ifdef MHD
  if ((emf1 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;

  if ((emf1_cc = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf2_cc = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((eta2 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((eta3 = (Real***)calloc_3d_array(Nx3, Nx2, Nx1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* H_CORRECTION */

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

  if ((B1_x1Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)
    goto on_error;

  if ((U1d =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;

  if ((Ul_x1Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D)))
    == NULL) goto on_error;
  if ((Ur_x1Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D)))
    == NULL) goto on_error;
  if ((Ul_x2Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D)))
    == NULL) goto on_error;
  if ((Ur_x2Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D)))
    == NULL) goto on_error;
  if ((Ul_x3Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;
  if ((Ur_x3Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;
  if ((x1Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;
  if ((x2Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;
  if ((x3Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) 
    == NULL) goto on_error;

  return;

  on_error:
  integrate_destruct();
  ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf1_corner
 * integrate_emf2_corner
 * integrate_emf3_corner
 *   Integrates face centered B-fluxes to compute corner EMFs.  Note:
 *   x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *   x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *   x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *   x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *   x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *   x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */

#ifdef MHD
static void integrate_emf1_corner(const Grid *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}


static void integrate_emf2_corner(const Grid *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k-1][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

static void integrate_emf3_corner(const Grid *pGrid)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  Real de3_l1, de3_r1, de3_l2, de3_r2;

  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1];
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i];
	else {
	  de3_l2 = 0.5*(x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1] + 
			x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = x2Flux[k][j][i].Bz - emf3_cc[k][j][i];
	else {
	  de3_r2 = 0.5*(x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1] + 
			x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ] );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}
#endif /* MHD */
#endif /* THREED_CTU */
