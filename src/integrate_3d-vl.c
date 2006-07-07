#include "copyright.h"
/*==============================================================================
 * FILE: integrate_3d-vl.c
 *
 * PURPOSE: 
 * Updates the input Grid structure pointed to by *pGrid by one timestep using
 * directionally unsplit van Leer scheme.
 *
 * For MHD,   uses (1 Gas + 9*Real + 9*Cons1D) 3D arrays = 80 3D arrays 
 * For hydro, uses (1 Gas + 3*Real + 9*Cons1D) 3D arrays = 53 3D arrays 
 *
 * The variables in Grid which are updated are:
 *    U.[d,M1,M2,M3,E,B1c,B2c,B3c] -- where U is of type Gas
 *    B1i, B2i, B3i  -- interface magnetic field
 *    time,dt,nstep
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static Gas ***Uhalf=NULL;
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

static ConsPotFun_t cons_pot_fun = NULL;

static void integrate_emf1_corner(const Grid *pGrid);
static void integrate_emf2_corner(const Grid *pGrid);
static void integrate_emf3_corner(const Grid *pGrid);

/*----------------------------------------------------------------------------*/
/* Function integrate_3d-vl
 */

void integrate_3d(Grid *pGrid)
{
  Real dtodx1,dtodx2,dtodx3,q1,q2,q3;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,il,iu;
  int j,jl,ju;
  int k,kl,ku;
  Real qa,cf2;
#ifdef MHD
  Real d, M1, M2, M3, B1c, B2c, B3c;
#endif

  dtodx1 = pGrid->dt/pGrid->dx1;
  il = is - 2;
  iu = ie + 2;

  dtodx2 = pGrid->dt/pGrid->dx2;
  jl = js - 2;
  ju = je + 2;

  dtodx3 = pGrid->dt/pGrid->dx3;
  kl = ks - 2;
  ku = ke + 2;

/*--- Step 1 ------------------------------------------------------------------
 * Compute first order fluxes in x1 direction
 */

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
	Ul_x1Face[k][j][i].d = pGrid->U[k][j][i-1].d;
	Ul_x1Face[k][j][i].Mx = pGrid->U[k][j][i-1].M1;
	Ul_x1Face[k][j][i].My = pGrid->U[k][j][i-1].M2;
	Ul_x1Face[k][j][i].Mz = pGrid->U[k][j][i-1].M3;
#ifndef ISOTHERMAL
	Ul_x1Face[k][j][i].E = pGrid->U[k][j][i-1].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        B1_x1Face[k][j][i] = pGrid->B1i[k][j][i];
	Ul_x1Face[k][j][i].By = pGrid->U[k][j][i-1].B2c;
	Ul_x1Face[k][j][i].Bz = pGrid->U[k][j][i-1].B3c;
#endif /* MHD */

	Ur_x1Face[k][j][i].d = pGrid->U[k][j][i].d;
	Ur_x1Face[k][j][i].Mx = pGrid->U[k][j][i].M1;
	Ur_x1Face[k][j][i].My = pGrid->U[k][j][i].M2;
	Ur_x1Face[k][j][i].Mz = pGrid->U[k][j][i].M3;
#ifndef ISOTHERMAL
	Ur_x1Face[k][j][i].E = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
	Ur_x1Face[k][j][i].By = pGrid->U[k][j][i].B2c;
	Ur_x1Face[k][j][i].Bz = pGrid->U[k][j][i].B3c;
#endif /* MHD */
      }
    }
  }

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
        GET_FLUXES(B1_x1Face[k][j][i]
          ,Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],&x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 2 ------------------------------------------------------------------
 * Compute first order fluxes in x2 direction
 */

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
	Ul_x2Face[k][j][i].d = pGrid->U[k][j-1][i].d;
	Ul_x2Face[k][j][i].Mx = pGrid->U[k][j-1][i].M2;
	Ul_x2Face[k][j][i].My = pGrid->U[k][j-1][i].M3;
	Ul_x2Face[k][j][i].Mz = pGrid->U[k][j-1][i].M1;
#ifndef ISOTHERMAL
	Ul_x2Face[k][j][i].E = pGrid->U[k][j-1][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        B2_x2Face[k][j][i] = pGrid->B2i[k][j][i];
	Ul_x2Face[k][j][i].By = pGrid->U[k][j-1][i].B3c;
	Ul_x2Face[k][j][i].Bz = pGrid->U[k][j-1][i].B1c;
#endif /* MHD */

	Ur_x2Face[k][j][i].d = pGrid->U[k][j][i].d;
	Ur_x2Face[k][j][i].Mx = pGrid->U[k][j][i].M2;
	Ur_x2Face[k][j][i].My = pGrid->U[k][j][i].M3;
	Ur_x2Face[k][j][i].Mz = pGrid->U[k][j][i].M1;
#ifndef ISOTHERMAL
	Ur_x2Face[k][j][i].E = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
	Ur_x2Face[k][j][i].By = pGrid->U[k][j][i].B3c;
	Ur_x2Face[k][j][i].Bz = pGrid->U[k][j][i].B1c;
#endif /* MHD */
      }
    }
  }

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
        GET_FLUXES(B2_x2Face[k][j][i]
          ,Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],&x2Flux[k][j][i]);
      }
    }
  }

/*--- Step 3 ------------------------------------------------------------------
 * Compute first order fluxes in x3 direction
 */

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
	Ul_x3Face[k][j][i].d = pGrid->U[k-1][j][i].d;
	Ul_x3Face[k][j][i].Mx = pGrid->U[k-1][j][i].M3;
	Ul_x3Face[k][j][i].My = pGrid->U[k-1][j][i].M1;
	Ul_x3Face[k][j][i].Mz = pGrid->U[k-1][j][i].M2;
#ifndef ISOTHERMAL
	Ul_x3Face[k][j][i].E = pGrid->U[k-1][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        B3_x3Face[k][j][i] = pGrid->B3i[k][j][i];
	Ul_x3Face[k][j][i].By = pGrid->U[k-1][j][i].B1c;
	Ul_x3Face[k][j][i].Bz = pGrid->U[k-1][j][i].B2c;
#endif /* MHD */

	Ur_x3Face[k][j][i].d = pGrid->U[k][j][i].d;
	Ur_x3Face[k][j][i].Mx = pGrid->U[k][j][i].M3;
	Ur_x3Face[k][j][i].My = pGrid->U[k][j][i].M1;
	Ur_x3Face[k][j][i].Mz = pGrid->U[k][j][i].M2;
#ifndef ISOTHERMAL
	Ur_x3Face[k][j][i].E = pGrid->U[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
	Ur_x3Face[k][j][i].By = pGrid->U[k][j][i].B1c;
	Ur_x3Face[k][j][i].Bz = pGrid->U[k][j][i].B2c;
#endif /* MHD */
      }
    }
  }

  for (k=ks-nghost+1; k<=ke+nghost; k++) {
    for (j=js-nghost+1; j<=je+nghost; j++) {
      for (i=is-nghost+1; i<=ie+nghost; i++) {
        GET_FLUXES(B3_x3Face[k][j][i]
          ,Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],&x3Flux[k][j][i]);
      }
    }
  }

/*--- Step 4 ------------------------------------------------------------------
 * Calculate the cell centered value of emf1 and integrate to corner.  Repeat
 * for emf2 and emf3
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

/* update face centered fields to dt/2 using CT */
  q1 = 0.5*dtodx1;
  q2 = 0.5*dtodx2;
  q3 = 0.5*dtodx3;
  for (k=ks-3; k<=ke+3; k++) {
    for (j=js-3; j<=je+3; j++) {
      for (i=is-3; i<=ie+3; i++) {
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
#endif /* MHD */

/*--- Step 5 ----------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x1-fluxes
 */

  q1 = 0.5*dtodx1;
  for (k=ks-nghost; k<=ke+nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost-1; i++) {
        Uhalf[k][j][i].d  = pGrid->U[k][j][i].d 
          - q1*(x1Flux[k][j][i+1].d -x1Flux[k][j][i].d );
        Uhalf[k][j][i].M1 = pGrid->U[k][j][i].M1
          - q1*(x1Flux[k][j][i+1].Mx-x1Flux[k][j][i].Mx);
        Uhalf[k][j][i].M2 = pGrid->U[k][j][i].M2
          - q1*(x1Flux[k][j][i+1].My-x1Flux[k][j][i].My);
        Uhalf[k][j][i].M3 = pGrid->U[k][j][i].M3
          - q1*(x1Flux[k][j][i+1].Mz-x1Flux[k][j][i].Mz);
#ifndef ISOTHERMAL
        Uhalf[k][j][i].E  = pGrid->U[k][j][i].E
          - q1*(x1Flux[k][j][i+1].E -x1Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        Uhalf[k][j][i].B1c= pGrid->U[k][j][i].B1c;
        Uhalf[k][j][i].B2c= pGrid->U[k][j][i].B2c
          - q1*(x1Flux[k][j][i+1].By-x1Flux[k][j][i].By);
        Uhalf[k][j][i].B3c= pGrid->U[k][j][i].B3c
          - q1*(x1Flux[k][j][i+1].Bz-x1Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }

/*--- Step 6 ----------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x2-fluxes
 */

  q2 = 0.5*dtodx2;
  for (k=ks-nghost; k<=ke+nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost-1; i++) {
        Uhalf[k][j][i].d -= q2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
        Uhalf[k][j][i].M1-= q2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
        Uhalf[k][j][i].M2-= q2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        Uhalf[k][j][i].M3-= q2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
#ifndef ISOTHERMAL
        Uhalf[k][j][i].E -= q2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        Uhalf[k][j][i].B3c-= q2*(x2Flux[k][j+1][i].By-x2Flux[k][j][i].By);
        Uhalf[k][j][i].B1c-= q2*(x2Flux[k][j+1][i].Bz-x2Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }

/*--- Step 7 ----------------------------------------------------------------
 * Update cell-centered variables to half-timestep using x3-fluxes
 */

  q3 = 0.5*dtodx3;
  for (k=ks-nghost; k<=ke+nghost-1; k++) {
    for (j=js-nghost; j<=je+nghost-1; j++) {
      for (i=is-nghost; i<=ie+nghost-1; i++) {
        Uhalf[k][j][i].d -= q3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );
        Uhalf[k][j][i].M1-= q3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        Uhalf[k][j][i].M2-= q3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
        Uhalf[k][j][i].M3-= q3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
#ifndef ISOTHERMAL
        Uhalf[k][j][i].E -= q3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
#endif /* ISOTHERMAL */
#ifdef MHD
        Uhalf[k][j][i].B1c-= q3*(x3Flux[k+1][j][i].By-x3Flux[k][j][i].By);
        Uhalf[k][j][i].B2c-= q3*(x3Flux[k+1][j][i].Bz-x3Flux[k][j][i].Bz);
#endif /* MHD */
      }
    }
  }


/*--- Step  1 ------------------------------------------------- 
 * Compute L and R states at X1-interfaces, store into 3D arrays.
 * U1d = (d, M1, M2, M3, E, B2c, B3c)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        U1d[i].d  = Uhalf[k][j][i].d;
        U1d[i].Mx = Uhalf[k][j][i].M1;
        U1d[i].My = Uhalf[k][j][i].M2;
        U1d[i].Mz = Uhalf[k][j][i].M3;
#ifndef ISOTHERMAL
        U1d[i].E  = Uhalf[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[i].By = Uhalf[k][j][i].B2c;
        U1d[i].Bz = Uhalf[k][j][i].B3c;
        Bxc[i] = Uhalf[k][j][i].B1c;
        Bxi[i] = B1_x1Face[k][j][i];
#endif /* MHD */
      }

      lr_states(U1d,Bxc,Bxi,0.0,0.0,is-2,ie+2,Ul_x1Face[k][j],Ur_x1Face[k][j]);

    }
  }
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        GET_FLUXES(B1_x1Face[k][j][i]
          ,Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],&x1Flux[k][j][i]);
      }
    }
  }

/*--- Step  2 ------------------------------------------------- 
 * Compute L and R states at X2-interfaces, store into 3D arrays.
 * U1d = (d, M1, M2, M3, E, B2c, B3c)
 */

  for (k=ks-2; k<=ke+2; k++) {
    for (i=is-2; i<=ie+2; i++) {
      for (j=js-nghost; j<=je+nghost; j++) {
        U1d[j].d  = Uhalf[k][j][i].d;
        U1d[j].Mx = Uhalf[k][j][i].M2;
        U1d[j].My = Uhalf[k][j][i].M3;
        U1d[j].Mz = Uhalf[k][j][i].M1;
#ifndef ISOTHERMAL
        U1d[j].E  = Uhalf[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[j].By = Uhalf[k][j][i].B3c;
        U1d[j].Bz = Uhalf[k][j][i].B1c;
        Bxc[j] = Uhalf[k][j][i].B2c;
        Bxi[j] = B2_x2Face[k][j][i];
#endif /* MHD */
      }

      lr_states(U1d,Bxc,Bxi,0.0,0.0,js-2,je+2,Ul,Ur);
     
      for (j=js-2; j<=je+2; j++) {
        Ul_x2Face[k][j][i] = Ul[j];
        Ur_x2Face[k][j][i] = Ur[j];
      }
    }
  }
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        GET_FLUXES(B2_x2Face[k][j][i]
          ,Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],&x2Flux[k][j][i]);
      }
    }
  }

/*--- Step x3_init_lr_states 1 ------------------------------------------------- * Compute L and R states at X3-interfaces, store into 3D arrays.
 * U1d = (d, M3, M1, M2, E, B1c, B2c)
 */

  for (j=js-2; j<=je+2; j++) {
    for (i=is-2; i<=ie+2; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
        U1d[k].d  = Uhalf[k][j][i].d;
        U1d[k].Mx = Uhalf[k][j][i].M3;
        U1d[k].My = Uhalf[k][j][i].M1;
        U1d[k].Mz = Uhalf[k][j][i].M2;
#ifndef ISOTHERMAL
        U1d[k].E  = Uhalf[k][j][i].E;
#endif /* ISOTHERMAL */
#ifdef MHD
        U1d[k].By = Uhalf[k][j][i].B1c;
        U1d[k].Bz = Uhalf[k][j][i].B2c;
        Bxc[k] = Uhalf[k][j][i].B3c;
        Bxi[k] = B3_x3Face[k][j][i];
#endif /* MHD */
      }

      lr_states(U1d,Bxc,Bxi,0.0,0.0,ks-2,ke+2,Ul,Ur);

      for (k=ks-2; k<=ke+2; k++) {
        Ul_x3Face[k][j][i] = Ul[k];
        Ur_x3Face[k][j][i] = Ur[k];
      }
    }
  }

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        GET_FLUXES(B3_x3Face[k][j][i],
           Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],&x3Flux[k][j][i]);
      }
    }
  }


/*--- Step XX ----------------------------------------------------------------
 * Calculate the cell centered value of emfs at the third-time-step, by
 * updating the cell-centered quantaties using only one flux gradient.
 */

#ifdef MHD
  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        d  = Uhalf[k][j][i].d ;
        M2 = Uhalf[k][j][i].M2;
        M3 = Uhalf[k][j][i].M3;

        B2c = 0.5*(B2_x2Face[k][j][i] + B2_x2Face[k  ][j+1][i]);
        B3c = 0.5*(B3_x3Face[k][j][i] + B3_x3Face[k+1][j  ][i]);

        emf1_cc[k][j][i] = (B2c*M3 - B3c*M2)/d;
      }
    }
  }
  integrate_emf1_corner(pGrid);

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        d  = Uhalf[k][j][i].d ;
        M1 = Uhalf[k][j][i].M1;
        M3 = Uhalf[k][j][i].M3;

        B1c = 0.5*(B1_x1Face[k][j][i] + B1_x1Face[k  ][j][i+1]);
        B3c = 0.5*(B3_x3Face[k][j][i] + B3_x3Face[k+1][j][i  ]);

        emf2_cc[k][j][i] = (B3c*M1 - B1c*M3)/d;
      }
    }
  }
  integrate_emf2_corner(pGrid);

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
        d  = Uhalf[k][j][i].d ;
        M1 = Uhalf[k][j][i].M1;
        M2 = Uhalf[k][j][i].M2;

        B1c = 0.5*(B1_x1Face[k][j][i] + B1_x1Face[k][j  ][i+1]);
        B2c = 0.5*(B2_x2Face[k][j][i] + B2_x2Face[k][j+1][i  ]);

        emf3_cc[k][j][i] = (B1c*M2 - B2c*M1)/d;
      }
    }
  }
  integrate_emf3_corner(pGrid);
#endif

/*--- Step 14(a) ---------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

#ifdef MHD
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


/*--- Step 16 ----------------------------------------------------------------
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

/*--- Step 17 -----------------------------------------------------------------
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


/*--- Step 18 ----------------------------------------------------------------
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

/*--- Step 19 -----------------------------------------------------------------
 * Update cell centered magnetic field using updated face centered fields.
 *
 * If the macro CT_CONSERVE_ENERGY is not defined (default), then internal
 * energy is conserved by subtracting off old magnetic energy, then adding
 * magnetic energy from updated face centered fields back in.  This makes the
 * scheme more robust at low beta, but violates conservation of total energy
 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#if !defined CT_CONSERVE_ENERGY && !defined ISOTHERMAL
        pGrid->U[k][j][i].E -= 0.5*( SQR(pGrid->U[k][j][i].B1c)
          + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c) );
#endif
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);

#if !defined CT_CONSERVE_ENERGY && !defined ISOTHERMAL
        pGrid->U[k][j][i].E += 0.5*( SQR(pGrid->U[k][j][i].B1c)
          + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c) );
#endif
      }
    }
  }
#endif /* MHD */

  pGrid->time += pGrid->dt;
  pGrid->nstep++;
}


/*----------------------------------------------------------------------------*/
/*  Functions integrate_emf*_corner
 *    Integrates face centered B-fluxes to compute corner EMFs.  Note:
 * x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX    */

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

/* Enroll a conservative potential function for integrating the total
   energy, including the potential energy. */
void cons_pot_fun_enroll_3d(ConsPotFun_t pfun)
{
  cons_pot_fun = pfun;
  return;
}

/*----------------------------------------------------------------------------*/
/*  Free temporary integration arrays */

void integrate_destruct_3d(void)
{
#ifdef MHD
  if (emf1 != NULL) free_3d_array((void***)emf1);
  if (emf2 != NULL) free_3d_array((void***)emf2);
  if (emf3 != NULL) free_3d_array((void***)emf3);
  if (emf1_cc != NULL) free_3d_array((void***)emf1_cc);
  if (emf2_cc != NULL) free_3d_array((void***)emf2_cc);
  if (emf3_cc != NULL) free_3d_array((void***)emf3_cc);
#endif /* MHD */
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

  if (Uhalf    != NULL) free_3d_array((void***)Uhalf);

  return;
}


/*----------------------------------------------------------------------------*/
/* Allocate temporary integration arrays */
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

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

  if ((B1_x1Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)    goto on_error;

  if ((B2_x2Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)    goto on_error;

  if ((B3_x3Face = (Real***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Real))) == NULL)    goto on_error;

  if ((U1d =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ul  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;
  if ((Ur  =      (Cons1D*)malloc(nmax*sizeof(Cons1D))) == NULL) goto on_error;

  if ((Ul_x1Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Ur_x1Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Ul_x2Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Ur_x2Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Ul_x3Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Ur_x3Face = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((x1Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((x2Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((x3Flux    = (Cons1D***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Cons1D))) == NULL)
    goto on_error;

  if ((Uhalf    = (Gas***)calloc_3d_array(Nx3,Nx2,Nx1, sizeof(Gas))) == NULL)
    goto on_error;

  return;

  on_error:
  integrate_destruct();
  ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}
