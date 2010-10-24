#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_2d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations in 2D using the directionally unsplit CTU
 *   method of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *      B1i, B2i -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically-thin cooling,
 *   and H-correction of Sanders et al.
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport", JCP, 205, 509 (2005)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_2d_radMHD()
 *   integrate_init_2d()
 *   integrate_destruct_2d()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(radMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS **Ul_x1Face=NULL, **Ur_x1Face=NULL;
static Cons1DS **Ul_x2Face=NULL, **Ur_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real **dhalf = NULL,**phalf = NULL;


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method in Gardiner & Stone (2005) 
 *============================================================================*/

#ifdef MHD
static void integrate_emf3_corner(GridS *pG);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_2d:  CTU integrator in 2D.
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */

void integrate_2d_radMHD(DomainS *pD)
{
	GridS *pG=(pD->Grid);
	Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;
	Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
	Real dx1 = pG->dx1, dx2 = pG->dx2;
	Real hdt = 0.5*pG->dt, dt = pG->dt;
  	int i,il,iu,is=pG->is, ie=pG->ie;
  	int j,jl,ju,js=pG->js, je=pG->je;
  	int ks=pG->ks, k, m, n;

	Real Bx = 0.0;
  	

#ifdef MHD
  	Real MHD_src,dbx,dby,B1,B2,B3,V3;
  	Real B1ch, B2ch, B3ch;
#endif

	Real temperature, velocity_x, velocity_y, pressure, density;
	Real Sigma_t, Sigma_a;
	Real SPP, alpha, Propa_44, SEE, SErho, SEmx, SEmy;

	Real Source_Inv[NVAR][NVAR], tempguess[NVAR], Uguess[NVAR], Source[NVAR], Source_guess[NVAR], Errort[NVAR];
	Real divFlux1[NVAR], divFlux2[NVAR];
	

	/* Initialize them to be zero */
	for(i=0; i<NVAR; i++){
		Source[i] = 0.0;
		Source_guess[i] = 0.0;
		Errort[i] = 0.0;
		tempguess[i] = 0.0;
		for(j=0; j<NVAR; j++) {
			Source_Inv[i][j] = 0.0;			
		if(i==j) {
		 Source_Inv[i][j] = 1.0;
		
		}
		}
	}

	il = is - 2;
  	iu = ie + 2;
  	jl = js - 2;
  	ju = je + 2;
	


/*----Step 1: Backward Euler is done in the main function for the whole mesh--------*/

	
/*=== STEP 2: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  	for (j=jl; j<=ju; j++) {
    		for (i=is-nghost; i<=ie+nghost; i++) {
      			U1d[i].d  = pG->U[ks][j][i].d;
      			U1d[i].Mx = pG->U[ks][j][i].M1;
      			U1d[i].My = pG->U[ks][j][i].M2;
      			U1d[i].Mz = pG->U[ks][j][i].M3;
      			U1d[i].E  = pG->U[ks][j][i].E;
			U1d[i].Er  = pG->U[ks][js][i].Er;
    			U1d[i].Fr1  = pG->U[ks][js][i].Fr1;
    			U1d[i].Fr2  = pG->U[ks][js][i].Fr2;
    			U1d[i].Fr3  = pG->U[ks][js][i].Fr3;
			U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;
			U1d[i].Edd_21  = pG->U[ks][js][i].Edd_21;
			U1d[i].Edd_22  = pG->U[ks][js][i].Edd_22;
#ifdef MHD
      			U1d[i].By = pG->U[ks][j][i].B2c;
      			U1d[i].Bz = pG->U[ks][j][i].B3c;
      			Bxc[i] = pG->U[ks][j][i].B1c;
      			Bxi[i] = pG->B1i[ks][j][i];
      			B1_x1Face[j][i] = pG->B1i[ks][j][i];
#endif /* MHD */
			}

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

    	for (i=is-nghost; i<=ie+nghost; i++) {
      		W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
    	}

	lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);

/*------Step 2c: Add source terms to the left and right state for 0.5*dt--------*/

	for (i=il+1; i<=iu; i++) {

	/* For left state */
		pressure = W[i-1].P;
		temperature = pressure / (U1d[i-1].d * R_ideal);
		velocity_x = U1d[i-1].Mx / U1d[i-1].d;
		velocity_y = U1d[i-1].My / U1d[i-1].d;
		Sigma_t    = pG->U[ks][j][i-1].Sigma_t;
		Sigma_a	   = pG->U[ks][j][i-1].Sigma_a;

		Source[1] = -Prat * (-Sigma_t * (U1d[i-1].Fr1/U1d[i-1].d 
			- ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y)* U1d[i-1].Er / (Crat * U1d[i-1].d))	
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - U1d[i-1].Er)/(Crat*U1d[i-1].d));
		Source[2] = -Prat * (-Sigma_t * (U1d[i-1].Fr2/U1d[i-1].d 
			- ((1.0 + U1d[i-1].Edd_22) * velocity_y + U1d[i-1].Edd_21 * velocity_x)* U1d[i-1].Er / (Crat * U1d[i-1].d))	
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - U1d[i-1].Er)/(Crat*U1d[i-1].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
			- U1d[i-1].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * (velocity_x
			* (U1d[i-1].Fr1 - ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y) * U1d[i-1].Er / Crat)
			+ velocity_y
			* (U1d[i-1].Fr2 - ((1.0 + U1d[i-1].Edd_22) * velocity_y + U1d[i-1].Edd_21 * velocity_x) * U1d[i-1].Er / Crat))/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[i-1].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
		Propa_44 = alpha;

		Wl[i].Vx += dt * Source[1] * 0.5;
		Wl[i].Vy += dt * Source[2] * 0.5;
		Wl[i].P += dt * Propa_44 * Source[4] * 0.5;

	/* For the right state */
	
	
		pressure = W[i].P;
		temperature = pressure / (U1d[i].d * R_ideal);
		velocity_x = U1d[i].Mx / U1d[i].d;
		velocity_y = U1d[i].My / U1d[i].d;
		Sigma_t    = pG->U[ks][j][i].Sigma_t;
		Sigma_a	   = pG->U[ks][j][i].Sigma_a;

		Source[1] = -Prat * (-Sigma_t * (U1d[i].Fr1/U1d[i].d 
			- ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y)* U1d[i].Er / (Crat * U1d[i].d))	
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - U1d[i].Er)/(Crat*U1d[i].d));
		Source[2] = -Prat * (-Sigma_t * (U1d[i].Fr2/U1d[i].d 
			- ((1.0 + U1d[i].Edd_22) * velocity_y + U1d[i].Edd_21 * velocity_x)* U1d[i].Er / (Crat * U1d[i].d))	
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - U1d[i].Er)/(Crat*U1d[i].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
			- U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * (velocity_x
			* (U1d[i].Fr1 - ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y) * U1d[i].Er / Crat)
			+ velocity_y
			* (U1d[i].Fr2 - ((1.0 + U1d[i].Edd_22) * velocity_y + U1d[i].Edd_21 * velocity_x) * U1d[i].Er / Crat))/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[i].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
		Propa_44 = alpha;

		Wr[i].Vx += dt * Source[1] * 0.5;
		Wr[i].Vy += dt * Source[2] * 0.5;
		Wr[i].P += dt * Propa_44 * Source[4] * 0.5;

	}
	
/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 2D array
 */

	
	 for (i=il+1; i<=iu; i++) {
     		Ul_x1Face[j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);
      		Ur_x1Face[j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

#ifdef MHD
      		Bx = B1_x1Face[j][i];
#endif

		x1Flux[j][i].d = dt;
		/* This is used to take dt to calculate alpha. 
                * x1Flux[i].d is recovered in the fluxes  */

		fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);

	} /* End to get the fluxes */

	}/* End big loop j */


/*=== STEP 3: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

	for (i=il; i<=iu; i++) {
    		for (j=js-nghost; j<=je+nghost; j++) {
      			U1d[j].d  = pG->U[ks][j][i].d;
      			U1d[j].Mx = pG->U[ks][j][i].M2;
      			U1d[j].My = pG->U[ks][j][i].M3;
      			U1d[j].Mz = pG->U[ks][j][i].M1;
      			U1d[j].E  = pG->U[ks][j][i].E;
			U1d[j].Er  = pG->U[ks][js][i].Er;
    			U1d[j].Fr1  = pG->U[ks][js][i].Fr1;
    			U1d[j].Fr2  = pG->U[ks][js][i].Fr2;
    			U1d[j].Fr3  = pG->U[ks][js][i].Fr3;
			U1d[j].Edd_11  = pG->U[ks][js][i].Edd_11;
			U1d[j].Edd_21  = pG->U[ks][js][i].Edd_21;
			U1d[j].Edd_22  = pG->U[ks][js][i].Edd_22;
#ifdef MHD
      			U1d[j].By = pG->U[ks][j][i].B3c;
      			U1d[j].Bz = pG->U[ks][j][i].B1c;
      			Bxc[j] = pG->U[ks][j][i].B2c;
      			Bxi[j] = pG->B2i[ks][j][i];
      			B2_x2Face[j][i] = pG->B2i[ks][j][i];
#endif /* MHD */
			}



/*--- Step 3b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add source terms for 0.5*dt
 */

    		for (j=js-nghost; j<=je+nghost; j++) {
      			W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
    			}

    		lr_states(pG,W,Bxc,pG->dt,dx2,jl+1,ju-1,Wl,Wr,2);


/*---------Add source terms----------------*/

	for (j=jl+1; j<=ju; j++) {

	/* For left state */
		pressure = W[j-1].P;
		temperature = pressure / (U1d[j-1].d * R_ideal);
		velocity_x = U1d[j-1].Mz / U1d[j-1].d;
		velocity_y = U1d[j-1].Mx / U1d[j-1].d;
		Sigma_t    = pG->U[ks][j-1][i].Sigma_t;
		Sigma_a	   = pG->U[ks][j-1][i].Sigma_a;

		Source[1] = -Prat * (-Sigma_t * (U1d[j-1].Fr1/U1d[j-1].d 
			- ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y)* U1d[j-1].Er / (Crat * U1d[j-1].d))	
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - U1d[j-1].Er)/(Crat*U1d[j-1].d));
		Source[2] = -Prat * (-Sigma_t * (U1d[j-1].Fr2/U1d[j-1].d 
			- ((1.0 + U1d[j-1].Edd_22) * velocity_y + U1d[j-1].Edd_21 * velocity_x)* U1d[j-1].Er / (Crat * U1d[j-1].d))	
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - U1d[j-1].Er)/(Crat*U1d[j-1].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
			- U1d[j-1].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * (velocity_x
			* (U1d[j-1].Fr1 - ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y) * U1d[j-1].Er / Crat)
			+ velocity_y
			* (U1d[j-1].Fr2 - ((1.0 + U1d[j-1].Edd_22) * velocity_y + U1d[j-1].Edd_21 * velocity_x) * U1d[j-1].Er / Crat))/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[j-1].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
		Propa_44 = alpha;

		/* "Vx" is actually vy, "vz" is actually vx; We stay with the correct meaning in source terms */
		Wl[j].Vx += dt * Source[2] * 0.5;
		Wl[j].Vz += dt * Source[1] * 0.5;
		Wl[j].P += dt * Propa_44 * Source[4] * 0.5;

	/* For the right state */
	
	
		pressure = W[j].P;
		temperature = pressure / (U1d[j].d * R_ideal);
		velocity_x = U1d[j].Mz / U1d[j].d;
		velocity_y = U1d[j].Mx / U1d[j].d;
		Sigma_t    = pG->U[ks][j][i].Sigma_t;
		Sigma_a	   = pG->U[ks][j][i].Sigma_a;

		Source[1] = -Prat * (-Sigma_t * (U1d[j].Fr1/U1d[j].d 
			- ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y)* U1d[j].Er / (Crat * U1d[j].d))	
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - U1d[j].Er)/(Crat*U1d[j].d));
		Source[2] = -Prat * (-Sigma_t * (U1d[j].Fr2/U1d[j].d 
			- ((1.0 + U1d[j].Edd_22) * velocity_y + U1d[j].Edd_21 * velocity_x)* U1d[j].Er / (Crat * U1d[j].d))	
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - U1d[j].Er)/(Crat*U1d[j].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
			- U1d[j].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * (velocity_x
			* (U1d[j].Fr1 - ((1.0 + U1d[j].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y) * U1d[j].Er / Crat)
			+ velocity_y
			* (U1d[j].Fr2 - ((1.0 + U1d[j].Edd_22) * velocity_y + U1d[j].Edd_21 * velocity_x) * U1d[j].Er / Crat))/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[j].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
		Propa_44 = alpha;

		/* "vx" is actually vy, "vy" is actually vx */
		Wr[j].Vx += dt * Source[2] * 0.5;
		Wr[j].Vz += dt * Source[1] * 0.5;
		Wr[j].P += dt * Propa_44 * Source[4] * 0.5;

	}

	for (j=jl+1; j<=ju; j++) {
      		Ul_x2Face[j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);
      		Ur_x2Face[j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);
#ifdef MHD
      		Bx = B2_x2Face[j][i];
#endif
		x2Flux[j][i].d = dt;
		/* Take the time step with x2Flux[][].d */	

      		fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],Bx,&x2Flux[j][i]);

		/* Note that x2Flux[][].Mx is actually momentum flux for vy ----
		-------------x2Flux[][].Mz is the actualy momentum flux for vx */
    	}/* End loop j to calculate the flux */

		}/*  End big loop i */
/*=== STEP 6: Correct x1-interface states with transverse flux gradients =====*/

/*--- Step 6a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS
 * The order is rotated-------------------------------------------
 */

/*--------If gravitational force is included, should correct flux due to graviational 
 * force. Gravitational force is not included in the Riemann problem But it is included 
 * in the flux as a conserved form */

/*------Radiation quantities are not calculated here-----------*/

	
  	for (j=jl+1; j<=ju-1; j++) {
    		for (i=il+1; i<=iu; i++) {

      		Ul_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i-1].d  - x2Flux[j][i-1].d );
      		Ul_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i-1].Mz - x2Flux[j][i-1].Mz);
      		Ul_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i-1].Mx - x2Flux[j][i-1].Mx);
      		Ul_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i-1].My - x2Flux[j][i-1].My);
		Ul_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i-1].E  - x2Flux[j][i-1].E );
#ifdef MHD
      		Ul_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i-1].By - x2Flux[j][i-1].By);
#endif

      		Ur_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i  ].d  - x2Flux[j][i  ].d );
      		Ur_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i  ].Mz - x2Flux[j][i  ].Mz);
      		Ur_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i  ].Mx - x2Flux[j][i  ].Mx);
      		Ur_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i  ].My - x2Flux[j][i  ].My);
		Ur_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i  ].E -  x2Flux[j][i  ].E );
#ifdef MHD
      		Ur_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i  ].By - x2Flux[j][i  ].By);
#endif

    		}
  	}


/*-------step 6b------------------------------*/
/* Correct x1 interface with x2 source gradient */
		




/*----Need extra flux gradient for self-gravity and gravitational source that is included as flux 
 * but not included in the Riemann problem---------*/

/*=== STEP 7: Correct x2-interface states with transverse flux gradients =====*/


	
  	for (j=jl+1; j<=ju; j++) {
    		for (i=il+1; i<=iu-1; i++) {

      		Ul_x2Face[j][i].d  -= hdtodx1*(x1Flux[j-1][i+1].d  - x1Flux[j-1][i].d );
      		Ul_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j-1][i+1].My - x1Flux[j-1][i].My);
      		Ul_x2Face[j][i].My -= hdtodx1*(x1Flux[j-1][i+1].Mz - x1Flux[j-1][i].Mz);
      		Ul_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j-1][i+1].Mx - x1Flux[j-1][i].Mx);
      		Ul_x2Face[j][i].E  -= hdtodx1*(x1Flux[j-1][i+1].E  - x1Flux[j-1][i].E );

#ifdef MHD
      		Ul_x2Face[j][i].By -= hdtodx1*(x1Flux[j-1][i+1].Bz - x1Flux[j-1][i].Bz);
#endif

      		Ur_x2Face[j][i].d  -= hdtodx1*(x1Flux[j  ][i+1].d  - x1Flux[j  ][i].d );
      		Ur_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j  ][i+1].My - x1Flux[j  ][i].My);
      		Ur_x2Face[j][i].My -= hdtodx1*(x1Flux[j  ][i+1].Mz - x1Flux[j  ][i].Mz);
      		Ur_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j  ][i+1].Mx - x1Flux[j  ][i].Mx);
      		Ur_x2Face[j][i].E  -= hdtodx1*(x1Flux[j  ][i+1].E  - x1Flux[j  ][i].E );

#ifdef MHD
      		Ur_x2Face[j][i].By -= hdtodx1*(x1Flux[j  ][i+1].Bz - x1Flux[j  ][i].Bz);
#endif

    }
  }





/*---------step 7b---------------*/
/* Correct x2 interface with x1 source gradient */




/*=== STEP 8: Compute 2D x1-Flux, x2-Flux ====================================*/


/*--- Step 8a ------------------------------------------------------------------
* Compute 2D x1-fluxes from corrected L/R states.
 */

	 for (j=js-1; j<=je+1; j++) {
    		for (i=is; i<=ie+1; i++) {

#ifdef MHD
      		Bx = B1_x1Face[j][i];
#endif
      		Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[j][i],&Bx);
      		Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[j][i],&Bx);

		/* Need parameter dt in radiation Riemann solver */
		x1Flux[j][i].d = dt;
	

      		fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);
    		}
  	}

	
/*--- Step 8b ------------------------------------------------------------------
 * Compute 2D x2-fluxes from corrected L/R states.
 */

  	for (j=js; j<=je+1; j++) {
    		for (i=is-1; i<=ie+1; i++) {

#ifdef MHD
      		Bx = B2_x2Face[j][i];
#endif
      		Wl[j] = Cons1D_to_Prim1D(&Ul_x2Face[j][i],&Bx);
      		Wr[j] = Cons1D_to_Prim1D(&Ur_x2Face[j][i],&Bx);

		/* take dt to calculate effective sound speed */
		x2Flux[j][i].d = dt;

      		fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],Bx,&x2Flux[j][i]);
    }
  }


/*-------Step 9: Predict and correct step after we get the flux------------- */
/*----------Radiation quantities are not updated in the modified Godunov corrector step */
	for (j=js; j<=je; j++) {
    		for (i=is; i<=ie; i++) {	
		density = pG->U[ks][j][i].d;
				
		pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 
				+ pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2) / density )	* (Gamma - 1);
		/* Should include magnetic energy for MHD */
		temperature = pressure / (density * R_ideal);
		velocity_x = pG->U[ks][j][i].M1 / density;
		velocity_y = pG->U[ks][j][i].M2 / density;
		Sigma_t    = pG->U[ks][j][i].Sigma_t;
		Sigma_a	   = pG->U[ks][j][i].Sigma_a;

	/* The Source term */
		SEE = 4.0 * Sigma_a * temperature * temperature * temperature * (Gamma - 1.0)/ (density * R_ideal);
		SErho = SEE * (-pG->U[ks][j][i].E / density + velocity_x * velocity_x + velocity_y * velocity_y);	
		SEmx = -SEE * velocity_x;
		SEmy = -SEE * velocity_y;
	
		Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][1] = -dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][2] = -dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	
		Source[1] = -Prat * (-Sigma_t * (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er)/Crat);
		Source[2] = -Prat * (-Sigma_t * (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er)/Crat);
		Source[4] = -Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er) 
			+ (Sigma_a - (Sigma_t - Sigma_a)) * velocity_x
			* (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)/Crat
			+ (Sigma_a - (Sigma_t - Sigma_a)) * velocity_y
			* (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)/Crat	); 


	/*--------Calculate the guess solution-------------*/
	
		divFlux1[0] = (x1Flux[j][i+1].d  - x1Flux[j][i].d ) / dx1;
		divFlux1[1] = (x1Flux[j][i+1].Mx - x1Flux[j][i].Mx) / dx1;
		divFlux1[2] = (x1Flux[j][i+1].My - x1Flux[j][i].My) / dx1;
		divFlux1[4] = (x1Flux[j][i+1].E  - x1Flux[j][i].E ) / dx1; 

		divFlux2[0] = (x2Flux[j+1][i].d  - x2Flux[j][i].d ) / dx2;
		divFlux2[1] = (x2Flux[j+1][i].Mz - x2Flux[j][i].Mz) / dx2;
		divFlux2[2] = (x2Flux[j+1][i].Mx - x2Flux[j][i].Mx) / dx2;
		divFlux2[4] = (x2Flux[j+1][i].E  - x2Flux[j][i].E ) / dx2; 

		for(n=0; n<5; n++) {
			tempguess[n] = 0.0;
			for(m=0; m<5; m++) {
				tempguess[n] += dt * Source_Inv[n][m] * (Source[m] - divFlux1[m] - divFlux2[m]);
			}
		}

		Uguess[0] = pG->U[ks][j][i].d + tempguess[0];
		Uguess[1] = pG->U[ks][j][i].M1 + tempguess[1];
		Uguess[2] = pG->U[ks][j][i].M2 + tempguess[2];
		Uguess[4] = pG->U[ks][j][i].E + tempguess[4];


	/*  Uguess[0] = d; Uguess[1]=Mx; Uguess[2]=My; Uguess[4]=E */

	/* Now calculate the source term due to the guess solution */

		
		density = Uguess[0];
		pressure = (Uguess[4] - 0.5 * (Uguess[1] * Uguess[1] 
				+ Uguess[2] * Uguess[2]) / density) * (Gamma - 1);
		/* Should include magnetic energy for MHD */
		temperature = pressure / (density * R_ideal);

		if(Opacity != NULL)
			Opacity(density,temperature, &Sigma_t, &Sigma_a);


		velocity_x = Uguess[1] / density;
		velocity_y = Uguess[2] / density;

	/* The Source term */
		SEE = 4.0 * Sigma_a * temperature * temperature * temperature * (Gamma - 1.0)/ (density * R_ideal);
		SErho = SEE * (-Uguess[4] / density + velocity_x * velocity_x + velocity_y * velocity_y);	
		SEmx = -SEE * velocity_x;
		SEmy = -SEE * velocity_y;
	
		Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][1] = -dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][2] = -dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE);
		Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	
		Source_guess[1] = -Prat * (-Sigma_t * (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)
			+ Sigma_a * velocity_x * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er)/Crat);
		Source_guess[2] = -Prat * (-Sigma_t * (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)
			+ Sigma_a * velocity_y * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er)/Crat);
		Source_guess[4] = -Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature - pG->U[ks][j][i].Er) 
			+ (Sigma_a - (Sigma_t - Sigma_a)) * velocity_x
			* (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)/Crat
			+ (Sigma_a - (Sigma_t - Sigma_a)) * velocity_y
			* (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)/Crat	); 

	/* Calculate the error term */
		Errort[0] = pG->U[ks][j][i].d + hdt * (Source[0] + Source_guess[0]) - dt * (divFlux1[0] + divFlux2[0]) - Uguess[0];
		Errort[1] = pG->U[ks][j][i].M1 + hdt * (Source[1] + Source_guess[1]) - dt * (divFlux1[1] + divFlux2[1]) - Uguess[1];
		Errort[2] = pG->U[ks][j][i].M2 + hdt * (Source[2] + Source_guess[2]) - dt * (divFlux1[2] + divFlux2[2]) - Uguess[2];
		Errort[4] = pG->U[ks][j][i].E + hdt * (Source[4] + Source_guess[4]) - dt * (divFlux1[4] + divFlux2[4]) - Uguess[4];

	/* Calculate the correction */
		for(n=0; n<5; n++) {
			tempguess[n] = 0.0;
			for(m=0; m<5; m++) {
				tempguess[n] += Source_Inv[n][m] * Errort[m];
			}
		}

	/* Apply the correction */
		
		pG->U[ks][j][i].d = Uguess[0] + tempguess[0];
		pG->U[ks][j][i].M1 = Uguess[1] + tempguess[1];
		pG->U[ks][j][i].M2 = Uguess[2] + tempguess[2];
		pG->U[ks][j][i].E = Uguess[4] + tempguess[4];		
	

		}/* End big loop i */
	}/* End big loop j */
		
	/* Boundary condition is applied in the main function */

	/* Update the opacity if Opacity function is set in the problem generator */
	if(Opacity != NULL){

		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++){
				
				density = pG->U[ks][j][i].d;
				
				pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 
				+ pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2) / density )	* (Gamma - 1);
			/* Should include magnetic energy for MHD */
				temperature = pressure / (density * R_ideal);
			
			Opacity(density,temperature,&(pG->U[ks][j][i].Sigma_t), &(pG->U[ks][j][i].Sigma_a));

			}
		}
	}

  return;

 
}

/*----------------------------------------------------------------------------*/
/* integrate_init_2d:    Allocate temporary integration arrays */

void integrate_init_2d(MeshS *pM)
{
  int nmax,size1=0,size2=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  nmax = MAX(size1,size2);

#ifdef MHD
  if ((emf3 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */


  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef MHD
  if ((B1_x1Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ul_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

#ifndef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  }



  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct_2d:  Free temporary integration arrays */

void integrate_destruct_2d(void)
{
#ifdef MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif /* MHD */

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_2d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_2d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_2d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_2d_array(Ur_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (dhalf     != NULL) free_2d_array(dhalf);
  if (phalf     != NULL) free_2d_array(phalf);



  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf3_corner:  */

#ifdef MHD
static void integrate_emf3_corner(GridS *pG)
{
  int i,is,ie,j,js,je;
  Real emf_l1, emf_r1, emf_l2, emf_r2;
  Real rsf=1.0,lsf=1.0;

  is = pG->is;   ie = pG->ie;
  js = pG->js;   je = pG->je;

/* NOTE: The x1-Flux of B2 is -E3.  The x2-Flux of B1 is +E3. */
  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
#ifdef CYLINDRICAL
      rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
      if (x1Flux[j-1][i].d > 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf;
      }
      else if (x1Flux[j-1][i].d < 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i].Bz - emf3_cc[j-1][i])*rsf;

      } else {
	emf_l2 = -x1Flux[j-1][i].By
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j-1][i  ])*rsf );
      }

      if (x1Flux[j][i].d > 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf;
      }
      else if (x1Flux[j][i].d < 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i].Bz - emf3_cc[j][i])*rsf;

      } else {
	emf_r2 = -x1Flux[j][i].By 
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j][i  ])*rsf );
      }

      if (x2Flux[j][i-1].d > 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]);
      }
      else if (x2Flux[j][i-1].d < 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz 
          + (-x1Flux[j][i].By - emf3_cc[j][i-1]);
      } else {
	emf_l1 = x2Flux[j][i-1].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i-1] );
      }

      if (x2Flux[j][i].d > 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i]);
      }
      else if (x2Flux[j][i].d < 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j][i].By - emf3_cc[j][i]);
      } else {
	emf_r1 = x2Flux[j][i].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i] );
      }

      emf3[j][i] = 0.25*(emf_l1 + emf_r1 + emf_l2 + emf_r2);
    }
  }

  return;
}
#endif /* MHD */

#endif /* radMHD_INTEGRATOR */
