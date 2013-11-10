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

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS **Ul_x1Face=NULL, **Ur_x1Face=NULL;
static Cons1DS **Ul_x2Face=NULL, **Ur_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;

#ifdef CONS_GRAVITY
static Real **x1Flux_grav=NULL;
static Real **x2Flux_grav=NULL;
static Real **density_old=NULL;
Real dotphil, dotgxl;
#endif

/* The interface magnetic fields and emfs */
#ifdef RADIATION_MHD
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

#ifdef RADIATION_MHD
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
  	int ks=pG->ks, m, n;

	int badcellflag = 0;

/* For static gravitational potential */
	Real x1, x2, x3, phicl, phicr, phifc, phic, phir, phil;

#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flux_m1l,flux_m1r,flux_m2l,flux_m2r;
#ifdef CONS_GRAVITY
  Real Tempswap;
#endif

#endif

	Real Bx = 0.0;
	Real M1h, M2h, M3h;
  	

#ifdef RADIATION_MHD
  	Real MHD_src,dbx,dby,B1,B2,B3,V3;
  	Real B1ch, B2ch, B3ch;
#endif
#ifdef SHEARING_BOX
	
	/* For radiation MHD code, we always assume shearing box is  in x-y plane */
/*	if(ShBoxCoord != xy)
		ath_error("[integrator_2d_radMHD]: Shearing box in radiation MHD code must be in x-y plane");	
*/	
	/* in XY shearing sheet 2=phi; in XZ shearing sheet 2=Z and 3=phi */
	Real Vphi, Mphi;
	Real M1n, dM2n, dM3n;   /* M1, dM2/3=(My+d*q*Omega_0*x) at time n */
	Real M1e, dM2e, dM3e;   /* M1, dM2/3 evolved by dt/2  */
	Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2;
	Real flx1_dM3, frx1_dM3, flx2_dM3, frx2_dM3;
	
/*	Real Source_M1h, Source_M2h, Source_M3h;
	Real ShearingSource_M1, ShearingSource_M2, ShearingSource_M3, ShearingSource_E;
	Real Shearingguess_M1, Shearingguess_M2, Shearingguess_M3, Shearingguess_E;
*/
	Real qom, om_dt = Omega_0*pG->dt;
	qom = qshear * Omega_0;
	Real fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
	
	Real ShearSource[4];

#endif /* SHEARING_BOX */
	
	

	Real temperature, velocity_x, velocity_y, velocity_z, pressure, density, Tguess, Fr0x, Fr0y, diffTEr, velocity, Smx0, Smy0;
	Real Prwork1, Prwork2, Prworksource, Ersource;
	Real Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE;
	Real SPP, alpha, Propa_44, SEE, SErho, SEmx, SEmy;
/*	Real dSigma[2*NOPACITY], dSigmadP[NOPACITY];
*/ 
	Real Sigma[NOPACITY];
	Cons1DS Usource;	
	/* for source term */
	PrimS Wopacity; /* temporary variable for opacity function */

	/* In case momentum becomes stiff */
	Real SFmx, SFmy, SVVx, SVVy, betax, betay;


	Real Source_Inv[NVAR][NVAR], Source_Invnew[NVAR][NVAR], tempguess[NVAR], Uguess[NVAR], Source[NVAR], Source_guess[NVAR], Errort[NVAR], SourceFlux[NVAR];
	Real divFlux1[NVAR], divFlux2[NVAR];
	
	/* SourceFlux is used to calculate flux due to other non-stiff source terms, such as static gravitational potential */

	/* Initialize them to be zero */
	for(i=0; i<NVAR; i++){
		Source[i] = 0.0;
		SourceFlux[i] = 0.0;
		Source_guess[i] = 0.0;
		Errort[i] = 0.0;
		divFlux1[i] = 0.0;
		divFlux2[i] = 0.0;
		tempguess[i] = 0.0;
		for(j=0; j<NVAR; j++) {
			Source_Inv[i][j] = 0.0;			
		if(i==j) {
		 Source_Inv[i][j] = 1.0;
		
		}
		}
	}
/*	for(i=0; i<2*NOPACITY; i++)
		dSigma[i] = 0.0;
*/

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
			U1d[i].Er  = pG->U[ks][j][i].Er;
    			U1d[i].Fr1  = pG->U[ks][j][i].Fr1;
    			U1d[i].Fr2  = pG->U[ks][j][i].Fr2;
    			U1d[i].Fr3  = pG->U[ks][j][i].Fr3;
			U1d[i].Edd_11  = pG->U[ks][j][i].Edd_11;
			U1d[i].Edd_21  = pG->U[ks][j][i].Edd_21;
			U1d[i].Edd_22  = pG->U[ks][j][i].Edd_22;
			
			for(m=0; m<NOPACITY;m++){
				U1d[i].Sigma[m] = pG->U[ks][j][i].Sigma[m];
			}
		
#ifdef RADIATION_MHD
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

	lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,2);

/*------Step 2c: Add source terms to the left and right state for 0.5*dt--------*/


	for (i=il+1; i<=iu; i++) {

	/* For left state */
		pressure = W[i-1].P;
		temperature = pressure / (U1d[i-1].d * R_ideal);
/*
		Tguess = pG->Tguess[ks][j][i-1];
*/
		Tguess = temperature;

		velocity_x = U1d[i-1].Mx / U1d[i-1].d;
		velocity_y = U1d[i-1].My / U1d[i-1].d;
		
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i-1,j,ks,&x1,&x2,&x3);
		velocity_y -= qom * x1;		
#endif

		velocity = velocity_x * velocity_x + velocity_y * velocity_y;
		
		
		Sigma_sF = U1d[i-1].Sigma[0];
		Sigma_aF = U1d[i-1].Sigma[1];
		Sigma_aP = U1d[i-1].Sigma[2];
		Sigma_aE = U1d[i-1].Sigma[3];


		Source[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[i-1].Fr1/U1d[i-1].d 
			- ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y)* U1d[i-1].Er / (Crat * U1d[i-1].d))	
			+ velocity_x * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[i-1].Er)/(Crat*U1d[i-1].d));
		Source[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[i-1].Fr2/U1d[i-1].d 
			- ((1.0 + U1d[i-1].Edd_22) * velocity_y + U1d[i-1].Edd_21 * velocity_x)* U1d[i-1].Er / (Crat * U1d[i-1].d))	
			+ velocity_y * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[i-1].Er)/(Crat*U1d[i-1].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_aP * (pow(Tguess, 4.0)
			- U1d[i-1].Er) + (Sigma_aF - Sigma_sF) * (velocity_x
			* (U1d[i-1].Fr1 - ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y) * U1d[i-1].Er / Crat)
			+ velocity_y
			* (U1d[i-1].Fr2 - ((1.0 + U1d[i-1].Edd_22) * velocity_y + U1d[i-1].Edd_21 * velocity_x) * U1d[i-1].Er / Crat))/Crat)
			- (Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2]) * U1d[i-1].d; 

/*
		if(Opacity != NULL) Opacity(U1d[i-1].d, temperature, NULL, dSigma);
		
		dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
		dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
		dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
		dSigmadP[3] =  dSigma[7] / (density * R_ideal); 
*/
		/*
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature /(U1d[i-1].d * R_ideal)
			-(Gamma - 1.0) * Prat * Crat * (dSigmadP[2] * pow(Tguess, 4.0) - dSigmadP[3] * U1d[i-1].Er)
		      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP[1] * (
			velocity_x * (U1d[i-1].Fr1 - ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y) * U1d[i-1].Er/Crat)
			+ velocity_y * (U1d[i-1].Fr2 - (U1d[i-1].Edd_21 * velocity_x + (1.0 + U1d[i-1].Edd_22) * velocity_y) * U1d[i-1].Er/Crat)
			);
		*/
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature * (1.0 - velocity/(Crat * Crat)) /(U1d[i-1].d * R_ideal);
		

		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
		SVVx = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[i-1].Edd_11) * W[i-1].Er / (W[i-1].d * Crat);
		
		if(fabs(SVVx * dt * 0.5) > 0.001)
		betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
		else 
		betax = 1.0 + 0.25 * SVVx * dt;

		SVVy = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[i-1].Edd_22) * W[i-1].Er / (W[i-1].d * Crat);
		
		if(fabs(SVVy * dt * 0.5) > 0.001)
		betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
		else 
		betay = 1.0 + 0.25 * SVVy * dt;
		/*===========================================================================*/
	
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
		Propa_44 = alpha;

        if(Erflag){
            Wl[i].Vx += dt * Source[1] * 0.5 * betax;
            Wl[i].Vy += dt * Source[2] * 0.5 * betay;
            Wl[i].P += dt * Propa_44 * Source[4] * 0.5;
        }

		if(Wl[i].P < TINY_NUMBER)
			Wl[i].P -= dt * Propa_44 * Source[4] * 0.5;
	
		for(m=0; m<NOPACITY; m++){
			Wl[i].Sigma[m] = U1d[i-1].Sigma[m];
		}


		Wl[i].Edd_11 = W[i-1].Edd_11;
		Wl[i].Edd_21 = W[i-1].Edd_21;
		Wl[i].Edd_22 = W[i-1].Edd_22;
		Wl[i].Edd_31 = W[i-1].Edd_31;
		Wl[i].Edd_32 = W[i-1].Edd_32;
		Wl[i].Edd_33 = W[i-1].Edd_33;

	/* For the right state */
	
	
		pressure = W[i].P;
		temperature = pressure / (U1d[i].d * R_ideal);
/*
		Tguess = pG->Tguess[ks][j][i];
*/
		Tguess = temperature;

		velocity_x = U1d[i].Mx / U1d[i].d;
		velocity_y = U1d[i].My / U1d[i].d;
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		velocity_y -= qom * x1;		
#endif
		
		velocity = velocity_x * velocity_x + velocity_y * velocity_y;
		
		Sigma_sF = U1d[i].Sigma[0];
		Sigma_aF = U1d[i].Sigma[1];
		Sigma_aP = U1d[i].Sigma[2];
		Sigma_aE = U1d[i].Sigma[3];

		Source[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[i].Fr1/U1d[i].d 
			- ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y)* U1d[i].Er / (Crat * U1d[i].d))	
			+ velocity_x * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[i].Er)/(Crat*U1d[i].d));
		Source[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[i].Fr2/U1d[i].d 
			- ((1.0 + U1d[i].Edd_22) * velocity_y + U1d[i].Edd_21 * velocity_x)* U1d[i].Er / (Crat * U1d[i].d))	
			+ velocity_y * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[i].Er)/(Crat*U1d[i].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * ((Sigma_aP * pow(Tguess, 4.0) 
			- Sigma_aE * U1d[i].Er) + (Sigma_aF - Sigma_sF) * (velocity_x
			* (U1d[i].Fr1 - ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y) * U1d[i].Er / Crat)
			+ velocity_y
			* (U1d[i].Fr2 - ((1.0 + U1d[i].Edd_22) * velocity_y + U1d[i].Edd_21 * velocity_x) * U1d[i].Er / Crat))/Crat)
			- (Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2]) * U1d[i].d; 

			
	/*	if(Opacity != NULL) Opacity(U1d[i].d, temperature, NULL, dSigma);
		
		dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
		dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
		dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
		dSigmadP[3] =  dSigma[7] / (density * R_ideal); 
	*/	
		/*
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature	* temperature /(U1d[i].d * R_ideal)
			-(Gamma - 1.0) * Prat * Crat * (dSigmadP[2] * pow(Tguess, 4.0) - dSigmadP[3] * U1d[i].Er)
		      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP[1] * (
			velocity_x * (U1d[i].Fr1 - ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y) * U1d[i].Er/Crat)
			+ velocity_y * (U1d[i].Fr2 - (U1d[i].Edd_21 * velocity_x + (1.0 + U1d[i].Edd_22) * velocity_y) * U1d[i].Er/Crat)
			);
		*/
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature	* temperature * (1.0 - velocity/(Crat * Crat))/(U1d[i].d * R_ideal);
		

		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
		SVVx = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[i].Edd_11) * W[i].Er / (W[i].d * Crat);
		
		if(fabs(SVVx * dt * 0.5) > 0.001)
		betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
		else 
		betax = 1.0 + 0.25 * SVVx * dt;

		SVVy = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[i].Edd_22) * W[i].Er / (W[i].d * Crat);
		
		if(fabs(SVVy * dt * 0.5) > 0.001)
		betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
		else 
		betay = 1.0 + 0.25 * SVVy * dt;
		/*===========================================================================*/
	


		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
		Propa_44 = alpha;
        
        if(Erflag){

            Wr[i].Vx += dt * Source[1] * 0.5 * betax;
            Wr[i].Vy += dt * Source[2] * 0.5 * betay;
            Wr[i].P += dt * Propa_44 * Source[4] * 0.5;

            if(Wr[i].P < TINY_NUMBER)
                Wr[i].P -= dt * Propa_44 * Source[4] * 0.5;
            
        }

		for(m=0; m<NOPACITY; m++){
			Wr[i].Sigma[m] = U1d[i].Sigma[m];
		}
	
		Wr[i].Edd_11 = W[i].Edd_11;
		Wr[i].Edd_21 = W[i].Edd_21;
		Wr[i].Edd_22 = W[i].Edd_22;
		Wr[i].Edd_31 = W[i].Edd_31;
		Wr[i].Edd_32 = W[i].Edd_32;
		Wr[i].Edd_33 = W[i].Edd_33;		

	}

/* For radiation_mhd, source term due to magnetic field part is also added */
#ifdef RADIATION_MHD
	for(i=il+1;i<=iu;i++){
		MHD_src = (pG->U[ks][j][i-1].M2/pG->U[ks][j][i-1].d)*
               		(pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1])/pG->dx1;
      		Wl[i].By += hdt*MHD_src;

		 MHD_src = (pG->U[ks][j][i].M2/pG->U[ks][j][i].d)*
               		(pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i])/pG->dx1;
      		Wr[i].By += hdt*MHD_src;
	}
#endif /* for radMHd */

	
 /* Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
	 if (StaticGravPot != NULL){
      		for (i=il+1; i<=iu; i++) {
        	cc_pos(pG,i,j,ks,&x1,&x2,&x3);

	        phicr = (*StaticGravPot)( x1             ,x2,x3);
        	phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
        	phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

        	Wl[i].Vx -= dtodx1*(phifc - phicl);
        	Wr[i].Vx -= dtodx1*(phicr - phifc);
      }
    }

#ifdef SELF_GRAVITY
    	for (i=il+1; i<=iu; i++) {
      		Wl[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      		Wr[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
    	}
#endif
		
		/*--- Step 1c (cont) -----------------------------------------------------------
		 * Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
		 * starting with tidal gravity terms added through the ShearingBoxPot
		 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
		 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
		 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
		 * (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
		 */		
		
#ifdef SHEARING_BOX
		if (ShearingBoxPot != NULL){
			for (i=il+1; i<=iu; i++) {
				cc_pos(pG,i,j,ks,&x1,&x2,&x3);
				phicr = (*ShearingBoxPot)( x1             ,x2,x3);
				phicl = (*ShearingBoxPot)((x1-    pG->dx1),x2,x3);
				phifc = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
				
				Wl[i].Vx -= dtodx1*(phifc - phicl);
				Wr[i].Vx -= dtodx1*(phicr - phifc);
			}
		}
		
		if (ShBoxCoord == xz){
			for (i=il+1; i<=iu; i++) {
				Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vz;
				Wr[i].Vx += pG->dt*Omega_0*W[i].Vz;
#ifdef FARGO
				Wl[i].Vz += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
				Wr[i].Vz += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
				Wl[i].Vz -= pG->dt*Omega_0*W[i-1].Vx;
				Wr[i].Vz -= pG->dt*Omega_0*W[i].Vx;
#endif
			}
		}
		
		if (ShBoxCoord == xy) {
			for (i=il+1; i<=iu; i++) {
				Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vy;
				Wr[i].Vx += pG->dt*Omega_0*W[i].Vy;
#ifdef FARGO
				Wl[i].Vy += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
				Wr[i].Vy += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
				Wl[i].Vy -= pG->dt*Omega_0*W[i-1].Vx;
				Wr[i].Vy -= pG->dt*Omega_0*W[i].Vx;
#endif
			}
		}
#endif /* SHEARING_BOX */	



	
/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 2D array
 */

	
	 for (i=il+1; i<=iu; i++) {
     		Ul_x1Face[j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);
      		Ur_x1Face[j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

#ifdef RADIATION_MHD
      		Bx = B1_x1Face[j][i];
#endif

		x1Flux[j][i].d = dt;
		x1Flux[j][i].Mx = 2; 
		/* This is used to take dt to calculate alpha. 
                * x1Flux[i].d is recovered in the fluxes  */

		fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);

	} /* End to get the fluxes */

	}/* End big loop j */


/*=== STEP 3: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 * NOTICE that momentum is rotated but flux and Eddington tensor are not *
 */

	for (i=il; i<=iu; i++) {
    		for (j=js-nghost; j<=je+nghost; j++) {
      			U1d[j].d  = pG->U[ks][j][i].d;
      			U1d[j].Mx = pG->U[ks][j][i].M2;
      			U1d[j].My = pG->U[ks][j][i].M3;
      			U1d[j].Mz = pG->U[ks][j][i].M1;
      			U1d[j].E  = pG->U[ks][j][i].E;
			U1d[j].Er  = pG->U[ks][j][i].Er;
    			U1d[j].Fr1  = pG->U[ks][j][i].Fr1;
    			U1d[j].Fr2  = pG->U[ks][j][i].Fr2;
    			U1d[j].Fr3  = pG->U[ks][j][i].Fr3;
			U1d[j].Edd_11  = pG->U[ks][j][i].Edd_11;
			U1d[j].Edd_21  = pG->U[ks][j][i].Edd_21;
			U1d[j].Edd_22  = pG->U[ks][j][i].Edd_22;
			for(m=0; m<NOPACITY;m++){
				U1d[j].Sigma[m] = pG->U[ks][j][i].Sigma[m];
			}


#ifdef RADIATION_MHD
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
/*
		Tguess = pG->Tguess[ks][j-1][i];
*/
		Tguess = temperature;
		velocity_x = U1d[j-1].Mz / U1d[j-1].d;
		velocity_y = U1d[j-1].Mx / U1d[j-1].d;
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i,j-1,ks,&x1,&x2,&x3);
		velocity_y -= qom * x1;		
#endif		


		velocity = velocity_x * velocity_x + velocity_y * velocity_y;
		
		Sigma_sF = U1d[j-1].Sigma[0];
		Sigma_aF = U1d[j-1].Sigma[1];
		Sigma_aP = U1d[j-1].Sigma[2];
		Sigma_aE = U1d[j-1].Sigma[3];


		Source[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[j-1].Fr1/U1d[j-1].d 
			- ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y)* U1d[j-1].Er / (Crat * U1d[j-1].d))	
			+ velocity_x * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[j-1].Er)/(Crat*U1d[j-1].d));
		Source[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[j-1].Fr2/U1d[j-1].d 
			- ((1.0 + U1d[j-1].Edd_22) * velocity_y + U1d[j-1].Edd_21 * velocity_x)* U1d[j-1].Er / (Crat * U1d[j-1].d))	
			+ velocity_y * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[j-1].Er)/(Crat*U1d[j-1].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * ((Sigma_aP * pow(Tguess, 4.0)
			- Sigma_aE * U1d[j-1].Er) + (Sigma_aF - Sigma_sF) * (velocity_x
			* (U1d[j-1].Fr1 - ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y) * U1d[j-1].Er / Crat)
			+ velocity_y
			* (U1d[j-1].Fr2 - ((1.0 + U1d[j-1].Edd_22) * velocity_y + U1d[j-1].Edd_21 * velocity_x) * U1d[j-1].Er / Crat))/Crat)
			- (Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2]) * U1d[j-1].d; 

		
	/*
		if(Opacity != NULL) Opacity(U1d[j-1].d, temperature, NULL, dSigma);
		
		dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
		dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
		dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
		dSigmadP[3] =  dSigma[7] / (density * R_ideal); 
	*/
		/*
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature /(U1d[j-1].d * R_ideal)
			-(Gamma - 1.0) * Prat * Crat * (dSigmadP[2] * pow(Tguess, 4.0) - dSigmadP[3] * U1d[j-1].Er)
		      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP[1] * (
			velocity_x * (U1d[j-1].Fr1 - ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y) * U1d[j-1].Er/Crat)
			+ velocity_y * (U1d[j-1].Fr2 - (U1d[j-1].Edd_21 * velocity_x + (1.0 + U1d[j-1].Edd_22) * velocity_y) * U1d[j-1].Er/Crat)
			);
		*/
		
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature * (1.0 - velocity /(Crat * Crat)) /(U1d[j-1].d * R_ideal);
		

		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */


		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
		SVVx = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[j-1].Edd_11) * W[j-1].Er / (W[j-1].d * Crat);
		
		if(fabs(SVVx * dt * 0.5) > 0.001)
		betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
		else 
		betax = 1.0 + 0.25 * SVVx * dt;

		SVVy = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[j-1].Edd_22) * W[j-1].Er / (W[j-1].d * Crat);
		
		if(fabs(SVVy * dt * 0.5) > 0.001)
		betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
		else 
		betay = 1.0 + 0.25 * SVVy * dt;
		/*===========================================================================*/


	
		/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
		Propa_44 = alpha;
        
        if(Erflag){

		/* "Vx" is actually vy, "vz" is actually vx; We stay with the correct meaning in source terms */
            Wl[j].Vx += dt * Source[2] * 0.5 * betay;
            Wl[j].Vz += dt * Source[1] * 0.5 * betax;
            Wl[j].P += dt * Propa_44 * Source[4] * 0.5;

            if(Wl[j].P < TINY_NUMBER)
                Wl[j].P -= dt * Propa_44 * Source[4] * 0.5;
            
        }

		for(m=0; m<NOPACITY; m++){
			Wl[j].Sigma[m] = U1d[j-1].Sigma[m];
		}

		
		Wl[j].Edd_11 = W[j-1].Edd_11;
		Wl[j].Edd_21 = W[j-1].Edd_21;
		Wl[j].Edd_22 = W[j-1].Edd_22;
		Wl[j].Edd_31 = W[j-1].Edd_31;
		Wl[j].Edd_32 = W[j-1].Edd_32;
		Wl[j].Edd_33 = W[j-1].Edd_33;

	/* For the right state */
	
	
		pressure = W[j].P;
		temperature = pressure / (U1d[j].d * R_ideal);
/*
		Tguess = pG->Tguess[ks][j][i];
*/
		Tguess = temperature;
		velocity_x = U1d[j].Mz / U1d[j].d;
		velocity_y = U1d[j].Mx / U1d[j].d;
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		velocity_y -= qom * x1;		
#endif

		velocity = velocity_x * velocity_x + velocity_y * velocity_y;
		
		
		Sigma_sF = U1d[j].Sigma[0];
		Sigma_aF = U1d[j].Sigma[1];
		Sigma_aP = U1d[j].Sigma[2];
		Sigma_aE = U1d[j].Sigma[3];


		Source[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[j].Fr1/U1d[j].d 
			- ((1.0 + U1d[j].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y)* U1d[j].Er / (Crat * U1d[j].d))	
			+ velocity_x * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[j].Er)/(Crat*U1d[j].d));
		Source[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (U1d[j].Fr2/U1d[j].d 
			- ((1.0 + U1d[j].Edd_22) * velocity_y + U1d[j].Edd_21 * velocity_x)* U1d[j].Er / (Crat * U1d[j].d))	
			+ velocity_y * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * U1d[j].Er)/(Crat*U1d[j].d));
		Source[4] = -(Gamma - 1.0) * Prat * Crat * ((Sigma_aP * pow(Tguess, 4.0) 
			- Sigma_aE * U1d[j].Er) + (Sigma_aF - Sigma_sF) * (velocity_x
			* (U1d[j].Fr1 - ((1.0 + U1d[j].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y) * U1d[j].Er / Crat)
			+ velocity_y
			* (U1d[j].Fr2 - ((1.0 + U1d[j].Edd_22) * velocity_y + U1d[j].Edd_21 * velocity_x) * U1d[j].Er / Crat))/Crat)
			- (Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2]) * U1d[j].d; 

	/*	
		if(Opacity != NULL) Opacity(U1d[j].d, temperature, NULL, dSigma);
		
		dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
		dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
		dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
		dSigmadP[3] =  dSigma[7] / (density * R_ideal); 
	*/
		/*
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature /(U1d[j].d * R_ideal)
			-(Gamma - 1.0) * Prat * Crat * (dSigmadP[2] * pow(Tguess,4.0) - dSigmadP[3] * U1d[j].Er)
		      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP[1] * (
			velocity_x * (U1d[j].Fr1 - ((1.0 + U1d[j].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y) * U1d[j].Er/Crat)
			+ velocity_y * (U1d[j].Fr2 - (U1d[j].Edd_21 * velocity_x + (1.0 + U1d[j].Edd_22) * velocity_y) * U1d[j].Er/Crat)
			);
		*/
		
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * temperature * temperature * temperature * (1.0 - velocity/(Crat * Crat))/(U1d[j].d * R_ideal);
		

		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */

		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
		SVVx = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[j].Edd_11) * W[j].Er / (W[j].d * Crat);
		
		if(fabs(SVVx * dt * 0.5) > 0.001)
		betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
		else 
		betax = 1.0 + 0.25 * SVVx * dt;

		SVVy = -Prat * (Sigma_aF + Sigma_sF) * (1.0 + W[j].Edd_22) * W[j].Er / (W[j].d * Crat);
		
		if(fabs(SVVy * dt * 0.5) > 0.001)
		betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
		else 
		betay = 1.0 + 0.25 * SVVy * dt;
		/*===========================================================================*/

	
		/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
		Propa_44 = alpha;
        
        if(Erflag){

            /* "vx" is actually vy, "vy" is actually vx */
            Wr[j].Vx += dt * Source[2] * 0.5 * betay;
            Wr[j].Vz += dt * Source[1] * 0.5 * betax;
            Wr[j].P += dt * Propa_44 * Source[4] * 0.5;

            if(Wr[j].P < TINY_NUMBER)
                Wr[j].P -= dt * Propa_44 * Source[4] * 0.5;
            
        }

		for(m=0; m<NOPACITY; m++){
			Wr[j].Sigma[m] = Sigma[m];
		}
	
		Wr[j].Edd_11 = W[j].Edd_11;
		Wr[j].Edd_21 = W[j].Edd_21;
		Wr[j].Edd_22 = W[j].Edd_22;
		Wr[j].Edd_31 = W[j].Edd_31;
		Wr[j].Edd_32 = W[j].Edd_32;
		Wr[j].Edd_33 = W[j].Edd_33;


	}


/********************************************
 * Add source term due to magnetic field */
/* notice that variables in y direction are rotated. "Bz" is actually Bx here */
#ifdef RADIATION_MHD
	for (j=jl+1; j<=ju; j++) {
      		MHD_src = (pG->U[ks][j-1][i].M1/pG->U[ks][j-1][i].d)*
        		(pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i])/dx2;
      		Wl[j].Bz += hdt*MHD_src;

      		MHD_src = (pG->U[ks][j][i].M1/pG->U[ks][j][i].d)*
        		(pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i])/dx2;
      		Wr[j].Bz += hdt*MHD_src;
    }

#endif

/*
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
		

	 if (StaticGravPot != NULL){
      		for (j=jl+1; j<=ju; j++) {
        		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        		phicr = (*StaticGravPot)(x1, x2             ,x3);
        		phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
        		phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        		Wl[j].Vx -= dtodx2*(phifc - phicl);
        		Wr[j].Vx -= dtodx2*(phicr - phifc);
      		}
    	}

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    	for (j=jl+1; j<=ju; j++) {
      		Wl[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      		Wr[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
    	}
#endif




	for (j=jl+1; j<=ju; j++) {
      		Ul_x2Face[j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);
      		Ur_x2Face[j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);
#ifdef RADIATION_MHD
      		Bx = B2_x2Face[j][i];
#endif
		x2Flux[j][i].d = dt;
		x2Flux[j][i].Mx = 2;
		/* Take the time step with x2Flux[][].d */	

      		fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[j],Wr[j],Bx,&x2Flux[j][i]);

		/* Note that x2Flux[][].Mx is actually momentum flux for vy ----
		-------------x2Flux[][].Mz is the actualy momentum flux for vx */
    	}/* End loop j to calculate the flux */

		}/*  End big loop i */



/******For magnetic field, calculate cell centered value of emf_3 at t and integrate to corner */
#ifdef RADIATION_MHD
	for (j=jl; j<=ju; j++) {
    	for (i=il; i<=iu; i++) {
      	emf3_cc[j][i] =
	(pG->U[ks][j][i].B1c*pG->U[ks][j][i].M2 -
	 pG->U[ks][j][i].B2c*pG->U[ks][j][i].M1 )/pG->U[ks][j][i].d;
    }
  }
  integrate_emf3_corner(pG);
/* update interface magnetic field */
	
  	for (j=jl+1; j<=ju-1; j++) {
    	for (i=il+1; i<=iu-1; i++) {

      	B1_x1Face[j][i] -= hdtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      	B2_x2Face[j][i] += hdtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }
    	B1_x1Face[j][iu] -= hdtodx2*(emf3[j+1][iu] - emf3[j][iu]);
  }
  	for (i=il+1; i<=iu-1; i++) {
    	B2_x2Face[ju][i] += hdtodx1*(emf3[ju][i+1] - emf3[ju][i]);
  }

#endif /* end MHD part */




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
#ifdef RADIATION_MHD
      		Ul_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i-1].By - x2Flux[j][i-1].By);
#endif

      		Ur_x1Face[j][i].d  -= hdtodx2*(x2Flux[j+1][i  ].d  - x2Flux[j][i  ].d );
      		Ur_x1Face[j][i].Mx -= hdtodx2*(x2Flux[j+1][i  ].Mz - x2Flux[j][i  ].Mz);
      		Ur_x1Face[j][i].My -= hdtodx2*(x2Flux[j+1][i  ].Mx - x2Flux[j][i  ].Mx);
      		Ur_x1Face[j][i].Mz -= hdtodx2*(x2Flux[j+1][i  ].My - x2Flux[j][i  ].My);
		Ur_x1Face[j][i].E  -= hdtodx2*(x2Flux[j+1][i  ].E -  x2Flux[j][i  ].E );
#ifdef RADIATION_MHD
      		Ur_x1Face[j][i].Bz -= hdtodx2*(x2Flux[j+1][i  ].By - x2Flux[j][i  ].By);
#endif

    		}
  	}


/*-------step 6b------------------------------*/
/* Correct x1 interface with x2 source gradient */
/* Add MHD source terms from x2 flux gradient to the conservative variable on x1Face */
#ifdef RADIATION_MHD
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
/* The left state */
      dbx = pG->B1i[ks][j][i] - pG->B1i[ks][j][i-1];

      B1 = pG->U[ks][j][i-1].B1c;
      B2 = pG->U[ks][j][i-1].B2c;
      B3 = pG->U[ks][j][i-1].B3c;
      V3 = pG->U[ks][j][i-1].M3/pG->U[ks][j][i-1].d;

      Ul_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ul_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ul_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ul_x1Face[j][i].Bz += hdtodx1*V3*dbx;

      Ul_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;

/* the right state */
      dbx = pG->B1i[ks][j][i+1] - pG->B1i[ks][j][i];

      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x1Face[j][i].Mx += hdtodx1*B1*dbx;
      Ur_x1Face[j][i].My += hdtodx1*B2*dbx;
      Ur_x1Face[j][i].Mz += hdtodx1*B3*dbx;
      Ur_x1Face[j][i].Bz += hdtodx1*V3*dbx;

      Ur_x1Face[j][i].E  += hdtodx1*B3*V3*dbx;

    }
  }

#endif		


/*----Need extra flux gradient for self-gravity and gravitational source that is included as flux 
 * but not included in the Riemann problem---------*/

/*
 * Add source terms for a static gravitational potential arising from x2-Flux
 * gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)(x1, x2             ,x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;

        Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                      x2Flux[j+1][i  ].d*(phir - phic));

        phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);


        Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;

        Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                      x2Flux[j+1][i-1].d*(phir - phic));

      }
    }
  }


#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu; i++) {
      	phic = pG->Phi[ks][j][i];
      	phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      	phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);

      	Ur_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;

      	Ur_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i  ].d*(phic - phil) +
                                    x2Flux[j+1][i  ].d*(phir - phic));

        phic = pG->Phi[ks][j][i-1];
        phir = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j+1][i-1]);
        phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j-1][i-1]);

        Ul_x1Face[j][i].My -= hdtodx2*(phir-phil)*pG->U[ks][j][i-1].d;

/* Adding source term still uses non-cnoservative form */
        Ul_x1Face[j][i].E -= hdtodx2*(x2Flux[j  ][i-1].d*(phic - phil) +
                                    x2Flux[j+1][i-1].d*(phir - phic));

    }
  }
#endif /* SELF_GRAVITY */




/*=== STEP 7: Correct x2-interface states with transverse flux gradients =====*/


	
  	for (j=jl+1; j<=ju; j++) {
    		for (i=il+1; i<=iu-1; i++) {

      		Ul_x2Face[j][i].d  -= hdtodx1*(x1Flux[j-1][i+1].d  - x1Flux[j-1][i].d );
      		Ul_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j-1][i+1].My - x1Flux[j-1][i].My);
      		Ul_x2Face[j][i].My -= hdtodx1*(x1Flux[j-1][i+1].Mz - x1Flux[j-1][i].Mz);
      		Ul_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j-1][i+1].Mx - x1Flux[j-1][i].Mx);
      		Ul_x2Face[j][i].E  -= hdtodx1*(x1Flux[j-1][i+1].E  - x1Flux[j-1][i].E );

#ifdef RADIATION_MHD
      		Ul_x2Face[j][i].By -= hdtodx1*(x1Flux[j-1][i+1].Bz - x1Flux[j-1][i].Bz);
#endif

      		Ur_x2Face[j][i].d  -= hdtodx1*(x1Flux[j  ][i+1].d  - x1Flux[j  ][i].d );
      		Ur_x2Face[j][i].Mx -= hdtodx1*(x1Flux[j  ][i+1].My - x1Flux[j  ][i].My);
      		Ur_x2Face[j][i].My -= hdtodx1*(x1Flux[j  ][i+1].Mz - x1Flux[j  ][i].Mz);
      		Ur_x2Face[j][i].Mz -= hdtodx1*(x1Flux[j  ][i+1].Mx - x1Flux[j  ][i].Mx);
      		Ur_x2Face[j][i].E  -= hdtodx1*(x1Flux[j  ][i+1].E  - x1Flux[j  ][i].E );

#ifdef RADIATION_MHD
      		Ur_x2Face[j][i].By -= hdtodx1*(x1Flux[j  ][i+1].Bz - x1Flux[j  ][i].Bz);
#endif

    }
  }


/*---------step 7b---------------*/
/* Correct x2 interface with x1 source gradient */
/* Add magnetic field source  term from x1-flux gradients */
#ifdef RADIATION_MHD

  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
/* for the left state */
      dby = pG->B2i[ks][j][i] - pG->B2i[ks][j-1][i];
      B1 = pG->U[ks][j-1][i].B1c;
      B2 = pG->U[ks][j-1][i].B2c;
      B3 = pG->U[ks][j-1][i].B3c;
      V3 = pG->U[ks][j-1][i].M3/pG->U[ks][j-1][i].d;

      Ul_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ul_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ul_x2Face[j][i].My += hdtodx2*B3*dby;
      Ul_x2Face[j][i].By += hdtodx2*V3*dby;

      Ul_x2Face[j][i].E  += hdtodx2*B3*V3*dby;
/* for the right state */

      dby = pG->B2i[ks][j+1][i] - pG->B2i[ks][j][i];
      B1 = pG->U[ks][j][i].B1c;
      B2 = pG->U[ks][j][i].B2c;
      B3 = pG->U[ks][j][i].B3c;
      V3 = pG->U[ks][j][i].M3/pG->U[ks][j][i].d;

      Ur_x2Face[j][i].Mz += hdtodx2*B1*dby;
      Ur_x2Face[j][i].Mx += hdtodx2*B2*dby;
      Ur_x2Face[j][i].My += hdtodx2*B3*dby;
      Ur_x2Face[j][i].By += hdtodx2*V3*dby;

      Ur_x2Face[j][i].E  += hdtodx2*B3*V3*dby;

    }
  }

#endif

/*-
 * Add source terms for a static gravitational potential arising from x1-Flux
 * gradients. To improve conservation of total energy, average the energy
 * source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);


        Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

        Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                      x1Flux[j  ][i+1].d*(phir - phic));


        phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);


        Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;

        Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                      x1Flux[j-1][i+1].d*(phir - phic));

      }
    }
  }


/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x1-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (j=jl+1; j<=ju; j++) {
    for (i=il+1; i<=iu-1; i++) {
      phic = pG->Phi[ks][j][i];
      phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);

      Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;

      Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
                                    x1Flux[j  ][i+1].d*(phir - phic));


      phic = pG->Phi[ks][j-1][i];
      phir = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i+1]);
      phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j-1][i-1]);

      Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;

/* Transverse flux correction still uses source term like format */
      Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
                                    x1Flux[j-1][i+1].d*(phir - phic));

    }
  }

#endif /* SELF_GRAVITY */	
	
	/*--- Step 6d (cont) -----------------------------------------------------------
	 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
	 * x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
	 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
	 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
	 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
	 * (x1,x2,x3) in code = (X,Z,Y) in shearing sheet
	 */
	
	/* No shearing source term for y component */
	
	
#ifdef SHEARING_BOX
	if (ShearingBoxPot != NULL){
		for (j=jl+1; j<=ju; j++) {
			for (i=il+1; i<=iu-1; i++) {
				cc_pos(pG,i,j,ks,&x1,&x2,&x3);
				phic = (*ShearingBoxPot)((x1            ),x2,x3);
				phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
				phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
				
				Ur_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
#ifndef BAROTROPIC
				Ur_x2Face[j][i].E -= hdtodx1*(x1Flux[j  ][i  ].d*(phic - phil) +
											  x1Flux[j  ][i+1].d*(phir - phic));
#endif
				
				phic = (*ShearingBoxPot)((x1            ),(x2-pG->dx2),x3);
				phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
				phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);
				
				Ul_x2Face[j][i].Mz -= hdtodx1*(phir-phil)*pG->U[ks][j-1][i].d;
#ifndef BAROTROPIC
				Ul_x2Face[j][i].E -= hdtodx1*(x1Flux[j-1][i  ].d*(phic - phil) +
											  x1Flux[j-1][i+1].d*(phir - phic));
#endif
			}
		}
	}
	
	if (ShBoxCoord == xz){
		for (j=jl+1; j<=ju; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j][i].M3;
				Ul_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j-1][i].M3;
#ifdef FARGO
				Ur_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
				Ul_x2Face[j][i].My += hdt*(qshear-2.)*Omega_0*pG->U[ks][j-1][i].M1;
#else
				Ur_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
				Ul_x2Face[j][i].My -= pG->dt*Omega_0*pG->U[ks][j-1][i].M1;
#endif
			}
		}
	}
	
	if (ShBoxCoord == xy){
		for (j=jl+1; j<=ju; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j][i].M2;
				Ul_x2Face[j][i].Mz += pG->dt*Omega_0*pG->U[ks][j-1][i].M2;
#ifdef FARGO
				Ur_x2Face[j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
				Ul_x2Face[j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[ks][j-1][i].M1;
#else
				Ur_x2Face[j][i].Mx -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
				Ul_x2Face[j][i].Mx -= pG->dt*Omega_0*pG->U[ks][j-1][i].M1;
#endif
			}
		}
	}
#endif /* SHEARING_BOX */
	




/************************************************************/
/* calculate d^{n+1/2}, this is needed to calculate emf3_cc^{n+1/2} */

 for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {

        dhalf[j][i] = pG->U[ks][j][i].d
          - hdtodx1*(x1Flux[j  ][i+1].d - x1Flux[j][i].d)
          - hdtodx2*(    x2Flux[j+1][i  ].d -     x2Flux[j][i].d);

      }
    }

#ifdef RADIATION_MHD
/* Now calculate cell centered emf3_cc^{n+1/2} */
      for (j=jl+1; j<=ju-1; j++) {
    for (i=il+1; i<=iu-1; i++) {


	/**************************************************/
	/* Update momentum using modified Godunov method **/
	divFlux1[1] = (x1Flux[j][i+1].Mx - x1Flux[j][i].Mx) / dx1;
	divFlux1[2] = (x1Flux[j][i+1].My - x1Flux[j][i].My) / dx1;

	divFlux2[1] = (x2Flux[j+1][i].Mz - x2Flux[j][i].Mz) / dx2;
	divFlux2[2] = (x2Flux[j+1][i].Mx - x2Flux[j][i].Mx) / dx2;

/* Add source terms for fixed gravitational potential */
        if (StaticGravPot != NULL){
        	cc_pos(pG,i,j,ks,&x1,&x2,&x3);

        	phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        	phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
		divFlux1[1] += (phir-phil)*pG->U[ks][j][i].d/dx1;


        	phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        	phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

		divFlux2[2] += (phir-phil)*pG->U[ks][j][i].d/dx2; 	
      	}

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      	phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i+1]);
      	phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j][i-1]);
	divFlux1[1] += (phir-phil)*pG->U[ks][j][i].d / dx1;

      	phir = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j+1][i]);
      	phil = 0.5*(pG->Phi[ks][j][i] + pG->Phi[ks][j-1][i]);
      	divFlux2[2] += (phir-phil)*pG->U[ks][j][i].d / dx2;

#endif /* SELF_GRAVITY */
		
		
/*****************************************************************
 * Add source  term from radiation *******************************/
/* We only need momentum source terms, not energy source  term */
	Usource.d  = pG->U[ks][j][i].d;
	Usource.Mx = pG->U[ks][j][i].M1;
    Usource.My = pG->U[ks][j][i].M2;
    Usource.Mz = pG->U[ks][j][i].M3;
    Usource.E  = pG->U[ks][j][i].E;
	Usource.Er  = pG->U[ks][j][i].Er;
    Usource.Fr1  = pG->U[ks][j][i].Fr1;
    Usource.Fr2  = pG->U[ks][j][i].Fr2;
    Usource.Fr3  = pG->U[ks][j][i].Fr3;
	Usource.Edd_11  = pG->U[ks][j][i].Edd_11;
	Usource.Edd_21  = pG->U[ks][j][i].Edd_21;
	Usource.Edd_22  = pG->U[ks][j][i].Edd_22;
	Usource.Edd_31	= pG->U[ks][j][i].Edd_31;
	Usource.Edd_32	= pG->U[ks][j][i].Edd_32;
	Usource.Edd_33	= pG->U[ks][j][i].Edd_33;
	for(m=0; m<NOPACITY; m++){
		Usource.Sigma[m] = pG->U[ks][j][i].Sigma[m];
	}

	Usource.By = pG->U[ks][j][i].B2c;
      	Usource.Bz = pG->U[ks][j][i].B3c;
	Bx = pG->U[ks][j][i].B1c;


	density = Usource.d;
	velocity_x = Usource.Mx / density;
	velocity_y = Usource.My / density;
	velocity_z = Usource.Mz / density;
	

	velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
		
	pressure = (Usource.E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */

	pressure -= 0.5 * (Bx * Bx + Usource.By * Usource.By + Usource.Bz * Usource.Bz) * (Gamma - 1.0);

	temperature = pressure / (density * R_ideal);
	
	diffTEr = Usource.Sigma[2] * pow(temperature, 4.0) - Usource.Sigma[3] * Usource.Er;
		
	/*************************************************/
	/* Note that we need to add background shearing after temperature is calculated */
		
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		/* Include background shearing in Usource, which is used in dSource */
		velocity_y -= qom * x1;	

#endif
			

	/* The Source term */
	dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, NULL, x1);

	/*=========================================================*/
	/* In case velocity is large and momentum source is stiff */
	SFmx = (Usource.Sigma[0] + Usource.Sigma[1]) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat); 
	/*	+ diffTEr / (density * Crat); */	

	SFmy = (Usource.Sigma[0] + Usource.Sigma[1]) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat); 
	/*	+ diffTEr / (density * Crat);*/

	
	Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
	Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
	
	/*=========================================================*/

	/* co-moving flux */
	Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y) * Usource.Er / Crat;
	Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x) * Usource.Er / Crat;
	

	/* Source term for momentum, not velocity*/
	Source[1] = -Prat * (-(Usource.Sigma[0] + Usource.Sigma[1]) * Fr0x + velocity_x * diffTEr / Crat);
	Source[2] = -Prat * (-(Usource.Sigma[0] + Usource.Sigma[1]) * Fr0y + velocity_y * diffTEr / Crat);


	/* Now we have source term and flux, update the momentum */
	/* Update momentum for half time step */
	M1h = pG->U[ks][j][i].M1 + hdt * Source_Inv[1][1] * (Source[1] - divFlux1[1] - divFlux2[1]);
	M2h = pG->U[ks][j][i].M2 + hdt * Source_Inv[2][2] * (Source[2] - divFlux1[2] - divFlux2[2]);


	/* We only need first order accuracy for momentum.*/
	/* Do not need to do correct step here */ 
/*
	velocity_x = Uguess[1] / density;
	velocity_y = Uguess[2] / density;
		
#ifdef FARGO

		velocity_y -= qom * x1;		
#endif
	
	Source_guess[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
	Source_guess[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);



	Errort[1] = pG->U[ks][j][i].M1 + 0.5 * hdt * (Source[1] + Source_guess[1]) - hdt * (divFlux1[1] + divFlux2[1]) - Uguess[1];
	Errort[2] = pG->U[ks][j][i].M2 + 0.5 * hdt * (Source[2] + Source_guess[2]) - hdt * (divFlux1[2] + divFlux2[2]) - Uguess[2];

	M1h = Uguess[1] + Source_Inv[1][1] * Errort[1];
	M2h = Uguess[2] + Source_Inv[2][2] * Errort[2];
*/	
	


  /****************************************************************/
		
/****************************************************************/
/* Add shearing box source terms, which are seperated from 
 * radiation field completely 
 */
		
		
		/* Add the tidal gravity and Coriolis terms for shearing box. */
#ifdef SHEARING_BOX
		if (ShearingBoxPot != NULL){
			cc_pos(pG,i,j,ks,&x1,&x2,&x3);
			phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
			phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
			M1h -= hdtodx1*(phir-phil)*pG->U[ks][j][i].d;
			
			phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
			phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
			M2h -= hdtodx2*(phir-phil)*pG->U[ks][j][i].d;
		}
		
		if (ShBoxCoord == xy) M1h += pG->dt*Omega_0*pG->U[ks][j][i].M2;
		if (ShBoxCoord == xz) M1h += pG->dt*Omega_0*pG->U[ks][j][i].M3;
#ifdef FARGO
		if (ShBoxCoord == xy) M2h += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
		if (ShBoxCoord == xz) M3h += hdt*(qshear-2.)*Omega_0*pG->U[ks][j][i].M1;
#else
		if (ShBoxCoord == xy) M2h -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
		if (ShBoxCoord == xz) M3h -= pG->dt*Omega_0*pG->U[ks][j][i].M1;
#endif
#endif /* SHEARING_BOX */		
		
		


/*
      Eh = pG->U[ks][j][i].E
        - hdtodx1*(x1Flux[j][i+1].E - x1Flux[j][i].E)
        - hdtodx2*(    x2Flux[j+1][i].E -     x2Flux[j][i].E);
*/

/* phalf is the gas pressure at half step, we do not calculate here. */
/* phalf is needed for cooling only */
/*
      phalf[j][i] = Eh - 0.5*(M1h*M1h + M2h*M2h + M3h*M3h)/dhalf[j][i];
*/
      B1ch = 0.5*(B1_x1Face[j][i] + B1_x1Face[j][i+1]);
      B2ch = 0.5*(    B2_x2Face[j][i] +     B2_x2Face[j+1][i]);
      B3ch = pG->U[ks][j][i].B3c 
        - hdtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz)
        - hdtodx2*(    x2Flux[j+1][i].By -     x2Flux[j][i].By);
      emf3_cc[j][i] = (B1ch*M2h - B2ch*M1h)/dhalf[j][i];

/*
      phalf[j][i] -= 0.5*(B1ch*B1ch + B2ch*B2ch + B3ch*B3ch);
      phalf[j][i] *= Gamma_1;
*/
	}
     }

#endif

/*=== STEP 8: Compute 2D x1-Flux, x2-Flux ====================================*/


/*--- Step 8a ------------------------------------------------------------------
* Compute 2D x1-fluxes from corrected L/R states.
 */

	 for (j=js-1; j<=je+1; j++) {
    		for (i=is; i<=ie+1; i++) {

#ifdef RADIATION_MHD
      		Bx = B1_x1Face[j][i];
#endif
      		Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[j][i],&Bx);
      		Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[j][i],&Bx);

		/* Need parameter dt in radiation Riemann solver */
		x1Flux[j][i].d = dt;
		x1Flux[j][i].Mx = 2;

      		fluxes(Ul_x1Face[j][i],Ur_x1Face[j][i],Wl[i],Wr[i],Bx,&x1Flux[j][i]);
    		}
  	}

	
/*--- Step 8b ------------------------------------------------------------------
 * Compute 2D x2-fluxes from corrected L/R states.
 */

  	for (j=js; j<=je+1; j++) {
    		for (i=is-1; i<=ie+1; i++) {

#ifdef RADIATION_MHD
      		Bx = B2_x2Face[j][i];
#endif
      		Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[j][i],&Bx);
      		Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[j][i],&Bx);

		/* take dt to calculate effective sound speed */
		x2Flux[j][i].d = dt;
		x2Flux[j][i].Mx = 2;

      		fluxes(Ul_x2Face[j][i],Ur_x2Face[j][i],Wl[i],Wr[i],Bx,&x2Flux[j][i]);
    }
  }





/******************************************/
/* integrate emf*^{n+1/2} to the grid cell corners */
#ifdef RADIATION_MHD
  integrate_emf3_corner(pG);

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->B1i[ks][j][i] -= dtodx2*(emf3[j+1][i  ] - emf3[j][i]);
      pG->B2i[ks][j][i] += dtodx1*(emf3[j  ][i+1] - emf3[j][i]);
    }

    pG->B1i[ks][j][ie+1] -= dtodx2*(emf3[j+1][ie+1] - emf3[j][ie+1]);
  }
  for (i=is; i<=ie; i++) {
    pG->B2i[ks][je+1][i] += dtodx1*(emf3[je+1][i+1] - emf3[je+1][i]);
  }

#endif /* end mhd part */




/*-------Step 9: Predict and correct step after we get the flux------------- */
	

#ifdef SELF_GRAVITY
/* Save mass fluxes in Grid structure for source term correction in main loop */
/* After we get the corrected fluxt */

  for (j=js; j<=je+1; j++) {
    for (i=is; i<=ie+1; i++) {
      pG->x1MassFlux[ks][j][i] = x1Flux[j][i].d;
      pG->x2MassFlux[ks][j][i] = x2Flux[j][i].d;
    }
  }

#ifdef CONS_GRAVITY

/* Need to update the density first so that we can calculate the new potential */
/* There is no source term for density from radiation */
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
		/* We need the density at time step n, for radiation source term and  */
		/* and self-gravity energy flux term */
		density_old[j][i] = pG->U[ks][j][i].d;
		pG->U[ks][j][i].d  -= dtodx1*(x1Flux[j][i+1].d -x1Flux[j][i].d );
		pG->U[ks][j][i].d  -= dtodx2*(x2Flux[j+1][i].d -x2Flux[j][i].d );

      }
    }
  
	bvals_mhd(pD);

/* With the updated flux, now we can calculate dphidt with Poisson Solver */
    	(*SelfGrav_cons)(pD);
	/* Need to apply boundary condition for the new calculated dphidt */
    	bvals_grav(pD);
/* Now calculate the energy flux */

  	for(j=js; j<=je+1;j++){	
		for (i=is; i<=ie+1; i++) {
			phil = 0.25*(pG->Phi[ks][j][i-1]+pG->Phi_old[ks][j][i]+pG->Phi_old[ks][j][i-1]+pG->Phi[ks][j][i]);
			gxl = 0.5 * (pG->Phi[ks][j][i-1] + pG->Phi_old[ks][j][i-1]  - pG->Phi[ks][j][i  ] - pG->Phi_old[ks][j][i ])/(pG->dx1);
			dotphil  = 0.5*(pG->dphidt[ks][j][i-1] + pG->dphidt[ks][j][i  ]);		
			dotgxl = (pG->dphidt[ks][j][i-1] - pG->dphidt[ks][j][i  ])/(pG->dx1);

			x1Flux_grav[j][i] =-0.5*(phil*dotgxl-dotphil*gxl)/four_pi_G + x1Flux[j][i].d*phil;

			phil = 0.25*(pG->Phi[ks][j-1][i]+pG->Phi_old[ks][j-1][i]+pG->Phi_old[ks][j][i]+pG->Phi[ks][j][i]);
			gxl = 0.5 * (pG->Phi[ks][j-1][i] + pG->Phi_old[ks][j-1][i]  - pG->Phi[ks][j][i  ] - pG->Phi_old[ks][j][i  ])/(pG->dx2);
			dotphil  = 0.5*(pG->dphidt[ks][j-1][i] + pG->dphidt[ks][j][i  ]);		
			dotgxl = (pG->dphidt[ks][j-1][i] - pG->dphidt[ks][j][i  ])/(pG->dx2);

			x2Flux_grav[j][i] =-0.5*(phil*dotgxl-dotphil*gxl)/four_pi_G + x2Flux[j][i].d*phil;

    		
		}
 	}
	
	/*-----------Now swap, desntiy_old actually save the new density  */


    	for (j=js; j<=je; j++) {
      		for (i=is; i<=ie; i++) {
		/* We need the density at time step n, for radiation source term and  */
		/* and self-gravity energy flux term */
		/* We do not need to update ghost zone here */
			Tempswap = pG->U[ks][j][i].d;
			pG->U[ks][j][i].d = density_old[j][i];
			density_old[j][i] = Tempswap;
      		}
    	}
  


#endif

#endif /* SELF_GRAVITY */	
	
	
/*----------Radiation quantities are not updated in the modified Godunov corrector step */
/*-----------NOTE that x1flux and x2flux are flux from Remann problem. If there is extra-----* 
 *----------source terms, flux due to those source terms should also be added in the-------*
 *---------- modified Godunov method -------------------------------------*/
	for (j=js; j<=je; j++) {
    		for (i=is; i<=ie; i++) {		
				
				
				
		/* Load 1D vector */
			Usource.d  = pG->U[ks][j][i].d;
      			Usource.Mx = pG->U[ks][j][i].M1;
      			Usource.My = pG->U[ks][j][i].M2;
      			Usource.Mz = pG->U[ks][j][i].M3;
      			Usource.E  = pG->U[ks][j][i].E;
			Usource.Er  = pG->U[ks][j][i].Er;
    			Usource.Fr1  = pG->U[ks][j][i].Fr1;
    			Usource.Fr2  = pG->U[ks][j][i].Fr2;
    			Usource.Fr3  = pG->U[ks][j][i].Fr3;
			Usource.Edd_11  = pG->U[ks][j][i].Edd_11;
			Usource.Edd_21  = pG->U[ks][j][i].Edd_21;
			Usource.Edd_22  = pG->U[ks][j][i].Edd_22;
			Usource.Edd_31	= pG->U[ks][j][i].Edd_31;
			Usource.Edd_32	= pG->U[ks][j][i].Edd_32;
			Usource.Edd_33	= pG->U[ks][j][i].Edd_33;
			for(m=0;m<NOPACITY;m++){
				Usource.Sigma[m] = pG->U[ks][j][i].Sigma[m];
			}		

	
#ifdef RADIATION_MHD
      			Usource.By = pG->U[ks][j][i].B2c;
      			Usource.Bz = pG->U[ks][j][i].B3c;
			Bx = pG->U[ks][j][i].B1c;
#else
			Bx = 0.0;
#endif /* MHD */

	
		density = pG->U[ks][j][i].d;
				
		pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 
				+ pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2 + pG->U[ks][j][i].M3 * pG->U[ks][j][i].M3) / density)	* (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c) * (Gamma - 1.0);
#endif
		temperature = pressure / (density * R_ideal);
		/* Tguess is uesed for source term T^4 - Er */
/*
		Tguess = pG->Tguess[ks][j][i];
*/
		Tguess = temperature;

		velocity_x = pG->U[ks][j][i].M1 / density;
		velocity_y = pG->U[ks][j][i].M2 / density;
#ifdef FARGO
	/* With FARGO, we should add background shearing to the source terms */

		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		/* Include background shearing in Usource, which is only used in dSource */
				
		velocity_y -= qom * x1;				
	
#endif				

				
		Sigma_sF = pG->U[ks][j][i].Sigma[0];
		Sigma_aF = pG->U[ks][j][i].Sigma[1];
		Sigma_aP = pG->U[ks][j][i].Sigma[2];
		Sigma_aE = pG->U[ks][j][i].Sigma[3];

		Sigma[0] = Sigma_sF;
		Sigma[1] = Sigma_aF;
		Sigma[2] = Sigma_aP;
		Sigma[3] = Sigma_aE;

		/* Part of momentum source term */
		Smx0 = -Prat * velocity_x * (Sigma_aP * pG->Tguess[ks][j][i] - Sigma_aE * Usource.Er) / Crat;
		Smy0 = -Prat * velocity_y * (Sigma_aP * pG->Tguess[ks][j][i] - Sigma_aE * Usource.Er) / Crat;

	/* The Source term */
		dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, NULL, x1);

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
		SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat); 
		/*	+ (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * Usource.Er) / (density * Crat);*/	

		SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat); 
		/*	+ (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * Usource.Er) / (density * Crat);*/	

		Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
		Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);

		/*=========================================================*/
		Source_Inv[0][0] = 1.0;
		Source_Inv[4][0] = 0.0;
		Source_Inv[4][1] = 0.0;
		Source_Inv[4][2] = 0.0;
		Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	
		Source[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat));
			/*+ velocity_x * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * pG->U[ks][j][i].Er)/Crat);*/
		Source[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat));
			/*+ velocity_y * (Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * pG->U[ks][j][i].Er)/Crat);*/
		if(Erflag){
			Source[4] = -Prat * Crat * ((Sigma_aP * pow(Tguess, 4.0) - Sigma_aE * pG->U[ks][j][i].Er) 
			+ Source_Inv[1][1] * (Sigma_aF - Sigma_sF) * velocity_x
			* (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)/Crat
			+ Source_Inv[2][2] * (Sigma_aF - Sigma_sF) * velocity_y
			* (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)/Crat); 
		}
		else{
			Source[4] = -Prat * Crat * (Source_Inv[1][1] * (Sigma_aF - Sigma_sF) * velocity_x
			* (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)/Crat
			+ Source_Inv[2][2] * (Sigma_aF - Sigma_sF) * velocity_y
			* (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat)/Crat); 
		}
		Prwork1 = -Prat * (Sigma_aF - Sigma_sF) * (Source_Inv[1][1] * velocity_x
			* (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat)		
			+ Source_Inv[2][2] * velocity_y
			* (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat));

		

	/*--------Calculate the guess solution-------------*/
	
		divFlux1[0] = (x1Flux[j][i+1].d  - x1Flux[j][i].d ) / dx1;
		divFlux1[1] = (x1Flux[j][i+1].Mx - x1Flux[j][i].Mx) / dx1;
		divFlux1[2] = (x1Flux[j][i+1].My - x1Flux[j][i].My) / dx1;
		divFlux1[4] = (x1Flux[j][i+1].E  - x1Flux[j][i].E ) / dx1; 

		divFlux2[0] = (x2Flux[j+1][i].d  - x2Flux[j][i].d ) / dx2;
		divFlux2[1] = (x2Flux[j+1][i].Mz - x2Flux[j][i].Mz) / dx2;
		divFlux2[2] = (x2Flux[j+1][i].Mx - x2Flux[j][i].Mx) / dx2;
		divFlux2[4] = (x2Flux[j+1][i].E  - x2Flux[j][i].E ) / dx2; 

	/*------Flux due to static gravitational potential ---------*/

		if (StaticGravPot != NULL){
    			cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        		phic = (*StaticGravPot)((x1            ),x2,x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

			SourceFlux[1] = dhalf[j][i]*(phir-phil)/dx1;
        		SourceFlux[4]= (x1Flux[j][i  ].d*(phic - phil) +
				 	x1Flux[j][i+1].d*(phir - phic))/dx1;

			divFlux1[1] += SourceFlux[1];
			divFlux1[4] += SourceFlux[4];
			
        		phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        		phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        		SourceFlux[2] = dhalf[j][i]*(phir-phil)/dx2;

        		SourceFlux[4] = (x2Flux[j  ][i].d*(phic - phil) +
                                	     x2Flux[j+1][i].d*(phir - phic))/dx2;
			divFlux2[2] += SourceFlux[2];
			divFlux2[4] += SourceFlux[4];      			
  		}

	/*----- Add flux due to self-gravity---------*/
#ifdef SELF_GRAVITY
		/*-----------x direction ------------*/
      		phic = pG->Phi[ks][j][i];
      		phil = 0.5*(pG->Phi[ks][j][i-1] + pG->Phi[ks][j][i  ]);
      		phir = 0.5*(pG->Phi[ks][j][i  ] + pG->Phi[ks][j][i+1]);

/* gx and gy centered at L and R x1-faces */
      		gxl = (pG->Phi[ks][j][i-1] - pG->Phi[ks][j][i  ])/(pG->dx1);
      		gxr = (pG->Phi[ks][j][i  ] - pG->Phi[ks][j][i+1])/(pG->dx1);

      		gyl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j+1][i-1]) +
                	  (pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) )/(pG->dx2);
      		gyr = 0.25*((pG->Phi[ks][j-1][i  ] - pG->Phi[ks][j+1][i  ]) +
                  	(pG->Phi[ks][j-1][i+1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx2);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
      		flux_m1l = 0.5*(gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
      		flux_m1r = 0.5*(gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

      		flux_m2l = gxl*gyl/four_pi_G;
      		flux_m2r = gxr*gyr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
      		divFlux1[1] += (flux_m1r - flux_m1l) / dx1;
      		divFlux1[2] += (flux_m2r - flux_m2l) / dx1;

#ifndef CONS_GRAVITY
      		divFlux1[4] += (x1Flux[j][i  ].d*(phic - phil) +
                                   x1Flux[j][i+1].d*(phir - phic)) / dx1;
#endif

		/*---------y direction ----------------*/
      		phil = 0.5*(pG->Phi[ks][j-1][i] + pG->Phi[ks][j  ][i]);
      		phir = 0.5*(pG->Phi[ks][j  ][i] + pG->Phi[ks][j+1][i]);

/* gx and gy centered at L and R x2-faces */
      		gxl = 0.25*((pG->Phi[ks][j-1][i-1] - pG->Phi[ks][j-1][i+1]) +
                  	(pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) )/(pG->dx1);
      		gxr = 0.25*((pG->Phi[ks][j  ][i-1] - pG->Phi[ks][j  ][i+1]) +
                  	(pG->Phi[ks][j+1][i-1] - pG->Phi[ks][j+1][i+1]) )/(pG->dx1);

      		gyl = (pG->Phi[ks][j-1][i] - pG->Phi[ks][j  ][i])/(pG->dx2);
      		gyr = (pG->Phi[ks][j  ][i] - pG->Phi[ks][j+1][i])/(pG->dx2);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
      		flux_m1l = gyl*gxl/four_pi_G;
      		flux_m1r = gyr*gxr/four_pi_G;

      		flux_m2l = 0.5*(gyl*gyl-gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      		flux_m2r = 0.5*(gyr*gyr-gxr*gxr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx1 terms  */
      		divFlux2[1] += (flux_m1r - flux_m1l) / dx2;
      		divFlux2[2] += (flux_m2r - flux_m2l) / dx2;

#ifndef CONS_GRAVITY
      		divFlux1[4] += (x2Flux[j  ][i].d*(phic - phil) +
                                   x2Flux[j+1][i].d*(phir - phic)) / dx2;
#endif			

	/* Now the energy flux due to self-gravity in cons-gravity */

#ifdef CONS_GRAVITY
		divFlux1[4] += (x1Flux_grav[j][i+1] - x1Flux_grav[j][i]) / dx1;
		divFlux2[4] += (x2Flux_grav[j+1][i] - x2Flux_grav[j][i]) / dx2;
		

#endif


#endif


	/*================== End static gravitational flux ================*/
		
		Source[0] = 0.0;
		for(n=0; n<5; n++) {
			tempguess[n] = dt * Source_Inv[n][n] * (Source[n] - divFlux1[n] - divFlux2[n]);
			
		}

		if(!Erflag){
			tempguess[4] += -Prat * pG->Ersource[ks][j][i];

		}

		Prworksource = dt * Source_Inv[4][4] * Prwork1;

		Uguess[0] = pG->U[ks][j][i].d + tempguess[0];
		Uguess[1] = pG->U[ks][j][i].M1 + tempguess[1] + dt * Smx0;
		Uguess[2] = pG->U[ks][j][i].M2 + tempguess[2] + dt * Smy0;
		Uguess[3] = pG->U[ks][j][i].M3;
		Uguess[4] = pG->U[ks][j][i].E + tempguess[4];

	

#ifdef CONS_GRAVITY
		/* density_old is now actually the updated density */
		Uguess[4] += 0.5*(pG->U[ks][j][i].d-grav_mean_rho)*pG->Phi_old[ks][j][i]-0.5*(density_old[j][i]-grav_mean_rho)*pG->Phi[ks][j][i];

#endif

			

	/*  Uguess[0] = d; Uguess[1]=Mx; Uguess[2]=My; Uguess[4]=E */


		/*====================================================================*/
		/* We need to add shearing box source in the guess solution */		

#ifdef SHEARING_BOX
			/* Initialize the shear source term to be zero */
			for(m=0; m<4; m++)
				ShearSource[m] = 0.0;

		
			cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		/* Store the current state */
			M1n  = pG->U[ks][j][i].M1;
#ifdef FARGO
			if (ShBoxCoord==xy) dM2n = pG->U[ks][j][i].M2;
			if (ShBoxCoord==xz) dM3n = pG->U[ks][j][i].M3;
#else
			if (ShBoxCoord==xy) dM2n = pG->U[ks][j][i].M2 + qom*x1*pG->U[ks][j][i].d;
			if (ShBoxCoord==xz) dM3n = pG->U[ks][j][i].M3 + qom*x1*pG->U[ks][j][i].d;
#endif
			
			/* Calculate the flux for the y-momentum fluctuation (M3 in 2D) */
			if (ShBoxCoord==xy){
				frx1_dM2 = x1Flux[j][i+1].My;
				flx1_dM2 = x1Flux[j][i  ].My;
				frx2_dM2 = x2Flux[j+1][i].Mx;
				flx2_dM2 = x2Flux[j  ][i].Mx;
			}
			if (ShBoxCoord==xz){
				frx1_dM3 = x1Flux[j][i+1].Mz;
				flx1_dM3 = x1Flux[j][i  ].Mz;
				frx2_dM3 = x2Flux[j+1][i].My;
				flx2_dM3 = x2Flux[j  ][i].My;
			}
#ifndef FARGO
			if (ShBoxCoord==xy){
				frx1_dM2 += qom*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
				flx1_dM2 += qom*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
				frx2_dM2 += qom*(x1            )*x2Flux[j+1][i].d;
				flx2_dM2 += qom*(x1            )*x2Flux[j  ][i].d;
			}
			if (ShBoxCoord==xz){
				frx1_dM3 += qom*(x1+0.5*pG->dx1)*x1Flux[j][i+1].d;
				flx1_dM3 += qom*(x1-0.5*pG->dx1)*x1Flux[j][i  ].d;
				frx2_dM3 += qom*(x1            )*x2Flux[j+1][i].d;
				flx2_dM3 += qom*(x1            )*x2Flux[j  ][i].d;
			}
#endif
			
			/* evolve M1n and dM3n by dt/2 using flux gradients */
			M1e = M1n - hdtodx1*(x1Flux[j][i+1].Mx - x1Flux[j][i].Mx)
			- hdtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
			
			if (ShBoxCoord==xy){
				dM2e = dM2n - hdtodx1*(frx1_dM2 - flx1_dM2)
				- hdtodx2*(frx2_dM2 - flx2_dM2);
			}
			if (ShBoxCoord==xz){
				dM3e = dM3n - hdtodx1*(frx1_dM3 - flx1_dM3)
				- hdtodx2*(frx2_dM3 - flx2_dM3);
			}

			/* Add radiation source term for half time step */
			
			M1e += 0.5 * tempguess[1];
			dM2e += 0.5 * tempguess[2];

			
			/* Update the 1- and 2-momentum (or 1- and 3-momentum in XZ 2D shearing box)
			 * for the Coriolis and tidal potential source terms using a Crank-Nicholson
			 * discretization for the momentum fluctuation equation. */
			
			if (ShBoxCoord==xy){
				ShearSource[0] = (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
				ShearSource[1] = 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;
#ifndef FARGO
				ShearSource[1] -=0.5*qshear*om_dt*(x1Flux[j][i].d+x1Flux[j][i+1].d);
#endif
			}
			if (ShBoxCoord==xz){
				ShearSource[0] = (4.0*dM3e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
				ShearSource[2] = 2.0*(qshear-2.)*(M1e + om_dt*dM3e)*fact;
#ifndef FARGO
				ShearSource[2] -=0.5*qshear*om_dt*(x1Flux[j][i].d+x1Flux[j][i+1].d);
#endif
			}
			
			/* Update the energy for a fixed potential.
			 * This update is identical to non-SHEARING_BOX below  */
			
			phic = (*ShearingBoxPot)((x1            ),x2,x3);
			phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
			phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
			ShearSource[3] -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
						 x1Flux[j][i+1].d*(phir - phic));
#endif
			
			phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
			phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
#ifndef BAROTROPIC
			ShearSource[3] -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
						 x2Flux[j+1][i].d*(phir - phic));
#endif


#endif

#ifdef SHEARING_BOX

			Uguess[1] += ShearSource[0];
			Uguess[2] += ShearSource[1];
			Uguess[3] += ShearSource[2];
			Uguess[4] += ShearSource[3];
#endif

	/* Now calculate the source term due to the guess solution */

		
		density = Uguess[0];
		pressure = (Uguess[4] - 0.5 * (Uguess[1] * Uguess[1] 
				+ Uguess[2] * Uguess[2] + Uguess[3] * Uguess[3]) / density) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		B1ch = 0.5*(    pG->B1i[ks][j][i] +     pG->B1i[ks][j][i+1]);
		B2ch = 0.5*(    pG->B2i[ks][j][i] +     pG->B2i[ks][j+1][i]);
		B3ch = pG->U[ks][j][i].B3c;

		pressure -= 0.5 * (B1ch * B1ch + B2ch * B2ch + B3ch * B3ch) * (Gamma - 1.0);
#endif

		if((pressure < TINY_NUMBER) || (pressure != pressure)){
			pressure = density * temperature * R_ideal;
			
			badcellflag = 1;

		}


		/* for bad cell, need to recalculate total energy */
		if(badcellflag){
			Uguess[4] =  pressure / (Gamma - 1.0) + 0.5 * (Uguess[1] * Uguess[1] 
				+ Uguess[2] * Uguess[2] + Uguess[3] * Uguess[3]) / density;
#ifdef RADIATION_MHD
                        Uguess[4] += 0.5 * (B1ch * B1ch + B2ch * B2ch + B3ch * B3ch);
#endif


		}	

		temperature = pressure / (density * R_ideal);
		Tguess = temperature;
		
		/* update source term */
		Usource.d = Uguess[0];
		Usource.Mx = Uguess[1];
		Usource.My = Uguess[2];
		Usource.E = Uguess[4];

		for(m=0; m<NOPACITY;m++){
			Usource.Sigma[m] = Sigma[m];
		}


		velocity_x = Uguess[1] / density;
		velocity_y = Uguess[2] / density;
				
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
		cc_pos(pG,i,j,ks,&x1,&x2,&x3);
		velocity_y -= qom * x1;		
#endif

	/* The Source term */
		if(pressure > TINY_NUMBER){
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, NULL, x1);
		}
		else{
			SEE = 0.0;			
			SErho = 0.0;
			SEmx = 0.0;	
			SEmy = 0.0;
		}



		Source_guess[1] = -Prat * (-(Sigma_aF + Sigma_sF) * (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
			+ pG->U[ks][j][i].Edd_21 * velocity_y)* pG->U[ks][j][i].Er / Crat));
		Source_guess[2] = -Prat * (-(Sigma_aF + Sigma_sF) * (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
			+ pG->U[ks][j][i].Edd_21 * velocity_x)* pG->U[ks][j][i].Er / Crat));
                
                
        /* Calculate the momentum source term before opacity is updated. It may cause trouble if opacity changes too much */
        Prwork2 = -Prat * (Sigma_aF - Sigma_sF) * (Source_Inv[1][1] * velocity_x
                 * (pG->U[ks][j][i].Fr1 - ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x
                 + pG->U[ks][j][i].Edd_21 * velocity_y) * pG->U[ks][j][i].Er / Crat)
                 + Source_Inv[2][2] * velocity_y
                 * (pG->U[ks][j][i].Fr2 - ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y
                  + pG->U[ks][j][i].Edd_21 * velocity_x) * pG->U[ks][j][i].Er / Crat));
        
		

		/* Prepare the Prims variable */
		Wopacity = Cons_to_Prim(&pG->U[ks][j][i]);
		/* Now update the density, pressure and velocity */
		Wopacity.d = density;
		Wopacity.P = pressure;
		Wopacity.V1 = velocity_x;
		Wopacity.V2 = velocity_y;
		/* background shearing should be included */


		if(Opacity != NULL){
			Opacity(&Wopacity, Sigma, NULL);

			Sigma_sF = Sigma[0];
			Sigma_aF = Sigma[1];
			Sigma_aP = Sigma[2];
			Sigma_aE = Sigma[3];
		}
		

		/* calculate the predict energy source term */
		ThermalRelaxation(temperature, pG->U[ks][j][i].Er, density, Sigma_aP, Sigma_aE, dt, NULL, &Ersource);
		Ersource = Ersource - pG->U[ks][j][i].Er;
				
#ifdef FLD
        pG->U[ks][j][i].Er += (1.0 - Eratio) * (0.5 * (Ersource + pG->Ersource[ks][j][i])) ;
           
        Uguess[4] += (-0.5 * Prat * (Ersource - pG->Ersource[ks][j][i]));
				
		pG->U[ks][j][i].d  = Uguess[0];
		pG->U[ks][j][i].M1 = Uguess[1];
		pG->U[ks][j][i].M2 = Uguess[2];			
		pG->U[ks][j][i].E  = Uguess[4];		
				
#else

		diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[ks][j][i].Er;

		

                
		Source_Invnew[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
                
                

		if(Erflag){
			Source_guess[4] = -Prat * Crat * diffTEr + Prwork2;
		}
		else{
			Source_guess[4] = Prwork2;
		}
		
				
				
/* Calculate the shearing source term due to the guess solution */				
		Errort[0] = 0.0;
		Errort[1] = hdt * (Source[1] + Source_guess[1]) 
						- dt * (divFlux1[1] + divFlux2[1]) 
						- dt * Source_Inv[1][1] * (Source[1] - divFlux1[1] - divFlux2[1]);
		
		Errort[2] = hdt * (Source[2] + Source_guess[2]) 
						- dt * (divFlux1[2] + divFlux2[2]) 
						- dt * Source_Inv[2][2] * (Source[2] - divFlux1[2] - divFlux2[2]);	

		for(m=1; m<3; m++) {
				tempguess[m] = Source_Inv[m][m] * Errort[m];
		}

		/* Correction to the energy source term needs special treatment */
		Errort[4] = hdt * (Source[4] + Source_guess[4]) - dt * (divFlux1[4] + divFlux2[4]) 
				- dt * Source_Inv[4][4] * (Source[4] - (divFlux1[4] + divFlux2[4]));

		tempguess[4] = Source_Invnew[4][4] * Errort[4];

		if(!Erflag){
			tempguess[4] += -0.5 * Prat * (Ersource - pG->Ersource[ks][j][i]);
		}



		Prworksource += Source_Invnew[4][4] * (hdt * (Prwork1 + Prwork2) - Prworksource);
		
		/* Estimate the added radiation source term  */
		
		if(Erflag){
		if(Prat > 0.0){
			pG->Ersource[ks][j][i] = Uguess[4] + tempguess[4] - 
				(pG->U[ks][j][i].E - dt * (divFlux1[4] + divFlux2[4]));

#ifdef SHEARING_BOX
			pG->Ersource[ks][j][i] -= ShearSource[3];

#endif

#ifdef CONS_GRAVITY
			pG->Ersource[ks][j][i] -= 0.5*(pG->U[ks][j][i].d-grav_mean_rho)*pG->Phi_old[ks][j][i]-0.5*(density_old[j][i]-grav_mean_rho)*pG->Phi[ks][j][i];
#endif

			/* Subtract the actual added work done by radiation force */
			/* This is added seperately for the radiation subsystem */
			
			pG->Ersource[ks][j][i] -= Prworksource;
				
			pG->Ersource[ks][j][i] /= -Prat;
			
		}
		else{
			
			pG->Ersource[ks][j][i] = 0.0;
			
		}
		}
		else{
			pG->Ersource[ks][j][i] = 0.5 * (pG->Ersource[ks][j][i] + Ersource);
		}
		
		pG->Eulersource[ks][j][i] = -Prworksource/Prat;

		if(badcellflag){

		/* Do not apply correction for bad cell */
			pG->Ersource[ks][j][i] = 0.0;
			pG->Eulersource[ks][j][i] = 0.0;
			for(n=0; n<5; n++){
				tempguess[n] = 0.0;
			}		

		}
		/* change due to radiation source term */
		/* This is used for shearing-box */
		pG->U[ks][j][i].d  = Uguess[0];
		pG->U[ks][j][i].M1 = Uguess[1] + tempguess[1];
		pG->U[ks][j][i].M2 = Uguess[2] + tempguess[2];			
		pG->U[ks][j][i].E  = Uguess[4] + tempguess[4];
				
				
#endif /* FLD */

	
/*
		if(!Erflag){
			pG->U[ks][j][i].d  = Uguess[0] + tempguess[0];
			pG->U[ks][j][i].M1 = Uguess[1] + tempguess[1];
			pG->U[ks][j][i].M2 = Uguess[2] + tempguess[2];			
			pG->U[ks][j][i].E  = Uguess[4];
			
		}
		else{
			
		}
		
*/		
		/*========================================================*/
		/* In 2D version , we always assume Fr_3 = 0. So we do not need source term for M3.*/
		/* otherwise, we have to go to 3D version of matrix to solve Fr_3 */
		pG->U[ks][j][i].M3 -= dtodx1 * (x1Flux[j][i+1].Mz - x1Flux[j][i].Mz);
		pG->U[ks][j][i].M3 -= dtodx2 * (x2Flux[j+1][i].My - x2Flux[j][i].My);
				
				
		}/* End big loop i */
	}/* End big loop j */



/* update the magnetic field */
/* magnetic field is independent of predict and correct step */
#ifdef RADIATION_MHD
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
	/* due to x1 flux */
	pG->U[ks][j][i].B2c -= dtodx1*(x1Flux[j][i+1].By - x1Flux[j][i].By);
      	pG->U[ks][j][i].B3c -= dtodx1*(x1Flux[j][i+1].Bz - x1Flux[j][i].Bz);


	/* due to x2 flux */
	pG->U[ks][j][i].B3c -= dtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
      	pG->U[ks][j][i].B1c -= dtodx2*(x2Flux[j+1][i].Bz - x2Flux[j][i].Bz);
    }
  }

/* set cell centered magnetic field to average of updated face centered field */
 for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {

      pG->U[ks][j][i].B1c =0.5*(pG->B1i[ks][j][i] + pG->B1i[ks][j][i+1]);
      pG->U[ks][j][i].B2c =0.5*(    pG->B2i[ks][j][i] +     pG->B2i[ks][j+1][i]);
/* Set the 3-interface magnetic field equal to the cell center field. */
      pG->B3i[ks][j][i] = pG->U[ks][j][i].B3c;
    }
  }


#endif /* radiation MHD part */



	/* calculate pG->Tguess after magnetic field is updated */
/*		for (j=js; j<=je; j++) {
      			for (i=is; i<=ie; i++) {

				pressure = pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 + pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2 + pG->U[ks][j][i].M3 * pG->U[ks][j][i].M3) / pG->U[ks][j][i].d;
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c);
#endif

				pressure *= (Gamma - 1.0);

				temperature = pressure / (pG->U[ks][j][i].d * R_ideal);
			
				pG->Tguess[ks][j][i] = pow(temperature, 4.0);
			
			}
		}
*/


/*=========================================================*/
/* Add non-stiff source term due to static gravity */
/*

  if (StaticGravPot != NULL){
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);


        pG->U[ks][j][i].M1 -= dtodx1*dhalf[j][i]*(phir-phil);


        pG->U[ks][j][i].E -= dtodx1*(x1Flux[j][i  ].d*(phic - phil) +
				 	x1Flux[j][i+1].d*(phir - phic));
			
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[ks][j][i].M2 -= dtodx2*dhalf[j][i]*(phir-phil);


        pG->U[ks][j][i].E -= dtodx2*(x2Flux[j  ][i].d*(phic - phil) +
                                     x2Flux[j+1][i].d*(phir - phic));

      }
    }
  }
*/


/* Check the pressure to make sure that it is positive */
	/* This will make the code more robust */
/*
		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++){
				density = pG->U[ks][j][i].d;
				pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 
				+ pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2 + pG->U[ks][j][i].M3 * pG->U[ks][j][i].M3) / density ) * (Gamma - 1);
			
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c) * (Gamma - 1.0);
#endif

				if(pressure < TINY_NUMBER){
					if(pG->U[ks][j][i].Er < TINY_NUMBER)
						pG->U[ks][j][i].Er = TINY_NUMBER;
					

				pressure = density * R_ideal * pow(pG->U[ks][j][i].Er,0.25);
				
				pG->U[ks][j][i].E = 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 + pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2 + pG->U[ks][j][i].M3 * pG->U[ks][j][i].M3) / density + pressure / (Gamma - 1.0);
#ifdef RADIATION_MHD
				pG->U[ks][j][i].E += 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c);
#endif	
				} 
				

			} 
		}
	
*/

		
	/* Boundary condition is applied in the main function */

	/* Update the opacity if Opacity function is set in the problem generator */
	if(Opacity != NULL){

		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++){

				Wopacity = Cons_to_Prim(&pG->U[ks][j][i]);
				
				/* Add background shearing */
#ifdef FARGO	
				cc_pos(pG,i,j,ks,&x1,&x2,&x3);
				Wopacity.V2 -= qom * x1;		
#endif
				
				if(Wopacity.P > TINY_NUMBER)
				{					
			
					Opacity(&Wopacity,Sigma,NULL);
						for(m=0;m<NOPACITY;m++){
							pG->U[ks][j][i].Sigma[m] = Sigma[m];
						}
				}

			}
		}
	}

  return;

 
}

/* ============================================================================================ *
 * ============================================================================================ *
 */



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

#ifdef RADIATION_MHD
  if ((emf3 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */


  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef RADIATION_MHD
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

#ifdef CONS_GRAVITY
  if ((x1Flux_grav   =(Real**)calloc_2d_array((size2+1),(size1+1),sizeof(Real)))==NULL)
    goto on_error;
  if ((x2Flux_grav   =(Real**)calloc_2d_array((size2+1),(size1+1),sizeof(Real)))==NULL)
    goto on_error;
  if ((density_old =(Real**)calloc_2d_array(size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif

  if ((dhalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;


  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct_2d:  Free temporary integration arrays */

void integrate_destruct_2d(void)
{
#ifdef RADIATION_MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif /* MHD */

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef RADIATION_MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

#ifdef CONS_GRAVITY
  if (x1Flux_grav    != NULL) free_2d_array(x1Flux_grav);
  if (x2Flux_grav    != NULL) free_2d_array(x2Flux_grav);
  if (density_old    != NULL) free_2d_array(density_old);

#endif

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

#ifdef RADIATION_MHD
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
