#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_3d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit CTU integrator of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *      B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically thin cooling,
 *   shearing box source terms, and the H-correction of Sanders et al.
 *     For adb hydro, requires (9*Cons1DS +  3*Real) = 48 3D arrays
 *     For adb mhd, requires   (9*Cons1DS + 10*Real) = 73 3D arrays
 *   The H-correction of Sanders et al. adds another 3 arrays.  
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport in three dimensions", JCP, 227, 4123 (2008)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 *   J.M. Stone et al., "Athena: A new code for astrophysical MHD", ApJS,
 *   178, 137 (2008)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_3d_radMHD()
 *   integrate_init_3d()
 *   integrate_destruct_3d()
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
static Cons1DS ***Ul_x1Face=NULL, ***Ur_x1Face=NULL;
static Cons1DS ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
static Cons1DS ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;

/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static Real **remapEyiib=NULL, **remapEyoib=NULL;
#endif

/* The interface magnetic fields and emfs */
#ifdef RADIATION_MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real ***dhalf = NULL, ***phalf=NULL;


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - the upwind CT method in GS05, for emf1
 *   integrate_emf2_corner() - the upwind CT method in GS05, for emf2
 *   integrate_emf3_corner() - the upwind CT method in GS05, for emf3
 *============================================================================*/

#ifdef RADIATION_MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif /* MHD */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/


void integrate_3d_radMHD(DomainS *pD)
{
	GridS *pG=(pD->Grid);
	Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;
	Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2, hdtodx3 = 0.5 * dtodx3;
	Real dx1 = pG->dx1, dx2 = pG->dx2, dx3 = pG->dx3;
	Real hdt = 0.5*pG->dt, dt = pG->dt;
  	int i,il,iu,is=pG->is, ie=pG->ie;
  	int j,jl,ju,js=pG->js, je=pG->je;
  	int k, kl, ku, ks=pG->ks, ke=pG->ke; 
	int m, n;

	/* For static gravitational potential */
	Real x1, x2, x3, phicl, phicr, phifc, phic, phir, phil;

	Real Bx = 0.0;
	Real M1h, M2h, M3h;
#ifdef RADIATION_MHD
  	Real MHD_src_By,MHD_src_Bz,mdb1,mdb2,mdb3;
  	Real db1,db2,db3,l1,l2,l3,B1,B2,B3,V1,V2,V3;
  	Real B1ch,B2ch,B3ch;
#endif

#if defined(RADIATION_MHD) || defined(SELF_GRAVITY)
  	Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;
#endif
	
/* Variables for shearing box */	
	
#ifdef SHEARING_BOX
	int my_iproc,my_jproc,my_kproc;
	Real M1n, dM2n; /* M1, dM2=(My+d*q*Omega_0*x) at time n */
	Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
	Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
	Real fact, qom, om_dt = Omega_0*pG->dt;
/*	Real ShearingSource_M1, ShearingSource_M2, ShearingSource_M3, ShearingSource_E;
	Real Shearingguess_M1, Shearingguess_M2, Shearingguess_M3, Shearingguess_E;
 */
#endif /* SHEARING_BOX */

	Real temperature, velocity_x, velocity_y, velocity_z, velocity, pressure, density;
	Real Fr0x, Fr0y, Fr0z, diffTEr; /* co-moving flux, Fr-vf E_r/C , diffTEr = T^4 - Er*/
	Real Sigma_t, Sigma_a, Sigma_s;
	Real SPP, alpha, Propa_44, SEE, SErho, SEmx, SEmy, SEmz, dSigmadP;
	Real dSigma[4];
	Cons1DS Usource;
	/* for source term */

	/* In case momentum becomes stiff */
	Real SFmx, SFmy, SFmz, SVVx, SVVy, SVVz, betax, betay, betaz;


	Real Source_Inv[NVAR][NVAR], tempguess[NVAR], Uguess[NVAR], Source[NVAR], Source_guess[NVAR], Errort[NVAR], SourceFlux[NVAR];
	Real divFlux1[NVAR], divFlux2[NVAR], divFlux3[NVAR];

	/* SourceFlux is used to calculate flux due to non-stiff source terms, such as gravity */

	/* Initialize them to be zero */
	/* NVAR is the same as normal hydro code. Rad variables are not included here */
	for(i=0; i<NVAR; i++){
		Source[i] = 0.0;
		SourceFlux[i] = 0.0;
		Source_guess[i] = 0.0;
		Errort[i] = 0.0;
		divFlux1[i] = 0.0;
		divFlux2[i] = 0.0;
		divFlux3[i] = 0.0;
		tempguess[i] = 0.0;
		for(j=0; j<NVAR; j++) {
			Source_Inv[i][j] = 0.0;			
			
			if(i==j) {
		 		Source_Inv[i][j] = 1.0;
		
			}
		}
	}

	for(i=0; i<4; i++)
		dSigma[i] = 0.0;


	il = is - 2;
  	iu = ie + 2;
  	jl = js - 2;
  	ju = je + 2;
  	kl = ks - 2;
  	ku = ke + 2;

/*=== Step 1: Backward Euler is done in the main function for the whole mesh ===*/
	
	
/*=== STEP 2: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */
	for (k=kl; k<=ku; k++){
		for (j=jl; j<=ju; j++) {
    			for (i=is-nghost; i<=ie+nghost; i++) {
      			
			U1d[i].d  = pG->U[k][j][i].d;
      			U1d[i].Mx = pG->U[k][j][i].M1;
      			U1d[i].My = pG->U[k][j][i].M2;
      			U1d[i].Mz = pG->U[k][j][i].M3;
      			U1d[i].E  = pG->U[k][j][i].E;
			U1d[i].Er  = pG->U[k][j][i].Er;
    			U1d[i].Fr1  = pG->U[k][j][i].Fr1;
    			U1d[i].Fr2  = pG->U[k][j][i].Fr2;
    			U1d[i].Fr3  = pG->U[k][j][i].Fr3;
			U1d[i].Edd_11  = pG->U[k][j][i].Edd_11;
			U1d[i].Edd_21  = pG->U[k][j][i].Edd_21;
			U1d[i].Edd_22  = pG->U[k][j][i].Edd_22;
			U1d[i].Edd_31  = pG->U[k][j][i].Edd_31;
			U1d[i].Edd_32  = pG->U[k][j][i].Edd_32;
			U1d[i].Edd_33  = pG->U[k][j][i].Edd_33;
			U1d[i].Sigma_a  = pG->U[k][j][i].Sigma_a;
                        U1d[i].Sigma_t  = pG->U[k][j][i].Sigma_t;
			

#ifdef RADIATION_MHD
      			U1d[i].By = pG->U[k][j][i].B2c;
      			U1d[i].Bz = pG->U[k][j][i].B3c;
      			Bxc[i] = pG->U[k][j][i].B1c;
      			Bxi[i] = pG->B1i[k][j][i];
      			B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */
			}


/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

			for (i=is-nghost; i<=ie+nghost; i++) {
      				W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
    			}

			lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,3);


/*------Step 2c: Add radiation source terms to the left and right state for 0.5*dt--------*/

		for (i=il+1; i<=iu; i++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			pressure = W[i-1].P;
			temperature = pressure / (U1d[i-1].d * R_ideal);
			
			diffTEr = pow(temperature, 4.0) - U1d[i-1].Er;	


			velocity_x = U1d[i-1].Mx / U1d[i-1].d;
			velocity_y = U1d[i-1].My / U1d[i-1].d;
			velocity_z = U1d[i-1].Mz / U1d[i-1].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k][j][i-1].Sigma_t;
			Sigma_a	   = pG->U[k][j][i-1].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;

			/* co-moving flux */
			Fr0x = U1d[i-1].Fr1 - ((1.0 + U1d[i-1].Edd_11) * velocity_x + U1d[i-1].Edd_21 * velocity_y + U1d[i-1].Edd_31 * velocity_z) * U1d[i-1].Er / Crat;
			Fr0y = U1d[i-1].Fr2 - ((1.0 + U1d[i-1].Edd_22) * velocity_y + U1d[i-1].Edd_21 * velocity_x + U1d[i-1].Edd_32 * velocity_z) * U1d[i-1].Er / Crat;
			Fr0z = U1d[i-1].Fr3 - ((1.0 + U1d[i-1].Edd_33) * velocity_z + U1d[i-1].Edd_31 * velocity_x + U1d[i-1].Edd_32 * velocity_y) * U1d[i-1].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[i-1].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[i-1].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[i-1].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
					-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[i-1].d;			


			if(Opacity != NULL) Opacity(U1d[i-1].d, temperature, NULL, NULL, dSigma);

			/* dSigma[0] = dSigmt/drho, dSigma[1] = dSigma/drho, dSigma[2]=dSigmat/dT, dsigma[3]=dSigmaa/dT */
		
			dSigmadP =  dSigma[3] / (U1d[i-1].d * R_ideal); 
			
			
			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[i-1].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[i-1].Edd_11) * W[i-1].Er + Sigma_a * diffTEr) / (W[i-1].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[i-1].Edd_22) * W[i-1].Er + Sigma_a * diffTEr) / (W[i-1].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[i-1].Edd_33) * W[i-1].Er + Sigma_a * diffTEr) / (W[i-1].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
		/*===========================================================================*/
	
			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
			Propa_44 = alpha;

			Wl[i].Vx += dt * Source[1] * 0.5 * betax;
			Wl[i].Vy += dt * Source[2] * 0.5 * betay;
			Wl[i].Vz += dt * Source[3] * 0.5 * betaz;
			Wl[i].P += dt * Propa_44 * Source[4] * 0.5;
	
			Wl[i].Sigma_a = Sigma_a;
			Wl[i].Sigma_t = Sigma_t;
			Wl[i].Edd_11 = W[i-1].Edd_11;
			Wl[i].Edd_21 = W[i-1].Edd_21;
			Wl[i].Edd_22 = W[i-1].Edd_22;
			Wl[i].Edd_31 = W[i-1].Edd_31;
			Wl[i].Edd_32 = W[i-1].Edd_32;
			Wl[i].Edd_33 = W[i-1].Edd_33;

		/* For the right state */
	
	
			pressure = W[i].P;
			temperature = pressure / (U1d[i].d * R_ideal);

			diffTEr = pow(temperature, 4.0) - U1d[i].Er;	

			velocity_x = U1d[i].Mx / U1d[i].d;
			velocity_y = U1d[i].My / U1d[i].d;
			velocity_z = U1d[i].Mz / U1d[i].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k][j][i].Sigma_t;
			Sigma_a	   = pG->U[k][j][i].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;


			/* co-moving flux */
			Fr0x = U1d[i].Fr1 - ((1.0 + U1d[i].Edd_11) * velocity_x + U1d[i].Edd_21 * velocity_y + U1d[i].Edd_31 * velocity_z) * U1d[i].Er / Crat;
			Fr0y = U1d[i].Fr2 - ((1.0 + U1d[i].Edd_22) * velocity_y + U1d[i].Edd_21 * velocity_x + U1d[i].Edd_32 * velocity_z) * U1d[i].Er / Crat;
			Fr0z = U1d[i].Fr3 - ((1.0 + U1d[i].Edd_33) * velocity_z + U1d[i].Edd_31 * velocity_x + U1d[i].Edd_32 * velocity_y) * U1d[i].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[i].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[i].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[i].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
					-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[i].d;			

			
			if(Opacity != NULL) Opacity(U1d[i].d, temperature, NULL, NULL, dSigma);
		
			dSigmadP =  dSigma[3] / (U1d[i].d * R_ideal);

			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[i].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
	
			/*===================================================================*/
			/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[i].Edd_11) * W[i].Er + Sigma_a * diffTEr) / (W[i].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[i].Edd_22) * W[i].Er + Sigma_a * diffTEr) / (W[i].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[i].Edd_33) * W[i].Er + Sigma_a * diffTEr) / (W[i].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
			/*===========================================================================*/
	


			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
			Propa_44 = alpha;

			Wr[i].Vx += dt * Source[1] * 0.5 * betax;
			Wr[i].Vy += dt * Source[2] * 0.5 * betay;
			Wr[i].Vz += dt * Source[3] * 0.5 * betaz;
			Wr[i].P += dt * Propa_44 * Source[4] * 0.5;

			Wr[i].Sigma_a = Sigma_a;
			Wr[i].Sigma_t = Sigma_t;	
			Wr[i].Edd_11 = W[i].Edd_11;
			Wr[i].Edd_21 = W[i].Edd_21;
			Wr[i].Edd_22 = W[i].Edd_22;
			Wr[i].Edd_31 = W[i].Edd_31;
			Wr[i].Edd_32 = W[i].Edd_32;
			Wr[i].Edd_33 = W[i].Edd_33;		
			
		}

/*---------For radiation MHD, source term due to magnetic field is also added---------*/
		
#ifdef RADIATION_MHD
      		for (i=il+1; i<=iu; i++) {
/* Source terms for left states in zone i-1 */

	        db1 = (    pG->B1i[k  ][j  ][i  ] -     pG->B1i[k][j][i-1])*dx1i;
       		db2 = (    pG->B2i[k  ][j+1][i-1] -     pG->B2i[k][j][i-1])*dx2i;
        	db3 = (    pG->B3i[k+1][j  ][i-1] -     pG->B3i[k][j][i-1])*dx3i;

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

	        MHD_src_By = (pG->U[k][j][i-1].M2/pG->U[k][j][i-1].d)*l2;
	        MHD_src_Bz = (pG->U[k][j][i-1].M3/pG->U[k][j][i-1].d)*l3;

	        Wl[i].By += hdt*MHD_src_By;
	        Wl[i].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone i */

	        db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        	db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        	db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;

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

        	MHD_src_By = (pG->U[k][j][i].M2/pG->U[k][j][i].d)*l2;
        	MHD_src_Bz = (pG->U[k][j][i].M3/pG->U[k][j][i].d)*l3;

        	Wr[i].By += hdt*MHD_src_By;
        	Wr[i].Bz += hdt*MHD_src_Bz;
      		}
#endif

/*-------Add source term due to static gravitational potential for left and right state--------*/
	
      		if (StaticGravPot != NULL){
        		for (i=il+1; i<=iu; i++) {
         			 cc_pos(pG,i,j,k,&x1,&x2,&x3);

	        	phicr = (*StaticGravPot)( x1             ,x2,x3);
          		phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
          		phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          		Wl[i].Vx -= dtodx1*(phifc - phicl);
          		Wr[i].Vx -= dtodx1*(phicr - phifc);

        		}
      		}
			
/* Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
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
			
#endif /* SHEARING_BOX */				


/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
	      	for (i=il+1; i<=iu; i++) {
        		Ul_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);
        		Ur_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

#ifdef RADIATION_MHD
        		Bx = B1_x1Face[k][j][i];
#endif
			/* take the time step with x1Flux */
			x1Flux[k][j][i].d = dt;	
				x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);
      		}
		

		} /* End j */

	}/* End k */	


/*=== STEP 3: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/
/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 * NOTICE that momentum is rotated but flux and Eddington tensor are not *
 */
	for (k=kl; k<=ku; k++){
		for (i=il; i<=iu; i++) {
    			for (j=js-nghost; j<=je+nghost; j++) {
      				U1d[j].d  = pG->U[k][j][i].d;
      				U1d[j].Mx = pG->U[k][j][i].M2;
      				U1d[j].My = pG->U[k][j][i].M3;
      				U1d[j].Mz = pG->U[k][j][i].M1;
      				U1d[j].E  = pG->U[k][j][i].E;
				U1d[j].Er  = pG->U[k][j][i].Er;
    				U1d[j].Fr1  = pG->U[k][j][i].Fr1;
				U1d[j].Fr2  = pG->U[k][j][i].Fr2;
    				U1d[j].Fr3  = pG->U[k][j][i].Fr3;
				U1d[j].Edd_11  = pG->U[k][j][i].Edd_11;
				U1d[j].Edd_21  = pG->U[k][j][i].Edd_21;
				U1d[j].Edd_22  = pG->U[k][j][i].Edd_22;
				U1d[j].Edd_31  = pG->U[k][j][i].Edd_31;
				U1d[j].Edd_32  = pG->U[k][j][i].Edd_32;
				U1d[j].Edd_33  = pG->U[k][j][i].Edd_33;
				U1d[j].Sigma_a  = pG->U[k][j][i].Sigma_a;
              	          	U1d[j].Sigma_t  = pG->U[k][j][i].Sigma_t;
#ifdef RADIATION_MHD
      				U1d[j].By = pG->U[k][j][i].B3c;
      				U1d[j].Bz = pG->U[k][j][i].B1c;
      				Bxc[j] = pG->U[k][j][i].B2c;
      				Bxi[j] = pG->B2i[k][j][i];
      				B2_x2Face[k][j][i] = pG->B2i[k][j][i];
#endif /* MHD */
			}


		for (j=js-nghost; j<=je+nghost; j++) {
        		W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx2,jl+1,ju-1,Wl,Wr,3);

/*------Step 3b: Add radiation source terms to the left and right state for 0.5*dt--------*/


		for (j=jl+1; j<=ju; j++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			pressure = W[j-1].P;
			temperature = pressure / (U1d[j-1].d * R_ideal);
			
			diffTEr = pow(temperature, 4.0) - U1d[j-1].Er;	


			velocity_x = U1d[j-1].Mz / U1d[j-1].d;
			velocity_y = U1d[j-1].Mx / U1d[j-1].d;
			velocity_z = U1d[j-1].My / U1d[j-1].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k][j-1][i].Sigma_t;
			Sigma_a	   = pG->U[k][j-1][i].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;

			/* co-moving flux */
			Fr0x = U1d[j-1].Fr1 - ((1.0 + U1d[j-1].Edd_11) * velocity_x + U1d[j-1].Edd_21 * velocity_y + U1d[j-1].Edd_31 * velocity_z) * U1d[j-1].Er / Crat;
			Fr0y = U1d[j-1].Fr2 - ((1.0 + U1d[j-1].Edd_22) * velocity_y + U1d[j-1].Edd_21 * velocity_x + U1d[j-1].Edd_32 * velocity_z) * U1d[j-1].Er / Crat;
			Fr0z = U1d[j-1].Fr3 - ((1.0 + U1d[j-1].Edd_33) * velocity_z + U1d[j-1].Edd_31 * velocity_x + U1d[j-1].Edd_32 * velocity_y) * U1d[j-1].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[j-1].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[j-1].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[j-1].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
			-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[j-1].d;			


			if(Opacity != NULL) Opacity(U1d[j-1].d, temperature, NULL, NULL, dSigma);

			/* dSigma[0] = dSigmt/drho, dSigma[1] = dSigma/drho, dSigma[2]=dSigmat/dT, dsigma[3]=dSigmaa/dT */
		
			dSigmadP =  dSigma[3] / (U1d[j-1].d * R_ideal); 
			
			
			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[j-1].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[j-1].Edd_11) * W[j-1].Er + Sigma_a * diffTEr) / (W[j-1].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[j-1].Edd_22) * W[j-1].Er + Sigma_a * diffTEr) / (W[j-1].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[j-1].Edd_33) * W[j-1].Er + Sigma_a * diffTEr) / (W[j-1].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
		/*===========================================================================*/
	
			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
			Propa_44 = alpha;

			/* "Vx" is actually vy, "vz" is actually vx and "vy" is actually vz. */
			/* We stay with the correct meaning in source terms */

			Wl[j].Vx += dt * Source[2] * 0.5 * betay;
			Wl[j].Vy += dt * Source[3] * 0.5 * betaz;
			Wl[j].Vz += dt * Source[1] * 0.5 * betax;
			Wl[j].P += dt * Propa_44 * Source[4] * 0.5;
	
			Wl[j].Sigma_a = Sigma_a;
			Wl[j].Sigma_t = Sigma_t;
			Wl[j].Edd_11 = W[j-1].Edd_11;
			Wl[j].Edd_21 = W[j-1].Edd_21;
			Wl[j].Edd_22 = W[j-1].Edd_22;
			Wl[j].Edd_31 = W[j-1].Edd_31;
			Wl[j].Edd_32 = W[j-1].Edd_32;
			Wl[j].Edd_33 = W[j-1].Edd_33;

		/* For the right state */
	
	
			pressure = W[j].P;
			temperature = pressure / (U1d[j].d * R_ideal);

			diffTEr = pow(temperature, 4.0) - U1d[j].Er;	

			velocity_x = U1d[j].Mz / U1d[j].d;
			velocity_y = U1d[j].Mx / U1d[j].d;
			velocity_z = U1d[j].My / U1d[j].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k][j][i].Sigma_t;
			Sigma_a	   = pG->U[k][j][i].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;


			/* co-moving flux */
			Fr0x = U1d[j].Fr1 - ((1.0 + U1d[j].Edd_11) * velocity_x + U1d[j].Edd_21 * velocity_y + U1d[j].Edd_31 * velocity_z) * U1d[j].Er / Crat;
			Fr0y = U1d[j].Fr2 - ((1.0 + U1d[j].Edd_22) * velocity_y + U1d[j].Edd_21 * velocity_x + U1d[j].Edd_32 * velocity_z) * U1d[j].Er / Crat;
			Fr0z = U1d[j].Fr3 - ((1.0 + U1d[j].Edd_33) * velocity_z + U1d[j].Edd_31 * velocity_x + U1d[j].Edd_32 * velocity_y) * U1d[j].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[j].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[j].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[j].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
			-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[j].d;			

			
			if(Opacity != NULL) Opacity(U1d[j].d, temperature, NULL, NULL, dSigma);
		
			dSigmadP =  dSigma[3] / (U1d[j].d * R_ideal);

			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[j].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
	
			/*===================================================================*/
			/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[j].Edd_11) * W[j].Er + Sigma_a * diffTEr) / (W[j].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[j].Edd_22) * W[j].Er + Sigma_a * diffTEr) / (W[j].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[j].Edd_33) * W[j].Er + Sigma_a * diffTEr) / (W[j].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
			/*===========================================================================*/
	


			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
			Propa_44 = alpha;

			Wr[j].Vx += dt * Source[2] * 0.5 * betay;
			Wr[j].Vy += dt * Source[3] * 0.5 * betaz;
			Wr[j].Vz += dt * Source[1] * 0.5 * betax;
			Wr[j].P += dt * Propa_44 * Source[4] * 0.5;

			Wr[j].Sigma_a = Sigma_a;
			Wr[j].Sigma_t = Sigma_t;	
			Wr[j].Edd_11 = W[j].Edd_11;
			Wr[j].Edd_21 = W[j].Edd_21;
			Wr[j].Edd_22 = W[j].Edd_22;
			Wr[j].Edd_31 = W[j].Edd_31;
			Wr[j].Edd_32 = W[j].Edd_32;
			Wr[j].Edd_33 = W[j].Edd_33;		
			
		}


#ifdef RADIATION_MHD

	      	for (j=jl+1; j<=ju; j++) {
/* Source terms for left states in zone j-1 */
        		db1 = (    pG->B1i[k  ][j-1][i+1] -     pG->B1i[k][j-1][i])*dx1i;
        		db2 = (    pG->B2i[k  ][j  ][i  ] -     pG->B2i[k][j-1][i])*dx2i;
        		db3 = (    pG->B3i[k+1][j-1][i  ] -     pG->B3i[k][j-1][i])*dx3i;

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

			MHD_src_By = (pG->U[k][j-1][i].M3/pG->U[k][j-1][i].d)*l3;
			MHD_src_Bz = (pG->U[k][j-1][i].M1/pG->U[k][j-1][i].d)*l1;

       	 		Wl[j].By += hdt*MHD_src_By;
        		Wl[j].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone j */
        		db1 = (    pG->B1i[k  ][j  ][i+1] -    pG->B1i[k][j][i])*dx1i;
        		db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        		db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;

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

        		MHD_src_By = (pG->U[k][j][i].M3/pG->U[k][j][i].d)*l3;
        		MHD_src_Bz = (pG->U[k][j][i].M1/pG->U[k][j][i].d)*l1;

        		Wr[j].By += hdt*MHD_src_By;
        		Wr[j].Bz += hdt*MHD_src_Bz;
      		}	
#endif


		
      		if (StaticGravPot != NULL){
        		for (j=jl+1; j<=ju; j++) {
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);

          			phicr = (*StaticGravPot)(x1, x2             ,x3);
         	 		phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          			phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          		 	Wl[j].Vx -= dtodx2*(phifc - phicl);
          			Wr[j].Vx -= dtodx2*(phicr - phifc);
        	
			}
      		}

/*-----Calculate flux along x2 direction ---------*/

		for (j=jl+1; j<=ju; j++) {
        		Ul_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);
        		Ur_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);

#ifdef RADIATION_MHD
        		Bx = B2_x2Face[k][j][i];
#endif
			/* Take the time step with x2Flux[][][].d */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[j],Wr[j],Bx,&x2Flux[k][j][i]);
      		}

		}/* END i */
	}/* end k */	




/*=== STEP 4: Compute L/R x3-interface states and 1D x3-Fluxes ===============*/
/*--- Step 4a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M3, M1, M2, E, B1c, B2c, s[n])
 * NOTICE that momentum is rotated but flux and Eddington tensor are not *
 */
	for (j=jl; j<=ju; j++){
		for (i=il; i<=iu; i++) {
    			for (k=ks-nghost; k<=ke+nghost; k++) {
      				U1d[k].d  = pG->U[k][j][i].d;
      				U1d[k].Mx = pG->U[k][j][i].M3;
      				U1d[k].My = pG->U[k][j][i].M1;
      				U1d[k].Mz = pG->U[k][j][i].M2;
      				U1d[k].E  = pG->U[k][j][i].E;
				U1d[k].Er  = pG->U[k][j][i].Er;
    				U1d[k].Fr1  = pG->U[k][j][i].Fr1;
				U1d[k].Fr2  = pG->U[k][j][i].Fr2;
    				U1d[k].Fr3  = pG->U[k][j][i].Fr3;
				U1d[k].Edd_11  = pG->U[k][j][i].Edd_11;
				U1d[k].Edd_21  = pG->U[k][j][i].Edd_21;
				U1d[k].Edd_22  = pG->U[k][j][i].Edd_22;
				U1d[k].Edd_31  = pG->U[k][j][i].Edd_31;
				U1d[k].Edd_32  = pG->U[k][j][i].Edd_32;
				U1d[k].Edd_33  = pG->U[k][j][i].Edd_33;
				U1d[k].Sigma_a  = pG->U[k][j][i].Sigma_a;
              	          	U1d[k].Sigma_t  = pG->U[k][j][i].Sigma_t;
#ifdef RADIATION_MHD
      				U1d[k].By = pG->U[k][j][i].B1c;
      				U1d[k].Bz = pG->U[k][j][i].B2c;
      				Bxc[k] = pG->U[k][j][i].B3c;
      				Bxi[k] = pG->B3i[k][j][i];
      				B3_x3Face[k][j][i] = pG->B3i[k][j][i];
#endif /* MHD */
			}


		for (k=ks-nghost; k<=ke+nghost; k++) {
        		W[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx3,kl+1,ku-1,Wl,Wr,3);

/*------Step 3b: Add radiation source terms to the left and right state for 0.5*dt--------*/


		for (k=kl+1; k<=ku; k++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			pressure = W[k-1].P;
			temperature = pressure / (U1d[k-1].d * R_ideal);
			
			diffTEr = pow(temperature, 4.0) - U1d[k-1].Er;	


			velocity_x = U1d[k-1].My / U1d[k-1].d;
			velocity_y = U1d[k-1].Mz / U1d[k-1].d;
			velocity_z = U1d[k-1].Mx / U1d[k-1].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k-1][j][i].Sigma_t;
			Sigma_a	   = pG->U[k-1][j][i].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;

			/* co-moving flux */
			Fr0x = U1d[k-1].Fr1 - ((1.0 + U1d[k-1].Edd_11) * velocity_x + U1d[k-1].Edd_21 * velocity_y + U1d[k-1].Edd_31 * velocity_z) * U1d[k-1].Er / Crat;
			Fr0y = U1d[k-1].Fr2 - ((1.0 + U1d[k-1].Edd_22) * velocity_y + U1d[k-1].Edd_21 * velocity_x + U1d[k-1].Edd_32 * velocity_z) * U1d[k-1].Er / Crat;
			Fr0z = U1d[k-1].Fr3 - ((1.0 + U1d[k-1].Edd_33) * velocity_z + U1d[k-1].Edd_31 * velocity_x + U1d[k-1].Edd_32 * velocity_y) * U1d[k-1].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[k-1].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[k-1].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[k-1].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
			-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[k-1].d;			


			if(Opacity != NULL) Opacity(U1d[k-1].d, temperature, NULL, NULL, dSigma);

			/* dSigma[0] = dSigmt/drho, dSigma[1] = dSigma/drho, dSigma[2]=dSigmat/dT, dsigma[3]=dSigmaa/dT */
		
			dSigmadP =  dSigma[3] / (U1d[k-1].d * R_ideal); 
			
			
			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[k-1].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[k-1].Edd_11) * W[k-1].Er + Sigma_a * diffTEr) / (W[k-1].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[k-1].Edd_22) * W[k-1].Er + Sigma_a * diffTEr) / (W[k-1].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[k-1].Edd_33) * W[k-1].Er + Sigma_a * diffTEr) / (W[k-1].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
		/*===========================================================================*/
	
			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d; */
			Propa_44 = alpha;

			/* "Vx" is actually vz, "vz" is actually vy and "vy" is actually vx. */
			/* We stay with the correct meaning in source terms */

			Wl[k].Vx += dt * Source[3] * 0.5 * betaz;
			Wl[k].Vy += dt * Source[1] * 0.5 * betax;
			Wl[k].Vz += dt * Source[2] * 0.5 * betay;
			Wl[k].P += dt * Propa_44 * Source[4] * 0.5;
	
			Wl[k].Sigma_a = Sigma_a;
			Wl[k].Sigma_t = Sigma_t;
			Wl[k].Edd_11 = W[k-1].Edd_11;
			Wl[k].Edd_21 = W[k-1].Edd_21;
			Wl[k].Edd_22 = W[k-1].Edd_22;
			Wl[k].Edd_31 = W[k-1].Edd_31;
			Wl[k].Edd_32 = W[k-1].Edd_32;
			Wl[k].Edd_33 = W[k-1].Edd_33;

		/* For the right state */
	
	
			pressure = W[k].P;
			temperature = pressure / (U1d[k].d * R_ideal);

			diffTEr = pow(temperature, 4.0) - U1d[k].Er;	

			velocity_x = U1d[k].My / U1d[k].d;
			velocity_y = U1d[k].Mz / U1d[k].d;
			velocity_z = U1d[k].Mx / U1d[k].d;
			velocity   = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			Sigma_t    = pG->U[k][j][i].Sigma_t;
			Sigma_a	   = pG->U[k][j][i].Sigma_a;
			Sigma_s	   = Sigma_t - Sigma_a;


			/* co-moving flux */
			Fr0x = U1d[k].Fr1 - ((1.0 + U1d[k].Edd_11) * velocity_x + U1d[k].Edd_21 * velocity_y + U1d[k].Edd_31 * velocity_z) * U1d[k].Er / Crat;
			Fr0y = U1d[k].Fr2 - ((1.0 + U1d[k].Edd_22) * velocity_y + U1d[k].Edd_21 * velocity_x + U1d[k].Edd_32 * velocity_z) * U1d[k].Er / Crat;
			Fr0z = U1d[k].Fr3 - ((1.0 + U1d[k].Edd_33) * velocity_z + U1d[k].Edd_31 * velocity_x + U1d[k].Edd_32 * velocity_y) * U1d[k].Er / Crat;

			/* Source term for velocity , not momentum*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat) / U1d[k].d;
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat) / U1d[k].d;
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat) / U1d[k].d;
			
			/* Source term for pressure */
			Source[4] = -(Gamma - 1.0) * Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat)
			-(Gamma - 1.0) * (velocity_x * Source[1] + velocity_y * Source[2] + velocity_z * Source[3]) * U1d[k].d;			

			
			if(Opacity != NULL) Opacity(U1d[k].d, temperature, NULL, NULL, dSigma);
		
			dSigmadP =  dSigma[3] / (U1d[k].d * R_ideal);

			SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * pow(temperature, 3.0) * (1.0 - velocity * velocity/(Crat * Crat)) /(U1d[k].d * R_ideal)
				-(Gamma - 1.0) * Prat * Crat * diffTEr * dSigmadP * (1.0 - velocity * velocity/(Crat * Crat))
			      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
	
			/*===================================================================*/
			/* In case velocity is large, momentum source term is also stiff */
			SVVx = -Prat * (Sigma_t * (1.0 + W[k].Edd_11) * W[k].Er + Sigma_a * diffTEr) / (W[k].d * Crat);
		
			if(fabs(SVVx * dt * 0.5) > 0.001)
			betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
			else 
			betax = 1.0 + 0.25 * SVVx * dt;

			SVVy = -Prat * (Sigma_t * (1.0 + W[k].Edd_22) * W[k].Er + Sigma_a * diffTEr) / (W[k].d * Crat);
		
			if(fabs(SVVy * dt * 0.5) > 0.001)
			betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
			else 
			betay = 1.0 + 0.25 * SVVy * dt;

			SVVz = -Prat * (Sigma_t * (1.0 + W[k].Edd_33) * W[k].Er + Sigma_a * diffTEr) / (W[k].d * Crat);
		
			if(fabs(SVVz * dt * 0.5) > 0.001)
			betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
			else 
			betaz = 1.0 + 0.25 * SVVz * dt;
			/*===========================================================================*/
	


			if(fabs(SPP * dt * 0.5) > 0.001)
			alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
			else 
			alpha = 1.0 + 0.25 * SPP * dt;
			/* In case SPP * dt  is small, use expansion expression */	
			/* Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d; */
			Propa_44 = alpha;

			Wr[k].Vx += dt * Source[3] * 0.5 * betaz;
			Wr[k].Vy += dt * Source[1] * 0.5 * betax;
			Wr[k].Vz += dt * Source[2] * 0.5 * betay;
			Wr[k].P += dt * Propa_44 * Source[4] * 0.5;

			Wr[k].Sigma_a = Sigma_a;
			Wr[k].Sigma_t = Sigma_t;	
			Wr[k].Edd_11 = W[k].Edd_11;
			Wr[k].Edd_21 = W[k].Edd_21;
			Wr[k].Edd_22 = W[k].Edd_22;
			Wr[k].Edd_31 = W[k].Edd_31;
			Wr[k].Edd_32 = W[k].Edd_32;
			Wr[k].Edd_33 = W[k].Edd_33;		
			
		}


#ifdef RADIATION_MHD

 		for (k=kl+1; k<=ku; k++) {
/* Source terms for left states in zone k-1 */
        		db1 = (    pG->B1i[k-1][j  ][i+1] -     pG->B1i[k-1][j][i])*dx1i;
        		db2 = (    pG->B2i[k-1][j+1][i  ] -     pG->B2i[k-1][j][i])*dx2i;
        		db3 = (    pG->B3i[k  ][j  ][i  ] -     pG->B3i[k-1][j][i])*dx3i;

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

				MHD_src_By = (pG->U[k-1][j][i].M1/pG->U[k-1][j][i].d)*l1;
				MHD_src_Bz = (pG->U[k-1][j][i].M2/pG->U[k-1][j][i].d)*l2;

        			Wl[k].By += hdt*MHD_src_By;
        			Wl[k].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone k */
        			db1 = (    pG->B1i[k][j][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k][j+1][i] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j][i] -     pG->B3i[k][j][i])*dx3i;

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

        		MHD_src_By = (pG->U[k][j][i].M1/pG->U[k][j][i].d)*l1;
       		 	MHD_src_Bz = (pG->U[k][j][i].M2/pG->U[k][j][i].d)*l2;

       		 	Wr[k].By += hdt*MHD_src_By;
        		Wr[k].Bz += hdt*MHD_src_Bz;
      }			
#endif


		
      		if (StaticGravPot != NULL){
        		for (k=kl+1; k<=ku; j++) {
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);

          			phicr = (*StaticGravPot)(x1,x2, x3             );
          			phicl = (*StaticGravPot)(x1,x2,(x3-    pG->dx3));
          			phifc = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          			Wl[k].Vx -= dtodx3*(phifc - phicl);
          			Wr[k].Vx -= dtodx3*(phicr - phifc);
        	
			}
      		}

/*-----Calculate flux along x3 direction ---------*/

		for (k=kl+1; k<=ku; k++) {
        		Ul_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wl[k],&Bxi[k]);
        		Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[k],&Bxi[k]);

#ifdef RADIATION_MHD
        		Bx = B3_x3Face[k][j][i];
#endif
			/* Take the time step with x3Flux[][][].d */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[k],Wr[k],Bx,&x3Flux[k][j][i]);
      		}

		}/* END i */
	}/* end j */	




/*=== STEP 5:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 5a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef RADIATION_MHD
/* emf1 */
  		for (k=kl; k<=ku; k++) {
    			for (j=jl; j<=ju; j++) {
      				for (i=il; i<=iu; i++) {
        				emf1_cc[k][j][i] = (pG->U[k][j][i].B2c*pG->U[k][j][i].M3 -
			    		pG->U[k][j][i].B3c*pG->U[k][j][i].M2)
                              		/pG->U[k][j][i].d;
        				emf2_cc[k][j][i] = (pG->U[k][j][i].B3c*pG->U[k][j][i].M1 -
			    		pG->U[k][j][i].B1c*pG->U[k][j][i].M3)
                              		/pG->U[k][j][i].d;
        				emf3_cc[k][j][i] = (pG->U[k][j][i].B1c*pG->U[k][j][i].M2 -
			    		pG->U[k][j][i].B2c*pG->U[k][j][i].M1)
                              		/pG->U[k][j][i].d;
      				}
    			}
  		}
  		
		integrate_emf1_corner(pG);
  		integrate_emf2_corner(pG);
  		integrate_emf3_corner(pG);

/*--- Step 5b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  		for (k=kl+1; k<=ku-1; k++) {
    			for (j=jl+1; j<=ju-1; j++) {
      				for (i=il+1; i<=iu-1; i++) {

        			B1_x1Face[k][j][i] += hdtodx3 * (emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                              		hdtodx2 * (emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        			B2_x2Face[k][j][i] += hdtodx1 * (emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                              		hdtodx3 * (emf1[k+1][j  ][i  ] - emf1[k][j][i]);

        			B3_x3Face[k][j][i] += hdtodx2 * (    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                              		hdtodx1 * (emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      			}

      			B1_x1Face[k][j][iu] += hdtodx3 * (emf2[k+1][j  ][iu]-emf2[k][j][iu]) -
                               hdtodx2 * (emf3[k  ][j+1][iu]-emf3[k][j][iu]);
    			}
    			
			for (i=il+1; i<=iu-1; i++) {
      				B2_x2Face[k][ju][i] += hdtodx1 * (emf3[k  ][ju][i+1]-emf3[k][ju][i]) -
                               		hdtodx3 * (emf1[k+1][ju][i  ]-emf1[k][ju][i]);
    				}
  			}

  			for (j=jl+1; j<=ju-1; j++) {
    				for (i=il+1; i<=iu-1; i++) {

      					B3_x3Face[ku][j][i] += hdtodx2 * (    emf1[ku][j+1][i  ] -     emf1[ku][j][i]) -
                             			hdtodx1 * (emf2[ku][j  ][i+1] - emf2[ku][j][i]);
    				}
  			}
#endif /* MHD */

/*=== STEP 6: Correct x1-interface states with transverse flux gradients =====*/
/*--- Step 6a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes and x3 fluxes  computed in previous steps.
 * Since the fluxes come from an -sweep, (x,y,z) on RHS -> (z,x,y) for y and (x,y,z) on RHS ->(z,x,y) on LHS
 * The order is rotated-------------------------------------------
 */

  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

        		Ul_x1Face[k][j][i].d -=hdtodx2*(x2Flux[k][j+1][i-1].d -x2Flux[k][j][i-1].d );
        		Ul_x1Face[k][j][i].Mx-=hdtodx2*(x2Flux[k][j+1][i-1].Mz-x2Flux[k][j][i-1].Mz);
        		Ul_x1Face[k][j][i].My-=hdtodx2*(x2Flux[k][j+1][i-1].Mx-x2Flux[k][j][i-1].Mx);
        		Ul_x1Face[k][j][i].Mz-=hdtodx2*(x2Flux[k][j+1][i-1].My-x2Flux[k][j][i-1].My);

        		Ul_x1Face[k][j][i].E -=hdtodx2*(x2Flux[k][j+1][i-1].E -x2Flux[k][j][i-1].E );

#ifdef RADIATION_MHD
/* Update B3 */
			Ul_x1Face[k][j][i].Bz+=hdtodx2*0.5*((emf1[k  ][j+1][i-1] - emf1[k  ][j][i-1]) 
					+ (emf1[k+1][j+1][i-1] - emf1[k+1][j][i-1]));
#endif

       			Ur_x1Face[k][j][i].d -=hdtodx2*(x2Flux[k][j+1][i  ].d -x2Flux[k][j][i  ].d );
        		Ur_x1Face[k][j][i].Mx-=hdtodx2*(x2Flux[k][j+1][i  ].Mz-x2Flux[k][j][i  ].Mz);
        		Ur_x1Face[k][j][i].My-=hdtodx2*(x2Flux[k][j+1][i  ].Mx-x2Flux[k][j][i  ].Mx);
        		Ur_x1Face[k][j][i].Mz-=hdtodx2*(x2Flux[k][j+1][i  ].My-x2Flux[k][j][i  ].My);
        		Ur_x1Face[k][j][i].E -=hdtodx2*(x2Flux[k][j+1][i  ].E -x2Flux[k][j][i  ].E );
#ifdef RADIATION_MHD
/* Update B3 */
			Ur_x1Face[k][j][i].Bz+=hdtodx2* 0.5 * ((emf1[k  ][j+1][i] - emf1[k  ][j][i]) + (emf1[k+1][j+1][i] - emf1[k+1][j][i]));
#endif


/*--- Step 6b ------------------------------------------------------------------
 * Correct x1-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */

        		Ul_x1Face[k][j][i].d  -= hdtodx3 * (x3Flux[k+1][j][i-1].d -x3Flux[k][j][i-1].d );
        		Ul_x1Face[k][j][i].Mx -= hdtodx3 * (x3Flux[k+1][j][i-1].My-x3Flux[k][j][i-1].My);
        		Ul_x1Face[k][j][i].My -= hdtodx3 * (x3Flux[k+1][j][i-1].Mz-x3Flux[k][j][i-1].Mz);
        		Ul_x1Face[k][j][i].Mz -= hdtodx3 * (x3Flux[k+1][j][i-1].Mx-x3Flux[k][j][i-1].Mx);

        		Ul_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k+1][j][i-1].E -x3Flux[k][j][i-1].E );
#ifdef RADIATION_MHD
/* Update B2 */
			Ul_x1Face[k][j][i].By -= hdtodx3 * 0.5* ((emf1[k+1][j  ][i-1] - emf1[k][j  ][i-1]) +
	   					(emf1[k+1][j+1][i-1] - emf1[k][j+1][i-1]));
#endif

        		Ur_x1Face[k][j][i].d  -= hdtodx3 * (x3Flux[k+1][j][i  ].d -x3Flux[k][j][i  ].d );
        		Ur_x1Face[k][j][i].Mx -= hdtodx3 * (x3Flux[k+1][j][i  ].My-x3Flux[k][j][i  ].My);
        		Ur_x1Face[k][j][i].My -= hdtodx3 * (x3Flux[k+1][j][i  ].Mz-x3Flux[k][j][i  ].Mz);
        		Ur_x1Face[k][j][i].Mz -= hdtodx3 * (x3Flux[k+1][j][i  ].Mx-x3Flux[k][j][i  ].Mx);
        		Ur_x1Face[k][j][i].E  -= hdtodx3 * (x3Flux[k+1][j][i  ].E -x3Flux[k][j][i  ].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ur_x1Face[k][j][i].By -= hdtodx3 * 0.5*((emf1[k+1][j  ][i] - emf1[k][j  ][i]) +
	   					(emf1[k+1][j+1][i] - emf1[k][j+1][i]));
#endif

      }
    }
  }


/*--- Step 6b ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2- and x3-flux-gradients to the
 * conservative variables on the x1Face.  Limiting is used as in GS (2007)
 */

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

        			db1 = (    pG->B1i[k  ][j  ][i  ] -     pG->B1i[k][j][i-1])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i-1] -     pG->B2i[k][j][i-1])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i-1] -     pG->B3i[k][j][i-1])*dx3i;
        			B1 = pG->U[k][j][i-1].B1c;
        			B2 = pG->U[k][j][i-1].B2c;
        			B3 = pG->U[k][j][i-1].B3c;
        			V2 = pG->U[k][j][i-1].M2/pG->U[k][j][i-1].d;
        			V3 = pG->U[k][j][i-1].M3/pG->U[k][j][i-1].d;

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

        			Ul_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
        			V2 = pG->U[k][j][i].M2/pG->U[k][j][i].d;
        			V3 = pG->U[k][j][i].M3/pG->U[k][j][i].d;

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

        			Ur_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );

      }
    }
  }
#endif /* RADIATION MHD */


/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x2-Flux
 * and x3-Flux gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {
        			cc_pos(pG,i,j,k,&x1,&x2,&x3);
        			
				phic = (*StaticGravPot)(x1, x2             ,x3);
        			phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        			phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

/* correct right states; x2 and x3 gradients */

        			Ur_x1Face[k][j][i].My -= hdtodx2 * (phir-phil)*pG->U[k][j][i].d;

        			Ur_x1Face[k][j][i].E -= hdtodx2 * (x2Flux[k][j  ][i  ].d*(phic - phil)
                                  + x2Flux[k][j+1][i  ].d*(phir - phic));


        			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        			phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));
        
        			Ur_x1Face[k][j][i].Mz -= hdtodx3 * (phir-phil)*pG->U[k][j][i].d;

        			Ur_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j][i  ].d*(phic - phil)
                                  	+ x3Flux[k+1][j][i  ].d*(phir - phic));


/* correct left states; x2 and x3 gradients */
        			phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        			phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        			phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);


        			Ul_x1Face[k][j][i].My -= hdtodx2 * (phir-phil)*pG->U[k][j][i-1].d;

        			Ul_x1Face[k][j][i].E -= hdtodx2 * (x2Flux[k][j  ][i-1].d*(phic - phil)
                                  	+ x2Flux[k][j+1][i-1].d*(phir - phic));


        			phir = (*StaticGravPot)((x1-pG->dx1),x2,(x3+0.5*pG->dx3));
        			phil = (*StaticGravPot)((x1-pG->dx1),x2,(x3-0.5*pG->dx3));
        
        			Ul_x1Face[k][j][i].Mz -= hdtodx3 * (phir-phil)*pG->U[k][j][i-1].d;

        			Ul_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j][i-1].d*(phic - phil)
                                  		+ x3Flux[k+1][j][i-1].d*(phir - phic));

      				}
    			}
  		}
	}


/*===========================================================================*/
/*=== STEP 7A: Correct x2-interface states with transverse flux gradients =====*/

	
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju; j++) {
      			for (i=il+1; i<=iu-1; i++) {
 
       				Ul_x2Face[k][j][i].d -=hdtodx1 * (x1Flux[k][j-1][i+1].d  - x1Flux[k][j-1][i].d );
        			Ul_x2Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k][j-1][i+1].My - x1Flux[k][j-1][i].My);
        			Ul_x2Face[k][j][i].My-=hdtodx1 * (x1Flux[k][j-1][i+1].Mz - x1Flux[k][j-1][i].Mz);
        			Ul_x2Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k][j-1][i+1].Mx - x1Flux[k][j-1][i].Mx);
       				Ul_x2Face[k][j][i].E -=hdtodx1 * (x1Flux[k][j-1][i+1].E - x1Flux[k][j-1][i].E );

#ifdef RADIATION_MHD
/* Update B3 */
				Ul_x2Face[k][j][i].By -=hdtodx1 * 0.5*
	  				((emf2[k  ][j-1][i+1] - emf2[k  ][j-1][i]) + 
	   				(emf2[k+1][j-1][i+1] - emf2[k+1][j-1][i]));
#endif

        			Ur_x2Face[k][j][i].d -=hdtodx1 * (x1Flux[k][j  ][i+1].d  - x1Flux[k][j  ][i].d );
        			Ur_x2Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k][j  ][i+1].My - x1Flux[k][j  ][i].My);
        			Ur_x2Face[k][j][i].My-=hdtodx1 * (x1Flux[k][j  ][i+1].Mz - x1Flux[k][j  ][i].Mz);
        			Ur_x2Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k][j  ][i+1].Mx - x1Flux[k][j  ][i].Mx);

        			Ur_x2Face[k][j][i].E -=hdtodx1 * (x1Flux[k][j  ][i+1].E  - x1Flux[k][j  ][i].E );
#ifdef RADIATION_MHD
/* Update B3 */
				Ur_x2Face[k][j][i].By-=hdtodx1 * 0.5*
	  				((emf2[k  ][j][i+1] - emf2[k  ][j][i]) + 
	   				(emf2[k+1][j][i+1] - emf2[k+1][j][i]));
#endif


/*--- Step 7b ------------------------------------------------------------------
 * Correct x2-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

        			Ul_x2Face[k][j][i].d -=hdtodx3 * (x3Flux[k+1][j-1][i].d -x3Flux[k][j-1][i].d );
        			Ul_x2Face[k][j][i].Mx-=hdtodx3 * (x3Flux[k+1][j-1][i].Mz-x3Flux[k][j-1][i].Mz);
        			Ul_x2Face[k][j][i].My-=hdtodx3 * (x3Flux[k+1][j-1][i].Mx-x3Flux[k][j-1][i].Mx);
        			Ul_x2Face[k][j][i].Mz-=hdtodx3 * (x3Flux[k+1][j-1][i].My-x3Flux[k][j-1][i].My);

        			Ul_x2Face[k][j][i].E -=hdtodx3 * (x3Flux[k+1][j-1][i].E -x3Flux[k][j-1][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
				Ul_x2Face[k][j][i].Bz+=hdtodx3 * 0.5*
	  				((emf2[k+1][j-1][i  ] - emf2[k][j-1][i  ]) +
	   				(emf2[k+1][j-1][i+1] - emf2[k][j-1][i+1]));
#endif

        			Ur_x2Face[k][j][i].d -=hdtodx3 * (x3Flux[k+1][j  ][i].d -x3Flux[k][j  ][i].d );
        			Ur_x2Face[k][j][i].Mx-=hdtodx3 * (x3Flux[k+1][j  ][i].Mz-x3Flux[k][j  ][i].Mz);
        			Ur_x2Face[k][j][i].My-=hdtodx3 * (x3Flux[k+1][j  ][i].Mx-x3Flux[k][j  ][i].Mx);
        			Ur_x2Face[k][j][i].Mz-=hdtodx3 * (x3Flux[k+1][j  ][i].My-x3Flux[k][j  ][i].My);

        			Ur_x2Face[k][j][i].E -=hdtodx3 * (x3Flux[k+1][j  ][i].E -x3Flux[k][j  ][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
				Ur_x2Face[k][j][i].Bz+=hdtodx3 * 0.5*
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

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju; j++) {
      			for (i=il+1; i<=iu-1; i++) {

        			db1 = (	   pG->B1i[k  ][j-1][i+1] -     pG->B1i[k][j-1][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j  ][i  ] -     pG->B2i[k][j-1][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j-1][i  ] -     pG->B3i[k][j-1][i])*dx3i;
        			B1 = pG->U[k][j-1][i].B1c;
        			B2 = pG->U[k][j-1][i].B2c;
        			B3 = pG->U[k][j-1][i].B3c;
        			V1 = pG->U[k][j-1][i].M1/pG->U[k][j-1][i].d;
        			V3 = pG->U[k][j-1][i].M3/pG->U[k][j-1][i].d;

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

        			Ul_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
        			V1 = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        			V3 = pG->U[k][j][i].M3/pG->U[k][j][i].d;

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
        			Ur_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );

      				}
    			}
  		}
#endif /* RADIATION MHD */

				
/*--- Step 7d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x3-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  	if (StaticGravPot != NULL){
  		for (k=kl+1; k<=ku-1; k++) {
    			for (j=jl+1; j<=ju; j++) {
      				for (i=il+1; i<=iu-1; i++) {
        				cc_pos(pG,i,j,k,&x1,&x2,&x3);
        
			phic = (*StaticGravPot)((x1            ),x2,x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* correct right states; x1 and x3 gradients */

        		Ur_x2Face[k][j][i].Mz -= hdtodx1 * (phir-phil)*pG->U[k][j][i].d;

        		Ur_x2Face[k][j][i].E -= hdtodx1 * (x1Flux[k][j  ][i  ].d*(phic - phil)
                                  + x1Flux[k][j  ][i+1].d*(phir - phic));

        		phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        		phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

        		Ur_x2Face[k][j][i].My -= hdtodx3 * (phir-phil)*pG->U[k][j][i].d;

        		Ur_x2Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j  ][i].d*(phic - phil)
                                  	+ x3Flux[k+1][j  ][i].d*(phir - phic));


/* correct left states; x1 and x3 gradients */
        		phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

        		Ul_x2Face[k][j][i].Mz -= hdtodx1 * (phir-phil)*pG->U[k][j-1][i].d;
        		Ul_x2Face[k][j][i].E -= hdtodx1 * (x1Flux[k][j-1][i  ].d*(phic - phil)
                                  + x1Flux[k][j-1][i+1].d*(phir - phic));

        		phir = (*StaticGravPot)(x1,(x2-pG->dx2),(x3+0.5*pG->dx3));
        		phil = (*StaticGravPot)(x1,(x2-pG->dx2),(x3-0.5*pG->dx3));

        		Ul_x2Face[k][j][i].My -= hdtodx3 * (phir-phil)*pG->U[k][j-1][i].d;

         		Ul_x2Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j-1][i].d*(phic - phil)
                                  + x3Flux[k+1][j-1][i].d*(phir - phic));

      				}
    			}
  		}
	}

/*========================Shearing box source term gradient *===============*/	
	
/* Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
* x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
*    Vx source term is (dt/2)( 2 Omega_0 Vy)
*    Vy source term is (dt/2)(-2 Omega_0 Vx)
*    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
*/
		
#ifdef SHEARING_BOX
		if (ShearingBoxPot != NULL){
			for (k=kl+1; k<=ku-1; k++) {
				for (j=jl+1; j<=ju; j++) {
					for (i=il+1; i<=iu-1; i++) {
						cc_pos(pG,i,j,k,&x1,&x2,&x3);
						phic = (*ShearingBoxPot)((x1            ),x2,x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
						
						/* correct right states; x1 and x3 gradients */
						Ur_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
						Ur_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j  ][i  ].d*(phic - phil)
													+ x1Flux[k][j  ][i+1].d*(phir - phic));
#endif
						
						phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
						
						Ur_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
						Ur_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j  ][i].d*(phic - phil)
													+ x3Flux[k+1][j  ][i].d*(phir - phic));
#endif
						
						/* correct left states; x1 and x3 gradients */
						phic = (*ShearingBoxPot)((x1            ),(x2-pG->dx2),x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);
						
						Ul_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j-1][i].d;
#ifndef BAROTROPIC
						Ul_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j-1][i  ].d*(phic - phil)
													+ x1Flux[k][j-1][i+1].d*(phir - phic));
#endif
						phir = (*ShearingBoxPot)(x1,(x2-pG->dx2),(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,(x2-pG->dx2),(x3-0.5*pG->dx3));
						
						Ul_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j-1][i].d;
#ifndef BAROTROPIC
						Ul_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j-1][i].d*(phic - phil)
													+ x3Flux[k+1][j-1][i].d*(phir - phic));
#endif
					}
				}
			}
		}
	
	for (k=kl+1; k<=ku-1; k++) {
		for (j=jl+1; j<=ju; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x2Face[k][j][i].Mz += pG->dt*Omega_0*pG->U[k][j][i].M2;
				Ul_x2Face[k][j][i].Mz += pG->dt*Omega_0*pG->U[k][j-1][i].M2;
#ifdef FARGO
				Ur_x2Face[k][j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
				Ul_x2Face[k][j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[k][j-1][i].M1;
#else
				Ur_x2Face[k][j][i].Mx -= pG->dt*Omega_0*pG->U[k][j][i].M1;
				Ul_x2Face[k][j][i].Mx -= pG->dt*Omega_0*pG->U[k][j-1][i].M1;
#endif
			}
		}
	}

#endif /* SHEARING_BOX */
	
/*========================Shearing box source term gradient *===============*/	
	
/*=== STEP 8: Correct x3-interface states with transverse flux gradients =====*/

/*--- Step 8a ------------------------------------------------------------------
 * Correct x3-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

  	for (k=kl+1; k<=ku; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {

        		Ul_x3Face[k][j][i].d -=hdtodx1 * (x1Flux[k-1][j][i+1].d -x1Flux[k-1][j][i].d );
        		Ul_x3Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k-1][j][i+1].Mz-x1Flux[k-1][j][i].Mz);
        		Ul_x3Face[k][j][i].My-=hdtodx1 * (x1Flux[k-1][j][i+1].Mx-x1Flux[k-1][j][i].Mx);
        		Ul_x3Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k-1][j][i+1].My-x1Flux[k-1][j][i].My);

        		Ul_x3Face[k][j][i].E -=hdtodx1 * (x1Flux[k-1][j][i+1].E -x1Flux[k-1][j][i].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ul_x3Face[k][j][i].Bz+=hdtodx1 * 0.5*
	  			((emf3[k-1][j  ][i+1] - emf3[k-1][j  ][i]) +
	   			(emf3[k-1][j+1][i+1] - emf3[k-1][j+1][i]));
#endif

        		Ur_x3Face[k][j][i].d -=hdtodx1 * (x1Flux[k  ][j][i+1].d  - x1Flux[k  ][j][i].d );
        		Ur_x3Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k  ][j][i+1].Mz - x1Flux[k  ][j][i].Mz);
        		Ur_x3Face[k][j][i].My-=hdtodx1 * (x1Flux[k  ][j][i+1].Mx - x1Flux[k  ][j][i].Mx);
        		Ur_x3Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k  ][j][i+1].My - x1Flux[k  ][j][i].My);

        		Ur_x3Face[k][j][i].E -=hdtodx1 * (x1Flux[k  ][j][i+1].E - x1Flux[k  ][j][i].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ur_x3Face[k][j][i].Bz+=hdtodx1 * 0.5*
	  			((emf3[k][j  ][i+1] - emf3[k][j  ][i]) +
	   			(emf3[k][j+1][i+1] - emf3[k][j+1][i]));
#endif


/*--- Step 8b ------------------------------------------------------------------
 * Correct x3-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (y,z,x) on LHS 
 */

        		Ul_x3Face[k][j][i].d -=hdtodx2 * (x2Flux[k-1][j+1][i].d -x2Flux[k-1][j][i].d );
        		Ul_x3Face[k][j][i].Mx-=hdtodx2 * (x2Flux[k-1][j+1][i].My-x2Flux[k-1][j][i].My);
        		Ul_x3Face[k][j][i].My-=hdtodx2 * (x2Flux[k-1][j+1][i].Mz-x2Flux[k-1][j][i].Mz);
        		Ul_x3Face[k][j][i].Mz-=hdtodx2 * (x2Flux[k-1][j+1][i].Mx-x2Flux[k-1][j][i].Mx);

        		Ul_x3Face[k][j][i].E -=hdtodx2 * (x2Flux[k-1][j+1][i].E -x2Flux[k-1][j][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
			Ul_x3Face[k][j][i].By-=hdtodx2 * 0.5*
	  			((emf3[k-1][j+1][i  ] - emf3[k-1][j][i  ]) +
	   			(emf3[k-1][j+1][i+1] - emf3[k-1][j][i+1]));
#endif

        		Ur_x3Face[k][j][i].d -=hdtodx2 * (x2Flux[k  ][j+1][i].d -x2Flux[k  ][j][i].d );
        		Ur_x3Face[k][j][i].Mx-=hdtodx2 * (x2Flux[k  ][j+1][i].My-x2Flux[k  ][j][i].My);
        		Ur_x3Face[k][j][i].My-=hdtodx2 * (x2Flux[k  ][j+1][i].Mz-x2Flux[k  ][j][i].Mz);
        		Ur_x3Face[k][j][i].Mz-=hdtodx2 * (x2Flux[k  ][j+1][i].Mx-x2Flux[k  ][j][i].Mx);

        		Ur_x3Face[k][j][i].E -=hdtodx2 * (x2Flux[k  ][j+1][i].E -x2Flux[k  ][j][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
			Ur_x3Face[k][j][i].By-=hdtodx2 * 0.5*
	  			((emf3[k][j+1][i  ] - emf3[k][j][i  ]) +
	   			(emf3[k][j+1][i+1] - emf3[k][j][i+1]));
#endif

      			}
    		}
  	}

		
/*--- Step 8c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x2-flux-gradients to the
 * conservative variables on the x3Face.  Limiting is used as in GS07.
 */

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {


        			db1 = (    pG->B1i[k-1][j  ][i+1] -     pG->B1i[k-1][j][i])*dx1i;
        			db2 = (    pG->B2i[k-1][j+1][i  ] -     pG->B2i[k-1][j][i])*dx2i;
        			db3 = (    pG->B3i[k  ][j  ][i  ] -     pG->B3i[k-1][j][i])*dx3i;
        			B1 = pG->U[k-1][j][i].B1c;
        			B2 = pG->U[k-1][j][i].B2c;
        			B3 = pG->U[k-1][j][i].B3c;
				V1 = pG->U[k-1][j][i].M1/pG->U[k-1][j][i].d;
				V2 = pG->U[k-1][j][i].M2/pG->U[k-1][j][i].d;

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

				Ul_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
				V1 = pG->U[k][j][i].M1/pG->U[k][j][i].d;
				V2 = pG->U[k][j][i].M2/pG->U[k][j][i].d;

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

				Ur_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );

      				}
    			}
  		}
#endif /* MHD */


/*--- Step 8d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x2-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  	if (StaticGravPot != NULL){
  		for (k=kl+1; k<=ku; k++) {
    			for (j=jl+1; j<=ju-1; j++) {
      				for (i=il+1; i<=iu-1; i++) {
        				cc_pos(pG,i,j,k,&x1,&x2,&x3);
        				phic = (*StaticGravPot)((x1            ),x2,x3);
        				phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        				phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* correct right states; x1 and x2 gradients */

        				Ur_x3Face[k][j][i].My -= hdtodx1 * (phir-phil)*pG->U[k][j][i].d;

        				Ur_x3Face[k][j][i].E -= hdtodx1 * (x1Flux[k  ][j][i  ].d*(phic - phil)
                                  			+ x1Flux[k  ][j][i+1].d*(phir - phic));


        				phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        				phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        				Ur_x3Face[k][j][i].Mz -= hdtodx2 * (phir-phil)*pG->U[k][j][i].d;

        				Ur_x3Face[k][j][i].E -= hdtodx2 * (x2Flux[k  ][j  ][i].d*(phic - phil)
                                  		+ x2Flux[k  ][j+1][i].d*(phir - phic));


/* correct left states; x1 and x2 gradients */
        				phic = (*StaticGravPot)((x1            ),x2,(x3-pG->dx3));
        				phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,(x3-pG->dx3));
        				phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,(x3-pG->dx3));


        				Ul_x3Face[k][j][i].My -= hdtodx1 * (phir-phil)*pG->U[k-1][j][i].d;

        				Ul_x3Face[k][j][i].E -= hdtodx1 * (x1Flux[k-1][j][i  ].d*(phic - phil)
                                  		+ x1Flux[k-1][j][i+1].d*(phir - phic));


        				phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),(x3-pG->dx3));
        				phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),(x3-pG->dx3));

        				Ul_x3Face[k][j][i].Mz -= hdtodx2 * (phir-phil)*pG->U[k-1][j][i].d;

        				Ul_x3Face[k][j][i].E -= hdtodx2 * (x2Flux[k-1][j  ][i].d*(phic - phil)
                                  		+ x2Flux[k-1][j+1][i].d*(phir - phic));

      						}
    					}
  				}
			}
	
/*=============== Shearing box source terms ================*/
	
	/*--- Step 7d (cont) -----------------------------------------------------------
	 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
	 * x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
	 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
	 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
	 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
	 */
	
#ifdef SHEARING_BOX
	if (ShearingBoxPot != NULL){
		for (k=kl+1; k<=ku; k++) {
			for (j=jl+1; j<=ju-1; j++) {
				for (i=il+1; i<=iu-1; i++) {
					cc_pos(pG,i,j,k,&x1,&x2,&x3);
					phic = (*ShearingBoxPot)((x1            ),x2,x3);
					phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
					phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
					
					/* correct right states; x1 and x2 gradients */
					Ur_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
					Ur_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k  ][j][i  ].d*(phic - phil)
												+ x1Flux[k  ][j][i+1].d*(phir - phic));
#endif
					
					phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
					phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
					
					Ur_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
					Ur_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k  ][j  ][i].d*(phic - phil)
												+ x2Flux[k  ][j+1][i].d*(phir - phic));
#endif
					
					/* correct left states; x1 and x2 gradients */
					phic = (*ShearingBoxPot)((x1            ),x2,(x3-pG->dx3));
					phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,(x3-pG->dx3));
					phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,(x3-pG->dx3));
					
					Ul_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k-1][j][i].d;
#ifndef BAROTROPIC
					Ul_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k-1][j][i  ].d*(phic - phil)
												+ x1Flux[k-1][j][i+1].d*(phir - phic));
#endif
					
					phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),(x3-pG->dx3));
					phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),(x3-pG->dx3));
					
					Ul_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k-1][j][i].d;
#ifndef BAROTROPIC
					Ul_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k-1][j  ][i].d*(phic - phil)
												+ x2Flux[k-1][j+1][i].d*(phir - phic));
#endif
				}
			}
		}
	}
	
	for (k=kl+1; k<=ku; k++) {
		for (j=jl+1; j<=ju-1; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x3Face[k][j][i].My += pG->dt*Omega_0*pG->U[k][j][i].M2;
				Ul_x3Face[k][j][i].My += pG->dt*Omega_0*pG->U[k-1][j][i].M2;
#ifdef FARGO
				Ur_x3Face[k][j][i].Mz += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
				Ul_x3Face[k][j][i].Mz += hdt*(qshear-2.)*Omega_0*pG->U[k-1][j][i].M1;
#else
				Ur_x3Face[k][j][i].Mz -= pG->dt*Omega_0*pG->U[k][j][i].M1;
				Ul_x3Face[k][j][i].Mz -= pG->dt*Omega_0*pG->U[k-1][j][i].M1;
#endif
			}
		}
	}
#endif /* SHEARING_BOX */	
	
	
	
	
	
	
	
	
/*===========================================================*/
	



/*=== STEP 9: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 9a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, cooling, or RADIATION_MHD)
 */

    	for (k=kl+1; k<=ku-1; k++) {
      		for (j=jl+1; j<=ju-1; j++) {
			for (i=il+1; i<=iu-1; i++) {

          			dhalf[k][j][i] = pG->U[k][j][i].d 
            				- hdtodx1 * (    x1Flux[k  ][j  ][i+1].d -     x1Flux[k][j][i].d)
            				- hdtodx2 * (    x2Flux[k  ][j+1][i  ].d -     x2Flux[k][j][i].d)
            				- hdtodx3 * (    x3Flux[k+1][j  ][i  ].d -     x3Flux[k][j][i].d);

			}
      		}
    	}
  

/*--- Step 9b ------------------------------------------------------------------
 * Calculate cell centered emf_cc^{n+1/2}
 */

#ifdef RADIATION_MHD
/* Update momentum for half time step. We also need to add radiation source term * 
 * using modified Godunov method */

  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {

		/* load 1D vector */
			Usource.d  = pG->U[k][j][i].d;
      			Usource.Mx = pG->U[k][j][i].M1;
      			Usource.My = pG->U[k][j][i].M2;
      			Usource.Mz = pG->U[k][j][i].M3;
      			Usource.E  = pG->U[k][j][i].E;
			Usource.Er  = pG->U[k][j][i].Er;
    			Usource.Fr1  = pG->U[k][j][i].Fr1;
    			Usource.Fr2  = pG->U[k][j][i].Fr2;
    			Usource.Fr3  = pG->U[k][j][i].Fr3;
			Usource.Edd_11  = pG->U[k][j][i].Edd_11;
			Usource.Edd_21  = pG->U[k][j][i].Edd_21;
			Usource.Edd_22  = pG->U[k][j][i].Edd_22;
			Usource.Edd_31	= pG->U[k][j][i].Edd_31;
			Usource.Edd_32	= pG->U[k][j][i].Edd_32;
			Usource.Edd_33	= pG->U[k][j][i].Edd_33;
			Usource.Sigma_a  = pG->U[k][j][i].Sigma_a;
                        Usource.Sigma_t  = pG->U[k][j][i].Sigma_t;
      			Usource.By = pG->U[k][j][i].B2c;
      			Usource.Bz = pG->U[k][j][i].B3c;
			Bx = pG->U[k][j][i].B1c;
			
			density = Usource.d;
			velocity_x = Usource.Mx / density;
			velocity_y = Usource.My / density;
			velocity_z = Usource.Mz / density;

			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
				
			pressure = (pG->U[k][j][i].E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */

			pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
	
			temperature = pressure / (density * R_ideal);

			diffTEr = pow(temperature, 4.0) - pG->U[k][j][i].Er;
		

			Sigma_t    = pG->U[k][j][i].Sigma_t;
			Sigma_a	   = pG->U[k][j][i].Sigma_a;

			/* The Source term */
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz);

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = Sigma_t * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	

			SFmy = Sigma_t * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);

			SFmz = Sigma_t * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/


			/* co-moving flux */
			Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
			Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
			Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;

			/* Source term for momentum, not velocity*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);
			

			/* Now calculate flux gradient */

			divFlux1[1] = (x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx) / dx1;
			divFlux1[2] = (x1Flux[k][j][i+1].My - x1Flux[k][j][i].My) / dx1;
			divFlux1[3] = (x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz) / dx1;
			
			divFlux2[1] = (x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz) / dx2;
			divFlux2[2] = (x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx) / dx2;
			divFlux2[3] = (x2Flux[k][j+1][i].My - x2Flux[k][j][i].My) / dx2;
			
			divFlux3[1] = (x3Flux[k+1][j][i].My - x3Flux[k][j][i].My) / dx3;
			divFlux3[2] = (x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz) / dx3;
			divFlux3[3] = (x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx) / dx3;
		
			
/* Add source terms for fixed gravitational potential */
        		if (StaticGravPot != NULL){
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);


          			phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          			phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
          			divFlux1[1] += (phir-phil)*pG->U[k][j][i].d/dx1;


          			phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          			phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
          			divFlux2[2] += (phir-phil)*pG->U[k][j][i].d/dx2;

          			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          			phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));
          			divFlux3[3] += (phir-phil)*pG->U[k][j][i].d/dx3;
        		}


			/* Guess solution */
			Uguess[1] = pG->U[k][j][i].M1 + hdt * Source_Inv[1][1] * (Source[1] - divFlux1[1] - divFlux2[1] - divFlux3[1]);
			Uguess[2] = pG->U[k][j][i].M2 + hdt * Source_Inv[2][2] * (Source[2] - divFlux1[2] - divFlux2[2] - divFlux3[2]);
			Uguess[3] = pG->U[k][j][i].M3 + hdt * Source_Inv[3][3] * (Source[3] - divFlux1[3] - divFlux2[3] - divFlux3[3]);

			/* Guess the source temrm */
			velocity_x = Uguess[1] / density;
			velocity_y = Uguess[2] / density;
			velocity_z = Uguess[3] / density;
			
	
			Source_guess[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
			Source_guess[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);
			Source_guess[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);




			/* do the predict step */
			Errort[1] = pG->U[k][j][i].M1 + 0.5 * hdt * (Source[1] + Source_guess[1]) - hdt * (divFlux1[1] + divFlux2[1] + divFlux3[1]) - Uguess[1];
			Errort[2] = pG->U[k][j][i].M2 + 0.5 * hdt * (Source[2] + Source_guess[2]) - hdt * (divFlux1[2] + divFlux2[2] + divFlux3[2]) - Uguess[2];
			Errort[3] = pG->U[k][j][i].M3 + 0.5 * hdt * (Source[3] + Source_guess[3]) - hdt * (divFlux1[3] + divFlux2[3] + divFlux3[3]) - Uguess[3];



       		M1h = Uguess[1] + Source_Inv[1][1] * Errort[1];
			M2h = Uguess[2] + Source_Inv[2][2] * Errort[2];
			M3h = Uguess[3] + Source_Inv[3][3] * Errort[3];			

/*
        		Eh = pG->U[k][j][i].E
           			- q1*(rsf*x1Flux[k  ][j  ][i+1].E - lsf*x1Flux[k][j][i].E)
           			- q2*(    x2Flux[k  ][j+1][i  ].E -     x2Flux[k][j][i].E)
           			- q3*(    x3Flux[k+1][j  ][i  ].E -     x3Flux[k][j][i].E);

*/
					
					
					
/* Add shearing box source term */
#ifdef SHEARING_BOX
					if (ShearingBoxPot != NULL){
						cc_pos(pG,i,j,k,&x1,&x2,&x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
						M1h -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
						
						phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
						phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
						M2h -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;
						
						phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
						M3h -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;
					}
					
					M1h += pG->dt*Omega_0*pG->U[k][j][i].M2;
#ifdef FARGO
					M2h += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
#else
					M2h -= pG->dt*Omega_0*pG->U[k][j][i].M1;
#endif
#endif /* SHEARING_BOX */				



        		B1ch = 0.5*(    B1_x1Face[k][j][i] +     B1_x1Face[k  ][j  ][i+1]);
        		B2ch = 0.5*(    B2_x2Face[k][j][i] +     B2_x2Face[k  ][j+1][i  ]);
        		B3ch = 0.5*(    B3_x3Face[k][j][i] +     B3_x3Face[k+1][j  ][i  ]);
        		emf1_cc[k][j][i] = (B2ch*M3h - B3ch*M2h)/dhalf[k][j][i];
        		emf2_cc[k][j][i] = (B3ch*M1h - B1ch*M3h)/dhalf[k][j][i];
        		emf3_cc[k][j][i] = (B1ch*M2h - B2ch*M1h)/dhalf[k][j][i];

      			}
    		}
  	}

#endif

/*=== STEP 10: Compute 3D x1-Flux, x2-Flux, x3-Flux ===========================*/

/*--- Step 10a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

/*--- Step 10b ------------------------------------------------------------------
 * Compute 3D x1-fluxes from corrected L/R states.
 */

  	for (k=ks-1; k<=ke+1; k++) {
    		for (j=js-1; j<=je+1; j++) {
      			for (i=is; i<=ie+1; i++) {
#ifdef RADIATION_MHD
        		Bx = B1_x1Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x1Flux[k][j][i].d = dt;
			x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);
      			}
    		}
  	}

/*--- Step 10c ------------------------------------------------------------------
 * Compute 3D x2-fluxes from corrected L/R states.
 */

  	for (k=ks-1; k<=ke+1; k++) {
    		for (j=js; j<=je+1; j++) {
      			for (i=is-1; i<=ie+1; i++) {

#ifdef RADIATION_MHD
        		Bx = B2_x2Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[i],Wr[i],Bx, &x2Flux[k][j][i]);
      			}
    		}
  	}

/*--- Step 10d ------------------------------------------------------------------
 * Compute 3D x3-fluxes from corrected L/R states.
 */

  	for (k=ks; k<=ke+1; k++) {
    		for (j=js-1; j<=je+1; j++) {
      			for (i=is-1; i<=ie+1; i++) {

#ifdef RADIATION_MHD
        		Bx = B3_x3Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x3Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x3Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[i],Wr[i],Bx, &x3Flux[k][j][i]);
      			}
    		}
  	}


	
/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef RADIATION_MHD
  	integrate_emf1_corner(pG);
  	integrate_emf2_corner(pG);
  	integrate_emf3_corner(pG);
	
	
	
	/* Remap Ey at is and ie+1 to conserve Bz in shearing box */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
	
	/* compute remapped Ey from opposite side of grid */
	
    if (my_iproc == 0) {
		RemapEy_ix1(pD, emf2, remapEyiib);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
		RemapEy_ox1(pD, emf2, remapEyoib);
    }
	
	/* Now average Ey and remapped Ey */
	
    if (my_iproc == 0) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2[k][j][is]  = 0.5*(emf2[k][j][is] + remapEyiib[k][j]);
			}
		}
    }
	
    if (my_iproc == (pD->NGrid[0]-1)) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2[k][j][ie+1]  = 0.5*(emf2[k][j][ie+1] + remapEyoib[k][j]);
			}
		}
    }
#endif /* SHEARING_BOX */	
	

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  	for (k=ks; k<=ke; k++) {
    		for (j=js; j<=je; j++) {
      			for (i=is; i<=ie; i++) {

        		pG->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                            		dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        		pG->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                            		dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);

        		pG->B3i[k][j][i] += dtodx2*(    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                            		dtodx1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      			}

      		pG->B1i[k][j][ie+1] +=
        		dtodx3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
        		dtodx2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    		}

    		for (i=is; i<=ie; i++) {
      			pG->B2i[k][je+1][i] +=
        		dtodx1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
        		dtodx3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    		}
  	}

  	for (j=js; j<=je; j++) {
    		for (i=is; i<=ie; i++) {

      			pG->B3i[ke+1][j][i] += 
        			dtodx2*(    emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
        			dtodx1*(    emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    		}
  	}
#endif /* RADIATION_MHD */


/*-------Step 11: Predict and correct step after we get the flux------------- */
	
	/* Add shearinb box source term first */
	
/* Add gravitational (or shearing box) source terms as a Static Potential.
*   A Crank-Nicholson update is used for shearing box terms.
*   The energy source terms computed at cell faces are averaged to improve
* conservation of total energy.
*    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
*/
	
#ifdef SHEARING_BOX
	fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
	qom = qshear*Omega_0;
	for(k=ks; k<=ke; k++) {
		for(j=js; j<=je; j++) {
			for(i=is; i<=ie; i++) {
				cc_pos(pG,i,j,k,&x1,&x2,&x3);
				
				/* Store the current state */
				M1n  = pG->U[k][j][i].M1;
#ifdef FARGO
				dM2n = pG->U[k][j][i].M2;
#else
				dM2n = pG->U[k][j][i].M2 + qom*x1*pG->U[k][j][i].d;
#endif
				
				/* Calculate the flux for the y-momentum fluctuation */
				frx1_dM2 = x1Flux[k][j][i+1].My;
				flx1_dM2 = x1Flux[k][j][i  ].My;
				frx2_dM2 = x2Flux[k][j+1][i].Mx;
				flx2_dM2 = x2Flux[k][j  ][i].Mx;
				frx3_dM2 = x3Flux[k+1][j][i].Mz;
				flx3_dM2 = x3Flux[k  ][j][i].Mz;
#ifndef FARGO
				frx1_dM2 += qom*(x1+0.5*pG->dx1)*x1Flux[k][j][i+1].d;
				flx1_dM2 += qom*(x1-0.5*pG->dx1)*x1Flux[k][j][i  ].d;
				frx2_dM2 += qom*(x1            )*x2Flux[k][j+1][i].d;
				flx2_dM2 += qom*(x1            )*x2Flux[k][j  ][i].d;
				frx3_dM2 += qom*(x1            )*x3Flux[k+1][j][i].d;
				flx3_dM2 += qom*(x1            )*x3Flux[k  ][j][i].d;
#endif
				
				/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
				M1e = M1n - hdtodx1*(x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx)
				- hdtodx2*(x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz)
				- hdtodx3*(x3Flux[k+1][j][i].My - x3Flux[k][j][i].My);
				
				dM2e = dM2n - hdtodx1*(frx1_dM2 - flx1_dM2)
	            - hdtodx2*(frx2_dM2 - flx2_dM2) 
				- hdtodx3*(frx3_dM2 - flx3_dM2);

				
				/* Update the 1- and 2-momentum for the Coriolis and tidal
				 * potential momentum source terms using a Crank-Nicholson
				 * discretization for the momentum fluctuation equation. */
				
				pG->U[k][j][i].M1 += (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
				pG->U[k][j][i].M2 += 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;
				
#ifndef FARGO
				pG->U[k][j][i].M2 -= 0.5*qshear*om_dt*
				(x1Flux[k][j][i].d + x1Flux[k][j][i+1].d);
#endif
				
				/* Update the energy for a fixed potential.
				 * This update is identical to non-SHEARING_BOX below  */
				
				phic = (*ShearingBoxPot)((x1            ),x2,x3);
				phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
				phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
				pG->U[k][j][i].E -= dtodx1*(x1Flux[k][j][i  ].d*(phic - phil) +
											x1Flux[k][j][i+1].d*(phir - phic));
#endif
				
				phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
				phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
#ifndef BAROTROPIC
				pG->U[k][j][i].E -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil) +
											x2Flux[k][j+1][i].d*(phir - phic));
#endif
				
				phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
				phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
#ifndef BAROTROPIC
				pG->U[k][j][i].E -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil) +
											x3Flux[k+1][j][i].d*(phir - phic));
#endif
			}
		}
	}
	
#endif /* SHEARING_BOX */
	
	
		
/*----------Radiation quantities are not updated in the modified Godunov corrector step */
/*-----------NOTE that x1flux, x2flux, x3flux are flux from Remann problem. If there is extra-----* 
 *----------source terms, flux due to those source terms should also be added in the-------*
 *---------- modified Godunov method -------------------------------------*/

	for (k=ks; k<=ke; k++){
		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++) {

			/* Load 1D vector */
			Usource.d  = pG->U[k][j][i].d;
      			Usource.Mx = pG->U[k][j][i].M1;
      			Usource.My = pG->U[k][j][i].M2;
      			Usource.Mz = pG->U[k][j][i].M3;
      			Usource.E  = pG->U[k][j][i].E;
			Usource.Er  = pG->U[k][j][i].Er;
    			Usource.Fr1  = pG->U[k][j][i].Fr1;
    			Usource.Fr2  = pG->U[k][j][i].Fr2;
    			Usource.Fr3  = pG->U[k][j][i].Fr3;
			Usource.Edd_11  = pG->U[k][j][i].Edd_11;
			Usource.Edd_21  = pG->U[k][j][i].Edd_21;
			Usource.Edd_22  = pG->U[k][j][i].Edd_22;
			Usource.Edd_31	= pG->U[k][j][i].Edd_31;
			Usource.Edd_32	= pG->U[k][j][i].Edd_32;
			Usource.Edd_33	= pG->U[k][j][i].Edd_33;
			Usource.Sigma_a  = pG->U[k][j][i].Sigma_a;
                        Usource.Sigma_t  = pG->U[k][j][i].Sigma_t;

#ifdef RADIATION_MHD
      			Usource.By = pG->U[k][j][i].B2c;
      			Usource.Bz = pG->U[k][j][i].B3c;
			Bx = pG->U[k][j][i].B1c;
#else
			Bx = 0.0;
#endif /* MHD */
			
			density = Usource.d;
			velocity_x = Usource.Mx / density;
			velocity_y = Usource.My / density;
			velocity_z = Usource.Mz / density;

			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
				
			pressure = (pG->U[k][j][i].E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
			pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
			temperature = pressure / (density * R_ideal);

			diffTEr = pow(temperature, 4.0) - pG->U[k][j][i].Er;
		

			Sigma_t    = pG->U[k][j][i].Sigma_t;
			Sigma_a	   = pG->U[k][j][i].Sigma_a;

			/* The Source term */
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz);

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = Sigma_t * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	

			SFmy = Sigma_t * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);

			SFmz = Sigma_t * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/

			Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
			Source_Inv[4][1] = (-dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[1][1];
			Source_Inv[4][2] = (-dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[2][2];
			Source_Inv[4][3] = (-dt * Prat * Crat * SEmz/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[3][3];
			Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	


			/* co-moving flux */
			Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
			Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
			Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;

			/* Source term for momentum, not velocity*/
			Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
			Source[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);
			Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);
			
			/* Source term for energy */
			Source[4] = -Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat);

			/*--------Calculate the guess solution-------------*/
	
			divFlux1[0] = (x1Flux[k][j][i+1].d  - x1Flux[k][j][i].d ) / dx1;
			divFlux1[1] = (x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx) / dx1;
			divFlux1[2] = (x1Flux[k][j][i+1].My - x1Flux[k][j][i].My) / dx1;
			divFlux1[3] = (x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz) / dx1;
			divFlux1[4] = (x1Flux[k][j][i+1].E  - x1Flux[k][j][i].E ) / dx1; 

			divFlux2[0] = (x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d ) / dx2;
			divFlux2[1] = (x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz) / dx2;
			divFlux2[2] = (x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx) / dx2;
			divFlux2[3] = (x2Flux[k][j+1][i].My - x2Flux[k][j][i].My) / dx2;
			divFlux2[4] = (x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E ) / dx2; 

			divFlux3[0] = (x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d ) / dx3;
			divFlux3[1] = (x3Flux[k+1][j][i].My - x3Flux[k][j][i].My) / dx3;
			divFlux3[2] = (x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz) / dx3;
			divFlux3[3] = (x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx) / dx3;
			divFlux3[4] = (x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E ) / dx3; 

	/*------Flux due to static gravitational potential ---------*/

		if (StaticGravPot != NULL){
    			cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        		phic = (*StaticGravPot)((x1            ),x2,x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

			/* x direction */

			SourceFlux[1] = dhalf[k][j][i]*(phir-phil)/dx1;
        		SourceFlux[4]= (x1Flux[k][j][i  ].d*(phic - phil) +
				 	x1Flux[k][j][i+1].d*(phir - phic))/dx1;

			divFlux1[1] += SourceFlux[1];
			divFlux1[4] += SourceFlux[4];

			/* y direction */
			
        		phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        		phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        		SourceFlux[2] = dhalf[k][j][i]*(phir-phil)/dx2;
        		SourceFlux[4] = (x2Flux[k][j  ][i].d*(phic - phil) +
                                	     x2Flux[k][j+1][i].d*(phir - phic))/dx2;
			divFlux2[2] += SourceFlux[2];
			divFlux2[4] += SourceFlux[4];      		


			/* z direction */
			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          		phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

			SourceFlux[3] = dhalf[k][j][i]*(phir-phil)/dx3;
        		SourceFlux[4] = (x3Flux[k][j  ][i].d*(phic - phil) +
                                	     x3Flux[k + 1][j][i].d*(phir - phic))/dx3;
			divFlux3[3] += SourceFlux[3];
			divFlux3[4] += SourceFlux[4];   
	
  		}

	/*================== End static gravitational flux ================*/

	/* cacluate guess solution */

		for(n=0; n<5; n++) {
			tempguess[n] = 0.0;
			for(m=0; m<5; m++) {
				tempguess[n] += dt * Source_Inv[n][m] * (Source[m] - divFlux1[m] - divFlux2[m] - divFlux3[m]);
			}
		}

		Uguess[0] = pG->U[k][j][i].d  + tempguess[0];
		Uguess[1] = pG->U[k][j][i].M1 + tempguess[1];
		Uguess[2] = pG->U[k][j][i].M2 + tempguess[2];
		Uguess[3] = pG->U[k][j][i].M3 + tempguess[3];
		Uguess[4] = pG->U[k][j][i].E  + tempguess[4];

		/*  Uguess[0] = d; Uguess[1]=Mx; Uguess[2]=My; Uguess[3]=Mz, Uguess[4]=E */
	
		/* Now calculate the source term due to the guess solution */
		/* Flux is not changed during the correction step */


		
		density    = Uguess[0];
		velocity_x = Uguess[1] / density;
		velocity_y = Uguess[2] / density;
		velocity_z = Uguess[3] / density;

		velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);



		pressure = (Uguess[4] - 0.5 * (density * velocity * velocity)) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
		temperature = pressure / (density * R_ideal);

		diffTEr = pow(temperature, 4.0) - pG->U[k][j][i].Er;

		if(Opacity != NULL)
			Opacity(density,temperature, &Sigma_t, &Sigma_a, NULL);

		/* update source term */
		Usource.d  = Uguess[0];
		Usource.Mx = Uguess[1];
		Usource.My = Uguess[2];
		Usource.Mz = Uguess[3];
		Usource.E  = Uguess[4];
		Usource.Sigma_a = Sigma_a;
		Usource.Sigma_t = Sigma_t;

		/* The Source term */
		dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz);
	
		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = Sigma_t * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	

			SFmy = Sigma_t * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);

			SFmz = Sigma_t * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ Sigma_a * diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/

			Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
			Source_Inv[4][1] = (-dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[1][1];
			Source_Inv[4][2] = (-dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[2][2];
			Source_Inv[4][3] = (-dt * Prat * Crat * SEmz/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[3][3];
			Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	


			/* co-moving flux */
			Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
			Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
			Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;

			/* Source term for momentum */
			Source_guess[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
			Source_guess[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);
			Source_guess[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);
			
			/* Source term for total Energy */
			Source_guess[4] = -Prat * Crat * (Sigma_a * diffTEr + (Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat);


			/* Calculate the error term */
			Errort[0] = pG->U[k][j][i].d  + hdt * (Source[0] + Source_guess[0]) 
					- dt * (divFlux1[0] + divFlux2[0] + divFlux3[0]) - Uguess[0];
			Errort[1] = pG->U[k][j][i].M1 + hdt * (Source[1] + Source_guess[1]) 
					- dt * (divFlux1[1] + divFlux2[1] + divFlux3[1]) - Uguess[1];
			Errort[2] = pG->U[k][j][i].M2 + hdt * (Source[2] + Source_guess[2]) 
					- dt * (divFlux1[2] + divFlux2[2] + divFlux3[2]) - Uguess[2];
			Errort[3] = pG->U[k][j][i].M3 + hdt * (Source[3] + Source_guess[3]) 
					- dt * (divFlux1[3] + divFlux2[3] + divFlux3[3]) - Uguess[3];
			Errort[4] = pG->U[k][j][i].E  + hdt * (Source[4] + Source_guess[4]) 
					- dt * (divFlux1[4] + divFlux2[4] + divFlux3[4]) - Uguess[4];

			/* Calculate the correction */
			for(n=0; n<5; n++) {
				tempguess[n] = 0.0;
				for(m=0; m<5; m++) {
					tempguess[n] += Source_Inv[n][m] * Errort[m];
				}
			}


			/* Apply the correction */
		
			pG->U[k][j][i].d  = Uguess[0] + tempguess[0];
			pG->U[k][j][i].M1 = Uguess[1] + tempguess[1];
			pG->U[k][j][i].M2 = Uguess[2] + tempguess[2];
			pG->U[k][j][i].M3 = Uguess[3] + tempguess[3];
			pG->U[k][j][i].E  = Uguess[4] + tempguess[4];



			} /* End i */
		}/* End j */
	}/* End k */



/*--- Step 12 -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef RADIATION_MHD
  	for (k=ks; k<=ke; k++) {
    		for (j=js; j<=je; j++) {
      			for (i=is; i<=ie; i++) {

        			pG->U[k][j][i].B1c = 0.5*(    pG->B1i[k][j][i] +     pG->B1i[k][j][i+1]);
        			pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
        			pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
      			}
    		}
  	}
#endif /* RADIATION MHD */


		
	/* Boundary condition is applied in the main function */

	/* Update the opacity if Opacity function is set in the problem generator */
	if(Opacity != NULL){
		for (k=ks; k<=ke; k++){
			for (j=js; j<=je; j++) {
    				for (i=is; i<=ie; i++){
				
				density = pG->U[k][j][i].d;
				
				pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
				+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density ) * (Gamma - 1);
			/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
				temperature = pressure / (density * R_ideal);
			
			Opacity(density,temperature,&(pG->U[k][j][i].Sigma_t), &(pG->U[k][j][i].Sigma_a),NULL);

				}
			}
		}
	}


  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_init_3d: Allocate temporary integration arrays 
*/

void integrate_init_3d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

#ifdef RADIATION_MHD
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;

  if ((emf1_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* MHD */



  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef RADIATION_MHD
  if ((B1_x1Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif /* MHD */

  if ((U1d=(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul =(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur =(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((x1Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((dhalf = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;

	
#ifdef SHEARING_BOX
	if ((remapEyiib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
		goto on_error;
	if ((remapEyoib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
		goto on_error;
#endif /* End SHEARING_BOX */

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_3d:  Free temporary integration arrays 
 */

void integrate_destruct_3d(void)
{

#ifdef RADIATION_MHD
  if (emf1    != NULL) free_3d_array(emf1);
  if (emf2    != NULL) free_3d_array(emf2);
  if (emf3    != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
#endif /* MHD */



  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef RADIATION_MHD
  if (B1_x1Face != NULL) free_3d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array(B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array(B3_x3Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_3d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_3d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_3d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_3d_array(Ur_x2Face);
  if (Ul_x3Face != NULL) free_3d_array(Ul_x3Face);
  if (Ur_x3Face != NULL) free_3d_array(Ur_x3Face);
  if (x1Flux    != NULL) free_3d_array(x1Flux);
  if (x2Flux    != NULL) free_3d_array(x2Flux);
  if (x3Flux    != NULL) free_3d_array(x3Flux);
  if (dhalf     != NULL) free_3d_array(dhalf);
  if (phalf     != NULL) free_3d_array(phalf);

#ifdef SHEARING_BOX
	if (remapEyiib != NULL) free_2d_array(remapEyiib);
	if (remapEyoib != NULL) free_2d_array(remapEyoib);
#endif

  return;
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

#ifdef RADIATION_MHD
static void integrate_emf1_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
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

static void integrate_emf2_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
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
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

static void integrate_emf3_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de3_l1, de3_r1, de3_l2, de3_r2;
  Real rsf=1.0,lsf=1.0;

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
#ifdef CYLINDRICAL
        rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
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
#endif /* RADIATION MHD */

#endif /* CTU_INTEGRATOR */
