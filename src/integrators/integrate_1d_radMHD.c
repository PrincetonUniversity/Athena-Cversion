#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_1d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations using modified_Godunov method
 *   Updates U.[d,M1,M2,M3,E,B2c,B3c,s] in Grid structure, where U is of type
 *   ConsS. Adds gravitational source terms, self-gravity, and optically-thin
 *   cooling.
 *   To be added later.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_1d_radMHD()
 *   integrate_init_1d()
 *   integrate_destruct_1d()
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
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
/* x1Flux is Cons1DS type, but it is the flux */
static Cons1DS *Ul_x1Face=NULL, *Ur_x1Face=NULL, *x1Flux=NULL;

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL;

static Real *dhalf = NULL;


/* Variables needed for cylindrical coordinates */
#ifdef CYLINDRICAL
static Real *geom_src=NULL;
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_1d: 1D version of integrator of radiation MHD.
 * 1D direction is assumed to be along the direction 1   
 *
 */

void integrate_1d_radMHD(DomainS *pD)
{
  	GridS *pG=(pD->Grid);
	/*Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1; */
	Real dt=pG->dt, dx=pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;;
	Real dtodx1 = dt / dx;
	int il,iu, is = pG->is, ie = pG->ie;
  	int i, j, m, n;
	int js = pG->js;
	int ks = pG->ks;
	Cons1DS Uguess, divFlux;
	Real *pUguess, *pdivFlux, *pfluxl, *pfluxr, *pU1d;

	Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;
	Real aeff;


	Real SEE, SErho, SEm;
	Real alpha;
	
	
	Real temperature, velocity, pressure, Tguess;

	/* magnetic field is not included in guess and predict step.
	 * So these variables have fixed size */
	Real Source_Inv[5][5], tempguess[5], Source[5];
	Real Propa_44, Det;
	Real Source_guess[5], Errort[5];

	
	Real Sigma_t, Sigma_a;

	/* Initialize them to be zero */
	for(i=0; i<5; i++){
		Source[i] = 0.0;
		Source_guess[i] = 0.0;
		Errort[i] = 0.0;
		tempguess[i] = 0.0;
		for(j=0; j<5; j++) {
			Source_Inv[i][j] = 0.0;			
		if(i==j) {
		 Source_Inv[i][j] = 1.0;
		
		}
	}
	}


 	 il = is - 1;
  	 iu = ie + 1;

		
	
/* Temperatory variables used to calculate the Matrix  */
  	
  
  	

/* Load 1D vector of conserved variables */
 	for (i=is-nghost; i<=ie+nghost; i++) {
    		U1d[i].d  = pG->U_old[ks][js][i].d;
    		U1d[i].Mx = pG->U_old[ks][js][i].M1;
    		U1d[i].My = pG->U_old[ks][js][i].M2;
    		U1d[i].Mz = pG->U_old[ks][js][i].M3;
    		U1d[i].E  = pG->U_old[ks][js][i].E;
    		U1d[i].Er  = pG->U[ks][js][i].Er;
    		U1d[i].Fr1  = pG->U[ks][js][i].Fr1;
    		U1d[i].Fr2  = pG->U[ks][js][i].Fr2;
    		U1d[i].Fr3  = pG->U[ks][js][i].Fr3;
		U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;
		U1d[i].Edd_21  = pG->U[ks][js][i].Edd_21;
		U1d[i].Edd_22  = pG->U[ks][js][i].Edd_22;
		U1d[i].Edd_31  = pG->U[ks][js][i].Edd_31;
		U1d[i].Edd_32  = pG->U[ks][js][i].Edd_32;
		U1d[i].Edd_33  = pG->U[ks][js][i].Edd_33;
		U1d[i].Sigma_a  = pG->U[ks][js][i].Sigma_a;
		U1d[i].Sigma_t  = pG->U[ks][js][i].Sigma_t;
#ifdef RADIATION_MHD
		U1d[i].By = pG->U_old[ks][js][i].B2c;
		U1d[i].Bz = pG->U_old[ks][js][i].B3c;
		Bxc[i] = pG->U_old[ks][js][i].B1c;
		Bxi[i] = pG->B1i[ks][js][i];

#endif
	}
	
	


/*----Step 1: Backward Euler is done in the main function for the whole mesh--------*/	



/*-----Step 2a----------------
 *  Calcualte the left and right state */
	 for (i=is-nghost; i<=ie+nghost; i++) {
    		W[i] = Cons1D_to_Prim1D(&U1d[i], &Bxc[i]);
	  }

	/* for radiation hydro and MHD, the last dir represents optical thin (0) or thick (1) */
  	lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,0);

	
/*------Step 2b: Add source terms to the left and right state--------*/
	for(i=il+1; i<=ie+1; i++){

	/* For left state */
		pressure = W[i-1].P;
		temperature = pressure / (U1d[i-1].d * R_ideal);
		velocity = U1d[i-1].Mx / U1d[i-1].d;

		Sigma_t = pG->U[ks][js][i-1].Sigma_t;
		Sigma_a = pG->U[ks][js][i-1].Sigma_a;

		Tguess = pG->Tguess[ks][js][i-1];


	Source[1] = -Prat * (-Sigma_t * (U1d[i-1].Fr1/U1d[i-1].d 
	- (1.0 + U1d[i-1].Edd_11) * velocity * U1d[i-1].Er / (Crat * U1d[i-1].d))	
	+ Sigma_a * velocity * (pow(Tguess, 4.0) - U1d[i-1].Er)/(Crat*U1d[i-1].d));
	Source[4] = -Source[1] * U1d[i-1].d * velocity * (Gamma - 1.0)
	-(Gamma - 1.0) * Prat * Crat * (Sigma_a * (pow(Tguess,4.0) - U1d[i-1].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i-1].Fr1 - (1.0 + U1d[i-1].Edd_11) * velocity * U1d[i-1].Er / Crat)/Crat); 

		aeff = eff_sound(W[i-1],dt,0);

		alpha = ((aeff * aeff * W[i-1].d / W[i-1].P) - 1.0) / (Gamma - 1.0);

		Propa_44 = alpha;

		Wl[i].Vx += dt * Source[1] * 0.5;
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
		velocity = U1d[i].Mx / U1d[i].d;
		Tguess = pG->Tguess[ks][js][i];
		
		Sigma_t = U1d[i].Sigma_t;
		Sigma_a = U1d[i].Sigma_a;



	Source[1] = -Prat * (-Sigma_t * (U1d[i].Fr1/U1d[i].d 
	- (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / (Crat * U1d[i].d))	
	+ Sigma_a * velocity * (pow(Tguess, 4.0) - U1d[i].Er)/(Crat*U1d[i].d));
	Source[4] = -Source[1] * U1d[i].d * velocity * (Gamma - 1.0)
	-(Gamma - 1.0) * Prat * Crat * (Sigma_a * (pow(Tguess, 4.0) - U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i].Fr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)/Crat); 

		
		aeff = eff_sound(W[i],dt,0);

		alpha = ((aeff * aeff * W[i].d / W[i].P) - 1.0) / (Gamma - 1.0);


		Propa_44 = alpha;

		Wr[i].Vx += dt * Source[1] * 0.5;
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

/*-----------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
	
  	if (StaticGravPot != NULL){
    		for (i=il+1; i<=iu; i++) {
      		cc_pos(pG,i,js,ks,&x1,&x2,&x3);

		phicr = (*StaticGravPot)( x1             ,x2,x3);
      		phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
      		phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      		Wl[i].Vx -= dtodx1*(phifc - phicl);
      		Wr[i].Vx -= dtodx1*(phicr - phifc);

    		}
  	}


		
/*---------Step 2c--------------*/ 
/* Calculate the flux */
		 for (i=il+1; i<=iu; i++) {

  		  Ul_x1Face[i] = Prim1D_to_Cons1D(&Wl[i], &Bxi[i]);
   		  Ur_x1Face[i] = Prim1D_to_Cons1D(&Wr[i], &Bxi[i]);


		  x1Flux[i].d = dt;
		/* This is used to take dt to calculate alpha. 
                * x1Flux[i].d is recovered in the fluxes function by using Wl */

   		 fluxes(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i]);
		
 		 }

/*---
 * Calculate d^{n+1/2} (needed with static potential, or cooling)
 */
		
  		if ((StaticGravPot != NULL))
  		{
    			for (i=il+1; i<=iu-1; i++) {
      			dhalf[i] = U1d[i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
    			}
  		}
	

/*----Step 3------------------
 * Modified Godunov Corrector Scheme   */
/*----------Radiation quantities are not updated in the modified Godunov corrector step */
	for(i=il+1; i<=iu-1; i++) {
	if(!(pG->Flagtau[ks][js][i]))   
{

	pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1.0);
#ifdef  RADIATION_MHD 
	pressure -= (Gamma - 1.0) * 0.5 * (Bxc[i] * Bxc[i] + U1d[i].By * U1d[i].By + U1d[i].Bz * U1d[i].Bz);
#endif
	/* Should include magnetic energy for MHD */
	temperature = pressure / (U1d[i].d * R_ideal);

	/* Tguess is uesed for source term T^4 - Er */
	Tguess = pG->Tguess[ks][js][i];

	velocity = U1d[i].Mx / U1d[i].d;

	Sigma_t = U1d[i].Sigma_t;
	Sigma_a = U1d[i].Sigma_a;

	
	

/* The Source term */
	dSource(U1d[i], Bxc[i], &SEE, &SErho, &SEm, NULL, NULL);

	Det = 1.0 + dt * Prat * Crat * SEE + dt * Crat * Sigma_a; 

	
	Source_Inv[4][0] = -dt * Prat * Crat * SErho/Det;
	Source_Inv[4][1] = -dt * Prat * Crat * SEm/Det;
	Source_Inv[4][4] = (1.0 + dt * Crat * Sigma_a) / Det;

	/* For term T^4 - Er, we use guess temperature */
	
	Source[1] = -(-Sigma_t * (U1d[i].Fr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)
	+ Sigma_a * velocity * (pow(Tguess, 4.0) - U1d[i].Er)/Crat);
	Source[4] = - Crat * (Sigma_a * (pow(Tguess, 4.0) - U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i].Fr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)/Crat); 

	pdivFlux = (Real*)&(divFlux);
	pfluxr = (Real*)&(x1Flux[i+1]);
	pfluxl = (Real*)&(x1Flux[i]);

	for(m=0; m<5; m++)
		pdivFlux[m] = (pfluxr[m] - pfluxl[m]) / dx;
		
	/* ATTENTION: magnetic field is not included in this predict and correction step */
	for(n=0; n<5; n++) {
		tempguess[n] = 0.0;
		for(m=0; m<5; m++){
		tempguess[n] += dt * Source_Inv[n][m] * (Prat * Source[m] - pdivFlux[m]);
		}
	}

	tempguess[4] += -Source[4] * dt * dt * Prat * Crat * Sigma_a / Det 
			+ Source[1] * dt * dt * Prat * (2.0 * Sigma_a - Sigma_t) * velocity / Det;
	
	/* copy U1d to Uguess, some of the variables will be updated here */
	Uguess = U1d[i];

	pUguess = (Real*)&(Uguess);
	pU1d = (Real*)&(U1d[i]);	

	for(m=0; m<5; m++)
		pUguess[m]= pU1d[m] + tempguess[m];




	pressure = (Uguess.E - 0.5 * Uguess.Mx * Uguess.Mx / Uguess.d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
#ifdef  RADIATION_MHD 
	pressure -= (Gamma - 1.0) * 0.5 * (Bxc[i] * Bxc[i] + U1d[i].By * U1d[i].By + U1d[i].Bz * U1d[i].Bz);
#endif

	temperature = pressure / (Uguess.d * R_ideal);
	
	velocity = Uguess.Mx / Uguess.d;

	/* Use the updated opacity due to the change of temperature and density */
	if(Opacity != NULL)
		Opacity(Uguess.d, temperature, &Sigma_t, &Sigma_a, NULL);

	Uguess.Sigma_t = Sigma_t;
	Uguess.Sigma_a = Sigma_a;

	dSource(Uguess, Bxc[i], &SEE, &SErho, &SEm, NULL, NULL);
	Det = 1.0 + dt * Prat * Crat * SEE;

	/* Guess solution doesn't include Er and Fr1, Er and Fr1 are in U1d */
	/* NOTICE that we do not include correction from radiation variables, which are always 1st order accuracy */
	Source_guess[1] = - (-Sigma_t * (U1d[i].Fr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)
	+ Sigma_a * velocity * (pow(Tguess, 4.0) - U1d[i].Er)/Crat);
	Source_guess[4] = - Crat * (Sigma_a * (pow(Tguess, 4.0) - U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i].Fr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)/Crat); 

	for(m=0; m<5; m++)
		Errort[m] = pU1d[m] + 0.5 * dt * (Source_guess[m] + Source[m]) * Prat 
			- dt * pdivFlux[m] - pUguess[m];


	/* (I - dt Div_U S(U))^-1 for the guess solution */
	
	Source_Inv[4][0] = -dt * Prat * Crat * SErho/ Det;
	Source_Inv[4][1] = -dt * Prat * Crat * SEm/ Det;
	Source_Inv[4][4] = 1.0 / Det;


	for(m=0; m<5; m++){
		tempguess[m]=0.0;
		for(n=0; n<5; n++){
			tempguess[m] += Source_Inv[m][n] * Errort[n];
		}
		pU1d[m] = pUguess[m] + tempguess[m];

	}
	
	/* Update the quantity in the Grids */
	pUguess = (Real*)&(pG->U[ks][js][i]);

/* Only update the cells when flag is 0 */
   
	for(m=0; m<5; m++) pUguess[m] = pU1d[m];	

#ifdef RADIATION_MHD
/* magnetic field is independent of modified Godunov method */
	pG->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    	pG->U[ks][js][i].B3c -= dtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
/* For consistency, set B2i and B3i to cell-centered values.  */
    	pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    	pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;

#endif

}

	/*Boundary condition is applied in the main.c function*/

/*-----------Finish---------------------*/
	} /* End big loop i */	


 /* Add gravitational source terms as a Static Potential.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

	  if (StaticGravPot != NULL){
    		for (i=is; i<=ie; i++) {
		if(!(pG->Flagtau[ks][js][i])){
      		cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      		phic = (*StaticGravPot)((x1            ),x2,x3);
      		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

	        pG->U[ks][js][i].M1 -= dtodx1*dhalf[i]*(phir-phil);


      		pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
    		}
		}
  	}


	/* Now update the reaction coefficient Sigma_t and Sigma_a in the conserved variable 
	 * If function Opacity is set in the problem generator 
	 *
	 */ 

	if(Opacity != NULL){
		for(i=il+1; i<=iu-1; i++) {
			if(!(pG->Flagtau[ks][js][i])){
			pressure = (pG->U[ks][js][i].E - 0.5 * pG->U[ks][js][i].M1 * pG->U[ks][js][i].M1 / pG->U[ks][js][i].d )
				* (Gamma - 1);
/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[ks][js][i].B1c * pG->U[ks][js][i].B1c + pG->U[ks][js][i].B2c * pG->U[ks][js][i].B2c + pG->U[ks][js][i].B3c * pG->U[ks][js][i].B3c) * (Gamma - 1.0);
#endif
		
			temperature = pressure / (pG->U[ks][js][i].d * R_ideal);
	
			/* We do not need to calculate dSigma here */
			Opacity(pG->U[ks][js][i].d, temperature,&(pG->U[ks][js][i].Sigma_t), &(pG->U[ks][js][i].Sigma_a), NULL);
		}	
		}
	}


  return;

}




/* This will be done first before backward Euler */
void integrate_1d_radMHD_thick(DomainS *pD)
{
  	GridS *pG=(pD->Grid);
	/*Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1; */
	Real dt=pG->dt, dx=pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;;
	Real dtodx1 = dt / dx;
	int il,iu, is = pG->is, ie = pG->ie;
  	int i, m;
	int js = pG->js;
	int ks = pG->ks;
	Cons1DS Uguess, divFlux;
	Real *pdivFlux, *pfluxl, *pfluxr;

	Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic;



	Real SFFr, InvFr;
	Real alpha;
	
	/* For correction step */ 
	Real Error, correction;
	
	Real temperature, velocity, pressure;

	/* magnetic field is not included in guess and predict step.
	 * So these variables have fixed size */
	Real SourceFr, SourceFrguess, Esum;

		
	Real Sigma_t, Sigma_a;

	/* Initialize them to be zero */

 	 il = is - 1;
  	 iu = ie + 1;

		
	
/* Temperatory variables used to calculate the Matrix  */
  	
  
  	

/* Load 1D vector of conserved variables */
 	for (i=is-nghost; i<=ie+nghost; i++) {
    		U1d[i].d  = pG->U[ks][js][i].d;
    		U1d[i].Mx = pG->U[ks][js][i].M1;
    		U1d[i].My = pG->U[ks][js][i].M2;
    		U1d[i].Mz = pG->U[ks][js][i].M3;
    		U1d[i].E  = pG->U[ks][js][i].E;
    		U1d[i].Er  = pG->U[ks][js][i].Er;
    		U1d[i].Fr1  = pG->U[ks][js][i].Fr1;
    		U1d[i].Fr2  = pG->U[ks][js][i].Fr2;
    		U1d[i].Fr3  = pG->U[ks][js][i].Fr3;
		U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;
		U1d[i].Edd_21  = pG->U[ks][js][i].Edd_21;
		U1d[i].Edd_22  = pG->U[ks][js][i].Edd_22;
		U1d[i].Edd_31  = pG->U[ks][js][i].Edd_31;
		U1d[i].Edd_32  = pG->U[ks][js][i].Edd_32;
		U1d[i].Edd_33  = pG->U[ks][js][i].Edd_33;
		U1d[i].Sigma_a  = pG->U[ks][js][i].Sigma_a;
		U1d[i].Sigma_t  = pG->U[ks][js][i].Sigma_t;
#ifdef RADIATION_MHD
		U1d[i].By = pG->U[ks][js][i].B2c;
		U1d[i].Bz = pG->U[ks][js][i].B3c;
		Bxc[i] = pG->U[ks][js][i].B1c;
		Bxi[i] = pG->B1i[ks][js][i];

#endif
	}
	
	


/*----Step 1: Backward Euler is done in the main function for the whole mesh--------*/	



/*-----Step 2a----------------
 *  Calcualte the left and right state */
	 for (i=is-nghost; i<=ie+nghost; i++) {
    		W[i] = Cons1D_to_Prim1D(&U1d[i], &Bxc[i]);
	  }

	/* Left and right states for Er, Fr are also calculated here */
  	lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);

	
/*------Step 2b: Add source terms to the left and right state--------*/
	for(i=il+1; i<=ie+1; i++){

	/* For left state */
		
		if(Opacity != NULL){
		
			/* We do not need to calculate dSigma here */
			Opacity(Wl[i].d, Wl[i].P / (Wl[i].d * R_ideal),&(Sigma_t), &(Sigma_a), NULL);
		}
		else{
			Sigma_a = W[i-1].Sigma_a;
			Sigma_t = W[i-1].Sigma_t;

		}
/*		
		temperature = Wl[i].P / (Wl[i].d * R_ideal);

		Wl[i].Er = pow(temperature, 4.0);
*/
		
		Wl[i].Sigma_a = Sigma_a;
		Wl[i].Sigma_t = Sigma_t;
		Wl[i].Edd_11 = W[i-1].Edd_11;
		Wl[i].Edd_21 = W[i-1].Edd_21;
		Wl[i].Edd_22 = W[i-1].Edd_22;
		Wl[i].Edd_31 = W[i-1].Edd_31;
		Wl[i].Edd_32 = W[i-1].Edd_32;
		Wl[i].Edd_33 = W[i-1].Edd_33;

		velocity = Wl[i].Vx;

		SFFr = -Crat * Sigma_t;

		alpha = (exp(SFFr * dt * 0.5) - 1.0)/(SFFr * dt * 0.5);

		SourceFr = -Crat * Sigma_t * (Wl[i].Fr1 - velocity * (1.0 + W[i-1].Edd_11) * Wl[i].Er / Crat);


		Wl[i].Fr1 += 0.5 * dt * alpha * SourceFr;
	
		

	/* For the right state */
	
	       
		if(Opacity != NULL){
		
			/* We do not need to calculate dSigma here */
			Opacity(Wr[i].d, Wr[i].P / (Wr[i].d * R_ideal),&(Sigma_t), &(Sigma_a), NULL);
		}
		else{
			Sigma_a = W[i].Sigma_a;
			Sigma_t = W[i].Sigma_t;

		}
/*
		temperature = Wl[i].P / (Wl[i].d * R_ideal);

		Wl[i].Er = pow(temperature, 4.0);
*/		
		Wr[i].Sigma_a = Sigma_a;
		Wr[i].Sigma_t = Sigma_t;
		Wr[i].Edd_11 = W[i].Edd_11;
		Wr[i].Edd_21 = W[i].Edd_21;
		Wr[i].Edd_22 = W[i].Edd_22;
		Wr[i].Edd_31 = W[i].Edd_31;
		Wr[i].Edd_32 = W[i].Edd_32;
		Wr[i].Edd_33 = W[i].Edd_33;
		


		velocity = Wr[i].Vx;

		SFFr = -Crat * Sigma_t;

		alpha = (exp(SFFr * dt * 0.5) - 1.0)/(SFFr * dt * 0.5);

		SourceFr = -Crat * Sigma_t * (Wr[i].Fr1 - velocity * (1.0 + W[i].Edd_11) * Wr[i].Er / Crat);


		Wr[i].Fr1 += 0.5 * dt * alpha * SourceFr;
	
		

	}

/*-----------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */
	
  	if (StaticGravPot != NULL){
    		for (i=il+1; i<=iu; i++) {
      		cc_pos(pG,i,js,ks,&x1,&x2,&x3);

		phicr = (*StaticGravPot)( x1             ,x2,x3);
      		phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
      		phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      		Wl[i].Vx -= dtodx1*(phifc - phicl);
      		Wr[i].Vx -= dtodx1*(phicr - phifc);

    		}
  	}


		
/*---------Step 2c--------------*/ 
/* Calculate the flux */

		 for (i=il+1; i<=iu; i++) {

  		  Ul_x1Face[i] = Prim1D_to_Cons1D(&Wl[i], &Bxi[i]);
   		  Ur_x1Face[i] = Prim1D_to_Cons1D(&Wr[i], &Bxi[i]);

		  /* Copy first order Er and Fr to Ul and Ur */

		  Ul_x1Face[i].Er = W[i-1].Er;
		  Ul_x1Face[i].Fr1 = W[i-1].Fr1;

		  Ur_x1Face[i].Er = W[i].Er;
                  Ur_x1Face[i].Fr1 = W[i].Fr1;
	


		  x1Flux[i].d = dt;
		/* This is used to take dt to calculate alpha. 
                * x1Flux[i].d is recovered in the fluxes function by using Wl */

   		 hlle_thick(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i]);
	
 		 }

/*---
 * Calculate d^{n+1/2} (needed with static potential, or cooling)
 */
		
  		if ((StaticGravPot != NULL))
  		{
    			for (i=il+1; i<=iu-1; i++) {
      			dhalf[i] = pG->U[ks][js][i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
    			}
  		}
	

/*----Step 3------------------
 * Modified Godunov Corrector Scheme   */

	for(i=il+1; i<=iu-1; i++) {
	if(pG->Flagtau[ks][js][i]){

	velocity = U1d[i].Mx / U1d[i].d;

	Sigma_t = U1d[i].Sigma_t;
	Sigma_a = U1d[i].Sigma_a;

	SFFr = -Crat * Sigma_t;
	

/* The Source term */
	SourceFr = -Crat * Sigma_t * (U1d[i].Fr1 - velocity * (1.0 + U1d[i].Edd_11) * U1d[i].Er / Crat);

	InvFr = 1.0 / (1.0 - dt * SFFr);

	pdivFlux = (Real*)&(divFlux);
	pfluxr = (Real*)&(x1Flux[i+1]);
	pfluxl = (Real*)&(x1Flux[i]);

	for(m=0; m<9; m++)
		pdivFlux[m] = (pfluxr[m] - pfluxl[m]) / dx;
		
	/* copy U1d to Uguess, some of the variables will be updated here */
	Uguess = U1d[i];

	Uguess.Fr1 = U1d[i].Fr1 + dt * InvFr * (SourceFr - divFlux.Fr1);
	Uguess.d   = U1d[i].d - dt * divFlux.d;
	Uguess.Mx  = U1d[i].Mx - dt * divFlux.Mx + (U1d[i].Fr1 - Uguess.Fr1) * Prat / Crat;
	Uguess.E   = U1d[i].E + Prat * U1d[i].Er - dt * divFlux.E;

	velocity   = Uguess.Mx / Uguess.d;
	Esum       = Uguess.E - 0.5 * Uguess.d * velocity * velocity;
	/*!< Esum should be only sum of thermal energy and radiation energy */ 
#ifdef RADIATION_MHD
	Esum  -=  0.5 * (Bxc[i] * Bxc[i] + U1d[i].By * U1d[i].By + U1d[i].Bz * U1d[i].Bz);
#endif
	
	temperature= EquState(Uguess.d, Esum, U1d[i].Er);
	/* Use the updated opacity due to the change of temperature and density */
	if(Opacity != NULL)
		Opacity(Uguess.d, temperature, &Sigma_t, &Sigma_a, NULL);

	Uguess.Er  = pow(temperature, 4.0);

	Uguess.E  -= Prat * Uguess.Er;
	

	/*!< Now do the correction step */
		
	SourceFrguess =  -Crat * Sigma_t * (Uguess.Fr1 - velocity * (1.0 + U1d[i].Edd_11) * Uguess.Er / Crat);

	Error = U1d[i].Fr1 + 0.5 * dt * (SourceFr + SourceFrguess) - dt * divFlux.Fr1 - Uguess.Fr1;

	correction = InvFr * Error;

	/*!< Now update the variables */

	pG->U[ks][js][i].Fr1 	= Uguess.Fr1 + correction;
	pG->U[ks][js][i].d 	= Uguess.d;
	pG->U[ks][js][i].M1 	= Uguess.Mx + Prat * (Uguess.Fr1 - pG->U[ks][js][i].Fr1) / Crat;

	Esum		       += 0.5 * Uguess.d * velocity * velocity 
				- 0.5 * pG->U[ks][js][i].M1 * pG->U[ks][js][i].M1 / Uguess.d;

	temperature= EquState(pG->U[ks][js][i].d, Esum, Uguess.Er);
		
	pG->U[ks][js][i].Er 	= pow(temperature, 4.0);

	pG->U[ks][js][i].E 	=  U1d[i].E + Prat * (U1d[i].Er - pG->U[ks][js][i].Er) - dt * divFlux.E;

#ifdef RADIATION_MHD
/* magnetic field is independent of modified Godunov method */
	pG->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    	pG->U[ks][js][i].B3c -= dtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
/* For consistency, set B2i and B3i to cell-centered values.  */
    	pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    	pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;

#endif


	/*Boundary condition is applied in the main.c function*/

/*-----------Finish---------------------*/
	} /* Finish optical thick flag*/
	} /* End big loop i */	


 /* Add gravitational source terms as a Static Potential.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

	  if (StaticGravPot != NULL){
    		for (i=is; i<=ie; i++) {
		if(pG->Flagtau[ks][js][i]){
      		cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      		phic = (*StaticGravPot)((x1            ),x2,x3);
      		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

	        pG->U[ks][js][i].M1 -= dtodx1*dhalf[i]*(phir-phil);


      		pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
    		}
		}
  	}


	/* Now update the reaction coefficient Sigma_t and Sigma_a in the conserved variable 
	 * If function Opacity is set in the problem generator 
	 *
	 */ 

	if(Opacity != NULL){
		for(i=il+1; i<=iu-1; i++) {
		if(pG->Flagtau[ks][js][i]){
			pressure = (pG->U[ks][js][i].E - 0.5 * pG->U[ks][js][i].M1 * pG->U[ks][js][i].M1 / pG->U[ks][js][i].d )
				* (Gamma - 1);
/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[ks][js][i].B1c * pG->U[ks][js][i].B1c + pG->U[ks][js][i].B2c * pG->U[ks][js][i].B2c + pG->U[ks][js][i].B3c * pG->U[ks][js][i].B3c) * (Gamma - 1.0);
#endif
		
			temperature = pressure / (pG->U[ks][js][i].d * R_ideal);
	
			/* We do not need to calculate dSigma here */
			Opacity(pG->U[ks][js][i].d, temperature,&(pG->U[ks][js][i].Sigma_t), &(pG->U[ks][js][i].Sigma_a), NULL);
		}
		}
	}


  return;

}





/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(MeshS *pM)
{
  int size1=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if ((pM->Domain[nl][nd].Grid->Nx[0]) > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
      }
    }
  }
/* We only need one boundary value at each side to solve the matrix */
   

  size1 = size1 + 2*nghost;

  if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;


  if ((U1d       =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
  if ((Ul_x1Face =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
  if ((Ur_x1Face =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
  if ((x1Flux    =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
 

  if ((W  = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;

  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
  {
    if ((dhalf  = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  }
  

#ifdef CYLINDRICAL
  if ((geom_src = (Real*)calloc_1d_array(size1, sizeof(Real))) == NULL)
    goto on_error;
#endif

  return;

/* Destruct is only used when error happens  */

  on_error:
    integrate_destruct();
    ath_error("[integrate_init_1d]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_1d: Free temporary integration arrays  */

void integrate_destruct_1d(void)
{

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);

  
  if (U1d != NULL) free(U1d);
  if (Ul_x1Face != NULL) free(Ul_x1Face);
  if (Ur_x1Face != NULL) free(Ur_x1Face);
  if (x1Flux != NULL) free(x1Flux);

  if (W  != NULL) free(W);
  if (Wl != NULL) free(Wl);
  if (Wr != NULL) free(Wr);

  if (dhalf != NULL) free(dhalf);

#ifdef CYLINDRICAL
  if (geom_src != NULL) free_1d_array((void*)geom_src);
#endif

  return;
}
#endif /* radMHD_INTEGRATOR */
