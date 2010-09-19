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

#if defined(radMHD_INTEGRATOR)
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

/* Radiation matter coupled speed */
static Real *Cspeeds;
/* The matrix coefficient */
static Real **MatrixEuler = NULL;
/* Right hand side of the Matrix equation */
static Real *RHSEuler = NULL;
/* Used to store the results from the Matrix solver */
static int *Ern1 = NULL;
static Real *Ern2 = NULL;

/* Variables at t^{n+1/2} used for source terms */
static Real *dhalf = NULL, *phalf = NULL;

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
	Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
	Real dt=pG->dt, dx=pG->dx1;
	int il,iu, is = pG->is, ie = pG->ie;
  	int i, j, m, n;
	int js = pG->js;
	int ks = pG->ks;
	int Nmatrix;
	Cons1DS Uguess, divFlux;
	Real *pUguess, *pdivFlux, *pfluxl, *pfluxr, *pU1d;


	Real SEE, SErho, SEm;
	Real SPP, alpha;
	Real temperature, velocity, pressure, Sigmas;

	Real Source_Inv[NVAR][NVAR], tempguess[NVAR], Source[NVAR];
	Real Propa_44;
//, Propa[NVAR][NVAR];
	Real Source_guess[NVAR], Errort[NVAR];

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


 	 il = is - 1;
  	 iu = ie + 1;

		

/* Allocate memory space for the Matrix calculation, just used for this grids */
/* Nmatrix is the number of active cells just in this grids */
/* Matrix size should be 2 * Nmatrix, ghost zones are included*/
   	Nmatrix = ie - is + 1 ;
	rad_hydro_init_1d(Nmatrix);


/* In principle, should load a routine to calculate the tensor f */

/* Temperatory variables used to calculate the Matrix  */
  	
  	Real temp1, temp2;
  	Real Matrixtheta[6];
  	Real Matrixphi[6];
  	Real Sigma_s = Sigma_t - Sigma_a;

/* Load 1D vector of conserved variables and calculate the source term */
 	for (i=is-nghost; i<=ie+nghost; i++) {
    		U1d[i].d  = pG->U[ks][js][i].d;
    		U1d[i].Mx = pG->U[ks][js][i].M1;
    		U1d[i].My = pG->U[ks][js][i].M2;
    		U1d[i].Mz = pG->U[ks][js][i].M3;
    		U1d[i].E  = pG->U[ks][js][i].E;
    		U1d[i].Er  = pG->U[ks][js][i].Er;
    		U1d[i].Fluxr1  = pG->U[ks][js][i].Fluxr1;
    		U1d[i].Fluxr2  = pG->U[ks][js][i].Fluxr2;
    		U1d[i].Fluxr3  = pG->U[ks][js][i].Fluxr3;
		U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;
	}
	

	
	
/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */


/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1 =size1, otherwise it is wrong */


	for(i=is; i<=ie+1; i++){
	 	Cspeeds[i-is] = (U1d[i].Edd_11 - U1d[i-1].Edd_11) 
				/ (U1d[i].Edd_11 + U1d[i-1].Edd_11); 		
	}


	for(i=is; i<=ie; i++){
/* E is the total energy. should subtract the kinetic energy and magnetic energy density */
    		pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1);
/* if MHD - 0.5 * Bx * Bx   */

    		temperature = pressure / (U1d[i].d * R_ideal);
      

    		RHSEuler[2*(i-is)]   = U1d[i].Er + Crat * pG->dt * Sigma_a 
				* temperature * temperature * temperature * temperature;
    		RHSEuler[2*(i-is)+1] = U1d[i].Fluxr1 + pG->dt *  Sigma_a
				* temperature * temperature * temperature * temperature * U1d[i].Mx / U1d[i].d;
	}

		

/*--------------------Note--------------------*.
/* Should judge the boundary condition. Here just use the perodic boundary condition first */

/* Step 1b: Setup the Matrix */
		

	for(i=is; i<=ie; i++){
		velocity = U1d[i].Mx / U1d[i].d; 
		Matrixtheta[0] = -Crat * hdtodx1 * (1.0 + Cspeeds[i-is]) * U1d[i-1].Edd_11;
		Matrixtheta[1] = -Crat * hdtodx1 * (1.0 + Cspeeds[i-is]);
		Matrixtheta[2] = 1.0 + Crat * hdtodx1 * (1.0 + Cspeeds[i-is+1]) * U1d[i].Edd_11 
			+ Crat * hdtodx1 * (1.0 - Cspeeds[i-is]) * U1d[i].Edd_11 + Crat * pG->dt * Sigma_a 
			+ pG->dt * (Sigma_a - Sigma_s) * (1.0 + U1d[i].Edd_11) * velocity * velocity / Crat;
		Matrixtheta[3] = Crat * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) 
				- pG->dt * (Sigma_a - Sigma_s) * velocity;
		Matrixtheta[4] = -Crat * hdtodx1 * (1.0 - Cspeeds[i-is+1]) * U1d[i+1].Edd_11;
		Matrixtheta[5] = Crat * hdtodx1 * (1.0 - Cspeeds[i-is+1]);

		Matrixphi[0]	= Matrixtheta[0];
		Matrixphi[1]	= Matrixtheta[1] * U1d[i-1].Edd_11;
		Matrixphi[2]	= Crat * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) * U1d[i].Edd_11 
				- pG->dt * Sigma_t * (1.0 + U1d[i].Edd_11) * velocity + pG->dt * Sigma_a * velocity;
		Matrixphi[3]	= 1.0 + Crat * hdtodx1 * (2.0 + Cspeeds[i-is+1] 
		- Cspeeds[i-is]) * U1d[i].Edd_11 + Crat * pG->dt * Sigma_t;
		Matrixphi[4]	= -Matrixtheta[4];
		Matrixphi[5]	= -Matrixtheta[5] * U1d[i+1].Edd_11;
	
/* For perodic boundary condition....................*/
		if(i == is) {
			MatrixEuler[0][0] = Matrixtheta[2];
			MatrixEuler[0][1] = Matrixtheta[3];
			MatrixEuler[0][3] = Matrixtheta[4];
			MatrixEuler[0][4] = Matrixtheta[5];
			MatrixEuler[0][2*Nmatrix-2] = Matrixtheta[0];
			MatrixEuler[0][2*Nmatrix-1] = Matrixtheta[1];
			MatrixEuler[1][0] = Matrixphi[2];
			MatrixEuler[1][1] = Matrixphi[3];
			MatrixEuler[1][3] = Matrixphi[4];
			MatrixEuler[1][4] = Matrixphi[5];
			MatrixEuler[1][2*Nmatrix-2] = Matrixphi[0];
			MatrixEuler[1][2*Nmatrix-1] = Matrixphi[1];
		}//end if
		else if(i == ie){
			MatrixEuler[2*Nmatrix-2][0] = Matrixtheta[4];
			MatrixEuler[2*Nmatrix-2][1] = Matrixtheta[5];
			MatrixEuler[2*Nmatrix-2][2*Nmatrix-4] = Matrixtheta[0];
			MatrixEuler[2*Nmatrix-2][2*Nmatrix-3] = Matrixtheta[1];
			MatrixEuler[2*Nmatrix-2][2*Nmatrix-2] = Matrixtheta[2];
			MatrixEuler[2*Nmatrix-2][2*Nmatrix-1] = Matrixtheta[3];

			MatrixEuler[2*Nmatrix-1][0] = Matrixphi[4];
			MatrixEuler[2*Nmatrix-1][1] = Matrixphi[5];
			MatrixEuler[2*Nmatrix-1][2*Nmatrix-4] = Matrixphi[0];
			MatrixEuler[2*Nmatrix-1][2*Nmatrix-3] = Matrixphi[1];
			MatrixEuler[2*Nmatrix-1][2*Nmatrix-2] = Matrixphi[2];
			MatrixEuler[2*Nmatrix-1][2*Nmatrix-1] = Matrixphi[3];
		}// end else if
		else {
			for(j=0; j<6; j++){
			MatrixEuler[2*(i-is)][j+2*(i-is-1)] = Matrixtheta[j];
			MatrixEuler[2*(i-is)+1][j+2*(i-is-1)] = Matrixphi[j];
			}// end for j loop
		}// end else
	}// end for i loop


		
	
	
/* Step 1c: Solve the Matrix equation */
		ludcmp(MatrixEuler,2*Nmatrix,Ern1,Ern2);
		
		lubksb(MatrixEuler,2*Nmatrix,Ern1,RHSEuler);
		
		
	for(i=is;i<=ie;i++){
		pG->U[ks][js][i].Er	= RHSEuler[2*(i-is)];
		pG->U[ks][js][i].Fluxr1	= RHSEuler[2*(i-is)+1];
		U1d[i].Er		= RHSEuler[2*(i-is)];
		U1d[i].Fluxr1		= RHSEuler[2*(i-is)+1];
	}
	/* May need to update Edd_11 */


/* Update the ghost zones for Periodic boundary condition to be used later */
/* Boundary condition for pG grids is applied after this loop is finished */

	 /* Inner boundary condition */
		U1d[is-1].Er = U1d[ie].Er;
		U1d[is-1].Fluxr1 = U1d[ie].Fluxr1;
		U1d[is-2].Er = U1d[ie-1].Er;
		U1d[is-2].Fluxr1 = U1d[ie-1].Fluxr1;

	/* Outer boundary condition */
		U1d[ie+1].Er = U1d[is].Er;
		U1d[ie+1].Fluxr1 = U1d[is].Fluxr1;
		U1d[ie+2].Er = U1d[is+1].Er;
		U1d[ie+2].Fluxr1 = U1d[is+1].Fluxr1;
		


/*-----Step 2a----------------
 *  Calcualte the left and right state */
	 for (i=is-nghost; i<=ie+nghost; i++) {
    		W[i] = Cons1D_to_Prim1D(&U1d[i], &Bxc[i]);
	  }

  	lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);

	
/*------Step 2b: Add source terms to the left and right state--------*/
	for(i=il+1; i<=ie+1; i++){

	/* For left state */
		pressure = W[i-1].P;
		temperature = pressure / (U1d[i-1].d * R_ideal);
		velocity = U1d[i-1].Mx / U1d[i-1].d;

	Source[1] = -Prat * (-Sigma_t * (U1d[i-1].Fluxr1/U1d[i-1].d 
	- (1.0 + U1d[i-1].Edd_11) * velocity * U1d[i-1].Er / (Crat * U1d[i-1].d))	
	+ Sigma_a * velocity * (temperature * temperature * temperature * temperature - U1d[i-1].Er)/(Crat*U1d[i-1].d));
	Source[4] = Source[1] * U1d[i-1].d * velocity * (Gamma - 1.0)
	-(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
		- U1d[i-1].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i-1].Fluxr1 - (1.0 + U1d[i-1].Edd_11) * velocity * U1d[i-1].Er / Crat)/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[i-1].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		//Propa[4][0] = (1.0 - alpha) * W[i-1].P / U1d[i-1].d;
		Propa_44 = alpha;

		Wl[i].P += dt * Propa_44 * Source[4] * 0.5;

	/* For the right state */
	
	
		pressure = W[i].P;
		temperature = pressure / (U1d[i].d * R_ideal);
		velocity = U1d[i].Mx / U1d[i].d;

	Source[1] = -Prat * (-Sigma_t * (U1d[i].Fluxr1/U1d[i].d 
	- (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / (Crat * U1d[i].d))	
	+ Sigma_a * velocity * (temperature * temperature * temperature * temperature - U1d[i].Er)/(Crat*U1d[i].d));
	Source[4] = Source[1] * U1d[i].d * velocity * (Gamma - 1.0)
	-(Gamma - 1.0) * Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
		- U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i].Fluxr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)/Crat); 
		SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_a * temperature * temperature 
			* temperature /(U1d[i].d * R_ideal);
		if(fabs(SPP * dt * 0.5) > 0.001)
		alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
		else 
		alpha = 1.0 + 0.25 * SPP * dt;
		/* In case SPP * dt  is small, use expansion expression */	
		//Propa[4][0] = (1.0 - alpha) * W[i].P / U1d[i].d;
		Propa_44 = alpha;

		Wr[i].P += dt * Propa_44 * Source[4] * 0.5;

	}

		
/*---------Step 2c--------------*/ 
		 for (i=il+1; i<=iu; i++) {
  		  Ul_x1Face[i] = Prim1D_to_Cons1D(&Wl[i], &Bxi[i]);
   		  Ur_x1Face[i] = Prim1D_to_Cons1D(&Wr[i], &Bxi[i]);
		  x1Flux[i].d = dt;
		/* This is used to take dt to calculate alpha. 
                * x1Flux[i].d is recovered in the fluxes function by using Wl */

   		 fluxes(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i]);
		
 		 }


	

/*----Step 3------------------
 * Modified Godunov Corrector Scheme   */
	for(i=il+1; i<=iu-1; i++) {

	pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (U1d[i].d * R_ideal);
	velocity = U1d[i].Mx / U1d[i].d;

/* The Source term */
	SEE = 4.0 * Sigma_a * temperature * temperature * temperature * (Gamma - 1.0)/ U1d[i].d;
	SErho = SEE * (-U1d[i].E/U1d[i].d + velocity * velocity);	
	SEm = -SEE * velocity;
	
	Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
	Source_Inv[4][1] = -dt * Prat * Crat * SEm/(1.0 + dt * Prat * Crat * SEE);
	Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	
	Source[1] = -Prat * (-Sigma_t * (U1d[i].Fluxr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)
	+ Sigma_a * velocity * (temperature * temperature * temperature * temperature - U1d[i].Er)/Crat);
	Source[4] = -Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
		- U1d[i].Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (U1d[i].Fluxr1 - (1.0 + U1d[i].Edd_11) * velocity * U1d[i].Er / Crat)/Crat); 

	pdivFlux = (Real*)&(divFlux);
	pfluxr = (Real*)&(x1Flux[i+1]);
	pfluxl = (Real*)&(x1Flux[i]);

	for(m=0; m<NVAR-4; m++)
		pdivFlux[m] = (pfluxr[m] - pfluxl[m]) / dx;
		
	
	for(n=0; n<NVAR-4; n++) {
		tempguess[n] = 0.0;
		for(m=0; m<NVAR-4; m++){
		tempguess[n] += dt * Source_Inv[n][m] * (Source[m] - pdivFlux[m]);
		}
	}
	
	pUguess = (Real*)&(Uguess);
	pU1d = (Real*)&(U1d[i]);	

	for(m=0; m<NVAR-4; m++)
		pUguess[m]= pU1d[m] + tempguess[m];


	pressure = (Uguess.E - 0.5 * Uguess.Mx * Uguess.Mx / Uguess.d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (Uguess.d * R_ideal);
	velocity = Uguess.Mx / Uguess.d;


	Source_guess[1] = -Prat * (-Sigma_t * (Uguess.Fluxr1 - (1.0 + U1d[i].Edd_11) * velocity * Uguess.Er / Crat)
	+ Sigma_a * velocity * (temperature * temperature * temperature * temperature - Uguess.Er)/Crat);
	Source_guess[4] = -Prat * Crat * (Sigma_a * (temperature * temperature * temperature * temperature 
		- Uguess.Er) + (Sigma_a - (Sigma_t - Sigma_a)) * velocity
		* (Uguess.Fluxr1 - (1.0 + U1d[i].Edd_11) * velocity * Uguess.Er / Crat)/Crat); 

	for(m=0; m<NVAR-4; m++)
		Errort[m] = pU1d[m] + 0.5 * dt * (Source_guess[m] + Source[m]) 
			- dt * pdivFlux[m] - pUguess[m];

	for(m=0; m<NVAR-4; m++){
		tempguess[m]=0.0;
		for(n=0; n<NVAR; n++){
			tempguess[m] += Source_Inv[m][n] * Errort[n];
		}
		pU1d[m] = pUguess[m] + tempguess[m];

	}
	
	/* Update the quantity in the Grids */
	pUguess = (Real*)&(pG->U[ks][js][i]);
	for(m=0; m<NVAR-4; m++) pUguess[m] = pU1d[m];	


	/*Boundary condition is applied in the main.c function*/

/*-----------Finish---------------------*/
	}	

  
	/* Free the temporary variables just used for this grid calculation*/
	rad_hydro_destruct_1d(Nmatrix);
	
  return;

	on_error:
	integrate_destruct();
	radiation_destruct_1d(Nmatrix);
	ath_error("[integrate_1d_radMHD]: malloc returned a NULL pointer\n");


}




/*-------------------------------------------------------------------------*/
/* radiation_init_1d(): function to allocate memory used just for radiation variables */
/* radiation_destruct_1d(): function to free memory */
void rad_hydro_init_1d(int Ngrids)
{

	int i, j;

	if ((Cspeeds = (Real*)malloc((Ngrids+1)*sizeof(Real))) == NULL) goto on_error;
	if ((RHSEuler = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;
	if ((Ern1 = (int*)malloc(2*Ngrids*sizeof(int))) == NULL) goto on_error;
	if ((Ern2 = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;

	if ((MatrixEuler = (Real**)malloc(2*Ngrids*sizeof(Real*))) == NULL) goto on_error;
	for(i=0; i<2*Ngrids; i++)
		if ((MatrixEuler[i] = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;

	 /* Initialize the matrix */
	for(i=0; i<2*Ngrids; i++)
		for(j=0; j<2*Ngrids; j++)
			MatrixEuler[i][j] = 0.0;
	
	return;

	on_error:
    	integrate_destruct();
	rad_hydro_destruct_1d(Ngrids);
	ath_error("[radiation_init_1d]: malloc returned a NULL pointer\n");
}


void rad_hydro_destruct_1d(int Ngrids)
{
	int i;
	if (Cspeeds 	!= NULL) free(Cspeeds);
	if (RHSEuler	!= NULL) free(RHSEuler);
	if (Ern1	!= NULL) free(Ern1);
	if (Ern2	!= NULL) free(Ern2);
	
	for(i=0; i<2*Ngrids; i++)
		if (MatrixEuler[i] != NULL) free(MatrixEuler[i]);

	if (MatrixEuler != NULL) free(MatrixEuler);
}


/*----------------------------------------------------------------------------*/
/* integrate_init_1d: Allocate temporary integration arrays */

void integrate_init_1d(MeshS *pM)
{
  int size1=0,nl,nd, i;

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

#ifdef CYLINDRICAL
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
  {
    if ((dhalf  = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  }
  if(CoolingFunc != NULL){
    if ((phalf  = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
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
  if (phalf != NULL) free(phalf);

#ifdef CYLINDRICAL
  if (geom_src != NULL) free_1d_array((void*)geom_src);
#endif

  return;
}
#endif /* radMHD_INTEGRATOR */
