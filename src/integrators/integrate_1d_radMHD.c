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
/* Tensor to relate radiation energy density and pressure */
/* To be modified later */
static Real *frad = NULL;
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
	Real temperature, velocity, pressure, Sigmas;

	Real Matrix_source[NWAVE], Matrix_source_Inv[NWAVE][NWAVE], tempguess[NWAVE];
	Real Matrix_source_guess[NWAVE], Errort[NWAVE];

	/* Initialize them to be zero */
	for(i=0; i<NWAVE; i++){
		Matrix_source[i] = 0.0;
		Matrix_source_guess[i] = 0.0;
		Errort[i] = 0.0;
		for(j=0; j<NWAVE; j++) {
			Matrix_source_Inv[i][j] = 0.0;
		if(i==j) {
		 Matrix_source_Inv[i][j] = 1.0;
		}
	}
	}


 	 il = is - 1;
  	 iu = ie + 1;

/* Allocate memory space for the Matrix calculation, just used for this grids */
/* Nmatrix is the number of active cells just in this grids */
/* Matrix size should be 2 * Nmatrix */
   	Nmatrix = ie - is + 1;
	radiation_init_1d(Nmatrix);


/* In principle, should load a routine to calculate the tensor f */

/* Temperatory variables used to calculate the Matrix  */
  	
  	Real temp1, temp2;
  	Real Matrixtheta[6];
  	Real Matrixphi[6];
  	Sigmas = Sigmat - Sigmaa;

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
	}
	



/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */


/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1=size1, otherwise it is wrong */


/* Now fra1D is a constant. To be improved later  */
	for(i=is-1; i<=ie+1; i++){
    		frad[i-is+1] = pG->fra1D;
	}

	for(i=is; i<=ie+1; i++){

    		Cspeeds[i-is] = (frad[i-is+1] - frad[i-is]) / (frad[i-is+1] + frad[i-is]); 	
	}


	for(i=is; i<=ie; i++){
/* E is the total energy. should subtract the kinetic energy and magnetic energy density */
    		pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1);
/* if MHD - 0.5 * Bx * Bx   */

    		temperature = pressure / (U1d[i].d * Ridealgas);
      

    		RHSEuler[2*(i-is)]   = U1d[i].Er + Cratio * pG->dt * Sigmaa 
				* temperature * temperature * temperature * temperature;
    		RHSEuler[2*(i-is)+1] = U1d[i].Fluxr1 + pG->dt * Cratio * Sigmaa
				* temperature * temperature * temperature * temperature * U1d[i].Mx / U1d[i].d;
	}

/*--------------------Note--------------------*.
/* Should judge the boundary condition. Here just use the perodic boundary condition first */

/* Step 1b: Setup the Matrix */
	for(i=is; i<=ie; i++){
		velocity = U1d[i].Mx / U1d[i].d; 
		Matrixtheta[0] = -Cratio * hdtodx1 * (1.0 + Cspeeds[i-is]) * frad[i-is];
		Matrixtheta[1] = -Cratio * hdtodx1 * (1.0 + Cspeeds[i-is]);
		Matrixtheta[2] = 1.0 + Cratio * hdtodx1 * (1.0 + Cspeeds[i-is+1]) * frad[i-is+1] 
			+ Cratio * hdtodx1 * (1.0 - Cspeeds[i-is]) * frad[i-is+1] + Cratio * pG->dt * Sigmaa 
			+ pG->dt * (Sigmaa - Sigmas) * (1.0 + frad[i-is+1]) * velocity * velocity / Cratio;
		Matrixtheta[3] = Cratio * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) 
				- pG->dt * (Sigmaa - Sigmas) * velocity;
		Matrixtheta[4] = -Cratio * hdtodx1 * (1.0 - Cspeeds[i-is+1]) * frad[i-is+2];
		Matrixtheta[5] = Cratio * hdtodx1 * (1.0 - Cspeeds[i-is+1]);

		Matrixphi[0]	= Matrixtheta[0];
		Matrixphi[1]	= Matrixtheta[1] * frad[i-is];
		Matrixphi[2]	= Cratio * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) * frad[i-is+1] 
				- pG->dt * Sigmat * (1.0 + frad[i-is+1]) * velocity + pG->dt * Sigmaa * velocity;
		Matrixphi[3]	= 1.0 + Cratio * hdtodx1 * (2.0 + Cspeeds[i-is+1] - Cspeeds[i-is]) * frad[i-is+1] 
				+ Cratio * pG->dt * Sigmat;
		Matrixphi[4]	= -Matrixtheta[4];
		Matrixphi[5]	= -Matrixtheta[5] * frad[i-is+2];
	
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
/* Update the ghost zones for Periodic boundary condition*/
	
	/* NOTE: Need to update the ghost zone for Er and Fluxr1. Not suer how now. */
	/* May call the boundary function  */

/*-----Step 2----------------
 *  Calcualte the left and right state and fluxes */
	lr_states_cons(pG,U1d,Bxc,pG->dt,pG->dx1,il+1,iu-1,Ul_x1Face,Ur_x1Face,1);
	

	 for (i=il+1; i<=iu; i++) {
	    Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[i], &Bxi[i]);
	    Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[i], &Bxi[i]);

	    rad_fluxes(Ul_x1Face[i],Ur_x1Face[i],Wl[i],Wr[i],Bxi[i],&x1Flux[i],dt);
  }

/*----Step 3------------------
 * Modified Godunov Corrector Scheme   */
	for(i=il+1; i<=iu-1; i++) {

	pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (U1d[i].d * Ridealgas);
	velocity = U1d[i].Mx / U1d[i].d;

/* The Source term */
	SEE = 4.0 * Sigmaa * temperature * temperature * temperature * (Gamma - 1.0)/ U1d[i].d;
	SErho = SEE * (-U1d[i].E/U1d[i].d + velocity * velocity);	
	SEm = -SEE * velocity;
	
	Matrix_source_Inv[4][0] = -dt * Pratio * Cratio * SErho/(1.0 + dt * Pratio * Cratio * SEE);
	Matrix_source_Inv[4][1] = -dt * Pratio * Cratio * SEm/(1.0 + dt * Pratio * Cratio * SEE);
	Matrix_source_Inv[4][4] = 1.0 / (1.0 + dt * Pratio * Cratio * SEE);
	
	Matrix_source[1] = -Pratio * (-Sigmat * (U1d[i].Fluxr1 - (1.0 + pG->fra1D) * velocity * U1d[i].Er / Cratio)
	+ Sigmaa * velocity * (temperature * temperature * temperature * temperature - U1d[i].Er)/Cratio);
	Matrix_source[4] = -Pratio * Cratio * (Sigmaa * (temperature * temperature * temperature * temperature 
		- U1d[i].Er) + (Sigmaa - (Sigmat - Sigmaa)) * velocity
		* (U1d[i].Fluxr1 - (1.0 + pG->fra1D) * velocity * U1d[i].Er / Cratio)/Cratio); 

	pdivFlux = (Real*)&(divFlux);
	pfluxr = (Real*)&(x1Flux[i+1]);
	pfluxl = (Real*)&(x1Flux[i]);

	for(m=0; m<NWAVE; m++)
		pdivFlux[m] = (pfluxr[m] - pfluxl[m]) / dx;
		
	
	for(n=0; n<NWAVE; n++) {
		tempguess[n] = 0.0;
		for(m=0; m<NWAVE; m++){
		tempguess[n] += dt * Matrix_source_Inv[n][m] * (Matrix_source[m] - pdivFlux[m]);
		}
	}
	
	pUguess = (Real*)&(Uguess);
	pU1d = (Real*)&(U1d[i]);	

	for(m=0; m<NWAVE; m++)
		pUguess[m]= pU1d[m] + tempguess[m];


	pressure = (Uguess.E - 0.5 * Uguess.Mx * Uguess.Mx / Uguess.d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (Uguess.d * Ridealgas);
	velocity = Uguess.Mx / Uguess.d;


	Matrix_source_guess[1] = -Pratio * (-Sigmat * (Uguess.Fluxr1 - (1.0 + pG->fra1D) * velocity * Uguess.Er / Cratio)
	+ Sigmaa * velocity * (temperature * temperature * temperature * temperature - Uguess.Er)/Cratio);
	Matrix_source_guess[4] = -Pratio * Cratio * (Sigmaa * (temperature * temperature * temperature * temperature 
		- Uguess.Er) + (Sigmaa - (Sigmat - Sigmaa)) * velocity
		* (Uguess.Fluxr1 - (1.0 + pG->fra1D) * velocity * Uguess.Er / Cratio)/Cratio); 

	for(m=0; m<NWAVE; m++)
		Errort[m] = pU1d[m] + 0.5 * dt * (Matrix_source_guess[m] + Matrix_source[m]) 
			- dt * pdivFlux[m] - pUguess[m];

	for(m=0; m<NWAVE; m++){
		tempguess[m]=0.0;
		for(n=0; n<NWAVE; n++){
			tempguess[m] += Matrix_source_Inv[m][n] * Errort[n];
		}
		pU1d[m] = pUguess[m] + tempguess[m];

	}
	
	/* Update the quantity in the Grids */
	pUguess = (Real*)&(pG->U[ks][js][i]);
	for(m=0; m<NWAVE; m++) pUguess[m] = pU1d[m];		

	/*Apply the boundary condition */

/*-----------Finish---------------------*/
}	

  
	/* Free the temporary variables just used for this grid calculation*/
	radiation_destruct_1d(Nmatrix);
	
  return;

	on_error:
	integrate_destruct();
	radiation_destruct_1d(Nmatrix);
	ath_error("[integrate_1d_radMHD]: malloc returned a NULL pointer\n");


}




/*-------------------------------------------------------------------------*/
/* radiation_init_1d(): function to allocate memory used just for radiation variables */
/* radiation_destruct_1d(): function to free memory */
void radiation_init_1d(int Ngrids)
{

	int i, j;

	if ((Cspeeds = (Real*)malloc((Ngrids+1)*sizeof(Real))) == NULL) goto on_error;
	if ((frad = (Real*)malloc((Ngrids+2)*sizeof(Real))) == NULL) goto on_error;
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
	radiation_destruct_1d(Ngrids);
	ath_error("[radiation_init_1d]: malloc returned a NULL pointer\n");
}


void radiation_destruct_1d(int Ngrids)
{
	int i;
	if (Cspeeds 	!= NULL) free(Cspeeds);
	if (frad	!= NULL) free(frad);
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
