#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_fullrad.c
 * copy from hydro_to_rad of the radiation-transfer
 * PURPOSE:  Contains functions for updating the RadGrid using the conserved 
 *           variables in Grid (hydro_to_rad) and for computing
 *           the radiation source term and updating the material energy in 
 *           Grid (rad_to_hydro).
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   hydro_to_rad()
 *   rad_to_hydro()
 *============================================================================*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER


void GetSource(const Real dt, const Real Tnew, const Real SigmaB, const Real SigmaI, const Real Jsource, const Real I0, Real *heatcool);
void GetTnew(const Real dt,const Real d,const Real Tgas, const Real J0, const Real Sigma[4], Real *heatcool, Real *Tnew);

/*----------------------------------------------------------------------------*/
/* hydro_to_rad:  */

void hydro_to_fullrad(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k,ifr,m;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;
  int nang = pRG->nang, noct = pRG->noct, l, n;
  Real d, etherm, Tgas, Tnew, J0, AngleV, Jsource;
  Real Sigma[4];
  Real mu[3];

/* Assumes ghost zone conserved variables have been set by
 * bvals routines.  These values are used to set B, chi, eps,
 * etc. so loops include RadGrid ghost zones*/
  if (pG->Nx[0] > 1) {
    ioff = nghost - Radghost; 
    il -= Radghost; 
    iu += Radghost;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - Radghost; 
    jl -= Radghost; 
    ju += Radghost; 
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - Radghost; 
    kl -= Radghost; 
    ku += Radghost;
  } else koff = 0;

/* Compute radiation variables from conserved variables */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
	ig = i + ioff;

	/* ------------------------------------*/
	/* First, update the opacity */
	for(ifr=0; ifr<nf; ifr++) {
	  get_full_opacity(pG,ifr,ig,jg,kg,&(pRG->R[ifr][k][j][i].Sigma[0]));
	}

	/*-------------------------------------------*/


	/* Compute gas temperature and store for later use */
	d = pG->U[kg][jg][ig].d;
	etherm = pG->U[kg][jg][ig].E - (0.5/d) * ( SQR(pG->U[kg][jg][ig].M1) +
		 SQR(pG->U[kg][jg][ig].M2) + SQR(pG->U[kg][jg][ig].M3) );
#ifdef MHD
	etherm -= 0.5 * (SQR(pG->U[kg][jg][ig].B1c) + SQR(pG->U[kg][jg][ig].B2c) + 
			 SQR(pG->U[kg][jg][ig].B3c) );
#endif
	pG->tgas[kg][jg][ig] = MAX(etherm * Gamma_1 / (d * R_ideal),0.0);
        Tgas = pG->tgas[kg][jg][ig];

	/* Calculate the frequency weighted J and opacity */

	J0 = 0.0;
	
	for(l=0; l<4; l++)
		Sigma[l] = 0.0;

	for(ifr=0; ifr<nf; ifr++){
		J0 += pRG->R[ifr][k][j][i].J * pRG->wnu[ifr];
		for(l=0; l<4; l++){
			Sigma[l] += pRG->R[ifr][k][j][i].Sigma[l];
		}

	}

	/* Calculate the new gas temperature */

	GetTnew(pG->dt,d,Tgas,J0,Sigma,&(pG->Radheat[kg][jg][ig]),&Tnew);

	/* Save the new gas temperature at pG->tgas */
	pG->tgas[kg][jg][ig] = Tnew;

	/* Set the momentum source term -Prat * v sigmaa(T^4 - Er)/C*/
	/* pG->Radheat is -Prat * Crat * sigmaa(T^4 - Er) */
	/* So the momentum source is Energy*v/Crat^2 */
	/* Also add momentum source term Prat Sigmaa Fr  */
	pG->Frsource[kg][jg][ig][0] = pG->U[kg][jg][ig].M1 * pG->Radheat[kg][jg][ig] / (d * Crat * Crat);
	pG->Frsource[kg][jg][ig][1] = pG->U[kg][jg][ig].M2 * pG->Radheat[kg][jg][ig] / (d * Crat * Crat);
	pG->Frsource[kg][jg][ig][2] = pG->U[kg][jg][ig].M3 * pG->Radheat[kg][jg][ig] / (d * Crat * Crat);

	/* The momentum source term is added in function UpdateRT when the new H is updated */
	/* We do not try to conserve the momentums because of the stiff scattering source terms */


       }/* Finish i */
      }/* Finish j */
    }/* Finish k */


	for(ifr=0; ifr<nf; ifr++) {	
		/* Now calculate the energy change due to for each array with this new temperature */
	  for(l=0; l<noct; l++){
		for(n=0; n<nang; n++){
		   for(k=kl; k<=ku; k++){
			kg = k + koff;
		   for(j=jl; j<=ju; j++){
			 jg = j + joff;
		   for(i=il; i<=iu; i++){
			ig = i + ioff;

			/* calculate AngleV */

#ifdef CYLINDRICAL
			for(m=0; m<3; m++)
                                mu[m] = pRG->Rphimu[l][n][k][j][i][m];
#else
			for(m=0; m<3; m++)
                                mu[m] = pRG->mu[l][n][k][j][i][m];	

#endif



			AngleV = (pG->U[kg][jg][ig].M1 * mu[0] + pG->U[kg][jg][ig].M2 * mu[1] + pG->U[kg][jg][ig].M3 * mu[2])/pG->U[kg][jg][ig].d;			

			 /* Now add the term 3n\dot vsigma (T^4/4pi - J) */
                        if(Prat > TINY_NUMBER)
                                Jsource = -pG->Radheat[kg][jg][ig] * 3.0 * AngleV/(Prat * Crat * 4.0 * PI);
			else
				Jsource = pG->dt * Crat * (pRG->R[ifr][k][j][i].Sigma[0] * pG->tgas[kg][jg][ig] * pG->tgas[kg][jg][ig] * pG->tgas[kg][jg][ig] * pG->tgas[kg][jg][ig]/(4.0*PI) - pRG->R[ifr][k][j][i].Sigma[1] * pRG->R[ifr][k][j][i].J)/(1.0 + pG->dt * Crat * pRG->R[ifr][k][j][i].Sigma[1]);

			GetSource(pG->dt, pG->tgas[kg][jg][ig], pRG->R[ifr][k][j][i].Sigma[0], pRG->R[ifr][k][j][i].Sigma[1], Jsource, pRG->imu[ifr][l][n][k][j][i],&(pRG->heatcool[ifr][l][n][k][j][i]));	
			/* Now calculate the momentum source term due to this partially updated I */
			for(m=0; m<3; m++){
				pG->Frsource[kg][jg][ig][m] += Prat * pRG->R[ifr][k][j][i].Sigma[1] * (pRG->imu[ifr][l][n][k][j][i] + pRG->heatcool[ifr][l][n][k][j][i] - Jsource) * pRG->mu[l][n][k][j][i][m] *  pRG->wmu[n][k][j][i] * 4.0 * PI;  				
			}
		   }/* Finish i */
                   }/* finish j */
                   }/* Finish k */

		}/* Finish number of angles */
	  }/* Finish number of octant */
	}/* Finish frequency */
      

  return;
}



/* function to calculate the energy exchange between J and B implicitly to handle the short thermalize time scale */
/* This function assumes that opacity does not change during this time step */
/* But the thermal function must be known in order to calculate implicitly */
/* The opacity is the frequency weighted opacity */
void GetTnew(const Real dt,const Real d,const Real Tgas, const Real J0, const Real Sigma[4], Real *heatcool, Real *Tnew)
{
	/* This implicit function assumes that thermal emission is Tgas^4/(4Pi), 
	*  isotropic in every direction */

	/* The implicit source terms solv the energy excahnge due to absorption */
 	/* As well as attentuiation of the source terms */

	/*******************************************/
	/* For absorption, we solve : 
	 * 4*pI * dJ/dt = CSigma_a(Tgas^4 - 4* Pi * J)
	 * de/dt = -PratC Sigma_a(Tgas^4 - 4*PI * J)
         */
	Real SigmaB, SigmaI, Sigmas;
	Real Tr;
	Real coef1, coef2, coef3, coef4, Ersum, pressure, Er0;

	SigmaB  = Sigma[0]; 
	SigmaI  = Sigma[1];
	Sigmas  = Sigma[2];

	/* SigmasJ and SigmasI must be the same in order to conserve energy locally */
	
	Er0 = J0 * 4.0 * PI;

	/* Negative radiation energy density */
	if(Er0 < 0.0 || Prat < TINY_NUMBER || ((SigmaB < TINY_NUMBER) && (SigmaI < TINY_NUMBER))){
		*heatcool = 0.0;
		*Tnew = Tgas;
		
		return;
	}
	else{

		/*----------------------------------------------------------*/
		/* First, energy source due to absorption opacity */

/*		The pow function can be very slow sometimes when Er0 is close to 1 
		Tr = pow(Er0, 0.25);
		
*/
		Tr = sqrt(Er0);
		Tr = sqrt(Tr);

		pressure = Tgas * d * R_ideal;
		Ersum = pressure / (Gamma - 1.0) + Prat * Er0;
	
		/* Here assume input gas temperature and Er0 is positive */
		if(Tgas < 0.0)
			ath_error("[FullRad_GetSource]: Negative gas temperature: %e!n\n",Tgas);
	
		   
		   coef1 = dt * Prat * Crat * SigmaB;
		   coef2 = d * R_ideal * (1.0 + dt * SigmaI * Crat) / (Gamma - 1.0);
		   coef3 = -pressure / (Gamma - 1.0) - dt * SigmaI * Crat * Ersum;
		   coef4 = 0.0;

		if(coef1 < 1.e-20){
			(*Tnew) = -coef3 / coef2;
		}
		else{
		   
		  
		   if(Tgas > Tr){			
			   (*Tnew) = rtsafe(Tequilibrium, Tr * (1.0 - 0.01), Tgas * (1.0 + 0.01), 1.e-12, coef1, coef2, coef3,coef4); 			   
		   }
		   else{
			   (*Tnew) = rtsafe(Tequilibrium, Tgas * (1.0 - 0.01), Tr * (1.0 + 0.01), 1.e-12, coef1, coef2, coef3, coef4);
		   }	
		}		


	
		*heatcool = d * ((*Tnew) - Tgas) * R_ideal/(Gamma - 1.0);
		

		
	}

}


void GetSource(const Real dt, const Real Tnew, const Real SigmaB, const Real SigmaI, const Real Jsource, const Real I0, Real *heatcool)
{
	Real Inew, Tnew4;
	
	Tnew4 = Tnew * Tnew * Tnew * Tnew;

	/* Negative radiation energy density */
	if(Prat < TINY_NUMBER){
		*heatcool = 0.0;		
		
		return;
	}
	else{


		Inew = (I0 + Jsource + dt * Crat * SigmaB * Tnew4 / (4.0 * PI))/(1.0 + dt * Crat * SigmaI);

		(*heatcool) = Inew - I0;	
		/*----------------------------------------------------------*/
		/* First, the energy source due to scattering opacity */
		/* We solve the equation dI/dt = csigma_s(J - I) */
		/* J is held to be a constant during this step */
		/* The formal solution is I(t) = (I0 - J) exp(-dt*csigmas) + J */
		/* Because for each cell, average of I0 over different angles is J */
		/* So total energy is conserved for the scattering process */

	/*	(*Scat) = (1.0 - exp(-Crat * Sigmas * dt)) * (J - I0);
	*/

	}

}

#endif /* FULL_RADIATION_TRANSFER */
