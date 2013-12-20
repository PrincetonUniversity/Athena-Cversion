#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_fullrad.c
 * copy from hydro_to_rad of the radiation-transfer
 * PURPOSE:  Calculate the opacity in the radiation grid and 
 *			 calculate the gas temperature for later use.
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
  int i,j,k,ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;
  Real d, etherm;


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
		  
		  
	  }/* end i */
	}/* end j */
  }/* end k */
		  
		  
		  



  return;
}


/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to absorption opacity */


void RadAsource(const int i, const int j, const int k, const int N, RadGridS *pRG, GridS *pG, Real **Tcoef, Real **Coefn, Real **sol)
{

	int l, n, ifr;
	
	
	Real Sigma, SigmaI;
	Real Er, Ersource, weightJ;
	Real Fr[3], Pr[6], Vel[3], Velguess[3], Msource[3], DeltaKin[3], CoFr[3]; /* velocity source terms and change of kinetic energy */
	Real dt = pG->dt;
	Real Tgas, Tgas4;
	Real rho;


    int nDim;
    
    nDim = 1;
	for (l=1; l<3; l++) if (pRG->Nx[l]>1) nDim++;
    
    int koff, joff, ioff;
    ioff = nghost - Radghost;

    if(nDim > 1)
        joff = ioff;
    else
        joff = 0;
        
    if(nDim > 2)
        koff = ioff;
    else
        koff = 0;

	
			
	/* loop through all l and n, nelements = nang * noct */
	/* weight for each ray is alerady included in sol[][][] */
	for(ifr=0; ifr<pRG->nf; ifr++){
		/* calculate Er and Fr used in the update */
		Er = 0.0;
		for(l=0; l<3; l++)
			Fr[l] = 0.0;
			
		for(l=0; l<6; l++)
			Pr[l] = 0.0;
	
		for(n=0; n<N; n++){
			weightJ = sol[ifr][n];
			
			Er += weightJ;
				
			for(l=0; l<3; l++)
				Fr[l] += weightJ * Coefn[n][l];
				
			for(l=0; l<6; l++)
				Pr[l] += weightJ * Coefn[n][3+l];
		}/* end n */
			
	/* multiple by 4 pi */
		Er *= 4.0 * PI;
		for(l=0; l<3; l++)
			Fr[l] *= 4.0 * PI;
			
		for(l=0; l<6; l++)
			Pr[l] *= 4.0 * PI;
			
		/* Now we have Er, Fr and Pr for this cell */
		rho = pG->U[k+koff][j+joff][i+ioff].d;
			
		Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
		Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
		Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;
			
		/* Note that velocity used in the source term is guess vel */
			
		for(l=0; l<3; l++)
			Velguess[l] = pG->Velguess[k+koff][j+joff][i+ioff][l];
		
		/* The new gas temperature is in sol[j][i][nelements], the last elements */
		Sigma = pRG->R[k][j][i][ifr].Sigma[0];
		SigmaI = pRG->R[k][j][i][ifr].Sigma[1];
			
		/* Only absorption opacity cares about gas temperature */
			
		Tgas = sol[ifr][N];
		Tgas4 = SQR(Tgas) * SQR(Tgas);
			
			
		CoFr[0] = Fr[0] - (Velguess[0] * Er + Velguess[0] * Pr[0] + Velguess[1] * Pr[1] + Velguess[2] * Pr[3])/Crat;
		CoFr[1] = Fr[1] - (Velguess[1] * Er + Velguess[0] * Pr[1] + Velguess[1] * Pr[2] + Velguess[2] * Pr[4])/Crat;
		CoFr[2] = Fr[2] - (Velguess[2] * Er + Velguess[0] * Pr[3] + Velguess[1] * Pr[4] + Velguess[2] * Pr[5])/Crat;
			
		/* The energy equation is self-consistent with specific intensity only requires \sum n\dot v=0 */
		/* However, the momentum equation requies \sum 3nn=1  in order to be consistent with specific intensity */
		for(l=0; l<3; l++){
			Msource[l] =  dt * Prat * SigmaI * CoFr[l] - dt * Prat * Velguess[l] * (Sigma * Tgas4 - SigmaI * Er)/Crat; /* This is actually momentum source term */
			DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]); /* the kinetic energy change related to the momentum change */
					
					
			if(DeltaKin[l] + 0.5 * rho * Vel[l] * Vel[l] < 0.0){
                 Msource[l] = -rho * Vel[l];
                 DeltaKin[l] = 0.5 * rho * Vel[l] * Vel[l];
        	}
					
		}/* end l for three directions */
				
				
		/* Add the change of internal energy */
		Ersource = Tcoef[ifr][0] * (Tgas - pG->tgas[k+koff][j+joff][i+ioff]) + DeltaKin[0] + DeltaKin[1] + DeltaKin[2];
			
			
		/* Now put the energy and momentum source term back */	
		/* This is initialized to be zero before it is called */
		pG->Radheat[k+koff][j+joff][i+ioff] += (pRG->wnu[ifr] * Ersource);
				
		/* Radiation source term for gas pressure alone */
		pG->Pgsource[k+koff][j+joff][i+ioff] += (pRG->wnu[ifr] * Tcoef[ifr][0] * (Tgas - pG->tgas[k+koff][j+joff][i+ioff]) * Gamma_1);
			
		for(l=0; l<3; l++)
			pG->Frsource[k+koff][j+joff][i+ioff][l] += (pRG->wnu[ifr] * Msource[l]);
			
			
	}/* end ifr */
	
	return;
	
}


/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to scattering opacity */

/* moments of the radiation are already updated before entering this function */
/* So no need to calculate the moments again */

void RadSsource(RadGridS *pRG, GridS *pG)
{
	int i, j, k, l, nf, ifr;
	int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
	int koff = 0, joff = 0, ioff = 0;
	int nDim; 
	
	
	Real Sigma;
	Real Er, Ersource;
	Real Fr[3], Pr[6], Vel[3], Velguess[3], Msource[3], DeltaKin[3], CoFr[3]; /* velocity source terms and change of kinetic energy */
	Real dt = pG->dt;
	Real rho;
	
	
	nDim = 1;
	for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;
	
	ioff = nghost - Radghost;
	
	/* Including the ghost zones */
	
	if(nDim > 1){
		jl -= Radghost;
		ju += Radghost;
		
		joff = nghost - Radghost;
	}
	
	if(nDim > 2){
		kl -= Radghost;
		ku += Radghost;
		
		koff = nghost - Radghost;
	}
	
	nf = pRG->nf;
	

	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
				for(ifr=0; ifr<nf; ifr++){
			
			
					/* calculate Er and Fr used in the update */
					Er = 4.0 * PI * pRG->R[k][j][i][ifr].J;
					for(l=0; l<3; l++)
						Fr[l] = 4.0 * PI * pRG->R[k][j][i][ifr].H[l];
			
					for(l=0; l<6; l++)
						Pr[l] = 4.0 * PI * pRG->R[k][j][i][ifr].K[l];
			
						
					/* Now we have Er, Fr and Pr for this cell */
					rho = pG->U[k+koff][j+joff][i+ioff].d;
			
					Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
					Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
					Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;
				
					for(l=0; l<3; l++)
						Velguess[l] = pG->Velguess[k+koff][j+joff][i+ioff][l];
			
			
					/* The new gas temperature is in sol[j][i][nelements], the last elements */
					Sigma = pRG->R[k][j][i][ifr].Sigma[2];
			
					CoFr[0] = Fr[0] - (Velguess[0] * Er + Velguess[0] * Pr[0] + Velguess[1] * Pr[1] + Velguess[2] * Pr[3])/Crat;
					CoFr[1] = Fr[1] - (Velguess[1] * Er + Velguess[0] * Pr[1] + Velguess[1] * Pr[2] + Velguess[2] * Pr[4])/Crat;
					CoFr[2] = Fr[2] - (Velguess[2] * Er + Velguess[0] * Pr[3] + Velguess[1] * Pr[4] + Velguess[2] * Pr[5])/Crat;
			
					/* The energy equation is self-consistent with specific intensity only requires \sum n\dot v=0 */
					/* However, the momentum equation requies \sum 3nn=1  in order to be consistent with specific intensity */
					

					/* For scattering opacity, use kinetic energy change to update velocity momentum so that gas temperature is unchanged */
					for(l=0; l<3; l++){
				
						Msource[l] = dt * Prat * Sigma * CoFr[l]; /* This is actually momentum source term */
						DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]); 
						
						
						if(DeltaKin[l] + 0.5 * rho * Vel[l] * Vel[l] < 0.0){
                 			Msource[l] = -rho * Vel[l];
                 			DeltaKin[l] = 0.5 * rho * Vel[l] * Vel[l];
        				}
					
														
					}/* end l for three directions */
				
					Ersource = DeltaKin[0] + DeltaKin[1] + DeltaKin[2];
				
				
		
					/* Now put the energy and momentum source term back */	
					/* This is initialized to be zero before it is called */
					pG->Radheat[k+koff][j+joff][i+ioff] += (pRG->wnu[ifr] * Ersource);
			
					for(l=0; l<3; l++)
						pG->Frsource[k+koff][j+joff][i+ioff][l] += (pRG->wnu[ifr] * Msource[l]);
					
					
				}/* End ifr */
			
			}/* end i */
		}/* end j */
	}/* end k */
	
	
	
}

/* Function to post processing gas T and Er */
void ComptTEr(DomainS *pD)
{
    
    RadGridS *pRG=(pD->RadGrid);
	GridS *pG = (pD->Grid);

    int i, j, k, ig, jg, kg;
	int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
	int koff = 0, joff = 0, ioff = 0;
	int nDim;
    int nf = pRG->nf;
    int ifr;

    Real dt =pG->dt;
    Real Tgas, Er, Ernew, Tr, Ersource, sigmas;
    Real coefA, coefK, coefB,coef1,coef2,coef3,coef4;
    PrimS Wtemp;

	
	nDim = 1;
	for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;
	
	ioff = nghost - Radghost;
	
	/* Including the ghost zones */
	
	if(nDim > 1){
		jl -= Radghost;
		ju += Radghost;
		
		joff = nghost - Radghost;
	}
	
	if(nDim > 2){
		kl -= Radghost;
		ku += Radghost;
		
		koff = nghost - Radghost;
	}
    
    
   
    
    for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
                kg = k + koff;
                jg = j + joff;
                ig = i + ioff;
                
                Ersource = 0.0;
                pG->Ercompt[kg][jg][ig] = 0.0;
                Wtemp = Cons_to_Prim(&(pG->U[kg][jg][ig]));
                Tgas = Wtemp.P / (R_ideal * Wtemp.d);
                
				for(ifr=0; ifr<nf; ifr++){
                    
                    Er = 4.0 * PI * pRG->R[k][j][i][ifr].J;
                    sigmas = pRG->R[k][j][i][ifr].Sigma[2];
                    
                    Tr = sqrt(Er);
                    Tr = sqrt(Tr);
                    
                    coefA = 4.0 * dt * Crat * sigmas / (T_e/Tunit);
                    coefK = (Gamma - 1.0) * Prat/(R_ideal * Wtemp.d);
                    coefB = Tgas + coefK * Er;
                    coef1 = coefA * coefK;
                    coef2 = coefA;
                    coef3 = 1.0 - coefA * coefB;
                    coef4 = -Er;
                    
                    if(Tr < Tgas){
                        Tr = rtsafe(Tcompton, Tr * (1.0 - 0.01), Tgas * (1.0 + 0.01), 1.e-10, coef1, coef2, coef3, coef4);
                    }
                    else{
                        
                        Tr = rtsafe(Tcompton, Tgas * (1.0 - 0.01), Tr * (1.0 + 0.01), 1.e-10, coef1, coef2, coef3, coef4);
                    }
                    
                    Ernew = SQR(SQR(Tr));
                    
                    pG->Ercompt[kg][jg][ig] += (pRG->wnu[ifr] * Ernew);
                    Ersource += (pRG->wnu[ifr] * (Ernew - Er));
                    
                }/* End ifr */
                
                pG->Tcompt[kg][jg][ig] = (Wtemp.P - Prat * Ersource * Gamma_1)/(R_ideal * Wtemp.d);
                
            }/* end i */
        }/* End j */
    }/* End k */
    
    
}



/* Update the opacity for the whole grid, including the ghost zones */
void UpdateOpacity(DomainS *pD)
{
	GridS *pG=(pD->Grid);
	RadGridS *pRG=(pD->RadGrid);
	int i,j,k,ifr;
	int il = pRG->is, iu = pRG->ie;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
	int nf = pRG->nf;
	int ig,jg,kg,ioff,joff,koff;
	
	
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
	
	for (k=kl; k<=ku; k++) {
		kg = k + koff;
		for (j=jl; j<=ju; j++) {
			jg = j + joff;
			for (i=il; i<=iu; i++) {
				ig = i + ioff;
				
				/* ------------------------------------*/
				/* update the opacity */
				for(ifr=0; ifr<nf; ifr++) {
					get_full_opacity(pG,ifr,ig,jg,kg,&(pRG->R[k][j][i][ifr].Sigma[0]));
				}/* end ifr */
			}/* end i */
		}/* end j */
	}/* end k */
	
	
}



/* Function to get the estimated velocity to get better 
 * momentum conservation */
/* The equation we used to estimate the velocity is *
 * dFr/dt = -C(sigmas+sigmaa)* (Fr-(1+f)vEr/C)
 * drhov/dt=P(sigmas+sigmaa)*(Fr-(1+f)vEr/C)
 * Here f is the diagonal component 
 * Er is kept constant during this estimate *
 */

void GetVelguess(DomainS *pD)
{
	
	GridS *pG=(pD->Grid);
	RadGridS *pRG=(pD->RadGrid);
	
	Real Er, sigma, rho, dt;
	Real Fr[3], Pr[3], Vel[3], M0[3];
	int i, j, k, ioff, joff, koff, nDim, ifr;
	int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
	int n;
	
	
	nDim = 1;
	for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;
	
	ioff = nghost - Radghost;
	joff = 0;
	koff = 0;
	
	/* Including the ghost zones */
	
	if(nDim > 1){
		jl -= Radghost;
		ju += Radghost;
		
		joff = nghost - Radghost;
	}
	
	if(nDim > 2){
		kl -= Radghost;
		ku += Radghost;
		
		koff = nghost - Radghost;
	}
	
	/* estimated the velocity at half time step */
	dt = 0.5 * pG->dt;
	
	
	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
				/* calculate Er and Fr used in the update */
				
				/* first, calculate frequency weighted Er, Pr and Fr */
				
				Er = 0.0;
				for(n=0; n<3; n++)
					Fr[n] = 0.0;
				for(n=0; n<3; n++)
					Pr[n] = 0.0;
				
				sigma = 0.0;
								
				for(ifr=0; ifr<pRG->nf; ifr++){
					Er += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].J);
					for(n=0; n<3; n++)
						Fr[n] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].H[n]);
					
					
					Pr[0] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].K[0]);
					Pr[1] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].K[2]);
					Pr[2] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].K[5]);
					
					
					sigma += (pRG->wnu[ifr] * (pRG->R[k][j][i][ifr].Sigma[1]+pRG->R[k][j][i][ifr].Sigma[2]));
					
				}
				
				Er *= (4.0 * PI);
				for(n=0; n<3; n++)
					Fr[n] *= (4.0 * PI);
				for(n=0; n<3; n++)
					Pr[n] *= (4.0 * PI);
				
					/* Now we have Er, Fr and Pr for this cell */
				
				M0[0] = Prat * Fr[0] / Crat + pG->U[k+koff][j+joff][i+ioff].M1;
				M0[1] = Prat * Fr[1] / Crat + pG->U[k+koff][j+joff][i+ioff].M2;
				M0[2] = Prat * Fr[2] / Crat + pG->U[k+koff][j+joff][i+ioff].M3;
				
				rho = pG->U[k+koff][j+joff][i+ioff].d;				
				
				
				Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
				Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
				Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;
				
				for(n=0; n<3; n++)
					pG->Velguess[k+koff][j+joff][i+ioff][n] = (rho * Vel[n] + dt * sigma * M0[n] * Crat)/(rho + dt * sigma * Crat * rho + dt * Prat * sigma * (Er + Pr[n])/Crat);
				
			}
		}
	}
	
	
	
}




/* Function to calculate the reduction factor of speed of light due to opacity */

void GetSpeedfactor(DomainS *pD)
{
	
	GridS *pG=(pD->Grid);
	RadGridS *pRG=(pD->RadGrid);
	
	
	int i, j, k, ifr, nDim, n;
	int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
    
    Real  sigmas, sigmaa, alpha;
    Real dS[3];
    dS[0] = pRG->dx1;
    dS[1] = pRG->dx2;
    dS[2] = pRG->dx3;
    
#ifdef CYLINDRICAL
	const Real *r=pG->r;
    int offset = nghost - Radghost;
#endif
	
	
	nDim = 1;
	for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;
	
		
	/* Including the ghost zones */
	
	if(nDim > 1){
		jl -= Radghost;
		ju += Radghost;
		
    }
	
	if(nDim > 2){
		kl -= Radghost;
		ku += Radghost;

	}
	

	
	
	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
#ifdef CYLINDRICAL
				/* The scale factor r[i] is the same for each i, for different angles j */
                dS[1] = pRG->dx2 * r[i+offset];
#endif
				for(ifr=0; ifr<pRG->nf; ifr++){
                    
                    sigmas = pRG->R[k][j][i][ifr].Sigma[2];
					/* The absorption opacity in front of I */
					sigmaa = pRG->R[k][j][i][ifr].Sigma[1];
                    
                    /* for three direction */
                    for(n=0; n<nDim; n++){
                        ReduceVelocity(sigmaa+sigmas, dS[n], &alpha);
                        
                        pRG->Speedfactor[k][j][i][ifr][n] = alpha;
                    
                    }/* end Dim */
                    
                }/* end ifr */
				
			}/* end i */
		}/* end j */
	}/* end k */
	
	
	return;
}


/* sol takes the initial guess, calculated without the velocity dependent terms */
/* inisol takes the solution at time step n */
/* Md and RHS are pre-allocated memory */
/* They are used as temporary memory */
/* This is for the 3D case */


void Absorption(const int nf, const int N, Real **sol, Real **inisol, Real ***Ma, Real **Mdcoef, Real **Tcoef, Real *Md, Real *RHS, int *flag)
{
	

	int n;
	const int MAXIte = 15;
	int count;
	const Real TOL = 1.e-16;
	Real Tgas, Tgas3, Tgas4;
	Real residual;
	int line=N-1;
	int ifr;
	

	
	
	/* Now solve the non-linear set of equations for each cell */
	for(ifr=0; ifr<nf; ifr++){
		
			
		Tgas = sol[ifr][N];
		Tgas3 = SQR(Tgas) * Tgas;
		Tgas4 = Tgas * Tgas3;
				
		*flag = 1;
			
			
		/* First, calculate the RHS for the guess solution */				
		for(n=0; n<N-1; n++){
			/* set Md, Md is used to invert the matrix */
			Md[n] = 4.0 * Mdcoef[ifr][n] * Tgas3;
				
			RHS[n] = Mdcoef[ifr][n] * Tgas4 - (inisol[ifr][n] - inisol[ifr][N-1]);
			RHS[n] += (Ma[ifr][n][0] * sol[ifr][n] - Ma[ifr][N-1][0] * sol[ifr][N-1]);
				
		}/* end n */
		/* Now the line N-1 */
		/* Now line=N-1 */
		Md[line] = 4.0 * Mdcoef[ifr][line] * Tgas3;
			
		RHS[line] = Mdcoef[ifr][line] * Tgas4 + (Ma[ifr][line][0] + Ma[ifr][line][1]) * sol[ifr][line] - inisol[ifr][line];
		for(n=0; n<N-1; n++){
			RHS[line] += (Ma[ifr][n][1] * sol[ifr][n]);
				
		}
			
		/* Now the last line nelements */
		Md[N] = Tcoef[ifr][0] + 4.0 * Tcoef[ifr][1] * Tgas3;
		RHS[N] = Tcoef[ifr][0] * (Tgas - inisol[ifr][N]) + Tcoef[ifr][1] * Tgas4;
		for(n=0; n<N; n++){
			RHS[N] += (Ma[ifr][n][2] * sol[ifr][n]);
		}
			
		/* now calculate the norm of residual */
		residual = 0.0;
		for(n=0; n<=N; n++)
			residual += SQR(RHS[n]);
			
		residual /= (N+1);
			
		residual = sqrt(residual);
		count = 0;
			
		/****************************************************/
		/* Do the iteration */
		while((residual > TOL) && (count < MAXIte)){
			count++;
			/* calculate Inverse(Jacobi) * RHS */
			AbsorptionMatrix(N, Ma[ifr], Md, RHS);
				
			/* The result is stored in RHS */
			/* Update the guess solution */
			for(n=0; n<=N; n++){
				sol[ifr][n] -= RHS[n];
			}
				
		/*------------------------------------------------*/
			Tgas = sol[ifr][N];
			Tgas3 = SQR(Tgas) * Tgas;
			Tgas4 = Tgas * Tgas3;
				
				
		/* update RHS and Md */
			for(n=0; n<N-1; n++){
				/* set Md, Md is used to invert the matrix */
				Md[n] = 4.0 * Mdcoef[ifr][n] * Tgas3;
					
				RHS[n] = Mdcoef[ifr][n] * Tgas4 - (inisol[ifr][n] - inisol[ifr][N-1]);
				RHS[n] += (Ma[ifr][n][0] * sol[ifr][n] - Ma[ifr][N-1][0] * sol[ifr][N-1]);
					
			}/* end n */
		/* Now the line N-1 */
		/* Now line=N-1 */
			Md[line] = 4.0 * Mdcoef[ifr][line] * Tgas3;				
			RHS[line] = Mdcoef[ifr][line] * Tgas4 + (Ma[ifr][line][0] + Ma[ifr][line][1]) * sol[ifr][line] - inisol[ifr][line];
			for(n=0; n<N-1; n++){
				RHS[line] += (Ma[ifr][n][1] * sol[ifr][n]);
						
			}
				
			/* Now the last line nelements */
			Md[N] = Tcoef[ifr][0] + 4.0 * Tcoef[ifr][1] * Tgas3;
			RHS[N] = Tcoef[ifr][0] * (Tgas - inisol[ifr][N]) + Tcoef[ifr][1] * Tgas4;
			for(n=0; n<N; n++){
				RHS[N] += (Ma[ifr][n][2] * sol[ifr][n]);
			}
				
			/*------------------------------------------------*/
			/* update the Residual */
			/* now calculate the norm of residual */
			residual = 0.0;
			for(n=0; n<=N; n++)
				residual += SQR(RHS[n]);
				
			residual /= (N+1);
				
			residual = sqrt(residual);
				
		}
				
		  /* When matrix does not converge, do not add radiation source term */
		 if(residual != residual){
            for(n=0; n<=N; n++){
                sol[ifr][n] = inisol[ifr][n];
         	}
		 	*flag = 0;
        }
			
		/* Now the solution is stored in sol[j][i] */
		if(residual > TOL)
			printf("Final residual: %e Iterations: %d\n",residual,count);
			
			
			
	}/* end ifr */
	
	
}




/*********************************************************************************************************************
 * The function below is actually not used anymore **************
 **********************************************************************************************************************
 */


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

		if(coef1 < 1.e-18){
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
