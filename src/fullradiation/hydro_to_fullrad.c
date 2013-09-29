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

void RadAsource2D(const int ifr, const int N, RadGridS *pRG, GridS *pG, Real **Tcoef, Real ****Coefn, Real ****Coefnn, Real ***sol)
{
	int i, j, n, l;
	int is, ie, js, je;
	int offset;

	
	Real Sigma;
	Real Er, Ersource, weightJ;
	Real Fr[3], Pr[6], Vel[3], Velguess[3], Msource[3], DeltaKin[3], CoFr[3]; /* velocity source terms and change of kinetic energy */
	Real dt = pG->dt;
	Real Tgas, Tgas4;
	Real rho;
	
	offset = nghost - Radghost;
	
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	
	for(j=js; j<=je; j++){
		for(i=is; i<=ie; i++){
			/* calculate Er and Fr used in the update */
			Er = 0.0;
			for(l=0; l<3; l++)
				Fr[l] = 0.0;
			
			for(l=0; l<6; l++)
				Pr[l] = 0.0;

			/* loop through all l and n, nelements = nang * noct */
			/* weight for each ray is alerady included in sol[][][] */
			for(n=0; n<N; n++){
				weightJ = sol[j][i][n];
				Er += weightJ;
				
				for(l=0; l<3; l++)
					Fr[l] += weightJ * Coefn[j][i][n][l];
				
				for(l=0; l<6; l++)
					Pr[l] += weightJ * Coefnn[j][i][n][l];
			}
			
			/* multiple by 4 pi */
			Er *= 4.0 * PI;
			for(l=0; l<3; l++)
				Fr[l] *= 4.0 * PI;
			
			for(l=0; l<6; l++)
				Pr[l] *= 4.0 * PI;
			
			/* Now we have Er, Fr and Pr for this cell */
			rho = pG->U[0][j+offset][i+offset].d;
			
			Vel[0] = pG->U[0][j+offset][i+offset].M1 / rho;
			Vel[1] = pG->U[0][j+offset][i+offset].M2 / rho;
			Vel[2] = pG->U[0][j+offset][i+offset].M3 / rho;
			
			/* Note that velocity used in the source term is guess vel */
			
			for(l=0; l<3; l++)
				Velguess[l] = pG->Velguess[0][j+offset][i+offset][l];
			
		/*	for(l=0; l<3; l++)
				Kin[l] = 0.5 * rho * SQR(Vel[l]);
		*/	
			/* The new gas temperature is in sol[j][i][nelements], the last elements */
			Sigma = pRG->R[ifr][0][j][i].Sigma[1];
			
			/* Only absorption opacity cares about gas temperature */
			
			Tgas = sol[j][i][N];
			Tgas4 = SQR(Tgas) * SQR(Tgas);
			
			
			CoFr[0] = Fr[0] - (Velguess[0] * Er + Velguess[0] * Pr[0] + Velguess[1] * Pr[1] + Velguess[2] * Pr[3])/Crat;
			CoFr[1] = Fr[1] - (Velguess[1] * Er + Velguess[0] * Pr[1] + Velguess[1] * Pr[2] + Velguess[2] * Pr[4])/Crat;
			CoFr[2] = Fr[2] - (Velguess[2] * Er + Velguess[0] * Pr[3] + Velguess[1] * Pr[4] + Velguess[2] * Pr[5])/Crat;
			
			/* The energy equation is self-consistent with specific intensity only requires \sum n\dot v=0 */
			/* However, the momentum equation requies \sum 3nn=1  in order to be consistent with specific intensity */
		/*	for(l=0; l<3; l++)
				DeltaKin[l] = dt * Prat * Sigma * (Velguess[l] * CoFr[l] - SQR(Velguess[l]) * (Tgas4 - Er)/Crat);
		*/	

			for(l=0; l<3; l++){
				Msource[l] =  dt * Prat * Sigma * CoFr[l] - dt * Prat * Velguess[l] * Sigma * (Tgas4 - Er)/Crat; /* This is actually momentum source term */
				DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]); /* the kinetic energy change related to the momentum change */
				
				
			/*	Kin[l] += DeltaKin[l];
				
				if(Kin[l] >= 0.0){
					Msource[l] = sqrt(2.0*rho*Kin[l]) - rho * Vel[l];
					Ersource += DeltaKin[l];
				}else {
					
					Msource[l] =  dt * Prat * Sigma * CoFr[l] - dt * Prat * Velguess[l] * Sigma * (Tgas4 - Er)/Crat; 
				 
					Ersource += (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]);
				}
			 */
				
			}/* end l for three directions */
			
						
			/* Add the change of internal energy */
			Ersource = Tcoef[j][i] * (Tgas - pG->tgas[0][j+offset][i+offset]) + DeltaKin[0] + DeltaKin[1] + DeltaKin[2];
		
			
			/* Now put the energy and momentum source term back */	
			/* This is initialized to be zero before it is called */
			pG->Radheat[0][j+offset][i+offset] += (pRG->wnu[ifr] * Ersource);
			
			for(l=0; l<3; l++)
				pG->Frsource[0][j+offset][i+offset][l] += (pRG->wnu[ifr] * Msource[l]);
			
		}/* end i */
	}/* end j */
	
	
	
}



/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to absorption opacity */
/* This is for 3D case */

void RadAsource3D(const int ifr, const int N, RadGridS *pRG, GridS *pG, Real ***Tcoef, Real *****Coefn, Real *****Coefnn, Real ****sol)
{
	int i, j, k, n, l;
	int is, ie, js, je, ks, ke;
	int offset;
	
	
	Real Sigma;
	Real Er, Ersource, weightJ;
	Real Fr[3], Pr[6], Vel[3], Velguess[3], Msource[3], DeltaKin[3], CoFr[3]; /* velocity source terms and change of kinetic energy */
	Real dt = pG->dt;
	Real Tgas, Tgas4;
	Real rho;
	
	offset = nghost - Radghost;
	
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	ks = pRG->ks; ke = pRG->ke;
	
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){
			for(i=is; i<=ie; i++){
				/* calculate Er and Fr used in the update */
				Er = 0.0;
				for(l=0; l<3; l++)
					Fr[l] = 0.0;
			
				for(l=0; l<6; l++)
					Pr[l] = 0.0;
			
				/* loop through all l and n, nelements = nang * noct */
				/* weight for each ray is alerady included in sol[][][] */
				for(n=0; n<N; n++){
					weightJ = sol[k][j][i][n];
					Er += weightJ;
				
					for(l=0; l<3; l++)
						Fr[l] += weightJ * Coefn[k][j][i][n][l];
				
					for(l=0; l<6; l++)
						Pr[l] += weightJ * Coefnn[k][j][i][n][l];
				}
			
				/* multiple by 4 pi */
				Er *= 4.0 * PI;
				for(l=0; l<3; l++)
					Fr[l] *= 4.0 * PI;
			
				for(l=0; l<6; l++)
					Pr[l] *= 4.0 * PI;
			
				/* Now we have Er, Fr and Pr for this cell */
				rho = pG->U[k+offset][j+offset][i+offset].d;
			
				Vel[0] = pG->U[k+offset][j+offset][i+offset].M1 / rho;
				Vel[1] = pG->U[k+offset][j+offset][i+offset].M2 / rho;
				Vel[2] = pG->U[k+offset][j+offset][i+offset].M3 / rho;
			
				/* Note that velocity used in the source term is guess vel */
			
				for(l=0; l<3; l++)
					Velguess[l] = pG->Velguess[k+offset][j+offset][i+offset][l];
				
				/* The new gas temperature is in sol[j][i][nelements], the last elements */
				Sigma = pRG->R[ifr][k][j][i].Sigma[1];
			
				/* Only absorption opacity cares about gas temperature */
			
				Tgas = sol[k][j][i][N];
				Tgas4 = SQR(Tgas) * SQR(Tgas);
			
			
				CoFr[0] = Fr[0] - (Velguess[0] * Er + Velguess[0] * Pr[0] + Velguess[1] * Pr[1] + Velguess[2] * Pr[3])/Crat;
				CoFr[1] = Fr[1] - (Velguess[1] * Er + Velguess[0] * Pr[1] + Velguess[1] * Pr[2] + Velguess[2] * Pr[4])/Crat;
				CoFr[2] = Fr[2] - (Velguess[2] * Er + Velguess[0] * Pr[3] + Velguess[1] * Pr[4] + Velguess[2] * Pr[5])/Crat;
			
				/* The energy equation is self-consistent with specific intensity only requires \sum n\dot v=0 */
				/* However, the momentum equation requies \sum 3nn=1  in order to be consistent with specific intensity */
				for(l=0; l<3; l++){
					Msource[l] =  dt * Prat * Sigma * CoFr[l] - dt * Prat * Velguess[l] * Sigma * (Tgas4 - Er)/Crat; /* This is actually momentum source term */
					DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]); /* the kinetic energy change related to the momentum change */
					
					
					/*	Kin[l] += DeltaKin[l];
					 
					 if(Kin[l] >= 0.0){
					 Msource[l] = sqrt(2.0*rho*Kin[l]) - rho * Vel[l];
					 Ersource += DeltaKin[l];
					 }else {
					 
					 Msource[l] =  dt * Prat * Sigma * CoFr[l] - dt * Prat * Velguess[l] * Sigma * (Tgas4 - Er)/Crat; 
					 
					 Ersource += (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]);
					 }
					 */
					
				}/* end l for three directions */
				
				
				/* Add the change of internal energy */
				Ersource = Tcoef[k][j][i] * (Tgas - pG->tgas[k+offset][j+offset][i+offset]) + DeltaKin[0] + DeltaKin[1] + DeltaKin[2];
			
			
				/* Now put the energy and momentum source term back */	
				/* This is initialized to be zero before it is called */
				pG->Radheat[k+offset][j+offset][i+offset] += (pRG->wnu[ifr] * Ersource);
			
				for(l=0; l<3; l++)
					pG->Frsource[k+offset][j+offset][i+offset][l] += (pRG->wnu[ifr] * Msource[l]);
			
			}/* end i */
		}/* end j */
	}/* end k */
	
	
	
}


/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to scattering opacity */

/* moments of the radiation are already updated before entering this function */
/* So no need to calculate the moments again */

void RadSsource(const int ifr, RadGridS *pRG, GridS *pG)
{
	int i, j, k, l;
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
	

	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
				/* calculate Er and Fr used in the update */
				Er = 4.0 * PI * pRG->R[ifr][k][j][i].J;
				for(l=0; l<3; l++)
					Fr[l] = 4.0 * PI * pRG->R[ifr][k][j][i].H[l];
			
				for(l=0; l<6; l++)
					Pr[l] = 4.0 * PI * pRG->R[ifr][k][j][i].K[l];
			
						
				/* Now we have Er, Fr and Pr for this cell */
				rho = pG->U[k+koff][j+joff][i+ioff].d;
			
				Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
				Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
				Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;
				
				for(l=0; l<3; l++)
					Velguess[l] = pG->Velguess[k+koff][j+joff][i+ioff][l];
			
			
				/* The new gas temperature is in sol[j][i][nelements], the last elements */
				Sigma = pRG->R[ifr][k][j][i].Sigma[2];
			
				CoFr[0] = Fr[0] - (Velguess[0] * Er + Velguess[0] * Pr[0] + Velguess[1] * Pr[1] + Velguess[2] * Pr[3])/Crat;
				CoFr[1] = Fr[1] - (Velguess[1] * Er + Velguess[0] * Pr[1] + Velguess[1] * Pr[2] + Velguess[2] * Pr[4])/Crat;
				CoFr[2] = Fr[2] - (Velguess[2] * Er + Velguess[0] * Pr[3] + Velguess[1] * Pr[4] + Velguess[2] * Pr[5])/Crat;
			
				/* The energy equation is self-consistent with specific intensity only requires \sum n\dot v=0 */
				/* However, the momentum equation requies \sum 3nn=1  in order to be consistent with specific intensity */
					

				/* For scattering opacity, use kinetic energy change to update velocity momentum so that gas temperature is unchanged */
				for(l=0; l<3; l++){
				
					Msource[l] = dt * Prat * Sigma * CoFr[l]; /* This is actually momentum source term */
					DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]); 
					
														
				}/* end l for three directions */
				
				Ersource = DeltaKin[0] + DeltaKin[1] + DeltaKin[2];
				
				
		
				/* Now put the energy and momentum source term back */	
				/* This is initialized to be zero before it is called */
				pG->Radheat[k+koff][j+joff][i+ioff] += (pRG->wnu[ifr] * Ersource);
			
				for(l=0; l<3; l++)
					pG->Frsource[k+koff][j+joff][i+ioff][l] += (pRG->wnu[ifr] * Msource[l]);
			
			}/* end i */
		}/* end j */
	}/* end k */
	
	
	
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
					get_full_opacity(pG,ifr,ig,jg,kg,&(pRG->R[ifr][k][j][i].Sigma[0]));
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
					Er += (pRG->wnu[ifr] * pRG->R[ifr][k][j][i].J);
					for(n=0; n<3; n++)
						Fr[n] += (pRG->wnu[ifr] * pRG->R[ifr][k][j][i].H[n]);
					
					
					Pr[0] += (pRG->wnu[ifr] * pRG->R[ifr][k][j][i].K[0]);
					Pr[1] += (pRG->wnu[ifr] * pRG->R[ifr][k][j][i].K[2]);
					Pr[2] += (pRG->wnu[ifr] * pRG->R[ifr][k][j][i].K[5]);
					
					
					sigma += (pRG->wnu[ifr] * (pRG->R[ifr][k][j][i].Sigma[1]+pRG->R[ifr][k][j][i].Sigma[2]));
					
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



/* sol takes the initial guess, calculated without the velocity dependent terms */
/* inisol takes the solution at time step n */
/* Md and RHS are pre-allocated memory */

void Absorption2D(const int N, RadGridS *pRG, Real ***sol, Real ***inisol, Real ***Ma, Real ***Mb, Real ***Mc, Real ***Mdcoef, Real **Tcoef, Real **T4coef, Real *Md, Real *RHS)
{
	
	int i, is, ie;
	int j, js, je;
	int n;
	const int MAXIte = 15;
	int count;
	const Real TOL = 1.e-12;
	Real Tgas, Tgas3, Tgas4;
	Real residual;
	int line=N-1;
	
	is=pRG->is; ie=pRG->ie;
	js=pRG->js; je=pRG->je;
	
	
	/* Now solve the non-linear set of equations for each cell */
	for(j=js; j<=je; j++){
		for(i=is; i<=ie; i++){
			
			Tgas = sol[j][i][N];
			Tgas3 = SQR(Tgas) * Tgas;
			Tgas4 = Tgas * Tgas3;
			
			
			/* First, calculate the RHS for the guess solution */				
			for(n=0; n<N-1; n++){
				/* set Md, Md is used to invert the matrix */
				Md[n] = 4.0 * Mdcoef[j][i][n] * Tgas3;
				
				RHS[n] = Mdcoef[j][i][n] * Tgas4 - (inisol[j][i][n] - inisol[j][i][N-1]);
				RHS[n] += (Ma[j][i][n] * sol[j][i][n] - Ma[j][i][N-1] * sol[j][i][N-1]);					
				
			}/* end n */
			/* Now the line N-1 */
			/* Now line=N-1 */
			Md[line] = 4.0 * Mdcoef[j][i][line] * Tgas3;
			
			RHS[line] = Mdcoef[j][i][line] * Tgas4 + (Ma[j][i][line] + Mb[j][i][line]) * sol[j][i][line] - inisol[j][i][line];
			for(n=0; n<N-1; n++){
				RHS[line] += (Mb[j][i][n] * sol[j][i][n]);
				
			}
			
			/* Now the last line nelements */
			Md[N] = Tcoef[j][i] + 4.0 * T4coef[j][i] * Tgas3;
			RHS[N] = Tcoef[j][i] * (Tgas - inisol[j][i][N]) + T4coef[j][i] * Tgas4;
			for(n=0; n<N; n++){
				RHS[N] += (Mc[j][i][n] * sol[j][i][n]);
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
				AbsorptionMatrix(N, Ma[j][i], Mb[j][i], Mc[j][i], Md, RHS);
				
				/* The result is stored in RHS */
				/* Update the guess solution */
				for(n=0; n<=N; n++){
					sol[j][i][n] -= RHS[n];	
				}
				
				/*------------------------------------------------*/
				Tgas = sol[j][i][N];
				Tgas3 = SQR(Tgas) * Tgas;
				Tgas4 = Tgas * Tgas3;
				
				
				/* update RHS and Md */
				for(n=0; n<N-1; n++){
					/* set Md, Md is used to invert the matrix */
					Md[n] = 4.0 * Mdcoef[j][i][n] * Tgas3;
					
					RHS[n] = Mdcoef[j][i][n] * Tgas4 - (inisol[j][i][n] - inisol[j][i][N-1]);
					RHS[n] += (Ma[j][i][n] * sol[j][i][n] - Ma[j][i][N-1] * sol[j][i][N-1]);					
					
				}/* end n */
				/* Now the line N-1 */
				/* Now line=N-1 */
				Md[line] = 4.0 * Mdcoef[j][i][line] * Tgas3;				
				RHS[line] = Mdcoef[j][i][line] * Tgas4 + (Ma[j][i][line] + Mb[j][i][line]) * sol[j][i][line] - inisol[j][i][line];
				for(n=0; n<N-1; n++){
					RHS[line] += (Mb[j][i][n] * sol[j][i][n]);
					
				}
				
				/* Now the last line nelements */
				Md[N] = Tcoef[j][i] + 4.0 * T4coef[j][i] * Tgas3;
				RHS[N] = Tcoef[j][i] * (Tgas - inisol[j][i][N]) + T4coef[j][i] * Tgas4;
				for(n=0; n<N; n++){
					RHS[N] += (Mc[j][i][n] * sol[j][i][n]);
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
			
			/* Now the solution is stored in sol[j][i] */
			if(residual > TOL)
				printf("Final residual: %e Iterations: %d\n",residual,count);
			
			
		}/* end i */
	}/* end j */
	
	
}





/* sol takes the initial guess, calculated without the velocity dependent terms */
/* inisol takes the solution at time step n */
/* Md and RHS are pre-allocated memory */
/* This is for the 3D case */

void Absorption3D(const int N, RadGridS *pRG, Real ****sol, Real ****inisol, Real ****Ma, Real ****Mb, Real ****Mc, Real ****Mdcoef, Real ***Tcoef, Real ***T4coef, Real *Md, Real *RHS)
{
	
	int i, is, ie;
	int j, js, je;
	int k, ks, ke;
	int n;
	const int MAXIte = 15;
	int count;
	const Real TOL = 1.e-10;
	Real Tgas, Tgas3, Tgas4;
	Real residual;
	int line=N-1;
	
	is=pRG->is; ie=pRG->ie;
	js=pRG->js; je=pRG->je;
	ks=pRG->ks; ke=pRG->ke;
	
	
	/* Now solve the non-linear set of equations for each cell */
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){
			for(i=is; i<=ie; i++){
			
				Tgas = sol[k][j][i][N];
				Tgas3 = SQR(Tgas) * Tgas;
				Tgas4 = Tgas * Tgas3;
			
			
				/* First, calculate the RHS for the guess solution */				
				for(n=0; n<N-1; n++){
					/* set Md, Md is used to invert the matrix */
					Md[n] = 4.0 * Mdcoef[k][j][i][n] * Tgas3;
				
					RHS[n] = Mdcoef[k][j][i][n] * Tgas4 - (inisol[k][j][i][n] - inisol[k][j][i][N-1]);
					RHS[n] += (Ma[k][j][i][n] * sol[k][j][i][n] - Ma[k][j][i][N-1] * sol[k][j][i][N-1]);					
				
				}/* end n */
				/* Now the line N-1 */
				/* Now line=N-1 */
				Md[line] = 4.0 * Mdcoef[k][j][i][line] * Tgas3;
			
				RHS[line] = Mdcoef[k][j][i][line] * Tgas4 + (Ma[k][j][i][line] + Mb[k][j][i][line]) * sol[k][j][i][line] - inisol[k][j][i][line];
				for(n=0; n<N-1; n++){
					RHS[line] += (Mb[k][j][i][n] * sol[k][j][i][n]);
				
				}
			
				/* Now the last line nelements */
				Md[N] = Tcoef[k][j][i] + 4.0 * T4coef[k][j][i] * Tgas3;
				RHS[N] = Tcoef[k][j][i] * (Tgas - inisol[k][j][i][N]) + T4coef[k][j][i] * Tgas4;
				for(n=0; n<N; n++){
					RHS[N] += (Mc[k][j][i][n] * sol[k][j][i][n]);
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
					AbsorptionMatrix(N, Ma[k][j][i], Mb[k][j][i], Mc[k][j][i], Md, RHS);
				
					/* The result is stored in RHS */
					/* Update the guess solution */
					for(n=0; n<=N; n++){
						sol[k][j][i][n] -= RHS[n];	
					}
				
					/*------------------------------------------------*/
					Tgas = sol[k][j][i][N];
					Tgas3 = SQR(Tgas) * Tgas;
					Tgas4 = Tgas * Tgas3;
				
				
					/* update RHS and Md */
					for(n=0; n<N-1; n++){
						/* set Md, Md is used to invert the matrix */
						Md[n] = 4.0 * Mdcoef[k][j][i][n] * Tgas3;
					
						RHS[n] = Mdcoef[k][j][i][n] * Tgas4 - (inisol[k][j][i][n] - inisol[k][j][i][N-1]);
						RHS[n] += (Ma[k][j][i][n] * sol[k][j][i][n] - Ma[k][j][i][N-1] * sol[k][j][i][N-1]);					
					
					}/* end n */
					/* Now the line N-1 */
					/* Now line=N-1 */
					Md[line] = 4.0 * Mdcoef[k][j][i][line] * Tgas3;				
					RHS[line] = Mdcoef[k][j][i][line] * Tgas4 + (Ma[k][j][i][line] + Mb[k][j][i][line]) * sol[k][j][i][line] - inisol[k][j][i][line];
					for(n=0; n<N-1; n++){
						RHS[line] += (Mb[k][j][i][n] * sol[k][j][i][n]);
						
					}
				
					/* Now the last line nelements */
					Md[N] = Tcoef[k][j][i] + 4.0 * T4coef[k][j][i] * Tgas3;
					RHS[N] = Tcoef[k][j][i] * (Tgas - inisol[k][j][i][N]) + T4coef[k][j][i] * Tgas4;
					for(n=0; n<N; n++){
						RHS[N] += (Mc[k][j][i][n] * sol[k][j][i][n]);
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
			
				/* Now the solution is stored in sol[j][i] */
				if(residual > TOL)
					printf("Final residual: %e Iterations: %d\n",residual,count);
			
			
			}/* end i */
		}/* end j */
	}/* end k */
	
	
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
