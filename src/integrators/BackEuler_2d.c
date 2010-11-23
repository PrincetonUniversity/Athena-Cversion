#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Use backward Euler method to update the radiation quantities
 * First set up the matrix
 * Then solve the matrix equations.
 * We need the flag for boundary condition.
 *
 * Backward Euler should be used for the whole mesh
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   BackEuler_2d()
 *
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






static int *IEuler = NULL;
static int *JEuler = NULL;
/*Store the index of all non-zero elements */
/* The matrix coefficient, this is oneD in GMRES solver */
static Real *Euler = NULL;


/* Right hand side of the Matrix equation */
static Real *RHSEuler = NULL;
/* RHSEuler is 0-------3*Nx*Ny-1; only for GMRES method */ 
static Real *INIguess = NULL;
/* Used for initial guess solution, also return the correct solution */







/********Public function****************/
/*-------BackEuler_2d(): Use back euler method to update E_r and Fluxr-----------*/



void BackEuler_2d(MeshS *pM)
{
/* Right now, only work for one domain. Modified later for SMR */


  	GridS *pG=(pM->Domain[0][0].Grid);
	Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
	Real dtodx2 = pG->dt/pG->dx2, hdtodx2 = 0.5 * pG->dt/pG->dx2;
	Real dt = pG->dt, dx = pG->dx1, dy = pG->dx2;
	int is = pG->is, ie = pG->ie;
  	int i, j, m, n, NoEr, NoFr1, NoFr2;
	int js = pG->js, je = pG->je;
	int ks = pG->ks;
	int Nx, Ny, Nmatrix, NZ_NUM;
	
	/* NZ_NUM is the number of non-zero elements in Matrix. It may change if periodic boundary condition is applied */

	Real SEE, SErho, SEm;
	Real temperature, velocity_x, velocity_y, pressure, Sigmas;

	Real Ci0, Ci1, Cj0, Cj1;
	/* This is equivilent to Cspeeds[] in 1D */

	Real temp1, temp2;
  	Real theta[11];
  	Real phi[11];
	Real psi[11];
  	Real Sigma_s, Sigma_t, Sigma_a;


	/* Boundary condition flag */
	int ix1, ox1, ix2, ox2, ix3, ox3;
	ix1 = pM->BCFlag_ix1;
	ox1 = pM->BCFlag_ox1;
	ix2 = pM->BCFlag_ix2;
	ox2 = pM->BCFlag_ox2;
	ix3 = pM->BCFlag_ix3;
	ox3 = pM->BCFlag_ox3;
	


/* Allocate memory space for the Matrix calculation, just used for this grids */
/* Nmatrix is the number of active cells just in this grids */
/* Matrix size should be 3 * Nmatrix, ghost zones are included*/
	Nx = ie - is + 1;
	Ny = je - js + 1;
   	Nmatrix = Ny * Nx;
	
	NZ_NUM = 31*Nmatrix;
	
	rad_hydro_init_2d(Nmatrix);




/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */


/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1 =size1, otherwise it is wrong */

/*
	for(i=is; i<=ie+1; i++){
	 	Cspeeds[i-is] = (U1d[i].Edd_11 - U1d[i-1].Edd_11) 
				/ (U1d[i].Edd_11 + U1d[i-1].Edd_11); 		
	}
*/
	for(j=js; j<=je; j++) {
		for(i=is; i<=ie; i++){
/* E is the total energy. should subtract the kinetic energy and magnetic energy density */
    		pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1 
			+ pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2)/pG->U[ks][j][i].d) * (Gamma - 1.0);
/* if MHD - 0.5 * Bx * Bx   */
    		temperature = pressure / (pG->U[ks][j][i].d * R_ideal);
		/* RHSEuler[0...N-1]  */
		Sigma_a = pG->U[ks][j][i].Sigma_a;
		
    		RHSEuler[3*(j-js)*Nx + 3*(i-is)]   = pG->U[ks][j][i].Er + Crat * dt * Sigma_a 
				* temperature * temperature * temperature * temperature;
    		RHSEuler[3*(j-js)*Nx + 3*(i-is) + 1] = pG->U[ks][j][i].Fr1 + dt * Sigma_a
				* temperature * temperature * temperature * temperature * pG->U[ks][j][i].M1 / pG->U[ks][j][i].d;
		RHSEuler[3*(j-js)*Nx + 3*(i-is) + 2] = pG->U[ks][j][i].Fr2 + dt * Sigma_a
				* temperature * temperature * temperature * temperature * pG->U[ks][j][i].M2 / pG->U[ks][j][i].d;

		/* For inflow boundary condition along x direction*/
		if((i == is) && (ix1 == 3)) {
			Ci0 = (sqrt(pG->U[ks][j][i].Edd_11) - sqrt(pG->U[ks][j][i-1].Edd_11)) 
				/ (sqrt(pG->U[ks][j][i].Edd_11) + sqrt(pG->U[ks][j][i-1].Edd_11));
			
			theta[2] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pG->U[ks][j][i-1].Edd_11);
			theta[3] = -Crat * hdtodx1 * (1.0 + Ci0);
			phi[2]	= theta[2] * sqrt(pG->U[ks][j][i-1].Edd_11);
			phi[3]	= theta[3] * sqrt(pG->U[ks][j][i-1].Edd_11);
			psi[2] = -Crat * hdtodx1 * (1.0 + Ci0) * pG->U[ks][j][i-1].Edd_21;
			psi[3] = phi[3];
			RHSEuler[3*(j-js)*Nx + 3*(i-is)] -= (theta[2] * pG->U[ks][j][i-1].Er + theta[3] * pG->U[ks][j][i-1].Fr1);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 1] -= (phi[2] * pG->U[ks][j][i-1].Er + phi[3] * pG->U[ks][j][i-1].Fr1);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 2] -= (psi[2] * pG->U[ks][j][i-1].Er + psi[3] * pG->U[ks][j][i-1].Fr2);
			
		}

		if((i == ie) && (ox1 == 3)) {
			Ci1 =  (sqrt(pG->U[ks][j][i+1].Edd_11) - sqrt(pG->U[ks][j][i].Edd_11)) 
				/ (sqrt(pG->U[ks][j][i+1].Edd_11) + sqrt(pG->U[ks][j][i].Edd_11));

			theta[7] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pG->U[ks][j][i+1].Edd_11);
			theta[8] = Crat * hdtodx1 * (1.0 - Ci1);
			phi[6]	= -theta[7] * sqrt(pG->U[ks][j][i+1].Edd_11);
			phi[7]	= -theta[8] * sqrt(pG->U[ks][j][i+1].Edd_11);
			psi[6]	= Crat * hdtodx1 * (1.0 - Ci1) * pG->U[ks][j][i+1].Edd_21;
			psi[7]	= phi[7];

			RHSEuler[3*(j-js)*Nx + 3*(i-is)] -= (theta[7] * pG->U[ks][j][i+1].Er + theta[8] * pG->U[ks][j][i+1].Fr1);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 1] -= (phi[6] * pG->U[ks][j][i+1].Er + phi[7] * pG->U[ks][j][i+1].Fr1);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 2] -= (psi[6] * pG->U[ks][j][i+1].Er + psi[7] * pG->U[ks][j][i+1].Fr2);
					
		}
		

		/* For inflow boundary condition along y direction*/
		if((j == js) && (ix2 == 3)) {
			Cj0 = (sqrt(pG->U[ks][j][i].Edd_22) - sqrt(pG->U[ks][j-1][i].Edd_22)) 
				/ (sqrt(pG->U[ks][j][i].Edd_22) + sqrt(pG->U[ks][j-1][i].Edd_22));
			
			theta[0] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pG->U[ks][j-1][i].Edd_22);
			theta[1] = -Crat * hdtodx2 * (1.0 + Cj0);
			phi[0]	= -Crat * hdtodx2 * (1.0 + Cj0) * pG->U[ks][j-1][i].Edd_21;
			phi[1]	= theta[1] * sqrt(pG->U[ks][j-1][i].Edd_22);
			psi[0] = theta[0] * sqrt(pG->U[ks][j-1][i].Edd_22);
			psi[1] = phi[1];

			RHSEuler[3*(j-js)*Nx] -= (theta[0] * pG->U[ks][j-1][i].Er + theta[1] * pG->U[ks][j-1][i].Fr2);
			RHSEuler[3*(j-js)*Nx + 1] -= (phi[0] * pG->U[ks][j-1][i].Er + phi[1] * pG->U[ks][j-1][i].Fr1);
			RHSEuler[3*(j-js)*Nx + 2] -= (psi[0] * pG->U[ks][j-1][i].Er + psi[1] * pG->U[ks][j-1][i].Fr2);
				
		}

		if((j == je) && (ox2 == 3)) {
			Cj1 =  (sqrt(pG->U[ks][j+1][i].Edd_22) - sqrt(pG->U[ks][j][i].Edd_22)) 
				/ (sqrt(pG->U[ks][j+1][i].Edd_22) + sqrt(pG->U[ks][j][i].Edd_22));

			theta[9] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pG->U[ks][j+1][i].Edd_22);
			theta[10] = Crat * hdtodx2 * (1.0 - Cj1);
			phi[8]	= theta[10] * pG->U[ks][j+1][i].Edd_21;
			phi[9]	= -theta[10] * sqrt(pG->U[ks][j+1][i].Edd_22);
			psi[8]	= -theta[9] * sqrt(pG->U[ks][j+1][i].Edd_22);
			psi[9]	= phi[9];

			RHSEuler[3*(j-js)*Nx + 3*(i-is)] -= (theta[9] * pG->U[ks][j+1][i].Er + theta[10] * pG->U[ks][j+1][i].Fr2);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 1] -= (phi[8] * pG->U[ks][j+1][i].Er + phi[9] * pG->U[ks][j+1][i].Fr1);
			RHSEuler[3*(j-js)*Nx + 3*(i-is) + 2] -= (psi[8] * pG->U[ks][j+1][i].Er + psi[9] * pG->U[ks][j+1][i].Fr2);
			
		}
				
		}
	}

	
/*--------------------Note--------------------*.


/* Step 1b: Setup the Matrix */
		
 	/* First, setup the guess solution. Guess solution is the solution from last time step */
	for(j=js; j<=je; j++){
		for(i=is; i<=ie; i++){
			INIguess[3*(j-js)*Nx + 3*(i-is)] = pG->U[ks][j][i].Er;
			INIguess[3*(j-js)*Nx + 3*(i-is)+1] = pG->U[ks][j][i].Fr1;
			INIguess[3*(j-js)*Nx + 3*(i-is)+2] = pG->U[ks][j][i].Fr2;
		}
	}	



	/*--------Now set the Euler matrix-----------*/
	for(j=js; j<=je; j++){
		for(i=is; i<=ie; i++){
			velocity_x = pG->U[ks][j][i].M1 / pG->U[ks][j][i].d;
			velocity_y = pG->U[ks][j][i].M2 / pG->U[ks][j][i].d;
			Sigma_a = pG->U[ks][j][i].Sigma_a;
			Sigma_t = pG->U[ks][j][i].Sigma_t;
			Sigma_s = Sigma_t - Sigma_a;
			Ci0 = (sqrt(pG->U[ks][j][i].Edd_11) - sqrt(pG->U[ks][j][i-1].Edd_11)) 
				/ (sqrt(pG->U[ks][j][i].Edd_11) + sqrt(pG->U[ks][j][i-1].Edd_11));
			Ci1 =  (sqrt(pG->U[ks][j][i+1].Edd_11) - sqrt(pG->U[ks][j][i].Edd_11)) 
				/ (sqrt(pG->U[ks][j][i+1].Edd_11) + sqrt(pG->U[ks][j][i].Edd_11));
			Cj0 = (sqrt(pG->U[ks][j][i].Edd_22) - sqrt(pG->U[ks][j-1][i].Edd_22)) 
				/ (sqrt(pG->U[ks][j][i].Edd_22) + sqrt(pG->U[ks][j-1][i].Edd_22));
			Cj1 =  (sqrt(pG->U[ks][j+1][i].Edd_22) - sqrt(pG->U[ks][j][i].Edd_22)) 
				/ (sqrt(pG->U[ks][j+1][i].Edd_22) + sqrt(pG->U[ks][j][i].Edd_22));
			theta[0] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pG->U[ks][j-1][i].Edd_22);
			theta[1] = -Crat * hdtodx2 * (1.0 + Cj0);
			theta[2] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pG->U[ks][j][i-1].Edd_11);
			theta[3] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[4] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pG->U[ks][j][i].Edd_11) 
				+ Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pG->U[ks][j][i].Edd_22)
				+ Crat * pG->dt * Sigma_a 
				+ pG->dt * (Sigma_a - Sigma_s) * ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x 
				+ velocity_y * pG->U[ks][j][i].Edd_21) * velocity_x / Crat
				+ pG->dt * (Sigma_a - Sigma_s) * ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y 
				+ velocity_x * pG->U[ks][j][i].Edd_21) * velocity_y / Crat;
			theta[5] = Crat * hdtodx1 * (Ci0 + Ci1)	- pG->dt * (Sigma_a - Sigma_s) * velocity_x;
			theta[6] = Crat * hdtodx2 * (Cj0 + Cj1)	- pG->dt * (Sigma_a - Sigma_s) * velocity_y;
			theta[7] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pG->U[ks][j][i+1].Edd_11);
			theta[8] = Crat * hdtodx1 * (1.0 - Ci1);
			theta[9] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pG->U[ks][j+1][i].Edd_22);
			theta[10] = Crat * hdtodx2 * (1.0 - Cj1);
			

			phi[0] = theta[1] * pG->U[ks][j-1][i].Edd_21;
			phi[1] = theta[0];
			phi[2] = theta[3] * pG->U[ks][j][i-1].Edd_11;
			phi[3] = theta[2]; 
			phi[4] = Crat * hdtodx1 * (Ci0 + Ci1) * pG->U[ks][j][i].Edd_11
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pG->U[ks][j][i].Edd_21   
			       - pG->dt * Sigma_t * ((1.0 + pG->U[ks][j][i].Edd_11) * velocity_x + pG->U[ks][j][i].Edd_21 * velocity_y) 
			       + pG->dt * Sigma_a * velocity_x;
			phi[5] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pG->U[ks][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pG->U[ks][j][i].Edd_22) 
				     + Crat * pG->dt * Sigma_t;
			phi[6] = theta[8] * pG->U[ks][j][i+1].Edd_11;
			phi[7] = theta[7];
			phi[8] = theta[10] * pG->U[ks][j+1][i].Edd_21;
			phi[9] = theta[9];


			psi[0] = theta[1] * pG->U[ks][j-1][i].Edd_22;
			psi[1] = theta[0];
			psi[2] = theta[3] * pG->U[ks][j][i-1].Edd_21;
			psi[3] = theta[2];			
			psi[4] = Crat * hdtodx1 * (Ci0 + Ci1) * pG->U[ks][j][i].Edd_21
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pG->U[ks][j][i].Edd_22   
			       - pG->dt * Sigma_t * ((1.0 + pG->U[ks][j][i].Edd_22) * velocity_y + pG->U[ks][j][i].Edd_21 * velocity_x) 
			       + pG->dt * Sigma_a * velocity_y;
			psi[5] = phi[5];
			psi[6] = theta[8] * pG->U[ks][j][i+1].Edd_21;
			psi[7] = theta[7];
			psi[8] = theta[10] * pG->U[ks][j+1][i].Edd_22;
			psi[9] = theta[9];

			NoEr = 31*((j-js)*Nx + (i-is)) + 4;
			NoFr1 = 31*((j-js)*Nx + (i-is)) + 15;
			NoFr2 = 31*((j-js)*Nx + (i-is)) + 25;

			if(i == is){
				/* Common elements for different boundary conditions */				

				Euler[NoEr+3]  = theta[7];				
				Euler[NoEr+4]  = theta[8];
				
				Euler[NoFr1+2]  = phi[6];				
				Euler[NoFr1+3]  = phi[7];
				
				
				Euler[NoFr2+2]  = psi[6];			
				Euler[NoFr2+3]  = psi[7];
				

				IEuler[NoEr]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoEr]   = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoEr+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				JEuler[NoEr+1] = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoEr+2] = 3*(j-js)*Nx + 3*(i-is) + 2;
				JEuler[NoEr+2] = 3*(j-js)*Nx + 3*(i-is);
				IEuler[NoEr+3] = 3*(j-js)*Nx + 3*(i-is) + 3;
				JEuler[NoEr+3] = 3*(j-js)*Nx + 3*(i-is);
				IEuler[NoEr+4] = 3*(j-js)*Nx + 3*(i-is) + 4;
				JEuler[NoEr+4] = 3*(j-js)*Nx + 3*(i-is);
								
				IEuler[NoFr1]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoFr1]   = 3*(j-js)*Nx + 3*(i-is) + 1;				
				IEuler[NoFr1+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				IEuler[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is) + 3;
				JEuler[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is) + 1;
				IEuler[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is) + 4;
				JEuler[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is) + 1;
							
				IEuler[NoFr2]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoFr2]   = 3*(j-js)*Nx + 3*(i-is) + 2;				
				IEuler[NoFr2+1] = 3*(j-js)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2+1] = 3*(j-js)*Nx + 3*(i-is) + 2;
				IEuler[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is) + 3;
				JEuler[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is) + 2;
				IEuler[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is) + 5;
				JEuler[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is) + 2;
				



				if(ix1 != 4){
				/* Nonperiodic boundary condition */
				
					Euler[NoEr-1]  = 0.0;
					IEuler[NoEr-1]   = 0;
					JEuler[NoEr-1]   = 0;
					Euler[NoEr-2]  = 0.0;
					IEuler[NoEr-2]   = 0;
					JEuler[NoEr-2]   = 0;

					Euler[NoFr1-1]  = 0.0;
					IEuler[NoFr1-1]   = 0;
					JEuler[NoFr1-1]   = 0;
					Euler[NoFr1-2]  = 0.0;
					IEuler[NoFr1-2]   = 0;
					JEuler[NoFr1-2]   = 0;

					Euler[NoFr2-1]  = 0.0;
					IEuler[NoFr2-1]   = 0;
					JEuler[NoFr2-1]   = 0;
					Euler[NoFr2-2]  = 0.0;
					IEuler[NoFr2-2]   = 0;
					JEuler[NoFr2-2]   = 0;

;				if(ix1 == 1){				
					Euler[NoEr]    = theta[2] + theta[4];
					Euler[NoEr+1]  = theta[5] - theta[3];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[2] + phi[4];
					Euler[NoFr1+1]  = phi[5] - phi[3];

					Euler[NoFr2]    = psi[2] + psi[4];
					Euler[NoFr2+1]  = psi[5] + psi[3];
				}
				else if(ix1 == 2){
					Euler[NoEr]    = theta[2] + theta[4];
					Euler[NoEr+1]  = theta[5] + theta[3];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[2] + phi[4];
					Euler[NoFr1+1]  = phi[5] + phi[3];

					Euler[NoFr2]    = psi[2] + psi[4];
					Euler[NoFr2+1]  = psi[5] + psi[3];
				}
				else if(ix1 == 3){
					Euler[NoEr]    = theta[4];
					Euler[NoEr+1]  = theta[5];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4];
					Euler[NoFr1+1]  = phi[5];

					Euler[NoFr2]    = psi[4];
					Euler[NoFr2+1]  = psi[5];
				}
				else 
					{goto on_error;}
				}
				else{
				/* For Periodic boundary condition ix1 == 4 */
					Euler[NoEr-1]  = theta[3];
					IEuler[NoEr-1]   = 3*(j-js)*Nx + 3*(ie-is) + 1;
					JEuler[NoEr-1]   = 3*(j-js)*Nx + 3*(i-is);
					Euler[NoEr-2]  = theta[2];
					IEuler[NoEr-2]   = 3*(j-js)*Nx + 3*(ie-is);
					JEuler[NoEr-2]   = 3*(j-js)*Nx + 3*(i-is);

					Euler[NoFr1-1]  = phi[3];
					IEuler[NoFr1-1]   =  3*(j-js)*Nx + 3*(ie-is) + 1;
					JEuler[NoFr1-1]   = 3*(j-js)*Nx + 3*(i-is) + 1;
					Euler[NoFr1-2]  = phi[2];
					IEuler[NoFr1-2]   = 3*(j-js)*Nx + 3*(ie-is);
					JEuler[NoFr1-2]   = 3*(j-js)*Nx + 3*(i-is) + 1;

					Euler[NoFr2-1]  = psi[3];
					IEuler[NoFr2-1]   = 3*(j-js)*Nx + 3*(ie-is) + 2;
					JEuler[NoFr2-1]   = 3*(j-js)*Nx + 3*(i-is) + 2;
					Euler[NoFr2-2]  = psi[2];
					IEuler[NoFr2-2]   = 3*(j-js)*Nx + 3*(ie-is);
					JEuler[NoFr2-2]   = 3*(j-js)*Nx + 3*(i-is) + 2;

					Euler[NoEr]    = theta[4];
					Euler[NoEr+1]  = theta[5];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4];
					Euler[NoFr1+1]  = phi[5];

					Euler[NoFr2]    = psi[4];
					Euler[NoFr2+1]  = psi[5];
				}
				
			}/* End i == is */
			else if(i == ie){
				/* Common elements for different boundary conditions */				

				Euler[NoEr-2]  = theta[2];				
				Euler[NoEr-1]  = theta[3];
				
				Euler[NoFr1-2]  = phi[2];				
				Euler[NoFr1-1]  = phi[3];
				
				
				Euler[NoFr2-2]  = psi[2];			
				Euler[NoFr2-1]  = psi[3];
				
				IEuler[NoEr-2] = 3*(j-js)*Nx + 3*(i-is) - 3;
				JEuler[NoEr-2] = 3*(j-js)*Nx + 3*(i-is);
				IEuler[NoEr-1] = 3*(j-js)*Nx + 3*(i-is) - 2;
				JEuler[NoEr-1] = 3*(j-js)*Nx + 3*(i-is);
				IEuler[NoEr]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoEr]   = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoEr+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				JEuler[NoEr+1] = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoEr+2] = 3*(j-js)*Nx + 3*(i-is) + 2;
				JEuler[NoEr+2] = 3*(j-js)*Nx + 3*(i-is);
				
				IEuler[NoFr1-2] = 3*(j-js)*Nx + 3*(i-is) - 3;
				JEuler[NoFr1-2] = 3*(j-js)*Nx + 3*(i-is) + 1;
				IEuler[NoFr1-1] = 3*(j-js)*Nx + 3*(i-is) - 2;
				JEuler[NoFr1-1] = 3*(j-js)*Nx + 3*(i-is) + 1;				
				IEuler[NoFr1]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoFr1]   = 3*(j-js)*Nx + 3*(i-is) + 1;				
				IEuler[NoFr1+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
				IEuler[NoFr2-2] = 3*(j-js)*Nx + 3*(i-is) - 3;
				JEuler[NoFr2-2] = 3*(j-js)*Nx + 3*(i-is) + 2;
				IEuler[NoFr2-1] = 3*(j-js)*Nx + 3*(i-is) - 1;
				JEuler[NoFr2-1] = 3*(j-js)*Nx + 3*(i-is) + 2;			
				IEuler[NoFr2]   = 3*(j-js)*Nx + 3*(i-is);
				JEuler[NoFr2]   = 3*(j-js)*Nx + 3*(i-is) + 2;				
				IEuler[NoFr2+1] = 3*(j-js)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2+1] = 3*(j-js)*Nx + 3*(i-is) + 2;
				
				



				if(ox1 != 4){
				/* Nonperiodic boundary condition */
				
					Euler[NoEr+3]  = 0.0;
					IEuler[NoEr+3]   = 0;
					JEuler[NoEr+3]   = 0;
					Euler[NoEr+4]  = 0.0;
					IEuler[NoEr+4]   = 0;
					JEuler[NoEr+4]   = 0;

					Euler[NoFr1+2]  = 0.0;
					IEuler[NoFr1+2]   = 0;
					JEuler[NoFr1+2]   = 0;
					Euler[NoFr1+3]  = 0.0;
					IEuler[NoFr1+3]   = 0;
					JEuler[NoFr1+3]   = 0;

					Euler[NoFr2+2]  = 0.0;
					IEuler[NoFr2+2]   = 0;
					JEuler[NoFr2+2]   = 0;
					Euler[NoFr2+3]  = 0.0;
					IEuler[NoFr2+3]   = 0;
					JEuler[NoFr2+3]   = 0;

;				if(ox1 == 1){				
					Euler[NoEr]    = theta[4] + theta[7];
					Euler[NoEr+1]  = theta[5] - theta[8];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4] + phi[6];
					Euler[NoFr1+1]  = phi[5] - phi[7];

					Euler[NoFr2]    = psi[4] + psi[6];
					Euler[NoFr2+1]  = psi[5] + psi[7];
				}
				else if(ox1 == 2){
					Euler[NoEr]    = theta[4] + theta[7];
					Euler[NoEr+1]  = theta[5] + theta[8];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4] + phi[6];
					Euler[NoFr1+1]  = phi[5] + phi[7];

					Euler[NoFr2]    = psi[4] + psi[6];
					Euler[NoFr2+1]  = psi[5] + psi[7];
				}
				else if(ox1 == 3){
					Euler[NoEr]    = theta[4];
					Euler[NoEr+1]  = theta[5];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4];
					Euler[NoFr1+1]  = phi[5];

					Euler[NoFr2]    = psi[4];
					Euler[NoFr2+1]  = psi[5];
				}
				else 
					{goto on_error;}
				}
				else{
				/* For Periodic boundary condition ox1 == 4 */
					Euler[NoEr+3]  = theta[7];
					IEuler[NoEr+3]   = 3*(j-js)*Nx;
					JEuler[NoEr+3]   = 3*(j-js)*Nx + 3*(i-is);
					Euler[NoEr+4]  = theta[8];
					IEuler[NoEr+4]   = 3*(j-js)*Nx + 1;
					JEuler[NoEr+4]   = 3*(j-js)*Nx + 3*(i-is);

					Euler[NoFr1+2]  = phi[6];
					IEuler[NoFr1+2]   =  3*(j-js)*Nx;
					JEuler[NoFr1+2]   = 3*(j-js)*Nx + 3*(i-is) + 1;
					Euler[NoFr1+3]  = phi[7];
					IEuler[NoFr1+3]   = 3*(j-js)*Nx + 1;
					JEuler[NoFr1+3]   = 3*(j-js)*Nx + 3*(i-is) + 1;

					Euler[NoFr2+2]  = psi[6];
					IEuler[NoFr2+2]   = 3*(j-js)*Nx;
					JEuler[NoFr2+2]   = 3*(j-js)*Nx + 3*(i-is) + 2;
					Euler[NoFr2+3]  = psi[7];
					IEuler[NoFr2+3]   = 3*(j-js)*Nx + 2;
					JEuler[NoFr2+3]   = 3*(j-js)*Nx + 3*(i-is) + 2;

					Euler[NoEr]    = theta[4];
					Euler[NoEr+1]  = theta[5];
					Euler[NoEr+2]  = theta[6];

					Euler[NoFr1]    = phi[4];
					Euler[NoFr1+1]  = phi[5];

					Euler[NoFr2]    = psi[4];
					Euler[NoFr2+1]  = psi[5];
				}
			}/* End i == ie */
			else{

				for(m=0; m<7; m++){
					Euler[NoEr -2+m] = theta[2+m];
					JEuler[NoEr-2+m] = 3*(j-js)*Nx + 3*(i-is);
				}
				for(m=0; m<6; m++){
					Euler[NoFr1-2 +m] = phi[2+m];
					JEuler[NoFr1-2+m] = 3*(j-js)*Nx + 3*(i-is) + 1;
					Euler[NoFr2-2 +m] = psi[2+m];
					JEuler[NoFr2-2+m] = 3*(j-js)*Nx + 3*(i-is) + 2;			
				}
				
				IEuler[NoEr-2] = 3*(j-js)*Nx + 3*(i-is) - 3;				
				IEuler[NoEr-1] = 3*(j-js)*Nx + 3*(i-is) - 2;
				IEuler[NoEr] = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoEr+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				IEuler[NoEr+2] = 3*(j-js)*Nx + 3*(i-is) + 2;				
				IEuler[NoEr+3] = 3*(j-js)*Nx + 3*(i-is) + 3;
				IEuler[NoEr+4] = 3*(j-js)*Nx + 3*(i-is) + 4;	

				IEuler[NoFr1-2] = 3*(j-js)*Nx + 3*(i-is) - 3;				
				IEuler[NoFr1-1] = 3*(j-js)*Nx + 3*(i-is) - 2;
				IEuler[NoFr1] = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoFr1+1] = 3*(j-js)*Nx + 3*(i-is) + 1;
				IEuler[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is) + 3;				
				IEuler[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is) + 4;

				IEuler[NoFr2-2] = 3*(j-js)*Nx + 3*(i-is) - 3;				
				IEuler[NoFr2-1] = 3*(j-js)*Nx + 3*(i-is) - 1;
				IEuler[NoFr2] = 3*(j-js)*Nx + 3*(i-is);				
				IEuler[NoFr2+1] = 3*(j-js)*Nx + 3*(i-is) + 2;
				IEuler[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is) + 3;				
				IEuler[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is) + 5;				
			}/* End i!= is && i!= ie */

/*-----------------------------------------------------------------------------------------------*/

			if(j == js){
				Euler[NoEr+5]  = theta[9];
				Euler[NoEr+6]  = theta[10];
				Euler[NoFr1+4] = phi[8];
				Euler[NoFr1+5] = phi[9];
				Euler[NoFr2+4] = psi[8];
				Euler[NoFr2+5] = psi[9];

				IEuler[NoEr+5] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoEr+6] = 3*(j-js+1)*Nx + 3*(i-is) + 2;
				JEuler[NoEr+5] = 3*(j-js)*Nx + 3*(i-is);
  				JEuler[NoEr+6] = 3*(j-js)*Nx + 3*(i-is);

				IEuler[NoFr1+4] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoFr1+5] = 3*(j-js+1)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  				JEuler[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
				IEuler[NoFr2+4] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoFr2+5] = 3*(j-js+1)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  				JEuler[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is) + 2;

				/* judge boundary condition for ix2 */
				if(ix2 != 4){
					/* non-periodic boundary condition */
					Euler[NoEr-4]  = 0.0;
					Euler[NoEr-3]  = 0.0;
					Euler[NoFr1-4] = 0.0;
					Euler[NoFr1-3] = 0.0;
					Euler[NoFr2-4] = 0.0;
					Euler[NoFr2-3] = 0.0;

					IEuler[NoEr-4] = 0;
					IEuler[NoEr-3] = 0;
					JEuler[NoEr-4] = 0;
  					JEuler[NoEr-3] = 0;

					IEuler[NoFr1-4] = 0;
					IEuler[NoFr1-3] = 0;
					JEuler[NoFr1-4] = 0;
  					JEuler[NoFr1-3] = 0;
				
					IEuler[NoFr2-4] = 0;
					IEuler[NoFr2-3] = 0;
					JEuler[NoFr2-4] = 0;
  					JEuler[NoFr2-3] = 0;

					if(ix2 == 1){
						Euler[NoEr]    += theta[0];
						Euler[NoEr+2]  -= theta[1];

						Euler[NoFr1]   += phi[0];
						Euler[NoFr1+1] += phi[1];

						Euler[NoFr2]    += psi[0];
						Euler[NoFr2+1]  -= psi[1];
					}
					else if(ix2 == 2){
						Euler[NoEr]    += theta[0];
						Euler[NoEr+2]  += theta[1];

						Euler[NoFr1]   += phi[0];
						Euler[NoFr1+1] += phi[1];

						Euler[NoFr2]    += psi[0];
						Euler[NoFr2+1]  += psi[1];

					}
					else if(ix2 == 3){
						/* Do nothing*/
					}
					else
						{goto on_error;}
				}
				else{
					/* for periodic boundary condition */
					Euler[NoEr-4]  = theta[0];
					Euler[NoEr-3]  = theta[1];
					Euler[NoFr1-4] = phi[0];
					Euler[NoFr1-3] = phi[1];
					Euler[NoFr2-4] = psi[0];
					Euler[NoFr2-3] = psi[1];

					/* Here, we assume je-js>1 , otherwise it is wrong */

					IEuler[NoEr-4] = 3*(je-js)*Nx + 3*(i-is);
					IEuler[NoEr-3] = 3*(je-js)*Nx + 3*(i-is) + 2; 
					JEuler[NoEr-4] = 3*(j-js)*Nx + 3*(i-is);
  					JEuler[NoEr-3] = 3*(j-js)*Nx + 3*(i-is);

					IEuler[NoFr1-4] = 3*(je-js)*Nx + 3*(i-is);
					IEuler[NoFr1-3] = 3*(je-js)*Nx + 3*(i-is) + 1;
					JEuler[NoFr1-4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  					JEuler[NoFr1-3] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
					IEuler[NoFr2-4] = 3*(je-js)*Nx + 3*(i-is);
					IEuler[NoFr2-3] = 3*(je-js)*Nx + 3*(i-is) + 2;
					JEuler[NoFr2-4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  					JEuler[NoFr2-3] = 3*(j-js)*Nx + 3*(i-is) + 2;
				}/* End periodic boundary condition */
			}/* End j==js */
			else if(j == je){
				Euler[NoEr-4]  = theta[0];
				Euler[NoEr-3]  = theta[1];
				Euler[NoFr1-4] = phi[0];
				Euler[NoFr1-3] = phi[1];
				Euler[NoFr2-4] = psi[0];
				Euler[NoFr2-3] = psi[1];

				IEuler[NoEr-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoEr-3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;
				JEuler[NoEr-4] = 3*(j-js)*Nx + 3*(i-is);
  				JEuler[NoEr-3] = 3*(j-js)*Nx + 3*(i-is);

				IEuler[NoFr1-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoFr1-3] = 3*(j-js-1)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1-4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  				JEuler[NoFr1-3] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
				IEuler[NoFr2-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoFr2-3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2-4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  				JEuler[NoFr2-3] = 3*(j-js)*Nx + 3*(i-is) + 2;

				/* judge boundary condition for ox2 */
				if(ox2 != 4){
					/* non-periodic boundary condition */
					Euler[NoEr+5]  = 0.0;
					Euler[NoEr+6]  = 0.0;
					Euler[NoFr1+4] = 0.0;
					Euler[NoFr1+5] = 0.0;
					Euler[NoFr2+4] = 0.0;
					Euler[NoFr2+5] = 0.0;

					IEuler[NoEr+5] = 0;
					IEuler[NoEr+6] = 0;
					JEuler[NoEr+5] = 0;
  					JEuler[NoEr+6] = 0;

					IEuler[NoFr1+4] = 0;
					IEuler[NoFr1+5] = 0;
					JEuler[NoFr1+4] = 0;
  					JEuler[NoFr1+5] = 0;
				
					IEuler[NoFr2+4] = 0;
					IEuler[NoFr2+5] = 0;
					JEuler[NoFr2+4] = 0;
  					JEuler[NoFr2+5] = 0;

					if(ox2 == 1){
						Euler[NoEr]    += theta[9];
						Euler[NoEr+2]  -= theta[10];

						Euler[NoFr1]   += phi[8];
						Euler[NoFr1+1] += phi[9];

						Euler[NoFr2]    += psi[8];
						Euler[NoFr2+1]  -= psi[9];
					}
					else if(ox2 == 2){
						Euler[NoEr]    += theta[9];
						Euler[NoEr+2]  += theta[10];

						Euler[NoFr1]   += phi[8];
						Euler[NoFr1+1] += phi[9];

						Euler[NoFr2]    += psi[8];
						Euler[NoFr2+1]  += psi[9];

					}
					else if(ox2 == 3){
						/* Do nothing*/
					}
					else
						{goto on_error;}
				}
				else{
					/* for periodic boundary condition */
					Euler[NoEr+5]  = theta[9];
					Euler[NoEr+6]  = theta[10];
					Euler[NoFr1+4] = phi[8];
					Euler[NoFr1+5] = phi[9];
					Euler[NoFr2+4] = psi[8];
					Euler[NoFr2+5] = psi[9];

					/* Here, we assume je-js>1 , otherwise it is wrong */

					IEuler[NoEr+5] = 3*(i-is);
					IEuler[NoEr+6] = 3*(i-is) + 2; 
					JEuler[NoEr+5] = 3*(j-js)*Nx + 3*(i-is);
  					JEuler[NoEr+6] = 3*(j-js)*Nx + 3*(i-is);

					IEuler[NoFr1+4] = 3*(i-is);
					IEuler[NoFr1+5] = 3*(i-is) + 1;
					JEuler[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  					JEuler[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
					IEuler[NoFr2+4] = 3*(i-is);
					IEuler[NoFr2+5] = 3*(i-is) + 2;
					JEuler[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  					JEuler[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is) + 2;
				}/* End periodic boundary condition */ 
			}/* End j==je */
			else{
				Euler[NoEr-4]  = theta[0];
				Euler[NoEr-3]  = theta[1];
				Euler[NoFr1-4] = phi[0];
				Euler[NoFr1-3] = phi[1];
				Euler[NoFr2-4] = psi[0];
				Euler[NoFr2-3] = psi[1];
				Euler[NoEr+5]  = theta[9];
				Euler[NoEr+6]  = theta[10];
				Euler[NoFr1+4] = phi[8];
				Euler[NoFr1+5] = phi[9];
				Euler[NoFr2+4] = psi[8];
				Euler[NoFr2+5] = psi[9];

				IEuler[NoEr-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoEr-3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;
				JEuler[NoEr-4] = 3*(j-js)*Nx + 3*(i-is);
  				JEuler[NoEr-3] = 3*(j-js)*Nx + 3*(i-is);

				IEuler[NoFr1-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoFr1-3] = 3*(j-js-1)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1-4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  				JEuler[NoFr1-3] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
				IEuler[NoFr2-4] = 3*(j-js-1)*Nx + 3*(i-is);
				IEuler[NoFr2-3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2-4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  				JEuler[NoFr2-3] = 3*(j-js)*Nx + 3*(i-is) + 2;


				IEuler[NoEr+5] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoEr+6] = 3*(j-js+1)*Nx + 3*(i-is) + 2;
				JEuler[NoEr+5] = 3*(j-js)*Nx + 3*(i-is);
  				JEuler[NoEr+6] = 3*(j-js)*Nx + 3*(i-is);

				IEuler[NoFr1+4] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoFr1+5] = 3*(j-js+1)*Nx + 3*(i-is) + 1;
				JEuler[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is) + 1;
  				JEuler[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is) + 1;
				
				IEuler[NoFr2+4] = 3*(j-js+1)*Nx + 3*(i-is);
				IEuler[NoFr2+5] = 3*(j-js+1)*Nx + 3*(i-is) + 2;
				JEuler[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is) + 2;
  				JEuler[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is) + 2;
			}/* End j!= js && j!= je*/
		}/* End loop i */
	}/* End loop j */

		/* Solve the matrix equation with retarded GMRES method */
		/* solution is input in INIguess */
		/* ITR_MAX, the maximum number of (outer) iterations to take.
    		 * MR, the maximum number of (inner) iterations to take.    MR must be less than N.*/

		int ITR_max = 50;
		int MR;
		double tol_abs = 1.0E-12;
    		double tol_rel = 1.0E-12;
		if(Nmatrix < 20) MR = Nmatrix - 1;
		else MR = 20;

		/* Note that in the C subroutine, J and I are actually I and J used in Athena */

		mgmres_st ( 3*Nmatrix, NZ_NUM, JEuler, IEuler, Euler, INIguess, RHSEuler, ITR_max, MR, tol_abs,  tol_rel );	

		
	/* update the radiation quantities in the mesh */	
	for(j=js;j<=je;j++){
		for(i=is; i<=ie; i++){
		pG->U[ks][j][i].Er	= INIguess[3*(j-js)*Nx + 3*(i-is)];
		if(pG->U[ks][j][i].Er < 0.0)
			fprintf(stderr,"[BackEuler_2d]: Negative Radiation energy: %e\n",pG->U[ks][j][i].Er);

		pG->U[ks][j][i].Fr1	= INIguess[3*(j-js)*Nx + 3*(i-is)+1];
		pG->U[ks][j][i].Fr2	= INIguess[3*(j-js)*Nx + 3*(i-is)+2];		
		}
	}
	/* May need to update Edd_11 */


/* Update the ghost zones for different boundary condition to be used later */
		bvals_radMHD(pM);


/*-----------Finish---------------------*/

  
	/* Free the temporary variables just used for this grid calculation*/
	rad_hydro_destruct_2d(Nmatrix);
	
  return;	


	on_error:
	
	rad_hydro_destruct_2d(Nmatrix);
	ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}



/*-------------------------------------------------------------------------*/
/* rad_hydro_init_1d: function to allocate memory used just for radiation variables */
/* rad_hydro_destruct_1d(): function to free memory */
void rad_hydro_init_2d(int Ngrids)
{

/* Ngrids = Nx * Ny

/* The matrix Euler is stored as a compact form. See $2.4 of numerical recipes */


	if ((RHSEuler = (Real*)malloc(3*Ngrids*sizeof(Real))) == NULL) goto on_error;
	if ((INIguess = (Real*)malloc(3*Ngrids*sizeof(Real))) == NULL) goto on_error;
	/* RHS_guess start from 1 */

	if ((Euler = (Real*)malloc(Ngrids*31*sizeof(Real))) == NULL) goto on_error;
	if ((IEuler = (int*)malloc(Ngrids*31*sizeof(int ))) == NULL) goto on_error;
	if ((JEuler = (int*)malloc(Ngrids*31*sizeof(int ))) == NULL) goto on_error;

	
	return;

	on_error:
    	rad_hydro_destruct_2d(Ngrids);
	ath_error("[BackEuler]: malloc returned a NULL pointer\n");
}


void rad_hydro_destruct_2d(int Ngrids)
{

	if (RHSEuler	!= NULL) free(RHSEuler);
	if (Euler	!= NULL) free(Euler);
	if (IEuler	!= NULL) free(IEuler);
	if (JEuler	!= NULL) free(JEuler);
	if (INIguess	!= NULL) free(INIguess);
	
}


#endif /* radMHD_INTEGRATOR */
