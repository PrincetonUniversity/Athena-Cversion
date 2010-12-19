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



/*================================*/
/* For the matrix solver */
/* we use lis library now */
#include <lis.h>

/*===============================*/

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */







/*Store the index of all non-zero elements */
/* The matrix coefficient, this is oneD in GMRES solver */
static LIS_MATRIX Euler;

static LIS_SOLVER solver;
/* Right hand side of the Matrix equation */
static LIS_VECTOR RHSEuler;
/* RHSEuler is 0-------3*Nx*Ny-1; only for GMRES method */ 
static LIS_VECTOR INIguess;
/* Used for initial guess solution, also return the correct solution */

static LIS_SCALAR *Value;
static int *indexValue;
static int *ptr;






/********Public function****************/
/*-------BackEuler_2d(): Use back euler method to update E_r and Fluxr-----------*/
/* we may need to use variables to represent the indices to simplify the code */



void BackEuler_2d(MeshS *pM)
{
/* Right now, only work for one domain. Modified later for SMR */


  	GridS *pG=(pM->Domain[0][0].Grid);
	Real hdtodx1 = 0.5*pG->dt/pG->dx1;
	Real hdtodx2 = 0.5 * pG->dt/pG->dx2;
	Real dt = pG->dt;
	int is = pG->is, ie = pG->ie;
  	int i, j, m, NoEr, NoFr1, NoFr2;
	int js = pG->js, je = pG->je;
	int ks = pG->ks;
	int Nx, Ny, Nmatrix, NZ_NUM, lines, count;
	
	/* NZ_NUM is the number of non-zero elements in Matrix. It may change if periodic boundary condition is applied */
	/* lines is the size of the matrix */
	/* count is the number of total non-zeros before that row */

	
	Real temperature, velocity_x, velocity_y, pressure;

	Real Ci0, Ci1, Cj0, Cj1;
	/* This is equivilent to Cspeeds[] in 1D */

	
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

	lines  = 3 * Nmatrix;

	NZ_NUM = 31 * Nmatrix; 
	
	NoEr = 0;/* Position of first non-zero element in row Er */
	NoFr1 = 0;/* Position of first non-zero element in row Fr1 */
	NoFr2 = 0;/* POsition of first non-zero element in row Fr2 */
	count = 0;
	/* For non-periodic boundary condition, this number will change */
	
/* setting for LIS library. Now this is noly for serial case. Need to change for parallal case */
/*	lis_initialize(0,0);
*/
	

	/* For temporary use only */
	int index,Matrixiter;
	Real tempvalue;
	
	
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
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c) * (Gamma - 1.0);
#endif

    		temperature = pressure / (pG->U[ks][j][i].d * R_ideal);
		/* RHSEuler[0...N-1]  */
		Sigma_a = pG->U[ks][j][i].Sigma_a;

		/*-----------------------------*/		
    		tempvalue   = pG->U[ks][j][i].Er + Crat * dt * Sigma_a 
				* temperature * temperature * temperature * temperature;
		index = 3*(j-js)*Nx + 3*(i-is);
		lis_vector_set_value(LIS_INS_VALUE,index,tempvalue,RHSEuler);

		/*----------------------------*/
    		tempvalue = pG->U[ks][j][i].Fr1 + dt * Sigma_a
				* temperature * temperature * temperature * temperature * pG->U[ks][j][i].M1 / pG->U[ks][j][i].d;
		++index;
		lis_vector_set_value(LIS_INS_VALUE,index,tempvalue,RHSEuler);
		
		/*-------------------------*/
		tempvalue = pG->U[ks][j][i].Fr2 + dt * Sigma_a
				* temperature * temperature * temperature * temperature * pG->U[ks][j][i].M2 / pG->U[ks][j][i].d;
		++index;
		lis_vector_set_value(LIS_INS_VALUE,index,tempvalue,RHSEuler);

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

			/* Subtract some value */
			tempvalue = -(theta[2] * pG->U[ks][j][i-1].Er + theta[3] * pG->U[ks][j][i-1].Fr1);
			index = 3*(j-js)*Nx + 3*(i-is);
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);


			tempvalue = -(phi[2] * pG->U[ks][j][i-1].Er + phi[3] * pG->U[ks][j][i-1].Fr1);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);

			tempvalue = -(psi[2] * pG->U[ks][j][i-1].Er + psi[3] * pG->U[ks][j][i-1].Fr2);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);
			
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

			tempvalue = -(theta[7] * pG->U[ks][j][i+1].Er + theta[8] * pG->U[ks][j][i+1].Fr1);
			index = 3*(j-js)*Nx + 3*(i-is);
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);

			tempvalue = -(phi[6] * pG->U[ks][j][i+1].Er + phi[7] * pG->U[ks][j][i+1].Fr1);
			++index; 
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);

			tempvalue = -(psi[6] * pG->U[ks][j][i+1].Er + psi[7] * pG->U[ks][j][i+1].Fr2);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);
					
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

			tempvalue = -(theta[0] * pG->U[ks][j-1][i].Er + theta[1] * pG->U[ks][j-1][i].Fr2);
			index = 3*(j-js)*Nx + 3 * (i - is);
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);


			tempvalue = -(phi[0] * pG->U[ks][j-1][i].Er + phi[1] * pG->U[ks][j-1][i].Fr1);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);


			tempvalue = -(psi[0] * pG->U[ks][j-1][i].Er + psi[1] * pG->U[ks][j-1][i].Fr2);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);
				
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

			tempvalue = -(theta[9] * pG->U[ks][j+1][i].Er + theta[10] * pG->U[ks][j+1][i].Fr2);
			index = 3*(j-js)*Nx + 3*(i-is);
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);

			tempvalue = -(phi[8] * pG->U[ks][j+1][i].Er + phi[9] * pG->U[ks][j+1][i].Fr1);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);


			tempvalue = -(psi[8] * pG->U[ks][j+1][i].Er + psi[9] * pG->U[ks][j+1][i].Fr2);
			++index;
			lis_vector_set_value(LIS_ADD_VALUE,index,tempvalue,RHSEuler);
			
		}
				
		}
	}

	
/*--------------------Note--------------------*/


/* Step 1b: Setup the Matrix */
		
 	/* First, setup the guess solution. Guess solution is the solution from last time step */
	for(j=js; j<=je; j++){
		for(i=is; i<=ie; i++){
			lis_vector_set_value(LIS_INS_VALUE,3*(j-js)*Nx + 3*(i-is),pG->U[ks][j][i].Er,INIguess);
			lis_vector_set_value(LIS_INS_VALUE,3*(j-js)*Nx + 3*(i-is)+1,pG->U[ks][j][i].Fr1,INIguess);
			lis_vector_set_value(LIS_INS_VALUE,3*(j-js)*Nx + 3*(i-is)+2,pG->U[ks][j][i].Fr2,INIguess);
			
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


		if(i == is){
			if(j == js){
				if(ix1 != 4){						
					if(ix2 !=4){
						NZ_NUM -= 12;
						NoEr = count;
						NoFr1 = count + 7;
						NoFr2 = count + 13;
						count += 19;
					}
					else{
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
				}/* Non periodic for x1 */
				else{
					if(ix2 !=4){
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
					else{
						NoEr = count;
						NoFr1 = count + 11;
						NoFr2 = count + 21;
						count += 31;
					}
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;
				/* For Er */
				for(m=0; m<5; m++)
					Value[NoEr+m] = theta[4+m];
				
				for(m=0;m<5;m++)
					indexValue[NoEr+m] = m;
				
				if(ix1 != 4){
						
					Value[NoEr+5]	= theta[9];
					Value[NoEr+6]	= theta[10];
					indexValue[NoEr+5] = 3 * Nx;
					indexValue[NoEr+6] = 3 * Nx + 2;
					
					if(ix2 == 4){
						Value[NoEr+7] = theta[0];
						Value[NoEr+8] = theta[1];
						indexValue[NoEr+7] = 3*(je-js)*Nx;
						indexValue[NoEr+8] = 3*(je-js)*Nx + 2;
					}
				}
				else {
					Value[NoEr+5]	= theta[2];
					Value[NoEr+6]	= theta[3];
					Value[NoEr+7]	= theta[9];
					Value[NoEr+8]	= theta[10];
					indexValue[NoEr+5] = 3 * (ie - is);
					indexValue[NoEr+6] = 3 * (ie - is) + 1;
					indexValue[NoEr+7] = 3 * Nx;
					indexValue[NoEr+8] = 3 * Nx + 2;
					
					if(ix2 == 4){
						Value[NoEr+9] = theta[0];
						Value[NoEr+10] = theta[1];
						indexValue[NoEr+9] = 3*(je-js)*Nx;
						indexValue[NoEr+10] = 3*(je-js)*Nx + 2;
					}
				}

				
						
				/* For Fr1 */
				for(m=0; m<4; m++)
					Value[NoFr1+m] = phi[4+m];

				indexValue[NoFr1+0] = 0;
				indexValue[NoFr1+1] = 1;
				indexValue[NoFr1+2] = 3;
				indexValue[NoFr1+3] = 4;
				
				if(ix1 != 4){
					Value[NoFr1+4]	= phi[8];
					Value[NoFr1+5]	= phi[9];
					indexValue[NoFr1+4] = 3 * Nx;
					indexValue[NoFr1+5] = 3 * Nx + 1;
					if(ix2 == 4){
						Value[NoFr1+6] = phi[0];
						Value[NoFr1+7] = phi[1];
						indexValue[NoFr1+6] = 3*(je-js)*Nx;
						indexValue[NoFr1+7] = 3*(je-js)*Nx + 1;
					}
				}
				else{
					Value[NoFr1+4]	= phi[2];
					Value[NoFr1+5]	= phi[3];
					Value[NoFr1+6]	= phi[8];
					Value[NoFr1+7]	= phi[9];
					indexValue[NoFr1+4] = 3 * (ie - is);
					indexValue[NoFr1+5] = 3 * (ie - is) + 1;
					indexValue[NoFr1+6] = 3 * Nx;
					indexValue[NoFr1+7] = 3 * Nx + 1;
					if(ix2 == 4){
						Value[NoFr1+8] = phi[0];
						Value[NoFr1+9] = phi[1];
						indexValue[NoFr1+8] = 3*(je-js)*Nx;
						indexValue[NoFr1+9] = 3*(je-js)*Nx + 1;
					}
				}				

					

				/* For Fr2 */
					
				for(m=0; m<4; m++)
					Value[NoFr2+m] = psi[4+m];

				indexValue[NoFr2+0] = 0;
				indexValue[NoFr2+1] = 2;
				indexValue[NoFr2+2] = 3;
				indexValue[NoFr2+3] = 5;
				
				if(ix1 != 4){
					Value[NoFr2+4]	= psi[8];
					Value[NoFr2+5]	= psi[9];
					indexValue[NoFr2+4] = 3 * Nx;
					indexValue[NoFr2+5] = 3 * Nx + 2;
					if(ix2 == 4){
						Value[NoFr2+6] = psi[0];
						Value[NoFr2+7] = psi[1];
						indexValue[NoFr2+6] = 3*(je-js)*Nx;
						indexValue[NoFr2+7] = 3*(je-js)*Nx + 2;
					}
				}
				else{
					Value[NoFr2+4]	= psi[2];
					Value[NoFr2+5]	= psi[3];
					Value[NoFr2+6]	= psi[8];
					Value[NoFr2+7]	= psi[9];
					indexValue[NoFr2+4] = 3 * (ie - is);
					indexValue[NoFr2+5] = 3 * (ie - is) + 2;
					indexValue[NoFr2+6] = 3 * Nx;
					indexValue[NoFr2+7] = 3 * Nx + 2;
					if(ix2 == 4){
						Value[NoFr2+8] = psi[0];
						Value[NoFr2+9] = psi[1];
						indexValue[NoFr2+8] = 3*(je-js)*Nx;
						indexValue[NoFr2+9] = 3*(je-js)*Nx + 2;
					}
				}


				/* other ix1 boundary condition */
				if(ix1 == 1 || ix1 == 5){
					
					Value[NoEr+0] += theta[2];
					Value[NoEr+1] -= theta[3];
			
					Value[NoFr1+0]+= phi[2];
					Value[NoFr1+1]-= phi[3];
				
					Value[NoFr2+0]+= psi[2];
					Value[NoFr2+1]+= psi[3];
				}
				else if(ix1 == 2){
					Value[NoEr+0] += theta[2];
					Value[NoEr+1] += theta[3];
			
					Value[NoFr1+0]+= phi[2];
					Value[NoFr1+1]+= phi[3];
				
					Value[NoFr2+0]+= psi[2];
					Value[NoFr2+1]+= psi[3];
				}
				else if(ix1 == 3){

					/* Do nothing */
				}
				else {
					goto on_error;
				}
				
				/* other ix2 boundary condition */	

				if(ix2 == 1 || ix2 == 5){
					
					Value[NoEr+0] += theta[0];
					Value[NoEr+2] -= theta[1];
			
					Value[NoFr1+0]+= phi[0];
					Value[NoFr1+1]+= phi[1];
				
					Value[NoFr2+0]+= psi[0];
					Value[NoFr2+1]-= psi[1];
				}
				else if(ix2 == 2){
					Value[NoEr+0] += theta[0];
					Value[NoEr+2] += theta[1];
			
					Value[NoFr1+0]+= phi[0];
					Value[NoFr1+1]+= phi[1];
				
					Value[NoFr2+0]+= psi[0];
					Value[NoFr2+1]+= psi[1];
				}
				else if(ix2 == 3){

					/* Do nothing */
				}
				else {
					goto on_error;
				}
				
			}/* End j == js */
			else if(j == je){
				if(ix1 != 4){						
					if(ox2 !=4){
						NZ_NUM -= 12;
						NoEr = count;
						NoFr1 = count + 7;
						NoFr2 = count + 13;
						count += 19;
					}
					else{
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
				}/* Non periodic for x1 */
				else{
					if(ox2 !=4){
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
					else{
						NoEr = count;
						NoFr1 = count + 11;
						NoFr2 = count + 21;
						count += 31;
					}
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;
				
				
				/* Now the important thing is ox2, which determines the first non-zero element */
				
				
				if(ox2 != 4){
					/* The following is true no matter ix1 == 4 or not */

					/* For Er */
					Value[NoEr] = theta[0];
					Value[NoEr+1] = theta[1];
							
					indexValue[NoEr] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is) + 2;

					/* For Fr1 */

					Value[NoFr1] = phi[0];
					Value[NoFr1+1] = phi[1];
							
					indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is) + 1;

					/* For Fr2 */
					
					Value[NoFr2] = psi[0];
					Value[NoFr2+1] = psi[1];
							
					indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is) + 2;

					
					/* for Er */
					for(m=0; m<5; m++){
						Value[NoEr+2+m] = theta[4+m];
						indexValue[NoEr+2+m] = 3*(j-js)*Nx + 3*(i-is) + m;
					}

					/* For Fr1 */
					Value[NoFr1+2] = phi[4];
					Value[NoFr1+3] = phi[5];
					Value[NoFr1+4] = phi[6];
					Value[NoFr1+5] = phi[7];

					indexValue[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is) + 1;
					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is) + 3;
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is) + 4;

					/* For Fr2 */
					Value[NoFr2+2] = psi[4];
					Value[NoFr2+3] = psi[5];
					Value[NoFr2+4] = psi[6];
					Value[NoFr2+5] = psi[7];

					indexValue[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is) + 2;
					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is) + 3;
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is) + 5;				

					
					if (ix1 == 4) {

						
						/* For Er */
						Value[NoEr+7] = theta[2];
						Value[NoEr+8] = theta[3];
							
						indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(ie-is) + 1;

						/* For Fr1 */

						Value[NoFr1+6] = phi[2];
						Value[NoFr1+7] = phi[3];
							
						indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(ie-is) + 1;

						/* For Fr2 */
					
						Value[NoFr2+6] = psi[2];
						Value[NoFr2+7] = psi[3];
							
						indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(ie-is) + 2;

					} /* for periodic boundary condition */
					else if(ix1 == 1 || ix1 == 5){
					
						Value[NoEr+2] += theta[2];
						Value[NoEr+3] -= theta[3];
			
						Value[NoFr1+2]+= phi[2];
						Value[NoFr1+3]-= phi[3];
				
						Value[NoFr2+2]+= psi[2];
						Value[NoFr2+3]+= psi[3];
					}
					else if(ix1 == 2){
						Value[NoEr+2] += theta[2];
						Value[NoEr+3] += theta[3];
			
						Value[NoFr1+2]+= phi[2];
						Value[NoFr1+3]+= phi[3];
				
						Value[NoFr2+2]+= psi[2];
						Value[NoFr2+3]+= psi[3];
					}
					else if(ix1 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}


					/* other ox2 boundary condition */
					if(ox2 == 1 || ox2 == 5){
					
						Value[NoEr+2] += theta[9];
						Value[NoEr+4] -= theta[10];
			
						Value[NoFr1+2]+= phi[8];
						Value[NoFr1+3]+= phi[9];
				
						Value[NoFr2+2]+= psi[8];
						Value[NoFr2+3]-= psi[9];
					}
					else if(ox2 == 2){
						Value[NoEr+2] += theta[9];
						Value[NoEr+4] += theta[10];
			
						Value[NoFr1+2]+= phi[8];
						Value[NoFr1+3]+= phi[9];
				
						Value[NoFr2+2]+= psi[8];
						Value[NoFr2+3]+= psi[9];
					}
					else if(ox2 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}
				}/* Non-periodic for x2 */
				else{

					/* The following is true no matter ix1 == 4 or not */
					/* For Er */
					Value[NoEr] = theta[9];
					Value[NoEr+1] = theta[10];
					Value[NoEr+2] = theta[0];
					Value[NoEr+3] = theta[1];
					
					indexValue[NoEr] = 3*(i-is);
					indexValue[NoEr+1] =3*(i-is) + 2;		
					indexValue[NoEr+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;

					/* For Fr1 */
					
					Value[NoFr1] = phi[8];
					Value[NoFr1+1] = phi[9];
					Value[NoFr1+2] = phi[0];
					Value[NoFr1+3] = phi[1];
					
					indexValue[NoFr1] = 3*(i-is);
					indexValue[NoFr1+1] = 3*(i-is) + 1;		
					indexValue[NoFr1+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+3] = 3*(j-js-1)*Nx + 3*(i-is) + 1;

					/* For Fr2 */
					
					Value[NoFr2] = psi[8];
					Value[NoFr2+1] = psi[9];
					Value[NoFr2+2] = psi[0];
					Value[NoFr2+3] = psi[1];
					
					indexValue[NoFr2] = 3*(i-is);
					indexValue[NoFr2+1] = 3*(i-is) + 2;		
					indexValue[NoFr2+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;

					
					/* for Er */
					for(m=0; m<5; m++){
						Value[NoEr+4+m] = theta[4+m];
						indexValue[NoEr+4+m] = 3*(j-js)*Nx + 3*(i-is) + m;
					}

					/* For Fr1 */
					Value[NoFr1+4] = phi[4];
					Value[NoFr1+5] = phi[5];
					Value[NoFr1+6] = phi[6];
					Value[NoFr1+7] = phi[7];

					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is) + 1;
					indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(i-is) + 3;
					indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(i-is) + 4;

					/* For Fr2 */
					Value[NoFr2+4] = psi[4];
					Value[NoFr2+5] = psi[5];
					Value[NoFr2+6] = psi[6];
					Value[NoFr2+7] = psi[7];

					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is) + 2;
					indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(i-is) + 3;
					indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(i-is) + 5;				

					
					if (ix1 == 4) {

						
						/* For Er */
						Value[NoEr+9] = theta[2];
						Value[NoEr+10] = theta[3];
							
						indexValue[NoEr+9] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoEr+10] = 3*(j-js)*Nx + 3*(ie-is) + 1;

						/* For Fr1 */

						Value[NoFr1+8] = phi[2];
						Value[NoFr1+9] = phi[3];
							
						indexValue[NoFr1+8] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoFr1+9] = 3*(j-js)*Nx + 3*(ie-is) + 1;

						/* For Fr2 */
					
						Value[NoFr2+8] = psi[2];
						Value[NoFr2+9] = psi[3];
							
						indexValue[NoFr2+8] = 3*(j-js)*Nx + 3*(ie-is);
						indexValue[NoFr2+9] = 3*(j-js)*Nx + 3*(ie-is) + 2;

					} /* for periodic boundary condition */
					else if(ix1 == 1 || ix1 == 5){
					
						Value[NoEr+4] += theta[2];
						Value[NoEr+5] -= theta[3];
			
						Value[NoFr1+4]+= phi[2];
						Value[NoFr1+5]-= phi[3];
				
						Value[NoFr2+4]+= psi[2];
						Value[NoFr2+5]+= psi[3];
					}
					else if(ix1 == 2){
						Value[NoEr+4] += theta[2];
						Value[NoEr+5] += theta[3];
			
						Value[NoFr1+4]+= phi[2];
						Value[NoFr1+5]+= phi[3];
				
						Value[NoFr2+4]+= psi[2];
						Value[NoFr2+5]+= psi[3];
					}
					else if(ix1 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}

				}/* periodic for x2 */
			} /* End j == je */
			else {
				if(ix1 != 4){						
					NZ_NUM -= 6;
					NoEr = count;
					NoFr1 = count + 9;
					NoFr2 = count + 17;
					count += 25;
					
				}/* Non periodic for x1 */
				else{
				
					NoEr = count;
					NoFr1 = count + 11;
					NoFr2 = count + 21;
					count += 31;
				
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;

				/* The following is true no matter ix1== 4 or not */

				/* For Er */

				Value[NoEr] = theta[0];
				Value[NoEr+1] = theta[1];

				Value[NoEr+2] = theta[4];
				Value[NoEr+3] = theta[5];
				Value[NoEr+4] = theta[6];
				Value[NoEr+5] = theta[7];
				Value[NoEr+6] = theta[8];

				indexValue[NoEr] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

				indexValue[NoEr+2] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoEr+3] = 3*(j-js)*Nx + 3*(i-is)+1;
				indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is)+2;
				indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is)+3;
				indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is)+4;

				/* For Fr1 */			
				
				Value[NoFr1] = phi[0];
				Value[NoFr1+1] = phi[1];

				Value[NoFr1+2] = phi[4];
				Value[NoFr1+3] = phi[5];
				Value[NoFr1+4] = phi[6];
				Value[NoFr1+5] = phi[7];

				indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is)+1;

				indexValue[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is)+1;
				indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is)+3;
				indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is)+4;

				/* For Fr2 */

				Value[NoFr2] = psi[0];
				Value[NoFr2+1] = psi[1];

				Value[NoFr2+2] = psi[4];
				Value[NoFr2+3] = psi[5];
				Value[NoFr2+4] = psi[6];
				Value[NoFr2+5] = psi[7];

				indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

				indexValue[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is)+2;
				indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is)+3;
				indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is)+5;


				if(ix1 != 4){
					
					/* For Er */
					
					Value[NoEr+7] = theta[9];
					Value[NoEr+8] = theta[10];

					
					
					indexValue[NoEr+7] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoEr+8] = 3*(j-js+1)*Nx + 3*(i-is)+2;


					/* For Fr1 */
					
					Value[NoFr1+6] = phi[8];
					Value[NoFr1+7] = phi[9];					
										
					indexValue[NoFr1+6] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr1+7] = 3*(j-js+1)*Nx + 3*(i-is)+1;

					/* For Fr2 */
		
					Value[NoFr2+6] = psi[8];
					Value[NoFr2+7] = psi[9];					
										
					indexValue[NoFr2+6] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr2+7] = 3*(j-js+1)*Nx + 3*(i-is)+2;			


				}/* no periodic boundary condition */
				else {
					/* For Er */
					
					Value[NoEr+7] = theta[2];
					Value[NoEr+8] = theta[3];
					Value[NoEr+9] = theta[9];
					Value[NoEr+10] = theta[10];
					
					
					indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(ie-is);
					indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(ie-is)+1;
					indexValue[NoEr+9] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoEr+10] = 3*(j-js+1)*Nx + 3*(i-is)+2;



					/* For Fr1 */
					
					Value[NoFr1+6] = phi[2];
					Value[NoFr1+7] = phi[3];	
					Value[NoFr1+8] = phi[8];
					Value[NoFr1+9] = phi[9];						
										
					indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(ie-is);
					indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(ie-is)+1;
					indexValue[NoFr1+8] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr1+9] = 3*(j-js+1)*Nx + 3*(i-is)+1;

					/* For Fr2 */
		
					Value[NoFr2+6] = psi[2];
					Value[NoFr2+7] = psi[3];
					Value[NoFr2+8] = psi[8];
					Value[NoFr2+9] = psi[9];						
										
					indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(ie-is);
					indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(ie-is)+2;
					indexValue[NoFr2+8] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr2+9] = 3*(j-js+1)*Nx + 3*(i-is)+2;


				}/* For periodic boundary condition */



				
				/* other ix1 boundary condition */

				if(ix1 == 1 || ix1 == 5){
					
					Value[NoEr+2] += theta[2];
					Value[NoEr+3] -= theta[3];
			
					Value[NoFr1+2]+= phi[2];
					Value[NoFr1+3]-= phi[3];
				
					Value[NoFr2+2]+= psi[2];
					Value[NoFr2+3]+= psi[3];
				}
				else if(ix1 == 2){
					Value[NoEr+2] += theta[2];
					Value[NoEr+3] += theta[3];
			
					Value[NoFr1+2]+= phi[2];
					Value[NoFr1+3]+= phi[3];
				
					Value[NoFr2+2]+= psi[2];
					Value[NoFr2+3]+= psi[3];
					}
				else if(ix1 == 3){

					/* Do nothing */
				}
				else {
					goto on_error;
				}

			} /* End j!= js & j != je */
		}/* End i==is */
		else if (i == ie){
			if(j == js){
				if(ox1 != 4){						
					if(ix2 !=4){
						NZ_NUM -= 12;
						NoEr = count;
						NoFr1 = count + 7;
						NoFr2 = count + 13;
						count += 19;
					}
					else{
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
				}/* Non periodic for x1 */
				else{
					if(ix2 !=4){
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
					else{
						NoEr = count;
						NoFr1 = count + 11;
						NoFr2 = count + 21;
						count += 31;
					}
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;


				if(ox1 !=4 ){
	
					/* For Er */
					for(m=0; m<5; m++)
						Value[NoEr+m] = theta[2+m];

					Value[NoEr+5] = theta[9];
					Value[NoEr+6] = theta[10];
				
					indexValue[NoEr] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoEr+1] = 3*(j-js)*Nx+3*(i-is-1)+1;

					for(m=0; m<3; m++)
						indexValue[NoEr+2+m] = 3*(j-js)*Nx+3*(i-is)+m;

					indexValue[NoEr+5] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoEr+6] = 3*(j-js+1)*Nx+3*(i-is)+2;

					/* For Fr1 */

					for(m=0; m<4; m++)
						Value[NoFr1+m] = phi[2+m];

					Value[NoFr1+4] = phi[8];
					Value[NoFr1+5] = phi[9];
				
					indexValue[NoFr1] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoFr1+1] = 3*(j-js)*Nx+3*(i-is-1)+1;

					for(m=0; m<2; m++)
						indexValue[NoFr1+2+m] = 3*(j-js)*Nx+3*(i-is)+m;

					indexValue[NoFr1+4] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoFr1+5] = 3*(j-js+1)*Nx+3*(i-is)+1;

					/* For Fr2 */

					for(m=0; m<4; m++)
						Value[NoFr2+m] = psi[2+m];

					Value[NoFr2+4] = psi[8];
					Value[NoFr2+5] = psi[9];
				
					indexValue[NoFr2] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoFr2+1] = 3*(j-js)*Nx+3*(i-is-1)+2;
					indexValue[NoFr2+2] = 3*(j-js)*Nx+3*(i-is);
					indexValue[NoFr2+3] = 3*(j-js)*Nx+3*(i-is)+2;					

					indexValue[NoFr2+4] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoFr2+5] = 3*(j-js+1)*Nx+3*(i-is)+2;		

					/* ix2 boundary condition */
					if(ix2 == 4){
	
						Value[NoEr+7] = theta[0];
						Value[NoEr+8] = theta[1];
					
						indexValue[NoEr+7] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoEr+8] = 3*(je-js)*Nx+3*(i-is)+2;
			
						Value[NoFr1+6] = phi[0];
						Value[NoFr1+7] = phi[1];

						indexValue[NoFr1+6] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoFr1+7] = 3*(je-js)*Nx+3*(i-is)+1;
				
						Value[NoFr2+6] = psi[0];
						Value[NoFr2+7] = psi[1];

						indexValue[NoFr2+6] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoFr2+7] = 3*(je-js)*Nx+3*(i-is)+2;				


					}
					else if(ix2 == 1 || ix2 == 5){
					
						Value[NoEr+2] += theta[0];
						Value[NoEr+4] -= theta[1];
			
						Value[NoFr1+2]+= phi[0];
						Value[NoFr1+3]+= phi[1];
				
						Value[NoFr2+2]+= psi[0];
						Value[NoFr2+3]-= psi[1];
					}
					else if(ix2 == 2){
						Value[NoEr+2] += theta[0];
						Value[NoEr+4] += theta[1];
			
						Value[NoFr1+2]+= phi[0];
						Value[NoFr1+3]+= phi[1];
				
						Value[NoFr2+2]+= psi[0];
						Value[NoFr2+3]+= psi[1];
					}
					else if(ix2 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}



				}/* Non periodic boundary condition */
				else {
					/* For Er */
					Value[NoEr] = theta[7];
					Value[NoEr+1] = theta[8];

					for(m=0; m<5; m++)
						Value[NoEr+m+2] = theta[2+m];

					Value[NoEr+7] = theta[9];
					Value[NoEr+8] = theta[10];
				
					indexValue[NoEr] = 3*(j-js)*Nx+3*(ie-is);
					indexValue[NoEr+1] = 3*(j-js)*Nx+3*(ie-is)+1;
					indexValue[NoEr+2] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoEr+3] = 3*(j-js)*Nx+3*(i-is-1)+1;

					for(m=0; m<3; m++)
						indexValue[NoEr+4+m] = 3*(j-js)*Nx+3*(i-is)+m;

					indexValue[NoEr+7] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoEr+8] = 3*(j-js+1)*Nx+3*(i-is)+2;

					/* For Fr1 */
					Value[NoFr1+1] = phi[6];
					Value[NoFr1+2] = phi[7];

					for(m=0; m<4; m++)
						Value[NoFr1+m+2] = phi[2+m];

					Value[NoFr1+6] = phi[8];
					Value[NoFr1+7] = phi[9];
				
					indexValue[NoFr1] = 3*(j-js)*Nx+3*(ie-is);
					indexValue[NoFr1+1] = 3*(j-js)*Nx+3*(ie-is)+1;
					indexValue[NoFr1+2] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoFr1+3] = 3*(j-js)*Nx+3*(i-is-1)+1;

					for(m=0; m<2; m++)
						indexValue[NoFr1+4+m] = 3*(j-js)*Nx+3*(i-is)+m;

					indexValue[NoFr1+6] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoFr1+7] = 3*(j-js+1)*Nx+3*(i-is)+1;

					/* For Fr2 */

					Value[NoFr2+1] = psi[6];
					Value[NoFr2+2] = psi[7];

					for(m=0; m<4; m++)
						Value[NoFr2+m+2] = psi[2+m];

					Value[NoFr2+6] = psi[8];
					Value[NoFr2+7] = psi[9];
				
					indexValue[NoFr2] = 3*(j-js)*Nx+3*(ie-is);
					indexValue[NoFr2+1] = 3*(j-js)*Nx+3*(ie-is)+2;
					indexValue[NoFr2+2] = 3*(j-js)*Nx+3*(i-is-1);
					indexValue[NoFr2+3] = 3*(j-js)*Nx+3*(i-is-1)+2;
					indexValue[NoFr2+4] = 3*(j-js)*Nx+3*(i-is);
					indexValue[NoFr2+5] = 3*(j-js)*Nx+3*(i-is)+2;			

					indexValue[NoFr2+6] = 3*(j-js+1)*Nx+3*(i-is);
					indexValue[NoFr2+7] = 3*(j-js+1)*Nx+3*(i-is)+2;


					/* ix2 boundary condition */
					if(ix2 == 4){
	
						Value[NoEr+9] = theta[0];
						Value[NoEr+10] = theta[1];
					
						indexValue[NoEr+9] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoEr+10] = 3*(je-js)*Nx+3*(i-is)+2;
			
						Value[NoFr1+8] = phi[0];
						Value[NoFr1+9] = phi[1];

						indexValue[NoEr+8] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoEr+9] = 3*(je-js)*Nx+3*(i-is)+1;
				
						Value[NoFr2+8] = psi[0];
						Value[NoFr2+9] = psi[1];

						indexValue[NoEr+8] = 3*(je-js)*Nx+3*(i-is);
						indexValue[NoEr+9] = 3*(je-js)*Nx+3*(i-is)+2;				


					}
					else if(ix2 == 1 || ix2 == 5){
					
						Value[NoEr+4] += theta[0];
						Value[NoEr+6] -= theta[1];
			
						Value[NoFr1+4]+= phi[0];
						Value[NoFr1+5]+= phi[1];
				
						Value[NoFr2+4]+= psi[0];
						Value[NoFr2+5]-= psi[1];
					}
					else if(ix2 == 2){
						Value[NoEr+4] += theta[0];
						Value[NoEr+6] += theta[1];
			
						Value[NoFr1+4]+= phi[0];
						Value[NoFr1+5]+= phi[1];
				
						Value[NoFr2+4]+= psi[0];
						Value[NoFr2+5]+= psi[1];
					}
					else if(ix2 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}

					
				}/* Periodic boundary condition */


				/* other ox1 boundary condition */
				if(ox1 == 1 || ox1 == 5){
					
					Value[NoEr+2] += theta[7];
					Value[NoEr+3] -= theta[8];
			
					Value[NoFr1+2]+= phi[6];
					Value[NoFr1+3]-= phi[7];
				
					Value[NoFr2+2]+= psi[6];
					Value[NoFr2+3]+= psi[7];
				}
				else if(ox1 == 2){
					Value[NoEr+2] += theta[7];
					Value[NoEr+3] += theta[8];
			
					Value[NoFr1+2]+= phi[6];
					Value[NoFr1+3]+= phi[7];
				
					Value[NoFr2+2]+= psi[6];
					Value[NoFr2+3]+= psi[7];
				}
				else if(ox1 == 3){

					/* Do nothing */
				}
				else {
					goto on_error;
				}
			} /* End j==js */
			else if(j == je){
				if(ox1 != 4){						
					if(ox2 !=4){
						NZ_NUM -= 12;
						NoEr = count;
						NoFr1 = count + 7;
						NoFr2 = count + 13;
						count += 19;
					}
					else{
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
				}/* Non periodic for x1 */
				else{
					if(ox2 !=4){
						NZ_NUM -= 6;
						NoEr = count;
						NoFr1 = count + 9;
						NoFr2 = count + 17;
						count += 25;
					}
					else{
						NoEr = count;
						NoFr1 = count + 11;
						NoFr2 = count + 21;
						count += 31;
					}
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;
				
				
				/* Now the important thing is ox2, which determines the first non-zero element */
				
				
				if(ox2 != 4){
					if(ox1 != 4){					
						/* For Er */
						for(m=0; m<7; m++)
							Value[NoEr+m] = theta[m];
				
						indexValue[NoEr] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoEr+1] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoEr+2] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoEr+3] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoEr+4] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoEr+5] = 3*(j-js)*Nx+3*(i-is)+1;
						indexValue[NoEr+6] = 3*(j-js)*Nx+3*(i-is)+2;				

						/* For Fr1 */

						for(m=0; m<6; m++)
							Value[NoFr1+m] = phi[m];
				
						indexValue[NoFr1] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr1+1] = 3*(j-js-1)*Nx+3*(i-is)+1;

						indexValue[NoFr1+2] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr1+3] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoFr1+4] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr1+5] = 3*(j-js)*Nx+3*(i-is)+1;
						
						/* For Fr2 */

						for(m=0; m<6; m++)
							Value[NoFr2+m] = psi[m];
				
						indexValue[NoFr2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr2+1] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoFr2+2] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr2+3] = 3*(j-js)*Nx+3*(i-is-1)+2;

						indexValue[NoFr2+4] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr2+5] = 3*(j-js)*Nx+3*(i-is)+2;


						/* other ox1 boundary condition */
						if(ox1 == 1 || ox1 == 5){
					
							Value[NoEr+4] += theta[7];
							Value[NoEr+5] -= theta[8];
			
							Value[NoFr1+4]+= phi[6];
							Value[NoFr1+5]-= phi[7];
				
							Value[NoFr2+4]+= psi[6];
							Value[NoFr2+5]+= psi[7];
						}
						else if(ox1 == 2){
							Value[NoEr+4] += theta[7];
							Value[NoEr+5] += theta[8];
			
							Value[NoFr1+4]+= phi[6];
							Value[NoFr1+5]+= phi[7];
				
							Value[NoFr2+4]+= psi[6];
							Value[NoFr2+5]+= psi[7];
						}
						else if(ox1 == 3){

							/* Do nothing */
						}
						else {
							goto on_error;
						}

						/* other x2 boundary condition */
						if(ox2 == 1 || ox2 == 5){
					
							Value[NoEr+4] += theta[9];
							Value[NoEr+6] -= theta[10];
			
							Value[NoFr1+4]+= phi[8];
							Value[NoFr1+5]+= phi[9];
				
							Value[NoFr2+4]+= psi[8];
							Value[NoFr2+5]-= psi[9];
						}
						else if(ox2 == 2){
							Value[NoEr+4] += theta[9];
							Value[NoEr+6] += theta[10];
			
							Value[NoFr1+4]+= phi[8];
							Value[NoFr1+5]+= phi[9];
				
							Value[NoFr2+4]+= psi[9];
							Value[NoFr2+5]+= psi[10];
						}
						else if(ox2 == 3){

							/* Do nothing */
						}
						else {
							goto on_error;
						}


					}/* non periodic for x1 */
					else{
						/* for Er */
						Value[NoEr] = theta[0];
						Value[NoEr+1] = theta[1];

						Value[NoEr+2] = theta[7];
						Value[NoEr+3] = theta[8];

						for(m=0; m<5; m++)
							Value[NoEr+4+m] = theta[2+m];
				
						indexValue[NoEr] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoEr+1] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoEr+2] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoEr+3] = 3*(j-js)*Nx+3*(ie-is)+1;

						indexValue[NoEr+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoEr+5] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoEr+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoEr+7] = 3*(j-js)*Nx+3*(i-is)+1;
						indexValue[NoEr+8] = 3*(j-js)*Nx+3*(i-is)+2;				

						/* For Fr1 */

						Value[NoFr1] = phi[0];
						Value[NoFr1+1] = phi[1];

						Value[NoFr1+2] = phi[6];
						Value[NoFr1+3] = phi[7];

						for(m=0; m<4; m++)
							Value[NoFr1+4+m] = phi[2+m];
				
						indexValue[NoFr1] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr1+1] = 3*(j-js-1)*Nx+3*(i-is)+1;

						indexValue[NoFr1+2] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoFr1+3] = 3*(j-js)*Nx+3*(ie-is)+1;

						indexValue[NoFr1+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr1+5] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoFr1+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr1+7] = 3*(j-js)*Nx+3*(i-is)+1;
						
						
						/* For Fr2 */
						Value[NoFr2] = psi[0];
						Value[NoFr2+1] = psi[1];

						Value[NoFr2+2] = psi[6];
						Value[NoFr2+3] = psi[7];

						for(m=0; m<4; m++)
							Value[NoFr2+4+m] = psi[2+m];
				
						indexValue[NoFr2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr2+1] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoFr2+2] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoFr2+3] = 3*(j-js)*Nx+3*(ie-is)+2;

						indexValue[NoFr2+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr2+5] = 3*(j-js)*Nx+3*(i-is-1)+2;

						indexValue[NoFr2+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr2+7] = 3*(j-js)*Nx+3*(i-is)+2;

						/* other x2 boundary condition */
						if(ox2 == 1 || ox2 == 5){
					
							Value[NoEr+6] += theta[9];
							Value[NoEr+8] -= theta[10];
			
							Value[NoFr1+6]+= phi[8];
							Value[NoFr1+7]+= phi[9];
				
							Value[NoFr2+6]+= psi[8];
							Value[NoFr2+7]-= psi[9];
						}
						else if(ox2 == 2){
							Value[NoEr+6] += theta[9];
							Value[NoEr+8] += theta[10];
			
							Value[NoFr1+6]+= phi[8];
							Value[NoFr1+7]+= phi[9];
				
							Value[NoFr2+6]+= psi[8];
							Value[NoFr2+7]+= psi[9];
						}
						else if(ox2 == 3){

							/* Do nothing */
						}
						else {
							goto on_error;
						}


					}/* periodic for x1 */				

				}/* Non-periodic for x2 */
				else{
					if(ox1 != 4){					
						/* For Er */
						Value[NoEr] = theta[9];
						Value[NoEr+1] = theta[10];

						for(m=0; m<7; m++)
							Value[NoEr+m+2] = theta[m];

						indexValue[NoEr] = 3*(i-is);
						indexValue[NoEr+1] = 3*(i-is)+2;
				
						indexValue[NoEr+2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoEr+3] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoEr+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoEr+5] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoEr+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoEr+7] = 3*(j-js)*Nx+3*(i-is)+1;
						indexValue[NoEr+8] = 3*(j-js)*Nx+3*(i-is)+2;				

						/* For Fr1 */

						Value[NoFr1] = phi[8];
						Value[NoFr1+1] = phi[9];

						for(m=0; m<6; m++)
							Value[NoFr1+m+2] = phi[m];

						indexValue[NoFr1] = 3*(i-is);
						indexValue[NoFr1+1] = 3*(i-is)+1;
				
						indexValue[NoFr1+2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr1+3] = 3*(j-js-1)*Nx+3*(i-is)+1;

						indexValue[NoFr1+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr1+5] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoFr1+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr1+7] = 3*(j-js)*Nx+3*(i-is)+1;
						
						/* For Fr2 */
	
						Value[NoFr2] = psi[8];
						Value[NoFr2+1] = psi[9];

						for(m=0; m<6; m++)
							Value[NoFr2+m+2] = psi[m];

						indexValue[NoFr2] = 3*(i-is);
						indexValue[NoFr2+1] = 3*(i-is)+2;
				
						indexValue[NoFr2+2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr2+3] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoFr2+4] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr2+5] = 3*(j-js)*Nx+3*(i-is-1)+2;

						indexValue[NoFr2+6] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr2+7] = 3*(j-js)*Nx+3*(i-is)+2;


						/* other ox1 boundary condition */
						if(ox1 == 1 || ox1 == 5){
					
							Value[NoEr+6] += theta[7];
							Value[NoEr+7] -= theta[8];
			
							Value[NoFr1+6]+= phi[6];
							Value[NoFr1+7]-= phi[7];
				
							Value[NoFr2+6]+= psi[6];
							Value[NoFr2+7]+= psi[7];
						}
						else if(ox1 == 2){
							Value[NoEr+6] += theta[7];
							Value[NoEr+7] += theta[8];
			
							Value[NoFr1+6]+= phi[6];
							Value[NoFr1+7]+= phi[7];
				
							Value[NoFr2+6]+= psi[6];
							Value[NoFr2+7]+= psi[7];
						}
						else if(ox1 == 3){

							/* Do nothing */
						}
						else {
							goto on_error;
						}

						


					}/* non periodic for x1 */
					else{
						/* for Er */
						Value[NoEr] = theta[9];
						Value[NoEr+1] = theta[10];

						Value[NoEr+2] = theta[0];
						Value[NoEr+3] = theta[1];

						Value[NoEr+4] = theta[7];
						Value[NoEr+5] = theta[8];

						for(m=0; m<5; m++)
							Value[NoEr+6+m] = theta[2+m];

						indexValue[NoEr] = 3*(i-is);
						indexValue[NoEr+1] = 3*(i-is)+2;
				
						indexValue[NoEr+1] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoEr+2] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoEr+3] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoEr+4] = 3*(j-js)*Nx+3*(ie-is)+1;

						indexValue[NoEr+5] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoEr+6] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoEr+7] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoEr+8] = 3*(j-js)*Nx+3*(i-is)+1;
						indexValue[NoEr+9] = 3*(j-js)*Nx+3*(i-is)+2;				

						/* For Fr1 */
						Value[NoFr1] = phi[8];
						Value[NoFr1+1] = phi[9];

						Value[NoFr1+2] = phi[0];
						Value[NoFr1+3] = phi[1];

						Value[NoFr1+4] = phi[6];
						Value[NoFr1+5] = phi[7];

						for(m=0; m<4; m++)
							Value[NoFr1+6+m] = phi[2+m];

						indexValue[NoFr1] = 3*(i-is);
						indexValue[NoFr1+1] = 3*(i-is)+1;
				
						indexValue[NoFr1+2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr1+3] = 3*(j-js-1)*Nx+3*(i-is)+1;

						indexValue[NoFr1+4] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoFr1+5] = 3*(j-js)*Nx+3*(ie-is)+1;

						indexValue[NoFr1+6] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr1+7] = 3*(j-js)*Nx+3*(i-is-1)+1;

						indexValue[NoFr1+8] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr1+9] = 3*(j-js)*Nx+3*(i-is)+1;
						
						
						/* For Fr2 */
						Value[NoFr2] = psi[8];
						Value[NoFr2+1] = psi[9];

						Value[NoFr2+2] = psi[0];
						Value[NoFr2+3] = psi[1];

						Value[NoFr2+4] = psi[6];
						Value[NoFr2+5] = psi[7];

						for(m=0; m<4; m++)
							Value[NoFr2+6+m] = psi[2+m];

						indexValue[NoFr2] = 3*(i-is);
						indexValue[NoFr2+1] = 3*(i-is)+2;
				
						indexValue[NoFr2+2] = 3*(j-js-1)*Nx+3*(i-is);
						indexValue[NoFr2+3] = 3*(j-js-1)*Nx+3*(i-is)+2;

						indexValue[NoFr2+4] = 3*(j-js)*Nx+3*(ie-is);
						indexValue[NoFr2+5] = 3*(j-js)*Nx+3*(ie-is)+2;

						indexValue[NoFr2+6] = 3*(j-js)*Nx+3*(i-is-1);
						indexValue[NoFr2+7] = 3*(j-js)*Nx+3*(i-is-1)+2;

						indexValue[NoFr2+8] = 3*(j-js)*Nx+3*(i-is);
						indexValue[NoFr2+9] = 3*(j-js)*Nx+3*(i-is)+2;					

					}/* periodic for x1 */	
				}/* periodic for x2 */
			} /* End j==je */
			else {
				if(ox1 != 4){						
					NZ_NUM -= 6;
					NoEr = count;
					NoFr1 = count + 9;
					NoFr2 = count + 17;
					count += 25;
					
				}/* Non periodic for x1 */
				else{
				
					NoEr = count;
					NoFr1 = count + 11;
					NoFr2 = count + 21;
					count += 31;
				
				}/* periodic for x1 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;

				if(ox1 == 4){

					/* For Er */
					Value[NoEr] = theta[0];
					Value[NoEr+1] =theta[1];

					Value[NoEr+2] = theta[7];
					Value[NoEr+3] = theta[8];

					for(m=0; m<5; m++)
						Value[NoEr+4+m] = theta[2+m];

					Value[NoEr+9] = theta[9];
					Value[NoEr+10] =theta[10];

					indexValue[NoEr] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

					indexValue[NoEr+2] = 3*(j-js)*Nx;
					indexValue[NoEr+3] = 3*(j-js)*Nx +1;

					indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(i-is)+1;
					indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoEr+9] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoEr+10] = 3*(j-js+1)*Nx + 3*(i-is)+2;

					/* For Fr1 */

					Value[NoFr1] = phi[0];
					Value[NoFr1+1] =phi[1];

					Value[NoFr1+2] = phi[6];
					Value[NoFr1+3] = phi[7];

					for(m=0; m<4; m++)
						Value[NoFr1+4+m] = phi[2+m];

					Value[NoFr1+8] = phi[8];
					Value[NoFr1+9] =phi[9];

					indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is)+1;

					indexValue[NoFr1+2] = 3*(j-js)*Nx;
					indexValue[NoFr1+3] = 3*(j-js)*Nx +1;

					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(i-is)+1;
					
					indexValue[NoFr1+8] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr1+9] = 3*(j-js+1)*Nx + 3*(i-is)+1;

					/* For Fr2 */

					Value[NoFr2] = psi[0];
					Value[NoFr2+1] =psi[1];

					Value[NoFr2+2] = psi[6];
					Value[NoFr2+3] = psi[7];

					for(m=0; m<4; m++)
						Value[NoFr2+4+m] = psi[2+m];

					Value[NoFr2+8] = psi[8];
					Value[NoFr2+9] =psi[9];

					indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

					indexValue[NoFr2+2] = 3*(j-js)*Nx;
					indexValue[NoFr2+3] = 3*(j-js)*Nx +2;

					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is-1)+2;

					indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoFr2+8] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr2+9] = 3*(j-js+1)*Nx + 3*(i-is)+2;


				}/* Periodic boundary condition */
				else{
					/* For Er */
					Value[NoEr] = theta[0];
					Value[NoEr+1] =theta[1];

					for(m=0; m<5; m++)
						Value[NoEr+2+m] = theta[2+m];

					Value[NoEr+7] = theta[9];
					Value[NoEr+8] =theta[10];

					indexValue[NoEr] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

					indexValue[NoEr+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoEr+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is)+1;
					indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoEr+7] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoEr+8] = 3*(j-js+1)*Nx + 3*(i-is)+2;

					/* For Fr1 */

					Value[NoFr1] = phi[0];
					Value[NoFr1+1] =phi[1];

					for(m=0; m<4; m++)
						Value[NoFr1+2+m] = phi[2+m];

					Value[NoFr1+6] = phi[8];
					Value[NoFr1+7] =phi[9];

					indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is)+1;

					indexValue[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is)+1;
					
					indexValue[NoFr1+6] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr1+7] = 3*(j-js+1)*Nx + 3*(i-is)+1;

					/* For Fr2 */

					Value[NoFr2] = psi[0];
					Value[NoFr2+1] =psi[1];


					for(m=0; m<4; m++)
						Value[NoFr2+2+m] = psi[2+m];

					Value[NoFr2+6] = psi[8];
					Value[NoFr2+7] =psi[9];

					indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

					indexValue[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is-1)+2;

					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoFr2+6] = 3*(j-js+1)*Nx + 3*(i-is);
					indexValue[NoFr2+7] = 3*(j-js+1)*Nx + 3*(i-is)+2;

					/* other ox1 boundary condition */
					if(ox1 == 1 || ox1 == 5){
					
						Value[NoEr+4] += theta[7];
						Value[NoEr+5] -= theta[8];
			
						Value[NoFr1+4]+= phi[6];
						Value[NoFr1+5]-= phi[7];
				
						Value[NoFr2+4]+= psi[6];
						Value[NoFr2+5]+= psi[7];
					}
					else if(ox1 == 2){
						Value[NoEr+4] += theta[7];
						Value[NoEr+5] += theta[8];
		
						Value[NoFr1+4]+= phi[6];
						Value[NoFr1+5]+= phi[7];
				
						Value[NoFr2+4]+= psi[6];
						Value[NoFr2+5]+= psi[7];
						}
						else if(ox1 == 3){

							/* Do nothing */
						}
						else {
							goto on_error;
						}

				}/* non-periodic boundary condition */
			} /* End j!=js & j!= je*/
		}/* End i==ie */
		else {
			if(j == js){
				if(ix2 != 4){						
					NZ_NUM -= 6;
					NoEr = count;
					NoFr1 = count + 9;
					NoFr2 = count + 17;
					count += 25;
					
				}/* Non periodic for x2 */
				else{
				
					NoEr = count;
					NoFr1 = count + 11;
					NoFr2 = count + 21;
					count += 31;
				
				}/* periodic for x2 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;

			/* The following is true no matter ix2==4 or not */
				/* For Er */
				for(m=0; m<9; m++)
					Value[NoEr+m] = theta[2+m];

				indexValue[NoEr] = 3*(i-is-1);
				indexValue[NoEr+1] = 3*(i-is-1) + 1;
		
				for(m=0; m<5; m++)
					indexValue[NoEr+2+m] = 3*(i-is)+m;

				indexValue[NoEr+7] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoEr+8] = 3*(j-js+1)*Nx + 3*(i-is)+2;

				/* For Fr1 */
				for(m=0; m<8; m++)
					Value[NoFr1+m] = phi[2+m];

				indexValue[NoFr1] = 3*(i-is-1);
				indexValue[NoFr1+1] = 3*(i-is-1) + 1;

				indexValue[NoFr1+2] = 3*(i-is);
				indexValue[NoFr1+3] = 3*(i-is) + 1;

				indexValue[NoFr1+4] = 3*(i-is)+3;
				indexValue[NoFr1+5] = 3*(i-is)+4;
				
				indexValue[NoFr1+6] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoFr1+7] = 3*(j-js+1)*Nx + 3*(i-is)+1;

				/* For Fr2 */

				for(m=0; m<8; m++)
					Value[NoFr2+m] = psi[2+m];

				indexValue[NoFr2] = 3*(i-is-1);
				indexValue[NoFr2+1] = 3*(i-is-1) + 2;

				indexValue[NoFr2+2] = 3*(i-is);
				indexValue[NoFr2+3] = 3*(i-is) + 2;

				indexValue[NoFr2+4] = 3*(i-is)+3;
				indexValue[NoFr2+5] = 3*(i-is)+5;
				
				indexValue[NoFr2+6] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoFr2+7] = 3*(j-js+1)*Nx + 3*(i-is)+2;

				if(ix2 == 4){
	
					Value[NoEr+9] = theta[0];
					Value[NoEr+10] = theta[1];
					
					indexValue[NoEr+9] = 3*(je-js)*Nx+3*(i-is);
					indexValue[NoEr+10] = 3*(je-js)*Nx+3*(i-is)+2;
			
					Value[NoFr1+8] = psi[0];
					Value[NoFr1+9] = psi[1];

					indexValue[NoFr1+8] = 3*(je-js)*Nx+3*(i-is);
					indexValue[NoFr1+9] = 3*(je-js)*Nx+3*(i-is)+1;
				
					Value[NoFr2+8] = psi[0];
					Value[NoFr2+9] = psi[1];

					indexValue[NoFr2+8] = 3*(je-js)*Nx+3*(i-is);
					indexValue[NoFr2+9] = 3*(je-js)*Nx+3*(i-is)+2;				


				}
				else if(ix2 == 1 || ix2 == 5){
					
					Value[NoEr+2] += theta[0];
					Value[NoEr+4] -= theta[1];
			
					Value[NoFr1+2]+= phi[0];
					Value[NoFr1+3]+= phi[1];
				
					Value[NoFr2+2]+= psi[0];
					Value[NoFr2+3]-= psi[1];
				}
				else if(ix2 == 2){
					Value[NoEr+2] += theta[0];
					Value[NoEr+4] += theta[1];
			
					Value[NoFr1+2]+= phi[0];
					Value[NoFr1+3]+= phi[1];
				
					Value[NoFr2+2]+= psi[0];
					Value[NoFr2+3]+= psi[1];
				}
				else if(ix2 == 3){

					/* Do nothing */
				}
				else {
					goto on_error;
				}
				
			}
			else if(j == je){

				if(ox2 != 4){						
					NZ_NUM -= 6;
					NoEr = count;
					NoFr1 = count + 9;
					NoFr2 = count + 17;
					count += 25;
					
				}/* Non periodic for x2 */
				else{
				
					NoEr = count;
					NoFr1 = count + 11;
					NoFr2 = count + 21;
					count += 31;
				
				}/* periodic for x2 */

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;

				if(ox2 == 4){
					/* For Er */
					Value[NoEr] = theta[9];
					Value[NoEr+1] = theta[10];
					
					for(m=0; m<9; m++)
						Value[NoEr+2+m] = theta[m];

					indexValue[NoEr] = 3*(i-is);
					indexValue[NoEr+1] = 3*(i-is) + 2;

					indexValue[NoEr+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;		

					indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(i-is)+1;
					indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(i-is)+2;

					indexValue[NoEr+9] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoEr+10] = 3*(j-js)*Nx + 3*(i-is+1)+1;

					/* For Fr1 */

					Value[NoFr1] = phi[8];
					Value[NoFr1+1] = phi[9];
					
					for(m=0; m<8; m++)
						Value[NoFr1+2+m] = phi[m];

					indexValue[NoFr1] = 3*(i-is);
					indexValue[NoFr1+1] = 3*(i-is) + 1;

					indexValue[NoFr1+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+3] = 3*(j-js-1)*Nx + 3*(i-is) + 1;		

					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(i-is)+1;
					
					indexValue[NoFr1+8] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoFr1+9] = 3*(j-js)*Nx + 3*(i-is+1)+1;

					/* For Fr2 */

					Value[NoFr2] = psi[8];
					Value[NoFr2+1] = psi[9];
					
					for(m=0; m<8; m++)
						Value[NoFr2+2+m] = psi[m];

					indexValue[NoFr2] = 3*(i-is);
					indexValue[NoFr2+1] = 3*(i-is) + 2;

					indexValue[NoFr2+2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+3] = 3*(j-js-1)*Nx + 3*(i-is) + 2;		

					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is-1)+2;

					indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoFr2+8] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoFr2+9] = 3*(j-js)*Nx + 3*(i-is+1)+2;
				}/* End periodic of x2 */
				else{
					/* For Er */
					
					for(m=0; m<9; m++)
						Value[NoEr+m] = theta[m];

					indexValue[NoEr+0] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is) + 2;		

					indexValue[NoEr+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoEr+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is)+1;
					indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is)+2;

					indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(i-is+1)+1;

					/* For Fr1 */

					
					for(m=0; m<8; m++)
						Value[NoFr1+m] = phi[m];


					indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is) + 1;		

					indexValue[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

					indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is)+1;
					
					indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(i-is+1)+1;

					/* For Fr2 */

					
					for(m=0; m<8; m++)
						Value[NoFr2+m] = psi[m];

					indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
					indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is) + 2;		

					indexValue[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is-1);
					indexValue[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is-1)+2;

					indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is);
					indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is)+2;
					
					indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(i-is+1);
					indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(i-is+1)+2;

					/* other x2 boundary condition */
					if(ox2 == 1 || ox2 == 5){
					
						Value[NoEr+4] += theta[9];
						Value[NoEr+6] -= theta[10];
			
						Value[NoFr1+4]+= phi[8];
						Value[NoFr1+5]+= phi[9];
				
						Value[NoFr2+4]+= psi[8];
						Value[NoFr2+5]-= psi[9];
					}
					else if(ox2 == 2){
						Value[NoEr+4] += theta[9];
						Value[NoEr+6] += theta[10];
			
						Value[NoFr1+4]+= phi[8];
						Value[NoFr1+5]+= phi[9];
				
						Value[NoFr2+4]+= psi[8];
						Value[NoFr2+5]+= psi[9];
					}
					else if(ox2 == 3){

						/* Do nothing */
					}
					else {
						goto on_error;
					}

				}/* End non-periodic of x2 */
			}
			else {
				/* no boundary cells on either direction */
				/*NZ_NUM;*/
				NoEr = count;
				NoFr1 = count + 11;
				NoFr2 = count + 21;
				count += 31;

				ptr[3*(j-js)*Nx+3*(i-is)] = NoEr;
				ptr[3*(j-js)*Nx+3*(i-is)+1] = NoFr1;
				ptr[3*(j-js)*Nx+3*(i-is)+2] = NoFr2;

				/* for Er */
				for(m=0; m<11; m++)
					Value[NoEr+m] = theta[m];

				indexValue[NoEr] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoEr+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

				indexValue[NoEr+2] = 3*(j-js)*Nx + 3*(i-is-1);
				indexValue[NoEr+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

				indexValue[NoEr+4] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoEr+5] = 3*(j-js)*Nx + 3*(i-is)+1;
				indexValue[NoEr+6] = 3*(j-js)*Nx + 3*(i-is)+2;

				indexValue[NoEr+7] = 3*(j-js)*Nx + 3*(i-is+1);
				indexValue[NoEr+8] = 3*(j-js)*Nx + 3*(i-is+1)+1;

				indexValue[NoEr+9] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoEr+10] = 3*(j-js+1)*Nx + 3*(i-is)+2;

				/* For Fr1 */
	
				for(m=0; m<10; m++)
					Value[NoFr1+m] = phi[m];

				indexValue[NoFr1] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoFr1+1] = 3*(j-js-1)*Nx + 3*(i-is)+1;

				indexValue[NoFr1+2] = 3*(j-js)*Nx + 3*(i-is-1);
				indexValue[NoFr1+3] = 3*(j-js)*Nx + 3*(i-is-1)+1;

				indexValue[NoFr1+4] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoFr1+5] = 3*(j-js)*Nx + 3*(i-is)+1;
				
				indexValue[NoFr1+6] = 3*(j-js)*Nx + 3*(i-is+1);
				indexValue[NoFr1+7] = 3*(j-js)*Nx + 3*(i-is+1)+1;

				indexValue[NoFr1+8] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoFr1+9] = 3*(j-js+1)*Nx + 3*(i-is)+1;

				/* For Fr2 */

				for(m=0; m<10; m++)
					Value[NoFr2+m] = psi[m];

				indexValue[NoFr2] = 3*(j-js-1)*Nx + 3*(i-is);
				indexValue[NoFr2+1] = 3*(j-js-1)*Nx + 3*(i-is)+2;

				indexValue[NoFr2+2] = 3*(j-js)*Nx + 3*(i-is-1);
				indexValue[NoFr2+3] = 3*(j-js)*Nx + 3*(i-is-1)+2;

				indexValue[NoFr2+4] = 3*(j-js)*Nx + 3*(i-is);
				indexValue[NoFr2+5] = 3*(j-js)*Nx + 3*(i-is)+2;
				
				indexValue[NoFr2+6] = 3*(j-js)*Nx + 3*(i-is+1);
				indexValue[NoFr2+7] = 3*(j-js)*Nx + 3*(i-is+1)+2;

				indexValue[NoFr2+8] = 3*(j-js+1)*Nx + 3*(i-is);
				indexValue[NoFr2+9] = 3*(j-js+1)*Nx + 3*(i-is)+2;
				
			}
		}/* End i!=is & i!=ie */
		}/* End loop i */
	}/* End loop j */

		/* The last element of ptr */
		ptr[lines] = NZ_NUM;


		/* Assemble the matrix and solve the matrix */
		lis_matrix_set_crs(NZ_NUM,ptr,indexValue,Value,Euler);
		lis_matrix_assemble(Euler);

		
		lis_solver_set_option("-i gmres -p none",solver);
		lis_solver_set_option("-tol 1.0e-12",solver);
		lis_solve(Euler,RHSEuler,INIguess,solver);
		
		/* check the iteration step to make sure 1.0e-12 is reached */
		lis_solver_get_iters(solver,&Matrixiter);

		ath_pout(0,"Matrix Iteration steps: %d\n",Matrixiter);

		
	/* update the radiation quantities in the mesh */	
	for(j=js;j<=je;j++){
		for(i=is; i<=ie; i++){

		lis_vector_get_value(INIguess,3*(j-js)*Nx + 3*(i-is),&(pG->U[ks][j][i].Er));
		lis_vector_get_value(INIguess,3*(j-js)*Nx + 3*(i-is)+1,&(pG->U[ks][j][i].Fr1));
		lis_vector_get_value(INIguess,3*(j-js)*Nx + 3*(i-is)+2,&(pG->U[ks][j][i].Fr2));

		if(pG->U[ks][j][i].Er < 0.0)
			fprintf(stderr,"[BackEuler_2d]: Negative Radiation energy: %e\n",pG->U[ks][j][i].Er);
		}
	
		
	}
	/* Eddington factor is updated in the integrator  */


/* Update the ghost zones for different boundary condition to be used later */
		bvals_radMHD(pM);


/*-----------Finish---------------------*/

  
	/* Free the temporary variables just used for this grid calculation*/
	/* This is done after main loop */
/*	rad_hydro_destruct_2d(Nmatrix);

*/
	
	lis_finalize();
	
  return;	


	on_error:
	
	BackEuler_destruct_2d(31*Nmatrix);
	ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}



/*-------------------------------------------------------------------------*/
/* BackEuler_init_2d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_2d(): function to free memory */
void BackEuler_init_2d(int Nelements)
{

/* Nelements = 31*(je-js+1)*(ie-is+1) */
	int lines;
	lines = 3*Nelements/31;

/* The matrix Euler is stored as a compact form.  */
	/*FOR LIS LIBRARY */
	lis_matrix_create(0,&Euler);	
	lis_matrix_set_size(Euler,0,lines);

	lis_vector_duplicate(Euler,&RHSEuler);
	lis_vector_duplicate(Euler,&INIguess);

	lis_solver_create(&solver);

	/* Allocate value and index array */
/*
	if ((Value = (Real*)malloc(Nelements*sizeof(Real))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");
*/
	Value = (LIS_SCALAR *)malloc(Nelements*sizeof(LIS_SCALAR));

	if ((indexValue = (int*)malloc(Nelements*sizeof(int))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	if ((ptr = (int*)malloc((lines+1)*sizeof(int))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");


	
	return;

}


void BackEuler_destruct_2d()
{

	lis_matrix_destroy(Euler);
	lis_solver_destroy(solver);	
	lis_vector_destroy(RHSEuler);
	lis_vector_destroy(INIguess);
	
	if(Value != NULL) free(Value);
	if(indexValue != NULL) free(indexValue);
	if(ptr != NULL) free(ptr);
}


#endif /* radMHD_INTEGRATOR */
