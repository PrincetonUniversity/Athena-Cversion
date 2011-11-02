#include "../../copyright.h"
/*==============================================================================
 * FILE: Jacobi3D.c
 *
 * PURPOSE: Use Jacobi method to solve matrix in 3D case. 
 * 
 *
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../prototypes.h"
#include "../../prototypes.h"

#ifdef MATRIX_MULTIGRID 

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */


/* Matrix boundary condition functions */
extern void bvals_Matrix(MatrixS *pMat);

/* In Gauss-Seidel scheme, we do not need a temporary array */
/* This is not exactly the original Gauss-Seidel scheme *
 * For ghost zones, we do not use the updated version */
/* So differenet CPUs do not need wait for each other to finish */

void GaussSeidel3D(MatrixS *pMat)
{

	int i, j, k, n;
	int is, ie, js, je, ks, ke;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;
	Real hdtodx2 = 0.5 * pMat->dt/pMat->dx2;
	Real hdtodx3 = 0.5 * pMat->dt/pMat->dx3;
	Real dt = pMat->dt;

	Real tempEr1, tempEr2, tempEr3;
	Real tempFr1, tempFr2, tempFr3, temp0;

	/* To store the coefficient */
	Real ****theta = NULL;
	Real ****phi = NULL;
	Real ****psi = NULL;
	Real ****varphi = NULL;

	if((theta = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,16,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	if((phi = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,16,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");
	
	if((psi = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,16,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	if((varphi = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,16,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");



	/* Temporary variables to setup the matrix */
	Real velocity_x, velocity_y, velocity_z, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1, Cj0, Cj1, Ck0, Ck1;

			
	/* First, Update the boundary cells */
	/* We do not set ghost zones after prolongation */
	/* velocity and T4 in the ghost zones are never used */
	/* We only need Er and Fr in the ghost zones */
	/* Use seperate loop to avoid if in next loop */

	/* Only need to calculate the coefficient once */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				
			velocity_x = pMat->U[k][j][i].V1;
			velocity_y = pMat->U[k][j][i].V2;
			velocity_z = pMat->U[k][j][i].V3;
			T4 = pMat->U[k][j][i].T4;
				
			/* Assuming the velocity is already the original velocity in case of FARGO */			
				
			Sigma_sF = pMat->U[k][j][i].Sigma[0];
			Sigma_aF = pMat->U[k][j][i].Sigma[1];
			Sigma_aP = pMat->U[k][j][i].Sigma[2];
			Sigma_aE = pMat->U[k][j][i].Sigma[3];

			Ci0 = (sqrt(pMat->U[k][j][i].Edd_11) - sqrt(pMat->U[k][j][i-1].Edd_11)) 
				/ (sqrt(pMat->U[k][j][i].Edd_11) + sqrt(pMat->U[k][j][i-1].Edd_11));
			Ci1 =  (sqrt(pMat->U[k][j][i+1].Edd_11) - sqrt(pMat->U[k][j][i].Edd_11)) 
				/ (sqrt(pMat->U[k][j][i+1].Edd_11) + sqrt(pMat->U[k][j][i].Edd_11));
			Cj0 = (sqrt(pMat->U[k][j][i].Edd_22) - sqrt(pMat->U[k][j-1][i].Edd_22)) 
				/ (sqrt(pMat->U[k][j][i].Edd_22) + sqrt(pMat->U[k][j-1][i].Edd_22));
			Cj1 =  (sqrt(pMat->U[k][j+1][i].Edd_22) - sqrt(pMat->U[k][j][i].Edd_22)) 
				/ (sqrt(pMat->U[k][j+1][i].Edd_22) + sqrt(pMat->U[k][j][i].Edd_22));
			Ck0 = (sqrt(pMat->U[k][j][i].Edd_33) - sqrt(pMat->U[k-1][j][i].Edd_33)) 
				/ (sqrt(pMat->U[k][j][i].Edd_33) + sqrt(pMat->U[k-1][j][i].Edd_33));
			Ck1 =  (sqrt(pMat->U[k+1][j][i].Edd_33) - sqrt(pMat->U[k][j][i].Edd_33)) 
				/ (sqrt(pMat->U[k+1][j][i].Edd_33) + sqrt(pMat->U[k][j][i].Edd_33));
			theta[k][j][i][0] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			theta[k][j][i][1] = -Crat * hdtodx3 * (1.0 + Ck0);
			theta[k][j][i][2] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			theta[k][j][i][3] = -Crat * hdtodx2 * (1.0 + Cj0);
			theta[k][j][i][4] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			theta[k][j][i][5] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[k][j][i][6] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				+ Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22)
				+ Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33);

/*				+ Crat * pMat->dt * Sigma_aE;
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_11) * velocity_x 
				+ velocity_y * pMat->U[k][j][i].Edd_21 + velocity_z * pMat->U[k][j][i].Edd_31) * velocity_x / Crat
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_22) * velocity_y 
				+ velocity_x * pMat->U[k][j][i].Edd_21 + velocity_z * pMat->U[k][j][i].Edd_32) * velocity_y / Crat
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_33) * velocity_z 
				+ velocity_x * pMat->U[k][j][i].Edd_31 + velocity_y * pMat->U[k][j][i].Edd_32) * velocity_z / Crat;
*/
			theta[k][j][i][7] = Crat * hdtodx1 * (Ci0 + Ci1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x;*/
			theta[k][j][i][8] = Crat * hdtodx2 * (Cj0 + Cj1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_y;*/
			theta[k][j][i][9] = Crat * hdtodx3 * (Ck0 + Ck1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_z;*/
			theta[k][j][i][10] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			theta[k][j][i][11] = Crat * hdtodx1 * (1.0 - Ci1);
			theta[k][j][i][12] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			theta[k][j][i][13] = Crat * hdtodx2 * (1.0 - Cj1);
			theta[k][j][i][14] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);
			theta[k][j][i][15] = Crat * hdtodx3 * (1.0 - Ck1);
			
			
			phi[k][j][i][0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_31;
			phi[k][j][i][1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			phi[k][j][i][2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_21;
			phi[k][j][i][3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			phi[k][j][i][4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_11;
			phi[k][j][i][5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			phi[k][j][i][6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_11
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_21   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_31   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_11) * velocity_x + pMat->U[k][j][i].Edd_21 * velocity_y + pMat->U[k][j][i].Edd_31 * velocity_z) 
			       + pMat->dt * Sigma_aE * velocity_x;
			phi[k][j][i][7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			phi[k][j][i][8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_11;
			phi[k][j][i][9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			phi[k][j][i][10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_21;
			phi[k][j][i][11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			phi[k][j][i][12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_31;
			phi[k][j][i][13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);



			psi[k][j][i][0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_32;
			psi[k][j][i][1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			psi[k][j][i][2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_22;
			psi[k][j][i][3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			psi[k][j][i][4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_21;
			psi[k][j][i][5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			psi[k][j][i][6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_21
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_22   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_32   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_22) * velocity_y + pMat->U[k][j][i].Edd_21 * velocity_x + pMat->U[k][j][i].Edd_32 * velocity_z) 
			       + pMat->dt * Sigma_aE * velocity_y;
			psi[k][j][i][7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			psi[k][j][i][8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_21;
			psi[k][j][i][9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			psi[k][j][i][10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_22;
			psi[k][j][i][11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			psi[k][j][i][12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_32;
			psi[k][j][i][13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);

			varphi[k][j][i][0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_33;
			varphi[k][j][i][1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			varphi[k][j][i][2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_32;
			varphi[k][j][i][3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			varphi[k][j][i][4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_31;
			varphi[k][j][i][5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			varphi[k][j][i][6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_31
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_32   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_33   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_33) * velocity_z + pMat->U[k][j][i].Edd_31 * velocity_x + pMat->U[k][j][i].Edd_32 * velocity_y) 
			       + pMat->dt * Sigma_aE * velocity_z;
			varphi[k][j][i][7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			varphi[k][j][i][8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_31;
			varphi[k][j][i][9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			varphi[k][j][i][10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_32;
			varphi[k][j][i][11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			varphi[k][j][i][12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_33;
			varphi[k][j][i][13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);

	}


/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
		/* The right hand side is stored in pMat */
		/* The coefficients are calculated according to the formula */

			/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
			/* For Er */
			pMat->U[k][j][i].Er  = pMat->RHS[k][j][i][0];

			tempEr3 = theta[k][j][i][0] * pMat->U[k-1][j][i].Er + theta[k][j][i][14] * pMat->U[k+1][j][i].Er;
			tempEr2 = theta[k][j][i][2] * pMat->U[k][j-1][i].Er + theta[k][j][i][12] * pMat->U[k][j+1][i].Er;
			tempEr1 = theta[k][j][i][4] * pMat->U[k][j][i-1].Er + theta[k][j][i][10] * pMat->U[k][j][i+1].Er;

			tempFr3 = theta[k][j][i][1] * pMat->U[k-1][j][i].Fr3 + theta[k][j][i][15] * pMat->U[k+1][j][i].Fr3;
			tempFr2 = theta[k][j][i][3] * pMat->U[k][j-1][i].Fr2 + theta[k][j][i][13] * pMat->U[k][j+1][i].Fr2;
			tempFr1 = theta[k][j][i][5] * pMat->U[k][j][i-1].Fr1 + theta[k][j][i][11] * pMat->U[k][j][i+1].Fr1;

			temp0 = theta[k][j][i][7] * pMat->U[k][j][i].Fr1 + theta[k][j][i][8] * pMat->U[k][j][i].Fr2 + theta[k][j][i][9] * pMat->U[k][j][i].Fr3;

			pMat->U[k][j][i].Er -= ((tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

			/* diagonal elements are not included */

			pMat->U[k][j][i].Er /= theta[k][j][i][6];

			/*****************************************************/
			/* For Fr1 */

			pMat->U[k][j][i].Fr1  = pMat->RHS[k][j][i][1];

			tempEr3 = phi[k][j][i][0] * pMat->U[k-1][j][i].Er + phi[k][j][i][12] * pMat->U[k+1][j][i].Er;
			tempEr2 = phi[k][j][i][2] * pMat->U[k][j-1][i].Er + phi[k][j][i][10] * pMat->U[k][j+1][i].Er;
			tempEr1 = phi[k][j][i][4] * pMat->U[k][j][i-1].Er + phi[k][j][i][8] * pMat->U[k][j][i+1].Er;

			tempFr3 = phi[k][j][i][1] * pMat->U[k-1][j][i].Fr1 + phi[k][j][i][13] * pMat->U[k+1][j][i].Fr1;
			tempFr2 = phi[k][j][i][3] * pMat->U[k][j-1][i].Fr1 + phi[k][j][i][11] * pMat->U[k][j+1][i].Fr1;
			tempFr1 = phi[k][j][i][5] * pMat->U[k][j][i-1].Fr1 + phi[k][j][i][9] * pMat->U[k][j][i+1].Fr1;

			temp0 = phi[k][j][i][6] * pMat->U[k][j][i].Er;

			
			pMat->U[k][j][i].Fr1 -= ((tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

			pMat->U[k][j][i].Fr1 /= phi[k][j][i][7];

			/**************************************************/

			/* For Fr2 */

			pMat->U[k][j][i].Fr2  = pMat->RHS[k][j][i][2];


			tempEr3 = psi[k][j][i][0] * pMat->U[k-1][j][i].Er + psi[k][j][i][12] * pMat->U[k+1][j][i].Er;
			tempEr2 = psi[k][j][i][2] * pMat->U[k][j-1][i].Er + psi[k][j][i][10] * pMat->U[k][j+1][i].Er;
			tempEr1 = psi[k][j][i][4] * pMat->U[k][j][i-1].Er + psi[k][j][i][8] * pMat->U[k][j][i+1].Er;

			tempFr3 = psi[k][j][i][1] * pMat->U[k-1][j][i].Fr2 + psi[k][j][i][13] * pMat->U[k+1][j][i].Fr2;
			tempFr2 = psi[k][j][i][3] * pMat->U[k][j-1][i].Fr2 + psi[k][j][i][11] * pMat->U[k][j+1][i].Fr2;
			tempFr1 = psi[k][j][i][5] * pMat->U[k][j][i-1].Fr2 + psi[k][j][i][9] * pMat->U[k][j][i+1].Fr2;

			temp0 = psi[k][j][i][6] * pMat->U[k][j][i].Er;

			
			pMat->U[k][j][i].Fr2 -= ((tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

			pMat->U[k][j][i].Fr2 /= psi[k][j][i][7];



			/***************************************************/
			/* For Fr3 */

			pMat->U[k][j][i].Fr3  = pMat->RHS[k][j][i][3];


			tempEr3 = varphi[k][j][i][0] * pMat->U[k-1][j][i].Er + varphi[k][j][i][12] * pMat->U[k+1][j][i].Er;
			tempEr2 = varphi[k][j][i][2] * pMat->U[k][j-1][i].Er + varphi[k][j][i][10] * pMat->U[k][j+1][i].Er;
			tempEr1 = varphi[k][j][i][4] * pMat->U[k][j][i-1].Er + varphi[k][j][i][8] * pMat->U[k][j][i+1].Er;

			tempFr3 = varphi[k][j][i][1] * pMat->U[k-1][j][i].Fr3 + varphi[k][j][i][13] * pMat->U[k+1][j][i].Fr3;
			tempFr2 = varphi[k][j][i][3] * pMat->U[k][j-1][i].Fr3 + varphi[k][j][i][11] * pMat->U[k][j+1][i].Fr3;
			tempFr1 = varphi[k][j][i][5] * pMat->U[k][j][i-1].Fr3 + varphi[k][j][i][9] * pMat->U[k][j][i+1].Fr3;

			temp0 = varphi[k][j][i][6] * pMat->U[k][j][i].Er;

			
			pMat->U[k][j][i].Fr3 -= ((tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

			pMat->U[k][j][i].Fr3 /= varphi[k][j][i][7];

	}
			
		
	/* Update the boundary cells */
	bvals_Matrix(pMat);


}	
  
	if(theta != NULL)
		free_4d_array(theta);

	if(phi != NULL)
		free_4d_array(phi);

	if(psi != NULL)
		free_4d_array(psi);

	if(varphi != NULL)
		free_4d_array(varphi);

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* matrix_multigrid */

