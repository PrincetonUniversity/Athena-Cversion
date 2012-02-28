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

	
	Real tempEr1, tempEr2, tempEr3;
	Real tempFr1, tempFr2, tempFr3, temp0;

	/* damped parameter */
	Real omega, Ert0, Frt0;
	omega = 0.5;

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

			
	/* First, Update the boundary cells */
	/* We do not set ghost zones after prolongation */
	/* velocity and T4 in the ghost zones are never used */
	/* We only need Er and Fr in the ghost zones */
	/* Use seperate loop to avoid if in next loop */

	/* Only need to calculate the coefficient once */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){				
			matrix_coef(pMat, NULL, 3, i, j, k, 0.0, &(theta[k][j][i][0]), &(phi[k][j][i][0]), &(psi[k][j][i][0]), &(varphi[k][j][i][0]));
							
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
			Ert0 = pMat->U[k][j][i].Er;
			pMat->U[k][j][i].Er  = pMat->RHS[k][j][i][0];

			tempEr3 = theta[k][j][i][0] * pMat->U[k-1][j][i].Er + theta[k][j][i][14] * pMat->U[k+1][j][i].Er;
			tempEr2 = theta[k][j][i][2] * pMat->U[k][j-1][i].Er + theta[k][j][i][12] * pMat->U[k][j+1][i].Er;
			tempEr1 = theta[k][j][i][4] * pMat->U[k][j][i-1].Er + theta[k][j][i][10] * pMat->U[k][j][i+1].Er;

			tempFr3 = theta[k][j][i][1] * pMat->U[k-1][j][i].Fr3 + theta[k][j][i][15] * pMat->U[k+1][j][i].Fr3;
			tempFr2 = theta[k][j][i][3] * pMat->U[k][j-1][i].Fr2 + theta[k][j][i][13] * pMat->U[k][j+1][i].Fr2;
			tempFr1 = theta[k][j][i][5] * pMat->U[k][j][i-1].Fr1 + theta[k][j][i][11] * pMat->U[k][j][i+1].Fr1;

			temp0 = theta[k][j][i][7] * pMat->U[k][j][i].Fr1 + theta[k][j][i][8] * pMat->U[k][j][i].Fr2 + theta[k][j][i][9] * pMat->U[k][j][i].Fr3;

			pMat->U[k][j][i].Er -= ((tempEr1 + tempEr2 + tempEr3));
			pMat->U[k][j][i].Er -= (tempFr1 + tempFr2 + tempFr3) + temp0;

			/* diagonal elements are not included */

			pMat->U[k][j][i].Er /= theta[k][j][i][6];
	
			pMat->U[k][j][i].Er = (1.0 - omega) * Ert0 + omega * pMat->U[k][j][i].Er;

			/*****************************************************/
			/* For Fr1 */

			Frt0 = pMat->U[k][j][i].Fr1;

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

			pMat->U[k][j][i].Fr1 = (1.0 - omega) * Frt0 + omega * pMat->U[k][j][i].Fr1;

			/**************************************************/

			/* For Fr2 */

			Frt0 = pMat->U[k][j][i].Fr2;

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

			pMat->U[k][j][i].Fr2 = (1.0 - omega) * Frt0 + omega * pMat->U[k][j][i].Fr2;
	
			/***************************************************/
			/* For Fr3 */

			Frt0 = pMat->U[k][j][i].Fr3;

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

			pMat->U[k][j][i].Fr3 = (1.0 - omega) * Frt0 + omega * pMat->U[k][j][i].Fr3;

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

