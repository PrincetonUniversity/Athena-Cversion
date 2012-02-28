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


void Jacobi3D(MatrixS *pMat)
{
/* Right now, only work for one domain. Modified later for SMR */


	

	int i, j, k, n;
	int is, ie, js, je, ks, ke;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	Real omega = 0.5;

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


	/* Only need to calculate the coefficient once */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){				
			matrix_coef(pMat, NULL, 3, i, j, k, 0.0, &(theta[k][j][i][0]), &(phi[k][j][i][0]), &(psi[k][j][i][0]), &(varphi[k][j][i][0]));
							
	}


	/* First, allocate memory for the temporary matrix */
	MatrixS *pMatnew;

	if((pMatnew = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");

	if((pMatnew->U = (RadMHDS***)calloc_3d_array(pMat->Nx[2]+2*Matghost,pMat->Nx[1]+2*Matghost, pMat->Nx[0]+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");




for(n=0; n<Ncycle; n++){

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){


		/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
		/* For Er */
		pMatnew->U[k][j][i].Er = pMat->RHS[k][j][i][0];
		pMatnew->U[k][j][i].Er -= theta[k][j][i][0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][1] * pMat->U[k-1][j][i].Fr3;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][3] * pMat->U[k][j-1][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][5] * pMat->U[k][j][i-1].Fr1;
		/* diagonal elements are not included */
		pMatnew->U[k][j][i].Er -= theta[k][j][i][7] * pMat->U[k][j][i].Fr1;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][8] * pMat->U[k][j][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][9] * pMat->U[k][j][i].Fr3;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][10] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][11] * pMat->U[k][j][i+1].Fr1;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][12] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][13] * pMat->U[k][j+1][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][14] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Er -= theta[k][j][i][15] * pMat->U[k+1][j][i].Fr3;
		pMatnew->U[k][j][i].Er /= theta[k][j][i][6];

		pMatnew->U[k][j][i].Er = (1.0 - omega) * pMat->U[k][j][i].Er + omega * pMatnew->U[k][j][i].Er;	

			/* For Fr1 */

		pMatnew->U[k][j][i].Fr1 =  pMat->RHS[k][j][i][1];
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][1] * pMat->U[k-1][j][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][3] * pMat->U[k][j-1][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][5] * pMat->U[k][j][i-1].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][6] * pMat->U[k][j][i].Er;
		/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][9] * pMat->U[k][j][i+1].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][11] * pMat->U[k][j+1][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[k][j][i][13] * pMat->U[k+1][j][i].Fr1;
	
		pMatnew->U[k][j][i].Fr1 /= phi[k][j][i][7];

		pMatnew->U[k][j][i].Fr1 = (1.0 - omega) * pMat->U[k][j][i].Fr1 + omega * pMatnew->U[k][j][i].Fr1;	

			/* For Fr2 */

		pMatnew->U[k][j][i].Fr2 = pMat->RHS[k][j][i][2];
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][1] * pMat->U[k-1][j][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][3] * pMat->U[k][j-1][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][5] * pMat->U[k][j][i-1].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][6] * pMat->U[k][j][i].Er;
			/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][9] * pMat->U[k][j][i+1].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][11] * pMat->U[k][j+1][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[k][j][i][13] * pMat->U[k+1][j][i].Fr2;
	
		pMatnew->U[k][j][i].Fr2 /= psi[k][j][i][7];

		pMatnew->U[k][j][i].Fr2 = (1.0 - omega) * pMat->U[k][j][i].Fr2 + omega * pMatnew->U[k][j][i].Fr2;	

			/* For Fr3 */

		pMatnew->U[k][j][i].Fr3 = pMat->RHS[k][j][i][3];
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][1] * pMat->U[k-1][j][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][3] * pMat->U[k][j-1][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][5] * pMat->U[k][j][i-1].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][6] * pMat->U[k][j][i].Er;
		/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][9] * pMat->U[k][j][i+1].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][11] * pMat->U[k][j+1][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[k][j][i][13] * pMat->U[k+1][j][i].Fr3;
	
		pMatnew->U[k][j][i].Fr3 /= varphi[k][j][i][7];

		pMatnew->U[k][j][i].Fr3 = (1.0 - omega) * pMat->U[k][j][i].Fr3 + omega * pMatnew->U[k][j][i].Fr3;	

	}
	
	/* Copy the new values to the old values for the cycle */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

				pMat->U[k][j][i].Er = pMatnew->U[k][j][i].Er;
				pMat->U[k][j][i].Fr1 = pMatnew->U[k][j][i].Fr1;
				pMat->U[k][j][i].Fr2 = pMatnew->U[k][j][i].Fr2;
				pMat->U[k][j][i].Fr3 = pMatnew->U[k][j][i].Fr3;
	} 

		
	/* Update the boundary cells */
	bvals_Matrix(pMat);

}	
  

	/* Free the temporary matrix */
	free_3d_array(pMatnew->U);
	free(pMatnew);

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

#endif /* MATRIX-MULTIGRID */
