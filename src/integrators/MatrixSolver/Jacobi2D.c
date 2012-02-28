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

#if defined(MATRIX_MULTIGRID) || defined(MATRIX_HYPRE) 

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */


/* Matrix boundary condition functions */
extern void bvals_Matrix(MatrixS *pMat);


void Jacobi2D(MatrixS *pMat)
{
/* Right now, only work for one domain. Modified later for SMR */


	

	int i, j, n;
	int is, ie, js, je, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	
	Real omega = 0.5;

	/* To store the coefficient */
	Real ***theta = NULL;
	Real ***phi = NULL;
	Real ***psi = NULL;
	

	if((theta = (Real***)calloc_3d_array(je-js+1+2*Matghost, ie-is+1+2*Matghost,11,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	if((phi = (Real***)calloc_3d_array(je-js+1+2*Matghost, ie-is+1+2*Matghost,11,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");
	
	if((psi = (Real***)calloc_3d_array(je-js+1+2*Matghost, ie-is+1+2*Matghost,11,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	/* Only need to calculate the coefficient once */
	for(j=js; j<=je; j++)
		for(i=is; i<=ie; i++){				
		matrix_coef(pMat, NULL, 2, i, j, ks, 0.0, &(theta[j][i][0]), &(phi[j][i][0]), &(psi[j][i][0]), NULL);
							
	}



	/* First, allocate memory for the temporary matrix */
	MatrixS *pMatnew;

	if((pMatnew = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");

	if((pMatnew->U = (RadMHDS***)calloc_3d_array(1,pMat->Nx[1]+2*Matghost, pMat->Nx[0]+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");




for(n=0; n<Ncycle; n++){

	for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		

		/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
		/* For Er */
		/* For Er */
			pMatnew->U[ks][j][i].Er  = pMat->RHS[ks][j][i][0];
			pMatnew->U[ks][j][i].Er -= theta[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Er -= theta[j][i][1] * pMat->U[ks][j-1][i].Fr2;
			pMatnew->U[ks][j][i].Er -= theta[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Er -= theta[j][i][3] * pMat->U[ks][j][i-1].Fr1;
			/* diagonal elements are not included */
			pMatnew->U[ks][j][i].Er -= theta[j][i][5] * pMat->U[ks][j][i].Fr1;
			pMatnew->U[ks][j][i].Er -= theta[j][i][6] * pMat->U[ks][j][i].Fr2;
			pMatnew->U[ks][j][i].Er -= theta[j][i][7] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Er -= theta[j][i][8] * pMat->U[ks][j][i+1].Fr1;
			pMatnew->U[ks][j][i].Er -= theta[j][i][9] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Er -= theta[j][i][10] * pMat->U[ks][j+1][i].Fr2;
			pMatnew->U[ks][j][i].Er /= theta[j][i][4];
		
			pMatnew->U[ks][j][i].Er = (1.0 - omega) * pMat->U[ks][j][i].Er + omega * pMatnew->U[ks][j][i].Er;

			/* For Fr1 */

			pMatnew->U[ks][j][i].Fr1  = pMat->RHS[ks][j][i][1];
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][1] * pMat->U[ks][j-1][i].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][3] * pMat->U[ks][j][i-1].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][6] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][7] * pMat->U[ks][j][i+1].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][8] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[j][i][9] * pMat->U[ks][j+1][i].Fr1;
			pMatnew->U[ks][j][i].Fr1 /= phi[j][i][5];

			pMatnew->U[ks][j][i].Fr1 = (1.0 - omega) * pMat->U[ks][j][i].Fr1 + omega * pMatnew->U[ks][j][i].Fr1;

			/* For Fr2 */

			pMatnew->U[ks][j][i].Fr2  = pMat->RHS[ks][j][i][2];
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][1] * pMat->U[ks][j-1][i].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][3] * pMat->U[ks][j][i-1].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][6] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][7] * pMat->U[ks][j][i+1].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][8] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[j][i][9] * pMat->U[ks][j+1][i].Fr2;
			pMatnew->U[ks][j][i].Fr2 /= psi[j][i][5];

			pMatnew->U[ks][j][i].Fr2 = (1.0 - omega) * pMat->U[ks][j][i].Fr2 + omega * pMatnew->U[ks][j][i].Fr2;

	}
	
	/* Copy the new values to the old values for the cycle */
	for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

				pMat->U[ks][j][i].Er = pMatnew->U[ks][j][i].Er;
				pMat->U[ks][j][i].Fr1 = pMatnew->U[ks][j][i].Fr1;
				pMat->U[ks][j][i].Fr2 = pMatnew->U[ks][j][i].Fr2;
	} 

		
	/* Update the boundary cells */
	bvals_Matrix(pMat);

}	
  

	/* Free the temporary matrix */
	free_3d_array(pMatnew->U);
	free(pMatnew);

	
	if(theta != NULL)
		free_3d_array(theta);

	if(phi != NULL)
		free_3d_array(phi);

	if(psi != NULL)
		free_3d_array(psi);



	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX-MULTIGRID */
