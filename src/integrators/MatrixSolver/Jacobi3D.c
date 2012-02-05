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


	/* To store the coefficient */
	Real theta[16];
	Real phi[16];
	Real psi[16];
	Real varphi[16];





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

		/* Only need to set the elements once, at the beginning */
		matrix_coef(pMat, NULL, 3, i, j, k, 0.0, &(theta[0]), &(phi[0]), &(psi[0]), &(varphi[0]));
		

		/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
		/* For Er */
		pMatnew->U[k][j][i].Er = pMat->RHS[k][j][i][0];
		pMatnew->U[k][j][i].Er -= theta[0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Er -= theta[1] * pMat->U[k-1][j][i].Fr3;
		pMatnew->U[k][j][i].Er -= theta[2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Er -= theta[3] * pMat->U[k][j-1][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Er -= theta[5] * pMat->U[k][j][i-1].Fr1;
		/* diagonal elements are not included */
		pMatnew->U[k][j][i].Er -= theta[7] * pMat->U[k][j][i].Fr1;
		pMatnew->U[k][j][i].Er -= theta[8] * pMat->U[k][j][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[9] * pMat->U[k][j][i].Fr3;
		pMatnew->U[k][j][i].Er -= theta[10] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Er -= theta[11] * pMat->U[k][j][i+1].Fr1;
		pMatnew->U[k][j][i].Er -= theta[12] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Er -= theta[13] * pMat->U[k][j+1][i].Fr2;
		pMatnew->U[k][j][i].Er -= theta[14] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Er -= theta[15] * pMat->U[k+1][j][i].Fr3;
		pMatnew->U[k][j][i].Er /= theta[6];
			/* For Fr1 */

		pMatnew->U[k][j][i].Fr1 =  pMat->RHS[k][j][i][1];
		pMatnew->U[k][j][i].Fr1 -= phi[0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[1] * pMat->U[k-1][j][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[3] * pMat->U[k][j-1][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[5] * pMat->U[k][j][i-1].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[6] * pMat->U[k][j][i].Er;
		/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr1 -= phi[8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[9] * pMat->U[k][j][i+1].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[11] * pMat->U[k][j+1][i].Fr1;
		pMatnew->U[k][j][i].Fr1 -= phi[12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr1 -= phi[13] * pMat->U[k+1][j][i].Fr1;
	
		pMatnew->U[k][j][i].Fr1 /= phi[7];

			/* For Fr2 */

		pMatnew->U[k][j][i].Fr2 = pMat->RHS[k][j][i][2];
		pMatnew->U[k][j][i].Fr2 -= psi[0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[1] * pMat->U[k-1][j][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[3] * pMat->U[k][j-1][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[5] * pMat->U[k][j][i-1].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[6] * pMat->U[k][j][i].Er;
			/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr2 -= psi[8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[9] * pMat->U[k][j][i+1].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[11] * pMat->U[k][j+1][i].Fr2;
		pMatnew->U[k][j][i].Fr2 -= psi[12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr2 -= psi[13] * pMat->U[k+1][j][i].Fr2;
	
		pMatnew->U[k][j][i].Fr2 /= psi[7];

			/* For Fr3 */

		pMatnew->U[k][j][i].Fr3 = pMat->RHS[k][j][i][3];
		pMatnew->U[k][j][i].Fr3 -= varphi[0] * pMat->U[k-1][j][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[1] * pMat->U[k-1][j][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[2] * pMat->U[k][j-1][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[3] * pMat->U[k][j-1][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[4] * pMat->U[k][j][i-1].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[5] * pMat->U[k][j][i-1].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[6] * pMat->U[k][j][i].Er;
		/* diagonal elements are not included */

		pMatnew->U[k][j][i].Fr3 -= varphi[8] * pMat->U[k][j][i+1].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[9] * pMat->U[k][j][i+1].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[10] * pMat->U[k][j+1][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[11] * pMat->U[k][j+1][i].Fr3;
		pMatnew->U[k][j][i].Fr3 -= varphi[12] * pMat->U[k+1][j][i].Er;
		pMatnew->U[k][j][i].Fr3 -= varphi[13] * pMat->U[k+1][j][i].Fr3;
	
		pMatnew->U[k][j][i].Fr3 /= varphi[7];

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

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX-MULTIGRID */
