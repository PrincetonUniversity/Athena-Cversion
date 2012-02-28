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


void Jacobi1D(MatrixS *pMat)
{
/* Right now, only work for one domain. Modified later for SMR */


	

	int i, n;
	int is, ie, js, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	ks = pMat->ks;
	
	Real omega = 0.5;

	/* To store the coefficient */
	Real theta[6];
	Real phi[6];
	





	/* First, allocate memory for the temporary matrix */
	MatrixS *pMatnew;

	if((pMatnew = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");

	if((pMatnew->U = (RadMHDS***)calloc_3d_array(1,1, pMat->Nx[0]+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");




for(n=0; n<Ncycle; n++){

		for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
	
			matrix_coef(pMat, NULL, 1, i, js, ks, 0.0, &(theta[0]), &(phi[0]), NULL, NULL);
		

		/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
		/* For Er */
		/* For Er */
			pMatnew->U[ks][js][i].Er  = pMat->RHS[ks][js][i][0];			
			pMatnew->U[ks][js][i].Er -= theta[0] * pMat->U[ks][js][i-1].Er;
			pMatnew->U[ks][js][i].Er -= theta[1] * pMat->U[ks][js][i-1].Fr1;
			/* diagonal elements are not included */
			pMatnew->U[ks][js][i].Er -= theta[3] * pMat->U[ks][js][i].Fr1;
			pMatnew->U[ks][js][i].Er -= theta[4] * pMat->U[ks][js][i+1].Er;
			pMatnew->U[ks][js][i].Er -= theta[5] * pMat->U[ks][js][i+1].Fr1;
			
			pMatnew->U[ks][js][i].Er /= theta[2];

			pMatnew->U[ks][js][i].Er = (1.0 - omega) * pMat->U[ks][js][i].Er + omega * pMatnew->U[ks][js][i].Er;
		
			/* For Fr1 */

			pMatnew->U[ks][js][i].Fr1  = pMat->RHS[ks][js][i][1];			
			pMatnew->U[ks][js][i].Fr1 -= phi[0] * pMat->U[ks][js][i-1].Er;
			pMatnew->U[ks][js][i].Fr1 -= phi[1] * pMat->U[ks][js][i-1].Fr1;
			pMatnew->U[ks][js][i].Fr1 -= phi[2] * pMat->U[ks][js][i].Er;
			/* diagonal elements are not included */

			pMatnew->U[ks][js][i].Fr1 -= phi[4] * pMat->U[ks][js][i+1].Er;
			pMatnew->U[ks][js][i].Fr1 -= phi[5] * pMat->U[ks][js][i+1].Fr1;
			
			pMatnew->U[ks][js][i].Fr1 /= phi[3];

			pMatnew->U[ks][js][i].Fr1 = (1.0 - omega) * pMat->U[ks][js][i].Fr1 + omega * pMatnew->U[ks][js][i].Fr1;
			

	}
	
	/* Copy the new values to the old values for the cycle */
		for(i=is; i<=ie; i++){

			pMat->U[ks][js][i].Er = pMatnew->U[ks][js][i].Er;
			pMat->U[ks][js][i].Fr1 = pMatnew->U[ks][js][i].Fr1;
			pMat->U[ks][js][i].Fr2 = pMatnew->U[ks][js][i].Fr2;
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
