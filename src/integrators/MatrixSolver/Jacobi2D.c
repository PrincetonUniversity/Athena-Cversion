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
	

	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;
	Real hdtodx2 = 0.5 * pMat->dt/pMat->dx2;

	Real dt = pMat->dt;

	/* To store the coefficient */
	Real theta[11];
	Real phi[11];
	Real psi[11];



	/* Temporary variables to setup the matrix */
	Real velocity_x, velocity_y, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1, Cj0, Cj1;


	/* First, allocate memory for the temporary matrix */
	MatrixS *pMatnew;

	if((pMatnew = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");

	if((pMatnew->U = (RadMHDS***)calloc_3d_array(1,pMat->Nx[1]+2*Matghost, pMat->Nx[0]+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[Jacobi3D]: malloc return a NULL pointer\n");




for(n=0; n<Ncycle; n++){

	for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
	
			velocity_x = pMat->U[ks][j][i].V1;
			velocity_y = pMat->U[ks][j][i].V2;

			T4 = pMat->U[ks][j][i].T4;
				
			/* Assuming the velocity is already the original velocity in case of FARGO */			
				
			Sigma_sF = pMat->U[ks][j][i].Sigma[0];
			Sigma_aF = pMat->U[ks][j][i].Sigma[1];
			Sigma_aP = pMat->U[ks][j][i].Sigma[2];
			Sigma_aE = pMat->U[ks][j][i].Sigma[3];

			Ci0 = (sqrt(pMat->U[ks][j][i].Edd_11) - sqrt(pMat->U[ks][j][i-1].Edd_11)) 
				/ (sqrt(pMat->U[ks][j][i].Edd_11) + sqrt(pMat->U[ks][j][i-1].Edd_11));
			Ci1 =  (sqrt(pMat->U[ks][j][i+1].Edd_11) - sqrt(pMat->U[ks][j][i].Edd_11)) 
				/ (sqrt(pMat->U[ks][j][i+1].Edd_11) + sqrt(pMat->U[ks][j][i].Edd_11));
			Cj0 = (sqrt(pMat->U[ks][j][i].Edd_22) - sqrt(pMat->U[ks][j-1][i].Edd_22)) 
				/ (sqrt(pMat->U[ks][j][i].Edd_22) + sqrt(pMat->U[ks][j-1][i].Edd_22));
			Cj1 =  (sqrt(pMat->U[ks][j+1][i].Edd_22) - sqrt(pMat->U[ks][j][i].Edd_22)) 
				/ (sqrt(pMat->U[ks][j+1][i].Edd_22) + sqrt(pMat->U[ks][j][i].Edd_22));
			
			theta[0] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[ks][j-1][i].Edd_22);
			theta[1] = -Crat * hdtodx2 * (1.0 + Cj0);
			theta[2] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][j][i-1].Edd_11);
			theta[3] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[4] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][j][i].Edd_11) 
				+ Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[ks][j][i].Edd_22);
/*
				+ Crat * pMat->dt * Sigma_aE 
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[ks][j][i].Edd_11) * velocity_x 
				+ velocity_y * pMat->U[ks][j][i].Edd_21) * velocity_x / Crat
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[ks][j][i].Edd_22) * velocity_y 
				+ velocity_x * pMat->U[ks][j][i].Edd_21) * velocity_y / Crat;
*/
			theta[5] = Crat * hdtodx1 * (Ci0 + Ci1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x; */
			theta[6] = Crat * hdtodx2 * (Cj0 + Cj1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_y; */
			theta[7] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[ks][j][i+1].Edd_11);
			theta[8] = Crat * hdtodx1 * (1.0 - Ci1);
			theta[9] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[ks][j+1][i].Edd_22);
			theta[10] = Crat * hdtodx2 * (1.0 - Cj1);

			
			phi[0] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[ks][j-1][i].Edd_21;
			phi[1] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[ks][j-1][i].Edd_22);
			phi[2] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[ks][j][i-1].Edd_11;
			phi[3] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][j][i-1].Edd_11);
			phi[4] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[ks][j][i].Edd_11
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[ks][j][i].Edd_21   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[ks][j][i].Edd_11) * velocity_x + pMat->U[ks][j][i].Edd_21 * velocity_y) 
			       + pMat->dt * Sigma_aE * velocity_x;
			phi[5] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[ks][j][i].Edd_22) 
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			phi[6] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[ks][j][i+1].Edd_11;
			phi[7] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[ks][j][i+1].Edd_11);
			phi[8] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[ks][j+1][i].Edd_21;
			phi[9] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[ks][j+1][i].Edd_22);
			


			psi[0] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[ks][j-1][i].Edd_22;
			psi[1] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[ks][j-1][i].Edd_22);
			psi[2] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[ks][j][i-1].Edd_21;
			psi[3] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][j][i-1].Edd_11);
			psi[4] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[ks][j][i].Edd_21
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[ks][j][i].Edd_22   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[ks][j][i].Edd_22) * velocity_y + pMat->U[ks][j][i].Edd_21 * velocity_x) 
			       + pMat->dt * Sigma_aE * velocity_y;
			psi[5] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[ks][j][i].Edd_22) 
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			psi[6] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[ks][j][i+1].Edd_21;
			psi[7] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[ks][j][i+1].Edd_11);
			psi[8] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[ks][j+1][i].Edd_22;
			psi[9] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[ks][j+1][i].Edd_22);
				

		

		/* The diagonal elements are theta[6], phi[7], psi[7], varphi[7] */
		
		/* For Er */
		/* For Er */
			pMatnew->U[ks][j][i].Er  = pMat->RHS[ks][j][i][0];
			pMatnew->U[ks][j][i].Er -= theta[0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Er -= theta[1] * pMat->U[ks][j-1][i].Fr2;
			pMatnew->U[ks][j][i].Er -= theta[2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Er -= theta[3] * pMat->U[ks][j][i-1].Fr1;
			/* diagonal elements are not included */
			pMatnew->U[ks][j][i].Er -= theta[5] * pMat->U[ks][j][i].Fr1;
			pMatnew->U[ks][j][i].Er -= theta[6] * pMat->U[ks][j][i].Fr2;
			pMatnew->U[ks][j][i].Er -= theta[7] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Er -= theta[8] * pMat->U[ks][j][i+1].Fr1;
			pMatnew->U[ks][j][i].Er -= theta[9] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Er -= theta[10] * pMat->U[ks][j+1][i].Fr2;
			pMatnew->U[ks][j][i].Er /= theta[4];
		
			/* For Fr1 */

			pMatnew->U[ks][j][i].Fr1  = pMat->RHS[ks][j][i][1];
			pMatnew->U[ks][j][i].Fr1 -= phi[0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[1] * pMat->U[ks][j-1][i].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[3] * pMat->U[ks][j][i-1].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMatnew->U[ks][j][i].Fr1 -= phi[6] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[7] * pMat->U[ks][j][i+1].Fr1;
			pMatnew->U[ks][j][i].Fr1 -= phi[8] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Fr1 -= phi[9] * pMat->U[ks][j+1][i].Fr1;
			pMatnew->U[ks][j][i].Fr1 /= phi[5];

			/* For Fr2 */

			pMatnew->U[ks][j][i].Fr2  = pMat->RHS[ks][j][i][2];
			pMatnew->U[ks][j][i].Fr2 -= psi[0] * pMat->U[ks][j-1][i].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[1] * pMat->U[ks][j-1][i].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[2] * pMat->U[ks][j][i-1].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[3] * pMat->U[ks][j][i-1].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMatnew->U[ks][j][i].Fr2 -= psi[6] * pMat->U[ks][j][i+1].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[7] * pMat->U[ks][j][i+1].Fr2;
			pMatnew->U[ks][j][i].Fr2 -= psi[8] * pMat->U[ks][j+1][i].Er;
			pMatnew->U[ks][j][i].Fr2 -= psi[9] * pMat->U[ks][j+1][i].Fr2;
			pMatnew->U[ks][j][i].Fr2 /= psi[5];

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

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX-MULTIGRID */
