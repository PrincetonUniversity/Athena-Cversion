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

/* In Gauss-Seidel scheme, we do not need a temporary array */
/* This is not exactly the original Gauss-Seidel scheme *
 * For ghost zones, we do not use the updated version */
/* So differenet CPUs do not need wait for each other to finish */

void GaussSeidel1D(MatrixS *pMat)
{

	int i, n;
	int is, ie, js, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	ks = pMat->ks;


	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;
	Real hdtodx2 = 0.5 * pMat->dt/pMat->dx2;

	Real dt = pMat->dt;

	/* To store the coefficient */
	Real theta[6];
	Real phi[6];


	/* Temporary variables to setup the matrix */
	Real velocity_x, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1;

			
	/* Update the boundary cells */


/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){


			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
		/* The right hand side is stored in pMat */
		/* The coefficients are calculated according to the formula */

		
			velocity_x = pMat->U[ks][js][i].V1;

			T4 = pMat->U[ks][js][i].T4;
				
				
			Sigma_sF = pMat->U[ks][js][i].Sigma[0];
			Sigma_aF = pMat->U[ks][js][i].Sigma[1];
			Sigma_aP = pMat->U[ks][js][i].Sigma[2];
			Sigma_aE = pMat->U[ks][js][i].Sigma[3];
	
			Ci0 = (sqrt(pMat->U[ks][js][i].Edd_11) - sqrt(pMat->U[ks][js][i-1].Edd_11)) 
				/ (sqrt(pMat->U[ks][js][i].Edd_11) + sqrt(pMat->U[ks][js][i-1].Edd_11));
			Ci1 =  (sqrt(pMat->U[ks][js][i+1].Edd_11) - sqrt(pMat->U[ks][js][i].Edd_11)) 
				/ (sqrt(pMat->U[ks][js][i+1].Edd_11) + sqrt(pMat->U[ks][js][i].Edd_11));
			
			
			
			theta[0] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][js][i-1].Edd_11);
			theta[1] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[2] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][js][i].Edd_11);
/* 
				+ Crat * pMat->dt * Sigma_aE 
				+ pMat->dt * (Sigma_aF - Sigma_sF) * (1.0 + pMat->U[ks][js][i].Edd_11) * velocity_x * velocity_x / Crat;
*/
			theta[3] = Crat * hdtodx1 * (Ci0 + Ci1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x; */
			theta[4] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[ks][js][i+1].Edd_11);
			theta[5] = Crat * hdtodx1 * (1.0 - Ci1);
			
			
			phi[0] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[ks][js][i-1].Edd_11;
			phi[1] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][js][i-1].Edd_11);
			phi[2] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[ks][js][i].Edd_11 
			       - pMat->dt * (Sigma_aF + Sigma_sF) * (1.0 + pMat->U[ks][js][i].Edd_11) * velocity_x 
			       + pMat->dt * Sigma_aE * velocity_x;
			phi[3] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][js][i].Edd_11) 
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			phi[4] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[ks][js][i+1].Edd_11;
			phi[5] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[ks][js][i+1].Edd_11);
		
			
		

			/* The diagonal elements are theta[4], phi[5], psi[5] */
			
			/* For Er */
			pMat->U[ks][js][i].Er  = pMat->RHS[ks][js][i][0];

			pMat->U[ks][js][i].Er -= theta[0] * pMat->U[ks][js][i-1].Er;
			pMat->U[ks][js][i].Er -= theta[1] * pMat->U[ks][js][i-1].Fr1;
			/* diagonal elements are not included */
			pMat->U[ks][js][i].Er -= theta[3] * pMat->U[ks][js][i].Fr1;
			pMat->U[ks][js][i].Er -= theta[4] * pMat->U[ks][js][i+1].Er;
			pMat->U[ks][js][i].Er -= theta[5] * pMat->U[ks][js][i+1].Fr1;

			pMat->U[ks][js][i].Er /= theta[2];
		
			/* For Fr1 */

			pMat->U[ks][js][i].Fr1  = pMat->RHS[ks][js][i][1];

			pMat->U[ks][js][i].Fr1 -= phi[0] * pMat->U[ks][js][i-1].Er;
			pMat->U[ks][js][i].Fr1 -= phi[1] * pMat->U[ks][js][i-1].Fr1;
			pMat->U[ks][js][i].Fr1 -= phi[2] * pMat->U[ks][js][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][js][i].Fr1 -= phi[4] * pMat->U[ks][js][i+1].Er;
			pMat->U[ks][js][i].Fr1 -= phi[5] * pMat->U[ks][js][i+1].Fr1;

			pMat->U[ks][js][i].Fr1 /= phi[3];
		

	}
			
		
	/* Update the boundary cells */
	bvals_Matrix(pMat);


}	
  

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */
