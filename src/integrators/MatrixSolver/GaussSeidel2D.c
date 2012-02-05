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

void GaussSeidel2D(MatrixS *pMat)
{

	int i, j, n;
	int is, ie, js, je, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;



	/* To store the coefficient */
	Real theta[11];
	Real phi[11];
	Real psi[11];


			
	/* Update the boundary cells */


/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){

		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
		/* The right hand side is stored in pMat */
		/* The coefficients are calculated according to the formula */

			matrix_coef(pMat, NULL, 2, i, j, ks, 0.0, &(theta[0]), &(phi[0]), &(psi[0]), NULL);
		

			/* The diagonal elements are theta[4], phi[5], psi[5] */
			
			/* For Er */
			pMat->U[ks][j][i].Er  = pMat->RHS[ks][j][i][0];
			pMat->U[ks][j][i].Er -= theta[0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Er -= theta[1] * pMat->U[ks][j-1][i].Fr2;
			pMat->U[ks][j][i].Er -= theta[2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Er -= theta[3] * pMat->U[ks][j][i-1].Fr1;
			/* diagonal elements are not included */
			pMat->U[ks][j][i].Er -= theta[5] * pMat->U[ks][j][i].Fr1;
			pMat->U[ks][j][i].Er -= theta[6] * pMat->U[ks][j][i].Fr2;
			pMat->U[ks][j][i].Er -= theta[7] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Er -= theta[8] * pMat->U[ks][j][i+1].Fr1;
			pMat->U[ks][j][i].Er -= theta[9] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Er -= theta[10] * pMat->U[ks][j+1][i].Fr2;
			pMat->U[ks][j][i].Er /= theta[4];
		
			/* For Fr1 */

			pMat->U[ks][j][i].Fr1  = pMat->RHS[ks][j][i][1];
			pMat->U[ks][j][i].Fr1 -= phi[0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Fr1 -= phi[1] * pMat->U[ks][j-1][i].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Fr1 -= phi[3] * pMat->U[ks][j][i-1].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][j][i].Fr1 -= phi[6] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Fr1 -= phi[7] * pMat->U[ks][j][i+1].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[8] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Fr1 -= phi[9] * pMat->U[ks][j+1][i].Fr1;
			pMat->U[ks][j][i].Fr1 /= phi[5];

			/* For Fr2 */

			pMat->U[ks][j][i].Fr2  = pMat->RHS[ks][j][i][2];
			pMat->U[ks][j][i].Fr2 -= psi[0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Fr2 -= psi[1] * pMat->U[ks][j-1][i].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Fr2 -= psi[3] * pMat->U[ks][j][i-1].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][j][i].Fr2 -= psi[6] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Fr2 -= psi[7] * pMat->U[ks][j][i+1].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[8] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Fr2 -= psi[9] * pMat->U[ks][j+1][i].Fr2;
			pMat->U[ks][j][i].Fr2 /= psi[5];

	}
			
		
	/* Update the boundary cells */
	bvals_Matrix(pMat);


}	
  

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */
