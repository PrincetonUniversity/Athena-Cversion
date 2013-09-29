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

void GaussSeidel2D(MatrixS *pMat, Real ***theta,  Real ***phi,  Real ***psi)
{

	int i, j, n;
	int is, ie, js, je, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;

	Real omega = 0.4;
#ifdef FLD
	omega = 1.0;
#endif
	Real Ert0, Frt0;


			
	/* Update the boundary cells */


/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){

		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
		/* The right hand side is stored in pMat */
		/* The coefficients are calculated according to the formula */
			/* The diagonal elements are theta[4], phi[5], psi[5] */
				
#ifdef FLD
				
			/* For Er */
			Ert0 = pMat->U[ks][j][i].Er;
			pMat->U[ks][j][i].Er  = pMat->RHS[ks][j][i][0];
			pMat->U[ks][j][i].Er -= theta[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Er -= theta[j][i][1] * pMat->U[ks][j][i-1].Er;
				
			pMat->U[ks][j][i].Er -= theta[j][i][3] * pMat->U[ks][j][i+1].Er;			
			pMat->U[ks][j][i].Er -= theta[j][i][4] * pMat->U[ks][j+1][i].Er;
				
			pMat->U[ks][j][i].Er /= theta[j][i][2];
				
			pMat->U[ks][j][i].Er = (1.0 - omega) * Ert0 + omega * pMat->U[ks][j][i].Er;
				
#else
			
			/* For Er */
			Ert0 = pMat->U[ks][j][i].Er;
			pMat->U[ks][j][i].Er  = pMat->RHS[ks][j][i][0];
			pMat->U[ks][j][i].Er -= theta[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Er -= theta[j][i][1] * pMat->U[ks][j-1][i].Fr2;
			pMat->U[ks][j][i].Er -= theta[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Er -= theta[j][i][3] * pMat->U[ks][j][i-1].Fr1;
			/* diagonal elements are not included */
			pMat->U[ks][j][i].Er -= theta[j][i][5] * pMat->U[ks][j][i].Fr1;
			pMat->U[ks][j][i].Er -= theta[j][i][6] * pMat->U[ks][j][i].Fr2;
			pMat->U[ks][j][i].Er -= theta[j][i][7] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Er -= theta[j][i][8] * pMat->U[ks][j][i+1].Fr1;
			pMat->U[ks][j][i].Er -= theta[j][i][9] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Er -= theta[j][i][10] * pMat->U[ks][j+1][i].Fr2;
			pMat->U[ks][j][i].Er /= theta[j][i][4];

			pMat->U[ks][j][i].Er = (1.0 - omega) * Ert0 + omega * pMat->U[ks][j][i].Er;
		
			/* For Fr1 */
			Frt0 = pMat->U[ks][j][i].Fr1;

			pMat->U[ks][j][i].Fr1  = pMat->RHS[ks][j][i][1];
			pMat->U[ks][j][i].Fr1 -= phi[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][1] * pMat->U[ks][j-1][i].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][3] * pMat->U[ks][j][i-1].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][j][i].Fr1 -= phi[j][i][6] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][7] * pMat->U[ks][j][i+1].Fr1;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][8] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Fr1 -= phi[j][i][9] * pMat->U[ks][j+1][i].Fr1;
			pMat->U[ks][j][i].Fr1 /= phi[j][i][5];

			pMat->U[ks][j][i].Fr1 = (1.0 - omega) * Frt0 + omega * pMat->U[ks][j][i].Fr1;

			/* For Fr2 */
			Frt0 = pMat->U[ks][j][i].Fr2;

			pMat->U[ks][j][i].Fr2  = pMat->RHS[ks][j][i][2];
			pMat->U[ks][j][i].Fr2 -= psi[j][i][0] * pMat->U[ks][j-1][i].Er;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][1] * pMat->U[ks][j-1][i].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][2] * pMat->U[ks][j][i-1].Er;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][3] * pMat->U[ks][j][i-1].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][4] * pMat->U[ks][j][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][j][i].Fr2 -= psi[j][i][6] * pMat->U[ks][j][i+1].Er;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][7] * pMat->U[ks][j][i+1].Fr2;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][8] * pMat->U[ks][j+1][i].Er;
			pMat->U[ks][j][i].Fr2 -= psi[j][i][9] * pMat->U[ks][j+1][i].Fr2;
			pMat->U[ks][j][i].Fr2 /= psi[j][i][5];

			pMat->U[ks][j][i].Fr2 = (1.0 - omega) * Frt0 + omega * pMat->U[ks][j][i].Fr2;
				
#endif /* FLD */

	}
			
		
	/* Update the boundary cells */
	bvals_Matrix(pMat);


}	
  	


	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */
