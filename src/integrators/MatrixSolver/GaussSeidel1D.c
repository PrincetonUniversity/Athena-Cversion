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

	Real omega= 0.5;
	
	Real tempEr, tempFr1;

	
	/* To store the coefficient */
	Real **theta = NULL;
	Real **phi = NULL;
	

	if((theta = (Real**)calloc_2d_array(ie-is+1+2*Matghost,6,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	if((phi   = (Real**)calloc_2d_array(ie-is+1+2*Matghost,6,sizeof(Real))) == NULL)
		ath_error("[GaussSeidel3D]: malloc return a NULL pointer\n");

	
	for(i=is; i<=ie; i++){				
		matrix_coef(pMat, NULL, 1, i, js, ks, 0.0, &(theta[i][0]), &(phi[i][0]), NULL, NULL);
							
	}



/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){


			for(i=is; i<=ie; i++){

		/* Only need to set the elements once, at the beginning */
		/* The right hand side is stored in pMat */
		/* The coefficients are calculated according to the formula */

			/* The diagonal elements are theta[4], phi[5], psi[5] */
			
			/* For Er */
			tempEr = pMat->U[ks][js][i].Er;
			pMat->U[ks][js][i].Er  = pMat->RHS[ks][js][i][0];

			pMat->U[ks][js][i].Er -= theta[i][0] * pMat->U[ks][js][i-1].Er;
			pMat->U[ks][js][i].Er -= theta[i][1] * pMat->U[ks][js][i-1].Fr1;
			/* diagonal elements are not included */
			pMat->U[ks][js][i].Er -= theta[i][3] * pMat->U[ks][js][i].Fr1;
			pMat->U[ks][js][i].Er -= theta[i][4] * pMat->U[ks][js][i+1].Er;
			pMat->U[ks][js][i].Er -= theta[i][5] * pMat->U[ks][js][i+1].Fr1;

			pMat->U[ks][js][i].Er /= theta[i][2];

			pMat->U[ks][js][i].Er = (1.0 - omega) * tempEr + omega * pMat->U[ks][js][i].Er;
		
			/* For Fr1 */

			tempFr1 = pMat->U[ks][js][i].Fr1;
			pMat->U[ks][js][i].Fr1  = pMat->RHS[ks][js][i][1];

			pMat->U[ks][js][i].Fr1 -= phi[i][0] * pMat->U[ks][js][i-1].Er;
			pMat->U[ks][js][i].Fr1 -= phi[i][1] * pMat->U[ks][js][i-1].Fr1;
			pMat->U[ks][js][i].Fr1 -= phi[i][2] * pMat->U[ks][js][i].Er;
			/* diagonal elements are not included */

			pMat->U[ks][js][i].Fr1 -= phi[i][4] * pMat->U[ks][js][i+1].Er;
			pMat->U[ks][js][i].Fr1 -= phi[i][5] * pMat->U[ks][js][i+1].Fr1;

			pMat->U[ks][js][i].Fr1 /= phi[i][3];

			pMat->U[ks][js][i].Fr1 = (1.0 - omega) * tempFr1 + omega * pMat->U[ks][js][i].Fr1;
			
		

	}
			
		
	/* Update the boundary cells */
	bvals_Matrix(pMat);


}	

	if(theta != NULL)
		free_2d_array(theta);

	if(phi != NULL)
		free_2d_array(phi);
  

	return;	
	
}


#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */
