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

	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;
	Real hdtodx2 = 0.5 * pMat->dt/pMat->dx2;
	Real hdtodx3 = 0.5 * pMat->dt/pMat->dx3;
	Real dt = pMat->dt;

	/* To store the coefficient */
	Real theta[16];
	Real phi[16];
	Real psi[16];
	Real varphi[16];

	/* Temporary variables to setup the matrix */
	Real velocity_x, velocity_y, velocity_z, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1, Cj0, Cj1, Ck0, Ck1;


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
	
		velocity_x = pMat->U[k][j][i].V1;
		velocity_y = pMat->U[k][j][i].V2;
		velocity_z = pMat->U[k][j][i].V3;
		T4 = pMat->U[k][j][i].T4;
				
		/* Assuming the velocity is already the original velocity in case of FARGO */			
				
			Sigma_sF = pMat->U[k][j][i].Sigma[0];
			Sigma_aF = pMat->U[k][j][i].Sigma[1];
			Sigma_aP = pMat->U[k][j][i].Sigma[2];
			Sigma_aE = pMat->U[k][j][i].Sigma[3];

			Ci0 = (sqrt(pMat->U[k][j][i].Edd_11) - sqrt(pMat->U[k][j][i-1].Edd_11)) 
				/ (sqrt(pMat->U[k][j][i].Edd_11) + sqrt(pMat->U[k][j][i-1].Edd_11));
			Ci1 =  (sqrt(pMat->U[k][j][i+1].Edd_11) - sqrt(pMat->U[k][j][i].Edd_11)) 
				/ (sqrt(pMat->U[k][j][i+1].Edd_11) + sqrt(pMat->U[k][j][i].Edd_11));
			Cj0 = (sqrt(pMat->U[k][j][i].Edd_22) - sqrt(pMat->U[k][j-1][i].Edd_22)) 
				/ (sqrt(pMat->U[k][j][i].Edd_22) + sqrt(pMat->U[k][j-1][i].Edd_22));
			Cj1 =  (sqrt(pMat->U[k][j+1][i].Edd_22) - sqrt(pMat->U[k][j][i].Edd_22)) 
				/ (sqrt(pMat->U[k][j+1][i].Edd_22) + sqrt(pMat->U[k][j][i].Edd_22));
			Ck0 = (sqrt(pMat->U[k][j][i].Edd_33) - sqrt(pMat->U[k-1][j][i].Edd_33)) 
				/ (sqrt(pMat->U[k][j][i].Edd_33) + sqrt(pMat->U[k-1][j][i].Edd_33));
			Ck1 =  (sqrt(pMat->U[k+1][j][i].Edd_33) - sqrt(pMat->U[k][j][i].Edd_33)) 
				/ (sqrt(pMat->U[k+1][j][i].Edd_33) + sqrt(pMat->U[k][j][i].Edd_33));
			theta[0] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			theta[1] = -Crat * hdtodx3 * (1.0 + Ck0);
			theta[2] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			theta[3] = -Crat * hdtodx2 * (1.0 + Cj0);
			theta[4] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			theta[5] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[6] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				+ Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22)
				+ Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)
				+ Eratio * (Crat * pMat->dt * Sigma_aE
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_11) * velocity_x 
				+ velocity_y * pMat->U[k][j][i].Edd_21 + velocity_z * pMat->U[k][j][i].Edd_31) * velocity_x / Crat
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_22) * velocity_y 
				+ velocity_x * pMat->U[k][j][i].Edd_21 + velocity_z * pMat->U[k][j][i].Edd_32) * velocity_y / Crat
				+ pMat->dt * (Sigma_aF - Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_33) * velocity_z 
				+ velocity_x * pMat->U[k][j][i].Edd_31 + velocity_y * pMat->U[k][j][i].Edd_32) * velocity_z / Crat);

			theta[7] = Crat * hdtodx1 * (Ci0 + Ci1)	- Eratio * pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x;
			theta[8] = Crat * hdtodx2 * (Cj0 + Cj1)	- Eratio * pMat->dt * (Sigma_aF - Sigma_sF) * velocity_y;
			theta[9] = Crat * hdtodx3 * (Ck0 + Ck1)	- Eratio * pMat->dt * (Sigma_aF - Sigma_sF) * velocity_z;
			theta[10] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			theta[11] = Crat * hdtodx1 * (1.0 - Ci1);
			theta[12] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			theta[13] = Crat * hdtodx2 * (1.0 - Cj1);
			theta[14] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);
			theta[15] = Crat * hdtodx3 * (1.0 - Ck1);
			
			
			phi[0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_31;
			phi[1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			phi[2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_21;
			phi[3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			phi[4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_11;
			phi[5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			phi[6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_11
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_21   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_31   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_11) * velocity_x + pMat->U[k][j][i].Edd_21 * velocity_y + pMat->U[k][j][i].Edd_31 * velocity_z) 
			       + pMat->dt * Sigma_aE * velocity_x;
			phi[7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			phi[8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_11;
			phi[9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			phi[10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_21;
			phi[11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			phi[12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_31;
			phi[13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);



			psi[0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_32;
			psi[1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			psi[2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_22;
			psi[3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			psi[4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_21;
			psi[5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			psi[6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_21
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_22   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_32   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_22) * velocity_y + pMat->U[k][j][i].Edd_21 * velocity_x + pMat->U[k][j][i].Edd_32 * velocity_z) 
			       + pMat->dt * Sigma_aE * velocity_y;
			psi[7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			psi[8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_21;
			psi[9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			psi[10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_22;
			psi[11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			psi[12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_32;
			psi[13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);

			varphi[0] = -Crat * hdtodx3 * (1.0 + Ck0) * pMat->U[k-1][j][i].Edd_33;
			varphi[1] = -Crat * hdtodx3 * (1.0 + Ck0) * sqrt(pMat->U[k-1][j][i].Edd_33);
			varphi[2] = -Crat * hdtodx2 * (1.0 + Cj0) * pMat->U[k][j-1][i].Edd_32;
			varphi[3] = -Crat * hdtodx2 * (1.0 + Cj0) * sqrt(pMat->U[k][j-1][i].Edd_22);
			varphi[4] = -Crat * hdtodx1 * (1.0 + Ci0) * pMat->U[k][j][i-1].Edd_31;
			varphi[5] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[k][j][i-1].Edd_11);
			varphi[6] = Crat * hdtodx1 * (Ci0 + Ci1) * pMat->U[k][j][i].Edd_31
			       + Crat * hdtodx2 * (Cj0 + Cj1) * pMat->U[k][j][i].Edd_32   
			       + Crat * hdtodx3 * (Ck0 + Ck1) * pMat->U[k][j][i].Edd_33   
			       - pMat->dt * (Sigma_aF + Sigma_sF) * ((1.0 + pMat->U[k][j][i].Edd_33) * velocity_z + pMat->U[k][j][i].Edd_31 * velocity_x + pMat->U[k][j][i].Edd_32 * velocity_y) 
			       + pMat->dt * Sigma_aE * velocity_z;
			varphi[7] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[k][j][i].Edd_11) 
				     + Crat * hdtodx2 * (2.0 + Cj1 - Cj0) * sqrt(pMat->U[k][j][i].Edd_22) 
				     + Crat * hdtodx3 * (2.0 + Ck1 - Ck0) * sqrt(pMat->U[k][j][i].Edd_33)	
				     + Crat * pMat->dt * (Sigma_aF + Sigma_sF);
			varphi[8] = Crat * hdtodx1 * (1.0 - Ci1) * pMat->U[k][j][i+1].Edd_31;
			varphi[9] = -Crat * hdtodx1 * (1.0 - Ci1) * sqrt(pMat->U[k][j][i+1].Edd_11);
			varphi[10] = Crat * hdtodx2 * (1.0 - Cj1) * pMat->U[k][j+1][i].Edd_32;
			varphi[11] = -Crat * hdtodx2 * (1.0 - Cj1) * sqrt(pMat->U[k][j+1][i].Edd_22);
			varphi[12] = Crat * hdtodx3 * (1.0 - Ck1) * pMat->U[k+1][j][i].Edd_33;
			varphi[13] = -Crat * hdtodx3 * (1.0 - Ck1) * sqrt(pMat->U[k+1][j][i].Edd_33);
				

		

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
