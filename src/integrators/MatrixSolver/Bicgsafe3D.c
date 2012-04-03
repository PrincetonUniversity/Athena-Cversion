#include "../../copyright.h"
/*==============================================================================
 * FILE: Bicgsafe3D.c
 *
 * PURPOSE: Use safe biconguate gradient safe convergece. 
 * From Fujino, Fujiwra, Yoshida
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

extern void bvals_matrix_vector(Real ****Apn, MatrixS *pMat);


void initial_step(Real *INInorm, MatrixS *pMat, Real ****theta, Real ****phi, Real ****psi, Real ****varphi, Real *r0rndot, Real *betan, Real ****Pn, Real ****Apn, Real ****Rnstar, Real ****un, Real ****Aun, Real ****yn, Real ****zn);

void iterations(Real *currentnorm, MatrixS *pMat, Real ****theta, Real ****phi, Real ****psi, Real ****varphi, Real *r0dotrn, Real *betan, Real ****Pn, Real ****Apn, Real ****Rnstar, Real ****un, Real ****Aun, Real ****yn, Real ****zn, Real ****Arn);

/* In Gauss-Seidel scheme, we do not need a temporary array */
/* This is not exactly the original Gauss-Seidel scheme *
 * For ghost zones, we do not use the updated version */
/* So differenet CPUs do not need wait for each other to finish */

void Bicgsafe3D(MatrixS *pMat, Real ****theta,  Real ****phi,  Real ****psi,  Real ****varphi)
{


	
	int n;
	int is, ie, js, je, ks, ke;
	
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;



	
	/* To store the temporary vectors used in Bicgsafe iteration */
	Real ****Pn = NULL;
	Real ****Apn = NULL;
	Real ****Rnstar = NULL;
	Real ****yn = NULL;
	Real ****un = NULL;
	Real ****Aun = NULL;	
	Real ****zn = NULL;
	Real ****Arn = NULL;
	
	Real r0dotrn, betan;
	Real INInorm, currentnorm;
	/* Used to decide whether we should stop iteration at this level */



	/* temporary vectors */
	if((Pn = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");

	if((Apn = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");


	if((Rnstar = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");
	
	if((yn = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");

	if((un = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");
	
	if((Aun = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");


	if((zn = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");

	if((Arn = (Real****)calloc_4d_array(ke-ks+1+2*Matghost,je-js+1+2*Matghost, ie-is+1+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[Bicgsafe3D]: malloc return a NULL pointer\n");




			
	/* First, Update the boundary cells */
	/* We do not set ghost zones after prolongation */
	
	/* We only need Er and Fr in the ghost zones */
	/* Use seperate loop to avoid if in next loop */

	
	
	

/* Hardware to Ncycle */
for(n=0; n<Ncycle; n++){
	

	if(n == 0){
		/* Initial betan is set to be zero */
		betan = 0.0;
		initial_step(&INInorm, pMat, theta, phi, psi, varphi, &r0dotrn, &betan, Pn, Apn, Rnstar, un, Aun, yn, zn);	
	/* Update the boundary cells for RHS*/

	}
	else{
		iterations(&currentnorm, pMat, theta, phi, psi, varphi, &r0dotrn, &betan, Pn, Apn, Rnstar, un, Aun, yn, zn, Arn);

		/* If it is converged at this level, do not need to do more */
		if(currentnorm/INInorm < TOL)
			n = Ncycle;	
	}


	/* update ghost zones for the residual */
	bvals_matrix_vector(pMat->RHS, pMat);
	

}

	/* Only need to update the boundary cells */
	/* bicgsafe doesn't need ghost zones of Er, Fr during the iterations */
	bvals_Matrix(pMat);
	
			
  
	if(Pn != NULL)
		free_4d_array(Pn);

	if(Apn != NULL)
		free_4d_array(Apn);

	if(Rnstar != NULL)
		free_4d_array(Rnstar);

	if(yn != NULL)
		free_4d_array(yn);

	if(un != NULL)
		free_4d_array(un);

	if(Aun != NULL)
		free_4d_array(Aun);

	if(zn != NULL)
		free_4d_array(zn);

	if(Arn != NULL)
		free_4d_array(Arn);

	return;	
	
}




void initial_step(Real *INInorm, MatrixS *pMat, Real ****theta, Real ****phi, Real ****psi, Real ****varphi, Real *r0dotrn, Real *betan, Real ****Pn, Real ****Apn, Real ****Rnstar, Real ****un, Real ****Aun, Real ****yn, Real ****zn)
{
	int i, j, k, n;
	int is, ie, js, je, ks, ke;
	
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	Real alphan, qsin;
	Real rndotrn, rndotApn, ApndotApn, rstardotApn, normtemp;
	

#ifdef MPI_PARALLEL
	Real norm[4];
 	Real totnorm[4];
#endif

	/* update ghost zones for the residual */
	bvals_matrix_vector(pMat->RHS, pMat);


	/* First, store RHS to Rnstar */
	/* first, initialize Pn and related vectors, assuming ghost zones of RHS have been updated */
	for(k=ks-Matghost; k<=ke+Matghost; k++)
		for(j=js-Matghost; j<=je+Matghost; j++)
			for(i=is-Matghost; i<=ie+Matghost; i++){
				for(n=0; n<4; n++){					
					Rnstar[k][j][i][n] = pMat->RHS[k][j][i][n];	
					Pn[k][j][i][n] = pMat->RHS[k][j][i][n];		
				}
	}

	/* now calculate Apn and related quantities */
	
		
	rndotrn     = 0.0;
	rndotApn    = 0.0;
	rstardotApn = 0.0;
	ApndotApn   = 0.0;
		
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		matrix_vector_product3D(theta[k][j][i], phi[k][j][i], psi[k][j][i], varphi[k][j][i], i, j, k, Pn, Apn[k][j][i]);
		
		vector_product(Rnstar[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
		rndotrn += normtemp;
		vector_product(Rnstar[k][j][i],Apn[k][j][i], 4 , &normtemp);
		rstardotApn += normtemp;
		vector_product(pMat->RHS[k][j][i],Apn[k][j][i], 4 , &normtemp);
		rndotApn += normtemp;
		vector_product(Apn[k][j][i],Apn[k][j][i], 4 , &normtemp);
		ApndotApn += normtemp;	
	}

#ifdef MPI_PARALLEL
	norm[0] = rndotrn;
	norm[1] = rndotApn;
	norm[2] = rstardotApn;
	norm[3] = ApndotApn;

	MPI_Allreduce(&norm[0],&totnorm[0],4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	rndotrn     = totnorm[0];
	rndotApn    = totnorm[1];
	rstardotApn = totnorm[2];
	ApndotApn   = totnorm[3];
#endif

	/* store the initial norm */
	*INInorm = rndotrn;

	/* stop the calculation if norm is close to zero */
	if(fabs(rndotApn) < TINY_NUMBER || fabs(ApndotApn) < TINY_NUMBER)
		return;

	alphan = rndotrn / rstardotApn;
	qsin   = rndotApn / ApndotApn;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
			for(n=0; n<4; n++){
				un[k][j][i][n] = qsin * Apn[k][j][i][n];
				zn[k][j][i][n] = qsin * Rnstar[k][j][i][n] - alphan * un[k][j][i][n];
			}
	}

	/* update boundary condition for un in order to calculate Aun */
	bvals_matrix_vector(un, pMat);
	
	/* now calculate Aun */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		matrix_vector_product3D(theta[k][j][i], phi[k][j][i], psi[k][j][i], varphi[k][j][i], i, j, k, un, Aun[k][j][i]);		
		
	}


	/* now calculate yn */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				for(n=0; n<4; n++){
					yn[k][j][i][n] = qsin * Apn[k][j][i][n] - alphan * Aun[k][j][i][n];
					pMat->RHS[k][j][i][n] -= alphan * Apn[k][j][i][n];
					pMat->RHS[k][j][i][n] -= yn[k][j][i][n];			
				}
			/* update the guess solution */
			pMat->U[k][j][i].Er  += alphan * Pn[k][j][i][0] + zn[k][j][i][0];
			pMat->U[k][j][i].Fr1 += alphan * Pn[k][j][i][1] + zn[k][j][i][1];
			pMat->U[k][j][i].Fr2 += alphan * Pn[k][j][i][2] + zn[k][j][i][2];
			pMat->U[k][j][i].Fr3 += alphan * Pn[k][j][i][3] + zn[k][j][i][3];
	}
	
	/* now calculate betan+1 */
	(*r0dotrn) = 0.0;
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		vector_product(Rnstar[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
		(*r0dotrn) += normtemp;			
	}
	
#ifdef MPI_PARALLEL
	norm[0] = (*r0dotrn);
	
	MPI_Allreduce(&norm[0],&totnorm[0],1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	

	(*r0dotrn) = totnorm[0];
#endif

	/* norm[0] is <r0, rn> */
	if(fabs(qsin) < TINY_NUMBER)
		*betan = 0.0;
	else
		*betan = (alphan/qsin) * (*r0dotrn) / rndotrn;


	return;

}

void iterations(Real *currentnorm, MatrixS *pMat, Real ****theta, Real ****phi, Real ****psi, Real ****varphi, Real *r0dotrn, Real *betan, Real ****Pn, Real ****Apn, Real ****Rnstar, Real ****un, Real ****Aun, Real ****yn, Real ****zn, Real ****Arn){

	Real ArndotArn, yndotyn, yndotrn, Arndotrn, Arndotyn, r0dotApn, r0dotrn1, normtemp, dottemp1, dottemp2;
	Real alphan1, qsin1, etan1;
	
#ifdef MPI_PARALLEL
	Real norm[6];
	Real totnorm[6];
#endif

	int i, j, k, m;
	int is, ie, js, je, ks, ke;
	
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;
	
	/* First, calculate Arn, where rn is right  hand side from last step */

		
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				matrix_vector_product3D(theta[k][j][i], phi[k][j][i], psi[k][j][i], varphi[k][j][i], i, j, k, pMat->RHS, Arn[k][j][i]);
				/* do not need to update ghost zones for Arn */
				for(m=0; m<4; m++){					
					Pn[k][j][i][m] = pMat->RHS[k][j][i][m] + (*betan) * (Pn[k][j][i][m] - un[k][j][i][m]);
					Apn[k][j][i][m] = Arn[k][j][i][m] + (*betan) * (Apn[k][j][i][m] - Aun[k][j][i][m]);
			}
	}

	/* now calculate inner product (r0,Apn), (yn, yn), (yn, rn), (yn, Arn), (Arn, Arn),  (Arn, rn)*/
	 ArndotArn = 0.0;
	 yndotyn   = 0.0;
	 yndotrn   = 0.0;
	 Arndotrn  = 0.0;
	 Arndotyn  = 0.0; 
	 r0dotApn  = 0.0;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				vector_product(Arn[k][j][i],Arn[k][j][i], 4 , &normtemp);
				ArndotArn += normtemp;
				vector_product(yn[k][j][i],yn[k][j][i], 4 , &normtemp);
				yndotyn += normtemp;
				vector_product(yn[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
				yndotrn += normtemp;	
				vector_product(Arn[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
				Arndotrn += normtemp;
				vector_product(Arn[k][j][i],yn[k][j][i], 4 , &normtemp);
				Arndotyn += normtemp;
				vector_product(Rnstar[k][j][i],Apn[k][j][i], 4 , &normtemp);
				r0dotApn += normtemp;				
	}

			

	
	/* now collect vector norm from other CPUs */

#ifdef MPI_PARALLEL
	norm[0] = ArndotArn;
	norm[1] = yndotyn;
	norm[2] = yndotrn;
	norm[3] = Arndotrn;
	norm[4] = Arndotyn;
	norm[5] = r0dotApn;
	

	MPI_Allreduce(&norm[0],&totnorm[0],6,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);


	ArndotArn  = totnorm[0];
	yndotyn    = totnorm[1];
	yndotrn    = totnorm[2];
	Arndotrn   = totnorm[3];
	Arndotyn   = totnorm[4];
	r0dotApn   = totnorm[5];
#endif

	

	/* stop the calculation if norm is close to zero */
	dottemp1 = ArndotArn * yndotyn - Arndotyn * Arndotyn;
	dottemp2 = yndotyn * Arndotrn - yndotrn * Arndotyn;

	if(fabs(r0dotApn) < TINY_NUMBER || fabs(dottemp1) < TINY_NUMBER)
		return;

	qsin1   = dottemp2 / dottemp1;

	dottemp2 = yndotrn * ArndotArn - Arndotyn * Arndotrn;

	etan1 = dottemp2 / dottemp1;

	alphan1 = (*r0dotrn) / r0dotApn;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				for(m=0; m<4; m++){
					un[k][j][i][m] = qsin1 * Apn[k][j][i][m] + etan1 * yn[k][j][i][m] + etan1 * (*betan) * un[k][j][i][m];
					zn[k][j][i][m] = qsin1 * pMat->RHS[k][j][i][m] + etan1 * zn[k][j][i][m] - alphan1 * un[k][j][i][m];
				}
	}
	/* need to update bounadry for Aun, as we need A.Aun to calculate un */
	/* need to pMat for information of boundary condition */
	bvals_matrix_vector(un, pMat);

	/* now calculate Aun */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

		matrix_vector_product3D(theta[k][j][i], phi[k][j][i], psi[k][j][i], varphi[k][j][i], i, j, k, un, Aun[k][j][i]);			
		
	}
	

	/* now calculate yn */
	r0dotrn1 = 0.0;
	*currentnorm = 0.0;
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				for(m=0; m<4; m++){
					yn[k][j][i][m] = qsin1 * Arn[k][j][i][m] + etan1 * yn[k][j][i][m] - alphan1 * Aun[k][j][i][m];
					pMat->RHS[k][j][i][m] -= alphan1 * Apn[k][j][i][m];
					pMat->RHS[k][j][i][m] -= yn[k][j][i][m];			
				}
			/* update the guess solution */
			pMat->U[k][j][i].Er  += alphan1 * Pn[k][j][i][0] + zn[k][j][i][0];
			pMat->U[k][j][i].Fr1 += alphan1 * Pn[k][j][i][1] + zn[k][j][i][1];
			pMat->U[k][j][i].Fr2 += alphan1 * Pn[k][j][i][2] + zn[k][j][i][2];
			pMat->U[k][j][i].Fr3 += alphan1 * Pn[k][j][i][3] + zn[k][j][i][3];

			vector_product(Rnstar[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
			r0dotrn1 += normtemp;
			vector_product(pMat->RHS[k][j][i],pMat->RHS[k][j][i], 4 , &normtemp);
			*currentnorm += normtemp;
	}

	
#ifdef MPI_PARALLEL
	norm[0] = r0dotrn1;
	norm[1] = *currentnorm;
	
	MPI_Allreduce(&norm[0],&totnorm[0],2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	

	r0dotrn1     = totnorm[0];
	*currentnorm = totnorm[1];
#endif
	

	if(fabs(r0dotrn1) < TINY_NUMBER)
		(*betan) = 0.0;
	else
		(*betan) = (alphan1/qsin1) * r0dotrn1/(*r0dotrn);

	*r0dotrn = r0dotrn1;
	
	return;

}

/* function to calculate matrix vector productt */
/* calculate A[4][n].vector[n]=result[4] */
/* the vector is stored according to the 3D grid */

void matrix_vector_product3D(Real *theta, Real *phi, Real *psi, Real *varphi, int i, int j, int k, Real ****vector, Real *result)
{

	Real tempEr3, tempEr2, tempEr1, tempFr3, tempFr2, tempFr1, temp0;

	/* The line Er,i,j,k */
	tempEr3 = theta[0] * vector[k-1][j][i][0] + theta[14] * vector[k+1][j][i][0];
	tempEr2 = theta[2] * vector[k][j-1][i][0] + theta[12] * vector[k][j+1][i][0];
	tempEr1 = theta[4] * vector[k][j][i-1][0] + theta[10] * vector[k][j][i+1][0];

	tempFr3 = theta[1] * vector[k-1][j][i][3] + theta[15] * vector[k+1][j][i][3];
	tempFr2 = theta[3] * vector[k][j-1][i][2] + theta[13] * vector[k][j+1][i][2];
	tempFr1 = theta[5] * vector[k][j][i-1][1] + theta[11] * vector[k][j][i+1][1];

	temp0 = theta[7] * vector[k][j][i][1] + theta[8] * vector[k][j][i][2] + theta[9] * vector[k][j][i][3];

	result[0] = ((tempEr1 + tempEr2 + tempEr3) + theta[6] * vector[k][j][i][0]) + tempFr1 + tempFr2 + tempFr3 + temp0;

	

	/* The line Fr1, i,j,k */

	tempEr3 = phi[0] * vector[k-1][j][i][0] + phi[12] * vector[k+1][j][i][0];
	tempEr2 = phi[2] * vector[k][j-1][i][0] + phi[10] * vector[k][j+1][i][0];
	tempEr1 = phi[4] * vector[k][j][i-1][0] + phi[8] * vector[k][j][i+1][0];

	tempFr3 = phi[1] * vector[k-1][j][i][1] + phi[13] * vector[k+1][j][i][1];
	tempFr2 = phi[3] * vector[k][j-1][i][1] + phi[11] * vector[k][j+1][i][1];
	tempFr1 = phi[5] * vector[k][j][i-1][1] + phi[9] * vector[k][j][i+1][1];

	temp0 = phi[6] * vector[k][j][i][0];

			
	result[1] = ((tempEr1 + tempEr2 + tempEr3) + temp0) + (tempFr1 + tempFr2 + tempFr3) + phi[7] * vector[k][j][i][1];

	/* The line Fr2, i,j,k */

	tempEr3 = psi[0] * vector[k-1][j][i][0] + psi[12] * vector[k+1][j][i][0];
	tempEr2 = psi[2] * vector[k][j-1][i][0] + psi[10] * vector[k][j+1][i][0];
	tempEr1 = psi[4] * vector[k][j][i-1][0] + psi[8] * vector[k][j][i+1][0];

	tempFr3 = psi[1] * vector[k-1][j][i][2] + psi[13] * vector[k+1][j][i][2];
	tempFr2 = psi[3] * vector[k][j-1][i][2] + psi[11] * vector[k][j+1][i][2];
	tempFr1 = psi[5] * vector[k][j][i-1][2] + psi[9] * vector[k][j][i+1][2];

	temp0 = psi[6] * vector[k][j][i][0];

			
	result[2] = ((tempEr1 + tempEr2 + tempEr3) + temp0) + (tempFr1 + tempFr2 + tempFr3) + psi[7] * vector[k][j][i][2];

	/* The line Fr3, i,j,k */

	tempEr3 = varphi[0] * vector[k-1][j][i][0] + varphi[12] * vector[k+1][j][i][0];
	tempEr2 = varphi[2] * vector[k][j-1][i][0] + varphi[10] * vector[k][j+1][i][0];
	tempEr1 = varphi[4] * vector[k][j][i-1][0] + varphi[8] * vector[k][j][i+1][0];

	tempFr3 = varphi[1] * vector[k-1][j][i][3] + varphi[13] * vector[k+1][j][i][3];
	tempFr2 = varphi[3] * vector[k][j-1][i][3] + varphi[11] * vector[k][j+1][i][3];
	tempFr1 = varphi[5] * vector[k][j][i-1][3] + varphi[9] * vector[k][j][i+1][3];

	temp0 = varphi[6] * vector[k][j][i][0];

			
	result[3] = ((tempEr1 + tempEr2 + tempEr3) + temp0) + (tempFr1 + tempFr2 + tempFr3) + varphi[7] * vector[k][j][i][3];
	


}


#endif /* radMHD_INTEGRATOR */

#endif /* matrix_multigrid */

