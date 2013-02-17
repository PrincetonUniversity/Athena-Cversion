#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Backward Euler for 3D with SMR. 
 * Solve the matrix for the whole mesh. 
 * Use V cycle, starting from the top level 
 * With SMR, we need to make sure Residual of each level 
 * be smaller than the tolerance
 * CONTAINS PUBLIC FUNCTIONS: 
 *   BackEuler_3d()
 * PRIVATE FUNCTION
 *	INI_matrix()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#ifdef PARTICLES
#include "../particles/particle.h"
#endif


#ifdef STATIC_MESH_REFINEMENT

#ifdef MATRIX_MULTIGRID 

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */


/* 2D array of matrix , corresponding to each domain at each level */
static MatrixS **Matrix;
static int *DomainNos; 

static int coefflag; /* To decide whether need to calculate coefficient or not */






/* The initial L2 norm of the right hand side   *
 * The convergence criterian is
 * ||r_n||_2 / ||r_0||_2  < Small number 	*
 */ 

static int Nlim = 4; /* the lim size of the coarsest grid in each CPU*/
static int Vcyclelim = 15; /* The limit of how many Wcycle we allow */ 

static int Matrixflag = 1; /* The matrix flag is used to choose Gauss_Seidel method or Jacobi method, or Bicgsafe */
			  /* 1 is GS method while 0 is Jacobi method , 2 is Bicgsafe method*/

static int RootLevel; /* Number of levels from below the root domain, log2(N)=Rootlevel */
static int DomainLevels; /* Number of levels above (include) the root domain */
static int TotLevels;

/********Private function ************/
static void Initialize_matrix(MatrixS *pMat, DomainS *pD); /* Function to copy data from grid to the matrix structure, used at the top of each tree */ 

static void RadMHD_multig_3D(MatrixS **Matrix, MeshS *pM);
static void RadMHD_multig_3D_first(MatrixS **Matrix, MeshS *pM);



static void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat);


static void Calculate_Coef(MatrixS *pMat);


/*
void Jacobi3D(MatrixS *pMat, MatrixS *pMatnew);
*/
extern void GaussSeidel3D(MatrixS *pMat, Real ****theta,  Real ****phi,  Real ****psi,  Real ****varphi);

extern void Jacobi3D(MatrixS *pMat, Real ****theta,  Real ****phi,  Real ****psi,  Real ****varphi);

extern void Bicgsafe3D(MatrixS *pMat, Real ****theta,  Real ****phi,  Real ****psi,  Real ****varphi);

/********Public function****************/



void BackEuler_3d(MeshS *pM)
{

	

	DomainS *pD;
	GridS *pG;
	MatrixS *pMat;
	int nd, nl;


	/* set the time */
	Real dt = pM->dt;
	
	/******************************/
	/* Matrix parameters */
	
	Real error, error0, norm;
	int Vcycle;

	int myID;
	
#ifdef MPI_PARALLEL
	myID = myID_Comm_world;
#else
	myID = 0;
#endif


#ifdef FARGO
	Real qom, x1, x2, x3;
	qom = qshear * Omega_0;
#endif

	/******************************/
	int i, j, k;
	int is, ie, js, je, ks, ke;
	int Mati, Matj, Matk;

	Real velocity_x, velocity_y, velocity_z;
	Real Fr0x, Fr0y, Fr0z;
	




	/*****************************/

		
	/* Do the multi-grid W cycle */
	/* reset the flags for each cycle */
	Vcycle = 0;
	error = 1.0;	
	coefflag = 1;
	/* Stop iteration if tolerance level is reached */
	/* Or the matrix doesn't converge in Wcycle limit */
	while((error > TOL) && (Vcycle < Vcyclelim)){
		
	/*	starting from the top Level *
	 * 	But we need go through each domina at this level *
	 */

		/* To advoid if statement in the loop, *
	 	 * First, go from top to bottom, and come back to top for one cycle first. *
		 * During this first cycle, we initialize the matrix and calculate the matrix coefficient */


		/*****************************************************************/
		/* At each level, always restrict the solution to coarse level as initial condition, */
		/* At the same time, restrict the right hand side and calculate the Residual for the coarse level */
		/* during each iteration, put the solution back to Grid, then the solution will be the final solution at the end */

		
		if(coefflag){
			coefflag = 0;
			RadMHD_multig_3D_first(Matrix, pM);

		}
		else{
			RadMHD_multig_3D(Matrix, pM);
		}


		/* Now sum the residual from Domain root level *
		 * to Domain Top level */

		error = 0.0;
		
		/* First, calculate Residual and update the solution for the top level */
		for(nd=0; nd<DomainNos[TotLevels-1]; nd++){

			if(Matrix[TotLevels-1][nd].CPUflag){
			
				RadSMR_Residual3D(&(Matrix[TotLevels-1][nd]),Matrix[TotLevels-1][nd].RHS,&(Matrix[TotLevels-1][nd].RHSnorm));
			
				pG = pM->Domain[TotLevels-1-RootLevel][nd].Grid;
				pMat = &(Matrix[TotLevels-1][nd]);
			
				for(k=pMat->ks; k<=pMat->ke; k++)
				for(j=pMat->js; j<=pMat->je; j++)
				for(i=pMat->is; i<=pMat->ie; i++){
						
					/* Now update the solution */
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Er  += pMat->U[k][j][i].Er;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr1 += pMat->U[k][j][i].Fr1;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr2 += pMat->U[k][j][i].Fr2;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr3 += pMat->U[k][j][i].Fr3;
						
				}/* End i, j, k */
			}/* End if this CPU works for this grid */
		}
				
		
		for(nl=TotLevels-1; nl>RootLevel; nl--){
			/* restrict the solution from top level downwards */
			/* RHS is updated in the Restriction step */
			Rad_Restriction(Matrix, nl, pM, 0);

			for(nd=0; nd<DomainNos[nl]; nd++){				

				if(Matrix[nl][nd].CPUflag){
					norm = Matrix[nl][nd].RHSnorm0;
					if(norm > TINY_NUMBER)
						error0 = Matrix[nl][nd].RHSnorm/norm;
					else
						error0 = Matrix[nl][nd].RHSnorm;

					error += error0;	
				}
			}
			
			
		} 
		/* Now include error from root domain */
		/* Residual for the RootLevel is already updated for te RootLevel in the Restriction step */
		for(nd=0; nd<DomainNos[RootLevel]; nd++){
			if(Matrix[RootLevel][nd].CPUflag){
				norm = Matrix[RootLevel][nd].RHSnorm0;
				if(norm > TINY_NUMBER)
					error0 = Matrix[RootLevel][nd].RHSnorm/norm;
				else
					error0 = Matrix[RootLevel][nd].RHSnorm;

				error += error0;
			}	
		}

		error /= DomainLevels;

		/* The error is the averaged error for the active levels */

		if(error != error)
			ath_error("[BackEuler3D]: NaN encountered!\n");
	
		Vcycle++;
		/* only need to calculate the coefficient in the first cycle */
		coefflag = 0;

		/* If error is smaller than tolerance level, the restricted solution is already at different levels */

	}

		/* Only output the residual for parent grid */	
		if(myID == 0)
			ath_pout(0,"Final residual: %e  Cycle No.: %d\n",error,Vcycle);

		
	/* add the correction back */


	for(nl=0; nl<pM->NLevels; nl++){
		for(nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
			if(Matrix[nl+RootLevel][nd].CPUflag){

				pMat = &(Matrix[nl+RootLevel][nd]);
				pD = &(pM->Domain[nl][nd]);
				pG = pD->Grid;

			
				is = pG->is;
				ie = pG->ie;
				js = pG->js;
				je = pG->je;
				ks = pG->ks;
				ke = pG->ke;

				for(k=ks; k<=ke; k++)
				for(j=js; j<=je; j++)
				for(i=is; i<= ie; i++){
				

					Mati = i - (nghost - Matghost);
					Matj = j - (nghost - Matghost);
					Matk = k - (nghost - Matghost);


					/* correction for work done by radiation force. May need this to reduce energy error */


					velocity_x = pMat->Ugas[Matk][Matj][Mati].V1;
					velocity_y = pMat->Ugas[Matk][Matj][Mati].V2;
					velocity_z = pMat->Ugas[Matk][Matj][Mati].V3;

					Fr0x =  pG->U[k][j][i].Fr1 - ((1.0 +  pG->U[k][j][i].Edd_11) * velocity_x +  pG->U[k][j][i].Edd_21 * velocity_y + pG->U[k][j][i].Edd_31 * velocity_z) *  pG->U[k][j][i].Er / Crat; 
					Fr0y =  pG->U[k][j][i].Fr2 - ((1.0 +  pG->U[k][j][i].Edd_22) * velocity_y +  pG->U[k][j][i].Edd_21 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_z) *  pG->U[k][j][i].Er / Crat;
					Fr0z =  pG->U[k][j][i].Fr3 - ((1.0 +  pG->U[k][j][i].Edd_33) * velocity_z +  pG->U[k][j][i].Edd_31 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_y) *  pG->U[k][j][i].Er / Crat; 

					/* Estimate the added energy source term */
					if(Prat > 0.0){
						pG->U[k][j][i].Er += (pG->Eulersource[k][j][i] - dt * (pG->U[k][j][i].Sigma[1] -  pG->U[k][j][i].Sigma[0]) * ( velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z));
					}

			
			/*			if(pG->U[k][j][i].Er < TINY_NUMBER){
                                			pG->U[k][j][i].Er = (pG->U[k][j][i-1].Er + pG->U[k][j][i+1].Er + pG->U[k][j-1][i].Er + pG->U[k][j+1][i].Er +  pG->U[k-1][j][i].Er + pG->U[k+1][j][i].Er) / 6.0;

                         			}
			*/
			
				}/* End i,j,k */
			}/* End if CPUflag */

		}/* End domain */
	}/* End level */

  return;	
	

}




void RadMHD_multig_3D_first(MatrixS **Matrix, MeshS *pM)
{

	/* We already allocate the memory for each level */
	/* This multigrid will finish one V cycle, starting from the top level */
	/* And return back to the top level */
	/* The array DomainNos stores the number of domains at each level */
	MatrixS *pMat;
	DomainS *pD;
	GridS *pG;

	int nd, nl;

	

	/* Do not do relaxation when going down for restriction */
	/* Do relaxation when going up after prolongation */
	/* [TotLevels-1] is the top level */
	
	/* First, for the levels above root, which has corresponding grid */
	for(nl=TotLevels-1; nl>=RootLevel; nl--){
		/*****************************************/

		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);

			pD = &(pM->Domain[nl-RootLevel][nd]);
			pG = pD->Grid;

			/* Initialize the matrix coefficient and right hand side */
			/* We need to do this for each domain at each level */
			/* Because there can be domains at coarse levels that do not have corresponding fine levels */
			
			/* In this initialize step, background state is subtracted and RHS is updated */
			/* No restriction is needed at this initialization step, as we calculate the */
			/*****************************************************/
			/* Do not need to do restriction in the first time */
			/* The information we need to calculate the matrix coefficient is */
			/* already restricted in the pG grid */
			/*****************************************/
			if(pMat->CPUflag)
				Initialize_matrix(pMat,pD);
		
		}

		
	}

	/* To be consistent with restriction and prolongation */
	/* We restrict the overlap region */
	/* We actually only restrict the RHS here */
	for(nl=TotLevels-1; nl>RootLevel; nl--){
		for(nd=0; nd<DomainNos[nl]; nd++){
			/* restrict the overlap region */
			/* Even this CPU does not work on nl, it still needs to go into this */
			/* As this CPU may need to get data */
			/* We will synchronize all CPUs for restriction */
			pMat = &(Matrix[nl][nd]);
			Rad_Restriction(Matrix, nl, pM, 0);
		}/* Finish domain at level nl */
	}/* Finish looping nl */

	/* Second, for the levels below the root level */
	/* There is only one Domain per Level below the root level */
	/* For these levels, we only need to calculate the matrix coefficient and restrict the matrix parameters */
	/* If this CPU works in this grid at the root domain, it will also work in the levels below that */

	if(Matrix[RootLevel][0].CPUflag){
		
		for(nl=RootLevel-1; nl>=0; nl--){
			/* First do restriction of the matrix coefficient for the whole level */
			Rad_Restriction(Matrix, nl+1, pM, 1);			
			pMat = &(Matrix[nl][0]);
			/* update the boundary for this level */
			/* Only need to update ghost zones for gas quantities */
			/* after first restriction, only right hand side changes */
			bvals_Matrix_gas_init(pMat);
#ifdef SHEARING_BOX
			bvals_Matrix_shear_gas_init(pMat);
#endif
		
			bvals_Matrix_gas(pMat);
		
#ifdef SHEARING_BOX
			bvals_Matrix_shear_gas_destruct();
#endif


			bvals_Matrix_gas_destruct(pMat);

		/* Then calculate the coefficient */
		
			Calculate_Coef(pMat);

		}	
	}

	/* After going down and initialize the matrix for the first time */
	/* Do relaxation and prolongate to the top levels*/

	for(nl=0; nl<TotLevels-1; nl++){
		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);
			if(pMat->CPUflag){
			
				bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_init(pMat);
#endif

		
				if(Matrixflag == 1){
					GaussSeidel3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 0){
					Jacobi3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 2){
					Bicgsafe3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else
					ath_error("[BackEuler3D: unknown matrix solver!\n]");

		
				bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_destruct();
#endif
			}/* End if CPU flag */
		}
		/* prolongate the whole level after relaxation */
		/* We only need to wait for the whole domain to be finished before prolongation */
		Rad_Prolongation(Matrix, nl, pM);
	}


	/* Do relaxation for the Top level */
		nl = TotLevels-1;
		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);
			if(pMat->CPUflag){
			
				bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_init(pMat);
#endif

		
				if(Matrixflag == 1){
					GaussSeidel3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 0){
					Jacobi3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 2){
					Bicgsafe3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else
					ath_error("[BackEuler3D: unknown matrix solver!\n]");

		
				bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_destruct();
#endif
			}/* End if CPU flag */
		}		


	return;
}




void RadMHD_multig_3D(MatrixS **Matrix, MeshS *pM)
{

	/* This is the normal multi-grid cycles. The matrix coefficients are already calculated  */
	/* We do the V cycle: starting from the top, do restriction, no relaxation during this process */
	/* relaxation and prolongation from the bottom level */
	/* Finish one V cycle and return to the top level */


	MatrixS *pMat;

	int nd, nl;



	/* Restriction to the bottom staring from the Root Domain */
	/* When we calculate the Error, we already restrict to the root domain */
	for(nl=RootLevel-1; nl>=0; nl--){
		
		Rad_Restriction(Matrix, nl+1, pM, 0);


	}

	/* Relaxation and prolongation when going up */
	for(nl=0; nl<TotLevels-1; nl++){
		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);
			/* If this CPU works on this grid */
			if(pMat->CPUflag){
				
				bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_init(pMat);
#endif

		
				if(Matrixflag == 1){
					GaussSeidel3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 0){
					Jacobi3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 2){
					Bicgsafe3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else
					ath_error("[BackEuler3D: unknown matrix solver!\n]");

			
				bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_destruct();
#endif



			} /* End if CPUflag */
		} /* End loop over domain in level nl */

		/* prolongate the whole level to the level above */
		Rad_Prolongation(Matrix,nl,pM);
	} /* End loop different levels */
	
	/* Do relaxation for the top level */
		nl = TotLevels-1;
		
		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);
			/* If this CPU works on this grid */
			if(pMat->CPUflag){
				
				bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_init(pMat);
#endif

		
				if(Matrixflag == 1){
					GaussSeidel3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 0){
					Jacobi3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else if(Matrixflag == 2){
					Bicgsafe3D(pMat,pMat->Ptheta,pMat->Pphi,pMat->Ppsi,pMat->Pvarphi);
				}
				else
					ath_error("[BackEuler3D: unknown matrix solver!\n]");

			
				bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
				bvals_Matrix_shear_destruct();
#endif



			} /* End if pMat->CPUflag */
		} /* end loop domain at level nl=TotLevels-1 */

			
			

	/* Exist when we come back to the btop */
	return;
}




/*  calculate the residual of right hand side from guess solution */
/* should allocate memory for newRHS before this function is called */
/* We also need MPI communicator for this Domain */

void RadSMR_Residual3D(MatrixS *pMat, Real ****newRHS, Real *error)
{

	int i, j, k, m;
	int is, ie, js, je, ks, ke;

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;


	
	Real tempEr1, tempEr2, tempEr3;
	Real tempFr1, tempFr2, tempFr3;
	Real temp0;
	Real norm, sum;

#ifdef MPI_PARALLEL
	int ierr;
#endif

	/* To store the coefficient */
	Real theta[16];
	Real phi[16];
	Real psi[16];
	Real varphi[16];

	Real tempRHS[4];


	sum = 0.0;
	
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				/* get the coefficient */
				for(m=0; m<16; m++){
					theta[m]  = pMat->Ptheta[k][j][i][m];
					phi[m] 	  = pMat->Pphi[k][j][i][m];
					psi[m]    = pMat->Ppsi[k][j][i][m];
					varphi[m] = pMat->Pvarphi[k][j][i][m];
				}

			tempRHS[0] = pMat->RHS[k][j][i][0];

			tempEr3 = theta[0] * pMat->U[k-1][j][i].Er + theta[14] * pMat->U[k+1][j][i].Er;
			tempEr2 = theta[2] * pMat->U[k][j-1][i].Er + theta[12] * pMat->U[k][j+1][i].Er;
			tempEr1 = theta[4] * pMat->U[k][j][i-1].Er + theta[10] * pMat->U[k][j][i+1].Er;

			tempFr3 = theta[1] * pMat->U[k-1][j][i].Fr3 + theta[15] * pMat->U[k+1][j][i].Fr3;
			tempFr2 = theta[3] * pMat->U[k][j-1][i].Fr2 + theta[13] * pMat->U[k][j+1][i].Fr2;
			tempFr1 = theta[5] * pMat->U[k][j][i-1].Fr1 + theta[11] * pMat->U[k][j][i+1].Fr1;

			temp0 = theta[7] * pMat->U[k][j][i].Fr1 + theta[8] * pMat->U[k][j][i].Fr2 + theta[9] * pMat->U[k][j][i].Fr3;

			tempRHS[0] -= (theta[6] * pMat->U[k][j][i].Er + (tempEr1 + tempEr2 + tempEr3));
			tempRHS[0] -= ((tempFr1 + tempFr2 + tempFr3) + temp0);

			/********************************************************/


			tempRHS[1] = pMat->RHS[k][j][i][1];


			tempEr3 = phi[0] * pMat->U[k-1][j][i].Er + phi[12] * pMat->U[k+1][j][i].Er;
			tempEr2 = phi[2] * pMat->U[k][j-1][i].Er + phi[10] * pMat->U[k][j+1][i].Er;
			tempEr1 = phi[4] * pMat->U[k][j][i-1].Er + phi[8] * pMat->U[k][j][i+1].Er;

			tempFr3 = phi[1] * pMat->U[k-1][j][i].Fr1 + phi[13] * pMat->U[k+1][j][i].Fr1;
			tempFr2 = phi[3] * pMat->U[k][j-1][i].Fr1 + phi[11] * pMat->U[k][j+1][i].Fr1;
			tempFr1 = phi[5] * pMat->U[k][j][i-1].Fr1 + phi[9] * pMat->U[k][j][i+1].Fr1;

			temp0 = phi[6] * pMat->U[k][j][i].Er;
			
			tempRHS[1] -= (phi[7] * pMat->U[k][j][i].Fr1 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);


			/********************************************************************/


			tempRHS[2] = pMat->RHS[k][j][i][2];


			tempEr3 = psi[0] * pMat->U[k-1][j][i].Er + psi[12] * pMat->U[k+1][j][i].Er;
			tempEr2 = psi[2] * pMat->U[k][j-1][i].Er + psi[10] * pMat->U[k][j+1][i].Er;
			tempEr1 = psi[4] * pMat->U[k][j][i-1].Er + psi[8] * pMat->U[k][j][i+1].Er;

			tempFr3 = psi[1] * pMat->U[k-1][j][i].Fr2 + psi[13] * pMat->U[k+1][j][i].Fr2;
			tempFr2 = psi[3] * pMat->U[k][j-1][i].Fr2 + psi[11] * pMat->U[k][j+1][i].Fr2;
			tempFr1 = psi[5] * pMat->U[k][j][i-1].Fr2 + psi[9] * pMat->U[k][j][i+1].Fr2;

			temp0 = psi[6] * pMat->U[k][j][i].Er;

			
			tempRHS[2] -= (psi[7] * pMat->U[k][j][i].Fr2 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);
			


			/******************************************************/


			tempRHS[3] = pMat->RHS[k][j][i][3];



			tempEr3 = varphi[0] * pMat->U[k-1][j][i].Er + varphi[12] * pMat->U[k+1][j][i].Er;
			tempEr2 = varphi[2] * pMat->U[k][j-1][i].Er + varphi[10] * pMat->U[k][j+1][i].Er;
			tempEr1 = varphi[4] * pMat->U[k][j][i-1].Er + varphi[8] * pMat->U[k][j][i+1].Er;

			tempFr3 = varphi[1] * pMat->U[k-1][j][i].Fr3 + varphi[13] * pMat->U[k+1][j][i].Fr3;
			tempFr2 = varphi[3] * pMat->U[k][j-1][i].Fr3 + varphi[11] * pMat->U[k][j+1][i].Fr3;
			tempFr1 = varphi[5] * pMat->U[k][j][i-1].Fr3 + varphi[9] * pMat->U[k][j][i+1].Fr3;

			temp0 = varphi[6] * pMat->U[k][j][i].Er;

			
			tempRHS[3] -= (varphi[7] * pMat->U[k][j][i].Fr3 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

			if(newRHS != NULL){
				/* Usually, newRHS will be pMat->RHS, so RHS is directly updated here */
				newRHS[k][j][i][0] = tempRHS[0];
				newRHS[k][j][i][1] = tempRHS[1];
				newRHS[k][j][i][2] = tempRHS[2];
				newRHS[k][j][i][3] = tempRHS[3];

			}

			sum += (tempRHS[0] * tempRHS[0] + tempRHS[1] * tempRHS[1] + tempRHS[2] * tempRHS[2] + tempRHS[3] * tempRHS[3]);

	}

		/* error = b - Ax */

	/* Note that this norm is  A * A, we do not take square root yet */

	if(error != NULL){
#ifdef MPI_PARALLEL
	
		ierr = MPI_Allreduce(&sum,&norm,1,MPI_DOUBLE,MPI_SUM,pMat->Comm_Domain);
#else
		norm = sum;
#endif	

		/* usually, *error is pMat->RHSnorm, which stores the norm of current RHS */
		*error = norm;
	}


	return;

}

/* Calculate the matrix coefficient for each level */
void Calculate_Coef(MatrixS *pMat)
{
	int i, j, k;
	int is, ie, js, je, ks, ke;
	

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

				matrix_coef(pMat, NULL, 3, i, j, k, 0.0, &(pMat->Ptheta[k][j][i][0]), &(pMat->Pphi[k][j][i][0]), &(pMat->Ppsi[k][j][i][0]), &(pMat->Pvarphi[k][j][i][0]));
		}

}


/* This function only sets the levels below the root level */
void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat)
{

#ifdef MPI_PARALLEL
        pMat_coarse->Comm_Domain = pMat->Comm_Domain;      /*!< MPI communicator between Grids on this Dom */
#endif

	pMat_coarse->dx1 = 2.0 * pMat->dx1;
	pMat_coarse->dx2 = 2.0 * pMat->dx2;
	pMat_coarse->dx3 = 2.0 * pMat->dx3;
	pMat_coarse->time = pMat->time;
	pMat_coarse->Lx = pMat->Lx;
	pMat_coarse->Ly = pMat->Ly;
	pMat_coarse->Lz = pMat->Lz;
	/* Grid numbers decrease by a factor of 2 */
	pMat_coarse->Nx[0] = pMat->Nx[0] / 2;
	pMat_coarse->Nx[1] = pMat->Nx[1] / 2;
	pMat_coarse->Nx[2] = pMat->Nx[2] / 2;
	/* Level label decreases when going down the stair */
	pMat_coarse->Level = pMat->Level - 1;
	/* top level grid number doesn't change */
	pMat_coarse->RootNx[0] = pMat->RootNx[0];
	pMat_coarse->RootNx[1] = pMat->RootNx[1];
	pMat_coarse->RootNx[2] = pMat->RootNx[2];
	pMat_coarse->MinX[0] = pMat->MinX[0];
	pMat_coarse->MinX[1] = pMat->MinX[1];
	pMat_coarse->MinX[2] = pMat->MinX[2];
	pMat_coarse->MaxX[0] = pMat->MaxX[0];
	pMat_coarse->MaxX[1] = pMat->MaxX[1];
	pMat_coarse->MaxX[2] = pMat->MaxX[2];
	pMat_coarse->NGrid[0] = pMat->NGrid[0];
	pMat_coarse->NGrid[1] = pMat->NGrid[1];
	pMat_coarse->NGrid[2] = pMat->NGrid[2];
	pMat_coarse->dt = pMat->dt;
	pMat_coarse->is = Matghost;
	pMat_coarse->ie = pMat_coarse->Nx[0] + Matghost - 1;
	pMat_coarse->js = Matghost;
	pMat_coarse->je = pMat_coarse->Nx[1] + Matghost - 1;
	pMat_coarse->ks = Matghost;
	pMat_coarse->ke = pMat_coarse->Nx[2] + Matghost - 1;
	
	/* For these levels, they are always solved with the same CPU */	
	pMat_coarse->ID = pMat->ID;
	pMat_coarse->my_iproc = pMat->my_iproc;
	pMat_coarse->my_jproc = pMat->my_jproc; 
	pMat_coarse->my_kproc = pMat->my_kproc;

	pMat_coarse->bgflag  = pMat->bgflag;
	pMat_coarse->CPUflag = pMat->CPUflag;

	pMat_coarse->rx1_id = pMat->rx1_id;
	pMat_coarse->lx1_id = pMat->lx1_id;
	pMat_coarse->rx2_id = pMat->rx2_id;
	pMat_coarse->lx2_id = pMat->lx2_id;
	pMat_coarse->lx3_id = pMat->lx3_id;
	pMat_coarse->rx3_id = pMat->rx3_id;

	/* Boundary flag */
	pMat_coarse->BCFlag_ix1 = pMat->BCFlag_ix1;
	pMat_coarse->BCFlag_ox1 = pMat->BCFlag_ox1;
	pMat_coarse->BCFlag_ix2 = pMat->BCFlag_ix2;
	pMat_coarse->BCFlag_ox2 = pMat->BCFlag_ox2;
	pMat_coarse->BCFlag_ix3 = pMat->BCFlag_ix3;
	pMat_coarse->BCFlag_ox3 = pMat->BCFlag_ox3;
	



	return;
}




/*-------------------------------------------------------------------------*/
/* BackEuler_init_3d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_3d(): function to free memory */
void BackEuler_init_3d(MeshS *pM)
{

	/* Initialization needs to be done for each domain and each level */
	/* Plus, the coarse grids are also included in the refinement levels */

	DomainS *pD;
	GridS *pG;
	MatrixS *pMat;
	
	int nd, nl, irefine, i;
	int Nx, Ny, Nz;
	Real temp;


	

	DomainLevels = pM->NLevels; /* Total levels of refinement above the root domain */

	/* Now consider the coarse grid below the root domain */
	/* First check that root level only has one domain */
	if(pM->DomainsPerLevel[0] > 1)
		ath_error("[BackEuler_init_3d]: Root Level can only have one domain!\n");
	/* We require grids at the root domain have the same size */
	pD = &(pM->Domain[0][0]);

	if(((pD->Nx[0] % pD->NGrid[0]) != 0) || ((pD->Nx[1] % pD->NGrid[1]) != 0) || ((pD->Nx[2] % pD->NGrid[2]) != 0))
		ath_error("[BackEuler_init_3d]: We require grids on the root domain have the same size!\n");

	Nx = pD->Nx[0] / pD->NGrid[0];
	Ny = pD->Nx[1] / pD->NGrid[1];
	Nz = pD->Nx[2] / pD->NGrid[2];

	/* Reach bottom first for the side with the smallest size*/
	RootLevel = Nx;
	if(Ny < Nx) RootLevel = Ny;
	if(Nz < Nx) RootLevel = Nz;

	RootLevel /= Nlim;
	
	temp = log10(RootLevel)/log10(2.0);

	RootLevel = (int)temp;
	if(fabs(temp-RootLevel) > 0.5) RootLevel++;

	TotLevels = DomainLevels + RootLevel;
	
	/* RootLevel is the number of levels below root domain, does not include root domain */
	/* Then Matrix[RootLevel] will be the root domain */
	/****************************************/

	/* The array to store number of domains at each level */
	if((DomainNos = (int*)calloc(DomainLevels+RootLevel,sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");
	
	for(nl=0; nl<DomainLevels; nl++){
		DomainNos[RootLevel+nl] = pM->DomainsPerLevel[nl];		
	}

	for(nl=0; nl<RootLevel; nl++){
		DomainNos[nl] = 1;
	}

	
	/****************************************/
	/* Memory for the arrays of matrix */
	if((Matrix = (MatrixS**)calloc(DomainLevels+RootLevel,sizeof(MatrixS *))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");
	

	/* Matrix[RootLevel] corresponds to the root domain Domain[0] */
	for(nl=0; nl<DomainLevels; nl++){
		nd = pM->DomainsPerLevel[nl];
		if((Matrix[nl+RootLevel] = (MatrixS*)calloc(nd,sizeof(MatrixS))) == NULL)
			ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");
	}
	/* Below the root level, one level only has one domain */
	for(nl=0; nl<RootLevel; nl++){
		if((Matrix[nl] = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
			ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");
	}

	/****************************************/



	/****************************************/
	/* Now initialize the Matrix array for each Grid */
	for(nl=0; nl<DomainLevels; nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			pD = &(pM->Domain[nl][nd]);
			pG = pD->Grid;
			/* Only do this if this CPU has grid */
			if(pG != NULL){
				irefine = 1; /* The ratio of resolution between this level and root level */
				for(i=1; i<=nl; i++)
					irefine *= 2;

				Nx = pG->ie - pG->is + 1;
				Ny = pG->je - pG->js + 1;
				Nz = pG->ke - pG->ks + 1;
		
				/* pMat will remain in the memory until the end of the simulation */
				pMat = &(Matrix[nl+RootLevel][nd]);
				

				if((pMat->U = (RadMHDS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadMHDS))) == NULL)
					ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

				if((pMat->Ugas = (RadCoefS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadCoefS))) == NULL)
					ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

				if((pMat->RHS = (Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,4,sizeof(Real))) == NULL)
					ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


				if((pMat->Ptheta=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
					ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

				if((pMat->Pphi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
					ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");


				if((pMat->Ppsi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
					ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

				if((pMat->Pvarphi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
					ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");




				/*==========================================================================*/
				/* now set the parameters */
#ifdef MPI_PARALLEL
  				pMat->Comm_Domain = pD->Comm_Domain;
#endif

				pMat->dx1 = pG->dx1;
				pMat->dx2 = pG->dx2;
				pMat->dx3 = pG->dx3;
				pMat->time = pG->time;
				/* dt in pG is not set at this time */
				/* pMat->dt = pG->dt;	*/
				pMat->Lx = pD->MaxX[0] - pD->MinX[0];
				pMat->Ly = pD->MaxX[1] - pD->MinX[1];
				pMat->Lz = pD->MaxX[2] - pD->MinX[2];
				pMat->CPUflag = 1;

				pMat->is = Matghost;
				pMat->ie = Nx + Matghost - 1;
				pMat->js = Matghost;
				pMat->je = Ny + Matghost - 1;
				pMat->ks = Matghost;
				pMat->ke = Nz + Matghost - 1;
				pMat->Nx[0] = Nx;
				pMat->Nx[1] = Ny;
				pMat->Nx[2] = Nz;				
				pMat->Level = nl + RootLevel;
				pMat->MinX[0] = pG->MinX[0];
				pMat->MinX[1] = pG->MinX[1];
				pMat->MinX[2] = pG->MinX[2];
				pMat->MaxX[0] = pG->MaxX[0];
				pMat->MaxX[1] = pG->MaxX[1];
				pMat->MaxX[2] = pG->MaxX[2];
				pMat->NGrid[0] = pD->NGrid[0];
				pMat->NGrid[1] = pD->NGrid[1];
				pMat->NGrid[2] = pD->NGrid[2];
				pMat->rx1_id = pG->rx1_id;
				pMat->lx1_id = pG->lx1_id;
				pMat->rx2_id = pG->rx2_id;
				pMat->lx2_id = pG->lx2_id;
				pMat->lx3_id = pG->lx3_id;
				pMat->rx3_id = pG->rx3_id;
#ifdef MPI_PARALLEL
				pMat->ID = myID_Comm_world;
				get_myGridIndex(pD, myID_Comm_world, &(pMat->my_iproc), &(pMat->my_jproc), &(pMat->my_kproc));
#else
				pMat->ID = 0;
				pMat->my_iproc = 0;
				pMat->my_jproc = 0;
				pMat->my_kproc = 0;
#endif


				/* Boundary flag */
				/* If the boundary does not touch the boundary of the root domain, ghost zones should */
				/* be set by prolongation, not by boundary function */
 
				if(pD->Disp[0] != 0)
					pMat->BCFlag_ix1 = 0;
				else
					pMat->BCFlag_ix1 = pM->BCFlag_ix1;

				if(((pD->Disp[0] + pD->Nx[0])/irefine) != pM->Nx[0])
					pMat->BCFlag_ox1 = 0;
				else
					pMat->BCFlag_ox1 = pM->BCFlag_ox1;	

				if(pD->Disp[1] != 0)
					pMat->BCFlag_ix2 = 0;	
				else
					pMat->BCFlag_ix2 = pM->BCFlag_ix2;

				if(((pD->Disp[1] + pD->Nx[1])/irefine) != pM->Nx[1])
					pMat->BCFlag_ox2 = 0;
				else
					pMat->BCFlag_ox2 = pM->BCFlag_ox2;

				if(pD->Disp[2] != 0)
					pMat->BCFlag_ix3 = 0;
				else
					pMat->BCFlag_ix3 = pM->BCFlag_ix3;


				if((pD->Disp[2] + pD->Nx[2])/irefine != pM->Nx[2])
					pMat->BCFlag_ox3 = 0;
				else
					pMat->BCFlag_ox3 = pM->BCFlag_ox3;

				/* To decide whether subtract background solution at top level or not */
				/* Default choice is not */
				pMat->bgflag = 1;

			/*==========================================================================*/

			}
			else{
				/* set the INIflag =-1, so that we know this CPU does not work on this domain at this level */
				Matrix[nl+RootLevel][nd].CPUflag = 0;
			}
		}
	}


	/****************************************/



	/****************************************/
	/* Now we need to set the levels below the root level */
	/* For these levels, each level has one domain */
	
	if(Matrix[RootLevel][0].CPUflag){
		for(nl=RootLevel; nl>0; nl--){
			/* First, we need to check whether this CPU works in this grid of the root domain */
			/* If this CPU works here, then it also works for all other coarse grids below the root level */



			/* set the coarse grid based on the fine grid */
			set_mat_level(&(Matrix[nl-1][0]), &(Matrix[nl][0]));

			/* Now allocate memory for these coarse levels */
			pMat = &(Matrix[nl-1][0]);

			Nx = pMat->Nx[0];
			Ny = pMat->Nx[1];
			Nz = pMat->Nx[2];


			if((pMat->U = (RadMHDS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadMHDS))) == NULL)
				ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

			if((pMat->Ugas = (RadCoefS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadCoefS))) == NULL)
				ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

			if((pMat->RHS = (Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,4,sizeof(Real))) == NULL)
				ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


			if((pMat->Ptheta=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
				ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

			if((pMat->Pphi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
				ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");


			if((pMat->Ppsi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
				ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

			if((pMat->Pvarphi=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,16,sizeof(Real))) == NULL)
				ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

		}
	}
	else{
		for(nl=RootLevel-1; nl>=0; nl--){
			/* set the coarse grid based on the fine grid */
			Matrix[nl][0].CPUflag = 0;
		}
	}


	/* Now initialize memory for Restriction and Prolongation */
	/* Static mesh refinement is already defined */
	SMR_Rad_init(pM, RootLevel); 


	return;
}


void BackEuler_destruct_3d(MeshS *pM)
{
	int nl, nd;
	MatrixS *pMat;

	/*********************************************/

	for(nl=0; nl<pM->NLevels; nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			pMat = &(Matrix[nl+RootLevel][nd]);
			
			if(pMat->U != NULL)
				free_3d_array(pMat->U);

			if(pMat->Ugas != NULL)
				free_3d_array(pMat->Ugas);

			if(pMat->RHS != NULL)
				free_4d_array(pMat->RHS);

			if(pMat->Ptheta != NULL)
				free_4d_array(pMat->Ptheta);

			
			if(pMat->Pphi != NULL)
				free_4d_array(pMat->Pphi);

			if(pMat->Ppsi != NULL)
				free_4d_array(pMat->Ppsi);

			
			if(pMat->Pvarphi != NULL)
				free_4d_array(pMat->Pvarphi);


		}

		/* Free this level */
		if(Matrix[nl+RootLevel] != NULL)
			free(Matrix[nl+RootLevel]);
	}
	/*****************************************/

	/****************************/
	/* Now free the root domain below */
	for(nl=0; nl<RootLevel; nl++){
		pMat = &(Matrix[nl][0]);
			
		if(pMat->U != NULL)
			free_3d_array(pMat->U);

		if(pMat->Ugas != NULL)
			free_3d_array(pMat->Ugas);

		if(pMat->RHS != NULL)
			free_4d_array(pMat->RHS);

		if(pMat->Ptheta != NULL)
			free_4d_array(pMat->Ptheta);

			
		if(pMat->Pphi != NULL)
			free_4d_array(pMat->Pphi);

		if(pMat->Ppsi != NULL)
			free_4d_array(pMat->Ppsi);

			
		if(pMat->Pvarphi != NULL)
			free_4d_array(pMat->Pvarphi);


		/* Free this level */
		if(Matrix[nl] != NULL)
			free(Matrix[nl]);

	}	
	/****************************/

	/* Now free the arrays of matrix */
		if(Matrix != NULL)
			free(Matrix);	

		if(DomainNos != NULL)
			free(DomainNos);

	


	/* Free the space for restriction and prolongation */
	SMR_Rad_destruct();

	

	return;

}


void Initialize_matrix(MatrixS *pMat, DomainS *pD)
{

	GridS *pG = pD->Grid;		
	
	int i, j, k, Nx, Ny, Nz;
	int is, ie, js, je, ks, ke;
	int Mati, Matj, Matk;

#ifdef FARGO
	double x1, x2, x3, qom;
	qom = qshear * Omega_0;

#endif
	Real velocity_x, velocity_y, velocity_z, T4, Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real AdvFx, AdvFy, AdvFz;

	/*Sigma_aF: flux mean absorption, Sigma_aP: Plank mean absorption, Sigma_aE: Er mean absorption opacity;*/

	Nx = pG->ie - pG->is + 1;
	Ny = pG->je - pG->js + 1;
	Nz = pG->ke - pG->ks + 1;


	Real dt = pG->dt;
	pMat->dt = dt;
	pMat->time = pG->time;


	/* Set the boundary */
	is = pG->is-Matghost;
	ie = pG->ie+Matghost;
	js = pG->js-Matghost;
	je = pG->je+Matghost;
	ks = pG->ks-Matghost;
	ke = pG->ke+Matghost;

	/* Step 1: Setup the parameters for all the domains at the top level  */
	/* Now copy the data */
	/* Including the ghost zones */
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){
			for(i=is; i<= ie; i++){
				

				Mati = i - (nghost - Matghost);
				Matj = j - (nghost - Matghost);
				Matk = k - (nghost - Matghost);

				velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
				velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
				velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;

				

#ifdef FARGO
				cc_pos(pG,i,j,k,&x1,&x2,&x3);
				velocity_y -= qom * x1;				
#endif


				/* Tguess is already T4 */
				T4 = pG->Tguess[k][j][i];

				Sigma_sF = pG->U[k][j][i].Sigma[0];
				Sigma_aF = pG->U[k][j][i].Sigma[1];
				Sigma_aP = pG->U[k][j][i].Sigma[2];
				Sigma_aE = pG->U[k][j][i].Sigma[3];


				pMat->U[Matk][Matj][Mati].Er  = pG->U[k][j][i].Er;				
				pMat->U[Matk][Matj][Mati].Fr1 = pG->U[k][j][i].Fr1;
				pMat->U[Matk][Matj][Mati].Fr2 = pG->U[k][j][i].Fr2;
				pMat->U[Matk][Matj][Mati].Fr3 = pG->U[k][j][i].Fr3;
				/* Store the background Er at time step n to each level in multigrid */
				pMat->Ugas[Matk][Matj][Mati].rho = pG->U[k][j][i].Er;
				pMat->Ugas[Matk][Matj][Mati].V1  = velocity_x;
				pMat->Ugas[Matk][Matj][Mati].V2  = velocity_y;
				pMat->Ugas[Matk][Matj][Mati].V3  = velocity_z;
				pMat->Ugas[Matk][Matj][Mati].T4  = T4;
				pMat->Ugas[Matk][Matj][Mati].Edd_11 = pG->U[k][j][i].Edd_11;
				pMat->Ugas[Matk][Matj][Mati].Edd_21 = pG->U[k][j][i].Edd_21;
				pMat->Ugas[Matk][Matj][Mati].Edd_22 = pG->U[k][j][i].Edd_22;
				pMat->Ugas[Matk][Matj][Mati].Edd_31 = pG->U[k][j][i].Edd_31;
				pMat->Ugas[Matk][Matj][Mati].Edd_32 = pG->U[k][j][i].Edd_32;
				pMat->Ugas[Matk][Matj][Mati].Edd_33 = pG->U[k][j][i].Edd_33;
				pMat->Ugas[Matk][Matj][Mati].Sigma[0] = Sigma_sF;
				pMat->Ugas[Matk][Matj][Mati].Sigma[1] = Sigma_aF;
				pMat->Ugas[Matk][Matj][Mati].Sigma[2] = Sigma_aP;
				pMat->Ugas[Matk][Matj][Mati].Sigma[3] = Sigma_aE;
			
				

		/* Now set the right hand side */
				Rad_Advection_Flux3D(pD, i, j, k, 1.0, &AdvFx, &AdvFy, &AdvFz);
						
				pMat->RHS[Matk][Matj][Mati][0] = pG->U[k][j][i].Er + dt * Sigma_aP * T4 * Crat * Eratio + (1.0 - Eratio) * pG->Ersource[k][j][i] + (AdvFx + AdvFy + AdvFz) + pG->Comp[k][j][i];
				pMat->RHS[Matk][Matj][Mati][1] = pG->U[k][j][i].Fr1 + Eratio * dt * Sigma_aP * T4 * velocity_x + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_x / Crat;
				pMat->RHS[Matk][Matj][Mati][2] = pG->U[k][j][i].Fr2 + Eratio * dt * Sigma_aP * T4 * velocity_y + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_y / Crat;
				pMat->RHS[Matk][Matj][Mati][3] = pG->U[k][j][i].Fr3 + Eratio * dt * Sigma_aP * T4 * velocity_z + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_z / Crat;	

					
			} /* End i */
		}/* End j */
	}/* End k */


	/* Step 2, calculate the coefficient for each matrix  */
	/* coefficient will not change until next time step */
	Calculate_Coef(pMat);



	/* Step 3: calculate the residual */
	/* With SMR, we always subtract the background state */
	RadSMR_Residual3D(pMat, pMat->RHS, &(pMat->RHSnorm0));

	/* Now set the matrix to be the residual */
	/* RHS is already updated with the new RHS in function RadSMR_Residual3D */
	for(k=pMat->ks - Matghost; k<=pMat->ke+Matghost; k++){
		for(j=pMat->js-Matghost; j<=pMat->je+Matghost; j++){
			for(i=pMat->is-Matghost; i<= pMat->ie+Matghost; i++){
				pMat->U[k][j][i].Er = 0.0;
				pMat->U[k][j][i].Fr1 = 0.0;
				pMat->U[k][j][i].Fr2 = 0.0;
				pMat->U[k][j][i].Fr3 = 0.0;

			}
		}
	}
	/* Right hand side in the ghost zones are never used */


}




#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */

#endif /* STATIC_MESH_REFINEMENT */





