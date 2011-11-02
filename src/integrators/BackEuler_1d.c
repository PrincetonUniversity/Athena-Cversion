#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Use backward Euler method to update the radiation quantities
 * First set up the matrix
 * Then solve the matrix equations.
 * We need the flag for boundary condition.
 *
 * Backward Euler should be used for the whole mesh
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   BackEuler_1d()
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

#if defined(MATRIX_MULTIGRID) || defined(MATRIX_HYPRE)

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */



static MatrixS *pMat;
static Real ****INIerror;


static int Nlim = 4; /* the lim size of the coarsest grid in each CPU*/
static int Wcyclelim = 10; /* The limit of how many Wcycle we allow */ 

static int Nlevel; /* Number of levels from top to bottom, log2(N)=Nlevel */
static int Wflag; /* To decide the position in the W flag */


/********Private function ************/
static Real CheckResidual(MatrixS *pMat, GridS *pG);
static void RadMHD_multig_1D(MatrixS *pMat);
static void prolongation1D(MatrixS *pMat_coarse, MatrixS *pMat_fine);
static void Restriction1D(MatrixS *pMat_fine, MatrixS *pMat_coarse);

static void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat);
static void RHSResidual1D(MatrixS *pMat, Real ****newRHS);

/* Matrix boundary function */
extern void bvals_Matrix_init(MatrixS *pMat);
extern void bvals_Matrix(MatrixS *pMat);
extern void bvals_Matrix_destruct(MatrixS *pMat);

/*
void Jacobi3D(MatrixS *pMat, MatrixS *pMatnew);
*/
extern void GaussSeidel1D(MatrixS *pMat);

/********Public function****************/
/*-------BackEuler_3d(): Use back euler method to update E_r and Fluxr-----------*/
/* Not work with SMR now. So there is only one Domain in the Mesh */


void BackEuler_1d(MeshS *pM)
{



	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	
	Real dt = pG->dt;

	/* Set the parameters that will change with time in pMat */
	pMat->dt = dt;
	pMat->time = pG->time;

	Real velocity_x, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF, pressure;

	Real error;
	int Wcycle;
	
	int i, j;
	int is, ie, js, ks;
	int Mati, Matj, Matk;
	/* Set the boundary */
	is = pG->is-Matghost;
	ie = pG->ie+Matghost;
	js = pG->js;
	ks = pG->ks;
	Matk = ks;
	Matj = js;

	
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
	

	/* Now copy the data */

			for(i=is; i<= ie; i++){
				

				Mati = i - (nghost - Matghost);

				
				velocity_x = pG->U[ks][js][i].M1 / pG->U[ks][js][i].d;

				/* Now Backward Euler step is done after the gas quantities are updated */
				pressure = (pG->U[ks][js][i].E - 0.5 * pG->U[ks][js][i].d * velocity_x * velocity_x) * (Gamma - 1.0);
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[ks][js][i].B1c * pG->U[ks][js][i].B1c) * (Gamma - 1.0);
#endif

				if(pressure < TINY_NUMBER){
					T4 = pG->U[ks][js][i].Er;
				}
				else{
					T4 = pow((pressure / (pG->U[ks][js][i].d * R_ideal)), 4.0);
				}

			
				
				Sigma_sF = pG->U[ks][js][i].Sigma[0];
				Sigma_aF = pG->U[ks][js][i].Sigma[1];
				Sigma_aP = pG->U[ks][js][i].Sigma[2];
				Sigma_aE = pG->U[ks][js][i].Sigma[3];

				pMat->U[Matk][Matj][Mati].Er  = pG->U[ks][js][i].Er;
				pMat->U[Matk][Matj][Mati].Fr1 = pG->U[ks][js][i].Fr1;
				pMat->U[Matk][Matj][Mati].V1  = velocity_x;

				pMat->U[Matk][Matj][Mati].T4  = T4;
				pMat->U[Matk][Matj][Mati].Edd_11 = pG->U[ks][js][i].Edd_11;


				pMat->U[Matk][Matj][Mati].Sigma[0] = Sigma_sF;
				pMat->U[Matk][Matj][Mati].Sigma[1] = Sigma_aF;
				pMat->U[Matk][Matj][Mati].Sigma[2] = Sigma_aP;
				pMat->U[Matk][Matj][Mati].Sigma[3] = Sigma_aE;

		/* Now set the right hand side */
				pMat->RHS[Matk][Matj][Mati][0] = pG->U[ks][js][i].Er + pG->Tguess[ks][js][i];
				pMat->RHS[Matk][Matj][Mati][1] = pG->U[ks][js][i].Fr1 + dt * Sigma_aP * T4 * velocity_x;
	
				
	} /* End i */

if(pMat->bgflag){
	/* calculate the residual */
	RHSResidual1D(pMat, INIerror);

	for(i=pMat->is-Matghost; i<= pMat->ie+Matghost; i++){
		pMat->U[ks][js][i].Er = 0.0;
		pMat->U[ks][js][i].Fr1 = 0.0;


		pMat->RHS[ks][js][i][0] = INIerror[ks][js][i][0];	
		pMat->RHS[ks][js][i][1] = INIerror[ks][js][i][1];	
	}
}
		
	/* Do the multi-grid W cycle */
	Wcycle = 0;
	error = 1.0;
	/* Stop iteration if tolerance level is reached */
	/* Or the matrix doesn't converge in Wcycle limit */
	while((error > TOL) && (Wcycle < Wcyclelim)){
		Wflag = 1;

		RadMHD_multig_1D(pMat);

		/* for parent grid, check residual */
		error = fabs(CheckResidual(pMat,pG));
		if(error != error)
			ath_error("[BackEuler3D]: NaN encountered!\n");
	
		Wcycle++;


	}
		/* Only output the residual for parent grid */	
		if(myID == 0)
			ath_pout(0,"Final residual: %e  Cycle No.: %d\n",error,Wcycle);

		
	/* copy the data back */
	/* Now copy the data */
			for(i=is; i<= ie; i++){
				

				Mati = i - (nghost - Matghost);

			if(pMat->bgflag){
				pG->U[ks][js][i].Er += pMat->U[Matk][Matj][Mati].Er;
				pG->U[ks][js][i].Fr1 += pMat->U[Matk][Matj][Mati].Fr1;
			}
			else{
				pG->U[ks][js][i].Er = pMat->U[Matk][Matj][Mati].Er;
				pG->U[ks][js][i].Fr1 = pMat->U[Matk][Matj][Mati].Fr1;
			}
				
	}

	/* Set the boundary condition */

	
/* Update the ghost zones for different boundary condition to be used later */
	for (i=0; i<pM->NLevels; i++){ 
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
        	if (pM->Domain[i][j].Grid != NULL){
  			bvals_radMHD(&(pM->Domain[i][j]));

        	}
      	     }
    	}

  return;	
	

}


void RadMHD_multig_1D(MatrixS *pMat)
{

	MatrixS pMat_coarse;
	int Nsmall;

	/* Once any dimension reaches size limit of Nlim, do Ncycle iteration and return */
	if(pMat->Nx[0] == Nlim){
		/* Create the temporary array */
		/* Need to create the temporary array for the boundary condition */		


		bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
		bvals_Matrix_shear_init(pMat);
#endif

		
		
		GaussSeidel1D(pMat);

		
		bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
		bvals_Matrix_shear_destruct();
#endif

		

	}
	else{

		/* First, do relaxation at the fine level during the process of going down */
		/*
		bvals_Matrix_init(pMat);
		
		GaussSeidel3D(pMat);
		
		bvals_Matrix_destruct(pMat);
		*/

		/* create the coarse grid */
		set_mat_level(&(pMat_coarse), pMat);

		
		/* project the data to coarse grid */
		Restriction1D(pMat, &pMat_coarse);

		bvals_Matrix_init(&pMat_coarse);
#ifdef SHEARING_BOX
		bvals_Matrix_shear_init(&pMat_coarse);
#endif
		
		bvals_Matrix(&pMat_coarse);
		
#ifdef SHEARING_BOX
		bvals_Matrix_shear_destruct();
#endif


		bvals_Matrix_destruct(&pMat_coarse);

		/* continue the V or W cycle recursively*/

		RadMHD_multig_1D(&pMat_coarse);

/* The following code is first reached after Ncycle iterations at the coarsest
 * level.  We then prolongate, do Ncycle iterations, and return.  This will return
 * execution to this same spot for the next coarsest level, so we will
 * prolongate, do Ncycle iterations, return, and so on.
 */

		/* Add the correction back to the fine grid */
		prolongation1D(&pMat_coarse, pMat);

		/* First destroy the coarse grid */

		free_3d_array(pMat_coarse.U);
		free_4d_array(pMat_coarse.RHS);

		/* Do relaxation when going up*/
			
		
		bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
		bvals_Matrix_shear_init(pMat);
#endif
		
		/* Update the ghost zones first */
		bvals_Matrix(pMat);

		GaussSeidel1D(pMat);
		
#ifdef SHEARING_BOX
		bvals_Matrix_shear_destruct();
#endif


		bvals_Matrix_destruct(pMat);
	
		/* To decide whether go W cycle */
		Nsmall = pMat->Nx[0];


		if((Nsmall > pow(2.0,Nlevel/2.0)) && Wflag){
			/* Set is Wcycle flag */
			/* This is Wcycle. We should go down again */
			Wflag = 0;	

			set_mat_level(&(pMat_coarse), pMat);

		
			/* project the data to coarse grid */
			Restriction1D(pMat, &pMat_coarse);

			/* We need to update the ghost zones after restriction */
			bvals_Matrix_init(&pMat_coarse);

#ifdef SHEARING_BOX
			bvals_Matrix_shear_init(&pMat_coarse);
#endif
		
			bvals_Matrix(&pMat_coarse);
		
#ifdef SHEARING_BOX
			bvals_Matrix_shear_destruct();
#endif


			bvals_Matrix_destruct(&pMat_coarse);


			/* continue the V or W cycle recursively*/

			RadMHD_multig_1D(&pMat_coarse);


			/* Add the correction back to the fine grid */
			prolongation1D(&pMat_coarse, pMat);

			/* First destroy the coarse grid */

			free_3d_array(pMat_coarse.U);
			free_4d_array(pMat_coarse.RHS);

			/* Do relaxation when going up*/
			
		
			bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
			bvals_Matrix_shear_init(pMat);
#endif

			/* Update the ghost zones first */
			bvals_Matrix(pMat);
		
			GaussSeidel1D(pMat);
		
#ifdef SHEARING_BOX
			bvals_Matrix_shear_destruct();
#endif


			bvals_Matrix_destruct(pMat);

		}


		


	}


	return;
}






/*-------------------------------------------------------------------------*/
/* BackEuler_init_1d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_1d(): function to free memory */
void BackEuler_init_1d(MeshS *pM)
{

	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	int Nx, Ny, Nz;
	

	Nx = pG->ie - pG->is + 1;
	Ny = 1;
	Nz = 1;

	/* Reach bottom first for the side with the smallest size*/
	Nlevel = Nx;



	Nlevel = (int)(log10(Nlevel)/log10(2.0));

	/* pMat will remain in the memory until the end of the simulation */

	if((pMat = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((pMat->U = (RadMHDS***)calloc_3d_array(Nz,Ny, Nx+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((pMat->RHS = (Real****)calloc_4d_array(Nz,Ny, Nx+2*Matghost,2,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


	if((INIerror=(Real****)calloc_4d_array(Nz,Ny,Nx+2*Matghost,2,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_1D]: malloc return a NULL pointer\n");


	/* now set the parameters */

	pMat->dx1 = pG->dx1;
	pMat->dx2 = pG->dx2;
	pMat->dx3 = pG->dx3;
	pMat->time = pG->time;
	/* dt in pG is not set at this time */
	/* pMat->dt = pG->dt;
	*/
	pMat->Lx = pD->RootMaxX[0] - pD->RootMinX[0];
	pMat->Ly = pD->RootMaxX[1] - pD->RootMinX[1];
	pMat->Lz = pD->RootMaxX[2] - pD->RootMinX[2];

	pMat->is = Matghost;
	pMat->ie = Nx + Matghost - 1;
	pMat->js = 0;
	pMat->je = 0;
	pMat->ks = 0;
	pMat->ke = 0;
	pMat->Nx[0] = Nx;
	pMat->Nx[1] = Ny;
	pMat->Nx[2] = Nz;
	pMat->RootNx[0] = Nx;
	pMat->RootNx[1] = Ny;
	pMat->RootNx[2] = Nz;
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
	pMat->BCFlag_ix1 = pM->BCFlag_ix1;
	pMat->BCFlag_ox1 = pM->BCFlag_ox1;
	pMat->BCFlag_ix2 = pM->BCFlag_ix2;
	pMat->BCFlag_ox2 = pM->BCFlag_ox2;
	pMat->BCFlag_ix3 = pM->BCFlag_ix3;
	pMat->BCFlag_ox3 = pM->BCFlag_ox3;
	
		
	/* To decide whether subtract background solution at top level or not */
	/* Default choice is not */
	pMat->bgflag = 0;

}


void BackEuler_destruct_1d()
{

	/* Free pMat and pMatnew */
	if(pMat->U != NULL)
		free_3d_array(pMat->U);

	if(pMat->RHS != NULL)
		free_4d_array(pMat->RHS);

	if(pMat != NULL)
		free(pMat);

	if(INIerror != NULL)
		free_4d_array(INIerror);
}



/*==========================================*/
/* Check the residual after one cycle */
/* This is the stop criterian */
Real CheckResidual(MatrixS *pMat, GridS *pG)
{
	Real Residual = 0.0;
	Real Norm = 0.0;


	int i;
	int is, ie, js, ks;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	ks = pMat->ks;

#ifdef MPI_PARALLEL
	int ierr;
	double tot_residual = 0.0, tot_norm = 0.0;
#endif


	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;

	Real dt = pMat->dt;
	int diffghost = nghost - Matghost;


	/* To store the coefficient */
	Real theta[6];
	Real phi[6];

	

	/* Temporary variables to setup the matrix */
	Real velocity_x, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1;


			for(i=is; i<=ie; i++){


			velocity_x = pMat->U[ks][js][i].V1;

			T4 = pMat->U[ks][js][i].T4;		
				
			/* Assuming the velocity is already the original velocity in case of FARGO */			
				
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
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x;
*/
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
		
		
			
						
			if(pMat->bgflag){
				Norm += fabs(pG->U[ks][js][i+diffghost].Er + pG->Tguess[ks][js][i+diffghost]);
			}
			else{
				Norm += fabs(pMat->RHS[ks][js][i][0]);
			}
			Residual += pMat->RHS[ks][js][i][0];

			Residual -= theta[0] * pMat->U[ks][js][i-1].Er;
			Residual -= theta[1] * pMat->U[ks][js][i-1].Fr1;
			Residual -= theta[2] * pMat->U[ks][js][i].Er;
			Residual -= theta[3] * pMat->U[ks][js][i].Fr1;
			Residual -= theta[4] * pMat->U[ks][js][i+1].Er;
			Residual -= theta[5] * pMat->U[ks][js][i+1].Fr1;



			if(pMat->bgflag){
				Norm += fabs(pG->U[ks][js][i+diffghost].Fr1 + dt * Sigma_aP * T4 * velocity_x);
			}
			else{
				Norm += fabs(pMat->RHS[ks][js][i][1]);
			}
			Residual += pMat->RHS[ks][js][i][1];

			Residual -= phi[0] * pMat->U[ks][js][i-1].Er;
			Residual -= phi[1] * pMat->U[ks][js][i-1].Fr1;
			Residual -= phi[2] * pMat->U[ks][js][i].Er;
			Residual -= phi[3] * pMat->U[ks][js][i].Fr1;
			Residual -= phi[4] * pMat->U[ks][js][i+1].Er;
			Residual -= phi[5] * pMat->U[ks][js][i+1].Fr1;




	

	}

		
	/* MPI call to collect error from different CPUs */


#ifdef MPI_PARALLEL
	ierr = MPI_Reduce(&Residual,&tot_residual,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	ierr = MPI_Reduce(&Norm,&tot_norm,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);

	/* If I'm the parent, copy the sum back to the total_error variable */
	Residual /= Norm;

    if(pMat->ID == 0){
	Norm = tot_norm;
	Residual = tot_residual;

	/* Relative residual */
	Residual /= Norm;

    }

    MPI_Bcast(&Residual,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

#else
	Residual /= Norm;
  
#endif	

	return Residual;

}



/* prolongation operator */


void prolongation1D(MatrixS *pMat_coarse, MatrixS *pMat_fine)
{
	/* Add the correction from the coarse grid back to the solution in fine grid */ 
	int i;
	int is, ie, js, ks;
	int num;

	Real *ptr_fine[2];
	Real *ptr_coarse[3];

	is = pMat_coarse->is;
	ie = pMat_coarse->ie;
	js = pMat_coarse->js;
	ks = pMat_coarse->ks;


			for(i=is; i<=ie; i++){
				/* Take the address */
				ptr_fine[0] = &(pMat_fine->U[ks][js][2*i  ].Er);
				ptr_fine[1] = &(pMat_fine->U[ks][js][2*i-1].Er);

				

				ptr_coarse[0] = &(pMat_coarse->U[ks][js ][i ].Er);
				ptr_coarse[1] = &(pMat_coarse->U[ks][js][i+1].Er);
				ptr_coarse[2] = &(pMat_coarse->U[ks][js][i-1].Er);

				

			/* Only needs to update Er, Frx, Fry, Frz. */
			/* Vx, Vy, Vz and T4 do not need to be updated */
				for(num=0; num<4; num++){
				/* We access a continuous array */
				ptr_fine[0][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[0][num] += 0.25*ptr_coarse[1][num];


      				ptr_fine[1][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[1][num] += 0.25*ptr_coarse[2][num];
 				

				}
	}

}


/* Restriction operator */

void Restriction1D(MatrixS *pMat_fine, MatrixS *pMat_coarse)
{
/* We actually should send the residual to the next level */

	int i;
	int is, ie, js, ks;

	is = pMat_fine->is;
	ie = pMat_fine->ie;
	js = pMat_fine->js;
	ks = pMat_fine->ks;


	Real hdtodx1 = 0.5 * pMat_fine->dt/pMat_fine->dx1;


	/* To store the coefficient */
	Real theta[6];
	Real phi[6];


	/* The right hand size is no-longer the original source terms. It is the residual  */
	/* But we cannot destroy the right hand size of the fine grid. We need that for the */
	/* relaxation step when we come back */ 


	/* Temporary variables to setup the matrix */
	Real velocity_x, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1;

	Real ****error;	
	if((error=(Real****)calloc_4d_array(1,1,pMat_fine->Nx[0]+2*Matghost,2,sizeof(Real))) == NULL)
			ath_error("[Restriction3D]: malloc return a NULL pointer\n");		

	Real *ptr_coarse;
	Real *ptr_fine[2];
	int num;	



			for(i=is; i<=ie; i++){


			velocity_x = pMat_fine->U[ks][js][i].V1;


			T4 = pMat_fine->U[ks][js][i].T4;	
				
			Sigma_sF = pMat->U[ks][js][i].Sigma[0];
			Sigma_aF = pMat->U[ks][js][i].Sigma[1];
			Sigma_aP = pMat->U[ks][js][i].Sigma[2];
			Sigma_aE = pMat->U[ks][js][i].Sigma[3];

			Ci0 = (sqrt(pMat_fine->U[ks][js][i].Edd_11) - sqrt(pMat_fine->U[ks][js][i-1].Edd_11)) 
				/ (sqrt(pMat_fine->U[ks][js][i].Edd_11) + sqrt(pMat_fine->U[ks][js][i-1].Edd_11));
			Ci1 =  (sqrt(pMat_fine->U[ks][js][i+1].Edd_11) - sqrt(pMat_fine->U[ks][js][i].Edd_11)) 
				/ (sqrt(pMat_fine->U[ks][js][i+1].Edd_11) + sqrt(pMat_fine->U[ks][js][i].Edd_11));


			theta[0] = -Crat * hdtodx1 * (1.0 + Ci0) * sqrt(pMat->U[ks][js][i-1].Edd_11);
			theta[1] = -Crat * hdtodx1 * (1.0 + Ci0);
			theta[2] = 1.0 + Crat * hdtodx1 * (2.0 + Ci1 - Ci0) * sqrt(pMat->U[ks][js][i].Edd_11);
/* 
				+ Crat * pMat->dt * Sigma_aE 
				+ pMat->dt * (Sigma_aF - Sigma_sF) * (1.0 + pMat->U[ks][js][i].Edd_11) * velocity_x * velocity_x / Crat;
*/
			theta[3] = Crat * hdtodx1 * (Ci0 + Ci1);
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x;
*/
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
		
		



			error[ks][js][i][0] += pMat_fine->RHS[ks][js][i][0];

			error[ks][js][i][0] -= theta[0] * pMat_fine->U[ks][js][i-1].Er;
			error[ks][js][i][0] -= theta[1] * pMat_fine->U[ks][js][i-1].Fr1;
			error[ks][js][i][0] -= theta[2] * pMat_fine->U[ks][js][i].Er;
			error[ks][js][i][0] -= theta[3] * pMat_fine->U[ks][js][i].Fr1;
			error[ks][js][i][0] -= theta[4] * pMat_fine->U[ks][js][i+1].Er;
			error[ks][js][i][0] -= theta[5] * pMat_fine->U[ks][js][i+1].Fr1;




			error[ks][js][i][1] += pMat_fine->RHS[ks][js][i][1];

			error[ks][js][i][1] -= phi[0] * pMat_fine->U[ks][js][i-1].Er;
			error[ks][js][i][1] -= phi[1] * pMat_fine->U[ks][js][i-1].Fr1;
			error[ks][js][i][1] -= phi[2] * pMat_fine->U[ks][js][i].Er;
			error[ks][js][i][1] -= phi[3] * pMat_fine->U[ks][js][i].Fr1;
			error[ks][js][i][1] -= phi[4] * pMat_fine->U[ks][js][i+1].Er;
			error[ks][js][i][1] -= phi[5] * pMat_fine->U[ks][js][i+1].Fr1;



	}

		/* error = b - Ax */
		/* Now restrict to the coarse grid */
	

			for(i=pMat_coarse->is; i<=pMat_coarse->ie; i++){
				ptr_coarse  = &(pMat_coarse->U[ks][js][i].Er);
				ptr_fine[0] = &(pMat_fine->U[ks][js][2*i ].Er);
				ptr_fine[1] = &(pMat_fine->U[ks][js][2*i-1].Er);


				for(num=0; num<14+NOPACITY; num++){

					ptr_coarse[num] =  (ptr_fine[0][num] + ptr_fine[1][num]) / 2.0;

				}
				/* The initial guess is taken to be zero */
				for(num=0; num<2; num++)
					ptr_coarse[num] = 0.0;
			/*
			pMat_coarse->U[ks][js][i] = (pMat_fine->U[2*k ][2*js ][2*i ]  + pMat_fine->U[2*k ][2*js ][2*i-1]
						+ pMat_fine->U[2*k ][2*js-1][2*i ] + pMat_fine->U[2*k ][2*js-1][2*i-1]
						+ pMat_fine->U[2*k-1][2*js ][2*i ] + pMat_fine->U[2*k-1][2*js ][2*i-1]
						+ pMat_fine->U[2*k-1][2*js-1][2*i ]+ pMat_fine->U[2*k-1][2*js-1][2*i-1]) / 8.0;

			*/
			pMat_coarse->RHS[ks][js][i][0] =  (error[ks][js ][2*i ][0] + error[ks][js][2*i-1][0]) / 2.0;	

			pMat_coarse->RHS[ks][js][i][1] =  (error[ks][js][2*i ][1] + error[ks][js][2*i-1][1]) / 2.0;	


	}

		/* Free the temporary array */
		free_4d_array(error);

	return;
}


void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat)
{

	pMat_coarse->dx1 = 2.0 * pMat->dx1;
	pMat_coarse->dx2 = pMat->dx2;
	pMat_coarse->dx3 = pMat->dx3;
	pMat_coarse->time = pMat->time;
	pMat_coarse->Lx = pMat->Lx;
	pMat_coarse->Ly = pMat->Ly;
	pMat_coarse->Lz = pMat->Lz;
	/* Grid numbers decrease by a factor of 2 */
	pMat_coarse->Nx[0] = pMat->Nx[0] / 2;
	pMat_coarse->Nx[1] = pMat->Nx[1];
	pMat_coarse->Nx[2] = pMat->Nx[2];
	/* top level grid number doesn't change */
	pMat_coarse->RootNx[0] = pMat->RootNx[0];
	pMat_coarse->RootNx[1] = pMat->RootNx[1];
	pMat_coarse->RootNx[2] = pMat->RootNx[2];
	pMat_coarse->NGrid[0] = pMat->NGrid[0];
	pMat_coarse->NGrid[1] = pMat->NGrid[1];
	pMat_coarse->NGrid[2] = pMat->NGrid[2];
	pMat_coarse->dt = pMat->dt;
	pMat_coarse->is = Matghost;
	pMat_coarse->ie = pMat_coarse->Nx[0] + Matghost - 1;
	pMat_coarse->js = 0;
	pMat_coarse->je = 0;
	pMat_coarse->ks = 0;
	pMat_coarse->ke = 0;
		
	pMat_coarse->ID = pMat->ID;
	pMat_coarse->my_iproc = pMat->my_iproc;
	pMat_coarse->my_jproc = pMat->my_jproc; 
	pMat_coarse->my_kproc = pMat->my_kproc;

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

	if((pMat_coarse->U = (RadMHDS***)calloc_3d_array(1,pMat_coarse->Nx[1]+2*Matghost, 	pMat_coarse->Nx[0]+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	/* Allocate memory for the right hand size */
	if((pMat_coarse->RHS = (Real****)calloc_4d_array(1,pMat_coarse->Nx[1]+2*Matghost, 	pMat_coarse->Nx[0]+2*Matghost,2,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


	return;
}

void RHSResidual1D(MatrixS *pMat, Real ****newRHS)
{

	
	int i;
	int is, ie, js, ks;

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	ks = pMat->ks;

	/* To store the coefficient */
	Real theta[6];
	Real phi[6];

	Real hdtodx1 = 0.5 * pMat->dt/pMat->dx1;


	/* Temporary variables to setup the matrix */
	Real velocity_x, T4;
	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	Real Ci0, Ci1;


		for(i=is; i<=ie; i++){


			velocity_x = pMat->U[ks][js][i].V1;

			T4 = pMat->U[ks][js][i].T4;		
				
			/* Assuming the velocity is already the original velocity in case of FARGO */			
				
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
/*	- pMat->dt * (Sigma_aF - Sigma_sF) * velocity_x;
*/
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
		
		
			
						

			
			newRHS[ks][js][i][0] = pMat->RHS[ks][js][i][0];

			newRHS[ks][js][i][0] -= theta[0] * pMat->U[ks][js][i-1].Er;
			newRHS[ks][js][i][0] -= theta[1] * pMat->U[ks][js][i-1].Fr1;
			newRHS[ks][js][i][0] -= theta[2] * pMat->U[ks][js][i].Er;
			newRHS[ks][js][i][0] -= theta[3] * pMat->U[ks][js][i].Fr1;
			newRHS[ks][js][i][0] -= theta[4] * pMat->U[ks][js][i+1].Er;
			newRHS[ks][js][i][0] -= theta[5] * pMat->U[ks][js][i+1].Fr1;




			
			newRHS[ks][js][i][1] = pMat->RHS[ks][js][i][1];

			newRHS[ks][js][i][1] -= phi[0] * pMat->U[ks][js][i-1].Er;
			newRHS[ks][js][i][1] -= phi[1] * pMat->U[ks][js][i-1].Fr1;
			newRHS[ks][js][i][1] -= phi[2] * pMat->U[ks][js][i].Er;
			newRHS[ks][js][i][1] -= phi[3] * pMat->U[ks][js][i].Fr1;
			newRHS[ks][js][i][1] -= phi[4] * pMat->U[ks][js][i+1].Er;
			newRHS[ks][js][i][1] -= phi[5] * pMat->U[ks][js][i+1].Fr1;

	}



	return;

}




#endif /* radMHD_INTEGRATOR */


#endif /* MATRIX_MULTIGRID */

