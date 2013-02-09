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



/* memory to send data */

static double **send_bufP= NULL; 
static double **send_bufRC=NULL; 
static double **recv_bufP= NULL;
#ifdef MPI_PARALLEL
static double **recv_bufRC=NULL;
static MPI_Request **recv_rq=NULL;
static MPI_Request  **send_rq=NULL;
#endif
static int maxND, *start_addrP;

static RadMHDS ***GZ;



/* The initial L2 norm of the right hand side   *
 * The convergence criterian is
 * ||r_n||_2 / ||r_0||_2  < Small number 	*
 */ 

static int Nlim = 4; /* the lim size of the coarsest grid in each CPU*/
static int Vcycle = 15; /* The limit of how many Wcycle we allow */ 

static int Matrixflag = 1; /* The matrix flag is used to choose Gauss_Seidel method or Jacobi method, or Bicgsafe */
			  /* 1 is GS method while 0 is Jacobi method , 2 is Bicgsafe method*/

static int RootLevel; /* Number of levels from below the root domain, log2(N)=Rootlevel */
static int DomainLevels; /* Number of levels above (include) the root domain */
static int TotLevels;
static int RootLevels;

/********Private function ************/
static void Initialize_matrix(MatrixS *pMat, DomainS *pD); /* Function to copy data from grid to the matrix structure, used at the top of each tree */ 

static void RadMHD_multig_3D(MatrixS **Matrix, MeshS *pM);
static void RadMHD_multig_3D_first(MatrixS **Matrix, MeshS *pM);

static void Restriction3D(MatrixS **Matrix, const int Level, const MeshS *pM, const int flag);

static void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat);


static void RHSResidual3D(MatrixS *pMat, Real ****newRHS, Real *error);

static void Calculate_Coef(MatrixS *pMat);


/* New restriction and prolongation function adopted from SMR */

static void Prolongation3D(MatrixS **Matrix, const int Level, const MeshS *pM);




static void ProU(const RadMHDS Uim1, const RadMHDS Ui, const RadMHDS Uip1, 
	  const RadMHDS Ujm1, const RadMHDS Ujp1,
	  const RadMHDS Ukm1, const RadMHDS Ukp1, RadMHDS PCon[2][2][2]);

#ifndef FIRST_ORDER
/* left: vl, right: vr, center: vc. Calculate the TVD slope */
static Real mcd_slope(const Real vl, const Real vc, const Real vr);
#endif /* FIRST_ORDER */




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
	int Wcycle;

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
	while((error > TOL) && (Wcycle < Vcycle)){
		
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
			
			RHSResidual3D(&(Matrix[TotLevels-1][nd]),Matrix[TotLevels-1][nd].RHS,&(Matrix[TotLevels-1][nd].RHSnorm));
			
			pG = pM->Domain[TotLevels-1-RootLevels][nd].Grid;
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
		
		}
				
		
		for(nl=TotLevels-1; nl>RootLevels; nl--){
			/* restrict the solution from top level downwards */
			/* RHS is updated in the Restriction step */
			Restriction3D(Matrix, nl, pM, 0);

			for(nd=0; nd<DomainNos[nl]; nd++){				

				norm = Matrix[nl][nd].RHSnorm0;
				if(norm > TINY_NUMBER)
					error0 = Matrix[nl][nd].RHSnorm/norm;
				else
					error0 = Matrix[nl][nd].RHSnorm;

				error += error0;	
			}
			
			
		} 
		/* Now include error from root domain */
		/* Residual for the RootLevel is already updated for te RootLevel in the Restriction step */
		for(nd=0; nd<DomainNos[RootLevels]; nd++){
			norm = Matrix[RootLevel][nd].RHSnorm0;
			if(norm > TINY_NUMBER)
				error0 = Matrix[RootLevel][nd].RHSnorm/norm;
			else
				error0 = Matrix[RootLevel][nd].RHSnorm;

			error += error0;	
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
			if(Matrix[nl+RootLevels][nd].CPUflag){

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

	int nd, nl, TotLevels;

	TotLevels = DomainLevels + RootLevels;

	

	/* Do not do relaxation when going down for restriction */
	/* Do relaxation when going up after prolongation */
	/* [TotLevels-1] is the top level */
	
	/* First, for the levels above root, which has corresponding grid */
	for(nl=TotLevels-1; nl>=RootLevels; nl--){
		/*****************************************/

		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);

			pD = &(pM->Domain[nl-RootLevels][nd]);
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
	for(nl=TotLevels-1; nl>RootLevels; nl--){
		for(nd=0; nd<DomainNos[nl]; nd++){
			/* restrict the overlap region */
			if(Matrix[nl][nd].CPUflag){
				pMat = &(Matrix[nl][nd]);
				Restriction3D(Matrix, nl, pM, 0);


			}

		}/* Finish domain at level nl */
	}/* Finish looping nl */

	/* Second, for the levels below the root level */
	/* There is only one Domain per Level below the root level */
	/* For these levels, we only need to calculate the matrix coefficient and restrict the matrix parameters */
	/* If this CPU works in this grid at the root domain, it will also work in the levels below that */

	if(Matrix[RootLevels][0].CPUflag){
		
		for(nl=RootLevels-1; nl>=0; nl--){
			/* First do restriction of the matrix coefficient for the whole level */
			Restriction3D(Matrix, nl+1, pM, 1);			
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

		}
		/* prolongate the whole level after relaxation */
		/* We only need to wait for the whole domain to be finished before prolongation */
		Prolongation3D(Matrix, nl,pM);
	}


	/* Do relaxation for the Top level */
		nl = TotLevels-1;
		for(nd=0; nd<DomainNos[nl]; nd++){
			pMat = &(Matrix[nl][nd]);
			
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

	int nd, nl, TotLevels;

	TotLevels = DomainLevels + RootLevels;

	

	/* Restriction to the bottom staring from the Root Domain */
	/* When we calculate the Error, we already restrict to the root domain */
	for(nl=RootLevels-1; nl>=0; nl--){
		
		Restriction3D(Matrix, nl+1, pM, 0);


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
		Prolongation3D(Matrix,nl,pM);
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



/* prolongation operator between Level and Level + 1 */
/* In the prolongation step, for levels above root, */



void Prolongation3D(MatrixS **Matrix, const int Level, const MeshS *pM)
{

	
  	int i, j, k, n, m, l, is, ie, js, je, ks, ke, ii, jj, kk, kfine, jfine, ifine;
	int ics, ice, jcs, jce, kcs, kce;
	int nd, nZeroP, npg, ncg, ngz1, ngz2, ngz3, igzs, igze, jgzs, jgze, kgzs, kgze, mend, nend;
	int nDim;
	int ips, ipe, jps, jpe, kps, kpe;	

	double *pRcv, *pSnd;
  	
	GridOvrlpS *pCO, *pPO;
	MatrixS *pMat, *pMatC;	
	RadMHDS Ptemp[2][2][2];	
	GridS *pG;
#ifdef MPI_PARALLEL
  	int mAddress, ierr, mIndex, mCount;
#endif

	/* number of dimensions in Grid. */
	/* First, determine the dimensionality */
  	nDim=1;
  	for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

	/*===========================================================================*/
	/* We only have two levels to deal with at each time */
	/* First, we need to judge whether the Level if below the RootLevel or not */ 
	/* Below the root level, we use the original way */
	if(Level < RootLevels){
		pMat  = &(Matrix[Level][0]); /* Matrix in the coarse Level */
		pMatC = &(Matrix[Level+1][0]);/* Matrix in the fine Level */

		is = pMat->is;
		ie = pMat->ie;
		js = pMat->js;
		je = pMat->je;
		ks = pMat->ks;
		ke = pMat->ke;

		for(k=ks; k<=ke; k++)
			for(j=js; j<=je; j++)
				for(i=is; i<=ie; i++){

					ProU(pMat->U[k][j][i-1],pMat->U[k][j][i],pMat->U[k][j][i+1],pMat->U[k][j-1][i],pMat->U[k][j+1][i],
						pMat->U[k-1][j][i],pMat->U[k+1][j][i],Ptemp);


			/* Now copy the data to the fine grid */
					ii = 2*(i-is) + pMatC->is;
					jj = 2*(j-js) + pMatC->js;
					kk = 2*(k-ks) + pMatC->ks;

				/* The coarse grid calculates the residual, but solution is already added to grid solution */
					for(kfine=0;kfine<2; kfine++)
						for(jfine=0; jfine<2; jfine++)
							for(ifine=0; ifine<2; ifine++){
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Er  = Ptemp[kfine][jfine][ifine].Er;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr1 = Ptemp[kfine][jfine][ifine].Fr1;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr2 = Ptemp[kfine][jfine][ifine].Fr2;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr3 = Ptemp[kfine][jfine][ifine].Fr3;
					}/* Finish copy solution back to child grid */				

		}/* Finish loop the whole grid i, j, k */

	}/* End if the Level is below the RootLevel */
	else{
		/* First, post non-blocking receives at Level + 1 for data from parent GridS */


#ifdef MPI_PARALLEL

		
    		for (nd=0; nd<(DomainNos[Level+1]); nd++){
      			if (pM->Domain[Level+1-RootLevels][nd].Grid != NULL) {
        			
				pG=pM->Domain[Level+1-RootLevels][nd].Grid;
        			nZeroP = 0;
        			mAddress = 0;        			
        			if (pG->NmyPGrid > 0) mAddress = pG->PGrid[0].Rad_nWordsP;

        		for (npg=(pG->NmyPGrid); npg<(pG->NPGrid); npg++){

			/* Skip if no prolongation needed for this child (only flux correction) */ 
			/* This is for the case when grids only touch */
          			if (pG->PGrid[npg].Rad_nWordsP == 0) { 
            				nZeroP += 1;
          			} else {

           				mIndex = npg - pG->NmyPGrid - nZeroP;
            				ierr = MPI_Irecv(&(recv_bufP[nd][mAddress]),
              				pG->PGrid[npg].Rad_nWordsP, MPI_DOUBLE, pG->PGrid[npg].ID,
              				pG->PGrid[npg].DomN, pM->Domain[Level+1-RootLevels][nd].Comm_Parent,
              					&(recv_rq[nd][mIndex]));
            				mAddress += pG->PGrid[npg].Rad_nWordsP;
          			}

        		}/* End loop all the parent grids */
      			}/* End if Grid != NULL */
    		}/* Finish looping all domains at Level +1 */
  
#endif /* MPI_PARALLEL */



		/*======Step 1, send the data, including the ghoze zones for Child grids==========*/
		/* Unlik the SMR in normal MHD part, we need to send the whole grid data */
		/* Now we works for Level */
		/* We need to send the data including ghost and active zones */
		for(nd=0; nd<DomainNos[Level]; nd++){
			if(Matrix[Level][nd].CPUflag){
				pG = pM->Domain[Level-RootLevels][nd].Grid;
				pMat = &(Matrix[Level][nd]);

				for(i=0; i<maxND; i++) start_addrP[i] = 0;
				nZeroP = 0;
				
				for(ncg=0; ncg<(pG->NCGrid); ncg++){
					if(pG->CGrid[ncg].Rad_nWordsP == 0){
						nZeroP += 1;
					}/* skip the grid that does need prolongation */
					else{
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);	/* ptr to child Grid overlap */

					/* index send_buf with DomN of child, since could be multiple child Domains on */
 					/* same processor.  Start address must be different for each DomN */


						pSnd = (double*)&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]);
						
						
						ics = pCO->ijks[0] - nghost + Matghost - (Matghost/2) - 1;
						ice = pCO->ijke[0] - nghost + Matghost + (Matghost/2) + 1;
						
						if(nDim > 1){
							jcs = pCO->ijks[1] - nghost + Matghost - (Matghost/2) - 1;
							jce = pCO->ijke[1] - nghost + Matghost + (Matghost/2) + 1;
						}
						else{
							jcs = pCO->ijks[1];
							jce = pCO->ijke[1];
						}

						if(nDim > 2){
							kcs = pCO->ijks[2] - nghost + Matghost - (Matghost/2) - 1;
							kce = pCO->ijks[2] - nghost + Matghost + (Matghost/2) + 1;
						}
						else{
							kcs = pCO->ijks[2];
							kce = pCO->ijks[2];
						}

							

						for(k=kcs; k<=kce; k++){
						for(j=jcs; j<=jce; j++){
						for(i=ics; i<=ice; i++){
							*(pSnd++) = pMat->U[k][j][i].Er;
							*(pSnd++) = pMat->U[k][j][i].Fr1;
							*(pSnd++) = pMat->U[k][j][i].Fr2;
							*(pSnd++) = pMat->U[k][j][i].Fr3;

						}/* end i */
						}/* end j */
						}/* end k */
						
					}/* End for the grids that needs prolongation */
				}/* End loop over the child grid */

			/* Step 1b: non-blocking send of data  to Child, using Domain number as tag */
#ifdef MPI_PARALLEL
      				if (ncg >= pG->NmyCGrid) {
        				mIndex = ncg - pG->NmyCGrid - nZeroP;
        				ierr = MPI_Isend(&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]),
          				pG->CGrid[ncg].nWordsP, MPI_DOUBLE, pG->CGrid[ncg].ID, nd,
          				pM->Domain[Level-RootLevels][nd].Comm_Children, &(send_rq[nd][mIndex]));
      				}
#endif /* MPI_PARALLEL */

      				start_addrP[pCO->DomN] += pG->CGrid[ncg].Rad_nWordsP;

			}/* End if CPUflag */
		}/* Finish looping all the domains at Level */



		/*==============================================================*/
		/* Because data is sent from Level, we need to clear the send_bufP for Level */
		/* Unlike the normal smr function for MHD, we only have two levels here and we do not need to loop over level */
		/* So for grids on the same CPU, we need set recv pointer to send buffer */

  	/******************************************************************************/
	/* This step is skipped as we just set Receive pointer to the send buffer if grids are on the same CPU 
	
		for (nd=0; nd<(pM->DomainsPerLevel[Level-RootLevels]); nd++){
    			if (pM->Domain[Level-RootLevels][nd].Grid != NULL) { 
     				 pG=pM->Domain[Level-RootLevels][nd].Grid; 
	

      				for (ncg=0; ncg<(pG->NmyCGrid); ncg++){
        				pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);   

        				for (i=0; i<pCO->Rad_nWordsP; i++) {
          					recv_bufP[pCO->DomN][i]=send_bufP[pCO->DomN][i];
        				}
      				}
    			}
  		}

 */

#ifdef MPI_PARALLEL
		/* For MPI jobs, wait for all non-blocking sends above to finish in order to continue to Level+1  */
		/* Ortherwise, the CPU may start to working while sent is not complete */

  		for (nd=0; nd<(pM->DomainsPerLevel[Level-RootLevels]); nd++){
    			if (pM->Domain[Level-RootLevels][nd].Grid != NULL) {
     				pG=pM->Domain[Level-RootLevels][nd].Grid;

      				nZeroP = 0;
      				for (i=0; i < pG->NCGrid; i++) if (pG->CGrid[i].Rad_nWordsP == 0) nZeroP++;

      				if (pG->NCGrid > pG->NmyCGrid) {
        				mCount = pG->NCGrid - pG->NmyCGrid - nZeroP;
        				ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      				}
    			}
  		}
#endif /* MPI_PARALLEL */



		/*=====================================================================*/
		/* Get solution from parent GridS and prolongation solution to ghost zones */
		/* Now we go back to Level + 1 */
		for(nd=0; nd<DomainNos[Level+1]; nd++){
			if(pM->Domain[Level+1-RootLevels][nd].Grid != NULL){
				pG = pM->Domain[Level+1-RootLevels][nd].Grid;
				pMat = &(Matrix[Level+1][nd]);

				/* Loop over number of parent grids with non-zero-size prolongation data */
				nZeroP = 0;
    				for (i=0; i < pG->NPGrid; i++) if (pG->PGrid[i].Rad_nWordsP == 0) nZeroP++;

				for (npg=0; npg<(pG->NPGrid - nZeroP); npg++){

				/* If parent Grid is on this processor, data is at start of recv buffer */

      					if (npg < pG->NmyPGrid) {
        					pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
						/* For grids on the same CPU, set the pointer to send buffer */
        					/* pRcv = (double*)&(recv_bufP[nd][0]); */
						pRcv = (double*)&(send_bufP[nd][0]);
	
      					} else {

#ifdef MPI_PARALLEL
					/* Check non-blocking receives posted above for data in ghost zone from parent
 						* Grids, sent in Step 1.  Accept messages in any order. */

        					mCount = pG->NPGrid - pG->NmyPGrid - nZeroP;
        					ierr = MPI_Waitany(mCount,recv_rq[nd],&mIndex,MPI_STATUS_IGNORE);
        					if(mIndex == MPI_UNDEFINED){
          						ath_error("[Prolong]: Invalid request index nl=%i nd=%i\n",Level+1,nd);
        					}

						/* mIndex returns the number that is completed */

					/* Recv buffer is addressed from PGrid[0].Rad_nWordsP for first MPI message
 					* if NmyPGrid>0.  Also must remove zero size messages from index. */

						/* re-build the value of mIndex to include NmyPGrid and grids that do not need prolongation */

        					mIndex += pG->NmyPGrid;
        					for (i=pG->NmyPGrid; i <= mIndex; i++) 
          						if (pG->PGrid[i].Rad_nWordsP == 0) mIndex++;

        					mAddress = 0;
        					for (i=0; i<mIndex; i++) mAddress += pG->PGrid[i].Rad_nWordsP;
        					pPO = (GridOvrlpS*)&(pG->PGrid[mIndex]); 
        					pRcv = (double*)&(recv_bufP[nd][mAddress]);
#else
				/* If not MPI_PARALLEL, and parent Grid not on this processor, then error */

        					ath_error("[Prolong]: no Parent Grid on Domain[%d][%d]\n",Level+1,nd);
#endif /* MPI_PARALLEL */
      					}/* End if npg > pG->NmyPGrid */

					
					/* Loop over 6 boundaries, set ghost zones */
					/* The difference between Matrix solver and normal MHD is that */
					/* We also need to prolongate data from each cell to the child grids */

					/* Get coordinates ON THIS GRID of ghost zones that overlap parent Grid */

					ngz1 = (pPO->ijke[0] - pPO->ijks[0] + 1)/2 + Matghost + 2;

					if(nDim > 1)
						ngz2 = (pPO->ijke[1] - pPO->ijks[1] + 1)/2 + Matghost + 2;
					else
						ngz2 = (pPO->ijke[1] - pPO->ijks[1] + 1)/2;


					if(nDim > 2)
						ngz3 = (pPO->ijke[2] - pPO->ijks[2] + 1)/2 + Matghost + 2;
					else
						ngz3 = (pPO->ijke[2] - pPO->ijks[2] + 1)/2;
	

					igzs = 0;
					igze = ngz1 - 1;
				
					if (pMat->Nx[1] > 1) {
						jgzs = 0;
						jgze = ngz2 - 1;
						mend = 1;
					} else {
						ngz2 = 1;
						jgzs = 1;
						jgze = 1;
						mend = 0;
					}
					if (pMat->Nx[2] > 1) {
						kgzs = 0;
						kgze = ngz3 - 1;
						nend = 1;
					} else {
						ngz3 = 1;
						kgzs = 1;
						kgze = 1;
						nend = 0;
					}

					/* Load the data */
					
					for(k=kgzs; k<=kgze; k++){
					for(j=jgzs; j<=jgze; j++){
					for(i=igzs; i<=igze; i++){
						GZ[k][j][i].Er = *(pRcv++);
						GZ[k][j][i].Fr1 = *(pRcv++);
						GZ[k][j][i].Fr2 = *(pRcv++);
						GZ[k][j][i].Fr3 = *(pRcv++);

					}/* End i */
					}/* End j */
					}/* End k */

					/* Set boundary conditions for GZ in 1D and 2D cases */
					/* This is needed for the prolongation below */
					 if (nDim == 1) {
            					for (i=igzs; i<=igze; i++) {
							GZ[1][0][i] = GZ[1][1][i];
              						GZ[1][2][i] = GZ[1][1][i];
              						GZ[0][1][i] = GZ[1][1][i];
              						GZ[2][1][i] = GZ[1][1][i];
							
            					}/* End for i*/
          				}/* End if nDim = 1 */

          				if (nDim == 2) {
            					for (j=jgzs; j<=jgze; j++) {
            					for (i=igzs; i<=igze; i++) {
							GZ[0][j][i] = GZ[1][j][i];
              						GZ[2][j][i] = GZ[1][j][i];								
            					}/* End for i */
						}/* End for j */
					}/* End if nDim = 2*/


					/* Now prolongate the ghost zones in array GZ */

					ips = pPO->ijks[0] - nghost;	/* This actually is - nghost + Matghost - Matghost */
					ipe = pPO->ijke[0] - nghost + Matghost + Matghost;
						
					if(pG->Nx[1] > 1){
						jps = pPO->ijks[1] - nghost;
						jpe = pPO->ijke[1] - nghost + Matghost + Matghost;
					}else{
						jps = pPO->ijks[1]; /* The value will be zero in this case */
						jpe = pPO->ijke[1];
					}

					if(pG->Nx[2] > 1){
						kps = pPO->ijks[2] - nghost;
						kpe = pPO->ijke[2] - nghost + Matghost + Matghost;
					}else{
						kps = pPO->ijks[2]; /* The value will be zero in this case */
						kpe = pPO->ijke[2];
					}

						/* Prolongate the GZ array data to Ptemp temporarily and then copy the data to pMat */
						/* i, j, k are for parent grids while kk, jj, ii are for child grids */
						

					for (k=kps, kk=1; k<=kpe; k+=2, kk++) {
					for (j=jps, jj=1; j<=jpe; j+=2, jj++) {
       					for (i=ips, ii=1; i<=ipe; i+=2, ii++) {
						ProU(GZ[kk][jj][ii-1], GZ[kk][jj][ii], GZ[kk][jj][ii+1], GZ[kk][jj-1][ii], 
								GZ[kk][jj+1][ii], GZ[kk-1][jj][ii], GZ[kk+1][jj][ii], Ptemp);

						/* Now set the solution */
						for(n=0; n<=nend; n++){
						for(m=0; m<=mend; m++){
						for(l=0; l<=1; l++){
							pMat->U[k+n][j+m][i+l].Er  = Ptemp[n][m][l].Er;
							pMat->U[k+n][j+m][i+l].Fr1 = Ptemp[n][m][l].Fr1;
							pMat->U[k+n][j+m][i+l].Fr2 = Ptemp[n][m][l].Fr2;
							pMat->U[k+n][j+m][i+l].Fr3 = Ptemp[n][m][l].Fr3;
						}/* End l */
						}/* End m */
						}/* End n */

					}/* End i */
					}/* End j */
					}/* End k */

				}/* End loop npg grid */



			}/* End if [Level+1][nd].Grid is not NULL */
		}/* Finish loop over domains at Level + 1*/

		/*===================================================================*/

	}/* End if the Level is above the RootLevel */




} /* End prolongation function */





/* This function only restricts the gas quantities from Level to Level-1 */
/* flag == 1: restrict the coefficient only; flag == 0:  restrict solution only */
/* We need pM here because we need to know the information about child and paraent grids */
/* Because the matrix coefficient needs primative variables, restrictions are done */
/* Even for levels above the root domain */


/* This function restricts the matrix coefficient, solution and RHS */
/* Assuming that RHS is already updated with the new solution */
/* We actually also copy the solution back from matrix to the grid */
/* Because this is just after the relaxation step */

void Restriction3D(MatrixS **Matrix, const int Level, const MeshS *pM, const int flag)
{
	int i, j, k, ips, ipe, jps, jpe, kps, kpe, ics, ice, jcs, jce, kcs, kce;
	int nd, npg, ncg, start_addr;
	MatrixS *pMat, *pMatP;	
	GridS *pG;

	double *pRcv, *pSnd;
	GridOvrlpS *pPO, *pCO;

	Real *ptrP;
	Real *ptr[8];
	int num;	
#ifdef MPI_PARALLEL
  	int ierr, mIndex, mCount, mAddress;
#endif
	

/*======== Step 3.: Restrict child solution and send =================*/
/* Loop over all Domains and parent GridS */
/* [Level] is the child level and [Level-1] is the parent level */
/* Separate levels above and below the RootLevel: Below the RootLevels, there is only one domain. */
/* Every child and parent grids are in the same CPU */
/* For levels below RootLevles, we do not need to restrict the solution */
/*==================================================================================*/
 	if((Level <= RootLevels) && (flag)){
		pMat  = &(Matrix[Level][0]);
		pMatP = &(Matrix[Level-1][0]);

		/*****************************************/
		/* Do not calculate the Residual here */	
		/* We wil do this before restriction function is called, if necessary */
		

		for(k=pMatP->ks; k<=pMatP->ke; k++)
			for(j=pMatP->js; j<=pMatP->je; j++)
				for(i=pMatP->is; i<=pMatP->ie; i++){
					if(flag){ /* We will never need to restrict the solution for level below RootLevels */
						ptrP  = &(pMatP->Ugas[k][j][i].rho);
						ptr[0] = &(pMat->Ugas[2*k ][2*j ][2*i ].rho);
						ptr[1] = &(pMat->Ugas[2*k ][2*j ][2*i-1].rho);
						ptr[2] = &(pMat->Ugas[2*k ][2*j-1][2*i ].rho);
						ptr[3] = &(pMat->Ugas[2*k ][2*j-1][2*i-1].rho);
						ptr[4] = &(pMat->Ugas[2*k-1][2*j ][2*i ].rho);
						ptr[5] = &(pMat->Ugas[2*k-1][2*j ][2*i-1].rho);
						ptr[6] = &(pMat->Ugas[2*k-1][2*j-1][2*i ].rho);	
						ptr[7] = &(pMat->Ugas[2*k-1][2*j-1][2*i-1].rho);

						for(num=0; num<11+NOPACITY; num++){

							ptrP[num] =  (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
						}
					}

					/* Now for the RHS */
					ptr[0] = &(pMat->RHS[2*k ][2*j ][2*i ][0]);
					ptr[1] = &(pMat->RHS[2*k ][2*j ][2*i-1][0]);
					ptr[2] = &(pMat->RHS[2*k ][2*j-1][2*i ][0]);
					ptr[3] = &(pMat->RHS[2*k ][2*j-1][2*i-1][0]);
					ptr[4] = &(pMat->RHS[2*k-1][2*j ][2*i ][0]);
					ptr[5] = &(pMat->RHS[2*k-1][2*j ][2*i-1][0]);
					ptr[6] = &(pMat->RHS[2*k-1][2*j-1][2*i ][0]);	
					ptr[7] = &(pMat->RHS[2*k-1][2*j-1][2*i-1][0]);
			
					for(num=0; num<4; num++){
						pMatP->RHS[k][j][i][num] = (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0; 
					}
					
	
						
		}/* End if k, j, i */

		/* Also set ghost zones to be zero */
		
		for(k=pMatP->ks-Matghost; k<=pMatP->ke+Matghost; k++)
			for(j=pMatP->js-Matghost; j<=pMatP->je+Matghost; j++)
				for(i=pMatP->is-Matghost; i<=pMatP->ie+Matghost; i++){

				
					pMatP->U[k][j][i].Er = 0.0;
					pMatP->U[k][j][i].Fr1 = 0.0;
					pMatP->U[k][j][i].Fr2 = 0.0;
					pMatP->U[k][j][i].Fr3 = 0.0;

		}
	
	}/* End if Level <= RootLevels */
/*==================================================================================*/
	else{
		/* We only work for Two levels at a time. When do restriction and receive, we should decide first, */ 
		/* Whether this CPU works in this grid */

		/* For MPI, first post non-block receiver */
		/* This will only post with MPI and for CPU work on this grid */
#ifdef MPI_PARALLEL
/* Post non-blocking receives at level Level-1 for data from child Grids at this
 * level (nl).  This data is sent in Step 3 below. */ 


    		for (nd=0; nd<pM->DomainsPerLevel[Level-RootLevels-1]; nd++){
      			if (pM->Domain[Level-RootLevels-1][nd].Grid != NULL) {
        			pG=pM->Domain[Level-RootLevels-1][nd].Grid;
        			mAddress = 0;
        			for (ncg=(pG->NmyCGrid); ncg<(pG->NCGrid); ncg++){
          				mIndex = ncg - pG->NmyCGrid;
          				ierr = MPI_Irecv(&(recv_bufRC[nd][mAddress]),
            				pG->CGrid[ncg].Rad_nWordsRC, MPI_DOUBLE, pG->CGrid[ncg].ID,
            				pG->CGrid[ncg].DomN, pM->Domain[Level-RootLevels-1][nd].Comm_Children,
            				&(recv_rq[nd][mIndex]));
          				mAddress += pG->CGrid[ncg].Rad_nWordsRC;
        			}/* End all the child grid */

      			}/* End if the grid is null */
    		}/* End loop all domains at Level -1 */
 
#endif /* MPI_PARALLEL */		

		/* First, at Level, restrict and send the data, this will be the first step anyway */
		for(nd=0; nd<DomainNos[Level]; nd++){
			/* If this CPU works for this Domain and this grid */
			if(Matrix[Level][nd].CPUflag){
				pG = pM->Domain[Level-RootLevels][nd].Grid; /* Level is gaurantee to be larger than RootLevels */
				pMat = &(Matrix[Level][nd]);
				start_addr = 0;
				/* Residual is calculated in multigrid main cycle. Here we assume residual is already calculated */

			
				/* There could be multiple parent grids overlap with this child grid */
				for(npg=0; npg<(pG->NPGrid); npg++){
					pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
					
					/* Get coordinates ON THIS fine MATRIX of overlap region of parent MATRIX */
					ips = pPO->ijks[0] - nghost + Matghost;
					ipe = pPO->ijke[0] - nghost + Matghost;
					jps = pPO->ijks[1] - nghost + Matghost;
					jpe = pPO->ijke[1] - nghost + Matghost;
					kps = pPO->ijks[2] - nghost + Matghost;
					kpe = pPO->ijke[2] - nghost + Matghost;

				
					/* Now restrict coefficient, RHS and solution */
					pSnd = (double*)&(send_bufRC[nd][start_addr]);

					for(k=kps; k<=kpe; k+=2){
					for(j=jps; j<=jpe; j+=2){
					for(i=ips; i<=ipe; i+=2){
						if(flag){
						/* Restrict gas quantities to calculate the coefficients */
							ptr[0] = &(pMat->Ugas[k][j][i].rho);
							ptr[1] = &(pMat->Ugas[k][j][i+1].rho);
							ptr[2] = &(pMat->Ugas[k][j+1][i].rho);
							ptr[3] = &(pMat->Ugas[k][j+1][i+1].rho);
							ptr[4] = &(pMat->Ugas[k+1][j ][i].rho);
							ptr[5] = &(pMat->Ugas[k+1][j][i+1].rho);
							ptr[6] = &(pMat->Ugas[k+1][j+1][i].rho);	
							ptr[7] = &(pMat->Ugas[k+1][j+1][i+1].rho);

							for(num=0; num<11+NOPACITY; num++){

								pSnd[num] =  (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
							}
							/* move the pointer  */
							pSnd += (11 + NOPACITY);				
						}/* end if flag = 1, end restricting matrix coefficient */
						/* Now restrict the new RHS hand side */
						
						
						for(num=0; num<4; num++){
							pSnd[num] = (pMat->RHS[k][j][i][num] + pMat->RHS[k][j][i+1][num]
						     		+  pMat->RHS[k][j+1][i][num] + pMat->RHS[k][j+1][i+1][num]
						     		+  pMat->RHS[k+1][j][i][num] + pMat->RHS[k+1][j][i+1][num]
						     		+  pMat->RHS[k+1][j+1][i][num]+ pMat->RHS[k+1][j+1][i+1][num]) / 8.0;	
						} /* End restricting RHS */

						/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
						/* than coarse grid solution */
					
							/* move the pointer */
						pSnd += 4;

						ptr[0] = &(pMat->U[k][j][i].Er);
						ptr[1] = &(pMat->U[k][j][i+1].Er);
						ptr[2] = &(pMat->U[k][j+1][i].Er);
						ptr[3] = &(pMat->U[k][j+1][i+1].Er);
						ptr[4] = &(pMat->U[k+1][j ][i].Er);
						ptr[5] = &(pMat->U[k+1][j][i+1].Er);
						ptr[6] = &(pMat->U[k+1][j+1][i].Er);	
						ptr[7] = &(pMat->U[k+1][j+1][i+1].Er);
						
						for(num=0; num<4; num++){
							pSnd[num] = (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
						 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
						}

					}/* end ips */		
					}/* end jps */
					}/* end kps */
					
#ifdef MPI_PARALLEL
			/* send the data for MPI case */
					 if (npg >= pG->NmyPGrid){
        					mIndex = npg - pG->NmyPGrid;
        					ierr = MPI_Isend(&(send_bufRC[nd][start_addr]), pG->PGrid[npg].Rad_nWordsRC,
          						MPI_DOUBLE, pG->PGrid[npg].ID, nd, pM->Domain[Level][nd].Comm_Parent,
          						&(send_rq[nd][mIndex]));
     					}

#endif /* MPI_PARALLEL */
					/* set the start address prepared for next parent grid */
					start_addr += pG->PGrid[npg].Rad_nWordsRC;

				}/* End npg parent grid */				
			}/* End if Matrix[Level][nd].CPUflag */
		} /* End loop all the domains at this level */

		/* At the end, we need to wait for send from all domains to be finished at Level */

#ifdef MPI_PARALLEL

  		for(nd=0; nd<(DomainNos[Level]); nd++){
    			if ((pM->Domain[Level-RootLevels][nd].Grid) != NULL) {
      				pG = pM->Domain[Level-RootLevels][nd].Grid;

      				if (pG->NPGrid > pG->NmyPGrid) {
        				mCount = pG->NPGrid - pG->NmyPGrid;
        				ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      				}
    			}
  		}/* Finish all Domains at Level */
#endif /* MPI_PARALLEL */			


/*===================================================================*/
/* After send the data, now we need to receive the data for Level -1 */

		/* Get Child solution */
		for(nd=0; nd<DomainNos[Level-1]; nd++){
			if(Matrix[Level-1][nd].CPUflag){
				pG = pM->Domain[Level-1-RootLevels][nd].Grid; /* Level is gaurantee to be larger than RootLevels */
				pMat = &(Matrix[Level-1][nd]);

					
			
				for(ncg=0; ncg<(pG->NCGrid); ncg++){
					if(ncg < pG->NmyCGrid){
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);
						/* send pointer to the beginning of the send buffer, if on the same process */
						pRcv = (double*)&(send_bufRC[pCO->DomN][0]);						

					} /* For the child grid on the same CPU */
					else {
#ifdef MPI_PARALLEL
						mCount = pG->NCGrid - pG->NmyCGrid;
        					ierr = MPI_Waitany(mCount,recv_rq[nd],&mIndex,MPI_STATUS_IGNORE);
        					if(mIndex == MPI_UNDEFINED){
          						ath_error("[RestCorr]: Invalid request index nl=%i nd=%i\n",Level-1,nd);
        					}
      
					/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0 */
        					mAddress = 0;
        					mIndex += pG->NmyCGrid;
        					for (i=pG->NmyCGrid; i<mIndex; i++) mAddress += pG->CGrid[i].Rad_nWordsRC;
        						pCO=(GridOvrlpS*)&(pG->CGrid[mIndex]);
        						pRcv = (double*)&(recv_bufRC[nd][mAddress]);
#else
				/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        					ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",Level-1,nd);
#endif /* MPI_PARALLEL */

					}

					/* shift the position with ghost zone difference */
					/* Now the index corresponds to position in the matrix structure */

					

					ics = pCO->ijks[0] - nghost + Matghost;
      					ice = pCO->ijke[0] - nghost + Matghost;
      					jcs = pCO->ijks[1] - nghost + Matghost;
      					jce = pCO->ijke[1] - nghost + Matghost;
      					kcs = pCO->ijks[2] - nghost + Matghost;
      					kce = pCO->ijke[2] - nghost + Matghost;

					

					/* First, get the solution or matrix coefficient */
					/* The restricted RHS is stored at a temporary first */
					/* The restricted RHS is b - AX, which is the residual */
					/* We need to first calculate b_2 - AX_2 for the whole coarse level, as it is not completely covered by fine level */
					/* Then we replace the RHS with the restricted value */
					
					for(k=kcs; k<=kce; k++){
					for(j=jcs; j<=jce; j++){
					for(i=ics; i<=ice; i++){
						if(flag){
						/* Restrict gas quantities to calculate the coefficients */
							ptrP  = &(pMat->Ugas[k][j][i].rho);							

							for(num=0; num<11+NOPACITY; num++){
								ptrP[num] =  pRcv[num];
							}
							/* move the pointer  */
							pRcv += (11 + NOPACITY);				
						}/* end if flag = 1, end restricting matrix coefficient */
						/* Now restrict the new RHS hand side */
						
						/* Do not replace RHS right now */	
					/*	ptrP = &(pMat->RHS[k][j][i][0]);

						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];	
						}
					*/
						/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
						/* than coarse grid solution */
					
						/* move the pointer */
						pRcv += 4;

						ptrP = &(pMat->U[k][j][i].Er);
						
						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];
						}
						/* move the pointer */
						pRcv += 4;

					}/* end ics */		
					}/* end jcs */
					}/* end kcs */
					/* Now update the solution in the overlap region, update the RHS for this domain */
				}/* Loop over all the child grids */

				/* Now we need to update the RHS and solution for this level */

				RHSResidual3D(pMat, pMat->RHS, &(pMat->RHSnorm));

				/* update RHS for the whole domain */
				/* We do not need ghost zones for RHS */
				for(k=pMat->ks; k<=pMat->ke; k++)
				for(j=pMat->js; j<=pMat->je; j++)
				for(i=pMat->is; i<=pMat->ie; i++){
					
					/* Now update the solution */
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Er  += pMat->U[k][j][i].Er;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr1 += pMat->U[k][j][i].Fr1;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr2 += pMat->U[k][j][i].Fr2;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr3 += pMat->U[k][j][i].Fr3;
				
				}


				/* Now we need to replace the RHS in the overlap region with restricted RHS from parent grids */
				/* The restricted data is already in the recv_buf, so we do not need to wait for MPI again */
 
		/*===================================================================================*/
				for(ncg=0; ncg<(pG->NCGrid); ncg++){
					if(ncg < pG->NmyCGrid){
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);
						/* send pointer to the beginning of the send buffer, if on the same process */
						pRcv = (double*)&(send_bufRC[pCO->DomN][0]);						

					} /* For the child grid on the same CPU */
					else {
#ifdef MPI_PARALLEL
      
					/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0 */
        					mAddress = 0;
        					mIndex += pG->NmyCGrid;
        					for (i=pG->NmyCGrid; i<mIndex; i++) mAddress += pG->CGrid[i].Rad_nWordsRC;
        						pCO=(GridOvrlpS*)&(pG->CGrid[mIndex]);
        						pRcv = (double*)&(recv_bufRC[nd][mAddress]);
#else
				/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        					ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",Level-1,nd);
#endif /* MPI_PARALLEL */

					}

					/* shift the position with ghost zone difference */
					/* Now the index corresponds to position in the matrix structure */
					ics = pCO->ijks[0] - nghost + Matghost;
      					ice = pCO->ijke[0] - nghost + Matghost;
      					jcs = pCO->ijks[1] - nghost + Matghost;
      					jce = pCO->ijke[1] - nghost + Matghost;
      					kcs = pCO->ijks[2] - nghost + Matghost;
      					kce = pCO->ijke[2] - nghost + Matghost;

					/* First, get the solution or matrix coefficient */
					/* The restricted RHS is stored at a temporary first */
					
					for(k=kcs; k<=kce; k++){
					for(j=jcs; j<=jce; j++){
					for(i=ics; i<=ice; i++){
						if(flag){
							pRcv += (11 + NOPACITY);				
						}/* end if flag = 1, end restricting matrix coefficient */
						
						/* Now replace the restrcited RHS */
						
						
						ptrP = &(pMat->RHS[k][j][i][0]);

						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];	
						}
					
						/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
						/* than coarse grid solution */
					
						/* move the pointer */
						pRcv += 8;

					}/* end ics */		
					}/* end jcs */
					}/* end kcs */
					/* Now update the solution in the overlap region, update the RHS for this domain */
				}/* Loop over all the child grids */

		/*==================================================================================*/

			}/* End if [Level-1][nd] CPU flag */
		}/* End loop all domains at Level-1 */


	}/* End for levels above RootLevel */




	return;
}





/*  calculate the residual of right hand side from guess solution */
/* should allocate memory for newRHS before this function is called */
/* We also need MPI communicator for this Domain */

void RHSResidual3D(MatrixS *pMat, Real ****newRHS, Real *error)
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
	int n;

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	/* current level */
	n = pMat->Level;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

				matrix_coef(pMat, NULL, 3, i, j, k, 0.0, &(pMat->Ptheta[k][j][i][0]), &(pMat->Pphi[k][j][i][0]), &(pMat->Ppsi[k][j][i][0]), &(pMat->Pvarphi[k][j][i][0]));
		}

}


/* Prolongation scheme adopted from SMR */
/* Only need to prolongate Er, Fr1,Fr2, Fr3 four variables */
/* from coarse grid to fine grid */
/* PU is used to take the data out */

void ProU(const RadMHDS Uim1, const RadMHDS Ui, const RadMHDS Uip1, 
	  const RadMHDS Ujm1, const RadMHDS Ujp1,
	  const RadMHDS Ukm1, const RadMHDS Ukp1, RadMHDS PCon[2][2][2])
{
  int i,j,k;

  Real dq1,dq2,dq3;


/* First order prolongation -- just copy values */
#ifdef FIRST_ORDER

  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Er  = Ui.Er;
    PCon[k][j][i].Fr1 = Ui.Fr1;
    PCon[k][j][i].Fr2 = Ui.Fr2;
    PCon[k][j][i].Fr3 = Ui.Fr3;
  }}}
/* second order prolongation -- apply limited slope reconstruction */
#else /* SECOND_ORDER or THIRD_ORDER */

/* Er */
  dq1 = mcd_slope(Uim1.Er, Ui.Er, Uip1.Er);
  dq2 = mcd_slope(Ujm1.Er, Ui.Er, Ujp1.Er);
  dq3 = mcd_slope(Ukm1.Er, Ui.Er, Ukp1.Er);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Er  = Ui.Er 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr1 */
  dq1 = mcd_slope(Uim1.Fr1, Ui.Fr1, Uip1.Fr1);
  dq2 = mcd_slope(Ujm1.Fr1, Ui.Fr1, Ujp1.Fr1);
  dq3 = mcd_slope(Ukm1.Fr1, Ui.Fr1, Ukp1.Fr1);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr1 = Ui.Fr1 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr2 */
  dq1 = mcd_slope(Uim1.Fr2, Ui.Fr2, Uip1.Fr2);
  dq2 = mcd_slope(Ujm1.Fr2, Ui.Fr2, Ujp1.Fr2);
  dq3 = mcd_slope(Ukm1.Fr2, Ui.Fr2, Ukp1.Fr2);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr2 = Ui.Fr2 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr3 */
  dq1 = mcd_slope(Uim1.Fr3, Ui.Fr3, Uip1.Fr3);
  dq2 = mcd_slope(Ujm1.Fr3, Ui.Fr3, Ujp1.Fr3);
  dq3 = mcd_slope(Ukm1.Fr3, Ui.Fr3, Ukp1.Fr3);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr3 = Ui.Fr3 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}
#endif /* FIRST_ORDER */
}




#ifndef FIRST_ORDER
static Real mcd_slope(const Real vl, const Real vc, const Real vr){

  Real dvl = (vc - vl), dvr = (vr - vc);
  Real dv, dvm;

  if(dvl > 0.0 && dvr > 0.0){
    dv = 2.0*(dvl < dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm < dv ? dvm : dv);
  }
  else if(dvl < 0.0 && dvr < 0.0){
    dv = 2.0*(dvl > dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm > dv ? dvm : dv);
  }

  return 0.0;
}
#endif /* FIRST_ORDER */




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
	
	int nd, nl, npg, ncg;
	int Nx, Ny, Nz;
	Real temp;

/*=========================================================*/
	int maxCG=1, max1=0, max2=0, max3=0;	
	int sendRC, recvRC, sendP, recvP;
	int max_sendRC=1, max_recvRC=1, max_sendP=1, max_recvP=1; /* maximum size of send buffer */
	
/*===========================================================*/
	

	DomainLevels = pM->NLevels; /* Total levels of refinement above the root domain */

	/* Now consider the coarse grid below the root domain */
	/* First check that root level only has one domain */
	if(pM->DomainsPerLevel[0] > 1)
		ath_error("[BackEuler_init_3d]: Root Level can only have one domain!\n");

	/****************************************/
	/* points to the root domain */
	pD = &(pM->Domain[0][0]);
	pG = pD->Grid;
	
	Nx = pG->ie - pG->is + 1;
	Ny = pG->je - pG->js + 1;
	Nz = pG->ke - pG->ks + 1;

	/* Reach bottom first for the side with the smallest size*/
	RootLevel = Nx;
	if(Ny < Nx) RootLevel = Ny;
	if(Nz < Nx) RootLevel = Nz;

	RootLevel /= Nlim;
	
	temp = log10(RootLevel)/log10(2.0);

	RootLevel = (int)temp;
	if(fabs(temp-RootLevel) > 0.5) RootLevel++;
	
	/* RootLevel is the number of levels below root domain, does not include root domain */
	/* Then Matrix[RootLevel] will be the root domain */
	/****************************************/

	/* The array to store number of domains at each level */
	if((DomainNos = (int*)calloc(DomainLevels+RootLevel,sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");
	
	maxND = 1;
	for(nl=0; nl<DomainLevels; nl++){
		DomainNos[RootLevel+nl] = pM->DomainsPerLevel[nl];
		maxND = MAX(maxND,pM->DomainsPerLevel[nl]);
	}

	if((start_addrP = (int*)calloc_1d_array(maxND,sizeof(int))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate start_addrP\n");


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
	/* Now initialize the Matrix array for each domain */
	for(nl=0; nl<DomainLevels; nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			pD = &(pM->Domain[nl][nd]);
			pG = pD->Grid;
			/* Only do this if this CPU has grid */
			if(pG != NULL){

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
				pMat->BCFlag_ix1 = pM->BCFlag_ix1;
				pMat->BCFlag_ox1 = pM->BCFlag_ox1;
				pMat->BCFlag_ix2 = pM->BCFlag_ix2;
				pMat->BCFlag_ox2 = pM->BCFlag_ox2;
				pMat->BCFlag_ix3 = pM->BCFlag_ix3;
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


	/*==========================================================*/
	/* Allocate memory and initialize parameters for restriction and prolongation */
	/* The maximum memory size required to send data from parent to child grids */
	for(nl=1; nl<DomainLevels; nl++){
		for(nd=0; nd<DomainNos[RootLevels+nl]; nd++){
			sendRC = 0;
			recvRC = 0;
			sendP = 1;
			recvP = 1;
			if(pM->Domain[nl][nd].Grid != NULL){
				pG = pM->Domain[nl][nd].Grid;				
				/* For parent grids; send data for all grids per domain together*/
				for(npg=0; npg<pG->NPGrid; npg++){
					sendRC += pG->PGrid[npg].Rad_nWordsRC;
					recvP  += pG->PGrid[npg].Rad_nWordsP;
				}
				/* For child grids; do all grids per domain together */
				for(ncg=0; ncg<pG->NCGrid; ncg++){
					recvRC += pG->CGrid[ncg].nWordsRC;
					sendP  += pG->CGrid[ncg].nWordsP;
				}
			
				max_sendRC = MAX(max_sendRC, sendRC);
				max_recvRC = MAX(max_recvRC, recvRC);
				max_sendP  = MAX(max_sendP, sendP);
				max_recvP  = MAX(max_recvP, recvP);
				maxCG = MAX(maxCG,pG->NCGrid);
				max1  = MAX(max1, (pG->Nx[0]+1));
				max2  = MAX(max2, (pG->Nx[1]+1));
				max3  = MAX(max3, (pG->Nx[2]+1));
			}/* End if grid is not NULL */

		}/* End domain at each level */
	}/* End level */

	/*======================================================================*/
	/* Allocate memory for variables used for restriction */

	if((send_bufRC = (double**)calloc_2d_array(maxND,max_sendRC,sizeof(double))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate send_bufRC\n");


#ifdef MPI_PARALLEL
	/* we do restriction between every two levels; We will finish one level first and then */
	/* go to another level */ 
	if((recv_bufRC =  (double**)calloc_2d_array(maxND,max_recvRC,sizeof(double))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv_bufRC\n");
	/* We only work for two Levels at a time, so the first index is 2, not total Levels */
	if((recv_rq = (MPI_Request**)calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv MPI_Request array\n");
  	if((send_rq = (MPI_Request**)calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate send MPI_Request array\n");
#endif /* MPI_PARALLEL */

	
	/*======================================================================*/
	
	if((send_bufP =(double**)calloc_2d_array(maxND,max_sendP,sizeof(double))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate send_bufP\n");

	if((recv_bufP =(double**)calloc_2d_array(maxND,max_recvP,sizeof(double))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv_bufP\n");


	max1 += 2*Matghost;
  	max2 += 2*Matghost;
  	max3 += 2*Matghost;

	/* Array to store the prolongation data temporary */
	/* Each GZ[k][j][i][Er-Fr?]: We only need to do the boundary for Er to Fr? */
	if((GZ=(RadMHDS***)calloc_3d_array(max3,max2,max1,sizeof(RadMHDS))) ==NULL) ath_error("[SMR_init]:Failed to allocate GZ[0]C\n");
  	
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

	if(start_addrP != NULL)
		free(start_addrP);

	/*==========================*/
	/* Free the MPI buffer */

	if(send_bufRC != NULL)
		free_2d_array(send_bufRC);

#ifdef MPI_PARALLEL

	if(send_rq != NULL)
		free_2d_array(send_rq);

	if(recv_rq != NULL)
		free_2d_array(recv_rq);

	if(recv_bufRC != NULL)
		free_2d_array(recv_bufRC);

#endif

	if(send_bufP != NULL)
		free_2d_array(send_bufP);

	if(recv_bufP != NULL)
		free_2d_array(recv_bufP);


	/* Free temporary array for prolongation ghost zones */
	
	if(GZ != NULL)
		free_3d_array(GZ);
		
	

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
	RHSResidual3D(pMat, pMat->RHS, &(pMat->RHSnorm0));

	/* Now set the matrix to be the residual */
	/* RHS is already updated with the new RHS in function RHSResidual3D */
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





