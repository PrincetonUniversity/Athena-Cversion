#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Use backward Euler method to update the radiation quantities
 * First set up the matrix
 * Then solve the matrix equations.
 * We need the flag for boundary condition.
 * There are three algorithms now: Jacobi, Gauss-Seidel and Bicgsafe 
 * Bicgsafe prefer to do more iterations during each level, just one W cycle may be enough  
 * 
 * Backward Euler should be used for the whole mesh
 *
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

#ifdef FLD

#ifdef MATRIX_MULTIGRID 

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */



static MatrixS *pMat;

/* memory at each coarse level */
/* Allocate the memory for U, Ugas and RHS in advance */
static RadMHDS ****U_coarse;
static RadCoefS ****Ugas_coarse;
static Real *****RHS_coarse;

static Real ****INIerror;

static Real *****Reshist; /* The residual after each cycle, we only store 5 cycles for iterant recombination */
static Real *****Solhist; /* Store the history of the solution */
static Real *Resnorm; /* store the norm of the stored residual */
static int MAXerror = 1; /* The maximum number of combined residual */
static int Resnum; /* actual number of residual stored, Resnum <= MAXerror */

/* The pointer array to store the address of the coefficient */
/* This make sure that matrix coefficient only calculate once */
static Real *****Ptheta;
static Real *****Pphi;
static Real *****Ppsi;
static Real *****Pvarphi;

static int coefflag; /* To decide whether need to calculate coefficient or not */


static Real INInorm;
/* The initial L2 norm of the right hand side   *
 * The convergence criterian is
 * ||r_n||_2 / ||r_0||_2  < Small number 	*
 */ 

static int Nlim = 4; /* the lim size of the coarsest grid in each CPU*/
static int Wcyclelim = 15; /* The limit of how many Wcycle we allow */ 

static int Matrixflag = 1; /* The matrix flag is used to choose Gauss_Seidel method or Jacobi method, or Bicgsafe */
			  /* 1 is GS method while 0 is Jacobi method , 2 is Bicgsafe method*/

static int Nlevel; /* Number of levels from top to bottom, log2(N)=Nlevel */
static int Wflag; /* To decide the position in the W flag */


/********Private function ************/
static Real CheckResidual(MatrixS *pMat, Real ****vector); /* This function basically calculates the L2 norm of vector */
static void RadMHD_multig_3D(MatrixS *pMat);

static void Restriction3D(MatrixS *pMat_fine, MatrixS *pMat_coarse);

static void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat);


static void RHSResidual3D(MatrixS *pMat, Real ****newRHS);

static void Calculate_Coef(MatrixS *pMat);

static void CopySolution(MatrixS *pMat, Real ****Solhist);

static void Recombination(MatrixS *pMat, Real *****Reshist, Real *****Solhist, const int num, const Real *Norms);

/* calculate inner product of two arbitrary vectors */
static void Inner_product(MatrixS *pMat, Real ****vector1, Real ****vector2, Real *result);

/* New restriction and prolongation function adopted from SMR */

static void prolongation3D(MatrixS *pMat_coarse, MatrixS *pMat_fine);



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
/*-------BackEuler_3d(): Use back euler method to update E_r and Fluxr-----------*/
/* Not work with SMR now. So there is only one Domain in the Mesh */


void BackEuler_3d(MeshS *pM)
{



	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	
	Real dt = pG->dt;

	/* Set the parameters that will change with time in pMat */
	pMat->dt = dt;
	pMat->time = pG->time;

	Real velocity_x, velocity_y, velocity_z, T4, Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF, pressure, density, temperature;
	Real vxi0, vxi1, vyj0, vyj1, vzk0, vzk1,dvxdx,dvydy,dvzdz;
	Real Fr0x, Fr0y, Fr0z;
	Real f11, f22, f33, f21, f32, f31;
	Real dErdx, dErdy, dErdz, divEr, limiter;
    Real AdvFx[2], AdvFy[2], AdvFz[2];

	/*Sigma_aF: flux mean absorption, Sigma_aP: Plank mean absorption, Sigma_aE: Er mean absorption opacity;*/
	Real Sigma[NOPACITY];
	/* source coefficient for radiation work term */	

	Real error;
	int Wcycle;
	
	int i, j, k, m;
	int is, ie, js, je, ks, ke;
	int Mati, Matj, Matk;
	/* Set the boundary */
	is = pG->is-Matghost;
	ie = pG->ie+Matghost;
	js = pG->js-Matghost;
	je = pG->je+Matghost;
	ks = pG->ks-Matghost;
	ke = pG->ke+Matghost;
	
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

/*===========================*/
/* variables used to subtract background state */
/*
	static int t0flag = 1;
	int bgflag;	

	bgflag = 0; 
*/
/*===================================*/

	/* Apply the boundary condition */
	/* Need to update the boundary condition of Er before calculating Fr */
/*        for (i=0; i<pM->NLevels; i++){
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
                if (pM->Domain[i][j].Grid != NULL){
                        bvals_radMHD(&(pM->Domain[i][j]));

                }                    }
        }  
*/

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
				pMat->Ugas[Matk][Matj][Mati].rho = pG->U[k][j][i].d;
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

				/* calculate the FLD limiter */
				dErdx = (pG->U[k][j][i+1].Er - pG->U[k][j][i-1].Er) / (2.0 * pG->dx1); 
				dErdy = (pG->U[k][j+1][i].Er - pG->U[k][j-1][i].Er) / (2.0 * pG->dx2);
				dErdz = (pG->U[k+1][j][i].Er - pG->U[k-1][j][i].Er) / (2.0 * pG->dx3);

				divEr = sqrt(dErdx * dErdx + dErdy * dErdy + dErdz * dErdz);

				FLD_limiter(divEr, pG->U[k][j][i].Er, Sigma_sF + Sigma_aF, &(pMat->Ugas[Matk][Matj][Mati].lambda)); 
				

		/* Now set the right hand side */
				Rad_Advection_Flux3D(pD, i, j, k, 1.0, AdvFx, AdvFy, AdvFz);
						
				pMat->RHS[Matk][Matj][Mati][0] = pG->U[k][j][i].Er + ((AdvFx[1] - AdvFx[0]) + (AdvFy[1] - AdvFy[0]) + (AdvFz[1] - AdvFz[0])) + dt * Sigma_aP * T4 * Crat * Eratio;
				/* Only need to calculate Er, one equation */
				
	} /* End i */
	}/* End j */
	}/* End k */

	/* Update boundary condition for gas quantities */
	/* The top level is always there, do not need to destruct every time */
	

/* First, calculate the coefficient for the top level  */
/* coefficient will not change until next time step */
	Calculate_Coef(pMat);


/* Only do this if background state is subtracted */
/* For inflow boundary condition, if background state is subtracted, ghost zones in the matrix should always set to be zero */
/* Bicgsafe always work with residual */
if((pMat->bgflag) || (Matrixflag == 2)){
	/* calculate the residual */
	RHSResidual3D(pMat, INIerror);

	/* Now set the matrix to be the residual */
	for(k=pMat->ks - Matghost; k<=pMat->ke+Matghost; k++){
		for(j=pMat->js-Matghost; j<=pMat->je+Matghost; j++){
			for(i=pMat->is-Matghost; i<= pMat->ie+Matghost; i++){
				pMat->U[k][j][i].Er = 0.0;
				

				pMat->RHS[k][j][i][0] = INIerror[k][j][i][0];	
				
			}
		}
	}
	/* Right hand side in the ghost zones are never used */
}

	/* calculate the initial norm of the right hand side */
	INInorm = CheckResidual(pMat, pMat->RHS);

	if(INInorm < TINY_NUMBER)
		return;

		
	/* Do the multi-grid W cycle */
	/* reset the flags for each cycle */
	Wcycle = 0;
	error = 1.0;
	Resnum = 0;
	coefflag = 1;
	/* Stop iteration if tolerance level is reached */
	/* Or the matrix doesn't converge in Wcycle limit */
	while((error > TOL) && (Wcycle < Wcyclelim)){
		Wflag = 1;

		/* Do the recombination before the multigrid cycle */
		/* num from 0... Resnum-1 */
		if(Resnum > 1){
			Recombination(pMat, Reshist, Solhist, Resnum, Resnorm);
			/* After recombination, replace solhist with improved solution */
		/*	Resnorm[Resnum-1] = CheckResidual(pMat,Reshist[Resnum-1]);
			
			CopySolution(pMat,Solhist[Resnum-1]);
		*/	
		}

		
		if(Resnum >= MAXerror)
			Resnum = 0;

		RadMHD_multig_3D(pMat);

		/* for parent grid, check residual */
		/* Bicgsafe works different from GS and Jacobi */
		/* pMat->RHS is changed during each iteration for Bicgsafe */
		if((Matrixflag == 0) || (Matrixflag == 1)){
			/* calculate the residual and store the residual in Reshist */
			RHSResidual3D(pMat, Reshist[Resnum]);			
			error = CheckResidual(pMat,Reshist[Resnum]);
			/* store the current solution and its norm */
			CopySolution(pMat,Solhist[Resnum]);
			Resnorm[Resnum] = error;
			
		}
		else if(Matrixflag == 2){
			error = CheckResidual(pMat,pMat->RHS);
		}	
		else{
			ath_error("[BackEuler3D: unknown matrix solver!\n]");
		}
		

		error /= INInorm;

		if(error != error)
			ath_error("[BackEuler3D]: NaN encountered!\n");
	
		Wcycle++;
		/* only need to calculate the coefficient in the first cycle */
		coefflag = 0;

	}
		/* Only output the residual for parent grid */	
		if(myID == 0)
			ath_pout(0,"Final residual: %e  Cycle No.: %d\n",error,Wcycle);

		
	/* add the correction back */
	/* Now copy the data */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<= ie; i++){
				

				Mati = i - (nghost - Matghost);
				Matj = j - (nghost - Matghost);
				Matk = k - (nghost - Matghost);

				T4 = pG->U[k][j][i].Er;

			if(pMat->bgflag){
				pG->U[k][j][i].Er += pMat->U[Matk][Matj][Mati].Er;
				
			}
			else{
				pG->U[k][j][i].Er = pMat->U[Matk][Matj][Mati].Er;
			
			}
					
			
/*			if(pG->U[k][j][i].Er < TINY_NUMBER){
                                pG->U[k][j][i].Er = (pG->U[k][j][i-1].Er + pG->U[k][j][i+1].Er + pG->U[k][j-1][i].Er + pG->U[k][j+1][i].Er +  pG->U[k-1][j][i].Er + pG->U[k+1][j][i].Er) / 6.0;

                         }
*/
			
	}
    
    /* Do the Div V P step */
    
    for(k=ks; k<=ke; k++){
        for(j=js; j<=je; j++){
            for(i=is; i<= ie; i++){
				f11 = pG->U[k][j][i].Edd_11;
                f22 = pG->U[k][j][i].Edd_22;
				f33 = pG->U[k][j][i].Edd_33;
                
                f21 = pG->U[k][j][i].Edd_21;
                f31 = pG->U[k][j][i].Edd_31;
                f32 = pG->U[k][j][i].Edd_32;
                
				vxi0 = pG->U[k][j][i-1].M1 / pG->U[k][j][i-1].d;
				vxi1 = pG->U[k][j][i+1].M1 / pG->U[k][j][i+1].d;
				vyj0 = pG->U[k][j-1][i].M2 / pG->U[k][j-1][i].d;
				vyj1 = pG->U[k][j+1][i].M2 / pG->U[k][j+1][i].d;
				vzk0 = pG->U[k-1][j][i].M3 / pG->U[k-1][j][i].d;
                vzk1 = pG->U[k+1][j][i].M3 / pG->U[k+1][j][i].d;
				
				dvxdx = (vxi1 - vxi0) / (2.0 * pG->dx1);
				dvydy = (vyj1 - vyj0) / (2.0 * pG->dx2);
				dvzdz = (vzk1 - vzk0) / (2.0 * pG->dx3);
				
				pG->U[k][j][i].Er /= (1.0 + pG->dt * ((f11 + f21 + f31) * dvxdx + (f21 + f22 + f32) * dvydy + (f31 + f32 + f33) * dvzdz));
			}
		}
	}


/* Need to update the boundary condition of Er before calculating Fr */
	for (i=0; i<pM->NLevels; i++){ 
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
        	if (pM->Domain[i][j].Grid != NULL){
  			bvals_radMHD(&(pM->Domain[i][j]));

        	}
      	     }
    	}

/* Update the Eddington tensor */
	/* calculate the Eddington tensor first */
	Eddington_FUN(&(pM->Domain[0][0]));

	/* Now update the cell centered  flux */
	for(k=pG->ks; k<=pG->ke; k++)
		for(j=pG->js; j<=pG->je; j++)
			for(i=pG->is; i<= pG->ie; i++){

				Mati = i - (nghost - Matghost);
				Matj = j - (nghost - Matghost);
				Matk = k - (nghost - Matghost);


	/* calculate the FLD limiter */
				dErdx = (pG->U[k][j][i+1].Er - pG->U[k][j][i-1].Er) / (2.0 * pG->dx1); 
				dErdy = (pG->U[k][j+1][i].Er - pG->U[k][j-1][i].Er) / (2.0 * pG->dx2);
				dErdz = (pG->U[k+1][j][i].Er - pG->U[k-1][j][i].Er) / (2.0 * pG->dx3);

				divEr = sqrt(dErdx * dErdx + dErdy * dErdy + dErdz * dErdz);

				Sigma_sF = pG->U[k][j][i].Sigma[0];
				Sigma_aF = pG->U[k][j][i].Sigma[1];
			

				FLD_limiter(divEr, pG->U[k][j][i].Er, Sigma_sF + Sigma_aF, &(limiter)); 
				
				velocity_x = pMat->Ugas[Matk][Matj][Mati].V1;
                		velocity_y = pMat->Ugas[Matk][Matj][Mati].V2;
                		velocity_z = pMat->Ugas[Matk][Matj][Mati].V3;

				/* First, the co-moving flux */
				pG->U[k][j][i].Fr1 = -limiter * dErdx / (Sigma_sF + Sigma_aF);
				pG->U[k][j][i].Fr2 = -limiter * dErdy / (Sigma_sF + Sigma_aF); 
				pG->U[k][j][i].Fr3 = -limiter * dErdz / (Sigma_sF + Sigma_aF);

				/* Now add the advection part */	
				pG->U[k][j][i].Fr1 += ((1.0 +  pG->U[k][j][i].Edd_11) * velocity_x +  pG->U[k][j][i].Edd_21 * velocity_y + pG->U[k][j][i].Edd_31 * velocity_z) *  pG->U[k][j][i].Er / Crat; 
				pG->U[k][j][i].Fr2 += ((1.0 +  pG->U[k][j][i].Edd_22) * velocity_y +  pG->U[k][j][i].Edd_21 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_z) *  pG->U[k][j][i].Er / Crat;
				pG->U[k][j][i].Fr3 += ((1.0 +  pG->U[k][j][i].Edd_33) * velocity_z +  pG->U[k][j][i].Edd_31 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_y) *  pG->U[k][j][i].Er / Crat; 
			}



/*
if(Opacity != NULL){
		for (k=pG->ks; k<=pG->ke; k++){
			for (j=pG->js; j<=pG->je; j++) {
    				for (i=pG->is; i<=pG->ie; i++){
				
				density = pG->U[k][j][i].d;
				
				pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
				+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density ) * (Gamma - 1);
			
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

				if(pressure > TINY_NUMBER)
				{
					temperature = pressure / (density * R_ideal);

						
					Opacity(density,temperature,Sigma,NULL);
					for(m=0;m<NOPACITY;m++){
						pG->U[k][j][i].Sigma[m] = Sigma[m];
					}

				}
				else{
					
					pG->Tguess[k][j][i] = pG->U[k][j][i].Er;
				}

			
				}
			}
		}
}
*/

	/* Set the boundary condition */

	
/* Update the ghost zones for different boundary condition to be used later */
/*	for (i=0; i<pM->NLevels; i++){ 
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
        	if (pM->Domain[i][j].Grid != NULL){
  			bvals_radMHD(&(pM->Domain[i][j]));

        	}
      	     }
    	}
*/
  return;	
	

}


void RadMHD_multig_3D(MatrixS *pMat)
{

	MatrixS pMat_coarse;

	/* Once any dimension reaches size limit of Nlim, do Ncycle iteration and return */
	if(pMat->Nx[0] <= Nlim || pMat->Nx[1] <= Nlim || pMat->Nx[2] <= Nlim){
		/* Create the temporary array */
		/* Need to create the temporary array for the boundary condition */		


		bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
		bvals_Matrix_shear_init(pMat);
#endif

		
		if(Matrixflag == 1){
			GaussSeidel3D(pMat,Ptheta[pMat->Level],NULL,NULL,NULL);
		}
		else if(Matrixflag == 0){
			Jacobi3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
		}
		else if(Matrixflag == 2){
			Bicgsafe3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
		}
		else
			ath_error("[BackEuler3D: unknown matrix solver!\n]");

		
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
		Restriction3D(pMat, &pMat_coarse);

		if(coefflag){

		/* Only need to update ghost zones for gas quantities */
		/* after first restriction, only right hand side changes */
		bvals_Matrix_gas_init(&pMat_coarse);
#ifdef SHEARING_BOX
		bvals_Matrix_shear_gas_init(&pMat_coarse);
#endif
		
		bvals_Matrix_gas(&pMat_coarse);
		
#ifdef SHEARING_BOX
		bvals_Matrix_shear_gas_destruct();
#endif


		bvals_Matrix_gas_destruct(&pMat_coarse);

		/* first time, calculate the coefficient */
		
			Calculate_Coef(&(pMat_coarse));

		}

		/* continue the V or W cycle recursively*/

		RadMHD_multig_3D(&pMat_coarse);

/* The following code is first reached after Ncycle iterations at the coarsest
 * level.  We then prolongate, do Ncycle iterations, and return.  This will return
 * execution to this same spot for the next coarsest level, so we will
 * prolongate, do Ncycle iterations, return, and so on.
 */

		/* Add the correction back to the fine grid */
		prolongation3D(&pMat_coarse, pMat);
		
		
		/* Do relaxation when going up*/
			
		
		bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
		bvals_Matrix_shear_init(pMat);
#endif

		/* Update the ghost zones first */
		bvals_Matrix(pMat);


		/* update the residual after prolongation*/
		/* It is no longer the original sequence after prolongation */
		if(Matrixflag == 2){
			RHSResidual3D(pMat, pMat->RHS);
		}

		
		if(Matrixflag == 1){
			GaussSeidel3D(pMat,Ptheta[pMat->Level],NULL,NULL,NULL);
		}
		else if(Matrixflag == 0){
			Jacobi3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
		}
		else if(Matrixflag == 2){
			Bicgsafe3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
		}
		else
			ath_error("[BackEuler3D: unknown matrix solver!\n]");
		
#ifdef SHEARING_BOX
		bvals_Matrix_shear_destruct();
#endif


		bvals_Matrix_destruct(pMat);
	
		/* To decide whether go W cycle */
		
		if((pMat->Level < Nlevel/2.0) && Wflag){
			/* Set is Wcycle flag */
			/* This is Wcycle. We should go down again */
			Wflag = 0;	

			set_mat_level(&(pMat_coarse), pMat);

		
			/* project the data to coarse grid */
			Restriction3D(pMat, &pMat_coarse);

			/* Do not need to update ghost zones inside W cycle */
			/* gas quantities are not changed, matrix coefficients are not changed */
			/* radiation quantities are always zero */
			
/*			bvals_Matrix_init(&pMat_coarse);

#ifdef SHEARING_BOX
			bvals_Matrix_shear_init(&pMat_coarse);
#endif
		
			bvals_Matrix(&pMat_coarse);
		
#ifdef SHEARING_BOX
			bvals_Matrix_shear_destruct();
#endif

			bvals_Matrix_destruct(&pMat_coarse);
*/
			/* continue the V or W cycle recursively*/

			RadMHD_multig_3D(&pMat_coarse);


			/* Add the correction back to the fine grid */
			prolongation3D(&pMat_coarse, pMat);

			/* Do relaxation when going up*/
			
		
			bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
			bvals_Matrix_shear_init(pMat);
#endif
		

			/* Update the ghost zones first */
			bvals_Matrix(pMat);

			/* update the residual after prolongation*/
			/* It is no longer the original sequence after prolongation */
			if(Matrixflag == 2){
				RHSResidual3D(pMat, pMat->RHS);
			}

			


			if(Matrixflag == 1){
				GaussSeidel3D(pMat,Ptheta[pMat->Level],NULL,NULL,NULL);
			}
			else if(Matrixflag == 0){
				Jacobi3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
			}
			else if(Matrixflag == 2){
				Bicgsafe3D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level],Pvarphi[pMat->Level]);
			}
			else
				ath_error("[BackEuler3D: unknown matrix solver!\n]");
		
#ifdef SHEARING_BOX
			bvals_Matrix_shear_destruct();
#endif


			bvals_Matrix_destruct(pMat);

		}


		


	}


	return;
}



/*==========================================*/
/* Check the residual after one cycle */
/* This is the stop criterian */
Real CheckResidual(MatrixS *pMat, Real ****vector)
{

	/* This function calculate the norm of the right hand side */
	Real Norm;

	Inner_product(pMat, vector, vector, &Norm);

	return Norm;

}



/* prolongation operator */


void prolongation3D(MatrixS *pMat_coarse, MatrixS *pMat_fine)
{
	
	int i, j, k, ifine, jfine, kfine;
	int ii, jj, kk;
	int is, ie, js, je, ks, ke;

	RadMHDS Ptemp[2][2][2];	

	is = pMat_coarse->is;
	ie = pMat_coarse->ie;
	js = pMat_coarse->js;
	je = pMat_coarse->je;
	ks = pMat_coarse->ks;
	ke = pMat_coarse->ke;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){

				ProU(pMat_coarse->U[k][j][i-1],pMat_coarse->U[k][j][i],pMat_coarse->U[k][j][i+1],
			 		pMat_coarse->U[k][j-1][i],pMat_coarse->U[k][j+1][i],
					pMat_coarse->U[k-1][j][i],pMat_coarse->U[k+1][j][i],Ptemp);


			/* Now copy the data to the fine grid */
				ii = 2*(i-is) + pMat_fine->is;
				jj = 2*(j-js) + pMat_fine->js;
				kk = 2*(k-ks) + pMat_fine->ks;

				/* The coarse grid calculates the residual, we add it back */
				for(kfine=0;kfine<2; kfine++)
					for(jfine=0; jfine<2; jfine++)
						for(ifine=0; ifine<2; ifine++){
							pMat_fine->U[kk+kfine][jj+jfine][ii+ifine].Er  += Ptemp[kfine][jfine][ifine].Er;
							
				}				

	}

}


/* Restriction operator */

void Restriction3D(MatrixS *pMat_fine, MatrixS *pMat_coarse)
{
/* We actually should send the residual to the next level */

/* After restriction step, we do not need to set the boundary values */
/* values of v, T, sigma in the ghost zones are not used */
/* Initial guess for Er , Fr1, Fr2, Fr3 are always zero */
/* Boundary values will be updated in GaussSeidal */

	int i, j, k;
	/* The right hand size is no-longer the original source terms. It is the residual  */
	/* But we cannot destroy the right hand size of the fine grid. We need that for the */
	/* relaxation step when we come back */ 


	Real ****error;	
	if((error=(Real****)calloc_4d_array(pMat_fine->Nx[2]+2*Matghost,pMat_fine->Nx[1]+2*Matghost,pMat_fine->Nx[0]+2*Matghost,1,sizeof(Real))) == NULL)
			ath_error("[Restriction3D]: malloc return a NULL pointer\n");		

	Real *ptr_coarse;
	Real *ptr_fine[8];
	int num;	

	/* calculate the residual */
	RHSResidual3D(pMat_fine, error);

		/* error = b - Ax */
		/* Now restrict to the coarse grid */
		/* only need to restrict the gas quantities once */
		/* But we need to restrict residual every time */
	
	for(k=pMat_coarse->ks; k<=pMat_coarse->ke; k++)
		for(j=pMat_coarse->js; j<=pMat_coarse->je; j++)
			for(i=pMat_coarse->is; i<=pMat_coarse->ie; i++){
			if(coefflag){
				ptr_coarse  = &(pMat_coarse->Ugas[k][j][i].rho);
				ptr_fine[0] = &(pMat_fine->Ugas[2*k ][2*j ][2*i ].rho);
				ptr_fine[1] = &(pMat_fine->Ugas[2*k ][2*j ][2*i-1].rho);
				ptr_fine[2] = &(pMat_fine->Ugas[2*k ][2*j-1][2*i ].rho);
				ptr_fine[3] = &(pMat_fine->Ugas[2*k ][2*j-1][2*i-1].rho);
				ptr_fine[4] = &(pMat_fine->Ugas[2*k-1][2*j ][2*i ].rho);
				ptr_fine[5] = &(pMat_fine->Ugas[2*k-1][2*j ][2*i-1].rho);
				ptr_fine[6] = &(pMat_fine->Ugas[2*k-1][2*j-1][2*i ].rho);	
				ptr_fine[7] = &(pMat_fine->Ugas[2*k-1][2*j-1][2*i-1].rho);

				/* 12 variables in RadCoef_S */ 
				for(num=0; num<12+NOPACITY; num++){

					ptr_coarse[num] =  (ptr_fine[0][num] + ptr_fine[1][num]
							 + ptr_fine[2][num] + ptr_fine[3][num]
							 + ptr_fine[4][num] + ptr_fine[5][num]
							 + ptr_fine[6][num] + ptr_fine[7][num]) / 8.0;

				}
			}
				
			/*
			pMat_coarse->U[k][j][i] = (pMat_fine->U[2*k ][2*j ][2*i ]  + pMat_fine->U[2*k ][2*j ][2*i-1]
						+ pMat_fine->U[2*k ][2*j-1][2*i ] + pMat_fine->U[2*k ][2*j-1][2*i-1]
						+ pMat_fine->U[2*k-1][2*j ][2*i ] + pMat_fine->U[2*k-1][2*j ][2*i-1]
						+ pMat_fine->U[2*k-1][2*j-1][2*i ]+ pMat_fine->U[2*k-1][2*j-1][2*i-1]) / 8.0;

			*/
			pMat_coarse->RHS[k][j][i][0] =  (error[2*k ][2*j ][2*i ][0] + error[2*k ][2*j ][2*i-1][0]
						     +	error[2*k ][2*j-1][2*i ][0]+ error[2*k ][2*j-1][2*i-1][0]
						     +  error[2*k-1][2*j ][2*i ][0]+ error[2*k-1][2*j ][2*i-1][0]
						     +	error[2*k-1][2*j-1][2*i ][0]+ error[2*k-1][2*j-1][2*i-1][0]) / 8.0;

		
	}

	/* Also set ghost zones to be zero */
	for(k=pMat_coarse->ks-Matghost; k<=pMat_coarse->ke+Matghost; k++)
		for(j=pMat_coarse->js-Matghost; j<=pMat_coarse->je+Matghost; j++)
			for(i=pMat_coarse->is-Matghost; i<=pMat_coarse->ie+Matghost; i++){

				/* The initial guess is taken to be zero */
				pMat_coarse->U[k][j][i].Er = 0.0;
				
		}
		

		/* Free the temporary array */
		free_4d_array(error);

	/* We need to update the Eddington tensor in the ghost zone */

	return;
}


/*  calculate the residual of right hand side from guess solution */
/* should allocate memory for newRHS before this function is called */
void RHSResidual3D(MatrixS *pMat, Real ****newRHS)
{

	int i, j, k, m, n;
	int is, ie, js, je, ks, ke;

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;


	
	Real tempEr1, tempEr2, tempEr3;	
	Real temp0;

	/* To store the coefficient */
	Real theta[7];
	

	/* current level */
	n = pMat->Level;

	
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				/* get the coefficient */
				for(m=0; m<7; m++){
					theta[m]  = Ptheta[n][k][j][i][m];					
				}

			newRHS[k][j][i][0] = pMat->RHS[k][j][i][0];

			tempEr3 = theta[0] * pMat->U[k-1][j][i].Er + theta[6] * pMat->U[k+1][j][i].Er;
			tempEr2 = theta[1] * pMat->U[k][j-1][i].Er + theta[5] * pMat->U[k][j+1][i].Er;
			tempEr1 = theta[2] * pMat->U[k][j][i-1].Er + theta[4] * pMat->U[k][j][i+1].Er;

			temp0 = theta[3] * pMat->U[k][j][i].Er;

			newRHS[k][j][i][0] -= (temp0 + (tempEr1 + tempEr2 + tempEr3));


	}

		/* error = b - Ax */

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

				matrix_coef_FLD(pMat, 3, i, j, k, &(Ptheta[n][k][j][i][0]));
		}

}
void Recombination(MatrixS *pMat, Real *****Reshist, Real *****Solhist, const int num, const Real *Norms)
{
	/* combine the latest 0 ... num-1 solution to minimize the current residual */
	/* Norms[i] = <Reshist[i], Reshist[i]> */

	int i, j, k, m;	
	int is, ie, js, je, ks, ke;
	int flag = 1; /* used to label whether matrix inversion sucesseed or not */

	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

	Real **Matrixcoef = NULL;
	Real *MatrixRHS = NULL;
	Real *crossNorm = NULL;

	/* temporary used for LU decomposition */
	int *index = NULL;
	Real *dtemp = NULL;
	

	Real tempNorm;

	if((Matrixcoef=(Real**)calloc_2d_array(num,num,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((MatrixRHS=(Real*)calloc(num,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");
	
	if((crossNorm=(Real*)calloc(num,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((index=(int*)calloc(num,sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((dtemp=(Real*)calloc(num,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	

	/* set the coefficient of H matrix and right hand side */
	/* crossNorm is <Res[num-1], Res[i]> */

	for(i=1; i<=num-1; i++){
		Inner_product(pMat, Reshist[num-1-i], Reshist[num-1], &(crossNorm[i]));
		MatrixRHS[i] = Norms[num-1] - crossNorm[i];
	}


	for(i=1; i<= num-1; i++){
		for(j=1; j<= num-1; j++){
			if(i != j){
				Inner_product(pMat, Reshist[num-1-i], Reshist[num-1-j], &(tempNorm));
				Matrixcoef[i][j] = tempNorm - crossNorm[i] - crossNorm[j] + Norms[num-1];
			}
			else{
				Matrixcoef[i][j] = Norms[num-1-i] - 2.0 * crossNorm[i] +  Norms[num-1];

			}
		}
	}

	/* now solve the (num-1)(num-1) matrix:  Matrixcoef alphas = MatrixRHS */ 
	/* This is smaller than 5 * 5 matrix, can use direct method to solve */
	/* LU decomposition */
	ludcmpnew(Matrixcoef,num-1, index, dtemp, &flag);	
	/* If matrix is ill conditinoed, do not do recombination */
	/* solve the matrix */
	if(flag){

	lubksb(Matrixcoef, num-1, index, MatrixRHS);
	/* solution is stored in MatrixRHS */
	
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
				for(m=1; m<num; m++){
					pMat->U[k][j][i].Er  += MatrixRHS[m] * (Solhist[num-1-m][k][j][i][0] - Solhist[num-1][k][j][i][0]);
					pMat->U[k][j][i].Fr1 += MatrixRHS[m] * (Solhist[num-1-m][k][j][i][1] - Solhist[num-1][k][j][i][1]);
					pMat->U[k][j][i].Fr2 += MatrixRHS[m] * (Solhist[num-1-m][k][j][i][2] - Solhist[num-1][k][j][i][2]);
					pMat->U[k][j][i].Fr3 += MatrixRHS[m] * (Solhist[num-1-m][k][j][i][3] - Solhist[num-1][k][j][i][3]);
				}
	}

	}
	/* The more accurate solution is now stored in pMat */




	if(Matrixcoef != NULL)
		free_2d_array(Matrixcoef);

	if(MatrixRHS != NULL)
		free(MatrixRHS);

	if(crossNorm != NULL)
		free(crossNorm);

	
	if(index != NULL)
		free(index);

	
	if(dtemp != NULL)
		free(dtemp);



}

/* Copy and store current estimate solution */
void CopySolution(MatrixS *pMat, Real ****Solhist)
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

				Solhist[k][j][i][0] = pMat->U[k][j][i].Er;
				Solhist[k][j][i][1] = pMat->U[k][j][i].Fr1;
				Solhist[k][j][i][2] = pMat->U[k][j][i].Fr2;
				Solhist[k][j][i][3] = pMat->U[k][j][i].Fr3;
	}

	return;
}

/* calculate the inner product of two vectors , handle MPI case */
void Inner_product(MatrixS *pMat, Real ****vector1, Real ****vector2, Real *result)
{


	Real Norm;
	Real normtemp;

	int i, j, k;
	int is, ie, js, je, ks, ke;
	is = pMat->is;
	ie = pMat->ie;
	js = pMat->js;
	je = pMat->je;
	ks = pMat->ks;
	ke = pMat->ke;

#ifdef MPI_PARALLEL
	int ierr;
	double tot_norm = 0.0;
#endif

	Norm = 0.0;
	/* There is only Er, no Fr */
	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){			
				vector_product(vector1[k][j][i],vector2[k][j][i], 1 , &normtemp);
				Norm += normtemp;							
			}	
		
	/* MPI call to collect norm from different CPUs */


#ifdef MPI_PARALLEL
	ierr = MPI_Allreduce(&Norm,&tot_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

	Norm = tot_norm;  
#endif	

	*result = Norm;

	return;

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





void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat)
{
    
#ifdef MPI_PARALLEL
	pMat_coarse->Comm_Domain = pMat->Comm_Domain;
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
	pMat_coarse->Level = pMat->Level + 1;
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
		
	pMat_coarse->ID = pMat->ID;
	pMat_coarse->my_iproc = pMat->my_iproc;
	pMat_coarse->my_jproc = pMat->my_jproc; 
	pMat_coarse->my_kproc = pMat->my_kproc;

	pMat_coarse->bgflag = pMat->bgflag;

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

	/* memory is already relocated. Just set the pointer */

	pMat_coarse->U    = U_coarse[pMat_coarse->Level];
	pMat_coarse->Ugas = Ugas_coarse[pMat_coarse->Level];
	pMat_coarse->RHS  = RHS_coarse[pMat_coarse->Level];



	return;
}




/*-------------------------------------------------------------------------*/
/* BackEuler_init_2d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_2d(): function to free memory */
void BackEuler_init_3d(MeshS *pM)
{

	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	int Nx, Ny, Nz;
	int Nx2, Ny2, Nz2;
	int i;
	Real temp;
	

	Nx = pG->ie - pG->is + 1;
	Ny = pG->je - pG->js + 1;
	Nz = pG->ke - pG->ks + 1;

	/* Reach bottom first for the side with the smallest size*/
	Nlevel = Nx;
	if(Ny < Nx) Nlevel = Ny;
	if(Nz < Nx) Nlevel = Nz;

	Nlevel /= Nlim;
	
	temp = log10(Nlevel)/log10(2.0);

	Nlevel = (int)temp;
	if(fabs(temp-Nlevel) > 0.5) Nlevel++;
	Nlevel++;
	

	/* pMat will remain in the memory until the end of the simulation */

	if((pMat = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((pMat->U = (RadMHDS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadMHDS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((pMat->Ugas = (RadCoefS***)calloc_3d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadCoefS))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((pMat->RHS = (Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost, Nx+2*Matghost,4,sizeof(Real))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

	if((INIerror=(Real****)calloc_4d_array(Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,4,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");	


	if((Reshist=(Real*****)calloc_5d_array(MAXerror,Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,4,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((Solhist=(Real*****)calloc_5d_array(MAXerror,Nz+2*Matghost,Ny+2*Matghost,Nx+2*Matghost,4,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((Resnorm=(Real*)calloc(MAXerror,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");	

	/*==================================================================*/
	/* The pointer to store the pointer of the coefficients at each level */
	if((Ptheta=(Real*****)calloc(Nlevel,sizeof(Real****))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");
	




	if((U_coarse=(RadMHDS****)calloc(Nlevel,sizeof(RadMHDS***))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((Ugas_coarse=(RadCoefS****)calloc(Nlevel,sizeof(RadCoefS***))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

	if((RHS_coarse=(Real*****)calloc(Nlevel,sizeof(Real****))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

#ifdef MPI_PARALLEL
	pMat->Comm_Domain = pD->Comm_Domain;
#endif
    

	/* allocate memory at each level */
	Nx2 = Nx;
	Ny2 = Ny;
	Nz2 = Nz;

	for(i=0; i<Nlevel; i++){

		if((Ptheta[i]=(Real****)calloc_4d_array(Nz2+2*Matghost,Ny2+2*Matghost,Nx2+2*Matghost,7,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");
		

		if((U_coarse[i] = (RadMHDS***)calloc_3d_array(Nz2+2*Matghost,Ny2+2*Matghost, Nx2+2*Matghost,sizeof(RadMHDS))) == NULL)
			ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

		if((Ugas_coarse[i] = (RadCoefS***)calloc_3d_array(Nz2+2*Matghost,Ny2+2*Matghost, Nx2+2*Matghost,sizeof(RadCoefS))) == NULL)
			ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

		if((RHS_coarse[i] = (Real****)calloc_4d_array(Nz2+2*Matghost,Ny2+2*Matghost, Nx2+2*Matghost,4,sizeof(Real))) == NULL)
			ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


		Nx2 /= 2;
		Ny2 /= 2;
		Nz2 /= 2;
	}
	/*==========================================================================*/
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
	pMat->js = Matghost;
	pMat->je = Ny + Matghost - 1;
	pMat->ks = Matghost;
	pMat->ke = Nz + Matghost - 1;
	pMat->Nx[0] = Nx;
	pMat->Nx[1] = Ny;
	pMat->Nx[2] = Nz;
	pMat->RootNx[0] = Nx;
	pMat->RootNx[1] = Ny;
	pMat->RootNx[2] = Nz;
	pMat->Level = 0;
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


	


}


void BackEuler_destruct_3d(MeshS *pM)
{
	int i;
	/* Free pMat and pMatnew */


	if(pMat->U != NULL)
		free_3d_array(pMat->U);

	if(pMat->Ugas != NULL)
		free_3d_array(pMat->Ugas);

	if(pMat->RHS != NULL)
		free_4d_array(pMat->RHS);

	if(pMat != NULL)
		free(pMat);

	if(INIerror != NULL)
		free_4d_array(INIerror);

	if(Reshist != NULL)
		free_5d_array(Reshist);

	if(Solhist != NULL)
		free_5d_array(Solhist);

	if(Resnorm != NULL)
		free(Resnorm);


	for(i=0; i<Nlevel; i++){
		free_4d_array(Ptheta[i]);
		

		free_3d_array(U_coarse[i]);
		free_3d_array(Ugas_coarse[i]);
		free_4d_array(RHS_coarse[i]);		
	}		

	free(Ptheta);
	
	free(U_coarse);
	free(Ugas_coarse);
	free(RHS_coarse);

}




#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */


#endif /* FLD */


/*

void prolongation3D(MatrixS *pMat_coarse, MatrixS *pMat_fine)
{
	
	int i, j, k;
	int is, ie, js, je, ks, ke;
	int num;

	Real *ptr_fine[8];
	Real *ptr_coarse[9];

	is = pMat_coarse->is;
	ie = pMat_coarse->ie;
	js = pMat_coarse->js;
	je = pMat_coarse->je;
	ks = pMat_coarse->ks;
	ke = pMat_coarse->ke;

	for(k=ks; k<=ke; k++)
		for(j=js; j<=je; j++)
			for(i=is; i<=ie; i++){
			
				ptr_fine[0] = &(pMat_fine->U[2*k  ][2*j  ][2*i  ].Er);
				ptr_fine[1] = &(pMat_fine->U[2*k  ][2*j  ][2*i-1].Er);
				ptr_fine[2] = &(pMat_fine->U[2*k  ][2*j-1][2*i  ].Er);
				ptr_fine[3] = &(pMat_fine->U[2*k-1][2*j  ][2*i  ].Er);
				ptr_fine[4] = &(pMat_fine->U[2*k  ][2*j-1][2*i-1].Er);
				ptr_fine[5] = &(pMat_fine->U[2*k-1][2*j-1][2*i  ].Er);
				ptr_fine[6] = &(pMat_fine->U[2*k-1][2*j  ][2*i-1].Er);
				ptr_fine[7] = &(pMat_fine->U[2*k-1][2*j-1][2*i-1].Er);

				ptr_coarse[0] = &(pMat_coarse->U[k ][j ][i ].Er);
				ptr_coarse[1] = &(pMat_coarse->U[k+1][j+1][i+1].Er);
				ptr_coarse[2] = &(pMat_coarse->U[k+1][j+1][i-1].Er);
				ptr_coarse[3] = &(pMat_coarse->U[k+1][j-1][i+1].Er);
				ptr_coarse[4] = &(pMat_coarse->U[k-1][j+1][i+1].Er);
				ptr_coarse[5] = &(pMat_coarse->U[k+1][j-1][i-1].Er);
				ptr_coarse[6] = &(pMat_coarse->U[k-1][j-1][i+1].Er);
				ptr_coarse[7] = &(pMat_coarse->U[k-1][j+1][i-1].Er);
				ptr_coarse[8] = &(pMat_coarse->U[k-1][j-1][i-1].Er);
				

				for(num=0; num<4; num++){
		
				ptr_fine[0][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[0][num] += 0.25*ptr_coarse[1][num];

      				ptr_fine[1][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[1][num] += 0.25*ptr_coarse[2][num];

      				ptr_fine[2][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[2][num] += 0.25*ptr_coarse[3][num];

      				ptr_fine[3][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[3][num] += 0.25*ptr_coarse[4][num];

      				ptr_fine[4][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[4][num] += 0.25*ptr_coarse[5][num];

      				ptr_fine[5][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[5][num] += 0.25*ptr_coarse[6][num];

      				ptr_fine[6][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[6][num] += 0.25*ptr_coarse[7][num];

      				ptr_fine[7][num] += 0.75*ptr_coarse[0][num];
      				ptr_fine[7][num] += 0.25*ptr_coarse[8][num];

				}
	}

}
*/



