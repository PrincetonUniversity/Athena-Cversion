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
 *   BackEuler()
 *
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

#ifdef MATRIX_LIS

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
/*================================*/
/* For the matrix solver */
/* we use lis library now */
#include <lis.h>

/*===============================*/

#endif


/* Radiation matter coupled speed */
static Real *Cspeeds;
/* The matrix coefficient */
static Real **Euler = NULL;

/* Used for initial guess solution, also return the correct solution */
/* Used for restarted GMRES method for periodic boundary condition */

static Real **lEuler = NULL;
/* Used to store the lower trianglular system in LU decompsion */

/* Right hand side of the Matrix equation */
static Real *RHSEuler = NULL;



/* Used to store the results from the Matrix solver */
static unsigned long *Ern1 = NULL;
static int *Ern1p = NULL;
static Real *Ern2 = NULL;

static Cons1DS *U1d=NULL;

/*========= variables used for periodic boundary condition =========*/
/* We use Lis library for periodic boundary condtion, similar to 2D case  */

/*Store the index of all non-zero elements */
/* The matrix coefficient, this is oneD in GMRES solver */
static LIS_MATRIX Eulerp;

static LIS_SOLVER solver;
/* Right hand side of the Matrix equation */
static LIS_VECTOR RHSEulerp;
/* RHSEuler is 0-------3*Nx*Ny-1; only for GMRES method */ 
static LIS_VECTOR INIguess;
/* Used for initial guess solution, also return the correct solution */

static LIS_SCALAR *Value;
static int *indexValue;
static int *ptr;



/********Public function****************/
/*-------BackEuler(): Use back euler method to update E_r and Fluxr-----------*/

/*---------1D right now. To be improved later----------*/


void BackEuler_1d(MeshS *pM)
{
/* Right now, only work for one domain. Modified later for SMR */


	
	DomainS *pD = &(pM->Domain[0][0]);
	GridS *pG = pD->Grid;
	Real hdtodx1 = 0.5*pG->dt/pG->dx1;
	Real dt = pG->dt;
	int il,iu, is = pG->is, ie = pG->ie;
  	int i, j, m;
	int js = pG->js;
	int ks = pG->ks;
	int Nmatrix, NZ_NUM, NoEr, NoFr, lines, count;
	int index, Matrixiter;
	Real tempvalue;
	
	Real tempEr;
	Real tempFr;
	Real temp0;
	
	Real temperature, velocity, Fr0x;
	Real AdvFx;
	
  	Real theta[7];
  	Real phi[7];
  	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF, pressure, density, T4;
	Real Sigma[NOPACITY];

	/* Boundary condition flag */
	int ix1, ox1, ix2, ox2, ix3, ox3;
	ix1 = pM->BCFlag_ix1;
	ox1 = pM->BCFlag_ox1;
	ix2 = pM->BCFlag_ix2;
	ox2 = pM->BCFlag_ox2;
	ix3 = pM->BCFlag_ix3;
	ox3 = pM->BCFlag_ox3;
	


/* Allocate memory space for the Matrix calculation, just used for this grids */
/* Nmatrix is the number of active cells just in this grids */
/* Matrix size should be 2 * Nmatrix, ghost zones are included*/

/* We keep similar form as 2D case , not all variables are necessary */
   	Nmatrix = ie - is + 1 ;
	NZ_NUM = 12 * Nmatrix;
	lines = 2 * Nmatrix;
	count = 0;
	NoEr = 0;
	NoFr = 0;

	


 	 il = is - 1;
  	 iu = ie + 1;


	int bgflag;		/* used to subtract whether subtract background or not */
	
	bgflag = 0;



/* In principle, should load a routine to calculate the tensor f */

/* Temperatory variables used to calculate the Matrix  */
  	

/* First, do the advection step and update bounary */
	


/* Load 1D vector of conserved variables and calculate the source term */
 	for (i=is-nghost; i<=ie+nghost; i++) {
    		U1d[i].d  = pG->U[ks][js][i].d;
    		U1d[i].Mx = pG->U[ks][js][i].M1;
    		U1d[i].My = pG->U[ks][js][i].M2;
    		U1d[i].Mz = pG->U[ks][js][i].M3;
    		U1d[i].E  = pG->U[ks][js][i].E;
    		U1d[i].Er  = pG->U[ks][js][i].Er;
    		U1d[i].Fr1  = pG->U[ks][js][i].Fr1;
    		U1d[i].Fr2  = pG->U[ks][js][i].Fr2;
    		U1d[i].Fr3  = pG->U[ks][js][i].Fr3;
		U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;		
	}
	
	
/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */

	/* calculate the guess temperature */
	/* Guess temperature is calculated in the main function */
/*	GetTguess(pM);
*/

/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1 =size1, otherwise it is wrong */



/* Step 1b: Setup the Matrix */
		
 	/* First, set the common elements for different boundary conditions */
	/* theta and phi are written according to the compact matrix form, not the order in the paper */
	for(i=is; i<=ie; i++){
		
		
		matrix_coef(NULL, pG, 1, i, js, ks, 0.0, &(theta[1]), &(phi[0]), NULL, NULL);
		theta[0] = 0.0;
		phi[6] = 0.0;


	/* non periodic boundary condition */
		/* Set the right hand side */
		T4 = pG->Tguess[ks][js][i];

		Sigma_sF = pG->U[ks][js][i].Sigma[0];
		Sigma_aF = pG->U[ks][js][i].Sigma[1];
		Sigma_aP = pG->U[ks][js][i].Sigma[2];
		Sigma_aE = pG->U[ks][js][i].Sigma[3];
	
		velocity = U1d[i].Mx / U1d[i].d;
	
		/* RHSEuler[0] is not used. RHSEuler[1...N]  */
		
		
		Rad_Advection_Flux1D(pD, i, js, ks, 1.0, &AdvFx);

    		RHSEuler[2*(i-is)+1]   = U1d[i].Er + dt * Sigma_aP * T4 * Crat * Eratio +  (1.0 - Eratio) * pG->Ersource[ks][js][i] + AdvFx;
    		RHSEuler[2*(i-is)+2] = U1d[i].Fr1 + Eratio * pG->dt *  Sigma_aP * T4 * velocity + (1.0 - Eratio) * pG->Ersource[ks][js][i] * velocity / Crat;


		
		/* For periodic boundary condition, we use Lis library */
		if((ix1 == 4) && (ox1 == 4)) {
				index = 2 * (i-is);
				tempvalue = RHSEuler[2*(i-is)+1];

				if(bgflag){
					tempEr = theta[1] * U1d[i-1].Er + theta[5] * U1d[i+1].Er;
					tempFr = theta[2] * U1d[i-1].Fr1 + theta[6] * U1d[i+1].Fr1;
					temp0 = theta[4] * U1d[i].Fr1;
					tempvalue -= theta[3] * U1d[i].Er + tempEr + tempFr + temp0;

				}

				lis_vector_set_value(LIS_INS_VALUE,index,tempvalue,RHSEulerp);

				++index;
				tempvalue = RHSEuler[2*(i-is)+2];


				if(bgflag){
					tempEr = phi[0] * U1d[i-1].Er + phi[4] * U1d[i+1].Er;
					tempFr = phi[1] * U1d[i-1].Fr1 + phi[5] * U1d[i+1].Fr1;
					temp0 = phi[3] * U1d[i].Fr1;
					tempvalue -= phi[2] * U1d[i].Er + tempEr + tempFr + temp0;

				}

				lis_vector_set_value(LIS_INS_VALUE,index,tempvalue,RHSEulerp);

				if(bgflag){
					lis_vector_set_value(LIS_INS_VALUE,2*(i-is),0.0,INIguess);
					lis_vector_set_value(LIS_INS_VALUE,2*(i-is)+1,0.0,INIguess);
				}
				else{
					lis_vector_set_value(LIS_INS_VALUE,2*(i-is),U1d[i].Er,INIguess);
					lis_vector_set_value(LIS_INS_VALUE,2*(i-is)+1,U1d[i].Fr1,INIguess);
				}
			
		

		}

		if(!bgflag){
		/* For inflow boundary condition */
		if((i == is) && (ix1 == 3)) {
			
			RHSEuler[1] -= (theta[1] * U1d[i-1].Er + theta[2] * U1d[i-1].Fr1);
			RHSEuler[2] -= (phi[0] * U1d[i-1].Er + phi[1] * U1d[i-1].Fr1);
		}

		if((i == ie) && (ox1 == 3)) {
			
			RHSEuler[2 * Nmatrix -1] -= (theta[5] * U1d[i+1].Er + theta[6] * U1d[i+1].Fr1);
			RHSEuler[2 * Nmatrix] -= (phi[4] * U1d[i+1].Er + phi[5] * U1d[i+1].Fr1);
		}
		}


	if((ix1 != 4) && (ox1 !=4)) {	
		for(j=1; j<=7; j++) {
			Euler[2*(i-is)+1][j] = theta[j-1];
			Euler[2*(i-is)+2][j] = phi[j-1];
		}	

		/* Judge the boundary condition */
		if(i == is) {
			/* clean the corner */
			Euler[1][1] = 0.0;
			Euler[1][2] = 0.0;
			Euler[1][3] = 0.0;
			Euler[2][1] = 0.0;
			Euler[2][2] = 0.0;
			
			if(ix1 == 2) {
				Euler[1][4] = theta[3] + theta[1];
				Euler[1][5] = theta[4] + theta[2];
				Euler[2][3] = phi[2] + phi[0];
				Euler[2][4] = phi[3] + phi[1];
			}/* outflow boundary condition */
			else if(ix1 == 1 || ix1 == 5) {
				Euler[1][4] = theta[3] + theta[1];
				Euler[1][5] = theta[4] - theta[2];
				Euler[2][3] = phi[2] + phi[0];
				Euler[2][4] = phi[3] - phi[1];
			}/* reflecting boundary condition and conducting boundary condition*/
			else if(ix1 == 3) ;/*inflow boundary condition, do nothing */
			else
			goto on_error;			
		}

		if(i == ie) {
			Euler[2*Nmatrix-1][6] = 0.0;
			Euler[2*Nmatrix-1][7] = 0.0;
			Euler[2*Nmatrix][5] = 0.0;
			Euler[2*Nmatrix][6] = 0.0;
			Euler[2*Nmatrix][7] = 0.0;
			
			if(ox1 == 2) {
				Euler[2*Nmatrix-1][4] = theta[3] + theta[5];
				Euler[2*Nmatrix-1][5] = theta[4] + theta[6];
				Euler[2*Nmatrix][3] = phi[2] + phi[4];
				Euler[2*Nmatrix][4] = phi[3] + phi[5];
			}/* outflow boundary condition */
			else if(ox1 == 1 || ox1 == 5) {
				Euler[2*Nmatrix-1][4] = theta[3] + theta[5];
				Euler[2*Nmatrix-1][5] = theta[4] - theta[6];
				Euler[2*Nmatrix][3] = phi[2] + phi[4];
				Euler[2*Nmatrix][4] = phi[3] - phi[5];
			}/* reflecting boundary condition and conducting bounary condition */
			else if(ox1 == 3) ;/* inflow boundary condition, do nothing*/
			else
			goto on_error;	
		}
	}
	/* End for non-periodic boundary condition */
	/* EulerLU is already initialized for zeros */
	if((ix1 == 4) && (ox1 == 4)) {

		NoEr = count;
		NoFr = count + 6;
		count += 12;

		ptr[2*(i-is)] = NoEr;
		ptr[2*(i-is)+1] = NoFr;
		
		if(i == is){			

			/* For Er */
			for(j=0; j<4; j++)
				Value[NoEr+j] = theta[3+j];

			/* For Fr */
			for(j=0; j<4; j++)
				Value[NoFr+j] = phi[2+j];

			/* For Er */
			for(j=0; j<4; j++)
				indexValue[NoEr+j] = j;

			/* For Fr */
			for(j=0; j<4; j++)
				indexValue[NoFr+j] = j;
			
			/* For Er */
			Value[NoEr + 4] = theta[1];
			Value[NoEr + 5] = theta[2];

			indexValue[NoEr + 4] = 2 * (ie - is);
			indexValue[NoEr + 5] = 2 * (ie - is) + 1;

			/* For Fr */

			Value[NoFr + 4] = phi[0];
			Value[NoFr + 5] = phi[1];

			indexValue[NoFr + 4] = 2 * (ie - is);
			indexValue[NoFr + 5] = 2 * (ie - is) + 1;
		} /* End i == is */
		else if(i == ie){
			/* For Er */
			Value[NoEr] = theta[5];
			Value[NoEr + 1] = theta[6];

			indexValue[NoEr] = 0;
			indexValue[NoEr + 1] = 1;

			/* For Fr */

			Value[NoFr] = phi[4];
			Value[NoFr + 1] = phi[5];

			indexValue[NoFr] = 0;
			indexValue[NoFr + 1] = 1;



			/* For Er */
			for(j=2; j<6; j++)
				Value[NoEr+j] = theta[j-1];

			/* For Fr */
			for(j=2; j<6; j++)
				Value[NoFr+j] = phi[j-2];

			/* For Er */
			for(j=2; j<6; j++)
				indexValue[NoEr+j] = 2 * (i - is - 1) + j - 2;

			/* For Fr */
			for(j=2; j<6; j++)
				indexValue[NoFr+j] = 2 * (i - is - 1) + j - 2;

		}/* End i == ie */
		else{
			/* For Er */
			for(j=0; j<6; j++)
				Value[NoEr+j] = theta[1+j];

			/* For Fr */
			for(j=0; j<6; j++)
				Value[NoFr+j] = phi[j];

			/* For Er */
			for(j=0; j<6; j++)
				indexValue[NoEr+j] = 2 * (i - is - 1) + j;

			/* For Fr */
			for(j=0; j<6; j++)
				indexValue[NoFr+j] = 2 * (i - is - 1) + j;

		}
		
	}
	/* End for periodic boundary condition */	
	}/* end for i loop */

	
	
	
/* Step 1c: Solve the Matrix equation for band diaagnoal system */
	if((ix1 != 4) && (ox1 !=4)){
		bandec(Euler, (unsigned long) (2*Nmatrix), 3, 3, lEuler, Ern1, Ern2);
		banbks(Euler, (unsigned long) (2*Nmatrix), 3, 3, lEuler, Ern1, RHSEuler);

	}
	else {
		ptr[lines] = NZ_NUM;
		/* Assemble the matrix and solve the matrix */
		lis_matrix_set_crs(NZ_NUM,ptr,indexValue,Value,Eulerp);
		lis_matrix_assemble(Eulerp);

		
		lis_solver_set_option("-i gmres -p ilu",solver);
		lis_solver_set_option("-tol 1.0e-12",solver);
		lis_solver_set_option("-maxiter 2000",solver);
		lis_solve(Eulerp,RHSEulerp,INIguess,solver);
		
		/* check the iteration step to make sure 1.0e-12 is reached */
		lis_solver_get_iters(solver,&Matrixiter);

		ath_pout(0,"Matrix Iteration steps: %d\n",Matrixiter);
	}
	
			
	for(i=is;i<=ie;i++){	
		temp0 = pG->U[ks][js][i].Er;	

		if((ix1==4)&&(ox1==4)){
			lis_vector_get_value(INIguess,2*(i-is),&(tempEr));
			lis_vector_get_value(INIguess,2*(i-is)+1,&(tempFr));
			if(bgflag){
				pG->U[ks][js][i].Er += tempEr;
				pG->U[ks][js][i].Fr1 += tempFr;				
			
			}
			else{
				pG->U[ks][js][i].Er = tempEr;
				pG->U[ks][js][i].Fr1 = tempFr;
				
			}			

			
		}
		else{
			pG->U[ks][js][i].Er	= RHSEuler[2*(i-is)+1];
			pG->U[ks][js][i].Fr1	= RHSEuler[2*(i-is)+2];
			U1d[i].Er		= RHSEuler[2*(i-is)+1];
			U1d[i].Fr1		= RHSEuler[2*(i-is)+2];
		}	

               	velocity = pG->U[ks][js][i].M1 /pG->U[ks][js][i].d;


		Fr0x =  pG->U[ks][js][i].Fr1 - (1.0 +  pG->U[ks][js][i].Edd_11) * velocity * pG->U[ks][js][i].Er / Crat; 
			
			/* Estimate the added energy source term */
		if(Prat > 0.0){
			
			pG->U[ks][js][i].Er += (pG->Eulersource[ks][js][i] - dt * (pG->U[ks][js][i].Sigma[1] -  pG->U[ks][js][i].Sigma[0]) * velocity * Fr0x);
			
		}			
			
		
	}
	/* May need to update Edd_11 */

/*
if(Opacity != NULL){
		for(i=il+1; i<=iu-1; i++) {

			pressure = (pG->U[ks][js][i].E - 0.5 * pG->U[ks][js][i].M1 * pG->U[ks][js][i].M1 / pG->U[ks][js][i].d )
				* (Gamma - 1);

#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[ks][js][i].B1c * pG->U[ks][js][i].B1c + pG->U[ks][js][i].B2c * pG->U[ks][js][i].B2c + pG->U[ks][js][i].B3c * pG->U[ks][js][i].B3c) * (Gamma - 1.0);
#endif
		
			temperature = pressure / (pG->U[ks][js][i].d * R_ideal);
	
		
			Opacity(pG->U[ks][js][i].d, temperature,Sigma, NULL);
			for(m=0;m<NOPACITY;m++){
				pG->U[ks][js][i].Sigma[m] = Sigma[m];
			}
	
		}
}	
*/

/* Update the ghost zones for different boundary condition to be used later */
/*
	for (i=0; i<pM->NLevels; i++){ 
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
        	if (pM->Domain[i][js].Grid != NULL){
  			bvals_radMHD(&(pM->Domain[i][js]));

        	}
      	     }
    	}
*/

/*-----------Finish---------------------*/


	
  return;	


	on_error:
	
	BackEuler_destruct_1d(Nmatrix);
	ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}




/*-------------------------------------------------------------------------*/
/* BackEuler_init_1d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_1d(): function to free memory */
void BackEuler_init_1d(MeshS *pM)
{
/* The matrix Euler is stored as a compact form. See $2.4 of numerical recipes */
/* Ngrids = ie - is + 1; */

	int size1=0,nl,nd, i, j;

	int Ngrids;

	Ngrids = pM->Domain[0][0].Grid->ie - pM->Domain[0][0].Grid->is + 1;

	/* Cycle over all Grids on this processor to find maximum Nx1 */
	  for (nl=0; nl<(pM->NLevels); nl++){
	    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
	      if (pM->Domain[nl][nd].Grid != NULL) {
	        if ((pM->Domain[nl][nd].Grid->Nx[0]) > size1){
        	  size1 = pM->Domain[nl][nd].Grid->Nx[0];
	        }
	      }
	    }
	  }

	size1 = size1 + 2*nghost;
	if ((U1d       =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;


	if ((Cspeeds = (Real*)malloc((Ngrids+1)*sizeof(Real))) == NULL) goto on_error;
	if ((RHSEuler = (Real*)malloc((2*Ngrids+1)*sizeof(Real))) == NULL) goto on_error;
	
	if ((Ern1 = (unsigned long *)malloc((2*Ngrids+1)*sizeof(unsigned long))) == NULL) goto on_error;
	if ((Ern1p = (int *)malloc((2*Ngrids+1)*sizeof(int))) == NULL) goto on_error;
	if ((Ern2 = (Real*)malloc((2*Ngrids+1)*sizeof(Real))) == NULL) goto on_error;

	if ((Euler = (Real**)malloc((2*Ngrids+1)*sizeof(Real*))) == NULL) goto on_error;
	
	if ((lEuler = (Real**)malloc((2*Ngrids+1)*sizeof(Real*))) == NULL) goto on_error;
	for(i=0; i<(2*Ngrids+1); i++) {
		if ((Euler[i] = (Real*)malloc(8*sizeof(Real))) == NULL) goto on_error;
		if ((lEuler[i] = (Real*)malloc(8*sizeof(Real))) == NULL) goto on_error;		
	}


	 /* Initialize the matrix */
	for(i=0; i<(2*Ngrids+1); i++)
		for(j=0; j<8; j++) {
			Euler[i][j] = 0.0;
			lEuler[i][j] = 0.0;
		}


/* ====== For Lis vectors and matrix ========*/
	int lines;
	lines = 2 * Ngrids;

	lis_matrix_create(0,&Eulerp);	
	lis_matrix_set_size(Eulerp,0,lines);

	lis_vector_duplicate(Eulerp,&RHSEulerp);
	lis_vector_duplicate(Eulerp,&INIguess);

	lis_solver_create(&solver);

	/* Allocate value and index array */
/*
	if ((Value = (Real*)malloc(Nelements*sizeof(Real))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");
*/
	Value = (LIS_SCALAR *)malloc(12*Ngrids*sizeof(LIS_SCALAR));

	if ((indexValue = (int*)malloc(12*Ngrids*sizeof(int))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	if ((ptr = (int*)malloc((lines+1)*sizeof(int))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	
	return;

	on_error:
    	BackEuler_destruct_1d(Ngrids);
	ath_error("[BackEuler]: malloc returned a NULL pointer\n");
}


void BackEuler_destruct_1d()
{
	int i;
	if (Cspeeds 	!= NULL) free(Cspeeds);
	if (RHSEuler	!= NULL) free(RHSEuler);
	if (Ern1	!= NULL) free(Ern1);
	if (Ern1p	!= NULL) free(Ern1p);
	if (Ern2	!= NULL) free(Ern2);
	if (U1d 	!= NULL) free(U1d);
	
	
	free_2d_array(Euler);
	free_2d_array(lEuler);

/* For Lis library */
	
	lis_matrix_destroy(Eulerp);
	lis_solver_destroy(solver);	
	lis_vector_destroy(RHSEulerp);
	lis_vector_destroy(INIguess);
	
/*	if(Value != NULL) free(Value);

	if(indexValue != NULL) free(indexValue);
	if(ptr != NULL) free(ptr);
*/
}




#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_LIS */
