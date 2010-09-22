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
 *   integrate_1d_radMHD()
 *   integrate_init_1d()
 *   integrate_destruct_1d()
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

#if defined(radMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */




/* Radiation matter coupled speed */
static Real *Cspeeds;
/* The matrix coefficient */
static Real **Euler = NULL;
/* Right hand side of the Matrix equation */
static Real *RHSEuler = NULL;
/* Used to store the results from the Matrix solver */
static int *Ern1 = NULL;
static Real *Ern2 = NULL;

static Cons1DS *U1d=NULL;



/********Public function****************/
/*-------BackEuler(): Use back euler method to update E_r and Fluxr-----------*/

/*---------1D right now. To be improved later----------*/


void BackEuler(MeshS *pM)
{
/* Right now, only work for one domain. Modified later for MPI */


  	GridS *pG=(pM->Domain[0][0].Grid);
	Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
	Real dt=pG->dt, dx=pG->dx1;
	int il,iu, is = pG->is, ie = pG->ie;
  	int i, j, m, n;
	int js = pG->js;
	int ks = pG->ks;
	int Nmatrix;
	

	Real SEE, SErho, SEm;
	Real temperature, velocity, pressure, Sigmas;

	Real temp1, temp2;
  	Real theta[6];
  	Real phi[6];
  	Real Sigma_s = Sigma_t - Sigma_a;

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
   	Nmatrix = ie - is + 1 ;
	rad_hydro_init_1d(Nmatrix,pM);


 	 il = is - 1;
  	 iu = ie + 1;

		


/* In principle, should load a routine to calculate the tensor f */

/* Temperatory variables used to calculate the Matrix  */
  	
  	

/* Load 1D vector of conserved variables and calculate the source term */
 	for (i=is-nghost; i<=ie+nghost; i++) {
    		U1d[i].d  = pG->U[ks][js][i].d;
    		U1d[i].Mx = pG->U[ks][js][i].M1;
    		U1d[i].My = pG->U[ks][js][i].M2;
    		U1d[i].Mz = pG->U[ks][js][i].M3;
    		U1d[i].E  = pG->U[ks][js][i].E;
    		U1d[i].Er  = pG->U[ks][js][i].Er;
    		U1d[i].Fluxr1  = pG->U[ks][js][i].Fluxr1;
    		U1d[i].Fluxr2  = pG->U[ks][js][i].Fluxr2;
    		U1d[i].Fluxr3  = pG->U[ks][js][i].Fluxr3;
		U1d[i].Edd_11  = pG->U[ks][js][i].Edd_11;
	}
	

	
	
/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */


/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1 =size1, otherwise it is wrong */


	for(i=is; i<=ie+1; i++){
	 	Cspeeds[i-is] = (U1d[i].Edd_11 - U1d[i-1].Edd_11) 
				/ (U1d[i].Edd_11 + U1d[i-1].Edd_11); 		
	}


	for(i=is; i<=ie; i++){
/* E is the total energy. should subtract the kinetic energy and magnetic energy density */
    		pressure = (U1d[i].E - 0.5 * U1d[i].Mx * U1d[i].Mx / U1d[i].d )
			* (Gamma - 1);
/* if MHD - 0.5 * Bx * Bx   */

    		temperature = pressure / (U1d[i].d * R_ideal);
      

    		RHSEuler[2*(i-is)]   = U1d[i].Er + Crat * pG->dt * Sigma_a 
				* temperature * temperature * temperature * temperature;
    		RHSEuler[2*(i-is)+1] = U1d[i].Fluxr1 + pG->dt *  Sigma_a
				* temperature * temperature * temperature * temperature * U1d[i].Mx / U1d[i].d;
	}

		

/*--------------------Note--------------------*.
/* Should judge the boundary condition. Here just use the perodic boundary condition first */

/* Step 1b: Setup the Matrix */
		
 	/* First, set the common elements for different boundary conditions */
	for(i=is; i<=ie; i++){
		velocity = U1d[i].Mx / U1d[i].d; 
		theta[0] = -Crat * hdtodx1 * (1.0 + Cspeeds[i-is]) * U1d[i-1].Edd_11;
		theta[1] = -Crat * hdtodx1 * (1.0 + Cspeeds[i-is]);
		theta[2] = 1.0 + Crat * hdtodx1 * (1.0 + Cspeeds[i-is+1]) * U1d[i].Edd_11 
			+ Crat * hdtodx1 * (1.0 - Cspeeds[i-is]) * U1d[i].Edd_11 + Crat * pG->dt * Sigma_a 
			+ pG->dt * (Sigma_a - Sigma_s) * (1.0 + U1d[i].Edd_11) * velocity * velocity / Crat;
		theta[3] = Crat * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) 
				- pG->dt * (Sigma_a - Sigma_s) * velocity;
		theta[4] = -Crat * hdtodx1 * (1.0 - Cspeeds[i-is+1]) * U1d[i+1].Edd_11;
		theta[5] = Crat * hdtodx1 * (1.0 - Cspeeds[i-is+1]);

		phi[0]	= theta[0];
		phi[1]	= theta[1] * U1d[i-1].Edd_11;
		phi[2]	= Crat * hdtodx1 * (Cspeeds[i-is] + Cspeeds[i-is+1]) * U1d[i].Edd_11 
				- pG->dt * Sigma_t * (1.0 + U1d[i].Edd_11) * velocity + pG->dt * Sigma_a * velocity;
		phi[3]	= 1.0 + Crat * hdtodx1 * (2.0 + Cspeeds[i-is+1] 
		- Cspeeds[i-is]) * U1d[i].Edd_11 + Crat * pG->dt * Sigma_t;
		phi[4]	= -theta[4];
		phi[5]	= -theta[5] * U1d[i+1].Edd_11;
	

		if(i == is) {
			if((ix1 == 4) && (ox1 == 4)) {
				Euler[0][0] = theta[2];
				Euler[0][1] = theta[3];
				Euler[0][2] = theta[4];
				Euler[0][3] = theta[5];
				Euler[0][2*Nmatrix-2] = theta[0];
				Euler[0][2*Nmatrix-1] = theta[1];
				Euler[1][0] = phi[2];
				Euler[1][1] = phi[3];
				Euler[1][2] = phi[4];
				Euler[1][3] = phi[5];
				Euler[1][2*Nmatrix-2] = phi[0];
				Euler[1][2*Nmatrix-1] = phi[1];
				}
		/* Periodic boundary condition */
			else if((ix1 == 2) && (ox2 == 2)){
				Euler[0][0] = theta[2] + theta[0];
				Euler[0][1] = theta[3] + theta[1];
				Euler[0][2] = theta[4];
				Euler[0][3] = theta[5];				
				Euler[1][0] = phi[2] + phi[0];
				Euler[1][1] = phi[3] + phi[1];
				Euler[1][2] = phi[4];
				Euler[1][3] = phi[5];
				}
		/* outflow boundary condition */
			else 
			goto on_error;

		}//end if
		else if(i == ie){
			if((ix1 == 4) && (ox1 == 4)) {
				Euler[2*Nmatrix-2][0] = theta[4];
				Euler[2*Nmatrix-2][1] = theta[5];
				Euler[2*Nmatrix-2][2*Nmatrix-4] = theta[0];
				Euler[2*Nmatrix-2][2*Nmatrix-3] = theta[1];
				Euler[2*Nmatrix-2][2*Nmatrix-2] = theta[2];
				Euler[2*Nmatrix-2][2*Nmatrix-1] = theta[3];

				Euler[2*Nmatrix-1][0] = phi[4];
				Euler[2*Nmatrix-1][1] = phi[5];
				Euler[2*Nmatrix-1][2*Nmatrix-4] = phi[0];
				Euler[2*Nmatrix-1][2*Nmatrix-3] = phi[1];
				Euler[2*Nmatrix-1][2*Nmatrix-2] = phi[2];
				Euler[2*Nmatrix-1][2*Nmatrix-1] = phi[3];
			}
			/* Periodic boundary condition */
			else if((ix1 == 2) && (ox2 == 2)){
				Euler[2*Nmatrix-2][2*Nmatrix-4] = theta[0];
				Euler[2*Nmatrix-2][2*Nmatrix-3] = theta[1];
				Euler[2*Nmatrix-2][2*Nmatrix-2] = theta[2] + theta[4];
				Euler[2*Nmatrix-2][2*Nmatrix-1] = theta[3] + theta[5];

				
				Euler[2*Nmatrix-1][2*Nmatrix-4] = phi[0];
				Euler[2*Nmatrix-1][2*Nmatrix-3] = phi[1];
				Euler[2*Nmatrix-1][2*Nmatrix-2] = phi[2] + phi[4];
				Euler[2*Nmatrix-1][2*Nmatrix-1] = phi[3] + phi[5];
			}
			/* outflow boundary condition */
			else 
			goto on_error;
			
		}// end else if
		else {
			for(j=0; j<6; j++){
			Euler[2*(i-is)][j+2*(i-is-1)] = theta[j];
			Euler[2*(i-is)+1][j+2*(i-is-1)] = phi[j];
			}// end for j loop
		}// end else
	}// end for i loop


		
	
	
/* Step 1c: Solve the Matrix equation */
		ludcmp(Euler,2*Nmatrix,Ern1,Ern2);
		
		lubksb(Euler,2*Nmatrix,Ern1,RHSEuler);
		
		
	for(i=is;i<=ie;i++){
		pG->U[ks][js][i].Er	= RHSEuler[2*(i-is)];
		pG->U[ks][js][i].Fluxr1	= RHSEuler[2*(i-is)+1];
		U1d[i].Er		= RHSEuler[2*(i-is)];
		U1d[i].Fluxr1		= RHSEuler[2*(i-is)+1];
	}
	/* May need to update Edd_11 */


/* Update the ghost zones for Periodic boundary condition to be used later */
/* Boundary condition for pG grids is applied after this loop is finished */

	 /* Inner boundary condition */
		if(ix1 == 4) {
		pG->U[ks][js][is-1].Er = U1d[ie].Er;
		pG->U[ks][js][is-1].Fluxr1 = U1d[ie].Fluxr1;
		pG->U[ks][js][is-2].Er = U1d[ie-1].Er;
		pG->U[ks][js][is-2].Fluxr1 = U1d[ie-1].Fluxr1;
		}
		else if(ix1 == 2) {
		pG->U[ks][js][is-1].Er = U1d[is].Er;
		pG->U[ks][js][is-1].Fluxr1 = U1d[is].Fluxr1;
		pG->U[ks][js][is-2].Er = U1d[is].Er;
		pG->U[ks][js][is-2].Fluxr1 = U1d[is].Fluxr1;
		}
		else 
		goto on_error;

	/* Outer boundary condition */
		if(ox1 == 4) {
		pG->U[ks][js][ie+1].Er = U1d[is].Er;
		pG->U[ks][js][ie+1].Fluxr1 = U1d[is].Fluxr1;
		pG->U[ks][js][ie+2].Er = U1d[is+1].Er;
		pG->U[ks][js][ie+2].Fluxr1 = U1d[is+1].Fluxr1;
		}
		else if(ox1 == 2) {
		pG->U[ks][js][ie+1].Er = U1d[ie].Er;
		pG->U[ks][js][ie+1].Fluxr1 = U1d[ie].Fluxr1;
		pG->U[ks][js][ie+2].Er = U1d[ie].Er;
		pG->U[ks][js][ie+2].Fluxr1 = U1d[ie].Fluxr1;
		}
		else
		goto on_error;
		


/*-----------Finish---------------------*/

  
	/* Free the temporary variables just used for this grid calculation*/
	rad_hydro_destruct_1d(Nmatrix);
	
  return;	


	on_error:
	
	rad_hydro_destruct_1d(Nmatrix);
	ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}




/*-------------------------------------------------------------------------*/
/* rad_hydro_init_1d: function to allocate memory used just for radiation variables */
/* rad_hydro_destruct_1d(): function to free memory */
void rad_hydro_init_1d(int Ngrids, MeshS *pM)
{

	int size1=0,nl,nd, i, j;

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
	if ((RHSEuler = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;
	if ((Ern1 = (int*)malloc(2*Ngrids*sizeof(int))) == NULL) goto on_error;
	if ((Ern2 = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;

	if ((Euler = (Real**)malloc(2*Ngrids*sizeof(Real*))) == NULL) goto on_error;
	for(i=0; i<2*Ngrids; i++)
		if ((Euler[i] = (Real*)malloc(2*Ngrids*sizeof(Real))) == NULL) goto on_error;

	 /* Initialize the matrix */
	for(i=0; i<2*Ngrids; i++)
		for(j=0; j<2*Ngrids; j++)
			Euler[i][j] = 0.0;
	
	return;

	on_error:
    	rad_hydro_destruct_1d(Ngrids);
	ath_error("[BackEuler]: malloc returned a NULL pointer\n");
}


void rad_hydro_destruct_1d(int Ngrids)
{
	int i;
	if (Cspeeds 	!= NULL) free(Cspeeds);
	if (RHSEuler	!= NULL) free(RHSEuler);
	if (Ern1	!= NULL) free(Ern1);
	if (Ern2	!= NULL) free(Ern2);
	if (U1d 	!= NULL) free(U1d);
	
	for(i=0; i<2*Ngrids; i++)
		if (Euler[i] != NULL) free(Euler[i]);

	if (Euler != NULL) free(Euler);
}


#endif /* radMHD_INTEGRATOR */
