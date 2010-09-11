#include "../copyright.h"
/*==============================================================================
 * FILE: lr_states_plm_radMHD.c
 *
 * PURPOSE: Second order (piecewise linear) spatial reconstruction in the
 *   conserved variables. With the CTU integrator, a time-evolution
 *   (characteristic tracing) step is used to interpolate interface values
 *   to the half time level {n+1/2}.
 *  
 *   But we do have to calculate the eigenvalues and eigenvectors
 *
 *   Limiting is performed in the conserved (rather than characteristic)
 *   variables.  
 *
 * NOTATION: 
 *   W_{L,i-1/2} is reconstructed value on the left-side of interface at i-1/2
 *   W_{R,i-1/2} is reconstructed value on the right-side of interface at i-1/2
 *
 *   The L- and R-states at the left-interface in each cell are indexed i.
 *   W_{L,i-1/2} is denoted by Wl[i  ];   W_{R,i-1/2} is denoted by Wr[i  ]
 *   W_{L,i+1/2} is denoted by Wl[i+1];   W_{R,i+1/2} is denoted by Wr[i+1]
 *
 *   Internally, in this routine, Wlv and Wrv are the reconstructed values on
 *   the left-and right-side of cell center.  Thus (see Step 8),
 *     W_{L,i-1/2} = Wrv(i-1);  W_{R,i-1/2} = Wlv(i)
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   lr_states()          - computes L/R states
 *   lr_states_init()     - initializes memory for static global arrays
 *   lr_states_destruct() - frees memory for static global arrays
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef RADIATION


static Real **pW=NULL;


/*----------------------------------------------------------------------------*/
/* lr_states:
 * Input Arguments:
 *   W = conserved variables at cell centers along 1-D slice
 *   Bxc = B in direction of slice at cell center
 *   dtodx = dt/dx
 *   il,iu = lower and upper indices of zone centers in slice
 *   il = is; iu = ie;
 *   NWAVE includes the radiation Er and Fluxr
 * W and Bxc must be initialized over [il-2:iu+2]
 *
 * Output Arguments:
 *   Wl,Wr = L/R-states of Conserved variables at interfaces over [il:iu+1]
 */

void lr_states_cons(const GridS *pG, const Cons1DS W[], const Real Bxc[], 
               const Real dt, const Real dx, const int il, const int iu, 
               Cons1DS Wl[], Cons1DS Wr[], const int dir)
{
	int i,n,m;
	Real lim_slope1,lim_slope2,qa,qx;
	Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
	Real Aeff[NWAVE][NWAVE], Propa[NWAVE][NWAVE], Source[NWAVE];
	Real dWc[NWAVE+NSCALARS],dWl[NWAVE+NSCALARS];
	Real dWr[NWAVE+NSCALARS],dWg[NWAVE+NSCALARS];
	Real dac[NWAVE+NSCALARS],dal[NWAVE+NSCALARS];
	Real dar[NWAVE+NSCALARS],dag[NWAVE+NSCALARS],da[NWAVE+NSCALARS];
	Real Wlv[NWAVE+NSCALARS],Wrv[NWAVE+NSCALARS];
	Real dW[NWAVE+NSCALARS],dWm[NWAVE+NSCALARS];
	Real *pWl, *pWr;
	Real qx1,qx2,C;

/***Temporary variables used for each cell */
	Real temperature, TEnergy, aeff, enthalpy, SEE, Alpha, velocity, pressure;
	Real Tempvari;
	const Real dtodx = dt / dx;

 

/* Zero eigenmatrices, set pointer to primitive variables */
	for (n=0; n<NWAVE; n++) {
		for (m=0; m<NWAVE; m++) {
			rem[n][m] = 0.0;
			lem[n][m] = 0.0;
			Aeff[n][m] = 0.0;
			Propa[n][m] = 0.0;
			if(n==m) Propa[n][m] = 1.0;
		}
	}

	for (n=0; n<NWAVE; n++) Source[n] = 0.0;


	for (i=il-2; i<=iu+2; i++) pW[i] = (Real*)&(W[i]);

/*========================== START BIG LOOP OVER i =======================*/
  for (i=il-1; i<=iu+1; i++) {

/*--- Step 1. ------------------------------------------------------------------
 * Compute eigensystem and eigenvalues for effective matrix  */

	TEnergy= W[i].E;
	pressure = (TEnergy - 0.5 * W[i].Mx * W[i].Mx / W[i].d )
			* (Gamma - 1);
	/* Should include magnetic energy for MHD */
	temperature = pressure / (W[i].d * Ridealgas);
	velocity = W[i].Mx / W[i].d;


	enthalpy = Gamma * TEnergy / W[i].d - (Gamma - 1.0) * velocity * velocity / 2.0;

	SEE = 4.0 * Sigmaa * temperature * temperature * temperature * (Gamma - 1.0)/ W[i].d;
	Alpha = (1.0 - exp(-Pratio * Cratio * SEE * dt/2.0))/(Pratio * Cratio* SEE * dt/2.0);
	aeff = -(Gamma - 1.0) * velocity * velocity/2.0 + Alpha * (Gamma - 1.0) * enthalpy 
			+ (1.0 - Alpha) * (temperature + (Gamma - 1.0) * velocity * velocity/2.0);
	aeff = sqrt(aeff);
	
	esys_rad_hyd(aeff, W[i].Mx/W[i].d, ev, rem, lem);
	
	/* Setup effective matrix Aeff, may be modified if we include radiation here */
	Aeff[0][1] = 1;
	Aeff[1][0] = (Gamma - 3.0) * velocity * velocity /2.0;
	Aeff[1][1] = -(Gamma - 3.0) * velocity;
	Aeff[1][4] = (Gamma - 1.0);
	Aeff[4][0] = velocity * ((Gamma - 1.0) * velocity * velocity/2.0 - Alpha * enthalpy 
			- (1.0 - Alpha) * (temperature / (Gamma - 1.0) + 0.5 * velocity * velocity));
	Aeff[4][1] = -(Gamma - 1.0) * velocity * velocity + Alpha * enthalpy 
			+ (1 - Alpha) * (temperature / (Gamma - 1.0) + 0.5 * velocity * velocity);
	Aeff[4][4] = Gamma * velocity;

	Propa[4][0] = (1.0 - Alpha) * (TEnergy / W[i].d - velocity * velocity);
	Propa[4][1] = (1.0 - Alpha) * velocity;
	Propa[4][4] = Alpha;

	Source[1] = -Pratio * (-Sigmat * (W[i].Fluxr1 - (1.0 + pG->fra1D) * velocity * W[i].Er / Cratio)
	+ Sigmaa * velocity * (temperature * temperature * temperature * temperature - W[i].Er)/Cratio);
	Source[4] = -Pratio * Cratio * (Sigmaa * (temperature * temperature * temperature * temperature - W[i].Er)
		+(Sigmaa - (Sigmat - Sigmaa)) * velocity
		* (W[i].Fluxr1 - (1.0 + pG->fra1D) * velocity * W[i].Er / Cratio)/Cratio);

/*--- Step 2. ------------------------------------------------------------------
 * Compute centered, L/R, and van Leer differences of primitive variables
 * Note we access contiguous array elements by indexing pointers for speed */

/* First try just the conserved variable. May use equation (51) in Miniati & Colella later */


    	for (n=0; n<(NWAVE+NSCALARS); n++) {
     		 dWc[n] = pW[i+1][n] - pW[i-1][n];
     		 dWl[n] = pW[i][n]   - pW[i-1][n];
     		 dWr[n] = pW[i+1][n] - pW[i][n];

      		if (dWl[n]*dWr[n] > 0.0) {
        		dWg[n] = 2.0*dWl[n]*dWr[n]/(dWl[n]+dWr[n]);
      		} else {
        		dWg[n] = 0.0;
      		}
    	}


/*--- Step 3. ------------------------------------------------------------------
 * Apply monotonicity constraints to characteristic projections */

    	for (n=0; n<(NWAVE+NSCALARS); n++) {
      		da[n] = 0.0;
     	 if (dWl[n]*dWr[n] > 0.0) {
        	lim_slope1 = MIN(    fabs(dWl[n]),fabs(dWr[n]));
        	lim_slope2 = MIN(0.5*fabs(dWc[n]),fabs(dWg[n]));
        	dWm[n] = SIGN(dWc[n])*MIN(2.0*lim_slope1,lim_slope2);
      		}
    	}



/*--- Step 4. ------------------------------------------------------------------
 * Compute L/R values, ensure they lie between neighboring cell-centered vals */

/* Note Wlv is left side with respect to W and Wrv is right side with respect to W. */

    	for (n=0; n<(NWAVE+NSCALARS); n++) {
		Wlv[n] = 0.0;
		Wrv[n] = 0.0;
		for (m=0; m<(NWAVE+NSCALARS); n++) {
			Wlv[n] +=  dtodx * Aeff[n][m] * dWm[m];
			Wrv[n] += dtodx * Aeff[n][m] * dWm[m];	
		}
		Wlv[n] = -0.5 * (dWm[n] + Wlv[n]);
		Wrv[n] = 0.5 * (dWm[n] - Wrv[n]);
		Wlv[n] += pW[i][n];
		Wrv[n] += pW[i][n];
    	}

/*----Step 5-----------------------------------
 * Add source terms */

	pWl = (Real *) &(Wl[i+1]);
   	pWr = (Real *) &(Wr[i]);


	
	for (n=0; n<(NWAVE+NSCALARS); n++){
	Tempvari = 0.0;
		for(m=0; m<(NWAVE+NSCALARS); m++){
		Tempvari += 0.5 * dt * Propa[n][m] * Source[m];
	}
		Wlv[n] += Tempvari;
		Wrv[n] += Tempvari;
	}

	
/*------Step 6--------------------------
 * change Wlv and Wrv to Wl and Wr  	  */
	for(n=0; n<(NWAVE+NSCALARS); n++){
		pWl[n] = Wrv[n];
		pWr[n] = Wlv[n];
	}

	

} /*===================== END BIG LOOP OVER i ===========================*/


  return;
}


/*----------------------------------------------------------------------------*/
/* lr_states_init:  Allocate enough memory for work arrays */

void lr_states_init(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

  if ((pW = (Real**)malloc(nmax*sizeof(Real*))) == NULL) goto on_error;

  return;
  on_error:
    lr_states_destruct();
    ath_error("[lr_states_init]: malloc returned a NULL pointer\n");
}


/*----------------------------------------------------------------------------*/
/* lr_states_destruct:  Free memory used by work arrays */

void lr_states_destruct(void)
{
  if (pW != NULL) free(pW);
  return;
}


#endif /* RADIATION */
