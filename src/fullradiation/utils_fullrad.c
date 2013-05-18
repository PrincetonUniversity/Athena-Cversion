#include "../copyright.h"
/*==============================================================================
 * FILE: utils_fullrad.c
 *
 * PURPOSE: contains misc. functions require for computation of rad. transfer
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   get_weights_linear()     - 
 *   get_weights_parabolic()  -
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"


#ifdef FULL_RADIATION_TRANSFER

void UpdateRT(DomainS *pD){

  RadGridS *pRG=(pD->RadGrid);
  int i,j,k ,ifr, l, n, m;

  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;

  Real wimu;
  Real mu[3]; /* cosins with respect to three axis */
  Real mu2[6]; 	/* products of two angles, used for radiation pressure */


for(ifr=0; ifr<pRG->nf; ifr++){
  for(k=ks; k<=ke; k++){
     for(j=js; j<=je; j++){	
	for(i=is; i<=ie; i++){	  
		/* set the momentums to be zero */
		pRG->R[ifr][k][j][i].J = 0.0;
		for(m=0; m<3; m++)
			pRG->R[ifr][k][j][i].H[m] = 0.0;

		for(m=0; m<6; m++)
			pRG->R[ifr][k][j][i].K[m] = 0.0;

		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){

				/* sum rays along different directions */
				wimu = pRG->imu[ifr][l][n][k][j][i] * pRG->wmu[n][k][j][i];
				for(m=0; m<3; m++)
					mu[m] = pRG->mu[l][n][k][j][i][m];
				/* for energy density */
				pRG->R[ifr][k][j][i].J += wimu;
	
				/* for flux */		
				for(m=0; m<3; m++)
					pRG->R[ifr][k][j][i].H[m] += mu[m] * wimu;

				/* for radiation pressure */
				mu2[0] = mu[0] * mu[0];
				mu2[1] = mu[0] * mu[1];
				mu2[2] = mu[1] * mu[1];
				mu2[3] = mu[0] * mu[2];
				mu2[4] = mu[1] * mu[2];
				mu2[5] = mu[2] * mu[2];

				for(m=0; m<6; m++)
					pRG->R[ifr][k][j][i].K[m] += mu2[m] * wimu;
			}/* end nang */
		}/* end octant l */
	}/* end i */
     }/* end j */
  }/* end k */
} /* end frequency */


	return;
}




#endif /* RADIATION_TRANSFER */


