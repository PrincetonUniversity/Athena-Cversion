
#include "copyright.h"
/*==============================================================================
 * FILE: bvals_mhd.c
 *
 * PURPOSE: Sets boundary conditions (quantities in ghost zones) for radiation 
 * energy density and radiation flux.  
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   bvals_rad()      - calls appropriate functions to set ghost cells
 *  
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* bvals_rad: set ghost zones for radiation quantities.
 * right now is only for one domain. To be extended later for multiple domains.
 */

/*----------------- Inflow boundary condition should be modified later */

#ifdef rad_hydro


void bvals_rad(MeshS *pM)
{
	GridS *pGrid=(pM->Domain[0][0].Grid);
	int is = pGrid->is, ie = pGrid->ie;
	int i, j, k;
	int js = pGrid->js;
	int ks = pGrid->ks;
	int je = pGrid->je;
	int ke = pGrid->ke;

	int ix1, ox1, ix2, ox2, ix3, ox3;
	ix1 = pM->BCFlag_ix1;
	ox1 = pM->BCFlag_ox1;
	ix2 = pM->BCFlag_ix2;
	ox2 = pM->BCFlag_ox2;
	ix3 = pM->BCFlag_ix3;
	ox3 = pM->BCFlag_ox3;
	
 

	/* Boundary condition for pG grids is applied after this loop is finished */

	 /* Inner boundary condition */
	if(ix1 == 1) {
		for (k=ks; k<=ke; k++) {
   			 for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er  =  pGrid->U[k][j][is+(i-1)].Er;
				pGrid->U[k][j][is-i].Fr1 = -pGrid->U[k][j][is+(i-1)].Fr1; /* reflect 1-flux. */
			      }
    			}
  		}
	}

	else if(ix1 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er  = pGrid->U[k][j][is].Er;
				pGrid->U[k][j][is-i].Fr1 = pGrid->U[k][j][is].Fr1;
      				}
    			}
  		}	
	}
	else if(ix1 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er = pGrid->U[k][j][ie-(i-1)].Er;
				pGrid->U[k][j][is-i].Fr1 = pGrid->U[k][j][ie-(i-1)].Fr1;
      				}
    			}
  		}	
	}
	else 
	goto on_error;

	/* Outer boundary condition */
	if(ox1 == 1) {

		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er    =  pGrid->U[k][j][ie-(i-1)].Er;
        			pGrid->U[k][j][ie+i].Fr1 = -pGrid->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */
      				}
    			}
  		}		
	}

	else if(ox1 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie].Er;
				pGrid->U[k][j][ie+i].Fr1 = pGrid->U[k][j][ie].Fr1;
      				}
    			}
  		}
	}
	else if(ox1 == 3) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er = 8.5 * 8.5 * 8.5 * 8.5;
				pGrid->U[k][j][ie+i].Fr1 = -30.0*0.333333*8.5*8.5*8.5/Sigma_t;
      				}
    			}
  		}
	}
	else if(ox1 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][is+(i-1)].Er;
				pGrid->U[k][j][ie+i].Fr1 = pGrid->U[k][j][is+(i-1)].Fr1;
      				}
    			}
  		}
	}
	else
	goto on_error;	

  	return;

		on_error:
		ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}

#endif

