
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

#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)


void bvals_radMHD(MeshS *pM)
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
	/* Boundary condition for x direction */

	 /* Inner boundary condition */
	/*reflecting boundary condition */
	if(pGrid->Nx[0] > 1) {
	if(ix1 == 1) {
		for (k=ks; k<=ke; k++) {
   			 for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er  =  pGrid->U[k][j][is+(i-1)].Er;
				pGrid->U[k][j][is-i].Fr2  =  pGrid->U[k][j][is+(i-1)].Fr2;
				pGrid->U[k][j][is-i].Fr1 = -pGrid->U[k][j][is+(i-1)].Fr1; /* reflect 1-flux. */
				pGrid->U[k][j][is-i].Sigma_a  =  pGrid->U[k][j][is+(i-1)].Sigma_a;
				pGrid->U[k][j][is-i].Sigma_t  =  pGrid->U[k][j][is+(i-1)].Sigma_t;
				pGrid->U[k][j][is-i].Edd_11  =  pGrid->U[k][j][is+(i-1)].Edd_11;
				pGrid->U[k][j][is-i].Edd_21  =  pGrid->U[k][j][is+(i-1)].Edd_21;
				pGrid->U[k][j][is-i].Edd_22  =  pGrid->U[k][j][is+(i-1)].Edd_22;
				pGrid->U[k][j][is-i].Edd_31  =  pGrid->U[k][j][is+(i-1)].Edd_31;
				pGrid->U[k][j][is-i].Edd_32  =  pGrid->U[k][j][is+(i-1)].Edd_32;
				pGrid->U[k][j][is-i].Edd_33  =  pGrid->U[k][j][is+(i-1)].Edd_33;
			      }
    			}
  		}
	}
	/* outflow boundary condition */
	else if(ix1 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er  = pGrid->U[k][j][is].Er;
				pGrid->U[k][j][is-i].Fr2  = pGrid->U[k][j][is].Fr2;
				pGrid->U[k][j][is-i].Fr1 = pGrid->U[k][j][is].Fr1;
				pGrid->U[k][j][is-i].Sigma_a  =  pGrid->U[k][j][is].Sigma_a;
				pGrid->U[k][j][is-i].Sigma_t  =  pGrid->U[k][j][is].Sigma_t;
				pGrid->U[k][j][is-i].Edd_11  =  pGrid->U[k][j][is].Edd_11;
				pGrid->U[k][j][is-i].Edd_21  =  pGrid->U[k][j][is].Edd_21;
				pGrid->U[k][j][is-i].Edd_22  =  pGrid->U[k][j][is].Edd_22;
				pGrid->U[k][j][is-i].Edd_31  =  pGrid->U[k][j][is].Edd_31;
				pGrid->U[k][j][is-i].Edd_32  =  pGrid->U[k][j][is].Edd_32;
				pGrid->U[k][j][is-i].Edd_33  =  pGrid->U[k][j][is].Edd_33;
      				}
    			}
  		}	
	}
	/* Inflow boundary condition */
	else if(ix1 == 3) {
		if(pM->Domain[0][0].rad_ix1_BCFun == NULL)
			goto on_error;
		else
			(*(pM->Domain[0][0].rad_ix1_BCFun))(pGrid);
	}
	/* periodic boundary condition */
	else if(ix1 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][is-i].Er = pGrid->U[k][j][ie-(i-1)].Er;
				pGrid->U[k][j][is-i].Fr2 = pGrid->U[k][j][ie-(i-1)].Fr2;
				pGrid->U[k][j][is-i].Fr1 = pGrid->U[k][j][ie-(i-1)].Fr1;
				pGrid->U[k][j][is-i].Sigma_a = pGrid->U[k][j][ie-(i-1)].Sigma_a;
				pGrid->U[k][j][is-i].Sigma_t = pGrid->U[k][j][ie-(i-1)].Sigma_t;
				pGrid->U[k][j][is-i].Edd_11  =  pGrid->U[k][j][ie-(i-1)].Edd_11;
				pGrid->U[k][j][is-i].Edd_21  =  pGrid->U[k][j][ie-(i-1)].Edd_21;
				pGrid->U[k][j][is-i].Edd_22  =  pGrid->U[k][j][ie-(i-1)].Edd_22;
				pGrid->U[k][j][is-i].Edd_31  =  pGrid->U[k][j][ie-(i-1)].Edd_31;
				pGrid->U[k][j][is-i].Edd_32  =  pGrid->U[k][j][ie-(i-1)].Edd_32;
				pGrid->U[k][j][is-i].Edd_33  =  pGrid->U[k][j][ie-(i-1)].Edd_33;
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
				pGrid->U[k][j][ie+i].Fr2    =  pGrid->U[k][j][ie-(i-1)].Fr2;
        			pGrid->U[k][j][ie+i].Fr1 = -pGrid->U[k][j][ie-(i-1)].Fr1; /* reflect 1-flux. */
				pGrid->U[k][j][ie+i].Sigma_a    =  pGrid->U[k][j][ie-(i-1)].Sigma_a;
				pGrid->U[k][j][ie+i].Sigma_t    =  pGrid->U[k][j][ie-(i-1)].Sigma_t;
				pGrid->U[k][j][ie+i].Edd_11  =  pGrid->U[k][j][ie-(i-1)].Edd_11;
				pGrid->U[k][j][ie+i].Edd_21  =  pGrid->U[k][j][ie-(i-1)].Edd_21;
				pGrid->U[k][j][ie+i].Edd_22  =  pGrid->U[k][j][ie-(i-1)].Edd_22;
				pGrid->U[k][j][ie+i].Edd_31  =  pGrid->U[k][j][ie-(i-1)].Edd_31;
				pGrid->U[k][j][ie+i].Edd_32  =  pGrid->U[k][j][ie-(i-1)].Edd_32;
				pGrid->U[k][j][ie+i].Edd_33  =  pGrid->U[k][j][ie-(i-1)].Edd_33;
      				}
    			}
  		}		
	}

	else if(ox1 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie].Er;
				pGrid->U[k][j][ie+i].Fr2 = pGrid->U[k][j][ie].Fr2;
				pGrid->U[k][j][ie+i].Fr1 = pGrid->U[k][j][ie].Fr1;
				pGrid->U[k][j][ie+i].Sigma_a = pGrid->U[k][j][ie].Sigma_a;
				pGrid->U[k][j][ie+i].Sigma_t = pGrid->U[k][j][ie].Sigma_t;
				pGrid->U[k][j][ie+i].Edd_11  =  pGrid->U[k][j][ie].Edd_11;
				pGrid->U[k][j][ie+i].Edd_21  =  pGrid->U[k][j][ie].Edd_21;
				pGrid->U[k][j][ie+i].Edd_22  =  pGrid->U[k][j][ie].Edd_22;
				pGrid->U[k][j][ie+i].Edd_31  =  pGrid->U[k][j][ie].Edd_31;
				pGrid->U[k][j][ie+i].Edd_32  =  pGrid->U[k][j][ie].Edd_32;
				pGrid->U[k][j][ie+i].Edd_33  =  pGrid->U[k][j][ie].Edd_33;
      				}
    			}
  		}
	}
	else if(ox1 == 3) {
		if(pM->Domain[0][0].rad_ox1_BCFun == NULL)
			goto on_error;
		else
			(*(pM->Domain[0][0].rad_ox1_BCFun))(pGrid);
	}
	else if(ox1 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {
        			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][is+(i-1)].Er;
				pGrid->U[k][j][ie+i].Fr2 = pGrid->U[k][j][is+(i-1)].Fr2;
				pGrid->U[k][j][ie+i].Fr1 = pGrid->U[k][j][is+(i-1)].Fr1;
				pGrid->U[k][j][ie+i].Sigma_a = pGrid->U[k][j][is+(i-1)].Sigma_a;
				pGrid->U[k][j][ie+i].Sigma_t = pGrid->U[k][j][is+(i-1)].Sigma_t;
				pGrid->U[k][j][ie+i].Edd_11  =  pGrid->U[k][j][is+(i-1)].Edd_11;
				pGrid->U[k][j][ie+i].Edd_21  =  pGrid->U[k][j][is+(i-1)].Edd_21;
				pGrid->U[k][j][ie+i].Edd_22  =  pGrid->U[k][j][is+(i-1)].Edd_22;
				pGrid->U[k][j][ie+i].Edd_31  =  pGrid->U[k][j][is+(i-1)].Edd_31;
				pGrid->U[k][j][ie+i].Edd_32  =  pGrid->U[k][j][is+(i-1)].Edd_32;
				pGrid->U[k][j][ie+i].Edd_33  =  pGrid->U[k][j][is+(i-1)].Edd_33;
      				}
    			}
  		}
	}
	else
	goto on_error;	
	}

	/* Boundary condition for y direction */
	/*------------------------------------------------------------*/
	/*reflecting boundary condition */
	if(pGrid->Nx[1] > 1){
	if(ix2 == 1) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][js-j][i].Er   =  pGrid->U[k][js+(j-1)][i].Er;
				pGrid->U[k][js-j][i].Fr2  = -pGrid->U[k][js+(j-1)][i].Fr2;
				pGrid->U[k][js-j][i].Fr1  =  pGrid->U[k][js+(j-1)][i].Fr1; /* reflect 1-flux. */
				pGrid->U[k][js-j][i].Sigma_a  = -pGrid->U[k][js+(j-1)][i].Sigma_a;
				pGrid->U[k][js-j][i].Sigma_t  = -pGrid->U[k][js+(j-1)][i].Sigma_t;
				pGrid->U[k][js-j][i].Edd_11  =  pGrid->U[k][js+(j-1)][i].Edd_11;
				pGrid->U[k][js-j][i].Edd_21  =  pGrid->U[k][js+(j-1)][i].Edd_21;
				pGrid->U[k][js-j][i].Edd_22  =  pGrid->U[k][js+(j-1)][i].Edd_22;
				pGrid->U[k][js-j][i].Edd_31  =  pGrid->U[k][js+(j-1)][i].Edd_31;
				pGrid->U[k][js-j][i].Edd_32  =  pGrid->U[k][js+(j-1)][i].Edd_32;
				pGrid->U[k][js-j][i].Edd_33  =  pGrid->U[k][js+(j-1)][i].Edd_33;
			      }
    			}
  		}
	}
	/* outflow boundary condition */
	else if(ix2 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][js-j][i].Er  = pGrid->U[k][js][i].Er;
				pGrid->U[k][js-j][i].Fr2  = pGrid->U[k][js][i].Fr2;
				pGrid->U[k][js-j][i].Fr1 = pGrid->U[k][js][i].Fr1;
				pGrid->U[k][js-j][i].Sigma_a = pGrid->U[k][js][i].Sigma_a;
				pGrid->U[k][js-j][i].Sigma_t = pGrid->U[k][js][i].Sigma_t;
				pGrid->U[k][js-j][i].Edd_11  =  pGrid->U[k][js][i].Edd_11;
				pGrid->U[k][js-j][i].Edd_21  =  pGrid->U[k][js][i].Edd_21;
				pGrid->U[k][js-j][i].Edd_22  =  pGrid->U[k][js][i].Edd_22;
				pGrid->U[k][js-j][i].Edd_31  =  pGrid->U[k][js][i].Edd_31;
				pGrid->U[k][js-j][i].Edd_32  =  pGrid->U[k][js][i].Edd_32;
				pGrid->U[k][js-j][i].Edd_33  =  pGrid->U[k][js][i].Edd_33;
      				}
    			}
  		}	
	}
	/* inflow boundary condition */
	else if(ix2 == 3) {
		if(pM->Domain[0][0].rad_ix2_BCFun == NULL)
			goto on_error;
		else
			(*(pM->Domain[0][0].rad_ix2_BCFun))(pGrid);
	}
	/* periodic boundary condition */
	else if(ix2 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][js-j][i].Er  = pGrid->U[k][je-(j-1)][i].Er;
				pGrid->U[k][js-j][i].Fr1 = pGrid->U[k][je-(j-1)][i].Fr1;
				pGrid->U[k][js-j][i].Fr2 = pGrid->U[k][je-(j-1)][i].Fr2;
				pGrid->U[k][js-j][i].Sigma_a = pGrid->U[k][je-(j-1)][i].Sigma_a;
				pGrid->U[k][js-j][i].Sigma_t = pGrid->U[k][je-(j-1)][i].Sigma_t;
				pGrid->U[k][js-j][i].Edd_11  =  pGrid->U[k][je-(j-1)][i].Edd_11;
				pGrid->U[k][js-j][i].Edd_21  =  pGrid->U[k][je-(j-1)][i].Edd_21;
				pGrid->U[k][js-j][i].Edd_22  =  pGrid->U[k][je-(j-1)][i].Edd_22;
				pGrid->U[k][js-j][i].Edd_31  =  pGrid->U[k][je-(j-1)][i].Edd_31;
				pGrid->U[k][js-j][i].Edd_32  =  pGrid->U[k][je-(j-1)][i].Edd_32;
				pGrid->U[k][js-j][i].Edd_33  =  pGrid->U[k][je-(j-1)][i].Edd_33;
      				}
    			}
  		}	
	}
	else 
	goto on_error;

	/* Outer boundary condition */
	if(ox2 == 1) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][je+j][i].Er   =  pGrid->U[k][je-(j-1)][i].Er;
				pGrid->U[k][je+j][i].Fr2  = -pGrid->U[k][je-(j-1)][i].Fr2;
				pGrid->U[k][je+j][i].Fr1  =  pGrid->U[k][je-(j-1)][i].Fr1; /* reflect 1-flux. */
				pGrid->U[k][je+j][i].Sigma_a  = -pGrid->U[k][je-(j-1)][i].Sigma_a;
				pGrid->U[k][je+j][i].Sigma_t  = -pGrid->U[k][je-(j-1)][i].Sigma_t;
				pGrid->U[k][je+j][i].Edd_11  =  pGrid->U[k][je-(j-1)][i].Edd_11;
				pGrid->U[k][je+j][i].Edd_21  =  pGrid->U[k][je-(j-1)][i].Edd_21;
				pGrid->U[k][je+j][i].Edd_22  =  pGrid->U[k][je-(j-1)][i].Edd_22;
				pGrid->U[k][je+j][i].Edd_31  =  pGrid->U[k][je-(j-1)][i].Edd_31;
				pGrid->U[k][je+j][i].Edd_32  =  pGrid->U[k][je-(j-1)][i].Edd_32;
				pGrid->U[k][je+j][i].Edd_33  =  pGrid->U[k][je-(j-1)][i].Edd_33;
			      }
    			}
  		}
	}
	/* outflow boundary condition */
	else if(ox2 == 2) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][je+j][i].Er  = pGrid->U[k][je][i].Er;
				pGrid->U[k][je+j][i].Fr2  = pGrid->U[k][je][i].Fr2;
				pGrid->U[k][je+j][i].Fr1 = pGrid->U[k][je][i].Fr1;
				pGrid->U[k][je+j][i].Sigma_a = pGrid->U[k][je][i].Sigma_a;
				pGrid->U[k][je+j][i].Sigma_t = pGrid->U[k][je][i].Sigma_t;
				pGrid->U[k][je+j][i].Edd_11  =  pGrid->U[k][je][i].Edd_11;
				pGrid->U[k][je+j][i].Edd_21  =  pGrid->U[k][je][i].Edd_21;
				pGrid->U[k][je+j][i].Edd_22  =  pGrid->U[k][je][i].Edd_22;
				pGrid->U[k][je+j][i].Edd_31  =  pGrid->U[k][je][i].Edd_31;
				pGrid->U[k][je+j][i].Edd_32  =  pGrid->U[k][je][i].Edd_32;
				pGrid->U[k][je+j][i].Edd_33  =  pGrid->U[k][je][i].Edd_33;
      				}
    			}
  		}	
	}
	/* inflow boundary condition */
	else if(ox2 == 3) {
	if(pM->Domain[0][0].rad_ox2_BCFun == NULL)
			goto on_error;
		else
			(*(pM->Domain[0][0].rad_ox2_BCFun))(pGrid);
	}
	/* periodic boundary condition */
	else if(ox2 == 4) {
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++) {
        			pGrid->U[k][je+j][i].Er  = pGrid->U[k][js+(j-1)][i].Er;
				pGrid->U[k][je+j][i].Fr1 = pGrid->U[k][js+(j-1)][i].Fr1;
				pGrid->U[k][je+j][i].Fr2 = pGrid->U[k][js+(j-1)][i].Fr2;
				pGrid->U[k][je+j][i].Sigma_a = pGrid->U[k][js+(j-1)][i].Sigma_a;
				pGrid->U[k][je+j][i].Sigma_t = pGrid->U[k][js+(j-1)][i].Sigma_t;
				pGrid->U[k][je+j][i].Edd_11  =  pGrid->U[k][js+(j-1)][i].Edd_11;
				pGrid->U[k][je+j][i].Edd_21  =  pGrid->U[k][js+(j-1)][i].Edd_21;
				pGrid->U[k][je+j][i].Edd_22  =  pGrid->U[k][js+(j-1)][i].Edd_22;
				pGrid->U[k][je+j][i].Edd_31  =  pGrid->U[k][js+(j-1)][i].Edd_31;
				pGrid->U[k][je+j][i].Edd_32  =  pGrid->U[k][js+(j-1)][i].Edd_32;
				pGrid->U[k][je+j][i].Edd_33  =  pGrid->U[k][js+(j-1)][i].Edd_33;
      				}
    			}
  		}	
	}
	else 
	goto on_error;
	}


  	return;

		on_error:
		ath_error("[BackEuler]: Boundary condition not allowed now!\n");

}

#endif

