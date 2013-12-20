#include "../copyright.h"
/*==============================================================================
 * FILE: formal_solution.c
 *
 * PURPOSE: integrator for the full radiation transfer in 2D 
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   formal_solution()  - interate formal solution until convergence
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER




static Real ****Divi = NULL;	/* temporary array to store flux for each array */
static Real ***FullAngleV = NULL;	/* temporary array to store vn */
static Real ***MatrixAngleV2 = NULL;	/* temporary array to store vKr/I */


static Real ***flux = NULL;

static Real ***tempimu = NULL;
/*static Real *vsource1 = NULL;
*/
static Real ***tempV = NULL;
/*
static Real *tempAdv = NULL; 
static Real *tempS = NULL;
static Real *tempSV= NULL;
*/
static Real **RHS = NULL; /* Right hand side of the local matrix for scattering opacity for each cell */

static Real ***Ma = NULL; /* The special diagonal elements for each angle (each line) */
/*static Real ***Mc = NULL; */ /* Matrix elements for the second special matrix */
/*static Real ***Mb = NULL;*/ /* The common elements for each angle */

static Real **Mdcoef = NULL;
static Real *Md = NULL;
static Real **inisol = NULL;
static Real **sol = NULL; /* temporary solution used for Newton-Raphon iteration */
static Real **Tcoef = NULL;  /* This is rho R/Gamma_1 for each cell [j][i] */
/*static Real **T4coef = NULL; *//* This is dt P C Sigma_a for each cell [j][i] */
/* prepare the coefficients required to calculate the momentums */
static Real **Coefn = NULL; /* the miux, miuy, miuz for each [j][i][Mi][] */
/*static Real ****Coefnn = NULL;*//* The nn tensor for each [j][i][Mi][] */

static Real *lN1 = NULL;
static Real *lN2 = NULL;
static Real *lN3 = NULL;
static Real *UN1 = NULL; /* Upper triangle matrix in LU decomposition */
static Real *UN2 = NULL; /* Upper triangle matrix in LU decomposition */
static Real *UN3 = NULL; /* Upper triangle matrix in LU decomposition */
static Real InvCrat = 0.0;
static Real QuaPI = 0.0;
static int nelements;


/* calculate the upwind and downwind intensity for k,j,i, frequency ifr, octant l and angle n */
/* This function is only used to calculate for one ray */
/* in 2D, the octant is numbered as *
 *           | 
 *	1    |    0 
 * ---------------------
 *	3    |    2
 *	     |
 ********************/
void Fluxloop2D(RadGridS *pRG, GridS *pG);
void Sourceloop2D(RadGridS *pRG, GridS *pG);


void fullRT_2d(DomainS *pD)
{

	RadGridS *pRG=(pD->RadGrid);
	GridS *pG = pD->Grid;

	int i, is, ie;
	int j, js, je;
	int ks;
	is = pG->is; ie = pG->ie;
	js = pG->js; je = pG->je;
	ks = pG->ks;
    
    int l;

	
	InvCrat=1.0/Crat;
	
  	/* calculate the flux for three directions */
    Fluxloop2D(pRG, pG);
    


	/* Now we have flux, now add the source terms due to absorption and scattering opacity seperately */ 
	
	/* First, Update the specific intensity with the absorption opacity related terms */
	/* solve the (T^4 - J) and velocity dependent terms together */
	/* This requires using Newton-Raphson to solve a set of non-linear systems */
	
	
	/* first set the gas energy and momentum source terms to be zero, this should be frequency averaged in principle */
	/* Including the ghost zone */
	/* ghost zones for absorption opacity source terms are set with boundary condition function */
	for(j=js-nghost; j<=je+nghost; j++){
		for(i=is-nghost; i<=ie+nghost; i++){
			pG->Radheat[ks][j][i] = 0.0;
			pG->Pgsource[ks][j][i] = 0.0;
			for(l=0; l<3; l++)
				pG->Frsource[ks][j][i][l] = 0.0;
		}
	}
    
    
    /****************************************************************/
	
	/* Add the absorption and scattering opacity terms */
    
        Sourceloop2D(pRG, pG);
	
    /****************************************************************/
    
	/* Moments are updated in the main loop */
    
	

  return;
}


void Fluxloop2D(RadGridS *pRG, GridS *pG)
{


#ifdef CYLINDRICAL
	const Real *r=pG->r, *ri=pG->ri;
#endif
	Real *Radr = NULL;
	int dir; 

	Real dt = pG->dt;
	int i, is, ie;
	int j, js, je;
	int ks;
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	ks = pRG->ks;

	int Mi;
	int offset, ifr;


	offset = nghost - Radghost;
#ifdef CYLINDRICAL
	Radr = (Real*)&(r[offset]);
#endif

	Real dx1, dx2, ds, alpha, dtods;
	dx1 = pRG->dx1;
	dx2 = pRG->dx2;

	Real sigmas, sigmaa, AngleV, vx, vy, vz, miux, miuy, miuz;
    Real temp;


	Real rsf, lsf;
	rsf = 1.0;
	lsf = 1.0;
    
    
    /* First, prepare the data for AngleV, AngleV2 and Vsource3 */
	for(j=0; j<=je+Radghost; j++){
		for(i=0; i<=ie+Radghost; i++){
                
                
			/* velocity is independent of the angles */
            /*
              vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
              vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
              vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
            */
                
            vx = pG->Velguess[0][j+offset][i+offset][0];
			vy = pG->Velguess[0][j+offset][i+offset][1];
			vz = pG->Velguess[0][j+offset][i+offset][2];
                
                
                

            for(Mi=0; Mi<nelements; Mi++){
                   /* First, prepare the array */
#ifdef CYLINDRICAL
						
				miux = pRG->Rphimu[ks][j][i][Mi][0];
				miuy = pRG->Rphimu[ks][j][i][Mi][1];
				miuz = pRG->Rphimu[ks][j][i][Mi][2];
#else
				miux = pRG->mu[ks][j][i][Mi][0];
				miuy = pRG->mu[ks][j][i][Mi][1];
				miuz = pRG->mu[ks][j][i][Mi][2];
#endif
                        
				FullAngleV[j][i][Mi] = miux * vx + miuy * vy + miuz * vz;
						
				MatrixAngleV2[j][i][Mi] = vx * vx * miux * miux + 2.0 * vx * vy * miux * miuy
                       + 2.0 * vx * vz * miux * miuz + vy * vy * miuy * miuy
                       + 2.0 * vy * vz * miuy * miuz + vz * vz * miuz * miuz;
                        
                        
             }/* end Mi */

		}/* end i */
	}/* end j */


    
    /* Now calculate the x flux */
	ds = dx1;
	dtods = dt/ds;
	dir = 1;

	for(j=js; j<=je; j++){
		/* first, prepare the temporary array */
		for(i=0; i<=ie+Radghost; i++){
                
           /* velocity is independent of the angles */
                
                
           /*		vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
             vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
           vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
                 */
                
            vx = pG->Velguess[0][j+offset][i+offset][0];
			vy = pG->Velguess[0][j+offset][i+offset][1];
			vz = pG->Velguess[0][j+offset][i+offset][2];
                
                
			for(ifr=0; ifr<pRG->nf; ifr++){
				/* First, prepare the array */
				sigmas = pRG->R[0][j][i][ifr].Sigma[2];
				/* The absorption opacity in front of I */
				sigmaa = pRG->R[0][j][i][ifr].Sigma[1];
                
                alpha = pRG->Speedfactor[ks][j][i][ifr][0];
                    
				for(Mi=0; Mi<nelements; Mi++){
						
                        
#ifdef CYLINDRICAL
						
					miux = pRG->Rphimu[ks][j][i][Mi][0];
					miuy = pRG->Rphimu[ks][j][i][Mi][1];
					miuz = pRG->Rphimu[ks][j][i][Mi][2];
#else
					miux = pRG->mu[ks][j][i][Mi][0];
					miuy = pRG->mu[ks][j][i][Mi][1];
					miuz = pRG->mu[ks][j][i][Mi][2];
#endif
						
                        
                    AngleV = FullAngleV[j][i][Mi];
                        
                        
                    if((sigmas + sigmaa) > TINY_NUMBER){
						temp = AngleV * (3.0 * pRG->R[ks][j][i][ifr].J);
							
                    }else{
						temp = 0.0;
					}
                        
					tempimu[i][ifr][2*Mi] = miux * (pRG->imu[ks][j][i][ifr][Mi] - temp * InvCrat);
					tempV[i][ifr][2*Mi] = Crat * alpha * SIGN(miux);
                        
                        
                    if((sigmas + sigmaa) > TINY_NUMBER){
                            
						tempimu[i][ifr][2*Mi+1] = miux * miux * (3.0 * pRG->R[ks][j][i][ifr].J);
						if(fabs(miux) > TINY_NUMBER)
							tempV[i][ifr][2*Mi+1] = vx + miuy * vy/miux + miuz * vz/miux;
						else
							tempV[i][ifr][2*Mi+1] = 0.0;
					}else{
                            
						tempimu[i][ifr][2*Mi+1] = 0.0;
						tempV[i][ifr][2*Mi+1] = 0.0;
                            
                            
					}
                }/* end Mi */
            }/* end ifr */
        }/* end i */
            
            
            /* first, calculate the advection part v(3J+I) */
            /*	flux_AdvJ(Radr, dir, tempS,	tempSV,	is, ie+1, ds, dt, tempAdv);
             */
            /* Third, calculate flux due to co-moving miu */
            /* The two parts are calculated together */
            flux_AdvJ(pRG->nf, 2*nelements, Radr, dir, tempimu,	tempV,	is, ie+1, ds, dt, flux);
            
            /* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
            for(i=is; i<=ie; i++){
#ifdef CYLINDRICAL
                rsf = ri[i+1+offset]/r[i+offset];  lsf = ri[i+offset]/r[i+offset];
#endif
                
                for(ifr=0; ifr<pRG->nf; ifr++){
                    for(Mi=0; Mi<nelements; Mi++){
                        Divi[j][i][ifr][Mi]  = Crat * (rsf * flux[i+1][ifr][2*Mi] - lsf * flux[i][ifr][2*Mi]) * dtods;
                        Divi[j][i][ifr][Mi] += 0.5 * (rsf * (tempV[i+1][ifr][2*Mi+1] + tempV[i][ifr][2*Mi+1]) * flux[i+1][ifr][2*Mi+1] - lsf * (tempV[i][ifr][2*Mi+1] + tempV[i-1][ifr][2*Mi+1]) * flux[i][ifr][2*Mi+1]) * dtods;
                        
                    }/* end angles */
                }/* end frequency */
            }/* End i */
            
            
    }/* End j */


 			
	/* Now calculate the flux along j direction */
	dir = 2;

    for(i=is; i<=ie; i++){
		ds = dx2;
#ifdef CYLINDRICAL
           /* The scale factor r[i] is the same for each i, for different angles j */
		ds *= r[i+offset];
#endif
            
		dtods = dt/ds;
            
		/* first save the data */
		for(j=0; j<=je+Radghost; j++){
                
         /*	vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
            vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
            vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
        */
                
			vx = pG->Velguess[0][j+offset][i+offset][0];
			vy = pG->Velguess[0][j+offset][i+offset][1];
			vz = pG->Velguess[0][j+offset][i+offset][2];
                
                
            for(ifr=0; ifr<pRG->nf; ifr++){
				/* First, prepare the array */
				sigmas = pRG->R[ks][j][i][ifr].Sigma[2];
				/* The absorption opacity in front of I */
				sigmaa = pRG->R[ks][j][i][ifr].Sigma[1];
                
                alpha = pRG->Speedfactor[ks][j][i][ifr][1];
                    
				for(Mi=0; Mi<nelements; Mi++){
                        
#ifdef CYLINDRICAL
						
                    miux = pRG->Rphimu[ks][j][i][Mi][0];
					miuy = pRG->Rphimu[ks][j][i][Mi][1];
					miuz = pRG->Rphimu[ks][j][i][Mi][2];
#else
					miux = pRG->mu[ks][j][i][Mi][0];
					miuy = pRG->mu[ks][j][i][Mi][1];
					miuz = pRG->mu[ks][j][i][Mi][2];
#endif
                        
                        
                        
                    AngleV = FullAngleV[j][i][Mi];
                        
                        
                        
                    if((sigmas + sigmaa) > TINY_NUMBER){
						temp = AngleV * (3.0 * pRG->R[ks][j][i][ifr].J);
							
                    }else{
						temp = 0.0;
					}
                        
					tempimu[j][ifr][2*Mi] = miuy * (pRG->imu[ks][j][i][ifr][Mi] - temp * InvCrat);
					tempV[j][ifr][2*Mi] = Crat * alpha * SIGN(miuy);
                        
                        
                    if((sigmas + sigmaa) > TINY_NUMBER){
                            
						tempimu[j][ifr][2*Mi+1] = miuy * miuy * (3.0 * pRG->R[ks][j][i][ifr].J);
						if(fabs(miuy) > TINY_NUMBER)
							tempV[j][ifr][2*Mi+1] = miux * vx/miuy + vy + miuz * vz/miuy;
						else
							tempV[j][ifr][2*Mi+1] = 0.0;
					}else{
                            
						tempimu[j][ifr][2*Mi+1] = 0.0;
						tempV[j][ifr][2*Mi+1] = 0.0;
                            
                            
					}
                }/* end Mi */
            }/* end ifr */
        }/* End j */
            
            
            
            /* first, calculate the advection part v(3J+I) */
            /*		flux_AdvJ(Radr, dir, tempS, tempSV,	js, je+1, ds, dt, tempAdv);
             */
            /* Third, calculate flux due to co-moving miu */
            flux_AdvJ(pRG->nf, 2*nelements, Radr, dir, tempimu,	tempV,	js, je+1, ds, dt, flux);
            
            /* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
            for(j=js; j<=je; j++){
                for(ifr=0; ifr<pRG->nf; ifr++){
                    for(Mi=0; Mi<nelements; Mi++){
                        Divi[j][i][ifr][Mi] += Crat * (flux[j+1][ifr][2*Mi] - flux[j][ifr][2*Mi]) * dtods;
                        Divi[j][i][ifr][Mi] += 0.5 * ((tempV[j+1][ifr][2*Mi+1] + tempV[j][ifr][2*Mi+1])	* flux[j+1][ifr][2*Mi+1] - (tempV[j][ifr][2*Mi+1] + tempV[j-1][ifr][2*Mi+1]) * flux[j][ifr][2*Mi+1]) * dtods;
                    }
                }/* end nf */						
            }/* End j */	
            
    } /* Finish i */

    


return;


}



void Sourceloop2D(RadGridS *pRG, GridS *pG)
{
    
    
	Real dt = pG->dt;
	int i, is, ie;
	int j, js, je;
	int ks;
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	ks = pRG->ks;
    
	int offset, ifr, Mi, flag;
	
	offset = nghost - Radghost;
	
	Real vx, vy, vz, vel2, AngleV, AngleV2, miux, miuy, miuz, AngleVN;
	Real sigmaa, sigmaaI, sigmas;
    Real Tgas4;
	

		for(j=js; j<=je; j++){
			for(i=is; i<=ie; i++){
                
				/* First, set the coefficient for gas temperature */
				/* This coefficients are needed for both RHS and Jacobi coefficient */
				/* We need Planck mean absorption opacity here */
                
				/*
                 vx = pG->U[0][j+offset][i+offset].M1 / pG->U[0][j+offset][i+offset].d;
                 vy = pG->U[0][j+offset][i+offset].M2 / pG->U[0][j+offset][i+offset].d;
                 vz = pG->U[0][j+offset][i+offset].M3 / pG->U[0][j+offset][i+offset].d;
				 */
                vx = pG->Velguess[0][j+offset][i+offset][0];
                vy = pG->Velguess[0][j+offset][i+offset][1];
                vz = pG->Velguess[0][j+offset][i+offset][2];
				
                
                vel2 = vx * vx + vy * vy + vz * vz;
                /* Prepare the matrix coefficient */
                /* First, calculate the last element */
                
                AngleVN = FullAngleV[j][i][nelements-1];
                
                /*-----------------------------------------------------------------------------*/
                /* First, the absorption opacity */
                
                
                for(ifr=0; ifr<pRG->nf; ifr++){
                    /* Matrix coefficient for absorption opacity */
                    sigmaa = pRG->R[ks][j][i][ifr].Sigma[0];
                    sigmaaI = pRG->R[ks][j][i][ifr].Sigma[1];
					
                    Tcoef[ifr][0] = pG->U[0][j+offset][i+offset].d * R_ideal/Gamma_1;
                    Tcoef[ifr][1] = dt * Prat * (Crat-vel2 * InvCrat) * sigmaa;
                    /* save the current gas temperature */
                    
                    for(Mi=0; Mi<nelements; Mi++){
#ifdef CYLINDRICAL
						
                        miux = pRG->Rphimu[ks][j][i][Mi][0];
                        miuy = pRG->Rphimu[ks][j][i][Mi][1];
                        miuz = pRG->Rphimu[ks][j][i][Mi][2];
#else
                        miux = pRG->mu[ks][j][i][Mi][0];
                        miuy = pRG->mu[ks][j][i][Mi][1];
                        miuz = pRG->mu[ks][j][i][Mi][2];
#endif
                        
                        /* calculate temporary array */
                        /* This is used to add radiation source term */
                        Coefn[Mi][0] = miux;
                        Coefn[Mi][1] = miuy;
                        Coefn[Mi][2] = miuz;
                        Coefn[Mi][3] = miux * miux;
                        Coefn[Mi][4] = miux * miuy;
                        Coefn[Mi][5] = miuy * miuy;
                        Coefn[Mi][6] = miux * miuz;
                        Coefn[Mi][7] = miuy * miuz;
                        Coefn[Mi][8] = miuz * miuz;
                        
						
                        AngleV = FullAngleV[j][i][Mi];
                        AngleV2 = MatrixAngleV2[j][i][Mi];
                        
                        
                        Ma[ifr][Mi][0] = (1.0 + dt * sigmaaI * (Crat - AngleV));
                        Ma[ifr][Mi][1] = dt * sigmaaI * (vel2 + AngleV2) * pRG->wmu[ks][j][i][Mi] * InvCrat;
                        Ma[ifr][Mi][2] = -4.0 * PI * Prat * dt * sigmaaI * (Crat - vel2 * InvCrat - 2.0 * AngleV) * pRG->wmu[ks][j][i][Mi] - Prat * 4.0 * PI * 2.0 * Ma[ifr][Mi][1];
                        /* The Mdcoef is for the original equations, not for the Jacobi matrix */
                        /* subtract line nelement from each line so that it becomes diagonal matrix */
                        if(Mi < nelements-1)
                            Mdcoef[ifr][Mi] = -dt * sigmaa * 3.0 * (AngleV - AngleVN) * QuaPI;
                        else {
                            Mdcoef[ifr][Mi] = -dt * sigmaa * (Crat + 3.0 * AngleV) * QuaPI;
                        }
						
                        
                        /* set the initial condition */
                        inisol[ifr][Mi] = pRG->imu[ks][j][i][ifr][Mi];
						
                        /* The absorption part is I itself, no weight */
                        sol[ifr][Mi] = inisol[ifr][Mi];
						
                    }/* end Mi */
					
					
                    inisol[ifr][nelements] = pG->tgas[0][j+offset][i+offset];
                    
                    /* Also set the first guess solution */
                    /* set the solution to be the one from last time step */
                    sol[ifr][nelements] = inisol[ifr][nelements];
                    
                }/* end ifr */
                
                
                /* Solve the absorption opacity related terms */
                /* Borrow the memory RHS */
                /* The energy and momentum source terms are already initialized */
                Absorption(pRG->nf, nelements, sol, inisol, Ma, Mdcoef, Tcoef, Md, RHS[0], &flag);
                
                
                /**********************************************************************/
                /* This part is only calculated in special cases */
                /* If the non-linear matrix does not converge, fix the temperature and invert a linear
                 * matrix as a rough estimate */
                if(flag == 0){
                    for(ifr=0; ifr<pRG->nf; ifr++){
                        sigmaa = pRG->R[ks][j][i][ifr].Sigma[0];
                        sigmaaI = pRG->R[ks][j][i][ifr].Sigma[1];
                        
                        Tgas4 = pG->tgas[ks][j+offset][i+offset];
                        Tgas4 = SQR(Tgas4);
                        Tgas4 = SQR(Tgas4);
                        /* save the current gas temperature */
                        
                        AngleV = FullAngleV[j][i][nelements-1];
                        
                        RHS[0][nelements-1] = pRG->imu[ks][j][i][ifr][nelements-1] + dt * Crat * sigmaa * Tgas4 * QuaPI + 3.0 * dt * AngleV * sigmaa * QuaPI;
                        
                        for(Mi=0; Mi<nelements; Mi++){
                            
                            
                            AngleV = FullAngleV[j][i][Mi];
                            AngleV2 = MatrixAngleV2[j][i][Mi];
                            
                            Ma[0][Mi][0] = (1.0 + dt * sigmaaI * (Crat - AngleV));
                            Ma[0][Mi][1] = dt * sigmaaI * (vel2 + AngleV2) * pRG->wmu[ks][j][i][Mi] * InvCrat;
                            
                            if(Mi < nelements - 1){
                                RHS[0][Mi] = pRG->imu[ks][j][i][ifr][Mi] + dt * Crat * sigmaa * Tgas4 * QuaPI + 3.0 * dt * AngleV * sigmaa * QuaPI;
                                RHS[0][Mi] -= RHS[0][nelements-1];
                            }
                            
                        }/* end Mi */
                        
                        /* Now solve the linear matrix */
                        LinearAbsorptionMatrix(nelements, Ma[0], RHS[0]);
                        
                        /* set the solution */
                        for(Mi=0; Mi<nelements; Mi++){
                            sol[ifr][Mi] = RHS[0][Mi];
                        }
                        
                        sol[ifr][nelements] = pG->tgas[ks][j+offset][i+offset];
                        
                        
                    }/* end ifr */
                    
                }/* end flag */

                
                /* Because the absorption opacity related terms are updated implicitly,
                 * use the update quantities to calculate the energy and momentum source terms
                 * for the gas
                 * The solution is in sol, no weight
                 */
                /* first, borrow the space of inisol to get the weighted sol */
                
                for(ifr=0; ifr<pRG->nf; ifr++){
                    for(Mi=0; Mi<nelements; Mi++){
                        inisol[ifr][Mi] = pRG->wmu[ks][j][i][Mi] * sol[ifr][Mi];
                    }/* end nelement */
                    
                    /* The last one is the gas temperature */
                    inisol[ifr][nelements] = sol[ifr][nelements];
                }/* end ifr */
                
				
                /* Now call the solution to add the energy and momentum source term */
                /* Only do this if absorption matrix is converged properly */
                if(flag > 0)
                    RadAsource(i, j, ks, nelements, pRG, pG, Tcoef, Coefn, inisol);
                
                
                
                /*--------------------------------------------------------------------------------------------------------------*/
                /* Now the scattering opacity */
                for(ifr=0; ifr<pRG->nf; ifr++){
                    
                    sigmas = pRG->R[ks][j][i][ifr].Sigma[2];
                    
                    
                    for(Mi=0; Mi<nelements; Mi++){
						
                        /* Set RHS */
                        RHS[ifr][Mi] = sol[ifr][Mi]- Divi[j][i][ifr][Mi];
						
                        AngleV = FullAngleV[j][i][Mi];
                        AngleV2 = MatrixAngleV2[j][i][Mi];
                        
                        Ma[ifr][Mi][0] = (1.0 + dt * sigmas * (Crat - AngleV))/pRG->wmu[ks][j][i][Mi];
                        Ma[ifr][Mi][1] = dt * sigmas * (2.0 * AngleV - (vel2 + AngleV2) * InvCrat);
                        Ma[ifr][Mi][2] = -dt * sigmas * (Crat + 3.0 * AngleV);
                    }/* end Mi */
                    
                    
                    /* Solve the scattering matrix for each frequency */
                    SpecialMatrix3(nelements, Ma[ifr], Md, RHS[ifr], lN1, lN2, lN3, UN1, UN2, UN3);
                    
                }/* end ifr */
                
                
   
                
                /*--------------------------------------------------------------------------------------------------------------*/
                /* Now set the solution */
                for(ifr=0; ifr<pRG->nf; ifr++){
                    for(Mi=0; Mi<nelements; Mi++){
                        pRG->imu[ks][j][i][ifr][Mi] = RHS[ifr][Mi]/pRG->wmu[ks][j][i][Mi];
                    }/* Mi */
                }/* ifr */
                
     
                
            }/* end i */
        }/* End j */

	
    
    return;	
    
    
}


void fullRT_2d_destruct(void)
{

  if(Divi != NULL) free_4d_array(Divi);
  if(FullAngleV != NULL) free_3d_array(FullAngleV);
	
  if(MatrixAngleV2 != NULL) free_3d_array(MatrixAngleV2);
		
  if(flux != NULL) free_3d_array(flux);
  if(tempimu != NULL) free_3d_array(tempimu);
/*  if(tempS != NULL) free_1d_array(tempS);
  if(tempSV != NULL) free_1d_array(tempSV);
  if(vsource1 != NULL) free_1d_array(vsource1);
*/
  if(tempV != NULL) free_3d_array(tempV);
/*  if(tempAdv != NULL) free_1d_array(tempAdv);
*/
  if(RHS != NULL) free_2d_array(RHS);
  if(Ma != NULL) free_3d_array(Ma);
/*  if(Mb != NULL) free_3d_array(Mb);
  if(Mc != NULL) free_3d_array(Mc);
*/
  if(Md != NULL) free_1d_array(Md);
  if(Coefn != NULL) free_2d_array(Coefn);
/*  if(Coefnn != NULL) free_4d_array(Coefnn);
 
  if(CoefW != NULL) free_3d_array(CoefW);
 */
  if(Mdcoef != NULL) free_2d_array(Mdcoef);
  if(sol != NULL) free_2d_array(sol);
  if(inisol != NULL) free_2d_array(inisol);
  if(Tcoef != NULL) free_2d_array(Tcoef);
/*  if(T4coef != NULL) free_2d_array(T4coef);
 */
  if(lN1 != NULL) free_1d_array(lN1);
  if(lN2 != NULL) free_1d_array(lN2);
  if(lN3 != NULL) free_1d_array(lN3);
  if(UN1 != NULL) free_1d_array(UN1);
  if(UN2 != NULL) free_1d_array(UN2);
  if(UN3 != NULL) free_1d_array(UN3);

  return;
}


void fullRT_2d_init(RadGridS *pRG)
{


	/* constants */
	QuaPI = 0.25/PI;
	
	int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1];
	int nfr = pRG->nf, noct = pRG->noct, nang = pRG->nang;
	int nmax;
	
    
	nelements = pRG->nang * pRG->noct;
	


	

	


	nmax = MAX(nx1,nx2);

	if ((flux = (Real ***)calloc_3d_array( nmax+2*Radghost, nfr, 2*nelements, sizeof(Real))) == NULL)
    		goto on_error;
	
	if ((Divi = (Real ****)calloc_4d_array(nx2+2*Radghost, nx1+2*Radghost, nfr, nelements, sizeof(Real))) == NULL)
    		goto on_error;
	
	if ((FullAngleV = (Real ***)calloc_3d_array(nx2+2*Radghost, nx1+2*Radghost, nelements, sizeof(Real))) == NULL)
		goto on_error;
	

	if ((MatrixAngleV2 = (Real ***)calloc_3d_array(nx2+2*Radghost, nx1+2*Radghost, nelements, sizeof(Real))) == NULL)
		goto on_error;
	
	
	if ((tempimu = (Real ***)calloc_3d_array(nmax+2*Radghost, nfr, 2*nelements, sizeof(Real))) == NULL)
		goto on_error;

/*    if ((tempS = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
                goto on_error;
	
	if ((tempSV = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
                goto on_error;

	if ((vsource1 = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
                goto on_error;
*/
    if ((tempV = (Real ***)calloc_3d_array(nmax+2*Radghost, nfr, 2*nelements, sizeof(Real))) == NULL)
                goto on_error;
/*
	if ((tempAdv = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
                goto on_error;
*/	
	/* Each column has noct * nang + 1 elements to make them work for absorption case */

    if ((RHS = (Real **)calloc_2d_array(nfr,noct*nang+1, sizeof(Real))) == NULL)
                goto on_error;

	if ((Ma = (Real ***)calloc_3d_array(nfr, noct*nang+1, 3, sizeof(Real))) == NULL)
                goto on_error;

/*	if ((Mb = (Real ***)calloc_3d_array(nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
                goto on_error;
	
	if ((Mc = (Real ***)calloc_3d_array(nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
                goto on_error;
*/	
	
	if ((Coefn = (Real **)calloc_2d_array(noct*nang, 9, sizeof(Real))) == NULL)
		goto on_error;
/*	
	if ((Coefnn = (Real ****)calloc_4d_array(nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, 6, sizeof(Real))) == NULL)
		goto on_error;
*/	

	
	
	if ((Tcoef = (Real **)calloc_2d_array(nfr, 2, sizeof(Real))) == NULL)
		goto on_error;

/*	if ((T4coef = (Real **)calloc_2d_array(nx2+2*Radghost,nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
*/

	if ((Md = (Real *)calloc_1d_array(noct*nang+1, sizeof(Real))) == NULL)
        goto on_error;
	
	if ((Mdcoef = (Real **)calloc_2d_array(nfr,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((sol = (Real **)calloc_2d_array(nfr, noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((inisol = (Real **)calloc_2d_array(nfr, noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;

	
	
	if ((lN1 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;

	if ((lN2 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;

	if ((lN3 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;

	if ((UN1 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;

	if ((UN2 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;

	if ((UN3 = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
                goto on_error;
	
	
	

	
	/*------------------------------------------------------*/

	return;

	on_error:
  	fullRT_2d_destruct();
  	ath_error("[fullRT_2d_init]: Error allocating memory\n");
  	return;




}
#endif /* FULL_RADIATION_TRANSFER */
