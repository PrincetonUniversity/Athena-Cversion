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



static Real ******Divi = NULL;	/* temporary array to store flux for each array */
static Real ******FullAngleV = NULL;	/* temporary array to store vn */
static Real ******MatrixAngleV2 = NULL;	/* temporary array to store vKr/I */


static Real *flux = NULL;

static Real *tempimu = NULL;
static Real *vsource1 = NULL;

static Real *tempV = NULL;
static Real *tempAdv = NULL; /* temporary array for advection flux */
static Real *tempS = NULL;
static Real *tempSV= NULL;
static Real ****RHS = NULL; /* Right hand side of the local matrix for scattering opacity for each cell */
static Real ****Mb = NULL; /* The common elements for each angle */
static Real ****Ma = NULL; /* The special diagonal elements for each angle (each line) */
static Real ****Mc = NULL; /* Matrix elements for the second special matrix */
static Real *Md = NULL;
static Real ****Mdcoef = NULL;
static Real ****inisol = NULL;
static Real ****sol = NULL;
static Real ***Tcoef = NULL;
static Real ***T4coef = NULL;
static Real *****Coefn = NULL;
static Real *****Coefnn = NULL;
static Real ****CoefW = NULL;
static Real *lN1 = NULL; /* lN[0][] */
static Real *lN2 = NULL;/* lN[Nzl][] */
static Real *lN3 = NULL; /* lN[N-1][] */
static Real *UN1 = NULL; /* UN[1][] Upper triangle matrix in LU decomposition */
static Real *UN2 = NULL;/* UN[Nzl][] */
static Real *UN3 = NULL;/* UN[N-1][] */



/* calculate the upwind and downwind intensity for k,j,i, frequency ifr, octant l and angle n */
/* This function is only used to calculate for one ray */
/* in 3D, the octant is numbered as *
 *           | 
 *	1    |    0 
 * ---------------------
 *	3    |    2
 *	     |
 *----------------------
 *---------------------
 *           | 
 *	5    |    4 
 * ---------------------
 *	7    |    6
 *	     |
 *----------------------
 ********************/



void fullRT_3d(DomainS *pD)
{

	RadGridS *pRG=(pD->RadGrid);
	GridS *pG=  pD->Grid;
#ifdef CYLINDRICAL
	const Real *r=pG->r, *ri=pG->ri;
#endif
	Real *Radr;
	int dir;
	int nelements = pRG->nang * pRG->noct;
	Real dt = pG->dt;
	int i, is, ie;
	int j, js, je;
	int k, ks, ke;
	is = pRG->is; ie = pRG->ie;
	js = pRG->js; je = pRG->je;
	ks = pRG->ks; ke = pRG->ke;

	int l, n, ifr,  Mi;
	int offset;
	
	offset = nghost - Radghost;
#ifdef CYLINDRICAL
	Radr = (Real*)&(r[offset]);
#endif
	
	Real dx1, dx2, dx3, ds, alpha, dtods;
	dx1 = pRG->dx1;
	dx2 = pRG->dx2;
	dx3 = pRG->dx3;
	
	Real sigmas, sigmaaI, sigmaa, AngleV, AngleV2, vx, vy, vz, vel2, miux, miuy, miuz;
	
	Real rsf, lsf;
	rsf = 1.0;
	lsf = 1.0;
	
	
	
	
	for(ifr=0; ifr<pRG->nf; ifr++){
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				
				/* First, prepare the data for AngleV, AngleV2 and Vsource3 */
				for(k=0; k<=ke+Radghost; k++){	
					for(j=0; j<=je+Radghost; j++){					
						for(i=0; i<=ie+Radghost; i++){
							/* First, prepare the array */
							
#ifdef CYLINDRICAL						
							miux = pRG->Rphimu[l][n][k][j][i][0];
                            				miuy = pRG->Rphimu[l][n][k][j][i][1];
							miuz = pRG->Rphimu[l][n][k][j][i][2];
#else
							miux = pRG->mu[l][n][k][j][i][0];
							miuy = pRG->mu[l][n][k][j][i][1];
							miuz = pRG->mu[l][n][k][j][i][2];	
#endif								
							
							/*
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							*/
							
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
							
							FullAngleV[ifr][l][n][k][j][i] = miux * vx + miuy * vy + miuz * vz;
							
							MatrixAngleV2[ifr][l][n][k][j][i] = vx * vx * miux * miux + 2.0 * vx * vy * miux * miuy 
												+ 2.0 * vx * vz * miux * miuz + vy * vy * miuy * miuy 
												+ 2.0 * vy * vz * miuy * miuz + vz * vz * miuz * miuz;
							
							
						}/* end i */
					}/* end j */
				}/* end k */
				
				

				/* Now calculate the x flux */
				ds = dx1;
				dtods = dt/ds;
				dir = 1;
				for(k=ks; k<=ke; k++){	
					for(j=js; j<=je; j++){
						
						/* first, prepare the temporary array */
						for(i=0; i<=ie+Radghost; i++){
							/* First, prepare the array */
							sigmas = pRG->R[ifr][k][j][i].Sigma[2];
							/* The absorption opacity in front of I */
							sigmaa = pRG->R[ifr][k][j][i].Sigma[1];
#ifdef CYLINDRICAL						
							miux = pRG->Rphimu[l][n][k][j][i][0];
                            miuy = pRG->Rphimu[l][n][k][j][i][1];
							miuz = pRG->Rphimu[l][n][k][j][i][2];
#else
							miux = pRG->mu[l][n][k][j][i][0];
							miuy = pRG->mu[l][n][k][j][i][1];
							miuz = pRG->mu[l][n][k][j][i][2];	
#endif								
							
							
						/*	vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
						*/
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
							
													
							AngleV = FullAngleV[ifr][l][n][k][j][i];
														
							if((sigmas + sigmaa) > TINY_NUMBER){
								vsource1[i] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
								
								
								tempS[i] = miux * miux * (3.0 * pRG->R[ifr][k][j][i].J);
								if(fabs(miux) > TINY_NUMBER)
									tempSV[i] = vx + miuy * vy/miux + miuz * vz/miux;
								else
									tempSV[i] = 0.0;				
							}else{
								vsource1[i] = 0.0;								
								tempS[i] = 0.0;
								tempSV[i] = 0.0;
								
								
							}
							tempimu[i] = miux * (pRG->imu[ifr][l][n][k][j][i] - vsource1[i]/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);	
							tempV[i] = Crat * alpha * SIGN(miux);		
						}
						
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS,	tempSV,	is, ie+1, ds, dt, tempAdv);
						
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	is, ie+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(i=is; i<=ie; i++){
#ifdef CYLINDRICAL
							rsf = ri[i+1+offset]/r[i+offset];  lsf = ri[i+offset]/r[i+offset];		
#endif						
							
							Divi[ifr][l][n][k][j][i] = 0.5 * (rsf * (tempSV[i+1]+tempSV[i]) * tempAdv[i+1] - lsf * (tempSV[i] + tempSV[i-1]) * tempAdv[i]) * dtods;
							Divi[ifr][l][n][k][j][i]+= Crat * (rsf * flux[i+1] - lsf * flux[i]) * dtods;
							
						}/* End i */
						
					
					}/* End j */
				}/* End k */

				/*---------------------------------------------------------------------*/

				
			/* Now calculate the flux along j direction */
				dir = 2;
				for(k=ks; k<=ke; k++){
					for(i=is; i<=ie; i++){
						ds = dx2;
#ifdef CYLINDRICAL
						/* The scale factor r[i] is the same for each i, for different angles j */
						ds *= r[i+offset];
#endif
						
						dtods = dt/ds;						
						
						/* first save the data */
						for(j=0; j<=je+Radghost; j++){
							/* First, prepare the array */
							sigmas = pRG->R[ifr][k][j][i].Sigma[2];
							/* The absorption opacity in front of I */
							sigmaa = pRG->R[ifr][k][j][i].Sigma[1];
#ifdef CYLINDRICAL						
							miux = pRG->Rphimu[l][n][k][j][i][0];
							miuy = pRG->Rphimu[l][n][k][j][i][1];
							miuz = pRG->Rphimu[l][n][k][j][i][2];
#else
							miux = pRG->mu[l][n][k][j][i][0];
							miuy = pRG->mu[l][n][k][j][i][1];
							miuz = pRG->mu[l][n][k][j][i][2];	
#endif						
							
						/*	vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
						*/
							
							
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
							
							
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							
							
							if((sigmas + sigmaa) > TINY_NUMBER){								
								vsource1[j] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
															
								tempS[j] = miuy * miuy * 3.0 * pRG->R[ifr][k][j][i].J;
								if(fabs(miuy) > TINY_NUMBER)
									tempSV[j] = miux * vx/miuy + vy + miuz * vz/miuy;
								else
									tempSV[j] = 0.0;								
							}else{
								vsource1[j] = 0.0;
								
								tempS[j] = 0.0;
								tempSV[j] = 0.0;							
								
							}
							tempimu[j] = miuy * (pRG->imu[ifr][l][n][k][j][i] - vsource1[j]/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);
							tempV[j] = Crat * alpha * SIGN(miuy);
						}/* End j */				
						
						
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS, tempSV,	js, je+1, ds, dt, tempAdv);
						
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	js, je+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(j=js; j<=je; j++){
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempSV[j+1]+tempSV[j])	* tempAdv[j+1] - (tempSV[j] + tempSV[j-1]) * tempAdv[j]) * dtods;
							Divi[ifr][l][n][k][j][i]+= Crat * (flux[j+1] - flux[j]) * dtods;						
						}/* End j */	
						
					} /* Finish i */
				}/* Finish k */


				/*---------------------------------------------------------------------*/

				/* Now calculate the flux along x3 direction */
				ds = dx3;
				dtods = dt/ds;
				dir = 3;
				for(j=js; j<=je; j++){
					for(i=is; i<=ie; i++){
						/* first save the data */
						for(k=0; k<=ke+Radghost; k++){
							/* First, prepare the array */
							sigmas = pRG->R[ifr][k][j][i].Sigma[2];
							/* The absorption opacity in front of I */
							sigmaa = pRG->R[ifr][k][j][i].Sigma[1];
#ifdef CYLINDRICAL						
							miux = pRG->Rphimu[l][n][k][j][i][0];
							miuy = pRG->Rphimu[l][n][k][j][i][1];
                            miuz = pRG->Rphimu[l][n][k][j][i][2];
#else
							miux = pRG->mu[l][n][k][j][i][0];
							miuy = pRG->mu[l][n][k][j][i][1];
							miuz = pRG->mu[l][n][k][j][i][2];	
#endif	
							
							/*
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							*/
							
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
							
							
							
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							
							
							if((sigmas + sigmaa) > TINY_NUMBER){								
								vsource1[k] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
								
								
								
								tempS[k] = miuz * miuz * 3.0 * pRG->R[ifr][k][j][i].J;
								if(fabs(miuz) > TINY_NUMBER)
									tempSV[k] = miux * vx/miuz + miuy * vy/miuz  + vz;
								else
									tempSV[k] = 0.0;
								
							}else{
								vsource1[k] = 0.0;
								
								tempS[k] = 0.0;
								tempSV[k] = 0.0;
								
								
							}
							tempimu[k] = miuz * (pRG->imu[ifr][l][n][k][j][i] - vsource1[k]/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);
							tempV[k] = Crat * alpha * SIGN(miuz);
						}/* End k */
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS,	tempSV,	ks, ke+1, ds, dt, tempAdv);
						
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	ks, ke+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(k=ks; k<=ke; k++){
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempSV[k+1]+tempSV[k])	* tempAdv[k+1] - (tempSV[k] + tempSV[k-1]) * tempAdv[k]) * dtods;							
							Divi[ifr][l][n][k][j][i]+= Crat * (flux[k+1] - flux[k]) * dtods;						
						}/* End k */	
				
					} /* Finish i */
				}/* Finish j */

			}/* end n */
		}/* end l */
	}/* end ifr */

	
	/****************************************************************/
	
	/* Now we have flux, now add the source terms due to absorption and scattering opacity seperately */ 
	
	/* First, Update the specific intensity with the absorption opacity related terms */
	/* solve the (T^4 - J) and velocity dependent terms together */
	/* This requires using Newton-Raphson to solve a set of non-linear systems */
	
	
	/* first set the gas energy and momentum source terms to be zero, this should be frequency averaged in principle */
	/* Including the ghost zone */
	/* ghost zones for absorption opacity source terms are set with boundary condition function */
	for(k=ks-Radghost; k<=ke+Radghost; k++){
		for(j=js-Radghost; j<=je+Radghost; j++){
			for(i=is-Radghost; i<=ie+Radghost; i++){
				pG->Radheat[k+offset][j+offset][i+offset] = 0.0;
				pG->Pgsource[k+offset][j+offset][i+offset] = 0.0;
				for(l=0; l<3; l++)
					pG->Frsource[k+offset][j+offset][i+offset][l] = 0.0;
			}
		}
	}
	/* first setup the commonly used coefficients */
	for(ifr=0; ifr<pRG->nf; ifr++){
		/* First, set the coefficient for gas temperature */
		/* This coefficients are needed for both RHS and Jacobi coefficient */
		/* We need Planck mean absorption opacity here */
		for(k=ks; k<=ke; k++){
			for(j=js; j<=je; j++){
				for(i=is; i<=ie; i++){
					/*
					 vx = pG->U[0][j+offset][i+offset].M1 / pG->U[0][j+offset][i+offset].d;
					 vy = pG->U[0][j+offset][i+offset].M2 / pG->U[0][j+offset][i+offset].d;
					 vz = pG->U[0][j+offset][i+offset].M3 / pG->U[0][j+offset][i+offset].d;
					 */
					vx = pG->Velguess[k+offset][j+offset][i+offset][0];
					vy = pG->Velguess[k+offset][j+offset][i+offset][1];
					vz = pG->Velguess[k+offset][j+offset][i+offset][2];
				
					vel2 = vx * vx + vy * vy  + vz * vz;
				
					Tcoef[k][j][i] = pG->U[k+offset][j+offset][i+offset].d * R_ideal/Gamma_1;
					T4coef[k][j][i] = dt * Prat * (Crat-vel2/Crat) * pRG->R[ifr][k][j][i].Sigma[0];	
					/* save the current gas temperature */
					inisol[k][j][i][nelements] = pG->tgas[k+offset][j+offset][i+offset];
				
					/* Also set the first guess solution */
					/* set the solution to be the one from last time step */
					sol[k][j][i][nelements] = inisol[k][j][i][nelements];
				
				
				}/* end i */
			}/* End j */
		}/* end k */
		
		
		
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=js; j<=je; j++){
						for(i=is; i<=ie; i++){
							/* Do not include the gas temperature in the Matrix coefficient so that it can be updated */
							Mi = l*(pRG->nang)+n;
							/* We need energy mean absorption opacity here */
							sigmaa = pRG->R[ifr][k][j][i].Sigma[0];
							sigmaaI = pRG->R[ifr][k][j][i].Sigma[1];
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							AngleV2 = MatrixAngleV2[ifr][l][n][k][j][i];
						
							/*
							 vx = pG->U[0][j+offset][i+offset].M1 / pG->U[0][j+offset][i+offset].d;
							 vy = pG->U[0][j+offset][i+offset].M2 / pG->U[0][j+offset][i+offset].d;
							 vz = pG->U[0][j+offset][i+offset].M3 / pG->U[0][j+offset][i+offset].d;
							 */
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
						
							vel2 = vx * vx + vy * vy  + vz * vz;
						
							/* We solve specific intensity directly, do not include the weight in the variable */
						
							Ma[k][j][i][Mi] = (1.0 + dt * sigmaaI * (Crat - AngleV));
							Mb[k][j][i][Mi] = dt * sigmaaI * (vel2 + AngleV2) * pRG->wmu[n][k][j][i]/Crat;
							Mc[k][j][i][Mi] = -4.0 * PI * Prat * dt * sigmaaI * (Crat - vel2/Crat - 2.0 * AngleV) * pRG->wmu[n][k][j][i] - Prat * 4.0 * PI * 2.0 * Mb[k][j][i][Mi];
							/* The Mdcoef is for the original equations, not for the Jacobi matrix */
							/* subtract line nelement from each line so that it becomes diagonal matrix */
							if(Mi < nelements-1)
								Mdcoef[k][j][i][Mi] = -dt * sigmaa * 3.0 * (AngleV - FullAngleV[ifr][pRG->noct-1][pRG->nang-1][k][j][i]) * 0.25 / PI;
							else {
								Mdcoef[k][j][i][Mi] = -dt * sigmaa * (Crat + 3.0 * AngleV) * 0.25 / PI;
							}
						
						
						
							/* set the initial condition */
							inisol[k][j][i][Mi] = pRG->imu[ifr][l][n][k][j][i];
						
							/* The absorption part is I itself, no weight */
							sol[k][j][i][Mi] = inisol[k][j][i][Mi];
						
						
						}/* end i */
					}/* end j */
				}/* end k */
			}/* end n */
		}/* end l */
		
		/* Solve the absorption opacity related terms */
		/* Borrow the memory RHS */
		/* The energy and momentum source terms are already initialized */
		Absorption3D(nelements, pRG, sol, inisol, Ma, Mb, Mc, Mdcoef, Tcoef, T4coef, Md, RHS[ks][js][is]);
		
		
		/* Because the absorption opacity related terms are updated implicitly, 
		 * use the update quantities to calculate the energy and momentum source terms 
		 * for the gas 
		 * The solution is in sol, no weight 
		 */
		/* first, borrow the space of inisol to get the weighted sol */
		for(k=ks; k<=ke; k++){
			for(j=js; j<=je; j++){
				for(i=is; i<=ie; i++){
				
					for(n=0; n<nelements; n++)
						inisol[k][j][i][n] = CoefW[k][j][i][n] * sol[k][j][i][n];
				
					/* The last elements is gas temperature */
				
					inisol[k][j][i][nelements] = sol[k][j][i][nelements];
				}/* end i */
			}/* end j */
		}/* end k */
		
		/* Now call the solution to add the energy and momentum source term */
		/* the flag 1 is sigma[1], for absorption opacity */
		RadAsource3D(ifr, nelements, pRG, pG, Tcoef, Coefn, Coefnn, inisol);

		
		
		
		/*====================================================================*/
	
	/*--------------------------------------------------------------------*/
	/* Now add the source terms with scattering opacity */
	/* All the other terms in the transfer equations are solved together
	 * The independent variable for the matrix solver is actually weight * I
	 * Ma[n] is the diagonal part for each I
	 * Mc[n] is the for the part ndot v (sigma_a+sigma_s)* 3J 
	 * Mb[n] is for the other momentums terms *
	 */
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=js; j<=je; j++){
						for(i=is; i<=ie; i++){
							/* first construct the Right hand side */
							/* We need an array for different angles for each cell. The index order of RHS is different */
							Mi = l*(pRG->nang)+n;
							RHS[k][j][i][Mi] = sol[k][j][i][Mi]- Divi[ifr][l][n][k][j][i];
							sigmas = pRG->R[ifr][k][j][i].Sigma[2];
							/*
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							 */
							vx = pG->Velguess[k+offset][j+offset][i+offset][0];
							vy = pG->Velguess[k+offset][j+offset][i+offset][1];
							vz = pG->Velguess[k+offset][j+offset][i+offset][2];
														
							
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							AngleV2 = MatrixAngleV2[ifr][l][n][k][j][i];
							
							Ma[k][j][i][Mi] = (1.0 + dt * sigmas * (Crat - AngleV))/pRG->wmu[n][k][j][i];
							Mc[k][j][i][Mi] = -dt * sigmas * (Crat + 3.0 * AngleV);
							Mb[k][j][i][Mi] = dt * sigmas * (2.0 * AngleV - (vx * vx + vy * vy + vz * vz + AngleV2)/Crat); 
						
						}/* end i */
					}/* end j */
				}/* end k */
			}/* end n */
		}/* end l */
		
		
		/* solve the first matrix */
		for(k=ks; k<=ke; k++)
			for(j=js; j<=je; j++)
				for(i=is; i<=ie; i++){
					SpecialMatrix3(nelements,Ma[k][j][i], Mb[k][j][i], Mc[k][j][i], Md, RHS[k][j][i], lN1, lN2, lN3, UN1, UN2, UN3);
			}
		
		/* Set the right hand side of the second matrix with the partially updated solution */
		/* The partially updated solution is RHS/wmu */
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=js; j<=je; j++){
						for(i=is; i<=ie; i++){
							Mi = l*(pRG->nang)+n;
							pRG->imu[ifr][l][n][k][j][i] = RHS[k][j][i][Mi]/pRG->wmu[n][k][j][i];
						
						}
					}
				}
			}
		}
		
	}/* End ifr */
	

	/* Moments are updated in the main loop */
	

  return;
}



void fullRT_3d_destruct(void)
{

 

	if(Divi != NULL) free_6d_array(Divi);
	if(FullAngleV != NULL) free_6d_array(FullAngleV);	
	
	if(MatrixAngleV2 != NULL) free_6d_array(MatrixAngleV2);	
	
	if(flux != NULL) free_1d_array(flux);
	

	if(tempimu != NULL) free_1d_array(tempimu);
	if(tempS != NULL) free_1d_array(tempS);
	if(tempSV != NULL) free_1d_array(tempSV);
	if(vsource1 != NULL) free_1d_array(vsource1);
	
	if(tempV != NULL) free_1d_array(tempV);
	if(tempAdv != NULL) free_1d_array(tempAdv);
	if(RHS != NULL) free_4d_array(RHS);
	if(Ma != NULL) free_4d_array(Ma);
	if(Mb != NULL) free_4d_array(Mb);
	if(Mc != NULL) free_4d_array(Mc);
	if(Md != NULL) free_1d_array(Md);
	if(Coefn != NULL) free_5d_array(Coefn);
	if(Coefnn != NULL) free_5d_array(Coefnn);
	if(CoefW != NULL) free_4d_array(CoefW);
	if(Mdcoef != NULL) free_4d_array(Mdcoef);
	if(sol != NULL) free_4d_array(sol);
	if(inisol != NULL) free_4d_array(inisol);
	if(Tcoef != NULL) free_3d_array(Tcoef);
	if(T4coef != NULL) free_3d_array(T4coef);	
	if(lN1 != NULL) free_1d_array(lN1);
	if(lN2 != NULL) free_1d_array(lN2);
	if(lN3 != NULL) free_1d_array(lN3);
	if(UN1 != NULL) free_1d_array(UN1);
	if(UN2 != NULL) free_1d_array(UN2);
	if(UN3 != NULL) free_1d_array(UN3);
	
  return;
}


void fullRT_3d_init(RadGridS *pRG)
{

	
	int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
	int nfr = pRG->nf, noct = pRG->noct, nang = pRG->nang;
	int nmax;

	int l, n, i, j, k, Mi;
	int ie = pRG->ie, je = pRG->je, ke = pRG->ke;
	
	Real miux, miuy, miuz;

	nmax = MAX(nx1,nx2);
	nmax = MAX(nmax,nx3);

	if ((flux = (Real *)calloc_1d_array( nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
		
	if ((Divi = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((FullAngleV = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	

	if ((MatrixAngleV2 = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
		
	if ((tempimu = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempS = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempSV = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((vsource1 = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	
	
	if ((tempV = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempAdv = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((RHS = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Ma = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Mb = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Mc = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	
	if ((Md = (Real *)calloc_1d_array(noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Mdcoef = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	
	if ((Coefn = (Real *****)calloc_5d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, 3, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Coefnn = (Real *****)calloc_5d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, 6, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((CoefW = (Real ****)calloc_4d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	
	if ((Tcoef = (Real ***)calloc_3d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((T4coef = (Real ***)calloc_3d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((sol = (Real ****)calloc_4d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((inisol = (Real ****)calloc_4d_array(nx3+2*Radghost,nx2+2*Radghost,nx1+2*Radghost,noct*nang+1, sizeof(Real))) == NULL)
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
	
	
	/* Store the angle and weight information to calculate the momentums */
	/* This information does not change during the simulation */
	/* So we can calculate it during the initialization after the memory is allocated */
	
	/*------------------------------------------------------*/
	
	
	
	for(l=0; l<noct; l++){
		for(n=0; n<nang; n++){
			/* First, prepare the data for AngleV, AngleV2 and Vsource3 */
			for(k=0; k<=ke+Radghost; k++){
				for(j=0; j<=je+Radghost; j++){	
					for(i=0; i<=ie+Radghost; i++){
					
						Mi = l * nang + n;
					
#ifdef CYLINDRICAL						
					
						miux = pRG->Rphimu[l][n][k][j][i][0]; 
						miuy = pRG->Rphimu[l][n][k][j][i][1];
						miuz = pRG->Rphimu[l][n][k][j][i][2];
#else
						miux = pRG->mu[l][n][k][j][i][0]; 
						miuy = pRG->mu[l][n][k][j][i][1];
						miuz = pRG->mu[l][n][k][j][i][2];
#endif
					
					
					
						Coefn[k][j][i][Mi][0] = miux;
						Coefn[k][j][i][Mi][1] = miuy;
						Coefn[k][j][i][Mi][2] = miuz;
					
						Coefnn[k][j][i][Mi][0] = miux * miux;
						Coefnn[k][j][i][Mi][1] = miux * miuy;
						Coefnn[k][j][i][Mi][2] = miuy * miuy;
						Coefnn[k][j][i][Mi][3] = miux * miuz;
						Coefnn[k][j][i][Mi][4] = miuy * miuz;
						Coefnn[k][j][i][Mi][5] = miuz * miuz;
					
						CoefW[k][j][i][Mi] = pRG->wmu[n][k][j][i];
					
					}/* end i */
				}/* end j */
			}/* end k */
		}/* end n */
	}/* end l */
	
	
	/*------------------------------------------------------*/
	

	return;

	on_error:
  	fullRT_3d_destruct();
  	ath_error("[fullRT_3d_init]: Error allocating memory\n");
  	return;

}
#endif /* FULL_RADIATION_TRANSFER */
