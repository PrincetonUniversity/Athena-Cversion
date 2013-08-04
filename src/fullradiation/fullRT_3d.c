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
static Real ******FullAngleV2 = NULL;	/* temporary array to store vPr */
static Real ******FullVsource3 = NULL;	/* temporary array to store vsource3 */

static Real *flux = NULL;
static Real *fluxsource3 = NULL;
static Real *tempV3V = NULL; /* The advection velocity for source 3 */
static Real *tempimu = NULL;
static Real *vsource1 = NULL;
static Real *vsource2 = NULL;
static Real *vsource3 = NULL;
static Real *tempV = NULL;
static Real *tempAdv = NULL; /* temporary array for advection flux */
static Real *tempS = NULL;
static Real *tempSV= NULL;
static Real ****RHS = NULL; /* Right hand side of the local matrix for scattering opacity for each cell */
static Real ****Mb = NULL; /* The common elements for each angle */
static Real ****Ma = NULL; /* The special diagonal elements for each angle (each line) */
static Real ****Mc = NULL; /* Matrix elements for the second special matrix */
static Real *Md = NULL;
static Real **lN = NULL;
static Real **UN = NULL; /* Upper triangle matrix in LU decomposition */



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
	
	Real sigmas, sigmaa, AngleV, AngleV2, vx, vy, vz, miux, miuy, miuz;
	
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
							
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							
							FullAngleV[ifr][l][n][k][j][i] = miux * vx + miuy * vy + miuz * vz;
							FullAngleV2[ifr][l][n][k][j][i] = vx * vx * pRG->R[ifr][k][j][i].K[0] + 2.0 * vx * vy * pRG->R[ifr][k][j][i].K[1] 
															+ 2.0 * vx * vz * pRG->R[ifr][k][j][i].K[3] + vy * vy * pRG->R[ifr][k][j][i].K[2] 
															+ 2.0 * vy * vz * pRG->R[ifr][k][j][i].K[4] + vz * vz * pRG->R[ifr][k][j][i].K[5];
							
							FullVsource3[ifr][l][n][k][j][i] = -2.0 * sigmas * (vx * pRG->R[ifr][k][j][i].H[0] + vy * pRG->R[ifr][k][j][i].H[1] 
								+ vz * pRG->R[ifr][k][j][i].H[2])/(sigmaa+sigmas)
								+ (sigmas-sigmaa)* ((vx*vx+vy*vy+vz*vz)*pRG->R[ifr][k][j][i].J+FullAngleV2[ifr][l][n][k][j][i])/(Crat*(sigmas+sigmaa));
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
							
							
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							AngleV = FullAngleV[ifr][l][n][k][j][i];
														
							if((sigmas + sigmaa) > TINY_NUMBER){
								vsource1[i] = AngleV * (pRG->imu[ifr][l][n][k][j][i]);
								vsource2[i] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
								vsource3[i] = FullVsource3[ifr][l][n][k][j][i];
								
								tempS[i] = miux * miux * (3.0 * pRG->R[ifr][k][j][i].J + pRG->imu[ifr][l][n][k][j][i]);
								if(fabs(miux) > TINY_NUMBER)
									tempSV[i] = vx + miuy * vy/miux + miuz * vz/miux;
								else
									tempSV[i] = 0.0;
								tempV3V[i] = miux;
							}else{
								vsource1[i] = 0.0;
								vsource2[i] = 0.0;
								vsource3[i] = 0.0;
								tempS[i] = 0.0;
								tempSV[i] = 0.0;
								tempV3V[i] = 0.0;
								
							}
							tempimu[i] = miux * (pRG->imu[ifr][l][n][k][j][i] - (vsource1[i] + vsource2[i] + vsource3[i])/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);	
							tempV[i] = Crat * alpha * SIGN(miux);		
						}
						
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS,	tempSV,	is, ie+1, ds, dt, tempAdv);
						/* second, calculate the flux due to vsource3 */
						flux_AdvJ(Radr, dir, vsource3, tempV3V,	is, ie+1, ds, dt, fluxsource3);
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	is, ie+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(i=is; i<=ie; i++){
#ifdef CYLINDRICAL
							rsf = ri[i+1+offset]/r[i+offset];  lsf = ri[i+offset]/r[i+offset];		
#endif						
							
							Divi[ifr][l][n][k][j][i] = 0.5 * (rsf * (tempSV[i+1]+tempSV[i]) * tempAdv[i+1] - lsf * (tempSV[i] + tempSV[i-1]) * tempAdv[i]) * dtods;
							Divi[ifr][l][n][k][j][i]+= 0.5 * (rsf * (tempV3V[i+1]+tempV3V[i]) * fluxsource3[i+1] - lsf * (tempV3V[i] + tempV3V[i-1]) * fluxsource3[i]) * dtods;
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
							
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							
							
							if((sigmas + sigmaa) > TINY_NUMBER){
								vsource1[j] = AngleV * (pRG->imu[ifr][l][n][k][j][i]);
								vsource2[j] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
								vsource3[j] = FullVsource3[ifr][l][n][k][j][i];
								
								
								tempS[j] = miuy * miuy * (3.0 * pRG->R[ifr][k][j][i].J + pRG->imu[ifr][l][n][k][j][i]);
								if(fabs(miuy) > TINY_NUMBER)
									tempSV[j] = miux * vx/miuy + vy + miuz * vz/miuy;
								else
									tempSV[j] = 0.0;
								tempV3V[j] = miuy;
							}else{
								vsource1[j] = 0.0;
								vsource2[j] = 0.0;
								vsource3[j] = 0.0;
								tempS[j] = 0.0;
								tempSV[j] = 0.0;
								tempV3V[j] = 0.0;
								
							}
							tempimu[j] = miuy * (pRG->imu[ifr][l][n][k][j][i] - (vsource1[j] + vsource2[j] + vsource3[j])/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);
							tempV[j] = Crat * alpha * SIGN(miuy);
						}/* End j */				
						
						
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS, tempSV,	js, je+1, ds, dt, tempAdv);
						/* second, calculate the flux due to vsource3 */
						flux_AdvJ(Radr, dir, vsource3, tempV3V,	js, je+1, ds, dt, fluxsource3);
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	js, je+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(j=js; j<=je; j++){
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempSV[j+1]+tempSV[j])	* tempAdv[j+1] - (tempSV[j] + tempSV[j-1]) * tempAdv[j]) * dtods;
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempV3V[j+1]+tempV3V[j]) * fluxsource3[j+1] - (tempV3V[j] + tempV3V[j-1]) * fluxsource3[j]) * dtods;
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
							
							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							
							
							if((sigmas + sigmaa) > TINY_NUMBER){
								vsource1[k] = AngleV * (pRG->imu[ifr][l][n][k][j][i]);
								vsource2[k] = AngleV * (3.0 * pRG->R[ifr][k][j][i].J);
								vsource3[k] = FullVsource3[ifr][l][n][k][j][i];
								
								
								tempS[k] = miuz * miuz * (3.0 * pRG->R[ifr][k][j][i].J + pRG->imu[ifr][l][n][k][j][i]);
								if(fabs(miuz) > TINY_NUMBER)
									tempSV[k] = miux * vx/miuz + miuy * vy/miuz  + vz;
								else
									tempSV[k] = 0.0;
								tempV3V[k] = miuz;
							}else{
								vsource1[k] = 0.0;
								vsource2[k] = 0.0;
								vsource3[k] = 0.0;
								tempS[k] = 0.0;
								tempSV[k] = 0.0;
								tempV3V[k] = 0.0;
								
							}
							tempimu[k] = miuz * (pRG->imu[ifr][l][n][k][j][i] - (vsource1[k] + vsource2[k] + vsource3[k])/Crat);
							ReduceVelocity(sigmaa+sigmas, ds, &alpha);
							tempV[k] = Crat * alpha * SIGN(miuz);
						}/* End k */
						
						/* first, calculate the advection part v(3J+I) */
						flux_AdvJ(Radr, dir, tempS,	tempSV,	ks, ke+1, ds, dt, tempAdv);
						/* second, calculate the flux due to vsource3 */
						flux_AdvJ(Radr, dir, vsource3, tempV3V,	ks, ke+1, ds, dt, fluxsource3);
						/* Third, calculate flux due to co-moving miu */
						flux_AdvJ(Radr, dir, tempimu,	tempV,	ks, ke+1, ds, dt, flux);
						
						/* Now save the flux difference. Note that flux_advJ only calculates the interface values, not the actual flux */
						for(k=ks; k<=ke; k++){
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempSV[k+1]+tempSV[k])	* tempAdv[k+1] - (tempSV[k] + tempSV[k-1]) * tempAdv[k]) * dtods;
							Divi[ifr][l][n][k][j][i]+= 0.5 * ((tempV3V[k+1]+tempV3V[k]) * fluxsource3[k+1] - (tempV3V[k] + tempV3V[k-1]) * fluxsource3[k]) * dtods;
							Divi[ifr][l][n][k][j][i]+= Crat * (flux[k+1] - flux[k]) * dtods;						
						}/* End k */	
				
					} /* Finish i */
				}/* Finish j */

			}/* end n */
		}/* end l */
	}/* end ifr */

	/* Now we have flux, now add the source terms due to absorption and scattering opacity seperately */ 
	
	
	/*--------------------------------------------------------------------*/
	/* The energy change due to T^4/4pi - I and T^4/4pi - J term is already calculated in pRG->heatcool */
	/* All the other terms in the transfer equations are solved together
	 * The independent variable for the matrix solver is actually weight * I
	 * Ma[n] is the diagonal part for each I
	 * Mc[n] is the for the part ndot v (sigma_a+sigma_s)* 3J 
	 * Mb[n] is for the other momentums terms *
	 */
	for(ifr=0; ifr<pRG->nf; ifr++){
		for(l=0; l<pRG->noct; l++){
			for(n=0; n<pRG->nang; n++){
				for(k=ks; k<=ke; k++){
					for(j=js; j<=je; j++){
						for(i=is; i<=ie; i++){
							/* first construct the Right hand side */
							/* We need an array for different angles for each cell. The index order of RHS is different */
							Mi = l*(pRG->nang)+n;
							RHS[k][j][i][Mi] = pRG->imu[ifr][l][n][k][j][i]- Divi[ifr][l][n][k][j][i] 
											+ pRG->heatcool[ifr][l][n][k][j][i];
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

							vx = pG->U[k+offset][j+offset][i+offset].M1 / pG->U[k+offset][j+offset][i+offset].d;
							vy = pG->U[k+offset][j+offset][i+offset].M2 / pG->U[k+offset][j+offset][i+offset].d;
							vz = pG->U[k+offset][j+offset][i+offset].M3 / pG->U[k+offset][j+offset][i+offset].d;
							AngleV = FullAngleV[ifr][l][n][k][j][i];
							AngleV2 = FullAngleV2[ifr][l][n][k][j][i];
							
							Ma[k][j][i][Mi] = (1.0 + dt * (sigmas * Crat - (sigmaa + sigmas) * AngleV))/pRG->wmu[n][k][j][i];
							Mc[k][j][i][Mi] = -dt * (sigmas * Crat + 3.0 * AngleV * (sigmas + sigmaa));
							Mb[k][j][i][Mi] = dt * 2.0 * sigmas * AngleV 
											+ dt * (sigmaa - sigmas) * (vx * vx + vy * vy + vz * vz + AngleV2)/Crat; 
						
						}/* end i */
					}/* end j */
				}/* end k */
			}/* end n */
		}/* end l */
		
		
		/* solve the first matrix */
		for(k=ks; k<=ke; k++)
			for(j=js; j<=je; j++)
				for(i=is; i<=ie; i++){
					SpecialMatrix3(Ma[k][j][i], Mb[k][j][i], Mc[k][j][i], Md, RHS[k][j][i], lN, UN,  pRG->nang*pRG->noct);
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
	if(FullAngleV2 != NULL) free_6d_array(FullAngleV2);	
	if(FullVsource3 != NULL) free_6d_array(FullVsource3);
	
	if(flux != NULL) free_1d_array(flux);
	if(fluxsource3 != NULL) free_1d_array(fluxsource3);
	if(tempV3V != NULL) free_1d_array(tempV3V);
	if(tempimu != NULL) free_1d_array(tempimu);
	if(tempS != NULL) free_1d_array(tempS);
	if(tempSV != NULL) free_1d_array(tempSV);
	if(vsource1 != NULL) free_1d_array(vsource1);
	if(vsource2 != NULL) free_1d_array(vsource2);
	if(vsource3 != NULL) free_1d_array(vsource3);
	if(tempV != NULL) free_1d_array(tempV);
	if(tempAdv != NULL) free_1d_array(tempAdv);
	if(RHS != NULL) free_4d_array(RHS);
	if(Ma != NULL) free_4d_array(Ma);
	if(Mb != NULL) free_4d_array(Mb);
	if(Mc != NULL) free_4d_array(Mc);
	if(Md != NULL) free_1d_array(Md);	
	if(lN != NULL) free_2d_array(lN);
	if(UN != NULL) free_2d_array(UN);
  return;
}


void fullRT_3d_init(RadGridS *pRG)
{

	
	int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
	int nfr = pRG->nf, noct = pRG->noct, nang = pRG->nang;
	int nmax;


	nmax = MAX(nx1,nx2);
	nmax = MAX(nmax,nx3);

	if ((flux = (Real *)calloc_1d_array( nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((fluxsource3 = (Real *)calloc_1d_array( nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempV3V = (Real *)calloc_1d_array( nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Divi = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((FullAngleV = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((FullAngleV2 = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((FullVsource3 = (Real ******)calloc_6d_array(nfr, noct, nang, nx3+2*Radghost, nx2+2*Radghost, nx1+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempimu = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempS = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempSV = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((vsource1 = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((vsource2 = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((vsource3 = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempV = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((tempAdv = (Real *)calloc_1d_array(nmax+2*Radghost, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((RHS = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Ma = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Mb = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Mc = (Real ****)calloc_4d_array(nx3+2*Radghost, nx2+2*Radghost,nx1+2*Radghost,noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((Md = (Real *)calloc_1d_array(noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	
	if ((lN = (Real **)calloc_2d_array(noct*nang, noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	
	if ((UN = (Real **)calloc_2d_array(noct*nang, noct*nang, sizeof(Real))) == NULL)
		goto on_error;
	

	return;

	on_error:
  	fullRT_3d_destruct();
  	ath_error("[fullRT_3d_init]: Error allocating memory\n");
  	return;

}
#endif /* FULL_RADIATION_TRANSFER */
