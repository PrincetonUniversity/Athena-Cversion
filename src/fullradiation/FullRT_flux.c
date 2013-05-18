#include "../copyright.h"
/*==============================================================================
 * FILE: FULLRT_flux.c
 *
 * PURPOSE: contain PLM and PPM functions to calculate flux 
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


/* calculate the specific intensity at i-1/2 and i+1/2 
   with second-order accurate van Leer slope, Stone & Mihalas (1992)
*/

void flux_PLM(const Real dt, const Real ds, const Real vel, Real imu[3], Real imhalf[1])
{
/* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
/* imup[0:2] is i-2, i-1, i*/
	
	Real dqi0, delq1, delq2;
	Real *imup;

	imup = &(imu[2]);


	/* The upwind slope */
	delq1 = (imup[0] - imup[-1]) / ds;
	delq2 = (imup[-1] - imup[-2]) / ds;
	if(delq1 * delq2 > 0.0){
		dqi0 = 2.0 * delq1 * delq2 / (delq1 + delq2);
	} 
	else{
		dqi0 = 0.0;
	}

	/* The upwind flux */
	imhalf[0] = imup[-1] + (ds - vel * dt) * (dqi0/2.0);
	
	return;
}



/* calculate the specific intensity at i-1/2 and i+1/2 
   with third-order accurate Colella & Woodward slope, Stone & Mihalas (1992)
*/

void flux_PPM(const Real dt, const Real ds, const Real vel, Real imu[5], Real imhalf[1])
{
/* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
/* imup[0:4] is i-3, i-2, i-1, i, i+1, from downwind to upwind */
	Real iLeft[1], iRight[1];
	Real xsi;
	Real *imup;
	Real Itemp;
	int n;
	/* First, we need to calculate left and right state for zone i and i-1, in order to calculate half state i+1/2 */

	lrstate_PPM(imu, &(iLeft[0]), &(iRight[0]));

	/* iLeft[0] and iRight[0] are for cell i-1 and iLeft[1] and iRight[1] are for cell i */

	/* The velocity is always positive, use equation 6 in SM to calculate the flux */
	xsi = vel * dt / ds;

	/* The upwind flux */
	for (n=0; n<1; n++){
		imup = &(imu[n+2]);
	/*	imhalf[n] = iRight[n] + xsi * (imup[n] - iRight[n]) + xsi * (1.0 - xsi) * (2.0 * imup[n] - iRight[n] - iLeft[n]);
	*/
		Itemp = 6.0 * imup[n] - 3.0 * (iLeft[n] + iRight[n]);
		imhalf[n] = iRight[n] - 0.5 * xsi * ((iRight[n] - iLeft[n]) - (1.0 - 2.0 * xsi/3.0) * Itemp);
		
	 }	

	return;
}

	/* Get left and right state with PPM for cell i, and i-1 */
	/* We need both left and right state */
	/* This assume uniform spacing */
void lrstate_PPM(Real imu[5], Real iLeft[1], Real iRight[1])
{
	/* imu[0:4] i-3, i-2, i-1, i, i+1 */
	Real *imup;
	Real Ihalf[2],lim_slope, IL, IR;
	Real qa, qb, qc, dI, dIc, dIl, dIr, dIlim;
	int n;
	
	/* first, calculate the half state */
	/*  i-3    | i-2   |    i -1     | i    | i + 1  */
	/*				  half[0]    half[1]
	 *	      dI[0]      dI[1]      dI[2]
	 */
	
	
	/* To calculate the values at i-3/2 and i-1/2 */
	/* We first calculate the gradient (equation 1.7 and 1.8 of CW) */
	for(n=0; n<2; n++){
		imup = &(imu[n+1]);
		Ihalf[n] = (7.0 * (imup[n] + imup[n+1]) - (imup[n+2] + imup[n-1])) / 12.0;
		dIc = 3.0 * (imup[n] - 2.0 * Ihalf[n] + imup[n+1]);
		dIl = (imup[n-1] - 2.0 * imup[n] + imup[n+1]);
		dIr = (imup[n] - 2.0 * imup[n+1] + imup[n+2]);
		dIlim = 0.0;
		lim_slope = MIN(fabs(dIl),fabs(dIr));
		
		if((dIc > 0.0) && (dIl > 0.0) && (dIr > 0.0)){
			dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
		}
		
		if((dIc < 0.0) && (dIl < 0.0) && (dIr < 0.0)){
			dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
		}
		
		Ihalf[n] = 0.5 * ((imup[n] + imup[n+1]) - dIlim/3.0);
	}
	
	for(n=0; n<1; n++){
		imup = &(imu[n+2]);
		
		IL = Ihalf[n];
		IR = Ihalf[n+1];
		
		qa = (IR - imup[n]) * (imup[n] - IL);
		qb = (imup[n-1] - imup[n]) * (imup[n] - imup[n+1]);
		
		if((qa <= 0.0) && (qb <= 0.0)){
			qc = 6.0 * (imup[n] - 0.5 * (IL + IR));
			dI = -2.0 * qc;
			dIc = imup[n-1] - 2.0 * imup[n] + imup[n+1];
			dIl = imup[n-2] - 2.0 * imup[n-1] + imup[n];
			dIr = imup[n] - 2.0 * imup[n+1] + imup[n+2];
			dIlim = 0.0;
			lim_slope = MIN(fabs(dIl), fabs(dIr));
			lim_slope = MIN(fabs(dIc),lim_slope);
			
			if((dIc > 0.0) && (dIl > 0.0) && (dIr > 0.0) && (dI > 0.0)){
				dIlim = SIGN(dI) * MIN(1.25 * lim_slope, fabs(dI));
				
			}
			
			if((dIc < 0.0) && (dIl < 0.0) && (dIr < 0.0) && (dI < 0.0)){
				dIlim = SIGN(dI) * MIN(1.25 * lim_slope, fabs(dI));
				
			}
			
			if(dI == 0.0){
				IL = imup[n];
				IR = imup[n];
			}else{
				IL = imup[n] + (IL - imup[n]) * dIlim/dI;
				IR = imup[n] + (IR - imup[n]) * dIlim/dI;
				
			}
		}
		
		
		/* constrain the values */
		qa = (IR - imup[n]) * (imup[n] - IL);
		qb = IR - IL;
		qc = 6.0 * (imup[n] - 0.5 * (IL + IR));
		
		if(qa <= 0.0){
			IL = imup[n];
			IR = imup[n];
		}else if((qb * qc) > (qb * qb)){
			IL = 3.0 * imup[n] - 2.0 * IR;
		} else if ((qb * qc) < -(qb * qb)){
			IR = 3.0 * imup[n] - 2.0 * IL;
		}
		
		IL = MAX(MIN(imup[n],imup[n-1]),IL);
		IL = MIN(MAX(imup[n],imup[n-1]),IL);
		
		IR = MAX(MIN(imup[n],imup[n+1]),IR);
		IR = MIN(MAX(imup[n],imup[n+1]),IR);
		
		/* Now set left and right state */
		iLeft[n] = IL;
		iRight[n] = IR;
		
	}/* Finish n */
	
	
	
	
	return;
	
}


#endif /* RADIATION_TRANSFER */


