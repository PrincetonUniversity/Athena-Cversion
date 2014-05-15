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
/* The r[3] and dir are needed for cylindrical coordinate */
/* dir label the direction along which the flux is calculated */
/* r[3] is coordinate center of each cell */

void flux_PLM(Real r[3]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[3], Real imhalf[1])
{
/* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
/* imup[0:2] is i-2, i-1, i*/

  Real dqi0, delq1, delq2;
  Real *imup;
  Real distance;
#ifdef CYLINDRICAL
  Real *pr;
  pr = &(r[2]);
  Real geom1 = 1.0, geom2 = 1.0;
#endif

  Real geom3 = 1.0, geom4 = 1.0; /* geometric weighting factor */

  imup = &(imu[2]);


/* The upwind slope */
  delq1 = (imup[0] - imup[-1]) / ds;
  delq2 = (imup[-1] - imup[-2]) / ds;
#ifdef CYLINDRICAL
/* Only need to apply the weighting factor when we want to calculate
 * flux along radial direction */
  if(dir == 1){
    geom1 = 1.0/(1.0 - ds * ds/(12.0 * pr[0] * pr[-1]));
    geom2 = 1.0/(1.0 - ds * ds/(12.0 * pr[-1] * pr[-2]));

    delq1 *= geom1;
    delq2 *= geom2;
  }
#endif


  if(delq1 * delq2 > 0.0){
    dqi0 = 2.0 * delq1 * delq2 / (delq1 + delq2);
  }
  else{
    dqi0 = 0.0;
  }

#ifdef CYLINDRICAL
/* Take into account the curvature effect for time averaged interface state */
/* See eq. 64 of skinner & ostriker 2010 */
  if(dir == 1){
    geom3 = 1.0 - ds /(6.0 * pr[-1]);
    geom4 = 1.0 - vel * dt/(6.0*(0.5 * (pr[-1] + pr[0]) - 0.5 * vel * dt));
  }
#endif

  distance = ds * geom3;
  distance -= ((vel * dt) * geom4);

#ifdef CYLINDRICAL
/* Take into account the curvature effect for time averaged interface state */
/*      if(dir == 1){

        if(distance < 0.0){
        ath_error("[FullRT_flux plm]: CFL condition should be satisfied for volume centered distance %e\n",distance);

        }
        }
*/
#endif

/* The upwind flux */
  imhalf[0] = imup[-1] + distance * (dqi0/2.0);

  return;
}



/* calculate the specific intensity at i-1/2 and i+1/2
   with third-order accurate Colella & Woodward slope, Stone & Mihalas (1992)
*/

void flux_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real dt, const Real ds, const Real vel, Real imu[5], Real imhalf[1])
{
/* for each ray, we always use the upwind specific intensity , therefore velocity is always positive */
/* imup[0:4] is i-3, i-2, i-1, i, i+1, from downwind to upwind */
  Real iLeft, iRight;
  Real xsi;
  Real imup;
  int n;

#ifdef CYLINDRICAL
  Real a6, gcurve, adiff, radius, radiusi, betacurve;

#endif

/* First, we need to calculate left and right state for zone i and i-1, in order to calculate half state i+1/2 */

  lrstate_PPM(r, dir, ds, imu, &(iLeft), &(iRight));

/* iLeft[0] and iRight[0] are for cell i-1 and iLeft[1] and iRight[1] are for cell i */

/* The velocity is always positive, use equation 6 in SM to calculate the flux */
  xsi = vel * dt / ds;

/* The upwind flux */
  for (n=0; n<1; n++){
    imup = imu[n+2];
/*      imhalf[n] = iRight[n] + xsi * (imup[n] - iRight[n]) + xsi * (1.0 - xsi) * (2.0 * imup[n] - iRight[n] - iLeft[n]);
 */
/*      Itemp = 6.0 * imup[n] - 3.0 * (iLeft[n] + iRight[n]);
        imhalf[n] = iRight[n] - 0.5 * xsi * ((iRight[n] - iLeft[n]) - (1.0 - 2.0 * xsi/3.0) * Itemp);
*/

#ifdef CYLINDRICAL
/* Only need to add the curvature factor along radial direction */
/* Following eq. 66a of Skinner & Ostriker 2010 */
/* If gcurve and betacurve are all zero, this formula is reduced to the normal one */
    if(dir == 1){
      radius = r[n+2];
      radiusi = 0.5 * (r[n+2] + r[n+3]);
      adiff = iRight - iLeft;
      gcurve = ds/(6.0 * radius);
      a6 = 6.0 * (imup - iLeft - 0.5 * adiff * (1.0 + gcurve));
      betacurve = vel * dt / (6.0 * (radiusi - 0.5 * vel * dt));
      imhalf[n] = iRight - 0.5 * xsi * (adiff - (1.0 - 2.0 * xsi / 3.0) * a6)
        + 0.5 * xsi * (adiff - (1.0 - xsi) * a6) * betacurve;

    }else {
      imhalf[n] = (xsi - 1.0) * (xsi - 1.0) * iRight + xsi * (xsi - 1.0) * iLeft + xsi * (3.0 - 2.0 * xsi) * imup;
    }



#else
    imhalf[n] = (xsi - 1.0) * (xsi - 1.0) * iRight + xsi * (xsi - 1.0) * iLeft + xsi * (3.0 - 2.0 * xsi) * imup;
#endif

  }/* end n=0 and 1 */

  return;
}

/* Get left and right state with PPM for cell i, and i-1 */
/* We need both left and right state */
/* This assume uniform spacing */
/* pG and dir are needed for cylindrical coordinate cases */
/* dir is the direction along which the flux is calculated */
void lrstate_PPM(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real ds __attribute__((unused)), Real imu[5], Real iLeft[1], Real iRight[1])
{
/* imu[0:4] i-3, i-2, i-1, i, i+1 */
  Real Ihalf0, Ihalf1, lim_slope, IL, IR;
  Real I03, I02, I01, I0, I1;
  Real qa, qb, qc, dI, dIc, dIl, dIr, dIlim;


/* first, calculate the half state */
/*  i-3    | i-2   |    i -1     | i    | i + 1  */
/*                half[0]    half[1]
 *            dI[0]      dI[1]      dI[2]
 */


/* To calculate the values at i-3/2 and i-1/2 */

/* Do not use pointer, save variables */
  I03 = imu[0]; I02 = imu[1]; I01 = imu[2]; I0 = imu[3]; I1 = imu[4];

/*----------------------------------------------*/

/* First, calculate Ihalf0 */
  Ihalf0 = (7.0 * (I02 + I01) - (I0 + I03)) / 12.0;
  dIc = 3.0 * (I02 - 2.0 * Ihalf0 + I01);
  dIl = (I03 - 2.0 * I02 + I01);
  dIr = (I02 - 2.0 * I01 + I0);
  dIlim = 0.0;

  lim_slope = MIN(fabs(dIl),fabs(dIr));

  if((dIc > 0.0) && (dIl > 0.0) && (dIr > 0.0)){
    dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
  }

  if((dIc < 0.0) && (dIl < 0.0) && (dIr < 0.0)){
    dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
  }

  Ihalf0 = 0.5 * ((I02 + I01) - dIlim/3.0);

/*----------------------------------------------*/
/* Now calculate Ihalf1 */


  Ihalf1 = (7.0 * (I01 + I0) - (I1 + I02)) / 12.0;
  dIc = 3.0 * (I01 - 2.0 * Ihalf1 + I0);
  dIl = (I02 - 2.0 * I01 + I0);
  dIr = (I01 - 2.0 * I0 + I1);
  dIlim = 0.0;

  lim_slope = MIN(fabs(dIl),fabs(dIr));

  if((dIc > 0.0) && (dIl > 0.0) && (dIr > 0.0)){
    dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
  }

  if((dIc < 0.0) && (dIl < 0.0) && (dIr < 0.0)){
    dIlim = SIGN(dIc) * MIN(1.25 * lim_slope, fabs(dIc));
  }

  Ihalf1 = 0.5 * ((I01 + I0) - dIlim/3.0);

/* Now calculate the left and right states */


  IL = Ihalf0;
  IR = Ihalf1;

  qa = (IR - I01) * (I01 - IL);
  qb = (I02 - I01) * (I01 - I0);

  if((qa <= 0.0) && (qb <= 0.0)){
    qc = 6.0 * (I01 - 0.5 * (IL + IR));
    dI = -2.0 * qc;
    dIc = I02 - 2.0 * I01 + I0;
    dIl = I03 - 2.0 * I02 + I01;
    dIr = I01 - 2.0 * I0 + I1;
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
      IL = I01;
      IR = I01;
    }else{
      IL = I01 + (IL - I01) * dIlim/dI;
      IR = I01 + (IR - I01) * dIlim/dI;

    }
  }


/* constrain the values */
  qa = (IR - I01) * (I01 - IL);
  qb = IR - IL;
  qc = 6.0 * (I01 - 0.5 * (IL + IR));

  if(qa <= 0.0){
    IL = I01;
    IR = I01;
  }else if((qb * qc) > (qb * qb)){
    IL = 3.0 * I01 - 2.0 * IR;
  } else if ((qb * qc) < -(qb * qb)){
    IR = 3.0 * I01 - 2.0 * IL;
  }

  IL = MAX(MIN(I01,I02),IL);
  IL = MIN(MAX(I01,I02),IL);

  IR = MAX(MIN(I01,I0),IR);
  IR = MIN(MAX(I01,I0),IR);

/* Now set left and right state */
  iLeft[0] = IL;
  iRight[0] = IR;


/*      for(n=0; n<2; n++){
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
*/
/*      for(n=0; n<1; n++){
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

*/
/* constrain the values */
/*              qa = (IR - imup[n]) * (imup[n] - IL);
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


                iLeft[n] = IL;
                iRight[n] = IR;

                }
*/
/* Finish n */




  return;

}

/* Try a new ppm scheme adopted from lr_state_ppm */
void lrstate_PPM2(Real r[5]  __attribute__((unused)), const int dir __attribute__((unused)), const Real ds __attribute__((unused)), Real imu[5], Real iLeft[1], Real iRight[1])
{
/* imu[0:4] i-3, i-2, i-1, i, i+1 */
  Real Ihalf0, Ihalf1, IL, IR;
  Real I03, I02, I01, I0, I1;
  Real qa, qb, qc, dI, dIl, dIr, dIlim;
#ifdef CYLINDRICAL
  Real *pr;
  pr = &(r[3]);
#endif
  Real gcurve = 0.0;


/* Do not use pointer, save variables */
  I03 = imu[0]; I02 = imu[1]; I01 = imu[2]; I0 = imu[3]; I1 = imu[4];

/* First, calculate the left state between i-2 and i-1 */
/* We need gradient for cell i-2 and i-1 */
#ifdef CYLINDRICAL
  if(dir == 1){
    dIl = (pr[-2] * I02 - pr[-3] * I03) / (0.5 * (pr[-2] + pr[-3]));
    dIr = (pr[-1] * I01 - pr[-2] * I02) / (0.5 * (pr[-1] + pr[-2]));
  }
  else
#endif
    {
      dIl = I02 - I03;
      dIr = I01 - I02;
    }

  if(dIl * dIr > 0.0)
    dIlim = 2.0 * dIl * dIr/(dIl + dIr);
  else {
    dIlim = 0.0;
  }

/*--------------------*/
#ifdef CYLINDRICAL
  if(dir == 1){
    dIl = (pr[-1] * I01 - pr[-2] * I02) / (0.5 * (pr[-1] + pr[-2]));
    dIr = (pr[0] * I0 - pr[-1] * I01) / (0.5 * (pr[0] + pr[-1]));
  }
  else
#endif
    {
      dIl = I01 - I02;
      dIr = I0 - I01;
    }

  if(dIl * dIr > 0.0)
    dI = 2.0 * dIl * dIr/(dIl + dIr);
  else {
    dI = 0.0;
  }

#ifdef CYLINDRICAL
  if(dir == 1){
    Ihalf0 = (0.5 * (I01 * pr[-1] + I02 * pr[-2]) - (dI * pr[-1] - dIlim * pr[-2])/6.0)/(0.5*(pr[-1]+pr[-2]));
  }
  else
#endif
    {
      Ihalf0 = 0.5 * (I01 + I02) - (dI - dIlim)/6.0;

    }
/*-------------------------------------------------*/
/* Second, calculate the Right state between i and i-1 */
/* We need gradient for cell i and i-1 */
#ifdef CYLINDRICAL
  if(dir == 1){
    dIl = (pr[-1] * I01 - pr[-2] * I02) / (0.5 * (pr[-1] + pr[-2]));
    dIr = (pr[0] * I0 - pr[-1] * I01) / (0.5 * (pr[0] + pr[-1]));
  }
  else
#endif
    {
      dIl = I01 - I02;
      dIr = I0 - I01;
    }


  if(dIl * dIr > 0.0)
    dIlim = 2.0 * dIl * dIr/(dIl + dIr);
  else {
    dIlim = 0.0;
  }


/*-------------------*/
#ifdef CYLINDRICAL
  if(dir == 1){
    dIl = (pr[0] * I0 - pr[-1] * I01) / (0.5 * (pr[0] + pr[-1]));
    dIr = (pr[1] * I1 - pr[0] * I0) / (0.5 * (pr[1] + pr[0]));
  }
  else
#endif
    {
      dIl = I0 - I01;
      dIr = I1 - I0;
    }

  if(dIl * dIr > 0.0)
    dI = 2.0 * dIl * dIr/(dIl + dIr);
  else {
    dI = 0.0;
  }



#ifdef CYLINDRICAL
  if(dir == 1){
    Ihalf1 = (0.5 * (I0 * pr[0] + I01 * pr[-1]) - (dI * pr[0] - dIlim * pr[-1])/6.0)/(0.5*(pr[0]+pr[-1]));
  }
  else
#endif
    {
      Ihalf1 = 0.5 * (I0 + I01) - (dI - dIlim)/6.0;
    }


/*================================================*/
/* Now make monotonic contrain */


  IL = Ihalf0;
  IR = Ihalf1;

#ifdef CYLINDRICAL
  if (dir == 1) gcurve = ds/(6.0*pr[-1]);
#endif


  qa = (IR - I01) * (I01 - IL);
  qb = IR - IL;
  qc = 6.0 * (I01 - 0.5 * (IL * (1.0 - gcurve) + IR * (1.0 + gcurve)));

  if(qa <= 0.0){
    IL = I01;
    IR = I01;
  }else if((qb * qc) > (qb * qb)){
    IL = (6.0 * I01 - IR * (4.0 + 3.0 * gcurve))/(2.0 - 3.0 * gcurve);
  } else if ((qb * qc) < -(qb * qb)){
    IR = (6.0 * I01 - IL * (4.0 - 3.0 * gcurve))/(2.0 + 3.0 * gcurve);
  }

  IL = MAX(MIN(I01,I02),IL);
  IL = MIN(MAX(I01,I02),IL);

  IR = MAX(MIN(I01,I0),IR);
  IR = MIN(MAX(I01,I0),IR);

/* Now set left and right state */
  iLeft[0] = IL;
  iRight[0] = IR;



}

/* This function is used to calculate the interface value, for the case the direction
 * of velocity is uncertain. So we need to determine the upwind direction first.
 * The returned value is the interface only */

void flux_AdvJ(const int nf, const int N, Real *r  __attribute__((unused)), const int dir __attribute__((unused)), Real ***tempJ, Real ***tempV, int nstart,int nend,Real ds, Real dt, Real ***tempAdv)
{
  int i, j, n, ifr;
  Real vel, meanJ, dtods;
  Real Jarray[5];
  Real tempr[5];
  dtods = dt/ds;


  for(i=nstart; i<=nend; i++){
    for(ifr=0; ifr<nf; ifr++){
      for(n=0; n<N; n++){
        vel = 0.5 * (tempV[i-1][ifr][n] + tempV[i][ifr][n]);
        if(vel > 0.0){
#ifdef SECOND_RAD_ORDER
          for(j=0; j<3; j++){
            Jarray[j] = tempJ[i-2+j][ifr][n];
#ifdef CYLINDRICAL
            tempr[j] = r[i-2+j];
#endif
          }

          flux_PLM(tempr, dir, dt, ds, vel, Jarray, &(meanJ));
#else

          for(j=0; j<5; j++){
            Jarray[j] = tempJ[i-3+j][ifr][n];
#ifdef CYLINDRICAL
            tempr[j] = r[i-3+j];
#endif

          }

          flux_PPM(tempr, dir, dt, ds, vel, Jarray, &(meanJ));

#endif
        }/* End if vel >0 */
        else{
#ifdef SECOND_RAD_ORDER
          for(j=0; j<3; j++){
            Jarray[j] = tempJ[i+1-j][ifr][n];
#ifdef CYLINDRICAL
            tempr[j] = r[i+1-j];
#endif

          }

          flux_PLM(tempr, dir, dt, ds, -vel, Jarray, &(meanJ));
#else

          for(j=0; j<5; j++){
            Jarray[j] = tempJ[i+2-j][ifr][n];
#ifdef CYLINDRICAL
            tempr[j] = r[i+2-j];
#endif

          }

          flux_PPM(tempr, dir, dt, ds, -vel, Jarray, &(meanJ));

#endif
        }
        tempAdv[i][ifr][n] = meanJ;
      }/* end n */
    }/* End if */
  }/* end i */



}

#endif /* RADIATION_TRANSFER */


