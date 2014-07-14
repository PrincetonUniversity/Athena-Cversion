#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_fullrad.c
 * copy from hydro_to_rad of the radiation-transfer
 * PURPOSE:  Calculate the opacity in the radiation grid and
 *                       calculate the gas temperature for later use.
 *
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   hydro_to_rad()
 *   rad_to_hydro()
 *============================================================================*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER


void GetSource(const Real dt, const Real Tnew, const Real SigmaB, const Real SigmaI, const Real Jsource, const Real I0, Real *heatcool);
void GetTnew(const Real dt,const Real d,const Real Tgas, const Real J0, const Real Sigma[4], Real *heatcool, Real *Tnew);

/*----------------------------------------------------------------------------*/
/* hydro_to_rad:  */

void hydro_to_fullrad(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k,ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;
  Real d, etherm;


/* Assumes ghost zone conserved variables have been set by
 * bvals routines.  These values are used to set B, chi, eps,
 * etc. so loops include RadGrid ghost zones*/
  if (pG->Nx[0] > 1) {
    ioff = nghost - Radghost;
    il -= Radghost;
    iu += Radghost;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - Radghost;
    jl -= Radghost;
    ju += Radghost;
  } else joff = 0;
  if (pG->Nx[2] > 1) {
    koff = nghost - Radghost;
    kl -= Radghost;
    ku += Radghost;
  } else koff = 0;

/* Compute radiation variables from conserved variables */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
        ig = i + ioff;

/* ------------------------------------*/
/* First, update the opacity */
        for(ifr=0; ifr<nf; ifr++) {
          get_full_opacity(pG,ifr,ig,jg,kg,&(pRG->R[k][j][i][ifr].Sigma[0]));
        }

/*-------------------------------------------*/


/* Compute gas temperature and store for later use */
        d = pG->U[kg][jg][ig].d;
        etherm = pG->U[kg][jg][ig].E - (0.5/d) * ( SQR(pG->U[kg][jg][ig].M1) +
                                                   SQR(pG->U[kg][jg][ig].M2) + SQR(pG->U[kg][jg][ig].M3) );
#ifdef MHD
        etherm -= 0.5 * (SQR(pG->U[kg][jg][ig].B1c) + SQR(pG->U[kg][jg][ig].B2c) +
                         SQR(pG->U[kg][jg][ig].B3c) );
#endif
        pG->tgas[kg][jg][ig] = MAX(etherm * Gamma_1 / (d * R_ideal),0.0);


      }/* end i */
    }/* end j */
  }/* end k */






  return;
}


/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to absorption opacity */


void RadAsource(const int i, const int j, const int k, const int N, const Real Tcoef, const Real Tnew, RadGridS *pRG, GridS *pG,  Real **Coefn, Real **sol)
{

  Real Fr0[3], Fr[3];
  Real weight, rho, Vel[3], Msource[3], DeltaKin[3], Ersource;
  Real Radheat, Frsource[3];

  int nDim, l, n, ifr;

  nDim = 1;
  for (l=1; l<3; l++) if (pRG->Nx[l]>1) nDim++;

  int koff, joff, ioff;
  ioff = nghost - Radghost;

  if(nDim > 1)
    joff = ioff;
  else
    joff = 0;

  if(nDim > 2)
    koff = ioff;
  else
    koff = 0;
  
    Radheat = 0.0;
    for(l=0; l<3; l++)
      Frsource[l] = 0.0;




/* loop through all l and n, nelements = nang * noct */
  for(ifr=0; ifr<pRG->nf; ifr++){
/* calculate Er and Fr used in the update */
    for(l=0; l<3; l++){
      Fr0[l] = 0.0;
      Fr[l] = 0.0;
    }



    for(n=0; n<N; n++){
      weight = pRG->wmu[k][j][i][n];

      for(l=0; l<3; l++){
        Fr0[l] += weight * Coefn[n][l] * sol[ifr][n];
        Fr[l] += weight * Coefn[n][l] * pRG->imu[k][j][i][ifr][n];
      }

    }/* end n */

/* multiple by 4 pi */

    for(l=0; l<3; l++){
      Fr0[l] *= 4.0 * PI;
      Fr[l] *= 4.0 * PI;
    }


    rho = pG->U[k+koff][j+joff][i+ioff].d;

    Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
    Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
    Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;



    for(l=0; l<3; l++){

/*   Msource[l] = dt * Prat * Sigma * CoFr[l];
 */
      Msource[l] = -Prat * (Fr[l] - Fr0[l])/Crat;

      DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]);


      if(DeltaKin[l] + 0.5 * rho * Vel[l] * Vel[l] < 0.0){
        Msource[l] = -rho * Vel[l];
        DeltaKin[l] = 0.5 * rho * Vel[l] * Vel[l];
      }


    }/* end l for three directions */

    Ersource = Tcoef * (Tnew - pG->tgas[k+koff][j+joff][i+ioff]) + DeltaKin[0] + DeltaKin[1] + DeltaKin[2];



/* Now put the energy and momentum source term back */
    Radheat += (pRG->wnu[ifr] * Ersource);

    for(l=0; l<3; l++)
      Frsource[l] += (pRG->wnu[ifr] * Msource[l]);


  }/* end ifr */
  
  pG->U[k+koff][j+joff][i+ioff].E += Radheat;
  pG->U[k+koff][j+joff][i+ioff].M1 += Frsource[0];
  pG->U[k+koff][j+joff][i+ioff].M2 += Frsource[1];
  pG->U[k+koff][j+joff][i+ioff].M3 += Frsource[2];
  
  
  
  return;

}



/* function to calculate the energy and momentum source terms from updated radiation quantities for a particular frequency band ifr*/
/* the updated sol already include the weight for each ray */
/* This function only calculate source terms due to scattering opacity */

/* moments of the radiation are already updated before entering this function */
/* So no need to calculate the moments again */
/* The old solution used is at sol while the updated solution is at pRG->imu */

void RadSsource(const int i, const int j, const int k, const int N, RadGridS *pRG, GridS *pG, Real **Coefn, Real **sol)
{

  Real Fr0[3], Fr[3];
  Real weight, rho, Vel[3], Msource[3], DeltaKin[3], Ersource;
  Real Frsource[3], Radheat;

  int nDim, l, n, ifr;

  nDim = 1;
  for (l=1; l<3; l++) if (pRG->Nx[l]>1) nDim++;

  int koff, joff, ioff;
  ioff = nghost - Radghost;

  if(nDim > 1)
    joff = ioff;
  else
    joff = 0;

  if(nDim > 2)
    koff = ioff;
  else
    koff = 0;

  for(l=0; l<3; l++)
    Frsource[l] = 0.0;
  
  Radheat = 0.0;



/* loop through all l and n, nelements = nang * noct */
  for(ifr=0; ifr<pRG->nf; ifr++){
/* calculate Er and Fr used in the update */
    for(l=0; l<3; l++){
      Fr0[l] = 0.0;
      Fr[l] = 0.0;
    }



    for(n=0; n<N; n++){
      weight = pRG->wmu[k][j][i][n];

      for(l=0; l<3; l++){
        Fr0[l] += weight * Coefn[n][l] * sol[ifr][n];
        Fr[l] += weight * Coefn[n][l] * pRG->imu[k][j][i][ifr][n];
      }

    }/* end n */

/* multiple by 4 pi */

    for(l=0; l<3; l++){
      Fr0[l] *= 4.0 * PI;
      Fr[l] *= 4.0 * PI;
    }


    rho = pG->U[k+koff][j+joff][i+ioff].d;

    Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
    Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
    Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;



    for(l=0; l<3; l++){

/*   Msource[l] = dt * Prat * Sigma * CoFr[l];
 */
      Msource[l] = -Prat * (Fr[l] - Fr0[l])/Crat;

      DeltaKin[l] = (SQR(Msource[l])*0.5/rho + Vel[l] * Msource[l]);


      if(DeltaKin[l] + 0.5 * rho * Vel[l] * Vel[l] < 0.0){
        Msource[l] = -rho * Vel[l];
        DeltaKin[l] = 0.5 * rho * Vel[l] * Vel[l];
      }


    }/* end l for three directions */

    Ersource = DeltaKin[0] + DeltaKin[1] + DeltaKin[2];



/* Now put the energy and momentum source term back */
    Radheat += (pRG->wnu[ifr] * Ersource);

    for(l=0; l<3; l++)
      Frsource[l] += (pRG->wnu[ifr] * Msource[l]);


  }/* end ifr */
  
  /* Add the energy and momentum to the gas */
  /* For scattering opacity, gas temperature should not change after this step */
  pG->U[k+koff][j+joff][i+ioff].E += Radheat;
  pG->U[k+koff][j+joff][i+ioff].M1 += Frsource[0];
  pG->U[k+koff][j+joff][i+ioff].M2 += Frsource[1];
  pG->U[k+koff][j+joff][i+ioff].M3 += Frsource[2];

}



/* Distribute the Compton scattering contribution to each rays in the lab frame */
/* There is no moment change to first order of v/c due to Compton scattering */
void ComptIntensity(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  GridS *pG = (pD->Grid);

  int i, j, k, ig, jg, kg;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int koff = 0, joff = 0, ioff = 0;
  int nDim;
  int nf = pRG->nf, nelements, noct, nang;
  int ifr, n;
  
  Real Tgas, Er, Ernew, Tr, Ersource;
  Real coefA, coefK, coefB,coef1,coef2,coef3,coef4;
  PrimS Wtemp;
  Real Jnew, Hnew[3], Knew[6], Vel[3];
  Real dtsigmas;
  Real miux, miuy, miuz;
  Real vdotn, vdotH, Hdotn, nnK, vnK, vsquar;
  
  
  noct = pRG->noct;
  nang = pRG->nang;
  nelements = noct * nang;




  nDim = 1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

  ioff = nghost - Radghost;

/* Including the ghost zones */

  if(nDim > 1){
    joff = nghost - Radghost;
  }

  if(nDim > 2){
    koff = nghost - Radghost;
  }
  
  

  for(k=kl; k<=ku; k++){
    for(j=jl; j<=ju; j++){
      for(i=il; i<=iu; i++){
        kg = k + koff;
        jg = j + joff;
        ig = i + ioff;
        
        Ersource = 0.0;
 
        
        Wtemp = Cons_to_Prim(&(pG->U[kg][jg][ig]));
        Tgas = Wtemp.P / (R_ideal * Wtemp.d);

        for(ifr=0; ifr<nf; ifr++){
          if(Comptflag){
              /* First, update zeroth moment J due to Compton process */
            
            Er = 4.0 * PI * pRG->R[k][j][i][ifr].J;
            dtsigmas = pG->dt * pRG->R[k][j][i][ifr].Sigma[2];

            Tr = sqrt(Er);
            Tr = sqrt(Tr);

            coefA = 4.0 * Crat * dtsigmas / (T_e/Tunit);
            coefK = (Gamma - 1.0) * Prat/(R_ideal * Wtemp.d);
            coefB = Tgas + coefK * Er;
            coef1 = coefA * coefK;
            coef2 = coefA;
            coef3 = 1.0 - coefA * coefB;
            coef4 = -Er;

            if(Tr < Tgas){
              Tr = rtsafe(Tcompton, Tr * (1.0 - 0.01), Tgas * (1.0 + 0.01), 1.e-10, coef1, coef2, coef3, coef4);
            }
            else{

              Tr = rtsafe(Tcompton, Tgas * (1.0 - 0.01), Tr * (1.0 + 0.01), 1.e-10, coef1, coef2, coef3, coef4);
            }

            Ernew = SQR(SQR(Tr));

            pRG->Ercompt[k][j][i][ifr] = (Ernew - Er)/(4.0*PI);
            Ersource += (pRG->wnu[ifr] * (Ernew - Er));
          
            Jnew = pRG->R[k][j][i][ifr].J + pRG->Ercompt[k][j][i][ifr];
            
            /* Use the actual velocity, not need to use guess velocity for Compton scattering */
            
            Vel[0] = pG->Velguess[k+koff][j+joff][i+ioff][0];
            Vel[1] = pG->Velguess[k+koff][j+joff][i+ioff][1];
            Vel[2] = pG->Velguess[k+koff][j+joff][i+ioff][2];
      
      

            
            /* update H due to bulk compton */
 /*           UpdateHcomp(&(Hnew[0]), pRG->R[k][j][i][ifr].H, dtsigmas, Vel);
  *         Do not do this now. There is no corresponding energy term for this bulk Compton momentum term */
  /*        for(n=0; n<3; n++)
                MomSource[n] += (pRG->wnu[ifr] * (Hnew[n] - pRG->R[k][j][i][ifr].H[n]));
    */
            for(n=0; n<3; n++){
                Hnew[n] = pRG->R[k][j][i][ifr].H[n];
            }
            
              
            /* update K due to Compton terms */
            UpdateKcomp(&(Knew[0]), pRG->R[k][j][i][ifr].K, Hnew, Vel, Jnew, dtsigmas, pRG->Ercompt[k][j][i][ifr]);
            
            
            /* With the updated J, H and K, we can update intensity due to Compton terms */
            for(n=0; n<nelements; n++){

#ifdef CYLINDRICAL
                miux = pRG->Rphimu[k][j][i][n][0];
                miuy = pRG->Rphimu[k][j][i][n][1];
                miuz = pRG->Rphimu[k][j][i][n][2];
#else
                miux = pRG->mu[k][j][i][n][0];
                miuy = pRG->mu[k][j][i][n][1];
                miuz = pRG->mu[k][j][i][n][2];
#endif
        
                vdotn = (Vel[0] * miux + Vel[1] * miuy + Vel[2] * miuz);
                vdotH = (Vel[0] * Hnew[0] + Vel[1] * Hnew[1] + Vel[2] * Hnew[2]);
                Hdotn = (Hnew[0] * miux + Hnew[1] * miuy + Hnew[2] * miuz);
                nnK = miux * miux * Knew[0] + 2.0 * miux * miuy * Knew[1] + miuy * miuy * Knew[2] + 2.0 * miux * miuz * Knew[3] + 2.0 * miuy * miuz * Knew[4] + miuz * miuz * Knew[5];
                vnK = Vel[0] * miux * Knew[0] + (Vel[0] * miuy + Vel[1] * miux) * Knew[1] + Vel[1] * miuy * Knew[2] + (Vel[0] * miuz + Vel[2] * miux) * Knew[3] + (Vel[1] * miuz + Vel[2] * miuy) * Knew[4] + Vel[2] * miuz * Knew[5];
              
      /*  Do not include (v/c)^2 terms right now */
          /*
                pRG->imu[k][j][i][ifr][Mi] = pRG->imu[k][j][i][ifr][Mi] + pRG->Ercompt[k][j][i][ifr] + dtsigmas * (-0.75 * vdotn + 5.25 * vdotn * vdotn/Crat - 1.75 * vsquar/Crat) * Jnew
                                              + dtsigmas * (0.5 * vdotH - 1.5 * vdotn * Hdotn - 7.5 * vdotn * vdotn * Hdotn/Crat + 1.5 * vsquar * Hdotn/Crat - 3.0 * vdotn * vdotH/Crat)
                                              + dtsigmas * (3.75 * vdotn * nnK - 1.5 * vnK);
          */
        
                /* Include the J and K terms, the first moment is not zero numerically, although it should be */
                pRG->ComptI[k][j][i][ifr][n] = pRG->Ercompt[k][j][i][ifr] + 0.0 * dtsigmas * (0.5 * vdotH - 1.5 * vdotn * Hdotn)
                                          + 0.0 * (dtsigmas * (3.75 * vdotn * nnK - 1.5 * vnK) - dtsigmas * 0.75 * vdotn * Jnew);


            }/* End nelements */
            
            
          }/* End comptflag */
          else{
            pRG->Ercompt[k][j][i][ifr] = 0.0;
            Ersource = 0.0;
            for(n=0; n<nelements; n++){
              pRG->ComptI[k][j][i][ifr][n] = 0.0;
            }
          }/* Comptflag */

        }/* End ifr */
  
        /* Add Compton energy source to gas end pressure source */
        pG->U[kg][jg][ig].E += -Prat * Ersource;


      }/* end i */
    }/* End j */
  }/* End k */
  
  
  /* Update boundary condition for gas quantity, as gas energy is changed */
  bvals_mhd(pD);
  /* Because gas temperature is updated here, we need to update opacity and tgas */
  /* scattering opaicty usually does not depends on temperature, otherwise we cannot update them here as *
   * we need to use the same opacity in Compton scattering and the other terms in scattering */
  hydro_to_fullrad(pD);

  return;

}


/* update H due to bulk Comptonization, given the old values of H, velocity and dtsigmas*/
/* This function basically inverts a 3x3 matrix */
void UpdateHcomp(Real *Hnew, Real Hold[3], const Real dtsigma, Real Vel[3])
{
  Real Ha1, Ha2, Ha3, Hb1, Hb2, Hc1;
  Real InvHa1, InvHa2, InvHa3, InvHb1, InvHb2, InvHc1;
  Real Delta;
  
  Ha1 = 1.0 + 2.0 * dtsigma * Vel[0] * Vel[0] / Crat;
  Ha2 = 2.0 * dtsigma * Vel[0] * Vel[1] / Crat;
  Ha3 = 2.0 * dtsigma * Vel[0] * Vel[2] /Crat;
  Hb1 = 1.0 + 2.0 * dtsigma * Vel[1] * Vel[1] / Crat;
  Hb2 = 2.0 * dtsigma * Vel[1] * Vel[2] / Crat;
  Hc1 = 1.0 + 2.0 * dtsigma * Vel[2] * Vel[2] / Crat;
  Delta = -Ha3 * Ha3 * Hb1 + 2.0 * Ha2 * Ha3 * Hb2 - Ha1 * Hb2 * Hb2 - Ha2 * Ha2 * Hc1 + Ha1 * Hb1 * Hc1;

  InvHa1 = Hb1 * Hc1 - Hb2 * Hb2;
  InvHa2 = Ha3 * Hb2 - Ha2 * Hc1;
  InvHa3 = -Ha3 * Hb1 + Ha2 * Hb2;
  InvHb1 = Ha1 * Hc1 - Ha3 * Ha3;
  InvHb2 = Ha2 * Ha3 - Ha1 * Hb2;
  InvHc1 = -Ha2 * Ha2 + Ha1 * Hb1;
  
  Hnew[0] = (InvHa1 * Hold[0] + InvHa2 * Hold[1] + InvHa3 * Hold[2]) / Delta;
  Hnew[1] = (InvHa2 * Hold[0] + InvHb1 * Hold[1] + InvHb2 * Hold[2]) / Delta;
  Hnew[2] = (InvHa3 * Hold[0] + InvHb2 * Hold[1] + InvHc1 * Hold[2]) / Delta;

  return;

}

/* update radiation pressure tensor according to  *
 * \partial K/\partial t = Delta S/3\delta t I + 1/6\sigma_s V\codt H I - 1/10 sigma_s(V H + H V + V/cdot H I) 
 * + sigma_s(-7 * v^2/(30Crat) I + 7 * VV/(10 Crat)) J
 * K[0] = K[0][0], K[1] = K[0][1], K[2] = K[1][1], K[3] = K[0][2], K[4] = K[1][2], K[5] = K[2][2]
 */
void UpdateKcomp(Real *Knew, Real Kold[6], Real Hnew[3], Real Vel[3], const Real Jnew, const Real dtsigma, const Real Source)
{
  /* Source is the term dt * Crat * Sigma_s * 4 (T - Tr)/T_e J */
  Real vsquar, vdotH;
  vsquar = Vel[0] * Vel[0] + Vel[1] * Vel[1] + Vel[2] * Vel[2];
  vdotH = Vel[0] * Hnew[0] + Vel[1] * Hnew[1] + Vel[2] * Hnew[2];
  /*
  Knew[0] = Kold[0] + Source/3.0 + dtsigma * (-7.0 * vsquar /(30.0 * Crat) + 7.0 * Vel[0] * Vel[0] /(10.0 * Crat)) * Jnew
            + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[0] * Hnew[0] * 2.0  + vdotH);
  Knew[1] = Kold[1] + dtsigma * (7.0 * Vel[0] * Vel[1] /(10.0 * Crat)) * Jnew
            - 0.1 * dtsigma * (Vel[0] * Hnew[1] + Hnew[0] * vel[1]);
  Knew[2] = Kold[2] + Source/3.0 + dtsigma * (-7.0 * vsquar /(30.0 * Crat) + 7.0 * Vel[1] * Vel[1] /(10.0 * Crat)) * Jnew
            + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[1] * Hnew[1] * 2.0  + vdotH);
  Knew[3] = Kold[3] + dtsigma * (7.0 * Vel[0] * Vel[2] /(10.0 * Crat)) * Jnew
            - 0.1 * dtsigma * (Vel[0] * Hnew[2] + Hnew[0] * vel[2]);
  Knew[4] = Kold[4] + dtsigma * (7.0 * Vel[1] * Vel[2] /(10.0 * Crat)) * Jnew
            - 0.1 * dtsigma * (Vel[1] * Hnew[2] + Hnew[1] * vel[2]);
  Knew[5] = Kold[5] + Source/3.0 + dtsigma * (-7.0 * vsquar /(30.0 * Crat) + 7.0 * Vel[2] * Vel[2] /(10.0 * Crat)) * Jnew
            + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[2] * Hnew[2] * 2.0  + vdotH);
*/
  /* Do not include v^2 term */
  Knew[0] = Kold[0] + Source/3.0 + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[0] * Hnew[0] * 2.0  + vdotH);
  Knew[1] = Kold[1] - 0.1 * dtsigma * (Vel[0] * Hnew[1] + Hnew[0] * Vel[1]);
  Knew[2] = Kold[2] + Source/3.0 + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[1] * Hnew[1] * 2.0  + vdotH);
  Knew[3] = Kold[3] - 0.1 * dtsigma * (Vel[0] * Hnew[2] + Hnew[0] * Vel[2]);
  Knew[4] = Kold[4] - 0.1 * dtsigma * (Vel[1] * Hnew[2] + Hnew[1] * Vel[2]);
  Knew[5] = Kold[5] + Source/3.0 + dtsigma * vdotH/6.0 - 0.1 * dtsigma * (Vel[2] * Hnew[2] * 2.0  + vdotH);
  
  return;
}





/* Update the opacity for the whole grid, including the ghost zones */
void UpdateOpacity(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k,ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;


/* Assumes ghost zone conserved variables have been set by
 * bvals routines.  These values are used to set B, chi, eps,
 * etc. so loops include RadGrid ghost zones*/
  if (pG->Nx[0] > 1) {
    ioff = nghost - Radghost;
    il -= Radghost;
    iu += Radghost;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - Radghost;
    jl -= Radghost;
    ju += Radghost;
  } else joff = 0;
  if (pG->Nx[2] > 1) {
    koff = nghost - Radghost;
    kl -= Radghost;
    ku += Radghost;
  } else koff = 0;

  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
        ig = i + ioff;

/* ------------------------------------*/
/* update the opacity */
        for(ifr=0; ifr<nf; ifr++) {
          get_full_opacity(pG,ifr,ig,jg,kg,&(pRG->R[k][j][i][ifr].Sigma[0]));
        }/* end ifr */
      }/* end i */
    }/* end j */
  }/* end k */


}



void GetVelguess(DomainS *pD)
{

  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);

  int nDim, i, ioff, joff, koff;
  int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;


  nDim = 1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;



/* Including the ghost zones */

  if(nDim > 1){
    jl -= Radghost;
    ju += Radghost;

  }

  if(nDim > 2){
    kl -= Radghost;
    ku += Radghost;

  }

  /* flag = 0 for scattering opacity only */
  EstimateVel(il, iu, jl, ju, kl, ku, 1, pRG, pG);

}

/* Estimate velocity from il to iu, jl to ju, kl to ku */
/* Flag is used to distinguish whether it is for scattering opacity or absorption opacity */

/* The estimated velocity will be the true solution if the system is in steady state */


void EstimateVel(const int il, const int iu, const int jl, const int ju, const int kl, const int ku, const int flag, RadGridS *pRG, GridS *pG)
{

  Real Er, sigma, rho, dt;
  Real Fr[3], Pr[6], M0[3], Vel[3];
  /* the matrix coefficient */
  Real Vx11, Vx12, Vx13, Vy11, Vy12, Vz11, RHS[3], factor;
  
  int i, j, k, ifr, n;
  int nDim, ioff, joff, koff;
  
  nDim = 1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

  ioff = nghost - Radghost;
  joff = 0;
  koff = 0;
  
  if(nDim > 1){
    joff = nghost - Radghost;
  }

  if(nDim > 2){
    koff = nghost - Radghost;
  }

  
/* estimated the velocity at half time step */
  dt = 0.5 * pG->dt;


  if(Vguessflag){

    for(k=kl; k<=ku; k++){
      for(j=jl; j<=ju; j++){
        for(i=il; i<=iu; i++){
/* calculate Er and Fr used in the update */

/* first, calculate frequency weighted Er, Pr and Fr */

          Er = 0.0;
          for(n=0; n<3; n++)
            Fr[n] = 0.0;
          for(n=0; n<6; n++)
            Pr[n] = 0.0;

          sigma = 0.0;

          for(ifr=0; ifr<pRG->nf; ifr++){
            Er += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].J);
          for(n=0; n<3; n++)
            Fr[n] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].H[n]);
          for(n=0; n<6; n++)
            Pr[n] += (pRG->wnu[ifr] * pRG->R[k][j][i][ifr].K[n]);


            /* flag =0 for absorption opacity, flag = 1 for scattering opacity */
            sigma += (pRG->wnu[ifr] * ((1-flag) * pRG->R[k][j][i][ifr].Sigma[1]+ flag * pRG->R[k][j][i][ifr].Sigma[2]));

          }

          Er *= (4.0 * PI);
          for(n=0; n<3; n++)
            Fr[n] *= (4.0 * PI);
          for(n=0; n<6; n++)
            Pr[n] *= (4.0 * PI);

    /* Pr[0] = Prxx, Pr[1] = Prxy, Pr[2] = Pryy, Pr[3] = Prxz, Pr[4] = Pryz, Pr[5] = Przz */
/* Now we have Er, Fr and Pr for this cell */

          M0[0] = Prat * Fr[0] / Crat + pG->U[k+koff][j+joff][i+ioff].M1;
          M0[1] = Prat * Fr[1] / Crat + pG->U[k+koff][j+joff][i+ioff].M2;
          M0[2] = Prat * Fr[2] / Crat + pG->U[k+koff][j+joff][i+ioff].M3;

          rho = pG->U[k+koff][j+joff][i+ioff].d;


         
          
          Vx11 = rho * (1.0 + dt * sigma * Crat) + Prat * dt * sigma * (Er + Pr[0])/Crat;
          Vx12 = dt * Prat * sigma * Pr[1] / Crat;
          Vx13 = dt * Prat * sigma * Pr[3] / Crat;
          Vy11 = rho * (1.0 + dt * sigma * Crat) + Prat * dt * sigma * (Er + Pr[2])/Crat;
          Vy12 = dt * Prat * sigma * Pr[4] / Crat;
          Vz11 = rho * (1.0 + dt * sigma * Crat) + Prat * dt * sigma * (Er + Pr[5])/Crat;
          
          RHS[0] = pG->U[k+koff][j+joff][i+ioff].M1 + dt * M0[0] * Crat * sigma;
          RHS[1] = pG->U[k+koff][j+joff][i+ioff].M2 + dt * M0[1] * Crat * sigma;
          RHS[2] = pG->U[k+koff][j+joff][i+ioff].M3 + dt * M0[2] * Crat * sigma;
          
          factor = Vx11 * Vy11 * Vz11 - Vy11 * Vx13 * Vx13 + 2.0 * Vx12 * Vx13 * Vy12 - Vx11 * Vy12 * Vy12 - Vx12 * Vx12 * Vz11;
          factor = 1.0/factor;
          
          pG->Velguess[k+koff][j+joff][i+ioff][0] = factor * (RHS[2] * (Vx12 * Vy12 - Vx13 * Vy11) + RHS[1] * (Vy12 * Vx13 - Vx12 * Vz11) + RHS[0] * (Vy11 * Vz11 - Vy12 * Vy12));
          
          pG->Velguess[k+koff][j+joff][i+ioff][1] = factor * (RHS[2] * (Vx12 * Vx13 - Vx11 * Vy12) + RHS[1] * (Vx11 * Vz11 - Vx13 * Vx13) + RHS[0] * (Vx13 * Vy12 - Vx12 * Vz11));
          
          pG->Velguess[k+koff][j+joff][i+ioff][2] = factor * (RHS[2] * (Vx11 * Vy11 - Vx12 * Vx12) + RHS[1] * (Vx12 * Vx13 - Vx11 * Vy12) + RHS[0] * (Vx12 * Vy12 - Vx13 * Vy11));
          
                  

 
        }
      }
    }/* end k */
  }else{
    for(k=kl; k<=ku; k++){
      for(j=jl; j<=ju; j++){
        for(i=il; i<=iu; i++){


          rho = pG->U[k+koff][j+joff][i+ioff].d;


          Vel[0] = pG->U[k+koff][j+joff][i+ioff].M1 / rho;
          Vel[1] = pG->U[k+koff][j+joff][i+ioff].M2 / rho;
          Vel[2] = pG->U[k+koff][j+joff][i+ioff].M3 / rho;

          for(n=0; n<3; n++)
            pG->Velguess[k+koff][j+joff][i+ioff][n] = Vel[n];

        }
      }
    }/* end k */


  }/* end else */



}




/* Function to calculate the reduction factor of speed of light due to opacity */

void GetSpeedfactor(DomainS *pD)
{


  RadGridS *pRG=(pD->RadGrid);


  int i, j, k, ifr, nDim, n;
  int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;

  Real  sigmas, sigmaa, alpha;
  Real dS[3];
  dS[0] = pRG->dx1;
  dS[1] = pRG->dx2;
  dS[2] = pRG->dx3;

#ifdef CYLINDRICAL
  const Real *r=(pD->Grid)->r;
  int offset = nghost - Radghost;
#endif


  nDim = 1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;


/* Including the ghost zones */

  if(nDim > 1){
    jl -= Radghost;
    ju += Radghost;

  }

  if(nDim > 2){
    kl -= Radghost;
    ku += Radghost;

  }




  for(k=kl; k<=ku; k++){
    for(j=jl; j<=ju; j++){
      for(i=il; i<=iu; i++){
#ifdef CYLINDRICAL
/* The scale factor r[i] is the same for each i, for different angles j */
        dS[1] = pRG->dx2 * r[i+offset];
#endif
        for(ifr=0; ifr<pRG->nf; ifr++){

          sigmas = pRG->R[k][j][i][ifr].Sigma[2];
/* The absorption opacity in front of I */
          sigmaa = pRG->R[k][j][i][ifr].Sigma[1];

/* for three direction */
          for(n=0; n<nDim; n++){
            ReduceVelocity(sigmaa+sigmas, dS[n], &alpha);

            pRG->Speedfactor[k][j][i][ifr][n] = alpha;

          }/* end Dim */

        }/* end ifr */

      }/* end i */
    }/* end j */
  }/* end k */


  return;
}


/* sol takes the initial guess, calculated without the velocity dependent terms */
/* inisol takes the solution at time step n */
/* Md and RHS are pre-allocated memory */
/* They are used as temporary memory */
/* This is for the 3D case */


void Absorption_old(const int nf, const int N, const int Absflag, Real **sol, Real **inisol, Real ***Ma, Real **Mdcoef, Real **Tcoef, Real *Md, Real *RHS, Real **Divi, int *flag)
{


  int n;
  const int MAXIte = 15;
  int count;
  const Real TOL = 1.e-15;
  Real Tgas, Tgas3, Tgas4;
  Real residual;
  int line=N-1;
  int ifr;




/* Now solve the non-linear set of equations for each cell */
  for(ifr=0; ifr<nf; ifr++){


    Tgas = sol[ifr][N];
    Tgas3 = SQR(Tgas) * Tgas;
    Tgas4 = Tgas * Tgas3;

    *flag = 1;


/* First, calculate the RHS for the guess solution */
    for(n=0; n<N-1; n++){
/* set Md, Md is used to invert the matrix */
      Md[n] = 4.0 * Mdcoef[ifr][n] * Tgas3;

      RHS[n] = Mdcoef[ifr][n] * Tgas4 - (inisol[ifr][n] - inisol[ifr][N-1]) + Absflag * (Divi[ifr][n] - Divi[ifr][N-1]);
      RHS[n] += (Ma[ifr][n][0] * sol[ifr][n] - Ma[ifr][N-1][0] * sol[ifr][N-1]);

    }/* end n */
/* Now the line N-1 */
/* Now line=N-1 */
    Md[line] = 4.0 * Mdcoef[ifr][line] * Tgas3;

    RHS[line] = Mdcoef[ifr][line] * Tgas4 + (Ma[ifr][line][0] + Ma[ifr][line][1]) * sol[ifr][line] - inisol[ifr][line] + Absflag * Divi[ifr][line];
    for(n=0; n<N-1; n++){
      RHS[line] += (Ma[ifr][n][1] * sol[ifr][n]);

    }

/* Now the last line nelements */
    Md[N] = Tcoef[ifr][0] + 4.0 * Tcoef[ifr][1] * Tgas3;
    RHS[N] = Tcoef[ifr][0] * (Tgas - inisol[ifr][N]) + Tcoef[ifr][1] * Tgas4;
    for(n=0; n<N; n++){
      RHS[N] += (Ma[ifr][n][2] * sol[ifr][n]);
    }

/* now calculate the norm of residual */
    residual = 0.0;
    for(n=0; n<=N; n++)
      residual += SQR(RHS[n]);

    residual /= (N+1);

    residual = sqrt(residual);
    count = 0;

/****************************************************/
/* Do the iteration */
    while((residual > TOL) && (count < MAXIte)){
      count++;
/* calculate Inverse(Jacobi) * RHS */
      AbsorptionMatrix(N, Ma[ifr], Md, RHS);

/* The result is stored in RHS */
/* Update the guess solution */
      for(n=0; n<=N; n++){
        sol[ifr][n] -= RHS[n];
      }

/*------------------------------------------------*/
      Tgas = sol[ifr][N];
      Tgas3 = SQR(Tgas) * Tgas;
      Tgas4 = Tgas * Tgas3;


/* update RHS and Md */
      for(n=0; n<N-1; n++){
/* set Md, Md is used to invert the matrix */
        Md[n] = 4.0 * Mdcoef[ifr][n] * Tgas3;

        RHS[n] = Mdcoef[ifr][n] * Tgas4 - (inisol[ifr][n] - inisol[ifr][N-1]) + Absflag * (Divi[ifr][n] - Divi[ifr][N-1]);
        RHS[n] += (Ma[ifr][n][0] * sol[ifr][n] - Ma[ifr][N-1][0] * sol[ifr][N-1]);

      }/* end n */
/* Now the line N-1 */
/* Now line=N-1 */
      Md[line] = 4.0 * Mdcoef[ifr][line] * Tgas3;
      RHS[line] = Mdcoef[ifr][line] * Tgas4 + (Ma[ifr][line][0] + Ma[ifr][line][1]) * sol[ifr][line] - inisol[ifr][line] + Absflag * Divi[ifr][line];
      for(n=0; n<N-1; n++){
        RHS[line] += (Ma[ifr][n][1] * sol[ifr][n]);

      }

/* Now the last line nelements */
      Md[N] = Tcoef[ifr][0] + 4.0 * Tcoef[ifr][1] * Tgas3;
      RHS[N] = Tcoef[ifr][0] * (Tgas - inisol[ifr][N]) + Tcoef[ifr][1] * Tgas4;
      for(n=0; n<N; n++){
        RHS[N] += (Ma[ifr][n][2] * sol[ifr][n]);
      }

/*------------------------------------------------*/
/* update the Residual */
/* now calculate the norm of residual */
      residual = 0.0;
      for(n=0; n<=N; n++)
        residual += SQR(RHS[n]);

      residual /= (N+1);

      residual = sqrt(residual);

    }

/* When matrix does not converge, do not add radiation source term */
    if(residual != residual){
      for(n=0; n<=N; n++){
        sol[ifr][n] = inisol[ifr][n];
      }
      *flag = 0;
    }

/* Now the solution is stored in sol[j][i] */
    if(residual > TOL)
      printf("Final residual: %e Iterations: %d\n",residual,count);



  }/* end ifr */


}



/* Function to solve the absorption opacity related terms, including all velocity dependent terms together */
/* Reduce the non-linear equation for a single variable */
/* The original solution is stored in imu, and the new solution will be returned in imu */

void Absorption(const int nelements, const Real rho, const Real dtsigma, const Real Told, const Real Jold, Real Vel[3], Real vdotn[], Real vvdotnn[], Real wmu[], Real imu[],Real **ABCm, Real *Tnew)
{
  Real Jfactor, vsquar, vcsquar, vNc, vnmc, J1, J2, J3;
  Real Coef1, Coef2, Coef3;
  Real T4coef, Tcoef, Tconst, Tmax;
  Real *Am, *Bm, *Cm;
  int n;
  
  Am = ABCm[0];
  Bm = ABCm[1];
  Cm = ABCm[2];
  
  vsquar = (SQR(Vel[0]) + SQR(Vel[1]) + SQR(Vel[2]));
  vcsquar = vsquar/(Crat * Crat);

  Jfactor = Prat * (2.0 + dtsigma * Crat * (1.0 + vcsquar));
  J1 = -rho * R_ideal * 0.25/(PI*(Gamma-1) * Jfactor);
  J2 = 0.25 * dtsigma * Prat * Crat * (1.0 + vcsquar)/(PI * Jfactor);
  J3 = 2.0 * Prat * Jold/Jfactor + 0.25 * rho * R_ideal * Told /(PI * (Gamma-1) * Jfactor);
  
  vNc = 1.0 - vdotn[nelements-1]/Crat;
  
  for(n=0; n<nelements-1; n++){
    vnmc = 1.0 - vdotn[n]/Crat;
  
    Am[n] = (1.0 + dtsigma * Crat * vNc)/(1.0 + dtsigma * Crat * vnmc);
    Bm[n] = 3.0 * dtsigma * (vdotn[n] - vdotn[nelements-1]) * 0.25/(PI * (1.0 + dtsigma * Crat * vnmc));
    Cm[n] = (imu[n] - imu[nelements-1])/(1.0 + dtsigma * Crat * vnmc);
  
  }
  
  Am[nelements-1] = 1.0;
  Bm[nelements-1] = 0.0;
  Cm[nelements-1] = 0.0;
  
  Coef1 = 1.0 + dtsigma * Crat * vNc;
  Coef2 = dtsigma * Crat * 0.25 * (1.0 + 3.0 * vdotn[nelements-1]/Crat)/PI;
  Coef3 = imu[nelements-1];
  
  for(n=0; n<nelements; n++){
    Coef1 += (dtsigma * (vsquar + vvdotnn[n]) * Am[n] * wmu[n]  / Crat);
    Coef2 += (-dtsigma * (vsquar + vvdotnn[n]) * Bm[n] * wmu[n] / Crat);
    Coef3 += (-dtsigma * (vsquar + vvdotnn[n]) * Cm[n] * wmu[n] / Crat);
  
  }
  
  Tcoef = Coef1 * J1;
  T4coef = Coef1 * J2;
  Tconst = Coef1 * J3;
  for(n=0; n<nelements; n++){
    T4coef += (-Coef1 * Bm[n] - Coef2 * Am[n]) * wmu[n];
    Tconst += (-Coef1 * Cm[n] - Coef3 * Am[n]) * wmu[n];
  
  }
  
  Tmax = (0.5 * rho * vsquar + Prat * 4.0 * PI * Jold) * (Gamma - 1.0)/(rho * R_ideal) + Told;
  
  /* Now we solve the equation T4coef * T^4 + Tcoef * T + Tconst= 0; */
  if(fabs(T4coef) < 1.e-20){
    *Tnew = -Tconst/Tcoef;
  }else{
    *Tnew =  rtsafe(Tequilibrium, 0.0, Tmax, 1.e-12, T4coef, Tcoef, Tconst, 0.0);
  }
  
  /* With Tnew, we can get other solutions directly */
  
  imu[nelements-1] = (Coef2 * SQR(SQR((*Tnew))) + Coef3)/Coef1;
  for(n=0; n<nelements-1; n++){
    imu[n] = Am[n] * imu[nelements-1] + Bm[n] * SQR(SQR((*Tnew))) + Cm[n];
  
  }
  
  return;

}



/*********************************************************************************************************************
 * The function below is actually not used anymore **************
 **********************************************************************************************************************
 */


/* function to calculate the energy exchange between J and B implicitly to handle the short thermalize time scale */
/* This function assumes that opacity does not change during this time step */
/* But the thermal function must be known in order to calculate implicitly */
/* The opacity is the frequency weighted opacity */
void GetTnew(const Real dt,const Real d,const Real Tgas, const Real J0, const Real Sigma[4], Real *heatcool, Real *Tnew)
{
/* This implicit function assumes that thermal emission is Tgas^4/(4Pi),
 *  isotropic in every direction */

/* The implicit source terms solv the energy excahnge due to absorption */
/* As well as attentuiation of the source terms */

/*******************************************/
/* For absorption, we solve :
 * 4*pI * dJ/dt = CSigma_a(Tgas^4 - 4* Pi * J)
 * de/dt = -PratC Sigma_a(Tgas^4 - 4*PI * J)
 */
  Real SigmaB, SigmaI, Sigmas;
  Real Tr;
  Real coef1, coef2, coef3, coef4, Ersum, pressure, Er0;

  SigmaB  = Sigma[0];
  SigmaI  = Sigma[1];
  Sigmas  = Sigma[2];

/* SigmasJ and SigmasI must be the same in order to conserve energy locally */

  Er0 = J0 * 4.0 * PI;

/* Negative radiation energy density */
  if(Er0 < 0.0 || Prat < TINY_NUMBER || ((SigmaB < TINY_NUMBER) && (SigmaI < TINY_NUMBER))){
    *heatcool = 0.0;
    *Tnew = Tgas;

    return;
  }
  else{

/*----------------------------------------------------------*/
/* First, energy source due to absorption opacity */

/*              The pow function can be very slow sometimes when Er0 is close to 1
                Tr = pow(Er0, 0.25);

*/
    Tr = sqrt(Er0);
    Tr = sqrt(Tr);

    pressure = Tgas * d * R_ideal;
    Ersum = pressure / (Gamma - 1.0) + Prat * Er0;

/* Here assume input gas temperature and Er0 is positive */
    if(Tgas < 0.0)
      ath_error("[FullRad_GetSource]: Negative gas temperature: %e!n\n",Tgas);


    coef1 = dt * Prat * Crat * SigmaB;
    coef2 = d * R_ideal * (1.0 + dt * SigmaI * Crat) / (Gamma - 1.0);
    coef3 = -pressure / (Gamma - 1.0) - dt * SigmaI * Crat * Ersum;
    coef4 = 0.0;

    if(coef1 < 1.e-18){
      (*Tnew) = -coef3 / coef2;
    }
    else{


      if(Tgas > Tr){
        (*Tnew) = rtsafe(Tequilibrium, Tr * (1.0 - 0.01), Tgas * (1.0 + 0.01), 1.e-12, coef1, coef2, coef3,coef4);
      }
      else{
        (*Tnew) = rtsafe(Tequilibrium, Tgas * (1.0 - 0.01), Tr * (1.0 + 0.01), 1.e-12, coef1, coef2, coef3, coef4);
      }
    }



    *heatcool = d * ((*Tnew) - Tgas) * R_ideal/(Gamma - 1.0);



  }

}




void GetSource(const Real dt, const Real Tnew, const Real SigmaB, const Real SigmaI, const Real Jsource, const Real I0, Real *heatcool)
{
  Real Inew, Tnew4;

  Tnew4 = Tnew * Tnew * Tnew * Tnew;

/* Negative radiation energy density */
  if(Prat < TINY_NUMBER){
    *heatcool = 0.0;

    return;
  }
  else{


    Inew = (I0 + Jsource + dt * Crat * SigmaB * Tnew4 / (4.0 * PI))/(1.0 + dt * Crat * SigmaI);

    (*heatcool) = Inew - I0;
/*----------------------------------------------------------*/
/* First, the energy source due to scattering opacity */
/* We solve the equation dI/dt = csigma_s(J - I) */
/* J is held to be a constant during this step */
/* The formal solution is I(t) = (I0 - J) exp(-dt*csigmas) + J */
/* Because for each cell, average of I0 over different angles is J */
/* So total energy is conserved for the scattering process */

/*      (*Scat) = (1.0 - exp(-Crat * Sigmas * dt)) * (J - I0);
 */

  }

}


#endif /* FULL_RADIATION_TRANSFER */
