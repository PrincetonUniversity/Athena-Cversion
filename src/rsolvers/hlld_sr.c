#include "../copyright.h"

/*=============================================================================
 * FILE: flux_hlld_rmhd.c
 *
 * PURPOSE: Compute 1D fluxes using the relativistic Riemann solver described
 * by Mignone, Ugliano, and Bodo.  For the equivalent hydro-only code, refer
 * to flux_hllc_rhd.c
 *
 * REFERENCES:
 *
 * A. Mignone, M. Ugliano and G. Bodo, "A five-wave HLL Riemann solver for
 * relativistic MHD", Mon. Not. R. Astron. Soc. 000, 1-15 (2007)
 *
 * V. Honkkila and P. Janhunen, "HLLC solver for ideal relativistic MHD",
 * Journal of Computational Physics, 233, 643 92007
 *
 *=============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef HLLD_FLUX
#ifdef SPECIAL_RELATIVITY

#if (NSCALARS > 0)
#error : The SR HLLD flux does not work with passive scalars.
#endif

#define MAX_ITER 20

typedef struct RIEMANN_STATE{
      int fail;
      Real vx, vy, vz;
      Real Bx, By, Bz;
      Real Kx, Ky, Kz, K2;
      Real w, sw, p, rho;
      Cons1D U, R;
      Real S, Sa;
} Riemann_State;

void printCons1D(const Cons1D *U);
void printPrim1D(const Prim1D *W);

void flux_LR(Cons1D U, Prim1D W, Cons1D *flux, Real Bx, Real* p);
Real ptot(Prim1D W, Real Bx);
Real Fstar(Riemann_State *PaL, Riemann_State *PaR, Real* Sc, Real p, Real Bx);
int get_Riemann_State(Riemann_State *Pv, Real p, Real Bx);
void get_astate(Riemann_State *Pa, Real p, Real Bx);
void get_cstate(Riemann_State *PaL, Riemann_State *PaR, Cons1D* Uc, Real p, Real Bx);
void getMaxSignalSpeeds(const Prim1D Wl, const Prim1D Wr,
                        const Real Bx, const Real error,
                        Real* low, Real* high);
int solveQuartic(double* a, double* root, double error);
int solveCubic(double* a, double* root);

#ifdef MHD
void printCons1D(const Cons1D *U){
   printf("d:  %.6e\n",U->d);
   printf("E:  %.6e\n",U->E);
   printf("Mx: %.6e\n",U->Mx);
   printf("My: %.6e\n",U->My);
   printf("Mz: %.6e\n",U->Mz);
   printf("By: %.6e\n",U->By);
   printf("Bz: %.6e\n",U->Bz);
   printf("\n");
}

void printPrim1D(const Prim1D *W){
   printf("d:  %.6e\n",W->d);
   printf("P:  %.6e\n",W->P);
   printf("Vx: %.6e\n",W->Vx);
   printf("Vy: %.6e\n",W->Vy);
   printf("Vz: %.6e\n",W->Vz);
   printf("By: %.6e\n",W->By);
   printf("Bz: %.6e\n",W->Bz);
   printf("\n");
}
#endif


/*----------------------------------------------------------------------------*/
/* flux_hlld_rmhd
 *   Input Arguments:
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *     Bx = B in direction of slice at cell interface
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1D Ul, const Cons1D Ur,
                    const Prim1D Wl, const Prim1D Wr, const Real Bx, 
                    Cons1D *pFlux)
{
   int k;
   int switch_to_hll;
   Real scrh;
   Cons1D Fl, Fr;
   Real Pl, Pr;
   Cons1D Uhll, Fhll;
   Prim1D Whll;
   Real Sl, Sr;
   
   Real p0, f0, p, f, dp, dS_1;
   Cons1D Uc;
   Riemann_State PaL, PaR;

   Real a,b,c;

   Real Sc;


   /* find min/max wave speeds */
   getMaxSignalSpeeds(Wl,Wr,Bx,1e-12,&Sl,&Sr);
   
   /* compute L/R fluxes */
   flux_LR(Ul,Wl,&Fl,Bx,&Pl);
   flux_LR(Ur,Wr,&Fr,Bx,&Pr);

   if(Sl >= 0.0){
      /*printf("Flux_L\n");*/
      pFlux->d  = Fl.d;
      pFlux->Mx = Fl.Mx;
      pFlux->My = Fl.My;
      pFlux->Mz = Fl.Mz;
      pFlux->By = Fl.By;
      pFlux->Bz = Fl.Bz;
      pFlux->E  = Fl.E;
   }
   else if(Sr <= 0.0){
      /*printf("Flux_R\n");*/
      pFlux->d  = Fr.d;
      pFlux->Mx = Fr.Mx;
      pFlux->My = Fr.My;
      pFlux->Mz = Fr.Mz;
      pFlux->By = Fr.By;
      pFlux->Bz = Fr.Bz;
      pFlux->E  = Fr.E;
   }
   else{
      
      /* Compute HLL average state */

      dS_1 = 1.0/(Sr - Sl);

      Uhll.d  = (Sr*Ur.d  - Sl*Ul.d  + Fl.d  - Fr.d ) * dS_1;
      Uhll.Mx = (Sr*Ur.Mx - Sl*Ul.Mx + Fl.Mx - Fr.Mx) * dS_1;
      Uhll.My = (Sr*Ur.My - Sl*Ul.My + Fl.My - Fr.My) * dS_1;
      Uhll.Mz = (Sr*Ur.Mz - Sl*Ul.Mz + Fl.Mz - Fr.Mz) * dS_1;
      Uhll.By = (Sr*Ur.By - Sl*Ul.By + Fl.By - Fr.By) * dS_1;
      Uhll.Bz = (Sr*Ur.Bz - Sl*Ul.Bz + Fl.Bz - Fr.Bz) * dS_1;
      Uhll.E  = (Sr*Ur.E  - Sl*Ul.E  + Fl.E  - Fr.E ) * dS_1;

      Fhll.d  = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
      Fhll.Mx = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
      Fhll.My = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
      Fhll.Mz = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
      Fhll.By = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
      Fhll.Bz = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;
      Fhll.E  = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;

      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      pFlux->E = Fhll.E;

      return;

      /* set up some variables */

      PaL.S = Sl;
      PaR.S = Sr;

      PaL.R.d  = Sl*Ul.d  - Fl.d;
      PaL.R.Mx = Sl*Ul.Mx - Fl.Mx;
      PaL.R.My = Sl*Ul.My - Fl.My;
      PaL.R.Mz = Sl*Ul.Mz - Fl.Mz;
      PaL.R.By = Sl*Ul.By - Fl.By;
      PaL.R.Bz = Sl*Ul.Bz - Fl.Bz;
      PaL.R.E  = Sl*Ul.E  - Fl.E;

      PaR.R.d  = Sr*Ur.d  - Fr.d;
      PaR.R.Mx = Sr*Ur.Mx - Fr.Mx;
      PaR.R.My = Sr*Ur.My - Fr.My;
      PaR.R.Mz = Sr*Ur.Mz - Fr.Mz;
      PaR.R.By = Sr*Ur.By - Fr.By;
      PaR.R.Bz = Sr*Ur.Bz - Fr.Bz;
      PaR.R.E  = Sr*Ur.E  - Fr.E;

      scrh = MAX(Pl,Pr);
      if(SQR(Bx)/scrh < 0.01){ /* -- try the B->0 limit */
         
         a = Sr - Sl;
         b = PaR.R.E - PaL.R.E + Sr*PaL.R.Mx - Sl*PaR.R.Mx;
         c = PaL.R.Mx*PaR.R.E - PaR.R.Mx*PaL.R.E;
         scrh = b*b - 4.0*a*c;
         scrh = MAX(scrh, 0.0);
         p0 = 0.5*(-b + sqrt(scrh))*dS_1;

      }
      else{
         Whll.d  = 0.5*(Wl.d  + Wr.d );
         Whll.Vx = 0.5*(Wl.Vx + Wr.Vx);
         Whll.Vy = 0.5*(Wl.Vy + Wr.Vy);
         Whll.Vz = 0.5*(Wl.Vz + Wr.Vz);
         Whll.By = 0.5*(Wl.By + Wr.By);
         Whll.Bz = 0.5*(Wl.Bz + Wr.Bz);
         Whll.P  = 0.5*(Wl.P  + Wr.P );
         Cons1D_to_Prim1D((const Cons1D*)(&Uhll),(Prim1D*)(&Whll),
                          (const Real*)(&Bx));
         p0 = ptot(Whll,Bx);
      }

      /* -- check if guess makes sense */

      switch_to_hll = 0;
      f0 = Fstar(&PaL, &PaR, &Sc, p0, Bx);
      if(f0 != f0 || PaL.fail) switch_to_hll = 1;

      /* -- Root finder -- */

      k = 0;
      if (fabs(f0) > 1.e-12 && !switch_to_hll){
         p = 1.025*p0; f = f0;
         for(k = 1; k < MAX_ITER; k++){
            
            f = Fstar(&PaL, &PaR, &Sc, p, Bx);
            if ( f != f || PaL.fail || (k > 9) ||
                (fabs(f) > fabs(f0) && k > 4) ) {
               switch_to_hll = 1;
               break;
            }

            dp = (p - p0)/(f - f0)*f;

            p0 = p; f0 = f;
            p -= dp;
            if (p < 0.0) p = 1.e-6;
            if (fabs(dp) < 1.e-6*p || fabs(f) < 1.e-6) break;
         }
      }else p = p0;

      /* too many iter? --> use HLL */

      if(PaL.fail) switch_to_hll = 1;
      if(switch_to_hll || 1){
         
         printf("Flux_HLL\n");

         pFlux->d  = Fhll.d;
         pFlux->Mx = Fhll.Mx;
         pFlux->My = Fhll.My;
         pFlux->Mz = Fhll.Mz;
         pFlux->By = Fhll.By;
         pFlux->Bz = Fhll.Bz;
         pFlux->E  = Fhll.E;

         return;
      }

      /* solution should be reliable */

      if(PaL.Sa >= -1.e-6){
         get_astate(&PaL, p, Bx);

         printf("Flux_aL\n");
         
         pFlux->d  = Fl.d  + Sl*(PaL.U.d  - Ul.d );
         pFlux->Mx = Fl.Mx + Sl*(PaL.U.Mx - Ul.Mx);
         pFlux->My = Fl.My + Sl*(PaL.U.My - Ul.My);
         pFlux->Mz = Fl.Mz + Sl*(PaL.U.Mz - Ul.Mz);
         pFlux->By = Fl.By + Sl*(PaL.U.By - Ul.By);
         pFlux->Bz = Fl.Bz + Sl*(PaL.U.Bz - Ul.Bz);
         pFlux->E  = Fl.E  + Sl*(PaL.U.E  - Ul.E );
      }
      else if(PaR.Sa <= 1.e-6){
         get_astate(&PaR, p, Bx);

         printf("Flux_aR\n");

         pFlux->d  = Fr.d  + Sr*(PaR.U.d  - Ur.d );
         pFlux->Mx = Fr.Mx + Sr*(PaR.U.Mx - Ur.Mx);
         pFlux->My = Fr.My + Sr*(PaR.U.My - Ur.My);
         pFlux->Mz = Fr.Mz + Sr*(PaR.U.Mz - Ur.Mz);
         pFlux->By = Fr.By + Sr*(PaR.U.By - Ur.By);
         pFlux->Bz = Fr.Bz + Sr*(PaR.U.Bz - Ur.Bz);
         pFlux->E  = Fr.E  + Sr*(PaR.U.E  - Ur.E );
      }
      else{
         get_cstate(&PaL,&PaR,&Uc,p,Bx);
         if(Sc > 0.0){
            printf("Flux_cL\n");
            
            pFlux->d  = Fl.d  + Sl*(PaL.U.d  - Ul.d )
                              + PaL.Sa*(Uc.d  - PaL.U.d );
            pFlux->Mx = Fl.Mx + Sl*(PaL.U.Mx - Ul.Mx)
                              + PaL.Sa*(Uc.Mx - PaL.U.Mx);
            pFlux->My = Fl.My + Sl*(PaL.U.My - Ul.My)
                              + PaL.Sa*(Uc.My - PaL.U.My);
            pFlux->Mz = Fl.Mz + Sl*(PaL.U.Mz - Ul.Mz)
                              + PaL.Sa*(Uc.Mz - PaL.U.Mz);
            pFlux->By = Fl.By + Sl*(PaL.U.By - Ul.By)
                              + PaL.Sa*(Uc.By - PaL.U.By);
            pFlux->Bz = Fl.Bz + Sl*(PaL.U.Bz - Ul.Bz)
                              + PaL.Sa*(Uc.Bz - PaL.U.Bz);
            pFlux->By = Fl.E  + Sl*(PaL.U.E  - Ul.E )
                              + PaL.Sa*(Uc.E  - PaL.U.E );
         }
         else{
            printf("Flux_cR\n");

            pFlux->d  = Fr.d  + Sr*(PaR.U.d  - Ur.d )
                              + PaR.Sa*(Uc.d  - PaR.U.d );
            pFlux->Mx = Fr.Mx + Sr*(PaR.U.Mx - Ur.Mx)
                              + PaR.Sa*(Uc.Mx - PaR.U.Mx);
            pFlux->My = Fr.My + Sr*(PaR.U.My - Ur.My)
                              + PaR.Sa*(Uc.My - PaR.U.My);
            pFlux->Mz = Fr.Mz + Sr*(PaR.U.Mz - Ur.Mz)
                              + PaR.Sa*(Uc.Mz - PaR.U.Mz);
            pFlux->By = Fr.By + Sr*(PaR.U.By - Ur.By)
                              + PaR.Sa*(Uc.By - PaR.U.By);
            pFlux->Bz = Fr.Bz + Sr*(PaR.U.Bz - Ur.Bz)
                              + PaR.Sa*(Uc.Bz - PaR.U.Bz);
            pFlux->By = Fr.E  + Sr*(PaR.U.E  - Ur.E )
                              + PaR.Sa*(Uc.E  - PaR.U.E );
         }
      }

   } /* end if for sL,sR */

}

void flux_LR(Cons1D U, Prim1D W, Cons1D *flux, Real Bx, Real* p){
   Real vB, b2, wtg2, Bmag2, pt;
   Real g, g2, g_2, h, gmmr, theta;
   Real bx, by, bz;

   /* calculate enthalpy */

   theta = W.P/W.d;
   gmmr = Gamma / Gamma_1;

   h = 1.0 + gmmr*theta;

   /* calculate gamma */

   g   = U.d/W.d;
   g2  = SQR(g);
   g_2 = 1.0/g2;

   vB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
   Bmag2 = SQR(Bx) + SQR(W.By) + SQR(W.Bz);

   bx = g*(  Bx*g_2 + vB*W.Vx);
   by = g*(W.By*g_2 + vB*W.Vy);
   bz = g*(W.Bz*g_2 + vB*W.Vz);

   b2 = Bmag2*g_2 + vB*vB;

   pt = W.P + 0.5*b2;

   wtg2 = (W.d*h + b2)*g2;

   flux->d  = U.d*W.Vx;
   flux->Mx = wtg2*W.Vx*W.Vx - bx*bx + pt;
   flux->My = wtg2*W.Vy*W.Vx - by*bx;
   flux->Mz = wtg2*W.Vz*W.Vx - bz*bx;
   flux->By = W.Vx*W.By - Bx*W.Vy;
   flux->Bz = W.Vx*W.Bz - Bx*W.Vz;
   flux->E  = U.Mx;

   *p = pt;
}

Real Fstar(Riemann_State *PaL, Riemann_State *PaR, Real* Sc, Real p, Real Bx){
   int success = 1;
   Real dK, Bxc, Byc, Bzc;
   Real sBx, fun;
   Real vxcL, KLBc;
   Real vxcR, KRBc;

   sBx = Bx > 0.0 ? 1.0 : -1.0;

   success *= get_Riemann_State(PaL, p, Bx);
   success *= get_Riemann_State(PaR, p, Bx);

   /* comptue B from average state */
   
   dK = PaR->Kx - PaL->Kx + 1.e-12;

   Bxc = Bx*dK;
   Byc =   PaR->By*(PaR->Kx - PaR->vx)
         - PaL->By*(PaL->Kx - PaL->vx)
         + Bx*(PaR->vy - PaL->vy);
   Bzc =   PaR->Bz*(PaR->Kx - PaR->vx)
         - PaL->Bz*(PaL->Kx - PaL->vx)
         + Bx*(PaR->vz - PaL->vz);

   KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
   KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;

   vxcL = PaL->Kx - dK*Bx*(1.0 - PaL->K2)/(PaL->sw*dK - KLBc);
   vxcR = PaR->Kx - dK*Bx*(1.0 - PaR->K2)/(PaR->sw*dK - KRBc);

   PaL->Sa = PaL->Kx;
   PaR->Sa = PaR->Kx;
   *Sc      = 0.5*(vxcL + vxcR);
   fun     = vxcL - vxcR;

   /* check if state makes physical sense */

   success *= (vxcL - PaL->Kx) > -1.e-6;
   success *= (PaR->Kx - vxcR) > -1.e-6;

   success *= (PaL->S - PaL->vx) < 0.0;
   success *= (PaR->S - PaR->vx) > 0.0;
   
   success *= (PaR->w - p) > 0.0;
   success *= (PaL->w - p) > 0.0;
   success *= (PaL->Sa - PaL->S) > -1.e-6;
   success *= (PaR->S  - PaR->Sa) > -1.e-6;

   PaL->fail = !success;

   return fun;
}

int get_Riemann_State(Riemann_State *Pv, Real p, Real Bx){
   Real A, C, G, X, s;
   Real vx, vy, vz, scrh, S;
   Cons1D R;

   S = Pv->S;
   R = Pv->R;

   A = R.Mx + p*(1.0 - S*S) - S*R.E;
   C = R.By*R.My + R.Bz*R.Mz;
   G = SQR(R.By) + SQR(R.Bz);
   X = Bx*(A*S*Bx + C) - (A + G)*(S*p + R.E);

   vx = ( Bx*(A*Bx + C*S) - (R.Mx + p)*(G + A) );
   vy = ( - (A + G - Bx*Bx*(1.0 - S*S))*R.My + R.By*(C + Bx*(S*R.Mx - R.E)) );
   vz = ( - (A + G - Bx*Bx*(1.0 - S*S))*R.Mz + R.Bz*(C + Bx*(S*R.Mx - R.E)) );

   scrh = vx*R.Mx + vy*R.My + vz*R.Mz;
   scrh = X*R.E - scrh;
   Pv->w = p + scrh/(X*S - vx);

   Pv->vx = vx/X;
   Pv->vy = vy/X;
   Pv->vz = vz/X;

   Pv->Bx = Bx;
   Pv->By = -(R.By*(S*p + R.E) - Bx*R.My)/A;
   Pv->Bz = -(R.Bz*(S*p + R.E) - Bx*R.Mz)/A;

   s = Bx > 0.0 ? 1.0 : -1.0;
   if(S < 0.0) s *= -1.0;
   
   if(Pv->w < 0.0){
      return (0); /* -- failure -- */
   }

   Pv->sw = s*sqrt(Pv->w);

   scrh = 1.0/(S*p + R.E + Bx*Pv->sw);
   Pv->Kx = scrh*(R.Mx + p +   Bx*Pv->sw);
   Pv->Ky = scrh*(R.My     + R.By*Pv->sw);
   Pv->Kz = scrh*(R.Mz     + R.Bz*Pv->sw);

   Pv->K2 = SQR(Pv->Kx) + SQR(Pv->Ky) + SQR(Pv->Kz);
   return (1); /* -- successful -- */
}

void get_astate(Riemann_State *Pa, Real p, Real Bx){
   Cons1D *Ua,*R;
   Real vB,S;

   Ua = &(Pa->U);
   S  = Pa->S;
   R  = &(Pa->R);

   Ua->d  = R->d/(S - Pa->vx);
   Ua->By = (R->By - Bx*Pa->vy)/(S - Pa->vx);
   Ua->Bz = (R->Bz - Bx*Pa->vz)/(S - Pa->vx);

   vB    = Pa->vx*Bx + Pa->vy*Ua->By + Pa->vz*Ua->Bz;
   Ua->E = (R->E + p*Pa->vx - vB*Bx)/(S - Pa->vx);

   Ua->Mx = (Ua->E + p)*Pa->vx - vB*Bx;
   Ua->My = (Ua->E + p)*Pa->vy - vB*Ua->By;
   Ua->Mz = (Ua->E + p)*Pa->vz - vB*Ua->Bz;
   
}

void get_cstate(Riemann_State *PaL, Riemann_State *PaR, Cons1D* Uc, Real p, Real Bx){
   Cons1D* Ua;
   Real dK;
   Real vxcL, vycL, vzcL, KLBc;
   Real vxcR, vycR, vzcR, KRBc;
   Real vxc, vyc, vzc, vBc;
   Real Bxc, Byc, Bzc, Sa, vxa;

   get_astate(PaL,p,Bx);
   get_astate(PaR,p,Bx);
   dK = (PaR->Kx - PaL->Kx) + 1.e-12;

   Bxc = Bx*dK;
   Byc =   PaR->By*(PaR->Kx - PaR->vx)
         - PaL->By*(PaL->Kx - PaL->vx)
         + Bx*(PaR->vy - PaL->vy);
   Bzc =   PaR->Bz*(PaR->Kx - PaR->vx)
         - PaL->Bz*(PaL->Kx - PaL->vx)
         + Bx*(PaR->vz - PaL->vz);

   Bxc = Bx;
   Byc /= dK;
   Bzc /= dK;

   Uc->By = Byc;
   Uc->Bz = Bzc;

   KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
   KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;

   vxcL = PaL->Kx -     Bx*(1.0 - PaL->K2)/(PaL->sw - KLBc);
   vxcR = PaR->Kx -     Bx*(1.0 - PaR->K2)/(PaR->sw - KRBc);

   vycL = PaL->Ky - Uc->By*(1.0 - PaL->K2)/(PaL->sw - KLBc);
   vycR = PaR->Ky - Uc->By*(1.0 - PaR->K2)/(PaR->sw - KRBc);

   vzcL = PaL->Kz - Uc->Bz*(1.0 - PaL->K2)/(PaL->sw - KLBc);
   vzcR = PaR->Kz - Uc->Bz*(1.0 - PaR->K2)/(PaR->sw - KRBc);

   vxc = 0.5*(vxcL + vxcR);
   vyc = 0.5*(vycL + vycR);
   vzc = 0.5*(vzcL + vzcR);

   if (vxc > 0.0){
      get_astate(PaL, p, Bx);
      Ua = &(PaL->U);
      Sa = PaL->Sa;
      vxa = PaL->vx;
   }
   else{
      get_astate(PaR, p, Bx);
      Ua = &(PaL->U);
      Sa = PaR->Sa;
      vxa = PaR->vx;
   }

   vBc = vxc*Bx + vyc*Uc->By + vzc*Uc->Bz;

   Uc->d = Ua->d*(Sa - vxa)/(Sa - vxc);
   Uc->E = (Sa*Ua->E - Ua->Mx + p*vxc - vBc*Bx)/(Sa - vxc);

   Uc->Mx = (Uc->E + p)*vxc - vBc*Bx;
   Uc->My = (Uc->E + p)*vyc - vBc*Uc->By;
   Uc->Mz = (Uc->E + p)*vzc - vBc*Uc->Bz;
   
}

Real ptot(Prim1D W, Real Bx){
   Real vel2, Bmag2, vB;
   double pt;

   vel2  = SQR(W.Vx) + SQR(W.Vy) + SQR(W.Vz);
   Bmag2 = SQR(Bx)   + SQR(W.By) + SQR(W.Bz);
   vB    = W.Vx*Bx   + W.Vy*W.By + W.Vz*W.Bz;

   pt = W.P + 0.5*(Bmag2*(1.0 - vel2) + vB*vB);
   return (pt);
}

void getMaxSignalSpeeds(const Prim1D Wl, const Prim1D Wr,
                        const Real Bx, const Real error,
                        Real* low, Real* high){

   Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
   Real lpl,lpr;        /* largest roots, Mignone Eq 55 */

   Real vlsq,vrsq;
   Real gammal, gammar;
   Real gammal2, gammar2;
   Real gammal4, gammar4;
   Real rhohl,rhohr;
   Real cslsq,csrsq;
   Real cslsq_1,csrsq_1;
   Real Blsq,Brsq;
   Real vDotBl,vDotBr;
   Real b0l,b0r;
   Real bxl,bxr;
   Real blsq,brsq;
   Real Ql,Qr;
   Real al[5],ar[5];
   Real rl[4],rr[4];
   Real discl,discr;

   int nl,nr;
   int i;

   /*printf("SPEED: Wl\n");
   printPrim1D(&Wl);
   printf("SPEED: Wr\n");
   printPrim1D(&Wr);*/

   /*printf("Bx: %f\n",Bx);*/


   rhohl = Wl.d + (Gamma/Gamma_1) * (Wl.P);
   rhohr = Wr.d + (Gamma/Gamma_1) * (Wr.P);

   /* Mignone RHD(2005) Eq 4 */
   cslsq = (Gamma * Wl.P) / (rhohl); 
   csrsq = (Gamma * Wr.P) / (rhohr);

   cslsq_1 = 1.0 - cslsq;
   csrsq_1 = 1.0 - csrsq;

   vlsq = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
   vrsq = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);

   gammal = 1.0 / sqrt(1 - vlsq);
   gammar = 1.0 / sqrt(1 - vrsq);

   gammal2 = SQR(gammal);
   gammar2 = SQR(gammar);

   gammal4 = SQR(gammal2);
   gammar4 = SQR(gammar2);

   Blsq = SQR(Bx) + SQR(Wl.By) + SQR(Wl.Bz);
   Brsq = SQR(Bx) + SQR(Wr.By) + SQR(Wr.Bz);

   vDotBl = Wl.Vx*Bx + Wl.Vy*Wl.By + Wl.Vz*Wl.Bz;
   vDotBr = Wr.Vx*Bx + Wr.Vy*Wr.By + Wr.Vz*Wr.Bz;

   b0l = gammal * vDotBl;
   b0r = gammar * vDotBr;

   bxl = Bx/gammal2 + Wl.Vx*vDotBl;
   bxr = Bx/gammar2 + Wr.Vx*vDotBr; 

   blsq = Blsq / gammal2 + SQR(vDotBl);
   brsq = Brsq / gammar2 + SQR(vDotBr);

   if( fabs(Bx) < error ){

      /*printf("Quadratic\n\n");*/

      /* Mignone Eq 58 */

      Ql = blsq - cslsq*vDotBl;
      Qr = brsq - csrsq*vDotBr;

      /*printf("Ql: %e\n",Ql);
        printf("Qr: %e\n\n",Qr);*/

      al[2] = rhohl*(cslsq + gammal2*cslsq_1) + Ql;
      ar[2] = rhohr*(csrsq + gammar2*csrsq_1) + Qr;

      al[1] = -2.0 * rhohl * gammal2 * Wl.Vx * cslsq_1;
      ar[1] = -2.0 * rhohr * gammar2 * Wr.Vx * csrsq_1;

      al[0] = rhohl*(gammal2*Wl.Vx*Wl.Vx*cslsq_1 - cslsq) - Ql;
      ar[0] = rhohr*(gammar2*Wr.Vx*Wr.Vx*csrsq_1 - csrsq) - Qr;

      /*printf("al[2]: %e\n",al[2]);
      printf("al[1]: %e\n",al[1]);
      printf("al[0]: %e\n\n",al[0]);

      printf("ar[2]: %e\n",ar[2]);
      printf("ar[1]: %e\n",ar[1]);
      printf("ar[0]: %e\n\n",ar[0]);*/

      discl = sqrt(al[1]*al[1] - 4.0*al[2]*al[0]);
      discr = sqrt(ar[1]*ar[1] - 4.0*ar[2]*ar[0]);

      lml = (-al[1] - discl) / (2.0*al[2]);
      lpl = (-al[1] + discl) / (2.0*al[2]);

      lmr = (-ar[1] - discr) / (2.0*ar[2]);
      lpr = (-ar[1] + discr) / (2.0*ar[2]);

      /*printf("lml: %f\n",lml);
      printf("lpl: %f\n",lpl);
      printf("lmr: %f\n",lmr);
      printf("lpr: %f\n\n",lpr);*/
   }
   else{
      /*printf("Quartic\n\n");*/

      al[4] = -SQR(b0l)*cslsq + cslsq_1*rhohl*gammal4 + 
         gammal2*(blsq + cslsq*rhohl);
      ar[4] = -SQR(b0r)*csrsq + csrsq_1*rhohr*gammar4 + 
         gammar2*(brsq + csrsq*rhohr);

      al[3] = 2*b0l*bxl*cslsq - 4*cslsq_1*rhohl*gammal4*Wl.Vx - 
         2*gammal2*(blsq + cslsq*rhohl)*Wl.Vx;
      ar[3] = 2*b0r*bxr*csrsq - 4*csrsq_1*rhohr*gammar4*Wr.Vx - 
         2*gammar2*(brsq + csrsq*rhohr)*Wr.Vx;

      al[2] = SQR(b0l)*cslsq - SQR(bxl)*cslsq - gammal2*(blsq + cslsq*rhohl) + 
         6*cslsq_1*rhohl*gammal4*SQR(Wl.Vx)
         + gammal2*SQR(Wl.Vx)*(blsq + cslsq*rhohl);
      ar[2] = SQR(b0r)*csrsq - SQR(bxr)*csrsq - gammar2*(brsq + csrsq*rhohr) + 
         6*csrsq_1*rhohr*gammar4*SQR(Wr.Vx)
         + gammar2*SQR(Wr.Vx)*(brsq + csrsq*rhohr);

      al[1] = -2*b0l*bxl*cslsq + 2*gammal2*Wl.Vx*(blsq + cslsq*rhohl) - 
         4*cslsq_1*rhohl*gammal4*SQR(Wl.Vx)*Wl.Vx;
      ar[1] = -2*b0r*bxr*csrsq + 2*gammar2*Wr.Vx*(brsq + csrsq*rhohr) - 
         4*csrsq_1*rhohr*gammar4*SQR(Wr.Vx)*Wr.Vx;

      al[0] = SQR(bxl)*cslsq - gammal2*SQR(Wl.Vx)*(blsq + cslsq*rhohl) + 
         cslsq_1*rhohl*gammal4*SQR(Wl.Vx)*SQR(Wl.Vx);
      ar[0] = SQR(bxr)*csrsq - gammar2*SQR(Wr.Vx)*(brsq + csrsq*rhohr) + 
         csrsq_1*rhohr*gammar4*SQR(Wr.Vx)*SQR(Wr.Vx);

/*      for(i=0; i<5; i++)
         printf("al[%d]: %f\n",i,al[i]);
      for(i=0; i<5; i++)
      printf("ar[%d]: %f\n",i,ar[i]);*/

      nl = solveQuartic(al,rl,1.0e-10);
      nr = solveQuartic(ar,rr,1.0e-10);

      /*printf("nl: %d\n",nl);
        printf("nr: %d\n",nr);*/

      if(nl == 0){
         lml = -1.0;
         lpl = 1.0;
      }
      else{
         lml = rl[0];
         lpl = rl[0];

         /* find smallest and largest roots */
         for(i=1; i<nl; i++){
            if(rl[i] < lml)
               lml = rl[i];
            if(rl[i] > lpl)
               lpl = rl[i];
         }
      }

      if(nr == 0){
         lmr = -1.0;
         lpr = 1.0;
      }
      else{
         lmr = rr[0];
         lpr = rr[0];

         /* find smallest and largest roots */
         for(i=1; i<nr; i++){
            if(rr[i] < lmr)
               lmr = rr[i];
            if(rr[i] > lpr)
               lpr = rr[i];
         }
      }
      
   }

   /* Mi */
   
   if(lml < lmr)
      *low = lml;
   else
      *low = lmr;

   if(lpl > lpr)
      *high = lpl;
   else
      *high = lpr;
}

/* algorithm from Abramowitz and Stegun (1972) */
int solveQuartic(double* a, double* root, double error){
   double a3,a2,a1,a0;
   double cubicA[4];
   double cubicR[3];
   double R,D,E;
   double rSq,dSq,eSq;
   double y1;
   double alpha,tmp;
   int i,nRoots;

   /* normalize coefficients */

   a3 = a[3]/a[4];
   a2 = a[2]/a[4];
   a1 = a[1]/a[4];
   a0 = a[0]/a[4];

   /* Reduce to cubic equation */

   cubicA[3] = 1.0;
   cubicA[2] = -a2;
   cubicA[1] = a1*a3 - 4.0*a0;
   cubicA[0] = 4.0*a2*a0 - SQR(a1) - SQR(a3)*a0;

   nRoots = solveCubic(cubicA,cubicR);

   /* select one of the cubic roots to serve as y1 */

   for(i=0; i<nRoots; i++){
      y1 = cubicR[i];
      rSq = SQR(a3)/4.0 - a2 + y1;
      tmp = SQR(y1) - 4.0*a0;

      if( (rSq > -error) && (tmp > -error) )
         break;
   }
   

   /* Reduce to solving two quadratic equations */

   nRoots = 0;

   /* if(R^2 == 0) */
   
   if(rSq < error){
      R = 0;
      if(tmp >= 0.0)
         tmp = sqrt(tmp);
      else
         tmp = 0.0;

      alpha = 3.0*SQR(a3)/4.0 - 2.0*a2;
         
      dSq = alpha + 2.0*tmp;
      eSq = alpha - 2.0*tmp;
   }
   
   /* if(R^2 > 0) */
   else{
      R = sqrt(rSq);
      tmp = (4.0*a3*a2 - 8.0*a1 - SQR(a3)*a3)/4.0;
      alpha = 3.0*SQR(a3)/4.0 - rSq - 2.0*a2;
      
      dSq = alpha + tmp/R;
      eSq = alpha - tmp/R;
   }
   
   if(dSq >= 0){
      D = sqrt(dSq);
      root[0] = -a3/4.0 + R/2.0 + D/2.0;
      root[1] = -a3/4.0 + R/2.0 - D/2.0;
      nRoots += 2;
   }

   /* incorrect but approximately right for situation */
   else{
      D = 0;
      root[0] = -a3/4.0 + R/2.0 + D/2.0;
      root[1] = -a3/4.0 + R/2.0 - D/2.0;
      nRoots += 2;
   }
   
   if(eSq >= 0){
      E = sqrt(eSq);
      root[nRoots] = -a3/4.0 - R/2.0 + E/2.0;
      root[nRoots+1] = -a3/4.0 - R/2.0 - E/2.0;
      nRoots += 2;
   }

   /* incorrect but approximately right for situation */
   else{
      E = 0;
      root[nRoots] = -a3/4.0 - R/2.0 + E/2.0;
      root[nRoots+1] = -a3/4.0 - R/2.0 - E/2.0;
      nRoots += 2;
   }
      
   return nRoots;
}   


/* Cubic equation solver */
int solveCubic(double* a, double* root){
   double a2,a1,a0;
   double q,r,d,e;
   double qCubed;
   double theta,sqrtQ;

   a2 = a[2]/a[3];
   a1 = a[1]/a[3];
   a0 = a[0]/a[3];

   q = (SQR(a2) - 3.0*a1)/9.0;
   r = (2*SQR(a2)*a2 - 9.0*a2*a1 + 27.0*a0)/54.0;
   qCubed = SQR(q)*q;
   d = qCubed - SQR(r);
   
   /* Three real roots */

   if(d >= 0){
      if(qCubed > 0)
         theta = acos(r/sqrt(qCubed));
      else
         theta = 0;

      sqrtQ = sqrt(q);

      root[0] = -2.0 * sqrtQ * cos(theta/3.0) - a2/3.0;
      root[1] = -2.0 * sqrtQ * cos((theta + 2.0*PI)/3.0) - a2/3.0;
      root[2] = -2.0 * sqrtQ * cos((theta + 4.0*PI)/3.0) - a2/3.0;

      return 3;
   }

   /* One real root */

   else{
      e = pow(sqrt(-d) + fabs(r),(ONE_3RD) );

      if(r > 0) e = -e;

      root[0] = (e + q/e) - a2/3.0;

      return 1;
   }
}

#undef MAX_ITER

#endif
#endif
