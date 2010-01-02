#include "../copyright.h"
/*=============================================================================
 * FILE: hlld_sr.c
 *
 * PURPOSE: Compute 1D fluxes using the relativistic Riemann solver described
 * by Mignone, Ugliano, and Bodo.  For the equivalent hydro-only code, refer
 * to hllc_sr.c
 *
 * REFERENCES:
 * A. Mignone, M. Ugliano and G. Bodo, "A five-wave HLL Riemann solver for
 * relativistic MHD", Mon. Not. R. Astron. Soc. 000, 1-15 (2007)
 *
 * V. Honkkila and P. Janhunen, "HLLC solver for ideal relativistic MHD",
 * Journal of Computational Physics, 233, 643 92007
 *============================================================================*/

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

/* computes left/right fluxes from left/right states */
void flux_LR(Cons1D U, Prim1D W, Cons1D *flux, Real Bx, Real* p);

/* computes total pressure */
Real ptot(Prim1D W, Real Bx);

/* want f = 0 */
Real Fstar(Riemann_State *PaL, Riemann_State *PaR, Real* Sc, Real p, Real Bx);

/* performs some computations and returns > 0 for success, 0 for failure */
int get_Riemann_State(Riemann_State *Pv, Real p, Real Bx);

/* computes left/right a states */
void get_astate(Riemann_State *Pa, Real p, Real Bx);

/* computes left/right c states */
void get_cstate(Riemann_State *PaL, Riemann_State *PaR, Cons1D* Uc, Real p, Real Bx);

/* computes min/max signal speeds */
void getMaxSignalSpeeds(const Prim1D Wl, const Prim1D Wr,const Real Bx, Real* low, Real* high);
void MaxChSpeed (const Prim1D W, const Real Bx, Real* cmin, Real* cmax);

/* solves quartic equation defined by a and returns roots in root
 * returns the number of real roots 
 * error specifies an accuracy
 * currently force four real solutions b/c it's physical */
int QUARTIC (real b, real c, real d, real e, real z[]);

/* solves cubic equation defined by a and stores roots in root
 * returns number of real roots */
int CUBIC (real b, real c, real d, real z[])

#ifdef MHD

/* functions for printing conserved/primitive vectors */
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
            const Prim1D Wl, const Prim1D Wr, const Real Bx, Cons1D *pFlux)
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
      if(switch_to_hll){
         
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
            pFlux->E  = Fl.E  + Sl*(PaL.U.E  - Ul.E )
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
            pFlux->E  = Fr.E  + Sr*(PaR.U.E  - Ur.E )
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
                        Real* low, Real* high)
{

  real sl_min,sl_max;
  real sr_min,sr_max;

  MaxChSpeed (Wl, Bx, &sl_min, &sl_max);
  MaxChSpeed (Wr, Bx, &sr_min, &sr_max);

  if (sl_min <= sr_min) {
     *low = sl_min;
  } else {
     *low = sr_min;
  }

  if (sl_max >= sr_max) {
     *high = sl_max;
  } else {
     *high = sr_max;
  }

}

void MaxChSpeed (const Prim1D W, const Real Bx, Real* cmin, Real* cmax)
{
  int   i, iflag;
  real  rhoh
  real  vB, vB2, w_1;
  real  eps2, one_m_eps2, a2_w;
  real  vx, vx2, u0, u02;
  real  a4, a3, a2, a1, a0;
  real  scrh1, scrh2, scrh;
  real  b2;
  real *vp, lambda[4];
  real  MAX_MACH_NUMBER
  static real *cs2, *h;

  rhoh = W.d + (Gamma/Gamma_1)* (W.P);
  cs2 = (Gamma * W.P) / (rhoh);
  
  scrh   = fabs(W.Vx)/sqrt(cs2);
  MAX_MACH_NUMBER = dmax(scrh, MAX_MACH_NUMBER);

  vB  = W.Vx*Bx + W.Vy*W.By + W.VZz*W.BZ);
  u02 = SQR(W.Vx) + SQR(W.Vy) + SQR(W.Vz);
  b2  = SQR(Bx) + SQR(W.By) + SQR(W.Bz);

  if (u02 >= 1.0){
     ath_error("[MaxChSpeed]:! |v|= %10.4 > 1\n",u02);
  }

  if (u02 < 1.e-14) {
    
    /* -----------------------------------------------------
        if total velocity is = 0 eigenvalue equation 
        reduces to a biquadratic:
         
            x^4 + a1*x^2 + a0 = 0
            
         with a0 = cs^2*bx*bx/W, a1 = - a0 - eps^2
       ----------------------------------------------------- */
       
      w_1  = 1.0/(rhoh + b2);   
      eps2 = cs2 + b2*w_1*(1.0 - cs2);
      a0   = cs2*SQR(Bx)*w_1;
      a1   = - a0 - eps2;
      scrh = a1*a1 - 4.0*a0;

      scrh = max(scrh, 0.0);
      scrh = sqrt(0.5*(-a1 + sqrt(scrh)));
      cmax =  scrh;
      cmin = -scrh;
      continue;      
    }

    vB2 = vB*vB;
    u02 = 1.0/(1.0 - u02);
    b2  = b2/u02 + vB2;
    u0  = sqrt(u02);
    w_1 = 1.0/(rhoh + b2);   
    vx  = W.Vx;
    vx2 = SQR(vx);

    if (fabs(Bx) < 1.e-14){

      eps2  = cs2 + b2*w_1*(1.0 - cs2);

      scrh1 = (1.0 - eps2)*u02;
      scrh2 = cs2*vB2*w_1 - eps2;
      
      a2  = scrh1 - scrh2;
      a1  = -2.0*vx*scrh1;
      a0  = vx2*scrh1 + scrh2;

      cmax = 0.5*(-a1 + sqrt(a1*a1 - 4.0*a2*a0))/a2;
      cmin = 0.5*(-a1 - sqrt(a1*a1 - 4.0*a2*a0))/a2;
      continue;

    }else{

      scrh1 = Bx/u02 + vB*vx;  /* -- this is bx/u0 -- */
      scrh2 = scrh1*scrh1;  
                     
      a2_w       = cs2*w_1;
      eps2       = (cs2*rhoh + b2)*w_1;
      one_m_eps2 = u02*rhoh*(1.0 - cs2)*w_1;

    /* ---------------------------------------
         Define coefficients for the quartic  
       --------------------------------------- */
    
      scrh = 2.0*(a2_w*vB*scrh1 - eps2*vx);
      a4 = one_m_eps2 - a2_w*vB2 + eps2;
      a3 = - 4.0*vx*one_m_eps2 + scrh;
      a2 =   6.0*vx2*one_m_eps2 + a2_w*(vB2 - scrh2) + eps2*(vx2 - 1.0);
      a1 = - 4.0*vx*vx2*one_m_eps2 - scrh;
      a0 = vx2*vx2*one_m_eps2 + a2_w*scrh2 - eps2*vx2;
   
      if (a4 < 1.e-12){
        ath_error ("[MaxChSpeed]: Can not divide by a4\n");
      }

      scrh = 1.0/a4;
     
      a3 *= scrh;
      a2 *= scrh;
      a1 *= scrh;
      a0 *= scrh;
      iflag = QUARTIC(a3, a2, a1, a0, lambda);
  
      if (iflag){
        printf ("Can not find max speed:\n");
        SHOW(uprim,i);
        print ("QUARTIC: f(x) = %12.6e + x*(%12.6e + x*(%12.6e ",
                 a0*a4, a1*a4, a2*a4);
        print ("+ x*(%12.6e + x*%12.6e)))\n", a3*a4, a4);
        ath_error ("[MaxChSpeed]: Can't find max speed\n");
      }

      cmax = max(lambda[3], lambda[2]);
      cmax = max(cmax, lambda[1]);
      cmax = max(cmax, lambda[0]);
      cmax = min(cmax, 1.0);

      cmin = min(lambda[3], lambda[2]);
      cmin = min(cmin, lambda[1]);
      cmin = min(cmin, lambda[0]);
      cmin = max(cmin, -1.0);

    }
  }
}

/* ******************************************** */
int QUARTIC (real b, real c, real d, 
             real e, real z[])
/* 
 *
 * PURPOSE:
 *
 *   Solve a quartic equation in the form 
 *
 *      z^4 + bz^3 + cz^2 + dz + e = 0
 *
 *   For its purpose, it is assumed that ALL 
 *   roots are real. This makes things faster.
 *
 *
 * ARGUMENTS
 *
 *   b, c,
 *   d, e  (IN)  = coefficient of the quartic
 *                 z^4 + bz^3 + cz^2 + dz + e = 0
 *
 *   z[]   (OUT) = a vector containing the 
 *                 (real) roots of the quartic
 *   
 *
 * REFERENCE:
 *
 *   http://www.1728.com/quartic2.htm 
 * 
 *
 *
 ********************************************** */
{
  int    n, ifail;
  real b2, f, g, h;
  real a2, a1, a0, u[4];
  real p, q, r, s;
  static real three_256 = 3.0/256.0;
  static real one_64 = 1.0/64.0;
  
  b2 = b*b;

  f = c - b2*0.375;
  g = d + b2*b*0.125 - b*c*0.5;
  h = e - b2*b2*three_256 + 0.0625*b2*c - 0.25*b*d;
  
  a2 = 0.5*f;
  a1 = (f*f - 4.0*h)*0.0625;
  a0 = -g*g*one_64;

  ifail = CUBIC(a2, a1, a0, u);

  if (ifail)return(1);

  if (u[1] < 1.e-14){

    p = sqrt(u[2]);
    s = 0.25*b;
    z[0] = z[2] = - p - s;
    z[1] = z[3] = + p - s;
    
  }else{
  
    p = sqrt(u[1]);
    q = sqrt(u[2]);
  
    r = -0.125*g/(p*q);
    s =  0.25*b;
     
    z[0] = - p - q + r - s;
    z[1] =   p - q - r - s;
    z[2] = - p + q - r - s;
    z[3] =   p + q + r - s;

  }  

  /* ----------------------------------------------
       verify that cmax and cmin satisfy original 
       equation
     ---------------------------------------------- */  

  for (n = 0; n < 4; n++){
    s = e + z[n]*(d + z[n]*(c + z[n]*(b + z[n])));
    if (s != s) {
      printf ("Nan found in QUARTIC \n");
      return(1);
    }
    if (fabs(s) > 1.e-6) {
      printf ("Solution does not satisfy f(z) = 0; f(z) = %12.6e\n",s);
      return(1);
    }
  }

  return(0);
/*  
  printf (" z: %f ; %f ; %f ; %f\n",z[0], z[1], z[2], z[3]);
  */
}
/* *************************************************** */
int CUBIC(real b, real c, real d, real z[])
/* 
 *
 * PURPOSE:
 *
 *   Solve a cubic equation in the form 
 *
 *      z^3 + bz^2 + cz + d = 0
 *
 *   For its purpose, it is assumed that ALL 
 *   roots are real. This makes things faster.
 *
 *
 * ARGUMENTS
 *
 *   b, c, d (IN)  = coefficient of the cubic
 *                    z^3 + bz^2 + cz + d = 0
 *
 *   z[]   (OUT)   = a vector containing the 
 *                   (real) roots of the cubic.
 *                   Roots should be sorted
 *                   in increasing order.
 *   
 *
 * REFERENCE:
 *
 *   http://www.1728.com/cubic2.htm 
 *
 *
 *
 ***************************************************** */
{
  real b2, g2;
  real f, g, h;
  real i, i2, j, k, m, n, p;
  static real one_3 = 1.0/3.0, one_27=1.0/27.0;

  b2 = b*b;
 
/*  ----------------------------------------------
     the expression for f should be 
     f = c - b*b/3.0; however, to avoid negative
     round-off making h > 0.0 or g^2/4 - h < 0.0
     we let c --> c(1- 1.1e-16)
    ---------------------------------------------- */

  f  = c*(1.0 - 1.e-16) - b2*one_3;
  g  = b*(2.0*b2 - 9.0*c)*one_27 + d; 
  g2 = g*g;
  i2 = -f*f*f*one_27;
  h  = g2*0.25 - i2;

/* --------------------------------------------
     Real roots are possible only when 
   
               h <= 0 
   -------------------------------------------- */

  if (h > 1.e-12){
    printf ("Only one real root (%12.6e)!\n", h);
  }
  if (i2 < 0.0){
/*
    printf ("i2 < 0.0 %12.6e\n",i2);
    return(1);
*/
    i2 = 0.0;
  }

/* --------------------------------------
       i^2 must be >= g2*0.25
   -------------------------------------- */
  
  i = sqrt(i2);       /*  > 0   */
  j = pow(i, one_3);  /*  > 0   */
  k = -0.5*g/i;

/*  this is to prevent unpleseant situation 
    where both g and i are close to zero       */

  k = (k < -1.0 ? -1.0:k);
  k = (k >  1.0 ?  1.0:k);
  
  k = acos(k)*one_3;       /*  pi/3 < k < 0 */
 
  m = cos(k);              /*   > 0   */
  n = sqrt(3.0)*sin(k);    /*   > 0   */
  p = -b*one_3;

  z[0] = -j*(m + n) + p;
  z[1] = -j*(m - n) + p;
  z[2] =  2.0*j*m + p;

/* ------------------------------------------------------
    Since j, m, n > 0, it should follow that from
    
      z0 = -jm - jn + p
      z1 = -jm + jn + p
      z2 = 2jm + p
      
    z2 is the greatest of the roots, while z0 is the 
    smallest one.
   ------------------------------------------------------ */
      
  return(0);
}

#undef MAX_ITER

#endif
#endif
