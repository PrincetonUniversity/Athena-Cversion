#include "../copyright.h"
/*=============================================================================
 * FILE: hlle_sr.c
 *
 * PURPOSE: Compute 1D fluxes using an HLLE-type relativistic Riemann solver.
 *   Works for both hydro and MHD, but is very diffusive.
 *
 * HISTORY:
 *   First version written by Kevin Tian, 2007
 *   Updated by Jonathan Fulton, February 2009
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef HLLE_FLUX
#ifdef SPECIAL_RELATIVITY

#if (NSCALARS > 0)
#error : The SR HLLE flux does not work with passive scalars.
#endif

void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p);
void getMaxSignalSpeeds(const Prim1DS Wl, const Prim1DS Wr,
                        const Real Bx, const Real error, Real* low, Real* high);
int solveQuartic(double* a, double* root, double error);
int solveCubic(double* a, double* root);

/*----------------------------------------------------------------------------*/
/* hlle_sr
 *   Input Arguments:
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *     Bx = B in direction of slice at cell interface
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pFlux)
{
   Cons1DS Fl, Fr;
   Cons1DS Uhll, Fhll;
   Real Pl, Pr;
   Real Sl, Sr;
   Real dS_1;

   /* find min/max wave speeds */
   getMaxSignalSpeeds(Wl,Wr,Bx,1.0e-6,&Sl,&Sr);
   /* compute L/R fluxes */
   flux_LR(Ul,Wl,&Fl,Bx,&Pl);
   flux_LR(Ur,Wr,&Fr,Bx,&Pr);

   if(Sl >= 0.0){
      /*printf("Flux_L\n");*/
      pFlux->d  = Fl.d;
      pFlux->Mx = Fl.Mx;
      pFlux->My = Fl.My;
      pFlux->Mz = Fl.Mz;
      pFlux->E  = Fl.E;
#ifdef MHD
      pFlux->By = Fl.By;
      pFlux->Bz = Fl.Bz;
#endif

      return;
   }
   else if(Sr <= 0.0){
      /*printf("Flux_R\n");*/
      pFlux->d  = Fr.d;
      pFlux->Mx = Fr.Mx;
      pFlux->My = Fr.My;
      pFlux->Mz = Fr.Mz;
      pFlux->E  = Fr.E;
#ifdef MHD
      pFlux->By = Fr.By;
      pFlux->Bz = Fr.Bz;
#endif

      return;
   }
   else{
      
      /* Compute HLL average state */

      dS_1 = 1.0/(Sr - Sl);

      Uhll.d  = (Sr*Ur.d  - Sl*Ul.d  + Fl.d  - Fr.d ) * dS_1;
      Uhll.Mx = (Sr*Ur.Mx - Sl*Ul.Mx + Fl.Mx - Fr.Mx) * dS_1;
      Uhll.My = (Sr*Ur.My - Sl*Ul.My + Fl.My - Fr.My) * dS_1;
      Uhll.Mz = (Sr*Ur.Mz - Sl*Ul.Mz + Fl.Mz - Fr.Mz) * dS_1;
      Uhll.E  = (Sr*Ur.E  - Sl*Ul.E  + Fl.E  - Fr.E ) * dS_1;
#ifdef MHD
      Uhll.By = (Sr*Ur.By - Sl*Ul.By + Fl.By - Fr.By) * dS_1;
      Uhll.Bz = (Sr*Ur.Bz - Sl*Ul.Bz + Fl.Bz - Fr.Bz) * dS_1;
#endif

      Fhll.d  = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
      Fhll.Mx = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
      Fhll.My = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
      Fhll.Mz = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
      Fhll.E  = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;
#ifdef MHD
      Fhll.By = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
      Fhll.Bz = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;
#endif

      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
#ifdef MHD
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
#endif

      return;
   }
}


void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p)
{
   Real wtg2, pt, g, g2, g_2, h, gmmr, theta;
#ifdef MHD
   Real bx, by, bz, vB, b2, Bmag2;
#endif

   /* calculate enthalpy */

   theta = W.P/W.d;
   gmmr = Gamma / Gamma_1;

   h = 1.0 + gmmr*theta;

   /* calculate gamma */

   g   = U.d/W.d;
   g2  = SQR(g);
   g_2 = 1.0/g2;

   pt = W.P;
   wtg2 = W.d*h*g2;

#ifdef MHD
   vB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
   Bmag2 = SQR(Bx) + SQR(W.By) + SQR(W.Bz);

   bx = g*(  Bx*g_2 + vB*W.Vx);
   by = g*(W.By*g_2 + vB*W.Vy);
   bz = g*(W.Bz*g_2 + vB*W.Vz);

   b2 = Bmag2*g_2 + vB*vB;

   pt += 0.5*b2;
   wtg2 += b2*g2;
#endif

   flux->d  = U.d*W.Vx;
   flux->Mx = wtg2*W.Vx*W.Vx + pt;
   flux->My = wtg2*W.Vy*W.Vx;
   flux->Mz = wtg2*W.Vz*W.Vx;
   flux->E  = U.Mx;
#ifdef MHD
   flux->Mx -= bx*bx;
   flux->My -= by*bx;
   flux->Mz -= bz*bx;
   flux->By = W.Vx*W.By - Bx*W.Vy;
   flux->Bz = W.Vx*W.Bz - Bx*W.Vz;
#endif

   *p = pt;
}


void getMaxSignalSpeeds(const Prim1DS Wl, const Prim1DS Wr,
                        const Real Bx, const Real error, Real* low, Real* high)
{
   Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
   Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
   Real rhohl, rhohr, cslsq, csrsq, cslsq_1, csrsq_1, vlsq, vrsq;
   Real gammal, gammar, gammal2, gammar2, gammal4, gammar4;
#ifdef MHD
   Real Blsq,Brsq,vDotBl,vDotBr,b0l,b0r,bxl,bxr,blsq,brsq,Ql,Qr;
   Real rl[4],rr[4];
#endif
   Real al[5],ar[5];
   Real discl,discr;
   int nl,nr;
   int i;

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

#ifdef MHD
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
#endif /* MHD */

      /* Mignone Eq 58.  These formulae are used in the case of hydro */

      al[0] = rhohl*(gammal2*Wl.Vx*Wl.Vx*cslsq_1 - cslsq);
      ar[0] = rhohr*(gammar2*Wr.Vx*Wr.Vx*csrsq_1 - csrsq);

      al[1] = -2.0 * rhohl * gammal2 * Wl.Vx * cslsq_1;
      ar[1] = -2.0 * rhohr * gammar2 * Wr.Vx * csrsq_1;

      al[2] = rhohl*(cslsq + gammal2*cslsq_1);
      ar[2] = rhohr*(csrsq + gammar2*csrsq_1);

#ifdef MHD
      Ql = blsq - cslsq*vDotBl;
      Qr = brsq - csrsq*vDotBr;
      al[2] += Ql;
      ar[2] += Qr;
      al[0] -= Ql;
      ar[0] -= Qr;
#endif 

      discl = sqrt(al[1]*al[1] - 4.0*al[2]*al[0]);
      discr = sqrt(ar[1]*ar[1] - 4.0*ar[2]*ar[0]);

      lml = (-al[1] - discl) / (2.0*al[2]);
      lpl = (-al[1] + discl) / (2.0*al[2]);

      lmr = (-ar[1] - discr) / (2.0*ar[2]);
      lpr = (-ar[1] + discr) / (2.0*ar[2]);

#ifdef MHD
   } else{
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

      nl = solveQuartic(al,rl,1.0e-10);
      nr = solveQuartic(ar,rr,1.0e-10);

      if(nl == 0){
         lml = -1.0;
         lpl = 1.0;
      } else{
         lml = rl[0];
         lpl = rl[0];

         /* find smallest and largest roots */
         for(i=1; i<nl; i++){
            if(rl[i] < lml) lml = rl[i];
            if(rl[i] > lpl) lpl = rl[i];
         }
      }

      if(nr == 0){
         lmr = -1.0;
         lpr = 1.0;
      } else{
         lmr = rr[0];
         lpr = rr[0];

         /* find smallest and largest roots */
         for(i=1; i<nr; i++){
            if(rr[i] < lmr) lmr = rr[i];
            if(rr[i] > lpr) lpr = rr[i];
         }
      }
      
   }
#endif /* MHD */

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

#endif /* SPECIAL_RELATIVITY */
#endif /* HLLE_FLUX */
