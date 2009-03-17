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
 *
 * HISTORY: Written by Jonathan Fulton, February 2009
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

#ifdef HLLD_FLUX
#ifdef SPECIAL_RELATIVITY
     
#ifdef HYDRO
#error : The HLLD flux only works for MHD.
#endif /* HYDRO */

#if (NSCALARS > 0)
#error : The SR HLLD flux does not work with passive scalars.
#endif

void printCons1D(const Cons1D *U);
void printPrim1D(const Prim1D *W);

/* stores slowest and fastest signal speeds into low and high */
void getMaxSignalSpeeds(const Prim1D Wl, const Prim1D Wr,
                        const Real Bx, const Real error,
                        Real* low, Real* high);

/* store solution to quartic equation defined by coefficients a[] in r[] */
int solveQuartic(double* a, double* r, const Real error);

/* stores solution to cubic equation defined by a[] in r[] */
int solveCubic(double* a, double* r);

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
/* fluxes
 *   Input Arguments:
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *     Wl,Wr = L/R-states of PRIMITIVE variables at cell interface
 *     Bx = B in direction of slice at cell interface
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

int count=0;

void fluxes(const Cons1D Ul, const Cons1D Ur,
            const Prim1D Wl, const Prim1D Wr, const Real Bx, 
            Cons1D *pFlux)
{
   Cons1D Fl,Fr;           /* left and right fluxes as in Mignone Eq 16 */
   Cons1D Uhll,Fhll;       /* hll-conserved quants, hll-flux */
   Prim1D Whll;            /* hll-primitive values */
   Cons1D Ual,Uar;         /* left-a and right-a conserved quantities */
   Cons1D Fal,Far;         /* left-a and right-a fluxes as in Mignone Eq 16 */
   Cons1D Ucl,Ucr;         /* conserved quantities as in Mignone Eq 16 */
   Cons1D Rl,Rr;           /* constants from Mignone Eq 12: R=lamda*U - F */
   Real spd[5];            /* wave speeds */
   Real pt,ptold;          /* total pressure, used to solve all fluxes */
   Real p0,phll;           /* pressure given Bx->0 and hll-pressure */
   Real ovlrMll;           /* 1.0 / (lamdar - lamdal) */
   Real lllr;              /* lamdal * lamdar */
   Real a,b,c;             /* quadratic equation coefficients */
   Real b0l,bxl,byl,bzl, 
        b0r,bxr,byr,bzr;   /* covariant magnetic field components */
   Real uxl,uyl,uzl,
      uxr,uyr,uzr;         /* four-velocity components of left/right states */
   Real gammal,vsql,vDotBl,
      gammar,vsqr,vDotBr;  /* gamma, v-squared, and dot product of v and B */
   Real wgl,bsql,wl,ptl,
      wgr,bsqr,wr,ptr;     /* gas enthalpy, b-squared, total enthalpy and total pressure */
   Real vxl,vyl,vzl,
      vxr,vyr,vzr;         /* velocity components defined in Mignone Eq 23-25 */
   Real Al,Gl,Cl,Ql,Xl,
      Ar,Gr,Cr,Qr,Xr;      /* constants defined in Mignone Eq 26-30 */
   Real ovspdMinusVxl,
      ovspdMinusVxr;       /* 1.0 / (spd - vx) */
   Real vDotRml,vDotRmr;   /* v dotted with Rm as in Mignone Eq 31 */
   Real wal,war;           /* total enthalpies */
   Real Kxl,Kyl,Kzl,
      Kxr,Kyr,Kzr;         /* K-vector components defined in Mignone Eq 36 */
   Real signBx;            /* the sign of Bx */
   Real etal,etar;         /* eta variable in Mignone Eq 43 */
   Real ovKldenom,
      ovKrdenom;           /* 1.0 / (denom) in Mignone Eq 43 */
   Real Bcy,Bcz;           /* Magnetic field components from Mignone Eq 45 */
   Real aRy,aRz,
      aLy,aLz;             /* Parts of numerator in Mignone Eq 45 */
   Real spdMinVxl,
      spdMinVxr;           /* lamda - Vx */
   Real ovs3Mns1;          /* 1.0 / (spd[3] - spd[1]) */
   Real dKx,Klsq,Krsq,
      KldotB,KrdotB,
      KldotBhc,KrdotBhc;   /* variables from Mignone Eq 49 */
   Real Yl,Yr;             /* Y(p) given in Mignone Eq 49 */
   Real fp,fpold;          /* f(p), Mignone Eq 48 */
   Real mKlsq,mKrsq,
      ovEtalmKdB,
      ovEtarmKdB;          /* variables from Mignone Eq 47 */
   Real vxcl,vycl,vzcl,
      vxcr,vycr,vzcr;      /* velocities given in Mignone Eq 47 */
   Real ovldalMvxcl,
      ovldarMvxcr;         /* denominators in Mignone Eqs 50-51 */
   Real vclDotBc,
      vcrDotBc;            /* V_c dotted with B_c, Mignone Eq 51-52 */
   Real ecrPpt,eclPpt;       /* E_c + pt, Mignone Eq 52 */
   Real jump,slope,error;  /* used for secant root finding */
   Real tmp;               /* temporary variable used to shorten equations and 
                              reduce redundant calculations */
   const Real 
      maxError = 1e-6;     /* Maximum error allowed */
   int status;             /* 0:first loop, 1:bounding,
                              2:iterating full secant method to find root */
   int i;

   /* initialize root-finding variables */
   status = 0;
   jump = 0.05;
   error = 1.0;

   /*printf("Wl\n");
   printf("%.3f %.3f %.3f %.3f\n",Wl.d,Wl.P,Wl.Vx,Wl.By);
   printf("Wr\n");
   printf("%.3f %.3f %.3f %.3f\n",Wr.d,Wr.P,Wr.Vx,Wr.By);
   printf("Ur\n");
   printf("%.3f %.3f %.3f %.3f\n",Ur.d,Ur.E,Ur.Mx,Ur.By);*/
   count++;
   /* printf("COUNT: %d\n",count); */

/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone
 */

   getMaxSignalSpeeds(Wl,Wr,Bx,maxError,spd,spd+4);

   printf("spd[0]: %f\n",spd[0]);
     printf("spd[4]: %f\n",spd[4]);


/*--- Step 2. ------------------------------------------------------------------
 * Compute left and right fluxes
 */

   /* compute covariant magnetic field components */
   
   vsql = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
   vsqr = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);

   gammal = 1.0 / sqrt(1 - vsql);  
   gammar = 1.0 / sqrt(1 - vsqr);

   vDotBl = Wl.Vx*Bx + Wl.Vy*Wl.By + Wl.Vz*Wl.Bz;
   vDotBr = Wr.Vx*Bx + Wr.Vy*Wr.By + Wr.Vz*Wr.Bz;

   b0l = gammal * vDotBl;
   b0r = gammar * vDotBr;

   bxl = (Bx / gammal) + (gammal * vDotBl * Wl.Vx);
   bxr = (Bx / gammar) + (gammar * vDotBr * Wr.Vx);

   byl = (Wl.By / gammal) + (gammal * vDotBl * Wl.Vy);
   byr = (Wr.By / gammar) + (gammar * vDotBr * Wr.Vy);

   bzl = (Wl.Bz / gammal) + (gammal * vDotBl * Wl.Vz);
   bzr = (Wr.Bz / gammar) + (gammar * vDotBr * Wr.Vz);

   /* Compute spacial four-velocity components */

   uxl = gammal * Wl.Vx;
   uxr = gammar * Wr.Vx;

   uyl = gammal * Wl.Vy;
   uyr = gammar * Wr.Vy;
   
   uzl = gammal * Wl.Vz;
   uzr = gammar * Wr.Vz;

   /* Compute total enthalpy and total pressure */

   bsql = SQR(bxl) + SQR(byl) + SQR(bzl) - SQR(b0l);
   bsqr = SQR(bxr) + SQR(byr) + SQR(bzr) - SQR(b0r);

   wgl = Wl.d + (Gamma / Gamma_1) * Wl.P;
   wgr = Wr.d + (Gamma / Gamma_1) * Wr.P;

   wl = wgl + bsql;
   wr = wgr + bsqr;

   ptl = Wl.P + 0.5*bsql;
   ptr = Wr.P + 0.5*bsqr;


   /* Now use Mignone Eq 7 to compute left/right fluxes */

   Fl.d = Ul.d * Wl.Vx;
   Fr.d = Ur.d * Wr.Vx;

   Fl.Mx = wl*uxl*uxl - bxl*bxl + ptl;
   Fr.Mx = wr*uxr*uxr - bxr*bxr + ptr;

   Fl.My = wl*uxl*uyl - bxl*byl;
   Fr.My = wr*uxr*uyr - bxr*byr;

   Fl.Mz = wl*uxl*uzl - bxl*bzl;
   Fr.Mz = wr*uxr*uzr - bxr*bzr;

   Fl.E = Ul.Mx;
   Fr.E = Ur.Mx;

   Fl.By = Ul.By*Wl.Vx - Bx*Wl.Vy;
   Fr.By = Ur.By*Wr.Vx - Bx*Wr.Vy;

   Fl.Bz = Ul.Bz*Wl.Vx - Bx*Wl.Vz;
   Fr.Bz = Ur.Bz*Wr.Vx - Bx*Wr.Vz;

   /*printf("F_L\n");
   printCons1D(&Fl);
   printf("F_R\n");
   printCons1D(&Fr);*/


/*--- Step 3. ------------------------------------------------------------------
 * If supersonic, return appropriate flux.
 */

   if(spd[0] >= 0.0){
      printf("Flux_L\n");
      *pFlux = Fl;
      return;
   }
   if(spd[4] <= 0.0){
      printf("Flux_R\n");
      *pFlux = Fr;
      return;
   }

/*--- Step 4. ------------------------------------------------------------------
 * Initialize total pressure variale according to Mignone Eq 53
 */

   /* compute hll-conserved vector (Eq 29 Mignone 2006)*/
   
   ovlrMll = 1.0 / (spd[4] - spd[0]);

   Uhll.d  = (spd[4]*Ur.d  - spd[0]*Ul.d  + Fl.d  - Fr.d ) * ovlrMll;
   Uhll.Mx = (spd[4]*Ur.Mx - spd[0]*Ul.Mx + Fl.Mx - Fr.Mx) * ovlrMll;
   Uhll.My = (spd[4]*Ur.My - spd[0]*Ul.My + Fl.My - Fr.My) * ovlrMll;
   Uhll.Mz = (spd[4]*Ur.Mz - spd[0]*Ul.Mz + Fl.Mz - Fr.Mz) * ovlrMll;
   Uhll.E  = (spd[4]*Ur.E  - spd[0]*Ul.E  + Fl.E  - Fr.E ) * ovlrMll;
   Uhll.By = (spd[4]*Ur.By - spd[0]*Ul.By + Fl.By - Fr.By) * ovlrMll;
   Uhll.Bz = (spd[4]*Ur.Bz - spd[0]*Ul.Bz + Fl.Bz - Fr.Bz) * ovlrMll;

   /* compute phll using variable converter; need initial guess for Whll */
   Whll.d = (Wl.d+Wr.d)/2.0;
   Whll.Vx = (Wl.Vx+Wr.Vx)/2.0;
   Whll.Vy = (Wl.Vy+Wr.Vy)/2.0;
   Whll.Vz = (Wl.Vz+Wr.Vz)/2.0;
   Whll.P = (Wl.P+Wr.P)/2.0;
   Whll.By = (Wl.By+Wr.By)/2.0;
   Whll.Bz = (Wl.Bz+Wr.Bz)/2.0;

   Cons1D_to_Prim1D((const Cons1D*)(&Uhll),(Prim1D*)(&Whll),
                    (const Real*)(&Bx));   

   /*printf("Uhll\n");
   printCons1D(&Uhll);
   printf("Whll\n");
   printPrim1D(&Whll);*/

   phll = Whll.P;

   /* printf("Phll: %f\n",phll); */

   /* compute hll-conserved flux (Eq 31 Mignone 2006) */
   
   lllr = spd[0]*spd[4];

   Fhll.d  = (spd[4]*Fl.d  - spd[0]*Fr.d  + lllr*(Ur.d  - Ul.d )) * ovlrMll;
   Fhll.Mx = (spd[4]*Fl.Mx - spd[0]*Fr.Mx + lllr*(Ur.Mx - Ul.Mx)) * ovlrMll;
   Fhll.My = (spd[4]*Fl.My - spd[0]*Fr.My + lllr*(Ur.My - Ul.My)) * ovlrMll;
   Fhll.Mz = (spd[4]*Fl.Mz - spd[0]*Fr.Mz + lllr*(Ur.Mz - Ul.Mz)) * ovlrMll;
   Fhll.E  = (spd[4]*Fl.E  - spd[0]*Fr.E  + lllr*(Ur.E  - Ul.E )) * ovlrMll;
   Fhll.By = (spd[4]*Fl.By - spd[0]*Fr.By + lllr*(Ur.By - Ul.By)) * ovlrMll;
   Fhll.Bz = (spd[4]*Fl.Bz - spd[0]*Fr.Bz + lllr*(Ur.Bz - Ul.Bz)) * ovlrMll;

   pFlux->d = Fhll.d;
   pFlux->Mx = Fhll.Mx;
   pFlux->My = Fhll.My;
   pFlux->Mz = Fhll.Mz;
   pFlux->E = Fhll.E;
   pFlux->By = Fhll.By;
   pFlux->Bz = Fhll.Bz;

   return;
   
   /* printf("Fhll\n");
      printCons1D(&Fhll); */

   /* set initial total pressure according to Mignone Eq 53 */

   if( SQR(Bx) / phll < 0.1){
      
      /* Now use Mignone Eq 55 */

      a = 1.0;
      b = Uhll.E - Fhll.Mx;
      c = Uhll.Mx*Fhll.E - Fhll.Mx*Uhll.E;

      /* positive root corresponds to physical solution */

      p0 = ( -b + sqrt(b*b - 4.0*a*c) )/(2.0*a);
      /*printf("P0\n\n");*/
   }
   else{
      /* printf("Phll\n\n"); */
      p0 = phll;
   }

/*--- Step 5. ------------------------------------------------------------------
 * Iterate according to Mignone section 3.4 to find the correct value for 
 * total pressure
 */

   /* set initial pressure */
   pt = p0;
   /* printf("p0: %f\n\n",pt); */

   while(error > maxError){
      /*printf("pt: %e\n",pt);
      printf("error: %e\n",error);
      printf("status: %d\n",status);*/

      /*--- Step 5.1. ----------------------------------------------------------
       * Solve jump conditions across the fast waves
       */

      Rl.d = spd[0]*Ul.d - Fl.d;
      Rr.d = spd[4]*Ur.d - Fr.d;

      Rl.Mx = spd[0]*Ul.Mx - Fl.Mx;
      Rr.Mx = spd[4]*Ur.Mx - Fr.Mx;

      Rl.My = spd[0]*Ul.My - Fl.My;
      Rr.My = spd[4]*Ur.My - Fr.My;

      Rl.Mz = spd[0]*Ul.Mz - Fl.Mz;
      Rr.Mz = spd[4]*Ur.Mz - Fr.Mz;

      Rl.E = spd[0]*Ul.E - Fl.E;
      Rr.E = spd[4]*Ur.E - Fr.E;

      Rl.By = spd[0]*Ul.By - Fl.By;
      Rr.By = spd[4]*Ur.By - Fr.By;

      Rl.Bz = spd[0]*Ul.Bz - Fl.Bz;
      Rr.Bz = spd[4]*Ur.Bz - Fr.Bz;

      /*** First solve for velocities using Mignone Eqs 22-30 ***/

      /* left side */
      Al = Rl.Mx - spd[0]*Rl.E + pt*(1.0 - SQR(spd[0]));
      Gl = Rl.By*Rl.By + Rl.Bz*Rl.Bz;
      Cl = Rl.My*Rl.By + Rl.Mz*Rl.Bz;
      Ql = -Al - Gl + SQR(Bx)*(1.0 - SQR(spd[0]));
      Xl = Bx*(Al*spd[0]*Bx + Cl) - (Al + Gl)*(spd[0]*pt + Rl.E);

      vxl = ( Bx*(Al*Bx + spd[0]*Cl) - (Al + Gl)*(pt + Rl.Mx) ) / Xl;

      tmp = Cl + Bx*(spd[0]*Rl.Mx - Rl.E);
      vyl = ( Ql*Rl.My + Rl.By*tmp ) / Xl;
      vzl = ( Ql*Rl.Mz + Rl.Bz*tmp ) / Xl;

      /* right side */
      Ar = Rr.Mx - spd[4]*Rr.E + pt*(1.0 - SQR(spd[4]));
      Gr = Rr.By*Rr.By + Rr.Bz*Rr.Bz;
      Cr = Rr.My*Rr.By + Rr.Mz*Rr.Bz;
      Qr = -Ar - Gr + SQR(Bx)*(1.0 - SQR(spd[4]));
      Xr = Bx*(Ar*spd[4]*Bx + Cr) - (Ar + Gr)*(spd[4]*pt + Rr.E);
      
      vxr = ( Bx*(Ar*Bx + spd[4]*Cr) - (Ar + Gr)*(pt + Rr.Mx) ) / Xr;

      tmp = Cr + Bx*(spd[4]*Rr.Mx - Rr.E);
      vyr = ( Qr*Rr.My + Rr.By*tmp ) / Xr;
      vzr = ( Qr*Rr.Mz + Rr.Bz*tmp ) / Xr;

      /*printf("vxl: %e\n",vxl);
        printf("vxr: %e\n",vxr); */

      /*** Mignone Eq 21 ***/

      /* left side */
      ovspdMinusVxl = 1.0 / (spd[0] - vxl);
      /* printf("spd[0]: %f\n",spd[0]);
         printf("vxl: %f\n",vxl); */

      Ual.By = (Rl.By - Bx*vyl) * ovspdMinusVxl;
      Ual.Bz = (Rl.Bz - Bx*vzl) * ovspdMinusVxl;

      /* right side */
      ovspdMinusVxr = 1.0 / (spd[4] - vxr);
      /* printf("spd[4]: %f\n",spd[4]);
         printf("vxr: %f\n",vxr); */

      Uar.By = (Rr.By - Bx*vyr) * ovspdMinusVxr;
      Uar.Bz = (Rr.Bz - Bx*vzr) * ovspdMinusVxr;

      /*** Mignone Eq 31 ***/
      
      /* left side */
      vDotRml = vxl*Rl.Mx + vyl*Rl.My + vzl*Rl.Mz;
      wal = pt + (Rl.E - vDotRml) * ovspdMinusVxl;

      /* right side */
      vDotRmr = vxr*Rr.Mx + vyr*Rr.My + vzr*Rr.Mz;
      war = pt + (Rr.E - vDotRmr) * ovspdMinusVxr;

      /*--- Step 5.2 -----------------------------------------------------------
       * Comput K-vectors according to Mignone Eq 43 and Bc-vector according to
       * Mignone Eq 45
       */
      
      if(Bx > 0)
         signBx = 1.0;
      else
         signBx = -1.0;

      /* Mignone Eq 35 */
      etal = -signBx * sqrt(wal);
      etar = signBx * sqrt(war);

      /* Mignone Eq 43 */

      ovKldenom = 1.0 / (spd[0]*pt + Rl.E + Bx*etal);
      ovKrdenom = 1.0 / (spd[4]*pt + Rr.E + Bx*etar);

      Kxl = (Rl.Mx + pt + (Bx*spd[0])*etal) * ovKldenom;
      Kxr = (Rr.Mx + pt + (Bx*spd[4])*etar) * ovKrdenom;

      Kyl = (Rl.My + Rl.By*etal) * ovKldenom;
      Kyr = (Rr.My + Rr.By*etar) * ovKrdenom;

      Kzl = (Rl.Mz + Rl.Bz*etal) * ovKldenom;
      Kzr = (Rr.Mz + Rr.Bz*etar) * ovKrdenom;

      /* Mignone states that lamda_a = K_a^x just before Eq 37 */

      spd[1] = Kxl;
      spd[3] = Kxr;

      
      /* Mignone Eq 45 */

      if(Kxl != Kxr){

         spdMinVxl = spd[1] - vxl;
         spdMinVxr = spd[3] - vxr;
         ovs3Mns1 = 1.0 / (spd[3] - spd[1]); 

         aLy = Ual.By*spdMinVxl + Bx*vyl;
         aRy = Uar.By*spdMinVxr + Bx*vyr;

         aLz = Ual.Bz*spdMinVxl + Bx*vzl;
         aRz = Uar.Bz*spdMinVxr + Bx*vzr;

         Bcy = (aRy - aLy) * ovs3Mns1;
         Bcz = (aRz - aLz) * ovs3Mns1;

         /*--- Step 5.3 -----------------------------------------------------------
          * Use Mignone Eq 48 to find next improved total pressure value
          */
      
         dKx  = Kxr - Kxl;

         Klsq = SQR(Kxl) + SQR(Kyl) + SQR(Kzl);
         Krsq = SQR(Kxr) + SQR(Kyr) + SQR(Kzr);

         KldotB = (Kxl*Bx + Kyl*Bcy + Kzl*Bcz);
         KrdotB = (Kxr*Bx + Kyr*Bcy + Kzr*Bcz);

         KldotBhc = dKx * KldotB;
         KrdotBhc = dKx * KrdotB;

         /* Mignone Eq 49 */
         Yl = (1.0 - Klsq) / (etal*dKx - KldotBhc);
         Yr = (1.0 - Krsq) / (etar*dKx - KrdotBhc);

         /*printf("dKx: %e\n",dKx);
            printf("Yl: %e\n",Yl);
            printf("Yr: %e\n",Yr); */
      }
      else{
         break;
      }

      /* Mignone Eq 48 */
      fp = dKx * (1.0 - Bx*(Yr - Yl));


      /* printf("Status: %d\n",status); */

      if(status == 0){ /* first loop */
         if(fabs(fp) < maxError)
            break;
         fpold = fp;
         ptold = pt;
         pt = (1.0 + jump)*pt;
         status = 1;
      }
      else if(status == 1){ /* second loop */

         /* if already bounded, then procede with secant */
         if( (fp>0 && fpold<0) || (fp<0 && fpold>0) ){
            tmp = pt;
            pt = pt - (pt - ptold)/(fp - fpold) * fp;
            ptold = tmp;
            fpold = fp;
            error = fabs(fp);
            status = 2;
         }

         /* if not bounded, bound the zero */
         else{
            slope = (fp - fpold)/(pt - ptold);
            ptold = pt;
            fpold = fp;

            /* sign of jump */
            if(fp*slope > 0)
               tmp = -1.0;
            else
               tmp = 1.0;
            
            pt = (1.0 + tmp*jump)*pt;
         }
      }
      else{ /* status == 2, secant method */
         tmp = pt;
         pt = pt - (pt - ptold)/(fp - fpold) * fp;
         ptold = tmp;
         fpold = fp;
         error = fabs(fp);
      }

      /* printf("pt: %f\n",pt);
      printf("fp: %f\n",fp);
      printf("error: %f\n\n",error); */
   }

/*--- Step 6.-------------------------------------------------------------------
 * Now begin back subsitution and returning fluxes
 */

   /* Mignone Eq 32-34 */

   /* left side */
   vDotBl = vxl*Bx + vyl*Ual.By + vzl*Ual.Bz;

   Ual.d  = Rl.d * ovspdMinusVxl;
   /*printf("Rl_D: %f\n",Rl.d);
     printf("den: %f\n",ovspdMinusVxl);*/
   Ual.E  = (Rl.E + pt*vxl - vDotBl*Bx) * ovspdMinusVxl;
   Ual.Mx = (Ual.E + pt)*vxl - vDotBl*Bx;
   Ual.My = (Ual.E + pt)*vyl - vDotBl*Ual.By;
   Ual.Mz = (Ual.E + pt)*vzl - vDotBl*Ual.Bz;

   /* right side */
   vDotBr = vxr*Bx + vyr*Uar.By + vzr*Uar.Bz;
      
   Uar.d  = Rr.d * ovspdMinusVxr;
   Uar.E  = (Rr.E + pt*vxr - vDotBr*Bx) * ovspdMinusVxr;
   Uar.Mx = (Uar.E + pt)*vxr - vDotBr*Bx;
   Uar.My = (Uar.E + pt)*vyr - vDotBr*Uar.By;
   Uar.Mz = (Uar.E + pt)*vzr - vDotBr*Uar.Bz;

   /* compute FaL and FaR according to Mignone Eq 16.5;
      Note, the equation as listed in paper is slightly
      incorrect.  lamda term should not have subscript */

   /* FaL */

   Fal.d  = Fl.d  + spd[0]*(Ual.d  - Ul.d );
   Fal.Mx = Fl.Mx + spd[0]*(Ual.Mx - Ul.Mx);
   Fal.My = Fl.My + spd[0]*(Ual.My - Ul.My);
   Fal.Mz = Fl.Mz + spd[0]*(Ual.Mz - Ul.Mz);
   Fal.E  = Fl.E  + spd[0]*(Ual.E  - Ul.E );
   Fal.By = Fl.By + spd[0]*(Ual.By - Ul.By);
   Fal.Bz = Fl.Bz + spd[0]*(Ual.Bz - Ul.Bz);

   /* FaR */

   Far.d  = Fr.d  + spd[4]*(Uar.d  - Ur.d );
   Far.Mx = Fr.Mx + spd[4]*(Uar.Mx - Ur.Mx);
   Far.My = Fr.My + spd[4]*(Uar.My - Ur.My);
   Far.Mz = Fr.Mz + spd[4]*(Uar.Mz - Ur.Mz);
   Far.E  = Fr.E  + spd[4]*(Uar.E  - Ur.E );
   Far.By = Fr.By + spd[4]*(Uar.By - Ur.By);
   Far.Bz = Fr.Bz + spd[4]*(Uar.Bz - Ur.Bz);


   /* now compute velocities in Mignone Eq 47 */

   mKlsq = 1.0 - Klsq;
   mKrsq = 1.0 - Krsq;

   ovEtalmKdB = 1.0 / (etal - KldotB);
   ovEtarmKdB = 1.0 / (etar - KrdotB);

   vxcl = Kxl - Bx * mKlsq * ovEtalmKdB;
   vycl = Kyl - Bcy * mKlsq * ovEtalmKdB;
   vzcl = Kzl - Bcz * mKlsq * ovEtalmKdB;

   vxcr = Kxr - Bx * mKrsq * ovEtarmKdB;
   vycr = Kxr - Bcy * mKrsq * ovEtarmKdB;
   vzcr = Kxr - Bcz * mKrsq * ovEtarmKdB;

   /*printf("Klsq: %f\n",Klsq);
   printf("Krsq: %f\n",Krsq);
   printf("Eta_L: %f\n",etal);
   printf("Eta_R: %f\n",etar);
   printf("KldotB: %f\n",KldotB);
   printf("KrdotB: %f\n\n",KrdotB);

   printf("Kxl: %f\n",Kxl);
   printf("Kxr: %f\n",Kxr);
   printf("Kyl: %f\n",Kyl);
   printf("Kyr: %f\n",Kyr);
   printf("Kzl: %f\n",Kzl);
   printf("Kzr: %f\n\n",Kzr);*/

   /*printf("Bx: %f\n",Bx);
   printf("Bcy: %f\n",Bcy);
   printf("Bcz: %f\n\n",Bcz);*/

   /*printf("vxcl: %e\n",vxcl);
   printf("vxcr: %e\n",vxcr);
   printf("vycl: %e\n",vycl);
   printf("vycr: %e\n",vycr);
   printf("vzcl: %e\n",vzcl);
   printf("vzcr: %e\n\n",vzcr);*/

   /* now compute UcL and UcR according to Mignone Eqs 50-52 */

   Ucl.By = Bcy;
   Ucr.By = Bcy;

   Ucl.Bz = Bcz;
   Ucr.Bz = Bcz;

   ovldalMvxcl = 1.0 / (spd[1] - vxcl);
   ovldarMvxcr = 1.0 / (spd[3] - vxcr);

   Ucl.d = Ual.d * (spd[1] - vxl) * ovldalMvxcl;
   Ucr.d = Uar.d * (spd[3] - vxr) * ovldarMvxcr;

   vclDotBc = vxcl*Bx + vycl*Ucl.By + vzcl*Ucl.Bz;
   vcrDotBc = vxcr*Bx + vycr*Ucr.By + vzcr*Ucr.Bz;

   Ucl.E = (spd[1]*Ual.E - Ual.Mx + pt*vxcl - vclDotBc*Bx) * ovldalMvxcl;
   Ucr.E = (spd[3]*Uar.E - Uar.Mx + pt*vxcr - vcrDotBc*Bx) * ovldarMvxcr;

   eclPpt = Ucl.E + pt;
   ecrPpt = Ucr.E + pt;

   Ucl.Mx = eclPpt*vxcl - vclDotBc*Bx;
   Ucr.Mx = ecrPpt*vxcr - vcrDotBc*Bx;

   Ucl.My = eclPpt*vycl - vclDotBc*Ucl.By;
   Ucr.My = ecrPpt*vycr - vcrDotBc*Ucr.By;

   Ucl.Mz = eclPpt*vzcl - vclDotBc*Ucl.Bz;
   Ucr.Mz = ecrPpt*vzcr - vcrDotBc*Ucr.Bz;

   /* Now use Mignone equation 47 to solve for lamda_c = Vx */

   spd[2] = (vxcl+vxcr)*0.5;

   /* printf("pt: %f\n",pt); */
   for(i=0; i<5; i++)
      printf("spd[%d]: %e\n",i,spd[i]);
   printf("\n");

   /* Check for consistancy according to Mignone Eq 54;
      If inconsistant, use HLL flux */

   if( (wl - pt <= 0) || (wr - pt <= 0) || 
       (vxcl - Kxl <= -1.e-6) || (Kxr - vxcr <= -1.e-6) ||
       (spd[1] - vxl >= 0.0) || (spd[3] - vxr <= 0.0) ||
       (spd[1] - spd[0] <= -1.e-6) || (spd[4] - spd[3] <= -1.e-6)){
      
      printf("Flux_HLL\n");

      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;

      return;
   }

   if(spd[1] >= -1.0e-6){ /* aL flux */

      printf("Flux_aL\n");

      /*printf("Ual\n");
        printCons1D(&Ual);*/

      pFlux->d = Fal.d;
      pFlux->Mx = Fal.Mx;
      pFlux->My = Fal.My;
      pFlux->Mz = Fal.Mz;
      pFlux->E = Fal.E;
      pFlux->By = Fal.By;
      pFlux->Bz = Fal.Bz;

      return;
   }
   else if(spd[3] <= 1.0e-6){ /* aR flux */

      printf("Flux_aR\n");
      
      pFlux->d = Far.d;
      pFlux->Mx = Far.Mx;
      pFlux->My = Far.My;
      pFlux->Mz = Far.Mz;
      pFlux->E = Far.E;
      pFlux->By = Far.By;
      pFlux->Bz = Far.Bz;
      
      return;
   }
   else{
      if(spd[2] > 0.0){ /* cL flux */

         printf("Flux_cL\n");

         pFlux->d  = Fal.d  + spd[2]*(Ucl.d  - Ual.d );
         pFlux->Mx = Fal.Mx + spd[2]*(Ucl.Mx - Ual.Mx);
         pFlux->My = Fal.My + spd[2]*(Ucl.My - Ual.My);
         pFlux->Mz = Fal.Mz + spd[2]*(Ucl.Mz - Ual.Mz);
         pFlux->E  = Fal.E  + spd[2]*(Ucl.E  - Ual.E );
         pFlux->By = Fal.By + spd[2]*(Ucl.By - Ual.By);
         pFlux->Bz = Fal.Bz + spd[2]*(Ucl.Bz - Ual.Bz);
         
         return;
      }
      else{ /* cR flux */

         printf("Flux_cR\n");

         pFlux->d  = Far.d  + spd[2]*(Ucr.d  - Uar.d );
         pFlux->Mx = Far.Mx + spd[2]*(Ucr.Mx - Uar.Mx);
         pFlux->My = Far.My + spd[2]*(Ucr.My - Uar.My);
         pFlux->Mz = Far.Mz + spd[2]*(Ucr.Mz - Uar.Mz);
         pFlux->E  = Far.E  + spd[2]*(Ucr.E  - Uar.E );
         pFlux->By = Far.By + spd[2]*(Ucr.By - Uar.By);
         pFlux->Bz = Far.Bz + spd[2]*(Ucr.Bz - Uar.Bz);

         return;
      }
   }
}

/* Use Mignone(2005) */
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
   blsq = Brsq / gammar2 + SQR(vDotBr);

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

      nl = solveQuartic(al,rl,1.0e-10);
      nr = solveQuartic(ar,rr,1.0e-10);

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
   
   if(eSq >= 0){
      E = sqrt(eSq);
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

#endif /* HLLD_FLUX */
#endif /* SPECIAL_RELATIVITY */
