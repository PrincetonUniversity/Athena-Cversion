#include "../copyright.h"
/*==============================================================================
 * FILE: hllc_sr.c
 *
 * PURPOSE: Computes 1D fluxes using the relativistic HLLC Riemann solver, 
 *   an extension of the HLLE fluxes to include the contact wave.  Currently 
 *   only works for hydrodynamics.  For an extension to MHD, see hlld_sr.c
 *
 * REFERENCES:
 *   A. Mignone and G. Bodo, "An HLLC Riemann solver for relativistic flows",
 *   Mon. Not. R. Astron. Soc. 364, 126-136 (2005)
 *
 *   A. Mignone and G. Bodo, "An HLLC Solver for Relativistic Flows - II:
 *   Magnetohydrodynamics", arxiv:astro-ph/0601640v1 (2006)
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

#ifdef HLLC_FLUX
#ifdef SPECIAL_RELATIVITY

#ifdef MHD
#error : The HLLC flux only works for hydro.
#endif

#if (NSCALARS > 0)
#error : The SR HLLC flux does not work with passive scalars.
#endif

void printCons1D(const Cons1D* c);
  void printPrim1D(const Prim1D* p);

void printCons1D(const Cons1D* c){
   printf("d:  %e\n",c->d);
   printf("E:  %e\n",c->E);
   printf("Mx: %e\n",c->Mx);
   printf("My: %e\n",c->My);
   printf("Mz: %e\n",c->Mz);
   printf("\n");
}

void printPrim1D(const Prim1D* p){
   printf("d:  %e\n",p->d);
   printf("P:  %e\n",p->P);
   printf("Vx: %e\n",p->Vx);
   printf("Vy: %e\n",p->Vy);
   printf("Vz: %e\n",p->Vz);
   printf("\n");
}

/*----------------------------------------------------------------------------*/
/* fluxes
 *   Input Arguments:
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *     Wl,Wr = L/R-states of PRIMITIVE variables at cell interface 
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1D Ul, const Cons1D Ur,
            const Prim1D Wl, const Prim1D Wr, Cons1D *pFlux)
{
  Cons1D Fl,Fr,Fhll,Uhll,Usl,Usr;
  Real rhl, rhr, csl, csr, cslsq, csrsq, vsql, vsqr, gammasql, gammasqr;
  Real ssl, ssr, radl, radr, lmdapl, lmdapr, lmdaml, lmdamr, lmdatlmda;
  Real lmdal,lmdar; /* Left and Right wave speeds */
  Real lmdas; /* Contact wave speed */
  Real ovlrmll;
  Real AL,AR,BL,BR,a,b,c;
  Real bsq,rad;
  Real den,ps; /* Pressure in inner region */

  /*printf("Wl\n");
  printPrim1D(&Wl);
  printf("Wr\n");
  printPrim1D(&Wr);
  printf("Ul\n");
  printCons1D(&Ul);
  printf("Ur\n");
  printCons1D(&Ur);*/

/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */

  rhl = Wl.d + Wl.P * Gamma / Gamma_1; /* Mignone Eq 3.5 */
  rhr = Wr.d + Wr.P * Gamma / Gamma_1;

  csl = sqrt(Gamma * Wl.P / rhl); /* Mignone Eq 4 */
  csr = sqrt(Gamma * Wr.P / rhr);

  cslsq = SQR(csl);
  csrsq = SQR(csr);

  vsql = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
  vsqr = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);

  gammasql = 1.0 / (1.0 - vsql);
  gammasqr = 1.0 / (1.0 - vsqr);

  ssl = cslsq / ( gammasql * (1.0 - cslsq) ); /* Mignone Eq 22.5 */
  ssr = csrsq / ( gammasqr * (1.0 - csrsq) );

  radl = sqrt( ssl*(1.0-SQR(Wl.Vx)+ssl) ); /* Mignone Eq 23 (radical part) */
  radr = sqrt( ssr*(1.0-SQR(Wr.Vx)+ssr) );

  lmdapl = (Wl.Vx + radl) / (1.0 + ssl); /* Mignone Eq 23 */
  lmdapr = (Wr.Vx + radr) / (1.0 + ssr);
  lmdaml = (Wl.Vx - radl) / (1.0 + ssl);
  lmdamr = (Wr.Vx - radr) / (1.0 + ssr);

  lmdal = MIN(lmdaml, lmdamr); /* Mignone Eq 21 */
  lmdar = MAX(lmdapl, lmdapr);
  

/*--- Step 2. ------------------------------------------------------------------
 * Compute L/R fluxes according to Mignone 2
 */

  Fl.d  = Ul.d * Wl.Vx;
  Fl.Mx = Ul.Mx * Wl.Vx + Wl.P;
  Fl.My = Ul.My * Wl.Vx;
  Fl.Mz = Ul.Mz * Wl.Vx;
  Fl.E  = Ul.Mx;

  Fr.d  = Ur.d * Wr.Vx;
  Fr.Mx = Ur.Mx * Wr.Vx + Wr.P;
  Fr.My = Ur.My * Wr.Vx;
  Fr.Mz = Ur.Mz * Wr.Vx;
  Fr.E  = Ur.Mx;

  /*printf("F_L\n");
  printCons1D(&Fl);
  printf("F_R\n");
  printCons1D(&Fr);*/

/*--- Step 3. ------------------------------------------------------------------
 * Compute HLL flux using Mignone Eq 11 (necessary for computing lmdas (Eq 18)
 * Compute HLL conserved quantities using Mignone eq 9
 */

  ovlrmll = 1.0 / ( lmdar - lmdal );
  lmdatlmda = lmdal*lmdar;

  Fhll.d  = (lmdar*Fl.d  - lmdal*Fr.d  + lmdatlmda * (Ur.d  - Ul.d) ) * ovlrmll;
  Fhll.Mx = (lmdar*Fl.Mx - lmdal*Fr.Mx + lmdatlmda * (Ur.Mx - Ul.Mx)) * ovlrmll;
  Fhll.My = (lmdar*Fl.My - lmdal*Fr.My + lmdatlmda * (Ur.My - Ul.My)) * ovlrmll;
  Fhll.Mz = (lmdar*Fl.Mz - lmdal*Fr.Mz + lmdatlmda * (Ur.Mz - Ul.Mz)) * ovlrmll;
  Fhll.E  = (lmdar*Fl.E  - lmdal*Fr.E  + lmdatlmda * (Ur.E  - Ul.E )) * ovlrmll;

  Uhll.d  = (lmdar * Ur.d  - lmdal * Ul.d  + Fl.d  - Fr.d ) * ovlrmll;
  Uhll.Mx = (lmdar * Ur.Mx - lmdal * Ul.Mx + Fl.Mx - Fr.Mx) * ovlrmll;
  Uhll.My = (lmdar * Ur.My - lmdal * Ul.My + Fl.My - Fr.My) * ovlrmll;
  Uhll.Mz = (lmdar * Ur.Mz - lmdal * Ul.Mz + Fl.Mz - Fr.Mz) * ovlrmll;
  Uhll.E  = (lmdar * Ur.E  - lmdal * Ul.E  + Fl.E  - Fr.E ) * ovlrmll;

/*  printf("Uhll\n");
  printCons1D(&Uhll);
  printf("Fhll\n");
  printCons1D(&Fhll);*/

/*--- Step 4. ------------------------------------------------------------------
 * Compute contact wave speed using larger root from Mignone Eq 18
 * Physical root is the root with the minus sign
 */

  /* quadratic formula calculation */

  if(fabs(Fhll.E) > (TINY_NUMBER)){
     b   = -(Uhll.E + Fhll.Mx);
     bsq = b * b;
     rad = sqrt(bsq - 4.0 * Fhll.E * Uhll.Mx);
     lmdas = (-b - rad) / (2.0 * Fhll.E);
     /*printf("%e\n",Fhll.E);*/
  }
  else{
     /*printf("HERE: %e\n",Fhll.E);*/
     lmdas = Uhll.Mx / (Uhll.E + Fhll.Mx);
  }

  /*********************************/
  /*printf("lamdal: %e\n",lmdal);
  printf("lamdas: %e\n",lmdas);
  printf("lamdar: %e\n\n",lmdar);*/

/*--- Step 5. ------------------------------------------------------------------
 * Determine intercell flux according to Mignone 13
 */

  if( lmdal >= 0.0){ /* Fl */
     /* intercell flux is left flux */
     pFlux->d  = Fl.d;
     pFlux->Mx = Fl.Mx;
     pFlux->My = Fl.My;
     pFlux->Mz = Fl.Mz;
     pFlux->E  = Fl.E;

     return;
  }
  else if( lmdas == 0.0 ){
     pFlux->d = 0.0;
     pFlux->Mx = 0.0;
     pFlux->My = 0.0;
     pFlux->Mz = 0.0;
     pFlux->E = 0.0;
  }
  else if( lmdas > 0.0){ /* Fls */

     /* Mignone 2006 Eq 48 */
     ps = -Fhll.E*lmdas + Fhll.Mx;

     /* now calculate Usl with Mignone Eq 16 */
     den = 1.0 / (lmdal - lmdas);

     Usl.d  = Ul.d * (lmdal - Wl.Vx) * den;
     Usl.Mx = (Ul.Mx * (lmdal - Wl.Vx) + ps - Wl.P) * den;
     Usl.My = Ul.My * (lmdal - Wl.Vx) * den;
     Usl.Mz = Ul.Mz * (lmdal - Wl.Vx) * den;
     Usl.E  = (Ul.E * (lmdal - Wl.Vx) + ps * lmdas - Wl.P * Wl.Vx) * den;

     /*printf("ps: %f\n",ps);
     printf("Fl\n");
     printCons1D(&Fl);
     printf("Usl\n");
     printCons1D(&Usl);*/

     /* now calculate Fsr using Mignone Eq 14 */

     pFlux->d  = lmdal*(Usl.d  - Ul.d ) + Fl.d;
     pFlux->Mx = lmdal*(Usl.Mx - Ul.Mx) + Fl.Mx;
     pFlux->My = lmdal*(Usl.My - Ul.My) + Fl.My;
     pFlux->Mz = lmdal*(Usl.Mz - Ul.Mz) + Fl.Mz;
     pFlux->E  = lmdal*(Usl.E  - Ul.E ) + Fl.E;

     return;
  }
  else if( lmdar >= 0.0){ /* Frs */

     /* Mignone 2006 Eq 48 */
     ps = -Fhll.E*lmdas + Fhll.Mx;

     /* now calculate Usr with Mignone Eq 16 */
     den = 1.0 / (lmdar - lmdas);

     Usr.d  = Ur.d * (lmdar - Wr.Vx) * den;
     Usr.Mx = (Ur.Mx * (lmdar - Wr.Vx) + ps - Wr.P) * den;
     Usr.My = Ur.My * (lmdar - Wr.Vx) * den;
     Usr.Mz = Ur.Mz * (lmdar - Wr.Vx) * den;
     Usr.E  = (Ur.E * (lmdar - Wr.Vx) + ps * lmdas - Wr.P * Wr.Vx) * den;

     /*printf("ps: %f\n",ps);
     printf("Fr\n");
     printCons1D(&Fr);
     printf("Usr\n");
     printCons1D(&Usr); */

     /* now calculate Fsr using Mignone Eq 14 */

     pFlux->d  = lmdar*(Usr.d  - Ur.d ) + Fr.d;
     pFlux->Mx = lmdar*(Usr.Mx - Ur.Mx) + Fr.Mx;
     pFlux->My = lmdar*(Usr.My - Ur.My) + Fr.My;
     pFlux->Mz = lmdar*(Usr.Mz - Ur.Mz) + Fr.Mz;
     pFlux->E  = lmdar*(Usr.E  - Ur.E ) + Fr.E;

     return;
  }
  else{ /* Fr */
     /* intercell flux is right flux */
     pFlux->d  = Fr.d;
     pFlux->Mx = Fr.Mx;
     pFlux->My = Fr.My;
     pFlux->Mz = Fr.Mz;
     pFlux->E  = Fr.E;

     return;
  }
  
  /* need to deal with scalar fluxes */
}



#endif /* HLLC_FLUX */
#endif /* SPECIAL_RELATIVITY */
