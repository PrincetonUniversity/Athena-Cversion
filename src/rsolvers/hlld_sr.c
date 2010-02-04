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

#if HYDRO
#error : The SR HLLD flux does not work for gas=hydro
#endif

#define MAX_ITER 20

typedef struct RIEMANN_STATE{
  int fail;
  Real vx, vy, vz;
  Real Bx, By, Bz;
  Real Kx, Ky, Kz, K2;
  Real w, eta, p, rho;
  Cons1DS U, R;
  Real S, Sa;
} Riemann_State;

void zeroCons1D(Cons1DS *U);
void printCons1D(const Cons1DS *U);
void printPrim1D(const Prim1DS *W);

/* computes left/right fluxes from left/right states */
void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p);
int getHLLCpress (const Prim1DS Wl,   const Prim1DS Wr,
		  const Cons1DS Uhll, const Cons1DS Fhll,
		  const Real Bxi,     Real *ps); 
void getHLLCflux (const Cons1DS Ul,   const Cons1DS Ur,
		  const Prim1DS Wl,   const Prim1DS Wr,
		  const Real Sl,      const Real Sr,
		  const Cons1DS Fl,   const Cons1DS Fr,
		  const Cons1DS Uhll, const Prim1DS Whll,
		  const Cons1DS Fhll, const Real Bxi, Cons1DS *pFlux);
Real Fstar (const Real Bx, const Real p, Riemann_State *PaL, Riemann_State *PaR, Real *Sc);
int getRiemannState (const int side, const Real Bx, const Real p, Riemann_State *Pv);
void getAstate (const Real Bx    , const Real p,
		Riemann_State *Pa, Cons1DS *U  );
void getCstate (const Real Bx     , const Real p,
		Riemann_State *PaL, Riemann_State *PaR,
		Cons1DS *UaL,       Cons1DS *UaR,       Cons1DS *Uc,
		Real *SaL,          Real *SaR,          Real *Sc);
void getPtot (const Real Bx, const Prim1DS W, Real *pt);
void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high);
void getMaxSignalSpeeds_echo(const Prim1DS Wl, const Prim1DS Wr,
			     const Real Bx, Real* low, Real* high);
void getVChar_echo(const Prim1DS W, const Real Bx, Real* lml, Real* lmr);
void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lml, Real* lmr);
void getSoundSpeed2(const Prim1DS W, Real *cs2, Real *h);
/* solves quartic equation defined by a and returns roots in root
 * returns the number of Real roots 
 * error specifies an accuracy
 * currently force four Real solutions b/c it's physical */
int QUARTIC (Real b, Real c, Real d, Real e, Real z[]);

/* solves cubic equation defined by a and stores roots in root
 * returns number of Real roots */
int CUBIC(Real b, Real c, Real d, Real z[]);

/* solves quartic equation defined by a and returns roots in root
 * returns the number of Real roots 
 * error specifies an accuracy
 * currently force four Real solutions b/c it's physical */
int QUARTIC (Real b, Real c, Real d, Real e, Real z[]);

/* solves cubic equation defined by a and stores roots in root
 * returns number of Real roots */
int CUBIC(Real b, Real c, Real d, Real z[]);

#ifdef MHD

void zeroCons1D(Cons1DS *U){
  U->d = 0.0;
  U->E = 0.0;
  U->Mx = 0.0;
  U->My = 0.0;
  U->Mz = 0.0;
  U->By = 0.0;
  U->Bz = 0.0;
}

/* functions for printing conserved/primitive vectors */
void printCons1D(const Cons1DS *U){
  printf("d:  %.6e\n",U->d);
  printf("E:  %.6e\n",U->E);
  printf("Mx: %.6e\n",U->Mx);
  printf("My: %.6e\n",U->My);
  printf("Mz: %.6e\n",U->Mz);
  printf("By: %.6e\n",U->By);
  printf("Bz: %.6e\n",U->Bz);
  printf("\n");
}

void printPrim1D(const Prim1DS *W){
  printf("d:  %.6e\n",W->d);
  printf("P:  %.6e\n",W->P);
  printf("Vx: %.6e\n",W->Vx);
  printf("Vy: %.6e\n",W->Vy);
  printf("Vz: %.6e\n",W->Vz);
  printf("By: %.6e\n",W->By);
  printf("Bz: %.6e\n",W->Bz);
  printf("\n");
}

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Fl,Fr;
  Cons1DS Fhll,Uhll;
  Cons1DS UaL,UaR,Uc;
  Prim1DS WaL,WaR,Wc,Whll;
  Riemann_State PaL, PaR;

  Real Sl, Sr, Pl, Pr;
  Real Sla, Sra;
  Real SaL, SaR, Sc;
  Real V2l, V2r, V2c;
  Real p,p0,f,f0;
  Real dPstep;
  Real dS_1, scrh;
  Real tol = 1.0e-3;

  int p_incs,switch_to_hllc,nr_success;
	
  /*--- Step 1a. ------------------------------------------------------------------
   * Compute the exact (estimated) max and min wave speeds via method from
   * PLUTO (ECHO).
   */
  switch_to_hllc = 0;
  getMaxSignalSpeeds_pluto(Wl,Wr,Bxi,&Sl,&Sr);
  getMaxSignalSpeeds_echo(Wl,Wr,Bxi,&Sla,&Sra);
  /*printf("[hlld_sr_mhd]: Exact Wave speeds     %10.4e %10.4e\n",Sl,Sr);
    printf("[hlld_sr_mhd]: Estimated Wave speeds %10.4e %10.4e\n",Sla,Sra);*/
	
  if (Sla != Sla) {
    printf("[hlld_sr_mhd]: NaN in estimated Sl %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra =  1.0;
  }
	
  if (Sra != Sra) {
    printf("[hlld_sr_mhd]: NaN in estimated Sr %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
	
  if (Sla < -1.0) {
    printf("[hlld_sr]: Superluminal estimated Sl %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
  if (Sra > 1.0) {
    printf("[hlld_sr]: Superluminal estimated Sr %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
	
  if (Sl != Sl) {
    switch_to_hllc = 1;
    printf("[hlld_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    switch_to_hllc = 1;
    printf("[hlld_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    switch_to_hllc = 1;
    printf("[hlld_sr]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    switch_to_hllc = 1;
    printf("[hlld_sr]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

  /*--- Step 1b. ------------------------------------------------------------------
   * If we've ended up with maximal wave speeds from the exact solution
   * then we try using the estimated wave speeds
   */
  if (Sl ==-1.0){
    switch_to_hllc = 1;
    Sl = Sla;
  }
  if (Sr == 1.0){
    switch_to_hllc = 1;
    Sr = Sra;
  }

  /*
    if (Sla < Sl) {
    Sl = Sla;
    }
	
    if (Sra > Sr) {
    Sl = Sla;
    }
  */

  /*--- Step 2a. ------------------------------------------------------------------
   * Comput L/R fluxes
   */
  flux_LR(Ul,Wl,&Fl,Bxi,&Pl);
  flux_LR(Ur,Wr,&Fr,Bxi,&Pr);
	
  /*--- Step 2b. ------------------------------------------------------------------
   * Construct HLL fluxes & average state to use as initial guesses and fall back
   */
  dS_1 = 1.0/(Sr - Sl);
	
  Uhll.d  = (Sr*Ur.d  - Sl*Ul.d  + Fl.d  - Fr.d ) * dS_1;
  Uhll.Mx = (Sr*Ur.Mx - Sl*Ul.Mx + Fl.Mx - Fr.Mx) * dS_1;
  Uhll.My = (Sr*Ur.My - Sl*Ul.My + Fl.My - Fr.My) * dS_1;
  Uhll.Mz = (Sr*Ur.Mz - Sl*Ul.Mz + Fl.Mz - Fr.Mz) * dS_1;
  Uhll.E  = (Sr*Ur.E  - Sl*Ul.E  + Fl.E  - Fr.E ) * dS_1;
  Uhll.By = (Sr*Ur.By - Sl*Ul.By + Fl.By - Fr.By) * dS_1;
  Uhll.Bz = (Sr*Ur.Bz - Sl*Ul.Bz + Fl.Bz - Fr.Bz) * dS_1;
	
  Whll = Cons1D_to_Prim1D (&Uhll,&Bxi);
  /*     V2l = SQR(Whll.Vx) + SQR(Whll.Vy) + SQR(Whll.Vz);
 	if (Whll.P < 0 || Whll.d < 0 || V2l > 1.0){
	printf("[hlld_sr_mhd]: hll average state has Pl = %10.4e, dl = %10.4e, Vsq = %10.4e\n",Whll.P,Whll.d,V2l);

	Sl = -1.0;
	Sr = 1.0;
	dS_1 = 1.0/(Sr - Sl);
		
	Fhll.d  = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
	Fhll.Mx = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
	Fhll.My = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
	Fhll.Mz = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
	Fhll.E  = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;
	Fhll.By = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
	Fhll.Bz = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;
		
	pFlux->d = Fhll.d;
	pFlux->Mx = Fhll.Mx;
	pFlux->My = Fhll.My;
	pFlux->Mz = Fhll.Mz;
	pFlux->E = Fhll.E;
	pFlux->By = Fhll.By;
	pFlux->Bz = Fhll.Bz;
		
	return;
	}
  */
	
  Fhll.d  = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
  Fhll.Mx = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
  Fhll.My = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
  Fhll.Mz = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
  Fhll.E  = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;
  Fhll.By = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
  Fhll.Bz = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;
	
  /* if the normal component of the field is zero, then use the hllc solver */
  if (fabs(Bxi) <= 1.0e-12){
    getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);
    return;
  }

  /*	if (switch_to_hllc){
	printf("[hlld_sr_mhd]: Switching to hllc flux\n");
	getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);
	return;
	}*/
	
  /*--- Step 2b. ------------------------------------------------------------------
   * If we've ended up with a maximal wave speed, then just compute the LF flux
   * and return
   */
  if (Sl == -1.0) {
    printf("[hlld_sr_mhd]: Switching to hll flux\n");
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;

    /*    if (pFlux->d != pFlux->d) {
      printf("[hlld_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
      printf("[hlld_sr_mhd]: NaN in hllc density flux\n");
      }*/			
		
    return;
  }
	
  if (Sr == 1.0) {
    /*printf("[hlld_sr_mhd]: Switching to hll flux\n");*/
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;

    /*if (pFlux->d != pFlux->d) {
      printf("[hlld_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
      printf("[hlld_sr_mhd]: NaN in hllc density flux\n");
      }*/			
		
    return;
  }
	
  if(Sl >= 0.0){
    /*printf("[hlld_sr_mhd]: State 1 \n");*/
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;
    pFlux->By = Fl.By;
    pFlux->Bz = Fl.Bz;
		
    /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      printf("[hlld_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
      printf("[hlld_sr_mhd]: NaN in hllc density flux\n");
      }*/
		
    return;
		
  }
  if(Sr <= 0.0){
    /*printf("[hlld_sr_mhd]: State 6 \n");*/
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;
    pFlux->By = Fr.By;
    pFlux->Bz = Fr.Bz;
		
    /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      printf("[hlld_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
      printf("[hlld_sr_mhd]: NaN in hllc density flux\n");
      }*/
		
    return;
  }

  /*--- Step 3b. -----------------------------------------------------------------
   * Construct the initial L & R states
   */
  PaL.S  = Sl;
  PaL.R.d  = Sl*Ul.d  - Fl.d ;
  PaL.R.Mx = Sl*Ul.Mx - Fl.Mx;
  PaL.R.My = Sl*Ul.My - Fl.My;
  PaL.R.Mz = Sl*Ul.Mz - Fl.Mz;
  PaL.R.E  = Sl*Ul.E  - Fl.E ;
  PaL.R.By = Sl*Ul.By - Fl.By;
  PaL.R.Bz = Sl*Ul.Bz - Fl.Bz;
  /*	printf("[hlld_sr_mhd]: Ul = \n");
	printCons1D(&Ul);
	printf("[hlld_sr_mhd]: Fl = \n");
	printCons1D(&Fl);
	printf("[hlld_sr_mhd]: PaL.R = \n");
	printCons1D(&PaL.R);*/

  PaR.S  = Sr;
  PaR.R.d  = (Sr*Ur.d  - Fr.d );
  PaR.R.Mx = (Sr*Ur.Mx - Fr.Mx);
  PaR.R.My = (Sr*Ur.My - Fr.My);
  PaR.R.Mz = (Sr*Ur.Mz - Fr.Mz);
  PaR.R.E  = (Sr*Ur.E  - Fr.E );
  PaR.R.By = (Sr*Ur.By - Fr.By);
  PaR.R.Bz = (Sr*Ur.Bz - Fr.Bz);
  /*	printf("[hlld_sr_mhd]: Ur = \n");
	printCons1D(&Ur);
	printf("[hlld_sr_mhd]: Fr = \n");
	printCons1D(&Fr);
	printf("[hlld_sr_mhd]: PaR.R = \n");
	printCons1D(&PaR.R);*/
	
  /*--- Step 3a. -----------------------------------------------------------------
   * Use the total pressure calculated from the hll average state as an initial
   * guess & make sure it's sensible. If not, revert to hllc flux.
   */
  /* ---- provide an initial guess ---- */
	
  /*switch_to_hllc != getHLLCpress(Wl,Wr,Uhll,Fhll,Bxi,&p0);*/
  /* ---- provide an initial guess ---- */
  switch_to_hllc = 0;
  scrh = MAX(Pl, Pr);
  if (Bxi*Bxi/scrh < 0.01) { /* -- try the B->0 limit -- */
		
    double a,b,c;
    a = Sr - Sl;
    b = PaR.R.E - PaL.R.E + Sr*PaL.R.Mx - Sl*PaR.R.Mx;
    c = PaL.R.Mx*PaR.R.E - PaR.R.Mx*PaL.R.E;
    scrh = b*b - 4.0*a*c;
    scrh = MAX(scrh,0.0);
    p0 = 0.5*(- b + sqrt(scrh))*dS_1;
		
  }else{  /* ----  use HLL average ---- */
    getPtot(Bxi,Whll,&p0);
  }
  p = 1.01*p0;
  /*printf("[hlld_sr_mhd]: p0 = %10.4e\n",p0);*/
	
  /*--- Step 3b. -----------------------------------------------------------------
   * Calculate states based on initial guesses
   */
  switch_to_hllc = 0;
  f0 = Fstar(Bxi, p0, &PaL, &PaR, &Sc);
  /*printf("[hlld_sr_mhd]: Result of initial guess for p0 %10.4e %i\n",f0,PaL.fail);*/
  if (f0 != f0 || PaL.fail) {
    /*printf("[hlld_sr_mhd]: Retrying initial guess\n");*/
    switch_to_hllc != getHLLCpress(Wl,Wr,Uhll,Fhll,Bxi,&p);
    if (p = p0) {
      p0 *= 0.1;
    } else {
      p0 = p;
    }
    switch_to_hllc = 0;
    f0 = Fstar(Bxi, p0, &PaL, &PaR, &Sc);
    if (f0 != f0 || PaL.fail) switch_to_hllc = 1;
    /*printf("[hlld_sr_mhd]: p0 = %10.4e\n",p0);*/
  }
	
  p = 1.025*p0;
  f = Fstar(Bxi, p, &PaL, &PaR, &Sc);
  /*printf("[hlld_sr_mhd]: Result of initial guess for p %10.4e %i\n",f0,PaL.fail);*/
  if (f0 != f0 || PaL.fail) {
    /*printf("[hlld_sr_mhd]: Retrying initial guess\n");*/
    switch_to_hllc != getHLLCpress(Wl,Wr,Uhll,Fhll,Bxi,&p);
    if (p = p0) {
      p0 *= 0.1;
    } else {
      p0 = p;
    }
    switch_to_hllc = 0;
    f0 = Fstar(Bxi, p0, &PaL, &PaR, &Sc);
    if (f0 != f0 || PaL.fail) switch_to_hllc = 1;
		
    p = 1.025*p0;
    f = Fstar(Bxi, p, &PaL, &PaR, &Sc);
    if (f0 != f0 || PaL.fail) switch_to_hllc = 1;
    /*printf("[hlld_sr_mhd]: p0 = %10.4e\n",p0);*/
  }
	
  /*--- Step 3c. -----------------------------------------------------------------
   * Newton-Raphson root finder to find the total pressure & consistent states
   */
  p_incs = 0;
  dPstep = 1.0;
  nr_success = 7*switch_to_hllc;
  if (fabs(f0) < tol && !switch_to_hllc) nr_success = 1;
  /*printf("[hlld_sr_mhd]: Root finder startup: p = %10.4e, f = %10.4e, nr_success = %i\n",p0,f0,nr_success)*/;
  while (nr_success == 0 && p_incs < 100) {
		
    if (fabs(dPstep) < tol || fabs(f) < tol) nr_success = 1;
		
    f  = Fstar(Bxi, p, &PaL, &PaR, &Sc);
    if (f != f) nr_success = 2;
		
    dPstep = (p - p0) / (f - f0)*f;
    p0 = p; f0 = f;
    p -= dPstep;
    dPstep *= 1.0/p;
    if (p != p) nr_success = 3;
    if (p < 0.0) nr_success = 4;	
    /*printf("[hlld_sr_mhd]: Root finder progress: %i, p = %10.4e, f = %10.4e, dP = %10.4e, nr_success = %i\n",p_incs,p,f,dPstep,nr_success);*/
		
    p_incs++;

  }

  if (p < 0.0) switch_to_hllc = 1;
  if (PaL.fail) switch_to_hllc = 1;
  if (nr_success != 1) switch_to_hllc = 1;
	
  if (switch_to_hllc) {
    /*printf("[hlld_sr_mhd]: Unable to find hlld state, error %i\n",nr_success);
    printf("[hlls_sr_mhd]: f0 = %10.4e, p0 = %10.4e, fail = %i\n",f0,p0,PaL.fail);
    printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
    /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;
    return;
  }

  if (fabs(Sc) < 1.0e-5) Sc = 0.0;
  /*printf("[hlld_sr_mhd]: Converged state has p = %10.4e\n",p);
    printf("[hlld_sr_mhd]: Wave speeds Sl = %10.4e, SaL = %10.4e, Sc = %10.4e, SaR = %10.4e, Sr = %10.4e\n",Sl,PaL.Sa,Sc,PaR.Sa,Sr);*/
	
  if (PaL.Sa > 0){      
		
    getAstate (Bxi,p,&PaL,&UaL);
    /*printf("[hlld_sr_mhd]: State 2 \n");*/
    WaL = Cons1D_to_Prim1D (&UaL,&Bxi);
    V2l = SQR(WaL.Vx) + SQR(WaL.Vy) + SQR(WaL.Vz);
    if (WaL.P < 0 || WaL.d < 0 || V2l > 1.0){
      /*printf("[hlld_sr_mhd]: hlld average state has Pl = %10.4e, dl = %10.4e, Vsq = %10.4e\n",WaL.P,WaL.d,V2l);
	printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }

    pFlux->d  = Fl.d  + Sl*(UaL.d  - Ul.d );
    pFlux->Mx = Fl.Mx + Sl*(UaL.Mx - Ul.Mx);
    pFlux->My = Fl.My + Sl*(UaL.My - Ul.My);
    pFlux->Mz = Fl.Mz + Sl*(UaL.Mz - Ul.Mz);
    pFlux->E  = Fl.E  + Sl*(UaL.E  - Ul.E );
    pFlux->By = Fl.By + Sl*(UaL.By - Ul.By);
    pFlux->Bz = Fl.Bz + Sl*(UaL.Bz - Ul.Bz);
		
    if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      /*printf("[hlld_sr_mhd] Sl, Usld, Uld, Fld, %15.8e %15.8e %15.8e %15.8e\n",Sl,UaL.d,Ul.d,Fl.d);
	printf("[hlld_sr_mhd]: NaN in hlld density flux, switching to hllc\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
    return;
		
  } else if (Sc >= 0.0) {
		
    getCstate (Bxi,p,&PaL,&PaR,&UaL,&UaR,&Uc,&SaL,&SaR,&Sc);
    /*		printf("[hlld_sr_mhd]: State 3 \n");*/
    WaL = Cons1D_to_Prim1D (&UaL,&Bxi);
    V2l = SQR(WaL.Vx) + SQR(WaL.Vy) + SQR(WaL.Vz);
    if (WaL.P < 0 || WaL.d < 0 || V2l > 1.0){
      /*printf("[hlld_sr_mhd]: hlld average state has Pl = %10.4e, dl = %10.4e, Vsq = %10.4e\n",WaL.P,WaL.d,V2l);
	printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
		
    Wc = Cons1D_to_Prim1D (&Uc,&Bxi);
    V2c = SQR(Wc.Vx) + SQR(Wc.Vy) + SQR(Wc.Vz);
    if (Wc.P < 0 || Wc.d < 0 || V2c > 1.0){
      /*			printf("[hlld_sr_mhd]: Uc = \n");
				printCons1D(&Uc);*/
      /*			printf("[hlld_sr_mhd]: Sc = %10.4e\n",Sc);
				printf("[hlld_sr_mhd]: hlld average state has Pc = %10.4e, dc = %10.4e, Vsq = %10.4e\n",Wc.P,Wc.d,V2c);*/
      zeroCons1D(&Uc);
      /*			printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      /*			pFlux->d = Fhll.d;
				pFlux->Mx = Fhll.Mx;
				pFlux->My = Fhll.My;
				pFlux->Mz = Fhll.Mz;
				pFlux->E = Fhll.E;
				pFlux->By = Fhll.By;
				pFlux->Bz = Fhll.Bz;
				return;*/
    }

    pFlux->d  = Fl.d  + PaL.Sa*Uc.d  - (PaL.Sa - PaL.S)*UaL.d  - PaL.S*Ul.d;
    pFlux->Mx = Fl.Mx + PaL.Sa*Uc.Mx - (PaL.Sa - PaL.S)*UaL.Mx - PaL.S*Ul.Mx;
    pFlux->My = Fl.My + PaL.Sa*Uc.My - (PaL.Sa - PaL.S)*UaL.My - PaL.S*Ul.My;
    pFlux->Mz = Fl.Mz + PaL.Sa*Uc.Mz - (PaL.Sa - PaL.S)*UaL.Mz - PaL.S*Ul.Mz;
    pFlux->E  = Fl.E  + PaL.Sa*Uc.E  - (PaL.Sa - PaL.S)*UaL.E  - PaL.S*Ul.E;
    pFlux->By = Fl.By + PaL.Sa*Uc.By - (PaL.Sa - PaL.S)*UaL.By - PaL.S*Ul.By;
    pFlux->Bz = Fl.Bz + PaL.Sa*Uc.Bz - (PaL.Sa - PaL.S)*UaL.Bz - PaL.S*Ul.Bz;
		
    if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      /*printf("[hlld_sr_mhd] Sl, SaL, Uald, Uld, Fld, %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,SaL,UaL.d,Ul.d,Fl.d);
	printf("[hlld_sr_mhd]: NaN in hlld density flux, switching to hllc\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
    return;
		
  } else if (Sc < 0.0 && PaR.Sa >= 0.0){
		
    getCstate (Bxi,p,&PaL,&PaR,&UaL,&UaR,&Uc,&SaL,&SaR,&Sc);
    /*		printf("[hlld_sr_mhd]: State 4 \n");*/
		
    WaR = Cons1D_to_Prim1D (&UaR,&Bxi);
    V2r = SQR(WaR.Vx) + SQR(WaR.Vy) + SQR(WaR.Vz);
    if (WaR.P < 0 || WaR.d < 0 || V2r > 1.0){
      /*printf("[hlld_sr_mhd]: hlld average state has Pr = %10.4e, dr = %10.4e, Vsq = %10.4e\n",WaR.P,WaR.d,V2r);*/
      /*printf("[hlld_sr_mhd]: Switching to hllc flux");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
		
    Wc = Cons1D_to_Prim1D (&Uc,&Bxi);
    V2c = SQR(Wc.Vx) + SQR(Wc.Vy) + SQR(Wc.Vz);
    if (Wc.P < 0 || Wc.d < 0 || V2c > 1.0){
      /*			printf("[hlld_sr_mhd]: Uc = \n");
				printCons1D(&Uc);*/
      /*			printf("[hlld_sr_mhd]: Sc = %10.4e\n",Sc);
				printf("[hlld_sr_mhd]: hlld average state has Pc = %10.4e, dc = %10.4e, Vsq = %10.4e\n",Wc.P,Wc.d,V2c);*/
      zeroCons1D(&Uc);
      /*			printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      /*			pFlux->d = Fhll.d;
				pFlux->Mx = Fhll.Mx;
				pFlux->My = Fhll.My;
				pFlux->Mz = Fhll.Mz;
				pFlux->E = Fhll.E;
				pFlux->By = Fhll.By;
				pFlux->Bz = Fhll.Bz;
				return;*/
    }
	
    pFlux->d  = Fr.d  + PaR.Sa*Uc.d  - (PaR.Sa - PaR.S)*UaR.d  - PaR.S*Ur.d;
    pFlux->Mx = Fr.Mx + PaR.Sa*Uc.Mx - (PaR.Sa - PaR.S)*UaR.Mx - PaR.S*Ur.Mx;
    pFlux->My = Fr.My + PaR.Sa*Uc.My - (PaR.Sa - PaR.S)*UaR.My - PaR.S*Ur.My;
    pFlux->Mz = Fr.Mz + PaR.Sa*Uc.Mz - (PaR.Sa - PaR.S)*UaR.Mz - PaR.S*Ur.Mz;
    pFlux->E  = Fr.E  + PaR.Sa*Uc.E  - (PaR.Sa - PaR.S)*UaR.E  - PaR.S*Ur.E;
    pFlux->By = Fr.By + PaR.Sa*Uc.By - (PaR.Sa - PaR.S)*UaR.By - PaR.S*Ur.By;
    pFlux->Bz = Fr.Bz + PaR.Sa*Uc.Bz - (PaR.Sa - PaR.S)*UaR.Bz - PaR.S*Ur.Bz;
			
		
    if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      /*printf("[hlld_sr_mhd] Sl, SaL, Uald, Uld, Fld, %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,SaL,UaL.d,Ul.d,Fl.d);
	printf("[hlld_sr_mhd]: NaN in hlld density flux, switching to hll\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
    return;
		
  } else if (PaR.Sa < 0.0){
		
    getAstate (Bxi,p,&PaR,&UaR);
    /*		printf("[hlld_sr_mhd]: State 5 \n");*/
    WaR = Cons1D_to_Prim1D (&UaR,&Bxi);
    V2r = SQR(WaR.Vx) + SQR(WaR.Vy) + SQR(WaR.Vz);
    if (WaR.P < 0 || WaR.d < 0 || V2r > 1.0){
      /*printf("[hlld_sr_mhd]: hlld average state has Pr = %10.4e, dr = %10.4e, Vsq = %10.4e\n",WaL.P,WaL.d,V2l);*/
      /*printf("[hlld_sr_mhd]: Switching to hllc flux\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
	
    pFlux->d  = Fr.d  + PaR.S*(UaR.d  - Ur.d );
    pFlux->Mx = Fr.Mx + PaR.S*(UaR.Mx - Ur.Mx);
    pFlux->My = Fr.My + PaR.S*(UaR.My - Ur.My);
    pFlux->Mz = Fr.Mz + PaR.S*(UaR.Mz - Ur.Mz);
    pFlux->E  = Fr.E  + PaR.S*(UaR.E  - Ur.E );
    pFlux->By = Fr.By + PaR.S*(UaR.By - Ur.By);
    pFlux->Bz = Fr.Bz + PaR.S*(UaR.Bz - Ur.Bz);
		
    if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      /*printf("[hlld_sr_mhd] Sr, Uard, Urd, Frd, %10.4e %10.4e %10.4e %10.4e\n",Sr,UaR.d,Ur.d,Fr.d);
	printf("[hlld_sr_mhd]: NaN in hlld density flux, switching to hll\n");*/
      /*getHLLCflux (Ul,Ur,Wl,Wr,Sl,Sr,Fl,Fr,Uhll,Whll,Fhll,Bxi,pFlux);*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
      return;
    }
		
    return;
	
  }
	
}

void hlle_fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pFlux)
{
  Cons1DS Fl, Fr;
  Cons1DS Uhll, Fhll;
  Real Pl, Pr;
  Real Sl, Sr;
  Real Sla, Sra;
  Real dS_1;

  /* find min/max wave speeds */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bx,&Sl,&Sr);
  getMaxSignalSpeeds_echo(Wl,Wr,Bx,&Sla,&Sra);
  /*printf("[hlld_sr_mhd]: Exact Wave speeds     %10.4e %10.4e\n",Sl,Sr);
    printf("[hlld_sr_mhd]: Estimated Wave speeds %10.4e %10.4e\n",Sla,Sra);*/
	
  if (Sla != Sla) {
    printf("[hlld_sr_mhd]: NaN in estimated Sl %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra =  1.0;
  }
	
  if (Sra != Sra) {
    printf("[hlld_sr_mhd]: NaN in estimated Sr %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
	
  if (Sla < -1.0) {
    printf("[hlld_sr]: Superluminal estimated Sl %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
  if (Sra > 1.0) {
    printf("[hlld_sr]: Superluminal estimated Sr %10.4e %10.4e\n",Sl,Sr);
    Sla = -1.0;
    Sra = 1.0;
  }
	
  if (Sl != Sl) {
    printf("[hlld_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    printf("[hlld_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    printf("[hlld_sr]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    printf("[hlld_sr]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
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

int getHLLCpress (const Prim1DS Wl,   const Prim1DS Wr,
		  const Cons1DS Uhll, const Cons1DS Fhll,
		  const Real Bxi,     Real *ps) 
{
  Real Bx, Bys, Bzs;
  Real BtFBt, Bt2, FBt2;
  Real a, b, c, vxl, vxr, scrh;
  Real vxs, vys, vzs, gammas_2, vBs, V2l, V2r;
	
  vxl = Wl.Vx;
  vxr = Wr.Vx;
	
  Bx  = Bxi;
  Bys = Uhll.By;
  Bzs = Uhll.Bz;
	
  if (fabs(Bx) < 1.0e-12) {
    a  = Fhll.E;
    b  = - (Fhll.Mx + Uhll.E);
    c  = Uhll.Mx;
  } else {
    BtFBt = Uhll.By*Fhll.By + Uhll.Bz*Fhll.Bz;
    Bt2 = Uhll.By*Uhll.By + Uhll.Bz*Uhll.Bz;
    FBt2 = Fhll.By*Fhll.By + Fhll.Bz*Fhll.Bz;                
		
    a  = Fhll.E - BtFBt;
    b  = Bt2 + FBt2 - (Fhll.Mx + Uhll.E);
    c  = Uhll.Mx - BtFBt;
  }
	
  if (fabs(a) > 1.e-12){       
    /*vxs = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;*/
    scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
    if (scrh != scrh) {
      /*printf("[hllc_sr_mhd]: NaN on scrh: a,b,c, 4*a*c/b/b %10.4e %10.4e %10.4e %10.4e\n",a,b,c,4.0*a*c/(b*b));*/
      return(0);
    }
    vxs  = - 2.0*c/(b*scrh);
  }else{
    vxs = -c/b;
  }
  if (vxs != vxs) {
    /*printf("[hllc_sr_mhd]: HLL fluxes %10.4e %10.4e\n",Fhll.E,BtFBt);
      printf("[hllc_sr_mhd]: NaN in vxs: Bx,a,b,c %10.4e %10.4e %10.4e %10.4e\n",Bx,a,b,c);*/
    return(0);
  }
	
  if (fabs(vxs) > 1.0) {
    /*printf("[hllc_sr_mhd]: vxs > 1: Bx,a,b,c %10.4e %10.4e %10.4e %10.4e\n",Bx,a,b,c);*/
    return(0);
  }
	
  if (fabs(Bx) < 1.0e-12) {
			
    /* -------------------------------
       the value of vy and vz
       is irrelevant in this case  
       ------------------------------- */
			
    *ps  = Fhll.Mx - Fhll.E*vxs;
			
    if (*ps < 0) {
      /*printf("[hllc_sr_mhd]: ps = %10.4e < 0\n",*ps);*/
      return(0);
    }
			
			
  }else{
			
    vys = (Bys*vxs - Fhll.By)/Bx;
    vzs = (Bzs*vxs - Fhll.Bz)/Bx;
			
    gammas_2 = vxs*vxs + vys*vys + vzs*vzs;
    gammas_2 = 1.0 - gammas_2;
    vBs = vxs*Bx + vys*Bys + vzs*Bzs;
			
    *ps = (Bx*vBs - Fhll.E)*vxs + (Bx*Bx*gammas_2) + Fhll.Mx;
    if (*ps < 0) {
      return(0);
      /*printf("[hllc_sr_mhd]: ps = %10.4e < 0\n",*ps);*/
    }
  }

  return(1);
	
}
	
void getHLLCflux (const Cons1DS Ul,   const Cons1DS Ur,
		  const Prim1DS Wl,   const Prim1DS Wr,
		  const Real Sl,      const Real Sr,
		  const Cons1DS Fl,   const Cons1DS Fr,
		  const Cons1DS Uhll, const Prim1DS Whll,
		  const Cons1DS Fhll, const Real Bxi, Cons1DS *pFlux)
{
  Cons1DS Usl,Usr;
  Prim1DS Wsl,Wsr;
  Real Pl, Pr;
  Real dS_1,scrh;
  Real Bx, Bys, Bzs;
  Real BtFBt, Bt2, FBt2;
  Real a, b, c;
  Real ps, vxs, vys, vzs, gammas_2, vBs, V2l, V2r;
  Real vxl, vxr, alpha_l, alpha_r;
  int ps_success,vxs_success,int_success;

  dS_1 = 1.0/(Sr - Sl);
	
  /*--- Step 3. ------------------------------------------------------------------
   * Compute fluxes based on wave speeds (Mignone et al. eqn 26)
   */
  if(Sl >= 0.0){
    /*printf("Flux_L\n");*/
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;
    pFlux->By = Fl.By;
    pFlux->Bz = Fl.Bz;
		
    /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      printf("[hllc_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
	prinf("[hllc_sr_mhd]: NaN in Fl density flux\n");
	}*/
		
    return;
		
  }
  else if(Sr <= 0.0){
    /*printf("Flux_R\n");*/
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;
    pFlux->By = Fr.By;
    pFlux->Bz = Fr.Bz;
		
    /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
      printf("[hllc_sr_mhd] Sl, Sr, dS_1, Fld, Frd, Uld, Urd %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e\n",Sl,Sr,dS_1,Fl.d,Fr.d,Ul.d,Ur.d);
	printf("[hllc_sr_mhd]: NaN in Fr density flux\n");
	}*/
		
    return;
  }
  else {		
    vxs_success = 1;
		
    /* Construct HLLC fluxes */
    /*printf("Flux_HLLC\n");*/
    vxl = Wl.Vx;
    vxr = Wr.Vx;
		
    Bx  = Bxi;/*(Sr*Bxi   - Sl*Bxi) * dS_1;*/
    Bys = Uhll.By;
    Bzs = Uhll.Bz;
		
    if (fabs(Bx) < 1.0e-12) {
      a  = Fhll.E;
      b  = - (Fhll.Mx + Uhll.E);
      c  = Uhll.Mx;
    } else {
      BtFBt = Uhll.By*Fhll.By + Uhll.Bz*Fhll.Bz;
      Bt2 = Uhll.By*Uhll.By + Uhll.Bz*Uhll.Bz;
      FBt2 = Fhll.By*Fhll.By + Fhll.Bz*Fhll.Bz;                
			
      a  = Fhll.E - BtFBt;
      b  = Bt2 + FBt2 - (Fhll.Mx + Uhll.E);
      c  = Uhll.Mx - BtFBt;
    }
		
    /*
      scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
      if (scrh != scrh) {
      hllc_success = 0;
      printf("[hllc_sr_mhd]: NaN on scrh: a,b,c, 4*a*c/b/b %10.4e %10.4e %10.4e %10.4e\n",a,b,c,4.0*a*c/(b*b));
      }
      vxs  = - 2.0*c/(b*scrh);
    */
    if (fabs(a) > 1.e-12){       
      /*vxs = 0.5*(- b - sqrt(b*b - 4.0*a*c))/a;*/
      scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
      if (scrh != scrh) {
	vxs_success = 0;
	/*printf("[hllc_sr_mhd]: NaN on scrh: a,b,c, 4*a*c/b/b %10.4e %10.4e %10.4e %10.4e\n",a,b,c,4.0*a*c/(b*b));*/
      }
      vxs  = - 2.0*c/(b*scrh);
    }else{
      vxs = -c/b;
    }
    if (vxs != vxs) {
      vxs_success = 0;
      /*printf("[hllc_sr_mhd]: HLL fluxes %10.4e %10.4e\n",Fhll.E,BtFBt);
	printf("[hllc_sr_mhd]: NaN in vxs: Bx,a,b,c %10.4e %10.4e %10.4e %10.4e\n",Bx,a,b,c);*/
    }
		
    if (fabs(vxs) > 1.0) {
      vxs_success = 0;
      /*printf("[hllc_sr_mhd]: vxs > 1: Bx,a,b,c %10.4e %10.4e %10.4e %10.4e\n",Bx,a,b,c);*/
    }
		
    ps_success = 0;
    if (vxs_success == 1) {
      ps_success = 1;
      if (fabs(Bx) < 1.0e-12) {
				
	/* -------------------------------
	   the value of vy and vz
	   is irrelevant in this case  
	   ------------------------------- */
				
	ps  = Fhll.Mx - Fhll.E*vxs;
				
	if (ps < 0) {
	  ps_success = 0;
	  /*printf("[hllc_sr_mhd]: ps = %10.4e < 0\n",ps);*/
	}
				
	alpha_l = (Sl - vxl)/(Sl - vxs);
	alpha_r = (Sr - vxr)/(Sr - vxs);
				
	Usl.d = Ul.d*alpha_l;
	Usr.d = Ur.d*alpha_r;
				
	Usl.E = (Sl*Ul.E - Fl.E + ps*vxs)/(Sl - vxs);
	Usr.E = (Sr*Ur.E - Fr.E + ps*vxs)/(Sr - vxs);
				
	Usl.Mx = (Usl.E + ps)*vxs; 
	Usr.Mx = (Usr.E + ps)*vxs;
	Usl.My = Ul.My*alpha_l; 
	Usr.My = Ur.My*alpha_r; 
	Usl.Mz = Ul.Mz*alpha_l; 
	Usr.Mz = Ur.Mz*alpha_r;
				
	Usl.By = Ul.By*alpha_l;
	Usr.By = Ur.By*alpha_r;
	Usl.Bz = Ul.Bz*alpha_l;
	Usr.Bz = Ur.Bz*alpha_r;
				
      }else{
				
	vys = (Bys*vxs - Fhll.By)/Bx;
	vzs = (Bzs*vxs - Fhll.Bz)/Bx;
				
	gammas_2 = vxs*vxs + vys*vys + vzs*vzs;
	gammas_2 = 1.0 - gammas_2;
	vBs = vxs*Bx + vys*Bys + vzs*Bzs;
				
	ps = (Bx*vBs - Fhll.E)*vxs + (Bx*Bx*gammas_2) + Fhll.Mx;
	if (ps < 0) {
	  ps_success = 0;
	  /*printf("[hllc_sr_mhd]: ps = %10.4e < 0\n",ps);*/
	}
				
	alpha_l = (Sl - vxl)/(Sl - vxs);
	alpha_r = (Sr - vxr)/(Sr - vxs);
				
	if (alpha_l != alpha_l) {
	  ps_success = 0;
	  /*printf("[hllc_sr_mhd]: NaN in alpha_l: Sl, vxl, vxs %10.4e %10.4e %10.4e\n",Sl,vxl,vxs);*/
	}
				
	if (alpha_r != alpha_r) {
	  ps_success = 0;
	  /*printf("[hllc_sr_mhd]: NaN in alpha_r: Sr, vxr, vxs %10.4e %10.4e %10.4e\n",Sl,vxl,vxs);*/
	}
				
	Usl.d = Ul.d*alpha_l;
	Usr.d = Ur.d*alpha_r;
				
	Usl.E = (Sl*Ul.E - Fl.E + ps*vxs - vBs*Bx)/(Sl - vxs);
	Usr.E = (Sr*Ur.E - Fr.E + ps*vxs - vBs*Bx)/(Sr - vxs);
				
	Usl.Mx = (Usl.E + ps)*vxs - vBs*Bx; 
	Usr.Mx = (Usr.E + ps)*vxs - vBs*Bx;
	Usl.My = (Sl*Ul.My - Fl.My - Bx*(Bys*gammas_2 + vBs*vys))/(Sl - vxs); 
	Usr.My = (Sr*Ur.My - Fr.My - Bx*(Bys*gammas_2 + vBs*vys))/(Sr - vxs);  
	Usl.Mz = (Sl*Ul.Mz - Fl.Mz - Bx*(Bzs*gammas_2 + vBs*vzs))/(Sl - vxs); 
	Usr.Mz = (Sr*Ur.Mz - Fr.Mz - Bx*(Bzs*gammas_2 + vBs*vzs))/(Sr - vxs);
				
	Usl.By = Usr.By = Bys;
	Usl.Bz = Usr.Bz = Bzs;
      }
    }
		
    /*  ----  Compute HLLC flux  ----  */
    int_success = 0;
    if (vxs_success && ps_success){
      int_success = 1;
      Wsl = Cons1D_to_Prim1D (&Usl,&Bx);
      Wsr = Cons1D_to_Prim1D (&Usr,&Bx);
      V2l = SQR(Wsl.Vx) + SQR(Wsl.Vy) + SQR(Wsl.Vz);
      V2r = SQR(Wsr.Vx) + SQR(Wsr.Vy) + SQR(Wsr.Vz);
      if (Wsl.P < 0.0) int_success = 0;
      if (Wsl.d < 0.0) int_success = 0;
      if (Wsr.P < 0.0) int_success = 0;
      if (Wsr.d < 0.0) int_success = 0;
      if (V2l   > 1.0) int_success = 0;
      if (V2r   > 1.0) int_success = 0;
      /*if (int_success == 0){
	printf("[hllc_sr_mhd]: Wave speeds            Sl = %10.4e, Sr = %10.4e, Vxs = %10.4e\n",Sl,Sr,vxs);
	printf("[hllc_sr_mhd]: Intermediate state has Pl = %10.4e, dl = %10.4e, V2l = %10.4e\n",Wsl.P,Wsl.d,V2l);
	printf("[hllc_sr_mhd]: Intermediate state has Pr = %10.4e, dr = %10.4e, V2l = %10.4e\n",Wsr.P,Wsr.d,V2r);
	}*/
    }
		
    if (int_success){
      if (vxs > 0.0) {
	pFlux->d  = Fl.d  + Sl*(Usl.d  - Ul.d );
	pFlux->Mx = Fl.Mx + Sl*(Usl.Mx - Ul.Mx);
	pFlux->My = Fl.My + Sl*(Usl.My - Ul.My);
	pFlux->Mz = Fl.Mz + Sl*(Usl.Mz - Ul.Mz);
	pFlux->E  = Fl.E  + Sl*(Usl.E  - Ul.E );
	pFlux->By = Fl.By + Sl*(Usl.By - Ul.By);
	pFlux->Bz = Fl.Bz + Sl*(Usl.Bz - Ul.Bz);
	/*pFlux->Mx+= Pl;*/
				
	if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
	  /*printf("[hllc_sr_mhd] Sl, Usld, Uld, Fld, %10.4e %10.4e %10.4e %10.4e\n",Sl,Usl.d,Ul.d,Fl.d);
	  printf("[hllc_sr_mhd]: NaN in hllc density flux\n");
	  printf("[hllc_sr_mhd]: Reverting to hll flux\n");*/
					
	  pFlux->d = Fhll.d;
	  pFlux->Mx = Fhll.Mx;
	  pFlux->My = Fhll.My;
	  pFlux->Mz = Fhll.Mz;
	  pFlux->E = Fhll.E;
	  pFlux->By = Fhll.By;
	  pFlux->Bz = Fhll.Bz;
					
	  /*if (pFlux->d != pFlux->d) {
	    printf("[hllc_sr_mhd] Sr, Usrd, Urd, Frd %10.4e %10.4e %10.4e %10.4e \n",Sl,Usr.d,Ur.d,Fr.d);
	    printf("[hllc_sr_mhd]: NaN in hll density flux\n");
	    }*/
					
	  return;
	}
				
	return;
				
      }else {
	pFlux->d  = Fr.d  + Sr*(Usr.d  - Ur.d );
	pFlux->Mx = Fr.Mx + Sr*(Usr.Mx - Ur.Mx);
	pFlux->My = Fr.My + Sr*(Usr.My - Ur.My);
	pFlux->Mz = Fr.Mz + Sr*(Usr.Mz - Ur.Mz);
	pFlux->E  = Fr.E  + Sr*(Usr.E  - Ur.E );
	pFlux->By = Fr.By + Sr*(Usr.By - Ur.By);
	pFlux->Bz = Fr.Bz + Sr*(Usr.Bz - Ur.Bz);
	/*pFlux->Mx+= Pr;*/
				
	if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
	  /*printf("[hllc_sr_mhd] Sr, Usrd, Urd, Frd %10.4e %10.4e %10.4e %10.4e \n",Sl,Usr.d,Ur.d,Fr.d);
	  printf("[hllc_sr_mhd]: NaN in hllc density flux\n");
	  printf("[hllc_sr_mhd]: Reverting to hll flux\n");*/
					
	  pFlux->d = Fhll.d;
	  pFlux->Mx = Fhll.Mx;
	  pFlux->My = Fhll.My;
	  pFlux->Mz = Fhll.Mz;
	  pFlux->E = Fhll.E;
	  pFlux->By = Fhll.By;
	  pFlux->Bz = Fhll.Bz;
					
	  /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
	    printf("[hllc_sr_mhd] Sr, Usrd, Urd, Frd %10.4e %10.4e %10.4e %10.4e \n",Sl,Usr.d,Ur.d,Fr.d);
	    printf("[hllc_sr_mhd]: NaN in hll density flux\n");
	    }*/
					
	  return;
	}
				
	return;
      }
			
    } else {
      /*printf("[hllc_sr_mhd]: Reverting to hll flux\n");*/
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
			
      /*if (pFlux->d != pFlux->d || pFlux->Mx != pFlux->Mx) {
	printf("[hllc_sr_mhd] Sr, Usrd, Urd, Frd %10.4e %10.4e %10.4e %10.4e \n",Sl,Usr.d,Ur.d,Fr.d);
	ath_error("[hllc_sr_mhd]: NaN in hllc density flux\n");
	}*/
			
      return;
			
    }
		
  }
	
}

void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p)
{
  Real wtg2, pt, g, g2, g_2, h, gmmr, theta;
  Real bx, by, bz, vB, b2, Bmag2;
	
  /* calcUlate enthalpy */
	
  theta = W.P/W.d;
  gmmr = Gamma / Gamma_1;
	
  h = 1.0 + gmmr*theta;
	
  /* calcUlate gamma */
  g   = U.d/W.d;
  g2  = SQR(g);
  g_2 = 1.0/g2;
	
  pt = W.P;
  wtg2 = W.d*h*g2;
	
  vB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  Bmag2 = SQR(Bx) + SQR(W.By) + SQR(W.Bz);
	
  bx = g*(  Bx*g_2 + vB*W.Vx);
  by = g*(W.By*g_2 + vB*W.Vy);
  bz = g*(W.Bz*g_2 + vB*W.Vz);
	
  b2 = Bmag2*g_2 + vB*vB;
	
  pt += 0.5*b2;
  wtg2 += b2*g2;
	
  flux->d  = U.d*W.Vx;
  flux->Mx = wtg2*W.Vx*W.Vx + pt;
  flux->My = wtg2*W.Vy*W.Vx;
  flux->Mz = wtg2*W.Vz*W.Vx;
  flux->E  = U.Mx;
	
  flux->Mx -= bx*bx;
  flux->My -= by*bx;
  flux->Mz -= bz*bx;
  flux->By = W.Vx*W.By - Bx*W.Vy;
  flux->Bz = W.Vx*W.Bz - Bx*W.Vz;
	
  *p = pt;
}

void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high)
{
	
  Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
  Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
  Real al,ar;
	
  getVChar_pluto(Wl,Bx,&lml,&lpl);
  /*printf("[hlle_sr]: Left characteristics %10.4e %10.4e\n",lml,lpl);*/
  getVChar_pluto(Wr,Bx,&lmr,&lpr);
  /*printf("[hlle_sr]: Right characteristics %10.4e %10.4e\n",lmr,lpr);*/
	
  *low =  MIN(lml, lmr);
  *high = MAX(lpl, lpr);
}

/* *********************************************************************** */
Real Fstar (const Real Bx, const Real p, Riemann_State *PaL, Riemann_State *PaR, Real *Sc)
/*
 *
 *
 *
 ************************************************************************* */
{
  int    success = 1;
  Real dK, Bxc, Byc, Bzc;
  Real sBx, fun, scrh;
  Real vxcL, KLBc, yspL;
  Real vxcR, KRBc, yspR;
	
	
  sBx = Bx > 0.0 ? 1.0:-1.0;
	
  success *= getRiemannState (-1, Bx, p, PaL);
  success *= getRiemannState ( 1, Bx, p, PaR);
	
  /* -- compute B from average state -- */
  if (success){
    dK  = PaR->Kx - PaL->Kx + 1.e-12;
			
    Bxc = dK*Bx;
    Byc =   (PaR->By*(PaR->Kx - PaR->vx) + Bx*PaR->vy) 
      - (PaL->By*(PaL->Kx - PaL->vx) + Bx*PaL->vy); 
    Bzc =   (PaR->Bz*(PaR->Kx - PaR->vx) + Bx*PaR->vz)
      - (PaL->Bz*(PaL->Kx - PaL->vx) + Bx*PaL->vz);
		
    KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
    KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;
		
    yspL = (1.0 - PaL->K2) / (dK*PaL->eta - KLBc);
    yspR = (1.0 - PaR->K2) / (dK*PaR->eta - KRBc);
		
    vxcL = PaL->Kx - Bx*(1.0 - PaL->K2)/(PaL->eta - KLBc/dK);
    vxcR = PaR->Kx - Bx*(1.0 - PaR->K2)/(PaR->eta - KRBc/dK);
    *Sc  = 0.5*(vxcL + vxcR);
	
    PaL->Sa = PaL->Kx;
    PaR->Sa = PaR->Kx;
    fun     = vxcL - vxcR;
    /*dK*(1.0 - Bx*(yspR - yspL));*/
    if (fun != fun){
      printf("[hlld_sr_mhd]: NaN in Fstar\n");
      printf("[hlld_sr_mhd]: dK = %10.4e, Bx = %10.4e, yspL = %10.4e yspR = %10.4e\n",dK,Bxc,yspL,yspR);
    }
		
    /* -- check if state makes physically sense -- */
    success  = vxcL >= PaL->Kx;
    success *= vxcR <= PaR->Kx;
    /*	if (!success){
	printf("[hlld_sr_mhd]: vxcL = %10.4e > PaL->Kx = %10.4e\n",vxcL,PaL->Kx);
	printf("[hlld_sr_mhd]: vxcR = %10.4e < PaR->Kx = %10.4e\n",vxcR,PaR->Kx);
	}*/
		
    success *= PaL->vx >= PaL->S;
    success *= PaR->vx <= PaR->S;
    /*	if (!success){
	printf("[hlld_sr_mhd]: SL   = %10.4e < PaL->Vx = %10.4e\n",PaL->S,PaL->vx);
	printf("[hlld_sr_mhd]: SR   = %10.4e > PaR->Vx = %10.4e\n",PaR->S,PaR->vx);
	}*/
	
    success *= (PaR->w >= p);
    success *= (PaL->w >= p);
    /*	if (!success){
	printf("[hlld_sr_mhd]: p    = %15.7e < PaL->w  = %15.7e\n",p,PaL->w);
	printf("[hlld_sr_mhd]: p    = %15.7e < PaR->w  = %15.7e\n",p,PaR->w);
	}*/


  }	
	
  PaL->fail = !success;
	
  return (fun);
}

/* *********************************************************************** */
int getRiemannState (const int side, const Real Bx, const Real p, Riemann_State *Pv)
/*
 *
 *  
 *  Return 1 if succesfull, 0 if w < 0 is encountered.
 *  side = -1 : means left
 *  side =  1 : means right
 *
 ************************************************************************* */
{
  Cons1DS R;
  Real A, C, G, X, Q, s;
  Real vx, vy, vz, scrh, S;
  Real vsq,gamma,gamma2;
  Real vB,E,sigma,b0;
	
  S = Pv->S;
  R = Pv->R;
	
  A = R.Mx + p*(1.0 - S*S) - S*R.E;
  C = R.By*R.My + R.Bz*R.Mz;
  G = SQR(R.By) + SQR(R.Bz);
  X = Bx*(A*S*Bx + C) - (A + G)*(S*p + R.E);
  Q = -A - G + Bx*Bx*(1.0 - S*S);
	
  vx = Bx*(A*Bx + C*S) - (R.Mx + p)*(G + A);
  vy = Q*R.My + R.By*(C + Bx*(S*R.Mx - R.E));
  vz = Q*R.Mz + R.Bz*(C + Bx*(S*R.Mx - R.E));
	
  scrh = 1.0/X;
  Pv->vx = vx*scrh; 
  Pv->vy = vy*scrh;
  Pv->vz = vz*scrh;
	
  Pv->Bx = Bx; 
  Pv->By = (R.By - Bx*Pv->vy)/(S - Pv->vx);
  Pv->Bz = (R.Bz - Bx*Pv->vz)/(S - Pv->vx);

  /*
    vsq = SQR(Pv->vx) + SQR(Pv->vy) + SQR(Pv->vz);
    if (vsq > 1.0) {
    printf("[hlld_sr_mhd]: |V|^2 = %10.4e > 1\n",vsq);
    return(0);
    }
  */
  /*	gamma2 = (1.0 - vsq);
	
	scrh = 1.0/(S - Pv->vx);
	vB  = Pv->vx*Bx + Pv->vy*Pv->By + Pv->vz*Pv->Bz;
	E = (R.E + p*Pv->vx - vB*Bx)*scrh;
	Pv->w = (E + p)*gamma2 + SQR(vB);*/
	
	
  scrh = Pv->vx*R.Mx + Pv->vy*R.My + Pv->vz*R.Mz;
  scrh = (R.E - scrh) / (S - Pv->vx);
  Pv->w = p + scrh;
  /*	if (Pv->w < p){
	printf("[hlld_sr_mhd]: Got w = %10.4e < p = %10.4e\n",Pv->w,p);
	printf("[hlld_sr_mhd]: Mx = %10.4e, My = %10.4e, Mz = %10.4e\n",R.Mx,R.My,R.Mz);
	printf("[hlld_sr_mhd]: Vx = %10.4e, Vy = %10.4e, Vz = %10.4e\n",Pv->vx,Pv->vy,Pv->vz);
	printf("[hlld_sr_mhd]: R.E = %10.4e, p = %10.4e, scrh = %10.e, S = %10.4e\n",R.E,p,scrh,S);
	}*/
	
  /* -- failure -- */
  if (Pv->w < 0.0){ 
    /*		printf("[hlld_sr_mhd]: Got a negative total entropy, w = %10.4e\n",Pv->w);
		printf("[hlld_sr_mhd]: Mx = %10.4e, My = %10.4e, Mz = %10.4e\n",R.Mx,R.My,R.Mz);
		printf("[hlld_sr_mhd]: Vx = %10.4e, Vy = %10.4e, Vz = %10.4e\n",vx,vy,vz);
		printf("[hlld_sr_mhd]: R.E = %10.4e, p = %10.4e, scrh = %10.e, S = %10.4e\n",R.E,p,scrh,S);
    */		return(0);  
  }
	
  s  = Bx > 0.0 ? 1.0:-1.0;
  if (side < 0) s *= -1.0;
	
  Pv->eta = s*sqrt(Pv->w);

  vsq = SQR(Pv->vx) + SQR(Pv->vy) + SQR(Pv->vz);
  if (vsq > 1.0) {
    /*printf("[hlld_sr_mhd]: |V|^2 = %10.4e > 1\n",vsq);*/
    return(0);
  }
  gamma = 1.0/sqrt(1.0 - vsq);
  vB  = Pv->vx*Bx + Pv->vy*Pv->By + Pv->vz*Pv->Bz;
  b0 = gamma*vB;
  sigma = Pv->eta*gamma + b0;
  scrh = 1.0/(gamma*sigma);
	
  Pv->Kx = Pv->vx + Pv->Bx*scrh;
  Pv->Ky = Pv->vy + Pv->By*scrh;
  Pv->Kz = Pv->vz + Pv->Bz*scrh;

  /*
    scrh = 1.0/(S*p +  R.E + Bx*Pv->eta);
    Pv->Kx = scrh*(R.Mx + p + Bx*Pv->eta);
    Pv->Ky = scrh*(R.My     + R.By*Pv->eta);
    Pv->Kz = scrh*(R.Mz     + R.Bz*Pv->eta);
  */
	
  Pv->K2 = SQR(Pv->Kx) + SQR(Pv->Ky) + SQR(Pv->Kz);
  return(1); /* -- success -- */
}

/* ************************************************************* */
void getAstate (const Real Bx    , const Real p,
		Riemann_State *Pa, Cons1DS *U  )
/*
 *        Compute states aL and aR (behind fast waves)
 *
 *************************************************************** */
{
  Cons1DS R;
  Real vB, S;
  Real scrh;
	
  S  = Pa->S;
  R  = Pa->R;
	
  scrh = 1.0/(S - Pa->vx);
	
  U->d =  R.d*scrh;

  U->By = (R.By - Bx*Pa->vy)*scrh;
  U->Bz = (R.Bz - Bx*Pa->vz)*scrh;
	
  vB  = Pa->vx*Bx + Pa->vy*U->By + Pa->vz*U->Bz;
  U->E = (R.E + p*Pa->vx - vB*Bx)*scrh;
	
  U->Mx = (U->E + p)*Pa->vx - vB*Bx  ;
  U->My = (U->E + p)*Pa->vy - vB*(U->By);
  U->Mz = (U->E + p)*Pa->vz - vB*(U->Bz);
}

void getCstate (const Real Bx     , const Real p,
		Riemann_State *PaL, Riemann_State *PaR,
		Cons1DS *UaL,       Cons1DS *UaR,       Cons1DS *Uc,
		Real *SaL,          Real *SaR,          Real *Sc)
/*
 *         Compute state  cL and cR across contact mode
 *
 *************************************************************** */
{
  Cons1DS R,Ua;
  Real vB, S, Sa, Vxa;
  Real Bxc, Byc, Bzc, dK;
  Real vxcL, vycL, vzcL, KLBc;
  Real vxcR, vycR, vzcR, KRBc;
  Real vxc,vyc,vzc;
  Real scrhL, scrhR;
	
  getAstate (Bx,p,PaL,UaL);
  getAstate (Bx,p,PaR,UaR);
	
  *SaL = PaL->Sa;
  *SaR = PaR->Sa;
	
  /* -- compute B from average state -- */
  dK  = PaR->Kx - PaL->Kx + 1.e-12;
	
  Byc =   (PaR->By*(PaR->Kx - PaR->vx) + Bx*PaR->vy) 
    - (PaL->By*(PaL->Kx - PaL->vx) + Bx*PaL->vy); 
  Byc*= 1.0/dK;
  Bzc =   (PaR->Bz*(PaR->Kx - PaR->vx) + Bx*PaR->vz)
    - (PaL->Bz*(PaL->Kx - PaL->vx) + Bx*PaL->vz);
  Bzc*= 1.0/dK;
	
  Uc->By = Byc;
  Uc->Bz = Bzc;
	
  KLBc = PaL->Kx*Bx + PaL->Ky*Byc + PaL->Kz*Bzc;
  KRBc = PaR->Kx*Bx + PaR->Ky*Byc + PaR->Kz*Bzc;

  scrhL = (1.0 - PaL->K2)/(PaL->eta - KLBc);
  scrhR = (1.0 - PaR->K2)/(PaR->eta - KRBc);

  vxcL = PaL->Kx - Bx*scrhL;
  vxcR = PaR->Kx - Bx*scrhR;
	
  vycL = PaL->Ky - Uc->By*scrhL;
  vycR = PaR->Ky - Uc->By*scrhR;
	
  vzcL = PaL->Kz - Uc->Bz*scrhL;
  vzcR = PaR->Kz - Uc->Bz*scrhR;
	
  vxc = 0.5 * (vxcL + vxcR);
  vyc = 0.5 * (vycL + vycR);
  vzc = 0.5 * (vzcL + vzcR);
	
  if (vxc > 0.0) {
    Ua = *UaL;
    Sa = *SaL;
    Vxa = PaL->vx;
  } else {
    Ua  = *UaR;
    Sa  = *SaR;
    Vxa = PaR->vx;
  }

  vB  = vxc*Bx + vyc*Uc->By + vzc*Uc->Bz;
	
  Uc->d =  Ua.d*(Sa - Vxa)/(Sa - vxc);
  if (Uc->d == 0.0 && vxc != 0.0) {
    printf("[hlld_sr_mhd]: Contact state has d = 0, Bx = %10.4e\n",Bx);
    printf("[hlld_sr_mhd]: Ua = \n");
    printCons1D(&Ua);
    printf("[hlld_sr_mhd]: Sa = %10.4e, Vxa = %10.4e, vxc = %10.4e\n",Sa,Vxa,vxc);
  }
	
  Uc->E = (Sa*Ua.E - Ua.Mx + p*vxc - vB*Bx)/(Sa - vxc);
	
  Uc->Mx = (Uc->E + p)*vxc - vB*Bx  ;
  Uc->My = (Uc->E + p)*vyc - vB*(Uc->By);
  Uc->Mz = (Uc->E + p)*vzc - vB*(Uc->Bz);
	
  *Sc      = 0.5*(vxcL + vxcR);
}

/* ************************************************************* */
void getPtot (const Real Bx, const Prim1DS W, Real *pt)
/*
 *
 * Compute total pressure
 *
 *************************************************************** */
{
  Real vel2, Bmag2, vB;
	
  vel2 = SQR(W.Vx) + SQR(W.Vy) + SQR(W.Vz);
  Bmag2 = SQR(Bx) + SQR(W.By) + SQR(W.Bz);
  vB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
	
  *pt = W.P + 0.5*(Bmag2*(1.0 - vel2) + vB*vB);
}

void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
{
  Real rhoh,vsq,bsq;
  Real cssq,vasq,asq;
  Real Vx2,gamma,gamma2;
  Real Bx2,Bsq,vDotB,vDotBsq,b0,bx;
  Real w_1,a0,a1,a2,a3,a4,Q;
  Real scrh,scrh1,scrh2,eps2;
  Real a2_w,one_m_eps2,lambda[5];
  int iflag;
	
  rhoh = W.d + (Gamma/Gamma_1) * (W.P);
	
  Vx2 = SQR(W.Vx);
  vsq = Vx2 + SQR(W.Vy) + SQR(W.Vz);
  if (vsq > 1.0){
    printf("[getVChar]: |v|= %f > 1\n",vsq);
		
    *lm = -1.0;
    *lp = 1.0;
		
    return;
		
  }
  gamma = 1.0 / sqrt(1 - vsq);    
  gamma2 = SQR(gamma);
	
  Bx2 = SQR(Bx);
  Bsq = Bx2 + SQR(W.By) + SQR(W.Bz);
  vDotB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  vDotBsq = SQR(vDotB);
  b0 = gamma * vDotB;
  bx = Bx/gamma2 + W.Vx*vDotB;
  bsq = Bsq / gamma2 + SQR(vDotB);
	
  cssq = (Gamma * W.P) / (rhoh);
  vasq = bsq / (rhoh + bsq);
	
  if (bsq < 0.0) bsq = 0.0;
  if (cssq < 0.0) cssq = 0.0;
  if (cssq > 1.0) cssq = 1.0;
  if (vasq > 1.0) bsq = rhoh + bsq;
  /*printf("[hlle_sr]: Sound, Alfven & Magneto sonic speeds %10.4e %10.4e %10.4e\n",cssq,vasq,asq);*/
	
  if (vsq < 1.0e-12) {
    w_1  = 1.0/(rhoh + bsq);   
    eps2 = cssq + bsq*w_1*(1.0 - cssq);
    a0   = cssq*Bx*Bx*w_1;
    a1   = - a0 - eps2;
    scrh = a1*a1 - 4.0*a0;
    if (scrh < 0.0) scrh = 0.0;
		
    scrh = sqrt(0.5*(-a1 + sqrt(scrh)));
    *lp =  scrh;
    *lm = -scrh;
    /*printf("[getVChar]: Zero velocity wave speeds %10.4e %10.4e %10.4e %10.4e %10.4e\n",*lm,*lp,a0,a1,scrh);*/
    return;
  }
	
  /*	
	vB2 = vB*vB;
	u02 = 1.0/(1.0 - u02);
	b2  = b2/u02 + vB2;
	u0  = sqrt(u02);
  */
  w_1 = 1.0/(rhoh + bsq);   
	
  if (Bx < 1.0e-14) {
		
    eps2  = cssq + bsq*w_1*(1.0 - cssq);
		
    scrh1 = (1.0 - eps2)*gamma2;
    scrh2 = cssq*vDotBsq*w_1 - eps2;
		
    a2  = scrh1 - scrh2;
    a1  = -2.0*W.Vx*scrh1;
    a0  = Vx2*scrh1 + scrh2;
		
    *lp = 0.5*(-a1 + sqrt(a1*a1 - 4.0*a2*a0))/a2;
    *lm = 0.5*(-a1 - sqrt(a1*a1 - 4.0*a2*a0))/a2;
		
    return;
  }
	
  scrh1 = bx;  /* -- this is bx/u0 -- */
  scrh2 = scrh1*scrh1;  
	
  a2_w       = cssq*w_1;
  eps2       = (cssq*rhoh + bsq)*w_1;
  one_m_eps2 = gamma2*rhoh*(1.0 - cssq)*w_1;
	
  /* ---------------------------------------
     Define coefficients for the quartic  
     --------------------------------------- */
	
  scrh = 2.0*(a2_w*vDotB*scrh1 - eps2*W.Vx);
  a4 = one_m_eps2 - a2_w*vDotBsq + eps2;
  a3 = - 4.0*W.Vx*one_m_eps2 + scrh;
  a2 =   6.0*Vx2*one_m_eps2 + a2_w*(vDotBsq - scrh2) + eps2*(Vx2 - 1.0);
  a1 = - 4.0*W.Vx*Vx2*one_m_eps2 - scrh;
  a0 = Vx2*Vx2*one_m_eps2 + a2_w*scrh2 - eps2*Vx2;
	
  if (a4 < 1.e-12){
    /*printPrim1D(W);*/
    printf("[MAX_CH_SPEED]: Can not divide by a4 in MAX_CH_SPEED\n");
		
    *lm = -1.0;
    *lp = 1.0;
		
    return;
  }
	
  scrh = 1.0/a4;
	
  a3 *= scrh;
  a2 *= scrh;
  a1 *= scrh;
  a0 *= scrh;
  iflag = QUARTIC(a3, a2, a1, a0, lambda);
	
  if (iflag){
    printf ("Can not find max speed:\n");
    /*SHOW(uprim,i);*/
    printf("QUARTIC: f(x) = %12.6e + x*(%12.6e + x*(%12.6e ",
	   a0*a4, a1*a4, a2*a4);
    printf("+ x*(%12.6e + x*%12.6e)))\n", a3*a4, a4);
    printf("[MAX_CH_SPEED]: Failed to find wave speeds");
		
    *lm = -1.0;
    *lp = 1.0;
		
    return;
  }

  /*
   *lp = MAX(lambda[3], lambda[2]);
   *lp = MAX(*lp, lambda[1]);
   *lp = MAX(*lp, lambda[0]);
   *lp = MIN(*lp, 1.0);
	
   *lm = MIN(lambda[3], lambda[2]);
   *lm = MIN(*lm, lambda[1]);
   *lm = MIN(*lm, lambda[0]);
   *lm = MAX(*lm, -1.0);
   */
  /*printf("[hlld_sr_mhd]: wave speeds %10.4e %10.4e %10.4e %10.4e %10.4e\n",lambda[0],lambda[1],lambda[2],lambda[3],lambda[4]);*/
	
  *lp = MIN(1.0,MAX(lambda[3], lambda[2]));
  *lp = MIN(1.0,MAX(*lp, lambda[1]));
  *lp = MIN(1.0,MAX(*lp, lambda[0]));

  *lm = MAX(-1.0,MIN(lambda[3], lambda[2]));
  *lm = MAX(-1.0,MIN(*lm, lambda[1]));
  *lm = MAX(-1.0,MIN(*lm, lambda[0]));
	
  return;
	
}

/* ******************************************** */
int QUARTIC (Real b, Real c, Real d, 
             Real e, Real z[])
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
  Real b2, f, g, h;
  Real a2, a1, a0, u[4];
  Real p, q, r, s;
  static Real three_256 = 3.0/256.0;
  static Real one_64 = 1.0/64.0;
	
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
int CUBIC(Real b, Real c, Real d, Real z[])
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
  Real b2, g2;
  Real f, g, h;
  Real i, i2, j, k, m, n, p;
  static Real one_3 = 1.0/3.0, one_27=1.0/27.0;
	
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

/* **************************************************************** */
void getSoundSpeed2  (const Prim1DS W, Real *cs2, Real *h)
/*
 *
 *    Define the square of the sound speed for adiabatic EOS
 *
 ****************************************************************** */
{
  Real  theta, Gamma_1;
	
  Gamma_1 = Gamma/(Gamma - 1.0);
  theta = W.P/W.d;
  *h = 1.0 + Gamma_1*theta;
  *cs2 = Gamma*theta/(*h);
}

void getMaxSignalSpeeds_echo (const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high)
{
	
  Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
  Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
  Real al,ar;
	
  getVChar_echo(Wl,Bx,&lml,&lpl);
  /*printf("[hlle_sr]: Left characteristics %10.4e %10.4e\n",lml,lpl);*/
  getVChar_echo(Wr,Bx,&lmr,&lpr);
  /*printf("[hlle_sr]: Right characteristics %10.4e %10.4e\n",lmr,lpr);*/
	
  *low =  MIN(lml, lmr);
  *high = MAX(lpl, lpr);
}

void getVChar_echo(const Prim1DS W, const Real Bx, Real* lm, Real* lp)
{
  Real rhoh,vsq,bsq;
  Real cssq,vasq,asq;
  Real gamma,gamma2;
  Real Bsq,vDotB,b0,bx;
  Real tmp1,tmp2,tmp3,tmp4,tmp5;
  Real vm,vp;
	
  rhoh = W.d + (Gamma/Gamma_1) * (W.P);
	
  vsq = SQR(W.Vx) + SQR(W.Vy) + SQR(W.Vz);
  gamma = 1.0 / sqrt(1 - vsq);    
  gamma2 = SQR(gamma);
	
  Bsq = SQR(Bx) + SQR(W.By) + SQR(W.Bz);
  vDotB = W.Vx*Bx + W.Vy*W.By + W.Vz*W.Bz;
  b0 = gamma * vDotB;
  bx = Bx/gamma2 + W.Vx*vDotB;
  bsq = Bsq / gamma2 + SQR(vDotB);
	
  cssq = (Gamma * W.P) / (rhoh);
  vasq = bsq / (rhoh + bsq);
  asq = cssq + vasq - (cssq*vasq);
	
  if (cssq < 0.0) cssq = 0.0;
  if (vasq > 0.0) vasq = 0.0;
  if (asq < 0.0) asq = 0.0;
  if (cssq > 1.0) cssq = 1.0;
  if (vasq > 1.0) vasq = 1.0;
  if (asq > 1.0) asq = 1.0;
  /*printf("[hlle_sr]: Sound, Alfven & Magneto sonic speeds %10.4e %10.4e %10.4e\n",cssq,vasq,asq);*/
	
  tmp1 = (1.0 - asq);
  tmp2 = (1.0 - vsq);
  tmp3 = (1.0 - vsq*asq);
  tmp4 = SQR(W.Vx);
  tmp5 = 1.0 / tmp3;
	
  vm = tmp1*W.Vx - sqrt(asq*tmp2*(tmp3 - tmp1*tmp4));
  vp = tmp1*W.Vx + sqrt(asq*tmp2*(tmp3 - tmp1*tmp4));
  vm *=tmp5;
  vp *=tmp5;
	
  if (vp > vm) {
    *lm = vm;
    *lp = vp;
  } else {
    *lm = vp;
    *lp = vm;
  }
}



#endif /* MHD */
#undef MAX_ITER

#endif
#endif
