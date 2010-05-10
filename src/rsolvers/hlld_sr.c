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

#define MAX_ITER 100

typedef struct CONS_STATE{
  Real DN,EN;
  Real M1,M2,M3;
  Real B1,B2,B3;
} CONS_STATE;

typedef struct RIEMANN_STATE{
  int fail;
  Real vx, vy, vz;
  Real Bx, By, Bz;
  Real Kx, Ky, Kz, K2;
  Real w, eta, p, rho;
  CONS_STATE u, R;
  Real S, Sa, sw;
} Riemann_State;

void zeroCons1D(Cons1DS *U);
void printCons1D(const Cons1DS *U);
void printPrim1D(const Prim1DS *W);

void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p);

void get_HLLC_flux (const Real Sl,  const Real Sr,
		    const Prim1DS Wl, const Prim1DS Wr,
		    const Cons1DS Ul, const Cons1DS Ur,
		    const Cons1DS Fl, const Cons1DS Fr,
		    const Cons1DS Uhll,const Cons1DS Fhll,
		    const Real Bxi, Cons1DS *pFlux);

Real Fstar (Riemann_State *PaL, Riemann_State *PaR, 
	    Real *Sc, Real p, const Real Bx);
int GET_RIEMANN_STATE (Riemann_State *Pv, Real p, int side, const Real Bx);
void GET_ASTATE (Riemann_State *Pa,  Real p, const Real Bx);
void GET_CSTATE (Riemann_State *PaL, Riemann_State *PaR, Real p,
		 CONS_STATE *Uc, const Real Bx);
int checkCurrentSheet (const Prim1DS Wl,const Prim1DS Wr,
			const Cons1DS Ul,const Cons1DS Ur,const Real Bxi);
void getPtot (const Real Bx, const Prim1DS W, Real *pt);
void getMaxSignalSpeeds_pluto(const Prim1DS Wl, const Prim1DS Wr,
                              const Real Bx, Real* low, Real* high);
void getMaxSignalSpeeds_echo(const Prim1DS Wl, const Prim1DS Wr,
                             const Real Bx, Real* low, Real* high);
void getVChar_echo(const Prim1DS W, const Real Bx, Real* lml, Real* lmr);
void getVChar_pluto(const Prim1DS W, const Real Bx, Real* lml, Real* lmr);

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

/*----------------------------------------------------------------------------*/
/* fluxes
 *   Input Arguments:
 *     Ul,Ur = L/R-states of CONSERVED variables at cell interface 
 *     Wl,Wr = L/R-states of PRIMITIVE variables at cell interface 
 *   Output Arguments:
 *     pFlux = pointer to fluxes of CONSERVED variables at cell interface 
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bxi, Cons1DS *pFlux)
{
  CONS_STATE Uhll,Fhll,Uc;
  Cons1DS Fl,Fr,Usl,Usr,Utmp,Ftmp;
  Prim1DS Wsl,Wsr,Whll,Wavg;
  Real Sl, Sr, Pl, Pr;
  Real Sla, Sra;
  Real dS_1,scrh,Sc;
  Real Bx,p0,pguess,f0,p,f,dp;
  Riemann_State PaL,PaR;
  int switch_to_hll,wave_speed_fail,k,ISTEP;
  int current_sheet_fail;

  wave_speed_fail = 0;
  switch_to_hll = 0;
	
/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bxi,&Sl,&Sr);
/*   Wavg.d = 0.5*(Wl.d + Wr.d); */
/*   Wavg.P = 0.5*(Wl.P + Wr.P); */
/*   Wavg.Vx = 0.5*(Wl.Vx + Wr.Vx); */
/*   Wavg.Vy = 0.5*(Wl.Vy + Wr.Vy); */
/*   Wavg.Vz = 0.5*(Wl.Vz + Wr.Vz); */
/*   Wavg.By = 0.5*(Wl.By + Wr.By); */
/*   Wavg.Bz = 0.5*(Wl.Bz + Wr.Bz); */
/*   getVChar_pluto(Wavg, Bxi, &Sl, &Sr); */
	
  if (Sl != Sl) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    wave_speed_fail = 1;
    printf("[hllc_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

/*--- Step 1a. -----------------------------------------------------------------
 * If PLUTO wavespeeds are bad, fall back to the estimate used in ECHO
 */
  if (wave_speed_fail){
    getMaxSignalSpeeds_echo (Wl,Wr,Bxi,&Sla,&Sra);
	
    if (Sla != Sla) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra =  1.0;
    }
	
    if (Sra != Sra) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
	
    if (Sla < -1.0) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
    if (Sra > 1.0) {
      switch_to_hll = 1;
      printf("[hllc_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }

    Sl = Sla;
    Sr = Sra;

  }

  /* compute L/R fluxes */
  flux_LR(Ul,Wl,&Fl,Bxi,&Pl);
  flux_LR(Ur,Wr,&Fr,Bxi,&Pr);
	
/*--- Step 2. ------------------------------------------------------------------
 * Construct HLL fluxes & average state
 */
  dS_1 = 1.0/(Sr - Sl);

  Uhll.DN = (Sr*Ur.d  - Sl*Ul.d  + Fl.d  - Fr.d ) * dS_1;
  Uhll.M1 = (Sr*Ur.Mx - Sl*Ul.Mx + Fl.Mx - Fr.Mx) * dS_1;
  Uhll.M2 = (Sr*Ur.My - Sl*Ul.My + Fl.My - Fr.My) * dS_1;
  Uhll.M3 = (Sr*Ur.Mz - Sl*Ul.Mz + Fl.Mz - Fr.Mz) * dS_1;
  Uhll.EN = (Sr*Ur.E  - Sl*Ul.E  + Fl.E  - Fr.E ) * dS_1;
  Uhll.B1 = (Sr*   Bxi- Sl*   Bxi               ) * dS_1;
  Uhll.B2 = (Sr*Ur.By - Sl*Ul.By + Fl.By - Fr.By) * dS_1;
  Uhll.B3 = (Sr*Ur.Bz - Sl*Ul.Bz + Fl.Bz - Fr.Bz) * dS_1;
		
  Fhll.DN = (Sr*Fl.d  - Sl*Fr.d  + Sl*Sr*(Ur.d  - Ul.d )) * dS_1;
  Fhll.M1 = (Sr*Fl.Mx - Sl*Fr.Mx + Sl*Sr*(Ur.Mx - Ul.Mx)) * dS_1;
  Fhll.M2 = (Sr*Fl.My - Sl*Fr.My + Sl*Sr*(Ur.My - Ul.My)) * dS_1;
  Fhll.M3 = (Sr*Fl.Mz - Sl*Fr.Mz + Sl*Sr*(Ur.Mz - Ul.Mz)) * dS_1;
  Fhll.EN = (Sr*Fl.E  - Sl*Fr.E  + Sl*Sr*(Ur.E  - Ul.E )) * dS_1;
  Fhll.B1 = (                      Sl*Sr*(   Bxi-   Bxi)) * dS_1;
  Fhll.B2 = (Sr*Fl.By - Sl*Fr.By + Sl*Sr*(Ur.By - Ul.By)) * dS_1;
  Fhll.B3 = (Sr*Fl.Bz - Sl*Fr.Bz + Sl*Sr*(Ur.Bz - Ul.Bz)) * dS_1;


  if (switch_to_hll) {
    pFlux->d = Fhll.DN;
    pFlux->Mx = Fhll.M1;
    pFlux->My = Fhll.M2;
    pFlux->Mz = Fhll.M3;
    pFlux->E = Fhll.EN;
    pFlux->By = Fhll.B2;
    pFlux->Bz = Fhll.B3;

    if (pFlux->d != pFlux->d) {
      printf("[hlld_sr_mhd]: NaN in hll density flux\n");
      printf("[hlld_sr_mhd]: Sl = %10.4e, Sr = %10.4e\n",Sl,Sr);
      printf("[hlld_sr_mhd]: dS_1 = %10.4e\n",dS_1);
      printf("[hlld_sr_mhd]: Fld = %10.4e, Frd = %10.4e\n",Fl.d,Fr.d);
      printf("[hlld_sr_mhd]: Uld = %10.4e, Urd = %10.4e\n",Ul.d,Ur.d);
    }
		
    return;
  }

/*--- Step 3. ------------------------------------------------------------------
 * Compute fluxes based on wave speeds (Mignone et al. eqn 26)
 */
  if(Sl >= 0.0){
    pFlux->d  = Fl.d;
    pFlux->Mx = Fl.Mx;
    pFlux->My = Fl.My;
    pFlux->Mz = Fl.Mz;
    pFlux->E  = Fl.E;
    pFlux->By = Fl.By;
    pFlux->Bz = Fl.Bz;

    if (pFlux->d != pFlux->d) {
      printf("[hlld_sr_mhd]: NaN in hllc density flux for Sl >= 0.0\n");
      printf("[hlld_sr_mhd]: Sl = %10.4e, Sr = %10.4e\n",Sl,Sr);
      printf("[hlld_sr_mhd]: dS_1 = %10.4e\n",dS_1);
      printf("[hlld_sr_mhd]: Fld = %10.4e, Frd = %10.4e\n",Fl.d,Fr.d);
      printf("[hlld_sr_mhd]: Uld = %10.4e, Urd = %10.4e\n",Ul.d,Ur.d);
    }
		
    return;
   
  }
  else if(Sr <= 0.0){
    pFlux->d  = Fr.d;
    pFlux->Mx = Fr.Mx;
    pFlux->My = Fr.My;
    pFlux->Mz = Fr.Mz;
    pFlux->E  = Fr.E;
    pFlux->By = Fr.By;
    pFlux->Bz = Fr.Bz;

    if (pFlux->d != pFlux->d) {
      printf("[hlld_sr_mhd]: NaN in hllc density flux for Sr <= 0.0\n");
      printf("[hlld_sr_mhd]: Sl = %10.4e, Sr = %10.4e\n",Sl,Sr);
      printf("[hlld_sr_mhd]: dS_1 = %10.4e\n",dS_1);
      printf("[hlld_sr_mhd]: Fld = %10.4e, Frd = %10.4e\n",Fl.d,Fr.d);
      printf("[hlld_sr_mhd]: Uld = %10.4e, Urd = %10.4e\n",Ul.d,Ur.d);
    }
		
    return;
  }
  else {

    PaL.S = Sl;
    PaR.S = Sr;
    Bx    = Uhll.B1;

    PaL.R.DN = Sl*Ul.d  - Fl.d;
    PaL.R.EN = Sl*Ul.E  - Fl.E;
    PaL.R.M1 = Sl*Ul.Mx - Fl.Mx;
    PaL.R.M2 = Sl*Ul.My - Fl.My;
    PaL.R.M3 = Sl*Ul.Mz - Fl.Mz;
    PaL.R.B1 = Sl*  Bxi        ;
    PaL.R.B2 = Sl*Ul.By - Fl.By;
    PaL.R.B3 = Sl*Ul.Bz - Fl.Bz;

    PaR.R.DN = Sr*Ur.d  - Fr.d;
    PaR.R.EN = Sr*Ur.E  - Fr.E;
    PaR.R.M1 = Sr*Ur.Mx - Fr.Mx;
    PaR.R.M2 = Sr*Ur.My - Fr.My;
    PaR.R.M3 = Sr*Ur.Mz - Fr.Mz;
    PaR.R.B1 = Sr*  Bxi        ;
    PaR.R.B2 = Sr*Ur.By - Fr.By;
    PaR.R.B3 = Sr*Ur.Bz - Fr.Bz;

    /* ---- provide an initial guess ---- */

    Utmp.d = Uhll.DN;
    Utmp.Mx = Uhll.M1;
    Utmp.My = Uhll.M2;
    Utmp.Mz = Uhll.M3;
    Utmp.E = Uhll.EN;
    Utmp.By = Uhll.B2;
    Utmp.Bz = Uhll.B3;

    Ftmp.d = Fhll.DN;
    Ftmp.Mx = Fhll.M1;
    Ftmp.My = Fhll.M2;
    Ftmp.Mz = Fhll.M3;
    Ftmp.E = Fhll.EN;
    Ftmp.By = Fhll.B2;
    Ftmp.Bz = Fhll.B3;

    /*Whll = check_Prim1D(&Utmp, &Bxi);
      getPtot(Bxi,Whll,&p0);*/

    scrh = MAX(Wl.P, Wr.P);
    if (Bx*Bx/scrh < 0.01) { /* -- try the B->0 limit -- */
	
      Real a,b,c;
      a = Sr - Sl;
      b = PaR.R.EN - PaL.R.EN + Sr*PaL.R.M1 - Sl*PaR.R.M1;
      c = PaL.R.M1*PaR.R.EN - PaR.R.M1*PaL.R.EN;
      scrh = b*b - 4.0*a*c;
      scrh = MAX(scrh,0.0);
      p0 = 0.5*(- b + sqrt(scrh))*dS_1;
      
    } else {  /* ----  use HLL average ---- */

	/*Utmp.d = Uhll.DN;
      Utmp.Mx = Uhll.M1;
      Utmp.My = Uhll.M2;
      Utmp.Mz = Uhll.M3;
      Utmp.E = Uhll.EN;
      Utmp.By = Uhll.B2;
      Utmp.Bz = Uhll.B3;*/
  
      Whll = check_Prim1D(&Utmp, &Bxi);                       
      getPtot(Bxi,Whll,&p0);
    }
      
    pguess = p0;
  /* ---- check if guess makes sense ---- */

    /*switch_to_hll = 0;
    f0 = Fstar(&PaL, &PaR, &Sc, p0, Bx);
    if (f0 != f0 || PaL.fail){

      Whll.d  = (Sr*Wr.d  - Sl*Wl.d ) * dS_1;
      Whll.Vx = (Sr*Wr.Vx - Sl*Wl.Vx) * dS_1;
      Whll.Vy = (Sr*Wr.Vy - Sl*Wl.Vy) * dS_1;
      Whll.Vz = (Sr*Wr.Vz - Sl*Wl.Vz) * dS_1;
      Whll.P  = (Sr*Wr.P  - Sl*Wl.P ) * dS_1;
      Whll.By = (Sr*Wr.By - Sl*Wl.By) * dS_1;
      Whll.Bz = (Sr*Wr.Bz - Sl*Wl.Bz) * dS_1;
      getPtot(Bxi,Whll,&p0);*/

      pguess = p0;
      switch_to_hll = 0;
      f0 = Fstar(&PaL, &PaR, &Sc, p0, Bx);
      if (f0 != f0 || PaL.fail) switch_to_hll = 1;
      /*}*/

    /* check whether state contains a \beta < 0.1 current sheet */
    /*current_sheet_fail = checkCurrentSheet (Wl,Wr,Ul,Ur,Bxi);
    if (current_sheet_fail){
      switch_to_hll = 1;
      printf("[hlld_sr]: Current sheet switch = %i\n",current_sheet_fail);
      }*/

    /* ---- Root finder ---- */

    k = 0;
    if (fabs(f0) > 1.e-12 && !switch_to_hll){
      p  = 1.025*p0; f  = f0;
      for (k = 1; k < MAX_ITER; k++){
	
	/*printf ("k = %i, Sc = %10.4e, p = %10.4e, f = %10.4e\n",k,Sc,p,f);*/
	
	f  = Fstar(&PaL, &PaR, &Sc, p, Bx);
	if ( f != f  || PaL.fail || (k > 7) || 
	     (fabs(f) > fabs(f0) && k > 4)) {
	  switch_to_hll = 1;
	  break;
	}
	/*if ( f != f ) {
	  switch_to_hll = 1;
	  break;
	  }*/
	
	dp = (p - p0)/(f - f0)*f;
	
	p0 = p; f0 = f;
	p -= dp;
	if (p < 0.0) p = 1.e-6;
	if (fabs(dp) < 1.e-5*p || fabs(f) < 1.e-6) break;
      }
    }else p = p0;
    
    if (PaL.fail) switch_to_hll = 1;
    
    if (switch_to_hll) {

      /*printf("[hlld_sr_mhd]: Failiure in non-linear root finder\n");
      printf("[hlld_sr_mhd]: k = %i, p = %10.4e, f = %10.4e, fail = %i\n",
	     k,p,f,PaL.fail);
	     printf("[hlld_sr_mhd]: Switching to HLL fluxes\n");*/

      /*get_HLLC_flux (Sl,Sr,Wl,Wr,Ul,Ur,Fl,Fr,Utmp,Ftmp,Bxi,pFlux);*/

      *pFlux = Ftmp;
      
      return;
    }

    if (PaL.Sa >= -1.e-6){

      GET_ASTATE (&PaL, p, Bx);

      pFlux->d  = Fl.d  + Sl*(PaL.u.DN  - Ul.d );
      pFlux->Mx = Fl.Mx + Sl*(PaL.u.M1  - Ul.Mx);
      pFlux->My = Fl.My + Sl*(PaL.u.M2  - Ul.My);
      pFlux->Mz = Fl.Mz + Sl*(PaL.u.M3  - Ul.Mz);
      pFlux->E  = Fl.E  + Sl*(PaL.u.EN  - Ul.E );
      pFlux->By = Fl.By + Sl*(PaL.u.B2  - Ul.By);
      pFlux->Bz = Fl.Bz + Sl*(PaL.u.B3  - Ul.Bz);

      if (pFlux->d  != pFlux->d || pFlux->E   != pFlux->E ||
          pFlux->Mx != pFlux->Mx || pFlux->My != pFlux->My ||
          pFlux->Mz != pFlux->Mz || pFlux->By != pFlux->By ||
          pFlux->Bz != pFlux->Bz) {
	printf("[hllc_sr]: NaN in hlld flux, PaL.Sa > 0.0\n");
	printf("[hlld_sr] F.d  = %10.4e, F.E  = %10.4e\n",pFlux->d,pFlux->E);
	printf("[hlld_sr] F.Mx = %10.4e, F.My = %10.4e\n",pFlux->Mx,pFlux->My);
        printf("[hlld_sr] F.Mz = %10.4e, F.By = %10.4e\n",pFlux->Mz,pFlux->By);
        printf("[hlld_sr] F.Bz = %10.4e\n",pFlux->Bz);
	printf("[hlld_sr]: Switching to hll fluxes\n");

	pFlux->d = Ftmp.d;
	pFlux->Mx = Ftmp.Mx;
	pFlux->My = Ftmp.My;
	pFlux->Mz = Ftmp.Mz;
	pFlux->E = Ftmp.E;
	pFlux->By = Ftmp.By;
	pFlux->Bz = Ftmp.Bz;

      }

      return;
      /*      for (nv = NFLX; nv--;   ) {
	state->flux[i][nv] = fL[nv] + SL[i]*(PaL.u[nv] - uL[nv]);
	}*/

    }else if (PaR.Sa <= 1.e-6){

      GET_ASTATE (&PaR, p, Bx);

      pFlux->d  = Fr.d  + Sr*(PaR.u.DN  - Ur.d );
      pFlux->Mx = Fr.Mx + Sr*(PaR.u.M1  - Ur.Mx);
      pFlux->My = Fr.My + Sr*(PaR.u.M2  - Ur.My);
      pFlux->Mz = Fr.Mz + Sr*(PaR.u.M3  - Ur.Mz);
      pFlux->E  = Fr.E  + Sr*(PaR.u.EN  - Ur.E );
      pFlux->By = Fr.By + Sr*(PaR.u.B2  - Ur.By);
      pFlux->Bz = Fr.Bz + Sr*(PaR.u.B3  - Ur.Bz);

      if (pFlux->d  != pFlux->d || pFlux->E   != pFlux->E ||
          pFlux->Mx != pFlux->Mx || pFlux->My != pFlux->My ||
          pFlux->Mz != pFlux->Mz || pFlux->By != pFlux->By ||
          pFlux->Bz != pFlux->Bz) {
	printf("[hllc_sr]: NaN in hlld flux, PaL.Sa > 0.0\n");
	printf("[hlld_sr] F.d  = %10.4e, F.E  = %10.4e\n",pFlux->d,pFlux->E);
	printf("[hlld_sr] F.Mx = %10.4e, F.My = %10.4e\n",pFlux->Mx,pFlux->My);
        printf("[hlld_sr] F.Mz = %10.4e, F.By = %10.4e\n",pFlux->Mz,pFlux->By);
        printf("[hlld_sr] F.Bz = %10.4e\n",pFlux->Bz);
	printf("[hlld_sr]: Switching to hll fluxes\n");

	pFlux->d = Ftmp.d;
	pFlux->Mx = Ftmp.Mx;
	pFlux->My = Ftmp.My;
	pFlux->Mz = Ftmp.Mz;
	pFlux->E = Ftmp.E;
	pFlux->By = Ftmp.By;
	pFlux->Bz = Ftmp.Bz;

      }

      return;

      /*for (nv = NFLX; nv--;   ) {
	state->flux[i][nv] = fR[nv] + SR[i]*(PaR.u[nv] - uR[nv]);
	}*/

    }else{

      GET_CSTATE (&PaL, &PaR, p, &Uc, Bx);

      if (Sc > 0.0){

	pFlux->d  = Fl.d  + Sl*(PaL.u.DN  - Ul.d ) + PaL.Sa*(Uc.DN - PaL.u.DN);
	pFlux->Mx = Fl.Mx + Sl*(PaL.u.M1  - Ul.Mx) + PaL.Sa*(Uc.M1 - PaL.u.M1);
	pFlux->My = Fl.My + Sl*(PaL.u.M2  - Ul.My) + PaL.Sa*(Uc.M2 - PaL.u.M2);
	pFlux->Mz = Fl.Mz + Sl*(PaL.u.M3  - Ul.Mz) + PaL.Sa*(Uc.M3 - PaL.u.M3);
	pFlux->E  = Fl.E  + Sl*(PaL.u.EN  - Ul.E ) + PaL.Sa*(Uc.EN - PaL.u.EN);
	pFlux->By = Fl.By + Sl*(PaL.u.B2  - Ul.By) + PaL.Sa*(Uc.B2 - PaL.u.B2);
	pFlux->Bz = Fl.Bz + Sl*(PaL.u.B3  - Ul.Bz) + PaL.Sa*(Uc.B3 - PaL.u.B3);

	if (pFlux->d  != pFlux->d || pFlux->E   != pFlux->E ||
	    pFlux->Mx != pFlux->Mx || pFlux->My != pFlux->My ||
	    pFlux->Mz != pFlux->Mz || pFlux->By != pFlux->By ||
	    pFlux->Bz != pFlux->Bz) {
	  printf("[hllc_sr]: NaN in hlld flux, PaL.Sa > 0.0\n");
	  printf("[hlld_sr] F.d  = %10.4e, F.E  = %10.4e\n",pFlux->d,pFlux->E);
	  printf("[hlld_sr] F.Mx = %10.4e, F.My = %10.4e\n",pFlux->Mx,pFlux->My);
	  printf("[hlld_sr] F.Mz = %10.4e, F.By = %10.4e\n",pFlux->Mz,pFlux->By);
	  printf("[hlld_sr] F.Bz = %10.4e\n",pFlux->Bz);
	  printf("[hlld_sr]: Switching to hll fluxes\n");

	  pFlux->d = Ftmp.d;
	  pFlux->Mx = Ftmp.Mx;
	  pFlux->My = Ftmp.My;
	  pFlux->Mz = Ftmp.Mz;
	  pFlux->E = Ftmp.E;
	  pFlux->By = Ftmp.By;
	  pFlux->Bz = Ftmp.Bz;

	}

	return;
	/*for (nv = NFLX; nv--;   ) {
	  state->flux[i][nv] = fL[nv] + SL[i]*(PaL.u[nv] - uL[nv]) 
                                        + PaL.Sa*(Uc[nv] - PaL.u[nv]);
					}*/


        }else{

	pFlux->d  = Fr.d  + Sr*(PaR.u.DN  - Ur.d ) + PaR.Sa*(Uc.DN - PaR.u.DN);
	pFlux->Mx = Fr.Mx + Sr*(PaR.u.M1  - Ur.Mx) + PaR.Sa*(Uc.M1 - PaR.u.M1);
	pFlux->My = Fr.My + Sr*(PaR.u.M2  - Ur.My) + PaR.Sa*(Uc.M2 - PaR.u.M2);
	pFlux->Mz = Fr.Mz + Sr*(PaR.u.M3  - Ur.Mz) + PaR.Sa*(Uc.M3 - PaR.u.M3);
	pFlux->E  = Fr.E  + Sr*(PaR.u.EN  - Ur.E ) + PaR.Sa*(Uc.EN - PaR.u.EN);
	pFlux->By = Fr.By + Sr*(PaR.u.B2  - Ur.By) + PaR.Sa*(Uc.B2 - PaR.u.B2);
	pFlux->Bz = Fr.Bz + Sr*(PaR.u.B3  - Ur.Bz) + PaR.Sa*(Uc.B3 - PaR.u.B3);

	if (pFlux->d  != pFlux->d || pFlux->E   != pFlux->E ||
	    pFlux->Mx != pFlux->Mx || pFlux->My != pFlux->My ||
	    pFlux->Mz != pFlux->Mz || pFlux->By != pFlux->By ||
	    pFlux->Bz != pFlux->Bz) {
	  printf("[hllc_sr]: NaN in hlld flux, PaL.Sa > 0.0\n");
	  printf("[hlld_sr] F.d  = %10.4e, F.E  = %10.4e\n",pFlux->d,pFlux->E);
	  printf("[hlld_sr] F.Mx = %10.4e, F.My = %10.4e\n",pFlux->Mx,pFlux->My);
	  printf("[hlld_sr] F.Mz = %10.4e, F.By = %10.4e\n",pFlux->Mz,pFlux->By);
	  printf("[hlld_sr] F.Bz = %10.4e\n",pFlux->Bz);
	  printf("[hlld_sr]: Switching to hll fluxes\n");

	  pFlux->d = Ftmp.d;
	  pFlux->Mx = Ftmp.Mx;
	  pFlux->My = Ftmp.My;
	  pFlux->Mz = Ftmp.Mz;
	  pFlux->E = Ftmp.E;
	  pFlux->By = Ftmp.By;
	  pFlux->Bz = Ftmp.Bz;

	}

	return;
	/*for (nv = NFLX; nv--;   ) {
            state->flux[i][nv] = fR[nv] + SR[i]*(PaR.u[nv] - uR[nv]) 
                                        + PaR.Sa*(Uc[nv] - PaR.u[nv]);
					}*/

      }  
    }
  }
}

void hlle_fluxes(const Cons1DS Ul, const Cons1DS Ur,
		 const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Cons1DS *pFlux)
{
  Cons1DS Fl, Fr;
  Cons1DS Uhll, Fhll;
  Prim1DS Whll;
  Real Pl, Pr;
  Real Sl, Sr;
  Real Sla, Sra;
  Real dS_1, V2l;
  int wave_speed_fail;

  wave_speed_fail = 0;
	
/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bx,&Sl,&Sr);
	
  if (Sl != Sl) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

/*--- Step 1a. -----------------------------------------------------------------
 * If PLUTO wavespeeds are bad, fall back to the estimate used in ECHO
 */
  if (wave_speed_fail){
    getMaxSignalSpeeds_echo (Wl,Wr,Bx,&Sla,&Sra);
	
    if (Sla != Sla) {
      printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra =  1.0;
    }
	
    if (Sra != Sra) {
      printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
	
    if (Sla < -1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
    if (Sra > 1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }

    Sl = Sla;
    Sr = Sra;

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
	
    /*Whll = check_Prim1D(&Uhll, &Bx);
    V2l = SQR(Whll.Vx) + SQR(Whll.Vy) + SQR(Whll.Vz);
    if (Whll.P < 0 || Whll.d < 0 || V2l > 1.0){
      printf("[hlle_sr_mhd]: Unphysical hll average state\n");
      printf("[hlle_sr_mhd]: Phll = %10.4e, dhll = %10.4e\n",Whll.P,Whll.d);
      printf("[hlle_sr_mhd]: V^2 hll = %10.4e\n",V2l);
      printf("[hlle_sr_mhd]: Reverting to LF fluxes\n");
      Sl = -1.0;
      Sr =  1.0;
      }*/

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

    if (pFlux->d != pFlux->d) {
      printf("[hlle_sr_mhd]: NaN in hll density flux\n");
      printf("[hlle_sr_mhd]: Sl = %10.4e, Sr = %10.4e\n",Sl,Sr);
      printf("[hlle_sr_mhd]: dS_1 = %10.4e\n",dS_1);
      printf("[hlle_sr_mhd]: Fld = %10.4e, Frd = %10.4e\n",Fl.d,Fr.d);
      printf("[hlle_sr_mhd]: Uld = %10.4e, Urd = %10.4e\n",Ul.d,Ur.d);
    }

    return;
  }
}

void get_HLLC_flux (const Real Sl,      const Real Sr,
		    const Prim1DS Wl,   const Prim1DS Wr,
		    const Cons1DS Ul,   const Cons1DS Ur,
		    const Cons1DS Fl,   const Cons1DS Fr,
		    const Cons1DS Uhll, const Cons1DS Fhll,
		    const Real Bxi,           Cons1DS *pFlux)
{
  Cons1DS Usl,Usr;
  Real Bx, Bys, Bzs;
  Real BtFBt, Bt2, FBt2;
  Real a, b, c, scrh;
  Real ps, vxs, vys, vzs, gammas_2, vBs, V2l, V2r;
  Real vxl, vxr, alpha_l, alpha_r;
  int switch_to_hll;

  switch_to_hll = 0;
		
  /* Construct HLLC fluxes */
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
    scrh = 1.0 + sqrt(1.0 - 4.0*a*c/(b*b));
    if (scrh != scrh) {
      switch_to_hll = 1;
    }
    vxs  = - 2.0*c/(b*scrh);
  } else {
    vxs = -c/b;
  }
  if ((vxs != vxs || vxs > 1.0) && (switch_to_hll == 0)) {
    switch_to_hll = 1;
  }
		
  if (switch_to_hll) {
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;
		
    return;
  }
		
  if (fabs(Bx) < 1.0e-12) {
			
    /* -------------------------------
       the value of vy and vz
       is irrelevant in this case  
       ------------------------------- */
			
    ps  = Fhll.Mx - Fhll.E*vxs;

    if (ps < 0) {
      switch_to_hll = 1;
    } else {			
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
    }
		
  } else {

    vys = (Bys*vxs - Fhll.By)/Bx;
    vzs = (Bzs*vxs - Fhll.Bz)/Bx;

    gammas_2 = vxs*vxs + vys*vys + vzs*vzs;
    gammas_2 = 1.0 - gammas_2;
    vBs = vxs*Bx + vys*Bys + vzs*Bzs;
    
    ps = (Bx*vBs - Fhll.E)*vxs + (Bx*Bx*gammas_2) + Fhll.Mx;
    if (ps < 0) {
      switch_to_hll = 1;
      /*printf("[hllc_sr_mhd]: ps = %10.4e < 0\n",ps);
	printf("[hllc_sr_mhd]: Switching to hll fluxes\n");*/
    } else {
			
      alpha_l = (Sl - vxl)/(Sl - vxs);
      alpha_r = (Sr - vxr)/(Sr - vxs);
			
      if (alpha_l != alpha_l) {
	switch_to_hll = 1;
	/*printf("[hllc_sr_mhd]: NaN in alpha_l:\n");
	  printf("[hllc_sr_mhd]: Sl = %10.4e\n",Sl);
	  printf("[hllc_sr_mhd]: vxl = %10.4e, vxs = %10.4e\n",vxl,vxs);
	  printf("[hllc_sr_mhd]: Switching to hll fluxes\n");*/
      }
			
      if (alpha_r != alpha_r) {
	switch_to_hll = 1;
	/*printf("[hllc_sr_mhd]: NaN in alpha_r:\n");
	  printf("[hllc_sr_mhd]: Sr = %10.4e\n",Sr);
	  printf("[hllc_sr_mhd]: vxr = %10.4e, vxs = %10.4e\n",vxr,vxs);
	  printf("[hllc_sr_mhd]: Switching to hll fluxes\n");*/
      }

      if (switch_to_hll == 0){
			
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
  }

  if (switch_to_hll) {
    pFlux->d = Fhll.d;
    pFlux->Mx = Fhll.Mx;
    pFlux->My = Fhll.My;
    pFlux->Mz = Fhll.Mz;
    pFlux->E = Fhll.E;
    pFlux->By = Fhll.By;
    pFlux->Bz = Fhll.Bz;
		
    return;
  }
		
    /*  ----  Compute HLLC flux  ----  */

  if (vxs > 0.0) {
    pFlux->d  = Fl.d  + Sl*(Usl.d  - Ul.d );
    pFlux->Mx = Fl.Mx + Sl*(Usl.Mx - Ul.Mx);
    pFlux->My = Fl.My + Sl*(Usl.My - Ul.My);
    pFlux->Mz = Fl.Mz + Sl*(Usl.Mz - Ul.Mz);
    pFlux->E  = Fl.E  + Sl*(Usl.E  - Ul.E );
    pFlux->By = Fl.By + Sl*(Usl.By - Ul.By);
    pFlux->Bz = Fl.Bz + Sl*(Usl.Bz - Ul.Bz);
			
    if (pFlux->d != pFlux->d) {
      pFlux->d = Fhll.d;
      pFlux->Mx = Fhll.Mx;
      pFlux->My = Fhll.My;
      pFlux->Mz = Fhll.Mz;
      pFlux->E = Fhll.E;
      pFlux->By = Fhll.By;
      pFlux->Bz = Fhll.Bz;
    }
			
    return;
			
  } else {
    pFlux->d  = Fr.d  + Sr*(Usr.d  - Ur.d );
    pFlux->Mx = Fr.Mx + Sr*(Usr.Mx - Ur.Mx);
    pFlux->My = Fr.My + Sr*(Usr.My - Ur.My);
    pFlux->Mz = Fr.Mz + Sr*(Usr.Mz - Ur.Mz);
    pFlux->E  = Fr.E  + Sr*(Usr.E  - Ur.E );
    pFlux->By = Fr.By + Sr*(Usr.By - Ur.By);
    pFlux->Bz = Fr.Bz + Sr*(Usr.Bz - Ur.Bz);
			
    if (pFlux->d != pFlux->d) {
	pFlux->d = Fhll.d;
	pFlux->Mx = Fhll.Mx;
	pFlux->My = Fhll.My;
	pFlux->Mz = Fhll.Mz;
	pFlux->E = Fhll.E;
	pFlux->By = Fhll.By;
	pFlux->Bz = Fhll.Bz;
    }
			
    return;

  }
}

/* *********************************************************************** */
Real Fstar (Riemann_State *PaL, Riemann_State *PaR, 
	    Real *Sc, Real p, const Real Bx)
/*
 *
 *
 *
 ************************************************************************* */
{
  int    success = 1;
  Real dK, Bxc, Byc, Bzc;
  Real sBx, fun;
  Real vxcL, KLBc;
  Real vxcR, KRBc;

  sBx = Bx > 0.0 ? 1.0:-1.0;

  success *= GET_RIEMANN_STATE (PaL, p, -1, Bx);
  success *= GET_RIEMANN_STATE (PaR, p,  1, Bx);

/* -- compute B from average state -- */

  dK  = PaR->Kx - PaL->Kx + 1.e-12;

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
  *Sc     = 0.5*(vxcL + vxcR);
  /*fun     = vxcL - vxcR;*/


fun = dK*(1.0 - Bx*(  (1.0 - PaR->K2)/(PaR->sw*dK - KRBc)
                    - (1.0 - PaL->K2)/(PaL->sw*dK - KLBc)) );


  /* -- check if state makes physically sense -- */

  success  = (vxcL - PaL->Kx)  > -1.e-6;
  success *= (PaR->Kx  - vxcR) > -1.e-6;

  success *= (PaL->S - PaL->vx) < 0.0;
  success *= (PaR->S - PaR->vx) > 0.0;

  success *= (PaR->w - p) > 0.0;
  success *= (PaL->w - p) > 0.0;
  success *= (PaL->Sa - PaL->S)  > -1.e-6;
  success *= (PaR->S  - PaR->Sa) > -1.e-6;

  PaL->fail = !success;

/*
scrh  = (1.0 - PaR->K2)*(PaL->sw*dK - KLBc);
scrh -= (1.0 - PaL->K2)*(PaR->sw*dK - KRBc);

PaL->fun1 = (PaR->sw*dK - KRBc)*(PaL->sw*dK - KLBc) - Bx*scrh; 
PaL->fun2 = scrh;
PaL->denL = (PaL->sw*dK - KLBc);
PaL->denR = (PaR->sw*dK - KRBc);
*/
  return (fun);
}

/* *********************************************************************** */
int GET_RIEMANN_STATE (Riemann_State *Pv, Real p, int side, const Real Bx)
/*
 *
 *  
 *  Return 1 if succesfull, 0 if w < 0 is encountered.
 *  side = -1 : means left
 *  side =  1 : means right
 *
 ************************************************************************* */
{
  double A, C, G, X, s;
  double vx, vy, vz, scrh;

  A = Pv->R.M1 + p*(1.0 - Pv->S*Pv->S) - Pv->S*Pv->R.EN;
  C = Pv->R.B2*Pv->R.M2 + Pv->R.B3*Pv->R.M3;
  G = Pv->R.B2*Pv->R.B2 + Pv->R.B3*Pv->R.B3;
  X = Bx*(A*Pv->S*Bx + C) - (A + G)*(Pv->S*p + Pv->R.EN);

  vx = ( Bx*(A*Bx + C*Pv->S) - (Pv->R.M1 + p)*(G + A) );
  vy = ( - (A + G - Bx*Bx*(1.0 - Pv->S*Pv->S))*Pv->R.M2     
	 + Pv->R.B2*(C + Bx*(Pv->S*Pv->R.M1 - Pv->R.EN)) );
  vz = ( - (A + G - Bx*Bx*(1.0 - Pv->S*Pv->S))*Pv->R.M3     
	 + Pv->R.B3*(C + Bx*(Pv->S*Pv->R.M1 - Pv->R.EN)) );

  scrh = vx*Pv->R.M1 + vy*Pv->R.M2 + vz*Pv->R.M3;
  scrh = X*Pv->R.EN - scrh;
  Pv->w = p + scrh/(X*Pv->S - vx);

  if (Pv->w < 0.0) return(0);  /* -- failure -- */

  Pv->vx = vx/X;
  Pv->vy = vy/X;
  Pv->vz = vz/X;
/*
  EXPAND(Pv->Bx = Bx;                            , 
         Pv->By = (R[B2]*X - Bx*vy)/(X*S - vx);  ,
         Pv->Bz = (R[B3]*X - Bx*vz)/(X*S - vx);)
*/
  Pv->Bx = Bx;                                   
  Pv->By = -(Pv->R.B2*(Pv->S*p + Pv->R.EN) - Bx*Pv->R.M2)/A;
  Pv->Bz = -(Pv->R.B3*(Pv->S*p + Pv->R.EN) - Bx*Pv->R.M3)/A;
  
  s  = Bx > 0.0 ? 1.0:-1.0;
  if (side < 0) s *= -1.0;

  Pv->sw = s*sqrt(Pv->w);

  scrh = 1.0/(Pv->S*p +  Pv->R.EN + Bx*Pv->sw);
  Pv->Kx = scrh*(Pv->R.M1 + p + Pv->R.B1*Pv->sw);
  Pv->Ky = scrh*(Pv->R.M2     + Pv->R.B2*Pv->sw);
  Pv->Kz = scrh*(Pv->R.M3     + Pv->R.B3*Pv->sw);

  Pv->K2 = Pv->Kx*Pv->Kx + Pv->Ky*Pv->Ky + Pv->Kz*Pv->Kz;
  return(1); /* -- success -- */
}
/* ************************************************************* */
void GET_ASTATE (Riemann_State *Pa,  Real p, const Real Bx)
/*
 *
 *        Compute states aL and aR (behind fast waves)
 *
 *************************************************************** */
{
  Real vB;
  Real scrh;

  scrh = 1.0/(Pa->S - Pa->vx);

  Pa->u.DN =  Pa->R.DN*scrh;
  Pa->u.B1 =  Bx;                       
  Pa->u.B2 = (Pa->R.B2 - Bx*Pa->vy)*scrh;
  Pa->u.B3 = (Pa->R.B3 - Bx*Pa->vz)*scrh;

  vB       = Pa->vx*Pa->u.B1 + Pa->vy*Pa->u.B2 + Pa->vz*Pa->u.B3;
  Pa->u.EN = (Pa->R.EN + p*Pa->vx - vB*Bx)*scrh;

  Pa->u.M1 = (Pa->u.EN + p)*Pa->vx - vB*Pa->u.B1;
  Pa->u.M2 = (Pa->u.EN + p)*Pa->vy - vB*Pa->u.B2;
  Pa->u.M3 = (Pa->u.EN + p)*Pa->vz - vB*Pa->u.B3;
}

/* ************************************************************* */
void GET_CSTATE (Riemann_State *PaL, Riemann_State *PaR, Real p,
		 CONS_STATE *Uc, const Real Bx)
/*
 *
 *     Compute state  cL and cR across contact mode
 *
 *************************************************************** */
{
  CONS_STATE ua;
  double dK;
  double vxcL, vycL, vzcL, KLBc;
  double vxcR, vycR, vzcR, KRBc;
  double vxc, vyc, vzc, vBc;
  double Bxc, Byc, Bzc, Sa, vxa;
  double scrhL, scrhR;

  GET_ASTATE (PaL, p, Bx);
  GET_ASTATE (PaR, p, Bx);
  dK = (PaR->Kx - PaL->Kx) + 1.e-12;

  Bxc = Bx*dK;
  Byc =   PaR->By*(PaR->Kx - PaR->vx) 
        - PaL->By*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vy - PaL->vy);
  Bzc =   PaR->Bz*(PaR->Kx - PaR->vx) 
        - PaL->Bz*(PaL->Kx - PaL->vx)
        + Bx*(PaR->vz - PaL->vz);
   
  Bxc  = Bx;
  Byc /= dK;
  Bzc /= dK;

  Uc->B1 = Bxc;
  Uc->B2 = Byc;
  Uc->B3 = Bzc;

  KLBc = PaL->Kx*Bxc + PaL->Ky*Byc + PaL->Kz*Bzc;
  KRBc = PaR->Kx*Bxc + PaR->Ky*Byc + PaR->Kz*Bzc;

  scrhL = (1.0 - PaL->K2)/(PaL->sw - KLBc);
  scrhR = (1.0 - PaR->K2)/(PaR->sw - KRBc);

  vxcL = PaL->Kx - Uc->B1*scrhL;  
  vxcR = PaR->Kx - Uc->B1*scrhR;

  vycL = PaL->Ky - Uc->B2*scrhL;  
  vycR = PaR->Ky - Uc->B2*scrhR;

  vzcL = PaL->Kz - Uc->B3*scrhL;
  vzcR = PaR->Kz - Uc->B3*scrhR;

  vxc = 0.5*(vxcL + vxcR);
  vyc = 0.5*(vycL + vycR);
  vzc = 0.5*(vzcL + vzcR);

/*
{
double scrh, Sx;
Sx   = Bx > 0.0 ? 1.0:-1.0;
scrh = vxc*Bxc + vyc*Byc + vzc*Bzc;
scrh = PaR->Ky*sqrt(PaR->w) + PaL->Ky*sqrt(PaL->w) + 
       Sx*scrh*(PaR->Ky - PaL->Ky) - (sqrt(PaR->w) + sqrt(PaL->w))*vyc;

if (fabs(scrh) > 1.e-5){
 printf ("! NOT satisfied, %f\n",scrh);
 exit(1);
}
}
*/

  if (vxc > 0.0) {
    GET_ASTATE (PaL, p, Bx);
    ua  = PaL->u;
    Sa  = PaL->Sa;
    vxa = PaL->vx;
  } else {
    GET_ASTATE (PaR, p, Bx);
    ua  = PaR->u;
    Sa  = PaR->Sa;
    vxa = PaR->vx;
  }
  
  vBc = vxc*Uc->B1 + vyc*Uc->B2 + vzc*Uc->B3;

  Uc->DN = ua.DN*(Sa - vxa)/(Sa - vxc);
  Uc->EN = (Sa*ua.EN - ua.M1 + p*vxc - vBc*Bx)/(Sa - vxc);

  Uc->M1 = (Uc->EN + p)*vxc - vBc*Bx;
  Uc->M2 = (Uc->EN + p)*vyc - vBc*Uc->B2;
  Uc->M3 = (Uc->EN + p)*vzc - vBc*Uc->B3;
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

int checkCurrentSheet (const Prim1DS Wl,const Prim1DS Wr,
			const Cons1DS Ul,const Cons1DS Ur,const Real Bxi)
{
  Real vel2, Bmag2, vB, b2l,b2r,b2,beta;
  Real sgn_dby,sgn_dbz;
  Real sgn_bl,sgn_br;
  int fail = 0;
        
  vel2 = SQR(Wl.Vx) + SQR(Wl.Vy) + SQR(Wl.Vz);
  Bmag2 = SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz);
  vB = Wl.Vx*Bxi + Wl.Vy*Wl.By + Wl.Vz*Wl.Bz;
  b2l = (Bmag2*(1.0 - vel2) + vB*vB);

  vel2 = SQR(Wr.Vx) + SQR(Wr.Vy) + SQR(Wr.Vz);
  Bmag2 = SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz);
  vB = Wr.Vx*Bxi + Wr.Vy*Wr.By + Wr.Vz*Wr.Bz;
  b2r = (Bmag2*(1.0 - vel2) + vB*vB);

  b2 = 0.5*(b2l + b2r);

  beta = 1.0e9;
  if (b2 > 0.0) beta = (Wl.P + Wr.P)/b2;

  sgn_bl = 1.0;
  sgn_br = 1.0;
  if (Wl.By < 0.0) sgn_bl = -1.0;
  if (Wr.By < 0.0) sgn_br = -1.0;
  sgn_dby = sgn_bl * sgn_br;

  sgn_bl = 1.0;
  sgn_br = 1.0;
  if (Wl.Bz < 0.0) sgn_bl = -1.0;
  if (Wr.Bz < 0.0) sgn_br = -1.0;
  sgn_dbz = sgn_bl * sgn_br;

  if (beta < 1.0){
    if (sgn_dby < 0.0) fail = 1;
    if (sgn_dbz < 0.0) fail = 1;
  }

  return(fail);
}

void entropy_flux (const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr, const Real Bx, Real *pFlux)
{
  Real Fl, Fr;
  Real USl, USr;
  Real WSl, WSr;
  Real Uhll, Fhll;
  Real Pl, Pr;
  Real Sl, Sr;
  Real Sla, Sra;
  Real dS_1;
  int wave_speed_fail;

  wave_speed_fail = 0;
	
/*--- Step 1. ------------------------------------------------------------------
 * Compute the max and min wave speeds used in Mignone 
 */
  getMaxSignalSpeeds_pluto(Wl,Wr,Bx,&Sl,&Sr);
	
  if (Sl != Sl) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr =  1.0;
  }
	
  if (Sr != Sr) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
	
  if (Sl < -1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }
  if (Sr > 1.0) {
    wave_speed_fail = 1;
    printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
    Sl = -1.0;
    Sr = 1.0;
  }

/*--- Step 1a. -----------------------------------------------------------------
 * If PLUTO wavespeeds are bad, fall back to the estimate used in ECHO
 */
  if (wave_speed_fail){
    getMaxSignalSpeeds_echo (Wl,Wr,Bx,&Sla,&Sra);
	
    if (Sla != Sla) {
      printf("[hlle_sr_mhd]: NaN in Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra =  1.0;
    }
	
    if (Sra != Sra) {
      printf("[hlle_sr_mhd]: NaN in Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
	
    if (Sla < -1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sl %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }
    if (Sra > 1.0) {
      printf("[hlle_sr_mhd]: Superluminal Sr %10.4e %10.4e\n",Sl,Sr);
      Sla = -1.0;
      Sra = 1.0;
    }

    Sl = Sla;
    Sr = Sra;

  }

  /* compute L/R fluxes */
  WSl = Wl.P*pow(Wl.d,1.0-Gamma);
  WSr = Wr.P*pow(Wr.d,1.0-Gamma);
  USl = WSl * Ul.d/Wl.d;
  USr = WSr * Ur.d/Wr.d;
  Fl = USl * Wl.Vx;
  Fr = USr * Wr.Vx;

  if(Sl >= 0.0){
    *pFlux = Fl;
    return;
  }
  else if(Sr <= 0.0){
    *pFlux = Fr;
    return;
  }
  else{
    /* Compute HLL average state */
    dS_1 = 1.0/(Sr - Sl);
    *pFlux = (Sr*Fl  - Sl*Fr  + Sl*Sr*(USr  - USl)) * dS_1;
    return;
  }
}

void flux_LR(Cons1DS U, Prim1DS W, Cons1DS *flux, Real Bx, Real* p)
{
  Real wtg2, pt, g, g2, g_1,g_2, h, gmmr, theta;
  Real bx, by, bz, vB, b2, Bmag2;
	
  /* calcUlate enthalpy */
	
  theta = W.P/W.d;
  gmmr = Gamma / Gamma_1;
	
  h = 1.0 + gmmr*theta;
	
  /* calcUlate gamma */
  g   = U.d/W.d;
  g2  = SQR(g);
  g_1 = 1.0/g;
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
  getVChar_pluto(Wr,Bx,&lmr,&lpr);
	
  *low =  MIN(lml, lmr);
  *high = MAX(lpl, lpr);
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
    /*printf("[getVChar]: |v|= %f > 1\n",vsq);*/	
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
    return;
  }
	
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
	
  *lp = MIN(1.0,MAX(lambda[3], lambda[2]));
  *lp = MIN(1.0,MAX(*lp, lambda[1]));
  *lp = MIN(1.0,MAX(*lp, lambda[0]));

  *lm = MAX(-1.0,MIN(lambda[3], lambda[2]));
  *lm = MAX(-1.0,MIN(*lm, lambda[1]));
  *lm = MAX(-1.0,MIN(*lm, lambda[0]));
	
  return;
	
}

void getMaxSignalSpeeds_echo (const Prim1DS Wl, const Prim1DS Wr,
			      const Real Bx, Real* low, Real* high)
{
	
  Real lml,lmr;        /* smallest roots, Mignone Eq 55 */
  Real lpl,lpr;        /* largest roots, Mignone Eq 55 */
  Real al,ar;
	
  getVChar_echo(Wl,Bx,&lml,&lpl);
  getVChar_echo(Wr,Bx,&lmr,&lpr);
	
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



#endif /* MHD */
#undef MAX_ITER

#endif
#endif
