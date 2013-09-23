#include "../copyright.h"
/*==============================================================================
 * FILE: utils_rad.c
 *
 * PURPOSE: contains misc. functions require for computation of rad. transfer
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *  interp_quad_chi()       -- 2nd order representation of opacity with slope lim.
 *  interp_quad_source_slope_lim() -- FS integ. weights with quad. interp. and
 *                                     slope limiting.
 *  interp_quad_source()    -- FS integ. weights with quad. interp. (obsolete)
 *  interp_linear_source()  -- FS integ. weights with linear interp. (obsolete)
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define TAUMIN 1.0E-6

#ifdef RADIATION_TRANSFER

/*=========================== PUBLIC FUNCTIONS ===============================*/


/*----------------------------------------------------------------------------*/
/*! \fn void interp_quad_chi(Real chi0, Real chi1, Real chi2, Real *chi)
 * Compute weights using parabolic interpolation of opacity.  Uses Bezier
 * interpolation scheme discussed in Auer (2003, ASP, 228, 3).  Method is
 * nearly equivalent to Hayek et al. (2010, A&A, 517, 49).  This version assumes
 * that the grid spacing is fixed. */
void interp_quad_chi(Real chi0, Real chi1, Real chi2, Real *chi)
{

  Real chimin, chimax, chic;

  chic = chi1 - 0.25 * (chi2 - chi0);

  /*chimax = MAX(chi0,chi1);
  chimin = MIN(chi0,chi1);
*/
  /* use standard interp if chimin < chic < chimax */
  /*  if ((chic >= chimin) && (chic <= chimax)) {
 */
   if ((chi0-chic)*(chi1-chic) <= 0.0) {
   (*chi) = 0.4166666666666667 * chi0 + 0.6666666666666667 * chi1 -  0.0833333333333333 * chi2;
    /*(*chi)= (5.0 * chi0 + 8.0 * chi1 - chi2) / 12.0;
  */
	 /* chic = chi1 */  
  } else {
    (*chi) = 0.3333333333333333 * chi0 + 0.6666666666666667 * chi1;
  }

}

/*----------------------------------------------------------------------------*/
/* \fn void interp_quad_source_slope_lim(Real dtaum, Real dtaup, Real *edtau, Real *a0,
 *			                 Real *a1, Real *a2, Real S0, Real S1, Real S2)
 * Copute weights using parabolic interpolation of source function.  Uses Bezier
 * interpolation scheme discussed in Auer (2003, ASP, 228, 3).  Method is
 * nearly equivalent to Hayek et al. (2010, A&A, 517, 49).  In addition to
 * improving stability, it ensures a positive intensity as long as source
 * terms are positive.*/
void interp_quad_source_slope_lim(Real dtaum, Real dtaup, Real *edtau, Real *a0,
			Real *a1, Real *a2, Real S0, Real S1, Real S2)
{
  Real c0, c1, c2;
  Real dtaus, dtaus1, dtaum1, dtaup1, dtausp1, dtaum2;
  Real Sc, taurat;

/* expand to 2nd order in dtaum */
  if (dtaum < TAUMIN) {
    (*edtau) = 1.0 - dtaum * (1.0 - 0.5* dtaum);
    taurat = dtaup / dtaum;
    Sc =  S1 - 0.5 / (1.0 + taurat) *  (taurat * (S1 - S0) + (S2 - S1) / taurat);
    if ((S0-Sc)*(S1-Sc) <= 0.0) {      
      (*a0) = dtaum * (taurat / (1.0 + taurat) - 0.5 * dtaum);
      (*a1) = dtaum * 0.5 * (1.0 + taurat) / taurat;
      (*a2) = - dtaum * 0.5 / (taurat * (1.0 + taurat));
    } else {
      (*a0) = 0.666666666666666667 * dtaum;
      (*a1) = 0.333333333333333333 * dtaum;
      (*a2) = 0.0;
    }
  } else {
/* use standard expresions */
    dtaus = dtaum + dtaup;
    dtaus1 = 1.0 / dtaus;
    dtaup1 = 1.0 / dtaup; 
    dtaum1 = 1.0 / dtaum;
    dtausp1 = dtaus1 * dtaup1;
    dtaum2 = dtaum * dtaum;
    
    (*edtau) = exp(-dtaum);
    
    c0 = 1.0 - (*edtau);
    c1 = dtaum - c0;
    c2 = dtaum2 - 2.0 * c1;

    Sc = S1 - 0.5 * (dtaup  * (S1 - S0) * dtaus1 + dtaum2 * (S2 - S1) * dtausp1);

/* use standard interp if Smin < Sc < Smax */
    if ((S0-Sc)*(S1-Sc) <= 0.0) {
      (*a0) = c0 + (c2 - (dtaus + dtaum) * c1) * dtaum1 * dtaus1;
      (*a1) = (dtaus * c1 - c2) * dtaum1 * dtaup1;
      (*a2) = (c2 - dtaum * c1) * dtausp1;
/* Sc = S1 if S1 is not an extremum */  
    } else {
      (*a1)  = c2 / dtaum2;
      (*a0)  = c0 - (*a1);
      (*a2)  = 0.0;
    }
/* Sc = S0 if S1 is an extremum */
  /*} else if (dSp * dSm < 0.0) {
    (*a0)  = c0 + (c2 - 2.0 * dtaum * c1) / dtaum2;
    (*a1)  = (2.0 * dtaum * c1 - c2) / dtaum2;
    (*a2)  = 0.0;
    }*/  
  }
}

/*----------------------------------------------------------------------------*/
/* \fn void interp_quad_source(Real dtaum, Real dtaup, Real *edtau, Real *a0,
 *       		       Real *a1, Real *a2)
 *  Computes quadrature weights for integration of the radiation transfer
 *  eq. along one ray in one gridzone using quadratic interpolation.  Retained only 
 *  for testing/debugging.  Replaced by interp_quad_source_slope_lim(). */
void interp_quad_source(Real dtaum, Real dtaup, Real *edtau, Real *a0,
			Real *a1, Real *a2)
{
  Real c0, c1, c2;
  Real dtaus, dtausp, dtausm, dtaum2;

  dtaus  = dtaum + dtaup;
  dtausp = dtaus * dtaup;
  dtausm = dtaus * dtaum;
  dtaum2 = dtaum * dtaum;


  (*edtau) = exp(-dtaum);

  c0 = 1.0 - (*edtau);
  c1 = dtaum - c0;
  c2 = dtaum2 - 2.0 * c1;

  (*a0) = c0 + (c2 - (dtaup + 2.0 * dtaum) * c1) / dtausm;
  (*a1) = (dtaus * c1 - c2) / (dtaum * dtaup);
  (*a2) = (c2 - dtaum * c1) / dtausp;
}


/*----------------------------------------------------------------------------*/
/* \fn void interp_linear_source(Real dtaum, Real *edtau, Real *a0, Real *a1)
 *  Computes quadrature weights for integration of the radiation transfer
 *  eq. along one ray in one gridzone using linear interpolation.  Retained only 
 *  for testing/debugging.  Replaced by interp_quad_source_slope_lim(). */
void interp_linear_source(Real dtaum, Real *edtau, Real *a0, Real *a1)
{
  Real c0, c1;

  (*edtau) = exp(-dtaum);

  c0 = 1.0 - (*edtau);
  c1 = dtaum - c0;
  
  (*a0) = c0 - c1 / dtaum;
  (*a1) = c1 / dtaum;

}

/*----------------------------------------------------------------------------*/
/*! \fn void interp_quad_chi_new(Real chi0, Real chi1, Real chi2, Real *chi12, Real *chi32)
 *  Experimental version of interp_quad_chi() */
void interp_quad_chi_new(Real chi0, Real chi1, Real chi2, Real *chi12, Real *chi32)
{

  Real chic12, chic32, dchi;

  dchi = 0.25 * (chi2 - chi0);
  chic12 = chi1 - dchi;
  chic32 = chi1 + dchi;

  /* use standard interp if chimin < chic < chimax */

  if ((chi0-chic12)*dchi <= 0.0) {
    (*chi12) = 0.4166666666666667 * chi0 + 0.6666666666666667 * chi1 -  0.0833333333333333 * chi2;
    /* chic = chi1 */  
  } else {
    (*chi12) = 0.3333333333333333 * chi0 + 0.6666666666666667 * chi1;
  }
  if ((chic32-chi2)*dchi <= 0.0) {
    (*chi32) = 0.4166666666666667 * chi2 + 0.6666666666666667 * chi1 -  0.0833333333333333 * chi0;
    /* chic = chi1 */  
  } else {
    (*chi32) = 0.3333333333333333 * chi2 + 0.6666666666666667 * chi1;
  }
}


#endif /* RADIATION_TRANSFER */


