#include "copyright.h"
/*==============================================================================
 * FILE: ionrad_chemistry.c
 *
 * PURPOSE: Contains functions to compute various chemistry updates
 * used in the Krumholz, Stone, & Gardiner (2007) ionizing radiation
 * transfer algorithm. Unless otherwise specified, an input of T is
 * temperature in Kelvin and x is ion fraction (H+ / (H+ + H).
 *
 *   Use of these routines requires that --enable-ion-radiation be set
 *   at compile time.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * recomb_rate_coef(T): e- + H+ -> H rate coefficient in cm^-3 s^-1
 * coll_ion_rate_coef(T): collisional H ionization rate coefficient in
 *      cm^3 s^-1
 * recomb_cool_rate_coef(T): radiative recombination cooling rate
 *      coefficient in erg cm^3 s^-1
 * dmc_cool_rate(x, T): Dalgarno-McCray (1973) atomic gas cooling rate
 *      in erg cm^3 s^-1
 * ki_cooling_rate(x, T): Koyama & Inutsuka (2002) cooling rate, including 
 *      H_2 and CO cooling, in erg cm^3 s^-1. Only applies to molecular
 *      gas.
 * ki_heating_rate(x, T): Koyama & Inutsuka (2002) heating rate, in erg 
 *      cm^-3 s^-1. Only applies to non-ionized gas.
 * osterbrock_cool_rate(T): cooling rate computed using transitions
 *      of first two ionized states of O, N, and Ne, as taken from
 *      Osterbrock & Ferland 2006, in erg cm^3 s^-1
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "prototypes.h"

/* Units and constants */
#define EV  1.602e-12
#define KB  1.38e-16
#define C   3.0e10
#define H   6.63e-27
#define IHI (13.6*EV)    /* H ionization potential */
#define GAMMAKI 2.0e-26  /* Koyama & Inutsuka (2002) heating rate, in erg/s */

/* MacDonald & Bailey cooling function table */
static Real xmat[] = {-0.133, 0.105, 0.452, 0.715, 0.901
		      , 1.030, 1.082, 1.174, 1.257, 1.362
		      , 1.448, 1.523, 1.569, 1.582, 1.539
		      , 1.430, 1.275, 1.168, 1.092, 1.019
		      , 1.000, 1.004, 1.008, 0.987, 0.905
		      , 0.738, 0.603, 0.555, 0.552, 0.554
		      , 0.552, 0.535, 0.425, 0.275, 0.251
		      , 0.232, 0.247, 0.283, 0.322, 0.363
		      , 0.397};

/* Element and ion abundances, and mean Gaunt factor,
   from Osterbrock & Ferland (2006). */
#define XO  7.0e-4
#define XNE 9.0e-5
#define XNI 9.0e-5
#define XSINGLE 0.8
#define XDOUBLE (1.0-XSINGLE)
#define GFF 1.3

/* Tables of collision strengths, wavelengths, and statistical weights
   of the ground state for transitions from Osterbrock & Ferland
   (2006). We use a different table for each species, and cgs units
   everywhere. */

/* OII, 4S->2D, 4S->2P */
static int nOII = 2;
static Real omegaOII[] = { 1.34, 0.4 };
static Real lambdaOII[] = { 3.73e-5, 2.47e-5 };
static int wgt1OII[] = { 4, 4 };

/* OIII, 3P->1D, 3P->1S, 3P0->3P1, 3P0->3P2 */
static int nOIII = 4;
static Real omegaOIII[] = { 2.29, 0.29, 0.55, 0.27 };
static Real lambdaOIII[] = { 5.0e-5, 2.33e-5, 8.84e-3, 3.27e-3 };
static int wgt1OIII[] = { 9, 9, 1, 1 };

/* NeII, 1P(1/2)->1P(3/2) */
static int nNeII = 1;
static Real omegaNeII[] = { 0.28 };
static Real lambdaNeII[] = { 1.28e-3 };
static int wgt1NeII[] = { 2 };

/* NeIII, 3P->1D, 3P->1S, 3P0->3P1, 3P0->3P2 */
static int nNeIII = 4;
static Real omegaNeIII[] = { 1.36, 0.15, 0.24, 0.21 };
static Real lambdaNeIII[] = { 3.95e-5, 1.80e-5, 3.60e-3, 1.09e-3 };
static int wgt1NeIII[] = { 9, 9, 1, 1 };

/* NII, 3P->1D, 3P->1S, 3P0->3P1, 3P0->3P2 */
static int nNII = 4;
static Real omegaNII[] = { 2.64, 0.29, 0.41, 0.27 };
static Real lambdaNII[] = { 6.55e-5, 3.07e-5, 2.06e-2, 7.65e-3 };
static int wgt1NII[] = { 9, 9, 1, 1 };

/* NIII, 1P(1/2)->1P(3/2) */
static int nNIII = 1;
static Real omegaNIII[] = { 1.45 };
static Real lambdaNIII[] = { 5.73e-3 };
static int wgt1NIII[] = { 2 };


Real recomb_rate_coef(Real T) {
  /* Osterbrock (1989), pg. 19, rough fit to table values. Note that
   * the table contains an error in the 20,000 K column for alpha_A or
   * alpha_B. Fit is taken from Rijkhorst, Plewa, Dewey, & Mellema
   * (2005), eqn. 18.
   */
#ifdef NORECOMBINATION
  return(0.0);
#else
#  ifndef FIXEDALPHAB
  return(2.59e-13*pow(T/1.0e4, -0.7));
#  else
  return(2.59e-13);
#  endif
#endif
}

Real coll_ion_rate_coef(Real T) {
  /* Tenorio-Tagle et al. 1986, eqn 8 */
#ifdef NORECOMBINATION
  return(0.0);
#else
  return(5.84e-11*sqrt(T)*exp(-IHI/(KB*T)));
#endif
}

Real recomb_cool_rate_coef(Real T) {
  /* Recombination cooling rate from Osterbrock (1989),
     pg. 50-51. Numerical values are a linear fit (in log space) to
     Table 3.2, column 3. */
#ifdef NORECOMBINATION
  return(0.0);
#else
  if (T < 100.0) {
    return(0.0);
  } else {
    return(6.11e-10*pow(T,-0.89)*KB*T);
  }
#endif
}

#define SCALEFACTOR 1.0e-23 /* Scaling to physical units */
Real dmc_cool_rate(Real x, Real T) {
#ifndef NORECOMBINATION
  /* Equilibrium cooling rate from Dalgarno & McCray (1972) */
  /* Implementation follows Jim Stone's dmc routine */
  Real le, lh, u, u2, om, dom, p1;
  Real qq2, qt1, qt2, qt3, qt4, xu1, xu2, xu3, xu4, tlost, tcool;
  int ipps, jaug;

  if (x < 1.0e-3) x = 1.0e-3;

  /* Electron impact excitation luminosity (eqn 3-10) */
  le = 0.0;
  if (T > 10) {
    le += 2.96e-23/sqrt(T) * exp(-92.0/T);
    if (T > 50) {
      le += 6.08e-23/sqrt(T) * exp(-413.0/T)
	+ 3.52e-23/sqrt(T) *
	(exp(-554.0/T) + 1.3* exp(-961.0/T));
      if (T > 2.0e4) {
	le = le + 4.14e-26*sqrt(T) * exp(-22700.0/T)
	  +  7.13e-26*sqrt(T)*(1.0-2.7e-9*T*T)
	  *exp(-27700.0/T);
      }
    }
  }

  /* Hydrogen cooling */
  if (T > 50.0) {
    lh = 2.37e-27*exp(-413/T)
      + 3.52e-27*(exp(-554.0/T)+1.4*exp(-961.0/T));
  } else {
    lh = 0.0;
  }

  /* Neutral cooling */
  u = T/157890. < 3.16 ? T/157890. : 3.16;
  u2 = u*u;
  om=.6098+1.489*u+.50755*u2-.38145*u*u2+.10196*u2*u2
    -.01007*u*u2*u2;
  dom=(1.489+2.*.50755*u-3.*.38145*u2+4.*.10196*u2*u
       -5.*.01007*u2*u2)/157890.;
  p1 = 0.0;
  if (T > 1.0e4)
    p1 = 0.5*1.41e-16*om*exp(-118000./T) / sqrt(T);

  /* Cooling in various regimes: < 100 K */
  if (T < 100.0)
    /* return(0.0); */
    return(x*le + lh + 
	   (1.0-x)*p1);
  if (T < 1.0e4)
    return(SCALEFACTOR*x*2.8347e-10*pow(T - 1.0e+02, 2.3562)
	   + x*le + lh
	   +(1.0-x)*p1);
  if (T > 1.27717e8)
    return(x*2.3988e-04*sqrt(T));

  tlost = log(T)/log(10.0);
  ipps = (int) floor(10.0*tlost) - 38;
  if (ipps > 41) ipps = 41;
  jaug = 2 > ipps ? 2 : ipps;
  qq2 = 3.8 + 0.1*jaug;
  qt2 = tlost - qq2;
  qt3 = qt2 - 0.1;

  if ((jaug == 2) || (jaug == 41)) {
    tcool = (xmat[jaug-1]*qt2 - xmat[jaug-2]*qt3)*10.0;
    return(SCALEFACTOR*pow(10.0, tcool)*x + (1.0-x)*p1);
  }

  qt1 = qt2 + 0.1;
  qt4 = qt3 - 0.1;

  xu1 = qt2*qt3*qt4/6.0e-03;
  xu2 = qt1*qt3*qt4/2.0e-03;
  xu3 = qt1*qt2*qt4/2.0e-03;
  xu4 = qt1*qt2*qt3/6.0e-03;

  tcool = -xmat[jaug-3]*xu1 + xmat[jaug-2]*xu2 -
    xmat[jaug-1]*xu3 + xmat[jaug]*xu4;
  return(SCALEFACTOR*pow(10.0,tcool)*x + (1.0-x)*p1);
#else
  return(0.0);
#endif
}
#undef SCALEFACTOR

Real ki_cool_rate(Real T) {
  return(GAMMAKI * (1.0e7 * exp(-118400.0/(T + 1000.0)) +
	  14*sqrt(T)*exp(-92.0/T)));
}

Real ki_heat_rate() {
  return(GAMMAKI);
}


Real osterbrock_cool_rate(Real T) {

  /* Osterbrock & Ferland cooling function table. The 101 values in the
     table give the log of the cooling rate in erg/s at temperatures
     from 1 K to 10^6 K, evenly spaced in log T. */
#define OBCOOLTABLO 0.0
#define OBCOOLTABHI 6.0
#define NOBCOOLTAB  101
#define LN10        2.3025851
  static Real obcooltab[] = {
    -26.733768, -26.703768, -26.673768, -26.643768, -26.613768,
    -26.583768, -26.553768, -26.523768, -26.493768, -26.463768,
    -26.433765, -26.403741, -26.373595, -26.342921, -26.310465,
    -26.273212, -26.225614, -26.160265, -26.071200, -25.958085,
    -25.826891, -25.686072, -25.542386, -25.399492, -25.258842,
    -25.121202, -24.987655, -24.859777, -24.739274, -24.627519,
    -24.525290, -24.432759, -24.349630, -24.275341, -24.209232,
    -24.150643, -24.098951, -24.053555, -24.013850, -23.979221,
    -23.949059, -23.922804, -23.899985, -23.880253, -23.863380,
    -23.849248, -23.837808, -23.829051, -23.822968, -23.819529,
    -23.818674, -23.820303, -23.824282, -23.830444, -23.838596,
    -23.848499, -23.859817, -23.871978, -23.883877, -23.893356,
    -23.896435, -23.886613, -23.855321, -23.794950, -23.703720,
    -23.587876, -23.458501, -23.326438, -23.199497, -23.082246,
    -22.976856, -22.884004, -22.803510, -22.734742, -22.676850,
    -22.628894, -22.589926, -22.559028, -22.535335, -22.518046,
    -22.506425, -22.499806, -22.497586, -22.499223, -22.504231,
    -22.512176, -22.522668, -22.535361, -22.549943, -22.566136,
    -22.583687, -22.602370, -22.621979, -22.642323, -22.663228,
    -22.684530, -22.706076, -22.727716, -22.749311, -22.770719,
    -22.791806 };
  int idx;
  Real idxrl, wgt, tabinterp;
  int n;
  Real energy, q12, cool = 0;

  /* Find index on table */
  idxrl = (NOBCOOLTAB-1) * (log10(T)-OBCOOLTABLO) / 
    (OBCOOLTABHI-OBCOOLTABLO);
  idx = (int) idxrl;

  if ((idx >= 0) && (idx < NOBCOOLTAB-1)) {

    /* We're on the table, so interpolate */
    wgt = idxrl - idx;
    tabinterp = (1.0-wgt)*obcooltab[idx] + wgt*obcooltab[idx+1];
    return(exp(LN10*tabinterp));

  } else {

    /* We're off the table, so do the calculation */

    /* Loop over ions, leaving out common factors for now */
    /* OII */
    for (n=0; n<nOII; n++) {
      energy = H*C/lambdaOII[n];
      q12 = omegaOII[n]/wgt1OII[n]*exp(-energy/(KB*T));
      cool += XO*XSINGLE*q12*energy;
    }
    /* OIII */
    for (n=0; n<nOIII; n++) {
      energy = H*C/lambdaOIII[n];
      q12 = omegaOIII[n]/wgt1OIII[n]*exp(-energy/(KB*T));
      cool += XO*XDOUBLE*q12*energy;
    }
    /* NeII */
    for (n=0; n<nNeII; n++) {
      energy = H*C/lambdaNeII[n];
      q12 = omegaNeII[n]/wgt1NeII[n]*exp(-energy/(KB*T));
      cool += XNE*XSINGLE*q12*energy;
    }
    /* NeIII */
    for (n=0; n<nNeIII; n++) {
      energy = H*C/lambdaNeIII[n];
      q12 = omegaNeIII[n]/wgt1NeIII[n]*exp(-energy/(KB*T));
      cool += XNE*XDOUBLE*q12*energy;
    }
    /* NII */
    for (n=0; n<nNII; n++) {
      energy = H*C/lambdaNII[n];
      q12 = omegaNII[n]/wgt1NII[n]*exp(-energy/(KB*T));
      cool += XNI*XSINGLE*q12*energy;
    }
    /* NIII */
    for (n=0; n<nNIII; n++) {
      energy = H*C/lambdaNIII[n];
      q12 = omegaNIII[n]/wgt1NIII[n]*exp(-energy/(KB*T));
      cool += XNI*XDOUBLE*q12*energy;
    }
    
    /* Multiply by common factors */
    cool *= 8.63e-6/sqrt(T);

    /* Add free-free cooling */
    cool += 1.42e-27*GFF*sqrt(T);

    return(cool);
  }
}
