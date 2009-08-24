#ifndef GLOBALS_H
#define GLOBALS_H  
/*==============================================================================
 * FILE: globals.h
 *
 * PURPOSE: Contains following global variables:
 *     CourNo       - The Courant, Friedrichs, & Lewy (CFL) Number
 *     Iso_csound   - Isothermal sound speed
 *     Iso_csound2  - Iso_csound^2
 *     Gamma       - Gamma C_p/C_v
 *     Gamma_1, Gamma_2 - Gamma-1, and Gamma-2
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else. 
 *============================================================================*/

#ifdef MAIN_C

Real CourNo;
#ifdef ISOTHERMAL
Real Iso_csound;
Real Iso_csound2;
#elif defined ADIABATIC
Real Gamma;
Real Gamma_1, Gamma_2;
#endif

GravPotFun_t StaticGravPot = NULL;
CoolingFun_t CoolingFunc = NULL;
#ifdef SELF_GRAVITY
Real four_pi_G, grav_mean_rho;
#endif

#ifdef SHEARING_BOX
Real Omega_0, qshear;
#endif

#ifdef PARTICLES
VGFun_t Integrate_Particles = NULL; /* function pointer to particle integrator */
TSFun_t     get_ts    = NULL;   /* get the stopping time */
WeightFun_t getweight = NULL;   /* get weight function */
#endif

#ifdef EXPLICIT_DIFFUSION
Real eta_Ohm=0.0, eta_Hall=0.0, nu_V=0.0, kappa_T=0.0, chi_C=0.0;
#endif

#ifdef CYLINDRICAL
StaticGravAcc_t x1GravAcc = NULL;
Real *r=NULL, *ri=NULL;
#endif

/*----------------------------------------------------------------------------*/
/* definitions included everywhere except main.c  */

#else /* MAIN_C */

extern Real CourNo;
#ifdef ISOTHERMAL
extern Real Iso_csound;
extern Real Iso_csound2;
#elif defined ADIABATIC
extern Real Gamma;
extern Real Gamma_1, Gamma_2;
#endif

extern GravPotFun_t StaticGravPot;
extern CoolingFun_t CoolingFunc;
#ifdef SELF_GRAVITY
extern Real four_pi_G, grav_mean_rho;
#endif

#ifdef SHEARING_BOX
extern Real Omega_0, qshear;
#endif

#ifdef PARTICLES
extern Real alamcoeff;
extern Real *grrhoa;

extern VGFun_t Integrate_Particles;
extern TSFun_t     get_ts;      /* get the stopping time */
extern WeightFun_t getweight;   /* get weight function */
#endif

#ifdef CONST_GRAVITY
extern Real g;
#endif

#ifdef EXPLICIT_DIFFUSION
extern Real eta_Ohm, eta_Hall, nu_V, kappa_T, chi_C;
#endif

#ifdef CYLINDRICAL
extern StaticGravAcc_t x1GravAcc;
extern Real *r, *ri;
#endif

#endif /* MAIN_C */
#endif /* GLOBALS_H */
