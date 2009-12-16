#ifndef GLOBALS_H
#define GLOBALS_H  
/*==============================================================================
 * FILE: globals.h
 *
 * PURPOSE: Contains global variables:
 *   The first occurence in this file is included in main.c and defines the
 *   variables.  The second is included everywhere else. 
 *============================================================================*/

#ifdef MAIN_C

Real CourNo;                 /* Courant, Friedrichs, & Lewy (CFL) number */
#ifdef ISOTHERMAL
Real Iso_csound;             /* isothermal sound speed */
Real Iso_csound2;            /* isothermal sound speed squared */
#elif defined ADIABATIC
Real Gamma;                  /* adiabatic index (ratio of specific heats) */
Real Gamma_1, Gamma_2;       /* (Gamma)-1 and (Gamma)-2 */
#endif
int myID_Comm_world;   /* Rank (proc ID) in MPI_COMM_WORLD, 0 for single proc */

GravPotFun_t StaticGravPot = NULL;
CoolingFun_t CoolingFunc = NULL;
#ifdef SELF_GRAVITY
Real four_pi_G, grav_mean_rho;    /* 4\pi G and mean density in domain */
#endif

#ifdef SHEARING_BOX
Real Omega_0, qshear; /* orbital frequency and shear parameter dln\Omega/dlnr */
#endif

#ifdef PARTICLES
TSFun_t     get_ts    = NULL;     /* get the stopping time */
WeightFun_t getweight = NULL;     /* get weight function */
#endif

#ifdef EXPLICIT_DIFFUSION
/* coefficients of Ohmic resistivity, Hall conduction, kinematic viscosity,
 * isotropic and anisotropic thermal conduction  */
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
extern Real Iso_csound, Iso_csound2;
#elif defined ADIABATIC
extern Real Gamma, Gamma_1, Gamma_2;
#endif
extern int myID_Comm_world;

extern GravPotFun_t StaticGravPot;
extern CoolingFun_t CoolingFunc;
#ifdef SELF_GRAVITY
extern Real four_pi_G, grav_mean_rho;
#endif

#ifdef SHEARING_BOX
extern Real Omega_0, qshear;
#endif

#ifdef PARTICLES
extern Real alamcoeff, *grrhoa;
extern TSFun_t     get_ts; 
extern WeightFun_t getweight; 
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
