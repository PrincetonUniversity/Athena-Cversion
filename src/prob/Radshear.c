#include "copyright.h"
/*==============================================================================
 * FILE: radMHD2d.c
 *
 * PURPOSE: Problem generator to test the radiation MHD code. 
 *  Development is underway. To be modified later 
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */
/*
void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
*/


/* For transfer module */
#ifdef RADIATION_TRANSFER

static Real eps0;

static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);

#endif

#ifdef SHEARING_BOX
static Real ShearingPot(const Real x1, const Real x2, const Real x3);

#endif


void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
	Real angle, x1, x2, x3, radius;
	angle = 3.0 * PI/ 180;

/* Parse global variables of unit ratio */

  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");

	/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
#ifdef SHEARING_BOX
	Omega_0 = par_getd_def("problem","omega",2.0);
	qshear  = par_getd_def("problem","qshear",1.5);
#endif

  Real alpha, miu, pressure, velocity_y;
  alpha = 20.0;
  miu = 0.0;
	

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

	iu = pGrid->ie;
	il = pGrid->is;

	ju = pGrid->je;
	jl = pGrid->js;


	


  for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {	
        for (i=il; i<=iu; i++) {
			cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
			radius=sqrt(x1*x1+x2*x2);
			
			
			
			
			pGrid->U[k][j][i].Er = 0.01 + 0.1 * exp(-(alpha*(x2-miu))*(alpha*(x2-miu)));
/*
			pGrid->U[k][j][i].Er = 1.0;
 */
 /*	pGrid->U[k][j][i].Fr2 = 2.0 * 3.0 * alpha * alpha * (x2-miu) * pGrid->U[k][j][i].Er / pGrid->U[k][j][i].Sigma_t;
			 */

			pGrid->U[k][j][i].Fr1 = 0.0;
			/* background Er */
			
		
			pGrid->U[k][j][i].d  = 1.0;
			pGrid->U[k][j][i].M1 = 0.0;
			pGrid->U[k][j][i].M2 = -0.0 * 1.5 * x1;
#ifdef SHEARING_BOX
#ifdef FARGO
			pGrid->U[k][j][i].M2 = 0.0;	
			velocity_y = -qshear * Omega_0 * x1;
#else
			pGrid->U[k][j][i].M2 = -qshear * Omega_0 * x1;
#endif
#endif
			pGrid->U[k][j][i].M3 = 0.0;
			
			pressure = pow(pGrid->U[k][j][i].Er, 0.25) * pGrid->U[k][j][i].d;

			
			pGrid->U[k][j][i].E = pressure/(Gamma - 1.0)+0.5 * (pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 
														   + pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1) / pGrid->U[k][j][i].d;
			
			
			pGrid->U[k][j][i].Edd_11 = 1.0/3.0; 
			pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
			pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
			pGrid->U[k][j][i].Sigma_t = 40.0;
			pGrid->U[k][j][i].Sigma_a = 40.0;


#ifdef RADIATION_MHD
	 /* interface magnetic field */


          pGrid->B1i[k][j][i] = Bx0;
          pGrid->B2i[k][j][i] = By0;
          pGrid->B3i[k][j][i] = 0.0;

	  pGrid->U[k][j][i].E += B0 * B0 / 2.0;
        
#endif
			
					pGrid->U[k][j][i].Fr2 = velocity_y * pGrid->U[k][j][i].Er * (1.0 + pGrid->U[k][j][i].Edd_22) / (Crat * pGrid->U[k][j][i].d);


	
	 	
        }
      }
    }
	
/* Boundary condition for trasnfer module */
	
#ifdef SHEARING_BOX
	ShearingBoxPot = ShearingPot;
	
#endif
	


/* data for radiation transfer method */
#ifdef RADIATION_TRANSFER
  RadGridS *pRG = (pDomain->RadGrid);


  int nf=pRG->nf, nang=pRG->nang;
  int ifr, l, m;

  
/* We do not need to iniatialize tau */

/* Read problem parameters. */


/* ------- Initialize boundary emission ---------------------------------- */

  for(ifr=0; ifr<nf; ifr++)
    for(l=0; l<4; l++) 
      for(m=0; m<nang; m++) {
	pRG->r1imu[ifr][pRG->ks][0][l][m] = 0.0;
	pRG->l1imu[ifr][pRG->ks][0][l][m] = 0.0;
      }
  for(j=pRG->js; j<=pRG->je; j++) {
/* incident radiation at left boundary */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<nang; m++) {
	  pRG->l1imu[ifr][pRG->ks][j][0][m] = 0.0;
	  pRG->l1imu[ifr][pRG->ks][j][2][m] = 0.0;
      }
/* incident radiation at right boundary */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<=nang; m++) {
	  pRG->r1imu[ifr][pRG->ks][j][1][m] = 0.0;
	  pRG->r1imu[ifr][pRG->ks][j][3][m] = 0.0;
      }
  }

  for(i=pRG->is; i<=pRG->ie; i++) {
/* incident radiation at upper and lower boundaries */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<nang; m++) {
/* lower boundary is tau=0, no irradiation */

	pRG->r2imu[ifr][pRG->ks][i][2][m] = 0.0;
	pRG->r2imu[ifr][pRG->ks][i][3][m] = 0.0;
      }
  }

get_thermal_source = Thermal_B;
get_thermal_fraction = const_eps;
get_total_opacity = transfer_opacity;

#endif


  return;
}




/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

#ifdef SHEARING_BOX

static Real ShearingPot(const Real x1, const Real x2, const Real x3)
{       
	Real phi=0.0;
#ifndef FARGO
	phi -= qshear*SQR(Omega_0*x1);
#endif  
	return phi;
}

#endif

void problem_write_restart(MeshS *pM, FILE *fp)
{

	


  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
	
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{

	


  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}

void Userwork_in_formal_solution(DomainS *pD)
{
	return;
}



/*=========================== PRIVATE FUNCTIONS ==============================*/

/*------------------------------------------------------------------------------
 * ran2: extracted from the Numerical Recipes in C (version 2) code.  Modified
 *   to use doubles instead of floats. -- T. A. Gardiner -- Aug. 12, 2003
 */

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define RNMX (1.0-DBL_EPSILON)

/* Long period (> 2 x 10^{18}) random number generator of L'Ecuyer
 * with Bays-Durham shuffle and added safeguards.  Returns a uniform
 * random deviate between 0.0 and 1.0 (exclusive of the endpoint
 * values).  Call with idum = a negative integer to initialize;
 * thereafter, do not alter idum between successive deviates in a
 * sequence.  RNMX should appriximate the largest floating point value
 * that is less than 1. 
 */

double ran2(long int *idum)
{
  int j;
  long int k;
  static long int idum2=123456789;
  static long int iy=0;
  static long int iv[NTAB];
  double temp;

  if (*idum <= 0) { /* Initialize */
    if (-(*idum) < 1) *idum=1; /* Be sure to prevent idum = 0 */
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) { /* Load the shuffle table (after 8 warm-ups) */
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;                 /* Start here when not initializing */
  *idum=IA1*(*idum-k*IQ1)-k*IR1; /* Compute idum=(IA1*idum) % IM1 without */
  if (*idum < 0) *idum += IM1;   /* overflows by Schrage's method */
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2; /* Compute idum2=(IA2*idum) % IM2 likewise */
  if (idum2 < 0) idum2 += IM2;
  j=(int)(iy/NDIV);              /* Will be in the range 0...NTAB-1 */
  iy=iv[j]-idum2;                /* Here idum is shuffled, idum and idum2 */
  iv[j] = *idum;                 /* are combined to generate output */
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX; /* No endpoint values */
  else return temp;
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef RNMX


/* Function for transfer module */
#ifdef RADIATION_TRANSFER

static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
	Real density, vx, vy, vz, energy, pressure, T, B;
	density = pG->U[k][j][i].d;
	vx = pG->U[k][j][i].M1/density;
	vy = pG->U[k][j][i].M2/density;
	vz = pG->U[k][j][i].M3/density;
	energy = pG->U[k][j][i].E;
	pressure = (energy - 0.5 * density * (vx* vx + vy * vy + vz * vz)) * (Gamma - 1.0);
#ifdef RADIATION_MHD

	pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
	T = pressure / (density * R_ideal);
	B = T * T * T * T / (4.0 * PI);

  return B;
}

static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k)
{
	Real eps;
	eps = pG->U[k][j][i].Sigma_a / pG->U[k][j][i].Sigma_t;

	return eps;
  
}

static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k)
{

  return (pG->U[k][j][i].Sigma_t);
  
}

#endif
