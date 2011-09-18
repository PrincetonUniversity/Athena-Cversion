#include "copyright.h"
/*==============================================================================
 * FILE: hgb.c
 *
 * PURPOSE:  Problem generator for 3D shearing sheet.  Based on the initial
 *   conditions described in "Local Three-dimensional Magnetohydrodynamic
 *   Simulations of Accretion Disks" by Hawley, Gammie & Balbus, or HGB.
 *
 * Several different field configurations and perturbations are possible:
 *
 *  ifield = 0 - uses field set by choice of ipert flag
 *  ifield = 1 - Bz=B0sin(kx*x1) field with zero-net-flux [default] (kx input)
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
 *
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *  ipert = 3 - J&G vortical shwave (hydro test)
 *  ipert = 4 - nonlinear density wave test of Fromang & Papaloizou
 *  ipert = 5 - 2nd MHD shwave test of JGG (2008) -- their figure 9
 *  ipert = 6 - 3rd MHD shwave test of JGG (2008) -- their figure 11
 *  ipert = 7 - nonlinear shearing wave test of Heinemann & Papaloizou (2008)
 *
 * To run simulations of stratified disks (including vertical gravity), use the
 * strat.c problem generator.
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Hawley, J. F. & Balbus, S. A., ApJ 400, 595-609 (1992).
 *            Johnson, Guan, & Gammie, ApJSupp, (2008)
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

Real Lx,Ly,Lz; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
#if defined(MHD) || defined(RADIATION_MHD)
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_dEw2(const GridS *pG, const int i, const int j, const int k);
static Real hst_dBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_EB(const GridS *pG, const int i, const int j, const int k);
#endif /* MHD */
static Real hst_T(const GridS *pG, const int i, const int j, const int k);


static Real kappaes;
static Real kappaff;

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
static void Thindiskopacity(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY]);
#endif
#ifndef SHEARING_BOX
static Real Omega_0 = 1.0;
static Real qshear = 1.5;
#endif

#if defined(RADIATION_MHD) || defined(MHD)
static Real betaz;
#endif
static Real pres;

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int ixs,jxs,kxs,i,j,k,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax;
  Real den = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real beta=1.0,B0,kx,ky,kz,amp;
  Real betay, B0z, B0y, Bamp;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */

	/* Parse global variables of unit ratio */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	Prat = par_getd("problem","Pratio");
	Crat = par_getd("problem","Cratio");
	R_ideal = par_getd("problem","R_ideal");
	Ncycle = par_getd_def("problem","Ncycle",5);
  	TOL  = par_getd_def("problem","TOL",1.e-8);
#endif
	
	
/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","Omega",1.0);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif

	betaz=40.0;
	betay=1.0;
	pres = 1.0;

  	B0z = sqrt((double)(2.0*pres/betaz));
  	B0y = 0.0;
	B0 = sqrt(B0z * B0z + B0y * B0y);	

	kappaes = 1953.93;
	kappaff = 0.105252;

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  kx = 16.0*PI/Lx;

	amp = 0.1;
	Bamp = 1.0;


/* Rescale amp to sound speed for ipert 2,3 */


  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

/* Initialize perturbations
 *  ipert = 1 - random perturbations to P and V [default, used by HGB]
 *  ipert = 2 - uniform Vx=amp (epicyclic wave test)
 *  ipert = 3 - vortical shwave (hydro test)
 *  ipert = 4 - Fromang & Papaloizou nonlinear density wave (hydro test)
 *  ipert = 5 & 6 - JGG MHD shwave tests
 *  ipert = 7 - Heinemann & Papaloizou (2008) nonlinear shwave (hydro test)
 */

        rval = amp*(ran2(&iseed) - 0.5);
#ifdef ADIABATIC
        rp = pres*(1.0 + 2.0*rval);
        rd = den;
#else
        rd = den*(1.0 + 2.0*rval);
#endif
/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
        rvx = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvy = 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
        rvz = 0.4*rval*sqrt(pres/den);
      

/* Initialize d, M, and P.  For 3D shearing box M1=Vx, M2=Vy, M3=Vz
 * With FARGO do not initialize the background shear */ 

      pGrid->U[k][j][i].d  = rd;
      pGrid->U[k][j][i].M1 = rd*rvx;
      pGrid->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pGrid->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pGrid->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pGrid->U[k][j][i].E = rp/Gamma_1
        + 0.5*(SQR(pGrid->U[k][j][i].M1) + SQR(pGrid->U[k][j][i].M2) 
             + SQR(pGrid->U[k][j][i].M3))/rd;
#endif

/* Initialize magnetic field.  For 3D shearing box B1=Bx, B2=By, B3=Bz
 *  ifield = 0 - 
 *  ifield = 1 - Bz=B0 sin(x1) field with zero-net-flux [default]
 *  ifield = 2 - uniform Bz
 *  ifield = 3 - B=(0,B0cos(kx*x1),B0sin(kx*x1))= zero-net flux w helicity
 *  ifield = 4 - B=(0,B0/sqrt(2),B0/sqrt(2))= net toroidal+vertical field
 */
#if defined(MHD) || defined(RADIATION_MHD)

        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0z;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0z;
     

#endif /* MHD */
		
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
		pGrid->U[k][j][i].Edd_11 = 1.0/3.0; 
		pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
		pGrid->U[k][j][i].Edd_33 = 1.0/3.0;
		pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
		pGrid->U[k][j][i].Edd_31 = 0.0;
		pGrid->U[k][j][i].Edd_32 = 0.0;
		pGrid->U[k][j][i].Sigma[0] = kappaes * rd;
		pGrid->U[k][j][i].Sigma[1] = kappaff * rd * rd * pow(rp/(rd*R_ideal),-3.5);
		pGrid->U[k][j][i].Sigma[2] = kappaff * rd * rd * pow(rp/(R_ideal*rd),-3.5);
		pGrid->U[k][j][i].Sigma[3] = kappaff * rd * rd * pow(rp/(R_ideal*rd),-3.5);
		
		pGrid->U[k][j][i].Er = 1.0;
		pGrid->U[k][j][i].Fr1 = 0.0;
		pGrid->U[k][j][i].Fr2 = 0.0;
#ifdef SHEARING_BOX
		pGrid->U[k][j][i].Fr2 = -qshear * Omega_0 * x1 * (1.0 + pGrid->U[k][j][i].Edd_22) * pGrid->U[k][j][i].Er/Crat;	
#else
		pGrid->U[k][j][i].Fr2 = -1.5 *  x1 * (1.0 + pGrid->U[k][j][i].Edd_22) * pGrid->U[k][j][i].Er/Crat;
#endif
		pGrid->U[k][j][i].Fr3 = 0.0;
		
		
#endif
		
		
		
		
		
    }
  }}
#if defined(MHD) || defined(RADIATION_MHD)
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
        pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
        pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
#ifdef ADIABATIC
      pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
         + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
      }
    }
  }
#endif /* MHD */

/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
#endif
	
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)	
	/* enroll the opacity function */
	Opacity = Thindiskopacity;
#endif
	
	
/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
    dump_history_enroll(hst_T,"<T>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif

#if defined(MHD) || defined(RADIATION_MHD)
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BxBy, "<-Bx By>");
    dump_history_enroll(hst_EB,"<B^2/2>");

#endif /* MHD */
    frst = 0;
  }

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/


void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ox1(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ix2(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ox2(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ix3(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ox3(VMatFun_t *Mat_BCFun)
{

} 



void problem_write_restart(MeshS *pM, FILE *fp)
{
	
	fwrite(&Gamma,sizeof(Real),1,fp);

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
	fwrite(&Prat,sizeof(Real),1,fp);
	fwrite(&Crat,sizeof(Real),1,fp);
	fwrite(&R_ideal,sizeof(Real),1,fp);
 	fwrite(&kappaes,sizeof(Real),1,fp);
	fwrite(&kappaff,sizeof(Real),1,fp);
#endif	
	
#ifdef SHEARING_BOX
	fwrite(&qshear,sizeof(Real),1,fp);
	fwrite(&Omega_0,sizeof(Real),1,fp);
#endif
#ifdef RADIATION_MHD
	fwrite(&betaz,sizeof(Real),1,fp);
	fwrite(&pres,sizeof(Real),1,fp);
#endif
#ifdef MATRIX_MULTIGRID 
	fwrite(&Ncycle,sizeof(Real),1,fp);
	fwrite(&TOL,sizeof(Real),1,fp);
#endif
	
  return;
}

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
	

/* Read Omega, and with viscosity and/or resistivity, read eta_Ohm and nu */
	fread(&Gamma,sizeof(Real),1,fp);
	
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)	
	fread(&Prat,sizeof(Real),1,fp);
	Prat = 376.255;
	fread(&Crat,sizeof(Real),1,fp);
	fread(&R_ideal,sizeof(Real),1,fp);
	fread(&kappaes,sizeof(Real),1,fp);
	fread(&kappaff,sizeof(Real),1,fp);
#endif
	
#ifdef SHEARING_BOX
    fread(&qshear,sizeof(Real),1,fp);
	fread(&Omega_0,sizeof(Real),1,fp);
#endif
	
#ifdef RADIATION_MHD
	fread(&betaz,sizeof(Real),1,fp);
	fread(&pres,sizeof(Real),1,fp);
#endif
#ifdef MATRIX_MULTIGRID 
	fread(&Ncycle,sizeof(Real),1,fp);
	fread(&TOL,sizeof(Real),1,fp);
#endif

/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
#endif
/* enroll new history variables */

	
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
    dump_history_enroll(hst_T,"<T>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
	
#if defined(MHD) || defined(RADIATION_MHD)
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BxBy, "<-Bx By>");
    dump_history_enroll(hst_EB,"<B^2/2>");
#endif

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
        /* enroll the opacity function */
        Opacity = Thindiskopacity;
#endif

	
	/* Increase the background magnetic field */
	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	GridS *pGrid = pD->Grid;
	int is = pGrid->is, ie = pGrid->ie;
	int js = pGrid->js, je = pGrid->je;
	int ks = pGrid->ks, ke = pGrid->ke;
	
	int i, j, k;
	Real betanew = 5.0;
	Real diffBz = sqrt((double)(2.0*pres/betanew)) - sqrt((double)(2.0*pres/betaz));

	Real pressure, averageP = 0.0;
	Real sumP = 0.0;
	Real Sigma[NOPACITY];
	int ierr, ID,m;
	int count = 0;
	Real x1, x2, x3, velocity_x, velocity_y, velocity_z;
#ifdef RADIATION_MHD
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=is; i<=ie; i++) {
				pressure = pGrid->U[k][j][i].E - 0.5 * (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 + pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2 
						+ pGrid->U[k][j][i].M3 * pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d
						- 0.5 * (pGrid->U[k][j][i].B1c * pGrid->U[k][j][i].B1c + pGrid->U[k][j][i].B2c * pGrid->U[k][j][i].B2c + pGrid->U[k][j][i].B3c * pGrid->U[k][j][i].B3c);
				pressure *= (Gamma - 1.0);

				if(pressure < TINY_NUMBER) pressure = TINY_NUMBER;

				averageP += pressure;
				count++;

			/* Update the opacity */
				Thindiskopacity(pGrid->U[k][j][i].d, pressure/(pGrid->U[k][j][i].d * R_ideal), Sigma, NULL);
				for(m=0;m<NOPACITY;m++){
					pGrid->U[k][j][i].Sigma[m] = Sigma[m];
				}
			


				pGrid->B3i[k][j][i] += diffBz;
				if (k==ke) pGrid->B3i[ke+1][j][i] += diffBz;
				/* Need to update total Energy when magnetic field is modified */
					
				pGrid->U[k][j][i].E += 0.5 * (2.0 * pGrid->U[k][j][i].B3c * diffBz + diffBz * diffBz);
				pGrid->U[k][j][i].B3c += diffBz;
				
				

			}
		}
	}

				averageP /= count;
	
#endif

#ifdef MPI_PARALLEL
		ierr = MPI_Reduce(&averageP,&sumP,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		ID = myID_Comm_world;
		if(ID == 0){
			sumP /= pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
		}
		MPI_Bcast(&sumP,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

		averageP = sumP;
#endif

	/* Initialize such that radiation pressure / gas pressure =125 */
	
	for (k=ks; k<=ke; k++) {
                for (j=js; j<=je; j++) {
	                  for (i=is; i<=ie; i++) {
				cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

				pGrid->U[k][j][i].Er = 125.0 * averageP * 3.0 / Prat;
				
				velocity_x = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
				velocity_y = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
				velocity_z = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;			

#ifdef FARGO
			      velocity_y -= qshear * Omega_0 * x1;

#endif
				pGrid->U[k][j][i].Fr1 = velocity_x * (1.0 + pGrid->U[k][j][i].Edd_11) * pGrid->U[k][j][i].Er/Crat;
				pGrid->U[k][j][i].Fr2 = velocity_y * (1.0 + pGrid->U[k][j][i].Edd_22) * pGrid->U[k][j][i].Er/Crat;
				pGrid->U[k][j][i].Fr3 = velocity_z * (1.0 + pGrid->U[k][j][i].Edd_33) * pGrid->U[k][j][i].Er/Crat;


			}
		}
	}
  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{
  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}

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

/*------------------------------------------------------------------------------
 * UnstratifiedDisk:
 */

static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
void Thindiskopacity(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY])
{
/* Sigma[0-NOPACITY] are: Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE respectively */
/* dSigma[0-2*NOPACITY] are: dSigma_sF/drho, dSigma_aF/drho, dSigma_aP/drho, dSigma_aE/drho */
/* 			     dSigma_sF/dT,   dSigma_aF/dT,   dSigma_aP/dT,   dSigma_aE/dT */
	
	Real Tpower, Tpower1;
	/* Tpower = T^3.5 , Tpower1 = T^4.5 */
	Tpower = 1.0 / (T * T * T * sqrt(T));
	Tpower1 = Tpower / T;
	
	if(Sigma != NULL){
		Sigma[0] =  kappaes * rho;
		Sigma[1] =  kappaff * rho * rho * Tpower;
		Sigma[2] =  kappaff * rho * rho * Tpower;
		Sigma[3] =  kappaff * rho * rho * Tpower;
	}
	
	if(dSigma != NULL){
		dSigma[0] = kappaes;
                dSigma[1] = 2.0 * rho * kappaff * Tpower;
		dSigma[2] = 2.0 * rho * kappaff * Tpower;
		dSigma[3] = 2.0 * rho * kappaff * Tpower;

                dSigma[4] = 0.0;
                dSigma[5] = -3.5 *kappaff * rho * rho * Tpower1;
		dSigma[6] = -3.5 *kappaff * rho * rho * Tpower1;
                dSigma[7] = -3.5 *kappaff * rho * rho * Tpower1;

		
	}
	
	
	return; 
	
}

#endif

/*------------------------------------------------------------------------------
 * expr_dV2: computes delta(Vy) 
 */

static Real expr_dV2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

/*------------------------------------------------------------------------------
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  return pG->U[k][j][i].M1*
    (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
}

static Real hst_rho_dVy2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  dVy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
  return pG->U[k][j][i].d*dVy*dVy;
}

static Real hst_T(const GridS *pG, const int i, const int j, const int k)
{

	Real pressure;
	pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#if defined(MHD) || defined(RADIATION_MHD)
	pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	pressure *= Gamma - 1.0;

	return pressure/pG->U[k][j][i].d;
}

#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */

/*------------------------------------------------------------------------------
 * MHD history variables
 * hst_Bx, etc.: Net flux, and Maxwell stress, added as history variables
 */

#if defined(MHD) || defined(RADIATION_MHD)
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B1c;
}

static Real hst_By(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B2c;
}

static Real hst_Bz(const GridS *pG, const int i, const int j, const int k)
{
  return pG->U[k][j][i].B3c;
}

static Real hst_BxBy(const GridS *pG, const int i, const int j, const int k)
{
  return -pG->U[k][j][i].B1c*pG->U[k][j][i].B2c;
}

static Real hst_dEw2(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVx,dVy,dVz,dBz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  dBz = pG->U[k][j][i].B3c-(sqrt(15.0/16.0))/(2.0*PI)/sqrt(4.*PI);
  dVx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  dVy = pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1;
  dVz = pG->U[k][j][i].M3/pG->U[k][j][i].d;
  
/*  return (dVx*dVx + dVy*dVy + dVz*dVz + pG->U[k][j][i].B1c*pG->U[k][j][i].B1c
    + pG->U[k][j][i].B2c*pG->U[k][j][i].B2c + dBz*dBz); */
  return (pG->U[k][j][i].B1c*pG->U[k][j][i].B1c
    + pG->U[k][j][i].B2c*pG->U[k][j][i].B2c + dBz*dBz); 
}

static Real hst_dBy(const GridS *pG, const int i, const int j, const int k)
{
  double fkx, fky, fkz; /* Fourier kx, ky */
  double dBy;
  Real x1,x2,x3;

/* Lx,Ly, and Lz are globals */

  fky = 2.0*PI/Ly;
  fkx = -4.0*PI/Lx + qshear*Omega_0*fky*pG->time;
  fkz = 2.0*PI/Lz;

/* compute real part of Fourier mode, for comparison to JGG fig 11 */
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  dBy = 2.0*(pG->U[k][j][i].B2c - (0.2-0.15*Omega_0*pG->time));
  dBy *= cos(fkx*x1 + fky*x2 + fkz*x3);

  return dBy;
}


static Real hst_EB(const GridS *pG, const int i, const int j, const int k)
{

	return 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
}

#endif /* MHD */

