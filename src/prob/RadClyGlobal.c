#include "copyright.h"
/*==============================================================================
 * FILE: RadGlobal.c
 *
 * PURPOSE:  Problem generator for global disk 
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
#ifdef SHEARING_BOX
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real grav_vertical(const Real x1, const Real x2, const Real x3);
#endif

void disk_ir(GridS *pG);
void disk_or(GridS *pG);
void disk_iz(GridS *pG);
void disk_oz(GridS *pG);


/*	Paczynski-Witt potential
*/ 
static Real PseudoNewton(const Real x1, const Real x2, const Real x3);
static Real PseudoNewtonAccR(const Real x1, const Real x2, const Real x3);
/* The keplerian rotation velocity along phi direction under the potential */
Real Vkep(const Real x1, const Real x2, const Real x3);

static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real expr_beta(const GridS *pG, const int i, const int j, const int k);
static Real expr_ME(const GridS *pG, const int i, const int j, const int k);
static Real expr_KE(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_VrVphi(const GridS *pG,const int i,const int j,const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
static Real hst_sigmas(const GridS *pG, const int i, const int j, const int k);
static Real hst_sigmaaP(const GridS *pG, const int i, const int j, const int k);
#endif

static Real hst_gravpot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Lx(const GridS *pG, const int i, const int j, const int k);
static Real hst_Ly(const GridS *pG, const int i, const int j, const int k);
static Real hst_Lz(const GridS *pG, const int i, const int j, const int k);

#if defined(MHD) || defined(RADIATION_MHD)
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BrBphi(const GridS *pG, const int i, const int j, const int k);
static Real hst_dEw2(const GridS *pG, const int i, const int j, const int k);
static Real hst_dBy(const GridS *pG, const int i, const int j, const int k);
static Real hst_EB(const GridS *pG, const int i, const int j, const int k);
#endif /* MHD */
static Real hst_T(const GridS *pG, const int i, const int j, const int k);

static Real hst_P(const GridS *pG, const int i, const int j, const int k);

static Real hst_rho2(const GridS *pG, const int i, const int j, const int k);

static Real hst_T2(const GridS *pG, const int i, const int j, const int k);


/* mass and energy flux at top and bottom */
static Real hst_rhofluxtop(const GridS *pG, const int i, const int j, const int k);
static Real hst_rhofluxbottom(const GridS *pG, const int i, const int j, const int k);



static Real kappaes;
static Real kappaffR;
static Real kappaffP;

static Real B0y;
static Real B0z;

static Real dz = 0.0; /* cell size along z direction */
static Real dx = 0.0;
static Real dy = 0.0;

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
static void Thindiskopacity(const PrimS *pW, Real Sigma[NOPACITY], Real dSigma[4]);
#endif
#ifndef SHEARING_BOX
static Real Omega_0 = 4.749736;
static Real qshear = 1.5;
#endif



#ifdef RESISTIVITY
/* To save the maximum J/rho from last step */
/* This is used as normalization factor */
static Real Jrhomax;
static Real Jmaxtemp = 0.0; /* to calculate the maximum J in current time step */

#endif


static Real ztop, zbtm;

static Real betay;
static Real betaz;

static Real pres;


/* Save the initial Er, Fr and rho at ke */

static Real Tfloor = 0.003;
static Real dfloor = 5.e-6;

/* The side of the mask region Rmas = 2 in unit of r_g */
static Real R0 = 30.0;
static Real Crat=800.0;
static Real R_ideal = 1.0;



/* For transfer module */
#ifdef RADIATION_TRANSFER


static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real Transfereps(const GridS *pG,  const int ifr, const int i, const int j, 
		      const int k);
static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
			  const int k);



#endif

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
/* If the problem generator is called, Domain.Grid must be not NULL */
void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;

#ifdef RADIATION_TRANSFER
  RadGridS *pRG = pDomain->RadGrid;
  int il,iu,jl,ju,kl,ku,ifr,m;
  int ioff,joff,koff,ig,jg,kg;
  int nf=pRG->nf, nang=pRG->nang;
  int niter0, dScnv0;
#endif

  int ixs,jxs,kxs,i,j,k,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax;
  Real Radius, distance, Vc, costheta, sintheta;
  Real den = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real density, pressure, Er, vx, vy, Frx, Fry, Frz;
  Real factor = 1.0;
  Real beta=1.0,B0,kx,ky,kz,amp;
  Real Bamp;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  ifield = par_geti_def("problem","ifield", 6);
  ifield = 6;
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */
  int lines;
  int ID, NGx, NGy, NGz, iproc, jproc, kproc;

/* Parameter used to read in the data */
#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	get_myGridIndex(pDomain, myID_Comm_world, &(iproc), &(jproc), &(kproc));	
#else
	ID = 0;
	
	iproc = 0;
	jproc = 0;
	kproc = 0;
#endif
	NGx = pDomain->NGrid[0];
	NGy = pDomain->NGrid[1];
	NGz = pDomain->NGrid[2];
	
	if((pDomain->Nx[0]) != (pDomain->Nx[1]))
		ath_error("[Problem] The global disk only cover one quarter and requires periodic boundary with respect to is and js!\n");


	




	/* Parse global variables of unit ratio */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	Prat = par_getd("problem","Pratio");
	Crat = par_getd("problem","Cratio");
	R_ideal = par_getd("problem","R_ideal");
	Eratio = par_getd_def("problem","Eratio",0.00);
	Erflag = par_getd_def("problem","Erflag",1);
	Ncycle = par_getd_def("problem","Ncycle",15);
  	TOL  = par_getd_def("problem","TOL",1.e-10);
#endif
	
	
/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
#ifdef SHEARING_BOX
  Omega_0 = par_getd_def("problem","Omega",4.749736);
  qshear  = par_getd_def("problem","qshear",1.5);
#endif

	betaz=0.0;
	betay=100.0;
	pres = 1.0;

  	B0z = 0.0;
  	B0y = sqrt((double)(2.0*pres/betay));
	B0 = sqrt(B0z * B0z + B0y * B0y);	

	kappaes = 6.48511e3;
        kappaffP = 177.58;
        kappaffR = 4.79946;

	/* Initialize temperature unit */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	T0= 1.e7;
#endif


/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs)) + ID;

/* Initialize boxsize */
  ztop = pDomain->RootMaxX[2];
  zbtm = pDomain->RootMinX[2];
  dz = pGrid->dx3;
  dx = pGrid->dx1;
  dy = pGrid->dx2;

  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  kx = 16.0*PI/Lx;

	amp = 0.01;
	Bamp = 1.0;

	/* Parameters to setup the initial torus */

	Real Lprofile = 0.0;
	Real vs0 = 8.0;
	Real rho0 = 5.0;
	Real L0 = sqrt(R0/2.0) * R0 * Crat/(R0 - 1.0);
	Real Effphi, Effphi0, tempphi;
	Real nindex, langular, vphi;
	nindex = 1.0/(Gamma-1.0);


/* Rescale amp to sound speed for ipert 2,3 */


  for (k=ks-nghost; k<=ke+nghost; k++) {
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
       x1 = x1vc(pGrid,i);

	Radius = x1;
	distance = sqrt(x1*x1+x3*x3);

	langular = L0*pow(Radius/R0,Lprofile);
	vphi = langular/Radius;

	vy = vphi;

	Effphi = -Crat*Crat/(2.0*(distance-1.0)) + (pow((langular/Radius),2.0))/(2.0*(1.0-Lprofile));
	Effphi0 = -Crat*Crat/(2.0*(R0-1.0)) + (pow((L0/R0),2.0))/(2.0*(1.0-Lprofile));
	tempphi = ((Effphi-Effphi0)/nindex)/(vs0*vs0);
	if((abs(tempphi) < 1.0) && (distance > 10.0)){
		density = rho0 * pow((1.0-tempphi),nindex);
		/* random perturbation */
		amp = 0.0;
		if(density < dfloor){
			density = dfloor;
		}
		pressure = rho0*vs0*vs0*pow(density/rho0,Gamma)/Gamma;
	}
	else{
		vy = Vkep(x1,x2,x3);
		density = dfloor;
		pressure = density * vs0 * vs0/Gamma;

/*		vx = 0.0;
		vy = 0.0;
		amp = 0.0;

		Effphi = -Crat*Crat/(2.0*(distance-1.0));
		Effphi0 = -Crat*Crat/(2.0*(R0-1.0));

		if(Radius < R0){
			
			tempphi = (Effphi-Effphi0)/(1.e5/Gamma);
			density = 1.e-5 * exp(-tempphi);
			pressure = sqrt(1.e5/Gamma) * density;
		}
		else{
			tempphi = (Effphi-Effphi0)/(1.e4/Gamma);
			density = 1.e-5 * exp(-tempphi);
			pressure = sqrt(1.e4/Gamma) * density;
		}
*/
		
	}


	
	

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
        rp = pressure*(1.0 + 2.0*rval);
        rd = density;
#else
        rd = density*(1.0 + 2.0*rval);
#endif
/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
	
        rvx = 0.0 * 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
	
        rvy = 0.0 * 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
	
        rvz = 0.0 * 0.4*rval*sqrt(pres/den);
      

/* Initialize d, M, and P.  For 3D shearing box M1=Vx, M2=Vy, M3=Vz
 * With FARGO do not initialize the background shear */ 

      pGrid->U[k][j][i].d  = rd;
      pGrid->U[k][j][i].M1 = rd*vx;
      pGrid->U[k][j][i].M2 = rd*vy;

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

	/* Apply a uniform azimuthal magnetic field */
        pGrid->B1i[k][j][i] = 0.0;
      	pGrid->B2i[k][j][i] = 0.0;
      	pGrid->B3i[k][j][i] = 0.0;
      	if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
      	if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
      	if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;

	if (ifield == 2) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = 0.0;
        pGrid->B3i[k][j][i] = B0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = 0.0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0;
      }	
      if (ifield == 3) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0*(cos((double)kx*x1));
        pGrid->B3i[k][j][i] = B0*(sin((double)kx*x1));
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0*(cos((double)kx*x1));
        if (k==ke) pGrid->B3i[ke+1][j][i] = B0*(sin((double)kx*x1));
      }
	/* Only apply magnetic field within a certain range */
	if (ifield == 4 && fabs(x3) < 1.0) {
        pGrid->B1i[k][j][i] = 0.0;
        pGrid->B2i[k][j][i] = B0;
        pGrid->B3i[k][j][i] = 0.0;
        if (i==ie) pGrid->B1i[k][j][ie+1] = 0.0;
        if (j==je) pGrid->B2i[k][je+1][i] = B0;
        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.0;

	}
	if(ifield == 4){

		 pGrid->B3i[k][j][i] = 0.1*B0*(sin((double)kx*x1));
                        if (k==ke) pGrid->B3i[ke+1][j][i] = 0.1*B0*(sin((double)kx*x1));
	}
	
     

#endif /* MHD */
		
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
		pGrid->U[k][j][i].Edd_11 = 1.0/3.0; 
		pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
		pGrid->U[k][j][i].Edd_33 = 1.0/3.0;
		pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
		pGrid->U[k][j][i].Edd_31 = 0.0;
		pGrid->U[k][j][i].Edd_32 = 0.0;
		pGrid->U[k][j][i].Sigma[0] = kappaes * rd;
		pGrid->U[k][j][i].Sigma[1] = kappaffR * rd * rd * pow(rp/(rd*R_ideal),-3.5);
		pGrid->U[k][j][i].Sigma[2] = kappaffP * rd * rd * pow(rp/(R_ideal*rd),-3.5);
		pGrid->U[k][j][i].Sigma[3] = kappaffP * rd * rd * pow(rp/(R_ideal*rd),-3.5);
		
		pGrid->U[k][j][i].Er = Er;
		pGrid->U[k][j][i].Fr1 = ((1.0 + pGrid->U[k][j][i].Edd_11) * vx + pGrid->U[k][j][i].Edd_21 * vy + pGrid->U[k][j][i].Edd_31 * rvz)* pGrid->U[k][j][i].Er / Crat;
		pGrid->U[k][j][i].Fr2 = ((1.0 + pGrid->U[k][j][i].Edd_22) * vy + pGrid->U[k][j][i].Edd_21 * vx + pGrid->U[k][j][i].Edd_32 * rvz)* pGrid->U[k][j][i].Er / Crat;

		pGrid->U[k][j][i].Fr3 = ((1.0 + pGrid->U[k][j][i].Edd_33) * rvz + pGrid->U[k][j][i].Edd_31 * vx + pGrid->U[k][j][i].Edd_32 * vy)* pGrid->U[k][j][i].Er / Crat;
		
		
#endif
		
		
		
		
		
    }
  }}
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
	/* Now add the comoving flux */
	for (k=ks; k<=ke; k++) {
  	for (j=js; j<=je; j++) {
    	for (i=is; i<=ie; i++) {
		Frx = -((pGrid->U[k][j][i+1].Edd_11 * pGrid->U[k][j][i+1].Er - pGrid->U[k][j][i-1].Edd_11 * pGrid->U[k][j][i-1].Er)/(2.0 * pGrid->dx1))/(pGrid->U[k][j][i].Sigma[0]+pGrid->U[k][j][i].Sigma[1]);
		Fry = -((pGrid->U[k][j+1][i].Edd_22 * pGrid->U[k][j+1][i].Er - pGrid->U[k][j-1][i].Edd_22 * pGrid->U[k][j-1][i].Er)/(2.0 * pGrid->dx2))/(pGrid->U[k][j][i].Sigma[0]+pGrid->U[k][j][i].Sigma[1]);
		Frz = -((pGrid->U[k+1][j][i].Edd_33 * pGrid->U[k+1][j][i].Er - pGrid->U[k-1][j][i].Edd_33 * pGrid->U[k-1][j][i].Er)/(2.0 * pGrid->dx3))/(pGrid->U[k][j][i].Sigma[0]+pGrid->U[k][j][i].Sigma[1]);

		pGrid->U[k][j][i].Fr1 += Frx;
		pGrid->U[k][j][i].Fr2 += Fry;
		pGrid->U[k][j][i].Fr3 += Frz;

	}
	}
	}
#endif
#if defined(MHD) || defined(RADIATION_MHD)
	if (ifield == 6) {
  /* flux tube of Hirose et al. We put this down here to break away from the
 *      for loops used above. */
                Real xc=0.0,zc=0.0,Bpratio=0.25,Bp0,rad0,rad,x1h,x3h;
		Real Aphi, Ay0, Ax0;

        Real ***Ay;
	Real ***Ax;
        if((Ay = (Real***)calloc_3d_array(pGrid->Nx[2]+2*nghost,pGrid->Nx[1]+2*nghost, pGrid->Nx[0]+2*nghost ,sizeof(Real))) == NULL)
                ath_error("[problem]: malloc return a NULL pointer\n");

	if((Ax = (Real***)calloc_3d_array(pGrid->Nx[2]+2*nghost,pGrid->Nx[1]+2*nghost, pGrid->Nx[0]+2*nghost ,sizeof(Real))) == NULL)
                ath_error("[problem]: malloc return a NULL pointer\n");


    for (k=ks; k<=ke+2; k++) {
      for (j=js; j<=je+2; j++) {
        for (i=is; i<=ie+2; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	Radius = sqrt(x1*x1+x2*x2);

	if(pGrid->U[k][j][i].d > dfloor)
		Aphi = B0 * pGrid->U[k][j][i].d;
	else
		Aphi = 0.0;

		/* Now convert to Ax and Ay */
		Ax[k][j][i] = Aphi * x1/Radius;
		Ay[k][j][i] = -Aphi * x2/Radius;
	
        }
      }
    }
  /* In this case, we calculate face fields from vect. potential */
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie+1; i++) {
          pGrid->B1i[k][j][i] = -(Ay[k+1][j][i]-Ay[k][j][i])/pGrid->dx3;
	  pGrid->B2i[k][j][i] = (Ax[k+1][j][i]-Ax[k][j][i])/pGrid->dx3;
          pGrid->B3i[k][j][i] =  (Ay[k][j][i+1]-Ay[k][j][i])/pGrid->dx1 - (Ax[k][j+1][i]-Ax[k][j][i])/pGrid->dx2;		
        }
      }
    }
  /* Sync poloidal centered fields to face fields for the next step */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
        }
      }
    }
  for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j+1][i]+pGrid->B2i[k][j][i]);      
      }
    }
 
  /* Finally, calculate cell centered By field from face fields */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
        
#ifdef ADIABATIC
          pGrid->U[k][j][i].E += 0.5*(SQR(pGrid->U[k][j][i].B1c)
               + SQR(pGrid->U[k][j][i].B2c) + SQR(pGrid->U[k][j][i].B3c));
#endif
        }
      }
    }
  } else {
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
  } /* End if ifield = 6 */
#endif /* MHD */

/* Check the magnetic field is divergence free */
#if defined(MHD) || defined(RADIATION_MHD)
  Real divb;
/* Finally, let's check that the field is divergenceless */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        divb = (pGrid->B1i[k][j][i+1]-pGrid->B1i[k][j][i])/pGrid->dx1 +
               (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/pGrid->dx2 +
               (pGrid->B3i[k+1][j][i]-pGrid->B3i[k][j][i])/pGrid->dx3;
        if (fabs(divb) >= 1.e-12) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
		printf("z: %f\n",x3);
          ath_error("[problem]: Nonzero divergence of initial magnetic field\n");
        }
      }
    }
  }

#endif

#ifdef RESISTIVITY
  eta_Ohm = par_getd_def("problem","eta_O",1.0);
  Jrhomax = 1.0;
#endif



/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
  StaticGravPot = grav_vertical;
#else
  StaticGravPot = PseudoNewton;
#endif

  x1GravAcc = PseudoNewtonAccR;
  bvals_mhd_fun(pDomain,left_x1,disk_ir);
  bvals_mhd_fun(pDomain,right_x1,disk_or);
  bvals_mhd_fun(pDomain,left_x3,disk_iz);
  bvals_mhd_fun(pDomain,right_x3,disk_oz);
	
#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)	
	/* enroll the opacity function */
	Opacity = Thindiskopacity;
#endif
	
	
/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_VrVphi, "<rho Vr Vphi>");
    dump_history_enroll(hst_T,"<T>");
    dump_history_enroll(hst_T2,"<T^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif

#if defined(MHD) || defined(RADIATION_MHD)
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BrBphi, "<-Br Bphi>");
    dump_history_enroll(hst_EB,"<B^2/2>");

#endif /* MHD */

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
   dump_history_enroll(hst_sigmas,"<Sigma_es>");
   dump_history_enroll(hst_sigmaaP,"<Sigma_aP>");
#endif

   dump_history_enroll(hst_rho2,"<rho^2>");
   dump_history_enroll(hst_P,"Pg");
   dump_history_enroll(hst_rhofluxtop,"rhov_top");
   dump_history_enroll(hst_rhofluxbottom,"rhov_bottom");
   dump_history_enroll(hst_gravpot,"GravPot");	
   dump_history_enroll(hst_Lx,"Lx");
   dump_history_enroll(hst_Ly,"Ly");
   dump_history_enroll(hst_Lz,"Lz");

    frst = 0;
  }




/* data for radiation transfer method */
#ifdef RADIATION_TRANSFER

  il = pRG->is-1, iu = pRG->ie+1;
  jl = pRG->js-1, ju = pRG->je+1;
  kl = pRG->ks-1, ku = pRG->ke+1;
  ioff = nghost - 1;
  joff = nghost - 1;
  koff = nghost - 1;
/* Initialize mean intensity */
  for(ifr=0; ifr<nf; ifr++) {
    for (k=kl; k<=ku; k++) {
      kg = k + koff;
      if (k == kl) kg++;
      if (k == ku) kg--;
      for (j=jl; j<=ju; j++) {
	jg = j + joff;
	if (j == jl) jg++;
	if (j == ju) jg--;
	for(i=il; i<=iu; i++) {
	  ig = i + ioff;
	  if (i == il) ig++;
	  if (i == iu) ig--;
	  pRG->R[ifr][k][j][i].J = Thermal_B(pGrid,ifr,ig,jg,kg);
	}}}
    
  }
 
/* ------- Initialize boundary emission ---------------------------------- */

/* Density gradient aligned with i3 */
/* For disk, zero intensity any direction */
    for(ifr=0; ifr<nf; ifr++) {
/* Initialize boundary intensity in x1 direction */      
      for(k=kl; k<=ku; k++) {
	kg = k + koff;
	for(j=jl; j<=ju; j++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][0][m] = 0.0;
	    pRG->Ghstr1i[ifr][k][j][2][m] = 0.0;
	    pRG->Ghstr1i[ifr][k][j][4][m] = 0.0;
	    pRG->Ghstr1i[ifr][k][j][6][m] = 0.0;
	    pRG->Ghstl1i[ifr][k][j][1][m] = 0.0;
	    pRG->Ghstl1i[ifr][k][j][3][m] = 0.0;
	    pRG->Ghstl1i[ifr][k][j][5][m] = 0.0;
	    pRG->Ghstl1i[ifr][k][j][7][m] = 0.0;	    
	  }}
/* Initialize boundary intensity in x2 direction */
	for(i=il; i<=iu; i++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl2i[ifr][k][i][0][m] = 0.0;
	    pRG->Ghstl2i[ifr][k][i][1][m] = 0.0;
	    pRG->Ghstl2i[ifr][k][i][4][m] = 0.0;
	    pRG->Ghstl2i[ifr][k][i][5][m] = 0.0;
	    pRG->Ghstr2i[ifr][k][i][2][m] = 0.0;
	    pRG->Ghstr2i[ifr][k][i][3][m] = 0.0;
	    pRG->Ghstr2i[ifr][k][i][6][m] = 0.0;
	    pRG->Ghstr2i[ifr][k][i][7][m] = 0.0;
	  }}
      }  
/* Initialize boundary intensity in x3 direction */
      for(j=jl; j<=ju; j++) {
	for(i=il; i<=iu; i++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][0][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][1][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][2][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][3][m] = 0.0;
	    pRG->Ghstr3i[ifr][j][i][4][m] = 0.0;
	    pRG->Ghstr3i[ifr][j][i][5][m] = 0.0;
	    pRG->Ghstr3i[ifr][j][i][6][m] = 0.0;
	    pRG->Ghstr3i[ifr][j][i][7][m] = 0.0;
	  }}}
    }

/* enroll radiation specification functions */
get_thermal_source = Thermal_B;
get_thermal_fraction = Transfereps;
get_total_opacity = transfer_opacity;

#endif



/* INItialize the Eddington tensor */
/* If initial data depends on Eddington tensor, we need to use the Eddington tensor calculated here */
#ifdef RADIATION_TRANSFER
/*
	hydro_to_rad(pDomain);
  

	formal_solution(pDomain);


	Eddington_FUN(pGrid, pRG);

*/

#endif


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




void problem_write_restart(MeshS *pM, FILE *fp)
{
	
	fwrite(&Gamma,sizeof(Real),1,fp);

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
	fwrite(&Prat,sizeof(Real),1,fp);
	fwrite(&Crat,sizeof(Real),1,fp);
	fwrite(&R_ideal,sizeof(Real),1,fp);
 	fwrite(&kappaes,sizeof(Real),1,fp);
	fwrite(&kappaffP,sizeof(Real),1,fp);
	fwrite(&kappaffR,sizeof(Real),1,fp);
#endif	
	
#ifdef SHEARING_BOX
	fwrite(&qshear,sizeof(Real),1,fp);
	fwrite(&Omega_0,sizeof(Real),1,fp);
#endif
#ifdef RADIATION_MHD
	fwrite(&betay,sizeof(Real),1,fp);
	fwrite(&pres,sizeof(Real),1,fp);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID
        fwrite(&Ncycle,sizeof(Real),1,fp);
        fwrite(&TOL,sizeof(Real),1,fp);
#endif
#endif

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
	fwrite(&Eratio,sizeof(Real),1,fp);
	fwrite(&Erflag,sizeof(int),1,fp);
#endif



	


#ifdef RESISTIVITY
	fwrite(&eta_Ohm,sizeof(Real),1,fp);
        fwrite(&Jrhomax,sizeof(Real),1,fp);
#endif

    fwrite(&zbtm,sizeof(Real),1,fp);
    fwrite(&ztop,sizeof(Real),1,fp);

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
    fwrite(&T0,sizeof(Real),1,fp);
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
	fread(&Crat,sizeof(Real),1,fp);
	fread(&R_ideal,sizeof(Real),1,fp);
	fread(&kappaes,sizeof(Real),1,fp);
	fread(&kappaffP,sizeof(Real),1,fp);
	fread(&kappaffR,sizeof(Real),1,fp);
	
#endif
	
#ifdef SHEARING_BOX
    	fread(&qshear,sizeof(Real),1,fp);
	fread(&Omega_0,sizeof(Real),1,fp);
#endif
	
#ifdef RADIATION_MHD
	fread(&betay,sizeof(Real),1,fp);
	fread(&pres,sizeof(Real),1,fp);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID
        fread(&Ncycle,sizeof(Real),1,fp);
        fread(&TOL,sizeof(Real),1,fp);
#endif
#endif

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
	fread(&Eratio,sizeof(Real),1,fp);
	fread(&Erflag,sizeof(int),1,fp);
#endif
/*
	Ncycle = 30;

	TOL = 1.e-10;
*/



#ifdef RESISTIVITY
	fread(&eta_Ohm,sizeof(Real),1,fp);
	fread(&Jrhomax,sizeof(Real),1,fp);
	eta_Ohm = 1;
#endif

       fread(&zbtm,sizeof(Real),1,fp);
       fread(&ztop,sizeof(Real),1,fp);

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
       fread(&T0,sizeof(Real),1,fp);
#endif
       

/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
  StaticGravPot = grav_vertical;
#else
  StaticGravPot = PseudoNewton;
#endif

  x1GravAcc = PseudoNewtonAccR;
 

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
        /* enroll the opacity function */
        Opacity = Thindiskopacity;

#endif
/* enroll new history variables */

	
    dump_history_enroll(hst_rho_VrVphi, "<rho Vx dVy>");
    dump_history_enroll(hst_T,"<T>");
    dump_history_enroll(hst_T2,"<T^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
	
#if defined(MHD) || defined(RADIATION_MHD)
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BrBphi, "<-Bx By>");
    dump_history_enroll(hst_EB,"<B^2/2>");
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
   dump_history_enroll(hst_sigmas,"<Sigma_es>");
   dump_history_enroll(hst_sigmaaP,"<Sigma_aP>");
#endif

   dump_history_enroll(hst_rho2,"<rho^2>");

   dump_history_enroll(hst_P,"Pg");
   dump_history_enroll(hst_rhofluxtop,"rhov_top");
   dump_history_enroll(hst_rhofluxbottom,"rhov_bottom");
   dump_history_enroll(hst_gravpot,"GravPot");
   dump_history_enroll(hst_Lx,"Lx");
   dump_history_enroll(hst_Ly,"Ly");
   dump_history_enroll(hst_Lz,"Lz");

	

	
	/* Increase the background magnetic field */

	int nl, nd;
	DomainS *pD;
	GridS *pGrid;

	for(nl=0; nl<(pM->NLevels); nl++){
	for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
	if(pM->Domain[nl][nd].Grid != NULL){


	pD= &(pM->Domain[nl][nd]);
	pGrid = pD->Grid;

	dz = pGrid->dx3;
	dx = pGrid->dx1;
	dy = pGrid->dx2;

	/* Enroll boundary condition */

	 bvals_mhd_fun(pD,left_x1,disk_ir);
  	 bvals_mhd_fun(pD,right_x1,disk_or);
         bvals_mhd_fun(pD,left_x3,disk_iz);
  	 bvals_mhd_fun(pD,right_x3,disk_oz);



	int is = pGrid->is, ie = pGrid->ie;
	int js = pGrid->js, je = pGrid->je;
	int ks = pGrid->ks, ke = pGrid->ke;
	
	int i, j, k;
	Real betanew = 10.0;
	Real diffBz = sqrt((double)(2.0*pres/betanew));

	Real pressure, averageP = 0.0;
	Real sumP = 0.0;
	int ierr, ID;
	int count = 0;
	Real x1, x2, x3, velocity_x, velocity_y, velocity_z;
#ifdef RADIATION_MHD
/*	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=is; i<=ie; i++) {
			
				 cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

                                if(fabs(x3) < 0.2){


                                        pGrid->B2i[k][j][i] += diffBz;
                                        if (j==je) pGrid->B2i[k][je+1][i] += diffBz;

                                        pGrid->U[k][j][i].E += 0.5 * (2.0 * pGrid->U[k][j][i].B2c * diffBz + diffBz * diffBz);
                                        pGrid->U[k][j][i].B2c += diffBz;

                                }
	
				

			}
		}
	}
*/
/*				averageP /= count;
*/	
#endif
/*
#ifdef MPI_PARALLEL
		ierr = MPI_Reduce(&averageP,&sumP,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
		ID = myID_Comm_world;
		if(ID == 0){
			sumP /= pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
		}
		MPI_Bcast(&sumP,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

		averageP = sumP;
#endif
*/
	/* Initialize such that radiation pressure / gas pressure =125 */
/*	
	for (k=ks; k<=ke; k++) {
                for (j=js; j<=je; j++) {
	                  for (i=is; i<=ie; i++) {
				cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

	
				velocity_x = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
				velocity_y = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
				velocity_z = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;

				if(pGrid->U[k][j][i].d > 2.0 * dfloor)
					pGrid->U[k][j][i].Er = 1.0;
				else
					pGrid->U[k][j][i].Er = pow(Tfloor, 4.0);

				pGrid->U[k][j][i].Fr1 = ((1.0 + pGrid->U[k][j][i].Edd_11) * velocity_x + pGrid->U[k][j][i].Edd_21 * velocity_y + pGrid->U[k][j][i].Edd_31 * velocity_z)* pGrid->U[k][j][i].Er / Crat;
				pGrid->U[k][j][i].Fr2 = ((1.0 + pGrid->U[k][j][i].Edd_22) * velocity_y + pGrid->U[k][j][i].Edd_21 * velocity_x + pGrid->U[k][j][i].Edd_32 * velocity_z)* pGrid->U[k][j][i].Er / Crat;

				pGrid->U[k][j][i].Fr3 = ((1.0 + pGrid->U[k][j][i].Edd_33) * velocity_z + pGrid->U[k][j][i].Edd_31 * velocity_x + pGrid->U[k][j][i].Edd_32 * velocity_y)* pGrid->U[k][j][i].Er / Crat;


			}
		}
	}

*/

/* define the radiation transfer module */
#ifdef RADIATION_TRANSFER
	
/*	RadGridS *pRG = pD->RadGrid;
	par_seti("radiation","niter0","%d",10,"# of initial iterations");
	par_setd("radiation","dScnv0","%e",5.0e-3,"threshold for initial convergence");
	par_seti("radiation","niter","%d",10,"# of initial iterations");
	par_setd("radiation","dScnv","%e",5.0e-3,"threshold for initial convergence");
*/
/*        pRG->nmu = 4.0;
 *                  pRG->nang =  pRG->nmu * (pRG->nmu + 1) / 2;
 *                            pRG->noct = 8;
 *                                      pRG->nf = 1;
 *                                      */
 /*       int nf=pRG->nf, nang=pRG->nang;
        int il = pRG->is-1, iu = pRG->ie+1;
        int jl = pRG->js-1, ju = pRG->je+1;
        int kl = pRG->ks-1, ku = pRG->ke+1;
        int ioff = nghost - 1;
        int joff = nghost - 1;
        int koff = nghost - 1;
        int ifr, kg, ig, jg;
        Real eps;
*/	
	/* Initialize mean intensity */
/*		for(ifr=0; ifr<nf; ifr++) {
	  for (k=kl; k<=ku; k++) {
	    kg = k + koff;
	    if (k == kl) kg++;
	    if (k == ku) kg--;
	    for (j=jl; j<=ju; j++) {
	      jg = j + joff;
	      if (j == jl) jg++;
	      if (j == ju) jg--;
	      for(i=il; i<=iu; i++) {
		ig = i + ioff;
                if (i == il) ig++;
                if (i == iu) ig--;
                eps = Transfereps(pGrid,ifr,ig,jg,kg);
		 pRG->R[ifr][k][j][i].J = eps * Thermal_B(pGrid,ifr,ig,jg,kg);
	      }}}

	      } 
*/

	/* Density gradient aligned with i3 */
/* For disk, zero intensity any direction */
 /*   for(ifr=0; ifr<nf; ifr++) {
      for(k=kl; k<=ku; k++) {
        kg = k + koff;
        for(j=jl; j<=ju; j++) {
          for(m=0; m<nang; m++) {
            pRG->Ghstr1i[ifr][k][j][0][m] = 0.0;
            pRG->Ghstr1i[ifr][k][j][2][m] = 0.0;
            pRG->Ghstr1i[ifr][k][j][4][m] = 0.0;
            pRG->Ghstr1i[ifr][k][j][6][m] = 0.0;
            pRG->Ghstl1i[ifr][k][j][1][m] = 0.0;
            pRG->Ghstl1i[ifr][k][j][3][m] = 0.0;
            pRG->Ghstl1i[ifr][k][j][5][m] = 0.0;
            pRG->Ghstl1i[ifr][k][j][7][m] = 0.0;
          }}
        for(i=il; i<=iu; i++) {
          for(m=0; m<nang; m++) {
            pRG->Ghstl2i[ifr][k][i][0][m] = 0.0;
            pRG->Ghstl2i[ifr][k][i][1][m] = 0.0;
            pRG->Ghstl2i[ifr][k][i][4][m] = 0.0;
            pRG->Ghstl2i[ifr][k][i][5][m] = 0.0;
            pRG->Ghstr2i[ifr][k][i][2][m] = 0.0;
            pRG->Ghstr2i[ifr][k][i][3][m] = 0.0;
            pRG->Ghstr2i[ifr][k][i][6][m] = 0.0;
            pRG->Ghstr2i[ifr][k][i][7][m] = 0.0;
          }}
      }

      for(j=jl; j<=ju; j++) {
        for(i=il; i<=iu; i++) {
          for(m=0; m<nang; m++) {
            pRG->Ghstl3i[ifr][j][i][0][m] = 0.0;
            pRG->Ghstl3i[ifr][j][i][1][m] = 0.0;
            pRG->Ghstl3i[ifr][j][i][2][m] = 0.0;
            pRG->Ghstl3i[ifr][j][i][3][m] = 0.0;
            pRG->Ghstr3i[ifr][j][i][4][m] = 0.0;
            pRG->Ghstr3i[ifr][j][i][5][m] = 0.0;
            pRG->Ghstr3i[ifr][j][i][6][m] = 0.0;
            pRG->Ghstr3i[ifr][j][i][7][m] = 0.0;
          }}}
    }

*/


        get_thermal_source = Thermal_B;
        get_thermal_fraction = Transfereps;
        get_total_opacity = transfer_opacity;


/*      hydro_to_rad(pDomain);
 *      formal_solution(pDomain);
 *      Eddington_FUN(pGrid, pRG);
 */
#endif
	}/* Grid != NULL */
	}/* nd */
	}/* nl */

  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
 return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
 
	return NULL;
}


#ifdef RESISTIVITY
/*void get_eta_user(GridS *pG, int i, int j, int k,
					  Real *eta_O, Real *eta_H, Real *eta_A)
{

	Real dl, eta0, lz;
        Real x1, x2, x3;
        Real pressure, T, Bp;
	dl = pG->dx1;
        if(pG->dx2 < dl) dl = pG->dx2;
        if(pG->dx3 < dl) dl = pG->dx3;

	if(pG->dt > TINY_NUMBER){	
		
        	eta0 = 0.03 * dl * dl/pG->dt;
	}	
	else{
		eta0 = 0.0;
	}
	cc_pos(pG, i, j,k, &x1, &x2, &x3);

	lz = 8.0;
	if(pG->time < 40.0){

		if(fabs(x3) < 2.0)
                        eta0 = 0.0;
	
	

       	*eta_O = 0.5 * eta0 * (sin(0.5 * PI * (1.0 * fabs(x3)  - 3.0 )) + 1.0);
	*eta_O = eta0;
	}
	else{
		*eta_O = 0.0;
	}

	*eta_O = 0.0;
	*eta_H = 0.0;
	*eta_A = 0.0;
		return;
}
*/

void get_eta_user(GridS *pG, int i, int j, int k,
					  Real *eta_O, Real *eta_H, Real *eta_A)
{

	int is = pG->is; 
	int ie = pG->ie;
	int js = pG->js;
	int je = pG->je;
	int ks = pG->ks;
	int ke = pG->ke;
	Real dl, eta0, lz;
	Real jx1, jx2, jx3, jx4, jy1, jy2, jy3, jy4, jz1, jz2, jz3, jz4, jx, jy, jz, J;
        Real x1, x2, x3, distance, Radius;
	dl = pG->dx1;
        if(pG->dx2 < dl) dl = pG->dx2;
        if(pG->dx3 < dl) dl = pG->dx3;

	if(pG->dt > TINY_NUMBER){	
		
        	eta0 = 0.03 * dl * dl/pG->dt;
	}	
	else{
		eta0 = 0.0;
	}
	Real factor = 0.01;
	/* First, calculate the current */
	/* current is on the edge */
	if((i != (is-4)) && (j != (js-4)) && (k != (ks-4)) && (i != (ie+4)) && (j != (je+4)) && (k != (ke+4))){
		/* use cell centered magnetic field */
		/* To get cell centered current */
		/* x component */
		jx1 = (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ])/pG->dx2 -
                        (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ])/pG->dx3;

		jx2 = (pG->B3i[k][j+1][i] - pG->B3i[k  ][j][i  ])/pG->dx2 -
                        (pG->B2i[k][j+1][i] - pG->B2i[k-1][j+1 ][i  ])/pG->dx3;

		jx3 = (pG->B3i[k+1][j][i] - pG->B3i[k+1][j-1][i  ])/pG->dx2 -
                        (pG->B2i[k+1][j][i] - pG->B2i[k][j  ][i  ])/pG->dx3;

		jx4 = (pG->B3i[k+1][j+1][i] - pG->B3i[k+1][j][i  ])/pG->dx2 -
                        (pG->B2i[k+1][j+1][i] - pG->B2i[k][j+1][i  ])/pG->dx3;

		jx = 0.25 * (jx1 + jx2 + jx3 + jx4);

		/* y component */
		jy1 = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])/pG->dx3 -
                        (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])/pG->dx1;

		jy2 = (pG->B1i[k][j][i+1] - pG->B1i[k-1][j  ][i+1])/pG->dx3 -
                        (pG->B3i[k][j][i+1] - pG->B3i[k  ][j  ][i])/pG->dx1;

		jy3 = (pG->B1i[k+1][j][i] - pG->B1i[k][j  ][i  ])/pG->dx3 -
                        (pG->B3i[k+1][j][i] - pG->B3i[k+1][j  ][i-1])/pG->dx1;

		jy4 = (pG->B1i[k+1][j][i+1] - pG->B1i[k][j  ][i+1])/pG->dx3 -
                        (pG->B3i[k+1][j][i+1] - pG->B3i[k+1][j ][i])/pG->dx1;

		jy = 0.25 * (jy1 + jy2 + jy3 + jy4);
		/* z component */
		jz1 = (pG->B2i[k][j][i] - pG->B2i[k  ][j  ][i-1])/pG->dx1 -
                        (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])/pG->dx2;

		jz2 = (pG->B2i[k][j][i+1] - pG->B2i[k  ][j  ][i])/pG->dx1 -
                        (pG->B1i[k][j][i+1] - pG->B1i[k  ][j-1][i+1  ])/pG->dx2;

		jz3 = (pG->B2i[k][j+1][i] - pG->B2i[k  ][j+1][i-1])/pG->dx1 -
                        (pG->B1i[k][j+1][i] - pG->B1i[k  ][j][i  ])/pG->dx2;

		jz4 = (pG->B2i[k][j+1][i+1] - pG->B2i[k  ][j+1][i])/pG->dx1 -
                        (pG->B1i[k][j+1][i+1] - pG->B1i[k  ][j][i+1])/pG->dx2;

		jz = 0.25 * (jz1 + jz2 + jz3 + jz4);
		
		J = sqrt(jx * jx + jy * jy + jz * jz);
	}
	else{
		J = 0.0;
	}

	J /= pG->U[k][j][i].d;

	if(J > Jmaxtemp)
		Jmaxtemp = J;

	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	
	
	if(x1 < 10.0){
		if(Jrhomax > 0.0){
			if((J > factor * Jrhomax) && (J > 0.01)){
				*eta_O = eta0 * J * J /( Jrhomax * Jrhomax);
			}
			else{
				*eta_O = 0.0;
			}
		}
		else{
			*eta_O = 0.0;
		}
	}
	else{
		*eta_O = 0.0;
	}

	if(*eta_O > eta0)
		*eta_O = eta0;

	

	*eta_H = 0.0;
	*eta_A = 0.0;

	return;
}
#endif




void disk_ir(GridS *pGrid) {
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
	Real Vkep,R,p,z, RBr, Lper;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	cc_pos(pGrid,is,j,k,&R,&p,&z);
#ifdef MHD
    	RBr = (R-0.5*pGrid->dx1)*pGrid->B1i[k][j][is];
#endif
    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = R*pGrid->U[k][j][is].M2;
#else
    	Vkep = sqrt(R/2.0)*Crat/(R-1.0);
	/* residual angular momentum with respect to the Keplerian angular momentum */
    	Lper = R*pGrid->U[k][j][is].M2 - R*pGrid->U[k][j][is].d*Vkep;
#endif

      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
				
	/* Calculate Keplerian velocity */
	cc_pos(pGrid,is-i,j,k,&R,&p,&z);
	Vkep = sqrt(R/2.0)*Crat/(R-1.0);
#ifdef FARGO
				
	pGrid->U[k][j][is-i].M2 = Lper/R;
				
#else
				
	pGrid->U[k][j][is-i].M2 = Lper/R + pGrid->U[k][j][is-i].d*Vkep;
				
#endif

      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	cc_pos(pGrid,is,j,k,&R,&p,&z);
    	RBr = (R-0.5*pGrid->dx1)*pGrid->B1i[k][j][is];
      for (i=1; i<=nghost-1; i++) {
      	cc_pos(pGrid,is,j,k,&R,&p,&z);
      	
      	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      	
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
  
      	pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      	
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      		pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      	

      }
    }
  }
#endif /* MHD */

  return;


}

void disk_or(GridS *pGrid) {

	int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
	Real Vkep, R,p,z, RBr, Lper;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
#ifdef MHD
    	RBr = pGrid->ri[ie+1]*pGrid->B1i[k][j][ie+1];
#endif
    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = pGrid->r[ie]*pGrid->U[k][j][ie].M2;
#else
	cc_pos(pGrid,ie,j,k,&R,&p,&z);
	Vkep = sqrt(R/2.0)*Crat/(R-1.0);
	Lper = R*pGrid->U[k][j][ie].M2 - R*pGrid->U[k][j][ie].d*Vkep;
#endif
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
				
	cc_pos(pGrid,ie+i,j,k,&R,&p,&z);
	Vkep = sqrt(R/2.0)*Crat/(R-1.0);
#ifdef FARGO
	
	pGrid->U[k][j][ie+i].M2 = Lper/R;
	
#else
	
	pGrid->U[k][j][ie+i].M2 = Lper/R + pGrid->U[k][j][ie+i].d*Vkep;
	
#endif

      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	RBr = pGrid->ri[ie+1]*pGrid->B1i[k][j][ie+1];
      for (i=2; i<=nghost; i++) {
	cc_pos(pGrid,ie+i,j,k,&R,&p,&z);

      	pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];

      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {

      	pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];

      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      		pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;

}


/* Vertical boundary */
void disk_oz(GridS *pGrid)
{
  	int i, je,j, js, m;
	int ks, is,ie, k, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real x1, x2, x3, pressure, dx1, dx2,dx3, velocity1, velocity2, velocity3, Bx, By, Bz, B, vx, vy, vz, T, density;
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	Real Sigma[NOPACITY];
#endif
	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;

	static Real x3t;
	x3t = ztop - 0.5 * pGrid->dx3;
	PrimS Wopacity;


#if defined(MHD) || defined(RADIATION_MHD)
	
	for(k=1;k<=nghost;k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=is-nghost;i<=ie+nghost;i++){
	

			if(pGrid->U[ke][j][i].M3 < 0.0){
                         	pGrid->B1i[ke+k][j][i] = 0.0 * pGrid->B1i[ke+k-1][j][i];
                                pGrid->B2i[ke+k][j][i] = 0.0 * pGrid->B2i[ke+k-1][j][i];
                        }
                        else{
                                pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
                                pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];
                        }


				if(k<nghost)
					pGrid->B3i[ke+k+1][j][i] = pGrid->B3i[ke+k][j][i];


			}
		}
	}

	/* Now update the cell centered values */
	for(k=1;k<=nghost;k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=is-nghost;i<=ie+nghost;i++){
				 if(pGrid->U[ke][j][i].M3 < 0.0){
                                         pGrid->U[ke+k][j][i].B1c = 0.0 * pGrid->U[ke][j][i].B1c;
                                         pGrid->U[ke+k][j][i].B2c = 0.0 * pGrid->U[ke][j][i].B2c;
                                 }
                                 else{
                                         pGrid->U[ke+k][j][i].B1c = pGrid->U[ke][j][i].B1c;
                                         pGrid->U[ke+k][j][i].B2c = pGrid->U[ke][j][i].B2c;
                                 }
				
				if(k<nghost){
					pGrid->U[ke+k][j][i].B3c = 0.5 * (pGrid->B3i[ke+k][j][i] + pGrid->B3i[ke+k+1][j][i]);
				}
				else{
					pGrid->U[ke+k][j][i].B3c = pGrid->B3i[ke+k][j][i];
				}

			}
		}
	}
		



#endif /* MHD */
	
	for (k=1;  k<=nghost;  k++) {
 		for(j=js-nghost;j<=je+nghost;j++){
    			for(i=is-nghost; i<=ie+nghost; i++){
	    
	

		cc_pos(pGrid, i, j,ke+k, &x1, &x2, &x3);
#if defined(MHD) || defined(RADIATION_MHD)
		Bx = pGrid->U[ke+k][j][i].B1c;
		By = pGrid->U[ke+k][j][i].B2c;
		Bz = pGrid->U[ke+k][j][i].B3c;
		
		B = sqrt(Bx * Bx + By * By + Bz * Bz);
#endif

		if(pGrid->U[ke][j][i].d < TINY_NUMBER)
			pGrid->U[ke][j][i].d = pGrid->U[ke][4][4].d;

	
		density  = pGrid->U[ke][j][i].d;

		
		/*
		pGrid->U[ke+k][j][i].d  = ghostrho[k-1] + pGrid->U[ke][j][i].d - inidata[0];
		*/
	
		
		/* boundary is optical thin. Gas temperature is independent of radiation temperature */

		velocity1 = pGrid->U[ke][j][i].M1 / pGrid->U[ke][j][i].d;
		velocity2 = pGrid->U[ke][j][i].M2 / pGrid->U[ke][j][i].d;
		velocity3 = pGrid->U[ke][j][i].M3 / pGrid->U[ke][j][i].d;

		pressure = (pGrid->U[ke][j][i].E - 0.5 * density * (velocity1 * velocity1 + velocity2 * velocity2 + velocity3 * velocity3)) * (Gamma - 1.0);
#if defined(MHD) || defined(RADIATION_MHD)

		pressure -= 0.5 * (pGrid->U[ke][j][i].B1c * pGrid->U[ke][j][i].B1c + pGrid->U[ke][j][i].B2c * pGrid->U[ke][j][i].B2c + pGrid->U[ke][j][i].B3c * pGrid->U[ke][j][i].B3c) * (Gamma - 1.0);
#endif

		T = pressure / (pGrid->U[ke][j][i].d * R_ideal);
/*
		if(T > Tfloor)
			T = Tfloor;
*/
		/* Tfloor is usually set as initial temperature in the ghost zoner */


		if(velocity3 < TINY_NUMBER){
			velocity3 = 0.0;
		}
		
		/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
		
/*        	pGrid->U[ke+k][j][i].d = MAX(pGrid->U[ke][j][i].d*exp(-(x3*x3-x3t*x3t)/(2.0*T/(Omega_0*Omega_0))), 1.e-8);
		if(pGrid->U[ke+k][j][i].d > pGrid->U[ke][j][i].d) pGrid->U[ke+k][j][i].d = pGrid->U[ke][j][i].d;
*/		
		pGrid->U[ke+k][j][i].d = pGrid->U[ke][j][i].d;

	
		/* set density upper limit in ghost zone */
	/*	if(pGrid->U[ke+k][j][i].d > dfloor)
			pGrid->U[ke+k][j][i].d = dfloor;				
	*/	

		
      		pGrid->U[ke+k][j][i].M1 = velocity1 * pGrid->U[ke+k][j][i].d;
		pGrid->U[ke+k][j][i].M2 = velocity2 * pGrid->U[ke+k][j][i].d;
			
		
		pGrid->U[ke+k][j][i].M3 = velocity3 * pGrid->U[ke+k][j][i].d;

		pressure = T * R_ideal * pGrid->U[ke+k][j][i].d;

		pGrid->U[ke+k][j][i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[ke+k][j][i].M1 * pGrid->U[ke+k][j][i].M1 + pGrid->U[ke+k][j][i].M2 * pGrid->U[ke+k][j][i].M2 + pGrid->U[ke+k][j][i].M3 * pGrid->U[ke+k][j][i].M3) / pGrid->U[ke+k][j][i].d;



#if defined(MHD) || defined(RADIATION_MHD)
		pGrid->U[ke+k][j][i].E += 0.5 * (pGrid->U[ke+k][j][i].B1c * pGrid->U[ke+k][j][i].B1c + pGrid->U[ke+k][j][i].B2c * pGrid->U[ke+k][j][i].B2c + pGrid->U[ke+k][j][i].B3c * pGrid->U[ke+k][j][i].B3c);

#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)

		Wopacity = Cons_to_Prim(&pGrid->U[ke+k][j][i]);
		/* Add background shearing */
#ifdef FARGO	
		cc_pos(pGrid,i,j,ke+k,&x1,&x2,&x3);
		Wopacity.V2 -= qshear * Omega_0 * x1;		
#endif
			
		Thindiskopacity(&Wopacity, Sigma, NULL);
		
		if(T < 2.0 * TINY_NUMBER)
			Sigma[0] = kappaes * pGrid->U[ke+k][j][i].d;

		for(m=0; m<NOPACITY;m++){
			pGrid->U[ke+k][j][i].Sigma[m] = Sigma[m];
		}	
	
#endif	
      		}
    		}
	}
  
}


void disk_iz(GridS *pGrid)
{
  	int i, je,j, ju, jl, js, m;
	int ks, is,ie, k, ke, ku;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real x1, x2, x3, pressure, dx1, dx2,dx3, velocity1, velocity2, velocity3, Bx, By, Bz, B, vx, vy, vz, T, density;
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	Real Sigma[NOPACITY];
#endif
	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;
	PrimS Wopacity;

	
	 static Real x3b;

  	x3b = zbtm+0.5*pGrid->dx3;



#if defined(MHD) || defined(RADIATION_MHD)
	for(k=1;k<=nghost;k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=is-nghost;i<=ie+nghost;i++){

				/* Just copy the magnetic field */
                          /*      if(fabs(pGrid->B1i[ks-k+2][j][i]) > fabs(pGrid->B1i[ks-k+1][j][i]) ){
                                         pGrid->B1i[ks-k][j][i] = 2.0 * pGrid->B1i[ks-k+1][j][i] - pGrid->B1i[ks-k+2][j][i];
                                        if(fabs(pGrid->B1i[ks-k][j][i]) > fabs(pGrid->B1i[ks-k+1][j][i])){
                                                pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks-k+1][j][i];
                                        }
                                }
                                else{
                                        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks-k+1][j][i];
                                }

                                if(fabs(pGrid->B2i[ks-k+2][j][i]) > fabs(pGrid->B2i[ks-k+1][j][i]) ){
                                         pGrid->B2i[ks-k][j][i] = 2.0 * pGrid->B2i[ks-k+1][j][i] - pGrid->B2i[ks-k+2][j][i];
                                        if(fabs(pGrid->B2i[ks-k][j][i]) > fabs(pGrid->B2i[ks-k+1][j][i])){
                                                pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks-k+1][j][i];
                                        }

                                }
                                else{
                                         pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks-k+1][j][i];
                                }
			*/
				if(pGrid->U[ks][j][i].M3 > 0.0){
                                        pGrid->B1i[ks-k][j][i] = 0.0 * pGrid->B1i[ks-k+1][j][i];
                                        pGrid->B2i[ks-k][j][i] = 0.0 * pGrid->B2i[ks-k+1][j][i];
                                }
                                else{
                                        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks-k+1][j][i];
                                        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks-k+1][j][i];
                                }
				
				pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks-k+1][j][i];	

			}
		}
	}
/*
		for(k=1;k<=nghost;k++){
                for(j=js-nghost;j<=je+nghost;j++){
                        for(i=is-nghost;i<=ie+nghost;i++){
				if(i< ie+nghost && j< je+nghost){
                                    pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks-k+1][j][i] + pGrid->dx3 * ((pGrid->B1i[ks-k][j][i+1] - pGrid->B1i[ks-k][j][i]) / pGrid->dx1
                                                                        + (pGrid->B2i[ks-k][j+1][i] - pGrid->B2i[ks-k][j][i]) / pGrid->dx2);
                                }
                                else{
                                        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks-k+1][j][i];
                                }


			}
		}
		}
*/
	/* Now update the cell centered values */
	for(k=1;k<=nghost;k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=is-nghost;i<=ie+nghost;i++){			

				if(pGrid->U[ks][j][i].M3 > 0.0){
					pGrid->U[ks-k][j][i].B1c = 0.0 * pGrid->U[ks][j][i].B1c;
                                        pGrid->U[ks-k][j][i].B2c = 0.0 * pGrid->U[ks][j][i].B2c;
				}
				else{
					pGrid->U[ks-k][j][i].B1c = pGrid->U[ks][j][i].B1c;
					pGrid->U[ks-k][j][i].B2c = pGrid->U[ks][j][i].B2c;
				}
				pGrid->U[ks-k][j][i].B3c = 0.5 * (pGrid->B3i[ks-k][j][i] + pGrid->B3i[ks-k+1][j][i]);

			}
		}
	}
		



#endif /* MHD */
	

 for(j=js-nghost;j<=je+nghost;j++)
    for(i=is-nghost; i<=ie+nghost; i++){
	    for (k=1;  k<=nghost;  k++) {
	

		cc_pos(pGrid, i, j,ks-k, &x1, &x2, &x3);
#if defined(MHD) || defined(RADIATION_MHD)
		Bx = pGrid->U[ks-k][j][i].B1c;
		By = pGrid->U[ks-k][j][i].B2c;
		Bz = pGrid->U[ks-k][j][i].B3c;
		
		B = sqrt(Bx * Bx + By * By + Bz * Bz);
#endif

		if(pGrid->U[ks][j][i].d < TINY_NUMBER)
			pGrid->U[ks][j][i].d = pGrid->U[ks][4][4].d;

		density  = pGrid->U[ks][j][i].d;

	/*	pGrid->U[ks-k][j][i].d  = ghostrho[k-1] + pGrid->U[ks][j][i].d - inidata[0];
	*/
		

		velocity1 = pGrid->U[ks][j][i].M1 / pGrid->U[ks][j][i].d;
		velocity2 = pGrid->U[ks][j][i].M2 / pGrid->U[ks][j][i].d;
		velocity3 = pGrid->U[ks][j][i].M3 / pGrid->U[ks][j][i].d;

		pressure = (pGrid->U[ks][j][i].E - 0.5 * density * (velocity1 * velocity1 + velocity2 * velocity2 + velocity3 * velocity3)) * (Gamma - 1.0);
#if defined(MHD) || defined(RADIATION_MHD)

		pressure -= 0.5 * (pGrid->U[ks][j][i].B1c * pGrid->U[ks][j][i].B1c + pGrid->U[ks][j][i].B2c * pGrid->U[ks][j][i].B2c + pGrid->U[ks][j][i].B3c * pGrid->U[ks][j][i].B3c) * (Gamma - 1.0);
#endif

		T = pressure / (pGrid->U[ks][j][i].d * R_ideal);
/*
		if(T > Tfloor)
			T = Tfloor;
*/			

		if(velocity3 > TINY_NUMBER){
			velocity3 = 0.0;
		}
		
		/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
  /*      	pGrid->U[ks-k][j][i].d = MAX(pGrid->U[ks][j][i].d*exp(-(x3*x3-x3b*x3b)/(2.0*T/(Omega_0*Omega_0))), 1.e-8);
		if(pGrid->U[ks-k][j][i].d > pGrid->U[ks][j][i].d) pGrid->U[ks-k][j][i].d = pGrid->U[ks][j][i].d;
*/
		pGrid->U[ks-k][j][i].d = pGrid->U[ks][j][i].d;

/*		if(pGrid->U[ks-k][j][i].d > dfloor)
			pGrid->U[ks-k][j][i].d = dfloor;
*/
		
      		pGrid->U[ks-k][j][i].M1 = velocity1 * pGrid->U[ks-k][j][i].d;
		pGrid->U[ks-k][j][i].M2 = velocity2 * pGrid->U[ks-k][j][i].d;
			
		
		pGrid->U[ks-k][j][i].M3 = velocity3 * pGrid->U[ks-k][j][i].d;
		
		pressure = T * R_ideal * pGrid->U[ks-k][j][i].d;

		pGrid->U[ks-k][j][i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[ks-k][j][i].M1 * pGrid->U[ks-k][j][i].M1 + pGrid->U[ks-k][j][i].M2 * pGrid->U[ks-k][j][i].M2 + pGrid->U[ks-k][j][i].M3 * pGrid->U[ks-k][j][i].M3) / pGrid->U[ks-k][j][i].d;



#if defined(MHD) || defined(RADIATION_MHD)
		pGrid->U[ks-k][j][i].E += 0.5 * (pGrid->U[ks-k][j][i].B1c * pGrid->U[ks-k][j][i].B1c + pGrid->U[ks-k][j][i].B2c * pGrid->U[ks-k][j][i].B2c + pGrid->U[ks-k][j][i].B3c * pGrid->U[ks-k][j][i].B3c);

#endif


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
		Wopacity = Cons_to_Prim(&pGrid->U[ks-k][j][i]);
		/* Add background shearing */
#ifdef FARGO	
		cc_pos(pGrid,i,j,ks-k,&x1,&x2,&x3);
		Wopacity.V2 -= qshear * Omega_0 * x1;		
#endif
			
		Thindiskopacity(&Wopacity, Sigma, NULL);
		

		if(T < 2.0 * TINY_NUMBER)
			Sigma[0] = kappaes * pGrid->U[ks-k][j][i].d;
		
		for(m=0; m<NOPACITY;m++){
			pGrid->U[ks-k][j][i].Sigma[m] = Sigma[m];
		}	
#endif	
      		
      }
    }
  
  
}



void Userwork_in_loop(MeshS *pM)
{

#ifdef RESISTIVITY
	int ierr;
	Real Jtemp;
	/* replace maximum J/rho with the maximum value from last time step */
	/* Find the maximum drift velocity over the whole mesh for current time step*/
#ifdef MPI_PARALLEL
	ierr = MPI_Allreduce(&Jmaxtemp,&Jtemp,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
	Jmaxtemp = Jtemp;
#endif

	Jrhomax = Jmaxtemp;
	Jmaxtemp = 0.0;
	/* Jmaxtemp only calculate maximum drift velocity for curent step */
#endif





       Real betafloor = 0.000;
	

	 GridS *pG;
        int i, j, k;
        int ie, is;
        int je, js;
        int ke, ks;
	int badcellflag;
        Real pressure, velocity, velocity_x, velocity_y, velocity_z, Bpre, beta;
	Real temprho, tempT, T0;
	PrimS Wtemp;
	Real Vmax = 0.9 * Crat;
	
        
	Real x1, x2, x3, distance, distance1;
	int m, n, l, count;
#if defined(MHD) || defined(RADIATION_MHD)
	Real Vmag;
	/* Limit the Alfven velocity in the mask region so that the time step will not be limited */
#endif


	int nl, nd;


	for(nl=0; nl<(pM->NLevels); nl++){
	for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
	if(pM->Domain[nl][nd].Grid != NULL){



	pG = pM->Domain[nl][nd].Grid;

        ie = pG->ie;
        is = pG->is;
        je = pG->je;
        js = pG->js;
        ke = pG->ke;
        ks = pG->ks;






        for(k=ks; k<=ke; k++) {
                for (j=js; j<=je; j++) {
                        for (i=is; i<=ie; i++) {



   			badcellflag = 0;

			velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
                         velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                         velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;


			/* Limit the maximum velocity */
			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);				
			if(velocity > Vmax){
				velocity_x *= Vmax/velocity;
				velocity_y *= Vmax/velocity;
				velocity_z *= Vmax/velocity;
				pG->U[k][j][i].M1 = pG->U[k][j][i].d * velocity_x;
				pG->U[k][j][i].M2 = pG->U[k][j][i].d * velocity_y;
				pG->U[k][j][i].M3 = pG->U[k][j][i].d * velocity_z;
			}


                         Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
#if defined(MHD) || defined(RADIATION_MHD)
                         Bpre = 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
			 if(Wtemp.P < 2.0 * TINY_NUMBER){
				if(pG->U[k][j][i].d < dfloor){
					 Wtemp.P = dfloor * R_ideal * Tfloor;
				}
				else{
					Wtemp.P = pG->U[k][j][i].d * R_ideal * Tfloor;
				}
				badcellflag = 1;
			}

			
			if(pG->U[k][j][i].d < dfloor){

                              pG->U[k][j][i].d = dfloor;

                              Wtemp.d =  dfloor;
			      pG->U[k][j][i].M1 =  pG->U[k][j][i].d * velocity_x;
			      pG->U[k][j][i].M2 =  pG->U[k][j][i].d * velocity_y;
			      pG->U[k][j][i].M3 =  pG->U[k][j][i].d * velocity_z;

			     badcellflag = 1;
			}

			if(badcellflag){
				 pG->U[k][j][i].E = Wtemp.P / (Gamma - 1.0) + 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
                                 pG->U[k][j][i].E += Bpre;
#endif


			}
#ifdef RADIATION_MHD
                                 if(fabs(Bpre) > 0.0){
                                        beta = Wtemp.P / Bpre;
                                        if(beta < betafloor){
                                                Wtemp.P *= betafloor / beta;
                                                pG->U[k][j][i].E = Wtemp.P / (Gamma - 1.0)  + 0.5 * pG->U[k][j][i].d * (velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
                                                pG->U[k][j][i].E += Bpre;

 					}/* beta < betafloor */

                                }/* Bpre > 0 */
#endif

                        } /* i */
                }/* j */
        }/* k */

	}/* Grid != NULL */
	}/* nd */
	} /* nl */
}

void Userwork_after_loop(MeshS *pM)
{
}
void Userwork_after_first_formal_solution(DomainS *pD)
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
#ifdef SHEARING_BOX
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

static Real grav_vertical(const Real x1, const Real x2, const Real x3)
{
	/* In the ghost zeons, set the gravitational potential to be flat */
	Real h;

	if(x3 > (ztop - 0.5 * dz))
		h = ztop - 0.5 * dz;
	else if(x3 < (zbtm + 0.5 * dz))
		h = zbtm + 0.5 * dz;
	else
		h = x3;		

  return 0.5 * Omega_0 * Omega_0 * h * h;
}

#endif


/* Paczynski - Witta potential */
static Real PseudoNewton(const Real x1, const Real x2, const Real x3)
{
	Real GM, distance, potential;
	/* The dimensionless number */
	GM = Crat*Crat*0.5;
	
	distance = sqrt(x1*x1+x3*x3);
	potential = -GM/(distance - 1.0);
	
	
	return potential;
}


/* Paczynski - Witta potential */
static Real PseudoNewtonAccR(const Real x1, const Real x2, const Real x3)
{
	Real GM, acc, distance;
	/* The dimensionless number */
	GM = Crat*Crat*0.5;
	distance = sqrt(x1*x1+x3*x3);
	
	acc = (x1*GM)/(distance*SQR(distance-1.0));
	
	return acc;
}


/* Paczynski - Witta potential */
Real Vkep(const Real x1, const Real x2, const Real x3)
{
	Real V;

	V = sqrt(x1/2.0)*Crat/(x1-1.0);

	return V;
}





static Real hst_gravpot(const GridS *pG,const int i,const int j, const int k)
{
	Real x1, x2, x3;
	cc_pos(pG,i,j,k,&x1, &x2, &x3);
#ifdef SHEARING_BOX
  return pG->U[k][j][i].d * grav_vertical(x1,x2,x3);
#else
  return pG->U[k][j][i].d * PseudoNewton(x1,x2,x3);
#endif
}


/* angular momentum is rho r X v
 */
static Real hst_Lx(const GridS *pG,const int i,const int j, const int k)
{
	Real x1, x2, x3, vx, vy, vz, rho, Lx;
	
	cc_pos(pG,i,j,k,&x1, &x2, &x3);
	rho = pG->U[k][j][i].d;
	vx = pG->U[k][j][i].M1 / rho;
	vy = pG->U[k][j][i].M2 / rho;
	vz = pG->U[k][j][i].M3 / rho;

	Lx = rho * (x2 * vz - x3 * vy); 

	return Lx;
}


/* angular momentum is rho r X v
 */
static Real hst_Ly(const GridS *pG,const int i,const int j, const int k)
{
	Real x1, x2, x3, vx, vy, vz, rho, Ly;
	
	cc_pos(pG,i,j,k,&x1, &x2, &x3);
	rho = pG->U[k][j][i].d;
	vx = pG->U[k][j][i].M1 / rho;
	vy = pG->U[k][j][i].M2 / rho;
	vz = pG->U[k][j][i].M3 / rho;

	Ly = rho * (x3 * vx - x1 * vz); 

	return Ly;
}


/* angular momentum is rho r X v
 */
static Real hst_Lz(const GridS *pG,const int i,const int j, const int k)
{
	Real x1, x2, x3, vx, vy, vz, rho, Lz;
	
	cc_pos(pG,i,j,k,&x1, &x2, &x3);
	rho = pG->U[k][j][i].d;
	vx = pG->U[k][j][i].M1 / rho;
	vy = pG->U[k][j][i].M2 / rho;
	vz = pG->U[k][j][i].M3 / rho;

	Lz = rho * (x1 * vy - x2 * vx); 

	return Lz;
}



#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
void Thindiskopacity(const PrimS *pW, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY])
{
/* Sigma[0-NOPACITY] are: Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE respectively */
/* dSigma[0-2*NOPACITY] are: dSigma_sF/drho, dSigma_aF/drho, dSigma_aP/drho, dSigma_aE/drho */
/* 			     dSigma_sF/dT,   dSigma_aF/dT,   dSigma_aP/dT,   dSigma_aE/dT */
	
/* When pressure becomes negative, we do not include radiation source term */
Real T, rho;

if((pW->P > TINY_NUMBER) && (pW->d > TINY_NUMBER)){	
	Real Tpower, Tpower1;
	rho = pW->d;
	T = pW->P / (rho * R_ideal);
	
	/* Tpower = T^3.5 , Tpower1 = T^4.5 */
	Tpower = 1.0 / (T * T * T * sqrt(T));
	Tpower1 = Tpower / T;
	
	if(Sigma != NULL){
		Sigma[0] =  kappaes * rho;
		Sigma[1] =  kappaffR * rho * rho * Tpower;
		Sigma[2] =  kappaffP * rho * rho * Tpower;
		Sigma[3] =  kappaffP * rho * rho * Tpower;
	}
	
	if(dSigma != NULL){
		dSigma[0] = kappaes;
                dSigma[1] = 2.0 * rho * kappaffR * Tpower;
		dSigma[2] = 2.0 * rho * kappaffP * Tpower;
		dSigma[3] = 2.0 * rho * kappaffP * Tpower;

                dSigma[4] = 0.0;
                dSigma[5] = -3.5 *kappaffR * rho * rho * Tpower1;
		dSigma[6] = -3.5 *kappaffP * rho * rho * Tpower1;
                dSigma[7] = -3.5 *kappaffP * rho * rho * Tpower1;

		
	}
}
else{

	if(Sigma != NULL){
		Sigma[0] =  0.0;
		Sigma[1] =  0.0;
		Sigma[2] =  0.0;
		Sigma[3] =  0.0;
	}
	
	if(dSigma != NULL){
		dSigma[0] = 0.0;
                dSigma[1] = 0.0;
		dSigma[2] = 0.0;
		dSigma[3] = 0.0;

                dSigma[4] = 0.0;
                dSigma[5] = 0.0;
		dSigma[6] = 0.0;
                dSigma[7] = 0.0;
		
	}
}
	
	
	return; 
	
}

#endif



/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_beta(const GridS *pG, const int i, const int j, 
 *			      const int k)
 *  \brief Computes beta=P/(B^2/2)  
 */
static Real expr_beta(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  Real pre;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#if defined(MHD) || defined(RADIATION_MHD)
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;

  pre = hst_T(pG, i, j, k) * R_ideal * pG->U[k][j][i].d;

  return pre/(B2*0.5);

#else
  return 0.0;
#endif /* MHD */
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_ME(const GridS *pG, const int i, const int j, 
 *			    const int k)
 *  \brief  Computes B^2/8pi
 */
static Real expr_ME(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#if defined(RADIATION_MHD) || defined(MHD)
  B2=pG->U[k][j][i].B1c*pG->U[k][j][i].B1c;
  B2+=pG->U[k][j][i].B2c*pG->U[k][j][i].B2c;
  B2+=pG->U[k][j][i].B3c*pG->U[k][j][i].B3c;
  return (B2*0.5);
#else
  return 0.0;
#endif
}
/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_KE(const GridS *pG, const int i, const int j, 
 *			    const int k)
 *  \brief Computes dens*(Vx^2+Vy^2+Vz^2)/2
 */
static Real expr_KE(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,Vy,Vx,Vz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
  Vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  Vz = pG->U[k][j][i].M3/pG->U[k][j][i].d;

  return pG->U[k][j][i].d*(Vx*Vx + Vy*Vy + Vz*Vz)/2.0;

}

/*------------------------------------------------------------------------------
 * Hydro history variables:
 * hst_rho_Vx_dVy: Reynolds stress, added as history variable.
 * hst_rho_dVy2: KE in y-velocity fluctuations
 * hst_E_total: total energy (including tidal potential).
 */

static Real hst_rho_VrVphi(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3, Radius, vr, vphi, vx, vy;
  Real cosphi, sinphi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  Radius = x1;
  cosphi = x1/Radius;
  sinphi = x2/Radius;
  vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  vy = pG->U[k][j][i].M2/pG->U[k][j][i].d;
  vr = cosphi * vx + sinphi * vy;
  vphi = -sinphi * vx + cosphi * vy;

  return (pG->U[k][j][i].d * vr * vphi);

}



static Real hst_T(const GridS *pG, const int i, const int j, const int k)
{

	Real pressure;
	pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#if defined(MHD) || defined(RADIATION_MHD)
	pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	pressure *= Gamma - 1.0;

	return pressure/(R_ideal * pG->U[k][j][i].d);
}

static Real hst_P(const GridS *pG, const int i, const int j, const int k)
{

	Real pressure;
	pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#if defined(MHD) || defined(RADIATION_MHD)
	pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	pressure *= Gamma - 1.0;

	return pressure;
}


#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi =  PseudoNewton(x1, x2, x3);
 
  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi;
}
#endif /* ADIABATIC */



#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
static Real hst_sigmas(const GridS *pG, const int i, const int j, const int k)
{
	return pG->U[k][j][i].Sigma[0];
}

static Real hst_sigmaaP(const GridS *pG, const int i, const int j, const int k)
{
	return pG->U[k][j][i].Sigma[2];
}

#endif

static Real hst_rho2(const GridS *pG, const int i, const int j, const int k)
{
	
	return (pG->U[k][j][i].d * pG->U[k][j][i].d);
}

static Real hst_T2(const GridS *pG, const int i, const int j, const int k)
{
	Real temperature;
	temperature = hst_T(pG, i, j, k);

	return temperature * temperature;
}




static Real hst_rhofluxtop(const GridS *pG, const int i, const int j, const int k)
{
	int ID, lx3_id, rx3_id;

#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	lx3_id = pG->lx3_id;
	rx3_id = pG->rx3_id;
#else
	ID = 0;
	lx3_id = -1;
	rx3_id = -1;
#endif

	/* If it is on the boundary */
	if(((k == pG->ke) && ((rx3_id < 0) || (rx3_id <= ID))))
	{
		return pG->U[k][j][i].M3;

	}
	else
		return 0.0;

}


static Real hst_rhofluxbottom(const GridS *pG, const int i, const int j, const int k)
{
	int ID, lx3_id, rx3_id;

#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	lx3_id = pG->lx3_id;
	rx3_id = pG->rx3_id;
#else
	ID = 0;
	lx3_id = -1;
	rx3_id = -1;
#endif

	/* If it is on the boundary */
	if(((k == pG->ks) && ((lx3_id < 0) || (lx3_id >= ID))))
	{
		return pG->U[k][j][i].M3;

	}
	else
		return 0.0;

}




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

static Real hst_BrBphi(const GridS *pG, const int i, const int j, const int k)
{

  Real x1,x2,x3, Radius, Br, Bphi;
  Real cosphi, sinphi;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  Radius = x1;
  cosphi = x1/Radius;
  sinphi = x2/Radius;
  Br = cosphi * pG->U[k][j][i].B1c + sinphi * pG->U[k][j][i].B2c;
  Bphi = -sinphi *  pG->U[k][j][i].B1c + cosphi * pG->U[k][j][i].B2c;

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

/* Function for transfer module */
#ifdef RADIATION_TRANSFER

static Real Thermal_B_old(const GridS *pG, const int ifr, const int i, const int j, 
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

	//*pG->U[k][j][i].Sigma[2] / (pG->U[k][j][i].Sigma[0] + 1.0 * pG->U[k][j][i].Sigma[2]) * 1000.;
 


  return B;
}

static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k)
{
  
  return pow(pG->tgas[k][j][i],4);

}

static Real Transfereps_old(const GridS *pG,  const int ifr, const int i, const int j, 
		      const int k)
{

  Real eps;

  eps = pG->U[k][j][i].Sigma[2] / (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[2]);

  return eps;
  
}

static Real Transfereps(const GridS *pG,  const int ifr, const int i, const int j, 
		      const int k)
{

  Real eps = 0.0, siga, sigs;
  Real epsmin = 1.0E-10;

  if(pG->tgas[k][j][i] > 0.0) {
    siga = kappaffP * SQR(pG->U[k][j][i].d) / pow(pG->tgas[k][j][i],3.5);
    sigs = kappaes * pG->U[k][j][i].d;
    eps = siga / (siga + sigs);
  }
  if (eps < epsmin) eps=epsmin;

  //eps = pG->U[k][j][i].Sigma[2] / (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[2]);

  return eps;
  
}

static Real transfer_opacity(const GridS *pG,  const int ifr, const int i, const int j, 
			  const int k)
{
  Real siga = 0.0, sigs;
  if(pG->tgas[k][j][i] > 0.0) {
    siga = kappaffP * SQR(pG->U[k][j][i].d) / pow(pG->tgas[k][j][i],3.5);
  }
  sigs = kappaes * pG->U[k][j][i].d;	

  return siga + sigs;
  //  return pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[2];
  
}

static Real transfer_opacity_old(const GridS *pG,  const int ifr, const int i, const int j, 
			  const int k)
{
  
  return pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[2];
  
}
void Userwork_in_formal_solution(DomainS *pD)
{


}


#endif



