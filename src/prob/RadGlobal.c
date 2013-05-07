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

/*	Paczynski-Witt potential
*/ 
static Real PseudoNewton(const Real x1, const Real x2, const Real x3);

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

static Real hst_Erfluxbottom(const GridS *pG, const int i, const int j, const int k);
static Real hst_Erfluxtop(const GridS *pG, const int i, const int j, const int k);

static void output_1d(MeshS *pM, OutputS *pOut);

static void output_2d_binary(MeshS *pM, OutputS *pOut);

/* function to determine position */
void Radial_pos(int *rpos, const Real *out_r, const int nrmx, const Real Radius);


/* Function for inflow boundary condition */

void radMHD_rad_outflowks(GridS *pGrid);
void radMHD_outflowks(GridS *pGrid);
void radMHD_rad_outflowke(GridS *pGrid);
void radMHD_outflowke(GridS *pGrid);

void radMHD_rad_outflowie(GridS *pGrid);
void radMHD_outflowie(GridS *pGrid);
void radMHD_rad_periodis(GridS *pGrid);
void radMHD_periodis(GridS *pGrid);



void radMHD_rad_outflowje(GridS *pGrid);
void radMHD_outflowje(GridS *pGrid);
void radMHD_rad_periodjs(GridS *pGrid);
void radMHD_periodjs(GridS *pGrid);


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

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID
void radMHD_Mat_outflowke(MatrixS *pMat);
void radMHD_Mat_outflowks(MatrixS *pMat);
void radMHD_Mat_outflowie(MatrixS *pMat);
void radMHD_Mat_periodis(MatrixS *pMat);
void radMHD_Mat_outflowje(MatrixS *pMat);
void radMHD_Mat_periodjs(MatrixS *pMat);
#endif
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
static Real Rmask = 3.0;
static Real R0 = 50.0;

static int Periodix1_id = -1; /* store the process id for left x and y boundaries */
static int Periodjx1_id = -1;



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
	printf("Eratio: %f  Erflag: %d\n",Eratio,Erflag);
	
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
	T0= 1.e7;


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
	Real vs0 = 2.0;
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

	Radius = sqrt(x1*x1+x2*x2);
	distance = sqrt(x1*x1+x2*x2+x3*x3);

	langular = L0*pow(distance/R0,Lprofile)*(distance/(distance-1.0))/(R0/(R0-1.0));
	vphi = -langular/distance;

	vx = -vphi * x2/Radius;
	vy = vphi * x1/Radius;

	Effphi = -Crat*Crat/(2.0*(distance-1.0)) + (pow((langular/Radius),2.0))/(2.0*(1.0-Lprofile));
	Effphi0 = -Crat*Crat/(2.0*(R0-1.0)) + (pow((L0/R0),2.0))/(2.0*(1.0-Lprofile));
	tempphi = ((Effphi-Effphi0)/nindex)/(vs0*vs0);
	if((abs(tempphi) < 1.0) && (distance > Rmask)){
		density = rho0 * pow((1.0-tempphi),nindex);
		/* random perturbation */
		amp = 0.01;
	}
	else{
		density = dfloor;
		langular = 0.0;
		vx = 0.0;
		vy = 0.0;
		amp = 0.0;

		if(Radius > Rmask){
			Vc = -Crat*sqrt(distance/2.0)*distance/(distance*(distance-1.0));
			costheta = x2/Radius;
			sintheta = x1/Radius;
			vx = -Vc * costheta;
			vy = Vc * sintheta;
		}
	}


	if(density < dfloor){
		density = dfloor;
		langular = 0.0;
		vx = 0.0;
		vy = 0.0;
		amp = 0.0;
	}
	pressure = rho0*vs0*vs0*pow(density/rho0,Gamma)/Gamma;


	if(density > 2.0 * dfloor){
		Er = pow(pressure/(density * R_ideal), 4.0);
		if(Er > 1.0)
			Er = 1.0;
	}
	else{

		Er = pow(0.1, 4.0);
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
   dump_history_enroll(hst_Erfluxtop,"Erflux_top");
   dump_history_enroll(hst_Erfluxbottom,"Erflux_bottom");
   dump_history_enroll(hst_gravpot,"GravPot");	
   dump_history_enroll(hst_Lx,"Lx");
   dump_history_enroll(hst_Ly,"Ly");
   dump_history_enroll(hst_Lz,"Lz");

    frst = 0;
  }

	/* set boundary functions */

	/* outflow boundary condition along vertical direction */

	/* boundary condition for root level */
	if(pDomain->Level == 0){

		bvals_mhd_fun(pDomain, right_x3, radMHD_outflowke);

		bvals_rad_fun(pDomain, right_x3, radMHD_rad_outflowke);

		bvals_mhd_fun(pDomain, left_x3, radMHD_outflowks);

		bvals_rad_fun(pDomain, left_x3, radMHD_rad_outflowks);

		bvals_mhd_fun(pDomain, right_x1, radMHD_outflowie);

		bvals_rad_fun(pDomain, right_x1, radMHD_rad_outflowie);
		
		bvals_mhd_fun(pDomain, right_x2, radMHD_outflowje);

		bvals_rad_fun(pDomain, right_x2, radMHD_rad_outflowje);

		
	}

	bvals_mhd_fun(pDomain, left_x1, radMHD_periodis);

	bvals_rad_fun(pDomain, left_x1, radMHD_rad_periodis);

	bvals_mhd_fun(pDomain, left_x2, radMHD_periodjs);

	bvals_rad_fun(pDomain, left_x2, radMHD_rad_periodjs);


#ifdef MPI_PARALLEL

	/* set the MPI flag */
	if(iproc == 0 && jproc > 0)
		Periodix1_id = pDomain->GData[kproc][iproc][jproc].ID_Comm_Domain;
	else
		Periodix1_id = -1;

	if(jproc == 0 && iproc > 0)
		Periodjx1_id = pDomain->GData[kproc][iproc][jproc].ID_Comm_Domain;
	else
		Periodjx1_id = -1;


#endif




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

	fwrite(&Eratio,sizeof(Real),1,fp);
	fwrite(&Erflag,sizeof(int),1,fp);




	


#ifdef RESISTIVITY
	fwrite(&eta_Ohm,sizeof(Real),1,fp);
        fwrite(&Jrhomax,sizeof(Real),1,fp);
#endif

    fwrite(&zbtm,sizeof(Real),1,fp);
    fwrite(&ztop,sizeof(Real),1,fp);

    fwrite(&T0,sizeof(Real),1,fp);
    fwrite(&Periodix1_id,sizeof(int),1,fp);
    fwrite(&Periodjx1_id,sizeof(int),1,fp);

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
	fread(&Eratio,sizeof(Real),1,fp);
	fread(&Erflag,sizeof(int),1,fp);

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

       fread(&T0,sizeof(Real),1,fp);
       fread(&Periodix1_id,sizeof(int),1,fp);
       fread(&Periodjx1_id,sizeof(int),1,fp);

/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
  StaticGravPot = grav_vertical;
#else
  StaticGravPot = PseudoNewton;
#endif

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
   dump_history_enroll(hst_Erfluxtop,"Erflux_top");
   dump_history_enroll(hst_Erfluxbottom,"Erflux_bottom");
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

	if(nl == 0){
	
	bvals_mhd_fun(pD, right_x3, radMHD_outflowke);

	bvals_rad_fun(pD, right_x3, radMHD_rad_outflowke);

	bvals_mhd_fun(pD, left_x3, radMHD_outflowks);

	bvals_rad_fun(pD, left_x3, radMHD_rad_outflowks);

	bvals_mhd_fun(pD, right_x1, radMHD_outflowie);

	bvals_rad_fun(pD, right_x1, radMHD_rad_outflowie);

	
	bvals_mhd_fun(pD, right_x2, radMHD_outflowje);

	bvals_rad_fun(pD, right_x2, radMHD_rad_outflowje);
	
	}


	bvals_mhd_fun(pD, left_x1, radMHD_periodis);

	bvals_rad_fun(pD, left_x1, radMHD_rad_periodis);

	bvals_mhd_fun(pD, left_x2, radMHD_periodjs);

	bvals_rad_fun(pD, left_x2, radMHD_rad_periodjs);



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
  if(strcmp(name,"1d")==0) return output_1d;
  if(strcmp(name,"2db")==0) return output_2d_binary; 	
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
	distance = sqrt(x1*x1+x2*x2+x3*x3);
	Radius = sqrt(x1*x1 + x2*x2);
	
	if(Radius < 10.0){
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

	
	/* Apply resistivity inside the mask */
	distance = sqrt(x1*x1+x2*x2+x3*x3);

	if(distance < Rmask)
		*eta_O = eta0;

	*eta_H = 0.0;
	*eta_A = 0.0;

	return;
}
#endif



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


			cc_pos(pG,i,j,k,&x1,&x2,&x3);

			distance = sqrt(x1*x1+x2*x2+x3*x3);
			/* Reset the mask region */
			velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
			velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
			velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
							
			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);				
			if(velocity > 0.9 * Crat){
				pG->U[k][j][i].M1 = pG->U[k][j][i].d * velocity_x * 0.5 * Crat / velocity;
				pG->U[k][j][i].M2 = pG->U[k][j][i].d * velocity_y * 0.5 * Crat / velocity;
				pG->U[k][j][i].M3 = pG->U[k][j][i].d * velocity_z * 0.5 * Crat / velocity;
			}
				
							
							
			if(distance < Rmask){

				count = 0;
				temprho = 0.0;
				tempT = 0.0;
				Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
				T0 = Wtemp.P/(R_ideal * Wtemp.d);
				pG->U[k][j][i].M1 = 0.0;
				pG->U[k][j][i].M2 = 0.0;
				pG->U[k][j][i].M3 = 0.0;
				for(m=k; m<=k+1; m++){
				for(n=j; n<=j+1; n++){
				for(l=i; l<=i+1; l++){	
					cc_pos(pG,l,n,m,&x1,&x2,&x3);

					distance1 = sqrt(x1*x1+x2*x2+x3*x3);
					Wtemp = Cons_to_Prim(&(pG->U[m][n][l]));

			/*		if(((Rmask - distance) < (2.0 * pG->dx1)) && (distance1 > distance)){					
						tempT += Wtemp.P/(R_ideal * Wtemp.d);
						temprho += pG->U[m][n][l].d;
						pG->U[k][j][i].M1 += pG->U[m][n][l].M1;
						pG->U[k][j][i].M2 += pG->U[m][n][l].M2;
						pG->U[k][j][i].M3 += pG->U[m][n][l].M3;
						count++;
					}
					else if(distance1 > distance){
						temprho += pG->U[m][n][l].d;
						tempT += Wtemp.P/(R_ideal * Wtemp.d);
						pG->U[k][j][i].M1 += exp(-(distance1 - distance)/(2.0 * pG->dx1)) * pG->U[m][n][l].M1;
						pG->U[k][j][i].M2 += exp(-(distance1 - distance)/(2.0* pG->dx1)) * pG->U[m][n][l].M2;
						pG->U[k][j][i].M3 += exp(-(distance1 - distance)/(2.0 * pG->dx1)) * pG->U[m][n][l].M3;
						count++;
					}
			*/
					if(distance1 > distance){
                                                temprho += pG->U[m][n][l].d;
                                                tempT += Wtemp.P/(R_ideal * Wtemp.d);
                                                pG->U[k][j][i].M1 += exp(-(distance1 - distance)/(4.0 * pG->dx1)) * pG->U[m][n][l].M1;
                                                pG->U[k][j][i].M2 += exp(-(distance1 - distance)/(4.0* pG->dx1)) * pG->U[m][n][l].M2;
                                                pG->U[k][j][i].M3 += exp(-(distance1 - distance)/(4.0 * pG->dx1)) * pG->U[m][n][l].M3;
                                                count++;
                                        }


				}
				}
				}

				if(count > 0){
					pG->U[k][j][i].d = temprho/count;
					tempT /= count;
					pG->U[k][j][i].M1 /= count;
					pG->U[k][j][i].M2 /= count;
					pG->U[k][j][i].M3 /= count;

					/* reconstruct E */
					pG->U[k][j][i].E = tempT * R_ideal * pG->U[k][j][i].d/(Gamma - 1.0) + 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3)/pG->U[k][j][i].d;
#if defined(MHD) || defined(RADIATION_MHD)
					pG->U[k][j][i].E += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
					
				}
				else{
					/* keep original T */
					/* reconstruct E */
					pG->U[k][j][i].E = T0 * R_ideal * pG->U[k][j][i].d/(Gamma - 1.0) + 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3)/pG->U[k][j][i].d;
#if defined(MHD) || defined(RADIATION_MHD)
					pG->U[k][j][i].E += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif


				}
							
				
				

			}/* Distance < Rmask */
			




   			badcellflag = 0;

			velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
                         velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                         velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;

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
	Real distance, GM, potential;
	/* The dimensionless number */
	GM = Crat*Crat*0.5;
	distance = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
	if(distance > Rmask)
		potential = -GM/(distance - 1.0);
	else
		potential = -GM/(2.0 - 1.0);
	
	return potential;
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
  Radius = sqrt(x1*x1+x2*x2);
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



static Real hst_Erfluxtop(const GridS *pG, const int i, const int j, const int k)
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
		return pG->U[k][j][i].Fr3;

	}
	else
		return 0.0;

}


static Real hst_Erfluxbottom(const GridS *pG, const int i, const int j, const int k)
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
		return pG->U[k][j][i].Fr3;

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
  Radius = sqrt(x1*x1+x2*x2);
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




/* Function for  boundary condition */

/* Vertical boundary */
void radMHD_outflowke(GridS *pGrid)
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
	Real Sigma[NOPACITY];
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

                                /* Just copy the magnetic field */
                /*                if(fabs(pGrid->B1i[ke+k-2][j][i]) > fabs(pGrid->B1i[ke+k-1][j][i]) ){
                                         pGrid->B1i[ke+k][j][i] = 2.0 * pGrid->B1i[ke+k-1][j][i] - pGrid->B1i[ke+k-2][j][i];
                                        if(fabs(pGrid->B1i[ke+k][j][i]) > fabs(pGrid->B1i[ke+k-1][j][i])){
                                                pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
                                        }
                                }
                                else{
                                        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
                                }

                                if(fabs(pGrid->B2i[ke+k-2][j][i]) > fabs(pGrid->B2i[ke+k-1][j][i]) ){
                                         pGrid->B2i[ke+k][j][i] = 2.0 * pGrid->B2i[ke+k-1][j][i] - pGrid->B2i[ke+k-2][j][i];
                                        if(fabs(pGrid->B2i[ke+k][j][i]) > fabs(pGrid->B2i[ke+k-1][j][i])){
                                                pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];
                                        }

                                }
                                else{
                                         pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];
                                }
		*/
			

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
	/* Then update B3 */
/*	for(k=1;k<=nghost;k++){
                for(j=js-nghost;j<=je+nghost;j++){
                        for(i=is-nghost;i<=ie+nghost;i++){
				if(k<nghost && i < ie+nghost && j < je+nghost){
                                        pGrid->B3i[ke+k+1][j][i] = pGrid->B3i[ke+k][j][i] - pGrid->dx3 * ((pGrid->B1i[ke+k][j][i+1] - pGrid->B1i[ke+k][j][i]) / pGrid->dx1
                                                                        + (pGrid->B2i[ke+k][j+1][i] - pGrid->B2i[ke+k][j][i]) / pGrid->dx2);
                                }
                                else if(k < nghost){
                                        pGrid->B3i[ke+k+1][j][i] = pGrid->B3i[ke+k][j][i];
                                }


			}
		}
	}
*/

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
#ifdef RADIATION_MHD

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



#ifdef RADIATION_MHD
		pGrid->U[ke+k][j][i].E += 0.5 * (pGrid->U[ke+k][j][i].B1c * pGrid->U[ke+k][j][i].B1c + pGrid->U[ke+k][j][i].B2c * pGrid->U[ke+k][j][i].B2c + pGrid->U[ke+k][j][i].B3c * pGrid->U[ke+k][j][i].B3c);

#endif

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
	
      		
      		}
    		}
	}
  
}


void radMHD_rad_outflowke(GridS *pGrid)
{
  	int i,j,js,je;
	int ks, is,ie,k, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real Fr0x, Fr0y, Fr0z;
	Real velocity_x, velocity_y, velocity_z;
	Real velocity_x1, velocity_y1, velocity_z1;
	Real Sigma_t, Sigma_t1;
	Real Eratio, reducefactor, tau;


	Real x1, x2, x3;
	

	 for(k=1;  k<=nghost;  k++) {
	for(j=js-nghost; j<=je+nghost;j++){
	for(i=is-nghost; i<=ie+nghost; i++){
	   
		cc_pos(pGrid, i, j,ke+k, &x1, &x2, &x3);

		pGrid->U[ke+k][j][i].Edd_11 = pGrid->U[ke][j][i].Edd_11;
		pGrid->U[ke+k][j][i].Edd_22 = pGrid->U[ke][j][i].Edd_22;
		pGrid->U[ke+k][j][i].Edd_21 = pGrid->U[ke][j][i].Edd_21;
		pGrid->U[ke+k][j][i].Edd_31 = pGrid->U[ke][j][i].Edd_31;
		pGrid->U[ke+k][j][i].Edd_32 = pGrid->U[ke][j][i].Edd_32;
		pGrid->U[ke+k][j][i].Edd_33 = pGrid->U[ke][j][i].Edd_33;
	
		matrix_alpha(0.0, pGrid->U[ke+k][j][i].Sigma, pGrid->dt, pGrid->U[ke+k][j][i].Edd_33, 0.0, &reducefactor, 0, dz);
	
	/*
		tau = pGrid->dt * Crat * (pGrid->U[ke+k][j][i].Sigma[0] + pGrid->U[ke+k][j][i].Sigma[1]);
		tau = tau * tau / (2.0 * pGrid->U[ke+k][j][i].Edd_33);

		if(tau > 0.001)
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
		else
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - 0.5 * tau));			
	*/
		Sigma_t = 0.5 * (pGrid->U[ke+k][j][i].Sigma[0] + pGrid->U[ke+k][j][i].Sigma[1] + pGrid->U[ke+k-1][j][i].Sigma[0] + pGrid->U[ke+k-1][j][i].Sigma[1]);
		Sigma_t1 = 0.5 * (pGrid->U[ke+k-1][j][i].Sigma[0] + pGrid->U[ke+k-1][j][i].Sigma[1] + pGrid->U[ke+k-2][j][i].Sigma[0] + pGrid->U[ke+k-2][j][i].Sigma[1]);
	
		velocity_x = pGrid->U[ke+k-1][j][i].M1 / pGrid->U[ke+k-1][j][i].d;
		velocity_y = pGrid->U[ke+k-1][j][i].M2 / pGrid->U[ke+k-1][j][i].d;
#ifdef FARGO
		velocity_y -= qshear * Omega_0 * x1;

#endif
		velocity_z = pGrid->U[ke+k-1][j][i].M3 / pGrid->U[ke+k-1][j][i].d;


		
		Fr0x = pGrid->U[ke+k-1][j][i].Fr1 - ((1.0 + pGrid->U[ke+k-1][j][i].Edd_11) * velocity_x + pGrid->U[ke+k-1][j][i].Edd_21 * velocity_y + pGrid->U[ke+k-1][j][i].Edd_31 * velocity_z)* pGrid->U[ke+k-1][j][i].Er / Crat;

		Fr0y = pGrid->U[ke+k-1][j][i].Fr2 - ((1.0 + pGrid->U[ke+k-1][j][i].Edd_22) * velocity_y + pGrid->U[ke+k-1][j][i].Edd_21 * velocity_x + pGrid->U[ke+k-1][j][i].Edd_32 * velocity_z)* pGrid->U[ke+k-1][j][i].Er / Crat;

		Fr0z = pGrid->U[ke+k-1][j][i].Fr3 - ((1.0 + pGrid->U[ke+k-1][j][i].Edd_33) * velocity_z + pGrid->U[ke+k-1][j][i].Edd_31 * velocity_x + pGrid->U[ke+k-1][j][i].Edd_32 * velocity_y)* pGrid->U[ke+k-1][j][i].Er / Crat;

		velocity_x1 = pGrid->U[ke+k][j][i].M1 / pGrid->U[ke+k][j][i].d;
		velocity_y1 = pGrid->U[ke+k][j][i].M2 / pGrid->U[ke+k][j][i].d;
		velocity_z1 = pGrid->U[ke+k][j][i].M3 / pGrid->U[ke+k][j][i].d;
#ifdef FARGO
		velocity_y1 -= qshear * Omega_0 * x1;

#endif

/*
		pGrid->U[ke+k][j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ke+k][j][i].Edd_11) * velocity_x1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_y1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ke+k][j][i].Edd_22) * velocity_y1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;



		pGrid->U[ke+k][j][i].Fr1 = pGrid->U[ke][j][i].Fr1;
		pGrid->U[ke+k][j][i].Fr2 = pGrid->U[ke][j][i].Fr2;
		if(pGrid->U[ke][j][i].Fr3 > 0.0)
			pGrid->U[ke+k][j][i].Fr3 = pGrid->U[ke][j][i].Fr3;
		else
			pGrid->U[ke+k][j][i].Fr3 = 0.0;
*/
/*
		if(Fr0z > 0.0){			
			pGrid->U[ke+k][j][i].Er = pGrid->U[ke+k-1][j][i].Er - pGrid->dx3 * Sigma_t * Fr0z / pGrid->U[ke+k][j][i].Edd_33;			

			if(pGrid->U[ke+k][j][i].Er < 0.0){
				pGrid->U[ke+k][j][i].Er = pGrid->U[ke+k-1][j][i].Er;

			}

			
		}	
		else{
			pGrid->U[ke+k][j][i].Er = pGrid->U[ke+k-1][j][i].Er;
			
		}

		Eratio = 1.0 + Sigma_t * (1.0 - pGrid->U[ke+k-2][j][i].Edd_33 * pGrid->U[ke+k-2][j][i].Er / (pGrid->U[ke+k-1][j][i].Edd_33 * pGrid->U[ke+k-1][j][i].Er)) / Sigma_t1;
		if(Eratio < 0.0 || Eratio > 1.0){
				Eratio = 1.0;
		}

		pGrid->U[ke+k][j][i].Er = Eratio * pGrid->U[ke+k-1][j][i].Er;		
*/
	
		Eratio = pGrid->U[ke+k][j][i].Edd_33 + 0.5 * pGrid->dx3 * reducefactor * Sigma_t;

		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
		pGrid->U[ke+k][j][i].Er = (pGrid->U[ke+k-1][j][i].Edd_33 * pGrid->U[ke+k-1][j][i].Er - 0.5 * pGrid->dx3 * Sigma_t * Fr0z ) / Eratio;

		if((pGrid->U[ke+k][j][i].Er > pGrid->U[ke+k-1][j][i].Er) || (pGrid->U[ke+k][j][i].Er < 0.0)){
			Eratio = pGrid->U[ke+k][j][i].Edd_33 + pGrid->dx3 * reducefactor * Sigma_t;
			pGrid->U[ke+k][j][i].Er = pGrid->U[ke+k-1][j][i].Edd_33 * pGrid->U[ke+k-1][j][i].Er / Eratio;
		}	
			
		
		Fr0z = reducefactor * pGrid->U[ke+k][j][i].Er;

		pGrid->U[ke+k][j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ke+k][j][i].Edd_11) * velocity_x1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_y1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ke+k][j][i].Edd_22) * velocity_y1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
		pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;


/*
		if(Fr0z > 0.0)
	                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;
		else
			pGrid->U[ke+k][j][i].Fr3 = ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;

*/			


		


    		}
		}
	}
  
}



void radMHD_outflowks(GridS *pGrid)
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
	Real Sigma[NOPACITY];
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
#ifdef RADIATION_MHD

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



#ifdef RADIATION_MHD
		pGrid->U[ks-k][j][i].E += 0.5 * (pGrid->U[ks-k][j][i].B1c * pGrid->U[ks-k][j][i].B1c + pGrid->U[ks-k][j][i].B2c * pGrid->U[ks-k][j][i].B2c + pGrid->U[ks-k][j][i].B3c * pGrid->U[ks-k][j][i].B3c);

#endif

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
	
      		
      }
    }
  
  
}


void radMHD_rad_outflowks(GridS *pGrid)
{
  	int i,j,js,je;
	int ks, is,ie,k, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real Fr0x, Fr0y, Fr0z;
	Real velocity_x, velocity_y, velocity_z;
	Real velocity_x1, velocity_y1, velocity_z1;
	Real Sigma_t, Sigma_t1;
	Real Eratio, reducefactor, tau;


	Real x1, x2, x3;


	 for(k=1;  k<=nghost;  k++) {
	for(j=js-nghost; j<=je+nghost;j++){
	for(i=is-nghost; i<=ie+nghost; i++){
	   
		cc_pos(pGrid, i, j,ks-k, &x1, &x2, &x3);

		pGrid->U[ks-k][j][i].Edd_11 = pGrid->U[ks][j][i].Edd_11;
		pGrid->U[ks-k][j][i].Edd_22 = pGrid->U[ks][j][i].Edd_22;
		pGrid->U[ks-k][j][i].Edd_21 = pGrid->U[ks][j][i].Edd_21;
		pGrid->U[ks-k][j][i].Edd_31 = pGrid->U[ks][j][i].Edd_31;
		pGrid->U[ks-k][j][i].Edd_32 = pGrid->U[ks][j][i].Edd_32;
		pGrid->U[ks-k][j][i].Edd_33 = pGrid->U[ks][j][i].Edd_33;

		matrix_alpha(0.0, pGrid->U[ks-k][j][i].Sigma, pGrid->dt, pGrid->U[ks-k][j][i].Edd_33, 0.0, &reducefactor, 0, dz);
	
	
	/*	tau = pGrid->dt * Crat * (pGrid->U[ks-k][j][i].Sigma[0] + pGrid->U[ks-k][j][i].Sigma[1]);
		tau = tau * tau / (2.0 * pGrid->U[ks-k][j][i].Edd_33);

		if(tau > 0.001)
			reducefactor = sqrt(pGrid->U[ks-k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
		else
			reducefactor = sqrt(pGrid->U[ks-k][j][i].Edd_33 * (1.0 - 0.5 * tau));	
	*/
		Sigma_t = 0.5 * (pGrid->U[ks-k][j][i].Sigma[0] + pGrid->U[ks-k][j][i].Sigma[1] + pGrid->U[ks-k+1][j][i].Sigma[0] + pGrid->U[ks-k+1][j][i].Sigma[1]);
		Sigma_t1 = 0.5 * (pGrid->U[ks-k+1][j][i].Sigma[0] + pGrid->U[ks-k+1][j][i].Sigma[1] + pGrid->U[ks-k+2][j][i].Sigma[0] + pGrid->U[ks-k+2][j][i].Sigma[1]);
	
		velocity_x = pGrid->U[ks-k+1][j][i].M1 / pGrid->U[ks-k+1][j][i].d;
		velocity_y = pGrid->U[ks-k+1][j][i].M2 / pGrid->U[ks-k+1][j][i].d;
#ifdef FARGO
         	velocity_y -= qshear * Omega_0 * x1;

#endif
		velocity_z = pGrid->U[ks-k+1][j][i].M3 / pGrid->U[ks-k+1][j][i].d;

		

		Fr0x = pGrid->U[ks-k+1][j][i].Fr1 - ((1.0 + pGrid->U[ks-k+1][j][i].Edd_11) * velocity_x + pGrid->U[ks-k+1][j][i].Edd_21 * velocity_y + pGrid->U[ks-k+1][j][i].Edd_31 * velocity_z)* pGrid->U[ks-k+1][j][i].Er / Crat;

		Fr0y = pGrid->U[ks-k+1][j][i].Fr2 - ((1.0 + pGrid->U[ks-k+1][j][i].Edd_22) * velocity_y + pGrid->U[ks-k+1][j][i].Edd_21 * velocity_x + pGrid->U[ks-k+1][j][i].Edd_32 * velocity_z)* pGrid->U[ks-k+1][j][i].Er / Crat;

		Fr0z = pGrid->U[ks-k+1][j][i].Fr3 - ((1.0 + pGrid->U[ks-k+1][j][i].Edd_33) * velocity_z + pGrid->U[ks-k+1][j][i].Edd_31 * velocity_x + pGrid->U[ks-k+1][j][i].Edd_32 * velocity_y)* pGrid->U[ks-k+1][j][i].Er / Crat;

		velocity_x1 = pGrid->U[ks-k][j][i].M1 / pGrid->U[ks-k][j][i].d;
		velocity_y1 = pGrid->U[ks-k][j][i].M2 / pGrid->U[ks-k][j][i].d;
#ifdef FARGO
		velocity_y1 -= qshear * Omega_0 * x1;

#endif

		velocity_z1 = pGrid->U[ks-k][j][i].M3 / pGrid->U[ks-k][j][i].d;


/*
		pGrid->U[ks-k][j][i].Fr1 = pGrid->U[ks-k+1][j][i].Fr1;
		pGrid->U[ks-k][j][i].Fr2 = pGrid->U[ks-k+1][j][i].Fr2;
		if(pGrid->U[ks][j][i].Fr3 < 0.0 )
			pGrid->U[ks-k][j][i].Fr3 = pGrid->U[ks][j][i].Fr3;
		else
			pGrid->U[ks-k][j][i].Fr3 = 0.0;

		if(Fr0z < 0.0){			
			pGrid->U[ks-k][j][i].Er = pGrid->U[ks-k+1][j][i].Er + pGrid->dx3 * Sigma_t * Fr0z / pGrid->U[ks-k][j][i].Edd_33;			

			if(pGrid->U[ks-k][j][i].Er < 0.0){
				pGrid->U[ks-k][j][i].Er = pGrid->U[ks-k+1][j][i].Er;

			}

			
		}
		else{
			pGrid->U[ks-k][j][i].Er = pGrid->U[ks-k+1][j][i].Er;
			
		}
*/
/*
		 Eratio = 1.0 + Sigma_t * (1.0 - pGrid->U[ks-k+2][j][i].Edd_33 * pGrid->U[ks-k+2][j][i].Er / (pGrid->U[ks-k+1][j][i].Edd_33 * pGrid->U[ks-k+1][j][i].Er)) / Sigma_t1;
                if(Eratio < 0.0 || Eratio > 1.0){
                                Eratio = 1.0;
                }

                pGrid->U[ks-k][j][i].Er = Eratio * pGrid->U[ks-k+1][j][i].Er;

*/

		Eratio = pGrid->U[ks-k][j][i].Edd_33 + 0.5 * pGrid->dx3 * reducefactor * Sigma_t;

		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
		pGrid->U[ks-k][j][i].Er = (pGrid->U[ks-k+1][j][i].Edd_33 * pGrid->U[ks-k+1][j][i].Er + 0.5 * pGrid->dx3 * Sigma_t * Fr0z ) / Eratio;

		if((pGrid->U[ks-k][j][i].Er > pGrid->U[ks-k+1][j][i].Er) || (pGrid->U[ks-k][j][i].Er < 0.0)){
			Eratio = pGrid->U[ks-k][j][i].Edd_33 + pGrid->dx3 * reducefactor * Sigma_t;
			pGrid->U[ks-k][j][i].Er = pGrid->U[ks-k+1][j][i].Edd_33 * pGrid->U[ks-k+1][j][i].Er / Eratio;
		}
		
		Fr0z = -reducefactor * pGrid->U[ks-k][j][i].Er;
		
		

		pGrid->U[ks-k][j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ks-k][j][i].Edd_11) * velocity_x1 + pGrid->U[ks-k][j][i].Edd_21 * velocity_y1 + pGrid->U[ks-k][j][i].Edd_31 * velocity_z1)* pGrid->U[ks-k][j][i].Er / Crat;
		pGrid->U[ks-k][j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ks-k][j][i].Edd_22) * velocity_y1 + pGrid->U[ks-k][j][i].Edd_21 * velocity_x1 + pGrid->U[ks-k][j][i].Edd_32 * velocity_z1)* pGrid->U[ks-k][j][i].Er / Crat;
		pGrid->U[ks-k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ks-k][j][i].Edd_33) * velocity_z1 + pGrid->U[ks-k][j][i].Edd_31 * velocity_x1 + pGrid->U[ks-k][j][i].Edd_32 * velocity_y1)* pGrid->U[ks-k][j][i].Er / Crat;



/*		if(Fr0z < 0.0)
	                pGrid->U[ks-k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ks-k][j][i].Edd_33) * velocity_z1 + pGrid->U[ks-k][j][i].Edd_31 * velocity_x1 + pGrid->U[ks-k][j][i].Edd_32 * velocity_y1)* pGrid->U[ks-k][j][i].Er / Crat;
		else
			pGrid->U[ks-k][j][i].Fr3 = ((1.0 + pGrid->U[ks-k][j][i].Edd_33) * velocity_z1 + pGrid->U[ks-k][j][i].Edd_31 * velocity_x1 + pGrid->U[ks-k][j][i].Edd_32 * velocity_y1)* pGrid->U[ks-k][j][i].Er / Crat;
*/
    		}
		}
	}
  
}


/* x boundary */
/* Assume to be keplerian motion  *
 * vx = v_c sin theta *
 * vy = -v_c cos theta */

void radMHD_outflowie(GridS *pGrid)
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
	Real Vc, distance, Radius; /* The Keplerian circular velocity */
	Real Sigma[NOPACITY];
	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;

	static Real x3t;
	x3t = ztop - 0.5 * pGrid->dx3;
	PrimS Wopacity;


#if defined(MHD) || defined(RADIATION_MHD)
	for(k=ks-nghost; k<=ke+nghost; k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=1; i<=nghost; i++){
                                /* Just copy the magnetic field */
                /*                if(fabs(pGrid->B1i[ke+k-2][j][i]) > fabs(pGrid->B1i[ke+k-1][j][i]) ){
                                         pGrid->B1i[ke+k][j][i] = 2.0 * pGrid->B1i[ke+k-1][j][i] - pGrid->B1i[ke+k-2][j][i];
                                        if(fabs(pGrid->B1i[ke+k][j][i]) > fabs(pGrid->B1i[ke+k-1][j][i])){
                                                pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
                                        }
                                }
                                else{
                                        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
                                }

                                if(fabs(pGrid->B2i[ke+k-2][j][i]) > fabs(pGrid->B2i[ke+k-1][j][i]) ){
                                         pGrid->B2i[ke+k][j][i] = 2.0 * pGrid->B2i[ke+k-1][j][i] - pGrid->B2i[ke+k-2][j][i];
                                        if(fabs(pGrid->B2i[ke+k][j][i]) > fabs(pGrid->B2i[ke+k-1][j][i])){
                                                pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];
                                        }

                                }
                                else{
                                         pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];
                                }
		*/
			

			if(pGrid->U[k][j][ie].M1 < 0.0){
				pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie+i-1];
                                pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie+i-1];
                        }
                        else{
                                pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie+i-1];
                                pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie+i-1];
                        }


				if(i<nghost)
					pGrid->B1i[k][j][ie+i+1] = pGrid->B1i[k][j][ie+i];


			}
		}
	}


	/* Now update the cell centered values */
	for(k=ks-nghost;k<=ke+nghost;k++){
		for(j=js-nghost;j<=je+nghost;j++){
			for(i=1;i<=nghost;i++){
				 if(pGrid->U[k][j][ie].M1 < 0.0){
                                         pGrid->U[k][j][ie+i].B3c = pGrid->U[k][j][ie+i-1].B3c;
                                         pGrid->U[k][j][ie+i].B2c = pGrid->U[k][j][ie+i-1].B2c;
                                 }
                                 else{
                                         pGrid->U[k][j][ie+i].B3c = pGrid->U[k][j][ie+i].B3c;
                                         pGrid->U[k][j][ie+i].B2c = pGrid->U[k][j][ie+i].B2c;
                                 }
				
				if(i<nghost){
					pGrid->U[k][j][ie+i].B1c = 0.5 * (pGrid->B1i[k][j][ie+i] + pGrid->B1i[k][j][ie+i+1]);
				}
				else{
					pGrid->U[k][j][ie+i].B1c = pGrid->B1i[k][j][ie+i];
				}

			}
		}
	}
		



#endif /* MHD */
	
	for (k=ks-nghost;  k<=ke+nghost;  k++) {
 		for(j=js-nghost;j<=je+nghost;j++){
    			for(i=1; i<=nghost; i++){
	    
	

		cc_pos(pGrid, ie+i, j,k, &x1, &x2, &x3);
#if defined(MHD) || defined(RADIATION_MHD)
		Bx = pGrid->U[k][j][ie+i].B1c;
		By = pGrid->U[k][j][ie+i].B2c;
		Bz = pGrid->U[k][j][ie+i].B3c;
		
		B = sqrt(Bx * Bx + By * By + Bz * Bz);
#endif

		if(pGrid->U[k][j][ie].d < TINY_NUMBER)
			pGrid->U[k][j][ie].d = pGrid->U[4][4][ie].d;

		distance = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
		Radius = sqrt(x1 * x1 + x2 * x2);
		Vc = Crat*sqrt(1.0/(2.0*distance))/(distance-1.0);

		density  = pGrid->U[k][j][ie].d;

		velocity1 = pGrid->U[k][j][ie].M1 / pGrid->U[k][j][ie].d;
		velocity2 = pGrid->U[k][j][ie].M2 / pGrid->U[k][j][ie].d;
		velocity3 = pGrid->U[k][j][ie].M3 / pGrid->U[k][j][ie].d;

		pressure = (pGrid->U[k][j][ie].E - 0.5 * density * (velocity1 * velocity1 + velocity2 * velocity2 + velocity3 * velocity3)) * (Gamma - 1.0);
#ifdef RADIATION_MHD

		pressure -= 0.5 * (pGrid->U[k][j][ie].B1c * pGrid->U[k][j][ie].B1c + pGrid->U[k][j][ie].B2c * pGrid->U[k][j][ie].B2c + pGrid->U[k][j][ie].B3c * pGrid->U[k][j][ie].B3c) * (Gamma - 1.0);
#endif

		T = pressure / (pGrid->U[k][j][ie].d * R_ideal);
/*
		if(T > Tfloor)
			T = Tfloor;
*/
		/* Tfloor is usually set as initial temperature in the ghost zoner */

		/* now reset the Keplerian velocity */
				

		velocity1 = Vc * x2;
		velocity2 = -Vc * x1;
		
		/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
		
	
		pGrid->U[k][j][ie+i].d = pGrid->U[k][j][ie].d;



		
      		pGrid->U[k][j][ie+i].M1 = velocity1 * pGrid->U[k][j][ie+i].d;
		pGrid->U[k][j][ie+i].M2 = velocity2 * pGrid->U[k][j][ie+i].d;
			
		
		pGrid->U[k][j][ie+i].M3 = velocity3 * pGrid->U[k][j][ie+i].d;

		pressure = T * R_ideal * pGrid->U[k][j][ie+i].d;

		pGrid->U[k][j][ie+i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[k][j][ie+i].M1 * pGrid->U[k][j][ie+i].M1 + pGrid->U[k][j][ie+i].M2 * pGrid->U[k][j][ie+i].M2 + pGrid->U[k][j][ie+i].M3 * pGrid->U[k][j][ie+i].M3) / pGrid->U[k][j][ie+i].d;



#ifdef RADIATION_MHD
		pGrid->U[k][j][ie+i].E += 0.5 * (pGrid->U[k][j][ie+i].B1c * pGrid->U[k][j][ie+i].B1c + pGrid->U[k][j][ie+i].B2c * pGrid->U[k][j][ie+i].B2c + pGrid->U[k][j][ie+i].B3c * pGrid->U[k][j][ie+i].B3c);

#endif

		Wopacity = Cons_to_Prim(&pGrid->U[k][j][ie+i]);
		/* Add background shearing */
#ifdef FARGO	
		cc_pos(pGrid,ie+i,j,k,&x1,&x2,&x3);
		Wopacity.V2 -= qshear * Omega_0 * x1;		
#endif
			
		Thindiskopacity(&Wopacity, Sigma, NULL);
		
		if(T < 2.0 * TINY_NUMBER)
			Sigma[0] = kappaes * pGrid->U[k][j][ie+i].d;

		for(m=0; m<NOPACITY;m++){
			pGrid->U[k][j][ie+i].Sigma[m] = Sigma[m];
		}	
	
      		
      		}
    		}
	}
  
}


void radMHD_rad_outflowie(GridS *pGrid)
{
  	int i,j,js,je;
	int ks, is,ie,k, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real Fr0x, Fr0y, Fr0z;
	Real velocity_x, velocity_y, velocity_z;
	Real velocity_x1, velocity_y1, velocity_z1;
	Real Sigma_t, Sigma_t1;
	Real Eratio, reducefactor, tau;


	Real x1, x2, x3;
	

	for(k=ks-nghost;  k<=ke+nghost;  k++) {
	for(j=js-nghost; j<=je+nghost;j++){
	for(i=1; i<=nghost; i++){
	   
		cc_pos(pGrid, ie+i, j,k, &x1, &x2, &x3);

		pGrid->U[k][j][ie+i].Edd_11 = pGrid->U[k][j][ie].Edd_11;
		pGrid->U[k][j][ie+i].Edd_22 = pGrid->U[k][j][ie].Edd_22;
		pGrid->U[k][j][ie+i].Edd_21 = pGrid->U[k][j][ie].Edd_21;
		pGrid->U[k][j][ie+i].Edd_31 = pGrid->U[k][j][ie].Edd_31;
		pGrid->U[k][j][ie+i].Edd_32 = pGrid->U[k][j][ie].Edd_32;
		pGrid->U[k][j][ie+i].Edd_33 = pGrid->U[k][j][ie].Edd_33;
	
		matrix_alpha(0.0, pGrid->U[k][j][ie+i].Sigma, pGrid->dt, pGrid->U[k][j][ie+i].Edd_11, 0.0, &reducefactor, 0, dx);
	
	/*
		tau = pGrid->dt * Crat * (pGrid->U[ke+k][j][i].Sigma[0] + pGrid->U[ke+k][j][i].Sigma[1]);
		tau = tau * tau / (2.0 * pGrid->U[ke+k][j][i].Edd_33);

		if(tau > 0.001)
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
		else
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - 0.5 * tau));			
	*/
		Sigma_t = 0.5 * (pGrid->U[k][j][ie+i].Sigma[0] + pGrid->U[k][j][ie+i].Sigma[1] + pGrid->U[k][j][ie+i-1].Sigma[0] + pGrid->U[k][j][ie+i-1].Sigma[1]);
		Sigma_t1 = 0.5 * (pGrid->U[k][j][ie+i-1].Sigma[0] + pGrid->U[k][j][ie+i-1].Sigma[1] + pGrid->U[k][j][ie+i-2].Sigma[0] + pGrid->U[k][j][ie+i-2].Sigma[1]);
	
		velocity_x = pGrid->U[k][j][ie+i-1].M1 / pGrid->U[k][j][ie+i-1].d;
		velocity_y = pGrid->U[k][j][ie+i-1].M2 / pGrid->U[k][j][ie+i-1].d;
#ifdef FARGO
		velocity_y -= qshear * Omega_0 * x1;

#endif
		velocity_z = pGrid->U[k][j][ie+i-1].M3 / pGrid->U[k][j][ie+i-1].d;


		
		Fr0x = pGrid->U[k][j][ie+i-1].Fr1 - ((1.0 + pGrid->U[k][j][ie+i-1].Edd_11) * velocity_x + pGrid->U[k][j][ie+i-1].Edd_21 * velocity_y + pGrid->U[k][j][ie+i-1].Edd_31 * velocity_z)* pGrid->U[k][j][ie+i-1].Er / Crat;

		Fr0y = pGrid->U[k][j][ie+i-1].Fr2 - ((1.0 + pGrid->U[k][j][ie+i-1].Edd_22) * velocity_y + pGrid->U[k][j][ie+i-1].Edd_21 * velocity_x + pGrid->U[k][j][ie+i-1].Edd_32 * velocity_z)* pGrid->U[k][j][ie+i-1].Er / Crat;

		Fr0z = pGrid->U[k][j][ie+i-1].Fr3 - ((1.0 + pGrid->U[k][j][ie+i-1].Edd_33) * velocity_z + pGrid->U[k][j][ie+i-1].Edd_31 * velocity_x + pGrid->U[k][j][ie+i-1].Edd_32 * velocity_y)* pGrid->U[k][j][ie+i-1].Er / Crat;

		velocity_x1 = pGrid->U[k][j][ie+i].M1 / pGrid->U[k][j][ie+i].d;
		velocity_y1 = pGrid->U[k][j][ie+i].M2 / pGrid->U[k][j][ie+i].d;
		velocity_z1 = pGrid->U[k][j][ie+i].M3 / pGrid->U[k][j][ie+i].d;

/*
		pGrid->U[ke+k][j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ke+k][j][i].Edd_11) * velocity_x1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_y1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ke+k][j][i].Edd_22) * velocity_y1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;



		pGrid->U[ke+k][j][i].Fr1 = pGrid->U[ke][j][i].Fr1;
		pGrid->U[ke+k][j][i].Fr2 = pGrid->U[ke][j][i].Fr2;
		if(pGrid->U[ke][j][i].Fr3 > 0.0)
			pGrid->U[ke+k][j][i].Fr3 = pGrid->U[ke][j][i].Fr3;
		else
			pGrid->U[ke+k][j][i].Fr3 = 0.0;
*/

		if(Fr0x > 0.0){			
			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie+i-1].Er - pGrid->dx1 * Sigma_t * Fr0x / pGrid->U[k][j][ie+i-1].Edd_11;			

			if(pGrid->U[k][j][ie+i].Er < 0.0){
				pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie+i-1].Er;

			}

			
		}	
		else{

			Fr0x = 0.0;
			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie+i-1].Er;
			
		}

		

	
/*		Eratio = pGrid->U[k][j][ie+i].Edd_11 + 0.5 * pGrid->dx1 * reducefactor * Sigma_t;
*/
		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
/*		pGrid->U[k][j][ie+i].Er = (pGrid->U[k][j][ie+i-1].Edd_11 * pGrid->U[k][j][ie+i-1].Er - 0.5 * pGrid->dx1 * Sigma_t * Fr0x ) / Eratio;

		if((pGrid->U[k][j][ie+i].Er > pGrid->U[k][j][ie+i-1].Er) || (pGrid->U[k][j][ie+i].Er < 0.0)){
			Eratio = pGrid->U[k][j][ie+i].Edd_11 + pGrid->dx1 * reducefactor * Sigma_t;
			pGrid->U[k][j][ie+i].Er = pGrid->U[k][j][ie+i-1].Edd_11 * pGrid->U[k][j][ie+i-1].Er / Eratio;
		}	
			
		
		Fr0x = reducefactor * pGrid->U[k][j][ie+i].Er;
*/
		pGrid->U[k][j][ie+i].Fr1 = Fr0x + ((1.0 + pGrid->U[k][j][ie+i].Edd_11) * velocity_x1 + pGrid->U[k][j][ie+i].Edd_21 * velocity_y1 + pGrid->U[k][j][ie+i].Edd_31 * velocity_z1)* pGrid->U[k][j][ie+i].Er / Crat;
                pGrid->U[k][j][ie+i].Fr2 = Fr0y + ((1.0 + pGrid->U[k][j][ie+i].Edd_22) * velocity_y1 + pGrid->U[k][j][ie+i].Edd_21 * velocity_x1 + pGrid->U[k][j][ie+i].Edd_32 * velocity_z1)* pGrid->U[k][j][ie+i].Er / Crat;
		pGrid->U[k][j][ie+i].Fr3 = Fr0z + ((1.0 + pGrid->U[k][j][ie+i].Edd_33) * velocity_z1 + pGrid->U[k][j][ie+i].Edd_31 * velocity_x1 + pGrid->U[k][j][ie+i].Edd_32 * velocity_y1)* pGrid->U[k][j][ie+i].Er / Crat;


/*
		if(Fr0z > 0.0)
	                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;
		else
			pGrid->U[ke+k][j][i].Fr3 = ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;

*/			


		


    		}
		}
	}
  
}


/* Need to copy data from [k][js][i] to [k][j][is] */

void radMHD_periodis(GridS *pGrid)
{
  	int i, j, k, is, ie, js, je, ks, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;
#if defined(MHD) || defined(RADIATION_MHD)
  	int ju, ku; /* j-upper, k-upper */
#endif

#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif


	/* Do not need MPI communication */
	if(Periodix1_id < 0){
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {

        				pGrid->U[k][j][is-i] = pGrid->U[k][js+i-1][j-js+is];

					/* reflect velocity */
					pGrid->U[k][j][is-i].M2 = pGrid->U[k][js+i-1][j-js+is].M1;
					pGrid->U[k][j][is-i].M1 = -pGrid->U[k][js+i-1][j-js+is].M2;

#if defined(MHD) || defined(RADIATION_MHD)
					pGrid->U[k][j][is-i].B2c = pGrid->U[k][js+i-1][j-js+is].B1c;
					pGrid->U[k][j][is-i].B1c = -pGrid->U[k][js+i-1][j-js+is].B2c;
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost-1; i++) {
        				pGrid->B1i[k][j][is-i] = -pGrid->B2i[k][js+i-1][j-js+is];
      				}
    			}
  		}

  		if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  			for (k=ks; k<=ke; k++) {
    				for (j=js; j<=ju; j++) {
      					for (i=1; i<=nghost; i++) {
        					pGrid->B2i[k][j][is-i] = pGrid->B1i[k][js+i-1][j-js+is];
      					}
    				}
  			}

  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  			for (k=ks; k<=ku; k++) {
    				for (j=js; j<=je; j++) {
      					for (i=1; i<=nghost; i++) {
        					pGrid->B3i[k][j][is-i] = pGrid->B3i[k][js+i-1][j-js+is];
      					}
    				}
  			}
#endif /* MHD */


	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
		/* First, count the number of data expected to receive */
		cntR = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(NVAR);
#if defined(MHD) || defined(RADIATION_MHD)
    		cnt2 = (pGrid->Nx[1] > 1) ? (pGrid->Nx[1] + 1) : 1;
    		cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    		cntR += (nghost-1)*(pGrid->Nx[1])*(pGrid->Nx[2]);	/* add B1i */
    		cntR += nghost*cnt2*(pGrid->Nx[2]);			/* add B2i */
    		cntR += nghost*(pGrid->Nx[1])*cnt3;			/* add B3i */
#endif	

		/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");

		/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodix1_id,1,
			pGrid->Comm_Domain, &(recv_rq));


		
		/* Now count the number of data needs to send */
		cntS = (pGrid->Nx[1] + 2*nghost)*nghost*(pGrid->Nx[2])*(NVAR);
#if defined(MHD) || defined(RADIATION_MHD)
    		cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    		cntS += (pGrid->Nx[0] + 2*nghost - 1)*nghost*(pGrid->Nx[2]);
    		cntS += (pGrid->Nx[0] + 2*nghost)*(nghost-1)*(pGrid->Nx[2]);
    		cntS += (pGrid->Nx[0] + 2*nghost)*nghost*cnt3;
#endif
		/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");

		

		
		/*--------------------------------------------------------------*/
		/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
		/* send data in the order that will be needed */
		
  		for (k=ks; k<=ke; k++) {
			for (i=is+(nghost-1); i>=is; i--) {
    				for (j=js-nghost; j<=je+nghost; j++) {
					/* send the data according to the order that will be received */				

        				*(pSnd++) = pGrid->U[k][j][i].d;
        				*(pSnd++) = pGrid->U[k][j][i].M2;
        				*(pSnd++) = -pGrid->U[k][j][i].M1;
        				*(pSnd++) = pGrid->U[k][j][i].M3;
#ifndef BAROTROPIC
        				*(pSnd++) = pGrid->U[k][j][i].E;
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
        				*(pSnd++) = pGrid->U[k][j][i].B2c;
        				*(pSnd++) = -pGrid->U[k][j][i].B1c;
        				*(pSnd++) = pGrid->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        				for (n=0; n<NSCALARS; n++) *(pSnd++) = pGrid->U[k][j][i].s[n];
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
				/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
			for(i=is+(nghost-1); i>=is; i--){
				for(j=js-(nghost-1); j<=je+nghost; j++){			
        				*(pSnd++) = pGrid->B2i[k][j][i];
      				}
    			}
  		}

/* B2i at j=js maps to B2i at j=je+1 and is not passed */
  		for (k=ks; k<=ke; k++) {
			for(i=is+(nghost-1); i>=is+1; i--){
				for(j=js-nghost; j<=je+nghost; j++){
			  				
        				*(pSnd++) = -pGrid->B1i[k][j][i];
      				}
    			}
  		}

  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  		for (k=ks; k<=ku; k++) {
			for(i=is+(nghost-1); i>=is; i--){
				for(j=js-nghost; j<=je+nghost; j++){
        				*(pSnd++) = pGrid->B3i[k][j][i];
      				}
    			}
  		}
#endif /* MHD */



		/* now actually send the  data */
		 ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodix1_id,1,
        		pGrid->Comm_Domain, &(send_rq));

		 /* check non-blocking sends have completed. */
	      	ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);

		/*now get the receive data  */

		 /* check non-blocking receive have finished. */
     		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);

		pRcv = (double*)&(recv_buf[0]);

		/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
    			for (j=js; j<=je; j++){
      				for (i=is-nghost; i<=is-1; i++){
        				pGrid->U[k][j][i].d  = *(pRcv++);
        				pGrid->U[k][j][i].M1 = *(pRcv++);
        				pGrid->U[k][j][i].M2 = *(pRcv++);
        				pGrid->U[k][j][i].M3 = *(pRcv++);
#ifndef BAROTROPIC
        				pGrid->U[k][j][i].E  = *(pRcv++);
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
        				pGrid->U[k][j][i].B1c = *(pRcv++);
        				pGrid->U[k][j][i].B2c = *(pRcv++);
        				pGrid->U[k][j][i].B3c = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        				for (n=0; n<NSCALARS; n++) pGrid->U[k][j][i].s[n] = *(pRcv++);
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=is-(nghost-1); i<=is-1; i++){
        				pGrid->B1i[k][j][i] = *(pRcv++);
      				}
    			}
  		}

  		if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=ju; j++) {
      				for (i=is-nghost; i<=is-1; i++){
        				pGrid->B2i[k][j][i] = *(pRcv++);
      				}
    			}
  		}

  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  		for (k=ks; k<=ku; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=is-nghost; i<=is-1; i++){
        				pGrid->B3i[k][j][i] = *(pRcv++);
      				}
    			}
  		}
#endif /* MHD */



		/* Free the memory */
		free(recv_buf);
		free(send_buf);

#else /* if not MPI_PARALLEL */
		ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodix1_id);
#endif
	}/* End MPI case */


  
}



/* Need to copy data from [k][js][i] to [k][j][is] */

void radMHD_rad_periodis(GridS *pGrid)
{
  	int i, j, k, is, ie, js, je, ks, ke, m;
	int NRad = 10 + NOPACITY;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;
#if defined(MHD) || defined(RADIATION_MHD)
  	int ju, ku; /* j-upper, k-upper */
#endif

#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif


	/* Do not need MPI communication */
	if(Periodix1_id < 0){
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=nghost; i++) {

        				pGrid->U[k][j][is-i].Er = pGrid->U[k][js+i-1][j-js+is].Er;
					pGrid->U[k][j][is-i].Fr1 = -pGrid->U[k][js+i-1][j-js+is].Fr2;
					pGrid->U[k][j][is-i].Fr2 = pGrid->U[k][js+i-1][j-js+is].Fr1;
					pGrid->U[k][j][is-i].Fr3 = pGrid->U[k][js+i-1][j-js+is].Fr3;
					pGrid->U[k][j][is-i].Edd_11 = pGrid->U[k][js+i-1][j-js+is].Edd_22;
					pGrid->U[k][j][is-i].Edd_21 = pGrid->U[k][js+i-1][j-js+is].Edd_21;
					pGrid->U[k][j][is-i].Edd_22 = pGrid->U[k][js+i-1][j-js+is].Edd_11;
					pGrid->U[k][j][is-i].Edd_31 = -pGrid->U[k][js+i-1][j-js+is].Edd_32;
					pGrid->U[k][j][is-i].Edd_32 = pGrid->U[k][js+i-1][j-js+is].Edd_31;
					pGrid->U[k][j][is-i].Edd_33 = pGrid->U[k][js+i-1][j-js+is].Edd_33;
				
					for(m=0; m<NOPACITY; m++){
						pGrid->U[k][j][is-i].Sigma[m] = pGrid->U[k][js+i-1][j-js+is].Sigma[m];
					}
      				}
    			}
  		}


	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
		/* First, count the number of data expected to receive */
		cntR = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(NRad);

		/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");

		/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodix1_id,1,
			pGrid->Comm_Domain, &(recv_rq));
		
		/* Now count the number of data needs to send */
		cntS = (pGrid->Nx[1] + 2*nghost)*nghost*(pGrid->Nx[2])*(NRad);

		/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");

		/*--------------------------------------------------------------*/
		/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
		/* send data in the order that will be needed */
		
  		for (k=ks; k<=ke; k++) {
			for (i=is+(nghost-1); i>=is; i--) {
    				for (j=js-nghost; j<=je+nghost; j++) {
      				
        				*(pSnd++) = pGrid->U[k][j][i].Er;
        				*(pSnd++) = pGrid->U[k][j][i].Fr2;
        				*(pSnd++) = -pGrid->U[k][j][i].Fr1;
        				*(pSnd++) = pGrid->U[k][j][i].Fr3;
					*(pSnd++) = pGrid->U[k][j][i].Edd_22;
        				*(pSnd++) = pGrid->U[k][j][i].Edd_21;
        				*(pSnd++) = pGrid->U[k][j][i].Edd_11;
        				*(pSnd++) = pGrid->U[k][j][i].Edd_32;
					*(pSnd++) = -pGrid->U[k][j][i].Edd_31;
        				*(pSnd++) = pGrid->U[k][j][i].Edd_33;
        				for(m=0; m<NOPACITY; m++){
						*(pSnd++) = pGrid->U[k][j][i].Sigma[m];
					}



      				}
    			}
  		}

		/* now actually send the  data */
		 ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodix1_id,1,
        		pGrid->Comm_Domain, &(send_rq));

		 /* check non-blocking sends have completed. */
	      	ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);

		/*now get the receive data  */

		 /* check non-blocking receive have finished. */
     		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);

		pRcv = (double*)&(recv_buf[0]);

		/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
    			for (j=js; j<=je; j++){
      				for (i=is-nghost; i<=is-1; i++){
        				pGrid->U[k][j][i].Er  = *(pRcv++);
        				pGrid->U[k][j][i].Fr1 = *(pRcv++);
        				pGrid->U[k][j][i].Fr2 = *(pRcv++);
        				pGrid->U[k][j][i].Fr3 = *(pRcv++);
					pGrid->U[k][j][i].Edd_11 = *(pRcv++);
        				pGrid->U[k][j][i].Edd_21 = *(pRcv++);
        				pGrid->U[k][j][i].Edd_22 = *(pRcv++);
        				pGrid->U[k][j][i].Edd_31 = *(pRcv++);
					pGrid->U[k][j][i].Edd_32 = *(pRcv++);
        				pGrid->U[k][j][i].Edd_33 = *(pRcv++);

					for(m=0; m<NOPACITY; m++){
						pGrid->U[k][j][i].Sigma[m] = *(pRcv++);
					}

      				}
    			}
  		}



		/* Free the memory */
		free(recv_buf);
		free(send_buf);

#else /* if not MPI_PARALLEL */
		ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodix1_id);
#endif
	}/* End MPI case */


  
}


/* y boundary */

void radMHD_periodjs(GridS *pGrid)
{
  	int i, j, k, is, ie, js, je, ks, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;
#if defined(MHD) || defined(RADIATION_MHD)
  	int ju, ku; /* j-upper, k-upper */
#endif

#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif


	/* Do not need MPI communication */
	if(Periodjx1_id < 0){
		
		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is; i<=ie; i++) {
        				pGrid->U[k][js-j][i] = pGrid->U[k][i-is+js][is+j-1];

					/* reflect velocity */
					pGrid->U[k][js-j][i].M2 = -pGrid->U[k][i-is+js][is+j-1].M1;
					pGrid->U[k][js-j][i].M1 =  pGrid->U[k][i-is+js][is+j-1].M2;
#if defined(MHD) || defined(RADIATION_MHD)
					pGrid->U[k][js-j][i].B2c = -pGrid->U[k][i-is+js][is+j-1].B1c;
					pGrid->U[k][js-j][i].B1c =  pGrid->U[k][i-is+js][is+j-1].B2c;
#endif
      				}
    			}
  		}

  /* Now the corner */
  /* The corner will copy some data in the just updated ghost zones */
 		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<is; i++) {
        				pGrid->U[k][js-j][i] = pGrid->U[k][i-is+js][is+j-1];

					pGrid->U[k][js-j][i].M2 = -pGrid->U[k][i-is+js][is+j-1].M1;
					pGrid->U[k][js-j][i].M1 =  pGrid->U[k][i-is+js][is+j-1].M2;
#if defined(MHD) || defined(RADIATION_MHD)
					pGrid->U[k][js-j][i].B2c = -pGrid->U[k][i-is+js][is+j-1].B1c;
					pGrid->U[k][js-j][i].B1c =  pGrid->U[k][i-is+js][is+j-1].B2c;
#endif
      				}
    			}
  		}

 		for (k=ks; k<=ke; k++) {	
    			for (j=1; j<=nghost; j++) {
      				for (i=ie+1; i<=ie+nghost; i++) {
        				pGrid->U[k][js-j][i] = pGrid->U[k][i-is+js][is+j-1];

					pGrid->U[k][js-j][i].M2 = -pGrid->U[k][i-is+js][is+j-1].M1;
					pGrid->U[k][js-j][i].M1 =  pGrid->U[k][i-is+js][is+j-1].M2;

#if defined(MHD) || defined(RADIATION_MHD)
					pGrid->U[k][js-j][i].B2c = -pGrid->U[k][i-is+js][is+j-1].B1c;
					pGrid->U[k][js-j][i].B1c =  pGrid->U[k][i-is+js][is+j-1].B2c;
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-(nghost-1); i<=ie+nghost; i++) {
        				pGrid->B1i[k][js-j][i] = pGrid->B2i[k][i-is+js][is+j-1];
      				}
    			}
  		}

  /* Now the corner */
   		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-(nghost-1); i<is; i++) {
        				pGrid->B1i[k][js-j][i] = pGrid->B2i[k][i-is+js][is+j-1];
      				}
    			}
  		}

  		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=ie+1; i<=ie+nghost; i++) {
        				pGrid->B1i[k][js-j][i] = pGrid->B2i[k][i-is+js][is+j-1];
      				}
    			}
  		}

/* B2i is not set at j=js-nghost */
  		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost-1; j++) {
      				for (i=is; i<=ie; i++) {
        				pGrid->B2i[k][js-j][i] = -pGrid->B1i[k][i-is+js][is+j-1];
      				}
    			}
  		}

 
  /* Now the corner */

  		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost-1; j++) {
      				for (i=is-nghost; i<is; i++) {
        				pGrid->B2i[k][js-j][i] = -pGrid->B1i[k][i-is+js][is+j-1];
      				}
    			}
  		}

 		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost-1; j++) {
      				for (i=ie+1; i<=ie+nghost; i++) {
        				pGrid->B2i[k][js-j][i] = -pGrid->B1i[k][i-is+js][is+j-1];
      				}
    			}
  		}


  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  		for (k=ks; k<=ku; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is; i<=ie; i++) {
        				pGrid->B3i[k][js-j][i] = pGrid->B3i[k][i-is+js][is+j-1];
      				}
    			}
  		}

  
  		/* Now the corner */
  		for (k=ks; k<=ku; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<is; i++) {
        				pGrid->B3i[k][js-j][i] = pGrid->B3i[k][i-is+js][is+j-1];
      				}
    			}
  		}

  		for (k=ks; k<=ku; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=ie+1; i<=ie+nghost; i++) {
        				pGrid->B3i[k][js-j][i] = pGrid->B3i[k][i-is+js][is+j-1];
      				}
    			}
  		}



#endif /* MHD */


	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
		/* Now count the number of data needs to receive */
		cntR = (pGrid->Nx[1] + 2*nghost)*nghost*(pGrid->Nx[2])*(NVAR);
#if defined(MHD) || defined(RADIATION_MHD)
    		cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    		cntR += (pGrid->Nx[0] + 2*nghost - 1)*nghost*(pGrid->Nx[2]);
    		cntR += (pGrid->Nx[0] + 2*nghost)*(nghost-1)*(pGrid->Nx[2]);
    		cntR += (pGrid->Nx[0] + 2*nghost)*nghost*cnt3;
#endif

		/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");

		/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodjx1_id,1,
			pGrid->Comm_Domain, &(recv_rq));



		/* Count the number of data needs to send */
		cntS = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(NVAR);
#if defined(MHD) || defined(RADIATION_MHD)
    		cnt2 = (pGrid->Nx[1] > 1) ? (pGrid->Nx[1] + 1) : 1;
    		cnt3 = (pGrid->Nx[2] > 1) ? (pGrid->Nx[2] + 1) : 1;
    		cntS += (nghost-1)*(pGrid->Nx[1])*(pGrid->Nx[2]);	/* add B1i */
    		cntS += nghost*cnt2*(pGrid->Nx[2]);			/* add B2i */
    		cntS += nghost*(pGrid->Nx[1])*cnt3;			/* add B3i */
#endif	

		/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");

		

		
		/*--------------------------------------------------------------*/
		/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
		/* send data in the order that will be needed */
		
  		for (k=ks; k<=ke; k++) {			
			for(i=is; i<=ie; i++){    				
				for(j=js+(nghost-1); j>=js; j--){
			
        				*(pSnd++) = pGrid->U[k][j][i].d;
        				*(pSnd++) = -pGrid->U[k][j][i].M2;
        				*(pSnd++) = pGrid->U[k][j][i].M1;
        				*(pSnd++) = pGrid->U[k][j][i].M3;
#ifndef BAROTROPIC
        				*(pSnd++) = pGrid->U[k][j][i].E;
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
        				*(pSnd++) = -pGrid->U[k][j][i].B2c;
        				*(pSnd++) = pGrid->U[k][j][i].B1c;
        				*(pSnd++) = pGrid->U[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
        				for (n=0; n<NSCALARS; n++) *(pSnd++) = pGrid->U[k][j][i].s[n];
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
				/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
			for(i=is; i<=ie; i++){				
				for(j=js+(nghost-2);j>=js; j--){			
        				*(pSnd++) = -pGrid->B2i[k][j][i];
      				}
    			}
  		}

/* B2i at j=js maps to B2i at j=je+1 and is not passed */
  		for (k=ks; k<=ke; k++) {
			for(i=is; i<=ie+1; i++){		
				for(j=js+(nghost-1); j>=js; j--){		
			  				
        				*(pSnd++) = pGrid->B1i[k][j][i];
      				}
    			}
  		}

  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  		for (k=ks; k<=ku; k++) {
			for(i=is; i<=ie; i++){			
				for(j=js+(nghost-1); j>=js; j--){				
        				*(pSnd++) = pGrid->B3i[k][j][i];
      				}
    			}
  		}
#endif /* MHD */



		/* now actually send the  data */
		 ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodjx1_id,1,
        		pGrid->Comm_Domain, &(send_rq));

		 /* check non-blocking sends have completed. */
	      	ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);

		/*now get the receive data  */

		 /* check non-blocking receive have finished. */
     		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);

		pRcv = (double*)&(recv_buf[0]);

		/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
    			for (j=js-nghost; j<js; j++){
      				for (i=is-nghost; i<=ie+nghost; i++){
        				pGrid->U[k][j][i].d  = *(pRcv++);
        				pGrid->U[k][j][i].M1 = *(pRcv++);
        				pGrid->U[k][j][i].M2 = *(pRcv++);
        				pGrid->U[k][j][i].M3 = *(pRcv++);
#ifndef BAROTROPIC
        				pGrid->U[k][j][i].E  = *(pRcv++);
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
        				pGrid->U[k][j][i].B1c = *(pRcv++);
        				pGrid->U[k][j][i].B2c = *(pRcv++);
        				pGrid->U[k][j][i].B3c = *(pRcv++);
#endif /* MHD */
#if (NSCALARS > 0)
        				for (n=0; n<NSCALARS; n++) pGrid->U[k][j][i].s[n] = *(pRcv++);
#endif
      				}
    			}
  		}

#if defined(MHD) || defined(RADIATION_MHD)
/* B1i is not set at i=is-nghost */
  		for (k=ks; k<=ke; k++) {
    			for (j=js-nghost; j<js; j++) {
      				for (i=is-(nghost-1); i<=ie+nghost; i++){
        				pGrid->B1i[k][j][i] = *(pRcv++);
      				}
    			}
  		}

  		
  		for (k=ks; k<=ke; k++) {
    			for (j=js-(nghost-1); j<js; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++){
        				pGrid->B2i[k][j][i] = *(pRcv++);
      				}
    			}
  		}

  		if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  		for (k=ks; k<=ku; k++) {
    			for (j=js-nghost; j<js; j++) {
      				for (i=is-nghost; i<=ie+nghost; i++){
        				pGrid->B3i[k][j][i] = *(pRcv++);
      				}
    			}
  		}
#endif /* MHD */



		/* Free the memory */
		free(recv_buf);
		free(send_buf);

#else /* if not MPI_PARALLEL */
		ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodjx1_id);
#endif
	}/* End MPI case */


  
}

												  
												  
void radMHD_rad_periodjs(GridS *pGrid)
{
	int i, j, k, is, ie, js, je, ks, ke, m;
	int NRad = 10 + NOPACITY;
	je = pGrid->je;
	js = pGrid->js;
	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;
#if defined(MHD) || defined(RADIATION_MHD)
	int ju, ku; /* j-upper, k-upper */
#endif
												  
#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif
												  
												  
	/* Do not need MPI communication */
	if(Periodjx1_id < 0){												  
		for (k=ks; k<=ke; k++) {
			for (j=1; j<=nghost; j++) {
				for (i=is; i<=ie; i++) {
					pGrid->U[k][js-j][i].Er = pGrid->U[k][i-is+js][is+j-1].Er;
					pGrid->U[k][js-j][i].Fr1 = pGrid->U[k][i-is+js][is+j-1].Fr2;
					pGrid->U[k][js-j][i].Fr2 = -pGrid->U[k][i-is+js][is+j-1].Fr1;
					pGrid->U[k][js-j][i].Fr3 = pGrid->U[k][i-is+js][is+j-1].Fr3;
					pGrid->U[k][js-j][i].Edd_11 = pGrid->U[k][i-is+js][is+j-1].Edd_22;
					pGrid->U[k][js-j][i].Edd_21 = pGrid->U[k][i-is+js][is+j-1].Edd_21;
					pGrid->U[k][js-j][i].Edd_22 = pGrid->U[k][i-is+js][is+j-1].Edd_11;
					pGrid->U[k][js-j][i].Edd_31 = pGrid->U[k][i-is+js][is+j-1].Edd_32;
					pGrid->U[k][js-j][i].Edd_32 = -pGrid->U[k][i-is+js][is+j-1].Edd_31;
					pGrid->U[k][js-j][i].Edd_33 = pGrid->U[k][i-is+js][is+j-1].Edd_33;
					for(m=0; m<NOPACITY; m++){
						pGrid->U[k][js-j][i].Sigma[m] = pGrid->U[k][i-is+js][is+j-1].Sigma[m];						
												  
					}
				}
			}
		}
			
		 /* Now the corner */
  /* The corner will copy some data in the just updated ghost zones */
 		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=nghost; j++) {
      				for (i=is-nghost; i<is; i++) {
        				pGrid->U[k][js-j][i].Er = pGrid->U[k][i-is+js][is+j-1].Er;
					pGrid->U[k][js-j][i].Fr1 = pGrid->U[k][i-is+js][is+j-1].Fr2;
					pGrid->U[k][js-j][i].Fr2 = -pGrid->U[k][i-is+js][is+j-1].Fr1;
					pGrid->U[k][js-j][i].Fr3 = pGrid->U[k][i-is+js][is+j-1].Fr3;
					pGrid->U[k][js-j][i].Edd_11 = pGrid->U[k][i-is+js][is+j-1].Edd_22;
					pGrid->U[k][js-j][i].Edd_21 = pGrid->U[k][i-is+js][is+j-1].Edd_21;
					pGrid->U[k][js-j][i].Edd_22 = pGrid->U[k][i-is+js][is+j-1].Edd_11;
					pGrid->U[k][js-j][i].Edd_31 = pGrid->U[k][i-is+js][is+j-1].Edd_32;
					pGrid->U[k][js-j][i].Edd_32 = -pGrid->U[k][i-is+js][is+j-1].Edd_31;
					pGrid->U[k][js-j][i].Edd_33 = pGrid->U[k][i-is+js][is+j-1].Edd_33;
					for(m=0; m<NOPACITY; m++){
						pGrid->U[k][js-j][i].Sigma[m] = pGrid->U[k][i-is+js][is+j-1].Sigma[m];						
												  
					}
      				}
    			}
  		}

 		for (k=ks; k<=ke; k++) {	
    			for (j=1; j<=nghost; j++) {
      				for (i=ie+1; i<=ie+nghost; i++) {
        				pGrid->U[k][js-j][i].Er = pGrid->U[k][i-is+js][is+j-1].Er;
					pGrid->U[k][js-j][i].Fr1 = pGrid->U[k][i-is+js][is+j-1].Fr2;
					pGrid->U[k][js-j][i].Fr2 = -pGrid->U[k][i-is+js][is+j-1].Fr1;
					pGrid->U[k][js-j][i].Fr3 = pGrid->U[k][i-is+js][is+j-1].Fr3;
					pGrid->U[k][js-j][i].Edd_11 = pGrid->U[k][i-is+js][is+j-1].Edd_22;
					pGrid->U[k][js-j][i].Edd_21 = pGrid->U[k][i-is+js][is+j-1].Edd_21;
					pGrid->U[k][js-j][i].Edd_22 = pGrid->U[k][i-is+js][is+j-1].Edd_11;
					pGrid->U[k][js-j][i].Edd_31 = pGrid->U[k][i-is+js][is+j-1].Edd_32;
					pGrid->U[k][js-j][i].Edd_32 = -pGrid->U[k][i-is+js][is+j-1].Edd_31;
					pGrid->U[k][js-j][i].Edd_33 = pGrid->U[k][i-is+js][is+j-1].Edd_33;
					for(m=0; m<NOPACITY; m++){
						pGrid->U[k][js-j][i].Sigma[m] = pGrid->U[k][i-is+js][is+j-1].Sigma[m];						
												  
					}
      				}
    			}
  		}

									  
												  
	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
	/* Now count the number of data needs to receive */
		cntR = (pGrid->Nx[1] + 2*nghost)*nghost*(pGrid->Nx[2])*(NRad);

												  
/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");
												  
	/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodjx1_id,1,
				 pGrid->Comm_Domain, &(recv_rq));
												  
												  
												  
	/* Count the number of data needs to send */
		cntS = nghost*(pGrid->Nx[1])*(pGrid->Nx[2])*(NRad);

												  
	/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");
											  
												  
 /*--------------------------------------------------------------*/
/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
/* send data in the order that will be needed */
												  
		for (k=ks; k<=ke; k++) {			
			for(i=is; i<=ie; i++){    				
				for(j=js+(nghost-1); j>=js; j--){
												  
					*(pSnd++) = pGrid->U[k][j][i].Er;
					*(pSnd++) = -pGrid->U[k][j][i].Fr2;
					*(pSnd++) = pGrid->U[k][j][i].Fr1;
					*(pSnd++) = pGrid->U[k][j][i].Fr3;
					*(pSnd++) = pGrid->U[k][j][i].Edd_22;
					*(pSnd++) = pGrid->U[k][j][i].Edd_21;
				    	*(pSnd++) = pGrid->U[k][j][i].Edd_11;
					*(pSnd++) = -pGrid->U[k][j][i].Edd_32;
					*(pSnd++) = pGrid->U[k][j][i].Edd_31;
					*(pSnd++) = pGrid->U[k][j][i].Edd_33;
					for(m=0; m<NOPACITY; m++){
							*(pSnd++) = pGrid->U[k][j][i].Sigma[m];
					}
												  


				}
			}
		}

																			
																			
																			
/* now actually send the  data */
		ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodjx1_id,1,
				pGrid->Comm_Domain, &(send_rq));
																			
/* check non-blocking sends have completed. */
		ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);
																			
/*now get the receive data  */
																			
/* check non-blocking receive have finished. */
		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);
																			
		pRcv = (double*)&(recv_buf[0]);
																			
/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
			for (j=js-nghost; j<js; j++){
				for (i=is-nghost; i<=ie+nghost; i++){
					pGrid->U[k][j][i].Er  = *(pRcv++);
					pGrid->U[k][j][i].Fr1 = *(pRcv++);
					pGrid->U[k][j][i].Fr2 = *(pRcv++);
					pGrid->U[k][j][i].Fr3 = *(pRcv++);
				    	pGrid->U[k][j][i].Edd_11 = *(pRcv++);
					pGrid->U[k][j][i].Edd_21 = *(pRcv++);
				    	pGrid->U[k][j][i].Edd_22 = *(pRcv++);
					pGrid->U[k][j][i].Edd_31 = *(pRcv++);
					pGrid->U[k][j][i].Edd_32 = *(pRcv++);
				    	pGrid->U[k][j][i].Edd_33 = *(pRcv++);
												  
				    for(m=0; m<NOPACITY; m++){
					   pGrid->U[k][j][i].Sigma[m] = *(pRcv++);
					}												  

				}
			}
		}
																			
																			
																			
	/* Free the memory */
	free(recv_buf);
	free(send_buf);
																			
#else /* if not MPI_PARALLEL */
	ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodjx1_id);
#endif
	}/* End MPI case */										
}
																			
																	
void radMHD_outflowje(GridS *pGrid)
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
	Real Sigma[NOPACITY];
	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;

	static Real x3t;
	x3t = ztop - 0.5 * pGrid->dx3;
	PrimS Wopacity;
	Real distance, Radius, Vc;

#if defined(MHD) || defined(RADIATION_MHD)
	for(k=ks-nghost; k<=ke+nghost; k++){
		for(i=is-nghost; i<=ie+nghost; i++){
			for(j=1; j<=nghost; j++){
		
			if(pGrid->U[k][je][i].M2 < 0.0){
				pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-1][i];
                                pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je+j-1][i];
                        }
                        else{
                                pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-1][i];
                                pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je+j-1][i];
                        }


				if(j<nghost)
					pGrid->B2i[k][je+j+1][i] = pGrid->B2i[k][je+j][i];


			}
		}
	}


	/* Now update the cell centered values */
	for(k=ks-nghost; k<=ke+nghost; k++){
		for(i=is-nghost; i<=ie+nghost; i++){
			for(j=1; j<=nghost; j++){
				 if(pGrid->U[k][je][i].M2 < 0.0){
                                         pGrid->U[k][je+j][i].B3c = pGrid->U[k][je+j-1][i].B3c;
                                         pGrid->U[k][je+j][i].B1c = pGrid->U[k][je+j-1][i].B1c;
                                 }
                                 else{
                                         pGrid->U[k][je+j][i].B3c = pGrid->U[k][je+j-1][i].B3c;
                                         pGrid->U[k][je+j][i].B1c = pGrid->U[k][je+j-1][i].B1c;
                                 }
				
				if(j<nghost){
					pGrid->U[k][je+j][i].B2c = 0.5 * (pGrid->B2i[k][je+j][i] + pGrid->B2i[k][je+j+1][i]);
				}
				else{
					pGrid->U[k][je+j][i].B2c = pGrid->B2i[k][je+j][i];
				}

			}
		}
	}
		



#endif /* MHD */

	for(k=ks-nghost; k<=ke+nghost; k++){
		for(i=is-nghost; i<=ie+nghost; i++){
			for(j=1; j<=nghost; j++){
	    
	

		cc_pos(pGrid, i, je+j,k, &x1, &x2, &x3);
#if defined(MHD) || defined(RADIATION_MHD)
		Bx = pGrid->U[k][je+j][i].B1c;
		By = pGrid->U[k][je+j][i].B2c;
		Bz = pGrid->U[k][je+j][i].B3c;
		
		B = sqrt(Bx * Bx + By * By + Bz * Bz);
#endif

		if(pGrid->U[k][je][i].d < TINY_NUMBER)
			pGrid->U[k][je][i].d = pGrid->U[4][je][4].d;

		
		distance = sqrt(x1 * x1 + x2 * x2 + x3 * x3);
		Radius = sqrt(x1 * x1 + x2 * x2);
		Vc = Crat*sqrt(1.0/(2.0*distance))/(distance-1.0);

	
		density  = pGrid->U[k][je][i].d;

		velocity1 = pGrid->U[k][je][i].M1 / pGrid->U[k][je][i].d;
		velocity2 = pGrid->U[k][je][i].M2 / pGrid->U[k][je][i].d;
		velocity3 = pGrid->U[k][je][i].M3 / pGrid->U[k][je][i].d;

		pressure = (pGrid->U[k][je][i].E - 0.5 * density * (velocity1 * velocity1 + velocity2 * velocity2 + velocity3 * velocity3)) * (Gamma - 1.0);
#ifdef RADIATION_MHD

		pressure -= 0.5 * (pGrid->U[k][je][i].B1c * pGrid->U[k][je][i].B1c + pGrid->U[k][je][i].B2c * pGrid->U[k][je][i].B2c + pGrid->U[k][je][i].B3c * pGrid->U[k][je][i].B3c) * (Gamma - 1.0);
#endif

		T = pressure / (pGrid->U[k][je][i].d * R_ideal);

		/* Now reset the velocity */
		velocity1 = Vc * x2;
		velocity2 = -Vc * x1;

/*
		if(T > Tfloor)
			T = Tfloor;
*/
		/* Tfloor is usually set as initial temperature in the ghost zoner */

		/* Now extrapolate the density to balance gravity assuming a constant temperature in the ghost zones */
		
	
		pGrid->U[k][je+j][i].d = pGrid->U[k][je][i].d;



		
      		pGrid->U[k][je+j][i].M1 = velocity1 * pGrid->U[k][je+j][i].d;
		pGrid->U[k][je+j][i].M2 = velocity2 * pGrid->U[k][je+j][i].d;
			
		
		pGrid->U[k][je+j][i].M3 = velocity3 * pGrid->U[k][je+j][i].d;

		pressure = T * R_ideal * pGrid->U[k][je+j][i].d;

		pGrid->U[k][je+j][i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[k][je+j][i].M1 * pGrid->U[k][je+j][i].M1 + pGrid->U[k][je+j][i].M2 * pGrid->U[k][je+j][i].M2 + pGrid->U[k][je+j][i].M3 * pGrid->U[k][je+j][i].M3) / pGrid->U[k][je+j][i].d;



#ifdef RADIATION_MHD
		pGrid->U[k][je+j][i].E += 0.5 * (pGrid->U[k][je+j][i].B1c * pGrid->U[k][je+j][i].B1c + pGrid->U[k][je+j][i].B2c * pGrid->U[k][je+j][i].B2c + pGrid->U[k][je+j][i].B3c * pGrid->U[k][je+j][i].B3c);

#endif

		Wopacity = Cons_to_Prim(&pGrid->U[k][je+j][i]);
		/* Add background shearing */
#ifdef FARGO	
		cc_pos(pGrid,i,je+j,k,&x1,&x2,&x3);
		Wopacity.V2 -= qshear * Omega_0 * x1;		
#endif
			
		Thindiskopacity(&Wopacity, Sigma, NULL);
		
		if(T < 2.0 * TINY_NUMBER)
			Sigma[0] = kappaes * pGrid->U[k][je+j][i].d;

		for(m=0; m<NOPACITY;m++){
			pGrid->U[k][je+j][i].Sigma[m] = Sigma[m];
		}	
	
      		
      		}
    		}
	}
  
}


void radMHD_rad_outflowje(GridS *pGrid)
{
  	int i,j,js,je;
	int ks, is,ie,k, ke;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	Real Fr0x, Fr0y, Fr0z;
	Real velocity_x, velocity_y, velocity_z;
	Real velocity_x1, velocity_y1, velocity_z1;
	Real Sigma_t, Sigma_t1;
	Real Eratio, reducefactor, tau;


	Real x1, x2, x3;
	
	for(k=ks-nghost; k<=ke+nghost; k++){
		for(i=is-nghost; i<=ie+nghost; i++){
			for(j=1; j<=nghost; j++){
	
		cc_pos(pGrid, i, je+j, k, &x1, &x2, &x3);

		pGrid->U[k][je+j][i].Edd_11 = pGrid->U[k][je][i].Edd_11;
		pGrid->U[k][je+j][i].Edd_22 = pGrid->U[k][je][i].Edd_22;
		pGrid->U[k][je+j][i].Edd_21 = pGrid->U[k][je][i].Edd_21;
		pGrid->U[k][je+j][i].Edd_31 = pGrid->U[k][je][i].Edd_31;
		pGrid->U[k][je+j][i].Edd_32 = pGrid->U[k][je][i].Edd_32;
		pGrid->U[k][je+j][i].Edd_33 = pGrid->U[k][je][i].Edd_33;
	
		matrix_alpha(0.0, pGrid->U[k][je+j][i].Sigma, pGrid->dt, pGrid->U[k][je+j][i].Edd_22, 0.0, &reducefactor, 0, dy);
	
	/*
		tau = pGrid->dt * Crat * (pGrid->U[ke+k][j][i].Sigma[0] + pGrid->U[ke+k][j][i].Sigma[1]);
		tau = tau * tau / (2.0 * pGrid->U[ke+k][j][i].Edd_33);

		if(tau > 0.001)
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
		else
			reducefactor = sqrt(pGrid->U[ke+k][j][i].Edd_33 * (1.0 - 0.5 * tau));			
	*/
		Sigma_t = 0.5 * (pGrid->U[k][je+j][i].Sigma[0] + pGrid->U[k][je+j][i].Sigma[1] + pGrid->U[k][je+j-1][i].Sigma[0] + pGrid->U[k][je+j-1][i].Sigma[1]);
		Sigma_t1 = 0.5 * (pGrid->U[k][je+j-1][i].Sigma[0] + pGrid->U[k][je+j-1][i].Sigma[1] + pGrid->U[k][je+j-2][i].Sigma[0] + pGrid->U[k][je+j-2][i].Sigma[1]);
	
		velocity_x = pGrid->U[k][je+j-1][i].M1 / pGrid->U[k][je+j-1][i].d;
		velocity_y = pGrid->U[k][je+j-1][i].M2 / pGrid->U[k][je+j-1][i].d;
#ifdef FARGO
		velocity_y -= qshear * Omega_0 * x1;

#endif
		velocity_z = pGrid->U[k][je+j-1][i].M3 / pGrid->U[k][je+j-1][i].d;


		
		Fr0x = pGrid->U[k][je+j-1][i].Fr1 - ((1.0 + pGrid->U[k][je+j-1][i].Edd_11) * velocity_x + pGrid->U[k][je+j-1][i].Edd_21 * velocity_y + pGrid->U[k][je+j-1][i].Edd_31 * velocity_z)* pGrid->U[k][je+j-1][i].Er / Crat;

		Fr0y = pGrid->U[k][je+j-1][i].Fr2 - ((1.0 + pGrid->U[k][je+j-1][i].Edd_22) * velocity_y + pGrid->U[k][je+j-1][i].Edd_21 * velocity_x + pGrid->U[k][je+j-1][i].Edd_32 * velocity_z)* pGrid->U[k][je+j-1][i].Er / Crat;

		Fr0z = pGrid->U[k][je+j-1][i].Fr3 - ((1.0 + pGrid->U[k][je+j-1][i].Edd_33) * velocity_z + pGrid->U[k][je+j-1][i].Edd_31 * velocity_x + pGrid->U[k][je+j-1][i].Edd_32 * velocity_y)* pGrid->U[k][j][ie+i-1].Er / Crat;

		velocity_x1 = pGrid->U[k][je+j][i].M1 / pGrid->U[k][je+j][i].d;
		velocity_y1 = pGrid->U[k][je+j][i].M2 / pGrid->U[k][je+j][i].d;
 		velocity_z1 = pGrid->U[k][je+j][i].M3 / pGrid->U[k][je+j][i].d;
#ifdef FARGO
		velocity_y1 -= qshear * Omega_0 * x1;

#endif

/*
		pGrid->U[ke+k][j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ke+k][j][i].Edd_11) * velocity_x1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_y1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ke+k][j][i].Edd_22) * velocity_y1 + pGrid->U[ke+k][j][i].Edd_21 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_z1)* pGrid->U[ke+k][j][i].Er / Crat;
                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;



		pGrid->U[ke+k][j][i].Fr1 = pGrid->U[ke][j][i].Fr1;
		pGrid->U[ke+k][j][i].Fr2 = pGrid->U[ke][j][i].Fr2;
		if(pGrid->U[ke][j][i].Fr3 > 0.0)
			pGrid->U[ke+k][j][i].Fr3 = pGrid->U[ke][j][i].Fr3;
		else
			pGrid->U[ke+k][j][i].Fr3 = 0.0;
*/

		if(Fr0y > 0.0){			
			pGrid->U[k][je+j][i].Er = pGrid->U[k][je+j-1][i].Er - pGrid->dx2 * Sigma_t * Fr0y / pGrid->U[k][je+j-1][i].Edd_22;			

			if(pGrid->U[k][je+j][i].Er < 0.0){
				pGrid->U[k][je+j][i].Er = pGrid->U[k][je+j-1][i].Er;

			}

			
		}	
		else{
			Fr0y = 0.0;
			pGrid->U[k][je+j][i].Er = pGrid->U[k][je+j-1][i].Er;
			
		}


/*	
		Eratio = pGrid->U[k][je+j][i].Edd_22 + 0.5 * pGrid->dx2 * reducefactor * Sigma_t;
*/
		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
/*		pGrid->U[k][je+j][i].Er = (pGrid->U[k][je+j-1][i].Edd_22 * pGrid->U[k][je+j-1][i].Er - 0.5 * pGrid->dx2 * Sigma_t * Fr0y ) / Eratio;

		if((pGrid->U[k][je+j][i].Er > pGrid->U[k][je+j-1][i].Er) || (pGrid->U[k][je+j][i].Er < 0.0)){
			Eratio = pGrid->U[k][je+j][i].Edd_22 + pGrid->dx2 * reducefactor * Sigma_t;
			pGrid->U[k][je+j][i].Er = pGrid->U[k][je+j-1][i].Edd_22 * pGrid->U[k][je+j-1][i].Er / Eratio;
		}	
			
		
		Fr0y = reducefactor * pGrid->U[k][je+j][i].Er;
*/
		pGrid->U[k][je+j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[k][je+j][i].Edd_11) * velocity_x1 + pGrid->U[k][je+j][i].Edd_21 * velocity_y1 + pGrid->U[k][je+j][i].Edd_31 * velocity_z1)* pGrid->U[k][je+j][i].Er / Crat;
                pGrid->U[k][je+j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[k][je+j][i].Edd_22) * velocity_y1 + pGrid->U[k][je+j][i].Edd_21 * velocity_x1 + pGrid->U[k][je+j][i].Edd_32 * velocity_z1)* pGrid->U[k][je+j][i].Er / Crat;
		pGrid->U[k][je+j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[k][je+j][i].Edd_33) * velocity_z1 + pGrid->U[k][je+j][i].Edd_31 * velocity_x1 + pGrid->U[k][je+j][i].Edd_32 * velocity_y1)* pGrid->U[k][je+j][i].Er / Crat;


/*
		if(Fr0z > 0.0)
	                pGrid->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;
		else
			pGrid->U[ke+k][j][i].Fr3 = ((1.0 + pGrid->U[ke+k][j][i].Edd_33) * velocity_z1 + pGrid->U[ke+k][j][i].Edd_31 * velocity_x1 + pGrid->U[ke+k][j][i].Edd_32 * velocity_y1)* pGrid->U[ke+k][j][i].Er / Crat;

*/			


		


    		}
		}
	}
  
}

		


#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
#ifdef MATRIX_MULTIGRID

void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_periodis;
} 
void bvals_mat_fun_ox1(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_outflowie;
} 
void bvals_mat_fun_ix2(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_periodjs;
} 
void bvals_mat_fun_ox2(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_outflowje;
} 
void bvals_mat_fun_ix3(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_outflowks;
} 
void bvals_mat_fun_ox3(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_outflowke;
} 




void radMHD_Mat_outflowke(MatrixS *pMat)
{
  	int i, j, k, je, ju, jl, js, ks, ke;
	int is,ie, il, iu;
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;

	ju = pMat->je + Matghost;
  	jl = pMat->js - Matghost;

	il = pMat->is - Matghost;
	iu = pMat->ie + Matghost;

	Real x1, x2, x3;
	Real vx, vy, vz, Fr0z, Sigma_t;
	Real reducefactor, Eratio, tau;

	for (k=1;  k<=Matghost;  k++) {
		for(j=js-Matghost; j<=je+Matghost; j++){
			for(i=is-Matghost; i<=ie+Matghost; i++){
	
			/*	pMat->Ugas[ke+k][j][i].Edd_11 = pMat->Ugas[ke][j][i].Edd_11;
				pMat->Ugas[ke+k][j][i].Edd_22 = pMat->Ugas[ke][j][i].Edd_22;
				pMat->Ugas[ke+k][j][i].Edd_21 = pMat->Ugas[ke][j][i].Edd_21;
				pMat->Ugas[ke+k][j][i].Edd_31 = pMat->Ugas[ke][j][i].Edd_31;
				pMat->Ugas[ke+k][j][i].Edd_32 = pMat->Ugas[ke][j][i].Edd_32;
				pMat->Ugas[ke+k][j][i].Edd_33 = pMat->Ugas[ke][j][i].Edd_33;
			*/
				pMat->Ugas[ke+k][j][i] = pMat->Ugas[ke][j][i];
			
				pMat->U[ke+k][j][i].Er  = 0.0;
				pMat->U[ke+k][j][i].Fr1 = 0.0;
				pMat->U[ke+k][j][i].Fr2 = 0.0;
				pMat->U[ke+k][j][i].Fr3 = 0.0;
			


    			}
		}
	}

  	return;
}

void radMHD_Mat_outflowks(MatrixS *pMat)
{
  	int i, j, k, je, ju, jl, js, ks, ke;
	int is,ie, il, iu;
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;

	ju = pMat->je + Matghost;
  	jl = pMat->js - Matghost;

	il = pMat->is - Matghost;
	iu = pMat->ie + Matghost;

	Real x1, x2, x3;
	Real vx, vy, vz, Fr0z, Sigma_t;
	Real reducefactor, Eratio, tau;

	for (k=1;  k<=Matghost;  k++) {
		for(j=js-Matghost; j<=je+Matghost; j++){
			for(i=is-Matghost; i<=ie+Matghost; i++){

			/*	pMat->Ugas[ks-k][j][i].Edd_11 = pMat->Ugas[ks][j][i].Edd_11;
				pMat->Ugas[ks-k][j][i].Edd_22 = pMat->Ugas[ks][j][i].Edd_22;
				pMat->Ugas[ks-k][j][i].Edd_21 = pMat->Ugas[ks][j][i].Edd_21;
				pMat->Ugas[ks-k][j][i].Edd_31 = pMat->Ugas[ks][j][i].Edd_31;
				pMat->Ugas[ks-k][j][i].Edd_32 = pMat->Ugas[ks][j][i].Edd_32;
				pMat->Ugas[ks-k][j][i].Edd_33 = pMat->Ugas[ks][j][i].Edd_33;
			*/
				pMat->Ugas[ks-k][j][i] = pMat->Ugas[ks][j][i];
		
				pMat->U[ks-k][j][i].Er  = 0.0;
				pMat->U[ks-k][j][i].Fr1 = 0.0;
				pMat->U[ks-k][j][i].Fr2 = 0.0;
				pMat->U[ks-k][j][i].Fr3 = 0.0;
				
    			}
		}
	}

  	return;
}




void radMHD_Mat_outflowie(MatrixS *pMat)
{
  	int i, j, k, je, ju, jl, js, ks, ke;
	int is,ie, il, iu;
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;

	ju = pMat->je + Matghost;
  	jl = pMat->js - Matghost;

	il = pMat->is - Matghost;
	iu = pMat->ie + Matghost;

	Real x1, x2, x3;
	Real vx, vy, vz, Fr0z, Sigma_t;
	Real reducefactor, Eratio, tau;

	for(k=ks-Matghost; k<=ke+Matghost; k++){	
		for(j=js-Matghost; j<=je+Matghost; j++){
			for(i=1; i<=Matghost; i++){
			
	
			/*	pMat->Ugas[ke+k][j][i].Edd_11 = pMat->Ugas[ke][j][i].Edd_11;
				pMat->Ugas[ke+k][j][i].Edd_22 = pMat->Ugas[ke][j][i].Edd_22;
				pMat->Ugas[ke+k][j][i].Edd_21 = pMat->Ugas[ke][j][i].Edd_21;
				pMat->Ugas[ke+k][j][i].Edd_31 = pMat->Ugas[ke][j][i].Edd_31;
				pMat->Ugas[ke+k][j][i].Edd_32 = pMat->Ugas[ke][j][i].Edd_32;
				pMat->Ugas[ke+k][j][i].Edd_33 = pMat->Ugas[ke][j][i].Edd_33;
			*/
				pMat->Ugas[k][j][ie+i] = pMat->Ugas[ke][j][i];
				
				pMat->U[k][j][ie+i].Er  = 0.0;
				pMat->U[k][j][ie+i].Fr1 = 0.0;
				pMat->U[k][j][ie+i].Fr2 = 0.0;
				pMat->U[k][j][ie+i].Fr3 = 0.0;
			


    			}
		}
	}

  	return;
}




void radMHD_Mat_periodis(MatrixS *pMat)
{
	int i, j, k, is, ie, js, je, ks, ke, m;
#ifdef MPI_PARALLEL
	int NRad = 15 + NOPACITY;	/* Need to transfer 11+NOPACITY variables for Ugas, plus 4 variables for U */
#endif
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;


#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif


	/* Do not need MPI communication */
	if(Periodix1_id < 0){
		for (k=ks; k<=ke; k++) {
    			for (j=js; j<=je; j++) {
      				for (i=1; i<=Matghost; i++) {

					pMat->Ugas[k][j][is-i] = pMat->Ugas[k][js+i-1][j-js+is];

					pMat->Ugas[k][j][is-i].V1 = -pMat->Ugas[k][js+i-1][j-js+is].V2;
					pMat->Ugas[k][j][is-i].V2 = pMat->Ugas[k][js+i-1][j-js+is].V1;

					pMat->Ugas[k][j][is-i].Edd_11 = pMat->Ugas[k][js+i-1][j-js+is].Edd_22;
					pMat->Ugas[k][j][is-i].Edd_22 = pMat->Ugas[k][js+i-1][j-js+is].Edd_11;
					pMat->Ugas[k][j][is-i].Edd_31 = -pMat->Ugas[k][js+i-1][j-js+is].Edd_32;
					pMat->Ugas[k][j][is-i].Edd_32 = pMat->Ugas[k][js+i-1][j-js+is].Edd_31;

					pMat->U[k][j][is-i] = pMat->U[k][js+i-1][j-js+is];
					pMat->U[k][j][is-i].Fr1 = -pMat->U[k][js+i-1][j-js+is].Fr2;
					pMat->U[k][j][is-i].Fr2 = pMat->U[k][js+i-1][j-js+is].Fr1;

      				}
    			}
  		}


	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
		/* First, count the number of data expected to receive */
		cntR = Matghost*(pMat->Nx[1])*(pMat->Nx[2])*(NRad);

		/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");

		/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodix1_id,1,
			pMat->Comm_Domain, &(recv_rq));
		
		/* Now count the number of data needs to send */
		cntS = (pMat->Nx[1] + 2*Matghost)*Matghost*(pMat->Nx[2])*(NRad);

		/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
      			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");

		/*--------------------------------------------------------------*/
		/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
		/* send data in the order that will be needed */
		
  		for (k=ks; k<=ke; k++) {
			for (i=is+(Matghost-1); i>=is; i--) {
    				for (j=js-Matghost; j<=je+Matghost; j++) {
					/* First Ugas, then U */
      				
        				*(pSnd++) = pMat->Ugas[k][j][i].rho;
					*(pSnd++) = pMat->Ugas[k][j][i].V2;
					*(pSnd++) = -pMat->Ugas[k][j][i].V1;
					*(pSnd++) = pMat->Ugas[k][j][i].V3;
					*(pSnd++) = pMat->Ugas[k][j][i].T4;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_22;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_21;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_11;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_32;
					*(pSnd++) = -pMat->Ugas[k][j][i].Edd_31;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_33;
					for(m=0;m<NOPACITY;m++){
						*(pSnd++) = pMat->Ugas[k][j][i].Sigma[m];
					}
					*(pSnd++) = pMat->U[k][j][i].Er;
        				*(pSnd++) = pMat->U[k][j][i].Fr2;
        				*(pSnd++) = -pMat->U[k][j][i].Fr1;
        				*(pSnd++) = pMat->U[k][j][i].Fr3;
      				}
    			}
  		}

		/* now actually send the  data */
		 ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodix1_id,1,
        		pMat->Comm_Domain, &(send_rq));

		 /* check non-blocking sends have completed. */
	      	ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);

		/*now get the receive data  */

		 /* check non-blocking receive have finished. */
     		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);

		pRcv = (double*)&(recv_buf[0]);

		/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
    			for (j=js; j<=je; j++){
      				for (i=is-Matghost; i<=is-1; i++){
					/* First Ugas then U */
					pMat->Ugas[k][j][i].rho 	= *(pRcv++);
        				pMat->Ugas[k][j][i].V1  	= *(pRcv++);
        				pMat->Ugas[k][j][i].V2 		= *(pRcv++);
        				pMat->Ugas[k][j][i].V3 		= *(pRcv++);
        				pMat->Ugas[k][j][i].T4 		= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_11  	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_21 	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_22 	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_31 	= *(pRcv++);
					pMat->Ugas[k][j][i].Edd_32 	= *(pRcv++);
					pMat->Ugas[k][j][i].Edd_33 	= *(pRcv++);
					for(m=0;m<NOPACITY;m++){	
						pMat->Ugas[k][j][i].Sigma[m]	= *(pRcv++);
					}

        				pMat->U[k][j][i].Er = *(pRcv++);
					pMat->U[k][j][i].Fr1 = *(pRcv++);
					pMat->U[k][j][i].Fr2 = *(pRcv++);
					pMat->U[k][j][i].Fr3 = *(pRcv++);

      				}
    			}
  		}



		/* Free the memory */
		free(recv_buf);
		free(send_buf);

#else /* if not MPI_PARALLEL */
		ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodix1_id);
#endif
	}/* End MPI case */



}




void radMHD_Mat_outflowje(MatrixS *pMat)
{
  	int i, j, k, je, ju, jl, js, ks, ke;
	int is,ie, il, iu;
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;

	ju = pMat->je + Matghost;
  	jl = pMat->js - Matghost;

	il = pMat->is - Matghost;
	iu = pMat->ie + Matghost;

	Real x1, x2, x3;
	Real vx, vy, vz, Fr0z, Sigma_t;
	Real reducefactor, Eratio, tau;

	for(k=ks-Matghost; k<=ke+Matghost; k++){
		for(i=is-Matghost; i<=ie+Matghost; i++){
			for(j=1; j<=Matghost; j++){
	
			/*	pMat->Ugas[ke+k][j][i].Edd_11 = pMat->Ugas[ke][j][i].Edd_11;
				pMat->Ugas[ke+k][j][i].Edd_22 = pMat->Ugas[ke][j][i].Edd_22;
				pMat->Ugas[ke+k][j][i].Edd_21 = pMat->Ugas[ke][j][i].Edd_21;
				pMat->Ugas[ke+k][j][i].Edd_31 = pMat->Ugas[ke][j][i].Edd_31;
				pMat->Ugas[ke+k][j][i].Edd_32 = pMat->Ugas[ke][j][i].Edd_32;
				pMat->Ugas[ke+k][j][i].Edd_33 = pMat->Ugas[ke][j][i].Edd_33;
			*/
				pMat->Ugas[k][je+j][i] = pMat->Ugas[k][je][i];
				
				pMat->U[k][je+j][i].Er  = 0.0;
				pMat->U[k][je+j][i].Fr1 = 0.0;
				pMat->U[k][je+j][i].Fr2 = 0.0;
				pMat->U[k][je+j][i].Fr3 = 0.0;
			


    			}
		}
	}

  	return;
}




void radMHD_Mat_periodjs(MatrixS *pMat)
{
	int i, j, k, is, ie, js, je, ks, ke, m;
	int NRad = 15 + NOPACITY;
	je = pMat->je;
	js = pMat->js;
	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;
												  
#ifdef MPI_PARALLEL
	int cntR, cntS, cnt2, cnt3, ierr, mIndex;
	double *recv_buf;
	double *send_buf;
	double *pSnd, *pRcv;
	MPI_Request recv_rq;
	MPI_Request send_rq;
#endif
												  
												  
	/* Do not need MPI communication */
	if(Periodjx1_id < 0){												  
		for (k=ks; k<=ke; k++) {
			for (j=1; j<=Matghost; j++) {
				for (i=is; i<=ie; i++) {
					
					pMat->Ugas[k][js-j][i] = pMat->Ugas[k][i-is+js][is+j-1];

					pMat->Ugas[k][js-j][i].V1 = pMat->Ugas[k][i-is+js][j+is-1].V2;
					pMat->Ugas[k][js-j][i].V2 = -pMat->Ugas[k][i-is+js][j+is-1].V1;

					pMat->Ugas[k][js-j][i].Edd_11 = pMat->Ugas[k][i-is+js][j+is-1].Edd_22;
					pMat->Ugas[k][js-j][i].Edd_22 = pMat->Ugas[k][i-is+js][j+is-1].Edd_11;
					pMat->Ugas[k][js-j][i].Edd_31 = pMat->Ugas[k][i-is+js][j+is-1].Edd_32;
					pMat->Ugas[k][js-j][i].Edd_32 = -pMat->Ugas[k][i-is+js][j+is-1].Edd_31;

					pMat->U[k][js-j][i] = pMat->U[k][i-is+js][is+j-1];

					pMat->U[k][js-j][i].Fr1 = pMat->U[k][i-is+js][j+is-1].Fr2;
					pMat->U[k][js-j][i].Fr2 = -pMat->U[k][i-is+js][j+is-1].Fr1;


				
					

					
				}
			}
		}
			
		 /* Now the corner */
  /* The corner will copy some data in the just updated ghost zones */
 		for (k=ks; k<=ke; k++) {
    			for (j=1; j<=Matghost; j++) {
      				for (i=is-Matghost; i<is; i++) {
					pMat->Ugas[k][js-j][i] = pMat->Ugas[k][i-is+js][is+j-1];

					pMat->Ugas[k][js-j][i].V1 = pMat->Ugas[k][i-is+js][j+is-1].V2;
					pMat->Ugas[k][js-j][i].V2 = -pMat->Ugas[k][i-is+js][j+is-1].V1;

					pMat->Ugas[k][js-j][i].Edd_11 = pMat->Ugas[k][i-is+js][j+is-1].Edd_22;
					pMat->Ugas[k][js-j][i].Edd_22 = pMat->Ugas[k][i-is+js][j+is-1].Edd_11;
					pMat->Ugas[k][js-j][i].Edd_31 = pMat->Ugas[k][i-is+js][j+is-1].Edd_32;
					pMat->Ugas[k][js-j][i].Edd_32 = -pMat->Ugas[k][i-is+js][j+is-1].Edd_31;

					pMat->U[k][js-j][i] = pMat->U[k][i-is+js][is+j-1];

					pMat->U[k][js-j][i].Fr1 = pMat->U[k][i-is+js][j+is-1].Fr2;
					pMat->U[k][js-j][i].Fr2 = -pMat->U[k][i-is+js][j+is-1].Fr1;

      				}
    			}
  		}

 		for (k=ks; k<=ke; k++) {	
    			for (j=1; j<=Matghost; j++) {
      				for (i=ie+1; i<=ie+Matghost; i++) {
        				pMat->Ugas[k][js-j][i] = pMat->Ugas[k][i-is+js][is+j-1];

					pMat->Ugas[k][js-j][i].V1 = pMat->Ugas[k][i-is+js][j+is-1].V2;
					pMat->Ugas[k][js-j][i].V2 = -pMat->Ugas[k][i-is+js][j+is-1].V1;

					pMat->Ugas[k][js-j][i].Edd_11 = pMat->Ugas[k][i-is+js][j+is-1].Edd_22;
					pMat->Ugas[k][js-j][i].Edd_22 = pMat->Ugas[k][i-is+js][j+is-1].Edd_11;
					pMat->Ugas[k][js-j][i].Edd_31 = pMat->Ugas[k][i-is+js][j+is-1].Edd_32;
					pMat->Ugas[k][js-j][i].Edd_32 = -pMat->Ugas[k][i-is+js][j+is-1].Edd_31;

					pMat->U[k][js-j][i] = pMat->U[k][i-is+js][is+j-1];

					pMat->U[k][js-j][i].Fr1 = pMat->U[k][i-is+js][j+is-1].Fr2;
					pMat->U[k][js-j][i].Fr2 = -pMat->U[k][i-is+js][j+is-1].Fr1;

      				}
    			}
  		}

									  
												  
	}/* End if no MPI call is needed */
	else{
#ifdef MPI_PARALLEL
	/* Now count the number of data needs to receive */
		cntR = (pMat->Nx[1] + 2*Matghost)*Matghost*(pMat->Nx[2])*(NRad);

												  
/* allocate memory for recv_buf and recv_rq */
		if((recv_buf = (double*)calloc_1d_array(cntR,sizeof(double))) == NULL)
			ath_error("[radMHD_periodix1]: Failed to allocate recv buffer\n");
												  
	/* post non-blocking receives for data from L grid */
		ierr = MPI_Irecv(&(recv_buf[0]),cntR,MPI_DOUBLE,Periodjx1_id,1,
				 pMat->Comm_Domain, &(recv_rq));
												  
												  
												  
	/* Count the number of data needs to send */
		cntS = Matghost*(pMat->Nx[1])*(pMat->Nx[2])*(NRad);

												  
	/* allocate memory for send buff */
		if((send_buf = (double*)calloc_1d_array(cntS,sizeof(double))) == NULL)
			ath_error("[radMHD_periodix1]: Failed to allocate send buffer\n");
											  
												  
 /*--------------------------------------------------------------*/
/* Now prepare the data to send */
		pSnd = (double*)&(send_buf[0]);
/* send data in the order that will be needed */
												  
		for (k=ks; k<=ke; k++) {			
			for(i=is; i<=ie; i++){    				
				for(j=js+(Matghost-1); j>=js; j--){
												  
					*(pSnd++) = pMat->Ugas[k][j][i].rho;
					*(pSnd++) = -pMat->Ugas[k][j][i].V2;
					*(pSnd++) = pMat->Ugas[k][j][i].V1;
					*(pSnd++) = pMat->Ugas[k][j][i].V3;
					*(pSnd++) = pMat->Ugas[k][j][i].T4;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_22;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_21;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_11;
					*(pSnd++) = -pMat->Ugas[k][j][i].Edd_32;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_31;
					*(pSnd++) = pMat->Ugas[k][j][i].Edd_33;
					for(m=0;m<NOPACITY;m++){
						*(pSnd++) = pMat->Ugas[k][j][i].Sigma[m];
					}
					*(pSnd++) = pMat->U[k][j][i].Er;
        				*(pSnd++) = -pMat->U[k][j][i].Fr2;
        				*(pSnd++) = pMat->U[k][j][i].Fr1;
        				*(pSnd++) = pMat->U[k][j][i].Fr3;					  


				}
			}
		}

																			
																			
																			
/* now actually send the  data */
		ierr = MPI_Isend(&(send_buf[0]),cntS,MPI_DOUBLE,Periodjx1_id,1,
				pMat->Comm_Domain, &(send_rq));
																			
/* check non-blocking sends have completed. */
		ierr = MPI_Waitall(1, &(send_rq), MPI_STATUS_IGNORE);
																			
/*now get the receive data  */
																			
/* check non-blocking receive have finished. */
		ierr = MPI_Waitany(1,&(recv_rq),&mIndex,MPI_STATUS_IGNORE);
																			
		pRcv = (double*)&(recv_buf[0]);
																			
/* Now have the receive data , unpack */
		for (k=ks; k<=ke; k++){
			for (j=js-Matghost; j<js; j++){
				for (i=is-Matghost; i<=ie+Matghost; i++){
					/* First Ugas then U */
					pMat->Ugas[k][j][i].rho 	= *(pRcv++);
        				pMat->Ugas[k][j][i].V1  	= *(pRcv++);
        				pMat->Ugas[k][j][i].V2 		= *(pRcv++);
        				pMat->Ugas[k][j][i].V3 		= *(pRcv++);
        				pMat->Ugas[k][j][i].T4 		= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_11  	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_21 	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_22 	= *(pRcv++);
        				pMat->Ugas[k][j][i].Edd_31 	= *(pRcv++);
					pMat->Ugas[k][j][i].Edd_32 	= *(pRcv++);
					pMat->Ugas[k][j][i].Edd_33 	= *(pRcv++);
					for(m=0;m<NOPACITY;m++){	
						pMat->Ugas[k][j][i].Sigma[m]	= *(pRcv++);
					}

        				pMat->U[k][j][i].Er = *(pRcv++);
					pMat->U[k][j][i].Fr1 = *(pRcv++);
					pMat->U[k][j][i].Fr2 = *(pRcv++);
					pMat->U[k][j][i].Fr3 = *(pRcv++);					  

				}
			}
		}
																			
																			
																			
	/* Free the memory */
	free(recv_buf);
	free(send_buf);
																			
#else /* if not MPI_PARALLEL */
	ath_error("[radMHD_periodis]: left id is : %d but no MPI is defined\n",Periodjx1_id);
#endif
	}/* End MPI case */
																	
}



#endif /* end multi_grid */

#endif /* end radMHD or radhydro */





/*! \fn static void output_1d(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D azimuthally averaged
    and vertically averaged radial profiles
     Extend to the SMR case. calculate Vertical average for each level */
/* Transform cartesian coordinate to spherical coordinate */
/* vr = cosphi * sintheta * vx + sintheta * sinphi * vy + costheta Vz *
 * vtheta = costheta * cosphi * vx + costheta * sinphi * vy - sintheta * vz
 * vphi = - sinphi * vx + cosphi * vy 	*/

/* Transform to Cynlindrical coordinate */
/**************************************
 * [vr, vphi] = {cosphi | sinphi, -sinphi |  cosphi} [vx, vy]
 *
 ***************************************/

static void output_1d(MeshS *pM, OutputS *pOut)
{
 int nl, nd;
 
 FILE *p_1dfile;
 char *fname, *plev=NULL, *pdom=NULL;
 char levstr[8], domstr[8];
 int big_end = ath_big_endian();
 
 int i,j,k;
 int tot1d,i1d,nrmx;
 int dnum = pOut->num;

  double **out1d;
  double x1,x2,x3,press;
  double cosphi, sinphi, distance, Radius;
  double *out_r;
  int *count_r;
  double vx, vy, vz, Fr01, Fr02, Fr03, Fr0r, Fr0phi;
  double vr, vphi, Fr, Fphi;
#if defined(MHD) || defined(RADIATION_MHD)
  double Br, Bphi;
#endif


#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int *my_count;
  int *g_count;
  int zproc;
  int ierr,myID_Comm_Domain;
  int nproc;
#endif

 /* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap, fxx, fyy, fzz, fxy, fxz, fyz */

#if defined(MHD) || defined(RADIATION_MHD)
  tot1d=51;
#elif defined(RADIATION_HYDRO)
  tot1d=41;
#else
  tot1d=19; 
#endif /* MHD */



  
  int is, ie, js, je, ks, ke, rpos;
/* Loop over all Domain and level */

  DomainS *pD;
  GridS *pGrid;

 for(nl=0; nl<(pM->NLevels); nl++){
	for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
		if(pM->Domain[nl][nd].Grid != NULL){

		if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          		(pOut->ndomain == -1 || pOut->ndomain == nd)){

			pD = &(pM->Domain[nl][nd]);
 			pGrid = pD->Grid;


			/* total number of radial grid zones for this Domain only */			
			nrmx = pD->Nx[0]; 


			/* index of this grid */
			is = pGrid->is; 
			ie = pGrid->ie;
			js = pGrid->js;
			je = pGrid->je;
			ks = pGrid->ks;
			ke = pGrid->ke;

#ifdef MPI_PARALLEL
  			nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  			ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  			if(ierr != MPI_SUCCESS)
    				ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  			/* We need to do it every time as it can be different for each domain at each level */
			/* count_r count how many cells are binned to a certain radial bins */
      			out_r = (double *) calloc_1d_array(nrmx,sizeof(double));
			count_r = (int *) calloc_1d_array(nrmx,sizeof(int));

  			

  			out1d = (double **) calloc_2d_array(nrmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  			my_out1d = (double *) calloc_1d_array(nrmx,sizeof(double));
  			g_out1d = (double *) calloc_1d_array(nrmx,sizeof(double));

			my_count = (int *) calloc_1d_array(nrmx,sizeof(int));
  			g_count = (int *) calloc_1d_array(nrmx,sizeof(int));
#endif
  			for (i=0; i<nrmx; i++) {
				count_r[i] = 0;
    				for (i1d=0; i1d<tot1d; i1d++) {
      					out1d[i][i1d] = 0.0;
    				}
  			}
  			

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */

			/* There will be some grids that not covered in the radial grids */
			/* This needs to be done for every CPU because we need this array to determine the position */				
			for (i=0; i<nrmx; i++) {
			/* Vertical coordinate for this whole domain */
				x1 = pD->MinX[0] + (i + 0.5)*pGrid->dx1;
      				out_r[i] = x1;
    			}
	
			/* Compute 1d averaged variables */
			/* The position of this grid in the radial array is determined by comparing the distance */
			 
  			for (k=ks; k<=ke; k++) {    				
    				for (j=js; j<=je; j++) {
      					for (i=is; i<=ie; i++) {
						/* First, calculate the transformation parameters */
						cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
						distance = sqrt(x1*x1+x2*x2+x3*x3);
						Radius = sqrt(x1*x1+x2*x2);						
						cosphi = x1/Radius;
						sinphi = x2/Radius;

						/* Determine the position of this cell in the array of radial distance */
						Radial_pos(&rpos, out_r, nrmx, Radius);
					if(rpos >= 0){
						/* counter increases */
						count_r[rpos]++;

						/* Now this cell should fill [rpos] in out_r and out1d */

						vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
						vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
						vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;

						vr   =  cosphi * vx + sinphi * vy;
						vphi = -sinphi * vx + cosphi * vy;
#if defined(MHD) || defined(RADIATION_MHD)
						Br   =  cosphi * pGrid->U[k][j][i].B1c + sinphi * pGrid->U[k][j][i].B2c;
						Bphi = -sinphi * pGrid->U[k][j][i].B1c + cosphi * pGrid->U[k][j][i].B2c;
#endif
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
						Fr = cosphi * pGrid->U[k][j][i].Fr1 + sinphi * pGrid->U[k][j][i].Fr2;
						Fphi = -sinphi * pGrid->U[k][j][i].Fr1 + cosphi * pGrid->U[k][j][i].Fr2;
#endif
						/* density */
        					i1d=0;
        					out1d[rpos][i1d] += pGrid->U[k][j][i].d;
        					i1d++;
#ifdef ISOTHERMAL
        					out1d[rpos][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        					press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
                                 				- expr_ME(pGrid,i,j,k)
#endif
                                				),TINY_NUMBER);
						/* pressure */
        					out1d[rpos][i1d] += press;
#endif
#ifdef ADIABATIC
        					i1d++;
						/* temperature */
        					out1d[rpos][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
        					i1d++;
						/* Total E */
        					out1d[rpos][i1d] += pGrid->U[k][j][i].E;
        					i1d++;
        					out1d[rpos][i1d] += hst_E_total(pGrid,i,j,k);
#endif
						/* momentum */
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].d * vr;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].d * vphi;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].d * vz;
						/* kinetic energy */
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(vr)*pGrid->U[k][j][i].d;
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(vphi)*pGrid->U[k][j][i].d;
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(vz)*pGrid->U[k][j][i].d;
        					i1d++;
        					out1d[rpos][i1d] += expr_KE(pGrid,i,j,k);
        					i1d++;
        					out1d[rpos][i1d] += hst_rho_VrVphi(pGrid,i,j,k);
						/* vr */
						i1d++;
						out1d[rpos][i1d] += vr;
						/* vphi */
						i1d++;
						out1d[rpos][i1d] += vphi;
						/* vz */
						i1d++;
						out1d[rpos][i1d] += vz;
#if defined(MHD) || defined(RADIATION_MHD)
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
        					i1d++;
        					out1d[rpos][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
        					i1d++;
        					out1d[rpos][i1d] += expr_ME(pGrid,i,j,k);
						/* Alfven velocity */
						i1d++;
						out1d[rpos][i1d] += sqrt(expr_ME(pGrid,i,j,k)/pGrid->U[k][j][i].d);
						/* Alfven velocity for Bz */
						i1d++;
						out1d[rpos][i1d] += sqrt(0.5*SQR(pGrid->U[k][j][i].B3c)/pGrid->U[k][j][i].d);
        					i1d++;
        					out1d[rpos][i1d] += hst_Bx(pGrid,i,j,k);
        					i1d++;
        					out1d[rpos][i1d] += hst_By(pGrid,i,j,k);
        					i1d++;
        					out1d[rpos][i1d] += hst_Bz(pGrid,i,j,k);
        					i1d++;
        					out1d[rpos][i1d] += hst_BrBphi(pGrid,i,j,k);
#endif
						/* rho v mass flux */
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].d * vr * out_r[rpos];
						
						/* angular momentum */
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].d * vphi * Radius;
						/* specific angular momentum */
						i1d++;
						out1d[rpos][i1d] += vphi * Radius;

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
						i1d++;
        					out1d[rpos][i1d] += pGrid->U[k][j][i].Er;
						i1d++;
        					out1d[rpos][i1d] += Fr;
						i1d++;
        					out1d[rpos][i1d] += Fphi;
						i1d++;
        					out1d[rpos][i1d] += pGrid->U[k][j][i].Fr3;
						/* co-moving flux */
						Fr01 = pGrid->U[k][j][i].Fr1 -(vx * (1.0 + pGrid->U[k][j][i].Edd_11) + vy * pGrid->U[k][j][i].Edd_21 + vz * pGrid->U[k][j][i].Edd_31) * pGrid->U[k][j][i].Er/Crat;
						Fr02 = pGrid->U[k][j][i].Fr2 -(vy * (1.0 + pGrid->U[k][j][i].Edd_22) + vx * pGrid->U[k][j][i].Edd_21 + vz * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;
						Fr03 = pGrid->U[k][j][i].Fr3 -(vz * (1.0 + pGrid->U[k][j][i].Edd_33) + vx * pGrid->U[k][j][i].Edd_31 + vy * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;
						i1d++;
						out1d[rpos][i1d] += Fr01;
						i1d++;
						out1d[rpos][i1d] += Fr02;
						i1d++;
						out1d[rpos][i1d] += Fr03;
						i1d++;
						/* Fr0r */
						out1d[rpos][i1d] += (cosphi * Fr01 + sinphi * Fr02);
						/* Fr0phi */
						i1d++;
						out1d[rpos][i1d] += (-sinphi * Fr01 + cosphi * Fr02);

						/* advection Er */
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Er * vr;
						/* kappes and kappaff */

						i1d++;
						out1d[rpos][i1d] += hst_sigmas(pGrid,i,j,k);
						i1d++;
        					out1d[rpos][i1d] += hst_sigmaaP(pGrid,i,j,k);
						/* radiation work term */
						i1d++;
						out1d[rpos][i1d] += (pGrid->U[k][j][i].Sigma[0] - pGrid->U[k][j][i].Sigma[1]) * (Fr01 * vx + Fr02 * vy + Fr03 * vz); 
						/* dErdx */
						i1d++;
						out1d[rpos][i1d] += (pGrid->U[k][j][i+1].Er - pGrid->U[k][j][i-1].Er)/(2.0 * pGrid->dx1);
						i1d++;
						out1d[rpos][i1d] += (pGrid->U[k][j+1][i].Er - pGrid->U[k][j-1][i].Er)/(2.0 * pGrid->dx2);
						i1d++;
						out1d[rpos][i1d] += (pGrid->U[k+1][j][1].Er - pGrid->U[k-1][j][i].Er)/(2.0 * pGrid->dx3);
						/* Eddington tensor */
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_11;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_22;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_33;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_21;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_31;
						i1d++;
						out1d[rpos][i1d] += pGrid->U[k][j][i].Edd_32;
	
#endif
					} /* if rpos > 0 */
      					}/* end i */
    				}/* end j */
  			}/* end k*/

 

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  			for(i1d=0; i1d<tot1d; i1d++){
    				for (k=0; k<nrmx; k++) {
      					my_out1d[k] = out1d[k][i1d];
    			}
    				ierr = MPI_Reduce(my_out1d, g_out1d, nrmx,
                      			MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    				if(ierr)
      					ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    				for (k=0; k<nrmx; k++) {
      					out1d[k][i1d] = g_out1d[k];
    				}
  			}

			/* sum the counter */
			for (k=0; k<nrmx; k++) {
      				my_count[k] = count_r[k];
    			}
    			ierr = MPI_Reduce(my_count, g_count, nrmx,
                      			MPI_INTEGER, MPI_SUM, 0, pD->Comm_Domain);
    			if(ierr)
      				ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    			for (k=0; k<nrmx; k++) {
      				count_r[k] = g_count[k];
    			}
#endif


/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  			if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  				
  				for (k=0; k<nrmx; k++) {
    					for (i1d=0; i1d<tot1d; i1d++) {
						if(count_r[k] > 0)
      						out1d[k][i1d] /= count_r[k];
    					}
  				}

/* Generate filename */
/* Adopt the way file name is generated from output_vtk_3d */
			if(nl > 0){
				plev = &levstr[0];
				sprintf(plev,"lev%d",nl);
			}
			if(nd > 0){
				pdom = &domstr[0];
				sprintf(pdom,"dom%d",nd);
			}




			if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
      					dnum,NULL,"1d")) == NULL){
    				ath_error("[output_1d]: Error constructing filename\n");
  			}



/* open filename */
  			p_1dfile = fopen(fname,"w");
  			if (p_1dfile == NULL) {
    				ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    				return;
  			}

/* Write out data */

  			for (k=0; k<nrmx; k++) {

#ifdef RADIATION_MHD
    				if (k == 0) {
      					fprintf(p_1dfile,"# [1]r	[2]dens	[3]pressure	[4]temperature	[5]E	[6]Etot	[7]Mr	[8]Mphi	[9]Mz	[10]KEr	[11]KEphi	[12]KEz	[13]KE	[14]rhoVrVphi	[15]Vr	[16]Vphi	[17]Vz	[18]MEx	[19]MEy	[20]MEz	[21]ME	[22]Va	[23]Vaz	[24]Bx	[25]By	[26]Bz	[27]BrBphi	[28]rhoVrR	[29]rhoVphiR	[30]VphiR	[31]Er	[32]Frr	[33]Fphi	[34]Fz	[35]Fr0x	[36]Fr0y	[37]Fr0z	[38]Fr0r	[39]Fr0phi	[40]VrEr	[41]kappaes	[42]kappaff	[43]Frwork	[44]dErdx	[45]dErdy	[46]dErdz	[47]f_xx	[48]f_yy	[49]f_zz	[50]f_xy	[51]f_xz	[52]f_yz\n");
    				}
    				fprintf(p_1dfile,"%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G\n",out_r[k],out1d[k][0],out1d[k][1],out1d[k][2],
            			out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            			out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28],out1d[k][29],out1d[k][30],out1d[k][31],out1d[k][32],out1d[k][33],out1d[k][34],out1d[k][35],out1d[k][36],out1d[k][37],out1d[k][38],out1d[k][39],out1d[k][40],out1d[k][41],out1d[k][42],out1d[k][43],out1d[k][44],out1d[k][45],out1d[k][46],out1d[k][47],out1d[k][48],out1d[k][49],out1d[k][50]);
#endif /* RADIATION_MHD */

  }

  				fclose(p_1dfile);
  				free(fname);
#ifdef MPI_PARALLEL
  			}/* End write 1D files */
#endif

			/* These temporary arrays are generated for each level and each domain */
  			free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  			free_1d_array(my_out1d); /* Free the memory we malloc'd */
  			free_1d_array(g_out1d); /* Free the memory we malloc'd */

			free_1d_array(my_count);
			free_1d_array(g_count);
#endif
 
			  			/* We need to do it every time as it can be different for each domain at each level */

				/* only need to free it for myID == 0 */
      				free_1d_array(out_r);
				free_1d_array(count_r);


		}/* End if level and domain match for this domain */
		}/* End if this CPU works on this domain */
	} /* End loop over domains at level nl */
 }/* end loop all levels */

	

return;
}





/*! \fn static void output_2d_binary(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate the azimuthal averaged 
 * 2D profiles: radial and vertical */
 
 
/* Transform to Cynlindrical coordinate */
/**************************************
 * [vr, vphi] = {cosphi | sinphi, -sinphi |  cosphi} [vx, vy]
 *
 ***************************************/

static void output_2d_binary(MeshS *pM, OutputS *pOut)
{
	int nl, nd;
	
	FILE *p_1dfile;
	char *fname, *plev=NULL, *pdom=NULL;
	char levstr[8], domstr[8];
	int big_end = ath_big_endian();
	
	int i,j,k;
	int tot1d,i1d,nrmx,nzmx,kg,kdispG,kdispD,index;
	int dnum = pOut->num;
	
	double **out1d;
	double x1,x2,x3,press;
	double cosphi, sinphi, distance, Radius;
	double *out_r, *out_x3;
	int *count_r;
	double vx, vy, vz, Fr01, Fr02, Fr03, Fr0r, Fr0phi;
	double vr, vphi, Fr, Fphi;
#if defined(MHD) || defined(RADIATION_MHD)
	double Br, Bphi;
#endif
	
	
#ifdef MPI_PARALLEL
	double *my_out1d;
	double *g_out1d;
	int *my_count;
	int *g_count;
	int zproc;
	int ierr,myID_Comm_Domain;
	int nproc;
#endif
	
	/* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap, fxx, fyy, fzz, fxy, fxz, fyz */
	
#if defined(MHD) || defined(RADIATION_MHD)
  	tot1d=51;
#elif defined(RADIATION_HYDRO)
  	tot1d=41;
#else
  	tot1d=19; 
#endif /* MHD */
	
	
	
	
	int is, ie, js, je, ks, ke, rpos;
	/* Loop over all Domain and level */
	
	DomainS *pD;
	GridS *pGrid;
	
	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				
				if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
					(pOut->ndomain == -1 || pOut->ndomain == nd)){
					
					pD = &(pM->Domain[nl][nd]);
					pGrid = pD->Grid;
					
					
					/* total number of radial grid zones for this Domain only */			
					nrmx = pD->Nx[0]; 
					nzmx = pD->Nx[2];
					
					
					/* index of this grid */
					is = pGrid->is; 
					ie = pGrid->ie;
					js = pGrid->js;
					je = pGrid->je;
					ks = pGrid->ks;
					ke = pGrid->ke;
					
#ifdef MPI_PARALLEL
					nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif
					
#ifdef MPI_PARALLEL
					ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
					if(ierr != MPI_SUCCESS)
						ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
					/* We need to do it every time as it can be different for each domain at each level */
					/* count_r count how many cells are binned to a certain radial bins */
					out_r = (double *) calloc_1d_array(nrmx,sizeof(double));
					
					count_r = (int *) calloc_1d_array(nrmx*nzmx,sizeof(int));
					
					/* We need to do it every time as it can be different for each domain at each level */
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
					}/* End myID */
#endif
					

									
					out1d = (double **) calloc_2d_array(nrmx*nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
					my_out1d = (double *) calloc_1d_array(nrmx*nzmx,sizeof(double));
					g_out1d = (double *) calloc_1d_array(nrmx*nzmx,sizeof(double));
					
					my_count = (int *) calloc_1d_array(nrmx*nzmx,sizeof(int));
					g_count = (int *) calloc_1d_array(nrmx*nzmx,sizeof(int));
#endif
					
					kdispD = pD->Disp[2];
					kdispG = pGrid->Disp[2];
					
					for (i=0; i<nrmx*nzmx; i++) {
						count_r[i] = 0;
						for (i1d=0; i1d<tot1d; i1d++) {
							out1d[i][i1d] = 0.0;
						}
					}
					
					
					/* First calculate the x3 coordinate and save it to be dumped
					 by root in every 1d file */
					
					/* There will be some grids that not covered in the radial grids */
					/* This needs to be done for every CPU because we need this array to determine the position */				
					for (i=0; i<nrmx; i++) {
						/* Vertical coordinate for this whole domain */
						x1 = pD->MinX[0] + (i + 0.5)*pGrid->dx1;
						out_r[i] = x1;
					}
					
					
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						for (k=0; k<nzmx; k++) {
							/* Vertical coordinate for this whole domain */
							x3 = pD->MinX[2] + (k + 0.5)*pGrid->dx3;
							out_x3[k] = x3;
						}
#ifdef MPI_PARALLEL
					}/* End my ID */
#endif
					
					
					
					/* Compute 1d averaged variables */
					/* The position of this grid in the radial array is determined by comparing the distance */
					
					for (k=ks; k<=ke; k++) {
						kg=k+kdispG-kdispD-nghost;
						for (j=js; j<=je; j++) {
							for (i=is; i<=ie; i++) {
								/* First, calculate the transformation parameters */
								cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
								distance = sqrt(x1*x1+x2*x2+x3*x3);
								Radius = sqrt(x1*x1+x2*x2);						
								cosphi = x1/Radius;
								sinphi = x2/Radius;
								
								/* Determine the position of this cell in the array of radial distance */
								Radial_pos(&rpos, out_r, nrmx, Radius);
								if(rpos >= 0){
									index = rpos * nzmx + kg;
									/* counter increases */
									count_r[index]++;
									
									/* Now this cell should fill [rpos] in out_r and out1d */
									
									vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
									vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
									vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;
									
									vr   =  cosphi * vx + sinphi * vy;
									vphi = -sinphi * vx + cosphi * vy;
#if defined(MHD) || defined(RADIATION_MHD)
									Br   =  cosphi * pGrid->U[k][j][i].B1c + sinphi * pGrid->U[k][j][i].B2c;
									Bphi = -sinphi * pGrid->U[k][j][i].B1c + cosphi * pGrid->U[k][j][i].B2c;
#endif
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
									Fr = cosphi * pGrid->U[k][j][i].Fr1 + sinphi * pGrid->U[k][j][i].Fr2;
									Fphi = -sinphi * pGrid->U[k][j][i].Fr1 + cosphi * pGrid->U[k][j][i].Fr2;
#endif
									i1d=0;
									out1d[index][i1d] += pGrid->U[k][j][i].d;
									i1d++;
#ifdef ISOTHERMAL
									out1d[index][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
									press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
																   - expr_ME(pGrid,i,j,k)
#endif
																   ),TINY_NUMBER);
									out1d[index][i1d] += press;
#endif
#ifdef ADIABATIC
									i1d++;
									out1d[index][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].E;
									i1d++;
									out1d[index][i1d] += hst_E_total(pGrid,i,j,k);
#endif
									/* momentum */
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].d * vr;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].d * vphi;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].d * vz;
									/* kinetic */
									i1d++;
									out1d[index][i1d] += 0.5*SQR(vr)*pGrid->U[k][j][i].d;
									i1d++;									
									out1d[index][i1d] += 0.5*SQR(vphi)*pGrid->U[k][j][i].d;		
									i1d++;
									out1d[index][i1d] += 0.5*SQR(vz)*pGrid->U[k][j][i].d;
									i1d++;
									out1d[index][i1d] += expr_KE(pGrid,i,j,k);
									i1d++;
									out1d[index][i1d] += hst_rho_VrVphi(pGrid,i,j,k);
									/* vr */
									i1d++;
									out1d[index][i1d] += vr;
									/* vphi */
									i1d++;
									out1d[index][i1d] += vphi;
									/* vz */
									i1d++;
									out1d[index][i1d] += vz;
#if defined(MHD) || defined(RADIATION_MHD)
									i1d++;
									out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
									i1d++;
									out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
									i1d++;
									out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
									i1d++;
									out1d[index][i1d] += expr_ME(pGrid,i,j,k);
									/* Alfven velocity */
									i1d++;
									out1d[index][i1d] += sqrt(expr_ME(pGrid,i,j,k)/pGrid->U[k][j][i].d);
									/* Alfven velocity for Bz */
									i1d++;
									out1d[index][i1d] += sqrt(0.5*SQR(pGrid->U[k][j][i].B3c)/pGrid->U[k][j][i].d);
									i1d++;
									out1d[index][i1d] += hst_Bx(pGrid,i,j,k);
									i1d++;
									out1d[index][i1d] += hst_By(pGrid,i,j,k);
									i1d++;
									out1d[index][i1d] += hst_Bz(pGrid,i,j,k);
									i1d++;
									out1d[index][i1d] += hst_BrBphi(pGrid,i,j,k);
#endif
									/* rho v mass flux */
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].d * vr * out_r[rpos];
									
									/* angular momentum */
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].d * vphi * Radius;
									/* specific angular momentum */
									i1d++;
									out1d[index][i1d] += vphi * Radius;
									
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Er;
									i1d++;
									out1d[index][i1d] += Fr;
									i1d++;
									out1d[index][i1d] += Fphi;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Fr3;
									/* co-moving flux */
									Fr01 = pGrid->U[k][j][i].Fr1 -(vx * (1.0 + pGrid->U[k][j][i].Edd_11) + vy * pGrid->U[k][j][i].Edd_21 + vz * pGrid->U[k][j][i].Edd_31) * pGrid->U[k][j][i].Er/Crat;
									Fr02 = pGrid->U[k][j][i].Fr2 -(vy * (1.0 + pGrid->U[k][j][i].Edd_22) + vx * pGrid->U[k][j][i].Edd_21 + vz * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;
									Fr03 = pGrid->U[k][j][i].Fr3 -(vz * (1.0 + pGrid->U[k][j][i].Edd_33) + vx * pGrid->U[k][j][i].Edd_31 + vy * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;
									i1d++;
									out1d[index][i1d] += Fr01;
									i1d++;
									out1d[index][i1d] += Fr02;
									i1d++;
									out1d[index][i1d] += Fr03;
									i1d++;
									/* Fr0r */
									out1d[index][i1d] += (cosphi * Fr01 + sinphi * Fr02);
									/* Fr0phi */
									i1d++;
									out1d[index][i1d] += (-sinphi * Fr01 + cosphi * Fr02);

									/* advection Er */
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Er * vr;
									/* kappes and kappaff */
									
									i1d++;
									out1d[index][i1d] += hst_sigmas(pGrid,i,j,k);
									i1d++;
									out1d[index][i1d] += hst_sigmaaP(pGrid,i,j,k);
									/* radiation work term */
									i1d++;
									out1d[index][i1d] += (pGrid->U[k][j][i].Sigma[0] - pGrid->U[k][j][i].Sigma[1]) * (Fr01 * vx + Fr02 * vy + Fr03 * vz); 
									/* dErdx */
									i1d++;
									out1d[index][i1d] += (pGrid->U[k][j][i+1].Er - pGrid->U[k][j][i-1].Er)/(2.0 * pGrid->dx1);
									i1d++;
									out1d[index][i1d] += (pGrid->U[k][j+1][i].Er - pGrid->U[k][j-1][i].Er)/(2.0 * pGrid->dx2);
									i1d++;
									out1d[index][i1d] += (pGrid->U[k+1][j][1].Er - pGrid->U[k-1][j][i].Er)/(2.0 * pGrid->dx3);
									/* Eddington tensor */
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_11;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_22;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_33;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_21;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_31;
									i1d++;
									out1d[index][i1d] += pGrid->U[k][j][i].Edd_32;
									
#endif
								} /* if rpos > 0 */
							}/* end i */
						}/* end j */
					}/* end k*/
					
					
					
					/* The parent sums the scal[] array.
					 * Note that this assumes (dx1,dx2,dx3) = const. */
					
#ifdef MPI_PARALLEL 
					for(i1d=0; i1d<tot1d; i1d++){
						for (k=0; k<nrmx*nzmx; k++) {
							my_out1d[k] = out1d[k][i1d];
						}
						ierr = MPI_Reduce(my_out1d, g_out1d, nrmx*nzmx,
										  MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
						if(ierr)
							ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
						for (k=0; k<nrmx*nzmx; k++) {
							out1d[k][i1d] = g_out1d[k];
						}
					}
					
					/* sum the counter */
					for (k=0; k<nrmx*nzmx; k++) {
						my_count[k] = count_r[k];
					}
					ierr = MPI_Reduce(my_count, g_count, nrmx*nzmx,
									  MPI_INTEGER, MPI_SUM, 0, pD->Comm_Domain);
					if(ierr)
						ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
					for (k=0; k<nrmx*nzmx; k++) {
						count_r[k] = g_count[k];
					}
#endif
					
					
					/* For parallel calculations, only the parent computes the average
					 * and writes the output. */
#ifdef MPI_PARALLEL
					if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif
						
						
						for (k=0; k<nrmx*nzmx; k++) {
							for (i1d=0; i1d<tot1d; i1d++) {
								if(count_r[k] > 0)
									out1d[k][i1d] /= count_r[k];
							}
						}
						
						/* Generate filename */
						/* Adopt the way file name is generated from output_vtk_3d */
						if(nl > 0){
							plev = &levstr[0];
							sprintf(plev,"lev%d",nl);
						}
						if(nd > 0){
							pdom = &domstr[0];
							sprintf(pdom,"dom%d",nd);
						}
						
						
						
						
						if((fname = ath_fname(plev,pM->outfilename,plev,pdom,num_digit,
											  dnum,NULL,"2db")) == NULL){
							ath_error("[output_2d_binary]: Error constructing filename\n");
						}
						
						
						
						/* open filename */
						p_1dfile = fopen(fname,"wb");
						if (p_1dfile == NULL) {
							ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
							return;
						}
						
						/* Write out data */
						/* First, wirte tot1d, nrmx, nzmx */
						fwrite(&tot1d,sizeof(int),1,p_1dfile);
						fwrite(&nrmx,sizeof(int),1,p_1dfile);
						fwrite(&nzmx,sizeof(int),1,p_1dfile);
						
						/* Then write vertical coordinate and radial coordinate */		
						fwrite(out_r,sizeof(Real),nrmx,p_1dfile);	
					
						fwrite(out_x3,sizeof(Real),nzmx,p_1dfile);	
						
						/* The name of each colum is */
						/*
						 fprintf(p_1dfile,"# [1]r	[2]dens	[3]pressure	[4]temperature	[5]E	[6]Etot	[7]KEr	[8]KEphi	[9]KEz	[10]KE	[11]rhoVrVphi	[12]MEx	[13]MEy	[14]MEz	[15]ME	[16]Va	[17]Vaz	[18]Bx	[19]By	[20]Bz	[21]BrBphi	[22]rhoVr	[23]rhoVphiR	[24]Er	[25]Fr	[26]Fphi	[27]Fz	[28]Fr0z	[29]VrEr	[30]kappaes	[31]kappaff	[32]f_xx	[33]f_yy	[34]f_zz	[35]f_xy	[36]f_xz	[37]f_yz\n");
						*/
						
						/* Now write the data */
						for (k=0; k<nrmx*nzmx; k++) {
								fwrite(out1d[k],sizeof(Real),tot1d,p_1dfile);
						}
							
						
						fclose(p_1dfile);
						free(fname);
#ifdef MPI_PARALLEL
					}/* End write 1D files */
#endif
					
					/* These temporary arrays are generated for each level and each domain */
					free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
					free_1d_array(my_out1d); /* Free the memory we malloc'd */
					free_1d_array(g_out1d); /* Free the memory we malloc'd */
					
					free_1d_array(my_count);
					free_1d_array(g_count);
#endif
					
					/* We need to do it every time as it can be different for each domain at each level */
					
					/* only need to free it for myID == 0 */
      				free_1d_array(out_r);
					free_1d_array(count_r);
					
					
				}/* End if level and domain match for this domain */
			}/* End if this CPU works on this domain */
		} /* End loop over domains at level nl */
	}/* end loop all levels */
	
	
	
	return;
}


void Radial_pos(int *rpos, const Real *out_r, const int nrmx, const Real Radius)
{
	int i;
	if(Radius > out_r[nrmx-1]){
		*rpos = -1;
		
	}
	else if(Radius < out_r[0]){
		*rpos = 0;
	}
	else{
		/* Now is within range */
		i = 0;
		while((Radius > 0.5 * (out_r[i] + out_r[i+1])) && (i < nrmx-1)){
			i++;
		}
		*rpos = i;
	}


	return;
}



