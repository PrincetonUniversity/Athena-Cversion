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
static Real grav_vertical(const Real x1, const Real x2, const Real x3);
static Real expr_dV2(const GridS *pG, const int i, const int j, const int k);
static Real expr_beta(const GridS *pG, const int i, const int j, const int k);
static Real expr_ME(const GridS *pG, const int i, const int j, const int k);
static Real expr_KE(const GridS *pG, const int i, const int j, const int k);
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real hst_rho_dVy2(const GridS *pG,const int i, const int j, const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
static Real hst_sigmas(const GridS *pG, const int i, const int j, const int k);
static Real hst_sigmaaP(const GridS *pG, const int i, const int j, const int k);
#endif

static Real hst_gravpot(const GridS *pG, const int i, const int j, const int k);


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

static Real hst_P(const GridS *pG, const int i, const int j, const int k);

static Real hst_rho2(const GridS *pG, const int i, const int j, const int k);

static Real hst_T2(const GridS *pG, const int i, const int j, const int k);

/* Integral of the stress over the radial surface */
static Real hst_stressL(const GridS *pG, const int i, const int j, const int k);
static Real hst_stressR(const GridS *pG, const int i, const int j, const int k);

/* mass and energy flux at top and bottom */
static Real hst_rhofluxtop(const GridS *pG, const int i, const int j, const int k);
static Real hst_rhofluxbottom(const GridS *pG, const int i, const int j, const int k);

static Real hst_Erfluxbottom(const GridS *pG, const int i, const int j, const int k);
static Real hst_Erfluxtop(const GridS *pG, const int i, const int j, const int k);

static void output_1d(MeshS *pM, OutputS *pOut);
static void output_1dx(MeshS *pM, OutputS *pOut);

/* Function for inflow boundary condition */

void radMHD_rad_inflowks(GridS *pGrid);
void radMHD_inflowks(GridS *pGrid);
void radMHD_rad_inflowke(GridS *pGrid);
void radMHD_inflowke(GridS *pGrid);


static Real kappaes;
static Real kappaffR;
static Real kappaffP;

static Real B0y;
static Real B0z;

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
static void Thindiskopacity(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[4]);
#endif
#ifndef SHEARING_BOX
static Real Omega_0 = 4.749736;
static Real qshear = 1.5;
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID
void radMHD_Mat_inflowke(MatrixS *pMat);
void radMHD_Mat_inflowks(MatrixS *pMat);

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

static Real ghostEr[4];
static Real ghostFr[4];
static Real ghostrho[4];

static Real inidata[3];
/* Save the initial Er, Fr and rho at ke */

static Real Tfloor = 0.03;
static Real dfloor = 1.e-5;

#define Dataline 256


/* For transfer module */
#ifdef RADIATION_TRANSFER


static Real Thermal_B(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
		    const int k);
static Real Transfereps(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
		      const int k);
static Real transfer_opacity(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
			  const int k);



#endif

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
  Real den = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real density, pressure, Er, Fr;
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

	if((pGrid->Nx[2] * NGz) != Dataline)
		ath_error("[Problem]: Input Data line: %d doesn't match vertical grid number: %d\n",Dataline,pGrid->Nx[2] * NGz);

	/* y direction cannot be periodic boundary condition */
	

	lines = kproc * Dataline/NGz;
	/* skip lines which are read by other CPUs */
	FILE *frho, *fEr, *fFr, *fp;
	if ( (fFr=fopen("./data/Fr_out.txt","r"))==NULL )
	{   
		printf("Open input fFr file error");
		return;
	
	}
	if ( (frho=fopen("./data/rho_out.txt","r"))==NULL )
	{   
		printf("Open input frho file error");
		return;
	
	}
	if ( (fEr=fopen("./data/Er_out.txt","r"))==NULL )
	{   
		printf("Open input fEr file error");
		return;
	
	}
	if ( (fp=fopen("./data/pressure_out.txt","r"))==NULL )
	{   
		printf("Open input fp file error");
		return;
	
	}


	/* First, read in ghost zones */
	for(i=0;i<4;i++){
		fscanf(frho,"%lf",&(ghostrho[3-i]));
	 	fscanf(fEr,"%lf",&(ghostEr[3-i]));
		fscanf(fFr,"%lf",&(ghostFr[3-i]));
		ghostFr[3-i] = factor * ghostFr[3-i];
	}



	for(i=0; i<lines;i++){
		fscanf(frho,"%lf",&(density));
	  	fscanf(fp,"%lf",&(pressure));
	 	fscanf(fEr,"%lf",&(Er));
		fscanf(fFr,"%lf",&(Fr));
	}



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
	betay=35.0;
	pres = 1.0;

  	B0z = 0.0;
  	B0y = sqrt((double)(2.0*pres/betay));
	B0 = sqrt(B0z * B0z + B0y * B0y);	

	kappaes = 2.719701e4;
	kappaffP = 182.926;
	kappaffR = 4.94395;

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs)) + ID;

/* Initialize boxsize */
  ztop = pDomain->RootMaxX[2];
  zbtm = pDomain->RootMinX[2];

  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
  kx = 16.0*PI/Lx;

	amp = 0.01;
	Bamp = 1.0;


/* Rescale amp to sound speed for ipert 2,3 */


  for (k=ks; k<=ke; k++) {
	fscanf(frho,"%lf",&(density));
	fscanf(fp,"%lf",&(pressure));
	fscanf(fEr,"%lf",&(Er));
	fscanf(fFr,"%lf",&(Fr));
	pressure = density * pow(Er, 0.25) * R_ideal;
	Fr = factor * Fr;
	

/*
	B0y = sqrt((double)(2.0*pressure/betay));
	B0 = sqrt(B0z * B0z + B0y * B0y);	
*/
	if(((k == ke) && (kproc == (NGz-1))) || ((k == ks) && kproc == 0)){
		inidata[0] = density;
		inidata[1] = Er;
		inidata[2] = -fabs(Fr);
	}

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
	if(fabs(x3) > 2.0)
                rval = 0.0;
#ifdef ADIABATIC
        rp = pressure*(1.0 + 2.0*rval);
        rd = density;
#else
        rd = density*(1.0 + 2.0*rval);
#endif
/* To conform to HGB, the perturbations to V/Cs are (1/5)amp/sqrt(Gamma)  */
        rval = amp*(ran2(&iseed) - 0.5);
	if(fabs(x3) > 2.0)
                rval = 0.0;

        rvx = 0.0 * 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
	if(fabs(x3) > 2.0)
                rval = 0.0;

        rvy = 0.0 * 0.4*rval*sqrt(pres/den);

        rval = amp*(ran2(&iseed) - 0.5);
	if(fabs(x3) > 2.0)
                rval = 0.0;

        rvz = 0.0 * 0.4*rval*sqrt(pres/den);
      

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
		pGrid->U[k][j][i].Fr1 = 0.0;
		pGrid->U[k][j][i].Fr2 = 0.0;
#ifdef SHEARING_BOX
		pGrid->U[k][j][i].Fr2 = -qshear * Omega_0 * x1 * (1.0 + pGrid->U[k][j][i].Edd_22) * pGrid->U[k][j][i].Er/Crat;	
#else
		pGrid->U[k][j][i].Fr2 = -1.5 *  x1 * (1.0 + pGrid->U[k][j][i].Edd_22) * pGrid->U[k][j][i].Er/Crat;
#endif
		pGrid->U[k][j][i].Fr3 = Fr;
		
		
#endif
		
		
		
		
		
    }
  }}

#if defined(MHD) || defined(RADIATION_MHD)

	if (ifield == 6) {
  /* flux tube of Hirose et al. We put this down here to break away from the
     for loops used above. */
    		Real xc=0.0,zc=0.0,Bpratio=0.125,Bp0,rad0,rad,Ay0,x1h,x3h;
    	
	 Real ***Ay;
  	if((Ay = (Real***)calloc_3d_array(pGrid->Nx[2]+2*nghost,pGrid->Nx[1]+2*nghost, pGrid->Nx[0]+2*nghost ,sizeof(Real))) == NULL)
		ath_error("[problem]: malloc return a NULL pointer\n");


    for (k=ks; k<=ke+2; k++) {
      for (j=js; j<=je+2; j++) {
        for (i=is; i<=ie+2; i++) {
          cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
/*
	 x1 = pGrid->MinX[0] + ((Real)(i - pGrid->is) + 0.5)*pGrid->dx1;
 	 x2 = pGrid->MinX[1] + ((Real)(j - pGrid->js) + 0.5)*pGrid->dx2;
	 x3 = pGrid->MinX[2] + ((Real)(k - pGrid->ks) + 0.5)*pGrid->dx3;
*/
          x1h = x1-0.5*pGrid->dx1;
          x3h = x3-0.5*pGrid->dx3;
          Bp0 = B0*Bpratio;
  /* We are assuming that the scale height is H = 1 as usual */
	 rad0 = Lx/4.0;
          Ay0 = Bp0*rad0/PI;
         if(x3h > 0.0)
            rad = sqrt(SQR(x1h-xc)+SQR(x3h-Lx/4.0));
         else
            rad = sqrt(SQR(x1h-xc)+SQR(x3h+Lx/4.0));


        /* Only limited to one scale height */
        /*   rad = fabs((x1h-xc)/(0.5) + fabs(x3h-zc)/(0.25));
        */
  /* Calculate vector potential */
          if (rad < rad0) {
              Ay[k][j][i]=-Ay0*(1.0+cos(PI*rad/rad0));
             if(x3h < 0.0)
                Ay[k][j][i]=Ay0*(1.0+cos(PI*rad/rad0));
          } else {
            Ay[k][j][i]=0.0;
          }
	
        }
      }
    }
  /* In this case, we calculate face fields from vect. potential */
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie+1; i++) {
          pGrid->B1i[k][j][i] = -(Ay[k+1][j][i]-Ay[k][j][i])/pGrid->dx3;
          pGrid->B3i[k][j][i] =  (Ay[k][j][i+1]-Ay[k][j][i])/pGrid->dx1;		
        }
      }
    }
  /* Sync poloidal centered fields to face fields for the next step */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i]+pGrid->B1i[k][j][i+1]);
          pGrid->U[k][j][i].B3c = 0.5*(pGrid->B3i[k][j][i]+pGrid->B3i[k+1][j][i]);
        }
      }
    }
  /* If mag. of total poloidal field (defined at cell center) is zero, then By is zero */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
          if (SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c) > TINY_NUMBER) {
		if(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c) > SQR(B0))
			pGrid->B2i[k][j][i] = 0.0;
		else
            		pGrid->B2i[k][j][i] = sqrt(SQR(B0)-(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c)));
		/* Make the net azimuthal flux to be zero */
/*		if(x3 < TINY_NUMBER)
			pGrid->B2i[k][j][i] = -pGrid->B2i[k][j][i];
*/
          } else {
            pGrid->B2i[k][j][i] = 0.0;

          }

		/* Fill in the mid-plane with uniform By */
		if(fabs(x3) < 0.8)
			pGrid->B2i[k][j][i] = sqrt(SQR(B0)-(SQR(pGrid->U[k][j][i].B1c)+SQR(pGrid->U[k][j][i].B3c)));

        }
      }
  /* Finally, calculate cell centered By field from face fields */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i]+pGrid->B2i[k][j+1][i]);
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
    dump_history_enroll(hst_T2,"<T^2>");
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

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
   dump_history_enroll(hst_sigmas,"<Sigma_es>");
   dump_history_enroll(hst_sigmaaP,"<Sigma_aP>");
#endif

   dump_history_enroll(hst_rho2,"<rho^2>");
   dump_history_enroll(hst_stressL,"StressL");
   dump_history_enroll(hst_stressR,"StressR");
   dump_history_enroll(hst_P,"Pg");
   dump_history_enroll(hst_rhofluxtop,"rhov_top");
   dump_history_enroll(hst_rhofluxbottom,"rhov_bottom");
   dump_history_enroll(hst_Erfluxtop,"Erflux_top");
   dump_history_enroll(hst_Erfluxbottom,"Erflux_bottom");
   dump_history_enroll(hst_gravpot,"GravPot");	


    frst = 0;
  }

	/* set boundary functions */

	/* use periodic boundary condition for gas variables */
	bvals_mhd_fun(pDomain, right_x3, radMHD_inflowke);

	bvals_rad_fun(pDomain, right_x3, radMHD_rad_inflowke);

	bvals_mhd_fun(pDomain, left_x3, radMHD_inflowks);

	bvals_rad_fun(pDomain, left_x3, radMHD_rad_inflowks);



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
	  pRG->R[ifr][k][j][i].J = Thermal_B(pGrid,pRG,ifr,ig,jg,kg);
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



	/* Write ghost zone */
	int m;
	for(m=0;m<nghost;m++){
		fwrite(&ghostrho[m],sizeof(Real),1,fp);
		fwrite(&ghostEr[m],sizeof(Real),1,fp);
		fwrite(&ghostFr[m],sizeof(Real),1,fp);
	}

	/* Write initial data */
	
	for(m=0;m<3;m++){
		fwrite(&inidata[m],sizeof(Real),1,fp);
		
	}


#ifdef RESISTIVITY
	fwrite(&eta_Ohm,sizeof(Real),1,fp);
        fwrite(&Jrhomax,sizeof(Real),1,fp);
#endif

    fwrite(&zbtm,sizeof(Real),1,fp);
    fwrite(&ztop,sizeof(Real),1,fp);

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



	/* read ghost zone */
	int m;
	for(m=0;m<nghost;m++){
		fread(&ghostrho[m],sizeof(Real),1,fp);
		fread(&ghostEr[m],sizeof(Real),1,fp);
		fread(&ghostFr[m],sizeof(Real),1,fp);
	}

	/* read  initial data */
	
	for(m=0;m<3;m++){
		fread(&inidata[m],sizeof(Real),1,fp);
		
	}

#ifdef RESISTIVITY
	fread(&eta_Ohm,sizeof(Real),1,fp);
	fread(&Jrhomax,sizeof(Real),1,fp);
	eta_Ohm = 1;
#endif

       fread(&zbtm,sizeof(Real),1,fp);
       fread(&ztop,sizeof(Real),1,fp);

/* enroll gravitational potential function */
#ifdef SHEARING_BOX
  ShearingBoxPot = UnstratifiedDisk;
  StaticGravPot = grav_vertical;
#endif

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
        /* enroll the opacity function */
        Opacity = Thindiskopacity;

#endif
/* enroll new history variables */

	
    dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
    dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
    dump_history_enroll(hst_T,"<T>");
    dump_history_enroll(hst_T2,"<T^2>");
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

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
   dump_history_enroll(hst_sigmas,"<Sigma_es>");
   dump_history_enroll(hst_sigmaaP,"<Sigma_aP>");
#endif

   dump_history_enroll(hst_rho2,"<rho^2>");

   dump_history_enroll(hst_stressL,"StressL");
   dump_history_enroll(hst_stressR,"StressR");
   dump_history_enroll(hst_P,"Pg");
   dump_history_enroll(hst_rhofluxtop,"rhov_top");
   dump_history_enroll(hst_rhofluxbottom,"rhov_bottom");
   dump_history_enroll(hst_Erfluxtop,"Erflux_top");
   dump_history_enroll(hst_Erfluxbottom,"Erflux_bottom");
   dump_history_enroll(hst_gravpot,"GravPot");


	

	
	/* Increase the background magnetic field */
	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	GridS *pGrid = pD->Grid;


	
	bvals_mhd_fun(pD, right_x3, radMHD_inflowke);

	bvals_rad_fun(pD, right_x3, radMHD_rad_inflowke);

	bvals_mhd_fun(pD, left_x3, radMHD_inflowks);

	bvals_rad_fun(pD, left_x3, radMHD_rad_inflowks);



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
*/
  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
  if(strcmp(expr,"dVy")==0) return expr_dV2;
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
  if(strcmp(name,"1d")==0) return output_1d;
  if(strcmp(name,"1dx")==0) return output_1dx;
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
        Real x1, x2, x3;
	dl = pG->dx1;
        if(pG->dx2 < dl) dl = pG->dx2;
        if(pG->dx3 < dl) dl = pG->dx3;

	if(pG->dt > TINY_NUMBER){	
		
        	eta0 = 0.03 * dl * dl/pG->dt;
	}	
	else{
		eta0 = 0.0;
	}
	Real factor = 0.2;
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
	
	if(fabs(x3) > 1.5){
		if(Jrhomax > 0.0){
			if(J > factor * Jrhomax){
				*eta_O = eta0 * J * J /(Jrhomax * Jrhomax);
			}
			else{
				*eta_O = eta0;
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
#endif



/*

 Real betafloor = 0.000;
	

	 GridS *pG;
        int i, j, k;
        int ie, is;
        int je, js;
        int ke, ks;
	int badcellflag;
        Real pressure, velocity, velocity_x, velocity_y, velocity_z, Bpre, beta;
	PrimS Wtemp;      
	pG = pM->Domain[0][0].Grid;
      

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

                         Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));

                         Bpre = 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);

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

 					}

                                }
#endif

                        }
                }
        }

*/
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

static Real grav_vertical(const Real x1, const Real x2, const Real x3)
{
	/* In the ghost zeons, set the gravitational potential to be flat */
	Real h;

	if(x3 > ztop)
		h = ztop;
	else if(x3 < zbtm)
		h = zbtm;
	else
		h = x3;		

  return 0.5 * Omega_0 * Omega_0 * h * h;
}

static Real hst_gravpot(const GridS *pG,const int i,const int j, const int k)
{
	Real x1, x2, x3;
	cc_pos(pG,i,j,k,&x1, &x2, &x3);

  return pG->U[k][j][i].d * grav_vertical(x1,x2,x3);
}


#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
void Thindiskopacity(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY])
{
/* Sigma[0-NOPACITY] are: Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE respectively */
/* dSigma[0-2*NOPACITY] are: dSigma_sF/drho, dSigma_aF/drho, dSigma_aP/drho, dSigma_aE/drho */
/* 			     dSigma_sF/dT,   dSigma_aF/dT,   dSigma_aP/dT,   dSigma_aE/dT */
	
/* When pressure becomes negative, we do not include radiation source term */
if((rho * T * R_ideal > 2.0 * TINY_NUMBER) && (rho > 0.0)){	
	Real Tpower, Tpower1;
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
#ifdef FARGO
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#else
  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d + qshear*Omega_0*x1);
#endif
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
  Real x1,x2,x3,phi, Kvy, Vy0, dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
  phi = UnstratifiedDisk(x1, x2, x3);
  Vy0 = -qshear * Omega_0 * x1;
  Kvy = 0.0; 
  dVy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
#ifdef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
  Kvy = 0.5 * pG->U[k][j][i].d * (2.0 * dVy * Vy0 + Vy0 * Vy0);
#endif

  return pG->U[k][j][i].E + pG->U[k][j][i].d*phi + Kvy;
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

static Real hst_stressL(const GridS *pG, const int i, const int j, const int k)
{
	Real rhov, BxBy;
	rhov = 0.0;
	BxBy = 0.0;

	int ID, lx1_id, rx1_id;

#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	lx1_id = pG->lx1_id;
	rx1_id = pG->rx1_id;
#else
	ID = 0;
	lx1_id = -1;
	rx1_id = -1;
#endif

	/* If it is on the boundary */
	if(((i == pG->is) && ((lx1_id < 0) || (lx1_id >= ID))))
	{
		rhov = hst_rho_Vx_dVy(pG,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
		BxBy =  hst_BxBy(pG, i, j, k);
#endif	

		return (qshear*Omega_0*(rhov + BxBy)* pG->dx2 * pG->dx3);
	}	
	else
		return 0.0;
}


static Real hst_stressR(const GridS *pG, const int i, const int j, const int k)
{
	Real rhov, BxBy;
	rhov = 0.0;
	BxBy = 0.0;

	int ID, lx1_id, rx1_id;

#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	lx1_id = pG->lx1_id;
	rx1_id = pG->rx1_id;
#else
	ID = 0;
	lx1_id = -1;
	rx1_id = -1;
#endif

	/* If it is on the boundary */
	if(((i == pG->ie) && ((rx1_id < 0) || (rx1_id <= ID))))
	{
		rhov = hst_rho_Vx_dVy(pG,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
		BxBy =  hst_BxBy(pG, i, j, k);
#endif	

		return (qshear*Omega_0*(rhov + BxBy)* pG->dx2 * pG->dx3);

	}
	else
		return 0.0;
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

/* Function for transfer module */
#ifdef RADIATION_TRANSFER

static Real Thermal_B(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
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

static Real Transfereps(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
		      const int k)
{
	Real eps;
	eps = pG->U[k][j][i].Sigma[1] / (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[1]);

  return eps;
  
}

static Real transfer_opacity(const GridS *pG, const RadGridS *pRG, const int ifr, const int i, const int j, 
			  const int k)
{
	
  return (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[1]);
  
}


void Userwork_in_formal_solution(DomainS *pD)
{


}

#endif




/* Function for inflow boundary condition */


void radMHD_inflowke(GridS *pGrid)
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
			pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke+k-1][j][i];
			pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke+k-1][j][i];

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
				pGrid->U[ke+k][j][i].B1c = pGrid->U[ke][j][i].B1c;
				pGrid->U[ke+k][j][i].B2c = pGrid->U[ke][j][i].B2c;
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
		if(pGrid->U[ke+k][j][i].d > dfloor)
			pGrid->U[ke+k][j][i].d = dfloor;				
		

		
      		pGrid->U[ke+k][j][i].M1 = velocity1 * pGrid->U[ke+k][j][i].d;
		pGrid->U[ke+k][j][i].M2 = velocity2 * pGrid->U[ke+k][j][i].d;
			
		
		pGrid->U[ke+k][j][i].M3 = velocity3 * pGrid->U[ke+k][j][i].d;

		pressure = T * R_ideal * pGrid->U[ke+k][j][i].d;

		pGrid->U[ke+k][j][i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[ke+k][j][i].M1 * pGrid->U[ke+k][j][i].M1 + pGrid->U[ke+k][j][i].M2 * pGrid->U[ke+k][j][i].M2 + pGrid->U[ke+k][j][i].M3 * pGrid->U[ke+k][j][i].M3) / pGrid->U[ke+k][j][i].d;



#ifdef RADIATION_MHD
		pGrid->U[ke+k][j][i].E += 0.5 * (pGrid->U[ke+k][j][i].B1c * pGrid->U[ke+k][j][i].B1c + pGrid->U[ke+k][j][i].B2c * pGrid->U[ke+k][j][i].B2c + pGrid->U[ke+k][j][i].B3c * pGrid->U[ke+k][j][i].B3c);

#endif

	
			
		Thindiskopacity(pGrid->U[ke+k][j][i].d, T, Sigma, NULL);
		
		for(m=0; m<NOPACITY;m++){
			pGrid->U[ke+k][j][i].Sigma[m] = Sigma[m];
		}	
	
      		
      		}
    		}
	}
  
}


void radMHD_rad_inflowke(GridS *pGrid)
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
	Real dz = pGrid->dx3;

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



void radMHD_inflowks(GridS *pGrid)
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
				pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks-k+1][j][i];
				pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks-k+1][j][i];
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
				pGrid->U[ks-k][j][i].B1c = pGrid->U[ks][j][i].B1c;
				pGrid->U[ks-k][j][i].B2c = pGrid->U[ks][j][i].B2c;
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

		if(pGrid->U[ks-k][j][i].d > dfloor)
			pGrid->U[ks-k][j][i].d = dfloor;

		
      		pGrid->U[ks-k][j][i].M1 = velocity1 * pGrid->U[ks-k][j][i].d;
		pGrid->U[ks-k][j][i].M2 = velocity2 * pGrid->U[ks-k][j][i].d;
			
		
		pGrid->U[ks-k][j][i].M3 = velocity3 * pGrid->U[ks-k][j][i].d;
		
		pressure = T * R_ideal * pGrid->U[ks-k][j][i].d;

		pGrid->U[ks-k][j][i].E =  pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[ks-k][j][i].M1 * pGrid->U[ks-k][j][i].M1 + pGrid->U[ks-k][j][i].M2 * pGrid->U[ks-k][j][i].M2 + pGrid->U[ks-k][j][i].M3 * pGrid->U[ks-k][j][i].M3) / pGrid->U[ks-k][j][i].d;



#ifdef RADIATION_MHD
		pGrid->U[ks-k][j][i].E += 0.5 * (pGrid->U[ks-k][j][i].B1c * pGrid->U[ks-k][j][i].B1c + pGrid->U[ks-k][j][i].B2c * pGrid->U[ks-k][j][i].B2c + pGrid->U[ks-k][j][i].B3c * pGrid->U[ks-k][j][i].B3c);

#endif

			
		Thindiskopacity(pGrid->U[ks-k][j][i].d, T, Sigma, NULL);
		
		for(m=0; m<NOPACITY;m++){
			pGrid->U[ks-k][j][i].Sigma[m] = Sigma[m];
		}	
	
      		
      }
    }
  
  
}


void radMHD_rad_inflowks(GridS *pGrid)
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
	Real dz = pGrid->dx3;

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

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
#ifdef MATRIX_MULTIGRID

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

	*Mat_BCFun = radMHD_Mat_inflowks;

} 
void bvals_mat_fun_ox3(VMatFun_t *Mat_BCFun)
{

	*Mat_BCFun = radMHD_Mat_inflowke;

} 




void radMHD_Mat_inflowke(MatrixS *pMat)
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
	    		
				Sigma_t = 0.5 * (pMat->Ugas[ke+k][j][i].Sigma[0] + pMat->Ugas[ke+k][j][i].Sigma[1] + pMat->Ugas[ke+k-1][j][i].Sigma[0] + pMat->Ugas[ke+k-1][j][i].Sigma[1]);

				pMat->Ugas[ke+k][j][i].Edd_11 = pMat->Ugas[ke][j][i].Edd_11;
				pMat->Ugas[ke+k][j][i].Edd_22 = pMat->Ugas[ke][j][i].Edd_22;
				pMat->Ugas[ke+k][j][i].Edd_21 = pMat->Ugas[ke][j][i].Edd_21;
				pMat->Ugas[ke+k][j][i].Edd_31 = pMat->Ugas[ke][j][i].Edd_31;
				pMat->Ugas[ke+k][j][i].Edd_32 = pMat->Ugas[ke][j][i].Edd_32;
				pMat->Ugas[ke+k][j][i].Edd_33 = pMat->Ugas[ke][j][i].Edd_33;
			if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0]) ){
				pMat->U[ke+k][j][i].Er  = 0.0;
				pMat->U[ke+k][j][i].Fr1 = 0.0;
				pMat->U[ke+k][j][i].Fr2 = 0.0;
				pMat->U[ke+k][j][i].Fr3 = 0.0;
			}
			else{
				pMat->U[ke+k][j][i].Fr1 = pMat->U[ke][j][i].Fr1;
				pMat->U[ke+k][j][i].Fr2 = pMat->U[ke][j][i].Fr2;
			/*	if(pMat->U[ke][j][i].Fr3 > 0.0)
					pMat->U[ke+k][j][i].Fr3 = pMat->U[ke][j][i].Fr3;
				else
					pMat->U[ke+k][j][i].Fr3 = 0.0;
			*/
				vx = pMat->Ugas[ke+k-1][j][i].V1;
				vy = pMat->Ugas[ke+k-1][j][i].V2;
				vz = pMat->Ugas[ke+k-1][j][i].V3;

				matrix_alpha(0.0, pMat->Ugas[ke+k][j][i].Sigma, pMat->dt, pMat->Ugas[ke+k][j][i].Edd_33, 0.0, &reducefactor, 0, pMat->dx3);

/*				tau = pMat->dt * Crat * (pMat->Ugas[ke+k][j][i].Sigma[0] + pMat->Ugas[ke+k][j][i].Sigma[1]);
				tau = tau * tau / (2.0 * pMat->Ugas[ke+k][j][i].Edd_33);

				if(tau > 0.001)
					reducefactor = sqrt(pMat->Ugas[ke+k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
				else
					reducefactor = sqrt(pMat->Ugas[ke+k][j][i].Edd_33 * (1.0 - 0.5 * tau));	
*/

				Fr0z = pMat->U[ke+k-1][j][i].Fr3 - ((1.0 + pMat->Ugas[ke+k-1][j][i].Edd_33) * vz + pMat->Ugas[ke+k-1][j][i].Edd_31 * vx + pMat->Ugas[ke+k-1][j][i].Edd_32 * vy)* pMat->U[ke+k-1][j][i].Er / Crat;

				Eratio = pMat->Ugas[ke+k][j][i].Edd_33 + 0.5 * pMat->dx3 * reducefactor * Sigma_t;

		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
				pMat->U[ke+k][j][i].Er = (pMat->Ugas[ke+k-1][j][i].Edd_33 * pMat->U[ke+k-1][j][i].Er - 0.5 * pMat->dx3 * Sigma_t * Fr0z ) / Eratio;

				if(pMat->U[ke+k][j][i].Er < 0.0){
					Eratio = pMat->Ugas[ke+k][j][i].Edd_33 + pMat->dx3 * reducefactor * Sigma_t;
				pMat->U[ke+k][j][i].Er = pMat->Ugas[ke+k-1][j][i].Edd_33 * pMat->U[ke+k-1][j][i].Er / Eratio;
				}
		
				Fr0z = reducefactor * pMat->U[ke+k][j][i].Er;
				
				pMat->U[ke+k][j][i].Fr3 = Fr0z + ((1.0 + pMat->Ugas[ke+k][j][i].Edd_33) * pMat->Ugas[ke+k][j][i].V3 + pMat->Ugas[ke+k][j][i].Edd_31 * pMat->Ugas[ke+k][j][i].V1 + pMat->Ugas[ke+k][j][i].Edd_32 * pMat->Ugas[ke+k][j][i].V2)* pMat->U[ke+k][j][i].Er / Crat;


			/*		

				if(Fr0z > 0.0){			
					pMat->U[ke+k][j][i].Er = pMat->U[ke+k-1][j][i].Er - pMat->dx3 * Sigma_t * Fr0z / pMat->U[ke+k][j][i].Edd_33;			

					if(pMat->U[ke+k][j][i].Er < 0.0){
					pMat->U[ke+k][j][i].Er = pMat->U[ke+k-1][j][i].Er;

				}

			
				}
				else{
					pMat->U[ke+k][j][i].Er = pMat->U[ke+k-1][j][i].Er;
			
				}
		*/



			}


    			}
		}
	}

  	return;
}

void radMHD_Mat_inflowks(MatrixS *pMat)
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
	    		
				Sigma_t = 0.5 * (pMat->Ugas[ks-k][j][i].Sigma[0] + pMat->Ugas[ks-k][j][i].Sigma[1] + pMat->Ugas[ks-k+1][j][i].Sigma[0] + pMat->Ugas[ks-k+1][j][i].Sigma[1]);

				pMat->Ugas[ks-k][j][i].Edd_11 = pMat->Ugas[ks][j][i].Edd_11;
				pMat->Ugas[ks-k][j][i].Edd_22 = pMat->Ugas[ks][j][i].Edd_22;
				pMat->Ugas[ks-k][j][i].Edd_21 = pMat->Ugas[ks][j][i].Edd_21;
				pMat->Ugas[ks-k][j][i].Edd_31 = pMat->Ugas[ks][j][i].Edd_31;
				pMat->Ugas[ks-k][j][i].Edd_32 = pMat->Ugas[ks][j][i].Edd_32;
				pMat->Ugas[ks-k][j][i].Edd_33 = pMat->Ugas[ks][j][i].Edd_33;

				
			if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){
				pMat->U[ks-k][j][i].Er  = 0.0;
				pMat->U[ks-k][j][i].Fr1 = 0.0;
				pMat->U[ks-k][j][i].Fr2 = 0.0;
				pMat->U[ks-k][j][i].Fr3 = 0.0;
			}
			else{
				pMat->U[ks-k][j][i].Fr1 = pMat->U[ks][j][i].Fr1;
				pMat->U[ks-k][j][i].Fr2 = pMat->U[ks][j][i].Fr2;
			/*	if(pMat->U[ke][j][i].Fr3 > 0.0)
					pMat->U[ks-k][j][i].Fr3 = pMat->U[ke][j][i].Fr3;
				else
					pMat->U[ks-k][j][i].Fr3 = 0.0;
			*/
				vx = pMat->Ugas[ks-k+1][j][i].V1;
				vy = pMat->Ugas[ks-k+1][j][i].V2;
				vz = pMat->Ugas[ks-k+1][j][i].V3;

				matrix_alpha(0.0, pMat->Ugas[ks-k][j][i].Sigma, pMat->dt, pMat->Ugas[ks-k][j][i].Edd_33, 0.0, &reducefactor, 0, pMat->dx3);
			
		/*		tau = pMat->dt * Crat * (pMat->Ugas[ks-k][j][i].Sigma[0] + pMat->Ugas[ks-k][j][i].Sigma[1]);
				tau = tau * tau / (2.0 * pMat->Ugas[ks-k][j][i].Edd_33);

				if(tau > 0.001)
					reducefactor = sqrt(pMat->Ugas[ks-k][j][i].Edd_33 * (1.0 - exp(- tau)) / tau);
				else
					reducefactor = sqrt(pMat->Ugas[ks-k][j][i].Edd_33 * (1.0 - 0.5 * tau));	
		*/

				Fr0z = pMat->U[ks-k+1][j][i].Fr3 - ((1.0 + pMat->Ugas[ks-k+1][j][i].Edd_33) * vz + pMat->Ugas[ks-k+1][j][i].Edd_31 * vx + pMat->Ugas[ks-k+1][j][i].Edd_32 * vy)* pMat->U[ks-k+1][j][i].Er / Crat;

				Eratio = pMat->Ugas[ks-k][j][i].Edd_33 + 0.5 * pMat->dx3 * reducefactor * Sigma_t;

		/* vacuum boundary condition. Fr0z = R * Er, where R is the reduced speed of light */
		/* Assume this relation holds in the ghost zones */
		/*  (f1 * Er(K+1, n+1) - f * Er(k,n))/(Sigmat * dx3) = -0.5*(Fr0z(k) + Fr0z(k+1)), and Fr0z(k+1)= R * Er(k+1) */
		
				pMat->U[ks-k][j][i].Er = (pMat->Ugas[ks-k+1][j][i].Edd_33 * pMat->U[ks-k+1][j][i].Er - 0.5 * pMat->dx3 * Sigma_t * Fr0z ) / Eratio;

				if(pMat->U[ks-k][j][i].Er < 0.0){
					Eratio = pMat->Ugas[ks-k][j][i].Edd_33 + pMat->dx3 * reducefactor * Sigma_t;
				pMat->U[ks-k][j][i].Er = pMat->Ugas[ks-k+1][j][i].Edd_33 * pMat->U[ks-k+1][j][i].Er / Eratio;
				}
		
				Fr0z = reducefactor * pMat->U[ks-k][j][i].Er;
				
				pMat->U[ks-k][j][i].Fr3 = Fr0z + ((1.0 + pMat->Ugas[ks-k][j][i].Edd_33) * pMat->Ugas[ks-k][j][i].V3 + pMat->Ugas[ks-k][j][i].Edd_31 * pMat->Ugas[ks-k][j][i].V1 + pMat->Ugas[ks-k][j][i].Edd_32 * pMat->Ugas[ks-k][j][i].V2)* pMat->U[ks-k][j][i].Er / Crat;

					
			}


    			}
		}
	}

  	return;
}


#endif /* end multi_grid */

#endif /* end radMHD or radhydro */




/*! \fn static void output_1d(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
    averaged quantities.  Currently, only outputs at lowest
    refinement level */

static void output_1d(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
  int i,j,k;
  int tot1d,i1d,nzmx,my_nz,kg,kdisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lx,Ly,press, press1, press3, Bpre1, Bpre3;
  static double *out_x3;
  double vx, vy, vz, Fr01,Fr02,Fr03;

  FILE *p_1dfile;
  char *fname;
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int zproc;
  int ierr,myID_Comm_Domain;
#endif

/* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap */

#if defined(MHD) || defined(RADIATION_MHD)
  tot1d=15+6+5+1;
#else
  tot1d=7+6+5+1;
#endif /* MHD */
#ifdef ADIABATIC
  tot1d=tot1d+3+1;
#endif /* ADIABATIC */

  Lx = pM->RootMaxX[0] - pM->RootMinX[0];
  Ly = pM->RootMaxX[1] - pM->RootMinX[1];
  nzmx = pM->Nx[2];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef MPI_PARALLEL
  int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->Disp[2];

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (k=0; k<nzmx; k++) {
      x3 = pM->RootMinX[2] + (k + 0.5)*pGrid->dx3;
      out_x3[k] = x3;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (k=ks; k<=ke; k++) {
    kg=k+kdisp-nghost;
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        i1d=0;
        out1d[kg][i1d] += pGrid->U[k][j][i].d;
        i1d++;
#ifdef ISOTHERMAL
        out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k)
#endif
                                ),TINY_NUMBER);
        out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[kg][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
        i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].E;
        i1d++;
        out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        i1d++;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
#ifdef FARGO
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
#else
        out1d[kg][i1d] += 0.5*pGrid->U[k][j][i].d*SQR(pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d + qshear*Omega_0*x1);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
        i1d++;
        out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
        i1d++;
        out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bx(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_By(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bz(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Er;
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Fr1;
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Fr2;
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Fr3;
	i1d++;
	vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
	vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
#ifdef FARGO
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;


	Fr02 = pGrid->U[k][j][i].Fr3 -(vz * (1.0 + pGrid->U[k][j][i].Edd_33) + vx * pGrid->U[k][j][i].Edd_31 + vy * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;

	out1d[kg][i1d] += Fr02;
	i1d++;

	/* dFr0/dz */
	vx = pGrid->U[k-1][j][i].M1 / pGrid->U[k-1][j][i].d;
	vy = pGrid->U[k-1][j][i].M2 / pGrid->U[k-1][j][i].d;
#ifdef FARGO
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k-1][j][i].M3 / pGrid->U[k-1][j][i].d;
	Fr01 = pGrid->U[k-1][j][i].Fr3 -(vz * (1.0 + pGrid->U[k-1][j][i].Edd_33) + vx * pGrid->U[k-1][j][i].Edd_31 + vy * pGrid->U[k-1][j][i].Edd_32) * pGrid->U[k-1][j][i].Er/Crat;

	vx = pGrid->U[k+1][j][i].M1 / pGrid->U[k+1][j][i].d;
	vy = pGrid->U[k+1][j][i].M2 / pGrid->U[k+1][j][i].d;
#ifdef FARGO
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k+1][j][i].M3 / pGrid->U[k+1][j][i].d;
	Fr03 = pGrid->U[k+1][j][i].Fr3 -(vz * (1.0 + pGrid->U[k+1][j][i].Edd_33) + vx * pGrid->U[k+1][j][i].Edd_31 + vy * pGrid->U[k+1][j][i].Edd_32) * pGrid->U[k+1][j][i].Er/Crat;
	out1d[kg][i1d] += (Fr03 - Fr01) * 0.5 / pGrid->dx3;
	i1d++;
	
	/* Er * vz */
	out1d[kg][i1d] += pGrid->U[k][j][i].Er * pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d;
	i1d++;

	/* dPg/dz / rho */
	press1           = MAX(Gamma_1*(pGrid->U[k-1][j][i].E - expr_KE(pGrid,i,j,k-1)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k-1)
#endif
                                ),TINY_NUMBER);
	press3           = MAX(Gamma_1*(pGrid->U[k+1][j][i].E - expr_KE(pGrid,i,j,k+1)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k+1)
#endif
                                ),TINY_NUMBER);
	out1d[kg][i1d] += (press3 - press1) * 0.5 / (pGrid->dx3 * pGrid->U[k][j][i].d);
	i1d++;

	/* dBpre/dz / rho */
	Bpre1 = expr_ME(pGrid,i,j,k-1);
	Bpre3 = expr_ME(pGrid,i,j,k+1);
	out1d[kg][i1d] += (Bpre3 - Bpre1) * 0.5 / (pGrid->dx3 * pGrid->U[k][j][i].d);
	i1d++;

        out1d[kg][i1d] += hst_sigmas(pGrid,i,j,k);
	i1d++;
        out1d[kg][i1d] += hst_sigmaaP(pGrid,i,j,k);
	/* dPrdz/sigm */
	i1d++;
	out1d[kg][i1d] += (pGrid->U[k+1][j][i].Er * pGrid->U[k+1][j][i].Edd_33 - pGrid->U[k-1][j][i].Er * pGrid->U[k-1][j][i].Edd_33) / (2.0* pGrid->dx3 * (pGrid->U[k][j][i].Sigma[0] + pGrid->U[k][j][i].Sigma[1]));

#endif

      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lx*Ly/(pGrid->dx1*pGrid->dx2);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1d");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
#ifdef MHD
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6]);
#endif /* MHD */
#else
#ifdef RADIATION_MHD
    if (k == 0) {
      fprintf(p_1dfile,"# [1]x3     [2]dens    [3]pressure    [4]temperature  [5]E     [6]Etot     [7]KEx         [8]KEy        [9] KEz       [10] KE        [11]Reynolds   [12]MEx        [13]MEy        [14]MEz        [15]ME         [16]Bx          [17]By         [18]Bz         [19]Maxwell     [20]Er      [21]Frx      [22]Fry       [23]Frz     [24]Frz0	[25]dFr0dz	[26]ErV		[27]dP/dz/rho		[28]dBpre/dz/rho	[29]kappaes          [30]kappaff	[31]dPrdz/sigma\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28],out1d[k][29]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens    pressure    temperature  E     Etot     KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9]);
#endif /* RADIATION_MHD */
#endif /* ISOTHERMAL */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
  if (FIRST == 0) {
    FIRST = 1;
  }

return;
}



/*! \fn static void output_1dx(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
    averaged quantities.  Currently, only outputs at lowest
    refinement level */

static void output_1dx(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
  int i,j,k;
  int tot1d,i1d,nzmx,my_nz,kg,kdisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lx,Ly,Lz,press, press1, press3, Bpre1, Bpre3;
  static double *out_x3;
  double vx, vy, vz, Fr01,Fr02,Fr03;
  int flag;

  FILE *p_1dfile;
  char *fname;
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int zproc;
  int ierr,myID_Comm_Domain;
#endif

/* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap */

#if defined(MHD) || defined(RADIATION_MHD)
  tot1d=15+6+5;
#else
  tot1d=7+6+5;
#endif /* MHD */
#ifdef ADIABATIC
  tot1d=tot1d+3;
#endif /* ADIABATIC */

  Lx = pM->RootMaxX[0] - pM->RootMinX[0];
  Ly = pM->RootMaxX[1] - pM->RootMinX[1];
  Lz = pM->RootMaxX[2] - pM->RootMinX[2];
  nzmx = pM->Nx[0];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef MPI_PARALLEL
  int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->Disp[0];

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (k=0; k<nzmx; k++) {
      x1 = pM->RootMinX[0] + (k + 0.5)*pGrid->dx1;
      out_x3[k] = x1;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (i=is; i<=ie; i++) {
    kg=i+kdisp-nghost;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        i1d=0;
        out1d[kg][i1d] += pGrid->U[k][j][i].d;
        i1d++;
#ifdef ISOTHERMAL
        out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k)
#endif
                                ),TINY_NUMBER);
        out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[kg][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
        i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].E;
        i1d++;
        out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        i1d++;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	if(x3 >= 0.0)	flag = 1;
	else	flag = -1;
#ifdef FARGO
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
#else
        out1d[kg][i1d] += 0.5*pGrid->U[k][j][i].d*SQR(pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d + qshear*Omega_0*x1);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
        i1d++;
        out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
        i1d++;
        out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bx(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_By(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bz(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Er;
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Fr1;
	i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].Fr2;
	i1d++;
	/* To avoid cancel */
	
        out1d[kg][i1d] += flag * pGrid->U[k][j][i].Fr3;
	
	i1d++;
	vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
	vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
#ifdef FARGO
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;


	Fr02 = pGrid->U[k][j][i].Fr3 -(vz * (1.0 + pGrid->U[k][j][i].Edd_33) + vx * pGrid->U[k][j][i].Edd_31 + vy * pGrid->U[k][j][i].Edd_32) * pGrid->U[k][j][i].Er/Crat;

	
	out1d[kg][i1d] += flag * Fr02;
	
	i1d++;

	/* dFr0/dz */
	vx = pGrid->U[k-1][j][i].M1 / pGrid->U[k-1][j][i].d;
	vy = pGrid->U[k-1][j][i].M2 / pGrid->U[k-1][j][i].d;
#ifdef FARGO
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k-1][j][i].M3 / pGrid->U[k-1][j][i].d;
	Fr01 = pGrid->U[k-1][j][i].Fr3 -(vz * (1.0 + pGrid->U[k-1][j][i].Edd_33) + vx * pGrid->U[k-1][j][i].Edd_31 + vy * pGrid->U[k-1][j][i].Edd_32) * pGrid->U[k-1][j][i].Er/Crat;

	vx = pGrid->U[k+1][j][i].M1 / pGrid->U[k+1][j][i].d;
	vy = pGrid->U[k+1][j][i].M2 / pGrid->U[k+1][j][i].d;
#ifdef FARGO
	vy -= qshear * Omega_0 * x1;
#endif
	vz = pGrid->U[k+1][j][i].M3 / pGrid->U[k+1][j][i].d;
	Fr03 = pGrid->U[k+1][j][i].Fr3 -(vz * (1.0 + pGrid->U[k+1][j][i].Edd_33) + vx * pGrid->U[k+1][j][i].Edd_31 + vy * pGrid->U[k+1][j][i].Edd_32) * pGrid->U[k+1][j][i].Er/Crat;
	
	out1d[kg][i1d] += flag * (Fr03 - Fr01) * 0.5 / pGrid->dx3;
	
	i1d++;
	
	/* Er * vz */
	out1d[kg][i1d] += flag * pGrid->U[k][j][i].Er * pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d;
	i1d++;

	/* dPg/dz / rho */
	press1           = MAX(Gamma_1*(pGrid->U[k-1][j][i].E - expr_KE(pGrid,i,j,k-1)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k-1)
#endif
                                ),TINY_NUMBER);
	press3           = MAX(Gamma_1*(pGrid->U[k+1][j][i].E - expr_KE(pGrid,i,j,k+1)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k+1)
#endif
                                ),TINY_NUMBER);
	out1d[kg][i1d] += flag * (press3 - press1) * 0.5 / (pGrid->dx3 * pGrid->U[k][j][i].d);
	i1d++;

	/* dBpre/dz / rho */
	Bpre1 = expr_ME(pGrid,i,j,k-1);
	Bpre3 = expr_ME(pGrid,i,j,k+1);
	out1d[kg][i1d] += flag * (Bpre3 - Bpre1) * 0.5 / (pGrid->dx3 * pGrid->U[k][j][i].d);
	i1d++;

        out1d[kg][i1d] += hst_sigmas(pGrid,i,j,k);
	i1d++;
        out1d[kg][i1d] += hst_sigmaaP(pGrid,i,j,k);

#endif

      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lz*Ly/(pGrid->dx3*pGrid->dx2);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dx");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dx");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
#ifdef MHD
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6]);
#endif /* MHD */
#else
#ifdef RADIATION_MHD
    if (k == 0) {
      fprintf(p_1dfile,"# [1]x3     [2]dens    [3]pressure    [4]temperature  [5]E     [6]Etot     [7]KEx         [8]KEy        [9] KEz       [10] KE        [11]Reynolds   [12]MEx        [13]MEy        [14]MEz        [15]ME         [16]Bx          [17]By         [18]Bz         [19]Maxwell     [20]Er      [21]Frx      [22]Fry       [23]Frz     [24]Frz0	[25]dFr0dz	[26]ErV		[27]dP/dz/rho		[28]dBpre/dz/rho	[29]kappaes          [30]kappaff\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens    pressure    temperature  E     Etot     KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9]);
#endif /* RADIATION_MHD */
#endif /* ISOTHERMAL */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
  if (FIRST == 0) {
    FIRST = 1;
  }

return;
}





