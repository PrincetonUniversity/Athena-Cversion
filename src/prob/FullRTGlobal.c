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

#ifndef CYLINDRICAL

#error FullRTGlobal.c requires CYLINDRICAL coordinate enabled!

#endif


Real Lx,Ly,Lz; /* root grid size, global to share with output functions */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static double ran2(long int *idum);


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


static Real expr_beta(GridS *pG, const int i, const int j, const int k);
static Real expr_ME(GridS *pG, const int i, const int j, const int k);
static Real expr_KE(GridS *pG, const int i, const int j, const int k);
static Real hst_rho_VrVphi(const GridS *pG,const int i,const int j,const int k);
#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif

#ifdef FULL_RADIATION_TRANSFER
static Real hst_sigmas(const GridS *pG, const int i, const int j, const int k);
static Real hst_sigmaaP(const GridS *pG, const int i, const int j, const int k);
#endif

static Real hst_gravpot(const GridS *pG, const int i, const int j, const int k);
static Real hst_Lx(const GridS *pG, const int i, const int j, const int k);
static Real hst_Ly(const GridS *pG, const int i, const int j, const int k);
static Real hst_Lz(const GridS *pG, const int i, const int j, const int k);

#ifdef FULL_RADIATION_TRANSFER
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_By(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bz(const GridS *pG, const int i, const int j, const int k);
static Real hst_BrBphi(const GridS *pG, const int i, const int j, const int k);
static Real hst_EB(const GridS *pG, const int i, const int j, const int k);
#endif /* MHD */
static Real hst_T(const GridS *pG, const int i, const int j, const int k);

static Real hst_P(const GridS *pG, const int i, const int j, const int k);

static Real hst_rho2(const GridS *pG, const int i, const int j, const int k);

static Real hst_T2(const GridS *pG, const int i, const int j, const int k);


/* Functions to write azimuthal averaged 1D and 2D functions */
static void output_1d(MeshS *pM, OutputS *pOut);
static void output_2d_binary(MeshS *pM, OutputS *pOut);


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

#ifdef FULL_RADIATION_TRANSFER
static void Thindiskopacity(GridS *pG, const int ifr, const int i, const int j, const int k, Real *Sigma);
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

static Real Tfloor = 1.e-4;
static Real dfloor = 5.e-6;

/* The side of the mask region Rmas = 2 in unit of r_g */
static Real R0 = 30.0;



/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
/* If the problem generator is called, Domain.Grid must be not NULL */
void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  RadGridS *pRG = (pDomain->RadGrid);
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;

  int isr = pRG->is, ier = pRG->ie;
  int jsr = pRG->js, jer = pRG->je;
  int ksr = pRG->ks, ker = pRG->ke; 
  int nf=pRG->nf, nang=pRG->nang;
  int noct = pRG->noct;
  int i, j, k, ifr, l, n, m;
  int ig, jg, kg;
  int offset;
	
  const Real *r = pGrid->r, *ri=pGrid->ri;
  Real rsf, lsf;
	
  int torusflag = 1;
  /* when torusflag =1, use adiabatic initial condition */
  /* When torusflag =0, use isothermal initial condition */
	

  int ixs,jxs,kxs,ipert,ifield;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax;
  Real Radius, distance, Vc, costheta, sintheta;
  Real den = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real density, pressure, Er, vx, vy, Frx, Fry, Frz, temperature;
  Real T0,Tr, coef1, coef2, coef3;
  Real factor = 1.0;
  Real beta=1.0,B0,kx,ky,kz,amp;
  Real Bamp;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  ifield = par_geti_def("problem","ifield", 6);
  ifield = 6;
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */

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
	

	/* Parse global variables of unit ratio */
#ifdef FULL_RADIATION_TRANSFER
	Prat = par_getd("problem","Pratio");
	Crat = par_getd("problem","Cratio");
	R_ideal = par_getd("problem","R_ideal");	

#ifdef MPI_PARALLEL 
  if(myID_Comm_world == 0){
#endif
	/* Print out the parameters for consistency check */
     printf("Parameters: Prat %G Crat %G R_ideal %G\n",Prat,Crat,R_ideal);
#ifdef MPI_PARALLEL 
  }
#endif

#endif
	
	betaz=0.0;
	betay=40.0;
	pres = 1.0;

  	B0z = 0.0;
  	B0y = sqrt((double)(2.0*pres/betay));
	B0 = sqrt(B0z * B0z + B0y * B0y);	

	kappaes = 6.48511e3;
        kappaffP = 177.58;
        kappaffR = 4.79946;

	/* Initialize temperature unit */
/*
   This Tunit is only needed when Compton term is added
*/
#ifdef FULL_RADIATION_TRANSFER
	Tunit= 1.e7;
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
	Real vs0 = 3.0;
	Real rho0 = 30.0;
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
		
	if(torusflag == 1){	
		tempphi = ((Effphi-Effphi0)/nindex)/(vs0*vs0);
		if((abs(tempphi) < 1.0) && (distance > 10.0)){
			amp = 0.01;
		
		
			density = rho0 * pow((1.0-tempphi),nindex);
			/* random perturbation */
			
			if(density < dfloor){
				density = dfloor;
			}
			
			pressure = rho0*vs0*vs0*pow(density/rho0,Gamma)/Gamma;
		
		
			

		}
		else{
                vx = 0.0;
                vy = 0.0;
                amp = 0.0;

                Effphi = -Crat*Crat/(2.0*(distance-1.0));
                Effphi0 = -Crat*Crat/(2.0*(sqrt(SQR(Lx)+SQR(Lz))-1.0));
                density = dfloor;
                temperature = -(Effphi - Effphi0) + 10.0 * Tfloor;

				temperature = 1.0;
                pressure = density  * temperature;

		
		}
	}/* end torusflag = 1 */
	else {
		tempphi = ((Effphi - Effphi0))/(vs0*vs0);
		density = rho0 * exp(-tempphi);
		
		if(density < dfloor)
			density = dfloor;
		
				
		pressure = density * vs0 * vs0;
	}

		
	T0 = pressure/(density*R_ideal);
	coef1 = Prat/3.0;
	coef2 = density * R_ideal;
	coef3 = -pressure;
	
/*	if(pressure < 5.02e-6)
		printf("P: %e Prat: %e T: %e\n",pressure,Prat,temperature);
*/

	
	temperature = rtsafe(Tequilibrium, 0.0, T0, 1.e-12, coef1, coef2, coef3,0.0);
	pressure = temperature * density * R_ideal;	
		
		

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
        rd = density*(1.0 + 2.0*rval);
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
#ifdef MHD

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
		
    }
  }}

#ifdef MHD
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

	Radius = x1; 

	if((pGrid->U[k][j][i].d > 10.0 *dfloor) && (Radius > 10.0))
		Aphi = B0 * pGrid->U[k][j][i].d;
	else
		Aphi = 0.0;

		/* Now convert to Ax and Ay */
		Ay[k][j][i] = Aphi;
	
        }
      }
    }
  /* In this case, we calculate face fields from vect. potential */
    for (k=ks; k<=ke+1; k++) {
      for (j=js; j<=je+1; j++) {
        for (i=is; i<=ie+1; i++) {
          pGrid->B1i[k][j][i] = -(Ay[k+1][j][i]-Ay[k][j][i])/pGrid->dx3;
	  pGrid->B2i[k][j][i] = 0.0;
	  rsf = ri[i+1]/r[i];
	  lsf = ri[i]/r[i];
          pGrid->B3i[k][j][i] =  (rsf*Ay[k][j][i+1]-lsf*Ay[k][j][i])/pGrid->dx1;		
        }
      }
    }
  /* Sync poloidal centered fields to face fields for the next step */
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
	  rsf = ri[i+1]/r[i];
          lsf = ri[i]/r[i];
          pGrid->U[k][j][i].B1c = 0.5*(lsf*pGrid->B1i[k][j][i]+rsf*pGrid->B1i[k][j][i+1]);
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
	  rsf = ri[i+1]/r[i];
          lsf = ri[i]/r[i];
          pGrid->U[k][j][i].B1c = 0.5*(lsf*pGrid->B1i[k][j][i]+rsf*pGrid->B1i[k][j][i+1]);
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
#ifdef MHD
  Real divb;
/* Finally, let's check that the field is divergenceless */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	  rsf = ri[i+1]/r[i];
          lsf = ri[i]/r[i];
        divb = (rsf*pGrid->B1i[k][j][i+1]-lsf*pGrid->B1i[k][j][i])/pGrid->dx1 +
               (pGrid->B2i[k][j+1][i]-pGrid->B2i[k][j][i])/(pGrid->dx2*r[i]) +
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

  StaticGravPot = PseudoNewton;


  x1GravAcc = PseudoNewtonAccR;
  bvals_mhd_fun(pDomain,left_x1,disk_ir);
  bvals_mhd_fun(pDomain,right_x1,disk_or);
  bvals_mhd_fun(pDomain,left_x3,disk_iz);
  bvals_mhd_fun(pDomain,right_x3,disk_oz);

/* The above functions only enroll boundary conditions for the MHD part */
/* Boundary conditinos for specific intensity can be enrolled with bvals_fullrad_trans_fun */
	
#ifdef FULL_RADIATION_TRANSFER
	/* enroll the opacity function */
	get_full_opacity = Thindiskopacity;
#endif
	
	
/* enroll new history variables, only once  */

  if (frst == 1) {
    dump_history_enroll(hst_rho_VrVphi, "<rho Vr Vphi>");
    dump_history_enroll(hst_T,"<T>");
    dump_history_enroll(hst_T2,"<T^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif

#ifdef MHD
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BrBphi, "<-Br Bphi>");
    dump_history_enroll(hst_EB,"<B^2/2>");

#endif /* MHD */

#ifdef FULL_RADIATION_TRANSFER
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


 /* Set initial specific intensity  *
  * As we run MHD without radiation first, 
  * we set specific intensity to be zero *
  * proper values needed to be setup when restart the simulation */

	PrimS Wtemp;
        offset = nghost - Radghost;

	/* Now initialize the radiation quantities */

      for (k=ksr; k<=ker; k++) {
        for (j=jsr; j<=jer; j++) {
           for (i=isr; i<=ier; i++) {
               for(ifr=0; ifr<nf; ifr++){
                   for(l=0; l<noct; l++){
                       for(n=0; n<nang; n++){
                           if(ker > ksr)
                               cc_pos(pGrid,i-Radghost+nghost,j-Radghost+nghost,k-Radghost+nghost,&x1,&x2,&x3);
                           else
                               cc_pos(pGrid,i-Radghost+nghost,j-Radghost+nghost,0,&x1,&x2,&x3);

                           kg = k + offset;
                           jg = j + offset;
                           ig = i + offset;
                           Wtemp = Cons_to_Prim(&(pGrid->U[kg][jg][ig]));
                           temperature = Wtemp.P/(Wtemp.d * R_ideal);

                           pRG->imu[k][j][i][ifr][l*pRG->nang+n] = temperature * temperature * temperature * temperature/(4.0 * PI);
			
                       }/* n */
                   }/* l */
               }/* ifr */
           }/* i */
        }/* j */
    }/* k */

	for(k=ksr; k<=ker; k++){
		for(j=jsr; j<=jer; j++){
			for(i=isr; i<=ier; i++){
                for(ifr=0; ifr<nf; ifr++){

                    Thindiskopacity(pGrid, ifr, i-Radghost+nghost, j-Radghost+nghost, k, &(pRG->R[k][j][i][ifr].Sigma[0]));
                }/* ifr */
			}/* i */
		}/* j */
	}/* k */

    CalMoment(isr, ier, jsr, jer, ksr, ker, pRG);

   



	/* set boundary condition */
/*	bvals_fullrad_trans_fun(pDomain, left_x2, TwoBeam_ix2);
*/

	/* Enroll the user defined function to dump history of specific intensity */
/*	Intensity_history_enroll(0, 0, 0, "<I_l0_n0>");
*/



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

#ifdef FULL_RADIATION_TRANSFER
	fwrite(&Prat,sizeof(Real),1,fp);
	fwrite(&Crat,sizeof(Real),1,fp);
	fwrite(&R_ideal,sizeof(Real),1,fp);
 	fwrite(&kappaes,sizeof(Real),1,fp);
	fwrite(&kappaffP,sizeof(Real),1,fp);
	fwrite(&kappaffR,sizeof(Real),1,fp);
#endif	
	

#ifdef MHD
	fwrite(&betay,sizeof(Real),1,fp);
	fwrite(&pres,sizeof(Real),1,fp);
#endif


	


#ifdef RESISTIVITY
	fwrite(&eta_Ohm,sizeof(Real),1,fp);
        fwrite(&Jrhomax,sizeof(Real),1,fp);
#endif

   	fwrite(&zbtm,sizeof(Real),1,fp);
    	fwrite(&ztop,sizeof(Real),1,fp);


#ifdef FULL_RADIATION_TRANSFER
    	fwrite(&Tunit,sizeof(Real),1,fp);
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
	
#ifdef FULL_RADIATION_TRANSFER
	fread(&Prat,sizeof(Real),1,fp);
	fread(&Crat,sizeof(Real),1,fp);
	fread(&R_ideal,sizeof(Real),1,fp);
	fread(&kappaes,sizeof(Real),1,fp);
	fread(&kappaffP,sizeof(Real),1,fp);
	fread(&kappaffR,sizeof(Real),1,fp);
	
#endif
	
#ifdef MHD
	fread(&betay,sizeof(Real),1,fp);
	fread(&pres,sizeof(Real),1,fp);
#endif



#ifdef RESISTIVITY
	fread(&eta_Ohm,sizeof(Real),1,fp);
	fread(&Jrhomax,sizeof(Real),1,fp);
	
#endif

       fread(&zbtm,sizeof(Real),1,fp);
       fread(&ztop,sizeof(Real),1,fp);


#ifdef FULL_RADIATION_TRANSFER
       fread(&Tunit,sizeof(Real),1,fp);
#endif


#ifdef MPI_PARALLEL 
  if(myID_Comm_world == 0){
#endif
	/* Print out the parameters for consistency check */
     printf("Parameters: Prat %G Crat %G R_ideal %G kappaes %G kappaffP %G kappaffR %G\n",Prat,Crat,R_ideal,kappaes,kappaffP,kappaffR);
#ifdef MPI_PARALLEL 
  }
#endif


  /* Enroll necessary functions */




  StaticGravPot = PseudoNewton;


  x1GravAcc = PseudoNewtonAccR;
 /* Vacuume boundary condition for MHD part */


#ifdef FULL_RADIATION_TRANSFER
	/* enroll the opacity function */
 get_full_opacity = Thindiskopacity;
#endif
	

  StaticGravPot = PseudoNewton;
  x1GravAcc = PseudoNewtonAccR;
 

/* enroll new history variables */

	
    dump_history_enroll(hst_rho_VrVphi, "<rho Vx dVy>");
    dump_history_enroll(hst_T,"<T>");
    dump_history_enroll(hst_T2,"<T^2>");
#ifdef ADIABATIC
    dump_history_enroll(hst_E_total, "<E + rho Phi>");
#endif
	
#ifdef MHD
    dump_history_enroll(hst_Bx, "<Bx>");
    dump_history_enroll(hst_By, "<By>");
    dump_history_enroll(hst_Bz, "<Bz>");
    dump_history_enroll(hst_BrBphi, "<-Bx By>");
    dump_history_enroll(hst_EB,"<B^2/2>");
#endif

#ifdef FULL_RADIATION_TRANSFER
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
  RadGridS *pRG;

  Real x1, x2, x3;
  Real T, rho;
  PrimS Wtemp;
    Real T0, coef1, coef2, coef3;

 for(nl=0; nl<(pM->NLevels); nl++){
  for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
    if(pM->Domain[nl][nd].Grid != NULL){


	pD= &(pM->Domain[nl][nd]);
	pRG = (pD->RadGrid);
	pGrid = pD->Grid;

	/* Enroll boundary condition */

	bvals_mhd_fun(pD,left_x1,disk_ir);
  	bvals_mhd_fun(pD,right_x1,disk_or);
  	bvals_mhd_fun(pD,left_x3,disk_iz);
  	bvals_mhd_fun(pD,right_x3,disk_oz);
	


	int isr = pRG->is, ier = pRG->ie;
  	int jsr = pRG->js, jer = pRG->je;
  	int ksr = pRG->ks, ker = pRG->ke; 
	
	int i, j, k, ifr, l, n, m, offset;
	int ig, jg, kg;
	offset = nghost - Radghost;

	int nf = pRG->nf, nang = pRG->nang;
	int noct = pRG->noct;	

/*
    for (k=ksr; k<=ker; k++) {
       for (j=jsr; j<=jer; j++) {
         for (i=isr; i<=ier; i++) {
            for(ifr=0; ifr<nf; ifr++){
                for(l=0; l<noct; l++){
                    for(n=0; n<nang; n++){
            
                        ig = i + offset;
                        jg = j + offset;
                        kg = k + offset;

                        cc_pos(pGrid,ig,jg,kg,&x1,&x2,&x3);

                         Wtemp = Cons_to_Prim(&(pGrid->U[kg][jg][ig]));
                        
                        rho = Wtemp.d;
                        T0 = Wtemp.P/(rho * R_ideal);
                        if(T0 <Tfloor){
                          Wtemp.P = rho * R_ideal * Tfloor;
						  T = Tfloor;
                        }
                        else if(T0 < 3.0){
                            T = T0;
                        }
                        else{

                                          
                            coef1 = Prat/3.0;
                            coef2 = rho * R_ideal;
                            coef3 = -Wtemp.P;
                                          
                                        
                            T = rtsafe(Tequilibrium, 0.0, T0, 1.e-12, coef1, coef2, coef3,0.0);
                            if(T < 1.0)
			    	T = 1.0;
			 
			    Wtemp.P = T * rho * R_ideal;
                        
                        }

                         pGrid->U[kg][jg][ig] = Prim_to_Cons(&(Wtemp));
                         pRG->imu[k][j][i][ifr][l*nang+n] = T * T * T * T /(4.0 * PI);

                    }
                }
              }
            }
          }
        }

	   for(k=ksr; k<=ker; k++){
	      for(j=jsr; j<=jer; j++){
            for(i=isr; i<=ier; i++){
                for(ifr=0; ifr<nf; ifr++){
                    Thindiskopacity(pGrid, ifr, i+offset, j+offset, k+offset, &(pRG->R[k][j][i][ifr].Sigma[0]));
                }
            }
	      }
	   }
 
		CalMoment(isr, ier, jsr, jer, ksr, ker, pRG);

*/



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
	const Real *r=pG->r;
	const Real *ri=pG->ri;
	Real dx2i, dx1i, dx3i;
	Real rsf, lsf;
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
		dx2i = 1.0/(r[i] * pG->dx2);
		dx3i = 1.0/(pG->dx3);
		dx1i = 1.0/(pG->dx1);
		rsf = r[i]/ri[i];
		lsf = r[i-1]/ri[i];

		jx1 = dx2i * (pG->B3i[k][j][i] - pG->B3i[k  ][j-1][i  ]) -
                       dx3i * (pG->B2i[k][j][i] - pG->B2i[k-1][j  ][i  ]);

		jx2 = dx2i * (pG->B3i[k][j+1][i] - pG->B3i[k  ][j][i  ]) -
                       dx3i *  (pG->B2i[k][j+1][i] - pG->B2i[k-1][j+1 ][i  ]);

		jx3 = (pG->B3i[k+1][j][i] - pG->B3i[k+1][j-1][i  ])*dx2i -
                        (pG->B2i[k+1][j][i] - pG->B2i[k][j  ][i  ])*dx3i;

		jx4 = (pG->B3i[k+1][j+1][i] - pG->B3i[k+1][j][i  ])*dx2i -
                        (pG->B2i[k+1][j+1][i] - pG->B2i[k][j+1][i  ])dx3i;

		jx = 0.25 * (jx1 + jx2 + jx3 + jx4);

		/* y component */
		jy1 = (pG->B1i[k][j][i] - pG->B1i[k-1][j  ][i  ])*dx3i -
                        (pG->B3i[k][j][i] - pG->B3i[k  ][j  ][i-1])*dx1i;

		jy2 = (pG->B1i[k][j][i+1] - pG->B1i[k-1][j  ][i+1])*dx3i -
                        (pG->B3i[k][j][i+1] - pG->B3i[k  ][j  ][i])*dx1i;

		jy3 = (pG->B1i[k+1][j][i] - pG->B1i[k][j  ][i  ])*dx3i -
                        (pG->B3i[k+1][j][i] - pG->B3i[k+1][j  ][i-1])*dx1i;

		jy4 = (pG->B1i[k+1][j][i+1] - pG->B1i[k][j  ][i+1])*dx3i -
                        (pG->B3i[k+1][j][i+1] - pG->B3i[k+1][j ][i])*dx1i;

		jy = 0.25 * (jy1 + jy2 + jy3 + jy4);
		/* z component */
		rsf = r[i]/ri[i];
                lsf = r[i-1]/ri[i];

		jz1 = (rsf * pG->B2i[k][j][i] - lsf * pG->B2i[k  ][j  ][i-1])*dx1i -
                        (pG->B1i[k][j][i] - pG->B1i[k  ][j-1][i  ])*dx2i;

		rsf = r[i+1]/ri[i+1];
                lsf = r[i]/ri[i+1];

		jz2 = (rsf * pG->B2i[k][j][i+1] -lsf * pG->B2i[k  ][j  ][i])*dx1i -
                        (pG->B1i[k][j][i+1] - pG->B1i[k  ][j-1][i+1  ])*dx2i;

		rsf = r[i]/ri[i];
                lsf = r[i-1]/ri[i];

		jz3 = (rsf * pG->B2i[k][j+1][i] -lsf * pG->B2i[k  ][j+1][i-1])*dx1i -
                        (pG->B1i[k][j+1][i] - pG->B1i[k  ][j][i  ])*dx2i;

		rsf = r[i+1]/ri[i+1];
                lsf = r[i]/ri[i+1];

		jz4 = (rsf * pG->B2i[k][j+1][i+1] -lsf * pG->B2i[k  ][j+1][i])*dx1i -
                        (pG->B1i[k][j+1][i+1] - pG->B1i[k  ][j][i+1])*dx2i;

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
	
	
	if((fabs(x3) > 30.0) || (x1 < 50.0)){
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

	if(x1 < 4 )
		*eta_O = eta0;

	*eta_H = 0.0;
	*eta_A = 0.0;

	return;
}
#endif



/* Do not allow inflow from inner boundary */
/* keep temperature and density to be the same */

void disk_ir(GridS *pGrid) {
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
	Real Vk,R,p,z, rho, Lper;

	PrimS Wtemp;
	
	Real rsf = 1.0;
	Real lsf = 1.0;

#ifdef CYLINDRICAL
  	const Real *r = pGrid->r, *ri=pGrid->ri;
#endif
	
#ifdef MHD
	/* B1i is not set at i=is-nghost */
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=1; i<=nghost-1; i++) {				
				pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
				
			}
		}
	}
	
	if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=ju; j++) {
			for (i=1; i<=nghost; i++) {
				if(pGrid->U[k][j][is].M1 < 0.0){				
					pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
				}
				else {
					pGrid->B2i[k][j][is-i] = 0.0 * pGrid->B2i[k][j][is];
				}
			}
		}
	}
	
	if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
	for (k=ks; k<=ku; k++) {
		for (j=js; j<=je; j++) {
			for (i=1; i<=nghost; i++) {
				if(pGrid->U[k][j][is].M1 < 0.0){
					pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
				}
				else {
					pGrid->B3i[k][j][is-i] = 0.0;
				}

				
				
			}
		}
	}
	/* Now set the cell centered value */
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=1; i<=nghost; i++) {
#ifdef CYLINDRICAL
				rsf = ri[is-i+1]/r[is-i];
				lsf = ri[is-i]/r[is-i];
	
#endif
				
				pGrid->U[k][j][is-i].B1c = 0.5 * (lsf * pGrid->B1i[k][j][is-i] + rsf * pGrid->B1i[k][j][is]);
				pGrid->U[k][j][is-i].B2c = 0.5 * (pGrid->B2i[k][j][is-i] + pGrid->B2i[k][j+1][is-i]);
				pGrid->U[k][j][is-i].B3c = 0.5 * (pGrid->B3i[k][j][is-i] + pGrid->B3i[k+1][j][is-i]);
				
			}
		}
	}
#endif /* MHD */
	

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	cc_pos(pGrid,is,j,k,&R,&p,&z);
		Wtemp = Cons_to_Prim(&(pGrid->U[k][j][is]));
		rho = Wtemp.d;
		
    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = R*pGrid->U[k][j][is].M2/rho;
#else
    	Vk = Vkep(R,p,z);
	/* residual angular momentum with respect to the Keplerian angular momentum */
    	Lper = R*pGrid->U[k][j][is].M2/rho - R*Vk;
#endif

		if(Wtemp.V1 > 0.0){
			Wtemp.V1 = 0.0;			
		}
		
      for (i=1; i<=nghost; i++) {

				
		  /* Calculate Keplerian velocity */
		  cc_pos(pGrid,is-i,j,k,&R,&p,&z);
		  Vk = Vkep(R,p,z);
	
		  Wtemp.V2 = Lper/R + Vk;

		  
		  pGrid->U[k][j][is-i].d = rho;
		  pGrid->U[k][j][is-i].M1 = Wtemp.V1 * rho;
		  pGrid->U[k][j][is-i].M2 = Wtemp.V2 * rho;
		  pGrid->U[k][j][is-i].M3 = Wtemp.V3 * rho;
		  
		  pGrid->U[k][j][is-i].E = Wtemp.P/Gamma_1 + 0.5 * rho * (Wtemp.V1 * Wtemp.V1 + Wtemp.V2 * Wtemp.V2 + Wtemp.V3 * Wtemp.V3);
		  
#ifdef MHD
		  pGrid->U[k][j][is-i].E += 0.5 * (SQR(pGrid->U[k][j][is-i].B1c) + SQR(pGrid->U[k][j][is-i].B2c) + SQR(pGrid->U[k][j][is-i].B3c));	
#endif
		  
		  /* No inflow from inner boundary */
		 


      }/* end i */
    }/* end j */
  }/* end k */


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
	Real Vk, R,p,z, rho, Lper;
	PrimS Wtemp;
	
	Real rsf = 1.0;
	Real lsf = 1.0;

#ifdef CYLINDRICAL
        const Real *r = pGrid->r, *ri=pGrid->ri;
#endif	
	
#ifdef MHD
	/* i=ie+1 is not a boundary condition for the interface field B1i */
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=2; i<=nghost; i++) {
				pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie+1];
				
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
	
	/* Now set the cell centered value */
	for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=1; i<=nghost; i++) {
				if(i<nghost){
				
#ifdef CYLINDRICAL
					rsf = ri[ie+i+1]/r[ie+i];
					lsf = ri[ie+i]/r[ie+i];
				
#endif
					
					pGrid->U[k][j][ie+i].B1c = 0.5 * (lsf * pGrid->B1i[k][j][ie+i] + rsf * pGrid->B1i[k][j][ie+i+1]);
				}
				else {
					pGrid->U[k][j][ie+i].B1c = pGrid->B1i[k][j][ie+i];
				}

				
				pGrid->U[k][j][ie+i].B2c = 0.5 * (pGrid->B2i[k][j][ie+i] + pGrid->B2i[k][j+1][ie+i]);
				pGrid->U[k][j][ie+i].B3c = 0.5 * (pGrid->B3i[k][j][ie+i] + pGrid->B3i[k+1][j][ie+i]);
				
			}
		}
	}
	
#endif /* MHD */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
		Wtemp = Cons_to_Prim(&(pGrid->U[k][j][ie]));
		rho = Wtemp.d;
		
		cc_pos(pGrid,ie,j,k,&R,&p,&z);

    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = R * Wtemp.V2;
#else
	
		Vk = Vkep(R,p,z);
		Lper = R*Wtemp.V2 - R*Vk;
#endif
      for (i=1; i<=nghost; i++) {
        				
		  cc_pos(pGrid,ie+i,j,k,&R,&p,&z);
		  Vk = Vkep(R,p,z);

		/*
 		  Wtemp.V2 = Lper/R + Vk;
		*/
		 if(Wtemp.V1 < 0.0)
			Wtemp.V1 = 0.0;
		  
		  pGrid->U[k][j][ie+i].d = rho;
		  pGrid->U[k][j][ie+i].M1 = rho * Wtemp.V1;
		  pGrid->U[k][j][ie+i].M2 = rho * Wtemp.V2;
		  pGrid->U[k][j][ie+i].M3 = rho * Wtemp.V3;
		  
		  pGrid->U[k][j][ie+i].E = Wtemp.P/Gamma_1 + 0.5 * rho * (Wtemp.V1 * Wtemp.V1 + Wtemp.V2 * Wtemp.V2 + Wtemp.V3 * Wtemp.V3);
		  
#ifdef MHD
		  pGrid->U[k][j][ie+i].E += 0.5 * (SQR(pGrid->U[k][j][ie+i].B1c) + SQR(pGrid->U[k][j][ie+i].B2c) + SQR(pGrid->U[k][j][ie+i].B3c));	
#endif
		  

      }
    }
  }


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

	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;

	static Real x3t;
	x3t = ztop - 0.5 * pGrid->dx3;
	PrimS Wopacity;


#ifdef MHD
	
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
#ifdef MHD
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
#ifdef MHD

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



#ifdef MHD
		pGrid->U[ke+k][j][i].E += 0.5 * (pGrid->U[ke+k][j][i].B1c * pGrid->U[ke+k][j][i].B1c + pGrid->U[ke+k][j][i].B2c * pGrid->U[ke+k][j][i].B2c + pGrid->U[ke+k][j][i].B3c * pGrid->U[ke+k][j][i].B3c);

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

	Real u0;
	u0 = 0.0;
	dx1 = pGrid->dx1;
	dx2 = pGrid->dx2;
	dx3 = pGrid->dx3;
	PrimS Wopacity;

	
	 static Real x3b;

  	x3b = zbtm+0.5*pGrid->dx3;



#ifdef MHD
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
#ifdef MHD
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
#ifdef MHD

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



#ifdef MHD
		pGrid->U[ks-k][j][i].E += 0.5 * (pGrid->U[ks-k][j][i].B1c * pGrid->U[ks-k][j][i].B1c + pGrid->U[ks-k][j][i].B2c * pGrid->U[ks-k][j][i].B2c + pGrid->U[ks-k][j][i].B3c * pGrid->U[ks-k][j][i].B3c);

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
	RadGridS *pRG;
        int i, j, k;
        int ie, is;
        int je, js;
        int ke, ks;
	int isr, ier, jsr, jer, ksr, ker;
	int l, n, m, ifr, nf, nang, noct;

	int badcellflag;
        Real pressure, velocity, velocity_x, velocity_y, velocity_z, Bpre, beta;
	Real temprho, tempT, T0, temperature;
	PrimS Wtemp;
	Real Vmax = 0.8 * Crat;
	
        
	Real x1, x2, x3, distance, distance1;
#ifdef MHD
	Real Vmag;
	/* Limit the Alfven velocity in the mask region so that the time step will not be limited */
#endif


	int nl, nd;


	for(nl=0; nl<(pM->NLevels); nl++){
	for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
	if(pM->Domain[nl][nd].Grid != NULL){



	pG = pM->Domain[nl][nd].Grid;
	pRG= pM->Domain[nl][nd].RadGrid;

        ie = pG->ie;
        is = pG->is;
        je = pG->je;
        js = pG->js;
        ke = pG->ke;
        ks = pG->ks;

	isr = pRG->is; ier = pRG->ie;
	jsr = pRG->js; jer = pRG->je;
	ksr = pRG->ks; ker = pRG->ke;
	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;




        for(k=ks; k<=ke; k++) {
                for (j=js; j<=je; j++) {
                        for (i=is; i<=ie; i++) {



   			badcellflag = 0;

			velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
                         velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                         velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;

			Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));

			/* Limit the maximum velocity */
			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);				
			if(velocity > Vmax){
				velocity_x *= Vmax/velocity;
				velocity_y *= Vmax/velocity;
				velocity_z *= Vmax/velocity;
				
				Wtemp.V1 = velocity_x;
				Wtemp.V2 = velocity_y;
				Wtemp.V3 = velocity_z; 

				badcellflag = 1;
			}
			
		
			if(Wtemp.d < dfloor){
				Wtemp.d = dfloor;
				badcellflag = 1;
			}

			temperature = Wtemp.P/(Wtemp.d*R_ideal);

			if(temperature > 0.5 * Crat * Crat){
				temperature =  0.5 * Crat * Crat;
				Wtemp.P = Wtemp.d * temperature * R_ideal;

                                badcellflag = 1;
			}

			if(temperature < 0.99*Tfloor){
/*
#ifndef MPI_PARALLEL
				printf("Below Tfloor i: %d j: %d k: %d T: %e rho: %e\n",i,j,k,temperature,Wtemp.d);
#else
				printf("Below Tfloor i: %d j: %d k: %d T: %e rho: %e ID: %d\n",i,j,k,temperature,Wtemp.d,myID_Comm_world);
#endif
*/
				Wtemp.P = Wtemp.d * Tfloor * R_ideal;

				badcellflag = 1;
			}	


			/* Reset the temperature */
			/* This is only used for isothermal runs */

/*
			 temperature = pG->tgas[k][j][i];
                        Wtemp.P = Wtemp.d * temperature * R_ideal;
                        badcellflag = 1;

*/


			if(badcellflag){
				pG->U[k][j][i] = Prim_to_Cons(&(Wtemp));	
			}



                        } /* i */
                }/* j */
        }/* k */


	/* Check specific intensity to make sure it is non-negative */

    for (k=ksr; k<=ker; k++) {
      	for (j=jsr; j<=jer; j++) {
            for (i=isr; i<=ier; i++) {
                for(ifr=0; ifr<nf; ifr++){
                    for(l=0; l<noct; l++){
                        for(n=0; n<nang; n++){
                            
                            if(pRG->imu[k][j][i][ifr][l*nang+n] < TINY_NUMBER)
                                pRG->imu[k][j][i][ifr][l*nang+n] = TINY_NUMBER;
		
                        }/* end n */
                    }/* end l */
                }/* end ifr */
		  }/* end i */
		 }/* end j */
		}/* end k */
	       
	}/* Grid != NULL */
	}/* nd */
	} /* nl */
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



/* Paczynski - Witta potential */
static Real PseudoNewton(const Real x1, const Real x2, const Real x3)
{
	Real GM, distance, potential, h;
	/* The dimensionless number */
	GM = Crat*Crat*0.5;
	if(x3 > (ztop - 0.5 * dz)){
		h = ztop - 0.5 * dz;
	}
	else if(x3 < (zbtm + 0.5 * dz)){
		h = zbtm + 0.5 * dz;
	}
	else{
		h = x3;
	}

	
	distance = sqrt(x1*x1+h*h);
	potential = -GM/(distance - 1.0);
	
	
	return potential;
}


/* Paczynski - Witta potential */
static Real PseudoNewtonAccR(const Real x1, const Real x2, const Real x3)
{
	Real GM, acc, distance, h;
	/* The dimensionless number */
	GM = Crat*Crat*0.5;
	if(x3 > (ztop - 0.5 * dz)){
		h = ztop - 0.5 * dz;
	}
	else if(x3 < (zbtm + 0.5 * dz)){
		h = zbtm + 0.5 * dz;
	}
	else{
		h = x3;
	}

	distance = sqrt(x1*x1+h*h);
	
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




static void Thindiskopacity(GridS *pG, const int ifr, const int i, const int j, const int k, Real *Sigma)
{
/* Sigma[0-NOPACITY] are: Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE respectively */
/* dSigma[0-2*NOPACITY] are: dSigma_sF/drho, dSigma_aF/drho, dSigma_aP/drho, dSigma_aE/drho */
/* 			     dSigma_sF/dT,   dSigma_aF/dT,   dSigma_aP/dT,   dSigma_aE/dT */
	
/* When pressure becomes negative, we do not include radiation source term */
	Real T, rho;


	rho = pG->U[k][j][i].d;
	T = pG->U[k][j][i].E - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) 
        	+ SQR(pG->U[k][j][i].M3))/rho;
#ifdef MHD
	T -= 0.5*(SQR(pG->U[k][j][i].B1c)
        	+ SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
	T *= Gamma_1;
	T /= (rho * R_ideal); 

	if(T < Tfloor)
		T = Tfloor;


if((T > TINY_NUMBER) && (rho > TINY_NUMBER)){	
	Real Tpower, Tpower1;
	
	/* Tpower = T^3.5 , Tpower1 = T^4.5 */
	 /* For FULL_Radiation_Transfer *
          * Sigma[0] and Sigma[1] are absorption opacity *
          * Sigma[2] and Sigma[3] are scattering opacity */

	Tpower = 1.0 / (T * T * T * sqrt(T));
	Tpower1 = Tpower / T;
	
	if(Sigma != NULL){
		Sigma[0] =  kappaffP * rho * rho * Tpower;
		Sigma[1] =  kappaffP * rho * rho * Tpower;
		Sigma[2] =  kappaes * rho;
                Sigma[3] =  kappaes * rho;
	}	
	
}
else{

	if(Sigma != NULL){
		Sigma[0] =  0.0;
		Sigma[1] =  0.0;
		Sigma[2] =  0.0;
		Sigma[3] =  0.0;
	}	
	
}
	
	return; 
	
}





/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_beta(const GridS *pG, const int i, const int j, 
 *			      const int k)
 *  \brief Computes beta=P/(B^2/2)  
 */
static Real expr_beta(GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  Real pre;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
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
static Real expr_ME(GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,B2;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);
#ifdef MHD
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
static Real expr_KE(GridS *pG, const int i, const int j, const int k)
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
  Real vr, vphi;

  vr = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  vphi = pG->U[k][j][i].M2/pG->U[k][j][i].d;
 
  return (pG->U[k][j][i].d * vr * vphi);

}



static Real hst_T(const GridS *pG, const int i, const int j, const int k)
{

	Real pressure;
	pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef MHD
	pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	pressure *= Gamma - 1.0;

	return pressure/(R_ideal * pG->U[k][j][i].d);
}

static Real hst_P(const GridS *pG, const int i, const int j, const int k)
{

	Real pressure;
	pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef MHD
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




static Real hst_sigmas(const GridS *pG, const int i, const int j, const int k)
{
	Real Sigma[4];
	Thindiskopacity( (GridS *)pG, 0, i, j,  k, &(Sigma[0]));
	return Sigma[0];
}

static Real hst_sigmaaP(const GridS *pG, const int i, const int j, const int k)
{
	Real Sigma[4];
	Thindiskopacity( (GridS *)pG, 0, i, j,  k, &(Sigma[0]));
	return Sigma[2];
}



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

#ifdef MHD
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

  return (-pG->U[k][j][i].B1c*pG->U[k][j][i].B2c);
}


static Real hst_EB(const GridS *pG, const int i, const int j, const int k)
{

	return 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
}


#endif /* MHD */




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
/* The above is for cartesian coordinate */
/* For cylindrical coordinate, x2 is already the azimuthal component */

/* This function calculate the radial profile, azimuthally and vertically averaged */
/* For each fixed r, the volume for each cell is the same */
/* So, volume average is the same as averaging the values */

static void output_1d(MeshS *pM, OutputS *pOut)
{
	int nl, nd;
	
	FILE *p_1dfile;
	char *fname, *plev=NULL, *pdom=NULL;
	char levstr[8], domstr[8];
	int big_end = ath_big_endian();
	
	int i,j,k,ig,n;
	int tot1d,i1d,nrmx, idispG, idispD;
	int dnum = pOut->num;
	
	double **out1d;
	double *out_r;
	double x1,x2,x3,press, darea;
	double distance;
	double vx, vy, vz, Fr01, Fr02, Fr03;
#ifdef MHD
	double Br, Bphi;
 	double PB1, PB2;
#endif
	double Pg1, Pg2;
	
#ifdef MPI_PARALLEL
	double *my_out1d;
	double *g_out1d;	
	int zproc;
	int ierr,myID_Comm_Domain;
	int nproc;
#endif
	
	/* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap, fxx, fyy, fzz, fxy, fxz, fyz */
	
#ifdef MHD
	tot1d=51;
#else
	tot1d=41;
#endif /* MHD */
	
	
	
	
	int is, ie, js, je, ks, ke;
	/* Loop over all Domain and level */
	int offset, nf, ifr;
	int ir, jr, kr;
	offset = Radghost - nghost;
	
	/* The frequency averaged quantities */
	Real Fr[3], Pr[6], Edd[6], Sigma[4];
	Real Er, Erk1, Erk0, Erj1, Erj0, Eri0, Eri1;
	
	DomainS *pD;
	GridS *pGrid;
	RadGridS *pRG;
	
	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				
				if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
					(pOut->ndomain == -1 || pOut->ndomain == nd)){
					
					pD = &(pM->Domain[nl][nd]);
					pGrid = pD->Grid;
					pRG = pD->RadGrid;
					
					nf = pRG->nf;
										
					
					
					/* total number of radial grid zones for this Domain only */			
					nrmx = pD->Nx[0]; 
					darea = 1.0/(pD->Nx[1] * pD->Nx[2]);
					
					
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
					/* This is for static mesh refinement */
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						out_r = (double *) calloc_1d_array(nrmx,sizeof(double));
#ifdef MPI_PARALLEL
					}/* End myID */
#endif
					
										
					out1d = (double **) calloc_2d_array(nrmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
					my_out1d = (double *) calloc_1d_array(nrmx,sizeof(double));
					g_out1d = (double *) calloc_1d_array(nrmx,sizeof(double));
					
#endif
					for (i=0; i<nrmx; i++) {						
						for (i1d=0; i1d<tot1d; i1d++) {
							out1d[i][i1d] = 0.0;
						}
					}
					
					
					idispD = pD->Disp[0];
					idispG = pGrid->Disp[0];
					
					
					/* First calculate the x1 coordinate and save it to be dumped
					 by root in every 1d file */
					
					/* There will be some grids that not covered in the radial grids */
					/* This needs to be done for every CPU because we need this array to determine the position */
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif					
						for (i=0; i<nrmx; i++) {
							/* Vertical coordinate for this whole domain */
							x1 = pD->MinX[0] + (i + 0.5)*pGrid->dx1;
							out_r[i] = x1;
						}
#ifdef MPI_PARALLEL
					}/* End my ID */
#endif						
					
					/* Compute 1d averaged variables */
					/* The position of this grid in the radial array is determined by comparing the distance */
					
					for (k=ks; k<=ke; k++) {    				
						for (j=js; j<=je; j++) {
							for (i=is; i<=ie; i++) {
								/* First, calculate the transformation parameters */
								/* For cylindrical coordinate, x1 is radius and x3 is vertical height */
								ig=i+idispG-idispD-nghost;
								
								ir = i + offset;
								jr = j + offset;
								kr = k + offset;

								/* Average the radiation quantities over frequencies */
								
								for(n=0; n<3; n++)
									Fr[n] = 0.0;
								for(n=0; n<4; n++)
									Sigma[n] = 0.0;
								for(n=0; n<6; n++)
									Pr[n] = 0.0;
								
								Er = 0.0;
								Erk1 = 0.0;
								Erk0 = 0.0;
								Erj1 = 0.0;
								Erj0 = 0.0;
								Eri1 = 0.0;
								Eri0 = 0.0;
								
								for(ifr=0; ifr<nf; ifr++){
									for(n=0; n<3; n++)
										Fr[n] += 4.0 * PI * pRG->R[kr][jr][ir][ifr].H[n] * pRG->wnu[ifr];
									for(n=0; n<4; n++)
										Sigma[n] += pRG->R[kr][jr][ir][ifr].Sigma[n] * pRG->wnu[ifr];
									for(n=0; n<6; n++)
										Pr[n] += 4.0 * PI * pRG->R[kr][jr][ir][ifr].K[n] * pRG->wnu[ifr];
									
									Er += 4.0 * PI * pRG->R[kr][jr][ir][ifr].J * pRG->wnu[ifr];
									Erk1 += 4.0 * PI * pRG->R[kr+1][jr][ir][ifr].J * pRG->wnu[ifr];
									Erk0 += 4.0 * PI * pRG->R[kr-1][jr][ir][ifr].J * pRG->wnu[ifr];
									Erj1 += 4.0 * PI * pRG->R[kr][jr+1][ir][ifr].J * pRG->wnu[ifr];
									Erj0 += 4.0 * PI * pRG->R[kr][jr-1][ir][ifr].J * pRG->wnu[ifr];
									Eri1 += 4.0 * PI * pRG->R[kr][jr][ir+1][ifr].J * pRG->wnu[ifr];
									Eri0 += 4.0 * PI * pRG->R[kr][jr][ir-1][ifr].J * pRG->wnu[ifr];
								}
								
								for(n=0; n<6; n++)
									Edd[n] = Pr[n] / Er;
								
								
								cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
								distance = x1;								
								
								
									
								vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
								vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
								vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;

								/* density */
								i1d=0;
								out1d[ig][i1d] += pGrid->U[k][j][i].d;
								i1d++;
#ifdef ISOTHERMAL
								out1d[ig][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
								press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#ifdef MHD
																   - expr_ME(pGrid,i,j,k)
#endif
																   ),TINY_NUMBER);

								Pg1           = MAX(Gamma_1*(pGrid->U[k][j][i-1].E - expr_KE(pGrid,i-1,j,k)
#ifdef MHD
																   - expr_ME(pGrid,i-1,j,k)
#endif
																   ),TINY_NUMBER);

								Pg2           = MAX(Gamma_1*(pGrid->U[k][j][i+1].E - expr_KE(pGrid,i+1,j,k)
#ifdef MHD
																   - expr_ME(pGrid,i+1,j,k)
#endif
																   ),TINY_NUMBER);



									/* pressure */
								out1d[ig][i1d] += press;
#endif



#ifdef ADIABATIC
								i1d++;
									/* temperature */
								out1d[ig][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
								i1d++;
									/* Total E */
								out1d[ig][i1d] += pGrid->U[k][j][i].E;
								i1d++;
								out1d[ig][i1d] += hst_E_total(pGrid,i,j,k);
#endif
									/* momentum */
								i1d++;
								out1d[ig][i1d] += pGrid->U[k][j][i].M1;
								i1d++;
								out1d[ig][i1d] += pGrid->U[k][j][i].M2;
								i1d++;
								out1d[ig][i1d] += pGrid->U[k][j][i].M3;
								/* kinetic energy */
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
								i1d++;
								out1d[ig][i1d] += expr_KE(pGrid,i,j,k);
								i1d++;
								out1d[ig][i1d] += hst_rho_VrVphi(pGrid,i,j,k);
								/* vr */
								i1d++;
								out1d[ig][i1d] += vx;
								/* vphi */
								i1d++;
								out1d[ig][i1d] += vy;
								/* vz */
								i1d++;
								out1d[ig][i1d] += vz;
#ifdef MHD
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
								i1d++;
								out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
								i1d++;
								out1d[ig][i1d] += expr_ME(pGrid,i,j,k);
								/* Alfven velocity */
								i1d++;
								out1d[ig][i1d] += sqrt(expr_ME(pGrid,i,j,k)/pGrid->U[k][j][i].d);
								/* Alfven velocity for Bz */
								i1d++;
								out1d[ig][i1d] += sqrt(0.5*SQR(pGrid->U[k][j][i].B3c)/pGrid->U[k][j][i].d);
								i1d++;
								out1d[ig][i1d] += hst_Bx(pGrid,i,j,k);
								i1d++;
								out1d[ig][i1d] += hst_By(pGrid,i,j,k);
								i1d++;
								out1d[ig][i1d] += hst_Bz(pGrid,i,j,k);
								i1d++;
								out1d[ig][i1d] += hst_BrBphi(pGrid,i,j,k);

								PB1 = expr_ME(pGrid,i-1,j,k);
								PB2 = expr_ME(pGrid,i+1,j,k);
#endif
								/* rho v mass flux */
								i1d++;
								out1d[ig][i1d] += pGrid->U[k][j][i].M1 * x1;
									
									/* angular momentum */
								i1d++;
								out1d[ig][i1d] += pGrid->U[k][j][i].M2 * distance;
								/* specific angular momentum */
								i1d++;
								out1d[ig][i1d] += vy * distance;
								i1d++;
								out1d[ig][i1d] += Er;
								i1d++;
								out1d[ig][i1d] += Fr[0];
								i1d++;
								out1d[ig][i1d] += Fr[1];
								i1d++;
								out1d[ig][i1d] += Fr[2];
									/* co-moving flux */
								Fr01 = Fr[0] -(vx * (1.0 + Edd[0]) + vy * Edd[1] + vz * Edd[3]) * Er/Crat;
								Fr02 = Fr[1] -(vy * (1.0 + Edd[2]) + vx * Edd[1] + vz * Edd[4]) * Er/Crat;
								Fr03 = Fr[2] -(vz * (1.0 + Edd[5]) + vx * Edd[3] + vy * Edd[4]) * Er/Crat;								
								i1d++;
								/* Fr0r */
								out1d[ig][i1d] += Fr01;
								/* Fr0phi */
								i1d++;
								out1d[ig][i1d] += Fr02;

								/* Fr0z */
								i1d++;
								out1d[ig][i1d] += Fr03;
								
								/* advection Er */
								i1d++;
								out1d[ig][i1d] += Er * vx;
								/* kappes and kappaff */
									
								i1d++;
								out1d[ig][i1d] += hst_sigmas(pGrid,i,j,k);
								i1d++;
								out1d[ig][i1d] += hst_sigmaaP(pGrid,i,j,k);
								/* radiation work term */
								i1d++;
								out1d[ig][i1d] += (Sigma[1] - Sigma[2]) * (Fr01 * vx + Fr02 * vy + Fr03 * vz); 
								/* dErdx */
								i1d++;
								out1d[ig][i1d] += (Eri1 - Eri0)/(2.0 * pGrid->dx1);
								i1d++;
								out1d[ig][i1d] += (Erj1 - Erj0)/(2.0 * pGrid->dx2);
								i1d++;
								out1d[ig][i1d] += (Erk1 - Erk0)/(2.0 * pGrid->dx3);
								/* Eddington tensor */
								i1d++;
								out1d[ig][i1d] += Edd[0];
								i1d++;
								out1d[ig][i1d] += Edd[2];
								i1d++;
								out1d[ig][i1d] += Edd[5];
								i1d++;
								out1d[ig][i1d] += Edd[1];
								i1d++;
								out1d[ig][i1d] += Edd[3];
								i1d++;
								out1d[ig][i1d] += Edd[4];
								
								/* Acceleration due to gas pressure gradient */
								i1d++;
								out1d[ig][i1d] += (Pg2-Pg1)/(2.0 * pGrid->dx1 * pGrid->U[k][j][i].d);
#ifdef MHD
								/* Acceleration due to magnetic pressure gradient */
								i1d++;
								out1d[ig][i1d] += (PB2-PB1)/(2.0 * pGrid->dx1 * pGrid->U[k][j][i].d);
#endif
									

								
							}/* end i */
						}/* end j */
					}/* end k*/
					
					
					
					/* The parent sums the scal[] array.
					 * Note that this assumes (dx1,dx2,dx3) = const. */
					
#ifdef MPI_PARALLEL 
					for(i1d=0; i1d<tot1d; i1d++){
						for (i=0; i<nrmx; i++) {
							my_out1d[i] = out1d[i][i1d];
						}
						ierr = MPI_Reduce(my_out1d, g_out1d, nrmx,
										  MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
						if(ierr)
							ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
						for (i=0; i<nrmx; i++) {
							out1d[i][i1d] = g_out1d[i];
						}
					}
					
					
#endif
					
					
					/* For parallel calculations, only the parent computes the average
					 * and writes the output. */
#ifdef MPI_PARALLEL
					if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif
						
						/* When average over the azimuthal direction, */
						/* The area for each cell is the same */
						for (i=0; i<nrmx; i++) {
							for (i1d=0; i1d<tot1d; i1d++) {
								out1d[i][i1d] *= darea;
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
							
#ifdef MHD
							if (k == 0) {
								fprintf(p_1dfile,"# [1]r	[2]dens	[3]pressure	[4]temperature	[5]E	[6]Etot	[7]Mr	[8]Mphi	[9]Mz	[10]KEr	[11]KEphi	[12]KEz	[13]KE	[14]rhoVrVphi	[15]Vr	[16]Vphi	[17]Vz	[18]MEx	[19]MEy	[20]MEz	[21]ME	[22]Va	[23]Vaz	[24]Bx	[25]By	[26]Bz	[27]BrBphi	[28]rhoVrR	[29]rhoVphiR	[30]VphiR	[31]Er	[32]Frr	[33]Fphi	[34]Fz	[35]Fr0r	[36]Fr0phi	[37]Fr0z	[38]VrEr	[39]kappaes	[40]kappaff	[41]Frwork	[42]dErdx	[43]dErdy	[44]dErdz	[45]f_xx	[46]f_yy	[47]f_zz	[48]f_xy	[49]f_xz	[50]f_yz	[51]dPgdr/rho	[52]dPB/dr/rho\n");
							}
							fprintf(p_1dfile,"%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G\n",out_r[k],out1d[k][0],out1d[k][1],out1d[k][2],
									out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
									out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28],out1d[k][29],out1d[k][30],out1d[k][31],out1d[k][32],out1d[k][33],out1d[k][34],out1d[k][35],out1d[k][36],out1d[k][37],out1d[k][38],out1d[k][39],out1d[k][40],out1d[k][41],out1d[k][42],out1d[k][43],out1d[k][44],out1d[k][45],out1d[k][46],out1d[k][47],out1d[k][48],out1d[k][49],out1d[k][50]);
#else

							if (k == 0) {
								fprintf(p_1dfile,"# [1]r	[2]dens	[3]pressure	[4]temperature	[5]E	[6]Etot	[7]Mr	[8]Mphi	[9]Mz	[10]KEr	[11]KEphi	[12]KEz	[13]KE	[14]rhoVrVphi	[15]Vr	[16]Vphi	[17]Vz		[18]rhoVrR	[19]rhoVphiR	[20]VphiR	[21]Er	[22]Frr	[23]Fphi	[24]Fz	[25]Fr0r	[26]Fr0phi	[27]Fr0z	[28]VrEr	[29]kappaes	[30]kappaff	[31]Frwork	[32]dErdx	[33]dErdy	[34]dErdz	[35]f_xx	[36]f_yy	[37]f_zz	[38]f_xy	[39]f_xz	[40]f_yz	[41]dPgdr/rho\n");
							}
							fprintf(p_1dfile,"%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G	%G\n",out_r[k],out1d[k][0],out1d[k][1],out1d[k][2],
									out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
									out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28],out1d[k][29],out1d[k][30],out1d[k][31],out1d[k][32],out1d[k][33],out1d[k][34],out1d[k][35],out1d[k][36],out1d[k][37],out1d[k][38],out1d[k][39]);


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
					
#endif
					
					/* We need to do it every time as it can be different for each domain at each level */
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						/* only need to free it for myID == 0 */
						free_1d_array(out_r);
#ifdef MPI_PARALLEL
					}/* End myID */
#endif
										
										
				}/* End if level and domain match for this domain */
			}/* End if this CPU works on this domain */
		} /* End loop over domains at level nl */
	}/* end loop all levels */
	
	
	
	return;
}





/*! \fn static void output_2d_binary(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate the azimuthal averaged 
 * 2D profiles: radial and vertical */
/* Keep the information along radial and vertical directions */


/* Transform to Cynlindrical coordinate */
/* No need to do this transformation if in cylindrical coordinate */
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
	
	int i,j,k, ig, kg,n;
	int tot1d,i1d,nrmx,nzmx,kdispG,kdispD,idispG,idispD,index;
	int dnum = pOut->num;
	
	double **out1d;
	double x1,x2,x3,press,darea;
	double distance;
	double *out_r, *out_x3;
	double vx, vy, vz, Fr01, Fr02, Fr03;
	
	double Pg1, Pg2;

	
	
#ifdef MPI_PARALLEL
	double *my_out1d;
	double *g_out1d;	
	int zproc;
	int ierr,myID_Comm_Domain;
	int nproc;
#endif
	
	/* For radiation case, we add, Er, Frx, Fry, Frz, Frz0, dFrz0/dz, Er*vz, dP/dz/rho, dB2/dz/rho, kappaes, kappap, fxx, fyy, fzz, fxy, fxz, fyz */
	
#ifdef MHD
  	tot1d=51;
#else
  	tot1d=41; 
#endif /* MHD */
	
	
	
	
	int is, ie, js, je, ks, ke;
	/* Loop over all Domain and level */
	int offset, nf, ifr;
	int ir, jr, kr;
	offset = Radghost - nghost;
	
	/* The frequency averaged quantities */
	Real Fr[3], Pr[6], Edd[6], Sigma[4];
	Real Er, Erk1, Erk0, Erj1, Erj0, Eri0, Eri1;
	
	DomainS *pD;
	GridS *pGrid;
	RadGridS *pRG;
	
	for(nl=0; nl<(pM->NLevels); nl++){
		for(nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
			if(pM->Domain[nl][nd].Grid != NULL){
				
				if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
					(pOut->ndomain == -1 || pOut->ndomain == nd)){
					
					pD = &(pM->Domain[nl][nd]);
					pGrid = pD->Grid;
					pRG = pD->RadGrid;

					nf = pRG->nf;					
					
					/* total number of radial grid zones for this Domain only */			
					nrmx = pD->Nx[0]; 
					nzmx = pD->Nx[2];
					darea = 1.0/(pD->Nx[1]);
					
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
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
						out_r = (double *) calloc_1d_array(nrmx,sizeof(double));
#ifdef MPI_PARALLEL
					}/* End myID */
#endif
					
					
					
					out1d = (double **) calloc_2d_array(nrmx*nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
					my_out1d = (double *) calloc_1d_array(nrmx*nzmx,sizeof(double));
					g_out1d = (double *) calloc_1d_array(nrmx*nzmx,sizeof(double));
					
					
#endif
					idispD = pD->Disp[0];
					idispG = pGrid->Disp[0];
					kdispD = pD->Disp[2];
					kdispG = pGrid->Disp[2];
					
					for (i=0; i<nrmx*nzmx; i++) {						
						for (i1d=0; i1d<tot1d; i1d++) {
							out1d[i][i1d] = 0.0;
						}
					}
					
					
					/* First calculate the x3 coordinate and save it to be dumped
					 by root in every 1d file */
					
					/* There will be some grids that not covered in the radial grids */
					/* This needs to be done for every CPU because we need this array to determine the position */				
					
					
					
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						
						for (i=0; i<nrmx; i++) {
							/* Vertical coordinate for this whole domain */
							x1 = pD->MinX[0] + (i + 0.5)*pGrid->dx1;
							out_r[i] = x1;
						}

						for (k=0; k<nzmx; k++) {
							/* Vertical coordinate for this whole domain */
							x3 = pD->MinX[2] + (k + 0.5)*pGrid->dx3;
							out_x3[k] = x3;
						}
#ifdef MPI_PARALLEL
					}/* End my ID */
#endif
					
					
					
					/* Compute 2d averaged variables */
				
					
					for (k=ks; k<=ke; k++) {
						kg=k+kdispG-kdispD-nghost;
						for (j=js; j<=je; j++) {
							for (i=is; i<=ie; i++) {
								ig=i+idispG-idispD-nghost;
								
								ir = i + offset;
								jr = j + offset;
								kr = k + offset;

								/* Average the radiation quantities over frequencies */
								
								for(n=0; n<3; n++)
									Fr[n] = 0.0;
								for(n=0; n<4; n++)
									Sigma[n] = 0.0;
								for(n=0; n<6; n++)
									Pr[n] = 0.0;
								
								Er = 0.0;
								Erk1 = 0.0;
								Erk0 = 0.0;
								Erj1 = 0.0;
								Erj0 = 0.0;
								Eri1 = 0.0;
								Eri0 = 0.0;
								
								for(ifr=0; ifr<nf; ifr++){
									for(n=0; n<3; n++)
										Fr[n] += 4.0 * PI * pRG->R[kr][jr][ir][ifr].H[n] * pRG->wnu[ifr];
									for(n=0; n<4; n++)
										Sigma[n] += pRG->R[kr][jr][ir][ifr].Sigma[n] * pRG->wnu[ifr];
									for(n=0; n<6; n++)
										Pr[n] += 4.0 * PI * pRG->R[kr][jr][ir][ifr].K[n] * pRG->wnu[ifr];
									
									Er += 4.0 * PI * pRG->R[kr][jr][ir][ifr].J * pRG->wnu[ifr];
									Erk1 += 4.0 * PI * pRG->R[kr+1][jr][ir][ifr].J * pRG->wnu[ifr];
									Erk0 += 4.0 * PI * pRG->R[kr-1][jr][ir][ifr].J * pRG->wnu[ifr];
									Erj1 += 4.0 * PI * pRG->R[kr][jr+1][ir][ifr].J * pRG->wnu[ifr];
									Erj0 += 4.0 * PI * pRG->R[kr][jr-1][ir][ifr].J * pRG->wnu[ifr];
									Eri1 += 4.0 * PI * pRG->R[kr][jr][ir+1][ifr].J * pRG->wnu[ifr];
									Eri0 += 4.0 * PI * pRG->R[kr][jr][ir-1][ifr].J * pRG->wnu[ifr];
								}
								
								for(n=0; n<6; n++)
									Edd[n] = Pr[n] / Er;
								
								
								cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
								distance = x1;								
								index = kg * nrmx + ig;

								/* Now this cell should fill [rpos] in out_r and out1d */
									
								vx = pGrid->U[k][j][i].M1 / pGrid->U[k][j][i].d;
								vy = pGrid->U[k][j][i].M2 / pGrid->U[k][j][i].d;
								vz = pGrid->U[k][j][i].M3 / pGrid->U[k][j][i].d;
								
								i1d=0;
								out1d[index][i1d] += pGrid->U[k][j][i].d;
								i1d++;
#ifdef ISOTHERMAL
								out1d[index][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
								press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#ifdef MHD
									   - expr_ME(pGrid,i,j,k)
#endif
									   ),TINY_NUMBER);


								
								Pg1           = MAX(Gamma_1*(pGrid->U[k-1][j][i].E - expr_KE(pGrid,i,j,k-1)
#ifdef MHD
														   - expr_ME(pGrid,i,j,k-1)
#endif
														   ),TINY_NUMBER);

								Pg2           = MAX(Gamma_1*(pGrid->U[k+1][j][i].E - expr_KE(pGrid,i,j,k+1)
#ifdef MHD
														   - expr_ME(pGrid,i,j,k+1)
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
								out1d[index][i1d] += pGrid->U[k][j][i].M1;
								i1d++;
								out1d[index][i1d] += pGrid->U[k][j][i].M2;
								i1d++;
								out1d[index][i1d] += pGrid->U[k][j][i].M3;
								/* kinetic */
								i1d++;
								out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
								i1d++;									
								out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;		
								i1d++;
								out1d[index][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
								i1d++;
								out1d[index][i1d] += expr_KE(pGrid,i,j,k);
								i1d++;
								out1d[index][i1d] += hst_rho_VrVphi(pGrid,i,j,k);
								/* vr */
								i1d++;
								out1d[index][i1d] += vx;
								/* vphi */
								i1d++;
								out1d[index][i1d] += vy;
								/* vz */
								i1d++;
								out1d[index][i1d] += vz;
#ifdef MHD
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
								out1d[index][i1d] += pGrid->U[k][j][i].M1 * x1;
									
								/* angular momentum */
								i1d++;
								out1d[index][i1d] += pGrid->U[k][j][i].M2 * distance;
								/* specific angular momentum */
								i1d++;
								out1d[index][i1d] += vy * distance;
									
								i1d++;
								out1d[index][i1d] += Er;
								i1d++;
								out1d[index][i1d] += Fr[0];
								i1d++;
								out1d[index][i1d] += Fr[1];
								i1d++;
								out1d[index][i1d] += Fr[2];
									/* co-moving flux */
								Fr01 = Fr[0] -(vx * (1.0 + Edd[0]) + vy * Edd[1] + vz * Edd[3]) * Er/Crat;
								Fr02 = Fr[1] -(vy * (1.0 + Edd[2]) + vx * Edd[1] + vz * Edd[4]) * Er/Crat;
								Fr03 = Fr[2] -(vz * (1.0 + Edd[5]) + vx * Edd[3] + vy * Edd[4]) * Er/Crat;
								
								i1d++;
								out1d[index][i1d] += Fr01;
								i1d++;
								out1d[index][i1d] += Fr02;
								i1d++;
								out1d[index][i1d] += Fr03;
									
								/* advection Er */
								i1d++;
								out1d[index][i1d] += Er * vx;
								/* kappes and kappaff */
									
								i1d++;
								out1d[index][i1d] += hst_sigmas(pGrid,i,j,k);
								i1d++;
								out1d[index][i1d] += hst_sigmaaP(pGrid,i,j,k);
								/* radiation work term */
								i1d++;
								out1d[index][i1d] += (Sigma[1] - Sigma[2]) * (Fr01 * vx + Fr02 * vy + Fr03 * vz); 
								/* dErdx */
								i1d++;
								out1d[index][i1d] += (Eri1 - Eri0)/(2.0 * pGrid->dx1);
								i1d++;
								out1d[index][i1d] += (Erj1 - Erj0)/(2.0 * pGrid->dx2);
								i1d++;
								out1d[index][i1d] += (Erk1 - Erk0)/(2.0 * pGrid->dx3);
								/* Eddington tensor */
								i1d++;
								out1d[index][i1d] += Edd[0];
								i1d++;
								out1d[index][i1d] += Edd[2];
								i1d++;
								out1d[index][i1d] += Edd[5];
								i1d++;
								out1d[index][i1d] += Edd[1];
								i1d++;
								out1d[index][i1d] += Edd[3];
								i1d++;
								out1d[index][i1d] += Edd[4];

								i1d++;
								out1d[index][i1d] += (Pg2-Pg1)/(2.0 * pGrid->dx3 * pGrid->U[k][j][i].d);
#ifdef MHD
								i1d++;
								out1d[index][i1d] += (expr_ME(pGrid,i,j,k+1) - expr_ME(pGrid,i,j,k-1))/(2.0 * pGrid->dx3 * pGrid->U[k][j][i].d);
#endif
									

								
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
					
	
#endif
					
					
					/* For parallel calculations, only the parent computes the average
					 * and writes the output. */
#ifdef MPI_PARALLEL
					if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif
						
						
						for (k=0; k<nrmx*nzmx; k++) {
							for (i1d=0; i1d<tot1d; i1d++) {
								
								out1d[k][i1d] *= darea;
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
#endif

					/* We need to do it every time as it can be different for each domain at each level */
#ifdef MPI_PARALLEL
					if (myID_Comm_Domain == 0) {
#endif
						/* only need to free it for myID == 0 */
						free_1d_array(out_r);
						free_1d_array(out_x3);
#ifdef MPI_PARALLEL
					}/* End myID */
#endif
					
									
					
				}/* End if level and domain match for this domain */
			}/* End if this CPU works on this domain */
		} /* End loop over domains at level nl */
	}/* end loop all levels */
	
	
	
	return;
}





