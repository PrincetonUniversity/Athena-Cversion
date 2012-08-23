#include "copyright.h"
/*==============================================================================
 * FILE: rt.c
 *
 * PURPOSE: Problem generator for galaxy model in a static ICM.  
 *   Gravitational pot. 
 *   is hardwired to match Roediger & Bruggen 2006. 
 *   
 *   
 * 
 *
 * FOR 3D:
 * Problem domain should be -1.4 < x < 1.4; -1.4 < y < 1.4, -.5 < z < 0.5
 * Use gamma=5/3
 * 
 * 
 *   
 *
 *
 * 
 *============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2() - random number generator from NR
 * reflect_ix3() - sets BCs on L-x3 (left edge) of grid used in 3D
 * reflect_ox3() - sets BCs on R-x3 (right edge) of grid used in 3D
 * grav_pot3() - gravitational potential for 3D problem
 * rps_ikb() -sets BC on L-x3 of grid in 3D
 *============================================================================*/

static Real ran2(long int *idum);
//static void reflect_ix3(GridS *pGrid);
//static void reflect_ox3(GridS *pGrid);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);
static void rps_ikb(GridS *pGrid);

/*===========================GALAXY DEFINING FUNCTIONS========================*/

Real av_den(Real cellwidth, Real z,Real xpos, Real ypos, Real zpos);
Real gas_vel(Real cellwidth, Real z, Real xpos, Real ypos, Real zpos, Real *temperature);
Real func1(Real zint, Real zicmP);
Real func2(Real zint,Real zicmP);
Real func3(Real zint, Real zicmP);
Real func4(Real zint,Real zicmP);
Real grav_potcalc(Real drcylin, Real x3g);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

static Real DensityUnits, LengthUnits,TemperatureUnits, TimeUnits, VelocityUnits;
static Real DiskRadius, mu;
static Real DiskPositionx, DiskPositiony, DiskPositionz, rBulge, MBulge, MgasScale, gScaleHeightR;
static Real gScaleHeightz, MSDisk, SDiskScaleHeightR, SDiskScaleHeightz, rDMConst, densDMConst;
static Real GravConst, Mpc, pi, SolarMass, mh, kboltz;
static Real densicm, ticm, Picm;
static Real AngularMomentumx, AngularMomentumy, AngularMomentumz;
static Real r2, drcyl, cellwidth;
static Real xuse, yuse,zuse;
static Real densinflow,tinflow,xvel_inflow,yvel_inflow,zvel_inflow;

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,iprob;
  long int iseed = -1;
  Real x1,x2,x3,lx,ly,lz; ///change these.
  Real DiskDensity, DiskVelocityMag;
  Real density, dens1, temp1;
  Real temperature, temperature0;
  			
  Real xpos, ypos, zpos, xpos1, ypos1, zpos1, zheight, drad; 
	
#ifdef MHD
  Real b0,angle;
#endif
  int ixs, jxs, kxs;

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

/* Read all inputs */
  densicm  = par_getd("problem","densicm");
  ticm = par_getd("problem","ticm");
//  Gamma = par_getd("problem","gamma");
  cellwidth = par_getd("problem","cellwidth");
	
/* Inputs with a default value */
  AngularMomentumx = par_getd_def("problem","AngularMomentumx",0.0);
  AngularMomentumy = par_getd_def("problem","AngularMomentumy",0.0);
  AngularMomentumz = par_getd_def("problem","AngularMomentumz",1.0);
  DiskPositionx = par_getd_def("problem","DiskPositionx",0.0);
  DiskPositiony = par_getd_def("problem","DiskPositiony",0.0);
  DiskPositionz = par_getd_def("problem","DiskPositionz",0.0);
  MgasScale = par_getd_def("problem","MgasScale",1.0e10);
  gScaleHeightR = par_getd_def("problem","gScaleHeightR",7.0e-3);
  gScaleHeightz = par_getd_def("problem","gScaleHeightz",4.0e-4);
  MSDisk = par_getd_def("problem","MSDisk",1.0e11);
  SDiskScaleHeightR = par_getd_def("problem","SDiskScaleHeightR",4.0e-3);
  SDiskScaleHeightz = par_getd_def("problem","SDiskScaleHeightz",2.5e-4);
  MBulge = par_getd_def("problem","MBulge",1.0e10);
  rBulge = par_getd_def("problem","rBulge",4.0e-4);
  rDMConst = par_getd_def("problem","rDMConst",2.3e-2);
  densDMConst	= par_getd_def("problem","densDMConst",3.81323e-25);
  
  densinflow = par_getd_def("problem","densinflow",densicm);
  tinflow = par_getd_def("problem","tinflow",ticm);
  xvel_inflow = par_getd_def("problem","xvel_inflow",0.0);	
  yvel_inflow = par_getd_def("problem","yvel_inflow",0.0);
  zvel_inflow = par_getd_def("problem","zvel_inflow",0.0);

  printf("rDMConst %g\n",rDMConst);
  DiskRadius = 1.0;
  mu = 0.6;
  DensityUnits = 1.0e-29;
  TimeUnits = 3.086e14;
  LengthUnits = 8.0236e22;
  VelocityUnits = LengthUnits/TimeUnits;
  TemperatureUnits = LengthUnits*LengthUnits/TimeUnits/TimeUnits*DensityUnits; //1.0;
  GravConst = 6.67e-8;
  SolarMass = 1.99e33;
  mh = 1.6733e-24;
  kboltz = 1.380658e-16;
  Mpc = 3.086e24;
  pi = 3.1415926536;
  Picm = densicm * DensityUnits / (mu*mh) * kboltz * ticm;
  printf("Picm %g Units %g\n",Picm,Picm/TemperatureUnits);
  temperature0 = temperature = ticm;
	
/* Read magnetic field strength, angle [should be in degrees, 0 is along +ve
 * X-axis (no rotation)] */
#ifdef MHD
  b0 = par_getd("problem","b0");
  angle = par_getd("problem","angle");
  angle = (angle/180.)*PI;
#endif


/* 3D PROBLEM ----------------------------------------------------------------*/
	
	
  if (pGrid->Nx[2] > 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	pGrid->U[k][j][i].d = densicm;
	pGrid->U[k][j][i].E = kboltz * ticm/TemperatureUnits/((Gamma-1.0)*mh*mu) * densicm*DensityUnits;
	//printf("E ICM: %g\n",pGrid->U[k][j][i].E);
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
	pGrid->U[k][j][i].M3 = 0.0;
	
		  /*If we are in the galaxy*/
	if (sqrt(pow((x1-DiskPositionx),2.0) + pow((x2-DiskPositiony),2.0) + pow((x3-DiskPositionz),2.0)) < DiskRadius) {

//			extern Real drcyl;
//		    extern Real xuse, yuse,zuse;
			/*compute position*/
			
			xpos = x1 - DiskPositionx;
			ypos = x2 - DiskPositiony;
			zpos = x3 - DiskPositionz;
			
			zheight = AngularMomentumx * xpos + AngularMomentumy * ypos + AngularMomentumz * zpos;
			
			xpos1 = xpos - zheight*AngularMomentumx;
			ypos1 = ypos - zheight*AngularMomentumy;
			zpos1 = zpos - zheight*AngularMomentumz;
			drad = sqrt(xpos1*xpos1 + ypos1*ypos1 + zpos1*zpos1);
			drcyl = drad;
			
			/*Normalize the vectors pointing along plane of disk*/
			
			xpos1 = xpos1/drad;
			ypos1 = ypos1/drad;
			zpos1 = zpos1/drad;
		    xuse = x1;
		    yuse = x2;
		    zuse = x3;
			DiskVelocityMag = gas_vel(cellwidth,zheight*LengthUnits,xpos1*drad*LengthUnits,
									  ypos1*drad*LengthUnits,zpos1*drad*LengthUnits, &temperature);
			
			temperature0 = temperature;
//			if (gScaleHeightz*Mpc/LengthUnits > CellWidth[0][0])/*I think that this does not need an if statement, because there is no AMR*/
//		    {
				
			dens1=av_den(cellwidth*LengthUnits, zheight*LengthUnits,xpos1*drad*LengthUnits,ypos1*drad*LengthUnits,
							 zpos1*drad*LengthUnits);
		    //printf("dens1 %g\n",dens1);

		    if (dens1 <= densicm){
				temperature0 = ticm;
				dens1 = densicm;
				density = dens1;
			}
		  
		    if (dens1 > densicm) {
			  if (drcyl <= 0.026*Mpc/LengthUnits) {
				  density = dens1;
				  //printf("temperature %g %g %g, dens1 %g, diskvel %g\n", temperature0,temperature,&temperature, dens1, DiskVelocityMag);
				  /*this is where I would insert the "color" parameter if I was doing so*/
			  }
			//  printf("dens1 %g density %g vrot %g\n",dens1,density, DiskVelocityMag);
				pGrid->U[k][j][i].d = density; //dens1; /*pretty sure that this is in code units*/
		  /* Compute velocity magnitude (divided by drad). 
		   This assumes PointSourceGravityPosition and Disk center are the same. */
		  
		  /* Compute velocty: L x r_perp. */  /*to get to Momentum units I will have to multiply by U[k][j][i].d*/
		  
		  pGrid->U[k][j][i].M1 = DiskVelocityMag*(AngularMomentumy*zpos1 -AngularMomentumz*ypos1);
		  
		  pGrid->U[k][j][i].M2 = DiskVelocityMag*(AngularMomentumz*xpos1 -   AngularMomentumx*zpos1);
		  
		  pGrid->U[k][j][i].M3 = DiskVelocityMag*(AngularMomentumx*ypos1 -  AngularMomentumy*xpos1);

		  pGrid->U[k][j][i].M1 *= density;
		  pGrid->U[k][j][i].M2 *= density;
		  pGrid->U[k][j][i].M3 *= density;
		
			
	
		pGrid->U[k][j][i].E = kboltz * temperature0/TemperatureUnits/((Gamma-1.0)*mh*mu) * pGrid->U[k][j][i].d*DensityUnits;
	
		  
		  pGrid->U[k][j][i].E += 0.5*(pow(pGrid->U[k][j][i].M3,2)+pow(pGrid->U[k][j][i].M2,2)+pow(pGrid->U[k][j][i].M1,2))/pGrid->U[k][j][i].d;
		  if (pGrid->U[k][j][i].E < kboltz * ticm/TemperatureUnits/((Gamma-1.0)*mh*mu) * densicm*DensityUnits) {
		    pGrid->U[k][j][i].E = kboltz * ticm/TemperatureUnits/((Gamma-1.0)*mh*mu) * densicm*DensityUnits;
		    printf("had a low energy!");
		  }

		 
			}
	}
#ifdef MHD
        switch(iprob){
        case 3: /* B only in light fluid, do not add B^2 to E, total P const */
          if (x3 <= 0.0) {
            pGrid->B1i[k][j][i] = b0;
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0;
            pGrid->U[k][j][i].B1c = b0;
          }
          break;
        case 4: /* discontinuous rotation of B by angle at interface */
          if (x3 <= 0.0) {
            pGrid->B1i[k][j][i] = b0;
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0;
            pGrid->U[k][j][i].B1c = b0;
            pGrid->U[k][j][i].E += 0.5*b0*b0;
          }
          else {
            pGrid->B1i[k][j][i] = b0*cos(angle);
            pGrid->B2i[k][j][i] = b0*sin(angle);
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0*cos(angle);
            if (j == je) pGrid->B2i[k][je+1][i] = b0*sin(angle);
            pGrid->U[k][j][i].B1c = b0*cos(angle);
            pGrid->U[k][j][i].B2c = b0*sin(angle);
            pGrid->U[k][j][i].E += 0.5*b0*b0;
          }
          break;
        case 5: /* rotation of B by angle over distance L_rot at interface */
          if (x3 <= (-L_rot/2.0)) {
            pGrid->B1i[k][j][i] = b0;
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0;
            pGrid->U[k][j][i].B1c = b0;
            pGrid->U[k][j][i].E += 0.5*b0*b0;
          }
          else if (x3 >= (L_rot/2.0)) {
            pGrid->B1i[k][j][i] = b0*cos(angle);
            pGrid->B2i[k][j][i] = b0*sin(angle);
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0*cos(angle);
            if (j == je) pGrid->B2i[k][je+1][i] = b0*sin(angle);
            pGrid->U[k][j][i].B1c = b0*cos(angle);
            pGrid->U[k][j][i].B2c = b0*sin(angle);
            pGrid->U[k][j][i].E += 0.5*b0*b0;
          }
          else {
            fact = ((L_rot/2.0)+x3)/L_rot;
            pGrid->B1i[k][j][i] = b0*cos(fact*angle);
            pGrid->B2i[k][j][i] = b0*sin(fact*angle);
            if (i == ie) pGrid->B1i[k][j][ie+1] = b0*cos(fact*angle);
            if (j == je) pGrid->B2i[k][je+1][i] = b0*sin(fact*angle);
            pGrid->U[k][j][i].B1c = b0*cos(fact*angle);
            pGrid->U[k][j][i].B2c = b0*sin(fact*angle);
            pGrid->U[k][j][i].E += 0.5*b0*b0;
          }

          break;
        default:
          pGrid->B1i[k][j][i] = b0;
          if (i == ie) pGrid->B1i[k][j][ie+1] = b0;
          pGrid->U[k][j][i].B1c = b0;
          pGrid->U[k][j][i].E += 0.5*b0*b0;
        }
#endif
      }
    }
  }

/* Enroll gravitational potential to give accn in z-direction for 3D
 * Use special boundary condition routines.  In 3D, gravity is in the
 * z-direction, so special boundary conditions needed for x3
 */

  StaticGravPot = grav_pot3;

  if (pDomain->Disp[2] == 0) bvals_mhd_fun(pDomain, left_x3,  rps_ikb);
  /*  if (pDomain->MaxX[2] == pDomain->RootMaxX[2])
      bvals_mhd_fun(pDomain, right_x3, reflect_ox3);*/
  
	
  } /* end of 3D initialization */

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

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/*
 * 'problem_read_restart' must enroll special boundary value functions,
 *    and initialize gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  int nl,nd;

  if (pM->Nx[2] > 1) {
    StaticGravPot = grav_pot3;
    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        bvals_mhd_fun(&(pM->Domain[nl][nd]), left_x3,  rps_ikb);
	//        bvals_mhd_fun(&(pM->Domain[nl][nd]), right_x3, reflect_ox3);
      }
    }
  }

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
}

void Userwork_after_loop(MeshS *pM)
{
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*-----------------------------------------------------------------------------
 * rps_ikb:  i:Inflow, k: x3, b: boundary--inflow of stripping ICM
 */

static void rps_ikb(GridS *pGrid)
{

  int i = 0, j = 0;
  int is,ie,js,je,ks,k,il,iu,jl,ju;
  Real d0,e0,u0,v0,w0;

  d0 = densinflow;
  e0 =  kboltz * tinflow/TemperatureUnits/((Gamma-1.0)*mh*mu) * densinflow*DensityUnits;
  u0 = xvel_inflow/VelocityUnits;
  v0 = yvel_inflow/VelocityUnits;
  w0 = zvel_inflow/VelocityUnits;
  
  ks = pGrid->ks;
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;
  
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        //dmr.c calls cc_pos(pGrid,i,j,ks,&x1,&x2,&x3);  I don't think I need it!
	pGrid->U[ks-k][j][i].d  = d0;
	pGrid->U[ks-k][j][i].M1 = d0*u0;
        pGrid->U[ks-k][j][i].M2 = d0*v0;
        pGrid->U[ks-k][j][i].M3 = d0*w0;
        pGrid->U[ks-k][j][i].E  = e0 + 0.5*d0*(u0*u0 + v0*v0 + w0*w0);
      }
    }
  }
}

/*------------------------------------------------------------------------------
 * reflect_ix3: special reflecting boundary functions in x3 for 2D sims
 */

static void reflect_ix3(GridS *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ks-k][j][i]    =  pGrid->U[ks+(k-1)][j][i];
        pGrid->U[ks-k][j][i].M3 = -pGrid->U[ks-k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ks-k][j][i].E +=
          pGrid->U[ks+(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks+(k-1)][j][i];
      }
    }
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks+(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks+(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

/*------------------------------------------------------------------------------
 * reflect_ox3: special reflecting boundary functions in x3 for 3D sims
 */

static void reflect_ox3(GridS *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k ,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;
  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i]    =  pGrid->U[ke-(k-1)][j][i];
        pGrid->U[ke+k][j][i].M3 = -pGrid->U[ke+k][j][i].M3; /* reflect 3-mom. */
        pGrid->U[ke+k][j][i].E -=
          pGrid->U[ke-(k-1)][j][i].d*0.1*(2*k-1)*pGrid->dx3/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke-(k-1)][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke-(k-1)][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke-(k-1)][j][i];
      }
    }
  }
#endif

  return;
}

/*------------------------------------------------------------------------------
 * grav_pot: Gravitational potentials 
 */

static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
//	extern Real AngularMomentumx, AngularMomentumy, AngularMomentumz;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst;
//	extern Real SDiskScaleHeightR, SDiskScaleHeightz;
//	extern Real LengthUnits,DensityUnits,TimeUnits;
//	extern Real GravConst, Mpc, pi,SolarMass;
	Real x1p,x2p,x3p,x1c,x2c,x3c,rcyl;  /*position in plane of disk*/
	Real zheight,radiusgp; /*spherical position*/
	Real Epotbulge, Epotdisk, Epotdm,Epottotal;
	/*AngularMomentumx,y,z x1c,x2c,x3c*/
	/*MBulge, rBulge, MSDisk,SDiskScaleHeightR, SDiskScaleHeightz, densDMConst,rDMConst*/
	x1c = DiskPositionx;
	x2c = DiskPositiony;
	x3c = DiskPositionz;
	zheight = AngularMomentumx*(x1-x1c) + AngularMomentumy*(x2-x2c) + AngularMomentumz*(x3-x3c);
	x1p = (x1-x1c) - zheight*AngularMomentumx;
	x2p = (x2-x2c) - zheight*AngularMomentumy;
	x3p = (x3-x3c) - zheight*AngularMomentumz;
	radiusgp = sqrt(x1p*x1p + x2p*x2p + x3p*x3p + zheight*zheight)*LengthUnits/Mpc;
	rcyl = sqrt(x1p*x1p + x2p*x2p + x3p*x3p)*LengthUnits/Mpc;
	
	Epotbulge = -MBulge*SolarMass / ((radiusgp + rBulge)*Mpc);
	Epotdisk = -MSDisk*SolarMass / (pow(pow(rcyl,(2.0)) + pow(SDiskScaleHeightR + pow((pow(zheight*LengthUnits/Mpc,(2.0)) + pow(SDiskScaleHeightz,(2.0))),(0.5)),(2.0)),(0.5))*Mpc);
	Epotdm = -pi * densDMConst * pow(rDMConst*Mpc,(2.0)) * (pi - 2.0*(1.0 + (rDMConst/radiusgp)) * atan(radiusgp/rDMConst) + 2.0*(1.0 + (rDMConst/radiusgp)) * log(1.0 + (radiusgp/rDMConst)) - (1-(rDMConst/radiusgp)) * log(1.0 + pow((radiusgp/rDMConst),(2.0))));
	if (radiusgp == 0.0) {
		Real radiusgpu = 0.000001;
		Epotdm = -pi * densDMConst * pow(rDMConst*Mpc,(2.0)) * (pi - 2.0*(1.0 + (rDMConst/radiusgpu)) * atan(radiusgpu/rDMConst) 
						+ 2.0*(1.0 + (rDMConst/radiusgpu)) * log(1.0 + (radiusgpu/rDMConst)) - (1-(rDMConst/radiusgpu)) * log(1.0 + pow((radiusgpu/rDMConst),(2.0))));
		//Epotdm = -pi*pi * densDMConst * pow(rDMConst*Mpc,(2.0));
	}
	Epottotal = Epotbulge + Epotdisk + Epotdm;
	if (Epotbulge > 0.0) printf("Epotbulge > 0.0!");
	if (Epotdisk > 0.0) printf("Epotdisk > 0.0!");
	if (Epotdm > 0.0) printf("Epotdm > 0.0!");
	//if (Epotbulge+Epotdisk+Epotdm == nan) {
	//printf("Epotbulge %g, Epotdisk %g, Epotdm %g, x1 %g, x2 %g x3 %g\n",Epotbulge,Epotdisk,Epotdm, x1,x2,x3);
//	}
    return GravConst*Epottotal/pow(LengthUnits,2)/pow(TimeUnits,-2); ///DensityUnits/pow(LengthUnits,5)/pow(TimeUnits,-2);  //codeunits
	//using the "-" here gives the wrong direction for gravity!  huh.
}
Real grav_potcalc(Real drcylin, Real x3g)
{
//	extern Real AngularMomentumx, AngularMomentumy, AngularMomentumz;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst;
//	extern Real SDiskScaleHeightR, SDiskScaleHeightz;
//	extern Real LengthUnits,DensityUnits,TimeUnits;
//	extern Real GravConst, Mpc, pi,SolarMass;
	Real x1p,x2p,x3p,x1c,x2c,x3c,rcylg;  /*position in plane of disk*/
	Real zheight,radiusggp; /*spherical position*/
	Real Epotbulge, Epotdisk, Epotdm,Epottotal;
	/*AngularMomentumx,y,z x1c,x2c,x3c*/
	/*MBulge, rBulge, MSDisk,SDiskScaleHeightR, SDiskScaleHeightz, densDMConst,rDMConst*/
	radiusggp = sqrt(drcylin*drcylin+x3g*x3g)*LengthUnits/Mpc;
	rcylg = drcylin*LengthUnits/Mpc;
	
	Epotbulge = -MBulge*SolarMass / ((radiusggp + rBulge)*Mpc);
	Epotdisk = -MSDisk*SolarMass / (pow(pow(rcylg,(2.0)) + pow(SDiskScaleHeightR + pow((pow(x3g*LengthUnits/Mpc,(2.0)) 
																						+ pow(SDiskScaleHeightz,(2.0))),(0.5)),(2.0)),(0.5))*Mpc);
	Epotdm = -pi * densDMConst * pow(rDMConst*Mpc,(2.0)) * (pi - 2.0*(1.0 + (rDMConst/radiusggp))
			* atan(radiusggp/rDMConst) + 2.0*(1.0 + (rDMConst/radiusggp)) * log(1.0 + (radiusggp/rDMConst)) - (1-(rDMConst/radiusggp)) 
															* log(1.0 + pow((radiusggp/rDMConst),(2.0))));
	if (radiusggp == 0.0) {
		Real radiusggpu = 0.000001;
		Epotdm = -pi * densDMConst * pow(rDMConst*Mpc,(2.0)) * (pi - 2.0*(1.0 + (rDMConst/radiusggpu))
																* atan(radiusggpu/rDMConst) + 2.0*(1.0 + (rDMConst/radiusggpu)) * log(1.0 + (radiusggpu/rDMConst)) - (1-(rDMConst/radiusggpu)) 
																* log(1.0 + pow((radiusggpu/rDMConst),(2.0))));
		//Epotdm = -pi*pi * densDMConst * pow(rDMConst*Mpc,(2.0));
	}
	Epottotal = Epotbulge + Epotdisk + Epotdm;
	if (Epotbulge > 0.0) printf("Epotbulge > 0.0!");
	if (Epotdisk > 0.0) printf("Epotdisk > 0.0!");
	if (Epotdm > 0.0) printf("Epotdm > 0.0!");
//	if (rcylg < 1e-2 && x3g *LengthUnits/Mpc < 4e-4) {
//		printf("Epotbulge %g Epotdisk %g Epotdm %g drcyl %g x3g %g\n",GravConst*Epotbulge,GravConst*Epotdisk,GravConst*Epotdm, rcylg,x3g*LengthUnits/Mpc);
//	}
    return GravConst*Epottotal/pow(LengthUnits,2)/pow(TimeUnits,-2);  //codeunits
	//using the "-" here gives the wrong direction for gravity!  huh.
}

/*_______________________________________________________________________*/
//Written by:  Stephanie Tonnesen
//Date:  November 2006
/*copied over to solve for density, temperature and velocity of galactic gas*/

//void nrerror(char error_text[])

//Real *vector(int nl,int nh)

//void free_vector(Real *v,int nl,int nh)


//Will be called by qromb to find the pressure at every point in disk.

//#define FUNC(x) ((*func)(x))

//Real trapzdd(Real (*func)(Real), Real a, Real b, int n)

//also called by qromb.

//void polint(Real xa[],Real ya[],int n,Real x,Real *y,Real *dy)
//main integration routine called in vel_gas.

//#define EPS 1.0e-6
//#define JMAX 60
//#define JMAXP JMAX+1
//#define K 12  //was 7 //30 never finished

//Real qromb(Real (*func)(Real), Real a, Real b)


//the two functions integrated by qromb

Real func1(Real zint, Real zicmP)
{
//	extern Real gScaleHeightR, gScaleHeightz;
//	extern Real LengthUnits,TimeUnits;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst, MgasScale;
//	extern Real GravConst, Mpc, pi, SolarMass;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real drcyl,xuse,yuse,zuse;
//	extern Real cellwidth,Picm;
	Real grav_potcalc(Real drcylin, Real x3g);
	Real numstepsP = 2000.;
	Real stepsizeP = (zicmP-zint)/numstepsP;
	//printf("stepsize %g\n",stepsizeP);
	Real Pressuresolve = 0.0;
	Real iP = 0.;
	for (iP=0.;iP<=numstepsP;iP++) {
		Pressuresolve += (MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)
			/cosh(fabs(zicmP-iP*stepsizeP-0.5*stepsizeP)/gScaleHeightz/Mpc)
			* (grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
						  *pow(LengthUnits,2)*pow(TimeUnits,-2));///fabs(stepsizeP)); //*fabs(zint)/sqrt(pow(zint,2)+pow(drcyl*LengthUnits,2)) 
		/*Pressuresolve +=(grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
		*pow(LengthUnits,2)*pow(TimeUnits,-2);*/
	}
	Pressuresolve = Pressuresolve;//*MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)
								  // /cosh(fabs(zint)/gScaleHeightz/Mpc);// * (zicmP-zint);
	//printf("Pressure %g, zicm %g, zin %g, zfinal %g\n",Pressuresolve,zicmP,zint,zicmP-50.*stepsizeP);
	return (Pressuresolve);
}	
Real func2(Real zint, Real zicmP)
{
//	extern Real gScaleHeightR, gScaleHeightz;
//	extern Real LengthUnits,TimeUnits;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst, MgasScale;
//	extern Real GravConst, Mpc, pi, SolarMass;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real drcyl,xuse,yuse,zuse;
//	extern Real cellwidth,Picm;
	Real grav_potcalc(Real drcylin, Real x3g);
	Real numstepsP = 2000.;
	Real stepsizeP = (zicmP-zint)/numstepsP;
	//printf("stepsize %g\n",stepsizeP);
	Real Pressuresolve = 0.0;
	Real iP = 0;
	for (iP=0;iP<=numstepsP;iP++) {
		Pressuresolve += (MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)
						  /cosh(fabs(zicmP-iP*stepsizeP-0.5*stepsizeP)/gScaleHeightz/Mpc)*(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc))))
						  * (grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
						  *pow(LengthUnits,2)*pow(TimeUnits,-2));///fabs(stepsizeP)); //*fabs(zint)/sqrt(pow(zint,2)+pow(drcyl*LengthUnits,2)) 
		/*Pressuresolve += (grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(drcyl,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
								   *pow(LengthUnits,2)*pow(TimeUnits,-2);///fabs(stepsizeP)); */
	}
	Pressuresolve = Pressuresolve;//*MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)
									// /cosh(fabs(zint)/gScaleHeightz/Mpc)*(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc))));// * (zicmP-zint);
	//printf("Pressure %g, zicm %g, zin %g, zfinal %g\n",Pressuresolve,zicmP,zint,zicmP-50.*stepsizeP);
return (Pressuresolve);
}

//need to repeat the same functions but with a new r

Real func3(Real zint, Real zicmP)
{
//	extern Real gScaleHeightR, gScaleHeightz;
//	extern Real LengthUnits,TimeUnits;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst, MgasScale;
//	extern Real GravConst, Mpc, pi, SolarMass;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real r2,xuse,yuse;
//	extern Real cellwidth,Picm;	
	Real grav_potcalc(Real drcylin, Real x3g);
	Real numstepsP = 2000;
	Real stepsizeP = (zicmP-zint)/numstepsP;
	Real Pressuresolve = 0.0;
	Real iP = 0;
	for (iP=0;iP<=numstepsP;iP++) {
		Pressuresolve += (MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)/
			cosh(fabs(zicmP-iP*stepsizeP-0.5*stepsizeP)/gScaleHeightz/Mpc)
			* (grav_potcalc(r2/LengthUnits,fabs(zicmP-iP*stepsizeP)/LengthUnits) -grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
						  *pow(LengthUnits,2)*pow(TimeUnits,-2));///fabs(stepsizeP));//*fabs(zint)/sqrt(pow(zint,2)+pow(r2,2))
		/*Pressuresolve +=(grav_potcalc(r2/LengthUnits,fabs(zicmP-iP*stepsizeP)/LengthUnits) -grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
		*pow(LengthUnits,2)*pow(TimeUnits,-2);*/
	}
	Pressuresolve = Pressuresolve;//*MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)/
									// cosh(fabs(zint)/gScaleHeightz/Mpc);// * (zicmP-zint);
	return (Pressuresolve);
}

Real func4(Real zint, Real zicmP)
{
//	extern Real gScaleHeightR, gScaleHeightz;
//	extern Real LengthUnits,TimeUnits;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst, MgasScale;
//	extern Real GravConst, Mpc, pi, SolarMass;
//	extern Real DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real r2,xuse,yuse,zuse;
//	extern Real cellwidth,Picm;
	Real grav_potcalc(Real drcylin, Real x3g);
	Real numstepsP = 2000.;
	Real stepsizeP = (zicmP-zint)/numstepsP;
	//printf("stepsize %g\n",stepsizeP);
	Real Pressuresolve = 0.0;
	Real iP = 0;
	for (iP=0;iP<=numstepsP;iP++) {
		Pressuresolve += (MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)
						  /cosh(fabs(zicmP-iP*stepsizeP-0.5*stepsizeP)/gScaleHeightz/Mpc)*(0.5*(1.0+cos(pi*(r2-0.02*Mpc)/(0.006*Mpc))))
						  * (grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
						  *pow(LengthUnits,2)*pow(TimeUnits,-2));///fabs(stepsizeP)); //*fabs(zint)/sqrt(pow(zint,2)+pow(drcyl*LengthUnits,2))
		/*Pressuresolve += (grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP)/LengthUnits)) - grav_potcalc(r2/LengthUnits,fabs((zicmP-iP*stepsizeP-stepsizeP)/LengthUnits)))
		*pow(LengthUnits,2)*pow(TimeUnits,-2);*/
	}
	Pressuresolve = Pressuresolve;//*MgasScale*SolarMass/(2*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)
									// /cosh(fabs(zint)/gScaleHeightz/Mpc)*(0.5*(1.0+cos(pi*(r2-0.02*Mpc)/(0.006*Mpc))));// * (zicmP-zint);
	//printf("Pressure %g, zicm %g, zin %g, zfinal %g\n",Pressuresolve,zicmP,zint,zicmP-50.*stepsizeP);
return (Pressuresolve);
}

Real gas_vel(Real cellwidth, Real z, Real xpos, Real ypos, Real zpos, Real *temperature)
{

/* Five steps:  1) Solve for zicm (at two nearby radii), to use as a limit to integral in 2) solve integral for Pressure at (r,z) (and at a nearby r). 
 3) Use P to solve for Temp(r,z), send to an array. 4) solve dP/dr using the two P values.  
 5) solve for vrot by balancing radial gravity force, dP/dr, and centrifugal force */

	//extern float AngularMomentumx, AngularMomentumxy, AngularMomentumxz;
	//extern float DiskPositionx, DiskPositiony, DiskPositionz;

//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst, MgasScale;
//	extern Real SDiskScaleHeightR, SDiskScaleHeightz, gScaleHeightR, gScaleHeightz;
//	extern Real LengthUnits,DensityUnits,TimeUnits, VelocityUnits;
//	extern Real GravConst, Mpc, pi,SolarMass, mh, kboltz;
//	extern Real DiskPositionx, DiskPositiony,DiskPositionz;
//	extern Real drcyl,xuse,yuse;
//	extern Real densicm, Picm;
	Real zicm,zicm2,zicmf=0.0,zsmall=0.0,zicm2f=0.0;
	Real Pressure,Pressure2,Pressure0,Pressure20;
	r2=(drcyl+0.01*cellwidth)*LengthUnits;
	// printf("r2,drcyl,cellwidth,LU=%g %g %g %g\n", r2, drcyl, cellwidth, LengthUnits);
	Real func1(Real zint,Real zicmP);       //(density times Stellar bulge force)
	Real func2(Real zint,Real zicmP);		 //(density times stellar disk force)
	Real func3(Real zint, Real zicmP);       //func1 but for r2
	Real func4(Real zint,Real zicmP);      //func2 but for r2
/* static Real grav_pot3(const Real x1, const Real x2, const Real x3); */
	Real grav_potcalc(Real drcylin, Real x3g);
	Real zint; 		 //zint variable to be integrated
	Real FdPdR;
	Real FtotR;
	Real FtotRpot;
	Real FDMR;
	Real denuse;
	Real rsph;
	rsph=sqrt(pow(drcyl*LengthUnits,2)+pow(z,2));
	Real vrot;
	
	Pressure = 0.0;
	Pressure2 = 0.0;
	
	if (fabs(drcyl*LengthUnits/Mpc) <= 0.02) {
		zicm=densicm*DensityUnits/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/
								   cosh(drcyl*LengthUnits/gScaleHeightR/Mpc));
		zicm=log(1.0/zicm+sqrt((1.0/pow(zicm,2))-1.0));
		zicm=fabs(zicm*gScaleHeightz*Mpc);
		
		//printf("zicm = %g, drcyl = %g\n", zicm/Mpc, drcyl*LengthUnits/Mpc);
		zicm2=densicm*DensityUnits/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/
									cosh(r2/gScaleHeightR/Mpc));
		zicm2=log(1.0/zicm2+sqrt((1.0/pow(zicm2,2))-1.0));
		zicm2=fabs(zicm2*gScaleHeightz*Mpc);
		
//		printf("zicm = %g, zicm2 = %g, drcyl = %g\n", zicm/Mpc, zicm2/Mpc, drcyl*LengthUnits/Mpc);
		
//		Pressure0 = qromb(func1,fabs(0.0),fabs(zicm));
//		Pressure20 = qromb(func1,fabs(0.0),fabs(zicm2));
//		printf("Pressure0 %g, Pressure20 %g \n",Pressure0,Pressure20);
		
		
		if (fabs(z) <= fabs(zicm)) {
			Pressure= func1(fabs(z), fabs(zicm)); // + qromb(func2, fabs(zicm), fabs(z));
			Pressure2= func3(fabs(z), fabs(zicm2)); // + qromb(func4, fabs(zicm2), fabs(z));
		}
	}
	else {
		if (fabs(drcyl*LengthUnits/Mpc) <= 0.026) {
			zicmf=densicm*DensityUnits/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)*(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc)))));
			zicm=log(1.0/zicmf+sqrt((1.0/pow(zicmf,2))-1.0));
			zicm=fabs(zicm*gScaleHeightz*Mpc);
			if (zicmf > 1.0) {
				zicm = 0.0;
			}
			
			//printf("zicm = %g, drcyl = %g\n", zicm/Mpc, drcyl*LengthUnits/Mpc);
			zicm2f=densicm*DensityUnits/(MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(r2/gScaleHeightR/Mpc)*(0.5*(1.0+cos(pi*(r2-0.02*Mpc)/(0.006*Mpc)))));
			zicm2=log(1.0/zicm2f+sqrt((1.0/pow(zicm2f,2))-1.0));
			zicm2=fabs(zicm2*gScaleHeightz*Mpc);
			if (zicm2f > 1.0) {
				zicm2 = 0.0;
			}
			
			if (densicm*DensityUnits >= (MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(z)/gScaleHeightz/Mpc)*(0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc))))) && fabs(z) < zicm) {
				printf("small density zicm = %g, z = %g\n", zicm/Mpc, z/Mpc);
			}
			if (fabs(z) <= fabs(zicm)) {
		//	if ((fabs(zicm)/Mpc-fabs(z)/Mpc) > 8.e-5) {
			//	printf("zicm %g, zicm2 %g, z %g, diff %g\n",zicm/Mpc,zicm2/Mpc,z/Mpc,(fabs(zicm)/Mpc-fabs(z)/Mpc));
				Pressure= (func2(fabs(z), fabs(zicm)));
				Pressure2 = 0.0;
			    if (fabs(z) <= fabs(zicm2)) {
				    Pressure2=(func4(fabs(z), fabs(zicm2)));
				}
				
			}
		}
	}
	
	denuse = av_den(cellwidth*LengthUnits, z, xpos, ypos, zpos)*DensityUnits;
	if (Pressure < 0.0 && fabs(drcyl)*LengthUnits/Mpc <= 0.026 && fabs(z) <= fabs(zicm)) {
		printf("neg pressure:  P = %g, z = %g, zicm = %g, P2 = %g, r = %g, x1 = %g, x2 = %g\n", Pressure, z/Mpc, zicm/Mpc, Pressure2, drcyl*LengthUnits/Mpc,xuse,yuse);
	}
	if (fabs(drcyl)*LengthUnits/Mpc >= 0.026 || fabs(zicm) <= fabs(z)){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm*DensityUnits;
	}
	if (Pressure2 <= 0.0 && Pressure <= 0.0){
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm*DensityUnits;
	}
	if (Pressure2 <= 0.0){
		Pressure2 = 0.0;
	}
	if (Pressure <= 0.0) {
		Pressure = 0.0;
		Pressure2 = 0.0;
		denuse = densicm*DensityUnits;
	}
	//if (Pressure2 < Pressure) {
	//   Pressure2 = Pressure;
	//  }
	if (denuse < densicm*DensityUnits) {
		printf("denuse small:  %g\n", denuse);
	}
	*temperature=0.6*mh*(Picm+Pressure)/(kboltz*denuse);
	//if (*temperature < 5e4) {
	//	*temperature=mh*(Picm+Pressure)/(kboltz*denuse);
	  //printf("temp %g, density %g, Picm %g, Pressure %g, Pressure2 %g, drcyl %g, r2 %g, cellwidth %g\n",0.6*mh*(Picm+Pressure)/(kboltz*denuse),denuse, Picm, Pressure,Pressure2, drcyl, r2/LengthUnits, cellwidth);
	//}
	FdPdR = (Pressure2 - Pressure)/(r2-drcyl*LengthUnits)/(av_den(cellwidth*LengthUnits, z, xpos, ypos, zpos)*DensityUnits);

	/*to get the total force in the R (cylindrical) direction, include stars and DM*/
	
	FtotRpot = -(grav_potcalc(r2/LengthUnits,fabs(z)/LengthUnits)-grav_potcalc(drcyl,fabs(z)/LengthUnits))*pow(LengthUnits,2)*pow(TimeUnits,-2)/(r2-drcyl*LengthUnits);
	//if (FtotR != FtotRpot) {
	//	printf("FtotR %g, FtotRpot %g,drcyl %g, x3 %g \n",FtotR,FtotRpot,drcyl*LengthUnits/Mpc, z/Mpc);
	//}
	if (FtotRpot > 0.0) {
		printf("FDMR = %g  FtotR = %g\n", FDMR, FtotR);
	}
	//if (drcyl*LengthUnits/Mpc < 1e3 && fabs(z)/Mpc < 4e3 && fabs(z)/Mpc > 2e3) {
	//	printf("FtotR %g, drcyl %g, x3 %g\n",FtotR, drcyl*LengthUnits/Mpc, z/Mpc);
	//}
	/*Roedigger claims that centrifugal force balances both grav and pressure force*/
	/* if (drcyl*LengthUnits*(FtotR+FdPdR) > 0.0) {
     if (denuse > densicm) {
     printf("v! FtotR = %g  dPdR = %g  FdPdR = %g  R = %g  z = %g  d = %g\n", FtotR, (Pressure2-Pressure)/(r2-drcyl*LengthUnits),FdPdR, drcyl*LengthUnits/Mpc, z/Mpc, denuse);
	 }
	 }*/
	
	if (*temperature < 0.0) {
		printf("temp = %g, P = %g, z = %g, zicm = %g, zicmf=%g, zsmall=%g, drcyl = %g\n", temperature, Pressure, z/Mpc, zicm/Mpc, zicmf, zsmall, drcyl*LengthUnits/Mpc);
	}
	if ((FtotRpot - FdPdR) > 0.0) {
		printf("FtotR = %g, FdPdR = %g, P = %g,P2 = %g, Picm = %g, dr = %g, drcyl = %g, z = %g\n", FtotRpot, FdPdR, Pressure, Pressure2, Picm, r2-drcyl*LengthUnits, drcyl*LengthUnits/Mpc, z/Mpc);
		FdPdR = 0.0;
	}
	
	vrot=sqrt(-drcyl*LengthUnits*(FtotRpot-FdPdR));
	
	//printf("vrot:  %g\n",vrot);
	if ((denuse == densicm*DensityUnits)) { // && (fabs(drcyl)*LengthUnits/Mpc <= 0.026)) {
		//printf("vrot original: %g\n", vrot);
		//vrot = vrot*exp(-fabs(z)/(0.030*Mpc));
		//printf("vrot new: %g, z: %g\n", vrot, z/Mpc);
		//if (fabs(z) > 0.060*Mpc) {
		vrot = 0.0;
		//}
	}
	
	//printf("vrot:  %g\n",vrot);
	//printf("FtotR = %g FdPdR = %g drcyl = %g r2 = %g\n", FtotR, FdPdR, drcyl, r2/LengthUnits);
	//if ((denuse == densicm) && (fabs(drcyl)*LengthUnits/Mpc) > 0.026) {
	//  vrot = 0.0;
	//  }
	
	
	return (vrot/VelocityUnits); //code units
	
}

Real av_den(Real cellwidth, Real z, Real xpos, Real ypos, Real zpos)
{
	// r and z are cylindrical coords.
	// routine to return the average gas density in a grid cell
	// Routine samples density in r plane of grid and averages
	// Assumes all input units are CGS!!
	
	int i,points;
	Real den,r1,nx,ny,nz;
	//	extern float AngularMomentumx, AngularMomentumxy, AngularMomentumxz;
	//	extern float DiskPositionx, DiskPositiony, DiskPositionz;
//	extern Real MBulge, rBulge, MSDisk, densDMConst, rDMConst,MgasScale;
//	extern Real gScaleHeightR, gScaleHeightz;
//	extern Real SDiskScaleHeightR, SDiskScaleHeightz;
//	extern Real LengthUnits,DensityUnits,TimeUnits;
//	extern Real GravConst, Mpc, pi,SolarMass;
//	extern Real densicm;
//	extern Real drcyl;
	Real avden;
	//points = 100;
	points = 10;
	//if (drcyl*LengthUnits/Mpc <= 0.01 && fabs(z) < 0.0001) {
	//printf("in density loop, r = %g, z = %g\n", drcyl*LengthUnits/Mpc, z/Mpc);
	//}
	if (fabs(drcyl)*LengthUnits <= 0.02*Mpc) 
	{
		den = MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2.0)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(z)/gScaleHeightz/Mpc);
		
		
		//avden = den/points;
		avden=den;
	}
	else 
	{
		//printf("using cosine density dist.  r = %g\n",drcyl*LengthUnits/Mpc);
		den = MgasScale*SolarMass/(2.0*pi*pow(gScaleHeightR*Mpc,2.0)*gScaleHeightz*Mpc)*0.25/cosh(drcyl*LengthUnits/gScaleHeightR/Mpc)/cosh(fabs(z)/gScaleHeightz/Mpc)*0.5*(1.0+cos(pi*(drcyl*LengthUnits-0.02*Mpc)/(0.006*Mpc)));
		
		
		//avden = den/points;
		 avden=den;
	}
	
	if (avden < densicm*DensityUnits || fabs(drcyl)*LengthUnits >= 0.026*Mpc) 
	{
		avden = densicm*DensityUnits;
	}
	//if (fabs(z) > 0.004*Mpc)
	//  {
	//    avden = densicm;
	//  }
	
	if (fabs(drcyl)*LengthUnits <= 0.026*Mpc) {
		if (avden < densicm*DensityUnits) {
			printf("avden = %g  drcyl = %g\n", avden, drcyl*LengthUnits/Mpc);
		}
	}
	//if (fabs(z)/Mpc < 0.0001 && drcyl*LengthUnits/Mpc <= 0.026) {
	//  printf("%g  %g\n",drcyl*LengthUnits/Mpc, avden);
	//} 
	//printf("avden %g\n", avden);
	return (avden/DensityUnits); //code units 
}


