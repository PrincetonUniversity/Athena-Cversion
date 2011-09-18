#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_3d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit CTU integrator of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *      B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically thin cooling,
 *   shearing box source terms, and the H-correction of Sanders et al.
 *     For adb hydro, requires (9*Cons1DS +  3*Real) = 48 3D arrays
 *     For adb mhd, requires   (9*Cons1DS + 10*Real) = 73 3D arrays
 *   The H-correction of Sanders et al. adds another 3 arrays.  
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport in three dimensions", JCP, 227, 4123 (2008)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 *   J.M. Stone et al., "Athena: A new code for astrophysical MHD", ApJS,
 *   178, 137 (2008)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_3d_radMHD()
 *   integrate_init_3d()
 *   integrate_destruct_3d()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS ***Ul_x1Face=NULL, ***Ur_x1Face=NULL;
static Cons1DS ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
static Cons1DS ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;


#ifdef CONS_GRAVITY
static Real ***x1Flux_grav=NULL;
static Real ***x2Flux_grav=NULL;
static Real ***x3Flux_grav=NULL;
static Real ***density_old=NULL;
Real dotphil, dotgxl;
#endif

/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static Real **remapEyiib=NULL, **remapEyoib=NULL;
#endif

/* The interface magnetic fields and emfs */
#ifdef RADIATION_MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;
static Real ****Source=NULL;
static Real ***Alpha = NULL;
static Real ****Beta = NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real ***dhalf = NULL, ***phalf=NULL;


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - the upwind CT method in GS05, for emf1
 *   integrate_emf2_corner() - the upwind CT method in GS05, for emf2
 *   integrate_emf3_corner() - the upwind CT method in GS05, for emf3
 *============================================================================*/

#ifdef RADIATION_MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif /* MHD */

void updatesource(GridS *pG);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/


void integrate_3d_radMHD(DomainS *pD)
{
	GridS *pG=(pD->Grid);
	Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;
	Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2, hdtodx3 = 0.5 * dtodx3;
	Real dx1 = pG->dx1, dx2 = pG->dx2, dx3 = pG->dx3;
	Real hdt = 0.5*pG->dt, dt = pG->dt;
  	int i,il,iu,is=pG->is, ie=pG->ie;
  	int j,jl,ju,js=pG->js, je=pG->je;
  	int k, kl, ku, ks=pG->ks, ke=pG->ke; 
	int m, n;

	/* For static gravitational potential */
	Real x1, x2, x3, phicl, phicr, phifc, phic, phir, phil;

	Real Bx = 0.0;
	Real M1h, M2h, M3h;
#ifdef RADIATION_MHD
  	Real MHD_src_By,MHD_src_Bz,mdb1,mdb2,mdb3;
  	Real db1,db2,db3,l1,l2,l3,B1,B2,B3,V1,V2,V3;
  	Real B1ch,B2ch,B3ch;
#endif

#if defined(RADIATION_MHD) || defined(SELF_GRAVITY)
  	Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;
#endif


#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#ifdef CONS_GRAVITY
  Real Tempswap;
#endif
#endif

	
/* Variables for shearing box */	
	
#ifdef SHEARING_BOX
	int my_iproc,my_jproc,my_kproc;
	Real M1n, dM2n; /* M1, dM2=(My+d*q*Omega_0*x) at time n */
	Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
	Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
	Real fact, qom, om_dt = Omega_0*pG->dt;
	fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
	qom = qshear*Omega_0;
	Real ShearSource[4];
/*	Real ShearingSource_M1, ShearingSource_M2, ShearingSource_M3, ShearingSource_E;
	Real Shearingguess_M1, Shearingguess_M2, Shearingguess_M3, Shearingguess_E;
 */
#endif /* SHEARING_BOX */

	Real temperature, velocity_x, velocity_y, velocity_z, velocity, pressure, density;
	Real Fr0x, Fr0y, Fr0z, diffTEr; /* co-moving flux, Fr-vf E_r/C , diffTEr = T^4 - Er*/
	Real Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE;
	/* The opacity function is:0-3 Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE */
	Real alpha, Propa_44, SEE, SErho, SEmx, SEmy, SEmz;
	Real dSigma[2*NOPACITY];
	Real Sigma[NOPACITY];
	/* dSigma for dSigma?/drho and dSigma?/dT */
	Cons1DS Usource;
	/* for source term */

	/* In case momentum becomes stiff */
	Real SFmx, SFmy, SFmz;


	Real Source_Inv[NVAR][NVAR], tempguess[NVAR], Uguess[NVAR], Source_guess[NVAR], Errort[NVAR], SourceFlux[NVAR];
	Real Psource;
	Real divFlux1[NVAR], divFlux2[NVAR], divFlux3[NVAR];

	/* SourceFlux is used to calculate flux due to non-stiff source terms, such as gravity */

	/* Initialize them to be zero */
	/* NVAR is the same as normal hydro code. Rad variables are not included here */
	for(i=0; i<NVAR; i++){
		SourceFlux[i] = 0.0;
		Source_guess[i] = 0.0;
		Errort[i] = 0.0;
		divFlux1[i] = 0.0;
		divFlux2[i] = 0.0;
		divFlux3[i] = 0.0;
		tempguess[i] = 0.0;
		for(j=0; j<NVAR; j++) {
			Source_Inv[i][j] = 0.0;			
			
			if(i==j) {
		 		Source_Inv[i][j] = 1.0;
		
			}
		}
	}

	for(i=0; i<4; i++)
		dSigma[i] = 0.0;


	il = is - 2;
  	iu = ie + 2;
  	jl = js - 2;
  	ju = je + 2;
  	kl = ks - 2;
  	ku = ke + 2;

/*=== Step 1: Backward Euler is done in the main function for the whole mesh ===*/
	
	/* First, calculate the source term */
	updatesource(pG);

	
/*=== STEP 2: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */
	for (k=kl; k<=ku; k++){
		for (j=jl; j<=ju; j++) {
    			for (i=is-nghost; i<=ie+nghost; i++) {
      			
			U1d[i].d  = pG->U[k][j][i].d;
			U1d[i].Mx = pG->U[k][j][i].M1;
			U1d[i].My = pG->U[k][j][i].M2;
			U1d[i].Mz = pG->U[k][j][i].M3;
			U1d[i].E  = pG->U[k][j][i].E;
#ifdef RADIATION_MHD
			U1d[i].By = pG->U[k][j][i].B2c;
			U1d[i].Bz = pG->U[k][j][i].B3c;
			Bxc[i] = pG->U[k][j][i].B1c;
			Bxi[i] = pG->B1i[k][j][i];
			B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */	
					
			U1d[i].Er  = pG->U[k][j][i].Er;
			U1d[i].Fr1  = pG->U[k][j][i].Fr1;
			U1d[i].Fr2  = pG->U[k][j][i].Fr2;
			U1d[i].Fr3  = pG->U[k][j][i].Fr3;
			U1d[i].Edd_11  = pG->U[k][j][i].Edd_11;
			U1d[i].Edd_21  = pG->U[k][j][i].Edd_21;
			U1d[i].Edd_22  = pG->U[k][j][i].Edd_22;
			U1d[i].Edd_31  = pG->U[k][j][i].Edd_31;
			U1d[i].Edd_32  = pG->U[k][j][i].Edd_32;
			U1d[i].Edd_33  = pG->U[k][j][i].Edd_33;
			for(m=0; m<NOPACITY;m++){
				U1d[i].Sigma[m] = pG->U[k][j][i].Sigma[m];

			}


			}


/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

			for (i=is-nghost; i<=ie+nghost; i++) {
      				W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
    			}

			lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,3);


/*------Step 2c: Add radiation source terms to the left and right state for 0.5*dt--------*/

		for (i=il+1; i<=iu; i++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			
			Propa_44 = Alpha[k][j][i-1];

			velocity_x = U1d[i-1].Mx / U1d[i-1].d;
			velocity_y = U1d[i-1].My / U1d[i-1].d;
			velocity_z = U1d[i-1].Mz / U1d[i-1].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i-1,j,k,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif

			/* The velocity here is the perturbed velocity */
			/* We do not need to add background shearing in */

			Psource =  (Gamma - 1.0) * Source[k][j][i-1][4]	- (Gamma - 1.0) * (velocity_x * Source[k][j][i-1][1] + velocity_y * Source[k][j][i-1][2] + velocity_z * Source[k][j][i-1][3]);



			Wl[i].Vx += dt * Source[k][j][i-1][1] * 0.5 * Beta[k][j][i-1][0] / U1d[i-1].d;
			Wl[i].Vy += dt * Source[k][j][i-1][2] * 0.5 * Beta[k][j][i-1][1] / U1d[i-1].d;
			Wl[i].Vz += dt * Source[k][j][i-1][3] * 0.5 * Beta[k][j][i-1][2] / U1d[i-1].d;
			Wl[i].P += dt * Propa_44 * Psource * 0.5;
			
			if(Wl[i].P < TINY_NUMBER) 
				Wl[i].P -= dt * Propa_44 * Psource * 0.5;	
	
			for(m=0; m<NOPACITY; m++)
				Wl[i].Sigma[m] = U1d[i-1].Sigma[m];
			
			Wl[i].Edd_11 = W[i-1].Edd_11;
			Wl[i].Edd_21 = W[i-1].Edd_21;
			Wl[i].Edd_22 = W[i-1].Edd_22;
			Wl[i].Edd_31 = W[i-1].Edd_31;
			Wl[i].Edd_32 = W[i-1].Edd_32;
			Wl[i].Edd_33 = W[i-1].Edd_33;

		/* For the right state */
	
	
			Propa_44 = alpha;

			velocity_x = U1d[i].Mx / U1d[i].d;
			velocity_y = U1d[i].My / U1d[i].d;
			velocity_z = U1d[i].Mz / U1d[i].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif
			
			Psource =  (Gamma - 1.0) * Source[k][j][i][4]	- (Gamma - 1.0) * (velocity_x * Source[k][j][i][1] + velocity_y * Source[k][j][i][2] + velocity_z * Source[k][j][i][3]);

			
			Wr[i].Vx += dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[i].d;
			Wr[i].Vy += dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[i].d;
			Wr[i].Vz += dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[i].d;
			Wr[i].P += dt * Propa_44 * Psource * 0.5;
			
			if(Wr[i].P < TINY_NUMBER)
				Wr[i].P -= dt * Propa_44 * Psource * 0.5;

			for(m=0; m<NOPACITY; m++)
				Wr[i].Sigma[m] = U1d[i].Sigma[m];

			Wr[i].Edd_11 = W[i].Edd_11;
			Wr[i].Edd_21 = W[i].Edd_21;
			Wr[i].Edd_22 = W[i].Edd_22;
			Wr[i].Edd_31 = W[i].Edd_31;
			Wr[i].Edd_32 = W[i].Edd_32;
			Wr[i].Edd_33 = W[i].Edd_33;		
			
		}

/*---------For radiation MHD, source term due to magnetic field is also added---------*/
		
#ifdef RADIATION_MHD
      		for (i=il+1; i<=iu; i++) {
/* Source terms for left states in zone i-1 */

	        db1 = (    pG->B1i[k  ][j  ][i  ] -     pG->B1i[k][j][i-1])*dx1i;
       		db2 = (    pG->B2i[k  ][j+1][i-1] -     pG->B2i[k][j][i-1])*dx2i;
        	db3 = (    pG->B3i[k+1][j  ][i-1] -     pG->B3i[k][j][i-1])*dx3i;

		if(db1 >= 0.0){
	  		l3 = db1 < -db3 ? db1 : -db3;
	  		l3 = l3 > 0.0 ? l3 : 0.0;

	  		l2 = db1 < -db2 ? db1 : -db2;
	  		l2 = l2 > 0.0 ? l2 : 0.0;
		}
		else{
	  		l3 = db1 > -db3 ? db1 : -db3;
	  		l3 = l3 < 0.0 ? l3 : 0.0;

	  		l2 = db1 > -db2 ? db1 : -db2;
	  		l2 = l2 < 0.0 ? l2 : 0.0;
		}

	        MHD_src_By = (pG->U[k][j][i-1].M2/pG->U[k][j][i-1].d)*l2;
	        MHD_src_Bz = (pG->U[k][j][i-1].M3/pG->U[k][j][i-1].d)*l3;

	        Wl[i].By += hdt*MHD_src_By;
	        Wl[i].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone i */

	        db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        	db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        	db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;

        	if(db1 >= 0.0){
          		l3 = db1 < -db3 ? db1 : -db3;
          		l3 = l3 > 0.0 ? l3 : 0.0;

          		l2 = db1 < -db2 ? db1 : -db2;
          		l2 = l2 > 0.0 ? l2 : 0.0;
        	}
        	else{
          		l3 = db1 > -db3 ? db1 : -db3;
          		l3 = l3 < 0.0 ? l3 : 0.0;

          		l2 = db1 > -db2 ? db1 : -db2;
          		l2 = l2 < 0.0 ? l2 : 0.0;
        	}

        	MHD_src_By = (pG->U[k][j][i].M2/pG->U[k][j][i].d)*l2;
        	MHD_src_Bz = (pG->U[k][j][i].M3/pG->U[k][j][i].d)*l3;

        	Wr[i].By += hdt*MHD_src_By;
        	Wr[i].Bz += hdt*MHD_src_Bz;
      		}
#endif

/*-------Add source term due to static gravitational potential for left and right state--------*/
	
      		if (StaticGravPot != NULL){
        		for (i=il+1; i<=iu; i++) {
         			 cc_pos(pG,i,j,k,&x1,&x2,&x3);

	        	phicr = (*StaticGravPot)( x1             ,x2,x3);
          		phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
          		phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          		Wl[i].Vx -= dtodx1*(phifc - phicl);
          		Wr[i].Vx -= dtodx1*(phicr - phifc);

        		}
      		}

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
      for (i=il+1; i<=iu; i++) {
        Wl[i].Vx -= hdtodx1*(pG->Phi[k][j][i] - pG->Phi[k][j][i-1]);
        Wr[i].Vx -= hdtodx1*(pG->Phi[k][j][i] - pG->Phi[k][j][i-1]);
      }
#endif
			
/* Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
 * starting with tidal gravity terms added through the ShearingBoxPot
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
 */		
				
#ifdef SHEARING_BOX
		if (ShearingBoxPot != NULL){
			for (i=il+1; i<=iu; i++) {
						cc_pos(pG,i,j,k,&x1,&x2,&x3);
						phicr = (*ShearingBoxPot)( x1             ,x2,x3);
						phicl = (*ShearingBoxPot)((x1-    pG->dx1),x2,x3);
						phifc = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
						
						Wl[i].Vx -= dtodx1*(phifc - phicl);
						Wr[i].Vx -= dtodx1*(phicr - phifc);
					}
				}
			
		for (i=il+1; i<=iu; i++) {
				Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vy;
				Wr[i].Vx += pG->dt*Omega_0*W[i].Vy;
#ifdef FARGO
				Wl[i].Vy += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
				Wr[i].Vy += hdt*(qshear-2.)*Omega_0*W[i].Vx;
#else
				Wl[i].Vy -= pG->dt*Omega_0*W[i-1].Vx;
				Wr[i].Vy -= pG->dt*Omega_0*W[i].Vx;
#endif
		}
			
#endif /* SHEARING_BOX */				


/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
	      	for (i=il+1; i<=iu; i++) {
        		Ul_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);
        		Ur_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

#ifdef RADIATION_MHD
        		Bx = B1_x1Face[k][j][i];
#endif
			/* take the time step with x1Flux */
			x1Flux[k][j][i].d = dt;	
				x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);
      		}
		

		} /* End j */

	}/* End k */	


/*=== STEP 3: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/
/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 * NOTICE that momentum is rotated but flux and Eddington tensor are not *
 */
	for (k=kl; k<=ku; k++){
		for (i=il; i<=iu; i++) {
    			for (j=js-nghost; j<=je+nghost; j++) {
      				U1d[j].d  = pG->U[k][j][i].d;
      				U1d[j].Mx = pG->U[k][j][i].M2;
      				U1d[j].My = pG->U[k][j][i].M3;
      				U1d[j].Mz = pG->U[k][j][i].M1;
      				U1d[j].E  = pG->U[k][j][i].E;
				U1d[j].Er  = pG->U[k][j][i].Er;
    				U1d[j].Fr1  = pG->U[k][j][i].Fr1;
				U1d[j].Fr2  = pG->U[k][j][i].Fr2;
    				U1d[j].Fr3  = pG->U[k][j][i].Fr3;
				U1d[j].Edd_11  = pG->U[k][j][i].Edd_11;
				U1d[j].Edd_21  = pG->U[k][j][i].Edd_21;
				U1d[j].Edd_22  = pG->U[k][j][i].Edd_22;
				U1d[j].Edd_31  = pG->U[k][j][i].Edd_31;
				U1d[j].Edd_32  = pG->U[k][j][i].Edd_32;
				U1d[j].Edd_33  = pG->U[k][j][i].Edd_33;
				for(m=0; m<NOPACITY;m++)
					U1d[j].Sigma[m] = pG->U[k][j][i].Sigma[m];
		
				
#ifdef RADIATION_MHD
      				U1d[j].By = pG->U[k][j][i].B3c;
      				U1d[j].Bz = pG->U[k][j][i].B1c;
      				Bxc[j] = pG->U[k][j][i].B2c;
      				Bxi[j] = pG->B2i[k][j][i];
      				B2_x2Face[k][j][i] = pG->B2i[k][j][i];
#endif /* MHD */
			}


		for (j=js-nghost; j<=je+nghost; j++) {
        		W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx2,jl+1,ju-1,Wl,Wr,3);

/*------Step 3b: Add radiation source terms to the left and right state for 0.5*dt--------*/


		for (j=jl+1; j<=ju; j++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			Propa_44 = Alpha[k][j-1][i];

			velocity_x = U1d[j-1].Mz / U1d[j-1].d;
			velocity_y = U1d[j-1].Mx / U1d[j-1].d;
			velocity_z = U1d[j-1].My / U1d[j-1].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j-1,k,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif

			Psource =  (Gamma - 1.0) * Source[k][j-1][i][4]	- (Gamma - 1.0) * (velocity_x * Source[k][j-1][i][1] + velocity_y * Source[k][j-1][i][2] + velocity_z * Source[k][j-1][i][3]);

			
			Wl[j].Vx += dt * Source[k][j-1][i][2] * 0.5 * Beta[k][j-1][i][1] / U1d[j-1].d;
			Wl[j].Vy += dt * Source[k][j-1][i][3] * 0.5 * Beta[k][j-1][i][2] / U1d[j-1].d;
			Wl[j].Vz += dt * Source[k][j-1][i][1] * 0.5 * Beta[k][j-1][i][0] / U1d[j-1].d;
			Wl[j].P += dt * Propa_44 * Psource * 0.5;
			
			if(Wl[j].P < TINY_NUMBER)
				Wl[j].P -= dt * Propa_44 * Psource * 0.5;

			
			for(m=0; m<NOPACITY;m++){
				Wl[j].Sigma[m] = U1d[j-1].Sigma[m];
			}

			Wl[j].Edd_11 = W[j-1].Edd_11;
			Wl[j].Edd_21 = W[j-1].Edd_21;
			Wl[j].Edd_22 = W[j-1].Edd_22;
			Wl[j].Edd_31 = W[j-1].Edd_31;
			Wl[j].Edd_32 = W[j-1].Edd_32;
			Wl[j].Edd_33 = W[j-1].Edd_33;



		/* For the right state */
	
	
			Propa_44 = Alpha[k][j][i];

			velocity_x = U1d[j].Mz / U1d[j].d;
			velocity_y = U1d[j].Mx / U1d[j].d;
			velocity_z = U1d[j].My / U1d[j].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif

			Psource =  (Gamma - 1.0) * Source[k][j][i][4]	- (Gamma - 1.0) * (velocity_x * Source[k][j][i][1] + velocity_y * Source[k][j][i][2] + velocity_z * Source[k][j][i][3]);

			
			Wr[j].Vx += dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[j].d;
			Wr[j].Vy += dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[j].d;
			Wr[j].Vz += dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[j].d;
			Wr[j].P += dt * Propa_44 * Psource * 0.5;

			if(Wr[j].P < TINY_NUMBER)
				Wr[j].P -= dt * Propa_44 * Psource * 0.5;
			
			for(m=0; m<NOPACITY; m++){
				Wr[j].Sigma[m] = U1d[j].Sigma[m];

			}


			Wr[j].Edd_11 = W[j].Edd_11;
			Wr[j].Edd_21 = W[j].Edd_21;
			Wr[j].Edd_22 = W[j].Edd_22;
			Wr[j].Edd_31 = W[j].Edd_31;
			Wr[j].Edd_32 = W[j].Edd_32;
			Wr[j].Edd_33 = W[j].Edd_33;		
			
		}


#ifdef RADIATION_MHD

	      	for (j=jl+1; j<=ju; j++) {
/* Source terms for left states in zone j-1 */
        		db1 = (    pG->B1i[k  ][j-1][i+1] -     pG->B1i[k][j-1][i])*dx1i;
        		db2 = (    pG->B2i[k  ][j  ][i  ] -     pG->B2i[k][j-1][i])*dx2i;
        		db3 = (    pG->B3i[k+1][j-1][i  ] -     pG->B3i[k][j-1][i])*dx3i;

			if(db2 >= 0.0){
	  			l1 = db2 < -db1 ? db2 : -db1;
	  			l1 = l1 > 0.0 ? l1 : 0.0;

	  			l3 = db2 < -db3 ? db2 : -db3;
	  			l3 = l3 > 0.0 ? l3 : 0.0;
			}
			else{
				l1 = db2 > -db1 ? db2 : -db1;
	  			l1 = l1 < 0.0 ? l1 : 0.0;

				l3 = db2 > -db3 ? db2 : -db3;
	  			l3 = l3 < 0.0 ? l3 : 0.0;
			}

			MHD_src_By = (pG->U[k][j-1][i].M3/pG->U[k][j-1][i].d)*l3;
			MHD_src_Bz = (pG->U[k][j-1][i].M1/pG->U[k][j-1][i].d)*l1;

       	 		Wl[j].By += hdt*MHD_src_By;
        		Wl[j].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone j */
        		db1 = (    pG->B1i[k  ][j  ][i+1] -    pG->B1i[k][j][i])*dx1i;
        		db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        		db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;

        		if(db2 >= 0.0){
          			l1 = db2 < -db1 ? db2 : -db1;
          			l1 = l1 > 0.0 ? l1 : 0.0;

          			l3 = db2 < -db3 ? db2 : -db3;
          			l3 = l3 > 0.0 ? l3 : 0.0;
        		}
        		else{
          			l1 = db2 > -db1 ? db2 : -db1;
          			l1 = l1 < 0.0 ? l1 : 0.0;

          			l3 = db2 > -db3 ? db2 : -db3;
          			l3 = l3 < 0.0 ? l3 : 0.0;
        		}

        		MHD_src_By = (pG->U[k][j][i].M3/pG->U[k][j][i].d)*l3;
        		MHD_src_Bz = (pG->U[k][j][i].M1/pG->U[k][j][i].d)*l1;

        		Wr[j].By += hdt*MHD_src_By;
        		Wr[j].Bz += hdt*MHD_src_Bz;
      		}	
#endif


		
      		if (StaticGravPot != NULL){
        		for (j=jl+1; j<=ju; j++) {
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);

          			phicr = (*StaticGravPot)(x1, x2             ,x3);
         	 		phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          			phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          		 	Wl[j].Vx -= dtodx2*(phifc - phicl);
          			Wr[j].Vx -= dtodx2*(phicr - phifc);
        	
			}
      		}


 /* Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
      		for (j=jl+1; j<=ju; j++) {
        		Wl[j].Vx -= hdtodx2*(pG->Phi[k][j][i] - pG->Phi[k][j-1][i]);
        		Wr[j].Vx -= hdtodx2*(pG->Phi[k][j][i] - pG->Phi[k][j-1][i]);
      		}
#endif



/*-----Calculate flux along x2 direction ---------*/

		for (j=jl+1; j<=ju; j++) {
        		Ul_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);
        		Ur_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);

#ifdef RADIATION_MHD
        		Bx = B2_x2Face[k][j][i];
#endif
			/* Take the time step with x2Flux[][][].d */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[j],Wr[j],Bx,&x2Flux[k][j][i]);
      		}

		}/* END i */
	}/* end k */	




/*=== STEP 4: Compute L/R x3-interface states and 1D x3-Fluxes ===============*/
/*--- Step 4a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M3, M1, M2, E, B1c, B2c, s[n])
 * NOTICE that momentum is rotated but flux and Eddington tensor are not *
 */
	for (j=jl; j<=ju; j++){
		for (i=il; i<=iu; i++) {
    			for (k=ks-nghost; k<=ke+nghost; k++) {
      				U1d[k].d  = pG->U[k][j][i].d;
      				U1d[k].Mx = pG->U[k][j][i].M3;
      				U1d[k].My = pG->U[k][j][i].M1;
      				U1d[k].Mz = pG->U[k][j][i].M2;
      				U1d[k].E  = pG->U[k][j][i].E;
				U1d[k].Er  = pG->U[k][j][i].Er;
    				U1d[k].Fr1  = pG->U[k][j][i].Fr1;
				U1d[k].Fr2  = pG->U[k][j][i].Fr2;
    				U1d[k].Fr3  = pG->U[k][j][i].Fr3;
				U1d[k].Edd_11  = pG->U[k][j][i].Edd_11;
				U1d[k].Edd_21  = pG->U[k][j][i].Edd_21;
				U1d[k].Edd_22  = pG->U[k][j][i].Edd_22;
				U1d[k].Edd_31  = pG->U[k][j][i].Edd_31;
				U1d[k].Edd_32  = pG->U[k][j][i].Edd_32;
				U1d[k].Edd_33  = pG->U[k][j][i].Edd_33;

				for(m=0; m<NOPACITY;m++){
					U1d[k].Sigma[m] = pG->U[k][j][i].Sigma[m];
				}

#ifdef RADIATION_MHD
      				U1d[k].By = pG->U[k][j][i].B1c;
      				U1d[k].Bz = pG->U[k][j][i].B2c;
      				Bxc[k] = pG->U[k][j][i].B3c;
      				Bxi[k] = pG->B3i[k][j][i];
      				B3_x3Face[k][j][i] = pG->B3i[k][j][i];
#endif /* MHD */
			}


		for (k=ks-nghost; k<=ke+nghost; k++) {
        		W[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx3,kl+1,ku-1,Wl,Wr,3);

/*------Step 3b: Add radiation source terms to the left and right state for 0.5*dt--------*/


		for (k=kl+1; k<=ku; k++) {
		/* NOTICE that we use primitive variables as left and right states to calculate flux */
		/* For left state */
			Propa_44 = Alpha[k-1][j][i];

			velocity_x = U1d[k-1].My / U1d[k-1].d;
			velocity_y = U1d[k-1].Mz / U1d[k-1].d;
			velocity_z = U1d[k-1].Mx / U1d[k-1].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j,k-1,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif

			Psource =  (Gamma - 1.0) * Source[k-1][j][i][4]	- (Gamma - 1.0) * (velocity_x * Source[k-1][j][i][1] + velocity_y * Source[k-1][j][i][2] + velocity_z * Source[k-1][j][i][3]);

			
			Wl[k].Vx += dt * Source[k-1][j][i][3] * 0.5 * Beta[k-1][j][i][2] / U1d[k-1].d;
			Wl[k].Vy += dt * Source[k-1][j][i][1] * 0.5 * Beta[k-1][j][i][0] / U1d[k-1].d;
			Wl[k].Vz += dt * Source[k-1][j][i][2] * 0.5 * Beta[k-1][j][i][1] / U1d[k-1].d;
			Wl[k].P += dt * Propa_44 * Psource * 0.5;
			
			if(Wl[k].P < TINY_NUMBER)
				Wl[k].P -= dt * Propa_44 * Psource * 0.5;
			
	
			for(m=0;m<NOPACITY;m++){
				Wl[k].Sigma[m] = U1d[k-1].Sigma[m];
			}

			Wl[k].Edd_11 = W[k-1].Edd_11;
			Wl[k].Edd_21 = W[k-1].Edd_21;
			Wl[k].Edd_22 = W[k-1].Edd_22;
			Wl[k].Edd_31 = W[k-1].Edd_31;
			Wl[k].Edd_32 = W[k-1].Edd_32;
			Wl[k].Edd_33 = W[k-1].Edd_33;

		/* For the right state */
	
	
			
			Propa_44 = Alpha[k][j][i];

			velocity_x = U1d[k].My / U1d[k].d;
			velocity_y = U1d[k].Mz / U1d[k].d;
			velocity_z = U1d[k].Mx / U1d[k].d;

			
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			velocity_y -= qom * x1;		
#endif

			Psource =  (Gamma - 1.0) * Source[k][j][i][4]	- (Gamma - 1.0) * (velocity_x * Source[k][j][i][1] + velocity_y * Source[k][j][i][2] + velocity_z * Source[k][j][i][3]);

			
			Wr[k].Vx += dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[k].d;
			Wr[k].Vy += dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[k].d;
			Wr[k].Vz += dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[k].d;
			Wr[k].P += dt * Propa_44 * Psource * 0.5;
			
			if(Wr[k].P < TINY_NUMBER)
				Wr[k].P -= dt * Propa_44 * Psource * 0.5;

			for(m=0; m<NOPACITY;m++){
				Wr[k].Sigma[m] = U1d[i].Sigma[m];
			}
	
			Wr[k].Edd_11 = W[k].Edd_11;
			Wr[k].Edd_21 = W[k].Edd_21;
			Wr[k].Edd_22 = W[k].Edd_22;
			Wr[k].Edd_31 = W[k].Edd_31;
			Wr[k].Edd_32 = W[k].Edd_32;
			Wr[k].Edd_33 = W[k].Edd_33;		
			
		}


#ifdef RADIATION_MHD

 		for (k=kl+1; k<=ku; k++) {
/* Source terms for left states in zone k-1 */
        		db1 = (    pG->B1i[k-1][j  ][i+1] -     pG->B1i[k-1][j][i])*dx1i;
        		db2 = (    pG->B2i[k-1][j+1][i  ] -     pG->B2i[k-1][j][i])*dx2i;
        		db3 = (    pG->B3i[k  ][j  ][i  ] -     pG->B3i[k-1][j][i])*dx3i;

			if(db3 >= 0.0){
	  			l1 = db3 < -db1 ? db3 : -db1;
	  			l1 = l1 > 0.0 ? l1 : 0.0;

	  			l2 = db3 < -db2 ? db3 : -db2;
	  			l2 = l2 > 0.0 ? l2 : 0.0;
			}
			else{
	  			l1 = db3 > -db1 ? db3 : -db1;
	  			l1 = l1 < 0.0 ? l1 : 0.0;

	  			l2 = db3 > -db2 ? db3 : -db2;
	  			l2 = l2 < 0.0 ? l2 : 0.0;
			}

				MHD_src_By = (pG->U[k-1][j][i].M1/pG->U[k-1][j][i].d)*l1;
				MHD_src_Bz = (pG->U[k-1][j][i].M2/pG->U[k-1][j][i].d)*l2;

        			Wl[k].By += hdt*MHD_src_By;
        			Wl[k].Bz += hdt*MHD_src_Bz;

/* Source terms for right states in zone k */
        			db1 = (    pG->B1i[k][j][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k][j+1][i] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j][i] -     pG->B3i[k][j][i])*dx3i;

			if(db3 >= 0.0){
         			 l1 = db3 < -db1 ? db3 : -db1;
          			l1 = l1 > 0.0 ? l1 : 0.0;

          			l2 = db3 < -db2 ? db3 : -db2;
          			l2 = l2 > 0.0 ? l2 : 0.0;
        		}
        		else{
          			l1 = db3 > -db1 ? db3 : -db1;
          			l1 = l1 < 0.0 ? l1 : 0.0;

          			l2 = db3 > -db2 ? db3 : -db2;
          			l2 = l2 < 0.0 ? l2 : 0.0;
        		}

        		MHD_src_By = (pG->U[k][j][i].M1/pG->U[k][j][i].d)*l1;
       		 	MHD_src_Bz = (pG->U[k][j][i].M2/pG->U[k][j][i].d)*l2;

       		 	Wr[k].By += hdt*MHD_src_By;
        		Wr[k].Bz += hdt*MHD_src_Bz;
      }			
#endif


		
      		if (StaticGravPot != NULL){
        		for (k=kl+1; k<=ku; j++) {
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);

          			phicr = (*StaticGravPot)(x1,x2, x3             );
          			phicl = (*StaticGravPot)(x1,x2,(x3-    pG->dx3));
          			phifc = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          			Wl[k].Vx -= dtodx3*(phifc - phicl);
          			Wr[k].Vx -= dtodx3*(phicr - phifc);
        	
			}
      		}

#ifdef SELF_GRAVITY
      		for (k=kl+1; k<=ku; k++) {
        		Wl[k].Vx -= hdtodx3*(pG->Phi[k][j][i] - pG->Phi[k-1][j][i]);
        		Wr[k].Vx -= hdtodx3*(pG->Phi[k][j][i] - pG->Phi[k-1][j][i]);
      		}
#endif



/*-----Calculate flux along x3 direction ---------*/

		for (k=kl+1; k<=ku; k++) {
        		Ul_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wl[k],&Bxi[k]);
        		Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[k],&Bxi[k]);

#ifdef RADIATION_MHD
        		Bx = B3_x3Face[k][j][i];
#endif
			/* Take the time step with x3Flux[][][].d */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[k],Wr[k],Bx,&x3Flux[k][j][i]);
      		}

		}/* END i */
	}/* end j */	




/*=== STEP 5:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 5a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */

#ifdef RADIATION_MHD
/* emf1 */
  		for (k=kl; k<=ku; k++) {
    			for (j=jl; j<=ju; j++) {
      				for (i=il; i<=iu; i++) {
        				emf1_cc[k][j][i] = (pG->U[k][j][i].B2c*pG->U[k][j][i].M3 -
			    		pG->U[k][j][i].B3c*pG->U[k][j][i].M2)
                              		/pG->U[k][j][i].d;
        				emf2_cc[k][j][i] = (pG->U[k][j][i].B3c*pG->U[k][j][i].M1 -
			    		pG->U[k][j][i].B1c*pG->U[k][j][i].M3)
                              		/pG->U[k][j][i].d;
        				emf3_cc[k][j][i] = (pG->U[k][j][i].B1c*pG->U[k][j][i].M2 -
			    		pG->U[k][j][i].B2c*pG->U[k][j][i].M1)
                              		/pG->U[k][j][i].d;
      				}
    			}
  		}
  		
		integrate_emf1_corner(pG);
  		integrate_emf2_corner(pG);
  		integrate_emf3_corner(pG);

/*--- Step 5b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */

  		for (k=kl+1; k<=ku-1; k++) {
    			for (j=jl+1; j<=ju-1; j++) {
      				for (i=il+1; i<=iu-1; i++) {

        			B1_x1Face[k][j][i] += hdtodx3 * (emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                              		hdtodx2 * (emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        			B2_x2Face[k][j][i] += hdtodx1 * (emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                              		hdtodx3 * (emf1[k+1][j  ][i  ] - emf1[k][j][i]);

        			B3_x3Face[k][j][i] += hdtodx2 * (    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                              		hdtodx1 * (emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      			}

      			B1_x1Face[k][j][iu] += hdtodx3 * (emf2[k+1][j  ][iu]-emf2[k][j][iu]) -
                               hdtodx2 * (emf3[k  ][j+1][iu]-emf3[k][j][iu]);
    			}
    			
			for (i=il+1; i<=iu-1; i++) {
      				B2_x2Face[k][ju][i] += hdtodx1 * (emf3[k  ][ju][i+1]-emf3[k][ju][i]) -
                               		hdtodx3 * (emf1[k+1][ju][i  ]-emf1[k][ju][i]);
    				}
  			}

  			for (j=jl+1; j<=ju-1; j++) {
    				for (i=il+1; i<=iu-1; i++) {

      					B3_x3Face[ku][j][i] += hdtodx2 * (    emf1[ku][j+1][i  ] -     emf1[ku][j][i]) -
                             			hdtodx1 * (emf2[ku][j  ][i+1] - emf2[ku][j][i]);
    				}
  			}
#endif /* MHD */

/*=== STEP 6: Correct x1-interface states with transverse flux gradients =====*/
/*--- Step 6a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes and x3 fluxes  computed in previous steps.
 * Since the fluxes come from an -sweep, (x,y,z) on RHS -> (z,x,y) for y and (x,y,z) on RHS ->(z,x,y) on LHS
 * The order is rotated-------------------------------------------
 */

  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

        		Ul_x1Face[k][j][i].d -=hdtodx2*(x2Flux[k][j+1][i-1].d -x2Flux[k][j][i-1].d );
        		Ul_x1Face[k][j][i].Mx-=hdtodx2*(x2Flux[k][j+1][i-1].Mz-x2Flux[k][j][i-1].Mz);
        		Ul_x1Face[k][j][i].My-=hdtodx2*(x2Flux[k][j+1][i-1].Mx-x2Flux[k][j][i-1].Mx);
        		Ul_x1Face[k][j][i].Mz-=hdtodx2*(x2Flux[k][j+1][i-1].My-x2Flux[k][j][i-1].My);

        		Ul_x1Face[k][j][i].E -=hdtodx2*(x2Flux[k][j+1][i-1].E -x2Flux[k][j][i-1].E );

#ifdef RADIATION_MHD
/* Update B3 */
			Ul_x1Face[k][j][i].Bz+=hdtodx2*0.5*((emf1[k  ][j+1][i-1] - emf1[k  ][j][i-1]) 
					+ (emf1[k+1][j+1][i-1] - emf1[k+1][j][i-1]));
#endif

       			Ur_x1Face[k][j][i].d -=hdtodx2*(x2Flux[k][j+1][i  ].d -x2Flux[k][j][i  ].d );
        		Ur_x1Face[k][j][i].Mx-=hdtodx2*(x2Flux[k][j+1][i  ].Mz-x2Flux[k][j][i  ].Mz);
        		Ur_x1Face[k][j][i].My-=hdtodx2*(x2Flux[k][j+1][i  ].Mx-x2Flux[k][j][i  ].Mx);
        		Ur_x1Face[k][j][i].Mz-=hdtodx2*(x2Flux[k][j+1][i  ].My-x2Flux[k][j][i  ].My);
        		Ur_x1Face[k][j][i].E -=hdtodx2*(x2Flux[k][j+1][i  ].E -x2Flux[k][j][i  ].E );
#ifdef RADIATION_MHD
/* Update B3 */
			Ur_x1Face[k][j][i].Bz+=hdtodx2* 0.5 * ((emf1[k  ][j+1][i] - emf1[k  ][j][i]) + (emf1[k+1][j+1][i] - emf1[k+1][j][i]));
#endif


/*--- Step 6b ------------------------------------------------------------------
 * Correct x1-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */

        		Ul_x1Face[k][j][i].d  -= hdtodx3 * (x3Flux[k+1][j][i-1].d -x3Flux[k][j][i-1].d );
        		Ul_x1Face[k][j][i].Mx -= hdtodx3 * (x3Flux[k+1][j][i-1].My-x3Flux[k][j][i-1].My);
        		Ul_x1Face[k][j][i].My -= hdtodx3 * (x3Flux[k+1][j][i-1].Mz-x3Flux[k][j][i-1].Mz);
        		Ul_x1Face[k][j][i].Mz -= hdtodx3 * (x3Flux[k+1][j][i-1].Mx-x3Flux[k][j][i-1].Mx);

        		Ul_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k+1][j][i-1].E -x3Flux[k][j][i-1].E );
#ifdef RADIATION_MHD
/* Update B2 */
			Ul_x1Face[k][j][i].By -= hdtodx3 * 0.5* ((emf1[k+1][j  ][i-1] - emf1[k][j  ][i-1]) +
	   					(emf1[k+1][j+1][i-1] - emf1[k][j+1][i-1]));
#endif

        		Ur_x1Face[k][j][i].d  -= hdtodx3 * (x3Flux[k+1][j][i  ].d -x3Flux[k][j][i  ].d );
        		Ur_x1Face[k][j][i].Mx -= hdtodx3 * (x3Flux[k+1][j][i  ].My-x3Flux[k][j][i  ].My);
        		Ur_x1Face[k][j][i].My -= hdtodx3 * (x3Flux[k+1][j][i  ].Mz-x3Flux[k][j][i  ].Mz);
        		Ur_x1Face[k][j][i].Mz -= hdtodx3 * (x3Flux[k+1][j][i  ].Mx-x3Flux[k][j][i  ].Mx);
        		Ur_x1Face[k][j][i].E  -= hdtodx3 * (x3Flux[k+1][j][i  ].E -x3Flux[k][j][i  ].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ur_x1Face[k][j][i].By -= hdtodx3 * 0.5*((emf1[k+1][j  ][i] - emf1[k][j  ][i]) +
	   					(emf1[k+1][j+1][i] - emf1[k][j+1][i]));
#endif

      }
    }
  }


/*--- Step 6b ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2- and x3-flux-gradients to the
 * conservative variables on the x1Face.  Limiting is used as in GS (2007)
 */

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

        			db1 = (    pG->B1i[k  ][j  ][i  ] -     pG->B1i[k][j][i-1])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i-1] -     pG->B2i[k][j][i-1])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i-1] -     pG->B3i[k][j][i-1])*dx3i;
        			B1 = pG->U[k][j][i-1].B1c;
        			B2 = pG->U[k][j][i-1].B2c;
        			B3 = pG->U[k][j][i-1].B3c;
        			V2 = pG->U[k][j][i-1].M2/pG->U[k][j][i-1].d;
        			V3 = pG->U[k][j][i-1].M3/pG->U[k][j][i-1].d;

/* Calculate mdb2 = min_mod(-db1,db2) */
        			if(db1 > 0.0 && db2 < 0.0){
          				mdb2 = db2 > -db1 ? db2 : -db1;
        			}
        			else if(db1 < 0.0 && db2 > 0.0){
          				mdb2 = db2 < -db1 ? db2 : -db1;
        			}
        			else mdb2 = 0.0;

/* Calculate mdb3 = min_mod(-db1,db3) */
        			if(db1 > 0.0 && db3 < 0.0){
          				mdb3 = db3 > -db1 ? db3 : -db1;
        			}
        			else if(db1 < 0.0 && db3 > 0.0){
          				mdb3 = db3 < -db1 ? db3 : -db1;
        			}
        			else mdb3 = 0.0;

        			Ul_x1Face[k][j][i].Mx += hdt*B1*db1;
        			Ul_x1Face[k][j][i].My += hdt*B2*db1;
        			Ul_x1Face[k][j][i].Mz += hdt*B3*db1;
        			Ul_x1Face[k][j][i].By += hdt*V2*(-mdb3);
        			Ul_x1Face[k][j][i].Bz += hdt*V3*(-mdb2);

        			Ul_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
        			V2 = pG->U[k][j][i].M2/pG->U[k][j][i].d;
        			V3 = pG->U[k][j][i].M3/pG->U[k][j][i].d;

/* Calculate mdb2 = min_mod(-db1,db2) */
        			if(db1 > 0.0 && db2 < 0.0){
          				mdb2 = db2 > -db1 ? db2 : -db1;
        			}
        			else if(db1 < 0.0 && db2 > 0.0){
          				mdb2 = db2 < -db1 ? db2 : -db1;
        			}
        			else mdb2 = 0.0;

/* Calculate mdb3 = min_mod(-db1,db3) */
        			if(db1 > 0.0 && db3 < 0.0){
          				mdb3 = db3 > -db1 ? db3 : -db1;
        			}
        			else if(db1 < 0.0 && db3 > 0.0){
         	 			mdb3 = db3 < -db1 ? db3 : -db1;
        			}
        			else mdb3 = 0.0;

        			Ur_x1Face[k][j][i].Mx += hdt*B1*db1;
        			Ur_x1Face[k][j][i].My += hdt*B2*db1;
        			Ur_x1Face[k][j][i].Mz += hdt*B3*db1;
        			Ur_x1Face[k][j][i].By += hdt*V2*(-mdb3);
        			Ur_x1Face[k][j][i].Bz += hdt*V3*(-mdb2);

        			Ur_x1Face[k][j][i].E  += hdt*(B2*V2*(-mdb3) + B3*V3*(-mdb2) );

      }
    }
  }
#endif /* RADIATION MHD */


/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x2-Flux
 * and x3-Flux gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  if (StaticGravPot != NULL){
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {
        			cc_pos(pG,i,j,k,&x1,&x2,&x3);
        			
				phic = (*StaticGravPot)(x1, x2             ,x3);
        			phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        			phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

/* correct right states; x2 and x3 gradients */

        			Ur_x1Face[k][j][i].My -= hdtodx2 * (phir-phil)*pG->U[k][j][i].d;

        			Ur_x1Face[k][j][i].E -= hdtodx2 * (x2Flux[k][j  ][i  ].d*(phic - phil)
                                  + x2Flux[k][j+1][i  ].d*(phir - phic));


        			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        			phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));
        
        			Ur_x1Face[k][j][i].Mz -= hdtodx3 * (phir-phil)*pG->U[k][j][i].d;

        			Ur_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j][i  ].d*(phic - phil)
                                  	+ x3Flux[k+1][j][i  ].d*(phir - phic));


/* correct left states; x2 and x3 gradients */
        			phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        			phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        			phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);


        			Ul_x1Face[k][j][i].My -= hdtodx2 * (phir-phil)*pG->U[k][j][i-1].d;

        			Ul_x1Face[k][j][i].E -= hdtodx2 * (x2Flux[k][j  ][i-1].d*(phic - phil)
                                  	+ x2Flux[k][j+1][i-1].d*(phir - phic));


        			phir = (*StaticGravPot)((x1-pG->dx1),x2,(x3+0.5*pG->dx3));
        			phil = (*StaticGravPot)((x1-pG->dx1),x2,(x3-0.5*pG->dx3));
        
        			Ul_x1Face[k][j][i].Mz -= hdtodx3 * (phir-phil)*pG->U[k][j][i-1].d;

        			Ul_x1Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j][i-1].d*(phic - phil)
                                  		+ x3Flux[k+1][j][i-1].d*(phir - phic));

      				}
    			}
  		}
	}


/*--- Step 5d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x2-Flux and x3-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (k=kl+1; k<=ku-1; k++) {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        phic = pG->Phi[k][j][i];
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j+1][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j-1][i]);

/* correct right states; x2 and x3 gradients */
        Ur_x1Face[k][j][i].My -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;

        Ur_x1Face[k][j][i].E -= hdtodx2*(x2Flux[k][j  ][i  ].d*(phic - phil)
                                  + x2Flux[k][j+1][i  ].d*(phir - phic));


        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k+1][j][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k-1][j][i]);

        Ur_x1Face[k][j][i].Mz -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;

        Ur_x1Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j][i  ].d*(phic - phil)
                                  + x3Flux[k+1][j][i  ].d*(phir - phic));


/* correct left states; x2 and x3 gradients */
        phic = pG->Phi[k][j][i-1];
        phir = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j+1][i-1]);
        phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j-1][i-1]);

        Ul_x1Face[k][j][i].My -= hdtodx2*(phir-phil)*pG->U[k][j][i-1].d;

        Ul_x1Face[k][j][i].E -= hdtodx2*(x2Flux[k][j  ][i-1].d*(phic - phil)
                                  + x2Flux[k][j+1][i-1].d*(phir - phic));


        phir = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k+1][j][i-1]);
        phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k-1][j][i-1]);

        Ul_x1Face[k][j][i].Mz -= hdtodx3*(phir-phil)*pG->U[k][j][i-1].d;

        Ul_x1Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j][i-1].d*(phic - phil)
                                  + x3Flux[k+1][j][i-1].d*(phir - phic));

      }
    }
  }
#endif /* SELF_GRAVITY */




/*============================================================================*/
/*== Add radiation momentum source terms along the perpendicular direction */
/*
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {


				
				density = pG->U[k][j][i].d;
				velocity_x = pG->U[k][j][i].M1 / density;
				velocity_y = pG->U[k][j][i].M2 / density;
				velocity_z = pG->U[k][j][i].M3 / density;

				
				Sigma_t    = pG->U[k][j][i].Sigma_t;
				Sigma_a	   = pG->U[k][j][i].Sigma_a;
					
#ifdef FARGO

					
				cc_pos(pG,i,j,k,&x1,&x2,&x3);

				velocity_y -= qom * x1;						
					
#endif			
			
				
				Fr0y = pG->U[k][j][i].Fr2 - ((1.0 + pG->U[k][j][i].Edd_22) * velocity_y + pG->U[k][j][i].Edd_21 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_z) * pG->U[k][j][i].Er / Crat;
				Fr0z = pG->U[k][j][i].Fr3 - ((1.0 + pG->U[k][j][i].Edd_33) * velocity_z + pG->U[k][j][i].Edd_31 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_y) * pG->U[k][j][i].Er / Crat;


				
				Source[2] = -Prat * (-Sigma_t * Fr0y );
				Source[3] = -Prat * (-Sigma_t * Fr0z);


				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_y * Fr0y * velocity_z * Fr0z)/Crat);

				Ur_x1Face[k][j][i].My += hdt * Source[2];
				Ur_x1Face[k][j][i].Mz += hdt * Source[3];
				
				Ur_x1Face[k][j][i].E  += hdt * Source[4];



				
				density = pG->U[k][j][i-1].d;
				velocity_x = pG->U[k][j][i-1].M1 / density;
				velocity_y = pG->U[k][j][i-1].M2 / density;
				velocity_z = pG->U[k][j][i-1].M3 / density;

				
				Sigma_t    = pG->U[k][j][i-1].Sigma_t;
				Sigma_a	   = pG->U[k][j][i-1].Sigma_a;
					
#ifdef FARGO
			
					
				cc_pos(pG,i-1,j,k,&x1,&x2,&x3);
		

				velocity_y -= qom * x1;						
					
#endif			
				
				Fr0y = pG->U[k][j][i-1].Fr2 - ((1.0 + pG->U[k][j][i-1].Edd_22) * velocity_y + pG->U[k][j][i-1].Edd_21 * velocity_x + pG->U[k][j][i-1].Edd_32 * velocity_z) * pG->U[k][j][i-1].Er / Crat;
				Fr0z = pG->U[k][j][i-1].Fr3 - ((1.0 + pG->U[k][j][i-1].Edd_33) * velocity_z + pG->U[k][j][i-1].Edd_31 * velocity_x + pG->U[k][j][i-1].Edd_32 * velocity_y) * pG->U[k][j][i-1].Er / Crat;


				
				Source[2] = -Prat * (-Sigma_t * Fr0y);
				Source[3] = -Prat * (-Sigma_t * Fr0z);

				
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_y * Fr0y * velocity_z * Fr0z)/Crat);

				Ul_x1Face[k][j][i].My += hdt * Source[2];
				Ul_x1Face[k][j][i].Mz += hdt * Source[3];
				
				Ul_x1Face[k][j][i].E  += hdt * Source[4];

			}
		}
	}
#endif

*/
/*===========================================================================*/
/*=== STEP 7A: Correct x2-interface states with transverse flux gradients =====*/

	
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju; j++) {
      			for (i=il+1; i<=iu-1; i++) {
 
       				Ul_x2Face[k][j][i].d -=hdtodx1 * (x1Flux[k][j-1][i+1].d  - x1Flux[k][j-1][i].d );
        			Ul_x2Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k][j-1][i+1].My - x1Flux[k][j-1][i].My);
        			Ul_x2Face[k][j][i].My-=hdtodx1 * (x1Flux[k][j-1][i+1].Mz - x1Flux[k][j-1][i].Mz);
        			Ul_x2Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k][j-1][i+1].Mx - x1Flux[k][j-1][i].Mx);
       				Ul_x2Face[k][j][i].E -=hdtodx1 * (x1Flux[k][j-1][i+1].E - x1Flux[k][j-1][i].E );

#ifdef RADIATION_MHD
/* Update B3 */
				Ul_x2Face[k][j][i].By -=hdtodx1 * 0.5*
	  				((emf2[k  ][j-1][i+1] - emf2[k  ][j-1][i]) + 
	   				(emf2[k+1][j-1][i+1] - emf2[k+1][j-1][i]));
#endif

        			Ur_x2Face[k][j][i].d -=hdtodx1 * (x1Flux[k][j  ][i+1].d  - x1Flux[k][j  ][i].d );
        			Ur_x2Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k][j  ][i+1].My - x1Flux[k][j  ][i].My);
        			Ur_x2Face[k][j][i].My-=hdtodx1 * (x1Flux[k][j  ][i+1].Mz - x1Flux[k][j  ][i].Mz);
        			Ur_x2Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k][j  ][i+1].Mx - x1Flux[k][j  ][i].Mx);

        			Ur_x2Face[k][j][i].E -=hdtodx1 * (x1Flux[k][j  ][i+1].E  - x1Flux[k][j  ][i].E );
#ifdef RADIATION_MHD
/* Update B3 */
				Ur_x2Face[k][j][i].By-=hdtodx1 * 0.5*
	  				((emf2[k  ][j][i+1] - emf2[k  ][j][i]) + 
	   				(emf2[k+1][j][i+1] - emf2[k+1][j][i]));
#endif


/*--- Step 7b ------------------------------------------------------------------
 * Correct x2-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

        			Ul_x2Face[k][j][i].d -=hdtodx3 * (x3Flux[k+1][j-1][i].d -x3Flux[k][j-1][i].d );
        			Ul_x2Face[k][j][i].Mx-=hdtodx3 * (x3Flux[k+1][j-1][i].Mz-x3Flux[k][j-1][i].Mz);
        			Ul_x2Face[k][j][i].My-=hdtodx3 * (x3Flux[k+1][j-1][i].Mx-x3Flux[k][j-1][i].Mx);
        			Ul_x2Face[k][j][i].Mz-=hdtodx3 * (x3Flux[k+1][j-1][i].My-x3Flux[k][j-1][i].My);

        			Ul_x2Face[k][j][i].E -=hdtodx3 * (x3Flux[k+1][j-1][i].E -x3Flux[k][j-1][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
				Ul_x2Face[k][j][i].Bz+=hdtodx3 * 0.5*
	  				((emf2[k+1][j-1][i  ] - emf2[k][j-1][i  ]) +
	   				(emf2[k+1][j-1][i+1] - emf2[k][j-1][i+1]));
#endif

        			Ur_x2Face[k][j][i].d -=hdtodx3 * (x3Flux[k+1][j  ][i].d -x3Flux[k][j  ][i].d );
        			Ur_x2Face[k][j][i].Mx-=hdtodx3 * (x3Flux[k+1][j  ][i].Mz-x3Flux[k][j  ][i].Mz);
        			Ur_x2Face[k][j][i].My-=hdtodx3 * (x3Flux[k+1][j  ][i].Mx-x3Flux[k][j  ][i].Mx);
        			Ur_x2Face[k][j][i].Mz-=hdtodx3 * (x3Flux[k+1][j  ][i].My-x3Flux[k][j  ][i].My);

        			Ur_x2Face[k][j][i].E -=hdtodx3 * (x3Flux[k+1][j  ][i].E -x3Flux[k][j  ][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
				Ur_x2Face[k][j][i].Bz+=hdtodx3 * 0.5*
	  				((emf2[k+1][j][i  ] - emf2[k][j][i  ]) +
	   				(emf2[k+1][j][i+1] - emf2[k][j][i+1]));
#endif

      				}
    			}
  		}

		
/*--- Step 7c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x3-flux-gradients to the
 * conservative variables on the x2Face.  Limiting is used as in GS (2007)
 */

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju; j++) {
      			for (i=il+1; i<=iu-1; i++) {

        			db1 = (	   pG->B1i[k  ][j-1][i+1] -     pG->B1i[k][j-1][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j  ][i  ] -     pG->B2i[k][j-1][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j-1][i  ] -     pG->B3i[k][j-1][i])*dx3i;
        			B1 = pG->U[k][j-1][i].B1c;
        			B2 = pG->U[k][j-1][i].B2c;
        			B3 = pG->U[k][j-1][i].B3c;
        			V1 = pG->U[k][j-1][i].M1/pG->U[k][j-1][i].d;
        			V3 = pG->U[k][j-1][i].M3/pG->U[k][j-1][i].d;

/* Calculate mdb1 = min_mod(-db2,db1) */
        			if(db2 > 0.0 && db1 < 0.0){
          				mdb1 = db1 > -db2 ? db1 : -db2;
        			}
        			else if(db2 < 0.0 && db1 > 0.0){
          				mdb1 = db1 < -db2 ? db1 : -db2;
        			}
        			else mdb1 = 0.0;

/* Calculate mdb3 = min_mod(-db2,db3) */
        			if(db2 > 0.0 && db3 < 0.0){
          				mdb3 = db3 > -db2 ? db3 : -db2;
        			}
        			else if(db2 < 0.0 && db3 > 0.0){
          				mdb3 = db3 < -db2 ? db3 : -db2;
        			}
        			else mdb3 = 0.0;

        			Ul_x2Face[k][j][i].Mz += hdt*B1*db2;
        			Ul_x2Face[k][j][i].Mx += hdt*B2*db2;
        			Ul_x2Face[k][j][i].My += hdt*B3*db2;
        			Ul_x2Face[k][j][i].By += hdt*V3*(-mdb1);
        			Ul_x2Face[k][j][i].Bz += hdt*V1*(-mdb3);

        			Ul_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
        			V1 = pG->U[k][j][i].M1/pG->U[k][j][i].d;
        			V3 = pG->U[k][j][i].M3/pG->U[k][j][i].d;

/* Calculate mdb1 = min_mod(-db2,db1) */
        			if(db2 > 0.0 && db1 < 0.0){
          				mdb1 = db1 > -db2 ? db1 : -db2;
        			}
        			else if(db2 < 0.0 && db1 > 0.0){
          				mdb1 = db1 < -db2 ? db1 : -db2;
        			}
       		 		else mdb1 = 0.0;

/* Calculate mdb3 = min_mod(-db2,db3) */
        			if(db2 > 0.0 && db3 < 0.0){
          				mdb3 = db3 > -db2 ? db3 : -db2;
        			}
        			else if(db2 < 0.0 && db3 > 0.0){
          				mdb3 = db3 < -db2 ? db3 : -db2;
        			}
        			else mdb3 = 0.0;

        			Ur_x2Face[k][j][i].Mz += hdt*B1*db2;
        			Ur_x2Face[k][j][i].Mx += hdt*B2*db2;
        			Ur_x2Face[k][j][i].My += hdt*B3*db2;
        			Ur_x2Face[k][j][i].By += hdt*V3*(-mdb1);
        			Ur_x2Face[k][j][i].Bz += hdt*V1*(-mdb3);
        			Ur_x2Face[k][j][i].E  += hdt*(B3*V3*(-mdb1) + B1*V1*(-mdb3) );

      				}
    			}
  		}
#endif /* RADIATION MHD */

				
/*--- Step 7d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x3-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  	if (StaticGravPot != NULL){
  		for (k=kl+1; k<=ku-1; k++) {
    			for (j=jl+1; j<=ju; j++) {
      				for (i=il+1; i<=iu-1; i++) {
        				cc_pos(pG,i,j,k,&x1,&x2,&x3);
        
			phic = (*StaticGravPot)((x1            ),x2,x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* correct right states; x1 and x3 gradients */

        		Ur_x2Face[k][j][i].Mz -= hdtodx1 * (phir-phil)*pG->U[k][j][i].d;

        		Ur_x2Face[k][j][i].E -= hdtodx1 * (x1Flux[k][j  ][i  ].d*(phic - phil)
                                  + x1Flux[k][j  ][i+1].d*(phir - phic));

        		phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        		phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

        		Ur_x2Face[k][j][i].My -= hdtodx3 * (phir-phil)*pG->U[k][j][i].d;

        		Ur_x2Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j  ][i].d*(phic - phil)
                                  	+ x3Flux[k+1][j  ][i].d*(phir - phic));


/* correct left states; x1 and x3 gradients */
        		phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

        		Ul_x2Face[k][j][i].Mz -= hdtodx1 * (phir-phil)*pG->U[k][j-1][i].d;
        		Ul_x2Face[k][j][i].E -= hdtodx1 * (x1Flux[k][j-1][i  ].d*(phic - phil)
                                  + x1Flux[k][j-1][i+1].d*(phir - phic));

        		phir = (*StaticGravPot)(x1,(x2-pG->dx2),(x3+0.5*pG->dx3));
        		phil = (*StaticGravPot)(x1,(x2-pG->dx2),(x3-0.5*pG->dx3));

        		Ul_x2Face[k][j][i].My -= hdtodx3 * (phir-phil)*pG->U[k][j-1][i].d;

         		Ul_x2Face[k][j][i].E -= hdtodx3 * (x3Flux[k  ][j-1][i].d*(phic - phil)
                                  + x3Flux[k+1][j-1][i].d*(phir - phic));

      				}
    			}
  		}
	}




/*--- Step 6d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x1-Flux and x3-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (k=kl+1; k<=ku-1; k++) {
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        phic = pG->Phi[k][j][i];
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i+1]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i-1]);

/* correct right states; x1 and x3 gradients */
        Ur_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;

        Ur_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j][i  ].d*(phic - phil)
                                  + x1Flux[k][j][i+1].d*(phir - phic));


        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k+1][j][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k-1][j][i]);

        Ur_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;

        Ur_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j][i].d*(phic - phil)
                                  + x3Flux[k+1][j][i].d*(phir - phic));

/* correct left states; x1 and x3 gradients */
        phic = pG->Phi[k][j-1][i];
        phir = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k][j-1][i+1]);
        phil = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k][j-1][i-1]);

        Ul_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j-1][i].d;

        Ul_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j-1][i  ].d*(phic - phil)
                                  + x1Flux[k][j-1][i+1].d*(phir - phic));
        phir = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k+1][j-1][i]);
        phil = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k-1][j-1][i]);

        Ul_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j-1][i].d;

        Ul_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j-1][i].d*(phic - phil)
                                  + x3Flux[k+1][j-1][i].d*(phir - phic));

      }
    }
  }
#endif /* SELF_GRAVITY */



/*========================Shearing box source term gradient *===============*/	
	
/* Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
* x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
*    Vx source term is (dt/2)( 2 Omega_0 Vy)
*    Vy source term is (dt/2)(-2 Omega_0 Vx)
*    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
*/
		
#ifdef SHEARING_BOX
		if (ShearingBoxPot != NULL){
			for (k=kl+1; k<=ku-1; k++) {
				for (j=jl+1; j<=ju; j++) {
					for (i=il+1; i<=iu-1; i++) {
						cc_pos(pG,i,j,k,&x1,&x2,&x3);
						phic = (*ShearingBoxPot)((x1            ),x2,x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
						
						/* correct right states; x1 and x3 gradients */
						Ur_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
						Ur_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j  ][i  ].d*(phic - phil)
													+ x1Flux[k][j  ][i+1].d*(phir - phic));
#endif
						
						phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
						
						Ur_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
						Ur_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j  ][i].d*(phic - phil)
													+ x3Flux[k+1][j  ][i].d*(phir - phic));
#endif
						
						/* correct left states; x1 and x3 gradients */
						phic = (*ShearingBoxPot)((x1            ),(x2-pG->dx2),x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);
						
						Ul_x2Face[k][j][i].Mz -= hdtodx1*(phir-phil)*pG->U[k][j-1][i].d;
#ifndef BAROTROPIC
						Ul_x2Face[k][j][i].E -= hdtodx1*(x1Flux[k][j-1][i  ].d*(phic - phil)
													+ x1Flux[k][j-1][i+1].d*(phir - phic));
#endif
						phir = (*ShearingBoxPot)(x1,(x2-pG->dx2),(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,(x2-pG->dx2),(x3-0.5*pG->dx3));
						
						Ul_x2Face[k][j][i].My -= hdtodx3*(phir-phil)*pG->U[k][j-1][i].d;
#ifndef BAROTROPIC
						Ul_x2Face[k][j][i].E -= hdtodx3*(x3Flux[k  ][j-1][i].d*(phic - phil)
													+ x3Flux[k+1][j-1][i].d*(phir - phic));
#endif
					}
				}
			}
		}
	
	for (k=kl+1; k<=ku-1; k++) {
		for (j=jl+1; j<=ju; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x2Face[k][j][i].Mz += pG->dt*Omega_0*pG->U[k][j][i].M2;
				Ul_x2Face[k][j][i].Mz += pG->dt*Omega_0*pG->U[k][j-1][i].M2;
#ifdef FARGO
				Ur_x2Face[k][j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
				Ul_x2Face[k][j][i].Mx += hdt*(qshear-2.)*Omega_0*pG->U[k][j-1][i].M1;
#else
				Ur_x2Face[k][j][i].Mx -= pG->dt*Omega_0*pG->U[k][j][i].M1;
				Ul_x2Face[k][j][i].Mx -= pG->dt*Omega_0*pG->U[k][j-1][i].M1;
#endif
			}
		}
	}

#endif /* SHEARING_BOX */
	
/*========================Shearing box source term gradient *===============*/	


/*============================================================================*/
/*== Add radiation momentum source terms along the perpendicular direction */
/*
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

				density = pG->U[k][j][i].d;
				velocity_x = pG->U[k][j][i].M1 / density;
				velocity_y = pG->U[k][j][i].M2 / density;
				velocity_z = pG->U[k][j][i].M3 / density;

				
				Sigma_t    = pG->U[k][j][i].Sigma_t;
				Sigma_a	   = pG->U[k][j][i].Sigma_a;
					
#ifdef FARGO
				
				cc_pos(pG,i,j,k,&x1,&x2,&x3);
	

				velocity_y -= qom * x1;						
					
#endif			
				
				
				Fr0x = pG->U[k][j][i].Fr1 - ((1.0 + pG->U[k][j][i].Edd_11) * velocity_x + pG->U[k][j][i].Edd_21 * velocity_y + pG->U[k][j][i].Edd_31 * velocity_z) * pG->U[k][j][i].Er / Crat;
				Fr0z = pG->U[k][j][i].Fr3 - ((1.0 + pG->U[k][j][i].Edd_33) * velocity_z + pG->U[k][j][i].Edd_31 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_y) * pG->U[k][j][i].Er / Crat;


				
				Source[1] = -Prat * (-Sigma_t * Fr0x);
				Source[3] = -Prat * (-Sigma_t * Fr0z);

				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x * velocity_z * Fr0z)/Crat);

				Ur_x1Face[k][j][i].Mz += hdt * Source[1];
				Ur_x1Face[k][j][i].My += hdt * Source[3];
				
				Ur_x1Face[k][j][i].E  += hdt * Source[4];


		
				
				density = pG->U[k][j-1][i].d;
				velocity_x = pG->U[k][j-1][i].M1 / density;
				velocity_y = pG->U[k][j-1][i].M2 / density;
				velocity_z = pG->U[k][j-1][i].M3 / density;

				
				Sigma_t    = pG->U[k][j-1][i].Sigma_t;
				Sigma_a	   = pG->U[k][j-1][i].Sigma_a;
					
#ifdef FARGO
			
					
				cc_pos(pG,i,j-1,k,&x1,&x2,&x3);
		

				velocity_y -= qom * x1;						
					
#endif			
			
				Fr0x = pG->U[k][j-1][i].Fr1 - ((1.0 + pG->U[k][j-1][i].Edd_11) * velocity_x + pG->U[k][j-1][i].Edd_21 * velocity_y + pG->U[k][j-1][i].Edd_31 * velocity_z) * pG->U[k][j-1][i].Er / Crat;
				Fr0z = pG->U[k][j-1][i].Fr3 - ((1.0 + pG->U[k][j-1][i].Edd_33) * velocity_z + pG->U[k][j-1][i].Edd_31 * velocity_x + pG->U[k][j-1][i].Edd_32 * velocity_y) * pG->U[k][j-1][i].Er / Crat;


				
				Source[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
				Source[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);

	
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x * velocity_z * Fr0z)/Crat);

				Ur_x1Face[k][j-1][i].Mz += hdt * Source[1];
				Ur_x1Face[k][j-1][i].My += hdt * Source[3];
				
				Ur_x1Face[k][j-1][i].E  += hdt * Source[4];


			}
		}
	}
#endif

*/
	
/*=== STEP 8: Correct x3-interface states with transverse flux gradients =====*/

/*--- Step 8a ------------------------------------------------------------------
 * Correct x3-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */

  	for (k=kl+1; k<=ku; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {

        		Ul_x3Face[k][j][i].d -=hdtodx1 * (x1Flux[k-1][j][i+1].d -x1Flux[k-1][j][i].d );
        		Ul_x3Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k-1][j][i+1].Mz-x1Flux[k-1][j][i].Mz);
        		Ul_x3Face[k][j][i].My-=hdtodx1 * (x1Flux[k-1][j][i+1].Mx-x1Flux[k-1][j][i].Mx);
        		Ul_x3Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k-1][j][i+1].My-x1Flux[k-1][j][i].My);

        		Ul_x3Face[k][j][i].E -=hdtodx1 * (x1Flux[k-1][j][i+1].E -x1Flux[k-1][j][i].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ul_x3Face[k][j][i].Bz+=hdtodx1 * 0.5*
	  			((emf3[k-1][j  ][i+1] - emf3[k-1][j  ][i]) +
	   			(emf3[k-1][j+1][i+1] - emf3[k-1][j+1][i]));
#endif

        		Ur_x3Face[k][j][i].d -=hdtodx1 * (x1Flux[k  ][j][i+1].d  - x1Flux[k  ][j][i].d );
        		Ur_x3Face[k][j][i].Mx-=hdtodx1 * (x1Flux[k  ][j][i+1].Mz - x1Flux[k  ][j][i].Mz);
        		Ur_x3Face[k][j][i].My-=hdtodx1 * (x1Flux[k  ][j][i+1].Mx - x1Flux[k  ][j][i].Mx);
        		Ur_x3Face[k][j][i].Mz-=hdtodx1 * (x1Flux[k  ][j][i+1].My - x1Flux[k  ][j][i].My);

        		Ur_x3Face[k][j][i].E -=hdtodx1 * (x1Flux[k  ][j][i+1].E - x1Flux[k  ][j][i].E );

#ifdef RADIATION_MHD
/* Update B2 */
			Ur_x3Face[k][j][i].Bz+=hdtodx1 * 0.5*
	  			((emf3[k][j  ][i+1] - emf3[k][j  ][i]) +
	   			(emf3[k][j+1][i+1] - emf3[k][j+1][i]));
#endif


/*--- Step 8b ------------------------------------------------------------------
 * Correct x3-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (y,z,x) on LHS 
 */

        		Ul_x3Face[k][j][i].d -=hdtodx2 * (x2Flux[k-1][j+1][i].d -x2Flux[k-1][j][i].d );
        		Ul_x3Face[k][j][i].Mx-=hdtodx2 * (x2Flux[k-1][j+1][i].My-x2Flux[k-1][j][i].My);
        		Ul_x3Face[k][j][i].My-=hdtodx2 * (x2Flux[k-1][j+1][i].Mz-x2Flux[k-1][j][i].Mz);
        		Ul_x3Face[k][j][i].Mz-=hdtodx2 * (x2Flux[k-1][j+1][i].Mx-x2Flux[k-1][j][i].Mx);

        		Ul_x3Face[k][j][i].E -=hdtodx2 * (x2Flux[k-1][j+1][i].E -x2Flux[k-1][j][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
			Ul_x3Face[k][j][i].By-=hdtodx2 * 0.5*
	  			((emf3[k-1][j+1][i  ] - emf3[k-1][j][i  ]) +
	   			(emf3[k-1][j+1][i+1] - emf3[k-1][j][i+1]));
#endif

        		Ur_x3Face[k][j][i].d -=hdtodx2 * (x2Flux[k  ][j+1][i].d -x2Flux[k  ][j][i].d );
        		Ur_x3Face[k][j][i].Mx-=hdtodx2 * (x2Flux[k  ][j+1][i].My-x2Flux[k  ][j][i].My);
        		Ur_x3Face[k][j][i].My-=hdtodx2 * (x2Flux[k  ][j+1][i].Mz-x2Flux[k  ][j][i].Mz);
        		Ur_x3Face[k][j][i].Mz-=hdtodx2 * (x2Flux[k  ][j+1][i].Mx-x2Flux[k  ][j][i].Mx);

        		Ur_x3Face[k][j][i].E -=hdtodx2 * (x2Flux[k  ][j+1][i].E -x2Flux[k  ][j][i].E );

#ifdef RADIATION_MHD
/* Update B1 */
			Ur_x3Face[k][j][i].By-=hdtodx2 * 0.5*
	  			((emf3[k][j+1][i  ] - emf3[k][j][i  ]) +
	   			(emf3[k][j+1][i+1] - emf3[k][j][i+1]));
#endif

      			}
    		}
  	}

		
/*--- Step 8c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x2-flux-gradients to the
 * conservative variables on the x3Face.  Limiting is used as in GS07.
 */

#ifdef RADIATION_MHD
  	for (k=kl+1; k<=ku; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {


        			db1 = (    pG->B1i[k-1][j  ][i+1] -     pG->B1i[k-1][j][i])*dx1i;
        			db2 = (    pG->B2i[k-1][j+1][i  ] -     pG->B2i[k-1][j][i])*dx2i;
        			db3 = (    pG->B3i[k  ][j  ][i  ] -     pG->B3i[k-1][j][i])*dx3i;
        			B1 = pG->U[k-1][j][i].B1c;
        			B2 = pG->U[k-1][j][i].B2c;
        			B3 = pG->U[k-1][j][i].B3c;
				V1 = pG->U[k-1][j][i].M1/pG->U[k-1][j][i].d;
				V2 = pG->U[k-1][j][i].M2/pG->U[k-1][j][i].d;

/* Calculate mdb1 = min_mod(-db3,db1) */
				if(db3 > 0.0 && db1 < 0.0){
	  				mdb1 = db1 > -db3 ? db1 : -db3;
				}
				else if(db3 < 0.0 && db1 > 0.0){
	  				mdb1 = db1 < -db3 ? db1 : -db3;
				}
				else mdb1 = 0.0;

/* Calculate mdb2 = min_mod(-db3,db2) */
				if(db3 > 0.0 && db2 < 0.0){
	  				mdb2 = db2 > -db3 ? db2 : -db3;
				}
				else if(db3 < 0.0 && db2 > 0.0){
	  				mdb2 = db2 < -db3 ? db2 : -db3;
				}
				else mdb2 = 0.0;

        			Ul_x3Face[k][j][i].My += hdt*B1*db3;
        			Ul_x3Face[k][j][i].Mz += hdt*B2*db3;
        			Ul_x3Face[k][j][i].Mx += hdt*B3*db3;
				Ul_x3Face[k][j][i].By += hdt*V1*(-mdb2);
				Ul_x3Face[k][j][i].Bz += hdt*V2*(-mdb1);

				Ul_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );


        			db1 = (    pG->B1i[k  ][j  ][i+1] -     pG->B1i[k][j][i])*dx1i;
        			db2 = (    pG->B2i[k  ][j+1][i  ] -     pG->B2i[k][j][i])*dx2i;
        			db3 = (    pG->B3i[k+1][j  ][i  ] -     pG->B3i[k][j][i])*dx3i;
        			B1 = pG->U[k][j][i].B1c;
        			B2 = pG->U[k][j][i].B2c;
        			B3 = pG->U[k][j][i].B3c;
				V1 = pG->U[k][j][i].M1/pG->U[k][j][i].d;
				V2 = pG->U[k][j][i].M2/pG->U[k][j][i].d;

/* Calculate mdb1 = min_mod(-db3,db1) */
				if(db3 > 0.0 && db1 < 0.0){
	  				mdb1 = db1 > -db3 ? db1 : -db3;
				}
				else if(db3 < 0.0 && db1 > 0.0){
	  				mdb1 = db1 < -db3 ? db1 : -db3;
				}
				else mdb1 = 0.0;	

/* Calculate mdb2 = min_mod(-db3,db2) */
				if(db3 > 0.0 && db2 < 0.0){
	  				mdb2 = db2 > -db3 ? db2 : -db3;
				}
				else if(db3 < 0.0 && db2 > 0.0){
	  				mdb2 = db2 < -db3 ? db2 : -db3;
				}
				else mdb2 = 0.0;

        			Ur_x3Face[k][j][i].My += hdt*B1*db3;
        			Ur_x3Face[k][j][i].Mz += hdt*B2*db3;
        			Ur_x3Face[k][j][i].Mx += hdt*B3*db3;
				Ur_x3Face[k][j][i].By += hdt*V1*(-mdb2);
				Ur_x3Face[k][j][i].Bz += hdt*V2*(-mdb1);

				Ur_x3Face[k][j][i].E  += hdt*(B1*V1*(-mdb2) + B2*V2*(-mdb1) );

      				}
    			}
  		}
#endif /* MHD */


/*--- Step 8d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x2-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

  	if (StaticGravPot != NULL){
  		for (k=kl+1; k<=ku; k++) {
    			for (j=jl+1; j<=ju-1; j++) {
      				for (i=il+1; i<=iu-1; i++) {
        				cc_pos(pG,i,j,k,&x1,&x2,&x3);
        				phic = (*StaticGravPot)((x1            ),x2,x3);
        				phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        				phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

/* correct right states; x1 and x2 gradients */

        				Ur_x3Face[k][j][i].My -= hdtodx1 * (phir-phil)*pG->U[k][j][i].d;

        				Ur_x3Face[k][j][i].E -= hdtodx1 * (x1Flux[k  ][j][i  ].d*(phic - phil)
                                  			+ x1Flux[k  ][j][i+1].d*(phir - phic));


        				phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        				phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        				Ur_x3Face[k][j][i].Mz -= hdtodx2 * (phir-phil)*pG->U[k][j][i].d;

        				Ur_x3Face[k][j][i].E -= hdtodx2 * (x2Flux[k  ][j  ][i].d*(phic - phil)
                                  		+ x2Flux[k  ][j+1][i].d*(phir - phic));


/* correct left states; x1 and x2 gradients */
        				phic = (*StaticGravPot)((x1            ),x2,(x3-pG->dx3));
        				phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,(x3-pG->dx3));
        				phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,(x3-pG->dx3));


        				Ul_x3Face[k][j][i].My -= hdtodx1 * (phir-phil)*pG->U[k-1][j][i].d;

        				Ul_x3Face[k][j][i].E -= hdtodx1 * (x1Flux[k-1][j][i  ].d*(phic - phil)
                                  		+ x1Flux[k-1][j][i+1].d*(phir - phic));


        				phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),(x3-pG->dx3));
        				phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),(x3-pG->dx3));

        				Ul_x3Face[k][j][i].Mz -= hdtodx2 * (phir-phil)*pG->U[k-1][j][i].d;

        				Ul_x3Face[k][j][i].E -= hdtodx2 * (x2Flux[k-1][j  ][i].d*(phic - phil)
                                  		+ x2Flux[k-1][j+1][i].d*(phir - phic));

      						}
    					}
  				}
			}




/*--- Step 7d (cont) -----------------------------------------------------------
 * Add source terms for self gravity arising from x1-Flux and x2-Flux gradients
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

#ifdef SELF_GRAVITY
  for (k=kl+1; k<=ku; k++) {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {
        phic = pG->Phi[k][j][i];
        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i+1]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i-1]);

/* correct right states; x1 and x2 gradients */
        Ur_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;

        Ur_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k][j][i  ].d*(phic - phil)
                                  + x1Flux[k][j][i+1].d*(phir - phic));

        phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j+1][i]);
        phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j-1][i]);

        Ur_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;

        Ur_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k][j  ][i].d*(phic - phil)
                                  + x2Flux[k][j+1][i].d*(phir - phic));


/* correct left states; x1 and x2 gradients */
        phic = pG->Phi[k-1][j][i];
        phir = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k-1][j][i+1]);
        phil = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k-1][j][i-1]);

        Ul_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k-1][j][i].d;

        Ul_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k-1][j][i  ].d*(phic - phil)
                                  + x1Flux[k-1][j][i+1].d*(phir - phic));


        phir = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k-1][j+1][i]);
        phil = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k-1][j-1][i]);

        Ul_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k-1][j][i].d;

        Ul_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k-1][j  ][i].d*(phic - phil)
                                  + x2Flux[k-1][j+1][i].d*(phir - phic));

      }
    }
  }
#endif /* SELF_GRAVITY */


	
/*=============== Shearing box source terms ================*/
	
	/*--- Step 7d (cont) -----------------------------------------------------------
	 * Add source terms for shearing box (Coriolis forces) for 0.5*dt arising from
	 * x1-Flux gradient.  The tidal gravity terms are added via ShearingBoxPot
	 *    Vx source term is (dt/2)( 2 Omega_0 Vy)
	 *    Vy source term is (dt/2)(-2 Omega_0 Vx)
	 *    Vy source term is (dt/2)((q-2) Omega_0 Vx) (with FARGO)
	 */
	
#ifdef SHEARING_BOX
	if (ShearingBoxPot != NULL){
		for (k=kl+1; k<=ku; k++) {
			for (j=jl+1; j<=ju-1; j++) {
				for (i=il+1; i<=iu-1; i++) {
					cc_pos(pG,i,j,k,&x1,&x2,&x3);
					phic = (*ShearingBoxPot)((x1            ),x2,x3);
					phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
					phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
					
					/* correct right states; x1 and x2 gradients */
					Ur_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
					Ur_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k  ][j][i  ].d*(phic - phil)
												+ x1Flux[k  ][j][i+1].d*(phir - phic));
#endif
					
					phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
					phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
					
					Ur_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
					Ur_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k  ][j  ][i].d*(phic - phil)
												+ x2Flux[k  ][j+1][i].d*(phir - phic));
#endif
					
					/* correct left states; x1 and x2 gradients */
					phic = (*ShearingBoxPot)((x1            ),x2,(x3-pG->dx3));
					phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,(x3-pG->dx3));
					phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,(x3-pG->dx3));
					
					Ul_x3Face[k][j][i].My -= hdtodx1*(phir-phil)*pG->U[k-1][j][i].d;
#ifndef BAROTROPIC
					Ul_x3Face[k][j][i].E -= hdtodx1*(x1Flux[k-1][j][i  ].d*(phic - phil)
												+ x1Flux[k-1][j][i+1].d*(phir - phic));
#endif
					
					phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),(x3-pG->dx3));
					phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),(x3-pG->dx3));
					
					Ul_x3Face[k][j][i].Mz -= hdtodx2*(phir-phil)*pG->U[k-1][j][i].d;
#ifndef BAROTROPIC
					Ul_x3Face[k][j][i].E -= hdtodx2*(x2Flux[k-1][j  ][i].d*(phic - phil)
												+ x2Flux[k-1][j+1][i].d*(phir - phic));
#endif
				}
			}
		}
	}
	
	for (k=kl+1; k<=ku; k++) {
		for (j=jl+1; j<=ju-1; j++) {
			for (i=il+1; i<=iu-1; i++) {
				Ur_x3Face[k][j][i].My += pG->dt*Omega_0*pG->U[k][j][i].M2;
				Ul_x3Face[k][j][i].My += pG->dt*Omega_0*pG->U[k-1][j][i].M2;
#ifdef FARGO
				Ur_x3Face[k][j][i].Mz += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
				Ul_x3Face[k][j][i].Mz += hdt*(qshear-2.)*Omega_0*pG->U[k-1][j][i].M1;
#else
				Ur_x3Face[k][j][i].Mz -= pG->dt*Omega_0*pG->U[k][j][i].M1;
				Ul_x3Face[k][j][i].Mz -= pG->dt*Omega_0*pG->U[k-1][j][i].M1;
#endif
			}
		}
	}
#endif /* SHEARING_BOX */	
	


/*===========================================================*/



/*============================================================================*/
/*== Add radiation momentum source terms along the perpendicular direction */

/*
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu; i++) {

				
				
				density = pG->U[k][j][i].d;
				velocity_x = pG->U[k][j][i].M1 / density;
				velocity_y = pG->U[k][j][i].M2 / density;
				velocity_z = pG->U[k][j][i].M3 / density;

				
				Sigma_t    = pG->U[k][j][i].Sigma_t;
				Sigma_a	   = pG->U[k][j][i].Sigma_a;
					
#ifdef FARGO
				
					
				cc_pos(pG,i,j,k,&x1,&x2,&x3);
				

				velocity_y -= qom * x1;						
					
#endif			
				
				
				Fr0x = pG->U[k][j][i].Fr1 - ((1.0 + pG->U[k][j][i].Edd_11) * velocity_x + pG->U[k][j][i].Edd_21 * velocity_y + pG->U[k][j][i].Edd_31 * velocity_z) * pG->U[k][j][i].Er / Crat;
				Fr0y = pG->U[k][j][i].Fr2 - ((1.0 + pG->U[k][j][i].Edd_22) * velocity_y + pG->U[k][j][i].Edd_21 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_z) * pG->U[k][j][i].Er / Crat;


				
				Source[1] = -Prat * (-Sigma_t * Fr0x);
				Source[2] = -Prat * (-Sigma_t * Fr0y);

				
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x * velocity_y * Fr0y)/Crat);

				Ur_x1Face[k][j][i].My += hdt * Source[1];
				Ur_x1Face[k][j][i].Mz += hdt * Source[2];
				
				Ur_x1Face[k][j][i].E  += hdt * Source[4];


				
				density = pG->U[k-1][j][i].d;
				velocity_x = pG->U[k-1][j][i].M1 / density;
				velocity_y = pG->U[k-1][j][i].M2 / density;
				velocity_z = pG->U[k-1][j][i].M3 / density;

				
				Sigma_t    = pG->U[k-1][j][i].Sigma_t;
				Sigma_a	   = pG->U[k-1][j][i].Sigma_a;
					
#ifdef FARGO
			
					
				cc_pos(pG,i,j,k-1,&x1,&x2,&x3);
			

				velocity_y -= qom * x1;						
					
#endif			
			
				
				Fr0x = pG->U[k-1][j][i].Fr1 - ((1.0 + pG->U[k-1][j][i].Edd_11) * velocity_x + pG->U[k-1][j][i].Edd_21 * velocity_y + pG->U[k-1][j][i].Edd_31 * velocity_z) * pG->U[k-1][j][i].Er / Crat;
				Fr0y = pG->U[k-1][j][i].Fr2 - ((1.0 + pG->U[k-1][j][i].Edd_22) * velocity_y + pG->U[k-1][j][i].Edd_21 * velocity_x + pG->U[k-1][j][i].Edd_32 * velocity_z) * pG->U[k-1][j][i].Er / Crat;


				
				Source[1] = -Prat * (-Sigma_t * Fr0x);
				Source[2] = -Prat * (-Sigma_t * Fr0y);

		
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x * velocity_y * Fr0y)/Crat);

				Ur_x1Face[k-1][j][i].My += hdt * Source[1];
				Ur_x1Face[k-1][j][i].Mz += hdt * Source[2];
				
				Ur_x1Face[k-1][j][i].E  += hdt * Source[4];

			}
		}
	}
#endif
*/
/* End adding radiation source terms */	



/*=== STEP 9: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 9a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, cooling, or RADIATION_MHD)
 */

    	for (k=kl+1; k<=ku-1; k++) {
      		for (j=jl+1; j<=ju-1; j++) {
			for (i=il+1; i<=iu-1; i++) {

          			dhalf[k][j][i] = pG->U[k][j][i].d 
            				- hdtodx1 * (    x1Flux[k  ][j  ][i+1].d -     x1Flux[k][j][i].d)
            				- hdtodx2 * (    x2Flux[k  ][j+1][i  ].d -     x2Flux[k][j][i].d)
            				- hdtodx3 * (    x3Flux[k+1][j  ][i  ].d -     x3Flux[k][j][i].d);

			}
      		}
    	}
  

/*--- Step 9b ------------------------------------------------------------------
 * Calculate cell centered emf_cc^{n+1/2}
 */

#ifdef RADIATION_MHD
/* Update momentum for half time step. We also need to add radiation source term * 
 * using modified Godunov method */

  	for (k=kl+1; k<=ku-1; k++) {
    		for (j=jl+1; j<=ju-1; j++) {
      			for (i=il+1; i<=iu-1; i++) {

		/* load 1D vector */
			Usource.d  = pG->U[k][j][i].d;
      			Usource.Mx = pG->U[k][j][i].M1;
      			Usource.My = pG->U[k][j][i].M2;
      			Usource.Mz = pG->U[k][j][i].M3;
      			Usource.E  = pG->U[k][j][i].E;
			Usource.Er  = pG->U[k][j][i].Er;
    			Usource.Fr1  = pG->U[k][j][i].Fr1;
    			Usource.Fr2  = pG->U[k][j][i].Fr2;
    			Usource.Fr3  = pG->U[k][j][i].Fr3;
			Usource.Edd_11  = pG->U[k][j][i].Edd_11;
			Usource.Edd_21  = pG->U[k][j][i].Edd_21;
			Usource.Edd_22  = pG->U[k][j][i].Edd_22;
			Usource.Edd_31	= pG->U[k][j][i].Edd_31;
			Usource.Edd_32	= pG->U[k][j][i].Edd_32;
			Usource.Edd_33	= pG->U[k][j][i].Edd_33;
			for(m=0; m<NOPACITY; m++){
				Usource.Sigma[m] =  pG->U[k][j][i].Sigma[m];
			}
      			Usource.By = pG->U[k][j][i].B2c;
      			Usource.Bz = pG->U[k][j][i].B3c;
			Bx = pG->U[k][j][i].B1c;
			
			density = Usource.d;
			velocity_x = Usource.Mx / density;
			velocity_y = Usource.My / density;
			velocity_z = Usource.Mz / density;

			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
				
			pressure = (pG->U[k][j][i].E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */

			pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
	
			temperature = pressure / (density * R_ideal);


			Sigma_sF = Usource.Sigma[0];
			Sigma_aF = Usource.Sigma[1];
			Sigma_aP = Usource.Sigma[2];
			Sigma_aE = Usource.Sigma[3];

			diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;
		

					
#ifdef FARGO
					/* With FARGO, we should add background shearing to the source terms */
					
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			/* Include background shearing in Usource, which is only used in dSource */
					
			velocity_y -= qom * x1;		
					
#endif	

			/* The Source term */
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	

			SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);

			SFmz = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/
		

			/* Now calculate flux gradient */

			divFlux1[1] = (x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx) / dx1;
			divFlux1[2] = (x1Flux[k][j][i+1].My - x1Flux[k][j][i].My) / dx1;
			divFlux1[3] = (x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz) / dx1;
			
			divFlux2[1] = (x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz) / dx2;
			divFlux2[2] = (x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx) / dx2;
			divFlux2[3] = (x2Flux[k][j+1][i].My - x2Flux[k][j][i].My) / dx2;
			
			divFlux3[1] = (x3Flux[k+1][j][i].My - x3Flux[k][j][i].My) / dx3;
			divFlux3[2] = (x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz) / dx3;
			divFlux3[3] = (x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx) / dx3;
		
			
/* Add source terms for fixed gravitational potential */
        		if (StaticGravPot != NULL){
          			cc_pos(pG,i,j,k,&x1,&x2,&x3);


          			phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
          			phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
          			divFlux1[1] += (phir-phil)*pG->U[k][j][i].d/dx1;


          			phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
          			phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);
          			divFlux2[2] += (phir-phil)*pG->U[k][j][i].d/dx2;

          			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          			phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));
          			divFlux3[3] += (phir-phil)*pG->U[k][j][i].d/dx3;
        		}

/* Add source terms due to self_gravity */
#ifdef SELF_GRAVITY
        		phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i+1]);
        		phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j][i-1]);
			divFlux1[1] += (phir-phil)*pG->U[k][j][i].d/dx1;


        		phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j+1][i]);
        		phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k][j-1][i]);
        		divFlux2[2] += (phir-phil)*pG->U[k][j][i].d/dx2;

        		phir = 0.5*(pG->Phi[k][j][i] + pG->Phi[k+1][j][i]);
        		phil = 0.5*(pG->Phi[k][j][i] + pG->Phi[k-1][j][i]);
        		divFlux3[3] += (phir-phil)*pG->U[k][j][i].d/dx3;
#endif /* SELF_GRAVITY */


			M1h = pG->U[k][j][i].M1 + hdt * Source_Inv[1][1] * (Source[k][j][i][1] - divFlux1[1] - divFlux2[1] - divFlux3[1]);
			M2h = pG->U[k][j][i].M2 + hdt * Source_Inv[2][2] * (Source[k][j][i][2] - divFlux1[2] - divFlux2[2] - divFlux3[2]);
			M3h = pG->U[k][j][i].M3 + hdt * Source_Inv[3][3] * (Source[k][j][i][3] - divFlux1[3] - divFlux2[3] - divFlux3[3]);
					
			/* Guess solution */
					
					
/*			Uguess[1] = pG->U[k][j][i].M1 + hdt * Source_Inv[1][1] * (Source[1] - divFlux1[1] - divFlux2[1] - divFlux3[1]);
			Uguess[2] = pG->U[k][j][i].M2 + hdt * Source_Inv[2][2] * (Source[2] - divFlux1[2] - divFlux2[2] - divFlux3[2]);
			Uguess[3] = pG->U[k][j][i].M3 + hdt * Source_Inv[3][3] * (Source[3] - divFlux1[3] - divFlux2[3] - divFlux3[3]);

			 Guess the source temrm /
			velocity_x = Uguess[1] / density;
			velocity_y = Uguess[2] / density;
			velocity_z = Uguess[3] / density;
					
#ifdef FARGO
					 With FARGO, we should add background shearing to the source terms 
					cc_pos(pG,i,j,k,&x1,&x2,&x3);
					velocity_y -= qom * x1;		
#endif
					
			
	
			Source_guess[1] = -Prat * (-Sigma_t * Fr0x + Sigma_a * velocity_x * diffTEr / Crat);
			Source_guess[2] = -Prat * (-Sigma_t * Fr0y + Sigma_a * velocity_y * diffTEr / Crat);
			Source_guess[3] = -Prat * (-Sigma_t * Fr0z + Sigma_a * velocity_z * diffTEr / Crat);




			 do the predict step 
			Errort[1] = pG->U[k][j][i].M1 + 0.5 * hdt * (Source[1] + Source_guess[1]) - hdt * (divFlux1[1] + divFlux2[1] + divFlux3[1]) - Uguess[1];
			Errort[2] = pG->U[k][j][i].M2 + 0.5 * hdt * (Source[2] + Source_guess[2]) - hdt * (divFlux1[2] + divFlux2[2] + divFlux3[2]) - Uguess[2];
			Errort[3] = pG->U[k][j][i].M3 + 0.5 * hdt * (Source[3] + Source_guess[3]) - hdt * (divFlux1[3] + divFlux2[3] + divFlux3[3]) - Uguess[3];



       		M1h = Uguess[1] + Source_Inv[1][1] * Errort[1];
			M2h = Uguess[2] + Source_Inv[2][2] * Errort[2];
			M3h = Uguess[3] + Source_Inv[3][3] * Errort[3];		
					
					*/

/*
        		Eh = pG->U[k][j][i].E
           			- q1*(rsf*x1Flux[k  ][j  ][i+1].E - lsf*x1Flux[k][j][i].E)
           			- q2*(    x2Flux[k  ][j+1][i  ].E -     x2Flux[k][j][i].E)
           			- q3*(    x3Flux[k+1][j  ][i  ].E -     x3Flux[k][j][i].E);

*/
					
					
					
/* Add shearing box source term */
#ifdef SHEARING_BOX
					if (ShearingBoxPot != NULL){
						cc_pos(pG,i,j,k,&x1,&x2,&x3);
						phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
						phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);
						M1h -= hdtodx1*(phir-phil)*pG->U[k][j][i].d;
						
						phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
						phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);
						M2h -= hdtodx2*(phir-phil)*pG->U[k][j][i].d;
						
						phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
						phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));
						M3h -= hdtodx3*(phir-phil)*pG->U[k][j][i].d;
					}
					
					M1h += pG->dt*Omega_0*pG->U[k][j][i].M2;
#ifdef FARGO
					M2h += hdt*(qshear-2.)*Omega_0*pG->U[k][j][i].M1;
#else
					M2h -= pG->dt*Omega_0*pG->U[k][j][i].M1;
#endif
#endif /* SHEARING_BOX */				



        		B1ch = 0.5*(    B1_x1Face[k][j][i] +     B1_x1Face[k  ][j  ][i+1]);
        		B2ch = 0.5*(    B2_x2Face[k][j][i] +     B2_x2Face[k  ][j+1][i  ]);
        		B3ch = 0.5*(    B3_x3Face[k][j][i] +     B3_x3Face[k+1][j  ][i  ]);
        		emf1_cc[k][j][i] = (B2ch*M3h - B3ch*M2h)/dhalf[k][j][i];
        		emf2_cc[k][j][i] = (B3ch*M1h - B1ch*M3h)/dhalf[k][j][i];
        		emf3_cc[k][j][i] = (B1ch*M2h - B2ch*M1h)/dhalf[k][j][i];

      			}
    		}
  	}

#endif

/*=== STEP 10: Compute 3D x1-Flux, x2-Flux, x3-Flux ===========================*/

/*--- Step 10a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

/*--- Step 10b ------------------------------------------------------------------
 * Compute 3D x1-fluxes from corrected L/R states.
 */

  	for (k=ks-1; k<=ke+1; k++) {
    		for (j=js-1; j<=je+1; j++) {
      			for (i=is; i<=ie+1; i++) {
#ifdef RADIATION_MHD
        		Bx = B1_x1Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x1Flux[k][j][i].d = dt;
			x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);
      			}
    		}
  	}

/*--- Step 10c ------------------------------------------------------------------
 * Compute 3D x2-fluxes from corrected L/R states.
 */

  	for (k=ks-1; k<=ke+1; k++) {
    		for (j=js; j<=je+1; j++) {
      			for (i=is-1; i<=ie+1; i++) {

#ifdef RADIATION_MHD
        		Bx = B2_x2Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[i],Wr[i],Bx, &x2Flux[k][j][i]);
      			}
    		}
  	}

/*--- Step 10d ------------------------------------------------------------------
 * Compute 3D x3-fluxes from corrected L/R states.
 */

  	for (k=ks; k<=ke+1; k++) {
    		for (j=js-1; j<=je+1; j++) {
      			for (i=is-1; i<=ie+1; i++) {

#ifdef RADIATION_MHD
        		Bx = B3_x3Face[k][j][i];
#endif
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x3Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x3Face[k][j][i],&Bx);

			/* Need parameter dt in radiation Riemann solver */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[i],Wr[i],Bx, &x3Flux[k][j][i]);
      			}
    		}
  	}


	
/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef RADIATION_MHD
  	integrate_emf1_corner(pG);
  	integrate_emf2_corner(pG);
  	integrate_emf3_corner(pG);
	
	
	
	/* Remap Ey at is and ie+1 to conserve Bz in shearing box */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
	
	/* compute remapped Ey from opposite side of grid */
	
    if (my_iproc == 0) {
		RemapEy_ix1(pD, emf2, remapEyiib);
    }
    if (my_iproc == (pD->NGrid[0]-1)) {
		RemapEy_ox1(pD, emf2, remapEyoib);
    }
	
	/* Now average Ey and remapped Ey */
	
    if (my_iproc == 0) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2[k][j][is]  = 0.5*(emf2[k][j][is] + remapEyiib[k][j]);
			}
		}
    }
	
    if (my_iproc == (pD->NGrid[0]-1)) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2[k][j][ie+1]  = 0.5*(emf2[k][j][ie+1] + remapEyoib[k][j]);
			}
		}
    }
#endif /* SHEARING_BOX */	
	

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */

  	for (k=ks; k<=ke; k++) {
    		for (j=js; j<=je; j++) {
      			for (i=is; i<=ie; i++) {

        		pG->B1i[k][j][i] += dtodx3*(emf2[k+1][j  ][i  ] - emf2[k][j][i]) -
                            		dtodx2*(emf3[k  ][j+1][i  ] - emf3[k][j][i]);
        		pG->B2i[k][j][i] += dtodx1*(emf3[k  ][j  ][i+1] - emf3[k][j][i]) -
                            		dtodx3*(emf1[k+1][j  ][i  ] - emf1[k][j][i]);

        		pG->B3i[k][j][i] += dtodx2*(    emf1[k  ][j+1][i  ] -     emf1[k][j][i]) -
                            		dtodx1*(emf2[k  ][j  ][i+1] - emf2[k][j][i]);
      			}

      		pG->B1i[k][j][ie+1] +=
        		dtodx3*(emf2[k+1][j  ][ie+1] - emf2[k][j][ie+1]) -
        		dtodx2*(emf3[k  ][j+1][ie+1] - emf3[k][j][ie+1]);
    		}

    		for (i=is; i<=ie; i++) {
      			pG->B2i[k][je+1][i] +=
        		dtodx1*(emf3[k  ][je+1][i+1] - emf3[k][je+1][i]) -
        		dtodx3*(emf1[k+1][je+1][i  ] - emf1[k][je+1][i]);
    		}
  	}

  	for (j=js; j<=je; j++) {
    		for (i=is; i<=ie; i++) {

      			pG->B3i[ke+1][j][i] += 
        			dtodx2*(    emf1[ke+1][j+1][i  ] - emf1[ke+1][j][i]) -
        			dtodx1*(    emf2[ke+1][j  ][i+1] - emf2[ke+1][j][i]);
    		}
  	}
#endif /* RADIATION_MHD */


#ifdef SELF_GRAVITY
/* Save mass fluxes in Grid structure for source term correction in main loop */
/* After we get the corrected fluxt */

  for (k=ks; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
        pG->x1MassFlux[k][j][i] = x1Flux[k][j][i].d;
        pG->x2MassFlux[k][j][i] = x2Flux[k][j][i].d;
        pG->x3MassFlux[k][j][i] = x3Flux[k][j][i].d;
      }
    }
  }


#ifdef CONS_GRAVITY

/* Need to update the density first so that we can calculate the new potential */
/* There is no source term for density from radiation */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
		/* We need the density at time step n, for radiation source term and  */
		/* and self-gravity energy flux term */
		density_old[k][j][i] = pG->U[k][j][i].d;
		pG->U[k][j][i].d  -= dtodx1*(x1Flux[k][j][i+1].d -x1Flux[k][j][i].d );
		pG->U[k][j][i].d  -= dtodx2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
		pG->U[k][j][i].d  -= dtodx3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );

      }
    }
  }
	bvals_mhd(pD);

/* With the updated flux, now we can calculate dphidt with Poisson Solver */
    	(*SelfGrav_cons)(pD);
	/* Need to apply boundary condition for the new calculated dphidt */
    	bvals_grav(pD);
/* Now calculate the energy flux */
	for(k=ks; k<=ke+1; k++){
  		for(j=js; j<=je+1;j++){	
			for (i=is; i<=ie+1; i++) {
				phil = 0.25*(pG->Phi[k][j][i-1]+pG->Phi_old[k][j][i]+pG->Phi_old[k][j][i-1]+pG->Phi[k][j][i]);
				gxl = 0.5 * (pG->Phi[k][j][i-1] + pG->Phi_old[k][j][i-1]  - pG->Phi[k][j][i  ] - pG->Phi_old[k][j][i  ])/(pG->dx1);
				dotphil  = 0.5*(pG->dphidt[k][j][i-1] + pG->dphidt[k][j][i  ]);		
				dotgxl = (pG->dphidt[k][j][i-1] - pG->dphidt[k][j][i  ])/(pG->dx1);

				x1Flux_grav[k][j][i] =-0.5*(phil*dotgxl-dotphil*gxl)/four_pi_G + x1Flux[k][j][i].d*phil;

				phil = 0.25*(pG->Phi[k][j-1][i]+pG->Phi_old[k][j-1][i]+pG->Phi_old[k][j][i]+pG->Phi[k][j][i]);
				gxl = 0.5 * (pG->Phi[k][j-1][i] + pG->Phi_old[k][j-1][i]  - pG->Phi[k][j][i  ] - pG->Phi_old[k][j][i  ])/(pG->dx2);
				dotphil  = 0.5*(pG->dphidt[k][j-1][i] + pG->dphidt[k][j][i  ]);		
				dotgxl = (pG->dphidt[k][j-1][i] - pG->dphidt[k][j][i  ])/(pG->dx2);

				x2Flux_grav[k][j][i] =-0.5*(phil*dotgxl-dotphil*gxl)/four_pi_G + x2Flux[k][j][i].d*phil;

				phil = 0.25*(pG->Phi[k-1][j][i]+pG->Phi_old[k-1][j][i]+pG->Phi_old[k][j][i]+pG->Phi[k][j][i]);
				gxl = 0.5 * (pG->Phi[k-1][j][i] + pG->Phi_old[k-1][j][i]  - pG->Phi[k][j][i  ] - pG->Phi_old[k][j][i  ])/(pG->dx3);
				dotphil  = 0.5*(pG->dphidt[k-1][j][i] + pG->dphidt[k][j][i  ]);		
				dotgxl = (pG->dphidt[k-1][j][i] - pG->dphidt[k][j][i  ])/(pG->dx3);

				x3Flux_grav[k][j][i] =-0.5*(phil*dotgxl-dotphil*gxl)/four_pi_G + x3Flux[k][j][i].d*phil;		
			}
 		}
	}
	
	/*-----------Now swap, desntiy_old actually save the new density  */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
		/* We need the density at time step n, for radiation source term and  */
		/* and self-gravity energy flux term */
		/* We do not need to update ghost zone here */
		Tempswap = pG->U[k][j][i].d;
		pG->U[k][j][i].d = density_old[k][j][i];
		density_old[k][j][i] = Tempswap;
      }
    }
  }


#endif

#endif /* SELF_GRAVITY */


/*-------Step 11: Predict and correct step after we get the flux------------- */
	
		
/*----------Radiation quantities are not updated in the modified Godunov corrector step */
/*-----------NOTE that x1flux, x2flux, x3flux are flux from Remann problem. If there is extra-----* 
 *----------source terms, flux due to those source terms should also be added in the-------*
 *---------- modified Godunov method -------------------------------------*/

	for (k=ks; k<=ke; k++){
		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++) {

			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			/* Load 1D vector */
			Usource.d  = pG->U[k][j][i].d;
      			Usource.Mx = pG->U[k][j][i].M1;
      			Usource.My = pG->U[k][j][i].M2;
      			Usource.Mz = pG->U[k][j][i].M3;
      			Usource.E  = pG->U[k][j][i].E;
			Usource.Er  = pG->U[k][j][i].Er;
    			Usource.Fr1  = pG->U[k][j][i].Fr1;
    			Usource.Fr2  = pG->U[k][j][i].Fr2;
    			Usource.Fr3  = pG->U[k][j][i].Fr3;
			Usource.Edd_11  = pG->U[k][j][i].Edd_11;
			Usource.Edd_21  = pG->U[k][j][i].Edd_21;
			Usource.Edd_22  = pG->U[k][j][i].Edd_22;
			Usource.Edd_31	= pG->U[k][j][i].Edd_31;
			Usource.Edd_32	= pG->U[k][j][i].Edd_32;
			Usource.Edd_33	= pG->U[k][j][i].Edd_33;
			for(m=0;m<NOPACITY;m++){
				Usource.Sigma[m] = pG->U[k][j][i].Sigma[m];
			}


#ifdef RADIATION_MHD
      			Usource.By = pG->U[k][j][i].B2c;
      			Usource.Bz = pG->U[k][j][i].B3c;
			Bx = pG->U[k][j][i].B1c;
#else
			Bx = 0.0;
#endif /* MHD */
			
			/* Initialize source to be zero */
			/* Because this is also used in other place */

		
#ifdef SHEARING_BOX
			for(m=0; m<4; m++)
				ShearSource[m] = 0.0;
#endif




			density = Usource.d;
			velocity_x = Usource.Mx / density;
			velocity_y = Usource.My / density;
			velocity_z = Usource.Mz / density;

			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
				
			pressure = (pG->U[k][j][i].E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
			pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
			temperature = pressure / (density * R_ideal);

			Sigma_sF = Usource.Sigma[0];
			Sigma_aF = Usource.Sigma[1];
			Sigma_aP = Usource.Sigma[2];
			Sigma_aE = Usource.Sigma[3];

		
			diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;
			

					
#ifdef FARGO
			
			/* Include background shearing in Usource, which is only used in dSource */

			velocity_y -= qom * x1;						
					
#endif			
			
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);
			

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	

			SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);

			SFmz =(Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/

			Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
			Source_Inv[4][1] = (-dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[1][1];
			Source_Inv[4][2] = (-dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[2][2];
			Source_Inv[4][3] = (-dt * Prat * Crat * SEmz/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[3][3];
			Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	


			/*--------Calculate the guess solution-------------*/
	
			divFlux1[0] = (x1Flux[k][j][i+1].d  - x1Flux[k][j][i].d ) / dx1;
			divFlux1[1] = (x1Flux[k][j][i+1].Mx - x1Flux[k][j][i].Mx) / dx1;
			divFlux1[2] = (x1Flux[k][j][i+1].My - x1Flux[k][j][i].My) / dx1;
			divFlux1[3] = (x1Flux[k][j][i+1].Mz - x1Flux[k][j][i].Mz) / dx1;
			divFlux1[4] = (x1Flux[k][j][i+1].E  - x1Flux[k][j][i].E ) / dx1; 

			divFlux2[0] = (x2Flux[k][j+1][i].d  - x2Flux[k][j][i].d ) / dx2;
			divFlux2[1] = (x2Flux[k][j+1][i].Mz - x2Flux[k][j][i].Mz) / dx2;
			divFlux2[2] = (x2Flux[k][j+1][i].Mx - x2Flux[k][j][i].Mx) / dx2;
			divFlux2[3] = (x2Flux[k][j+1][i].My - x2Flux[k][j][i].My) / dx2;
			divFlux2[4] = (x2Flux[k][j+1][i].E  - x2Flux[k][j][i].E ) / dx2; 

			divFlux3[0] = (x3Flux[k+1][j][i].d  - x3Flux[k][j][i].d ) / dx3;
			divFlux3[1] = (x3Flux[k+1][j][i].My - x3Flux[k][j][i].My) / dx3;
			divFlux3[2] = (x3Flux[k+1][j][i].Mz - x3Flux[k][j][i].Mz) / dx3;
			divFlux3[3] = (x3Flux[k+1][j][i].Mx - x3Flux[k][j][i].Mx) / dx3;
			divFlux3[4] = (x3Flux[k+1][j][i].E  - x3Flux[k][j][i].E ) / dx3; 

	/*------Flux due to static gravitational potential ---------*/

		if (StaticGravPot != NULL){
    			
        		phic = (*StaticGravPot)((x1            ),x2,x3);
        		phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        		phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

			/* x direction */

			SourceFlux[1] = dhalf[k][j][i]*(phir-phil)/dx1;
        		SourceFlux[4]= (x1Flux[k][j][i  ].d*(phic - phil) +
				 	x1Flux[k][j][i+1].d*(phir - phic))/dx1;

			divFlux1[1] += SourceFlux[1];
			divFlux1[4] += SourceFlux[4];

			/* y direction */
			
        		phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        		phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        		SourceFlux[2] = dhalf[k][j][i]*(phir-phil)/dx2;
        		SourceFlux[4] = (x2Flux[k][j  ][i].d*(phic - phil) +
                                	     x2Flux[k][j+1][i].d*(phir - phic))/dx2;
			divFlux2[2] += SourceFlux[2];
			divFlux2[4] += SourceFlux[4];      		


			/* z direction */
			phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
          		phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

			SourceFlux[3] = dhalf[k][j][i]*(phir-phil)/dx3;
        		SourceFlux[4] = (x3Flux[k][j  ][i].d*(phic - phil) +
                                	     x3Flux[k + 1][j][i].d*(phir - phic))/dx3;
			divFlux3[3] += SourceFlux[3];
			divFlux3[4] += SourceFlux[4];   
	
  		}

	/*--------Flux due to self-gravity ------------------------*/
#ifdef SELF_GRAVITY
		/************************************************/
		/*------------x direction ----------------------*/

		phic = pG->Phi[k][j][i];
        	phil = 0.5*(pG->Phi[k][j][i-1] + pG->Phi[k][j][i  ]);
        	phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+1]);

/* gx, gy and gz centered at L and R x1-faces */
        	gxl = (pG->Phi[k][j][i-1] - pG->Phi[k][j][i  ])*(dx1i);
        	gxr = (pG->Phi[k][j][i  ] - pG->Phi[k][j][i+1])*(dx1i);

        	gyl = 0.25*((pG->Phi[k][j-1][i-1] - pG->Phi[k][j+1][i-1]) +
                	    (pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]) )*(dx2i);
        	gyr = 0.25*((pG->Phi[k][j-1][i  ] - pG->Phi[k][j+1][i  ]) +
                	    (pG->Phi[k][j-1][i+1] - pG->Phi[k][j+1][i+1]) )*(dx2i);

        	gzl = 0.25*((pG->Phi[k-1][j][i-1] - pG->Phi[k+1][j][i-1]) +
                	    (pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) )*(dx3i);
        	gzr = 0.25*((pG->Phi[k-1][j][i  ] - pG->Phi[k+1][j][i  ]) +
                	    (pG->Phi[k-1][j][i+1] - pG->Phi[k+1][j][i+1]) )*(dx3i);

/* momentum fluxes in x1.  2nd term is needed only if Jean's swindle used */
        	flx_m1l = 0.5*(gxl*gxl-gyl*gyl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        	flx_m1r = 0.5*(gxr*gxr-gyr*gyr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        	flx_m2l = gxl*gyl/four_pi_G;
        	flx_m2r = gxr*gyr/four_pi_G;

        	flx_m3l = gxl*gzl/four_pi_G;
        	flx_m3r = gxr*gzr/four_pi_G;

/* Update momenta and energy with d/dx1 terms  */
        	
		SourceFlux[1] = (flx_m1r - flx_m1l)/dx1;
		divFlux1[1] += SourceFlux[1];
	
		SourceFlux[2] = (flx_m2r - flx_m2l)/dx1;
		divFlux1[2] += SourceFlux[2];
		
		SourceFlux[3] = (flx_m3r - flx_m3l)/dx1;
		divFlux1[3] += SourceFlux[3];

#ifndef CONS_GRAVITY
		SourceFlux[4] = (x1Flux[k][j][i  ].d*(phic - phil) +
                	      	x1Flux[k][j][i+1].d*(phir - phic)) / dx1;
		divFlux1[4] += SourceFlux[4];
#endif

		/************************************************/
		/*------------y direction ----------------------*/

  		phil = 0.5*(pG->Phi[k][j-1][i] + pG->Phi[k][j  ][i]);
        	phir = 0.5*(pG->Phi[k][j  ][i] + pG->Phi[k][j+1][i]);

/* gx, gy and gz centered at L and R x2-faces */
        	gxl = 0.25*((pG->Phi[k][j-1][i-1] - pG->Phi[k][j-1][i+1]) +
                	    (pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]) )*(dx1i);
        	gxr = 0.25*((pG->Phi[k][j  ][i-1] - pG->Phi[k][j  ][i+1]) +
                	    (pG->Phi[k][j+1][i-1] - pG->Phi[k][j+1][i+1]) )*(dx1i);

        	gyl = (pG->Phi[k][j-1][i] - pG->Phi[k][j  ][i])*(dx2i);
        	gyr = (pG->Phi[k][j  ][i] - pG->Phi[k][j+1][i])*(dx2i);

        	gzl = 0.25*((pG->Phi[k-1][j-1][i] - pG->Phi[k+1][j-1][i]) +
                    	(pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]) )*(dx3i);
        	gzr = 0.25*((pG->Phi[k-1][j  ][i] - pG->Phi[k+1][j  ][i]) +
                	    (pG->Phi[k-1][j+1][i] - pG->Phi[k+1][j+1][i]) )*(dx3i);

/* momentum fluxes in x2.  2nd term is needed only if Jean's swindle used */
        	flx_m1l = gyl*gxl/four_pi_G;
        	flx_m1r = gyr*gxr/four_pi_G;

        	flx_m2l = 0.5*(gyl*gyl-gxl*gxl-gzl*gzl)/four_pi_G + grav_mean_rho*phil;
        	flx_m2r = 0.5*(gyr*gyr-gxr*gxr-gzr*gzr)/four_pi_G + grav_mean_rho*phir;

        	flx_m3l = gyl*gzl/four_pi_G;
        	flx_m3r = gyr*gzr/four_pi_G;


/* Update momenta and energy with d/dx1 terms  */
        	
		SourceFlux[1] = (flx_m1r - flx_m1l)/dx2;
		divFlux2[1] += SourceFlux[1];
	
		SourceFlux[2] = (flx_m2r - flx_m2l)/dx2;
		divFlux2[2] += SourceFlux[2];
		
		SourceFlux[3] = (flx_m3r - flx_m3l)/dx2;
		divFlux2[3] += SourceFlux[3];

#ifndef CONS_GRAVITY
		SourceFlux[4] = (x2Flux[k][j  ][i].d*(phic - phil) +
                                    x2Flux[k][j+1][i].d*(phir - phic)) / dx2;
		divFlux2[4] += SourceFlux[4];
#endif

		/************************************************/
		/*------------z direction ----------------------*/
		
      		phil = 0.5*(pG->Phi[k-1][j][i] + pG->Phi[k  ][j][i]);
        	phir = 0.5*(pG->Phi[k  ][j][i] + pG->Phi[k+1][j][i]);

/* gx, gy and gz centered at L and R x3-faces */
        	gxl = 0.25*((pG->Phi[k-1][j][i-1] - pG->Phi[k-1][j][i+1]) +
                	    (pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]) )*(dx1i);
        	gxr = 0.25*((pG->Phi[k  ][j][i-1] - pG->Phi[k  ][j][i+1]) +
                	    (pG->Phi[k+1][j][i-1] - pG->Phi[k+1][j][i+1]) )*(dx1i);

        	gyl = 0.25*((pG->Phi[k-1][j-1][i] - pG->Phi[k-1][j+1][i]) +
                	    (pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]) )*(dx2i);
        	gyr = 0.25*((pG->Phi[k  ][j-1][i] - pG->Phi[k  ][j+1][i]) +
                	    (pG->Phi[k+1][j-1][i] - pG->Phi[k+1][j+1][i]) )*(dx2i);

        	gzl = (pG->Phi[k-1][j][i] - pG->Phi[k  ][j][i])*(dx3i);
        	gzr = (pG->Phi[k  ][j][i] - pG->Phi[k+1][j][i])*(dx3i);

/* momentum fluxes in x3.  2nd term is needed only if Jean's swindle used */
        	flx_m1l = gzl*gxl/four_pi_G;
        	flx_m1r = gzr*gxr/four_pi_G;

        	flx_m2l = gzl*gyl/four_pi_G;
        	flx_m2r = gzr*gyr/four_pi_G;

        	flx_m3l = 0.5*(gzl*gzl-gxl*gxl-gyl*gyl)/four_pi_G + grav_mean_rho*phil;
        	flx_m3r = 0.5*(gzr*gzr-gxr*gxr-gyr*gyr)/four_pi_G + grav_mean_rho*phir;

/* Update momenta and energy with d/dx1 terms  */
        	
		SourceFlux[1] = (flx_m1r - flx_m1l)/dx3;
		divFlux3[1] += SourceFlux[1];
	
		SourceFlux[2] = (flx_m2r - flx_m2l)/dx3;
		divFlux3[2] += SourceFlux[2];
		
		SourceFlux[3] = (flx_m3r - flx_m3l)/dx3;
		divFlux3[3] += SourceFlux[3];

#ifndef CONS_GRAVITY
		SourceFlux[4] = (x3Flux[k  ][j][i].d*(phic - phil) +
                                    x3Flux[k+1][j][i].d*(phir - phic)) / dx3;
		divFlux3[4] += SourceFlux[4];
#endif

	/* Now the energy flux due to self-gravity in cons-gravity */

#ifdef CONS_GRAVITY
		divFlux1[4] += (x1Flux_grav[k][j][i+1] - x1Flux_grav[k][j][i]) / dx1;
		divFlux2[4] += (x2Flux_grav[k][j+1][i] - x2Flux_grav[k][j][i]) / dx2;
		divFlux3[4] += (x3Flux_grav[k+1][j][i] - x3Flux_grav[k][j][i]) / dx3;

#endif

#endif /* SELF_GRAVITY */


	/*================== End static gravitational flux ================*/

	/* cacluate guess solution */

		for(n=0; n<5; n++) {
			tempguess[n] = 0.0;
			for(m=0; m<5; m++) {
				tempguess[n] += dt * Source_Inv[n][m] * (Source[k][j][i][m] - divFlux1[m] - divFlux2[m] - divFlux3[m]);
			}
		}

		Uguess[0] = pG->U[k][j][i].d  + tempguess[0];
		Uguess[1] = pG->U[k][j][i].M1 + tempguess[1];
		Uguess[2] = pG->U[k][j][i].M2 + tempguess[2];
		Uguess[3] = pG->U[k][j][i].M3 + tempguess[3];
		Uguess[4] = pG->U[k][j][i].E  + tempguess[4];

#ifdef CONS_GRAVITY
		/* density_old is now actually the updated density */
		Uguess[4] += 0.5*(pG->U[k][j][i].d-grav_mean_rho)*pG->Phi_old[k][j][i]-0.5*(density_old[k][j][i]-grav_mean_rho)*pG->Phi[k][j][i];

#endif

		/*  Uguess[0] = d; Uguess[1]=Mx; Uguess[2]=My; Uguess[3]=Mz, Uguess[4]=E */
	
		/* Now calculate the source term due to the guess solution */
		/* Flux is not changed during the correction step */

/*================================================================*/
/* Guess solution also needs to add shearing_box source term */
/* Just use forward Euler to add the shearing source term in guess solution */
					
#ifdef SHEARING_BOX
		phic = (*ShearingBoxPot)((x1            ),x2,x3);
		phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
		phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

		ShearSource[3] = dtodx1*(x1Flux[k][j][i  ].d*(phic - phil) +
				x1Flux[k][j][i+1].d*(phir - phic));
					
		phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
		phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);

		ShearSource[3] -= dtodx2*(x2Flux[k][j  ][i].d*(phic - phil) +
				x2Flux[k][j+1][i].d*(phir - phic));

							
		phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
		phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));

		ShearSource[3] -= dtodx3*(x3Flux[k  ][j][i].d*(phic - phil) +
				x3Flux[k+1][j][i].d*(phir - phic));

		/* Use Crank Nichson but do not include radiation source term for guess solution */
		/* Store the current state */
		M1n  = pG->U[k][j][i].M1;
#ifdef FARGO
		dM2n = pG->U[k][j][i].M2;
#else
		dM2n = pG->U[k][j][i].M2 + qom*x1*pG->U[k][j][i].d;
#endif
					

#ifndef FARGO
		frx1_dM2 = qom*(x1+0.5*pG->dx1)*x1Flux[k][j][i+1].d;
		flx1_dM2 = qom*(x1-0.5*pG->dx1)*x1Flux[k][j][i  ].d;
		frx2_dM2 = qom*(x1            )*x2Flux[k][j+1][i].d;
		flx2_dM2 = qom*(x1            )*x2Flux[k][j  ][i].d;
		frx3_dM2 = qom*(x1            )*x3Flux[k+1][j][i].d;
		flx3_dM2 = qom*(x1            )*x3Flux[k  ][j][i].d;
#endif
/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
		M1e  = M1n + 0.5 * tempguess[1];
		dM2e = dM2n + 0.5 * tempguess[2];
#ifndef FARGO
		dM2e -= hdtodx1*(frx1_dM2 - flx1_dM2)
			- hdtodx2*(frx2_dM2 - flx2_dM2) 
			- hdtodx3*(frx3_dM2 - flx3_dM2);
#endif

		ShearSource[0] = (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact; 
		ShearSource[1] = 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;
#ifndef FARGO
		ShearSource[1] -= 0.5*qshear*om_dt*(x1Flux[k][j][i].d + x1Flux[k][j][i+1].d);
#endif		
		/* Now add the source term to guess solution */


		Uguess[1] += ShearSource[0];
		Uguess[2] += ShearSource[1];
		Uguess[3] += ShearSource[2];
		Uguess[4] += ShearSource[3];
				
#endif
/* End SHEARING_BOX */
/*==========================================================*/
		
		density    = Uguess[0];
		velocity_x = Uguess[1] / density;
		velocity_y = Uguess[2] / density;
		velocity_z = Uguess[3] / density;

		velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);



		pressure = (Uguess[4] - 0.5 * (density * velocity * velocity)) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif
		temperature = pressure / (density * R_ideal);
		
		if(pressure < TINY_NUMBER){
			pressure = density * R_ideal * pow(pG->U[k][j][i].Er, 0.25);
			Uguess[4] = 0.5 * (density * velocity * velocity) + pressure / (Gamma - 1.0);
#ifdef RADIATION_MHD
			Uguess[4] += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
		}

		temperature = pressure / (density * R_ideal);

		/* If Opacity is not set, Sigma_? will not be changed. */

		if(Opacity != NULL){
			Opacity(density,temperature, Sigma, NULL);
		
			Sigma_sF = Sigma[0];
			Sigma_aF = Sigma[1];
			Sigma_aP = Sigma[2];
			Sigma_aE = Sigma[3];
		}
		else{
			Sigma[0] = Sigma_sF;
			Sigma[1] = Sigma_aF;
			Sigma[2] = Sigma_aP;
			Sigma[3] = Sigma_aE;
		}

		
	

		diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;
		
		

		/* update source term */
		Usource.d  = Uguess[0];
		Usource.Mx = Uguess[1];
		Usource.My = Uguess[2];
		Usource.Mz = Uguess[3];
		Usource.E  = Uguess[4];

		for(m=0; m<NOPACITY;m++)
			Usource.Sigma[m] = Sigma[m];

					
#ifdef FARGO
		/* With FARGO, we should add background shearing to the source terms */
					
		/* Include background shearing in Usource, which is only used in dSource */
					
		velocity_y -= qom * x1;	
					
#endif	

		/* The Source term */
		
		dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);

		
	
		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	

			SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);

			SFmz = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat) 
				+ diffTEr / (density * Crat);	


			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);

		/*=========================================================*/

			Source_Inv[4][0] = -dt * Prat * Crat * SErho/(1.0 + dt * Prat * Crat * SEE);
			Source_Inv[4][1] = (-dt * Prat * Crat * SEmx/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[1][1];
			Source_Inv[4][2] = (-dt * Prat * Crat * SEmy/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[2][2];
			Source_Inv[4][3] = (-dt * Prat * Crat * SEmz/(1.0 + dt * Prat * Crat * SEE)) * Source_Inv[3][3];
			Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	


			/* co-moving flux */
			Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
			Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
			Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;

			/* Source term for momentum */
			Source_guess[1] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0x + velocity_x * diffTEr / Crat);
			Source_guess[2] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0y + velocity_y * diffTEr / Crat);
			Source_guess[3] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0z + velocity_z * diffTEr / Crat);
			
			/* Source term for total Energy */
			Source_guess[4] = -Prat * Crat * (diffTEr + (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat);


			/* Calculate the error term */
			/* Error term should include shearing source term as Guess solution include shearing source term */
			Errort[0] = pG->U[k][j][i].d  + hdt * (Source[k][j][i][0] + Source_guess[0]) 
					- dt * (divFlux1[0] + divFlux2[0] + divFlux3[0]) - Uguess[0];
			Errort[1] = pG->U[k][j][i].M1 + hdt * (Source[k][j][i][1] + Source_guess[1]) 
					- dt * (divFlux1[1] + divFlux2[1] + divFlux3[1]) - Uguess[1];
			Errort[2] = pG->U[k][j][i].M2 + hdt * (Source[k][j][i][2] + Source_guess[2]) 
					- dt * (divFlux1[2] + divFlux2[2] + divFlux3[2]) - Uguess[2];
			Errort[3] = pG->U[k][j][i].M3 + hdt * (Source[k][j][i][3] + Source_guess[3]) 
					- dt * (divFlux1[3] + divFlux2[3] + divFlux3[3]) - Uguess[3];
			Errort[4] = pG->U[k][j][i].E  + hdt * (Source[k][j][i][4] + Source_guess[4]) 
					- dt * (divFlux1[4] + divFlux2[4] + divFlux3[4]) - Uguess[4];

#ifdef CONS_GRAVITY
			Errort[4] += 0.5*(pG->U[k][j][i].d-grav_mean_rho)*pG->Phi_old[k][j][i]-0.5*(density_old[k][j][i]-grav_mean_rho)*pG->Phi[k][j][i];
#endif

		/* substract shearing source term */
#ifdef SHEARING_BOX
			Errort[1] += ShearSource[0];
			Errort[2] += ShearSource[1];
			Errort[4] += ShearSource[3];
#endif


			/* Calculate the correction */
			for(n=0; n<5; n++) {
				tempguess[n] = 0.0;
				for(m=0; m<5; m++) {
					tempguess[n] += Source_Inv[n][m] * Errort[m];
				}
			}

				
			/* Apply the correction */
					
		
			pG->U[k][j][i].d  = Uguess[0] + tempguess[0];
			pG->U[k][j][i].M1 = Uguess[1] + tempguess[1];
			pG->U[k][j][i].M2 = Uguess[2] + tempguess[2];
			pG->U[k][j][i].M3 = Uguess[3] + tempguess[3];
			pG->U[k][j][i].E  = Uguess[4] + tempguess[4];



			} /* End i */
		}/* End j */
	}/* End k */



/*--- Step 12 -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef RADIATION_MHD
  	for (k=ks; k<=ke; k++) {
    		for (j=js; j<=je; j++) {
      			for (i=is; i<=ie; i++) {

        			pG->U[k][j][i].B1c = 0.5*(    pG->B1i[k][j][i] +     pG->B1i[k][j][i+1]);
        			pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
        			pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
					
      			}
    		}
  	}
#endif /* RADIATION MHD */


		
	/* Boundary condition is applied in the main function */
	/* Check the pressure to make sure that it is positive */
	/* This will make the code more robust */
/*
	for (k=ks; k<=ke; k++){
		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++){
				if(pG->U[k][j][i].d < TINY_NUMBER)
					pG->U[k][j][i].d = 1.e-3;

				density = pG->U[k][j][i].d;

				pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
				+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density ) * (Gamma - 1);
		
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

				if(pressure < TINY_NUMBER){
					if(pG->U[k][j][i].Er < TINY_NUMBER)
						pG->U[k][j][i].Er = TINY_NUMBER;
					

				pressure = density * R_ideal * pow(pG->U[k][j][i].Er,0.25);
				
				pG->U[k][j][i].E = 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density + pressure / (Gamma - 1.0);
#ifdef RADIATION_MHD
				pG->U[k][j][i].E += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif	
				} 
				

			} 
		}
	}
*/


	/* Update the opacity if Opacity function is set in the problem generator */
	if(Opacity != NULL){
		for (k=ks; k<=ke; k++){
			for (j=js; j<=je; j++) {
    				for (i=is; i<=ie; i++){
				
				density = pG->U[k][j][i].d;
				
				pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
				+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density ) * (Gamma - 1);
			/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

				if(pressure > TINY_NUMBER)
				{
					temperature = pressure / (density * R_ideal);

						
					Opacity(density,temperature,Sigma,NULL);
					for(m=0;m<NOPACITY;m++){
						pG->U[k][j][i].Sigma[m] = Sigma[m];
					}

				}

			
				}
			}
		}
	}


  return;
}


/*----------------------------------------------------------------------------*/
/* integrate_init_3d: Allocate temporary integration arrays 
*/

void integrate_init_3d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

#ifdef RADIATION_MHD
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;

  if ((emf1_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* MHD */



  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef RADIATION_MHD
  if ((B1_x1Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif /* MHD */

  if ((U1d=(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul =(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur =(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((x1Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;


#ifdef CONS_GRAVITY
  if ((x1Flux_grav   =(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((x2Flux_grav   =(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((x3Flux_grav   =(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    == NULL) goto on_error;
  if ((density_old   =(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))
    ==NULL)  goto on_error;
#endif


  if ((dhalf = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  
	/* For source term */
  if ((Source = (Real****)calloc_4d_array(size3, size2, size1,5, sizeof(Real))) == NULL)
    goto on_error;

  if ((Alpha = (Real***)calloc_3d_array(size3, size2, size1,sizeof(Real))) == NULL)
    goto on_error;

  if ((Beta = (Real****)calloc_4d_array(size3, size2, size1,3,sizeof(Real))) == NULL)
    goto on_error;

	
#ifdef SHEARING_BOX
	if ((remapEyiib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
		goto on_error;
	if ((remapEyoib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
		goto on_error;
#endif /* End SHEARING_BOX */

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_3d:  Free temporary integration arrays 
 */

void integrate_destruct_3d(void)
{

#ifdef RADIATION_MHD
  if (emf1    != NULL) free_3d_array(emf1);
  if (emf2    != NULL) free_3d_array(emf2);
  if (emf3    != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
#endif /* MHD */



  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef RADIATION_MHD
  if (B1_x1Face != NULL) free_3d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array(B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array(B3_x3Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_3d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_3d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_3d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_3d_array(Ur_x2Face);
  if (Ul_x3Face != NULL) free_3d_array(Ul_x3Face);
  if (Ur_x3Face != NULL) free_3d_array(Ur_x3Face);
  if (x1Flux    != NULL) free_3d_array(x1Flux);
  if (x2Flux    != NULL) free_3d_array(x2Flux);
  if (x3Flux    != NULL) free_3d_array(x3Flux);
  if (dhalf     != NULL) free_3d_array(dhalf);
  if (phalf     != NULL) free_3d_array(phalf);

  if(Source	!= NULL) free_4d_array(Source);

  if(Alpha	!= NULL) free_3d_array(Alpha);

  if(Beta	!= NULL) free_4d_array(Beta);

#ifdef SHEARING_BOX
	if (remapEyiib != NULL) free_2d_array(remapEyiib);
	if (remapEyoib != NULL) free_2d_array(remapEyoib);
#endif

#ifdef CONS_GRAVITY
  if (x1Flux_grav    != NULL) free_3d_array(x1Flux_grav);
  if (x2Flux_grav    != NULL) free_3d_array(x2Flux_grav);
  if (x3Flux_grav    != NULL) free_3d_array(x3Flux_grav);
  if (density_old    != NULL) free_3d_array(density_old);
#endif


  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf1_corner
 * integrate_emf2_corner
 * integrate_emf3_corner
 *   Integrates face centered B-fluxes to compute corner EMFs.  Note:
 *   x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *   x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *   x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *   x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *   x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *   x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */

#ifdef RADIATION_MHD
static void integrate_emf1_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

static void integrate_emf2_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

static void integrate_emf3_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de3_l1, de3_r1, de3_l2, de3_r2;
  Real rsf=1.0,lsf=1.0;

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
#ifdef CYLINDRICAL
        rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}
#endif /* RADIATION MHD */


/*========================================*/
/* Private function to calculate the source term */
/* We only need to calculate the source term once */
/* This is source term for momentum and energy */
/* When used for velocity and pressure, should change */
void updatesource(GridS *pG)
{

	int i,il,iu,is=pG->is, ie=pG->ie;
  	int j,jl,ju,js=pG->js, je=pG->je;
  	int k, kl, ku, ks=pG->ks, ke=pG->ke, m;
	double x1, x2, x3;
	double dt = pG->dt;

	double density, velocity_x, velocity_y, velocity_z, velocity, pressure, temperature;
	double Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE, Bx, Fr0x, Fr0y, Fr0z;
	Real SPP, dSigmadP[4], diffTEr, diffTErdP;
	Real dSigma[8];
#ifdef SHEARING_BOX
	Real qom = qshear*Omega_0;
#endif
	Cons1DS Usource;
	/* for source term */

	/* In case momentum becomes stiff */
	Real SVVx, SVVy, SVVz, betax, betay, betaz, alpha;
 
	il = is - 2;
  	iu = ie + 2;
  	jl = js - 2;
  	ju = je + 2;
  	kl = ks - 2;
  	ku = ke + 2;


	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){
			for(i=il; i<=iu; i++){
	
			cc_pos(pG,i,j,k,&x1,&x2,&x3);

	/* Load 1D vector */
			Usource.d  = pG->U[k][j][i].d;
      			Usource.Mx = pG->U[k][j][i].M1;
      			Usource.My = pG->U[k][j][i].M2;
      			Usource.Mz = pG->U[k][j][i].M3;
      			Usource.E  = pG->U[k][j][i].E;
			Usource.Er  = pG->U[k][j][i].Er;
    			Usource.Fr1  = pG->U[k][j][i].Fr1;
    			Usource.Fr2  = pG->U[k][j][i].Fr2;
    			Usource.Fr3  = pG->U[k][j][i].Fr3;
			Usource.Edd_11  = pG->U[k][j][i].Edd_11;
			Usource.Edd_21  = pG->U[k][j][i].Edd_21;
			Usource.Edd_22  = pG->U[k][j][i].Edd_22;
			Usource.Edd_31	= pG->U[k][j][i].Edd_31;
			Usource.Edd_32	= pG->U[k][j][i].Edd_32;
			Usource.Edd_33	= pG->U[k][j][i].Edd_33;
			for(m=0;m<NOPACITY;m++){
				Usource.Sigma[m] = pG->U[k][j][i].Sigma[m];
			}


#ifdef RADIATION_MHD
      			Usource.By = pG->U[k][j][i].B2c;
      			Usource.Bz = pG->U[k][j][i].B3c;
			Bx = pG->U[k][j][i].B1c;
#else
			Bx = 0.0;
#endif /* MHD */
			

			density = Usource.d;
			velocity_x = Usource.Mx / density;
			velocity_y = Usource.My / density;
			velocity_z = Usource.Mz / density;

			velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
				
			pressure = (pG->U[k][j][i].E - 0.5 * density * velocity * velocity) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
			pressure -= 0.5 * (Bx * Bx + Usource.By * Usource.By + Usource.Bz * Usource.Bz) * (Gamma - 1.0);
#endif
			/* in case pressure is negative */

			if(pressure > TINY_NUMBER){

				temperature = pressure / (density * R_ideal);

				Sigma_sF = pG->U[k][j][i].Sigma[0];
				Sigma_aF = pG->U[k][j][i].Sigma[1];
				Sigma_aP = pG->U[k][j][i].Sigma[2];
				Sigma_aE = pG->U[k][j][i].Sigma[3];


				/* Opacity is already included in source term */
				diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;


		
					
#ifdef FARGO
					/* With FARGO, we should add background shearing to the source terms */
				/* Include background shearing in Usource, which is only used in dSource */

				velocity_y -= qom * x1;						
					
#endif			

				/* co-moving flux */
				Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
				Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
				Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;

			/* Source term for momentum, not velocity*/
				Source[k][j][i][0] = 0.0;
				Source[k][j][i][1] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0x + velocity_x * diffTEr / Crat);
				Source[k][j][i][2] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0y + velocity_y * diffTEr / Crat);
				Source[k][j][i][3] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0z + velocity_z * diffTEr / Crat);
			
			/* Source term for energy */
				Source[k][j][i][4] = -Prat * Crat * (diffTEr + (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x + velocity_y * Fr0y * velocity_z * Fr0z)/Crat);

			

				if(Opacity != NULL){
					 Opacity(density, temperature, NULL, dSigma);
				}
				else{
					for(m=0;m<2*NOPACITY;m++)
						dSigma[m] = 0.0;
				}
				

				/* dSigma[0] = dSigma_sF/drho, dSigma[1] = dSigma_aF/drho, dSigma[2]=dSigma_aP/drho, dSigma[3]= dSigma_aE/drho */
				/* dSigma[4] = dSigma_sF/dT, dSigma[5] = dSigma_aF/dT, dSigma[6]=dSigma_aP/dT, dSigma[7]= dSigma_aE/dT */
						

				dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
				dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
				dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
				dSigmadP[3] =  dSigma[7] / (density * R_ideal); 

				diffTErdP = dSigmadP[2] * pow(temperature, 4.0) - dSigmadP[3] * pG->U[k][j][i].Er;

				/* The velocity used to convert primitive variable and conservative variables are not the same as velocity in source term */


				velocity = velocity_x * velocity_x + velocity_z * velocity_z + velocity_y * velocity_y;

				/* If  opacity depends on density, temperature, SPP may becomes positive and unstable */
		/*
				SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * pow(temperature, 3.0) * (1.0 - velocity/(Crat * Crat)) /(density * R_ideal)
					-(Gamma - 1.0) * Prat * Crat * diffTErdP * (1.0 - velocity/(Crat * Crat))
				      -(Gamma - 1.0) * Prat * 2.0 * dSigmadP[1] * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
		*/
				/* If SPP is positive, then this is numerical unstable. */
				/* We need to subtract the unstable mode and gurantee that it is negative */ 
				
				SPP = -4.0 * (Gamma - 1.0) * Prat * Crat * Sigma_aP * pow(temperature, 3.0) * (1.0 - velocity/(Crat * Crat)) /(density * R_ideal);					
				

		/*===================================================================*/
		/* In case velocity is large, momentum source term is also stiff */
				SVVx = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er + diffTEr) / (density * Crat);
		
				if(fabs(SVVx * dt * 0.5) > 0.001)
				betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
				else 
				betax = 1.0 + 0.25 * SVVx * dt;

				SVVy = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er + diffTEr) / (density * Crat);
		
				if(fabs(SVVy * dt * 0.5) > 0.001)
				betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
				else 
				betay = 1.0 + 0.25 * SVVy * dt;

				SVVz = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er + diffTEr) / (density * Crat);
		
				if(fabs(SVVz * dt * 0.5) > 0.001)
				betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
				else 	
				betaz = 1.0 + 0.25 * SVVz * dt;
		/*===========================================================================*/
	
				if(fabs(SPP * dt * 0.5) > 0.001)
				alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
				else 
				alpha = 1.0 + 0.25 * SPP * dt;


				Alpha[k][j][i] = alpha;
				Beta[k][j][i][0] = betax;
				Beta[k][j][i][1] = betay;
				Beta[k][j][i][2] = betaz;

			}
			else{
	
				/* If pressure becomes negative, we do not add radiation source term */
				Source[k][j][i][0] = 0.0;
				Source[k][j][i][1] = 0.0;
				Source[k][j][i][2] = 0.0;
				Source[k][j][i][3] = 0.0;
				Source[k][j][i][4] = 0.0;

				Alpha[k][j][i] = 1.0;
				Beta[k][j][i][0] = 1.0;
				Beta[k][j][i][1] = 1.0;
				Beta[k][j][i][2] = 1.0;

			}

			} /* End i */
		} /* End j */
	} /* End k */

}


#endif /* CTU_INTEGRATOR */
