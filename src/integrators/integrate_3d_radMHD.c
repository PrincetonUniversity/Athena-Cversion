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


#ifdef FIRST_ORDER_FLUX_CORRECTION
static Cons1DS ***x1FluxP=NULL, ***x2FluxP=NULL, ***x3FluxP=NULL;
static Real ***dhalfP = NULL;
static Real ***emf1P=NULL, ***emf2P=NULL, ***emf3P=NULL;
static Cons1DS x1FD_i, x1FD_ip1, x2FD_j, x2FD_jp1, x3FD_k, x3FD_kp1;
static Real ****INIdata;
#if defined(MHD) || defined(RADIATION_MHD)
static Real emf1D_kj, emf1D_kp1j, emf1D_kjp1, emf1D_kp1jp1;
static Real emf2D_ki, emf2D_kip1, emf2D_kp1i, emf2D_kp1ip1;
static Real emf3D_ji, emf3D_jip1, emf3D_jp1i, emf3D_jp1ip1;
#endif
#endif

#ifdef CONS_GRAVITY
static Real ***x1Flux_grav=NULL;
static Real ***x2Flux_grav=NULL;
static Real ***x3Flux_grav=NULL;
static Real ***density_old=NULL;
Real dotphil, dotgxl;
#endif


/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static ConsS **Flxiib=NULL, **Flxoib=NULL;
static ConsS **rFlxiib=NULL, **rFlxoib=NULL;
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

static Real zmin;
static Real Zmax = 1.5;
static Real Vmax = 10.0;
static Real betafloor = 0.0;
static Real dfloor = 1.e-20;
static Real Tfloor = 1.e-20;
/* electron equilivalent rest mass temperature, used in compton scattering */
static Real T_e = 5.94065e9; 
/*static Real T0 = 2.63375e7; *
 * Now defined in global */
 


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

#ifdef FIRST_ORDER_FLUX_CORRECTION
static void integrate_emf1_corner_FOFC(const GridS *pG);
static void integrate_emf2_corner_FOFC(const GridS *pG);
static void integrate_emf3_corner_FOFC(const GridS *pG);
#endif

#endif /* MHD */

#ifdef FIRST_ORDER_FLUX_CORRECTION
static void FixCell(GridS *pG, Int3Vect);
static void ApplyCorr(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3);
static void ApplyCorr2(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3);
static void FixMHD(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3);

void FOFC_Flux(const Cons1DS Ul, const Cons1DS Ur,
                   const Prim1DS Wl, const Prim1DS Wr,
                   const Real Bxi, Cons1DS *pFlux);
#endif


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

#ifdef STATIC_MESH_REFINEMENT
  	int ncg,npg,dim;
  	int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif

	/* For static gravitational potential */
	Real x1, x2, x3, phicl, phicr, phifc, phic, phir, phil;
	/* parameters used in Compton Scattering */
	Real Tr, coefA, coefK, coefB, coef1, coef2, coef3, coef4;
	Real Ersum;
	/* decide the way to add source term */
	int Sourceflag = 1;
	int Sourceflag2 = 1;
	int badcellflag = 0;
	
	/* The safe factor used in the Riemann solver to decide whether adopt first order 
	* left and right states or not */
	const Real fluxfactor = 2.0;

	/* To correct the negative pressure */
	zmin = 4.0;
	Real beta, Bpre;
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

#ifdef FIRST_ORDER_FLUX_CORRECTION
  int flag_cell=0,negd=0,negP=0,superl=0,NaNFlux=0;
  Real Vsq;
  Int3Vect BadCell;
  PrimS Wtemp1, Wtemp2, Wtemp3, Wtemp4, Wtemp5, Wtemp6;
  int Np; /* number of cells used to average pressure */

#endif 
	
  PrimS Wtemp;
  PrimS Wopacity;
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
	Real Smx0, Smy0, Smz0; /* the momentum source term v(T^4-E_r)* sigma/Crat */
	Real Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE;
	/* The opacity function is:0-3 Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE */
	Real Propa_44, SEE, SErho, SEmx, SEmy, SEmz;
	Real dSigma[2*NOPACITY];
	Real Sigma[NOPACITY];
	/* dSigma for dSigma?/drho and dSigma?/dT */
	Cons1DS Usource;
	/* for source term */

	/* In case momentum becomes stiff */
	Real SFmx, SFmy, SFmz;
	/* for radiation work term */
	Real Prwork1, Prworksource, Prwork2;


	Real Source_Inv[NVAR][NVAR], Source_Invnew[NVAR][NVAR], tempguess[NVAR], Uguess[NVAR], Source_guess[NVAR], Errort[NVAR], SourceFlux[NVAR];
	Real Psource, Ersource;
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
			Source_Invnew[i][j] = 0.0;
			
			if(i==j) {
		 		Source_Inv[i][j] = 1.0;
				Source_Invnew[i][j] = 1.0;
		
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


	/* Check negative pressure */
	/* This is particular for the ghost zones */
	
 	for (k=ks-nghost; k<=ke+nghost; k++){
                for (j=js-nghost; j<=je+nghost; j++) {
                        for (i=is-nghost; i<=ie+nghost; i++) {

                                Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
			if(Wtemp.d < dfloor){
				Wtemp.d = dfloor;
				pG->U[k][j][i] = Prim_to_Cons(&(Wtemp));
			}

			if(Wtemp.P/(R_ideal * Wtemp.d) < Tfloor){
					Wtemp.P = Tfloor * Wtemp.d * R_ideal;
                                        pG->U[k][j][i] = Prim_to_Cons(&(Wtemp));


                                } /* End if wtemp negative */
                        }/* end i */
                }/* end j*/
        }/* end k*/




/*=== Step 1: Backward Euler is done in the main function for the whole mesh ===*/
	
	/* First, calculate the source term */
	updatesource(pG);

/* Store the initial data */
#ifdef FIRST_ORDER_FLUX_CORRECTION
	for (k=ks-nghost; k<=ke+nghost; k++){
		for (j=js-nghost; j<=je+nghost; j++) {
    			for (i=is-nghost; i<=ie+nghost; i++) {
				INIdata[k][j][i][0] = pG->U[k][j][i].d;
				INIdata[k][j][i][1] = pG->U[k][j][i].M1;
				INIdata[k][j][i][2] = pG->U[k][j][i].M2;
				INIdata[k][j][i][3] = pG->U[k][j][i].M3;
				INIdata[k][j][i][4] = pG->U[k][j][i].E;

				Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));

				INIdata[k][j][i][5] = Wtemp.P;


			}
		}
	}	

#endif


	
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

/* Do not include source terms for the left and right states */
/* Using first order reconstruction and do not include source terms to the left and right states */
#ifdef FIRST_ORDER_FLUX_CORRECTION
			
			for (i=il+1; i<=iu; i++) {
        			Ul_x1Face[k][j][i] = U1d[i-1];
        			Ur_x1Face[k][j][i] = U1d[i];
				
				Wl[i] = Cons1D_to_Prim1D(&U1d[i-1],&Bxc[i-1]);
				Wr[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);

#ifdef RADIATION_MHD
        			Bx = B1_x1Face[k][j][i];
#endif
			/* take the time step with x1Flux */
				x1FluxP[k][j][i].d = dt;	
				x1FluxP[k][j][i].Mx = 3;

        		FOFC_Flux(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1FluxP[k][j][i]);
      		}


#endif /* End first order flux */

			for (i=is-nghost; i<=ie+nghost; i++) {
      				W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
    			}

			lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,3);
			
/*			for (i=il+1; i<=iu; i++) {
				velocity = sqrt(Wl[i].Vx * Wl[i].Vx + Wl[i].Vy * Wl[i].Vy + Wl[i].Vz * Wl[i].Vz);
				velocity_x = sqrt(Wr[i].Vx * Wr[i].Vx + Wr[i].Vy * Wr[i].Vy + Wr[i].Vz * Wr[i].Vz);
				if( (velocity > Vmax) || (velocity_x > Vmax) || (velocity!= velocity) || (velocity_x != velocity_x) || (Wl[i].P/(Wl[i].d * R_ideal) < Tfloor) || (Wr[i].P/(Wr[i].d * R_ideal) < Tfloor)){
					Wl[i] = W[i-1];
					Wr[i] = W[i];
				}
			}
*/

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
			
			if(Wl[i].P < TINY_NUMBER){ 
				Wl[i].P -= dt * Propa_44 * Psource * 0.5;
				Wl[i].Vx -= dt * Source[k][j][i-1][1] * 0.5 * Beta[k][j][i-1][0] / U1d[i-1].d;
				Wl[i].Vy -= dt * Source[k][j][i-1][2] * 0.5 * Beta[k][j][i-1][1] / U1d[i-1].d;
				Wl[i].Vz -= dt * Source[k][j][i-1][3] * 0.5 * Beta[k][j][i-1][2] / U1d[i-1].d;
			}	
	
			for(m=0; m<NOPACITY; m++)
				Wl[i].Sigma[m] = U1d[i-1].Sigma[m];
			
			Wl[i].Edd_11 = W[i-1].Edd_11;
			Wl[i].Edd_21 = W[i-1].Edd_21;
			Wl[i].Edd_22 = W[i-1].Edd_22;
			Wl[i].Edd_31 = W[i-1].Edd_31;
			Wl[i].Edd_32 = W[i-1].Edd_32;
			Wl[i].Edd_33 = W[i-1].Edd_33;

		/* For the right state */
	
	
			Propa_44 = Alpha[k][j][i];

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
			
			if(Wr[i].P < TINY_NUMBER){
				Wr[i].P -= dt * Propa_44 * Psource * 0.5;
				Wr[i].Vx -= dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[i].d;
				Wr[i].Vy -= dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[i].d;
				Wr[i].Vz -= dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[i].d;
			}

			for(m=0; m<NOPACITY; m++)
				Wr[i].Sigma[m] = U1d[i].Sigma[m];

			Wr[i].Edd_11 = W[i].Edd_11;
			Wr[i].Edd_21 = W[i].Edd_21;
			Wr[i].Edd_22 = W[i].Edd_22;
			Wr[i].Edd_31 = W[i].Edd_31;
			Wr[i].Edd_32 = W[i].Edd_32;
			Wr[i].Edd_33 = W[i].Edd_33;	

/*			velocity = sqrt(Wl[i].Vx * Wl[i].Vx + Wl[i].Vy * Wl[i].Vy + Wl[i].Vz * Wl[i].Vz);
                        velocity_x = sqrt(Wr[i].Vx * Wr[i].Vx + Wr[i].Vy * Wr[i].Vy + Wr[i].Vz * Wr[i].Vz);	
			if(velocity != velocity || velocity_x != velocity_x){
				printf("lvx: %e lvy: %e lvz: %e rvx: %e rvy: %e rvz: %e ls: %e lBeta: %e rs: %e rBeta: %e\n",Wl[i].Vx,Wl[i].Vy,Wl[i].Vz,Wr[i].Vx,Wr[i].Vz,Wr[i].Vz,Source[k][j][i-1][1],Beta[k][j][i-1][0],Source[k][j][i][1],Beta[k][j][i][0]);

			}
*/			
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


				if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor) || (Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor)){
                                        /* copy left primitive states to be the left sate */
                                        Wl[i] = Cons1D_to_Prim1D(&U1d[i-1],&Bxc[i-1]);
                                        Ul_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bxi[i]);

                                        /* copy right primitive states to be the right sate */
                                        Wr[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i]);
                                        Ur_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bxi[i]);

                                }


#ifdef RADIATION_MHD
        		Bx = B1_x1Face[k][j][i];
#endif
			/* take the time step with x1Flux */
			x1Flux[k][j][i].d = dt;	
				x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);

			 /* revert to predictor flux if this flux Nan'ed */
/*
                                if ((x1Flux[k][j][i].d  != x1Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
                                        (x1Flux[k][j][i].E  != x1Flux[k][j][i].E)  ||
#endif
#ifdef RADIATION_MHD
                                        (x1Flux[k][j][i].By != x1Flux[k][j][i].By) ||
                                        (x1Flux[k][j][i].Bz != x1Flux[k][j][i].Bz) ||
#endif
                                        (x1Flux[k][j][i].Mx != x1Flux[k][j][i].Mx) ||
                                        (x1Flux[k][j][i].My != x1Flux[k][j][i].My) ||
                                        (x1Flux[k][j][i].Mz != x1Flux[k][j][i].Mz)) {
                                        printf("1nd flux1: ld: %e lp: %e lvx: %e lvy: %e lvz: %e lBx: %e lby: %e lbz: %e Uld: %e ulE: %e UlMx: %e Ulmy: %e Ulmz: %e rd: %e rp: %e rvx: %e rvy: %e rvz: %e rby: %e rbz: %e Urd: %e urE: %e UrMx: %e Urmy: %e Urmz: %e \n",Wl[i].d,Wl[i].P, Wl[i].Vx,Wl[i].Vy,Wl[i].Vz,Bx,Wl[i].By,Wl[i].Bz,pG->U[k][j][i-1].d,pG->U[k][j][i-1].E,pG->U[k][j][i-1].M1,pG->U[k][j][i-1].M2,pG->U[k][j][i-1].M3,Wr[i].d,Wr[i].P,Wr[i].Vx,Wr[i].Vy,Wr[i].Vz,Wr[i].By,Wr[i].Bz,pG->U[k][j][i].d,pG->U[k][j][i].E,pG->U[k][j][i].M1,pG->U[k][j][i].M2,pG->U[k][j][i].M3);
                                         printf("i: %d j: %d k: %d ID: %d\n",i,j,k,myID_Comm_world);
                                }
*/

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

/* Do not include source terms for the left and right states */
#ifdef FIRST_ORDER_FLUX_CORRECTION
			
			for (j=jl+1; j<=ju; j++) {
        			Ul_x2Face[k][j][i] = U1d[j-1];
        			Ur_x2Face[k][j][i] = U1d[j];

				Wl[j] = Cons1D_to_Prim1D(&U1d[j-1],&Bxc[j-1]);
				Wr[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);


#ifdef RADIATION_MHD
        			Bx = B2_x2Face[k][j][i];
#endif
				/* Take the time step with x2Flux[][][].d */
				x2FluxP[k][j][i].d = dt;
				x2FluxP[k][j][i].Mx = 3;
			
        			FOFC_Flux(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[j],Wr[j],Bx,&x2FluxP[k][j][i]);
      			}


#endif /* End first order flux */


		for (j=js-nghost; j<=je+nghost; j++) {
        		W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx2,jl+1,ju-1,Wl,Wr,3);

/*		for (j=jl+1; j<=ju; j++) {
				velocity = sqrt(Wl[j].Vx * Wl[j].Vx + Wl[j].Vy * Wl[j].Vy + Wl[j].Vz * Wl[j].Vz);
				velocity_x = sqrt(Wr[j].Vx * Wr[j].Vx + Wr[j].Vy * Wr[j].Vy + Wr[j].Vz * Wr[j].Vz);
				if((velocity > Vmax) || (velocity_x > Vmax) || (velocity_x != velocity_x) || (Wl[j].P/(Wl[j].d * R_ideal) < Tfloor) || (Wr[j].P/(Wr[j].d * R_ideal) < Tfloor)){
					Wl[j] = W[j-1];
					Wr[j] = W[j];
				}
		}
*/

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
			
			if(Wl[j].P < TINY_NUMBER){
				Wl[j].P -= dt * Propa_44 * Psource * 0.5;
				Wl[j].Vx -= dt * Source[k][j-1][i][2] * 0.5 * Beta[k][j-1][i][1] / U1d[j-1].d;
				Wl[j].Vy -= dt * Source[k][j-1][i][3] * 0.5 * Beta[k][j-1][i][2] / U1d[j-1].d;
				Wl[j].Vz -= dt * Source[k][j-1][i][1] * 0.5 * Beta[k][j-1][i][0] / U1d[j-1].d;
			}

			
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

			if(Wr[j].P < TINY_NUMBER){
				Wr[j].P -= dt * Propa_44 * Psource * 0.5;
				Wr[j].Vx -= dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[j].d;
				Wr[j].Vy -= dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[j].d;
				Wr[j].Vz -= dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[j].d;		
			}
			
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
				
 				if((Wl[j].P < 2.0 * TINY_NUMBER) || (Wl[j].d < dfloor) || (Wr[j].P < 2.0 * TINY_NUMBER) || (Wr[j].d < dfloor) ){
                                        /* copy left primitive states to be the left sate */
                                        Wl[j] = Cons1D_to_Prim1D(&U1d[j-1],&Bxc[j-1]);
                                        Ul_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wl[j],&Bxi[j]);

                                        /* copy right primitive states to be the right sate */
                                        Wr[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j]);
                                        Ur_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wr[j],&Bxi[j]);

                                }


#ifdef RADIATION_MHD
        		Bx = B2_x2Face[k][j][i];
#endif
			/* Take the time step with x2Flux[][][].d */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[j],Wr[j],Bx,&x2Flux[k][j][i]);

/*
                        if ((x2Flux[k][j][i].d  != x2Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
                                (x2Flux[k][j][i].E  != x2Flux[k][j][i].E)  ||
#endif
#if defined(MHD) || defined(RADIATION_MHD)
                                (x2Flux[k][j][i].By != x2Flux[k][j][i].By) ||
                                (x2Flux[k][j][i].Bz != x2Flux[k][j][i].Bz) ||
#endif
                                (x2Flux[k][j][i].Mx != x2Flux[k][j][i].Mx) ||
                                (x2Flux[k][j][i].My != x2Flux[k][j][i].My) ||
                                (x2Flux[k][j][i].Mz != x2Flux[k][j][i].Mz)) {
                                printf("1nd flux2: ld: %e lp: %e lvx: %e lvy: %e lvz: %e lBx: %e lby: %e lbz: %e Uld: %e ulE: %e UlMx: %e Ulmy: %e Ulmz: %e rd: %e rp: %e rvx: %e rvy: %e rvz: %e rby: %e rbz: %e Urd: %e urE: %e UrMx: %e Urmy: %e Urmz: %e \n",Wl[j].d,Wl[j].P, Wl[j].Vx,Wl[j].Vy,Wl[j].Vz,Bx,Wl[j].By,Wl[j].Bz,pG->U[k][j-1][i].d,pG->U[k][j-1][i].E,pG->U[k][j-1][i].M1,pG->U[k][j-1][i].M2,pG->U[k][j-1][i].M3,Wr[j].d,Wr[j].P,Wr[j].Vx,Wr[j].Vy,Wr[j].Vz,Wr[j].By,Wr[j].Bz,pG->U[k][j][i].d,pG->U[k][j][i].E,pG->U[k][j][i].M1,pG->U[k][j][i].M2,pG->U[k][j][i].M3);
                        }
*/

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


/* Do not include source terms for the left and right states */
#ifdef FIRST_ORDER_FLUX_CORRECTION
			
			for (k=kl+1; k<=ku; k++) {
        			Ul_x3Face[k][j][i] = U1d[k-1];
        			Ur_x3Face[k][j][i] = U1d[k];

				Wl[k] = Cons1D_to_Prim1D(&U1d[k-1],&Bxc[k-1]);
				Wr[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);


#ifdef RADIATION_MHD
        			Bx = B3_x3Face[k][j][i];
#endif
			/* Take the time step with x3Flux[][][].d */
				x3FluxP[k][j][i].d = dt;
				x3FluxP[k][j][i].Mx = 3;
			
        			FOFC_Flux(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[k],Wr[k],Bx,&x3FluxP[k][j][i]);
      		}

#endif /* End first order flux */




		for (k=ks-nghost; k<=ke+nghost; k++) {
        		W[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
      		}

      		lr_states(pG,W,Bxc,pG->dt,dx3,kl+1,ku-1,Wl,Wr,3);

/*		for (k=kl+1; k<=ku; k++) {
				velocity = sqrt(Wl[k].Vx * Wl[k].Vx + Wl[k].Vy * Wl[k].Vy + Wl[k].Vz * Wl[k].Vz);
				velocity_x = sqrt(Wr[k].Vx * Wr[k].Vx + Wr[k].Vy * Wr[k].Vy + Wr[k].Vz * Wr[k].Vz);
				if((velocity > Vmax) || (velocity_x > Vmax) || (velocity_x != velocity_x) || (Wl[k].P/(Wl[k].d * R_ideal) < Tfloor) || (Wr[k].P/(Wr[k].d * R_ideal) < Tfloor)){
					Wl[k] = W[k-1];
					Wr[k] = W[k];
				}
		}
*/
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
			
			if(Wl[k].P < TINY_NUMBER){
				Wl[k].P -= dt * Propa_44 * Psource * 0.5;
				Wl[k].Vx -= dt * Source[k-1][j][i][3] * 0.5 * Beta[k-1][j][i][2] / U1d[k-1].d;
				Wl[k].Vy -= dt * Source[k-1][j][i][1] * 0.5 * Beta[k-1][j][i][0] / U1d[k-1].d;
				Wl[k].Vz -= dt * Source[k-1][j][i][2] * 0.5 * Beta[k-1][j][i][1] / U1d[k-1].d;
			}
			
	
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
			
			if(Wr[k].P < TINY_NUMBER){
				Wr[k].P -= dt * Propa_44 * Psource * 0.5;
				Wr[k].Vx -= dt * Source[k][j][i][3] * 0.5 * Beta[k][j][i][2] / U1d[k].d;
				Wr[k].Vy -= dt * Source[k][j][i][1] * 0.5 * Beta[k][j][i][0] / U1d[k].d;
				Wr[k].Vz -= dt * Source[k][j][i][2] * 0.5 * Beta[k][j][i][1] / U1d[k].d;
			}

			for(m=0; m<NOPACITY;m++){
				Wr[k].Sigma[m] = U1d[k].Sigma[m];
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
        		for (k=kl+1; k<=ku; k++) {
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

				 if((Wl[k].P < 2.0 * TINY_NUMBER) || (Wl[k].d < dfloor) || (Wr[k].P < 2.0 * TINY_NUMBER) || (Wr[k].d < dfloor)){
                                        /* copy left primitive states to be the left sate */
                                        Wl[k] = Cons1D_to_Prim1D(&U1d[k-1],&Bxc[k-1]);
                                        Ul_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wl[k],&Bxi[k]);

                                        /* copy right primitive states to be the right sate */
                                        Wr[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k]);
                                        Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[k],&Bxi[k]);

                                }
			

#ifdef RADIATION_MHD
        		Bx = B3_x3Face[k][j][i];
#endif
			/* Take the time step with x3Flux[][][].d */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[k],Wr[k],Bx,&x3Flux[k][j][i]);

/*
                        if ((x3Flux[k][j][i].d  != x3Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
                                (x3Flux[k][j][i].E  != x3Flux[k][j][i].E)  ||
#endif
#if defined(MHD) || defined(RADIATION_MHD)
                                (x3Flux[k][j][i].By != x3Flux[k][j][i].By) ||
                                (x3Flux[k][j][i].Bz != x3Flux[k][j][i].Bz) ||
#endif
                                (x3Flux[k][j][i].Mx != x3Flux[k][j][i].Mx) ||
                                (x3Flux[k][j][i].My != x3Flux[k][j][i].My) ||
                                (x3Flux[k][j][i].Mz != x3Flux[k][j][i].Mz)) {
                                printf("1nd flux3: ld: %e lp: %e lvx: %e lvy: %e lvz: %e lBx: %e lby: %e lbz: %e Uld: %e ulE: %e UlMx: %e Ulmy: %e Ulmz: %e rd: %e rp: %e rvx: %e rvy: %e rvz: %e rby: %e rbz: %e Urd: %e urE: %e UrMx: %e Urmy: %e Urmz: %e \n",Wl[k].d,Wl[k].P,Wl[k]
.Vx,Wl[k].Vy,Wl[k].Vz,Bx,Wl[k].By,Wl[k].Bz,pG->U[k-1][j][i].d,pG->U[k-1][j][i].E,pG->U[k-1][j][i].M1,pG->U[k-1][j][i].M2,pG->U[k-1][j][i].M3,Wr[k].d,Wr[k].P,Wr[k].Vx,Wr[k].Vy,Wr[k].Vz,Wr[k].By,Wr[k].Bz,pG->U[k][j][i].d,pG->U[k][j][i].E,pG->U[k][j][i].M1,pG->U[k][j][i].M2,pG->U[
k][j][i].M3);
                        }
*/

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




#ifdef FIRST_ORDER_FLUX_CORRECTION
/* Save first order flux */
#if defined(MHD) || defined(RADIATION_MHD)
	/* calculate the EMF with the first order flux */
	/* emf1,2, 3 P are updated there */
  	integrate_emf1_corner_FOFC(pG);
  	integrate_emf2_corner_FOFC(pG);
  	integrate_emf3_corner_FOFC(pG);
#endif
#if defined (RADIATION_MHD) && defined (SHEARING_BOX)
  get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);

/* compute remapped Ey from opposite side of grid */
  if (my_iproc == 0) {
      RemapEyFlux_ix1(pD, emf2P, remapEyiib, x1FluxP, remapx1Fiib);
  }
  if (my_iproc == (pD->NGrid[0]-1)) {
      RemapEyFlux_ox1(pD, emf2P, remapEyoib, x1FluxP, remapx1Foib);
  }


/* Now average Ey and remapped Ey */

	
    if (my_iproc == 0) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2P[k][j][is]  = 0.5*(emf2P[k][j][is] + remapEyiib[k][j]);
				x1FluxP[k][j][is].d = 0.5*(x1FluxP[k][j][is].d + remapx1Fiib[k][j]);
			}
		}
    }
	
    if (my_iproc == (pD->NGrid[0]-1)) {
		for(k=ks; k<=ke+1; k++) {
			for(j=js; j<=je; j++){
				emf2P[k][j][ie+1]  = 0.5*(emf2P[k][j][ie+1] + remapEyoib[k][j]);
				x1FluxP[k][j][ie+1].d = 0.5*(x1FluxP[k][j][ie+1].d + remapx1Foib[k][j]);
			}
		}
    }


#endif /* MHD & SHEARING_BOX */


#endif /* End FOFC */


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


				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_y * Fr0y + velocity_z * Fr0z)/Crat);

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

				
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_y * Fr0y + velocity_z * Fr0z)/Crat);

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

				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_z * Fr0z)/Crat);

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

				
				Source[4] = -Prat * Crat * ((Sigma_a - Sigma_s) * (velocity_x * Fr0x + velocity_y * Fr0y)/Crat);

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
  


/*--- Step 9a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (with FOFC flux)
 */
#ifdef FIRST_ORDER_FLUX_CORRECTION
    	for (k=kl+1; k<=ku-1; k++) {
      		for (j=jl+1; j<=ju-1; j++) {
			for (i=il+1; i<=iu-1; i++) {

          			dhalfP[k][j][i] = pG->U[k][j][i].d 
            				- hdtodx1 * (    x1FluxP[k  ][j  ][i+1].d -     x1FluxP[k][j][i].d)
            				- hdtodx2 * (    x2FluxP[k  ][j+1][i  ].d -     x2FluxP[k][j][i].d)
            				- hdtodx3 * (    x3FluxP[k+1][j  ][i  ].d -     x3FluxP[k][j][i].d);

			}
      		}
    	}
#endif /* FOFC */

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

			/* Part of momentum source term */
			/* For this term, use Tguess, not temperature ^4 */
			Smx0 = -Prat * velocity_x * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;
			Smy0 = -Prat * velocity_y * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;
			Smz0 = -Prat * velocity_z * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;

			/* The Source term */
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);

		/*=========================================================*/
		/* In case velocity is large and momentum source is stiff */
			SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat); 
			/*	+ diffTEr / (density * Crat);
			*/	

			SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat); 
			/*	+ diffTEr / (density * Crat);
			*/

			SFmz = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat); 
			/*	+ diffTEr / (density * Crat);
			*/	


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
					
			/* Part of radiation momentum source term */
			M1h += hdt * Smx0;
			M2h += hdt * Smy0;
			M3h += hdt * Smz0;

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

			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[k][j][i],&Bx);
			Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[k][j][i],&Bx);
/*

			velocity = sqrt(Wl[i].Vx * Wl[i].Vx + Wl[i].Vy * Wl[i].Vy + Wl[i].Vz * Wl[i].Vz);
                        velocity_x = sqrt(Wr[i].Vx * Wr[i].Vx + Wr[i].Vy * Wr[i].Vy + Wr[i].Vz * Wr[i].Vz);

			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor) ||(Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor) || (velocity> Vmax) || (velocity_x > Vmax) || (fabs(x3) > Zmax)){
*/


			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor) ||(Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor)){
#ifdef RADIATION_MHD
				B1_x1Face[k][j][i] = pG->B1i[k][j][i];
				Bx = B1_x1Face[k][j][i];
#endif
				 Wl[i].P = pG->U[k][j][i-1].E - 0.5 * (pG->U[k][j][i-1].M1 * pG->U[k][j][i-1].M1 + pG->U[k][j][i-1].M2 * pG->U[k][j][i-1].M2 + pG->U[k][j][i-1].M3 * pG->U[k][j][i-1].M3) / pG->U[k][j][i-1].d;
#ifdef RADIATION_MHD
                                Wl[i].P -= 0.5 * (pG->U[k][j][i-1].B1c * pG->U[k][j][i-1].B1c + pG->U[k][j][i-1].B2c * pG->U[k][j][i-1].B2c + pG->U[k][j][i-1].B3c * pG->U[k][j][i-1].B3c);
#endif

                                Wl[i].P *= (Gamma - 1.0);

                                Wl[i].d = pG->U[k][j][i-1].d;
                                Wl[i].Vx = pG->U[k][j][i-1].M1 / pG->U[k][j][i-1].d;
                                Wl[i].Vy = pG->U[k][j][i-1].M2 / pG->U[k][j][i-1].d;
                                Wl[i].Vz = pG->U[k][j][i-1].M3 / pG->U[k][j][i-1].d;

                                Ul_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bx);


                                Wr[i].P = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
                                Wr[i].P -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif

                                Wr[i].P *= (Gamma - 1.0);

                                Wr[i].d = pG->U[k][j][i].d;
                                Wr[i].Vx = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
                                Wr[i].Vy = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                                Wr[i].Vz = pG->U[k][j][i].M3 / pG->U[k][j][i].d;

                                Ur_x1Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bx);
#ifdef RADIATION_MHD

				emf1_cc[k][j][i-1] = (pG->U[k][j][i-1].B2c*pG->U[k][j][i-1].M3 -
			    		pG->U[k][j][i-1].B3c*pG->U[k][j][i-1].M2)/pG->U[k][j][i-1].d;
        			emf2_cc[k][j][i-1] = (pG->U[k][j][i-1].B3c*pG->U[k][j][i-1].M1 -
			    		pG->U[k][j][i-1].B1c*pG->U[k][j][i-1].M3)/pG->U[k][j][i-1].d;
        			emf3_cc[k][j][i-1] = (pG->U[k][j][i-1].B1c*pG->U[k][j][i-1].M2 -
			    		pG->U[k][j][i-1].B2c*pG->U[k][j][i-1].M1)/pG->U[k][j][i-1].d;

				emf1_cc[k][j][i] = (pG->U[k][j][i].B2c*pG->U[k][j][i].M3 -
			    		pG->U[k][j][i].B3c*pG->U[k][j][i].M2)/pG->U[k][j][i].d;
        			emf2_cc[k][j][i] = (pG->U[k][j][i].B3c*pG->U[k][j][i].M1 -
			    		pG->U[k][j][i].B1c*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
        			emf3_cc[k][j][i] = (pG->U[k][j][i].B1c*pG->U[k][j][i].M2 -
			    		pG->U[k][j][i].B2c*pG->U[k][j][i].M1)/pG->U[k][j][i].d;
#endif

			}
			
		
			/* Need parameter dt in radiation Riemann solver */
			x1Flux[k][j][i].d = dt;
			x1Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx, &x1Flux[k][j][i]);



#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux Nan'ed */
        if ((x1Flux[k][j][i].d  != x1Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x1Flux[k][j][i].E  != x1Flux[k][j][i].E)  ||
#endif
#ifdef RADIATION_MHD
            (x1Flux[k][j][i].By != x1Flux[k][j][i].By) ||
            (x1Flux[k][j][i].Bz != x1Flux[k][j][i].Bz) ||
#endif
            (x1Flux[k][j][i].Mx != x1Flux[k][j][i].Mx) ||
            (x1Flux[k][j][i].My != x1Flux[k][j][i].My) ||
            (x1Flux[k][j][i].Mz != x1Flux[k][j][i].Mz)) {
          x1Flux[k][j][i] = x1FluxP[k][j][i];
          NaNFlux++;
        }
#endif /* End FOFC */
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

			cc_pos(pG,i,j,k,&x1,&x2,&x3);
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[k][j][i],&Bx);


/*

			velocity = sqrt(Wl[i].Vx * Wl[i].Vx + Wl[i].Vy * Wl[i].Vy + Wl[i].Vz * Wl[i].Vz);
                        velocity_x = sqrt(Wr[i].Vx * Wr[i].Vx + Wr[i].Vy * Wr[i].Vy + Wr[i].Vz * Wr[i].Vz);


			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor)  || (Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor) || (velocity > Vmax) || (velocity_x > Vmax) || (fabs(x3) > Zmax)){


*/

			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor)  || (Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor)){
#ifdef RADIATION_MHD			
				B2_x2Face[k][j][i] = pG->B2i[k][j][i];
				Bx = B2_x2Face[k][j][i];
#endif	
 				Wl[i].P = pG->U[k][j-1][i].E - 0.5 * (pG->U[k][j-1][i].M1 * pG->U[k][j-1][i].M1 + pG->U[k][j-1][i].M2 * pG->U[k][j-1][i].M2 + pG->U[k][j-1][i].M3 * pG->U[k][j-1][i].M3) / pG->U[k][j-1][i].d;
#ifdef RADIATION_MHD
                                Wl[i].P -= 0.5 * (pG->U[k][j-1][i].B1c * pG->U[k][j-1][i].B1c + pG->U[k][j-1][i].B2c * pG->U[k][j-1][i].B2c + pG->U[k][j-1][i].B3c * pG->U[k][j-1][i].B3c);
#endif

                                Wl[i].P *= (Gamma - 1.0);

                                Wl[i].d = pG->U[k][j-1][i].d;
                                Wl[i].Vx = pG->U[k][j-1][i].M2 / pG->U[k][j-1][i].d;
                                Wl[i].Vy = pG->U[k][j-1][i].M3 / pG->U[k][j-1][i].d;
                                Wl[i].Vz = pG->U[k][j-1][i].M1 / pG->U[k][j-1][i].d;

                                Ul_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bx);

                                Wr[i].P = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
                                Wr[i].P -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif

                               Wr[i].P *= (Gamma - 1.0);

                               Wr[i].d = pG->U[k][j][i].d;
                               Wr[i].Vx = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                               Wr[i].Vy = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
                               Wr[i].Vz = pG->U[k][j][i].M1 / pG->U[k][j][i].d;

                               Ur_x2Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bx);
#ifdef RADIATION_MHD
				emf1_cc[k][j-1][i] = (pG->U[k][j-1][i].B2c*pG->U[k][j-1][i].M3 -
			    		pG->U[k][j-1][i].B3c*pG->U[k][j-1][i].M2)/pG->U[k][j-1][i].d;
        			emf2_cc[k][j-1][i] = (pG->U[k][j-1][i].B3c*pG->U[k][j-1][i].M1 -
			    		pG->U[k][j-1][i].B1c*pG->U[k][j-1][i].M3)/pG->U[k][j-1][i].d;
        			emf3_cc[k][j-1][i] = (pG->U[k][j-1][i].B1c*pG->U[k][j-1][i].M2 -
			    		pG->U[k][j-1][i].B2c*pG->U[k][j-1][i].M1)/pG->U[k][j-1][i].d;

				emf1_cc[k][j][i] = (pG->U[k][j][i].B2c*pG->U[k][j][i].M3 -
			    		pG->U[k][j][i].B3c*pG->U[k][j][i].M2)/pG->U[k][j][i].d;
        			emf2_cc[k][j][i] = (pG->U[k][j][i].B3c*pG->U[k][j][i].M1 -
			    		pG->U[k][j][i].B1c*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
        			emf3_cc[k][j][i] = (pG->U[k][j][i].B1c*pG->U[k][j][i].M2 -
			    		pG->U[k][j][i].B2c*pG->U[k][j][i].M1)/pG->U[k][j][i].d;
#endif

                       }


			/* Need parameter dt in radiation Riemann solver */
			x2Flux[k][j][i].d = dt;
			x2Flux[k][j][i].Mx = 3;

        		fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[i],Wr[i],Bx, &x2Flux[k][j][i]);


#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x2Flux[k][j][i].d  != x2Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x2Flux[k][j][i].E  != x2Flux[k][j][i].E)  ||
#endif
#if defined(MHD) || defined(RADIATION_MHD)
            (x2Flux[k][j][i].By != x2Flux[k][j][i].By) ||
            (x2Flux[k][j][i].Bz != x2Flux[k][j][i].Bz) ||
#endif
            (x2Flux[k][j][i].Mx != x2Flux[k][j][i].Mx) ||
            (x2Flux[k][j][i].My != x2Flux[k][j][i].My) ||
            (x2Flux[k][j][i].Mz != x2Flux[k][j][i].Mz)) {
          x2Flux[k][j][i] = x2FluxP[k][j][i];
          NaNFlux++;
        }

#endif /* End FOFC */
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
					
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			
        		Wl[i] = Cons1D_to_Prim1D(&Ul_x3Face[k][j][i],&Bx);
        		Wr[i] = Cons1D_to_Prim1D(&Ur_x3Face[k][j][i],&Bx);

/*

			velocity = sqrt(Wl[i].Vx * Wl[i].Vx + Wl[i].Vy * Wl[i].Vy + Wl[i].Vz * Wl[i].Vz);
                        velocity_x = sqrt(Wr[i].Vx * Wr[i].Vx + Wr[i].Vy * Wr[i].Vy + Wr[i].Vz * Wr[i].Vz);

			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor) || (Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor) || (velocity > Vmax) || (velocity_x > Vmax) || (fabs(x3) > Zmax)){

*/

			if((Wl[i].P < 2.0 * TINY_NUMBER) || (Wl[i].d < dfloor) || (Wr[i].P < 2.0 * TINY_NUMBER) || (Wr[i].d < dfloor)){
#ifdef RADIATION_MHD		
				B3_x3Face[k][j][i] = pG->B3i[k][j][i];
				Bx = B3_x3Face[k][j][i];
#endif
				Wl[i].P = pG->U[k-1][j][i].E - 0.5 * (pG->U[k-1][j][i].M1 * pG->U[k-1][j][i].M1 + pG->U[k-1][j][i].M2 * pG->U[k-1][j][i].M2 + pG->U[k-1][j][i].M3 * pG->U[k-1][j][i].M3) / pG->U[k-1][j][i].d;
#ifdef RADIATION_MHD
                                Wl[i].P -= 0.5 * (pG->U[k-1][j][i].B1c * pG->U[k-1][j][i].B1c + pG->U[k-1][j][i].B2c * pG->U[k-1][j][i].B2c + pG->U[k-1][j][i].B3c * pG->U[k-1][j][i].B3c);
#endif

                                Wl[i].P *= (Gamma - 1.0);

                                Wl[i].d  = pG->U[k-1][j][i].d;
                                Wl[i].Vx = pG->U[k-1][j][i].M3 / pG->U[k-1][j][i].d;
                                Wl[i].Vy = pG->U[k-1][j][i].M1 / pG->U[k-1][j][i].d;
                                Wl[i].Vz = pG->U[k-1][j][i].M2 / pG->U[k-1][j][i].d;

                                Ul_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wl[i],&Bx);


                                Wr[i].P = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
                                Wr[i].P -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif

                                Wr[i].P *= (Gamma - 1.0);

                                Wr[i].d = pG->U[k][j][i].d;
                                Wr[i].Vx = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
                                Wr[i].Vy = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
                                Wr[i].Vz = pG->U[k][j][i].M2 / pG->U[k][j][i].d;


                                 Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[i],&Bx);


#ifdef RADIATION_MHD
				emf1_cc[k-1][j][i] = (pG->U[k-1][j][i].B2c*pG->U[k-1][j][i].M3 -
			    		pG->U[k-1][j][i].B3c*pG->U[k-1][j][i].M2)/pG->U[k-1][j][i].d;
        			emf2_cc[k-1][j][i] = (pG->U[k-1][j][i].B3c*pG->U[k-1][j][i].M1 -
			    		pG->U[k-1][j][i].B1c*pG->U[k-1][j][i].M3)/pG->U[k-1][j][i].d;
        			emf3_cc[k-1][j][i] = (pG->U[k-1][j][i].B1c*pG->U[k-1][j][i].M2 -
			    		pG->U[k-1][j][i].B2c*pG->U[k-1][j][i].M1)/pG->U[k-1][j][i].d;

				emf1_cc[k][j][i] = (pG->U[k][j][i].B2c*pG->U[k][j][i].M3 -
			    		pG->U[k][j][i].B3c*pG->U[k][j][i].M2)/pG->U[k][j][i].d;
        			emf2_cc[k][j][i] = (pG->U[k][j][i].B3c*pG->U[k][j][i].M1 -
			    		pG->U[k][j][i].B1c*pG->U[k][j][i].M3)/pG->U[k][j][i].d;
        			emf3_cc[k][j][i] = (pG->U[k][j][i].B1c*pG->U[k][j][i].M2 -
			    		pG->U[k][j][i].B2c*pG->U[k][j][i].M1)/pG->U[k][j][i].d;
#endif

			}
	
	
			/* Need parameter dt in radiation Riemann solver */
			x3Flux[k][j][i].d = dt;
			x3Flux[k][j][i].Mx = 3;
			
        		fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[i],Wr[i],Bx, &x3Flux[k][j][i]);

#ifdef FIRST_ORDER_FLUX_CORRECTION
/* revert to predictor flux if this flux NaN'ed */
        if ((x3Flux[k][j][i].d  != x3Flux[k][j][i].d)  ||
#ifndef BAROTROPIC
            (x3Flux[k][j][i].E  != x3Flux[k][j][i].E)  ||
#endif
#if defined(MHD) || defined(RADIATION_MHD)
            (x3Flux[k][j][i].By != x3Flux[k][j][i].By) ||
            (x3Flux[k][j][i].Bz != x3Flux[k][j][i].Bz) ||
#endif
            (x3Flux[k][j][i].Mx != x3Flux[k][j][i].Mx) ||
            (x3Flux[k][j][i].My != x3Flux[k][j][i].My) ||
            (x3Flux[k][j][i].Mz != x3Flux[k][j][i].Mz)) {
          x3Flux[k][j][i] = x3FluxP[k][j][i];
          NaNFlux++;
        }
#endif /* End FOFC */

      			}
    		}
  	}

#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (NaNFlux != 0) {
    printf("[Step10] %i second-order fluxes replaced\n",NaNFlux);
    NaNFlux=0;
  }
#endif

	
/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */

#ifdef RADIATION_MHD
  	integrate_emf1_corner(pG);
  	integrate_emf2_corner(pG);
  	integrate_emf3_corner(pG);
#endif	
	
	
	/* Remap Ey at is and ie+1 to conserve Bz in shearing box */
#ifdef SHEARING_BOX
    get_myGridIndex(pD, myID_Comm_world, &my_iproc, &my_jproc, &my_kproc);
	
/* initialize remapped Fluxes */


    for(k=ks; k<=ke+1; k++) {
      for(j=js; j<=je+1; j++){
        Flxiib[k][j].d = x1Flux[k][j][is].d;
        Flxiib[k][j].M1 = x1Flux[k][j][is].Mx;
        Flxiib[k][j].M2 = x1Flux[k][j][is].My;
        Flxiib[k][j].M3 = x1Flux[k][j][is].Mz;

        Flxoib[k][j].d = x1Flux[k][j][ie+1].d;
        Flxoib[k][j].M1 = x1Flux[k][j][ie+1].Mx;
        Flxoib[k][j].M2 = x1Flux[k][j][ie+1].My;
        Flxoib[k][j].M3 = x1Flux[k][j][ie+1].Mz;
#ifndef BAROTROPIC
        Flxiib[k][j].E = x1Flux[k][j][is].E;
        Flxoib[k][j].E = x1Flux[k][j][ie+1].E;
#endif
#ifdef RADIATION_MHD
        Flxiib[k][j].B1c = emf1[k][j][is];
        Flxiib[k][j].B2c = emf2[k][j][is];
        Flxiib[k][j].B3c = emf3[k][j][is];

        Flxoib[k][j].B1c = emf1[k][j][ie+1];
        Flxoib[k][j].B2c = emf2[k][j][ie+1];
        Flxoib[k][j].B3c = emf3[k][j][ie+1];
#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) {
          Flxiib[k][j].s[n] = x1Flux[k][j][is].s[n];
          Flxoib[k][j].s[n] = x1Flux[k][j][ie+1].s[n];
        }
#endif
      }
    }


/* compute remapped Fluxes from opposite side of grid */

    if (my_iproc == 0) {

      RemapFlx_ix1(pD, Flxiib, Flxoib, rFlxiib);
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      RemapFlx_ox1(pD, Flxiib, Flxoib, rFlxoib);
    }

/* Now average fluxes and remapped fluxes */
	

    if (my_iproc == 0) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
#ifdef RADIATION_MHD
          emf2[k][j][is] = 0.5*(emf2[k][j][is] + rFlxiib[k][j].B2c);
#endif
         x1Flux[k][j][is].d  = 0.5*(x1Flux[k][j][is].d  + rFlxiib[k][j].d);
          x1Flux[k][j][is].Mx = 0.5*(x1Flux[k][j][is].Mx + rFlxiib[k][j].M1);
          x1Flux[k][j][is].My = 0.5*(x1Flux[k][j][is].My + rFlxiib[k][j].M2);
          x1Flux[k][j][is].Mz = 0.5*(x1Flux[k][j][is].Mz + rFlxiib[k][j].M3);
        }
      }
    }

    if (my_iproc == (pD->NGrid[0]-1)) {
      for(k=ks; k<=ke+1; k++) {
        for(j=js; j<=je; j++){
#ifdef RADIATION_MHD
          emf2[k][j][ie+1] = 0.5*(emf2[k][j][ie+1] + rFlxoib[k][j].B2c);
#endif
          x1Flux[k][j][ie+1].d =0.5*(x1Flux[k][j][ie+1].d  + rFlxoib[k][j].d);
          x1Flux[k][j][ie+1].Mx=0.5*(x1Flux[k][j][ie+1].Mx + rFlxoib[k][j].M1);
          x1Flux[k][j][ie+1].My=0.5*(x1Flux[k][j][ie+1].My + rFlxoib[k][j].M2);
          x1Flux[k][j][ie+1].Mz=0.5*(x1Flux[k][j][ie+1].Mz + rFlxoib[k][j].M3);
        }
      }
    }

#endif /* SHEARING_BOX */	
	
	

	
	
/*=======================================================*/	
	

/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */
#ifdef RADIATION_MHD
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

			/* First, assume this is good cell */
			badcellflag = 0;

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

			Sigma[0] = Sigma_sF;
			Sigma[1] = Sigma_aF;
			Sigma[2] = Sigma_aP;
			Sigma[3] = Sigma_aE;


			if(pressure > 2.0 * TINY_NUMBER){		
				diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;
			}
			else{
				diffTEr = 0.0;
			}

			/* This is the prefactor used in SEE, if this is smaller than 1, may cause trouble */
			if(4.0 * Prat * temperature * temperature * temperature * (Gamma - 1.0)/(density * R_ideal) < 1.0)
				Sourceflag = 0;
			else
				Sourceflag = 1;
			

					
#ifdef FARGO
			
			/* Include background shearing in Usource, which is only used in dSource */

			velocity_y -= qom * x1;						
					
#endif	
					
			/* co-moving flux */
			Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
			Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
			Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;
					
		
		
			/* Part of momentum source term */
			Smx0 = -Prat * velocity_x * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;
			Smy0 = -Prat * velocity_y * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;
			Smz0 = -Prat * velocity_z * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * Usource.Er) / Crat;
			
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);
	
		
			/*=========================================================*/
			/* In case velocity is large and momentum source is stiff */
			SFmx = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er / (density * Crat);
			/*      + diffTEr / (density * Crat);
			*/
					
			SFmy = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er / (density * Crat);
			/*      + diffTEr / (density * Crat);
			*/
					
			SFmz = (Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er / (density * Crat);
			/*      + diffTEr / (density * Crat);
			*/
					
					
			Source_Inv[0][0] = 1.0;
			Source_Inv[1][1] = 1.0 / (1.0 + dt * Prat * SFmx);
			Source_Inv[2][2] = 1.0 / (1.0 + dt * Prat * SFmy);
			Source_Inv[3][3] = 1.0 / (1.0 + dt * Prat * SFmz);
					
			/*=========================================================*/
					
	
			Source_Inv[4][0] = 0.0;
			Source_Inv[4][1] = 0.0;
			Source_Inv[4][2] = 0.0;
			Source_Inv[4][3] = 0.0;
			Source_Inv[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
	

			/*=======================================================*/
			/****************************************/
			/* Modify the energy source term to include the stiffness of the momentum source term */
			if(Erflag && Sourceflag){
					Source[k][j][i][4] = -Prat * Crat * (diffTEr + (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x * Source_Inv[1][1] + velocity_y * Fr0y * Source_Inv[2][2] + velocity_z * Fr0z * Source_Inv[3][3])/Crat);
			}
			else{
					Source[k][j][i][4] = - Prat * (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x * Source_Inv[1][1] + velocity_y * Fr0y * Source_Inv[2][2] + velocity_z * Fr0z * Source_Inv[3][3]);
			}
					
			Prwork1 = -Prat * (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x * Source_Inv[1][1] + velocity_y * Fr0y * Source_Inv[2][2] + velocity_z * Fr0z * Source_Inv[3][3]);
					
					
					
			/*****************************************************/
					



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

		Source[k][j][i][0] = 0.0;
		for(n=0; n<5; n++) {
				tempguess[n] = dt * Source_Inv[n][n] * (Source[k][j][i][n] - divFlux1[n] - divFlux2[n] - divFlux3[n]);			
		}
					
		if((!Erflag) || (!Sourceflag)){
			tempguess[4] += -Prat * pG->Ersource[k][j][i];
						
		}
					
				
			Prworksource = dt * Source_Inv[4][4] * Prwork1;
					
		Uguess[0] = pG->U[k][j][i].d  + tempguess[0];
		Uguess[1] = pG->U[k][j][i].M1 + tempguess[1] + dt * Smx0;
		Uguess[2] = pG->U[k][j][i].M2 + tempguess[2] + dt * Smy0;
		Uguess[3] = pG->U[k][j][i].M3 + tempguess[3] + dt * Smz0;
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

		ShearSource[3] = -(dtodx1*(x1Flux[k][j][i  ].d*(phic - phil) +
				x1Flux[k][j][i+1].d*(phir - phic)));
					
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
/* Finish the predict step. All the source terms are added */
		
		density    = Uguess[0];
		velocity_x = Uguess[1] / density;
		velocity_y = Uguess[2] / density;
		velocity_z = Uguess[3] / density;

		velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);



		pressure = (Uguess[4] - 0.5 * (density * velocity * velocity)) * (Gamma - 1.0);
		/* Should include magnetic energy for MHD */
#ifdef RADIATION_MHD
		/* Should use the updated cell centered magnetic field to calculate the temperature */
		B1ch = 0.5*(    pG->B1i[k][j][i] +     pG->B1i[k][j][i+1]);
		B2ch = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
		B3ch = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
		
		pressure -= 0.5 * (B1ch * B1ch + B2ch * B2ch + B3ch * B3ch) * (Gamma - 1.0);
#endif

		if((pressure < TINY_NUMBER) || (pressure != pressure)){
			if(density < dfloor){
				/* keep original temperature */
				pressure = pG->U[k][j][i].d * temperature * R_ideal;
			}
			else{
				pressure = density * temperature * R_ideal;
			}
			badcellflag = 1;

		}

		
		/* temperature and density are original temperature and density */
		if((density < dfloor) || (density != density)){
			Uguess[0] = pG->U[k][j][i].d;
			density = pG->U[k][j][i].d;
			Uguess[1] = density * velocity_x;
			Uguess[2] = density * velocity_y;
			Uguess[3] = density * velocity_z;
			
			
			badcellflag = 1;
		}

		/* for bad cell, need to recalculate total energy */
		if(badcellflag){
			Uguess[4] =  pressure / (Gamma - 1.0) 
				+ 0.5 * density * (velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
#ifdef RADIATION_MHD
                        Uguess[4] += 0.5 * (B1ch * B1ch + B2ch * B2ch + B3ch * B3ch);
#endif


		}		

		
		temperature = pressure / (density * R_ideal);
					
#ifdef FARGO
			/* With FARGO, we should add background shearing to the source terms */
					
			/* Include background shearing in Usource, which is only used in dSource */
					
			velocity_y -= qom * x1;	
					
#endif	
			
					
		/* co-moving flux */
		/* Correct the radiation work term with original opacity, but with updated velocity */			
		Fr0x = Usource.Fr1 - ((1.0 + Usource.Edd_11) * velocity_x + Usource.Edd_21 * velocity_y + Usource.Edd_31 * velocity_z) * Usource.Er / Crat;
		Fr0y = Usource.Fr2 - ((1.0 + Usource.Edd_22) * velocity_y + Usource.Edd_21 * velocity_x + Usource.Edd_32 * velocity_z) * Usource.Er / Crat;
		Fr0z = Usource.Fr3 - ((1.0 + Usource.Edd_33) * velocity_z + Usource.Edd_31 * velocity_x + Usource.Edd_32 * velocity_y) * Usource.Er / Crat;
					
		Prwork2 = -Prat * (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x * Source_Inv[1][1] + velocity_y * Fr0y * Source_Inv[2][2] + velocity_z * Fr0z * Source_Inv[3][3]);
					
					/* Source term for momentum */
		Source_guess[1] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0x); /* + velocity_x * diffTEr / Crat); */
		Source_guess[2] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0y); /* + velocity_y * diffTEr / Crat); */
		Source_guess[3] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0z); /* + velocity_z * diffTEr / Crat); */
					
					

		/* If Opacity is not set, Sigma_? will not be changed. */
		/* Negative pressur is handled in the opacity function */
					
		/* update opacity for the thermalization term, not for momentum term, which can cause trouble */	

		/* Prepare the Prims variable */
		Wopacity = Cons_to_Prim(&pG->U[k][j][i]);
		/* Now update the density, pressure and velocity */
		Wopacity.d = density;
		Wopacity.P = pressure;
		Wopacity.V1 = velocity_x;
		Wopacity.V2 = velocity_y;
		Wopacity.V3 = velocity_z;
		/* background shearing should be included */
		
	
		if(Opacity != NULL){
			Opacity(&Wopacity, Sigma, NULL);
		
			Sigma_sF = Sigma[0];
			Sigma_aF = Sigma[1];
			Sigma_aP = Sigma[2];
			Sigma_aE = Sigma[3];
			
			
		}
		
		/* Update the sourceflag, It can happen that initially it is radiation dominated, but become gas *
		  * pressure dominated in the predict step */
		if(4.0 * Prat * temperature * temperature * temperature * (Gamma - 1.0)/(density * R_ideal) < 1.0)
				Sourceflag2 = 0;
		else
				Sourceflag2 = 1;
			
			
		
		/* calculate the predict energy source term */
		ThermalRelaxation(temperature, pG->U[k][j][i].Er, density, Sigma_aP, Sigma_aE, dt, NULL, &Ersource);
		Ersource = Ersource - pG->U[k][j][i].Er;

		diffTEr = Sigma_aP * pow(temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er;
		
		

		/* update source term */
		Usource.d  = Uguess[0];
		Usource.Mx = Uguess[1];
		Usource.My = Uguess[2];
		Usource.Mz = Uguess[3];
		Usource.E  = Uguess[4];

		for(m=0; m<NOPACITY;m++)
			Usource.Sigma[m] = Sigma[m];

					


		/* The Source term */
		/* Only do this if density and pressure are normal */
		if((Usource.d > 0.0) && (pressure > 2.0 * TINY_NUMBER)){
			dSource(Usource, Bx, &SEE, &SErho, &SEmx, &SEmy, &SEmz, x1);
		}
		
		/* Do not need to update Source_Inv[1][1], Source_Inv[2][2], Source_Inv[3][3] */
		/* We only need SEE here, in principle should use the updated Bx */	
		Source_Invnew[4][4] = 1.0 / (1.0 + dt * Prat * Crat * SEE);
		
		if(!Sourceflag2)
                       Source_Invnew[4][4] = Source_Inv[4][4];

	
		/* Source term for total Energy */
		if(Erflag && Sourceflag2){
			Source_guess[4] = -Prat * Crat * diffTEr + Prwork2;
		}
		else{
			Source_guess[4] =  Prwork2;
		}
			
			
			
			/* Calculate the error term */
			/* Operator split terms don't need to be included here */
			Errort[1] = hdt * (Source[k][j][i][1] + Source_guess[1]) 
						- dt * (divFlux1[1] + divFlux2[1] + divFlux3[1]) 
						- dt * Source_Inv[1][1] * (Source[k][j][i][1] - divFlux1[1] - divFlux2[1] - divFlux3[1]);
			Errort[2] = hdt * (Source[k][j][i][2] + Source_guess[2]) 
						- dt * (divFlux1[2] + divFlux2[2] + divFlux3[2]) 
						- dt * Source_Inv[2][2] * (Source[k][j][i][2] - divFlux1[2] - divFlux2[2] - divFlux3[2]);
			Errort[3] = hdt * (Source[k][j][i][3] + Source_guess[3]) 
						- dt * (divFlux1[3] + divFlux2[3] + divFlux3[3]) 
						- dt * Source_Inv[3][3] * (Source[k][j][i][3] - divFlux1[3] - divFlux2[3] - divFlux3[3]);
			
			for(m=1; m<4; m++) {
				tempguess[m] = Source_Inv[m][m] * Errort[m];
			}
			
			/* Correction to the energy source term needs special treatment */
			Errort[4] = hdt * (Source[k][j][i][4] + Source_guess[4]) - dt * (divFlux1[4] + divFlux2[4] + divFlux3[4]) 
						- dt * Source_Inv[4][4] * (Source[k][j][i][4] - (divFlux1[4] + divFlux2[4] + divFlux3[4]));
			tempguess[4] = Source_Invnew[4][4] * Errort[4];
			
			if(Erflag){
				if((Sourceflag) && (!Sourceflag2)){
					tempguess[4] += -0.5 * Prat * Ersource;
				}
				else if((!Sourceflag) && (!Sourceflag2)){
					tempguess[4] += -0.5 * Prat * (Ersource - pG->Ersource[k][j][i]);
				}
				else if((!Sourceflag) && (Sourceflag2)){
					tempguess[4] += 0.5 * Prat * pG->Ersource[k][j][i];	
				}
				
			}
			else {
					tempguess[4] += -0.5 * Prat * (Ersource - pG->Ersource[k][j][i]);
			}

			
			/* This is the actual added radiation work term */
			Prworksource += Source_Invnew[4][4] * (hdt * (Prwork1 + Prwork2) - Prworksource);
			
			
				
			/* Apply the correction */
			/* Estimate the added radiation source term  */
			
			/* If Erflag == 0 , only add source term as first order */
			if(Erflag){
				if(Prat > 0.0){
					pG->Ersource[k][j][i] = Uguess[4] + tempguess[4] 
										- (pG->U[k][j][i].E - dt * (divFlux1[4] + divFlux2[4] + divFlux3[4]));
					/* Energy error seems smaller with Errot term added */
#ifdef SHEARING_BOX
					pG->Ersource[k][j][i] -= ShearSource[3];

#endif

#ifdef CONS_GRAVITY
			        pG->Ersource[k][j][i] -= 0.5*(pG->U[k][j][i].d-grav_mean_rho)*pG->Phi_old[k][j][i]-0.5*(density_old[k][j][i]-grav_mean_rho)*pG->Phi[k][j][i];
#endif

					/* subtract the radiation work term, which is not seperated */
				
					/* Radiation work term is not seperated. Will correct error due to this term later */
					pG->Ersource[k][j][i] -= Prworksource;
				

					pG->Ersource[k][j][i] /= -Prat;
							
				

				}
				else{				

					pG->Ersource[k][j][i] = 0.0;
				}/* End prat = 0 */
			} /* End Erflag */
			else {
				
				pG->Ersource[k][j][i] = 0.5 * (pG->Ersource[k][j][i] + Ersource);
			}


			pG->Eulersource[k][j][i] = -Prworksource/Prat;

			if(badcellflag){

				/* Do not apply correction for bad cell */
				pG->Ersource[k][j][i] = 0.0;
				pG->Eulersource[k][j][i] = 0.0;
				for(n=0; n<5; n++){
					tempguess[n] = 0.0;
				}		


			}
					
			/* The momentum source term Prat * v * Sigma(T^4 - Er)/Crat is added explicitly */
			/* In case overshoot make negative pressure */
			/* Apply the correction */

			pG->U[k][j][i].d  = Uguess[0];
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
					
/* 
				if(pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c < TINY_NUMBER)
				{
						
					pG->U[k][j][i].M1 = 0.0;
					pG->U[k][j][i].M2 = 0.0;
					pG->U[k][j][i].M3 = 0.0;
						
				 }
*/	
					
      			}
    		}
  	}
#endif /* RADIATION MHD */



#ifdef STATIC_MESH_REFINEMENT
/* ======== This is adopted from integrate_3d_ctu.c ================== 
 * With SMR, store fluxes at boundaries of child and parent grids.  */
/* Loop over all child grids ------------------*/

  	for (ncg=0; ncg<pG->NCGrid; ncg++) {

/* x1-boundaries of child Grids (interior to THIS Grid) */

    		for (dim=0; dim<2; dim++){
      			if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        			if (dim==0) i = pG->CGrid[ncg].ijks[0];
        			if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;
        			jcs = pG->CGrid[ncg].ijks[1];
        			jce = pG->CGrid[ncg].ijke[1];
        			kcs = pG->CGrid[ncg].ijks[2];
        			kce = pG->CGrid[ncg].ijke[2];

        			for (k=kcs, kk=0; k<=kce; k++, kk++){
          			for (j=jcs, jj=0; j<=jce; j++, jj++){
            				pG->CGrid[ncg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d; 
            				pG->CGrid[ncg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx; 
            				pG->CGrid[ncg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            				pG->CGrid[ncg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz; 
#ifndef BAROTROPIC
            				pG->CGrid[ncg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->CGrid[ncg].myFlx[dim][kk][jj].B1c = 0.0;
            				pG->CGrid[ncg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By; 
            				pG->CGrid[ncg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
            				for (n=0; n<NSCALARS; n++)
              					pG->CGrid[ncg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n]; 
#endif
          			}/* end j */
        			}/* end k */
#if defined(MHD) || defined(RADIATION_MHD)
        			for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          			for (j=jcs, jj=0; j<=jce; j++, jj++){
            				pG->CGrid[ncg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          			}/* end j */
        			}/* end k */

        			for (k=kcs, kk=0; k<=kce; k++, kk++){
          			for (j=jcs, jj=0; j<=jce+1; j++, jj++){
            				pG->CGrid[ncg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          			}/* end j */
        			}/* end k */
#endif /* MHD */
      			}/* end if flux is not null */
    		}/* end dim from 0 to 2 */

/* x2-boundaries of child Grids (interior to THIS Grid) */

    		for (dim=2; dim<4; dim++){
      			if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        			ics = pG->CGrid[ncg].ijks[0];
        			ice = pG->CGrid[ncg].ijke[0];
        			if (dim==2) j = pG->CGrid[ncg].ijks[1];
        			if (dim==3) j = pG->CGrid[ncg].ijke[1] + 1;
        			kcs = pG->CGrid[ncg].ijks[2];
        			kce = pG->CGrid[ncg].ijke[2];

        			for (k=kcs, kk=0; k<=kce; k++, kk++){
          			for (i=ics, ii=0; i<=ice; i++, ii++){
            				pG->CGrid[ncg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d; 
            				pG->CGrid[ncg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz; 
            				pG->CGrid[ncg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            				pG->CGrid[ncg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My; 
#ifndef BAROTROPIC
            				pG->CGrid[ncg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->CGrid[ncg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz; 
            				pG->CGrid[ncg].myFlx[dim][kk][ii].B2c = 0.0;
            				pG->CGrid[ncg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By; 
#endif /* MHD */
#if (NSCALARS > 0)
            			for (n=0; n<NSCALARS; n++)
              				pG->CGrid[ncg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n]; 
#endif
          			}/* end i */
        			}/* end k */
#if defined(MHD) || defined(RADIATION_MHD)
        			for (k=kcs, kk=0; k<=kce+1; k++, kk++){
          			for (i=ics, ii=0; i<=ice; i++, ii++){
            				pG->CGrid[ncg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          			}/* end i */
        			}/* end k */
        			for (k=kcs, kk=0; k<=kce; k++, kk++){
          			for (i=ics, ii=0; i<=ice+1; i++, ii++){
            				pG->CGrid[ncg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          			}/* end i */
        			}/* end k */
#endif /* MHD */
      			}/* end if mflux is not NULL */
    		}/* end m from 2 to 3 */

/* x3-boundaries of child Grids (interior to THIS Grid) */

    		for (dim=4; dim<6; dim++){
      			if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        			ics = pG->CGrid[ncg].ijks[0];
        			ice = pG->CGrid[ncg].ijke[0];
        			jcs = pG->CGrid[ncg].ijks[1];
        			jce = pG->CGrid[ncg].ijke[1];
        			if (dim==4) k = pG->CGrid[ncg].ijks[2];
        			if (dim==5) k = pG->CGrid[ncg].ijke[2] + 1;

        			for (j=jcs, jj=0; j<=jce; j++, jj++){
          			for (i=ics, ii=0; i<=ice; i++, ii++){
           	 			pG->CGrid[ncg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d; 
            				pG->CGrid[ncg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My; 
            				pG->CGrid[ncg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            				pG->CGrid[ncg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx; 
#ifndef BAROTROPIC
            				pG->CGrid[ncg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->CGrid[ncg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By; 
            				pG->CGrid[ncg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz; 
            				pG->CGrid[ncg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            			for (n=0; n<NSCALARS; n++)
              				pG->CGrid[ncg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n]; 
#endif
          			}/* end i */
        			}/* end j */
#if defined(MHD) || defined(RADIATION_MHD)
        			for (j=jcs, jj=0; j<=jce+1; j++, jj++){
          			for (i=ics, ii=0; i<=ice; i++, ii++){
            				pG->CGrid[ncg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          			}/* end i */
        			}/* end j */
        			for (j=jcs, jj=0; j<=jce; j++, jj++){
          			for (i=ics, ii=0; i<=ice+1; i++, ii++){
            				pG->CGrid[ncg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          			}/* end i */
        			}/* end j */
#endif /* MHD */
      			}/* end if myFlux is not null */
    		}/* end dim 4 and 5 */
  	} /* end loop over child Grids */

/* Loop over all parent grids ------------------*/

  	for (npg=0; npg<pG->NPGrid; npg++) {

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */

    		for (dim=0; dim<2; dim++){
      			if (pG->PGrid[npg].myFlx[dim] != NULL) {

        			if (dim==0) i = pG->PGrid[npg].ijks[0];
        			if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;
        			jps = pG->PGrid[npg].ijks[1];
        			jpe = pG->PGrid[npg].ijke[1];
        			kps = pG->PGrid[npg].ijks[2];
        			kpe = pG->PGrid[npg].ijke[2];

        			for (k=kps, kk=0; k<=kpe; k++, kk++){
          			for (j=jps, jj=0; j<=jpe; j++, jj++){
            				pG->PGrid[npg].myFlx[dim][kk][jj].d  = x1Flux[k][j][i].d; 
            				pG->PGrid[npg].myFlx[dim][kk][jj].M1 = x1Flux[k][j][i].Mx; 
            				pG->PGrid[npg].myFlx[dim][kk][jj].M2 = x1Flux[k][j][i].My;
            				pG->PGrid[npg].myFlx[dim][kk][jj].M3 = x1Flux[k][j][i].Mz; 
#ifndef BAROTROPIC
            				pG->PGrid[npg].myFlx[dim][kk][jj].E  = x1Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->PGrid[npg].myFlx[dim][kk][jj].B1c = 0.0;
            				pG->PGrid[npg].myFlx[dim][kk][jj].B2c = x1Flux[k][j][i].By; 
            				pG->PGrid[npg].myFlx[dim][kk][jj].B3c = x1Flux[k][j][i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
            				for (n=0; n<NSCALARS; n++)
              					pG->PGrid[npg].myFlx[dim][kk][jj].s[n]  = x1Flux[k][j][i].s[n]; 
#endif
          			}/* end j */
        			}/* end k */
#if defined(MHD) || defined(RADIATION_MHD)
        			for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          			for (j=jps, jj=0; j<=jpe; j++, jj++){
            				pG->PGrid[npg].myEMF2[dim][kk][jj] = emf2[k][j][i];
          			}/* end j */
        			}/* end k */
        			for (k=kps, kk=0; k<=kpe; k++, kk++){
          			for (j=jps, jj=0; j<=jpe+1; j++, jj++){
            				pG->PGrid[npg].myEMF3[dim][kk][jj] = emf3[k][j][i];
          			}/* end j*/
        			}/* end  k */
#endif /* MHD */
      			}/* end mflux is not null */
    		}/* end dim 0 to 1 */

/* x2-boundaries of parent Grids (at boundaries of THIS Grid)  */

    		for (dim=2; dim<4; dim++){
      			if (pG->PGrid[npg].myFlx[dim] != NULL) {

        			ips = pG->PGrid[npg].ijks[0];
        			ipe = pG->PGrid[npg].ijke[0];
        			if (dim==2) j = pG->PGrid[npg].ijks[1];
        			if (dim==3) j = pG->PGrid[npg].ijke[1] + 1;
        			kps = pG->PGrid[npg].ijks[2];
        			kpe = pG->PGrid[npg].ijke[2];

        			for (k=kps, kk=0; k<=kpe; k++, kk++){
          			for (i=ips, ii=0; i<=ipe; i++, ii++){
            				pG->PGrid[npg].myFlx[dim][kk][ii].d  = x2Flux[k][j][i].d; 
            				pG->PGrid[npg].myFlx[dim][kk][ii].M1 = x2Flux[k][j][i].Mz; 
            				pG->PGrid[npg].myFlx[dim][kk][ii].M2 = x2Flux[k][j][i].Mx;
            				pG->PGrid[npg].myFlx[dim][kk][ii].M3 = x2Flux[k][j][i].My; 
#ifndef BAROTROPIC
            				pG->PGrid[npg].myFlx[dim][kk][ii].E  = x2Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->PGrid[npg].myFlx[dim][kk][ii].B1c = x2Flux[k][j][i].Bz; 
            				pG->PGrid[npg].myFlx[dim][kk][ii].B2c = 0.0;
            				pG->PGrid[npg].myFlx[dim][kk][ii].B3c = x2Flux[k][j][i].By; 
#endif /* MHD */
#if (NSCALARS > 0)
            			for (n=0; n<NSCALARS; n++)
              				pG->PGrid[npg].myFlx[dim][kk][ii].s[n]  = x2Flux[k][j][i].s[n]; 
#endif
          			}/* end i*/
        			}/* end k*/
#if defined(MHD) || defined(RADIATION_MHD)
        			for (k=kps, kk=0; k<=kpe+1; k++, kk++){
          			for (i=ips, ii=0; i<=ipe; i++, ii++){
            				pG->PGrid[npg].myEMF1[dim][kk][ii] = emf1[k][j][i];
          			}/* end i*/
        			}/* end k */
        			for (k=kps, kk=0; k<=kpe; k++, kk++){
          			for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            				pG->PGrid[npg].myEMF3[dim][kk][ii] = emf3[k][j][i];
          			}/* end i */
        			}/* end k*/
#endif /* MHD */
      			}/* end if myflux not null */
    		}/* end dim 2 to 3 */

/* x3-boundaries of parent Grids (at boundaries of THIS Grid)  */

    		for (dim=4; dim<6; dim++){
      			if (pG->PGrid[npg].myFlx[dim] != NULL) {

        			ips = pG->PGrid[npg].ijks[0];
        			ipe = pG->PGrid[npg].ijke[0];
        			jps = pG->PGrid[npg].ijks[1];
        			jpe = pG->PGrid[npg].ijke[1];
        			if (dim==4) k = pG->PGrid[npg].ijks[2];
        			if (dim==5) k = pG->PGrid[npg].ijke[2] + 1;

        			for (j=jps, jj=0; j<=jpe; j++, jj++){
         	 		for (i=ips, ii=0; i<=ipe; i++, ii++){
            				pG->PGrid[npg].myFlx[dim][jj][ii].d  = x3Flux[k][j][i].d; 
            				pG->PGrid[npg].myFlx[dim][jj][ii].M1 = x3Flux[k][j][i].My; 
            				pG->PGrid[npg].myFlx[dim][jj][ii].M2 = x3Flux[k][j][i].Mz;
            				pG->PGrid[npg].myFlx[dim][jj][ii].M3 = x3Flux[k][j][i].Mx; 
#ifndef BAROTROPIC
            				pG->PGrid[npg].myFlx[dim][jj][ii].E  = x3Flux[k][j][i].E; 
#endif /* BAROTROPIC */
#if defined(MHD) || defined(RADIATION_MHD)
            				pG->PGrid[npg].myFlx[dim][jj][ii].B1c = x3Flux[k][j][i].By; 
            				pG->PGrid[npg].myFlx[dim][jj][ii].B2c = x3Flux[k][j][i].Bz; 
            				pG->PGrid[npg].myFlx[dim][jj][ii].B3c = 0.0;
#endif /* MHD */
#if (NSCALARS > 0)
            				for (n=0; n<NSCALARS; n++)
              					pG->PGrid[npg].myFlx[dim][jj][ii].s[n]  = x3Flux[k][j][i].s[n]; 
#endif
          			}/* end i */
        			}/* end j*/
#if defined(MHD) || defined(RADIATION_MHD)
        			for (j=jps, jj=0; j<=jpe+1; j++, jj++){
          			for (i=ips, ii=0; i<=ipe; i++, ii++){
            				pG->PGrid[npg].myEMF1[dim][jj][ii] = emf1[k][j][i];
          			}/* end i*/
        			}/* end j*/
        			for (j=jps, jj=0; j<=jpe; j++, jj++){
          			for (i=ips, ii=0; i<=ipe+1; i++, ii++){
            				pG->PGrid[npg].myEMF2[dim][jj][ii] = emf2[k][j][i];
          			}/* end i*/
        			}/* end j*/
#endif /* MHD */
      			}/* end flux is not null */
    		}/* end dim 4 to 5 */
  	}/* end loop parent grids */

#endif /* STATIC_MESH_REFINEMENT */
	
	
/* Add density floor and beta floor */
/*	for(k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=is; i<=ie; i++) {
	 
				badcellflag = 0;
	 
				velocity_x = pG->U[k][j][i].M1 / pG->U[k][j][i].d;
				velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
				velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
 
				velocity = sqrt(velocity_x * velocity_x + velocity_y * velocity_y + velocity_z * velocity_z);
			 
				Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
	 
				temperature = Wtemp.P / (Wtemp.d * R_ideal);
#ifdef RADIATION_MHD	 
				Bpre = 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif	 
				if((Wtemp.P < 2.0 * TINY_NUMBER) || (temperature < Tfloor)){
					if(pG->U[k][j][i].d < dfloor){
						Wtemp.P = dfloor * R_ideal * Tfloor;
					}
					else{
						Wtemp.P = pG->U[k][j][i].d * R_ideal * Tfloor;
					}
						temperature = Tfloor;
						badcellflag = 1;
					}
	 
	 
				if(pG->U[k][j][i].d < dfloor){
	 
					pG->U[k][j][i].d = dfloor;
	 
					Wtemp.d =  dfloor;
					Wtemp.P = Wtemp.d * temperature * R_ideal;
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
	
	/*====================================================================*/
	/* Add Compton scattering source term at the end of the integrator */	
	
	/* calculate pG->Tguess after magnetic field is updated */
	
		
/*	 for (k=ks; k<=ke; k++) {
		for (j=js; j<=je; j++) {
			for (i=is; i<=ie; i++) {
	 
				pressure = pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	 
				pressure *= (Gamma - 1.0);
	 
				temperature = pressure / (pG->U[k][j][i].d * R_ideal);
	 
	 
			if(fabs(pow(temperature, 4.0) - pG->U[k][j][i].Er) < TINY_NUMBER){
				pG->Comp[k][j][i] = 0.0;
			}
			else{
				if(temperature > TINY_NUMBER){
					Tr = pow(pG->U[k][j][i].Er, 0.25);
					coefA = 4.0 * dt * Crat * pG->U[k][j][i].Sigma[0] / (T_e/T0);
					coefK = (Gamma - 1.0) * Prat / (R_ideal * pG->U[k][j][i].d);
					coefB = temperature + coefK * pG->U[k][j][i].Er;
					coef1 = coefA * coefK;
					coef2 = coefA;
					coef3 = 1.0 - coefA * coefB;
					coef4 = -pG->U[k][j][i].Er;
	 
					if(Tr < temperature){
						Tr = rtsafe(Tcompton, Tr * (1.0 - 0.01), temperature * (1.0 + 0.01), 1.e-14, coef1, coef2, coef3, coef4);
					}
					else{
	 
						Tr = rtsafe(Tcompton, temperature * (1.0 - 0.01), Tr * (1.0 + 0.01), 1.e-14, coef1, coef2, coef3, coef4);
					}
	 
					pG->Comp[k][j][i] = pow(Tr, 4.0) - pG->U[k][j][i].Er;
					pG->U[k][j][i].E += -Prat * pG->Comp[k][j][i];
	 
	 
				}
				else{
					pressure = pow(pG->U[k][j][i].Er, 0.25) * pG->U[k][j][i].d * R_ideal;
	 
					pG->U[k][j][i].E = pressure / (Gamma - 1.0) +  0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 + pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / pG->U[k][j][i].d;
#ifdef RADIATION_MHD
					pG->U[k][j][i].E += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif
	 
					pG->Comp[k][j][i] = 0.0;
		
				} 
			}
		}
	 }
	 }
*/
	/* Boundary condition is applied in the main function */
	/* Check the pressure to make sure that it is positive */
	/* This will make the code more robust */
/* If cell-centered d or P have gone negative, or if v^2 > 1 in SR, correct
 * by using 1st order predictor fluxes */
#ifdef FIRST_ORDER_FLUX_CORRECTION

	for (k=ks; k<=ke; k++){
		for (j=js; j<=je; j++) {
    			for (i=is; i<=ie; i++){
 		Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
	if (Wtemp.d < 0.0) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negd++;
        }

        if (Wtemp.P < 2.0 * TINY_NUMBER) {
          flag_cell = 1;
          BadCell.i = i;
          BadCell.j = j;
          BadCell.k = k;
          negP++;
	
        }

	
	if (flag_cell != 0 && (negd >0 || negP > 0)) {
		cc_pos(pG,i,j,k,&x1,&x2,&x3);
		if(fabs(x3) < fabs(zmin))
			zmin = x3;


	  if(Prat > 0.0){
		/*
	  	pG->Tguess[k][j][i] = -pG->dt * Source[k][j][i][4] / Prat;	
		*/
		/* For bad cells, do not include radiation source term */
		pG->Tguess[k][j][i] = pG->U[k][j][i].Er;

	  }	  
          FixCell(pG, BadCell);
          flag_cell=0;
        
	/* Check whether fix this issue */
	 Wtemp = Cons_to_Prim(&(pG->U[k][j][i]));
	if(Wtemp.P < 2.0 * TINY_NUMBER) {
		/* Do average */
		Np = 0;
		Wtemp1 = Cons_to_Prim(&(pG->U[k][j][i-1]));
			if(Wtemp1.P > 2.0 * TINY_NUMBER)	Np++;
		Wtemp2 = Cons_to_Prim(&(pG->U[k][j][i+1]));
			if(Wtemp2.P > 2.0 * TINY_NUMBER)	Np++;
		Wtemp3 = Cons_to_Prim(&(pG->U[k][j-1][i]));
			if(Wtemp3.P > 2.0 * TINY_NUMBER)	Np++;
		Wtemp4 = Cons_to_Prim(&(pG->U[k][j+1][i]));
			if(Wtemp4.P > 2.0 * TINY_NUMBER)	Np++;
		Wtemp5 = Cons_to_Prim(&(pG->U[k-1][j][i]));
			if(Wtemp5.P > 2.0 * TINY_NUMBER)	Np++;
		Wtemp6 = Cons_to_Prim(&(pG->U[k+1][j][i]));
			if(Wtemp6.P > 2.0 * TINY_NUMBER)	Np++;
		
		if(Np > 0){
			Wtemp.P = (Wtemp1.P + Wtemp2.P + Wtemp3.P + Wtemp4.P + Wtemp5.P + Wtemp6.P) / Np;
		}
		else{
			Wtemp.P = INIdata[k][j][i][5];		
		}
		pG->U[k][j][i].E = Wtemp.P / (Gamma - 1.0) + 0.5 * Wtemp.d * (Wtemp.V1 * Wtemp.V1 + Wtemp.V2 * Wtemp.V2 + Wtemp.V3 * Wtemp.V3);
#ifdef RADIATION_MHD
		pG->U[k][j][i].E += 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c);
#endif

		}

		}

	}/* end i */
	}/* end j */
	}/* end k */
	if(fabs(zmin) < 4.0 && (negd >0 || negP > 0)) printf("[Bad Cell]: Maximum height: %f\n",zmin);
	  if (negd > 0 || negP > 0)
    printf("[Bad Cells]: at t=%g; %i cells had d<0; %i cells had P<0\n",pG->time,negd,negP);

#endif /* FIRST_ORDER_FLUX_CORRECTION */



	/* Update the opacity if Opacity function is set in the problem generator */
/* Opacity is updated in the BackEuler step so that they are consistent */
	if(Opacity != NULL){
		for (k=ks; k<=ke; k++){
			for (j=js; j<=je; j++) {
    				for (i=is; i<=ie; i++){
				
				Wopacity = Cons_to_Prim(&pG->U[k][j][i]);
				
				/* Add background shearing */
#ifdef FARGO	
				cc_pos(pG,i,j,k,&x1,&x2,&x3);
				Wopacity.V2 -= qom * x1;		
#endif

				if(Wopacity.P > TINY_NUMBER)
				{
					Opacity(&Wopacity,Sigma,NULL);

					for(m=0;m<NOPACITY;m++){
						pG->U[k][j][i].Sigma[m] = Sigma[m];
					}

				}
				else{
					
					pG->Tguess[k][j][i] = pG->U[k][j][i].Er;
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


#ifdef FIRST_ORDER_FLUX_CORRECTION
  if ((x1FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3FluxP =(Cons1DS***)calloc_3d_array(size3,size2,size1, sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((INIdata =(Real****)calloc_4d_array(size3,size2,size1, 6, sizeof(Real)))
    == NULL) goto on_error;
  if ((dhalfP = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#if defined(MHD) || defined(RADIATION_MHD)
  if ((emf1P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf2P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((emf3P = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */



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
  if ((Flxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((Flxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
#endif


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


#ifdef FIRST_ORDER_FLUX_CORRECTION
  if (x1FluxP != NULL) free_3d_array(x1FluxP);
  if (x2FluxP != NULL) free_3d_array(x2FluxP);
  if (x3FluxP != NULL) free_3d_array(x3FluxP);
  if (INIdata != NULL) free_4d_array(INIdata);
  if (dhalfP != NULL) free_4d_array(dhalfP);
#if defined(MHD) || defined(RADIATION_MHD)
  if (emf1P   != NULL) free_3d_array(emf1P);
  if (emf2P   != NULL) free_3d_array(emf2P);
  if (emf3P   != NULL) free_3d_array(emf3P);
#endif
#endif /* FIRST_ORDER_FLUX_CORRECTION */





  if (dhalf     != NULL) free_3d_array(dhalf);
  if (phalf     != NULL) free_3d_array(phalf);

  if(Source	!= NULL) free_4d_array(Source);

  if(Alpha	!= NULL) free_3d_array(Alpha);

  if(Beta	!= NULL) free_4d_array(Beta);

#ifdef SHEARING_BOX
  if (Flxiib != NULL) free_2d_array(Flxiib);
  if (Flxoib != NULL) free_2d_array(Flxoib);
  if (rFlxiib != NULL) free_2d_array(rFlxiib);
  if (rFlxoib != NULL) free_2d_array(rFlxoib);
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




#ifdef FIRST_ORDER_FLUX_CORRECTION

static void integrate_emf1_corner_FOFC(const GridS *pG)
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
	if (x2FluxP[k-1][j][i].d > 0.0)
	  de1_l3 = x3FluxP[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2FluxP[k-1][j][i].d < 0.0)
	  de1_l3 = x3FluxP[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3FluxP[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3FluxP[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2FluxP[k][j][i].d > 0.0)
	  de1_r3 = x3FluxP[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2FluxP[k][j][i].d < 0.0)
	  de1_r3 = x3FluxP[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3FluxP[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3FluxP[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3FluxP[k][j-1][i].d > 0.0)
	  de1_l2 = -x2FluxP[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3FluxP[k][j-1][i].d < 0.0)
	  de1_l2 = -x2FluxP[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2FluxP[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2FluxP[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3FluxP[k][j][i].d > 0.0)
	  de1_r2 = -x2FluxP[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3FluxP[k][j][i].d < 0.0)
	  de1_r2 = -x2FluxP[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2FluxP[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2FluxP[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1P[k][j][i] = 0.25*(  x3FluxP[k][j][i].Bz + x3FluxP[k][j-1][i].Bz
                              - x2FluxP[k][j][i].By - x2FluxP[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

static void integrate_emf2_corner_FOFC(const GridS *pG)
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
	if (x1FluxP[k-1][j][i].d > 0.0)
	  de2_l3 = -x3FluxP[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1FluxP[k-1][j][i].d < 0.0)
	  de2_l3 = -x3FluxP[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3FluxP[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3FluxP[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1FluxP[k][j][i].d > 0.0)
	  de2_r3 = -x3FluxP[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1FluxP[k][j][i].d < 0.0)
	  de2_r3 = -x3FluxP[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3FluxP[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3FluxP[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3FluxP[k][j][i-1].d > 0.0)
	  de2_l1 = x1FluxP[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3FluxP[k][j][i-1].d < 0.0)
	  de2_l1 = x1FluxP[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1FluxP[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1FluxP[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3FluxP[k][j][i].d > 0.0)
	  de2_r1 = x1FluxP[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3FluxP[k][j][i].d < 0.0)
	  de2_r1 = x1FluxP[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1FluxP[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1FluxP[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2P[k][j][i] = 0.25*(  x1FluxP[k][j][i].Bz + x1FluxP[k-1][j][i  ].Bz
                              - x3FluxP[k][j][i].By - x3FluxP[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

static void integrate_emf3_corner_FOFC(const GridS *pG)
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
	if (x1FluxP[k][j-1][i].d > 0.0)
	  de3_l2 = (x2FluxP[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1FluxP[k][j-1][i].d < 0.0)
	  de3_l2 = (x2FluxP[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2FluxP[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf + 
			(x2FluxP[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1FluxP[k][j][i].d > 0.0)
	  de3_r2 = (x2FluxP[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1FluxP[k][j][i].d < 0.0)
	  de3_r2 = (x2FluxP[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2FluxP[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf + 
			(x2FluxP[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
	}

	if (x2FluxP[k][j][i-1].d > 0.0)
	  de3_l1 = -x1FluxP[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2FluxP[k][j][i-1].d < 0.0)
	  de3_l1 = -x1FluxP[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1FluxP[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1FluxP[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2FluxP[k][j][i].d > 0.0)
	  de3_r1 = -x1FluxP[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2FluxP[k][j][i].d < 0.0)
	  de3_r1 = -x1FluxP[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1FluxP[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1FluxP[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3P[k][j][i] = 0.25*(  x2FluxP[k][j  ][i-1].Bz + x2FluxP[k][j][i].Bz
			      - x1FluxP[k][j-1][i  ].By - x1FluxP[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}
#endif


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
	Real SPP, diffTEr;
/*	Real dSigma[8];
*/
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
				/* The momentum source term v * diffTEr/ Crat is added explicitly */
				Source[k][j][i][0] = 0.0;
				Source[k][j][i][1] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0x); /* + velocity_x * diffTEr / Crat); */
				Source[k][j][i][2] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0y); /* + velocity_y * diffTEr / Crat); */
				Source[k][j][i][3] = -Prat * (-(Sigma_aF + Sigma_sF) * Fr0z); /* + velocity_z * diffTEr / Crat); */
			
			/* Source term for energy */
				Source[k][j][i][4] = -Prat * Crat * (diffTEr + (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x + velocity_y * Fr0y 
												+ velocity_z * Fr0z)/Crat);

/*
				if(Opacity != NULL){
					 Opacity(density, temperature, NULL, dSigma);
				}
				else{
					for(m=0;m<2*NOPACITY;m++)
						dSigma[m] = 0.0;
				}
				
*/
				/* dSigma[0] = dSigma_sF/drho, dSigma[1] = dSigma_aF/drho, dSigma[2]=dSigma_aP/drho, dSigma[3]= dSigma_aE/drho */
				/* dSigma[4] = dSigma_sF/dT, dSigma[5] = dSigma_aF/dT, dSigma[6]=dSigma_aP/dT, dSigma[7]= dSigma_aE/dT */
						

/*				dSigmadP[0] =  dSigma[4] / (density * R_ideal); 
				dSigmadP[1] =  dSigma[5] / (density * R_ideal); 
				dSigmadP[2] =  dSigma[6] / (density * R_ideal); 
				dSigmadP[3] =  dSigma[7] / (density * R_ideal); 

				diffTErdP = dSigmadP[2] * pow(temperature, 4.0) - dSigmadP[3] * pG->U[k][j][i].Er;
*/
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
				SVVx = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_11) * Usource.Er) / (density * Crat);
			/* + diffTEr) / (density * Crat);*/
	
				
				if(fabs(SVVx * dt * 0.5) > 50.0)
					betax = -1.0/(SVVx * dt * 0.5);	
				else if(fabs(SVVx * dt * 0.5) > 0.001)
					betax = (exp(SVVx * dt * 0.5) - 1.0)/(SVVx * dt * 0.5);
				else 
					betax = 1.0 + 0.25 * SVVx * dt;

				SVVy = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_22) * Usource.Er) / (density * Crat);
				/* + diffTEr) / (density * Crat); */
		
				
				if(fabs(SVVy * dt * 0.5) > 50.0)
					betay = -1.0/(SVVy * dt * 0.5);
				else if(fabs(SVVy * dt * 0.5) > 0.001)
					betay = (exp(SVVy * dt * 0.5) - 1.0)/(SVVy * dt * 0.5);
				else 
					betay = 1.0 + 0.25 * SVVy * dt;

				SVVz = -Prat * ((Sigma_aF + Sigma_sF) * (1.0 + Usource.Edd_33) * Usource.Er)/ (density * Crat);
				/* + diffTEr) / (density * Crat);*/
		
				
				if(fabs(SVVz * dt * 0.5) > 50.0)
					betaz = -1.0/(SVVz * dt * 0.5);
				else if(fabs(SVVz * dt * 0.5) > 0.001)
					betaz = (exp(SVVz * dt * 0.5) - 1.0)/(SVVz * dt * 0.5);
				else 	
					betaz = 1.0 + 0.25 * SVVz * dt;
		/*===========================================================================*/
	
				if(fabs(SPP * dt * 0.5) > 50.0)
					alpha = -1.0/(SPP * dt * 0.5);
				else if(fabs(SPP * dt * 0.5) > 0.001)
					alpha = (exp(SPP * dt * 0.5) - 1.0)/(SPP * dt * 0.5);
				else 
					alpha = 1.0 + 0.25 * SPP * dt;

				
				Alpha[k][j][i] = alpha;
				Beta[k][j][i][0] = betax;
				Beta[k][j][i][1] = betay;
				Beta[k][j][i][2] = betaz;
			
				/* In gas pressure dominated case, do something special */
				if((4.0 * (Gamma - 1.0) * Prat * temperature * temperature * temperature / (density * R_ideal)) < 1.0){
					/* Source term for energy */
					Source[k][j][i][4] = -Prat * (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z)
					-Prat * pG->Ersource[k][j][i] / (dt * alpha);
				}		

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




/* During this correction step, do not include Radiation source term */
#ifdef FIRST_ORDER_FLUX_CORRECTION
/*----------------------------------------------------------------------------*/
/*! \fn static void FixCell(GridS *pG, Int3Vect ix)
 *  \brief Uses first order fluxes to fix negative d,P or superluminal v
 */ 
static void FixCell(GridS *pG, Int3Vect ix)
{
  int n;
#if defined(MHD) || defined(RADIATION_MHD)
  int i,j,k;
#endif
  Real rsf=1.0,lsf=1.0;


  x1FD_i.d = x1Flux[ix.k][ix.j][ix.i].d - x1FluxP[ix.k][ix.j][ix.i].d;
  x2FD_j.d = x2Flux[ix.k][ix.j][ix.i].d - x2FluxP[ix.k][ix.j][ix.i].d;
  x3FD_k.d = x3Flux[ix.k][ix.j][ix.i].d - x3FluxP[ix.k][ix.j][ix.i].d;

  x1FD_ip1.d = x1Flux[ix.k][ix.j][ix.i+1].d - x1FluxP[ix.k][ix.j][ix.i+1].d;
  x2FD_jp1.d = x2Flux[ix.k][ix.j+1][ix.i].d - x2FluxP[ix.k][ix.j+1][ix.i].d;
  x3FD_kp1.d = x3Flux[ix.k+1][ix.j][ix.i].d - x3FluxP[ix.k+1][ix.j][ix.i].d;

  x1FD_i.Mx = x1Flux[ix.k][ix.j][ix.i].Mx - x1FluxP[ix.k][ix.j][ix.i].Mx;
  x2FD_j.Mx = x2Flux[ix.k][ix.j][ix.i].Mx - x2FluxP[ix.k][ix.j][ix.i].Mx;
  x3FD_k.Mx = x3Flux[ix.k][ix.j][ix.i].Mx - x3FluxP[ix.k][ix.j][ix.i].Mx;

  x1FD_ip1.Mx = x1Flux[ix.k][ix.j][ix.i+1].Mx - x1FluxP[ix.k][ix.j][ix.i+1].Mx;
  x2FD_jp1.Mx = x2Flux[ix.k][ix.j+1][ix.i].Mx - x2FluxP[ix.k][ix.j+1][ix.i].Mx;
  x3FD_kp1.Mx = x3Flux[ix.k+1][ix.j][ix.i].Mx - x3FluxP[ix.k+1][ix.j][ix.i].Mx;

  x1FD_i.My = x1Flux[ix.k][ix.j][ix.i].My - x1FluxP[ix.k][ix.j][ix.i].My;
  x2FD_j.My = x2Flux[ix.k][ix.j][ix.i].My - x2FluxP[ix.k][ix.j][ix.i].My;
  x3FD_k.My = x3Flux[ix.k][ix.j][ix.i].My - x3FluxP[ix.k][ix.j][ix.i].My;

  x1FD_ip1.My = x1Flux[ix.k][ix.j][ix.i+1].My - x1FluxP[ix.k][ix.j][ix.i+1].My;
  x2FD_jp1.My = x2Flux[ix.k][ix.j+1][ix.i].My - x2FluxP[ix.k][ix.j+1][ix.i].My;
  x3FD_kp1.My = x3Flux[ix.k+1][ix.j][ix.i].My - x3FluxP[ix.k+1][ix.j][ix.i].My;

  x1FD_i.Mz = x1Flux[ix.k][ix.j][ix.i].Mz - x1FluxP[ix.k][ix.j][ix.i].Mz;
  x2FD_j.Mz = x2Flux[ix.k][ix.j][ix.i].Mz - x2FluxP[ix.k][ix.j][ix.i].Mz;
  x3FD_k.Mz = x3Flux[ix.k][ix.j][ix.i].Mz - x3FluxP[ix.k][ix.j][ix.i].Mz;

  x1FD_ip1.Mz = x1Flux[ix.k][ix.j][ix.i+1].Mz - x1FluxP[ix.k][ix.j][ix.i+1].Mz;
  x2FD_jp1.Mz = x2Flux[ix.k][ix.j+1][ix.i].Mz - x2FluxP[ix.k][ix.j+1][ix.i].Mz;
  x3FD_kp1.Mz = x3Flux[ix.k+1][ix.j][ix.i].Mz - x3FluxP[ix.k+1][ix.j][ix.i].Mz;


  x1FD_i.E = x1Flux[ix.k][ix.j][ix.i].E - x1FluxP[ix.k][ix.j][ix.i].E;
  x2FD_j.E = x2Flux[ix.k][ix.j][ix.i].E - x2FluxP[ix.k][ix.j][ix.i].E;
  x3FD_k.E = x3Flux[ix.k][ix.j][ix.i].E - x3FluxP[ix.k][ix.j][ix.i].E;

  x1FD_ip1.E = x1Flux[ix.k][ix.j][ix.i+1].E - x1FluxP[ix.k][ix.j][ix.i+1].E;
  x2FD_jp1.E = x2Flux[ix.k][ix.j+1][ix.i].E - x2FluxP[ix.k][ix.j+1][ix.i].E;
  x3FD_kp1.E = x3Flux[ix.k+1][ix.j][ix.i].E - x3FluxP[ix.k+1][ix.j][ix.i].E;




/* We just re-do the flux divergence for hydro part */
/* For MHD part, we add the difference of flux divergence */

#if defined(MHD) || defined(RADIATION_MHD)
  emf1D_kj     = emf1[ix.k  ][ix.j  ][ix.i] - emf1P[ix.k  ][ix.j  ][ix.i];
  emf1D_kjp1   = emf1[ix.k  ][ix.j+1][ix.i] - emf1P[ix.k  ][ix.j+1][ix.i];
  emf1D_kp1j   = emf1[ix.k+1][ix.j  ][ix.i] - emf1P[ix.k+1][ix.j  ][ix.i];
  emf1D_kp1jp1 = emf1[ix.k+1][ix.j+1][ix.i] - emf1P[ix.k+1][ix.j+1][ix.i];

  emf2D_ki     = emf2[ix.k  ][ix.j][ix.i  ] - emf2P[ix.k  ][ix.j][ix.i  ];
  emf2D_kip1   = emf2[ix.k  ][ix.j][ix.i+1] - emf2P[ix.k  ][ix.j][ix.i+1];
  emf2D_kp1i   = emf2[ix.k+1][ix.j][ix.i  ] - emf2P[ix.k+1][ix.j][ix.i  ];
  emf2D_kp1ip1 = emf2[ix.k+1][ix.j][ix.i+1] - emf2P[ix.k+1][ix.j][ix.i+1];

  emf3D_ji     = emf3[ix.k][ix.j  ][ix.i  ] - emf3P[ix.k][ix.j  ][ix.i  ];
  emf3D_jip1   = emf3[ix.k][ix.j  ][ix.i+1] - emf3P[ix.k][ix.j  ][ix.i+1];
  emf3D_jp1i   = emf3[ix.k][ix.j+1][ix.i  ] - emf3P[ix.k][ix.j+1][ix.i  ];
  emf3D_jp1ip1 = emf3[ix.k][ix.j+1][ix.i+1] - emf3P[ix.k][ix.j+1][ix.i+1];
#endif /* MHD */

/* Use flux differences to correct bad cell-centered quantities */
  ApplyCorr(pG,ix.i,ix.j,ix.k,1,1,1,1,1,1);

#ifdef SELF_GRAVITY
/* Save mass fluxes in Grid structure for source term correction in main loop */
  pG->x1MassFlux[ix.k  ][ix.j  ][ix.i  ] = x1FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x2MassFlux[ix.k  ][ix.j  ][ix.i  ] = x2FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x3MassFlux[ix.k  ][ix.j  ][ix.i  ] = x3FluxP[ix.k  ][ix.j  ][ix.i  ].d;
  pG->x1MassFlux[ix.k  ][ix.j  ][ix.i+1] = x1FluxP[ix.k  ][ix.j  ][ix.i+1].d;
  pG->x2MassFlux[ix.k  ][ix.j+1][ix.i  ] = x2FluxP[ix.k  ][ix.j+1][ix.i  ].d;
  pG->x3MassFlux[ix.k+1][ix.j  ][ix.i  ] = x3FluxP[ix.k+1][ix.j  ][ix.i  ].d;
#endif /* SELF_GRAVITY */



	/* correct neighbor cells */
  if (ix.i > pG->is) ApplyCorr2(pG,ix.i-1,ix.j,ix.k,-1,0,0,0,0,0);
  if (ix.i < pG->ie) ApplyCorr2(pG,ix.i+1,ix.j,ix.k,0,-1,0,0,0,0);

  if (ix.j > pG->js) ApplyCorr2(pG,ix.i,ix.j-1,ix.k,0,0,-1,0,0,0);
  if (ix.j < pG->je) ApplyCorr2(pG,ix.i,ix.j+1,ix.k,0,0,0,-1,0,0);

  if (ix.k > pG->ks) ApplyCorr2(pG,ix.i,ix.j,ix.k-1,0,0,0,0,-1,0);
  if (ix.k < pG->ke) ApplyCorr2(pG,ix.i,ix.j,ix.k+1,0,0,0,0,0,-1);


#if defined(MHD) || defined(RADIATION_MHD)
  FixMHD(pG,ix.i,ix.j,ix.k,1,1,1,1,1,1);

/* Use flux differences to correct cell-centered values at i-1 and i+1 */
  
  if (ix.i > pG->is) FixMHD(pG,ix.i-1,ix.j,ix.k,-1,0,0,0,0,0);
  if (ix.i < pG->ie) FixMHD(pG,ix.i+1,ix.j,ix.k,0,-1,0,0,0,0);

/* Use flux differences to correct cell-centered values at j-1 and j+1 */
  
  if (ix.j > pG->js) FixMHD(pG,ix.i,ix.j-1,ix.k,0,0,-1,0,0,0);
  if (ix.j < pG->je) FixMHD(pG,ix.i,ix.j+1,ix.k,0,0,0,-1,0,0);

/* Use flux differences to correct cell-centered values at k-1 and k+1 */
 
  if (ix.k > pG->ks) FixMHD(pG,ix.i,ix.j,ix.k-1,0,0,0,0,-1,0);
  if (ix.k < pG->ke) FixMHD(pG,ix.i,ix.j,ix.k+1,0,0,0,0,0,-1);

/* Compute new cell-centered fields */
  for (k=(ix.k-1); k<=(ix.k+1); k++) {
  for (j=(ix.j-1); j<=(ix.j+1); j++) {
  for (i=(ix.i-1); i<=(ix.i+1); i++) {
#ifdef CYLINDRICAL
    rsf = pG->ri[i+1]/pG->r[i];  lsf = pG->ri[i]/pG->r[i];
#endif
    pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
    pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
    pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
  }}}
#endif /* MHD */


/* Must replace higher-order fluxes/emfs with predict fluxes/emfs in case two
 * adjacent cells are corrected.  Otherwise, flux/emf differences may get
 * corrected more than once.  By overwriting the higher-order fluxes/emfs,
 * the second time through, the differences will be zero. */
  x1Flux[ix.k  ][ix.j  ][ix.i  ] = x1FluxP[ix.k  ][ix.j  ][ix.i  ];
  x2Flux[ix.k  ][ix.j  ][ix.i  ] = x2FluxP[ix.k  ][ix.j  ][ix.i  ];
  x3Flux[ix.k  ][ix.j  ][ix.i  ] = x3FluxP[ix.k  ][ix.j  ][ix.i  ];
  x1Flux[ix.k  ][ix.j  ][ix.i+1] = x1FluxP[ix.k  ][ix.j  ][ix.i+1];
  x2Flux[ix.k  ][ix.j+1][ix.i  ] = x2FluxP[ix.k  ][ix.j+1][ix.i  ];
  x3Flux[ix.k+1][ix.j  ][ix.i  ] = x3FluxP[ix.k+1][ix.j  ][ix.i  ];

#if defined(MHD) || defined(RADIATION_MHD)
  emf1[ix.k  ][ix.j  ][ix.i  ] = emf1P[ix.k  ][ix.j  ][ix.i  ];
  emf1[ix.k  ][ix.j+1][ix.i  ] = emf1P[ix.k  ][ix.j+1][ix.i  ];
  emf1[ix.k+1][ix.j  ][ix.i  ] = emf1P[ix.k+1][ix.j  ][ix.i  ];
  emf1[ix.k+1][ix.j+1][ix.i  ] = emf1P[ix.k+1][ix.j+1][ix.i  ];
  emf2[ix.k  ][ix.j  ][ix.i  ] = emf2P[ix.k  ][ix.j  ][ix.i  ];
  emf2[ix.k  ][ix.j  ][ix.i+1] = emf2P[ix.k  ][ix.j  ][ix.i+1];
  emf2[ix.k+1][ix.j  ][ix.i  ] = emf2P[ix.k+1][ix.j  ][ix.i  ];
  emf2[ix.k+1][ix.j  ][ix.i+1] = emf2P[ix.k+1][ix.j  ][ix.i+1];
  emf3[ix.k  ][ix.j  ][ix.i  ] = emf3P[ix.k  ][ix.j  ][ix.i  ];
  emf3[ix.k  ][ix.j  ][ix.i+1] = emf3P[ix.k  ][ix.j  ][ix.i+1];
  emf3[ix.k  ][ix.j+1][ix.i  ] = emf3P[ix.k  ][ix.j+1][ix.i  ];
  emf3[ix.k  ][ix.j+1][ix.i+1] = emf3P[ix.k  ][ix.j+1][ix.i+1];
#endif
}



/* Radiation source terms are not added */
static void ApplyCorr(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3)
{
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real x1,x2,x3,phil,phir,phic;
  int n;
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  Real M1e, dM2e, M1n, dM2n; /* M1, dM2 evolved by dt/2 */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
  Real fact, qom, om_dt = Omega_0*pG->dt;
  Real ShearSource[4];
#endif /* SHEARING_BOX */
  Real lsf=1.0,rsf=1.0;
  Real velocity_x, velocity_y, velocity_z, diffTEr, Fr0x, Fr0y, Fr0z, Sigma_aP, Sigma_aE, Sigma_aF, Sigma_sF, pressure, Temperature, Prwork, TEr, density;
  Real SFmx, SFmy, SFmz;
  Real Smx0, Smy0, Smz0;
  Real coef1, coef2, coef3, Ersum;

#ifdef CYLINDRICAL
  rsf = (rx1>0) ? pG->ri[i+1]/pG->r[i] : pG->ri[i  ]/pG->r[i];
  lsf = (lx1>0) ? pG->ri[i  ]/pG->r[i] : pG->ri[i+1]/pG->r[i];
  dtodx2 = pG->dt/(pG->r[i]*pG->dx2);
#endif


 pG->U[k][j][i].d  = INIdata[k][j][i][0] - (dtodx1 * (x1FluxP[k][j][i+1].d  - x1FluxP[k][j][i].d)
					 +  dtodx2 * (x2FluxP[k][j+1][i].d  - x2FluxP[k][j][i].d)
					 +  dtodx3 * (x3FluxP[k+1][j][i].d  - x3FluxP[k][j][i].d));


 pG->U[k][j][i].M1 = INIdata[k][j][i][1] - (dtodx1 * (x1FluxP[k][j][i+1].Mx  - x1FluxP[k][j][i].Mx)
					 +  dtodx2 * (x2FluxP[k][j+1][i].Mz  - x2FluxP[k][j][i].Mz)
					 +  dtodx3 * (x3FluxP[k+1][j][i].My  - x3FluxP[k][j][i].My));


 pG->U[k][j][i].M2 = INIdata[k][j][i][2] - (dtodx1 * (x1FluxP[k][j][i+1].My  - x1FluxP[k][j][i].My)
					 +  dtodx2 * (x2FluxP[k][j+1][i].Mx  - x2FluxP[k][j][i].Mx)
					 +  dtodx3 * (x3FluxP[k+1][j][i].Mz  - x3FluxP[k][j][i].Mz));



 pG->U[k][j][i].M3 = INIdata[k][j][i][3] - (dtodx1 * (x1FluxP[k][j][i+1].Mz  - x1FluxP[k][j][i].Mz)
					 +  dtodx2 * (x2FluxP[k][j+1][i].My  - x2FluxP[k][j][i].My)
					 +  dtodx3 * (x3FluxP[k+1][j][i].Mx  - x3FluxP[k][j][i].Mx));



 pG->U[k][j][i].E  = INIdata[k][j][i][4] - (dtodx1 * (x1FluxP[k][j][i+1].E  - x1FluxP[k][j][i].E)
					 +  dtodx2 * (x2FluxP[k][j+1][i].E  - x2FluxP[k][j][i].E)
					 +  dtodx3 * (x3FluxP[k+1][j][i].E  - x3FluxP[k][j][i].E));


#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  phic = (*ShearingBoxPot)((x1            ),x2,x3);
		phir = (*ShearingBoxPot)((x1+0.5*pG->dx1),x2,x3);
		phil = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

		ShearSource[3] = dtodx1*(x1FluxP[k][j][i  ].d*(phic - phil) +
				x1FluxP[k][j][i+1].d*(phir - phic));
					
		phir = (*ShearingBoxPot)(x1,(x2+0.5*pG->dx2),x3);
		phil = (*ShearingBoxPot)(x1,(x2-0.5*pG->dx2),x3);

		ShearSource[3] -= dtodx2*(x2FluxP[k][j  ][i].d*(phic - phil) +
				x2FluxP[k][j+1][i].d*(phir - phic));

							
		phir = (*ShearingBoxPot)(x1,x2,(x3+0.5*pG->dx3));
		phil = (*ShearingBoxPot)(x1,x2,(x3-0.5*pG->dx3));

		ShearSource[3] -= dtodx3*(x3FluxP[k  ][j][i].d*(phic - phil) +
				x3FluxP[k+1][j][i].d*(phir - phic));

		/* Use Crank Nichson but do not include radiation source term for guess solution */
		/* Store the current state */
		M1n  = INIdata[k][j][i][1];
#ifdef FARGO
		dM2n = INIdata[k][j][i][2];
#else
		dM2n = INIdata[k][j][i][2] + qom*x1*INIdata[k][j][i][0];
#endif
					

#ifndef FARGO
		frx1_dM2 = qom*(x1+0.5*pG->dx1)*x1FluxP[k][j][i+1].d;
		flx1_dM2 = qom*(x1-0.5*pG->dx1)*x1FluxP[k][j][i  ].d;
		frx2_dM2 = qom*(x1            )*x2FluxP[k][j+1][i].d;
		flx2_dM2 = qom*(x1            )*x2FluxP[k][j  ][i].d;
		frx3_dM2 = qom*(x1            )*x3FluxP[k+1][j][i].d;
		flx3_dM2 = qom*(x1            )*x3FluxP[k  ][j][i].d;
#endif
/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
		M1e  = M1n + 0.5 * (pG->U[k][j][i].M1 - INIdata[k][j][i][1]);
		dM2e = dM2n + 0.5 * (pG->U[k][j][i].M2 - INIdata[k][j][i][2]);
#ifndef FARGO
		dM2e -= hdtodx1*(frx1_dM2 - flx1_dM2)
			- hdtodx2*(frx2_dM2 - flx2_dM2) 
			- hdtodx3*(frx3_dM2 - flx3_dM2);
#endif

		ShearSource[0] = (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact; 
		ShearSource[1] = 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;
#ifndef FARGO
		ShearSource[1] -= 0.5*qshear*om_dt*(x1FluxP[k][j][i].d + x1FluxP[k][j][i+1].d);
#endif		
		/* Now add the source term to guess solution */


		pG->U[k][j][i].M1 += ShearSource[0];
		pG->U[k][j][i].M2 += ShearSource[1];
		pG->U[k][j][i].M3 += ShearSource[2];
		pG->U[k][j][i].E += ShearSource[3];


#endif /* SHEARING_BOX */

  if (StaticGravPot != NULL){
    	cc_pos(pG,i,j,k,&x1,&x2,&x3);
   	phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

	/* x direction */

	pG->U[k][j][i].M1 -= dtodx1 * dhalfP[k][j][i]*(phir-phil);
        pG->U[k][j][i].E  -= dtodx1 * (x1FluxP[k][j][i  ].d*(phic - phil) +
		x1FluxP[k][j][i+1].d*(phir - phic));			

	/* y direction */
			
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        pG->U[k][j][i].M2 -= dtodx2 * dhalfP[k][j][i]*(phir-phil);
        pG->U[k][j][i].E  -= dtodx2 * (x2FluxP[k][j  ][i].d*(phic - phil) +
                               	     x2FluxP[k][j+1][i].d*(phir - phic));
	/* z direction */
	phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

	pG->U[k][j][i].M3 -= dtodx3 * dhalfP[k][j][i]*(phir-phil);
        pG->U[k][j][i].E  -= dtodx3 * (x3FluxP[k][j  ][i].d*(phic - phil) +
        	     x3FluxP[k + 1][j][i].d*(phir - phic));	
  }
/* Add radiation momentum source term seperately as first orer */

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	/* Part of momentum source term */
	/* Use unpdated gas quantities */
	density = INIdata[k][j][i][0];

	velocity_x = INIdata[k][j][i][1] / density;
	velocity_y = INIdata[k][j][i][2] / density;
#ifdef FARGO
				
		/* With FARGO, we should add background shearing to the source terms */
	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	velocity_y -= qom * x1;			
				
#endif	
	velocity_z = INIdata[k][j][i][3] / density;

	pressure = INIdata[k][j][i][5];
	Temperature = pressure / (density * R_ideal);

	Sigma_sF = pG->U[k][j][i].Sigma[0];
	Sigma_aF = pG->U[k][j][i].Sigma[1];
	Sigma_aP = pG->U[k][j][i].Sigma[2];
	Sigma_aE = pG->U[k][j][i].Sigma[3];

	diffTEr = Sigma_aP * pow(Temperature, 4.0) - Sigma_aE * pG->U[k][j][i].Er; 

	SFmx = (Sigma_aF + Sigma_sF) * (1.0 + pG->U[k][j][i].Edd_11) * pG->U[k][j][i].Er / (density * Crat); 
	SFmy = (Sigma_aF + Sigma_sF) * (1.0 + pG->U[k][j][i].Edd_22) * pG->U[k][j][i].Er / (density * Crat); 
	SFmz = (Sigma_aF + Sigma_sF) * (1.0 + pG->U[k][j][i].Edd_33) * pG->U[k][j][i].Er / (density * Crat); 

	Smx0 = -Prat * velocity_x * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * pG->U[k][j][i].Er) / Crat;
	Smy0 = -Prat * velocity_y * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * pG->U[k][j][i].Er) / Crat;
	Smz0 = -Prat * velocity_z * (Sigma_aP * pG->Tguess[k][j][i] - Sigma_aE * pG->U[k][j][i].Er) / Crat;
	
	/* Part of radiation momentum source term is not included in Source */
	/* We only need first order, we do not do the correct step */
	pG->U[k][j][i].M1 += pG->dt * Source[k][j][i][1] / (1.0 + pG->dt * Prat * SFmx) + 0.5 * pG->dt * Smx0;
	pG->U[k][j][i].M2 += pG->dt * Source[k][j][i][2] / (1.0 + pG->dt * Prat * SFmy) + 0.5 * pG->dt * Smy0;
	pG->U[k][j][i].M3 += pG->dt * Source[k][j][i][3] / (1.0 + pG->dt * Prat * SFmz) + 0.5 * pG->dt * Smz0;
	
	/* First add radiation work term explicitly */
	Fr0x = pG->U[k][j][i].Fr1 - ((1.0 + pG->U[k][j][i].Edd_11) * velocity_x + pG->U[k][j][i].Edd_21 * velocity_y
							+ pG->U[k][j][i].Edd_31 * velocity_z) * pG->U[k][j][i].Er/Crat;
	Fr0y = pG->U[k][j][i].Fr2 - ((1.0 + pG->U[k][j][i].Edd_22) * velocity_y + pG->U[k][j][i].Edd_21 * velocity_x
							+ pG->U[k][j][i].Edd_32 * velocity_z) * pG->U[k][j][i].Er/Crat;
	Fr0z = pG->U[k][j][i].Fr3 - ((1.0 + pG->U[k][j][i].Edd_33) * velocity_z + pG->U[k][j][i].Edd_31 * velocity_x
							+ pG->U[k][j][i].Edd_32 * velocity_y) * pG->U[k][j][i].Er/Crat;

	Prwork = -Prat * (Sigma_aF - Sigma_sF) * (velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z);
	
	/* First, add Prwork term explicitly */
	pG->U[k][j][i].E  += pG->dt * Prwork;

	/* Then add Sigma(T4 - Er) implicitly */
	TEr = pow(pG->U[k][j][i].Er, 0.25);
	Ersum = pressure / (Gamma - 1.0) + Prat * pG->U[k][j][i].Er;
	coef1 = pG->dt * Prat * Crat * Sigma_aP;
	coef2 = density * R_ideal * (1.0 + pG->dt * Sigma_aE * Crat) / (Gamma - 1.0);
	coef3 = -pressure / (Gamma - 1.0) - pG->dt * Sigma_aE * Crat * Ersum;
	if(fabs(TEr - Temperature) < TINY_NUMBER){
		pG->Ersource[k][j][i] = 0.0;
		pG->Tguess[k][j][i] = pG->U[k][j][i].Er;
	}
	else{
		if(TEr < Temperature){
			Temperature = rtsafe(Tequilibrium, TEr * (1.0 - 0.01), Temperature * (1.0 + 0.01), 1.e-14, coef1, coef2, coef3, 0.0);
		}
		else{
		
			Temperature = rtsafe(Tequilibrium, Temperature * (1.0 - 0.01), TEr * (1.0 + 0.01), 1.e-14, coef1, coef2, coef3, 0.0);
		}
		
		 pG->U[k][j][i].E += (density * R_ideal * Temperature - pressure) / (Gamma - 1.0);

		 pG->Ersource[k][j][i] = -(density * R_ideal * Temperature - pressure) / (Prat * (Gamma - 1.0));

		
		pG->Tguess[k][j][i] = pow(Temperature, 4.0);
		
	
	}
#endif

}


/* Use to correct flux of neighbor cells */
static void ApplyCorr2(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3)
{
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real x1,x2,x3,phil,phir,phic;
  int n;
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */
  Real lsf=1.0,rsf=1.0;

#ifdef CYLINDRICAL
  rsf = (rx1>0) ? pG->ri[i+1]/pG->r[i] : pG->ri[i  ]/pG->r[i];
  lsf = (lx1>0) ? pG->ri[i  ]/pG->r[i] : pG->ri[i+1]/pG->r[i];
  dtodx2 = pG->dt/(pG->r[i]*pG->dx2);
#endif
  pG->U[k][j][i].d  += dtodx1*(rsf*rx1*x1FD_ip1.d  - lsf*lx1*x1FD_i.d );
  pG->U[k][j][i].M1 += dtodx1*(rsf*rx1*x1FD_ip1.Mx - lsf*lx1*x1FD_i.Mx);
  pG->U[k][j][i].M2 += dtodx1*(SQR(rsf)*rx1*x1FD_ip1.My - SQR(lsf)*lx1*x1FD_i.My);
  pG->U[k][j][i].M3 += dtodx1*(rsf*rx1*x1FD_ip1.Mz - lsf*lx1*x1FD_i.Mz);

  pG->U[k][j][i].E  += dtodx1*(rsf*rx1*x1FD_ip1.E  - lsf*lx1*x1FD_i.E );

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++)
    pG->U[k][j][i].s[n] += dtodx1*(rsf*rx1*x1FD_ip1.s[n] - lsf*lx1*x1FD_i.s[n]);
#endif

  pG->U[k][j][i].d  += dtodx2*(rx2*x2FD_jp1.d  - lx2*x2FD_j.d );
  pG->U[k][j][i].M1 += dtodx2*(rx2*x2FD_jp1.Mz - lx2*x2FD_j.Mz);
  pG->U[k][j][i].M2 += dtodx2*(rx2*x2FD_jp1.Mx - lx2*x2FD_j.Mx);
  pG->U[k][j][i].M3 += dtodx2*(rx2*x2FD_jp1.My - lx2*x2FD_j.My);

  pG->U[k][j][i].E  += dtodx2*(rx2*x2FD_jp1.E  - lx2*x2FD_j.E );

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++)
    pG->U[k][j][i].s[n] += dtodx2*(rx2*x2FD_jp1.s[n] - lx2*x2FD_j.s[n]);
#endif

  pG->U[k][j][i].d  += dtodx3*(rx3*x3FD_kp1.d  - lx3*x3FD_k.d );
  pG->U[k][j][i].M1 += dtodx3*(rx3*x3FD_kp1.My - lx3*x3FD_k.My);
  pG->U[k][j][i].M2 += dtodx3*(rx3*x3FD_kp1.Mz - lx3*x3FD_k.Mz);
  pG->U[k][j][i].M3 += dtodx3*(rx3*x3FD_kp1.Mx - lx3*x3FD_k.Mx);

  pG->U[k][j][i].E  += dtodx3*(rx3*x3FD_kp1.E  - lx3*x3FD_k.E );



#ifdef SHEARING_BOX
  fact = om_dt/(2. + (2.-qshear)*om_dt*om_dt);
  qom = qshear*Omega_0;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

/* Calculate the flux for the y-momentum fluctuation */
  frx1_dM2 = x1FD_ip1.My;
  flx1_dM2 = x1FD_i.My;
  frx2_dM2 = x2FD_jp1.Mx;
  flx2_dM2 = x2FD_j.Mx;
  frx3_dM2 = x3FD_kp1.Mz;
  flx3_dM2 = x3FD_k.Mz;
#ifndef FARGO
  frx1_dM2 += qom*(x1+rx1*0.5*pG->dx1)*x1FD_ip1.d;
  flx1_dM2 += qom*(x1-lx1*0.5*pG->dx1)*x1FD_i.d;
  frx2_dM2 += qom*(x1            )*x2FD_jp1.d;
  flx2_dM2 += qom*(x1            )*x2FD_j.d;
  frx3_dM2 += qom*(x1            )*x3FD_kp1.d;
  flx3_dM2 += qom*(x1            )*x3FD_k.d;
#endif

/* Now evolve M1n and dM2n by dt/2 using Forward Euler */
  M1e = q1*(rx1*x1FD_ip1.Mx - lx1*x1FD_i.Mx)
      + q2*(rx2*x2FD_jp1.Mz - lx2*x2FD_j.Mz)
      + q3*(rx3*x3FD_kp1.My - lx3*x3FD_k.My);

  dM2e = q1*(rx1*frx1_dM2 - lx1*flx1_dM2)
       + q2*(rx2*frx2_dM2 - lx2*flx2_dM2) 
       + q3*(rx3*frx3_dM2 - lx3*flx3_dM2);

/* Update the 1- and 2-momentum for the Coriolis and tidal
 * potential momentum source terms using a Crank-Nicholson
 * discretization for the momentum fluctuation equation. */

  pG->U[k][j][i].M1 += (4.0*dM2e + 2.0*(qshear-2.)*om_dt*M1e)*fact;
  pG->U[k][j][i].M2 += 2.0*(qshear-2.)*(M1e + om_dt*dM2e)*fact;

#ifndef FARGO
  pG->U[k][j][i].M2 += 0.5*qshear*om_dt*(ABS(lx1)*x1FD_i.d + ABS(rx1)*x1FD_ip1.d);
#endif

/* Update the energy for a fixed potential.
 * This update is identical to non-SHEARING_BOX below  */

  phic = (*ShearingBoxPot)((x1            ),x2,x3);
  phir = (*ShearingBoxPot)((x1+rx1*0.5*pG->dx1),x2,x3);
  phil = (*ShearingBoxPot)((x1-lx1*0.5*pG->dx1),x2,x3);
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx1*(lx1*x1FD_i.d*(phic - phil) +
                              rx1*x1FD_ip1.d*(phir - phic));
#endif

  phir = (*ShearingBoxPot)(x1,(x2+rx2*0.5*pG->dx2),x3);
  phil = (*ShearingBoxPot)(x1,(x2-lx2*0.5*pG->dx2),x3);
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                              rx2*x2FD_jp1.d*(phir - phic));
#endif

  phir = (*ShearingBoxPot)(x1,x2,(x3+rx3*0.5*pG->dx3));
  phil = (*ShearingBoxPot)(x1,x2,(x3-lx3*0.5*pG->dx3));
#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                              rx3*x3FD_kp1.d*(phir - phic));
#endif
#endif /* SHEARING_BOX */

  if (StaticGravPot != NULL){
    cc_pos(pG,i,j,k,&x1,&x2,&x3);
    phic = (*StaticGravPot)(x1,x2,x3);
    phir = (*StaticGravPot)((x1+rx1*0.5*pG->dx1),x2,x3);
    phil = (*StaticGravPot)((x1-lx1*0.5*pG->dx1),x2,x3);

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx1*(lsf*lx1*x1FD_i.d*(phic - phil) +
                                rsf*rx1*x1FD_ip1.d*(phir - phic));
#endif

    phir = (*StaticGravPot)(x1,(x2+rx2*0.5*pG->dx2),x3);
    phil = (*StaticGravPot)(x1,(x2-lx2*0.5*pG->dx2),x3);

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                                rx2*x2FD_jp1.d*(phir - phic));
#endif

    phir = (*StaticGravPot)(x1,x2,(x3+rx3*0.5*pG->dx3));
    phil = (*StaticGravPot)(x1,x2,(x3-lx3*0.5*pG->dx3));

#ifndef BAROTROPIC
    pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                                rx3*x3FD_kp1.d*(phir - phic));
#endif
  }

#ifdef SELF_GRAVITY
  phic = pG->Phi[k][j][i];
  phil = 0.5*(pG->Phi[k][j][i-lx1] + pG->Phi[k][j][i  ]);
  phir = 0.5*(pG->Phi[k][j][i  ] + pG->Phi[k][j][i+rx1]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx1*(lx1*x1FD_i.d*(phic - phil) +
                              rx1*x1FD_ip1.d*(phir - phic));
#endif

  phil = 0.5*(pG->Phi[k][j-lx2][i] + pG->Phi[k][j  ][i]);
  phir = 0.5*(pG->Phi[k][j  ][i] + pG->Phi[k][j+rx2][i]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx2*(lx2*x2FD_j.d*(phic - phil) +
                              rx2*x2FD_jp1.d*(phir - phic));
#endif

  phil = 0.5*(pG->Phi[k-lx3][j][i] + pG->Phi[k  ][j][i]);
  phir = 0.5*(pG->Phi[k  ][j][i] + pG->Phi[k+rx3][j][i]);

#ifndef BAROTROPIC
  pG->U[k][j][i].E += dtodx3*(lx3*x3FD_k.d*(phic - phil) +
                              rx3*x3FD_kp1.d*(phir - phic));
#endif

#endif /* SELF_GRAVITY */


}


#if defined(MHD) || defined(RADIATION_MHD)
static void FixMHD(GridS *pG, int i, int j, int k, 
                      int lx1, int rx1, int lx2, int rx2, int lx3, int rx3)
{
 Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real lsf=1.0,rsf=1.0;

#ifdef CYLINDRICAL
  dtodx2 = pG->dt/(pG->ri[i]*pG->dx2);
#endif
  pG->B1i[k][j][i  ] -= dtodx3*(rx3*emf2D_kp1i   - lx3*emf2D_ki) -
                        dtodx2*(rx2*emf3D_jp1i   - lx2*emf3D_ji);
#ifdef CYLINDRICAL
  dtodx2 = pG->dt/(pG->ri[i+1]*pG->dx2);
#endif
  pG->B1i[k][j][i+1] -= dtodx3*(rx3*emf2D_kp1ip1 - lx3*emf2D_kip1) -
                        dtodx2*(rx2*emf3D_jp1ip1 - lx2*emf3D_jip1);
  pG->B2i[k][j  ][i] -= dtodx1*(rx1*emf3D_jip1   - lx1*emf3D_ji) -
                        dtodx3*(rx3*emf1D_kp1j   - lx3*emf1D_kj);
  pG->B2i[k][j+1][i] -= dtodx1*(rx1*emf3D_jp1ip1 - lx1*emf3D_jp1i) -
                        dtodx3*(rx3*emf1D_kp1jp1 - lx3*emf1D_kjp1);
#ifdef CYLINDRICAL
  rsf = (rx1>0) ? pG->ri[i+1]/pG->r[i] : pG->ri[i  ]/pG->r[i];
  lsf = (lx1>0) ? pG->ri[i  ]/pG->r[i] : pG->ri[i+1]/pG->r[i];
  dtodx2 = pG->dt/(pG->r[i]*pG->dx2);
#endif
  pG->B3i[k  ][j][i] -= dtodx2*(rx2*emf1D_kjp1   - lx2*emf1D_kj) -
                        dtodx1*(rsf*rx1*emf2D_kip1   - lsf*lx1*emf2D_ki);
  pG->B3i[k+1][j][i] -= dtodx2*(rx2*emf1D_kp1jp1 - lx2*emf1D_kp1j) -
                        dtodx1*(rsf*rx1*emf2D_kp1ip1 - lsf*lx1*emf2D_kp1i);
}

#endif /* End MHD */



void FOFC_Flux(const Cons1DS Ul, const Cons1DS Ur,
                   const Prim1DS Wl, const Prim1DS Wr,
                   const Real Bxi, Cons1DS *pFlux)
{
 Real sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe,pbl=0.0,pbr=0.0;
  Real asq,vaxsq=0.0,qsq,cfsq,cfl,cfr,bp,bm,ct2=0.0,tmp;
#ifndef ISOTHERMAL
  Real hroe;
#endif
#if  defined(MHD) || defined(RADIATION_MHD)
  Real b2roe,b3roe,x,y;
#endif
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Real dt=pFlux->d;
  Real Proe, Sigma_roe[NOPACITY];
  Real Erroe, Frroe[3], Edd[6];
  int DIM = (int)pFlux->Mx;
  int m;
#endif

  Real ev[NWAVE],al,ar;
  Real *pFl, *pFr, *pF;
  Cons1DS Fl,Fr;
  int n;


/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute Roe-averaged data from left- and right-states
 */

  sqrtdl = sqrt((double)Wl.d);
  sqrtdr = sqrt((double)Wr.d);
  isdlpdr = 1.0/(sqrtdl + sqrtdr);

  droe  = sqrtdl*sqrtdr;
  v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
  v2roe = (sqrtdl*Wl.Vy + sqrtdr*Wr.Vy)*isdlpdr;
  v3roe = (sqrtdl*Wl.Vz + sqrtdr*Wr.Vz)*isdlpdr;

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Proe = (sqrtdl*Wl.P + sqrtdr*Wr.P)*isdlpdr;
  for(m=0;m<NOPACITY;m++){
	Sigma_roe[m] = (sqrtdl*Wl.Sigma[m] + sqrtdr*Wr.Sigma[m] )*isdlpdr;
  }

  Erroe = (sqrtdl*Wl.Er + sqrtdr*Wr.Er)*isdlpdr;
  Frroe[0] = (sqrtdl*Wl.Fr1 + sqrtdr*Wr.Fr1)*isdlpdr;
  Frroe[1] = (sqrtdl*Wl.Fr2 + sqrtdr*Wr.Fr2)*isdlpdr;
  Frroe[2] = (sqrtdl*Wl.Fr3 + sqrtdr*Wr.Fr3)*isdlpdr;
  Edd[0] = (sqrtdl*Wl.Edd_11 + sqrtdr*Wr.Edd_11)*isdlpdr;
  Edd[1] = (sqrtdl*Wl.Edd_21 + sqrtdr*Wr.Edd_21)*isdlpdr;
  Edd[2] = (sqrtdl*Wl.Edd_22 + sqrtdr*Wr.Edd_22)*isdlpdr;
  Edd[3] = (sqrtdl*Wl.Edd_31 + sqrtdr*Wr.Edd_31)*isdlpdr;
  Edd[4] = (sqrtdl*Wl.Edd_32 + sqrtdr*Wr.Edd_32)*isdlpdr;
  Edd[5] = (sqrtdl*Wl.Edd_33 + sqrtdr*Wr.Edd_33)*isdlpdr;
#endif

/* The Roe average of the magnetic field is defined differently.  */

#if defined(MHD) || defined(RADIATION_MHD)
  b2roe = (sqrtdr*Wl.By + sqrtdl*Wr.By)*isdlpdr;
  b3roe = (sqrtdr*Wl.Bz + sqrtdl*Wr.Bz)*isdlpdr;
  x = 0.5*(SQR(Wl.By - Wr.By) + SQR(Wl.Bz - Wr.Bz))/(SQR(sqrtdl + sqrtdr));
  y = 0.5*(Wl.d + Wr.d)/droe;
  pbl = 0.5*(SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz));
#endif

/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */

#ifndef ISOTHERMAL
  hroe  = ((Ul.E + Wl.P + pbl)/sqrtdl + (Ur.E + Wr.P + pbr)/sqrtdr)*isdlpdr;
#endif

/*--- Step 3. ------------------------------------------------------------------
 * Compute eigenvalues using Roe-averaged values, needed in step 4.
 */

#ifdef HYDRO
#ifdef ISOTHERMAL
  esys_roe_iso_hyd(v1roe, v2roe, v3roe,       ev, NULL, NULL);
#else
  esys_roe_adb_hyd(v1roe, v2roe, v3roe, hroe, ev, NULL, NULL);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#if defined(MHD) || defined(RADIATION_MHD)
#ifdef ISOTHERMAL
 esys_roe_iso_mhd(droe,v1roe,v2roe,v3roe,     Bxi,b2roe,b3roe,x,y,ev,NULL,NULL);
#else
 esys_roe_adb_mhd(droe,v1roe,v2roe,v3roe,hroe,Bxi,b2roe,b3roe,x,y,ev,NULL,NULL);
#endif /* ISOTHERMAL */
#endif /* MHD */

/****************************
 radiation hydro and radiation mhd 
 *****************************/

#ifdef RADIATION_HYDRO
 esys_roe_rad_hyd(v1roe, v2roe, v3roe, hroe, dt, Proe, Erroe, Frroe, Edd, Sigma_roe, DIM, ev, NULL, NULL);
#endif

#ifdef RADIATION_MHD
esys_roe_adb_mhd(droe,v1roe,v2roe,v3roe,hroe,Bxi,b2roe,b3roe,x,y,ev,NULL,NULL);
#endif

/*--- Step 4. ------------------------------------------------------------------
 * Compute the max and min wave speeds
 */

/* left state */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  
 asq = eff_sound(Wl,dt,DIM);
/*
  asq = sqrt(Gamma * Wl.P / Wl.d);
*/
  asq = asq * asq;

#else 

#ifdef ISOTHERMAL
  asq = Iso_csound2;
#else
  asq = Gamma*Wl.P/Wl.d;
#endif

#endif /* radiation hydro and mhd */


#if defined(MHD) || defined(RADIATION_MHD)
  vaxsq = Bxi*Bxi/Wl.d;
  ct2 = (Ul.By*Ul.By + Ul.Bz*Ul.Bz)/Wl.d;
#endif
  qsq = vaxsq + ct2 + asq;
  tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((double)(tmp*tmp + 4.0*asq*ct2)));
  cfl = sqrt((double)cfsq);

/* right state */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  
 asq = eff_sound(Wr,dt,DIM);
/*
  asq = sqrt(Gamma * Wr.P / Wr.d);
*/	
  asq = asq * asq;

#else

#ifdef ISOTHERMAL
  asq = Iso_csound2;
#else
  asq = Gamma*Wr.P/Wr.d; 
#endif

#endif /* radiation hydro and mhd */



#if defined(MHD) || defined(RADIATION_MHD)
  vaxsq = Bxi*Bxi/Wr.d;
  ct2 = (Ur.By*Ur.By + Ur.Bz*Ur.Bz)/Wr.d;
#endif
  qsq = vaxsq + ct2 + asq;
  tmp = vaxsq + ct2 - asq;
  cfsq = 0.5*(qsq + sqrt((double)(tmp*tmp + 4.0*asq*ct2)));
  cfr = sqrt((double)cfsq);

/* take max/min of Roe eigenvalues and L/R state wave speeds */

  ar = MAX(ev[NWAVE-1],(Wr.Vx + cfr));
  al = MIN(ev[0]      ,(Wl.Vx - cfl));

  bp = MAX(ar, 0.0);
  bm = MIN(al, 0.0);

/*--- Step 5. ------------------------------------------------------------------
 * Compute L/R fluxes along the lines bm/bp: F_{L}-S_{L}U_{L}; F_{R}-S_{R}U_{R}
 */

  Fl.d  = Ul.Mx - bm*Ul.d;
  Fr.d  = Ur.Mx - bp*Ur.d;

  Fl.Mx = Ul.Mx*(Wl.Vx - bm);
  Fr.Mx = Ur.Mx*(Wr.Vx - bp);

  Fl.My = Ul.My*(Wl.Vx - bm);
  Fr.My = Ur.My*(Wr.Vx - bp);

  Fl.Mz = Ul.Mz*(Wl.Vx - bm);
  Fr.Mz = Ur.Mz*(Wr.Vx - bp);

#ifdef ISOTHERMAL
  Fl.Mx += Wl.d*Iso_csound2;
  Fr.Mx += Wr.d*Iso_csound2;
#else
  Fl.Mx += Wl.P;
  Fr.Mx += Wr.P;

  Fl.E  = Ul.E*(Wl.Vx - bm) + Wl.P*Wl.Vx;
  Fr.E  = Ur.E*(Wr.Vx - bp) + Wr.P*Wr.Vx;
#endif /* ISOTHERMAL */

#if defined(MHD) || defined(RADIATION_MHD)
  Fl.Mx -= 0.5*(Bxi*Bxi - SQR(Wl.By) - SQR(Wl.Bz));
  Fr.Mx -= 0.5*(Bxi*Bxi - SQR(Wr.By) - SQR(Wr.Bz));

  Fl.My -= Bxi*Wl.By;
  Fr.My -= Bxi*Wr.By;
    
  Fl.Mz -= Bxi*Wl.Bz;
  Fr.Mz -= Bxi*Wr.Bz;

#ifndef ISOTHERMAL
  Fl.E += (pbl*Wl.Vx - Bxi*(Bxi*Wl.Vx + Wl.By*Wl.Vy + Wl.Bz*Wl.Vz));
  Fr.E += (pbr*Wr.Vx - Bxi*(Bxi*Wr.Vx + Wr.By*Wr.Vy + Wr.Bz*Wr.Vz));
#endif /* ISOTHERMAL */

  Fl.By = Wl.By*(Wl.Vx - bm) - Bxi*Wl.Vy;
  Fr.By = Wr.By*(Wr.Vx - bp) - Bxi*Wr.Vy;

  Fl.Bz = Wl.Bz*(Wl.Vx - bm) - Bxi*Wl.Vz;
  Fr.Bz = Wr.Bz*(Wr.Vx - bp) - Bxi*Wr.Vz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

#ifdef CYLINDRICAL
#ifndef ISOTHERMAL
  Fl.Pflux = Wl.P;
  Fr.Pflux = Wr.P;
#if defined(MHD) || defined(RADIATION_MHD)
  Fl.Pflux += pbl;
  Fr.Pflux += pbr;
#endif /* MHD */
#endif /* ISOTHERMAL */
#endif /* CYLINDRICAL */

/*--- Step 6. ------------------------------------------------------------------
 * Compute the HLLE flux at interface.
 */

  pFl = (Real *)&(Fl);
  pFr = (Real *)&(Fr);
  pF  = (Real *)pFlux;
  tmp = 0.5*(bp + bm)/(bp - bm);
  for (n=0; n<(NWAVE+NSCALARS); n++){
    pF[n] = 0.5*(pFl[n] + pFr[n]) + (pFl[n] - pFr[n])*tmp;
  }

#ifdef CYLINDRICAL
  n = NWAVE+NSCALARS;
  pF[n] = 0.5*(pFl[n] + pFr[n]) + (pFl[n] - pFr[n])*tmp;
#endif 

  return;
}



#endif /* FIRST_ORDER_FLUX_CORRECTION */



#endif /* CTU_INTEGRATOR */
