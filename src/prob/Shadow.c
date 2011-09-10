#include "copyright.h"
/*==============================================================================
 * FILE: Shadow.c
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


static Real sigma0;
static Real T0;
static Real T1;
static Real d0;
static Real d1;

void constopa(const Real rho, const Real T, Real *Sigma_t, Real *Sigma_a, Real dSigma[4]);

void radMHD_Mat_inflow(MatrixS *pMat);


void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);

void radMHD_inflow2(GridS *pGrid);
void radMHD_rad_inflow2(GridS *pGrid);


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


void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  Real x1, x2, x3, lx, ly, lz, density;

/* Parse global variables of unit ratio */

  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
  Ncycle = par_getd_def("problem","Ncycle",20);
  TOL  = par_getd_def("problem","TOL",1.e-16);

  int ixs, jxs, kxs;
  /* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];

	

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




  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

  Real ybottom, xbottom;
  ybottom = pDomain->RootMinX[1];
  xbottom = pDomain->RootMinX[0];
	

	
	Real aaxis, baxis, delta,pressure,temperature;
	aaxis = lx/10.0;
	baxis = ly/10.0;

	sigma0 = 1.0;
	
	d0 = 1.0;
	d1 = 10.0;
	
	
	T0 = 1.0;
	T1 = 6.0;


/* Initialize the grid including the ghost cells.  */

	


  for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
		    for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
				delta = 10.0 * (x1 * x1/(aaxis * aaxis) + x2 * x2 / (baxis * baxis) - 1.0);
	
        	pGrid->U[k][j][i].d  = d0 + (d1-d0)/(1.0+exp(delta));
				
				density = pGrid->U[k][j][i].d;
				
				pressure = d0 * R_ideal * T0;
				
				temperature = pressure / pGrid->U[k][j][i].d;


	
	pGrid->U[k][j][i].M1  = 0.0;
	pGrid->U[k][j][i].M2  = 0.0;

	

        pGrid->U[k][j][i].M3 = 0.0;

        pGrid->U[k][j][i].E = pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[k][j][i].M1 * pGrid->U[k][j][i].M1 + 
																pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2) / pGrid->U[k][j][i].d;


/*	 pGrid->U[k][j][i].Edd_11 = tempEdd11; 
	 pGrid->U[k][j][i].Edd_22 = tempEdd22; 
*/
	 pGrid->U[k][j][i].Edd_11 = 1.0/3.0; 
	 pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
	 pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
	 pGrid->U[k][j][i].Sigma_t = sigma0 * density * density * pow(temperature, -3.5);
	 pGrid->U[k][j][i].Sigma_a = sigma0 * density * density * pow(temperature, -3.5);


#ifdef RADIATION_MHD
	 /* interface magnetic field */


          pGrid->B1i[k][j][i] = Bx0;
          pGrid->B2i[k][j][i] = By0;
          pGrid->B3i[k][j][i] = 0.0;

	  pGrid->U[k][j][i].E += B0 * B0 / 2.0;
        
#endif



	  pGrid->U[k][j][i].Er = pow(temperature,4.0);
	  pGrid->U[k][j][i].Fr1 = 0.0;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;


	  /* background Er */
	 
	
	 	
        }
      }
    }
	

#ifdef RADIATION_MHD

	 for (k=pGrid->ks; k<=pGrid->ke; k++) {
      for (j=pGrid->js; j<=pGrid->je; j++) {
        for (i=pGrid->is; i<=pGrid->ie; i++) {
	  pGrid->U[k][j][i].B1c = Bx0;
          pGrid->U[k][j][i].B2c = By0;
          pGrid->U[k][j][i].B3c = 0.0;


	}
	}
      }

	for(j=jl-nghost; j<=ju+nghost; j++){
		for(i=1; i<=nghost; i++){
			pGrid->B1i[0][j][il-i] = Bx0;
			pGrid->B1i[0][j][iu+i] = Bx0;
			pGrid->U[0][j][il-i].B1c = Bx0;
			pGrid->U[0][j][iu+i].B1c = Bx0;

			pGrid->B2i[0][j][il-i] = By0;
                        pGrid->B2i[0][j][iu+i] = By0;
                        pGrid->U[0][j][il-i].B2c = By0;
                        pGrid->U[0][j][iu+i].B2c = By0;

			
		}
	}

	
	for(i=il-nghost; i<=iu+nghost; i++){
		for(j=1; j<=nghost; j++){
			pGrid->B1i[0][jl-j][i] = Bx0;
			pGrid->B1i[0][ju+j][i] = Bx0;
			pGrid->U[0][jl-j][i].B1c = Bx0;
			pGrid->U[0][ju+j][i].B1c = Bx0;

			pGrid->B2i[0][jl-j][i] = By0;
                        pGrid->B2i[0][ju+j][i] = By0;
                        pGrid->U[0][jl-j][i].B2c = By0;
                        pGrid->U[0][ju+j][i].B2c = By0;

		}
	}


#endif



	Opacity = constopa;


	bvals_mhd_fun(pDomain, left_x1, radMHD_inflow);
	bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow);

/*
	bvals_mhd_fun(pDomain, left_x2, radMHD_inflow2);
	bvals_rad_fun(pDomain, left_x2, radMHD_rad_inflow2);
	
	bvals_mhd_fun(pDomain, right_x2, radMHD_inflow2);
	bvals_rad_fun(pDomain, right_x2, radMHD_rad_inflow2);

*/

/* Boundary condition for trasnfer module */


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
 
  for(j=pRG->js-1; j<=pRG->je+1; j++) {
/* incident radiation at left boundary */
    for(ifr=0; ifr<nf; ifr++)
      for(m=0; m<nang/2; m++) {
		  /* Corresponding to radiation temperature */
		  if(m==5){
				pRG->l1imu[ifr][pRG->ks][j][0][m] = pow(6.0, 4.0) * Thermal_B(pGrid, ifr, pGrid->is-1,j+nghost-1, pGrid->ks);
			    pRG->l1imu[ifr][pRG->ks][j][2][m] = pow(6.0, 4.0) * Thermal_B(pGrid, ifr, pGrid->is-1,j+nghost-1, pGrid->ks);
		  }
		  else{
				pRG->l1imu[ifr][pRG->ks][j][0][m] = 0.0;
				pRG->l1imu[ifr][pRG->ks][j][2][m] = 0.0;
		  }
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

	if(m == 5){
		pRG->r2imu[ifr][pRG->ks][i][2][m] = pow(6.0, 4.0) * Thermal_B(pGrid, ifr, pGrid->is-1,nghost, pGrid->ks);
		pRG->r2imu[ifr][pRG->ks][i][3][m] = 0.0;
		 
		 	  
		pRG->l2imu[ifr][pRG->ks][i][0][m] = pow(6.0, 4.0) * Thermal_B(pGrid, ifr, pGrid->is-1,nghost, pGrid->ks);
		pRG->l2imu[ifr][pRG->ks][i][1][m] = 0.0;
	}
	else{
		pRG->r2imu[ifr][pRG->ks][i][2][m] = 0.0;
		pRG->r2imu[ifr][pRG->ks][i][3][m] = 0.0;
		
		
		pRG->l2imu[ifr][pRG->ks][i][0][m] = 0.0;
		pRG->l2imu[ifr][pRG->ks][i][1][m] = 0.0;
		
			  
	}
      }
  }

get_thermal_source = Thermal_B;
get_thermal_fraction = const_eps;
get_total_opacity = transfer_opacity;

#endif





  return;
}



void constopa(const Real rho, const Real T, Real *Sigma_t, Real *Sigma_a, Real dSigma[4]){
	if(Sigma_t != NULL)
		*Sigma_t = sigma0 * rho * rho * pow(T, -3.5);
	
	if(Sigma_a != NULL)
		*Sigma_a = sigma0 * rho * rho * pow(T, -3.5);

	if(dSigma != NULL){
		dSigma[0] = 2.0 * sigma0 * rho * pow(T, -3.5);
		dSigma[1] = 2.0 * sigma0 * rho * pow(T,-3.5);
		dSigma[2] = -3.5 * sigma0 * rho * rho * pow(T, -4.5);
		dSigma[3] = -3.5 * sigma0 * rho * rho * pow(T, -4.5);
	}
	

 return; 

}


/* right hand side of x2 direction */

void radMHD_inflow(GridS *pGrid)
{
  	int i, je,j, ju, jl, js;
	int ks, is,ie, k, ke, ku;
	Real pressure, x1, x2, x3, density;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;


	
	
	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;


#if defined(MHD) || defined(RADIATION_MHD)
/* B1i is not set at i=is-nghost */
/* Assuming x periodic boundary is already applied */
/* Boundary condition for j=0, assuming B3i = 0 in 2D case */
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
	/* first update B1i */
      for (i=is-(nghost-1); i<=ie+nghost; i++) {
 /*       pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-1][i] + dxratio * (pGrid->B2i[k][je+j][i] - pGrid->B2i[k][je+j][i-1]);
*/
/*
	 pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-1][i];
*/
	
	pGrid->B1i[k][je+j][i] = Bx0;
      }
	
	/* Then update B2i */
	/* j=je+1 is not a boundary condition for the interface field B2i */
	/* bounary condition from div B=0 */

	if(j<nghost){
		for (i=is-nghost; i<=ie+nghost-1; i++) {
/*        		pGrid->B2i[k][je+j+1][i] = pGrid->B2i[k][je+j][i] - dxratio * (pGrid->B1i[k][je+j][i+1] - pGrid->B1i[k][je+j][i]);
*/
/*			pGrid->B2i[k][je+j+1][i] = pGrid->B2i[k][je+1][i];
*/
			pGrid->B2i[k][je+j+1][i] = By0;
	
		
      		}
	}

    }
  }


  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je][i];
      }
    }
  }

	/* Now update cell centered magnetic field */
	 for(i=is-nghost; i<=ie+nghost; i++){
	    for (j=1;  j<nghost;  j++) {

		pGrid->U[ks][je+j][i].B1c = 0.5 * (pGrid->B1i[ks][je+j][i] + pGrid->B1i[ks][je+j][i+1]);
		pGrid->U[ks][je+j][i].B2c = 0.5 * (pGrid->B2i[ks][je+j+1][i] + pGrid->B2i[ks][je+j][i]);
		pGrid->U[ks][je+j][i].B3c = 0.5 * (pGrid->B3i[ks][je+j][i] + pGrid->B3i[ks][je+j][i]);
		
		}
	}


#endif /* MHD */
	

 
    for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
	

		cc_pos(pGrid, is-i, j,ks, &x1, &x2, &x3);

		pGrid->U[ks][j][is-i].d  = d0;
			
			density = d0;
		
     	pGrid->U[ks][j][is-i].M1 = 0.0;
		pGrid->U[ks][j][is-i].M2 = 0.0;

		
		
		pGrid->U[ks][j][is-i].M3 = 0.0;
			
		pressure = d0 * T0 * R_ideal;

		pGrid->U[ks][j][is-i].E = pressure / (Gamma - 1.0) + 0.5 * (pGrid->U[ks][j][is-i].M1 * pGrid->U[ks][j][is-i].M1 + 
																		pGrid->U[ks][j][is-i].M2 * pGrid->U[ks][j][is-i].M2) / pGrid->U[ks][j][is-i].d;
			
#ifdef RADIATION_MHD
		pGrid->U[ks][j][is-i].E += 0.5 * (pGrid->U[ks][j][is-i].B1c * pGrid->U[ks][j][is-i].B1c + pGrid->U[ks][j][is-i].B2c * pGrid->U[ks][j][is-i].B2c + pGrid->U[ks][j][is-i].B3c * pGrid->U[ks][j][is-i].B3c);

#endif

		pGrid->U[ks][j][is-i].Sigma_t = sigma0 * density * density * pow(T0, -3.5);
		pGrid->U[ks][j][is-i].Sigma_a = sigma0 * density * density * pow(T0, -3.5);
      		
      }
    }
  
}


void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, je,j, ju, jl, js;
	int ks, is,ie,k, ke, ku;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	ke = pGrid->ke;
	is = pGrid->is;
	ie = pGrid->ie;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	Real x1, x2, x3;

	
	for(j=js-nghost; j<=je+nghost; j++){
	    for (i=1;  i<=nghost;  i++) {
			
		cc_pos(pGrid, is-i, j,ks, &x1, &x2, &x3);

		pGrid->U[ks][j][is-i].Edd_11 = pGrid->U[ks][j][is].Edd_11;
		pGrid->U[ks][j][is-i].Edd_22 = pGrid->U[ks][j][is].Edd_22;
		pGrid->U[ks][j][is-i].Edd_21 = pGrid->U[ks][j][is].Edd_21;
	
		pGrid->U[ks][j][is-i].Er  = pow(T1,4.0);
		pGrid->U[ks][j][is-i].Fr1 = 0.0;
		pGrid->U[ks][j][is-i].Fr2 = 0.0;

		pGrid->U[ks][j][is-i].Fr3 = 0.0;		

    }
	}

  
}


void radMHD_inflow2(GridS *pGrid)
{
}

void radMHD_rad_inflow2(GridS *pGrid)
{
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



void radMHD_Mat_inflow(MatrixS *pMat)
{
  	int i, je,j, ju, jl, js;
	int ks, is,ie,k, ke, ku;
	je = pMat->je;
	js = pMat->js;
  	ks = pMat->ks;
	ke = pMat->ke;
	is = pMat->is;
	ie = pMat->ie;

	ju = pMat->je + Matghost;
  	jl = pMat->js - Matghost;

	Real x1, x2, x3;

	
	for(j=js-Matghost; j<=je+Matghost; j++){
	    for (i=1;  i<=Matghost;  i++) {
			

		pMat->U[ks][j][is-i].Edd_11 = pMat->U[ks][j][is].Edd_11;
		pMat->U[ks][j][is-i].Edd_22 = pMat->U[ks][j][is].Edd_22;
		pMat->U[ks][j][is-i].Edd_21 = pMat->U[ks][j][is].Edd_21;
	
		pMat->U[ks][j][is-i].Er  = pow(T1,4.0);
		pMat->U[ks][j][is-i].Fr1 = 0.0;
		pMat->U[ks][j][is-i].Fr2 = 0.0;

		pMat->U[ks][j][is-i].Fr3 = 0.0;		

    }
	}

  
}

void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{

	*Mat_BCFun = radMHD_Mat_inflow;

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


void Userwork_in_formal_solution(DomainS *pD)
{
  return;
}

void problem_write_restart(MeshS *pM, FILE *fp)
{


/*
	fwrite(&Gamma,sizeof(Real),1,fp);
 	fwrite(&Prat,sizeof(Real),1,fp);
	fwrite(&Crat,sizeof(Real),1,fp);
	fwrite(&R_ideal,sizeof(Real),1,fp);
 	fwrite(&kappae,sizeof(Real),1,fp);
	fwrite(&kappa0,sizeof(Real),1,fp);
	fwrite(&grav,sizeof(Real),1,fp);	
	fwrite(&consFr,sizeof(Real),1,fp);
	fwrite(&B0,sizeof(Real),1,fp);
	fwrite(&Bx0,sizeof(Real),1,fp);
	fwrite(&By0,sizeof(Real),1,fp);
	fwrite(&rho_up[4],sizeof(Real),1,fp);
	fwrite(&Energy_up[4],sizeof(Real),1,fp);
	fwrite(&Er_up[4],sizeof(Real),1,fp);
	fwrite(&rho_lower[4],sizeof(Real),1,fp);
	fwrite(&Energy_lower[4],sizeof(Real),1,fp);
	fwrite(&Er_lower[4],sizeof(Real),1,fp);
	
*/

  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
	int i;
/*

	fread(&Gamma,sizeof(Real),1,fp);
	fread(&Prat,sizeof(Real),1,fp);
	fread(&Crat,sizeof(Real),1,fp);
	fread(&R_ideal,sizeof(Real),1,fp);
	fread(&kappae,sizeof(Real),1,fp);
	fread(&kappa0,sizeof(Real),1,fp);
	fread(&grav,sizeof(Real),1,fp);	
	fread(&consFr,sizeof(Real),1,fp);
	fread(&B0,sizeof(Real),1,fp);
	fread(&Bx0,sizeof(Real),1,fp);
	fread(&By0,sizeof(Real),1,fp);
	fread(&rho_up[4],sizeof(Real),1,fp);
	fread(&Energy_up[4],sizeof(Real),1,fp);
	fread(&Er_up[4],sizeof(Real),1,fp);
	fread(&rho_lower[4],sizeof(Real),1,fp);
	fread(&Energy_lower[4],sizeof(Real),1,fp);
	fread(&Er_lower[4],sizeof(Real),1,fp);

	


	Opacity = constopa;


	bvals_mhd_fun(&(pM->Domain[0][0]), left_x1, radMHD_inflow);
	bvals_rad_fun(&(pM->Domain[0][0]), left_x1, radMHD_rad_inflow);

*/

	

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
