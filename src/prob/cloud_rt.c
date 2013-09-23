#include "copyright.h"
/*==============================================================================
 * FILE: cloud_rt.c
 *
 * PURPOSE: Problem generator to model an irradiated cloud in 2D using the ray
 * tracing algorithm.  This version is adapted by S.W. Davis from the generator 
 * written by Y.-F. Jiang to model irradiation without using ray tracing
 *
 * configure:
 *  --with-problem=cloud_rt --enable-radiation-transfer  --with-gas=rad_hydro
 *  --enable-ray-tracing 
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
#if defined(RADIATION_TRANSFER) && defined(RAY_TRACING)

static Real sigma0;
static Real T0;
static Real T1;
static Real d0;
static Real d1;
static Real kappa;
static Real consFr;

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
void constopa(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY]);
void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
#ifdef MATRIX_MULTIGRID
void radMHD_Mat_inflowis(MatrixS *pMat);
#endif
#endif

static Real const_B(const GridS *pG, const int ifr, const int i,
		    const int j, const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);
static Real const_chi(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);

/*=========================== PUBLIC FUNCTIONS =================================
 *============================================================================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  RadGridS *pRG = (pDomain->RadGrid);

  int nf=pRG->nf, nang=pRG->nang;
  int ifr, l, m;
  int i, j, k, iu, il, ju, jl, ku, kl;
  Real x1, x2, x3, lx, ly, lz, density;
  Real ybottom, xbottom;
  Real aaxis, baxis, delta, pressure, temperature;
  int ixs, jxs, kxs;
  Real rt_rat, flux_sc, taub, taut, x10;

/* Parse global variables of unit ratio */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  Eratio = par_getd_def("problem","Eratio",1.00);
  Erflag = par_getd_def("problem","Erflag",1);
  Ncycle = par_getd_def("problem","Ncycle",6);
  TOL  = par_getd_def("problem","TOL",1.e-10);
  sigma0 = par_getd_def("problem","sigma0",1.0);
  consFr = par_getd_def("problem","flux",10.0);
  if (myID_Comm_world == 0) 
    printf("Eratio: %f  Erflag: %d\n",Eratio,Erflag);
#endif


  R_ideal = par_getd("problem","R_ideal");
  kappa = par_getd_def("problem","kappa",0.0);
  d0 = par_getd("problem","drare");
  d1 = par_getd("problem","dcloud");
  T0 = par_getd("problem","Trare");
  rt_rat = par_getd_def("problem","rt_rat",1.0);
  flux_sc = consFr * (1.0 - rt_rat);
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


  ybottom = pDomain->RootMinX[1];
  xbottom = pDomain->RootMinX[0];

  aaxis = lx/20.0;
  baxis = ly/10.0;
  
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

        pGrid->U[k][j][i].E = pressure / (Gamma_1) + 0.5 * (SQR(pGrid->U[k][j][i].M1) +
			      SQR(pGrid->U[k][j][i].M2)) / pGrid->U[k][j][i].d;
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	pGrid->U[k][j][i].Edd_11 = 1.0/3.0; 
	pGrid->U[k][j][i].Edd_22 = 1.0/3.0;
	pGrid->U[k][j][i].Edd_21 = 0.0; /* Set to be a constant in 1D. To be modified later */
	pGrid->U[k][j][i].Sigma[0] = kappa * density;
	pGrid->U[k][j][i].Sigma[1] = sigma0 * SQR(density) * pow(temperature, -3.5);
	pGrid->U[k][j][i].Sigma[2] = sigma0 * SQR(density) * pow(temperature, -3.5);
	pGrid->U[k][j][i].Sigma[3] = sigma0 * SQR(density) * pow(temperature, -3.5);

	pGrid->U[k][j][i].Er = pow(temperature,4.0);
	pGrid->U[k][j][i].Fr1 = 0.0;
	pGrid->U[k][j][i].Fr2 = 0.0;
	pGrid->U[k][j][i].Fr3 = 0.0;
#endif
#ifdef RADIATION_TRANSFER
	pGrid->tgas[k][j][i] = temperature; 
#endif
      }
    }
  }
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Opacity = constopa;
  
  bvals_mhd_fun(pDomain, left_x1, radMHD_inflow);
  bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow);
#endif

/* ------- Initialize boundary emission ---------------------------------- */

  for(ifr=0; ifr<nf; ifr++) {
    for(j=pRG->js-1; j<=pRG->je+1; j++) {
      for(m=0; m<nang; m++) {
/* incident thermal radiation at left and right boundaries */
/*	pRG->Ghstl1i[ifr][pRG->ks][j][0][m] = const_B(pGrid,ifr,pGrid->is-1,j+nghost-1,ku);
	pRG->Ghstl1i[ifr][pRG->ks][j][2][m] = const_B(pGrid,ifr,pGrid->is-1,j+nghost-1,ku);
	pRG->Ghstr1i[ifr][pRG->ks][j][1][m] = const_B(pGrid,ifr,pGrid->ie+1,j+nghost-1,ku);
	pRG->Ghstr1i[ifr][pRG->ks][j][3][m] = const_B(pGrid,ifr,pGrid->ie+1,j+nghost-1,ku);*/
	pRG->Ghstl1i[ifr][pRG->ks][j][0][m] = 2.0 * flux_sc;
	pRG->Ghstl1i[ifr][pRG->ks][j][2][m] = 2.0 * flux_sc;
      }
    }

    if (pRG->lx1_id == -1) {
      for(m=0; m<nang; m++) {
	pRG->r2imu[ifr][pRG->ks][0][0][m] = 2.0 * flux_sc;
	pRG->l2imu[ifr][pRG->ks][0][2][m] = 2.0 * flux_sc;
      }
    }
  
    for(i=pRG->is-1; i<=pRG->ie+1; i++) {
/* incident thermal radiation at upper and lower boundaries */
      for(m=0; m<nang; m++) {
	/*	pRG->Ghstl2i[ifr][pRG->ks][i][0][m] = const_B(pGrid,ifr,i+nghost-1,pGrid->js-1,ku);
	pRG->Ghstl2i[ifr][pRG->ks][i][1][m] = const_B(pGrid,ifr,i+nghost-1,pGrid->js-1,ku);
	pRG->Ghstr2i[ifr][pRG->ks][i][2][m] = const_B(pGrid,ifr,i+nghost-1,pGrid->je+1,ku);
	pRG->Ghstr2i[ifr][pRG->ks][i][3][m] = const_B(pGrid,ifr,i+nghost-1,pGrid->je+1,ku);*/

	/*cc_pos(pGrid, i+nghost-1, nghost,0, &x1, &x2, &x3);
	if (pRG->lx2_id == -1) {
	  taub = (pGrid->U[0][nghost-1][i+nghost-1].Sigma[0] + pGrid->U[0][nghost-1][i+nghost-1].Sigma[1]) *
	  (x1-xbottom);
	  pRG->Ghstl2i[ifr][pRG->ks][i][0][m] = 2.0 * flux_sc/ exp(taub);
	}
	if (pRG->rx2_id == -1) {
	  taut = (pGrid->U[0][pRG->je+nghost][i+nghost-1].Sigma[0] + pGrid->U[0][pRG->je+nghost][i+nghost-1].Sigma[1]) *
	  (x1-xbottom);
	  pRG->Ghstr2i[ifr][pRG->ks][i][2][m] = 2.0 * flux_sc/ exp(taut);
	}*/


      }
    }
  }

  get_thermal_source = const_B;
  get_thermal_fraction = const_eps;
  get_total_opacity = const_chi;

/* ------- Initialize ray tracing boundary ---------------------------------- */

/* Redefine jl, ju, kl, ku for RadGrid */
  jl = pRG->js, ju = pRG->je;
  kl = pRG->ks, ku = pRG->ke;
  if (pRG->Nx[1] > 1) {
    jl--; ju++; 
  }
  if (pRG->Nx[2] > 1) {
    kl--; ku++;
  }

  for(ifr=0; ifr<pRG->nf_rt; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++) {	
	pRG->H[ifr][k][j][pRG->is-1] = consFr * rt_rat;
      }}}
/* enroll ray tracing opacity functions */
  //get_raytrace_thermal_fraction = const_eps;
  //get_raytrace_opacity = const_chi;

  return;
}

#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)

void constopa(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY])
{

/* Sigma[0-NOPACITY] are: Sigma_sF, Sigma_aF, Sigma_aP, Sigma_aE respectively */
/* dSigma[0-2*NOPACITY] are: dSigma_sF/drho, dSigma_aF/drho, dSigma_aP/drho, dSigma_aE/drho */
/* 			     dSigma_sF/dT,   dSigma_aF/dT,   dSigma_aP/dT,   dSigma_aE/dT */

/* When pressure becomes negative, we do not include radiation source term */
  Real Tpower, Tpower1;

  if((rho * T * R_ideal > 2.0 * TINY_NUMBER) && (rho > 0.0)) {	

    /* Tpower = T^3.5 , Tpower1 = T^4.5 */
    Tpower = 1.0 / (T * T * T * sqrt(T));
    Tpower1 = Tpower / T;
	
    if(Sigma != NULL) {
      Sigma[0] =  kappa * rho;
      Sigma[1] =  sigma0 * rho * rho * Tpower;
      Sigma[2] =  sigma0 * rho * rho * Tpower;
      Sigma[3] =  sigma0 * rho * rho * Tpower;
    }	
    if(dSigma != NULL) {
      dSigma[0] = kappa;
      dSigma[1] = 2.0 * rho * sigma0 * Tpower;
      dSigma[2] = 2.0 * rho * sigma0 * Tpower;
      dSigma[3] = 2.0 * rho * sigma0 * Tpower;

      dSigma[4] = 0.0;
      dSigma[5] = -3.5 * sigma0 * rho * rho * Tpower1;
      dSigma[6] = -3.5 * sigma0 * rho * rho * Tpower1;
      dSigma[7] = -3.5 * sigma0 * rho * rho * Tpower1;		
    }
  } else {

    if(Sigma != NULL) {
      Sigma[0] =  0.0;
      Sigma[1] =  0.0;
      Sigma[2] =  0.0;
      Sigma[3] =  0.0;
    }
	
    if(dSigma != NULL) {
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

  for(j=js-nghost; j<=je+nghost; j++){
    for (i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, is-i, j,ks, &x1, &x2, &x3);

      pGrid->U[ks][j][is-i].d  = d0;
      density = d0;		
      pGrid->U[ks][j][is-i].M1 = 0.0;
      pGrid->U[ks][j][is-i].M2 = 0.0;
      pGrid->U[ks][j][is-i].M3 = 0.0;			
      pressure = d0 * T0 * R_ideal;
      pGrid->U[ks][j][is-i].E = pressure / (Gamma - 1.0) + 0.5 * (SQR(pGrid->U[ks][j][is-i].M1) +
				SQR(pGrid->U[ks][j][is-i].M2)) / pGrid->U[ks][j][is-i].d;

#ifdef RADIATION_MHD
      pGrid->U[ks][j][is-i].E += 0.5 * (pGrid->U[ks][j][is-i].B1c * pGrid->U[ks][j][is-i].B1c + pGrid->U[ks][j][is-i].B2c * pGrid->U[ks][j][is-i].B2c + pGrid->U[ks][j][is-i].B3c * pGrid->U[ks][j][is-i].B3c);

#endif
      pGrid->U[ks][j][is-i].Sigma[0] = 0.0;
      pGrid->U[ks][j][is-i].Sigma[1] = sigma0 * density * density * pow(T0, -3.5);
      pGrid->U[ks][j][is-i].Sigma[2] = sigma0 * density * density * pow(T0, -3.5);
      pGrid->U[ks][j][is-i].Sigma[1] = sigma0 * density * density * pow(T0, -3.5);      		
    }
  }
  
}


void radMHD_rad_inflow(GridS *pGrid)
{
  int i, je,j, ju, jl, js;
  int ks, is,ie,k, ke, ku;
  Real x1, x2, x3;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  ke = pGrid->ke;
  is = pGrid->is;
  ie = pGrid->ie;
  
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;
  
  for(j=js-nghost; j<=je+nghost; j++){
    for (i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, is-i, j,ks, &x1, &x2, &x3);

      pGrid->U[ks][j][is-i].Edd_11 = pGrid->U[ks][j][is].Edd_11;
      pGrid->U[ks][j][is-i].Edd_22 = pGrid->U[ks][j][is].Edd_22;
      pGrid->U[ks][j][is-i].Edd_21 = pGrid->U[ks][j][is].Edd_21;
      
      pGrid->U[ks][j][is-i].Er  = consFr;
      //pGrid->U[ks][j][is-i].Er  = pow(T1,4.0);
      pGrid->U[ks][j][is-i].Fr1 =  0.0;
      pGrid->U[ks][j][is-i].Fr2 = 0.0;
      pGrid->U[ks][j][is-i].Fr3 = 0.0;		
    }
  }
  
  return;
}




#ifdef MATRIX_MULTIGRID

void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{

  *Mat_BCFun = radMHD_Mat_inflowis;

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

void radMHD_Mat_inflowis(MatrixS *pMat)
{
  int i, j, k, je, ju, jl, js, ks, ke;
  int is,ie, il, iu;
  Real x1, x2, x3;
  Real vx, vy, vz, Fr0z, Sigma_t;

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

  for(j=js-Matghost; j<=je+Matghost; j++){
    for(i=1; i<=Matghost; i++){

      pMat->Ugas[ks][j][is-i].Edd_11 = pMat->Ugas[ks][j][is].Edd_11;
      pMat->Ugas[ks][j][is-i].Edd_22 = pMat->Ugas[ks][j][is].Edd_22;
      pMat->Ugas[ks][j][is-i].Edd_21 = pMat->Ugas[ks][j][is].Edd_21;
      pMat->Ugas[ks][j][is-i].Edd_31 = pMat->Ugas[ks][j][is].Edd_31;
      pMat->Ugas[ks][j][is-i].Edd_32 = pMat->Ugas[ks][j][is].Edd_32;
      pMat->Ugas[ks][j][is-i].Edd_33 = pMat->Ugas[ks][j][is].Edd_33;
				
      pMat->Ugas[ks][j][is-i].V1 =  pMat->Ugas[ks][j][is].V1;
      pMat->Ugas[ks][j][is-i].V2 = -pMat->Ugas[ks][j][is].V2;
      pMat->Ugas[ks][j][is-i].T4 =  pMat->Ugas[ks][j][is].T4;
      pMat->Ugas[ks][j][is-i].Sigma[0] = pMat->Ugas[ks][j][is].Sigma[0];
      pMat->Ugas[ks][j][is-i].Sigma[1] = pMat->Ugas[ks][j][is].Sigma[1];
      pMat->Ugas[ks][j][is-i].Sigma[2] = pMat->Ugas[ks][j][is].Sigma[2];
      pMat->Ugas[ks][j][is-i].Sigma[3] = pMat->Ugas[ks][j][is].Sigma[3];

      //if((pMat->bgflag) || (pMat->Nx[0] < pMat->RootNx[0])){
	pMat->U[ks][j][is-i].Er  = 0.0;
	pMat->U[ks][j][is-i].Fr1 = 0.0;
	pMat->U[ks][j][is-i].Fr2 = 0.0;
	pMat->U[ks][j][is-i].Fr3 = 0.0;
	//} else {
	//	pMat->U[ks][j][is-i].Fr1 = pMat->U[ks][j][is].Fr1;
	//pMat->U[ks][j][is-i].Fr2 = pMat->U[ks][j][is].Fr2;	
	//pMat->U[ks][j][is-i].Er = pow(T1, 4.0);
	//}
    }
  }

  return;
}

#endif /* end multi_grid */

#endif /* end radMHD or radhydro */

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * Userwork_after_formal_solution  - problem specific work after formal solution
 *----------------------------------------------------------------------------*/

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

void Userwork_after_formal_solution(DomainS *pD)
{

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

static Real const_B(const GridS *pG, const int ifr, const int i,
		    const int j, const int k)
{
  return pow(pG->tgas[k][j][i], 4) / (4.0 * PI);

}

static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k)
{
  Real eps;
  eps = pG->U[k][j][i].Sigma[1] / (pG->U[k][j][i].Sigma[1] + 
        pG->U[k][j][i].Sigma[0]);

  return eps;
}

static Real const_chi(const GridS *pG, const int ifr, const int i,
		      const int j, const int k)
{
  return (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[1]);

}

#endif /* RADIATION_TRANSFER */
