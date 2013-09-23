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



/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/



static Real grav_pot(const Real x1, const Real x2, const Real x3);

static Real kappaes = 1.0;
static Real kappaff = 0.0;
static Real grav = 16.815;
static Real consFr;
static Real Lx, Ly, Lz;
static Real d1, d2, T1;
static Real pressure;
static Real Ertop, Erbottom, Ermid;
static Real ytop, ybtm, yintt, yintb;

void constopa(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY]);
static double ran2(long int *idum);

void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
void radMHD_inflow2(GridS *pGrid);
void radMHD_rad_inflow2(GridS *pGrid);

void radMHD_Mat_inflowi2(MatrixS *pMat);
void radMHD_Mat_inflowo2(MatrixS *pMat);


#ifdef RADIATION_TRANSFER

void const_H_ix2(RadGridS *pRG, int ifr);
void const_J_ox2(RadGridS *pRG, int ifr);
static Real Thermal_B(const GridS *pG, const int ifr, const int i, const int j, 
		    const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i, const int j, 
		      const int k);
static Real transfer_opacity(const GridS *pG, const int ifr, const int i, const int j, 
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
  int ixs,jxs,kxs,i,j,k,ipert,ifield;
  int nx3;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax;
  Real den = 1.0, rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real beta=1.0,B0,kx,ky,kz,amp=0.04;
  Real Itop, Ibtm, Irat;
  Real alpha;
  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  double rval;
  static int frst=1;  /* flag so new history variables enrolled only once */
  int iprob = 2;
  int ID, lx2, rx2, NGx, lines, NGy;
  Real tempEdd22, tempEdd11;
#ifdef RADIATION_TRANSFER
  RadGridS *pRG = (pDomain->RadGrid);
  int nang=pRG->nang;
  int noct = pRG->noct;
  int ifr = 0, l, m;
  int il, iu, jl, ju, kl, ku;
#endif

/* Parse global variables of unit ratio */
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  Ncycle = par_getd_def("problem","Ncycle",60);
  TOL  = par_getd_def("problem","TOL",1.e-10);
  Ncycle = 120.0;
#endif

/* obtain problem paramters from input file */
  R_ideal = par_getd("problem","R_ideal");	
  alpha  = par_getd_def("problem","alpha",1.0);
  yintt = par_getd("problem","shl_top");
  yintb = par_getd("problem","shl_btm");
  Irat = par_getd_def("problem","Irat",0.0);
  
  d1 = par_getd("problem","drare");
  d2 = par_getd("problem","dshell");
  T1 = par_getd("problem","Trare");
  
  consFr = alpha * grav / (Prat * (kappaes + kappaff));
  pressure = d1 * T1;
	
  Itop = Irat * consFr;
  ytop = pDomain->RootMaxX[1];
  ybtm = pDomain->RootMinX[1];	
  
  Ertop = ytop + (0.57735 + 0.33333 * Itop / consFr) / ((kappaes + kappaff) * d1);
  Ermid = yintt + (Ertop - yintt) * d1 / d2;
  Erbottom = (yintt - yintb) * d2 / d1 + Ertop + yintb - yintt;
  
  Ibtm = d1 * (kappaes + kappaff) * consFr * (Erbottom - ybtm) * 3.0;

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));
  printf("%d %d %d %d %ld\n",myID_Comm_world,ixs,jxs,kxs,iseed);
/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
 

/* Rescale amp to sound speed for ipert 2,3 */
  
  tempEdd11 = 1.0/3.0;
  tempEdd22 = 1.0/3.0;
  for (k=ks; k<=ke; k++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
    for (j=js; j<=je; j++) {
      //      for (i=is-nghost; i<=ie+nghost; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	pressure = d1 * T1;

	pGrid->U[k][j][i].d  = d1;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
	if ((x2 < yintt) && (x2 > yintb)) {
	  pGrid->U[k][j][i].d = d2;
          pGrid->U[k][j][i].M2 *= d2/d1;
	}

        pGrid->U[k][j][i].E = pressure/(Gamma - 1.0)+0.5*((pGrid->U[k][j][i].M1 * 
			      pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d + 
			      (pGrid->U[k][j][i].M2 * pGrid->U[k][j][i].M2) / 
			      pGrid->U[k][j][i].d);
	pGrid->U[k][j][i].Edd_11 = tempEdd11;
	pGrid->U[k][j][i].Edd_22 = tempEdd22;
	pGrid->U[k][j][i].Edd_21 = 0.0;
	pGrid->U[k][j][i].Sigma[0] = kappaes * pGrid->U[k][j][i].d;
	pGrid->U[k][j][i].Sigma[1] = kappaff * pGrid->U[k][j][i].d * pGrid->U[k][j][i].d * pow(pressure/pGrid->U[k][j][i].d,-3.5);
	pGrid->U[k][j][i].Sigma[2] = kappaff * pGrid->U[k][j][i].d * pGrid->U[k][j][i].d * pow(pressure/pGrid->U[k][j][i].d,-3.5);
	pGrid->U[k][j][i].Sigma[3] = kappaff * pGrid->U[k][j][i].d * pGrid->U[k][j][i].d * pow(pressure/pGrid->U[k][j][i].d,-3.5);

	pGrid->U[k][j][i].Fr1 = 0.0;
	pGrid->U[k][j][i].Fr2 = consFr;
	pGrid->U[k][j][i].Fr3 = 0.0;	
			
	if(x2 < yintb) {	  
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
                                 (Erbottom - x2) / pGrid->U[k][j][i].Edd_22;
	} else if(x2 < yintt) {
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                         (Ermid - x2) / pGrid->U[k][j][i].Edd_22;
	} else {
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                         (Ertop - x2) / pGrid->U[k][j][i].Edd_22;	
	}

	/* perturb density */

	if(fabs(x2) < 0.25 * Ly){

	  if (iprob == 1) {
	    pGrid->U[k][j][i].d *= (1.0 + amp/4.0*
	      (1.0+cos(2.0*PI*x1/Lx))*(1.0+cos(PI*x2)));
	  }
	  else if(iprob == 2){
	    pGrid->U[k][j][i].d *= (1.0 + 0.25 * amp*(ran2(&iseed) - 0.5)*(1.0+cos(PI*x2)));
	  }
	  else{
	   
	  }				
	}


      }
    }
  }

  /* Now ghost zone in the up J direction */
  for(i=is-nghost; i<=ie+nghost; i++){
    for (j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, je+j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	pressure = d1 * T1;
      pGrid->U[ks][je+j][i].d  = d1;
      pGrid->U[ks][je+j][i].M1  = 0.0;
      pGrid->U[ks][je+j][i].M2  = 0.0;
      pGrid->U[ks][je+j][i].M3  = 0.0;
      pGrid->U[ks][je+j][i].E = pressure/(Gamma - 1.0);

      pGrid->U[ks][je+j][i].Edd_11 = pGrid->U[ks][je][i].Edd_11;
      pGrid->U[ks][je+j][i].Edd_22 = pGrid->U[ks][je][i].Edd_22; 
      pGrid->U[ks][je+j][i].Edd_21 = 0.0;
      pGrid->U[ks][je+j][i].Sigma[0] = kappaes *  pGrid->U[ks][je][i].d;
      pGrid->U[ks][je+j][i].Sigma[1] = kappaff * pGrid->U[ks][je][i].d * pGrid->U[ks][je][i].d * pow(pressure/pGrid->U[ks][je][i].d,-3.5);
      pGrid->U[ks][je+j][i].Sigma[2] = kappaff * pGrid->U[ks][je][i].d * pGrid->U[ks][je][i].d * pow(pressure/pGrid->U[ks][je][i].d,-3.5);
      pGrid->U[ks][je+j][i].Sigma[3] = kappaff * pGrid->U[ks][je][i].d * pGrid->U[ks][je][i].d * pow(pressure/pGrid->U[ks][je][i].d,-3.5);

      pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je+j][i].d * (kappaes + kappaff) * consFr * 
	                         (Ertop - x2) / pGrid->U[ks][je+j][i].Edd_22;
      pGrid->U[ks][je+j][i].Fr1 = 0.0;
      pGrid->U[ks][je+j][i].Fr2 = consFr;
      pGrid->U[ks][je+j][i].Fr3 = 0.0;
    }
  }

  for(i=is-nghost; i<=ie+nghost; i++){
    for (j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, js-j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */

      pGrid->U[ks][js-j][i].d  = d1;
      pGrid->U[ks][js-j][i].M1  = 0.0;
      pGrid->U[ks][js-j][i].M2  = 0.0;
      pGrid->U[ks][js-j][i].M3  = 0.0;      
      pGrid->U[ks][js-j][i].E = pressure/(Gamma - 1.0);

      pGrid->U[ks][js-j][i].Edd_11 = pGrid->U[ks][js][i].Edd_11;
      pGrid->U[ks][js-j][i].Edd_22 = pGrid->U[ks][js][i].Edd_22;
      pGrid->U[ks][js-j][i].Edd_21 = 0.0;

      pGrid->U[ks][js-j][i].Sigma[0] = kappaes * pGrid->U[ks][js][i].d;
      pGrid->U[ks][js-j][i].Sigma[1] = kappaff * pGrid->U[ks][js][i].d * pGrid->U[ks][js][i].d * pow(pressure/pGrid->U[ks][js][i].d,-3.5);
      pGrid->U[ks][js-j][i].Sigma[2] = kappaff * pGrid->U[ks][js][i].d * pGrid->U[ks][js][i].d * pow(pressure/pGrid->U[ks][js][i].d,-3.5);
      pGrid->U[ks][js-j][i].Sigma[3] = kappaff * pGrid->U[ks][js][i].d * pGrid->U[ks][js][i].d * pow(pressure/pGrid->U[ks][js][i].d,-3.5);

      pGrid->U[ks][js-j][i].Er = pGrid->U[ks][js-j][i].d * (kappaes + kappaff) * consFr * 
	                         (Erbottom - x2) / pGrid->U[ks][js-j][i].Edd_22;
      pGrid->U[ks][js-j][i].Fr1 = 0.0;
      pGrid->U[ks][js-j][i].Fr2 = consFr;
      pGrid->U[ks][js-j][i].Fr3 = 0.0;
    }
  }

  Opacity = constopa;
  StaticGravPot = grav_pot;

/* set boundary functions */
  bvals_mhd_fun(pDomain, right_x2, radMHD_inflow);
  bvals_rad_fun(pDomain, right_x2, radMHD_rad_inflow);
  bvals_mhd_fun(pDomain, left_x2, radMHD_inflow2);
  bvals_rad_fun(pDomain, left_x2, radMHD_rad_inflow2);

/* initialization for radiation transfer method */
#ifdef RADIATION_TRANSFER
 
  il = pRG->is-1, iu = pRG->ie+1;
  jl = pRG->js-1, ju = pRG->je+1;
  kl = pRG->ks,   ku = pRG->ke;
  //if (pRG->Nx[1] > 1) { jl -= 1; ju += 1; }
  if (pRG->Nx[2] > 1) { kl -= 1; ku += 1; }

/* Initialize mean intensity */
  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for(i=il; i<=iu; i++) {
	pRG->R[ifr][k][j][i].J = pGrid->U[ks][j+nghost-1][i+nghost-1].Er;
	pRG->R[ifr][k][j][i].H[1] = consFr;
      }


/* ------- Initialize boundary emission ---------------------------------- */
  for(k=kl; k<=ku; k++) {
    for(j=jl+1; j<=ju-1; j++) {
      for(m=0; m<nang; m++) {
	pRG->Ghstl1i[ifr][k][j][0][m] = pGrid->U[ks][j+nghost-1][il+nghost-1].Er;
	if (noct > 2) {
	  pRG->Ghstl1i[ifr][k][j][2][m] = pGrid->U[ks][j+nghost-1][il+nghost-1].Er;
	  if (noct == 8) {
	    pRG->Ghstl1i[ifr][k][j][4][m] = pGrid->U[ks][j+nghost-1][il+nghost-1].Er;
	    pRG->Ghstl1i[ifr][k][j][6][m] = pGrid->U[ks][j+nghost-1][il+nghost-1].Er;
	  }
	}
	pRG->Ghstr1i[ifr][k][j][1][m] = pGrid->U[ks][j+nghost-1][iu+nghost-1].Er;
	if (noct > 2) {
	  pRG->Ghstr1i[ifr][k][j][3][m] = pGrid->U[ks][j+nghost-1][iu+nghost-1].Er;
	  if (noct == 8) {
	    pRG->Ghstr1i[ifr][k][j][5][m] = pGrid->U[ks][j+nghost-1][iu+nghost-1].Er;
	    pRG->Ghstr1i[ifr][k][j][7][m] = pGrid->U[ks][j+nghost-1][iu+nghost-1].Er;
	  }
	}
      }}
      for(m=0; m<nang; m++) {
	pRG->Ghstl1i[ifr][k][jl][0][m] = Ibtm;
	pRG->Ghstr1i[ifr][k][jl][1][m] = Ibtm;
	pRG->Ghstl1i[ifr][k][ju][2][m] = Itop;
	pRG->Ghstr1i[ifr][k][ju][3][m] = Itop;
	if (noct == 8) {
	  pRG->Ghstl1i[ifr][k][jl][4][m] = Ibtm;
	  pRG->Ghstr1i[ifr][k][jl][5][m] = Ibtm;
	  pRG->Ghstl1i[ifr][k][ju][6][m] = Itop;
	  pRG->Ghstr1i[ifr][k][ju][7][m] = Itop;
	}
      }
    if (noct > 2) {

/* Initialize boundary intensity in x2 direction */
      /* for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][k][il][l][m] = Itop; 
	  pRG->Ghstl2i[ifr][k][il][l][m] =  Ibtm;
	  
	}}
      for(l=0; l<noct; l++) {
	for(m=0; m<nang; m++) {
	  pRG->Ghstr2i[ifr][k][iu][l][m] = Itop; 
	  pRG->Ghstl2i[ifr][k][iu][l][m] = Ibtm;
	  }}*/
      for(i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = Ibtm;
	  pRG->Ghstl2i[ifr][k][i][1][m] = Ibtm;
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = Ibtm;
	    pRG->Ghstl2i[ifr][k][i][5][m] = Ibtm;	    
	  }
	}
	for(m=0; m<=nang; m++) {
	  pRG->Ghstr2i[ifr][k][i][2][m] = Itop;
	  pRG->Ghstr2i[ifr][k][i][3][m] = Itop;
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = Itop;
	    pRG->Ghstr2i[ifr][k][i][7][m] = Itop;
	  }
	}
      }
    }
  }

/* enrol user-defined  boundary functions */
  bvals_rad_trans_fun(pDomain, right_x2, const_J_ox2);
  bvals_rad_trans_fun(pDomain, left_x2, const_H_ix2);

 
/* enroll radiation specification functions */
get_thermal_source = Thermal_B;
get_thermal_fraction = const_eps;
get_total_opacity = transfer_opacity;

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
 * Userwork_in_formal_solution  - problem specific work in formal solution loop
 * Userwork_after_formal_solution  - problem specific work after formal solution
 *----------------------------------------------------------------------------*/

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ox1(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ix2(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_inflowi2;
} 
void bvals_mat_fun_ox2(VMatFun_t *Mat_BCFun)
{
	*Mat_BCFun = radMHD_Mat_inflowo2;
} 
void bvals_mat_fun_ix3(VMatFun_t *Mat_BCFun)
{

} 
void bvals_mat_fun_ox3(VMatFun_t *Mat_BCFun)
{

} 
#endif

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
	


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID 
	fwrite(&Ncycle,sizeof(Real),1,fp);
	fwrite(&TOL,sizeof(Real),1,fp);
#endif
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
	fread(&kappaff,sizeof(Real),1,fp);
#endif
	


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
#ifdef MATRIX_MULTIGRID 
	fread(&Ncycle,sizeof(Real),1,fp);
	fread(&TOL,sizeof(Real),1,fp);
#endif
#endif



  return;
}

/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
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


#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)

void constopa(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY])
{

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


static Real grav_pot(const Real x1, const Real x2, const Real x3)
{

  return grav*x2;
}


/* right hand side of x2 direction */

void radMHD_inflow(GridS *pGrid)
{

  int i, je,j, ju, jl, js;
  int ks, is,ie;
  Real u0, opac;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  u0 = 0.0;
	
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;
 
  for(i=is-nghost; i<=ie+nghost; i++) {
    /* this looks wrong: */
/*    opac = pGrid->U[ks][je][i].d * pGrid->U[ks][je+j][i].d * pow(pressure/pGrid->U[ks][je+j][i].d,-3.5);
*/
    for (j=1;  j<=nghost;  j++) {	
#ifdef RADIATION_MHD
      pGrid->B1i[ks][je+j][i] = sqrt(fabs(u0));
      pGrid->U[ks][je+j][i].B1c =  sqrt(fabs(u0));
#endif
      pGrid->U[ks][je+j][i].M1 = pGrid->U[ks][je-j+1][i].M1;
      pGrid->U[ks][je+j][i].M2 = -0.0 * pGrid->U[ks][je-j+1][i].M2;
      pGrid->U[ks][je+j][i].M3 = 0.0;
      /*pGrid->U[ks][je+j][i].d  = pGrid->U[ks][je+j-1][i].d;
      pGrid->U[ks][je+j][i].Sigma[0] = kappaes *  pGrid->U[ks][je+j][i].d;
      pGrid->U[ks][je+j][i].Sigma[1] = kappaff * opac;
      pGrid->U[ks][je+j][i].Sigma[2] = kappaff * opac;
      pGrid->U[ks][je+j][i].Sigma[3] = kappaff * opac;
      
      pGrid->U[ks][je+j][i].E  = pressure/(Gamma_1)+0.5*(SQR(pGrid->U[ks][je+j][i].M1) + 
      SQR(pGrid->U[ks][je+j][i].M2)) / pGrid->U[ks][je+j][i].d;*/

      pGrid->U[ks][je+j][i].d  = d1;
      pGrid->U[ks][je+j][i].Sigma[0] = pGrid->U[ks][je+j-1][i].Sigma[0];
      pGrid->U[ks][je+j][i].Sigma[1] = pGrid->U[ks][je+j-1][i].Sigma[1];
      pGrid->U[ks][je+j][i].Sigma[2] = pGrid->U[ks][je+j-1][i].Sigma[2];  
      pGrid->U[ks][je+j][i].Sigma[3] = pGrid->U[ks][je+j-1][i].Sigma[3];
      pGrid->U[ks][je+j][i].E = d1 * T1/(Gamma_1)+0.5*(SQR(pGrid->U[ks][je+j][i].M1) + 
				 SQR(pGrid->U[ks][je+j][i].M2)) / pGrid->U[ks][je+j][i].d;
    }
  }
  return;
}

void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, je,j, ju, jl, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	Real x1, x2, x3;

	Real Fr0x, Fr0y;
	Real velocity_x, velocity_y;
	Real velocity_x1, velocity_y1;
	Real Sigma_t, Sigma_t1;
	Real Eratio, reducefactor, tau;

	Real dz = pGrid->dx3;

	
	for(i=is-nghost; i<=ie+nghost; i++){
	    for(j=1;  j<=nghost;  j++) {
		cc_pos(pGrid, i, je+j,ks, &x1, &x2, &x3);

		pGrid->U[ks][je+j][i].Edd_11 = pGrid->U[ks][je][i].Edd_11;
		pGrid->U[ks][je+j][i].Edd_22 = pGrid->U[ks][je][i].Edd_22;
		pGrid->U[ks][je+j][i].Edd_21 = pGrid->U[ks][je][i].Edd_21;

		matrix_alpha(0.0, pGrid->U[ks][je+j][i].Sigma, pGrid->dt, pGrid->U[ks][je+j][i].Edd_33, 0.0, &reducefactor, 0, dz);

		Sigma_t = 0.5 * (pGrid->U[ks][je+j][i].Sigma[0] + pGrid->U[ks][je+j][i].Sigma[1] + pGrid->U[ks][je+j-1][i].Sigma[0] + pGrid->U[ks][je+j-1][i].Sigma[1]);

		Sigma_t1 = 0.5 * (pGrid->U[ks][je+j-1][i].Sigma[0] + pGrid->U[ks][je+j-1][i].Sigma[1] + pGrid->U[ks][je+j-2][i].Sigma[0] + pGrid->U[ks][je+j-2][i].Sigma[1]);

		velocity_x = pGrid->U[ks][je+j-1][i].M1 / pGrid->U[ks][je+j-1][i].d;
		velocity_y = pGrid->U[ks][je+j-1][i].M2 / pGrid->U[ks][je+j-1][i].d;

		Fr0x = pGrid->U[ks][je+j-1][i].Fr1 - ((1.0 + pGrid->U[ks][je+j-1][i].Edd_11) * velocity_x + pGrid->U[ks][je+j-1][i].Edd_21 * velocity_y)* pGrid->U[ks][je+j-1][i].Er / Crat;

		Fr0y = pGrid->U[ks][je+j-1][i].Fr2 - ((1.0 + pGrid->U[ks][je+j-1][i].Edd_22) * velocity_y + pGrid->U[ks][je+j-1][i].Edd_21 * velocity_x)* pGrid->U[ks][je+j-1][i].Er / Crat;

		Fr0y = consFr;

		velocity_x1 = pGrid->U[ks][je+j][i].M1 / pGrid->U[ks][je+j][i].d;
		velocity_y1 = pGrid->U[ks][je+j][i].M2 / pGrid->U[ks][je+j][i].d;

/*		pGrid->U[ks][je+j][i].Er  = pGrid->U[ks][je][i].d * (kappaes + kappaff) * consFr * (Ertop - x2) / pGrid->U[ks][je][i].Edd_22;
*/
		/* use diffusion boundary */
		if(Fr0y > 0.0)
			pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je+j-1][i].Er - pGrid->dx2 * Sigma_t1 * Fr0y / pGrid->U[ks][je+j][i].Edd_22;
		
		if(pGrid->U[ks][je+j][i].Er < 0.0)
			pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je][i].Er;


/*		pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je+j][i].d * (kappaes + kappaff) * consFr * (Ertop - x2) / pGrid->U[ks][je][i].Edd_22;
*/
/*		Eratio = pGrid->U[ks][je+j][i].Edd_22 + 0.5 * pGrid->dx2 * reducefactor * Sigma_t;

		pGrid->U[ks][je+j][i].Er = (pGrid->U[ks][je+j-1][i].Edd_22 * pGrid->U[ks][je+j-1][i].Er - 0.5 * pGrid->dx2 * Sigma_t * Fr0y ) / Eratio;

		if((pGrid->U[ks][je+j][i].Er > pGrid->U[ks][je+j-1][i].Er) || (pGrid->U[ks][je+j][i].Er < 0.0)){
			Eratio = pGrid->U[ks][je+j][i].Edd_22 + pGrid->dx2 * reducefactor * Sigma_t;
			pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je+j-1][i].Edd_22 * pGrid->U[ks][je+j-1][i].Er / Eratio;
		}

		Fr0y = reducefactor * pGrid->U[ks][je+j][i].Er;
*/
		pGrid->U[ks][je+j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ks][je+j][i].Edd_11) * velocity_x1 + pGrid->U[ks][je+j][i].Edd_21 * velocity_y1)* pGrid->U[ks][je+j][i].Er / Crat;

                pGrid->U[ks][je+j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ks][je+j][i].Edd_22) * velocity_y1 + pGrid->U[ks][je+j][i].Edd_21 * velocity_x1)* pGrid->U[ks][je+j][i].Er / Crat;


		pGrid->U[ks][je+j][i].Fr3 = 0.0;

    }
	}

  
}

void radMHD_rad_inflow_old(GridS *pGrid)
{
  int i, je,j, ju, jl, js;
  int ks, is,ie;
  Real x1, x2, x3, xt;
  Real Fr0x, Fr0y;
  Real velocity_x, velocity_y;
  Real velocity_x1, velocity_y1;
  Real Sigma_t;
  
  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;
  
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  cc_pos(pGrid, is, je,ks, &x1, &xt, &x3);
  for(i=is-nghost; i<=ie+nghost; i++){
    for(j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, je+j,ks, &x1, &x2, &x3);
      
      pGrid->U[ks][je+j][i].Edd_11 = pGrid->U[ks][je][i].Edd_11;
      pGrid->U[ks][je+j][i].Edd_22 = pGrid->U[ks][je][i].Edd_22;
      pGrid->U[ks][je+j][i].Edd_21 = pGrid->U[ks][je][i].Edd_21;
      
      /* Sigma_t = 0.5 * (pGrid->U[ks][je+j  ][i].Sigma[0] + pGrid->U[ks][je+j  ][i].Sigma[1] + 
                       pGrid->U[ks][je+j-1][i].Sigma[0] + pGrid->U[ks][je+j-1][i].Sigma[1]);

      velocity_x = pGrid->U[ks][je+j-1][i].M1 / pGrid->U[ks][je+j-1][i].d;
      velocity_y = pGrid->U[ks][je+j-1][i].M2 / pGrid->U[ks][je+j-1][i].d;

      Fr0x = pGrid->U[ks][je+j-1][i].Fr1 - ((1.0 + pGrid->U[ks][je+j-1][i].Edd_11) * velocity_x + 
             pGrid->U[ks][je+j-1][i].Edd_21 * velocity_y)* pGrid->U[ks][j][i].Er / Crat;
      Fr0y = pGrid->U[ks][je+j-1][i].Fr2 - ((1.0 + pGrid->U[ks][je+j-1][i].Edd_22) * velocity_y + 
      pGrid->U[ks][je+j-1][i].Edd_21 * velocity_x)* pGrid->U[ks][j][i].Er / Crat;*/
     
      /* use diffusion boundary */
      /*if(Fr0y > 0.0)
	pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je+j-1][i].Er - pGrid->dx2 * Sigma_t * Fr0y / pGrid->U[ks][je+j][i].Edd_22;		
      if(pGrid->U[ks][je+j][i].Er < 0.0)
	pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je][i].Er;*/


      pGrid->U[ks][je+j][i].Er = pGrid->U[ks][je][i].Er + pGrid->U[ks][je+j][i].d * (kappaes + kappaff) * 
	                         consFr * (xt - x2) / pGrid->U[ks][je][i].Edd_22;

      pGrid->U[ks][je+j][i].Fr1 = pGrid->U[ks][je][i].Fr1;
      pGrid->U[ks][je+j][i].Fr2 = consFr;
      pGrid->U[ks][je+j][i].Fr3 = 0.0;
      
    }
  }

  return;
}


void radMHD_inflow2(GridS *pGrid)
{
  int i, je,j, ju, jl, js;
  int ks, is,ie;
  Real u0, opac;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;
  
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  u0 = 0.0;
 
  for(i=is-nghost; i<=ie+nghost; i++) {
    opac = SQR(pGrid->U[ks][js][i].d) * pow(pressure/pGrid->U[ks][js][i].d,-3.5);
    for (j=1;  j<=nghost;  j++) {	
#ifdef RADIATION_MHD
      pGrid->B1i[ks][js-j][i] = sqrt(fabs(u0));
      pGrid->U[ks][js-j][i].B1c =  sqrt(fabs(u0));
#endif
		
      pGrid->U[ks][js-j][i].M1 = pGrid->U[ks][js+j-1][i].M1;
      pGrid->U[ks][js-j][i].M2 = -0.0 * pGrid->U[ks][js+j-1][i].M2;
      pGrid->U[ks][js-j][i].M3 = 0.0;
      pGrid->U[ks][js-j][i].d  = d1;
      pGrid->U[ks][js-j][i].Sigma[0] = kappaes *  pGrid->U[ks][js-j][i].d;
      pGrid->U[ks][js-j][i].Sigma[1] = kappaff * opac;
      pGrid->U[ks][js-j][i].Sigma[2] = kappaff * opac;
      pGrid->U[ks][js-j][i].Sigma[3] = kappaff * opac;      
      pGrid->U[ks][js-j][i].E  = d1 * T1/(Gamma_1)+0.5*(SQR(pGrid->U[ks][js-j][i].M1) + 
				 SQR(pGrid->U[ks][js-j][i].M2)) / pGrid->U[ks][js-j][i].d;
/*pGrid->U[ks][js+j-1][i].E;
*/ 
/*
*/
    }
  }
  return;
}

void radMHD_FeqH_ix2(GridS *pGrid)
{
  int i, j, k;
  int is,ie,js,ks, ke;
  Real x1, x2, x3, Fy;

  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  for(i=is-nghost; i<=ie+nghost; i++){
    for(j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, js-j,ks, &x1, &x2, &x3);

      pGrid->U[ks][js-j][i].Edd_11 = pGrid->U[ks][js][i].Edd_11;
      pGrid->U[ks][js-j][i].Edd_22 = pGrid->U[ks][js][i].Edd_22;
      pGrid->U[ks][js-j][i].Edd_21 = pGrid->U[ks][js][i].Edd_21;
      
      pGrid->U[ks][js-j][i].Er  = pGrid->U[ks][js-j][i].d * (kappaes + kappaff) * Fy * (Erbottom - x2) / 
	                          pGrid->U[ks][js][i].Edd_22;
      pGrid->U[ks][js-j][i].Fr1 = pGrid->U[ks][js][i].Fr1;
      pGrid->U[ks][js-j][i].Fr2 = Fy;      
      pGrid->U[ks][js-j][i].Fr3 = 0.0;
    }
  }
  return;
}

void radMHD_FeqH_ox2(GridS *pGrid)
{
  int i, j, k;
  int is,ie,je,ks, ke;
  Real x1, x2, x3, Fy;

  je = pGrid->je;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  for(i=is-nghost; i<=ie+nghost; i++){
    for(j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, je+j, ks, &x1, &x2, &x3);

      pGrid->U[ks][je+j][i].Edd_11 = pGrid->U[ks][je][i].Edd_11;
      pGrid->U[ks][je+j][i].Edd_22 = pGrid->U[ks][je][i].Edd_22;
      pGrid->U[ks][je+j][i].Edd_21 = pGrid->U[ks][je][i].Edd_21;
      
      pGrid->U[ks][je+j][i].Er  = pGrid->U[ks][je+j][i].d * (kappaes + kappaff) * Fy * (Erbottom - x2) / 
	                          pGrid->U[ks][je][i].Edd_22;
      pGrid->U[ks][je+j][i].Fr1 = pGrid->U[ks][je][i].Fr1;
      pGrid->U[ks][je+j][i].Fr2 = Fy;      
      pGrid->U[ks][je+j][i].Fr3 = 0.0;
    }
  }
  return;
}

void radMHD_rad_inflow2(GridS *pGrid)
{
  	int i, je,j, ju, jl, js;
	int ks, is,ie;
	je = pGrid->je;
	js = pGrid->js;
  	ks = pGrid->ks;
	is = pGrid->is;
	ie = pGrid->ie;

	ju = pGrid->je + nghost;
  	jl = pGrid->js - nghost;

	Real x1, x2, x3;
	Real Sigma_t1;

	Real velocity_x, velocity_y, velocity_x1, velocity_y1;
	Real Fr0x, Fr0y;
	
	for(i=is-nghost; i<=ie+nghost; i++){
	    for(j=1;  j<=nghost;  j++) {
		cc_pos(pGrid, i, js-j,ks, &x1, &x2, &x3);

		pGrid->U[ks][js-j][i].Edd_11 = pGrid->U[ks][js][i].Edd_11;
		pGrid->U[ks][js-j][i].Edd_22 = pGrid->U[ks][js][i].Edd_22;
		pGrid->U[ks][js-j][i].Edd_21 = pGrid->U[ks][js][i].Edd_21;

		Sigma_t1 = 0.5 * (pGrid->U[ks][js-j+1][i].Sigma[0] + pGrid->U[ks][js-j+1][i].Sigma[1] + pGrid->U[ks][js-j+2][i].Sigma[0] + pGrid->U[ks][js-j+2][i].Sigma[1]);
	
		
		velocity_x = pGrid->U[ks][js-j+1][i].M1 / pGrid->U[ks][js-j+1][i].d;
		velocity_y = pGrid->U[ks][js-j+1][i].M2 / pGrid->U[ks][js-j+1][i].d;

		Fr0x = pGrid->U[ks][js-j+1][i].Fr1 - ((1.0 + pGrid->U[ks][js-j+1][i].Edd_11) * velocity_x + pGrid->U[ks][js-j+1][i].Edd_21 * velocity_y)* pGrid->U[ks][js-j+1][i].Er / Crat;

		Fr0y = consFr;

		velocity_x1 = pGrid->U[ks][js-j][i].M1 / pGrid->U[ks][js-j][i].d;
		velocity_y1 = pGrid->U[ks][js-j][i].M2 / pGrid->U[ks][js-j][i].d;

/*		pGrid->U[ks][je+j][i].Er  = pGrid->U[ks][je][i].d * (kappaes + kappaff) * consFr * (Ertop - x2) / pGrid->U[ks][je][i].Edd_22;
*/
		/* use diffusion boundary */
		if(Fr0y > 0.0)
			pGrid->U[ks][js-j][i].Er = pGrid->U[ks][js-j+1][i].Er + pGrid->dx2 * Sigma_t1 * Fr0y / pGrid->U[ks][js-j][i].Edd_22;
		
		pGrid->U[ks][js-j][i].Fr1 = Fr0x + ((1.0 + pGrid->U[ks][js-j][i].Edd_11) * velocity_x1 + pGrid->U[ks][js-j][i].Edd_21 * velocity_y1)* pGrid->U[ks][js-j][i].Er / Crat;

                pGrid->U[ks][js-j][i].Fr2 = Fr0y + ((1.0 + pGrid->U[ks][js-j][i].Edd_22) * velocity_y1 + pGrid->U[ks][js-j][i].Edd_21 * velocity_x1)* pGrid->U[ks][js-j][i].Er / Crat;


/*
		pGrid->U[ks][js-j][i].Er  = pGrid->U[ks][js-j][i].d * (kappaes + kappaff) * consFr * (Erbottom - x2) / pGrid->U[ks][js][i].Edd_22;
      		pGrid->U[ks][js-j][i].Fr1 = pGrid->U[ks][js][i].Fr1;
		pGrid->U[ks][js-j][i].Fr2 = consFr;
*/
		pGrid->U[ks][js-j][i].Fr3 = 0.0;

    }
	}

  
}



void radMHD_rad_inflow2_old(GridS *pGrid)
{
  int i, je,j, ju, jl, js;
  int ks, is,ie;
  Real x1, x2, x3, xb;
  
  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;
  
  ju = pGrid->je + nghost;
  jl = pGrid->js - nghost;

  cc_pos(pGrid, is, js, ks, &x1, &xb, &x3);
  for(i=is-nghost; i<=ie+nghost; i++){
    for(j=1;  j<=nghost;  j++) {
      cc_pos(pGrid, i, js-j,ks, &x1, &x2, &x3);
      
      pGrid->U[ks][js-j][i].Edd_11 = pGrid->U[ks][js][i].Edd_11;
      pGrid->U[ks][js-j][i].Edd_22 = pGrid->U[ks][js][i].Edd_22;
      pGrid->U[ks][js-j][i].Edd_21 = pGrid->U[ks][js][i].Edd_21;
      
      pGrid->U[ks][js-j][i].Er  = pGrid->U[ks][js][i].Er + pGrid->U[ks][js-j][i].d * (kappaes + kappaff) * 
	                          consFr * (xb - x2) / pGrid->U[ks][js][i].Edd_22;
      pGrid->U[ks][js-j][i].Fr1 = pGrid->U[ks][js][i].Fr1;
      pGrid->U[ks][js-j][i].Fr2 = consFr;      
      pGrid->U[ks][js-j][i].Fr3 = 0.0;
    }
  }
  return;
}



void radMHD_Mat_inflowi2(MatrixS *pMat)
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

	
	for(j=1; j<=Matghost; j++){
	    for (i=is-Matghost;  i<=ie+Matghost;  i++) {
			

	      	pMat->Ugas[ks][js-j][i].Edd_11 = pMat->Ugas[ks][js][i].Edd_11;
		pMat->Ugas[ks][js-j][i].Edd_22 = pMat->Ugas[ks][js][i].Edd_22;
		pMat->Ugas[ks][js-j][i].Edd_21 = pMat->Ugas[ks][js][i].Edd_21;
	
		pMat->Ugas[ks][js-j][i].V1 = pMat->Ugas[ks][js+j-1][i].V1;
                pMat->Ugas[ks][js-j][i].V2 = -pMat->Ugas[ks][js+j-1][i].V2;
		pMat->Ugas[ks][js-j][i].T4 = pMat->Ugas[ks][js][i].T4;
		pMat->Ugas[ks][js-j][i].Sigma[0] = pMat->Ugas[ks][js][i].Sigma[0];
		pMat->Ugas[ks][js-j][i].Sigma[1] = pMat->Ugas[ks][js][i].Sigma[1];
		pMat->Ugas[ks][js-j][i].Sigma[2] = pMat->Ugas[ks][js][i].Sigma[2];
		pMat->Ugas[ks][js-j][i].Sigma[3] = pMat->Ugas[ks][js][i].Sigma[3];

		pMat->U[ks][js-j][i].Er  = 0.0;
		pMat->U[ks][js-j][i].Fr1 = 0.0;
		pMat->U[ks][js-j][i].Fr2 = 0.0;

		pMat->U[ks][js-j][i].Fr3 = 0.0;		

    }
	}

  
}



void radMHD_Mat_inflowo2(MatrixS *pMat)
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

	
	for(j=1; j<=Matghost; j++){
	    for (i=is-Matghost;  i<=ie+Matghost;  i++) {
			

	      	pMat->Ugas[ks][je+j][i].Edd_11 = pMat->Ugas[ks][je][i].Edd_11;
		pMat->Ugas[ks][je+j][i].Edd_22 = pMat->Ugas[ks][je][i].Edd_22;
		pMat->Ugas[ks][je+j][i].Edd_21 = pMat->Ugas[ks][je][i].Edd_21;
	
		pMat->Ugas[ks][je+j][i].V1 = pMat->Ugas[ks][je][i].V1;
                pMat->Ugas[ks][je+j][i].V2 = -pMat->Ugas[ks][je-j+1][i].V2;
                pMat->Ugas[ks][je+j][i].T4 = pMat->Ugas[ks][je][i].T4;
                pMat->Ugas[ks][je+j][i].Sigma[0] = pMat->Ugas[ks][je][i].Sigma[0];
                pMat->Ugas[ks][je+j][i].Sigma[1] = pMat->Ugas[ks][je][i].Sigma[1];
                pMat->Ugas[ks][je+j][i].Sigma[2] = pMat->Ugas[ks][je][i].Sigma[2];
                pMat->Ugas[ks][je+j][i].Sigma[3] = pMat->Ugas[ks][je][i].Sigma[3];

		pMat->U[ks][je+j][i].Er  = 0.0;
		pMat->U[ks][je+j][i].Fr1 = 0.0;
		pMat->U[ks][je+j][i].Fr2 = 0.0;

		pMat->U[ks][je+j][i].Fr3 = 0.0;		

    }
	}

  
}

/* Function for transfer module */
#ifdef RADIATION_TRANSFER

static Real Thermal_B(const GridS *pG, const int ifr, const int i, 
		      const int j, const int k)
{
  
  /* no thermal emission -- scattering dominated */
  return 0.0;
}


static Real const_eps(const GridS *pG, const int ifr, const int i, 
		      const int j, const int k)
{
	Real eps;
/* eps =0 for pure scattering opacity */

	eps = pG->U[k][j][i].Sigma[1] / (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[1]);

  return eps;
  
}

static Real transfer_opacity(const GridS *pG, const int ifr, const int i, 
			     const int j, const int k)
{
	
  return (pG->U[k][j][i].Sigma[0] + pG->U[k][j][i].Sigma[1]);
  
}

void const_H_ix2(RadGridS *pRG, int ifr)
{
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js-1;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, k, l, m, n;
  Real I0, Jm, H, Hm, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

/* update Ghstl2i using l2imu */
  for (k=kl; k<=ku; k++) {
    for (i=il; i<=iu; i++) {
      Hm = 0.0;
      Jm = 0.0;
      for (m=0; m<nang; m++) {
	Hm += pRG->l2imu[ifr][k][i][2][m] * pRG->mu[2][m][1] * pRG->wmu[m];
	Hm += pRG->l2imu[ifr][k][i][3][m] * pRG->mu[3][m][1] * pRG->wmu[m];
	Jm += pRG->l2imu[ifr][k][i][2][m] * pRG->wmu[m];
	Jm += pRG->l2imu[ifr][k][i][3][m] * pRG->wmu[m];	  
	if (noct == 8) {
	  Hm += pRG->l2imu[ifr][k][i][6][m] * pRG->mu[6][m][1] * pRG->wmu[m];
	  Hm += pRG->l2imu[ifr][k][i][7][m] * pRG->mu[7][m][1] * pRG->wmu[m];
	  Jm += pRG->l2imu[ifr][k][i][6][m] * pRG->wmu[m];
	  Jm += pRG->l2imu[ifr][k][i][7][m] * pRG->wmu[m];
	}
	}	
      H = consFr;
      I0 = (H - Hm) / gamma;
      if (I0 < 0.0) I0 = 0.0; 
      pRG->R[ifr][k][jl][i].J = 0.5 * I0 + Jm;
      for (m=0; m<nang; m++) {
	pRG->Ghstl2i[ifr][k][i][0][m] = I0;
	pRG->Ghstl2i[ifr][k][i][1][m] = I0;
	if (noct == 8) {
	  pRG->Ghstl2i[ifr][k][i][4][m] = I0;
	  pRG->Ghstl2i[ifr][k][i][5][m] = I0;
	}
      }
    }
/* update Ghstr1i and Ghstl1i so corner intensities are correct w/ periodic bcs */
    for(l=0; l<noct; l++) {
      for (m=0; m<nang; m++) {
	pRG->Ghstl2i[ifr][k][il-1][l][m] = pRG->Ghstl2i[ifr][k][iu][l][m];
	pRG->Ghstl2i[ifr][k][iu+1][l][m] = pRG->Ghstl2i[ifr][k][il][l][m];
      }}
    for (m=0; m<nang; m++) {
      pRG->Ghstl1i[ifr][k][jl][0][m] = pRG->Ghstl2i[ifr][k][iu][0][m];
      pRG->Ghstr1i[ifr][k][jl][1][m] = pRG->Ghstl2i[ifr][k][il][1][m];
    }
  }

  return;
}

void const_J_ox2(RadGridS *pRG, int ifr)
{
  int il = pRG->is-1, iu = pRG->ie+1;
  int ju = pRG->je+1;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, k, l, m, n;
  int ig, jg, kg, io, ko;
  Real I0, Jp, J, Hp, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][1] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;
  jg = ju + nghost - 1;
  io = nghost - 1;
  if (pRG->pG->Nx[2] > 1) {
    ko = nghost - 1;
  } else ko = 0;

/* update Ghstr2i using r2imu */
  for (k=kl; k<=ku; k++) {
    kg = k + ko;
    for (i=il; i<=iu; i++) {
      ig = i + io;
      Hp = 0.0;
      Jp = 0.0;
      for (m=0; m<nang; m++) {
	Hp += pRG->r2imu[ifr][k][i][0][m] * pRG->mu[0][m][1] * pRG->wmu[m];
	Hp += pRG->r2imu[ifr][k][i][1][m] * pRG->mu[1][m][1] * pRG->wmu[m];
	Jp += pRG->r2imu[ifr][k][i][0][m] * pRG->wmu[m];
	Jp += pRG->r2imu[ifr][k][i][1][m] * pRG->wmu[m];
	if (noct == 8) {
	  Hp += pRG->r2imu[ifr][k][i][4][m] * pRG->mu[4][m][1] * pRG->wmu[m];
	  Hp += pRG->r2imu[ifr][k][i][5][m] * pRG->mu[5][m][1] * pRG->wmu[m];
	  Jp += pRG->r2imu[ifr][k][i][4][m] * pRG->wmu[m];
	  Jp += pRG->r2imu[ifr][k][i][5][m] * pRG->wmu[m];
	}
      }
      //J = pRG->R[ifr][k][ju][i].J;
      J = pRG->pG->U[kg][jg][ig].Er;
      I0 = 2.0 * (J - Jp); 
      if (I0 < 0.0) I0 = 0.0; 
      for (m=0; m<nang; m++) {
	pRG->Ghstr2i[ifr][k][i][2][m] = I0;
	pRG->Ghstr2i[ifr][k][i][3][m] = I0;
	if (noct == 8) {
	  pRG->Ghstr2i[ifr][k][i][6][m] = I0;
	  pRG->Ghstr2i[ifr][k][i][7][m] = I0;
	}
      }
    }}

  return;
}

void Userwork_in_formal_solution(DomainS *pD)
{
  /*GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i, j, k, DIM;
  int is, ie, js, je, ks, ke;
  int ri, rj, rk;
  int ioff, joff, koff;
  int ifr = 0;
 
  DIM = 0;
  is = pG->is;
  ie = pG->ie;
  js = pG->js;
  je = pG->je;
  ks = pG->ks;
  ke = pG->ke;

 for (i=0; i<3; i++) if(pG->Nx[i] > 1) ++DIM;
  ioff = 1 - nghost;
  if(DIM > 1) { 
    joff = 1 - nghost;
    if (DIM == 3)
      koff = 1 - nghost;
    else
      koff = 0;
  } else
    joff = 0;
  
  for(k=ks; k<=ke; k++) {
    rk = k + koff;
    for(i=is; i<=ie; i++) {
      ri = i + ioff;      
      H[k][i][0] = pRG->R[ifr][rk][js+joff][ri].H[1];
      H[k][i][1] = pRG->R[ifr][rk][je+joff][ri].H[1];
    }}
  */
  return;
}

void Userwork_after_formal_solution(DomainS *pD)
{

  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k, ifr;
  int il = pG->is, iu = pG->ie;
  int jl = pG->js, ju = pG->je;
  int kl = pG->ks, ku = pG->ke;
  int nDim;
  int ir,jr,kr,ioff,joff,koff;
  Real J, H, K;
  Real x1, x2, x3;
#ifdef MPI_PARALLEL
  Real Kloc, Hloc, Jloc;
#endif
  static int fstflag = 1;

  if (fstflag != 1) return;

  ioff = nghost - 1;
  nDim = 1;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1;
    nDim = 2;
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1;
    nDim = 3;
  } else koff = 0;

/* Compute new value of Erad at surface */
 
#ifdef MPI_PARALLEL
  Kloc = 0.0; Hloc = 0.0; 
  Jloc =0.0;
  for (i = pRG->is; i <= pRG->ie; i++) {
    if(pG->rx2_id == -1) {
      Kloc += pRG->R[0][pRG->ks][pRG->je][i].K[2];
      Hloc += pRG->R[0][pRG->ks][pRG->je][i].H[1];
      Jloc += pRG->R[0][pRG->ks][pRG->je][i].J;
    }
  }
  MPI_Allreduce(&Kloc, &K, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Hloc, &H, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Jloc, &J, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  K = 0.0; H = 0.0;
  for (i = pRG->is; i <= pRG->ie; i++) {
    K += pRG->R[0][pRG->ks][pRG->je][i].K[2];
    H += pRG->R[0][pRG->ks][pRG->je][i].H[1];
  }
#endif

  Ertop = ytop + K / ((kappaes + kappaff) * d1 * H);
  Ermid = yintt + (Ertop - yintt) * d1 / d2;
  Erbottom = (yintt - yintb) * d2 / d1 + Ertop + yintb - yintt;

/* Update Erad and f22 as in intialization but using results of radiative transfer scheme */
  for (k=kl; k<=ku; k++) {
    kr = k - koff;
    for (j=jl; j<=ju; j++) {
      jr = j - joff;
      for (i=il; i<=iu; i++) {
	ir = i - ioff;
	J = pRG->R[0][kr][jr][ir].J;
	if((fabs(J) > TINY_NUMBER) && (J > 0.0)) {
	  pG->U[k][j][i].Edd_22 = pRG->R[0][kr][jr][ir].K[2] / J;
	}
	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	if(x2 < yintb) {
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Erbottom - x2) / pG->U[k][j][i].Edd_22;
	} else if(x2 < yintt) {
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Ermid - x2) / pG->U[k][j][i].Edd_22;
	} else{
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Ertop - x2) / pG->U[k][j][i].Edd_22;
	}
      }}}
  fstflag = 0;

  return;
}

#endif /* End radiation transfer */


