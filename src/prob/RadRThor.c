#include "copyright.h"
/*==============================================================================
 * FILE: RadRThor.c
 *
 * PURPOSE:  Problem generator for 2D Rayleigh-Taylor simulation of a
 *   thin shell.  This problem generator assumes gravity points in the x
 *   direction so that it can be used with the ray tracing algorithm.
 *   Adapted by S.W. Davis from a similar problem generator written by
 *   Y.-F. Jiang.
 *
 * configure:
 *  --with-problem=RadRThor --enable-radiation-transfer --with-gas=rad_hydro 
 *  --enable-ray-tracing
 *============================================================================*/

#include <float.h>
#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"


static Real kappaes = 1.0;
static Real kappaff = 0.0;
static Real grav = 16.815;
static Real consFr;
static Real Lx, Ly, Lz;
static Real d1, d2, T1;
static Real pressure;
static Real Ertop, Erbottom, Ermid;
static Real ytop, ybtm, yintt, yintb;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * ran2()          - random number generator from NR
 * UnstratifiedDisk() - tidal potential in 3D shearing box
 * expr_dV2()       - computes delta(Vy)
 * hst_*            - new history variables
 *============================================================================*/

static Real grav_pot(const Real x1, const Real x2, const Real x3);
void constopa(const Real rho, const Real T, Real Sigma[NOPACITY], Real dSigma[2*NOPACITY]);
static double ran2(long int *idum);

void radMHD_inflow(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
void radMHD_inflow2(GridS *pGrid);
void radMHD_rad_inflow2(GridS *pGrid);

void radMHD_Mat_inflowi2(MatrixS *pMat);
void radMHD_Mat_inflowo2(MatrixS *pMat);


#ifdef RADIATION_TRANSFER
void const_H_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife);
void const_J_ox1(GridS *pG, RadGridS *pRG, int ifs, int ife);

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

#ifdef 	RAY_TRACING
  Itop = 0.0;
#else
  Itop = Irat * consFr;
#endif
  ytop = pDomain->RootMaxX[0];
  ybtm = pDomain->RootMinX[0];	
  
  Ertop = ytop + (0.57735 + 0.33333 * Itop / consFr) / ((kappaes + kappaff) * d1);
  Ermid = yintt + (Ertop - yintt) * d1 / d2;
  Erbottom = (yintt - yintb) * d2 / d1 + Ertop + yintb - yintt;

#ifdef RAY_TRACING 
  Ibtm = 0.0;
#else
  Ibtm = d1 * (kappaes + kappaff) * consFr * (Erbottom - ybtm) * 3.0;
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pGrid->Disp[0];
  jxs = pGrid->Disp[1];
  kxs = pGrid->Disp[2];
  iseed = -1 - (jxs + pDomain->Nx[1]*(ixs + pDomain->Nx[0]*kxs));
  printf("%d %d %d %d %d\n",myID_Comm_world,ixs,jxs,kxs,iseed);
/* Initialize boxsize */
  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  Lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];
 

/* Rescale amp to sound speed for ipert 2,3 */
  
  tempEdd11 = 1.0/3.0;
  tempEdd22 = 1.0/3.0;
  for (k=ks; k<=ke; k++) {
    for (j=js-nghost; j<=je+nghost; j++) {
      for (i=is; i<=ie; i++) {
	cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	pressure = d1 * T1;

	pGrid->U[k][j][i].d  = d1;
	pGrid->U[k][j][i].M1 = 0.0;
	pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
	if ((x1 < yintt) && (x1 > yintb)) {
	  pGrid->U[k][j][i].d = d2;
          pGrid->U[k][j][i].M1 *= d2/d1;
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

	pGrid->U[k][j][i].Fr1 = consFr;
	pGrid->U[k][j][i].Fr2 = 0.0;
	pGrid->U[k][j][i].Fr3 = 0.0;	
			
	if(x1 < yintb) {	  
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
                                 (Erbottom - x1) / pGrid->U[k][j][i].Edd_11;
	} else if(x2 < yintt) {
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                         (Ermid - x1) / pGrid->U[k][j][i].Edd_11;
	} else {
	  pGrid->U[k][j][i].Er = pGrid->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                         (Ertop - x1) / pGrid->U[k][j][i].Edd_11;	
	}

	/* perturb density */

	if(fabs(x1) < 0.25 * Lx){

	  if (iprob == 1) {
	    pGrid->U[k][j][i].d *= (1.0 + amp/4.0*
	      (1.0+cos(2.0*PI*x2/Ly))*(1.0+cos(PI*x1)));
	  }
	  else if(iprob == 2){
	    pGrid->U[k][j][i].d *= (1.0 + 0.25 * amp*(ran2(&iseed) - 0.5)*(1.0+cos(PI*x1)));
	  }
	  else{
	   
	  }				
	}


      }
    }
  }

  /* Now ghost zone in the up J direction */
  for(j=js-nghost; j<=je+nghost; j++){
    for (i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, ie+i, j, k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	pressure = d1 * T1;
      pGrid->U[ks][j][ie+i].d  = d1;
      pGrid->U[ks][j][ie+i].M1  = 0.0;
      pGrid->U[ks][j][ie+i].M2  = 0.0;
      pGrid->U[ks][j][ie+i].M3  = 0.0;
      pGrid->U[ks][j][ie+i].E = pressure/(Gamma - 1.0);

      pGrid->U[ks][j][ie+i].Edd_11 = pGrid->U[ks][j][ie].Edd_11;
      pGrid->U[ks][j][ie+i].Edd_22 = pGrid->U[ks][j][ie].Edd_22; 
      pGrid->U[ks][j][ie+i].Edd_21 = 0.0;
      pGrid->U[ks][j][ie+i].Sigma[0] = kappaes *  pGrid->U[ks][j][ie].d;
      pGrid->U[ks][j][ie+i].Sigma[1] = kappaff * pGrid->U[ks][j][ie].d * pGrid->U[ks][j][ie].d * pow(pressure/pGrid->U[ks][j][ie].d,-3.5);
      pGrid->U[ks][j][ie+i].Sigma[2] = kappaff * pGrid->U[ks][j][ie].d * pGrid->U[ks][j][ie].d * pow(pressure/pGrid->U[ks][j][ie].d,-3.5);
      pGrid->U[ks][j][ie+i].Sigma[3] = kappaff * pGrid->U[ks][j][ie].d * pGrid->U[ks][j][ie].d * pow(pressure/pGrid->U[ks][j][ie].d,-3.5);

      pGrid->U[ks][j][ie+i].Er = pGrid->U[ks][j][ie+i].d * (kappaes + kappaff) * consFr * 
	                         (Ertop - x1) / pGrid->U[ks][j][ie+i].Edd_11;
      pGrid->U[ks][j][ie+i].Fr1 = consFr;
      pGrid->U[ks][j][ie+i].Fr2 = 0.0;
      pGrid->U[ks][j][ie+i].Fr3 = 0.0;
    }
  }

  for(j=js-nghost; j<=je+nghost; j++){
    for (i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, is-i, j, k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */

      pGrid->U[ks][j][is-i].d  = d1;
      pGrid->U[ks][j][is-i].M1  = 0.0;
      pGrid->U[ks][j][is-i].M2  = 0.0;
      pGrid->U[ks][j][is-i].M3  = 0.0;      
      pGrid->U[ks][j][is-i].E = pressure/(Gamma - 1.0);

      pGrid->U[ks][j][is-i].Edd_11 = pGrid->U[ks][js][is].Edd_11;
      pGrid->U[ks][j][is-i].Edd_22 = pGrid->U[ks][js][is].Edd_22;
      pGrid->U[ks][j][is-i].Edd_21 = 0.0;

      pGrid->U[ks][j][is-i].Sigma[0] = kappaes * pGrid->U[ks][js][is].d;
      pGrid->U[ks][j][is-i].Sigma[1] = kappaff * pGrid->U[ks][js][is].d * pGrid->U[ks][js][is].d * pow(pressure/pGrid->U[ks][js][is].d,-3.5);
      pGrid->U[ks][j][is-i].Sigma[2] = kappaff * pGrid->U[ks][js][is].d * pGrid->U[ks][js][is].d * pow(pressure/pGrid->U[ks][js][is].d,-3.5);
      pGrid->U[ks][j][is-i].Sigma[3] = kappaff * pGrid->U[ks][js][is].d * pGrid->U[ks][js][is].d * pow(pressure/pGrid->U[ks][js][is].d,-3.5);

      pGrid->U[ks][j][is-i].Er = pGrid->U[ks][j][is-i].d * (kappaes + kappaff) * consFr * 
	                         (Erbottom - x1) / pGrid->U[ks][j][is-i].Edd_11;
      pGrid->U[ks][j][is-i].Fr1 = consFr;
      pGrid->U[ks][j][is-i].Fr2 = 0.0;
      pGrid->U[ks][j][is-i].Fr3 = 0.0;
    }
  }

  Opacity = constopa;
  StaticGravPot = grav_pot;

/* set boundary functions */
  bvals_mhd_fun(pDomain, right_x1, radMHD_inflow);
  bvals_rad_fun(pDomain, right_x1, radMHD_rad_inflow);
  bvals_mhd_fun(pDomain, left_x1, radMHD_inflow2);
  bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow2);


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
#ifdef RAY_TRACING
	pRG->R[ifr][k][j][i].H[0] = 0.0;
#else
	pRG->R[ifr][k][j][i].H[0] = consFr;
#endif
      }


/* ------- Initialize boundary emission ---------------------------------- */
  for(k=kl; k<=ku; k++) {
    for(j=jl; j<=ju; j++) {
      for(m=0; m<nang; m++) {
	pRG->Ghstl1i[ifr][k][j][0][m] = Ibtm;
	if (noct > 2) {
	  pRG->Ghstl1i[ifr][k][j][2][m] = Ibtm;
	  if (noct == 8) {
	    pRG->Ghstl1i[ifr][k][j][4][m] = Ibtm;
	    pRG->Ghstl1i[ifr][k][j][6][m] = Ibtm;
	  }
	}
	pRG->Ghstr1i[ifr][k][j][1][m] = Itop;
	if (noct > 2) {
	  pRG->Ghstr1i[ifr][k][j][3][m] = Itop;
	  if (noct == 8) {
	    pRG->Ghstr1i[ifr][k][j][5][m] = Itop;
	    pRG->Ghstr1i[ifr][k][j][7][m] = Itop;
	  }
	}
      }}

    if (noct > 2) {
/* Initialize boundary intensity in x2 direction */

      for(i=il+1; i<=iu-1; i++) {
	/* periodic radiation at left boundary */
	for(m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = pGrid->U[ks][jl+nghost][i+nghost-1].Er;
	  pRG->Ghstl2i[ifr][k][i][1][m] = pGrid->U[ks][jl+nghost][i+nghost-1].Er;
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = pGrid->U[ks][jl+nghost][i+nghost-1].Er;
	    pRG->Ghstl2i[ifr][k][i][5][m] = pGrid->U[ks][jl+nghost][i+nghost-1].Er;	    
	  }
	}
	/* periodic radiation at right boundary */
	for(m=0; m<=nang; m++) {
	  pRG->Ghstr2i[ifr][k][i][2][m] = pGrid->U[ks][ju+nghost][i+nghost-1].Er;
	  pRG->Ghstr2i[ifr][k][i][3][m] = pGrid->U[ks][ju+nghost][i+nghost-1].Er;
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = pGrid->U[ks][ju+nghost][i+nghost-1].Er;
	    pRG->Ghstr2i[ifr][k][i][7][m] = pGrid->U[ks][ju+nghost][i+nghost-1].Er;
	  }
	}
      }

      for(m=0; m<nang; m++) {
	pRG->Ghstl2i[ifr][k][il][0][m] = Ibtm;
	pRG->Ghstr2i[ifr][k][il][2][m] = Ibtm;
	pRG->Ghstl2i[ifr][k][iu][1][m] = Itop;
	pRG->Ghstr2i[ifr][k][iu][3][m] = Itop;
	if (noct == 8) {
	  pRG->Ghstl2i[ifr][k][il][4][m] = Ibtm;
	  pRG->Ghstr2i[ifr][k][il][6][m] = Ibtm;
	  pRG->Ghstl2i[ifr][k][iu][5][m] = Itop;
	  pRG->Ghstr2i[ifr][k][iu][7][m] = Itop;
	}
      }
    }
  }

/* enrol user-defined  boundary functions */
  bvals_rad_trans_fun(pDomain, right_x1, const_J_ox1);
  bvals_rad_trans_fun(pDomain, left_x1, const_H_ix1);

/* enroll radiation specification functions */
  get_thermal_source = Thermal_B;
  get_thermal_fraction = const_eps;
  get_total_opacity = transfer_opacity;

/* ------- Initialize ray tracing boundary ----------------------- */
#ifdef RAY_TRACING
  for(ifr=0; ifr<pRG->nf_rt; ifr++) {
    for(k=pRG->ks; k<=pRG->ke; k++) {
      for(j=pRG->js; j<=pRG->je; j++) {	
	pRG->H[ifr][k][j][pRG->is-1] = consFr;
      }}}

  get_raytrace_thermal_fraction = const_eps;
  get_raytrace_opacity = transfer_opacity;
#endif /* RAY_TRACING */

#endif /* RADIATION_TRANSFER */

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

#if defined(RADIATION_MHD) || defined(RADIATION_HYDRO)
void bvals_mat_fun_ix1(VMatFun_t *Mat_BCFun)
{
  *Mat_BCFun = radMHD_Mat_inflowi2;
} 
void bvals_mat_fun_ox1(VMatFun_t *Mat_BCFun)
{
  *Mat_BCFun = radMHD_Mat_inflowo2;
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

  return grav*x1;
}


/* right hand side of x2 direction */

void radMHD_inflow(GridS *pGrid)
{

  int i, j;
  int ks, is, ie, js, je;
  Real u0;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  u0 = 0.0;	
 
  for(j=js-nghost; j<=je+nghost; j++) {
    for (i=1;  i<=nghost;  i++) {	
#ifdef RADIATION_MHD
      pGrid->B1i[ks][j][ie+i] = sqrt(fabs(u0));
      pGrid->U[ks][j][ie+i].B1c =  sqrt(fabs(u0));
#endif
      pGrid->U[ks][j][ie+i].M1 = - 0.0 * pGrid->U[ks][j][ie-i+1].M1;
      pGrid->U[ks][j][ie+i].M2 = pGrid->U[ks][j][ie-i+1].M2;
      pGrid->U[ks][j][ie+i].M3 = 0.0;
 

      pGrid->U[ks][j][ie+i].d  = d1;
      pGrid->U[ks][j][ie+i].Sigma[0] = pGrid->U[ks][j][ie+i-1].Sigma[0];
      pGrid->U[ks][j][ie+i].Sigma[1] = pGrid->U[ks][j][ie+i-1].Sigma[1];
      pGrid->U[ks][j][ie+i].Sigma[2] = pGrid->U[ks][j][ie+i-1].Sigma[2];  
      pGrid->U[ks][j][ie+i].Sigma[3] = pGrid->U[ks][j][ie+i-1].Sigma[3];
      pGrid->U[ks][j][ie+i].E = d1 * T1/(Gamma_1)+0.5*(SQR(pGrid->U[ks][j][ie+i].M1) + 
				 SQR(pGrid->U[ks][j][ie+i].M2)) / pGrid->U[ks][j][ie+i].d;
    }
  }
  return;
}

void radMHD_rad_inflow(GridS *pGrid)
{
  int i, j;
  int ks, is, ie, js, je;
  Real x1, x2, x3;
  Real Fr0x, Fr0y;
  Real velocity_x, velocity_y;
  Real velocity_x1, velocity_y1;
  Real Sigma_t, Sigma_t1;
  Real Eratio, reducefactor, tau;
  Real dz = pGrid->dx3;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;
	
  for(j=js-nghost; j<=je+nghost; j++){
    for(i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, ie+i, j,ks, &x1, &x2, &x3);

      pGrid->U[ks][j][ie+i].Edd_11 = pGrid->U[ks][j][ie].Edd_11;
      pGrid->U[ks][j][ie+i].Edd_22 = pGrid->U[ks][j][ie].Edd_22;
      pGrid->U[ks][j][ie+i].Edd_21 = pGrid->U[ks][j][ie].Edd_21;
      
      matrix_alpha(0.0, pGrid->U[ks][j][ie+i].Sigma, pGrid->dt, pGrid->U[ks][j][ie+i].Edd_33, 0.0, &reducefactor, 0, dz);
      
      Sigma_t = 0.5 * (pGrid->U[ks][j][ie+i].Sigma[0] + pGrid->U[ks][j][ie+i].Sigma[1] + pGrid->U[ks][j][ie+i-1].Sigma[0] + pGrid->U[ks][j][ie+i-1].Sigma[1]);

      Sigma_t1 = 0.5 * (pGrid->U[ks][j][ie+i-1].Sigma[0] + pGrid->U[ks][j][ie+i-1].Sigma[1] + pGrid->U[ks][j][ie+i-2].Sigma[0] + pGrid->U[ks][j][ie+i-2].Sigma[1]);

      velocity_x = pGrid->U[ks][j][ie+i-1].M1 / pGrid->U[ks][j][ie+i-1].d;
      velocity_y = pGrid->U[ks][j][ie+i-1].M2 / pGrid->U[ks][j][ie+i-1].d;

      Fr0x = pGrid->U[ks][j][ie+i-1].Fr1 - ((1.0 + pGrid->U[ks][j][ie+i-1].Edd_11) * velocity_x + pGrid->U[ks][j][ie+i-1].Edd_21 * velocity_y)* pGrid->U[ks][j][ie+i-1].Er / Crat;

      Fr0y = pGrid->U[ks][j][ie+i-1].Fr2 - ((1.0 + pGrid->U[ks][j][ie+i-1].Edd_22) * velocity_y + pGrid->U[ks][j][ie+i-1].Edd_21 * velocity_x)* pGrid->U[ks][j][ie+i-1].Er / Crat;

      Fr0x = consFr;

      velocity_x1 = pGrid->U[ks][j][ie+i].M1 / pGrid->U[ks][j][ie+i].d;
      velocity_y1 = pGrid->U[ks][j][ie+i].M2 / pGrid->U[ks][j][ie+i].d;

      /* use diffusion boundary */
      if(Fr0x > 0.0)
	pGrid->U[ks][j][ie+i].Er = pGrid->U[ks][j][ie+i-1].Er - pGrid->dx1 * Sigma_t1 * Fr0x / pGrid->U[ks][j][ie+i].Edd_11;
		
      if(pGrid->U[ks][j][ie+i].Er < 0.0)
	pGrid->U[ks][j][ie+i].Er = pGrid->U[ks][j][ie].Er;

      pGrid->U[ks][j][ie+i].Fr1 = Fr0x + ((1.0 + pGrid->U[ks][j][ie+i].Edd_11) * velocity_x1 + pGrid->U[ks][j][ie+i].Edd_21 * velocity_y1)* pGrid->U[ks][j][ie+i].Er / Crat;

      pGrid->U[ks][j][ie+i].Fr2 = Fr0y + ((1.0 + pGrid->U[ks][j][ie+i].Edd_22) * velocity_y1 + pGrid->U[ks][j][ie+i].Edd_21 * velocity_x1)* pGrid->U[ks][j][ie+i].Er / Crat;

      pGrid->U[ks][j][ie+i].Fr3 = 0.0;
    }
  }

  return;
}


void radMHD_inflow2(GridS *pGrid)
{
  int i, j;
  int ks, is, ie, js, je;
  Real u0, opac;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  u0 = 0.0;
 
  for(j=js-nghost; j<=je+nghost; j++) {
    opac = SQR(pGrid->U[ks][j][is].d) * pow(pressure/pGrid->U[ks][j][is].d,-3.5);
    for (i=1;  i<=nghost;  i++) {	
#ifdef RADIATION_MHD
      pGrid->B1i[ks][j][is-i] = sqrt(fabs(u0));
      pGrid->U[ks][j][is-i].B1c =  sqrt(fabs(u0));
#endif
		
      pGrid->U[ks][j][is-i].M1 = pGrid->U[ks][j][is+i-1].M1;
      pGrid->U[ks][j][is-i].M2 = -0.0 * pGrid->U[ks][j][is+i-1].M2;
      pGrid->U[ks][j][is-i].M3 = 0.0;
      pGrid->U[ks][j][is-i].d  = d1;
      pGrid->U[ks][j][is-i].Sigma[0] = kappaes *  pGrid->U[ks][j][is-i].d;
      pGrid->U[ks][j][is-i].Sigma[1] = kappaff * opac;
      pGrid->U[ks][j][is-i].Sigma[2] = kappaff * opac;
      pGrid->U[ks][j][is-i].Sigma[3] = kappaff * opac;      
      pGrid->U[ks][j][is-i].E  = d1 * T1/(Gamma_1)+0.5*(SQR(pGrid->U[ks][j][is-i].M1) + 
				 SQR(pGrid->U[ks][j][is-i].M2)) / pGrid->U[ks][j][is-i].d;

    }
  }
  return;
}


void radMHD_rad_inflow2(GridS *pGrid)
{
  int i, j;
  int ks, is, ie, js, je;
  Real x1, x2, x3;
  Real Sigma_t1;
  Real velocity_x, velocity_y, velocity_x1, velocity_y1;
  Real Fr0x, Fr0y;

  je = pGrid->je;
  js = pGrid->js;
  ks = pGrid->ks;
  is = pGrid->is;
  ie = pGrid->ie;

  for(j=js-nghost; j<=je+nghost; j++){
    for(i=1;  i<=nghost;  i++) {
      cc_pos(pGrid, is-i, j,ks, &x1, &x2, &x3);

      pGrid->U[ks][j][is-i].Edd_11 = pGrid->U[ks][j][is].Edd_11;
      pGrid->U[ks][j][is-i].Edd_22 = pGrid->U[ks][j][is].Edd_22;
      pGrid->U[ks][j][is-i].Edd_21 = pGrid->U[ks][j][is].Edd_21;
      
      Sigma_t1 = 0.5 * (pGrid->U[ks][j][is-i+1].Sigma[0] + pGrid->U[ks][j][is-i+1].Sigma[1] + pGrid->U[ks][j][is-i+2].Sigma[0] + pGrid->U[ks][j][is-i+2].Sigma[1]);	
		
      velocity_x = pGrid->U[ks][j][is-i+1].M1 / pGrid->U[ks][j][is-i+1].d;
      velocity_y = pGrid->U[ks][j][is-i+1].M2 / pGrid->U[ks][j][is-i+1].d;
      
      Fr0x = consFr;

      //      Fr0x = pGrid->U[ks][j][is-i+1].Fr1 - ((1.0 + pGrid->U[ks][j][is-i+1].Edd_11) * velocity_x + pGrid->U[ks][j][is-i+1].Edd_21 * velocity_y)* pGrid->U[ks][j][is-i+1].Er / Crat;

      Fr0y = pGrid->U[ks][j][is-i+1].Fr2 - ((1.0 + pGrid->U[ks][j][is-i+1].Edd_22) * velocity_y + pGrid->U[ks][j][is-i+1].Edd_21 * velocity_x)* pGrid->U[ks][j][is-i+1].Er / Crat;

      velocity_x1 = pGrid->U[ks][j][is-i].M1 / pGrid->U[ks][j][is-i].d;
      velocity_y1 = pGrid->U[ks][j][is-i].M2 / pGrid->U[ks][j][is-i].d;

      /* use diffusion boundary */
      if(Fr0x > 0.0)
	pGrid->U[ks][j][is-i].Er = pGrid->U[ks][j][is-i+1].Er + pGrid->dx1 * Sigma_t1 * Fr0x / pGrid->U[ks][j][is-i].Edd_11;
		
      pGrid->U[ks][j][is-i].Fr1 = Fr0x + ((1.0 + pGrid->U[ks][j][is-i].Edd_11) * velocity_x1 + pGrid->U[ks][j][is-i].Edd_21 * velocity_y1)* pGrid->U[ks][j][is-i].Er / Crat;

      pGrid->U[ks][j][is-i].Fr2 = Fr0y + ((1.0 + pGrid->U[ks][j][is-i].Edd_22) * velocity_y1 + pGrid->U[ks][j][is-i].Edd_21 * velocity_x1)* pGrid->U[ks][j][is-i].Er / Crat;

      pGrid->U[ks][j][is-i].Fr3 = 0.0;
    }
  }

  return;
}


void radMHD_Mat_inflowi2(MatrixS *pMat)
{
  int i, j;
  int ks, is, ie, js, je;
  Real x1, x2, x3;

  je = pMat->je;
  js = pMat->js;
  ks = pMat->ks;
  is = pMat->is;
  ie = pMat->ie;

  for (j=js-Matghost;  j<=je+Matghost;  j++) {	
    for(i=1; i<=Matghost; i++) {
      pMat->Ugas[ks][j][is-i].Edd_11 = pMat->Ugas[ks][j][is].Edd_11;
      pMat->Ugas[ks][j][is-i].Edd_22 = pMat->Ugas[ks][j][is].Edd_22;
      pMat->Ugas[ks][j][is-i].Edd_21 = pMat->Ugas[ks][j][is].Edd_21;
	
      pMat->Ugas[ks][j][is-i].V1 =  pMat->Ugas[ks][j][is+i-1].V1;
      pMat->Ugas[ks][j][is-i].V2 = -pMat->Ugas[ks][j][is+i-1].V2;
      pMat->Ugas[ks][j][is-i].T4 = pMat->Ugas[ks][j][is].T4;
      pMat->Ugas[ks][j][is-i].Sigma[0] = pMat->Ugas[ks][j][is].Sigma[0];
      pMat->Ugas[ks][j][is-i].Sigma[1] = pMat->Ugas[ks][j][is].Sigma[1];
      pMat->Ugas[ks][j][is-i].Sigma[2] = pMat->Ugas[ks][j][is].Sigma[2];
      pMat->Ugas[ks][j][is-i].Sigma[3] = pMat->Ugas[ks][j][is].Sigma[3];

      pMat->U[ks][j][is-i].Er  = 0.0;
      pMat->U[ks][j][is-i].Fr1 = 0.0;
      pMat->U[ks][j][is-i].Fr2 = 0.0;
      pMat->U[ks][j][is-i].Fr3 = 0.0;		
      
    }
  }

  return;
}



void radMHD_Mat_inflowo2(MatrixS *pMat)
{
  int i, j;
  int ks, is, ie, js, je;
  Real x1, x2, x3;

  je = pMat->je;
  js = pMat->js;
  ks = pMat->ks;
  is = pMat->is;
  ie = pMat->ie;
	
  for (j=js-Matghost;  j<=je+Matghost;  j++) {
    for(i=1; i<=Matghost; i++){
      pMat->Ugas[ks][j][ie+i].Edd_11 = pMat->Ugas[ks][j][ie].Edd_11;
      pMat->Ugas[ks][j][ie+i].Edd_22 = pMat->Ugas[ks][j][ie].Edd_22;
      pMat->Ugas[ks][j][ie+i].Edd_21 = pMat->Ugas[ks][j][ie].Edd_21;
      
      pMat->Ugas[ks][j][ie+i].V1 = pMat->Ugas[ks][j][ie].V1;
      pMat->Ugas[ks][j][ie+i].V2 = -pMat->Ugas[ks][j][ie-i+1].V2;
      pMat->Ugas[ks][j][ie+i].T4 = pMat->Ugas[ks][j][ie].T4;
      pMat->Ugas[ks][j][ie+i].Sigma[0] = pMat->Ugas[ks][j][ie].Sigma[0];
      pMat->Ugas[ks][j][ie+i].Sigma[1] = pMat->Ugas[ks][j][ie].Sigma[1];
      pMat->Ugas[ks][j][ie+i].Sigma[2] = pMat->Ugas[ks][j][ie].Sigma[2];
      pMat->Ugas[ks][j][ie+i].Sigma[3] = pMat->Ugas[ks][j][ie].Sigma[3];

      pMat->U[ks][j][ie+i].Er  = 0.0;
      pMat->U[ks][j][ie+i].Fr1 = 0.0;
      pMat->U[ks][j][ie+i].Fr2 = 0.0;

      pMat->U[ks][j][ie+i].Fr3 = 0.0;		
    }
  }

  return;
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

void const_H_ix1(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int il = pRG->is-1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, k, l, m, n, ifr;
  int ig,kg,ioff,joff,koff;
  Real I0, Jm, H, Hm, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][0] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  ioff = nghost - 1;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1;
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1;
  } else koff = 0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstl1i using l1imu */
    for (k=kl; k<=ku; k++) {
      kg = k + koff;
      for (j=jl; j<=ju; j++) {
	Hm = 0.0;
	Jm = 0.0;
	for (m=0; m<nang; m++) {
	  Hm += pRG->l1imu[ifr][k][j][1][m] * pRG->mu[1][m][0] * pRG->wmu[m];
	  Hm += pRG->l1imu[ifr][k][j][3][m] * pRG->mu[3][m][0] * pRG->wmu[m];
	  Jm += pRG->l1imu[ifr][k][j][1][m] * pRG->wmu[m];
	  Jm += pRG->l1imu[ifr][k][j][3][m] * pRG->wmu[m];	  
	  if (noct == 8) {
	    Hm += pRG->l1imu[ifr][k][j][5][m] * pRG->mu[5][m][0] * pRG->wmu[m];
	    Hm += pRG->l1imu[ifr][k][j][7][m] * pRG->mu[7][m][0] * pRG->wmu[m];
	    Jm += pRG->l1imu[ifr][k][j][5][m] * pRG->wmu[m];
	    Jm += pRG->l1imu[ifr][k][j][7][m] * pRG->wmu[m];
	  }
	}
#ifdef RAY_TRACING 
	H = 0.0;
#else
	H = consFr;
#endif
	I0 = (H - Hm) / gamma;
	if (I0 < 0.0) I0 = 0.0; 
	pRG->R[ifr][k][j][il].J = 0.5 * I0 + Jm;
	for (m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = I0;
	  pRG->Ghstl1i[ifr][k][j][2][m] = I0;
	  if (noct == 8) {
	    pRG->Ghstl1i[ifr][k][j][4][m] = I0;
	    pRG->Ghstl1i[ifr][k][j][6][m] = I0;
	  }
	}
      }
/* update r2imu and l2imu so corner intensities are correct w/ periodic bcs*/
      for (m=0; m<nang; m++) {
	pRG->r2imu[ifr][k][il][0][m] = pRG->Ghstl1i[ifr][k][ju][0][m];
	pRG->l2imu[ifr][k][il][2][m] = pRG->Ghstl1i[ifr][k][jl][2][m];
      }
    }
  }
  return;
}

void const_J_ox1(GridS *pG, RadGridS *pRG, int ifs, int ife)
{
  int iu = pRG->ie+1;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int i, j, k, l, m, n, ifr;
  Real I0, Jp, J, Hp, gamma = 0.0;

  /* gamma ~ 1/4 */
  for (m=0; m<nang; m++) {
    gamma += pRG->mu[0][m][0] * pRG->wmu[m];
  }
  if (noct == 8) gamma *= 4.0; else gamma *= 2.0;

  for(ifr=ifs; ifr<=ife; ifr++) {
/* update Ghstr2i using r2imu */
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	Hp = 0.0;
	Jp = 0.0;
	for (m=0; m<nang; m++) {
	  Hp += pRG->r1imu[ifr][k][j][0][m] * pRG->mu[0][m][0] * pRG->wmu[m];
	  Hp += pRG->r1imu[ifr][k][j][2][m] * pRG->mu[2][m][0] * pRG->wmu[m];
	  Jp += pRG->r1imu[ifr][k][j][0][m] * pRG->wmu[m];
	  Jp += pRG->r1imu[ifr][k][j][2][m] * pRG->wmu[m];
	  if (noct == 8) {
	    Hp += pRG->r1imu[ifr][k][j][4][m] * pRG->mu[4][m][0] * pRG->wmu[m];
	    Hp += pRG->r1imu[ifr][k][j][6][m] * pRG->mu[6][m][0] * pRG->wmu[m];
	    Jp += pRG->r1imu[ifr][k][j][4][m] * pRG->wmu[m];
	    Jp += pRG->r1imu[ifr][k][j][6][m] * pRG->wmu[m];
	  }
	}
#ifdef RAY_TRACING
	Jp += pRG->H[ifr][k][j][iu];
#endif
	J = pRG->R[ifr][k][j][iu].J;
	I0 = 2.0 * (J - Jp); 
	if (I0 < 0.0) I0 = 0.0; 
	for (m=0; m<nang; m++) {
	  pRG->Ghstr1i[ifr][k][j][1][m] = I0;
	  pRG->Ghstr1i[ifr][k][j][3][m] = I0;
	  if (noct == 8) {
	    pRG->Ghstr1i[ifr][k][j][5][m] = I0;
	    pRG->Ghstr1i[ifr][k][j][7][m] = I0;
	  }
	}
      }
/* update r2imu and l2imu so corner intensities are correct w/ periodic bcs*/
      for (m=0; m<nang; m++) {
	pRG->r2imu[ifr][k][iu][1][m] = pRG->Ghstr1i[ifr][k][ju][1][m];
	pRG->l2imu[ifr][k][iu][3][m] = pRG->Ghstr1i[ifr][k][jl][3][m];
      }
    }
  }

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

void Userwork_after_first_formal_solution(DomainS *pD)
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
 
  Kloc = 0.0; Hloc = 0.0; 
  for (j = pRG->js; j <= pRG->je; j++) {
    if(pG->rx1_id == -1) {
      Kloc += pRG->R[0][pRG->ks][j][pRG->ie].K[0];
      Hloc += pRG->R[0][pRG->ks][j][pRG->ie].H[0];
#ifdef RAY_TRACING
      Kloc += pRG->H[0][pRG->ks][j][pRG->ie];
      Hloc += pRG->H[0][pRG->ks][j][pRG->ie];
#endif
    }
  }

#ifdef MPI_PARALLEL
  MPI_Allreduce(&Kloc, &K, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&Hloc, &H, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
  K = Kloc;
  H = Hloc;
#endif /* MPI_PARALLEL */

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
	K = pRG->R[0][kr][jr][ir].K[0];
#ifdef RAY_TRACING
	J += pRG->H[0][kr][jr][ir];
	K += pRG->H[0][kr][jr][ir];
#endif
	if((fabs(J) > TINY_NUMBER) && (J > 0.0)) {
	  pG->U[k][j][i].Edd_11 = K / J;
	}
	cc_pos(pG,i,j,k,&x1,&x2,&x3);
	if(x1 < yintb) {
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Erbottom - x1) / pG->U[k][j][i].Edd_11;
	} else if(x1 < yintt) {
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Ermid - x1) / pG->U[k][j][i].Edd_11;
	} else{
	  pG->U[k][j][i].Er = pG->U[k][j][i].d * (kappaes + kappaff) * consFr * 
	                      (Ertop - x1) / pG->U[k][j][i].Edd_11;
	}
      }}}
}

#endif /* End radiation transfer */


