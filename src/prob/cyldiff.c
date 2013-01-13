#include "copyright.h"
/*============================================================================*/
/*! \file cyldiff.c
 *  \brief A simple magnetic diffusion test in cylindrical coordinates
 *
 * PURPOSE:  A simple magnetic diffusion test in cylindrical coordinates
 *   withno pressure or tension forces.
 *
 * AUTHOR:  Chang-Goo Kim
 *
 * IMPORTANT NOTICE:
 *   To perform this test, the main integrator MUST be turned off!
 *   iprob = 1; Eigenmode test;
 *     IC Bz(r,phi,0)=J_m(lambda_m,n*R)*cos(m*phi);
 *     BC Bz(0,phi,t)<infty, Bz(1,phi,t)=0 for R and periodic for phi;
 *     Solution: Bz(r,phi,t)=Bz(r,phi,0)*exp(-eta_Ohm/lambda_m,n^2*t)
 *   iprob = 2; 1D Cylindrical test with IC Bz(x,0)=1-x^2
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef MHD
#error : The diffusion test only works for mhd.
#endif

#ifndef RESISTIVITY
#error : The diffusion test only works with resistivity.
#endif


/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * cyldiff_ix1() - Inner-R boundary conditions for diffusion problem
 * cyldiff_ox1() - Outer-R boundary conditions for diffusion problem
 * cyldiff_ix2() - Inner-phi boundary conditions for diffusion problem
 * cyldiff_ox2() - Outer-phi boundary conditions for diffusion problem
 * ROOTJ(), SECANT() - CALCULATE THE FIRST NK ZEROES OF BESSEL FUNCTION J(N,X)
 *                     Obtained from http://jean-pierre.moreau.pagesperso-orange.fr/c_bessel.html
 *============================================================================*/

static int iprob,nn,mm;
static Real lambdamn,om_diff;
static Real ***RootSoln=NULL;
void cyldiff_ix1(GridS *pG);
void cyldiff_ox1(GridS *pG);
void cyldiff_ix2(GridS *pG);
void cyldiff_ox2(GridS *pG);
void SECANT(int N,int NITMX, double TOL, double *ZEROJ, int *IER);
void ROOTJ(int N, int NK, double *JZERO, int *IER);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */
void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3;
  int *ierr;
  double *lambda;
  Real x1c,x2c,x3c;
  Real x1f,x2f,x3f;
#ifndef CYLINDRICAL
  Real r,ri;
#endif

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

  /* allocate memory for solution */
  if ((RootSoln = (Real***)calloc_3d_array(nx3,nx2,nx1,sizeof(Real))) == NULL)
    ath_error("[cyldiff]: Error allocating memory for solution\n");

  /* parse input file */
  iprob  = par_geti("problem", "iprob");
  mm = par_geti("problem","m");
  nn = par_geti("problem","n");
  
  if ((lambda = (double *)malloc(nn*sizeof(double))) == NULL)
    ath_error("[problem] Error initializing lambda array");
  if ((ierr = (int *)malloc(nn*sizeof(int))) == NULL)
    ath_error("[problem] Error initializing ierr array");

  ROOTJ(mm,nn,lambda,ierr);
  lambdamn=lambda[nn];
  printf("m=%i, n=%i, lambda_m,n=%g\n",mm,nn,lambdamn);
  par_setd("problem","lambda","%e",lambdamn,"lambda_m,n");

  /* Set density and phi-momentum */
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1c,&x2c,&x3c);
        fc_pos(pG,i,j,k,&x1f,&x2f,&x3f);
        memset(&(pG->U[k][j][i]),0.0,sizeof(ConsS));

        pG->U[k][j][i].d  = 1.0;
        switch (iprob) {
          case 1:
#ifndef CYLINDRICAL
  ath_error("[cyldiff]: iprob=1 only works in cylindrical!\n");
#endif
#ifdef MHD
                  pG->U[k][j][i].B3c = jn(mm,lambdamn*x1c)*cos(mm*x2c);
                  pG->B3i[k][j][i]   = jn(mm,lambdamn*x1f)*cos(mm*x2f);
#endif

                  break;
          case 2:
#ifndef CYLINDRICAL
  ath_error("[cyldiff]: iprob=2 only works in cylindrical!\n");
#endif

#ifdef MHD
                  pG->U[k][j][i].B3c = 1.0-SQR(x1c);
                  pG->B3i[k][j][i]   = 1.0-SQR(x1f);
#endif
                  break;
          default:
                  ath_error("[cyldiff]:  Not an accepted problem number\n");
        }

        RootSoln[k][j][i] = pG->U[k][j][i].B3c;
      }
    }
  }

  bvals_mhd_fun(pDomain, left_x1,  cyldiff_ix1);
  bvals_mhd_fun(pDomain, right_x1, cyldiff_ox1);
/*
  if (iprob == 2) {
    bvals_mhd_fun(pDomain, left_x2,  cyldiff_ix2);
    bvals_mhd_fun(pDomain, right_x2, cyldiff_ox2);
  }
*/

#ifdef RESISTIVITY 
  eta_Ohm = par_getd("problem","eta_Ohm");
  Q_AD    = 0.0;
  Q_Hall  = 0.0;
  d_ind   = 0.0;
#endif

  om_diff=eta_Ohm*SQR(lambdamn); // diffusion rate


  return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
#ifdef RESISTIVITY 
  eta_Ohm = par_getd("problem","eta_Ohm");
  Q_AD    = 0.0;
  Q_Hall  = 0.0;
  d_ind   = 0.0;
#endif

  return;
}

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
#ifdef MHD
  printf("Max divB = %1.10e\n", compute_div_b(pM->Domain[0][0].Grid));
#endif
}

void Userwork_after_loop(MeshS *pM)
{
  GridS *pG = pM->Domain[0][0].Grid;
  int i,j,k;
  int is, ie, js, je, ks, ke;
  int Nx1, Nx2, Nx3, count;
  FILE *fp;
  char *fname;
  Real error=0.0,tfact;

  tfact=exp(-om_diff*pG->time);

  is = pG->is;  ie = pG->ie;
  js = pG->js;  je = pG->je;
  ks = pG->ks;  ke = pG->ke;

  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;
  count = Nx1*Nx2*Nx3;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        error += fabs(pG->U[k][j][i].B3c - RootSoln[k][j][i]*tfact);
      }
    }
  }

  error /= (Real)count;

  fname = ath_fname(NULL,"cyldiff-errors",NULL,NULL,1,10*mm+nn,NULL,"dat");
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_perr(-1,"[Userwork_after_loop]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_perr(-1,"[Userwork_after_loop]: Unable to open file.\n");
      free(fname);
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  Error");
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,error);
  fprintf(fp,"\n");
  fclose(fp);
  free(fname);

  return;


}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void cyldiff_ix1(GridS *pG)
 *  \brief Inner-R boundary conditions for diffusion problem
 */

void cyldiff_ix1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1c,x2c,x3c;
  Real x1f,x2f,x3f;
  Real Eint,Emag,Ekin;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        cc_pos(pG,is-i,j,k,&x1c,&x2c,&x3c);
        fc_pos(pG,is-i,j,k,&x1f,&x2f,&x3f);
        switch(iprob){
          case 1:        
#ifdef MHD

            pG->U[k][j][is-i].B3c = jn(mm,lambdamn*x1c)*cos(mm*x2c)*exp(-om_diff*pG->time);
            pG->B3i[k][j][is-i]   = jn(mm,lambdamn*x1f)*cos(mm*x2f)*exp(-om_diff*pG->time);
#endif
            break;
          default:
#ifdef MHD
            pG->U[k][j][is-i].B3c = pG->U[k][j][is].B3c;
            pG->B3i[k][j][is-i]   = pG->B3i[k][j][is];
#endif
        }

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void cyldiff_ox1(GridS *pG)
 *  \brief Outer-R boundary conditions for diffusion problem
 */

void cyldiff_ox1(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3;
#ifdef CYLINDRICAL
  const Real *r=pG->r, *ri=pG->ri;
#endif
  Real rsf=1., lsf=1.;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        switch(iprob){
          case 1:        
#ifdef MHD
            pG->U[k][j][ie+i].B3c = 0.;
            pG->B3i[k][j][ie+i]   = 0.;
#endif
            break;
          default:
#ifdef MHD
            pG->U[k][j][ie+i].B3c = pG->U[k][j][ie].B3c;
            pG->B3i[k][j][ie+i]   = pG->B3i[k][j][ie];
#endif
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void cyldiff_ix2(GridS *pG)
 *  \brief Inner-phi boundary conditions for diffusion problem
 */

void cyldiff_ix2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3;

  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      for (j=1; j<=nghost; j++) {
        cc_pos(pG,i,js-j,k,&x1,&x2,&x3);
        switch(iprob){
          case 1:        
#ifdef MHD
            pG->U[k][js-j][i].B3c = 0.;
            pG->B3i[k][js-j][i]   = 0.;
#endif
            break;
          default:
#ifdef MHD
            pG->U[k][js-j][i].B3c = pG->U[k][js][i].B3c;
            pG->B3i[k][js-j][i]   = pG->B3i[k][js][i];
#endif
        }

      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void cyldiff_ox2(GridS *pG)
 *  \brief Outer-phi boundary conditions for diffusion problem
 */

void cyldiff_ox2(GridS *pG)
{
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int i,j,k;
  Real x1,x2,x3;

  for (k=ks; k<=ke; k++) {
    for (i=is; i<=ie; i++) {
      for (j=1; j<=nghost; j++) {
        cc_pos(pG,i,je+j,k,&x1,&x2,&x3);
        switch(iprob){
          case 1:        
#ifdef MHD
            pG->U[k][je+j][i].B3c = 0.;
            pG->B3i[k][je+j][i]   = 0.;
#endif
            break;
          default:
#ifdef MHD
            pG->U[k][je+j][i].B3c = pG->U[k][je][i].B3c;
            pG->B3i[k][je+j][i]   = pG->B3i[k][je][i];
#endif
        }
      }
    }
  }

  return;
}

/* -------------------------------------------------------------------------
! Reference: From Numath Library By Tuan Dang Trong in Fortran 77
!            [BIBLI 18].
!
!                               C++ Release 1.0 By J-P Moreau, Paris
! ------------------------------------------------------------------------ */


//----------------------------------------------------------------------
void ROOTJ(int N, int NK, double *JZERO, int *IER)  {
/*----------------------------------------------------------------------
!     CALCULATE THE FIRST NK ZEROES OF BESSEL FUNCTION J(N,X)
!
!     INPUTS:
!       N    ORDER OF FUNCTION J (INTEGER >= 0)                  I*4
!       NK   NUMBER OF FIRST ZEROES  (INTEGER > 0)               I*4
!     OUTPUTS:
!       JZERO(NK)  TABLE OF FIRST ZEROES (ABCISSAS)              R*8
!       IER(NK)    TABLE OF ERROR CODES (MUST BE ZEROES)         I*4
!
!     REFERENCE :
!     ABRAMOWITZ M. & STEGUN IRENE A.
!     HANDBOOK OF MATHEMATICAL FUNCTIONS
! -------------------------------------------------------------------- */
  double ZEROJ,B0,B1,B2,B3,B5,B7,T0,T1,T3,T5,T7,FN,FK;
  double C1,C2,C3,C4,C5,F1,F2,F3,TOL,ERRJ;
  int IERROR,K,NITMX;

  TOL=1E-8; NITMX=10;
  C1=1.8557571; C2=1.033150; C3=0.00397; C4=0.0908; C5=0.043;
  FN = 1.0*N;

//    FIRST ZERO
  if (N==0) {
    ZEROJ = C1+C2-C3-C4+C5;
    SECANT(N,NITMX,TOL,&ZEROJ,&IERROR);
    IER[1]=IERROR;
    JZERO[1]=ZEROJ;
  }
  else {
    F1 = pow(FN,(1.0/3.0));
    F2 = F1*F1*FN;
    F3 = F1*FN*FN;
    ZEROJ = FN+C1*F1+(C2/F1)-(C3/FN)-(C4/F2)+(C5/F3);
    SECANT(N,NITMX,TOL,&ZEROJ,&IERROR);
    IER[1]=IERROR;
    JZERO[1]=ZEROJ;
  }

  T0 = 4.0*FN*FN;
  T1 = T0-1.0;
  T3 = 4.0*T1*(7.0*T0-31.0);
  T5 = 32.0*T1*((83.0*T0-982.0)*T0+3779.0);
  T7 = 64.0*T1*(((6949.0*T0-153855.0)*T0+1585743.0)*T0-6277237.0);

//    OTHER ZEROES
  for (K = 2; K <= NK; K++) {
    JZERO[K] = 0.0;
    FK = 1.0*K;

//      MAC MAHON'S SERIES FOR K>>N

    B0 = (FK+0.5*FN-0.25)*PI;
    B1 = 8.0*B0;
    B2 = B1*B1;
    B3 = 3.0*B1*B2;
    B5 = 5.0*B3*B2;
    B7 = 7.0*B5*B2;

    ZEROJ = B0-(T1/B1)-(T3/B3)-(T5/B5)-(T7/B7);

    ERRJ=fabs(jn(N,ZEROJ));

//      IMPROVE SOLUTION USING PROCEDURE SECANT

    if (ERRJ > TOL) SECANT(N,NITMX,TOL,&ZEROJ,&IERROR);
    JZERO[K]=ZEROJ;
    IER[K]=IERROR;
  }
}
// ------------------------------------------------------------------------------
void SECANT(int N,int NITMX, double TOL, double *ZEROJ, int *IER)  {
    //Labels: e5,e10,e15,e20
    
  double P0,P1,Q0,Q1,DP,P;
  double C[3];
  int IT,NEV,NTRY;

  C[1]=0.95; C[2]=0.999;
  NTRY=1;
  *IER=0;

e5:   P0 = C[NTRY]*(*ZEROJ);

  P1 = *ZEROJ;
  NEV = 2;
  Q0 = jn(N,P0);
  Q1 = jn(N,P1);
  for (IT = 1; IT <= NITMX; IT++) {
    if (Q1==Q0) goto e15;
    P = P1-Q1*(P1-P0)/(Q1-Q0);
    DP = P-P1;
    if (IT==1) goto e10;
    if (fabs(DP) < TOL) goto e20;

e10:    NEV = NEV+1;
    P0 = P1;
    Q0 = Q1;
    P1 = P;
    Q1 = jn(N,P1);
  }
e15:  NTRY++;
  if (NTRY <= 2) goto e5;
  *IER = NTRY;
e20:  *ZEROJ = P;
}
