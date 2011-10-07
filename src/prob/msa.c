#include "copyright.h"
/*==============================================================================
 * FILE: msa.c
 *
 * PURPOSE:  Problem generator for MSA test
 *   3D shearing box code. 
 *
 * Code must be configured using --enable-shearing-box
 *
 * REFERENCE: Kim, W.-T and Ostriker, E. C. (2001)
 * Wriiten by Chang-Goo Kim
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
 *============================================================================*/
static void slab_grav_bc_lower(GridS *pG);
static void slab_grav_bc_upper(GridS *pG);
static void slab_mhd_bc_lower(GridS *pG);
static void slab_mhd_bc_upper(GridS *pG);
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3);
static Real hst_sigma(const GridS *pG, const int i, const int j, const int k);
static Real hst_ux(const GridS *pG, const int i, const int j, const int k);
static Real hst_uy(const GridS *pG, const int i, const int j, const int k);
static Real hst_dSigma(const GridS *pG, const int i, const int j, const int k);
static Real hst_Vx(const GridS *pG, const int i, const int j, const int k);
static Real hst_dVy(const GridS *pG, const int i, const int j, const int k);
#ifdef MHD
static Real hst_m1(const GridS *pG, const int i, const int j, const int k);
static Real hst_m2(const GridS *pG, const int i, const int j, const int k);
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k);
static Real hst_dBy(const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef SELF_GRAVITY
static Real hst_Phi(const GridS *pG, const int i, const int j, const int k);
static Real hst_dPhi(const GridS *pG, const int i, const int j, const int k);
#endif
#ifdef ADIABATIC
static Real hst_dE(const GridS *pG, const int i, const int j, const int k);
#endif
static Real dVol;

/*=========================== PUBLIC FUNCTIONS =================================
 * Contains the usual, plus:
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pG = pDomain->Grid;
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int ixs,jxs,kxs,i,j,k;
  long int iseed = -1; /* Initialize on the first call to ran2 */
  Real x1,x2,x3,xmin,xmax,Lx,Ly,Lz;
  Real rd, rp, rvx, rvy, rvz, rbx, rby, rbz;
  Real beta,B0,P0,kx,ky,kz,amp,press;
  Real Q,nJ,cs,cs2;
  Real time0,kxt;
#ifdef SELF_GRAVITY
  Real Gcons;
#endif

  int nwx,nwy,nwz;  /* input number of waves per Lx,Ly,Lz [default=1] */
  double rval;

  if(pG->Nx[2] == 1) ShBoxCoord = xy; /* 2D xy-plane */

/* Read problem parameters. */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd("problem","qshear");
  amp = par_getd("problem","amp");

/* Read parameters for magnetic field */
  beta = par_getd("problem","beta"); 

/* Read parameters for self gravity */
  Q=par_getd("problem","Q");
  nJ= par_getd("problem","nJ");

  time0=par_getd_def("problem","time0",0.0);

  cs=sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  cs2=SQR(cs);

#ifdef SELF_GRAVITY
  Gcons = nJ*cs2;
  grav_mean_rho = 1.0;
#ifndef SELF_GRAVITY_USING_FFT_DISK
  if(pG->Nx[2] >1) grav_mean_rho = 1.0;
#endif

/* Set gravity constant*/
  four_pi_G = 4.0*PI*Gcons;
#endif /* SELF_GRAVITY */

  B0 = cs/sqrt(beta);
#ifndef BAROTROPIC
  P0 = cs2/Gamma;
#endif

/* Ensure a different initial random seed for each process in an MPI calc. */
  ixs = pG->Disp[0];
  jxs = pG->Disp[1];
  kxs = pG->Disp[2];
  iseed = -1 - (ixs + pDomain->Nx[0]*(jxs + pDomain->Nx[1]*kxs));

  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];

/* initialize wavenumbers, given input number of waves per L */
  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  ky = nwy*2.0*PI;
  kx = nwx*2.0*PI;
  kxt = kx+qshear*Omega_0*ky*time0;

  pG->time=time0;

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,j,k,&x1,&x2,&x3);
      if (((i-pG->Disp[0]) == 58) && ((j-pG->Disp[1]) == 16)) 
        printf("i=%d j=%d k=%d x1=%e x2=%e\n",i,j,k,x1,x2);

      rd  = 1.0+amp*cos(kxt*x1+ky*x2);
      rvx = amp*kx/ky*sin(kxt*x1+ky*x2);
      rvy = amp*sin(kxt*x1+ky*x2);
      rvz = 0.0;
      rp  = cs2*(rd-1.0);

      rbx = amp*nwy*cos(kxt*(x1-0.5*pG->dx1)+ky*x2);
      rby = -amp*nwx*cos(kxt*x1+ky*(x2-0.5*pG->dx2));
      rbz = 0.0;

      pG->U[k][j][i].d  = rd;
      pG->U[k][j][i].M1 = rd*rvx;
      pG->U[k][j][i].M2 = rd*rvy;
#ifndef FARGO
      pG->U[k][j][i].M2 -= rd*(qshear*Omega_0*x1);
#endif
      pG->U[k][j][i].M3 = rd*rvz;
#ifdef ADIABATIC
      pG->U[k][j][i].E = (P0+rp)/Gamma_1
         + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) 
         + SQR(pG->U[k][j][i].M3))/rd;
#endif

#ifdef MHD
        pG->B1i[k][j][i] = rbx;
        pG->B2i[k][j][i] = B0+rby;
        pG->B3i[k][j][i] = 0.0;

        if (i==ie) cc_pos(pG,ie+1,j,k,&x1,&x2,&x3);
        rbx = amp*nwy*cos(kx*(x1-0.5*pG->dx1)+ky*x2);
        if (j==je) cc_pos(pG,i,je+1,k,&x1,&x2,&x3);
        rby = -amp*nwx*cos(kx*x1+ky*(x2-0.5*pG->dx2));
        if (i==ie) pG->B1i[k][j][ie+1] = rbx;
        if (j==je) pG->B2i[k][je+1][i] = B0+rby;
        if (pG->Nx[2] > 1 && k==ke) pG->B3i[ke+1][j][i] = 0.0;
#endif /* MHD */
    }
  }}
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].B1c = 0.5*(pG->B1i[k][j][i]+pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(pG->B2i[k][j][i]+pG->B2i[k][j+1][i]);
        if (pG->Nx[2] >1) pG->U[k][j][i].B3c = 0.5*(pG->B3i[k][j][i]+pG->B3i[k+1][j][i]); else pG->U[k][j][i].B3c =pG->B3i[k][j][i];
#ifdef ADIABATIC
        pG->U[k][j][i].E += 0.5*(SQR(pG->U[k][j][i].B1c)
         + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
#endif
      }
    }
  }
#endif /* MHD */

/* enroll gravitational potential function */

  ShearingBoxPot = UnstratifiedDisk;

/* enroll new history variables, only once with SMR  */

  dVol = pDomain->Nx[0]*pDomain->Nx[1]*pDomain->Nx[2];

/* history dump for linear perturbation amplitude. See Kim & Ostriker 2001 */
  dump_history_enroll(hst_sigma, "<sigma>");
  dump_history_enroll(hst_ux, "<ux>");
  dump_history_enroll(hst_uy, "<uy>");
#ifdef MHD
  dump_history_enroll(hst_m1, "<m1>");
  dump_history_enroll(hst_m2, "<m2>");
#endif

/* history dump for peturbed quantities at a specific grid point */
  dump_history_enroll(hst_dSigma, "<dSigma>");
  dump_history_enroll(hst_Vx, "<Vx>");
  dump_history_enroll(hst_dVy, "<dVy>");

#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_dBy, "<dBy>");
#endif /* MHD */
#ifdef SELF_GRAVITY
  dump_history_enroll(hst_Phi, "<Phi>");
  dump_history_enroll(hst_dPhi, "<dPhi>");
#endif
#ifdef ADIABATIC
  dump_history_enroll(hst_dE, "<dE>");
#endif
  
  printf("=== end of problem setting ===\n");
  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
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

/*
 * 'problem_read_restart' must enroll gravity on restarts
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  Real beta, Q, nJ, cs, cs2, Gcons;

/* Read problem parameters. */
  Omega_0 = par_getd("problem","omega");
  qshear = par_getd("problem","qshear");

/* Read parameters for magnetic field */
  beta = par_getd("problem","beta"); 

/* Read parameters for self gravity */
  Q=par_getd("problem","Q");
  nJ= par_getd("problem","nJ");

  cs=sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  cs2=SQR(cs);

#ifdef SELF_GRAVITY
  Gcons = nJ*cs2;
  grav_mean_rho = 1.0;
#ifndef SELF_GRAVITY_USING_FFT_DISK
  if(pM->Nx[2] >1) grav_mean_rho = 1.0;
#endif

/* Set gravity constant*/
  four_pi_G = 4.0*PI*Gcons;
#endif /* SELF_GRAVITY */

/* enroll gravitational potential function */

  ShearingBoxPot = UnstratifiedDisk;

/* history dump for linear perturbation amplitude. See Kim & Ostriker 2001 */
  dump_history_enroll(hst_sigma, "<sigma>");
  dump_history_enroll(hst_ux, "<ux>");
  dump_history_enroll(hst_uy, "<uy>");
#ifdef MHD
  dump_history_enroll(hst_m1, "<m1>");
  dump_history_enroll(hst_m2, "<m2>");
#endif


/* history dump for peturbed quantities at a specific grid point */
  dump_history_enroll(hst_dSigma, "<dSigma>");
  dump_history_enroll(hst_Vx, "<Vx>");
  dump_history_enroll(hst_dVy, "<dVy>");

#ifdef MHD
  dump_history_enroll(hst_Bx, "<Bx>");
  dump_history_enroll(hst_dBy, "<dBy>");
#endif /* MHD */
#ifdef SELF_GRAVITY
  dump_history_enroll(hst_Phi, "<Phi>");
  dump_history_enroll(hst_dPhi, "<dPhi>");
#endif

/* should be modified for SMR version */
  dVol = 1.0;
  if (pM->dx[0] > 0.0) dVol /= pM->dx[0];
  if (pM->dx[1] > 0.0) dVol /= pM->dx[1];
  if (pM->dx[2] > 0.0) dVol /= pM->dx[2];


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

/*==============================================================================
 * PRIVATE FUNCTION:
 *============================================================================*/
static Real UnstratifiedDisk(const Real x1, const Real x2, const Real x3)
{
  Real phi=0.0;
#ifndef FARGO
  phi -= qshear*Omega_0*Omega_0*x1*x1;
#endif
  return phi;
}

static Real hst_sigma(const GridS *pG, const int i, const int j, const int k)
{
  Real kx,kxt,ky;
  Real x1,x2,x3;
  int nwx,nwy;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  return (pG->U[k][j][i].d-1)/cos(kxt*x1+ky*x2);
}

static Real hst_ux(const GridS *pG, const int i, const int j, const int k)
{
  Real kx,kxt,ky;
  Real x1,x2,x3;
  int nwx,nwy;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  return (pG->U[k][j][i].M1/pG->U[k][j][i].d-1)/sin(kxt*x1+ky*x2);
}

static Real hst_uy(const GridS *pG, const int i, const int j, const int k)
{
  Real kx,kxt,ky;
  Real x1,x2,x3;
  int nwx,nwy;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  return (pG->U[k][j][i].M2/pG->U[k][j][i].d-1)/sin(kxt*x1+ky*x2);
}

#ifdef MHD
static Real hst_m1(const GridS *pG, const int i, const int j, const int k)
{
  Real kx,kxt,ky;
  Real x1,x2,x3;
  int nwx,nwy;

  Real nJ,Q,beta,cs,B0;

  nJ = par_getd("problem","nJ");
  Q = par_getd("problem","Q");
  beta = par_getd("problem","beta");

  cs = sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  B0 = cs/sqrt(beta);

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  return (pG->U[k][j][i].B1c/ky/B0)/cos(kxt*x1+ky*x2);
}

static Real hst_m2(const GridS *pG, const int i, const int j, const int k)
{
  Real kx,kxt,ky;
  Real x1,x2,x3;
  int nwx,nwy;

  Real nJ,Q,beta,cs,B0;

  nJ = par_getd("problem","nJ");
  Q = par_getd("problem","Q");
  beta = par_getd("problem","beta");

  cs = sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  B0 = cs/sqrt(beta);

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  return -(pG->U[k][j][i].B2c-B0)/kxt/B0/cos(kxt*x1+ky*x2);
}
#endif

static Real hst_dSigma(const GridS *pG, const int i, const int j, const int k)
{
  Real dSigma=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  if((ii == 58) && (jj == 16) && (kk == ks)) dSigma=pG->U[kk][jj][ii].d-1.0;
  return dSigma*dVol;
}

static Real hst_Vx(const GridS *pG, const int i, const int j, const int k)
{
  Real Vx=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  if(ii == 58 && jj == 16 && kk == ks) Vx=pG->U[kk][jj][ii].M1/pG->U[kk][jj][ii].d;
  return Vx*dVol;
}

static Real hst_dVy(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,dVy=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  if(ii == 58 && jj == 16 && kk == ks){
    cc_pos(pG,ii,jj,kk,&x1,&x2,&x3);
#ifndef FARGO
    dVy = (pG->U[kk][jj][ii].M2/pG->U[kk][jj][ii].d + qshear*Omega_0*x1);
#else
    dVy = (pG->U[kk][jj][ii].M2/pG->U[kk][jj][ii].d);
#endif
  }
  return dVy*dVol;
}

#ifdef MHD
static Real hst_Bx(const GridS *pG, const int i, const int j, const int k)
{
  Real B1c=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  if(ii == 58 && jj == 16 && kk == ks) B1c= pG->U[kk][jj][ii].B1c;
  return B1c*dVol;
}

static Real hst_dBy(const GridS *pG, const int i, const int j, const int k)
{
  Real nJ,Q,beta,cs,B0,dBy=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  nJ = par_getd("problem","nJ");
  Q = par_getd("problem","Q");
  beta = par_getd("problem","beta");

  cs = sqrt(4.0-2.0*qshear)/PI/nJ/Q;
  B0 = cs/sqrt(beta);

  if(ii == 58 && jj == 16 && kk == ks) dBy = pG->U[kk][jj][ii].B2c - B0;

  return dBy*dVol;
}
#endif

#ifdef SELF_GRAVITY
static Real hst_Phi(const GridS *pG, const int i, const int j, const int k)
{
  Real Phi=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;
  if(ii == 58 && jj == 16 && kk == ks) Phi=pG->Phi[kk][jj][ii];
  return Phi*dVol;
}


static Real hst_dPhi(const GridS *pG, const int i, const int j, const int k)
{
  Real nJ,Q,cs,dPhi,Phi1;
  Real kx,kxt,ky,k2,k20;
  Real x1,x2,x3,zmax;
  int nwx,nwy;
  int ks = pG->ks, ke = pG->ke;

  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  nJ = par_getd("problem","nJ");
  Q = par_getd("problem","Q");
  nwx = par_geti_def("problem","nwx",-6);
  nwy = par_geti_def("problem","nwy",1);
  zmax = par_getd("domain1","x3max");

  kx = nwx*2*PI;
  ky = nwy*2*PI;
  kxt = kx+qshear*ky*pG->time;
  k2 =kxt*kxt+ky*ky;
  k20 =kx*kx+ky*ky;

  cs = sqrt(4.0-2.0*qshear)/(PI*nJ*Q);
  Phi1 = -4.0*PI*nJ*cs*cs/k20*(pG->U[k][j][i].d-1.0);
#ifdef SELF_GRAVITY_USING_FFT_DISK
  Phi1 *= 1-0.5*(exp(-sqrt(k20)*(zmax-fabs(x3)))+exp(-sqrt(k20)*(zmax+fabs(x3))));
#endif

  dPhi = (pG->Phi[k][j][i]-Phi1)/(pG->U[k][j][i].d-1)*k2;
  dPhi /=  1-0.5*(exp(-sqrt(k2)*(zmax-fabs(x3)))+exp(-sqrt(k2)*(zmax+fabs(x3))));
  return dPhi;
}

#endif

#ifdef ADIABATIC
static Real hst_dE(const GridS *pG, const int i, const int j, const int k)
{
  Real nJ,Q,cs,cs2,E0,dE=0.;
  int ii=i-pG->Disp[0],jj=j-pG->Disp[1],kk=k-pG->Disp[2];
  int ks = pG->ks, ke = pG->ke;

  nJ = par_getd("problem","nJ");
  Q = par_getd("problem","Q");

  cs = sqrt(4.0-2.0*qshear)/(PI*nJ*Q);
  cs2 = SQR(cs);
  E0 = cs2/Gamma/Gamma_1;
  if(ii == 58 && jj == 16 && kk == ks) {
  dE = pG->U[kk][jj][ii].E 
       - 0.5*(SQR(pG->U[kk][jj][ii].M1) + SQR(pG->U[kk][jj][ii].M2) 
         + SQR(pG->U[kk][jj][ii].M3))/pG->U[kk][jj][ii].d;
#ifdef MHD
  dE -= 0.5*(SQR(pG->U[kk][jj][ii].B1c)
           + SQR(pG->U[kk][jj][ii].B2c) + SQR(pG->U[kk][jj][ii].B3c));
#endif
  dE -= E0;

  }
  return dE*dVol;
}

#endif
