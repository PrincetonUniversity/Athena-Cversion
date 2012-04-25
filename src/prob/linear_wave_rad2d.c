#include "copyright.h"
/*==============================================================================
 * FILE: linear_wave_rad2d.c
 *
 * PURPOSE: Problem generator for plane-parallel, non grid-aligned linear wave
 *   tests with the radiative transfer module in 2D. 
 *
 *   Configure --with-problem=linear_wave_rad2d --enable-radiation-transfer --with-gas=hydro
 *
 *   Can be used for either standing (problem/vflow=1.0) or travelling
 *   (problem/vflow=0.0) waves.
 *
 * USERWORK_AFTER_LOOP function computes L1 error norm in solution by computing
 * evolution assuming the initial state is an eigenvector.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include <complex.h>

Real kappa, B0, E0, T0, vph, rdamp, Bo, tau, amp, vflow, V0R, V0I, E0R, E0I;
Real sin_a, cos_a,lambda;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * const_B()
 * const_eps()
 * const_chi()
 *============================================================================*/
static void gauleg(Real x1, Real x2,  Real *x, Real *w, int n);
static void zroots(double complex *a, int m, double complex *roots, int polish);
static void laguer(double complex *ad, int m, double complex *x, int *its);
static Real grey_B(const GridS *pG, const int ifr, const int i,
		   const int j, const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);
static Real const_chi(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);
static void acoustic_wave_rad(Real Bo, Real tau, Real cs, Real d0, Real *vph, Real *rdamp,
                       Real *V0R, Real *V0I, Real *E0R, Real *E0I);
/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  RadGridS *pRG=(pDomain->RadGrid);
  ConsS ***Soln;
  int i=0,j=0,k=0,ifr;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3;
  Real d0,p0,u0,v0,w0;
  Real x1,x2,x3,r;
  int nang = pRG->nang, nf=pRG->nf, noct= pRG->noct;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  FILE *fp;
  char *fname;
  Real angle,x1size,x2size;

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

/* Initialization of ghost zones needed for first call to formal solution */
  is -= nghost;
  ie += nghost;
  if (pGrid->Nx[1] > 1) {
    js -= nghost;
    je += nghost;
  }
  if (pGrid->Nx[2] > 1) {
    ks -= nghost;
    ke += nghost;
  }

/* Read initial conditions  */

  amp = par_getd("problem","amp");
  vflow = par_getd("problem","vflow");
  Bo = par_getd("problem","Bo");
  tau = par_getd("problem","tau");

  x1size = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  x2size = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  angle = atan(x1size/x2size);
  sin_a = sin(angle);
  cos_a = cos(angle);

/* Use the larger angle to determine the wavelength */
  if (cos_a >= sin_a) {
    lambda = x1size*cos_a;
  } else {
    lambda = x2size*sin_a;
  }

/* Get eigenmatrix, where the quantities u0 and bx0 are parallel to the
 * wavevector, and v0,w0 are perpendicular. */

  R_ideal = 1.0;
  d0 = 1.0;
  p0 = d0/Gamma;  /*  c_s=1 */
  u0 = vflow*sqrt(Gamma*p0/d0);
  v0 = 0.0;
  w0 = 0.0;

/* get initial conditions for damped radiating wave */
  acoustic_wave_rad(Bo,tau,1.0,d0,&vph,&rdamp,&V0R,&V0I,&E0R,&E0I);
  printf("E, V: %g %g %g %g\n",E0R,E0I,V0R,V0I);

/* Now initialize 1D solution vector  */
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {

/* Set background state */

    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
    r = (x1*cos_a + x2*sin_a)/lambda;

    pGrid->U[k][j][i].d = d0;
    pGrid->U[k][j][i].E = p0/Gamma_1 + 0.5*d0*u0*u0;


    pGrid->U[k][j][i].d += amp*sin(2.0*PI*r);
    pGrid->U[k][j][i].E += amp*(sin(2.0*PI*r)*E0R - cos(2.0*PI*r)*E0I);
    pGrid->U[k][j][i].M1 = cos_a*d0*vflow;
    pGrid->U[k][j][i].M2 = sin_a*d0*vflow;
    pGrid->U[k][j][i].M3 = 0.0;
    pGrid->U[k][j][i].M1 += cos_a*amp*(sin(2.0*PI*r)*V0R - cos(2.0*PI*r)*V0I);
    pGrid->U[k][j][i].M2 += sin_a*amp*(sin(2.0*PI*r)*V0R - cos(2.0*PI*r)*V0I);
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pGrid->U[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*r));
#endif

  }}}



/* ------- Initialize boundary emission ---------------------------------- */
/* Compute Planck function given Bo, tau */
  E0 = 1.0 / (Gamma*Gamma_1);
  T0 = E0 * Gamma_1 / (d0 * R_ideal);
  kappa = tau * 2.0 * PI;
  B0 = Gamma * E0 / (Bo*PI);
  B00 = B0;
  printf("E0, B0 kappa: %g %g %g\n",E0,B0,kappa);

  il = pRG->is-1, iu = pRG->ie+1;
  jl = pRG->js,   ju = pRG->je;
  kl = pRG->ks,   ku = pRG->ke;
  if (pRG->Nx[1] > 1) { jl -= 1; ju += 1; }
  if (pRG->Nx[2] > 1) { kl -= 1; ku += 1; }

  for(ifr=0; ifr<nf; ifr++) {
    for(k=kl; k<=ku; k++) {
      for(j=jl; j<=ju; j++) {
	for(m=0; m<nang; m++) {
	  pRG->Ghstl1i[ifr][k][j][0][m] = B0;
	  pRG->Ghstl1i[ifr][k][j][2][m] = B0;
	  if (noct == 8) {
	    pRG->Ghstl1i[ifr][k][j][4][m] = B0;
	    pRG->Ghstl1i[ifr][k][j][6][m] = B0;
	  }	  
	  pRG->Ghstr1i[ifr][k][j][1][m] = B0;
	  pRG->Ghstr1i[ifr][k][j][3][m] = B0;
	  if (noct == 8) {
	    pRG->Ghstr1i[ifr][k][j][5][m] = B0;
	    pRG->Ghstr1i[ifr][k][j][7][m] = B0;
	  }
	}}}

    for(k=kl; k<=ku; k++) {
      for(i=il; i<=iu; i++) {
	for(m=0; m<nang; m++) {
	  pRG->Ghstl2i[ifr][k][i][0][m] = B0;
	  pRG->Ghstl2i[ifr][k][i][1][m] = B0;
	  if (noct == 8) {
	    pRG->Ghstl2i[ifr][k][i][4][m] = B0;
	    pRG->Ghstl2i[ifr][k][i][5][m] = B0;
	  }	  
	  pRG->Ghstr2i[ifr][k][i][2][m] = B0;
	  pRG->Ghstr2i[ifr][k][i][3][m] = B0;
	  if (noct == 8) {
	    pRG->Ghstr2i[ifr][k][i][6][m] = B0;
	    pRG->Ghstr2i[ifr][k][i][7][m] = B0;
	  }
	}}}   

    if (noct == 8) {
      for(j=jl; j<=ju; j++) {
	for(i=il; i<=iu; i++) {
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl3i[ifr][j][i][0][m] = B0;
	    pRG->Ghstl3i[ifr][j][i][1][m] = B0;
	    pRG->Ghstl3i[ifr][j][i][2][m] = B0;
	    pRG->Ghstl3i[ifr][j][i][3][m] = B0;
	    pRG->Ghstr3i[ifr][j][i][4][m] = B0;
	    pRG->Ghstr3i[ifr][j][i][5][m] = B0;
	    pRG->Ghstr3i[ifr][j][i][6][m] = B0;
	    pRG->Ghstr3i[ifr][j][i][7][m] = B0;
	  }}}
    }
   }

/* enroll radiation specification functions */
get_thermal_source = grey_B;
get_thermal_fraction = const_eps;
get_total_opacity = const_chi;

/* --------------computed formal solution ---------------------------- */
/* Called to intitialze radiation variables for first set of outputs,  */
/* not needed for evolution. */
  hydro_to_rad(pDomain);
  formal_solution(pDomain);

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
 * Userwork_in_formal_solution  - problem specific work in formal solution loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
#ifdef SELF_GRAVITY
  Real d0 = 1.0;
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = d0;
#endif /* SELF_GRAVITY */
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
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

void Userwork_in_formal_solution(DomainS *pD)
{
  return;
}

void Userwork_in_loop(MeshS *pM)
{
}

/*------------------------------------------------------------------------------
 * Userwork_after_loop: computes L1-error in linear waves,
 * Assumes intitial state is eigenvector.  Does not work
 * with MHD.
 */

void Userwork_after_loop(MeshS *pM)
{
  GridS *pGrid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
#if (NSCALARS > 0)
   int n;
#endif

  Real rms_error=0.0;
  ConsS error,total_error;
  ConsS ***Soln;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count;
  Real time, amp_t;
  Real d0,p0,u0,v0,w0;
  Real x1,x2,x3,r;
  int nx1,nx2,nx3;

  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
  total_error.E = 0.0;
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = 0.0;
#endif

/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

  time = pGrid->time;
  d0 = 1.0;
  p0 = d0/Gamma;  /*  c_s=1 */
  u0 = vflow*sqrt(Gamma*p0/d0);
  amp_t = amp * exp(-time*rdamp);
  printf("%14.10g %14.10g\n",rdamp,vph);

/* compute L1 error in each variable, and rms total error */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))==NULL)
    ath_error("[userwork_after_loop]: Error allocating memory for Soln\n");

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
    error.E = 0.0;
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) error.s[n] = 0.0;
#endif

    for (i=is; i<=ie; i++) {


      cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
      //x1 -= vph*time;
      //x2 -= vph*time;
      r = (x1*cos_a + x2*sin_a)/lambda;

      Soln[k][j][i].d = d0;
      Soln[k][j][i].E = p0/Gamma_1 + 0.5*d0*u0*u0;
      
/* Select appropriate solution based on direction of wavevector */

      //amp_t=amp;
      r -= vph*time;
      //printf("v amp: %g %g\n",vph*time,amp_t);
      Soln[k][j][i].d += amp_t*sin(2.0*PI*r);
      Soln[k][j][i].E += amp_t*(sin(2.0*PI*r)*E0R - cos(2.0*PI*r)*E0I);
      Soln[k][j][i].M1 = cos_a*d0*vflow;
      Soln[k][j][i].M2 = sin_a*d0*vflow;
      Soln[k][j][i].M3 = 0.0;
      Soln[k][j][i].M1 += cos_a*amp_t*(sin(2.0*PI*r)*V0R - cos(2.0*PI*r)*V0I);
      Soln[k][j][i].M2 += sin_a*amp_t*(sin(2.0*PI*r)*V0R - cos(2.0*PI*r)*V0I);
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
	Soln[k][j][i].s[n] = amp_t*(1.0 + sin(2.0*PI*r));
#endif

      error.d   += fabs(pGrid->U[k][j][i].d   - Soln[k][j][i].d );
      error.M1  += fabs(pGrid->U[k][j][i].M1  - Soln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - Soln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - Soln[k][j][i].M3); 
      error.E   += fabs(pGrid->U[k][j][i].E   - Soln[k][j][i].E );
      if(i==is) printf("%d %g %g %g\n",j,r,pGrid->U[k][j][i].d-1.0,Soln[k][j][i].d-1.0);
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) 
        error.s[n] += fabs(pGrid->U[k][j][i].s[n] - Soln[k][j][i].s[n]);;
#endif
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
    total_error.E += error.E;
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) total_error.s[n] += error.s[n];
#endif
  }}

  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;

  count = Nx1*Nx2*Nx3;

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3) + SQR(total_error.E);
  rms_error = sqrt(rms_error)/(double)count;

/* Print error to file "LinWave-errors.#.dat", where #=wave_flag  */

  fname = ath_fname(NULL,"LinWave-errors",NULL,NULL,1,0,NULL,"dat");
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
    fprintf(fp,"  E");
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  S[ %d ]",n);
    }
#endif
    fprintf(fp,"\n#\n");
  }
 
  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double) count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count));
  fprintf(fp,"  %e",(total_error.E/(double)count));
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  %e",total_error.s[n]/(double)count);
    }
#endif
  fprintf(fp,"\n");

  fclose(fp);

  free_3d_array(Soln);

  return;
}

static Real grey_B(const GridS *pG, const int ifr, const int i,
		   const int j, const int k)
{
  return B0 * pow(pG->tgas[k][j][i]/T0,4);
}

static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k)
{
  return 1.0;
  
}

static Real const_chi(const GridS *pG, const int ifr, const int i,
		      const int j, const int k)
{
  return kappa; 
}

static void acoustic_wave_rad_old(Real Bo, Real tau, Real cs, Real d0, Real *vph, Real *rdamp,
                       Real *V0R, Real *V0I, Real *E0R, Real *E0I) 
{
  
  double complex omega, delta, delta13, delta23;
  double complex omg, V0, E0;
  Real tau2, tau3, tau4;
  Real Gamma2, Gamma3;
  Real Bo2, Bo4;
  Real phi, phi2, phi4, xi;

  Bo2=Bo*Bo;
  Bo4=Bo2*Bo2;
  tau2=tau*tau;
  tau3=tau*tau2;
  tau4=tau2*tau2;
  Gamma2=Gamma*Gamma;
  Gamma3=Gamma*Gamma2;

  phi=1.0+3.0*tau2;
  phi2=phi*phi;
  phi4=phi2*phi2;
  xi=18.0*Gamma+Gamma2-27.0;

  delta = (72.0*(Gamma - 3.0)*tau*Bo2*phi2 - 4096.0*Gamma3*tau3) * I + 3.0*sqrt(3.0) *
          csqrt(Bo2*phi2*(64.0*Bo2*xi*tau2*phi2 - 65536.0*Gamma3*tau4 - Bo4*phi4));
  delta13 = cpow(delta,ONE_3RD);
  delta23 = cpow(delta,TWO_3RDS);

  omega = (3.0*Bo2 + 18.0*Bo2*tau2 - 256.0*Gamma2*tau2 + 27.0*Bo2*tau4 + 
	   I*16.0*Gamma*tau*delta13 + delta23) / (3.0*Bo*phi*delta13);
  
  omg = fabs(creal(omega)) + I * cimag(omega);
  V0 = cs * omega / d0;
  E0 = cs * cs * omega * omega / Gamma_1;

  (*V0R) = creal(V0); (*V0I) = cimag(V0);
  (*E0R) = creal(E0); (*E0I) = cimag(E0);
  (*vph) = creal(omg);
  (*rdamp) = 2.0 * PI * cimag(omg);

  return;
}

static void acoustic_wave_rad(Real Bo, Real tau, Real cs, Real d0, Real *vph, Real *rdamp,
                       Real *V0R, Real *V0I, Real *E0R, Real *E0I) 
{
  
  double complex *roots = NULL, *coeff = NULL;
  double complex omega, V0, E0;
  Real theta, mu;
  int i;
  
  if ((coeff = (double complex*)calloc_1d_array(4,sizeof(double complex)))==NULL)
    ath_error("[problem]: Error allocating memory for coeff\n");

  if ((roots = (double complex*)calloc_1d_array(3,sizeof(double complex)))==NULL)
    ath_error("[problem]: Error allocating memory for roots\n");

  mu = 1.0 - tau * atan(1.0 / tau);
  theta = 16.0 * tau * mu / Bo;
  coeff[0] = I * theta;
  coeff[1] = -1.0;
  coeff[2] = -I * theta * Gamma;
  coeff[3] = 1.0;
  zroots(coeff,3,roots,1);
  for(i=0; i<3; i++) 
    if(creal(roots[i]) > 0.5) {
      omega = roots[i];
    }

  V0 = cs * omega / d0;
  E0 = cs * cs * omega * omega / Gamma_1;
  (*V0R) = creal(V0); (*V0I) = cimag(V0);
  (*E0R) = creal(E0); (*E0I) = cimag(E0);
  (*vph) = creal(omega);
  (*rdamp) = 2.0 * PI * cimag(omega);

  free_1d_array(coeff);
  free_1d_array(roots);

  return;
}

static void zroots(double complex *a, int m, double complex *roots, int polish)
{

  double complex ad[4],x,b,c;
  Real eps = 1.0E-8;
  int i,j,jj,its;

  for(i=0; i<=m; i++) 
    ad[i] = a[i];

  for(j=m; j>0; j--) {
    x = 0.0 + I * 0.0;
    laguer(ad,j,&x,&its);
    if(fabs(cimag(x)) <= 2.0*eps*fabs(creal(x))) x = creal(x);
    roots[j-1] = x;
    b = ad[j];
    for(jj=j-1; jj >=0; jj--) {
      c = ad[jj];
      ad[jj] = b;
      b = x * b + c;
    }
  }
  if (polish == 1) 
    for(j=0; j<m; j++) {
      laguer(a,m,&(roots[j]),&its);
    }
  return;
}

static void laguer(double complex *a, int m, double complex *x, int *its)
{
  Real epss = 2.0e-7;
  int iter, j;
  Real abx,abp,abm,err;
  double complex dx,x1,b,d,f,g,h,sq,gp,gm,g2;
  static Real frac[9] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};

  for(iter=1; iter <=18; iter++) {
    (*its) = iter;
    b = a[m];
    err = cabs(b);
    d = f = 0.0 + I * 0.0;
    abx = cabs(*x);
    for(j=m-1; j>=0; j--) {
      f = (*x) * f + d;
      d = (*x) * d + b;
      b = (*x) * b + a[j];
      err = cabs(b) + abx * err;
    }
    err *= epss;
    if (cabs(b) <= err) return;
    g = d / b;
    g2 = g * g;
    h = g2 - 2.0 * f / b;
    sq = csqrt((m-1) * (m * h - g2));
    gp = g + sq;
    gm = g - sq;
    abp = cabs(gp);
    abm = cabs(gm);
    if (abp < abm) gp = gm;
    if (MAX(abp,abm) > 0.0) 
      dx = m / gp;
    else 
      dx = cexp( log(1.0+abx) + I * (float)iter);
    x1 = (*x) - dx;
    if ((*x) == x1) return;
    if (iter % 10) 
      (*x) = x1;
    else
      (*x) -= frac[iter/10] * dx;
  }
  printf("too many laguer iterations\n");
  return;
} 
