#include "copyright.h"
/*============================================================================*/
/*! \file shkset1d.c 
 *  \brief Problem generator for 1-D Riemann problems.  
 *
 * PURPOSE: Problem generator for 1-D Riemann problems.  Initial discontinuity
 *   is located so there are equal numbers of cells to the left and right (at
 *   center of grid based on integer index).  Initializes plane-parallel shock
 *   along x1 (in 1D, 2D, 3D), along x2 (in 2D, 3D), and along x3 (in 3D).
 *
 * If error_test=1 in the <problem> block, then the L1 error in the final
 * solution will be computed for the Sod shocktube (hydrodynamics) or the RJ2a
 * test (MHD).  This is useful for regression tests.
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* Analytic solution at stopping time, shared with Userwork_after_loop to
 * compute L1 error */
static ConsS ***RootSoln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:    */

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i,il,iu,j,jl,ju,k,kl,ku;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;

  int shk_dir; /* Shock direction: {1,2,3} -> {x1,x2,x3} */
  Real ang_2, ang_3; /* Rotation angles about the y and z' axis */
  Real sin_a2, cos_a2, sin_a3, cos_a3;

  Real x1,x2,x3;
  Prim1DS Wl, Wr;
  Cons1DS U1d, Ul, Ur;
  Real Bxl=0.0, Bxr=0.0;
/* speeds of shock, contact, head and foot of rarefaction for Sod test */
/* speeds of slow/fast shocks, Alfven wave and contact in RJ2a test */
  Real tlim;
  int err_test;
  Real r,xs,xc,xf,xh,vs,vc,vf,vh;
  Real xfp,xrp,xsp,xsm,xrm,xfm,vfp,vrp,vsp,vsm,vrm,vfm;
  Real d0,v0,Mx,My,Mz,E0,r0,Bx,By,Bz;
#if (NSCALARS > 0)
  int n;
#endif

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }

/* Parse left state read from input file: dl,pl,ul,vl,wl,bxl,byl,bzl */

  Wl.d = par_getd("problem","dl");
#ifdef ADIABATIC
  Wl.P = par_getd("problem","pl");
#endif
  Wl.Vx = par_getd("problem","v1l");
  Wl.Vy = par_getd("problem","v2l");
  Wl.Vz = par_getd("problem","v3l");
#ifdef MHD
  Bxl = par_getd("problem","b1l");
  Wl.By = par_getd("problem","b2l");
  Wl.Bz = par_getd("problem","b3l");
#endif
#if (NSCALARS > 0)
  Wl.r[0] = par_getd("problem","r0l");
#endif

/* Parse right state read from input file: dr,pr,ur,vr,wr,bxr,byr,bzr */

  Wr.d = par_getd("problem","dr");
#ifdef ADIABATIC
  Wr.P = par_getd("problem","pr");
#endif
  Wr.Vx = par_getd("problem","v1r");
  Wr.Vy = par_getd("problem","v2r");
  Wr.Vz = par_getd("problem","v3r");
#ifdef MHD
  Bxr = par_getd("problem","b1r");
  Wr.By = par_getd("problem","b2r");
  Wr.Bz = par_getd("problem","b3r");
  if (Bxr != Bxl) ath_error(0,"[shkset1d] L/R values of Bx not the same\n");
#endif
#if (NSCALARS > 0)
  Wr.r[0] = par_getd("problem","r0r");
#endif

  Ul = Prim1D_to_Cons1D(&Wl, &Bxl);
  Ur = Prim1D_to_Cons1D(&Wr, &Bxr);

/* Parse shock direction */
  shk_dir = par_geti("problem","shk_dir");
  if (shk_dir != 1 && shk_dir != 2 && shk_dir != 3) {
    ath_error("[problem]: shk_dir = %d must be either 1,2 or 3\n",shk_dir);
  }

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

/* Initialize the grid including the ghost cells.  Discontinuity is always
 * located at x=0, so xmin/xmax in input file must be set appropriately. */

  switch(shk_dir) {
/*--- shock in 1-direction ---------------------------------------------------*/
  case 1:  /* shock in 1-direction  */
    ang_2 = 0.0;
    ang_3 = 0.0;

    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive and conserved variables to be L or R state */
          if (x1 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = U1d.Mx;
          pGrid->U[k][j][i].M2 = U1d.My;
          pGrid->U[k][j][i].M3 = U1d.Mz;
#ifdef MHD
          pGrid->B1i[k][j][i] = Bxl;
          pGrid->B2i[k][j][i] = U1d.By;
          pGrid->B3i[k][j][i] = U1d.Bz;
          pGrid->U[k][j][i].B1c = Bxl;
          pGrid->U[k][j][i].B2c = U1d.By;
          pGrid->U[k][j][i].B3c = U1d.Bz;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
    break;

/*--- shock in 2-direction ---------------------------------------------------*/
  case 2:  /* shock in 2-direction  */
    ang_2 = 0.0;
    ang_3 = PI/2.0;
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
          if (x2 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = -U1d.My;
          pGrid->U[k][j][i].M2 = U1d.Mx;
          pGrid->U[k][j][i].M3 = U1d.Mz;
#ifdef MHD
          pGrid->B1i[k][j][i] = -U1d.By;
          pGrid->B2i[k][j][i] = Bxl;
          pGrid->B3i[k][j][i] = U1d.Bz;
          pGrid->U[k][j][i].B1c = -U1d.By;
          pGrid->U[k][j][i].B2c = Bxl;
          pGrid->U[k][j][i].B3c = U1d.Bz;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
    break;

/*--- shock in 3-direction ---------------------------------------------------*/
  case 3:  /* shock in 3-direction  */
    ang_2 = PI/2.0;
    ang_3 = 0.0;
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          cc_pos(pGrid, i, j, k, &x1, &x2, &x3);

/* set primitive variables to be L or R state */
          if (x3 <= 0.0) {
            U1d = Ul;
          } else {
            U1d = Ur;
          }

/* Initialize conserved (and with SR the primitive) variables in Grid */
          pGrid->U[k][j][i].d  = U1d.d;
          pGrid->U[k][j][i].M1 = -U1d.Mz;
          pGrid->U[k][j][i].M2 = U1d.My;
          pGrid->U[k][j][i].M3 = U1d.Mx;
#ifdef MHD
          pGrid->B1i[k][j][i] = -U1d.Bz;
          pGrid->B2i[k][j][i] = U1d.By;
          pGrid->B3i[k][j][i] = Bxl;
          pGrid->U[k][j][i].B1c = -U1d.Bz;
          pGrid->U[k][j][i].B2c = U1d.By;
          pGrid->U[k][j][i].B3c = Bxl;
#endif
#ifdef ADIABATIC
          pGrid->U[k][j][i].E = U1d.E;
#endif
#if (NSCALARS > 0)
          pGrid->U[k][j][i].s[0] = U1d.s[0];
#endif
        }
      }
    }
  break;
  default:
    ath_error("[shkset1d]: invalid shk_dir = %i\n",shk_dir);
  }

/* Compute Analytic solution for Sod and RJ4a tests, if required */

  tlim = par_getd("time","tlim");
  err_test = par_getd_def("problem","error_test",0);
  if (err_test == 1) {

    sin_a3 = sin(ang_3);
    cos_a3 = cos(ang_3);
    sin_a2 = sin(ang_2);
    cos_a2 = cos(ang_2);

/* wave speeds for Sod test */
#ifdef HYDRO
    vs = 1.7522; xs = vs*tlim;
    vc = 0.92745; xc = vc*tlim;
    vf = -0.07027; xf = vf*tlim;
    vh = -1.1832; xh = vh*tlim;
#endif /* HYDRO */

/* wave speeds for RJ2a test */
#ifdef MHD
    vfp = 2.2638; xfp = vfp*tlim;
    vrp = (0.53432 + 1.0/sqrt(PI*1.309)); xrp = vrp*tlim;
    vsp = (0.53432 + 0.48144/1.309); xsp = vsp*tlim;
    vc = 0.57538; xc = vc*tlim;
    vsm = (0.60588 - 0.51594/1.4903); xsm = vsm*tlim;
    vrm = (0.60588 - 1.0/sqrt(PI*1.4903)); xrm = vrm*tlim;
    vfm = (1.2 - 2.3305/1.08); xfm = vfm*tlim;
#endif /* MHD */

    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        r = cos_a2*(x1*cos_a3 + x2*sin_a3) + x3*sin_a2;

/* Sod solution */
#ifdef HYDRO
        My = Mz = 0.0;
        if (r > xs) {
          d0 = 0.125;
          Mx = 0.0;
          E0 = 0.25;
          r0 = 0.0;
        } else if (r > xc) {
          d0 = 0.26557;
          Mx = 0.92745*d0;
          E0 = 0.87204;
          r0 = 0.0;
        } else if (r > xf) {
          d0 = 0.42632;
          Mx = 0.92745*d0;
          E0 = 0.94118;
          r0 = 1.0;
        } else if (r > xh) {
          v0 = 0.92745*(r-xh)/(xf-xh);
          d0 = 0.42632*pow((1.0+0.20046*(0.92745-v0)),5);
          E0 = (0.30313*pow((1.0+0.20046*(0.92745-v0)),7))/0.4 + 0.5*d0*v0*v0;
          r0 = 1.0;
          Mx = v0*d0;
        } else {
          d0 = 1.0;
          Mx = 0.0;
          E0 = 2.5;
          r0 = 1.0;
        }
#endif /* HYDRO */
/* RJ2a solution (Dai & Woodward 1994 Tables Ia and Ib) */
#ifdef MHD
        Bx = 2.0/sqrt(4.0*PI);
        if (r > xfp) {
          d0 = 1.0;
          Mx = 0.0;
          My = 0.0;
          Mz = 0.0;
          By = 4.0/sqrt(4.0*PI);
          Bz = 2.0/sqrt(4.0*PI);
          E0 = 1.0/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 0.0;
        } else if (r > xrp) {
          d0 = 1.3090;
          Mx = 0.53432*d0;
          My = -0.094572*d0;
          Mz = -0.047286*d0;
          By = 5.3452/sqrt(4.0*PI);
          Bz = 2.6726/sqrt(4.0*PI);
          E0 = 1.5844/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 0.0;
        } else if (r > xsp) {
          d0 = 1.3090;
          Mx = 0.53432*d0;
          My = -0.18411*d0;
          Mz = 0.17554*d0;
          By = 5.7083/sqrt(4.0*PI);
          Bz = 1.7689/sqrt(4.0*PI);
          E0 = 1.5844/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 0.0;
        } else if (r > xc) {
          d0 = 1.4735;
          Mx = 0.57538*d0;
          My = 0.047601*d0;
          Mz = 0.24734*d0;
          By = 5.0074/sqrt(4.0*PI);
          Bz = 1.5517/sqrt(4.0*PI);
          E0 = 1.9317/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 0.0;
        } else if (r > xsm) {
          d0 = 1.6343;
          Mx = 0.57538*d0;
          My = 0.047601*d0;
          Mz = 0.24734*d0;
          By = 5.0074/sqrt(4.0*PI);
          Bz = 1.5517/sqrt(4.0*PI);
          E0 = 1.9317/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 1.0;
        } else if (r > xrm) {
          d0 = 1.4903;
          Mx = 0.60588*d0;
          My = 0.22157*d0;
          Mz = 0.30125*d0;
          By = 5.5713/sqrt(4.0*PI);
          Bz = 1.7264/sqrt(4.0*PI);
          E0 = 1.6558/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 1.0;
        } else if (r > xfm) {
          d0 = 1.4903;
          Mx = 0.60588*d0;
          My = 0.11235*d0;
          Mz = 0.55686*d0;
          By = 5.0987/sqrt(4.0*PI);
          Bz = 2.8326/sqrt(4.0*PI);
          E0 = 1.6558/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 1.0;
        } else {
          d0 = 1.08;
          Mx = 1.2*d0;
          My = 0.01*d0;
          Mz = 0.5*d0;
          By = 3.6/sqrt(4.0*PI);
          Bz = 2.0/sqrt(4.0*PI);
          E0 = 0.95/Gamma_1 + 0.5*((Mx*Mx+My*My+Mz*Mz)/d0 + (Bx*Bx+By*By+Bz*Bz));
          r0 = 1.0;
        }
#endif /* MHD */
 
        RootSoln[k][j][i].d = d0;

        RootSoln[k][j][i].M1 = Mx*cos_a2*cos_a3 - My*sin_a3 - Mz*sin_a2*cos_a3;
        RootSoln[k][j][i].M2 = Mx*cos_a2*sin_a3 + My*cos_a3 - Mz*sin_a2*sin_a3;
        RootSoln[k][j][i].M3 = Mx*sin_a2                    + Mz*cos_a2;

#ifdef MHD
        RootSoln[k][j][i].B1c = Bx*cos_a2*cos_a3 - By*sin_a3 - Bz*sin_a2*cos_a3;
        RootSoln[k][j][i].B2c = Bx*cos_a2*sin_a3 + By*cos_a3 - Bz*sin_a2*sin_a3;
        RootSoln[k][j][i].B3c = Bx*sin_a2                    + Bz*cos_a2;
#endif /* MHD */

#ifndef ISOTHERMAL
        RootSoln[k][j][i].E = E0;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) RootSoln[k][j][i].s[n] = r0*d0;
#endif

      }
    }}

  } /* end calculation of analytic (root) solution */

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
  GridS *pGrid;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke;
  Real rms_error=0.0;
  ConsS error,total_error;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count, min_zones;
#if defined MPI_PARALLEL
  double err[8+(NSCALARS)], tot_err[8+(NSCALARS)];
  int ierr,myID;
#endif
#if (NSCALARS > 0)
  int n;
#endif

  int err_test = par_getd_def("problem","error_test",0);
  if (err_test == 0) return;

  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
#ifdef MHD
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
  total_error.E = 0.0;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = 0.0;
#endif

/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

/* compute L1 error in each variable, and rms total error */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
#ifdef MHD
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
    error.E = 0.0;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) error.s[n] = 0.0;
#endif

    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[k][j][i].d   - RootSoln[k][j][i].d );
      error.M1  += fabs(pGrid->U[k][j][i].M1  - RootSoln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - RootSoln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
#ifdef MHD
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[k][j][i].E   - RootSoln[k][j][i].E );
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        error.s[n] += fabs(pGrid->U[k][j][i].s[n] - RootSoln[k][j][i].s[n]);
#endif
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
#ifdef MHD
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E += error.E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) total_error.s[n] += error.s[n];
#endif
  }}

#ifdef MPI_PARALLEL
  Nx1 = pM->Domain[0][0].Nx[0];
  Nx2 = pM->Domain[0][0].Nx[1];
  Nx3 = pM->Domain[0][0].Nx[2];
#else
  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;
#endif
  count = Nx1*Nx2*Nx3;

#ifdef MPI_PARALLEL 
/* Now we have to use an All_Reduce to get the total error over all the MPI
 * grids.  Begin by copying the error into the err[] array */

  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
#ifdef MHD
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) err[8+n] = total_error.s[n];
#endif

  /* Sum up the Computed Error */
  ierr = MPI_Reduce(err,tot_err,(8+(NSCALARS)),MPI_DOUBLE,MPI_SUM,0,
    pM->Domain[0][0].Comm_Domain);

/* If I'm the parent, copy the sum back to the total_error variable */

  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
#ifdef MHD
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) total_error.s[n] = tot_err.s[8+n];
#endif
  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
#ifdef MHD
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) rms_error += SQR(total_error.s[n]);
#endif
  rms_error = sqrt(rms_error)/(double)count;

/* Print warning to stdout if rms_error exceeds estimate of 1st-order conv */
/* For 1D, assume shock propagates along direction with MAX number of zones */

  min_zones = Nx1;
  if (Nx2 > 1) min_zones = MAX(min_zones,Nx2);
  if (Nx3 > 1) min_zones = MAX(min_zones,Nx3);
  if (rms_error > 8.0/min_zones)
    printf("WARNING: rms_error=%e exceeds estimate\n",rms_error);

/* Print error to file "LinWave-errors.#.dat", where #=wave_flag  */

#ifdef MPI_PARALLEL
  fname = "../shock-errors.dat";
#else
  fname = "shock-errors.dat";
#endif

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
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
#ifdef MHD
    fprintf(fp,"  B1c  B2c  B3c");
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
      fprintf(fp,"  S[ %d ]",n);
    }
#endif
    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double)count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count) );

#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(total_error.E/(double)count) );
#endif /* ISOTHERMAL */

#ifdef MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    fprintf(fp,"  %e",total_error.s[n]/(double)count);
  }
#endif

  fprintf(fp,"\n");

  fclose(fp);

  return;
}
