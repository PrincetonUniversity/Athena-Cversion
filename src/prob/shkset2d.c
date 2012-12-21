#include "copyright.h"
/*============================================================================*/
/*! \file shkset2d.c
 *  \brief Sets up shock at angle to grid to test multidimensional algorithm.
 *
 * PURPOSE: Sets up shock at ang_3=atan(Ly/Lx) to grid to test multidimensional
 *   algorithm.  Nx1/Nx2 must be the same ratio as Lx/Ly.  Uses the angle of the
 *   shock to remap ghost cells to the equivalent active grid cells, which
 *   requires Nx1>32, using special function pointers.  The shock is initialized
 *   with reference to a coordinate system (x,y,z) with transformation rules to
 *   the code coordinate system (x1,x2,x3)
 *   -  x =  x1*cos(ang_3) + x2*sin(ang_3)
 *   -  y = -x1*sin(ang_3) + x2*cos(ang_3)
 *   -  z = x3

 *   This inverts to:
 *   -  x1 = x*cos(ang_3) - y*sin(ang_3)
 *   -  x2 = x*sin(ang_3) + y*cos(ang_3)
 *   -  x3 = z								  
 *
 * If error_test=1 in the <problem> block, then the L1 error in the final
 * solution will be computed for the Sod shocktube (hydrodynamics) or the RJ2a
 * test (MHD).  This is useful for regression tests.
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * - shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * - shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * - shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.		      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * shkset2d_iib() - sets BCs on L-x1 (left edge) of grid.
 * shkset2d_oib() - sets BCs on R-x1 (right edge) of grid.
 * shkset2d_ijb() - sets BCs on L-x2 (bottom edge) of grid.
 * shkset2d_ojb() - sets BCs on R-x2 (top edge) of grid.
 *============================================================================*/

void shkset2d_iib(GridS *pGrid);
void shkset2d_oib(GridS *pGrid);
void shkset2d_ijb(GridS *pGrid);
void shkset2d_ojb(GridS *pGrid);

/* Make size of box and dimension of unit cell (r1 x r2) static globals so they
 * can be accessed by boundary value functions */
static Real Lx,Ly;
static int r1,r2;

/* Analytic solution at stopping time, shared with Userwork_after_loop to
 * compute L1 error */
static ConsS ***RootSoln=NULL;

/*----------------------------------------------------------------------------*/
/* problem:   */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  int kl,ku,irefine,ir,ix1,ix2,nx1,nx2,nx3,gcd;
  Real ang_2, ang_3; /* Rotation angles about the y and z' axis */
  Real sin_a2, cos_a2, sin_a3, cos_a3;

/* speeds of shock, contact, head and foot of rarefaction for Sod test */
/* speeds of slow/fast shocks, Alfven wave and contact in RJ2a test */
  Real tlim;
  int err_test;
  Real x1,x2,x3,r,xs,xc,xf,xh,vs,vc,vf,vh;
  Real xfp,xrp,xsp,xsm,xrm,xfm,vfp,vrp,vsp,vsm,vrm,vfm;
  Real d0,v0,Mx,My,Mz,E0,r0,Bx,By,Bz;

  Real rootdx1, rootdx2;
  Prim1DS Wl, Wr;
  Cons1DS Ul, Ur;
  ConsS ql, qr;
  Real Bxl=0.0,Bxr=0.0;
  div_t id;   /* structure containing remainder and quotient */

/* Following are used to compute volume of cell crossed by initial interface
 * that is assigned to left/right states */
  int dll, dlr, drr, drl;
  Real afl_lx, afr_lx, afl_rx, afr_rx;
  Real afl_ly, afr_ly, afl_ry, afr_ry;
  Real vfl, vfr, B1r, B2r;

  if (pGrid->Nx[1] == 1)
    ath_error("[shkset2d]: This problem can only be run in 2D or 3D\n");

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;

  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }

/* Find number of cells on root grid */

  irefine = 1;
  for (ir=1;ir<=pDomain->Level;ir++) irefine *= 2;

  Lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  Ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];

  rootdx1 = pGrid->dx1*((double)(irefine)); 
  rootdx2 = pGrid->dx2*((double)(irefine)); 

  nx1 = (int)(Lx/rootdx1);
  nx2 = (int)(Ly/rootdx2);

/* Compute greatest common divisor of nx1,nx2.  The size of the "unit cell"
 * is nx1/gcd by nx2/gcd */

  if((gcd = ath_gcd(nx1,nx2)) < 10)
    ath_error("[shkset2d]: Greatest Common Divisor (nx1,nx2) = %d\n",gcd);

  id = div(nx1,gcd);
  r1 = id.quot;
  if(id.rem != 0)
    ath_error("[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx1,gcd,id.rem);

  id = div(nx2,gcd);
  r2 = id.quot;
  if(id.rem != 0)
    ath_error("[shkset2d]: GCD failed, Remainder of %d / %d is %d\n",
	      nx2,gcd,id.rem);

  ath_pout(1,"The unit cell is (%d,1,%d) grid cells in size\n",r1,r2);

/* Compute angles initial interface makes to the grid */

  ang_2 = 0.0;
  cos_a2=cos(ang_2);
  sin_a2=sin(ang_2);
  if(Lx == Ly){
    cos_a3 = sin_a3 = sqrt(0.5);
  }
  else{
    ang_3 = atan((double)(Lx/Ly));
    sin_a3 = sin(ang_3);
    cos_a3 = cos(ang_3);
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
  if (Bxr != Bxl) ath_error(0,"[shkset2d] L/R values of Bx not the same\n");
#endif

  Ul = Prim1D_to_Cons1D(&Wl, &Bxl);
  Ur = Prim1D_to_Cons1D(&Wr, &Bxr);

/* Initialize ql rotated to the (x1,x2,x3) coordinate system */
  ql.d   = Ul.d;
  ql.M1  = Ul.Mx*cos_a3 - Ul.My*sin_a3;
  ql.M2  = Ul.Mx*sin_a3 + Ul.My*cos_a3;
  ql.M3  = Ul.Mz;
#ifdef MHD
  ql.B1c = Bxl*cos_a3 - Ul.By*sin_a3;
  ql.B2c = Bxl*sin_a3 + Ul.By*cos_a3;
  ql.B3c = Ul.Bz;
#endif
#ifdef ADIABATIC
  ql.E   = Ul.E;
#endif

/* Initialize qr rotated to the (x1,x2,x3) coordinate system */
  qr.d   = Ur.d;
  qr.M1  = Ur.Mx*cos_a3 - Ur.My*sin_a3;
  qr.M2  = Ur.Mx*sin_a3 + Ur.My*cos_a3;
  qr.M3  = Ur.Mz;
#ifdef MHD
  qr.B1c = Bxr*cos_a3 - Ur.By*sin_a3;
  qr.B2c = Bxr*sin_a3 + Ur.By*cos_a3;
  qr.B3c = Ur.Bz;
#endif
#ifdef ADIABATIC
  qr.E   = Ur.E;
#endif

/* Initialize the grid */

  for (k=kl; k<=ku; k++) {
    for (j=0; j<=je+nghost; j++) {
      ix2 = j + pGrid->Disp[1];
      for (i=0; i<=ie+nghost; i++) {
	ix1 = i + pGrid->Disp[0];

/* cell is completely in the left state */
	if((drr = r2*(ix1) + r1*(ix2) - gcd*r1*r2) <= 0){
	  pGrid->U[k][j][i] = ql;
#ifdef MHD
	  pGrid->B1i[k][j][i] = ql.B1c;
	  pGrid->B2i[k][j][i] = ql.B2c;
	  pGrid->B3i[k][j][i] = ql.B3c;
#endif /* MHD */
	}
/* cell is completely in the right state */
	else if((dll = r2*(ix1-1) + r1*(ix2-1) - gcd*r1*r2) >= 0){
	  pGrid->U[k][j][i] = qr;
#ifdef MHD
	  pGrid->B1i[k][j][i] = qr.B1c;
	  pGrid->B2i[k][j][i] = qr.B2c;
	  pGrid->B3i[k][j][i] = qr.B3c;
#endif /* MHD */
	}
/* The more complicated case of a cell  split by the interface boundary */
	else{
	  dlr = r2*(ix1-1) + r1*(ix2) - gcd*r1*r2;

	  if(dlr < 0){ /* The boundary hits the right y-face */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = (Real)(-dlr)/(Real)(r2);
	    afr_ry = 1.0 - afl_ry;
	  }
	  else if(dlr > 0){ /* The boundary hits the left x-face */
	    afl_lx = (Real)(-dll)/(Real)(r1);
	    afr_lx = 1.0 - afl_lx;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }
	  else{ /* dlr == 0.0, The boundary hits the grid cell corner */
	    afl_lx = 1.0;
	    afr_lx = 0.0;
	    afl_ry = 0.0;
	    afr_ry = 1.0;
	  }

	  drl = r2*(ix1) + r1*(ix2-1) - gcd*r1*r2;

	  if(drl < 0){ /* The boundary hits the right x-face */
	    afl_rx = (Real)(-drl)/(Real)(r1);
	    afr_rx = 1.0 - afl_rx;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }
	  else if(drl > 0){ /* The boundary hits the left y-face */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = (Real)(-dll)/(Real)(r2);
	    afr_ly = 1.0 - afl_ly;
	  }
	  else{ /* drl == 0.0, The boundary hits the grid cell corner */
	    afl_rx = 0.0;
	    afr_rx = 1.0;
	    afl_ly = 1.0;
	    afr_ly = 0.0;
	  }

/* The boundary hits both x-interfaces */
	  if(dlr > 0 && drl < 0){ 
	    vfl = 0.5*(afl_lx + afl_rx);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both y-interfaces */
	  else if(dlr < 0 && drl > 0){ 
	    vfl = 0.5*(afl_ly + afl_ry);
	    vfr = 1.0 - vfl;
	  }
/* The boundary hits both grid cell corners */
	  else if(dlr == 0 && drl == 0){ 
	    vfl = vfr = 0.5;
	  }
/* The boundary hits the left x- and left y-interface */
	  else if(dlr > 0 && drl > 0){
	    vfl = 0.5*afl_lx*afl_ly;
	    vfr = 1.0 - vfl;
	  }
/* dlr<0 && drl<0:  The boundary hits the right x- and right y-interface */
	  else{ 
	    vfr = 0.5*afr_rx*afr_ry;
	    vfl = 1.0 - vfr;
	  }

/* Initialize the x- and y-interface magnetic fields */
#ifdef MHD
	  pGrid->B1i[k][j][i] = afl_lx*ql.B1c + afr_lx*qr.B1c;
	  B1r              = afl_rx*ql.B1c + afr_rx*qr.B1c;

	  pGrid->B2i[k][j][i] = afl_ly*ql.B2c + afr_ly*qr.B2c;
	  B2r              = afl_ry*ql.B2c + afr_ry*qr.B2c;

	  pGrid->B3i[k][j][i] = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */

/* Initialize the volume averaged quantities */
	  pGrid->U[k][j][i].d  = vfl*ql.d + vfr*qr.d;
	  pGrid->U[k][j][i].M1 = vfl*ql.M1 + vfr*qr.M1;
	  pGrid->U[k][j][i].M2 = vfl*ql.M2 + vfr*qr.M2;
	  pGrid->U[k][j][i].M3 = vfl*ql.M3 + vfr*qr.M3;
#ifdef MHD
	  pGrid->U[k][j][i].B1c = 0.5*(pGrid->B1i[k][j][i] + B1r);
	  pGrid->U[k][j][i].B2c = 0.5*(pGrid->B2i[k][j][i] + B2r);
	  pGrid->U[k][j][i].B3c = vfl*ql.B3c + vfr*qr.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
	  pGrid->U[k][j][i].E  = vfl*ql.E + vfr*qr.E;
#endif
	}
      }
    }
  }

/* Set boundary value function pointers */

  bvals_mhd_fun(pDomain, left_x1,  shkset2d_iib);
  bvals_mhd_fun(pDomain, left_x2,  shkset2d_ijb);
  bvals_mhd_fun(pDomain, right_x1, shkset2d_oib);
  bvals_mhd_fun(pDomain, right_x2, shkset2d_ojb);

/* Compute Analytic solution for Sod and RJ4a tests, if required */

  tlim = par_getd("time","tlim");
  err_test = par_getd_def("problem","error_test",0);
  if (err_test == 1) {

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
          r0 = 1.0;
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
          r0 = 0.0;
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
  double err[8], tot_err[8];
  int ierr,myID;
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

  /* Sum up the Computed Error */
  ierr = MPI_Reduce(err,tot_err,8,MPI_DOUBLE,MPI_SUM,0,
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
  rms_error = sqrt(rms_error)/(double)count;

/* Print warning to stdout if rms_error exceeds estimate of 1st-order conv */

  min_zones = Nx1;
  if (Nx2 > 1) min_zones = MIN(min_zones,Nx2);
#ifdef HYDRO
  if (rms_error > 8.0/min_zones)
    printf("WARNING: rms_error=%e exceeds estimate\n",rms_error);
#endif
#ifdef MHD
  if (rms_error > 12.0/min_zones)
    printf("WARNING: rms_error=%e exceeds estimate\n",rms_error);
#endif

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

  fprintf(fp,"\n");

  fclose(fp);

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_iib(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_iib(GridS *pGrid)
{
  const int is = pGrid->is;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost + r2;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U  [k][j][is-i] = pGrid->U  [k][j-r2][is-i+r1];
#ifdef MHD
	pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j-r2][is-i+r1];
	pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j-r2][is-i+r1];
	pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j-r2][is-i+r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_oib(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_oib(GridS *pGrid)
{
  const int ie = pGrid->ie;
  int i, j, k, ju, jl, kl, ku; /* j-upper, j-lower */

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost - r2;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that i=ie+1 is not a boundary condition for the interface field B1i */

  for (k=kl; k<=ku; k++) {
    for (i=1; i<=nghost; i++) { /* Do NOT Change this loop ordering! */
      for (j=jl; j<=ju; j++) {
	pGrid->U[k][j][ie+i] = pGrid->U[k][j+r2][ie+i-r1];
#ifdef MHD
	if(i>1) pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j+r2][ie+i-r1];
	pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j+r2][ie+i-r1];
	pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j+r2][ie+i-r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_ijb(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_ijb(GridS *pGrid)
{
  const int js = pGrid->js;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost + r1;

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U  [k][js-j][i] = pGrid->U  [k][js-j+r2][i-r1];
#ifdef MHD
	pGrid->B1i[k][js-j][i] = pGrid->B1i[k][js-j+r2][i-r1];
	pGrid->B2i[k][js-j][i] = pGrid->B2i[k][js-j+r2][i-r1];
	pGrid->B3i[k][js-j][i] = pGrid->B3i[k][js-j+r2][i-r1];
#endif
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void shkset2d_ojb(GridS *pGrid)
 *  \brief Sets ghost zones using the nearest computational grid
 * cells implied by the size of the unit cell (r1xr2).
 */

void shkset2d_ojb(GridS *pGrid)
{
  const int je = pGrid->je;
  int i, j, k, iu, il, kl, ku; /* i-upper, i-lower */

  iu = pGrid->ie + nghost - r1;
  il = pGrid->is - nghost;

  if (pGrid->Nx[2] > 1){
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  } else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Note that j=je+1 is not a boundary condition for the interface field B2i */

  for (k=kl; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
	pGrid->U[k][je+j][i] = pGrid->U[k][je+j-r2][i+r1];
#ifdef MHD
	pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je+j-r2][i+r1];
	if(j>1) pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je+j-r2][i+r1];
	pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je+j-r2][i+r1];
#endif
      }
    }
  }
  return;
}
