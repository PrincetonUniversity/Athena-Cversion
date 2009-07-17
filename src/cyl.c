#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "cyl.h"
//#include "debug.h"

#ifdef CYLINDRICAL


/*============================================================================
 * CONVERSION FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*  FUNCTION Gas_to_Prim
 *
 *    THIS FUNCTION WAS DEPRECATED FOR SOME REASON, SO WE RESTORE IT HERE. 
 *    IT IS MOSTLY USED FOR DEBUGGING ANYWAY.
 */
void Gas_to_Prim(const Gas *pU, Prim *pW)
{
#if (NSCALARS > 0)
  int n;
#endif
  Real di = 1.0/pU->d;

  pW->d  = pU->d;
  pW->V1 = pU->M1*di;
  pW->V2 = pU->M2*di;
  pW->V3 = pU->M3*di;

#ifndef ISOTHERMAL
  pW->P =  pU->E - 0.5*(SQR(pU->M1) + SQR(pU->M2) + SQR(pU->M3))*di;
#ifdef MHD
  pW->P -= 0.5*(SQR(pU->B1c) + SQR(pU->B2c) + SQR(pU->B3c));
#endif
  pW->P *= Gamma_1;
// if (ppW->P < 0.0) printf("WARNING:  P WAS NEGATIVE\n");
  pW->P = MAX(pW->P,TINY_NUMBER);
#endif /* ISOTHERMAL */

#ifdef MHD
  pW->B1c = pU->B1c;
  pW->B2c = pU->B2c;
  pW->B3c = pU->B3c;
#endif

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) pW->r[n] = pU->s[n]*di;
#endif

  return;
}


/*----------------------------------------------------------------------------*/
/*  FUNCTION Prim_to_Gas
 *
 *    THIS FUNCTION WAS DEPRECATED FOR SOME REASON, SO WE RESTORE IT HERE. 
 *    IT IS MOSTLY USED FOR DEBUGGING ANYWAY.
 */
void Prim_to_Gas(Gas *pU, const Prim *pW)
{
#if (NSCALARS > 0)
  int n;
#endif

  pU->d  = pW->d;
  pU->M1 = pW->V1*pW->d;
  pU->M2 = pW->V2*pW->d;
  pU->M3 = pW->V3*pW->d;

#ifndef ISOTHERMAL
  pU->E = MAX(pW->P,TINY_NUMBER)/Gamma_1 + 0.5*pW->d*(SQR(pW->V1) + SQR(pW->V2) + SQR(pW->V3));
#ifdef MHD
  pU->E += 0.5*(SQR(pW->B1c) + SQR(pW->B2c) + SQR(pW->B3c));
#endif
#endif /* ISOTHERMAL */

#ifdef MHD
  pU->B1c = pW->B1c;
  pU->B2c = pW->B2c;
  pU->B3c = pW->B3c;
#endif

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) pU->s[n] = pW.r[n]*pW->d;
#endif

  return;
}


/*============================================================================
 * BOUNDARY CONDITION FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*  FUNCTION do_nothing_bc
 *
 *  DOES ABSOLUTELY NOTHING!  THUS, WHATEVER THE BOUNDARY ARE SET TO INITIALLY,
 *  THEY REMAIN FOR ALL TIME.
 */
void do_nothing_bc(Grid *pG)
{
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions (w/diode condition), Inner x1 boundary
 */

void strict_outflow_ix1(Grid *pGrid)
{
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3,L;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      cc_pos(pGrid,is,j,k,&x1,&x2,&x3);
      L = x1*pGrid->U[k][j][is].M1;
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
        // HOLD ANGULAR MOMENTUM CONSTANT
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        pGrid->U[k][j][is-i].M1 = L/x1;
        pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      }
    }
  }


  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions (w/diode condition), Outer x1 boundary
 */

void strict_outflow_ox1(Grid *pGrid)
{
  int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
  Real x1,x2,x3,L;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      cc_pos(pGrid,ie,j,k,&x1,&x2,&x3);
      L = x1*pGrid->U[k][j][ie].M1;
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
        // HOLD ANGULAR MOMENTUM CONSTANT
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        pGrid->U[k][j][ie+i].M1 = L/x1;
        pGrid->U[k][j][ie+i].M1 = MAX(pGrid->U[k][j][ie+i].M1,0.0);
      }
    }
  }

#ifdef MHD
/* Note that i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=2; i<=nghost; i++) {
        pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx2 > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      }
    }
  }

  if (pGrid->Nx3 > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
        pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions w/diode condition, Inner x3 boundary (ibc_x3=2)
 * Phi set by multipole expansion external to this function
 */
void diode_outflow_ix3(Grid *pGrid)
{
  int ks = pGrid->ks;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U  [ks-k][j][i] = pGrid->U  [ks][j][i];
        pGrid->U[ks-k][j][i].M3 = MIN(pGrid->U[ks-k][j][i].M3,0.0);
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ks-k][j][i] = pGrid->B1i[ks][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ks-k][j][i] = pGrid->B2i[ks][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ks-k][j][i] = pGrid->B3i[ks][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}

/*----------------------------------------------------------------------------*/
/* OUTFLOW boundary conditions w/diode condition, Outer x3 boundary (obc_x3=2)
 * Phi set by multipole expansion external to this function
 */
void diode_outflow_ox3(Grid *pGrid)
{
  int ke = pGrid->ke;
  int i,j,k,il,iu,jl,ju; /* i-lower/upper;  j-lower/upper */

  if (pGrid->Nx1 > 1){
    iu = pGrid->ie + nghost;
    il = pGrid->is - nghost;
  } else {
    iu = pGrid->ie;
    il = pGrid->is;
  }
  if (pGrid->Nx2 > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[ke+k][j][i] = pGrid->U[ke][j][i];
        pGrid->U[ke+k][j][i].M3 = MAX(pGrid->U[ke+k][j][i].M3,0.0);
      }
    }
  }

#ifdef MHD
  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[ke+k][j][i] = pGrid->B1i[ke][j][i];
      }
    }
  }

  for (k=1; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[ke+k][j][i] = pGrid->B2i[ke][j][i];
      }
    }
  }

/* Note that k=ke+1 is not a boundary condition for the interface field B3i */
  for (k=2; k<=nghost; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[ke+k][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

  return;
}


/*============================================================================
 * ERROR-ANALYSIS FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* FUNCTION compute_div_B
 *
 *  COMPUTE THE DIVERGENCE OF THE MAGNETIC FIELD USING FACE-CENTERED FIELDS
 *  OVER THE ENTIRE ACTIVE GRID.  RETURNS THE MAXIMUM OF DIV B IN ABSOLUTE VALUE.
 */
Real compute_div_b(Grid *pG)
{
#ifdef MHD
  int i,j,k,is,ie,js,je,ks,ke;
  Real x1,x2,x3,divB,maxdivB=0.0;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        divB = (pG->dx2)*(pG->dx3)*((x1+0.5*pG->dx1)*pG->B1i[k][j][i+1] - (x1-0.5*pG->dx1)*pG->B1i[k][j][i])
             + (pG->dx1)*(pG->dx3)*(pG->B2i[k][j+1][i] - pG->B2i[k][j][i]);
        if (ke > ks)
          divB += (pG->dx1)*(x1*pG->dx2)*(pG->B3i[k+1][j][i] - pG->B3i[k][j][i]);
        maxdivB = MAX(maxdivB,fabs(divB));
      }
    }
  }

  return maxdivB;

#else
  fprintf(stderr,"[compute_div_b]: This only works for MHD!\n");
  exit(EXIT_FAILURE);
  return 0.0;
#endif /* MHD */
}

/*----------------------------------------------------------------------------*/
/* FUNCTION compute_l1_error
 *
 *  COMPUTE THE L1-ERRORS IN ALL VARIABLES AT THE CURRENT (USUALLY THE FINAL)
 *  TIMESTEP USING THE INITIAL SOLUTION.  THIS MEANS THAT THE SOLUTION MUST 
 *  EITHER BE STATIC (STEADY-STATE) OR MUST HAVE COMPLETED A FULL PERIOD OF
 *  ROTATION OR WHATEVER.  FOR THE ERRORTYPE FLAG, 0 MEANS ABSOLUTE ERROR, AND
 *  1 MEANS RELATIVE ERROR.  DO NOT USE RELATIVE ERROR IF THE SOLUTION IS 
 *  INITIALLY ZERO!!!
 */
void compute_l1_error(char *problem, Grid *pG, Domain *pDomain, Gas ***Soln, const int errortype)
{
  int i=0,j=0,k=0;
#if (NSCALARS > 0)
   int n;
#endif
  int is,ie,js,je,ks,ke;
  Real rms_error=0.0;
  Real x1,x2,x3,x1min,x1max,x2min,x2max,x3min,x3max,grid_vol,dvol;
  Gas error,total_error;
  FILE *fp;
  char *fname, fnamestr[256];
  int Nx1, Nx2, Nx3, count;
#if defined MPI_PARALLEL
  double err[8+NSCALARS], tot_err[8+NSCALARS];
  int mpi_err;
#endif

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

/* compute L1 error in each variable, and rms total error */

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;

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
      cc_pos(pG,i,j,k,&x1,&x2,&x3);

      if (errortype == 0) {
        error.d   += x1*fabs(pG->U[k][j][i].d   - Soln[k][j][i].d );
        error.M1  += x1*fabs(pG->U[k][j][i].M1  - Soln[k][j][i].M1);
        error.M2  += x1*fabs(pG->U[k][j][i].M2  - Soln[k][j][i].M2);
        error.M3  += x1*fabs(pG->U[k][j][i].M3  - Soln[k][j][i].M3); 
#ifdef MHD
        error.B1c += x1*fabs(pG->U[k][j][i].B1c - Soln[k][j][i].B1c);
        error.B2c += x1*fabs(pG->U[k][j][i].B2c - Soln[k][j][i].B2c);
        error.B3c += x1*fabs(pG->U[k][j][i].B3c - Soln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
        error.E   += x1*fabs(pG->U[k][j][i].E   - Soln[k][j][i].E );
#endif /* ISOTHERMAL */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          error.s[n] += x1*fabs(pG->U[k][j][i].s[n] - Soln[k][j][i].s[n]);;
#endif

      }
      else {
        error.d   += Soln[k][j][i].d  != 0.0 ? x1*fabs(pG->U[k][j][i].d   - Soln[k][j][i].d )/fabs(Soln[k][j][i].d ) 
                    : x1*fabs(pG->U[k][j][i].d   - Soln[k][j][i].d );
        error.M1  += Soln[k][j][i].M1 != 0.0 ? x1*fabs(pG->U[k][j][i].M1  - Soln[k][j][i].M1)/fabs(Soln[k][j][i].M1) 
                    : x1*fabs(pG->U[k][j][i].M1  - Soln[k][j][i].M1 );
        error.M2  += Soln[k][j][i].M2 != 0.0 ? x1*fabs(pG->U[k][j][i].M2  - Soln[k][j][i].M2)/fabs(Soln[k][j][i].M2) 
                    : x1*fabs(pG->U[k][j][i].M2  - Soln[k][j][i].M2 );
        error.M3  += Soln[k][j][i].M3 != 0.0 ? x1*fabs(pG->U[k][j][i].M3  - Soln[k][j][i].M3)/fabs(Soln[k][j][i].M3) 
                    : x1*fabs(pG->U[k][j][i].M3  - Soln[k][j][i].M3 );
#ifndef ISOTHERMAL
        error.E   += Soln[k][j][i].E  != 0.0 ? x1*fabs(pG->U[k][j][i].E   - Soln[k][j][i].E )/fabs(Soln[k][j][i].E ) 
                    : x1*fabs(pG->U[k][j][i].E   - Soln[k][j][i].E );
#endif /* ISOTHERMAL */
      }
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

#if defined MPI_PARALLEL
  Nx1 = pDomain->ixe - pDomain->ixs + 1;
  Nx2 = pDomain->jxe - pDomain->jxs + 1;
  Nx3 = pDomain->kxe - pDomain->kxs + 1;
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
  mpi_err = MPI_Reduce(err, tot_err, (8+NSCALARS),
		       MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(mpi_err)
    ath_error("[compute_l1_error]: MPI_Reduce call returned error = %d\n",
	      mpi_err);

/* If I'm the parent, copy the sum back to the total_error variable */
  if(pG->my_id == 0){ /* I'm the parent */
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
  for (n=0; n<NSCALARS; n++) total_error.s[n] = err[8+n];
#endif

  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */


  // VOLUME AVERAGE
  x1min = pG->x1_0 + (is + pG->idisp    )*pG->dx1;
  x1max = pG->x1_0 + (ie + pG->idisp + 1)*pG->dx1;
  x2min = pG->x2_0 + (js + pG->jdisp    )*pG->dx2;
  x2max = pG->x2_0 + (je + pG->jdisp + 1)*pG->dx2;
  x3min = pG->x3_0 + (ks + pG->kdisp    )*pG->dx3;
  x3max = pG->x3_0 + (ke + pG->kdisp + 1)*pG->dx3;
  grid_vol = 0.5*(x1max + x1min)*(x1max - x1min)*(x2max - x2min);
  if (pG->Nx3 > 1)
    grid_vol *= (x3max - x3min);
  dvol = pG->dx1*pG->dx2/grid_vol;
  if (pG->Nx3 > 1)
    dvol *= pG->dx3;


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
//   rms_error = sqrt(rms_error)/(double)count;
  rms_error = sqrt(rms_error)*dvol;


/* Print error to file "BLAH-errors.#.dat"  */
   sprintf(fnamestr,"%s-errors",problem);
   fname = ath_fname("",fnamestr,1,0,NULL,"dat");

/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[compute_l1_error]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[compute_l1_error]: Unable to open file.\n");
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
/*	  (total_error.d/(double) count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count));*/
	  (total_error.d*dvol),
	  (total_error.M1*dvol),
	  (total_error.M2*dvol),
	  (total_error.M3*dvol));

#ifndef ISOTHERMAL
//   fprintf(fp,"  %e",(total_error.E/(double)count));
  fprintf(fp,"  %e",(total_error.E*dvol));
#endif /* ISOTHERMAL */

#ifdef MHD
  fprintf(fp,"  %e  %e  %e",
/*	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));*/
	  (total_error.B1c*dvol),
	  (total_error.B2c*dvol),
	  (total_error.B3c*dvol));
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) {
//       fprintf(fp,"  %e",total_error.s[n]/(double)count);
      fprintf(fp,"  %e",total_error.s[n]*dvol);
    }
#endif

  fprintf(fp,"\n");

  fclose(fp);
  free(fname);

  return;
}


/*============================================================================
 * ROOT-FINDING FUNCTIONS
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/* FUNCTION sign_change
 *
 *  SEARCH FOR A SIGN CHANGE.  THIS FUNCTION PARTITIONS THE INTERVAL (a0,b0) INTO
 *  2^k EQUALLY SPACED GRID POINTS, EVALUATES THE FUNCTION f AT THOSE POINTS,
 *  AND THEN SEARCHES FOR A SIGN CHANGE IN f BETWEEN ADJACENT GRID POINTS.  THE
 *  FIRST SUCH INTERVAL FOUND, (a,b), IS RETURNED.
 */
int sign_change(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *a, Real *b) {
  const int kmax=15;
  int k, n, i;
  Real delta, fk, fkp1;

  for (k=1; k<=kmax; k++) {
    n = pow(2,k);
    delta = (b0-a0)/(n-1);
    *a = a0;
    fk = func(x,*a);
    for (i=1; i<n; i++) {
      *b = *a + delta;
      fkp1 = func(x,*b);
      if (fkp1*fk < 0)
        return 1;
      *a = *b;
      fk = fkp1;
    }
  }
  ath_error("[sign_change]: No sign change was detected in (%f,%f) for x=%f!\n",a0,b0,x);
  return 0;
} 


/*----------------------------------------------------------------------------*/
/* FUNCTION bisection
 *
 *  THIS FUNCTION IMPLEMENTS THE BISECTION METHOD FOR ROOT FINDING.
 */
int bisection(Real (*func)(const Real,const Real), const Real a0, const Real b0, const Real x, Real *root) 
{
  const Real tol = 1.0E-10;
  const int maxiter = 400;
  Real a=a0, b=b0, c, fa, fb, fc;
  int i;

  fa = func(x,a);
  fb = func(x,b);
  if (fabs(fa) < tol) {
    *root = a;
    return 1;
  }
  if (fabs(fb) < tol) {
    *root = b;
    return 1;
  }
// printf("fa = %f, fb = %f\n", fa, fb);

  for (i = 0; i < maxiter; i++) {
    c = 0.5*(a+b);
// printf("x = %f, a = %f, b = %f, c = %f\n", x,a,b,c);
#ifdef MYDEBUG
    printf("c = %f\n", c);
#endif
    if (fabs((b-a)/c) < tol) {
#ifdef MYDEBUG
      printf("Bisection converged within tolerance of %f!\n", eps);
#endif
      *root = c;
      return 1;
    }
    fc = func(x,c);
    if (fa*fc < 0) {
      b = c;
      fb = fc;
    }
    else if (fc*fb < 0) {
      a = c;
      fa = fc;
    }
    else if (fc == 0) {
      *root = c;
      return 1;
    }
    else {
      printf("Error:  There is no single root in (%f,%f) for x = %13.10f!!\n", a, b,x);
      *root = c;
      return 0;
    }
  }

  printf("Error:  Bisection did not converge in %d iterations for x = %13.10f!!\n", maxiter,x);
  *root = c;
  return 0;
}


#endif /* CYLINDRICAL */

