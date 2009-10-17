#include "copyright.h"
/*==============================================================================
 * FILE: dump_tab.c
 *
 * PURPOSE: Functions to write a dump as a formatted table.  Note that all the
 *   data is output over the whole Grid using formatted writes, so the
 *   resulting output files can be extremely large.  Useful for 1D calculations,
 *   some 2D calculations, and for very small 3D runs.
 *
 * REMINDER: use the slicing option available in output_tab() to write selected
 *   variables as a formatted table along any arbitrary 1D slice, or in any
 *   sub-volume, of 2D or 3D calculations.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   dump_tab_cons - writes conserved variables + P as formatted table
 *   dump_tab_prim - writes primitive variables     as formatted table
 *============================================================================*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#ifdef PARTICLES
#include "particles/particle.h"
#endif

/*----------------------------------------------------------------------------*/
/* dump_tab_cons: output CONSERVED variables  */

void dump_tab_cons(Grid *pG, Domain *pD, Output *pOut)
{
  int dnum = pOut->num;
  FILE *pfile;
  Gas *pGas;
  Real x1,x2,x3;
#ifndef BAROTROPIC
  Real KE,ME=0.0;
#endif
  char zone_fmt[20], fmt[80];
  int col_cnt=1, nmax;
#if (NSCALARS > 0)
  int n;
#endif
/* Upper and Lower bounds on i,j,k for data dump */
  int i, il = pG->is; int iu = pG->ie;
  int j, jl = pG->js; int ju = pG->je;
  int k, kl = pG->ks; int ku = pG->ke;

  nmax =  pG->Nx1 > pG->Nx2  ? pG->Nx1 : pG->Nx2;
  nmax = (pG->Nx3 > nmax ? pG->Nx3 : nmax);

#ifdef WRITE_GHOST_CELLS
  if(pG->Nx1 > 1) {
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  }

  if(pG->Nx2 > 1) {
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  }

  if(pG->Nx3 > 1) {
    ku = pG->ke + nghost;
    kl = pG->ks - nghost;
  }
  nmax += 2*nghost;
#endif

/* Open the output file */
  if((pfile = ath_fopen(NULL,pG->outfilename,num_digit,dnum,NULL,"tab","w")) 
     == NULL){
    ath_error("[dump_tab_cons]: File Open Error Occured");
    return;
  }

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }
  sprintf(zone_fmt,"%%%dd", (int)(1+log10((double)(nmax))));

/* Write out some header information */
  if (pG->Nx1 > 1) {
    fprintf(pfile,"# Nx1 = %d\n",iu-il+1);
    fprintf(pfile,"# x1-size = %g\n",(iu-il+1)*pG->dx1);
  }
  if (pG->Nx2 > 1) {
    fprintf(pfile,"# Nx2 = %d\n",ju-jl+1);
    fprintf(pfile,"# x2-size = %g\n",(ju-jl+1)*pG->dx2);
  }
  if (pG->Nx3 > 1) {
    fprintf(pfile,"# Nx3 = %d\n",ku-kl+1);
    fprintf(pfile,"# x3-size = %g\n",(ku-kl+1)*pG->dx3);
  }
  fprintf(pfile,"# Dump of CONSERVED variables at Time = %g\n",pG->time);
  fprintf(pfile,"#\n#");

/* write out i,j,k column headers.  Note column number is embedded in header */
  if (pG->Nx1 > 1) {
    fprintf(pfile," [%d]=i-zone",col_cnt);
    col_cnt++;
  }
  if (pG->Nx2 > 2) {
    fprintf(pfile," [%d]=j-zone",col_cnt);
    col_cnt++;
  }
  if (pG->Nx3 > 3) {
    fprintf(pfile," [%d]=k-zone",col_cnt);
    col_cnt++;
  }

/* write out x1,x2,x3 column headers.  */
  if (pG->Nx1 > 1) {
    fprintf(pfile," [%d]=x1",col_cnt);
    col_cnt++;
  }
  if (pG->Nx2 > 2) {
    fprintf(pfile," [%d]=x2",col_cnt);
    col_cnt++;
  }
  if (pG->Nx3 > 3) {
    fprintf(pfile," [%d]=x3",col_cnt);
    col_cnt++;
  }

/* write out d,M1,M2,M3 column headers */
  fprintf(pfile," [%d]=d",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=M1",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=M2",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=M3",col_cnt);
  col_cnt++;

/* write out P,E column headers, if not barotropic */
#ifndef BAROTROPIC
  fprintf(pfile," [%d]=P",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=E",col_cnt);
  col_cnt++;
#endif /* BAROTROPIC */

/* write out magnetic field component column headers, if mhd */
#ifdef MHD
  fprintf(pfile," [%d]=B1c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B2c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B3c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B1i",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B2i",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B3i",col_cnt);
  col_cnt++;
#endif /* MHD */

/* write out column header for gravitational potential (self-gravity) */
#ifdef SELF_GRAVITY
  fprintf(pfile," [%d]=Phi",col_cnt);
  col_cnt++;
#endif

/* write out column headers for particles */
#ifdef PARTICLES
  if (pOut->out_pargrid) {
    fprintf(pfile," [%d]=dpar",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=M1par",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=M2par",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=M3par",col_cnt);
    col_cnt++;
  }
#endif

/* write out column headers for passive scalars */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    fprintf(pfile," [%d]=s%d",col_cnt,n);
    col_cnt++;
  }
#endif

  fprintf(pfile,"\n#\n");

/* Write out data */
  for(k=kl; k<=ku; k++){
    for(j=jl; j<=ju; j++){
      for(i=il; i<=iu; i++){
	pGas = &(pG->U[k][j][i]);

/* Calculate the cell center position of the cell i,j,k */
	cc_pos(pG,i,j,k,&x1,&x2,&x3);

	if (pG->Nx1 > 1) fprintf(pfile,zone_fmt,i);
	if (pG->Nx2 > 1) fprintf(pfile,zone_fmt,j);
	if (pG->Nx3 > 1) fprintf(pfile,zone_fmt,k);
	if (pG->Nx1 > 1) fprintf(pfile,fmt,x1);
	if (pG->Nx2 > 1) fprintf(pfile,fmt,x2);
	if (pG->Nx3 > 1) fprintf(pfile,fmt,x3);

/* Dump all variables */
	fprintf(pfile,fmt,pGas->d);
	fprintf(pfile,fmt,pGas->M1);
	fprintf(pfile,fmt,pGas->M2);
	fprintf(pfile,fmt,pGas->M3);

#ifndef BAROTROPIC
	KE = 0.5*(pGas->M1*pGas->M1 + pGas->M2*pGas->M2 + pGas->M3*pGas->M3)/pGas->d;
#ifdef MHD
	ME = 0.5*(pGas->B1c*pGas->B1c + pGas->B2c*pGas->B2c + pGas->B3c*pGas->B3c); 
#endif /* MHD */
	fprintf(pfile,fmt,(pGas->E - KE - ME)*Gamma_1);
	fprintf(pfile,fmt,pGas->E);
#endif /* BAROTROPIC */

#ifdef MHD
	fprintf(pfile,fmt,pGas->B1c);
	fprintf(pfile,fmt,pGas->B2c);
	fprintf(pfile,fmt,pGas->B3c);
	fprintf(pfile,fmt,pG->B1i[k][j][i]);
	fprintf(pfile,fmt,pG->B2i[k][j][i]);
	fprintf(pfile,fmt,pG->B3i[k][j][i]);
#endif

#ifdef SELF_GRAVITY
        fprintf(pfile,fmt,pG->Phi[k][j][i]);
#endif

#ifdef PARTICLES
        if (pOut->out_pargrid) {
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_d);
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v1);
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v2);
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v3);
        }
#endif

#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) fprintf(pfile,fmt,pGas->s[n]);
#endif

	fprintf(pfile,"\n");
      }
    }
  }
  fclose(pfile);

  return;
}

/*----------------------------------------------------------------------------*/
/* dump_tab_prim: output PRIMITIVE variables  */

void dump_tab_prim(Grid *pG, Domain *pD, Output *pOut)
{
  int dnum = pOut->num;
  FILE *pfile;
  Prim Primitives;
  Real x1,x2,x3,d1;
#ifndef BAROTROPIC
  Real KE,ME=0.0;
#endif
  char zone_fmt[20], fmt[80];
  int col_cnt=1, nmax;
#if (NSCALARS > 0)
  int n;
#endif
/* Upper and Lower bounds on i,j,k for data dump */
  int i, il = pG->is; int iu = pG->ie;
  int j, jl = pG->js; int ju = pG->je;
  int k, kl = pG->ks; int ku = pG->ke;

  nmax =  pG->Nx1 > pG->Nx2  ? pG->Nx1 : pG->Nx2;
  nmax = (pG->Nx3 > nmax ? pG->Nx3 : nmax);

#ifdef WRITE_GHOST_CELLS
  if(pG->Nx1 > 1) {
    iu = pG->ie + nghost;
    il = pG->is - nghost;
  }

  if(pG->Nx2 > 1) {
    ju = pG->je + nghost;
    jl = pG->js - nghost;
  }

  if(pG->Nx3 > 1) {
    ku = pG->ke + nghost;
    kl = pG->ks - nghost;
  }
  nmax += 2*nghost;
#endif

/* Open the output file */
  if((pfile = ath_fopen(NULL,pG->outfilename,num_digit,dnum,NULL,"tab","w")) 
     == NULL){
    ath_error("[dump_tab_prim]: File Open Error Occured");
    return;
  }

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }
  sprintf(zone_fmt,"%%%dd", (int)(1+log10((double)(nmax))));

/* Write out some header information */
  if (pG->Nx1 > 1) {
    fprintf(pfile,"# Nx1 = %d\n",iu-il+1);
    fprintf(pfile,"# x1-size = %g\n",(iu-il+1)*pG->dx1);
  }
  if (pG->Nx2 > 1) {
    fprintf(pfile,"# Nx2 = %d\n",ju-jl+1);
    fprintf(pfile,"# x2-size = %g\n",(ju-jl+1)*pG->dx2);
  }
  if (pG->Nx3 > 1) {
    fprintf(pfile,"# Nx3 = %d\n",ku-kl+1);
    fprintf(pfile,"# x3-size = %g\n",(ku-kl+1)*pG->dx3);
  }
  fprintf(pfile,"# Dump of PRIMITIVE variables at Time = %g\n",pG->time);
  fprintf(pfile,"#\n#");

/* write out i,j,k column headers.  Note column number is embedded in header */
  if (pG->Nx1 > 1) {
    fprintf(pfile," [%d]=i-zone",col_cnt);
    col_cnt++;
  }
  if (pG->Nx2 > 2) {
    fprintf(pfile," [%d]=j-zone",col_cnt);
    col_cnt++;
  }
  if (pG->Nx3 > 3) {
    fprintf(pfile," [%d]=k-zone",col_cnt);
    col_cnt++;
  }

/* write out x1,x2,x3 column headers.  */
  if (pG->Nx1 > 1) {
    fprintf(pfile," [%d]=x1",col_cnt);
    col_cnt++;
  }
  if (pG->Nx2 > 2) {
    fprintf(pfile," [%d]=x2",col_cnt);
    col_cnt++;
  }
  if (pG->Nx3 > 3) {
    fprintf(pfile," [%d]=x3",col_cnt);
    col_cnt++;
  }

/* write out d,V1,V2,V3 column headers */
  fprintf(pfile," [%d]=d",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=V1",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=V2",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=V3",col_cnt);
  col_cnt++;

/* write out P,E column headers, if not barotropic */
#ifndef BAROTROPIC
  fprintf(pfile," [%d]=P",col_cnt);
  col_cnt++;
#endif /* BAROTROPIC */

/* write out magnetic field component column headers, if mhd */
#ifdef MHD
  fprintf(pfile," [%d]=B1c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B2c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B3c",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B1i",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B2i",col_cnt);
  col_cnt++;
  fprintf(pfile," [%d]=B3i",col_cnt);
  col_cnt++;
#endif /* MHD */

/* write out column header for gravitational potential (self-gravity) */
#ifdef SELF_GRAVITY
  fprintf(pfile," [%d]=Phi",col_cnt);
  col_cnt++;
#endif

/* write out column headers for particles */
#ifdef PARTICLES
  if (pOut->out_pargrid) {
    fprintf(pfile," [%d]=dpar",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=V1par",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=V2par",col_cnt);
    col_cnt++;
    fprintf(pfile," [%d]=V3par",col_cnt);
    col_cnt++;
  }
#endif

/* write out column headers for passive scalars */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    fprintf(pfile," [%d]=s%d",col_cnt,n);
    col_cnt++;
  }
#endif

  fprintf(pfile,"\n#\n");

/* Write out data */
  for(k=kl; k<=ku; k++){
    for(j=jl; j<=ju; j++){
      for(i=il; i<=iu; i++){
#ifdef SPECIAL_RELATIVITY
	Primitives = pG->W[k][j][i];
#else
        Primitives.d  = pG->U[k][j][i].d;
        d1 = 1.0/pG->U[k][j][i].d;
        Primitives.V1 = pG->U[k][j][i].M1*d1;
        Primitives.V2 = pG->U[k][j][i].M2*d1;
        Primitives.V3 = pG->U[k][j][i].M3*d1;
#ifndef BAROTROPIC
	KE = 0.5*(pG->U[k][j][i].M1*pG->U[k][j][i].M1 
                + pG->U[k][j][i].M2*pG->U[k][j][i].M2
                + pG->U[k][j][i].M3*pG->U[k][j][i].M3)*d1;
#ifdef MHD
	ME = 0.5*(pG->U[k][j][i].B1c*pG->U[k][j][i].B1c +
                  pG->U[k][j][i].B2c*pG->U[k][j][i].B2c +
                  pG->U[k][j][i].B3c*pG->U[k][j][i].B3c); 
#endif /* MHD */
	Primitives.P = (pG->U[k][j][i].E - KE - ME)*Gamma_1;
#endif /* BAROTROPIC */
#ifdef MHD
        Primitives.B1c = pG->U[k][j][i].B1c;
        Primitives.B2c = pG->U[k][j][i].B2c;
        Primitives.B3c = pG->U[k][j][i].B3c;
#endif /* MHD */
#endif /* SPECIAL_RELATIVITY */

/* Calculate the cell center position of the cell i,j,k */
	cc_pos(pG,i,j,k,&x1,&x2,&x3);

	if (pG->Nx1 > 1) fprintf(pfile,zone_fmt,i);
	if (pG->Nx2 > 1) fprintf(pfile,zone_fmt,j);
	if (pG->Nx3 > 1) fprintf(pfile,zone_fmt,k);
	if (pG->Nx1 > 1) fprintf(pfile,fmt,x1);
	if (pG->Nx2 > 1) fprintf(pfile,fmt,x2);
	if (pG->Nx3 > 1) fprintf(pfile,fmt,x3);

/* Dump all variables */
	fprintf(pfile,fmt,Primitives.d);
	fprintf(pfile,fmt,Primitives.V1);
	fprintf(pfile,fmt,Primitives.V2);
	fprintf(pfile,fmt,Primitives.V3);

#ifndef BAROTROPIC
	fprintf(pfile,fmt,Primitives.P);
#endif /* BAROTROPIC */

#ifdef MHD
	fprintf(pfile,fmt,Primitives.B1c);
	fprintf(pfile,fmt,Primitives.B2c);
	fprintf(pfile,fmt,Primitives.B3c);
	fprintf(pfile,fmt,pG->B1i[k][j][i]);
	fprintf(pfile,fmt,pG->B2i[k][j][i]);
	fprintf(pfile,fmt,pG->B3i[k][j][i]);
#endif

#ifdef SELF_GRAVITY
        fprintf(pfile,fmt,pG->Phi[k][j][i]);
#endif

#ifdef PARTICLES
        if (pOut->out_pargrid) {
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_d);
          if (pG->Coup[k][j][i].grid_d>0.0)
            d1 = 1.0/pG->Coup[k][j][i].grid_d;
          else
            d1 = 0.0;
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v1*d1);
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v2*d1);
          fprintf(pfile,fmt,pG->Coup[k][j][i].grid_v3*d1);
        }
#endif

#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) fprintf(pfile,fmt,Primitives.r[n]);
#endif

	fprintf(pfile,"\n");
      }
    }
  }
  fclose(pfile);

  return;
}
