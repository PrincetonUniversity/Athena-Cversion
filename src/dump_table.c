#include "copyright.h"
/*==============================================================================
 * FILE: dump_table.c
 *
 * PURPOSE: Function to write a dump as a formatted table.  Note that it writes
 *   all the data over the whole output grid using formatted writes, so the
 *   resulting output files can be extremely large.  Useful for 1D calculations,
 *   some 2D calculations, and for de-bugging small 3D runs.
 *
 * REMINDER: use the slicing option available in output_tab() to write selected
 *   variables as a formatted table along any arbitrary 1D slice, or in any
 *   sub-volume, of 2D or 3D calculations.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   dump_table - writes conserved variables + P as formatted table
 *============================================================================*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* dump_table  */

void dump_table(Grid *pG, Output *pOut)
{
  int dnum = pOut->num;
  FILE *pfile;
  Gas *pq;
  Real KE,ME=0.0;
  Real x1,x2,x3;
  char zone_fmt[20];
  char fmt[80];
/* Upper and Lower bounds on i,j,k for data dump */
  int i, il = pG->is; int iu = pG->ie;
  int j, jl = pG->js; int ju = pG->je;
  int k, kl = pG->ks; int ku = pG->ke;

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
#endif

/* Open the output file */
  if((pfile = ath_fopen(pG->outfilename,num_digit,dnum,NULL,"tab","w")) 
     == NULL){
    ath_error("[dump_table]: File Open Error Occured");
    return;
  }

/* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }
  sprintf(zone_fmt,"%%%dd %%%dd %%%dd",
	  (int)(1+log10((double)(iu-il+1))),
	  (int)(1+log10((double)(ju-jl+1))),
	  (int)(1+log10((double)(ku-kl+1))));

/* Write out some header information */
  fprintf(pfile,"# Nx1 = %d\n",iu-il+1);
  fprintf(pfile,"# x1-size = %g\n",(iu-il+1)*pG->dx1);
  fprintf(pfile,"# Nx2 = %d\n",ju-jl+1);
  fprintf(pfile,"# x2-size = %g\n",(ju-jl+1)*pG->dx2);
  fprintf(pfile,"# Nx3 = %d\n",ku-kl+1);
  fprintf(pfile,"# x3-size = %g\n",(ku-kl+1)*pG->dx3);
  fprintf(pfile,"# Time = %g\n",pG->time);
  fprintf(pfile,"#\n");
  fprintf(pfile,"# i-zone# j-zone# k-zone# x1 x2 x3 d M1 M2 M3");

#ifndef ISOTHERMAL
  fprintf(pfile," P E");
#endif /* ISOTHERMAL */

#ifdef MHD
  fprintf(pfile," B1c B2c B3c B1i B2i B3i");
#endif /* MHD */
  fprintf(pfile,"\n#\n");

/* Write out the table */
  for(k=kl; k<=ku; k++){
    for(j=jl; j<=ju; j++){
      for(i=il; i<=iu; i++){
	pq = &(pG->U[k][j][i]);

/* Calculate the cell center position of the cell i,j,k */
	cc_pos(pG,i,j,k,&x1,&x2,&x3);

	fprintf(pfile,zone_fmt,i,j,k);
	fprintf(pfile,fmt,x1);
	fprintf(pfile,fmt,x2);
	fprintf(pfile,fmt,x3);
	fprintf(pfile,fmt,pq->d);
	fprintf(pfile,fmt,pq->M1);
	fprintf(pfile,fmt,pq->M2);
	fprintf(pfile,fmt,pq->M3);

#ifndef ISOTHERMAL
	KE = 0.5*(pq->M1*pq->M1 + pq->M2*pq->M2 + pq->M3*pq->M3)/pq->d;
#ifdef MHD
	ME = 0.5*(pq->B1c*pq->B1c + pq->B2c*pq->B2c + pq->B3c*pq->B3c); 
#endif /* MHD */
	fprintf(pfile,fmt,(pq->E - KE - ME)*Gamma_1);
	fprintf(pfile,fmt,pq->E);
#endif /* ISOTHERMAL */

#ifdef MHD
	fprintf(pfile,fmt,pq->B1c);
	fprintf(pfile,fmt,pq->B2c);
	fprintf(pfile,fmt,pq->B3c);
	fprintf(pfile,fmt,pG->B1i[k][j][i]);
	fprintf(pfile,fmt,pG->B2i[k][j][i]);
	fprintf(pfile,fmt,pG->B3i[k][j][i]);
#endif
	fprintf(pfile,"\n");
      }
    }
  }
  fclose(pfile);

  return;
}
