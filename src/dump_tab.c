#include "copyright.h"
/*==============================================================================
 * FILE: dump_tab.c
 *
 * PURPOSE: Functions to write a dump as a formatted table.  Note that all the
 *   data is output for all Grids over the whole Mesh hierarchy on this
 *   processor using formatted writes, so the resulting output files can be
 *   extremely large.  Useful for 1D calculations, some 2D calculations, and
 *   for very small 3D runs.  
 *
 * REMINDER: use the slicing option available in output_tab() to write selected
 *   variables as a formatted table along any arbitrary 1D slice, or in any
 *   sub-volume, of 2D or 3D calculations.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   dump_tab_cons - writes conserved variables as formatted table
 *   dump_tab_prim - writes primitive variables as formatted table
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

void dump_tab_cons(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  int dnum = pOut->num;
  int nl,nd,i,j,k,il,iu,jl,ju,kl,ku,nGrid,nvar,nscal,ngrav;
  FILE *pfile;
  Real x1,x2,x3;
  char zone_fmt[20], fmt[80];
  int col_cnt, nmax;
#if (NSCALARS > 0)
  int n;
#endif

/* Open the output file */
  if((pfile = ath_fopen(NULL,pM->outfilename,num_digit,dnum,NULL,"tab","w")) 
     == NULL){
    ath_error("[dump_tab_cons]: File Open Error Occured");
    return;
  }

/* Write out some header information about Mesh hierarchy */

  fprintf(pfile,"# Dump of CONSERVED variables at Time = %g\n",pM->time);
  fprintf(pfile,"# NLevels in Mesh hierarchy\n%d\n",pM->NLevels);
  fprintf(pfile,"# DomainsPerLevel in Mesh hierarchy\n");
  for (nl=0;nl<pM->NLevels;nl++) fprintf(pfile,"  %d",pM->DomainsPerLevel[nl]);

/* Find and write number of Grids output in this file (on this processor) */

  nGrid=0;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) nGrid++;
    }
  }
  fprintf(pfile,"\n# nGrids in this file\n%d\n",nGrid);

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG = pM->Domain[nl][nd].Grid;
        col_cnt = 1;

/* Upper and Lower bounds on i,j,k for data dump */
        il = pG->is; iu = pG->ie;
        jl = pG->js; ju = pG->je;
        kl = pG->ks; ku = pG->ke;

        nmax =  pG->Nx[0] > pG->Nx[1]  ? pG->Nx[0] : pG->Nx[1];
        nmax = (pG->Nx[2] > nmax ? pG->Nx[2] : nmax);

#ifdef WRITE_GHOST_CELLS
        iu = pG->ie + nghost;
        il = pG->is - nghost;

        if(pG->Nx[1] > 1) {
          ju = pG->je + nghost;
          jl = pG->js - nghost;
        }

        if(pG->Nx[2] > 1) {
          ku = pG->ke + nghost;
          kl = pG->ks - nghost;
        }
        nmax += 2*nghost;
#endif

/* Add a white space to the format, setup format for integer zone columns */
        if(pOut->dat_fmt == NULL){
          sprintf(fmt," %%12.8e"); /* Use a default format */
        }
        else{
          sprintf(fmt," %s",pOut->dat_fmt);
        }
        sprintf(zone_fmt,"%%%dd", (int)(2+log10((double)(nmax))));

/* Write info about this Grid */

        ngrav = 0;
        nvar = NVAR;
        nscal = NSCALARS;
#ifdef SELF_GRAVITY
        ngrav = 1;
#endif
        fprintf(pfile,"\n#  Level Domain Nx1  Nx2  Nx3 nvar nscal ngrav\n");
        fprintf(pfile,"   %d     %d      %d  %d  %d     %d    %d    %d\n"
          ,nl,nd,(iu-il+1),(ju-jl+1),(ku-kl+1),nvar,nscal,ngrav);

/* write out i,j,k column headers.  Note column number is embedded in header */

        fprintf(pfile,"# [%d]=i-zone",col_cnt);
        col_cnt++;
        if (pG->Nx[1] > 2) {
          fprintf(pfile," [%d]=j-zone",col_cnt);
          col_cnt++;
        }
        if (pG->Nx[2] > 3) {
          fprintf(pfile," [%d]=k-zone",col_cnt);
          col_cnt++;
        }

/* write out x1,x2,x3 column headers.  */

        if (pG->Nx[0] > 1) {
          fprintf(pfile," [%d]=x1",col_cnt);
          col_cnt++;
        }
        if (pG->Nx[1] > 2) {
          fprintf(pfile," [%d]=x2",col_cnt);
          col_cnt++;
        }
        if (pG->Nx[2] > 3) {
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

/* write out E column header, if not barotropic */
#ifndef BAROTROPIC
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
        fprintf(pfile,"\n");

/* Write out data */

        for(k=kl; k<=ku; k++){
          for(j=jl; j<=ju; j++){
            for(i=il; i<=iu; i++){
              cc_pos(pG,i,j,k,&x1,&x2,&x3);

              if (pG->Nx[0] > 1) fprintf(pfile,zone_fmt,i);
              if (pG->Nx[1] > 1) fprintf(pfile,zone_fmt,j);
              if (pG->Nx[2] > 1) fprintf(pfile,zone_fmt,k);
              if (pG->Nx[0] > 1) fprintf(pfile,fmt,x1);
              if (pG->Nx[1] > 1) fprintf(pfile,fmt,x2);
              if (pG->Nx[2] > 1) fprintf(pfile,fmt,x3);

/* Dump all variables */

              fprintf(pfile,fmt,pG->U[k][j][i].d);
              fprintf(pfile,fmt,pG->U[k][j][i].M1);
              fprintf(pfile,fmt,pG->U[k][j][i].M2);
              fprintf(pfile,fmt,pG->U[k][j][i].M3);

#ifndef BAROTROPIC
              fprintf(pfile,fmt,pG->U[k][j][i].E);
#endif /* BAROTROPIC */

#ifdef MHD
              fprintf(pfile,fmt,pG->U[k][j][i].B1c);
              fprintf(pfile,fmt,pG->U[k][j][i].B2c);
              fprintf(pfile,fmt,pG->U[k][j][i].B3c);
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
              for (n=0; n<NSCALARS; n++) fprintf(pfile,fmt,pG->U[k][j][i].s[n]);
#endif

      	      fprintf(pfile,"\n");
            }
          }
        }
      }
    } /* end loop over domains */
  } /* end loop over levels */

  fclose(pfile);

  return;
}

/*----------------------------------------------------------------------------*/
/* dump_tab_prim: output PRIMITIVE variables  */

void dump_tab_prim(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  int dnum = pOut->num;
  int nl,nd,i,j,k,il,iu,jl,ju,kl,ku,nGrid,nvar,nscal,ngrav;
  FILE *pfile;
  PrimVarS W;
  Real x1,x2,x3,d1;
  char zone_fmt[20], fmt[80];
  int col_cnt, nmax;
#if (NSCALARS > 0)
  int n;
#endif

/* Open the output file */
  if((pfile = ath_fopen(NULL,pM->outfilename,num_digit,dnum,NULL,"tab","w")) 
     == NULL){
    ath_error("[dump_tab_prim]: File Open Error Occured");
    return;
  }

/* Write out some header information about Mesh hierarchy */

  fprintf(pfile,"# Dump of PRIMITIVE variables at Time = %g\n",pM->time);
  fprintf(pfile,"# NLevels in Mesh hierarchy\n%d\n",pM->NLevels);
  fprintf(pfile,"# DomainsPerLevel in Mesh hierarchy\n");
  for (nl=0; nl<pM->NLevels; nl++) fprintf(pfile,"  %d",pM->DomainsPerLevel[nl]);

/* Find and write number of Grids output in this file (on this processor) */

  nGrid=0;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) nGrid++;
    }
  }
  fprintf(pfile,"\n# nGrids in this file\n%d\n",nGrid);

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        pG = pM->Domain[nl][nd].Grid;
        col_cnt = 1;

/* Upper and Lower bounds on i,j,k for data dump */

        il = pG->is; iu = pG->ie;
        jl = pG->js; ju = pG->je;
        kl = pG->ks; ku = pG->ke;

        nmax =  pG->Nx[0] > pG->Nx[1]  ? pG->Nx[0] : pG->Nx[1];
        nmax = (pG->Nx[2] > nmax ? pG->Nx[2] : nmax);

#ifdef WRITE_GHOST_CELLS
        iu = pG->ie + nghost;
        il = pG->is - nghost;

        if(pG->Nx[1] > 1) {
          ju = pG->je + nghost;
          jl = pG->js - nghost;
        }

        if(pG->Nx[2] > 1) {
          ku = pG->ke + nghost;
          kl = pG->ks - nghost;
        }
        nmax += 2*nghost;
#endif

/* Add a white space to the format, setup format for integer zone columns */
        if(pOut->dat_fmt == NULL){
          sprintf(fmt," %%12.8e"); /* Use a default format */
        }
        else{
          sprintf(fmt," %s",pOut->dat_fmt);
        }
        sprintf(zone_fmt,"%%%dd", (int)(2+log10((double)(nmax))));

/* Write info about this Grid */

        ngrav = 0;
        nvar = NVAR;
        nscal = NSCALARS;
#ifdef SELF_GRAVITY
        ngrav = 1;
#endif
        fprintf(pfile,"\n#  Level Domain Nx1  Nx2  Nx3 nvar nscal ngrav\n");
        fprintf(pfile,"   %d     %d      %d  %d  %d     %d    %d    %d\n"
          ,nl,nd,(iu-il+1),(ju-jl+1),(ku-kl+1),nvar,nscal,ngrav);

/* write out i,j,k column headers.  Note column number is embedded in header */

        fprintf(pfile,"# [%d]=i-zone",col_cnt);
        col_cnt++;
        if (pG->Nx[1] > 2) {
          fprintf(pfile," [%d]=j-zone",col_cnt);
          col_cnt++;
        }
        if (pG->Nx[2] > 3) {
          fprintf(pfile," [%d]=k-zone",col_cnt);
          col_cnt++;
        }

/* write out x1,x2,x3 column headers.  */

        fprintf(pfile," [%d]=x1",col_cnt);
        col_cnt++;
        if (pG->Nx[1] > 2) {
          fprintf(pfile," [%d]=x2",col_cnt);
          col_cnt++;
        }
        if (pG->Nx[2] > 3) {
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

/* write out P column header, if not barotropic */
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
        fprintf(pfile,"\n");

/* Write out data */

        for(k=kl; k<=ku; k++){
          for(j=jl; j<=ju; j++){
            for(i=il; i<=iu; i++){
              cc_pos(pG,i,j,k,&x1,&x2,&x3);
              W = Cons_to_Prim(&(pG->U[k][j][i])); 

              if (pG->Nx[0] > 1) fprintf(pfile,zone_fmt,i);
              if (pG->Nx[1] > 1) fprintf(pfile,zone_fmt,j);
              if (pG->Nx[2] > 1) fprintf(pfile,zone_fmt,k);
              if (pG->Nx[0] > 1) fprintf(pfile,fmt,x1);
              if (pG->Nx[1] > 1) fprintf(pfile,fmt,x2);
              if (pG->Nx[2] > 1) fprintf(pfile,fmt,x3);

/* Dump all variables */

              fprintf(pfile,fmt,W.d);
              fprintf(pfile,fmt,W.V1);
              fprintf(pfile,fmt,W.V2);
              fprintf(pfile,fmt,W.V3);

#ifndef BAROTROPIC
              fprintf(pfile,fmt,W.P);
#endif /* BAROTROPIC */

#ifdef MHD
              fprintf(pfile,fmt,W.B1c);
              fprintf(pfile,fmt,W.B2c);
              fprintf(pfile,fmt,W.B3c);
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
              for (n=0; n<NSCALARS; n++) fprintf(pfile,fmt,W.r[n]);
#endif
              fprintf(pfile,"\n");
            }
          }
        }
      }
    } /* end loop over domains */
  } /* end loop over levels */

  fclose(pfile);

  return;
}
