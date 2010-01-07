#include "copyright.h"
/*==============================================================================
 * FILE: dump_vtk.c
 *
 * PURPOSE: Function to write a dump in VTK "legacy" format.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   dump_vtk() - writes VTK dump (all variables).
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#ifdef PARTICLES
#include "particles/particle.h"
#endif

/*----------------------------------------------------------------------------*/
/* dump_vtk:   */

void dump_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  FILE *pfile;
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku;
  int big_end = ath_big_endian();
  int ndata0;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3;
#if (NSCALARS > 0)
  int n;
#endif

/* Return if Grid is not on this processor */

  pGrid = pM->Domain[pOut->nlevel][pOut->ndomain].Grid;
  if (pGrid == NULL) return;

  il = pGrid->is, iu = pGrid->ie;
  jl = pGrid->js, ju = pGrid->je;
  kl = pGrid->ks, ku = pGrid->ke;

#ifdef WRITE_GHOST_CELLS
  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if(pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }

  if(pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
#endif /* WRITE_GHOST_CELLS */

/* Open output file, constructing filename in-line */

  pfile = ath_fopen(NULL,pM->outfilename,num_digit,pOut->num,NULL,"vtk","w");
  if(pfile == NULL){
    ath_error("[dump_vtk]: File Open Error Occured");
    return;
  }

/* Allocate memory for temporary array of floats */

  ndata0 = iu-il+1;
  if((data = (float *)malloc(3*ndata0*sizeof(float))) == NULL){
    ath_error("[dump_vtk]: malloc failed for temporary array\n");
    return;
  }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

  fprintf(pfile,"Really cool Athena data at time = %e\n",pGrid->time);

/*  3. File format */

  fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

  x1 = pGrid->MinX[0];
  x2 = pGrid->MinX[1];
  x3 = pGrid->MinX[2];

  fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,ku-kl+2);
  fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
  fprintf(pfile,"SPACING %e %e %e \n",pGrid->dx1,pGrid->dx2,pGrid->dx3);

/*  5. Data  */

  fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ju-jl+1)*(ku-kl+1));

/* Write density */

  fprintf(pfile,"SCALARS density float\n");
  fprintf(pfile,"LOOKUP_TABLE default\n");
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        data[i-il] = (float)pGrid->U[k][j][i].d;
      }
      if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
      fwrite(data,sizeof(float),(size_t)ndata0,pfile);
    }
  }

/* Write Velocity */

  fprintf(pfile,"\nVECTORS velocity float\n");
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
	data[3*(i-il)] = (float)(pGrid->U[k][j][i].M1/pGrid->U[k][j][i].d);
	data[3*(i-il)+1] = (float)(pGrid->U[k][j][i].M2/pGrid->U[k][j][i].d);
	data[3*(i-il)+2] = (float)(pGrid->U[k][j][i].M3/pGrid->U[k][j][i].d);
      }
      if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
      fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
    }
  }

/* Write total energy */

#ifndef BAROTROPIC
  fprintf(pfile,"\nSCALARS total_energy float\n");
  fprintf(pfile,"LOOKUP_TABLE default\n");
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        data[i-il] = (float)pGrid->U[k][j][i].E;
      }
      if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
      fwrite(data,sizeof(float),(size_t)ndata0,pfile);
    }
  }
#endif

/* Write cell centered B */

#ifdef MHD
  fprintf(pfile,"\nVECTORS cell_centered_B float\n");
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        data[3*(i-il)] = (float)pGrid->U[k][j][i].B1c;
        data[3*(i-il)+1] = (float)pGrid->U[k][j][i].B2c;
        data[3*(i-il)+2] = (float)pGrid->U[k][j][i].B3c;
      }
      if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
      fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
    }
  }
#endif

/* Write gravitational potential */

#ifdef SELF_GRAVITY
  fprintf(pfile,"\nSCALARS gravitational_potential float\n");
  fprintf(pfile,"LOOKUP_TABLE default\n");
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        data[i-il] = (float)pGrid->Phi[k][j][i];
      }
      if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
      fwrite(data,sizeof(float),(size_t)ndata0,pfile);
    }
  }
#endif

/* Write binned particle grid */

#ifdef PARTICLES
  if (pOut->out_pargrid) {
    fprintf(pfile,"\nSCALARS particle_density float\n");
    fprintf(pfile,"LOOKUP_TABLE default\n");
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          data[i-il] = pGrid->Coup[k][j][i].grid_d;
        }
        if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
        fwrite(data,sizeof(float),(size_t)ndata0,pfile);
      }
    }
    fprintf(pfile,"\nVECTORS particle_momentum float\n");
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
	  data[3*(i-il)] = pGrid->Coup[k][j][i].grid_v1;
	  data[3*(i-il)+1] = pGrid->Coup[k][j][i].grid_v2;
	  data[3*(i-il)+2] = pGrid->Coup[k][j][i].grid_v3;
        }
        if(!big_end) ath_bswap(data,sizeof(float),3*(iu-il+1));
        fwrite(data,sizeof(float),(size_t)(3*ndata0),pfile);
      }
    }
  }
#endif

/* Write passive scalars */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++){
    fprintf(pfile,"\nSCALARS scalar[%d] float\n",n);
    fprintf(pfile,"LOOKUP_TABLE default\n");
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {
          data[i-il] = (float)pGrid->U[k][j][i].s[n];
        }
        if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
        fwrite(data,sizeof(float),(size_t)ndata0,pfile);
      }
    }
  }
#endif

/* close file and free memory */

  fclose(pfile);
  free(data);
  return;
}
