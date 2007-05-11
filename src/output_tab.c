#include "copyright.h"
/*==============================================================================
 * FILE: output_tab.c
 *
 * PURPOSE: Functions for writing output in tabular format.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_tab() - opens file and calls appropriate 1D/2D/3D output function
 *     Uses subset1,2,3() to extract appropriate section to be output.
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   output_tab_1d() - write tab file for 1D slice of data
 *   output_tab_2d() - write tab file for 2D plane of data
 *   output_tab_3d() - write tab file for 3D section of data
 *============================================================================*/

void output_tab_1d(Grid *pGrid, Output *pOut, FILE *pFile);
void output_tab_2d(Grid *pGrid, Output *pOut, FILE *pFile);
void output_tab_3d(Grid *pGrid, Output *pOut, FILE *pFile);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_tab:  open file, call 1D/2D/3D writer; called by data_ouput  */

void output_tab(Grid *pGrid, Domain *pD, Output *pOut)
{
  FILE *pFile;
  char *fname;
  int dnum = pOut->num;

/* construct output filename */
  fname = fname_construct(pGrid->outfilename,num_digit,dnum,pOut->id,"tab");
  if (fname == NULL) {
    ath_error("[output_tab]: Error constructing output filename\n");
    return;
  }

/* open filename */
  pFile = fopen(fname,"w");
  if (pFile == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
    return;
  }

/* call appropriate 1D, 2D or 3D output routine */
  if (pOut->ndim == 3) {
    output_tab_3d(pGrid,pOut,pFile);
  } else if (pOut->ndim == 2) {
    output_tab_2d(pGrid,pOut,pFile);
  } else if (pOut->ndim == 1) {
    output_tab_1d(pGrid,pOut,pFile);
  }

  fclose(pFile);
  free(fname);
  return;
}

/*----------------------------------------------------------------------------*/
/* output_tab_1d: writes 1D data  */

void output_tab_1d(Grid *pGrid, Output *pOut, FILE *pFile)
{
  int i,nx1;
  float *data, dmin, dmax, xworld;

  data = subset1(pGrid,pOut);
  nx1 = pOut->Nx1;     /* we know it's 1dim data */
  minmax1(data,nx1,&dmin,&dmax);

  for (i=0; i<pOut->Nx1; i++) {
    xworld = pOut->x1_0  + pOut->dx1*(float)(i);
    fprintf(pFile,"%g %g\n",xworld,data[i]);
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free_1d_array((void *)data); /* Free the memory we malloc'd */
}

/*----------------------------------------------------------------------------*/
/* output_tab_2d: writes 2D data  */

void output_tab_2d(Grid *pGrid, Output *pOut, FILE *pFile)
{
  int i,j,nx1,nx2;
  float **data, dmin, dmax, xworld, yworld;

  data = subset2(pGrid,pOut);
  nx1 = pOut->Nx1;     /* we know it's 2dim data */
  nx2 = pOut->Nx2;
  minmax2(data,nx2,nx1,&dmin,&dmax);

  for (j=0; j<pOut->Nx2; j++) {
    for (i=0; i<pOut->Nx1; i++) {
      xworld = pOut->x1_0  + pOut->dx1*(float)(i);
      yworld = pOut->x2_0  + pOut->dx2*(float)(j);
      fprintf(pFile,"%g %g %g\n",xworld,yworld,data[j][i]);
    }
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free_2d_array(data); /* Free the memory we malloc'd */
}

/*----------------------------------------------------------------------------*/
/* output_tab_3d: writes 3D data   */

void output_tab_3d(Grid *pGrid, Output *pOut, FILE *pFile)
{
  int i,j,k,nx1,nx2,nx3;
  float ***data, dmin, dmax, xworld, yworld, zworld;

  data = subset3(pGrid,pOut);
  nx1 = pOut->Nx1;     /* we know it's 3dim data */
  nx2 = pOut->Nx2;
  nx3 = pOut->Nx3;
  minmax3(data,nx3,nx2,nx1,&dmin,&dmax);

  for (k=0; k<pOut->Nx3; k++) {
    for (j=0; j<pOut->Nx2; j++) {
      for (i=0; i<pOut->Nx1; i++) {
        xworld = pOut->x1_0  + pOut->dx1*(float)(i);
        yworld = pOut->x2_0  + pOut->dx2*(float)(j);
        zworld = pOut->x3_0  + pOut->dx3*(float)(k);
        fprintf(pFile,"%g %g %g %g\n",xworld,yworld,zworld,data[k][j][i]);
      }
    }
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free_3d_array(data); /* Free the memory we malloc'd */
}
