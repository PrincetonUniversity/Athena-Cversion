#include "copyright.h"
/*==============================================================================
 * FILE: output_vtk.c
 *
 * PURPOSE: Function to write a single variable in VTK "legacy" format.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_vtk() - writes VTK file (single variable).
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   output_vtk_2d() - write vtk file for 2D data
 *   output_vtk_3d() - write vtk file for 3D data
 *============================================================================*/

static void output_vtk_2d(MeshS *pM, OutputS *pOut);
static void output_vtk_3d(MeshS *pM, OutputS *pOut);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_vtk:   */

void output_vtk(MeshS *pM, OutputS *pOut)
{
  if (pM->Domain[pOut->nlevel][pOut->ndomain].Grid == NULL) return;

  if (pOut->ndim == 3) {
    output_vtk_3d(pM, pOut);
  } else if (pOut->ndim == 2) {
    output_vtk_2d(pM, pOut);
  } else {
    ath_error("[output_vtk]: Only able to output 2D or 3D");
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* output_vtk_2d: writes 2D data  */

static void output_vtk_2d(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid=pM->Domain[pOut->nlevel][pOut->ndomain].Grid;
  FILE *pfile;
/* Upper and Lower bounds on i,j,k for data dump */
  int big_end = ath_big_endian();
  int ndata0;
  float **data2d=NULL;
  double x1, x2, x3, dx1, dx2, dx3;

/* Open output file, constructing filename in-line */

  pfile=ath_fopen(NULL,pM->outfilename,num_digit,pOut->num,pOut->id,"vtk","w");
  if(pfile == NULL){
    ath_error("[output_vtk]: File Open Error Occured");
    return;
  }

/* Allocate memory for temporary array of floats */
  data2d = subset2(pGrid,pOut);

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

  fprintf(pfile,"Really cool Athena data at time = %e\n",pGrid->time);

/*  3. File format */

  fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
  x1 = pGrid->x1min;
  x2 = pGrid->x2min;
  x3 = pGrid->x3min;
  ndata0 = pOut->Nx1 * pOut->Nx2 * pOut->Nx3;
  dx1 = (pOut->Nx1 == 1 ? pGrid->dx1 * pOut->Nx1 : pGrid->dx1);
  dx2 = (pOut->Nx2 == 1 ? pGrid->dx2 * pOut->Nx2 : pGrid->dx2);
  dx3 = (pOut->Nx3 == 1 ? pGrid->dx3 * pOut->Nx3 : pGrid->dx3);

  fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",pOut->Nx1+1,pOut->Nx2+1,pOut->Nx3+1);
  fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
  fprintf(pfile,"SPACING %e %e %e \n",dx1,dx2,dx3);

/*  5. Data  */

  fprintf(pfile,"CELL_DATA %d \n", ndata0);

/* Write density */

  fprintf(pfile,"SCALARS %s float\n", pOut->id);
  fprintf(pfile,"LOOKUP_TABLE default\n");
  if(!big_end) ath_bswap(data2d[0],sizeof(float),ndata0);
  fwrite(data2d[0],sizeof(float),(size_t)ndata0,pfile);

/* close file and free memory */

  fclose(pfile);
  free_2d_array(data2d);
  return;
}

/*----------------------------------------------------------------------------*/
/* output_vtk_3d: writes 3D data  */

static void output_vtk_3d(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid=pM->Domain[pOut->nlevel][pOut->ndomain].Grid;
  FILE *pfile;
/* Upper and Lower bounds on i,j,k for data dump */
  int big_end = ath_big_endian();
  int ndata0, k;
  float ***data3d=NULL;
  double x1, x2, x3;

/* Open output file, constructing filename in-line */

  pfile=ath_fopen(NULL,pM->outfilename,num_digit,pOut->num,pOut->id,"vtk","w");
  if(pfile == NULL){
    ath_error("[output_vtk]: File Open Error Occured");
    return;
  }

/* Allocate memory for temporary array of floats */
  data3d = subset3(pGrid,pOut);

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

  fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */

  fprintf(pfile,"Really cool Athena data at time = %e\n",pGrid->time);

/*  3. File format */

  fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
  x1 = pGrid->x1min;
  x2 = pGrid->x2min;
  x3 = pGrid->x3min;
  ndata0 = pOut->Nx1 * pOut->Nx2;

  fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
  fprintf(pfile,"DIMENSIONS %d %d %d\n",pOut->Nx1+1,pOut->Nx2+1,pOut->Nx3+1);
  fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
  fprintf(pfile,"SPACING %e %e %e \n",pGrid->dx1,pGrid->dx2,pGrid->dx3);

/*  5. Data  */

  fprintf(pfile,"CELL_DATA %d \n", pOut->Nx1 * pOut->Nx2 * pOut->Nx3);

/* Write density */

  fprintf(pfile,"SCALARS %s float\n", pOut->id);
  fprintf(pfile,"LOOKUP_TABLE default\n");
  for (k=0; k<pOut->Nx3; k++) {
    if(!big_end) ath_bswap(data3d[k][0],sizeof(float),ndata0);
    fwrite(data3d[k][0],sizeof(float),(size_t)ndata0,pfile);
  }

/* close file and free memory */

  fclose(pfile);
  free_3d_array(data3d);
  return;
}
