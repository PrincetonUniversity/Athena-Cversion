#include "copyright.h"
/*==============================================================================
 * FILE: output_ppm.c
 *
 * PURPOSE: Writes single variable as a PPM image with color table.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_ppm()
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS:
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static float **data=NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   compute_rgb()  
 *============================================================================*/

static void compute_rgb(double data, double min, double max, int *pR, int *pG,
                    int *pB, Output *p);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_ppm:  output PPM image */

void output_ppm(Grid *pGrid, Domain *pD, Output *pOut)
{
  FILE *pfile;
  char *fname;
  int nx1,nx2,i,j;
  float dmin, dmax;
  int dnum = pOut->num;
  int red,green,blue;

/* check output data is 2D (output must be a 2D slice for 3D runs) */
  if (pOut->ndim != 2) {
    ath_error("[output_ppm:] Data must be 2D\n");
    return;
  }

/* construct output filename.  pOut->id will either be name of variable,
 * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
 * <output> block.  */
  if((fname = fname_construct(pGrid->outfilename,num_digit,dnum,pOut->id,"ppm"))
     == NULL){
    ath_error("[output_ppm]: Error constructing filename\n");
    return;
  }

/* open output file */
  if((pfile = fopen(fname,"w")) == NULL){
    ath_error("[output_ppm]: Unable to open ppm file %s\n",fname);
    return;
  }

/* Extract 2D data from 3D data,  Can either be slice or average along axis,
 * depending on range of ix1,ix2,ix3 in <ouput> block */

  data = subset2(pGrid,pOut);

/*Set the dimensions of the array corresponding to the sliced axis*/
/*A slice perpendicular to the ix3-axis*/
  if(pOut->Nx3 == 1){
    nx1 = pOut->Nx1;
    nx2 = pOut->Nx2;
  }
/*A slice perpendicular to the ix2-axis*/
  if(pOut->Nx2 == 1){
    nx1 = pOut->Nx1;
    nx2 = pOut->Nx3;
  }
/*A slice perpendicular to the ix1-axis*/
  if(pOut->Nx1 == 1){
    nx1 = pOut->Nx2;
    nx2 = pOut->Nx3;
  }

/* Store the global min / max, for output at end of run */
  minmax2(data,nx2,nx1,&dmin,&dmax);
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  fprintf(pfile,"P6\n");
  fprintf(pfile,"# dmin = %.7e, dmax = %.7e, gmin = %.7e, gmax = %.7e\n",
	  dmin,dmax,pOut->gmin,pOut->gmax);
  fprintf(pfile,"%d %d\n255\n",nx1,nx2);

/* Override auto-scale? */
  if (pOut->sdmin != 0) dmin = pOut->dmin;
  if (pOut->sdmax != 0) dmax = pOut->dmax;

  for (j=nx2-1; j>=0; j--) {
    for (i=0; i<nx1; i++) {
      compute_rgb(data[j][i],dmin,dmax,&red,&green,&blue,pOut);
      fputc(red,pfile);
      fputc(green,pfile);
      fputc(blue,pfile);
    }
  }

/* Close the file, free memory */
  fclose(pfile);
  free_2d_array((void **)data);
  free(fname);
}

/*----------------------------------------------------------------------------*/
/* compute_rgb: converts data into RGB values using palette in Output struct  */

static void compute_rgb(double data, double min, double max,
  int *R, int *G, int *B, Output *pOut)
{
  int i;
  float x, *rgb = pOut->rgb, *der = pOut->der;

  if (rgb == 0) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }

  if (min==max) {
    *R = *G = *B = (data > max ? 255 : 0);    
    return;
  }
#if 1
  x = (data-min)*255.0/(max-min);
  if (x<=0.0 || x>=255.0) {         /* out of bounds */
    i=  (x <= 0.0 ? 0 : 255);
    *R = (int) (rgb[i*3+0] * 255.0);
    *G = (int) (rgb[i*3+1] * 255.0);
    *B = (int) (rgb[i*3+2] * 255.0);
    return;
  }
  i = (int) x;
  *R = (int)  ((rgb[i*3+0] + (x-i)*der[i*3+0])*255.0);
  *G = (int)  ((rgb[i*3+1] + (x-i)*der[i*3+1])*255.0);
  *B = (int)  ((rgb[i*3+2] + (x-i)*der[i*3+2])*255.0);
#else
  i = (int) ((data-min)*255.0/(max-min));
  if (i<0) i=0;
  if (i>255) i=255;
  *R = (int) (rgb[i*3+0] * 255.0);
  *G = (int) (rgb[i*3+1] * 255.0);
  *B = (int) (rgb[i*3+2] * 255.0);
#endif
}
