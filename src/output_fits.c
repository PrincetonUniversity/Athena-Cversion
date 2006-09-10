#include "copyright.h"
/*==============================================================================
 * FILE: output_fits.c
 *
 * PURPOSE: Functions for writing FITS files.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_fits() - opens file and calls appropriate 2D/3D output function
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   output_fits_2d() - write fits file for 2D data
 *   output_fits_3d() - write fits file for 3D data
 *============================================================================*/

void output_fits_2d(Grid *pGrid, Output *pOut, FILE *pFile);
void output_fits_3d(Grid *pGrid, Output *pOut, FILE *pFile);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* output_fits:  open file, call 2D/3D writer; called by data_ouput  */

void output_fits(Grid *pGrid, Domain *pD, Output *pOut)
{
  FILE *pFile;
  char *fname;
  int dnum = pOut->num;

/* construct output filename */
  fname = fname_construct(pGrid->outfilename,num_digit,dnum,pOut->id,"fits");
  if (fname == NULL) {
    ath_error("[output_fits]: Error constructing output filename\n");
    return;
  }

/* open filename */
  pFile = fopen(fname,"w");
  if (pFile == NULL) {
    ath_error("[output_fits]: Unable to open fits file %s\n",fname);
    return;
  }

/* call appropriate 2D or 3D output routine */
  if (pOut->ndim == 3) {
    output_fits_3d(pGrid,pOut,pFile);
  } else if (pOut->ndim == 2) {
    output_fits_2d(pGrid,pOut,pFile);
  } else {
    ath_error("output_fits: ouput data must be 2D or 3D\n");
    return;
  }

  fclose(pFile);
  free(fname);
  return;
}

/*----------------------------------------------------------------------------*/
/* output_fits_2d: writes 2D data  */

void output_fits_2d(Grid *pGrid, Output *pOut, FILE *pFile)
{
  int n, nx1, nx2;
  float **data, dmin, dmax, zero = 0.0;
  int dnum = pOut->num;
  int recsize, ncard = 0;
  char card[81];

  data = subset2(pGrid,pOut);
  nx1 = pOut->Nx1;     /* we know it's a 2dim image */
  nx2 = pOut->Nx2;
  minmax2(data,nx2,nx1,&dmin,&dmax);

  ncard = 0;
  sprintf(card,"SIMPLE  = T");
  fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"BITPIX  = -32");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS   = 2");             fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS1  = %d",nx1);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS2  = %d",nx2);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRPIX1  = %g",1.0);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRPIX2  = %g",1.0);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CDELT1  = %g",pGrid->dx1); fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CDELT2  = %g",pGrid->dx2); fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRVAL1  = %g",pGrid->x1_0);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRVAL2  = %g",pGrid->x2_0);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CTYPE1  = 'X'");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CTYPE2  = 'Y'");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"DATAMIN = %g",dmin);       fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"DATAMAX = %g",dmax);       fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"TIME    = %g",pGrid->time);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"END");                     fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card," ");

/* better make sure we never wrote more than 2880/80 here ... */
  for ( ; ncard < 2880/80; ncard++)   /* write the trailing end of the header */
    fprintf(pFile,"%-80s",card);

  n = nx1*nx2;
/* careful: we're swapping bytes, don't use data anymore */
  if (!ath_big_endian())     
    ath_bswap(&data[0][0],sizeof(float),n);
  ath_fwrite(&data[0][0],n,sizeof(float),pFile);

  recsize = 2880 / sizeof(float);   /* this never has a remainder */
  n %= recsize;                   /* check if any trailing data to be written */
  if (n) {                /* if so, write a bunch of zero's until 2880 filled */
    n = recsize - n;
    zero = 0.0;
    ath_fwrite(&zero,n,sizeof(float),pFile);
  }
  
/* Compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free_2d_array((void **)data); /* Free the memory we malloc'd */
}


/*----------------------------------------------------------------------------*/
/* output_fits_3d: writes 3D data   */

void output_fits_3d(Grid *pGrid, Output *pOut, FILE *pFile)
{
  int ndata[3],i,j,k,n;
  float ***data, dmin, dmax, zero = 0.0;
  int dnum = pOut->num;
  int recsize, ncard = 0;
  char card[81];

/* no re-ordering of axes for 3d->3d mapping - for now */
  ndata[0] = pOut->Nx1;     
  ndata[1] = pOut->Nx2;
  ndata[2] = pOut->Nx3;

  data = subset3(pGrid,pOut);
  minmax3(data, ndata[2], ndata[1], ndata[0], &dmin, &dmax);

  ncard = 0;
  sprintf(card,"SIMPLE  = T");             fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"BITPIX  = -32");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS   = 3");             fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS1  = %d",ndata[0]);   fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS2  = %d",ndata[1]);   fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"NAXIS3  = %d",ndata[2]);   fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRPIX1  = %g",1.0);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRPIX2  = %g",1.0);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRPIX3  = %g",1.0);        fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CDELT1  = %g",pGrid->dx1); fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CDELT2  = %g",pGrid->dx2); fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CDELT3  = %g",pGrid->dx2); fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRVAL1  = %g",pGrid->x1_0);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRVAL2  = %g",pGrid->x2_0);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CRVAL3  = %g",pGrid->x3_0);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CTYPE1  = 'X'");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CTYPE2  = 'Y'");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"CTYPE3  = 'Z'");           fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"TIME    = %g",pGrid->time);fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"DATAMIN = %g",dmin);       fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"DATAMAX = %g",dmax);       fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card,"END");                     fprintf(pFile,"%-80s",card); ncard++;
  sprintf(card," ");

/* better make sure we never wrote more than 2880/80 here ... */
  for ( ; ncard < 2880/80; ncard++)   /* write the trailing end of the header */
    fprintf(pFile,"%-80s",card);


  if (pOut->sdmin != 0 || pOut->sdmax != 0) { /* override auto-scale */
    if (pOut->sdmin != 0) dmin = pOut->dmin;
    if (pOut->sdmax != 0) dmax = pOut->dmax;
/* construct the data matrix */
    for (k=0; k<ndata[2]; k++) {      
      for (j=0; j<ndata[1]; j++) { 
	for (i=0; i<ndata[0]; i++) {
	  if (data[k][j][i] > dmax) data[k][j][i] = dmax;
	  if (data[k][j][i] < dmin) data[k][j][i] = dmin;
	}
      }
    }
  }

  n = ndata[0]*ndata[1]*ndata[2];
/* careful: we're swapping bytes, don't use data anymore */
  if (!ath_big_endian())      
    ath_bswap(&data[0][0][0],sizeof(float),n);
  ath_fwrite(&data[0][0][0],n,sizeof(float),pFile);

  recsize = 2880 / sizeof(float);   /* this never has a remainder */
  n %= recsize;                   /* check if any trailing data to be written */
  if (n) {                /* if so, write a bunch of zero's until 2880 filled */
    n = recsize - n;
    zero = 0.0;
    ath_fwrite(&zero,n,sizeof(float),pFile);
  }
  
/* compute and store global min/max, for output at end of run */
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free_3d_array((void ***)data); /* Free the memory we malloc'd */
}
