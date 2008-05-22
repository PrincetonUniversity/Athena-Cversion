#include "copyright.h"
/*==============================================================================
 * FILE: output_pgm.c
 *
 * PURPOSE: Writes Portable Gray Map (PGM) outputs.  These are extremely simple
 *   grayscale 2D images, see e.g. sourceforge for documentation.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_pgm() -  outputs 2D PGM images
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* output_pgm: output 2D PGM image   */

void output_pgm(Grid *pGrid, Domain *pD, Output *pOut)
{
  FILE *pfile;
  char *fname;
  int nx1,nx2,gray,i,j;
  float **data, dmin, dmax, max_min, sfact;
  int dnum = pOut->num;

/* check output data is 2D (output must be a 2D slice for 3D runs) */
  if (pOut->ndim != 2) {
    ath_error("[output_pgm]: Output must be a 2D slice\n");
    return;
  }

/* construct output filename.  pOut->id will either be name of variable,
 * if 'id=...' was included in <ouput> block, or 'outN' where N is number of
 * <output> block.  */
  if((fname = ath_fname(NULL,pGrid->outfilename,num_digit,dnum,pOut->id,"pgm"))
     == NULL){
    ath_error("[output_pgm]: Error constructing filename\n");
    return;
  }

/* open output file */
  if((pfile = fopen(fname,"w")) == NULL){
    ath_error("[output_pgm]: Unable to open pgm file %s\n",fname);
    return;
  }

/* Extract 2D data from 3D data,  Can either be slice or average along axis,
 * depending on range of ix1,ix2,ix3 in <ouput> block */
  data = subset2(pGrid,pOut);
  nx1 = pOut->Nx1; /* we know it's a 2dim image */
  nx2 = pOut->Nx2;
  fprintf(pfile,"P5\n%d %d\n255\n",nx1,nx2);

/* Store the global min / max, for output at end of run */
  minmax2(data,nx2,nx1,&dmin,&dmax);
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

/* Override auto-scale? */
  if (pOut->sdmin != 0) dmin = pOut->dmin;
  if (pOut->sdmax != 0) dmax = pOut->dmax;

  max_min = (dmax - dmin)*(1.0 + FLT_EPSILON);

/* map the data which satisfies [min <= data <= max] to the range 
 * [0.0 , 256.0] -- Not inclusive of 256 */

  if(max_min > 0.0) {
    sfact = 256.0/max_min;
    for (j=nx2-1; j>=0; j--) {
      for (i=0; i<nx1; i++) {
/* Map the data to an 8 bit int, i.e. 0 -> 255 */
	gray = (int)(sfact*(data[j][i] - dmin));
/* Out of bounds data is mapped to the min or max integer value */
	gray = gray >   0 ? gray :   0;
	gray = gray < 255 ? gray : 255;

	fputc(gray,pfile);
      }
    }

/* else, if max=min set image to constant */

  } else {
    gray = 0;
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
	fputc(gray,pfile);
      }
    }
  }

/* Close the file, free memory */

  fclose(pfile); 
  free_2d_array(data);
  free(fname);
}
