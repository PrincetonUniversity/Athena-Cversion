#include "copyright.h"
/*==============================================================================
 * FILE: output_tab.c
 *
 * PURPOSE: Output formatted tables of data.
 *
 * Writes only one single data point.  Needs to write 1D slice, and N vars.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_tab() 
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* output_tab:  */

void output_tab(Grid *pGrid, Output *pOut)
{
  FILE *pfile;
  char *fname;
  int ndata[3],i,iu,il,j,ju,jl;
  float data, dmin, dmax, xworld, yworld;
  int dnum = pOut->num;

#ifdef WRITE_GHOST_CELLS
  if(pGrid->Nx1 > 1){
    il = pGrid->is - nghost;
    iu = pGrid->ie + nghost;
  }
  else{
    il = pGrid->is;
    iu = pGrid->ie;
  }

  if(pGrid->Nx2 > 1){
    jl = pGrid->js - nghost;
    ju = pGrid->je + nghost;
  }
  else{
    jl = pGrid->js;
    ju = pGrid->je;
  }
#else
  il = pGrid->is;
  iu = pGrid->ie;
  jl = pGrid->js;
  ju = pGrid->je;
#endif

/* construct output filename */
  fname = fname_construct(pGrid->outfilename,num_digit,dnum,pOut->id,"tab");
  if (fname == NULL) {
    ath_error("[output_tab]: Error constructing filename\n");
    return;
  }

/* open output file */
  if ((pfile = fopen(fname,"w")) == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
    return;
  }

  ndata[0] = iu-il+1;
  ndata[1] = ju-jl+1;


  for (j=0; j<ndata[1]; j++) {              /* construct the data matrix */
    yworld = pGrid->x2_0 +  pGrid->dx2*(j-pGrid->js);
    for (i=0; i<ndata[0]; i++) {
      xworld = pGrid->x1_0  + pGrid->dx1*(i-pGrid->is);
      data = (*pOut->expr)(pGrid,i+il,j+jl,pGrid->ks);
      fprintf(pfile,"%g %g  %g\n",xworld,yworld,data);
      if (i>0 || j>0) {
	dmin = MIN(dmin,data);
	dmax = MAX(dmax,data);
      } else
	dmin = dmax = data;
    }
  }

  fclose(pfile);
  
  if (pOut->num == 0) {
    pOut->gmin = dmin;
    pOut->gmax = dmax;
  } else {
    pOut->gmin = MIN(dmin,pOut->gmin);
    pOut->gmax = MAX(dmax,pOut->gmax);
  }

  free(fname);
}
