#include "copyright.h"
/*==============================================================================
 * FILE: dump_binary.c
 *
 * PURPOSE: Function to write an unformatted dump of the field variables that
 *   can be read, e.g., by IDL scripts.  Also called by dump_dx().
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   dump_binary -
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* dump_binary:  */

void dump_binary(Grid *pGrid, Output *pOut)
{
  int dnum = pOut->num;
  FILE *p_binfile;
  char *fname;
  int n,ndata[4];
/* Upper and Lower bounds on i,j,k for data dump */
  int i, il = pGrid->is, iu = pGrid->ie;
  int j, jl = pGrid->js, ju = pGrid->je;
  int k, kl = pGrid->ks, ku = pGrid->ke;
  float eos[2],*data;
  Real *pq;

#ifdef WRITE_GHOST_CELLS
  if(pGrid->Nx1 > 1){
    il = pGrid->is - nghost;
    iu = pGrid->ie + nghost;
  }

  if(pGrid->Nx2 > 1){
    jl = pGrid->js - nghost;
    ju = pGrid->je + nghost;
  }

  if(pGrid->Nx3 > 1){
    kl = pGrid->ks - nghost;
    ku = pGrid->ke + nghost;
  }
#endif /* WRITE_GHOST_CELLS */

  if((fname = fname_construct(pGrid->outfilename,num_digit,dnum,NULL,"bin")) 
     == NULL){
    ath_error("[dump_binary]: Error constructing filename\n");
    return;
  }

  if((p_binfile = fopen(fname,"wb")) == NULL){
    ath_error("[dump_binary]: Unable to open binary dump file\n");
    return;
  }

/* Write number of zones and variables */
  ndata[0] = iu-il+1;
  ndata[1] = ju-jl+1;
  ndata[2] = ku-kl+1;
  ndata[3] = NVAR;
  fwrite(ndata,sizeof(int),4,p_binfile);

/* Write (gamma-1) and isothermal sound speed */

#ifdef ISOTHERMAL
  eos[0] = (float)0.0;
  eos[1] = (float)Iso_csound;
#else
  eos[0] = (float)Gamma_1 ;
  eos[1] = (float)0.0;
#endif
  fwrite(eos,sizeof(float),2,p_binfile);

/* Allocate Memory */

  if((data = (float *)malloc(ndata[0]*sizeof(float))) == NULL){
    ath_error("[dump_binary]: malloc failed for temporary array\n");
    fclose(p_binfile);
    return;
  }

/* Write cell-centered date in Gas array pGrid->U[n] */

  for (n=0;n<NVAR; n++) {
    for (k=0; k<ndata[2]; k++) {
    for (j=0; j<ndata[1]; j++) {
      for (i=0; i<ndata[0]; i++) {
        pq = ((Real *) &(pGrid->U[k+kl][j+jl][i+il])) + n;
        data[i] = (float)(*pq);
      }
      fwrite(data,sizeof(float),(size_t)ndata[0],p_binfile);
    }}
  }

/* close file and free memory */
  fclose(p_binfile); 
  free(data); 
  free(fname);
}
