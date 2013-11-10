#include "../copyright.h"
/*==============================================================================
 * FILE: output_spec.c
 *
 * PURPOSE: Functions for writing output of radiation grid information:
 * angles and frequencies to be used with intensity dumps and outputs.
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_spec() - opens file and writes frequency grid and angular grid
 *   in ascii 
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER

void output_spec(MeshS *pM)
{

  RadGridS *pRG;
  DomainS *pD;
  FILE *pfile;
  char *fname, *plev=NULL, *pdom=NULL;
  char levstr[8], domstr[8];
  int nl, nd, i, j;
  int nf, noct, nang;
  int myID_Comm_Domain=1;
#ifdef MPI_PARALLEL
  int ierr;
#endif
/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	pD = (DomainS*)&(pM->Domain[nl][nd]);
	/* Only the parent (rank=0) process writes output.
	 * For single-processor jobs, myID_Comm_world is always zero. */
#ifdef MPI_PARALLEL
        ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
#endif
	if((myID_Comm_Domain==0) || (myID_Comm_world==0)){ 
	
/* write files if domain and level match input, or are not specified (-1) */
	  pRG = pM->Domain[nl][nd].RadGrid;
	  nf = pRG->nf;
	  nang = pRG->nang;
	  noct = pRG->noct;

/* construct filename and  open file */
	  if (nl>0) {
	    plev = &levstr[0];
	    sprintf(plev,"lev%d",nl);
	  }
	  if (nd>0) {
	    pdom = &domstr[0];
	    sprintf(pdom,"dom%d",nd);
	  }
	  if((fname = ath_fname(plev,pM->outfilename,plev,pdom,0,
				0,NULL,"spc")) == NULL){
	    ath_error("[output_spec]: Error constructing filename\n");
	  }
	  if((pfile = fopen(fname,"w")) == NULL){
	    ath_error("[output_spec]: Unable to open file\n");
	    return;
	  }
	  fprintf(pfile,"%d %d %d\n",nf,noct,nang);
	  for (i=0; i<nf; i++) {
	    fprintf(pfile,"%d %12.8e %12.8e\n",i,pRG->nu[i],pRG->wnu[i]);
	  }
	  for (i=0; i<nang; i++) {
	    fprintf(pfile,"%d ",i);
	    for(j=0; j<3; j++) 
	      fprintf(pfile,"%12.8e ",pRG->mu[0][0][0][i][j]);
	    fprintf(pfile,"%12.8e\n",pRG->wmu[0][0][0][i]);
	  }
	  /* close file */
	  fclose(pfile);	
	}
      }
    }}

  return;
}

#endif
