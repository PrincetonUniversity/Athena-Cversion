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
 *                   in ascii 
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* \fn void output_rad_mesh(MeshS *pM)
 *  Writes ascii table with information about the the frequency and angular
 *  grids and corresponding quadratures.  Called by data_output(). */
void output_rad_mesh(MeshS *pM)
{

  RadGridS *pRG;
  DomainS *pD;
  FILE *pfile;
  char *fname, *plev=NULL, *pdom=NULL;
  char levstr[8], domstr[8];
  int nl, nd, i, j;
  int nf, noct, nang;
  int mhst, myID_Comm_Domain=1;
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
				0,NULL,"rad")) == NULL){
	    ath_error("[output_rad_mesh]: Error constructing filename\n");
	  }
	  if((pfile = fopen(fname,"w")) == NULL){
	    ath_error("[output_rad_mesh]: Unable to open file\n");
	    return;
	  }

	  if ( (radt_mode == 0) || (radt_mode == 2) ) {
	    pRG = pM->Domain[nl][nd].RadGrid;
	    nf = pRG->nf;
	    nang = pRG->nang;
	    noct = pRG->noct;
	    fprintf(pfile,"RadGrid:\n");
	    fprintf(pfile,"%d %d %d\n",nf,noct,nang);
	    for (i=0; i<nf; i++) {
	      fprintf(pfile,"%d %12.8e %12.8e\n",i,pRG->nu[i],pRG->wnu[i]);
	    }
	    for (i=0; i<nang; i++) {
	      fprintf(pfile,"%d ",i);
	      for(j=0; j<3; j++) 
		fprintf(pfile,"%12.8e ",pRG->mu[0][i][j]);
	      fprintf(pfile,"%12.8e\n",pRG->wmu[i]);
	    }
	  }
	  if ( (radt_mode == 1) || (radt_mode == 2) ) {
	    pRG = pM->Domain[nl][nd].RadOutGrid;
	    nf = pRG->nf;
	    nang = pRG->nang;
	    noct = pRG->noct;
	    fprintf(pfile,"RadOutGrid:\n");
	    fprintf(pfile,"%d %d %d\n",nf,noct,nang);
	    for (i=0; i<nf; i++) {
	      fprintf(pfile,"%d %12.8e %12.8e\n",i,pRG->nu[i],pRG->wnu[i]);
	    }
	    for (i=0; i<nang; i++) {
	      fprintf(pfile,"%d ",i);
	      for(j=0; j<3; j++) 
		fprintf(pfile,"%12.8e ",pRG->mu[0][i][j]);
	      fprintf(pfile,"%12.8e\n",pRG->wmu[i]);
	    }
	  }
	  /* close file */
	  fclose(pfile);	
	}
      }
    }}

  return;
}

#endif
