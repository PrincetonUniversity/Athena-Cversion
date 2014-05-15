#include "../copyright.h"
/*==============================================================================
 * FILE: output_intensity_vtk.c
 *
 * PURPOSE: Function to write a dump of the specific intensity on specific
 *   faces in VTK "legacy" format.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   dump_???_vtk() - writes VTK dump of intensity at ?? boundary
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"

#ifdef FULL_RADIATION_TRANSFER

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/* \fn dump_ix1_vtk(MeshS *pM, OutputS *pOut)  */
void dump_intensity_vtk(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  RadGridS *pRG;

  int irl,iru,jrl,jru,krl,kru, i, j, k;
  int ioff, joff, koff;
  Real x1, x2, x3;

  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  int slen;
  char *basename = NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int nl,nd;
  int ifr, nf, n, l, nang, noct;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */




/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){

/* write files if domain and level match input, or are not specified (-1) */
        if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
            (pOut->ndomain == -1 || pOut->ndomain == nd)){

          pG = pM->Domain[nl][nd].Grid;

          pRG = pM->Domain[nl][nd].RadGrid;
          nf = pRG->nf;

          nang = pRG->nang;
          noct = pRG->noct;

          irl = pRG->is, iru = pRG->ie;
          jrl = pRG->js, jru = pRG->je;
          krl = pRG->ks, kru = pRG->ke;

          ioff = nghost - Radghost;

          joff = 0;
          if(pRG->Nx[1] > 1) {
            joff = ioff;
          }

          koff = 0;
          if(pRG->Nx[2] > 1) {
            koff = ioff;
          }



#ifdef WRITE_GHOST_CELLS
          iru += Radghost;
          irl -= Radghost;

          if(pRG->Nx[1] > 1) {
            jru += Radghost;
            jrl -= Radghost;
          }

          if(pRG->Nx[2] > 1) {
            kru += Radghost;
            krl -= Radghost;
          }
#endif /* WRITE_GHOST_CELLS */

          ndata0 = iru-irl+1;
          ndata1 = jru-jrl+1;
          ndata2 = kru-krl+1;


/* calculate primitive variables, if needed */

/* construct filename, open file */
          if (nl>0) {
            plev = &levstr[0];
            sprintf(plev,"lev%d",nl);
          }

          if (nd>0) {
            pdom = &domstr[0];
            sprintf(pdom,"dom%d",nd);
          }
/* Add .intensity to basename */

          slen = 10 + strlen(pM->outfilename);
          if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
            ath_error("[output_ix1_vtk]: Error constructing basename\n");
          }
          strcpy(basename,pM->outfilename);
          strcpy(&(basename[slen-10]),".intensity");


          if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
                                pOut->num,NULL,"vtk")) == NULL){
            ath_error("[dump_vtk]: Error constructing filename\n");
          }

          free(basename);


          if((pfile = fopen(fname,"w")) == NULL){
            ath_error("[dump_vtk]: Unable to open vtk dump file\n");
            return;
          }
          free(fname);

/* Allocate memory for temporary array of floats */

          if((data = (float *)malloc(2*ndata0*sizeof(float))) == NULL){
            ath_error("[dump_vtk]: malloc failed for temporary array\n");
            return;
          }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

          fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


          fprintf(pfile,"Intensity vars at time= %e, level= %i, domain= %i\n",
                  pG->time,nl,nd);


/*  3. File format */

          fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

          fc_pos(pG, irl+ioff, jrl+joff, krl+koff, &x1, &x2, &x3);

          fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
          if (pRG->Nx[1] == 1) {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",iru-irl+2,1,1);
          } else {
            if (pRG->Nx[2] == 1) {
              fprintf(pfile,"DIMENSIONS %d %d %d\n",iru-irl+2,jru-jrl+2,1);
            } else {
              fprintf(pfile,"DIMENSIONS %d %d %d\n",iru-irl+2,jru-jrl+2,kru-krl+2);
            }
          }
          fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
          fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

          fprintf(pfile,"CELL_DATA %d \n", (iru-irl+1)*(jru-jrl+1)*(kru-krl+1));

/* Write intensity for every ifr, l, n */
          for(ifr=0; ifr<nf; ifr++){
            for(l=0; l<noct; l++){
              for(n=0; n<nang; n++){

                fprintf(pfile,"SCALARS imu_%i_%i_%i float\n",ifr,l,n);
                fprintf(pfile,"LOOKUP_TABLE default\n");

                for (k=krl; k<=kru; k++) {
                  for (j=jrl; j<=jru; j++) {
                    for (i=irl; i<=iru; i++) {
                      data[i-irl] = (float)pRG->imu[k][j][i][ifr][nang * l + n];
                    }
                    if(!big_end) ath_bswap(data,sizeof(float),iru-irl+1);
                    fwrite(data,sizeof(float),(size_t)ndata0,pfile);
                  }
                }


              }/* end n */
            }/* end l */
          }/* end ifr */



/* close file and free memory */

          fclose(pfile);
          free(data);
        }/* ==nl and ==nd */
      }/* Grid != NULL */
    }/* end nl */
  }/* end nd */

  return;
}

#endif /* FULL_RADIATION_TRANSFER */


