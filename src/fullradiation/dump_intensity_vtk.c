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
void dump_ix1_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    il = pRG->is;
	jl = pRG->js, ju = pRG->je;
	kl = pRG->ks, ku = pRG->ke;
	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    il -=1;
          
	if(pRG->Nx[1] > 1) {
	  jl -= 1; ju += 1;
	}
	if(pRG->Nx[2] > 1) {
	  kl -= 1; ku += 1;
        }
#endif /* WRITE_GHOST_CELLS */
	
	ndata1 = ju-jl+1;
	ndata2 = ku-kl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ix1 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
          ath_error("[dump_ix1_vtk]: Error constructing basename\n");
	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ix1");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ix1_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ix1_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata1*sizeof(float))) == NULL){
          ath_error("[dump_ix1_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
	x1 = pRG->MinX[0];
        x2 = pRG->MinX[1];
        x3 = pRG->MinX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
        if (pRG->Nx[1] == 1) {
          fprintf(pfile,"DIMENSIONS %d %d %d\n",1,1,1);
        } else {
          if (pRG->Nx[2] == 1) {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",1,ju-jl+2,1);
          } else {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",1,ju-jl+2,ku-kl+2);
          }
        }
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (ju-jl+1)*(ku-kl+1));

      
/* Write left intensities for each frequency and angle */

	for (ifr=0; ifr<nf; ifr++) {
	  for(l=0; l<noct; l++) {
/* Only write outgoing intensities */
	    if((l == 1) || (l == 3) || (l == 5) || (l == 7)) {
	      for(m=0; m<nang; m++) {
              n=nang * l + m;
		fprintf(pfile,"SCALARS l1imu_%i_%i_%i float\n",ifr,l,m);
		fprintf(pfile,"LOOKUP_TABLE default\n");
		for (k=kl; k<=ku; k++) {
		  for (j=jl; j<=ju; j++) {
		    data[j-jl] = (float)pRG->imu[k][j][il][ifr][n];
		  }
		  if(!big_end) ath_bswap(data,sizeof(float),ju-jl+1);
		  fwrite(data,sizeof(float),(size_t)ndata1,pfile);
		}
	      }
	    }
	  }}
/* close file and free memory */
        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* \fn dump_ox1_vtk(MeshS *pM, OutputS *pOut)  */
void dump_ox1_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    iu = pRG->ie;
	jl = pRG->js, ju = pRG->je;
	kl = pRG->ks, ku = pRG->ke;
	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    iu += 1;
	if(pRG->Nx[1] > 1) {
	  jl -= 1; ju += 1;
	}
	if(pRG->Nx[2] > 1) {
	  kl -= 1; ku += 1;
        }
#endif /* WRITE_GHOST_CELLS */
	
	ndata1 = ju-jl+1;
	ndata2 = ku-kl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ox1 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
          ath_error("[dump_ox1_vtk]: Error constructing basename\n");

	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ox1");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ox1_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ox1_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata1*sizeof(float))) == NULL){
          ath_error("[dump_ox1_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
	x1 = pRG->MaxX[0];
        x2 = pRG->MinX[1];
        x3 = pRG->MinX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
        if (pRG->Nx[1] == 1) {
          fprintf(pfile,"DIMENSIONS %d %d %d\n",1,1,1);
        } else {
          if (pRG->Nx[2] == 1) {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",1,ju-jl+2,1);
          } else {
            fprintf(pfile,"DIMENSIONS %d %d %d\n",1,ju-jl+2,ku-kl+2);
          }
        }
	fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (ju-jl+1)*(ku-kl+1));

      
/* Write right intensities for each frequency and angle */


	for (ifr=0; ifr<nf; ifr++) {
	  for(l=0; l<noct; l++) {
/* Only write outgoing intensities */
	    if((l == 0) || (l == 2) || (l == 4) || (l == 6)) {
	      for(m=0; m<nang; m++) {
              n = nang * l + m;
		fprintf(pfile,"SCALARS r1imu_%i_%i_%i float\n",ifr,l,m);
		fprintf(pfile,"LOOKUP_TABLE default\n");
		for (k=kl; k<=ku; k++) {
		  for (j=jl; j<=ju; j++) {
		    data[j-jl] = (float)pRG->imu[k][j][iu][ifr][n];
		  }
		  if(!big_end) ath_bswap(data,sizeof(float),ju-jl+1);
		  fwrite(data,sizeof(float),(size_t)ndata1,pfile);
		}
	      }
	    }
	  }}
/* close file and free memory */
        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* \fn dump_ix2_vtk(MeshS *pM, OutputS *pOut)  */
void dump_ix2_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    jl = pRG->js;
	il = pRG->is, iu = pRG->ie;
	kl = pRG->ks, ku = pRG->ke;
	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    jl -= 1;
	il -= 1; iu += 1;
	if(pRG->Nx[2] > 1) {
	  kl -= 1; ku += 1;
        }
#endif /* WRITE_GHOST_CELLS */
	
	ndata0 = iu-il+1;
	ndata2 = ku-kl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ix2 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
	  ath_error("[dump_ix2_vtk]: Error constructing basename\n");
	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ix2");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ix2_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ix2_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata0*sizeof(float))) == NULL){
          ath_error("[dump_ix2_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

        x1 = pRG->MinX[0];
	x2 = pRG->MinX[1];
        x3 = pRG->MinX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
	if (pRG->Nx[2] == 1) {
	  fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,1);
	} else {
	  fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,ku-kl+2);
	}
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ku-kl+1));

      
/* Write left intensities for each frequency and angle */

	for (ifr=0; ifr<nf; ifr++) {
	  for(l=0; l<noct; l++) {
/* Only write outgoing intensities */
	    if((l == 2) || (l == 3) || (l == 6) || (l == 7)) {
	      for(m=0; m<nang; m++) {
              n=l*nang + m;
		fprintf(pfile,"SCALARS l2imu_%i_%i_%i float\n",ifr,l,m);
		fprintf(pfile,"LOOKUP_TABLE default\n");
		for (k=kl; k<=ku; k++) {
		  for (i=il; i<=iu; i++) {
		    data[i-il] = (float)pRG->imu[k][jl][i][ifr][n];
		  }
		  if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
		  fwrite(data,sizeof(float),(size_t)ndata0,pfile);
		}
	      }
	    }
	  }}	 
/* close file and free memory */
        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* \fn dump_ox2_vtk(MeshS *pM, OutputS *pOut)  */
void dump_ox2_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    ju = pRG->je;
	il = pRG->is, iu = pRG->ie;
	kl = pRG->ks, ku = pRG->ke;
	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    ju += 1;
	il -= 1; iu += 1;
	if(pRG->Nx[2] > 1) {
	  kl -= 1; ku += 1;
        }
#endif /* WRITE_GHOST_CELLS */
	
	ndata0 = iu-il+1;
	ndata2 = ku-kl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ox2 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
	  ath_error("[dump_ox2_vtk]: Error constructing basename\n");
	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ox2");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ox2_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ox2_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata0*sizeof(float))) == NULL){
          ath_error("[dump_ox2_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */

        x1 = pRG->MinX[0];
	x2 = pRG->MaxX[1];
        x3 = pRG->MinX[2];

	fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
	if (pRG->Nx[2] == 1) {
	  fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,1);
	} else {
	  fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,1,ku-kl+2);
	}
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);
/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ku-kl+1));

      
/* Write right intensities for each frequency and angle */

	for (ifr=0; ifr<nf; ifr++) {
	  for(l=0; l<noct; l++) {
/* Only write outgoing intensities */
	    if((l == 0) || (l == 1) || (l == 4) || (l == 5)) {
	      for(m=0; m<nang; m++) {
              n = l * nang + m;
		fprintf(pfile,"SCALARS r2imu_%i_%i_%i float\n",ifr,l,m);
		fprintf(pfile,"LOOKUP_TABLE default\n");
		for (k=kl; k<=ku; k++) {
		  for (i=il; i<=iu; i++) {
		    data[i-il] = (float)pRG->imu[k][ju][i][ifr][n];
		  }
		  if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
		  fwrite(data,sizeof(float),(size_t)ndata0,pfile);
		}
	      }
	    }
	  }}	    
/* close file and free memory */
        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* \fn dump_ix3_vtk(MeshS *pM, OutputS *pOut)  */
void dump_ix3_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    kl = pRG->ks;
	il = pRG->is, iu = pRG->ie;
	jl = pRG->js, ju = pRG->je;

	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    kl -= 1;
	il -= 1; iu += 1;
	jl -= 1; ju += 1;
#endif /* WRITE_GHOST_CELLS */

	ndata0 = iu-il+1;	
	ndata1 = ju-jl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ix3 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
          ath_error("[dump_ix3_vtk]: Error constructing basename\n");
	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ix3");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ix3_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ix3_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata0*sizeof(float))) == NULL){
          ath_error("[dump_ix3_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
	x1 = pRG->MinX[0];
        x2 = pRG->MinX[1];
	x3 = pRG->MinX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
	fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,1);
        fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ju-jl+1));
      
/* Write left intensities for each frequency and angle */
	for (ifr=0; ifr<nf; ifr++) {
/* Only write outgoing intensities */
	  for(l=4; l<=7; l++) {
	    for(m=0; m<nang; m++) {
            n=l * nang +m;
	      fprintf(pfile,"SCALARS l3imu_%i_%i_%i float\n",ifr,l,m);
	      fprintf(pfile,"LOOKUP_TABLE default\n");
	      for (j=jl; j<=ju; j++) {
		for (i=il; i<=iu; i++) {
		  data[i-il] = (float)pRG->imu[kl][j][i][ifr][n];
		}
		if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
		fwrite(data,sizeof(float),(size_t)ndata0,pfile);
	      }
	    }}}	    
/* close file and free memory */

        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/* \fn dump_ox3_vtk(MeshS *pM, OutputS *pOut)  */
void dump_ox3_vtk(MeshS *pM, OutputS *pOut)
{
  RadGridS *pRG;
  FILE *pfile;
  char *fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
/* Upper and Lower bounds on i,j,k for data dump */
  int i,j,k,il,iu,jl,ju,kl,ku,nl,nd,l,m,n,ifr;
  int big_end = ath_big_endian();
  int ndata0,ndata1,ndata2;
  float *data;   /* points to 3*ndata0 allocated floats */
  double x1, x2, x3, time;
  int nf, nang, noct;
  char *basename = NULL;
  int slen;

/* Loop over all Domains in Mesh, and output Grid data */

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].RadGrid != NULL){
	time = pM->Domain[nl][nd].Grid->time;

/* write files if domain and level match input, or are not specified (-1) */
      if ((pOut->nlevel == -1 || pOut->nlevel == nl) &&
          (pOut->ndomain == -1 || pOut->ndomain == nd)){
        pRG = pM->Domain[nl][nd].RadGrid;

    ku = pRG->ke;
	il = pRG->is, iu = pRG->ie;
	jl = pRG->js, ju = pRG->je;

	nf = pRG->nf;
	nang = pRG->nang;
	noct = pRG->noct;

#ifdef WRITE_GHOST_CELLS
    ku += 1;
	il -= 1; iu += 1;
	jl -= 1; ju += 1;
#endif /* WRITE_GHOST_CELLS */

	ndata0 = iu-il+1;	
	ndata1 = ju-jl+1;

/* construct filename, open file */
        if (nl>0) {
          plev = &levstr[0];
          sprintf(plev,"lev%d",nl);
        }
        if (nd>0) {
          pdom = &domstr[0];
          sprintf(pdom,"dom%d",nd);
	}
/* Add .ox3 to basename */
	slen = 4 + strlen(pM->outfilename);
	if((basename = (char*)malloc(slen*sizeof(char))) == NULL){
          ath_error("[dump_ox3_vtk]: Error constructing basename\n");
	}
	strcpy(basename,pM->outfilename);
	strcpy(&(basename[slen-4]),".ox3");

        if((fname = ath_fname(plev,basename,plev,pdom,num_digit,
            pOut->num,NULL,"vtk")) == NULL){
          ath_error("[dump_ox3_vtk]: Error constructing filename\n");
        }
	free(basename);

        if((pfile = fopen(fname,"w")) == NULL){
          ath_error("[dump_ox3_vtk]: Unable to open vtk dump file\n");
          return;
        }

/* Allocate memory for temporary array of floats */

        if((data = (float *)malloc(ndata0*sizeof(float))) == NULL){
          ath_error("[dump_ox3_vtk]: malloc failed for temporary array\n");
          return;
        }

/* There are five basic parts to the VTK "legacy" file format.  */
/*  1. Write file version and identifier */

        fprintf(pfile,"# vtk DataFile Version 2.0\n");

/*  2. Header */


	fprintf(pfile,"intensities vars at time= %e, level= %i, domain= %i\n",
            time,nl,nd); 

/*  3. File format */

        fprintf(pfile,"BINARY\n");

/*  4. Dataset structure */

/* Set the Grid origin */
	x1 = pRG->MinX[0];
        x2 = pRG->MinX[1];
	x3 = pRG->MaxX[2];

        fprintf(pfile,"DATASET STRUCTURED_POINTS\n");
	fprintf(pfile,"DIMENSIONS %d %d %d\n",iu-il+2,ju-jl+2,1);
	fprintf(pfile,"ORIGIN %e %e %e \n",x1,x2,x3);
        fprintf(pfile,"SPACING %e %e %e \n",pRG->dx1,pRG->dx2,pRG->dx3);

/*  5. Data  */

        fprintf(pfile,"CELL_DATA %d \n", (iu-il+1)*(ju-jl+1));
      
/* Write right intensities for each frequency and angle */
	for (ifr=0; ifr<nf; ifr++) {
/* Only write outgoing intensities */
	  for(l=0; l<=3; l++) {
	    for(m=0; m<nang; m++) {
            n=l * nang + m;
	      fprintf(pfile,"SCALARS r3imu_%i_%i_%i float\n",ifr,l,m);
	      fprintf(pfile,"LOOKUP_TABLE default\n");
	      for (j=jl; j<=ju; j++) {
		for (i=il; i<=iu; i++) {
		  data[i-il] = (float)pRG->imu[ku][j][i][ifr][n];
		}
		if(!big_end) ath_bswap(data,sizeof(float),iu-il+1);
		fwrite(data,sizeof(float),(size_t)ndata0,pfile);
	      }
	    }}}	 	    
/* close file and free memory */

        fclose(pfile);
        free(data);
      }}
    }
  }
  return;
}

#endif /* FULL_RADIATION_TRANSFER */


