#include "copyright.h"
/*==============================================================================
 * FILE: output_pdf.c
 *
 * PURPOSE: Outputs Probability Distribution Functions of selected variables
 *   in formatted tabular form.  Fully MPI enabled, which requires passing
 *   lots of global sums and means (only the parent process produces output).
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   output_pdf() -
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

static int size_dat=0; /* Number of elements in the data[] array */
static double *data=NULL; /* Computed data array: data[size_dat] */

static int size_pdf=0; /* Number of elements in the pdf[] array */
static int *pdf=NULL; /* (non-normalized) PDF */
#ifdef MPI_PARALLEL
static int *cg_pdf=NULL; /* (non-normalized) complete grid PDF */
#endif /* MPI_PARALLEL */

static char def_fmt[]="%21.15e"; /* A default tabular dump data format */

/*----------------------------------------------------------------------------*/
/* output_pdf:   */

void output_pdf(Grid *pG, Domain *pD, Output *pout)
{
  FILE *pfile;
  char fmt[80];
  char fid[80]; /* File "id" for the statistics table */
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int n, data_cnt;
  double dmin, dmax, delta, dpdf, dat, scl;
  double mean=0.0, var=0.0; /* mean and variance of the distribution */
  double adev=0.0, sdev=0.0; /* average & standard deviation */
  double skew=0.0, kurt=0.0; /* skewness and kurtosis of the distribution */
  double r, s, ep=0.0; /* Temp. variables for calculating the variance, etc. */
#ifdef MPI_PARALLEL
  int err, cg_data_cnt;
  cg_data_cnt = (pD->ixe - pD->ixs +1)*(pD->jxe - pD->jxs +1)*(pD->kxe - pD->kxs +1);
#endif /* MPI_PARALLEL */

  data_cnt = (ie - is + 1)*(je - js + 1)*(ke - ks + 1);

/* Are the requisite arrays allocated? */
  if(data == NULL){
    size_dat = data_cnt;
    data = (double *)calloc(size_dat,sizeof(double));
    if(data == NULL)
      ath_error("[output_pdf]: Failed to allocate data array\n");

/* This choice for size_pdf represents a balance between
 * resolution in the PDF and "shot noise" in the data binning. */
#ifdef MPI_PARALLEL
    size_pdf = (int)sqrt((double)cg_data_cnt);
#else /* MPI_PARALLEL */
    size_pdf = (int)sqrt((double)size_dat);
#endif /* MPI_PARALLEL */
    pdf = (int *)calloc(size_pdf,sizeof(int));
    if(pdf == NULL)
      ath_error("[output_pdf]: Failed to allocate pdf array\n");

#ifdef MPI_PARALLEL
    if(pG->my_id == 0){ /* I'm the parent */
      cg_pdf = (int *)calloc(size_pdf,sizeof(int));
      if(cg_pdf == NULL)
	ath_error("[output_pdf]: Failed to allocate cg_pdf array\n");
    }
#endif /* MPI_PARALLEL */
  }


/* Initialize dmin, dmax */
  dmin = dmax = (*pout->expr)(pG,is,js,ks);

/* Fill the data array */
  n=0;
  for(k = ks; k<=ke; k++){
    for(j = js; j<=je; j++){
      for(i = is; i<=ie; i++){
	data[n] = (double)(*pout->expr)(pG,i,j,k);
	dmin = data[n] < dmin ? data[n] : dmin;
	dmax = data[n] > dmax ? data[n] : dmax;
	mean += data[n];
	n++;
      }
    }
  }

#ifdef MPI_PARALLEL

  dat = dmin;
  err = MPI_Allreduce(&dat, &dmin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if(err) ath_error("[dump_pdf]: MPI_Allreduce (dmin) error = %d\n",err);

  dat = dmax;
  err = MPI_Allreduce(&dat, &dmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(err) ath_error("[dump_pdf]: MPI_Allreduce (dmax) error = %d\n",err);

  dat = mean;
  err = MPI_Allreduce(&dat, &mean, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if(err) ath_error("[dump_pdf]: MPI_Allreduce (mean) error = %d\n",err);

  mean /= (double)cg_data_cnt; /* Complete the calc. of the mean */

#else /* MPI_PARALLEL */

  mean /= (double)data_cnt; /* Complete the calc. of the mean */

#endif /* MPI_PARALLEL */

  if(data_cnt > 1){
/* Calculate the variance, etc. with the corrected 2-pass formula */
    for(n=0; n<data_cnt; n++){
      s = data[n] - mean;
      adev += fabs(s);
      ep += s;
      var += (r = s*s);
      skew += (r *= s);
      kurt += (r *= s);
    }

#ifdef MPI_PARALLEL

    dat = ep;
    err = MPI_Allreduce(&dat, &ep, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(err) ath_error("[dump_pdf]: MPI_Allreduce (ep) error = %d\n",err);

    dat = var;
    err = MPI_Allreduce(&dat, &var, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(err) ath_error("[dump_pdf]: MPI_Allreduce (var) error = %d\n",err);

    dat = skew;
    err = MPI_Allreduce(&dat, &skew, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(err) ath_error("[dump_pdf]: MPI_Allreduce (skew) error = %d\n",err);

    dat = kurt;
    err = MPI_Allreduce(&dat, &kurt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(err) ath_error("[dump_pdf]: MPI_Allreduce (kurt) error = %d\n",err);

    adev /= (double)cg_data_cnt;
    var = (var - ep*ep/(double)cg_data_cnt)/(double)(cg_data_cnt-1);
    sdev = sqrt(var);
    if(sdev > 0.0){
      skew /= var*sdev*cg_data_cnt;
      kurt = kurt/(var*var*cg_data_cnt) - 3.0;
    }

#else /* MPI_PARALLEL */

    adev /= (double)data_cnt;
    var = (var - ep*ep/(double)data_cnt)/(double)(data_cnt-1);
    sdev = sqrt(var);
    if(sdev > 0.0){
      skew /= var*sdev*data_cnt;
      kurt = kurt/(var*var*data_cnt) - 3.0;
    }

#endif /* MPI_PARALLEL */
  }

/* Store the global maximum and minimum of the quantity */
  if(pout->num > 0){
    pout->gmin = dmin < pout->gmin ? dmin : pout->gmin;
    pout->gmax = dmax > pout->gmax ? dmax : pout->gmax;
  }
  else{
    pout->gmin = dmin;
    pout->gmax = dmax;
  }

/* Compute the pdf directly using sampling. Define size_pdf bins, each of equal
 * size, and fill them with the number of cells whose data value falls in the
 * range spanned by the bin. */

  if(dmax - dmin > 0.0){
/* Initialize pdf[] to zero */
    for(n=0; n<size_pdf; n++) pdf[n] = 0;
/* Calculate the number of cells whose data falls in each bin */
    scl = (double)size_pdf/(dmax - dmin);
    for(n=0; n<data_cnt; n++){
      i = (int)(scl*(data[n] - dmin));
      i = i < size_pdf ? i : size_pdf - 1;
      pdf[i]++;
    }
  }

#ifdef MPI_PARALLEL

/* Sum up the pdf in the array cg_pdf */
  err = MPI_Reduce(pdf, cg_pdf, size_pdf, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(err) ath_error("[dump_pdf]: MPI_Reduce (pdf) error = %d\n",err);

#endif /* MPI_PARALLEL */

#ifdef MPI_PARALLEL
/* For parallel calculations, only the parent writes the output. */
  if(pG->my_id != 0) return;
#endif /* MPI_PARALLEL */

/* Open the output file */
  pfile = ath_fopen(pG->outfilename,num_digit,pout->num,pout->id,"prb","w");
  if(pfile == NULL){
    fprintf(stderr,"[output_pdf]: File Open Error Occured");
    return;
  }

/* Write out some extra information in a header */
  fprintf(pfile,"# Time = %21.15e\n",pG->time);
  fprintf(pfile,"# nstep = %d\n",pG->nstep);
  fprintf(pfile,"# expr = \"%s\"\n",pout->out);
  fprintf(pfile,"# Nbin = %d\n",((dmax - dmin) > 0.0 ? size_pdf : 1));
  fprintf(pfile,"# dmin = %21.15e\n",dmin); 
  fprintf(pfile,"# dmax = %21.15e\n",dmax); 
  fprintf(pfile,"# mean = %21.15e\n",mean); 
  fprintf(pfile,"# variance = %21.15e\n",var); 
  fprintf(pfile,"# std. dev. = %21.15e\n",sdev); 
  fprintf(pfile,"# avg. dev. = %21.15e\n",adev); 
  fprintf(pfile,"# skewness = %21.15e\n",skew); 
  fprintf(pfile,"# kurtosis = %21.15e\n#\n",kurt); 

/* Add a white space to the format */
  if(pout->dat_fmt == NULL)
    sprintf(fmt,"%s  %s\n",def_fmt, def_fmt);
  else
    sprintf(fmt,"%s  %s\n",pout->dat_fmt,pout->dat_fmt);

/* write out the normalized Proabability Distribution Function */
  if(dmax - dmin > 0.0){
    delta = (dmax - dmin)/(double)(size_pdf);
#ifdef MPI_PARALLEL
    scl = (double)size_pdf/(double)(cg_data_cnt*(dmax - dmin));
#else
    scl = (double)size_pdf/(double)(data_cnt*(dmax - dmin));
#endif /* MPI_PARALLEL */
    for(n=0; n<size_pdf; n++){
/* Calculate the normalized Prob. Dist. Fun. */
      dat = dmin + (n + 0.5)*delta;
#ifdef MPI_PARALLEL
      dpdf = (double)(cg_pdf[n])*scl;
#else
      dpdf = (double)(pdf[n])*scl;
#endif /* MPI_PARALLEL */
      fprintf(pfile, fmt, dat, dpdf);
    }
  }
  else
    fprintf(pfile,fmt,dmax,1.0);

  fclose(pfile);

/* Also write a history type file on the statistics */
  sprintf(fid,"prb_stat.%s",pout->id);
  pfile = ath_fopen(pG->outfilename,0,0,fid,"tab","a");
  if(pfile == NULL){
    fprintf(stderr,"[output_pdf]: File Open Error Occured");
    return;
  }

  if(pout->num == 0){
    fprintf(pfile,"# expr = \"%s\"\n#\n",pout->out);
    fprintf(pfile,"# time  nstep  dmin  dmax  mean  variance  \"std. dev.\"  ");
    fprintf(pfile,"\"avg. dev.\"  skewness  kurtosis\n#\n");
  }

/* Add a white space to the format */
  if(pout->dat_fmt == NULL) sprintf(fmt," %s",def_fmt);
  else                      sprintf(fmt," %s",pout->dat_fmt);

  fprintf(pfile,"%21.15e %d",pG->time,pG->nstep);

/* write out the table of statistics */
  fprintf(pfile,fmt,dmin);
  fprintf(pfile,fmt,dmax);
  fprintf(pfile,fmt,mean);
  fprintf(pfile,fmt,var);
  fprintf(pfile,fmt,sdev);
  fprintf(pfile,fmt,adev);
  fprintf(pfile,fmt,skew);
  fprintf(pfile,fmt,kurt);
  fprintf(pfile,"\n");

  fclose(pfile);

  return;
}
