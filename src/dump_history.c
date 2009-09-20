#include "copyright.h"
/*==============================================================================
 * FILE: dump_history.c
 *
 * PURPOSE: Functions to write dumps of scalar "history" variables in a
 *   formatted table.  "History" dumps are scalars (usually volume averages)
 *   written periodically, giving the time-evolution of that quantitity.
 *   Default (hardwired) values are:
 *     scal[0] = time
 *     scal[1] = dt
 *     scal[2] = mass
 *     scal[3] = total energy
 *     scal[4] = d*v1
 *     scal[5] = d*v2
 *     scal[6] = d*v3
 *     scal[7] = 0.5*d*v1**2
 *     scal[8] = 0.5*d*v2**2
 *     scal[9] = 0.5*d*v3**2
 *     scal[10] = 0.5*b1**2
 *     scal[11] = 0.5*b2**2
 *     scal[12] = 0.5*b3**2
 *     scal[13] = d*Phi
 *     scal[14+NSCALARS] = passively advected scalars
 *     scal[15+NSCALARS] = angular momentum
 * More variables can be hardwired by increasing NSCAL=number of variables, and
 * adding calculation of desired quantities below.
 *
 * Alternatively, up to MAX_USR_H_COUNT new history variables can be added using
 * dump_history_enroll() in the problem generator.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   dump_history()        - Writes variables as formatted table
 *   dump_history_enroll() - Adds new user-defined history variables
 *============================================================================*/

#include <stdio.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/* Maximum Number of default history dump columns. */
#define NSCAL 14

/* Maximum number of history dump columns that the user routine can add. */
#define MAX_USR_H_COUNT 10

/* Array of strings / labels for the user added history column. */
static char *usr_label[MAX_USR_H_COUNT];

/* Array of history dump function pointers for user added history columns. */
static Gasfun_t phst_fun[MAX_USR_H_COUNT];

static int usr_hst_cnt = 0; /* User History Counter <= MAX_USR_H_COUNT */

/*----------------------------------------------------------------------------*/
/* dump_history:  */

void dump_history(Grid *pGrid, Domain *pD, Output *pOut)
{
  int i, is = pGrid->is, ie = pGrid->ie;
  int j, js = pGrid->js, je = pGrid->je;
  int k, ks = pGrid->ks, ke = pGrid->ke;
  double dvol, scal[NSCAL + NSCALARS + MAX_USR_H_COUNT], d1;
  FILE *p_hstfile;
  char fmt[80];
  int vol_rat; /* (Grid Volume)/(dx1*dx2*dx3) */
  int n, total_hst_cnt, mhst;
#ifdef MPI_PARALLEL
  double my_scal[NSCAL + NSCALARS + MAX_USR_H_COUNT]; /* My Volume averages */
  int my_vol_rat; /* My (Grid Volume)/(dx1*dx2*dx3) */
  int err;
#endif

  Real x1=1.0;
#ifdef CYLINDRICAL
  Real x2, x3;
  double Rmin, Rmax;
#endif 

  total_hst_cnt = 9 + NSCALARS + usr_hst_cnt;
#ifdef ADIABATIC
  total_hst_cnt++;
#endif
#ifdef MHD
  total_hst_cnt += 3;
#endif
#ifdef SELF_GRAVITY
  total_hst_cnt += 1;
#endif
#ifdef CYLINDRICAL
  total_hst_cnt++;  /* for angular momentum */
#endif

/* Add a white space to the format */
  if(pOut->dat_fmt == NULL){
    sprintf(fmt," %%14.6e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

  for (i=2; i<total_hst_cnt; i++) {
    scal[i] = 0.0;
  }
 
/* Compute history variables */

  scal[0] = pGrid->time;
  scal[1] = pGrid->dt;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        /* WEIGHT THE VARIABLES IN EACH CELL BY x1 */
#ifdef CYLINDRICAL
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
#endif

        mhst = 2;
	scal[mhst] += x1*pGrid->U[k][j][i].d;
	d1 = 1.0/pGrid->U[k][j][i].d;
#ifndef BAROTROPIC
        mhst++;
	scal[mhst] += x1*pGrid->U[k][j][i].E;
#endif
        mhst++;
	scal[mhst] += x1*pGrid->U[k][j][i].M1;
        mhst++;
	scal[mhst] += x1*pGrid->U[k][j][i].M2;
        mhst++;
	scal[mhst] += x1*pGrid->U[k][j][i].M3;
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].M1)*d1;
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].M2)*d1;
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].M3)*d1;
#ifdef MHD
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].B1c);
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].B2c);
        mhst++;
	scal[mhst] += x1*0.5*SQR(pGrid->U[k][j][i].B3c);
#endif
#ifdef SELF_GRAVITY
        mhst++;
        scal[mhst] += x1*pGrid->U[k][j][i].d*pGrid->Phi[k][j][i];
#endif
#if (NSCALARS > 0)
	for(n=0; n<NSCALARS; n++){
          mhst++;
	  scal[mhst] += x1*pGrid->U[k][j][i].s[n];
	}
#endif

#ifdef CYLINDRICAL
        mhst++;
        scal[mhst] += x1*(x1*pGrid->U[k][j][i].M2);
#endif

/* Calculate the user defined history variables */
	for(n=0; n<usr_hst_cnt; n++){
          mhst++;
	  scal[mhst] += x1*(*phst_fun[n])(pGrid, i, j, k);
	}
      }
    }
  }

/* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  vol_rat = (pD->ide - pD->ids + 1)*(pD->jde - pD->jds + 1)*
            (pD->kde - pD->kds + 1);
#ifdef CYLINDRICAL
  cc_pos(pGrid,pD->ids,pD->jds,pD->kds,&Rmin,&x2,&x3);
  cc_pos(pGrid,pD->ide,pD->jde,pD->kde,&Rmax,&x2,&x3);
  Rmin -= 0.5*pGrid->dx1;
  Rmax += 0.5*pGrid->dx1;
#endif

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i=2; i<total_hst_cnt; i++){
    my_scal[i] = scal[i];
  }
  err = MPI_Reduce(&(my_scal[2]), &(scal[2]), total_hst_cnt - 2,
		   MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(err)
    ath_error("[dump_history]: MPI_Reduce call returned error = %d\n",err);
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(pGrid->my_id == 0){ /* I'm the parent */
#endif

  dvol = 1.0/(double)vol_rat;
#ifdef CYLINDRICAL
  dvol /= 0.5*(Rmin+Rmax);
#endif

  for (i=2; i<total_hst_cnt; i++) {
    scal[i] *= dvol;
  }

/* Write history file */
#ifdef MPI_PARALLEL
  p_hstfile = ath_fopen("../",pGrid->outfilename,0,0,NULL,"hst","a");
#else
  p_hstfile = ath_fopen(NULL,pGrid->outfilename,0,0,NULL,"hst","a");
#endif
  if(p_hstfile == NULL){
    ath_perr(-1,"[dump_history]: Unable to open the history file\n");
    return;
  }

/* Write out column headers */

  mhst = 0;
  if(pOut->num == 0){
    mhst++;
    fprintf(p_hstfile,"#   [%i]=time   ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=dt      ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=mass    ",mhst);
#ifdef ADIABATIC
    mhst++;
    fprintf(p_hstfile,"   [%i]=total E ",mhst);
#endif
    mhst++;
    fprintf(p_hstfile,"   [%i]=x1 Mom. ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x2 Mom. ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x3 Mom. ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x1-KE   ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x2-KE   ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x3-KE   ",mhst);
#ifdef MHD
    mhst++;
    fprintf(p_hstfile,"   [%i]=x1-ME   ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x2-ME   ",mhst);
    mhst++;
    fprintf(p_hstfile,"   [%i]=x3-ME   ",mhst);
#endif
#ifdef SELF_GRAVITY
    mhst++;
    fprintf(p_hstfile,"   [%i]=grav PE ",mhst);
#endif
#if (NSCALARS > 0)
    for(n=0; n<NSCALARS; n++){
      mhst++;
      fprintf(p_hstfile,"  [%i]=scalar %i",mhst,n);
    }
#endif

#ifdef CYLINDRICAL
    mhst++;
    fprintf(p_hstfile,"   [%i]=Ang.Mom.",mhst);
#endif

    for(n=0; n<usr_hst_cnt; n++){
      mhst++;
      fprintf(p_hstfile,"  [%i]=%s",mhst,usr_label[n]);
    }
    fprintf(p_hstfile,"\n#\n");
  }

/* Write out data */

  for (i=0; i<total_hst_cnt; i++) {
    fprintf(p_hstfile,fmt,scal[i]);
  }
  fprintf(p_hstfile,"\n");
  fclose(p_hstfile);

#ifdef MPI_PARALLEL
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* dump_history_enroll:  */

void dump_history_enroll(const Gasfun_t pfun, const char *label){

  if(usr_hst_cnt >= MAX_USR_H_COUNT)
    ath_error("[dump_history_enroll]: MAX_USR_H_COUNT = %d exceeded\n",
	      MAX_USR_H_COUNT);

/* Copy the label string */
  if((usr_label[usr_hst_cnt] = ath_strdup(label)) == NULL)
    ath_error("[dump_history_enroll]: Error on sim_strdup(\"%s\")\n",label);

/* Store the function pointer */
  phst_fun[usr_hst_cnt] = pfun;

  usr_hst_cnt++;

  return;

}

#undef NSCAL
#undef MAX_USR_H_COUNT
