#include "copyright.h"
/*==============================================================================
 * FILE: dump_history.c
 *
 * PURPOSE: Functions to write dumps of scalar "history" variables in a
 *   formatted table.  "History" variables are usually volume integrals: for
 *   example S = \Sum_{ijk} q[k][j][i]*dx1[i]*dx2[j]*dx3[k], although other
 *   quantities like time and dt are also output.  Note history variables are
 *   TOTAL, VOLUME INTEGRATED quantities, NOT volume averages (just divide by
 *   the Domain volume to compute the latter).  History data is written
 *   periodically, giving the time-evolution of that quantitity.
 * Default (hardwired) values are:
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
 * Data is integrated over a single Domain, by default the root (level=0)
 * Domain.  If <output>level or <output>domain != 0, then the data integrated
 * over this Domain will be output to a seperate file with <output>id added to
 * the filename (this can be set by user to identify file).
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
#define MAX_USR_H_COUNT 30

/* Array of strings / labels for the user added history column. */
static char *usr_label[MAX_USR_H_COUNT];

/* Array of history dump function pointers for user added history columns. */
static GasFun_t phst_fun[MAX_USR_H_COUNT];

static int usr_hst_cnt = 0; /* User History Counter <= MAX_USR_H_COUNT */

/*----------------------------------------------------------------------------*/
/* dump_history:  */

void dump_history(MeshS *pM, OutputS *pOut)
{
  GridS *pG;
  DomainS *pD;
  int i,j,k,is,ie,js,je,ks,ke;
  double dVol, scal[NSCAL + NSCALARS + MAX_USR_H_COUNT], d1;
  FILE *pfile;
  char fmt[80];
  int n, total_hst_cnt, mhst;
#ifdef MPI_PARALLEL
  double my_scal[NSCAL + NSCALARS + MAX_USR_H_COUNT]; /* My Volume averages */
  int ierr, myID_Comm_Domain;
#endif
#ifdef CYLINDRICAL
  Real x1,x2,x3;
#endif 

/* Return if Grid is not on this processor */

  pG = pM->Domain[pOut->nlevel][pOut->ndomain].Grid;
  if (pG == NULL) return;

#ifdef MPI_PARALLEL
  pD = (DomainS*)&(pM->Domain[pOut->nlevel][pOut->ndomain]);
#endif
  is = pG->is, ie = pG->ie;
  js = pG->js, je = pG->je;
  ks = pG->ks, ke = pG->ke;

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

  scal[0] = pG->time;
  scal[1] = pG->dt;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        dVol = x1*(pG->dx1)*(pG0->dx2)*(pG->dx3); 
#else
        dVol = (pG->dx1)*(pG->dx2)*(pG->dx3); 
#endif

        mhst = 2;
	scal[mhst] += dVol*pG->U[k][j][i].d;
	d1 = 1.0/pG->U[k][j][i].d;
#ifndef BAROTROPIC
        mhst++;
	scal[mhst] += dVol*pG->U[k][j][i].E;
#endif
        mhst++;
	scal[mhst] += dVol*pG->U[k][j][i].M1;
        mhst++;
	scal[mhst] += dVol*pG->U[k][j][i].M2;
        mhst++;
	scal[mhst] += dVol*pG->U[k][j][i].M3;
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M1)*d1;
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M2)*d1;
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].M3)*d1;
#ifdef MHD
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B1c);
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B2c);
        mhst++;
	scal[mhst] += dVol*0.5*SQR(pG->U[k][j][i].B3c);
#endif
#ifdef SELF_GRAVITY
        mhst++;
        scal[mhst] += dVol*pG->U[k][j][i].d*pG->Phi[k][j][i];
#endif
#if (NSCALARS > 0)
	for(n=0; n<NSCALARS; n++){
          mhst++;
	  scal[mhst] += dVol*pG->U[k][j][i].s[n];
	}
#endif

#ifdef CYLINDRICAL
        mhst++;
        scal[mhst] += dVol*(x1*pG->U[k][j][i].M2);
#endif

/* Calculate the user defined history variables */
	for(n=0; n<usr_hst_cnt; n++){
          mhst++;
	  scal[mhst] += dVol*(*phst_fun[n])(pG, i, j, k);
	}
      }
    }
  }

/* Compute the sum over all Grids in Domain */

#ifdef MPI_PARALLEL 
  for(i=2; i<total_hst_cnt; i++){
    my_scal[i] = scal[i];
  }
  ierr = MPI_Reduce(&(my_scal[2]), &(scal[2]), (total_hst_cnt - 2),
    MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
#endif

/* Only the parent (rank=0) process computes the average and writes output. */
#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

/* Write history file */
#ifdef MPI_PARALLEL
  if (pOut->nlevel > 0) {
    pfile = ath_fopen("../",pM->outfilename,0,0,pOut->id,"hst","a");
  } else {
    pfile = ath_fopen("../",pM->outfilename,0,0,NULL,"hst","a");
  }
#else
  if (pOut->nlevel > 0) {
    pfile = ath_fopen(NULL,pM->outfilename,0,0,pOut->id,"hst","a");
  } else {
    pfile = ath_fopen(NULL,pM->outfilename,0,0,NULL,"hst","a");
  }
#endif
  if(pfile == NULL){
    ath_perr(-1,"[dump_history]: Unable to open the history file\n");
    return;
  }

/* Write out column headers */

  mhst = 0;
  if(pOut->num == 0){
    mhst++;
    fprintf(pfile,"#   [%i]=time   ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=dt      ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=mass    ",mhst);
#ifdef ADIABATIC
    mhst++;
    fprintf(pfile,"   [%i]=total E ",mhst);
#endif
    mhst++;
    fprintf(pfile,"   [%i]=x1 Mom. ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x2 Mom. ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x3 Mom. ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x1-KE   ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x2-KE   ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x3-KE   ",mhst);
#ifdef MHD
    mhst++;
    fprintf(pfile,"   [%i]=x1-ME   ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x2-ME   ",mhst);
    mhst++;
    fprintf(pfile,"   [%i]=x3-ME   ",mhst);
#endif
#ifdef SELF_GRAVITY
    mhst++;
    fprintf(pfile,"   [%i]=grav PE ",mhst);
#endif
#if (NSCALARS > 0)
    for(n=0; n<NSCALARS; n++){
      mhst++;
      fprintf(pfile,"  [%i]=scalar %i",mhst,n);
    }
#endif

#ifdef CYLINDRICAL
    mhst++;
    fprintf(pfile,"   [%i]=Ang.Mom.",mhst);
#endif

    for(n=0; n<usr_hst_cnt; n++){
      mhst++;
      fprintf(pfile,"  [%i]=%s",mhst,usr_label[n]);
    }
    fprintf(pfile,"\n#\n");
  }

/* Write out data */

  for (i=0; i<total_hst_cnt; i++) {
    fprintf(pfile,fmt,scal[i]);
  }
  fprintf(pfile,"\n");
  fclose(pfile);

#ifdef MPI_PARALLEL
  }
#endif

  return;
}

/*----------------------------------------------------------------------------*/
/* dump_history_enroll:  */

void dump_history_enroll(const GasFun_t pfun, const char *label){

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
