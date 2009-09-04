#include "copyright.h"
/*==============================================================================
 * FILE: restart.c
 *
 * PURPOSE: Functions for writing and reading restart files.  Restart files
 *   begin with the entire input parameter (athinput.XX) file used to start the
 *   run in text format, with nstep, time, and dt appended, followed by the Gas
 *   array and face-centered B in binary, followed by any problem-specific data
 *   written by a user function in the problem generator.  The restart input
 *   file is parsed by main() using the standard par_* functions, with
 *   the parameters used by init_grid_block() to initialize the grid; values
 *   can be superceded by input from the command line, or another input file.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   restart_grid_block() - reads nstep,time,dt,Gas and B from restart file 
 *   dump_restart()       - writes a restart file
 *
 * VARIABLE TYPE AND STRUCTURE DEFINITIONS: none
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "particles/particle.h"

/*----------------------------------------------------------------------------*/
/* restart_grid_block: Reads nstep,time,dt,arrays of Gas and face-centered B in 
 *   Grid structure from restart file.  The remaining data is initialized
 *   in init_grid_block using input parameters at start of restart file,
 *   the command line, or from a new input file.
 */

void restart_grid_block(char *res_file, Grid *pG, Domain *pD)
{
  FILE *fp;
  char line[MAXLEN];
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
#ifdef MHD
  int ib, jb, kb;
#endif
#if (NSCALARS > 0)
  int n;
  char scalarstr[16];
#endif
#ifdef ION_RADPOINT
  Real x1rad, x2rad, x3rad, srad;
  int nradpoint, rebuildctr;
  Rad_Ran2_State ranstate;
  float rotation[3][3];
#endif
#ifdef ION_RADPLANE
  int dir, nradplane;
  Real flux;
#endif
#ifdef PARTICLES
  long p;
#endif

/* Open the restart file */
  if((fp = fopen(res_file,"r")) == NULL)
    ath_error("[restart_grid_block]: Error opening the restart file\nIf this is a MPI job, make sure each file from each processor is in the same directory.\n");

/* Skip over the parameter file at the start of the restart file */
  do{
    fgets(line,MAXLEN,fp);
  }while(strncmp(line,"<par_end>",9) != 0);

/* read nstep */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"N_STEP",6) != 0)
    ath_error("[restart_grid_block]: Expected N_STEP, found %s",line);
  fread(&(pG->nstep),sizeof(int),1,fp);

/* read time */
  fgets(line,MAXLEN,fp);   /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME",4) != 0)
    ath_error("[restart_grid_block]: Expected TIME, found %s",line);
  fread(&(pG->time),sizeof(Real),1,fp);

/* read dt */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME_STEP",9) != 0)
    ath_error("[restart_grid_block]: Expected TIME_STEP, found %s",line);
  fread(&(pG->dt),sizeof(Real),1,fp);

/* Read the density */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"DENSITY",7) != 0)
    ath_error("[restart_grid_block]: Expected DENSITY, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->U[k][j][i].d),sizeof(Real),1,fp);
      }
    }
  }

/* Read the x1-momentum */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"1-MOMENTUM",10) != 0)
    ath_error("[restart_grid_block]: Expected 1-MOMENTUM, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->U[k][j][i].M1),sizeof(Real),1,fp);
      }
    }
  }

/* Read the x2-momentum */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"2-MOMENTUM",10) != 0)
    ath_error("[restart_grid_block]: Expected 2-MOMENTUM, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->U[k][j][i].M2),sizeof(Real),1,fp);
      }
    }
  }

/* Read the x3-momentum */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"3-MOMENTUM",10) != 0)
    ath_error("[restart_grid_block]: Expected 3-MOMENTUM, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->U[k][j][i].M3),sizeof(Real),1,fp);
      }
    }
  }

#ifndef BAROTROPIC
/* Read energy density */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"ENERGY",6) != 0)
    ath_error("[restart_grid_block]: Expected ENERGY, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->U[k][j][i].E),sizeof(Real),1,fp);
      }
    }
  }
#endif

#ifdef MHD
/* if there is more than one cell in each dimension, need to read one more
 * face-centered field component than the number of cells.  [ijbk]b is
 * the number of extra cells to be read  */

  ib = ie > is ? 1 : 0;
  jb = je > js ? 1 : 0;
  kb = ke > ks ? 1 : 0;

/* Read the x1-field */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"1-FIELD",7) != 0)
    ath_error("[restart_grid_block]: Expected 1-FIELD, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+ib; i++) {
	fread(&(pG->B1i[k][j][i]),sizeof(Real),1,fp);
      }
    }
  }

/* Read the x2-field */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"2-FIELD",7) != 0)
    ath_error("[restart_grid_block]: Expected 2-FIELD, found %s",line);
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+jb; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->B2i[k][j][i]),sizeof(Real),1,fp);
      }
    }
  }

/* Read the x3-field */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"3-FIELD",7) != 0)
    ath_error("[restart_grid_block]: Expected 3-FIELD, found %s",line);
  for (k=ks; k<=ke+kb; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	fread(&(pG->B3i[k][j][i]),sizeof(Real),1,fp);
      }
    }
  }

/* initialize the cell center magnetic fields as either the average of the face
 * centered field if there is more than one cell in that dimension, or just
 * the face centered field if not  */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	pG->U[k][j][i].B1c = 
	  ib ? 0.5*(pG->B1i[k][j][i] + pG->B1i[k][j][i+1]) : pG->B1i[k][j][i];

	pG->U[k][j][i].B2c = 
	  jb ? 0.5*(pG->B2i[k][j][i] + pG->B2i[k][j+1][i]) : pG->B2i[k][j][i];

	pG->U[k][j][i].B3c = 
	  kb ? 0.5*(pG->B3i[k][j][i] + pG->B3i[k+1][j][i]) : pG->B3i[k][j][i];
      }
    }
  }
#endif

#ifdef ION_RADPOINT
  /* Read list of radiators and call routine to add them */
  fgets(line,MAXLEN,fp);/* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"RADIATOR LIST",13) != 0)
    ath_error("[restart_grid_block]: Expected RADIATOR LIST, found %s",line);
  fread(&nradpoint,sizeof(int),1,fp);
  for (n=0; n<nradpoint; n++) {
    fread(&x1rad,sizeof(Real),1,fp);
    fread(&x2rad,sizeof(Real),1,fp);
    fread(&x3rad,sizeof(Real),1,fp);
    fread(&srad,sizeof(Real),1,fp);
    fread(&rebuildctr,sizeof(int),1,fp);
    fread(rotation,sizeof(float),9,fp);
    restore_radpoint_3d(pG, x1rad, x2rad, x3rad, srad, rebuildctr, rotation);
  }

  /* Read the current random number generator state and set it */
  fgets(line,MAXLEN,fp);/* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"RAN2 STATE",10) != 0)
    ath_error("[restart_grid_block]: Expected RAN2 STATE, found %s",line);
  fread(&ranstate,sizeof(Rad_Ran2_State),1,fp);
  ion_radpoint_set_ranstate(ranstate);
#endif /* ION_RADPOINT */

#ifdef ION_RADPLANE
  /* Read list of radiator planes */
  fgets(line,MAXLEN,fp);/* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"RADIATOR PLANE LIST",19) != 0)
    ath_error("[restart_grid_block]: Expected RADIATOR PLANE LIST, found %s",line);
  fread(&nradplane,sizeof(int),1,fp);
  for (n=0; n<nradplane; n++) {
    fread(&dir,sizeof(int),1,fp);
    fread(&flux,sizeof(Real),1,fp);
    add_radplane_3d(pG, dir, flux);
  }
#endif /* ION_RADPLANE */

/* Read any passively advected scalars */
/* Following code only works if NSCALARS < 10 */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    fgets(line,MAXLEN,fp);/* Read the '\n' preceeding the next string */
    fgets(line,MAXLEN,fp);
    sprintf(scalarstr, "SCALAR %d", n);
    if(strncmp(line,scalarstr,8) != 0)
      ath_error("[restart_grid_block]: Expected %s, found %s",scalarstr,line);
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          fread(&(pG->U[k][j][i].s[n]),sizeof(Real),1,fp);
        }
      }
    }
  }
#endif

/* Read particle properties and the complete particle list */
#ifdef PARTICLES
  fgets(line,MAXLEN,fp);/* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE LIST",13) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE LIST, found %s",line);
  fread(&(pG->nparticle),sizeof(long),1,fp);
  fread(&(pG->partypes),sizeof(int),1,fp);
  for (i=0; i<pG->partypes; i++) {          /* particle property list */
#ifdef FEEDBACK
    fread(&(pG->grproperty[i].m),sizeof(Real),1,fp);
#endif
    fread(&(pG->grproperty[i].rad),sizeof(Real),1,fp);
    fread(&(pG->grproperty[i].rho),sizeof(Real),1,fp);
    fread(&(tstop0[i]),sizeof(Real),1,fp);
    fread(&(grrhoa[i]),sizeof(Real),1,fp);
  }
  fread(&(alamcoeff),sizeof(Real),1,fp);   /* coefficient to calculate the Reynolds number */

  for (i=0; i<pG->partypes; i++)
    fread(&(pG->grproperty[i].integrator),sizeof(short),1,fp);

/* Read the x1-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE X1",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE X1, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].x1),sizeof(Real),1,fp);
  }

/* Read the x2-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE X2",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE X2, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].x2),sizeof(Real),1,fp);
  }

/* Read the x3-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE X3",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE X3, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].x3),sizeof(Real),1,fp);
  }

/* Read the v1-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE V1",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE V1, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].v1),sizeof(Real),1,fp);
  }

/* Read the v2-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE V2",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE V2, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].v2),sizeof(Real),1,fp);
  }

/* Read the v3-coordinate */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE V3",11) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE V3, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].v3),sizeof(Real),1,fp);
  }

/* Read particle properties */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE PROPERTY",17) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE PROPERTY, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].property),sizeof(int),1,fp);
    pG->particle[p].pos = 1;	/* grid particle */
  }

/* Read particle my_id */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE MY_ID",14) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE MY_ID, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].my_id),sizeof(long),1,fp);
  }

#ifdef MPI_PARALLEL
/* Read particle init_id */
  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"PARTICLE INIT_ID",16) != 0)
    ath_error("[restart_particle_block]: Expected PARTICLE INIT_ID, found %s",line);
  for (p=0; p<pG->nparticle; p++) {
    fread(&(pG->particle[p].init_id),sizeof(int),1,fp);
  }
#endif

/* count the number of particles with different types */
  for (i=0; i<pG->partypes; i++)
    pG->grproperty[i].num = 0;
  for (p=0; p<pG->nparticle; p++)
    pG->grproperty[pG->particle[p].property].num += 1;

#endif /* PARTICLES */

  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"USER_DATA",9) != 0)
    ath_error("[restart_grid_block]: Expected USER_DATA, found %s",line);

/* Call a user function to read his/her problem-specific data! */
  problem_read_restart(pG, pD, fp);

  fclose(fp);

  return;
}

/*----------------------------------------------------------------------------*/
/* dump_restart: writes a restart file, including problem-specific data from
 *   a user defined function  */

void dump_restart(Grid *pG, Domain *pD, Output *pout)
{
  FILE *fp;
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
#ifdef MHD
  int ib, jb, kb;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef ION_RADPOINT
  Rad_Ran2_State ranstate;
#endif
#ifdef PARTICLES
  int nprop, *ibuf = NULL, ibufsize, lbufsize, sbufsize, nibuf, nlbuf, nsbuf;
  long np, p, *lbuf = NULL;
  short *sbuf = NULL;
  nibuf = 0;	nsbuf = 0;	nlbuf = 0;
#endif
  int bufsize, nbuf = 0;
  Real *buf = NULL;

/* Allocate memory for buffer */
  bufsize = 262144 / sizeof(Real);  /* 256 KB worth of Reals */
  if ((buf = (Real*)calloc_1d_array(bufsize, sizeof(Real))) == NULL) {
    ath_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
#ifdef PARTICLES
  ibufsize = 262144 / sizeof(int);  /* 256 KB worth of ints */
  if ((ibuf = (int*)calloc_1d_array(ibufsize, sizeof(int))) == NULL) {
    ath_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
  lbufsize = 262144 / sizeof(long);  /* 256 KB worth of longs */
  if ((lbuf = (long*)calloc_1d_array(lbufsize, sizeof(long))) == NULL) {
    ath_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
  sbufsize = 262144 / sizeof(short);  /* 256 KB worth of longs */
  if ((sbuf = (short*)calloc_1d_array(sbufsize, sizeof(short))) == NULL) {
    ath_perr(-1,"[dump_restart]: Error allocating memory for buffer\n");
    return;  /* Right now, we just don't write instead of aborting completely */
  }
#endif

/* Open the output file */
  fp = ath_fopen(NULL,pG->outfilename,num_digit,pout->num,NULL,"rst","wb");
  if(fp == NULL){
    ath_perr(-1,"[dump_restart]: Error opening the restart file\n");
    return;
  }

/* Add the current time & nstep to the parameter file */
  par_setd("time","time","%e",pG->time,"Current Simulation Time");
  par_seti("time","nstep","%d",pG->nstep,"Current Simulation Time Step");

/* Write the current state of the parameter file */
  par_dump(2,fp);

/* Write out the current simulation step number */
  fprintf(fp,"N_STEP\n");
  if(fwrite(&(pG->nstep),sizeof(int),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time */
  fprintf(fp,"\nTIME\n");
  if(fwrite(&(pG->time),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time step */
  fprintf(fp,"\nTIME_STEP\n");
  if(fwrite(&(pG->dt),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");

/* Write the density */
  fprintf(fp,"\nDENSITY\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->U[k][j][i].d;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

/* Write the x1-momentum */
  fprintf(fp,"\n1-MOMENTUM\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->U[k][j][i].M1;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

/* Write the x2-momentum */
  fprintf(fp,"\n2-MOMENTUM\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->U[k][j][i].M2;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

/* Write the x3-momentum */
  fprintf(fp,"\n3-MOMENTUM\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->U[k][j][i].M3;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

#ifndef BAROTROPIC
/* Write energy density */
  fprintf(fp,"\nENERGY\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->U[k][j][i].E;
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }
#endif

#ifdef MHD
/* see comments in restart_grid_block() for use of [ijk]b */
  ib = ie > is ? 1 : 0;
  jb = je > js ? 1 : 0;
  kb = ke > ks ? 1 : 0;

/* Write the x1-field */
  fprintf(fp,"\n1-FIELD\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie+ib; i++) {
        buf[nbuf++] = pG->B1i[k][j][i];
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

/* Write the x2-field */
  fprintf(fp,"\n2-FIELD\n");
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je+jb; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->B2i[k][j][i];
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

/* Write the x3-field */
  fprintf(fp,"\n3-FIELD\n");
  for (k=ks; k<=ke+kb; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        buf[nbuf++] = pG->B3i[k][j][i];
        if ((nbuf+1) > bufsize) {
          fwrite(buf,sizeof(Real),nbuf,fp);
          nbuf = 0;
        }
      }
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }
#endif

#ifdef ION_RADPOINT
  /* Write out properties of radiators */
  fprintf(fp,"\nRADIATOR LIST\n");
  fwrite(&(pG->nradpoint),sizeof(int),1,fp);
  for (n=0; n<pG->nradpoint; n++) {
    fwrite(&(pG->radpointlist[n].x1),sizeof(Real),1,fp);
    fwrite(&(pG->radpointlist[n].x2),sizeof(Real),1,fp);
    fwrite(&(pG->radpointlist[n].x3),sizeof(Real),1,fp);
    fwrite(&(pG->radpointlist[n].s),sizeof(Real),1,fp);
    fwrite(&(pG->radpointlist[n].tree.rebuild_ctr),sizeof(int),1,fp);
    fwrite(&(pG->radpointlist[n].tree.rotation),sizeof(float),9,fp);
  }

  /* Write out the state of the random number generator */
  fprintf(fp,"\nRAN2 STATE\n");
  ranstate=ion_radpoint_get_ranstate();
  fwrite(&ranstate, sizeof(Rad_Ran2_State), 1, fp);

#endif /* ION_RADPOINT */

#ifdef ION_RADPLANE
  /* Write out properties of radiator planes */
  fprintf(fp,"\nRADIATOR PLANE LIST\n");
  fwrite(&(pG->nradplane),sizeof(int),1,fp);
  for (n=0; n<pG->nradplane; n++) {
    fwrite(&(pG->radplanelist[n].dir),sizeof(int),1,fp);
    fwrite(&(pG->radplanelist[n].flux),sizeof(Real),1,fp);
  }
#endif /* ION_RADPLANE */

/* Write out passively advected scalars */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    fprintf(fp,"\nSCALAR %d\n", n);
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          buf[nbuf++] = pG->U[k][j][i].s[n];
          if ((nbuf+1) > bufsize) {
            fwrite(buf,sizeof(Real),nbuf,fp);
            nbuf = 0;
          }
        }
      }
    }
    if (nbuf > 0) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
#endif

#ifdef PARTICLES
  /* Write out the number of particles */
  fprintf(fp,"\nPARTICLE LIST\n");
  np = 0;
  for (p=0; p<pG->nparticle; p++)
    if (pG->particle[p].pos == 1) np += 1;
  fwrite(&(np),sizeof(long),1,fp);

  /* Write out the particle properties */
#ifdef FEEDBACK
  nprop = 5;
#else
  nprop = 4;
#endif
  fwrite(&(pG->partypes),sizeof(int),1,fp); /* number of particle types */
  for (i=0; i<pG->partypes; i++) {          /* particle property list */
#ifdef FEEDBACK
    buf[nbuf++] = pG->grproperty[i].m;
#endif
    buf[nbuf++] = pG->grproperty[i].rad;
    buf[nbuf++] = pG->grproperty[i].rho;
    buf[nbuf++] = tstop0[i];
    buf[nbuf++] = grrhoa[i];
    if ((nbuf+nprop) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }
  fwrite(&(alamcoeff),sizeof(Real),1,fp);  /* coefficient to calculate the Reynolds number */

  for (i=0; i<pG->partypes; i++) {          /* particle integrator type */
    sbuf[nsbuf++] = pG->grproperty[i].integrator;
    if ((nsbuf+1) > sbufsize) {
      fwrite(sbuf,sizeof(short),nsbuf,fp);
      nsbuf = 0;
    }
  }
  if (nsbuf > 0) {
    fwrite(sbuf,sizeof(short),nsbuf,fp);
    nsbuf = 0;
  }

  /* Write x1 */
  fprintf(fp,"\nPARTICLE X1\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].x1;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write x2 */
  fprintf(fp,"\nPARTICLE X2\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].x2;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write x3 */
  fprintf(fp,"\nPARTICLE X3\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].x3;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write v1 */
  fprintf(fp,"\nPARTICLE V1\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].v1;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write v2 */
  fprintf(fp,"\nPARTICLE V2\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].v2;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write v3 */
  fprintf(fp,"\nPARTICLE V3\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    buf[nbuf++] = pG->particle[p].v3;
    if ((nbuf+1) > bufsize) {
      fwrite(buf,sizeof(Real),nbuf,fp);
      nbuf = 0;
    }
  }
  if (nbuf > 0) {
    fwrite(buf,sizeof(Real),nbuf,fp);
    nbuf = 0;
  }

  /* Write properties */
  fprintf(fp,"\nPARTICLE PROPERTY\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    ibuf[nibuf++] = pG->particle[p].property;
    if ((nibuf+1) > ibufsize) {
      fwrite(ibuf,sizeof(int),nibuf,fp);
      nibuf = 0;
    }
  }
  if (nibuf > 0) {
    fwrite(ibuf,sizeof(int),nibuf,fp);
    nibuf = 0;
  }

  /* Write my_id */
  fprintf(fp,"\nPARTICLE MY_ID\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    lbuf[nlbuf++] = pG->particle[p].my_id;
    if ((nlbuf+1) > lbufsize) {
      fwrite(lbuf,sizeof(long),nlbuf,fp);
      nlbuf = 0;
    }
  }
  if (nlbuf > 0) {
    fwrite(lbuf,sizeof(long),nlbuf,fp);
    nlbuf = 0;
  }

#ifdef MPI_PARALLEL
  /* Write init_id */
  fprintf(fp,"\nPARTICLE INIT_ID\n");
  for (p=0;p<pG->nparticle;p++)
  if (pG->particle[p].pos == 1){
    ibuf[nibuf++] = pG->particle[p].init_id;
    if ((nibuf+1) > ibufsize) {
      fwrite(ibuf,sizeof(int),nibuf,fp);
      nibuf = 0;
    }
  }
  if (nibuf > 0) {
    fwrite(ibuf,sizeof(int),nibuf,fp);
    nibuf = 0;
  }
#endif

#endif /*PARTICLES*/

  fprintf(fp,"\nUSER_DATA\n");

/* call a user function to write his/her problem-specific data! */
  problem_write_restart(pG, pD, fp);

  fclose(fp);

  free_1d_array(buf);

  return;
}
