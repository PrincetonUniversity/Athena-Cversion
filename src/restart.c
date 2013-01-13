#include "copyright.h"
/*============================================================================*/
/*! \file restart.c
 *  \brief Functions for writing and reading restart files.
 *
 * PURPOSE: Functions for writing and reading restart files.  Restart files
 *   begin with the input parameter file (athinput.XX) used to start the run in
 *   text format, followed by binary format data for selected quantities (time,
 *   dt, etc) for each of the Grid structures being updated by this processor.
 *   Lastly, any problem-specific data in binary format is included, which must
 *   be written by a user function in the problem generator.  Since all of the
 *   mesh data (NLevels, number of Domains, size of Grids, etc.) can be
 *   recalculated from the athinput file, only the dependent variables (time,
 *   dt, nstep, U, B, etc.) are dumped in restart files.
 *
 *   The athinput file included at the start of the restart file is parsed by
 *   main() in Step 2, using the standard par_* functions, with the parameters
 *   used by init_mesh() and init_grid() to re-initialize the grid.  Parameter
 *   values read from the athinput file contained in the resfile can be
 *   superceded by input from the command line, or another input file.
 *
 * MPI parallel jobs must be restarted on the same number of processors as they
 * were run originally.
 *
 * With SMR, restart files contain ALL levels and domains being updated by each
 * processor in one file, written in the default directory for the process.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - restart_grids() - reads nstep,time,dt,ConsS and B from restart file 
 * - dump_restart()  - writes a restart file
 *									      */
/*============================================================================*/

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
/*! \fn void restart_grids(char *res_file, MeshS *pM)
 *  \brief Reads nstep, time, dt, and arrays of ConsS and interface B
 *   for each of the Grid structures in the restart file.  
 *
 *    By the time this
 *   function is called (in Step 6 of main()), the Mesh hierarchy has already
 *   been re-initialized by init_mesh() and init_grid() in Step 4 of main()
 *   using parameters in the athinput file at the start of this restart file,
 *   the command line, or from a new input file.
 */

void restart_grids(char *res_file, MeshS *pM)
{
  GridS *pG;
  FILE *fp;
  char line[MAXLEN];
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
#if defined(MHD) || defined(RADIATION_MHD)
  int ib=0,jb=0,kb=0;
#endif
#if (NSCALARS > 0)
  int n;
  char scalarstr[16];
#endif
#ifdef PARTICLES
  long p;
#endif
#ifdef RADIATION_TRANSFER
  RadGridS *pRG;
  int nf, nang, noct, l, m, ifr;
#endif

/* Open the restart file */

  if((fp = fopen(res_file,"r")) == NULL)
    ath_error("[restart_grids]: Error opening the restart file\nIf this is a MPI job, make sure each file from each processor is in the same directory.\n");

/* Skip over the parameter file at the start of the restart file */

  do{
    fgets(line,MAXLEN,fp);
  }while(strncmp(line,"<par_end>",9) != 0);

/* read nstep */

  fgets(line,MAXLEN,fp);
  if(strncmp(line,"N_STEP",6) != 0)
    ath_error("[restart_grids]: Expected N_STEP, found %s",line);
  fread(&(pM->nstep),sizeof(int),1,fp);

/* read time */

  fgets(line,MAXLEN,fp);   /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME",4) != 0)
    ath_error("[restart_grids]: Expected TIME, found %s",line);
  fread(&(pM->time),sizeof(Real),1,fp);

/* read dt */

  fgets(line,MAXLEN,fp);    /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"TIME_STEP",9) != 0)
    ath_error("[restart_grids]: Expected TIME_STEP, found %s",line);
  fread(&(pM->dt),sizeof(Real),1,fp);
#ifdef STS
  fread(&(pM->diff_dt),sizeof(Real),1,fp);
  fread(&(N_STS),sizeof(int),1,fp);
  fread(&(nu_STS),sizeof(Real),1,fp);
#endif

/* Now loop over all Domains containing a Grid on this processor */

  for (nl=0; nl<=(pM->NLevels)-1; nl++){
  for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;
      is = pG->is;
      ie = pG->ie;
      js = pG->js;
      je = pG->je;
      ks = pG->ks;
      ke = pG->ke;

/* propagate time and dt to all Grids */

      pG->time = pM->time;
      pG->dt   = pM->dt;

/* Read the density */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"DENSITY",7) != 0)
        ath_error("[restart_grids]: Expected DENSITY, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].d),sizeof(Real),1,fp);
          }
        }
      }

/* Read the x1-momentum */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"1-MOMENTUM",10) != 0)
        ath_error("[restart_grids]: Expected 1-MOMENTUM, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].M1),sizeof(Real),1,fp);
          }
        }
      }

/* Read the x2-momentum */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"2-MOMENTUM",10) != 0)
        ath_error("[restart_grids]: Expected 2-MOMENTUM, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].M2),sizeof(Real),1,fp);
          }
        }
      }

/* Read the x3-momentum */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"3-MOMENTUM",10) != 0)
        ath_error("[restart_grids]: Expected 3-MOMENTUM, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].M3),sizeof(Real),1,fp);
          }
        }
      }

#ifndef BAROTROPIC
/* Read energy density */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"ENERGY",6) != 0)
        ath_error("[restart_grids]: Expected ENERGY, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].E),sizeof(Real),1,fp);
          }
        }
      }
#endif


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
/* Read radiation energy density */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"radENERGY",9) != 0)
        ath_error("[restart_grids]: Expected radENERGY, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Er),sizeof(Real),1,fp);
          }
        }
      }


/* Read the 1-radiation flux */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"1-radFLUX",9) != 0)
        ath_error("[restart_grids]: Expected 1-radFLUX, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Fr1),sizeof(Real),1,fp);
          }
        }
      }

/* Read the 2-radiation Flux */

     fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"2-radFLUX",9) != 0)
        ath_error("[restart_grids]: Expected 2-radFLUX, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Fr2),sizeof(Real),1,fp);
          }
        }
      }
/* Read the 3-radiation Flux */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"3-radFLUX",9) != 0)
        ath_error("[restart_grids]: Expected 3-radFLUX, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Fr3),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Eddington tensor 11 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"11-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 11-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_11),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Eddington tensor 21 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"21-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 21-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_21),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Eddington tensor 22 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"22-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 22-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_22),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Eddington tensor 31 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"31-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 31-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_31),sizeof(Real),1,fp);
          }
        }
      }


/* Read the Eddington tensor 32 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"32-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 32-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_32),sizeof(Real),1,fp);
          }
        }
      }


/* Read the Eddington tensor 33 */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"33-Eddtensor",12) != 0)
        ath_error("[restart_grids]: Expected 33-Eddtensor, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Edd_33),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Flux mean scattering opacity */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Sigma_sF",8) != 0)
        ath_error("[restart_grids]: Expected Sigma_sF, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Sigma[0]),sizeof(Real),1,fp);
          }
        }
      }


/* Read the flux mean absorption opacity */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Sigma_aF",8) != 0)
        ath_error("[restart_grids]: Expected Sigma_aF, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Sigma[1]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Plank mean absorption opacity */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Sigma_aP",8) != 0)
        ath_error("[restart_grids]: Expected Sigma_aP, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Sigma[2]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the Energy mean absorption opacity */

    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Sigma_aE",8) != 0)
        ath_error("[restart_grids]: Expected Sigma_aE, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->U[k][j][i].Sigma[3]),sizeof(Real),1,fp);
          }
        }
      }


#endif 
/* End read radiation quantity */

#if defined(MHD) || defined(RADIATION_MHD)
/* if there is more than one cell in each dimension, need to read one more
 * face-centered field component than the number of cells.  [ijk]b is
 * the number of extra cells to be read  */

      if (ie > is) ib=1;
      if (je > js) jb=1;
      if (ke > ks) kb=1;

/* Read the face-centered x1 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"1-FIELD",7) != 0)
        ath_error("[restart_grids]: Expected 1-FIELD, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie+ib; i++) {
            fread(&(pG->B1i[k][j][i]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the face-centered x2 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"2-FIELD",7) != 0)
        ath_error("[restart_grids]: Expected 2-FIELD, found %s",line);
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je+jb; j++) {
          for (i=is; i<=ie; i++) {
            fread(&(pG->B2i[k][j][i]),sizeof(Real),1,fp);
          }
        }
      }

/* Read the face-centered x3 B-field */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"3-FIELD",7) != 0)
        ath_error("[restart_grids]: Expected 3-FIELD, found %s",line);
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
        pG->U[k][j][i].B1c = pG->B1i[k][j][i];
        pG->U[k][j][i].B2c = pG->B2i[k][j][i];
        pG->U[k][j][i].B3c = pG->B3i[k][j][i];
        if(ib==1) pG->U[k][j][i].B1c=0.5*(pG->B1i[k][j][i] +pG->B1i[k][j][i+1]);
        if(jb==1) pG->U[k][j][i].B2c=0.5*(pG->B2i[k][j][i] +pG->B2i[k][j+1][i]);
        if(kb==1) pG->U[k][j][i].B3c=0.5*(pG->B3i[k][j][i] +pG->B3i[k+1][j][i]);
      }}}
#endif

#if (NSCALARS > 0)
/* Read any passively advected scalars */
/* Following code only works if NSCALARS < 10 */

      for (n=0; n<NSCALARS; n++) {
        fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
        fgets(line,MAXLEN,fp);
        sprintf(scalarstr, "SCALAR %d", n);
        if(strncmp(line,scalarstr,8) != 0)
          ath_error("[restart_grids]: Expected %s, found %s",scalarstr,line);
        for (k=ks; k<=ke; k++) {
          for (j=js; j<=je; j++) {
            for (i=is; i<=ie; i++) {
              fread(&(pG->U[k][j][i].s[n]),sizeof(Real),1,fp);
            }
          }
        }
      }
#endif

#ifdef PARTICLES
/* Read particle properties and the complete particle list */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE LIST",13) != 0)
        ath_error("[restart_grids]: Expected PARTICLE LIST, found %s",line);
      fread(&(pG->nparticle),sizeof(long),1,fp);

      if (pG->nparticle > pG->arrsize-2)
        particle_realloc(pG, pG->nparticle+2);

      fread(&(npartypes),sizeof(int),1,fp);
      for (i=0; i<npartypes; i++) {          /* particle property list */
#ifdef FEEDBACK
        fread(&(grproperty[i].m),sizeof(Real),1,fp);
#endif
        fread(&(grproperty[i].rad),sizeof(Real),1,fp);
        fread(&(grproperty[i].rho),sizeof(Real),1,fp);
        fread(&(tstop0[i]),sizeof(Real),1,fp);
        fread(&(grrhoa[i]),sizeof(Real),1,fp);
      }
      fread(&(alamcoeff),sizeof(Real),1,fp);  /* coef to calc Reynolds number */

      for (i=0; i<npartypes; i++)
        fread(&(grproperty[i].integrator),sizeof(short),1,fp);

/* Read the x1-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X1",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE X1, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x1),sizeof(Real),1,fp);
      }

/* Read the x2-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X2",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE X2, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x2),sizeof(Real),1,fp);
      }

/* Read the x3-positions */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE X3",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE X3, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].x3),sizeof(Real),1,fp);
      }

/* Read the v1 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V1",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE V1, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v1),sizeof(Real),1,fp);
      }

/* Read the v2 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V2",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE V2, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v2),sizeof(Real),1,fp);
      }

/* Read the v3 velocity */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE V3",11) != 0)
        ath_error("[restart_grids]: Expected PARTICLE V3, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].v3),sizeof(Real),1,fp);
      }

/* Read particle properties */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE PROPERTY",17) != 0)
        ath_error("[restart_grids]: Expected PARTICLE PROPERTY, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].property),sizeof(int),1,fp);
        pG->particle[p].pos = 1;	/* grid particle */
      }

/* Read particle my_id */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE MY_ID",14) != 0)
        ath_error("[restart_grids]: Expected PARTICLE MY_ID, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].my_id),sizeof(long),1,fp);
      }

#ifdef MPI_PARALLEL
/* Read particle init_id */

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"PARTICLE INIT_ID",16) != 0)
        ath_error("[restart_grids]: Expected PARTICLE INIT_ID, found %s",line);
      for (p=0; p<pG->nparticle; p++) {
        fread(&(pG->particle[p].init_id),sizeof(int),1,fp);
      }
#endif

/* count the number of particles with different types */

      for (i=0; i<npartypes; i++)
        grproperty[i].num = 0;
      for (p=0; p<pG->nparticle; p++)
        grproperty[pG->particle[p].property].num += 1;

#endif /* PARTICLES */

    }


    //#ifdef NOT_DEFINED_EVER
#ifdef RADIATION_TRANSFER
    pRG=pM->Domain[nl][nd].RadGrid;
    is = pRG->is;
    ie = pRG->ie;
    js = pRG->js;
    je = pRG->je;
    ks = pRG->ks;
    ke = pRG->ke;
    nf = pRG->nf;
    nang = pRG->nang;
    noct = pRG->noct;
    printf("noct, nang: %d %d\n",nang,noct);

/* Read the mean intensity */
    fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
    fgets(line,MAXLEN,fp);
    if(strncmp(line,"rad_J",5) != 0)
      ath_error("[restart_grids]: Expected rad_J, found %s",line);
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
	for (i=is; i<=ie; i++) {
	  for (ifr=0; ifr<nf; ifr++) {
	    fread(&(pRG->R[ifr][k][j][i].J),sizeof(Real),1,fp);
	  }}}}

    if (pRG->Nx[0] > 1) {
/* Read intensities on ix1 grid face */
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstl1i",7) != 0)
	ath_error("[restart_grids]: Expected Ghstl1i, found %s",line);
      for (ifr=0; ifr<nf; ifr++) {
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (l=0; l<noct; l++) {
	      for (m=0; m<nang; m++) {
		fread(&(pRG->Ghstl1i[ifr][k][j][l][m]),sizeof(Real),1,fp);
	      }}}}}

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"l1imu",5) != 0)
	ath_error("[restart_grids]: Expected l1imu, found %s",line);
      for (ifr=0; ifr<nf; ifr++) {
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (l=0; l<noct; l++) {
	      for (m=0; m<nang; m++) {
		fread(&(pRG->l1imu[ifr][k][j][l][m]),sizeof(Real),1,fp);
	      }}}}}

/* Read intensities on ox1 grid face */
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */     
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstr1i",7) != 0)
       ath_error("[restart_grids]: Expected Ghstr1i, found %s",line);
      for (ifr=0; ifr<nf; ifr++) {
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (l=0; l<noct; l++) {
	      for (m=0; m<nang; m++) {
		fread(&(pRG->Ghstr1i[ifr][k][j][l][m]),sizeof(Real),1,fp);
	      }}}}}

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */     
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"r1imu",5) != 0)
        ath_error("[restart_grids]: Expected r1imu, found %s",line);
      for (ifr=0; ifr<nf; ifr++) {
	for (k=ks; k<=ke; k++) {
	  for (j=js; j<=je; j++) {
	    for (l=0; l<noct; l++) {
	      for (m=0; m<nang; m++) {
		fread(&(pRG->r1imu[ifr][k][j][l][m]),sizeof(Real),1,fp);
	      }}}}}
    }

    if (pRG->Nx[1] > 1) {
/* Read intensities on ix2 grid face */
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstl2i",7) != 0)
        ath_error("[restart_grids]: Expected Ghstl2i, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->Ghstl2i[ifr][k][i][l][m]),sizeof(Real),1,fp);
		}}}}}

      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"l2imu",5) != 0)
        ath_error("[restart_grids]: Expected l2imu, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->l2imu[ifr][k][i][l][m]),sizeof(Real),1,fp);
		}}}}}

/* Read intensities on ox2 grid face */
     fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstr2i",7) != 0)
        ath_error("[restart_grids]: Expected Ghstr2i, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->Ghstr2i[ifr][k][i][l][m]),sizeof(Real),1,fp);
		}}}}}

     fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"r2imu",5) != 0)
        ath_error("[restart_grids]: Expected r2imu, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->r2imu[ifr][k][i][l][m]),sizeof(Real),1,fp);
		}}}}}
      }

      if (pRG->Nx[2] > 1) {
/* Read intensities on ix3 grid face */
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstl3i",7) != 0)
        ath_error("[restart_grids]: Expected Ghstl3i, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->Ghstl3i[ifr][j][i][l][m]),sizeof(Real),1,fp);
		}}}}}
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"l3imu",5) != 0)
        ath_error("[restart_grids]: Expected l3imu, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->l3imu[ifr][j][i][l][m]),sizeof(Real),1,fp);
		}}}}}

/* Read intensities on ox3 grid face */
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"Ghstr3i",7) != 0)
        ath_error("[restart_grids]: Expected Ghstr3i, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->Ghstr3i[ifr][j][i][l][m]),sizeof(Real),1,fp);
		}}}}}
      fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
      fgets(line,MAXLEN,fp);
      if(strncmp(line,"r3imu",5) != 0)
        ath_error("[restart_grids]: Expected r3imu, found %s",line);
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  fread(&(pRG->r3imu[ifr][j][i][l][m]),sizeof(Real),1,fp);
		}}}}}
      }
#endif /*RADIATION_TRANSFER*/

  }} /* End loop over all Domains --------------------------------------------*/

/* Call a user function to read his/her problem-specific data! */

  fgets(line,MAXLEN,fp); /* Read the '\n' preceeding the next string */
  fgets(line,MAXLEN,fp);
  if(strncmp(line,"USER_DATA",9) != 0)
    ath_error("[restart_grids]: Expected USER_DATA, found %s",line);
  problem_read_restart(pM, fp);

  fclose(fp);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_restart(MeshS *pM, OutputS *pout)
 *  \brief Writes a restart file, including problem-specific data from
 *   a user defined function  */

void dump_restart(MeshS *pM, OutputS *pout)
{
  GridS *pG;
  FILE *fp;
  char *fname;
  int i,j,k,is,ie,js,je,ks,ke,nl,nd;
#if defined(MHD) || defined(RADIATION_MHD)
  int ib=0,jb=0,kb=0;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef PARTICLES
  int nprop, *ibuf = NULL, ibufsize, lbufsize, sbufsize, nibuf, nlbuf, nsbuf;
  long np, p, *lbuf = NULL;
  short *sbuf = NULL;
  nibuf = 0;	nsbuf = 0;	nlbuf = 0;
#endif
  int bufsize, nbuf = 0;
  Real *buf = NULL;

#ifdef RADIATION_TRANSFER
  RadGridS *pRG;
  int nf, nang, noct, l, m, ifr;
#endif

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

/* Create filename and Open the output file */

  if((fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,
      pout->num,NULL,"rst")) == NULL){
    ath_error("[dump_restart]: Error constructing filename\n");
  }

  if((fp = fopen(fname,"wb")) == NULL){
    ath_error("[dump_restart]: Unable to open restart file\n");
    return;
  }
  free(fname);

/* Add the current time & nstep to the parameter file */

  par_setd("time","time","%e",pM->time,"Current Simulation Time");
  par_seti("time","nstep","%d",pM->nstep,"Current Simulation Time Step");

/* Write the current state of the parameter file */

  par_dump(2,fp);

/* Write out the current simulation step number */

  fprintf(fp,"N_STEP\n");
  if(fwrite(&(pM->nstep),sizeof(int),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time */

  fprintf(fp,"\nTIME\n");
  if(fwrite(&(pM->time),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");

/* Write out the current simulation time step */

  fprintf(fp,"\nTIME_STEP\n");
  if(fwrite(&(pM->dt),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");
#ifdef STS
  if(fwrite(&(pM->diff_dt),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");
  if(fwrite(&(N_STS),sizeof(int),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");
  if(fwrite(&(nu_STS),sizeof(Real),1,fp) != 1)
    ath_error("[dump_restart]: fwrite() error\n");
#endif

/* Now loop over all Domains containing a Grid on this processor */

  for (nl=0; nl<=(pM->NLevels)-1; nl++){
  for (nd=0; nd<=(pM->DomainsPerLevel[nl])-1; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pG=pM->Domain[nl][nd].Grid;
      is = pG->is;
      ie = pG->ie;
      js = pG->js;
      je = pG->je;
      ks = pG->ks;
      ke = pG->ke;

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


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)

/* Write radiation energy density Er */

      fprintf(fp,"\nradENERGY\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Er;
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


/* Write the x1-radFlux */

      fprintf(fp,"\n1-radFLUX\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Fr1;
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

/* Write the x2-radFlux */

      fprintf(fp,"\n2-radFLUX\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Fr2;
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
    
/* Write the x3-radFlux */

      fprintf(fp,"\n3-radFLUX\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Fr3;
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

/* Write the Eddington tensor 11 */

      fprintf(fp,"\n11-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_11;
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

/* Write the Eddington tensor 21 */

      fprintf(fp,"\n21-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_21;
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

/* Write the Eddington tensor 22 */

      fprintf(fp,"\n22-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_22;
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

/* Write the Eddington tensor 31 */

      fprintf(fp,"\n31-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_31;
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

/* Write the Eddington tensor 32 */

      fprintf(fp,"\n32-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_32;
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


/* Write the Eddington tensor 33 */

      fprintf(fp,"\n33-Eddtensor\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Edd_33;
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

/* Write the flux mean scattering opacity */

      fprintf(fp,"\nSigma_sF\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Sigma[0];
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

/* Write flux mean absorption opacity */

      fprintf(fp,"\nSigma_aF\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Sigma[1];
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

/* Write Planck mean absorption opacity */

      fprintf(fp,"\nSigma_aP\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Sigma[2];
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

/* Write Energy mean absorption opacity */

      fprintf(fp,"\nSigma_aE\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            buf[nbuf++] = pG->U[k][j][i].Sigma[3];
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
/* Finish output radiation quantity */


#if defined(MHD) || defined(RADIATION_MHD)
/* see comments in restart_grid_block() for use of [ijk]b */

  if (ie > is) ib = 1;
  if (je > js) jb = 1;
  if (ke > ks) kb = 1;

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
      fwrite(&(npartypes),sizeof(int),1,fp); /* number of particle types */
      for (i=0; i<npartypes; i++) {          /* particle property list */
#ifdef FEEDBACK
        buf[nbuf++] = grproperty[i].m;
#endif
        buf[nbuf++] = grproperty[i].rad;
        buf[nbuf++] = grproperty[i].rho;
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
      fwrite(&(alamcoeff),sizeof(Real),1,fp);  /* coef for Reynolds number */
    
      for (i=0; i<npartypes; i++) {         /* particle integrator type */
        sbuf[nsbuf++] = grproperty[i].integrator;
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


#ifdef RADIATION_TRANSFER
      pRG=pM->Domain[nl][nd].RadGrid;
      is = pRG->is;
      ie = pRG->ie;
      js = pRG->js;
      je = pRG->je;
      ks = pRG->ks;
      ke = pRG->ke;
      nf = pRG->nf;
      nang = pRG->nang;
      noct = pRG->noct;

/* Write the mean intensity */

     fprintf(fp,"\nrad_J\n");
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
	    for (ifr=0; ifr<nf; ifr++) {
	      buf[nbuf++] = pRG->R[ifr][k][j][i].J;
	      if ((nbuf+1) > bufsize) {
		fwrite(buf,sizeof(Real),nbuf,fp);
		nbuf = 0;
	      }
	    }}}}
      if (nbuf > 0) {
        fwrite(buf,sizeof(Real),nbuf,fp);
        nbuf = 0;
      }

      if (pRG->Nx[0] > 1) {
/* Write intensities on ix1 grid face */
	fprintf(fp,"\nGhstl1i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstl1i[ifr][k][j][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}

	fprintf(fp,"\nl1imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->l1imu[ifr][k][j][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}

/* Write intensities on ox1 grid face */
	fprintf(fp,"\nGhstr1i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstr1i[ifr][k][j][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}

	fprintf(fp,"\nr1imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->r1imu[ifr][k][j][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
      }

      if (pRG->Nx[1] > 1) {
/* Write intensities on ix2 grid face */
	fprintf(fp,"\nGhstl2i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstl2i[ifr][k][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}

	fprintf(fp,"\nl2imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->l2imu[ifr][k][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
/* Write intensities on ox2 grid face */
	fprintf(fp,"\nGhstr2i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstr2i[ifr][k][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}

	fprintf(fp,"\nr2imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (k=ks; k<=ke; k++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->r2imu[ifr][k][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
      }

      if (pRG->Nx[2] > 1) {
/* Write intensities on ix3 grid face */
	fprintf(fp,"\nGhstl3i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstl3i[ifr][j][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
	fprintf(fp,"\nl3imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->l3imu[ifr][j][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
/* Write intensities on ox3 grid face */
	fprintf(fp,"\nGhstr3i\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->Ghstr3i[ifr][j][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
	fprintf(fp,"\nr3imu\n");
	for (ifr=0; ifr<nf; ifr++) {
	  for (j=js; j<=je; j++) {
	    for (i=is; i<=ie; i++) {
	      for (l=0; l<noct; l++) {
		for (m=0; m<nang; m++) {
		  buf[nbuf++] = pRG->r3imu[ifr][j][i][l][m];
		  if ((nbuf+1) > bufsize) {
		    fwrite(buf,sizeof(Real),nbuf,fp);
		    nbuf = 0;
		  }
		}}}}}
	if (nbuf > 0) {
	  fwrite(buf,sizeof(Real),nbuf,fp);
	  nbuf = 0;
	}
      }

#endif /*RADIATION_TRANSFER*/

    }
  }}  /*---------- End loop over all Domains ---------------------------------*/
    
/* call a user function to write his/her problem-specific data! */
    
  fprintf(fp,"\nUSER_DATA\n");
  problem_write_restart(pM, fp);

  fclose(fp);

  free_1d_array(buf);

  return;
}
