#include "../copyright.h"
/*===========================================================================*/
/*! \file output_particle.c
 *  \brief Contains routines necessary for outputting particles.
 *
 * PURPOSE: contains routines necessary for outputting particles.
 *   There are two basic formats:
 *  - 1. Bin particles to the grid, then output particles as a grid.
 *  - 2. Dump the particle list directly.
 *
 *   For particle binning, there can be many choices since particles may have
 *   different properties. We provide a default (and trivial) particle selection
 *   function, which select all the particles with different properties. The 
 *   user can define their own particle selection functions in the problem
 *   generator and pass them to the main code.
 *
 *   The output quantities include, density, momentum density and velocity of
 *   the selected particles averaged in one grid cell. The binned data are
 *   saved in arrays dpar, and grid_v. The latter is borrowed from particle.c
 *   to save memory. The expression functions expr_??? are used to pick the 
 *   relevant quantities, which is part of the output data structure. The way
 *   to output these binned particle quantities are then exactly the same as
 *   other gas quantities.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - particle_to_grid();
 * - dump_particle_binary();
 * - property_all();
 * 
 *============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

#ifdef PARTICLES         /* endif at the end of the file */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   expr_*
 *============================================================================*/

Real expr_dpar (const GridS *pG, const int i, const int j, const int k);
Real expr_M1par(const GridS *pG, const int i, const int j, const int k);
Real expr_M2par(const GridS *pG, const int i, const int j, const int k);
Real expr_M3par(const GridS *pG, const int i, const int j, const int k);
Real expr_V1par(const GridS *pG, const int i, const int j, const int k);
Real expr_V2par(const GridS *pG, const int i, const int j, const int k);
Real expr_V3par(const GridS *pG, const int i, const int j, const int k);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*! \fn void particle_to_grid(Grid *pG, PropFun_t par_prop)
 *  \brief Bin the particles to grid cells
 */
void particle_to_grid(DomainS *pD, PropFun_t par_prop)
{
  GridS *pG = pD->Grid;
  int i,j,k, is,js,ks, i0,j0,k0, i1,j1,k1, i2,j2,k2;
  int n0 = ncell-1;
  long p;
  Real drho;
  Real weight[3][3][3];
  Real3Vect cell1;
  GrainS *gr;

  /* Get grid limit related quantities */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else                cell1.x1 = 0.0;

  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else                cell1.x2 = 0.0;

  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else                cell1.x3 = 0.0;

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        pG->Coup[k][j][i].grid_d = 0.0;
        pG->Coup[k][j][i].grid_v1 = 0.0;
        pG->Coup[k][j][i].grid_v2 = 0.0;
        pG->Coup[k][j][i].grid_v3 = 0.0;
      }

  /* bin the particles */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);

    /* judge if the particle should be selected */
    if ((*par_prop)(gr, &(pG->parsub[p]))) {/* 1: true; 0: false */

      getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);

      /* distribute particles */
      k1 = MAX(ks, klp);    k2 = MIN(ks+n0, kup);
      j1 = MAX(js, jlp);    j2 = MIN(js+n0, jup);
      i1 = MAX(is, ilp);    i2 = MIN(is+n0, iup);

      for (k=k1; k<=k2; k++) {
        k0 = k-k1;
        for (j=j1; j<=j2; j++) {
          j0 = j-j1;
          for (i=i1; i<=i2; i++) {
            i0 = i-i1;
            /* interpolate the particles to the grid */
#ifdef FEEDBACK
            drho = grproperty[gr->property].m;
#else
            drho = 1.0;
#endif
            pG->Coup[k][j][i].grid_d  += weight[k0][j0][i0]*drho;
            pG->Coup[k][j][i].grid_v1 += weight[k0][j0][i0]*drho*gr->v1;
            pG->Coup[k][j][i].grid_v2 += weight[k0][j0][i0]*drho*gr->v2;
            pG->Coup[k][j][i].grid_v3 += weight[k0][j0][i0]*drho*gr->v3;

          }
        }
      }
    }
  }

/* deposit ghost zone values into the boundary zones */
  exchange_gpcouple(pD, 0);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void dump_particle_binary(MeshS *pM, OutputS *pOut)
 *  \brief Dump unbinned particles in binary format
 */
void dump_particle_binary(MeshS *pM, OutputS *pOut)
{
  DomainS *pD = (DomainS*)&(pM->Domain[0][0]);  
  GridS   *pG = pD->Grid; 
  int dnum = pOut->num;
  FILE *pfile;
  char *fname;
  long p, nout, my_id;
  int i,is,js,ks,h,init_id = 0;
  short pos;
  Real3Vect cell1;
  Real weight[3][3][3];         /* weight function */
  Real dpar,u1,u2,u3,cs;
#ifdef FEEDBACK
  Real stiffness;
#endif
  GrainS *gr;
  float fdata[12];  /* coordinate of grid and domain boundary */

  if((fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,
      pOut->num,pOut->id,"lis")) == NULL){
    ath_error("[dump_particle_binary]: Error constructing filename\n");
  }

  /* open output file */
  if((pfile = fopen(fname,"wb")) == NULL){
    ath_error("[dump_particle_binary]: Unable to open lis file %s\n",fname);
  }

  /* bin all the particles to the grid */
  particle_to_grid(pD, property_all);

  /* Get grid limit related quantities */
  if (pG->Nx[0] > 1)  cell1.x1 = 1.0/pG->dx1;
  else               cell1.x1 = 0.0;
  if (pG->Nx[1] > 1)  cell1.x2 = 1.0/pG->dx2;
  else               cell1.x2 = 0.0;
  if (pG->Nx[2] > 1)  cell1.x3 = 1.0/pG->dx3;
  else               cell1.x3 = 0.0;

  /* update the particle auxilary array */
  for (p=0; p<pG->nparticle; p++)
  {
    gr = &(pG->particle[p]);

    /* get the local particle density */
    getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);
#ifdef FEEDBACK
    h = getvalues(pG, weight, is, js, ks, &dpar,&u1,&u2,&u3,&cs,&stiffness);
#else
    h = getvalues(pG, weight, is, js, ks, &dpar, &u1, &u2, &u3, &cs);
#endif

    pG->parsub[p].dpar = dpar;
  }

  /* find out how many particles is to be output */
  nout = 0;
  for (p=0; p<pG->nparticle; p++)
  if ((*(pOut->par_prop))(&(pG->particle[p]), &(pG->parsub[p])))
    nout += 1;

/* write the basic information */

  /* Write the grid and domain boundary */
  fdata[0]  = (float)(pG->MinX[0]);
  fdata[1]  = (float)(pG->MaxX[0]);
  fdata[2]  = (float)(pG->MinX[1]);
  fdata[3]  = (float)(pG->MaxX[1]);
  fdata[4]  = (float)(pG->MinX[2]);
  fdata[5]  = (float)(pG->MaxX[2]);
  fdata[6]  = (float)(pM->RootMinX[0]);
  fdata[7]  = (float)(pM->RootMaxX[0]);
  fdata[8]  = (float)(pM->RootMinX[1]);
  fdata[9]  = (float)(pM->RootMaxX[1]);
  fdata[10] = (float)(pM->RootMinX[2]);
  fdata[11] = (float)(pM->RootMaxX[2]);

  fwrite(fdata,sizeof(float),12,pfile);

  /* Write particle property information */
  fwrite(&(npartypes),sizeof(int),1,pfile);

  for (i=0; i<npartypes; i++)
  {
    fdata[0] = (float)(grproperty[i].rad);
    fwrite(fdata,sizeof(float),1,pfile);
  }

  /* Write time, dt */
  fdata[0] = (float)pG->time;
  fdata[1] = (float)pG->dt;
  fwrite(fdata,sizeof(float),2,pfile);

/* Write all the selected particles */

  /* Write particle number */
  fwrite(&nout,sizeof(long),1,pfile);

  /* Write particle information */
  for (p=0; p<pG->nparticle; p++)
  {
    gr = &(pG->particle[p]);
    if ((*(pOut->par_prop))(gr,&(pG->parsub[p]))) { /* 1: true; 0: false */

      /* collect data */
      fdata[0] = (float)(gr->x1);
      fdata[1] = (float)(gr->x2);
      fdata[2] = (float)(gr->x3);
      fdata[3] = (float)(gr->v1);
      fdata[4] = (float)(gr->v2);
      fdata[5] = (float)(gr->v3);
//      fdata[6] = (float)(pG->grproperty[gr->property].rad);
      fdata[6] = (float)(pG->parsub[p].dpar);
      my_id = gr->my_id;
#ifdef MPI_PARALLEL
      init_id = gr->init_id;
#endif
      pos = gr->pos;

      fwrite(fdata,sizeof(float),7,pfile);
      fwrite(&(gr->property),sizeof(int),1,pfile);
      fwrite(&(my_id),sizeof(long),1,pfile);
      fwrite(&(init_id),sizeof(int),1,pfile);

    }
  }

  fclose(pfile);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn int property_all(const Grain *gr, const GrainAux *grsub)
 *  \brief Default choice for binning particles to the grid: 
 * All the particles are binned, return true for any value.
 */
int property_all(const GrainS *gr, const GrainAux *grsub)
{
//  if (gr->pos == 1)
    return 1;  /* always true */
//  else
//    return 0;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */

/* expr_*: where * are variables d,M1,M2,M3,V1,V2,V3 for particles */

/*! \fn Real expr_dpar(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle density */
Real expr_dpar(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_d;
}
/*! \fn Real expr_M1par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 1-momentum */
Real expr_M1par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_v1;
}

/*! \fn Real expr_M2par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 2-momentum */
Real expr_M2par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_v2;
}
/*! \fn Real expr_M3par(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle 3-momentum */
Real expr_M3par(const GridS *pG, const int i, const int j, const int k) {
  return pG->Coup[k][j][i].grid_v3;
}
/*! \fn Real expr_V1par(const Grid *pG, const int i, const int j, const int k) 
 *  \brief Wrapper for particle 1-velocity */
Real expr_V1par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_v1/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}
/*! \fn Real expr_V2par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 2-velocity */
Real expr_V2par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_v2/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}
/*! \fn Real expr_V3par(const Grid *pG, const int i, const int j, const int k)
 *  \brief Wrapper for particle 3-velocity */
Real expr_V3par(const GridS *pG, const int i, const int j, const int k) {
  if (pG->Coup[k][j][i].grid_d>0.0)
    return pG->Coup[k][j][i].grid_v3/pG->Coup[k][j][i].grid_d;
  else return 0.0;
}

#endif /*PARTICLES*/
