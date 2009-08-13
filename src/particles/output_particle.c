#include "../copyright.h"
/*=============================================================================
FILE: output_particle.c
PURPOSE: contains all the routines necessary for outputting particles. There are
  two basic formats:
  1. Bin particles to the grid, then output particles as a grid.
  2. Dump the particle list directly.

  For particle binning, there can be many choices since particles may have different
  properties. We provide a default (and trivial) particle selection function, which
  select all the particles with different properties. The user can define their own
  particle selection functions in the problem generator and pass them to the main code.

  The output quantities include, density, momentum density and velocity of the selected
  particles averaged in one grid cell. The binned data are saved in arrays dpar, and grid_v.
  The latter is borrowed from particle.c to save memory. The expression functions expr_???
  are used to pick the relevant quantities, which is part of the output data structure.
  The way to output these binned particle quantities are then exactly the same as other
  gas quantities.

  Dumping particle list has not been developed yet since we need to figure out how to
  do visualization.

CONTAINS PUBLIC FUNCTIONS:
  void particle_to_grid(Grid *pG, Domain *pD, PropFun_t par_prop);
  void dump_particle_binary(Grid *pG, Domain *pD);
  int property_all(Grain *gr);

History:
  Written by Xuening Bai, Mar. 2009

==============================================================================*/
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
 *   property_all
 *============================================================================*/

Real expr_dpar (const Grid *pG, const int i, const int j, const int k);
Real expr_M1par(const Grid *pG, const int i, const int j, const int k);
Real expr_M2par(const Grid *pG, const int i, const int j, const int k);
Real expr_M3par(const Grid *pG, const int i, const int j, const int k);
Real expr_V1par(const Grid *pG, const int i, const int j, const int k);
Real expr_V2par(const Grid *pG, const int i, const int j, const int k);
Real expr_V3par(const Grid *pG, const int i, const int j, const int k);
int property_all(Grain *gr);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/* Bin the particles to grid cells */
void particle_to_grid(Grid *pG, Domain *pD, PropFun_t par_prop)
{
  int i,j,k, is,js,ks, i1,j1,k1;
  long p;
  Real drho;
  Real weight[3][3][3];
  Vector cell1;
  Grain *gr;

  /* Get grid limit related quantities */

  if (pG->Nx1 > 1)  cell1.x1 = 1.0/pG->dx1;
  else              cell1.x1 = 0.0;

  if (pG->Nx2 > 1)  cell1.x2 = 1.0/pG->dx2;
  else              cell1.x2 = 0.0;

  if (pG->Nx3 > 1)  cell1.x3 = 1.0/pG->dx3;
  else              cell1.x3 = 0.0;

  /* initialization */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        grid_d[k][j][i] = 0.0;
        grid_v[k][j][i].x1 = 0.0;
        grid_v[k][j][i].x2 = 0.0;
        grid_v[k][j][i].x3 = 0.0;
      }

  /* bin the particles */
  for (p=0; p<pG->nparticle; p++) {
    gr = &(pG->particle[p]);
    /* judge if the particle should be selected */
    if ((*par_prop)(gr)) {/* 1: true; 0: false */

      getweight(pG, gr->x1, gr->x2, gr->x3, cell1, weight, &is, &js, &ks);

      /* distribute feedback force */
      for (k=0; k<ncell; k++) {
        k1 = k+ks;
        if ((k1 <= kup) && (k1 >= klp)) {
          for (j=0; j<ncell; j++) {
            j1 = j+js;
            if ((j1 <= jup) && (j1 >= jlp)) {
              for (i=0; i<ncell; i++) {
                i1 = i+is;
                if ((i1 <= iup) && (i1 >= ilp)) {
                  /* interpolate the particles to the grid */
#ifdef FEEDBACK
                  drho = pG->grproperty[gr->property].m;
#else
                  drho = 1.0;
#endif
                  grid_d[k1][j1][i1] += weight[k][j][i]*drho;
                  grid_v[k1][j1][i1].x1 += weight[k][j][i]*drho*gr->v1;
                  grid_v[k1][j1][i1].x2 += weight[k][j][i]*drho*gr->v2;
                  grid_v[k1][j1][i1].x3 += weight[k][j][i]*drho*gr->v3;
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}

/* dump unbinned particles in binary format */
void dump_particle_binary(Grid *pG, Domain *pD, Output *pOut)
{
  int dnum = pOut->num;
  FILE *p_binfile;
  char *fname;
  long p, nout, my_id;
  int init_id = 0;
  Grain *gr;
  float fdata[12];  /* coordinate of grid and domain boundary */

  /* open the binary file */
  if((fname = ath_fname(NULL,pG->outfilename,num_digit,dnum,pOut->id,"lis"))
     == NULL){
    ath_error("[dump_particle_binary]: Error constructing filename\n");
    return;
  }

  if((p_binfile = fopen(fname,"wb")) == NULL){
    ath_error("[dump_particle_binary]: Unable to open binary dump file\n");
    return;
  }

  /* find out how many particles is to be output */
  nout = 0;
  for (p=0; p<pG->nparticle; p++)
    if ((*(pOut->par_prop))(&(pG->particle[p]))) /* 1: true; 0: false */
      nout += 1;

/* write the basic information */

  /* write the grid and domain boundary */
  fdata[0]  = (float)(pG->x1_0 + (pG->is + pG->idisp)*pG->dx1);
  fdata[1]  = (float)(pG->x1_0 + (pG->ie+1 + pG->idisp)*pG->dx1);
  fdata[2]  = (float)(pG->x2_0 + (pG->js + pG->jdisp)*pG->dx2);
  fdata[3]  = (float)(pG->x2_0 + (pG->je+1 + pG->jdisp)*pG->dx2);
  fdata[4]  = (float)(pG->x3_0 + (pG->ks + pG->kdisp)*pG->dx3);
  fdata[5]  = (float)(pG->x3_0 + (pG->ke+1 + pG->kdisp)*pG->dx3);
  fdata[6]  = (float)(par_getd("grid","x1min"));
  fdata[7]  = (float)(par_getd("grid","x1max"));
  fdata[8]  = (float)(par_getd("grid","x2min"));
  fdata[9]  = (float)(par_getd("grid","x2max"));
  fdata[10] = (float)(par_getd("grid","x3min"));
  fdata[11] = (float)(par_getd("grid","x3max"));

  fwrite(fdata,sizeof(float),12,p_binfile);

  /* Write time, dt */
  fdata[0] = (float)pG->time;
  fdata[1] = (float)pG->dt;
  fwrite(fdata,sizeof(float),2,p_binfile);

/* Write all the selected particles */

  /* Write particle number */
  fwrite(&nout,sizeof(long),1,p_binfile);

  /* Write particle information */
  for (p=0; p<pG->nparticle; p++)
  {
    gr = &(pG->particle[p]);
    if ((*(pOut->par_prop))(gr)) { /* 1: true; 0: false */

      fdata[0] = (float)(gr->x1);
      fdata[1] = (float)(gr->x2);
      fdata[2] = (float)(gr->x3);
      fdata[3] = (float)(gr->v1);
      fdata[4] = (float)(gr->v2);
      fdata[5] = (float)(gr->v3);
      fdata[6] = (float)(pG->grproperty[gr->property].rad);
#ifdef FEEDBACK
      fdata[7] = (float)(pG->grproperty[gr->property].m);
#else
      fdata[7] = (float)(1.0);
#endif
      my_id = gr->my_id;
#ifdef MPI_PARALLEL
      init_id = gr->init_id;
#endif

      fwrite(fdata,sizeof(float),8,p_binfile);
      fwrite(&(my_id),sizeof(long),1,p_binfile);
      fwrite(&(init_id),sizeof(int),1,p_binfile);

    }
  }

  return;
}

/* default choice for binning particles to the grid: 
   All the particles are binned, return true for any value.
*/
int property_all(Grain *gr)
{
  return 1;  /* always true */
}

/*=========================== PRIVATE FUNCTIONS ==============================*/
/*--------------------------------------------------------------------------- */
/* expr_*: where * are variables d,M1,M2,M3,V1,V2,V3 for particles */

Real expr_dpar(const Grid *pG, const int i, const int j, const int k) {
  return grid_d[k][j][i];
}

Real expr_M1par(const Grid *pG, const int i, const int j, const int k) {
  return grid_v[k][j][i].x1;
}

Real expr_M2par(const Grid *pG, const int i, const int j, const int k) {
  return grid_v[k][j][i].x2;
}

Real expr_M3par(const Grid *pG, const int i, const int j, const int k) {
  return grid_v[k][j][i].x3;
}

Real expr_V1par(const Grid *pG, const int i, const int j, const int k) {
  if (grid_d[k][j][i]>0.0)
    return grid_v[k][j][i].x1/grid_d[k][j][i];
  else return 0.0;
}

Real expr_V2par(const Grid *pG, const int i, const int j, const int k) {
  if (grid_d[k][j][i]>0.0)
    return grid_v[k][j][i].x2/grid_d[k][j][i];
  else return 0.0;
}

Real expr_V3par(const Grid *pG, const int i, const int j, const int k) {
  if (grid_d[k][j][i]>0.0)
    return grid_v[k][j][i].x3/grid_d[k][j][i];
  else return 0.0;
}

#endif /*PARTICLES*/
