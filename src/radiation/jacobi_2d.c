#include "../copyright.h"
/*==============================================================================
 * FILE: jacobi_2d.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 2D grid using jacobi's method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *          Uses quadratic interpolation to compute intensities at edges.  This
 *          requires dividing the sweeps for x-aligned and y-aligned rays.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   formal_solution_2d.c()
 *   formal_solution_2d_destruct()
 *   formal_solution_2d_init()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"


#ifdef RADIATION_TRANSFER
#if defined(JACOBI) && defined(QUADRATIC_INTENSITY)

/* Working arrays used in formal solution */
static Real ***lamstr = NULL;
static Real *****imuo = NULL;
static Real ***Jold = NULL;
static Real *am0;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   update_cell_x() - update grid cell variables for x-aligned rays
 *   update_cell_y() - update grid cell variables for y-aligned rays
 *   update_sfunc() - updates source function after compute of mean intensity
 *   sweep_2d_forward_x() - sweep x aligned rays from lower left to upper right
 *   sweep_2d_backward_x() - sweep x aligned rays from upper right to lower left
 *   sweep_2d_forward_y() - sweep y aligned rays from lower left to upper right
 *   sweep_2d_backward_y() - sweep y aligned rays from upper right to lower left
 *============================================================================*/

static void update_cell_x(RadGridS *pRG, Real *****imuo, int ifr, int k, int j,
                          int i, int l);
static void update_cell_y(RadGridS *pRG, Real *****imuo, int ifr, int k, int j,
                          int i, int l);
static void update_sfunc(RadS *R, Real *dSr, Real lamstr);
static void sweep_2d_forward_y(RadGridS *pRG, int ifr);
static void sweep_2d_backward_y(RadGridS *pRG, int ifr);
static void sweep_2d_forward_x(RadGridS *pRG, int ifr);
static void sweep_2d_backward_x(RadGridS *pRG, int ifr);

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_2d(RadGridS *pRG, Real *dSrmax, int ifr)
 *  \brief formal solution for single freq. in 2D using Jacobi method
 *  with quadratic interpolation of intensities */
void formal_solution_2d(RadGridS *pRG, Real *dSrmax, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf;
  int ismx, jsmx;
  Real dSr, dJ, dJmax;

  if (radt_mode == 2) {
/* Initialize am0 on each call to formal soution when
   both passive and feedback based radiation transfer schemes are enabled */
    for(i=0; i<pRG->nang; i++) {
      am0[i] = fabs(pRG->dx2 * pRG->mu[0][i][0] /(pRG->dx1 * pRG->mu[0][i][1]));
    }
  }

/* if LTE then store J values from previous iteration */
  if(lte != 0) {
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++)
        Jold[ifr][j][i] = pRG->R[ifr][ks][j][i].J;
  }

/* initialize mean intensities at all depths to zero */
  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++) {
      pRG->R[ifr][ks][j][i].J = 0.0;
      pRG->R[ifr][ks][j][i].H[0] = 0.0;
      pRG->R[ifr][ks][j][i].H[1] = 0.0;
      pRG->R[ifr][ks][j][i].K[0] = 0.0;
      pRG->R[ifr][ks][j][i].K[1] = 0.0;
      pRG->R[ifr][ks][j][i].K[2] = 0.0;
      lamstr[ifr][j][i] = 0.0;
    }

/* Compute formal solution and for all rays in each gridzone and
 * update boundary emission*/

  sweep_2d_forward_x(pRG,ifr);
  sweep_2d_forward_y(pRG,ifr);

  sweep_2d_backward_x(pRG,ifr);
  sweep_2d_backward_y(pRG,ifr);

  if(lte == 0) {
/* Update source function */
    (*dSrmax) = 0.0;
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) {
        update_sfunc(&(pRG->R[ifr][ks][j][i]),&dSr,lamstr[ifr][j][i]);
        if( dSr > (*dSrmax)) {
          (*dSrmax) = dSr; ismx=i; jsmx=j;
        }
      }
  } else {
/* Use delta J / J as convergence criterion */
    (*dSrmax) = 0.0;
    dJmax = 0.0;
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) {
        dJ = fabs(pRG->R[ifr][ks][j][i].J - Jold[ifr][j][i]);
        if(dJ > dJmax) dJmax = dJ;
        if (Jold[ifr][j][i] > 0.0)
          dSr = dJ / Jold[ifr][j][i];
        else
          dSr = 0;
        if( dSr > (*dSrmax)) {
          (*dSrmax) = dSr; ismx=i; jsmx=j;
        }
      }
    if(((*dSrmax) == 0.0) && (dJmax > 0.0)) (*dSrmax) = 1.0;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_2d_destruct(void)
 *  \brief free temporary working arrays */
void formal_solution_2d_destruct(void)
{
  if (lamstr != NULL) free_3d_array(lamstr);
  if (imuo   != NULL) free_5d_array(imuo);
  if (Jold   != NULL) free_3d_array(Jold);
  if (am0    != NULL) free_1d_array(am0);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_2d_init(DomainS *pD)
 *  \brief Allocate memory for working arrays */
void formal_solution_2d_init(DomainS *pD)
{

  RadGridS *pRG, *pROG;
  int nx1, nx2, nmx;
  int i, nf, nang;

  if (radt_mode == 0) { /* integration only */
    pRG = pD->RadGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1];
    nf = pRG->nf;
    nang = pRG->nang;
  } else if (radt_mode == 1) { /* output only */
    pRG = pD->RadOutGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1];
    nf = pRG->nf;
    nang = pRG->nang;
  } else if (radt_mode == 2) { /* integration and output */
    pRG = pD->RadGrid;
    pROG = pD->RadOutGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1];
    nf = MAX(pRG->nf,pROG->nf);
    nang = MAX(pRG->nang,pROG->nang);
  }

  nmx = MAX(nx1,nx2);
  if ((imuo = (Real *****)calloc_5d_array(nf,nmx+2,4,nang,2,
                                          sizeof(Real))) == NULL)
    goto on_error;

  if ((am0 = (Real *)calloc_1d_array(nang,sizeof(Real))) == NULL)
    goto on_error;

  if ((lamstr = (Real ***)calloc_3d_array(nf,nx2+2,nx1+2,sizeof(Real))) == NULL)
    goto on_error;

  if (lte != 0)
    if ((Jold = (Real ***)calloc_3d_array(nf,nx2+2,nx1+2,sizeof(Real))) == NULL)
      goto on_error;

  if ( (radt_mode == 0) ||  (radt_mode == 1) ) {
/* Initialize am0 array only once if radiation transfer is only
   used for output or integration but not both */
    for(i=0; i<pRG->nang; i++) {
      am0[i] = fabs(pRG->dx2 * pRG->mu[0][i][0] /(pRG->dx1 * pRG->mu[0][i][1]));
    }
  }

  return;

  on_error:
  formal_solution_2d_destruct();
  ath_error("[formal_solution__2d_init]: Error allocating memory\n");
  return;

}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_2d_forward_y(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from lower left to upper right
 *  for x-aligned rays */
static void sweep_2d_forward_y(RadGridS *pRG, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;
  Real am;

/* Account for ix2 boundary intensities */
  for(i=is-1; i<=ie+1; i++) {
    for(l=0; l<=1; l++)  {
      for(m=0; m<nang; m++) {
        if(am0[m] <= 1.0) {
          imuo[ifr][i][l][m][0] = pRG->Ghstl2i[ifr][ks][i][l][m];
        }
      }}}

  /* sweep forward in x2 */
  for(j=js; j<=je; j++) {

    /* Account for ix1 boundary intensities */
    for(l=0; l<=1; l++)  {
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        if(am0[m] <= 1.0) {
          imuo[ifr][is-1][l][m][1] = imuo[ifr][is-1][l][m][0];
          imuo[ifr][ie+1][l][m][1] = imuo[ifr][ie+1][l][m][0];
          imuo[ifr][is-1][l][m][0] = pRG->Ghstl1i[ifr][ks][j][l][m];
          imuo[ifr][ie+1][l][m][0] = pRG->Ghstr1i[ifr][ks][j][l][m];
        }
      }}

    /* Sweep forward in x1 */
    update_cell_y(pRG,imuo,ifr,ks,j,is,0);
    /* Update intensity at the ix1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->l1imu[ifr][ks][j][0][m] = imuo[ifr][is][0][m][0];
      }
    }
    for(i=is+1; i<=ie; i++)
      update_cell_y(pRG,imuo,ifr,ks,j,i,0);

    /* Update intensity at the ox1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
            pRG->r1imu[ifr][ks][j][0][m] = imuo[ifr][ie][0][m][0];
      }
    }

    /* Sweep backward in x1 */
    update_cell_y(pRG,imuo,ifr,ks,j,ie,1);
    /* Update intensity at the ox1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->r1imu[ifr][ks][j][1][m] = imuo[ifr][ie][1][m][0];
      }
    }
    for(i=ie-1; i>=is; i--)
      update_cell_y(pRG,imuo,ifr,ks,j,i,1);

    /* Update intensity at the ix1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->l1imu[ifr][ks][j][1][m] = imuo[ifr][is][1][m][0];
      }
    }
    if (j == js) {
      /* Update intensity at the ix2 boundary */
      for(i=is; i<=ie; i++) {
        for(l=0; l<=1; l++) {
          for(m=0; m<nang; m++) {
            if(am0[m] <= 1.0) {
              pRG->l2imu[ifr][ks][i][l][m] = imuo[ifr][i][l][m][0];
            }
          }}}
    }
  }
  /* Update intensity at the ox2 boundary */
  for(i=is; i<=ie; i++) {
    for(l=0; l<=1; l++) {
      for(m=0; m<nang; m++) {
        if(am0[m] <= 1.0) {
          pRG->r2imu[ifr][ks][i][l][m] = imuo[ifr][i][l][m][0];
        }
      }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_2d_backward_y(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from upper right to lower left
 *   for y-aligned rays*/
static void sweep_2d_backward_y(RadGridS *pRG, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ox2 boundary intensities */
  for(i=is-1; i<=ie+1; i++) {
    for(l=2; l<=3; l++)  {
      for(m=0; m<nang; m++) {
        if(am0[m] <= 1.0) {
          imuo[ifr][i][l][m][0] = pRG->Ghstr2i[ifr][ks][i][l][m];
        }
      }}}

  /* sweep backward in x2 */
  for(j=je; j>=js; j--) {

    /* Account for ix1 boundary intensities */
    for(l=2; l<=3; l++)  {
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        if(am0[m] <= 1.0) {
          imuo[ifr][is-1][l][m][1] = imuo[ifr][is-1][l][m][0];
          imuo[ifr][ie+1][l][m][1] = imuo[ifr][ie+1][l][m][0];
          imuo[ifr][is-1][l][m][0] = pRG->Ghstl1i[ifr][ks][j][l][m];
          imuo[ifr][ie+1][l][m][0] = pRG->Ghstr1i[ifr][ks][j][l][m];
        }
      }}

    /* Sweep forward in x1 */
    update_cell_y(pRG,imuo,ifr,ks,j,is,2);
    /* Update intensity at the ix1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->l1imu[ifr][ks][j][2][m] = imuo[ifr][is][2][m][0];
      }
    }
    for(i=is+1; i<=ie; i++)
      update_cell_y(pRG,imuo,ifr,ks,j,i,2);

    /* Update intensity at the ox1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->r1imu[ifr][ks][j][2][m] = imuo[ifr][ie][2][m][0];
      }
    }

    /* Sweep backward in x1 */
    update_cell_y(pRG,imuo,ifr,ks,j,ie,3);
    /* Update intensity at the ox1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->r1imu[ifr][ks][j][3][m] = imuo[ifr][ie][3][m][0];
      }
    }
    for(i=ie-1; i>=is; i--)
      update_cell_y(pRG,imuo,ifr,ks,j,i,3);

    /* Update intensity at the ix1 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] <= 1.0) {
        pRG->l1imu[ifr][ks][j][3][m] = imuo[ifr][is][3][m][0];
      }
    }
    if (j == je) {
      /* Update intensity at the ox2 boundary */
      for(i=is; i<=ie; i++) {
        for(l=2; l<=3; l++) {
          for(m=0; m<nang; m++) {
            if(am0[m] <= 1.0) {
              pRG->r2imu[ifr][ks][i][l][m] = imuo[ifr][i][l][m][0];
            }
          }}}
    }
  }

  /* Update intensity at the ix2 boundary */
  for(i=is; i<=ie; i++) {
    for(l=2; l<=3; l++) {
      for(m=0; m<nang; m++) {
        if(am0[m] <= 1.0) {
          pRG->l2imu[ifr][ks][i][l][m] = imuo[ifr][i][l][m][0];
        }
      }}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_2d_forward_x(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from lower left to upper right
 *   for x-aligned rays */
static void sweep_2d_forward_x(RadGridS *pRG, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ix1 boundary intensities */
  for(j=js-1; j<=je+1; j++) {
    for(m=0; m<nang; m++) {
      if(am0[m] > 1.0) {
        imuo[ifr][j][0][m][0] = pRG->Ghstl1i[ifr][ks][j][0][m];
        imuo[ifr][j][2][m][0] = pRG->Ghstl1i[ifr][ks][j][2][m];
      }
    }}

  /* sweep forward in x1 */
  for(i=is; i<=ie; i++) {

    /* Account for ix2 boundary intensities */
    for(m=0; m<nang; m++) {
      /* ix2/ox2 boundary conditions*/
      if(am0[m] > 1.0) {
        imuo[ifr][js-1][0][m][1] = imuo[ifr][js-1][0][m][0];
        imuo[ifr][js-1][2][m][1] = imuo[ifr][js-1][2][m][0];
        imuo[ifr][je+1][0][m][1] = imuo[ifr][je+1][0][m][0];
        imuo[ifr][je+1][2][m][1] = imuo[ifr][je+1][2][m][0];
        imuo[ifr][js-1][0][m][0] = pRG->Ghstl2i[ifr][ks][i][0][m];
        imuo[ifr][js-1][2][m][0] = pRG->Ghstl2i[ifr][ks][i][2][m];
        imuo[ifr][je+1][0][m][0] = pRG->Ghstr2i[ifr][ks][i][0][m];
        imuo[ifr][je+1][2][m][0] = pRG->Ghstr2i[ifr][ks][i][2][m];
      }
    }

    /* Sweep forward in x2 */
    update_cell_x(pRG,imuo,ifr,ks,js,i,0);
    /* Update intensity at the ix2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->l2imu[ifr][ks][i][0][m] = imuo[ifr][js][0][m][0];
      }
    }
    for(j=js+1; j<=je; j++)
      update_cell_x(pRG,imuo,ifr,ks,j,i,0);

    /* Update intensity at the ox2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->r2imu[ifr][ks][i][0][m] = imuo[ifr][je][0][m][0];
      }
    }

    /* Sweep backward in x2 */
    update_cell_x(pRG,imuo,ifr,ks,je,i,2);
    /* Update intensity at the ox2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->r2imu[ifr][ks][i][2][m] = imuo[ifr][je][2][m][0];
      }
    }
    for(j=je-1; j>=js; j--)
      update_cell_x(pRG,imuo,ifr,ks,j,i,2);

    /* Update intensity at the ix2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->l2imu[ifr][ks][i][2][m] = imuo[ifr][js][2][m][0];
      }
    }
    if(i == is) {
      /* Update intensity at the ix1 boundary */
      for(j=js; j<=je; j++) {
        for(m=0; m<nang; m++) {
          if(am0[m] > 1.0) {
            pRG->l1imu[ifr][ks][j][0][m] = imuo[ifr][j][0][m][0];
            pRG->l1imu[ifr][ks][j][2][m] = imuo[ifr][j][2][m][0];
          }
        }}
    }
  }
  /* Update intensity at the ox1 boundary */
  for(j=js; j<=je; j++) {
    for(m=0; m<nang; m++) {
      if(am0[m] > 1.0) {
        pRG->r1imu[ifr][ks][j][0][m] = imuo[ifr][j][0][m][0];
        pRG->r1imu[ifr][ks][j][2][m] = imuo[ifr][j][2][m][0];
      }
    }}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_2d_backward(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from upper right to lower left
 *   for x-aligned rays*/
static void sweep_2d_backward_x(RadGridS *pRG, int ifr)
{
  int i, j, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ox1 boundary intensities */
  for(j=js-1; j<=je+1; j++) {
    for(m=0; m<nang; m++) {
      if(am0[m] > 1.0) {
        imuo[ifr][j][1][m][0] = pRG->Ghstr1i[ifr][ks][j][1][m];
        imuo[ifr][j][3][m][0] = pRG->Ghstr1i[ifr][ks][j][3][m];
      }
    }}

  /* sweep backward in x1 */
  for(i=ie; i>=is; i--) {

    /* Account for ix2 boundary intensities */
    for(m=0; m<nang; m++) {
      /* ix2/ox2 boundary conditions*/
      if(am0[m] > 1.0) {
        imuo[ifr][js-1][1][m][1] = imuo[ifr][js-1][1][m][0];
        imuo[ifr][js-1][3][m][1] = imuo[ifr][js-1][3][m][0];
        imuo[ifr][je+1][1][m][1] = imuo[ifr][je+1][1][m][0];
        imuo[ifr][je+1][3][m][1] = imuo[ifr][je+1][3][m][0];
        imuo[ifr][js-1][1][m][0] = pRG->Ghstl2i[ifr][ks][i][1][m];
        imuo[ifr][js-1][3][m][0] = pRG->Ghstl2i[ifr][ks][i][3][m];
        imuo[ifr][je+1][1][m][0] = pRG->Ghstr2i[ifr][ks][i][1][m];
        imuo[ifr][je+1][3][m][0] = pRG->Ghstr2i[ifr][ks][i][3][m];
      }
    }

    /* Sweep forward in x2 */
    update_cell_x(pRG,imuo,ifr,ks,js,i,1);
    /* Update intensity at the ix2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->l2imu[ifr][ks][i][1][m] = imuo[ifr][js][1][m][0];
      }
    }
    for(j=js+1; j<=je; j++)
      update_cell_x(pRG,imuo,ifr,ks,j,i,1);

    /* Update intensity at the ox2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->r2imu[ifr][ks][i][1][m] = imuo[ifr][je][1][m][0];
      }
    }

    /* Sweep backward in x2 */
    update_cell_x(pRG,imuo,ifr,ks,je,i,3);
    /* Update intensity at the ox2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->r2imu[ifr][ks][i][3][m] = imuo[ifr][je][3][m][0];
      }
    }
    for(j=je-1; j>=js; j--)
      update_cell_x(pRG,imuo,ifr,ks,j,i,3);

    /* Update intensity at the ix2 boundary */
    for(m=0; m<nang; m++)  {
      if(am0[m] > 1.0) {
        pRG->l2imu[ifr][ks][i][3][m] = imuo[ifr][js][3][m][0];
      }
    }
    if(i == ie) {
      /* Update intensity at the ox1 boundary */
      for(j=js; j<=je; j++) {
        for(m=0; m<nang; m++) {
          if(am0[m] > 1.0) {
            pRG->r1imu[ifr][ks][j][1][m] = imuo[ifr][j][1][m][0];
            pRG->r1imu[ifr][ks][j][3][m] = imuo[ifr][j][3][m][0];
          }
        }}
    }
  }

  /* Update intensity at the ix1 boundary */
  for(j=js; j<=je; j++) {
    for(m=0; m<nang; m++) {
      if(am0[m] > 1.0) {
        pRG->l1imu[ifr][ks][j][1][m] = imuo[ifr][j][1][m][0];
        pRG->l1imu[ifr][ks][j][3][m] = imuo[ifr][j][3][m][0];
      }
    }}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_cell_y(RadGridS *pRG, Real *****imuo, int ifr, int k,
 *                                int j, int i, int l)
 *  \brief Update radiation variables in a single cell for y-aligned 
 *         intensities */
static void update_cell_y(RadGridS *pRG, Real *****imuo, int ifr, int k, int j,
                          int i, int l)
{

  int im, ip, jm, jp;
  int m, nf = pRG->nf, nang = pRG->nang;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real am, am1, bm, bm1;
  Real w0, w1, w2;
  Real maxint, minint;
  Real dx = pRG->dx1, dy = pRG->dx2;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;

/* initialize stencil base on quadrant*/
  if(l == 0) {
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 1) {
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 2) {
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else {
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  }

  for(m=0; m<nang; m++) {
    chi1 = pRG->R[ifr][k][j][i].chi;
/* --------- Interpolate intensity and source functions at endpoints --------
 * --------- of characteristics                                      -------- */
    am = am0[m];
    if (am <= 1.0) {
      am1 = 1.0 - am;
      /* Use linear interpolation for source functions */
      S0 = am  * pRG->R[ifr][k][jm][im].S +
           am1 * pRG->R[ifr][k][jm][i ].S;
      S2 = am  * pRG->R[ifr][k][jp][ip].S +
           am1 * pRG->R[ifr][k][jp][i ].S;
        /* Use quadratic interpolation for intensity */
      w0 = 0.5 * am * (1.0 + am);
      w1 = am1 * (1.0 + am);
      w2 = -0.5 * am * am1;
      imu0 = w0 * imuo[ifr][im][l][m][1] + w1 * imuo[ifr][i][l][m][0] +
             w2 * imuo[ifr][ip][l][m][0];
      maxint = MAX(imuo[ifr][im][l][m][1],imuo[ifr][i][l][m][0]);
      minint = MIN(imuo[ifr][im][l][m][1],imuo[ifr][i][l][m][0]);
      if(imu0 > maxint) imu0 = maxint;
      if(imu0 < minint) imu0 = minint;
/* -------  compute intensity at grid center and add to mean intensity ------ */

      chi0 = am  * pRG->R[ifr][k][jm][im].chi +
             am1 * pRG->R[ifr][k][jm][i ].chi;
      chi2 = am  * pRG->R[ifr][k][jp][ip].chi +
             am1 * pRG->R[ifr][k][jp][i ].chi;
      /*dtaum = 0.5 * (chi0 + chi1);
        dtaup = 0.5 * (chi2 + chi1); */
      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dy / pRG->mu[0][m][1];
      dtaup *= dy / pRG->mu[0][m][1];
      interp_quad_source_slope_lim(dtaum, dtaup, &edtau, &a0, &a1, &a2,
                                   S0, pRG->R[ifr][k][j][i].S, S2);
      imu = a0 * S0 + a1 * pRG->R[ifr][k][j][i].S + a2 * S2 + edtau * imu0;
      lamstr[ifr][j][i] += pRG->wmu[m] * a1;
/* Add to radiation moments and save for next iteration */
      wimu = pRG->wmu[m] * imu;
      pRG->R[ifr][k][j][i].J += wimu;
      pRG->R[ifr][k][j][i].H[0] += pRG->mu[l][m][0] * wimu;
      pRG->R[ifr][k][j][i].H[1] += pRG->mu[l][m][1] * wimu;
      pRG->R[ifr][k][j][i].K[0] += pRG->mu[l][m][0] * pRG->mu[l][m][0] * wimu;
      pRG->R[ifr][k][j][i].K[1] += pRG->mu[l][m][0] * pRG->mu[l][m][1] * wimu;
      pRG->R[ifr][k][j][i].K[2] += pRG->mu[l][m][1] * pRG->mu[l][m][1] * wimu;
      /* Update intensity workspace */
      imuo[ifr][i][l][m][1] = imuo[ifr][i][l][m][0];
      imuo[ifr][i][l][m][0] = imu;
    }
  }
  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_cell_x(RadGridS *pRG, Real *****imuo, int ifr, int k,
 *                                int j, int i, int l)
 *  \brief Update radiation variables in a single cell for x-aligned
           intensities */
static void update_cell_x(RadGridS *pRG, Real *****imuo, int ifr, int k, int j,
                          int i, int l)
{

  int im, ip, jm, jp;
  int m, nf = pRG->nf, nang = pRG->nang;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real am, am1, bm, bm1;
  Real w0, w1, w2;
  Real maxint, minint;
  Real dx = pRG->dx1, dy = pRG->dx2;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;

/* initialize stencil base on quadrant*/
  if(l == 0) {
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 1) {
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 2) {
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else {
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  }

  for(m=0; m<nang; m++) {
    chi1 = pRG->R[ifr][k][j][i].chi;
/* -------- Interpolate intensity and source functions at endpoints  --------
 * -------- of characteristics                                       -------- */
    am = am0[m];
    if (am > 1.0) {
      bm = 1.0 / am;
      bm1 = 1.0 - bm;
      /* Use linear interpolation for source functions */
      S0 = bm  * pRG->R[ifr][k][jm][im].S +
           bm1 * pRG->R[ifr][k][j ][im].S;
      S2 = bm  * pRG->R[ifr][k][jp][ip].S +
           bm1 * pRG->R[ifr][k][j ][ip].S;
      /* Use linear interpolation for intensity */
      w0 = 0.5 * bm * (1.0 + bm);
      w1 = bm1 * (1.0 + bm);
      w2 = -0.5 * bm * bm1;
      imu0 = w0 * imuo[ifr][jm][l][m][1] + w1 * imuo[ifr][j][l][m][0] +
             w2 * imuo[ifr][jp][l][m][0];
      maxint = MAX(imuo[ifr][jm][l][m][1],imuo[ifr][j][l][m][0]);
      minint = MIN(imuo[ifr][jm][l][m][1],imuo[ifr][j][l][m][0]);
      if(imu0 > maxint) imu0 = maxint;
      if(imu0 < minint) imu0 = minint;
/* -------  compute intensity at grid center and add to mean intensity ------ */

      chi0 = bm  * pRG->R[ifr][k][jm][im].chi +
             bm1 * pRG->R[ifr][k][j ][im].chi;
      chi2 = bm  * pRG->R[ifr][k][jp][ip].chi +
             bm1 * pRG->R[ifr][k][j ][ip].chi;
      /*dtaum = 0.5 * (chi0 + chi1);
        dtaup = 0.5 * (chi2 + chi1); */
      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dx / pRG->mu[0][m][0];
      dtaup *= dx / pRG->mu[0][m][0];
      interp_quad_source_slope_lim(dtaum, dtaup, &edtau, &a0, &a1, &a2,
                                   S0, pRG->R[ifr][k][j][i].S, S2);
      imu = a0 * S0 + a1 * pRG->R[ifr][k][j][i].S + a2 * S2 + edtau * imu0;
      lamstr[ifr][j][i] += pRG->wmu[m] * a1;

/* Add to radiation moments and save for next iteration */
      wimu = pRG->wmu[m] * imu;
      pRG->R[ifr][k][j][i].J += wimu;
      pRG->R[ifr][k][j][i].H[0] += pRG->mu[l][m][0] * wimu;
      pRG->R[ifr][k][j][i].H[1] += pRG->mu[l][m][1] * wimu;
      pRG->R[ifr][k][j][i].K[0] += pRG->mu[l][m][0] * pRG->mu[l][m][0] * wimu;
      pRG->R[ifr][k][j][i].K[1] += pRG->mu[l][m][0] * pRG->mu[l][m][1] * wimu;
      pRG->R[ifr][k][j][i].K[2] += pRG->mu[l][m][1] * pRG->mu[l][m][1] * wimu;
 /* Update intensity workspace */
      imuo[ifr][j][l][m][1] = imuo[ifr][j][l][m][0];
      imuo[ifr][j][l][m][0] = imu;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_sfunc(RadS *R, Real *dSr, Real lamstr)
 *  \brief Jacobi update of source function with new mean intensity */
static void update_sfunc(RadS *R, Real *dSr, Real lamstr)
{
  Real Snew, dS;

  Snew = (1.0 - R->eps) * R->J + R->eps * R->B + R->Snt;
  dS = (Snew - R->S) / (1.0 - (1.0 - R->eps) * lamstr);
  if (R->S > 0.0) (*dSr) = fabs(dS / R->S);
  R->S += dS;

  return;
}


#endif /* JACOBI  && QUADTRATIC_INTENSITY */
#endif /* RADIATION_TRANSFER */
