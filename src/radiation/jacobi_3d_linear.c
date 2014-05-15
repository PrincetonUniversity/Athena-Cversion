
#include "../copyright.h"
/*==============================================================================
 * FILE: jacobi_3d_linear.c
 *
 * PURPOSE: Solves a single iteration of the formal solution of radiative
 *          transfer on a 3D grid using jacobi's method.  The basic algorithm
 *          is described in Trujillo Bueno and Fabiani Benedicho, ApJ, 455, 646.
 *          Uses linear interpolation to compuet intensities at edges.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   formal_solution_3d.c()
 *   formal_solution_3d_destruct()
 *   formal_solution_3d_init()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"


#ifdef RADIATION_TRANSFER
#ifdef JACOBI


/* Working arrays used in formal solution */
static Real ****lamstr = NULL;
static Real ******imuo = NULL;
static int *face = NULL;
static Real **coeff = NULL;
static Real ****Jold = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   update_cell() - update variables in a single grid cell
 *   update_sfunc() - updates source function after compute of mean intensity
 *   sweep_3d_forward() - sweep grid from lower left to upper right
 *   sweep_3d_backward() - sweep grid from upper right to lower left
 *============================================================================*/

static void update_cell(RadGridS *pRG, Real ******imuo, int ifr, int k, int j,
                        int i, int l);
static void update_sfunc(RadS *R, Real *dSr, Real lam);
static void sweep_3d_forward(RadGridS *pRG, int ifr);
static void sweep_3d_backward(RadGridS *pRG, int ifr);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_3d(RadGridS *pRG, Real *dSrmax, int ifr)
 *  \brief formal solution for single freq. in 3D using Jacobi method
 *  with linear interpolation of intensities */
void formal_solution_3d(RadGridS *pRG, Real *dSrmax, int ifr)
{
  int i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf;
  int ismx, jsmx, ksmx;
  Real dSr, dJ, dJmax;
  Real gdSrmax;
  Real lx, ly, lz, lmin;
  Real dx, dy, dz;
  Real am, bm;

#ifdef QUADRATIC_INTENSITY
  ath_error("[jacob_3d_linear.c]: quadratic intensity not currently supported "
            "on 3D domains.\n");
#endif

  if (radt_mode == 2) {
/* Initialize face and coeff arrays on each call to formal soution when
   both passive and feedback based radiation transfer schemes are enabled */
    for(i=0; i<pRG->nang; i++) {
      lx = dx / pRG->mu[0][i][0];
      ly = dy / pRG->mu[0][i][1];
      lz = dz / pRG->mu[0][i][2];
      lmin = MIN(MIN(lx,ly),lz);
      if (lz == lmin) {
        face[i]=2;
        am = lmin/lx; bm = lmin/ly;
      } else if (ly == lmin) {
        face[i]=1;
        am = lmin/lx; bm = lmin/lz;
      } else {
        face[i]=0;
        am = lmin/ly; bm = lmin/lz;
      }
      coeff[i][0] = (1.0 - am)*(1.0 - bm);
      coeff[i][1] = (1.0 - am)*       bm;
      coeff[i][2] =        am *       bm;
      coeff[i][3] =        am *(1.0 - bm);
    }
  }

/* if LTE then store J values from previous iteration */
  if(lte != 0) {
    for(k=ks; k<=ke; k++)
      for(j=js; j<=je; j++)
        for(i=is; i<=ie; i++)
          Jold[k][j][i][ifr] = pRG->R[ifr][k][j][i].J;
  }

/* initialize mean intensities at all depths to zero */
  for(k=ks; k<=ke; k++)
    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++) {
        pRG->R[ifr][k][j][i].J = 0.0;
        pRG->R[ifr][k][j][i].H[0] = 0.0;
        pRG->R[ifr][k][j][i].H[1] = 0.0;
        pRG->R[ifr][k][j][i].H[2] = 0.0;
        pRG->R[ifr][k][j][i].K[0] = 0.0;
        pRG->R[ifr][k][j][i].K[1] = 0.0;
        pRG->R[ifr][k][j][i].K[2] = 0.0;
        pRG->R[ifr][k][j][i].K[3] = 0.0;
        pRG->R[ifr][k][j][i].K[4] = 0.0;
        pRG->R[ifr][k][j][i].K[5] = 0.0;
        lamstr[ifr][k][j][i] = 0.0;
      }

/* Compute formal solution and for all rays in each gridzone and
 * update boundary emission*/
  sweep_3d_forward(pRG,ifr);

  sweep_3d_backward(pRG,ifr);

  if(lte == 0) {
/* Update source function */
    (*dSrmax) = 0.0;
    for(k=ks; k<=ke; k++)
      for(j=js; j<=je; j++)
        for(i=is; i<=ie; i++) {
          update_sfunc(&(pRG->R[ifr][k][j][i]),&dSr,lamstr[ifr][k][j][i]);
          if( dSr > (*dSrmax)) {
            (*dSrmax) = dSr; ismx=i; jsmx=j; ksmx=k;
          }
        }
  } else {
/* Use delta J / J as convergence criterion */
    (*dSrmax) = 0.0;
    dJmax = 0.0;
    for(k=ks; k<=ke; k++)
      for(j=js; j<=je; j++)
        for(i=is; i<=ie; i++) {
          dJ = fabs(pRG->R[ifr][k][j][i].J - Jold[k][j][i][ifr]);
          if(dJ > dJmax) dJmax = dJ;
          if (Jold[k][j][i][ifr] > 0.0)
            dSr = dJ / Jold[k][j][i][ifr];
          else
            dSr = 0;
          if( dSr > (*dSrmax)) {
            (*dSrmax) = dSr; ismx=i; jsmx=j; ksmx=k;
          }
        }
    if(((*dSrmax) == 0.0) && (dJmax > 0.0)) (*dSrmax) = 1.0;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_3d_destruct(void)
 *  \brief free temporary working arrays */
void formal_solution_3d_destruct(void)
{

  if (lamstr != NULL) free_4d_array(lamstr);
  if (imuo   != NULL) free_6d_array(imuo);
  if (face   != NULL) free_1d_array(face);
  if (coeff  != NULL) free_2d_array(coeff);
  if (Jold   != NULL) free_4d_array(Jold);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_3d_init(DomainS *pD)
 *  \brief Allocate memory for working arrays */
void formal_solution_3d_init(DomainS *pD)
{

  RadGridS *pRG, *pROG;
  int nx1, nx2, nx3;
  Real dx, dy, dz;
  int nf, nang;
  int i;
  Real am, bm, lx, ly, lz, lmin;

  if (radt_mode == 0) { /* integration only */
    pRG = pD->RadGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1]; nx3 = pRG->Nx[2];
    dx = pRG->dx1; dy = pRG->dx2; dz = pRG->dx3;
    nf = pRG->nf;
    nang = pRG->nang;
  } else if (radt_mode == 1) { /* output only */
    pRG = pD->RadOutGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1]; nx3 = pRG->Nx[2];
    dx = pRG->dx1; dy = pRG->dx2; dz = pRG->dx3;
    nf = pRG->nf;
    nang = pRG->nang;
  } else if (radt_mode == 2) { /* integration and output */
    pRG = pD->RadGrid;
    pROG = pD->RadOutGrid;
    nx1 = pRG->Nx[0]; nx2 = pRG->Nx[1]; nx3 = pRG->Nx[2];
    dx = pRG->dx1; dy = pRG->dx2; dz = pRG->dx3;
    nf = MAX(pRG->nf,pROG->nf);
    nang = MAX(pRG->nang,pROG->nang);
  }

  if ((lamstr = (Real ****)calloc_4d_array(nf,nx3+2,nx2+2,nx1+2,
                                           sizeof(Real))) == NULL)
   goto on_error;

  if ((imuo = (Real ******)calloc_6d_array(nf,nx2+2,nx1+2,8,nang,2,
                                           sizeof(Real))) == NULL)
    goto on_error;

  if ((face = (int *)calloc_1d_array(nang,sizeof(int))) == NULL)
    goto on_error;

  if ((coeff = (Real **)calloc_2d_array(nang,4,sizeof(Real))) == NULL)
    goto on_error;

  if ( (radt_mode == 0) ||  (radt_mode == 1) ) {
/* Initialize face and coeff arrays only once if radiation transfer is only
   used for output or integration but not both */

    for(i=0; i<nang; i++) {
      lx = dx / pRG->mu[0][i][0];
      ly = dy / pRG->mu[0][i][1];
      lz = dz / pRG->mu[0][i][2];
      lmin = MIN(MIN(lx,ly),lz);
      if (lz == lmin) {
        face[i]=2;
        am = lmin/lx; bm = lmin/ly;
      } else if (ly == lmin) {
        face[i]=1;
        am = lmin/lx; bm = lmin/lz;
      } else {
        face[i]=0;
        am = lmin/ly; bm = lmin/lz;
      }
      coeff[i][0] = (1.0 - am)*(1.0 - bm);
      coeff[i][1] = (1.0 - am)*       bm;
      coeff[i][2] =        am *       bm;
      coeff[i][3] =        am *(1.0 - bm);
    }
  }

  if (lte != 0)
    if ((Jold = (Real ****)calloc_4d_array(nx3+2,nx2+2,nx1+2,nf,
                                           sizeof(Real))) == NULL)
      goto on_error;

  return;

  on_error:
  formal_solution_3d_destruct();
  ath_error("[formal_solution__3d_init]: Error allocating memory\n");
  return;

}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_3d_forward(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from lower left to upper right */
static void sweep_3d_forward(RadGridS *pRG, int ifr)
{
  int i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ix3 boundary intensities */
  for(j=js-1; j<=je+1; j++) {
    for(i=is-1; i<=ie+1; i++) {
      for(l=0; l<=3; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][j][i][l][m][0] = pRG->Ghstl3i[ifr][j][i][l][m];
        }}}}

  /* sweep forward in x3 */
  for(k=ks; k<=ke; k++) {

    /* Account for ix2 boundary intensities.  Note that this uses
     * l2imu to initialize imuo on edge */
    for(i=is-1; i<=ie+1; i++) {
      for(l=0; l<=1; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][js-1][i][l][m][1] = imuo[ifr][js-1][i][l][m][0];
          imuo[ifr][js-1][i][l][m][0] = pRG->Ghstl2i[ifr][k][i][l][m];
        }}}

    /* Sweep forward in x2 */
    for(j=js; j<=je; j++) {

      /* Account for ix1 boundary intensities */
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        imuo[ifr][j][is-1][0][m][1] = imuo[ifr][j][is-1][0][m][0];
        imuo[ifr][j][ie+1][1][m][1] = imuo[ifr][j][ie+1][1][m][0];
        imuo[ifr][j][is-1][0][m][0] = pRG->Ghstl1i[ifr][k][j][0][m];
        imuo[ifr][j][ie+1][1][m][0] = pRG->Ghstr1i[ifr][k][j][1][m];
      }

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++)
        update_cell(pRG,imuo,ifr,k,j,i,0);

      /* Update intensity at the ox1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->r1imu[ifr][k][j][0][m] = imuo[ifr][j][ie][0][m][0];
      }

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--)
        update_cell(pRG,imuo,ifr,k,j,i,1);

      /* Update intensity at the ix1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->l1imu[ifr][k][j][1][m] = imuo[ifr][j][is][1][m][0];
      }
    }

    /* Update intensity at the ox2 boundary */
    for(i=is; i<=ie; i++) {
      for(l=0; l<=1; l++) {
        for(m=0; m<nang; m++) {
          pRG->r2imu[ifr][k][i][l][m] = imuo[ifr][je][i][l][m][0];
        }}}

/* ----------------  Start of reverse sweep ---------------------- */

    /* Account for ox2 boundary intensities */
    for(i=is-1; i<=ie+1; i++) {
      for(l=2; l<=3; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][je+1][i][l][m][1] = imuo[ifr][je+1][i][l][m][0];
          imuo[ifr][je+1][i][l][m][0] = pRG->Ghstr2i[ifr][k][i][l][m];
        }}}

    /* sweep backward in x2 */
    for(j=je; j>=js; j--) {

      /* Account for ix1 boundary intensities */
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        imuo[ifr][j][is-1][2][m][1] = imuo[ifr][j][is-1][2][m][0];
        imuo[ifr][j][ie+1][3][m][1] = imuo[ifr][j][ie+1][3][m][0];
        imuo[ifr][j][is-1][2][m][0] = pRG->Ghstl1i[ifr][k][j][2][m];
        imuo[ifr][j][ie+1][3][m][0] = pRG->Ghstr1i[ifr][k][j][3][m];
      }

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++)
        update_cell(pRG,imuo,ifr,k,j,i,2);

      /* Update intensity at the ox1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->r1imu[ifr][k][j][2][m] = imuo[ifr][j][ie][2][m][0];
      }

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--)
        update_cell(pRG,imuo,ifr,k,j,i,3);

      /* Update intensity at the ix1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->l1imu[ifr][k][j][3][m] = imuo[ifr][j][is][3][m][0];
      }
    }

    /* Update intensity at the ix2 boundary */
    for(i=is; i<=ie; i++) {
      for(l=2; l<=3; l++) {
        for(m=0; m<nang; m++) {
          pRG->l2imu[ifr][k][i][l][m] = imuo[ifr][js][i][l][m][0];
        }}}
  }

   /* Update intensity at the ox3 boundary */
  for(j=js; j<=je; j++) {
    for(i=is; i<=ie; i++) {
      for(l=0; l<=3; l++) {
        for(m=0; m<nang; m++) {
          pRG->r3imu[ifr][j][i][l][m] = imuo[ifr][j][i][l][m][0];
        }}}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void sweep_3d_backward(RadGridS *pRG, int ifr)
 *  \brief Perform Jacobi sweep from upper right to lower left*/
static void sweep_3d_backward(RadGridS *pRG, int ifr)
{
  int i, j, k, l, m;
  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js, je = pRG->je;
  int ks = pRG->ks, ke = pRG->ke;
  int nf = pRG->nf, nang = pRG->nang;

/* Account for ix3 boundary intensities */
  for(j=js-1; j<=je+1; j++) {
    for(i=is-1; i<=ie+1; i++) {
      for(l=4; l<=7; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][j][i][l][m][0] = pRG->Ghstr3i[ifr][j][i][l][m];
        }}}}

  /* sweep forward in x3 */
  for(k=ke; k>=ks; k--) {

    /* Account for ix2 boundary intensities.  Note that this uses
     * l2imu to initialize imuo on edge */
    for(i=is-1; i<=ie+1; i++) {
      for(l=4; l<=5; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][js-1][i][l][m][1] = imuo[ifr][js-1][i][l][m][0];
          imuo[ifr][js-1][i][l][m][0] = pRG->Ghstl2i[ifr][k][i][l][m];
        }}}

    /* Sweep forward in x2 */
    for(j=js; j<=je; j++) {

      /* Account for ix1 boundary intensities */
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        imuo[ifr][j][is-1][4][m][1] = imuo[ifr][j][is-1][4][m][0];
        imuo[ifr][j][ie+1][5][m][1] = imuo[ifr][j][ie+1][5][m][0];
        imuo[ifr][j][is-1][4][m][0] = pRG->Ghstl1i[ifr][k][j][4][m];
        imuo[ifr][j][ie+1][5][m][0] = pRG->Ghstr1i[ifr][k][j][5][m];
      }

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++)
        update_cell(pRG,imuo,ifr,k,j,i,4);

      /* Update intensity at the ox1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->r1imu[ifr][k][j][4][m] = imuo[ifr][j][ie][4][m][0];
      }

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--)
        update_cell(pRG,imuo,ifr,k,j,i,5);

      /* Update intensity at the ix1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->l1imu[ifr][k][j][5][m] = imuo[ifr][j][is][5][m][0];
      }
    }

    /* Update intensity at the ox2 boundary */
    for(i=is; i<=ie; i++) {
      for(l=4; l<=5; l++) {
        for(m=0; m<nang; m++) {
          pRG->r2imu[ifr][k][i][l][m] = imuo[ifr][je][i][l][m][0];
        }}}

/* ----------------  Start of reverse sweep ---------------------- */

    /* Account for ox2 boundary intensities */
    for(i=is-1; i<=ie+1; i++) {
      for(l=6; l<=7; l++)  {
        for(m=0; m<nang; m++) {
          imuo[ifr][je+1][i][l][m][1] = imuo[ifr][je+1][i][l][m][0];
          imuo[ifr][je+1][i][l][m][0] = pRG->Ghstr2i[ifr][k][i][l][m];
        }}}

    /* sweep backward in x2 */
    for(j=je; j>=js; j--) {

      /* Account for ix1 boundary intensities */
      for(m=0; m<nang; m++) {
        /* ix1/ox1 boundary conditions*/
        imuo[ifr][j][is-1][6][m][1] = imuo[ifr][j][is-1][6][m][0];
        imuo[ifr][j][ie+1][7][m][1] = imuo[ifr][j][ie+1][7][m][0];
        imuo[ifr][j][is-1][6][m][0] = pRG->Ghstl1i[ifr][k][j][6][m];
        imuo[ifr][j][ie+1][7][m][0] = pRG->Ghstr1i[ifr][k][j][7][m];
      }

      /* Sweep forward in x1 */
      for(i=is; i<=ie; i++)
        update_cell(pRG,imuo,ifr,k,j,i,6);

      /* Update intensity at the ox1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->r1imu[ifr][k][j][6][m] = imuo[ifr][j][ie][6][m][0];
      }

      /* Sweep backward in x1 */
      for(i=ie; i>=is; i--)
        update_cell(pRG,imuo,ifr,k,j,i,7);

      /* Update intensity at the ix1 boundary */
      for(m=0; m<nang; m++)  {
        pRG->l1imu[ifr][k][j][7][m] = imuo[ifr][j][is][7][m][0];
      }
    }

    /* Update intensity at the ix2 boundary */
    for(i=is; i<=ie; i++) {
      for(l=6; l<=7; l++) {
        for(m=0; m<nang; m++) {
          pRG->l2imu[ifr][k][i][l][m] = imuo[ifr][js][i][l][m][0];
        }}}
  }

   /* Update intensity at the ox3 boundary */
  for(j=js; j<=je; j++) {
    for(i=is; i<=ie; i++) {
      for(l=4; l<=7; l++) {
        for(m=0; m<nang; m++) {
          pRG->l3imu[ifr][j][i][l][m] = imuo[ifr][j][i][l][m][0];
        }}}}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_cell(RadGridS *pRG, Real ******imuo, int ifr, int k,
 *                              int j, int i, int l)
 *  \brief Update radiation variables in a single cell */
static void update_cell(RadGridS *pRG, Real ******imuo, int ifr, int k, int j, int i, int l)
{

  int im, ip, jm, jp, km, kp;
  int m, nf = pRG->nf, nang = pRG->nang;
  Real imu, imu0, wimu;
  Real S0, S2;
  Real w0, w1, w2;
  Real dx = pRG->dx1, dy = pRG->dx2, dz = pRG->dx3;
  Real chi0, chi1, chi2, dtaum, dtaup;
  Real edtau, a0, a1, a2;

/* initialize stencil base on quadrant*/
  if(l == 0) {
    kp = k + 1;  km = k - 1;
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 1) {
    kp = k + 1;  km = k - 1;
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 2) {
    kp = k + 1;  km = k - 1;
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 3) {
    kp = k + 1;  km = k - 1;
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 4) {
    kp = k - 1;  km = k + 1;
    jp = j + 1;  jm = j - 1;
    ip = i + 1;  im = i - 1;
  } else if (l == 5) {
    kp = k - 1;  km = k + 1;
    jp = j + 1;  jm = j - 1;
    ip = i - 1;  im = i + 1;
  } else if (l == 6) {
    kp = k - 1;  km = k + 1;
    jp = j - 1;  jm = j + 1;
    ip = i + 1;  im = i - 1;
  } else {
    kp = k - 1;  km = k + 1;
    jp = j - 1;  jm = j + 1;
    ip = i - 1;  im = i + 1;
  }


  chi1 = pRG->R[ifr][k][j][i].chi;
  for(m=0; m<nang; m++) {
/* --------- Interpolate intensity and source functions at endpoints --------
 * --------- of characteristics                                      -------- */
    if (face[m] == 0) {
      /* interpolation in x2-x3 plane */

      S0    = coeff[m][0] * pRG->R[ifr][k ][j ][im].S;
      chi0  = coeff[m][0] * pRG->R[ifr][k ][j ][im].chi;
      S0   += coeff[m][1] * pRG->R[ifr][km][j ][im].S;
      chi0 += coeff[m][1] * pRG->R[ifr][km][j ][im].chi;
      S0   += coeff[m][2] * pRG->R[ifr][km][jm][im].S;
      chi0 += coeff[m][2] * pRG->R[ifr][km][jm][im].chi;
      S0   += coeff[m][3] * pRG->R[ifr][k ][jm][im].S;
      chi0 += coeff[m][3] * pRG->R[ifr][k ][jm][im].chi;

      S2    = coeff[m][0] * pRG->R[ifr][k ][j ][ip].S;
      chi2  = coeff[m][0] * pRG->R[ifr][k ][j ][ip].chi;
      S2   += coeff[m][1] * pRG->R[ifr][kp][j ][ip].S;
      chi2 += coeff[m][1] * pRG->R[ifr][kp][j ][ip].chi;
      S2   += coeff[m][2] * pRG->R[ifr][kp][jp][ip].S;
      chi2 += coeff[m][2] * pRG->R[ifr][kp][jp][ip].chi;
      S2   += coeff[m][3] * pRG->R[ifr][k ][jp][ip].S;
      chi2 += coeff[m][3] * pRG->R[ifr][k ][jp][ip].chi;

      imu0  = coeff[m][0] * imuo[ifr][j ][im][l][m][0] +
              coeff[m][1] * imuo[ifr][j ][im][l][m][1] +
              coeff[m][2] * imuo[ifr][jm][im][l][m][1] +
              coeff[m][3] * imuo[ifr][jm][im][l][m][0];

      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dx / pRG->mu[0][m][0];
      dtaup *= dx / pRG->mu[0][m][0];
    } else if(face[m] == 1) {
      /* interpolation in x1-x3 plane */
      S0    = coeff[m][0] * pRG->R[ifr][k ][jm][i ].S;
      chi0  = coeff[m][0] * pRG->R[ifr][k ][jm][i ].chi;
      S0   += coeff[m][1] * pRG->R[ifr][km][jm][i ].S;
      chi0 += coeff[m][1] * pRG->R[ifr][km][jm][i ].chi;
      S0   += coeff[m][2] * pRG->R[ifr][km][jm][im].S;
      chi0 += coeff[m][2] * pRG->R[ifr][km][jm][im].chi;
      S0   += coeff[m][3] * pRG->R[ifr][k ][jm][im].S;
      chi0 += coeff[m][3] * pRG->R[ifr][k ][jm][im].chi;

      S2    = coeff[m][0] * pRG->R[ifr][k ][jp][i ].S;
      chi2  = coeff[m][0] * pRG->R[ifr][k ][jp][i ].chi;
      S2   += coeff[m][1] * pRG->R[ifr][kp][jp][i ].S;
      chi2 += coeff[m][1] * pRG->R[ifr][kp][jp][i ].chi;
      S2   += coeff[m][2] * pRG->R[ifr][kp][jp][ip].S;
      chi2 += coeff[m][2] * pRG->R[ifr][kp][jp][ip].chi;
      S2   += coeff[m][3] * pRG->R[ifr][k ][jp][ip].S;
      chi2 += coeff[m][3] * pRG->R[ifr][k ][jp][ip].chi;

      imu0  = coeff[m][0] * imuo[ifr][jm][i ][l][m][0] +
              coeff[m][1] * imuo[ifr][jm][i ][l][m][1] +
              coeff[m][2] * imuo[ifr][jm][im][l][m][1] +
              coeff[m][3] * imuo[ifr][jm][im][l][m][0];

      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dy / pRG->mu[0][m][1];
      dtaup *= dy / pRG->mu[0][m][1];
    } else  {
      /* interpolation in x1-x2 plane */
      S0    = coeff[m][0] * pRG->R[ifr][km][j ][i ].S;
      chi0  = coeff[m][0] * pRG->R[ifr][km][j ][i ].chi;
      S0   += coeff[m][1] * pRG->R[ifr][km][jm][i ].S;
      chi0 += coeff[m][1] * pRG->R[ifr][km][jm][i ].chi;
      S0   += coeff[m][2] * pRG->R[ifr][km][jm][im].S;
      chi0 += coeff[m][2] * pRG->R[ifr][km][jm][im].chi;
      S0   += coeff[m][3] * pRG->R[ifr][km][j ][im].S;
      chi0 += coeff[m][3] * pRG->R[ifr][km][j ][im].chi;

      S2    = coeff[m][0] * pRG->R[ifr][kp][j ][i ].S;
      chi2  = coeff[m][0] * pRG->R[ifr][kp][j ][i ].chi;
      S2   += coeff[m][1] * pRG->R[ifr][kp][jp][i ].S;
      chi2 += coeff[m][1] * pRG->R[ifr][kp][jp][i ].chi;
      S2   += coeff[m][2] * pRG->R[ifr][kp][jp][ip].S;
      chi2 += coeff[m][2] * pRG->R[ifr][kp][jp][ip].chi;
      S2   += coeff[m][3] * pRG->R[ifr][kp][j ][ip].S;
      chi2 += coeff[m][3] * pRG->R[ifr][kp][j ][ip].chi;

      imu0  = coeff[m][0] * imuo[ifr][j ][i ][l][m][0] +
              coeff[m][1] * imuo[ifr][jm][i ][l][m][1] +
              coeff[m][2] * imuo[ifr][jm][im][l][m][1] +
              coeff[m][3] * imuo[ifr][j ][im][l][m][1];

      interp_quad_chi(chi0,chi1,chi2,&dtaum);
      interp_quad_chi(chi2,chi1,chi0,&dtaup);
      dtaum *= dz / pRG->mu[0][m][2];
      dtaup *= dz / pRG->mu[0][m][2];
    }
/* -------  compute intensity at grid center and add to mean intensity ------ */
    interp_quad_source_slope_lim(dtaum, dtaup, &edtau, &a0, &a1, &a2,
                                 S0, pRG->R[ifr][k][j][i].S, S2);
    imu = a0 * S0 + a1 * pRG->R[ifr][k][j][i].S + a2 * S2 + edtau * imu0;
    /*  if (imu != imu) {
      printf("%d %d %d %d %g %g %g %g\n",k,j,i,m,imu0,dtaum,dtaup,pRG->R[ifr][k][j][i].S);
      printf("S %g %g\n",S0,S2);
      printf("chi %g %g %g\n",chi0,chi1,chi2);
      printf("mu %d %g %g %g\n",face[m],pRG->mu[0][m][0],pRG->mu[0][m][1],pRG->mu[0][m][2]);
      printf("coeff %g %g %g %g\n",coeff[m][0],coeff[m][1],coeff[m][2],coeff[m][3]);
      if (face[m] == 2) {
        printf("%g %g %g %g\n", pRG->R[ifr][km][j ][i ].S,pRG->R[ifr][km][jm][i ].S,
                pRG->R[ifr][km][jm][im].S,pRG->R[ifr][km][j ][im].S);
      }
      }*/
    lamstr[ifr][k][j][i] += pRG->wmu[m] * a1;
/* Add to radiation moments and save for next iteration */
    wimu = pRG->wmu[m] * imu;
    pRG->R[ifr][k][j][i].J += wimu;
    pRG->R[ifr][k][j][i].H[0] += pRG->mu[l][m][0] * wimu;
    pRG->R[ifr][k][j][i].H[1] += pRG->mu[l][m][1] * wimu;
    pRG->R[ifr][k][j][i].H[2] += pRG->mu[l][m][2] * wimu;
    pRG->R[ifr][k][j][i].K[0] += pRG->mu[l][m][0] * pRG->mu[l][m][0] * wimu;
    pRG->R[ifr][k][j][i].K[1] += pRG->mu[l][m][0] * pRG->mu[l][m][1] * wimu;
    pRG->R[ifr][k][j][i].K[2] += pRG->mu[l][m][1] * pRG->mu[l][m][1] * wimu;
    pRG->R[ifr][k][j][i].K[3] += pRG->mu[l][m][0] * pRG->mu[l][m][2] * wimu;
    pRG->R[ifr][k][j][i].K[4] += pRG->mu[l][m][1] * pRG->mu[l][m][2] * wimu;
    pRG->R[ifr][k][j][i].K[5] += pRG->mu[l][m][2] * pRG->mu[l][m][2] * wimu;
/* Update intensity workspace */
    imuo[ifr][j][i][l][m][1] = imuo[ifr][j][i][l][m][0];
    imuo[ifr][j][i][l][m][0] = imu;
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn static void update_sfunc(RadS *R, Real *dSr, Real lamstr)
 *  \brief Jacobi update of source function with new mean intensity */
static void update_sfunc(RadS *R, Real *dSr, Real lam)
{
  Real Snew, dS;

  Snew = (1.0 - R->eps) * R->J + R->eps * R->B + R->Snt;
  dS = (Snew - R->S) / (1.0 - (1.0 - R->eps) * lam);
  if (R->S > 0.0) (*dSr) = fabs(dS / R->S);
  R->S += dS;

  return;
}

#endif /* JACOBI */
#endif /* RADIATION_TRANSFER */
