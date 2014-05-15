#include "../copyright.h"
/*==============================================================================
 * FILE: angles.c
 *
 * PURPOSE:  Contains all routiens for specifing and intializing angular
 *           grid and quadrature.  Called by init_radiation().
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   init_angles()
 *============================================================================*/

#include <math.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * carlson() - Implements Carlson symmetric S_N quadrature
 * legendre_equal_weight() - Legendre weights for polar angles
 * singe_z_ray() - single Legendre weight in z
 * InverseMatrix() -  Compute matrix inverse
 * MatrixMult() - Matrix multiplication
 * permutation() - Checks for permutations of previously assigned vector
 * gauleg() -  Gauss-Legendre quadrature for Numerical Recipes
 * ludcmp_nr() - LU decomposition from Numerical Recipes
 * lubksb_nr() - Backward substitution from Numerical Recipies
 *============================================================================*/

void carlson(RadGridS *pRG, const int nmu);
void legendre_equal_weight(RadGridS *pRG, const int nmu, const int nphi);
void single_z_ray(RadGridS *pRG, const int nmu);
void InverseMatrix(Real **a, int n, Real **b);
void MatrixMult(Real **a, Real *b, int m, int n, Real *c);
int permutation(int i, int j, int k, int **pl, int np);
void gauleg(Real x1, Real x2,  Real *x, Real *w, int n);
void ludcmp_nr(Real **a, int n, int *indx, Real *d);
void lubksb_nr(Real **a, int n, int *indx, Real b[]);

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int qmeth, const int nmu)
 *  \brief Call appropriate function to initialize angles/quadratures */
void init_angles(RadGridS *pRG, const int qmeth, const int outflag)
{
  int nmu, nphi;
  char *radbl;

  if (outflag == 0) radbl = "radiation"; else radbl = "radiation_output";

  switch(qmeth) {

  case 1: /* Carlson symmetric S_N method (default) */
    nmu = par_geti_def(radbl,"nmu",0);
    carlson(pRG,nmu);
    break;

  case 2: /* Legendre Equal Weight */
    nmu = par_geti_def(radbl,"nmu",0);
    nphi = par_geti_def(radbl,"nphi",1);
    legendre_equal_weight(pRG,nmu,nphi);
    break;

  case 10: /* single mu_z */
    nmu = par_geti_def(radbl,"nmu",0);
    single_z_ray(pRG,nmu);
    break;

  default:
    ath_error("[init_agnles]: ang_quad = %d.\n",qmeth);
  }
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int qmeth, const int nmu)
 *  \brief Initialize angles/angle quadratures using Carlson's symmetric
 *  S_N method */
void carlson(RadGridS *pRG, const int nmu)
{

  int nDim;
  int np, ip, iang;
  int i,j,k,l,m;
  Real deltamu, Wsum, wsum, W2;
  Real *mu2tmp = NULL, **mutmp = NULL, *mutmp1d = NULL;
  Real *Wtmp = NULL, *wtmp = NULL;
  Real **pmat = NULL, **pinv = NULL, *wpf = NULL;
  int **pl = NULL, *plab = NULL;

/* number of dimensions in RadGrid. */
  nDim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

  //printf("ndim, nmu: %d %d\n",nDim,nmu);
  if(nDim == 1) {
/* 1D  compute angles using gaussian quadarature */
    mutmp1d = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
    if (mutmp1d == NULL) goto on_error1;
    wtmp = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
    if (wtmp == NULL) goto on_error2;

    gauleg(-1.0, 1.0, mutmp1d, wtmp, 2*nmu);

    for(i=nmu; i<2*nmu; i++) {
      pRG->wmu[i-nmu] = 0.5 * wtmp[i];
      pRG->mu[0][i-nmu][0] = mutmp1d[i];
      pRG->mu[1][i-nmu][0] = -mutmp1d[i];
    }
    free_1d_array(wtmp);
    free_1d_array(mutmp1d);
  } else {
/* 2D and 3D:  compute angles and weights for angular quadratures following
 * the algorithm described in Bruls et al. 1999, A&A, 348, 233 */

    mu2tmp = (Real *)calloc_1d_array(nmu,sizeof(Real));
    if (mu2tmp == NULL) goto on_error3;
    mutmp = (Real **)calloc_2d_array(pRG->nang,3,sizeof(Real));
    if (mutmp == NULL) goto on_error4;
    Wtmp = (Real *)calloc_1d_array(nmu-1,sizeof(Real));
    if (Wtmp == NULL) goto on_error5;
    wtmp = (Real *)calloc_1d_array(nmu,sizeof(Real));
    if (wtmp == NULL) goto on_error6;

/* first compute polar weights and angles */
    if (nmu <= 6) {
      deltamu = 2.0 / (2 * nmu - 1);
      mu2tmp[0] = 1.0 / (3.0 * (2 * nmu - 1));
      for (i=1; i<nmu; i++) {
        mu2tmp[i] = mu2tmp[i-1] + deltamu;
      }
    } else {
      mu2tmp[0] = 1.0 / SQR((Real)nmu-1.0);
      deltamu = (1.0 - 3.0 * mu2tmp[0]) / ((Real)nmu-1.0);
      for (i=1; i<nmu; i++) {
        mu2tmp[i] = mu2tmp[i-1] + deltamu;
      }
    }

    W2 = 4.0 * mu2tmp[0];
    Wsum = Wtmp[0] = sqrt(W2);
    for (i=1; i<nmu-2; i++) {
      W2 += deltamu;
      Wsum += Wtmp[i] = sqrt(W2);
    }
    if (nmu > 2) Wtmp[nmu-2] = 2.0*(nmu-1)/3.0 - Wsum;

    wsum = wtmp[0] = Wtmp[0];
    for (i=1; i<nmu-1; i++) {
      wsum += wtmp[i] = Wtmp[i] - Wtmp[i-1];
    }
    wtmp[nmu-1] = 1.0 - wsum;

/* Next, set up system of equations for determining how polar weights
 * are distributed in azimuth (along circles of section), subject to
 * the constraint that members of permutation families have identical
 * weights */

    pmat = (Real **)calloc_2d_array(nmu,nmu,sizeof(Real));
    if (pmat == NULL) goto on_error7;
    pinv = (Real **)calloc_2d_array(nmu-1,nmu-1,sizeof(Real));
    if (pinv == NULL) goto on_error8;
    plab = (int *)calloc_1d_array(pRG->nang,sizeof(int));
    if (plab == NULL) goto on_error9;
    pl = (int **)calloc_2d_array(nmu,3,sizeof(int));
    if (pl == NULL) goto on_error10;
    wpf = (Real *)calloc_1d_array(nmu-1,sizeof(Real));
    if (wpf == NULL) goto on_error11;

    np = 0;
    iang = 0;
    for (i=0; i<nmu; i++) {
      for (j=0; j<nmu; j++) {
        for (k=0; k<nmu; k++) {
          if (i + j + k == nmu - 1) {
/* assign cosines to temporary array grid */
            mutmp[iang][0] = sqrt(mu2tmp[j]);
            mutmp[iang][1] = sqrt(mu2tmp[k]);
            mutmp[iang][2] = sqrt(mu2tmp[i]);
            if (nmu <= 6) {
              ip=permutation(i,j,k,pl,np);
              if (ip == -1) {
                pl[np][0] = i;
                pl[np][1] = j;
                pl[np][2] = k;
                pmat[i][np] += 1.0;
                plab[iang] = np;
                np++;
              } else {
                pmat[i][ip] += 1.0;
                plab[iang] = ip;
              }
            }
            iang++;
          }
        }
      }
    }

    if (nmu <= 6) {
/* Use Bruls/Carlsson formulation */
      if (nmu > 1) {
/*  Invert matrix of permutations families */
        InverseMatrix(pmat,nmu-1,pinv);
/* Solve for and assign weights for each permutation family */
        MatrixMult(pinv,wtmp,nmu-1,nmu-1,wpf);
        for (i=0; i<pRG->nang; i++)
          pRG->wmu[i] = wpf[plab[i]];
      } else
        pRG->wmu[0] = 1.0;
    } else {
/* Use equal weights for all angles */
      for (i=0; i<pRG->nang; i++)
        pRG->wmu[i] = 1.0/(Real)pRG->nang;
    }

/*  assign angles to RadGrid elements */
    if (nDim == 2) {
      for (i=0; i<pRG->nang; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<2; k++) {
            l=2*j+k;
            if (k == 0)
              pRG->mu[l][i][0] =  mutmp[i][0];
            else
              pRG->mu[l][i][0] = -mutmp[i][0];
            if (j == 0)
              pRG->mu[l][i][1] =  mutmp[i][1];
            else
              pRG->mu[l][i][1] = -mutmp[i][1];
            pRG->mu[l][i][2] =  mutmp[i][2];
          }
        }
        pRG->wmu[i] *= 0.25;
      }
    } else if (nDim == 3) {
      for (i=0; i<pRG->nang; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<2; k++) {
            for (l=0; l<2; l++) {
              m=4*j+2*k+l;
              if (l == 0)
                pRG->mu[m][i][0] =  mutmp[i][0];
              else
                pRG->mu[m][i][0] = -mutmp[i][0];
              if (k == 0)
                pRG->mu[m][i][1] =  mutmp[i][1];
              else
                pRG->mu[m][i][1] = -mutmp[i][1];
              if (j == 0)
                pRG->mu[m][i][2] =  mutmp[i][2];
              else
                pRG->mu[m][i][2] = -mutmp[i][2];
            }
          }
        }
        pRG->wmu[i] *= 0.125;
      }
    }
/* deallocate temporary arrays */
    free_1d_array(wpf);
    free_2d_array(pl);
    free_1d_array(plab);
    free_2d_array(pinv);
    free_2d_array(pmat);
    free_1d_array(wtmp);
    free_1d_array(Wtmp);
    free_2d_array(mutmp);
    free_1d_array(mu2tmp);
  }
/* end of multidimensional quadratures */


  return;

/*--- Error messages ---------------------------------------------------------*/

 on_error11:
  if (nDim > 1) free_1d_array(wpf);
 on_error10:
  if (nDim > 1) free_2d_array(pl);
 on_error9:
  if (nDim > 1) free_1d_array(plab);
 on_error8:
  if (nDim > 1) free_2d_array(pinv);
 on_error7:
  if (nDim > 1) free_2d_array(pmat);
 on_error6:
  if (nDim > 1) free_1d_array(wtmp);
 on_error5:
  if (nDim > 1) free_1d_array(Wtmp);
 on_error4:
  if (nDim > 1) free_2d_array(mutmp);
 on_error3:
  if (nDim > 1) free_1d_array(mu2tmp);
 on_error2:
  if (nDim ==  1) free_1d_array(wtmp);
 on_error1:
  if (nDim ==  1) free_1d_array(mutmp1d);
  ath_error("[carlson]: Error allocating memory\n");

}

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int qmeth, const int nmu)
 *  \brief Initialize angles/angle quadratures using Gauss-Legendre
 *  quadrature */

void legendre_equal_weight(RadGridS *pRG, const int nmu, const int nphi)
{

  int nDim;
  int i,j,k,l,m;
  Real **mutmp = NULL, *mutmp1d = NULL;
  Real *wtmp1d = NULL;
  Real dphi, sintheta;

/* number of dimensions in RadGrid. */
  nDim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

/* 1D  compute angles using gaussian quadarature */
  mutmp1d = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
  if (mutmp1d == NULL) goto on_error1;
  wtmp1d = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
  if (wtmp1d == NULL) goto on_error2;

  gauleg(-1.0, 1.0, mutmp1d, wtmp1d, 2*nmu);

  if (nDim == 1) {
    for(i=0; i<nmu; i++) {
      pRG->wmu[i] = 0.5 * wtmp1d[i+nmu];
      pRG->mu[0][i][0] = mutmp1d[i+nmu];
      pRG->mu[1][i][0] = -mutmp1d[i+nmu];
    }
  } else {
    mutmp = (Real **)calloc_2d_array(pRG->nang,3,sizeof(Real));
    if (mutmp == NULL) goto on_error3;

    k=0;
    for(i=0; i<nmu; i++) {
      dphi = 0.5 * PI / (Real) (2*(nphi+i));
      l=2*nmu-1-i;
      sintheta = sqrt(1.0 - mutmp1d[l] * mutmp1d[l]);
      for(j=0; j<=(i+nphi-1); j++) {
        mutmp[k][0] = sintheta * cos( dphi * (Real)(2*j+1));
        mutmp[k][1] = sintheta * sin( dphi * (Real)(2*j+1));
        mutmp[k][2] = mutmp1d[l];
        pRG->wmu[k] = wtmp1d[l] / (Real)(i+nphi);
        //      printf("wmu %d %d %g %g\n",k,i,pRG->wmu[k],wtmp1d[l]);
        //      printf("mu: %d %d %g %g %g %g\n",k,i,mutmp[k][0],mutmp[k][1],mutmp[k][2],dphi * (Real)(2*j+1));
        k++;
      }
    }

/*  assign angles to RadGrid elements */
    if (nDim == 2) {
      for (i=0; i<pRG->nang; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<2; k++) {
            l=2*j+k;
            if (k == 0)
              pRG->mu[l][i][0] =  mutmp[i][0];
            else
              pRG->mu[l][i][0] = -mutmp[i][0];
            if (j == 0)
              pRG->mu[l][i][1] =  mutmp[i][1];
            else
              pRG->mu[l][i][1] = -mutmp[i][1];
            pRG->mu[l][i][2] =  mutmp[i][2];
          }
        }
        pRG->wmu[i] *= 0.25;
      }
    } else if (nDim == 3) {
      for (i=0; i<pRG->nang; i++) {
        for (j=0; j<2; j++) {
          for (k=0; k<2; k++) {
            for (l=0; l<2; l++) {
              m=4*j+2*k+l;
              if (l == 0)
                pRG->mu[m][i][0] =  mutmp[i][0];
              else
                pRG->mu[m][i][0] = -mutmp[i][0];
              if (k == 0)
                pRG->mu[m][i][1] =  mutmp[i][1];
              else
                pRG->mu[m][i][1] = -mutmp[i][1];
              if (j == 0)
                pRG->mu[m][i][2] =  mutmp[i][2];
              else
                pRG->mu[m][i][2] = -mutmp[i][2];
            }
          }
        }
        pRG->wmu[i] *= 0.125;
      }
    }
    free_2d_array(mutmp);
  }

/* Free temporary arrays */
  free_1d_array(wtmp1d);
  free_1d_array(mutmp1d);

  return;

/*--- Error messages ---------------------------------------------------------*/

 on_error3:
  free_2d_array(mutmp);
 on_error2:
  free_1d_array(wtmp1d);
 on_error1:
  free_1d_array(mutmp1d);
  ath_error("[legendre_equal_weight]: Error allocating memory\n");

}

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int qmeth, const int nmu)
 *  \brief Initialize angles/angle quadratures using Carlson's symmetric
 *  S_N method */
void single_z_ray(RadGridS *pRG, const int nmu)
{

  int nDim;
  int i,j,k,l;
  Real dphi, mu, sintheta;
  Real **mutmp;

  /* number of dimensions in RadGrid. */
  nDim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

  if (nDim == 1) {
    ath_error("[single_z_ray]: ang_quad = 10 should be used only for 2D "
              "problems.\n");
  } else if (nDim == 2) {
    mutmp = (Real **)calloc_2d_array(pRG->nang,3,sizeof(Real));
    if (mutmp == NULL)
      ath_error("[single_z_ray]: Error allocating memory\n");

    mu=1.0/sqrt((Real)3.0);
    sintheta = sqrt((Real)2.0)/sqrt((Real)3.0);
    dphi = 0.5 * PI / (Real) (2*nmu);
    for(i=0; i<nmu; i++) {
      mutmp[i][0] = sintheta * cos( dphi * (Real)(2*i+1));
      mutmp[i][1] = sintheta * sin( dphi * (Real)(2*i+1));
      mutmp[i][2] = mu;
      pRG->wmu[i] = 1 / (Real)nmu;
    }
    for (i=0; i<pRG->nang; i++) {
      for (j=0; j<2; j++) {
        for (k=0; k<2; k++) {
          l=2*j+k;
          if (k == 0)
            pRG->mu[l][i][0] =  mutmp[i][0];
          else
            pRG->mu[l][i][0] = -mutmp[i][0];
          if (j == 0)
            pRG->mu[l][i][1] =  mutmp[i][1];
          else
            pRG->mu[l][i][1] = -mutmp[i][1];
          pRG->mu[l][i][2] =  mutmp[i][2];
        }
      }
      pRG->wmu[i] *= 0.25;
    }
    free_2d_array(mutmp);
  } else if (nDim == 3) {
    ath_error("[single_z_ray]: ang_quad = 10 should be used only for 2D "
              "problems.\n");
  }

}


/*----------------------------------------------------------------------------*/
/*! \fn int permutation(int i, int j, int k, int **pl, int np)
 * Checks if an indicies triplet is a permutation of vectors in pl[][]
 * (previously loaded triplets) and returns index of permutation if matched
 * or -1 if no match */
int permutation(int i, int j, int k, int **pl, int np)
{
  int ip=-1;
  int l,m,n,o;

/*  This routine is only called at initialization so brute force
 *  algorithm is fine */


  for(l=0; l<np; l++) {
/* check each permutation in the table */
    for(m=0; m<3; m++)
      if(i == pl[l][m])
        for(n=0; n<3; n++)
          if(n != m)
            if(j == pl[l][n])
              for(o=0;o<3;o++)
                if((o != m) && (o != n))
                  if(k == pl[l][o])
                    ip = l;
  }

  return ip;
}

/*----------------------------------------------------------------------------*/
/*! \fn void ludcmp_nr(Real **a, int n, int *indx, Real *d)
 * LU decomposition from Numerical Recipes
 * Using Crout's method with partial pivoting
 * a is the input matrix, and is returned with LU decomposition readily made,
 * n is the matrix size, indx records the history of row permutation,
 * whereas d =1(-1) for even(odd) number of permutations.
 */
void ludcmp_nr(Real **a, int n, int *indx, Real *d)
{
  int i,imax,j,k;
  Real big,dum,sum,temp;
  Real *rowscale;  /* the implicit scaling of each row */

  rowscale = (Real*)calloc_1d_array(n, sizeof(Real));
  *d=1.0;  /* No row interchanges yet */

  for (i=0;i<n;i++)
    { /* Loop over rows to get the implicit scaling information */
      big=0.0;
      for (j=0;j<n;j++)
        if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0) ath_error("[LUdecomp]:Input matrix is singular!");
      rowscale[i]=1.0/big;  /* Save the scaling */
    }

  for (j=0;j<n;j++) { /* Loop over columns of Crout's method */
    /* Calculate the upper block */
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    /* Calculate the lower block (first step) */
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      /* search for the largest pivot element */
      if ( (dum=rowscale[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    /* row interchange */
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      rowscale[imax]=rowscale[j];
    }
    indx[j]=imax; /* record row interchange history */
    /* Calculate the lower block (second step) */
    if (a[j][j] == 0.0) a[j][j]=TINY_NUMBER;
    dum=1.0/(a[j][j]);
    for (i=j+1;i<n;i++) a[i][j] *= dum;
  }
  free(rowscale);
}

/*----------------------------------------------------------------------------*/
/*! \fn void lubksb_nr(Real **a, int n, int *indx, Real b[])
 *  Backward substitution (from numerical recipies)
 *  a is the input matrix done with LU decomposition, n is the matrix size
 *  indx id the history of row permutation
 *  b is the vector on the right (AX=b), and is returned with the solution
 */
void lubksb_nr(Real **a, int n, int *indx, Real b[])
{
  int i,ii=-1,ip,j;
  Real sum;
  /* Solve L*y=b */
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  /* Solve U*x=y */
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void InverseMatrix(Real **a, int n, Real **b)
 *  Inverse matrix solver
 *  a: input matrix; n: matrix size, b: return matrix
 *  Note: the input matrix will be DESTROYED
 */
void InverseMatrix(Real **a, int n, Real **b)
{
  int i,j,*indx;
  Real *col,d;

  indx = (int*)calloc_1d_array(n, sizeof(int));
  col = (Real*)calloc_1d_array(n, sizeof(Real));

  ludcmp_nr(a,n,indx,&d);

  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_nr(a, n, indx, col);
    for (i=0; i<n; i++)    b[i][j] = col[i];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MatrixMult(Real **a, Real *b, int m, int n, Real *c)
 *  Matrix multiplication: a(m*n) * b(n) = c(m) */
void MatrixMult(Real **a, Real *b, int m, int n, Real *c)
{
  int i, j;
  for (i=0; i<m; i++) {
    c[i] = 0.0;
    for (j=0; j<n; j++) c[i] += a[i][j] * b[j];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
 * gauss-legendre weight routine from numerical recipes */
void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
{

  Real eps = 3.0e-14;
  Real xm, xl, z, z1;
  Real p1, p2, p3, pp;
  int m, i, j;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=1; i<=m; i++) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * (Real)j - 1.0) * z * p2 - ((Real)j - 1.0) * p3) / (Real)j;
      }
      pp = (Real)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}

#endif /* RADIATION_TRANSFER */
