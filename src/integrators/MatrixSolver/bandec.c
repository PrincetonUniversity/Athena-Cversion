#include "../../copyright.h"
/*==============================================================================
 * FILE: bandec.c
 *
 * PURPOSE: This file is adopted from Numerical Recipes to solve a set of
 * linear equations with LU decomposition. This is designed especially
 * for band diagnol matrix equations
 * Input matrix is stored as compact form.
 * bandec do the decomposition and banbks do the backsubstitution
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   bandec()
 *   banbks()
 *============================================================================*/

#define SWAP(a,b) {dum=(a);(a)=(b);(b)=dum;}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../prototypes.h"
#include "../../prototypes.h"

#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)

void bandec(Real **a, unsigned long n, int m1, int m2, Real **al,
            unsigned long indx[], Real *d)
{
  unsigned long i,j,k,l;
  int mm;
  Real dum;

  mm=m1+m2+1;
  l=m1;
  for (i=1;i<=m1;i++) {
    for (j=m1+2-i;j<=mm;j++) a[i][j-l]=a[i][j];
    l--;
    for (j=mm-l;j<=mm;j++) a[i][j]=0.0;
  }
  *d=1.0;
  l=m1;
  for (k=1;k<=n;k++) {
    dum=a[k][1];
    i=k;
    if (l < n) l++;
    for (j=k+1;j<=l;j++) {
      if (fabs(a[j][1]) > fabs(dum)) {
        dum=a[j][1];
        i=j;
      }
    }
    indx[k]=i;
    if (dum == 0.0) {
      a[k][1]=TINY_NUMBER;
      fprintf(stderr,"[bandec.c]: Matrix is singular!\n");
    }
    if (i != k) {
      *d = -(*d);
      for (j=1;j<=mm;j++) SWAP(a[k][j],a[i][j])
                            }
    for (i=k+1;i<=l;i++) {
      dum=a[i][1]/a[k][1];
      al[k][i-k]=dum;
      for (j=2;j<=mm;j++) a[i][j-1]=a[i][j]-dum*a[k][j];
      a[i][mm]=0.0;
    }
  }
}


void banbks(Real **a, unsigned long n, int m1, int m2, Real **al,
            unsigned long indx[], Real b[])
{
  unsigned long i,k,l;
  int mm;
  Real dum;

  mm=m1+m2+1;
  l=m1;
  for (k=1;k<=n;k++) {
    i=indx[k];
    if (i != k) SWAP(b[k],b[i])
                  if (l < n) l++;
    for (i=k+1;i<=l;i++) b[i] -= al[k][i-k]*b[k];
  }
  l=1;
  for (i=n;i>=1;i--) {
    dum=b[i];
    for (k=2;k<=l;k++) dum -= a[i][k]*b[k+i-1];
    b[i]=dum/a[i][1];
    if (l < mm) l++;
  }
}
#undef SWAP

#endif

