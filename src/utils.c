#include "copyright.h"
/*==============================================================================
 * FILE: utils.c
 *
 * PURPOSE: A variety of useful utility functions.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ath_strdup()     - not supplied by fancy ANSI C, but ok in C89 
 *   ath_gcd()        - computes greatest common divisor by Euler's method
 *   ath_big_endian() - run-time detection of endianism of the host cpu
 *   ath_bswap()      - fast byte swapping routine
 *   ath_error()      - fatal error routine
 *   minmax1()        - fast Min/Max for a 1d array using registers
 *   minmax2()        - fast Min/Max for a 2d array using registers
 *   minmax3()        - fast Min/Max for a 3d array using registers
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* ath_strdup: this is really strdup(), but strdup is not available in 
 *   ANSI  (-pendantic or -ansi will leave it undefined in gcc)
 *   much like allocate.
 */

char *ath_strdup(const char *in)
{
  char *out = (char *)malloc((1+strlen(in))*sizeof(char));
  if(out == NULL) {
    fprintf(stderr,"ath_strdup: failed to alloc %d\n",(int)(1+strlen(in)));
    return NULL; /* malloc failed */
  }
  return strcpy(out,in);
}

/*----------------------------------------------------------------------------*/
/* ath_gcd: Calculate the Greatest Common Divisor by Euler's method
 */

int ath_gcd(int a, int b)
{
  int c;
  if(b>a) {c=a; a=b; b=c;} 
  while((c=a%b)) {a=b; b=c;}
  return b;
}

/*----------------------------------------------------------------------------*/
/* ath_big_endian:  return 1 if the machine is big endian (e.g. Sun, PowerPC)
 * return 0 if not (e.g. Intel)
 */

int ath_big_endian(void)
{
  short int n = 1;
  char *ep = (char *)&n;

  return (*ep == 0); /* Returns 1 on a big endian machine */
}

/*----------------------------------------------------------------------------*/
/* ath_bswap: swap bytes, code stolen from NEMO  
 */
 
void ath_bswap(void *vdat, int len, int cnt)
{
  char tmp, *dat = (char *) vdat;
  int k;
 
  if (len==1)
    return;
  else if (len==2)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[1];  dat[1] = tmp;
      dat += 2;
    }
  else if (len==4)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[3];  dat[3] = tmp;
      tmp = dat[1];  dat[1] = dat[2];  dat[2] = tmp;
      dat += 4;
    }
  else if (len==8)
    while (cnt--) {
      tmp = dat[0];  dat[0] = dat[7];  dat[7] = tmp;
      tmp = dat[1];  dat[1] = dat[6];  dat[6] = tmp;
      tmp = dat[2];  dat[2] = dat[5];  dat[5] = tmp;
      tmp = dat[3];  dat[3] = dat[4];  dat[4] = tmp;
      dat += 8;
    }
  else {  /* the general SLOOOOOOOOOW case */
    for(k=0; k<len/2; k++) {
      tmp = dat[k];
      dat[k] = dat[len-1-k];
      dat[len-1-k] = tmp;
    }
  }
}

/*----------------------------------------------------------------------------*/
/* ath_error: Terminate execution and output error message
 *  Uses variable-length argument lists provided in <stdarg.h>
 */

void ath_error(char *fmt, ...)
{
  va_list ap;
 
  fprintf(stderr,"### Fatal error: ");   /* prefix */
  va_start(ap, fmt);              /* ap starts with string 'fmt' */
  vfprintf(stderr, fmt, ap);      /* print out on stderr */
  fflush(stderr);                 /* flush it NOW */
  va_end(ap);                     /* end varargs */

#ifdef MPI_PARALLEL
  MPI_Abort(MPI_COMM_WORLD, 1);
#endif

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*/
/* minmax1,2,3: return the min and max of a 1D, 2D or 3D array using registers
 *  Works on data of type float, not Real.
 */

void minmax1(float *data, int nx1, float *dmino, float *dmaxo)
{
  int i;
  register float dmin, dmax;

  dmin = dmax = data[0];
  for (i=0; i<nx1; i++) {
    dmin = MIN(dmin,data[i]);
    dmax = MAX(dmax,data[i]);
  }
  *dmino = dmin;
  *dmaxo = dmax;
}

void minmax2(float **data, int nx2, int nx1, float *dmino, float *dmaxo)
{
  int i,j;
  register float dmin, dmax;

  dmin = dmax = data[0][0];
  for (j=0; j<nx2; j++) {
    for (i=0; i<nx1; i++) {
      dmin = MIN(dmin,data[j][i]);
      dmax = MAX(dmax,data[j][i]);
    }
  }
  *dmino = dmin;
  *dmaxo = dmax;
}

void minmax3(float ***data,int nx3,int nx2,int nx1, float *dmino, float *dmaxo)
{
  int i,j,k;
  register float dmin, dmax;

  dmin = dmax = data[0][0][0];
  for (k=0; k<nx3; k++) {
    for (j=0; j<nx2; j++) {
      for (i=0; i<nx1; i++) {
	dmin = MIN(dmin,data[k][j][i]);
	dmax = MAX(dmax,data[k][j][i]);
      }
    }
  }
  *dmino = dmin;
  *dmaxo = dmax;
}
