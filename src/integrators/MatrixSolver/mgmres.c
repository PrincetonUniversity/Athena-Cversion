
/* Matrix solver for general sparse matrix with restarted GMRES method.
 * Adopted and modified from website:
 * http://people.sc.fsu.edu/~jburkardt/c_src/mgmres/mgmres.html
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../../defs.h"
#include "../../athena.h"
#include "../../globals.h"
#include "../prototypes.h"
#include "../../prototypes.h"

#ifdef rad_hydro

void mgmres_st ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel )

/******************************************************************************/
/*
  Purpose:

    MGMRES_ST applies restarted GMRES to a matrix in sparse triplet form.

  Discussion:

    The matrix A is assumed to be stored in sparse triplet form.  Only 
    the nonzero entries of A are stored.  For instance, the 
    K-th nonzero entry in the matrix is stored by:

      A(K) = value of entry,
      IA(K) = row of entry,
      JA(K) = column of entry.

    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
    corrections to the code on 31 May 2007.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    18 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994.
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 2003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the linear system.

    Input, int NZ_NUM, the number of nonzero matrix values.

    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices of 
    the matrix values.

    Input, double A[], the matrix values.

    Input/output, double X[N]; on input, an approximation to
    the solution.  On output, an improved approximation.

    Input, double RHS[N], the right hand side of the linear system.

    Input, int ITR_MAX, the maximum number of (outer) iterations to take.

    Input, int MR, the maximum number of (inner) iterations to take.
    MR must be less than N.

    Input, double TOL_ABS, an absolute tolerance applied to the
    current residual.

    Input, double TOL_REL, a relative tolerance comparing the
    current residual to the initial residual.
*/
{
  double av;
  double *c;
  double delta = 1.0e-03;
  double *g;
  double **h;
  double htmp;
  int i;
  int itr;
  int itr_used;
  int j;
  int k;
  int k_copy;
  double mu;
  double *r;
  double rho;
  double rho_tol;
  double *s;
  double **v;
  int verbose = 1;
  double *y;

  if ( n < mr )
      ath_error("n<mr in matrix solver mgmres!\n");

  itr_used = 0;

  if((c = ( double * ) malloc ( mr * sizeof ( double ) ))==NULL)
	{goto on_error;}

  if((g = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) ))==NULL)
	{goto on_error;}
  
  if((h = (double **) malloc ((mr + 1) * sizeof (double * ) ))==NULL)
	{goto on_error;}

  for(i=0; i< mr+1; i++)
	if((h[i] = (double *) malloc ( mr * sizeof (double ) ))==NULL)
		{goto on_error;}
  
  if((r = ( double * ) malloc ( n * sizeof ( double ) ))==NULL)
	{goto on_error;}

  if((s = ( double * ) malloc ( mr * sizeof ( double ) ))==NULL)
	{goto on_error;}

  if((v = (double **) malloc ((mr + 1) * sizeof (double * ) ))==NULL)
	{goto on_error;}

  for(i=0; i< mr + 1; i++)
	if((v[i] = (double *) malloc ( n * sizeof (double ) ))==NULL)
		{goto on_error;}


  if((y = ( double * ) malloc ( ( mr + 1 ) * sizeof ( double ) ))==NULL)
	{goto on_error;}

  for ( itr = 0; itr < itr_max; itr++ ) 
  {
    ax_st ( n, nz_num, ia, ja, a, x, r );

    for ( i = 0; i < n; i++ )
    {
      r[i] = rhs[i] - r[i];
    }

    rho = sqrt ( r8vec_dot ( n, r, r ) );

/* * do not print residual in Athena 
    if ( verbose )
    {
      printf ( "  ITR = %8d  Residual = %e\n", itr, rho );
    }
*/
    if ( itr == 0 )
    {
      rho_tol = rho * tol_rel;
    }

    for ( i = 0; i < n; i++ )
    {
      v[0][i] = r[i] / rho;
    }

    g[0] = rho;
    for ( i = 1; i < mr + 1; i++ )
    {
      g[i] = 0.0;
    }

    for ( i = 0; i < mr + 1; i++ )
    {
      for ( j = 0; j < mr; j++ ) 
      {
        h[i][j] = 0.0;
      }
    }

    for ( k = 0; k < mr; k++ )
    {
      k_copy = k;

      ax_st ( n, nz_num, ia, ja, a, v[k], v[k+1] );

      av = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      for ( j = 0; j < k+1; j++ )
      {
        h[j][k] = r8vec_dot ( n, v[k+1], v[j] );
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] - h[j][k] * v[j][i];
        }
      }

      h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );

      if ( ( av + delta * h[k+1][k] ) == av )
      {
        for ( j = 0; j < k+1; j++ )
        {
          htmp = r8vec_dot ( n, v[k+1], v[j] );
          h[j][k] = h[j][k] + htmp;
          for ( i = 0; i < n; i++ ) 
          {
            v[k+1][i] = v[k+1][i] - htmp * v[j][i];
          }
        }
        h[k+1][k] = sqrt ( r8vec_dot ( n, v[k+1], v[k+1] ) );
      }

      if ( h[k+1][k] != 0.0 )
      {
        for ( i = 0; i < n; i++ ) 
        {
          v[k+1][i] = v[k+1][i] / h[k+1][k];
        }
      }

      if ( 0 < k )
      {
        for ( i = 0; i < k + 2; i++ )
        {
          y[i] = h[i][k];
        }
        for ( j = 0; j < k; j++ ) 
        {
          mult_givens ( c[j], s[j], j, y );
        }
        for ( i = 0; i < k + 2; i++ ) 
        {
          h[i][k] = y[i];
        }
      }

      mu = sqrt ( h[k][k] * h[k][k] + h[k+1][k] * h[k+1][k] );
      c[k] = h[k][k] / mu;
      s[k] = -h[k+1][k] / mu;
      h[k][k] = c[k] * h[k][k] - s[k] * h[k+1][k];
      h[k+1][k] = 0.0;
      mult_givens ( c[k], s[k], k, g );

      rho = fabs ( g[k+1] );

      itr_used = itr_used + 1;
/*
      if ( verbose )
      {
        printf ( "  K =   %8d  Residual = %e\n", k, rho );
      }
*/
      if ( rho <= rho_tol && rho <= tol_abs )
      {
        break;
      }
    }

    k = k_copy;

    y[k] = g[k] / h[k][k];
    for ( i = k - 1; 0 <= i; i-- )
    {
      y[i] = g[i];
      for ( j = i+1; j < k + 1; j++ ) 
      {
        y[i] = y[i] - h[i][j] * y[j];
      }
      y[i] = y[i] / h[i][i];
    }

    for ( i = 0; i < n; i++ )
    {
      for ( j = 0; j < k + 1; j++ )
      {
        x[i] = x[i] + v[j][i] * y[j];
      }
    }

    if ( rho <= rho_tol && rho <= tol_abs ) 
    {
      break;
    }
  }
/*
  if ( verbose )
  {
    printf ( "\n" );
    printf ( "MGMRES_ST:\n" );
    printf ( "  Iterations = %d\n", itr_used );
    printf ( "  Final residual = %e\n", rho );
  }
*/

  free ( c );
  free ( g );

  for(i=0; i<mr+1; i++)
	free(h[i]);

  free(h);

  free ( r );
  free ( s );

  for(i=0; i< mr+1; i++)
  free(v[i]);

  free(v);

  free ( y );

  return;

  on_error:
    	
    	ath_error("[mgmres]: malloc returned a NULL pointer\n");
	
} 

/******************************************************************************/

void ax_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] )

/******************************************************************************/
/*
  Purpose:

    AX_ST computes A*x for a matrix stored in sparse triplet form.

  Discussion:

    The matrix A is assumed to be stored in sparse triplet format.  Only
    the nonzero entries of A are stored.  For instance, the K-th nonzero
    entry in the matrix is stored by:

      A(K) = value of entry,
      IA(K) = row of entry,
      JA(K) = column of entry.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    17 July 2007

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, int N, the order of the system.

    Input, int NZ_NUM, the number of nonzeros.

    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
    of the matrix values.

    Input, double A[NZ_NUM], the matrix values.

    Input, double X[N], the vector to be multiplied by A.

    Output, double W[N], the value of A*X.
*/
{
  int i;
  int j;
  int k;

  for ( i = 0; i < n; i++ )
  {
    w[i] = 0.0;
  }

  for ( k = 0; k < nz_num; k++ )
  {
    i = ia[k];
    j = ja[k];
    w[i] = w[i] + a[k] * x[j];
  }

  return;
}
/******************************************************************************/




double r8vec_dot ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:

    R8VEC_DOT computes the dot product of a pair of R8VEC's.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    26 July 2007

  Author:

    John Burkardt

  Parameters:

    Input, int N, the number of entries in the vectors.

    Input, double A1[N], A2[N], the two vectors to be considered.

    Output, double R8VEC_DOT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/



void mult_givens ( double c, double s, int k, double *g )

/******************************************************************************/
/*
  Purpose:

    MULT_GIVENS applies a Givens rotation to two vector elements.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    10 August 2006

  Author:

    Lili Ju

  Reference:

    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
    Charles Romine, Henk van der Vorst,
    Templates for the Solution of Linear Systems:
    Building Blocks for Iterative Methods,
    SIAM, 1994,
    ISBN: 0898714710,
    LC: QA297.8.T45.

    Tim Kelley,
    Iterative Methods for Linear and Nonlinear Equations,
    SIAM, 2004,
    ISBN: 0898713528,
    LC: QA297.8.K45.

    Yousef Saad,
    Iterative Methods for Sparse Linear Systems,
    Second Edition,
    SIAM, 20003,
    ISBN: 0898715342,
    LC: QA188.S17.

  Parameters:

    Input, double C, S, the cosine and sine of a Givens
    rotation.

    Input, int K, indicates the location of the first vector entry.

    Input/output, double G[K+2], the vector to be modified.  On output,
    the Givens rotation has been applied to entries G(K) and G(K+1).
*/
{
  double g1;
  double g2;

  g1 = c * g[k] - s * g[k+1];
  g2 = s * g[k] + c * g[k+1];

  g[k]   = g1;
  g[k+1] = g2;

  return;
}
/******************************************************************************/

#endif

