#include "../copyright.h"
/*==============================================================================
 * FILE: selfg_fft_obc.c
 *
 * PURPOSE: Contains functions to solve Poisson's equation for self-gravity in
 *   3D using FFTs, using OPEN BCs in all three directions 
 *
 *
 *   The function uses FFTW3.x, and for MPI parallel use Steve Plimpton's
 *   block decomposition routines added by N. Lemaster to /athena/fftsrc.
 *   This means to use these fns the code must be
 *     (1) configured with --with-gravity=fft_obc --enable-fft
 *     (2) compiled with links to FFTW libraries (may need to edit Makeoptions)
 *
 *   For PERIODIC BCs, use selfg_fft() functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_fft_obc_3d() - 3D Poisson solver using FFTs
 *   selfg_fft_obc_3d_init() - initializes FFT plans for 3D
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SELF_GRAVITY_USING_FFT_OBC

#ifndef FFT_ENABLED
#error self gravity with FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL;


#ifdef STATIC_MESH_REFINEMENT
#error self gravity with FFT not yet implemented to work with SMR
#endif


/*----------------------------------------------------------------------------*/
/* selfg_fft_3d:
 *   Only works for uniform grid, periodic boundary conditions
 */

void selfg_fft_obc_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ip, jp, kp;
  int hNx1=pD->Nx[0]/2,hNx2=pD->Nx[1]/2,hNx3=pD->Nx[2]; 
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3);
  Real dkx,dky,dkz;
  Real den; 
  Real den_parent=0.;

/* To compute kx,ky,kz, note that indices relative to whole Domain are needed */
  dkx = 2.0*PI/(double)(pD->Nx[0]);
  dky = 2.0*PI/(double)(pD->Nx[1]);
  dkz = 2.0*PI/(double)(pD->Nx[2]);


// NOTE: the following is done only once and could be moved to 
//   selfg_fft_3d_init:
/* coefficients for Poisson kernel */
  static int coeff_set=0;
  static  Real ***Geee=NULL,***Gooo=NULL,***Goee=NULL,***Geoo=NULL,***Geoe=NULL, ***Goeo=NULL,***Geeo=NULL,***Gooe=NULL;
  if (!coeff_set){
/* first time through: compute coefficients of poisson kernel */

/*   allocates memory for Poisson kernel coefficient arrays */
  if ((Geee = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Gooo = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Goee = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Geoo = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Geoe = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Goeo = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Geeo = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");
  if ((Gooe = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[poiss_coeff]: malloc returned a NULL pointer\n");

/* Compute potential coeffs in k space. Zero wavenumber is special
   case; need to avoid divide by zero.  Since we only do this once in setup, 
   if statement in loop is cleaner than multiple loops.    */
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
   if ((k-ks+pG->Disp[2])==0 && (j-js+pG->Disp[1])==0 && (i-is+pG->Disp[0])==0) 
	   Geee[0][0][0] =0.0;
          else{
	   Geee[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos((     (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos((     (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos((     (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	    }
	   Gooo[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos(( 0.5+(i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos(( 0.5+(j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos(( 0.5+(k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Goee[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos(( 0.5+(i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos((     (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos((     (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Geoo[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos((     (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos(( 0.5+(j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos(( 0.5+(k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Geoe[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos((     (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos(( 0.5+(j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos((     (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Goeo[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos(( 0.5+(i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos((     (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos(( 0.5+(k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Geeo[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos((     (i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos((     (j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos(( 0.5+(k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
	   Gooe[i-is][j-js][k-ks] = 1.0/ 
        (((2.0*cos(( 0.5+(i-is)+pG->Disp[0] )*dkx)-2.0)/dx1sq) +
         ((2.0*cos(( 0.5+(j-js)+pG->Disp[1] )*dky)-2.0)/dx2sq) +
         ((2.0*cos((     (k-ks)+pG->Disp[2] )*dkz)-2.0)/dx3sq));
      }
    }
  }
  coeff_set=1; /* done computing coeffs */
  } /* end of one-time-only setup */

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
   for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
    }
  }
}


/* There are eight different terms.  For each term, need to do the 
   following steps:

  (1) fourier transform of 

      (density[i j k] - parent_density[i j k] )* 
  [1 or  (cos(pi i /Nx)+i sin(pi i /Nx)]  for [even or odd] in i
  [1 or  (cos(pi j /Ny)+i sin(pi j /Ny)]  for [even or odd] in j
  [1 or  (cos(pi k /Nx)+i sin(pi k /Nx)]  for [even or odd] in k

  (2) multiply by appropriate Geee Gooo Goee etc.

  (3) take inverse transform

  (4) multiply by 

  [1 or  (cos(pi i /Nx) -i sin(pi i /Nx))]  for [even or odd] in i 
  [1 or  (cos(pi j /Ny) -i sin(pi j /Ny))]  for [even or odd] in j 
  [1 or  (cos(pi k /Nx) -i sin(pi k /Nx))]  for [even or odd] in k


After these eight terms are done, result is added up and multiplied by 
  4 pi G/(8 Nx Ny Nz) 

*/

// NOTE:  MAY NEED TO CHANGE THE SIGN OF ALL THE sin TERMS; NEED TO CHECK THIS


/* First term: eee  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
   }
  }
}
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Geee[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Geee[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential in real space. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] = 
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
        *four_pi_G/(8.*bplan3d->gcnt);
    }
  }
}


/* Second term: ooo  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*(((i-is)+pG->Disp[0])*dkx+
		 ((j-js)+pG->Disp[1])*dky+ 
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*(((i-is)+pG->Disp[0])*dkx+
		 ((j-js)+pG->Disp[1])*dky+ 
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Gooo[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Gooo[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((j-js)+pG->Disp[1])*dky+ 
		  ((k-ks)+pG->Disp[2])*dkz))
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((j-js)+pG->Disp[1])*dky+ 
		  ((k-ks)+pG->Disp[2])*dkz))
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}

/* Third term: oee  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*((i-is)+pG->Disp[0])*dkx)*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*((i-is)+pG->Disp[0])*dkx)*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Goee[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Goee[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*((i-is)+pG->Disp[0])*dkx)
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*((i-is)+pG->Disp[0])*dkx)
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}

/* Fourth term: eoo  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*(((j-js)+pG->Disp[1])*dky+ 
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*(((j-js)+pG->Disp[1])*dky+ 
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Geoo[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Geoo[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*(((j-js)+pG->Disp[1])*dky+ 
		  ((k-ks)+pG->Disp[2])*dkz))
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*(((j-js)+pG->Disp[1])*dky+ 
		  ((k-ks)+pG->Disp[2])*dkz))
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}


/* Fifth term: eoe  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*((j-js)+pG->Disp[1])*dky)*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*((j-js)+pG->Disp[1])*dky)*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Geoe[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Geoe[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*((j-js)+pG->Disp[1])*dky)
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*((j-js)+pG->Disp[1])*dky)
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}


/* Sixth term: ooo  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*(((i-is)+pG->Disp[0])*dkx+		  
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*(((i-is)+pG->Disp[0])*dkx+ 
		 ((k-ks)+pG->Disp[2])*dkz))*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Goeo[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Goeo[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((k-ks)+pG->Disp[2])*dkz))
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((k-ks)+pG->Disp[2])*dkz))
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}
		  

/* Seventh term: eeo  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*((k-ks)+pG->Disp[2])*dkz)*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*((k-ks)+pG->Disp[2])*dkz)*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Geeo[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Geeo[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*((k-ks)+pG->Disp[2])*dkz)
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*((k-ks)+pG->Disp[2])*dkz)
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}


/* Eighth term: ooo  */
/* Subtract off background density, and set up real and imaginary arrays
   based on multiplication by [even or odd] exponential factors for i,j,k */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
	cos(0.5*(((i-is)+pG->Disp[0])*dkx+
		 ((j-js)+pG->Disp[1])*dky))*
        (pG->U[k][j][i].d - den_parent);
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
	sin(0.5*(((i-is)+pG->Disp[0])*dkx+
		 ((j-js)+pG->Disp[1])*dky))*
        (pG->U[k][j][i].d - den_parent);
   }
  }
 }
/* take forward transform */  
  ath_3d_fft(fplan3d, work);
/* compute potential contribution in Fourier space, using pre-computed 
   coefficient*/  
  for (i=is; i<=ie; i++){
   for (j=js; j<=je; j++){
    for (k=ks; k<=ke; k++){
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= 
    Gooe[i-is][j-js][k-ks];
      work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= 
    Gooe[i-is][j-js][k-ks];
    }
  }
}
/* Backward FFT */ 
  ath_3d_fft(bplan3d, work);
/* Multiply by [even or odd] exponential factors for i,j,k and 
   apply constant factor and normalization over total number of cells in Domain.
   Add contribution to potential, keeping just the real part. */  
 */
  for (k=ks; k<=ke; k++){
   for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[k][j][i] += (
       work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
	*cos(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((j-js)+pG->Disp[1])*dky))
     + work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1]
	*sin(0.5*(((i-is)+pG->Disp[0])*dkx+
		  ((j-js)+pG->Disp[1])*dky))
			   )*four_pi_G/(8.*bplan3d->gcnt);

    }
  }
}

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_fft_3d_init:
 *   Initializes plans for forward/backward FFTs, and allocates memory needed
 *   by FFTW 
 *
 *   NOTE: for SMR, allocation of memory for Poisson kernel arrays 
 *      Geee, Gooo, Goee, Geoo, Geoe, Goeo, Geeo, Gooe 
 *   could also be done here
 *
 */

void selfg_fft_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* SELF_GRAVITY_USING_FFT */
