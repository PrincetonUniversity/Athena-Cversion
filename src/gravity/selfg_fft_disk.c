#include "../copyright.h"
/*=============================================================================*/
/*! \file selfg_fft_disk.c
 *  \brief Contains functions to solve Poisson's equation for self-gravity 
 *   in disk symmetry, in 1D, 2D and 3D using FFTs 
 *
 *   For 1D, x1 is perpendicular to the plane
 *   For 2D, x1 is in plane and periodic, and x2 is perpendicular to the plane
 *   For 3D, x1 and x2 are in plane and periodic, and x3 is perpendicular 
 *   to the plane 
 *   
 *
 *   The FFT's use the FFTW3.x libraries, and for MPI parallel use 
 *   Steve Plimpton's block decomposition routines added by N. Lemaster 
 *   to /athena/fftsrc.
 *   This means to use these fns the code must be
 *      (1) configured with --enable-fft
 *      (2) compiled with links to FFTW libraries
 *
 *   For NON-PERIODIC BCs, use selfg_multig() functions.
 *   For FULLY-PERIODIC BCs, use selfg_fft functions
 *
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_by_fft_disk_1d() - actually uses recursion; single processor only
 *   selfg_by_fft_disk_2d() - 2D Poisson solver using FFTs
 *   selfg_by_fft_disk_3d() - 3D Poisson solver using FFTs
 *   selfg_by_fft_disk_2d_init() - initializes FFT plans for 2D
 *   selfg_by_fft_disk_3d_init() - initializes FFT plans for 3D
 *
 *  NOTE:     The functions in selfg_fft assume PERIODIC BC in ALL directions.
 *            The functions here implement OPEN BC in ONE direction and 
 *            PERIODIC BC in the other direction(s). */
/*============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SELF_GRAVITY_USING_FFT_DISK

#ifndef FFT_ENABLED
#error self gravity with FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_2d_fft_plan *fplan2d, *bplan2d;
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL, *work2=NULL;

#ifdef STATIC_MESH_REFINEMENT
#error self gravity with FFT in DISK not yet implemented to work with SMR
#endif

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_disk_1d(DomainS *pD)
 *  \brief Actually uses recursion formula.  
 *  ONLY WORKS FOR SINGLE PROCESSOR! 
 */

void selfg_fft_disk_1d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real total_Phi=0.0,dx1sq = (pG->dx1*pG->dx1);
/* Copy current potential into old */

  for (i=is-nghost; i<=ie+nghost; i++){
    pG->Phi_old[ks][js][i] = pG->Phi[ks][js][i];
  }

/* Compute new potential */

  pG->Phi[ks][js][is] = 0.0;
  for (i=is; i<=ie; i++) {
    pG->Phi[ks][js][is] += pG->U[ks][js][i].d; 
  }

  pG->Phi[ks][js][is  ] *= 0.25*four_pi_G*dx1sq*(float)((pG->Nx[0])-1);
  pG->Phi[ks][js][is+1] = pG->Phi[ks][js][is] + 
           four_pi_G*dx1sq*pG->U[ks][js][is].d - 
           2.*pG->Phi[ks][js][is]/(float)((pG->Nx[0])-1);
  for (i=is+2; i<=ie; i++) {
    pG->Phi[ks][js][i] = four_pi_G*dx1sq*pG->U[ks][js][i-1].d 
      + 2.0*pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i-2];
  }
/* apply open BC in x1 direction to obtain values in ghost zones */
      pG->Phi[ks][js][ie+1] = 2.0*pG->Phi[ks][js][ie] - pG->Phi[ks][js][ie-1] +
	dx1sq*four_pi_G*pG->U[ks][js][ie].d;
      pG->Phi[ks][js][is-1] = 2.0*pG->Phi[ks][js][is] - pG->Phi[ks][js][is+1] +
	dx1sq*four_pi_G*pG->U[ks][js][is].d;
}




/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_disk_2d(DomainS *pD)
 *  \brief Periodic boundary conditions in x1; open bc in x2
 */

void selfg_fft_disk_2d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real dkx;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2);
  Real xmin,xmax,Lperp; 
  static int coeff_set=0;
  static Real **Acoeff=NULL,**Bcoeff=NULL; 
  static Real dky=0.;
  int ip;

/* first time through: compute coefficients of poisson kernel and dky*/

  if (!coeff_set){

/*   allocates memory for Acoeff and Bcoeff arrays */
  if ((Acoeff = (Real**)calloc_2d_array(pG->Nx[0],pG->Nx[1],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

  if ((Bcoeff = (Real**)calloc_2d_array(pG->Nx[0],pG->Nx[1],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");


/* To compute kx,ky,kz, note that indices relative to whole Domain are needed */  
    dkx = 2.0*PI/(double)(pD->Nx[0]);
    dky = 2.0*PI/(double)(pD->Nx[1]);
/* This is size of whole Domain perpendicular to the plane (=disk thickness)*/
    xmin = pD->RootMinX[1];
    xmax = pD->RootMaxX[1];
    Lperp = xmax-xmin;

/* Compute potential coeffs in k space. Zero wavenumber is special
   case; need to avoid divide by zero */
    for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      ip=KCOMP(i-is,pG->Disp[0],pD->Nx[0]);
      if (((j-js)+pG->Disp[1])==0 && ((i-is)+pG->Disp[0])==0) 
        Acoeff[0][0]=0.0;
      else{
        Acoeff[i-is][j-js]= 0.5*(1.0-exp(-fabs(ip*dkx/pG->dx1)*Lperp))/ 
	  (((2.0*cos(((i-is)+pG->Disp[0])*dkx)-2.0)/dx1sq) + 
	   ((2.0*cos(((j-js)+pG->Disp[1])*dky)-2.0)/dx2sq));
      }
      Bcoeff[i-is][j-js]= 0.5*(1.0+exp(-fabs(ip*dkx/pG->dx1)*Lperp))/ 
        (((2.0*cos((    (i-is)+pG->Disp[0])*dkx)-2.0)/dx1sq) + 
	 ((2.0*cos((0.5+(j-js)+pG->Disp[1])*dky)-2.0)/dx2sq));
    }
    }
  coeff_set=1; /* done computing coeffs */
  }
  
/* Copy current potential into old */

  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[ks][j][i] = pG->Phi[ks][j][i];
    }
  }

/* Fill complex work arrays with 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      /* real part */
      work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] = four_pi_G*pG->U[ks][j][i].d;
      /* imaginary part */
      work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] = 0.0;
      /* real part */
      work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0]= four_pi_G*pG->U[ks][j][i].d*
          cos(0.5*((j-js)+pG->Disp[1])*dky) ;
      /* imaginary part */
      work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1]= four_pi_G*pG->U[ks][j][i].d*
         -sin(0.5*((j-js)+pG->Disp[1])*dky);
    }
  }

/* Forward FFT of 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */

  ath_2d_fft(fplan2d, work);
  ath_2d_fft(fplan2d, work2);

/* Compute potential in Fourier space, using pre-computed coefficients. */

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] *=Acoeff[i-is][j-js]; 
      work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] *=Acoeff[i-is][j-js]; 
      work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] *=Bcoeff[i-is][j-js]; 
      work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] *=Bcoeff[i-is][j-js]; 
    }
  }

  /* Backward FFT */ 

  ath_2d_fft(bplan2d, work);
  ath_2d_fft(bplan2d, work2);

  /* Set potential in real space.  Normalization of Phi is over
      total number of cells in Domain */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pG->Phi[ks][j][i] = (work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0]
        +cos(0.5*((j-js)+pG->Disp[1])*dky)*work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] 
        -sin(0.5*((j-js)+pG->Disp[1])*dky)*work2[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1])/
        bplan2d->gcnt;
    }
  }

  return;
}



/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_disk_3d(DomainS *pD)
 *  \brief Periodic boundary conditions in x1 and x2; open bc in x3
 */

void selfg_fft_disk_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int ip, jp;
  Real kxtdx,kydy;
  Real dkx,dky,dkz;
  Real dx1sq=(pG->dx1*pG->dx1),dx2sq=(pG->dx2*pG->dx2),dx3sq=(pG->dx3*pG->dx3); 
  Real xmin,xmax;
  Real Lperp,den; 
  Real ***Acoeff=NULL,***Bcoeff=NULL; 

#ifdef SHEARING_BOX
  int nx3=pG->Nx[2]+2*nghost;
  int nx2=pG->Nx[1]+2*nghost;
  int nx1=pG->Nx[0]+2*nghost;
  Real ***RollDen=NULL, ***UnRollPhi=NULL;
  Real Lx,Ly,qomt,dt;

  if((RollDen=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");
  if((UnRollPhi=(Real***)calloc_3d_array(nx3,nx1,nx2,sizeof(Real)))==NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

  xmin = pD->RootMinX[0];
  xmax = pD->RootMaxX[0];
  Lx = xmax - xmin;

  xmin = pD->RootMinX[1];
  xmax = pD->RootMaxX[1];
  Ly = xmax - xmin;

  dt = pG->time-((int)(qshear*Omega_0*pG->time*Lx/Ly))*Ly/(qshear*Omega_0*Lx);
  qomt = qshear*Omega_0*dt;
#endif

/* allocates memory for Acoeff and Bcoeff arrays */
  if ((Acoeff = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

  if ((Bcoeff = (Real***)calloc_3d_array(pG->Nx[0],pG->Nx[1],pG->Nx[2],sizeof(Real))) == NULL)
    ath_error("[selfg_fft_disk]: malloc returned a NULL pointer\n");

/* To compute kx,ky,kz, note indices relative to whole Domain are needed */
  dkx = 2.0*PI/(double)(pD->Nx[0]);
  dky = 2.0*PI/(double)(pD->Nx[1]);
  dkz = 2.0*PI/(double)(pD->Nx[2]);

/* This is size of whole Domain perpendicular to the plane (=disk thickness)*/
  xmin = pD->RootMinX[2];
  xmax = pD->RootMaxX[2];
  Lperp = xmax-xmin;

/* Compute potential coeffs in k space. Zero wavenumber is special
   case; need to avoid divide by zero */
  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        ip=KCOMP(i-is,pG->Disp[0],pD->Nx[0]);
        jp=KCOMP(j-js,pG->Disp[1],pD->Nx[1]);
#ifdef SHEARING_BOX
        kxtdx = (ip+qomt*Lx/Ly*jp)*dkx;
#else
        kxtdx = ip*dkx;
#endif
        kydy = jp*dky;
	if (((k-ks)+pG->Disp[2])==0 && ((j-js)+pG->Disp[1])==0 && ((i-is)+pG->Disp[0])==0) 
          Acoeff[0][0][0] = 0.0;
        else{
          Acoeff[i-is][j-js][k-ks] = 0.5*  
            (1.0-exp(-sqrt(SQR(kxtdx)/dx1sq+SQR(kydy)/dx2sq)*Lperp))/ 
            (((2.0*cos(  kxtdx                 )-2.0)/dx1sq) + 
	     ((2.0*cos(  kydy                  )-2.0)/dx2sq) +
             ((2.0*cos(((k-ks)+pG->Disp[2])*dkz)-2.0)/dx3sq));
	}
        Bcoeff[i-is][j-js][k-ks] = 0.5*  
          (1.0+exp(-sqrt(SQR(kxtdx)/dx1sq+SQR(kydy)/dx2sq)*Lperp))/ 
          (((2.0*cos(      kxtdx                 )-2.0)/dx1sq) + 
	   ((2.0*cos(      kydy                  )-2.0)/dx2sq) +
           ((2.0*cos((0.5+(k-ks)+pG->Disp[2])*dkz)-2.0)/dx3sq));
      }
    }
  }

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
#ifdef SHEARING_BOX
        RollDen[k][i][j] = pG->U[k][j][i].d;
/* should add star particle density to RollDen using assign_starparticles_3d(pD,work), where work is the 1D
version of grid.  Note that assign_starparticles_3d only fills active zones.  Does RemapVar really need
the ghost zones? */
#endif
      }
    }
  }

#ifdef SHEARING_BOX
  RemapVar(pD,RollDen,-dt);
#endif

/* Fill arrays of 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */


  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        den=RollDen[k][i][j];
#else
        den=pG->U[k][j][i].d;
#endif
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = den;
      }
    }
  }

#ifndef SHEARING_BOX 
#ifdef STAR_PARTICLE
   assign_starparticles_3d(pD,work); 
#endif
#endif

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=four_pi_G;
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;

        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
          cos(0.5*((k-ks)+pG->Disp[2])*dkz)*
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
         -sin(0.5*((k-ks)+pG->Disp[2])*dkz)*
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];
      }
    }
  }

/*
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        den=RollDen[k][i][j]-grav_mean_rho;
#else
        den=pG->U[k][j][i].d-grav_mean_rho;
#endif
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
          four_pi_G*den;
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;

        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 
          four_pi_G*den*cos(0.5*((k-ks)+pG->Disp[2])*dkz);
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 
         -four_pi_G*den*sin(0.5*((k-ks)+pG->Disp[2])*dkz);
      }
    }
  }
*/

/* Forward FFT of 4\piG*d and 4\piG*d *exp(-i pi x2/Lperp) */

  ath_3d_fft(fplan3d, work);
  ath_3d_fft(fplan3d, work2);

/* Compute potential in Fourier space, using pre-computed coefficients */

  for (i=is; i<=ie; i++){
    for (j=js; j<=je; j++){
      for (k=ks; k<=ke; k++){
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=
          Acoeff[i-is][j-js][k-ks]; 
        work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *=
          Acoeff[i-is][j-js][k-ks]; 
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=
          Bcoeff[i-is][j-js][k-ks]; 
        work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *=
          Bcoeff[i-is][j-js][k-ks]; 
      }
    }
  }


/* Backward FFT and set potential in real space.  Normalization of Phi is over
 * total number of cells in Domain */

  ath_3d_fft(bplan3d, work);
  ath_3d_fft(bplan3d, work2);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
#ifdef SHEARING_BOX
        UnRollPhi[k][i][j] = 
#else
        pG->Phi[k][j][i] =
#endif
                           (work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
                         + cos(0.5*((k-ks)+pG->Disp[2])*dkz)*
		           work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] 
                         - sin(0.5*((k-ks)+pG->Disp[2])*dkz)*
		           work2[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1])/
                           bplan3d->gcnt;
      }
    }
  }

#ifdef SHEARING_BOX
  RemapVar(pD,UnRollPhi,dt);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
         pG->Phi[k][j][i] = UnRollPhi[k][i][j];
      }
    }
  }

  free_3d_array(RollDen);
  free_3d_array(UnRollPhi);
#endif
  free_3d_array(Acoeff);
  free_3d_array(Bcoeff);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_disk_2d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *  needed by FFTW.  
 */

void selfg_fft_disk_2d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_2d_fft_malloc(fplan2d);
        work2 = ath_2d_fft_malloc(fplan2d);
      }
    }
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_disk_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *  needed by FFTW.
 */

void selfg_fft_disk_3d_init(MeshS *pM)
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
        work2 = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* SELF_GRAVITY_USING_FFT_DISK */
