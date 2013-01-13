#include "copyright.h"
/*============================================================================*/
/*! \file fft_test.c
 *  \brief Problem generator for testing FFT rountines with SMR and MPI
 *
 *  - Forward and backward FFT is done for input (Gaussian) profile.
 *  - nlim should be set 1 (one-step calculation).
 *  - Use input file tst/3D-hydro/athinput.fft_test
 *  - OUTPUT: each variable doesn't represent its own physical meaning. 
 *      input data -> d
 *      real part of FFT(d) -> M1
 *      imaginary part of FFT(d) -> M2
 *      IFFT(FFT(d)) -> M3
 *
 *  - Last updated Mar 29, 2012
 *
 *  WRITTEN BY Chang-Goo Kim 						      */
/*============================================================================*/

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "prototypes.h"
#include "globals.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif /* MPI_PARALLEL */

#ifndef FFT_ENABLED
#error FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

static ath_fft_data *work=NULL;
void do_fft(DomainS *pD);

void problem(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is=pG->is, ie = pG->ie;
  int j, js=pG->js, je = pG->je;
  int k, ks=pG->ks, ke = pG->ke;
  Real r2,x1,x2,x3;

  if(pD->Nx[2] > 1) {
/* 3D case */
    pD->fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
    pD->bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
  } else if (pD->Nx[1] > 1) {
/* 2D case */
    pD->fplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
    pD->bplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
  } else {
/* 1D case */
  }

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        r2=SQR(x1)+SQR(x2)+SQR(x3);
        pG->U[k][j][i].d = exp(-r2);
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M2 = 0.0;
        pG->U[k][j][i].M3 = 0.0;
      }
    }
  }

  return;
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

/* problem restart needs to enroll cooling function and boundary conditions. 
 */

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef RESISTIVITY
void get_eta_user(GridS *pG, int i, int j, int k,
                             Real *eta_O, Real *eta_H, Real *eta_A)
{

  *eta_O = 0.0;
  *eta_H = 0.0;
  *eta_A = 0.0;

  return;
}
#endif

void Userwork_in_loop(MeshS *pM)
{
  int nl,nd;
  int cntmax=0;
  int me, nproc;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        if(pM->Domain[nl][nd].fplan2d != NULL) cntmax=MAX(cntmax,pM->Domain[nl][nd].fplan2d->cnt);
        if(pM->Domain[nl][nd].fplan3d != NULL) cntmax=MAX(cntmax,pM->Domain[nl][nd].fplan3d->cnt);
#ifdef MPI_PARALLEL
        MPI_Comm_rank(pM->Domain[nl][nd].Comm_Domain,&me);
        MPI_Comm_size(pM->Domain[nl][nd].Comm_Domain,&nproc);
        printf("nl=%d nd=%d rank=%d nproc=%d\n\n",nl,nd,me,nproc);
#endif
      }
    }
  }

  work = (ath_fft_data *)fftw_malloc(sizeof(ath_fft_data) * cntmax);

  printf("id=%d work array is allocated with cnt=%d!\n",myID_Comm_world,cntmax);

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        printf("id=%d nl=%d nd=%d do FFT\n",myID_Comm_world,nl,nd);
        do_fft(&(pM->Domain[nl][nd]));
        printf("id=%d FFT is done.\n",myID_Comm_world);
      }
    }
  }

  data_output(pM, 0);

  return;
}

void Userwork_after_loop(MeshS *pM)
{
  int nl, nd;

  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL){
        if(pM->Nx[2] > 1) {
/* 3D case */
          ath_3d_fft_destroy_plan(pM->Domain[nl][nd].fplan3d);
          ath_3d_fft_destroy_plan(pM->Domain[nl][nd].bplan3d);
          if (work != NULL) ath_3d_fft_free(work);
        } else if (pM->Nx[1] > 1) {
/* 2D case */
          ath_2d_fft_destroy_plan(pM->Domain[nl][nd].fplan2d);
          ath_2d_fft_destroy_plan(pM->Domain[nl][nd].bplan2d);
          if (work != NULL) ath_2d_fft_free(work);
        } else {
/* 1D case */
        }
      }
    }
  }

  return;
}

void do_fft(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i, is=pG->is, ie = pG->ie;
  int j, js=pG->js, je = pG->je;
  int k, ks=pG->ks, ke = pG->ke;
  Real r2,x1,x2,x3;


  if(pD->Nx[2] > 1) {
/* 3D case */

    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          r2=SQR(x1)+SQR(x2)+SQR(x3);
          pG->U[k][j][i].d = exp(-r2);
          work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = pG->U[k][j][i].d;
        }
      }
    }

    ath_3d_fft(pD->fplan3d, work);
   
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pG->U[k][j][i].M1 = work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];
          pG->U[k][j][i].M2 = work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1];
        }
      }
    }

    ath_3d_fft(pD->bplan3d, work);

    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pG->U[k][j][i].M3 = work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]/pD->bplan3d->gcnt;
        }
      }
    }

  } else if (pD->Nx[1] > 1) {
/* 2D case */

    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          r2=SQR(x1)+SQR(x2);
          pG->U[k][j][i].d = exp(-r2);
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] = pG->U[k][j][i].d;
        }
      }
    }

    ath_2d_fft(pD->fplan2d, work);
   
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pG->U[k][j][i].M1 = work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0];
          pG->U[k][j][i].M2 = work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1];
        }
      }
    }

    ath_2d_fft(pD->bplan2d, work);

    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          pG->U[k][j][i].M3 = work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0]/pD->bplan2d->gcnt;
        }
      }
    }

  } else {
/* 1D case */
  }

  return;
}
