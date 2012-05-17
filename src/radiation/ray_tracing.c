#include "../copyright.h"
/*==============================================================================
 * FILE: ray_tracing.c
 *
 * PURPOSE: Contatins functions for handling ray tracing on Cartesian,
 * grid-alinged rays.  Models attenuation of "external" radiation that can
 * be approximated by parallel rays (e.g. a distant point source).  When
 * scattering is present, computes source term for the "diffuse" emission.
 * Current implementation assumes rays aligned with x1 direction with source
 * of radiation at ix1 boundary.  Currently assumes nu and nu_rt grids are
 * equivalent.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ray_trace()  - compute external radiation
 *   init_ray_tracing()  -memory allocation for ray tracing module
 *   destruct_ray_tracing() -free up memmory used by ray tracing module
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define TAUMIN 1.0E-6

#ifdef RAY_TRACING


#ifdef MPI_PARALLEL
/* MPI send and receive buffers */
static double *send_buf = NULL, *recv_buf = NULL;
static  int buf_size;
#endif
static Real ****S_rt = NULL;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   attenuate()     - computes attenuation of external radiation
 *============================================================================*/
void attenuate(DomainS *pD);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
void ray_trace(DomainS *pD)
{
  RadGridS *pRG = (pD->RadGrid);
  int nf = pRG->nf_rt, ifr;
  int il = pRG->is-1;
  int is = pRG->is, ie = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int i, j, k;
#ifdef MPI_PARALLEL
  int cnt, ierr, mIndex;
  double *pSnd,*pRcv;
  MPI_Status status;
#endif

  if (pRG->Nx[1] > 1) {
    jl--; ju++; 
  }
  if (pRG->Nx[2] > 1) {
    kl--; ku++;
  }

#ifdef MPI_PARALLEL
  if(pRG->lx1_id >= 0) {
/* set LHS values to 1.0 */
    for(ifr=0; ifr<nf; ifr++) {
      for(k=kl; k<=ku; k++) { 
	for(j=jl; j<=ju; j++) {
	  pRG->H[ifr][k][j][il] = 1.0;
	}}}
  }

/* compute the attenuation across the grid */
   attenuate(pD);

/* if left most grid send data immediately */  
   if (pRG->lx1_id < 0) {
     if (pRG->rx1_id >= 0) {
/* pack flux into send buffer */
       pSnd = send_buf;
       for(ifr=0; ifr<nf; ifr++) {
	 for (k=kl; k<=ku; k++) {
	   for (j=jl; j<=ju; j++) {
	     *(pSnd++) = pRG->H[ifr][k][j][ie];  
	   }}}
/* send data to right grid */
       ierr = MPI_Send(send_buf,buf_size,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
		       pD->Comm_Domain);
     }
   } else {
   
/* receive data from left grid */
     ierr = MPI_Recv(recv_buf,buf_size,MPI_DOUBLE,pRG->lx1_id,LtoR_tag,
		     pD->Comm_Domain,&status);

/* unpack flux from receive buffer */
     pRcv = recv_buf;
     for(ifr=0; ifr<nf; ifr++) {
       for (k=kl; k<=ku; k++) {
	 for (j=jl; j<=ju; j++) {
	   pRG->H[ifr][k][j][il] = *(pRcv++);
	 }}}

/* if grid on right, pack and send data */ 
     if (pRG->rx1_id >= 0) {

/* rescale right boundary flux and pack into send buffer */
       pSnd = send_buf;
       for(ifr=0; ifr<nf; ifr++) {
	 for(k=kl; k<=ku; k++) { 
	   for(j=jl; j<=ju; j++) {
	     *(pSnd++) = pRG->H[ifr][k][j][ie] * pRG->H[ifr][k][j][il];
	 }}}

/* send data to right grid */
       ierr = MPI_Send(send_buf,buf_size,MPI_DOUBLE,pRG->rx1_id,LtoR_tag,
		       pD->Comm_Domain);
     }

/* rescale flux throughout the volume */
     for(ifr=0; ifr<nf; ifr++) {
       for(k=kl; k<=ku; k++) { 
	 for(j=jl; j<=ju; j++) {
	   for(i=is; i<=ie; i++) {
	     pRG->H[ifr][k][j][i] *= pRG->H[ifr][k][j][il];
	     S_rt[ifr][k][j][i] *= pRG->H[ifr][k][j][il];
	   }
	   S_rt[ifr][k][j][il] *= pRG->H[ifr][k][j][il];
	   S_rt[ifr][k][j][ie+1] *= pRG->H[ifr][k][j][il];
	 }}}
   }
#else
/* Only need to compute the attenuation across the grid for serial computation */
   attenuate(pD);
#endif

/* Update S and S_nt with scattered component */
   for(ifr=0; ifr<nf; ifr++) {
     for(k=kl; k<=ku; k++) { 
       for(j=jl; j<=ju; j++) {
	 for(i=is; i<=ie; i++) {
	   pRG->R[ifr][k][j][i].S += S_rt[ifr][k][j][i];
	   pRG->R[ifr][k][j][i].Snt = S_rt[ifr][k][j][i];
	 }}}}

  return;
}

/* ------------------------------------------------------------ 
 * computes attenuation of external irradiation
 */
void attenuate(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG = (pD->RadGrid);
  int i, j, k;
  int is = pRG->is, ie = pRG->ie; 
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf_rt, ifr;
  int jg,kg,ioff,joff,koff;
  Real dx = pRG->dx1;
  Real chit, chitp, chitm, dtau, edtau;
  Real dH, eps;

  if (pG->Nx[0] > 1) {
    ioff = nghost - 1;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1; jl--; ju++; 
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1; kl--; ku++;
  } else koff = 0;

  for(ifr=0; ifr<nf; ifr++) 
    for(k=kl; k<=ku; k++) { 
      kg = k + koff;
      for(j=jl; j<=ju; j++) {
	jg = j + joff;
	chitm = get_raytrace_opacity(pG,ifr,is+ioff-1,jg,kg);
	chit = get_raytrace_opacity(pG,ifr,is+ioff,jg,kg);
	for(i=is; i<=ie; i++) {
	  chitp = get_raytrace_opacity(pG,ifr,i+ioff+1,jg,kg);
	  interp_quad_chi(chitm,chit,chitp,&dtau);
	  dtau *= dx;
	  if (dtau < TAUMIN) {
	    dH = (1.0 - 0.5 * dtau);
	    edtau = 1.0 - dH;
	  } else {
	    edtau = exp(-dtau);
	    dH = (1.0 - edtau) / dtau;	      
	  }
	  eps = get_raytrace_thermal_fraction(pG,ifr,i+ioff+1,jg,kg);
	  S_rt[ifr][k][j][i] = dH * (1.0 - eps) * pRG->H[ifr][k][j][i-1];
	  pRG->H[ifr][k][j][i] = pRG->H[ifr][k][j][i-1] * edtau;
	  chitm = chit;
	  chit = chitp;	  
	}
	S_rt[ifr][k][j][is-1] = S_rt[ifr][k][j][is];
	S_rt[ifr][k][j][ie+1] = S_rt[ifr][k][j][ie];
      }
    }

  return;
}

/* ------------------------------------------------------------
 * initalizes communication buffers and allocates memory for 
 * working arrays used by ray_tracing algorithms 
 */
void init_ray_tracing(RadGridS *pRG)
{
  int nx1 = pRG->Nx[0], nx2 = pRG->Nx[1], nx3 = pRG->Nx[2];
  int nf = pRG->nf;

#ifdef MPI_PARALLEL
/* allocate memory for send/receive buffers */
  buf_size = (pRG->Nx[1] + 2) * (pRG->Nx[2] + 2) * pRG->nf_rt;

  if((send_buf = (double*)calloc_1d_array(buf_size,sizeof(double))) == NULL)
    ath_error("[init_ray_tracing]: Failed to allocate send buffer\n");
  if((recv_buf = (double*)calloc_1d_array(buf_size,sizeof(double))) == NULL)
    ath_error("[init_ray_tracing]: Failed to allocate recv buffer\n");
#endif

/* allocate memory for source function array */
  if ((S_rt = (Real ****)calloc_4d_array(nf,nx3+2,nx2+2,nx1+2,sizeof(Real))) == NULL)
    ath_error("[init_ray_tracing]: error allocating memory for S_rt\n");

  return;

}

/* ------------------------------------------------------------
 * deallocates memory for arrays used by ray_tracing
 */
void destruct_ray_tracing(void)
{
#ifdef MPI_PARALLEL
  if (send_buf != NULL) free(send_buf);
  if (recv_buf != NULL) free(recv_buf);
#endif
  if (S_rt   != NULL) free_4d_array(S_rt);

  return;
}

#endif /* RAY_TRACING */
    
