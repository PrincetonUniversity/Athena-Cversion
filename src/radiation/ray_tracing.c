#include "../copyright.h"
/*==============================================================================
 * FILE: ray_tracing.c
 *
 * PURPOSE: Contatins functions for handling ray tracing on Cartesian,
 * grid-alinged rays.  Models attenuation of "external" radiation that can
 * be approximated by parallel rays (e.g. a distant point source).  When
 * scattering is present, computes source term for the "diffuse" emission.
 * Current implementation assumes rays aligned with x1 direction with source
 * of radiation at ix1 boundary.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   ray_trace()  - compute external radiation
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define TAUMIN 1.0E-3

#ifdef RAY_TRACING

#ifdef MPI_PARALLEL
static Real ***H0 = NULL;
#endif

#ifdef MPI_PARALLEL
void allocate_working_array_ray_tracing(RayGridS *pRayG);
void destruct_working_array_ray_tracing(void);
#endif
void attenuate(DomainS *pD);

void ray_trace(DomainS *pD)
{
  
#if MPI_PARRALLEL
  /* store values of H on lhs and reset to 1.0 */

#endif
  attenuate(pD);
  return;
}

void attenuate(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RayGridS *pRayG=(pD->RayGrid);
  int i, j, k;
  int is = pRayG->is, ie = pRayG->ie; 
  int js = pRayG->js, je = pRayG->je;
  int ks = pRayG->ks, ke = pRayG->ke;
  int nf = pRayG->nf, ifr;
  int jg,kg,ioff,joff,koff;
  Real dx = pRayG->dx1;
  Real chit, chitp, chitm, dtau, edtau;
  Real dH;

  if (pG->Nx[0] > 1) {
    ioff = nghost - 1;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1;
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1;
  } else koff = 0;

  for(ifr=0; ifr<nf; ifr++) 
    for(k=ks; k<=ke; k++) { 
      kg = k + koff;
      for(j=js; j<=je; j++) {
	jg = j + joff;
	chitm = get_raytrace_opacity(pG,ifr,is+ioff-1,jg,kg);
	chit = get_raytrace_opacity(pG,ifr,is+ioff,jg,kg);
	for(i=is; i<=ie; i++) {
	  chitp = get_raytrace_opacity(pG,ifr,i+ioff+1,jg,kg);
	  interp_quad_chi(chitm,chit,chitp,&dtau);
	  dtau *= dx;
	  if (dtau < TAUMIN) {
	    dH = SQR(dtau) * (1.0 - 0.5 * dtau);
	    edtau = 1.0 - dH;
	  } else {
	    edtau = exp(-dtau);
	    dH = 1.0 - edtau;	      
	  }
	  dH *= pRayG->H[ifr][k][j][i];
	  pRayG->H[ifr][k][j][i] = pRayG->H[ifr][k][j][i-1] * edtau;
	  // if (j == 128) {
	  //  printf("H: %d %g %g %g\n",i,pRayG->H[ifr][k][j][i],dtau,edtau);
	  // }
	  chitm = chit;
	  chit = chitp;	  
	}
      }
    }

  return;
}

#ifdef MPI_PARALLEL
/* allocates memory for working arrays used by ray_tracing */
void allocate_working_array_ray_tracing(RayGridS *pRayG)
{
  int nx2 = pRayG->Nx[1], nx3 = pRayG->Nx[2];
  int nf = pRayG->nf;
  
  if ((H0 = (Real ***)calloc_3d_array(nf,nx3,nx2,sizeof(Real))) == NULL)
    goto on_error;

  return;

  on_error:
  formal_solution_2d_destruct();
  ath_error("[allocate_working_array_ray_tracing]: Error allocating memory\n");
  return;
}

/* allocates memory for working arrays used by ray_tracing */
void destruct_working_array_ray_tracing(void)
{

   if (H0 != NULL) free_3d_array(H0);
  return;
}
#endif /* MPI_PARALLEL */


#endif /* RAY_TRACING */
    
