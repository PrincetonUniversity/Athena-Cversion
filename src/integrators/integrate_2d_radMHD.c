#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_2d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations in 2D using the directionally unsplit CTU
 *   method of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *      B1i, B2i -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically-thin cooling,
 *   and H-correction of Sanders et al.
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport", JCP, 205, 509 (2005)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_2d_radMHD()
 *   integrate_init_2d()
 *   integrate_destruct_2d()
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(radMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS **Ul_x1Face=NULL, **Ur_x1Face=NULL;
static Cons1DS **Ul_x2Face=NULL, **Ur_x2Face=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL, *Ul=NULL, *Ur=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real **dhalf = NULL,**phalf = NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/* VARIABLES NEEDED FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
static Real **geom_src=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf3_corner() - the upwind CT method in Gardiner & Stone (2005) 
 *============================================================================*/

#ifdef MHD
static void integrate_emf3_corner(GridS *pG);
#endif

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_2d:  CTU integrator in 2D.
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */

void integrate_2d_radMHD(DomainS *pD)
{
  

  return;

 
}

/*----------------------------------------------------------------------------*/
/* integrate_init_2d:    Allocate temporary integration arrays */

void integrate_init_2d(MeshS *pM)
{
  int nmax,size1=0,size2=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  nmax = MAX(size1,size2);

#ifdef MHD
  if ((emf3 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((emf3_cc = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((eta2 = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* H_CORRECTION */

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef MHD
  if ((B1_x1Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((B2_x2Face = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif /* MHD */

  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ul = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Ur = (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x1Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ul_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((Ur_x2Face=(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

#ifndef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  }

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if ((geom_src = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL) 
    goto on_error;
#endif

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}


/*----------------------------------------------------------------------------*/
/* integrate_destruct_2d:  Free temporary integration arrays */

void integrate_destruct_2d(void)
{
#ifdef MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_2d_array(eta1);
  if (eta2 != NULL) free_2d_array(eta2);
#endif /* H_CORRECTION */
  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_2d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_2d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_2d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_2d_array(Ur_x2Face);
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (dhalf     != NULL) free_2d_array(dhalf);
  if (phalf     != NULL) free_2d_array(phalf);

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_2d_array(geom_src);
#endif

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf3_corner:  */

#ifdef MHD
static void integrate_emf3_corner(GridS *pG)
{
  int i,is,ie,j,js,je;
  Real emf_l1, emf_r1, emf_l2, emf_r2;
  Real rsf=1.0,lsf=1.0;

  is = pG->is;   ie = pG->ie;
  js = pG->js;   je = pG->je;

/* NOTE: The x1-Flux of B2 is -E3.  The x2-Flux of B1 is +E3. */
  for (j=js-1; j<=je+2; j++) {
    for (i=is-1; i<=ie+2; i++) {
#ifdef CYLINDRICAL
      rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
      if (x1Flux[j-1][i].d > 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf;
      }
      else if (x1Flux[j-1][i].d < 0.0) {
	emf_l2 = -x1Flux[j-1][i].By
	  + (x2Flux[j][i].Bz - emf3_cc[j-1][i])*rsf;

      } else {
	emf_l2 = -x1Flux[j-1][i].By
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j-1][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j-1][i  ])*rsf );
      }

      if (x1Flux[j][i].d > 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf;
      }
      else if (x1Flux[j][i].d < 0.0) {
	emf_r2 = -x1Flux[j][i].By 
	  + (x2Flux[j][i].Bz - emf3_cc[j][i])*rsf;

      } else {
	emf_r2 = -x1Flux[j][i].By 
	  + 0.5*((x2Flux[j][i-1].Bz - emf3_cc[j][i-1])*lsf + 
		 (x2Flux[j][i  ].Bz - emf3_cc[j][i  ])*rsf );
      }

      if (x2Flux[j][i-1].d > 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]);
      }
      else if (x2Flux[j][i-1].d < 0.0) {
	emf_l1 = x2Flux[j][i-1].Bz 
          + (-x1Flux[j][i].By - emf3_cc[j][i-1]);
      } else {
	emf_l1 = x2Flux[j][i-1].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i-1]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i-1] );
      }

      if (x2Flux[j][i].d > 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j-1][i].By - emf3_cc[j-1][i]);
      }
      else if (x2Flux[j][i].d < 0.0) {
	emf_r1 = x2Flux[j][i].Bz
	  + (-x1Flux[j][i].By - emf3_cc[j][i]);
      } else {
	emf_r1 = x2Flux[j][i].Bz
	  + 0.5*(-x1Flux[j-1][i].By - emf3_cc[j-1][i]
		 -x1Flux[j  ][i].By - emf3_cc[j  ][i] );
      }

      emf3[j][i] = 0.25*(emf_l1 + emf_r1 + emf_l2 + emf_r2);
    }
  }

  return;
}
#endif /* MHD */

#endif /* radMHD_INTEGRATOR */
