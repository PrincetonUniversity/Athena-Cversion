#include "../copyright.h"
/*==============================================================================
 * FILE: integrate_3d_radMHD.c
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   unsplit CTU integrator of Colella (1990).  The variables updated are:
 *      U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *      B1i, B2i, B3i  -- interface magnetic field
 *   Also adds gravitational source terms, self-gravity, optically thin cooling,
 *   shearing box source terms, and the H-correction of Sanders et al.
 *     For adb hydro, requires (9*Cons1DS +  3*Real) = 48 3D arrays
 *     For adb mhd, requires   (9*Cons1DS + 10*Real) = 73 3D arrays
 *   The H-correction of Sanders et al. adds another 3 arrays.  
 *
 * REFERENCES:
 *   P. Colella, "Multidimensional upwind methods for hyperbolic conservation
 *   laws", JCP, 87, 171 (1990)
 *
 *   T. Gardiner & J.M. Stone, "An unsplit Godunov method for ideal MHD via
 *   constrained transport in three dimensions", JCP, 227, 4123 (2008)
 *
 *   R. Sanders, E. Morano, & M.-C. Druguet, "Multidimensinal dissipation for
 *   upwind schemes: stability and applications to gas dynamics", JCP, 145, 511
 *   (1998)
 *
 *   J.M. Stone et al., "Athena: A new code for astrophysical MHD", ApJS,
 *   178, 137 (2008)
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   integrate_3d_radMHD()
 *   integrate_init_3d()
 *   integrate_destruct_3d()
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
static Cons1DS ***Ul_x1Face=NULL, ***Ur_x1Face=NULL;
static Cons1DS ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
static Cons1DS ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL;

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real ***dhalf = NULL, ***phalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real ***eta1=NULL, ***eta2=NULL, ***eta3=NULL;
#endif

/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static Real **remapEyiib=NULL, **remapEyoib=NULL;
#endif

/* VARIABLES NEED FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
static Real ***geom_src=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - the upwind CT method in GS05, for emf1
 *   integrate_emf2_corner() - the upwind CT method in GS05, for emf2
 *   integrate_emf3_corner() - the upwind CT method in GS05, for emf3
 *============================================================================*/

#ifdef MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif /* MHD */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* integrate_3d: 3D CTU integrator for MHD using 6-solve method */

void integrate_3d_radMHD(DomainS *pD)
{
  
  return;
}

/*----------------------------------------------------------------------------*/
/* integrate_init_3d: Allocate temporary integration arrays 
*/

void integrate_init_3d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

#ifdef MHD
  if ((emf1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;

  if ((emf1_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf2_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((emf3_cc=(Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* MHD */
#ifdef H_CORRECTION
  if ((eta1 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta2 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((eta3 = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
#endif /* H_CORRECTION */

  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

#ifdef MHD
  if ((B1_x1Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B2_x2Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
  if ((B3_x3Face = (Real***)calloc_3d_array(size3,size2,size1, sizeof(Real)))
    == NULL) goto on_error;
#endif /* MHD */

  if ((U1d=(Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wl =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr =(Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Ul_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x1Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x2Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ul_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((Ur_x3Face=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x1Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

#ifdef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((phalf = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  }

#ifdef SHEARING_BOX
  if ((remapEyiib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
    goto on_error;
  if ((remapEyoib = (Real**)calloc_2d_array(size3,size2, sizeof(Real))) == NULL)
    goto on_error;
#endif

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if ((geom_src = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
    goto on_error;
#endif

  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/* integrate_destruct_3d:  Free temporary integration arrays 
 */

void integrate_destruct_3d(void)
{

#ifdef MHD
  if (emf1    != NULL) free_3d_array(emf1);
  if (emf2    != NULL) free_3d_array(emf2);
  if (emf3    != NULL) free_3d_array(emf3);
  if (emf1_cc != NULL) free_3d_array(emf1_cc);
  if (emf2_cc != NULL) free_3d_array(emf2_cc);
  if (emf3_cc != NULL) free_3d_array(emf3_cc);
#endif /* MHD */
#ifdef H_CORRECTION
  if (eta1 != NULL) free_3d_array(eta1);
  if (eta2 != NULL) free_3d_array(eta2);
  if (eta3 != NULL) free_3d_array(eta3);
#endif /* H_CORRECTION */

  if (Bxc != NULL) free(Bxc);
  if (Bxi != NULL) free(Bxi);
#ifdef MHD
  if (B1_x1Face != NULL) free_3d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_3d_array(B2_x2Face);
  if (B3_x3Face != NULL) free_3d_array(B3_x3Face);
#endif /* MHD */

  if (U1d      != NULL) free(U1d);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Ul_x1Face != NULL) free_3d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_3d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_3d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_3d_array(Ur_x2Face);
  if (Ul_x3Face != NULL) free_3d_array(Ul_x3Face);
  if (Ur_x3Face != NULL) free_3d_array(Ur_x3Face);
  if (x1Flux    != NULL) free_3d_array(x1Flux);
  if (x2Flux    != NULL) free_3d_array(x2Flux);
  if (x3Flux    != NULL) free_3d_array(x3Flux);
  if (dhalf     != NULL) free_3d_array(dhalf);
  if (phalf     != NULL) free_3d_array(phalf);
#ifdef SHEARING_BOX
  if (remapEyiib != NULL) free_2d_array(remapEyiib);
  if (remapEyoib != NULL) free_2d_array(remapEyoib);
#endif

  /* DATA STRUCTURES FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_3d_array(geom_src);
#endif

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/* integrate_emf1_corner
 * integrate_emf2_corner
 * integrate_emf3_corner
 *   Integrates face centered B-fluxes to compute corner EMFs.  Note:
 *   x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 *   x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 *   x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 *   x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 *   x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 *   x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */

#ifdef MHD
static void integrate_emf1_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de1_l2, de1_r2, de1_l3, de1_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-2; i<=ie+2; i++) {
/* NOTE: The x2-Flux of By is -E1. */
/*       The x3-Flux of Bz is +E1. */
	if (x2Flux[k-1][j][i].d > 0.0)
	  de1_l3 = x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i];
	else if (x2Flux[k-1][j][i].d < 0.0)
	  de1_l3 = x3Flux[k][j][i].Bz - emf1_cc[k-1][j][i];
	else {
	  de1_l3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k-1][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k-1][j  ][i] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de1_r3 = x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de1_r3 = x3Flux[k][j][i].Bz - emf1_cc[k][j][i];
	else {
	  de1_r3 = 0.5*(x3Flux[k][j-1][i].Bz - emf1_cc[k][j-1][i] +
			x3Flux[k][j  ][i].Bz - emf1_cc[k][j  ][i] );
	}

	if (x3Flux[k][j-1][i].d > 0.0)
	  de1_l2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i];
	else if (x3Flux[k][j-1][i].d < 0.0)
	  de1_l2 = -x2Flux[k][j][i].By - emf1_cc[k][j-1][i];
	else {
	  de1_l2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j-1][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j-1][i] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de1_r2 = -x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de1_r2 = -x2Flux[k][j][i].By - emf1_cc[k][j][i];
	else {
	  de1_r2 = 0.5*(-x2Flux[k-1][j][i].By - emf1_cc[k-1][j][i]
			-x2Flux[k  ][j][i].By - emf1_cc[k  ][j][i] );
	}

        emf1[k][j][i] = 0.25*(  x3Flux[k][j][i].Bz + x3Flux[k][j-1][i].Bz
                              - x2Flux[k][j][i].By - x2Flux[k-1][j][i].By 
			      + de1_l2 + de1_r2 + de1_l3 + de1_r3);
      }
    }
  }

  return;
}

static void integrate_emf2_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de2_l1, de2_r1, de2_l3, de2_r3;

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-2; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
/* NOTE: The x1-Flux of Bz is +E2. */
/*       The x3-Flux of By is -E2. */
	if (x1Flux[k-1][j][i].d > 0.0)
	  de2_l3 = -x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1];
	else if (x1Flux[k-1][j][i].d < 0.0)
	  de2_l3 = -x3Flux[k][j][i].By - emf2_cc[k-1][j][i];
	else {
	  de2_l3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k-1][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k-1][j][i  ] );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de2_r3 = -x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1];
	else if (x1Flux[k][j][i].d < 0.0)
	  de2_r3 = -x3Flux[k][j][i].By - emf2_cc[k][j][i];
	else {
	  de2_r3 = 0.5*(-x3Flux[k][j][i-1].By - emf2_cc[k][j][i-1] 
			-x3Flux[k][j][i  ].By - emf2_cc[k][j][i  ] );
	}

	if (x3Flux[k][j][i-1].d > 0.0)
	  de2_l1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1];
	else if (x3Flux[k][j][i-1].d < 0.0)
	  de2_l1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i-1];
	else {
	  de2_l1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i-1] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i-1] );
	}

	if (x3Flux[k][j][i].d > 0.0)
	  de2_r1 = x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i];
	else if (x3Flux[k][j][i].d < 0.0)
	  de2_r1 = x1Flux[k][j][i].Bz - emf2_cc[k][j][i];
	else {
	  de2_r1 = 0.5*(x1Flux[k-1][j][i].Bz - emf2_cc[k-1][j][i] +
			x1Flux[k  ][j][i].Bz - emf2_cc[k  ][j][i] );
	}

	emf2[k][j][i] = 0.25*(  x1Flux[k][j][i].Bz + x1Flux[k-1][j][i  ].Bz
                              - x3Flux[k][j][i].By - x3Flux[k  ][j][i-1].By
			      + de2_l1 + de2_r1 + de2_l3 + de2_r3);
      }
    }
  }

  return;
}

static void integrate_emf3_corner(const GridS *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real de3_l1, de3_r1, de3_l2, de3_r2;
  Real rsf=1.0,lsf=1.0;

  for (k=ks-2; k<=ke+2; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+2; i++) {
/* NOTE: The x1-Flux of By is -E3. */
/*       The x2-Flux of Bx is +E3. */
#ifdef CYLINDRICAL
        rsf = pG->ri[i]/pG->r[i];  lsf = pG->ri[i]/pG->r[i-1];
#endif
	if (x1Flux[k][j-1][i].d > 0.0)
	  de3_l2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf;
	else if (x1Flux[k][j-1][i].d < 0.0)
	  de3_l2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j-1][i])*rsf;
	else {
	  de3_l2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j-1][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j-1][i  ])*rsf );
	}

	if (x1Flux[k][j][i].d > 0.0)
	  de3_r2 = (x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf;
	else if (x1Flux[k][j][i].d < 0.0)
	  de3_r2 = (x2Flux[k][j][i].Bz - emf3_cc[k][j][i])*rsf;
	else {
	  de3_r2 = 0.5*((x2Flux[k][j][i-1].Bz - emf3_cc[k][j][i-1])*lsf + 
			(x2Flux[k][j][i  ].Bz - emf3_cc[k][j][i  ])*rsf );
	}

	if (x2Flux[k][j][i-1].d > 0.0)
	  de3_l1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1];
	else if (x2Flux[k][j][i-1].d < 0.0)
	  de3_l1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i-1];
	else {
	  de3_l1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i-1]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i-1] );
	}

	if (x2Flux[k][j][i].d > 0.0)
	  de3_r1 = -x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i];
	else if (x2Flux[k][j][i].d < 0.0)
	  de3_r1 = -x1Flux[k][j][i].By - emf3_cc[k][j][i];
	else {
	  de3_r1 = 0.5*(-x1Flux[k][j-1][i].By - emf3_cc[k][j-1][i]
			-x1Flux[k][j  ][i].By - emf3_cc[k][j  ][i] );
	}

	emf3[k][j][i] = 0.25*(  x2Flux[k][j  ][i-1].Bz + x2Flux[k][j][i].Bz
			      - x1Flux[k][j-1][i  ].By - x1Flux[k][j][i].By
			      + de3_l1 + de3_r1 + de3_l2 + de3_r2);
      }
    }
  }

  return;
}
#endif /* MHD */

#endif /* CTU_INTEGRATOR */
