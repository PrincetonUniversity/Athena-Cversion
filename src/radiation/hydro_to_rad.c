#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_rad.c
 *
 * PURPOSE:  Contains functions for updating the RadGrid using the conserved
 *           variables in Grid (hydro_to_rad) and for computing
 *           the radiation source term and updating the material energy in
 *           Grid (rad_to_hydro).
 *
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   hydro_to_rad()
 *   rad_to_hydro()
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void hydro_to_rad(DomainS *pD, const int outflag)
 *  Computes opacities, thermalization paraters, thermal source function
 *  from primative variables using user defined functions */
void hydro_to_rad(DomainS *pD, const int outflag)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG;
  int i,j,k,ifr;
  int il, iu, jl, ju, kl, ku;
  int nf;
  int ig,jg,kg,ioff,joff,koff;
  Real eps;
  Real etherm, ekin, B, d;

  if (outflag == 0)
    pRG = pD->RadGrid;    /* set ptr to RadGrid */
  else
    pRG = pD->RadOutGrid; /* set ptr to RadOutGrid */
  il = pRG->is; iu = pRG->ie;
  jl = pRG->js; ju = pRG->je;
  kl = pRG->ks; ku = pRG->ke;
  nf = pRG->nf;


/* Assumes ghost zone conserved variables have been set by
 * bvals routines.  These values are used to set B, chi, eps,
 * etc. so loops include RadGrid ghost zones*/
  if (pG->Nx[0] > 1) {
    ioff = nghost - 1; il--; iu++;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1; jl--; ju++;
  } else joff = 0;
  if (pG->Nx[2] > 1) {
    koff = nghost - 1; kl--; ku++;
  } else koff = 0;

/* Compute radiation variables from conserved variables */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
        ig = i + ioff;

        /* Compute gas temperature and store for later use */
        d = pG->U[kg][jg][ig].d;
        etherm = pG->U[kg][jg][ig].E - (0.5/d) * ( SQR(pG->U[kg][jg][ig].M1) +
                 SQR(pG->U[kg][jg][ig].M2) + SQR(pG->U[kg][jg][ig].M3) );
#if defined(MHD) || defined(RADIATION_MHD)
        etherm -= 0.5 * (SQR(pG->U[kg][jg][ig].B1c)+SQR(pG->U[kg][jg][ig].B2c)+
                         SQR(pG->U[kg][jg][ig].B3c));
#endif
        if (outflag == 0)
          pG->tgas[kg][jg][ig] = MAX(etherm * Gamma_1 / (d * R_ideal),0.0);

        for(ifr=0; ifr<nf; ifr++) {
          eps = get_thermal_fraction(pG,ifr,ig,jg,kg);
          pRG->R[ifr][k][j][i].B = get_thermal_source(pG,ifr,ig,jg,kg);
          pRG->R[ifr][k][j][i].eps = eps;
          pRG->R[ifr][k][j][i].S = (1.0 - eps) * pRG->R[ifr][k][j][i].J +
                                          eps  * pRG->R[ifr][k][j][i].B;
          pRG->R[ifr][k][j][i].chi = get_total_opacity(pG,ifr,ig,jg,kg);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void rad_to_hydro(DomainS *pD)
 *  Performs operator split update of the energy equation using radiation
 *  transfer solution directly.  Uses differential form for high optical
 *  depth and integral form for low optical depth */
void rad_to_hydro(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k, ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int nDim;
  int ig,jg,kg,ioff,joff,koff;
  Real esource, kappa;
  Real dx1=0.5/pRG->dx1, dx2=0.5/pRG->dx2, dx3=0.5/pRG->dx3;
  Real dxmin;
  int flag = 0;
  static Real dt = 1.0e-3;

  dxmin = pD->dx[0];
  if (pD->Nx[1] > 1) dxmin = MIN( dxmin, (pD->dx[1]) );
  if (pD->Nx[2] > 1) dxmin = MIN( dxmin, (pD->dx[2]) );

  ioff = nghost - 1;
  nDim = 1;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1;
    nDim = 2;
  } else joff = 0;
  if (pG->Nx[2] > 1) {
    koff = nghost - 1;
    nDim = 3;
  } else koff = 0;

/* Update thermal energy */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
        ig = i + ioff;
        esource = 0.0;
        for(ifr=0; ifr<nf; ifr++) {
          if(pRG->R[ifr][k][j][i].chi*dxmin <= 1.0) {

            esource += pRG->wnu[ifr] * pRG->R[ifr][k][j][i].eps *
                       pRG->R[ifr][k][j][i].chi * (pRG->R[ifr][k][j][i].J -
                       pRG->R[ifr][k][j][i].B);
          } else {
            esource += pRG->wnu[ifr] * dx1 * (pRG->R[ifr][k][j][i-1].H[0] -
                                              pRG->R[ifr][k][j][i+1].H[0]);
            if (nDim > 1) {
              esource += pRG->wnu[ifr] * dx2 * (pRG->R[ifr][k][j-1][i].H[1] -
                                                pRG->R[ifr][k][j+1][i].H[1]);
              if (nDim == 3) {
                esource += pRG->wnu[ifr] * dx3 * (pRG->R[ifr][k-1][j][i].H[2] -
                                                  pRG->R[ifr][k+1][j][i].H[2]);
              }
            }
          }

        }
        pG->U[kg][jg][ig].E += pG->dt * 4.0 * PI * esource * CPrat;
        /*pG->U[kg][jg][ig].E += dt * 4.0 * PI * esource * CPrat;*/
      }}}


  return;
}

#endif /* RADIATION_TRANSFER */
