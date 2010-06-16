#include "../copyright.h"
/*==============================================================================
 * FILE: init_radiation.c
 *
 * PURPOSE: Initializes most variables in the RadGrid structure.  Modelled
 *          after init_grid.c.  Also contains radiation_destruct for
 *          memory deallocation.
 *
 *          For the moment, it does not work with SMR.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_radiation()
 *   radiation_temp_array_init
 *   radiation_desrtuct()
 *   radgrid_destruct()
 *============================================================================*/

#include <math.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION

/*----------------------------------------------------------------------------*/
/* init_radiation:  */

void init_radiation(MeshS *pM)
{
  DomainS *pD;
  RadGridS *pRG;
  int nDim,nl,nd,myL,myM,myN;
  int i,l,m,n;

/* number of dimensions in Grid. */
  nDim=1;
  for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */
      pRG = pM->Domain[nl][nd].RadGrid;          /* set ptr to RadGrid */

/* Initialize nf,nmu,ng */
      pRG->nf = par_geti("problem","nf");
      pRG->nmu = par_geti("problem","nmu");
      pRG->ng = par_geti("problem","ng");

/* get (l,m,n) coordinates of Grid being updated on this processor */

      get_myGridIndex(pD, myID_Comm_world, &myL, &myM, &myN);
 

/* ---------------------  Intialize grid in 1-direction --------------------- */
/* Initialize is,ie,dx1
 * Compute Disp, MinX[0], and MaxX[0] using displacement of Domain and Grid
 * location within Domain */

      pRG->Nx[0] = pD->GData[myN][myM][myL].Nx[0];

      if(pRG->Nx[0] > 1) {
        pRG->is = 1;
        pRG->ie = pRG->Nx[0];
      }
      else
        pRG->is = pRG->ie = 0;

      pRG->dx1 = pD->dx[0];
    
      pRG->Disp[0] = pD->Disp[0];
      pRG->MinX[0] = pD->MinX[0];
      for (l=1; l<=myL; l++) {
        pRG->Disp[0] +=        pD->GData[myN][myM][l-1].Nx[0];
        pRG->MinX[0] += (Real)(pD->GData[myN][myM][l-1].Nx[0])*pRG->dx1;
      }
      pRG->MaxX[0] = pRG->MinX[0] + (Real)(pRG->Nx[0])*pRG->dx1;
    
/* ---------------------  Intialize grid in 2-direction --------------------- */
/* Initialize js,je,dx2
 * Compute Disp, MinX[1], and MaxX[1] using displacement of Domain and Grid
 * location within Domain */

      pRG->Nx[1] = pD->GData[myN][myM][myL].Nx[1];
    
      if(pRG->Nx[1] > 1) {
        pRG->js = 1;
        pRG->je = pRG->Nx[1];
      }
      else
        pRG->js = pRG->je = 0;

      pRG->dx2 = pD->dx[1];

      pRG->Disp[1] = pD->Disp[1];
      pRG->MinX[1] = pD->MinX[1];
      for (m=1; m<=myM; m++) {
        pRG->Disp[1] +=        pD->GData[myN][m-1][myL].Nx[1];
        pRG->MinX[1] += (Real)(pD->GData[myN][m-1][myL].Nx[1])*pRG->dx2;
      }
      pRG->MaxX[1] = pRG->MinX[1] + (Real)(pRG->Nx[1])*pRG->dx2;

/* ---------------------  Intialize grid in 3-direction --------------------- */
/* Initialize ks,ke,dx3
 * Compute Disp, MinX[2], and MaxX[2] using displacement of Domain and Grid
 * location within Domain */

      pRG->Nx[2] = pD->GData[myN][myM][myL].Nx[2];

      if(pRG->Nx[2] > 1) {
        pRG->ks = 1;
        pRG->ke = pRG->Nx[2];
      }
      else
        pRG->ks = pRG->ke = 0;

      pRG->dx3 = pD->dx[2];

      pRG->Disp[2] = pD->Disp[2];
      pRG->MinX[2] = pD->MinX[2];
      for (n=1; n<=myN; n++) {
        pRG->Disp[2] +=        pD->GData[n-1][myM][myL].Nx[2];
        pRG->MinX[2] += (Real)(pD->GData[n-1][myM][myL].Nx[2])*pRG->dx3;
      }
      pRG->MaxX[2] = pRG->MinX[2] + (Real)(pRG->Nx[2])*pRG->dx3;
      

/*  Allocate memory for array of RadS */
      pRG->R = (RadS ****)calloc_4d_array(pRG->Nx[2]+2,pRG->Nx[1]+2,
        pRG->Nx[0]+2,pRG->nf,sizeof(RadS));
      if (pRG->R == NULL) goto on_error1;
      
/* Allocate memory for angles and weights for angular quadratures */

      pRG->mu = (Real *)calloc_1d_array(pRG->nmu,sizeof(Real));
      if (pRG->mu ==NULL) goto on_error2;

      pRG->w = (Real **)calloc_2d_array(pRG->nmu,pRG->ng,sizeof(Real));
      if (pRG->w ==NULL) goto on_error3;

      if(nDim > 1) {
	pRG->gamma = (Real *)calloc_1d_array(pRG->ng,sizeof(Real));
	if (pRG->gamma ==NULL) goto on_error4;
      }

/* Allocate memory for intensity at boundaries */

      if (pRG->Nx[0] > 1) {
	pRG->r1imu = (Real *****)calloc_5d_array(pRG->Nx[2]+2,pRG->Nx[1]+2,
          pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->r1imu == NULL) goto on_error5;

	pRG->l1imu = (Real *****)calloc_5d_array(pRG->Nx[2]+2,pRG->Nx[1]+2,
	  pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->l1imu == NULL) goto on_error6;
      }

      if (pRG->Nx[1] > 1) {
	pRG->r2imu = (Real *****)calloc_5d_array(pRG->Nx[2]+2,pRG->Nx[0]+2,
          pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->r2imu == NULL) goto on_error7;

	pRG->l2imu = (Real *****)calloc_5d_array(pRG->Nx[2]+2,pRG->Nx[0]+2,
          pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->l2imu == NULL) goto on_error8;
      }

      if (pRG->Nx[2] > 1) {
	pRG->r3imu = (Real *****)calloc_5d_array(pRG->Nx[1]+2,pRG->Nx[0]+2,
          pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->r3imu == NULL) goto on_error9;

	pRG->l3imu = (Real *****)calloc_5d_array(pRG->Nx[1]+2,pRG->Nx[0]+2,
          pRG->nf,pRG->nmu,pRG->ng,sizeof(Real));
	if (pRG->l3imu == NULL) goto on_error10;
      }

/*-- Get IDs of neighboring Grids in Domain communicator ---------------------*/
/* If Grid is at the edge of the Domain (so it is either a physical boundary,
 * or an internal boundary between fine/coarse grids), then ID is set to -1
 */

/* Left-x1 */
      if(myL > 0) pRG->lx1_id = pD->GData[myN][myM][myL-1].ID_Comm_Domain;
      else pRG->lx1_id = -1;

/* Right-x1 */
      if(myL <(pD->NGrid[0])-1)
        pRG->rx1_id = pD->GData[myN][myM][myL+1].ID_Comm_Domain;
      else pRG->rx1_id = -1;

/* Left-x2 */
      if(myM > 0) pRG->lx2_id = pD->GData[myN][myM-1][myL].ID_Comm_Domain;
      else pRG->lx2_id = -1;

/* Right-x2 */
      if(myM <(pD->NGrid[1])-1)
        pRG->rx2_id = pD->GData[myN][myM+1][myL].ID_Comm_Domain;
      else pRG->rx2_id = -1;

/* Left-x3 */
      if(myN > 0) pRG->lx3_id = pD->GData[myN-1][myM][myL].ID_Comm_Domain;
      else pRG->lx3_id = -1;

/* Right-x3 */
      if(myN <(pD->NGrid[2])-1)
        pRG->rx3_id = pD->GData[myN+1][myM][myL].ID_Comm_Domain;
      else pRG->rx3_id = -1;

#ifdef RAD_MULTIG
/* set maximum number of multigrid refinements */
  if(nDim == 1)
    nmgrid = pRG->Nx[0] / 16;
  else if(nDim == 2)
    nmgrid = pRG->Nx[1] / 16;
#endif


    }
  }
  }

  return;

/*--- Error messages ---------------------------------------------------------*/

 on_error10:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->l3imu);
 on_error9:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->r3imu);
 on_error8:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->l2imu);
 on_error7:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->r2imu);
 on_error6:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->l1imu);
 on_error5:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->r1imu);
 on_error4:
  if(nDim > 1) free_1d_array(pRG->gamma);
 on_error3:
  free_2d_array(pRG->w);
 on_error2:
  free_1d_array(pRG->mu);
 on_error1:
  free_4d_array(pRG->R);
  ath_error("[init_radiation]: Error allocating memory\n");

}

void radiation_temp_array_init(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  int i, dim;

/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pRG->Nx[i] > 1) dim++;

/* set function pointer to appropriate integrator based on dimensions */
  /*switch(dim){

  case 1:
    formal_solution_1d_alloc(pRG);
  case 2:
    formal_solution_2d_alloc(pRG);
    }*/

  return;
}

void radiation_destruct(MeshS *pM)
{
 
  RadGridS *pRG;
  int nl,nd;

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pRG = pM->Domain[nl][nd].RadGrid;          /* set ptr to RadGrid */
      radgrid_destruct(pRG);
    }
  }}
  return;
}

void radgrid_destruct(RadGridS *pRG)
{
 
  if (pRG->R != NULL) free_4d_array(pRG->R);  
  if (pRG->w != NULL) free_2d_array(pRG->w);
  if (pRG->mu != NULL) free(pRG->mu);
  if (pRG->gamma != NULL) free(pRG->gamma);
  if (pRG->r3imu != NULL) free_5d_array(pRG->r3imu);
  if (pRG->l3imu != NULL) free_5d_array(pRG->l3imu);
  if (pRG->r2imu != NULL) free_5d_array(pRG->r2imu);
  if (pRG->l2imu != NULL) free_5d_array(pRG->l2imu);
  if (pRG->r1imu != NULL) free_5d_array(pRG->r1imu);
  if (pRG->l1imu != NULL) free_5d_array(pRG->l1imu);

  return;
}

#endif /* RADIATION */
