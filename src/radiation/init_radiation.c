#include "../copyright.h"
/*==============================================================================
 * FILE: init_radiation.c
 *
 * PURPOSE: Initializes most variables in the RadGrid structure.  Modelled
 *          after init_grid.c.  Also contains radiation_destruct for
 *          memory deallocation.
 *
 *          It is not compatible with SMR.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   init_radiation()
 *   radiation_desrtuct()
 *============================================================================*/

#include <math.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"

#ifdef RADIATION_TRANSFER

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * init_radiation_grid() - Allocate memory and intialize RadGrid
 * radgrid_destruct() - Free up memory allocated to RadGrid structures
 * formal_solution_init() - Calls appropriate formal_solution_*d_init()
 * formal_solution_destruct() -  Calls appropriate formal_solution_*d_destruct()
 * init_angles() - Initialize angles/quadratures with default method
 * InverseMatrix() -  Compute matrix inverse
 * MatrixMult() - Matrix multiplication
 * permutation() - Checks for permutations of previously assigned vector
 * gauleg() -  Gauss-Legendre quadrature for Numerical Recipes
 * ludcmp_nr() - LU decomposition from Numerical Recipes
 * lubksb_nr() - Backward substitution from Numerical Recipies
 *============================================================================*/

void init_radiation_grid(DomainS *pD, int outflag);
void radgrid_destruct(RadGridS *pRG);
void formal_solution_init(DomainS *pD);
void formal_solution_destruct(DomainS *pD);
void init_angles(RadGridS *pRG, const int nmu);
void InverseMatrix(Real **a, int n, Real **b);
void MatrixMult(Real **a, Real *b, int m, int n, Real *c);
int permutation(int i, int j, int k, int **pl, int np);
void gauleg(Real x1, Real x2,  Real *x, Real *w, int n);
void ludcmp_nr(Real **a, int n, int *indx, Real *d);
void lubksb_nr(Real **a, int n, int *indx, Real b[]);

/*=========================== PUBLIC FUNCTIONS ===============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void init_radiation(MeshS *pM)
 *  \brief Call routines to intialize RadGrid/RadOutGrid */
void init_radiation(MeshS *pM)
{
  DomainS *pD;
  RadGridS *pRG;
  int nl,nd;

/* initialize global variable lte which controls whether J or S is tested
   for convergence */
  lte = par_geti("radiation","lte");
/* Initialize global variable CPrat to 1 */
  CPrat = 1.0;

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
      pD = (DomainS*)&(pM->Domain[nl][nd]);  /* set ptr to Domain */

/* Initialize RadGrid used for integration of the hydrodynamics */
      if (radt_mode == 0) { /* integration only */
	init_radiation_grid(pD,0);
	pD->RadOutGrid = NULL;
      } else if (radt_mode == 1) { /* output only */
	pD->RadGrid = NULL;
	init_radiation_grid(pD,1);
      } else if (radt_mode == 2) { /* integration and output */
	init_radiation_grid(pD,0);
	init_radiation_grid(pD,1);
      } else {
	ath_error("[init_radiation]: radiation mode must be 0,1, or 2 but is set to %d\n");
      }

/* Initialize working arrays for formal solution */
      formal_solution_init(pD);
    }
  }
  }

  return;

}

/*----------------------------------------------------------------------------*/
/*! \fn void radiation_destruct(MeshS *pM)
 *  \brief Free memory allocated to RadGrid/RadOutGrid */
void radiation_destruct(MeshS *pM)
{
 
  RadGridS *pRG;
  DomainS *pD;
  int nl,nd;

  for (nl=0; nl<pM->NLevels; nl++){
  for (nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
    if (pM->Domain[nl][nd].Grid != NULL) {
/* Destruct RadGrid used for integration of the hydrodynamics */
      if ( (radt_mode == 0) || (radt_mode == 2) ) {
	pRG = pM->Domain[nl][nd].RadGrid;  /* set ptr to RadGrid */
	radgrid_destruct(pRG);
      }
/* Destruct RadGrid used for output diagnostics */
     if ( (radt_mode == 1) || (radt_mode == 2) ) {
       pRG = pM->Domain[nl][nd].RadOutGrid; /* set ptr to RadGrid */
       radgrid_destruct(pRG);
     }
/* Call subroutine to deallocate working arrays in formal solution */
      pD = (DomainS*)&(pM->Domain[nl][nd]);
      formal_solution_destruct(pD);
    }
  }}
  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn void init_radiation_grid(DomainS *pD, int outflag)
 *  \brief Allocates memory and intializes RadGrid structure.
 *  calls init_angles(). */
void init_radiation_grid(DomainS *pD, int outflag)
{
/* Allocates memory and sets up angular grid for pRG */

  RadGridS *pRG;
  int nDim,myL,myM,myN;;
  int i,j,k,l,m,n;
  int nmu;

  if (outflag == 0) {
/* -------------------------- Initialize RadGrid array  -------------------------- */
      pRG = pD->RadGrid;    /* set ptr to RadGrid */
      pRG->outflag = 0;
      pRG->pG = pD->Grid;  /* set correspond Grids pointer */
      pRG->nf  = par_geti("radiation","nf");
      nmu = par_geti_def("radiation","nmu",0);
      if (nmu == 0) {
/* if nmu == 0, nang must be set in input block and angles initalized in problem generator */
	pRG->nang = par_geti_def("radiation","nang",0);
	if (pRG->nang == 0) 
	  ath_error("[init_radiation]: neither nang or nmu are set in radiation block.\n");
      }
#ifdef RAY_TRACING
      pRG->nf_rt = par_geti_def("radiation","nf_rt",pRG->nf);
#endif
  } else {
/* -------------------------- Initialize RadOutGrid array  -------------------------- */  
      pRG = pD->RadOutGrid;          /* set ptr to RadOutGrid */
      pRG->outflag = 1;
      pRG->pG = pD->Grid;  /* set correspond Grids pointer */
      pRG->nf  = par_geti("radiation_output","nf");
      nmu = par_geti_def("radiation_output","nmu",0);
      if (nmu == 0) {
/* if nmu == 0, nang must be set in input block and angles initalized in problem generator */
	pRG->nang = par_geti_def("radiation_output","nang",0);
	if (pRG->nang == 0) 
	  ath_error("[init_radiation]: neither nang or nmu are set in radiation output block.\n");
      }
#ifdef RAY_TRACING
      pRG->nf_rt = par_geti_def("radiation_output","nf_rt",pRG->nf);
#endif
  }

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

/* number of dimensions in RadGrid. */
  nDim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

/* Initialize nf, nang, and noct members of RadGrid */

  if (nDim == 1) {
    pRG->noct = 2;
  } else  if (nDim == 2) {
    pRG->noct = 4;
  } else if (nDim == 3) {
    pRG->noct = 8;
  }
  if (nmu != 0) {
    if (nDim == 1) {
      pRG->nang = nmu;
    } else  if (nDim == 2) {
      pRG->nang = nmu * (nmu + 1) / 2;
    } else if (nDim == 3) {
      pRG->nang = nmu * (nmu + 1) / 2;
    }
  }

/* -------------------- Allocate memory for RadGrid arrays --------------------*/
/*  Allocate memory for array of RadS */
  pRG->R = (RadS ****)calloc_4d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[1]+2,
				      pRG->Nx[0]+2,sizeof(RadS));
  if (pRG->R == NULL) goto on_error1;
/* set non-thermal source function to zero everywhere */
  for(l=0; l<pRG->nf; l++) {
    for(k=0; k<pRG->Nx[2]+2; k++) { 
      for(j=0; j<pRG->Nx[1]+2; j++) {
	for(i=0; i<pRG->Nx[0]+2; i++) {
	  pRG->R[l][k][j][i].Snt = 0.0;
	}}}}
/* Allocate memory for intensities, angles and weights for angular quadratures */

  pRG->mu = (Real ***)calloc_3d_array(pRG->noct,pRG->nang,3,sizeof(Real));
  if (pRG->mu == NULL) goto on_error2;

  pRG->wmu = (Real *)calloc_1d_array(pRG->nang,sizeof(Real));
  if (pRG->wmu == NULL) goto on_error3;     

 
/* Allocate memory for intensity at boundaries */ 
  if (pRG->Nx[0] > 1) {
    pRG->r1imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[1]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->r1imu == NULL) goto on_error4;

    pRG->l1imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[1]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->l1imu == NULL) goto on_error5;
  } else {
    pRG->r1imu = NULL;
    pRG->l1imu = NULL;	
  }

  if (pRG->Nx[1] > 1) {
    pRG->r2imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[0]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->r2imu == NULL) goto on_error6;
    pRG->l2imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[0]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->l2imu == NULL) goto on_error7;
  } else {
    pRG->r2imu = NULL;
    pRG->l2imu = NULL;	
  }

  if (pRG->Nx[2] > 1) {
    pRG->r3imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[1]+2,pRG->Nx[0]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->r3imu == NULL) goto on_error8;
    
    pRG->l3imu = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[1]+2,pRG->Nx[0]+2,
					     pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->l3imu == NULL) goto on_error9;
  } else {
    pRG->r3imu = NULL;
    pRG->l3imu = NULL;	
  }

/* Allocate memory for frequency and quadrature arrays */ 
  pRG->nu = (Real *)calloc_1d_array(pRG->nf,sizeof(Real));
  if (pRG->nu == NULL) goto on_error10;

  pRG->wnu = (Real *)calloc_1d_array(pRG->nf,sizeof(Real));
  if (pRG->wnu == NULL) goto on_error11;
/* initialize wnu to unity if nf=1 */
  if (pRG->nf == 1) pRG->wnu[0] = 1.0;

/* Allocate memory for intensity Ghost Zones */ 
  if (pRG->Nx[0] > 1) {
    pRG->Ghstr1i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[1]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstr1i == NULL) goto on_error12;
    
    pRG->Ghstl1i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[1]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstl1i == NULL) goto on_error13;
  } else {
    pRG->Ghstr1i = NULL;
    pRG->Ghstl1i = NULL;	
  }

  if (pRG->Nx[1] > 1) {
    pRG->Ghstr2i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[0]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstr2i == NULL) goto on_error14;
    pRG->Ghstl2i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[2]+2,pRG->Nx[0]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstl2i == NULL) goto on_error15;
  } else {
    pRG->Ghstr2i = NULL;
    pRG->Ghstl2i = NULL;	
  }

  if (pRG->Nx[2] > 1) {
    pRG->Ghstr3i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[1]+2,pRG->Nx[0]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstr3i == NULL) goto on_error16;
    
    pRG->Ghstl3i = (Real *****)calloc_5d_array(pRG->nf,pRG->Nx[1]+2,pRG->Nx[0]+2,
					       pRG->noct,pRG->nang,sizeof(Real));
    if (pRG->Ghstl3i == NULL) goto on_error17;
  } else {
    pRG->Ghstr3i = NULL;
	pRG->Ghstl3i = NULL;	
  }
  
  if (nmu != 0) {
/* Call init_angles() to initialize default angle grid and angular quadrature
 * If nmu = 0, angles/quadratures must be initialized in problem generator */
    init_angles(pRG,nmu);
  }

#ifdef RAY_TRACING
/* allocate memory for ray tracing arrays */
      init_ray_tracing(pRG);
#endif 

  return;

/*--- Error messages ---------------------------------------------------------*/

 on_error17:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->Ghstl3i);
 on_error16:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->Ghstr3i);
 on_error15:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->Ghstl2i);
 on_error14:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->Ghstr2i);
 on_error13:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->Ghstl1i);
 on_error12:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->Ghstr1i);
 on_error11:
  free_1d_array(pRG->wnu);  
 on_error10:
  free_1d_array(pRG->nu);
 on_error9:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->l3imu);
 on_error8:
  if (pRG->Nx[2] > 1) free_5d_array(pRG->r3imu);
 on_error7:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->l2imu);
 on_error6:
  if (pRG->Nx[1] > 1) free_5d_array(pRG->r2imu);
 on_error5:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->l1imu);
 on_error4:
  if (pRG->Nx[0] > 1) free_5d_array(pRG->r1imu);
 on_error3:
  free_1d_array(pRG->wmu);
 on_error2:
  free_3d_array(pRG->mu);
 on_error1:
  free_4d_array(pRG->R);
  ath_error("[init_radiation]: Error allocating memory\n");

}

/*----------------------------------------------------------------------------*/
/*! \fn void radgrid_destruct(RadGridS *pRG)
 *  \brief free up memory allocated to RadGrid structures  */
void radgrid_destruct(RadGridS *pRG)
{
  if (pRG->R != NULL) free_4d_array(pRG->R);  
  if (pRG->wmu != NULL) free_1d_array(pRG->wmu);
  if (pRG->mu != NULL) free_3d_array(pRG->mu);
  if (pRG->wnu != NULL) free_1d_array(pRG->wnu);
  if (pRG->nu != NULL) free_1d_array(pRG->nu);
  if (pRG->r3imu != NULL) free_5d_array(pRG->r3imu);
  if (pRG->l3imu != NULL) free_5d_array(pRG->l3imu);
  if (pRG->r2imu != NULL) free_5d_array(pRG->r2imu);
  if (pRG->l2imu != NULL) free_5d_array(pRG->l2imu);
  if (pRG->r1imu != NULL) free_5d_array(pRG->r1imu);
  if (pRG->l1imu != NULL) free_5d_array(pRG->l1imu);
  if (pRG->Ghstr3i != NULL) free_5d_array(pRG->Ghstr3i);
  if (pRG->Ghstl3i != NULL) free_5d_array(pRG->Ghstl3i);
  if (pRG->Ghstr2i != NULL) free_5d_array(pRG->Ghstr2i);
  if (pRG->Ghstl2i != NULL) free_5d_array(pRG->Ghstl2i);
  if (pRG->Ghstr1i != NULL) free_5d_array(pRG->Ghstr1i);
  if (pRG->Ghstl1i != NULL) free_5d_array(pRG->Ghstl1i);
#ifdef RAY_TRACING
      destruct_ray_tracing(pRG);
#endif /* RAY_TRACING */

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_init(DomainS *pD)
 *  \brief Call appropriate formal_solution_*d_init() */
void formal_solution_init(DomainS *pD)
{

  int i, dim;

/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pD->Nx[i] > 1) dim++;

/* set function pointer to appropriate initalization routine based on dimensions */
  if (dim == 1) 
    formal_solution_1d_init(pD);
  else if (dim == 2) 
    formal_solution_2d_init(pD);
  else if (dim == 3)
    formal_solution_3d_init(pD);
  else
    ath_error("[formal_solution_init]: incorrect number of dim: %d.\n",dim);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void formal_solution_destruct(DomainS *pD)
 *  \brief Call appropriate formal_solution_*d_destruct() */
void formal_solution_destruct(DomainS *pD)
{

  int i, dim;

/* Calculate the dimensions (using root Domain)  */
  dim = 0;
  for (i=0; i<3; i++) if(pD->Nx[i] > 1) dim++;

/* set function pointer to appropriate destruct routine based on dimensions */
  if (dim == 1) 
    formal_solution_1d_destruct();
  else if (dim == 2) 
    formal_solution_2d_destruct();
  else if (dim == 3)
    formal_solution_3d_destruct();
  else
    ath_error("[formal_solution_destruct]: incorrect number of dim: %d.\n",dim);

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void init_angles(RadGridS *pRG, const int nmu)
 *  \brief Initialize angles/angle quadratures using default methods */
void init_angles(RadGridS *pRG, const int nmu)
{
  int nDim;
  int np, ip, iang;
  int i,j,k,l,m;
  Real deltamu, Wsum, wsum, W2;
  Real *mu2tmp = NULL, **mutmp = NULL, *mutmp1d = NULL;
  Real *Wtmp = NULL, *wtmp = NULL;
  Real **pmat = NULL, **pinv = NULL, *wpf = NULL;
  int **pl = NULL, *plab = NULL;

/* number of dimensions in RadGrid. */
  nDim=1;
  for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;

 if(nDim == 1) {
/* 1D  compute angles using gaussian quadarature */
    mutmp1d = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
    if (mutmp1d == NULL) goto on_error1; 
    wtmp = (Real *)calloc_1d_array(2.0*nmu,sizeof(Real));
    if (wtmp == NULL) goto on_error2; 

    gauleg(-1.0, 1.0, mutmp1d, wtmp, 2*nmu);

    for(i=nmu; i<2*nmu; i++) {
      pRG->wmu[i-nmu] = 0.5 * wtmp[i];
      pRG->mu[0][i-nmu][0] = mutmp1d[i];
      pRG->mu[1][i-nmu][0] = -mutmp1d[i];
    }
    free_1d_array(wtmp);
    free_1d_array(mutmp1d);
  } else {
/* 2D and 3D:  compute angles and weights for angular quadratures following 
   the algorithm described in Bruls et al. 1999, A&A, 348, 233 */

    mu2tmp = (Real *)calloc_1d_array(nmu,sizeof(Real));
    if (mu2tmp == NULL) goto on_error3;
    mutmp = (Real **)calloc_2d_array(pRG->nang,3,sizeof(Real));
    if (mutmp == NULL) goto on_error4;
    Wtmp = (Real *)calloc_1d_array(nmu-1,sizeof(Real));
    if (Wtmp == NULL) goto on_error5;
    wtmp = (Real *)calloc_1d_array(nmu,sizeof(Real));
    if (wtmp == NULL) goto on_error6;

/* first compute polar weights and angles */
    if (nmu <= 6) {
      deltamu = 2.0 / (2 * nmu - 1);
      mu2tmp[0] = 1.0 / (3.0 * (2 * nmu - 1));
      for (i=1; i<nmu; i++) {
	mu2tmp[i] = mu2tmp[i-1] + deltamu;
      }
    } else {
      mu2tmp[0] = 1.0 / SQR((Real)nmu-1.0);
      deltamu = (1.0 - 3.0 * mu2tmp[0]) / ((Real)nmu-1.0);
      for (i=1; i<nmu; i++) {
	mu2tmp[i] = mu2tmp[i-1] + deltamu;
      }
    }
	
    W2 = 4.0 * mu2tmp[0];
    Wsum = Wtmp[0] = sqrt(W2);
    for (i=1; i<nmu-2; i++) {
      W2 += deltamu;
      Wsum += Wtmp[i] = sqrt(W2);
    }
    if (nmu > 2) Wtmp[nmu-2] = 2.0*(nmu-1)/3.0 - Wsum;
    
    wsum = wtmp[0] = Wtmp[0];
    for (i=1; i<nmu-1; i++) {
      wsum += wtmp[i] = Wtmp[i] - Wtmp[i-1];
    }    
    wtmp[nmu-1] = 1.0 - wsum;

/* Next, set up system of equations for determining how polar weights
   are distributed in azimuth (along circles of section), subject to
   the constraint that members of permutation families have identical
   weights */

    pmat = (Real **)calloc_2d_array(nmu,nmu,sizeof(Real));
    if (pmat == NULL) goto on_error7;
    pinv = (Real **)calloc_2d_array(nmu-1,nmu-1,sizeof(Real));
    if (pinv == NULL) goto on_error8;
    plab = (int *)calloc_1d_array(pRG->nang,sizeof(int));
    if (plab == NULL) goto on_error9;
    pl = (int **)calloc_2d_array(nmu,3,sizeof(int));
    if (pl == NULL) goto on_error10;
    wpf = (Real *)calloc_1d_array(nmu-1,sizeof(Real));
    if (wpf == NULL) goto on_error11;

    np = 0;
    iang = 0;
    for (i=0; i<nmu; i++) {
      for (j=0; j<nmu; j++) {
	for (k=0; k<nmu; k++) {
	  if (i + j + k == nmu - 1) {
/* assign cosines to temporary array grid */
	    mutmp[iang][0] = sqrt(mu2tmp[j]);
	    mutmp[iang][1] = sqrt(mu2tmp[k]);
	    mutmp[iang][2] = sqrt(mu2tmp[i]);
	    if (nmu <= 6) {
	      ip=permutation(i,j,k,pl,np); 
	      if (ip == -1) {
		pl[np][0] = i;
		pl[np][1] = j;
		pl[np][2] = k;		  
		pmat[i][np] += 1.0;
		plab[iang] = np;
		np++;		
	      } else {
		pmat[i][ip] += 1.0;
		plab[iang] = ip;
	      }
	    }
	    iang++;
	  }
	}
      }
    }

    if (nmu <= 6) {	  
/* Use Bruls/Carlsson formulation */
      if (nmu > 1) {
/*  Invert matrix of permutations families */
	InverseMatrix(pmat,nmu-1,pinv);
/* Solve for and assign weights for each permutation family */
	MatrixMult(pinv,wtmp,nmu-1,nmu-1,wpf);
	for (i=0; i<pRG->nang; i++) 
	  pRG->wmu[i] = wpf[plab[i]];
      } else 
	pRG->wmu[0] = 1.0;
    } else {
/* Use equal weights for all angles */
      for (i=0; i<pRG->nang; i++) 
	pRG->wmu[i] = 1.0/(Real)pRG->nang;
    }

/*  assign angles to RadGrid elements */
    if (nDim == 2) {
      for (i=0; i<pRG->nang; i++) {
	for (j=0; j<2; j++) {
	  for (k=0; k<2; k++) {
	    l=2*j+k;
	    if (k == 0)
	      pRG->mu[l][i][0] =  mutmp[i][0];
	    else
	      pRG->mu[l][i][0] = -mutmp[i][0];
	    if (j == 0)
	      pRG->mu[l][i][1] =  mutmp[i][2];
	    else
	      pRG->mu[l][i][1] = -mutmp[i][2];	    
	  }
	}
	pRG->wmu[i] *= 0.25;
      }
    } else if (nDim == 3) {
      for (i=0; i<pRG->nang; i++) {
	for (j=0; j<2; j++) {
	  for (k=0; k<2; k++) {
	    for (l=0; l<2; l++) {
	      m=4*j+2*k+l;
	      if (l == 0)
		pRG->mu[m][i][0] =  mutmp[i][0];
	      else
		pRG->mu[m][i][0] = -mutmp[i][0];
	      if (k == 0)
		pRG->mu[m][i][1] =  mutmp[i][1];
	      else
		pRG->mu[m][i][1] = -mutmp[i][1];
	      if (j == 0)
		pRG->mu[m][i][2] =  mutmp[i][2];
	      else
		pRG->mu[m][i][2] = -mutmp[i][2];	      
	    }
	  }
	}
	pRG->wmu[i] *= 0.125;	    
      }
    }
/* deallocate temporary arrays */
    free_1d_array(wpf);
    free_2d_array(pl);
    free_1d_array(plab);
    free_2d_array(pinv);
    free_2d_array(pmat);
    free_1d_array(wtmp);
    free_1d_array(Wtmp);
    free_2d_array(mutmp);
    free_1d_array(mu2tmp);
  } 
/* end of multidimensional quadratures */


  return;

/*--- Error messages ---------------------------------------------------------*/

 on_error11:
  if (nDim > 1) free_1d_array(wpf);
 on_error10:
  if (nDim > 1) free_2d_array(pl);
 on_error9:
  if (nDim > 1) free_1d_array(plab);
 on_error8:
  if (nDim > 1) free_2d_array(pinv);
 on_error7:
  if (nDim > 1) free_2d_array(pmat);
 on_error6:
  if (nDim > 1) free_1d_array(wtmp);
 on_error5:
  if (nDim > 1) free_1d_array(Wtmp);
 on_error4:
  if (nDim > 1) free_2d_array(mutmp);
 on_error3:
  if (nDim > 1) free_1d_array(mu2tmp);
 on_error2:
  if (nDim ==  1) free_1d_array(wtmp);
 on_error1:
  if (nDim ==  1) free_1d_array(mutmp1d);
  ath_error("[init_angles]: Error allocating memory\n");

}

/*----------------------------------------------------------------------------*/
/*! \fn int permutation(int i, int j, int k, int **pl, int np)
 * Checks if an indicies triplet is a permutation of vectors in pl[][]
 * (previously loaded triplets) and returns index of permutation if matched
 * or -1 if no match */
int permutation(int i, int j, int k, int **pl, int np)
{
  int ip=-1;
  int l,m,n,o;

/*  This routine is only called at initialization so brute force
 *  algorithm is fine
 */


  for(l=0; l<np; l++) {
/* check each permutation in the table */
    for(m=0; m<3; m++)
      if(i == pl[l][m])
	for(n=0; n<3; n++)
	  if(n != m)
	    if(j == pl[l][n])
	      for(o=0;o<3;o++)
		if((o != m) && (o != n))
		  if(k == pl[l][o]) 
		    ip = l;
  }

  return ip;
}

/*----------------------------------------------------------------------------*/
/*! \fn void ludcmp_nr(Real **a, int n, int *indx, Real *d)
 * LU decomposition from Numerical Recipes
 * Using Crout's method with partial pivoting
 * a is the input matrix, and is returned with LU decomposition readily made,
 * n is the matrix size, indx records the history of row permutation,
 * whereas d =1(-1) for even(odd) number of permutations.
 */
void ludcmp_nr(Real **a, int n, int *indx, Real *d)
{
  int i,imax,j,k;
  Real big,dum,sum,temp;
  Real *rowscale;  /* the implicit scaling of each row */

  rowscale = (Real*)calloc_1d_array(n, sizeof(Real));
  *d=1.0;  /* No row interchanges yet */

  for (i=0;i<n;i++)
  { /* Loop over rows to get the implicit scaling information */
    big=0.0;
    for (j=0;j<n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) ath_error("[LUdecomp]:Input matrix is singular!");
    rowscale[i]=1.0/big;  /* Save the scaling */
  }

  for (j=0;j<n;j++) { /* Loop over columns of Crout's method */
    /* Calculate the upper block */
    for (i=0;i<j;i++) {
      sum=a[i][j];
      for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    /* Calculate the lower block (first step) */
    big=0.0;
    for (i=j;i<n;i++) {
      sum=a[i][j];
      for (k=0;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      /* search for the largest pivot element */
      if ( (dum=rowscale[i]*fabs(sum)) >= big) {
        big=dum;
        imax=i;
      }
    }
    /* row interchange */
    if (j != imax) {
      for (k=0;k<n;k++) {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      rowscale[imax]=rowscale[j];
    }
    indx[j]=imax; /* record row interchange history */
    /* Calculate the lower block (second step) */
    if (a[j][j] == 0.0) a[j][j]=TINY_NUMBER;
    dum=1.0/(a[j][j]);
    for (i=j+1;i<n;i++) a[i][j] *= dum;
  }
  free(rowscale);
}

/*----------------------------------------------------------------------------*/
/*! \fn void lubksb_nr(Real **a, int n, int *indx, Real b[])
 *  Backward substitution (from numerical recipies)
 *  a is the input matrix done with LU decomposition, n is the matrix size
 *  indx id the history of row permutation
 *  b is the vector on the right (AX=b), and is returned with the solution
 */
void lubksb_nr(Real **a, int n, int *indx, Real b[])
{
  int i,ii=-1,ip,j;
  Real sum;
  /* Solve L*y=b */
  for (i=0;i<n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii>=0)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  /* Solve U*x=y */
  for (i=n-1;i>=0;i--) {
    sum=b[i];
    for (j=i+1;j<n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void InverseMatrix(Real **a, int n, Real **b)
 *  Inverse matrix solver
 *  a: input matrix; n: matrix size, b: return matrix
 *  Note: the input matrix will be DESTROYED
 */
void InverseMatrix(Real **a, int n, Real **b)
{
  int i,j,*indx;
  Real *col,d;

  indx = (int*)calloc_1d_array(n, sizeof(int));
  col = (Real*)calloc_1d_array(n, sizeof(Real));

  ludcmp_nr(a,n,indx,&d);

  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) col[i]=0.0;
    col[j]=1.0;
    lubksb_nr(a, n, indx, col);
    for (i=0; i<n; i++)    b[i][j] = col[i];
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void MatrixMult(Real **a, Real *b, int m, int n, Real *c)
 *  Matrix multiplication: a(m*n) * b(n) = c(m) */
void MatrixMult(Real **a, Real *b, int m, int n, Real *c)
{
  int i, j;
  for (i=0; i<m; i++) {
    c[i] = 0.0;
    for (j=0; j<n; j++) c[i] += a[i][j] * b[j];
  }
}

/*----------------------------------------------------------------------------*/
/*! \fn void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
 * gauss-legendre weight routine from numerical recipes */
void gauleg(Real x1, Real x2,  Real *x, Real *w, int n)
{

  Real eps = 3.0e-14;
  Real xm, xl, z, z1;
  Real p1, p2, p3, pp;
  int m, i, j;

  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);

  for (i=1; i<=m; i++) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(j=1; j<=n; j++) {
	p3 = p2;
	p2 = p1;
	p1 = ((2.0 * (Real)j - 1.0) * z * p2 - ((Real)j - 1.0) * p3) / (Real)j;
      }
      pp = (Real)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x[i-1] = xm - xl * z;
    x[n-i] = xm + xl * z;
    w[i-1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w[n-i] = w[i-1];
  }

}

#endif /* RADIATION_TRANSFER */
