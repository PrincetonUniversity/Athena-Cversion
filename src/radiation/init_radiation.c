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
  int nmu, qmeth;

  if (outflag == 0) { /* Initialize RadGrid array */
      pRG = pD->RadGrid;    /* set ptr to RadGrid */
      pRG->outflag = 0;
  } else {   /* Initialize RadOutGrid array */
      pRG = pD->RadOutGrid;          /* set ptr to RadOutGrid */
      pRG->outflag = 1;
  }
  pRG->pG = pD->Grid;  /* set correspond Grids pointer */

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
  if (outflag == 0) {
/* -------------------------- Initialize RadGrid array  -------------------------- */
    pRG->nf  = par_geti("radiation","nf");
    qmeth = par_geti_def("radiation","ang_quad",1);
    switch(qmeth) {
      
    case 0: /* User defined angular quadrature */
      pRG->nang = par_geti_def("radiation","nang",0);
      if (pRG->nang == 0) 
	ath_error("[init_radiation]: ang_quad = %d with nang = %d.\n",qmeth,pRG->nang);	
      break;
      
    case 1: /* Carlson symmetric S_N method (default) */
      nmu = par_geti_def("radiation","nmu",0);
      if (nmu <= 0) 
	ath_error("[init_radiation]: ang_quad = %d with nmu = %d.\n",qmeth,nmu);
    if (nDim == 1) {
      pRG->nang = nmu;
    } else {
      pRG->nang = nmu * (nmu + 1) / 2;
    }
    break;
      
    case 2: /* Legendre Equal Weight */
      nmu = par_geti_def("radiation","nmu",0);
      if (nmu == 0) 
	ath_error("[init_radiation]: ang_quad = %d with nmu = %d.\n",qmeth,nmu);
      if (nDim == 1) {
	pRG->nang = nmu;
      } else {
	pRG->nang = nmu * (nmu + 1) / 2;
      }
      break;

    case 10: /* single mu_z */
      nmu = par_geti_def("radiation","nmu",0);
      if (nmu == 0) 
	ath_error("[init_radiation]: ang_quad = %d with nmu = %d.\n",qmeth,nmu);
	pRG->nang = nmu;
      break;

    default:
      ath_error("[init_radiation]: ang_quad = %d.\n",qmeth);
    }
#ifdef RAY_TRACING
    pRG->nf_rt = par_geti_def("radiation","nf_rt",pRG->nf);
#endif
  } else {
/* -------------------------- Initialize RadOutGrid array  -------------------------- */  
    pRG->nf  = par_geti("radiation_output","nf");
    qmeth = par_geti_def("radiation_output","ang_quad",1);
    switch(qmeth) {
      
    case 0: /* User defined angular quadrature */
      pRG->nang = par_geti_def("radiation_output","nang",0);
      if (pRG->nang == 0) 
	ath_error("[init_radiation]: radiation_output: ang_quad = %d with nang = %d.\n",qmeth,pRG->nang);
      break;
      
    case 1: /* Carlson symmetric S_N method (default) */
      nmu = par_geti_def("radiation_output","nmu",0);
      if (nmu == 0) 
	ath_error("[init_radiation]: radiation_output: ang_quad = %d with nmu = %d.\n",qmeth,nmu);
      if (nDim == 1) {
	pRG->nang = nmu;
      } else {
	pRG->nang = nmu * (nmu + 1) / 2;
      }
      break;
      
    case 2: /* Legendre Equal Weight */
      nmu = par_geti_def("radiation_output","nmu",0);
      if (nmu == 0) 
	ath_error("[init_radiation]: radiation_output: ang_quad = %d with nmu = %d.\n",qmeth,nmu);
      if (nDim == 1) {
	pRG->nang = nmu;
      } else {
	pRG->nang = nmu * (nmu + 1) / 2;
      }
      break;
      
    default:
      ath_error("[init_radiation]: radiation_output: ang_quad = %d.\n",qmeth);
    }
#ifdef RAY_TRACING
    pRG->nf_rt = par_geti_def("radiation_output","nf_rt",pRG->nf);
#endif
  }

  if (nDim == 1) {
    pRG->noct = 2;
  } else  if (nDim == 2) {
    pRG->noct = 4;
  } else if (nDim == 3) {
    pRG->noct = 8;
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
  
  if (qmeth != 0) {
/* Call init_angles() to initialize default angle grid and angular quadrature
 * If qmeth = 0, angles/quadratures must be initialized in problem generator */
    init_angles(pRG,qmeth,outflag);
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


#endif /* RADIATION_TRANSFER */
