#include "../copyright.h"
/*==============================================================================
 * FILE: BackEuler.c
 *
 * PURPOSE: Use backward Euler method to update the radiation quantities
 * First set up the matrix
 * Then solve the matrix equations.
 * We need the flag for boundary condition.
 *
 * Backward Euler should be used for the whole mesh
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   BackEuler_2d()
 * PRIVATE FUNCTION
 *      INI_matrix()
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

#ifndef FLD

#if defined(MATRIX_MULTIGRID) || defined(MATRIX_HYPRE)

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */



static MatrixS *pMat;

/* memory at each coarse level */
/* Allocate the memory for U, Ugas and RHS in advance */
static RadMHDS ****U_coarse;
static RadCoefS ****Ugas_coarse;
static Real *****RHS_coarse;


static Real ****INIerror;

static Real ****Reshist; /* The residual after each cycle, we only store 5 cycles for iterant recombination */
static Real ****Solhist; /* Store the history of the solution */
static Real *Resnorm; /* store the norm of the stored residual */
static int MAXerror = 1; /* The maximum number of combined residual */
static int Resnum; /* actual number of residual stored, Resnum <= MAXerror */

/* The pointer array to store the address of the coefficient */
/* This make sure that matrix coefficient only calculate once */
static Real ****Ptheta;
static Real ****Pphi;
static Real ****Ppsi;

static int coefflag; /* To decide whether need to calculate coefficient or not */

static Real INInorm;


static int Nlim = 4; /* the lim size of the coarsest grid in each CPU*/
static int Wcyclelim = 15; /* The limit of how many Wcycle we allow */

static int Matrixflag = 1; /* This is used to choose Gauss-Seidel or Jacobi method */
/* 1 is GS method while 0 is Jacobi */

static int Nlevel; /* Number of levels from top to bottom, log2(N)=Nlevel */
static int Wflag; /* To decide the position in the W flag */

/* subtract background state, which is initial condition */


/********Private function ************/
static Real CheckResidual(MatrixS *pMat, Real ***vector);
static void RadMHD_multig_2D(MatrixS *pMat);
static void prolongation2D(MatrixS *pMat_coarse, MatrixS *pMat_fine);
static void Restriction2D(MatrixS *pMat_fine, MatrixS *pMat_coarse);

static void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat);
static void RHSResidual2D(MatrixS *pMat, Real ***newRHS);


static void Calculate_Coef2D(MatrixS *pMat);

static void CopySolution2D(MatrixS *pMat, Real ***Solhist);

static void Recombination2D(MatrixS *pMat, Real ****Reshist, Real ****Solhist, const int num, const Real *Norms);

/* calculate inner product of two arbitrary vectors */
static void Inner_product2D(MatrixS *pMat, Real ***vector1, Real ***vector2, Real *result);



/*
  void Jacobi3D(MatrixS *pMat, MatrixS *pMatnew);
*/
extern void GaussSeidel2D(MatrixS *pMat,Real ***theta,  Real ***phi,  Real ***psi);
extern void Jacobi2D(MatrixS *pMat, Real ***theta,  Real ***phi,  Real ***psi);

/********Public function****************/
/*-------BackEuler_3d(): Use back euler method to update E_r and Fluxr-----------*/
/* Not work with SMR now. So there is only one Domain in the Mesh */


void BackEuler_2d(MeshS *pM)
{



  DomainS *pD;
  pD= &(pM->Domain[0][0]);

  GridS *pG=pD->Grid;

  Real dt = pG->dt;

/* Set the parameters that will change with time in pMat */
  pMat->dt = dt;
  pMat->time = pG->time;

  Real velocity_x, velocity_y, T4, pressure, density, temperature, Fr0x, Fr0y;
  Real AdvFx[2], AdvFy[2];
  Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
  Real Sigma[NOPACITY];

  Real error;
  int Wcycle;

  int i, j, m;
  int is, ie, js, je, ks;
  int Mati, Matj, Matk;
/* Set the boundary */
  is = pG->is-Matghost;
  ie = pG->ie+Matghost;
  js = pG->js-Matghost;
  je = pG->je+Matghost;
  ks = pG->ks;
  Matk = ks;


  int myID;

#ifdef MPI_PARALLEL
  myID = myID_Comm_world;
#else
  myID = 0;
#endif


#ifdef FARGO
  Real qom, x1, x2, x3;
  qom = qshear * Omega_0;
#endif


/* First, do the advection step and update bounary */


/* Now copy the data */
  for(j=js; j<=je; j++){
    for(i=is; i<= ie; i++){


      Mati = i - (nghost - Matghost);
      Matj = j - (nghost - Matghost);

      velocity_x = pG->U[ks][j][i].M1 / pG->U[ks][j][i].d;
      velocity_y = pG->U[ks][j][i].M2 / pG->U[ks][j][i].d;

#ifdef FARGO
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);
      velocity_y -= qom * x1;
#endif


      T4 = pG->Tguess[ks][j][i];


      Sigma_sF = pG->U[ks][j][i].Sigma[0];
      Sigma_aF = pG->U[ks][j][i].Sigma[1];
      Sigma_aP = pG->U[ks][j][i].Sigma[2];
      Sigma_aE = pG->U[ks][j][i].Sigma[3];



      pMat->U[Matk][Matj][Mati].Er  = pG->U[ks][j][i].Er;
      pMat->U[Matk][Matj][Mati].Fr1 = pG->U[ks][j][i].Fr1;
      pMat->U[Matk][Matj][Mati].Fr2 = pG->U[ks][j][i].Fr2;
/* Store the background Er */
      pMat->Ugas[Matk][Matj][Mati].rho = pG->U[ks][j][i].Er;
      pMat->Ugas[Matk][Matj][Mati].V1  = velocity_x;
      pMat->Ugas[Matk][Matj][Mati].V2  = velocity_y;
      pMat->Ugas[Matk][Matj][Mati].T4  = T4;
      pMat->Ugas[Matk][Matj][Mati].Edd_11 = pG->U[ks][j][i].Edd_11;
      pMat->Ugas[Matk][Matj][Mati].Edd_21 = pG->U[ks][j][i].Edd_21;
      pMat->Ugas[Matk][Matj][Mati].Edd_22 = pG->U[ks][j][i].Edd_22;

      pMat->Ugas[Matk][Matj][Mati].Sigma[0] = Sigma_sF;
      pMat->Ugas[Matk][Matj][Mati].Sigma[1] = Sigma_aF;
      pMat->Ugas[Matk][Matj][Mati].Sigma[2] = Sigma_aP;
      pMat->Ugas[Matk][Matj][Mati].Sigma[3] = Sigma_aE;

      Rad_Advection_Flux2D(pD, i, j, ks, 1.0, AdvFx, AdvFy);

      pMat->RHS[Matk][Matj][Mati][0] = pG->U[ks][j][i].Er + dt * Sigma_aP * T4 * ReduceC * Crat * Eratio +  (1.0 - Eratio) * pG->Ersource[ks][j][i] + ((AdvFx[1] - AdvFx[0]) + (AdvFy[1] - AdvFy[0]));
      pMat->RHS[Matk][Matj][Mati][1] = pG->U[ks][j][i].Fr1 + ReduceC * Eratio * dt * Sigma_aP * T4 * velocity_x + (1.0 - Eratio) * pG->Ersource[ks][j][i] * velocity_x / Crat;
      pMat->RHS[Matk][Matj][Mati][2] = pG->U[ks][j][i].Fr2 + ReduceC * Eratio * dt * Sigma_aP * T4 * velocity_y + (1.0 - Eratio) * pG->Ersource[ks][j][i] * velocity_y / Crat;



    } /* End i */
  }/* End j */


/* First, calculate the coefficient for the top level  */
/* coefficient will not change until next time step */
  Calculate_Coef2D(pMat);

  if(pMat->bgflag){
/* calculate the residual */
    RHSResidual2D(pMat, &(INIerror[0][0]));

/* Now set the matrix to be the residual */

    for(j=pMat->js-Matghost; j<=pMat->je+Matghost; j++){
      for(i=pMat->is-Matghost; i<= pMat->ie+Matghost; i++){
        pMat->U[ks][j][i].Er = 0.0;
        pMat->U[ks][j][i].Fr1 = 0.0;
        pMat->U[ks][j][i].Fr2 = 0.0;


        pMat->RHS[ks][j][i][0] = INIerror[ks][j][i][0];
        pMat->RHS[ks][j][i][1] = INIerror[ks][j][i][1];
        pMat->RHS[ks][j][i][2] = INIerror[ks][j][i][2];


      }
    }
  }
/* Right hand side in the ghost zones are never used */

/* calculate the initial norm of the right hand side */
  INInorm = CheckResidual(pMat, &(pMat->RHS[0][0]));
  if(INInorm < TINY_NUMBER)
    return;

/* Do the multi-grid W cycle */
  Wcycle = 0;
  error = 1.0;
  Resnum = 0;
  coefflag = 1;
/* Stop iteration if tolerance level is reached */
/* Or the matrix doesn't converge in Wcycle limit */
  while((error > TOL) && (Wcycle < Wcyclelim)){
    Wflag = 1;

    if(Resnum > 1){
      Recombination2D(pMat, Reshist, Solhist, Resnum, Resnorm);
/* After recombination, replace solhist with improved solution */
/*      Resnorm[Resnum-1] = CheckResidual(pMat,Reshist[Resnum-1]);

        CopySolution2D(pMat,Solhist[Resnum-1]);
*/
    }

    if(Resnum >= MAXerror)
      Resnum = 0;

    RadMHD_multig_2D(pMat);

/* for parent grid, check residual */
    RHSResidual2D(pMat, Reshist[Resnum]);
    error = CheckResidual(pMat,Reshist[Resnum]);
/* store the current solution and its norm */
    CopySolution2D(pMat,Solhist[Resnum]);
    Resnorm[Resnum] = error;
    Resnum++;

    if(error != error)
      ath_error("[BackEuler3D]: NaN encountered!\n");

    error /= INInorm;

    Wcycle++;

    coefflag = 0;


  }
/* Only output the residual for parent grid */
  if(myID == 0)
    ath_pout(0,"Final residual: %e  Cycle No.: %d\n",error,Wcycle);


/* copy the data back */
/* Now copy the data */
  for(j=js; j<=je; j++)
    for(i=is; i<= ie; i++){


      Mati = i - (nghost - Matghost);
      Matj = j - (nghost - Matghost);

      T4 = pG->U[ks][j][i].Er;

      if(pMat->bgflag){
        pG->U[ks][j][i].Er += pMat->U[Matk][Matj][Mati].Er;
        pG->U[ks][j][i].Fr1 += pMat->U[Matk][Matj][Mati].Fr1;
        pG->U[ks][j][i].Fr2 += pMat->U[Matk][Matj][Mati].Fr2;
      }
      else{
        pG->U[ks][j][i].Er = pMat->U[Matk][Matj][Mati].Er;
        pG->U[ks][j][i].Fr1 = pMat->U[Matk][Matj][Mati].Fr1;
        pG->U[ks][j][i].Fr2 = pMat->U[Matk][Matj][Mati].Fr2;
      }

      velocity_x = pMat->Ugas[Matk][Matj][Mati].V1;
      velocity_y = pMat->Ugas[Matk][Matj][Mati].V2;


      Fr0x =  pG->U[ks][j][i].Fr1 - ((1.0 +  pG->U[ks][j][i].Edd_11) * velocity_x +  pG->U[ks][j][i].Edd_21 * velocity_y) *  pG->U[ks][j][i].Er / Crat;
      Fr0y =  pG->U[ks][j][i].Fr2 - ((1.0 +  pG->U[ks][j][i].Edd_22) * velocity_y +  pG->U[ks][j][i].Edd_21 * velocity_x) *  pG->U[ks][j][i].Er / Crat;

/* Estimate the added energy source term */
      if(Prat > 0.0){

        pG->U[ks][j][i].Er += (pG->Eulersource[ks][j][i] - ReduceC * dt * (pG->U[ks][j][i].Sigma[1] -  pG->U[ks][j][i].Sigma[0]) * (velocity_x * Fr0x + velocity_y * Fr0y));


      }
    }


/* Set the boundary condition */
/*
  if(Opacity != NULL){

  for (j=pG->js; j<=pG->je; j++) {
  for (i=pG->is; i<=pG->ie; i++){

  density = pG->U[ks][j][i].d;

  pressure = (pG->U[ks][j][i].E - 0.5 * (pG->U[ks][j][i].M1 * pG->U[ks][j][i].M1
  + pG->U[ks][j][i].M2 * pG->U[ks][j][i].M2 +  pG->U[ks][j][i].M3 * pG->U[ks][j][i].M3) / density )     * (Gamma - 1);

  #ifdef RADIATION_MHD
  pressure -= 0.5 * (pG->U[ks][j][i].B1c * pG->U[ks][j][i].B1c + pG->U[ks][j][i].B2c * pG->U[ks][j][i].B2c + pG->U[ks][j][i].B3c * pG->U[ks][j][i].B3c) * (Gamma - 1.0);
  #endif
  temperature = pressure / (density * R_ideal);

  Opacity(density,temperature,Sigma,NULL);
  for(m=0;m<NOPACITY;m++){
  pG->U[ks][j][i].Sigma[m] = Sigma[m];
  }


  }
  }
  }
*/

/* Update the ghost zones for different boundary condition to be used later */
/*      for (i=0; i<pM->NLevels; i++){
        for (j=0; j<pM->DomainsPerLevel[i]; j++){
        if (pM->Domain[i][j].Grid != NULL){
        bvals_radMHD(&(pM->Domain[i][j]));

        }
        }
        }
*/
  return;


}


void RadMHD_multig_2D(MatrixS *pMat)
{

  MatrixS pMat_coarse;
  int Nsmall;

/* Once any dimension reaches size limit of Nlim, do Ncycle iteration and return */
  if(pMat->Nx[0] <= Nlim || pMat->Nx[1] <= Nlim){
/* Create the temporary array */
/* Need to create the temporary array for the boundary condition */


    bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
    bvals_Matrix_shear_init(pMat);
#endif


    if(Matrixflag){
      GaussSeidel2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
    }
    else{
      Jacobi2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
    }


    bvals_Matrix_destruct(pMat);

#ifdef SHEARING_BOX
    bvals_Matrix_shear_destruct();
#endif



  }
  else{

/* First, do relaxation at the fine level during the process of going down */
/*
  bvals_Matrix_init(pMat);

  GaussSeidel3D(pMat);

  bvals_Matrix_destruct(pMat);
*/

/* create the coarse grid */
    set_mat_level(&(pMat_coarse), pMat);


/* project the data to coarse grid */
    Restriction2D(pMat, &pMat_coarse);
/* Need to set boundary conditions after restriction */

    if(coefflag){

      bvals_Matrix_gas_init(&pMat_coarse);
#ifdef SHEARING_BOX
      bvals_Matrix_shear_gas_init(&pMat_coarse);
#endif

      bvals_Matrix_gas(&pMat_coarse);

#ifdef SHEARING_BOX
      bvals_Matrix_shear_gas_destruct();
#endif


      bvals_Matrix_gas_destruct(&pMat_coarse);

/* first time, calculate the coefficient */

      Calculate_Coef2D(&(pMat_coarse));
    }


/* continue the V or W cycle recursively*/
    RadMHD_multig_2D(&pMat_coarse);

/* The following code is first reached after Ncycle iterations at the coarsest
 * level.  We then prolongate, do Ncycle iterations, and return.  This will return
 * execution to this same spot for the next coarsest level, so we will
 * prolongate, do Ncycle iterations, return, and so on.
 */

/* Add the correction back to the fine grid */
    prolongation2D(&pMat_coarse, pMat);



/* Do relaxation when going up*/


    bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
    bvals_Matrix_shear_init(pMat);
#endif

/* Update the ghost zones first */
    bvals_Matrix(pMat);

    if(Matrixflag){
      GaussSeidel2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
    }
    else{
      Jacobi2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
    }

#ifdef SHEARING_BOX
    bvals_Matrix_shear_destruct();
#endif


    bvals_Matrix_destruct(pMat);

/* To decide whether go W cycle */
    Nsmall = pMat->Nx[0];
    if(pMat->Nx[1] < pMat->Nx[0]) Nsmall = pMat->Nx[1];


    if((Nsmall > pow(2.0,Nlevel/2.0)) && Wflag){
/* Set is Wcycle flag */
/* This is Wcycle. We should go down again */
      Wflag = 0;

      set_mat_level(&(pMat_coarse), pMat);


/* project the data to coarse grid */
      Restriction2D(pMat, &pMat_coarse);

/* Do not need to update ghost zones inside W cycle */
/* gas quantities are not changed, matrix coefficients are not changed */
/* radiation quantities are always zero */
/*
  bvals_Matrix_init(&pMat_coarse);

  #ifdef SHEARING_BOX
  bvals_Matrix_shear_init(&pMat_coarse);
  #endif

  bvals_Matrix(&pMat_coarse);

  #ifdef SHEARING_BOX
  bvals_Matrix_shear_destruct();
  #endif


  bvals_Matrix_destruct(&pMat_coarse);

*/

/* continue the V or W cycle recursively*/

      RadMHD_multig_2D(&pMat_coarse);


/* Add the correction back to the fine grid */
      prolongation2D(&pMat_coarse, pMat);


/* Do relaxation when going up*/


      bvals_Matrix_init(pMat);

#ifdef SHEARING_BOX
      bvals_Matrix_shear_init(pMat);
#endif
/* Update the ghost zones first */
      bvals_Matrix(pMat);

      if(Matrixflag){
        GaussSeidel2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
      }
      else{
        Jacobi2D(pMat,Ptheta[pMat->Level],Pphi[pMat->Level],Ppsi[pMat->Level]);
      }

#ifdef SHEARING_BOX
      bvals_Matrix_shear_destruct();
#endif


      bvals_Matrix_destruct(pMat);

    }





  }


  return;
}





/*==========================================*/
/* Check the residual after one cycle */
/* This is the stop criterian */
Real CheckResidual(MatrixS *pMat, Real ***vector)
{
/* This function calculate the norm of the right hand side */
  Real Norm;

  Inner_product2D(pMat, vector, vector, &Norm);

  return Norm;

}



/* prolongation operator */


void prolongation2D(MatrixS *pMat_coarse, MatrixS *pMat_fine)
{
/* Add the correction from the coarse grid back to the solution in fine grid */
  int i, j;
  int is, ie, js, je, ks;
  int num;

  Real *ptr_fine[4];
  Real *ptr_coarse[5];

  is = pMat_coarse->is;
  ie = pMat_coarse->ie;
  js = pMat_coarse->js;
  je = pMat_coarse->je;
  ks = pMat_coarse->ks;


  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++){
/* Take the address */
      ptr_fine[0] = &(pMat_fine->U[ks][2*j  ][2*i  ].Er);
      ptr_fine[1] = &(pMat_fine->U[ks][2*j  ][2*i-1].Er);
      ptr_fine[2] = &(pMat_fine->U[ks][2*j-1][2*i  ].Er);
      ptr_fine[3] = &(pMat_fine->U[ks][2*j-1][2*i-1].Er);


      ptr_coarse[0] = &(pMat_coarse->U[ks][j ][i ].Er);
      ptr_coarse[1] = &(pMat_coarse->U[ks][j+1][i+1].Er);
      ptr_coarse[2] = &(pMat_coarse->U[ks][j+1][i-1].Er);
      ptr_coarse[3] = &(pMat_coarse->U[ks][j-1][i+1].Er);
      ptr_coarse[4] = &(pMat_coarse->U[ks][j-1][i-1].Er);



/* Only needs to update Er, Frx, Fry, Frz. */
/* Vx, Vy, Vz and T4 do not need to be updated */
      for(num=0; num<4; num++){
/* We access a continuous array */
        ptr_fine[0][num] += 0.75*ptr_coarse[0][num];
        ptr_fine[0][num] += 0.25*ptr_coarse[1][num];


        ptr_fine[1][num] += 0.75*ptr_coarse[0][num];
        ptr_fine[1][num] += 0.25*ptr_coarse[2][num];

        ptr_fine[2][num] += 0.75*ptr_coarse[0][num];
        ptr_fine[2][num] += 0.25*ptr_coarse[3][num];

        ptr_fine[3][num] += 0.75*ptr_coarse[0][num];
        ptr_fine[3][num] += 0.25*ptr_coarse[4][num];

      }
    }

}


/* Restriction operator */

void Restriction2D(MatrixS *pMat_fine, MatrixS *pMat_coarse)
{
/* We actually should send the residual to the next level */

  int i, j;
  int is, ie, js, je, ks;

  is = pMat_fine->is;
  ie = pMat_fine->ie;
  js = pMat_fine->js;
  je = pMat_fine->je;
  ks = pMat_fine->ks;


/* The right hand size is no-longer the original source terms. It is the residual  */
/* But we cannot destroy the right hand size of the fine grid. We need that for the */
/* relaxation step when we come back */




  Real ****error;
  if((error=(Real****)calloc_4d_array(1,pMat_fine->Nx[1]+2*Matghost,pMat_fine->Nx[0]+2*Matghost,3,sizeof(Real))) == NULL)
    ath_error("[Restriction3D]: malloc return a NULL pointer\n");

  Real *ptr_coarse;
  Real *ptr_fine[4];
  int num;


/* calculate the residual */
  RHSResidual2D(pMat_fine, &(error[0][0]));




/* error = b - Ax */
/* Now restrict to the coarse grid */

  for(j=pMat_coarse->js; j<=pMat_coarse->je; j++)
    for(i=pMat_coarse->is; i<=pMat_coarse->ie; i++){
      if(coefflag){
        ptr_coarse  = &(pMat_coarse->Ugas[ks][j][i].rho);
        ptr_fine[0] = &(pMat_fine->Ugas[ks][2*j ][2*i ].rho);
        ptr_fine[1] = &(pMat_fine->Ugas[ks][2*j ][2*i-1].rho);
        ptr_fine[2] = &(pMat_fine->Ugas[ks][2*j-1][2*i ].rho);
        ptr_fine[3] = &(pMat_fine->Ugas[ks][2*j-1][2*i-1].rho);

        for(num=0; num<11+NOPACITY; num++){

          ptr_coarse[num] =  (ptr_fine[0][num] + ptr_fine[1][num]
                              + ptr_fine[2][num] + ptr_fine[3][num] ) / 4.0;

        }
      }
/*
  pMat_coarse->U[ks][j][i] = (pMat_fine->U[2*k ][2*j ][2*i ]  + pMat_fine->U[2*k ][2*j ][2*i-1]
  + pMat_fine->U[2*k ][2*j-1][2*i ] + pMat_fine->U[2*k ][2*j-1][2*i-1]
  + pMat_fine->U[2*k-1][2*j ][2*i ] + pMat_fine->U[2*k-1][2*j ][2*i-1]
  + pMat_fine->U[2*k-1][2*j-1][2*i ]+ pMat_fine->U[2*k-1][2*j-1][2*i-1]) / 8.0;

*/
      pMat_coarse->RHS[ks][j][i][0] =  (error[ks][2*j ][2*i ][0] + error[ks][2*j ][2*i-1][0]
                                        +       error[ks][2*j-1][2*i ][0]+ error[ks][2*j-1][2*i-1][0]) / 4.0;

      pMat_coarse->RHS[ks][j][i][1] =  (error[ks][2*j ][2*i ][1] + error[ks][2*j ][2*i-1][1]
                                        +       error[ks][2*j-1][2*i ][1]+ error[ks][2*j-1][2*i-1][1]) / 4.0;

      pMat_coarse->RHS[ks][j][i][2] =  (error[ks][2*j ][2*i ][2] + error[ks][2*j ][2*i-1][2]
                                        +       error[ks][2*j-1][2*i ][2]+ error[ks][2*j-1][2*i-1][2]) / 4.0;

    }

  for(j=pMat_coarse->js-Matghost; j<=pMat_coarse->je+Matghost; j++)
    for(i=pMat_coarse->is-Matghost; i<=pMat_coarse->ie+Matghost; i++){
/* The initial guess is taken to be zero */
      pMat_coarse->U[ks][j][i].Er = 0.0;
      pMat_coarse->U[ks][j][i].Fr1 = 0.0;
      pMat_coarse->U[ks][j][i].Fr2 = 0.0;

    }


/* Free the temporary array */
  free_4d_array(error);

  return;
}



/*  calculate the residual of right hand side from guess solution */
/* should allocate memory for newRHS before this function is called */
void RHSResidual2D(MatrixS *pMat, Real ***newRHS)
{

  int i, j, m, n;
  int is, ie, js, je, ks;

  is = pMat->is;
  ie = pMat->ie;
  js = pMat->js;
  je = pMat->je;
  ks = pMat->ks;



/* To store the coefficient */
  Real theta[11];
  Real phi[11];
  Real psi[11];

/* current level */
  n = pMat->Level;

  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++){

/* get the coefficient */
      for(m=0; m<11; m++){
        theta[m]  = Ptheta[n][j][i][m];
        phi[m]    = Pphi[n][j][i][m];
        psi[m]    = Ppsi[n][j][i][m];

      }


      newRHS[j][i][0] = pMat->RHS[ks][j][i][0];
      newRHS[j][i][0] -= theta[0] * pMat->U[ks][j-1][i].Er;
      newRHS[j][i][0] -= theta[1] * pMat->U[ks][j-1][i].Fr2;
      newRHS[j][i][0] -= theta[2] * pMat->U[ks][j][i-1].Er;
      newRHS[j][i][0] -= theta[3] * pMat->U[ks][j][i-1].Fr1;
      newRHS[j][i][0] -= theta[4] * pMat->U[ks][j][i].Er;
      newRHS[j][i][0] -= theta[5] * pMat->U[ks][j][i].Fr1;
      newRHS[j][i][0] -= theta[6] * pMat->U[ks][j][i].Fr2;
      newRHS[j][i][0] -= theta[7] * pMat->U[ks][j][i+1].Er;
      newRHS[j][i][0] -= theta[8] * pMat->U[ks][j][i+1].Fr1;
      newRHS[j][i][0] -= theta[9] * pMat->U[ks][j+1][i].Er;
      newRHS[j][i][0] -= theta[10] * pMat->U[ks][j+1][i].Fr2;



      newRHS[j][i][1] = pMat->RHS[ks][j][i][1];
      newRHS[j][i][1] -= phi[0] * pMat->U[ks][j-1][i].Er;
      newRHS[j][i][1] -= phi[1] * pMat->U[ks][j-1][i].Fr1;
      newRHS[j][i][1] -= phi[2] * pMat->U[ks][j][i-1].Er;
      newRHS[j][i][1] -= phi[3] * pMat->U[ks][j][i-1].Fr1;
      newRHS[j][i][1] -= phi[4] * pMat->U[ks][j][i].Er;
      newRHS[j][i][1] -= phi[5] * pMat->U[ks][j][i].Fr1;
      newRHS[j][i][1] -= phi[6] * pMat->U[ks][j][i+1].Er;
      newRHS[j][i][1] -= phi[7] * pMat->U[ks][j][i+1].Fr1;
      newRHS[j][i][1] -= phi[8] * pMat->U[ks][j+1][i].Er;
      newRHS[j][i][1] -= phi[9] * pMat->U[ks][j+1][i].Fr1;


      newRHS[j][i][2] = pMat->RHS[ks][j][i][2];
      newRHS[j][i][2] -= psi[0] * pMat->U[ks][j-1][i].Er;
      newRHS[j][i][2] -= psi[1] * pMat->U[ks][j-1][i].Fr2;
      newRHS[j][i][2] -= psi[2] * pMat->U[ks][j][i-1].Er;
      newRHS[j][i][2] -= psi[3] * pMat->U[ks][j][i-1].Fr2;
      newRHS[j][i][2] -= psi[4] * pMat->U[ks][j][i].Er;
      newRHS[j][i][2] -= psi[5] * pMat->U[ks][j][i].Fr2;
      newRHS[j][i][2] -= psi[6] * pMat->U[ks][j][i+1].Er;
      newRHS[j][i][2] -= psi[7] * pMat->U[ks][j][i+1].Fr2;
      newRHS[j][i][2] -= psi[8] * pMat->U[ks][j+1][i].Er;
      newRHS[j][i][2] -= psi[9] * pMat->U[ks][j+1][i].Fr2;



    }

/* error = b - Ax */

  return;

}



/* Calculate the matrix coefficient for each level */
void Calculate_Coef2D(MatrixS *pMat)
{
  int i, j;
  int is, ie, js, je, ks;
  int n;

  is = pMat->is;
  ie = pMat->ie;
  js = pMat->js;
  je = pMat->je;
  ks = pMat->ks;

/* current level */
  n = pMat->Level;

  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++){

      matrix_coef(pMat, NULL, 2, i, j, ks, 0.0, &(Ptheta[n][j][i][0]), &(Pphi[n][j][i][0]), &(Ppsi[n][j][i][0]), NULL);
    }

}
void Recombination2D(MatrixS *pMat, Real ****Reshist, Real ****Solhist, const int num, const Real *Norms)
{
/* combine the latest 0 ... num-1 solution to minimize the current residual */
/* Norms[i] = <Reshist[i], Reshist[i]> */

  int i, j, m;
  int is, ie, js, je, ks;
  int flag = 1; /* used to label whether matrix inversion sucesseed or not */

  is = pMat->is;
  ie = pMat->ie;
  js = pMat->js;
  je = pMat->je;
  ks = pMat->ks;

  Real **Matrixcoef = NULL;
  Real *MatrixRHS = NULL;
  Real *crossNorm = NULL;

/* temporary used for LU decomposition */
  int *index = NULL;
  Real *dtemp = NULL;


  Real tempNorm;

  if((Matrixcoef=(Real**)calloc_2d_array(num,num,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((MatrixRHS=(Real*)calloc(num,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((crossNorm=(Real*)calloc(num,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((index=(int*)calloc(num,sizeof(int))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((dtemp=(Real*)calloc(num,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");



/* set the coefficient of H matrix and right hand side */
/* crossNorm is <Res[num-1], Res[i]> */

  for(i=1; i<=num-1; i++){
    Inner_product2D(pMat, Reshist[num-1-i], Reshist[num-1], &(crossNorm[i]));
    MatrixRHS[i] = Norms[num-1] - crossNorm[i];
  }


  for(i=1; i<= num-1; i++){
    for(j=1; j<= num-1; j++){
      if(i != j){
        Inner_product2D(pMat, Reshist[num-1-i], Reshist[num-1-j], &(tempNorm));
        Matrixcoef[i][j] = tempNorm - crossNorm[i] - crossNorm[j] + Norms[num-1];
      }
      else{
        Matrixcoef[i][j] = Norms[num-1-i] - 2.0 * crossNorm[i] +  Norms[num-1];

      }
    }
  }

/* now solve the (num-1)(num-1) matrix:  Matrixcoef alphas = MatrixRHS */
/* This is smaller than 5 * 5 matrix, can use direct method to solve */
/* LU decomposition */
  ludcmpnew(Matrixcoef,num-1, index, dtemp, &flag);
/* If matrix is ill conditinoed, do not do recombination */
/* solve the matrix */
  if(flag){

    lubksb(Matrixcoef, num-1, index, MatrixRHS);
/* solution is stored in MatrixRHS */

    for(j=js; j<=je; j++)
      for(i=is; i<=ie; i++){
        for(m=1; m<num; m++){
          pMat->U[ks][j][i].Er  += MatrixRHS[m] * (Solhist[num-1-m][j][i][0] - Solhist[num-1][j][i][0]);
          pMat->U[ks][j][i].Fr1 += MatrixRHS[m] * (Solhist[num-1-m][j][i][1] - Solhist[num-1][j][i][1]);
          pMat->U[ks][j][i].Fr2 += MatrixRHS[m] * (Solhist[num-1-m][j][i][2] - Solhist[num-1][j][i][2]);
        }
      }

  }
/* The more accurate solution is now stored in pMat */




  if(Matrixcoef != NULL)
    free_2d_array(Matrixcoef);

  if(MatrixRHS != NULL)
    free(MatrixRHS);

  if(crossNorm != NULL)
    free(crossNorm);


  if(index != NULL)
    free(index);


  if(dtemp != NULL)
    free(dtemp);



}


/* Copy and store current estimate solution */
void CopySolution2D(MatrixS *pMat, Real ***Solhist)
{
  int i, j;
  int is, ie, js, je, ks;

  is = pMat->is;
  ie = pMat->ie;
  js = pMat->js;
  je = pMat->je;
  ks = pMat->ks;




  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++){

      Solhist[j][i][0] = pMat->U[ks][j][i].Er;
      Solhist[j][i][1] = pMat->U[ks][j][i].Fr1;
      Solhist[j][i][2] = pMat->U[ks][j][i].Fr2;
    }

  return;
}

/* calculate the inner product of two vectors , handle MPI case */
void Inner_product2D(MatrixS *pMat, Real ***vector1, Real ***vector2, Real *result)
{


  Real Norm;
  Real normtemp;

  int i, j;
  int is, ie, js, je, ks;
  is = pMat->is;
  ie = pMat->ie;
  js = pMat->js;
  je = pMat->je;
  ks = pMat->ks;

#ifdef MPI_PARALLEL
  int ierr;
  double tot_norm = 0.0;
#endif

  Norm = 0.0;


  for(j=js; j<=je; j++)
    for(i=is; i<=ie; i++){
      vector_product(vector1[j][i],vector2[j][i], 3 , &normtemp);
      Norm += normtemp;
    }

/* MPI call to collect norm from different CPUs */


#ifdef MPI_PARALLEL
  ierr = MPI_Allreduce(&Norm,&tot_norm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  Norm = tot_norm;
#endif

  *result = Norm;

  return;

}






void set_mat_level(MatrixS *pMat_coarse, MatrixS *pMat)
{

#ifdef MPI_PARALLEL
  pMat_coarse->Comm_Domain = pMat->Comm_Domain;
#endif

  pMat_coarse->dx1 = 2.0 * pMat->dx1;
  pMat_coarse->dx2 = 2.0 * pMat->dx2;
  pMat_coarse->dx3 = pMat->dx3;
  pMat_coarse->time = pMat->time;
  pMat_coarse->Lx = pMat->Lx;
  pMat_coarse->Ly = pMat->Ly;
  pMat_coarse->Lz = pMat->Lz;
/* Grid numbers decrease by a factor of 2 */
  pMat_coarse->Nx[0] = pMat->Nx[0] / 2;
  pMat_coarse->Nx[1] = pMat->Nx[1] / 2;
  pMat_coarse->Nx[2] = pMat->Nx[2];
  pMat_coarse->Level = pMat->Level + 1;
/* top level grid number doesn't change */
  pMat_coarse->RootNx[0] = pMat->RootNx[0];
  pMat_coarse->RootNx[1] = pMat->RootNx[1];
  pMat_coarse->RootNx[2] = pMat->RootNx[2];
  pMat_coarse->NGrid[0] = pMat->NGrid[0];
  pMat_coarse->NGrid[1] = pMat->NGrid[1];
  pMat_coarse->NGrid[2] = pMat->NGrid[2];
  pMat_coarse->dt = pMat->dt;
  pMat_coarse->is = Matghost;
  pMat_coarse->ie = pMat_coarse->Nx[0] + Matghost - 1;
  pMat_coarse->js = Matghost;
  pMat_coarse->je = pMat_coarse->Nx[1] + Matghost - 1;
  pMat_coarse->ks = 0;
  pMat_coarse->ke = 0;

  pMat_coarse->ID = pMat->ID;
  pMat_coarse->my_iproc = pMat->my_iproc;
  pMat_coarse->my_jproc = pMat->my_jproc;
  pMat_coarse->my_kproc = pMat->my_kproc;

  pMat_coarse->rx1_id = pMat->rx1_id;
  pMat_coarse->lx1_id = pMat->lx1_id;
  pMat_coarse->rx2_id = pMat->rx2_id;
  pMat_coarse->lx2_id = pMat->lx2_id;
  pMat_coarse->lx3_id = pMat->lx3_id;
  pMat_coarse->rx3_id = pMat->rx3_id;

/* Boundary flag */
  pMat_coarse->BCFlag_ix1 = pMat->BCFlag_ix1;
  pMat_coarse->BCFlag_ox1 = pMat->BCFlag_ox1;
  pMat_coarse->BCFlag_ix2 = pMat->BCFlag_ix2;
  pMat_coarse->BCFlag_ox2 = pMat->BCFlag_ox2;
  pMat_coarse->BCFlag_ix3 = pMat->BCFlag_ix3;
  pMat_coarse->BCFlag_ox3 = pMat->BCFlag_ox3;

/* memory is already relocated. Just set the pointer */

  pMat_coarse->U    = U_coarse[pMat_coarse->Level];
  pMat_coarse->Ugas = Ugas_coarse[pMat_coarse->Level];
  pMat_coarse->RHS  = RHS_coarse[pMat_coarse->Level];


  return;
}




/*-------------------------------------------------------------------------*/
/* BackEuler_init_2d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_2d(): function to free memory */
void BackEuler_init_2d(MeshS *pM)
{

  DomainS *pD;
  pD= &(pM->Domain[0][0]);

  GridS *pG=pD->Grid;
  int Nx, Ny, Nz;
  int Nx2, Ny2;
  int i;
  Real temp;


  Nx = pG->ie - pG->is + 1;
  Ny = pG->je - pG->js + 1;
  Nz = 1;

/* Reach bottom first for the side with the smallest size*/
  Nlevel = Nx;
  if(Ny < Nx) Nlevel = Ny;


  Nlevel /= Nlim;

  temp = log10(Nlevel)/log10(2.0);

  Nlevel = (int)temp;
  if(fabs(temp-Nlevel) > 0.5) Nlevel++;
  Nlevel++;

/* pMat will remain in the memory until the end of the simulation */

  if((pMat = (MatrixS*)calloc(1,sizeof(MatrixS))) == NULL)
    ath_error("[BackEuler_init_2d]: malloc return a NULL pointer\n");

  if((pMat->U = (RadMHDS***)calloc_3d_array(Nz,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadMHDS))) == NULL)
    ath_error("[BackEuler_init_2d]: malloc return a NULL pointer\n");

  if((pMat->Ugas = (RadCoefS***)calloc_3d_array(Nz,Ny+2*Matghost, Nx+2*Matghost,sizeof(RadCoefS))) == NULL)
    ath_error("[BackEuler_init_2d]: malloc return a NULL pointer\n");


  if((pMat->RHS = (Real****)calloc_4d_array(Nz,Ny+2*Matghost, Nx+2*Matghost,3,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2d]: malloc return a NULL pointer\n");

  if((INIerror=(Real****)calloc_4d_array(Nz,Ny+2*Matghost,Nx+2*Matghost,3,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");


  if((Reshist=(Real****)calloc_4d_array(MAXerror,Ny+2*Matghost,Nx+2*Matghost,3,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((Solhist=(Real****)calloc_4d_array(MAXerror,Ny+2*Matghost,Nx+2*Matghost,3,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((Resnorm=(Real*)calloc(MAXerror,sizeof(Real))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

/*==================================================================*/
/* The pointer to store the pointer of the coefficients at each level */
  if((Ptheta=(Real****)calloc(Nlevel,sizeof(Real***))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");


  if((Pphi=(Real****)calloc(Nlevel,sizeof(Real***))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((Ppsi=(Real****)calloc(Nlevel,sizeof(Real***))) == NULL)
    ath_error("[BackEuler_init_2D]: malloc return a NULL pointer\n");

  if((U_coarse=(RadMHDS****)calloc(Nlevel,sizeof(RadMHDS***))) == NULL)
    ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

  if((Ugas_coarse=(RadCoefS****)calloc(Nlevel,sizeof(RadCoefS***))) == NULL)
    ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

  if((RHS_coarse=(Real*****)calloc(Nlevel,sizeof(Real****))) == NULL)
    ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

#ifdef MPI_PARALLEL
  pMat->Comm_Domain = pD->Comm_Domain;
#endif
/* allocate memory at each level */
  Nx2 = Nx;
  Ny2 = Ny;

  for(i=0; i<Nlevel; i++){

    if((Ptheta[i]=(Real***)calloc_3d_array(Ny2+2*Matghost,Nx2+2*Matghost,11,sizeof(Real))) == NULL)
      ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");
    if((Pphi[i]=(Real***)calloc_3d_array(Ny2+2*Matghost,Nx2+2*Matghost,11,sizeof(Real))) == NULL)
      ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");
    if((Ppsi[i]=(Real***)calloc_3d_array(Ny2+2*Matghost,Nx2+2*Matghost,11,sizeof(Real))) == NULL)
      ath_error("[BackEuler_init_3D]: malloc return a NULL pointer\n");

    if((U_coarse[i] = (RadMHDS***)calloc_3d_array(Nz,Ny2+2*Matghost, Nx2+2*Matghost,sizeof(RadMHDS))) == NULL)
      ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

    if((Ugas_coarse[i] = (RadCoefS***)calloc_3d_array(Nz,Ny2+2*Matghost, Nx2+2*Matghost,sizeof(RadCoefS))) == NULL)
      ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");

    if((RHS_coarse[i] = (Real****)calloc_4d_array(Nz,Ny2+2*Matghost, Nx2+2*Matghost,4,sizeof(Real))) == NULL)
      ath_error("[BackEuler_init_3d]: malloc return a NULL pointer\n");


    Nx2 /= 2;
    Ny2 /= 2;
  }
/*==========================================================================*/



/* now set the parameters */

  pMat->dx1 = pG->dx1;
  pMat->dx2 = pG->dx2;
  pMat->dx3 = pG->dx3;
  pMat->time = pG->time;
/* dt in pG is not set at this time */
/* pMat->dt = pG->dt;
 */
  pMat->Lx = pD->RootMaxX[0] - pD->RootMinX[0];
  pMat->Ly = pD->RootMaxX[1] - pD->RootMinX[1];
  pMat->Lz = pD->RootMaxX[2] - pD->RootMinX[2];

  pMat->is = Matghost;
  pMat->ie = Nx + Matghost - 1;
  pMat->js = Matghost;
  pMat->je = Ny + Matghost - 1;
  pMat->ks = 0;
  pMat->ke = 0;
  pMat->Nx[0] = Nx;
  pMat->Nx[1] = Ny;
  pMat->Nx[2] = Nz;
  pMat->Level = 0;
  pMat->RootNx[0] = Nx;
  pMat->RootNx[1] = Ny;
  pMat->RootNx[2] = Nz;
  pMat->NGrid[0] = pD->NGrid[0];
  pMat->NGrid[1] = pD->NGrid[1];
  pMat->NGrid[2] = pD->NGrid[2];
  pMat->rx1_id = pG->rx1_id;
  pMat->lx1_id = pG->lx1_id;
  pMat->rx2_id = pG->rx2_id;
  pMat->lx2_id = pG->lx2_id;
  pMat->lx3_id = pG->lx3_id;
  pMat->rx3_id = pG->rx3_id;
#ifdef MPI_PARALLEL
  pMat->ID = myID_Comm_world;
  get_myGridIndex(pD, myID_Comm_world, &(pMat->my_iproc), &(pMat->my_jproc), &(pMat->my_kproc));
#else
  pMat->ID = 0;
  pMat->my_iproc = 0;
  pMat->my_jproc = 0;
  pMat->my_kproc = 0;
#endif


/* Boundary flag */
  pMat->BCFlag_ix1 = pM->BCFlag_ix1;
  pMat->BCFlag_ox1 = pM->BCFlag_ox1;
  pMat->BCFlag_ix2 = pM->BCFlag_ix2;
  pMat->BCFlag_ox2 = pM->BCFlag_ox2;
  pMat->BCFlag_ix3 = pM->BCFlag_ix3;
  pMat->BCFlag_ox3 = pM->BCFlag_ox3;


/* To decide whether subtract background solution at top level or not */
/* Default choice is not */
  pMat->bgflag = 1;


}


void BackEuler_destruct_2d(MeshS *pM)
{
  int i;

/* Free pMat and pMatnew */
  if(pMat->U != NULL)
    free_3d_array(pMat->U);

  if(pMat->Ugas != NULL)
    free_3d_array(pMat->Ugas);

  if(pMat->RHS != NULL)
    free_4d_array(pMat->RHS);

  if(pMat != NULL)
    free(pMat);

  if(INIerror != NULL)
    free_4d_array(INIerror);

  if(Reshist != NULL)
    free_4d_array(Reshist);

  if(Solhist != NULL)
    free_4d_array(Solhist);

  if(Resnorm != NULL)
    free(Resnorm);


  for(i=0; i<Nlevel; i++){
    free_3d_array(Ptheta[i]);
    free_3d_array(Pphi[i]);
    free_3d_array(Ppsi[i]);

    free_3d_array(U_coarse[i]);
    free_3d_array(Ugas_coarse[i]);
    free_4d_array(RHS_coarse[i]);
  }

  free(Ptheta);
  free(Pphi);
  free(Ppsi);

  free(U_coarse);
  free(Ugas_coarse);
  free(RHS_coarse);

}



#endif /* radMHD_INTEGRATOR */

#endif /* MATRIX_MULTIGRID */


#endif /* FLD */

