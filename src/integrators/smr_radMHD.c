#include "../copyright.h"
/*============================================================================*/
/*! \file smr_radMHD.c
 *  \brief The restriction and prolongation function for the matrix solver with SMR
 *
 * PURPOSE: 
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef STATIC_MESH_REFINEMENT
#ifdef RADIATIONMHD_INTEGRATOR
#ifdef MATRIX_MULTIGRID 

/* memory to send data */

static double **send_bufP= NULL; 
static double **send_bufRC=NULL; 
static double **recv_bufP= NULL;
#ifdef MPI_PARALLEL
static double **recv_bufRC=NULL;
static MPI_Request **recv_rq=NULL;
static MPI_Request  **send_rq=NULL;
#endif
static int maxND, *start_addrP;

static RadMHDS ***GZ[3];
static RadMHDS ***Pro_buf; /* temporary array to load data during prolongation */

static int RootLevel;

/*===================================*/
/* Private Function */
/*=====================================*/


static void ProU(const RadMHDS Uim1, const RadMHDS Ui, const RadMHDS Uip1, 
	  const RadMHDS Ujm1, const RadMHDS Ujp1,
	  const RadMHDS Ukm1, const RadMHDS Ukp1, RadMHDS PCon[2][2][2]);

#ifndef FIRST_ORDER
/* left: vl, right: vr, center: vc. Calculate the TVD slope */
static Real mcd_slope(const Real vl, const Real vc, const Real vr);
#endif /* FIRST_ORDER */



/*==============================================*/
/********Public functions **********************/
/*===============================================*/



/* This function only restricts the gas quantities from Level to Level-1 */
/* flag == 1: restrict the coefficient only; flag == 0:  restrict solution only */
/* We need pM here because we need to know the information about child and paraent grids */
/* Because the matrix coefficient needs primative variables, restrictions are done */
/* Even for levels above the root domain */


/* This function restricts the matrix coefficient, solution and RHS */
/* Assuming that RHS is already updated with the new solution */
/* We actually also copy the solution back from matrix to the grid */
/* Because this is just after the relaxation step */




void Rad_Restriction(MatrixS **Matrix, const int Level, const MeshS *pM, const int flag)
{
	/* When there is no data for restriction, we shouldn't initialize a MPI call for this case */
	/* That is what nZeroRC for */
	int i, j, k, ips, ipe, jps, jpe, kps, kpe, ics, ice, jcs, jce, kcs, kce, nDim;
	int nd, npg, ncg, start_addr, nZeroRC;
	MatrixS *pMat, *pMatP;	
	GridS *pG;

	double *pRcv, *pSnd;
	GridOvrlpS *pPO, *pCO;

	Real *ptrP;
	Real *ptr[8];
	int num;	
#ifdef MPI_PARALLEL
  	int ierr, mIndex, mCount, mAddress;
	int **RHStemp; /* temporary array used to store index and address for each child grid */
#endif

	/* judge the dimensionality of the problem */
	nDim = 1;	
	for(i=1; i<3; i++)
		if(pM->Nx[i] > 1) nDim++;
	

/*======== Step 3.: Restrict child solution and send =================*/
/* Loop over all Domains and parent GridS */
/* [Level] is the child level and [Level-1] is the parent level */
/* Separate levels above and below the RootLevel: Below the RootLevel, there is only one domain. */
/* Every child and parent grids are in the same CPU */
/* For levels below RootLevles, we do not need to restrict the solution */
/*==================================================================================*/
	/* Only do this if this CPU works here */
 	if(Level <= RootLevel){
	if(Matrix[0][0].CPUflag){
		pMat  = &(Matrix[Level][0]);
		pMatP = &(Matrix[Level-1][0]);

		/*****************************************/
		/* Do not calculate the Residual here */	
		/* We wil do this before restriction function is called, if necessary */
		

		for(k=pMatP->ks; k<=pMatP->ke; k++)
			for(j=pMatP->js; j<=pMatP->je; j++)
				for(i=pMatP->is; i<=pMatP->ie; i++){
					if(flag){ /* We will never need to restrict the solution for level below RootLevel */
						ptrP  = &(pMatP->Ugas[k][j][i].rho);
						ptr[0] = &(pMat->Ugas[2*k ][2*j ][2*i ].rho);
						ptr[1] = &(pMat->Ugas[2*k ][2*j ][2*i-1].rho);
						ptr[2] = &(pMat->Ugas[2*k ][2*j-1][2*i ].rho);
						ptr[3] = &(pMat->Ugas[2*k ][2*j-1][2*i-1].rho);
						ptr[4] = &(pMat->Ugas[2*k-1][2*j ][2*i ].rho);
						ptr[5] = &(pMat->Ugas[2*k-1][2*j ][2*i-1].rho);
						ptr[6] = &(pMat->Ugas[2*k-1][2*j-1][2*i ].rho);	
						ptr[7] = &(pMat->Ugas[2*k-1][2*j-1][2*i-1].rho);

						for(num=0; num<11+NOPACITY; num++){

							ptrP[num] =  (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
						}
					}

					/* Now for the RHS */
					ptr[0] = &(pMat->RHS[2*k ][2*j ][2*i ][0]);
					ptr[1] = &(pMat->RHS[2*k ][2*j ][2*i-1][0]);
					ptr[2] = &(pMat->RHS[2*k ][2*j-1][2*i ][0]);
					ptr[3] = &(pMat->RHS[2*k ][2*j-1][2*i-1][0]);
					ptr[4] = &(pMat->RHS[2*k-1][2*j ][2*i ][0]);
					ptr[5] = &(pMat->RHS[2*k-1][2*j ][2*i-1][0]);
					ptr[6] = &(pMat->RHS[2*k-1][2*j-1][2*i ][0]);	
					ptr[7] = &(pMat->RHS[2*k-1][2*j-1][2*i-1][0]);
			
					for(num=0; num<4; num++){
						pMatP->RHS[k][j][i][num] = (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0; 
					}
					
	
						
		}/* End if k, j, i */

		/* Also set ghost zones to be zero */
		
		for(k=pMatP->ks-Matghost; k<=pMatP->ke+Matghost; k++)
			for(j=pMatP->js-Matghost; j<=pMatP->je+Matghost; j++)
				for(i=pMatP->is-Matghost; i<=pMatP->ie+Matghost; i++){

				
					pMatP->U[k][j][i].Er = 0.0;
					pMatP->U[k][j][i].Fr1 = 0.0;
					pMatP->U[k][j][i].Fr2 = 0.0;
					pMatP->U[k][j][i].Fr3 = 0.0;

		}
	}/* End if this CPU works on this grid */
	}/* End if Level <= RootLevel */
/*==================================================================================*/
	else{
		/* We only work for Two levels at a time. When do restriction and receive, we should decide first, */ 
		/* Whether this CPU works in this grid */

		/* For MPI, first post non-block receiver */
		/* This will only post with MPI and for CPU work on this grid */
#ifdef MPI_PARALLEL
/* Post non-blocking receives at level Level-1 for data from child Grids at this
 * level (nl).  This data is sent in Step 3 below. */ 


    		for (nd=0; nd<pM->DomainsPerLevel[Level-RootLevel-1]; nd++){
      			if (pM->Domain[Level-RootLevel-1][nd].Grid != NULL) {
        			pG=pM->Domain[Level-RootLevel-1][nd].Grid;
				nZeroRC = 0;
        			mAddress = 0; 
				/* mAddress always starts from Zero, even if NmyCGrid  > 0 */
				/* This is because for grids on the same process, data is not stored in Recv buffer. */
				/* We just change the receive pointer to the send buffer */		
				
        			for (ncg=(pG->NmyCGrid); ncg<(pG->NCGrid); ncg++){
          				
					/* Only do this if we expect data to come */
					if(pG->CGrid[ncg].Rad_nWordsRC == 0){
						nZeroRC += 1;
					}
					else{
						mIndex = ncg - pG->NmyCGrid - nZeroRC;
          					ierr = MPI_Irecv(&(recv_bufRC[nd][mAddress]),
            					pG->CGrid[ncg].Rad_nWordsRC, MPI_DOUBLE, pG->CGrid[ncg].ID,
            					pG->CGrid[ncg].DomN, pM->Domain[Level-RootLevel-1][nd].Comm_Children,
            					&(recv_rq[nd][mIndex]));
          					mAddress += pG->CGrid[ncg].Rad_nWordsRC;
					}
					
        			}/* End all the child grid */

      			}/* End if the grid is null */
    		}/* End loop all domains at Level -1 */
 
#endif /* MPI_PARALLEL */		

		/* First, at Level, restrict and send the data, this will be the first step anyway */
		for(nd=0; nd<pM->DomainsPerLevel[Level-RootLevel]; nd++){
			/* If this CPU works for this Domain and this grid */
			pG = pM->Domain[Level-RootLevel][nd].Grid;
			if((Matrix[Level][nd].CPUflag) && (pG->NPGrid > 0)){
				 /* Level is gaurantee to be larger than RootLevel */
				pMat = &(Matrix[Level][nd]);
				start_addr = 0;
				nZeroRC = 0;
				/* Residual is calculated in multigrid main cycle. Here we assume residual is already calculated */


				/* There could be multiple parent grids overlap with this child grid */
				for(npg=0; npg<(pG->NPGrid); npg++){
					if(pG->PGrid[npg].Rad_nWordsRC == 0){
						nZeroRC += 1;
					}
					else{
					/* Only do  restriction if we actually have data to send, otherwise send_bufRC will not have enough space */
						pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
					
						/* Get coordinates ON THIS fine MATRIX of overlap region of parent MATRIX */
						/* If two grids just touch, *pe < *ps, there is no restriction or prolongation */

						ips = pPO->ijks[0] - nghost + Matghost;
						ipe = pPO->ijke[0] - nghost + Matghost;
						jps = pPO->ijks[1] - nghost + Matghost;
						jpe = pPO->ijke[1] - nghost + Matghost;
						kps = pPO->ijks[2] - nghost + Matghost;
						kpe = pPO->ijke[2] - nghost + Matghost;

				
						/* Now restrict coefficient, RHS and solution */
						pSnd = (double*)&(send_bufRC[nd][start_addr]);

						for(k=kps; k<=kpe; k+=2){
						for(j=jps; j<=jpe; j+=2){
						for(i=ips; i<=ipe; i+=2){
							if(flag){
						/* Restrict gas quantities to calculate the coefficients */
								ptr[0] = &(pMat->Ugas[k][j][i].rho);
								ptr[1] = &(pMat->Ugas[k][j][i+1].rho);
								ptr[2] = &(pMat->Ugas[k][j+1][i].rho);
								ptr[3] = &(pMat->Ugas[k][j+1][i+1].rho);
								ptr[4] = &(pMat->Ugas[k+1][j ][i].rho);
								ptr[5] = &(pMat->Ugas[k+1][j][i+1].rho);
								ptr[6] = &(pMat->Ugas[k+1][j+1][i].rho);	
								ptr[7] = &(pMat->Ugas[k+1][j+1][i+1].rho);

								for(num=0; num<11+NOPACITY; num++){

									pSnd[num] =  (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   		+ ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
								}
								/* move the pointer  */
								pSnd += (11 + NOPACITY);				
							}/* end if flag = 1, end restricting matrix coefficient */
						/* Now restrict the new RHS hand side */
						
						
							for(num=0; num<4; num++){
								pSnd[num] = (pMat->RHS[k][j][i][num] + pMat->RHS[k][j][i+1][num]
						     			+  pMat->RHS[k][j+1][i][num] + pMat->RHS[k][j+1][i+1][num]
						     			+  pMat->RHS[k+1][j][i][num] + pMat->RHS[k+1][j][i+1][num]
						     			+  pMat->RHS[k+1][j+1][i][num]+ pMat->RHS[k+1][j+1][i+1][num]) / 8.0;	
							} /* End restricting RHS */

							/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
							/* than coarse grid solution */
					
							/* move the pointer */
							pSnd += 4;

							ptr[0] = &(pMat->U[k][j][i].Er);
							ptr[1] = &(pMat->U[k][j][i+1].Er);
							ptr[2] = &(pMat->U[k][j+1][i].Er);
							ptr[3] = &(pMat->U[k][j+1][i+1].Er);
							ptr[4] = &(pMat->U[k+1][j ][i].Er);
							ptr[5] = &(pMat->U[k+1][j][i+1].Er);
							ptr[6] = &(pMat->U[k+1][j+1][i].Er);	
							ptr[7] = &(pMat->U[k+1][j+1][i+1].Er);
						
							for(num=0; num<4; num++){
								pSnd[num] = (ptr[0][num] + ptr[1][num]  + ptr[2][num] + ptr[3][num]
							 	   + ptr[4][num] + ptr[5][num]	+ ptr[6][num] + ptr[7][num]) / 8.0;
							}

							/* Move the pointer for Er, Fr? */
							pSnd += 4;

						}/* end ips */		
						}/* end jps */
						}/* end kps */
					
#ifdef MPI_PARALLEL
			/* send the data for MPI case */			
					 	if (npg >= pG->NmyPGrid ){
        						mIndex = npg - pG->NmyPGrid - nZeroRC;
        						ierr = MPI_Isend(&(send_bufRC[nd][start_addr]), pG->PGrid[npg].Rad_nWordsRC,
          							MPI_DOUBLE, pG->PGrid[npg].ID, nd, pM->Domain[Level-RootLevel][nd].Comm_Parent,
          							&(send_rq[nd][mIndex]));
     						}
#endif /* MPI_PARALLEL */
						/* set the start address prepared for next parent grid */
						start_addr += pG->PGrid[npg].Rad_nWordsRC;
					}/* End if Rad_nWordsRc != 0 */

				}/* End npg parent grid */				
			}/* End if Matrix[Level][nd].CPUflag */
		} /* End loop all the domains at this level */

		/* At the end, we need to wait for send from all domains to be finished at Level */

#ifdef MPI_PARALLEL

  		for(nd=0; nd<pM->DomainsPerLevel[Level-RootLevel]; nd++){
    			if ((pM->Domain[Level-RootLevel][nd].Grid) != NULL) {
      				pG = pM->Domain[Level-RootLevel][nd].Grid;
				nZeroRC = 0;
				for(i=0; i<pG->NPGrid; i++)
					if(pG->PGrid[i].Rad_nWordsRC == 0)	nZeroRC++;

      				if (pG->NPGrid > pG->NmyPGrid) {
        				mCount = pG->NPGrid - pG->NmyPGrid - nZeroRC;
        				ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      				}
    			}
  		}/* Finish all Domains at Level */
#endif /* MPI_PARALLEL */			


/*===================================================================*/
/* After send the data, now we need to receive the data for Level -1 */

		/* Get Child solution */
		for(nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel-1]); nd++){
			pG = pM->Domain[Level-1-RootLevel][nd].Grid; /* Level is gaurantee to be larger than RootLevel */
			if((Matrix[Level-1][nd].CPUflag) && (pG->NCGrid > 0)){				
				pMat = &(Matrix[Level-1][nd]);
				nZeroRC = 0;
				for(i=0; i< pG->NCGrid; i++)
					if(pG->CGrid[i].Rad_nWordsRC == 0)
						nZeroRC++;

#ifdef MPI_PARALLEL
					 if((RHStemp = (int**)calloc_2d_array(pG->NCGrid-nZeroRC,2,sizeof(double))) == NULL)
    						ath_error("[Restriction3D]:Failed to allocate RHStemp\n");
#endif
					
			
				for(ncg=0; ncg<(pG->NCGrid-nZeroRC); ncg++){
					if(ncg < pG->NmyCGrid){
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);
						/* send pointer to the beginning of the send buffer, if on the same process */
						pRcv = (double*)&(send_bufRC[pCO->DomN][0]);						

					} /* For the child grid on the same CPU */
					else {
#ifdef MPI_PARALLEL
						mCount = pG->NCGrid - pG->NmyCGrid - nZeroRC;
        					ierr = MPI_Waitany(mCount,recv_rq[nd],&mIndex,MPI_STATUS_IGNORE);
        					if(mIndex == MPI_UNDEFINED){
          						ath_error("[RestCorr]: Invalid request index nl=%i nd=%i\n",Level-1,nd);
        					}
      
					/* Recv buffer is addressed from 0 for first MPI message, even if NmyCGrid>0 */
        					
        					mIndex += pG->NmyCGrid;
						for(i=pG->NmyCGrid; i<mIndex; i++)
							if(pG->CGrid[i].Rad_nWordsRC == 0) mIndex++;

						mAddress = 0;
        					for (i=pG->NmyCGrid; i<mIndex; i++) mAddress += pG->CGrid[i].Rad_nWordsRC;
        						pCO=(GridOvrlpS*)&(pG->CGrid[mIndex]);
        						pRcv = (double*)&(recv_bufRC[nd][mAddress]);

						/* store the address information in RHStemp array */
						RHStemp[ncg][0] = mIndex;
						RHStemp[ncg][1] = mAddress;
#else
				/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        					ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",Level-1,nd);
#endif /* MPI_PARALLEL */

					}

					/* shift the position with ghost zone difference */
					/* Now the index corresponds to position in the matrix structure */

					

					ics = pCO->ijks[0] - nghost + Matghost;
      					ice = pCO->ijke[0] - nghost + Matghost;
      					jcs = pCO->ijks[1] - nghost + Matghost;
      					jce = pCO->ijke[1] - nghost + Matghost;
      					kcs = pCO->ijks[2] - nghost + Matghost;
      					kce = pCO->ijke[2] - nghost + Matghost;

					

					/* First, get the solution or matrix coefficient */
					/* The restricted RHS is stored at a temporary first */
					/* The restricted RHS is b - AX, which is the residual */
					/* We need to first calculate b_2 - AX_2 for the whole coarse level, as it is not completely covered by fine level */
					/* Then we replace the RHS with the restricted value */
					
					for(k=kcs; k<=kce; k++){
					for(j=jcs; j<=jce; j++){
					for(i=ics; i<=ice; i++){
						if(flag){
						/* Restrict gas quantities to calculate the coefficients */
							ptrP  = &(pMat->Ugas[k][j][i].rho);							

							for(num=0; num<11+NOPACITY; num++){
								ptrP[num] =  pRcv[num];
							}
							/* move the pointer  */
							pRcv += (11 + NOPACITY);				
						}/* end if flag = 1, end restricting matrix coefficient */
						/* Now restrict the new RHS hand side */
						
						/* Do not replace RHS right now */	
					/*	ptrP = &(pMat->RHS[k][j][i][0]);

						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];	
						}
					*/
						/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
						/* than coarse grid solution */
					
						/* move the pointer */
						pRcv += 4;

						ptrP = &(pMat->U[k][j][i].Er);
						
						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];
						}
						/* move the pointer */
						pRcv += 4;

					}/* end ics */		
					}/* end jcs */
					}/* end kcs */
					/* Now update the solution in the overlap region, update the RHS for this domain */
				}/* Loop over all the child grids */

				/* Now we need to update the RHS and solution for this level */

				/* We do not update RHS before calculate Residual for the whole level, this is required */
				/* For the un-refined region. Then update Residual for from fine levels */
				if(nDim == 3)
					RadSMR_Residual3D(pMat, pMat->RHS, &(pMat->RHSnorm));
				else
					ath_error("The function to calculate Residual only works for 3D now!");

				/* update RHS for the whole domain */
				/* We do not need ghost zones for RHS */
				for(k=pMat->ks; k<=pMat->ke; k++)
				for(j=pMat->js; j<=pMat->je; j++)
				for(i=pMat->is; i<=pMat->ie; i++){
					
					/* Now update the solution */
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Er  += pMat->U[k][j][i].Er;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr1 += pMat->U[k][j][i].Fr1;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr2 += pMat->U[k][j][i].Fr2;
					pG->U[k-Matghost+nghost][j-Matghost+nghost][i-Matghost+nghost].Fr3 += pMat->U[k][j][i].Fr3;
				
				}/* end i, j, k */


				/* Now we need to replace the RHS in the overlap region with restricted RHS from parent grids */
				/* The restricted data is already in the recv_buf, so we do not need to wait for MPI again */
 
		/*===================================================================================*/
				for(ncg=0; ncg<(pG->NCGrid-nZeroRC); ncg++){
					if(ncg < pG->NmyCGrid){
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);
						/* send pointer to the beginning of the send buffer, if on the same process */
						pRcv = (double*)&(send_bufRC[pCO->DomN][0]);						

					} /* For the child grid on the same CPU */
					else {
#ifdef MPI_PARALLEL
       						pCO=(GridOvrlpS*)&(pG->CGrid[RHStemp[ncg][0]]);
       						pRcv = (double*)&(recv_bufRC[nd][RHStemp[ncg][1]]);
#else
				/* If not MPI_PARALLEL, and child Grid not on this processor, then error */

        					ath_error("[RestCorr]: no Child grid on Domain[%d][%d]\n",Level-1,nd);
#endif /* MPI_PARALLEL */

					}

					/* shift the position with ghost zone difference */
					/* Now the index corresponds to position in the matrix structure */
					ics = pCO->ijks[0] - nghost + Matghost;
      					ice = pCO->ijke[0] - nghost + Matghost;
      					jcs = pCO->ijks[1] - nghost + Matghost;
      					jce = pCO->ijke[1] - nghost + Matghost;
      					kcs = pCO->ijks[2] - nghost + Matghost;
      					kce = pCO->ijke[2] - nghost + Matghost;

					/* First, get the solution or matrix coefficient */
					/* The restricted RHS is stored at a temporary first */
					
					for(k=kcs; k<=kce; k++){
					for(j=jcs; j<=jce; j++){
					for(i=ics; i<=ice; i++){
						if(flag){
							pRcv += (11 + NOPACITY);				
						}/* end if flag = 1, end restricting matrix coefficient */
						
						/* Now replace the restrcited RHS */
						
						
						ptrP = &(pMat->RHS[k][j][i][0]);

						for(num=0; num<4; num++){
							ptrP[num] = pRcv[num];	
						}
					
						/* Always restrict the solution, as we assume fine grid solution is ALWAYS better */ 
						/* than coarse grid solution */
					
						/* move the pointer */
						pRcv += 8;

					}/* end ics */		
					}/* end jcs */
					}/* end kcs */
					/* Now update the solution in the overlap region, update the RHS for this domain */
				}/* Loop over all the child grids */

		/*==================================================================================*/
				/* Free the temporary array */
#ifdef MPI_PARALLEL
					 if(RHStemp != NULL)
    						free_2d_array(RHStemp);
#endif

			}/* End if [Level-1][nd] CPU flag */
		}/* End loop all domains at Level-1 */


	}/* End for levels above RootLevel */




	return;
}


/* prolongation operator between Level and Level + 1 */
/* In the prolongation step, for levels above root, */



void Rad_Prolongation(MatrixS **Matrix, const int Level, const MeshS *pM)
{



	
  	int i, j, k, n, m, l, is, ie, js, je, ks, ke, ii, jj, kk, kfine, jfine, ifine;
	int ics, ice, jcs, jce, kcs, kce;
	int nd, nZeroP, npg, ncg, ngz1, ngz2, ngz3, igzs, igze, jgzs, jgze, kgzs, kgze, mend, nend, lend;
	int nDim, dim, id;
	int ips, ipe, jps, jpe, kps, kpe;
	

	double *pRcv, *pSnd;
  	
	GridOvrlpS *pCO, *pPO;
	MatrixS *pMat, *pMatC;	
	RadMHDS Ptemp[2][2][2];	
	GridS *pG;
#ifdef MPI_PARALLEL
  	int mAddress, ierr, mIndex, mCount;
#endif

	/* number of dimensions in Grid. */
	/* First, determine the dimensionality */
  	nDim=1;
  	for (i=1; i<3; i++) if (pM->Nx[i]>1) nDim++;

	/*===========================================================================*/
	/* We only have two levels to deal with at each time */
	/* First, we need to judge whether the Level if below the RootLevel or not */ 
	/* Below the root level, we use the original way */
	/* Only do this if this CPU works in this grid */
	/* Grids below the root work under the same CPU */
	if(Level < RootLevel){
	if(Matrix[RootLevel][0].CPUflag){
		pMat  = &(Matrix[Level][0]); /* Matrix in the coarse Level */
		pMatC = &(Matrix[Level+1][0]);/* Matrix in the fine Level */

		is = pMat->is;
		ie = pMat->ie;
		js = pMat->js;
		je = pMat->je;
		ks = pMat->ks;
		ke = pMat->ke;

		for(k=ks; k<=ke; k++)
			for(j=js; j<=je; j++)
				for(i=is; i<=ie; i++){

					ProU(pMat->U[k][j][i-1],pMat->U[k][j][i],pMat->U[k][j][i+1],pMat->U[k][j-1][i],pMat->U[k][j+1][i],
						pMat->U[k-1][j][i],pMat->U[k+1][j][i],Ptemp);


			/* Now copy the data to the fine grid */
					ii = 2*(i-is) + pMatC->is;
					jj = 2*(j-js) + pMatC->js;
					kk = 2*(k-ks) + pMatC->ks;

				/* The coarse grid calculates the residual, but solution is already added to grid solution */
					for(kfine=0;kfine<2; kfine++)
						for(jfine=0; jfine<2; jfine++)
							for(ifine=0; ifine<2; ifine++){
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Er  = Ptemp[kfine][jfine][ifine].Er;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr1 = Ptemp[kfine][jfine][ifine].Fr1;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr2 = Ptemp[kfine][jfine][ifine].Fr2;
								pMatC->U[kk+kfine][jj+jfine][ii+ifine].Fr3 = Ptemp[kfine][jfine][ifine].Fr3;
					}/* Finish copy solution back to child grid */				

		}/* Finish loop the whole grid i, j, k */
	}/* End if this CPU works in this grid */
	}/* End if the Level is below the RootLevel */
	else{
		/* First, post non-blocking receives at Level + 1 for data from parent GridS */


#ifdef MPI_PARALLEL

		
    		for (nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel+1]); nd++){
      			if (pM->Domain[Level+1-RootLevel][nd].Grid != NULL) {
        			
				pG=pM->Domain[Level+1-RootLevel][nd].Grid;
        			nZeroP = 0;
        			mAddress = 0;        			
        			if (pG->NmyPGrid > 0) mAddress = pG->PGrid[0].Rad_nWordsP;

        		for (npg=(pG->NmyPGrid); npg<(pG->NPGrid); npg++){

			/* Skip if no prolongation needed for this child (only flux correction) */ 
			/* This is for the case when grids only touch */
          			if (pG->PGrid[npg].Rad_nWordsP == 0) { 
            				nZeroP += 1;
          			} else {

           				mIndex = npg - pG->NmyPGrid - nZeroP;
            				ierr = MPI_Irecv(&(recv_bufP[nd][mAddress]),
              				pG->PGrid[npg].Rad_nWordsP, MPI_DOUBLE, pG->PGrid[npg].ID,
              				pG->PGrid[npg].DomN, pM->Domain[Level+1-RootLevel][nd].Comm_Parent,
              					&(recv_rq[nd][mIndex]));
            				mAddress += pG->PGrid[npg].Rad_nWordsP;
          			}

        		}/* End loop all the parent grids */
      			}/* End if Grid != NULL */
    		}/* Finish looping all domains at Level +1 */
  
#endif /* MPI_PARALLEL */



		/*======Step 1, send the data, including the ghoze zones for Child grids==========*/
		/* Unlik the SMR in normal MHD part, we need to send the whole grid data */
		/* Now we works for Level */
		/* We need to send the data including ghost and active zones */
		/* But ghost zones will only be prolongated if this boundary is interior of the coarse grid */

		for(nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel]); nd++){
			pG = pM->Domain[Level-RootLevel][nd].Grid;
			if((Matrix[Level][nd].CPUflag) && (pG->NCGrid > 0)){
				/* Check that there is child grid at this level, otherwise we shouldn't come here. */
				

				pMat = &(Matrix[Level][nd]);

				for(i=0; i<maxND; i++) start_addrP[i] = 0;
				nZeroP = 0;
				
				for(ncg=0; ncg<(pG->NCGrid); ncg++){
					if(pG->CGrid[ncg].Rad_nWordsP == 0){
						nZeroP += 1;
					}/* skip the grid that does need prolongation */
					else{
						pCO = (GridOvrlpS*)&(pG->CGrid[ncg]);	/* ptr to child Grid overlap */

					/* index send_buf with DomN of child, since could be multiple child Domains on */
 					/* same processor.  Start address must be different for each DomN */

						/* First prolongate the region that overlaps */
						
						pSnd = (double*)&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]);

						ics = pCO->ijks[0] - nghost + Matghost;
						ice = pCO->ijke[0] - nghost + Matghost;
						jcs = pCO->ijks[1] - nghost + Matghost;
						jce = pCO->ijke[1] - nghost + Matghost;
						kcs = pCO->ijks[2] - nghost + Matghost;
						kce = pCO->ijke[2] - nghost + Matghost;
					
						
						for(k=kcs; k<=kce; k++){
						for(j=jcs; j<=jce; j++){
						for(i=ics; i<=ice; i++){
							*(pSnd++) = pMat->U[k][j][i].Er;
							*(pSnd++) = pMat->U[k][j][i].Fr1;
							*(pSnd++) = pMat->U[k][j][i].Fr2;
							*(pSnd++) = pMat->U[k][j][i].Fr3;

						}/* end i */
						}/* end j */
						}/* end k */
						
						/* Now send the data for the necessary boundary */
						/* This only works for Matghost == 1 */
						for(dim=0; dim<(2*nDim); dim++){
							if(pCO->AdvEr[dim]){
								ics = pCO->ijks[0] - nghost + Matghost - 1;
								ice = pCO->ijke[0] - nghost + Matghost + 1;
						
								if(nDim > 1){
									jcs = pCO->ijks[1] - nghost + Matghost - 1;
									jce = pCO->ijke[1] - nghost + Matghost + 1;
								}
								else{
									jcs = pCO->ijks[1];
									jce = pCO->ijke[1];
								}		

								if(nDim > 2){
									kcs = pCO->ijks[2] - nghost + Matghost - 1;
									kce = pCO->ijke[2] - nghost + Matghost + 1;
								}
								else{
									kcs = pCO->ijks[2];
									kce = pCO->ijke[2];
								}


								if (dim == 0) (ice = pCO->ijks[0] - nghost + Matghost);
          							if (dim == 1) (ics = pCO->ijke[0] - nghost + Matghost);
          							if (dim == 2) (jce = pCO->ijks[1] - nghost + Matghost);
          							if (dim == 3) (jcs = pCO->ijke[1] - nghost + Matghost);
         							if (dim == 4) (kce = pCO->ijks[2] - nghost + Matghost);
          							if (dim == 5) (kcs = pCO->ijke[2] - nghost + Matghost);

	
								for(k=kcs; k<=kce; k++){
								for(j=jcs; j<=jce; j++){
								for(i=ics; i<=ice; i++){
									*(pSnd++) = pMat->U[k][j][i].Er;
									*(pSnd++) = pMat->U[k][j][i].Fr1;
									*(pSnd++) = pMat->U[k][j][i].Fr2;
									*(pSnd++) = pMat->U[k][j][i].Fr3;

								}/* end i */
								}/* end j */
								}/* end k */


							}/* This dim is a boundary */
						}/* end 6 faces */



			/* Step 1b: non-blocking send of data  to Child, using Domain number as tag */
#ifdef MPI_PARALLEL
				/* Only do this if we actually have data to send */
      						if (ncg >= pG->NmyCGrid) {
        						mIndex = ncg - pG->NmyCGrid - nZeroP;
        						ierr = MPI_Isend(&(send_bufP[pCO->DomN][start_addrP[pCO->DomN]]),
          						pG->CGrid[ncg].Rad_nWordsP, MPI_DOUBLE, pG->CGrid[ncg].ID, nd,
          						pM->Domain[Level-RootLevel][nd].Comm_Children, &(send_rq[nd][mIndex]));
      						}
#endif /* MPI_PARALLEL */

      						start_addrP[pCO->DomN] += pG->CGrid[ncg].Rad_nWordsP;
					}/* End for the grids that needs prolongation */
				}/* End loop over the child grid */

			}/* End if CPUflag */
		}/* Finish looping all the domains at Level */



		/*==============================================================*/
		/* Because data is sent from Level, we need to clear the send_bufP for Level */
		/* Unlike the normal smr function for MHD, we only have two levels here and we do not need to loop over level */
		/* So for grids on the same CPU, we need set recv pointer to send buffer */

  	/******************************************************************************/
	/* This step is skipped as we just set Receive pointer to the send buffer if grids are on the same CPU 
	
		for (nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel]); nd++){
    			if (pM->Domain[Level-RootLevel][nd].Grid != NULL) { 
     				 pG=pM->Domain[Level-RootLevel][nd].Grid; 
	

      				for (ncg=0; ncg<(pG->NmyCGrid); ncg++){
        				pCO=(GridOvrlpS*)&(pG->CGrid[ncg]);   

        				for (i=0; i<pCO->Rad_nWordsP; i++) {
          					recv_bufP[pCO->DomN][i]=send_bufP[pCO->DomN][i];
        				}
      				}
    			}
  		}

 */

#ifdef MPI_PARALLEL
		/* For MPI jobs, wait for all non-blocking sends above to finish in order to continue to Level+1  */
		/* Ortherwise, the CPU may start to working while sent is not complete */

  		for (nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel]); nd++){
    			if (pM->Domain[Level-RootLevel][nd].Grid != NULL) {
     				pG=pM->Domain[Level-RootLevel][nd].Grid;

      				nZeroP = 0;
      				for (i=0; i < pG->NCGrid; i++) if (pG->CGrid[i].Rad_nWordsP == 0) nZeroP++;

      				if (pG->NCGrid > pG->NmyCGrid) {
        				mCount = pG->NCGrid - pG->NmyCGrid - nZeroP;
        				ierr = MPI_Waitall(mCount, send_rq[nd], MPI_STATUS_IGNORE);
      				}
    			}
  		}
#endif /* MPI_PARALLEL */



		/*=====================================================================*/
		/* Get solution from parent GridS and prolongation solution to ghost zones */
		/* Now we go back to Level + 1 */
		for(nd=0; nd<(pM->DomainsPerLevel[Level-RootLevel+1]); nd++){
			if(pM->Domain[Level+1-RootLevel][nd].Grid != NULL){
				pG = pM->Domain[Level+1-RootLevel][nd].Grid;
				pMat = &(Matrix[Level+1][nd]);

				/* Loop over number of parent grids with non-zero-size prolongation data */
				nZeroP = 0;
    				for (i=0; i < pG->NPGrid; i++) if (pG->PGrid[i].Rad_nWordsP == 0) nZeroP++;

				for (npg=0; npg<(pG->NPGrid - nZeroP); npg++){

				/* If parent Grid is on this processor, data is at start of recv buffer */

      					if (npg < pG->NmyPGrid) {
        					pPO = (GridOvrlpS*)&(pG->PGrid[npg]);
						/* For grids on the same CPU, set the pointer to send buffer */
        					/* pRcv = (double*)&(recv_bufP[nd][0]); */
						/* One CPU only works for one grid in one Domain */
						/* So, the address at buffer for domain nd must start from 0 */
						pRcv = (double*)&(send_bufP[nd][0]);
	
      					} else {

#ifdef MPI_PARALLEL
					/* Check non-blocking receives posted above for data in ghost zone from parent
 						* Grids, sent in Step 1.  Accept messages in any order. */

        					mCount = pG->NPGrid - pG->NmyPGrid - nZeroP;
        					ierr = MPI_Waitany(mCount,recv_rq[nd],&mIndex,MPI_STATUS_IGNORE);
        					if(mIndex == MPI_UNDEFINED){
          						ath_error("[Prolong]: Invalid request index nl=%i nd=%i\n",Level+1,nd);
        					}

						/* mIndex returns the number that is completed */

					/* Recv buffer is addressed from PGrid[0].Rad_nWordsP for first MPI message
 					* if NmyPGrid>0.  Also must remove zero size messages from index. */

						/* re-build the value of mIndex to include NmyPGrid and grids that do not need prolongation */

        					mIndex += pG->NmyPGrid;
        					for (i=pG->NmyPGrid; i <= mIndex; i++) 
          						if (pG->PGrid[i].Rad_nWordsP == 0) mIndex++;

        					mAddress = 0;
        					for (i=0; i<mIndex; i++) mAddress += pG->PGrid[i].Rad_nWordsP;
        					pPO = (GridOvrlpS*)&(pG->PGrid[mIndex]); 
        					pRcv = (double*)&(recv_bufP[nd][mAddress]);
#else
				/* If not MPI_PARALLEL, and parent Grid not on this processor, then error */

        					ath_error("[Prolong]: no Parent Grid on Domain[%d][%d]\n",Level+1,nd);
#endif /* MPI_PARALLEL */
      					}/* End if npg > pG->NmyPGrid */

					
					/* Loop over 6 boundaries, set ghost zones */
					/* The difference between Matrix solver and normal MHD is that */
					/* We also need to prolongate data from each cell to the child grids */
					
					/* Only do this if prolongation data is non-zero */
					/* This is especially for NmyPgrid */
					if(pPO->Rad_nWordsP > 0){

					/* Get coordinates ON THIS GRID of ghost zones that overlap parent Grid */
						/* first prolongate the overlap region */
						/* To use this prolongated solution as initial guess in fine level */

						/*---------------------------------------------------*/
						ips = 0;
						ipe = (pPO->ijke[0] - pPO->ijks[0] + 1)/2 - 1;
						if(nDim > 1){
							jps = 0;
							jpe = (pPO->ijke[1] - pPO->ijks[1] + 1)/2 - 1;
						}
						else{
							jps = 1;
							jpe = 1;
						}
						if(nDim > 2){
							kps = 0;
							kpe = (pPO->ijke[2] - pPO->ijks[2] + 1)/2 - 1;
						}
						else{
							kps = 1;
							kpe = 1;
						}
						/* Load the data */
					
						for(k=kps; k<=kpe; k++){
						for(j=jps; j<=jpe; j++){
						for(i=ips; i<=ipe; i++){
							Pro_buf[k][j][i].Er = *(pRcv++);
							Pro_buf[k][j][i].Fr1 = *(pRcv++);
							Pro_buf[k][j][i].Fr2 = *(pRcv++);
							Pro_buf[k][j][i].Fr3 = *(pRcv++);

						}/* End i */
						}/* End j */
						}/* End k */

						/* Fill the junk zones for 1D and 2D cases */
						if(nDim == 1){
							for(i=ips; i<=ipe; i++){
								Pro_buf[1][0][i] = Pro_buf[1][1][i];
								Pro_buf[1][2][i] = Pro_buf[1][1][i];
								Pro_buf[0][1][i] = Pro_buf[1][1][i];
								Pro_buf[2][1][i] = Pro_buf[1][1][i];
							}

						}else if(nDim == 2){
							for(j=jps; j<=jpe; j++){
							for(i=ips; i<=ipe; i++){
								Pro_buf[0][j][i] = Pro_buf[1][j][i];
								Pro_buf[2][j][i] = Pro_buf[1][j][i];
							}
							}

						}/* End if nDim == 2 */


						/* Do prolongation */
						ips = pPO->ijks[0] - nghost + Matghost;
						ipe = pPO->ijke[0] - nghost + Matghost;
						lend = 1;
						if(nDim > 1){
							jps = pPO->ijks[1] - nghost + Matghost;
							jpe = pPO->ijke[1] - nghost + Matghost;
							mend = 1;
						}
						else{
							jps = pPO->ijks[1];
							jpe = pPO->ijke[1];
							mend = 0;
						}

						if(nDim > 2){
							kps = pPO->ijks[2] - nghost + Matghost;
							kpe = pPO->ijke[2] - nghost + Matghost;
							nend = 1;
						}
						else{
							kps = pPO->ijks[2];
							kpe = pPO->ijke[2];
							nend = 0;
						}

							/* Prolongate the Pro_buf array data to Ptemp temporarily and then copy the data to pMat */
							/* i, j, k are for parent grids while kk, jj, ii are for child grids */
						for (k=kps, kk=1; k<=kpe; k+=2, kk++) {
						for (j=jps, jj=1; j<=jpe; j+=2, jj++) {
       						for (i=ips, ii=1; i<=ipe; i+=2, ii++) {
							ProU(Pro_buf[kk][jj][ii-1], Pro_buf[kk][jj][ii], Pro_buf[kk][jj][ii+1], Pro_buf[kk][jj-1][ii], 
									Pro_buf[kk][jj+1][ii], Pro_buf[kk-1][jj][ii], Pro_buf[kk+1][jj][ii], Ptemp);

							/* Now set the solution */
							for(n=0; n<=nend; n++){
							for(m=0; m<=mend; m++){
							for(l=0; l<=1; l++){
								pMat->U[k+n][j+m][i+l].Er  = Ptemp[n][m][l].Er;
								pMat->U[k+n][j+m][i+l].Fr1 = Ptemp[n][m][l].Fr1;
								pMat->U[k+n][j+m][i+l].Fr2 = Ptemp[n][m][l].Fr2;
								pMat->U[k+n][j+m][i+l].Fr3 = Ptemp[n][m][l].Fr3;
							}/* End l */
							}/* End m */
							}/* End n */

						}/* End i */
						}/* End j */
						}/* End k */


						/*-----------------------------------------------------*/
						/* Now set ghost zones for required faces */
						for(dim=0; dim<(2*nDim); dim++){
							if(pPO->AdvEr[dim]){

								if(dim == 0 || dim == 1){
									ngz1 = (Matghost/2) + 2;
									id = 0;
								}else{
									ngz1 = (pPO->ijke[0] - pPO->ijks[0] + 1)/2 + 2;
								}

								if(dim == 2 || dim == 3){
									ngz2 = (Matghost/2) + 2;
									id = 1;
								}else{
									ngz2 = (pPO->ijke[1] - pPO->ijks[1] + 1)/2 + 2;
								}

								if(dim == 4 || dim == 5){
									ngz3 = (Matghost/2) + 2;
									id = 2;
								}else{
									ngz3 = (pPO->ijke[2] - pPO->ijks[2] + 1)/2 + 2;
								}

								igzs = 0;
								igze = ngz1 - 1;
				
								if (nDim > 1) {
									jgzs = 0;
									jgze = ngz2 - 1;							
								} else {
									ngz2 = 1;
									jgzs = 1;
									jgze = 1;							
								}
								if (nDim > 2) {
									kgzs = 0;
									kgze = ngz3 - 1;							
								} else {
									ngz3 = 1;
									kgzs = 1;
									kgze = 1;							
								}

						/* Load the data */
					
								for(k=kgzs; k<=kgze; k++){
								for(j=jgzs; j<=jgze; j++){
								for(i=igzs; i<=igze; i++){
									GZ[id][k][j][i].Er = *(pRcv++);
									GZ[id][k][j][i].Fr1 = *(pRcv++);
									GZ[id][k][j][i].Fr2 = *(pRcv++);
									GZ[id][k][j][i].Fr3 = *(pRcv++);
	
								}/* End i */
								}/* End j */
								}/* End k */

						/* Set boundary conditions for GZ in 1D and 2D cases */
						/* This is needed for the prolongation below */
								 if (nDim == 1) {
            								for (i=igzs; i<=igze; i++) {
										GZ[id][1][0][i] = GZ[id][1][1][i];
              									GZ[id][1][2][i] = GZ[id][1][1][i];
              									GZ[id][0][1][i] = GZ[id][1][1][i];
              									GZ[id][2][1][i] = GZ[id][1][1][i];
							
            								}/* End for i*/
          							}/* End if nDim = 1 */
								else if (nDim == 2) {
            								for (j=jgzs; j<=jgze; j++) {
            								for (i=igzs; i<=igze; i++) {
										GZ[id][0][j][i] = GZ[id][1][j][i];
              									GZ[id][2][j][i] = GZ[id][1][j][i];								
            								}/* End for i */
									}/* End for j */
								}/* End if nDim = 2*/


						/* Now prolongate the ghost zones in array GZ */

								ips = pPO->ijks[0] - nghost;	/* This actually is - nghost + Matghost - Matghost */
								ipe = pPO->ijke[0] - nghost + Matghost + Matghost;
								
								if(nDim > 1){
									jps = pPO->ijks[1] - nghost;
									jpe = pPO->ijke[1] - nghost + Matghost + Matghost;
								}else{
									jps = pPO->ijks[1]; /* The value will be zero in this case */
									jpe = pPO->ijke[1];
								}

								if(nDim > 2){
									kps = pPO->ijks[2] - nghost;
									kpe = pPO->ijke[2] - nghost + Matghost + Matghost;
								}else{
									kps = pPO->ijks[2]; /* The value will be zero in this case */
									kpe = pPO->ijke[2];
								}

								
								/* In case Matghost=1, we only need to prolongate one cell */
								if (dim == 0) {
									ipe = pPO->ijks[0] - 1 - nghost + Matghost;
									lend = MIN((Matghost-1),1);
								}
          							if (dim == 1) {
									ips = pPO->ijke[0] + 1 - nghost + Matghost;
									lend = MIN((Matghost-1),1);
								}
		
          							if (dim == 2) {
									jpe = pPO->ijks[1] - 1 - nghost + Matghost;
									mend = MIN((Matghost-1),mend);
								}
          							if (dim == 3) {
									jps = pPO->ijke[1] + 1 - nghost + Matghost;
									mend = MIN((Matghost-1),mend);
								}
          							if (dim == 4) {
									kpe = pPO->ijks[2] - 1 - nghost + Matghost;
									nend = MIN((Matghost-1),nend);
								}
          							if (dim == 5) {
									kps = pPO->ijke[2] + 1 - nghost + Matghost;
									nend = MIN((Matghost-1),nend);
								}

							/* Prolongate the GZ array data to Ptemp temporarily and then copy the data to pMat */
							/* i, j, k are for parent grids while kk, jj, ii are for child grids */
						

								for (k=kps, kk=1; k<=kpe; k+=2, kk++) {
								for (j=jps, jj=1; j<=jpe; j+=2, jj++) {
       								for (i=ips, ii=1; i<=ipe; i+=2, ii++) {
									ProU(GZ[id][kk][jj][ii-1], GZ[id][kk][jj][ii], GZ[id][kk][jj][ii+1], GZ[id][kk][jj-1][ii], 
									GZ[id][kk][jj+1][ii], GZ[id][kk-1][jj][ii], GZ[id][kk+1][jj][ii], Ptemp);

									/* Now set the solution */
									for(n=0; n<=nend; n++){
									for(m=0; m<=mend; m++){
									for(l=0; l<=lend; l++){
										pMat->U[k+n][j+m][i+l].Er  = Ptemp[n][m][l].Er;
										pMat->U[k+n][j+m][i+l].Fr1 = Ptemp[n][m][l].Fr1;
										pMat->U[k+n][j+m][i+l].Fr2 = Ptemp[n][m][l].Fr2;
										pMat->U[k+n][j+m][i+l].Fr3 = Ptemp[n][m][l].Fr3;
									}/* End l */
									}/* End m */
									}/* End n */

								}/* End i */
								}/* End j */
								}/* End k */

							}/* End if dim */
						}/* Finish looping all six faces */

					} /* End if there are data to prolongate */
				}/* End loop npg grid */



			}/* End if [Level+1][nd].Grid is not NULL */
		}/* Finish loop over domains at Level + 1*/

		/*===================================================================*/

	}/* End if the Level is above the RootLevel */







}





/* Prolongation scheme adopted from SMR */
/* Only need to prolongate Er, Fr1,Fr2, Fr3 four variables */
/* from coarse grid to fine grid */
/* PU is used to take the data out */

void ProU(const RadMHDS Uim1, const RadMHDS Ui, const RadMHDS Uip1, 
	  const RadMHDS Ujm1, const RadMHDS Ujp1,
	  const RadMHDS Ukm1, const RadMHDS Ukp1, RadMHDS PCon[2][2][2])
{
  int i,j,k;

  Real dq1,dq2,dq3;


/* First order prolongation -- just copy values */
#ifdef FIRST_ORDER

  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Er  = Ui.Er;
    PCon[k][j][i].Fr1 = Ui.Fr1;
    PCon[k][j][i].Fr2 = Ui.Fr2;
    PCon[k][j][i].Fr3 = Ui.Fr3;
  }}}
/* second order prolongation -- apply limited slope reconstruction */
#else /* SECOND_ORDER or THIRD_ORDER */

/* Er */
  dq1 = mcd_slope(Uim1.Er, Ui.Er, Uip1.Er);
  dq2 = mcd_slope(Ujm1.Er, Ui.Er, Ujp1.Er);
  dq3 = mcd_slope(Ukm1.Er, Ui.Er, Ukp1.Er);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Er  = Ui.Er 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr1 */
  dq1 = mcd_slope(Uim1.Fr1, Ui.Fr1, Uip1.Fr1);
  dq2 = mcd_slope(Ujm1.Fr1, Ui.Fr1, Ujp1.Fr1);
  dq3 = mcd_slope(Ukm1.Fr1, Ui.Fr1, Ukp1.Fr1);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr1 = Ui.Fr1 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr2 */
  dq1 = mcd_slope(Uim1.Fr2, Ui.Fr2, Uip1.Fr2);
  dq2 = mcd_slope(Ujm1.Fr2, Ui.Fr2, Ujp1.Fr2);
  dq3 = mcd_slope(Ukm1.Fr2, Ui.Fr2, Ukp1.Fr2);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr2 = Ui.Fr2 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}

/* Fr3 */
  dq1 = mcd_slope(Uim1.Fr3, Ui.Fr3, Uip1.Fr3);
  dq2 = mcd_slope(Ujm1.Fr3, Ui.Fr3, Ujp1.Fr3);
  dq3 = mcd_slope(Ukm1.Fr3, Ui.Fr3, Ukp1.Fr3);
  for (k=0; k<2; k++){
  for (j=0; j<2; j++){
  for (i=0; i<2; i++){
    PCon[k][j][i].Fr3 = Ui.Fr3 
      + (0.5*i - 0.25)*dq1 + (0.5*j - 0.25)*dq2 + (0.5*k - 0.25)*dq3;
  }}}
#endif /* FIRST_ORDER */
}




#ifndef FIRST_ORDER
static Real mcd_slope(const Real vl, const Real vc, const Real vr){

  Real dvl = (vc - vl), dvr = (vr - vc);
  Real dv, dvm;

  if(dvl > 0.0 && dvr > 0.0){
    dv = 2.0*(dvl < dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm < dv ? dvm : dv);
  }
  else if(dvl < 0.0 && dvr < 0.0){
    dv = 2.0*(dvl > dvr ? dvl : dvr);
    dvm = 0.5*(dvl + dvr);
    return (dvm > dv ? dvm : dv);
  }

  return 0.0;
}
#endif /* FIRST_ORDER */





void SMR_Rad_init(MeshS *pM, const int Root)
{

	int nl, nd, DomainLevels, npg, ncg;

	GridS *pG;
	
	/* Root is the number of Levels below the root */
	/* RootLevel is a local global variable */
	RootLevel = Root;

/*=========================================================*/
	int maxCG=1, max1=0, max2=0, max3=0;	
	int sendRC, recvRC, sendP, recvP;
	int max_sendRC=1, max_recvRC=1, max_sendP=1, max_recvP=1; /* maximum size of send buffer */
	
/*===========================================================*/
	/* first, calculate DomainLevels and RootLevel */
	DomainLevels = pM->NLevels;

	maxND = 1;
	for(nl=0; nl<DomainLevels; nl++)
		maxND = MAX(maxND,pM->DomainsPerLevel[nl]);

	if((start_addrP = (int*)calloc_1d_array(maxND,sizeof(int))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate start_addrP\n");


	/*==========================================================*/
	/* Allocate memory and initialize parameters for restriction and prolongation */
	/* The maximum memory size required to send data from parent to child grids */
	for(nl=0; nl<DomainLevels; nl++){
		for(nd=0; nd<pM->DomainsPerLevel[nl]; nd++){
			sendRC = 0;
			recvRC = 0;
			sendP = 1;
			recvP = 1;
			if(pM->Domain[nl][nd].Grid != NULL){
				pG = pM->Domain[nl][nd].Grid;				
				/* For parent grids; send data for all grids per domain together*/
				for(npg=0; npg<pG->NPGrid; npg++){
					sendRC += pG->PGrid[npg].Rad_nWordsRC;
					recvP  += pG->PGrid[npg].Rad_nWordsP;
				}
				/* For child grids; do all grids per domain together */
				for(ncg=0; ncg<pG->NCGrid; ncg++){
					recvRC += pG->CGrid[ncg].Rad_nWordsRC;
					sendP  += pG->CGrid[ncg].Rad_nWordsP;
				}
			
				max_sendRC = MAX(max_sendRC, sendRC);
				max_recvRC = MAX(max_recvRC, recvRC);
				max_sendP  = MAX(max_sendP, sendP);
				max_recvP  = MAX(max_recvP, recvP);
				maxCG = MAX(maxCG,pG->NCGrid);
				max1  = MAX(max1, (pG->Nx[0]+1));
				max2  = MAX(max2, (pG->Nx[1]+1));
				max3  = MAX(max3, (pG->Nx[2]+1));
			}/* End if grid is not NULL */

		}/* End domain at each level */
	}/* End level */

	/*======================================================================*/
	/* Allocate memory for variables used for restriction */

	if((send_bufRC = (double**)calloc_2d_array(maxND,max_sendRC,sizeof(double))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate send_bufRC\n");


#ifdef MPI_PARALLEL
	/* we do restriction between every two levels; We will finish one level first and then */
	/* go to another level */ 
	if((recv_bufRC =  (double**)calloc_2d_array(maxND,max_recvRC,sizeof(double))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv_bufRC\n");
	/* We only work for two Levels at a time, so the first index is 2, not total Levels */
	if((recv_rq = (MPI_Request**)calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv MPI_Request array\n");
  	if((send_rq = (MPI_Request**)calloc_2d_array(maxND,maxCG,sizeof(MPI_Request))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate send MPI_Request array\n");
#endif /* MPI_PARALLEL */

	
	/*======================================================================*/
	
	if((send_bufP =(double**)calloc_2d_array(maxND,max_sendP,sizeof(double))) == NULL)
    		ath_error("[SMR_init]:Failed to allocate send_bufP\n");

	if((recv_bufP =(double**)calloc_2d_array(maxND,max_recvP,sizeof(double))) == NULL)
    		ath_error("[SMR_init]: Failed to allocate recv_bufP\n");


	max1 += 2*Matghost;
  	max2 += 2*Matghost;
  	max3 += 2*Matghost;

	/* Array to store the prolongation data temporary */
	/* Each GZ[k][j][i][Er-Fr?]: We only need to do the boundary for Er to Fr? */
	if((GZ[0]=(RadMHDS***)calloc_3d_array(max3,max2,nghost,sizeof(RadMHDS))) ==NULL) 
		ath_error("[SMR_Rad_init]:Failed to allocate GZ[0]C\n");

	if((GZ[1]=(RadMHDS***)calloc_3d_array(max3,nghost,max1,sizeof(RadMHDS))) ==NULL) 
		ath_error("[SMR_Rad_init]:Failed to allocate GZ[0]C\n");

	if((GZ[2]=(RadMHDS***)calloc_3d_array(nghost,max2,max1,sizeof(RadMHDS))) ==NULL) 
		ath_error("[SMR_Rad_init]:Failed to allocate GZ[0]C\n");
  	
  	if((Pro_buf=(RadMHDS***)calloc_3d_array((max3-2*Matghost),(max2-2*Matghost),(max1-2*Matghost),sizeof(RadMHDS))) ==NULL) 
		ath_error("[SMR_Rad_init]:Failed to allocate GZ[0]C\n");
  	

	return;
}



void SMR_Rad_destruct()
{

	int i;
			
	
	if(start_addrP != NULL){
		free(start_addrP);
		start_addrP = NULL;
	}

/*==========================*/
	/* Free the MPI buffer */

	if(send_bufRC != NULL){
		free_2d_array(send_bufRC);
		send_bufRC = NULL;
	}

#ifdef MPI_PARALLEL

	if(send_rq != NULL){
		free_2d_array(send_rq);
		send_rq = NULL;
	}

	if(recv_rq != NULL){
		free_2d_array(recv_rq);
		recv_rq = NULL;
	}

	if(recv_bufRC != NULL){
		free_2d_array(recv_bufRC);
		recv_bufRC = NULL;
	}

#endif

	if(send_bufP != NULL){
		free_2d_array(send_bufP);
		send_bufP = NULL;
	}

	if(recv_bufP != NULL){
		free_2d_array(recv_bufP);
		recv_bufP = NULL;
	}


	/* Free temporary array for prolongation ghost zones */
	for(i=0; i<3; i++){
		if(GZ[i] != NULL){
			free_3d_array(GZ[i]);
			GZ[i] = NULL;
		}

	}

	if(Pro_buf != NULL){
		free_3d_array(Pro_buf);
		Pro_buf = NULL;
	}
	return;


}












#endif /* only for the radiation integrator */
#endif /* End matrix multigrid */
#endif /* STATIC_MESH_REFINEMENT */
