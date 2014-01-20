#include "../copyright.h"
/*==============================================================================
 * FILE: utils_fullrad.c
 *
 * PURPOSE: contains misc. functions require for computation of rad. transfer
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   get_weights_linear()     - 
 *   get_weights_parabolic()  -
 *============================================================================*/

#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"


#ifdef FULL_RADIATION_TRANSFER

/* This function is called to update the momentums in the ghost zones only */


void UpdateRT(DomainS *pD){
	
	RadGridS *pRG=(pD->RadGrid);

	
	int i;
	int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
	int jl = pRG->js, ju = pRG->je;
	int kl = pRG->ks, ku = pRG->ke;
	int koff = 0, joff = 0, ioff = 0;
	int nDim; 
	
	
	nDim = 1;
	for (i=1; i<3; i++) if (pRG->Nx[i]>1) nDim++;
	
	ioff = nghost - Radghost;
	
	/* Including the ghost zones */
	
	if(nDim > 1){
		jl -= Radghost;
		ju += Radghost;
		
		joff = nghost - Radghost;
	}
	
	if(nDim > 2){
		kl -= Radghost;
		ku += Radghost;
		
		koff = nghost - Radghost;
	}
	
	

		CalMoment(il, iu, jl, ju, kl, ku, pRG);
		
    return;
}


/* Calculate the radiation moments for the specific index range:
 * il to iu, jl to ju, kl to ku 
 */


void CalMoment(const int il, const int iu, const int jl, const int ju, const int kl, const int ku, RadGridS *pRG)
{
	int i, j, k, Mi, m;
	Real wimu;
	Real mu[3]; /* cosins with respect to three axis */
	Real mu2[6]; 	/* products of two angles, used for radiation pressure */
	int nelements = pRG->noct * pRG->nang;
	int nf = pRG->nf;
	int ifr;
	
	/* First, initialize to be zero */
	
	for(k=kl; k<=ku; k++){
		for(j=jl; j<=ju; j++){	
			for(i=il; i<=iu; i++){	
				for(ifr=0; ifr<nf; ifr++){
			
					pRG->R[k][j][i][ifr].J = 0.0;
					for(m=0; m<3; m++)
						pRG->R[k][j][i][ifr].H[m] = 0.0;
				
					for(m=0; m<6; m++)
						pRG->R[k][j][i][ifr].K[m] = 0.0;
						
					for(Mi=0; Mi<nelements; Mi++){
						/* First, check specific intensity is not negative */	  
						if(pRG->imu[k][j][i][ifr][Mi] < TINY_NUMBER)
						 	pRG->imu[k][j][i][ifr][Mi] = TINY_NUMBER;
						 
						/* sum rays along different directions */
						wimu = pRG->imu[k][j][i][ifr][Mi] * pRG->wmu[k][j][i][Mi];
						
						/* For cylindrical coordinate, we nned to convert the angle */
#ifdef CYLINDRICAL
						for(m=0; m<3; m++)
							mu[m] = pRG->Rphimu[k][j][i][Mi][m];
#else
						for(m=0; m<3; m++)
							mu[m] = pRG->mu[k][j][i][Mi][m];
#endif
						
						/* for energy density */
						pRG->R[k][j][i][ifr].J += wimu;
						
						/* for flux */		
						for(m=0; m<3; m++)
							pRG->R[k][j][i][ifr].H[m] += mu[m] * wimu;
						
						/* for radiation pressure */
						mu2[0] = mu[0] * mu[0];
						mu2[1] = mu[0] * mu[1];
						mu2[2] = mu[1] * mu[1];
						mu2[3] = mu[0] * mu[2];
						mu2[4] = mu[1] * mu[2];
						mu2[5] = mu[2] * mu[2];
						
						for(m=0; m<6; m++)
							pRG->R[k][j][i][ifr].K[m] += mu2[m] * wimu;
										
					}/* end Mi */	
				}/* end ifr */		
			}/* end i*/
		}/* end j */
	}/* end k */
	

}




/* use LU decomposition to solve a special matrix, 
 * resulting from the coupling between different angles 
 * due to scattering opacity 
 * The matrix form is 
 * Ma1+Mb1  Mb2     Mb3  ... ... Mbn   *
 * Mb1     Ma2+Mb2  Mb3  ... ... Mbn   *
 * Mb1     Mb2    Ma3+Mb3        Mbn   *
 * ...   ...     ...		       *
 * Mb1     Mb2     Mb3 ...  ... Man+Mbn*
 
 * We can do the LU decomposition by hand 
 * for this special matrix 
 */

void SpecialMatrix(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N)
{
    Real swap, Mas, Mbs;
    Real temp;
    int i, Nzl, j;



   /* first, search Mb, find the first non-zero element */
   i = 0;
   while(((fabs(Mb[i]) < TINY_NUMBER))&&(i < N))
	i++;
  
   Nzl = i;

  if(Nzl == N){
    /* All the elements of Mb are zero *
     * The angles are decoupled */
     for(i=0; i<N; i++)
	RHS[i] /= Ma[i];
  }
  else{
	
   	/* switch the RHS between the last line and the first none zero line */
	/* The order of the solution does not change */
  	swap = RHS[Nzl];
  	RHS[Nzl] = RHS[N-1];
  	RHS[N-1] = swap;

	Mas = Ma[Nzl];
	Mbs = Mb[Nzl];
   /*-----------------------------------*/
   /* ------------------------------------ */
   /* The lower triangle matrix is *
    * 1   0   0  ... ... 0 *
    * 1   1   0 ... ...  0 *
    * 1   0   1 ... ...  0 *
    * l1 l2   l3 ... ... ln *
    */

   	if(Nzl < N-1){
		/* If Nzl = N-1, the L matrix is trival */
		/* No need to calculate in this case */
		lN[Nzl] = (Mas + Mbs)/Mbs;
	
	
		for(i=Nzl+1; i<N-1; i++){
			lN[i] = -Mas * Mb[i]/(Mbs * Ma[i]);
		}
			lN[N-1] = 1.0;

		/* Now multiple the Inversion of L and RHS */
		/* The 0 --- Nzl-1 lines are unchanged */
		temp = -lN[Nzl];
		tempRHS[N-1] = 0.0;
		for(i=Nzl+1; i< N-1; i++){
			tempRHS[i] = RHS[i] - RHS[Nzl];
			temp += lN[i];
			tempRHS[N-1] += (-lN[i] * RHS[i]);
		}
			tempRHS[N-1] += (temp * RHS[Nzl] + RHS[N-1]);
 
		/* Now copy tempRHS to RHS */
		for(i=Nzl+1; i<N; i++)
			RHS[i] = tempRHS[i];

	}/* end if Nzl < N-1 */

	/*-------------------------------------------*/
	/*--------------------------------------------*/
	/* Now handle the upper triangle matrix */
	/* The first element */

	/* the last element */

	if(Nzl < N-1){
		UN[N-1] = Ma[N-1] + Mb[N-1];
		for(i=N-2; i>=Nzl; i--){
			UN[N-1] *= Ma[i];
			temp = Mb[i];
			for(j=i+1; j<N; j++)
				temp *= Ma[j];
			UN[N-1] += temp;
		}
		UN[N-1] *= (-1.0);	

		temp = Mb[Nzl];
		for(i=Nzl+1; i<N-1; i++)
			temp *= Ma[i];

		UN[N-1] /= temp;

	}/* end the line N-1 */

	/* Now we have the Upper matrix, solve the equation */

	RHS[N-1] /= UN[N-1];
 
	for(i=N-2; i>Nzl; i--){
		RHS[i] = (RHS[i] + RHS[N-1] * Ma[N-1]) / Ma[i];
	}

	/* The line Nzl */
	for(i=Nzl+1; i<N; i++)
		RHS[Nzl] = RHS[Nzl] - Mb[i] * RHS[i];

		RHS[Nzl] = RHS[Nzl] - Ma[N-1] * RHS[N-1];
		RHS[Nzl] /= Mb[Nzl];

	/* The lines between 0 and Nzl-1 */
	for(j=0; j<Nzl; j++){
		for(i=j+1; i<N; i++){
			RHS[j] = RHS[j] - Mb[i] * RHS[i];
		}

		RHS[j] /= Ma[j];
	}


  }/* end Nzl < N */



   return;

}


/* use LU decomposition to solve a special matrix, 
 * resulting from the coupling between different angles 
 * due to scattering opacity 
 * The matrix form is 
 * Ma1+Mb1  Mb1     Mb1  ... ... Mb1   *
 * Mb2     Ma2+Mb2  Mb2  ... ... Mb2   *
 * Mb3     Mb3    Ma3+Mb3        Mb3   *
 * ...   ...     ...		       *
 * Mbn     Mbn     Mbn ...  ... Man+Mbn*
 
 * We can do the LU decomposition by hand 
 * for this special matrix
 * This matrix is the transpose of previous matrix
 * But we do the decomposition slightly different  
 */

void SpecialMatrix2(Real *Ma, Real *Mb, Real *RHS, Real *lN, Real *tempRHS, Real *UN, const int N)
{
    Real swap;
    Real temp;
    int i, Nzl, j;



   /* first, search Mb, find the first non-zero element */
   i = N-1;
   while(((fabs(Mb[i]) < TINY_NUMBER))&&(i >= 0))
	i--;
  
   Nzl = i;

  if(Nzl < 0){
    /* All the elements of Mb are zero *
     * The angles are decoupled */
     for(i=0; i<N; i++)
	RHS[i] /= Ma[i];
  }
  else{
	
   	/* switch the RHS between the first line and the first none zero line */
	/* The order of the solution does not change */
  	swap = RHS[Nzl];
  	RHS[Nzl] = RHS[0];
  	RHS[0] = swap;

   /*-----------------------------------*/
   /* ------------------------------------ */
   /* The lower triangle matrix is *
    * 1   0   0  ... ... 0 *
    *     1   0 ... ...  0 *
    * 1      1 ... ...  0 *
    * l1 l2   l3 ... ... ln *
    */

   	if(Nzl > 0){
		/* If Nzl = 0, the L matrix is trival */
		/* No need to calculate in this case */
		lN[0] = (Ma[0]+Mb[0])/Mb[Nzl];
	
	
		for(i=1; i<Nzl; i++){
			lN[i] = -Ma[0]/Ma[i];
		}
			lN[Nzl] = 1.0;

		/* Now multiple the Inversion of L and RHS */
		/* The 0 --- Nzl-1 lines are unchanged */
		temp = -lN[0];
		tempRHS[Nzl] = RHS[Nzl];
		for(i=1; i< Nzl; i++){
			tempRHS[i] = RHS[i] - RHS[0] * Mb[i]/Mb[Nzl];
			temp += lN[i] * Mb[i]/Mb[Nzl];
			tempRHS[Nzl] += (-lN[i] * RHS[i]);
		}
			tempRHS[Nzl] += (temp * RHS[0]);
 
		/* Now copy tempRHS to RHS */
		for(i=1; i<=Nzl; i++)
			RHS[i] = tempRHS[i];

	}/* end if Nzl > 0 */

	/*-------------------------------------------*/
	/*--------------------------------------------*/
	/* Now handle the upper triangle matrix */
		
	/* The first element */
	UN[0] = Mb[Nzl];

	for(i=1; i<Nzl; i++)
		UN[i] = -Ma[Nzl] * Mb[i]/Mb[Nzl];

	/* now the element Nzl */
	if(Nzl > 0){
		UN[Nzl] = Ma[Nzl] + Mb[Nzl];
		for(i=Nzl-1; i>=0; i--){
			UN[Nzl] *= Ma[i];
			temp = Mb[i];
                        for(j=Nzl; j>i; j--)
                                temp *= Ma[j];

                        UN[Nzl] += temp;
		}

		UN[Nzl] *= -1.;
		temp = Mb[Nzl];
		for(i=1; i<Nzl; i++)
			temp *= Ma[i];
		UN[Nzl] /= temp;
	}

	/* Now we have the upper triangle matrix */
	/* The lines between Nzl +1 to N-1 */
	for(i=N-1; i>Nzl; i--)
		RHS[i] /= Ma[i];

	/* The line Nzl */
	if(Nzl > 0){
		for(i=Nzl+1; i<N; i++)
			RHS[Nzl] += Ma[0] * RHS[i];

		RHS[Nzl] /= UN[Nzl];
	}

	/* The lines between 1 and Nzl-1 */
	for(i=1; i<Nzl; i++)
		RHS[i] = (RHS[i] - UN[i] * RHS[Nzl])/Ma[i];

	/* The first element */
	for(i=1; i<N; i++)
		RHS[0] = RHS[0] - Mb[Nzl] * RHS[i];

		RHS[0] = RHS[0] - Ma[Nzl] * RHS[Nzl];
		RHS[0] /= Mb[Nzl];



  }/* end Nzl < N */



   return;


}

/* use LU decomposition to solve a special matrix, 
 * resulting from the resulting from full 
 * radiation transfer equation with velocity dependent terms  
 * The matrix form is 
 * Ma1+Mc1+Mb1  Mc1+Mb2     Mc1+Mb3  ... ... Mc1+Mbn   *
 * Mc2+Mb1     Ma2+Mc2+Mb2  Mc2+Mb3  ... ... Mc2+Mbn   *
 * Mc3+Mb1     Mc3+Mb2    Ma3+Mc3+Mb3        Mc3+Mbn   *
 * ...   ...     ...                   *
 * Mcn+Mb1     Mcn+Mb2     Mcn+Mb3 ...  ... Man+Mcn+Mbn*
 
 * We can do the LU decomposition by hand 
 * for this special matrix
 * To solve this matrix, we first subtract the last line from each line *
 * Then we add the first line to the last line *
 * Then we redefine Mc = Mc-Mb and Md = Mc - Mcn
 */

/* We only need three lines for lN: lN[0], lN[Nl], lN[N-1] */
/* We also only need three lines for UN */
/* lN1[] is lN[0][] *
 * lN2[] is lN[Nl][] *
 * lN3[] is lN[N-1][] */

/* UN1[] is UN[1][]    *
 * UN2[] is UN[Nl][]  *
 * UN3[] is UN[N-1][]
 */

void SpecialMatrix3(const int N, Real **Ma, Real *Md, Real *RHS,  Real *lN1, Real *lN2, Real *lN3,  Real *UN1, Real *UN2, Real *UN3)
{

	Real swap;
	Real tempMc,temp,tempup, tempdown, MCB0, MCB1, MCB2, BInv, CInv, norm;
	int i, Nzl, j; 
	
	/* First, subtract the last line from lines 1 to N-1 */
	/* RHS is also changed during this process but Xi unchanged */
	tempMc = Ma[0][2];
	
	
	/* Redfined Mb */
	Ma[0][1] += Ma[0][0] + tempMc;
	for(i=1; i<N; i++)
		Ma[i][1] += tempMc; 
		
		
	/* Redefine Mc */	
		
	tempMc = Ma[N-1][2];
    temp = RHS[0];

	for(i=0; i<N-1; i++){
		RHS[i] -= RHS[N-1];
		/* Redefine Mc */
		Ma[i][2] -= tempMc;
	}
	RHS[N-1] = temp;
	

	for(i=0; i<N-1; i++)
		Md[i] = Ma[i][2] - Ma[N-1][0];

	/* Mc[N-1] is no longer used anymore */
	/* Ma[i] and Mb[0] are always no-zero */

	/* Find the first non-zero Mc */
	i=N-2;
    
    if(fabs(tempMc) > 1.e-15)
        norm = tempMc;
    else
        norm = 1.0;

    
	while((fabs(Ma[i][2]/norm) < 1.e-13)){
        	i--;
        if( i<0 ) break;
    }

   	Nzl = i;
	
       
	if(Nzl < 0){
        /* For the isotripic case, do something special to conserve energy */
        
     		/* First, Multiple L^-1 and RHS */
		for(i=0; i<N-1; i++){
			RHS[N-1] -= Ma[i][1] * RHS[i]/Ma[i][0];
		}
		/* Then, Multiple U^-1 and RHS */
		/* First RHS[N-1] */
		temp = Ma[N-1][1];
		for(i=0; i<N-1; i++)
			temp += Ma[i][1] * Ma[N-1][0]/Ma[i][0];

		RHS[N-1] /= temp;

		/* Now solve RHS[0] to RHS[n-2]	*/
		for(i=0; i<N-1; i++)
			RHS[i] = (RHS[i] + Ma[N-1][0] * RHS[N-1])/Ma[i][0];		

	}/* End if Nzl < 0 */
	else if(Nzl == 0){

        	/* Exchange the first and last line */

		swap = RHS[0];
		RHS[0] = RHS[N-1];
		RHS[N-1] = swap;

		/* calculate L^-1 RHS */

	/*	lN[0][0] = (Ma[0] + Mc[0])/Mb[0];
	*/
		MCB0 = (Ma[0][0] + Ma[0][2])/Ma[0][1];
		lN1[0] = MCB0;
/*
		for(i=1; i<N-1; i++)
			lN[0][i] = (Mc[0] - Mb[i] * (Ma[0] + Mc[0])/Mb[0])/Ma[i];
*/
		for(i=1; i<N-1; i++)
			lN1[i] = (Ma[0][2] - Ma[i][1] * MCB0)/Ma[i][0];

		for(i=0; i<N-1; i++)
			RHS[N-1] -= lN1[i] * RHS[i];

		/* calculate U^-1 RHS */
		temp = -Ma[N-1][0] + Ma[0][2] - Ma[N-1][1] * MCB0;
		for(i=1; i<N-1; i++)
			temp += Ma[N-1][0] * (Ma[0][2] - Ma[i][1] * MCB0)/Ma[i][0];

		RHS[N-1] /= temp;

		for(i=1; i<N-1; i++){
			RHS[i] = (RHS[i] + Ma[N-1][0] * RHS[N-1]) / Ma[i][0];
		}

		for(i=1; i<N; i++)
			RHS[0] -= Ma[i][1] * RHS[i];

		RHS[0] /= Ma[0][1];


	} /* end Nzl ==0 */	
	else{
		/* For Nzl +1 to N-2 , The following equation is true */
		/* Ma[i] * xi + Md[i] xn-1  = RHS[i] * 
		 * Therefore  xi =  RHS[i]/Ma[i] - Md[i]/Ma[i] xn-1 *
		 * The elements added to line j  are Mc[j] * RHS[i] /Ma[i] - Mc[j] * Md[i]/Ma[i] xn-1 */
		
		/* Replace all xi with xn-1 for lines between Nzl+1 to N-2 */
		for(i=Nzl+1; i<N-1; i++){
			for(j=0; j<=Nzl; j++){
				/* Redefine RHS */
				RHS[j] -= (Ma[j][2] * RHS[i] / Ma[i][0]);
				/* Redefine Md[j] */
				Md[j] -= (Ma[j][2] * Md[i] / Ma[i][0]);
			}
			/* Also replace the line N-1 */
			RHS[N-1] -= (Ma[i][1] * RHS[i] / Ma[i][0]);
			Ma[N-1][1] -= (Ma[i][1] * Md[i]/Ma[i][0]);
		}

		/* Now we only need to solve lines between 0 to Nzl and N-1 lines */
		/* First, swap the first and last line */
		swap = RHS[0];
                RHS[0] = RHS[N-1];
                RHS[N-1] = swap;
		
		/* If Nzl > 1, also need to swap line 1 and Nzl */
		if(Nzl > 1){
			swap = RHS[1];
	                RHS[1] = RHS[Nzl];
        	        RHS[Nzl] = swap;
		}
	
		/*-----------------------------------------*/

		/* Now calculate the L matrix */
		/* For line 1 to Nzl-1, the first element is Mc[Nzl]/Mb[0] *
		 * The second element (starting from line 2 ) is Mc[Nzl-i+1]/Mc[Nzl]
		*/
		/* Now the line Nzl */
		MCB0 = (Ma[0][0] + Ma[0][2])/Ma[0][1];
		MCB1 = Ma[1][0]/(Ma[0][1] - Ma[1][1]);
		MCB2 = Ma[0][0]/(Ma[0][1] - Ma[1][1]);
		
	/*	lN[Nzl][0] = Mc[1]/Mb[0];
		lN[N-1][0] = (Ma[0] + Mc[0])/Mb[0];
	*/
		lN2[0] = Ma[1][2]/Ma[0][1];
		lN3[0] = MCB0;
		/* The second element */
		if(Nzl > 1){
			lN2[1] = (Ma[1][0] * Ma[0][1] + (Ma[0][1] - Ma[1][1]) * Ma[1][2])/((Ma[0][1] - Ma[1][1]) * Ma[Nzl][2]);
			lN3[1] = -(Ma[0][0] * Ma[1][1] + Ma[0][2] * (Ma[1][1] - Ma[0][1]))/((Ma[0][1] - Ma[1][1]) * Ma[Nzl][2]);
		}

		for(i=2; i<Nzl; i++){
			lN2[i] = MCB1 * (Ma[i][1] - Ma[0][1])/Ma[i][0];
			lN3[i] = MCB2 * (Ma[1][1] - Ma[i][1])/Ma[i][0];
		}
			


		/* Now the element [N-1][Nzl] */
		/* To avoid many Ma[i] especially when Ma[i] is huge */
		/* Divide a common factor for up and down */
		tempup = Ma[0][2] * (Ma[0][1] - Ma[1][1])/(Ma[0][0]*Ma[1][0]);
		temp = Ma[1][1]/Ma[1][0];
		tempup -= temp;

		for(j=2; j<=Nzl; j++){
			temp = Ma[j][2] * (Ma[1][1] - Ma[j][1])/(Ma[1][0]*Ma[j][0]);
			tempup -= temp;
		}

		tempdown = Ma[1][2] * (Ma[0][1] - Ma[1][1])/(Ma[0][0]*Ma[1][0]);
		temp = Ma[0][1]/Ma[0][0];
		tempdown += temp;

		for(j=2; j<=Nzl; j++){
			temp = Ma[j][2] * (Ma[0][1] - Ma[j][1])/(Ma[0][0]*Ma[j][0]);
			tempdown += temp;
		}


		lN3[Nzl] = tempup/tempdown;
		
		/*-----------------------------------------*/

		/* Now calculate L^-1 RHS */
		/* The first element is unchanged */
		/* The second element */
		BInv = 1.0/Ma[0][1];
		CInv = 1.0/Ma[Nzl][2];
		
		RHS[1] -= (Ma[Nzl][2] * RHS[0] * BInv);
		for(j=2; j<Nzl; j++){
			RHS[j] -= (Ma[j][2] * RHS[0] * BInv);
			RHS[j] -= (Ma[j][2] * RHS[1] * CInv);
		}
		/* The line Nzl */
		if(Nzl > 1){
			for(i=0; i<Nzl; i++)
				RHS[Nzl] -= (lN2[i] * RHS[i]);
		}

		/* The line N-1 */
		for(i=0; i<=Nzl; i++)
			RHS[N-1] -= (lN3[i] * RHS[i]);

		/*-----------------------------------------------*/
		/* Now construct the U matrix */
		/* UN1[] is UN[1][]    *
		 * UN2[] is UN[Nzl][]  *
		 * UN3[] is UN[N-1][]
		 */
		

		/* For line 1 */

	
		for(i=1; i<=Nzl; i++)
			UN1[i] = (Ma[0][1] - Ma[i][1]) * Ma[Nzl][2] * BInv;

		UN1[Nzl] += Ma[Nzl][0];
		
		UN1[N-1] = Md[Nzl] - Ma[N-1][1] * Ma[Nzl][2] * BInv;

		/* For the lines 2 to Nzl-1, the elements are a[i], ... -a[Nzl]*c[i]/c[Nzl];

		* now the line Nzl, there are only two elements in lines Nzl */
		if(Nzl > 1){
			UN2[Nzl] = -Ma[0][1];

			for(j=1; j<=Nzl; j++){
				temp = (Ma[0][1] - Ma[j][1]) * Ma[j][2] / Ma[j][0];
				UN2[Nzl] -= temp;
			}

			temp = (Ma[0][1] - Ma[1][1]) * Ma[Nzl][2] / (Ma[1][0] * Ma[Nzl][0]);
			UN2[Nzl] /= temp;
		

		/* Now the element [Nzl][N-1] */
			UN2[N-1] = Ma[N-1][1] * Ma[Nzl][2] - Ma[0][1] * Md[Nzl];
		
			for(j=1; j<Nzl; j++){
				temp = (Ma[Nzl][2] * Md[j] - Ma[j][2] * Md[Nzl]) * (Ma[0][1] - Ma[j][1])/Ma[j][0];
				UN2[N-1] += temp;
			}

			temp = (Ma[0][1] - Ma[1][1]) * Ma[Nzl][2]/Ma[1][0];

			UN2[N-1] /= temp;
		}/* End Nzl > 1 */
		else {
			/* In this case, array UN2 is the same as array UN1 */
			UN2[Nzl] = UN1[Nzl];
			UN2[N-1] = UN1[N-1];
		}


		/* now the last element UN[N-1][N-1] */
		/* sum LN[N-1][i] * UN[i][N-1] = Mb */
		UN3[N-1] = Md[0] - lN3[0] * Ma[N-1][1] - lN3[1] * UN1[N-1];
		for(i=2; i<Nzl; i++)
			UN3[N-1] -= (lN3[i] * (Md[i] - Ma[i][2] * Md[Nzl] * CInv));

		if(Nzl > 1){
			UN3[N-1] -= (lN3[Nzl] * UN2[N-1]);
		}
		
		/* Now calculate U^-1 RHS */
		RHS[N-1] /= UN3[N-1];
	
		RHS[Nzl] -= (UN2[N-1] * RHS[N-1]);
		RHS[Nzl] /= UN2[Nzl]; 

		for(i=Nzl-1; i>1; i--){
			RHS[i] -= (Md[i] - Ma[i][2] * Md[Nzl] * CInv) * RHS[N-1]; 
			RHS[i] += (Ma[Nzl][0] * Ma[i][2] * CInv) * RHS[Nzl];
			RHS[i] /= Ma[i][0];
		}
		
		if(Nzl > 1){	
			for(i=2; i<=Nzl; i++)
				RHS[1] -= (UN1[i] * RHS[i]);

			RHS[1] -= (UN1[N-1] * RHS[N-1]);
			RHS[1] /= UN1[1];
		}
			
		for(i=1; i<=Nzl; i++)
			RHS[0] -= Ma[i][1] * RHS[i];

		RHS[0] -= (Ma[N-1][1] * RHS[N-1]);
		RHS[0] /= Ma[0][1];


			/* now calculate the solution between Nzl+1 to N-2 */
			/* xi =  RHS[i]/Ma[i] - Md[i]/Ma[i] xn-1 */
		for(i=Nzl+1; i<N-1; i++)
			RHS[i] = (RHS[i] - Md[i] * RHS[N-1])/Ma[i][0];			

		
	}/* end if Nzl > 0 */

}



/* use LU decomposition to invert a special matrix, 
 * resulting from the all the terms, especially 
 * the velocity dependent terms related with absorption opacity
 * The matrix form is 
 * Ma1+Mb1  Mb2     Mb3  ... ... Mbn   Md1		*
 * Mb1     Ma2+Mb2  Mb3  ... ... Mbn   Md2		*
 * Mb1     Mb2    Ma3+Mb3        Mbn   Md3		*
 * ...   ...     ...							*
 * Mb1     Mb2     Mb3 ...  ... Man+Mbn Mdn		*
 * Mc1		Mc2		Mc3	... ...  Mcn	Mdn+1	*
 
 * We can do the LU decomposition by hand 
 * for this special matrix
 * To solve this matrix, we first subtract the n line from each line *
 * Then we redefine Md = Md-Mdn
 */




void AbsorptionMatrix(const int N, Real **Ma, Real *Md, Real *RHS)
{
	
	/* N always larger than 2 */
	/* The total line s of the matrix is N+1 */
	/* Lines 1 to N are for weight * I */
	/* Lines N+1 is for gas temperature */
	
	Real CoefN[2], CoefN1[2], tempRHS[2], temp;
	Real MaInv;
	int i; 
	
	/* Now the matrix becomes the form */
	/* Ma1 0    0   0  ... ...	 -Man Md1 *
	 * 0   Ma2  0   0  ... ...   -Man Md2 *
	 * 0   0    Ma3 0  ... ...   -Man Md3 *
	 * ... ...  ...	...
	 * ... ... ... ...
	 * Mb1 Mb2  Mb3 Mb4...  ...  Man+Mbn Mdn *
	 * Mc1 Mc2 Mc3  Mc4 ... ...  Mcn  Mdn+1 *
	 */
		
	/* All the variables 1 to N-1 can be expressed in terms of variable N and N+1 */
	CoefN[0] = Ma[N-1][0] + Ma[N-1][1];
	CoefN[1] = Md[N-1];
	
	CoefN1[0] = Ma[N-1][2];
	CoefN1[1] = Md[N];
	
	tempRHS[0] = RHS[N-1];
	tempRHS[1] = RHS[N];
	
	for(i=0; i<N-1; i++){
		MaInv = 1.0/Ma[i][0];
		/* line N */
		tempRHS[0] += (-Ma[i][1] * MaInv * RHS[i]);
		CoefN[0] += (Ma[i][1] * MaInv * Ma[N-1][0]);
		CoefN[1] += (-Ma[i][1] * MaInv * Md[i]);
		
		/* line N+1 */
		tempRHS[1] += (-Ma[i][2] * MaInv * RHS[i]);
		CoefN1[0] += (Ma[i][2] * MaInv * Ma[N-1][0]);
		CoefN1[1] += (-Ma[i][2] * MaInv * Md[i]);		
		
	}
	
	/* Now solve the 2*2 Matrix */
	/* CoefN[0] is always non zero */
	/* CoefN[0] * xN + CoefN[1] * xN1 = tempRHS[0]
	 * CoefN1[0] * xN + CoefN1[1] * xN1 = tempRHS[1]
	 */
	
	temp = CoefN[1] * CoefN1[0]/CoefN[0] - CoefN1[1];
	RHS[N] = (tempRHS[0] * CoefN1[0]/CoefN[0] - tempRHS[1])/temp;
	RHS[N-1] = (tempRHS[0] - CoefN[1] * RHS[N])/CoefN[0];
	
	
	temp = -Ma[N-1][0] * RHS[N-1];
	
	/* Now set the solution */
	for(i=0; i<N-1; i++){
		RHS[i] = (RHS[i] - (temp + Md[i] * RHS[N]))/Ma[i][0];		
	}
	
	/* Now the new solutions are stored in RHS[] */
}




/* This is the matrix for absorption opacity when 
 * gas temperature is taken fixed
 * This is used when the non-linear matrix cannot converge
 * The matrix form is
 * Ma1+Mb1  Mb2     Mb3  ... ... Mbn   		*
 * Mb1     Ma2+Mb2  Mb3  ... ... Mbn   		*
 * Mb1     Mb2    Ma3+Mb3        Mbn   		*
 * ...   ...     ...							*
 * Mb1     Mb2     Mb3 ...  ... Man+Mbn		*
 
 * We can do the LU decomposition by hand
 * for this special matrix
 * To solve this matrix, we first subtract the n line from each line *
 */




void LinearAbsorptionMatrix(const int N, Real **Ma, Real *RHS)
{
	
	/* N always larger than 2 */
	/* The total line s of the matrix is N */
	/* Lines 1 to N are for  I */
	
	Real CoefN, tempRHS, temp;
	Real MaInv;
	int i;
	
	/* Now the matrix becomes the form */
	/* Ma1 0    0   0  ... ...	 -Man*
	 * 0   Ma2  0   0  ... ...   -Man *
	 * 0   0    Ma3 0  ... ...   -Man *
	 * ... ...  ...	...
	 * ... ... ... ...
	 * Mb1 Mb2  Mb3 Mb4...  ...  Man+Mbn *
	 */
    
	/* All the variables 1 to N-1 can be expressed in terms of variable N  */
	CoefN = Ma[N-1][0] + Ma[N-1][1];
	
	tempRHS = RHS[N-1];

	
	for(i=0; i<N-1; i++){
		MaInv = 1.0/Ma[i][0];
		/* line N */
		tempRHS += (-Ma[i][1] * MaInv * RHS[i]);
		CoefN += (Ma[i][1] * MaInv * Ma[N-1][0]);
	}
    
    /* Now get solution N */
    RHS[N-1] = tempRHS / CoefN;
	
	
	temp = -Ma[N-1][0] * RHS[N-1];
	
	/* Now set the solution */
	for(i=0; i<N-1; i++){
		RHS[i] = (RHS[i] - temp)/Ma[i][0];
	}
	
	/* Now the new solutions are stored in RHS[] */
}


void ReduceVelocity(const Real sigma, const Real ds, Real *alpha)
{

	Real tau, y, y2, y3, y4, expression;
	tau = Taufactor * ds * sigma;
	tau = tau * tau;

	/* As this function is called many times, we approximate the function */
	/* (1-exp(-tau))/tau with different functions at different range to */
	/* speed up the calculation */
/*

	if(tau > 0.001)
		*alpha = sqrt((1.0 - exp(- tau)) / tau);
	else
		*alpha = sqrt(1.0 - 0.5 * tau);
*/
	if(tau < 0.574374765791614){
		y2 = tau * tau;
		y3 = y2 * tau;
		y4 = y3 * tau;
		expression = 1.0 - 0.5* tau + 0.16666666666666666 * y2 - 0.041666666666666664 * y3 + 0.008333333333333333 * y4;
	}else if(tau < 1.5){
		y = 1.0 - tau;
		y2 = y * y;
		y3 = y2 * y;
		y4 = y3 * y;
		expression = (0.6321205588285577 + 0.36787944117144233 * y - 0.18393972058572117 * y2 + 0.061313240195240384 * y3 - 0.015328310048810096 * y4)/tau;
	}else if(tau < 2.5){
		y = 1.0/tau - 0.5;
		y2 = y * y;
		y3 = y2 * y;
		y4 = y3 * y;
		expression = 0.43233235838169365 + 0.5939941502901619 * y - 0.5413411329464508 * y2 + 0.36089408863096717 * y3;

	}else if(tau < 5.268707747422607){
		y = 1.0/tau - 0.2857142857142857;
		y2 = y * y;
		y3 = y2 * y;
		y4 = y3 * y;
		expression = 0.27708646187933755 + 0.8641117745995668 * y -  0.6473564071159528 * y2 - 0.3776245708176391 * y3 + 2.4781612459907563 * y4;

	}else if(tau < 11.411075605604893){
		y=1.0/tau - 0.125;
		y2 = y * y;
		y3 = y2 * y;
		y4 = y3 * y;
		expression = 0.12495806717151219 + 0.9969808363488774 * y - 0.08587843274304303 * y2 - 1.1450457699072405 * y3 - 5.496219695554754 * y4;

	}else{
		expression = 1.0/tau;
	}
	
	*alpha = sqrt(expression);

	return;

}

/* Convert the cosphi values from the global cartesian coordinate to the local cylindrical coordinate */
/* cosphi0 = miux0, sinphi0 = miuy0 
 * miux = cos(phi0 - x2), miuy = sin(phi0-x2)
 */
void convert_angle(const Real x2, const Real miux0, const Real miuy0, Real *miux, Real *miuy)
{
	Real cosphi, sinphi;
	cosphi = cos(x2);
	sinphi = sin(x2);
	
	*miux = cosphi * miux0 + sinphi * miuy0;
	*miuy = cosphi * miuy0 - sinphi * miux0;


}


#endif /* RADIATION_TRANSFER */


