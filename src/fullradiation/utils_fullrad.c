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

void UpdateRT(DomainS *pD){

  RadGridS *pRG=(pD->RadGrid);
  GridS *pG = (pD->Grid);
  int i,j,k ,ifr, l, n, m;
#ifdef CYLINDRICAL
  Real x1, x2, x3;
#endif

  int il = pRG->is-Radghost, iu = pRG->ie+Radghost;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int kg, jg, ig;
  int koff = 0, joff = 0, ioff = 0;
  int nDim; 
  Real vel[3], Fr0[3], VKr[3], VEr[3];;

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

  Real wimu;
  Real mu[3]; /* cosins with respect to three axis */
  Real mu2[6]; 	/* products of two angles, used for radiation pressure */


for(ifr=0; ifr<pRG->nf; ifr++){
	/* First, initialize to be zero */
  for(k=kl; k<=ku; k++){
     for(j=jl; j<=ju; j++){	
	for(i=il; i<=iu; i++){	
		pRG->R[ifr][k][j][i].J = 0.0;
		for(m=0; m<3; m++)
			pRG->R[ifr][k][j][i].H[m] = 0.0;

		for(m=0; m<6; m++)
			pRG->R[ifr][k][j][i].K[m] = 0.0;
	}
     }
  }

  /* now recalculate the momentums */


  for(l=0; l<pRG->noct; l++){
  	for(n=0; n<pRG->nang; n++){
  		for(k=kl; k<=ku; k++){
     			for(j=jl; j<=ju; j++){	
				for(i=il; i<=iu; i++){	  
		

					/* sum rays along different directions */
					wimu = pRG->imu[ifr][l][n][k][j][i] * pRG->wmu[n][k][j][i];
					for(m=0; m<3; m++)
						mu[m] = pRG->mu[l][n][k][j][i][m];

					/* For cylindrical coordinate, we nned to convert the angle */
#ifdef CYLINDRICAL
					cc_pos(pG,i+ioff,j+joff,k+koff,&x1,&x2,&x3);
					convert_angle(x2,mu[0],mu[1],&mu[0],&mu[1]);
#endif

					/* for energy density */
					pRG->R[ifr][k][j][i].J += wimu;
	
					/* for flux */		
					for(m=0; m<3; m++)
						pRG->R[ifr][k][j][i].H[m] += mu[m] * wimu;

					/* for radiation pressure */
					mu2[0] = mu[0] * mu[0];
					mu2[1] = mu[0] * mu[1];
					mu2[2] = mu[1] * mu[1];
					mu2[3] = mu[0] * mu[2];
					mu2[4] = mu[1] * mu[2];
					mu2[5] = mu[2] * mu[2];

					for(m=0; m<6; m++)
						pRG->R[ifr][k][j][i].K[m] += mu2[m] * wimu;
				}/* end i  */
			}/* end j */
		}/* end k */
     	}/* end n */
  }/* end l */
} /* end frequency */

	/*-----------------------------------------------------*/
	/* Now calculate the frequency momentum source terms */
	/* Initialize to be zero */

	/* The momentum source term is set to be v(T^4-Er) in hydro_to_fullRT */
	/* So hydro_to_fullRT must be called before this momentum source term is calculated */ 
/*	for(k=kl; k<=ku; k++){
		kg = k + koff;
		for(j=jl; j<=ju; j++){
		jg = j + joff;	
			for(i=il; i<=iu; i++){
			ig = i + ioff;
				for(m=0; m<3; m++)
					pG->Frsource[kg][jg][ig][m] = 0.0;
				
			}
		}
	}
*/

	/* Now sum over frequency */
for(ifr=0; ifr<pRG->nf; ifr++){	
	for(k=kl; k<=ku; k++){
		kg = k + koff;
		for(j=jl; j<=ju; j++){
		jg = j + joff;	
			for(i=il; i<=iu; i++){
				ig = i + ioff;
				vel[0] = pG->U[kg][jg][ig].M1 / pG->U[kg][jg][ig].d;
				vel[1] = pG->U[kg][jg][ig].M2 / pG->U[kg][jg][ig].d;
				vel[2] = pG->U[kg][jg][ig].M3 / pG->U[kg][jg][ig].d;
				VKr[0] = vel[0] * pRG->R[ifr][k][j][i].K[0] + vel[1] * pRG->R[ifr][k][j][i].K[1]
                                                + vel[2] * pRG->R[ifr][k][j][i].K[3];
				VKr[1] = vel[0] * pRG->R[ifr][k][j][i].K[1] + vel[1] * pRG->R[ifr][k][j][i].K[2]
                                                + vel[2] * pRG->R[ifr][k][j][i].K[4];
				VKr[2] = vel[0] * pRG->R[ifr][k][j][i].K[3] + vel[1] * pRG->R[ifr][k][j][i].K[4]
                                                + vel[2] * pRG->R[ifr][k][j][i].K[5]; 

	
				for(m=0; m<3; m++){
					Fr0[m] = 4.0 * PI * (pRG->R[ifr][k][j][i].H[m] - vel[m] * pRG->R[ifr][k][j][i].J/Crat 
						- VKr[m]/Crat);

					VEr[m] = 4.0 * PI * (- vel[m] * pRG->R[ifr][k][j][i].J/Crat
                                                - VKr[m]/Crat);
 
					/* sigma_a Fr term is added in hydro_to_full when T^4-Er is calculated */
					pG->Frsource[kg][jg][ig][m] += (Prat * pRG->wnu[ifr] * (pRG->R[ifr][k][j][i].Sigma[2] * Fr0[m] 
									+ pRG->R[ifr][k][j][i].Sigma[1] * VEr[m]));
				}

				/* Also add the radiation work term to the energy source term */
				/* The energy source term includes dt, but momentum source term does not */
				pG->Radheat[kg][jg][ig] += (-pG->dt * Prat * (pRG->R[ifr][k][j][i].Sigma[1] - pRG->R[ifr][k][j][i].Sigma[2]) 
							* (vel[0] * Fr0[0] + vel[1] * Fr0[1] + vel[2] * Fr0[2])); 	
				
			}/* end i */
		}/* end j */
	}/* end k */
}/* End frequency */


	/* No Need to update the boundary condition for the source terms, they are calculate in the ghost zones */


	return;
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

void SpecialMatrix3(Real *Ma, Real *Mb, Real *Mc, Real *Md, Real *RHS,  Real **lN,  Real **UN, const int N)
{

	Real swap;
	Real tempMc,temp,tempup, tempdown;
    	int i, Nzl, j; 
	
	/* First, subtract the last line from lines 1 to N-1 */
	/* RHS is also changed during this process but Xi unchanged */
	tempMc = Mc[0];
	temp = RHS[0];
	for(i=0; i<N-1; i++){
		RHS[i] -= RHS[N-1];
		/* Redefine Mc */
		Mc[i] -= Mc[N-1];
	}
	RHS[N-1] = temp;
	/* Redfined Mb */
	Mb[0] += Ma[0] + tempMc;
	for(i=1; i<N; i++)
		Mb[i] += tempMc; 

	for(i=0; i<N-1; i++)
		Md[i] = Mc[i] - Ma[N-1];

	/* Mc[N-1] is no longer used anymore */
	/* Ma[i] and Mb[0] are always no-zero */

	/* Find the first non-zero Mc */
	i=N-2; 

	while(((fabs(Mc[i]) < TINY_NUMBER))&&(i >= 0))
        	i--;

   	Nzl = i;
       
        if(Nzl < 0){
     		/* First, Multiple L^-1 and RHS */
		for(i=0; i<N-1; i++){
			RHS[N-1] -= Mb[i] * RHS[i]/Ma[i];
		}
		/* Then, Multiple U^-1 and RHS */
		/* First RHS[N-1] */
		temp = Mb[N-1];
		for(i=0; i<N-1; i++)
			temp += Mb[i] * Ma[N-1]/Ma[i];

		RHS[N-1] /= temp;

		/* Now solve RHS[0] to RHS[n-2]	*/
		for(i=0; i<N-1; i++)
			RHS[i] = (RHS[i] + Ma[N-1] * RHS[N-1])/Ma[i];		

	}/* End if Nzl < 0 */
	else if(Nzl == 0){

        	/* Exchange the first and last line */

		swap = RHS[0];
		RHS[0] = RHS[N-1];
		RHS[N-1] = swap;

		/* calculate L^-1 RHS */

		lN[0][0] = (Ma[0] + Mc[0])/Mb[0];
		for(i=1; i<N-1; i++)
			lN[0][i] = (Mc[0] - Mb[i] * (Ma[0] + Mc[0])/Mb[0])/Ma[i];

		for(i=0; i<N-1; i++)
			RHS[N-1] -= lN[0][i] * RHS[i];

		/* calculate U^-1 RHS */
		temp = -Ma[N-1] + Mc[0] - Mb[N-1] * (Ma[0] + Mc[0]) / Mb[0];
		for(i=1; i<N-1; i++)
			temp += Ma[N-1] * (Mc[0] - Mb[i] * (Ma[0] + Mc[0])/Mb[0])/Ma[i];

		RHS[N-1] /= temp;

		for(i=1; i<N-1; i++){
			RHS[i] = (RHS[i] + Ma[N-1] * RHS[N-1]) / Ma[i];
		}

		for(i=1; i<N; i++)
			RHS[0] -= Mb[i] * RHS[i];

		RHS[0] /= Mb[0];


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
				RHS[j] -= (Mc[j] * RHS[i] / Ma[i]);
				/* Redefine Md[j] */
				Md[j] -= (Mc[j] * Md[i] / Ma[i]);
			}
			/* Also replace the line N-1 */
			RHS[N-1] -= (Mb[i] * RHS[i] / Ma[i]);
			Mb[N-1] -= (Mb[i] * Md[i]/Ma[i]);
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
		
		lN[Nzl][0] = Mc[1]/Mb[0];
		lN[N-1][0] = (Ma[0] + Mc[0])/Mb[0];
		/* The second element */
		if(Nzl > 1){
			lN[Nzl][1] = (Ma[1] * Mb[0] + (Mb[0] - Mb[1]) * Mc[1])/((Mb[0] - Mb[1]) * Mc[Nzl]);
			lN[N-1][1] = -(Ma[0] * Mb[1] + Mc[0] * (Mb[1] - Mb[0]))/((Mb[0] - Mb[1]) * Mc[Nzl]);
		}

		for(i=2; i<Nzl; i++){
			lN[Nzl][i] = Ma[1] * (Mb[i] - Mb[0])/(Ma[i] * (Mb[0] - Mb[1]));
			lN[N-1][i] = Ma[0] * (Mb[1] - Mb[i])/(Ma[i] * (Mb[0] - Mb[1]));
		}
			


		/* Now the element [N-1][Nzl] */
		tempup = Mc[0] * (Mb[0] - Mb[1]);
		temp = Ma[0] * Mb[1];
		for(i=2; i<=Nzl; i++){
			tempup *= Ma[i];
			temp *= Ma[i];
		}
		tempup -= temp;

		for(j=2; j<=Nzl; j++){
			temp = Ma[0] * (Mb[1] - Mb[j]);
			for(i=2; i<=Nzl; i++){
				if(i != j)
					temp *= Ma[i];
				else
					temp *= Mc[i];
			}
			tempup -= temp;
		}

		tempdown = Mc[1] * (Mb[0] - Mb[1]);
		temp = Ma[1] * Mb[0];
		for(i=2; i<=Nzl; i++){
			tempdown *= Ma[i];
			temp *= Ma[i];
		}
		tempdown += temp;

		for(j=2; j<=Nzl; j++){
			temp = Ma[1] * (Mb[0] - Mb[j]);
			for(i=2; i<=Nzl; i++){
				if(i != j)
					temp *= Ma[i];
				else
					temp *= Mc[i];
			}
			tempdown += temp;
		}


		lN[N-1][Nzl] = tempup/tempdown;
		
		/*-----------------------------------------*/

		/* Now calculate L^-1 RHS */
		/* The first element is unchanged */
		/* The second element */
		RHS[1] -= (Mc[Nzl] * RHS[0] / Mb[0]);
		for(j=2; j<Nzl; j++){
			RHS[j] -= (Mc[j] * RHS[0]/Mb[0]);
			RHS[j] -= (Mc[j] * RHS[1]/Mc[Nzl]);
		}
		/* The line Nzl */
		if(Nzl > 1){
			for(i=0; i<Nzl; i++)
				RHS[Nzl] -= (lN[Nzl][i] * RHS[i]);
		}

		/* The line N-1 */
		for(i=0; i<=Nzl; i++)
			RHS[N-1] -= (lN[N-1][i] * RHS[i]);

		/*-----------------------------------------------*/
		/* Now construct the U matrix */
		

		/* For line 1 */

	
		for(i=1; i<=Nzl; i++)
			UN[1][i] = (Mb[0] - Mb[i]) * Mc[Nzl]/Mb[0];

		UN[1][Nzl] += Ma[Nzl];
		
		UN[1][N-1] = Md[Nzl] - Mb[N-1] * Mc[Nzl]/Mb[0];

		/* For the lines 2 to Nzl-1, the elements are a[i], ... -a[Nzl]*c[i]/c[Nzl];

		* now the line Nzl, there are only two elements in lines Nzl */
		if(Nzl > 1){
			UN[Nzl][Nzl] = -Mb[0];
			for(i=1; i<=Nzl; i++)
				UN[Nzl][Nzl] *= Ma[i];

			for(j=1; j<=Nzl; j++){
				temp = Mb[0] - Mb[j];
				for(i=1; i<=Nzl; i++){
					if(i == j)
						temp *= Mc[i];
					else
						temp *= Ma[i];
				}
				UN[Nzl][Nzl] -= temp;
			}

			temp = Mb[0] - Mb[1];
			for(i=2; i<Nzl; i++)
				temp *= Ma[i];

			temp *= Mc[Nzl];
			UN[Nzl][Nzl] /= temp;
		

		/* Now the element [Nzl][N-1] */
			UN[Nzl][N-1] = Mb[N-1] * Mc[Nzl] - Mb[0] * Md[Nzl];
			for(i=1; i<Nzl; i++)
				UN[Nzl][N-1] *= Ma[i];
		
			for(j=1; j<Nzl; j++){
				temp = Mc[Nzl] * Md[j] - Mc[j] * Md[Nzl];
				for(i=1; i<Nzl; i++){
					if(i == j)
						temp *= (Mb[0] - Mb[i]);
					else
						temp *= Ma[i];
				}
				UN[Nzl][N-1] += temp;
			}

			temp = (Mb[0] - Mb[1]) * Mc[Nzl];
			for(i=2; i<Nzl; i++)
				temp *= Ma[i];

			UN[Nzl][N-1] /= temp;
		}/* End Nzl > 1 */

		/* now the last element UN[N-1][N-1] */
		/* sum LN[N-1][i] * UN[i][N-1] = Mb */
		UN[N-1][N-1] = Md[0] - lN[N-1][0] * Mb[N-1] - lN[N-1][1] * UN[1][N-1];
		for(i=2; i<Nzl; i++)
			UN[N-1][N-1] -= (lN[N-1][i] * (Md[i] - Mc[i] * Md[Nzl]/Mc[Nzl]));

		if(Nzl > 1){
			UN[N-1][N-1] -= (lN[N-1][Nzl] * UN[Nzl][N-1]);
		}
		
		/* Now calculate U^-1 RHS */
		RHS[N-1] /= UN[N-1][N-1];
	
		RHS[Nzl] -= (UN[Nzl][N-1] * RHS[N-1]);
		RHS[Nzl] /= UN[Nzl][Nzl]; 

		for(i=Nzl-1; i>1; i--){
			RHS[i] -= (Md[i] - Mc[i] * Md[Nzl]/Mc[Nzl]) * RHS[N-1]; 
			RHS[i] += (Ma[Nzl] * Mc[i]/Mc[Nzl]) * RHS[Nzl];
			RHS[i] /= Ma[i];
		}
		
		if(Nzl > 1){	
			for(i=2; i<=Nzl; i++)
				RHS[1] -= (UN[1][i] * RHS[i]);

			RHS[1] -= (UN[1][N-1] * RHS[N-1]);
			RHS[1] /= UN[1][1];
		}
			
		for(i=1; i<=Nzl; i++)
			RHS[0] -= Mb[i] * RHS[i];

		RHS[0] -= (Mb[N-1] * RHS[N-1]);
		RHS[0] /= Mb[0];


			/* now calculate the solution between Nzl+1 to N-2 */
			/* xi =  RHS[i]/Ma[i] - Md[i]/Ma[i] xn-1 */
		for(i=Nzl+1; i<N-1; i++)
			RHS[i] = (RHS[i] - Md[i] * RHS[N-1])/Ma[i];			

		
	}/* end if Nzl > 0 */

}

void ReduceVelocity(const Real sigma, const Real ds, Real *alpha)
{

	Real tau;
	tau = 10.0 * ds * sigma;
	tau = tau * tau;

	if(tau > 0.001)
		*alpha = sqrt((1.0 - exp(- tau)) / tau);
	else
		*alpha = sqrt(1.0 - 0.5 * tau);

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


