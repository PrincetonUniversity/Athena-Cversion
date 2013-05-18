#include "copyright.h"

/*==============================================================================
 * FILE: FullRTtest.c
 *
 * PURPOSE:  
 * Initial conditions available:
 *
 *============================================================================*/

#include <math.h>

#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifndef FULL_RADIATION_TRANSFER

#error FullRTtest.c requires FULL_RADIATION_TRANSFER enabled!

#endif




static Real eps0;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/

static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma);

void problem(DomainS *pDomain)
{
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
  int isr = pRG->is, ier = pRG->ie;
  int jsr = pRG->js, jer = pRG->je;
  int ksr = pRG->ks, ker = pRG->ke; 
  int nf=pRG->nf, nang=pRG->nang;
  int noct = pRG->noct;
  int i, j, k, ifr, l, n, m;

  Real x1, x2, x3;
 
  Real rho = 1.0, T;

  Real Jr;

/* Read problem parameters. */

  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  R_ideal = par_getd("problem","R_ideal");
#ifdef MPI_PARALLEL 
  if(myID_Comm_world == 0){
#endif
     printf("Parameters: Prat %G Crat %G R_ideal %G\n",Prat,Crat,R_ideal);
#ifdef MPI_PARALLEL 
  }
#endif

	/* First, initialize gas quantities */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
       for (i=is; i<=ie; i++) {
		cc_pos(pG,i,j,k,&x1,&x2,&x3);

		pG->U[k][j][i].d = rho;
		pG->U[k][j][i].M1 = 0.0;
		pG->U[k][j][i].M2 = 0.0;
		pG->U[k][j][i].M3 = 0.0;
		if(sqrt(SQR(x1-0.5)+SQR(x2-0.5))<0.2)
			T = 1.0;
		else
			T = 1.0;

		pG->U[k][j][i].E = rho*T/(Gamma-1.0)  + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) 
             + SQR(pG->U[k][j][i].M3))/rho;
       }
    }
   }


	/* Now initialize the radiation quantities */
for(ifr=0; ifr<nf; ifr++){
  for(l=0; l<noct; l++){
    for(n=0; n<nang; n++){
      for (k=ksr; k<=ker; k++) {
        for (j=jsr; j<=jer; j++) {
           for (i=isr; i<=ier; i++) {
		cc_pos(pG,i-Radghost+nghost,j-Radghost+nghost,0,&x1,&x2,&x3);
	  
		/*	if(sqrt(SQR(x1-0.5)+SQR(x2-0.5))<0.2)
				pRG->imu[k][j][i][ifr][l][n] = 10.0/(4.0*PI);
			else
				pRG->imu[k][j][i][ifr][l][n] = 1.0/(4.0*PI);
		*/

			if(sqrt(SQR(x1-0.5)+SQR(x2-0.5))<0.2)
				Jr = 10.0/(4.0*PI);
			else
				Jr = 1.0/(4.0*PI);

			if((n == 0 && l == 1) )
				pRG->imu[ifr][l][n][k][j][i] = Jr;
			else
				pRG->imu[ifr][l][n][k][j][i] = 0.0;
			
		}
	     }
	  }
       }
    }
   }

	get_full_opacity = const_opacity;

	/* set the momentums of the specific intensitites */
	UpdateRT(pDomain);


    return;

}



void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}


/* Get_user_expression computes dVy */
ConsFun_t get_usr_expr(const char *expr)
{
 /* if(strcmp(expr,"dVy")==0) return expr_dV2;
*/
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name)
{
/*  if(strcmp(name,"1d")==0) return output_1d;
  if(strcmp(name,"1dx")==0) return output_1dx;
*/
	return NULL;
}


void problem_write_restart(MeshS *pM, FILE *fp)
{

}


void problem_read_restart(MeshS *pM, FILE *fp)
{

}

static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma)
{

	Sigma[0] = 0.0;
	Sigma[1] = 0.0;
	Sigma[2] = 0.0;
	Sigma[3] = 0.0;

  return;
  
}

