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



static Real sigma0;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
#ifdef FULL_RADIATION_TRANSFER
static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma);
#endif
#ifndef FULL_RADIATION_TRANSFER
static Real R_ideal = 1.0;
#endif

/* Initial solution, shared with Userwork_after_loop to compute L1 error */
/* We only calculate the error from density, velocity and pressure */
/* Radiation quantities are not in cons, they are related to the gas quantities */
static ConsS ***RootSoln=NULL;


void problem(DomainS *pDomain)
{
#ifdef FULL_RADIATION_TRANSFER
  RadGridS *pRG = (pDomain->RadGrid);
#endif
  GridS *pG = (pDomain->Grid);
  int is = pG->is, ie = pG->ie;
  int js = pG->js, je = pG->je;
  int ks = pG->ks, ke = pG->ke;
#ifdef FULL_RADIATION_TRANSFER 
  int isr = pRG->is, ier = pRG->ie;
  int jsr = pRG->js, jer = pRG->je;
  int ksr = pRG->ks, ker = pRG->ke; 
  int nf=pRG->nf, nang=pRG->nang;
  int noct = pRG->noct;
#endif  
 int i, j, k, ifr, l, n, m;

  Real x1, x2, x3;
 
  Real rho0 = 1.e0, T0 = 1.0, P0, E0, amp=1.e-6;
  Real Jr, Hr, delv, dP, dEr, dFr, pre;
  Real knum, Omegareal, Omegaimg, theta, cos1, cos2, cos3, t = pG->time, weight, tlim;

/* Read problem parameters. */
#ifdef FULL_RADIATION_TRANSFER
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
#endif
	
  if (pDomain->Level == 0){
    if ((RootSoln = (ConsS***)calloc_3d_array(pG->Nx[2]+2*nghost,pG->Nx[1]+2*nghost,pG->Nx[0]+2*nghost,sizeof(ConsS)))
      == NULL) ath_error("[problem]: Error alloc memory for RootSoln\n");
  }


   P0 = rho0 * R_ideal * T0;
   E0 = P0/(Gamma - 1.0);
   Omegareal = 6.398479314825398; 
   Omegaimg = 0.7044806086435424;
   sigma0 = 1.0;

   tlim = 2.0 * PI / Omegareal;
	
#ifdef FULL_RADIATION_TRANSFER
   cos1 = pRG->mu[0][4][4][0][0];
   cos2 = pRG->mu[0][4][4][0][1];
   cos3 = pRG->mu[0][4][4][0][2];
#endif
	
	knum = 2.0 * PI;

   

	/* First, initialize gas quantities */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
       for (i=is; i<=ie; i++) {
		cc_pos(pG,i,j,k,&x1,&x2,&x3);

		theta = knum * x1 - Omegareal * t;
		delv =  amp* (1.0183496112257058 * cos(theta) + 0.1121215711780068 * sin(theta));
		dP = amp * (1.0220380692314723 * cos(theta) + 0.18993018794365163 * sin(theta));
		
		pre = P0 + dP;
	
		pG->U[k][j][i].d = rho0 + amp * cos(theta);
		pG->U[k][j][i].M1 = delv;
		pG->U[k][j][i].M2 = 0.0;
		pG->U[k][j][i].M3 = 0.0;
		

		pG->U[k][j][i].E = pre/(Gamma-1.0)  + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) 
             + SQR(pG->U[k][j][i].M3))/(pG->U[k][j][i].d);
       }
    }
   }

#ifdef FULL_RADIATION_TRANSFER

	/* Now initialize the radiation quantities */

  for (k=ksr; k<=ker; k++) {
        for (j=jsr; j<=jer; j++) {
           for (i=isr; i<=ier; i++) {
               for(ifr=0; ifr<nf; ifr++){
                   for(l=0; l<noct; l++){
                       for(n=0; n<nang; n++){
                           
		if(ker > ksr)
			cc_pos(pG,i-Radghost+nghost,j-Radghost+nghost,k-Radghost+nghost,&x1,&x2,&x3);
		else
	  		cc_pos(pG,i-Radghost+nghost,j-Radghost+nghost,0,&x1,&x2,&x3);

			theta = knum * x1 - Omegareal * t;

			dEr = amp * (-0.026018127896336885 * cos(theta) + 0.12095401964915764 * sin(theta));
			dFr = amp * (-0.10566859341556321 * cos(theta) + 0.030196412832965945 * sin(theta));

			
			Jr = (1.0 + dEr)/(4.0*PI);
			Hr = dFr/(4.0*PI*cos1);

			weight = pRG->wmu[k][j][i][l*nang+n];

			if(l == 0 || l == 2){
				pRG->imu[k][j][i][ifr][l*nang+n] = (Jr + Hr)/(4.0 * weight);
			}
			else{
				pRG->imu[k][j][i][ifr][l*nang+n] = (Jr - Hr)/(4.0 * weight);
			}
			

			
		}
	     }
	  }
       }
    }
   }

	get_full_opacity = const_opacity;

	/* set the momentums of the specific intensitites */
	UpdateRT(pDomain);
#endif
	/* set boundary condition */
/*	bvals_fullrad_trans_fun(pDomain, left_x2, TwoBeam_ix2);
*/

/* save solution on root grid */

  if (pDomain->Level == 0) {
    for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      RootSoln[k][j][i].d  = exp(-Omegaimg * tlim) * (pG->U[k][j][i].d - rho0) + rho0;
      RootSoln[k][j][i].M1 = exp(-Omegaimg * tlim) * pG->U[k][j][i].M1;
      RootSoln[k][j][i].M2 = exp(-Omegaimg * tlim) * pG->U[k][j][i].M2;
      RootSoln[k][j][i].M3 = exp(-Omegaimg * tlim) * pG->U[k][j][i].M3;
#ifndef ISOTHERMAL
      RootSoln[k][j][i].E  = exp(-Omegaimg * tlim) * (pG->U[k][j][i].E - E0) + E0;
#endif /* ISOTHERMAL */


    }}}
  }

    return;

}



void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  
  GridS *pGrid;
  int i=0,j=0,k=0;

  int is,ie,js,je,ks,ke;
  Real rms_error=0.0;
  ConsS error,total_error;
  FILE *fp;
  char *fname;
  int Nx1, Nx2, Nx3, count;
#if defined MPI_PARALLEL
  double err[12+NSCALARS], tot_err[12+NSCALARS];
  int ierr,myID;
#endif

  total_error.d = 0.0;
  total_error.M1 = 0.0;
  total_error.M2 = 0.0;
  total_error.M3 = 0.0;
#ifdef RADIATION_MHD
  total_error.B1c = 0.0;
  total_error.B2c = 0.0;
  total_error.B3c = 0.0;
#endif /* MHD */

#ifndef ISOTHERMAL
  total_error.E = 0.0;
#endif /* ISOTHERMAL */


/* Compute error only on root Grid, which is in Domain[0][0] */

  pGrid=pM->Domain[0][0].Grid;
  if (pGrid == NULL) return;

/* compute L1 error in each variable, and rms total error */

  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
    error.d = 0.0;
    error.M1 = 0.0;
    error.M2 = 0.0;
    error.M3 = 0.0;
#ifdef RADIATION_MHD
    error.B1c = 0.0;
    error.B2c = 0.0;
    error.B3c = 0.0;
#endif /* MHD */
#ifndef ISOTHERMAL
    error.E = 0.0;
#endif /* ISOTHERMAL */
   


    for (i=is; i<=ie; i++) {
      error.d   += fabs(pGrid->U[k][j][i].d   - RootSoln[k][j][i].d );
      error.M1  += fabs(pGrid->U[k][j][i].M1  - RootSoln[k][j][i].M1);
      error.M2  += fabs(pGrid->U[k][j][i].M2  - RootSoln[k][j][i].M2);
      error.M3  += fabs(pGrid->U[k][j][i].M3  - RootSoln[k][j][i].M3); 
#ifdef RADIATION_MHD
      error.B1c += fabs(pGrid->U[k][j][i].B1c - RootSoln[k][j][i].B1c);
      error.B2c += fabs(pGrid->U[k][j][i].B2c - RootSoln[k][j][i].B2c);
      error.B3c += fabs(pGrid->U[k][j][i].B3c - RootSoln[k][j][i].B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
      error.E   += fabs(pGrid->U[k][j][i].E   - RootSoln[k][j][i].E );
#endif /* ISOTHERMAL */      
    }

    total_error.d += error.d;
    total_error.M1 += error.M1;
    total_error.M2 += error.M2;
    total_error.M3 += error.M3;
#ifdef RADIATION_MHD
    total_error.B1c += error.B1c;
    total_error.B2c += error.B2c;
    total_error.B3c += error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E += error.E;
#endif /* ISOTHERMAL */
  
  }}

#ifdef MPI_PARALLEL
  Nx1 = pM->Domain[0][0].Nx[0];
  Nx2 = pM->Domain[0][0].Nx[1];
  Nx3 = pM->Domain[0][0].Nx[2];
#else
  Nx1 = ie - is + 1;
  Nx2 = je - js + 1;
  Nx3 = ke - ks + 1;
#endif
  count = Nx1*Nx2*Nx3;

#ifdef MPI_PARALLEL 
/* Now we have to use an All_Reduce to get the total error over all the MPI
 * grids.  Begin by copying the error into the err[] array */

  err[0] = total_error.d;
  err[1] = total_error.M1;
  err[2] = total_error.M2;
  err[3] = total_error.M3;
#ifdef RADIATION_MHD
  err[4] = total_error.B1c;
  err[5] = total_error.B2c;
  err[6] = total_error.B3c;
#endif /* MHD */
#ifndef ISOTHERMAL
  err[7] = total_error.E;
#endif /* ISOTHERMAL */


/* Sum up the Computed Error */
  ierr = MPI_Reduce(err,tot_err,(12+NSCALARS),MPI_DOUBLE,MPI_SUM,0,
    pM->Domain[0][0].Comm_Domain);

/* If I'm the parent, copy the sum back to the total_error variable */

  ierr = MPI_Comm_rank(pM->Domain[0][0].Comm_Domain, &myID);
  if(myID == 0){ /* I'm the parent */
    total_error.d   = tot_err[0];
    total_error.M1  = tot_err[1];
    total_error.M2  = tot_err[2];
    total_error.M3  = tot_err[3];
#ifdef RADIATION_MHD
    total_error.B1c = tot_err[4];
    total_error.B2c = tot_err[5];
    total_error.B3c = tot_err[6];
#endif /* MHD */
#ifndef ISOTHERMAL
    total_error.E   = tot_err[7];
#endif /* ISOTHERMAL */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) total_error.s[n] = err[8+n];
#endif

  }
  else return; /* The child grids do not do any of the following code */

#endif /* MPI_PARALLEL */

/* Compute RMS error over all variables, and print out */

  rms_error = SQR(total_error.d) + SQR(total_error.M1) + SQR(total_error.M2)
                + SQR(total_error.M3);
#ifdef RADIATION_MHD
  rms_error += SQR(total_error.B1c) + SQR(total_error.B2c) 
               + SQR(total_error.B3c);
#endif /* MHD */
#ifndef ISOTHERMAL
  rms_error += SQR(total_error.E);
#endif /* ISOTHERMAL */
 
  rms_error = sqrt(rms_error)/(double)count;


/* Print error to file "LinWave-errors.#.dat", where #=wave_flag  */

#ifdef MPI_PARALLEL
  fname = ath_fname("../","LinWave-errors",NULL,NULL,1,NULL,NULL,"dat");
#else
  fname = ath_fname(NULL,"LinWave-errors",NULL,NULL,1,NULL,NULL,"dat");
#endif

/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_after_loop]: Unable to reopen file.\n");
      free(fname);
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_after_loop]: Unable to open file.\n");
      free(fname);
      return;
    }
/* Now write out some header information */
    fprintf(fp,"# Nx1  Nx2  Nx3  RMS-Error  d  M1  M2  M3");
#ifndef ISOTHERMAL
    fprintf(fp,"  E");
#endif /* ISOTHERMAL */
#ifdef RADIATION_MHD
    fprintf(fp,"  B1c  B2c  B3c");
#endif /* MHD */
  

    fprintf(fp,"\n#\n");
  }

  fprintf(fp,"%d  %d  %d  %e",Nx1,Nx2,Nx3,rms_error);

  fprintf(fp,"  %e  %e  %e  %e",
	  (total_error.d/(double) count),
	  (total_error.M1/(double)count),
	  (total_error.M2/(double)count),
	  (total_error.M3/(double)count));

#ifndef ISOTHERMAL
  fprintf(fp,"  %e",(total_error.E/(double)count));
#endif /* ISOTHERMAL */

#ifdef RADIATION_MHD
  fprintf(fp,"  %e  %e  %e",
	  (total_error.B1c/(double)count),
	  (total_error.B2c/(double)count),
	  (total_error.B3c/(double)count));
#endif /* MHD */
  

  fprintf(fp,"\n");

  fclose(fp);
  free(fname);

  return;


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
	fwrite(&Gamma,sizeof(Real),1,fp);

#ifdef FULL_RADIATION_TRANSFER
	fwrite(&Prat,sizeof(Real),1,fp);
	fwrite(&Crat,sizeof(Real),1,fp);
	fwrite(&R_ideal,sizeof(Real),1,fp); 	
#endif	

}


void problem_read_restart(MeshS *pM, FILE *fp)
{

	fread(&Gamma,sizeof(Real),1,fp);

#ifdef FULL_RADIATION_TRANSFER
	fread(&Prat,sizeof(Real),1,fp);
	fread(&Crat,sizeof(Real),1,fp);
	fread(&R_ideal,sizeof(Real),1,fp); 	
#endif	

#ifdef FULL_RADIATION_TRANSFER
	get_full_opacity = const_opacity;
	
#endif
	/* Increase the background magnetic field */
	DomainS *pD;
	pD= &(pM->Domain[0][0]);

	/* set boundary condition */
/*	bvals_fullrad_trans_fun(pD, left_x2, TwoBeam_ix2);
*/

}

#ifdef FULL_RADIATION_TRANSFER
static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma)
{

	Sigma[0] = sigma0;
	Sigma[1] = sigma0;
	Sigma[2] = 0.0;
	Sigma[3] = 0.0;

  return;
  
}

#endif



