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
static Real kappaes;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/

static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma);

static void TwoBeam_ix2(GridS *pG, RadGridS *pRG);

static void output_1dx(MeshS *pM, OutputS *pOut);
static void output_1dy(MeshS *pM, OutputS *pOut);

#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k);
#endif
static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j,const int k);
static Real expr_KE(const GridS *pG, const int i, const int j, const int k);


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

  Real x1, x2, x3, ytop;
 
  kappaes = 4.e4;

  Real rho0 = 1.e-3, T;
  Real rho;
  Real Jr, slope1, slope2, dis1, dis2;
  Real Hr, cos1,cos2, weight, Er, Fr;

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

  ytop = pDomain->RootMaxX[1];
  cos1 = pRG->mu[0][0][0][4][4][0];
  cos2 = pRG->mu[0][0][0][4][4][1];
	/* First, initialize gas quantities */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
       for (i=is; i<=ie; i++) {
		cc_pos(pG,i,j,k,&x1,&x2,&x3);
		rho = rho0 * exp(fabs(x2-ytop));
		rho = 1.0;
		pG->U[k][j][i].d = rho;
		pG->U[k][j][i].M1 = 0.0;
		pG->U[k][j][i].M2 = 1.0;
		pG->U[k][j][i].M3 = 0.0;
		if(sqrt(SQR(x1-0.5)+SQR(x2-0.5)+SQR(x3-0.5))<0.2)
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
		if(ker > ksr)
			cc_pos(pG,i-Radghost+nghost,j-Radghost+nghost,k-Radghost+nghost,&x1,&x2,&x3);
		else
	  		cc_pos(pG,i-Radghost+nghost,j-Radghost+nghost,0,&x1,&x2,&x3);

			Er = exp(-40.0*x2*x2);
			Fr= 2.0*40.0*x2*exp(-40.0*x2*x2)/(3.0*kappaes);

			if(fabs(x2) > 0.5){
				Er = exp(-40.0*0.5*0.5);
				Fr = 0.0;
			}

			Jr = Er/(4.0*PI);
			Hr = Fr/(4.0*PI*cos2);

			weight = pRG->wmu[n][k][j][i];

			if(l == 0 || l == 1){
				pRG->imu[ifr][l][n][k][j][i] = (Jr + 0.0 * Hr)/(4.0 * weight);
			}
			else{
				pRG->imu[ifr][l][n][k][j][i] = (Jr - 0.0 * Hr)/(4.0 * weight);
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

	/* set boundary condition */
/*	bvals_fullrad_trans_fun(pDomain, left_x2, TwoBeam_ix2);
*/

	/* Enroll the user defined function to dump history of specific intensity */
	Intensity_history_enroll(0, 0, 0, "<I_l0_n0>");

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
  if(strcmp(name,"1dy")==0) return output_1dy;

  if(strcmp(name,"1dx")==0) return output_1dx;

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


	get_full_opacity = const_opacity;
	
	/* Increase the background magnetic field */
	DomainS *pD;
	pD= &(pM->Domain[0][0]);

	/* set boundary condition */
/*	bvals_fullrad_trans_fun(pD, left_x2, TwoBeam_ix2);
*/

}

static void const_opacity(GridS *pG, const int ifr, const int i,
			     const int j, const int k, Real *Sigma)
{
	Real rho;
	rho = pG->U[k][j][i].d;

	Sigma[0] = 0.0;
	Sigma[1] = 0.0;
	Sigma[2] = kappaes * rho;
	Sigma[3] = kappaes * rho;

  return;
  
}

static void TwoBeam_ix2(GridS *pG, RadGridS *pRG)
{


  int is = pRG->is, ie = pRG->ie;
  int js = pRG->js;
  int ks = pRG->ks, ke = pRG->ke;
  int nang = pRG->nang;
  int noct = pRG->noct;
  int nf = pRG->nf;
  int i, j, k, l, n, ifr;
  int ig, jg, kg;
  Real Jr, slope1, slope2, dis1, dis2;
  Real x1, x2, x3;
  int ioff, joff, koff;


  ioff = nghost - Radghost;
 
  if(pG->Nx[1] > 1) joff = nghost - Radghost;
  if(pG->Nx[2] > 1) koff = nghost - Radghost;

  

  Jr = 10.0/(4.0*PI);


for(ifr=0; ifr<nf; ifr++){
   for(l=0; l<noct; l++){
      for(n=0; n<nang; n++){
	 for(k=ks; k<=ke; k++){
		kg = k + koff;
       	    for(j=1; j<=Radghost; j++){
		jg = js-j + joff;
	       for(i=is-Radghost; i<=ie+Radghost; i++){
		ig = i + ioff;	
			cc_pos(pG,ig,jg,kg,&x1, &x2, &x3);
			slope1 = -pRG->mu[0][0][k][j][i][1]/pRG->mu[0][0][k][j][i][0];
			slope2 = -pRG->mu[1][0][k][j][i][1]/pRG->mu[1][0][k][j][i][0];
			dis1 = fabs(slope1 * (x1 - 0.1) + (x2 + 0.5));
			dis2 = fabs(slope2 * (x1 + 0.1) + (x2 + 0.5));
		

			if((k==(ks+ke)/2) && (((l == 0) && (n ==0) && (dis1 < pG->dx1)) || ((l == 1) && (n ==0) && (dis2 < pG->dx1)))){
				pRG->imu[ifr][l][n][k][js-j][i] = Jr;				
			}
			else
				pRG->imu[ifr][l][n][k][js-j][i] = 0.0;



	   	}/* end i */
       	     }/* end J */
     	 } /* End k */
       }/* end nang */
    }/* end noctant */  
}/* end ifr */


}



static Real hst_rho_Vx_dVy(const GridS *pG,const int i,const int j, const int k)
{
  Real x1,x2,x3;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  return pG->U[k][j][i].M1*(pG->U[k][j][i].M2/pG->U[k][j][i].d);

}

#ifdef ADIABATIC
static Real hst_E_total(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,phi, Kvy, Vy0, dVy;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  return pG->U[k][j][i].E;
}
#endif /* ADIABATIC */

/*----------------------------------------------------------------------------*/
/*! \fn static Real expr_KE(const GridS *pG, const int i, const int j, 
 *			    const int k)
 *  \brief Computes dens*(Vx^2+Vy^2+Vz^2)/2
 */
static Real expr_KE(const GridS *pG, const int i, const int j, const int k)
{
  Real x1,x2,x3,Vy,Vx,Vz;
  cc_pos(pG,i,j,k,&x1,&x2,&x3);

  Vy = (pG->U[k][j][i].M2/pG->U[k][j][i].d);

  Vx = pG->U[k][j][i].M1/pG->U[k][j][i].d;
  Vz = pG->U[k][j][i].M3/pG->U[k][j][i].d;

  return pG->U[k][j][i].d*(Vx*Vx + Vy*Vy + Vz*Vz)/2.0;

}


/*! \fn static void output_1dx(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
    averaged quantities.  Currently, only outputs at lowest
    refinement level */

static void output_1dx(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
#ifdef FULL_RADIATION_TRANSFER
  RadGridS *pRG;
#endif
  int i,j,k;
  int tot1d,i1d,nzmx,my_nz,kg,kdisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lx,Ly,Lz,press, press1, press3, Bpre1, Bpre3;
  static double *out_x3;
  double vx, vy, vz, Fr01,Fr02,Fr03;
  int flag;

  FILE *p_1dfile;
  char *fname;
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int zproc;
  int ierr,myID_Comm_Domain;
#endif

/* For radiation case, we add, Er, Frx, Fry, Frz, */

#ifdef FULL_RADIATION_TRANSFER
  int koff, joff, ioff;
  ioff = Radghost - nghost;
  if(pM->Nx[1] > 1) joff = ioff;
  else		    joff = 0;

  if(pM->Nx[2] > 1) koff = ioff;
  else	            koff = 0;

#endif

#if defined(MHD)
  tot1d=15+8-3;
#else
  tot1d=15-3;
#endif /* MHD */
	tot1d++;
#ifdef ADIABATIC
  tot1d=tot1d+3;
#endif /* ADIABATIC */

  Lx = pM->RootMaxX[0] - pM->RootMinX[0];
  Ly = pM->RootMaxX[1] - pM->RootMinX[1];
  Lz = pM->RootMaxX[2] - pM->RootMinX[2];
  nzmx = pM->Nx[0];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef FULL_RADIATION_TRANSFER
  pRG=pD->RadGrid;

#endif

#ifdef MPI_PARALLEL
  int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] = 0.0;
    }
  }
  kdisp=pGrid->Disp[0];

/* First calculate the x3 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (k=0; k<nzmx; k++) {
      x1 = pM->RootMinX[0] + (k + 0.5)*pGrid->dx1;
      out_x3[k] = x1;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (i=is; i<=ie; i++) {
    kg=i+kdisp-nghost;
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        i1d=0;
        out1d[kg][i1d] += pGrid->U[k][j][i].d;
        i1d++;
#ifdef ISOTHERMAL
        out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
        press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
                                 - expr_ME(pGrid,i,j,k)
#endif
                                ),TINY_NUMBER);
        out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[kg][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
        i1d++;
        out1d[kg][i1d] += pGrid->U[k][j][i].E;
        i1d++;
        out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
        i1d++;
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
	
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;

        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
        i1d++;
        out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
        i1d++;
        out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
        i1d++;
        out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bx(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_By(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_Bz(pGrid,i,j,k);
        i1d++;
        out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif

#ifdef FULL_RADIATION_TRANSFER
	i1d++;
        out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].J);
	i1d++;
        out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[0]);
	i1d++;
        out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[1]);
	i1d++;
	/* To avoid cancel */
	
        out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[2]);
		  i1d++;
		  /* To avoid cancel */
		  
		  out1d[kg][i1d] += 4.0*PI*(pRG->imu[0][0][0][k+koff][j+joff][i+ioff]);
		  

#endif

      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lz*Ly/(pGrid->dx3*pGrid->dx2);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (k=0; k<nzmx; k++) {
      my_out1d[k] = out1d[k][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
    for (k=0; k<nzmx; k++) {
      out1d[k][i1d] = g_out1d[k];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (k=0; k<nzmx; k++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[k][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dx");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dx");
#endif
  if (fname == NULL) {
    ath_error("[output_1d]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
#ifdef MHD
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6]);
#endif /* MHD */
#else
#ifdef RADIATION_MHD
    if (k == 0) {
      fprintf(p_1dfile,"# [1]x3     [2]dens    [3]pressure    [4]temperature  [5]E     [6]Etot     [7]KEx         [8]KEy        [9] KEz       [10] KE        [11]Reynolds   [12]MEx        [13]MEy        [14]MEz        [15]ME         [16]Bx          [17]By         [18]Bz         [19]Maxwell     [20]Er      [21]Frx      [22]Fry       [23]Frz     [24]Frz0	[25]dFr0dz	[26]ErV		[27]dP/dz/rho		[28]dBpre/dz/rho	[29]kappaes          [30]kappaff\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
            out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
            out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28]);
#else
    if (k == 0) {
      fprintf(p_1dfile,"# [1]x3		[2]dens		[3]pressure	[4]temperature	[5]E	[6]Etot		[7]KEx		[8]KEy		[9]KEz	[10]KE	[11]Reynolds	[12]Er	[13]Frx		[14]Fry		[15]Frz		[16]Il0n0\n");
    }
    fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
            out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],out1d[k][12],out1d[k][13],out1d[k][14]);
#endif /* RADIATION_MHD */
#endif /* ISOTHERMAL */
  }

  fclose(p_1dfile);
  free(fname);
#ifdef MPI_PARALLEL
  }
#endif

  free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
  free_1d_array(my_out1d); /* Free the memory we malloc'd */
  free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
  if (FIRST == 0) {
    FIRST = 1;
  }

return;
}




/*! \fn static void output_1dx(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D horizontally
 averaged quantities.  Currently, only outputs at lowest
 refinement level */

static void output_1dy(MeshS *pM, OutputS *pOut)
{
	GridS *pGrid;
	DomainS *pD;
#ifdef FULL_RADIATION_TRANSFER
	RadGridS *pRG;
#endif
	int i,j,k;
	int tot1d,i1d,nzmx,my_nz,kg,kdisp;
	int dnum = pOut->num,nl,nd;
	static int FIRST = 0;
	double darea,**out1d;
	double x1,x2,x3,Lx,Ly,Lz,press, press1, press3, Bpre1, Bpre3;
	static double *out_x3;
	double vx, vy, vz, Fr01,Fr02,Fr03;
	int flag;
	
	FILE *p_1dfile;
	char *fname;
	double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */
	
#ifdef MPI_PARALLEL
	double *my_out1d;
	double *g_out1d;
	int zproc;
	int ierr,myID_Comm_Domain;
#endif
	
	/* For radiation case, we add, Er, Frx, Fry, Frz, */
	
#ifdef FULL_RADIATION_TRANSFER
	int koff, joff, ioff;
	ioff = Radghost - nghost;
	if(pM->Nx[1] > 1) joff = ioff;
	else		    joff = 0;
	
	if(pM->Nx[2] > 1) koff = ioff;
	else	            koff = 0;
	
#endif
	
#if defined(MHD)
	tot1d=15+8-3;
#else
	tot1d=15-3;
#endif /* MHD */
	tot1d++;
#ifdef ADIABATIC
	tot1d=tot1d+3;
#endif /* ADIABATIC */
	
	Lx = pM->RootMaxX[0] - pM->RootMinX[0];
	Ly = pM->RootMaxX[1] - pM->RootMinX[1];
	Lz = pM->RootMaxX[2] - pM->RootMinX[2];
	nzmx = pM->Nx[1];
	
	/* At level=0, there is only one domain */
	
	pGrid = pM->Domain[0][0].Grid;
	int is = pGrid->is, ie = pGrid->ie;
	int js = pGrid->js, je = pGrid->je;
	int ks = pGrid->ks, ke = pGrid->ke;
	pD = (DomainS*)&(pM->Domain[0][0]);
	
#ifdef FULL_RADIATION_TRANSFER
	pRG=pD->RadGrid;
	
#endif
	
#ifdef MPI_PARALLEL
	int nproc = pD->NGrid[0]*pD->NGrid[1]*pD->NGrid[2];
#endif
	
#ifdef MPI_PARALLEL
	ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
	if(ierr != MPI_SUCCESS)
		ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
	if (FIRST == 0){
#ifdef MPI_PARALLEL
		if (myID_Comm_Domain == 0) {
#endif
			out_x3 = (double *) calloc_1d_array(nzmx,sizeof(double));
#ifdef MPI_PARALLEL
		}
#endif
	}
	
	out1d = (double **) calloc_2d_array(nzmx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
	my_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
	g_out1d = (double *) calloc_1d_array(nzmx,sizeof(double));
#endif
	for (k=0; k<nzmx; k++) {
		for (i1d=0; i1d<tot1d; i1d++) {
			out1d[k][i1d] = 0.0;
		}
	}
	kdisp=pGrid->Disp[1];
	
	/* First calculate the x3 coordinate and save it to be dumped
	 by root in every 1d file */
	if (FIRST == 0) {
#ifdef MPI_PARALLEL
		if (myID_Comm_Domain == 0) {
#endif
			for (k=0; k<nzmx; k++) {
				x1 = pM->RootMinX[1] + (k + 0.5)*pGrid->dx2;
				out_x3[k] = x1;
			}
#ifdef MPI_PARALLEL
		}
#endif
	}
	
	/* Compute 1d averaged variables */
	for (j=js; j<=je; j++) {
		kg=j+kdisp-nghost;
		for (k=ks; k<=ke; k++) {
			for (i=is; i<=ie; i++) {
				i1d=0;
				out1d[kg][i1d] += pGrid->U[k][j][i].d;
				i1d++;
#ifdef ISOTHERMAL
				out1d[kg][i1d] += pGrid->U[k][j][i].d*Iso_csound2;
#else
				press           = MAX(Gamma_1*(pGrid->U[k][j][i].E - expr_KE(pGrid,i,j,k)
#if defined(MHD) || defined(RADIATION_MHD)
											   - expr_ME(pGrid,i,j,k)
#endif
											   ),TINY_NUMBER);
				out1d[kg][i1d] += press;
#endif
#ifdef ADIABATIC
				i1d++;
				out1d[kg][i1d] += press/(R_ideal * pGrid->U[k][j][i].d);
				i1d++;
				out1d[kg][i1d] += pGrid->U[k][j][i].E;
				i1d++;
				out1d[kg][i1d] += hst_E_total(pGrid,i,j,k);
#endif
				i1d++;
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d;
				i1d++;
				cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
				
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;
				
				i1d++;
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d;
				i1d++;
				out1d[kg][i1d] += expr_KE(pGrid,i,j,k);
				i1d++;
				out1d[kg][i1d] += hst_rho_Vx_dVy(pGrid,i,j,k);
#if defined(MHD) || defined(RADIATION_MHD)
				i1d++;
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c);
				i1d++;
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c);
				i1d++;
				out1d[kg][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c);
				i1d++;
				out1d[kg][i1d] += expr_ME(pGrid,i,j,k);
				i1d++;
				out1d[kg][i1d] += hst_Bx(pGrid,i,j,k);
				i1d++;
				out1d[kg][i1d] += hst_By(pGrid,i,j,k);
				i1d++;
				out1d[kg][i1d] += hst_Bz(pGrid,i,j,k);
				i1d++;
				out1d[kg][i1d] += hst_BxBy(pGrid,i,j,k);
#endif
				
#ifdef FULL_RADIATION_TRANSFER
				i1d++;
				out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].J);
				i1d++;
				out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[0]);
				i1d++;
				out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[1]);
				i1d++;
				/* To avoid cancel */
				
				out1d[kg][i1d] += 4.0*PI*(pRG->R[0][k+koff][j+joff][i+ioff].H[2]);
				i1d++;
				/* To avoid cancel */
				
				out1d[kg][i1d] += 4.0*PI*(pRG->imu[0][0][0][k+koff][j+joff][i+ioff]);
				
				
#endif
				
			}
		}
	}
	
	/* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
	area_rat = Lz*Ly/(pGrid->dx3*pGrid->dx1);
	
	/* The parent sums the scal[] array.
	 * Note that this assumes (dx1,dx2,dx3) = const. */
	
#ifdef MPI_PARALLEL 
	for(i1d=0; i1d<tot1d; i1d++){
		for (k=0; k<nzmx; k++) {
			my_out1d[k] = out1d[k][i1d];
		}
		ierr = MPI_Reduce(my_out1d, g_out1d, nzmx,
						  MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
		if(ierr)
			ath_error("[output_1d]: MPI_Reduce call returned error = %d\n",ierr);
		for (k=0; k<nzmx; k++) {
			out1d[k][i1d] = g_out1d[k];
		}
	}
#endif
	
	/* For parallel calculations, only the parent computes the average
	 * and writes the output. */
#ifdef MPI_PARALLEL
	if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif
		
		darea = 1.0/(double)area_rat;
		for (k=0; k<nzmx; k++) {
			for (i1d=0; i1d<tot1d; i1d++) {
				out1d[k][i1d] *= darea;
			}
		}
		
		/* Generate filename */
#ifdef MPI_PARALLEL
		fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dy");
#else
		fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"1dy");
#endif
		if (fname == NULL) {
			ath_error("[output_1d]: Error constructing output filename\n");
			return;
		}
		
		/* open filename */
		p_1dfile = fopen(fname,"w");
		if (p_1dfile == NULL) {
			ath_error("[output_1d]: Unable to open 1d average file %s\n",fname);
			return;
		}
		
		/* Write out data */
		
		for (k=0; k<nzmx; k++) {
#ifdef ISOTHERMAL
#ifdef MHD
			if (k == 0) {
				fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds    MEx         MEy         MEz         ME          Bx           By           Bz          Maxwell\n");
			}
			fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
					out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
					out1d[k][12],out1d[k][13],out1d[k][14]);
#else
			if (k == 0) {
				fprintf(p_1dfile,"# x3     dens  pressure    KEx         KEy         KEz         KE          Reynolds\n");
			}
			fprintf(p_1dfile,"%G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
					out1d[k][5],out1d[k][6]);
#endif /* MHD */
#else
#ifdef RADIATION_MHD
			if (k == 0) {
				fprintf(p_1dfile,"# [1]x3     [2]dens    [3]pressure    [4]temperature  [5]E     [6]Etot     [7]KEx         [8]KEy        [9] KEz       [10] KE        [11]Reynolds   [12]MEx        [13]MEy        [14]MEz        [15]ME         [16]Bx          [17]By         [18]Bz         [19]Maxwell     [20]Er      [21]Frx      [22]Fry       [23]Frz     [24]Frz0	[25]dFr0dz	[26]ErV		[27]dP/dz/rho		[28]dBpre/dz/rho	[29]kappaes          [30]kappaff\n");
			}
			fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],
					out1d[k][3],out1d[k][4],out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],
					out1d[k][12],out1d[k][13],out1d[k][14],out1d[k][15],out1d[k][16],out1d[k][17],out1d[k][18],out1d[k][19],out1d[k][20],out1d[k][21],out1d[k][22],out1d[k][23],out1d[k][24],out1d[k][25],out1d[k][26],out1d[k][27],out1d[k][28]);
#else
			if (k == 0) {
				fprintf(p_1dfile,"# [1]x3		[2]dens		[3]pressure	[4]temperature	[5]E	[6]Etot		[7]KEx		[8]KEy		[9]KEz	[10]KE	[11]Reynolds	[12]Er	[13]Frx		[14]Fry		[15]Frz		[16]Il0n0\n");
			}
			fprintf(p_1dfile,"%G %G %G %G %G %G %G %G %G %G %G %G %G %G %G %G\n",out_x3[k],out1d[k][0],out1d[k][1],out1d[k][2],out1d[k][3],out1d[k][4],
					out1d[k][5],out1d[k][6],out1d[k][7],out1d[k][8],out1d[k][9],out1d[k][10],out1d[k][11],out1d[k][12],out1d[k][13],out1d[k][14]);
#endif /* RADIATION_MHD */
#endif /* ISOTHERMAL */
		}
		
		fclose(p_1dfile);
		free(fname);
#ifdef MPI_PARALLEL
	}
#endif
	
	free_2d_array(out1d); /* Free the memory we malloc'd */
#ifdef MPI_PARALLEL
	free_1d_array(my_out1d); /* Free the memory we malloc'd */
	free_1d_array(g_out1d); /* Free the memory we malloc'd */
#endif
	if (FIRST == 0) {
		FIRST = 1;
	}
	
	return;
}



