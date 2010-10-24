#include "copyright.h"
/*==============================================================================
 * FILE: Marshak.c
 *
 * PURPOSE: Problem generator to test the radiation MHD code. 
 *  Compare the numerical results from the semi-analytical solution fro Su & Olson 1996
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/
/* problem:    */
void radMHD_inflow(GridS *pGrid);
void radMHD_inflow2(GridS *pGrid);
void radMHD_rad_inflow(GridS *pGrid);
void radMHD_rad_inflow2(GridS *pGrid);
double qromo(double (*func)(double,double), double a, double b,
	double (*choose)(double(*)(double,double), double, double, int,double),double t);
double midpnt(double (*func)(double,double), double a, double b, int n, double t);
double func1(double x, double t);
double func2(double x, double t);
double func3(double x, double t);
double func4(double x, double t);

void problem(DomainS *pDomain)
{
  GridS *pGrid=(pDomain->Grid);
  int i, j, k, iu, il, ju, jl, ku, kl;
  int shift;

/* Parse global variables of unit ratio */
#ifdef rad_hydro
  Prat = par_getd("problem","Pratio");
  Crat = par_getd("problem","Cratio");
  Sigma_t = par_getd("problem","Sigma_t");
  Sigma_a = par_getd("problem","Sigma_a");
  R_ideal = par_getd("problem","R_ideal");
#endif

/* Set up the index bounds for initializing the grid */
  iu = pGrid->ie;
  il = pGrid->is;

  if (pGrid->Nx[1] > 1) {
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  }
  else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  if (pGrid->Nx[2] > 1) {
    ku = pGrid->ke + nghost;
    kl = pGrid->ks - nghost;
  }
  else {
    ku = pGrid->ke;
    kl = pGrid->ks;
  }

/* Initialize the grid including the ghost cells.  */
	Real d0, u0, T0, x1, x2, x3, temperature;
	d0 = 1.0;
	u0 = 0.0;
	T0 = 0.0;


    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
        for (i=il; i<=iu; i++) {

	cc_pos(pGrid, i, j,k, &x1, &x2, &x3);

/* Initialize conserved (and  the primitive) variables in Grid */
	
	  temperature = T0;
          pGrid->U[k][j][i].d  = d0;
          pGrid->U[k][j][i].M1 = d0 * u0;
          pGrid->U[k][j][i].M2 = 0.0;
          pGrid->U[k][j][i].M3 = 0.0;

          pGrid->U[k][j][i].E = 0.5 * d0 * u0 * u0 + d0 * temperature /(Gamma - 1.0);

	 pGrid->U[k][j][i].Edd_11 = 0.33333333; /* Set to be a constant in 1D. To be modified later */

#ifdef MHD
          pGrid->B1i[k][j][i] = 0.0;
          pGrid->B2i[k][j][i] = 0.0;
          pGrid->B3i[k][j][i] = 0.0;
          pGrid->U[k][j][i].B1c = 0.0;
          pGrid->U[k][j][i].B2c = 0.0;
          pGrid->U[k][j][i].B3c = 0.0;
#endif
#ifdef rad_hydro
	  pGrid->U[k][j][i].Er =  0.0;
	  pGrid->U[k][j][i].Fr1 = 0.0;
	  pGrid->U[k][j][i].Fr2 = 0.0;
	  pGrid->U[k][j][i].Fr3 = 0.0;

	 		
#endif
        }
      }
    }

	bvals_mhd_fun(pDomain, left_x1, radMHD_inflow);
	bvals_mhd_fun(pDomain, right_x1, radMHD_inflow2);
	bvals_rad_fun(pDomain, left_x1, radMHD_rad_inflow);
	bvals_rad_fun(pDomain, right_x1, radMHD_rad_inflow2);

  return;
}

void radMHD_inflow(GridS *pGrid)
{
 
	int i, is;
        int ks, js;
        is = pGrid->is;
        ks = pGrid->ks;
        js = pGrid->js;

        double dt=pGrid->dt;
        double epsilon = 0.1;
        double tempEr, tempT, updateT;
	double t=pGrid->time;
	double DimT;
	DimT = epsilon * Crat * Sigma_a * t;

	for (i=1;  i<=nghost;  i++) {


        	tempEr = pGrid->U[ks][js][is-i].Er;
/*	        tempT = pGrid->U[ks][js][is-i].E * (Gamma - 1.0); */
 		updateT = tempEr - qromo(func3,0.0,1.0,midpnt,DimT) + qromo(func4,0.0,1.0,midpnt,DimT);
		       	

	/*	updateT = (tempT * tempT * tempT * tempT + dt * epsilon * Crat * Sigma_a * tempEr) / (1.0 + dt * epsilon * Crat * Sigma_a); 
		updateT = 1.0 - 0.1667 * (4.0 * exp(-0.25 * DimT) * 0.596 + 0.771) - 0.16667 * exp(-21.0 * DimT) * 4.0 * 0.301 
			-(0.464 - 0.34 * DimT) + exp(-DimT) * (0.464 - 4.3 * DimT);
	
		if(updateT < 0.0) updateT=0.0;
	*/
      		pGrid->U[ks][js][is-i].d  = 1.0;
	      	pGrid->U[ks][js][is-i].E  =  pow(updateT,0.25)/(Gamma - 1.0);
		
    }



  
}


void radMHD_inflow2(GridS *pGrid)
{
        int i, ie;
        int ks, js;
        ie = pGrid->ie;
        ks = pGrid->ks;
        js = pGrid->js;

        double dt=pGrid->dt;
	double epsilon = 0.1;
	double tempEr, tempT, updateT;

	for (i=1;  i<=nghost;  i++) {

		tempEr = pGrid->U[ks][js][ie+i].Er; 
		tempT = pGrid->U[ks][js][ie+i].E * (Gamma - 1.0);
/*		updateT = (tempT * tempT * tempT * tempT + dt * epsilon * Crat * Sigma_a * tempEr) / (1.0 + dt * epsilon * Crat * Sigma_a); */
		updateT = 0.0;
	      	pGrid->U[ks][js][ie+i].d  = 1.0;
      		pGrid->U[ks][js][ie+i].E  =  pow(updateT,0.25)/(Gamma - 1.0);
	
    }

}


void radMHD_rad_inflow(GridS *pGrid)
{
  	int i, is;
	int ks, js;
	is = pGrid->is;
  	ks = pGrid->ks;
	js = pGrid->js;

	double t=pGrid->time;
	double epsilon = 0.1;
	double DimT;
	double Frin = 0.25;
	double tempEr;

	DimT = epsilon * Crat * Sigma_a * t;
	if(DimT < 0.001){
		for(i=1; i<= nghost; i++){
			pGrid->U[ks][js][is-i].Er  = (sqrt(3.0 * DimT/(3.1415926 * epsilon)) - 3.0 * DimT / (4.0 * epsilon) 
							+ DimT * sqrt(DimT/(3.0 * 3.1415926 * epsilon))/(2.0 * epsilon)) * 4.0 * Frin;
			pGrid->U[ks][js][is-i].Fr1 = (1.0 - pGrid->U[ks][js][is-i].Er) * 0.5;
		}
	}	
	else {

		/*	tempEr = 0.615 - 2.0 * 1.74 * 0.092 * exp(-DimT)*(0.143*exp(-DimT/(0.1*0.95)) + 4.0 * 0.542 * exp(-DimT/(0.1*0.5))); 
			tempEr = 1.0 - 0.1667 * (4.0 * exp(-0.25 * DimT) * 0.596 + 0.771) - 0.16667 * exp(-21.0 * DimT) * 4.0 * 0.301;
		*/

		tempEr = 1.0 - qromo(func1, 0.0, 1.0, midpnt,DimT) - qromo(func2, 0.0, 1.0, midpnt,DimT);

		for(i=1; i<= nghost; i++){
			pGrid->U[ks][js][is-i].Er = tempEr * 4.0 * Frin;
			pGrid->U[ks][js][is-i].Fr1 = (1.0 - tempEr) * 0.5;

		}
	}

		for(i=1; i<= nghost; i++){
                        pGrid->U[ks][js][is-i].Edd_11 = 0.3333;                

                }

  
}

void radMHD_rad_inflow2(GridS *pGrid)
{
        int i, ie;
        int ks, js;
        ie = pGrid->ie;
        ks = pGrid->ks;
        js = pGrid->js;


    for (i=1;  i<=nghost;  i++) {
      pGrid->U[ks][js][ie+i].Er  = 0.0;
      /*pGrid->U[ks][js][ie+i].Fr1 = Crat * 0.3333 * pGrid->U[ks][js][ie].Er /(Sigma_t * pGrid->dx1); */
     pGrid->U[ks][js][ie+i].Fr1 = 0.0;
	 pGrid->U[ks][js][ie+i].Edd_11 = 0.333333;
	 }
}


/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
	fprintf(fp,"%5.3e\n",Gamma);
	fprintf(fp,"%5.3e\n",Prat);
	fprintf(fp,"%5.3e\n",Crat);
	fprintf(fp,"%5.3e\n",Sigma_t);
	fprintf(fp,"%5.3e\n",Sigma_a);
	fprintf(fp,"%5.3e\n",R_ideal);
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{

	bvals_mhd_fun(&(pM->Domain[0][0]), right_x1, radMHD_inflow);
	fscanf(fp,"%lf",&Gamma);
	fscanf(fp,"%lf\n",&Prat);
	fscanf(fp,"%lf\n",&Crat);
	fscanf(fp,"%lf\n",&Sigma_t);
	fscanf(fp,"%lf\n",&Sigma_a);
	fscanf(fp,"%lf\n",&R_ideal);


  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  return;
}

void Userwork_after_loop(MeshS *pM)
{
  return;
}



/* For integrate improper function */


#define FUNC(x,y) ((*func)(x,y))

double midpnt(double (*func)(double,double), double a, double b, int n, double t)
{
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b),t));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x,t);
			x += ddel;
			sum += FUNC(x,t);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC

#define EPS 1.0e-6
#define JMAX 14
#define JMAXP (JMAX+1)
#define K 5


double qromo(double (*func)(double,double), double a, double b,
	double (*choose)(double(*)(double,double), double, double, int,double),double t) 
{
	void polint(double xa[], double ya[], int n, double x, double *y, double *dy);
	
	int j;
	double ss,dss,h[JMAXP+1],s[JMAXP];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j,t);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	ath_error("Too many steps in routing qromo");
	return 0.0;
}

#undef EPS
#undef JMAX
#undef JMAXP
#undef K



void polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=(double *)malloc((n+1)*sizeof(double));
	d=(double *)malloc((n+1)*sizeof(double));
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) ath_error("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free(d);
	free(c);
}

double func1(double x, double t)
{

	double epsi,gamma1,theta1,result;
	epsi = 0.1;
	gamma1 = x * sqrt(epsi + 1.0/(1.0-x * x));
	theta1=sqrt(1.0 - 3.0/(3.0 + 4.0 * x * x * (0.1 + 1.0/(1.0 - x * x))));
	result = 2.0 * 1.732 * 0.318 * exp(-t * x * x) * theta1/(x * sqrt(3.0 + 4.0 * gamma1 * gamma1));
	return result;
}

double func2(double x, double t)
{

	double epsi,gamma2,theta2,result;
	epsi = 0.1;
	gamma2 = sqrt((1.0 - x) * (epsi + 1.0/x));
	theta2=sqrt(1.0 - 3.0/(3.0 + 4.0 * (0.1 + 1.0/x) * (1.0 - x)));
	result = 1.732 * 0.318 * exp(-t) * exp(-t/(epsi * x)) * theta2/(x * (1.0+epsi * x)* sqrt(3.0 + 4.0 * gamma2 * gamma2));
	return result;
}


double func3(double x, double t)
{

	double epsi,gamma3,theta3,result;
	epsi = 0.1;
	gamma3 = sqrt((1.0 - x * x) * (epsi + 1.0/(x * x)));
	theta3=sqrt(1.0 - 3.0/(3.0 + 4.0 * (0.1 + 1.0/(x*x)) * (1.0 - x*x)));
	result = 2.0 * 1.732 * 0.318 * exp(-t * (1.0 -x * x)) * theta3/sqrt(4.0 - x * x + 4.0 * epsi * x * x * (1.0 - x * x));
	return result;
}

double func4(double x, double t)
{

	double epsi,gamma2,theta2,result;
	epsi = 0.1;
	gamma2 = sqrt((1.0 - x) * (epsi + 1.0/x));
	theta2=sqrt(1.0 - 3.0/(3.0 + 4.0 * (0.1 + 1.0/x) * (1.0 - x)));
	result = 1.732 * 0.318 * exp(-t) * exp(-t/(epsi * x)) * theta2/(x * sqrt(3.0 + 4.0 * gamma2 * gamma2));
	return result;
}


