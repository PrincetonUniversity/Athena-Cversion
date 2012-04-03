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
 *   BackEuler_3d()
 * PRIVATE FUNCTION
 *	INI_matrix()
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


#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
/*================================*/
/* For the matrix solver */
/* we use lis library now */
#include <HYPRE_utilities.h>
#include <HYPRE_krylov.h>
#include <HYPRE.h>
#include <HYPRE_parcsr_ls.h>

/*===============================*/

#endif

#if defined(RADIATIONMHD_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The radiation MHD integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */



/*===========================================================
 * Private functions used for different boundary conditions *
 * and MPI case *
 */

void INI_matrix(MeshS *pM);



/* No boundary case */
void i_j_k(const int i, const int j, const int k);

/* one boundary case */
void i_j_ke_phy(const int i, const int j);
void i_j_ke_MPI(const int i, const int j);
void i_j_ks_phy(const int i, const int j);
void i_j_ks_MPI(const int i, const int j);


void i_je_k_phy(const int i, const int k);
void i_je_k_MPI(const int i, const int k);
void i_js_k_phy(const int i, const int k);
void i_js_k_MPI(const int i, const int k);

void ie_j_k_phy(const int j, const int k);
void ie_j_k_MPI(const int j, const int k);
void is_j_k_phy(const int j, const int k);
void is_j_k_MPI(const int j, const int k);

/* two boundaries case */ 

/* for i */
void i_je_ke_phy_phy(const int i);
void i_je_ke_MPI_phy(const int i);
void i_je_ke_phy_MPI(const int i);
void i_je_ke_MPI_MPI(const int i);

void i_je_ks_phy_phy(const int i);
void i_je_ks_MPI_phy(const int i);
void i_je_ks_phy_MPI(const int i);
void i_je_ks_MPI_MPI(const int i);

void i_js_ke_phy_phy(const int i);
void i_js_ke_MPI_phy(const int i);
void i_js_ke_phy_MPI(const int i);
void i_js_ke_MPI_MPI(const int i);

void i_js_ks_phy_phy(const int i);
void i_js_ks_MPI_phy(const int i);
void i_js_ks_phy_MPI(const int i);
void i_js_ks_MPI_MPI(const int i);

/* for j */

void ie_j_ke_phy_phy(const int j);
void ie_j_ke_MPI_phy(const int j);
void ie_j_ke_phy_MPI(const int j);
void ie_j_ke_MPI_MPI(const int j);

void ie_j_ks_phy_phy(const int j);
void ie_j_ks_MPI_phy(const int j);
void ie_j_ks_phy_MPI(const int j);
void ie_j_ks_MPI_MPI(const int j);

void is_j_ke_phy_phy(const int j);
void is_j_ke_MPI_phy(const int j);
void is_j_ke_phy_MPI(const int j);
void is_j_ke_MPI_MPI(const int j);

void is_j_ks_phy_phy(const int j);
void is_j_ks_MPI_phy(const int j);
void is_j_ks_phy_MPI(const int j);
void is_j_ks_MPI_MPI(const int j);

/* for k */


void ie_je_k_phy_phy(const int k);
void ie_je_k_MPI_phy(const int k);
void ie_je_k_phy_MPI(const int k);
void ie_je_k_MPI_MPI(const int k);

void ie_js_k_phy_phy(const int k);
void ie_js_k_MPI_phy(const int k);
void ie_js_k_phy_MPI(const int k);
void ie_js_k_MPI_MPI(const int k);

void is_je_k_phy_phy(const int k);
void is_je_k_MPI_phy(const int k);
void is_je_k_phy_MPI(const int k);
void is_je_k_MPI_MPI(const int k);

void is_js_k_phy_phy(const int k);
void is_js_k_MPI_phy(const int k);
void is_js_k_phy_MPI(const int k);
void is_js_k_MPI_MPI(const int k);


/* Three boundaries case */

void ie_je_ke_phy_phy_phy();
void ie_je_ke_MPI_phy_phy();
void ie_je_ke_phy_MPI_phy();
void ie_je_ke_MPI_MPI_phy();
void ie_je_ke_phy_phy_MPI();
void ie_je_ke_MPI_phy_MPI();
void ie_je_ke_phy_MPI_MPI();
void ie_je_ke_MPI_MPI_MPI();


void is_je_ke_phy_phy_phy();
void is_je_ke_MPI_phy_phy();
void is_je_ke_phy_MPI_phy();
void is_je_ke_MPI_MPI_phy();
void is_je_ke_phy_phy_MPI();
void is_je_ke_MPI_phy_MPI();
void is_je_ke_phy_MPI_MPI();
void is_je_ke_MPI_MPI_MPI();

void ie_js_ke_phy_phy_phy();
void ie_js_ke_MPI_phy_phy();
void ie_js_ke_phy_MPI_phy();
void ie_js_ke_MPI_MPI_phy();
void ie_js_ke_phy_phy_MPI();
void ie_js_ke_MPI_phy_MPI();
void ie_js_ke_phy_MPI_MPI();
void ie_js_ke_MPI_MPI_MPI();


void is_js_ke_phy_phy_phy();
void is_js_ke_MPI_phy_phy();
void is_js_ke_phy_MPI_phy();
void is_js_ke_MPI_MPI_phy();
void is_js_ke_phy_phy_MPI();
void is_js_ke_MPI_phy_MPI();
void is_js_ke_phy_MPI_MPI();
void is_js_ke_MPI_MPI_MPI();




void ie_je_ks_phy_phy_phy();
void ie_je_ks_MPI_phy_phy();
void ie_je_ks_phy_MPI_phy();
void ie_je_ks_MPI_MPI_phy();
void ie_je_ks_phy_phy_MPI();
void ie_je_ks_MPI_phy_MPI();
void ie_je_ks_phy_MPI_MPI();
void ie_je_ks_MPI_MPI_MPI();


void is_je_ks_phy_phy_phy();
void is_je_ks_MPI_phy_phy();
void is_je_ks_phy_MPI_phy();
void is_je_ks_MPI_MPI_phy();
void is_je_ks_phy_phy_MPI();
void is_je_ks_MPI_phy_MPI();
void is_je_ks_phy_MPI_MPI();
void is_je_ks_MPI_MPI_MPI();

void ie_js_ks_phy_phy_phy();
void ie_js_ks_MPI_phy_phy();
void ie_js_ks_phy_MPI_phy();
void ie_js_ks_MPI_MPI_phy();
void ie_js_ks_phy_phy_MPI();
void ie_js_ks_MPI_phy_MPI();
void ie_js_ks_phy_MPI_MPI();
void ie_js_ks_MPI_MPI_MPI();


void is_js_ks_phy_phy_phy();
void is_js_ks_MPI_phy_phy();
void is_js_ks_phy_MPI_phy();
void is_js_ks_MPI_MPI_phy();
void is_js_ks_phy_phy_MPI();
void is_js_ks_MPI_phy_MPI();
void is_js_ks_phy_MPI_MPI();
void is_js_ks_MPI_MPI_MPI();




/*===========================================================*/



static HYPRE_IJMatrix Euler;
static HYPRE_ParCSRMatrix parcsr_Euler;
static HYPRE_IJVector RHSEuler2;
static HYPRE_ParVector par_RHSEuler2;
static HYPRE_IJVector INIguess2;
static HYPRE_ParVector par_INIguess2;

static HYPRE_Solver solver, precond;

static int solverflag = 1;
/* flag to choose different solver */
/* flag = 1 GMRES with BoomerAMG preconditioner , multigrid */
/* flag = 2 GMRES with ilu preconditioner, or hydro */
/* flag = 3 BoomerAMG solver */

static double *Value;
static int *indexValue;
static int *ptr;
static int *ncols;
static int *rows;
static double *RHSEuler;
static double *INIguess;
static int *NZrow;



/* parameters used to setup the matrix elements */
static int ilower;
static int iupper;
static int lines;
static int NoEr;
static int NoFr1;
static int NoFr2;
static int NoFr3;
static int is;
static int ie;
static int js;
static int je;
static int ks;
static int ke;
static Real *theta;
static Real *phi;
static Real *psi;
static Real *varphi;


/* subtract background state, which is initial condition */
/* For this to work, we assume no  initial perturbation in Er and Fr, initial v=0 */


/* boundary flag */
static int ix1;
static int ox1;
static int ix2;
static int ox2;
static int ix3;
static int ox3;

static int NGx; /* Number of Grid in x direction */
static int NGy; /* Number of Grid in y direction */
static int NGz; /* Number of Grid in z direction */
static int Nx; /* Number of zones in x direction of each Grid */
static int Ny; /* Number of Zones in y direction of each Grid */
static int Nz; /* Number of Zones in z direction of each Grid */
static int count_Grids;
static int ID;
static int lx1;
static int rx1;
static int lx2;
static int rx2;
static int lx3;
static int rx3;


/* index variables to setup the blocks */
/* those indices are updated in the top functions */
/* Those indices are used to record the order in each row */
static int indexk0Er; static int indexk0Fr1; static int indexk0Fr2; static int indexk0Fr3;
static int indexj0Er; static int indexj0Fr1; static int indexj0Fr2; static int indexj0Fr3;
static int indexi0Er; static int indexi0Fr1; static int indexi0Fr2; static int indexi0Fr3;
static int indexiEr; static int indexiFr1; static int indexiFr2; static int indexiFr3;
static int indexi1Er; static int indexi1Fr1; static int indexi1Fr2; static int indexi1Fr3;
static int indexj1Er; static int indexj1Fr1; static int indexj1Fr2; static int indexj1Fr3;
static int indexk1Er; static int indexk1Fr1; static int indexk1Fr2; static int indexk1Fr3;

/* Used to save the index of each column so that we do not need to calculate it again */
static int **colindex;

/* the followings are used for column number */
static int colk0;
static int colj0;
static int coli0;
static int coli;
static int coli1;
static int colj1;
static int colk1;



/* for shearing box boundary condition */
#ifdef SHEARING_BOX
static int joffset; /* total sheared j cells */
static int jshearing; /* j index in the sheared position */
static int shearing_grid; /* ID of the shifted grid */
static int col_shear;

static int iproc;
static int jproc;
static int kproc;

static int sheark0;
static int shearj0;
static int sheari0;
static int sheari;
static int sheari1;
static int shearj1;
static int sheark1;

static void shearing_matrix(MeshS *pM);

#endif


static int MPIcount1;
static int MPIcount2; 
static int MPIcount2F;
/* Used for MPI periodic boundary condition */



/********Public function****************/
/*-------BackEuler_3d(): Use back euler method to update E_r and Fluxr-----------*/
/* Not work with SMR now. So there is only one Domain in the Mesh */


void BackEuler_3d(MeshS *pM)
{
/* Cell number in root domain is Nx[0,1,2] */
	

	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	Real hdtodx1 = 0.5 * pG->dt/pG->dx1;
	Real hdtodx2 = 0.5 * pG->dt/pG->dx2;
	Real hdtodx3 = 0.5 * pG->dt/pG->dx3;
	Real dt = pG->dt;
	
	int i, j, k;
	is = pG->is;
	ie = pG->ie;
	js = pG->js;
	je = pG->je;
	ks = pG->ks;
	ke = pG->ke;
	
	
	/* ID of this grid */
#ifdef MPI_PARALLEL
	ID = myID_Comm_world;
	lx1 = pG->lx1_id;
	rx1 = pG->rx1_id;
	lx2 = pG->lx2_id;
	rx2 = pG->rx2_id;
	lx3 = pG->lx3_id;
	rx3 = pG->rx3_id;
#else
	ID = 0;
	lx1 = -1;
	rx1 = -1;
	lx2 = -1;
	rx2 = -1;
	lx3 = -1;
	rx3 = -1;
#endif	

	Nx = ie - is + 1;
	Ny = je - js + 1;
	Nz = ke - ks + 1;
	
	lines  = 4 * Nx * Ny * Nz; 
	
	count_Grids = ID * lines;
	
	ilower = count_Grids;
	iupper = ilower + lines - 1;
	
		/* count_Grids is the total number of lines before this Grids, which depends on the relative position of this grid */	
	
	/* Boundary condition flag */
	
	ix1 = pM->BCFlag_ix1;
	ox1 = pM->BCFlag_ox1;
	ix2 = pM->BCFlag_ix2;
	ox2 = pM->BCFlag_ox2;
	ix3 = pM->BCFlag_ix3;
	ox3 = pM->BCFlag_ox3;
	
	

	Real vshear = 0.0;
	Real qom = 0.0;
	/* For shearing box boundary condition */
#ifdef SHEARING_BOX
	Real xmin, xmax, Lx, Ly, yshear, deltay, qomL, dFrycoef;
	xmin = pD->RootMinX[0];
	xmax = pD->RootMaxX[0];
	Lx = xmax - xmin;
	
	xmin = pD->RootMinX[1];
	xmax = pD->RootMaxX[1];
	Ly = xmax - xmin;
	
	qomL = qshear*Omega_0*Lx;
	yshear = qomL*pG->time;
	
	
	/* Only need to apply shearing boundary condition in x direction for flux y */
	/* Others are periodic boundary condition */
	/* This is now set in each cell to include the factor (1+f_yy) */	
/*	deltay = fmod(yshear, Ly);
	
	joffset = (int)(deltay/pG->dx2);
	
	if((deltay - joffset * pG->dx2) > 0.5 * pG->dx2)
		joffset += 1;
	
	if(joffset >= NGy * Ny)
		joffset -= NGy * Ny;
*/	
	
	
#ifdef MPI_PARALLEL

	get_myGridIndex(pD, myID_Comm_world, &iproc, &jproc, &kproc);
	
	/* Total number of lines from level k=0 to k=kGrids -1 (if kGrids >= 1) */
#else
	
	iproc = 0;
	jproc = 0;
	kproc = 0;

	
#endif /* End MPI_ parallel */
	
	
	/* check the boundary flag is correct */
	/* Boundary condition along z direction is not necessary periodic */
	if((ix1 != 4) || (ix2 !=4)  || (ox1 != 4) || (ox2 != 4) )
		ath_error("[BackEuler_3d]: Shearing box must have boundary flag 4 along three directions!\n");
	

	/* For FARGO algorithm */
#ifdef FARGO
	
	Real x1, x2, x3;
	qom = qshear * Omega_0;
	
#endif
	
#endif /* End shearing box */
	
	

	Real temperature, velocity_x, velocity_y, velocity_z, pressure, T4, Fr0x, Fr0y, Fr0z, density;
	Real AdvFx, AdvFy, AdvFz;
	Real tempEr3, tempEr2, tempEr1;
	Real tempFr3, tempFr2, tempFr1;
	Real temp0;
	/* This is equivilent to Cspeeds[] in 1D */
#ifdef RADIATION_TRANSFER
	Real fluxi0, fluxi1, fluxj0, fluxj1, fluxk0, fluxk1;

#endif


  	Real Sigma_aF, Sigma_aP, Sigma_aE, Sigma_sF;
	int m;
	Real Sigma[NOPACITY];


/* variables used to subtract background state */
/* NOTICE Here we assume background state grad E_r = -Sigma_t Fr, a uniform flux and background velocity is zero */
/* We assume the background state, Eddington tensor is 1/3 */
/* Initial perturbations are not added in Er and Fr so that we can use variables at t0 as background state */

	static int firstflag = 1; 
	static int shearflag = 0;
	int bgflag;		/* used to subtract whether subtract background or not */
	
	bgflag = 0;	/* 1 means subtract background, 0 means not. Default is not */
	
	

	/* For temporary use only */
	int index,Matrixiter;
	Real tempvalue;

	
	/* *****************************************************/
/* Step 1 : Use Backward Euler to update the radiation energy density and flux */


	/* First, calculate the advection flux  */
	/* and subtract the advection flux */



/* Step 1a: Calculate the Matrix elements  */
/* ie-is+1 =size1, otherwise it is wrong */


/*--------------------Note--------------------*/


/* Step 1b: Setup the Matrix */
	
	
	/*--------Now set the Euler matrix-----------*/

	
	/* Set the positions of the elements. Only do this once */
	if(firstflag)
	{

		
		INI_matrix(pM);
		
		/* If the code starts with time nonzero, start from restart, we still need to *
		 * apply shearing boundary condition */
#ifdef SHEARING_BOX
		shearing_matrix(pM);
#endif		
		firstflag = 0;
	}
	
	/* For shearing-box, we should update the position of the elements for boundary cells */
	
#ifdef SHEARING_BOX
	if(shearflag){
		shearing_matrix(pM);
	
	
	/* For shearing_box, we need to destroy the matrix completely *
	 * because the patter of the matrix is changed */
		
		/* Create the matrix */
		HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &Euler);
		HYPRE_IJMatrixSetObjectType(Euler, HYPRE_PARCSR);
		
		/* set the maximum number of elements in each row */
		HYPRE_IJMatrixSetRowSizes(Euler, NZrow);
	
	
	}
	else{
		shearflag = 1;
	}
#endif

/* First, do the advection step and update bounary */

	
	
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){
			for(i=is; i<=ie; i++){
			vshear = 0.0;
#ifdef SHEARING_BOX

		/*	yshear *= (1.0 + pG->U[k][j][i].Edd_22);
		*/
			deltay = fmod(yshear, Ly);
	
			joffset = (int)(deltay/pG->dx2);
	
			if((deltay - joffset * pG->dx2) > 0.5 * pG->dx2)
				joffset += 1;
	
			if(joffset >= NGy * Ny)
				joffset -= NGy * Ny;
			
#ifdef FARGO				
			/* With FARGO, we should add background shearing to the source terms */
			cc_pos(pG,i,j,k,&x1,&x2,&x3);
			vshear = qom * x1;			
				
#endif				
#endif	

			/*==============================================*/
			/* if background state is subtracted, initial guess is zero */
			/* first, setup the initial guess */

			if(bgflag){
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is)] = 0.0;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 1] = 0.0;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 2] = 0.0;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 3] = 0.0;
			}
			else{
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is)] = pG->U[k][j][i].Er;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 1] = pG->U[k][j][i].Fr1;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 2] = pG->U[k][j][i].Fr2;
				INIguess[4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is) + 3] = pG->U[k][j][i].Fr3;

			}



			/*================================================*/





			matrix_coef(NULL, pG, 3, i, j, k, qom, theta, phi, psi, varphi);


			/* Set the right hand side */
			T4 = pG->Tguess[k][j][i];

		

		/* RHSEuler[0...N-1]  */
		Sigma_sF = pG->U[k][j][i].Sigma[0];
		Sigma_aF = pG->U[k][j][i].Sigma[1];
		Sigma_aP = pG->U[k][j][i].Sigma[2];
		Sigma_aE = pG->U[k][j][i].Sigma[3];

		velocity_x = pG->U[k][j][i].M1 /pG->U[k][j][i].d;
		velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
		velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;
				
#ifdef FARGO

		velocity_y -= vshear;			
				
#endif				

						  				
		/*-----------------------------*/	
		/* Tguess is now the energy source term, which should be added */
		/* calculate the advection term */
		Rad_Advection_Flux3D(pD, i, j, k, 1.0, &AdvFx, &AdvFy, &AdvFz);
    		tempvalue   = pG->U[k][j][i].Er + dt * Sigma_aP * T4 * Crat * Eratio + (1.0 - Eratio) * pG->Ersource[k][j][i] + (AdvFx + AdvFy + AdvFz) + pG->Comp[k][j][i];
			
		if(bgflag){
			tempEr3 = theta[0] * pG->U[k-1][j][i].Er + theta[14] * pG->U[k+1][j][i].Er;
			tempEr2 = theta[2] * pG->U[k][j-1][i].Er + theta[12] * pG->U[k][j+1][i].Er;
			tempEr1 = theta[4] * pG->U[k][j][i-1].Er + theta[10] * pG->U[k][j][i+1].Er;

			tempFr3 = theta[1] * pG->U[k-1][j][i].Fr3 + theta[15] * pG->U[k+1][j][i].Fr3;
			tempFr2 = theta[3] * pG->U[k][j-1][i].Fr2 + theta[13] * pG->U[k][j+1][i].Fr2;
			tempFr1 = theta[5] * pG->U[k][j][i-1].Fr1 + theta[11] * pG->U[k][j][i+1].Fr1;

			temp0 = theta[7] * pG->U[k][j][i].Fr1 + theta[8] * pG->U[k][j][i].Fr2 + theta[9] * pG->U[k][j][i].Fr3;

			tempvalue -= (theta[6] * pG->U[k][j][i].Er + (tempEr1 + tempEr2 + tempEr3));
			tempvalue -= ((tempFr1 + tempFr2 + tempFr3) + temp0);
		
		}



		index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is);
		RHSEuler[index] = tempvalue;

		/*----------------------------*/
    		tempvalue = pG->U[k][j][i].Fr1 + Eratio * dt * Sigma_aP * T4 * velocity_x + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_x / Crat;
		
		if(bgflag){

			tempEr3 = phi[0] * pG->U[k-1][j][i].Er + phi[12] * pG->U[k+1][j][i].Er;
			tempEr2 = phi[2] * pG->U[k][j-1][i].Er + phi[10] * pG->U[k][j+1][i].Er;
			tempEr1 = phi[4] * pG->U[k][j][i-1].Er + phi[8] * pG->U[k][j][i+1].Er;

			tempFr3 = phi[1] * pG->U[k-1][j][i].Fr1 + phi[13] * pG->U[k+1][j][i].Fr1;
			tempFr2 = phi[3] * pG->U[k][j-1][i].Fr1 + phi[11] * pG->U[k][j+1][i].Fr1;
			tempFr1 = phi[5] * pG->U[k][j][i-1].Fr1 + phi[9] * pG->U[k][j][i+1].Fr1;

			temp0 = phi[6] * pG->U[k][j][i].Er;
			tempvalue -= (phi[7] * pG->U[k][j][i].Fr1 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);


		}
	

		++index;
		RHSEuler[index] = tempvalue;
		
		/*-------------------------*/
		tempvalue = pG->U[k][j][i].Fr2 + Eratio * dt * Sigma_aP * T4 * velocity_y + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_y / Crat;
		
		if(bgflag){


			tempEr3 = psi[0] * pG->U[k-1][j][i].Er + psi[12] * pG->U[k+1][j][i].Er;
			tempEr2 = psi[2] * pG->U[k][j-1][i].Er + psi[10] * pG->U[k][j+1][i].Er;
			tempEr1 = psi[4] * pG->U[k][j][i-1].Er + psi[8] * pG->U[k][j][i+1].Er;

			tempFr3 = psi[1] * pG->U[k-1][j][i].Fr2 + psi[13] * pG->U[k+1][j][i].Fr2;
			tempFr2 = psi[3] * pG->U[k][j-1][i].Fr2 + psi[11] * pG->U[k][j+1][i].Fr2;
			tempFr1 = psi[5] * pG->U[k][j][i-1].Fr2 + psi[9] * pG->U[k][j][i+1].Fr2;

			temp0 = psi[6] * pG->U[k][j][i].Er;

			
			tempvalue -= (psi[7] * pG->U[k][j][i].Fr2 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);


		}

		++index;
		RHSEuler[index] = tempvalue;

		/*-------------------------*/
		tempvalue = pG->U[k][j][i].Fr3 + Eratio * dt * Sigma_aP * T4 * velocity_z + (1.0 - Eratio) * pG->Ersource[k][j][i] * velocity_z / Crat;

		if(bgflag){

			tempEr3 = varphi[0] * pG->U[k-1][j][i].Er + varphi[12] * pG->U[k+1][j][i].Er;
			tempEr2 = varphi[2] * pG->U[k][j-1][i].Er + varphi[10] * pG->U[k][j+1][i].Er;
			tempEr1 = varphi[4] * pG->U[k][j][i-1].Er + varphi[8] * pG->U[k][j][i+1].Er;

			tempFr3 = varphi[1] * pG->U[k-1][j][i].Fr3 + varphi[13] * pG->U[k+1][j][i].Fr3;
			tempFr2 = varphi[3] * pG->U[k][j-1][i].Fr3 + varphi[11] * pG->U[k][j+1][i].Fr3;
			tempFr1 = varphi[5] * pG->U[k][j][i-1].Fr3 + varphi[9] * pG->U[k][j][i+1].Fr3;

			temp0 = varphi[6] * pG->U[k][j][i].Er;

			
			tempvalue -= (varphi[7] * pG->U[k][j][i].Fr3 + (tempEr1 + tempEr2 + tempEr3) + (tempFr1 + tempFr2 + tempFr3) + temp0);

		}		


		++index;
		RHSEuler[index] = tempvalue;

		/* When background is substracted, boundary condition is included in the right hand side */
		if(!bgflag){
		/* For inflow boundary condition along x direction*/
		if((i == is) && (ix1 == 3) && (lx1 < 0)) {
			

			/* Subtract some value */
			tempvalue = -(theta[4] * pG->U[k][j][i-1].Er + theta[5] * pG->U[k][j][i-1].Fr1);
			index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is);
			RHSEuler[index] += tempvalue;



			tempvalue = -(phi[4] * pG->U[k][j][i-1].Er + phi[5] * pG->U[k][j][i-1].Fr1);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(psi[4] * pG->U[k][j][i-1].Er + psi[5] * pG->U[k][j][i-1].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[4] * pG->U[k][j][i-1].Er + varphi[5] * pG->U[k][j][i-1].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
			
		}

		if((i == ie) && (ox1 == 3) && (rx1 < 0)) {
			
			
			

			tempvalue = -(theta[10] * pG->U[k][j][i+1].Er + theta[11] * pG->U[k][j][i+1].Fr1);
			index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is); 
			RHSEuler[index] += tempvalue;

			tempvalue = -(phi[8] * pG->U[k][j][i+1].Er + phi[9] * pG->U[k][j][i+1].Fr1);
			++index; 
			RHSEuler[index] += tempvalue;

			tempvalue = -(psi[8] * pG->U[k][j][i+1].Er + psi[9] * pG->U[k][j][i+1].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[8] * pG->U[k][j][i+1].Er + varphi[9] * pG->U[k][j][i+1].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
					
		}
		

		/* For inflow boundary condition along y direction*/
		if((j == js) && (ix2 == 3) && (lx2 < 0)) {
			

			tempvalue = -(theta[2] * pG->U[k][j-1][i].Er + theta[3] * pG->U[k][j-1][i].Fr2);
			index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4 * (i - is);
			RHSEuler[index] += tempvalue;


			tempvalue = -(phi[2] * pG->U[k][j-1][i].Er + phi[3] * pG->U[k][j-1][i].Fr1);
			++index;
			RHSEuler[index] += tempvalue;


			tempvalue = -(psi[2] * pG->U[k][j-1][i].Er + psi[3] * pG->U[k][j-1][i].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[2] * pG->U[k][j-1][i].Er + varphi[3] * pG->U[k][j-1][i].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
				
		}

		if((j == je) && (ox2 == 3) && (rx2 < 0)) {
			
			

			tempvalue = -(theta[12] * pG->U[k][j+1][i].Er + theta[13] * pG->U[k][j+1][i].Fr2);
			index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is);
			RHSEuler[index] += tempvalue;

			tempvalue = -(phi[10] * pG->U[k][j+1][i].Er + phi[11] * pG->U[k][j+1][i].Fr1);
			++index;
			RHSEuler[index] += tempvalue;


			tempvalue = -(psi[10] * pG->U[k][j+1][i].Er + psi[11] * pG->U[k][j+1][i].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[10] * pG->U[k][j+1][i].Er + varphi[11] * pG->U[k][j+1][i].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
			
		}

		/* inflow boundary condition in z direction */
		if((k == ks) && (ix3 == 3) && (lx3 < 0)) {
			


			tempvalue = -(theta[0] * pG->U[k-1][j][i].Er + theta[1] * pG->U[k-1][j][i].Fr3);
			index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4 * (i - is);
			RHSEuler[index] += tempvalue;


			tempvalue = -(phi[0] * pG->U[k-1][j][i].Er + phi[1] * pG->U[k-1][j][i].Fr1);
			++index;
			RHSEuler[index] += tempvalue;


			tempvalue = -(psi[0] * pG->U[k-1][j][i].Er + psi[1] * pG->U[k-1][j][i].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[0] * pG->U[k-1][j][i].Er + varphi[1] * pG->U[k-1][j][i].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
				
		}

		if((k == ke) && (ox3 == 3) && (rx3 < 0)) {
			/* k + 1 */
			
			

			tempvalue = -(theta[14] * pG->U[k+1][j][i].Er + theta[15] * pG->U[k+1][j][i].Fr3);
			index =  4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4 * (i - is);
			RHSEuler[index] += tempvalue;

			tempvalue = -(phi[12] * pG->U[k+1][j][i].Er + phi[13] * pG->U[k+1][j][i].Fr1);
			++index;
			RHSEuler[index] += tempvalue;


			tempvalue = -(psi[12] * pG->U[k+1][j][i].Er + psi[13] * pG->U[k+1][j][i].Fr2);
			++index;
			RHSEuler[index] += tempvalue;

			tempvalue = -(varphi[12] * pG->U[k+1][j][i].Er + varphi[13] * pG->U[k+1][j][i].Fr3);
			++index;
			RHSEuler[index] += tempvalue;
			
		}

		} /* End bgflag */		
				/* Need to modify the matrix elements for shearing box boundary conditoin */
				/* Eddington tensor in the ghost zone is already set by boundary condition */
				/* Eddington tensor is assumed to be a constant tduring this step. It is explicit */
					
#ifdef SHEARING_BOX
				if((i == is) && ((pG->lx1_id < 0) || (pG->lx1_id > ID))){
					dFrycoef = (1.0 + pG->U[k][j][i-1].Edd_22) * qomL / Crat;
					psi[4] += psi[5] * dFrycoef;
					
				}else if((i == ie) && ((pG->rx1_id < 0) || (pG->rx1_id < ID))){
					dFrycoef = (1.0 + pG->U[k][j][i+1].Edd_22) * qomL / Crat;
					psi[8] -= psi[9] * dFrycoef;				
					
				}			
#endif
				

				
		 /* get the index */
				indexk0Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][0]; 
				indexk0Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][1]; 
				indexk0Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][2];
				indexk0Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][3];
				indexj0Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][4]; 
				indexj0Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][5]; 
				indexj0Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][6];  
				indexj0Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][7]; 
				indexi0Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][8]; 
				indexi0Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][9];  
				indexi0Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][10]; 
				indexi0Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][11]; 
				indexiEr	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][12];   
				indexiFr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][13];   
				indexiFr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][14];   
				indexiFr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][15]; 
				indexi1Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][16];  
				indexi1Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][17];  
				indexi1Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][18];  
				indexi1Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][19]; 
				indexj1Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][20];  
				indexj1Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][21];  
				indexj1Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][22];  
				indexj1Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][23]; 
				indexk1Er	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][24];  
				indexk1Fr1	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][25]; 
				indexk1Fr2	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][26];  
				indexk1Fr3	= colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][27]; 
				
		/* Need to modify the  matrix elements to account for different boundary 
		 * conditions. This needs to be done for every step before the matrix elements 
		 * are set.
		 */				
				
				


		/* Now indices are updated. We update the vectors */
			/* There are 7 blocks */
			if(indexk0Er >= 0){
				/* for Er */
				Value[indexk0Er] = theta[0];
				Value[indexk0Er+1] = theta[1];

				/* For Fr1 */
				Value[indexk0Fr1] = phi[0];
				Value[indexk0Fr1+1] = phi[1];
				
				/* For Fr2 */
				Value[indexk0Fr2] = psi[0];
				Value[indexk0Fr2+1] = psi[1];
	
				/* For Fr3 */
				Value[indexk0Fr3] = varphi[0];
				Value[indexk0Fr3+1] = varphi[1];
			
			}
			else{
				if(ix3 == 1 || ix3 == 5){
					
					theta[6] += theta[0];
					theta[9] -= theta[1];
					
					phi[6] += phi[0];
					phi[7] += phi[1];
					
					psi[6] += psi[0];
					psi[7] += psi[1];
					
					varphi[6] += varphi[0];
					varphi[7] -= varphi[1];
				}
				else if(ix3 == 2){
					theta[6] += theta[0];
					theta[9] += theta[1];
					
					phi[6] += phi[0];
					phi[7] += phi[1];
					
					psi[6] += psi[0];
					psi[7] += psi[1];
					
					varphi[6] += varphi[0];
					varphi[7] += varphi[1];
				}
				else if(ix3 == 3){
					/* Do nothing */
				}
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ix3\n");
				}
				
			}

			if(indexj0Er >= 0){
				/* for Er */
				Value[indexj0Er] = theta[2];
				Value[indexj0Er+1] = theta[3];

				/* For Fr1 */
				Value[indexj0Fr1] = phi[2];
				Value[indexj0Fr1+1] = phi[3];

				/* For Fr2 */
				Value[indexj0Fr2] = psi[2];
				Value[indexj0Fr2+1] = psi[3];
	
				/* For Fr3 */
				Value[indexj0Fr3] = varphi[2];
				Value[indexj0Fr3+1] = varphi[3];
	
			}
			else{
				if(ix2 == 1 || ix2 == 5){
					
					theta[6] += theta[2];
					theta[8] -= theta[3];
					
					phi[6] += phi[2];
					phi[7] += phi[3];
					
					psi[6] += psi[2];
					psi[7] -= psi[3];
					
					varphi[6] += varphi[2];
					varphi[7] += varphi[3];
				}
				else if(ix2 == 2){
					theta[6] += theta[2];
					theta[8] += theta[3];
					
					phi[6] += phi[2];
					phi[7] += phi[3];
					
					psi[6] += psi[2];
					psi[7] += psi[3];
					
					varphi[6] += varphi[2];
					varphi[7] += varphi[3];
				}
				else if(ix2 == 3){
					/* Do nothing */
				}
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ix2\n");
				}
			
			}

			if(indexi0Er >= 0){
				/* for Er */
				Value[indexi0Er] = theta[4];
				Value[indexi0Er+1] = theta[5];

				/* For Fr1 */
				Value[indexi0Fr1] = phi[4];
				Value[indexi0Fr1+1] = phi[5];
	
				/* For Fr2 */
				Value[indexi0Fr2] = psi[4];
				Value[indexi0Fr2+1] = psi[5];
	
				/* For Fr3 */
				Value[indexi0Fr3] = varphi[4];
				Value[indexi0Fr3+1] = varphi[5];

			}
			else{
				if(ix1 == 1 || ix1 == 5){
					
					theta[6] += theta[4];
					theta[7] -= theta[5];
					
					phi[6] += phi[4];
					phi[7] -= phi[5];
					
					psi[6] += psi[4];
					psi[7] += psi[5];
					
					varphi[6] += varphi[4];
					varphi[7] += varphi[5];
				}
				else if(ix1 == 2){
					theta[6] += theta[4];
					theta[7] += theta[5];
					
					phi[6] += phi[4];
					phi[7] += phi[5];
					
					psi[6] += psi[4];
					psi[7] += psi[5];
					
					varphi[6] += varphi[4];
					varphi[7] += varphi[5];
				}
				else if(ix1 == 3){
					/* Do nothing */
				}
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ix1\n");
				}
			}


			if(indexi1Er >= 0){
				/* for Er */
				Value[indexi1Er] = theta[10];
				Value[indexi1Er+1] = theta[11];


				/* For Fr1 */
				Value[indexi1Fr1] = phi[8];
				Value[indexi1Fr1+1] = phi[9];
				
				/* For Fr2 */
				Value[indexi1Fr2] = psi[8];
				Value[indexi1Fr2+1] = psi[9];
		
				/* For Fr3 */
				Value[indexi1Fr3] = varphi[8];
				Value[indexi1Fr3+1] = varphi[9];

			}		
			else{
				if(ox1 == 1 || ox1 == 5){
					
					theta[6] += theta[10];
					theta[7] -= theta[11];
					
					phi[6] += phi[8];
					phi[7] -= phi[9];
					
					psi[6] += psi[8];
					psi[7] += psi[9];
					
					varphi[6] += varphi[8];
					varphi[7] += varphi[9];
				}
				else if(ox1 == 2){
					theta[6] += theta[10];
					theta[7] += theta[11];
					
					phi[6] += phi[8];
					phi[7] += phi[9];
					
					psi[6] += psi[8];
					psi[7] += psi[9];
					
					varphi[6] += varphi[8];
					varphi[7] += varphi[9];
				}
				else if(ox1 == 3){
					/* Do nothing */
				}
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ox1\n");
				}
					
			}

			if(indexj1Er >= 0){
				/* for Er */
				Value[indexj1Er] = theta[12];
				Value[indexj1Er+1] = theta[13];


				/* For Fr1 */
				Value[indexj1Fr1] = phi[10];
				Value[indexj1Fr1+1] = phi[11];

				/* For Fr2 */
				Value[indexj1Fr2] = psi[10];
				Value[indexj1Fr2+1] = psi[11];

				/* For Fr3 */
				Value[indexj1Fr3] = varphi[10];
				Value[indexj1Fr3+1] = varphi[11];

			}		
			else{
				if(ox2 == 1 || ox2 == 5){
					
					theta[6] += theta[12];
					theta[8] -= theta[13];
					
					phi[6] += phi[10];
					phi[7] += phi[11];
					
					psi[6] += psi[10];
					psi[7] -= psi[11];
					
					varphi[6] += varphi[10];
					varphi[7] += varphi[11];
				}
				else if(ox2 == 2){
					theta[6] += theta[12];
					theta[8] += theta[13];
					
					phi[6] += phi[10];
					phi[7] += phi[11];
					
					psi[6] += psi[10];
					psi[7] += psi[11];
					
					varphi[6] += varphi[10];
					varphi[7] += varphi[11];
				}
				else if(ox2 == 3){
					/* Do nothing */
				}			
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ox2\n");
				}
				
				
			}
				
				

			if(indexk1Er >= 0){
				/* for Er */
				Value[indexk1Er] = theta[14];
				Value[indexk1Er+1] = theta[15];

				/* For Fr1 */
				Value[indexk1Fr1] = phi[12];
				Value[indexk1Fr1+1] = phi[13];
				
				/* For Fr2 */
				Value[indexk1Fr2] = psi[12];
				Value[indexk1Fr2+1] = psi[13];
	
				/* For Fr3 */
				Value[indexk1Fr3] = varphi[12];
				Value[indexk1Fr3+1] = varphi[13];

			}
			else {
				if(ox3 == 1 || ox3 == 5){
					
					theta[6] += theta[14];
					theta[9] -= theta[15];
					
					phi[6] += phi[12];
					phi[7] += phi[13];
					
					psi[6] += psi[12];
					psi[7] += psi[13];
					
					varphi[6] += varphi[12];
					varphi[7] -= varphi[13];
				}
				else if(ox3 == 2){
					theta[6] += theta[14];
					theta[9] += theta[15];
					
					phi[6] += phi[12];
					phi[7] += phi[13];
					
					psi[6] += psi[12];
					psi[7] += psi[13];
					
					varphi[6] += varphi[12];
					varphi[7] += varphi[13];
				}
				else if(ox3 == 3){
					/* Do nothing */
				}
				else{
					ath_error("BackEuler_3d: Unknown boundary condition flag for ox3\n");
				}
			}

				
				
			/* This block needs to be setup at the last
			 * in case matrix elements changed due to boundary 
			 * condition. indexiEr can never be negative *
			 */
				if(indexiEr >= 0){
					/* for Er */
					Value[indexiEr] = theta[6];
					Value[indexiEr+1] = theta[7];
					
					Value[indexiEr+2] = theta[8];
					Value[indexiEr+3] = theta[9];
					
					/* For Fr1 */
					Value[indexiFr1] = phi[6];
					Value[indexiFr1+1] = phi[7];
					
					/* For Fr2 */
					Value[indexiFr2] = psi[6];
					Value[indexiFr2+1] = psi[7];
					
					/* For Fr3 */
					Value[indexiFr3] = varphi[6];
					Value[indexiFr3+1] = varphi[7];
					
				}	
				

			}/* End loop i */
		}/* End loop j */
	}/* End loop k */
			

		/* Initialize the matrix and vectors */
		HYPRE_IJVectorInitialize(RHSEuler2);
		HYPRE_IJVectorInitialize(INIguess2);
		HYPRE_IJMatrixInitialize(Euler);


		/* Set the matrix element */
		HYPRE_IJMatrixSetValues(Euler, lines, ncols, rows, indexValue, Value);
	
	
		HYPRE_IJMatrixAssemble(Euler);
	
	
		HYPRE_IJMatrixGetObject(Euler, (void**) &parcsr_Euler);
	
		HYPRE_IJVectorSetValues(RHSEuler2, lines, rows, RHSEuler);

	
		HYPRE_IJVectorSetValues(INIguess2, lines, rows, INIguess);

	

	
	
	
		HYPRE_IJVectorAssemble(RHSEuler2);
		HYPRE_IJVectorGetObject(RHSEuler2, (void **) &par_RHSEuler2);
	 
		HYPRE_IJVectorAssemble(INIguess2);
		HYPRE_IJVectorGetObject(INIguess2, (void **) &par_INIguess2);


/* debug node */

/*		HYPRE_ParCSRFlexGMRESSetup(solver, parcsr_Euler, par_RHSEuler2, par_INIguess2);
		HYPRE_ParCSRFlexGMRESSolve(solver, parcsr_Euler, par_RHSEuler2, par_INIguess2);
	
		HYPRE_FlexGMRESGetNumIterations(solver, &Matrixiter);
*/


	
if(solverflag == 1 || solverflag == 2){
	HYPRE_ParCSRGMRESSetup(solver, parcsr_Euler, par_RHSEuler2, par_INIguess2);
	HYPRE_ParCSRGMRESSolve(solver, parcsr_Euler, par_RHSEuler2, par_INIguess2);
	
	
	HYPRE_GMRESGetNumIterations(solver, &Matrixiter);
}
else if(solverflag == 3){	
	HYPRE_BoomerAMGSetup(solver,parcsr_Euler,par_RHSEuler2,par_INIguess2);
	HYPRE_BoomerAMGSolve(solver,parcsr_Euler,par_RHSEuler2,par_INIguess2);

	HYPRE_BoomerAMGGetNumIterations(solver,&Matrixiter);	
}
	/* Get solution */
	HYPRE_IJVectorGetValues(INIguess2,lines,rows,INIguess);
	
		ath_pout(0,"Matrix Iteration steps: %d\n",Matrixiter);
		
	/* update the radiation quantities in the mesh */
	for(k=ks;k<=ke;k++){	
		for(j=js;j<=je;j++){
			for(i=is; i<=ie; i++){
			temp0 = pG->U[k][j][i].Er;

				index = 4*(k-ks)*Nx*Ny + 4*(j-js)*Nx + 4*(i-is);
				if(bgflag){
					pG->U[k][j][i].Er  += INIguess[index];
					pG->U[k][j][i].Fr1 += INIguess[index+1];
					pG->U[k][j][i].Fr2 += INIguess[index+2];
					pG->U[k][j][i].Fr3 += INIguess[index+3];
				}
				else{
					pG->U[k][j][i].Er  += INIguess[index];
					pG->U[k][j][i].Fr1 += INIguess[index+1];
					pG->U[k][j][i].Fr2 += INIguess[index+2];
					pG->U[k][j][i].Fr3 += INIguess[index+3];
				}


#ifdef FARGO  
                        cc_pos(pG,i,j,k,&x1,&x2,&x3);
#endif
                	velocity_x = pG->U[k][j][i].M1 /pG->U[k][j][i].d;
                	velocity_y = pG->U[k][j][i].M2 / pG->U[k][j][i].d;
                	velocity_z = pG->U[k][j][i].M3 / pG->U[k][j][i].d;

#ifdef FARGO

                	velocity_y -= qshear * Omega_0 * x1;
#endif

			Fr0x =  pG->U[k][j][i].Fr1 - ((1.0 +  pG->U[k][j][i].Edd_11) * velocity_x +  pG->U[k][j][i].Edd_21 * velocity_y + pG->U[k][j][i].Edd_31 * velocity_z) *  pG->U[k][j][i].Er / Crat; 
			Fr0y =  pG->U[k][j][i].Fr2 - ((1.0 +  pG->U[k][j][i].Edd_22) * velocity_y +  pG->U[k][j][i].Edd_21 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_z) *  pG->U[k][j][i].Er / Crat;
			Fr0z =  pG->U[k][j][i].Fr3 - ((1.0 +  pG->U[k][j][i].Edd_33) * velocity_z +  pG->U[k][j][i].Edd_31 * velocity_x + pG->U[k][j][i].Edd_32 * velocity_y) *  pG->U[k][j][i].Er / Crat; 

			
			/* Estimate the added energy source term */
			if(Prat > 0.0){
				
				pG->U[k][j][i].E += (pG->Eulersource[k][j][i] - dt * (pG->U[k][j][i].Sigma[1] -  pG->U[k][j][i].Sigma[0]) * ( velocity_x * Fr0x + velocity_y * Fr0y + velocity_z * Fr0z));
				
			}
/*
		if(pG->U[ks][j][i].Er < 0.0)
			fprintf(stderr,"[BackEuler_2d]: Negative Radiation energy: %e\n",pG->U[ks][j][i].Er);
*/	
			

			}
		}
	}



	/* Eddington factor is updated in the integrator  */
/*
if(Opacity != NULL){
		for (k=pG->ks; k<=pG->ke; k++){
			for (j=pG->js; j<=pG->je; j++) {
    				for (i=pG->is; i<=pG->ie; i++){
				
				density = pG->U[k][j][i].d;
				
				pressure = (pG->U[k][j][i].E - 0.5 * (pG->U[k][j][i].M1 * pG->U[k][j][i].M1 
				+ pG->U[k][j][i].M2 * pG->U[k][j][i].M2 + pG->U[k][j][i].M3 * pG->U[k][j][i].M3) / density ) * (Gamma - 1);
			
#ifdef RADIATION_MHD
				pressure -= 0.5 * (pG->U[k][j][i].B1c * pG->U[k][j][i].B1c + pG->U[k][j][i].B2c * pG->U[k][j][i].B2c + pG->U[k][j][i].B3c * pG->U[k][j][i].B3c) * (Gamma - 1.0);
#endif

				if(pressure > TINY_NUMBER)
				{
					temperature = pressure / (density * R_ideal);

						
					Opacity(density,temperature,Sigma,NULL);
					for(m=0;m<NOPACITY;m++){
						pG->U[k][j][i].Sigma[m] = Sigma[m];
					}

				}
				else{
					
					pG->Tguess[k][j][i] = pG->U[k][j][i].Er;
				}

			
				}
			}
		}
}

*/
	
/* Update the ghost zones for different boundary condition to be used later */
/*	for (i=0; i<pM->NLevels; i++){ 
            for (j=0; j<pM->DomainsPerLevel[i]; j++){  
        	if (pM->Domain[i][j].Grid != NULL){
  			bvals_radMHD(&(pM->Domain[i][j]));

        	}
      	     }
    	}
*/
	
#ifdef SHEARING_BOX
	/* clean the matrix */
	HYPRE_IJMatrixDestroy(Euler);
#endif


	
  return;	
	

}



/*-------------------------------------------------------------------------*/
/* BackEuler_init_2d: function to allocate memory used just for radiation variables */
/* BackEuler_destruct_2d(): function to free memory */
void BackEuler_init_3d(MeshS *pM)
{

/* Nelements = (je-js+1)*(ie-is+1)*(ke-ks+1) */
/* Nx = (ie-is+1);
   Ny = (je-js+1);
   Nz = (ke-ks+1);
*/
	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	int Nx, Ny, Nz;
	

	Nx = pG->ie - pG->is + 1;
	Ny = pG->je - pG->js + 1;
	Nz = pG->ke - pG->ks + 1;

	NGx = pD->NGrid[0];
	NGy = pD->NGrid[1];
	NGz = pD->NGrid[2];


	int line, Nmatrix, i;
	line = 4*Nx*Ny*Nz;
	Nmatrix = Nx * Ny * Nz;
/* In 3D, we have Er, Fr1, Fr2, Fr3 */


	int    restart = 30;


	/* Allocate value and index array */


/* For temporary vector theta, phi, psi */
	if ((theta = (Real*)malloc(16*sizeof(Real))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	if ((phi = (Real*)malloc(16*sizeof(Real))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	if ((psi = (Real*)malloc(16*sizeof(Real))) == NULL) 
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");

	if ((varphi = (Real*)malloc(16*sizeof(Real))) == NULL)
	ath_error("[BackEuler_init_2d]: malloc returned a NULL pointer\n");
	
/* For Variables used in the matrix solver */	
	
	if((Value = (double *)malloc(58*Nmatrix*sizeof(double))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	if((ncols = (int *)malloc(line*sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	if((rows = (int *)malloc(line*sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	if((RHSEuler = (double *)malloc(line*sizeof(double))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	if((INIguess = (double *)malloc(line*sizeof(double))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	if ((indexValue = (int*)malloc(58*Nmatrix*sizeof(int))) == NULL) 
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	
	if ((ptr = (int*)malloc(line*sizeof(int))) == NULL) 
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	
	if ((NZrow = (int*)malloc(line*sizeof(int))) == NULL) 
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	
	
	
	
	/* NZrow is the estimated size of each row */
	for(i=0; i<Nmatrix; i++){
		NZrow[4*i]   = 16;
		NZrow[4*i+1] = 14;
		NZrow[4*i+2] = 14;
		NZrow[4*i+3] = 14;
	}
	
	
	

	if((colindex = (int **)calloc_2d_array(Nmatrix, 28, sizeof(int))) == NULL)
		ath_error("[BackEuler_init_3d]: malloc returned a NULL pointer\n");
	
/* Objects for the matrix solver */	

/* BoomerAMG solver */
    if(solverflag == 3){
	HYPRE_BoomerAMGCreate(&solver);
	HYPRE_BoomerAMGSetPrintLevel(solver, 0); 
	HYPRE_BoomerAMGSetCoarsenType(solver, 3);
	HYPRE_BoomerAMGSetRelaxType(solver, 4); 
	HYPRE_BoomerAMGSetNumSweeps(solver, 1);
	HYPRE_BoomerAMGSetTol(solver, 1.e-10); 
	HYPRE_BoomerAMGSetMaxIter(solver, 2000); 
	HYPRE_BoomerAMGSetStrongThreshold(solver,0.5);
   }

	


/* GMRES Method */
   if(solverflag == 1 || solverflag == 2){	
	HYPRE_ParCSRGMRESCreate(MPI_COMM_WORLD, &solver);
	
	
	HYPRE_GMRESSetKDim(solver, restart);
	HYPRE_GMRESSetMaxIter(solver, 2000); 
	HYPRE_GMRESSetTol(solver, 1e-10); 
	
    	HYPRE_GMRESSetPrintLevel(solver, 0); 
	HYPRE_GMRESSetLogging(solver, 1);
    }
  
/* BoomerAMGC preconditioner */
  if(solverflag == 1){	
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 0); 
	HYPRE_BoomerAMGSetCoarsenType(precond, 3);
	HYPRE_BoomerAMGSetRelaxType(precond, 4); 
	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); 
	HYPRE_BoomerAMGSetMaxIter(precond, 2); 
	HYPRE_BoomerAMGSetStrongThreshold(precond,0.5);
	
	HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
						  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
	
   }
	
	/* Euclid preconditioner */
   if(solverflag == 2){
	HYPRE_EuclidCreate(MPI_COMM_WORLD, &precond);

	HYPRE_EuclidSetLevel(precond, 1);
 

	
	HYPRE_GMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_EuclidSolve,
						  (HYPRE_PtrToSolverFcn) HYPRE_EuclidSetup, precond);
	
  }
/************************************************/
/* Flexible GMRES method */
	
/*	HYPRE_ParCSRFlexGMRESCreate(MPI_COMM_WORLD, &solver);
	
	
	HYPRE_FlexGMRESSetKDim(solver, restart);
	HYPRE_FlexGMRESSetMaxIter(solver, 2000); 
	HYPRE_FlexGMRESSetTol(solver, 1e-8); 
	HYPRE_FlexGMRESSetPrintLevel(solver, 0); 
	HYPRE_FlexGMRESSetLogging(solver, 1); 
	
	
	HYPRE_BoomerAMGCreate(&precond);
	HYPRE_BoomerAMGSetPrintLevel(precond, 0); 
	HYPRE_BoomerAMGSetCoarsenType(precond, 6);
	HYPRE_BoomerAMGSetRelaxType(precond, 6); 

	HYPRE_BoomerAMGSetNumSweeps(precond, 1);
	HYPRE_BoomerAMGSetTol(precond, 0.0); 
	HYPRE_BoomerAMGSetMaxIter(precond, 1); 
	
	
	HYPRE_FlexGMRESSetPrecond(solver, (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSolve,
							  (HYPRE_PtrToSolverFcn) HYPRE_BoomerAMGSetup, precond);
	
*/	
	
	

	
	return;

}


void BackEuler_destruct_3d()
{

	


/* For temporary vector theta, phi, psi */
	if(theta != NULL) free(theta);
	if(phi != NULL) free(phi);
	if(psi != NULL) free(psi);
	if(varphi != NULL) free(varphi);


	
	if(Value != NULL) free(Value);
	if(indexValue != NULL) free(indexValue);
	if(ptr != NULL) free(ptr);
	if(ncols != NULL) free(ncols);
	if(rows != NULL) free(rows);
	if(RHSEuler != NULL) free(RHSEuler);
	if(INIguess != NULL) free(INIguess);
	if(NZrow != NULL) free(NZrow);
	
	free_2d_array(colindex);
	
	/* clean the matrix */
	
#ifndef SHEARING_BOX
	HYPRE_IJMatrixDestroy(Euler);
#endif
	
    HYPRE_IJVectorDestroy(RHSEuler2);
	HYPRE_IJVectorDestroy(INIguess2);
	
     if(solverflag == 1){
	HYPRE_ParCSRGMRESDestroy(solver);
     }

     if(solverflag == 2){  
	HYPRE_EuclidDestroy(precond);
    }
	
/*	HYPRE_ParCSRFlexGMRESDestroy(solver);
*/
   if(solverflag == 1){
	HYPRE_BoomerAMGDestroy(precond);
    }

     if(solverflag == 3){
	HYPRE_BoomerAMGDestroy(solver);	
     }	

}





/*********************************************************************************/
/* ========= Private Function =====================*/
/* ===============================================*/

/* This function is used to set the positions of different elements */
/* Which are unchanged most of the time */
/* For shearing boundary condition, the boundary cells need to be updated again */

void INI_matrix(MeshS *pM)
{
	
	/* NZ_NUM is the number of non-zero elements in Matrix. It may change if periodic boundary condition is applied */
	/* lines is the size of the partial matrix */
	/* count is the number of total non-zeros before that row */
	
	
	
	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;

	
	int i, j, k;
	is = pG->is;
	ie = pG->ie;
	js = pG->js;
	je = pG->je;
	ks = pG->ks;
	ke = pG->ke;
	int Nmatrix, NZ_NUM, count;
	
	/* Nmatrix is the number of active cells just in this grids */
	/* Matrix size should be 3 * Nmatrix, ghost zones are not included*/
   	Nmatrix = Ny * Nx * Nz;

	/* total number of lines in the grid. This is also the same for every grid */
	
	

	
	
	/* Number of non-zero elements is the same for both periodic 
	 * boundary condition and MPI boundary. 
	 * It is just that index will be different *
	 */
	
	NoEr = 0;	/* Position of first non-zero element in row Er */
	NoFr1 = 0;	/* Position of first non-zero element in row Fr1 */
	NoFr2 = 0;	/* POsition of first non-zero element in row Fr2 */
	NoFr3 = 0;
	count = 0;
	/* For non-periodic boundary condition, this number will change */
	
	
	/* total number of lines before this grids. This is also the starting line for this grid */ 
	/* ID first goes from x direction, then goes from y direction and finally goes from z direction */ 
	
	
	
	NZ_NUM = 58 * Nmatrix; 
	

	
	/* Create the right hand side and guess solution */
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&RHSEuler2);
	HYPRE_IJVectorSetObjectType(RHSEuler2, HYPRE_PARCSR);

	
	HYPRE_IJVectorCreate(MPI_COMM_WORLD, ilower, iupper,&INIguess2);
	HYPRE_IJVectorSetObjectType(INIguess2, HYPRE_PARCSR);
	
	
	/* Create the matrix */
	HYPRE_IJMatrixCreate(MPI_COMM_WORLD, ilower, iupper, ilower, iupper, &Euler);
	HYPRE_IJMatrixSetObjectType(Euler, HYPRE_PARCSR);
	
	/* set the maximum number of elements in each row */
	HYPRE_IJMatrixSetRowSizes(Euler, NZrow);
	
	/* set the maximum number of elements that will be set by other process */
/*	HYPRE_IJMatrixSetMaxOffProcElmts(Euler, 0);
*/
	
	for(i=0;i<lines; i++){
		rows[i] = ilower + i;
		
	}
	

	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){
			for(i=is; i<=ie; i++){


	
	
	/* Initialize the block pointer *
	 * For non-periodic nonMPI case, some blocks does not exist 
	 */

	/* We should do this for every line */

	indexk0Er = -1; indexk0Fr1 = -1; indexk0Fr2 = -1; indexk0Fr3 = -1;
	indexj0Er = -1; indexj0Fr1 = -1; indexj0Fr2 = -1; indexj0Fr3 = -1;
	indexi0Er = -1; indexi0Fr1 = -1; indexi0Fr2 = -1; indexi0Fr3 = -1;
	indexiEr = -1;  indexiFr1 = -1;  indexiFr2 = -1;  indexiFr3 = -1;
	indexi1Er = -1; indexi1Fr1 = -1; indexi1Fr2 = -1; indexi1Fr3 = -1;
	indexj1Er = -1; indexj1Fr1 = -1; indexj1Fr2 = -1; indexj1Fr3 = -1;
	indexk1Er = -1; indexk1Fr1 = -1; indexk1Fr2 = -1; indexk1Fr3 = -1;
	
				
				if(i == is){
					if(j == js){
						if(k == ks){
							
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for y direction*/
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}	
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for y direction */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								is_js_ks_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								is_js_ks_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								is_js_ks_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								is_js_ks_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								is_js_ks_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								is_js_ks_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								is_js_ks_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								is_js_ks_MPI_MPI_MPI();
						}/* End k == ks, j == js, i == is */
						
						else if(k == ke){
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for y direction */
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								} /* non periodic boundary for y direction */
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for y direction */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								is_js_ke_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								is_js_ke_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								is_js_ke_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								is_js_ke_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								is_js_ke_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								is_js_ke_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								is_js_ke_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								is_js_ke_MPI_MPI_MPI();
							
						}/* End k == ke, j = js, i = is */
						else{
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
									
								} /* either periodic or MPI boundary condition for y direction  */
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{						
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
									
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0))
								is_js_k_phy_phy(k);
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0))
								is_js_k_phy_MPI(k);				
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0))
								is_js_k_MPI_phy(k);
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0))
								is_js_k_MPI_MPI(k);
							
						}/* End k!=ks and k!=ke, j=js, i = is */	
					}/* End j == js , i = is */
					else if(j == je){
						if(k == ks){	
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}	
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for y direction */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition y direction */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								is_je_ks_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								is_je_ks_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								is_je_ks_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								is_je_ks_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								is_je_ks_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								is_je_ks_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								is_je_ks_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								is_je_ks_MPI_MPI_MPI();
						}/* End k == ks , j = je, i = is*/
						else if(k == ke){
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								} /* End non periodic boundary condition for y direction */
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition y direction */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}/* End non-periodic for y direction */
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition y direction */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								is_je_ke_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								is_je_ke_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								is_je_ke_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								is_je_ke_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								is_je_ke_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								is_je_ke_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								is_je_ke_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								is_je_ke_MPI_MPI_MPI();
							
						}/* End k == ke, j = je, i = is */
						/* begin k != ks and k != ke, j=je, i = is */
						else{
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
									
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{						
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
									
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0))
								is_je_k_phy_phy(k);
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0))
								is_je_k_phy_MPI(k);				
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0))
								is_je_k_MPI_phy(k);
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0))
								is_je_k_MPI_MPI(k);
							
							
						} /* End k != ke , k != ks*/
					}/* End j=je, Now begin j != js && j != je */			
					else{
						if(k == ks){	
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* Non periodic for x1 */
							else{
								
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx3_id < 0))
								is_j_ks_phy_phy(j);
							else if((pG->lx1_id < 0) && (pG->lx3_id >= 0))
								is_j_ks_phy_MPI(j);
							else if((pG->lx1_id >= 0) && (pG->lx3_id < 0))
								is_j_ks_MPI_phy(j);
							else if((pG->lx1_id >= 0) && (pG->lx3_id >= 0))
								is_j_ks_MPI_MPI(j);
						}/* End k == ks */
						else if(k == ke){
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* Non periodic for x1 */
							else{
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx3_id < 0))
								is_j_ke_phy_phy(j);
							else if((pG->lx1_id < 0) && (pG->rx3_id >= 0))
								is_j_ke_phy_MPI(j);
							else if((pG->lx1_id >= 0) && (pG->rx3_id < 0))
								is_j_ke_MPI_phy(j);
							else if((pG->lx1_id >= 0) && (pG->rx3_id >= 0))
								is_j_ke_MPI_MPI(j);
							
						}/* End k == ke */
						else{
							if((ix1 != 4) && (pG->lx1_id < 0)){
								/* physical boundary on the left */						
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
								
								
							}/* Non periodic for x1 */
							else{
								
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0))
								is_j_k_phy(j,k);
							else
								is_j_k_MPI(j,k);								
							
						} /* End k!=ks k!= ke */
					}/* End j!=js, j!=je, i=is */
				}/* End i==is and begin i=ie */
				else if (i == ie){
					if(j == js){
						if(k == ks){	
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition for z direction*/
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								ie_js_ks_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								ie_js_ks_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								ie_js_ks_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								ie_js_ks_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								ie_js_ks_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								ie_js_ks_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								ie_js_ks_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								ie_js_ks_MPI_MPI_MPI();
						}/* End k == ks */
						else if(k == ke){
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								ie_js_ke_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								ie_js_ke_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								ie_js_ke_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								ie_js_ke_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								ie_js_ke_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								ie_js_ke_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								ie_js_ke_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								ie_js_ke_MPI_MPI_MPI();
							
						}/* End k == ke */
						else{/* Now begin k!=ks && k!= ke */
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ix2 !=4) && (pG->lx2_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
									
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ix2 !=4) && (pG->lx2_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{						
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
									
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0))
								ie_js_k_phy_phy(k);
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0))
								ie_js_k_phy_MPI(k);				
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0))
								ie_js_k_MPI_phy(k);
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0))
								ie_js_k_MPI_MPI(k);
							
						}/* End k!=ks and k!=ke */	
					}/* End j == js */
					else if(j == je){
						if(k == ks){	
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}
								else{
									
									if((ix3 != 4) && (pG->lx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								ie_je_ks_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								ie_je_ks_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								ie_je_ks_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								ie_je_ks_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								ie_je_ks_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								ie_je_ks_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								ie_je_ks_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								ie_je_ks_MPI_MPI_MPI();
						}/* End k == ks */
						else if(k == ke){
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 24;
										NoEr = count;
										NoFr1 = count + 10;
										NoFr2 = count + 18;
										NoFr3 = count + 26;
										count += 34;
									}
									else{
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									} /* End for z direction */
									
								}
								else{
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 16;
										NoEr = count;
										NoFr1 = count + 12;
										NoFr2 = count + 22;
										NoFr3 = count + 32;
										count += 42;
									}
									else{
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									} /* End for z direction */
									
								}/* non periodic for je */
								else{
									
									if((ox3 != 4) && (pG->rx3_id < 0)){
										
										NZ_NUM -= 8;
										NoEr = count;
										NoFr1 = count + 14;
										NoFr2 = count + 26;
										NoFr3 = count + 38;
										count += 50;
									}
									else{
										/* NZ_NUM do not change */
										NoEr = count;
										NoFr1 = count + 16;
										NoFr2 = count + 30;
										NoFr3 = count + 44;
										count += 58;
									} /* End for z direction */
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								ie_je_ke_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								ie_je_ke_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								ie_je_ke_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								ie_je_ke_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								ie_je_ke_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								ie_je_ke_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								ie_je_ke_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								ie_je_ke_MPI_MPI_MPI();
							
						}/* End k == ke */
						else{/* begin k!= ks and k!= ke */
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								if((ox2 !=4) && (pG->rx2_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
									
								} /* either periodic or MPI boundary condition */
							}/* Non periodic for x1 */
							else{
								if((ox2 !=4) && (pG->rx2_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{						
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
									
								} /* either periodic or MPI boundary condition */
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0))
								ie_je_k_phy_phy(k);
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0))
								ie_je_k_phy_MPI(k);				
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0))
								ie_je_k_MPI_phy(k);
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0))
								ie_je_k_MPI_MPI(k);
							
						}/* End k!= ks and k!= ke */
					} /* End j == je */
					else{/* begin j!= js and j!= je */
						if(k == ks){	
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* Non periodic for x1 */
							else{
								
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx3_id < 0))
								ie_j_ks_phy_phy(j);
							else if((pG->rx1_id < 0) && (pG->lx3_id >= 0))
								ie_j_ks_phy_MPI(j);
							else if((pG->rx1_id >= 0) && (pG->lx3_id < 0))
								ie_j_ks_MPI_phy(j);
							else if((pG->rx1_id >= 0) && (pG->lx3_id >= 0))
								ie_j_ks_MPI_MPI(j);
						}/* End k == ks */
						else if(k == ke){
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* Non periodic for x1 */
							else{
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx3_id < 0))
								ie_j_ke_phy_phy(j);
							else if((pG->rx1_id < 0) && (pG->rx3_id >= 0))
								ie_j_ke_phy_MPI(j);
							else if((pG->rx1_id >= 0) && (pG->rx3_id < 0))
								ie_j_ke_MPI_phy(j);
							else if((pG->rx1_id >= 0) && (pG->rx3_id >= 0))
								ie_j_ke_MPI_MPI(j);
							
						}/* End k == ke */
						else{/* begin k!= ks && k!= ke */
							if((ox1 != 4) && (pG->rx1_id < 0)){
								/* physical boundary on the left */						
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
								
								
							}/* Non periodic for x1 */
							else{
								
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
								
							}/* either periodic for x1 or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0))
								ie_j_k_phy(j,k);
							else
								ie_j_k_MPI(j,k);								
						}/* End k!=ks && k!= ke */
					} /* End j!= js & j != je */
				}/* End i==ie */
				else{/* begin i!= is && i!= ie */
					
					if(j == js){
						if(k == ks){	
							if((ix2 !=4) && (pG->lx2_id < 0)){
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}
							else{
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
							} /* either periodic or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx2_id < 0) && (pG->lx3_id < 0))
								i_js_ks_phy_phy(i);
							else if((pG->lx2_id < 0) && (pG->lx3_id >= 0))
								i_js_ks_phy_MPI(i);
							else if((pG->lx2_id >= 0) && (pG->lx3_id < 0))
								i_js_ks_MPI_phy(i);
							else if( (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								i_js_ks_MPI_MPI(i);
						}/* End k == ks */
						else if(k == ke){
							
							if((ix2 !=4) && (pG->lx2_id < 0)){
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* non-periodic for x1 */
							else{
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
							} /* either periodic or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->lx2_id < 0) && (pG->rx3_id < 0))
								i_js_ke_phy_phy(i);
							else if((pG->lx2_id < 0) && (pG->rx3_id >= 0))
								i_js_ke_phy_MPI(i);
							else if((pG->lx2_id >= 0) && (pG->rx3_id < 0))
								i_js_ke_MPI_phy(i);
							else if((pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								i_js_ke_MPI_MPI(i);				
							
						}/* End k == ke */
						else{
							
							if((ix2 !=4) && (pG->lx2_id < 0)){
								
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
							}
							else{						
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
								
							} /* either periodic or MPI boundary condition */
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if(pG->lx2_id < 0)
								i_js_k_phy(i,k);
							else 
								i_js_k_MPI(i,k);					
							
						}/* End k!=ks and k!=ke */	
					}/* End j == js */
					else if(j == je){
						if(k == ks){	
							
							if((ox2 !=4) && (pG->rx2_id < 0)){
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}/* non periodic for x2 */
							else{
								
								if((ix3 != 4) && (pG->lx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
							} /* either periodic or MPI boundary condition */
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx2_id < 0) && (pG->lx3_id < 0))
								i_je_ks_phy_phy(i);
							else if((pG->rx2_id < 0) && (pG->lx3_id >= 0))
								i_je_ks_phy_MPI(i);
							else if((pG->rx2_id >= 0) && (pG->lx3_id < 0))
								i_je_ks_MPI_phy(i);
							else if((pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								i_je_ks_MPI_MPI(i);
							
						}/* End k == ks */
						else if(k == ke){
							
							if((ox2 !=4) && (pG->rx2_id < 0)){
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 16;
									NoEr = count;
									NoFr1 = count + 12;
									NoFr2 = count + 22;
									NoFr3 = count + 32;
									count += 42;
								}
								else{
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								} /* End for z direction */
								
							}
							else{
								
								if((ox3 != 4) && (pG->rx3_id < 0)){
									
									NZ_NUM -= 8;
									NoEr = count;
									NoFr1 = count + 14;
									NoFr2 = count + 26;
									NoFr3 = count + 38;
									count += 50;
								}
								else{
									/* NZ_NUM do not change */
									NoEr = count;
									NoFr1 = count + 16;
									NoFr2 = count + 30;
									NoFr3 = count + 44;
									count += 58;
								} /* End for z direction */
							} /* either periodic or MPI boundary condition */
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx2_id < 0) && (pG->rx3_id < 0))
								i_je_ke_phy_phy(i);
							else if((pG->rx2_id < 0) && (pG->rx3_id >= 0))
								i_je_ke_phy_MPI(i);
							else if((pG->rx2_id >= 0) && (pG->rx3_id < 0))
								i_je_ke_MPI_phy(i);
							else if((pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								i_je_ke_MPI_MPI(i);
							
							
						}/* End k == ke */
						else{
							
							if((ox2 !=4) && (pG->rx2_id < 0)){
								
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
							}
							else{						
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
								
							} /* either periodic or MPI boundary condition */
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if((pG->rx2_id < 0))
								i_je_k_phy(i,k);
							else
								i_je_k_MPI(i,k);				
							
						} /* End k!= ks && k!= ke */
					}/* End j== je */
					else {/* begin j!=js && j!= je */
						if(k == ks){	
							if((ix3 != 4) && (pG->lx3_id < 0)){
								
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
							}
							else{
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
							} /* End for z direction */
							
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if(pG->lx3_id < 0)
								i_j_ks_phy(i,j);
							else
								i_j_ks_MPI(i,j);
						}/* End k == ks */
						else if(k == ke){
							
							if((ox3 != 4) && (pG->rx3_id < 0)){
								
								NZ_NUM -= 8;
								NoEr = count;
								NoFr1 = count + 14;
								NoFr2 = count + 26;
								NoFr3 = count + 38;
								count += 50;
							}
							else{
								/* NZ_NUM do not change */
								NoEr = count;
								NoFr1 = count + 16;
								NoFr2 = count + 30;
								NoFr3 = count + 44;
								count += 58;
							} /* End for z direction */
							
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							/* judge MPI or physical boundary condition */
							if(pG->rx3_id < 0)
								i_j_ke_phy(i,j);
							else 
								i_j_ke_MPI(i,j);
							
							
						}/* End k == ke */
						else{
							
							
							/* NZ_NUM do not change */
							NoEr = count;
							NoFr1 = count + 16;
							NoFr2 = count + 30;
							NoFr3 = count + 44;
							count += 58;								
							
							
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoEr;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr1;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr2;
							ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = NoFr3;
							
							i_j_k(i,j,k);	
							
							
						}/* End k!=ks & k != ke */ 	
						
					} /* End j!= js & j != je */
				}/* End i!=is & i!=ie */
			
				
	/* There are 7 blocks */
	if(indexk0Er >= 0){
		/* for Er */
		indexValue[indexk0Er] = colk0;
					
		indexValue[indexk0Er+1] = colk0 + 3;
					
		/* For Fr1 */
		indexValue[indexk0Fr1] = colk0;
					
		indexValue[indexk0Fr1+1] = colk0 + 1;
					
					
		/* For Fr2 */

		indexValue[indexk0Fr2] = colk0;
					
		indexValue[indexk0Fr2+1] = colk0 + 2;
					
					
		/* For Fr3 */

		indexValue[indexk0Fr3] = colk0;
					

		indexValue[indexk0Fr3+1] = colk0 + 3;
	}
				
	if(indexj0Er >= 0){
		/* for Er */

		indexValue[indexj0Er] = colj0;
					

		indexValue[indexj0Er+1] = colj0 + 2;
					
		/* For Fr1 */
		
		indexValue[indexj0Fr1] = colj0;
					
		
		indexValue[indexj0Fr1+1] = colj0 + 1;
					
					
		/* For Fr2 */
		
		indexValue[indexj0Fr2] = colj0;
					
		
		indexValue[indexj0Fr2+1] = colj0 + 2;
					
					
		/* For Fr3 */
		
		indexValue[indexj0Fr3] = colj0;
					
		
		indexValue[indexj0Fr3+1] = colj0 + 3;
	}		
				
	if(indexi0Er >= 0){
		/* for Er */
		
	    indexValue[indexi0Er] = coli0;
					

		indexValue[indexi0Er+1] = coli0 + 1;
					
		/* For Fr1 */

		indexValue[indexi0Fr1] = coli0;
			
		indexValue[indexi0Fr1+1] = coli0 + 1;
					
					
		/* For Fr2 */

		indexValue[indexi0Fr2] = coli0;
					
		indexValue[indexi0Fr2+1] = coli0 + 2;
					
					
		/* For Fr3 */
		indexValue[indexi0Fr3] = coli0;
					
		indexValue[indexi0Fr3+1] = coli0 + 3;
	}					
				
				
	if(indexiEr >= 0){
		/* for Er */

		indexValue[indexiEr] = coli;
					
		indexValue[indexiEr+1] = coli + 1;
					

		indexValue[indexiEr+2] = coli + 2;
					
		indexValue[indexiEr+3] = coli + 3;
					
		/* For Fr1 */
		indexValue[indexiFr1] = coli;
					

		indexValue[indexiFr1+1] = coli + 1;
					
					
		/* For Fr2 */
		indexValue[indexiFr2] = coli;
					
		indexValue[indexiFr2+1] = coli + 2;
					
					
		/* For Fr3 */
		indexValue[indexiFr3] = coli;
					
		indexValue[indexiFr3+1] = coli + 3;
	}	
				
	if(indexi1Er >= 0){
		/* for Er */
		indexValue[indexi1Er] = coli1;
					

		indexValue[indexi1Er+1] = coli1 + 1;
					
		/* For Fr1 */
		indexValue[indexi1Fr1] = coli1;
					
		indexValue[indexi1Fr1+1] = coli1 + 1;
					
					
		/* For Fr2 */

		indexValue[indexi1Fr2] = coli1;
					
		indexValue[indexi1Fr2+1] = coli1 + 2;
					
					
		/* For Fr3 */

		indexValue[indexi1Fr3] = coli1;
					
		indexValue[indexi1Fr3+1] = coli1 + 3;
	}		
				
	if(indexj1Er >= 0){
	/* for Er */
		indexValue[indexj1Er] = colj1;
		
		indexValue[indexj1Er+1] = colj1 + 2;
					
		/* For Fr1 */
		indexValue[indexj1Fr1] = colj1;
					
		
		indexValue[indexj1Fr1+1] = colj1 + 1;
					
					
		/* For Fr2 */
		indexValue[indexj1Fr2] = colj1;
					

		indexValue[indexj1Fr2+1] = colj1 + 2;
					
					
		/* For Fr3 */
	
		indexValue[indexj1Fr3] = colj1;
					
		indexValue[indexj1Fr3+1] = colj1 + 3;
	}		
				
	if(indexk1Er >= 0){
		/* for Er */

		indexValue[indexk1Er] = colk1;
					

		indexValue[indexk1Er+1] = colk1 + 3;
					
		/* For Fr1 */
		indexValue[indexk1Fr1] = colk1;
					
		indexValue[indexk1Fr1+1] = colk1 + 1;
					
					
		/* For Fr2 */

		indexValue[indexk1Fr2] = colk1;
					
		indexValue[indexk1Fr2+1] = colk1 + 2;
					
					
		/* For Fr3 */

		indexValue[indexk1Fr3] = colk1;
					

		indexValue[indexk1Fr3+1] = colk1 + 3;
	}
				
				
				
	ncols[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)] = NoFr1 - NoEr;
	ncols[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1] = NoFr2 - NoFr1;
	ncols[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2] = NoFr3 - NoFr2;
	ncols[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3] = count - NoFr3;
				
		/* Save the index for furture use */
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][0] = indexk0Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][1] = indexk0Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][2] = indexk0Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][3] = indexk0Fr3;			
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][4] = indexj0Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][5] = indexj0Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][6] = indexj0Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][7] = indexj0Fr3;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][8] = indexi0Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][9] = indexi0Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][10] = indexi0Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][11] = indexi0Fr3;			
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][12] = indexiEr;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][13] = indexiFr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][14] = indexiFr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][15] = indexiFr3;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][16] = indexi1Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][17] = indexi1Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][18] = indexi1Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][19] = indexi1Fr3;			
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][20] = indexj1Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][21] = indexj1Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][22] = indexj1Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][23] = indexj1Fr3;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][24] = indexk1Er;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][25] = indexk1Fr1;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][26] = indexk1Fr2;
	colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][27] = indexk1Fr3;			

				
				
			}/* End loop i */
		}/* End loop j */
	}/* End loop k */
	

	
	return;
	
}

/* Update the boundary cells for shearing box boundary condition */
#ifdef SHEARING_BOX
void shearing_matrix(MeshS *pM)
{

	/* Only apply shearing periodic boundary condition for physical boundaries. */
	
	DomainS *pD;
	pD= &(pM->Domain[0][0]);
	
	GridS *pG=pD->Grid;
	
	
	int i, j, k;
	is = pG->is;
	ie = pG->ie;
	js = pG->js;
	je = pG->je;
	ks = pG->ks;
	ke = pG->ke;



	
	/* First, left x boundary */
	
	i = is;
if((pG->lx1_id < 0) || (pG->lx1_id > ID))
{	
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){

	/* We need to do this for every line */

	indexk0Er = -1; indexk0Fr1 = -1; indexk0Fr2 = -1; indexk0Fr3 = -1;
	indexj0Er = -1; indexj0Fr1 = -1; indexj0Fr2 = -1; indexj0Fr3 = -1;
	indexi0Er = -1; indexi0Fr1 = -1; indexi0Fr2 = -1; indexi0Fr3 = -1;
	indexiEr = -1;  indexiFr1 = -1;  indexiFr2 = -1;  indexiFr3 = -1;
	indexi1Er = -1; indexi1Fr1 = -1; indexi1Fr2 = -1; indexi1Fr3 = -1;
	indexj1Er = -1; indexj1Fr1 = -1; indexj1Fr2 = -1; indexj1Fr3 = -1;
	indexk1Er = -1; indexk1Fr1 = -1; indexk1Fr2 = -1; indexk1Fr3 = -1;





			
			
			/* We need NoEr, NoFr1, NoFr2, NoFr3 to determine the positions of the elements. */
			NoEr = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)];
			NoFr1 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1];
			NoFr2 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2];
			NoFr3 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3];
				

					if(j == js){
						if(k == ks){
														
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								is_js_ks_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								is_js_ks_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								is_js_ks_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								is_js_ks_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								is_js_ks_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								is_js_ks_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								is_js_ks_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								is_js_ks_MPI_MPI_MPI();
						}/* End k == ks, j == js, i == is */
						
						else if(k == ke){
														
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								is_js_ke_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								is_js_ke_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								is_js_ke_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								is_js_ke_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								is_js_ke_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								is_js_ke_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								is_js_ke_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								is_js_ke_MPI_MPI_MPI();
							
						}/* End k == ke, j = js, i = is */
						else{
													
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx2_id < 0))
								is_js_k_phy_phy(k);
							else if((pG->lx1_id < 0) && (pG->lx2_id >= 0))
								is_js_k_phy_MPI(k);				
							else if((pG->lx1_id >= 0) && (pG->lx2_id < 0))
								is_js_k_MPI_phy(k);
							else if((pG->lx1_id >= 0) && (pG->lx2_id >= 0))
								is_js_k_MPI_MPI(k);
							
						}/* End k!=ks and k!=ke, j=js, i = is */	
					}/* End j == js , i = is */
					else if(j == je){
						if(k == ks){	
														
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								is_je_ks_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								is_je_ks_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								is_je_ks_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								is_je_ks_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								is_je_ks_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								is_je_ks_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								is_je_ks_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								is_je_ks_MPI_MPI_MPI();
						}/* End k == ks , j = je, i = is*/
						else if(k == ke){
							
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								is_je_ke_phy_phy_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								is_je_ke_phy_phy_MPI();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								is_je_ke_phy_MPI_phy();
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								is_je_ke_phy_MPI_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								is_je_ke_MPI_phy_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								is_je_ke_MPI_phy_MPI();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								is_je_ke_MPI_MPI_phy();
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								is_je_ke_MPI_MPI_MPI();
							
						}/* End k == ke, j = je, i = is */
						/* begin k != ks and k != ke, j=je, i = is */
						else{
														
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx2_id < 0))
								is_je_k_phy_phy(k);
							else if((pG->lx1_id < 0) && (pG->rx2_id >= 0))
								is_je_k_phy_MPI(k);				
							else if((pG->lx1_id >= 0) && (pG->rx2_id < 0))
								is_je_k_MPI_phy(k);
							else if((pG->lx1_id >= 0) && (pG->rx2_id >= 0))
								is_je_k_MPI_MPI(k);
							
							
						} /* End k != ke , k != ks*/
					}/* End j=je, Now begin j != js && j != je */			
					else{
						if(k == ks){	
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->lx3_id < 0))
								is_j_ks_phy_phy(j);
							else if((pG->lx1_id < 0) && (pG->lx3_id >= 0))
								is_j_ks_phy_MPI(j);
							else if((pG->lx1_id >= 0) && (pG->lx3_id < 0))
								is_j_ks_MPI_phy(j);
							else if((pG->lx1_id >= 0) && (pG->lx3_id >= 0))
								is_j_ks_MPI_MPI(j);
						}/* End k == ks */
						else if(k == ke){
					
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0) && (pG->rx3_id < 0))
								is_j_ke_phy_phy(j);
							else if((pG->lx1_id < 0) && (pG->rx3_id >= 0))
								is_j_ke_phy_MPI(j);
							else if((pG->lx1_id >= 0) && (pG->rx3_id < 0))
								is_j_ke_MPI_phy(j);
							else if((pG->lx1_id >= 0) && (pG->rx3_id >= 0))
								is_j_ke_MPI_MPI(j);
							
						}/* End k == ke */
						else{
		
							
							/* judge MPI or physical boundary condition */
							if((pG->lx1_id < 0))
								is_j_k_phy(j,k);
							else
								is_j_k_MPI(j,k);								
							
						} /* End k!=ks k!= ke */
					}/* End j!=js, j!=je, i=is */


						
				/* There are 7 blocks */
				if(indexk0Er >= 0){
					/* for Er */
					indexValue[indexk0Er] = colk0;
					
					indexValue[indexk0Er+1] = colk0 + 3;
					
					/* For Fr1 */
					indexValue[indexk0Fr1] = colk0;
					
					indexValue[indexk0Fr1+1] = colk0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexk0Fr2] = colk0;
					
					indexValue[indexk0Fr2+1] = colk0 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexk0Fr3] = colk0;
					
					
					indexValue[indexk0Fr3+1] = colk0 + 3;
				}
				
				if(indexj0Er >= 0){
					/* for Er */
					
					indexValue[indexj0Er] = colj0;
					
					
					indexValue[indexj0Er+1] = colj0 + 2;
					
					/* For Fr1 */
					
					indexValue[indexj0Fr1] = colj0;
					
					
					indexValue[indexj0Fr1+1] = colj0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexj0Fr2] = colj0;
					
					
					indexValue[indexj0Fr2+1] = colj0 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexj0Fr3] = colj0;
					
					
					indexValue[indexj0Fr3+1] = colj0 + 3;
				}		
				
				if(indexi0Er >= 0){
					/* for Er */
					
					indexValue[indexi0Er] = coli0;
					
					
					indexValue[indexi0Er+1] = coli0 + 1;
					
					/* For Fr1 */
					
					indexValue[indexi0Fr1] = coli0;
					
					indexValue[indexi0Fr1+1] = coli0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexi0Fr2] = coli0;
					
					indexValue[indexi0Fr2+1] = coli0 + 2;
					
					
					/* For Fr3 */
					indexValue[indexi0Fr3] = coli0;
					
					indexValue[indexi0Fr3+1] = coli0 + 3;
				}					
				
				
				if(indexiEr >= 0){
					/* for Er */
					
					indexValue[indexiEr] = coli;
					
					indexValue[indexiEr+1] = coli + 1;
					
					
					indexValue[indexiEr+2] = coli + 2;
					
					indexValue[indexiEr+3] = coli + 3;
					
					/* For Fr1 */
					indexValue[indexiFr1] = coli;
					
					
					indexValue[indexiFr1+1] = coli + 1;
					
					
					/* For Fr2 */
					indexValue[indexiFr2] = coli;
					
					indexValue[indexiFr2+1] = coli + 2;
					
					
					/* For Fr3 */
					indexValue[indexiFr3] = coli;
					
					indexValue[indexiFr3+1] = coli + 3;
				}	
				
				if(indexi1Er >= 0){
					/* for Er */
					indexValue[indexi1Er] = coli1;
					
					
					indexValue[indexi1Er+1] = coli1 + 1;
					
					/* For Fr1 */
					indexValue[indexi1Fr1] = coli1;
					
					indexValue[indexi1Fr1+1] = coli1 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexi1Fr2] = coli1;
					
					indexValue[indexi1Fr2+1] = coli1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexi1Fr3] = coli1;
					
					indexValue[indexi1Fr3+1] = coli1 + 3;
				}		
				
				if(indexj1Er >= 0){
					/* for Er */
					indexValue[indexj1Er] = colj1;
					
					indexValue[indexj1Er+1] = colj1 + 2;
					
					/* For Fr1 */
					indexValue[indexj1Fr1] = colj1;
					
					
					indexValue[indexj1Fr1+1] = colj1 + 1;
					
					
					/* For Fr2 */
					indexValue[indexj1Fr2] = colj1;
					
					
					indexValue[indexj1Fr2+1] = colj1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexj1Fr3] = colj1;
					
					indexValue[indexj1Fr3+1] = colj1 + 3;
				}		
				
				if(indexk1Er >= 0){
					/* for Er */
					
					indexValue[indexk1Er] = colk1;
					
					
					indexValue[indexk1Er+1] = colk1 + 3;
					
					/* For Fr1 */
					indexValue[indexk1Fr1] = colk1;
					
					indexValue[indexk1Fr1+1] = colk1 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexk1Fr2] = colk1;
					
					indexValue[indexk1Fr2+1] = colk1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexk1Fr3] = colk1;
					
					
					indexValue[indexk1Fr3+1] = colk1 + 3;
				}
				
		

				/* Save the index for furture use */
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][0] = indexk0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][1] = indexk0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][2] = indexk0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][3] = indexk0Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][4] = indexj0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][5] = indexj0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][6] = indexj0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][7] = indexj0Fr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][8] = indexi0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][9] = indexi0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][10] = indexi0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][11] = indexi0Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][12] = indexiEr;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][13] = indexiFr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][14] = indexiFr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][15] = indexiFr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][16] = indexi1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][17] = indexi1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][18] = indexi1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][19] = indexi1Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][20] = indexj1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][21] = indexj1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][22] = indexj1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][23] = indexj1Fr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][24] = indexk1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][25] = indexk1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][26] = indexk1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][27] = indexk1Fr3;			
				
		}/* End loop j */
	}/* End loop k */
}	
						
/****************** End i = is ************/	

	/* Right x boundary */

if((pG->rx1_id < 0) || (pG->rx1_id < ID))
{
	i = ie;
	
	for(k=ks; k<=ke; k++){
		for(j=js; j<=je; j++){

	/* We need to do this for every line */
	indexk0Er = -1; indexk0Fr1 = -1; indexk0Fr2 = -1; indexk0Fr3 = -1;
	indexj0Er = -1; indexj0Fr1 = -1; indexj0Fr2 = -1; indexj0Fr3 = -1;
	indexi0Er = -1; indexi0Fr1 = -1; indexi0Fr2 = -1; indexi0Fr3 = -1;
	indexiEr = -1;  indexiFr1 = -1;  indexiFr2 = -1;  indexiFr3 = -1;
	indexi1Er = -1; indexi1Fr1 = -1; indexi1Fr2 = -1; indexi1Fr3 = -1;
	indexj1Er = -1; indexj1Fr1 = -1; indexj1Fr2 = -1; indexj1Fr3 = -1;
	indexk1Er = -1; indexk1Fr1 = -1; indexk1Fr2 = -1; indexk1Fr3 = -1;

	
				NoEr = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)];
				NoFr1 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+1];
				NoFr2 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+2];
				NoFr3 = ptr[4*(k-ks)*Nx*Ny+4*(j-js)*Nx+4*(i-is)+3];		
						
				
					if(j == js){
						if(k == ks){	
									
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								ie_js_ks_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								ie_js_ks_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								ie_js_ks_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								ie_js_ks_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id < 0))
								ie_js_ks_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->lx3_id >= 0))
								ie_js_ks_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id < 0))
								ie_js_ks_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->lx3_id >= 0))
								ie_js_ks_MPI_MPI_MPI();
						}/* End k == ks */
						else if(k == ke){
											
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								ie_js_ke_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								ie_js_ke_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								ie_js_ke_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								ie_js_ke_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id < 0))
								ie_js_ke_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0) && (pG->rx3_id >= 0))
								ie_js_ke_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id < 0))
								ie_js_ke_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0) && (pG->rx3_id >= 0))
								ie_js_ke_MPI_MPI_MPI();
							
						}/* End k == ke */
						else{/* Now begin k!=ks && k!= ke */
												
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx2_id < 0))
								ie_js_k_phy_phy(k);
							else if((pG->rx1_id < 0) && (pG->lx2_id >= 0))
								ie_js_k_phy_MPI(k);				
							else if((pG->rx1_id >= 0) && (pG->lx2_id < 0))
								ie_js_k_MPI_phy(k);
							else if((pG->rx1_id >= 0) && (pG->lx2_id >= 0))
								ie_js_k_MPI_MPI(k);
							
						}/* End k!=ks and k!=ke */	
					}/* End j == js */
					else if(j == je){
						if(k == ks)						
						{	
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								ie_je_ks_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								ie_je_ks_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								ie_je_ks_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								ie_je_ks_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id < 0))
								ie_je_ks_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->lx3_id >= 0))
								ie_je_ks_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id < 0))
								ie_je_ks_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->lx3_id >= 0))
								ie_je_ks_MPI_MPI_MPI();
						}/* End k == ks */
						else if(k == ke){
												
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								ie_je_ke_phy_phy_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								ie_je_ke_phy_phy_MPI();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								ie_je_ke_phy_MPI_phy();
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								ie_je_ke_phy_MPI_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id < 0))
								ie_je_ke_MPI_phy_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0) && (pG->rx3_id >= 0))
								ie_je_ke_MPI_phy_MPI();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id < 0))
								ie_je_ke_MPI_MPI_phy();
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0) && (pG->rx3_id >= 0))
								ie_je_ke_MPI_MPI_MPI();
							
						}/* End k == ke */
						else{/* begin k!= ks and k!= ke */
												
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx2_id < 0))
								ie_je_k_phy_phy(k);
							else if((pG->rx1_id < 0) && (pG->rx2_id >= 0))
								ie_je_k_phy_MPI(k);				
							else if((pG->rx1_id >= 0) && (pG->rx2_id < 0))
								ie_je_k_MPI_phy(k);
							else if((pG->rx1_id >= 0) && (pG->rx2_id >= 0))
								ie_je_k_MPI_MPI(k);
							
						}/* End k!= ks and k!= ke */
					} /* End j == je */
					else{/* begin j!= js and j!= je */
						if(k == ks){	
						
							
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->lx3_id < 0))
								ie_j_ks_phy_phy(j);
							else if((pG->rx1_id < 0) && (pG->lx3_id >= 0))
								ie_j_ks_phy_MPI(j);
							else if((pG->rx1_id >= 0) && (pG->lx3_id < 0))
								ie_j_ks_MPI_phy(j);
							else if((pG->rx1_id >= 0) && (pG->lx3_id >= 0))
								ie_j_ks_MPI_MPI(j);
						}/* End k == ks */
						else if(k == ke){
													
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0) && (pG->rx3_id < 0))
								ie_j_ke_phy_phy(j);
							else if((pG->rx1_id < 0) && (pG->rx3_id >= 0))
								ie_j_ke_phy_MPI(j);
							else if((pG->rx1_id >= 0) && (pG->rx3_id < 0))
								ie_j_ke_MPI_phy(j);
							else if((pG->rx1_id >= 0) && (pG->rx3_id >= 0))
								ie_j_ke_MPI_MPI(j);
							
						}/* End k == ke */
						else{/* begin k!= ks && k!= ke */
												
							
							/* judge MPI or physical boundary condition */
							if((pG->rx1_id < 0))
								ie_j_k_phy(j,k);
							else
								ie_j_k_MPI(j,k);								
						}/* End k!=ks && k!= ke */
					} /* End j!= js & j != je */
			

				
				
				
				/* There are 7 blocks */
				if(indexk0Er >= 0){
					/* for Er */
					indexValue[indexk0Er] = colk0;
					
					indexValue[indexk0Er+1] = colk0 + 3;
					
					/* For Fr1 */
					indexValue[indexk0Fr1] = colk0;
					
					indexValue[indexk0Fr1+1] = colk0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexk0Fr2] = colk0;
					
					indexValue[indexk0Fr2+1] = colk0 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexk0Fr3] = colk0;
					
					
					indexValue[indexk0Fr3+1] = colk0 + 3;
				}
				
				if(indexj0Er >= 0){
					/* for Er */
					
					indexValue[indexj0Er] = colj0;
					
					
					indexValue[indexj0Er+1] = colj0 + 2;
					
					/* For Fr1 */
					
					indexValue[indexj0Fr1] = colj0;
					
					
					indexValue[indexj0Fr1+1] = colj0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexj0Fr2] = colj0;
					
					
					indexValue[indexj0Fr2+1] = colj0 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexj0Fr3] = colj0;
					
					
					indexValue[indexj0Fr3+1] = colj0 + 3;
				}		
				
				if(indexi0Er >= 0){
					/* for Er */
					
					indexValue[indexi0Er] = coli0;
					
					
					indexValue[indexi0Er+1] = coli0 + 1;
					
					/* For Fr1 */
					
					indexValue[indexi0Fr1] = coli0;
					
					indexValue[indexi0Fr1+1] = coli0 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexi0Fr2] = coli0;
					
					indexValue[indexi0Fr2+1] = coli0 + 2;
					
					
					/* For Fr3 */
					indexValue[indexi0Fr3] = coli0;
					
					indexValue[indexi0Fr3+1] = coli0 + 3;
				}					
				
				
				if(indexiEr >= 0){
					/* for Er */
					
					indexValue[indexiEr] = coli;
					
					indexValue[indexiEr+1] = coli + 1;
					
					
					indexValue[indexiEr+2] = coli + 2;
					
					indexValue[indexiEr+3] = coli + 3;
					
					/* For Fr1 */
					indexValue[indexiFr1] = coli;
					
					
					indexValue[indexiFr1+1] = coli + 1;
					
					
					/* For Fr2 */
					indexValue[indexiFr2] = coli;
					
					indexValue[indexiFr2+1] = coli + 2;
					
					
					/* For Fr3 */
					indexValue[indexiFr3] = coli;
					
					indexValue[indexiFr3+1] = coli + 3;
				}	
				
				if(indexi1Er >= 0){
					/* for Er */
					indexValue[indexi1Er] = coli1;
					
					
					indexValue[indexi1Er+1] = coli1 + 1;
					
					/* For Fr1 */
					indexValue[indexi1Fr1] = coli1;
					
					indexValue[indexi1Fr1+1] = coli1 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexi1Fr2] = coli1;
					
					indexValue[indexi1Fr2+1] = coli1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexi1Fr3] = coli1;
					
					indexValue[indexi1Fr3+1] = coli1 + 3;
				}		
				
				if(indexj1Er >= 0){
					/* for Er */
					indexValue[indexj1Er] = colj1;
					
					indexValue[indexj1Er+1] = colj1 + 2;
					
					/* For Fr1 */
					indexValue[indexj1Fr1] = colj1;
					
					
					indexValue[indexj1Fr1+1] = colj1 + 1;
					
					
					/* For Fr2 */
					indexValue[indexj1Fr2] = colj1;
					
					
					indexValue[indexj1Fr2+1] = colj1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexj1Fr3] = colj1;
					
					indexValue[indexj1Fr3+1] = colj1 + 3;
				}		
				
				if(indexk1Er >= 0){
					/* for Er */
					
					indexValue[indexk1Er] = colk1;
					
					
					indexValue[indexk1Er+1] = colk1 + 3;
					
					/* For Fr1 */
					indexValue[indexk1Fr1] = colk1;
					
					indexValue[indexk1Fr1+1] = colk1 + 1;
					
					
					/* For Fr2 */
					
					indexValue[indexk1Fr2] = colk1;
					
					indexValue[indexk1Fr2+1] = colk1 + 2;
					
					
					/* For Fr3 */
					
					indexValue[indexk1Fr3] = colk1;
					
					
					indexValue[indexk1Fr3+1] = colk1 + 3;
				}
				
		

				/* Save the index for furture use */
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][0] = indexk0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][1] = indexk0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][2] = indexk0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][3] = indexk0Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][4] = indexj0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][5] = indexj0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][6] = indexj0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][7] = indexj0Fr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][8] = indexi0Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][9] = indexi0Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][10] = indexi0Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][11] = indexi0Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][12] = indexiEr;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][13] = indexiFr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][14] = indexiFr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][15] = indexiFr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][16] = indexi1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][17] = indexi1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][18] = indexi1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][19] = indexi1Fr3;			
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][20] = indexj1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][21] = indexj1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][22] = indexj1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][23] = indexj1Fr3;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][24] = indexk1Er;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][25] = indexk1Fr1;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][26] = indexk1Fr2;
				colindex[(k-ks)*Nx*Ny+(j-js)*Nx+(i-is)][27] = indexk1Fr3;			
				
				
				

		}/* End loop j */
	}/* End loop k */
	
}	
	
	return;	
	
}

#endif


/* Function to assign values of the vectors and matrix for different 
 * boundary condition. Judge whether physical boundary or MPI boundary */

/* No boundary case */
void i_j_k(const int i, const int j, const int k)
{
	/* There are seven blocks */
	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexj0Er  = NoEr + 2;
	indexj0Fr1 = NoFr1 + 2;
	indexj0Fr2 = NoFr2 + 2;
	indexj0Fr3 = NoFr3 + 2;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

	indexi0Er  = NoEr + 4;
	indexi0Fr1 = NoFr1 + 4;
	indexi0Fr2 = NoFr2 + 4;
	indexi0Fr3 = NoFr3 + 4;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6;
	indexiFr1 = NoFr1 + 6;
	indexiFr2 = NoFr2 + 6;
	indexiFr3 = NoFr3 + 6;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10;
	indexi1Fr1 = NoFr1 + 8;
	indexi1Fr2 = NoFr2 + 8;
	indexi1Fr3 = NoFr3 + 8;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 12;
	indexj1Fr1 = NoFr1 + 10;
	indexj1Fr2 = NoFr2 + 10;
	indexj1Fr3 = NoFr3 + 10;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 14;
	indexk1Fr1 = NoFr1 + 12;
	indexk1Fr2 = NoFr2 + 12;
	indexk1Fr3 = NoFr3 + 12;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
}

/* one boundary case */
/* boundary along z direction */
void i_j_ks_phy(const int i, const int j)
{
	int k; 
	k = ks;
/* No matter ix3 == 4 or not, it will not affect others */
	
	if(ix3 == 4){
	
		indexk0Er  = NoEr + 14;
		indexk0Fr1 = NoFr1 + 12;
		indexk0Fr2 = NoFr2 + 12;
		indexk0Fr3 = NoFr3 + 12;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}
	/* other x2 boundary condition */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}		


		indexj0Er  = NoEr;
		indexj0Fr1 = NoFr1;
		indexj0Fr2 = NoFr2;
		indexj0Fr3 = NoFr3;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

		indexi0Er  = NoEr + 2;
		indexi0Fr1 = NoFr1 + 2;
		indexi0Fr2 = NoFr2 + 2;
		indexi0Fr3 = NoFr3 + 2;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 4;
		indexiFr1 = NoFr1 + 4;
		indexiFr2 = NoFr2 + 4;
		indexiFr3 = NoFr3 + 4;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8;
		indexi1Fr1 = NoFr1 + 6;
		indexi1Fr2 = NoFr2 + 6;
		indexi1Fr3 = NoFr3 + 6;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 10;
		indexj1Fr1 = NoFr1 + 8;
		indexj1Fr2 = NoFr2 + 8;
		indexj1Fr3 = NoFr3 + 8;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12;
		indexk1Fr1 = NoFr1 + 10;
		indexk1Fr2 = NoFr2 + 10;
		indexk1Fr3 = NoFr3 + 10;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}


void i_j_ks_MPI(const int i, const int j)
{
	int k, shiftz;
	k = ks;

	shiftz = 4 * Nx * Ny * Nz * (lx3 - ID);

	if(lx3 > ID){	
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}
	else{
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */
		/*************************************/
		/* First, MPI part */
		indexk0Er  = NoEr + MPIcount2;
		indexk0Fr1 = NoFr1 + MPIcount2F;
		indexk0Fr2 = NoFr2 + MPIcount2F;
		indexk0Fr3 = NoFr3 + MPIcount2F;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
		/***********************************/
		
		/* non-MPI part */
		indexj0Er  = NoEr + MPIcount1;
		indexj0Fr1 = NoFr1 + MPIcount1;
		indexj0Fr2 = NoFr2 + MPIcount1;
		indexj0Fr3 = NoFr3 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

		indexi0Er  = NoEr + 2 + MPIcount1;
		indexi0Fr1 = NoFr1 + 2 + MPIcount1;
		indexi0Fr2 = NoFr2 + 2 + MPIcount1;
		indexi0Fr3 = NoFr3 + 2 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 4 + MPIcount1;
		indexiFr1 = NoFr1 + 4 + MPIcount1;
		indexiFr2 = NoFr2 + 4 + MPIcount1;
		indexiFr3 = NoFr3 + 4 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8 + MPIcount1;
		indexi1Fr1 = NoFr1 + 6 + MPIcount1;
		indexi1Fr2 = NoFr2 + 6 + MPIcount1;
		indexi1Fr3 = NoFr3 + 6 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 10 + MPIcount1;
		indexj1Fr1 = NoFr1 + 8 + MPIcount1;
		indexj1Fr2 = NoFr2 + 8 + MPIcount1;
		indexj1Fr3 = NoFr3 + 8 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12 + MPIcount1;
		indexk1Fr1 = NoFr1 + 10 + MPIcount1;
		indexk1Fr2 = NoFr2 + 10 + MPIcount1;
		indexk1Fr3 = NoFr3 + 10 + MPIcount1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;	



	return;

}


void i_j_ke_phy(const int i, const int j)
{
	int k, count1;
	k = ke;
	if(ox3 == 4){
				
		count1 = 2;
	}
	else{
			
		count1 = 0;
	}

	
		indexk0Er  = NoEr + count1;
		indexk0Fr1 = NoFr1 + count1;
		indexk0Fr2 = NoFr2 + count1;
		indexk0Fr3 = NoFr3 + count1;
		colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexj0Er  = NoEr + 2 + count1;
		indexj0Fr1 = NoFr1 + 2 + count1;
		indexj0Fr2 = NoFr2 + 2 + count1;
		indexj0Fr3 = NoFr3 + 2 + count1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

		indexi0Er  = NoEr + 4 + count1;
		indexi0Fr1 = NoFr1 + 4 + count1;
		indexi0Fr2 = NoFr2 + 4 + count1;
		indexi0Fr3 = NoFr3 + 4 + count1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 6 + count1;
		indexiFr1 = NoFr1 + 6 + count1;
		indexiFr2 = NoFr2 + 6 + count1;
		indexiFr3 = NoFr3 + 6 + count1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 10 + count1;
		indexi1Fr1 = NoFr1 + 8 + count1;
		indexi1Fr2 = NoFr2 + 8 + count1;
		indexi1Fr3 = NoFr3 + 8 + count1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 12 + count1;
		indexj1Fr1 = NoFr1 + 10 + count1;
		indexj1Fr2 = NoFr2 + 10 + count1;
		indexj1Fr3 = NoFr3 + 10 + count1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* Periodic bounary condition */
	/* other x2 boundary condition */
	else if(ox3 == 1 || ox3 == 5){
					
			theta[6] += theta[14];
			theta[9] -= theta[15];
			
			phi[6] += phi[12];
			phi[7] += phi[13];
	
			psi[6] += psi[12];
			psi[7] += psi[13];
			
			varphi[6] += varphi[12];
			varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
			theta[6] += theta[14];
			theta[9] += theta[15];
			
			phi[6] += phi[12];
			phi[7] += phi[13];
	
			psi[6] += psi[12];
			psi[7] += psi[13];
			
			varphi[6] += varphi[12];
			varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}
	
	return;
}



void i_j_ke_MPI(const int i, const int j)
{
	int k, shiftz;
	
	k = ke;

	shiftz = 4 * Nx * Ny * Nz * (rx3 - ID);

	if(rx3 < ID){
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */

				
		indexk0Er  = NoEr + MPIcount1;
		indexk0Fr1 = NoFr1 + MPIcount1;
		indexk0Fr2 = NoFr2 + MPIcount1;
		indexk0Fr3 = NoFr3 + MPIcount1;
		colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		
		indexj0Er  = NoEr + 2 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

		indexi0Er  = NoEr + 4 + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 6 + MPIcount1;
		indexiFr1 = NoFr1 + 6 + MPIcount1;
		indexiFr2 = NoFr2 + 6 + MPIcount1;
		indexiFr3 = NoFr3 + 6 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 10 + MPIcount1;
		indexi1Fr1 = NoFr1 + 8 + MPIcount1;
		indexi1Fr2 = NoFr2 + 8 + MPIcount1;
		indexi1Fr3 = NoFr3 + 8 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 12 + MPIcount1;
		indexj1Fr1 = NoFr1 + 10 + MPIcount1;
		indexj1Fr2 = NoFr2 + 10 + MPIcount1;
		indexj1Fr3 = NoFr3 + 10 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		/*******************************/
		/* MPI part */
		indexk1Er  = NoEr + MPIcount2;
		indexk1Fr1 = NoFr1 + MPIcount2F;
		indexk1Fr2 = NoFr2 + MPIcount2F;
		indexk1Fr3 = NoFr3 + MPIcount2F;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;	

		/**********************************/

	return;

}


void i_js_k_phy(const int i, const int k)
{
	int j, count1; 
	j = js;

	if(ix2 == 4)	count1 = 2;
	else		count1 = 0;



	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	if(ix2 == 4){
		indexj0Er  = NoEr + 12;
		indexj0Fr1 = NoFr1 + 10;
		indexj0Fr2 = NoFr2 + 10;
		indexj0Fr3 = NoFr3 + 10;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
		
	}
	/* other x2 boundary condition */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}		


		
		indexi0Er  = NoEr + 2;
		indexi0Fr1 = NoFr1 + 2;
		indexi0Fr2 = NoFr2 + 2;
		indexi0Fr3 = NoFr3 + 2;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 4;
		indexiFr1 = NoFr1 + 4;
		indexiFr2 = NoFr2 + 4;
		indexiFr3 = NoFr3 + 4;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8;
		indexi1Fr1 = NoFr1 + 6;
		indexi1Fr2 = NoFr2 + 6;
		indexi1Fr3 = NoFr3 + 6;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 10;
		indexj1Fr1 = NoFr1 + 8;
		indexj1Fr2 = NoFr2 + 8;
		indexj1Fr3 = NoFr3 + 8;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12 + count1;
		indexk1Fr1 = NoFr1 + 10 + count1;
		indexk1Fr2 = NoFr2 + 10 + count1;
		indexk1Fr3 = NoFr3 + 10 + count1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}

void i_js_k_MPI(const int i, const int k)
{
	int j, shifty;
	j = js;

	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 > ID){	
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}
	else{
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */

		
		indexk0Er  = NoEr + MPIcount1;
		indexk0Fr1 = NoFr1 + MPIcount1;
		indexk0Fr2 = NoFr2 + MPIcount1;
		indexk0Fr3 = NoFr3 + MPIcount1;
		colk0 = 4 * (k - ks  -1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/********************************************/	
		/* MPI part */
		indexj0Er  = NoEr + MPIcount2;
		indexj0Fr1 = NoFr1 + MPIcount2F;
		indexj0Fr2 = NoFr2 + MPIcount2F;
		indexj0Fr3 = NoFr3 + MPIcount2F;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/*******************************************/


		indexi0Er  = NoEr + 2 + MPIcount1;
		indexi0Fr1 = NoFr1 + 2 + MPIcount1;
		indexi0Fr2 = NoFr2 + 2 + MPIcount1;
		indexi0Fr3 = NoFr3 + 2 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 4 + MPIcount1;
		indexiFr1 = NoFr1 + 4 + MPIcount1;
		indexiFr2 = NoFr2 + 4 + MPIcount1;
		indexiFr3 = NoFr3 + 4 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8 + MPIcount1;
		indexi1Fr1 = NoFr1 + 6 + MPIcount1;
		indexi1Fr2 = NoFr2 + 6 + MPIcount1;
		indexi1Fr3 = NoFr3 + 6 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 10 + MPIcount1;
		indexj1Fr1 = NoFr1 + 8 + MPIcount1;
		indexj1Fr2 = NoFr2 + 8 + MPIcount1;
		indexj1Fr3 = NoFr3 + 8 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12 + MPIcount1;
		indexk1Fr1 = NoFr1 + 10 + MPIcount1;
		indexk1Fr2 = NoFr2 + 10 + MPIcount1;
		indexk1Fr3 = NoFr3 + 10 + MPIcount1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;	



	return;
}


void i_je_k_phy(const int i, const int k)
{
	int j, count1; 
	j = je;
	if(ox2 == 4){
				
		count1 = 2;
	}
	else{
			
		count1 = 0;
	}


	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	

	indexj0Er  = NoEr + 2 + count1;
	indexj0Fr1 = NoFr1 + 2 + count1;
	indexj0Fr2 = NoFr2 + 2 + count1;
	indexj0Fr3 = NoFr3 + 2 + count1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js  -1) * Nx + 4 * (i - is) + count_Grids;
		
		
	indexi0Er  = NoEr + 4 + count1;
	indexi0Fr1 = NoFr1 + 4 + count1;
	indexi0Fr2 = NoFr2 + 4 + count1;
	indexi0Fr3 = NoFr3 + 4 + count1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + count1;
	indexiFr1 = NoFr1 + 6 + count1;
	indexiFr2 = NoFr2 + 6 + count1;
	indexiFr3 = NoFr3 + 6 + count1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10 + count1;
	indexi1Fr1 = NoFr1 + 8 + count1;
	indexi1Fr2 = NoFr2 + 8 + count1;
	indexi1Fr3 = NoFr3 + 8 + count1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + 2;
		indexj1Fr1 = NoFr1 + 2;
		indexj1Fr2 = NoFr2 + 2;
		indexj1Fr3 = NoFr3 + 2;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}
	/* other x2 boundary condition */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}		

	

	indexk1Er  = NoEr + 12 + count1;
	indexk1Fr1 = NoFr1 + 10 + count1;
	indexk1Fr2 = NoFr2 + 10 + count1;
	indexk1Fr3 = NoFr3 + 10 + count1;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;


}


void i_je_k_MPI(const int i, const int k)
{
	int j, shifty;
	
	j = je;

	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */

		
		indexk0Er  = NoEr + MPIcount1;
		indexk0Fr1 = NoFr1 + MPIcount1;
		indexk0Fr2 = NoFr2 + MPIcount1;
		indexk0Fr3 = NoFr3 + MPIcount1;
		colk0 = 4 * (k - ks  -1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
		indexj0Er  = NoEr + 2 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	


		indexi0Er  = NoEr + 4 + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		indexiEr  = NoEr + 6 + MPIcount1;
		indexiFr1 = NoFr1 + 6 + MPIcount1;
		indexiFr2 = NoFr2 + 6 + MPIcount1;
		indexiFr3 = NoFr3 + 6 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 10 + MPIcount1;
		indexi1Fr1 = NoFr1 + 8 + MPIcount1;
		indexi1Fr2 = NoFr2 + 8 + MPIcount1;
		indexi1Fr3 = NoFr3 + 8 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		/********************************************/	
		/* MPI part */

		indexj1Er  = NoEr +  MPIcount2;
		indexj1Fr1 = NoFr1 + MPIcount2F;
		indexj1Fr2 = NoFr2 + MPIcount2F;
		indexj1Fr3 = NoFr3 + MPIcount2F;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

		/*******************************************/

		indexk1Er  = NoEr + 12 + MPIcount1;
		indexk1Fr1 = NoFr1 + 10 + MPIcount1;
		indexk1Fr2 = NoFr2 + 10 + MPIcount1;
		indexk1Fr3 = NoFr3 + 10 + MPIcount1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;	



	return;

}

void is_j_k_phy(const int j, const int k)
{
	int i, count1; 
	i = is;

	if(ix1 == 4) 	count1 = 2;
	else		count1 = 0;

	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	


	indexj0Er  = NoEr + 2;
	indexj0Fr1 = NoFr1 + 2;
	indexj0Fr2 = NoFr2 + 2;
	indexj0Fr3 = NoFr3 + 2;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix1 == 4){	
		indexi0Er  = NoEr + 10;
		indexi0Fr1 = NoFr1 + 8;
		indexi0Fr2 = NoFr2 + 8;
		indexi0Fr3 = NoFr3 + 8;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	/* other x2 boundary condition */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}		



		indexiEr  = NoEr + 4;
		indexiFr1 = NoFr1 + 4;
		indexiFr2 = NoFr2 + 4;
		indexiFr3 = NoFr3 + 4;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8;
		indexi1Fr1 = NoFr1 + 6;
		indexi1Fr2 = NoFr2 + 6;
		indexi1Fr3 = NoFr3 + 6;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

		indexj1Er  = NoEr + 10 + count1;
		indexj1Fr1 = NoFr1 + 8 + count1;
		indexj1Fr2 = NoFr2 + 8 + count1;
		indexj1Fr3 = NoFr3 + 8 + count1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12 + count1;
		indexk1Fr1 = NoFr1 + 10 + count1;
		indexk1Fr2 = NoFr2 + 10 + count1;
		indexk1Fr3 = NoFr3 + 10 + count1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	

	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	if(col_shear < colk0){
		sheark0 = 2;
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj0){
		sheari0 = 2;
		shearj0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < coli){
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj1){
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colk1){
		sheark1 = 2;
	}
	
	indexk0Er  = NoEr  + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	indexiEr  = NoEr + 4 + sheari;
	indexiFr1 = NoFr1 + 4 + sheari;
	indexiFr2 = NoFr2 + 4 + sheari;
	indexiFr3 = NoFr3 + 4 + sheari;
	
	
	indexi1Er  = NoEr + 8 + sheari;
	indexi1Fr1 = NoFr1 + 6 + sheari;
	indexi1Fr2 = NoFr2 + 6 + sheari;
	indexi1Fr3 = NoFr3 + 6 + sheari;
	
	indexj1Er  = NoEr + 10 + shearj1;
	indexj1Fr1 = NoFr1 + 8 + shearj1;
	indexj1Fr2 = NoFr2 + 8 + shearj1;
	indexj1Fr3 = NoFr3 + 8 + shearj1;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
#endif

	return;

}

void is_j_k_MPI(const int j, const int k)
{
	int i, shiftx;
	i = is;
	

	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */

		
		indexk0Er  = NoEr + MPIcount1;
		indexk0Fr1 = NoFr1 + MPIcount1;
		indexk0Fr2 = NoFr2 + MPIcount1;
		indexk0Fr3 = NoFr3 + MPIcount1;
		colk0 = 4 * (k - ks  -1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
		indexj0Er  = NoEr + 2 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
		/********************************************/	
		/* MPI part */

		indexi0Er  = NoEr + MPIcount2;
		indexi0Fr1 = NoFr1 + MPIcount2F;
		indexi0Fr2 = NoFr2 + MPIcount2F;
		indexi0Fr3 = NoFr3 + MPIcount2F;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;

		/*******************************************/

		indexiEr  = NoEr + 4 + MPIcount1;
		indexiFr1 = NoFr1 + 4 + MPIcount1;
		indexiFr2 = NoFr2 + 4 + MPIcount1;
		indexiFr3 = NoFr3 + 4 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		indexi1Er  = NoEr + 8 + MPIcount1;
		indexi1Fr1 = NoFr1 + 6 + MPIcount1;
		indexi1Fr2 = NoFr2 + 6 + MPIcount1;
		indexi1Fr3 = NoFr3 + 6 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;


		indexj1Er  = NoEr +  10 + MPIcount1;
		indexj1Fr1 = NoFr1 + 8  + MPIcount1;
		indexj1Fr2 = NoFr2 + 8  + MPIcount1;
		indexj1Fr3 = NoFr3 + 8  + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		

		indexk1Er  = NoEr + 12 + MPIcount1;
		indexk1Fr1 = NoFr1 + 10 + MPIcount1;
		indexk1Fr2 = NoFr2 + 10 + MPIcount1;
		indexk1Fr3 = NoFr3 + 10 + MPIcount1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID){
	/* Only use this for left boundary */	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;

	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	if(col_shear < colk0){
		sheark0 = 2;
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj0){
		sheari0 = 2;
		shearj0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < coli){
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj1){
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colk1){
		sheark1 = 2;
	}
	
	indexk0Er  = NoEr  + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	indexiEr  = NoEr + 4 + sheari;
	indexiFr1 = NoFr1 + 4 + sheari;
	indexiFr2 = NoFr2 + 4 + sheari;
	indexiFr3 = NoFr3 + 4 + sheari;
	
	
	indexi1Er  = NoEr + 8 + sheari;
	indexi1Fr1 = NoFr1 + 6 + sheari;
	indexi1Fr2 = NoFr2 + 6 + sheari;
	indexi1Fr3 = NoFr3 + 6 + sheari;
	
	indexj1Er  = NoEr + 10 + shearj1;
	indexj1Fr1 = NoFr1 + 8 + shearj1;
	indexj1Fr2 = NoFr2 + 8 + shearj1;
	indexj1Fr3 = NoFr3 + 8 + shearj1;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
	} /* End lx1 > ID */	
#endif


	return;
}

void ie_j_k_phy(const int j, const int k)
{
	int i, count1; 
	i = ie;

	if(ox1 == 4){
				
		count1 = 2;
	}
	else{
			
		count1 = 0;
	}


	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	indexj0Er  = NoEr + 2;
	indexj0Fr1 = NoFr1 + 2;
	indexj0Fr2 = NoFr2 + 2;
	indexj0Fr3 = NoFr3 + 2;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + count1;
	indexi0Fr1 = NoFr1 + 4 + count1;
	indexi0Fr2 = NoFr2 + 4 + count1;
	indexi0Fr3 = NoFr3 + 4 + count1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	
		
	indexiEr  = NoEr + 6 + count1;
	indexiFr1 = NoFr1 + 6 + count1;
	indexiFr2 = NoFr2 + 6 + count1;
	indexiFr3 = NoFr3 + 6 + count1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){	
		indexi1Er  = NoEr + 4;
		indexi1Fr1 = NoFr1 + 4;
		indexi1Fr2 = NoFr2 + 4;
		indexi1Fr3 = NoFr3 + 4;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}
		/* other x2 boundary condition */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

		indexj1Er  = NoEr + 10 + count1;
		indexj1Fr1 = NoFr1 + 8 + count1;
		indexj1Fr2 = NoFr2 + 8 + count1;
		indexj1Fr3 = NoFr3 + 8 + count1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		indexk1Er  = NoEr + 12 + count1;
		indexk1Fr1 = NoFr1 + 10 + count1;
		indexk1Fr2 = NoFr2 + 10 + count1;
		indexk1Fr3 = NoFr3 + 10 + count1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	if(col_shear < colk0){
		sheark0 = 2;
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj0){
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < coli){
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj1){
		
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colk1){
		sheark1 = 2;
	}
	
	indexk0Er  = NoEr  + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	indexi0Er  = NoEr + 4 + sheari;
	indexi0Fr1 = NoFr1 + 4 + sheari;
	indexi0Fr2 = NoFr2 + 4 + sheari;
	indexi0Fr3 = NoFr3 + 4 + sheari;
	
	indexiEr  = NoEr + 6 + sheari;
	indexiFr1 = NoFr1 + 6 + sheari;
	indexiFr2 = NoFr2 + 6 + sheari;
	indexiFr3 = NoFr3 + 6 + sheari;
	
	indexj1Er  = NoEr + 10 + shearj1;
	indexj1Fr1 = NoFr1 + 8 + shearj1;
	indexj1Fr2 = NoFr2 + 8 + shearj1;
	indexj1Fr3 = NoFr3 + 8 + shearj1;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;	
	
#endif

	return;

}

void ie_j_k_MPI(const int j, const int k)
{
	int i, shiftx;
	
	i = ie;
	

	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 14;
		MPIcount2F = 12;
	}

/* For MPI boundary, no difference for periodic and non-periodic boundary condition */

		
		indexk0Er  = NoEr + MPIcount1;
		indexk0Fr1 = NoFr1 + MPIcount1;
		indexk0Fr2 = NoFr2 + MPIcount1;
		indexk0Fr3 = NoFr3 + MPIcount1;
		colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
		indexj0Er  = NoEr + 2 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

		indexi0Er  = NoEr + 4 + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

		
		indexiEr  = NoEr + 6 + MPIcount1;
		indexiFr1 = NoFr1 + 6 + MPIcount1;
		indexiFr2 = NoFr2 + 6 + MPIcount1;
		indexiFr3 = NoFr3 + 6 + MPIcount1;
		coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

		
		/********************************************/	
		/* MPI part */

		indexi1Er  = NoEr + MPIcount2;
		indexi1Fr1 = NoFr1 + MPIcount2F;
		indexi1Fr2 = NoFr2 + MPIcount2F;
		indexi1Fr3 = NoFr3 + MPIcount2F;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
		/*******************************************/



		indexj1Er  = NoEr +  10 + MPIcount1;
		indexj1Fr1 = NoFr1 + 8  + MPIcount1;
		indexj1Fr2 = NoFr2 + 8  + MPIcount1;
		indexj1Fr3 = NoFr3 + 8  + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		

		indexk1Er  = NoEr + 12 + MPIcount1;
		indexk1Fr1 = NoFr1 + 10 + MPIcount1;
		indexk1Fr2 = NoFr2 + 10 + MPIcount1;
		indexk1Fr3 = NoFr3 + 10 + MPIcount1;
		colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	if(col_shear < colk0){
		sheark0 = 2;
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj0){
		shearj0 = 2;
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < coli){
		sheari0 = 2;
		sheari = 2;
		sheari1 = 2;
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colj1){
		
		shearj1 = 2;
		sheark1 = 2;
	}
	else if(col_shear < colk1){
		sheark1 = 2;
	}
	
	indexk0Er  = NoEr  + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	indexi0Er  = NoEr + 4 + sheari;
	indexi0Fr1 = NoFr1 + 4 + sheari;
	indexi0Fr2 = NoFr2 + 4 + sheari;
	indexi0Fr3 = NoFr3 + 4 + sheari;
	
	indexiEr  = NoEr + 6 + sheari;
	indexiFr1 = NoFr1 + 6 + sheari;
	indexiFr2 = NoFr2 + 6 + sheari;
	indexiFr3 = NoFr3 + 6 + sheari;
	
	indexj1Er  = NoEr + 10 + shearj1;
	indexj1Fr1 = NoFr1 + 8 + shearj1;
	indexj1Fr2 = NoFr2 + 8 + shearj1;
	indexj1Fr3 = NoFr3 + 8 + shearj1;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
}
	
#endif


	return;

}


/* two boundaries case */
/* for i first, i, js, ks */
void i_js_ks_phy_phy(const int i)
{
	int j, k, count1y;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + count1y;
		indexk0Fr1 = NoFr1 + 10 + count1y;
		indexk0Fr2 = NoFr2 + 10 + count1y;
		indexk0Fr3 = NoFr3 + 10 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10;
		indexj0Fr1 = NoFr1 + 8;
		indexj0Fr2 = NoFr2 + 8;
		indexj0Fr3 = NoFr3 + 8;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	indexi0Er  = NoEr;
	indexi0Fr1 = NoFr1;
	indexi0Fr2 = NoFr2;
	indexi0Fr3 = NoFr3;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 2;
	indexiFr1 = NoFr1 + 2;
	indexiFr2 = NoFr2 + 2;
	indexiFr3 = NoFr3 + 2;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6;
	indexi1Fr1 = NoFr1 + 4;
	indexi1Fr2 = NoFr2 + 4;
	indexi1Fr3 = NoFr3 + 4;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8;
	indexj1Fr1 = NoFr1 + 6;
	indexj1Fr2 = NoFr2 + 6;
	indexj1Fr3 = NoFr3 + 6;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + count1y;
	indexk1Fr1 = NoFr1 + 8 + count1y;
	indexk1Fr2 = NoFr2 + 8 + count1y;
	indexk1Fr3 = NoFr3 + 8 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;
}

void i_js_ks_MPI_phy(const int i)
{

	int j, k, count1z, shifty, MPIcount2F;
	j = js;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + MPIcount1;
		indexk0Fr1 = NoFr1 + 10 + MPIcount1;
		indexk0Fr2 = NoFr2 + 10 + MPIcount1;
		indexk0Fr3 = NoFr3 + 10 + MPIcount1;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	indexi0Er  = NoEr + MPIcount1;
	indexi0Fr1 = NoFr1 + MPIcount1;
	indexi0Fr2 = NoFr2 + MPIcount1;
	indexi0Fr3 = NoFr3 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}

void i_js_ks_phy_MPI(const int i)
{

	int j, k, count1y, shiftz, MPIcount2F;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/



	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	indexi0Er  = NoEr + MPIcount1;
	indexi0Fr1 = NoFr1 + MPIcount1;
	indexi0Fr2 = NoFr2 + MPIcount1;
	indexi0Fr3 = NoFr3 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}

void i_js_ks_MPI_MPI(const int i)
{
	
	int j, k, MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz, shiftz, shifty;
	j = js;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	

	indexi0Er  = NoEr + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1z + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;


}


/******************/
/* begin i, js, ke */

void i_js_ke_phy_phy(const int i)
{
	int j, k, count1z;
	j = js;
	k = ke;

		

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

		/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 12 + count1z;
		indexj0Fr1 = NoFr1 + 10 + count1z;
		indexj0Fr2 = NoFr2 + 10 + count1z;
		indexj0Fr3 = NoFr3 + 10 + count1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	indexi0Er  = NoEr + 2 + count1z;
	indexi0Fr1 = NoFr1 + 2 + count1z;
	indexi0Fr2 = NoFr2 + 2 + count1z;
	indexi0Fr3 = NoFr3 + 2 + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + count1z;
	indexiFr1 = NoFr1 + 4 + count1z;
	indexiFr2 = NoFr2 + 4 + count1z;
	indexiFr3 = NoFr3 + 4 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + count1z;
	indexi1Fr1 = NoFr1 + 6 + count1z; 
	indexi1Fr2 = NoFr2 + 6 + count1z; 
	indexi1Fr3 = NoFr3 + 6 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + count1z;
	indexj1Fr1 = NoFr1 + 8 + count1z;
	indexj1Fr2 = NoFr2 + 8 + count1z;
	indexj1Fr3 = NoFr3 + 8 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
	indexk1Er  = NoEr;
	indexk1Fr1 = NoFr1;
	indexk1Fr2 = NoFr2;
	indexk1Fr3 = NoFr3;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}	
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	


	return;
}

void i_js_ke_MPI_phy(const int i)
{

	int j, k, count1z, shifty;
	j = js;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	indexi0Er  = NoEr + 2  + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	return;

}

void i_js_ke_phy_MPI(const int i)
{

	int j, k, count1y, shiftz, MPIcount2F;
	j = js;
	k = ke;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix2 == 4){
		indexj0Er  = NoEr + 12 + MPIcount1;
		indexj0Fr1 = NoFr1 + 10 + MPIcount1;
		indexj0Fr2 = NoFr2 + 10 + MPIcount1;
		indexj0Fr3 = NoFr3 + 10 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	indexi0Er  = NoEr + 2 + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/


	return;

}

void i_js_ke_MPI_MPI(const int i)
{
	
	int j, k, MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz, shiftz, shifty;
	j = js;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}



	/* There are seven blocks */	
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1y;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	

	indexi0Er  = NoEr + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1z + MPIcount1y;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1y;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1y;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1y;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/

	return;
}

/* Now begin i, je, ks part */

void i_je_ks_phy_phy(const int i)
{
	int j, k, count1y;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + count1y;
		indexk0Fr1 = NoFr1 + 10 + count1y;
		indexk0Fr2 = NoFr2 + 10 + count1y;
		indexk0Fr3 = NoFr3 + 10 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	
	indexj0Er  = NoEr + count1y;
	indexj0Fr1 = NoFr1 + count1y;
	indexj0Fr2 = NoFr2 + count1y;
	indexj0Fr3 = NoFr3 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 2 + count1y;
	indexi0Fr1 = NoFr1 + 2 + count1y;
	indexi0Fr2 = NoFr2 + 2 + count1y;
	indexi0Fr3 = NoFr3 + 2 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + count1y;
	indexiFr1 = NoFr1 + 4 + count1y;
	indexiFr2 = NoFr2 + 4 + count1y;
	indexiFr3 = NoFr3 + 4 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + count1y;
	indexi1Fr1 = NoFr1 + 6 + count1y;
	indexi1Fr2 = NoFr2 + 6 + count1y;
	indexi1Fr3 = NoFr3 + 6 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr;
		indexj1Fr1 = NoFr1;
		indexj1Fr2 = NoFr2;
		indexj1Fr3 = NoFr3;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 10 + count1y;
	indexk1Fr1 = NoFr1 + 8 + count1y;
	indexk1Fr2 = NoFr2 + 8 + count1y;
	indexk1Fr3 = NoFr3 + 8 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;
}

void i_je_ks_MPI_phy(const int i)
{

	int j, k, count1z, shifty, MPIcount2F;
	j = je;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
   	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + MPIcount1;
		indexk0Fr1 = NoFr1 + 10 + MPIcount1;
		indexk0Fr2 = NoFr2 + 10 + MPIcount1;
		indexk0Fr3 = NoFr3 + 10 + MPIcount1;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 2 + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;


	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 10 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}

void i_je_ks_phy_MPI(const int i)
{

	int j, k, count1y, shiftz;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



	/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/

	
	indexj0Er  = NoEr + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	indexi0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + MPIcount1;
		indexj1Fr1 = NoFr1 + MPIcount1;
		indexj1Fr2 = NoFr2 + MPIcount1;
		indexj1Fr3 = NoFr3 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;

}

void i_je_ks_MPI_MPI(const int i)
{
	
	int j, k, MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz, shiftz, shifty;
	j = je;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	indexj0Er  = NoEr + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexi0Er  = NoEr + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	
	
	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 10 + MPIcount1z + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	return;


}


/* Now begin i, je, ke */


void i_je_ke_phy_phy(const int i)
{
	int j, k, count1y, count1z;
	j = je;
	k = ke;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox3 == 4) count1z = 2;
	else count1z = 0;


		/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary */
	
	indexj0Er  = NoEr + 2 + count1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + count1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + count1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + count1y + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + count1y + count1z;
	indexi0Fr1 = NoFr1 + 4 + count1y + count1z;
	indexi0Fr2 = NoFr2 + 4 + count1y + count1z;
	indexi0Fr3 = NoFr3 + 4 + count1y + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + count1y + count1z;
	indexiFr1 = NoFr1 + 6 + count1y + count1z;
	indexiFr2 = NoFr2 + 6 + count1y + count1z;
	indexiFr3 = NoFr3 + 6 + count1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10 + count1y + count1z;
	indexi1Fr1 = NoFr1 + 8 + count1y + count1z;
	indexi1Fr2 = NoFr2 + 8 + count1y + count1z;
	indexi1Fr3 = NoFr3 + 8 + count1y + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + count1z;
		indexj1Fr1 = NoFr1 + 2 + count1z;
		indexj1Fr2 = NoFr2 + 2 + count1z;
		indexj1Fr3 = NoFr3 + 2 + count1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	return;
}

void i_je_ke_MPI_phy(const int i)
{

	int j, k, count1z, shifty;
	j = je;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1 + count1z; 
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;


	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	return;

}

void i_je_ke_phy_MPI(const int i)
{

	int j, k, count1y, shiftz;
	j = je;
	k = ke;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1y;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + MPIcount1;
		indexj1Fr1 = NoFr1 +  2 + MPIcount1;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids  + shiftz;
	/******************************/


	return;

}

void i_je_ke_MPI_MPI(const int i)
{
	
	int j, k, MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz, shiftz, shifty;
	j = je;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1y;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 +  MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexi0Er  = NoEr + 4 + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 10 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	
	
	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	/***************************/
	/* for y boundary, MPI part */
	
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids  + shiftz;

	/******************************/


	return;


}

/* Now for j part */
/* begin, is, j, ks */

void is_j_ks_phy_phy(const int j)
{
	int i, k, count1x;
	i = is;
	k = ks;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + count1x;
		indexk0Fr1 = NoFr1 + 10 + count1x;
		indexk0Fr2 = NoFr2 + 10 + count1x;
		indexk0Fr3 = NoFr3 + 10 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	
	indexj0Er  = NoEr;
	indexj0Fr1 = NoFr1;
	indexj0Fr2 = NoFr2;
	indexj0Fr3 = NoFr3;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix1 == 4){
		indexi0Er  = NoEr + 8;
		indexi0Fr1 = NoFr1 + 6;
		indexi0Fr2 = NoFr2 + 6;
		indexi0Fr3 = NoFr3 + 6;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2;
	indexiFr1 = NoFr1 + 2;
	indexiFr2 = NoFr2 + 2;
	indexiFr3 = NoFr3 + 2;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6;
	indexi1Fr1 = NoFr1 + 4;
	indexi1Fr2 = NoFr2 + 4;
	indexi1Fr3 = NoFr3 + 4;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + count1x;
	indexj1Fr1 = NoFr1 + 6 + count1x;
	indexj1Fr2 = NoFr2 + 6 + count1x;
	indexj1Fr3 = NoFr3 + 6 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + count1x;
	indexk1Fr1 = NoFr1 + 8 + count1x;
	indexk1Fr2 = NoFr2 + 8 + count1x;
	indexk1Fr3 = NoFr3 + 8 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0;
	indexj0Fr1 = NoFr1 + shearj0;
	indexj0Fr2 = NoFr2 + shearj0;
	indexj0Fr3 = NoFr3 + shearj0;
	
	indexiEr  = NoEr + 2 + sheari;
	indexiFr1 = NoFr1 + 2 + sheari;
	indexiFr2 = NoFr2 + 2 + sheari;
	indexiFr3 = NoFr3 + 2 + sheari;
	
	
	indexi1Er  = NoEr + 6 + sheari;
	indexi1Fr1 = NoFr1 + 4 + sheari;
	indexi1Fr2 = NoFr2 + 4 + sheari;
	indexi1Fr3 = NoFr3 + 4 + sheari;
	
	indexj1Er  = NoEr + 8 + shearj1;
	indexj1Fr1 = NoFr1 + 6 + shearj1;
	indexj1Fr2 = NoFr2 + 6 + shearj1;
	indexj1Fr3 = NoFr3 + 6 + shearj1;
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	
	
#endif

	return;
}



void is_j_ks_MPI_phy(const int j)
{

	int i, k, count1z, shiftx;
	i = is;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + MPIcount1;
		indexk0Fr1 = NoFr1 + 10 + MPIcount1;
		indexk0Fr2 = NoFr2 + 10 + MPIcount1;
		indexk0Fr3 = NoFr3 + 10 + MPIcount1;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	

	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}

	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	

	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0;
	indexj0Fr1 = NoFr1 + shearj0;
	indexj0Fr2 = NoFr2 + shearj0;
	indexj0Fr3 = NoFr3 + shearj0;
	
	indexiEr  = NoEr + 2 + sheari;
	indexiFr1 = NoFr1 + 2 + sheari;
	indexiFr2 = NoFr2 + 2 + sheari;
	indexiFr3 = NoFr3 + 2 + sheari;
	
	
	indexi1Er  = NoEr + 6 + sheari;
	indexi1Fr1 = NoFr1 + 4 + sheari;
	indexi1Fr2 = NoFr2 + 4 + sheari;
	indexi1Fr3 = NoFr3 + 4 + sheari;
	
	indexj1Er  = NoEr + 8 + shearj1;
	indexj1Fr1 = NoFr1 + 6 + shearj1;
	indexj1Fr2 = NoFr2 + 6 + shearj1;
	indexj1Fr3 = NoFr3 + 6 + shearj1;
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	

	

	
}	
#endif
	

	return;

}

void is_j_ks_phy_MPI(const int j)
{

	int i, k, count1x, shiftz;
	i = is;
	k = ks;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/



	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	if(ix1 == 4){
	indexi0Er  = NoEr + 8 + MPIcount1;
	indexi0Fr1 = NoFr1 + 6 + MPIcount1;
	indexi0Fr2 = NoFr2 + 6 + MPIcount1;
	indexi0Fr3 = NoFr3 + 6 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(lx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount2 + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount2F + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount2F + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount2F + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0 + MPIcount1;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1;

	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
	
#endif
	

	return;

}

void is_j_ks_MPI_MPI(const int j)
{
	
	int i, k, MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz, shiftz, shiftx;
	i = is;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	

	indexj0Er  = NoEr + MPIcount1x + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for y boundary, MPI part */
	indexi0Er  = NoEr + MPIcount1z + MPIcount2x;
	indexi0Fr1 = NoFr1 + MPIcount1z + MPIcount2Fx;
	indexi0Fr2 = NoFr2 + MPIcount1z + MPIcount2Fx;
	indexi0Fr3 = NoFr3 + MPIcount1z + MPIcount2Fx;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1x;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1x;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1x;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1z + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)	
{	
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount2z + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount2Fz + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount2Fz + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount2Fz + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1z;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1z;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
	

} 	
#endif

	return;


}



/*=======================*/

/* begin, is, j, ke */

void is_j_ke_phy_phy(const int j)
{
	int i, k, count1x, count1z;
	i = is;
	k = ke;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

		/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	indexj0Er  = NoEr + 2 + count1z;
	indexj0Fr1 = NoFr1 + 2 + count1z;
	indexj0Fr2 = NoFr2 + 2 + count1z;
	indexj0Fr3 = NoFr3 + 2 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	/***************************/
	/* for x boundary */
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + count1z;
		indexi0Fr1 = NoFr1 + 8 + count1z;
		indexi0Fr2 = NoFr2 + 8 + count1z;
		indexi0Fr3 = NoFr3 + 8 + count1z;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 4 + count1z;
	indexiFr1 = NoFr1 + 4 + count1z;
	indexiFr2 = NoFr2 + 4 + count1z;
	indexiFr3 = NoFr3 + 4 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + count1z;
	indexi1Fr1 = NoFr1 + 6 + count1z; 
	indexi1Fr2 = NoFr2 + 6 + count1z; 
	indexi1Fr3 = NoFr3 + 6 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + count1z + count1x;
	indexj1Fr1 = NoFr1 + 8 + count1z + count1x;
	indexj1Fr2 = NoFr2 + 8 + count1z + count1x;
	indexj1Fr3 = NoFr3 + 8 + count1z + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
	indexk1Er  = NoEr;
	indexk1Fr1 = NoFr1;
	indexk1Fr2 = NoFr2;
	indexk1Fr3 = NoFr3;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}	
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	indexj0Er  = NoEr  + 2 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + count1z;
	
	indexiEr  = NoEr + 4 + sheari + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + count1z;
	
	
	indexi1Er  = NoEr + 8 + sheari + count1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + count1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + count1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + count1z;
	
	indexj1Er  = NoEr + 10 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 8 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 8 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 8 + shearj1 + count1z;
	
	
	
#endif
	
	

	return;
}

void is_j_ke_MPI_phy(const int j)
{

	int i, k, count1z, shiftx;
	i = is;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	


	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
#ifdef SHEARING_BOX
	if(lx1 > ID)
	{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
		
		indexj0Er  = NoEr  + 2 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 2 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 2 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 2 + shearj0 + count1z;
		
		indexiEr  = NoEr + 4 + sheari + count1z;
		indexiFr1 = NoFr1 + 4 + sheari + count1z;
		indexiFr2 = NoFr2 + 4 + sheari + count1z;
		indexiFr3 = NoFr3 + 4 + sheari + count1z;
		
		
		indexi1Er  = NoEr + 8 + sheari + count1z;
		indexi1Fr1 = NoFr1 + 6 + sheari + count1z;
		indexi1Fr2 = NoFr2 + 6 + sheari + count1z;
		indexi1Fr3 = NoFr3 + 6 + sheari + count1z;
		
		indexj1Er  = NoEr + 10 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 8 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 8 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 8 + shearj1 + count1z;
	
		
	}	
#endif

	
	return;

}

void is_j_ke_phy_MPI(const int j)
{

	int i, k, count1x, shiftz;
	i = is;
	k = ke;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + MPIcount1;
		indexi0Fr1 = NoFr1 + 8 + MPIcount1;
		indexi0Fr2 = NoFr2 + 8 + MPIcount1;
		indexi0Fr3 = NoFr3 + 8 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for z boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(rx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1 + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1 + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1 + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1 + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1;
	
	indexj1Er  = NoEr + 10 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 8 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 8 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 8 + shearj1 + MPIcount1;

	
	indexk1Er  = NoEr  + MPIcount2 + sheark1;
	indexk1Fr1 = NoFr1 + MPIcount2F + sheark1;
	indexk1Fr2 = NoFr2 + MPIcount2F + sheark1;
	indexk1Fr3 = NoFr3 + MPIcount2F + sheark1;
	
	
	
#endif	


	return;

}

void is_j_ke_MPI_MPI(const int j)
{
	
	int i, k, MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz, shiftz, shiftx;
	i = is;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



	/* There are seven blocks */	
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	

	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount1z + MPIcount2x;
	indexi0Fr1 = NoFr1 + MPIcount1z + MPIcount2Fx;
	indexi0Fr2 = NoFr2 + MPIcount1z + MPIcount2Fx;
	indexi0Fr3 = NoFr3 + MPIcount1z + MPIcount2Fx;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1z + MPIcount1x;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1x;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1x;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 10 + MPIcount1z + MPIcount1x;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1z + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1z + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1z + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1z + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1z;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 10 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 8 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 8 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 8 + shearj1 + MPIcount1z;
	
	
	indexk1Er  = NoEr  + MPIcount2z + sheark1;
	indexk1Fr1 = NoFr1 + MPIcount2Fz + sheark1;
	indexk1Fr2 = NoFr2 + MPIcount2Fz + sheark1;
	indexk1Fr3 = NoFr3 + MPIcount2Fz + sheark1;
	
}	
	
#endif		
	

	return;
}



/* begin ie, j, ks */


void ie_j_ks_phy_phy(const int j)
{
	int i, k, count1x;
	i = ie;
	k = ks;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + count1x;
		indexk0Fr1 = NoFr1 + 10 + count1x;
		indexk0Fr2 = NoFr2 + 10 + count1x;
		indexk0Fr3 = NoFr3 + 10 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	
	indexj0Er  = NoEr;
	indexj0Fr1 = NoFr1;
	indexj0Fr2 = NoFr2;
	indexj0Fr3 = NoFr3;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 2 + count1x;
	indexi0Fr1 = NoFr1 + 2 + count1x;
	indexi0Fr2 = NoFr2 + 2 + count1x;
	indexi0Fr3 = NoFr3 + 2 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + count1x;
	indexiFr1 = NoFr1 + 4 + count1x;
	indexiFr2 = NoFr2 + 4 + count1x;
	indexiFr3 = NoFr3 + 4 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 2;
		indexi1Fr1 = NoFr1 + 2;
		indexi1Fr2 = NoFr2 + 2;
		indexi1Fr3 = NoFr3 + 2;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	
	indexj1Er  = NoEr + 8 + count1x;
	indexj1Fr1 = NoFr1 + 6 + count1x;
	indexj1Fr2 = NoFr2 + 6 + count1x;
	indexj1Fr3 = NoFr3 + 6 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexk1Er  = NoEr + 10 + count1x;
	indexk1Fr1 = NoFr1 + 8 + count1x;
	indexk1Fr2 = NoFr2 + 8 + count1x;
	indexk1Fr3 = NoFr3 + 8 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0;
	indexj0Fr1 = NoFr1 + shearj0;
	indexj0Fr2 = NoFr2 + shearj0;
	indexj0Fr3 = NoFr3 + shearj0;
	
	
	indexi0Er  = NoEr + 2 + sheari;
	indexi0Fr1 = NoFr1 + 2 + sheari;
	indexi0Fr2 = NoFr2 + 2 + sheari;
	indexi0Fr3 = NoFr3 + 2 + sheari;
	
	indexiEr  = NoEr + 4 + sheari;
	indexiFr1 = NoFr1 + 4 + sheari;
	indexiFr2 = NoFr2 + 4 + sheari;
	indexiFr3 = NoFr3 + 4 + sheari;
	
	indexj1Er  = NoEr + 8 + shearj1;
	indexj1Fr1 = NoFr1 + 6 + shearj1;
	indexj1Fr2 = NoFr2 + 6 + shearj1;
	indexj1Fr3 = NoFr3 + 6 + shearj1;
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
#endif
		

	return;
}

void ie_j_ks_MPI_phy(const int j)
{

	int i, k, count1z, shiftx;
	i = ie;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
   	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 12 + MPIcount1;
		indexk0Fr1 = NoFr1 + 10 + MPIcount1;
		indexk0Fr2 = NoFr2 + 10 + MPIcount1;
		indexk0Fr3 = NoFr3 + 10 + MPIcount1;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 2 + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;

	/******************************/
	

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	

	indexk1Er  = NoEr + 10 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 12 + sheark0;
			indexk0Fr1 = NoFr1 + 10 + sheark0;
			indexk0Fr2 = NoFr2 + 10 + sheark0;
			indexk0Fr3 = NoFr3 + 10 + sheark0;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexj0Er  = NoEr  + shearj0;
		indexj0Fr1 = NoFr1 + shearj0;
		indexj0Fr2 = NoFr2 + shearj0;
		indexj0Fr3 = NoFr3 + shearj0;
		
		
		indexi0Er  = NoEr + 2 + sheari;
		indexi0Fr1 = NoFr1 + 2 + sheari;
		indexi0Fr2 = NoFr2 + 2 + sheari;
		indexi0Fr3 = NoFr3 + 2 + sheari;
	
		indexiEr  = NoEr + 4 + sheari;
		indexiFr1 = NoFr1 + 4 + sheari;
		indexiFr2 = NoFr2 + 4 + sheari;
		indexiFr3 = NoFr3 + 4 + sheari;
		
		indexj1Er  = NoEr + 8 + shearj1;
		indexj1Fr1 = NoFr1 + 6 + shearj1;
		indexj1Fr2 = NoFr2 + 6 + shearj1;
		indexj1Fr3 = NoFr3 + 6 + shearj1;
		
		indexk1Er = NoEr + 10 + sheark1;
		indexk1Fr1 = NoFr1 + 8 + sheark1;
		indexk1Fr2 = NoFr2 + 8 + sheark1;
		indexk1Fr3 = NoFr3 + 8 + sheark1;
	
}	
#endif


	return;

}

void ie_j_ks_phy_MPI(const int j)
{

	int i, k, count1x, shiftz;
	i = ie;
	k = ks;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/

	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	indexi0Er  = NoEr + 2 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + MPIcount1 + 2;
		indexi1Fr1 = NoFr1 + MPIcount1 + 2;
		indexi1Fr2 = NoFr2 + MPIcount1 + 2;
		indexi1Fr3 = NoFr3 + MPIcount1 + 2;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	
	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	
	
#ifdef SHEARING_BOX

		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(lx3 > ID){
		
			MPIcount2 = 12;
			MPIcount2F = 10;
		}
	
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexk0Er  = NoEr  + MPIcount2 + sheark0;
		indexk0Fr1 = NoFr1 + MPIcount2F + sheark0;
		indexk0Fr2 = NoFr2 + MPIcount2F + sheark0;
		indexk0Fr3 = NoFr3 + MPIcount2F + sheark0;
	
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexj0Er  = NoEr  + shearj0 + MPIcount1;
		indexj0Fr1 = NoFr1 + shearj0 + MPIcount1;
		indexj0Fr2 = NoFr2 + shearj0 + MPIcount1; 
		indexj0Fr3 = NoFr3 + shearj0 + MPIcount1;
		
		
		indexi0Er  = NoEr + 2 + sheari + MPIcount1;
		indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1;
		indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1;
		indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1;
		
		indexiEr  = NoEr + 4 + sheari + MPIcount1;
		indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
		indexiFr2 = NoFr2 + 4 + sheari + MPIcount1; 
		indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
		
		indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
		indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
		indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1;
		indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
		
		indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
		
#endif
	

	return;

}

void ie_j_ks_MPI_MPI(const int j)
{
	
	int i, k, MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz, shiftz, shiftx;
	i = ie;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	indexj0Er  = NoEr + MPIcount1x + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexi0Er  = NoEr + 2 + MPIcount1z + MPIcount1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount1z + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	
	/******************************/
	
	
	

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	indexk1Er  = NoEr + 10 + MPIcount1z + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er  = NoEr  + MPIcount2z + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount2Fz + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount2Fz + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount2Fz + sheark0;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1z; 
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1z;
	
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z; 
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1z;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
}
	
#endif

 return;


}


/* begin ie, j, ke part */



void ie_j_ke_phy_phy(const int j)
{
	int i, k, count1x, count1z;
	i = ie;
	k = ke;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox3 == 4) count1z = 2;
	else count1z = 0;


		/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + count1z;
	indexj0Fr1 = NoFr1 + 2 + count1z;
	indexj0Fr2 = NoFr2 + 2 + count1z;
	indexj0Fr3 = NoFr3 + 2 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + count1x + count1z;
	indexi0Fr1 = NoFr1 + 4 + count1x + count1z;
	indexi0Fr2 = NoFr2 + 4 + count1x + count1z;
	indexi0Fr3 = NoFr3 + 4 + count1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + count1x + count1z;
	indexiFr1 = NoFr1 + 6 + count1x + count1z;
	indexiFr2 = NoFr2 + 6 + count1x + count1z;
	indexiFr3 = NoFr3 + 6 + count1x + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + count1z + 4;
		indexi1Fr1 = NoFr1 + count1z + 4;
		indexi1Fr2 = NoFr2 + count1z + 4;
		indexi1Fr3 = NoFr3 + count1z + 4;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

		}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	
	indexj1Er  = NoEr + 10 + count1z + count1x;
	indexj1Fr1 = NoFr1 + 8 + count1z + count1x;
	indexj1Fr2 = NoFr2 + 8 + count1z + count1x;
	indexj1Fr3 = NoFr3 + 8 + count1z + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + 2 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + count1z;
	
	
	indexi0Er  = NoEr + 4 + sheari + count1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + count1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + count1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + count1z;
	
	indexiEr  = NoEr + 6 + sheari + count1z;
	indexiFr1 = NoFr1 + 6 + sheari + count1z;
	indexiFr2 = NoFr2 + 6 + sheari + count1z;
	indexiFr3 = NoFr3 + 6 + sheari + count1z;
	
	indexj1Er  = NoEr + 10 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 8 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 8 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 8 + shearj1 + count1z;
	
#endif

	return;
}

void ie_j_ke_MPI_phy(const int j)
{

	int i, k, count1z, shiftx;
	i = ie;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1z;
		MPIcount2F = 10 + count1z;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1 + count1z; 
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	/***************************/
	/* for x boundary, MPI part */


	indexi1Er  = NoEr + MPIcount2; 
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;

	
	/******************************/


	indexj1Er  = NoEr + 10 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
#ifdef SHEARING_BOX
	if(rx1 < ID)
	{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
	
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexj0Er  = NoEr  + 2 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 2 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 2 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 2 + shearj0 + count1z;
		
		
		indexi0Er  = NoEr + 4 + sheari + count1z;
		indexi0Fr1 = NoFr1 + 4 + sheari + count1z;
		indexi0Fr2 = NoFr2 + 4 + sheari + count1z;
		indexi0Fr3 = NoFr3 + 4 + sheari + count1z;
		
		indexiEr  = NoEr + 6 + sheari + count1z;
		indexiFr1 = NoFr1 + 6 + sheari + count1z;
		indexiFr2 = NoFr2 + 6 + sheari + count1z;
		indexiFr3 = NoFr3 + 6 + sheari + count1z;
		
		indexj1Er  = NoEr + 10 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 8 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 8 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 8 + shearj1 + count1z;
	
		
	}	
#endif
	

	return;

}

void ie_j_ke_phy_MPI(const int j)
{

	int i, k, count1x, shiftz;
	i = ie;
	k = ke;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;



	indexiEr  = NoEr + 6 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	

	if(ox1 == 4){
		indexi1Er  = NoEr + MPIcount1 + 4;
		indexi1Fr1 = NoFr1 + MPIcount1 + 4;
		indexi1Fr2 = NoFr2 + MPIcount1 + 4;
		indexi1Fr3 = NoFr3 + MPIcount1 + 4;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for x direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	indexj1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids  + shiftz;
	/******************************/
	
	
	
#ifdef SHEARING_BOX

		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(rx3 > ID){
		
			MPIcount2 = 12;
			MPIcount2F = 10;
		}
	
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexk0Er = NoEr + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1;
		
		
		indexi0Er  = NoEr + 4 + sheari + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1;
		
		indexiEr  = NoEr + 6 + sheari + MPIcount1;
		indexiFr1 = NoFr1 + 6 + sheari + MPIcount1;
		indexiFr2 = NoFr2 + 6 + sheari + MPIcount1;
		indexiFr3 = NoFr3 + 6 + sheari + MPIcount1;
		
		indexj1Er  = NoEr + 10 + shearj1 + MPIcount1;
		indexj1Fr1 = NoFr1 + 8 + shearj1 + MPIcount1;
		indexj1Fr2 = NoFr2 + 8 + shearj1 + MPIcount1;
		indexj1Fr3 = NoFr3 + 8 + shearj1 + MPIcount1;
	
		indexk1Er = NoEr + sheark1 + MPIcount2;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount2F;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount2F;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount2F;
	
		
		

#endif	
		

	return;

}

void ie_j_ke_MPI_MPI(const int j)
{
	
	int i, k, MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz, shiftz, shiftx;
	i = ie;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexi0Er  = NoEr + 4 + MPIcount1z + MPIcount1x;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1x;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1x;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1z + MPIcount1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1z + MPIcount1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1z + MPIcount1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1z + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount1z + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/
	
	
	

	indexj1Er  = NoEr + 10 + MPIcount1z + MPIcount1x;
	indexj1Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1x;
	indexj1Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	


	/***************************/
	/* for y boundary, MPI part */
	
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids  + shiftz;

	/******************************/
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)	
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1z;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1z;
	
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 10 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 8 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 8 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 8 + shearj1 + MPIcount1z;
	
	indexk1Er = NoEr + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
}		
	
#endif		


	return;


}


/* Now for k */

/* begin is, js, k */


void is_js_k_phy_phy(const int k)
{
	int i, j, count1x, count1y;
	i = is;
	j = js;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;


	/* There are seven blocks */
	
	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + count1x;
		indexj0Fr1 = NoFr1 + 8 + count1x;
		indexj0Fr2 = NoFr2+ 8 + count1x;
		indexj0Fr3 = NoFr3 + 8 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for y direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	

	if(ix1 == 4){
		indexi0Er  = NoEr + 8;
		indexi0Fr1 = NoFr1 + 6;
		indexi0Fr2 = NoFr2 + 6;
		indexi0Fr3 = NoFr3 + 6;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2;
	indexiFr1 = NoFr1 + 2;
	indexiFr2 = NoFr2 + 2;
	indexiFr3 = NoFr3 + 2;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6;
	indexi1Fr1 = NoFr1 + 4;
	indexi1Fr2 = NoFr2 + 4;
	indexi1Fr3 = NoFr3 + 4;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + count1x;
	indexj1Fr1 = NoFr1 + 6 + count1x;
	indexj1Fr2 = NoFr2 + 6 + count1x;
	indexj1Fr3 = NoFr3 + 6 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 8 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 8 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 8 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	/* y direction is periodic */
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexk0Er = NoEr + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	
	indexiEr  = NoEr + 2 + sheari;
	indexiFr1 = NoFr1 + 2 + sheari;
	indexiFr2 = NoFr2 + 2 + sheari;
	indexiFr3 = NoFr3 + 2 + sheari;
	
	
	indexi1Er  = NoEr + 6 + sheari;
	indexi1Fr1 = NoFr1 + 4 + sheari;
	indexi1Fr2 = NoFr2 + 4 + sheari;
	indexi1Fr3 = NoFr3 + 4 + sheari;
	
	indexj1Er  = NoEr + 8 + shearj1;
	indexj1Fr1 = NoFr1 + 6 + shearj1;
	indexj1Fr2 = NoFr2 + 6 + shearj1;
	indexj1Fr3 = NoFr3 + 6 + shearj1;
	
	indexj0Er  = NoEr  + 10 + shearj0;
	indexj0Fr1 = NoFr1 + 8 + shearj0;
	indexj0Fr2 = NoFr2 + 8 + shearj0;
	indexj0Fr3 = NoFr3 + 8 + shearj0;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
#endif

	return;
}



void is_js_k_MPI_phy(const int k)
{

	int i, j, count1y, shiftx;
	j = js;
	i = is;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



		/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
	{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		/* y direction is periodic */
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexk0Er = NoEr + sheark0;
		indexk0Fr1 = NoFr1 + sheark0;
		indexk0Fr2 = NoFr2 + sheark0;
		indexk0Fr3 = NoFr3 + sheark0;
		
		
		indexiEr  = NoEr + 2 + sheari;
		indexiFr1 = NoFr1 + 2 + sheari;
		indexiFr2 = NoFr2 + 2 + sheari;
		indexiFr3 = NoFr3 + 2 + sheari;
		
		
		indexi1Er  = NoEr + 6 + sheari;
		indexi1Fr1 = NoFr1 + 4 + sheari;
		indexi1Fr2 = NoFr2 + 4 + sheari;
		indexi1Fr3 = NoFr3 + 4 + sheari;
		
		indexj1Er  = NoEr + 8 + shearj1;
		indexj1Fr1 = NoFr1 + 6 + shearj1;
		indexj1Fr2 = NoFr2 + 6 + shearj1;
		indexj1Fr3 = NoFr3 + 6 + shearj1;
		
		indexj0Er  = NoEr  + 10 + shearj0;
		indexj0Fr1 = NoFr1 + 8 + shearj0;
		indexj0Fr2 = NoFr2 + 8 + shearj0;
		indexj0Fr3 = NoFr3 + 8 + shearj0;
		
		indexk1Er = NoEr + 12 + sheark1;
		indexk1Fr1 = NoFr1 + 10 + sheark1;
		indexk1Fr2 = NoFr2 + 10 + sheark1;
		indexk1Fr3 = NoFr3 + 10 + sheark1;
		
		
	}	
#endif
	

	return;

}

void is_js_k_phy_MPI(const int k)
{

	int i, j, count1x, shifty;
	i = is;
	j = js;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary, MPI part */
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(lx2 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1 + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1 + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1 + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1 + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + MPIcount2 + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
	
#endif

	return;

}

void is_js_k_MPI_MPI(const int k)
{
	
	int i, j, MPIcount1x, MPIcount1y, MPIcount2x, MPIcount2y, MPIcount2Fx, MPIcount2Fy, shifty, shiftx;
	i = is;
	j = js;

		
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 14;
		MPIcount2Fy = 12;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1x + MPIcount1y;
	indexk0Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
	indexk0Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
	indexk0Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;



	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids  + shifty;

	/******************************/
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount1y + MPIcount2x;
	indexi0Fr1 = NoFr1 + MPIcount1y + MPIcount2Fx;
	indexi0Fr2 = NoFr2 + MPIcount1y + MPIcount2Fx;
	indexi0Fr3 = NoFr3 + MPIcount1y + MPIcount2Fx;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(lx2 > ID){
		
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1y + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1y + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1y + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1y + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1y;
	
}	
	
#endif	
	
	return;


}

/* begin is, je, k */


void is_je_k_phy_phy(const int k)
{
	int i, j, count1y, count1x;
	i = is;
	j = je;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

		/* There are seven blocks */
	
	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	indexj0Er  = NoEr + 2 + count1y;
	indexj0Fr1 = NoFr1 + 2 + count1y;
	indexj0Fr2 = NoFr2 + 2 + count1y;
	indexj0Fr3 = NoFr3 + 2 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	/***************************/
	/* for x boundary */
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + count1y;
		indexi0Fr1 = NoFr1 + 8 + count1y;
		indexi0Fr2 = NoFr2 + 8 + count1y;
		indexi0Fr3 = NoFr3 + 8 + count1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 4 + count1y;
	indexiFr1 = NoFr1 + 4 + count1y;
	indexiFr2 = NoFr2 + 4 + count1y;
	indexiFr3 = NoFr3 + 4 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + count1y;
	indexi1Fr1 = NoFr1 + 6 + count1y; 
	indexi1Fr2 = NoFr2 + 6 + count1y; 
	indexi1Fr3 = NoFr3 + 6 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){

		indexj1Er  = NoEr + 2;
		indexj1Fr1 = NoFr1 + 2;
		indexj1Fr2 = NoFr2 + 2;
		indexj1Fr3 = NoFr3 + 2;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}	
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	



	
	indexk1Er  = NoEr + 10 + count1y + count1x;
	indexk1Fr1 = NoFr1 + 8 + count1y + count1x;
	indexk1Fr2 = NoFr2 + 8 + count1y + count1x;
	indexk1Fr3 = NoFr3 + 8 + count1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	/* y direction is periodic */
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexk0Er = NoEr + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexj1Er  = NoEr + 2 + shearj1;
	indexj1Fr1 = NoFr1 + 2 + shearj1;
	indexj1Fr2 = NoFr2 + 2 + shearj1;
	indexj1Fr3 = NoFr3 + 2 + shearj1;
	
	indexj0Er  = NoEr  + 4 + shearj0;
	indexj0Fr1 = NoFr1 + 4 + shearj0;
	indexj0Fr2 = NoFr2 + 4 + shearj0;
	indexj0Fr3 = NoFr3 + 4 + shearj0;
	
	
	indexiEr  = NoEr + 6 + sheari;
	indexiFr1 = NoFr1 + 6 + sheari;
	indexiFr2 = NoFr2 + 6 + sheari;
	indexiFr3 = NoFr3 + 6 + sheari;
	
	
	indexi1Er  = NoEr + 10 + sheari;
	indexi1Fr1 = NoFr1 + 8 + sheari;
	indexi1Fr2 = NoFr2 + 8 + sheari;
	indexi1Fr3 = NoFr3 + 8 + sheari;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
#endif


	return;
}

void is_je_k_MPI_phy(const int k)
{

	int i, j, count1y, shiftx;
	i = is;
	j = je;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + MPIcount1 + 2;
		indexj1Fr1 = NoFr1 + MPIcount1 + 2;
		indexj1Fr2 = NoFr2 + MPIcount1 + 2;
		indexj1Fr3 = NoFr3 + MPIcount1 + 2;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	



	
	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
	{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		/* y direction is periodic */
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexk0Er = NoEr + sheark0;
		indexk0Fr1 = NoFr1 + sheark0;
		indexk0Fr2 = NoFr2 + sheark0;
		indexk0Fr3 = NoFr3 + sheark0;
		
		indexj1Er  = NoEr + 2 + shearj1;
		indexj1Fr1 = NoFr1 + 2 + shearj1;
		indexj1Fr2 = NoFr2 + 2 + shearj1;
		indexj1Fr3 = NoFr3 + 2 + shearj1;
		
		indexj0Er  = NoEr  + 4 + shearj0;
		indexj0Fr1 = NoFr1 + 4 + shearj0;
		indexj0Fr2 = NoFr2 + 4 + shearj0;
		indexj0Fr3 = NoFr3 + 4 + shearj0;
		
		
		indexiEr  = NoEr + 6 + sheari;
		indexiFr1 = NoFr1 + 6 + sheari;
		indexiFr2 = NoFr2 + 6 + sheari;
		indexiFr3 = NoFr3 + 6 + sheari;
		
		
		indexi1Er  = NoEr + 10 + sheari;
		indexi1Fr1 = NoFr1 + 8 + sheari;
		indexi1Fr2 = NoFr2 + 8 + sheari;
		indexi1Fr3 = NoFr3 + 8 + sheari;
		
		indexk1Er = NoEr + 12 + sheark1;
		indexk1Fr1 = NoFr1 + 10 + sheark1;
		indexk1Fr2 = NoFr2 + 10 + sheark1;
		indexk1Fr3 = NoFr3 + 10 + sheark1;
		
		
	}	
#endif
	
	
	return;

}

void is_je_k_phy_MPI(const int k)
{

	int i, j, count1x, shifty;
	i = is;
	j = je;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + MPIcount1;
		indexi0Fr1 = NoFr1 + 8 + MPIcount1;
		indexi0Fr2 = NoFr2 + 8 + MPIcount1;
		indexi0Fr3 = NoFr3 + 8 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	

	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x; 
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(rx2 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1 + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1 + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1 + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1 + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	

	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1;
	

	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	indexj1Er  = NoEr  + MPIcount2 + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;
	
	
	
#endif	

	return;

}

void is_je_k_MPI_MPI(const int k)
{
	
	int i, j, MPIcount1x, MPIcount1y, MPIcount2x, MPIcount2y, MPIcount2Fx, MPIcount2Fy, shifty, shiftx;
	i = is;
	j = je;

		
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 14;
		MPIcount2Fy = 12;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



	/* There are seven blocks */	
	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	

	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount1y + MPIcount2x;
	indexi0Fr1 = NoFr1 + MPIcount1y + MPIcount2Fx;
	indexi0Fr2 = NoFr2 + MPIcount1y + MPIcount2Fx;
	indexi0Fr3 = NoFr3 + MPIcount1y + MPIcount2Fx;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	

	indexk1Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)	
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	if(rx2 > ID){
		
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}
	
	
	
	
	indexk0Er  = NoEr  + MPIcount1y + sheark0;
	indexk0Fr1 = NoFr1 + MPIcount1y + sheark0;
	indexk0Fr2 = NoFr2 + MPIcount1y + sheark0;
	indexk0Fr3 = NoFr3 + MPIcount1y + sheark0;
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1y;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1y;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1y;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1y;
	
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1y;
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
}	
	
#endif	
	

	return;
}

/* begin ie, js, k */


void ie_js_k_phy_phy(const int k)
{
	int i, j, count1x, count1y;
	i = ie;
	j = js;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;


	/* There are seven blocks */
	
	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + count1x;
		indexj0Fr1 = NoFr1 + 8 + count1x;
		indexj0Fr2 = NoFr2 + 8 + count1x;
		indexj0Fr3 = NoFr3 + 8 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	indexi0Er  = NoEr + 2 + count1x;
	indexi0Fr1 = NoFr1 + 2 + count1x;
	indexi0Fr2 = NoFr2 + 2 + count1x;
	indexi0Fr3 = NoFr3 + 2 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + count1x;
	indexiFr1 = NoFr1 + 4 + count1x;
	indexiFr2 = NoFr2 + 4 + count1x;
	indexiFr3 = NoFr3 + 4 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 2;
		indexi1Fr1 = NoFr1 + 2;
		indexi1Fr2 = NoFr2 + 2;
		indexi1Fr3 = NoFr3 + 2;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	
	indexj1Er  = NoEr + 8 + count1x;
	indexj1Fr1 = NoFr1 + 6 + count1x;
	indexj1Fr2 = NoFr2 + 6 + count1x;
	indexj1Fr3 = NoFr3 + 6 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexk1Er  = NoEr + 10 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 8 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 8 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 8 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexk0Er = NoEr + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexi0Er  = NoEr + 2 + sheari;
	indexi0Fr1 = NoFr1 + 2 + sheari;
	indexi0Fr2 = NoFr2 + 2 + sheari;
	indexi0Fr3 = NoFr3 + 2 + sheari;
	
	indexiEr  = NoEr + 4 + sheari;
	indexiFr1 = NoFr1 + 4 + sheari;
	indexiFr2 = NoFr2 + 4 + sheari;
	indexiFr3 = NoFr3 + 4 + sheari;
	
	indexj1Er  = NoEr + 8 + shearj1;
	indexj1Fr1 = NoFr1 + 6 + shearj1;
	indexj1Fr2 = NoFr2 + 6 + shearj1;
	indexj1Fr3 = NoFr3 + 6 + shearj1;
	
	indexj0Er  = NoEr  + 10 + shearj0;
	indexj0Fr1 = NoFr1 + 8 + shearj0;
	indexj0Fr2 = NoFr2 + 8 + shearj0;
	indexj0Fr3 = NoFr3 + 8 + shearj0;
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
#endif

	return;
}

void ie_js_k_MPI_phy(const int k)
{

	int i, j, count1y, shiftx;
	i = ie;
	j = js;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
   	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



		/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	indexi0Er  = NoEr + 2 + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1;
	indexiFr1 = NoFr1 + 4 + MPIcount1;
	indexiFr2 = NoFr2 + 4 + MPIcount1;
	indexiFr3 = NoFr3 + 4 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;

	/******************************/
	

	indexj1Er  = NoEr + 8 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	

	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
	{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		indexk0Er = NoEr + sheark0;
		indexk0Fr1 = NoFr1 + sheark0;
		indexk0Fr2 = NoFr2 + sheark0;
		indexk0Fr3 = NoFr3 + sheark0;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
				
		
		indexi0Er  = NoEr + 2 + sheari;
		indexi0Fr1 = NoFr1 + 2 + sheari;
		indexi0Fr2 = NoFr2 + 2 + sheari;
		indexi0Fr3 = NoFr3 + 2 + sheari;
		
		indexiEr  = NoEr + 4 + sheari;
		indexiFr1 = NoFr1 + 4 + sheari;
		indexiFr2 = NoFr2 + 4 + sheari;
		indexiFr3 = NoFr3 + 4 + sheari;
		
		indexj1Er  = NoEr + 8 + shearj1;
		indexj1Fr1 = NoFr1 + 6 + shearj1;
		indexj1Fr2 = NoFr2 + 6 + shearj1;
		indexj1Fr3 = NoFr3 + 6 + shearj1;
		
		indexj0Er  = NoEr  + 10 + shearj0;
		indexj0Fr1 = NoFr1 + 8 + shearj0;
		indexj0Fr2 = NoFr2 + 8 + shearj0;
		indexj0Fr3 = NoFr3 + 8 + shearj0;
		
		indexk1Er = NoEr + 12 + sheark1;
		indexk1Fr1 = NoFr1 + 10 + sheark1;
		indexk1Fr2 = NoFr2 + 10 + sheark1;
		indexk1Fr3 = NoFr3 + 10 + sheark1;
		
	}	
#endif
	

	return;

}

void ie_js_k_phy_MPI(const int k)
{

	int i, j, count1x, shifty;
	i = ie;
	j = js;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */


	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	

	indexi0Er  = NoEr + 2 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + MPIcount1 + 2;
		indexi1Fr1 = NoFr1 + MPIcount1 + 2;
		indexi1Fr2 = NoFr2 + MPIcount1 + 2;
		indexi1Fr3 = NoFr3 + MPIcount1 + 2;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	
	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	
#ifdef SHEARING_BOX
	
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
	
		if(lx2 > ID){
		
			MPIcount2 = 12;
			MPIcount2F = 10;
		}
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		indexk0Er = NoEr + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		
		indexi0Er  = NoEr + 2 + sheari + MPIcount1;
		indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1;
		indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1;
		indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1;
		
		indexiEr  = NoEr + 4 + sheari + MPIcount1;
		indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
		indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
		indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
		
		indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
		indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
		indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1; 
		indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
		
		indexj0Er  = NoEr  + MPIcount2 + shearj0;
		indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
		indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
		indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
		
		indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;

	
#endif	
	

	return;

}

void ie_js_k_MPI_MPI(const int k)
{
	
	int i, j, MPIcount1x, MPIcount1y, MPIcount2x, MPIcount2y, MPIcount2Fx, MPIcount2Fy, shifty, shiftx;
	i = ie;
	j = js;

		
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 14;
		MPIcount2Fy = 12;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary, MPI part */
	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	

	indexi0Er  = NoEr + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 4 + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	
	/******************************/
	
	
	

	indexj1Er  = NoEr + 8 + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	indexk1Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(lx2 > ID){
		
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1y;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y; 
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1y;
	
}
		
#endif	

	return;
}


/* begin ie, je, k */


void ie_je_k_phy_phy(const int k)
{
	int i, j, count1x, count1y;
	i = ie;
	j = je;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox2 == 4) count1y = 2;
	else count1y = 0;


		/* There are seven blocks */
	
	indexk0Er  = NoEr;
	indexk0Fr1 = NoFr1;
	indexk0Fr2 = NoFr2;
	indexk0Fr3 = NoFr3;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + count1y;
	indexj0Fr1 = NoFr1 + 2 + count1y;
	indexj0Fr2 = NoFr2 + 2 + count1y;
	indexj0Fr3 = NoFr3 + 2 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + count1x + count1y;
	indexi0Fr1 = NoFr1 + 4 + count1x + count1y;
	indexi0Fr2 = NoFr2 + 4 + count1x + count1y;
	indexi0Fr3 = NoFr3 + 4 + count1x + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + count1x + count1y;
	indexiFr1 = NoFr1 + 6 + count1x + count1y;
	indexiFr2 = NoFr2 + 6 + count1x + count1y;
	indexiFr3 = NoFr3 + 6 + count1x + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 4 + count1y;
		indexi1Fr1 = NoFr1 + 4 + count1y;
		indexi1Fr2 = NoFr2 + 4 + count1y;
		indexi1Fr3 = NoFr3 + 4 + count1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

		}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	if(ox2 == 4){
		indexj1Er  = NoEr + 2;
		indexj1Fr1 = NoFr1 + 2; 
		indexj1Fr2 = NoFr2 + 2;
		indexj1Fr3 = NoFr3 + 2;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	
	indexk1Er  = NoEr + 10 + count1y + count1x;
	indexk1Fr1 = NoFr1 + 8 + count1y + count1x;
	indexk1Fr2 = NoFr2 + 8 + count1y + count1x;
	indexk1Fr3 = NoFr3 + 8 + count1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexk0Er = NoEr + sheark0;
	indexk0Fr1 = NoFr1 + sheark0;
	indexk0Fr2 = NoFr2 + sheark0;
	indexk0Fr3 = NoFr3 + sheark0;
	
	indexj1Er  = NoEr + 2 + shearj1;
	indexj1Fr1 = NoFr1 + 2 + shearj1;
	indexj1Fr2 = NoFr2 + 2 + shearj1;
	indexj1Fr3 = NoFr3 + 2 + shearj1;
	
	indexj0Er  = NoEr  + 4 + shearj0;
	indexj0Fr1 = NoFr1 + 4 + shearj0;
	indexj0Fr2 = NoFr2 + 4 + shearj0;
	indexj0Fr3 = NoFr3 + 4 + shearj0;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexi0Er  = NoEr + 6 + sheari;
	indexi0Fr1 = NoFr1 + 6 + sheari;
	indexi0Fr2 = NoFr2 + 6 + sheari;
	indexi0Fr3 = NoFr3 + 6 + sheari;
	
	indexiEr  = NoEr + 8 + sheari;
	indexiFr1 = NoFr1 + 8 + sheari;
	indexiFr2 = NoFr2 + 8 + sheari;
	indexiFr3 = NoFr3 + 8 + sheari;
	
	
	indexk1Er = NoEr + 12 + sheark1;
	indexk1Fr1 = NoFr1 + 10 + sheark1;
	indexk1Fr2 = NoFr2 + 10 + sheark1;
	indexk1Fr3 = NoFr3 + 10 + sheark1;
	
	
#endif
	
	

	return;
}


void ie_je_k_MPI_phy(const int k)
{

	int i, j, count1y, shiftx;
	i = ie;
	j = je;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1y;
		MPIcount2F = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1y;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1 + count1y; 
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	/***************************/
	/* for x boundary, MPI part */


	indexi1Er  = NoEr + MPIcount2; 
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;

	
	/******************************/

	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + MPIcount1;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	
	indexk1Er  = NoEr + 10 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		indexk0Er = NoEr + sheark0;
		indexk0Fr1 = NoFr1 + sheark0;
		indexk0Fr2 = NoFr2 + sheark0;
		indexk0Fr3 = NoFr3 + sheark0;
	
		indexj1Er  = NoEr + 2 + shearj1;
		indexj1Fr1 = NoFr1 + 2 + shearj1;
		indexj1Fr2 = NoFr2 + 2 + shearj1;
		indexj1Fr3 = NoFr3 + 2 + shearj1;
	
		indexj0Er  = NoEr  + 4 + shearj0;
		indexj0Fr1 = NoFr1 + 4 + shearj0;
		indexj0Fr2 = NoFr2 + 4 + shearj0;
		indexj0Fr3 = NoFr3 + 4 + shearj0;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		
		indexi0Er  = NoEr + 6 + sheari;
		indexi0Fr1 = NoFr1 + 6 + sheari;
		indexi0Fr2 = NoFr2 + 6 + sheari;
		indexi0Fr3 = NoFr3 + 6 + sheari;
		
		indexiEr  = NoEr + 8 + sheari;
		indexiFr1 = NoFr1 + 8 + sheari;
		indexiFr2 = NoFr2 + 8 + sheari;
		indexiFr3 = NoFr3 + 8 + sheari;
		
	
		indexk1Er = NoEr + 12 + sheark1;
		indexk1Fr1 = NoFr1 + 10 + sheark1;
		indexk1Fr2 = NoFr2 + 10 + sheark1;
		indexk1Fr3 = NoFr3 + 10 + sheark1;
		
	}	
#endif	

	return;

}

void ie_je_k_phy_MPI(const int k)
{

	int i, j, count1x, shifty;
	i = ie;
	j = je;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 12 + count1x;
		MPIcount2F = 10 + count1x;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	indexi0Er  = NoEr + 4 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;



	indexiEr  = NoEr + 6 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	

	if(ox1 == 4){
		indexi1Er  = NoEr + MPIcount1 + 4;
		indexi1Fr1 = NoFr1 + MPIcount1 + 4;
		indexi1Fr2 = NoFr2 + MPIcount1 + 4;
		indexi1Fr3 = NoFr3 + MPIcount1 + 4;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for x direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for z boundary, MPI part */
	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids  + shifty;
	/******************************/
	
	indexk1Er  = NoEr + 10 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX

		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
	
		if(rx2 > ID){
		
			MPIcount2 = 12;
			MPIcount2F = 10;
		}
	
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		indexk0Er = NoEr + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
		

		
		indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1;
		indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1;
		indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1;
		indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1;
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		
		indexi0Er  = NoEr + 4 + sheari + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1;
		
		indexiEr  = NoEr + 6 + sheari + MPIcount1;
		indexiFr1 = NoFr1 + 6 + sheari + MPIcount1;
		indexiFr2 = NoFr2 + 6 + sheari + MPIcount1;
		indexiFr3 = NoFr3 + 6 + sheari + MPIcount1;
			
		
		indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
		indexj1Er  = NoEr + MPIcount2 + shearj1;
		indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
		indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
		indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;

#endif	
	


	return;

}

void ie_je_k_MPI_MPI(const int k)
{
	
	int i, j, MPIcount1x, MPIcount1y, MPIcount2x, MPIcount2y, MPIcount2Fx, MPIcount2Fy, shifty, shiftx;
	i = ie;
	j = je;

		
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 14;
		MPIcount2Fy = 12;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 12;
		MPIcount2Fx = 10;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexi0Er  = NoEr + 4 + MPIcount1y + MPIcount1x;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1y + MPIcount1x;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1y + MPIcount1x;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1y + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;

	indexiEr  = NoEr + 6 + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */

	indexi1Er  = NoEr + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/
	
	
	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	
	
	indexk1Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids ;

	
	
#ifdef SHEARING_BOX
if(rx1 < ID)	
{	
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(rx2 > ID){
		
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y;
	
	
	
	indexj0Er  = NoEr  + 2 + shearj0 + MPIcount1y;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1y;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1y;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1y;
	
	
	indexj1Er  = NoEr + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
}
	
#endif	
	
	

	return;


}


/*================================================================*/
/*================================================================*/

/* Three boundaries case */
/* begin is, js, ks */

void is_js_ks_phy_phy_phy()
{
	int i, j, k, count1y, count1x;
	i = is;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + count1y + count1x;
		indexk0Fr1 = NoFr1 + 8 + count1y + count1x;
		indexk0Fr2 = NoFr2 + 8 + count1y + count1x;
		indexk0Fr3 = NoFr3 + 8 + count1y + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + count1x;
		indexj0Fr1 = NoFr1 + 6 + count1x;
		indexj0Fr2 = NoFr2 + 6 + count1x;
		indexj0Fr3 = NoFr3 + 6 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for x boundary */
	if(ix1 == 4){
	indexi0Er  = NoEr + 6;
	indexi0Fr1 = NoFr1 + 4;
	indexi0Fr2 = NoFr2 + 4;
	indexi0Fr3 = NoFr3 + 4;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr;
	indexiFr1 = NoFr1;
	indexiFr2 = NoFr2;
	indexiFr3 = NoFr3;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4;
	indexi1Fr1 = NoFr1 + 2;
	indexi1Fr2 = NoFr2 + 2;
	indexi1Fr3 = NoFr3 + 2;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + count1x;
	indexj1Fr1 = NoFr1 + 4 + count1x;
	indexj1Fr2 = NoFr2 + 4 + count1x;
	indexj1Fr3 = NoFr3 + 4 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 6 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 6 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 6 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari;
	indexiFr1 = NoFr1 + sheari;
	indexiFr2 = NoFr2 + sheari;
	indexiFr3 = NoFr3 + sheari;
	
	
	indexi1Er  = NoEr + 4 + sheari;
	indexi1Fr1 = NoFr1 + 2 + sheari;
	indexi1Fr2 = NoFr2 + 2 + sheari;
	indexi1Fr3 = NoFr3 + 2 + sheari;
	
	indexj1Er  = NoEr + 6 + shearj1;
	indexj1Fr1 = NoFr1 + 4 + shearj1;
	indexj1Fr2 = NoFr2 + 4+ shearj1;
	indexj1Fr3 = NoFr3 + 4 + shearj1;
	
	indexj0Er  = NoEr  + 8 + shearj0;
	indexj0Fr1 = NoFr1 + 6 + shearj0;
	indexj0Fr2 = NoFr2 + 6 + shearj0;
	indexj0Fr3 = NoFr3 + 6 + shearj0;
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	
	
#endif
	

	return;
}

void is_js_ks_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = is;
	j = js;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;
	if(ix2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1y;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + MPIcount1;
	indexiFr1 = NoFr1 + MPIcount1;
	indexiFr2 = NoFr2 + MPIcount1;
	indexiFr3 = NoFr3 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 12 + sheark0;
			indexk0Fr1 = NoFr1 + 10 + sheark0;
			indexk0Fr2 = NoFr2 + 10 + sheark0;
			indexk0Fr3 = NoFr3 + 10 + sheark0;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
				
		indexiEr  = NoEr + sheari;
		indexiFr1 = NoFr1 + sheari;
		indexiFr2 = NoFr2 + sheari;
		indexiFr3 = NoFr3 + sheari;
		
		
		indexi1Er  = NoEr + 4 + sheari;
		indexi1Fr1 = NoFr1 + 2 + sheari;
		indexi1Fr2 = NoFr2 + 2 + sheari;
		indexi1Fr3 = NoFr3 + 2 + sheari;
		
		indexj1Er  = NoEr + 6 + shearj1;
		indexj1Fr1 = NoFr1 + 4 + shearj1;
		indexj1Fr2 = NoFr2 + 4+ shearj1;
		indexj1Fr3 = NoFr3 + 4 + shearj1;
		
		indexj0Er  = NoEr  + 8 + shearj0;
		indexj0Fr1 = NoFr1 + 6 + shearj0;
		indexj0Fr2 = NoFr2 + 6 + shearj0;
		indexj0Fr3 = NoFr3 + 6 + shearj0;
	
		indexk1Er = NoEr + 10 + sheark1;
		indexk1Fr1 = NoFr1 + 8 + sheark1;
		indexk1Fr2 = NoFr2 + 8 + sheark1;
		indexk1Fr3 = NoFr3 + 8 + sheark1;
		
	
		
	}	
#endif
	

	return;

}

void is_js_ks_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = is;
	j = js;
	k = ks;

	
	if(ix1 == 4) count1x = 2;
	else	count1x = 0;
	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for x boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	
	/******************************/
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 6 + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	
	

	indexiEr  = NoEr + MPIcount1;
	indexiFr1 = NoFr1 + MPIcount1;
	indexiFr2 = NoFr2 + MPIcount1;
	indexiFr3 = NoFr3 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
		
		
		if(lx2 > ID){
			
			MPIcount2 = 10 + count1z;
			MPIcount2F = 8 + count1z;
		}
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1;
			indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1;
			indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1;
			indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		indexj0Er  = NoEr  + MPIcount2 + shearj0;
		indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
		indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
		indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
	
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexiEr  = NoEr + sheari + MPIcount1;
		indexiFr1 = NoFr1 + sheari + MPIcount1;
		indexiFr2 = NoFr2 + sheari + MPIcount1;
		indexiFr3 = NoFr3 + sheari + MPIcount1;
		
		
		indexi1Er  = NoEr + 4 + sheari + MPIcount1;
		indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1;
		indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1;
		indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1;
		
		indexj1Er  = NoEr + 6 + shearj1 + MPIcount1;
		indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1;
		indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1;
		indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1;
		

		
		indexk1Er = NoEr + 8 + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1;
		
		
#endif

	return;

}

void is_js_ks_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = is;
	j = js;
	k = ks;

	if(ix3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/


	indexiEr  = NoEr + MPIcount1x + MPIcount1y;
	indexiFr1 = NoFr1 + MPIcount1x + MPIcount1y;
	indexiFr2 = NoFr2 + MPIcount1x + MPIcount1y;
	indexiFr3 = NoFr3 + MPIcount1x + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1x + MPIcount1y;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1x + MPIcount1y;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1y;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1y;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1y;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1y;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + sheari + MPIcount1y;
	
	
	indexi1Er  = NoEr + 4 + sheari + MPIcount1y;
	indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1y;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y;
	
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y;
	
}	
#endif

	return;


}




void is_js_ks_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = is;
	j = js;
	k = ks;
	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1 + count1x;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	if(ix1 == 4){
		indexi0Er  = NoEr + 6 + MPIcount1;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + MPIcount1;
	indexiFr1 = NoFr1 + MPIcount1;
	indexiFr2 = NoFr2 + MPIcount1;
	indexiFr3 = NoFr3 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2F;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2F;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2F;
	

	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1;
	indexiFr1 = NoFr1 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 4 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1;
	
	indexj0Er  = NoEr  + 8 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 6 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 6 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 6 + MPIcount1 + shearj0;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
#endif
	

	return;
}


void is_js_ks_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = is;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	/***************************/
	/* for z boundary, MPI part */
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + MPIcount1x + MPIcount1z;
	indexiFr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexiFr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexiFr3 = NoFr3 + MPIcount1x + MPIcount1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1x + MPIcount1z;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexi1Fr3 = NoFr3 + 2 +MPIcount1x + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 4 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1z;
	
	indexj0Er  = NoEr  + 8 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 6 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 6 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 6 + MPIcount1z + shearj0;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
}	
#endif

	return;

}

void is_js_ks_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = is;
	j = js;
	k = ks;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 6 + MPIcount1z + MPIcount1y;
		indexi0Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
		indexi0Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
		indexi0Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
	
#endif
	

	return;

}

void is_js_ks_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = is;
	j = js;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	/***************************/
	/* for y boundary, MPI part */
	
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
}	
#endif

	return;


}

/*======================================================*/

/* Now begin ie, js, ks */


void ie_js_ks_phy_phy_phy()
{
	int i, j, k, count1y, count1x;
	i = ie;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + count1y + count1x;
		indexk0Fr1 = NoFr1 + 8 + count1y + count1x;
		indexk0Fr2 = NoFr2 + 8 + count1y + count1x;
		indexk0Fr3 = NoFr3 + 8 + count1y + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + count1x;
		indexj0Fr1 = NoFr1 + 6 + count1x;
		indexj0Fr2 = NoFr2 + 6 + count1x;
		indexj0Fr3 = NoFr3 + 6 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	
	indexi0Er  = NoEr + count1x;
	indexi0Fr1 = NoFr1 + count1x;
	indexi0Fr2 = NoFr2 + count1x;
	indexi0Fr3 = NoFr3 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 2 + count1x;
	indexiFr1 = NoFr1 + 2 + count1x;
	indexiFr2 = NoFr2 + 2 + count1x;
	indexiFr3 = NoFr3 + 2 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary */
	if(ox1 == 4){
		indexi1Er  = NoEr;
		indexi1Fr1 = NoFr1;
		indexi1Fr2 = NoFr2;
		indexi1Fr3 = NoFr3;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	indexj1Er  = NoEr + 6 + count1x;
	indexj1Fr1 = NoFr1 + 4 + count1x;
	indexj1Fr2 = NoFr2 + 4 + count1x;
	indexj1Fr3 = NoFr3 + 4 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 6 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 6 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 6 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny ) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari;
	indexi0Fr1 = NoFr1 + sheari;
	indexi0Fr2 = NoFr2 + sheari;
	indexi0Fr3 = NoFr3 + sheari;
	
	indexiEr  = NoEr + 2 + sheari;
	indexiFr1 = NoFr1 + 2 + sheari;
	indexiFr2 = NoFr2 + 2 + sheari;
	indexiFr3 = NoFr3 + 2 + sheari;
	
	
	indexj1Er  = NoEr + 6 + shearj1;
	indexj1Fr1 = NoFr1 + 4 + shearj1;
	indexj1Fr2 = NoFr2 + 4+ shearj1;
	indexj1Fr3 = NoFr3 + 4 + shearj1;
	
	indexj0Er  = NoEr  + 8 + shearj0;
	indexj0Fr1 = NoFr1 + 6 + shearj0;
	indexj0Fr2 = NoFr2 + 6 + shearj0;
	indexj0Fr3 = NoFr3 + 6 + shearj0;
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	

#endif
	

	return;
}

void ie_js_ks_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = ie;
	j = js;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;
	if(ix2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1y;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	
	indexi0Er  = NoEr + MPIcount1;
	indexi0Fr1 = NoFr1 + MPIcount1;
	indexi0Fr2 = NoFr2 + MPIcount1;
	indexi0Fr3 = NoFr3 + MPIcount1;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids  + shiftx;
	/******************************/

	indexj1Er  = NoEr + 6 + MPIcount1;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny ) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 12 + sheark0;
			indexk0Fr1 = NoFr1 + 10 + sheark0;
			indexk0Fr2 = NoFr2 + 10 + sheark0;
			indexk0Fr3 = NoFr3 + 10 + sheark0;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli0)	sheari0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
		indexi0Er  = NoEr + sheari;
		indexi0Fr1 = NoFr1 + sheari;
		indexi0Fr2 = NoFr2 + sheari;
		indexi0Fr3 = NoFr3 + sheari;
		
		indexiEr  = NoEr + 2 + sheari;
		indexiFr1 = NoFr1 + 2 + sheari;
		indexiFr2 = NoFr2 + 2 + sheari;
		indexiFr3 = NoFr3 + 2 + sheari;
		
	
		indexj1Er  = NoEr + 6 + shearj1;
		indexj1Fr1 = NoFr1 + 4 + shearj1;
		indexj1Fr2 = NoFr2 + 4+ shearj1;
		indexj1Fr3 = NoFr3 + 4 + shearj1;
		
		indexj0Er  = NoEr  + 8 + shearj0;
		indexj0Fr1 = NoFr1 + 6 + shearj0;
		indexj0Fr2 = NoFr2 + 6 + shearj0;
		indexj0Fr3 = NoFr3 + 6 + shearj0;
		
		indexk1Er = NoEr + 10 + sheark1;
		indexk1Fr1 = NoFr1 + 8 + sheark1;
		indexk1Fr2 = NoFr2 + 8 + sheark1;
		indexk1Fr3 = NoFr3 + 8 + sheark1;
			
		
	}	
	
#endif

	return;

}

void ie_js_ks_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = ie;
	j = js;
	k = ks;

	
	if(ox1 == 4) count1x = 2;
	else	count1x = 0;
	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	
	/******************************/
	
	
	
	indexi0Er  = NoEr + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 2 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	if(ox1 == 4){
		indexi1Er  = NoEr + MPIcount1;
		indexi1Fr1 = NoFr1 + MPIcount1;
		indexi1Fr2 = NoFr2 + MPIcount1;
		indexi1Fr3 = NoFr3 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	
	

	indexj1Er  = NoEr + 6 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexj0Er  = NoEr  + MPIcount2 + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	

	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1;
	
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1;
	
	
#endif
	

	return;

}

void ie_js_ks_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = ie;
	j = js;
	k = ks;

	if(ix3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	
	indexi0Er  = NoEr + MPIcount1y + MPIcount1x;
	indexi0Fr1 = NoFr1 + MPIcount1y + MPIcount1x;
	indexi0Fr2 = NoFr2 + MPIcount1y + MPIcount1x;
	indexi0Fr3 = NoFr3 + MPIcount1y + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 6 + MPIcount1x + MPIcount1y;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1y;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1y;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1y;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1y;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari + MPIcount1y;
	indexi0Fr1 = NoFr1 + sheari + MPIcount1y;
	indexi0Fr2 = NoFr2 + sheari + MPIcount1y;
	indexi0Fr3 = NoFr3 + sheari + MPIcount1y;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	
	
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y;
	indexj1Fr2 = NoFr2 + 4 + shearj1 + MPIcount1y;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y;
	
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y;
	
}	
#endif

	return;


}



void ie_js_ks_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = ie;
	j = js;
	k = ks;
	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1 + count1x;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	
	indexi0Er  = NoEr + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 2 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){
		indexi1Er  = NoEr + MPIcount1;
		indexi1Fr1 = NoFr1 + MPIcount1;
		indexi1Fr2 = NoFr2 + MPIcount1;
		indexi1Fr3 = NoFr3 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	indexj1Er  = NoEr + 6 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	

	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2F;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2F;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2F;
	

	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	
	
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1;
	
	
	indexj0Er  = NoEr  + 8 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 6 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 6 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 6 + MPIcount1 + shearj0;
	
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
#endif

	return;
}



void ie_js_ks_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = ie;
	j = js;
	k = ks;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	/***************************/
	/* for z boundary, MPI part */
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z;
		indexj0Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z;
		indexj0Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	
	
	indexi0Er  = NoEr + MPIcount1x + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	
	
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4+ shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1z;
	
	
	indexj0Er  = NoEr  + 8 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 6 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 6 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 6 + MPIcount1z + shearj0;
	
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
}	
#endif
	

	return;

}

void ie_js_ks_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = ie;
	j = js;
	k = ks;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	/***************************/
	/* for y boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	
	
	
	indexi0Er  = NoEr + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr1 = NoFr1 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr2 = NoFr2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr3 = NoFr3 + MPIcount1z + MPIcount1y + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y + count1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + MPIcount1z + MPIcount1y;
		indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
		indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
		indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	indexj1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + sheari + MPIcount1z + MPIcount1y;
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
	
#endif
	

	return;

}

void ie_js_ks_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = ie;
	j = js;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	
	
	indexi0Er  = NoEr + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount1x + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount1z + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + sheari + MPIcount1z + MPIcount1y;
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1y + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 4 + shearj1 + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;

}		
#endif
	

	return;


}

/******************************************/
/* Now begin is, je, ks */


void is_je_ks_phy_phy_phy()
{
	int i, j, k, count1y, count1x;
	i = is;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + count1y + count1x;
		indexk0Fr1 = NoFr1 + 8 + count1y + count1x;
		indexk0Fr2 = NoFr2 + 8 + count1y + count1x;
		indexk0Fr3 = NoFr3 + 8 + count1y + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	indexj0Er  = NoEr + count1y;
	indexj0Fr1 = NoFr1 + count1y;
	indexj0Fr2 = NoFr2 + count1y;
	indexj0Fr3 = NoFr3 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for x boundary */
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + count1y;
		indexi0Fr1 = NoFr1 + 6 + count1y;
		indexi0Fr2 = NoFr2 + 6 + count1y;
		indexi0Fr3 = NoFr3 + 6 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2 + count1y;
	indexiFr1 = NoFr1 + 2 + count1y;
	indexiFr2 = NoFr2 + 2 + count1y;
	indexiFr3 = NoFr3 + 2 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + count1y;
	indexi1Fr1 = NoFr1 + 4 + count1y;
	indexi1Fr2 = NoFr2 + 4 + count1y;
	indexi1Fr3 = NoFr3 + 4 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr;
		indexj1Fr1 = NoFr1;
		indexj1Fr2 = NoFr2;
		indexj1Fr3 = NoFr3;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + count1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + count1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + count1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + count1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	indexj1Er  = NoEr + shearj1;
	indexj1Fr1 = NoFr1 + shearj1;
	indexj1Fr2 = NoFr2 + shearj1;
	indexj1Fr3 = NoFr3 + shearj1;
	
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari;
	indexiFr1 = NoFr1 + 4 + sheari;
	indexiFr2 = NoFr2 + 4 + sheari;
	indexiFr3 = NoFr3 + 4 + sheari;
	
	
	indexi1Er  = NoEr + 8 + sheari;
	indexi1Fr1 = NoFr1 + 6 + sheari;
	indexi1Fr2 = NoFr2 + 6 + sheari;
	indexi1Fr3 = NoFr3 + 6 + sheari;
	
	
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	
	
	
	
#endif

	return;
}


void is_je_ks_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = is;
	j = je;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;
	if(ox2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1y;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	
	indexj0Er  = NoEr + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;


	if(ox2 == 4){
		indexj1Er  = NoEr + MPIcount1;
		indexj1Fr1 = NoFr1 + MPIcount1;
		indexj1Fr2 = NoFr2 + MPIcount1;
		indexj1Fr3 = NoFr3 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 12 + sheark0;
			indexk0Fr1 = NoFr1 + 10 + sheark0;
			indexk0Fr2 = NoFr2 + 10 + sheark0;
			indexk0Fr3 = NoFr3 + 10 + sheark0;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari0 = 2;
	
	
		indexj1Er  = NoEr + shearj1;
		indexj1Fr1 = NoFr1 + shearj1;
		indexj1Fr2 = NoFr2 + shearj1;
		indexj1Fr3 = NoFr3 + shearj1;

	
		indexj0Er  = NoEr  + 2 + shearj0;
		indexj0Fr1 = NoFr1 + 2 + shearj0;
		indexj0Fr2 = NoFr2 + 2 + shearj0;
		indexj0Fr3 = NoFr3 + 2 + shearj0;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexiEr  = NoEr + 4 + sheari;
		indexiFr1 = NoFr1 + 4 + sheari;
		indexiFr2 = NoFr2 + 4 + sheari;
		indexiFr3 = NoFr3 + 4 + sheari;
		
		
		indexi1Er  = NoEr + 8 + sheari;
		indexi1Fr1 = NoFr1 + 6 + sheari;
		indexi1Fr2 = NoFr2 + 6 + sheari;
		indexi1Fr3 = NoFr3 + 6 + sheari;
		

		
		indexk1Er = NoEr + 10 + sheark1;
		indexk1Fr1 = NoFr1 + 8 + sheark1;
		indexk1Fr2 = NoFr2 + 8 + sheark1;
		indexk1Fr3 = NoFr3 + 8 + sheark1;
		
		
		
}	
#endif
	

	return;

}

void is_je_ks_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = is;
	j = je;
	k = ks;

	
	if(ix1 == 4) count1x = 2;
	else	count1x = 0;
	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	
	

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	

	
	indexj0Er  = NoEr + shearj0 + MPIcount1;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1;
	

	indexj1Er  = NoEr  + MPIcount2 + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1;
	
	
#endif

	return;

}

void is_je_ks_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = is;
	j = je;
	k = ks;

	if(ix3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	

	indexj0Er  = NoEr + MPIcount1y + MPIcount1x;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1x;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1x;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1x;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/


	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1x + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1y;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1y;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1y;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1y;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexj0Er  = NoEr + shearj0 + MPIcount1y;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1y;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1y;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1y;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y;
	
}	
#endif
	

	return;


}


void is_je_ks_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = is;
	j = je;
	k = ks;
	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	

	
	indexj0Er  = NoEr + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1 + count1y;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 2 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary */
	if(ox2 == 4){

		indexj1Er  = NoEr + MPIcount1;
		indexj1Fr1 = NoFr1 + MPIcount1;
		indexj1Fr2 = NoFr2 + MPIcount1;
		indexj1Fr3 = NoFr3 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	indexk1Er  = NoEr + 8 + MPIcount1 + count1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2F;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2F;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2F;
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount1;
	
	
	
	indexj0Er  = NoEr  + 2 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + shearj0;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
#endif
	

	return;
}



void is_je_ks_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = is;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	/***************************/
	/* for z boundary, MPI part */
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/



	
	indexj0Er  = NoEr + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1x + MPIcount1z + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1z + count1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr3 = NoFr3 + 4 +MPIcount1x + MPIcount1z + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + MPIcount1x + MPIcount1z;
		indexj1Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
		indexj1Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
		indexj1Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount1z;
	
	
	
	indexj0Er  = NoEr  + 2 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1z + shearj0;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
}
#endif
	

	return;

}


void is_je_ks_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = is;
	j = je;
	k = ks;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	
	
	indexj0Er  = NoEr + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount2y + MPIcount1z;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount2Fy + MPIcount1z;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount2Fy + MPIcount1z;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount2Fy + MPIcount1z;
	
	
	indexj0Er  = NoEr  + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1z + MPIcount1y;
	

	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;

	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
	
#endif
	

	return;

}

void is_je_ks_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = is;
	j = je;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	

	indexj0Er  = NoEr + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1z + MPIcount1x;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount2y + MPIcount1z;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount2Fy + MPIcount1z;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount2Fy + MPIcount1z;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount2Fy + MPIcount1z;
	
	
	indexj0Er  = NoEr  + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1z + MPIcount1y;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1z + MPIcount1y;
	
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
}	
#endif
	

	return;


}


/* begin ie, je, ks */

void ie_je_ks_phy_phy_phy()
{
	int i, j, k, count1y, count1x;
	i = ie;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;


		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + count1y + count1x;
		indexk0Fr1 = NoFr1 + 8 + count1y + count1x;
		indexk0Fr2 = NoFr2 + 8 + count1y + count1x;
		indexk0Fr3 = NoFr3 + 8 + count1y + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	indexj0Er  = NoEr + count1y;
	indexj0Fr1 = NoFr1 + count1y;
	indexj0Fr2 = NoFr2 + count1y;
	indexj0Fr3 = NoFr3 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexi0Er  = NoEr + 2 + count1x + count1y;
	indexi0Fr1 = NoFr1 + 2 + count1x + count1y;
	indexi0Fr2 = NoFr2 + 2 + count1x + count1y;
	indexi0Fr3 = NoFr3 + 2 + count1x + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + count1x + count1y;
	indexiFr1 = NoFr1 + 4 + count1x + count1y;
	indexiFr2 = NoFr2 + 4 + count1x + count1y;
	indexiFr3 = NoFr3 + 4 + count1x + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary */
	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + count1y;
		indexi1Fr1 = NoFr1 + 2 + count1y;
		indexi1Fr2 = NoFr2 + 2 + count1y;
		indexi1Fr3 = NoFr3 + 2 + count1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr;
		indexj1Fr1 = NoFr1;
		indexj1Fr2 = NoFr2;
		indexj1Fr3 = NoFr3;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ix2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 6 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 6 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 6 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny ) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 12 + sheark0;
		indexk0Fr1 = NoFr1 + 10 + sheark0;
		indexk0Fr2 = NoFr2 + 10 + sheark0;
		indexk0Fr3 = NoFr3 + 10 + sheark0;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexj1Er  = NoEr + shearj1;
	indexj1Fr1 = NoFr1 + shearj1;
	indexj1Fr2 = NoFr2 + shearj1;
	indexj1Fr3 = NoFr3 + shearj1;
	
	indexj0Er  = NoEr  + 2 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + shearj0;
	
	indexi0Er  = NoEr + 4 + sheari;
	indexi0Fr1 = NoFr1 + 4 + sheari;
	indexi0Fr2 = NoFr2 + 4 + sheari;
	indexi0Fr3 = NoFr3 + 4 + sheari;
	
	indexiEr  = NoEr + 6 + sheari;
	indexiFr1 = NoFr1 + 6 + sheari;
	indexiFr2 = NoFr2 + 6 + sheari;
	indexiFr3 = NoFr3 + 6 + sheari;
	
	
	
	
	indexk1Er = NoEr + 10 + sheark1;
	indexk1Fr1 = NoFr1 + 8 + sheark1;
	indexk1Fr2 = NoFr2 + 8 + sheark1;
	indexk1Fr3 = NoFr3 + 8 + sheark1;
	
#endif
	

	return;
}

void ie_je_ks_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = ie;
	j = je;
	k = ks;

	if(ix3 == 4) count1z = 2;
	else	count1z = 0;
	if(ox2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



		/* There are seven blocks */
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1y;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	
	indexj0Er  = NoEr + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids  + shiftx;
	/******************************/

	if(ox2 == 4){

		indexj1Er  = NoEr + MPIcount1;
		indexj1Fr1 = NoFr1 + MPIcount1;
		indexj1Fr2 = NoFr2 + MPIcount1;
		indexj1Fr3 = NoFr3 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny ) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ix3 == 4){
			if(col_shear < colk0)	sheark0 = 2;
			
			indexk0Er  = NoEr  + 12 + sheark0;
			indexk0Fr1 = NoFr1 + 10 + sheark0;
			indexk0Fr2 = NoFr2 + 10 + sheark0;
			indexk0Fr3 = NoFr3 + 10 + sheark0;
		}
		else {
			sheark0 = 2;
		}
		
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli0)	sheari0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < colk1)	sheark1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
		indexj1Er  = NoEr + shearj1;
		indexj1Fr1 = NoFr1 + shearj1;
		indexj1Fr2 = NoFr2 + shearj1;
		indexj1Fr3 = NoFr3 + shearj1;
	
		indexj0Er  = NoEr  + 2 + shearj0;
		indexj0Fr1 = NoFr1 + 2 + shearj0;
		indexj0Fr2 = NoFr2 + 2 + shearj0;
		indexj0Fr3 = NoFr3 + 2 + shearj0;
		
		indexi0Er  = NoEr + 4 + sheari;
		indexi0Fr1 = NoFr1 + 4 + sheari;
		indexi0Fr2 = NoFr2 + 4 + sheari;
		indexi0Fr3 = NoFr3 + 4 + sheari;
		
		indexiEr  = NoEr + 6 + sheari;
		indexiFr1 = NoFr1 + 6 + sheari;
		indexiFr2 = NoFr2 + 6 + sheari;
		indexiFr3 = NoFr3 + 6 + sheari;
		
		

		
		indexk1Er = NoEr + 10 + sheark1;
		indexk1Fr1 = NoFr1 + 8 + sheark1;
		indexk1Fr2 = NoFr2 + 8 + sheark1;
		indexk1Fr3 = NoFr3 + 8 + sheark1;
		
		
	}	
	
#endif

	return;

}

void ie_je_ks_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = ie;
	j = je;
	k = ks;

	
	if(ox1 == 4) count1x = 2;
	else	count1x = 0;
	if(ix3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	


	
	
	indexj0Er  = NoEr + MPIcount1;
	indexj0Fr1 = NoFr1 + MPIcount1;
	indexj0Fr2 = NoFr2 + MPIcount1;
	indexj0Fr3 = NoFr3 + MPIcount1;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;	
	
	
	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + MPIcount1;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	
	
	/***************************/
	/* for y boundary, MPI part */
	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexj1Er  = NoEr  + MPIcount2 + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;
	
		
	indexj0Er  = NoEr + shearj0 + MPIcount1;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1;
	
	
#endif

	return;

}

void ie_je_ks_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = ie;
	j = je;
	k = ks;

	if(ix3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */


	
	if(ix3 == 4){
		indexk0Er  = NoEr + 10 + MPIcount1y + MPIcount1x;
		indexk0Fr1 = NoFr1 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr2 = NoFr2 + 8 + MPIcount1y + MPIcount1x;
		indexk0Fr3 = NoFr3 + 8 + MPIcount1y + MPIcount1x;
		colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix3 == 1 || ix3 == 5){
					
		theta[6] += theta[0];
		theta[9] -= theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] -= varphi[1];
	}
	else if(ix3 == 2){
		theta[6] += theta[0];
		theta[9] += theta[1];
			
		phi[6] += phi[0];
		phi[7] += phi[1];
	
		psi[6] += psi[0];
		psi[7] += psi[1];
			
		varphi[6] += varphi[0];
		varphi[7] += varphi[1];
	}
	else if(ix3 == 3){
			/* Do nothing */
	}	
	

	

	indexj0Er  = NoEr + MPIcount1x + MPIcount1y;
	indexj0Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
	indexj0Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
	indexj0Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexi0Er  = NoEr + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
	
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ix3 == 4){
		if(col_shear < colk0)	sheark0 = 2;
		
		indexk0Er  = NoEr  + 10 + sheark0 + MPIcount1y;
		indexk0Fr1 = NoFr1 + 8 + sheark0 + MPIcount1y;
		indexk0Fr2 = NoFr2 + 8 + sheark0 + MPIcount1y;
		indexk0Fr3 = NoFr3 + 8 + sheark0 + MPIcount1y;
	}
	else {
		sheark0 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
	
	indexj0Er  = NoEr + shearj0 + MPIcount1y;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1y;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1y;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1y;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1y;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y;
	
}	
#endif

	return;


}




void ie_je_ks_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = ie;
	j = je;
	k = ks;
	
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(lx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount2;
	indexk0Fr1 = NoFr1 + MPIcount2F;
	indexk0Fr2 = NoFr2 + MPIcount2F;
	indexk0Fr3 = NoFr3 + MPIcount2F;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	


	
	indexj0Er  = NoEr + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1x + count1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x + count1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x + count1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1x + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + MPIcount1 + count1y;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr + MPIcount1;
		indexj1Fr1 = NoFr1 + MPIcount1;
		indexj1Fr2 = NoFr2 + MPIcount1;
		indexj1Fr3 = NoFr3 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + MPIcount1 + count1x + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1 + count1x + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1 + count1x + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1 + count1x + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2F;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2F;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2F;
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount1;
	
	
	indexj0Er  = NoEr  + 2 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + shearj0;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1;
	
	
#endif

	return;
}

void ie_je_ks_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = ie;
	j = je;
	k = ks;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	/***************************/
	/* for z boundary, MPI part */
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	
	indexj0Er  = NoEr + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr1 = NoFr1 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr2 = NoFr2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr3 = NoFr3 + MPIcount1x + MPIcount1z + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	indexi0Er  = NoEr + 2 +  MPIcount1x + MPIcount1z + count1y;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	
	if(ox2 == 4){

		indexj1Er  = NoEr + MPIcount1x + MPIcount1z;
		indexj1Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
		indexj1Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
		indexj1Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	indexk1Er  = NoEr + 8 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z + count1y;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	
#ifdef SHEARING_BOX
	
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	indexj1Er  = NoEr + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + shearj1 + MPIcount1z;
	
	
	indexj0Er  = NoEr  + 2 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1z + shearj0;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	
	indexk1Er = NoEr + 10 + sheark1 + MPIcount1z;
	indexk1Fr1 = NoFr1 + 8 + sheark1 + MPIcount1z;
	indexk1Fr2 = NoFr2 + 8 + sheark1 + MPIcount1z;
	indexk1Fr3 = NoFr3 + 8 + sheark1 + MPIcount1z;
	
}	
#endif	
	
	return;

}

void ie_je_ks_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = ie;
	j = je;
	k = ks;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	
	indexj0Er  = NoEr + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

	
	
	
	indexi0Er  = NoEr + 2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for y boundary, MPI part */
	
	indexj1Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids  + shifty;
	
	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
		
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
	indexj0Er  = NoEr + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	indexi0Er  = NoEr + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
	
#endif
	

	return;

}

void ie_je_ks_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = ie;
	j = je;
	k = ks;

		
	shiftz = 4 * Ny * Nx * Nz * (lx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	

	if(lx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	/***************************/
	/* for z boundary, MPI part */
	
	indexk0Er  = NoEr + MPIcount2z;
	indexk0Fr1 = NoFr1 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + MPIcount2Fz;
	colk0 = 4 * (ke - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/


	

	indexj0Er  = NoEr + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr1 = NoFr1 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr2 = NoFr2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr3 = NoFr3 + MPIcount1y + MPIcount1z + MPIcount1x;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexi0Er  = NoEr + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount1z + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/

	indexk1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colk1 = 4 * (k - ks + 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk0Er  = NoEr  + sheark0 + MPIcount2z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount2Fz;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount2Fz;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount2Fz;
	
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
	indexj0Er  = NoEr + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	
	indexk1Er = NoEr + 8 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr1 = NoFr1 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr2 = NoFr2 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	indexk1Fr3 = NoFr3 + 6 + sheark1 + MPIcount1y + MPIcount1z;
	
}	
#endif

	return;


}


/* begin is, js, ke */


void is_js_ke_phy_phy_phy()
{
	int i, j, k, count1y, count1x, count1z;
	i = is;
	j = js;
	k = ke;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;


		/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + count1x + count1z;
		indexj0Fr1 = NoFr1 + 8 + count1x + count1z;
		indexj0Fr2 = NoFr2 + 8 + count1x + count1z;
		indexj0Fr3 = NoFr3 + 8 + count1x + count1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for x boundary */
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + count1z;
		indexi0Fr1 = NoFr1 + 6 + count1z;
		indexi0Fr2 = NoFr2 + 6 + count1z;
		indexi0Fr3 = NoFr3 + 6 + count1z;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2 + count1z;
	indexiFr1 = NoFr1 + 2 + count1z;
	indexiFr2 = NoFr2 + 2 + count1z;
	indexiFr3 = NoFr3 + 2 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + count1z;
	indexi1Fr1 = NoFr1 + 4 + count1z;
	indexi1Fr2 = NoFr2 + 4 + count1z;
	indexi1Fr3 = NoFr3 + 4 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + count1z + count1x;
	indexj1Fr1 = NoFr1 + 6 + count1z + count1x;
	indexj1Fr2 = NoFr2 + 6 + count1z + count1x;
	indexj1Fr3 = NoFr3 + 6 + count1z + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){

		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	indexiEr  = NoEr + 2 + sheari + count1z;
	indexiFr1 = NoFr1 + 2 + sheari + count1z;
	indexiFr2 = NoFr2 + 2 + sheari + count1z;
	indexiFr3 = NoFr3 + 2 + sheari + count1z;
	
	
	indexi1Er  = NoEr + 6 + sheari + count1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + count1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + count1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + count1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + count1z;
	
	indexj0Er  = NoEr  + 10 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 8 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 8 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 8 + shearj0 + count1z;
	
#endif

	return;
}

void is_js_ke_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = is;
	j = js;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;
	if(ix2 == 4) count1y = 2;
	else	count1y = 0;
	
	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1 + count1z;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
		
		indexiEr  = NoEr + 2 + sheari + count1z;
		indexiFr1 = NoFr1 + 2 + sheari + count1z;
		indexiFr2 = NoFr2 + 2 + sheari + count1z;
		indexiFr3 = NoFr3 + 2 + sheari + count1z;
		
		
		indexi1Er  = NoEr + 6 + sheari + count1z;
		indexi1Fr1 = NoFr1 + 4 + sheari + count1z;
		indexi1Fr2 = NoFr2 + 4 + sheari + count1z;
		indexi1Fr3 = NoFr3 + 4 + sheari + count1z;
		
		indexj1Er  = NoEr + 8 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 6 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 6 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 6 + shearj1 + count1z;
		
		indexj0Er  = NoEr  + 10 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 8 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 8 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 8 + shearj0 + count1z;
		
	}	
#endif
	

	return;

}

void is_js_ke_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = is;
	j = js;
	k = ke;

	
	if(ix1 == 4) count1x = 2;
	else	count1x = 0;
	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	
	/******************************/
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1 + count1z;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	
	

	indexiEr  = NoEr + 2 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 2 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1 + count1z + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1z + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1z + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1z + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1 + count1z;
	
	
	indexj0Er  = NoEr  + MPIcount2 + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1 + count1z;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1 + count1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1 + count1z;
	

	
	
#endif

	return;

}

void is_js_ke_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = is;
	j = js;
	k = ke;

	if(ox3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



	/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/


	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	

	indexj1Er  = NoEr + 8 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

		
	
	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1x + MPIcount1y;
		indexk1Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
		indexk1Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
		indexk1Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}


	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1y;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1y;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1y;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1y;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + count1z;
	
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1y + count1z;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y + count1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y + count1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y + count1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y + count1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + count1z;
	
}
	
#endif	
	
	return;


}




void is_js_ke_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = is;
	j = js;
	k = ke;
	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 2 + MPIcount1;
	indexiFr1 = NoFr1 + 2 + MPIcount1;
	indexiFr2 = NoFr2 + 2 + MPIcount1;
	indexiFr3 = NoFr3 + 2 + MPIcount1;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	/* MPI boundary */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2F;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2F;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
	
	indexj0Er  = NoEr  + 10 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 8 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 8 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 8 + MPIcount1 + shearj0;
	

	
	
#endif
	

	return;
}

void is_js_ke_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = is;
	j = js;
	k = ke;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1x + MPIcount1z;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1x + MPIcount1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1x + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1x + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 +MPIcount1x + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1x + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for z boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1z;
	
	indexj0Er  = NoEr  + 10 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 8 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 8 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 8 + MPIcount1z + shearj0;
	
	
}	
	
#endif
	


	return;

}


void is_js_ke_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = is;
	j = js;
	k = ke;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1y;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
		indexi0Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
		indexi0Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
		indexi0Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	
	
	
#endif
	
	return;

}

void is_js_ke_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = is;
	j = js;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	/***************************/
	/* for y boundary, MPI part */
	
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z; 
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/

	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 2 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 2 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	
	
}	
#endif	
	
	return;


}

/* begin ie, js, ke */


void ie_js_ke_phy_phy_phy()
{
	int i, j, k, count1y, count1x, count1z;
	i = ie;
	j = js;
	k = ke;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;



	/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + count1x + count1z;
		indexj0Fr1 = NoFr1 + 8 + count1x + count1z;
		indexj0Fr2 = NoFr2 + 8 + count1x + count1z;
		indexj0Fr3 = NoFr3 + 8 + count1x + count1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	
	indexi0Er  = NoEr + 2 + count1x + count1z;
	indexi0Fr1 = NoFr1 + 2 + count1x + count1z;
	indexi0Fr2 = NoFr2 + 2 + count1x + count1z;
	indexi0Fr3 = NoFr3 + 2 + count1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + count1x + count1z;
	indexiFr1 = NoFr1 + 4 + count1x + count1z;
	indexiFr2 = NoFr2 + 4 + count1x + count1z;
	indexiFr3 = NoFr3 + 4 + count1x + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary */
	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + count1z;
		indexi1Fr1 = NoFr1 + 2 + count1z;
		indexi1Fr2 = NoFr2 + 2 + count1z;
		indexi1Fr3 = NoFr3 + 2 + count1z;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	
	indexj1Er  = NoEr + 8 + count1x + count1z;
	indexj1Fr1 = NoFr1 + 6 + count1x + count1z;
	indexj1Fr2 = NoFr2 + 6 + count1x + count1z;
	indexj1Fr3 = NoFr3 + 6 + count1x + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny ) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	indexi0Er  = NoEr + 2 + sheari + count1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + count1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + count1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + count1z;
	
	indexiEr  = NoEr + 4 + sheari + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + count1z;
	
	
	indexj1Er  = NoEr + 8 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + count1z;
	
	indexj0Er  = NoEr  + 10 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 8 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 8 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 8 + shearj0 + count1z;	
	
#endif
	

	return;
}


void ie_js_ke_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = ie;
	j = js;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;
	if(ix2 == 4) count1y = 2;
	else	count1y = 0;
	
	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z ;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;



	
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1 + count1z;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	

	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids  + shiftx;
	/******************************/

	indexj1Er  = NoEr + 8 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){

		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny ) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
		
	
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli0)	sheari0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
		
		indexi0Er  = NoEr + 2 + sheari + count1z;
		indexi0Fr1 = NoFr1 + 2 + sheari + count1z;
		indexi0Fr2 = NoFr2 + 2 + sheari + count1z;
		indexi0Fr3 = NoFr3 + 2 + sheari + count1z;
		
		indexiEr  = NoEr + 4 + sheari + count1z;
		indexiFr1 = NoFr1 + 4 + sheari + count1z;
		indexiFr2 = NoFr2 + 4 + sheari + count1z;
		indexiFr3 = NoFr3 + 4 + sheari + count1z;
		
		
		indexj1Er  = NoEr + 8 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 6 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 6 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 6 + shearj1 + count1z;
		
		indexj0Er  = NoEr  + 10 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 8 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 8 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 8 + shearj0 + count1z;

		
		
	}	
	
#endif

	return;

}

void ie_js_ke_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = ie;
	j = js;
	k = ke;

	
	if(ox1 == 4) count1x = 2;
	else	count1x = 0;
	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(lx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2;
	indexj0Fr1 = NoFr1 + MPIcount2F;
	indexj0Fr2 = NoFr2 + MPIcount2F;
	indexj0Fr3 = NoFr3 + MPIcount2F;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	
	/******************************/
	
	
	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1x + count1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x + count1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x + count1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1x + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + MPIcount1 + count1z;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	
	

	indexj1Er  = NoEr + 8 + MPIcount1 + count1x + count1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x + count1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x + count1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	

	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(lx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;	
	if(col_shear < coli)	sheari1 = 2;
	
	indexj0Er  = NoEr  + MPIcount2 + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2F + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2F + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2F + shearj0;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1 + count1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1 + count1z;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1 + count1z;
	
	
	
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1 + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1 + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1 + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1 + count1z;
	

	
	
#endif
	

	return;

}

void ie_js_ke_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = ie;
	j = js;
	k = ke;

	if(ox3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



	/* There are seven blocks */


	indexk0Er  = NoEr + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	

	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y;
	indexj0Fr1 = NoFr1 + MPIcount2Fy;
	indexj0Fr2 = NoFr2 + MPIcount2Fy;
	indexj0Fr3 = NoFr3 + MPIcount2Fy;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	
	indexi0Er  = NoEr + 2 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 8 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y + count1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	if(ox3 == 4){

		indexk1Er  = NoEr + MPIcount1x + MPIcount1y;
		indexk1Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
		indexk1Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
		indexk1Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(lx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1y;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1y;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1y;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1y;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + count1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1y + count1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y + count1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y + count1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y + count1z;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y + count1z;
	
	
	
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + count1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + count1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + count1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + count1z;
	
	
}	
	
#endif
	

	return;


}




void ie_js_ke_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = ie;
	j = js;
	k = ke;
	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary */
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1 + count1x;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1 + count1x;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1 + count1x;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1 + count1x;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	

	
	indexi0Er  = NoEr + 2 + MPIcount1 + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1 + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1 + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1 + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){
		indexi1Er  = NoEr + 2 + MPIcount1;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	indexj1Er  = NoEr + 8 + MPIcount1 + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1 + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1 + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1 + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/************************/
	/* MPI z direction */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/************************/
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2F;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2F;
	
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1;
	
	
	
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1;
	
	
	indexj0Er  = NoEr  + 10 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 8 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 8 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 8 + MPIcount1 + shearj0;
	
	
	
	
#endif

	return;
}

void ie_js_ke_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = ie;
	j = js;
	k = ke;

	if(ix2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	if(ix2 == 4){
		indexj0Er  = NoEr + 10 + MPIcount1x + MPIcount1z;
		indexj0Fr1 = NoFr1 + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr2 = NoFr2 + 8 + MPIcount1x + MPIcount1z;
		indexj0Fr3 = NoFr3 + 8 + MPIcount1x + MPIcount1z;
		colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix2 == 1 || ix2 == 5){
					
		theta[6] += theta[2];
		theta[8] -= theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] -= psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 2){
		theta[6] += theta[2];
		theta[8] += theta[3];
			
		phi[6] += phi[2];
		phi[7] += phi[3];
	
		psi[6] += psi[2];
		psi[7] += psi[3];
			
		varphi[6] += varphi[2];
		varphi[7] += varphi[3];
	}
	else if(ix2 == 3){
			/* Do nothing */
	}	
	
	
	indexi0Er  = NoEr + 2 + MPIcount1x + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 8 + MPIcount1x + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for z boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1z;
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z;
	
	
	
	
	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1z;
	
	
	indexj0Er  = NoEr  + 10 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 8 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 8 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 8 + MPIcount1z + shearj0;
	
}	
	
#endif

	return;

}

void ie_js_ke_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = ie;
	j = js;
	k = ke;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(lx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1z + MPIcount1y;
	indexk0Fr1 = NoFr1 + MPIcount1z + MPIcount1y;
	indexk0Fr2 = NoFr2 + MPIcount1z + MPIcount1y;
	indexk0Fr3 = NoFr3 + MPIcount1z + MPIcount1y;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary, MPI part */
	
	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	
	
	
	indexi0Er  = NoEr + 2 +  MPIcount1z + MPIcount1y + count1x;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr1 = NoFr1 + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr2 = NoFr2 + 2 + MPIcount1z + MPIcount1y;
		indexi1Fr3 = NoFr3 + 2 + MPIcount1z + MPIcount1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z + MPIcount1y;
	

	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	

	
	
#endif
	

	return;

}

void ie_js_ke_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = ie;
	j = js;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (lx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(lx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	indexk0Er  = NoEr + MPIcount1x + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + MPIcount1x + MPIcount1y + MPIcount1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	/***************************/
	/* for y boundary, MPI part */

	indexj0Er  = NoEr + MPIcount2y + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (je - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/

	
	
	
	indexi0Er  = NoEr + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount1z + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	indexj1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexj1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js + 1) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for z boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(lx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexj0Er  = NoEr  + MPIcount2y + shearj0 + MPIcount1z;
	indexj0Fr1 = NoFr1 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr2 = NoFr2 + MPIcount2Fy + shearj0 + MPIcount1z;
	indexj0Fr3 = NoFr3 + MPIcount2Fy + shearj0 + MPIcount1z;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	indexi0Er  = NoEr + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 2 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 2 + sheari + MPIcount1y + MPIcount1z;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z + MPIcount1y;
	

	indexj1Er  = NoEr + 8 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr1 = NoFr1 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr2 = NoFr2 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	indexj1Fr3 = NoFr3 + 6 + shearj1 + MPIcount1y + MPIcount1z;
	
	
}	
	
#endif

	return;


}


/* begin is, je, ke */


void is_je_ke_phy_phy_phy()
{
	int i, j, k, count1y, count1z;
	i = is;
	j = je;
	k = ke;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;



	/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
		


	
	indexj0Er  = NoEr + 2 + count1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + count1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + count1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + count1y + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	/***************************/
	/* for x boundary */
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + count1y + count1z;
		indexi0Fr1 = NoFr1 + 8 + count1y + count1z;
		indexi0Fr2 = NoFr2 + 8 + count1y + count1z;
		indexi0Fr3 = NoFr3 + 8 + count1y + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 4 + count1y + count1z;
	indexiFr1 = NoFr1 + 4 + count1y + count1z;
	indexiFr2 = NoFr2 + 4 + count1y + count1z;
	indexiFr3 = NoFr3 + 4 + count1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + count1y + count1z;
	indexi1Fr1 = NoFr1 + 6 + count1y + count1z;
	indexi1Fr2 = NoFr2 + 6 + count1y + count1z;
	indexi1Fr3 = NoFr3 + 6 + count1y + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + count1z;
		indexj1Fr1 = NoFr1 + 2 + count1z;
		indexj1Fr2 = NoFr2 + 2 + count1z;
		indexj1Fr3 = NoFr3 + 2 + count1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}
	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	
	indexj1Er  = NoEr + 2 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + count1z;
	
	indexj0Er  = NoEr  + 4 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 4 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 4 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 4 + shearj0 + count1z;
	
	indexiEr  = NoEr + 6 + sheari + count1z;
	indexiFr1 = NoFr1 + 6 + sheari + count1z;
	indexiFr2 = NoFr2 + 6 + sheari + count1z;
	indexiFr3 = NoFr3 + 6 + sheari + count1z;
	
	
	indexi1Er  = NoEr + 10 + sheari + count1z;
	indexi1Fr1 = NoFr1 + 8 + sheari + count1z;
	indexi1Fr2 = NoFr2 + 8 + sheari + count1z;
	indexi1Fr3 = NoFr3 + 8 + sheari + count1z;
	
	
#endif

	return;
}

void is_je_ke_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = is;
	j = je;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;
	if(ox2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);

	if(lx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
		


	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2;
	indexi0Fr1 = NoFr1 + MPIcount2F;
	indexi0Fr2 = NoFr2 + MPIcount2F;
	indexi0Fr3 = NoFr3 + MPIcount2F;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1 + count1y + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1y + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1y + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1y + count1z;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1y + count1z;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1y + count1z;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1y + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;


	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + MPIcount1 + count1z;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}
	
	
	
	
#ifdef SHEARING_BOX
if(lx1 > ID)
{
		
		jshearing = Ny * jproc + (j - js) - joffset;
		
		if(jshearing < 0) {
			jshearing += NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli0 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
		if(col_shear < colj1)	shearj1 = 2;
		if(col_shear < coli)	sheari0 = 2;
		
		
		
		indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	
		indexj1Er  = NoEr + 2 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 2 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 2 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 2 + shearj1 + count1z;
	
		indexj0Er  = NoEr  + 4 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 4 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 4 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 4 + shearj0 + count1z;
		
		indexiEr  = NoEr + 6 + sheari + count1z;
		indexiFr1 = NoFr1 + 6 + sheari + count1z;
		indexiFr2 = NoFr2 + 6 + sheari + count1z;
		indexiFr3 = NoFr3 + 6 + sheari + count1z;
		
		
		indexi1Er  = NoEr + 10 + sheari + count1z;
		indexi1Fr1 = NoFr1 + 8 + sheari + count1z;
		indexi1Fr2 = NoFr2 + 8 + sheari + count1z;
		indexi1Fr3 = NoFr3 + 8 + sheari + count1z;
		
	}	
#endif

	return;

}

void is_je_ke_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = is;
	j = je;
	k = ke;

	
	if(ix1 == 4) count1x = 2;
	else	count1x = 0;
	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}

	/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + MPIcount1 + count1z;
		indexi0Fr1 = NoFr1 + 8 + MPIcount1 + count1z;
		indexi0Fr2 = NoFr2 + 8 + MPIcount1 + count1z;
		indexi0Fr3 = NoFr3 + 8 + MPIcount1 + count1z;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	
	

	indexiEr  = NoEr + 4 + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1 + count1z;
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1 + count1z;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1 + count1z;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1 + count1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1 + count1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1 + count1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1 + count1z;
	

	indexj1Er  = NoEr  + MPIcount2 + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;
	
	
	
	
#endif

	return;

}

void is_je_ke_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = is;
	j = je;
	k = ke;

	if(ox3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */


	indexk0Er  = NoEr + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1x + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1x + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1x + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1x + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/


	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y + count1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/


	
	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1x + MPIcount1y;
		indexk1Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
		indexk1Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
		indexk1Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1y;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1y;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1y;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1y;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + count1z;
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + count1z;
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1y + count1z;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1y + count1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1y + count1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1y + count1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1y + count1z;
	
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
	
}	
	
#endif
	

	return;


}



void is_je_ke_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = is;
	j = je;
	k = ke;
	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + MPIcount1 + count1y;
		indexi0Fr1 = NoFr1 + 8 + MPIcount1 + count1y;
		indexi0Fr2 = NoFr2 + 8 + MPIcount1 + count1y;
		indexi0Fr3 = NoFr3 + 8 + MPIcount1 + count1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	

	indexiEr  = NoEr + 4 + MPIcount1 + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1 + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1 + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1 + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1 + count1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1 + count1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1 + count1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1 + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary */
	if(ox2 == 4){

		indexj1Er  = NoEr + 2 + MPIcount1;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	

	/********************************/
	/* MPI boundary in z direction */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/************************/
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2F;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2F;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
	
	indexj1Er  = NoEr + 2 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + MPIcount1;
	
	indexj0Er  = NoEr  + 4 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 4 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 4 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 4 + MPIcount1 + shearj0;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1;
	
	
	indexi1Er  = NoEr + 10 + sheari + MPIcount1;
	indexi1Fr1 = NoFr1 + 8 + sheari + MPIcount1;
	indexi1Fr2 = NoFr2 + 8 + sheari + MPIcount1;
	indexi1Fr3 = NoFr3 + 8 + sheari + MPIcount1;
	
	
	
	
#endif
	

	return;
}

void is_je_ke_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = is;
	j = je;
	k = ke;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(lx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1x + MPIcount1z;
	indexk0Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexk0Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexk0Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	/***************************/
	/* for x boundary, MPI part */
	indexi0Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexi1Fr3 = NoFr3 + 6 +MPIcount1x + MPIcount1z + count1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1z;
	
	indexj1Er  = NoEr + 2 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + MPIcount1z;
	
	indexj0Er  = NoEr  + 4 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 4 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 4 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 4 + MPIcount1z + shearj0;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	
	indexi1Er  = NoEr + 10 + sheari + MPIcount1z;
	indexi1Fr1 = NoFr1 + 8 + sheari + MPIcount1z;
	indexi1Fr2 = NoFr2 + 8 + sheari + MPIcount1z;
	indexi1Fr3 = NoFr3 + 8 + sheari + MPIcount1z;
	
	
}	
	
#endif



	return;

}

void is_je_ke_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = is;
	j = je;
	k = ke;

	if(ix1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */

	
	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	
	
	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

	
	
	if(ix1 == 4){
		indexi0Er  = NoEr + 10 + MPIcount1z + MPIcount1y;
		indexi0Fr1 = NoFr1 + 8 + MPIcount1z + MPIcount1y;
		indexi0Fr2 = NoFr2 + 8 + MPIcount1z + MPIcount1y;
		indexi0Fr3 = NoFr3 + 8 + MPIcount1z + MPIcount1y;
		coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids;
	}/* periodic for z direction */
	else if(ix1 == 1 || ix1 == 5){
					
		theta[6] += theta[4];
		theta[7] -= theta[5];
			
		phi[6] += phi[4];
		phi[7] -= phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 2){
		theta[6] += theta[4];
		theta[7] += theta[5];
			
		phi[6] += phi[4];
		phi[7] += phi[5];
	
		psi[6] += psi[4];
		psi[7] += psi[5];
			
		varphi[6] += varphi[4];
		varphi[7] += varphi[5];
	}
	else if(ix1 == 3){
			/* Do nothing */
	}	


	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1z + MPIcount1y;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/
	

	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
	
	
#endif


	return;

}

void is_je_ke_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = is;
	j = je;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (lx1 - ID);
	

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(lx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	/***************************/
	/* for x boundary, MPI part */
	
	indexi0Er  = NoEr + MPIcount2x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (ie - is) + count_Grids + shiftx;
	/******************************/

	indexiEr  = NoEr + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	indexi1Er  = NoEr + 8 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexi1Fr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is + 1) + count_Grids;

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	/***************************/
	/* for z boundary, MPI part */

	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;

	/******************************/

	
	
	
#ifdef SHEARING_BOX
	
if(lx1 > ID)
{	
	jshearing = Ny * jproc + (j - js) - joffset;
	
	if(jshearing < 0) {
		jshearing += NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (ie - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli0 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari1 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari0 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	
	indexi0Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi0Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 4 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 4 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexi1Er  = NoEr + 8 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr1 = NoFr1 + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr2 = NoFr2 + 6 + sheari + MPIcount1y + MPIcount1z;
	indexi1Fr3 = NoFr3 + 6 + sheari + MPIcount1y + MPIcount1z;
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
	
}	
#endif

	return;


}


/* begin ie, je, ke */


void ie_je_ke_phy_phy_phy()
{
	int i, j, k, count1y, count1x, count1z;
	i = ie;
	j = je;
	k = ke;

	
	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;


	/* There are seven blocks */
	
	indexk0Er  = NoEr + count1z;
	indexk0Fr1 = NoFr1 + count1z;
	indexk0Fr2 = NoFr2 + count1z;
	indexk0Fr3 = NoFr3 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
		


	
	indexj0Er  = NoEr + 2 + count1z + count1y;
	indexj0Fr1 = NoFr1 + 2 + count1z + count1y;
	indexj0Fr2 = NoFr2 + 2 + count1z + count1y;
	indexj0Fr3 = NoFr3 + 2 + count1z + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexi0Er  = NoEr + 4 + count1z + count1y + count1x;
	indexi0Fr1 = NoFr1 + 4 + count1z + count1y + count1x;
	indexi0Fr2 = NoFr2 + 4 + count1z + count1y + count1x;
	indexi0Fr3 = NoFr3 + 4 + count1z + count1y + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 6 + count1z + count1y + count1x;
	indexiFr1 = NoFr1 + 6 + count1z + count1y + count1x;
	indexiFr2 = NoFr2 + 6 + count1z + count1y + count1x;
	indexiFr3 = NoFr3 + 6 + count1z + count1y + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary */
	if(ox1 == 4){
		indexi1Er  = NoEr + 4 + count1y + count1z;
		indexi1Fr1 = NoFr1 + 4 + count1y + count1z;
		indexi1Fr2 = NoFr2 + 4 + count1y + count1z;
		indexi1Fr3 = NoFr3 + 4 + count1y + count1z;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	


	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + count1z;
		indexj1Fr1 = NoFr1 + 2 + count1z;
		indexj1Fr2 = NoFr2 + 2 + count1z;
		indexj1Fr3 = NoFr3 + 2 + count1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ix2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	if(ox3 == 4){
		indexk1Er  = NoEr;
		indexk1Fr1 = NoFr1;
		indexk1Fr2 = NoFr2;
		indexk1Fr3 = NoFr3;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}
	
	
	
	
#ifdef SHEARING_BOX
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny ) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1;
		indexk1Fr1 = NoFr1 + sheark1;
		indexk1Fr2 = NoFr2 + sheark1;
		indexk1Fr3 = NoFr3 + sheark1;
	}
	else {
		sheark1 = 2;
	}
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;		
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
	indexk0Er = NoEr + sheark0 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	
	indexj1Er  = NoEr + 2 + shearj1 + count1z;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + count1z;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + count1z;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + count1z;
	
	indexj0Er  = NoEr  + 4 + shearj0 + count1z;
	indexj0Fr1 = NoFr1 + 4 + shearj0 + count1z;
	indexj0Fr2 = NoFr2 + 4 + shearj0 + count1z;
	indexj0Fr3 = NoFr3 + 4 + shearj0 + count1z;
	
	
	indexi0Er  = NoEr + 6 + sheari + count1z;
	indexi0Fr1 = NoFr1 + 6 + sheari + count1z;
	indexi0Fr2 = NoFr2 + 6 + sheari + count1z;
	indexi0Fr3 = NoFr3 + 6 + sheari + count1z;
	
	indexiEr  = NoEr + 8 + sheari + count1z;
	indexiFr1 = NoFr1 + 8 + sheari + count1z;
	indexiFr2 = NoFr2 + 8 + sheari + count1z;
	indexiFr3 = NoFr3 + 8 + sheari + count1z;
	
	
#endif
	

	return;
}

void ie_je_ke_MPI_phy_phy()
{

	int i, j, k, count1z, shiftx, count1y;
	i = ie;
	j = je;
	k = ke;

	if(ox3 == 4) count1z = 2;
	else	count1z = 0;
	if(ox2 == 4) count1y = 2;
	else	count1y = 0;
	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);

	if(rx1 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1z + count1y;
		MPIcount2F = 8 + count1z + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	

	
	indexi0Er  = NoEr + 4 + MPIcount1 + count1y + count1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1y + count1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1y + count1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1y + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 6 + MPIcount1 + count1y + count1z;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1y + count1z;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1y + count1z;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2;
	indexi1Fr1 = NoFr1 + MPIcount2F;
	indexi1Fr2 = NoFr2 + MPIcount2F;
	indexi1Fr3 = NoFr3 + MPIcount2F;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids  + shiftx;
	/******************************/

	if(ox2 == 4){

		indexj1Er  = NoEr + 2 + MPIcount1 + count1z;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{
		
		jshearing = Ny * jproc + (j - js) + joffset;
		
		if(jshearing >= NGy * Ny ) {
			jshearing -= NGy * Ny;
		}		
		
		shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
		
		jshearing = jshearing - shearing_grid * Ny;
		
		shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
		
		
		col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
		
		
		
		sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
		
		coli1 = col_shear;
		
		if(ox3 == 4){
			if(col_shear < colk1)	sheark1 = 2;
			
			indexk1Er  = NoEr  + sheark1;
			indexk1Fr1 = NoFr1 + sheark1;
			indexk1Fr2 = NoFr2 + sheark1;
			indexk1Fr3 = NoFr3 + sheark1;
		}
		else {
			sheark1 = 2;
		}
	
		
		if(col_shear < colk0)	sheark0 = 2;
		if(col_shear < colj0)	shearj0 = 2;
		if(col_shear < coli0)	sheari0 = 2;
		if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
		if(col_shear < colj1)	shearj1 = 2;		
		if(col_shear < coli)	sheari1 = 2;
		
		
		
		indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
		indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	
		indexk0Er = NoEr + sheark0 + count1z;
		indexk0Fr1 = NoFr1 + sheark0 + count1z;
		indexk0Fr2 = NoFr2 + sheark0 + count1z;
		indexk0Fr3 = NoFr3 + sheark0 + count1z;
	
	
		indexj1Er  = NoEr + 2 + shearj1 + count1z;
		indexj1Fr1 = NoFr1 + 2 + shearj1 + count1z;
		indexj1Fr2 = NoFr2 + 2 + shearj1 + count1z;
		indexj1Fr3 = NoFr3 + 2 + shearj1 + count1z;
	
		indexj0Er  = NoEr  + 4 + shearj0 + count1z;
		indexj0Fr1 = NoFr1 + 4 + shearj0 + count1z;
		indexj0Fr2 = NoFr2 + 4 + shearj0 + count1z;
		indexj0Fr3 = NoFr3 + 4 + shearj0 + count1z;
		
				
		indexi0Er  = NoEr + 6 + sheari + count1z;
		indexi0Fr1 = NoFr1 + 6 + sheari + count1z;
		indexi0Fr2 = NoFr2 + 6 + sheari + count1z;
		indexi0Fr3 = NoFr3 + 6 + sheari + count1z;
		
		indexiEr  = NoEr + 8 + sheari + count1z;
		indexiFr1 = NoFr1 + 8 + sheari + count1z;
		indexiFr2 = NoFr2 + 8 + sheari + count1z;
		indexiFr3 = NoFr3 + 8 + sheari + count1z;
	
		
	}	
	
#endif
	

	return;

}

void ie_je_ke_phy_MPI_phy()
{

	int i, j, k, count1x, shifty, count1z;
	i = ie;
	j = je;
	k = ke;

	
	if(ox1 == 4) count1x = 2;
	else	count1x = 0;
	if(ox3 == 4) count1z = 2;
	else	count1z = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx2 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1x + count1z;
		MPIcount2F = 8 + count1x + count1z;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1 + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;	
	
	
	
	indexi0Er  = NoEr + 4 + MPIcount1 + count1x + count1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1x + count1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1x + count1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 6 + MPIcount1 + count1x + count1z;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1x + count1z;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1x + count1z;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1x + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	if(ox1 == 4){
		indexi1Er  = NoEr + 4 + MPIcount1 + count1z;
		indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1z;
		indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1z;
		indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1z;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	
	
	/***************************/
	/* for y boundary, MPI part */
	indexj1Er  = NoEr + MPIcount2;
	indexj1Fr1 = NoFr1 + MPIcount2F;
	indexj1Fr2 = NoFr2 + MPIcount2F;
	indexj1Fr3 = NoFr3 + MPIcount2F;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;
	/******************************/


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1;
		indexk1Fr1 = NoFr1 + MPIcount1;
		indexk1Fr2 = NoFr2 + MPIcount1;
		indexk1Fr3 = NoFr3 + MPIcount1;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx2 > ID){
		
		MPIcount2 = 10 + count1z;
		MPIcount2F = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	

	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1 + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1 + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1 + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1 + count1z;

	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1 + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1 + count1z;

	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1 + count1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1 + count1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1 + count1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1 + count1z;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1 + count1z;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1 + count1z;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1 + count1z;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1 + count1z;
	
	indexj1Er  = NoEr  + MPIcount2 + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2F + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2F + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2F + shearj1;
	
	
#endif


	return;

}

void ie_je_ke_MPI_MPI_phy()
{
	
	int i, j, k, MPIcount1y, MPIcount1x, MPIcount2y, MPIcount2x, MPIcount2Fy, MPIcount2Fx, shiftx, shifty, count1z;
	i = ie;
	j = je;
	k = ke;

	if(ox3 == 4) 	count1z = 2;
	else		count1z = 0;
		
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1z;
		MPIcount2Fx = 8 + count1z;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12 + count1z;
		MPIcount2Fy = 10 + count1z;
	}



		/* There are seven blocks */

	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1x + count1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1x + count1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	
	

	

	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1y + count1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexi0Er  = NoEr + 4 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1y + MPIcount1x + count1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1y + MPIcount1x + count1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 6 + MPIcount1x + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 6 + MPIcount1x + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 6 + MPIcount1x + MPIcount1y + count1z;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;


	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1y;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1y;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1y;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1y;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/


	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	if(ox3 == 4){
		indexk1Er  = NoEr + MPIcount1x + MPIcount1y;
		indexk1Fr1 = NoFr1 + MPIcount1x + MPIcount1y;
		indexk1Fr2 = NoFr2 + MPIcount1x + MPIcount1y;
		indexk1Fr3 = NoFr3 + MPIcount1x + MPIcount1y;
		colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox3 == 1 || ox3 == 5){
					
		theta[6] += theta[14];
		theta[9] -= theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] -= varphi[13];
	}
	else if(ox3 == 2){
		theta[6] += theta[14];
		theta[9] += theta[15];
			
		phi[6] += phi[12];
		phi[7] += phi[13];
	
		psi[6] += psi[12];
		psi[7] += psi[13];
			
		varphi[6] += varphi[12];
		varphi[7] += varphi[13];
	}
	else if(ox3 == 3){
			/* Do nothing */
	}	
	
	
	
	
	
#ifdef SHEARING_BOX
	
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx2 > ID){
		
		MPIcount2y = 10 + count1z;
		MPIcount2Fy = 8 + count1z;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	if(ox3 == 4){
		if(col_shear < colk1)	sheark1 = 2;
		
		indexk1Er  = NoEr  + sheark1 + MPIcount1y;
		indexk1Fr1 = NoFr1 + sheark1 + MPIcount1y;
		indexk1Fr2 = NoFr2 + sheark1 + MPIcount1y;
		indexk1Fr3 = NoFr3 + sheark1 + MPIcount1y;
	}
	else {
		sheark1 = 2;
	}
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;	
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + count1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + count1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + count1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + count1z;
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + count1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + count1z;
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1y + count1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1y + count1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1y + count1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1y + count1z;
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1y + count1z;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1y + count1z;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1y + count1z;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1y + count1z;
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1;
	
}	
#endif

	return;


}




void ie_je_ke_phy_phy_MPI()
{
	int i, j, k, count1y, count1x, shiftz;
	i = ie;
	j = je;
	k = ke;
	
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	if(rx3 < ID){	
		
		MPIcount1 = 2;
		MPIcount2 = 0;
		MPIcount2F = 0;
	}
	else{
		
		MPIcount1 = 0;
		MPIcount2 = 10 + count1y + count1x;
		MPIcount2F = 8 + count1y + count1x;
	}


	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1;
	indexk0Fr1 = NoFr1 + MPIcount1;
	indexk0Fr2 = NoFr2 + MPIcount1;
	indexk0Fr3 = NoFr3 + MPIcount1;
	colk0 = 4 * (k - ks -1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1 + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1 + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1 + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1 + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	indexi0Er  = NoEr + 4 + MPIcount1 + count1x + count1y;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1 + count1x + count1y;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1 + count1x + count1y;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1 + count1x + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 6 + MPIcount1 + count1x + count1y;
	indexiFr1 = NoFr1 + 6 + MPIcount1 + count1x + count1y;
	indexiFr2 = NoFr2 + 6 + MPIcount1 + count1x + count1y;
	indexiFr3 = NoFr3 + 6 + MPIcount1 + count1x + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){
		indexi1Er  = NoEr + 4 + MPIcount1 + count1y;
		indexi1Fr1 = NoFr1 + 4 + MPIcount1 + count1y;
		indexi1Fr2 = NoFr2 + 4 + MPIcount1 + count1y;
		indexi1Fr3 = NoFr3 + 4 + MPIcount1 + count1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;
	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for y boundary */
	if(ox2 == 4){
		indexj1Er  = NoEr + 2 + MPIcount1;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;
	}
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	


	/* MPI boundary z  direction */
	indexk1Er  = NoEr + MPIcount2;
	indexk1Fr1 = NoFr1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + MPIcount2F;
	indexk1Fr3 = NoFr3 + MPIcount2F;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2 = 12;
		MPIcount2F = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1;
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2F;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2F;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2F;
	
	
	
	
	indexj1Er  = NoEr + 2 + shearj1 + MPIcount1;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + MPIcount1;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + MPIcount1;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + MPIcount1;
	
	
	indexj0Er  = NoEr  + 4 + MPIcount1 + shearj0;
	indexj0Fr1 = NoFr1 + 4 + MPIcount1 + shearj0;
	indexj0Fr2 = NoFr2 + 4 + MPIcount1 + shearj0;
	indexj0Fr3 = NoFr3 + 4 + MPIcount1 + shearj0;
	
	

	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 6 + sheari + MPIcount1;
	indexi0Fr1 = NoFr1 + 6 + sheari + MPIcount1;
	indexi0Fr2 = NoFr2 + 6 + sheari + MPIcount1;
	indexi0Fr3 = NoFr3 + 6 + sheari + MPIcount1;
	
	indexiEr  = NoEr + 8 + sheari + MPIcount1;
	indexiFr1 = NoFr1 + 8 + sheari + MPIcount1;
	indexiFr2 = NoFr2 + 8 + sheari + MPIcount1;
	indexiFr3 = NoFr3 + 8 + sheari + MPIcount1;
	

	
	
#endif

	return;
}

void ie_je_ke_MPI_phy_MPI()
{

	int i, j, k, count1y, shiftz, shiftx;
	int MPIcount1x, MPIcount1z, MPIcount2x, MPIcount2z, MPIcount2Fx, MPIcount2Fz;
	i = ie;
	j = je;
	k = ke;

	if(ox2 == 4) count1y = 2;
	else	count1y = 0;

	
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx1 < ID){	
		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10 + count1y;
		MPIcount2Fx = 8 + count1y;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1y;
		MPIcount2Fz = 10 + count1y;
	}



	/* There are seven blocks */
	
	indexk0Er  = NoEr + MPIcount1x + MPIcount1z;
	indexk0Fr1 = NoFr1 + MPIcount1x + MPIcount1z;
	indexk0Fr2 = NoFr2 + MPIcount1x + MPIcount1z;
	indexk0Fr3 = NoFr3 + MPIcount1x + MPIcount1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z + count1y;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z + count1y;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	
	
	
	indexi0Er  = NoEr + 4 +  MPIcount1x + MPIcount1z + count1y;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1z + count1y;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1z + count1y;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 6 + MPIcount1x + MPIcount1z + count1y;
	indexiFr1 = NoFr1 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexiFr2 = NoFr2 + 6 + MPIcount1x + MPIcount1z + count1y;
	indexiFr3 = NoFr3 + 6 + MPIcount1x + MPIcount1z + count1y;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount2x + MPIcount1z;
	indexi1Fr1 = NoFr1 + MPIcount2Fx + MPIcount1z;
	indexi1Fr2 = NoFr2 + MPIcount2Fx + MPIcount1z;
	indexi1Fr3 = NoFr3 + MPIcount2Fx + MPIcount1z;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	
	if(ox2 == 4){

		indexj1Er  = NoEr + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr1 = NoFr1 + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr2 = NoFr2 + 2 + MPIcount1x + MPIcount1z;
		indexj1Fr3 = NoFr3 + 2 + MPIcount1x + MPIcount1z;
		colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox2 == 1 || ox2 == 5){
					
		theta[6] += theta[12];
		theta[8] -= theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] -= psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 2){
		theta[6] += theta[12];
		theta[8] += theta[13];
			
		phi[6] += phi[10];
		phi[7] += phi[11];
	
		psi[6] += psi[10];
		psi[7] += psi[11];
			
		varphi[6] += varphi[10];
		varphi[7] += varphi[11];
	}
	else if(ox2 == 3){
			/* Do nothing */
	}	



	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
	
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli0)	sheari0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1z;
	
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	
	
	indexj1Er  = NoEr + 2 + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + 2 + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + 2 + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + 2 + shearj1 + MPIcount1z;
	
	
	indexj0Er  = NoEr  + 4 + MPIcount1z + shearj0;
	indexj0Fr1 = NoFr1 + 4 + MPIcount1z + shearj0;
	indexj0Fr2 = NoFr2 + 4 + MPIcount1z + shearj0;
	indexj0Fr3 = NoFr3 + 4 + MPIcount1z + shearj0;
	
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexi0Er  = NoEr + 6 + sheari + MPIcount1z;
	indexi0Fr1 = NoFr1 + 6 + sheari + MPIcount1z;
	indexi0Fr2 = NoFr2 + 6 + sheari + MPIcount1z;
	indexi0Fr3 = NoFr3 + 6 + sheari + MPIcount1z;
	
	indexiEr  = NoEr + 8 + sheari + MPIcount1z;
	indexiFr1 = NoFr1 + 8 + sheari + MPIcount1z;
	indexiFr2 = NoFr2 + 8 + sheari + MPIcount1z;
	indexiFr3 = NoFr3 + 8 + sheari + MPIcount1z;
	
}
		
#endif

	return;

}

void ie_je_ke_phy_MPI_MPI()
{

	int i, j, k, count1x, shiftz, shifty;
	int MPIcount1y, MPIcount1z, MPIcount2y, MPIcount2z, MPIcount2Fy, MPIcount2Fz;
	i = ie;
	j = je;
	k = ke;

	if(ox1 == 4) count1x = 2;
	else	count1x = 0;

	
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);

	if(rx2 < ID){	
		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 10 + count1x;
		MPIcount2Fy = 8 + count1x;
	}

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 12 + count1x;
		MPIcount2Fz = 10 + count1x;
	}


	/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1z;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	


	
	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1z;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;

	
	
	
	indexi0Er  = NoEr + 4 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y + count1x;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y + count1x;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	


	indexiEr  = NoEr + 6 + MPIcount1z + MPIcount1y + count1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + count1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + count1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	if(ox1 == 4){

		indexi1Er  = NoEr + 4 + MPIcount1z + MPIcount1y;
		indexi1Fr1 = NoFr1 + 4 + MPIcount1z + MPIcount1y;
		indexi1Fr2 = NoFr2 + 4 + MPIcount1z + MPIcount1y;
		indexi1Fr3 = NoFr3 + 4 + MPIcount1z + MPIcount1y;
		coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids;

	}/* periodic for z direction */
	else if(ox1 == 1 || ox1 == 5){
					
		theta[6] += theta[10];
		theta[7] -= theta[11];
			
		phi[6] += phi[8];
		phi[7] -= phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 2){
		theta[6] += theta[10];
		theta[7] += theta[11];
			
		phi[6] += phi[8];
		phi[7] += phi[9];
	
		psi[6] += psi[8];
		psi[7] += psi[9];
			
		varphi[6] += varphi[8];
		varphi[7] += varphi[9];
	}
	else if(ox1 == 3){
			/* Do nothing */
	}	

	/***************************/
	/* for y boundary, MPI part */
	
	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids  + shifty;
	
	/******************************/



	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
	
#ifdef SHEARING_BOX
	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
		
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1z + MPIcount1y;


	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
	
	
#endif

	return;

}

void ie_je_ke_MPI_MPI_MPI()
{
	
	int i, j, k;
	int MPIcount1x, MPIcount1y, MPIcount1z, MPIcount2x, MPIcount2y, MPIcount2z, MPIcount2Fx, MPIcount2Fy, MPIcount2Fz, shiftx, shiftz, shifty;
	i = ie;
	j = je;
	k = ke;

		
	shiftz = 4 * Ny * Nx * Nz * (rx3 - ID);
	shifty = 4 * Ny * Nx * Nz * (rx2 - ID);
	shiftx = 4 * Ny * Nx * Nz * (rx1 - ID);
	

	if(rx3 < ID){	
		
		MPIcount1z = 2;
		MPIcount2z = 0;
		MPIcount2Fz = 0;
	}
	else{
		
		MPIcount1z = 0;
		MPIcount2z = 14;
		MPIcount2Fz = 12;
	}

	if(rx2 < ID){		
		MPIcount1y = 2;
		MPIcount2y = 0;
		MPIcount2Fy = 0;
	}
	else{
		
		MPIcount1y = 0;
		MPIcount2y = 12;
		MPIcount2Fy = 10;
	}

	if(rx1 < ID){		
		MPIcount1x = 2;
		MPIcount2x = 0;
		MPIcount2Fx = 0;
	}
	else{
		
		MPIcount1x = 0;
		MPIcount2x = 10;
		MPIcount2Fx = 8;
	}



		/* There are seven blocks */


	
	
	indexk0Er  = NoEr + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr1 = NoFr1 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr2 = NoFr2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexk0Fr3 = NoFr3 + MPIcount1y + MPIcount1z + MPIcount1x;
	colk0 = 4 * (k - ks - 1) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;
	

	

	indexj0Er  = NoEr + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr1 = NoFr1 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr2 = NoFr2 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	indexj0Fr3 = NoFr3 + 2 + MPIcount1y + MPIcount1z + MPIcount1x;
	colj0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js - 1) * Nx + 4 * (i - is) + count_Grids;
	

	
	
	indexi0Er  = NoEr + 4 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + MPIcount1x + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + MPIcount1x + MPIcount1y + MPIcount1z;
	coli0 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is - 1) + count_Grids;
	

	indexiEr  = NoEr + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr1 = NoFr1 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr2 = NoFr2 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	indexiFr3 = NoFr3 + 6 + MPIcount1z + MPIcount1y + MPIcount1x;
	coli = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids;

	/***************************/
	/* for x boundary, MPI part */
	indexi1Er  = NoEr + MPIcount1z + MPIcount1y + MPIcount2x;
	indexi1Fr1 = NoFr1 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr2 = NoFr2 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	indexi1Fr3 = NoFr3 + MPIcount1z + MPIcount1y + MPIcount2Fx;
	coli1 = 4 * (k - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (is - is) + count_Grids + shiftx;
	/******************************/

	/***************************/
	/* for y boundary, MPI part */

	indexj1Er  = NoEr + MPIcount1z + MPIcount2y;
	indexj1Fr1 = NoFr1 + MPIcount1z + MPIcount2Fy;
	indexj1Fr2 = NoFr2 + MPIcount1z + MPIcount2Fy;
	indexj1Fr3 = NoFr3 + MPIcount1z + MPIcount2Fy;
	colj1 = 4 * (k - ks) * Nx * Ny + 4 * (js - js) * Nx + 4 * (i - is) + count_Grids + shifty;

	/******************************/


	/***************************/
	/* for z boundary, MPI part */
	indexk1Er  = NoEr + MPIcount2z;
	indexk1Fr1 = NoFr1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + MPIcount2Fz;
	colk1 = 4 * (ks - ks) * Nx * Ny + 4 * (j - js) * Nx + 4 * (i - is) + count_Grids + shiftz;
	/******************************/
	
	
	
	
	
#ifdef SHEARING_BOX
if(rx1 < ID)
{	
	jshearing = Ny * jproc + (j - js) + joffset;
	
	if(jshearing >= NGy * Ny) {
		jshearing -= NGy * Ny;
	}		
	
	shearing_grid = (int)(jshearing/Ny); /* The integer part of grid */
	
	jshearing = jshearing - shearing_grid * Ny;
	
	shearing_grid = kproc * NGx * NGy + shearing_grid * NGx + iproc;
	
	
	col_shear = 4 * (k - ks) * Nx * Ny + 4 * jshearing * Nx + 4 * (is - is) + count_Grids + (shearing_grid - ID) * lines + shiftx;
	
	
	if(rx3 > ID){
		
		MPIcount2z = 12;
		MPIcount2Fz = 10;
	}
	
	if(rx2 > ID){
		
		MPIcount2y = 10;
		MPIcount2Fy = 8;
	}
	
	
	
	sheark0 = 0; shearj0 = 0; sheari0 = 0; sheari = 0; sheari1 = 0; shearj1 = 0; sheark1 = 0;
	
	coli1 = col_shear;
	
	
	if(col_shear < colk0)	sheark0 = 2;
	if(col_shear < colj0)	shearj0 = 2;
	if(col_shear < coli)	{sheari = 2; sheari0 = 2;}
	if(col_shear < colj1)	shearj1 = 2;
	if(col_shear < colk1)	sheark1 = 2;
	if(col_shear < coli)	sheari1 = 2;
	
	indexk1Er  = NoEr  + sheark1 + MPIcount2z;
	indexk1Fr1 = NoFr1 + sheark1 + MPIcount2Fz;
	indexk1Fr2 = NoFr2 + sheark1 + MPIcount2Fz;
	indexk1Fr3 = NoFr3 + sheark1 + MPIcount2Fz;
	
	
	indexk0Er = NoEr + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr1 = NoFr1 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr2 = NoFr2 + sheark0 + MPIcount1y + MPIcount1z;
	indexk0Fr3 = NoFr3 + sheark0 + MPIcount1y + MPIcount1z;
	
	
	indexj0Er  = NoEr + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr1 = NoFr1 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr2 = NoFr2 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	indexj0Fr3 = NoFr3 + 2 + shearj0 + MPIcount1y + MPIcount1z;
	
	
	indexi0Er  = NoEr + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr1 = NoFr1 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr2 = NoFr2 + 4 + sheari + MPIcount1y + MPIcount1z;
	indexi0Fr3 = NoFr3 + 4 + sheari + MPIcount1y + MPIcount1z;
	
	
	indexi1Er  = NoEr  + 14 - (sheark0 + shearj0 + sheari1 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr1 = NoFr1 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr2 = NoFr2 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	indexi1Fr3 = NoFr3 + 12 - (sheark0 + shearj0 + 2 * sheari + shearj1 + sheark1);
	
	indexiEr  = NoEr + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr1 = NoFr1 + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr2 = NoFr2 + 6 + sheari + MPIcount1z + MPIcount1y;
	indexiFr3 = NoFr3 + 6 + sheari + MPIcount1z + MPIcount1y;
	
	
	indexj1Er  = NoEr  + MPIcount2y + shearj1 + MPIcount1z;
	indexj1Fr1 = NoFr1 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr2 = NoFr2 + MPIcount2Fy + shearj1 + MPIcount1z;
	indexj1Fr3 = NoFr3 + MPIcount2Fy + shearj1 + MPIcount1z;
	
}	
	
#endif
	

	return;

}




#endif /* radMHD_INTEGRATOR */

