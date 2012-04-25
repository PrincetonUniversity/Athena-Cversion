#include "copyright.h"

/*==============================================================================
 * FILE: radtest.c
 *
 * PURPOSE:  Problem generator for a non-LTE test of radiative transfer routine
 *           assuming a 1D variation of the optical depth.  Can be run for
 *           any number of dimensions with periodic boundary conditions in
 *           the non-varying directions.  Directions of variation is set
 *           by vert_dir=1,2,3 in the problem block.  eps sets the degree of
 *           deviation from LTE.  See Fabiani Bendicho & Trujillo Bueno ApJ,
 *           455, 646..
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

static Real eps0;
/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *============================================================================*/
static Real const_B(const GridS *pG, const int ifr, const int i,
		    const int j, const int k);
static Real const_eps(const GridS *pG, const int ifr, const int i,
		      const int j, const int k);
static Real const_opacity(const GridS *pG, const int ifr, const int i,
			  const int j, const int k);

void problem(DomainS *pDomain)
{
  RadGridS *pRG = (pDomain->RadGrid);
  GridS *pG = (pDomain->Grid);
  int il = pG->is, iu = pG->ie;
  int jl = pG->js, ju = pG->je;
  int kl = pG->ks, ku = pG->ke;
  int nf=pRG->nf, nang=pRG->nang;
  int noct = pRG->noct;
  int i, j, k, ifr, l, m;
  int vdir;
  Real y, ytop, ybtm;  
  Real den = 1.0;
  Real *tau = NULL, taumax, taumin;

/* Read problem parameters. */

  eps0 = par_getd("problem","eps");
  vdir = par_geti("problem","vert_dir");
  taumax = par_getd("problem","taumax");
  taumin = par_getd("problem","taumin");
  R_ideal = 1.0;

/* ---------- Initialize Grid ------------ */

/* Setup density structure */ 
/* tau il used to initialize density grid */

  switch(vdir) {
  
  case 1:
    ytop = pDomain->RootMaxX[0];
    ybtm = pDomain->RootMinX[0];
    if ((tau = (Real *)calloc_1d_array(pG->Nx[0]+2*nghost,sizeof(Real))) == NULL) {
      ath_error("[problem]: Error allocating memory");
    }
    for(i=il; i<=iu+2; i++) {
      y = pG->MinX[0] + (Real)(i-il)*pG->dx1;
      tau[i] = pow(10.0,taumin + (taumax-taumin) * ((y-ybtm)/(ytop-ybtm)));
    }
    il -= 1; iu += 1;
    if (pG->Nx[1] > 1) {
      jl -= 1; ju += 1;
    } else if (pG->Nx[2] > 1) {
      kl -= 1; ku += 1;
    }
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (i=il; i<=iu; i++) {
	  pG->U[k][j][i].d  = (tau[i+1] - tau[i]) / pG->dx1;
	  pG->U[k][j][i].E = 1.0;  /* needed for tgas init */
	}}}

    break;

  case 2:
    ytop = pDomain->RootMaxX[1];
    ybtm = pDomain->RootMinX[1];
    if ((tau = (Real *)calloc_1d_array(pG->Nx[1]+2*nghost,sizeof(Real))) == NULL) {
      ath_error("[problem]: Error allocating memory");
    }
    for(j=jl; j<=ju+2; j++) {
      y = pG->MinX[1] + (Real)(j-jl)*pG->dx2;
      tau[j] = pow(10.0,taumin + (taumax-taumin) * ((y-ybtm)/(ytop-ybtm)));
    }
    il -= 1; iu += 1;
    jl -= 1; ju += 1;
    if (pG->Nx[2] > 1) {
      kl -= 1; ku += 1;
    }
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (i=il; i<=iu; i++) {
	  pG->U[k][j][i].d  = (tau[j+1] - tau[j]) / pG->dx2;
	  pG->U[k][j][i].E = 1.0;  /* needed for tgas init */
	}}}
    
    break;

 case 3:
    ytop = pDomain->RootMaxX[2];
    ybtm = pDomain->RootMinX[2];
    if ((tau = (Real *)calloc_1d_array(pG->Nx[2]+2*nghost,sizeof(Real))) == NULL) {
      ath_error("[problem]: Error allocating memory");
    }
    for(k=kl; k<=ku+2; k++) {
      y = pG->MinX[2] + (Real)(k-kl)*pG->dx3;
      tau[k] = pow(10.0,taumin + (taumax-taumin) * ((y-ybtm)/(ytop-ybtm)));
    }
    il -= 1; iu += 1;
    jl -= 1; ju += 1;
    kl -= 1; ku += 1;
    for (k=kl; k<=ku; k++) {
      for (j=jl; j<=ju; j++) {
	for (i=il; i<=iu; i++) {
	  pG->U[k][j][i].d  = (tau[k+1] - tau[k]) / pG->dx3;
	  pG->U[k][j][i].E = 1.0;  /* needed for tgas init */
	}}}
    break;

  default:
    ath_error("[rad2d]: direction vert_dir must be 1-3\n");
    break;
  }
/* Free up memory */
  free_1d_array(tau);

/* ---------- Initialize RadGrid ------------ */

  il = pRG->is-1, iu = pRG->ie+1;
  jl = pRG->js,   ju = pRG->je;
  kl = pRG->ks,   ku = pRG->ke;
  if (pRG->Nx[1] > 1) { jl -= 1; ju += 1; }
  if (pRG->Nx[2] > 1) { kl -= 1; ku += 1; }

/* Initialize mean intensity */
  for(ifr=0; ifr<nf; ifr++)
    for (k=kl; k<=ku; k++)
      for (j=jl; j<=ju; j++)
	for(i=il; i<=iu; i++) {
	  pRG->R[ifr][k][j][i].J = 1.0;
	}
 
/* ------- Initialize boundary emission ---------------------------------- */

  switch(vdir) {

  case 1:
/* Density gradient aligned with i3 */
    for(ifr=0; ifr<nf; ifr++) {

      /* Initialize J to zero in the top boundary gridzones */
      for(k=kl; k<=ku; k++) {
	for(j=jl; j<=ju; j++) {
	  pRG->R[ifr][k][j][0].J = 0.0;
	}}
      /* Initialize boundary intensity in x1 direction */
      for(k=kl; k<=ku; k++) {
	for(j=jl; j<=ju; j++) {
	  for(m=0; m<nang; m++) {
	    /* lower boundary is tau=0, no irradiation */
	    pRG->Ghstl1i[ifr][k][j][0][m] = 0.0;
	    if (noct > 2) {
	      pRG->Ghstl1i[ifr][k][j][2][m] = 0.0;
	      if (noct == 8) {
		pRG->Ghstl1i[ifr][k][j][4][m] = 0.0;
		pRG->Ghstl1i[ifr][k][j][6][m] = 0.0;
	      }
	    }
	    /* upper boundary is large tau, eps=1 */
	    pRG->Ghstr1i[ifr][k][j][1][m] = 1.0;
	    if (noct > 2) {
	      pRG->Ghstr1i[ifr][k][j][3][m] = 1.0;
	      if (noct == 8) {
		pRG->Ghstr1i[ifr][k][j][5][m] = 1.0;
		pRG->Ghstr1i[ifr][k][j][7][m] = 1.0;
	      }
	    }
	  }}

	if (noct > 2) {
/* Initialize boundary intensity in x2 direction */
	  /* lower boundary is tau=0, no irradiation */
	  for(l=0; l<noct; l++) {
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr2i[ifr][k][il][l][m] = 0.0; 
	      pRG->Ghstl2i[ifr][k][il][l][m] = 0.0;
	    }}
	  for(i=il+1; i<=iu-1; i++) {
	    /* periodic radiation at left boundary */
	    for(m=0; m<nang; m++) {
	      pRG->Ghstl2i[ifr][k][i][0][m] = 1.0;
	      pRG->Ghstl2i[ifr][k][i][1][m] = 1.0;
	      if (noct == 8) {
		pRG->Ghstl2i[ifr][k][i][4][m] = 1.0;
		pRG->Ghstl2i[ifr][k][i][5][m] = 1.0;
	      }
	    }
	    /* periodic radiation at right boundary */
	    for(m=0; m<=nang; m++) {
	      pRG->Ghstr2i[ifr][k][i][2][m] = 1.0;
	      pRG->Ghstr2i[ifr][k][i][3][m] = 1.0;
	      if (noct == 8) {
		pRG->Ghstr2i[ifr][k][i][6][m] = 1.0;
		pRG->Ghstr2i[ifr][k][i][7][m] = 1.0;
	      }
	    }
	  }
	  /* upper boundary is large tau, eps=1 */
	  for(l=0; l<noct; l++) {
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr2i[ifr][k][iu][l][m] = 1.0; 
	      pRG->Ghstl2i[ifr][k][iu][l][m] = 1.0;
	    }}
	}
      }
/* Initialize boundary intensity in x3 direction */
      if (noct == 8) {
	for(j=jl; j<=ju; j++) {
	  /* lower boundary is tau=0, no irradiation */
	  for(l=0; l<noct; l++) { 
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr3i[ifr][j][il][l][m] = 0.0; 
	      pRG->Ghstl3i[ifr][j][il][l][m] = 0.0;
	    }}
	  for(i=il+1; i<=iu-1; i++) {
	    for(m=0; m<nang; m++) {
	      /* periodic radiation at left boundary */
	      pRG->Ghstl3i[ifr][j][i][0][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][1][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][2][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][3][m] = 1.0;
	      /* periodic radiation at right boundary */
	      pRG->Ghstr3i[ifr][j][i][4][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][5][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][6][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][7][m] = 1.0;
	    }}
	  /* upper boundary is large tau, eps=1 */
	  for(l=0; l<noct; l++) { 
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr3i[ifr][j][iu][l][m] = 1.0; 
	      pRG->Ghstl3i[ifr][j][iu][l][m] = 1.0;
	    }}
	}
      }
    }
    break;

  case 2:
/* Density gradient aligned with i2 */
    for(ifr=0; ifr<nf; ifr++) {

/* Initialize J to zero in the top boundary gridzones */
      for(k=kl; k<=ku; k++) {
	for(i=il; i<=iu; i++) {
	  pRG->R[ifr][k][0][i].J = 0.0;
	}}

/* Initialize boundary intensity in x1 direction */
      for(k=kl; k<=ku; k++) {
	/* lower boundary is tau=0, no irradiation */
	for(l=0; l<noct; l++) 
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][jl][l][m] = 0.0;
	    pRG->Ghstl1i[ifr][k][jl][l][m] = 0.0;
	  }
	for(j=jl+1; j<=ju-1; j++) {
	  /* periodic radiation at left boundary */
	  for(m=0; m<nang; m++) {
	    pRG->Ghstl1i[ifr][k][j][0][m] = 1.0;
	    pRG->Ghstl1i[ifr][k][j][2][m] = 1.0;
	    if (noct == 8) {
	      pRG->Ghstl1i[ifr][k][j][4][m] = 1.0;
	      pRG->Ghstl1i[ifr][k][j][6][m] = 1.0;
	    }
	  }
	  /* periodic radiation at right boundary */
	  for(m=0; m<=nang; m++) {
	    pRG->Ghstr1i[ifr][k][j][1][m] = 1.0;
	    pRG->Ghstr1i[ifr][k][j][3][m] = 1.0;
	    if (noct == 8) {
	      pRG->Ghstr1i[ifr][k][j][5][m] = 1.0;
	      pRG->Ghstr1i[ifr][k][j][7][m] = 1.0;
	    }
	  }
	}
	/* upper boundary is large tau, eps=1 */
      	for(l=0; l<noct; l++) 
	  for(m=0; m<nang; m++) {
	    pRG->Ghstr1i[ifr][k][ju][l][m] = 1.0;
	    pRG->Ghstl1i[ifr][k][ju][l][m] = 1.0;
	  }

/* Initialize boundary intensity in x2 direction */
	for(i=il; i<=iu; i++) { 
	  for(m=0; m<nang; m++) {
	    /* lower boundary is tau=0, no irradiation */
	    pRG->Ghstl2i[ifr][k][i][0][m] = 0.0;
	    pRG->Ghstl2i[ifr][k][i][1][m] = 0.0;
	    if (noct == 8) {
	      pRG->Ghstl2i[ifr][k][i][4][m] = 0.0;
	      pRG->Ghstl2i[ifr][k][i][5][m] = 0.0;
	    }
	    /* upper boundary is large tau, eps=1 */
	    pRG->Ghstr2i[ifr][k][i][2][m] = 1.0;
	    pRG->Ghstr2i[ifr][k][i][3][m] = 1.0;
	    if (noct == 8) {
	      pRG->Ghstr2i[ifr][k][i][6][m] = 1.0;
	      pRG->Ghstr2i[ifr][k][i][7][m] = 1.0;
	    }
	  }
	}
      }
/* Initialize boundary intensity in x3 direction */
      if (noct == 8) {
	for(i=il; i<=iu; i++) {
	  /* lower boundary is tau=0, no irradiation */
	  for(l=0; l<noct; l++) { 
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr3i[ifr][jl][i][l][m] = 0.0; 
	      pRG->Ghstl3i[ifr][jl][i][l][m] = 0.0;
	    }}}
	for(j=jl+1; j<=ju-1; j++) {
	  for(i=il; i<=iu; i++) {
	    for(m=0; m<nang; m++) {
	      /* periodic radiation at left boundary */
	      pRG->Ghstl3i[ifr][j][i][0][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][1][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][2][m] = 1.0;
	      pRG->Ghstl3i[ifr][j][i][3][m] = 1.0;
	      /* periodic radiation at right boundary */
	      pRG->Ghstr3i[ifr][j][i][4][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][5][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][6][m] = 1.0;
	      pRG->Ghstr3i[ifr][j][i][7][m] = 1.0;
	    }}}
	for(i=il; i<=iu; i++) {
	  /* upper boundary is large tau, eps=1 */
	  for(l=0; l<noct; l++) { 
	    for(m=0; m<nang; m++) {
	      pRG->Ghstr3i[ifr][ju][i][l][m] = 1.0; 
	      pRG->Ghstl3i[ifr][ju][i][l][m] = 1.0;
	    }}}
      }
    }
    break;

 case 3:
/* Density gradient aligned with i3 */
   for(ifr=0; ifr<nf; ifr++) {
/* Initialize J to zero in the top boundary gridzones */
     for(j=jl; j<=ju; j++) {
       for(i=il; i<=iu; i++) {
	 pRG->R[ifr][0][j][i].J = 0.0;
       }}     
     /* Initialize boundary intensity in x1 direction */
     /* lower boundary is tau=0, no irradiation */
     for(j=jl; j<=ju; j++) {
       for(l=0; l<8; l++) { 
	 for(m=0; m<nang; m++) {
	   pRG->Ghstr1i[ifr][kl][j][l][m] = 0.0;
	   pRG->Ghstl1i[ifr][kl][j][l][m] = 0.0;
	 }}}
     for(k=kl+1; k<=ku-1; k++) {
       for(j=jl; j<=ju; j++) {
	 for(m=0; m<nang; m++) {
	   /* periodic radiation at left boundary */
	   pRG->Ghstl1i[ifr][k][j][0][m] = 1.0;
	   pRG->Ghstl1i[ifr][k][j][2][m] = 1.0;
	   pRG->Ghstl1i[ifr][k][j][4][m] = 1.0;
	   pRG->Ghstl1i[ifr][k][j][6][m] = 1.0;
	   /* periodic radiation at right boundary */
	   pRG->Ghstr1i[ifr][k][j][1][m] = 1.0;
	   pRG->Ghstr1i[ifr][k][j][3][m] = 1.0;
	   pRG->Ghstr1i[ifr][k][j][5][m] = 1.0;
	   pRG->Ghstr1i[ifr][k][j][7][m] = 1.0;
	 }}}
     /* upper boundary is large tau, eps=1 */
     for(j=jl; j<=ju; j++) {
       for(l=0; l<8; l++) { 
	 for(m=0; m<nang; m++) {
	   pRG->Ghstr1i[ifr][ku][j][l][m] = 1.0;
	   pRG->Ghstl1i[ifr][ku][j][l][m] = 1.0;
	 }}}

     /* Initialize boundary intensity in x2 direction */
     /* lower boundary is tau=0, no irradiation */
     for(i=il; i<=iu; i++) {
       for(l=0; l<8; l++) { 
	 for(m=0; m<nang; m++) {
	   pRG->Ghstr2i[ifr][kl][i][l][m] = 0.0;
	   pRG->Ghstl2i[ifr][kl][i][l][m] = 0.0;
	 }}}
     for(k=kl+1; k<=ku-1; k++) {
       for(i=il; i<=iu; i++) { 
	 for(m=0; m<nang; m++) {
	   /* periodic radiation at left boundary */     
	   pRG->Ghstl2i[ifr][k][i][0][m] = 1.0;
	   pRG->Ghstl2i[ifr][k][i][1][m] = 1.0;
	   pRG->Ghstl2i[ifr][k][i][4][m] = 1.0;
	   pRG->Ghstl2i[ifr][k][i][5][m] = 1.0;
	   /* periodic radiation at right boundary */
	   pRG->Ghstr2i[ifr][k][i][2][m] = 1.0;
	   pRG->Ghstr2i[ifr][k][i][3][m] = 1.0;
	   pRG->Ghstr2i[ifr][k][i][6][m] = 1.0;
	   pRG->Ghstr2i[ifr][k][i][7][m] = 1.0;
	 }}}
     /* upper boundary is large tau, eps=1 */
     for(i=il; i<=iu; i++) {
       for(l=0; l<8; l++) { 
	 for(m=0; m<nang; m++) {
	   pRG->Ghstr2i[ifr][ku][i][l][m] = 1.0;
	   pRG->Ghstl2i[ifr][ku][i][l][m] = 1.0;
	 }}}

     /* Initialize boundary intensity in x3 direction */
     for(j=jl; j<=ju; j++) {
       for(i=il; i<=iu; i++) {
	  for(m=0; m<nang; m++) {
	    /* lower boundary is tau=0, no irradiation */
	    pRG->Ghstl3i[ifr][j][i][0][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][1][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][2][m] = 0.0;
	    pRG->Ghstl3i[ifr][j][i][3][m] = 0.0;
	    /* upper boundary is large tau, eps=1 */
	    pRG->Ghstr3i[ifr][j][i][4][m] = 1.0;
	    pRG->Ghstr3i[ifr][j][i][5][m] = 1.0;
	    pRG->Ghstr3i[ifr][j][i][6][m] = 1.0;
	    pRG->Ghstr3i[ifr][j][i][7][m] = 1.0;
	  }
       }}
    }
    break;

  }
/* enroll radiation specification functions */
get_thermal_source = const_B;
get_thermal_fraction = const_eps;
get_total_opacity = const_opacity;

  return;
}

/*==============================================================================
 * PUBLIC PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 * Userwork_in_formal_solution  - problem specific work in formal solution loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}


void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_formal_solution(DomainS *pD)
{

  RadGridS *pRG=(pD->RadGrid);
  int i,j,k;
  int ifr = 0;
  int is=pRG->is, ie=pRG->ie;
  int js=pRG->js, je=pRG->je;
  int ks=pRG->ks, ke=pRG->ke;
  int ixmax, iymax, izmax;
  Real ds;
  static Real ***dst = NULL, *tau0 = NULL, ***sol = NULL;
  static int frstflag = 1, iter = 0;
  int vdir;
  Real jsol, chio, chim, chip, dtaum, dtaup;
  FILE *fp;
  char *fname;

  vdir = par_geti("problem","vert_dir");
  if (frstflag == 1) {
    if ((dst = (Real ***)calloc_3d_array(pRG->Nx[2]+2,pRG->Nx[1]+2,pRG->Nx[0]+2,
					 sizeof(Real))) == NULL) {
      ath_error("[Userwork_in_formal_solution]: Error allocating memory\n");
    }
    if ((sol = (Real ***)calloc_3d_array(pRG->Nx[2],pRG->Nx[1],pRG->Nx[0],
					 sizeof(Real))) == NULL) {
      ath_error("[Userwork_in_formal_solution]: Error allocating memory\n");
    }

    switch(vdir) {

    case 1:
      if ((tau0 = (Real *)calloc_1d_array(pRG->Nx[0]+2,sizeof(Real))) == NULL) {
	ath_error("[problem]: Error allocating memory");
      }
      tau0[0]= 0.0;
      for (i=is; i<=ie; i++) {
	chio = pRG->R[ifr][ks][js][i  ].chi;
	chim = pRG->R[ifr][ks][js][i-1].chi;
	chip = pRG->R[ifr][ks][js][i+1].chi;
	/*dtaum = 0.5 * (chim + chio);
	/*dtaup = 0.5 * (chip + chio);*/
	interp_quad_chi(chim,chio,chip,&dtaum);
	interp_quad_chi(chip,chio,chim,&dtaup);
	dtaum *= pRG->dx1; 
	dtaup *= pRG->dx1;
	tau0[i] = tau0[i-1] + dtaum;
      }
      for(k=ks; k<=ke; k++) {
	for(j=js; j<=je; j++) {
	  for(i=is; i<=ie; i++) {
	    jsol = 1.0 - exp(-sqrt(3.0 * eps0) * tau0[i]) / (1.0 + sqrt(eps0));
	    sol[k-ks][j-js][i-is] = eps0 + (1.0-eps0) * jsol;
	  }}}
      break;

    case 2:
      if ((tau0 = (Real *)calloc_1d_array(pRG->Nx[1]+2,sizeof(Real))) == NULL) {
	ath_error("[problem]: Error allocating memory");
      }
      tau0[0]= 0.0;
      for (j=js; j<=je; j++) {
	chio = pRG->R[ifr][ks][j  ][is].chi;
	chim = pRG->R[ifr][ks][j-1][is].chi;
	chip = pRG->R[ifr][ks][j+1][is].chi;
	/*dtaum = 0.5 * (chim + chio);
	/*dtaup = 0.5 * (chip + chio);*/
	interp_quad_chi(chim,chio,chip,&dtaum);
	interp_quad_chi(chip,chio,chim,&dtaup);
	dtaum *= pRG->dx2; 
	dtaup *= pRG->dx2;
	tau0[j] = tau0[j-1] + dtaum;
      }
      for(k=ks; k<=ke; k++) {
	for(j=js; j<=je; j++) {
	  for(i=is; i<=ie; i++) {
	    jsol = 1.0 - exp(-sqrt(3.0 * eps0) * tau0[j]) / (1.0 + sqrt(eps0));
	    sol[k-ks][j-js][i-is] = eps0 + (1.0-eps0) * jsol;
	  }}}
      break;
      
    case 3:
      if ((tau0 = (Real *)calloc_1d_array(pRG->Nx[2]+2,sizeof(Real))) == NULL) {
	ath_error("[problem]: Error allocating memory");
      }
      tau0[0]= 0.0;
      for (k=ks; k<=ke; k++) {
	chio = pRG->R[ifr][k  ][js][is].chi;
	chim = pRG->R[ifr][k-1][js][is].chi;
	chip = pRG->R[ifr][k+1][js][is].chi;
	/*dtaum = 0.5 * (chim + chio);
	/*dtaup = 0.5 * (chip + chio);*/
	interp_quad_chi(chim,chio,chip,&dtaum);
	interp_quad_chi(chip,chio,chim,&dtaup);
	dtaum *= pRG->dx3; 
	dtaup *= pRG->dx3;
	tau0[k] = tau0[k-1] + dtaum;
      }
      for(k=ks; k<=ke; k++) {
	for(j=js; j<=je; j++) {
	  for(i=is; i<=ie; i++) {
	    jsol = 1.0 - exp(-sqrt(3.0 * eps0) * tau0[k]) / (1.0 + sqrt(eps0));
	    sol[k-ks][j-js][i-is] = eps0 + (1.0-eps0) * jsol;
	  }}}
      break;

    }
      frstflag = 0;
  }

  iter++;
  printf("Iteration # %d \n",iter);
  ds = 0;
  for(k=ks; k<=ke; k++) {
    for(j=js; j<=je; j++) {
      for(i=is; i<=ie; i++) {
	dst[k][j][i] = fabs(pRG->R[ifr][k][j][i].S - sol[k-ks][j-js][i-is]) / 
	  sol[k-ks][j-js][i-is];
	if (dst[k][j][i] > ds) {ds = dst[k][j][i]; ixmax = i; iymax = j; izmax = k;}
      }}}

/* Print error to file "Radtest-diag.0.dat"  */

  fname = ath_fname(NULL,"RadTest-diag",NULL,NULL,1,0,NULL,"dat");
/* The file exists -- reopen the file in append mode */
  if((fp=fopen(fname,"r")) != NULL){
    if((fp = freopen(fname,"a",fp)) == NULL){
      ath_error("[Userwork_in_formal_solution]: Unable to reopen file.\n");
      return;
    }
  }
/* The file does not exist -- open the file in write mode */
  else{
    if((fp = fopen(fname,"w")) == NULL){
      ath_error("[Userwork_in_formal_solution]: Unable to open file.\n");
      return;
    }
/* Now write out some header information */
    fprintf(fp,"#  dSmax  ixmax  iymax  izmax\n");
  }
 
  fprintf(fp,"%d %g %d %d %d\n",iter,ds,ixmax,iymax,izmax);
  fclose(fp);

/* Print error to file "Radtest-error.0.dat"  */
  /*fname = ath_fname(NULL,"RadTest-error",NULL,NULL,1,0,NULL,"dat");
  if((fp=fopen(fname,"w")) != NULL){
    switch(vdir) {
  
    case 1:
      fprintf(fp,"%d %d\n",vdir,pRG->Nx[0]);
      for(i=pRG->is; i<=pRG->ie; i++) {    
	    fprintf(fp,"%g %g %g %g\n",tau0[i],dst[ks][js][i],pRG->R[ifr][ks][js][i].S,sol[0][0][i-is]);
      }
      break;
      
    case 2:
      fprintf(fp,"%d %d\n",vdir,pRG->Nx[1]);
      for(j=pRG->js; j<=pRG->je; j++) {    
	    fprintf(fp,"%g %g %g %g\n",tau0[j],dst[ks][j][is],pRG->R[ifr][ks][j][is].S,sol[0][j-js][0]);
      }
      break;

    case 3:
      fprintf(fp,"%d %d\n",vdir,pRG->Nx[2]);
      for(k=pRG->ks; k<=pRG->ke; k++) {    
	    fprintf(fp,"%g %g %g %g\n",tau0[k],dst[k][js][is],pRG->R[ifr][k][js][is].S,sol[k-ks][0][0]);
      }
      break;
    }

    fclose(fp);
  }*/

 
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

static Real const_B(const GridS *pG, const int ifr, const int i,
		    const int j, const int k)
{
  return 1.0;
}
static Real const_eps(const GridS *pG, const int ifr, const int i, 
		      const int j, const int k)
{

  return eps0;
  
}
static Real const_opacity(const GridS *pG, const int ifr, const int i,
			  const int j, const int k)
{

  return pG->U[k][j][i].d;
  
}

