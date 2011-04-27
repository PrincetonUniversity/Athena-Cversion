#include "../copyright.h"
/*==============================================================================
 * FILE: hydro_to_rad.c
 *
 * PURPOSE: Computes radiation variables (thermal source function,
 *          thermalization parameter, and total opacity) based on
 *          state of conserved variables after last integration
 *          step.  Currently uses simple constructions for testing
 *          purposes.
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 *   hydro_to_rad()
 *============================================================================*/

#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "../prototypes.h"
#define BFRAC 1.0E2

#ifdef RADIATION_TRANSFER
void output_rad_1d(RadGridS *pRG);
static char *construct_filename(char *basename,char *key,int dump,char *ext);
static int lte;

/*----------------------------------------------------------------------------*/
/* hydro_to_rad:  */

void hydro_to_rad(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k,ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;
  Real eps;
  Real etherm, ekin, B, d;

  lte = par_geti("radiation","lte");

  if (pG->Nx[0] > 1) {
    ioff = nghost - 1; il--; iu++;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1; jl--; ju++; 
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1; kl--; ku++;
  } else koff = 0;

/* Compute radiation variables from conserved variables */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
	ig = i + ioff;

	/* Comput gas temperature and store for later use */
	d = pG->U[kg][jg][ig].d;
	ekin =  pG->U[kg][jg][ig].M1 * pG->U[kg][jg][ig].M1;
	ekin += pG->U[kg][jg][ig].M2 * pG->U[kg][jg][ig].M2;
	ekin += pG->U[kg][jg][ig].M3 * pG->U[kg][jg][ig].M3;
	ekin *= 0.5 / d;
	etherm=pG->U[kg][jg][ig].E - ekin;
	pG->tgas[kg][jg][ig] = etherm * Gamma_1 / (d * R_ideal);

	for(ifr=0; ifr<nf; ifr++) {
#if defined(RADIATION_HYDRO) || defined(RADIATION_MHD)
	  if (lte == 1) 
	    pRG->R[k][j][i][ifr].J = pG->U[kg][jg][ig].Er / (4.0 * PI);
	  eps = get_thermal_fraction(pG,ifr,ig,jg,kg);
	  pRG->R[k][j][i][ifr].B = (1.0 - eps) * pRG->R[k][j][i][ifr].J +
	    eps  * get_thermal_source(pG,ifr,ig,jg,kg);
	  pRG->R[k][j][i][ifr].eps = 1.0;
	  pRG->R[k][j][i][ifr].S = pRG->R[k][j][i][ifr].B;
#else
	  eps = get_thermal_fraction(pG,ifr,ig,jg,kg);
	  pRG->R[k][j][i][ifr].B = get_thermal_source(pG,ifr,ig,jg,kg);
	  pRG->R[k][j][i][ifr].eps = eps;
	  pRG->R[k][j][i][ifr].S = (1.0 - eps) * pRG->R[k][j][i][ifr].J +
	                                  eps  * pRG->R[k][j][i][ifr].B;
#endif
	  pRG->R[k][j][i][ifr].chi = get_total_opacity(pG,ifr,ig,jg,kg);
	}
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* rad_to_hydro:  */

void rad_to_hydro(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  RadGridS *pRG=(pD->RadGrid);
  int i,j,k, ifr;
  int il = pRG->is, iu = pRG->ie;
  int jl = pRG->js, ju = pRG->je;
  int kl = pRG->ks, ku = pRG->ke;
  int nf = pRG->nf;
  int ig,jg,kg,ioff,joff,koff;
  Real esource, kappa;
  Real Delta, kDelta;
  Real dx = pRG->dx1;

  if (pG->Nx[0] > 1) {
    ioff = nghost - 1;
  } else ioff = 0;
  if (pG->Nx[1] > 1) {
    joff = nghost - 1; 
  } else joff = 0; 
  if (pG->Nx[2] > 1) {
    koff = nghost - 1;
  } else koff = 0;

/* Update thermal energy */
  for (k=kl; k<=ku; k++) {
    kg = k + koff;
    for (j=jl; j<=ju; j++) {
      jg = j + joff;
      for (i=il; i<=iu; i++) {
	ig = i + ioff;
	esource=0.0;
	for(ifr=0; ifr<nf; ifr++) {
	  kappa = pRG->R[k][j][i][ifr].eps * pRG->R[k][j][i][ifr].chi;
	  /* Must add frequency weights */
	  Delta = pRG->R[k][j][i][ifr].J - pRG->R[k][j][i][ifr].B;
	  /* if (fabs(Delta/pRG->R[k][j][i][ifr].B) >= BFRAC) { */
	  esource += kappa * Delta;
	  
	  /*if(i == il) 
	    esource += (pRG->R[k][j][i][ifr].H[0]-pRG->R[k][j][i+1][ifr].H[0])/dx; 
	  else if(i == iu)
	    esource += (pRG->R[k][j][i-1][ifr].H[0]-pRG->R[k][j][i][ifr].H[0])/dx;
	  else
	  esource += 0.5*(pRG->R[k][j][i-1][ifr].H[0]-pRG->R[k][j][i+1][ifr].H[0])/dx;*/

	  /* esource += kappa * (B00 - pRG->R[k][j][i][ifr].B) * (1.0 - kappa/(2.0*PI) * 
	     atan((2.0 * PI)/kappa)); */

	}
	pG->U[kg][jg][ig].E += pG->dt * 4.0 * PI * esource;
      }
    }
  }

  return;
}

void output_rad_1d(RadGridS *pRG)
{
  FILE *fp;
  char *fname;
  int i,l;
  static int itr;

  fname = construct_filename("rad", NULL, itr, "out");

  if ((fp = fopen(fname, "w")) == NULL) {
    ath_error("## Error opening radxxxx.out\n");
  }
  for(i=pRG->is; i<=pRG->ie; i++) {
    fprintf(fp,"%d %g %g %g %g %g\n",i,pRG->R[pRG->ks][pRG->js][i][0].B,
	    pRG->R[pRG->ks][pRG->js][i][0].J,pRG->R[pRG->ks][pRG->js][i][0].S,
	    pRG->R[pRG->ks][pRG->js][i][0].eps,pRG->R[pRG->ks][pRG->js][i][0].chi);
  }
  itr++;
  fclose(fp);
  return;
}

static char *construct_filename(char *basename,char *key,int dump,char *ext)
{
  char *fname = NULL;

  int namelen = strlen(basename)+1+4+1+strlen(ext)+1;
  if (key != NULL) namelen += 1+strlen(key);

  if ((fname = (char *)calloc(namelen,sizeof(char))) == NULL)
    ath_error("[problem]: Error allocating memory for filename\n");

  if (key != NULL) {
    sprintf(fname, "%s-%s.%04d.%s", basename, key, dump, ext);
  } else {
    sprintf(fname, "%s.%04d.%s", basename, dump, ext);
  }

  return fname;
}

#endif /* RADIATION_TRANSFER */
