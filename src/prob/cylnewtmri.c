
#include "copyright.h"
/*==============================================================================
 * FILE: newtmri.c
 *
 * Test of the evolution/saturation of the MRI in an unstratified/uniform
 * Newtonian disk
 *
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#define SEED 31337

static Real rho0,Amp,Beta, R0, Rm, Lb,rhomin, Flux, Hbc, Mbc;
Real Mollifier(Real Center, Real Scl, Real x);
void ScaleToBeta(GridS *pG, Real beta);
void disk_ir(GridS *pG);
void disk_or(GridS *pG);
static Real Mrp(const GridS *pG, const int i, const int j, const int k);
static Real Trp(const GridS *pG, const int i, const int j, const int k);
static Real Pb(const GridS *pG, const int i, const int j, const int k);
static Real Vaz(const GridS *pG, const int i, const int j, const int k);
static Real Br(const GridS *pG, const int i, const int j, const int k);
static Real Bp(const GridS *pG, const int i, const int j, const int k);
static Real Bz(const GridS *pG, const int i, const int j, const int k);
static Real Vr(const GridS *pG, const int i, const int j, const int k);
static Real Vp(const GridS *pG, const int i, const int j, const int k);
static Real Vz(const GridS *pG, const int i, const int j, const int k);
static Real MdotR1(const GridS *pG, const int i, const int j, const int k);
static Real MdotR2(const GridS *pG, const int i, const int j, const int k);
static Real MdotR3(const GridS *pG, const int i, const int j, const int k);
static Real MdotR4(const GridS *pG, const int i, const int j, const int k);



static Real grav_pot(const Real x1, const Real x2, const Real x3) {

	return (-1.0/x1);
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
	Real g;
	g = (1/SQR(x1));

	return g;

}

static Real Omega(const Real R) {
	Real Arg;
	Arg = pow(R,1.5);
	return (1.0/Arg);
}
static Real Shear(const Real R) {
	return 1.5;
}

static Real vphi(const Real x1, const Real x2, const Real x3) {
	return x1*Omega(x1);
}


static Real BzZero(const Real R, const Real phi, const Real z) {
	Real Arg, bz, H0, m;


	H0 = sqrt(2.0)*Iso_csound/Omega(R0); 
	Arg = 2.0*PI*(R-R0)/H0;
	m = (int)(R0/H0);
	bz = Omega(R)*cos(Arg)*sin(m*phi)*Mollifier(Rm,Lb,R);
	return bz;
}
static Real BzNet(const Real R, const Real phi, const Real z) {
	
	return 0.0;
}


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3,myid=0;
  Real x1,x2,x3,y1, r;
	GridS *pG = pDomain->Grid;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;


  rho0        = par_getd("problem", "rho0");
  Amp         = par_getd("problem", "Amp");
  Beta        = par_getd("problem", "Beta");
  R0          = par_getd("problem", "R0");
	Rm          = par_getd("problem", "Rm");
  Lb          = par_getd("problem","Lb");
  rhomin      = par_getd("problem","rhomin");
  Flux      = par_getd("problem","Flux");
  Hbc       = par_getd("problem","Hbc");
  Mbc       = par_getd("problem","Mbc");

  srand(SEED + myID_Comm_world);
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        r = 2.0*rand()/((double)RAND_MAX) - 1.0;
        pG->U[k][j][i].d = rho0;
#ifndef FARGO
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*avg1d(vphi,pG,i,j,k)*(1.0+Mollifier(Rm,Lb,x1)*Amp*r);
#else
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*avg1d(vphi,pG,i,j,k)*(Mollifier(Rm,Lb,x1)*Amp*r);
#endif
				r = 2.0*rand()/((double)RAND_MAX)-1.0;
        pG->U[k][j][i].M1 = 0.0;
        pG->U[k][j][i].M3 = rho0*Amp*r;
#ifdef MHD
        pG->U[k][j][i].B1c = 0.0;
        pG->U[k][j][i].B2c = 0.0;
        if (Flux == 0) {
          pG->U[k][j][i].B3c = avg2d(BzZero,pG,i,j,k);
        } else {
          pG->U[k][j][i].B3c = avg2d(BzNet,pG,i,j,k);
        }

        pG->B1i[k][j][i] = 0.0;
        pG->B2i[k][j][i] = 0.0;
        pG->B3i[k][j][i] = pG->U[k][j][i].B3c;

#endif

      }
    }
  }
#ifdef MHD
  ScaleToBeta(pG,Beta);
#endif /* MHD */

  StaticGravPot = grav_pot;
  x1GravAcc = grav_acc;
  bvals_mhd_fun(pDomain,left_x1,disk_ir);
  bvals_mhd_fun(pDomain,right_x1,disk_or);
#ifdef FARGO
  OrbitalProfile = Omega;
  ShearProfile = Shear;
#endif
	// Enroll history functions
	//
#ifdef MHD
	dump_history_enroll(Br, "<Br>");
	dump_history_enroll(Bp, "<Bp>");
	dump_history_enroll(Bz, "<Bz>");
	dump_history_enroll(Mrp, "<Mrp>");
	dump_history_enroll(Trp, "<Trp>");
	dump_history_enroll(MdotR1, "<MdotR1>");
	dump_history_enroll(MdotR2, "<MdotR2>");
	dump_history_enroll(MdotR3, "<MdotR3>");
	dump_history_enroll(MdotR4, "<MdotR4>");

#endif
  return;
}
/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/
static Real MdotR1(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=1.0, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;

}

static Real MdotR2(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=2.0, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;
}

static Real MdotR3(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=3.0, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;
}
static Real MdotR4(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=4.0, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;
}
static Real Pb(const GridS *pG, const int i, const int j, const int k)
{
#ifdef MHD
	return ( 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c)) );
#else
	return 0.0;
#endif
}

static Real Mrp(const GridS *pG, const int i, const int j, const int k)
{
#ifdef MHD
	return ( -1.0*pG->U[k][j][i].B1c*pG->U[k][j][i].B2c );  
#else
	return 0.0;
#endif
}

static Real Trp(const GridS *pG, const int i, const int j, const int k)
{
	// Note, this assumes Fargo to calculate perturbation Vp
	return ((pG->U[k][j][i].M1*pG->U[k][j][i].M2)/(pG->U[k][j][i].d));
}
static Real Vaz(const GridS *pG, const int i, const int j, const int k)
{
#ifdef MHD
	return (fabs(pG->U[k][j][i].B3c)/sqrt(pG->U[k][j][i].d));
#else
	return 0.0;
#endif
}
static Real Br(const GridS *pG, const int i, const int j, const int k) {
#ifdef MHD
	return pG->U[k][j][i].B1c;
#else
	return 0.0;
#endif

}
static Real Bp(const GridS *pG, const int i, const int j, const int k) {
#ifdef MHD
	return pG->U[k][j][i].B2c;
#else
	return 0.0;
#endif

}
static Real Bz(const GridS *pG, const int i, const int j, const int k) {
#ifdef MHD
	return pG->U[k][j][i].B3c;
#else
	return 0.0;
#endif

}
static Real Vr(const GridS *pG, const int i, const int j, const int k) {
	return (pG->U[k][j][i].M1/pG->U[k][j][i].d);

}
static Real Vp(const GridS *pG, const int i, const int j, const int k) {
	return (pG->U[k][j][i].M2/pG->U[k][j][i].d);

}
static Real Vz(const GridS *pG, const int i, const int j, const int k) {
	return (pG->U[k][j][i].M3/pG->U[k][j][i].d);

}
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
	if (strcmp(expr,"Pb")==0) return Pb;
	else if(strcmp(expr,"Mrp")==0) return Mrp;
	else if(strcmp(expr,"Trp")==0) return Trp;
	else if(strcmp(expr,"Vaz")==0) return Vaz;
	else if(strcmp(expr,"Br")==0) return Br;
	else if(strcmp(expr,"Bp")==0) return Bp;
	else if(strcmp(expr,"Bz")==0) return Bz;
	else if(strcmp(expr,"Vr")==0) return Vr;
	else if(strcmp(expr,"Vp")==0) return Vp;
	else if(strcmp(expr,"Vz")==0) return Vz;

  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
}

void Userwork_after_loop(MeshS *pM)
{
}
// Private functions
//
Real Mollifier(Real Center, Real Scl, Real z) {
        Real Arg, a, M;

        a = fabs( (z-Center)/Scl );

        if (a + TINY_NUMBER < 1.0) {
                Arg = 1.0 - SQR(a);
                M = exp(-1.0/Arg);
        } else {
                M = 0.0;
        }
        return M;
}

void ScaleToBeta(GridS *pG, Real beta) {
#ifdef MHD
  int i,j,k;
  int is,ie,js,je,ks,ke,nx1,nx2,nx3;
  int il,iu,jl,ju,kl,ku;
  Real Pgas = 0.0, Pb = 0.0, TotPgas = 0.0, TotPb = 0.0, scl = 0.0, CellVol;

  is = pG->is; ie = pG->ie;
  js = pG->js; je = pG->je;
  ks = pG->ks; ke = pG->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;
  il = is - nghost*(ie > is);
  jl = js - nghost*(je > js);
  kl = ks - nghost*(ke > ks);
  iu = ie + nghost*(ie > is);
  ju = je + nghost*(je > js);
  ku = ke + nghost*(ke > ks);

  /* Calculate total gas/magnetic pressure on local tile */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
				CellVol = (pG->r[i])*(pG->dx1)*(pG->dx2)*(pG->dx3);
#ifndef ISOTHERMAL
        Pgas += (Gamma-1)*pG->U[k][j][i].E*CellVol;
#else
        Pgas += Iso_csound2*pG->U[k][j][i].d*CellVol;
#endif
        Pb += 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c)
              + SQR(pG->U[k][j][i].B3c))*CellVol;
      }
    }
  }

#ifdef MPI_PARALLEL
  MPI_Reduce(&Pgas, &TotPgas, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce(&Pb, &TotPb, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotPgas, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&TotPb, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#else
  TotPgas = Pgas;
  TotPb = Pb;
#endif //PARALLEL
  //calculate and use scaling factor so that the correct beta is ensured
  scl = sqrt(TotPgas/(TotPb*beta));
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i].B1c *= scl;
        pG->U[k][j][i].B2c *= scl;
        pG->U[k][j][i].B3c *= scl;
        pG->B1i[k][j][i]   *= scl;
        pG->B2i[k][j][i]   *= scl;
        pG->B3i[k][j][i]   *= scl;
      }
    }
  }
#endif /* MHD */
}

void disk_ir(GridS *pGrid) {
  int is = pGrid->is;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
	Real Vkep,R,p,z, RBr, Lper;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	cc_pos(pGrid,is,j,k,&R,&p,&z);
#ifdef MHD
    	RBr = (R-0.5*pGrid->dx1)*pGrid->B1i[k][j][is];
#endif
    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = R*pGrid->U[k][j][is].M2;
#else
    	Vkep = sqrt(R*(*x1GravAcc)(R,p,z));
    	Lper = R*pGrid->U[k][j][is].M2 - R*pGrid->U[k][j][is].d*Vkep;
#endif
    	Lper = MIN(Lper,0.0);

      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
				// Enforce diode
				pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);
				// Calculate Keplerian velocity
				cc_pos(pGrid,is-i,j,k,&R,&p,&z);
				Vkep = sqrt(R*(*x1GravAcc)(R,p,z));
#ifdef FARGO
				if (Hbc == 1) {
					pGrid->U[k][j][is-i].M2 = 0.0;
				} else {
					pGrid->U[k][j][is-i].M2 = Lper/R;
				}
#else
				if (Hbc == 1) {
					pGrid->U[k][j][is-i].M2 = pGrid->U[k][j][is].d*Vkep;
				} else {
					pGrid->U[k][j][is-i].M2 = Lper/R + pGrid->U[k][j][is-i].d*Vkep;
				}
#endif
#ifdef MHD
				if (Mbc == 0) {
					pGrid->U[k][j][is-i].B1c = RBr/R;
					pGrid->U[k][j][is-i].B2c = 0.0;
					pGrid->U[k][j][is-i].B3c = 0.0;
				}
#endif
      }
    }
  }

#ifdef MHD
/* B1i is not set at i=is-nghost */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	cc_pos(pGrid,is,j,k,&R,&p,&z);
    	RBr = (R-0.5*pGrid->dx1)*pGrid->B1i[k][j][is];
      for (i=1; i<=nghost-1; i++) {
      	cc_pos(pGrid,is,j,k,&R,&p,&z);
      	if (Mbc == 0) {
      		pGrid->B1i[k][j][is-i] = RBr/(R-0.5*pGrid->dx1);
      	} else {
      		pGrid->B1i[k][j][is-i] = pGrid->B1i[k][j][is];
      	}
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
      	if (Mbc == 1) {
      		pGrid->B2i[k][j][is-i] = pGrid->B2i[k][j][is];
      	} else {
      		pGrid->B2i[k][j][is-i] = 0.0;
      	}
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      	if (Mbc == 1) {
      		pGrid->B3i[k][j][is-i] = pGrid->B3i[k][j][is];
      	} else {
      		pGrid->B3i[k][j][is-i] = 0.0;
      	}

      }
    }
  }
#endif /* MHD */

  return;


}

void disk_or(GridS *pGrid) {

	int ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  int i,j,k;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
	Real Vkep, R,p,z, RBr, Lper;

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
#ifdef MHD
    	RBr = pGrid->ri[ie+1]*pGrid->B1i[k][j][ie+1];
#endif
    	// Calculate angular momentum in rotating frame
#ifdef FARGO
    	Lper = pGrid->r[ie]*pGrid->U[k][j][ie].M2;
#else
			cc_pos(pGrid,ie,j,k,&R,&p,&z);
			Vkep = sqrt(R*(*x1GravAcc)(R,p,z));
			Lper = R*pGrid->U[k][j][ie].M2 - x1p*pGrid->U[k][j][ie]*Vkep;
#endif
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
				// Enforce diode
				pGrid->U[k][j][ie+i].M1 = MAX(pGrid->U[k][j][ie+i].M1,0.0);
				cc_pos(pGrid,ie+i,j,k,&R,&p,&z);
				Vkep = sqrt(R*(*x1GravAcc)(R,p,z));
#ifdef FARGO
				if (Hbc == 1) {
					pGrid->U[k][j][ie+i].M2 = 0.0;
				} else {
					pGrid->U[k][j][ie+i].M2 = Lper/R;
				}
#else
				if (Hbc == 1) {
					pGrid->U[k][j][ie+i].M2 = pGrid->U[k][j][ie].d*Vkep;
				} else {
					pGrid->U[k][j][ie+i].M2 = Lper/R + pGrid->U[k][j][ie+i].d*Vkep;
				}
#endif
#ifdef MHD
				if (Mbc == 0) {
					pGrid->U[k][j][ie+i].B1c = RBr/R;
					pGrid->U[k][j][ie+i].B2c = 0.0;
					pGrid->U[k][j][ie+i].B3c = 0.0;
				}
#endif
	
      }
    }
  }

#ifdef MHD
/* i=ie+1 is not a boundary condition for the interface field B1i */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
    	RBr = pGrid->ri[ie+1]*pGrid->B1i[k][j][ie+1];
      for (i=2; i<=nghost; i++) {
				cc_pos(pGrid,ie+i,j,k,&R,&p,&z);
      	if (Mbc == 0) {
      		pGrid->B1i[k][j][ie+i] = RBr/R;
      	} else {
      		pGrid->B1i[k][j][ie+i] = pGrid->B1i[k][j][ie];
      	}
      }
    }
  }

  if (pGrid->Nx[1] > 1) ju=je+1; else ju=je;
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=ju; j++) {
      for (i=1; i<=nghost; i++) {
      	if (Mbc == 0) {
      		pGrid->B2i[k][j][ie+i] = 0.0;
      	} else {
      		pGrid->B2i[k][j][ie+i] = pGrid->B2i[k][j][ie];
      	}
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++) {
      	if (Mbc == 0) {
      		pGrid->B3i[k][j][ie+i] = 0.0;
      	} else {
      		pGrid->B3i[k][j][ie+i] = pGrid->B3i[k][j][ie];
      	}
      }
    }
  }
#endif /* MHD */

  return;

}


