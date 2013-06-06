
#include "copyright.h"
/*============================================================================*/
/*! \file cylcvmri.c
 * \brief Test of the evolution/saturation of the MRI in an unstratified/uniform
 * Newtonian disk in Cataclysmic Variable (CV) system.
 * Features: Rotating frame; Roche lobe overflow through L1; Cylindrical Coordinate.
 * Last modified: by Wendy Ju, June 04, 2013
 * Last added: sub regions(r<2) and sub2 regions(r<3) 
 * Filename in previous versions: cylnewtmri-mhd.c
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#define SEED 31337
#define Pi   3.14159

static Real rho0, d_floor, Amp,Beta, R0, Field, Hbc, Mbc, Mc, IC;
static Real q_m, R_comp; /* q_m: M_dwarf/M_companion, R_comp: distance between companion star and white dwarf, by Wendy*/
#ifndef ISOTHERMAL
static Real pgas, Iso_csound, Tmin, presmin;
#endif
#ifdef L1_INFLOW
static Real Vr_l1;
#endif
Real Mollifier(Real Center, Real Scl, Real x);
void ScaleToBeta(GridS *pG, Real beta);
void disk_ir(GridS *pG);
void disk_or(GridS *pG);
static Real ChiMag(Real R);
static Real ChiSub(Real R);
static Real vphi(const Real x1, const Real x2, const Real x3);
static Real Br2(const GridS *pG, const int i, const int j, const int k);
static Real Mrp(const GridS *pG, const int i, const int j, const int k);
static Real Trp(const GridS *pG, const int i, const int j, const int k);
static Real rdVp(const GridS *pG, const int i, const int j, const int k);
static Real Pb(const GridS *pG, const int i, const int j, const int k);
static Real Pg(const GridS *pG, const int i, const int j, const int k);
static Real Vaz(const GridS *pG, const int i, const int j, const int k);
static Real Va(const GridS *pG, const int i, const int j, const int k);
static Real Br(const GridS *pG, const int i, const int j, const int k);
static Real Bp(const GridS *pG, const int i, const int j, const int k);
static Real Bz(const GridS *pG, const int i, const int j, const int k);
static Real Vr(const GridS *pG, const int i, const int j, const int k);
static Real Vp(const GridS *pG, const int i, const int j, const int k);
static Real Vz(const GridS *pG, const int i, const int j, const int k);
static Real L_AM(const GridS *pG, const int i, const int j, const int k);
static Real MR1(const GridS *pG, const int i, const int j, const int k);
static Real MR2(const GridS *pG, const int i, const int j, const int k);
static Real MR3(const GridS *pG, const int i, const int j, const int k);
static Real Mdot(const GridS *pG, const int i, const int j, const int k);/*to produce Mdot(r)*/
static Real MdotR1(const GridS *pG, const int i, const int j, const int k);/*to produce Mdot(r=0.04,t)*/
static Real MdotR2(const GridS *pG, const int i, const int j, const int k);/*to produce Mdot(r=0.1,t)*/ 
static Real MdotR3(const GridS *pG, const int i, const int j, const int k);/*to produce Mdot(r=0.3,t)*/ 
static Real MdotR4(const GridS *pG, const int i, const int j, const int k);/*to produce Mdot(r=0.5,t)*/ 
static Real Msub(const GridS *pG, const int i, const int j, const int k);
static Real Msub2(const GridS *pG, const int i, const int j, const int k);
static Real Trpsub(const GridS *pG, const int i, const int j, const int k);
static Real Trpsub2(const GridS *pG, const int i, const int j, const int k);
static Real Mrpsub(const GridS *pG, const int i, const int j, const int k);
static Real Mrpsub2(const GridS *pG, const int i, const int j, const int k);
static Real Bzsub(const GridS *pG, const int i, const int j, const int k);
static Real Bpsub(const GridS *pG, const int i, const int j, const int k);
static Real Pbsub(const GridS *pG, const int i, const int j, const int k);
static Real Pbsub2(const GridS *pG, const int i, const int j, const int k);
static Real Pgsub(const GridS *pG, const int i, const int j, const int k);
static Real Pgsub2(const GridS *pG, const int i, const int j, const int k);
static Real Temp(const GridS *pG, const int i, const int j, const int k);
static void dump_vtksub(MeshS *pM, OutputS *pOut);
void out_ktab(MeshS *pM, OutputS *pOut);
static void out_jtab(MeshS *pM, OutputS *pOut); /*tabular file output by Wendy Ju*/

static Real grav_pot(const Real x1, const Real x2, const Real x3) {
  #ifdef ROTATING_FRAME
  Real r2, Pot;
  r2 = sqrt(SQR(x1) + SQR(R_comp) - 2.0*x1*R_comp*cos(x2));
  Pot = (-1.0/x1 - q_m/r2)/(1.0+q_m) - 0.5*SQR(Omega_0*x1); /* Scaled unit: GM0 = G(M1+M2) =1 */
  return Pot;
  #else
  return (-1.0/x1)/(1.0+q_m); /* Scaled unit: M0 = M1+M2*/
  #endif
}

static Real grav_acc(const Real x1, const Real x2, const Real x3) {
  Real g, r2;
  #ifdef ROTATING_FRAME
  r2 = sqrt(SQR(x1) + SQR(R_comp) - 2.0*x1*R_comp*cos(x2));
  g = (1.0/SQR(x1) + q_m*(x1-R_comp*cos(x2))/(r2*SQR(r2)))/(1.0+q_m) - SQR(Omega_0)*x1; /* Scaled unit: M0 = M1+M2*/
  #else 
  g = (1.0/SQR(x1)/(1.0+q_m)); /* Scaled unit: M0 = M1+M2*/
  #endif
  return g;
}

static Real Omega(const Real R) {
	Real Arg;
	Arg = pow(R,1.5)*sqrt(1.0+q_m); /* Scaled unit: M0 = M1+M2*/
	return (1.0/Arg);
}
static Real Shear(const Real R) {
	return 1.5;
}

/*Kepler velocity*/
static Real vphi(const Real x1, const Real x2, const Real x3) {
	#ifdef ROTATING_FRAME
	return x1*(Omega(x1)-Omega_0);
	#else
	return x1*Omega(x1);
	#endif
}


static Real BzZero(const Real R, const Real phi, const Real z) {
	Real Arg, bz, H0, n, Rs, I;

	H0 = sqrt(2.0)*Iso_csound/Omega(R0);
	Arg = 2.0*PI*(R-R0)/H0;
	n = floor( (R-1.2)/H0 );
	Rs = 1.3 + n*H0;
	I = ChiMag(R);
	bz = I*sin(Arg)*(Rs/R)*Omega(Rs);
	return bz;
}

static Real BzNet(const Real R, const Real phi, const Real z) {
	Real bz, I;

	I = ChiMag(R);
	bz = I*Omega(R);
	return bz;
}

static Real BpNet(const Real R, const Real phi, const Real z) {
	Real Bp, I;
	I = ChiMag(R);
	Bp = sqrt(rho0/R)*I/Mc;
	return Bp;
}

static Real ChiMag(Real R) {
	Real I;
/*	if ((R >= 1.2) && (R <= 3.8)) {
		I = 1.0;
	} else {
		I = 0.0;
	}
*/	I = 1.0;
	return I;
}
static Real ChiSub(Real R) {
	Real I;
	if (R <= 0.2) {
		I = 1.0;
	} else {
		I = 0.0;
	}
	return I;
}
static Real ChiSub2(Real R) {
        Real I;
        if (R <= 0.3) {
                I = 1.0;
        } else {
                I = 0.0;
        }
        return I;
}

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  int i,j,k;
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  int nx1,nx2,nx3,myid=0;
  Real x1M;
  Real x1,x2,x3,y1, r;
	GridS *pG = pDomain->Grid;

  is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
  js = pG->js;  je = pG->je;  nx2 = je-js+1;
  ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

  il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);  nx1 = iu-il+1;
  jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);  nx2 = ju-jl+1;
  kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);  nx3 = ku-kl+1;

#ifndef CYLINDRICAL
ath_error("[problem]: This problem generator only applies to cylindrical coordinnates.\n");
#endif

  rho0        = par_getd_def("problem", "rho0", 100.0);
  Amp         = par_getd_def("problem", "Amp", 1.0e-2);
  Beta        = par_getd_def("problem", "Beta",100.0);
  d_floor     = par_getd_def("problem","d_floor",1.0e-3);
  d_MIN	      = 0.1*d_floor;
  presmin     = par_getd_def("problem","presmin",1.0e-8);
  Field      = par_getd_def("problem","Field",0);
  Hbc       = par_getd_def("problem","Hbc",1);
  Mbc       = par_getd_def("problem","Mbc",1);
  Mc        = par_getd_def("problem","Mc",20.0);
  IC        = par_getd_def("problem","IC",1); 		/*flag for Initial Conditions*/
#ifndef ISOTHERMAL
  pgas      = par_getd_def("problem","pgas", 1.0e-4);
  Iso_csound = par_getd_def("problem","iso_csound", 0.01);
  Tmin	    = par_getd_def("problem", "Tmin", 1.0e-4);
#endif
#ifdef ROTATING_FRAME
  q_m       = par_getd_def("problem","q_m", 1.0);
  R_comp    = 1.0; 
  Omega_0 = sqrt(1.0/pow(R_comp,3.0)); /* Scaled unit: M0 = M1+M2*/
  Rc = q_m/(1.0+q_m)*R_comp;
#endif
#ifdef L1_INFLOW
  Vr_l1 =  par_getd_def("problem","Vr_l1", 0.01);
#endif

  srand(SEED + myID_Comm_world);
  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        x1 = x1vc(pG,i);
        r = 2.0*rand()/((double)RAND_MAX) - 1.0;

	/*Initial Condition 1: for no inflow*/
	if(IC == 1.0) {
	pG->U[k][j][i].d = rho0;
	pG->U[k][j][i].M1 = 0.0;
#ifdef FARGO
        pG->U[k][j][i].M2 = 0.0;
#else
	pG->U[k][j][i].M2 = pG->U[k][j][i].d*vphi(x1,x2,x3);
#endif
	pG->U[k][j][i].M2 += pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
                                r = 2.0*rand()/((double)RAND_MAX)-1.0;
	pG->U[k][j][i].M3 = pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
#ifndef ISOTHERMAL
	pG->U[k][j][i].E  = pgas/(Gamma*Gamma_1)
                            + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#endif
	} /*IC 1*/

        /*Initial Condition 2: for L1-inflow*/
        if(IC == 2.0) {
        x1M = pG->MaxX[0];
        if(x1<=0.2*x1M) pG->U[k][j][i].d = 0.2*rho0;
        else pG->U[k][j][i].d = d_floor;

        pG->U[k][j][i].M1 = 0.0;
#ifdef FARGO
        pG->U[k][j][i].M2 = 0.0;
#else
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*vphi(x1,x2,x3);
#endif
        pG->U[k][j][i].M2 += pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
                                r = 2.0*rand()/((double)RAND_MAX)-1.0;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
#ifndef ISOTHERMAL
        if(x1<=0.2*x1M) pG->U[k][j][i].E  = 0.2*pgas/(Gamma*Gamma_1)
                                        + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
        else           pG->U[k][j][i].E  = 0.01*pgas/(Gamma*Gamma_1)
                                        + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#endif
        } /*IC 2*/

        /*Initial Condition 3:*/
        if(IC == 3) {
	R0 = pG->MinX[0];
        pG->U[k][j][i].d = 0.2*rho0*pow((x1/R0), -1.5);
        pG->U[k][j][i].M1 = 0.0;
#ifdef FARGO
        pG->U[k][j][i].M2 = 0.0;
#else
        pG->U[k][j][i].M2 = pG->U[k][j][i].d*vphi(x1,x2,x3);
#endif
        pG->U[k][j][i].M2 += pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
                                r = 2.0*rand()/((double)RAND_MAX)-1.0;
        pG->U[k][j][i].M3 = pG->U[k][j][i].d*Amp*r*Iso_csound*ChiMag(x1);
#ifndef ISOTHERMAL
	pG->U[k][j][i].E  = 2.0*pgas/(Gamma*Gamma_1)*pow((x1/R0), -1.5)
		            + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
#endif
        } /*IC 3*/

	if((IC!=1)&&(IC!=2)&&(IC!=3)) ath_error("The initial condition flag is not defined."); 

#ifdef MHD
        pG->U[k][j][i].B1c = 0.0;

        if (Field == 2) {
        	pG->U[k][j][i].B2c = BpNet(x1,x2,x3);
        } else {
          pG->U[k][j][i].B2c = 0.0;
        }
        if (Field == 0) {
          pG->U[k][j][i].B3c = BzZero(x1,x2,x3);
        } else if (Field == 1) {
          pG->U[k][j][i].B3c = BzNet(x1,x2,x3);
        } else {
        	pG->U[k][j][i].B3c = 0.0;
        }

        pG->B1i[k][j][i] = 0.0;
        pG->B2i[k][j][i] = pG->U[k][j][i].B2c;
        pG->B3i[k][j][i] = pG->U[k][j][i].B3c;
	
	/* Initial Condition for Initial Vacuum Disk */
	if (Field == 3){
		pG->U[k][j][i].B1c = 0.0;
		pG->U[k][j][i].B2c = 0.0;
		pG->U[k][j][i].B3c = 0.0;
		pG->B1i[k][j][i] = 0.0;
        	pG->B2i[k][j][i] = 0.0;
        	pG->B3i[k][j][i] = 0.0;
	}

#ifndef ISOTHERMAL
        pG->U[k][j][i].E  += 0.5* (SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c)) ;
#endif
#endif /* MHD */

      }
    }
  }
#ifdef MHD
  if (Field != 2 && Field != 3) {
  	ScaleToBeta(pG,Beta);
  }
#endif /* MHD */

  StaticGravPot = grav_pot;
  bvals_mhd_fun(pDomain,left_x1,disk_ir);
  bvals_mhd_fun(pDomain,right_x1,disk_or);
#ifdef FARGO
  OrbitalProfile = Omega;
  ShearProfile = Shear;
#endif
        dump_history_enroll(MR1, "<MR1>");
        dump_history_enroll(MR2, "<MR2>");
        dump_history_enroll(MR3, "<MR3>");
        dump_history_enroll(MdotR1, "<MdotR1>");
        dump_history_enroll(MdotR2, "<MdotR2>");
        dump_history_enroll(MdotR3, "<MdotR3>");
        dump_history_enroll(MdotR4, "<MdotR4>");
	dump_history_enroll(L_AM, "<L_AngularM>");
        dump_history_enroll(Trp, "<Trp>");
        dump_history_enroll(Trpsub, "<Trpsub>");
        dump_history_enroll(Trpsub2, "<Trpsub2>");
#ifdef MHD
	dump_history_enroll(Br, "<Br>");
	dump_history_enroll(Bp, "<Bp>");
	dump_history_enroll(Bz, "<Bz>");
	dump_history_enroll(Br2, "<Br2>");
	dump_history_enroll(Mrp, "<Mrp>");
        dump_history_enroll(Mrpsub, "<Mrpsub>");
        dump_history_enroll(Mrpsub2, "<Mrpsub2>");
	dump_history_enroll(Pb, "<Pb>");
        dump_history_enroll(Pbsub, "<Pbsub>");
        dump_history_enroll(Pbsub2, "<Pbsub2>");
#endif
        dump_history_enroll(Pg, "<Pg>");
        dump_history_enroll(Pgsub, "<Pgsub>");
        dump_history_enroll(Pgsub2, "<Pgsub2>");
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
/*Mass enclosed within R0=0.1*/
static Real MR1(const GridS *pG, const int i, const int j, const int k) {
        Real R0=0.1, R, phi, z, Mass;

        cc_pos(pG,i,j,k,&R,&phi,&z);
        if (R<R0) {
                Mass = pG->U[k][j][i].d ; 
        } else {
                Mass = 0.0;
        }
        return Mass;
}
/*Mass enclosed within R0=0.5*/
static Real MR2(const GridS *pG, const int i, const int j, const int k) {
        Real R0=0.3, R, phi, z, Mass;

        cc_pos(pG,i,j,k,&R,&phi,&z);
        if (R<R0) {
                Mass = pG->U[k][j][i].d ;
        } else {
                Mass = 0.0;
        }
        return Mass;
}
/*Mass enclosed within R0=0.75*/
static Real MR3(const GridS *pG, const int i, const int j, const int k) {
        Real R0=0.6, R, phi, z, Mass;

        cc_pos(pG,i,j,k,&R,&phi,&z); 
        if (R<R0) {
                Mass = pG->U[k][j][i].d ;
        } else {
                Mass = 0.0;
        }
        return Mass;
}

/*Mdot used to calculated phi-integrated value as a function of r*/
static Real Mdot(const GridS *pG, const int i, const int j, const int k) {
        Real dphi=pG->dx2, dz=pG->dx3, R, phi, z, Md;

        cc_pos(pG,i,j,k,&R,&phi,&z);
	Md = -1.0*pG->U[k][j][i].M1 * R;
        if (dphi>0) Md *= dphi;
	if (dz>0) Md *= dz;
	Md *= (pG->je-pG->js+1)*(pG->ke-pG->ks+1);
        return Md;
}

static Real MdotR1(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=0.07, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;

}

static Real MdotR2(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=0.2, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;
}

static Real MdotR3(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=0.3, R, phi, z, Md;

	cc_pos(pG,i,j,k,&R,&phi,&z);
	if ((R0 > (R-dR)) && (R0 < (R+dR))) {
		Md = -1.0*pG->U[k][j][i].M1/dR;
	} else {
		Md = 0.0;
	}
	return Md;
}
static Real MdotR4(const GridS *pG, const int i, const int j, const int k) {
	Real dR=pG->dx1, R0=0.62, R, phi, z, Md;

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
static Real Pg(const GridS *pG, const int i, const int j, const int k)
{
#ifndef ISOTHERMAL
    Real Pb, Pg, Ek;
        Ek = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#ifdef MHD
        Pb = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
	Pg = (pG->U[k][j][i].E - Pb - Ek) * (Gamma-1);
#else
	Pg = (pG->U[k][j][i].E - Ek) * (Gamma-1);
#endif
#else
    Real Pg;
	Pg = Iso_csound2*pG->U[k][j][i].d; 
#endif
	return Pg;
}

static Real Br2(const GridS *pG, const int i, const int j, const int k)
{
#ifdef MHD
        return ( pG->U[k][j][i].B1c * pG->U[k][j][i].B1c );
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
	Real R,p,z;
	cc_pos(pG,i,j,k,&R,&p,&z);
	return ((pG->U[k][j][i].M1*(pG->U[k][j][i].M2-vphi(R,p,z)*pG->U[k][j][i].d))/(pG->U[k][j][i].d));
}
/* rho * (Vphi - Vkep)*/
static Real rdVp(const GridS *pG, const int i, const int j, const int k)
{
        Real R,p,z;
        cc_pos(pG,i,j,k,&R,&p,&z);
        return (pG->U[k][j][i].M2-vphi(R,p,z)*pG->U[k][j][i].d);
}

static Real Va(const GridS *pG, const int i, const int j, const int k)
{
#ifdef MHD
        return (sqrt((SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c))/pG->U[k][j][i].d));
#else
        return 0.0;
#endif
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
	Real R=x1vc(pG, i);
#ifdef ROTATING_FRAME
	return (pG->U[k][j][i].M2/pG->U[k][j][i].d + Omega_0*R);
#else
	return (pG->U[k][j][i].M2/pG->U[k][j][i].d);
#endif

}
static Real Vz(const GridS *pG, const int i, const int j, const int k) {
	return (pG->U[k][j][i].M3/pG->U[k][j][i].d);

}

static Real L_AM(const GridS *pG, const int i, const int j, const int k) {
	Real L, R;
	
	R = x1vc(pG, i);
#ifdef ROTATING_FRAME
	L = (pG->U[k][j][i].M2 + pG->U[k][j][i].d*R*Omega_0) * R;
#else
	L = pG->U[k][j][i].M2 * R;
#endif
	return(L);
}


static Real Msub(const GridS *pG, const int i, const int j, const int k) {
	Real I, R;

	R = x1vc(pG,i);
	I = ChiSub(R);
	return (pG->U[k][j][i].d*I);
}
static Real Msub2(const GridS *pG, const int i, const int j, const int k) {
        Real I, R;

        R = x1vc(pG,i);
        I = ChiSub2(R);
        return (pG->U[k][j][i].d*I);
}

static Real Trpsub(const GridS *pG, const int i, const int j, const int k){
	Real I, R0;
        R0 = x1vc(pG,i);
        I = ChiSub(R0);
#ifdef ROTATING_FRAME
        Real R,p,z;
        cc_pos(pG,i,j,k,&R,&p,&z);
        return ((pG->U[k][j][i].M1*(pG->U[k][j][i].M2-vphi(R,p,z)*pG->U[k][j][i].d))/(pG->U[k][j][i].d) * I);
#else
        return ((pG->U[k][j][i].M1*pG->U[k][j][i].M2)/(pG->U[k][j][i].d) * I);
#endif
}
static Real Trpsub2(const GridS *pG, const int i, const int j, const int k){
        Real I, R0;
        R0 = x1vc(pG,i);
        I = ChiSub2(R0);
#ifdef ROTATING_FRAME
        Real R,p,z;
        cc_pos(pG,i,j,k,&R,&p,&z);
        return ((pG->U[k][j][i].M1*(pG->U[k][j][i].M2-vphi(R,p,z)*pG->U[k][j][i].d))/(pG->U[k][j][i].d) * I);
#else
        return ((pG->U[k][j][i].M1*pG->U[k][j][i].M2)/(pG->U[k][j][i].d) * I);
#endif
}

static Real Mrpsub(const GridS *pG, const int i, const int j, const int k) {
	Real I, R;

	R = x1vc(pG,i);
	I = ChiSub(R);
#ifdef MHD
	return (-1.0*pG->U[k][j][i].B1c*pG->U[k][j][i].B2c*I);
#else
	return 0.0;
#endif

}
static Real Mrpsub2(const GridS *pG, const int i, const int j, const int k) {
        Real I, R;

        R = x1vc(pG,i);
        I = ChiSub2(R);
#ifdef MHD
        return (-1.0*pG->U[k][j][i].B1c*pG->U[k][j][i].B2c*I);
#else
        return 0.0;
#endif

}

static Real Bpsub(const GridS *pG, const int i, const int j, const int k) {
	Real I, R;

	R = x1vc(pG,i);
	I = ChiSub(R);
#ifdef MHD
	return (pG->U[k][j][i].B2c*I);
#else
	return 0.0;
#endif

}

static Real Bzsub(const GridS *pG, const int i, const int j, const int k) {
	Real I, R;

	R = x1vc(pG,i);
	I = ChiSub(R);
#ifdef MHD
	return (pG->U[k][j][i].B3c*I);
#else
	return 0.0;
#endif

}

static Real Pbsub(const GridS *pG, const int i, const int j, const int k) {
	Real I, R;

	R = x1vc(pG,i);
	I = ChiSub(R);
#ifdef MHD
	return (0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c))*I);
#else
	return 0.0;
#endif
}
static Real Pbsub2(const GridS *pG, const int i, const int j, const int k) {
        Real I, R;

        R = x1vc(pG,i);
        I = ChiSub2(R);
#ifdef MHD
        return (0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c))*I);
#else
        return 0.0;
#endif
}

static Real Pgsub(const GridS *pG, const int i, const int j, const int k)
{
    Real Pg, I, R;
	R = x1vc(pG,i);
        I = ChiSub(R);
#ifndef ISOTHERMAL
    Real Pb, Ek;
        Ek = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#ifdef MHD
        Pb = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        Pg = (pG->U[k][j][i].E - Pb - Ek) * (Gamma-1);
#else
        Pg = (pG->U[k][j][i].E - Ek) * (Gamma-1);
#endif
#else
	Pg = Iso_csound2*pG->U[k][j][i].d;
#endif
	return (Pg*I);
}
static Real Pgsub2(const GridS *pG, const int i, const int j, const int k)
{
    Real Pg, I, R;
        R = x1vc(pG,i);
        I = ChiSub2(R);
#ifndef ISOTHERMAL
    Real Pb, Ek;
        Ek = 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
#ifdef MHD
        Pb = 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
        Pg = (pG->U[k][j][i].E - Pb - Ek) * (Gamma-1);
#else
        Pg = (pG->U[k][j][i].E - Ek) * (Gamma-1);
#endif
#else
        Pg = Iso_csound2*pG->U[k][j][i].d;
#endif
        return (Pg*I);
}


static Real Temp(const GridS *pG, const int i, const int j, const int k) {
        Real T;

#ifdef ISOTHERMAL
	T = Iso_csound2;
#else
	#ifdef MHD
        T = (pG->U[k][j][i].E - 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c)))/pG->U[k][j][i].d
            - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/SQR(pG->U[k][j][i].d);
	#else
	T = pG->U[k][j][i].E/pG->U[k][j][i].d
            - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/SQR(pG->U[k][j][i].d);
	#endif
#endif
	return T;
}


void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{

  rho0        = par_getd_def("problem", "rho0", 100.0);
  Amp         = par_getd_def("problem", "Amp", 1.0e-2);
  Beta        = par_getd_def("problem", "Beta",100.0);
  R0          = par_getd_def("problem", "R0",2.0);
  d_floor     = par_getd_def("problem","d_floor",1.0e-3);
  d_MIN       = 0.1*d_floor;
  presmin     = par_getd_def("problem","presmin",1.0e-8); 
  Field      = par_getd_def("problem","Field",0);
  Hbc       = par_getd_def("problem","Hbc",1);
  Mbc       = par_getd_def("problem","Mbc",1);
  Mc        = par_getd_def("problem","Mc",20.0);
  IC	     = par_getd_def("problem","IC",1);

#ifndef ISOTHERMAL
  pgas      = par_getd_def("problem","pgas", 1.0e-4);
  Iso_csound = par_getd_def("problem","iso_csound", 0.01);
  Tmin      = par_getd_def("problem", "Tmin", 1.0e-4);
#endif

#ifdef ROTATING_FRAME
  q_m       = par_getd_def("problem","q_m", 1.0);
  R_comp    = 1.0; /*Length Unit*/
/*  Omega_0 = sqrt((1.0+q_m)/pow(R_comp,3.0));
*/  Omega_0 = sqrt(1.0/pow(R_comp,3.0));
  Rc = q_m/(1.0+q_m)*R_comp;
#endif

#ifdef L1_INFLOW
  Vr_l1 =  par_getd_def("problem","Vr_l1", 0.01);
#endif

  StaticGravPot = grav_pot;
  bvals_mhd_fun(&(pM->Domain[0][0]),left_x1,disk_ir);
  bvals_mhd_fun(&(pM->Domain[0][0]),right_x1,disk_or);
#ifdef FARGO
  OrbitalProfile = Omega;
  ShearProfile = Shear;
#endif
        dump_history_enroll(MR1, "<MR1>");
        dump_history_enroll(MR2, "<MR2>");
        dump_history_enroll(MR3, "<MR3>");
        dump_history_enroll(MdotR1, "<MdotR1>");
        dump_history_enroll(MdotR2, "<MdotR2>");
        dump_history_enroll(MdotR3, "<MdotR3>");
        dump_history_enroll(MdotR4, "<MdotR4>");
        dump_history_enroll(L_AM, "<L_AngularM>");
        dump_history_enroll(Trp, "<Trp>");
        dump_history_enroll(Trpsub, "<Trpsub>");
        dump_history_enroll(Trpsub2, "<Trpsub2>");
#ifdef MHD
        dump_history_enroll(Br, "<Br>");
        dump_history_enroll(Bp, "<Bp>");
        dump_history_enroll(Bz, "<Bz>");
        dump_history_enroll(Br2, "<Br2>");
        dump_history_enroll(Mrp, "<Mrp>");
        dump_history_enroll(Mrpsub, "<Mrpsub>");
        dump_history_enroll(Mrpsub2, "<Mrpsub2>");
        dump_history_enroll(Pb, "<Pb>");
        dump_history_enroll(Pbsub, "<Pbsub>");
        dump_history_enroll(Pbsub2, "<Pbsub2>");
#endif
        dump_history_enroll(Pg, "<Pg>");
        dump_history_enroll(Pgsub, "<Pgsub>");
        dump_history_enroll(Pgsub2, "<Pgsub2>");

  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
	if (strcmp(expr,"Pb")==0) return Pb;
        else if(strcmp(expr,"Pg")==0) return Pg;
	else if(strcmp(expr,"Br2")==0) return Br2;
	else if(strcmp(expr,"Mrp")==0) return Mrp;
	else if(strcmp(expr,"Trp")==0) return Trp;
	else if(strcmp(expr,"rdVp")==0) return rdVp;
        else if(strcmp(expr,"Va")==0) return Va;
	else if(strcmp(expr,"Vaz")==0) return Vaz;
	else if(strcmp(expr,"Br")==0) return Br;
	else if(strcmp(expr,"Bp")==0) return Bp;
	else if(strcmp(expr,"Bz")==0) return Bz;
	else if(strcmp(expr,"Vr")==0) return Vr;
	else if(strcmp(expr,"Vp")==0) return Vp;
	else if(strcmp(expr,"Vz")==0) return Vz;
        else if(strcmp(expr,"L_AM")==0) return L_AM;
        else if(strcmp(expr,"Mdot")==0) return Mdot;
	else if(strcmp(expr,"Msub")==0) return Msub;
        else if(strcmp(expr,"Msub2")==0) return Msub2;
        else if(strcmp(expr,"Trpsub")==0) return Trpsub;
        else if(strcmp(expr,"Trpsub2")==0) return Trpsub2;
	else if(strcmp(expr,"Mrpsub")==0) return Mrpsub;
        else if(strcmp(expr,"Mrpsub2")==0) return Mrpsub2;
	else if(strcmp(expr,"Bpsub")==0) return Bpsub;
	else if(strcmp(expr,"Bzsub")==0) return Bzsub;
	else if(strcmp(expr,"Pbsub")==0) return Pbsub;
        else if(strcmp(expr,"Pbsub2")==0) return Pbsub2;
        else if(strcmp(expr,"Pgsub")==0) return Pgsub;
        else if(strcmp(expr,"Pgsub2")==0) return Pgsub2;
        else if(strcmp(expr,"Temp")==0) return Temp;

  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
	if(strcmp(name,"vtksub")==0) {
		return dump_vtksub;
	} else if(strcmp(name,"ktab")==0) {
		return out_ktab;
	} else if(strcmp(name,"jtab")==0) {
                return out_jtab;
	}
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{
  GridS *pG;
  int i,j,k,nl,nd;
  int nx1,nx2,nx3; 
  int is,ie,il,iu,js,je,jl,ju,ks,ke,kl,ku;
  Real T, Pres;

    for (nl=0; nl<(pM->NLevels); nl++){
      for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
        if (pM->Domain[nl][nd].Grid != NULL){
		pG = pM->Domain[nl][nd].Grid;

		is = pG->is;  ie = pG->ie;  nx1 = ie-is+1;
		js = pG->js;  je = pG->je;  nx2 = je-js+1;
		ks = pG->ks;  ke = pG->ke;  nx3 = ke-ks+1;

                il = is-nghost*(nx1>1);  iu = ie+nghost*(nx1>1);
                jl = js-nghost*(nx2>1);  ju = je+nghost*(nx2>1);
                kl = ks-nghost*(nx3>1);  ku = ke+nghost*(nx3>1);

		  for (k=kl; k<=ku; k++) {
		    for (j=jl; j<=ju; j++) {
		      for (i=il; i<=iu; i++) {     
			if(pG->U[k][j][i].d < d_MIN) pG->U[k][j][i].d = d_MIN;
			#ifndef BAROTROPIC
			#ifdef MHD
			Pres = pG->U[k][j][i].E - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
						- 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
			if(Pres < presmin)
			pG->U[k][j][i].E = presmin + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d
						   + 0.5*(SQR(pG->U[k][j][i].B1c) + SQR(pG->U[k][j][i].B2c) + SQR(pG->U[k][j][i].B3c));
			#else
			Pres = pG->U[k][j][i].E - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
			if(Pres < presmin)
			pG->U[k][j][i].E = presmin + 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d;
			#endif	

/*			T = pG->U[k][j][i].E/pG->U[k][j][i].d 
			- 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/SQR(pG->U[k][j][i].d);
			if(T < Tmin)
			pG->U[k][j][i].E = pG->U[k][j][i].d*Tmin 
			+ 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d; 
*/			
			#endif /*BAROTROPIC*/
		      }
		    }	
		  }
	}
      }
    }

}

void Userwork_after_loop(MeshS *pM)
{
}
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
        Pgas += (Gamma-1)*(pG->U[k][j][i].E - 0.5*(SQR(pG->U[k][j][i].M1) + SQR(pG->U[k][j][i].M2) + SQR(pG->U[k][j][i].M3))/pG->U[k][j][i].d)*CellVol;
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
#endif /*PARALLEL*/
  /*calculate and use scaling factor so that the correct beta is ensured*/
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
	/* Calculate angular momentum in rotating frame */
#ifdef FARGO
    	Lper = R*pGrid->U[k][j][is].M2;
#else
	/*Vkep = sqrt(R*(*x1GravAcc)(R,p,z));*/
	Vkep = vphi(R,p,z); /*different with or without ROTATING_FRAME, added by Wendy Ju*/
    	Lper = R*pGrid->U[k][j][is].M2 - R*pGrid->U[k][j][is].d*Vkep;
#endif /* FARGO */
    	Lper = MIN(Lper,0.0);

      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][is-i] = pGrid->U[k][j][is];
				/* Enforce diode */
				pGrid->U[k][j][is-i].M1 = MIN(pGrid->U[k][j][is-i].M1,0.0);
				/* Calculate Keplerian velocity */
				cc_pos(pGrid,is-i,j,k,&R,&p,&z);
				/* Vkep = sqrt(R*(*x1GravAcc)(R,p,z)); */
				Vkep = vphi(R,p,z); /*different with or without ROTATING_FRAME, added by Wendy Ju */
#ifdef FARGO
				if (Hbc == 1) {
					pGrid->U[k][j][is-i].M2 = 0.0;
				} else if (Hbc == 2){  /*copy the inner velocity*/
                                        pGrid->U[k][j][is-i].M2 = pGrid->U[k][j][is].M2;
                                }else{
					pGrid->U[k][j][is-i].M2 = Lper/R;
				}
#else
				if (Hbc == 1) {  /*local kepler velocity*/
					pGrid->U[k][j][is-i].M2 = pGrid->U[k][j][is].d*Vkep;
				} else if (Hbc == 2){  /*copy the inner velocity*/
				        pGrid->U[k][j][is-i].M2 = pGrid->U[k][j][is].M2;
				}else{  /*keep velocity perturbation constant*/
					pGrid->U[k][j][is-i].M2 = Lper/R + pGrid->U[k][j][is-i].d*Vkep;
				}
#endif
#ifdef MHD
				if (Mbc == 0) {
					pGrid->U[k][j][is-i].B1c = RBr/R;
					pGrid->U[k][j][is-i].B2c = 0.0;
					pGrid->U[k][j][is-i].B3c = 0.0;
				}
				else{
					pGrid->U[k][j][is-i].B1c = pGrid->U[k][j][is].B1c;
                                        pGrid->U[k][j][is-i].B2c = pGrid->U[k][j][is].B2c;
                                        pGrid->U[k][j][is-i].B3c = pGrid->U[k][j][is].B3c;
				}
#endif
#ifndef ISOTHERMAL
                                pGrid->U[k][j][is-i].E = pGrid->U[k][j][is].E;
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
  int jmid;
#ifdef MHD
  int ju, ku; /* j-upper, k-upper */
#endif
  Real Vkep, R,p,z, RBr, Lper;
#ifdef L1_INFLOW
  Real Lx, Ly, A0, Rmov; 
#endif

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
#ifdef MHD
    	RBr = pGrid->ri[ie+1]*pGrid->B1i[k][j][ie+1];
#endif
	/* Calculate angular momentum in rotating frame */
#ifdef FARGO
    	Lper = pGrid->r[ie]*pGrid->U[k][j][ie].M2;
#else
			cc_pos(pGrid,ie,j,k,&R,&p,&z);
			Vkep = vphi(R,p,z); /*different with or without ROTATING_FRAME, added by Wendy Ju*/
			Lper = R*pGrid->U[k][j][ie].M2 - R*pGrid->U[k][j][ie].d*Vkep;
#endif
      for (i=1; i<=nghost; i++) {
        pGrid->U[k][j][ie+i] = pGrid->U[k][j][ie];
				/* Enforce diode */
				pGrid->U[k][j][ie+i].M1 = MAX(pGrid->U[k][j][ie+i].M1,0.0);
				cc_pos(pGrid,ie+i,j,k,&R,&p,&z);
/*				Vkep = sqrt(R*(*x1GravAcc)(R,p,z));

*/
				Vkep = vphi(R,p,z); /*different with or without ROTATING_FRAME, added by Wendy Ju*/
#ifdef FARGO
				if (Hbc == 1) {
					pGrid->U[k][j][ie+i].M2 = 0.0;
				} else if (Hbc == 2){  /*copy the inner velocity*/
                                        pGrid->U[k][j][ie+i].M2 = pGrid->U[k][j][ie].M2;
                                }else{
					pGrid->U[k][j][ie+i].M2 = Lper/R;
				}
#else
				if (Hbc == 1) {  /*local kepler velocity*/
					pGrid->U[k][j][ie+i].M2 = pGrid->U[k][j][ie].d*Vkep;
				} else if (Hbc == 2){  /*copy the inner velocity*/
				        pGrid->U[k][j][ie+i].M2 = pGrid->U[k][j][ie].M2;
				}else {  /*keep velocity perturbation constant*/
					pGrid->U[k][j][ie+i].M2 = Lper/R + pGrid->U[k][j][ie+i].d*Vkep;
				}
#endif
#ifdef MHD
				if (Mbc == 0) {
					pGrid->U[k][j][ie+i].B1c = RBr/R;
					pGrid->U[k][j][ie+i].B2c = 0.0;
					pGrid->U[k][j][ie+i].B3c = 0.0;
				}
				else{
					pGrid->U[k][j][ie+i].B1c = pGrid->U[k][j][ie].B1c;
                                        pGrid->U[k][j][ie+i].B2c = pGrid->U[k][j][ie].B2c;
                                        pGrid->U[k][j][ie+i].B3c = pGrid->U[k][j][ie].B3c;
				}
#endif
#ifndef ISOTHERMAL
                                pGrid->U[k][j][ie+i].E = pGrid->U[k][j][ie].E;
#endif
	
      }/*i*/
    }/*j*/
  }/*k*/

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


#ifdef L1_INFLOW
  cc_pos(pGrid,ie,js,ks,&R,&p,&z);
  Vkep = R * Omega(R);

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) { 
      for (i=1; i<=nghost; i++) {
        cc_pos(pGrid,ie+i,j,k,&R,&p,&z);

	/* To let the geometry move with the flow*/
	Rmov = pGrid->r[ie+i] + Vr_l1*pGrid->time;
	if((p>-Pi/60.0)&&(p<Pi/60.0)) { /*about 6 degrees*/
	pGrid->U[k][j][ie+i].d = 0.5 * rho0;
	pGrid->U[k][j][ie+i].M1 = -1.0 * pGrid->U[k][j][ie+i].d * Vr_l1;
#ifdef ROTATING_FRAME
	pGrid->U[k][j][ie+i].M2 = 0.0;
#else
	pGrid->U[k][j][ie+i].M2 = R * Omega(R);
#endif
#ifndef ISOTHERMAL
        pGrid->U[k][j][ie+i].E  = 10.0*pgas/(Gamma*Gamma_1) 
				+ 0.5*(SQR(pGrid->U[k][j][ie+i].M1) + SQR(pGrid->U[k][j][ie+i].M2) + SQR(pGrid->U[k][j][ie+i].M3))/pGrid->U[k][j][ie+i].d;
#endif

#ifdef MHD /** Toroidal Inflow **/
	/*Diameter of the toroidal*/
	Lx = 4.0*pGrid->dx1;
	Ly = Lx;
	#ifndef BAROTROPIC
	A0 = sqrt(2.0*10.0*pgas/Gamma/Beta);
	#else
	A0 = sqrt(2.0*pGrid->U[k][j][ie+i].d*Iso_csound2/Beta);
	#endif /*BAROTROPIC*/
        pGrid->U[k][j][ie+i].B1c = A0 * cos(Pi/Lx*Rmov*cos(p)) * sin(Pi/Ly*Rmov*sin(p)) * sin(p)
				 + A0 * sin(Pi/Lx*Rmov*cos(p)) * cos(Pi/Ly*Rmov*sin(p)) * cos(p);
        pGrid->U[k][j][ie+i].B2c = A0 * cos(Pi/Lx*Rmov*cos(p)) * sin(Pi/Ly*Rmov*sin(p)) * cos(p)
                                 + A0 * sin(Pi/Lx*Rmov*cos(p)) * cos(Pi/Ly*Rmov*sin(p)) * sin(p);
	pGrid->U[k][j][ie+i].B3c = 0.0;

	if(i!=1){
        pGrid->B1i[k][j][ie+i] = A0 * cos(Pi/Lx*(Rmov-0.5*pGrid->dx1)*cos(p)) * sin(Pi/Ly*(Rmov-0.5*pGrid->dx1)*sin(p)) * sin(p)
                                 + A0 * sin(Pi/Lx*(Rmov-0.5*pGrid->dx1)*cos(p)) * cos(Pi/Ly*(Rmov-0.5*pGrid->dx1)*sin(p)) * cos(p);
        pGrid->B2i[k][j][ie+i] = A0 * cos(Pi/Lx*(Rmov-0.5*pGrid->dx1)*cos(p)) * sin(Pi/Ly*(Rmov-0.5*pGrid->dx1)*sin(p)) * cos(p)
                                 + A0 * sin(Pi/Lx*(Rmov-0.5*pGrid->dx1)*cos(p)) * cos(Pi/Ly*(Rmov-0.5*pGrid->dx1)*sin(p)) * sin(p);
        pGrid->B3i[k][j][ie+i] = 0.0;
	}
        #ifndef ISOTHERMAL
           pGrid->U[k][j][ie+i].E  += 0.5*(SQR(pGrid->U[k][j][ie+i].B1c) + SQR(pGrid->U[k][j][ie+i].B2c) + SQR(pGrid->U[k][j][ie+i].B3c));
        #endif

#endif /*MHD*/
	}/*if p*/
      }/*i*/
    }/*j*/
  }/*k*/

#endif /*L1_INFLOW*/

  return;

}

static void dump_vtksub(MeshS *pM, OutputS *pOut) {
	GridS *pGrid;
	FILE *pfile;
	char *fname, *plev = NULL, *pdom = NULL;
	char levstr[8], domstr[8];
	/* Upper and Lower bounds on i,j,k for data dump */
	int i, j, k, il, iu, jl, ju, kl, ku;
	int big_end = ath_big_endian();
	int ndata0;
	float *data; /* points to 3*ndata0 allocated floats */
	double x1m, x2m, x3m, x1M, x1, x2, x3, dR;

	/* Loop over all Domains in Mesh, and output Grid data */

	pGrid = pM->Domain[0][0].Grid;

	il = pGrid->is, iu = pGrid->ie;
	jl = pGrid->js, ju = pGrid->je;
	kl = pGrid->ks, ku = pGrid->ke;
	x1m = pGrid->MinX[0];
	x2m = pGrid->MinX[1];
	x3m = pGrid->MinX[2];

	x1M = pGrid->MaxX[0];
	dR = pGrid->dx1;
	/* Exclude tiles not in the subvolume */
	if ((x1M <= 1.6) || (x1m >= 2.4)) {
		return;
	}
	il = -1;
	for (i = pGrid->is; i <= pGrid->ie; i++) {
		cc_pos(pGrid, i, jl, kl, &x1, &x2, &x3);
		/* Get first index in range */
		if ((x1 + 0.5 * dR >= 1.6) && (il < 0)) {
			il = i;
			x1m = x1 - 0.5 * dR;
		}
		if (x1 - 0.5 * dR <= 2.4) {
			iu = i;
		}
	}

	ndata0 = iu - il + 1;

	/* construct filename, open file */
	if ((fname = ath_fname(plev, pM->outfilename, plev, pdom, num_digit,
			pOut->num, "sub", "vtk")) == NULL) {
		ath_error("[dump_vtk]: Error constructing filename\n");
	}

	if ((pfile = fopen(fname, "w")) == NULL) {
		ath_error("[dump_vtk]: Unable to open vtk dump file\n");
		return;
	}

	/* Allocate memory for temporary array of floats */

	if ((data = (float *) malloc(3 * ndata0 * sizeof(float))) == NULL) {
		ath_error("[dump_vtk]: malloc failed for temporary array\n");
		return;
	}

	/* There are five basic parts to the VTK "legacy" file format.  */
	/*  1. Write file version and identifier */

	fprintf(pfile, "# vtk DataFile Version 2.0\n");

	/*  2. Header */

	fprintf(pfile, "Subvolume variables at time= %e, level= %i, domain= %i\n",
			pGrid->time, 0, 0);

	/*  3. File format */

	fprintf(pfile, "BINARY\n");

	/*  4. Dataset structure */

	/* Set the Grid origin */

	fprintf(pfile, "DATASET STRUCTURED_POINTS\n");
	if (pGrid->Nx[1] == 1) {
		fprintf(pfile, "DIMENSIONS %d %d %d\n", iu - il + 2, 1, 1);
	} else {
		if (pGrid->Nx[2] == 1) {
			fprintf(pfile, "DIMENSIONS %d %d %d\n", iu - il + 2, ju - jl + 2, 1);
		} else {
			fprintf(pfile, "DIMENSIONS %d %d %d\n", iu - il + 2, ju - jl + 2, ku - kl
					+ 2);
		}
	}
	fprintf(pfile, "ORIGIN %e %e %e \n", x1m, x2m, x3m);
	fprintf(pfile, "SPACING %e %e %e \n", pGrid->dx1, pGrid->dx2, pGrid->dx3);

	/*  5. Data  */

	fprintf(pfile, "CELL_DATA %d \n", (iu - il + 1) * (ju - jl + 1) * (ku - kl
			+ 1));

	/* Write density */

	fprintf(pfile, "SCALARS density float\n");
	fprintf(pfile, "LOOKUP_TABLE default\n");
	for (k = kl; k <= ku; k++) {
		for (j = jl; j <= ju; j++) {
			for (i = il; i <= iu; i++) {
				data[i - il] = (float) pGrid->U[k][j][i].d;
			}
			if (!big_end)
				ath_bswap(data, sizeof(float), iu - il + 1);
			fwrite(data, sizeof(float), (size_t) ndata0, pfile);
		}
	}

	/* Write momentum or velocity */

	fprintf(pfile, "\nVECTORS momentum float\n");
	for (k = kl; k <= ku; k++) {
		for (j = jl; j <= ju; j++) {
			for (i = il; i <= iu; i++) {
				data[3 * (i - il)] = (float) pGrid->U[k][j][i].M1;
				data[3 * (i - il) + 1] = (float) pGrid->U[k][j][i].M2;
				data[3 * (i - il) + 2] = (float) pGrid->U[k][j][i].M3;
			}
			if (!big_end)
				ath_bswap(data, sizeof(float), 3 * (iu - il + 1));
			fwrite(data, sizeof(float), (size_t) (3 * ndata0), pfile);
		}
	}

	/* Write cell centered B */

#ifdef MHD
	fprintf(pfile, "\nVECTORS cell_centered_B float\n");
	for (k = kl; k <= ku; k++) {
		for (j = jl; j <= ju; j++) {
			for (i = il; i <= iu; i++) {
				data[3 * (i - il)] = (float) pGrid->U[k][j][i].B1c;
				data[3 * (i - il) + 1] = (float) pGrid->U[k][j][i].B2c;
				data[3 * (i - il) + 2] = (float) pGrid->U[k][j][i].B3c;
			}
			if (!big_end)
				ath_bswap(data, sizeof(float), 3 * (iu - il + 1));
			fwrite(data, sizeof(float), (size_t) (3 * ndata0), pfile);
		}
	}

#endif

	/* close file and free memory */

	fclose(pfile);
	free(data);

	return;
}

void out_ktab(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid=pM->Domain[0][0].Grid;
  int i,nx1;
  FILE *pFile;
  char fmt[80],*fname,*plev=NULL,*pdom=NULL;
  char levstr[8],domstr[8];
  Real *data=NULL;
  Real dmin, dmax, xworld;

printf("out_ktab used\n");

/* Add a white space to the format, setup format for integer zone columns */
  if(pOut->dat_fmt == NULL){
     sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

/* compute 1D array of data */
  data = OutData1(pGrid,pOut,&nx1);
  if (data == NULL) { printf("error: slice not in range of Grid\n"); return;}  /* slice not in range of Grid */

  minmax1(data,nx1,&dmin,&dmax);


  if((fname = ath_fname(plev,pM->outfilename,0,0,0,0,
      pOut->id,"tab")) == NULL){
    ath_error("[output_tab]: Error constructing filename\n");
  }

  pFile = fopen(fname,"a");
/* open filename */
  if (pFile == NULL) {
    ath_error("[output_tab]: Unable to open tab file %s\n",fname);
  }

/* write data */

  if (pOut->num == 0) {
  	for (i=0; i<nx1; i++) {
  		fprintf(pFile,"%12d\t",pGrid->Disp[0]+i);
  	}
  	fprintf(pFile,"\n");
  }
  for (i=0; i<nx1; i++) {
  	fprintf(pFile,fmt,data[i]);
  	fprintf(pFile,"\t");
  }
  fprintf(pFile,"\n");
/* Compute and store global min/max, for output at end of run */
  pOut->gmin = MIN(dmin,pOut->gmin);
  pOut->gmax = MAX(dmax,pOut->gmax);

  fclose(pFile);
  free_1d_array(data); /* Free the memory we malloc'd */
}


/*! \fn static void out_jtab(MeshS *pM, OutputS *pOut)
 *  \brief output routine to calculate 1D azimuthally & vertically
    averaged quantities. */

static void out_jtab(MeshS *pM, OutputS *pOut)
{
  GridS *pGrid;
  DomainS *pD;
  int i,j,k;
  int tot1d,i1d,nx,ig,idisp;
  int dnum = pOut->num,nl,nd;
  static int FIRST = 0;
  double darea,**out1d;
  double x1,x2,x3,Lphi,Lz;
  static double *out_r;

  FILE *p_1dfile;
  char *fname, fmt[80];
  double area_rat; /* (Grid Volume)/(dx1*dx2*dx3) */

  if(pOut->dat_fmt == NULL){
     sprintf(fmt," %%12.8e"); /* Use a default format */
  }
  else{
    sprintf(fmt," %s",pOut->dat_fmt);
  }

#ifdef MPI_PARALLEL
  double *my_out1d;
  double *g_out1d;
  int ierr,myID_Comm_Domain;
#endif

#ifdef MHD
  tot1d=15;
#else
  tot1d=7;
#endif /* MHD */
#ifdef ADIABATIC
  tot1d=tot1d+2;
#endif /* ADIABATIC */

  Lphi = pM->RootMaxX[1] - pM->RootMinX[1];
  Lz = pM->RootMaxX[2] - pM->RootMinX[2];
  nx = pM->Nx[0];

/* At level=0, there is only one domain */

  pGrid = pM->Domain[0][0].Grid;
  int is = pGrid->is, ie = pGrid->ie;
  int js = pGrid->js, je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke;
  pD = (DomainS*)&(pM->Domain[0][0]);

#ifdef MPI_PARALLEL
  ierr = MPI_Comm_rank(pD->Comm_Domain, &myID_Comm_Domain);
  if(ierr != MPI_SUCCESS)
    ath_error("[change_rundir]: MPI_Comm_rank error = %d\n",ierr);
#endif
  if (FIRST == 0){
#ifdef MPI_PARALLEL
    if (myID_Comm_Domain == 0) {
#endif
      out_r = (double *) calloc_1d_array(nx,sizeof(double));
#ifdef MPI_PARALLEL
    }
#endif
  }

  out1d = (double **) calloc_2d_array(nx,tot1d,sizeof(double));
#ifdef MPI_PARALLEL
  my_out1d = (double *) calloc_1d_array(nx,sizeof(double));
  g_out1d = (double *) calloc_1d_array(nx,sizeof(double));
#endif
  for (i=0; i<nx; i++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[i][i1d] = 0.0;
    }
  }
  idisp=pGrid->Disp[0]; 

/* First calculate the x1 coordinate and save it to be dumped
   by root in every 1d file */
  if (FIRST == 0) {
#ifdef MPI_PARALLEL
  if (myID_Comm_Domain == 0) {
#endif
    for (i=0; i<nx; i++) {
      x1 = pM->RootMinX[0] + (i + 0.5)*pGrid->dx1;
      out_r[i] = x1;
    }
#ifdef MPI_PARALLEL
  }
#endif
  }

/* Compute 1d averaged variables */
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
	ig=i+idisp-nghost;
        i1d=0;
        out1d[ig][i1d] += pGrid->U[k][j][i].d; /*density*/
        i1d++;
        out1d[ig][i1d] += Pg(pGrid,i,j,k); /*gas pressure*/
#ifdef MHD
        i1d++;
        out1d[ig][i1d] += Pb(pGrid,i,j,k); /*magnetic pressure*/
#endif
#ifdef ADIABATIC
        i1d++;
        out1d[ig][i1d] += Temp(pGrid,i,j,k); /*gas temperature*/
        i1d++;
        out1d[ig][i1d] += pGrid->U[k][j][i].E; /*total energy*/
#endif
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].M1)/pGrid->U[k][j][i].d; /*KEr*/
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(Vp(pGrid,i,j,k))*pGrid->U[k][j][i].d; /*KEphi*/
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].M3)/pGrid->U[k][j][i].d; /*KEz*/
        i1d++;
        out1d[ig][i1d] += Trp(pGrid,i,j,k); /*Reynolds stress*/
#ifdef MHD
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B1c); /*MEr*/
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B2c); /*MEphi*/
        i1d++;
        out1d[ig][i1d] += 0.5*SQR(pGrid->U[k][j][i].B3c); /*MEz*/
        i1d++;
        out1d[ig][i1d] += Br(pGrid,i,j,k); /*Br*/
        i1d++;
        out1d[ig][i1d] += Bp(pGrid,i,j,k); /*Bphi*/
        i1d++;
        out1d[ig][i1d] += Bz(pGrid,i,j,k); /*Bz*/
        i1d++;
        out1d[ig][i1d] += Mrp(pGrid,i,j,k); /*Mrp*/
#endif
        i1d++;
        out1d[ig][i1d] += Mdot(pGrid,i,j,k); /*Mdot*/
      }
    }
  }

  /* Calculate the (Grid Volume) / (Grid Cell Volume) Ratio */
  area_rat = Lphi*Lz/(pGrid->dx2*pGrid->dx3);

/* The parent sums the scal[] array.
 * Note that this assumes (dx1,dx2,dx3) = const. */

#ifdef MPI_PARALLEL 
  for(i1d=0; i1d<tot1d; i1d++){
    for (i=0; i<nx; i++) {
      my_out1d[i] = out1d[i][i1d];
    }
    ierr = MPI_Reduce(my_out1d, g_out1d, nx,
                      MPI_DOUBLE, MPI_SUM, 0, pD->Comm_Domain);
    if(ierr)
      ath_error("[output_jtab]: MPI_Reduce call returned error = %d\n",ierr);
    for (i=0; i<nx; i++) {
      out1d[i][i1d] = g_out1d[i];
    }
  }
#endif

/* For parallel calculations, only the parent computes the average
 * and writes the output. */
#ifdef MPI_PARALLEL
  if(myID_Comm_Domain == 0){ /* I'm the parent */
#endif

  darea = 1.0/(double)area_rat;
  for (i=0; i<nx; i++) {
    for (i1d=0; i1d<tot1d; i1d++) {
      out1d[i][i1d] *= darea;
    }
  }

/* Generate filename */
#ifdef MPI_PARALLEL
  fname = ath_fname("../",pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"jtab");
#else
  fname = ath_fname(NULL,pM->outfilename,NULL,NULL,num_digit,dnum,NULL,"jtab");
#endif
  if (fname == NULL) {
    ath_error("[output_jtab]: Error constructing output filename\n");
    return;
  }

/* open filename */
  p_1dfile = fopen(fname,"w");
  if (p_1dfile == NULL) {
    ath_error("[output_jtab]: Unable to open 1d average file %s\n",fname);
    return;
  }

/* Write out data */

  for (i=0; i<nx; i++) {
#ifdef ISOTHERMAL
#ifdef MHD
    if (i == 0) {
      fprintf(p_1dfile,"# r               dens            Pg              Pb              KEr             KEphi           KEz             Trp             MEr             MEphi           MEz             Br              Bphi            Bz              Mrp             Mdot\n");
    }
#else
    if (i == 0) {
      fprintf(p_1dfile,"# r               dens            Pg              KEr             KEphi           KEz             Reynolds        Mdot\n");
    }
#endif /* MHD */
#else
#ifdef MHD
    if (i == 0) {
      fprintf(p_1dfile,"# r               dens            Pg              Pb              Temp            Etot            KEr             KEphi           KEz             Trp             MEr             MEphi           MEz              Br              Bphi            Bz              Mrp             Mdot\n");
    }
#else
    if (i == 0) {
      fprintf(p_1dfile,"# r               dens            Pg              Temp             Etot            KEr             KEphi          KEz             Reynolds        Mdot\n");
    }
#endif /* MHD */
#endif /* ISOTHERMAL */

    fprintf(p_1dfile, fmt, out_r[i]);
    fprintf(p_1dfile, "\t");
    for(i1d=0; i1d<tot1d; i1d++){
        fprintf(p_1dfile, fmt, out1d[i][i1d]);
        fprintf(p_1dfile, "\t");
    }
    fprintf(p_1dfile, "\n");

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

