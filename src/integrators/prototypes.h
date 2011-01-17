#ifndef INTEGRATORS_PROTOTYPES_H
#define INTEGRATORS_PROTOTYPES_H 
#include "../copyright.h"
/*==============================================================================
 * FILE: prototypes.h
 *
 * PURPOSE: Prototypes for all public functions in the /src/integrators dir
 *============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* integrate.c */
VDFun_t integrate_init(MeshS *pM);
void integrate_destruct(void);

/* Only used for rad_hydro integrators */
#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)
void BackEuler_init_1d(int Ngrids, MeshS *pM);
void BackEuler_init_2d(int Ngrids);
void BackEuler_destruct_1d(int Ngrids);
void BackEuler_destruct_2d();
VMFun_t BackEuler_init(MeshS *pM);
/* General LU decomposition functions */
void ludcmp(Real **a, int n, int *indx, Real *d);
void lubksb(Real **a, int n, int *indx, Real b[]);
/* Matrix solver for band diagnol equations */
void bandec(Real **a, unsigned long n, int m1, int m2, Real **al,
	unsigned long indx[], Real *d);
void banbks(Real **a, unsigned long n, int m1, int m2, Real **al,
	unsigned long indx[], Real b[]);

/* Matrix solver using restarted GMRES method to sovle sparse matrix */

void ax_st ( int n, int nz_num, int ia[], int ja[], double a[], double x[],
  double w[] );

void mult_givens ( double c, double s, int k, double *g );

double r8vec_dot ( int n, double a1[], double a2[] );

void mgmres_st ( int n, int nz_num, int ia[], int ja[], double a[], 
  double x[], double rhs[], int itr_max, int mr, double tol_abs, 
  double tol_rel );


#endif

/* integrate_1d_ctu.c and integrate_1d_vl.c */
void integrate_destruct_1d(void);
void integrate_init_1d(MeshS *pM);
void integrate_1d_ctu(DomainS *pD);
void integrate_1d_vl(DomainS *pD);

#if defined (RADIATION_HYDRO) || defined (RADIATION_MHD)
void integrate_1d_radMHD(DomainS *pD);
void integrate_2d_radMHD(DomainS *pD);
void integrate_3d_radMHD(DomainS *pD);

void integrate_1d_radMHD_thick(DomainS *pD);
void integrate_2d_radMHD_thick(DomainS *pD);
void integrate_3d_radMHD_thick(DomainS *pD);
void BackEuler_1d(MeshS *pM);
void BackEuler_2d(MeshS *pM);
#endif

/* integrate_2d_ctu.c and integrate_2d_vl.c */
void integrate_destruct_2d(void);
void integrate_init_2d(MeshS *pM);
void integrate_2d_ctu(DomainS *pD);
void integrate_2d_vl(DomainS *pD);


/* integrate_3d_ctu.c and integrate_3d_vl.c */
void integrate_destruct_3d(void);
void integrate_init_3d(MeshS *pM);
void integrate_3d_ctu(DomainS *pD);
void integrate_3d_vl(DomainS *pD);


#endif /* INTEGRATORS_PROTOTYPES_H */
