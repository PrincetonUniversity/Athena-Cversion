#include "copyright.h"
/*==============================================================================
 * FILE: selfg_multigrid.c
 *
 * PURPOSE: Contains functions to solve Poisson's equation for self-gravity in
 *   3D using multigrid.
 *
 *   These functions work for non-periodic domains.  A low-order multipole
 *   expansion is used to compute the potential on the boundaries.
 *
 * HISTORY:
 *   june-2007 - 2D and 3D solvers written by Irene Balmes
 *   july-2007 - routines incorporated into Athena by JMS and IB
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   selfg_by_multig_2d() - 2D Poisson solver using multigrid
 *   selfg_by_multig_3d() - 3D Poisson solver using multigrid
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/* MGrid structure; holds RHS, potential, and information about grid
 * size for a given level in the multi-grid hierarchy  */
typedef struct MGrid_s{
  Real ***rhs,***Phi;  /* RHS of elliptic equation, and solution */
  Real dx1,dx2,dx3;
  int Nx1,Nx2,Nx3;
  int is,ie;
  int js,je;
  int ks,ke; 
}MGrid;

/* 3D temporary array needed for restriction of errors  */
Real ***error;

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   multig_3d() -
 *============================================================================*/

void multig_3d(MGrid *pMG);
void Jacobi(MGrid *pMG);
void Restriction_3d(MGrid *pMG_fine, MGrid *pMG_coarse);
void Prolongation_3d(MGrid *pMG_coarse, MGrid *pMG_fine);

#ifdef MPI_PARALLEL
void set_iterate_bvals(MGrid *pMG);
void swap_iterate_ix1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_iterate_ox1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_iterate_ix2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_iterate_ox2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_iterate_ix3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
void swap_iterate_ox3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq);
#endif


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* selfg_by_multig_1d:
 */

void selfg_by_multig_1d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_2d:
 */

void selfg_by_multig_2d(Grid *pG, Domain *pD)
{
  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_3d:  Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */

void selfg_by_multig_3d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY_USING_MULTIGRID
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int Nx1z, Nx2z, Nx3z;
  MGrid Root_grid;
  Real mass = 0.0, dVol, rad, x1, x2, x3;
  Real Grav_const = four_pi_G/(4.0*PI);

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
      }
    }
  }

/* Allocate memory for error array used during restriction */

  Nx1z = pG->Nx1 + 2*nghost;
  Nx2z = pG->Nx2 + 2*nghost;
  Nx3z = pG->Nx3 + 2*nghost;
  error = (Real ***) calloc_3d_array(Nx3z, Nx2z, Nx1z, sizeof(Real));
  if (error == NULL) {
    ath_error("[multig_3d]: Error allocating memory for error array\n");
  }

/* Compute solution at boundaries using monopole expansion */

  dVol = pG->dx1*pG->dx2*pG->dx3;
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        mass += pG->U[k][j][i].d*dVol;
      }
    }
  }

/*  Inner and outer x1 boundaries */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=1; i<=nghost; i++){
        cc_pos(pG,is-i,j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][j][is-i] = -Grav_const*mass/rad;

        cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][j][ie+i] = -Grav_const*mass/rad;
      }
    }
  }

/*  Inner and outer x2 boundaries */

  for (k=ks; k<=ke; k++){
    for (j=1; j<=nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        cc_pos(pG,i,js-j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][js-j][i] = -Grav_const*mass/rad;

        cc_pos(pG,i,je+j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][je+j][i] = -Grav_const*mass/rad;
      }
    }
  }

/*  Inner and outer x3 boundaries */

  for (k=1; k<=nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[ks-k][j][i] = -Grav_const*mass/rad;

        cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[ke+k][j][i] = -Grav_const*mass/rad;
      }
    }
  }

/* Initialize MGrid structure for top level (the root, or finest, level) */

  Nx1z = pG->Nx1 + 2;
  Nx2z = pG->Nx2 + 2;
  Nx3z = pG->Nx3 + 2;

  Root_grid.Nx1 = pG->Nx1;
  Root_grid.Nx2 = pG->Nx2;
  Root_grid.Nx3 = pG->Nx3;
  Root_grid.is = 1;  Root_grid.ie = pG->Nx1;
  Root_grid.js = 1;  Root_grid.je = pG->Nx2;
  Root_grid.ks = 1;  Root_grid.ke = pG->Nx3;
  Root_grid.dx1 = pG->dx1;
  Root_grid.dx2 = pG->dx2;
  Root_grid.dx3 = pG->dx3;

  Root_grid.rhs = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
  Root_grid.Phi = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
  if (Root_grid.rhs == NULL) {
    ath_error("[selfg_by_multig_3d]: Error allocating memory\n");
  }
  if (Root_grid.Phi == NULL) {
    ath_error("[selfg_by_multig_3d]: Error allocating memory\n");
  }

  for (k=ks-1; k<=ke+1; k++){
    for (j=js-1; j<=je+1; j++){
      for (i=is-1; i<=ie+1; i++){
        Root_grid.rhs[k-ks+1][j-js+1][i-is+1] = four_pi_G*pG->U[k][j][i].d;
        Root_grid.Phi[k-ks+1][j-js+1][i-is+1] = pG->Phi[k][j][i];
      }
    }
  }

/* Compute new potential.  Note multig_3d calls itself recursively. */

  multig_3d(&Root_grid);

/* copy solution for potential from MGrid into Grid structure */

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        pG->Phi[k][j][i] = Root_grid.Phi[k-ks+1][j-js+1][i-is+1];
      }
    }
  }

#endif /* SELF_GRAVITY_USING_MULTIGRID */
  return;
}

/*----------------------------------------------------------------------------*/
/* Functions needed for the multigrid solver in 3D
 */
#ifdef SELF_GRAVITY_USING_MULTIGRID

void multig_3d(MGrid *pMG)
{
  int i, is = pMG->is, ie = pMG->ie;
  int j, js = pMG->js, je = pMG->je;
  int k, ks = pMG->ks, ke = pMG->ke;
  int Nx1z, Nx2z, Nx3z;
  MGrid Coarse_grid;

/* If we are down to 4 cells in any dimension do 10 iterations and return */

  if (pMG->Nx1==4 || pMG->Nx2==4 || pMG->Nx3==4)
    Jacobi(pMG);

/* Else, do 10 iterations at this level, restrict to a coarser grid, and call
 * multig_3d again with this coarse grid */

  else { 
    Jacobi(pMG);

/* Allocate and initialize MGrid at next coarsest level */

    Nx1z = (pMG->Nx1)/2 + 2;
    Nx2z = (pMG->Nx2)/2 + 2;
    Nx3z = (pMG->Nx3)/2 + 2;

    Coarse_grid.Nx1 = pMG->Nx1/2;
    Coarse_grid.Nx2 = pMG->Nx2/2;
    Coarse_grid.Nx3 = pMG->Nx3/2;
    Coarse_grid.is = 1;  Coarse_grid.ie = Coarse_grid.Nx1;
    Coarse_grid.js = 1;  Coarse_grid.je = Coarse_grid.Nx2;
    Coarse_grid.ks = 1;  Coarse_grid.ke = Coarse_grid.Nx3;
    Coarse_grid.dx1 = 2.0*pMG->dx1;
    Coarse_grid.dx2 = 2.0*pMG->dx2;
    Coarse_grid.dx3 = 2.0*pMG->dx3;

    Coarse_grid.rhs = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
    Coarse_grid.Phi = (Real ***) calloc_3d_array(Nx3z,Nx2z,Nx1z,sizeof(Real));
    if (Coarse_grid.rhs == NULL) {
      ath_error("[multig_3d]: Error allocating memory for some level\n");
    }
    if (Coarse_grid.Phi == NULL) {
      ath_error("[multig_3d]: Error allocating memory for some level\n");
    }

    Restriction_3d(pMG, &Coarse_grid);
#ifdef MPI_PARALLEL
    set_iterate_bvals(pMG);
#endif

    multig_3d(&Coarse_grid);

/* The following code is first reached after 10 iterations at the coarsest
 * level.  We then prolongate, do 10 iterations, and return.  This will return
 * execution to this same spot for the next coarsest level, so we will
 * prolongate, do 10 iterations, return, and so on.
 */

    Prolongation_3d(&Coarse_grid, pMG);
#ifdef MPI_PARALLEL
    set_iterate_bvals(pMG);
#endif

/* End with 10 iterations at the finest level */

    Jacobi(pMG);
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_3d:  Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */
void Jacobi(MGrid *pMG)
{
  int i, is = pMG->is, ie = pMG->ie;
  int j, js = pMG->js, je = pMG->je;
  int k, ks = pMG->ks, ke = pMG->ke;
  int n;
  Real trhs;
  Real dx1sq = (pMG->dx1*pMG->dx1);
  Real dx2sq = (pMG->dx2*pMG->dx2);
  Real dx3sq = (pMG->dx3*pMG->dx3);

/* Jacobi iterations in 3D */

  if (pMG->Nx3 > 1) {
    for (n=0; n<=10; n++){  /* hardwired to do 10 iterations */
      for (k=ks; k<=ke; k++){
        for (j=js; j<=je; j++){
          for (i=is; i<=ie; i++){
            trhs = -pMG->rhs[k][j][i];
            trhs += (pMG->Phi[k][j][i-1] + pMG->Phi[k][j][i+1])/dx1sq;
            trhs += (pMG->Phi[k][j+1][i] + pMG->Phi[k][j-1][i])/dx2sq;
            trhs += (pMG->Phi[k+1][j][i] + pMG->Phi[k-1][j][i])/dx3sq;
            pMG->Phi[k][j][i] = trhs/(2.0/dx1sq + 2.0/dx2sq + 2.0/dx3sq);
          }
        }
      }
#ifdef MPI_PARALLEL
      set_iterate_bvals(pMG);
#endif
    }
  }

/* Jacobi iterations in 2D (x-y plane) */

  if (pMG->Nx3 == 1) {
    for (n=0; n<=10; n++){  /* hardwired to do 10 iterations */
      for (j=js; j<=je; j++){
        for (i=is; i<=ie; i++){
          trhs = -pMG->rhs[ks][j][i];
          trhs += (pMG->Phi[ks][j][i-1] + pMG->Phi[ks][j][i+1])/dx1sq;
          trhs += (pMG->Phi[ks][j+1][i] + pMG->Phi[ks][j-1][i])/dx2sq;
          pMG->Phi[ks][j][i] = trhs/(2.0/dx1sq + 2.0/dx2sq);
        }
      }
#ifdef MPI_PARALLEL
      set_iterate_bvals(pMG);
#endif
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_3d:  Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */

void Restriction_3d(MGrid *pMG_fine, MGrid *pMG_coarse)
{
  int i, is = 1, ie = pMG_fine->ie;
  int j, js = 1, je = pMG_fine->je;
  int k, ks = 1, ke = pMG_fine->ke;
  Real dx1sq = (pMG_fine->dx1*pMG_fine->dx1);
  Real dx2sq = (pMG_fine->dx2*pMG_fine->dx2);
  Real dx3sq = (pMG_fine->dx3*pMG_fine->dx3);

  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        error[k][j][i] = pMG_fine->rhs[k][j][i];
        error[k][j][i] -= (pMG_fine->Phi[k][j][i+1] + pMG_fine->Phi[k][j][i-1]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx1sq;
        error[k][j][i] -= (pMG_fine->Phi[k][j+1][i] + pMG_fine->Phi[k][j-1][i]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx2sq;
        error[k][j][i] -= (pMG_fine->Phi[k+1][j][i] + pMG_fine->Phi[k-1][j][i]
          - 2.0*pMG_fine->Phi[k][j][i]) / dx3sq;
      }
    }
  }

  for(k=ks; k<=pMG_coarse->ke; k++){
    for (j=js; j<=pMG_coarse->je; j++){
      for (i=is; i<=pMG_coarse->ie; i++){
        pMG_coarse->rhs[k][j][i] =
           (error[2*k  ][2*j  ][2*i] + error[2*k  ][2*j  ][2*i-1]
          + error[2*k  ][2*j-1][2*i] + error[2*k  ][2*j-1][2*i-1]
          + error[2*k-1][2*j  ][2*i] + error[2*k-1][2*j  ][2*i-1]
          + error[2*k-1][2*j-1][2*i] + error[2*k-1][2*j-1][2*i-1])/8.0 ;
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_3d:  Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */
void Prolongation_3d(MGrid *pMG_coarse, MGrid *pMG_fine)
{
  int i, is = 1, ie = pMG_coarse->ie;
  int j, js = 1, je = pMG_coarse->je;
  int k, ks = 1, ke = pMG_coarse->ke;

  for (k=ks; k<=ke; k++){
  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      pMG_fine->Phi[2*k  ][2*j  ][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j  ][2*i  ] += 0.25*pMG_coarse->Phi[k+1][j+1][i+1];

      pMG_fine->Phi[2*k  ][2*j  ][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j  ][2*i-1] += 0.25*pMG_coarse->Phi[k+1][j+1][i-1];

      pMG_fine->Phi[2*k  ][2*j-1][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j-1][2*i  ] += 0.25*pMG_coarse->Phi[k+1][j-1][i+1];

      pMG_fine->Phi[2*k-1][2*j  ][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j  ][2*i  ] += 0.25*pMG_coarse->Phi[k-1][j+1][i+1];

      pMG_fine->Phi[2*k  ][2*j-1][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k  ][2*j-1][2*i-1] += 0.25*pMG_coarse->Phi[k+1][j-1][i-1];

      pMG_fine->Phi[2*k-1][2*j-1][2*i  ] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j-1][2*i  ] += 0.25*pMG_coarse->Phi[k-1][j-1][i+1];

      pMG_fine->Phi[2*k-1][2*j]  [2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j]  [2*i-1] += 0.25*pMG_coarse->Phi[k-1][j+1][i-1];

      pMG_fine->Phi[2*k-1][2*j-1][2*i-1] += 0.75*pMG_coarse->Phi[k  ][j  ][i  ];
      pMG_fine->Phi[2*k-1][2*j-1][2*i-1] += 0.25*pMG_coarse->Phi[k-1][j-1][i-1];
    }
  }}

  return;
}
#endif /* SELF_GRAVITY_USING_MULTIGRID */

/*----------------------------------------------------------------------------*/
/* set_iterate_bvals:  sets BC for Jacobi iterates for MPI parallel jobs.
 *   With self-gravity using multigrid, the boundary conditions at the edge of
 *   the Domain are held fixed.  So only ghostzones associated with internal
 *   boundaries between MPI grids need to be passed.
 *
 * This routine is largely a copy of set_bvals().
 * Order for updating boundary conditions must always be x1-x2-x3 in order to
 * fill the corner cells properly
 */

#ifdef SELF_GRAVITY_USING_MULTIGRID
#ifdef MPI_PARALLEL
void set_iterate_bvals(MGrid *pMG)
{
  int cnt3, cnt, err;
  MPI_Status stat;
  MPI_Request rq;

/*--- Step 1. ------------------------------------------------------------------
 * Boundary Conditions in x1-direction */

  cnt3 = pMG->Nx3 > 1 ? pMG->Nx3 + 1 : 1;
  cnt = nghost*(pMG->Nx2 + 1)*cnt3;

/* MPI blocks to both left and right */
  if (pMG->rx1_id >= 0 && pMG->lx1_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ox1(pMG,cnt,0,&rq);  /* send R */
    swap_iterate_ix1(pMG,cnt,1,&rq);  /* listen L */

    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ix1(pMG,cnt,0,&rq);  /* send L */
    swap_iterate_ox1(pMG,cnt,1,&rq);  /* listen R */
  }

/* Physical boundary on left, MPI block on right */
  if (pMG->rx1_id >= 0 && pMG->lx1_id < 0) {
    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ox1(pMG,cnt,0,&rq);  /* send R */
    swap_iterate_ox1(pMG,cnt,1,&rq);  /* listen R */
  }

/* MPI block on left, Physical boundary on right */
  if (pMG->rx1_id < 0 && pMG->lx1_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ix1(pMG,cnt,0,&rq);  /* send L */
    swap_iterate_ix1(pMG,cnt,1,&rq);  /* listen L */
  }

/*--- Step 2. ------------------------------------------------------------------
 * Boundary Conditions in x2-direction */

  cnt3 = pMG->Nx3 > 1 ? pMG->Nx3 + 1 : 1;
  cnt = nghost*(pMG->Nx1 + 2*nghost)*cnt3;

/* MPI blocks to both left and right */
  if (pMG->rx2_id >= 0 && pMG->lx2_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ox2(pMG,cnt,0,&rq);  /* send R */
    swap_iterate_ix2(pMG,cnt,1,&rq);  /* listen L */

    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ix2(pMG,cnt,0,&rq);  /* send L */
    swap_iterate_ox2(pMG,cnt,1,&rq);  /* listen R */
  }

/* Physical boundary on left, MPI block on right */
  if (pMG->rx2_id >= 0 && pMG->lx2_id < 0) {
    /* Post a non-blocking receive for the input data from the right grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ox2(pMG,cnt,0,&rq);  /* send R */
    swap_iterate_ox2(pMG,cnt,1,&rq);  /* listen R */
  }

/* MPI block on left, Physical boundary on right */
  if (pMG->rx2_id < 0 && pMG->lx2_id >= 0) {
    /* Post a non-blocking receive for the input data from the left grid */
    err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
      boundary_cells_tag, MPI_COMM_WORLD, &rq);
    if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

    swap_iterate_ix2(pMG,cnt,0,&rq);  /* send L */
    swap_iterate_ix2(pMG,cnt,1,&rq);  /* listen L */
  }

/*--- Step 3. ------------------------------------------------------------------
 * Boundary Conditions in x3-direction */

  if (pMG->Nx3 > 1){

    cnt = nghost*(pMG->Nx1 + 2*nghost)*(pMG->Nx2 + 2*nghost);

/* MPI blocks to both left and right */
    if (pMG->rx3_id >= 0 && pMG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      swap_iterate_ox3(pMG,cnt,0,&rq);  /* send R */
      swap_iterate_ix3(pMG,cnt,1,&rq);  /* listen L */

      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      swap_iterate_ix3(pMG,cnt,0,&rq);  /* send L */
      swap_iterate_ox3(pMG,cnt,1,&rq);  /* listen R */
    }

/* Physical boundary on left, MPI block on right */
    if (pMG->rx3_id >= 0 && pMG->lx3_id < 0) {
      /* Post a non-blocking receive for the input data from the right grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      swap_iterate_ox3(pMG,cnt,0,&rq);  /* send R */
      swap_iterate_ox3(pMG,cnt,1,&rq);  /* listen R */
    }

/* MPI block on left, Physical boundary on right */
    if (pMG->rx3_id < 0 && pMG->lx3_id >= 0) {
      /* Post a non-blocking receive for the input data from the left grid */
      err = MPI_Irecv(recv_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
		      boundary_cells_tag, MPI_COMM_WORLD, &rq);
      if(err) ath_error("[set_bvals]: MPI_Irecv error = %d\n",err);

      swap_iterate_ix3(pMG,cnt,0,&rq);  /* send L */
      swap_iterate_ix3(pMG,cnt,1,&rq);  /* listen L */
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Inner x1 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix1 and set_bvals/receive_ix1
 */

void swap_iterate_ix1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  jl = pMG->js;
  ju = pMG->je + 1;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke + 1;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack iterate into send buffer */

  if (swap_flag == 0) {
    il = pMG->is;
    iu = pMG->is + nghost - 1;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on L-x1 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx1_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ix1]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    il = pMG->is - nghost;
    iu = pMG->is - 1;

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ix1]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Outer x1 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ox1 and set_bvals/receive_ox1
 */

void swap_ox1(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  jl = pMG->js;
  ju = pMG->je + 1;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke + 1;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack iterate into send buffer */

  if (swap_flag == 0) {
    il = pMG->ie - nghost + 1;
    iu = pMG->ie;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on R-x1 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx1_id,
      boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ox1]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    il = pMG->ie + 1;
    iu = pMG->ie + nghost;

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ox1]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Inner x2 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix2 and set_bvals/receive_ix2
 */

void swap_iterate_ix2(MGrid *pG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  il = pMG->is - nghost;
  iu = pMG->ie + nghost;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke + 1;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack iterate into send buffer */

  if (swap_flag == 0) {
    jl = pMG->js;
    ju = pMG->js + nghost - 1;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on L-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx2_id,
       boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ix2]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    jl = pMG->js - nghost;
    ju = pMG->js - 1;

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ix2]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Outer x2 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix2 and set_bvals/receive_ix2
 */

void swap_iterate_ox2(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  il = pMG->is - nghost;
  iu = pMG->ie + nghost;

  if(pMG->Nx3 > 1){
    kl = pMG->ks;
    ku = pMG->ke + 1;
  } else {
    kl = ku = pMG->ks;
  }

/* Pack iterate into send buffer */

  if (swap_flag == 0) {
    jl = pMG->je - nghost + 1;
    ju = pMG->je;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on R-x2 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx2_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ox2]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    jl = pMG->je + 1;
    ju = pMG->je + nghost;

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ox2]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Inner x3 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix3 and set_bvals/receive_ix3
 */

void swap_iterate_ix3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  il = pMG->is - nghost;
  iu = pMG->ie + nghost;
  jl = pMG->js - nghost;
  ju = pMG->je + nghost;

/* Pack iterate into send buffer */

  if (swap_flag == 0) {
    kl = pMG->ks;
    ku = pMG->ks + nghost - 1;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on L-x3 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->lx3_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ix3]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    kl = pMG->ks - nghost;
    ku = pMG->ks - 1;

    /* Wait to receive the input data from the left grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ix3]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/* MPI_SWAP of boundary conditions, Outer x3 boundary
 *   This function does either a send (swap_flag=0), or receive (swap_flag=1)
 *   Largely copied from set_bvals/send_ix1 and set_bvals/receive_ix1
 */

void swap_iterate_ox3(MGrid *pMG, int cnt, int swap_flag, MPI_Request *prq)
{
  int i,il,iu,j,jl,ju,k,kl,ku,err;
  double *pd = send_buf;

  il = pMG->is - nghost;
  iu = pMG->ie + nghost;
  jl = pMG->js - nghost;
  ju = pMG->je + nghost;

/* Pack data in Gas structure into send buffer */

  if (swap_flag == 0) {
    kl = pMG->ke - nghost + 1;
    ku = pMG->ke;
    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          *(pd++) = pMG->Phi[k][j][i];
        }
      }
    }

    /* send contents of buffer to the neighboring grid on R-x3 */
    err = MPI_Send(send_buf, cnt, MPI_DOUBLE, pMG->rx3_id,
                   boundary_cells_tag, MPI_COMM_WORLD);
    if(err) ath_error("[swap_iterate_ox3]: MPI_Send error = %d\n",err);
  }

/* Receive message and unpack iterate */

  if (swap_flag == 1) {
    kl = pMG->ke + 1;
    ku = pMG->ke + nghost;

    /* Wait to receive the input data from the right grid */
    err = MPI_Wait(prq, &stat);
    if(err) ath_error("[swap_iterate_ox3]: MPI_Wait error = %d\n",err);

    for (k=kl; k<=ku; k++){
      for (j=jl; j<=ju; j++){
        for (i=il; i<=iu; i++){
          pMG->Phi[k][j][i] = *(pd++);
        }
      }
    }
  }

  return;
}
#endif /* MPI_PARALLEL */
#endif /* SELF_GRAVITY_USING_MULTIGRID */
