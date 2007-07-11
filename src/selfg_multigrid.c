#include "copyright.h"
/*==============================================================================
 * FILE: selfg_multigrid.c
 *
 * PURPOSE: Contains functions to solve Poisson's equation for self-gravity in
 *   2D and 3D using multigrid.
 *
 *   These functions work for non-periodic domains.  A low-order multipole
 *   expansion is used to compute the potential on the boundaries.
 *
 * HISTORY:
 *   june-2007 - 2D and 3D solvers written by Irene Balmes
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

float G=6.67;

void selfg_by_multig_2d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY
  int i, is = pG->is, ie = pG->ie, im = (is+ie)/2;
  int j, js = pG->js, je = pG->je, jm = (js+je)/2;
  int k, ks = pG->ks, ke = pG->ke;
  float mass, error, r, solution[pG->Nx3+2*nghost][pG->Nx2+2*nghost][pG->Nx1+2*nghost];

/* Copy current potential into old */

  for (j=js-nghost; j<=je+nghost; j++){
    for (i=is-nghost; i<=ie+nghost; i++){
      pG->Phi_old[ks][j][i] = pG->Phi[ks][j][i];
      pG->Phi[ks][j][i] = 0;
    }
  }

/* Apply boundary conditions : first order */

  for (j=js; j<=je; j++){
    for (i=is; i<=ie; i++){
      mass += pG->U[ks][j][i].d;
    }
  }

  mass = mass*pG->dx1*pG->dx2;

  for (j=js-nghost; j<=je+nghost; j++)
    for (i=is-nghost; i<=ie+nghost; i++){
      r = pow(pow((i-im)*pG->dx1, 2) + pow((j-jm)*pG->dx2, 2), .5);
      if(r<=.3)
        solution[ks][j][i] = G*PI*(pow(r,2)-pow(.3,2))/2 + G*mass*log(.3);
      else
        solution[ks][j][i] = G*mass*log(r);
    }

  for (j=js-nghost; j<js; j++)
    for (i=is-nghost; i<=ie+nghost; i++){
      r = pow(pow((i-im)*pG->dx1, 2) + pow((j-jm)*pG->dx2, 2), .5);
      pG->Phi[ks][j][i] = G*mass*log(r);
    }

  for (j=js; j<=je; j++){
    for (i=is-nghost; i<is; i++){
      r = pow(pow((i-im)*pG->dx1, 2) + pow((j-jm)*pG->dx2, 2), .5);
      pG->Phi[ks][j][i] = G*mass*log(r);
    }
    for (i=ie+1; i<=ie+nghost; i++){
      r = pow(pow((i-im)*pG->dx1, 2) + pow((j-jm)*pG->dx2, 2), .5);
      pG->Phi[ks][j][i] = G*mass*log(r);
    }
  }

  for (j=je+1; j<=je+nghost; j++)
    for (i=is-nghost; i<=ie+nghost; i++){
      r = pow(pow((i-im)*pG->dx1, 2) + pow((j-jm)*pG->dx2, 2), .5);
      pG->Phi[ks][j][i] = G*mass*log(r);
    }

  for (j=js; j<=je; j++)
    for (i=is; i<=ie; i++)
      pG->Phi[ks][j][i] = -10;

/* Compute new potential */

  multi_2d(pG);

  for (j=js-nghost; j<=je+nghost; j++)
    for (i=is-nghost; i<=ie+nghost; i++){
//       printf("%f\t%f\t%f\n", pG->Phi[ks][j][i], solution[ks][j][i], pG->Phi[ks][j][i]-solution[ks][j][i]);
      error += fabs(pG->Phi[ks][j][i]-solution[ks][j][i]);
    }

  printf("error is %f\n", error/pow(256,2));

#endif
  return;
}
/* Functions for the 2d solver */

#ifdef SELF_GRAVITY
void multi_2d(Grid *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Grid pG_coarse;
  Real ***fine_Phi;

  if (pG->Nx1==4 || pG->Nx2==4)
    Jacobi_2d(pG);

  else{
    Jacobi_2d(pG);

    pG_coarse.U= (Gas ***) calloc_3d_array(pG->Nx3+2*nghost, pG->Nx2/2+2*nghost, pG->Nx1/2+2*nghost, sizeof(Gas));

    pG_coarse.Nx1 = pG->Nx1/2;
    pG_coarse.Nx2 = pG->Nx2/2;
    pG_coarse.Nx3 = pG->Nx3;

    pG_coarse.is = is;
    pG_coarse.js = js;
    pG_coarse.ks = ks;
    pG_coarse.ie = is+pG_coarse.Nx1-1;
    pG_coarse.je = js+pG_coarse.Nx2-1;

    pG_coarse.dx1 = pG->dx1*2;
    pG_coarse.dx2 = pG->dx2*2;

    Restriction_2d(pG, pG_coarse);

    pG_coarse.Phi=(Real ***) calloc_3d_array(pG_coarse.Nx3+2*nghost, pG_coarse.Nx2+2*nghost, pG_coarse.Nx1+2*nghost, sizeof(Real));

    multi_2d(&pG_coarse);

    fine_Phi=(Real ***) calloc_3d_array(pG->Nx3+2*nghost, pG->Nx2+2*nghost, pG->Nx1+2*nghost, sizeof(Real));

//     if (pG->Nx1==256)
//       for (j=pG_coarse.js; j<=pG_coarse.je; j++)
//         for (i=pG_coarse.is; i<=pG_coarse.ie; i++)
//           printf("%f\n", pG_coarse.Phi[ks][j][i]);

    Prolongation_2d(&pG_coarse, fine_Phi);

    for (j=js; j<=je; j++)
      for (i=is; i<=ie; i++)
        pG->Phi[ks][j][i] = pG->Phi[ks][j][i] + fine_Phi[ks][j][i];

    Jacobi_2d(pG);

//     Restriction_2d(pG, pG_coarse);
// 
//     multi_2d(&pG_coarse);
// 
//     Prolongation_2d(&pG_coarse, fine_Phi);
// 
//     for (j=js; j<=je; j++)
//       for (i=is; i<=ie; i++)
//         pG->Phi[ks][j][i] = pG->Phi[ks][j][i] + fine_Phi[ks][j][i];
// 
//     Jacobi_2d(pG);
    }

  return;
}

void Jacobi_2d(Grid *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  int n;
  Real temp[pG->Nx3+2*nghost][pG->Nx2+2*nghost][pG->Nx1+2*nghost];

  for (n=0; n<10; n++){
    for (j=js; j<=je; j++)
      for (i=is; i<=ie; i++){
        temp[ks][j][i] = 4*PI*G*pG->U[ks][j][i].d - (pG->Phi[ks][j][i-1] + pG->Phi[ks][j][i+1]) /pow(pG->dx1,2) - (pG->Phi[ks][j+1][i] + pG->Phi[ks][j-1][i]) /pow(pG->dx2,2);
        pG->Phi[ks][j][i] = -temp[ks][j][i]/(2*(1/pow(pG->dx1,2) + 1/pow(pG->dx2,2)));
        }
    }

  return;
}

void Restriction_2d(Grid *pG, Grid pG_coarse)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real ***temp;

  temp= (Real ***) calloc_3d_array(pG->Nx3+2*nghost, pG->Nx2+2*nghost, pG->Nx1+2*nghost, sizeof(Real));
  for (j=js; j<=je; j++)
    for (i=is; i<=ie; i++)
      temp[ks][j][i] = 4*PI*G*pG->U[ks][j][i].d - (pG->Phi[ks][j][i+1] + pG->Phi[ks][j][i-1] - 2*pG->Phi[ks][j][i]) /pow(pG->dx1,2) - (pG->Phi[ks][j+1][i] + pG->Phi[ks][j-1][i] - 2*pG->Phi[ks][j][i]) /pow(pG->dx2,2);

  for (j=js; j<=pG_coarse.ie; j++)
    for (i=is; i<=pG_coarse.je; i++)
      (pG_coarse.U[ks][j][i]).d = ((temp[ks][2*j-nghost][2*i-nghost] + (temp[ks][2*j-nghost+1][2*i-nghost] + temp[ks][2*j-nghost-1][2*i-nghost] + temp[ks][2*j-nghost][2*i-nghost+1] + temp[ks][2*j-nghost][2*i-nghost-1]) /2) /3)/(4*PI*G);

  return;
}

void Prolongation_2d(Grid *pG, Real ***fine_grid)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int ks = pG->ks;

  for (j=js; j<=je; j++)
    for (i=is; i<=ie; i++){
      if (i%2==0 && j%2==0)
        fine_grid[ks][2*j-nghost][2*i-nghost] = pG->Phi[ks][j][i];
      else if (i%2==1 && j%2==0)
        fine_grid[ks][2*j-nghost][2*i-nghost] = (pG->Phi[ks][j][i+1] + pG->Phi[ks][j][i-1]) /2;
      else if (i%2==0 && j%2==1)
        fine_grid[ks][2*j-nghost][2*i-nghost] = (pG->Phi[ks][j+1][i] + pG->Phi[ks][j-1][i]) /2;
      else if (i%2==1 && j%2==1)
        fine_grid[ks][2*j-nghost][2*i-nghost] = (pG->Phi[ks][j-1][i-1] + pG->Phi[ks][j-1][i+1] + pG->Phi[ks][j+1][i-1]) /3;
      }

  return;
}

#endif

/*----------------------------------------------------------------------------*/
/* selfg_by_multig_3d: Do not use with periodic BCs, uses multipole expansion
 *   to compute potential at boundary
 */

void selfg_by_multig_3d(Grid *pG, Domain *pD)
{
#ifdef SELF_GRAVITY_USING_MULTIGRID
  int i, is = pG->is, ie = pG->ie, im = (ie+is)/2;
  int j, js = pG->js, je = pG->je, jm = (je+js)/2;
  int k, ks = pG->ks, ke = pG->ke, km = (ke+ks)/2;
  float mass = 0.0, dVol, r, error;
  Real Grav_const = four_pi_G/(4.0*PI);

/* Copy current potential into old */

  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
      }
    }
  }

/* Compute solution at boundaries using monopole expansion */

  dvol = pG->dx1*pG->dx2*pG->dx3;
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
        pG->Phi[k][j][is-i]= Grav_const*mass/rad;

        cc_pos(pG,ie+i,j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][j][ie+i]= Grav_const*mass/rad;
      }
    }
  }

/*  Inner and outer x2 boundaries */

  for (k=ks; k<=ke; k++){
    for (j=1; j<=nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        cc_pos(pG,i,js-j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][js-j][i]= Grav_const*mass/rad;

        cc_pos(pG,i,je+j,k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[k][je+j][i]= Grav_const*mass/rad;
      }
    }
  }

/*  Inner and outer x3 boundaries */

  for (k=1; k<=nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        cc_pos(pG,i,j,ks-k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[ks-k][j][i]= Grav_const*mass/rad;

        cc_pos(pG,i,j,ke+k,&x1,&x2,&x3);
        rad = sqrt(x1*x1 + x2*x2 + x3*x3);
        pG->Phi[ke+k][j][i]= Grav_const*mass/rad;
      }
    }
  }

/* Compute new potential.  Note multi_3d calls itself recursively. */

  multi_3d(pG);

#endif
  return;
}

/*----------------------------------------------------------------------------*/
/* Functions needed for the multigrid solver in 3D
 */
#ifdef SELF_GRAVITY

void multi_3d(Grid *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Grid pG_coarse;
  Real ***fine_Phi;

  if (pG->Nx1==4 || pG->Nx2==4 || pG->Nx3==4)

    Jacobi_3d(pG);

  else{ 
    Jacobi_3d(pG);

    pG_coarse.U= (Gas ***) calloc_3d_array(pG->Nx3/2+2*nghost, pG->Nx2/2+2*nghost, pG->Nx1/2+2*nghost, sizeof(Gas));

    pG_coarse.Nx1 = pG->Nx1/2;
    pG_coarse.Nx2 = pG->Nx2/2;
    pG_coarse.Nx3 = pG->Nx3/2;

    pG_coarse.is = is;
    pG_coarse.js = js;
    pG_coarse.ks = ks;
    pG_coarse.ie = is+pG_coarse.Nx1-1;
    pG_coarse.je = js+pG_coarse.Nx2-1;
    pG_coarse.ke = ks+pG_coarse.Nx3-1;

    pG_coarse.dx1 = pG->dx1*2;
    pG_coarse.dx2 = pG->dx2*2;
    pG_coarse.dx3 = pG->dx3*2;

    Restriction_3d(pG, pG_coarse);

    pG_coarse.Phi=(Real ***) calloc_3d_array(pG_coarse.Nx3+2*nghost, pG_coarse.Nx2+2*nghost, pG_coarse.Nx1+2*nghost, sizeof(Real));

    multi_3d(&pG_coarse);
    multi_3d(&pG_coarse);

//     for (k=ks; k<=pG_coarse.ke; k++)
//       for (j=js; j<=pG_coarse.je; j++)
//         for (i=is; i<=pG_coarse.ie; i++)
//           printf("%f\n", pG_coarse.Phi[k][j][i]);

    fine_Phi=(Real ***) calloc_3d_array(pG->Nx3+2*nghost, pG->Nx2+2*nghost, pG->Nx1+2*nghost, sizeof(Real));

    Prolongation_3d(&pG_coarse, fine_Phi);

//     for (k=ks; k<=ke; k++)
//       for (j=js; j<=je; j++)
//         for (i=is; i<=ie; i++)
//           printf("%f\n", fine_Phi[k][j][i]);

    for (k=ks; k<=ke; k++)
      for (j=js; j<=je; j++)
        for (i=is; i<=ie; i++)
          pG->Phi[k][j][i] = pG->Phi[k][j][i] + fine_Phi[k][j][i];

    Jacobi_3d(pG);

//     Restriction_3d(pG, pG_coarse);
// 
//     multi_3d(&pG_coarse);
// 
//     Prolongation_3d(&pG_coarse, fine_Phi);
// 
//     for (k=ks; k<=ke; k++)
//       for (j=js; j<=je; j++)
//         for (i=is; i<=ie; i++)
//           pG->Phi[k][j][i] = pG->Phi[k][j][i] + fine_Phi[k][j][i];
// 
//     Jacobi_3d(pG);
// 
//     Restriction_3d(pG, pG_coarse);
// 
//     multi_3d(&pG_coarse);
// 
//     Prolongation_3d(&pG_coarse, fine_Phi);
// 
//     for (k=ks; k<=ke; k++)
//       for (j=js; j<=je; j++)
//         for (i=is; i<=ie; i++)
//           pG->Phi[k][j][i] = pG->Phi[k][j][i] + fine_Phi[k][j][i];
// 
//     Jacobi_3d(pG);
    }

  return;
}

void Jacobi_3d(Grid *pG)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  int n;
  Real bla;
  Real temp[pG->Nx3+2*nghost][pG->Nx2+2*nghost][pG->Nx1+2*nghost];

  for (n=0; n<=10; n++){
    for (k=ks; k<=ke; k++)
      for (j=js; j<=je; j++)
        for (i=is; i<=ie; i++){
          temp[k][j][i] = four_pi_G*pG->U[k][j][i].d;
          temp -= (pG->Phi[k][j][i-1] + pG->Phi[k][j][i+1]) /pow(pG->dx1,2)
          temp -= (pG->Phi[k][j+1][i] + pG->Phi[k][j-1][i]) /pow(pG->dx2,2)
          temp -= (pG->Phi[k+1][j][i] + pG->Phi[k-1][j][i]) /pow(pG->dx3,2) ;
          pG->Phi[k][j][i] = -temp[k][j][i]/(2*(pow(pG->dx1,-2) + pow(pG->dx2,-2) + pow(pG->dx3,-2)));
          }
    }

  return;
}

void Restriction_3d(Grid *pG, Grid pG_coarse)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;
  Real ***temp;

  temp= (Real ***) calloc_3d_array(pG->Nx3+2*nghost, pG->Nx2+2*nghost, pG->Nx1+2*nghost, sizeof(Real));
  for (k=ks; k<=ke; k++)
    for (j=js; j<=je; j++)
      for (i=is; i<=ie; i++)
        temp[k][j][i] = 4*PI*G*pG->U[k][j][i].d - (pG->Phi[k][j][i+1] + pG->Phi[k][j][i-1] - 2*pG->Phi[k][j][i]) /pow(pG->dx1,2) - (pG->Phi[k][j+1][i] + pG->Phi[k][j-1][i] - 2*pG->Phi[k][j][i]) /pow(pG->dx2,2) - (pG->Phi[k+1][j][i] + pG->Phi[k-1][j][i] - 2*pG->Phi[k][j][i]) /pow(pG->dx3,2);

  for(k=ks; k<=pG_coarse.ke; k++)
    for (j=js; j<=pG_coarse.je; j++)
      for (i=is; i<=pG_coarse.ie; i++)
        (pG_coarse.U[k][j][i]).d = ((temp[2*k-nghost][2*j-nghost][2*i-nghost] + (temp[2*k-nghost+1][2*j-nghost][2*i-nghost] + temp[2*k-nghost-1][2*j-nghost][2*i-nghost] + temp[2*k-nghost][2*j-nghost+1][2*i-nghost] + temp[2*k-nghost][2*j-nghost-1][2*i-nghost] + temp[2*k-nghost][2*j-nghost][2*i-nghost+1] + temp[2*k-nghost][2*j-nghost][2*i-nghost-1]) /2) /4)/(4*PI*G);

  return;
}

void Prolongation_3d(Grid *pG, Real ***fine_grid)
{
  int i, is = pG->is, ie = pG->ie;
  int j, js = pG->js, je = pG->je;
  int k, ks = pG->ks, ke = pG->ke;

  for (k=ks; k<=ke; k++)
    for (j=js; j<=je; j++)
      for (i=is; i<=ie; i++){
        if (i%2==0 && j%2==0 && k%2==0)
          fine_grid[k-nghost][j-nghost][i-nghost] = pG->Phi[k/2][j/2][i/2];
        else if (i%2!=0 && j%2==0 && k%2==0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[k/2][j/2][(i-1)/2] + pG->Phi[k/2][j/2][(i+1)/2]) /2;
        else if (i%2==0 && j%2!=0 && k%2==0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[k/2][(j-1)/2][i/2] + pG->Phi[k/2][(j+1)/2][i/2]) /2;
        else if (i%2==0 && j%2==0 && k%2!=0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[(k-1)/2][j/2][i/2] + pG->Phi[(k+1)/2][j/2][i/2]) /2;
        else if (i%2!=0 && j%2!=0 && k%2==0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[k/2][(j-1)/2][(i-1)/2] + pG->Phi[k/2][(j-1)/2][(i+1)/2] + pG->Phi[k/2][(j+1)/2][(i-1)/2]) /3;
        else if (i%2!=0 && j%2==0 && k%2!=0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[(k-1)/2][j/2][(i-1)/2] + pG->Phi[(k-1)/2][j/2][(i+1)/2] + pG->Phi[(k+1)/2][j/2][(i-1)/2]) /3;
        else if (i%2==0 && j%2!=0 && k%2!=0)
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[(k-1)/2][(j-1)/2][i/2] + pG->Phi[(k+1)/2][(j-1)/2][i/2] + pG->Phi[(k-1)/2][(j+1)/2][i/2]) /3;
        else
          fine_grid[k-nghost][j-nghost][i-nghost] = (pG->Phi[(k-1)/2][(j-1)/2][(i-1)/2] + pG->Phi[(k+1)/2][(j-1)/2][(i-1)/2] + pG->Phi[(k-1)/2][(j+1)/2][(i-1)/2] + pG->Phi[(k-1)/2][(j-1)/2][(i+1)/2]) /4;
       }

  return;
}
#endif
