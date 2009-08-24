#include "copyright.h"
/*==============================================================================
 * FILE: jet.c
 *
 * PURPOSE: Sets up a jet on x1-L boundary (left edge)
 *
 *============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 * jet_iib() - sets BCs on L-x1 (left edge) of grid.
 *============================================================================*/

void jet_iib(Grid *pGrid);

/* Make radius of jet and jet variables static global so they can be accessed
   by BC functions */
static Real rjet;
static Prim1D Wjet;
static Cons1D Ujet;
static Real x1_mid,x2_mid, x3_mid;

#ifdef MHD
static Real bxjet;
#endif

/*----------------------------------------------------------------------------*/
/* problem */

void problem(Grid *pGrid, Domain *pDomain){
   int i, is = pGrid->is, ie = pGrid->ie;
   int j, js = pGrid->js, je = pGrid->je;
   int k, ks = pGrid->ks, ke = pGrid->ke;
   int il,iu,jl,ju,kl,ku;
   Prim1D W;
   Cons1D U;
   Real x1_min, x1_max;
   Real x2_min, x2_max;
   Real x3_min, x3_max;
#ifdef MHD
   Real bx;
#endif

/* read parameters from input file */

   W.d  = par_getd("problem", "d");
   W.P  = par_getd("problem", "p");
   W.Vx = par_getd("problem", "vx");
   W.Vy = par_getd("problem", "vy");
   W.Vz = par_getd("problem", "vz");
#ifdef MHD
   bx   = par_getd("problem", "bx");
   W.By = par_getd("problem", "by");
   W.Bz = par_getd("problem", "bz");
#endif

   Wjet.d  = par_getd("problem", "djet");
   Wjet.P  = par_getd("problem", "pjet");
   Wjet.Vx = par_getd("problem", "vxjet");
   Wjet.Vy = par_getd("problem", "vyjet");
   Wjet.Vz = par_getd("problem", "vzjet");
#ifdef MHD
   bxjet   = par_getd("problem", "bxjet");
   Wjet.By = par_getd("problem", "byjet");
   Wjet.Bz = par_getd("problem", "bzjet");
#endif

   rjet = par_getd("problem", "rjet");
   
   Prim1D_to_Cons1D(&U,&W MHDARG( , &bx));
   Prim1D_to_Cons1D(&Ujet,&Wjet MHDARG( , &bxjet));

   x1_min = par_getd("grid", "x1min");
   x1_max = par_getd("grid", "x1max");
   x2_min = par_getd("grid", "x2min");
   x2_max = par_getd("grid", "x2max");
   x3_min = par_getd("grid", "x3min");
   x3_max = par_getd("grid", "x3max");

   x1_mid = 0.5 * (x1_max + x1_min);
   x2_mid = 0.5 * (x2_max + x2_min);
   x3_mid = 0.5 * (x3_max + x3_min);

/* initialize index bounds assuming problem in xy plane */
   iu = ie + nghost;
   il = is - nghost;
   ju = je + nghost;
   jl = js - nghost;
   if(pGrid->Nx3 > 1){
      ku = pGrid->ke + nghost;
      kl = pGrid->ks - nghost;
   }
   else{
      ku = pGrid->ke;
      kl = pGrid->ks;
   }

/* initialize conserved variables */
   
   for(k=kl; k<=ku; k++){
      for(j=jl; j<=ju; j++){
         for(i=il; i<=iu; i++){
            pGrid->U[k][j][i].d  = U.d;
            pGrid->U[k][j][i].M1 = U.Mx;
            pGrid->U[k][j][i].M2 = U.My;
            pGrid->U[k][j][i].M3 = U.Mz;
            pGrid->U[k][j][i].E  = U.E;
#ifdef MHD
            pGrid->U[k][j][i].B1c = bx;
            pGrid->U[k][j][i].B2c = U.By;
            pGrid->U[k][j][i].B3c = U.Bz;
            pGrid->B1i[k][j][i] = bx;
            pGrid->B2i[k][j][i] = U.By;
            pGrid->B3i[k][j][i] = U.Bz;
#endif

#ifdef SPECIAL_RELATIVITY
            pGrid->W[k][j][i].d  = W.d;
            pGrid->W[k][j][i].V1 = W.Vx;
            pGrid->W[k][j][i].V2 = W.Vy;
            pGrid->W[k][j][i].V3 = W.Vz;
            pGrid->W[k][j][i].P  = W.P;
#ifdef MHD
            pGrid->W[k][j][i].B1c = bx;
            pGrid->W[k][j][i].B2c = W.By;
            pGrid->W[k][j][i].B3c = W.Bz;
#endif
#endif
         }
      }
   }

/* Set boundary value function pointers */

   set_bvals_mhd_fun(left_x1,jet_iib);

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

void problem_write_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, Domain *pD, FILE *fp)
{
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VGFunout_t get_usr_out_fun(const char *name){
  return NULL;
}

#ifdef PARTICLES
PropFun_t get_usr_par_prop(const char *name)
{
  return NULL;
}

void gasvshift(const Real x1, const Real x2, const Real x3, Real *u1, Real *u2, Real *u3)
{
  return;
}

void Userforce_particle(Vector *ft, const Real x1, const Real x2, const Real x3, Real *w1, Real *w2, Real *w3)
{
  return;
}
#endif

void Userwork_in_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}

void Userwork_after_loop(Grid *pGrid, Domain *pDomain)
{
  return;
}


/*===================== PRIVATE FUNCTIONS ====================================*/

/*******************************************************************************
 * FUNCTION: jet_iib: sets ghost zones to either outflow BC or
 * if cell is within jet radius, to jet values
 ******************************************************************************/

void jet_iib(Grid *pGrid){
   int i, is = pGrid->is;
   int j, js = pGrid->js, je = pGrid->je;
   int k, ks = pGrid->ks, ke = pGrid->ke;
   Real rad,x1,x2,x3;

   for(k=ks; k<=ke; k++){
      for(j=js; j<=je; j++){
         for(i=1; i<=nghost; i++){
            cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
            rad = sqrt(SQR(x2 - x2_mid) + SQR(x3 - x3_mid));
            
            if(rad <= rjet){
               pGrid->U[k][j][i].d  = Ujet.d;
               pGrid->U[k][j][i].M1 = Ujet.Mx;
               pGrid->U[k][j][i].M2 = Ujet.My;
               pGrid->U[k][j][i].M3 = Ujet.Mz;
               pGrid->U[k][j][i].E  = Ujet.E;
#ifdef MHD
               pGrid->U[k][j][i].B1c = bxjet;
               pGrid->U[k][j][i].B2c = Ujet.By;
               pGrid->U[k][j][i].B3c = Ujet.Bz;
#endif

#ifdef SPECIAL_RELATIVITY
               pGrid->W[k][j][i].d  = Wjet.d;
               pGrid->W[k][j][i].V1 = Wjet.Vx;
               pGrid->W[k][j][i].V2 = Wjet.Vy;
               pGrid->W[k][j][i].V3 = Wjet.Vz;
               pGrid->W[k][j][i].P  = Wjet.P;
#ifdef MHD
               pGrid->W[k][j][i].B1c = bxjet;
               pGrid->W[k][j][i].B2c = Wjet.By;
               pGrid->W[k][j][i].B3c = Wjet.Bz;
#endif
#endif
            }
            else{
               pGrid->U[k][j][i] = pGrid->U[k][j][is];
               pGrid->W[k][j][i] = pGrid->W[k][j][is];
            }
         }
      }
   }
}
