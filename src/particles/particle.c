#include "../copyright.h"
/*=============================================================================
FILE: particle.c
PURPOSE: This file contains the main particle integrator. The integrator is 2nd
  order and fully implicit, which has absolute stability and is capable of integrate
  particles of any size with any stopping time. The current version implements
  linear interpolation of fluid quantities (density, velocity, sound speed). The
  feedback of the particle drag to the gas is also calculated in the integrator.
  The additional feedback routine calculates the feedback for the gas's predictor
  step. Other routines include the resort the particles using the quick sort algorithm,
  and removal of particles in the ghost zone.

CONTAINS PUBLIC FUNCTIONS:
  void integrate_particle(Grid* pG);
  void init_particle(Grid *pG, Domain *pD);
  void particle_destruct(Grid *pG);
  void particle_realloc(Grid *pG, long n);
  void feedback(Grid* pG);
  void shuffle(Grid* pG);
  void remove_ghost_particle(Grid *pG);

History:
 Created:	Emmanuel Jacquet	May 2008
 Rewritten:	Xuening Bai		Feb. 2009

==============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "../globals.h"

#ifdef PARTICLES         /* endif at the end of the file */

/* Fluid information used for drag force and feedback calculation */
typedef struct Gas_Info_s{
  Real v1, v2, v3;	
#ifndef ISOTHERMAL	
  Real cs;
#endif /* ISOTHERMAL */
}Gas_Info;

/* Filewide global variables */
int il,iu, jl,ju, kl,ku;	/* left and right limit of grid indices */
Real x1lpar,x1upar,x2lpar,x2upar,x3lpar,x3upar;	/* left and right limit of grid boundary */
static Gas_Info ***mygas;	/* additional gas information */

static GVDFun_t gasvshift = NULL;/* the gas velocity difference from Keplerian due to pressure gradient */
Real vsc1, vsc2;		/* coefficients for velocity difference */


/*=========================== PROTOTYPES OF PRIVATE FUNCTIONS ===============================*/

void getwei_linear(Grid *pG, Real x1, Real x2, Real x3, Real dx11, Real dx21, Real dx31, Real weight[2][2][2], int *is, int *js, int *ks);
int getval_linear(Grid *pG, Real weight[2][2][2], int is, int js, int ks, Real *rho, Real *u1, Real *u2, Real *u3, Real *cs);

Real get_ts(Grid *pG, Grain *curG, Real rho, Real cs, Real vd);
void grid_limit(Grid *pG, Domain *pD);
void get_gasinfo(Grid *pG);

#ifdef FEEDBACK
void distF_linear(Grid *pG, Real weight[2][2][2], Real is, Real js, Real ks, Real fb1, Real fb2, Real fb3);
#endif

void gasvshift_quad(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3);

int compare_gr(Grid *pG, Real dx11, Real dx21, Real dx31, Grain gr1, Grain gr2);
void quicksort_particle(Grid *pG, Real dx11, Real dx21, Real dx31, long start, long end);

/*=========================== PUBLIC FUNCTIONS ===============================*/

/* --------------Main Particle Integrator (corrector step) -----------------------
   Evolve particles for one full step using 2nd order implicit mid-
   point method.
   Input: 
     pG: grid which is already evolved in the predictor step. The
         particles are unevolved.
   Output:
     vmax: maximum of the updated particle velocity.
     pG: particle velocity updated with one full time step.
         feedback force is also calculated for the corrector step.
   Note: numerical experiment shows that velocity evolution may be unstable if t_stop/dt<0.006.
         therefore, we record the stopping time.
*/
void integrate_particle(Grid *pG)
{
  /* loca variables */
  Grain *curG, *curP, mygr;	/* pointer of the current working position */
  Real weight[2][2][2];		/* weight function */
  int is, js, ks;		/* starting location of interpolation */
  long p;			/* particle index */
  Real rho, u1, u2, u3;		/* density and velocity of the fluid from interpolation */
  Real cs, tstop;		/* sound speek, stopping time */
  Real vd1, vd2, vd3, vd;	/* velocity difference between particle and fluid */
  Real dv1, dv2, dv3;		/* amount of velocity update */
  Real fc1, fc2, fc3, fp1, fp2, fp3;/* force at current and predictor position */
  Real ft1, ft2, ft3;		/* total force */
  Real dx11, dx21, dx31, cellvol1;/* one over dx1, dx2, dx3 and cell volumn */
  Real ts11, ts12;		/* 1/stopping time */
  Real b0, b1, b2, oh;		/* other shortcut expressions */
  Real x1n, x2n, x3n;		/* first order new position at half a time step */
  Real omg, omg2;		/* Omega, Omega^2 */
#ifdef FEEDBACK
  Real fb1, fb2, fb3;		/* feedback force */
#endif

  /*-------------------- Initialization --------------------*/
  get_gasinfo(pG);		/* calculate gas information */

#ifdef FEEDBACK
  for (ks=kl; ks<=ku; ks++)	/* refresh feedback array */
    for (js=jl; js<=ju; js++)
      for (is=il; is<=iu; is++) {
        pG->feedback[ks][js][is].x1 = 0.0;
        pG->feedback[ks][js][is].x2 = 0.0;
        pG->feedback[ks][js][is].x3 = 0.0;
      }
#endif /* FEEDBACK */

  curP = &(mygr);	/* temperory particle */

  /* calculate some shortcut expressions */
#ifdef SHEARING_BOX
  omg = Omega;         omg2 = SQR(omg);
#else
  omg = 0.0;		   omg2 = 0.0;
#endif

  cellvol1 = 1.0;
  /* dx11, dx21, dx31 are shortcut expressions as well as dimension indicators */
  if (pG->Nx1 > 1)  { dx11 = 1.0/pG->dx1; cellvol1 *= dx11; }
  else dx11 = 0.0;

  if (pG->Nx2 > 1)  { dx21 = 1.0/pG->dx2; cellvol1 *= dx21; }
  else dx21 = 0.0;

  if (pG->Nx3 > 1)  { dx31 = 1.0/pG->dx3; cellvol1 *= dx31; }
  else dx31 = 0.0;

  /*---------------Main loop over all the particles--------------*/

  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particles */
    curG = &(pG->particle[p]);

    /* step 1: predictor of the particle position after one time step */
    x1n = curG->x1+curG->v1*pG->dt;
    x2n = curG->x2+curG->v2*pG->dt;
    x3n = curG->x3+curG->v3*pG->dt;

    /* step 2: calculate the force at current position */
    /* interpolation to get fluid density, velocity and the sound speed */
    getwei_linear(pG, curG->x1, curG->x2, curG->x3, dx11, dx21, dx31, weight, &is, &js, &ks);
    if (getval_linear(pG, weight, is, js, ks, &rho, &u1, &u2, &u3, &cs) == 0)
    { /* particle in the grid */
      /* apply gas velocity shift due to pressure gradient */
      gasvshift(curG->x1, curG->x2, curG->x3, &u1, &u2, &u3);
      /* velocity difference */
      vd1 = curG->v1-u1;
      vd2 = curG->v2-u2;
      vd3 = curG->v3-u3;
      vd = sqrt(SQR(vd1) + SQR(vd2) + SQR(vd3)); /* dimension independent */

      /* particle stopping time */
      tstop = get_ts(pG, curG, rho, cs, vd);
      ts11 = 1.0/tstop;
    }
    else { /* particle out of the grid, free motion, with warning sign */
      vd1 = 0.0;        vd2 = 0.0;      vd3 = 0.0;      ts11 = 0.0;
      ath_perr(1, "Particle move out of grid %d!\n", pG->my_id); /* level = ? */
    }

    /* Drag force */
    fc1 = -ts11*vd1;
    fc2 = -ts11*vd2;
    fc3 = -ts11*vd3;

    /* Force due to rotation */
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
  #ifdef FARGO
    fc1 += 2.0*curG->v3*omg;
    fc3 += -0.5*curG->v1*omg;
  #else
    fc1 += 3.0*omg2*curG->x1 + 2.0*curG->v3*omg;
    fc3 += -2.0*curG->v1*omg;
  #endif /* FARGO */
  #ifdef VERTICAL_GRAVITY
    fc2 += -omg2*curG->x2;
  #endif /* VERTICAL_GRAVITY */
#endif /* SHEARING_BOX */

    /* step 3: calculate the force at the predicted positoin */
    /* interpolation to get fluid density, velocity and the sound speed */
    getwei_linear(pG, x1n, x2n, x3n, dx11, dx21, dx31, weight, &is, &js, &ks);
    if (getval_linear(pG, weight, is, js, ks, &rho, &u1, &u2, &u3, &cs) == 0)
    { /* particle in the grid */
      /* apply gas velocity shift due to pressure gradient */
      gasvshift(x1n, x2n, x3n, &u1, &u2, &u3);
      /* velocity difference */
      vd1 = curG->v1-u1;
      vd2 = curG->v2-u2;
      vd3 = curG->v3-u3;
      vd = sqrt(SQR(vd1) + SQR(vd2) + SQR(vd3)); /* dimension independent */

      /* particle stopping time */
      tstop = get_ts(pG, curG, rho, cs, vd);
      ts12 = 1.0/tstop;
    }
    else { /* particle out of the grid, free motion, with warning sign */
      vd1 = 0.0;	vd2 = 0.0;	vd3 = 0.0;	ts12 = 0.0;
      ath_perr(1, "Particle move out of grid %d!\n", pG->my_id); /* level = ? */
    }

    /* Drag force */
    fp1 = -ts12*vd1;
    fp2 = -ts12*vd2;
    fp3 = -ts12*vd3;

    /* Force due to rotation */
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
  #ifdef FARGO
    fp1 += 2.0*curG->v3*omg;
    fp3 += -0.5*curG->v1*omg;
  #else
    fp1 += 3.0*omg2*x1n + 2.0*curG->v3*omg;
    fp3 += -2.0*curG->v1*omg;
  #endif /* FARGO */
  #ifdef VERTICAL_GRAVITY
    fp2 += -omg2*x2n;
  #endif /* VERTICAL_GRAVITY */
#endif /* SHEARING_BOX */

    /* step 4: calculate the velocity update */
    /* shortcut expressions */
    b0 = 1.0+pG->dt*ts11;
    oh = omg*pG->dt;

    /* Total force */
    ft1 = 0.5*(fc1+b0*fp1-2.0*oh*fp3);
    ft2 = 0.5*(fc2+b0*fp2);
#ifdef FARGO
    ft3 = 0.5*(fc3+b0*fp3+0.5*oh*fp1);
#else
    ft3 = 0.5*(fc3+b0*fp3+2.0*oh*fp1);
#endif

    /* calculate the inverse matrix elements */
    b1 = 1.0+0.5*pG->dt*(ts11 + ts12 + pG->dt*ts11*ts12);
    b2 = pG->dt/(b1*b1);

    /* velocity update */
    dv1 = b2*(b1*ft1+2.0*oh*ft3);
    dv2 = b2*b1*ft2;
#ifdef FARGO
    dv3 = b2*(b1*ft3-0.5*oh*ft1);
#else
    dv3 = b2*(b1*ft3-2.0*oh*ft1);
#endif

    /* Step 5: particle update to curP */
    /* velocity update */
    curP->v1 = curG->v1 + dv1;
    curP->v2 = curG->v2 + dv2;
    curP->v3 = curG->v3 + dv3;

    /* position update */
    if (pG->Nx1 > 1)
      curP->x1 = curG->x1 + 0.5*pG->dt*(curG->v1+curP->v1);
    else /* do not move if this dimension collapses */
      curP->x1 = curG->x1;

    if (pG->Nx2 > 1)
      curP->x2 = curG->x2 + 0.5*pG->dt*(curG->v2+curP->v2);
    else /* do not move if this dimension collapses */
      curP->x2 = curG->x2;

    if (pG->Nx3 > 1)
      curP->x3 = curG->x3 + 0.5*pG->dt*(curG->v3+curP->v3);
    else /* do not move if this dimension collapses */
      curP->x3 = curG->x3;

#ifdef FARGO
    /* shift = -3/2 * Omega * x * dt */
    curG->shift = -0.75*omg*(curG->x1+curP->x1)*pG->dt; /* Question: WHAT IS THE SHIFT CONVENSION FOR MHD FARGO? */
#endif

    curP->property = curG->property;

    /* step 6: calculate feedback force to the gas */
#ifdef FEEDBACK
    /* Force other than drag force */
    fb1 = 0.0;	fb2 = 0.0;	fb3 = 0.0;
    x1n = 0.5*(curG->x1+curP->x1);
    x2n = 0.5*(curG->x2+curP->x2);
    x3n = 0.5*(curG->x3+curP->x3);
#ifdef SHEARING_BOX
/* Note in shearing sheet, (X,Y,Z)=(x1,x3,x2) */
  #ifdef FARGO
    fb1 += 1.0*omg*(curG->v3+curP->v3);
    fb3 += -0.25*omg*(curG->v1+curP->v1);
  #else
    fb1 += 3.0*omg2*x1n + 1.0*omg*(curG->v3+curP->v3);
    fb3 += -1.0*omg*(curG->v1+curP->v1);
  #endif /* FARGO */
  #ifdef VERTICAL_GRAVITY
    fb2 += -omg2*x2n;
  #endif /* VERTICAL_GRAVITY */
#endif /* SHEARING_BOX */

    /* Velocity change due to the gas drag */
    fb1 = dv1 - pG->dt*fb1;
    fb2 = dv2 - pG->dt*fb2;
    fb3 = dv3 - pG->dt*fb3;

    /* Drag force density */
    fb1 = pG->grproperty[curG->property].m * fb1 * cellvol1;
    fb2 = pG->grproperty[curG->property].m * fb2 * cellvol1;
    fb3 = pG->grproperty[curG->property].m * fb3 * cellvol1;

    /* distribute the drag force (density) to the grid */
    getwei_linear(pG, x1n, x2n, x3n, dx11, dx21, dx31, weight, &is, &js, &ks);
    distF_linear(pG, weight, is, js, ks, fb1, fb2, fb3);
#endif /* FEEDBACK */

    /* step 7: update particle in pG */
    curG->x1 = curP->x1;
    curG->x2 = curP->x2;
    curG->x3 = curP->x3;
    curG->v1 = curP->v1;
    curG->v2 = curP->v2;
    curG->v3 = curP->v3;

  }/* end of the for loop */

  /* remove all ghost particles */
  remove_ghost_particle(pG);

  /* output the status */
  ath_pout(0, "In processor %d, there are %ld particles.\n", pG->my_id, pG->nparticle);

  return;
}

/* Initialization for particles.
   We assume to have "partypes" types of particles, each type has "parnum" particles.
   We enforce that each type has equal number of particles to ensure equal resolution.
   Allocate memory for the mygas array, feedback array.
*/
void init_particle(Grid *pG, Domain *pD)
{
  int N1T, N2T, N3T;
  Grain *GrArray;

  /* get coordinate limit */
  grid_limit(pG, pD);
  N1T = iu-il+1;
  N2T = ju-jl+1;
  N3T = ku-kl+1;

  /* check particle types */
  pG->partypes = par_geti("particle","partypes");

  if (pG->partypes < 0)
    ath_error("[init_particle]: Particle types must not be negative!\n");

  /* initialize the particle array */
  pG->arrsize = (long)(2*pG->partypes*par_geti("particle","parnum"));

  if (pG->arrsize < 0)
    ath_error("[init_particle]: Particle number must not be negative!\n");

  if (pG->arrsize<100) pG->arrsize = 100;

  pG->particle = (Grain*)calloc_1d_array(pG->arrsize, sizeof(Grain));
  if (pG->particle == NULL) goto on_error;

  /* allocate memory for particle properties */
  pG->grproperty = (Grain_Property*)calloc_1d_array(pG->partypes, sizeof(Grain_Property));
  if (pG->grproperty == NULL) goto on_error;

  grrhoa = (Real*)calloc_1d_array(pG->partypes, sizeof(Real));
  if (grrhoa == NULL) goto on_error;

  /* set gas velocity shift function pointer */
  gasvshift = gasvshift_quad;	/* by default will adopt the quadratic velocity shift */
  /* set velocity shift coefficients (for gas pressure gradient) */
  vsc1 = par_getd_def("problem","vsc1",0.0);
  vsc2 = par_getd_def("problem","vsc2",0.0);

  /* allocate the memory for gas and feedback arrays */
  mygas = (Gas_Info***)calloc_3d_array(N3T, N2T, N1T, sizeof(Gas_Info));
  if (mygas == NULL) goto on_error;

#ifdef FEEDBACK
  pG->feedback = (Vector***)calloc_3d_array(N3T, N2T, N1T, sizeof(Vector));
  if (pG->feedback == NULL) goto on_error;
#endif

  return;

  on_error:
    ath_error("[init_particle]: Error allocating memory.\n");
}

/* Finalization for particles */
void particle_destruct(Grid *pG)
{
  free_1d_array(pG->particle);

  free_1d_array(pG->grproperty);
  free_1d_array(grrhoa);

  /* free memory for gas and feedback arrays */
  if (mygas != NULL) free_3d_array(mygas);
#ifdef FEEDBACK
  if (pG->feedback != NULL) free_3d_array(pG->feedback);
#endif

  return;
}

/* Enlarge the particle array */
void particle_realloc(Grid *pG, long n)
{
  pG->arrsize = MAX((long)(1.2*pG->arrsize), n);
  if ((pG->particle = (Grain*)realloc(pG->particle, pG->arrsize*sizeof(Grain))) == NULL)
    ath_error("[init_particle]: Error re-allocating memory with array size %ld.\n", n);
}

/* Calculate the feedback of the drag force from the particle to the fluid
   Input: pG: grid with particles
   Output: pG: the array of drag forces exerted by the particle is updated
   This subroutine is used ONLY in the predictor step.
*/
#ifdef FEEDBACK
void feedback(Grid* pG)
{
  int i,j,k, is,js,ks;
  long p;			/* particle index */
  Real weight[2][2][2];		/* weight function */
  Real rho, u1, u2, u3;		/* density and velocity of the fluid from interpolation */
  Real cs, tstop;		/* sound speed, stopping time */
  Real vd1, vd2, vd3, vd;	/* velocity difference between particle and fluid */
  Real f1, f2, f3;		/* feedback force */
  Real m, ts1;			/* grain mass, 1/tstop */
  Real cellvol1;		/* 1/(dx1*dx2*dx3) */
  Real dx11, dx21, dx31;	/* one over dx1, dx2, dx3 */
  Grain *cur;			/* pointer of the current working position */

  /* initialization */
  get_gasinfo(pG);		/* calculate gas information */

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (i=il; i<=iu; i++) {
        pG->feedback[k][j][i].x1 = 0.0;
        pG->feedback[k][j][i].x2 = 0.0;
        pG->feedback[k][j][i].x3 = 0.0;
      }

  /* convenient expressions */
  cellvol1 = 1.0;
  if (pG->Nx1 > 1)  { dx11 = 1.0/pG->dx1; cellvol1 *= dx11; }
  else  dx11 = 0.0;

  if (pG->Nx2 > 1)  { dx21 = 1.0/pG->dx2;  cellvol1 *= dx21; }
  else dx21 = 0.0;

  if (pG->Nx3 > 1)  { dx31 = 1.0/pG->dx3;   cellvol1 *= dx31; }
  else dx31 = 0.0;

  /* loop over all particles to calculate the drag force */
  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particle */
    cur = &(pG->particle[p]);
    /* interpolation to get fluid density and velocity */
    getwei_linear(pG, cur->x1, cur->x2, cur->x3, dx11, dx21, dx31, weight, &is, &js, &ks);
    if (getval_linear(pG, weight, is, js, ks, &rho, &u1, &u2, &u3, &cs) == 0)
    { /* particle is in the grid */
      /* apply gas velocity shift due to pressure gradient */
      gasvshift(cur->x1, cur->x2, cur->x3, &u1, &u2, &u3);
      /* velocity difference */
      vd1 = cur->v1-u1;
      vd2 = cur->v2-u2;
      vd3 = cur->v3-u3;
      vd = sqrt(vd1*vd1 + vd2*vd2 + vd3*vd3);

      /* calculate particle stopping time */
      tstop = MAX(get_ts(pG, cur, rho, cs, vd), pG->dt); /* to avoid the stiff dependence on tstop */
      ts1 = 1.0/tstop;

      /* Drag force density */
      m = pG->grproperty[cur->property].m;
      f1 = m * vd1 * ts1 * cellvol1;
      f2 = m * vd2 * ts1 * cellvol1;
      f3 = m * vd3 * ts1 * cellvol1;

      /* distribute the drag force (density) to the grid */
      distF_linear(pG, weight, is, js, ks, f1, f2, f3);
    }
  }/* end of the for loop */

  return;
}

/* Temporary routine to implement feedback to the grid
   This routine is just first order accurate, and is just used
   for test purposes.
*/
void appy_feedback(Grid *pG)
{
  int i,j,k;

  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (i=il; i<=iu; i++) {
        pG->U[k][j][i].M1 -= pG->feedback[k][j][i].x1;
        pG->U[k][j][i].M2 -= pG->feedback[k][j][i].x2;
        pG->U[k][j][i].M3 -= pG->feedback[k][j][i].x3;
      }
}

#endif /* FEEDBACK */

/* Shuffle the particles
   Input: pG: grid with particles;
   Output: pG: particles in the linked list are rearranged by the order of their
           locations that are consistent with grid cell storage.
*/
void shuffle(Grid *pG)
{
  Grain *cur, *allgr;
  Real dx11,dx21,dx31;

  if (pG->Nx1 > 1) dx11 = 1.0/pG->dx1;  else  dx11 = 0.0;
  if (pG->Nx2 > 1) dx21 = 1.0/pG->dx2;  else  dx21 = 0.0;
  if (pG->Nx3 > 1) dx31 = 1.0/pG->dx3;  else  dx31 = 0.0;

  /* output status */
  ath_pout(0, "Resorting particles...\n");
  /* sort the particles according to their positions */
  quicksort_particle(pG, dx11, dx21, dx31, 0, pG->nparticle-1);

  return;
}

/* remove all the particle in the ghost cells */
void remove_ghost_particle(Grid *pG)
{
  Grain *cur;
  long p;

  /* loop over all particle to remove ghosts */
  p = 0;
  while (p<pG->nparticle) {
    cur = &(pG->particle[p]);
    if ((cur->x1 >= x1upar) || (cur->x1 < x1lpar) || (cur->x2 >= x2upar) || (cur->x2 < x2lpar) || (cur->x3 >= x3upar) || (cur->x3 < x3lpar))
    { /* if the particle is outside the box */
      pG->nparticle -= 1;
      pG->grproperty[cur->property].num -= 1;
      pG->particle[p] = pG->particle[pG->nparticle];
    }
    else
      p += 1;
  }

  return;
}

/*=========================== PRIVATE FUNCTIONS ===============================*/

/* Linear interpolation step 1: get weight
   Input: pG: grid; x1,x2,x3: global coordinate; dx11,dx21,dx31: 1 over dx1,dx2,dx3
   Output: weight: weight function; is,js,ks: starting cell indices in the grid.
   Note: this interpolation works in any 1-3 dimensions.
*/

void getwei_linear(Grid *pG, Real x1, Real x2, Real x3, Real dx11, Real dx21, Real dx31, Real weight[2][2][2], int *is, int *js, int *ks)
{
  int i, j, k, i1, j1, k1;
  Real a, b, c;				/* grid coordinate for the position (x1,x2,x3) */
  Real wei1[2], wei2[2], wei3[2];	/* weight function in x1,x2,x3 directions */

  /* find cell locations and calculate 1D weight */
  /* x1 direction */
  if (dx11 > 0.0) {
    i = celli(pG, x1, dx11, &i1, &a);		/* x1 index */
    i1 = i+i1-1;	*is = i1;		/* starting x1 index */
    wei1[1] = a - i1 - 0.5;			/* one direction weight */
    wei1[0] = 1.0 - wei1[1];			/* 0: left; 1: right */
  }
  else { /* x1 dimension collapses */
    *is = pG->is;				/* starting x1 index */
    wei1[1] = 0.0;				/* one direction weight */
    wei1[0] = 1.0;				/* 0: left; 1: right */
  }

  /* x2 direction */
  if (dx21 > 0.0) {
    j = cellj(pG, x2, dx21, &j1, &b);		/* x2 index */
    j1 = j+j1-1;	*js = j1;		/* starting x2 index */
    wei2[1] = b - j1 - 0.5;			/* one direction weight */
    wei2[0] = 1.0 - wei2[1];			/* 0: left; 1: right */
  }
  else { /* x2 dimension collapses */
    *js = pG->js;				/* starting x2 index */
    wei2[1] = 0.0;				/* one direction weight */
    wei2[0] = 1.0;				/* 0: left; 1: right */
  }

  /* x3 direction */
  if (dx31 > 0.0) {
    k = cellk(pG, x3, dx31, &k1, &c);		/* x3 index */
    k1 = k+k1-1;	*ks = k1;		/* starting x3 index */
    wei3[1] = c - k1 - 0.5;			/* one direction weight */
    wei3[0] = 1.0 - wei3[1];			/* 0: left; 1: right */
  }
  else { /* x3 dimension collapses */
    *ks = pG->ks;				/* starting x3 index */
    wei3[1] = 0.0;				/* one direction weight */
    wei3[0] = 1.0;				/* 0: left; 1: right */
  }

  /* calculate 3D weight */
  for (k=0; k<2; k++)
    for (j=0; j<2; j++)
      for (i=0; i<2; i++)
        weight[k][j][i] = wei1[i] * wei2[j] * wei3[k];

  return;
}

/* Linear interpolation step 2: get interpolated value
   Input:
     pG: grid; weight: weight function;
     is,js,ks: starting cell indices in the grid.
   Output:
     interpolated values of density, velocity and sound speed of the fluid
   Return: 0: normal exit;  -1: particle lie out of the grid, cannot interpolate!
   Note: this interpolation works in any 1-3 dimensions.
*/
int getval_linear(Grid *pG, Real weight[2][2][2], int is, int js, int ks, Real *rho, Real *u1, Real *u2, Real *u3, Real *cs)
{
  int i, j, k, i1, j1, k1;
  Real D, v1, v2, v3;		/* density and velocity of the fluid */
#ifndef ISOTHERMAL
  Real C = 0.0;			/* fluid sound speed */
#endif
  Real totwei, totwei1;		/* total weight (in case of edge cells) */

  /* linear interpolation */
  D = 0.0; v1 = 0.0; v2 = 0.0; v3 = 0.0;
  totwei = 0.0;		totwei1 = 1.0;
  /* Interpolate density, velocity and sound speed */
  /* Note: in lower dimensions only wei[0] is non-zero, which ensures the validity */
  for (k=0; k<2; k++) {
    k1 = k+ks;
    if ((k1 <= ku) && (k1 >= kl)) {
      for (j=0; j<2; j++) {
        j1 = j+js;
        if ((j1 <= ju) && (j1 >= jl)) {
          for (i=0; i<2; i++) {
            i1 = i+is;
            if ((i1 <= iu) && (i1 >= il)) {
              D += weight[k][j][i] * pG->U[k1][j1][i1].d;
              v1 += weight[k][j][i] * mygas[k1][j1][i1].v1;
              v2 += weight[k][j][i] * mygas[k1][j1][i1].v2;
              v3 += weight[k][j][i] * mygas[k1][j1][i1].v3;
#ifndef ISOTHERMAL
              C += weight[k][j][i] * mygas[k1][j1][i1].cs;
#endif
              totwei += weight[k][j][i];
            }
          }
        }
      }
    }
  }
  if (totwei < TINY_NUMBER) {
    ath_perr(0, "[particle]: Particle lies out of the grid: (is, js, ks)=(%d, %d, %d).\n",is, js, ks);
    return -1;
  }

  totwei1 = 1.0/totwei;
  *rho = D*totwei1;
  *u1 = v1*totwei1;	*u2 = v2*totwei1;	*u3 = v3*totwei1;
#ifdef ISOTHERMAL
  *cs = Iso_csound;
#else
  *cs = C*totwei1;
#endif /* ISOTHERMAL */

  return 0;
}

/* Distribute the feedback force to grid cells
   Input: 
     pG: grid;   weight: weight function; 
     is,js,ks: starting cell indices in the grid.
     f1, f2, f3: feedback force from one particle.
   Output:
     pG: feedback array is updated.
*/
#ifdef FEEDBACK
void distF_linear(Grid *pG, Real weight[2][2][2], Real is, Real js, Real ks, Real fb1, Real fb2, Real fb3)
{
  int i,j,k,i1,j1,k1;
  /* distribute feedback force */
  for (k=0; k<2; k++) {
    k1 = k+ks;
    if ((k1 <= ku) && (k1 >= kl)) {
      for (j=0; j<2; j++) {
        j1 = j+js;
        if ((j1 <= ju) && (j1 >= jl)) {
          for (i=0; i<2; i++) {
            i1 = i+is;
            if ((i1 <= iu) && (i1 >= il)) {
              pG->feedback[k1][j1][i1].x1 += weight[k][j][i] * fb1;
              pG->feedback[k1][j1][i1].x2 += weight[k][j][i] * fb2;
              pG->feedback[k1][j1][i1].x3 += weight[k][j][i] * fb3;
            }
          }
        }
      }
    }
  }

  return;
}
#endif /* FEEDBACK */

/* Calculate the stopping time 
   The relavent scale to calculate is: 
   1. a/lambda_m == alam
   2. rho_s*a in normalized unit == rhoa
*/
Real get_ts(Grid *pG, Grain *cur, Real rho, Real cs, Real vd)
{
  Real tstop;		/* stopping time */
  int type;		/* type of the particle (to look for properties) */
  Real a, rhos;		/* primitive properties: size, solid density in cgs unit */
  Real alam;		/* a/lambda: particle size/gas mean free path */
  Real rhoa;		/* rhoa: rho_s*a in normalized unit (density, velocity, time) */
  Real Re, CD;		/* Reynolds number and drag coefficient */

  /* particle properties */
  type = cur->property;
  a = pG->grproperty[type].rad;
  rhos = pG->grproperty[type].rho;
  rhoa = grrhoa[type];

  /* calculate particle size/gas mean free path */
  alam = alamcoeff * a * rho;  /* alamcoeff is defined in global.h */

  /* calculate the stopping time */
  if (alam < 2.25) {		/* Epstein regime */
    tstop = rhoa/(rho*cs);
  }

  else {
    Re = 4.0*alam*vd/cs;	/* the Reynolds number */

    if (Re < 1.) CD = 24.0/Re;	/* Stokes regime */
    else if (Re < 800.0) CD = 24.0*exp(-0.6*log(Re));
    else CD = 0.44;

    tstop = rhoa/(rho*vd*CD);
  } /* endif */

  return tstop;
}

/* Calculate the left and right grid limit
   Input: pG: grid;
   Output: il,iu,jl,ju,kl,ku: grid limit indices;
           x1lpar,x1upar,x2lpar,x2upar,x3lpar,x3upar: grid boundary coordinates
*/
void grid_limit(Grid *pG, Domain *pD)
{
  int m1, m2, m3;	/* dimension flags */
  int my_iproc, my_jproc, my_kproc;

  if (pG->Nx1 > 1) m1 = 1;
  else m1 = 0;

  if (pG->Nx2 > 1) m2 = 1;
  else m2 = 0;

  if (pG->Nx3 > 1) m3 = 1;
  else m3 = 0;

  /* set left and right grid indices */
  il = pG->is - m1*nghost;
  iu = pG->ie + m1*nghost;

  jl = pG->js - m2*nghost;
  ju = pG->je + m2*nghost;

  kl = pG->ks - m3*nghost;
  ku = pG->ke + m3*nghost;

  /* set left and right boundary for removing particles */
  /* Note: for outflow B.C. (ibc=2), we only remove the particles in
   * the outermost layer of the ghost cells */
  get_myGridIndex(pD, pG->my_id, &my_iproc, &my_jproc, &my_kproc);

  if ((par_geti("grid","ibc_x1") == 2) && (my_iproc == 0))
    x1lpar = pG->x1_0 + (il+m1 + pG->idisp)*pG->dx1;
  else
    x1lpar = pG->x1_0 + (pG->is + pG->idisp)*pG->dx1;

  if ((par_geti("grid","obc_x1") == 2) && (my_iproc == pD->NGrid_x1-1))
    x1upar = pG->x1_0 + (iu + pG->idisp)*pG->dx1;
  else
    x1upar = pG->x1_0 + (pG->ie + 1 + pG->idisp)*pG->dx1;

  if ((par_geti("grid","ibc_x2") == 2) && (my_jproc == 0))
    x2lpar = pG->x2_0 + (jl+m2 + pG->jdisp)*pG->dx2;
  else
    x2lpar = pG->x2_0 + (pG->js + pG->jdisp)*pG->dx2;

  if ((par_geti("grid","obc_x2") == 2) && (my_jproc == pD->NGrid_x2-1))
    x2upar = pG->x2_0 + (ju + pG->jdisp)*pG->dx2;
  else
    x2upar = pG->x2_0 + (pG->je + 1 + pG->jdisp)*pG->dx2;

  if ((par_geti("grid","ibc_x3") == 2) && (my_kproc == 0))
    x3lpar = pG->x3_0 + (kl+m3 + pG->kdisp)*pG->dx3;
  else
    x3lpar = pG->x3_0 + (pG->ks + pG->kdisp)*pG->dx3;

  if ((par_geti("grid","obc_x3") == 2) && (my_kproc == pD->NGrid_x3-1))
    x3upar = pG->x3_0 + (ku + pG->kdisp)*pG->dx3;
  else
    x3upar = pG->x3_0 + (pG->ke + 1 + pG->kdisp)*pG->dx3;

  if (pG->Nx1 == 1) x1upar += MAX(0.1*fabs(pG->x1_0), 1.);
  if (pG->Nx2 == 1) x2upar += MAX(0.1*fabs(pG->x2_0), 1.);
  if (pG->Nx3 == 1) x3upar += MAX(0.1*fabs(pG->x3_0), 1.);

  return;
}

/* Calculate the gas information from conserved variables
   Input: pG: grid (now already evolved in the predictor step).
   Output: calculate 3D array mygas in the grid structure.
           Calculated are gas velocity and sound speed.
*/
void get_gasinfo(Grid *pG)
{
  int i,j,k;
  Real rho1;
#ifdef ADIABATIC
  Real P;
#endif

  /* get gas information */
  for (k=kl; k<=ku; k++)
    for (j=jl; j<=ju; j++)
      for (i=il; i<=iu; i++) {
        rho1 = 1.0/(pG->U[k][j][i].d);
        mygas[k][j][i].v1 = pG->U[k][j][i].M1 * rho1;
        mygas[k][j][i].v2 = pG->U[k][j][i].M2 * rho1;
        mygas[k][j][i].v3 = pG->U[k][j][i].M3 * rho1;
#ifndef ISOTHERMAL
  #ifdef ADIABATIC
        /* E = P/(gamma-1) + rho*v^2/2 + B^2/2 */
        P = pG->U[k][j][i].E - 0.5*pG->U[k][j][i].d*(SQR(mygas[k][j][i].v1) \
             + SQR(mygas[k][j][i].v2) + SQR(mygas[k][j][i].v3));
    #ifdef MHD
        P = P - 0.5*(SQR(pG->U[k][j][i].B1c)+SQR(pG->U[k][j][i].B2c)+SQR(pG->U[k][j][i].B3c));
    #endif /* MHD */
        P = MAX(Gamma_1*P, TINY_NUMBER);
        mygas[k][j][i].cs = sqrt(Gamma*P/pG->U[k][j][i].d);
  #else
        ath_error("[get_gasinfo] can not calculate the sound speed!\n");
  #endif /* ADIABATIC */
#endif /* ISOTHERMAL */
      }

  return;
}

/*----------------Subroutines for pressure gradient -----------------*/

/* Calculate the gas velocity difference to Keplerian velocity as a
   function of position and apply the shift to the velocity (u1,u2,u3).
   This routine assumes a quadratic dependence on the x2 (or Z in shearing box).
   vsc1, vsc2 are coefficients.
*/
void gasvshift_quad(Real x1, Real x2, Real x3, Real *u1, Real *u2, Real *u3)
{
  *u2 += vsc1+vsc2*SQR(x2);
  return;
}

/*---------------Subroutines for the shuffle algorithm---------------*/

/* Compare the order of two particles according to their positions in the grid
   Input: pG: grid; 
          dx11,dx21,dx31: 1/dx`,1/dx2,1/dx3, or 0 if that dimension collapses.
          gr1,gr2: pointers of the two particles to be compared.
   Output: pointer of the particle that should be put in front of the other.
*/
int compare_gr(Grid *pG, Real dx11, Real dx21, Real dx31, Grain gr1, Grain gr2)
{
  int i1,j1,k1, i2,j2,k2;

  k1 = (int)((gr1.x3 - pG->x3_0) * dx31);	/* x3 index of gr1 */
  k2 = (int)((gr2.x3 - pG->x3_0) * dx31);	/* x3 index of gr2 */
  if (k1 < k2) return 1;
  if (k1 > k2) return 2;

  j1 = (int)((gr1.x2 - pG->x2_0) * dx21);	/* x2 index of gr1 */
  j2 = (int)((gr2.x2 - pG->x2_0) * dx21);	/* x2 index of gr2 */
  if (j1 < j2) return 1;
  if (j1 > j2) return 2;

  i1 = (int)((gr1.x1 - pG->x1_0) * dx11);	/* x1 index of gr1 */
  i2 = (int)((gr2.x1 - pG->x1_0) * dx11);	/* x1 index of gr2 */
  if (i1 < i2) return 1;
  if (i1 > i2) return 2;

  /* if they have equal indices, arbitrarily choose gr1 */
  return 1;
}

/* Quick sort algorithm to shuffle the particles
   Input: pG, dx11,dx21,dx31: for compare_gr subroutine only. See above.
          head, rear: head and rear of the linked list.
                      They do not contain data, or equal the pivot in the recursion.
          length: length of the linked list (does not contain head or rear).
   Output: *head: linked list with shuffling finished.
*/
void quicksort_particle(Grid *pG, Real dx11, Real dx21, Real dx31, long start, long end)
{
  long i, pivot;
  Grain gr;
  if (end <= start) return;	/* automatically sorted already */

  /* location of the pivot at half chain length */
  pivot = (long)((start+end+1)/2);

  /* move the pivot to the start */
  gr = pG->particle[pivot];
  pG->particle[pivot] = pG->particle[start];
  pG->particle[start] = gr;
  /* initial configuration */
  pivot = start;
  i = start + 1;
  /* move the particles that are "smaller" than the pivot before it */
  while (i <= end) {
    if (compare_gr(pG, dx11, dx21, dx31, pG->particle[pivot], pG->particle[i]) == 2)
    {/* the ith particle is smaller, move it before the pivot */
      gr = pG->particle[pivot];
      pG->particle[pivot] = pG->particle[i];
      pG->particle[i] = pG->particle[pivot+1];
      pG->particle[pivot+1] = gr;
      pivot += 1;
    }
    i += 1;
  }

  /* recursively call this routine to complete sorting */
  quicksort_particle(pG, dx11, dx21, dx31, start, pivot-1);
  quicksort_particle(pG, dx11, dx21, dx31, pivot+1, end);

  return;
}

#endif /*PARTICLES*/
