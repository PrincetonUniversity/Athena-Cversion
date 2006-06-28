#ifndef ATHENA_H
#define ATHENA_H 
/*==============================================================================
 * FILE: athena.h
 *
 * PURPOSE: Contains definitions of the following data types and structures:
 *   Real   - either float or double, depending on configure option
 *   Gas    - cell-centered conserved variables: d,M1,M2,M3,[E],[B1c,B2c,B3c] 
 *   Prim   - cell-centered primitive variables: d,V1,V2,V3,[P],[B1c,B2c,B3c]
 *   Cons1D - conserved variables in 1D: d,Mx,My,Mz,[E],[By,Bz]
 *   Prim1D - primitive variables in 1D: d,Vx,Vy,Vz,[E],[By,Bz]
 *   Grid   - arrays of Gas and face-centered B, coordinates, time, etc.
 *   Output
 * Also contains the following global variables:
 *   CourNo, Iso_csound, Iso_csound2, Gamma, Gamma_1, Gamma_2;
 *============================================================================*/
#include "defs.h"

/*----------------------------------------------------------------------------*/
/* variable type Real
 *   depends on macro set by configure
 */ 

#if defined(SINGLE_PREC)
typedef float  Real;
#elif defined(DOUBLE_PREC)
typedef double Real;
#else
# error "Not a valid precision flag"
#endif

/*----------------------------------------------------------------------------*/
/* Some global variables
 */

extern Real CourNo; /* The Courant, Friedrichs, & Lewy (CFL) Number */

#ifdef ISOTHERMAL
extern Real Iso_csound,Iso_csound2; /* Isothermal sound speed, and its square */
#else
extern Real Gamma, Gamma_1, Gamma_2; /* adiabatic index, and g-1, g-2 */
#endif

/*----------------------------------------------------------------------------*/
/* structure Gas: 
 *  IMPORTANT!! The order of the elements in Gas CANNOT be changed.
 */

typedef struct{
  Real d;			/* density */
  Real M1;			/* Momenta in 1,2,3.  Use 1,2,3 to label */
  Real M2;                      /* directions in anticipation of         */
  Real M3;                      /* covariant coordinate in future        */
#ifndef ISOTHERMAL
  Real E;			/* Total energy density */
#endif /* ISOTHERMAL */
#ifdef MHD
  Real B1c;			/* cell centered magnetic fields in 1,2,3 */
  Real B2c;
  Real B3c;
#endif /* MHD */
}Gas;

/*----------------------------------------------------------------------------*/
/* structure Prim:  
 *  IMPORTANT!! The order of the elements in Prim CANNOT be changed.
 */

typedef struct{
  Real d;			/* density  */
  Real V1;			/* Velocity in 1,2,3 */
  Real V2;
  Real V3;
#ifndef ISOTHERMAL
  Real P;			/* pressure */
#endif /* ISOTHERMAL */
#ifdef MHD
  Real B1c;                     /* cell centered magnetic fields in 1,2,3 */
  Real B2c;
  Real B3c;
#endif /* MHD */
}Prim;

/*----------------------------------------------------------------------------*/
/* structure Cons1D:  Note 1D vectors do not contain Bx.
 *  IMPORTANT!! The order of the elements in Cons1D CANNOT be changed.
 */

typedef struct{
  Real d;			/* density */
  Real Mx;			/* Momenta in X,Y,Z.  Use X,Y,Z now instead */
  Real My;                      /* 1,2,3 since this structure can contain a */
  Real Mz;                      /* slice in any dimension, 1,2,or 3 */
#ifndef ISOTHERMAL
  Real E;			/* Total energy density */
#endif /* ISOTHERMAL */
#ifdef MHD
  Real By;			/* cell centered magnetic fields in Y,Z */
  Real Bz;
#endif /* MHD */
}Cons1D;

/*----------------------------------------------------------------------------*/
/* structure Prim1D:  Note 1D vectors do not contain Bx.
 *  IMPORTANT!! The order of the elements in Prim1D CANNOT be changed.
 */

typedef struct{
  Real d;			/* density */
  Real Vx;			/* Velocity in X,Y,Z */
  Real Vy;
  Real Vz;
#ifndef ISOTHERMAL
  Real P;			/* pressure */
#endif /* ISOTHERMAL */
#ifdef MHD
  Real By;			/* cell centered magnetic fields in Y,Z */
  Real Bz;
#endif /* MHD */
}Prim1D;

/*----------------------------------------------------------------------------*/
/* structure Grid:  
 *   contains array of Gas and face-centered B, coordinates, integration time
 *   and timestep, basename for outputs, and various other data associated
 *   with an individual grid block.  Initialized by init_grid_block().
 *   By using an array of Gas, rather than arrays of each variable, we
 *   guarantee all the variables for each cell are contiguous in memory.
 */

typedef struct{
  Gas ***U;			/* pointer to a 3D array of Gas'es */
#ifdef MHD
  Real ***B1i,***B2i,***B3i;    /* pointer to a 3D array of interface B's */
#endif /* MHD */
  Real x1_0;		        /* x1-position of coordinate ix = 0 */
  Real x2_0;		        /* x2-position of coordinate jx = 0 */
  Real x3_0;		        /* x3-position of coordinate kx = 0 */
  Real dx1,dx2,dx3;	       	/* cell size */
  Real dt,time;			/* time step, absolute time */
  int nstep;			/* number of integration steps taken */
  int Nx1,Nx2,Nx3;		/* number of zones in x1, x2, x3 direction */
  int is,ie;			/* start/end cell index in x1 direction */
  int js,je;			/* start/end cell index in x2 direction */
  int ks,ke;			/* start/end cell index in x3 direction */
  int idisp;                    /* coordinate ix = index i + idisp */
  int jdisp;                    /* coordinate jx = index j + jdisp */
  int kdisp;                    /* coordinate kx = index k + kdisp */
  char *outfilename;		/* basename for output files */
#ifdef MPI_PARALLEL
  int my_id;                    /* In MPI this is also called the rank */
  int nproc;                    /* Number of processes in the Calculation */
/* The following are the id's (in MPI the rank) of tasks which are
 *  working on neighboring sections of the complete grid.  These are
 *  intended solely for use by the boundary condition routines. In
 *  particular, for periodic boundary conditions, these id's may look
 *  like a circular linked list.
 */
  int rx1_id, lx1_id; /* Right- and Left-x1 neighbor grid's id */
  int rx2_id, lx2_id; /* Right- and Left-x2 neighbor grid's id */
  int rx3_id, lx3_id; /* Right- and Left-x3 neighbor grid's id */
#endif
}Grid;

typedef struct Domain_s{
  int ixs, jxs, kxs; /* Minimum coordinate of computational grid cells */
  int ixe, jxe, kxe; /* Maximum coordinate of computational grid cells */
  int my_id; /* The task ID in PVM or the rank in MPI */
}Domain;


struct Output_s;

typedef void (*VGFunout_t)(Grid *pGrid, struct Output_s *pout);

typedef Real (*Gasfun_t)(const Grid *pG, const int i, const int j, const int k);


/*----------------------------------------------------------------------------*/
/* structure Output: */
  
typedef struct Output_s{

/* variables which describe output time/number/file */
  int n;          /* which N (from outN) it is in the <output> block */
  Real dt;        /* diag_dt: */
  Real t;         /* next time to output */
  int num;        /* dump number (0=first) */
  char *out;      /* out: "fun(var)[dmin:dmax]{x=,y=,z=,amode}" */
  char *id;       /* filename is of the form <basename>[.idump][.id].<ext> */

/* variables which describe data min/max */
  Real dmin;      /* user defined min for scaling data */
  Real dmax;      /* user defined max for scaling data */
  Real gmin;      /* computed global min (over all data output so far) */
  Real gmax;      /* computed global max (over all data output so far) */
  int sdmin;      /* 0 = auto scale, otherwise use dmin */
  int sdmax;      /* 0 = auto scale, otherwise use dmax */

/* variables which describe coordinates of output data volume */
  int ndim;       /* 3=cube 2=slice 1=vector 0=scalar */
  int ix1l, ix1u; /* lower/upper x1 indices for data slice  -1 = all data */
  int ix2l, ix2u; /* lower/upper x2 indices for data slice  -1 = all data */
  int ix3l, ix3u; /* lower/upper x3 indices for data slice  -1 = all data */
  int Nx1;        /* number of grid points to be output in x1 */
  int Nx2;        /* number of grid points to be output in x2 */
  int Nx3;        /* number of grid points to be output in x3 */
  Real x1_0, dx1; /* origin and grid spacing of output slice in x1 */
  Real x2_0, dx2; /* origin and grid spacing of output slice in x2 */
  Real x3_0, dx3; /* origin and grid spacing of output slice in x3 */

/* variables which describe output format */
  char *out_fmt;  /* output format = {bin, tab, hdf, hst, pgm, ppm, ...} */
  char *dat_fmt;  /* format string for tabular type output, e.g. "%10.5e" */
  char *palette;  /* name of palette for RGB conversions */
  float *rgb;     /* array of RGB[256*3] values derived from palette */
  float *der;     /* helper array of derivatives for interpolation into RGB */

/* pointers to output functions; data expressions */
  VGFunout_t fun; /* function pointer */
  Gasfun_t expr;  /* expression for output */
} Output;


typedef void (*VGFun_t)(Grid *pG); /* VGFun_t -> void grid function type */

/* Primitive Source Function Prototype */
typedef void (*PrimSrcFun_t)(const Grid *pG, const Prim *pW,
			     const int i, const int j, const int k, Prim *psrc);

/* Conservative Source Function Prototype */
typedef void (*ConsSrcFun_t)(const Grid *pG, const Gas *pU,
			     const int i, const int j, const int k, Gas *psrc);

/* Conservative Potential Function Prototype -- An example of this
   type of potential is a time independent gravitational potential. */
typedef Real (*ConsPotFun_t)(const Real x1, const Real x2, const Real x3);

/* Directions for the set_bvals_fun() function */
enum Direction {left_x1, right_x1, left_x2, right_x2, left_x3, right_x3};

#include "prototypes.h"

#endif /* ATHENA_H */
