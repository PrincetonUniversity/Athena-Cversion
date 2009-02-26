#ifndef ATHENA_H
#define ATHENA_H 
/*==============================================================================
 * FILE: athena.h
 *
 * PURPOSE: Contains definitions of the following data types and structures:
 *   Real  - either float or double, depending on configure option
 *   Gas   - cell-centered conserved variables
 *   Prim  - cell-centered primitive variables
 *   Cons1D - conserved variables in 1D: same as Gas minus Bx
 *   Prim1D - primitive variables in 1D: same as Prim minus Bx
 *   Grain  - basic properties of particles
 *   Ray,Rad_Ran2_State,Ray_Tree,Radpoint,Radplane - for ionizing rad transport
 *   Grid   - everything needed by a Grid: arrays of Gas, B, indices, time, etc.
 *   Grid_Indices - indices and ID of one Grid in a Domain
 *   Domain - info on array of Grids covering computational domain 
 *   Output - everything associated with an individual output: time, type, etc.
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

/* general 3-vector */
typedef struct Vector_s{
  Real f1, f2, f3;		
}Vector;

/*----------------------------------------------------------------------------*/
/* structure Gas: conserved variables 
 *  IMPORTANT!! The order of the elements in Gas CANNOT be changed.
 */

typedef struct Gas_s{
  Real d;			/* density */
  Real M1;			/* Momenta in 1,2,3.  Use 1,2,3 to label */
  Real M2;                      /* directions in anticipation of         */
  Real M3;                      /* covariant coordinates in future       */
#ifndef BAROTROPIC
  Real E;			/* Total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;			/* cell centered magnetic fields in 1,2,3 */
  Real B2c;
  Real B3c;
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];              /* passively advected scalars */
#endif
}Gas;

/*----------------------------------------------------------------------------*/
/* structure Prim: primitive variables, used with special relativity 
 *  IMPORTANT!! The order of the elements in Prim CANNOT be changed.
 */

typedef struct Prim_s{
  Real d;			/* density  */
  Real V1;			/* Velocity in 1,2,3 */
  Real V2;
  Real V3;
#ifndef BAROTROPIC
  Real P;			/* pressure */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;                     /* cell centered magnetic fields in 1,2,3 */
  Real B2c;
  Real B3c;
#endif /* MHD */
#if (NSCALARS > 0)
  Real r[NSCALARS];              /* density-normalized advected scalars */
#endif
}Prim;

/*----------------------------------------------------------------------------*/
/* structure Cons1D:  conserved variables in 1D (does not contain Bx)
 *  IMPORTANT!! The order of the elements in Cons1D CANNOT be changed.
 */

typedef struct Cons1D_s{
  Real d;			/* density */
  Real Mx;			/* Momenta in X,Y,Z.  Use X,Y,Z now instead  */
  Real My;                      /* of 1,2,3 since this structure can contain */
  Real Mz;                      /* a slice in any dimension: 1,2,or 3        */
#ifndef BAROTROPIC
  Real E;			/* Total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real By;			/* cell centered magnetic fields in Y,Z */
  Real Bz;
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];              /* passively advected scalars */
#endif
}Cons1D;

/*----------------------------------------------------------------------------*/
/* structure Prim1D:  primitive variables in 1D (does not contain Bx)
 *  IMPORTANT!! The order of the elements in Prim1D CANNOT be changed.
 */

typedef struct Prim1D_s{
  Real d;			/* density */
  Real Vx;			/* Velocity in X,Y,Z */
  Real Vy;
  Real Vz;
#ifndef BAROTROPIC
  Real P;			/* pressure */
#endif /* BAROTROPIC */
#ifdef MHD
  Real By;			/* cell centered magnetic fields in Y,Z */
  Real Bz;
#endif /* MHD */
#if (NSCALARS > 0)
  Real r[NSCALARS];              /* density-normalized advected scalars */
#endif
}Prim1D;

/*--------------------------------------------------------------------------*/
/* Grain structure: Basic quantities for one particle.
 * Note: One particle here represents a collection of billions of real particles
 * Created: Emmanuel Jacquet, May 2008;  Modified: Xuening Bai, Dec. 2008
 */

#ifdef PARTICLES

/* Physical quantities of a grain particle */
typedef struct Grain_s{
  Real x1,x2,x3;		/* coordinate in X,Y,Z */
  Real v1,v2,v3;		/* velocity in X,Y,Z */
  int property;			/* index of grain properties */
#ifdef FARGO
  Real shift;			/* amount of shift in x2 direction */
#endif
}Grain;

/* List of physical grain properties */
typedef struct Grain_Property_s{
#ifdef FEEDBACK
  Real m;		/* mass of this type of particle (g) */
#endif
  Real rad;		/* radius of this type of particle (cm) */
  Real rho;		/* solid density of this type of particle (g/cm^3) */
  long num;		/* number of particles with this property */
}Grain_Property;

#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/* structure Ray: a structure for storing information about a ray
 * passing through the grid, which we use to do radaitive transfer.
 */
#ifdef ION_RADPOINT

typedef struct Ray_s {
  Real n1, n2, n3;         /* Unit vector specifying ray direction */
  Real x1_0, x2_0, x3_0;   /* Physical position of ray origin */
  Real x1_1, x2_1, x3_1;   /* Physical position of ray end */
#ifdef MPI_PARALLEL
  int nproc;               /* Number of processor domains the ray intersects */
  int *proc_list;          /* List of processor domains through which the ray
			    * passes
			    */
  Real *x1_dom;            /* Coordinates where ray enters and exits
			    * processor domains
			    */
  Real *x2_dom;
  Real *x3_dom;
  Real *proc_length_list;  /* Length of ray passing through proc domain */
  int *ncell;
  int **i_list;            /* List of cells through which this ray passes, 
			    * indexed by processor, then by cell number
			    */
  int **j_list;
  int **k_list;
  Real **cell_length_list; /* Length of ray passing through cell */
#else
  int ncell;
  int *i_list;             /* List of cells through which this ray passes, 
			    * indexed by cell number
			    */
  int *j_list;
  int *k_list;
  Real *cell_length_list;  /* Length of ray passing through cell */
#endif
  Real flux;               /* Ionizing flux */
  int leaf;                /* Is this ray a leaf? */
  Real length;             /* Length of this ray segment, start to end */
  Real cum_length;         /* Length of ray up to this point from free base */
  struct Ray_s *child[4];  /* Pointer to 4 children */
  struct Ray_s *parent;    /* Pointer to parent */
  int lev;                 /* Level of the ray */
  int mynum, levnum;       /* Number of the ray, relative to parent
			      and on its level as a whole */
} Ray;

/*----------------------------------------------------------------------------*/
/* structure Rad_Ran2_State: The internal state of the random number
 * generator used by the ionizing radiative transfer routines
 */
#define NTAB_RAN2 32
typedef struct Rad_Ran2_State_s {
  long seed, idum2, iy;
  long iv[NTAB_RAN2];
} Rad_Ran2_State;

/*----------------------------------------------------------------------------*/
/* structure Ray_Tree: a tree of rays
 */
typedef struct Ray_Tree_s {
  Real x1, x2, x3;         /* Physical position of tree center */
  Real x1blo, x2blo, x3blo; /* Bounding box around tree */
  Real x1bhi, x2bhi, x3bhi;
  float rotation[3][3];    /* Rotation matrix describing tree orientation */
  int max_level;           /* Maximum level of the tree */
  int rebuild_ctr;         /* Counter to count time to rebuild tree */
  Rad_Ran2_State ranstate; /* State of random number generator
			      associated with this tree */
  Ray *rays;               /* List of rays */
} Ray_Tree;

/*----------------------------------------------------------------------------*/
/* structure Radpoint:  An object that describes a point source of ionizing
 * radiation.
 */
typedef struct Radpoint_s {
  Real x1, x2, x3;           /* Physical position of radiator */
  Real s;                    /* Ionizing luminosity of source, in units
				of 1/time (i.e. photons / sec, not erg /sec */
  Ray_Tree tree;             /* The ray tree associated with this
				radiator */
} Radpoint;

#endif /* ION_RADPOINT */

#ifdef ION_RADPLANE
/*----------------------------------------------------------------------------*/
/* structure Radplane:  An object that describes source of
 * plane-parallel ionizing radiation entering the grid on one side.
 */
typedef struct Radplane_s {
  int dir;                   /* Direction of radiation propagation: -1
				= -x direction, +1 = +x direction, -2
				= -y direction, +2 = +y direction,
				etc. */
  Real flux;                 /* Ionizing flux, in units
				of 1/time/area (i.e. photons / sec /
				cm^2, not erg /sec / cm^2 */
} Radplane;

#endif /* ION_RADPLANE */

/*----------------------------------------------------------------------------*/
/* structure Grid: All data needed by a single processor to integrate equations
 *   Initialized by init_grid().  By using an array of Gas, rather than arrays
 *   of each variable, we guarantee data for each cell are contiguous in memory.
 */

typedef struct Grid_s{
  Gas ***U;			/* conserved variables */
#ifdef MHD
  Real ***B1i,***B2i,***B3i;    /* interface magnetic fields */
#endif /* MHD */
#ifdef SPECIAL_RELATIVITY
  Prim ***W;                    /* primitive variables, needed with SR */
#endif /* SPECIAL_RELATIVITY */
#ifdef SELF_GRAVITY
  Real ***Phi, ***Phi_old;      /* gravitational potential */
  Real ***x1MassFlux;           /* x1 mass flux for source term correction */
  Real ***x2MassFlux;           /* x2 mass flux for source term correction */
  Real ***x3MassFlux;           /* x3 mass flux for source term correction */
#endif /* GRAVITY */
  Real x1_0;	            /* x1-position of coordinate ix = 0 */
  Real x2_0;	            /* x2-position of coordinate jx = 0 */
  Real x3_0;	            /* x3-position of coordinate kx = 0 */
  Real dx1,dx2,dx3;         /* cell size */
  Real dt,time;		    /* time step, absolute time */
  int nstep;		    /* number of integration steps taken */
  int Nx1,Nx2,Nx3;          /* number of zones in each direction in this Grid */
  int is,ie;		    /* start/end cell index in x1 direction */
  int js,je;		    /* start/end cell index in x2 direction */
  int ks,ke;		    /* start/end cell index in x3 direction */
  int idisp;                /* coordinate ix = index i + idisp */
  int jdisp;                /* coordinate jx = index j + jdisp */
  int kdisp;                /* coordinate kx = index k + kdisp */
  char *outfilename;        /* basename for output files */

#ifdef PARTICLES
  Grain *particle; 		/* linked list of all particles */
  long arrsize;			/* size of the particle array */
  long nparticle;		/* number of particles */
  Grain_Property *grproperty;	/* array of particle properties of all types */
  int partypes;			/* number of particle types types (size, density, mass) */
#ifdef FEEDBACK
  Vector ***feedback;		/* array of feedback force to grid */
#endif /* FEEDBACK */
#endif /* PARTICLES */

#ifdef ION_RADPOINT
  int nradpoint;            /* number of point radiatiors */
  Radpoint *radpointlist;   /* list of radiators */
#endif
#ifdef ION_RADPLANE
  int nradplane;            /* number of planar rad fronts */
  Radplane *radplanelist;   /* list of radiation plane sources */
#endif
  int my_id;                /* process ID (or rank in MPI) updating this Grid */
  int rx1_id, lx1_id;       /* ID of grids to R/L in x1-dir (default = -1) */
  int rx2_id, lx2_id;       /* ID of grids to R/L in x2-dir (default = -1) */
  int rx3_id, lx3_id;       /* ID of grids to R/L in x3-dir (default = -1) */
}Grid;

/*----------------------------------------------------------------------------*/
/* structure GridIndices: indices and IDs of each Grid in Domain
 * structure Domain: 3D array of GridIndices, and all other information about
 *   this Domain that might be needed by an individual processor updating an
 *   individual Grid
 */

typedef struct Grid_Indices_s{
  int igs, jgs, kgs;             /* Minimum coordinate of cells in this Grid */
  int ige, jge, kge;             /* Maximum coordinate of cells in this Grid */
  int id;                        /* process ID (rank in MPI) */
}Grid_Indices;

typedef struct Domain_s{
  Grid_Indices ***GridArray;     /* 3D array of Grids tiling this Domain */
  int ids, jds, kds;        /* Minimum coordinate of cells over entire Domain */
  int ide, jde, kde;        /* Maximum coordinate of cells over entire Domain */
  int NGrid_x1;             /* Number of Grids in x1 direction for this Domain */
  int NGrid_x2;            /* Number of Grids in x2 direction for this Domain */
  int NGrid_x3;            /* Number of Grids in x3 direction for this Domain */
  int Nx1,Nx2,Nx3;  /* total number of zones in each direction over all Grids */
}Domain;

/*----------------------------------------------------------------------------*/
/* structure Output: everything for outputs */
  
struct Output_s;
typedef void (*VGFunout_t)(Grid *pGrid, Domain *pD, struct Output_s *pout);
typedef Real (*Gasfun_t)(const Grid *pG, const int i, const int j, const int k);

typedef struct Output_s{
  int n;          /* the N from the <outputN> block of this output */
  Real dt;        /* time interval between outputs  */
  Real t;         /* next time to output */
  int num;        /* dump number (0=first) */
  char *out;      /* variable (or user fun) to be output */
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


/* typedefs for functions that compute static gravitational potential and
 * cooling, set in problem generator, and used by integrators */
typedef Real (*GravPotFun_t)(const Real x1, const Real x2, const Real x3);
typedef Real (*CoolingFun_t)(const Real d, const Real p, const Real dt);

/* Directions for the set_bvals_fun() function */
enum Direction {left_x1, right_x1, left_x2, right_x2, left_x3, right_x3};

/* Definitions of various functions */
typedef void (*VBCFun_t)(Grid *pG);    /* void boundary cond fn */
typedef void (*VGFun_t) (Grid *pG);    /* void grid function */
typedef void (*VGDFun_t)(Grid *pG, Domain *pD);     /*void grid + domain func */
#ifdef PARTICLES
/* describe the gas velocity difference due to the pressure gradient */
typedef void (*GVDFun_t)(const Real x1, const Real x2, const Real x3, Real *u1, Real *u2, Real *u3);
#endif /* PARTICLES */

#endif /* ATHENA_H */
