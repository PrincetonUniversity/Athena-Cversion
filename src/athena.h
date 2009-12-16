#ifndef ATHENA_H
#define ATHENA_H 
/*==============================================================================
 * FILE: athena.h
 *
 * PURPOSE: Contains definitions of the following data types and structures:
 *   Real   - either float or double, depending on configure option
 *   Gas    - cell-centered conserved variables
 *   Prim   - cell-centered primitive variables
 *   Cons1D - conserved variables in 1D: same as Gas minus Bx
 *   Prim1D - primitive variables in 1D: same as Prim minus Bx
 *   Grain  - basic properties of particles
 *   GridS   - everything in a single Grid: arrays of Gas, B field, etc.
 *   DomainS - everything in a single Domain (potentially many Grids)
 *   MeshS   - everything across whole Mesh (potentially many Domains)
 *   OutputS - everything associated with an individual output
 *============================================================================*/
#include "defs.h"

#ifdef MPI_PARALLEL
#include "mpi.h"
#endif

/* variable type Real:   depends on macro set by configure
 */ 
#if defined(SINGLE_PREC)
typedef float  Real;
#elif defined(DOUBLE_PREC)
typedef double Real;
#else
# error "Not a valid precision flag"
#endif

/* general 3-vectors of Reals and integers 
 */
typedef struct Real3Vect_s{
  Real x, y, z;
}Real3Vect;
typedef struct Int3Vect_s{
  int i, j, k;
}Int3Vect;

/* sides of a cube, used to find overlaps between Grids at different levels 
 */
typedef struct Side_s{
  int ijkl[3];    /* indices of left-sides  in each dir [0,1,2]=[i,j,k] */ 
  int ijkr[3];    /* indices of right-sides in each dir [0,1,2]=[i,j,k] */ 
}SideS;

/* number of zones in, and identifying information about, a Grid
 */
typedef struct GridsData_s{
  int Nx[3];                /* number of zones in each dir [0,1,2]=[x1,x2,x3] */
  int Disp[3];     /* i,j,k displacements from origin of root [0,1,2]=[i,j,k] */
  int ID_Comm_world;      /* ID of process for this Grid in MPI_COMM_WORLD */
  int ID_Comm_Domain;     /* ID of process for this Grid in Comm_Domain    */
#ifdef STATIC_MESH_REFINEMENT
  int ID_Comm_Children;     /* ID updating this Grid in Comm_Domain    */
  int ID_Comm_Parent;     /* ID updating this Grid in Comm_Domain    */
#endif
}GridsDataS;

/*----------------------------------------------------------------------------*/
/* Gas structure: conserved variables 
 *  IMPORTANT!! The order of the elements in Gas CANNOT be changed.
 */

typedef struct Gas_s{
  Real d;			/* density */
  Real M1;			/* momentum density in 1,2,3 directions */
  Real M2;
  Real M3;
#ifndef BAROTROPIC
  Real E;			/* total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real B1c;			/* cell centered magnetic fields in 1,2,3 */
  Real B2c;
  Real B3c;
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];             /* passively advected scalars */
#endif
}Gas;

/*----------------------------------------------------------------------------*/
/* Prim structure: primitive variables, used with special relativity 
 *  IMPORTANT!! The order of the elements in Prim CANNOT be changed.
 */

typedef struct Prim_s{
  Real d;			/* density */
  Real V1;			/* velocity in 1,2,3 */
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
  Real r[NSCALARS];             /* density-normalized advected scalars */
#endif
}Prim;

/*----------------------------------------------------------------------------*/
/* Cons1D structure:  conserved variables in 1D (does not contain Bx)
 *  IMPORTANT!! The order of the elements in Cons1D CANNOT be changed.
 */

typedef struct Cons1D_s{
  Real d;			/* density */
  Real Mx;			/* momentum density in X,Y,Z; where X is     */
  Real My;                      /* direction longitudinal to 1D slice; which */
  Real Mz;                      /* can be in any dimension: 1,2,or 3         */
#ifndef BAROTROPIC
  Real E;			/* total energy density */
#endif /* BAROTROPIC */
#ifdef MHD
  Real By;			/* cell centered magnetic fields in Y,Z */
  Real Bz;
#endif /* MHD */
#if (NSCALARS > 0)
  Real s[NSCALARS];             /* passively advected scalars */
#endif
#ifdef CYLINDRICAL
  Real Pflux;	 		/* pressure component of flux */
#endif
}Cons1D;

/*----------------------------------------------------------------------------*/
/* Prim1D structure:  primitive variables in 1D (does not contain Bx)
 *  IMPORTANT!! The order of the elements in Prim1D CANNOT be changed.
 */

typedef struct Prim1D_s{
  Real d;			/* density */
  Real Vx;			/* velocity in X,Y,Z */
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
  Real r[NSCALARS];             /* density-normalized advected scalars */
#endif
}Prim1D;

/*----------------------------------------------------------------------------*/
/* Grain structure: Basic quantities for one particle.
 * Note: One particle here represents a collection of billions of real particles
 */

#ifdef PARTICLES

/* Physical quantities of a dust particle */
typedef struct Grain_s{
  Real x1,x2,x3;	/* coordinate in X,Y,Z */
  Real v1,v2,v3;	/* velocity in X,Y,Z */
  int property;		/* index of grain properties */
  short pos;		/* position: 0: ghost; 1: grid; >=10: cross out/in; */
  long my_id;		/* particle id */
#ifdef MPI_PARALLEL
  int init_id;          /* particle's initial host processor id */
#endif
#ifdef FARGO
  Real shift;           /* amount of shift in x2 direction */
#endif
}Grain;

/* Auxilary quantities for a dust particle */
typedef struct GrainAux_s{
  Real dpar;            /* local particle density */
#ifdef FARGO
  Real shift;           /* amount of shift in x2 direction */
#endif
}GrainAux;

/* List of physical grain properties */
typedef struct Grain_Property_s{
#ifdef FEEDBACK
  Real m;		/* mass of this type of particle */
#endif
  Real rad;		/* radius of this type of particle (cm) */
  Real rho;		/* solid density of this type of particle (g/cm^3) */
  long num;		/* number of particles with this property */
  short integrator;	/* integrator type: exp (1), semi (2) or full (3) */
}Grain_Property;

/* Grid elements for gas-particle coupling */
typedef struct GPCouple_s{
  Real grid_d;		/* gas density (at 1/2 step) */
  Real grid_v1;		/* gas velocities (at 1/2 step) */
  Real grid_v2;
  Real grid_v3;
#ifndef BAROTROPIC
  Real grid_cs;		/* gas sound speed */
#endif
#ifdef FEEDBACK
  Real fb1;             /* momentum feedback to the grid */
  Real fb2;
  Real fb3;
  Real FBstiff;         /* stiffness of the feedback term */
  Real Eloss;           /* energy dissipation */
#endif
}GPCouple;

#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/* Grid overlap structures, used for SMR
 */

#ifdef STATIC_MESH_REFINEMENT
typedef struct GridOvrlp_s{
  int ijks[3];           /* start ijk on this Grid of overlap [0,1,2]=[i,j,k] */
  int ijke[3];           /* end   ijk on this Grid of overlap [0,1,2]=[i,j,k] */
  int ID, DomN;          /* processor ID, and Domain #, of OVERLAP Grid */
  int nWordsRC, nWordsP; /* # of words communicated for Rest/Corr and Prol */
  Gas **myFlx[6];        /* fluxes of conserved variables at 6 boundaries */
#ifdef MHD
  Real **myEMF1[6];      /* fluxes of magnetic field (EMF1) at 6 boundaries */
  Real **myEMF2[6];      /* fluxes of magnetic field (EMF2) at 6 boundaries */
  Real **myEMF3[6];      /* fluxes of magnetic field (EMF3) at 6 boundaries */
#endif
}GridOvrlp;
#endif /* STATIC_MESH_REFINEMENT */

/*----------------------------------------------------------------------------*/
/* GridS: 3D arrays of dependent variables, plus grid data, plus particle data,
 *   plus data about child and parent Grids, plus MPI rank information for a
 *   Grid, where a Grid is defined to be the region of a Domain at some
 *   refinement level being updated by a single processor.  Uses an array of
 *   Gas, rather than arrays of each variable, to increase locality of data for
 *   a given cell in memory.  */

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
  Real x1min,x1max;        /* min/max of x1 on this Grid */
  Real x2min,x2max;        /* min/max of x2 on this Grid */
  Real x3min,x3max;        /* min/max of x3 on this Grid */
  Real dx1,dx2,dx3;        /* cell size on this Grid */
  Real time, dt;           /* current time and timestep  */
  int is,ie;		   /* start/end cell index in x1 direction */
  int js,je;		   /* start/end cell index in x2 direction */
  int ks,ke;		   /* start/end cell index in x3 direction */
  int Nx[3];       /* # of zones in each dir on Grid [0,1,2]=[x1,x2,x3] */
  int Disp[3];     /* i,j,k displacements of Grid from origin [0,1,2]=[i,j,k] */

  int rx1_id, lx1_id;   /* ID of Grid to R/L in x1-dir (default=-1; no Grid) */
  int rx2_id, lx2_id;   /* ID of Grid to R/L in x2-dir (default=-1; no Grid) */
  int rx3_id, lx3_id;   /* ID of Grid to R/L in x3-dir (default=-1; no Grid) */

#ifdef PARTICLES
  int partypes;              /* number of particle types */
  Grain_Property *grproperty;/* array of particle properties of all types */
  long nparticle;            /* number of particles */
  long arrsize;              /* size of the particle array */
  Grain *particle;           /* array of all particles */
  GrainAux *parsub;          /* supplemental particle information */
  GPCouple ***Coup;          /* array of gas-particle coupling */
#endif /* PARTICLES */

#ifdef STATIC_MESH_REFINEMENT
  int NCGrid;         /* # of child  Grids that overlap this Grid */
  int NPGrid;         /* # of parent Grids that this Grid overlaps */
  int NmyCGrid;       /* # of child  Grids on same processor as this Grid */
  int NmyPGrid;       /* # of parent Grids on same processor (either 0 or 1) */

  GridOvrlp *CGrid;     /* 1D array of data for NCGrid child  overlap regions */
  GridOvrlp *PGrid;     /* 1D array of data for NPGrid parent overlap regions */
/* NB: The order of the Grids in these two arrays is such that the first
 * NmyCGrid[NmyPGrid] elements contain overlap regions being updated by the
 * same processor as this Grid */
#endif /* STATIC_MESH_REFINEMENT */

}GridS;

typedef void (*VGFun_t)(GridS *pG);    /* generic void function of Grid */

/*----------------------------------------------------------------------------*/
/* DomainS: information about one region of Mesh at some particular level.
 *
 * Contains pointer to a single Grid, even though the Domain may contain many
 * Grids, because for any general parallelization mode, no more than one Grid
 * can exist per Domain per processor.
 *
 * The i,j,k displacements are measured in units of grid cells on this Domain
 */

typedef struct Domain_s{
  Real RootMinX[3];   /* min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3];   /* max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real MinX[3];       /* min(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real MaxX[3];       /* max(x) in each dir on this Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];                  /* cell size in this Domain [0,1,2]=[x1,x2,x3] */
  int Nx[3];      /* # of zones in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int NGrid[3];   /* # of Grids in each dir in this Domain [0,1,2]=[x1,x2,x3] */
  int Disp[3];   /* i,j,k displacements of Domain from origin [0,1,2]=[i,j,k] */
  int Level,DomNumber;   /* level and ID number of this Domain */
  int InputBlock;        /* # of <domain> block in input file for this Domain */
  GridS *Grid;       /* pointer to Grid in this Dom updated on this processor */

  GridsDataS ***GData; /* size,location, & processor IDs of Grids in this Dom */

  VGFun_t ix1_BCFun, ox1_BCFun;  /* ix1/ox1 BC function pointers for this Dom */
  VGFun_t ix2_BCFun, ox2_BCFun;  /* ix1/ox1 BC function pointers for this Dom */
  VGFun_t ix3_BCFun, ox3_BCFun;  /* ix1/ox1 BC function pointers for this Dom */

#ifdef MPI_PARALLEL
  MPI_Comm Comm_Domain;        /* MPI communicator between Grids on this Dom */
  MPI_Group Group_Domain;      /* MPI group for Domain communicator */
#ifdef STATIC_MESH_REFINEMENT
  MPI_Comm Comm_Parent;        /* MPI communicator to Grids in parent Domain  */
  MPI_Comm Comm_Children;      /* MPI communicator to Grids in  child Domains */
  MPI_Group Group_Children;    /* MPI group for Children communicator */
#endif /* STATIC_MESH_REFINEMENT */
#endif /* MPI_PARALLEL */
}DomainS;

typedef void (*VDFun_t)(DomainS *pD);  /* generic void function of Domain */

/*----------------------------------------------------------------------------*/
/* MeshS: information about entire mesh hierarchy, including array of Domains
 */

typedef struct Mesh_s{
  Real RootMinX[3];   /* min(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real RootMaxX[3];   /* max(x) in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  Real dx[3];     /* cell size on root Domain [0,1,2]=[x1,x2,x3] */
  Real time, dt;  /* current time and timestep for entire Mesh */
  int Nx[3];      /* # of zones in each dir on root Domain [0,1,2]=[x1,x2,x3] */
  int nstep;                 /* number of integration steps taken */
  int BCFlag_ix1, BCFlag_ox1;  /* BC flag on root domain for inner/outer x1 */
  int BCFlag_ix2, BCFlag_ox2;  /* BC flag on root domain for inner/outer x2 */
  int BCFlag_ix3, BCFlag_ox3;  /* BC flag on root domain for inner/outer x3 */
  int NLevels;               /* overall number of refinement levels in mesh */
  int *DomainsPerLevel;      /* number of Domains per level (DPL) */
  DomainS **Domain;          /* array of Domains, indexed over levels and DPL */
  char *outfilename;           /* basename for output files containing -id#  */
}MeshS;

/*----------------------------------------------------------------------------*/
/* OutputS: everything for outputs */
  
struct Output_s;
typedef void (*VOutFun_t)(MeshS *pM, struct Output_s *pout);
typedef void (*VResFun_t)(MeshS *pM, struct Output_s *pout);
typedef Real (*GasFun_t)(const GridS *pG, const int i,const int j,const int k);
#ifdef PARTICLES
typedef int (*PropFun_t)(const Grain *gr, const GrainAux *grsub);
typedef Real (*Parfun_t)(const Grid *pG, const Grain *gr);
typedef Real (*Parfun_t)(const Grid *pG, const Grain *gr);
#endif

typedef struct Output_s{
  int n;          /* the N from the <outputN> block of this output */
  Real dt;        /* time interval between outputs  */
  Real t;         /* next time to output */
  int num;        /* dump number (0=first) */
  char *out;      /* variable (or user fun) to be output */
  char *id;       /* filename is of the form <basename>[.idump][.id].<ext> */
#ifdef PARTICLES
  int out_pargrid;    /* output grid binned particles (=1) or not (=0) */
  PropFun_t par_prop; /* particle property selection function */
#endif

/* level and domain number of output (default = [0,0] = root level) */

  int nlevel, ndomain;

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
  VOutFun_t out_fun; /* output function pointer */
  VResFun_t res_fun; /* restart function pointer */
  GasFun_t expr;     /* pointer to expression that computes quant for output */

}OutputS;


/*----------------------------------------------------------------------------*/
/* typedefs for functions:
 */
/* for static gravitational potential and cooling, set in problem generator,
 * and used by integrators */

typedef Real (*GravPotFun_t)(const Real x1, const Real x2, const Real x3);
#ifdef CYLINDRICAL
typedef Real (*StaticGravAcc_t)(const Real x1, const Real x2, const Real x3);
#endif
typedef Real (*CoolingFun_t)(const Real d, const Real p, const Real dt);

#ifdef PARTICLES
/* function types for interpolation schemes and stopping time */
typedef void (*WeightFun_t)(GridS *pG, Real x1, Real x2, Real x3,
  Real3Vector cell1, Real weight[3][3][3], int *is, int *js, int *ks);
typedef Real (*TSFun_t)(GridS *pG, int type, Real rho, Real cs, Real vd);
#endif /* PARTICLES */

/*----------------------------------------------------------------------------*/
/* Directions for the set_bvals_fun() function */
enum BCDirection {left_x1, right_x1, left_x2, right_x2, left_x3, right_x3};

#endif /* ATHENA_H */
