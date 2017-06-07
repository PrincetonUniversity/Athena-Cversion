---
title: Grid Structure
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Grid Structure

We organize the information about each Grid in the Domain, as well as the dependent
variables on that Grid, into a C structure called `GridS`.  The conserved
variables *at cell center* are organized into another C structure `ConsS`.
Thus, a GridS contains a 3D array of the ConsS structure, a 3D array of each component of the face centered
magnetic field, plus information about coordinates.
A portion of the GridS structure is reproduced below.

	typedef struct Grid_s{
	  ConsS ***U;                /* conserved variables */
	#ifdef MHD
	  Real ***B1i,***B2i,***B3i;    /* interface magnetic fields */
	#endif /* MHD */
	#ifdef SELF_GRAVITY
	  Real ***Phi, ***Phi_old;      /* gravitational potential */
	  Real ***x1MassFlux;           /* x1 mass flux for source term correction */
	  Real ***x2MassFlux;           /* x2 mass flux for source term correction */
	  Real ***x3MassFlux;           /* x3 mass flux for source term correction */
	#endif /* GRAVITY */
	  Real MinX[3];         /* min(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
	  Real MaxX[3];         /* max(x) in each dir on this Grid [0,1,2]=[x1,x2,x3] */
	  Real dx1,dx2,dx3;        /* cell size on this Grid */
	  Real time, dt;           /* current time and timestep  */
	  int is,ie;               /* start/end cell index in x1 direction */
	  int js,je;               /* start/end cell index in x2 direction */
	  int ks,ke;               /* start/end cell index in x3 direction */
	  int Nx[3];       /* # of zones in each dir on Grid [0,1,2]=[x1,x2,x3] */
	  int Disp[3];     /* i,j,k displacements of Grid from origin [0,1,2]=[i,j,k] */
	
	  int rx1_id, lx1_id;   /* ID of Grid to R/L in x1-dir (default=-1; no Grid) */
	  int rx2_id, lx2_id;   /* ID of Grid to R/L in x2-dir (default=-1; no Grid) */
	  int rx3_id, lx3_id;   /* ID of Grid to R/L in x3-dir (default=-1; no Grid) */
	
	#ifdef CYLINDRICAL
	  Real *r,*ri;                  /* cylindrical scaling factors */
	#endif /* CYLINDRICAL */
	
	}GridS;

Some variables relating to SMR and particles have been deleted in the code reproduced above for clarity.

The advantage of making the cell centered variables a 3D
array of type `ConsS`, rather than making the `ConsS` structure itself contain
3D arrays of each variable, is that we guarantee that different variables at the same
grid cell will be stored contiguously in memory, which can improve cache
performance.
