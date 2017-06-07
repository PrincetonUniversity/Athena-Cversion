---
title: Domain Structure
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[ProgrammerGuide]({{site.baseurl}}/AthenaDocsPG)/Domain Structure

Another important data structure defined in `/athena/src/athena.h` is `DomainS`

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

This structure stores the edges of the Domain in the (x1,x2,x3) coordinate system,
the integer displacement of the Domain from the origin of the root Domain, and information
needed for SMR and parallelization with MPI.  It also includes a pointer to a GridS structure,
which stores the Grid on this Domain.  With MPI and SMR, it is possible that there is no Grid
on this Domain being updated by this processor, in which case `*GridS` will be the
`NULL` pointer.  In any case, we also store information about *all* Grids on this Domain in a 3D
array with dimensions equal to the domain decomposition used with MPI.
