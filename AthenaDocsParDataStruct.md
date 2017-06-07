---
title: Particle Data Structures
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Data Structure

The information for a single particle is stored in the structure named `GrainS`, which records the particle position, velocity, property, and identifier. The declaration of the GrainS structure is shown as below:


	typedef struct Grain_s{
	  Real x1,x2,x3;        /* coordinate in X,Y,Z */
	  Real v1,v2,v3;        /* velocity in X,Y,Z */
	  int property;         /* index of grain properties */
	  short pos;            /* position: 0: ghost; 1: grid; >=10: cross out/in; */
	  long my_id;           /* particle id */
	#ifdef MPI_PARALLEL
	  int init_id;          /* particle's initial host processor id */
	#endif
	}GrainS;


Particles are associated with the grid. An array of the above particle structure is contained in the structure `GridS`, which stores all the particles whose positions are within the grid. Particles are passed to other processors via MPI calls when they cross the grid boundaries.

Additional information about the particles are recorded in the structure `Grain_Property` defined as below: 


	typedef struct Grain_Property_s{
	#ifdef FEEDBACK
	  Real m;               /* mass of this type of particle */
	#endif
	  Real rad;             /* radius of this type of particle (cm) */
	  Real rho;             /* solid density of this type of particle (g/cm^3) */
	  long num;             /* number of particles with this property */
	  short integrator;     /* integrator type: exp (1), semi (2) or full (3) */
	}Grain_Property;


The properties include the quantities that are needed to compute the particle stopping time (`rad` and `rho`), particle mass for calculating the back reaction, as well as the index of particle integrator used to integrate this type of particles. See [Particle Input Block]({{site.baseurl}}/AthenaDocsParBlock) and [Particle Problem Generator]({{site.baseurl}}/AthenaDocsParProb) for more details.

An arbitrary number particle types are allowed, with each type having different particle properties. The number of particle types is recorded in the global variable `npartypes`, and global array 

	Grain_Property *grproperty;

stores the information of each type of particles. The element `property` in the structure `GrainS` records the property type index that the particle belongs.

In numerical simulations, it is common practice that constant particle stopping time for each particle type is used. For this purpose, we further define a global array 

	Real *tstop0;

which records the stopping time of each particle types.
