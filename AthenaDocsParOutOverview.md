---
title: Overview of Particle Data Output
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Output Overview

There are two choices for the particle output: 1). Bin the particle-related quantities to the grid and output arrays of grid-binned quantities; 2). Output a list of individual particles. Both choices are supported in Athena.

Enrolling particle output is similar to [that in the MHD code]({{site.baseurl}}/AthenaDocsUGOutputOverview): one creates an output block in the input file specifying the output format, variable, time interval and so on. With particles, two additional items, namely, particle binning switch and particle property selection function, are to be specified in the output block (in many cases they are set by default). An example of their usage is shown as follows. 


	<output1>
	  ...
	  pargrid	= 1		/* bin particles to grid (=1) or not (=0) */
	  par_prop	= limit		/* particle property selection function */


**pargrid**: whether quantities associated with the particles are to be binned to the grid (yes=1, no=0). Particle binning is automatically assigned in most situations based on the output format and output variable in the block, but it needs to be set explicitly if the output is user-defined.

**par_prop**: the particle property selection function. By default, all the particles are selected (`par_prop=all`), but one can properly define the function by selecting a subsample of the particles for the output. This is particularly useful if one wants to track the trajectories of a small number of particles in a big simulation.

One example of particle property selection function is shown below


	static int property_limit(const GrainS *gr, const GrainAux *grsub)
	{
	  if ((gr->pos == 1) && (gr->my_id<nlis))
	    return 1;
	  else
	    return 0;
	}


In this function, the parameter `nlis` is previously assigned as a local variable in the problem generator. The function select particles that are within the grid (`pos==1`) and whose particle identifier number `my_id` is less than `nlis`.

To enroll this selection function, one chooses a name (here "limit") for this selection function and set `par_prop` in the output block with this name. Then, in the problem generator, a function pointer is to be returned to this function if the string in the argument (read from the output block) matches "limit".


	PropFun_t get_usr_par_prop(const char *name)
	{
	  if (strcmp(name,"limit")==0)    return property_limit;
	
	  return NULL;
	}
