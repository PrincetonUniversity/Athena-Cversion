---
title: Particle List
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Particle List Output

The binary format "`lis`" is designed to output a list of individual particles without any loss of information. One example of setting the `lis` output is the following


	<output1>
	out_fmt  = lis               # particle list data output
	dt       = 3.0               # time step between output
	id       = ds
	par_prop = limit             # user defined particle selection function


Here the user-defined [particle property selection function]({{site.baseurl}}/AthenaDocsParOutOverview) "`limit`" is used to track only a subsample of all the particles (hence the `id` name "`ds`" means down sampling). The output contains the total number of particles in the file, as well as complete information about all particle types and all individual particles (position, velocity, id, etc.) in the file.

A matlab script `/atheta/vis/particle/Read_par_standard.m` can be used to read the binary file and make simple plots, based on which the users can make their own plotting scripts.

The script join_lis.c can be used to combine the lis files from different processors in MPI jobs into a single lis file. The script sort_lis.c can be used to sort the particles in the lis file based on their identifiers. Both C scripts are located in `/athena/vis/particle`.
