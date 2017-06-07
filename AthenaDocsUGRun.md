---
title: Running Athena
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Running

After editing the input file, the Athena executable can be run with the `-i` option to
specify the name of
the input file.  For example, to run the Brio & Wu shocktube use

	% athena -i ../tst/1D-mhd/athinput.brio-wu

The code will first echo the values of all input parameters to stdout.
During the main integration loop
it will print the cycle number, current and next timestep, and
when it concludes it will print final diagnostic information.

A variety of command line options have been implemented in Athena.
A list is given by the `-h` switch:


	% athena -h
	Athena version 4.0 - 01-JUN-2010
	  Last configure: Thu Mar 25 22:15:28 EDT 2010
	
	Usage: athena [options] [block/par=value ...]
	
	Options:
	  -i <file>       Alternate input file [athinput]
	  -r <file>       Restart a simulation with this file
	  -d <directory>  Alternate run dir [current dir]
	  -h              This Help, and configuration settings
	  -n              Parse input, but don't run program
	  -c              Show Configuration details and quit
	  -t hh:mm:ss     With MPI, wall time limit for final output
	
	Configuration details:
	
	 Problem:                 linear_wave1d
	 Gas properties:          MHD
	 Equation of State:       ADIABATIC
	 Passive scalars:         0
	 Self-gravity:            OFF
	 Ohmic resistivity:       OFF
	 Viscosity:               OFF
	 Thermal conduction:      OFF
	 Particles:               OFF
	 Coordinate System:       Cartesian
	 Special Relativity:      OFF
	 Order of Accuracy:       2 (SECOND_ORDER_CHAR)
	 Flux:                    roe
	 Unsplit integrator:      ctu
	 Precision:               DOUBLE_PREC
	 Ghost cell Output:       OFF
	 Parallel Modes: MPI:     OFF
	 H-correction:            OFF
	 FFT:                     OFF
	 Shearing Box:            OFF
	 FARGO:                   OFF
	 All-wave integration:    OFF
	 Static mesh refinement:  OFF


The `-d` option can be used to create a new directory in which
Athena will run and write the output files.  The `-n` option is
useful for debugging any parsing errors, as it will dump the contents
of all parsed block/parameters.

A value for any of the valid parameter names in the input file can
also be input from the command line, this over-rides the values in the
input file.  This, in combination with the `-d` option, is useful
for parameter surveys.  The `-c` option is useful for checking the
configuration parameters with which the executable was compiled.
