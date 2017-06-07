---
title: Command Line Interface
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[Tutorial]({{site.baseurl}}/AthenaDocsTut)/Command Line

Modifying Input Parameters
==========================

Any existing parameter in the input file can also be modified at run time using the command line interface.

For example, to change the grid resolution in the 2D Orszag-Tang vortex test, simply run Athena using

	% athena -i ../tst/2D-mhd/athinput.orszag-tang domain1/Nx1=128 domain1/Nx2=128

This will have exactly the same result as changing the parameters by editing the input file, as described in the first
section of the [Editing the Input File]({{site.baseurl}}/AthenaDocsTutIF) tutorial.

Similarly, to change the format of outputs from binary to vtk in this same test, just use

	% athena -i ../tst/2D-mhd/athinput.orszag-tang output2/out_fmt=vtk

This will have the same effect as editing the `<output2>` block in the input file.

Finally, only parameters that currently exist in the input file can be modified using the command line.  Thus, it is not
possible to add an entirely new output block (for example, the `<output5>` block described in the 
[Editing the Input File]({{site.baseurl}}/AthenaDocsTutIF) tutorial), since this block does not exist in the input file to begin with.
Thus, *adding new parameters requires editing the input file.*

Other Command Line Options
==========================

The command line interface can be used to control more than just changing parameter values in the input file.  To see the
available options, use

	% athena -h

The code will output a variety of diagnostic information, including the configuration details, as well as a list of
available command line options

	Usage: athena [options] [block/par=value ...]
	
	Options:
	  -i <file>       Alternate input file [athinput]
	  -r <file>       Restart a simulation with this file
	  -d <directory>  Alternate run dir [current dir]
	  -h              This Help, and configuration settings
	  -n              Parse input, but don't run program
	  -c              Show Configuration details and quit
	  -t hh:mm:ss     With MPI, wall time limit for final output

The `-d` option is particularly useful for outputting data to a directory other than the current.
