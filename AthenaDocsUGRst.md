---
title: Restart Files
layout: page
---
[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Restart Files

Restart files (sometimes called *checkpoints*) are useful when a
calculation must be continued from a previous point.  The files contain
enough information, and with the necessary accuracy, that a restart calculation
generates *identical* data (to all significant digits) to a calculation run
continuously.  Athena defines its
own format for restart files.  The file `restart.c` contains all
the functions needed to read and write these files.

The following example shows an `<output>` block in an
input file that generates restart files:

        <output2>
        out_fmt = rst            # Restart dump
        dt      = 1.0            # time increment between outputs

The time increment `<output>/dt` is measured in problem time,
and should be set to give the desired output frequency of files
(usually writing one restart dump every 6 hours of wall clock time is useful).
For jobs run in parallel with MPI, there will be one restart file per process,
and the restarted job must use the same number of processors as there are
restart files.

If the problem contains special, user-defined data, these must be added to the
restart dumps.  Athena provides a mechanism for automatically adding
such data.  In the problem generator, two functions are provided:

        void problem_write_restart(MeshS *pM, FILE *fp)
        {
          return;
        }
        
        void problem_read_restart(MeshS *pM, FILE *fp)
        {
          return;
        }

Generally these functions are empty, but if necessary they can be used
to read and write extra parameters using `fwrite` and `fread` statements, and the file
pointer `*fp` passed as the second argument to these functions.  They can also set problem-specific boundary
conditions on restart, and enroll user-defined outputs and variables.  The problem generator `src/prob/rt.c`
contains an example of usage.  See the [Programmer Guide]({{site.baseurl}}/AthenaDocsPG) for more
information about the structures in the argument list to these functions.

To read a restart file, the `-r` command is used on the command line:

        % athena -r myfile.rst

Note that an input file, specified by `-i myinput`, is not needed for
restarts.  This is because the restart file contains the original input file,
in ASCII format, at the beginning, from which all the necessary parameters
are read by `par.c` on restart.  
This also makes restart files self-documenting: the values of input parameters used in the calculation that generated
the restart file can be read with an editor.  If an input file is specified
along with a restart,

        % athena -r myfile.rst -i myinput

then the values in `myinput` overwrite the values stored in the restart file
itself.  Alternatively, values in the input file can be overwritten using the
command line, for example:

        % athena -r myfile.rst time/tlim=20.0 

Usually the `time/tlim` parameter needs to be modified on restart.
For parallel jobs run with MPI, only the name of the restart file for the
root (rank 0) process needs to be specified, all other processes will
create their own appropriate restart filename based on this name.  See [Running Athena on Parallel Processors]({{site.baseurl}}/AthenaDocsUGParallel).
