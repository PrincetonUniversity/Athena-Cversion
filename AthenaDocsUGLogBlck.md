---
title: The <Log> Block
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Log Block

Parameters in this optional block control output to stdout and stderr.  On many
machines, output to stderr and stdout does not scale well (to $10^4-10^5$
processors).  However, having output from all MPI processes can be useful
for debugging.  Thus, some kind of runtime control that allows output to
stderr and stdout be turned off for production runs, or redirected to
files for debugging, is very helpful.  

This control is providing by adding integer flags called *out_level* and *err_level* which
control how verbose output will be from the root (rank=0) and child processes
with MPI.  The following table describes the behavior of output with different values of these
flags:

| out_level or err_level | behavior |
|------------------------|----------|
| -1 | only fatal errors produce output |
| 0  | only root (rank=0) process produces output |
| 1  | all processes produce output |

In the example below, each of the parameters
in this block are set to their default values.

	<log>
	file_open = 0
	iflush = 0
	lazy = 1
	out_level = 0
	err_level = 0
	child_out_level = -1
	child_err_level = -1


**file_open:** Set to one to direct output from stderr to a file named *problem_id*.err, and output from stdout to a file
named *problem_id*.out, where *problem_id* is the problem basename set in the `<job>` block.

**iflush:**  Number of cycles between flush of buffers.  Useful to guarantee error output every cycle.

**lazy:**  Set to zero to force opening of .err and .out files if `file_open=1`,
even if no output has been generated.  If lazy is not
equal to zero, these files are only opened if there is output to be written.  This is useful for
parallel jobs where the children, if made sufficiently quiet, would
otherwise generate a large number of empty files.

**out_level, err_level:** Out and err level of root process.

**child_out_level, child_err_level:**  Out and err level of child processes.

The `<log>` block is optional, and if not included in the input file, then the default values for all parameters
given above will be used.  Most input files in the ./athena4.0/tst directories do not include it.
