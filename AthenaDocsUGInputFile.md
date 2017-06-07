---
title: The Input File
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Input File

Before running the code, run time parameters must be set in an *Input File*.

Usually this file is given the name
`athinput.`*problem-name*, where *problem-name* is a
string identifier, often the same as the name of the
problem-generator, i.e. the function in `./athena/src/prob` which
is used to initialize the data.  For example, the 1D MHD Brio & Wu shocktube test uses the input file *athinput.brio-wu*.
A large number of input files for different tests are provided in `./athena/tst`.

Within the input file, parameters are grouped into named *blocks*, with the name of each
block appearing on a single line within angle brackets, for example

	<job>
	problem_id      = Strat      # problem ID: basename of output filenames
	maxout          = 5          # Output blocks number from 1 -> maxout
	num_domains     = 1          # number of Domains in Mesh

Block names must always appear in angle brackets on a separate line
(blank lines above and below the block names are not required, but can be used for clarity).

Below each block name is a list of parameters, with syntax

	parameter-name = value [# comments]

White space after the parameter name, after the `=', and before the
# character is ignored.  Everything after (and including) the #
character is also ignored.  Only one parameter value can appear per
line.  Comment lines (i.e. a line beginning with #) are allowed for
documentation purposes.  A maximum number of 256 characters is permitted
per line in the input file.  Both block names and parameter names are
case sensitive.

The input file is read by a very flexible parser written for Athena
`./athena/src/par.c`.  The entire input file is read at the
very beginning of the main program, and the parameter names and their
values stored in memory.  Thereafter, these values can be accessed as
necessary by any function at any time during execution.  The parser
allows the parameter names to appear in any order within a named
block, extra (or misspelled) parameter names will be parsed and never
used.  If a value is requested from the parser but its name does
not exist, the parser will print an error message and terminate the program. 
The parameters may be
integers, floating point numbers, or strings, but the input values must be of the right type (for example, strings
cannot be used to enter integers).  The parser will do
automatic type conversion, for example converting floating point
numbers to double precision if necessary.  Parameter values can also be set at run time through the
command line, see [Tutorial/Command Line]({{site.baseurl}}/AthenaDocsTutCL).

Functions defined in `par.c` to read the input file are described in the
section on [Problem Generators]({{site.baseurl}}/AthenaDocsUGProbGens).
