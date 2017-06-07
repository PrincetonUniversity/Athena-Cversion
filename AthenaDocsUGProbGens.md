---
title: Problem Generators
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Problem Generators

What is a Problem Generator?
============================

The purposes of the problem generator are to
 * set initial conditions for all variables, in a function that must be called `problem` and has the following prototype

	void problem(DomainS *pDomain)

 * enroll any problem-specific boundary conditions, if required.  See [Boundary Conditions]({{site.baseurl}}/AthenaDocsUGBCs).
 * enroll any problem-specific user-defined history outputs, if required.  See [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress).
 * enroll any problem-specific physics controlled by function pointers, like forces due to a static gravitational potential, or optically-thin cooling.
It will probably be necessary to read the [Programmer Guide]({{site.baseurl}}/AthenaDocsPG) in order to understand the data structures and Mesh in Athena well
enough to write a new problem generator.  The existing files in the `/src/prob` directory can be used as starting templates.

In addition, the file containing the `problem()` function must also contain a number of other required functions

	* problem_write_restart() - writes problem-specific user data to restart files
	* problem_read_restart()  - reads problem-specific user data from restart files
	* get_usr_expr()          - sets pointer to expression for special output data
	* get_usr_out_fun()       - returns a user defined output function pointer
	* Userwork_in_loop        - problem specific work IN     main loop
	* Userwork_after_loop     - problem specific work AFTER  main loop

In particular, see [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress) for a description of how to use the function
`get_usr_expr()` to add new user-defined output variables using one of the existing file formats, and see
[User-defined Output Formats]({{site.baseurl}}/AthenaDocsUGUserFormats) for a description of how to use the function
`get_usr_out_fun()` to add new user-defined output formats.

It may also be necessary to include special user-defined functions in the same file that contains the problem generator if
new output variables, or new output formats, are used.

A large number of problem generators are included in the 
[./athena/src/prob](https://github.com/PrincetonUniversity/Athena-Cversion/tree/master/src/prob)
directory.

Parsing the Input File
======================

As metioned in the section on [Input Files]({{site.baseurl}}/AthenaDocsUGInputFile), data in the input file must be read using functions
defined in a parser written for Athena located in the file `src/par.c`.  It is quite likely that every problem generator
will require the input of at least a few parameters from the `<problem>` block in the input file.  The following functions can be
used in the problem generator to read data from the input file.

	char  *par_gets(char *block, char *name);  /* reads a string called "name" from input block "block" */
	int    par_geti(char *block, char *name);  /* reads an integer called "name" from input block "block" */
	double par_getd(char *block, char *name);  /* reads a double called "name" from input block "block" */

The following three functions do the same thing, except set the value to that given in the `def`, if the name
cannot be found in the input block.  This is useful for setting default values to parameters without having to
always include them in the input file.

	char  *par_gets_def(char *block, char *name, char   *def);
	int    par_geti_def(char *block, char *name, int    def);
	double par_getd_def(char *block, char *name, double def);

Most of the problem generators in the `src/prob` have examples of the usage of the above functions.  Below are some examples.

	/* Read problem parameters.  Note Omega_0 set to 10^{-3} by default */
	Omega_0 = par_getd_def("problem","Omega",1.0e-3);
	qshear  = par_getd_def("problem","qshear",1.5);
	amp = par_getd("problem","amp");
	filename = par_gets("problem","fname");

Finally, sometimes it is useful to set parameters that are not already defined a specific `<input block>`.
The following functions can be used for this purpose.

	void par_sets()       - sets/adds a string
	void par_seti()       - sets/adds an integer
	void par_setd()       - sets/adds a Real

