---
title: New Problems
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/NewProblems

The real utility of the Athena code is as a solver for new problems
(i.e. problems that 
are not initialized by the set of problem generators included in the source
code distribution).  For new problems, the following steps are required.

1. Write a new function that initializes the problem, contained in a file in the `./athena/src/prob` directory, with the following prototype

        void problem(DomainS *pDomain)

    The main purpose of this function is to set the initial conditions for all dependent variables over the entire Domain.  Values in the ghost
    zones do not need to be set, since Athena calls the boundary condition functions immediately after the problem generator.  However, values for
    the magnetic field on the faces of the Domain (`B1i` at `is` and `ie+1`; `B2i` at `js` and `je+1`; and
    `B3i` at `ks` and `ke+1`) must be set, since these values are *evolved* by Athena.

2. Write new functions called `problem_write_restart` and `problem_read_restart` (which may be no-ops if needed) and include them in the
   problem generator file.  These functions write and read any problem-dependent parameters to the restart files, or read them from the input file on
   restarts, and also set any problem-specific outputs or boundary conditions on restarts.  As an example, the following code reads some parameters
   and sets new history variables, but does not write any new parameters to the restart file

        void problem_write_restart(MeshS *pM, FILE *fp)
        {
          return;
        }
        
        void problem_read_restart(MeshS *pM, FILE *fp)
        {
          Omega_0 = par_getd_def("problem","omega",1.0e-3);
        
          dump_history_enroll(hst_rho_Vx_dVy, "<rho Vx dVy>");
          dump_history_enroll(hst_rho_dVy2, "<rho dVy^2>");
        
          return;
        }

    **These functions must be present in the problem generator file, even if they are just no-ops.**

3. Write new functions called `Userwork_in_loop` and `Userwork_after_loop` (which may be no-ops if not needed) and include them in the file containing `problem`.  As the names suggest, these functions can be used to perform special problem-dependent work in or after the main loop (see `./athena/src/prob/linear_wave.c` for an example).    **These functions must be present in the problem generator file, even if they are just no-ops.**

4. If special purpose boundary conditions are needed, write special functions that implement them, and enroll them using the function `bvals_mhd_fun` see [Boundary Conditions]({{site.baseurl}}/AthenaDocsUGBCs).

5. If special purpose data output is needed, write special functions that implement them (see [User-defined Output Variables]({{site.baseurl}}/AthenaDocsUGUserExpress) and [User-defined Output Formats]({{site.baseurl}}/AthenaDocsUGUserFormats)).  Adding new variables requires the use of the `get_usr_expr()` function, while adding new formats requires the `get_usr_out_fun()` function.  Both of these functions **must be present in the problem generator file, even if they are just no-ops.**

6. Once the above is complete, configure and compile the code using the appropriate physics options, and including the new problem generator using `--with-problem=`*new-name*, where *new-name* is the file created in step (1).


It is likely the [Programmer's Guide]({{site.baseurl}}/AthenaDocsPG) will be needed to write a
new problem generator to understand the data structures, variable names, and array indexing used
in Athena.  As a start, the problem generators in `./athena/src/prob` can be used
as templates.

