---
title: Boundary Conditions
layout: page
---

[Documentation]({{site.baseurl}}/AthenaDocs)/[UserGuide]({{site.baseurl}}/AthenaDocsUG)/Boundary Conditions

Boundary conditions are applied using *ghost zones*, a set of `nghost` cells which extend beyond the grid on both
sides in each dimension.  By default `nghost=4`, but this value can be changed in `./athena/src/defs.h.in`.

Default Boundary Conditions
===========================

In the [Domain Blocks]({{site.baseurl}}/AthenaDocsUGDomainBlck) in the input file, integer flags can be set to specify a limited set
of boundary conditions automatically in Athena.  The actual implementation of
the boundary conditions uses function pointers.  The
flags are used to enroll the appropriate default
functions from the complete list in `./athena/src/bvals_mhd.c`, which set values in the ghost zones
according to the following values for the flags:

	1 = reflecting
	2 = flow out
	4 = periodic

Flow in boundary conditions, in which the values in the ghost zones are set to specific values, are set using
function pointers as described below.

New Problem-specific Boundary Conditions
========================================

The use of function pointers makes adding new boundary conditions for
specific problems quite easy.  For example, to add a new problem-specific
boundary condition along the inner X1 boundary, the user would
 1. Write a new function which sets the values in the ghost zones and include it in the same file as the problem generator, and
 2. enroll this new function by adding the following line at the end of the problem generator

	bvals_mhd_fun(pDomain, right_x1, bc_function_name);

where the arguments are

**pDomain** pointer to the `DomainS` structure.  See [Programmer Guide]({{site.baseurl}}/AthenaDocsPG).

**right_x1** specifies the boundary on which the special function is enrolled; use
` left_x1` or `right_x1` for the inner or outer x1-boundary,
and similarly for the x2- and x3-boundaries.

**bc_function_name** is the name of the special function written in step (1).

Examples of problem generators that enroll special boundary conditions include
`dmr.c`, `noh.c`, and `shkset3d.c`.

Users should also note the following:

 1. Boundary condition flags in the input file are only required for directions in which the grid is integrated.  That is, if `Nx1>1` and `Nx2=Nx3=1`, then only bc_ix1 and bc_ox1 are required in the input file.  The parameters bc_ix2, bc_ox2, bc_ix3, and bc_ox3 may be present, but their value will not be checked.

 2. If the user enrolls a boundary condition routine for say the inner x1-boundary, the boundary condition flag bc_ix1 in the parameter file is not required.  Again it may be in the parameter file, but its value will not be checked.

New Problem-specific Boundary Conditions with SMR
=================================================

With SMR, boundary conditions must be specified only on the root Domain, and only on higher level Domains that touch the edge of the root Domain
(edge of the physical region of the problem).  The boundary conditions for Domains whose edges lie entirely inside the root Domain are handled 
by the prolongation operators in the SMR algorithm.

Thus, when specifying new problem-specific BCs with SMR, the user must check whether the Domain touches the edge of the root Domain before setting
the function pointers.  For example, calls to the `bvals_mhd_fun()` should now be of the form:

	if (pDomain->Disp[1] == 0)                    bvals_mhd_fun(pDomain, left_x2,  reflect_ix2);
	if (pDomain->MaxX[1] == pDomain->RootMaxX[1]) bvals_mhd_fun(pDomain, right_x2, reflect_ox2);


In this example, boundary conditions on the left and right x2-boundaries are set, but only if the edge of the Domain is equal to the
edge of the root Domain.

Outputting Values in the Ghost Zones
====================================

Sometimes it is useful to query the values of the conserved variables in the ghost zones themselves.  This is possible by enabling the ghost-zone feature
during the configure step (before compiling), for example

	% configure --enable-ghost

After the code is compiled when configured for ghost-zone outputs, virtually all output types specified in the input file will output data over 
the ghost zones as well as the active zones.
